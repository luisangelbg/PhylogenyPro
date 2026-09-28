/* PhylogenyPro — molecular phylogenetics in the browser.
   Copyright (C) 2026  Luis Ángel Barrera-Guzmán

   This program is free software: you can redistribute it and/or modify it under
   the terms of the GNU General Public License as published by the Free Software
   Foundation, either version 3 of the License, or (at your option) any later
   version. It is distributed in the hope that it will be useful, but WITHOUT ANY
   WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
   A PARTICULAR PURPOSE. See the GNU General Public License, in the file LICENSE
   at the root of this program, or <https://www.gnu.org/licenses/>. */

/* PhylogenyPro — comparing trees, and admitting they disagree.

   Every earlier block produces a tree. This one asks what to do when two of
   them are not the same, which in real data is always.

   How far apart two trees are
     Robinson–Foulds counts the splits one has and the other has not. It is the
     standard and it is also the bluntest: moving one tip can change it by the
     maximum. The branch score of Kuhner & Felsenstein uses the lengths, the
     path difference of Steel & Penny uses the distances between tips, and the
     quartet distance counts the four-taxon statements the two trees disagree
     on — which is the one that degrades gracefully.

   Where they disagree, and how much
     Concordance factors (Minh et al. 2020) put a number on every branch of a
     reference tree: what percentage of the gene trees contains it (gCF), and
     what percentage of the sites supports it (sCF). A branch with 100 %
     bootstrap and 35 % concordance is a real thing, and a common one.

   What the disagreement might be
     A species tree estimated from quartet frequencies, which is what is left
     when incomplete lineage sorting makes concatenation lie; a split network,
     which draws the conflict instead of resolving it; and the ABBA-BABA test,
     which asks whether the conflict is symmetric, as lineage sorting predicts,
     or lopsided, as introgression does.

   Wrapped in a named function so a Web Worker can be built from its own source
   text (js/pool.js). Nothing here touches the DOM. */

function TreeCmpCore(g) {
const Cmp = {};

  /* ================================================================
     1 · splits, as membership vectors
     ================================================================ */
  /* Tree.splits keys a split by the sorted smaller side, which is what set
     operations want. For the quartet work a membership vector is faster, so
     both are built once and carried together. */
  /* The lengths are accumulated by split rather than read off one node, because
     a rooted tree has two edges under its root and they are one edge of the
     unrooted tree: their split is the same and their length is the sum. Reading
     only one of them is a quiet halving that every weighted measure inherits. */
  function splitTable(tree, n, opts) {
    opts = opts || {};
    const nTaxa = n || g.Tree.tips(tree).length;
    const F = g.Tree.flatten(tree);
    const below = new Array(F.n);
    const acc = new Map();
    for (let i = 0; i < F.post.length; i++) {
      const k = F.post[i];
      below[k] = F.isTip[k] ? [F.tipRow[k]] : [].concat.apply([], F.kids[k].map(c => below[c]));
      if (k === 0) continue;
      const size = below[k].length;
      if (size === nTaxa) continue;
      /* A split is trivial when *either* side is a single tip. Testing only the
         side below the node misses the case of a tree rooted on a tip, where
         the branch above its sister holds n−1 tips and its split is that same
         tip against the rest: counted as internal it inflates every comparison
         with a tree rooted elsewhere by one. */
      if (Math.min(size, nTaxa - size) === 1 && !opts.trivial) continue;
      const key = canonKey(below[k], nTaxa);
      acc.set(key, (acc.get(key) || 0) + (F.len[k] || 0));
    }
    const out = { keys: [], sets: [], lens: [], nTaxa };
    acc.forEach((len, key) => {
      const set = new Uint8Array(nTaxa);
      key.split(',').forEach(v => { set[+v] = 1; });
      out.keys.push(key);
      out.sets.push(set);
      out.lens.push(len);
    });
    out.byKey = new Map(out.keys.map((k, i) => [k, i]));
    return out;
  }
  /* the length of the branch a split stands for, or zero when the tree does
     not have it */
  const lenOf = (tab, key) => (tab.byKey.has(key) ? tab.lens[tab.byKey.get(key)] : 0);

  /* ================================================================
     2 · how far apart
     ================================================================ */
  function rf(t1, t2, n) {
    const a = splitTable(t1, n), b = splitTable(t2, n);
    let only1 = 0, only2 = 0;
    a.keys.forEach(k => { if (!b.byKey.has(k)) only1++; });
    b.keys.forEach(k => { if (!a.byKey.has(k)) only2++; });
    const nTaxa = a.nTaxa;
    /* the largest the distance can be for two trees of this size, which is what
       phangorn normalises by: every internal split different on both sides */
    const maxRF = 2 * (nTaxa - 3);
    return {
      rf: only1 + only2, only1, only2,
      shared: a.keys.length - only1,
      n1: a.keys.length, n2: b.keys.length,
      maxRF, normalised: maxRF > 0 ? (only1 + only2) / maxRF : 0,
    };
  }
  /* the weighted version: how much branch length the two trees put on splits
     the other does not have, plus how differently they weight the shared ones */
  function weightedRF(t1, t2, n) {
    const a = splitTable(t1, n, { trivial: true }), b = splitTable(t2, n, { trivial: true });
    const keys = new Set([...a.keys, ...b.keys]);
    let s = 0;
    keys.forEach(k => { s += Math.abs(lenOf(a, k) - lenOf(b, k)); });
    return s;
  }
  /* Kuhner & Felsenstein's (1994) branch score: the same differences, squared */
  function branchScore(t1, t2, n) {
    const a = splitTable(t1, n, { trivial: true }), b = splitTable(t2, n, { trivial: true });
    const keys = new Set([...a.keys, ...b.keys]);
    let s = 0;
    keys.forEach(k => { const d = lenOf(a, k) - lenOf(b, k); s += d * d; });
    return Math.sqrt(s);
  }
  /* Steel & Penny's (1993) path difference: for every pair of tips, how far
     apart they are in the tree — in edges, or in branch length — and then the
     Euclidean distance between the two vectors of those numbers */
  function pathVector(tree, n, weighted) {
    const nTaxa = n || g.Tree.tips(tree).length;
    /* Counted in edges, the path between two tips depends on where the root
       sits: a root inside one branch adds an edge to some paths and not to
       others. The measure is about the unrooted tree, so the root goes first. */
    const F = g.Tree.flatten(g.Tree.unroot(tree));
    const below = new Array(F.n);
    const depth = new Float64Array(F.n);
    for (let i = F.post.length - 1; i >= 0; i--) {
      const k = F.post[i], p = F.parent[k];
      if (p >= 0) depth[k] = depth[p] + (weighted ? (F.len[k] || 0) : 1);
    }
    for (let i = 0; i < F.post.length; i++) {
      const k = F.post[i];
      below[k] = F.isTip[k] ? [F.tipRow[k]] : [].concat.apply([], F.kids[k].map(c => below[c]));
    }
    const m = nTaxa * (nTaxa - 1) / 2;
    const v = new Float64Array(m);
    const at = (i, j) => (i < j ? i * nTaxa - (i * (i + 1)) / 2 + (j - i - 1)
      : j * nTaxa - (j * (j + 1)) / 2 + (i - j - 1));
    const tipDepth = new Float64Array(nTaxa);
    for (let k = 0; k < F.n; k++) if (F.isTip[k]) tipDepth[F.tipRow[k]] = depth[k];
    /* every pair meets at exactly one node: the one where they first end up in
       different daughters */
    for (let i = 0; i < F.post.length; i++) {
      const k = F.post[i];
      if (F.isTip[k]) continue;
      const kids = F.kids[k];
      for (let a = 0; a < kids.length; a++) {
        for (let b = a + 1; b < kids.length; b++) {
          below[kids[a]].forEach(u => below[kids[b]].forEach(w => {
            v[at(u, w)] = tipDepth[u] + tipDepth[w] - 2 * depth[k];
          }));
        }
      }
    }
    return v;
  }
  function pathDistance(t1, t2, n, opts) {
    opts = opts || {};
    const a = pathVector(t1, n, !!opts.weighted), b = pathVector(t2, n, !!opts.weighted);
    let s = 0;
    for (let i = 0; i < a.length; i++) { const d = a[i] - b[i]; s += d * d; }
    return Math.sqrt(s);
  }

  /* ================================================================
     3 · quartets
     ================================================================ */
  /* Four tips, three possible statements about them. For a tree the answer is
     read off its splits: the quartet is ab|cd when some split puts a and b on
     one side and c and d on the other. */
  const QUARTET = ['12|34', '13|24', '14|23'];
  function quartetOf(tab, a, b, c, d) {
    const sets = tab.sets;
    for (let s = 0; s < sets.length; s++) {
      const S = sets[s];
      const pa = S[a], pb = S[b], pc = S[c], pd = S[d];
      const k = pa + pb + pc + pd;
      if (k !== 2) continue;
      if (pa === pb) return 0;                    // ab | cd
      if (pa === pc) return 1;                    // ac | bd
      return 2;                                   // ad | bc
    }
    return -1;                                    // the tree does not resolve it
  }
  function quartetDistance(t1, t2, n) {
    const a = splitTable(t1, n), b = splitTable(t2, n);
    const nTaxa = a.nTaxa;
    let same = 0, different = 0, unresolved = 0, total = 0;
    for (let i = 0; i < nTaxa - 3; i++) {
      for (let j = i + 1; j < nTaxa - 2; j++) {
        for (let k = j + 1; k < nTaxa - 1; k++) {
          for (let l = k + 1; l < nTaxa; l++) {
            total++;
            const x = quartetOf(a, i, j, k, l), y = quartetOf(b, i, j, k, l);
            if (x < 0 || y < 0) { unresolved++; continue; }
            if (x === y) same++; else different++;
          }
        }
      }
    }
    return {
      same, different, unresolved, total,
      distance: different,
      normalised: total > 0 ? different / total : 0,
    };
  }

  /* Over a set of trees, how often each of the three statements about each
     quartet is made. That table is all a quartet species tree needs. */
  function quartetCounts(trees, n, opts) {
    opts = opts || {};
    const tabs = trees.map(t => splitTable(t, n));
    const nTaxa = tabs[0].nTaxa;
    const counts = new Map();
    const keyOf = (i, j, k, l) => i * 1e6 + j * 1e4 + k * 1e2 + l;
    for (let i = 0; i < nTaxa - 3; i++) {
      for (let j = i + 1; j < nTaxa - 2; j++) {
        for (let k = j + 1; k < nTaxa - 1; k++) {
          for (let l = k + 1; l < nTaxa; l++) {
            const c = new Float64Array(3);
            let dec = 0;
            tabs.forEach(t => {
              const q = quartetOf(t, i, j, k, l);
              if (q >= 0) { c[q]++; dec++; }
            });
            if (dec > 0) counts.set(keyOf(i, j, k, l), { c, dec, tips: [i, j, k, l] });
          }
        }
      }
    }
    return { counts, nTaxa, nTrees: trees.length, keyOf };
  }
  /* how many of those statements a candidate tree agrees with */
  function quartetScore(tree, qc, n) {
    const tab = splitTable(tree, n || qc.nTaxa);
    let score = 0, total = 0;
    qc.counts.forEach(v => {
      const [i, j, k, l] = v.tips;
      const q = quartetOf(tab, i, j, k, l);
      total += v.dec;
      if (q >= 0) score += v.c[q];
    });
    return { score, total, proportion: total > 0 ? score / total : 0 };
  }
  /* and a hill climb over nearest-neighbour interchanges to improve it. This is
     a heuristic, not ASTRAL's exact dynamic program over a constrained set of
     clades, and the block says so. */
  function quartetSearch(start, qc, n, opts) {
    opts = opts || {};
    const maxRounds = opts.maxRounds || 30;
    let best = g.Tree.clone(start);
    let bestScore = quartetScore(best, qc, n).score;
    const visited = [bestScore];
    for (let round = 0; round < maxRounds; round++) {
      let improved = false;
      /* Tree.nniMoves is a generator, so the first move that improves the score
         is taken and the rest are never built: on a twenty-taxon tree that is
         the difference between a second and a minute */
      for (const cand of g.Tree.nniMoves(best)) {
        const s = quartetScore(cand, qc, n).score;
        if (s > bestScore + 1e-9) { best = cand; bestScore = s; improved = true; break; }
      }
      visited.push(bestScore);
      if (!improved) break;
    }
    const final = quartetScore(best, qc, n);
    return { tree: best, score: final.score, total: final.total,
      proportion: final.proportion, rounds: visited.length - 1, trace: visited };
  }

  /* ================================================================
     4 · concordance factors
     ================================================================ */
  /* An internal branch of the reference tree divides the tips in two; each side
     divides again at the nodes the branch joins, so every internal branch has
     four groups around it. That quartet is what both concordance factors are
     about: a gene tree or a site can only speak about the branch if it has a
     tip in each of the four groups, and then it says one of three things. */
  function branchQuartets(tree, n) {
    const nTaxa = n || g.Tree.tips(tree).length;
    const F = g.Tree.flatten(tree);
    const below = new Array(F.n);
    for (let i = 0; i < F.post.length; i++) {
      const k = F.post[i];
      below[k] = F.isTip[k] ? [F.tipRow[k]] : [].concat.apply([], F.kids[k].map(c => below[c]));
    }
    const all = [];
    for (let i = 0; i < nTaxa; i++) all.push(i);
    const out = [];
    const rootPair = (F.kids[0] || []).filter(c => !F.isTip[c]);
    for (let k = 1; k < F.n; k++) {
      if (F.isTip[k]) continue;
      /* the two branches under a bifurcating root are one branch of the
         unrooted tree, so only one of them is counted */
      if (F.parent[k] === 0 && rootPair.length === 2 && k === rootPair[1]) continue;
      const kids = F.kids[k];
      if (kids.length !== 2) continue;
      const X1 = below[kids[0]], X2 = below[kids[1]];
      const inside = new Set(below[k]);
      const outside = all.filter(x => !inside.has(x));
      /* The other side divides at the parent: the sister and everything above.
         Except directly under a bifurcating root, where "everything above" is
         empty — the root is not a node of the unrooted tree — and the division
         happens at the sister's own two daughters instead. Missing that case
         costs exactly one branch, and it is always the deepest one. */
      const p = F.parent[k];
      const sibs = F.kids[p].filter(c => c !== k);
      let Y1 = [], Y2 = [];
      if (p === 0 && sibs.length === 1) {
        const sib = sibs[0];
        if (F.isTip[sib] || F.kids[sib].length < 2) continue;
        Y1 = below[F.kids[sib][0]];
        Y2 = below[F.kids[sib][1]];
      } else if (sibs.length >= 1) {
        Y1 = below[sibs[0]];
        Y2 = outside.filter(x => Y1.indexOf(x) < 0);
      }
      if (!X1.length || !X2.length || !Y1.length || !Y2.length) continue;
      out.push({
        node: k, key: below[k].slice().sort((a, b) => a - b).join(','),
        clade: below[k].slice().sort((a, b) => a - b),
        groups: [X1, X2, Y1, Y2], nTips: below[k].length,
        length: F.len[k] || 0,
      });
    }
    return out;
  }
  /* the split key in Tree.splits' canonical form, so the two can be compared */
  function canonKey(members, nTaxa) {
    const a = members.slice().sort((x, y) => x - y);
    const b = [];
    for (let k = 0; k < nTaxa; k++) if (a.indexOf(k) < 0) b.push(k);
    return (a.length < b.length || (a.length === b.length && a[0] < b[0]) ? a : b).join(',');
  }

  /* gCF: among the gene trees decisive for a branch, the percentage that have
     it, and the two alternative resolutions */
  function gcf(refTree, geneTrees, n, opts) {
    opts = opts || {};
    const nTaxa = n || g.Tree.tips(refTree).length;
    const branches = branchQuartets(refTree, nTaxa);
    const tabs = geneTrees.map(t => splitTable(t, nTaxa));
    /* which tips each gene tree actually has, so that a gene missing taxa is
       counted as undecided rather than as disagreeing */
    const present = geneTrees.map(t => {
      const v = new Uint8Array(nTaxa);
      g.Tree.tips(t).forEach(x => { v[x.tip] = 1; });
      return v;
    });
    const rows = branches.map(br => {
      let conc = 0, d1 = 0, d2 = 0, decisive = 0, paraphyletic = 0;
      tabs.forEach((tab, gi) => {
        const has = present[gi];
        /* one tip from each of the four groups, chosen among those the gene
           tree has; if any group is empty here the gene cannot speak */
        const pick = br.groups.map(grp => grp.filter(x => has[x]));
        if (pick.some(p => !p.length)) return;
        decisive++;
        /* the gene tree's statement is read on one representative quartet per
           group combination, and the majority over a sample of them is taken */
        const reps = opts.reps || 1;
        const tally = [0, 0, 0];
        for (let r = 0; r < reps; r++) {
          const a = pick[0][r % pick[0].length], b = pick[1][r % pick[1].length];
          const c = pick[2][r % pick[2].length], d = pick[3][r % pick[3].length];
          const srt = [a, b, c, d].slice().sort((x, y) => x - y);
          const q = quartetOf(tab, srt[0], srt[1], srt[2], srt[3]);
          if (q < 0) continue;
          /* translate the quartet's answer into "with the reference or not" */
          const pos = x => srt.indexOf(x);
          const together = q === 0 ? [[srt[0], srt[1]], [srt[2], srt[3]]]
            : q === 1 ? [[srt[0], srt[2]], [srt[1], srt[3]]]
              : [[srt[0], srt[3]], [srt[1], srt[2]]];
          const isRef = together.some(pr => (pr[0] === a && pr[1] === b) || (pr[0] === b && pr[1] === a));
          if (isRef) tally[0]++;
          else {
            /* which of the two alternatives: a with c, or a with d */
            const withC = together.some(pr => (pr[0] === a && pr[1] === c) || (pr[0] === c && pr[1] === a));
            if (withC) tally[1]++; else tally[2]++;
          }
          void pos;
        }
        const m = tally.indexOf(Math.max.apply(null, tally));
        if (tally[m] === 0) { paraphyletic++; return; }
        if (m === 0) conc++; else if (m === 1) d1++; else d2++;
      });
      const tot = decisive || 1;
      return {
        node: br.node, clade: br.clade, key: canonKey(br.clade, nTaxa), nTips: br.nTips,
        length: br.length, decisive,
        gCF: 100 * conc / tot, gDF1: 100 * d1 / tot, gDF2: 100 * d2 / tot,
        gDFP: 100 * paraphyletic / tot,
        nConcordant: conc, nDiscordant1: d1, nDiscordant2: d2,
      };
    });
    return rows;
  }

  /* sCF: the same question asked of the sites. For each branch a number of
     quartets is drawn, one tip from each of its four groups, and every site
     that is informative for that quartet votes for one of the three
     resolutions. */
  function scf(refTree, seqs, n, opts) {
    opts = opts || {};
    const nTaxa = n || g.Tree.tips(refTree).length;
    const nQuartets = opts.quartets || 100;
    const r = g.rng(opts.seed || 12345);
    const branches = branchQuartets(refTree, nTaxa);
    const L = seqs[0] ? seqs[0].length : 0;
    const isBase = ch => ch === 'A' || ch === 'C' || ch === 'G' || ch === 'T';
    return branches.map(br => {
      let sum0 = 0, sum1 = 0, sum2 = 0, used = 0, totalSites = 0;
      for (let q = 0; q < nQuartets; q++) {
        const pick = br.groups.map(grp => grp[Math.floor(r() * grp.length)]);
        const [a, b, c, d] = pick;
        if (new Set(pick).size < 4) continue;
        const sa = seqs[a], sb = seqs[b], sc = seqs[c], sd = seqs[d];
        if (!sa || !sb || !sc || !sd) continue;
        let n0 = 0, n1 = 0, n2 = 0;
        for (let i = 0; i < L; i++) {
          const A = sa[i], B = sb[i], C = sc[i], D = sd[i];
          if (!isBase(A) || !isBase(B) || !isBase(C) || !isBase(D)) continue;
          /* a site is informative for a quartet when it has exactly two
             characters, two of each: that is the only pattern that says
             anything about which pair goes together */
          if (A === B && C === D && A !== C) n0++;
          else if (A === C && B === D && A !== B) n1++;
          else if (A === D && B === C && A !== B) n2++;
        }
        const tot = n0 + n1 + n2;
        if (!tot) continue;
        used++;
        totalSites += tot;
        sum0 += n0 / tot; sum1 += n1 / tot; sum2 += n2 / tot;
      }
      const u = used || 1;
      return {
        node: br.node, clade: br.clade, key: canonKey(br.clade, nTaxa), nTips: br.nTips,
        length: br.length, quartets: used, sites: totalSites,
        sCF: 100 * sum0 / u, sDF1: 100 * sum1 / u, sDF2: 100 * sum2 / u,
      };
    });
  }

  /* ================================================================
     5 · split networks
     ================================================================ */
  /* A consensus network keeps every split that enough of the trees have,
     including splits that contradict each other — which is the point: a
     consensus tree has to throw those away, and a network does not. */
  function consensusNetwork(trees, n, opts) {
    opts = opts || {};
    const threshold = opts.threshold == null ? 0.3 : opts.threshold;
    const nTaxa = n || g.Tree.tips(trees[0]).length;
    const seen = new Map();
    trees.forEach(t => {
      const tab = splitTable(t, nTaxa);
      tab.keys.forEach((k, i) => {
        if (!seen.has(k)) seen.set(k, { key: k, count: 0, lengthSum: 0, set: tab.sets[i] });
        const e = seen.get(k);
        e.count++;
        e.lengthSum += tab.lens[i];
      });
    });
    const out = [];
    seen.forEach(e => {
      const freq = e.count / trees.length;
      if (freq + 1e-12 < threshold) return;
      out.push({
        key: e.key, set: e.set, frequency: freq, count: e.count,
        weight: e.lengthSum / e.count,
        members: e.key.split(',').map(Number),
      });
    });
    out.sort((a, b) => b.frequency - a.frequency || a.members.length - b.members.length);
    /* which splits contradict which: two splits are compatible when one of the
       four intersections of their sides is empty */
    const conflicts = out.map(() => []);
    for (let i = 0; i < out.length; i++) {
      for (let j = i + 1; j < out.length; j++) {
        if (!compatible(out[i].set, out[j].set, nTaxa)) { conflicts[i].push(j); conflicts[j].push(i); }
      }
    }
    out.forEach((s, i) => { s.conflictsWith = conflicts[i]; });
    return { splits: out, nTaxa, nTrees: trees.length, threshold,
      compatible: out.every((s, i) => conflicts[i].length === 0) };
  }
  function compatible(A, B, n) {
    let ab = 0, aB = 0, Ab = 0, AB = 0;
    for (let i = 0; i < n; i++) {
      if (A[i]) { if (B[i]) AB++; else Ab++; } else { if (B[i]) aB++; else ab++; }
    }
    return ab === 0 || aB === 0 || Ab === 0 || AB === 0;
  }

  /* A circular ordering of the tips: the order they come out in when the tree
     is drawn in a plane without crossings. Any tree gives one, and it is what
     turns a distance matrix into a split system that can be drawn. */
  function circularOrder(tree) {
    const out = [];
    (function walk(nd) {
      if (nd.tip != null) { out.push(nd.tip); return; }
      nd.children.forEach(c => walk(c.node));
    })(tree);
    return out;
  }
  /* The circular decomposition (Bandelt & Dress 1992): given a circular
     ordering, the weight of every interval split follows in closed form from
     the distances, and the reconstruction is exact when the distance really is
     a sum of those splits — which it is, for instance, when it comes from a
     tree whose drawing gave the ordering. Negative weights mean the data are
     not circular in that ordering; they are clipped to zero and counted, so
     that the figure never draws an edge of negative length and the caller can
     say how much was dropped. */
  function circularSplits(D, order, opts) {
    opts = opts || {};
    const n = order.length;
    const pos = new Array(n);
    order.forEach((t, i) => { pos[t] = i; });
    /* the indices wrap round the circle, and JavaScript's % keeps the sign of
       the dividend, so −1 % 16 is −1 and not 15 */
    const mod = i => ((i % n) + n) % n;
    const d = (i, j) => D[order[mod(i)]][order[mod(j)]];
    const out = [];
    let clipped = 0, clippedMass = 0;
    const seenKey = new Set();
    /* the scale below which a weight is rounding noise rather than a split */
    let scale = 0;
    for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) if (D[i][j] > scale) scale = D[i][j];
    const eps = opts.minWeight == null ? Math.max(1e-12, 1e-9 * scale) : opts.minWeight;
    for (let i = 0; i < n; i++) {
      for (let len = 1; len < n; len++) {
        const j = i + len - 1;
        if (len > n - 1) continue;
        /* the split {x_i .. x_j}: only one of each complementary pair is kept */
        const members = [];
        for (let k = i; k <= j; k++) members.push(order[k % n]);
        if (members.length > n - members.length) continue;
        /* when a split is exactly half the tips, both its sides are intervals
           of the circle and the loop meets it twice; the canonical key is what
           tells them apart */
        const kk = canonKey(members, n);
        if (seenKey.has(kk)) continue;
        /* The split {x_i … x_j} is bounded by the gap before x_i and the gap
           after x_j, and its weight is the four-point difference across those
           two gaps: the two distances that cross both, minus the two that cross
           one each. Getting the four terms in the wrong pairs gives a number
           that looks like a weight and is not one. */
        const w = 0.5 * (d(i - 1, j) + d(i, j + 1) - d(i - 1, j + 1) - d(i, j));
        /* A weight of −3·10⁻¹⁷ is a zero that rounded the wrong way, not a
           distance that fails to be circular; counting those as clipped would
           report fifteen failures on a tree that decomposes perfectly. */
        if (w < -eps) { clipped++; clippedMass += -w; continue; }
        if (w < eps) continue;
        const set = new Uint8Array(n);
        members.forEach(t => { set[t] = 1; });
        out.push({ members: members.slice().sort((a, b) => a - b), set, weight: w,
          start: i % n, end: j % n, size: members.length,
          key: kk });
        seenKey.add(kk);
      }
    }
    out.sort((a, b) => b.weight - a.weight);
    return { splits: out, order, nTaxa: n, clipped, clippedMass };
  }
  /* the distance a split system implies, so that the decomposition can be
     checked against the distances it came from */
  function splitDistance(splits, n) {
    const D = Array.from({ length: n }, () => new Float64Array(n));
    splits.forEach(s => {
      for (let i = 0; i < n; i++) {
        for (let j = i + 1; j < n; j++) {
          if (s.set[i] !== s.set[j]) { D[i][j] += s.weight; D[j][i] += s.weight; }
        }
      }
    });
    return D;
  }

  /* The outer cycle of the network. Every split is a band of parallel edges;
     walking once round the tips in circular order crosses each band exactly
     twice, which is what makes the conflicting splits show up as boxes rather
     than as a decision. */
  function networkOuterCycle(cs) {
    const n = cs.nTaxa, order = cs.order;
    const splits = cs.splits;
    /* each split spans an arc of the circle; its edges run perpendicular to the
       line that separates the two sides */
    const angleOf = s => {
      let sx = 0, sy = 0;
      for (let k = 0; k < n; k++) {
        if (!s.set[order[k]]) continue;
        const a = 2 * Math.PI * k / n;
        sx += Math.cos(a); sy += Math.sin(a);
      }
      return Math.atan2(sy, sx);
    };
    splits.forEach(s => { s.angle = angleOf(s); });
    /* at each gap between two consecutive tips, the splits whose arc begins or
       ends there are crossed */
    const gaps = Array.from({ length: n }, () => []);
    splits.forEach((s, i) => {
      const k1 = (s.start - 1 + n) % n;              // entering the arc
      const k2 = s.end;                               // leaving it
      gaps[k1].push({ i, sign: +1 });
      gaps[k2].push({ i, sign: -1 });
    });
    const pts = [];
    const tipAt = new Array(n);
    let x = 0, y = 0;
    for (let k = 0; k < n; k++) {
      tipAt[order[k]] = pts.length;
      pts.push({ x, y, tip: order[k] });
      /* cross the bands in angular order, so the polygon stays convex */
      const here = gaps[k].slice().sort((a, b) => {
        const aa = splits[a.i].angle + (a.sign > 0 ? 0 : Math.PI);
        const bb = splits[b.i].angle + (b.sign > 0 ? 0 : Math.PI);
        return aa - bb;
      });
      here.forEach(e => {
        const s = splits[e.i];
        const dir = s.angle + Math.PI / 2;
        x += e.sign * s.weight * Math.cos(dir);
        y += e.sign * s.weight * Math.sin(dir);
        pts.push({ x, y, split: e.i, sign: e.sign });
      });
    }
    /* the walk closes on itself up to rounding; the residue is reported so a
       drawing that does not close says so instead of looking like a network */
    const closure = Math.hypot(x, y);
    return { points: pts, tipAt, closure, splits };
  }

  /* ================================================================
     6 · ABBA-BABA
     ================================================================ */
  /* Patterson's D asks whether the two kinds of site that lineage sorting
     produces in equal numbers really are equal. With the tree (((P1,P2),P3),O),
     a site where P2 and P3 share a derived state (ABBA) and one where P1 and P3
     do (BABA) are both incongruent, and under sorting alone they are equally
     likely. A lopsided count is what introgression leaves behind — or what a
     badly chosen outgroup, or unequal rates, can imitate. */
  function dStatistic(seqs, p1, p2, p3, out, opts) {
    opts = opts || {};
    const blocks = opts.blocks || 20;
    const A = seqs[p1], B = seqs[p2], C = seqs[p3], O = seqs[out];
    const L = Math.min(A.length, B.length, C.length, O.length);
    const isBase = ch => ch === 'A' || ch === 'C' || ch === 'G' || ch === 'T';
    const abbaAt = new Float64Array(L), babaAt = new Float64Array(L);
    let abba = 0, baba = 0, used = 0;
    for (let i = 0; i < L; i++) {
      const a = A[i], b = B[i], c = C[i], o = O[i];
      if (!isBase(a) || !isBase(b) || !isBase(c) || !isBase(o)) continue;
      used++;
      if (a === o && b === c && a !== b) { abba++; abbaAt[i] = 1; }
      else if (b === o && a === c && a !== b) { baba++; babaAt[i] = 1; }
    }
    const D = (abba + baba) > 0 ? (abba - baba) / (abba + baba) : 0;
    /* a block jackknife, because neighbouring sites are not independent: drop
       one contiguous block at a time and see how much the answer moves */
    const nb = Math.max(2, Math.min(blocks, L));
    const edges = [];
    for (let k = 0; k <= nb; k++) edges.push(Math.round(k * L / nb));
    const partial = [];
    for (let k = 0; k < nb; k++) {
      let ab = abba, ba = baba;
      for (let i = edges[k]; i < edges[k + 1]; i++) { ab -= abbaAt[i]; ba -= babaAt[i]; }
      partial.push((ab + ba) > 0 ? (ab - ba) / (ab + ba) : 0);
    }
    const mean = partial.reduce((a, b) => a + b, 0) / nb;
    let s2 = 0;
    partial.forEach(v => { s2 += (v - mean) * (v - mean); });
    const se = Math.sqrt((nb - 1) / nb * s2);
    const z = se > 0 ? D / se : 0;
    return {
      ABBA: abba, BABA: baba, D, se, z, blocks: nb, sites: used,
      p: 2 * (1 - normalCdf(Math.abs(z))),
      f4: used > 0 ? (abba - baba) / used : 0,
    };
  }
  function normalCdf(x) {
    /* Abramowitz & Stegun 26.2.17, which is plenty for a z score */
    const t = 1 / (1 + 0.2316419 * Math.abs(x));
    const d = 0.3989422804014327 * Math.exp(-x * x / 2);
    const p = d * t * (0.319381530 + t * (-0.356563782 + t * (1.781477937 +
      t * (-1.821255978 + t * 1.330274429))));
    return x > 0 ? 1 - p : p;
  }
  /* every trio against one outgroup, which is how the test is usually run when
     nobody has a particular hypothesis yet */
  function dScan(seqs, taxa, outgroup, opts) {
    opts = opts || {};
    const rows = [];
    const idx = taxa.map((t, i) => i).filter(i => i !== outgroup);
    for (let a = 0; a < idx.length; a++) {
      for (let b = a + 1; b < idx.length; b++) {
        for (let c = 0; c < idx.length; c++) {
          if (c === a || c === b) continue;
          const r = dStatistic(seqs, idx[a], idx[b], idx[c], outgroup, opts);
          /* `|| 20` would turn a deliberate zero into twenty, which is the
             difference between "show me every trio" and "show me none" */
          const minSites = opts.minSites == null ? 20 : opts.minSites;
          if (r.ABBA + r.BABA < minSites) continue;
          rows.push(Object.assign({ p1: idx[a], p2: idx[b], p3: idx[c] }, r));
        }
      }
    }
    rows.sort((x, y) => Math.abs(y.z) - Math.abs(x.z));
    return rows;
  }

  /* ================================================================
     7 · tanglegrams
     ================================================================ */
  /* Two trees facing each other with a line per tip. The lines cross wherever
     the trees disagree, and a tanglegram is only readable once the rotations
     that do not change either tree have been used to make them cross as little
     as possible. */
  function leafOrder(tree) {
    const out = [];
    (function walk(nd) {
      if (nd.tip != null) { out.push(nd.tip); return; }
      nd.children.forEach(c => walk(c.node));
    })(tree);
    return out;
  }
  function crossings(order1, order2) {
    const pos = new Map(order1.map((t, i) => [t, i]));
    const seq = order2.map(t => (pos.has(t) ? pos.get(t) : -1)).filter(v => v >= 0);
    let c = 0;
    for (let i = 0; i < seq.length; i++) {
      for (let j = i + 1; j < seq.length; j++) if (seq[i] > seq[j]) c++;
    }
    return c;
  }
  /* rotate the daughters of the second tree, node by node, keeping whichever
     way round crosses less: the standard step-wise untangling, which is not
     guaranteed optimal and is what every program does */
  function untangle(fixed, movable, opts) {
    opts = opts || {};
    const rounds = opts.rounds || 4;
    const t2 = g.Tree.clone(movable);
    const o1 = leafOrder(fixed);
    let best = crossings(o1, leafOrder(t2));
    for (let r = 0; r < rounds; r++) {
      let improved = false;
      (function walk(nd) {
        if (nd.tip != null) return;
        if (nd.children.length > 1) {
          nd.children.reverse();
          const c = crossings(o1, leafOrder(t2));
          if (c < best) { best = c; improved = true; } else nd.children.reverse();
        }
        nd.children.forEach(c => walk(c.node));
      })(t2);
      if (!improved) break;
    }
    return { tree: t2, crossings: best,
      before: crossings(o1, leafOrder(movable)),
      order1: o1, order2: leafOrder(t2) };
  }

  Object.assign(Cmp, {
    splitTable, rf, weightedRF, branchScore, pathVector, pathDistance,
    QUARTET, quartetOf, quartetDistance, quartetCounts, quartetScore, quartetSearch,
    branchQuartets, canonKey, gcf, scf,
    consensusNetwork, compatible, circularOrder, circularSplits, splitDistance,
    networkOuterCycle,
    dStatistic, dScan, normalCdf,
    leafOrder, crossings, untangle,
  });
  g.Cmp = Cmp;
}
TreeCmpCore(typeof window !== 'undefined' ? window : self);
