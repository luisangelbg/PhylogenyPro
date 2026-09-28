/* PhylogenyPro — molecular phylogenetics in the browser.
   Copyright (C) 2026  Luis Ángel Barrera-Guzmán

   This program is free software: you can redistribute it and/or modify it under
   the terms of the GNU General Public License as published by the Free Software
   Foundation, either version 3 of the License, or (at your option) any later
   version. It is distributed in the hope that it will be useful, but WITHOUT ANY
   WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
   A PARTICULAR PURPOSE. See the GNU General Public License, in the file LICENSE
   at the root of this program, or <https://www.gnu.org/licenses/>. */

/* PhylogenyPro — parsimony, tree search and consensus.

   Parsimony asks one question: which tree needs the fewest changes? It is the
   oldest criterion, the most transparent, and the one that fails in the
   Felsenstein zone (the first lab of the home page shows exactly how). It is
   still the right tool for morphology, for data with little change, and as a
   fast starting tree for maximum likelihood.

   What is here
     · Fitch (1971) for unordered characters with equal cost, on compressed
       patterns and with per-pattern weights;
     · Sankoff (1975) for a general cost matrix — ordered characters, step
       matrices, transversion weighting;
     · a heuristic search: stepwise addition from several random orders, then
       branch swapping by NNI, SPR or TBR, keeping every tree that ties for the
       best score;
     · the consistency, retention, rescaled consistency and homoplasy indices;
     · Bremer support (how many extra steps it costs to lose a clade);
     · bootstrap and jackknife, both reporting split frequencies;
     · strict, majority-rule, semistrict and greedy consensus trees.

   Validated against phangorn in validation/block4. */

/* Wrapped in a named function so that a Web Worker can be built from its own
   source text: see js/pool.js. Nothing here touches the DOM. */
function ParsCore(g) {
const Pars = {};
const Consensus = {};

  /* ================================================================
     1 · Fitch
     ================================================================ */
  /* Down-pass only: the number of steps is complete after one post-order
     traversal, which is all a search needs. States are bit masks. */
  function fitch(F, A) {
    const nPat = A.nPat, S = A.nStates;
    const sets = new Int32Array(F.n * nPat);
    const full = (1 << S) - 1;
    /* tips */
    for (let k = 0; k < F.n; k++) {
      if (!F.isTip[k]) continue;
      const row = A.tipMask[F.tipRow[k]];
      const off = k * nPat;
      for (let p = 0; p < nPat; p++) sets[off + p] = row[p];
    }
    let steps = 0;
    const perPattern = new Float64Array(nPat);
    for (let i = 0; i < F.post.length; i++) {
      const k = F.post[i];
      if (F.isTip[k]) continue;
      const kids = F.kids[k], off = k * nPat;
      if (!kids.length) continue;
      for (let p = 0; p < nPat; p++) {
        let cur = sets[kids[0] * nPat + p];
        let add = 0;
        for (let c = 1; c < kids.length; c++) {
          const other = sets[kids[c] * nPat + p];
          const inter = cur & other;
          if (inter) cur = inter;
          else { cur = cur | other; add++; }
        }
        sets[off + p] = cur;
        if (add) { perPattern[p] += add; steps += add * A.weights[p]; }
      }
    }
    return { steps, perPattern, sets };
  }

  /* ================================================================
     2 · Sankoff
     ================================================================ */
  /* A general cost matrix: cost[i][j] is what it costs to change from state i
     to state j along a branch. Ordered characters use |i − j|; transversion
     weighting uses a 4 × 4 matrix with transitions cheaper. */
  function sankoff(F, A, cost) {
    const nPat = A.nPat, S = A.nStates;
    const BIG = 1e12;
    const g = new Float64Array(F.n * nPat * S);
    for (let k = 0; k < F.n; k++) {
      if (!F.isTip[k]) continue;
      const mask = A.tipMask[F.tipRow[k]];
      for (let p = 0; p < nPat; p++) {
        const off = (k * nPat + p) * S, m = mask[p];
        for (let s = 0; s < S; s++) g[off + s] = (m & (1 << s)) ? 0 : BIG;
      }
    }
    for (let i = 0; i < F.post.length; i++) {
      const k = F.post[i];
      if (F.isTip[k]) continue;
      const kids = F.kids[k];
      for (let p = 0; p < nPat; p++) {
        const off = (k * nPat + p) * S;
        for (let s = 0; s < S; s++) {
          let tot = 0;
          for (let c = 0; c < kids.length; c++) {
            const co = (kids[c] * nPat + p) * S;
            let best = BIG;
            for (let j = 0; j < S; j++) {
              const v = cost[s][j] + g[co + j];
              if (v < best) best = v;
            }
            tot += best;
          }
          g[off + s] = tot;
        }
      }
    }
    let steps = 0;
    const perPattern = new Float64Array(nPat);
    for (let p = 0; p < nPat; p++) {
      const off = (0 * nPat + p) * S;
      let best = Infinity;
      for (let s = 0; s < S; s++) if (g[off + s] < best) best = g[off + s];
      perPattern[p] = best;
      steps += best * A.weights[p];
    }
    return { steps, perPattern };
  }

  /* cost matrices ready to use */
  function costMatrix(kind, S, opts) {
    opts = opts || {};
    const M = Array.from({ length: S }, () => new Float64Array(S));
    for (let i = 0; i < S; i++) for (let j = 0; j < S; j++) {
      if (i === j) M[i][j] = 0;
      else if (kind === 'ordered') M[i][j] = Math.abs(i - j);
      else if (kind === 'transversion') M[i][j] = isTransition(i, j) ? (opts.tsCost == null ? 1 : opts.tsCost) : (opts.tvCost == null ? 2 : opts.tvCost);
      else M[i][j] = 1;
    }
    return M;
  }

  /* ================================================================
     3 · data preparation
     ================================================================ */
  /* Parsimony works on the same compressed patterns as the likelihood engine,
     but it needs the state masks rather than the partial likelihoods. */
  function prepare(A) {
    if (A.tipMask) return A;
    const S = A.nStates, nPat = A.nPat;
    const tipMask = A.tips.map(arr => {
      const m = new Int32Array(nPat);
      for (let p = 0; p < nPat; p++) {
        let bits = 0;
        for (let s = 0; s < S; s++) if (arr[p * S + s] > 0) bits |= 1 << s;
        m[p] = bits;
      }
      return m;
    });
    A.tipMask = tipMask;
    /* which patterns can tell two topologies apart */
    const informative = new Uint8Array(nPat);
    for (let p = 0; p < nPat; p++) {
      const counts = {};
      for (let i = 0; i < tipMask.length; i++) {
        const m = tipMask[i][p];
        /* only unambiguous states count towards informativeness */
        if (m && (m & (m - 1)) === 0) counts[m] = (counts[m] || 0) + 1;
      }
      informative[p] = Object.values(counts).filter(v => v >= 2).length >= 2 ? 1 : 0;
    }
    A.informative = informative;
    return A;
  }

  function score(tree, A, opts) {
    opts = opts || {};
    prepare(A);
    const F = Tree.flatten(tree);
    return opts.cost ? sankoff(F, A, opts.cost) : fitch(F, A);
  }

  /* ================================================================
     4 · tree search
     ================================================================ */
  /* Stepwise addition: taxa are added one at a time, each into the branch that
     costs the fewest extra steps. The order is shuffled between replicates,
     which is what makes several starts worth doing. */
  function stepwise(A, order, opts) {
    opts = opts || {};
    const n = A.nSeq;
    const ord = order || Array.from({ length: n }, (_, i) => i);
    let tree = {
      children: [
        { node: { tip: ord[0], children: [] }, len: 1 },
        { node: { tip: ord[1], children: [] }, len: 1 },
        { node: { tip: ord[2], children: [] }, len: 1 },
      ],
    };
    for (let i = 3; i < n; i++) {
      const t = ord[i];
      let best = null;
      /* try inserting on every branch */
      const branches = [];
      (function collect(node, parent) {
        node.children.forEach(c => { branches.push({ parent: node, child: c }); collect(c.node, node); });
      })(tree, null);
      branches.forEach(b => {
        const cand = Tree.clone(tree);
        /* find the same branch in the clone by position */
        const idx = branches.indexOf(b);
        const cbranches = [];
        (function collect(node) {
          node.children.forEach(c => { cbranches.push({ parent: node, child: c }); collect(c.node); });
        })(cand);
        const target = cbranches[idx];
        const newNode = { children: [{ node: target.child.node, len: 1 }, { node: { tip: t, children: [] }, len: 1 }] };
        target.child.node = newNode;
        const s = score(cand, A, opts).steps;
        if (!best || s < best.steps) best = { tree: cand, steps: s };
      });
      tree = best.tree;
    }
    return tree;
  }

  /* --- branch swapping --- */
  /* every internal branch of the tree, as {node, parent} pairs */
  function internalBranches(tree) {
    const out = [];
    (function walk(node) {
      node.children.forEach(c => {
        if (c.node.children.length) out.push({ parent: node, edge: c });
        walk(c.node);
      });
    })(tree);
    return out;
  }
  /* Nearest-neighbour interchange: every sibling of an internal branch swapped
     against every child of it. On a binary tree that is the classical pair of
     rearrangements per branch; at a node of degree three there are four, and
     generating only two leaves the search stuck in a local optimum. */
  function* nniNeighbours(tree) {
    const branches = internalBranches(tree);
    for (let b = 0; b < branches.length; b++) {
      const probe = internalBranches(tree)[b];
      if (!probe) continue;
      const nSibs = probe.parent.children.length - 1;
      const nGrand = probe.edge.node.children.length;
      for (let si = 0; si < nSibs; si++) {
        for (let gi = 0; gi < nGrand; gi++) {
          const cand = Tree.clone(tree);
          const cb = internalBranches(cand)[b];
          if (!cb) continue;
          const sibs = cb.parent.children.filter(c => c !== cb.edge);
          if (si >= sibs.length || gi >= cb.edge.node.children.length) continue;
          const sib = sibs[si], grand = cb.edge.node.children[gi];
          const tmp = sib.node;
          sib.node = grand.node;
          grand.node = tmp;
          yield cand;
        }
      }
    }
  }

  /* subtree pruning and regrafting: cut a subtree and try it everywhere else.
     With reroot = true this becomes TBR, which also tries every rooting of the
     pruned subtree. */
  function* sprNeighbours(tree, reroot, maxMoves) {
    const all = [];
    (function walk(node) {
      node.children.forEach(c => { all.push(c); walk(c.node); });
    })(tree);
    let done = 0;
    for (let i = 0; i < all.length; i++) {
      /* prune the i-th branch */
      const base = Tree.clone(tree);
      const baseEdges = [];
      (function walk(node) { node.children.forEach(c => { baseEdges.push({ parent: node, edge: c }); walk(c.node); }); })(base);
      const cut = baseEdges[i];
      if (!cut) continue;
      const pruned = cut.edge.node;
      const parent = cut.parent;
      /* remove the branch and splice out the parent if it is left with one child */
      parent.children = parent.children.filter(c => c !== cut.edge);
      let rest = base;
      if (parent.children.length === 1 && parent !== base) {
        /* find the grandparent and splice */
        const gp = findParent(base, parent);
        if (gp) {
          const e = gp.children.find(c => c.node === parent);
          e.node = parent.children[0].node;
          e.len = (e.len || 0) + (parent.children[0].len || 0);
        }
      } else if (base.children.length === 1) {
        rest = base.children[0].node;
      }
      if (Tree.tips(rest).length < 2) continue;
      /* the rootings of the pruned subtree that TBR tries */
      const rootings = reroot && pruned.children.length ? tbrRootings(pruned) : [pruned];
      for (let r = 0; r < rootings.length; r++) {
        const restEdges = [];
        (function walk(node) { node.children.forEach(c => { restEdges.push({ parent: node, edge: c }); walk(c.node); }); })(rest);
        for (let j = 0; j < restEdges.length; j++) {
          const cand = Tree.clone(rest);
          const candEdges = [];
          (function walk(node) { node.children.forEach(c => { candEdges.push({ parent: node, edge: c }); walk(c.node); }); })(cand);
          const at = candEdges[j];
          const attach = Tree.clone(rootings[r]);
          const newNode = { children: [{ node: at.edge.node, len: 1 }, { node: attach, len: 1 }] };
          at.edge.node = newNode;
          yield cand;
          if (maxMoves && ++done >= maxMoves) return;
        }
      }
    }
  }
  function findParent(root, target) {
    let found = null;
    (function walk(node) {
      node.children.forEach(c => { if (c.node === target) found = node; else walk(c.node); });
    })(root);
    return found;
  }
  /* every way of rooting a pruned subtree, for TBR */
  function tbrRootings(sub) {
    const out = [sub];
    const idx = Tree.index(sub);
    Tree.nodes(sub).forEach(nd => {
      if (nd === sub) return;
      try { out.push(Tree.rerootAbove(sub, nd, 0.5)); } catch (e) { /* skip */ }
    });
    return out.slice(0, 12);          // a cap, or TBR explodes on big trees
  }

  /* The search itself. Keeps every tree that ties for the best score, up to
     maxTrees, which is what a parsimony analysis reports. */
  function search(A, opts) {
    opts = opts || {};
    prepare(A);
    const n = A.nSeq;
    const starts = opts.starts || 5;
    const swap = opts.swap || 'SPR';
    const maxTrees = opts.maxTrees || 100;
    const r = rng(opts.seed || 1);
    const t0 = (typeof performance !== 'undefined' ? performance.now() : Date.now());
    let best = Infinity, keep = [];
    const seen = new Set();
    let rearrangements = 0;

    for (let s = 0; s < starts; s++) {
      const order = Array.from({ length: n }, (_, i) => i);
      if (s > 0) shuffle(order, r);
      let tree = stepwise(A, order, opts);
      let cur = score(tree, A, opts).steps;
      /* swap until nothing improves */
      let improved = true;
      while (improved) {
        improved = false;
        const gen = swap === 'NNI' ? nniNeighbours(tree)
          : sprNeighbours(tree, swap === 'TBR', opts.maxMoves || 4000);
        for (const cand of gen) {
          rearrangements++;
          const sc = score(cand, A, opts).steps;
          if (sc < cur - 1e-9) { tree = cand; cur = sc; improved = true; break; }
          if (Math.abs(sc - cur) < 1e-9 && keep.length < maxTrees) {
            const key = Tree.writeNewick(sortTree(cand), { lengths: false, support: false });
            if (!seen.has(key) && Math.abs(cur - best) < 1e-9) { seen.add(key); keep.push(cand); }
          }
          if (opts.cancelled && opts.cancelled()) break;
        }
        if (opts.cancelled && opts.cancelled()) break;
      }
      if (cur < best - 1e-9) {
        best = cur; keep = [tree]; seen.clear();
        seen.add(Tree.writeNewick(sortTree(tree), { lengths: false, support: false }));
      } else if (Math.abs(cur - best) < 1e-9) {
        const key = Tree.writeNewick(sortTree(tree), { lengths: false, support: false });
        if (!seen.has(key) && keep.length < maxTrees) { seen.add(key); keep.push(tree); }
      }
      if (opts.onProgress) opts.onProgress(s + 1, starts, best);
      if (opts.cancelled && opts.cancelled()) break;
    }
    return {
      steps: best, trees: keep, nTrees: keep.length, rearrangements,
      ms: (typeof performance !== 'undefined' ? performance.now() : Date.now()) - t0,
      swap, starts,
    };
  }
  /* a canonical ordering, so the same topology always writes the same string */
  function sortTree(tree) {
    const t = Tree.clone(tree);
    (function walk(n) {
      n.children.forEach(c => walk(c.node));
      n.children.sort((a, b) => minTip(a.node) - minTip(b.node));
    })(t);
    return t;
  }
  function minTip(n) { return n.tip != null ? n.tip : Math.min(...n.children.map(c => minTip(c.node))); }

  /* ================================================================
     5 · the indices
     ================================================================ */
  /* m = the minimum possible number of steps for each character (states − 1),
     g = the maximum (on a star tree), s = the steps on this tree.
       CI = m/s · RI = (g − s)/(g − m) · RC = CI·RI · HI = 1 − CI            */
  function indices(tree, A, opts) {
    prepare(A);
    const sc = score(tree, A, opts);
    const nPat = A.nPat, S = A.nStates;
    let m = 0, g = 0, s = sc.steps;
    let mInf = 0, gInf = 0, sInf = 0;
    for (let p = 0; p < nPat; p++) {
      /* observed states in this pattern */
      const counts = {};
      let states = 0;
      for (let i = 0; i < A.nSeq; i++) {
        const mask = A.tipMask[i][p];
        if (mask && (mask & (mask - 1)) === 0) { counts[mask] = (counts[mask] || 0) + 1; states |= mask; }
      }
      let k = 0;
      for (let b = 0; b < S; b++) if (states & (1 << b)) k++;
      const mp = Math.max(0, k - 1);
      /* the maximum: every tip keeps its own state on a star, so the steps are
         the total minus the commonest state */
      const vals = Object.values(counts);
      const tot = vals.reduce((a, b2) => a + b2, 0);
      const gp = tot > 0 ? tot - Math.max(...vals, 0) : 0;
      m += mp * A.weights[p];
      g += gp * A.weights[p];
      if (A.informative[p]) { mInf += mp * A.weights[p]; gInf += gp * A.weights[p]; sInf += sc.perPattern[p] * A.weights[p]; }
    }
    const CI = s > 0 ? m / s : 1;
    const RI = (g - m) > 0 ? (g - s) / (g - m) : 1;
    return {
      steps: s, minSteps: m, maxSteps: g,
      CI, RI, RC: CI * RI, HI: 1 - CI,
      CIinf: sInf > 0 ? mInf / sInf : 1,
      RIinf: (gInf - mInf) > 0 ? (gInf - sInf) / (gInf - mInf) : 1,
    };
  }

  /* ================================================================
     6 · Bremer support
     ================================================================ */
  /* For every clade of the best tree: the number of extra steps of the shortest
     tree that does NOT contain that clade. Computed by searching among the
     trees visited while swapping, which is the practical way — a full
     constrained search per clade would take far longer than it is worth. */
  function bremer(bestTree, A, opts) {
    opts = opts || {};
    prepare(A);
    const n = A.nSeq;
    const bestScore = score(bestTree, A, opts).steps;
    const bestSplits = Tree.splits(bestTree, n);
    const support = new Map();
    bestSplits.forEach((node, key) => support.set(key, Infinity));
    const consider = cand => {
      const sc = score(cand, A, opts).steps;
      const sp = Tree.splits(cand, n);
      support.forEach((v, key) => {
        if (!sp.has(key)) {
          const extra = sc - bestScore;
          if (extra < v) support.set(key, extra);
        }
      });
    };
    /* every NNI and SPR neighbour of the best tree, and of the trees one step
       away from it */
    let count = 0;
    const limit = opts.limit || 20000;
    for (const cand of sprNeighbours(bestTree, false, limit)) {
      consider(cand);
      if (++count >= limit) break;
      if (opts.cancelled && opts.cancelled()) break;
    }
    for (const cand of nniNeighbours(bestTree)) consider(cand);
    const out = [];
    support.forEach((v, key) => out.push({ split: key, bremer: isFinite(v) ? v : null, node: bestSplits.get(key) }));
    return { support: out, bestScore, examined: count };
  }

  /* ================================================================
     7 · resampling
     ================================================================ */
  /* Bootstrap resamples the columns with replacement; the jackknife deletes a
     fraction of them (37 % is the usual choice, the one that makes it
     comparable to the bootstrap). Each replicate is searched from the previous
     best tree with NNI, which is what makes a hundred replicates affordable. */
  function resample(A, opts) {
    opts = opts || {};
    prepare(A);
    const reps = opts.reps || 100;
    const kind = opts.kind || 'bootstrap';
    const r = rng(opts.seed || 7);
    const n = A.nSeq;
    const counts = new Map();
    const t0 = (typeof performance !== 'undefined' ? performance.now() : Date.now());
    const starts = opts.starts || 1;
    let done = 0;
    for (let rep = 0; rep < reps; rep++) {
      /* new weights for the patterns */
      const w = new Float64Array(A.nPat);
      if (kind === 'jackknife') {
        const del = opts.deleteFraction == null ? 0.37 : opts.deleteFraction;
        for (let p = 0; p < A.nPat; p++) for (let k = 0; k < A.weights[p]; k++) if (r() >= del) w[p]++;
      } else {
        const total = A.nSites;
        /* draw `total` sites, each landing in a pattern with its own probability */
        const cum = new Float64Array(A.nPat);
        let acc = 0;
        for (let p = 0; p < A.nPat; p++) { acc += A.weights[p]; cum[p] = acc; }
        for (let i = 0; i < total; i++) {
          const u = r() * acc;
          let lo = 0, hi = A.nPat - 1;
          while (lo < hi) { const mid = (lo + hi) >> 1; if (u <= cum[mid]) hi = mid; else lo = mid + 1; }
          w[lo]++;
        }
      }
      const A2 = Object.assign({}, A, { weights: w });
      const res = search(A2, Object.assign({}, opts, { starts, seed: (opts.seed || 7) + rep, swap: opts.repSwap || 'NNI', maxTrees: 1 }));
      const tree = res.trees[0];
      if (tree) {
        const sp = Tree.splits(tree, n);
        sp.forEach((node, key) => counts.set(key, (counts.get(key) || 0) + 1));
      }
      done++;
      if (opts.onProgress) opts.onProgress(done, reps);
      if (opts.cancelled && opts.cancelled()) break;
    }
    const freq = new Map();
    counts.forEach((v, k) => freq.set(k, v / done));
    return { freq, reps: done, kind, ms: (typeof performance !== 'undefined' ? performance.now() : Date.now()) - t0 };
  }

  /* put support values on the nodes of a tree */
  function applySupport(tree, freq, nTaxa, asPercent) {
    const sp = Tree.splits(tree, nTaxa);
    sp.forEach((node, key) => {
      const v = freq.get(key);
      node.support = v == null ? 0 : (asPercent === false ? v : Math.round(v * 100));
    });
    return tree;
  }

  Object.assign(Pars, {
    fitch, sankoff, costMatrix, prepare, score, stepwise, search, indices, bremer,
    resample, applySupport, nniNeighbours, sprNeighbours, sortTree, internalBranches,
  });

  /* ================================================================
     8 · consensus trees
     ================================================================ */
  /* Splits are counted across the trees, then added to a star tree from the
     most frequent downwards, skipping any that conflicts with what is already
     there. That single routine gives every kind of consensus: strict (p = 1),
     majority rule (p = 0.5), and greedy (p = 0, add whatever is compatible). */
  function splitFrequencies(trees, nTaxa) {
    const counts = new Map();
    trees.forEach(t => {
      Tree.splits(t, nTaxa).forEach((node, key) => counts.set(key, (counts.get(key) || 0) + 1));
    });
    const freq = new Map();
    counts.forEach((v, k) => freq.set(k, v / trees.length));
    return freq;
  }
  function compatible(a, b, nTaxa) {
    const A2 = new Set(a.split(',').map(Number)), B2 = new Set(b.split(',').map(Number));
    let inter = 0;
    A2.forEach(x => { if (B2.has(x)) inter++; });
    const onlyA = A2.size - inter, onlyB = B2.size - inter;
    const outside = nTaxa - A2.size - onlyB;
    /* compatible if one of the four intersections of the two bipartitions is empty */
    return inter === 0 || onlyA === 0 || onlyB === 0 || outside === 0;
  }
  function buildFromSplits(splits, nTaxa, labels) {
    /* start from a star */
    const tips = Array.from({ length: nTaxa }, (_, i) => ({ tip: i, label: labels ? labels[i] : String(i), children: [] }));
    const root = { children: tips.map(t => ({ node: t, len: 1 })) };
    /* biggest splits first, so that nesting works */
    const sorted = splits.slice().sort((a, b) => b.set.length - a.set.length || b.freq - a.freq);
    sorted.forEach(sp => {
      const want = new Set(sp.set);
      /* the node whose tips contain the split and is deepest */
      let host = root;
      let changed = true;
      while (changed) {
        changed = false;
        for (const c of host.children) {
          const t = Tree.tips(c.node).map(x => x.tip);
          if (t.length >= want.size && sp.set.every(x => t.indexOf(x) >= 0)) { host = c.node; changed = true; break; }
        }
      }
      const inside = host.children.filter(c => Tree.tips(c.node).every(x => want.has(x.tip)));
      if (inside.length < 2 || inside.length === host.children.length) return;
      const node = { children: inside.map(c => ({ node: c.node, len: c.len })), support: sp.freq };
      host.children = host.children.filter(c => inside.indexOf(c) < 0);
      host.children.push({ node, len: 1 });
    });
    return root;
  }
  function consensus(trees, opts) {
    opts = opts || {};
    const nTaxa = opts.nTaxa || Tree.tips(trees[0]).length;
    const p = opts.p == null ? 0.5 : opts.p;
    const freq = splitFrequencies(trees, nTaxa);
    const chosen = [];
    [...freq.entries()].sort((a, b) => b[1] - a[1]).forEach(([key, f]) => {
      if (f < p - 1e-9) return;
      const set = key.split(',').map(Number);
      if (opts.greedy || p < 0.5) {
        if (chosen.some(c => !compatible(c.key, key, nTaxa))) return;
      }
      chosen.push({ key, set, freq: f });
    });
    const tree = buildFromSplits(chosen, nTaxa, opts.labels);
    return { tree, freq, splits: chosen, nTrees: trees.length };
  }
  const strict = (trees, opts) => consensus(trees, Object.assign({}, opts, { p: 1 }));
  const majority = (trees, opts) => consensus(trees, Object.assign({}, opts, { p: 0.5 }));
  const greedy = (trees, opts) => consensus(trees, Object.assign({}, opts, { p: 0, greedy: true }));
  /* semistrict (combinable components): keep every split that no tree contradicts */
  function semistrict(trees, opts) {
    opts = opts || {};
    const nTaxa = opts.nTaxa || Tree.tips(trees[0]).length;
    const all = new Map();
    const perTree = trees.map(t => Tree.splits(t, nTaxa));
    perTree.forEach(sp => sp.forEach((node, key) => all.set(key, (all.get(key) || 0) + 1)));
    const chosen = [];
    all.forEach((count, key) => {
      const conflicts = perTree.some(sp => {
        if (sp.has(key)) return false;
        return [...sp.keys()].some(k2 => !compatible(k2, key, nTaxa));
      });
      if (!conflicts) chosen.push({ key, set: key.split(',').map(Number), freq: count / trees.length });
    });
    return { tree: buildFromSplits(chosen, nTaxa, opts.labels), splits: chosen, nTrees: trees.length };
  }

  Object.assign(Consensus, { splitFrequencies, compatible, buildFromSplits, consensus, strict, majority, greedy, semistrict });
  g.Pars = Pars;
  g.Consensus = Consensus;
}
ParsCore(typeof window !== 'undefined' ? window : self);
