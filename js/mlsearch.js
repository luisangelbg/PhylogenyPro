/* PhylogenyPro — molecular phylogenetics in the browser.
   Copyright (C) 2026  Luis Ángel Barrera-Guzmán

   This program is free software: you can redistribute it and/or modify it under
   the terms of the GNU General Public License as published by the Free Software
   Foundation, either version 3 of the License, or (at your option) any later
   version. It is distributed in the hope that it will be useful, but WITHOUT ANY
   WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
   A PARTICULAR PURPOSE. See the GNU General Public License, in the file LICENSE
   at the root of this program, or <https://www.gnu.org/licenses/>. */

/* PhylogenyPro — maximum likelihood: searching tree space, and measuring how
   much the data support what it finds.

   The search follows what PhyML and RAxML do, and for the same reason: a full
   reoptimisation of every branch for every candidate topology is unaffordable,
   so a candidate is scored by reoptimising only the branches its rearrangement
   touched — five for an NNI — and the whole tree is reoptimised only when a
   round ends. The engine's incremental recomputation (js/like.js) makes each of
   those local optimisations cost a fraction of a full likelihood.

   Support
     · bootstrap ......... resample the columns, redo the search, count splits;
                           correct and slow, spread over Workers when available
     · ultrafast bootstrap  resample the per-site log-likelihoods of the trees
       (UFBoot) ..........  visited during the search and see which one wins
                            (the RELL idea of Kishino & Hasegawa 1989, as used
                            by Minh et al. 2013); nearly free, and read with a
                            different threshold (95 %, not 70 %)
     · SH-aLRT .......... Guindon et al. (2010): compare each internal branch
                          with its two NNI alternatives
     · aBayes ........... Anisimova et al. (2011): the same three likelihoods
                          read as an approximate posterior

   Topology tests: KH, SH and AU, all built on RELL resampling of the per-site
   log-likelihoods, with AU following the multiscale idea of Shimodaira (2002).

   Validated against phangorn in validation/block5. */

/* Wrapped in a named function so that a Web Worker can be built from its own
   source text: see js/pool.js. Nothing here touches the DOM. */
function MLCore(g) {
const ML = {};

  /* ================================================================
     1 · scoring a candidate topology
     ================================================================ */
  /* Branches whose node carries __hot are the ones the rearrangement touched;
     everything else keeps the length it had. The mark has two levels: 2 for the
     one or two branches the move actually created — the only ones PhyML
     reoptimises when it is merely *ranking* candidates — and 1 for the
     neighbours, which are worth optimising once a candidate is taken seriously. */
  function hotIndices(F, level) {
    const out = [];
    const min = level || 1;
    for (let k = 1; k < F.n; k++) if ((F.nodes[k].__hot || 0) >= min) out.push(k);
    return out;
  }
  function clearHot(tree) { Tree.nodes(tree).forEach(n => { delete n.__hot; }); }

  /* optimise a subset of branches, in place, with the incremental engine */
  function optimiseSome(lik, lens, which, opts) {
    opts = opts || {};
    const maxLen = opts.maxLen || 5;
    const tol = opts.tol || 1e-5;
    const iters = opts.iters || 20;
    let lnL = lik.full(lens);
    for (let pass = 0; pass < (opts.passes || 2); pass++) {
      const before = lnL;
      which.forEach(k => {
        const res = Like.brentMin(t => {
          lens[k] = t;
          const v = -lik.propose([k], lens); lik.accept(); return v;
        }, 1e-9, maxLen, tol, iters);
        lens[k] = res.x;
      });
      lnL = lik.full(lens);
      if (Math.abs(lnL - before) < (opts.stop || 1e-4)) break;
    }
    return lnL;
  }

  /* score one candidate: inherit the branch lengths, optimise the hot ones.

     opts.pool  buffers from Like.makePool(), reused between candidates. Without
                them every candidate allocates a few megabytes, which is what
                turns a search into a frozen page. Only one pooled engine may be
                alive at a time, so the `lik` returned here must not be kept.
     opts.lazy  score with the inherited lengths and no optimisation: the cheap
                first pass that decides which candidates deserve the expensive
                one (the "lazy SPR" of RAxML).
     opts.core  optimise only the branches the move created, not their
                neighbours: enough to rank candidates, and several times faster. */
  function scoreCandidate(tree, A, model, opts) {
    opts = opts || {};
    const F = Tree.flatten(tree);
    const lik = Like.engine(F, A, model, opts.pool);
    const lens = Float64Array.from(F.len, v => Math.max(1e-7, v || 0.02));
    const hot = hotIndices(F, opts.core ? 2 : 1);
    const lnL = (hot.length && !opts.lazy) ? optimiseSome(lik, lens, hot, opts) : lik.full(lens);
    return { lnL, lens, F, lik };
  }

  /* ================================================================
     2 · rearrangements that remember what they touched
     ================================================================ */
  function* nniHot(tree) {
    const n = Tree.internalBranches(tree).length;
    for (let b = 0; b < n; b++) {
      const probe = Tree.internalBranches(tree)[b];
      if (!probe) continue;
      const nSibs = Tree.nniPartners(tree, probe).length;
      const nGrand = probe.edge.node.children.length;
      for (let si = 0; si < nSibs; si++) {
        for (let gi = 0; gi < nGrand; gi++) {
          const cand = Tree.clone(tree);
          clearHot(cand);
          const cb = Tree.internalBranches(cand)[b];
          if (!cb) continue;
          const sibs = Tree.nniPartners(cand, cb);
          if (si >= sibs.length || gi >= cb.edge.node.children.length) continue;
          const sib = sibs[si], grand = cb.edge.node.children[gi];
          const tmp = sib.node; sib.node = grand.node; grand.node = tmp;
          /* the five branches an NNI changes: the internal one, the two that
             moved, and the two that stayed */
          cb.edge.node.children.forEach(c => { c.node.__hot = 1; });
          cb.parent.children.forEach(c => { c.node.__hot = 1; });
          sib.node.__hot = 1;
          grand.node.__hot = 1;
          cb.edge.node.__hot = 2;          // the internal branch: the one that matters
          yield cand;
        }
      }
    }
  }

  /* how far each edge of `rest` lies from the place the subtree was cut,
     counted in nodes along the unrooted tree. Regrafting far away is what makes
     SPR expensive and is almost never what improves a tree, so the sweep is
     capped at a radius, exactly as PhyML and RAxML do. */
  function edgesWithin(rest, anchor, radius) {
    const edges = Tree.allEdges(rest);
    if (!radius || !anchor) return edges.map((_, j) => j);
    const parentOf = new Map();
    edges.forEach(e => parentOf.set(e.edge.node, e.parent));
    const dist = new Map([[anchor, 0]]);
    let front = [anchor];
    for (let d = 1; d <= radius + 1 && front.length; d++) {
      const next = [];
      front.forEach(node => {
        const nbrs = node.children.map(c => c.node);
        const up = parentOf.get(node);
        if (up) nbrs.push(up);
        nbrs.forEach(m => { if (!dist.has(m)) { dist.set(m, d); next.push(m); } });
      });
      front = next;
    }
    const out = [];
    edges.forEach((e, j) => {
      const a = dist.has(e.parent) ? dist.get(e.parent) : Infinity;
      const b = dist.has(e.edge.node) ? dist.get(e.edge.node) : Infinity;
      if (Math.min(a, b) <= radius) out.push(j);
    });
    return out;
  }

  function* sprHot(tree, maxMoves, radius) {
    let done = 0;
    const nEdges = Tree.allEdges(tree).length;
    for (let i = 0; i < nEdges; i++) {
      const base = Tree.clone(tree);
      clearHot(base);
      const edges = Tree.allEdges(base);
      const cut = edges[i];
      if (!cut) continue;
      const pruned = cut.edge.node;
      const parent = cut.parent;
      parent.children = parent.children.filter(c => c !== cut.edge);
      let rest = base, anchor = parent;
      if (parent.children.length === 1 && parent !== base) {
        const gp = findParent(base, parent);
        if (gp) {
          const e = gp.children.find(c => c.node === parent);
          e.node = parent.children[0].node;
          e.len = (e.len || 0) + (parent.children[0].len || 0);
          e.node.__hot = 1;
          anchor = e.node;
        }
      } else if (base.children.length === 1) { rest = base.children[0].node; anchor = rest; }
      if (Tree.tips(rest).length < 2) continue;
      const allowed = edgesWithin(rest, anchor, radius);
      for (const j of allowed) {
        const cand = Tree.clone(rest);
        const at = Tree.allEdges(cand)[j];
        if (!at) continue;
        const attach = Tree.clone(pruned);
        const node = { children: [{ node: at.edge.node, len: at.edge.len || 0.02 }, { node: attach, len: cut.edge.len || 0.02 }] };
        at.edge.node = node;
        node.children.forEach(c => { c.node.__hot = 1; });
        node.__hot = 2;                    // the branch the regraft created
        attach.__hot = 2;                  // and the one that holds the subtree
        yield cand;
        if (maxMoves && ++done >= maxMoves) return;
      }
    }
  }
  function findParent(root, target) {
    let found = null;
    (function walk(node) { node.children.forEach(c => { if (c.node === target) found = node; else walk(c.node); }); })(root);
    return found;
  }

  /* ================================================================
     3 · the search
     ================================================================ */
  function search(A, spec, opts) {
    opts = opts || {};
    const t0 = (typeof performance !== 'undefined' ? performance.now() : Date.now());
    const model = Like.model(spec, A);
    /* the starting tree: whatever the caller gives, or BIONJ on JC distances */
    let tree = opts.start ? Tree.clone(opts.start) : null;
    if (!tree) {
      const D = Dist.matrix(opts.seqs || [], 'jc', { type: spec.type }).D;
      tree = Tree.bionj(D, opts.labels);
    }
    clearHot(tree);
    /* a proper fit of the starting tree, so the comparison is fair */
    let fit = Like.fit(tree, A, Object.assign({}, spec), { passes: opts.initPasses || 6, tol: opts.roughTol || 1e-3, fixModel: opts.fixModel, optFreqs: opts.optFreqs });
    tree = fit.tree;
    let best = fit.lnL;
    let bestSpec = fit.spec;
    const visited = [];              // for the ultrafast bootstrap
    let rounds = 0, evaluated = 0;
    const eps = opts.eps == null ? 0.01 : opts.eps;

    const recordVisited = (t, lnL) => {
      if (!opts.collect) return;
      if (visited.length >= (opts.maxCollect || 300)) return;
      visited.push({ tree: Tree.clone(t), lnL });
    };
    recordVisited(tree, best);

    /* the buffers every candidate borrows, and the model every candidate uses:
       both built once, because building them per candidate costs more than the
       likelihood itself */
    const nNodes = Tree.nodes(tree).length + 2;
    const pool = Like.makePool(nNodes, A, Like.model(bestSpec, A));
    let curModel = Like.model(bestSpec, A);
    const accepts = opts.filter || null;
    const keep = opts.shortlist || 12;

    let improved = true;
    while (improved && rounds < (opts.maxRounds || 40)) {
      improved = false;
      rounds++;
      /* --- NNI sweep: rank on the central branch alone, as PhyML does, and
             reoptimise the neighbourhood only for the one that wins --- */
      let bestCand = null;
      for (const cand of nniHot(tree)) {
        if (accepts && !accepts(cand)) continue;
        const sc = scoreCandidate(cand, A, curModel, { core: true, passes: 1, iters: 10, pool });
        evaluated++;
        recordVisited(cand, sc.lnL);
        if (!bestCand || sc.lnL > bestCand.lnL) bestCand = { tree: cand, lnL: sc.lnL };
        if (opts.cancelled && opts.cancelled()) break;
      }
      if (bestCand) {
        const sc = scoreCandidate(bestCand.tree, A, curModel, { passes: 2, iters: 16, pool });
        if (sc.lnL > best + eps) {
          tree = Tree.unflatten(sc.F, sc.lens);
          best = sc.lnL;
          improved = true;
        }
      }
      /* --- SPR sweep, when NNI has nothing left --- */
      if (!improved && opts.spr !== false) {
        /* first pass: every candidate within the radius, scored with the
           lengths it inherits — one likelihood each, no optimisation */
        const shortlist = [];
        for (const cand of sprHot(tree, opts.maxMoves || 1500, opts.radius == null ? 5 : opts.radius)) {
          if (accepts && !accepts(cand)) continue;
          const sc = scoreCandidate(cand, A, curModel, { lazy: true, pool });
          evaluated++;
          recordVisited(cand, sc.lnL);
          if (shortlist.length < keep) {
            shortlist.push({ tree: cand, lnL: sc.lnL });
            shortlist.sort((a, b) => a.lnL - b.lnL);
          } else if (sc.lnL > shortlist[0].lnL) {
            shortlist[0] = { tree: cand, lnL: sc.lnL };
            shortlist.sort((a, b) => a.lnL - b.lnL);
          }
          if (opts.cancelled && opts.cancelled()) break;
        }
        /* second pass: only the best few get their branches reoptimised */
        let bestSpr = null;
        for (let s = shortlist.length - 1; s >= 0; s--) {
          const sc = scoreCandidate(shortlist[s].tree, A, curModel, { passes: 2, iters: 14, pool });
          if (!bestSpr || sc.lnL > bestSpr.lnL) bestSpr = { lnL: sc.lnL, lens: sc.lens, F: sc.F };
        }
        if (bestSpr && bestSpr.lnL > best + eps) {
          tree = Tree.unflatten(bestSpr.F, bestSpr.lens);
          best = bestSpr.lnL;
          improved = true;
        }
      }
      /* --- every few rounds, reoptimise everything --- */
      if (improved && (rounds % (opts.refitEvery || 2) === 0)) {
        clearHot(tree);
        fit = Like.fit(tree, A, Object.assign({}, bestSpec), { passes: opts.refitPasses || 4, tol: opts.roughTol || 1e-3, fixModel: opts.fixModel, optFreqs: opts.optFreqs });
        tree = fit.tree; best = fit.lnL; bestSpec = fit.spec;
        curModel = fit.model || Like.model(bestSpec, A);
      }
      if (opts.onProgress) opts.onProgress(rounds, best, evaluated);
      if (opts.cancelled && opts.cancelled()) break;
    }
    /* a full, careful fit of the winner */
    clearHot(tree);
    fit = Like.fit(tree, A, Object.assign({}, bestSpec), { passes: opts.finalPasses || 10, tol: opts.finalTol, fixModel: opts.fixModel, optFreqs: opts.optFreqs });
    return {
      tree: fit.tree, lnL: fit.lnL, spec: fit.spec, model: fit.model, lens: fit.lens,
      k: fit.k, AIC: fit.AIC, AICc: fit.AICc, BIC: fit.BIC,
      rounds, evaluated, visited,
      ms: (typeof performance !== 'undefined' ? performance.now() : Date.now()) - t0,
    };
  }

  /* per-site log-likelihoods of a tree under a model, on the compressed
     patterns — the raw material of every resampling test here */
  function siteLnL(tree, A, model, pool) {
    const F = Tree.flatten(tree);
    const lik = Like.engine(F, A, model, pool);
    lik.full(F.len);
    return lik.siteLnL();
  }

  /* ================================================================
     4 · ultrafast bootstrap (RELL)
     ================================================================ */
  /* Instead of realigning and researching, each replicate draws new weights for
     the site patterns and asks which of the trees already visited has the
     highest resampled log-likelihood. That is the RELL approximation; it costs
     one dot product per tree per replicate. The split frequencies it gives are
     less biased than the classical bootstrap and must be read against 95 %,
     not 70 % (Minh, Nguyen & von Haeseler 2013). */
  function ufboot(A, trees, model, opts) {
    opts = opts || {};
    const reps = opts.reps || 1000;
    const r = rng(opts.seed || 13);
    const nTrees = trees.length;
    if (!nTrees) return { freq: new Map(), reps: 0 };
    const pool = Like.makePool(Tree.nodes(trees[0].tree || trees[0]).length + 2, A, model);
    const S = trees.map(t => siteLnL(t.tree || t, A, model, pool));
    const nPat = A.nPat;
    const w = A.weights;
    const counts = new Map();
    const nTaxa = Tree.tips(trees[0].tree || trees[0]).length;
    const splitsOf = trees.map(t => [...Tree.splits(t.tree || t, nTaxa).keys()]);
    /* cumulative pattern weights, to draw sites fast */
    const cum = new Float64Array(nPat);
    let acc = 0;
    for (let p = 0; p < nPat; p++) { acc += w[p]; cum[p] = acc; }
    const draw = new Float64Array(nPat);
    for (let rep = 0; rep < reps; rep++) {
      draw.fill(0);
      for (let i = 0; i < A.nSites; i++) {
        const u = r() * acc;
        let lo = 0, hi = nPat - 1;
        while (lo < hi) { const mid = (lo + hi) >> 1; if (u <= cum[mid]) hi = mid; else lo = mid + 1; }
        draw[lo]++;
      }
      let bestI = 0, bestV = -Infinity;
      for (let t = 0; t < nTrees; t++) {
        const s = S[t];
        let v = 0;
        for (let p = 0; p < nPat; p++) if (draw[p]) v += draw[p] * s[p];
        if (v > bestV) { bestV = v; bestI = t; }
      }
      splitsOf[bestI].forEach(key => counts.set(key, (counts.get(key) || 0) + 1));
      if (opts.onProgress && rep % 50 === 0) opts.onProgress(rep, reps);
    }
    const freq = new Map();
    counts.forEach((v, k) => freq.set(k, v / reps));
    return { freq, reps, nTrees };
  }

  /* ================================================================
     5 · SH-aLRT and aBayes
     ================================================================ */
  /* For every internal branch, the likelihood of the best tree and of its two
     NNI alternatives (with the five affected branches reoptimised). The
     SH-like statistic of Guindon et al. (2010) compares the best with the
     second best; aBayes turns the same three numbers into a posterior. */
  function branchTests(tree, A, model, opts) {
    opts = opts || {};
    const nTaxa = Tree.tips(tree).length;
    const pool = Like.makePool(Tree.nodes(tree).length + 2, A, model);
    const base = scoreCandidate(Tree.clone(tree), A, model, { passes: 3, iters: 20, pool });
    const lnL0 = base.lnL;
    const branches = Tree.internalBranches(tree);
    const out = [];
    for (let b = 0; b < branches.length; b++) {
      const alts = [];
      const probe = Tree.internalBranches(tree)[b];
      const nSibs = Tree.nniPartners(tree, probe).length;
      const nGrand = probe.edge.node.children.length;
      for (let si = 0; si < nSibs; si++) {
        for (let gi = 0; gi < nGrand; gi++) {
          const cand = Tree.clone(tree);
          clearHot(cand);
          const cb = Tree.internalBranches(cand)[b];
          if (!cb) continue;
          const sibs = Tree.nniPartners(cand, cb);
          if (si >= sibs.length || gi >= cb.edge.node.children.length) continue;
          const sib = sibs[si], grand = cb.edge.node.children[gi];
          const tmp = sib.node; sib.node = grand.node; grand.node = tmp;
          cb.edge.node.__hot = true; sib.node.__hot = true; grand.node.__hot = true;
          cb.edge.node.children.forEach(c => { c.node.__hot = true; });
          cb.parent.children.forEach(c => { c.node.__hot = true; });
          const sc = scoreCandidate(cand, A, model, { passes: 2, iters: 16, pool });
          alts.push(sc.lnL);
        }
      }
      alts.sort((a, b2) => b2 - a);
      const best2 = alts.length ? alts[0] : -Infinity;
      const diff = lnL0 - best2;
      /* aBayes: the posterior of the best of the three, with a flat prior */
      const all = [lnL0].concat(alts.slice(0, 2));
      const mx = Math.max(...all);
      const ws = all.map(v => Math.exp(v - mx));
      const sum = ws.reduce((a, b2) => a + b2, 0);
      const aBayes = ws[0] / sum;
      /* SH-aLRT: the statistic is 2Δ, and its null distribution is well
         approximated by a mixture that puts half its mass at zero; the p-value
         below follows Guindon et al. (2010, eq. 5) */
      const stat = 2 * diff;
      const p = stat <= 0 ? 1 : 0.5 * (1 - Dist.pchisq(stat, 1)) + 0.5 * (1 - Dist.pchisq(stat, 2));
      out.push({
        branch: b, lnL: lnL0, best2, delta: diff,
        shAlrt: Math.max(0, 100 * (1 - p)), p, aBayes,
        split: splitKeyOf(tree, b, nTaxa),
      });
      if (opts.onProgress) opts.onProgress(b + 1, branches.length);
      if (opts.cancelled && opts.cancelled()) break;
    }
    return { tests: out, lnL: lnL0 };
  }
  function splitKeyOf(tree, b, nTaxa) {
    const br = Tree.internalBranches(tree)[b];
    if (!br) return null;
    const set = Tree.tips(br.edge.node).map(t => t.tip).sort((x, y) => x - y);
    const other = [];
    for (let k = 0; k < nTaxa; k++) if (set.indexOf(k) < 0) other.push(k);
    return (set.length < other.length || (set.length === other.length && set[0] < other[0]) ? set : other).join(',');
  }

  /* ================================================================
     6 · topology tests: KH, SH and AU
     ================================================================ */
  /* All three work on the matrix of per-site log-likelihoods of the candidate
     trees, resampled with RELL. KH compares two trees; SH compares many while
     controlling the multiple comparison; AU corrects the selection bias that
     makes SH conservative, by resampling at several sample sizes and reading
     the slope (Shimodaira 2002). */
  function topologyTests(siteMatrix, weights, opts) {
    opts = opts || {};
    const reps = opts.reps || 1000;
    const nT = siteMatrix.length, nPat = siteMatrix[0].length;
    const r = rng(opts.seed || 23);
    const nSites = weights.reduce((a, b) => a + b, 0);
    const lnL = siteMatrix.map(s => {
      let v = 0;
      for (let p = 0; p < nPat; p++) v += weights[p] * s[p];
      return v;
    });
    const bestIdx = lnL.indexOf(Math.max(...lnL));
    const cum = new Float64Array(nPat);
    let acc = 0;
    for (let p = 0; p < nPat; p++) { acc += weights[p]; cum[p] = acc; }

    /* one RELL replicate at sample size `scale` × nSites */
    function replicate(scale) {
      const draw = new Float64Array(nPat);
      const n = Math.round(nSites * scale);
      for (let i = 0; i < n; i++) {
        const u = r() * acc;
        let lo = 0, hi = nPat - 1;
        while (lo < hi) { const mid = (lo + hi) >> 1; if (u <= cum[mid]) hi = mid; else lo = mid + 1; }
        draw[lo]++;
      }
      const v = new Float64Array(nT);
      for (let t = 0; t < nT; t++) {
        const s = siteMatrix[t];
        let x = 0;
        for (let p = 0; p < nPat; p++) if (draw[p]) x += draw[p] * s[p];
        v[t] = x / scale;             // rescaled to the original sample size
      }
      return v;
    }

    /* --- KH: pairwise, against the best tree --- */
    const khP = new Array(nT).fill(1);
    const shCount = new Array(nT).fill(0);
    const centred = [];
    const diffsKH = Array.from({ length: nT }, () => []);
    const maxCentred = [];
    for (let rep = 0; rep < reps; rep++) {
      const v = replicate(1);
      /* KH: the difference between each tree and the best, centred */
      for (let t = 0; t < nT; t++) diffsKH[t].push((v[bestIdx] - v[t]));
      /* SH: centre every tree on its own bootstrap mean, then take the max */
      centred.push(v);
    }
    /* KH p-values: how often the resampled difference exceeds the observed one
       after centring */
    for (let t = 0; t < nT; t++) {
      if (t === bestIdx) { khP[t] = 1; continue; }
      const obs = lnL[bestIdx] - lnL[t];
      const arr = diffsKH[t];
      const mean = arr.reduce((a, b) => a + b, 0) / arr.length;
      let ge = 0;
      arr.forEach(d => { if (d - mean >= obs) ge++; });
      khP[t] = ge / arr.length;
    }
    /* SH: the classical centring on the bootstrap means of every tree */
    const means = new Array(nT).fill(0);
    centred.forEach(v => { for (let t = 0; t < nT; t++) means[t] += v[t] / centred.length; });
    centred.forEach(v => {
      let mx = -Infinity;
      for (let t = 0; t < nT; t++) { const d = v[t] - means[t]; if (d > mx) mx = d; }
      for (let t = 0; t < nT; t++) {
        const obs = lnL[bestIdx] - lnL[t];
        if (mx - (v[t] - means[t]) >= obs) shCount[t]++;
      }
    });
    const shP = shCount.map(c => c / reps);

    /* --- AU: multiscale bootstrap --- */
    const scales = opts.scales || [0.5, 0.6, 0.8, 1, 1.4, 2];
    const auP = new Array(nT).fill(null);
    if (opts.au !== false) {
      const bpByScale = scales.map(sc => {
        const win = new Array(nT).fill(0);
        const n = Math.max(50, Math.round(reps / 3));
        for (let rep = 0; rep < n; rep++) {
          const v = replicate(sc);
          let bi = 0, bv = -Infinity;
          for (let t = 0; t < nT; t++) if (v[t] > bv) { bv = v[t]; bi = t; }
          win[bi]++;
        }
        return win.map(w => w / n);
      });
      for (let t = 0; t < nT; t++) {
        /* z = −Φ⁻¹(BP) fitted against √r and 1/√r: AU = 1 − Φ(d − c) */
        const xs = [], ys = [];
        scales.forEach((sc, i) => {
          const bp = Math.min(1 - 1e-6, Math.max(1e-6, bpByScale[i][t]));
          xs.push(sc);
          ys.push(-qnorm(bp));
        });
        /* fit y = d·√r + c/√r by least squares */
        let s11 = 0, s12 = 0, s22 = 0, t1 = 0, t2 = 0;
        xs.forEach((sc, i) => {
          const a = Math.sqrt(sc), b = 1 / Math.sqrt(sc);
          s11 += a * a; s12 += a * b; s22 += b * b;
          t1 += a * ys[i]; t2 += b * ys[i];
        });
        const det = s11 * s22 - s12 * s12;
        if (Math.abs(det) < 1e-12) { auP[t] = null; continue; }
        const d = (t1 * s22 - t2 * s12) / det;
        const c = (t2 * s11 - t1 * s12) / det;
        auP[t] = 1 - pnorm(d - c);
      }
    }
    return {
      lnL, best: bestIdx, diff: lnL.map(v => lnL[bestIdx] - v),
      kh: khP, sh: shP, au: auP, reps,
    };
  }
  function pnorm(z) {
    /* Abramowitz & Stegun 26.2.17 */
    const t = 1 / (1 + 0.2316419 * Math.abs(z));
    const d = 0.3989422804014327 * Math.exp(-z * z / 2);
    const p = d * t * (0.319381530 + t * (-0.356563782 + t * (1.781477937 + t * (-1.821255978 + t * 1.330274429))));
    return z >= 0 ? 1 - p : p;
  }
  function qnorm(p) {
    if (p <= 0) return -Infinity;
    if (p >= 1) return Infinity;
    let lo = -8, hi = 8;
    for (let i = 0; i < 200; i++) { const m = (lo + hi) / 2; if (pnorm(m) < p) lo = m; else hi = m; }
    return (lo + hi) / 2;
  }

  /* ================================================================
     7 · constrained search
     ================================================================ */
  /* Force a group to be monophyletic and find the best tree under that
     constraint: the honest way to test a hypothesis of monophyly, because the
     constrained tree is then compared with the unconstrained one by the
     topology tests above. */
  function searchConstrained(A, spec, constraintTips, opts) {
    opts = opts || {};
    const want = new Set(constraintTips);
    const ok = tree => {
      const nTaxa = Tree.tips(tree).length;
      if (want.size < 2 || want.size >= nTaxa) return true;
      const key = [...want].sort((a, b) => a - b).join(',');
      const other = [];
      for (let k = 0; k < nTaxa; k++) if (!want.has(k)) other.push(k);
      const alt = other.join(',');
      const sp = Tree.splits(tree, nTaxa);
      return sp.has(key) || sp.has(alt);
    };
    /* start from a tree that already obeys the constraint */
    let start = opts.start ? Tree.clone(opts.start) : null;
    if (!start || !ok(start)) start = buildConstrained(A, constraintTips, opts);
    const res = search(A, spec, Object.assign({}, opts, {
      start, collect: false,
      filter: ok,
    }));
    /* the generic search does not know about constraints, so the result is
       checked and, if it broke the constraint, rebuilt from the constrained
       start with NNI moves that keep it */
    if (ok(res.tree)) return Object.assign(res, { constrained: true });
    return constrainedNNI(A, spec, start, ok, opts);
  }
  function buildConstrained(A, constraintTips, opts) {
    const labels = opts.labels || Array.from({ length: A.nSeq }, (_, i) => 'T' + i);
    const D = opts.D || Dist.matrix(opts.seqs, 'jc', { type: opts.type || 'dna' }).D;
    const inSet = constraintTips.slice().sort((a, b) => a - b);
    const outSet = [];
    for (let i = 0; i < A.nSeq; i++) if (inSet.indexOf(i) < 0) outSet.push(i);
    const sub = idxs => {
      if (idxs.length === 1) return { tip: idxs[0], label: labels[idxs[0]], children: [] };
      const d = idxs.map(i => Float64Array.from(idxs.map(j => D[i][j])));
      const t = Tree.nj(d.map(r => Array.from(r)), idxs.map(i => labels[i]));
      /* renumber the tips back to the global indices */
      Tree.tips(t).forEach(tp => { tp.tip = idxs[tp.tip]; tp.label = labels[tp.tip]; });
      return t;
    };
    return { children: [{ node: sub(inSet), len: 0.05 }, { node: sub(outSet), len: 0.05 }] };
  }
  function constrainedNNI(A, spec, start, ok, opts) {
    const model = Like.model(spec, A);
    let tree = Tree.clone(start);
    const pool = Like.makePool(Tree.nodes(tree).length + 2, A, model);
    let fit = Like.fit(tree, A, Object.assign({}, spec), { passes: 6 });
    tree = fit.tree;
    let best = fit.lnL;
    let improved = true, rounds = 0;
    while (improved && rounds < (opts.maxRounds || 25)) {
      improved = false; rounds++;
      let bestCand = null;
      for (const cand of nniHot(tree)) {
        if (!ok(cand)) continue;
        const sc = scoreCandidate(cand, A, model, { passes: 2, iters: 14, pool });
        if (!bestCand || sc.lnL > bestCand.lnL) bestCand = { F: sc.F, lens: sc.lens, lnL: sc.lnL };
      }
      if (bestCand && bestCand.lnL > best + 0.01) {
        tree = Tree.unflatten(bestCand.F, bestCand.lens);
        best = bestCand.lnL;
        improved = true;
      }
    }
    clearHot(tree);
    fit = Like.fit(tree, A, Object.assign({}, spec), { passes: 8 });
    return Object.assign({}, fit, { constrained: true, rounds });
  }

  Object.assign(ML, {
    search, searchConstrained, buildConstrained, scoreCandidate, optimiseSome,
    nniHot, sprHot, edgesWithin, siteLnL, ufboot, branchTests, topologyTests, pnorm, qnorm, clearHot,
  });
  g.ML = ML;
}
MLCore(typeof window !== 'undefined' ? window : self);
