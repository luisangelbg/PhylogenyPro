/* PhylogenyPro — molecular phylogenetics in the browser.
   Copyright (C) 2026  Luis Ángel Barrera-Guzmán

   This program is free software: you can redistribute it and/or modify it under
   the terms of the GNU General Public License as published by the Free Software
   Foundation, either version 3 of the License, or (at your option) any later
   version. It is distributed in the hope that it will be useful, but WITHOUT ANY
   WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
   A PARTICULAR PURPOSE. See the GNU General Public License, in the file LICENSE
   at the root of this program, or <https://www.gnu.org/licenses/>. */

/* PhylogenyPro — Bayesian inference: Metropolis–Hastings over trees.

   Maximum likelihood answers "which tree makes the data most probable"; this
   answers "given the data and what I assumed beforehand, how probable is each
   tree". The machinery is the one MrBayes and BEAST use:

     · a state = (topology, branch lengths, model parameters);
     · a prior on every one of those — independent exponentials on the branch
       lengths, a uniform distribution over labelled topologies, exponentials on
       the relative rates and on the shape of the gamma, a Dirichlet on the base
       frequencies, a uniform on the proportion of invariable sites;
     · proposals with their Hastings ratios written down, not guessed: the
       multiplier move of Larget & Simon (1999), a stochastic NNI, an SPR whose
       Jacobian is T/(a+b) because it merges two branches and splits one, and a
       Dirichlet proposal on the frequencies;
     · Metropolis coupling (MC³, Geyer 1991): several chains at different
       temperatures, only the cold one sampled, with swaps between neighbours,
       so the sampler can cross the valleys between tree islands;
     · marginal likelihoods by stepping stone (Xie et al. 2011) when Bayes
       factors are wanted.

   The proposals are checked the way they should be: run the sampler with the
   likelihood switched off and the samples must reproduce the prior — uniform
   over topologies, exponential over branch lengths. A wrong Hastings ratio
   fails that test immediately. See validation/block6.

   Wrapped in a named function so a Web Worker can be built from its own source
   text (js/pool.js). Nothing here touches the DOM. */

function McmcCore(g) {
const Mcmc = {};

  /* ================================================================
     1 · priors
     ================================================================ */
  function defaultPriors(spec) {
    return {
      branchMean: 0.1,                 // exponential, mean 0.1 substitutions/site
      alphaMean: 1,                    // exponential on the gamma shape
      rateMean: 1,                     // exponential on each free relative rate
      freqsDirichlet: 1,               // Dirichlet(1,…,1): uniform on the simplex
      pInv: [0, 1],                    // uniform
      topology: 'uniform',
    };
  }

  const logExp = (x, mean) => (x <= 0 ? -Infinity : -Math.log(mean) - x / mean);

  function logBranchPrior(lens, nEdges, pr) {
    let s = 0;
    for (let k = 1; k < nEdges; k++) s += logExp(lens[k], pr.branchMean);
    return s;
  }

  function logModelPrior(spec, pr, DNA_SCHEMES) {
    let s = 0;
    if (spec.type === 'dna' && spec.rates) {
      const scheme = DNA_SCHEMES[spec.model] || DNA_SCHEMES.GTR;
      const nFree = Math.max.apply(null, scheme);
      for (let r = 1; r <= nFree; r++) s += logExp(spec.rates[r], pr.rateMean);
    }
    if (spec.alpha != null) s += logExp(spec.alpha, pr.alphaMean);
    if (spec.pInv != null) {
      if (spec.pInv < pr.pInv[0] || spec.pInv > pr.pInv[1]) return -Infinity;
      s += -Math.log(pr.pInv[1] - pr.pInv[0]);
    }
    /* a Dirichlet(a,…,a) on the frequencies; with a = 1 it is a constant, kept
       explicit so that a different concentration can be used */
    if (spec.freqs && pr.freqsDirichlet !== 1) {
      const a = pr.freqsDirichlet;
      for (let i = 0; i < spec.freqs.length; i++) s += (a - 1) * Math.log(spec.freqs[i]);
    }
    return s;
  }

  /* ================================================================
     2 · one chain
     ================================================================ */
  /* beta is the heating power: 1 for the cold chain, less for the heated ones.
     Both the likelihood and the prior are raised to it, as MrBayes does. */
  function makeChain(A, spec, tree, pr, opts) {
    const c = {
      A, pr, beta: opts.beta == null ? 1 : opts.beta,
      spec: g.Like.cloneSpec(spec),
      /* the state is an unrooted tree: a trifurcating root, 2n − 3 branches */
      tree: g.Tree.unroot(g.Tree.clone(tree)),
      noLikelihood: !!opts.noLikelihood,
      pool: opts.pool || null,
    };
    rebuild(c);
    c.lnL = c.noLikelihood ? 0 : c.lik.full(c.lens);
    c.lnPrior = logBranchPrior(c.lens, c.F.n, pr) + logModelPrior(c.spec, pr, g.Like.DNA_SCHEMES);
    return c;
  }
  /* the flattened tree, the engine and the lengths, after a topology change */
  function rebuild(c) {
    c.F = g.Tree.flatten(c.tree);
    c.lens = Float64Array.from(c.F.len, v => Math.max(1e-8, v || 0.01));
    c.M = g.Like.model(c.spec, c.A);
    if (!c.noLikelihood) c.lik = g.Like.engine(c.F, c.A, c.M, c.pool);
  }
  function target(c) { return c.beta * (c.lnL + c.lnPrior); }

  /* ================================================================
     3 · the proposals
     ================================================================ */
  /* Every move returns { logHastings, undo } or null when it cannot apply.
     `undo` puts the chain back exactly as it was, which is cheaper and safer
     than cloning the whole state on every generation. */

  /* the multiplier of Larget & Simon (1999): m = exp(lambda (u − 1/2)),
     Hastings ratio m, so log m */
  function multiplier(r, lambda) { return Math.exp(lambda * (r() - 0.5)); }

  /* one branch, multiplied */
  function moveBranch(c, r, tune) {
    const n = c.F.n;
    if (n < 3) return null;
    const k = 1 + Math.floor(r() * (n - 1));
    const m = multiplier(r, tune.branchLambda);
    const old = c.lens[k];
    const nu = Math.min(20, Math.max(1e-8, old * m));
    c.lens[k] = nu;
    const dPrior = logExp(nu, c.pr.branchMean) - logExp(old, c.pr.branchMean);
    let lnLnew = c.lnL;
    if (!c.noLikelihood) { lnLnew = c.lik.propose([k], c.lens); }
    return {
      logHastings: Math.log(nu / old),
      dLnL: lnLnew - c.lnL, dPrior,
      accept() { if (!c.noLikelihood) c.lik.accept(); },
      undo() { c.lens[k] = old; if (!c.noLikelihood) c.lik.reject(); },
    };
  }

  /* every branch, multiplied by the same factor: the move that lets the tree
     grow or shrink as a whole, which single-branch moves reach only slowly */
  function moveTreeScale(c, r, tune) {
    const n = c.F.n;
    const m = multiplier(r, tune.scaleLambda);
    const old = Float64Array.from(c.lens);
    let dPrior = 0;
    for (let k = 1; k < n; k++) {
      const nu = Math.min(20, Math.max(1e-8, old[k] * m));
      c.lens[k] = nu;
      dPrior += logExp(nu, c.pr.branchMean) - logExp(old[k], c.pr.branchMean);
    }
    let lnLnew = c.lnL;
    if (!c.noLikelihood) lnLnew = c.lik.full(c.lens);
    return {
      logHastings: (n - 1) * Math.log(m),
      dLnL: lnLnew - c.lnL, dPrior,
      accept() {},
      undo() { c.lens.set(old); if (!c.noLikelihood) c.lik.full(c.lens); },
    };
  }

  /* The tree the sampler carries is UNROOTED, written with a trifurcating root:
     the root is a real internal node of degree three, every other internal node
     has two children plus its parent, and the number of branches is exactly
     2n − 3. That matters for more than tidiness — with a bifurcating root the
     two basal branches would be two halves of one real branch, so the prior
     would count 2n − 2 exponentials for 2n − 3 branches and the posterior on
     tree length would be wrong.

     Every topology move works on a CLONE and, if rejected, simply puts the old
     tree back. Cloning thirty-odd nodes costs nothing next to one likelihood. */

  /* The chain's branch lengths live in c.lens, not in the node objects of
     c.tree: the incremental engine works on the flat array and never writes
     back. So a topology move must rebuild the tree FROM the flat state, or it
     silently resets every branch to whatever the last topology move left in the
     node objects — which is what made the sampled branch lengths come out with
     the right mean and far too little spread. */
  function currentTree(c) { return g.Tree.unflatten(c.F, c.lens); }

  /* A stochastic NNI, done on the FLAT tree and scored incrementally.

     Rebuilding the tree and recomputing every partial would cost a full
     likelihood — five times an incremental one — and the NNI is the move the
     sampler makes most. It does not need any of that. An NNI across the branch
     (u, v) swaps a neighbour a of u with a neighbour b of v; in the flattened
     arrays that is four assignments, and the postorder stays valid, because
     every node that was below u still is. Only v, u and the nodes above them
     change, which is exactly what the engine's incremental path recomputes.

     Branch lengths travel with the subtrees they belong to, so the proposal is
     symmetric and the Hastings ratio is one. */
  function moveNNI(c, r) {
    const F = c.F;
    const inner = [];
    for (let v = 1; v < F.n; v++) if (!F.isTip[v]) inner.push(v);
    if (!inner.length) return null;
    const v = inner[Math.floor(r() * inner.length)];
    const u = F.parent[v];
    if (u < 0) return null;
    const sibs = F.kids[u].filter(k => k !== v);
    const kids = F.kids[v];
    if (!sibs.length || !kids.length) return null;
    const a = sibs[Math.floor(r() * sibs.length)];
    const b = kids[Math.floor(r() * kids.length)];

    const oldPost = F.post;
    F.kids[u] = F.kids[u].map(k => (k === a ? b : k));
    F.kids[v] = F.kids[v].map(k => (k === b ? a : k));
    F.parent[a] = v; F.parent[b] = u;
    /* The postorder has to be rebuilt. It survives the swap only when the moved
       subtree already came before its new parent, and whether it did depends on
       the order the children happen to be stored in — so half the time it did
       not, and the engine then combined a child that had not been recomputed
       yet. The likelihood came out ABOVE the maximum, which is how the bug was
       caught. Rebuilding it is a walk over thirty-odd nodes, nothing next to
       one likelihood. */
    F.post = postorder(F);

    const lnLnew = c.noLikelihood ? 0 : c.lik.propose([a, b], c.lens);
    return {
      logHastings: 0, dLnL: lnLnew - c.lnL, dPrior: 0,
      accept() { if (!c.noLikelihood) c.lik.accept(); },
      undo() {
        F.kids[u] = F.kids[u].map(k => (k === b ? a : k));
        F.kids[v] = F.kids[v].map(k => (k === a ? b : k));
        F.parent[a] = u; F.parent[b] = v;
        F.post = oldPost;
        if (!c.noLikelihood) c.lik.reject();
      },
    };
  }
  function postorder(F) {
    const out = [];
    (function walk(k) { F.kids[k].forEach(walk); out.push(k); })(0);
    return out;
  }

  /* An SPR that keeps every branch length. Pruning leaves a node of degree two
     behind, which is spliced out and its two branches merged into one of length
     a + b; regrafting splits the target branch T into uT and (1−u)T with u
     uniform. Forward and reverse proposals have the same density — the same
     number of branches to cut and the same number to graft onto — so the whole
     Hastings ratio is the Jacobian of
         (a, b, T, u) → (a+b, uT, (1−u)T, a/(a+b)),   which is T/(a+b). */
  function moveSPR(c, r) {
    const tree = currentTree(c);
    if (g.Tree.tips(tree).length < 4) return null;
    const edges = g.Tree.allEdges(tree);
    const cut = edges[Math.floor(r() * edges.length)];
    const parent = cut.parent, pruned = cut.edge;
    /* what stays behind has to keep at least three taxa, or there is nothing
       left to graft onto */
    const nPruned = g.Tree.tips(pruned.node).length;
    if (g.Tree.tips(tree).length - nPruned < 3) return null;

    parent.children = parent.children.filter(e => e !== pruned);
    let root = tree, merged = null;
    if (parent !== tree) {
      /* an ordinary internal node: it had two children, now one — splice it */
      const gp = findParent(tree, parent);
      if (!gp) return null;
      const e = gp.children.find(x => x.node === parent);
      merged = (e.len || 0) + (parent.children[0].len || 0);
      e.len = merged;
      e.node = parent.children[0].node;
    } else {
      /* the trifurcating root: it had three children, now two, so it is a node
         of degree two in the unrooted tree and has to disappear as well. One of
         its two children becomes the new root and the other hangs from it. */
      const a = tree.children[0], b = tree.children[1];
      const host = a.node.children.length ? a : (b.node.children.length ? b : null);
      if (!host) return null;
      const other = host === a ? b : a;
      merged = (a.len || 0) + (b.len || 0);
      host.node.children.push({ node: other.node, len: merged });
      root = host.node;
    }

    const rests = g.Tree.allEdges(root);
    if (!rests.length) return null;
    const at = rests[Math.floor(r() * rests.length)];
    const T = Math.max(1e-8, at.edge.len || 1e-8);
    const u = r();
    at.edge.node = {
      children: [
        { node: at.edge.node, len: Math.max(1e-9, u * T) },
        { node: pruned.node, len: pruned.len },
      ],
    };
    at.edge.len = Math.max(1e-9, (1 - u) * T);
    return topologyChanged(c, root, Math.log(T) - Math.log(Math.max(1e-12, merged)));
  }
  function findParent(root, target2) {
    let found = null;
    (function walk(n) { n.children.forEach(ch => { if (ch.node === target2) found = n; else walk(ch.node); }); })(root);
    return found;
  }

  /* what every topology move does afterwards: adopt the new tree, reflatten,
     rebuild the engine and recompute from scratch */
  function topologyChanged(c, newTree, logJ) {
    const old = { tree: c.tree, F: c.F, lens: c.lens, lik: c.lik };
    c.tree = newTree;
    c.F = g.Tree.flatten(newTree);
    c.lens = Float64Array.from(c.F.len, v => Math.max(1e-9, v || 0.01));
    if (!c.noLikelihood) c.lik = g.Like.engine(c.F, c.A, c.M, c.pool);
    const lnLnew = c.noLikelihood ? 0 : c.lik.full(c.lens);
    const dPrior = logBranchPrior(c.lens, c.F.n, c.pr) - logBranchPrior(old.lens, old.F.n, c.pr);
    return {
      logHastings: logJ,
      dLnL: lnLnew - c.lnL, dPrior,
      accept() {},
      undo() {
        c.tree = old.tree; c.F = old.F; c.lens = old.lens; c.lik = old.lik;
        if (!c.noLikelihood && c.lik) c.lik.full(c.lens);
      },
    };
  }

  /* ---- model parameters ---- */
  function moveRate(c, r, tune) {
    if (c.spec.type !== 'dna' || !c.spec.rates) return null;
    const scheme = g.Like.DNA_SCHEMES[c.spec.model] || g.Like.DNA_SCHEMES.GTR;
    const nFree = Math.max.apply(null, scheme);
    if (nFree < 1) return null;
    const i = 1 + Math.floor(r() * nFree);
    const old = c.spec.rates[i];
    const m = multiplier(r, tune.rateLambda);
    const nu = Math.min(500, Math.max(1e-4, old * m));
    c.spec.rates[i] = nu;
    return modelChanged(c, Math.log(nu / old),
      logExp(nu, c.pr.rateMean) - logExp(old, c.pr.rateMean),
      () => { c.spec.rates[i] = old; });
  }
  function moveAlpha(c, r, tune) {
    if (c.spec.alpha == null) return null;
    const old = c.spec.alpha;
    const m = multiplier(r, tune.alphaLambda);
    const nu = Math.min(200, Math.max(0.01, old * m));
    c.spec.alpha = nu;
    return modelChanged(c, Math.log(nu / old),
      logExp(nu, c.pr.alphaMean) - logExp(old, c.pr.alphaMean),
      () => { c.spec.alpha = old; });
  }
  function movePInv(c, r, tune) {
    if (c.spec.pInv == null) return null;
    const old = c.spec.pInv;
    let nu = old + (r() - 0.5) * tune.pInvWindow;
    /* reflect at the boundaries, which keeps the proposal symmetric */
    while (nu < 0 || nu > 0.99) { if (nu < 0) nu = -nu; if (nu > 0.99) nu = 1.98 - nu; }
    c.spec.pInv = nu;
    return modelChanged(c, 0, 0, () => { c.spec.pInv = old; });
  }
  /* a Dirichlet proposal centred on the current frequencies */
  function moveFreqs(c, r, tune) {
    if (!c.spec.freqs || c.spec.freqs.length < 2) return null;
    const old = Array.from(c.spec.freqs);
    const n = old.length, a0 = tune.freqsAlpha;
    const alphaF = old.map(v => Math.max(1e-6, v * a0));
    const draw = alphaF.map(a => gammaDraw(r, a));
    let sum = 0; draw.forEach(v => { sum += v; });
    if (!(sum > 0)) return null;
    const nu = draw.map(v => Math.max(1e-6, v / sum));
    let s2 = 0; nu.forEach(v => { s2 += v; });
    for (let i = 0; i < n; i++) nu[i] /= s2;
    c.spec.freqs = nu;
    /* Hastings: q(old | new) / q(new | old), both Dirichlet densities */
    const alphaB = nu.map(v => Math.max(1e-6, v * a0));
    const logH = dirichletLogPdf(old, alphaB) - dirichletLogPdf(nu, alphaF);
    const a = c.pr.freqsDirichlet;
    let dPrior = 0;
    if (a !== 1) for (let i = 0; i < n; i++) dPrior += (a - 1) * (Math.log(nu[i]) - Math.log(old[i]));
    return modelChanged(c, logH, dPrior, () => { c.spec.freqs = old; });
  }
  function modelChanged(c, logH, dPrior, restore) {
    const oldM = c.M;
    c.M = g.Like.model(c.spec, c.A);
    let lnLnew = c.lnL;
    if (!c.noLikelihood) { c.lik.setModel(c.M); lnLnew = c.lik.full(c.lens); }
    return {
      logHastings: logH, dLnL: lnLnew - c.lnL, dPrior,
      accept() {},
      undo() {
        restore();
        c.M = oldM;
        if (!c.noLikelihood) { c.lik.setModel(c.M); c.lik.full(c.lens); }
      },
    };
  }

  /* gamma variate (Marsaglia & Tsang 2000), for the Dirichlet proposal */
  function gammaDraw(r, shape) {
    if (shape < 1) return gammaDraw(r, shape + 1) * Math.pow(Math.max(1e-300, r()), 1 / shape);
    const d = shape - 1 / 3, cc = 1 / Math.sqrt(9 * d);
    for (;;) {
      let x, v;
      do { x = normalDraw(r); v = 1 + cc * x; } while (v <= 0);
      v = v * v * v;
      const u = r();
      if (u < 1 - 0.0331 * x * x * x * x) return d * v;
      if (Math.log(u) < 0.5 * x * x + d * (1 - v + Math.log(v))) return d * v;
    }
  }
  function normalDraw(r) {
    let u = 0, v = 0;
    while (u === 0) u = r();
    while (v === 0) v = r();
    return Math.sqrt(-2 * Math.log(u)) * Math.cos(2 * Math.PI * v);
  }
  function lgamma(x) {
    const cof = [76.18009172947146, -86.50532032941677, 24.01409824083091,
      -1.231739572450155, 0.1208650973866179e-2, -0.5395239384953e-5];
    let y = x, tmp = x + 5.5;
    tmp -= (x + 0.5) * Math.log(tmp);
    let ser = 1.000000000190015;
    for (let j = 0; j < 6; j++) ser += cof[j] / ++y;
    return -tmp + Math.log(2.5066282746310005 * ser / x);
  }
  function dirichletLogPdf(x, alpha) {
    let s = 0, sa = 0;
    for (let i = 0; i < x.length; i++) { s += (alpha[i] - 1) * Math.log(Math.max(1e-300, x[i])); sa += alpha[i]; }
    let lb = lgamma(sa);
    for (let i = 0; i < alpha.length; i++) lb -= lgamma(alpha[i]);
    return lb + s;
  }

  /* ================================================================
     4 · the sampler
     ================================================================ */
  const MOVES = [
    /* The weights are chosen by what each move COSTS as much as by what it
       does. A branch multiplier and an NNI are incremental — they recompute
       only the path to the root — while a tree rescaling, an SPR and every
       model parameter force a full likelihood, five times dearer. So the two
       cheap moves carry most of the weight and the expensive ones are made
       often enough to mix and no more. */
    { name: 'branch', weight: 40, fn: moveBranch, tunes: 'branchLambda', target: 0.35 },
    { name: 'nni', weight: 12, fn: moveNNI, tunes: null },
    { name: 'treeScale', weight: 2, fn: moveTreeScale, tunes: 'scaleLambda', target: 0.35 },
    { name: 'spr', weight: 4, fn: moveSPR, tunes: null },
    { name: 'rate', weight: 3, fn: moveRate, tunes: 'rateLambda', target: 0.35 },
    { name: 'alpha', weight: 2, fn: moveAlpha, tunes: 'alphaLambda', target: 0.35 },
    { name: 'pInv', weight: 2, fn: movePInv, tunes: 'pInvWindow', target: 0.35 },
    { name: 'freqs', weight: 2, fn: moveFreqs, tunes: 'freqsAlpha', target: 0.35, inverse: true },
  ];

  function defaultTune() {
    return { branchLambda: 1.2, scaleLambda: 0.6, rateLambda: 1.0, alphaLambda: 1.0, pInvWindow: 0.2, freqsAlpha: 300 };
  }

  /* one Metropolis–Hastings step on one chain.
      lets a caller switch a move off (weight 0) or favour another; the
     audit in run() uses it to find which move is responsible when the running
     log-likelihood drifts away from a fresh recomputation. */
  function step(c, r, tune, stats, weights) {
    const w = m => (weights && weights[m.name] != null ? weights[m.name] : m.weight);
    const total = MOVES.reduce((s, m) => s + w(m), 0);
    if (!(total > 0)) return null;
    let u = r() * total, pick = MOVES[0];
    for (let i = 0; i < MOVES.length; i++) { u -= w(MOVES[i]); if (u <= 0) { pick = MOVES[i]; break; } }
    const prop = pick.fn(c, r, tune);
    if (!prop) return null;
    const st = stats[pick.name] || (stats[pick.name] = { tried: 0, accepted: 0 });
    st.tried++;
    const logAlpha = c.beta * (prop.dLnL + prop.dPrior) + prop.logHastings;
    if (logAlpha >= 0 || Math.log(r()) < logAlpha) {
      c.lnL += prop.dLnL;
      c.lnPrior += prop.dPrior;
      prop.accept();
      st.accepted++;
      return pick.name;
    }
    prop.undo();
    return null;
  }

  /* auto-tuning during burn-in, towards a target acceptance rate */
  function tuneStep(tune, stats) {
    MOVES.forEach(m => {
      if (!m.tunes) return;
      const st = stats[m.name];
      if (!st || st.tried < 30) return;
      const rate = st.accepted / st.tried;
      const f = rate > m.target ? 1.1 : 1 / 1.1;
      const dir = m.inverse ? 1 / f : f;
      tune[m.tunes] = Math.min(m.inverse ? 1e5 : 10, Math.max(m.inverse ? 5 : 0.02, tune[m.tunes] * dir));
      st.tried = 0; st.accepted = 0;
    });
  }

  /* ---------------------------------------------------------------
     run: one independent run, with as many heated chains as asked for
     --------------------------------------------------------------- */
  /* opts: { tree, labels, spec, generations, burnin, sampleEvery, chains, deltaT,
             seed, priors, noLikelihood, beta (for stepping stone), onProgress,
             cancelled } */
  function run(A, opts) {
    opts = opts || {};
    const r = g.rng(opts.seed || 1);
    const pr = Object.assign(defaultPriors(opts.spec), opts.priors || {});
    const nChains = Math.max(1, opts.chains || 4);
    const deltaT = opts.deltaT == null ? 0.1 : opts.deltaT;
    const gens = opts.generations || 100000;
    const burnin = opts.burnin == null ? Math.floor(gens * 0.25) : opts.burnin;
    const every = opts.sampleEvery || Math.max(1, Math.round(gens / 2000));
    /* ONE POOL PER CHAIN. The buffers a pooled engine borrows are its working
       memory, and in Metropolis coupling several chains are alive at once: a
       shared pool means every chain computes its likelihood in another chain's
       partials. The symptom was spectacular and easy to miss — the sampled
       log-likelihoods came out four thousand units ABOVE the maximum, which is
       impossible, and the model parameters never moved at all. */
    const nNodes = g.Tree.nodes(opts.tree).length + 4;
    const M0 = g.Like.model(opts.spec, A);

    const chains = [];
    for (let i = 0; i < nChains; i++) {
      const beta = (opts.beta == null ? 1 : opts.beta) / (1 + deltaT * i);
      const pool = opts.noLikelihood ? null : g.Like.makePool(nNodes, A, M0);
      chains.push(makeChain(A, opts.spec, opts.tree, pr, { beta, noLikelihood: opts.noLikelihood, pool }));
    }
    const tune = defaultTune();
    const stats = {};
    const swaps = { tried: 0, accepted: 0 };
    const audit = { worst: null, byMove: {} };

    const samples = [];                 // one row per sample of the cold chain
    const trees = [];                   // its topology, as Newick
    const splitCounts = new Map();
    let nSampled = 0;
    const t0 = (typeof performance !== 'undefined' ? performance.now() : Date.now());

    for (let gen = 1; gen <= gens; gen++) {
      let lastMove = null;
      for (let i = 0; i < chains.length; i++) {
        const took = step(chains[i], r, tune, i === 0 ? stats : {}, opts.weights);
        if (i === 0) lastMove = took;
      }
      /* The audit: every `audit` generations the running log-likelihood is
         compared with one computed from scratch. They must agree to the last
         bit; anything else means a move is leaving the engine in a state that
         does not match what it reported, and the culprit is the move that ran
         just before the discrepancy appeared. */
      if (opts.audit && !opts.noLikelihood && gen % opts.audit === 0) {
        const c0 = chains[0];
        const fresh = g.Like.engine(c0.F, c0.A, c0.M).full(c0.lens);
        const d = Math.abs(fresh - c0.lnL);
        if (!audit.worst || d > audit.worst.d) audit.worst = { gen, d, move: lastMove, running: c0.lnL, fresh };
        audit.byMove[lastMove || 'ninguno'] = Math.max(audit.byMove[lastMove || 'ninguno'] || 0, d);
        if (opts.auditFix !== false) c0.lnL = fresh;
      }
      /* a swap between two neighbouring temperatures */
      if (chains.length > 1) {
        const i = Math.floor(r() * (chains.length - 1));
        const a = chains[i], b = chains[i + 1];
        swaps.tried++;
        const logR = (a.beta - b.beta) * ((b.lnL + b.lnPrior) - (a.lnL + a.lnPrior));
        if (logR >= 0 || Math.log(r()) < logR) {
          const tb = a.beta; a.beta = b.beta; b.beta = tb;
          chains[i] = b; chains[i + 1] = a;
          swaps.accepted++;
        }
      }
      if (gen <= burnin && gen % 100 === 0) tuneStep(tune, stats);
      if (gen > burnin && gen % every === 0) {
        const c = chains[0];
        const row = {
          gen, lnL: c.lnL, lnPrior: c.lnPrior, lnPosterior: c.lnL + c.lnPrior,
          treeLength: totalLen(c.lens, c.F.n),
        };
        if (c.spec.alpha != null) row.alpha = c.spec.alpha;
        if (c.spec.pInv != null) row.pInv = c.spec.pInv;
        if (c.spec.type === 'dna' && c.spec.rates) {
          const scheme = g.Like.DNA_SCHEMES[c.spec.model] || g.Like.DNA_SCHEMES.GTR;
          const nFree = Math.max.apply(null, scheme);
          for (let k = 1; k <= nFree; k++) row['r' + k] = c.spec.rates[k];
        }
        if (c.spec.freqs) c.spec.freqs.forEach((v, i2) => { row['pi' + 'ACGT'[i2]] = v; });
        samples.push(row);
        const t = g.Tree.unflatten(c.F, c.lens);
        trees.push(g.Tree.writeNewick(t, { labels: opts.labels, support: false }));
        g.Tree.splits(t, opts.labels.length).forEach((node, key) => {
          splitCounts.set(key, (splitCounts.get(key) || 0) + 1);
        });
        nSampled++;
      }
      if (opts.onProgress && gen % Math.max(1, Math.floor(gens / 100)) === 0) {
        opts.onProgress(gen, gens, chains[0].lnL);
      }
      if (opts.cancelled && opts.cancelled()) break;
    }

    const splitFreq = new Map();
    splitCounts.forEach((v, k) => splitFreq.set(k, v / Math.max(1, nSampled)));
    const acceptance = {};
    Object.keys(stats).forEach(k => { acceptance[k] = stats[k].tried ? stats[k].accepted / stats[k].tried : 0; });
    return {
      samples, trees, splitFreq, acceptance, tune, audit,
      swapRate: swaps.tried ? swaps.accepted / swaps.tried : 0,
      generations: gens, burnin, sampleEvery: every, chains: nChains,
      ms: (typeof performance !== 'undefined' ? performance.now() : Date.now()) - t0,
    };
  }
  function totalLen(lens, n) { let s = 0; for (let k = 1; k < n; k++) s += lens[k]; return s; }

  /* ================================================================
     5 · diagnostics
     ================================================================ */
  /* Effective sample size: n divided by the integrated autocorrelation time,
     summed while the autocorrelation stays positive (the initial positive
     sequence of Geyer 1992, which is what Tracer reports). */
  function ess(x) {
    const n = x.length;
    if (n < 10) return n;
    let mean = 0;
    for (let i = 0; i < n; i++) mean += x[i];
    mean /= n;
    let var0 = 0;
    for (let i = 0; i < n; i++) var0 += (x[i] - mean) * (x[i] - mean);
    var0 /= n;
    if (var0 <= 0) return n;
    const maxLag = Math.min(n - 1, Math.floor(n / 3));
    let sum = 0;
    for (let lag = 1; lag <= maxLag; lag++) {
      let s = 0;
      for (let i = 0; i < n - lag; i++) s += (x[i] - mean) * (x[i + lag] - mean);
      const rho = s / (n * var0);
      if (rho < 0.05 && lag > 3) break;
      sum += rho;
    }
    const act = 1 + 2 * sum;
    return act <= 0 ? n : Math.min(n, n / act);
  }

  /* the highest posterior density interval: the shortest interval that holds
     the requested mass */
  function hpd(x, prob) {
    const s = Array.from(x).sort((a, b) => a - b);
    const n = s.length;
    const p = prob == null ? 0.95 : prob;
    const m = Math.max(1, Math.floor(n * p));
    let best = 0, width = Infinity;
    for (let i = 0; i + m < n; i++) {
      const w = s[i + m] - s[i];
      if (w < width) { width = w; best = i; }
    }
    return { lower: s[best], upper: s[best + m] };
  }
  function quantile(x, p) {
    const s = Array.from(x).sort((a, b) => a - b);
    const i = (s.length - 1) * p;
    const lo = Math.floor(i), hi = Math.ceil(i);
    return s[lo] + (s[hi] - s[lo]) * (i - lo);
  }
  function mean(x) { let s = 0; for (let i = 0; i < x.length; i++) s += x[i]; return s / x.length; }
  function sd(x) {
    const m = mean(x);
    let s = 0;
    for (let i = 0; i < x.length; i++) s += (x[i] - m) * (x[i] - m);
    return Math.sqrt(s / Math.max(1, x.length - 1));
  }

  /* Gelman & Rubin's potential scale reduction factor, over several runs of the
     same length. Values near 1 mean the runs cannot be told apart. */
  function psrf(chainsX) {
    const m = chainsX.length;
    if (m < 2) return null;
    const n = Math.min.apply(null, chainsX.map(c => c.length));
    if (n < 10) return null;
    const means = chainsX.map(c => mean(c.slice(0, n)));
    const vars = chainsX.map(c => { const v = c.slice(0, n); const mu = mean(v); let s = 0; v.forEach(z => { s += (z - mu) * (z - mu); }); return s / (n - 1); });
    const W = mean(vars);
    const grand = mean(means);
    let B = 0; means.forEach(mu => { B += (mu - grand) * (mu - grand); });
    B = B * n / (m - 1);
    if (W <= 0) return 1;
    const varHat = (n - 1) / n * W + B / n;
    return Math.sqrt(varHat / W);
  }

  /* The average standard deviation of split frequencies: the number MrBayes
     prints while it runs. Splits seen in fewer than `minFreq` of the samples of
     every run are left out, as MrBayes does, because rare splits are noise. */
  function asdsf(freqMaps, minFreq) {
    const cut = minFreq == null ? 0.1 : minFreq;
    const keys = new Set();
    freqMaps.forEach(f => f.forEach((v, k) => { if (v >= cut) keys.add(k); }));
    if (!keys.size) return 0;
    let sum = 0, n = 0;
    keys.forEach(k => {
      const vals = freqMaps.map(f => f.get(k) || 0);
      const mu = mean(vals);
      let s = 0;
      vals.forEach(v => { s += (v - mu) * (v - mu); });
      sum += Math.sqrt(s / vals.length);
      n++;
    });
    return n ? sum / n : 0;
  }

  /* ================================================================
     6 · summarising the posterior
     ================================================================ */
  /* the majority-rule consensus with posterior probabilities, plus the tree
     that was sampled most often (the maximum a posteriori topology) and the
     95 % credible set of topologies */
  function summarise(res, labels) {
    const nTaxa = labels.length;
    const counts = new Map();
    res.trees.forEach(nw => {
      const key = topologyKey(nw, labels, nTaxa);
      const e = counts.get(key);
      if (e) e.n++; else counts.set(key, { n: 1, newick: nw });
    });
    const sorted = [...counts.values()].sort((a, b) => b.n - a.n);
    const total = res.trees.length;
    let acc = 0;
    const credible = [];
    for (const t of sorted) {
      credible.push({ newick: t.newick, p: t.n / total });
      acc += t.n / total;
      if (acc >= 0.95) break;
    }
    /* mean branch lengths of every split, over the samples that contain it */
    const lenSum = new Map(), lenN = new Map();
    res.trees.forEach(nw => {
      const t = g.Tree.parseNewick(nw, labels);
      g.Tree.splits(t, nTaxa).forEach((node, key) => {
        const parent = findParent(t, node);
        if (!parent) return;
        const e = parent.children.find(x => x.node === node);
        if (!e) return;
        lenSum.set(key, (lenSum.get(key) || 0) + (e.len || 0));
        lenN.set(key, (lenN.get(key) || 0) + 1);
      });
    });
    return {
      mapTopology: sorted.length ? sorted[0].newick : null,
      mapProbability: sorted.length ? sorted[0].n / total : 0,
      nTopologies: sorted.length,
      credibleSet: credible,
      splitFreq: res.splitFreq,
      meanLength: new Map([...lenSum.keys()].map(k => [k, lenSum.get(k) / lenN.get(k)])),
    };
  }
  /* two Newick strings describe the same topology when their split sets match */
  function topologyKey(newick, labels, nTaxa) {
    const t = g.Tree.parseNewick(newick, labels);
    return [...g.Tree.splits(t, nTaxa).keys()].sort().join('|');
  }

  /* ================================================================
     7 · the marginal likelihood, by stepping stone
     ================================================================ */
  /* Xie et al. (2011): a ladder of powers beta_k between 0 and 1; at each one a
     chain samples from the tempered posterior, and the ratio of consecutive
     normalising constants is estimated from the sampled likelihoods. Summing
     the logs gives log P(data | model), which is what a Bayes factor compares.
     The powers follow a beta(0.3, 1) spacing, which puts more of them near the
     prior where the integrand changes fastest. */
  function steppingStone(A, opts) {
    const K = opts.steps || 16;
    const betas = [];
    for (let k = 0; k < K; k++) betas.push(Math.pow(k / (K - 1), 1 / 0.3));
    let logML = 0;
    const perStep = [];
    for (let k = K - 1; k > 0; k--) {
      const bHigh = betas[k], bLow = betas[k - 1];
      const res = run(A, Object.assign({}, opts, {
        beta: bLow, chains: 1, seed: (opts.seed || 1) + k * 101,
        generations: opts.stepGenerations || 4000,
        burnin: Math.floor((opts.stepGenerations || 4000) * 0.3),
        sampleEvery: opts.stepSampleEvery || 10,
      }));
      const d = bHigh - bLow;
      const ls = res.samples.map(s => s.lnL);
      const maxL = Math.max.apply(null, ls);
      let sum = 0;
      ls.forEach(v => { sum += Math.exp(d * (v - maxL)); });
      const term = d * maxL + Math.log(sum / ls.length);
      logML += term;
      perStep.push({ betaLow: bLow, betaHigh: bHigh, term, n: ls.length });
      if (opts.onStep) opts.onStep(K - k, K - 1, logML);
      if (opts.cancelled && opts.cancelled()) break;
    }
    return { logMarginalLikelihood: logML, steps: perStep, betas };
  }

  /* a log file in the format Tracer reads */
  function tracerLog(samples) {
    if (!samples.length) return '';
    const keys = Object.keys(samples[0]).filter(k => k !== 'gen');
    const head = ['state'].concat(keys).join('\t');
    const rows = samples.map(s => [s.gen].concat(keys.map(k => (s[k] == null ? 'NA' : Number(s[k]).toPrecision(8)))).join('\t'));
    return head + '\n' + rows.join('\n') + '\n';
  }

  Object.assign(Mcmc, {
    run, steppingStone, summarise, tracerLog,
    ess, hpd, quantile, mean, sd, psrf, asdsf,
    defaultPriors, defaultTune, logBranchPrior, logModelPrior,
    gammaDraw, normalDraw, dirichletLogPdf, lgamma, topologyKey, MOVES,
  });
  g.Mcmc = Mcmc;
}
McmcCore(typeof window !== 'undefined' ? window : self);
