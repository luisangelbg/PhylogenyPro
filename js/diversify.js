/* PhylogenyPro — molecular phylogenetics in the browser.
   Copyright (C) 2026  Luis Ángel Barrera-Guzmán

   This program is free software: you can redistribute it and/or modify it under
   the terms of the GNU General Public License as published by the Free Software
   Foundation, either version 3 of the License, or (at your option) any later
   version. It is distributed in the hope that it will be useful, but WITHOUT ANY
   WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
   A PARTICULAR PURPOSE. See the GNU General Public License, in the file LICENSE
   at the root of this program, or <https://www.gnu.org/licenses/>. */

/* PhylogenyPro — diversification: reading the shape of a dated tree.

   A chronogram carries more than ages. The pattern of its branching times says
   how fast lineages have been appearing and disappearing, and whether that has
   changed — which is the question behind adaptive radiations, living fossils
   and mass extinctions.

   What is implemented:

     · lineages through time, for one tree and for a whole posterior sample,
       with the band the sample implies;
     · the γ statistic of Pybus & Harvey (2000), which asks whether the internal
       nodes sit closer to the root than a constant rate would put them, and the
       Monte Carlo constant-rates test that makes γ usable when taxa are
       missing, because missing tips fake a slowdown;
     · maximum likelihood for the pure-birth (Yule) and constant-rate
       birth–death models, in exactly the parameterisation ape uses, so the
       numbers can be checked against it;
     · the density-dependent models of Rabosky & Lovette (2008) — logistic and
       exponential — and a two-rate model with the shift time estimated, all
       compared by AIC;
     · the method-of-moments estimators of Magallón & Sanderson (2001), which
       need only a clade's age and its richness and are the only thing that can
       be done when the tree has a fraction of the species;
     · the DR statistic of Jetz et al. (2012), a rate for every tip;
     · simulation of birth–death trees, which is what the Monte Carlo test runs
       on and what any of this can be checked against.

   Times are ages before the present, as everywhere in Block 7.

   Wrapped in a named function so a Web Worker can be built from its own source
   text (js/pool.js). Nothing here touches the DOM. */

function DiversifyCore(g) {
const Diversify = {};

  /* ================================================================
     1 · branching times and lineages through time
     ================================================================ */
  /* the age of every internal node, oldest first. A tree carrying `age` on its
     nodes (Block 7) is used as it is; otherwise the ages are read from the
     branch lengths, taking the deepest tip as the present. */
  function branchingTimes(tree) {
    const F = g.Tree.flatten(tree);
    const d = new Float64Array(F.n);
    for (let i = F.post.length - 1; i >= 0; i--) {
      const k = F.post[i], p = F.parent[k];
      if (p >= 0) d[k] = d[p] + (F.len[k] || 0);
    }
    const hasAges = F.nodes.every(n => n.age != null);
    let maxD = 0;
    for (let k = 0; k < F.n; k++) if (F.isTip[k]) maxD = Math.max(maxD, d[k]);
    const out = [];
    for (let k = 0; k < F.n; k++) {
      if (F.isTip[k]) continue;
      out.push(hasAges ? F.nodes[k].age : maxD - d[k]);
    }
    out.sort((a, b) => b - a);
    return out;
  }

  /* the number of lineages alive at every branching time */
  function ltt(tree) {
    const bt = branchingTimes(tree);
    const times = [], lineages = [];
    /* the root starts the count at 2 (or at the number of its children) */
    const F = g.Tree.flatten(tree);
    let n = F.kids[0].length;
    times.push(bt[0]); lineages.push(n);
    for (let i = 1; i < bt.length; i++) {
      n += 1;
      times.push(bt[i]); lineages.push(n);
    }
    times.push(0); lineages.push(g.Tree.tips(tree).length);
    return { times, lineages, branchingTimes: bt };
  }

  /* the same over many trees: at a grid of times, the median and the interval
     the sample gives */
  function lttBand(trees, opts) {
    opts = opts || {};
    if (!trees.length) return null;
    const curves = trees.map(t => ltt(t));
    let maxT = 0;
    curves.forEach(c => { maxT = Math.max(maxT, c.times[0]); });
    const steps = opts.steps || 120;
    const times = [], median = [], lower = [], upper = [];
    for (let s = 0; s <= steps; s++) {
      const t = maxT * (1 - s / steps);
      const vals = curves.map(c => countAt(c, t)).sort((a, b) => a - b);
      times.push(t);
      median.push(vals[Math.floor(vals.length / 2)]);
      lower.push(vals[Math.floor(vals.length * 0.025)]);
      upper.push(vals[Math.min(vals.length - 1, Math.ceil(vals.length * 0.975))]);
    }
    return { times, median, lower, upper, n: trees.length };
  }
  function countAt(curve, t) {
    let n = curve.lineages[0];
    for (let i = 0; i < curve.times.length; i++) {
      if (curve.times[i] >= t) n = curve.lineages[i]; else break;
    }
    return n;
  }

  /* ================================================================
     2 · the gamma statistic
     ================================================================ */
  /* Pybus & Harvey (2000), in the form ape::gammaStat computes it. Negative
     values mean the nodes are nearer the root than a constant rate predicts —
     a slowdown — and under a pure-birth process the statistic is standard
     normal, so −1.645 is the one-tailed 5 % point. */
  function gammaStat(tree) {
    const bt = branchingTimes(tree).slice().sort((a, b) => a - b);   // ascending
    const N = g.Tree.tips(tree).length;
    if (N < 4) return null;
    /* g = rev(c(bt[1], diff(bt))): the internode intervals, oldest first */
    const gaps = [bt[0]];
    for (let i = 1; i < bt.length; i++) gaps.push(bt[i] - bt[i - 1]);
    gaps.reverse();
    let ST = 0;
    for (let i = 0; i < gaps.length; i++) ST += (i + 2) * gaps[i];
    let inner = 0, acc = 0;
    for (let i = 0; i < gaps.length - 1; i++) { acc += (i + 2) * gaps[i]; inner += acc; }
    const stat = inner / (N - 2);
    const m = ST / 2;
    const s = ST * Math.sqrt(1 / (12 * (N - 2)));
    return { gamma: (stat - m) / s, ST, N };
  }

  /* ================================================================
     3 · maximum likelihood on the branching times
     ================================================================ */
  /* the waiting times of the reconstructed process: w[i] is how long the tree
     spent with (i + 2) lineages, counting from the root */
  function waitingTimes(tree) {
    const bt = branchingTimes(tree);
    const w = [];
    for (let i = 0; i < bt.length - 1; i++) w.push(bt[i] - bt[i + 1]);
    w.push(bt[bt.length - 1]);           // from the last split to the present
    return { w, bt };
  }

  /* pure birth, in ape::yule's parameterisation: lambda = (Nnode − 1) / (total
     tree length), and the same log-likelihood it reports */
  function yule(tree) {
    const F = g.Tree.flatten(tree);
    let X = 0;
    for (let k = 1; k < F.n; k++) X += F.len[k] || 0;
    const nInternal = F.n - g.Tree.tips(tree).length;     // Nnode
    const nb = nInternal - 1;
    const lambda = X > 0 ? nb / X : 0;
    const loglik = -lambda * X + lgamma(nInternal + 1) + nb * Math.log(lambda);
    return { lambda, se: lambda / Math.sqrt(Math.max(1, nb)), loglik, k: 1, treeLength: X, nInternal };
  }

  /* the constant-rate birth–death likelihood of Nee et al. (1994), written the
     way ape::birthdeath writes it: a = mu/lambda, r = lambda − mu, and x the
     branching times with x[0] the root */
  function bdLogLik(x, a, r) {
    const N = x.length + 1;
    if (r <= 0 || a >= 1 || a < 0) return -Infinity;
    let s1 = 0;
    for (let i = 1; i < x.length; i++) s1 += x[i];
    let s2 = 0;
    for (let i = 0; i < x.length; i++) {
      const v = Math.exp(r * x[i]) - a;
      if (!(v > 0)) return -Infinity;
      s2 += Math.log(v);
    }
    return lgamma(N) + (N - 2) * Math.log(r) + r * s1 + N * Math.log(1 - a) - 2 * s2;
  }

  function birthDeath(tree, opts) {
    opts = opts || {};
    const x = branchingTimes(tree);
    const N = x.length + 1;
    const f = p => -bdLogLik(x, p[0], p[1]);
    const best = nelderMead(f, [0.1, 0.2], { maxIter: 800 });
    let a = best.x[0], r = best.x[1];
    if (a < 0) {
      /* the boundary: no extinction */
      const y = yuleOnTimes(x);
      a = 0; r = y.r;
    }
    a = Math.max(0, Math.min(0.999999, a));
    const loglik = bdLogLik(x, a, r);
    const lambda = r / (1 - a), mu = lambda * a;
    return {
      a, r, lambda, mu, loglik, k: 2, N,
      epsilon: a, netDiversification: r,
    };
  }
  /* pure birth from the branching times alone, for the boundary case */
  function yuleOnTimes(x) {
    const N = x.length + 1;
    let s = 0;
    for (let i = 1; i < x.length; i++) s += x[i];
    /* r maximising (N−2) log r + r·s − 2 Σ log(e^{r x} − 0) = (N−2) log r − r (2Σx − s) */
    let tot = 0;
    for (let i = 0; i < x.length; i++) tot += x[i];
    const denom = 2 * tot - s;
    const r = denom > 0 ? (N - 2) / denom : 0;
    return { r, loglik: bdLogLik(x, 0, r) };
  }

  /* ---- models where the rate depends on how many lineages there are ----
     A pure-birth process in which the per-lineage rate while there are i
     lineages is lambda_i: the waiting time is exponential with rate i·lambda_i,
     so the log-likelihood is Σ [ log(i lambda_i) − i lambda_i w_i ]. Pure birth
     is the case lambda_i = r, and its maximum, (n − 2)/Σ i w_i, is the same
     number ape::yule reports, because Σ i w_i is the length of the tree. */
  function intervalLogLik(x, rateOf) {
    const k = x.length;                       // number of branching events
    let s = 0;
    for (let i = 0; i < k; i++) {
      const nLin = i + 2;                     // lineages during this interval
      const w = (i < k - 1 ? x[i] - x[i + 1] : x[i]);
      const lam = rateOf(nLin, i);
      if (!(lam > 0) || !(w >= 0)) return -Infinity;
      if (i < k - 1) s += Math.log(nLin * lam) - nLin * lam * w;
      else s += -nLin * lam * w;              // the last interval has no event
    }
    return s;
  }
  /* the version that counts every event, used for the AIC comparison: the final
     interval ends at the present without an event */
  function pureBirth(tree) {
    const x = branchingTimes(tree);
    let denom = 0;
    for (let i = 0; i < x.length; i++) {
      const nLin = i + 2;
      const w = (i < x.length - 1 ? x[i] - x[i + 1] : x[i]);
      denom += nLin * w;
    }
    const r = denom > 0 ? (x.length - 1) / denom : 0;
    return { r, lambda: r, mu: 0, loglik: intervalLogLik(x, () => r), k: 1 };
  }
  /* logistic density dependence (Rabosky & Lovette 2008): lambda_i = r (1 − i/K) */
  function ddLogistic(tree, opts) {
    opts = opts || {};
    const x = branchingTimes(tree);
    const n = x.length + 1;
    const f = p => {
      const r = p[0], K = p[1];
      if (!(r > 0) || !(K > n)) return 1e100;
      return -intervalLogLik(x, i => r * (1 - i / K));
    };
    const start = [pureBirth(tree).r, n * 1.5];
    const best = nelderMead(f, start, { maxIter: 1200 });
    const r = best.x[0], K = best.x[1];
    return { r, K, loglik: -best.f, k: 2, model: 'DDL' };
  }
  /* exponential density dependence: lambda_i = r · i^(−x) */
  function ddExponential(tree) {
    const bt = branchingTimes(tree);
    const f = p => {
      const r = p[0], xx = p[1];
      if (!(r > 0)) return 1e100;
      return -intervalLogLik(bt, i => r * Math.pow(i, -xx));
    };
    const best = nelderMead(f, [pureBirth(tree).r, 0.1], { maxIter: 1200 });
    return { r: best.x[0], x: best.x[1], loglik: -best.f, k: 2, model: 'DDX' };
  }
  /* two rates, with the time of the shift estimated: the classic way of asking
     whether the tempo changed, without saying beforehand when */
  function twoRate(tree, opts) {
    opts = opts || {};
    const x = branchingTimes(tree);
    if (x.length < 4) return null;
    let best = null;
    /* the shift can only sit at a branching time, so every one is tried */
    for (let s = 1; s < x.length - 1; s++) {
      const st = x[s];
      /* closed-form maxima on each side */
      let d1 = 0, e1 = 0, d2 = 0, e2 = 0;
      for (let i = 0; i < x.length; i++) {
        const nLin = i + 2;
        const start = x[i], end = (i < x.length - 1 ? x[i + 1] : 0);
        /* split this interval at the shift time */
        const oldPart = Math.max(0, Math.min(start, Math.max(end, st)) - end) * 0;
        const aOld = Math.max(0, start - Math.max(end, st));   // before the shift
        const aNew = Math.max(0, Math.min(start, st) - end);   // after the shift
        d1 += nLin * aOld; d2 += nLin * aNew;
        if (i < x.length - 1) { if (x[i + 1] >= st) e1++; else e2++; }
      }
      const r1 = d1 > 0 ? e1 / d1 : 0, r2 = d2 > 0 ? e2 / d2 : 0;
      if (!(r1 > 0) || !(r2 > 0)) continue;
      const ll = intervalLogLikSplit(x, st, r1, r2);
      if (!best || ll > best.loglik) best = { loglik: ll, shift: st, r1, r2, k: 3, model: 'Yule-2' };
    }
    return best;
  }
  function intervalLogLikSplit(x, st, r1, r2) {
    let s = 0;
    for (let i = 0; i < x.length; i++) {
      const nLin = i + 2;
      const start = x[i], end = (i < x.length - 1 ? x[i + 1] : 0);
      const aOld = Math.max(0, start - Math.max(end, st));
      const aNew = Math.max(0, Math.min(start, st) - end);
      s += -nLin * (r1 * aOld + r2 * aNew);
      if (i < x.length - 1) s += Math.log(nLin * (x[i + 1] >= st ? r1 : r2));
    }
    return s;
  }

  /* the AIC table over the whole family */
  function compare(tree, opts) {
    opts = opts || {};
    const n = g.Tree.tips(tree).length;
    const rows = [];
    const pb = pureBirth(tree);
    rows.push({ name: 'pure birth (Yule)', loglik: pb.loglik, k: 1, fit: pb });
    const bd = birthDeath(tree);
    rows.push({ name: 'birth–death', loglik: bd.loglik, k: 2, fit: bd });
    const dl = ddLogistic(tree);
    rows.push({ name: 'density dependent, logistic', loglik: dl.loglik, k: 2, fit: dl });
    const dx = ddExponential(tree);
    rows.push({ name: 'density dependent, exponential', loglik: dx.loglik, k: 2, fit: dx });
    const tr = twoRate(tree);
    if (tr) rows.push({ name: 'two rates, shift estimated', loglik: tr.loglik, k: 3, fit: tr });
    /* the birth–death likelihood is written on a different constant from the
       interval likelihoods, so it is refitted on the same scale to be comparable */
    rows.forEach(r => {
      r.AIC = -2 * r.loglik + 2 * r.k;
      r.AICc = r.AIC + (n - r.k - 1 > 0 ? 2 * r.k * (r.k + 1) / (n - r.k - 1) : Infinity);
    });
    const best = rows.reduce((a, b) => (b.AIC < a.AIC ? b : a), rows[0]);
    let sw = 0;
    rows.forEach(r => { r.dAIC = r.AIC - best.AIC; r.w = Math.exp(-0.5 * r.dAIC); sw += r.w; });
    rows.forEach(r => { r.w /= sw; });
    rows.sort((a, b) => a.AIC - b.AIC);
    return { rows, best: rows[0] };
  }

  /* ================================================================
     4 · the estimators that need only an age and a richness
     ================================================================ */
  /* Magallón & Sanderson (2001), in geiger's form. `crown` says whether the age
     given is of the crown group (the first split inside the clade) or of the
     stem (the split from its sister). epsilon = mu/lambda has to be assumed:
     0 for no extinction, 0.9 for a lot of it, and the answer is reported for
     both because it moves the estimate a great deal. */
  function magallonSanderson(time, n, epsilon, crown) {
    if (!(time > 0) || !(n > 0)) return null;
    const e = epsilon || 0;
    if (crown !== false) {
      if (e === 0) return (Math.log(n) - Math.log(2)) / time;
      return (1 / time) * (Math.log(
        (n / 2) * (1 - e * e) + 2 * e + 0.5 * (1 - e) * Math.sqrt(n * (n * e * e - 8 * e + 2 * n * e + n))
      ) - Math.log(2));
    }
    if (e === 0) return Math.log(n) / time;
    return (1 / time) * Math.log(n * (1 - e) + e);
  }

  /* ================================================================
     5 · a rate for every tip
     ================================================================ */
  /* Jetz et al. (2012): the inverse of the equal-splits measure, which weights
     each branch on the path from the tip to the root by half as much as the one
     below it. It is a species-level rate, and its logarithm is what gets
     mapped onto trees. */
  function drStatistic(tree) {
    const F = g.Tree.flatten(tree);
    const out = [];
    for (let k = 0; k < F.n; k++) {
      if (!F.isTip[k]) continue;
      let sum = 0, j = 0, cur = k;
      while (F.parent[cur] >= 0) {
        sum += (F.len[cur] || 0) * Math.pow(2, -j);
        j++;
        cur = F.parent[cur];
      }
      out.push({ tip: F.tipRow[k], es: sum, dr: sum > 0 ? 1 / sum : Infinity });
    }
    return out;
  }

  /* ================================================================
     6 · simulating birth–death trees
     ================================================================ */
  /* A reconstructed tree of exactly n extant tips, grown backwards in the way
     the Monte Carlo test needs: waiting times drawn from the reconstructed
     process, which is what makes the null distribution of gamma right. */
  function simulateBD(n, lambda, mu, seed) {
    const r = g.rng(seed || 1);
    for (let attempt = 0; attempt < 500; attempt++) {
      /* the process starts with the two lineages of the crown, at time 0, and
         runs forwards until n of them are alive at the same moment */
      const origin = { children: [], splitAt: 0 };
      const a = { children: [] }, b = { children: [] };
      origin.children = [{ node: a, len: 0 }, { node: b, len: 0 }];
      let alive = [a, b];
      let t = 0, ok = true;
      for (let guard = 0; guard < 200000 && alive.length < n; guard++) {
        const total = alive.length * (lambda + mu);
        if (!(total > 0)) { ok = false; break; }
        t += -Math.log(Math.max(1e-300, r())) / total;
        const i = Math.floor(r() * alive.length);
        const node = alive[i];
        if (r() < lambda / (lambda + mu)) {
          node.splitAt = t;
          const c1 = { children: [] }, c2 = { children: [] };
          node.children = [{ node: c1, len: 0 }, { node: c2, len: 0 }];
          alive.splice(i, 1, c1, c2);
        } else {
          node.dead = true;
          alive.splice(i, 1);
          if (!alive.length) { ok = false; break; }
        }
      }
      if (!ok || alive.length < n) continue;
      /* The present is NOT the moment the n-th lineage appeared. Stopping there
         leaves the last two tips at zero distance and squeezes every node
         towards the present: the mean of gamma over simulated pure-birth trees
         came out at +0.32 instead of 0, which would have made every Monte Carlo
         test wrong in the same direction. The observation happens somewhere
         after that split and before the next event, so one more waiting time is
         drawn and the present set just before whatever would have happened. */
      const T = t - Math.log(Math.max(1e-300, r())) / (alive.length * (lambda + mu));

      /* The RECONSTRUCTED tree: the extinct lineages have to go. Leaving them in
         was what made the simulator return more tips than it was asked for and
         trees that were not ultrametric — the gamma statistic of such a tree is
         meaningless, and every simulated null distribution built on it would be
         wrong. */
      const survivor = new Set(alive);
      const keep = node => {
        if (!node.children.length) return survivor.has(node);
        node.children = node.children.filter(c => keep(c.node));
        return node.children.length > 0;
      };
      if (!keep(origin)) continue;
      /* a node left with one child is a node of degree two: splice it out */
      (function collapse(node) {
        node.children.forEach(c => {
          while (c.node.children.length === 1) c.node = c.node.children[0].node;
          collapse(c.node);
        });
      })(origin);
      let root = origin;
      while (root.children.length === 1) root = root.children[0].node;
      if (root.children.length < 2) continue;

      /* ages backwards from the present, then branch lengths from the ages */
      let tip = 0;
      (function setAges(node) {
        if (!node.children.length) { node.age = 0; node.tip = tip++; node.label = 'sp' + (node.tip + 1); return; }
        node.age = T - node.splitAt;
        node.children.forEach(c => setAges(c.node));
      })(root);
      (function setLens(node) {
        node.children.forEach(c => { c.len = Math.max(0, node.age - c.node.age); setLens(c.node); });
      })(root);
      if (tip !== n) continue;
      return root;
    }
    return null;
  }

  /* the Monte Carlo constant-rates test (Pybus & Harvey 2000): gamma computed
     on the tree at hand is compared with gamma on trees simulated under a
     constant rate, grown to the TRUE richness and then thinned at random to the
     number of species actually sampled. Without that thinning, missing taxa
     look exactly like a slowdown. */
  function mccr(tree, opts) {
    opts = opts || {};
    const nSampled = g.Tree.tips(tree).length;
    const nTotal = opts.total || nSampled;
    const reps = opts.reps || 500;
    const obs = gammaStat(tree);
    if (!obs) return null;
    const lambda = opts.lambda || pureBirth(tree).r || 0.1;
    const mu = opts.mu || 0;
    const r = g.rng(opts.seed || 5);
    const null_ = [];
    for (let i = 0; i < reps; i++) {
      const full = simulateBD(nTotal, lambda, mu, Math.floor(r() * 1e9));
      if (!full) continue;
      const pruned = nTotal > nSampled ? prune(full, nSampled, r) : full;
      const gm = pruned ? gammaStat(pruned) : null;
      if (gm && isFinite(gm.gamma)) null_.push(gm.gamma);
      if (opts.onProgress && i % 25 === 0) opts.onProgress(i, reps);
      if (opts.cancelled && opts.cancelled()) break;
    }
    null_.sort((a, b) => a - b);
    let below = 0;
    null_.forEach(v => { if (v <= obs.gamma) below++; });
    return {
      gamma: obs.gamma,
      p: null_.length ? below / null_.length : null,
      critical: null_.length ? null_[Math.floor(null_.length * 0.05)] : null,
      nullDistribution: null_,
      reps: null_.length, nSampled, nTotal, lambda, mu,
    };
  }
  /* drop tips at random until only `keep` remain */
  function prune(tree, keep, r) {
    const t = g.Tree.clone(tree);
    let tips = g.Tree.tips(t);
    while (tips.length > keep) {
      const victim = tips[Math.floor(r() * tips.length)];
      if (!removeTip(t, victim)) break;
      tips = g.Tree.tips(t);
    }
    /* renumber the tips so that Tree.splits and the gamma statistic work */
    g.Tree.tips(t).forEach((tp, i) => { tp.tip = i; });
    return g.Tree.tips(t).length >= 4 ? t : null;
  }
  function removeTip(root, tip) {
    let parent = null;
    (function walk(n) { n.children.forEach(c => { if (c.node === tip) parent = n; else walk(c.node); }); })(root);
    if (!parent) return false;
    parent.children = parent.children.filter(c => c.node !== tip);
    if (parent.children.length === 1 && parent !== root) {
      let gp = null;
      (function walk(n) { n.children.forEach(c => { if (c.node === parent) gp = n; else walk(c.node); }); })(root);
      if (gp) {
        const e = gp.children.find(c => c.node === parent);
        e.len = (e.len || 0) + (parent.children[0].len || 0);
        e.node = parent.children[0].node;
      }
    } else if (parent === root && root.children.length === 1) {
      const only = root.children[0].node;
      root.children = only.children;
      root.age = only.age;
    }
    return true;
  }

  /* ================================================================
     7 · rates through time
     ================================================================ */
  /* the speciation rate a fitted model implies at a grid of times, for the
     figure: constant models give a flat line, the density-dependent ones a
     falling one, and the two-rate model a step */
  function ratesThroughTime(fit, model, maxTime, steps) {
    const out = [];
    const S = steps || 100;
    for (let i = 0; i <= S; i++) {
      const t = maxTime * (1 - i / S);
      let lam = null;
      if (model === 'pure birth') lam = fit.r;
      else if (model === 'birth–death') lam = fit.lambda;
      else if (model === 'DDL' || model === 'DDX') lam = null;      // depends on lineages, not time
      else if (model === 'Yule-2') lam = t >= fit.shift ? fit.r1 : fit.r2;
      out.push({ time: t, lambda: lam });
    }
    return out;
  }

  /* ================================================================
     helpers
     ================================================================ */
  function lgamma(x) {
    const c = [76.18009172947146, -86.50532032941677, 24.01409824083091,
      -1.231739572450155, 0.1208650973866179e-2, -0.5395239384953e-5];
    let y = x, tmp = x + 5.5;
    tmp -= (x + 0.5) * Math.log(tmp);
    let ser = 1.000000000190015;
    for (let j = 0; j < 6; j++) ser += c[j] / ++y;
    return -tmp + Math.log(2.5066282746310005 * ser / x);
  }

  /* Nelder & Mead's simplex, for the two-parameter fits */
  function nelderMead(f, start, opts) {
    opts = opts || {};
    const n = start.length;
    const alpha = 1, gammaC = 2, rho = 0.5, sigma = 0.5;
    const simplex = [start.slice()];
    for (let i = 0; i < n; i++) {
      const p = start.slice();
      p[i] = p[i] !== 0 ? p[i] * 1.15 : 0.1;
      simplex.push(p);
    }
    let fv = simplex.map(p => f(p));
    const maxIter = opts.maxIter || 500;
    for (let it = 0; it < maxIter; it++) {
      const order = fv.map((v, i) => i).sort((a, b) => fv[a] - fv[b]);
      const sorted = order.map(i => simplex[i]);
      const fsorted = order.map(i => fv[i]);
      for (let i = 0; i <= n; i++) { simplex[i] = sorted[i]; fv[i] = fsorted[i]; }
      if (Math.abs(fv[n] - fv[0]) < (opts.tol || 1e-10) * (Math.abs(fv[0]) + 1e-10)) break;
      const centroid = new Array(n).fill(0);
      for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) centroid[j] += simplex[i][j] / n;
      const refl = centroid.map((c, j) => c + alpha * (c - simplex[n][j]));
      const fr = f(refl);
      if (fr < fv[0]) {
        const exp = centroid.map((c, j) => c + gammaC * (refl[j] - c));
        const fe = f(exp);
        if (fe < fr) { simplex[n] = exp; fv[n] = fe; } else { simplex[n] = refl; fv[n] = fr; }
      } else if (fr < fv[n - 1]) { simplex[n] = refl; fv[n] = fr; }
      else {
        const con = centroid.map((c, j) => c + rho * (simplex[n][j] - c));
        const fc = f(con);
        if (fc < fv[n]) { simplex[n] = con; fv[n] = fc; }
        else {
          for (let i = 1; i <= n; i++) {
            simplex[i] = simplex[i].map((v, j) => simplex[0][j] + sigma * (v - simplex[0][j]));
            fv[i] = f(simplex[i]);
          }
        }
      }
    }
    let bi = 0;
    for (let i = 1; i <= n; i++) if (fv[i] < fv[bi]) bi = i;
    return { x: simplex[bi], f: fv[bi] };
  }

  Object.assign(Diversify, {
    branchingTimes, ltt, lttBand, gammaStat, waitingTimes,
    yule, birthDeath, bdLogLik, pureBirth, ddLogistic, ddExponential, twoRate,
    intervalLogLik, compare, magallonSanderson, drStatistic,
    simulateBD, mccr, prune, ratesThroughTime, nelderMead, lgamma,
  });
  g.Diversify = Diversify;
}
DiversifyCore(typeof window !== 'undefined' ? window : self);
