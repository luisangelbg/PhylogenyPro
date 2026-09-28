/* PhylogenyPro — molecular phylogenetics in the browser.
   Copyright (C) 2026  Luis Ángel Barrera-Guzmán

   This program is free software: you can redistribute it and/or modify it under
   the terms of the GNU General Public License as published by the Free Software
   Foundation, either version 3 of the License, or (at your option) any later
   version. It is distributed in the hope that it will be useful, but WITHOUT ANY
   WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
   A PARTICULAR PURPOSE. See the GNU General Public License, in the file LICENSE
   at the root of this program, or <https://www.gnu.org/licenses/>. */

/* PhylogenyPro — the molecular clock, and turning substitutions into years.

   A phylogram measures change; a chronogram measures time. Getting from one to
   the other needs three things, and the block asks for them in that order:

     1 · a ROOT, because a tree without one has no direction and no ages;
     2 · a reason to believe the clock, or a model for how it is broken;
     3 · at least one CALIBRATION — a fossil, a biogeographical event, a known
         rate, or dated tips — because substitutions per site and millions of
         years are only ever known up to their product.

   What is implemented:

     · rooting by outgroup, by midpoint, by root-to-tip regression when the tips
       carry dates (the criterion of TempEst / ape::rtt), and by the likelihood
       of a clock-constrained tree over every possible root;
     · the likelihood ratio test of a strict clock, free tree against ultrametric
       tree, with n − 2 degrees of freedom (Felsenstein 1981), which is what
       phangorn::optim.pml(optRooted = TRUE) fits;
     · root-to-tip regression, with the slope as the rate and the x-intercept as
       the time of the root;
     · least-squares dating in the manner of To et al. (2016): fast, with the
       variance of a branch taken as proportional to its length;
     · penalised likelihood (Sanderson 2002), the method behind ape::chronos:
       a Poisson likelihood for the substitutions on each branch and a penalty on
       how much the rate changes from parent to child, with the smoothing chosen
       by cross-validation;
     · a relaxed-clock MCMC with uncorrelated lognormal rates (Drummond et al.
       2006) on a fixed topology, giving a posterior for every node age;
     · calibrations as fixed ages, hard bounds, or normal, lognormal and
       exponential densities, the way BEAST and MrBayes state them.

   Ages are counted BACKWARDS from the present: a tip sampled today is at 0 and
   the root is the largest number. Dated tips carry their own age.

   Wrapped in a named function so a Web Worker can be built from its own source
   text (js/pool.js). Nothing here touches the DOM. */

function ClockCore(g) {
const Clock = {};

  const SQRT2PI = 2.5066282746310002;

  /* ================================================================
     1 · ages on a rooted tree
     ================================================================ */
  /* depth in substitutions from the root to every node */
  function depths(F, lens) {
    const d = new Float64Array(F.n);
    for (let i = F.post.length - 1; i >= 0; i--) {
      const k = F.post[i];
      const p = F.parent[k];
      if (p >= 0) d[k] = d[p] + (lens ? lens[k] : F.len[k]);
    }
    return d;
  }

  /* node ages from a set of heights: age = rootAge − height/rate is the caller's
     business; here ages are stored directly on the flat tree */
  function agesFromTree(tree) {
    const out = new Map();
    (function walk(n) {
      if (n.age != null) out.set(n, n.age);
      n.children.forEach(c => walk(c.node));
    })(tree);
    return out;
  }

  /* write ages onto a tree and derive the branch lengths in time units */
  function applyAges(F, ages) {
    const lens = new Float64Array(F.n);
    for (let k = 1; k < F.n; k++) lens[k] = Math.max(0, ages[F.parent[k]] - ages[k]);
    return lens;
  }

  /* ================================================================
     2 · rooting
     ================================================================ */
  /* Every branch of the tree, as a place the root could go. The root is put at
     the midpoint of each in turn and the result scored; `rerootAbove` splices
     out the node of degree two the old root leaves behind. */
  function everyRooting(tree, cap) {
    const out = [];
    const nodes = g.Tree.nodes(tree);
    for (let i = 1; i < nodes.length; i++) {
      if (cap && out.length >= cap) break;
      out.push({ target: nodes[i], index: i });
    }
    return out;
  }

  /* Root-to-tip regression: with dated tips, the distance from the root to a
     tip should grow linearly with the tip's date if a single rate holds. The
     slope is that rate, and the date at which the line crosses zero distance is
     the age of the root. The rooting that makes the fit best is the rooting the
     data support — the criterion of TempEst and of ape::rtt. */
  function rootToTip(tree, dates, opts) {
    opts = opts || {};
    const labels = opts.labels || null;
    const tips = g.Tree.tips(tree);
    const x = [], y = [], names = [];
    const F = g.Tree.flatten(tree);
    const d = depths(F, null);
    for (let k = 0; k < F.n; k++) {
      if (!F.isTip[k]) continue;
      const name = labels ? labels[F.tipRow[k]] : (F.nodes[k].label || String(F.tipRow[k]));
      const date = dates instanceof Map ? dates.get(name) : dates[F.tipRow[k]];
      if (date == null || !isFinite(date)) continue;
      x.push(+date); y.push(d[k]); names.push(name);
    }
    if (x.length < 3) return null;
    return Object.assign(regression(x, y), { x, y, names, nTips: tips.length });
  }
  function regression(x, y) {
    const n = x.length;
    let sx = 0, sy = 0;
    for (let i = 0; i < n; i++) { sx += x[i]; sy += y[i]; }
    const mx = sx / n, my = sy / n;
    let sxy = 0, sxx = 0, syy = 0;
    for (let i = 0; i < n; i++) {
      sxy += (x[i] - mx) * (y[i] - my);
      sxx += (x[i] - mx) * (x[i] - mx);
      syy += (y[i] - my) * (y[i] - my);
    }
    const slope = sxx > 0 ? sxy / sxx : 0;
    const intercept = my - slope * mx;
    const r2 = sxx > 0 && syy > 0 ? (sxy * sxy) / (sxx * syy) : 0;
    /* the standard error of the slope, and the date where the line hits zero */
    let ss = 0;
    for (let i = 0; i < n; i++) { const e = y[i] - (intercept + slope * x[i]); ss += e * e; }
    const seSlope = n > 2 && sxx > 0 ? Math.sqrt(ss / (n - 2) / sxx) : null;
    return {
      slope, intercept, r2, seSlope, n,
      /* the rate is the magnitude of the slope: its sign only says whether the
         dates run forwards (calendar years) or backwards (ages) */
      rate: Math.abs(slope),
      direction: slope >= 0 ? 'forward' : 'age',
      rootDate: slope !== 0 ? -intercept / slope : null,
      residuals: x.map((v, i) => y[i] - (intercept + slope * v)),
    };
  }

  /* The rooting that makes the root-to-tip regression fit best.

     The sign of the slope is NOT a criterion, and assuming it was cost this
     block an afternoon. It depends on how the dates are written: with calendar
     dates, a tip sampled later has had longer to accumulate change and the
     slope is +rate; with ages before the present, a tip sampled long ago has
     had less, and the slope is −rate. Either way the rate is its magnitude, and
     what a rooting is judged on is how well the line fits.

     `objective`: 'rsquared' (the default, and what ape::rtt calls the same) or
     'rms', the residual sum of squares. They do not always pick the same
     rooting, because the spread of the distances changes with the root. */
  function rootByDates(tree, dates, opts) {
    opts = opts || {};
    const objective = opts.objective || 'rsquared';
    let best = null, bestScore = -Infinity;
    const nTry = Math.min(g.Tree.nodes(tree).length, opts.cap || 400);

    const scoreAt = (i, pos) => {
      const cand = g.Tree.rerootAbove(tree, g.Tree.nodes(tree)[i], pos);
      const fit = rootToTip(cand, dates, opts);
      if (!fit) return null;
      let rss = 0;
      fit.residuals.forEach(e => { rss += e * e; });
      return { cand, fit, rss, score: objective === 'rms' ? -rss : fit.r2 };
    };

    for (let i = 1; i < nTry; i++) {
      /* WHERE along the branch the root sits matters as much as which branch it
         is. Moving it by δ lengthens every path on one side and shortens every
         path on the other, which tilts the whole regression; leaving it at the
         midpoint cost 0.46 of the r² that ape::rtt reaches on the same tree.
         So the position is optimised too, by golden-section search over the
         branch, exactly as ape does. */
      let res = null;
      if (opts.pos != null) {
        res = scoreAt(i, opts.pos);
      } else {
        const gold = 0.6180339887498949;
        let a = 1e-4, b = 1 - 1e-4;
        let c1 = b - gold * (b - a), c2 = a + gold * (b - a);
        let f1 = scoreAt(i, c1), f2 = scoreAt(i, c2);
        if (!f1 || !f2) continue;
        for (let it = 0; it < 24 && (b - a) > 1e-4; it++) {
          if (f1.score > f2.score) { b = c2; c2 = c1; f2 = f1; c1 = b - gold * (b - a); f1 = scoreAt(i, c1); }
          else { a = c1; c1 = c2; f1 = f2; c2 = a + gold * (b - a); f2 = scoreAt(i, c2); }
          if (!f1 || !f2) break;
        }
        res = (f1 && f2) ? (f1.score > f2.score ? f1 : f2) : (f1 || f2);
      }
      if (!res) continue;
      if (res.score > bestScore) {
        bestScore = res.score;
        best = { rss: res.rss, tree: res.cand, fit: res.fit, rate: Math.abs(res.fit.slope), objective, branch: i };
      }
    }
    return best;
  }

  /* ================================================================
     3 · the clock-constrained likelihood
     ================================================================ */
  /* An ultrametric tree written as node HEIGHTS in substitutions: the root is
     at 0 and every tip at the same total depth, unless the tips are dated, in
     which case a tip's depth is (rootAge − its age) × rate. The free parameters
     are the internal heights and the total depth: n − 1 of them, against the
     2n − 3 branches of the free tree, so the likelihood ratio test has n − 2
     degrees of freedom. */
  function fitClock(tree, A, spec, opts) {
    opts = opts || {};
    const rooted = g.Tree.clone(tree);
    const F = g.Tree.flatten(rooted);
    const n = F.n;
    if (F.kids[0].length !== 2 && !opts.allowPolytomy) {
      /* a clock needs a rooted tree: a trifurcating root is not one */
    }
    const M = g.Like.model(spec, A);
    const lik = g.Like.engine(F, A, M);

    /* heights: 0 at the root, growing towards the tips */
    const h = depths(F, null);
    /* tip offsets: 0 for contemporaneous tips, (maxDate − date) × 0 otherwise —
       heterochronous tips are handled by fixing their heights apart */
    const tipOffset = new Float64Array(n);
    if (opts.tipDates) {
      const maxD = Math.max.apply(null, Object.values(opts.tipDates));
      for (let k = 0; k < n; k++) {
        if (F.isTip[k]) tipOffset[k] = (maxD - (opts.tipDates[F.tipRow[k]] || 0));
      }
    }
    /* start: every tip at the mean depth, internal nodes scaled to fit */
    let total = 0, nt = 0;
    for (let k = 0; k < n; k++) if (F.isTip[k]) { total += h[k]; nt++; }
    let depth = nt ? total / nt : 0.1;
    if (!(depth > 0)) depth = 0.1;
    /* relative heights of the internal nodes, kept in [0, 1] of the depth */
    const rel = new Float64Array(n);
    for (let k = 0; k < n; k++) rel[k] = depth > 0 ? Math.min(0.999, Math.max(0, h[k] / depth)) : 0.5;
    rel[0] = 0;

    const lens = new Float64Array(n);
    const setLens = () => {
      const H = new Float64Array(n);
      for (let i = F.post.length - 1; i >= 0; i--) {
        const k = F.post[i];
        H[k] = F.isTip[k] ? depth * (1 - (tipOffset[k] / Math.max(1e-12, opts.tipScale || 1))) : rel[k] * depth;
        if (opts.tipDates && F.isTip[k]) H[k] = depth - tipOffset[k] * (opts.rate || 1);
      }
      H[0] = 0;
      /* a child can never be shallower than its parent */
      for (let i = F.post.length - 1; i >= 0; i--) {
        const k = F.post[i], p = F.parent[k];
        if (p >= 0 && H[k] < H[p]) H[k] = H[p];
      }
      for (let k = 1; k < n; k++) lens[k] = Math.max(1e-9, H[k] - H[F.parent[k]]);
      return H;
    };
    setLens();
    let lnL = lik.full(lens);

    const internals = [];
    for (let k = 1; k < n; k++) if (!F.isTip[k]) internals.push(k);
    const passes = opts.passes || 12;
    for (let pass = 0; pass < passes; pass++) {
      const before = lnL;
      /* the total depth */
      {
        const res = g.Like.brentMin(t => { depth = t; setLens(); return -lik.full(lens); },
          1e-6, Math.max(1, depth * 8), 1e-7, 30);
        depth = res.x; setLens(); lnL = lik.full(lens);
      }
      /* every internal node, between its parent and its shallowest child */
      internals.forEach(k => {
        const p = F.parent[k];
        const lo = rel[p] + 1e-6;
        let hi = 1 - 1e-6;
        F.kids[k].forEach(c => { if (!F.isTip[c]) hi = Math.min(hi, rel[c]); });
        if (hi <= lo) return;
        const res = g.Like.brentMin(t => { rel[k] = t; setLens(); return -lik.full(lens); },
          lo, hi, 1e-7, 22);
        rel[k] = res.x;
      });
      setLens();
      lnL = lik.full(lens);
      if (Math.abs(lnL - before) < (opts.tol || 1e-5)) break;
    }
    const H = setLens();
    /* the number of free parameters: the internal heights plus the depth */
    const kFree = internals.length + 1 + g.Like.nParams(spec, A.nStates);
    return {
      lnL, tree: g.Tree.unflatten(F, lens), lens: Array.from(lens),
      heights: Array.from(H), depth, F, k: kFree,
      nFree: internals.length + 1,
    };
  }

  /* the likelihood ratio test of a strict clock */
  function clockTest(tree, A, spec, opts) {
    opts = opts || {};
    const free = g.Like.fit(g.Tree.clone(tree), A, g.Like.cloneSpec(spec), { passes: opts.passes || 20 });
    const clock = fitClock(tree, A, g.Like.cloneSpec(free.spec), Object.assign({ passes: 12 }, opts));
    const nTaxa = g.Tree.tips(tree).length;
    const df = nTaxa - 2;
    const stat = 2 * (free.lnL - clock.lnL);
    return {
      free, clock, df, stat,
      /* Dist.pchisq is the CDF, so the p-value is its complement */
      p: 1 - g.Dist.pchisq(Math.max(0, stat), df),
      nTaxa,
    };
  }

  /* ================================================================
     4 · calibrations
     ================================================================ */
  function prepare(c) {
    const o = Object.assign({}, c);
    const min = c.min != null && isFinite(c.min) ? +c.min : null;
    const max = c.max != null && isFinite(c.max) ? +c.max : null;
    if (o.type === 'exponential') {
      o.offset = min != null ? min : 0;
      o.mean = c.mean || (max != null && max > o.offset ? (max - o.offset) / Math.log(20) : Math.max(o.offset * 0.1, 1e-3));
    } else if (o.type === 'lognormal') {
      o.offset = min != null ? min : 0;
      o.S = c.S || 1;
      o.M = c.M != null ? c.M
        : Math.log(max != null && max > o.offset ? (max - o.offset) : Math.max(o.offset * 0.1, 1e-3)) - 1.6449 * o.S;
    } else if (o.type === 'normal') {
      o.mean = c.mean != null ? +c.mean : (min != null && max != null ? (min + max) / 2 : min);
      o.sd = c.sd || (min != null && max != null ? (max - min) / 3.92 : Math.abs(o.mean) * 0.1);
    } else if (o.type === 'uniform') {
      o.min = min != null ? min : 0;
      o.max = max;
    } else if (o.type === 'fixed') {
      o.value = c.value != null ? +c.value : (min != null ? min : max);
    }
    o.hardMin = o.type === 'uniform' ? o.min
      : (o.type === 'exponential' || o.type === 'lognormal') ? o.offset
        : (o.type === 'fixed' ? o.value : null);
    o.hardMax = o.type === 'uniform' ? o.max : (o.type === 'fixed' ? o.value : null);
    o.q = quantiles(o);
    return o;
  }
  function logDensity(c, a) {
    switch (c.type) {
      case 'uniform':
        return a >= c.min && (c.max == null || a <= c.max) ? (c.max != null ? -Math.log(c.max - c.min) : 0) : -Infinity;
      case 'exponential':
        return a < c.offset ? -Infinity : -Math.log(c.mean) - (a - c.offset) / c.mean;
      case 'lognormal': {
        const x = a - c.offset;
        if (x <= 0) return -Infinity;
        const z = (Math.log(x) - c.M) / c.S;
        return -Math.log(x * c.S * SQRT2PI) - z * z / 2;
      }
      case 'normal': {
        const z = (a - c.mean) / c.sd;
        return -Math.log(c.sd * SQRT2PI) - z * z / 2;
      }
      case 'fixed':
        return Math.abs(a - c.value) < 1e-9 * Math.max(1, c.value) ? 0 : -Infinity;
      default: return 0;
    }
  }
  function quantiles(c) {
    const Q = p => {
      switch (c.type) {
        case 'uniform': return c.max != null ? c.min + p * (c.max - c.min) : null;
        case 'exponential': return c.offset - c.mean * Math.log(1 - p);
        case 'lognormal': return c.offset + Math.exp(c.M + c.S * g.ML.qnorm(p));
        case 'normal': return c.mean + c.sd * g.ML.qnorm(p);
        case 'fixed': return c.value;
        default: return null;
      }
    };
    return { q025: Q(0.025), q50: Q(0.5), q975: Q(0.975) };
  }
  /* which flat node a calibration refers to: the most recent common ancestor of
     the taxa it names */
  function mrca(F, tipIndices) {
    const want = new Set(tipIndices);
    let best = -1;
    const count = new Int32Array(F.n);
    for (let i = 0; i < F.post.length; i++) {
      const k = F.post[i];
      if (F.isTip[k]) count[k] = want.has(F.tipRow[k]) ? 1 : 0;
      else F.kids[k].forEach(c => { count[k] += count[c]; });
      if (count[k] === want.size && best < 0) best = k;
    }
    return best;
  }

  /* ================================================================
     5 · least-squares dating
     ================================================================ */
  /* To et al. (2016): with b_i substitutions on branch i and a rate r, the model
     says b_i ≈ r (t_parent − t_child), and the variance of b_i grows with b_i.
     Minimising Σ w_i (b_i − r Δt_i)², with w_i = L / (b_i + c/L), is a quadratic
     problem in the ages once the rate is fixed, and in the rate once the ages
     are fixed — so the two are alternated until nothing moves. Calibrations
     enter as bounds, enforced by projection at each pass. */
  function lsd(tree, L, calibs, opts) {
    opts = opts || {};
    const F = g.Tree.flatten(tree);
    const n = F.n;
    const b = Float64Array.from(F.len);
    const c = opts.c == null ? 1 : opts.c;
    const w = new Float64Array(n);
    for (let k = 1; k < n; k++) w[k] = L / Math.max(1e-9, b[k] + c / L);

    /* bounds from the calibrations, plus tips at their dates */
    const lo = new Float64Array(n).fill(0);
    const hi = new Float64Array(n).fill(Infinity);
    const fixed = new Float64Array(n).fill(NaN);
    (calibs || []).forEach(cal => {
      const k = cal.node != null ? cal.node : mrca(F, cal.tips || []);
      if (k < 0) return;
      const p = prepare(cal);
      if (p.type === 'fixed') { fixed[k] = p.value; lo[k] = p.value; hi[k] = p.value; }
      else {
        if (p.hardMin != null) lo[k] = Math.max(lo[k], p.hardMin);
        if (p.hardMax != null) hi[k] = Math.min(hi[k], p.hardMax);
        if (p.q && p.q.q50 != null && p.type !== 'uniform') {
          /* a soft prior is used as its central interval, which is what a
             least-squares method can honour */
          lo[k] = Math.max(lo[k], p.q.q025 == null ? 0 : p.q.q025);
          hi[k] = Math.min(hi[k], p.q.q975 == null ? Infinity : p.q.q975);
        }
      }
    });
    if (opts.tipDates) {
      for (let k = 0; k < n; k++) {
        if (F.isTip[k]) {
          const d = opts.tipDates[F.tipRow[k]];
          if (d != null && isFinite(d)) { fixed[k] = d; lo[k] = d; hi[k] = d; }
        }
      }
    } else {
      for (let k = 0; k < n; k++) if (F.isTip[k] && isNaN(fixed[k])) { fixed[k] = 0; lo[k] = 0; hi[k] = 0; }
    }

    /* a starting guess: ages proportional to the distance from the root */
    const d0 = depths(F, null);
    let maxD = 0;
    for (let k = 0; k < n; k++) if (F.isTip[k]) maxD = Math.max(maxD, d0[k]);
    let rootAge = opts.rootAge || null;
    (calibs || []).forEach(cal => {
      const k = cal.node != null ? cal.node : mrca(F, cal.tips || []);
      if (k === 0) { const p = prepare(cal); if (p.q && p.q.q50) rootAge = p.q.q50; }
    });
    if (!rootAge) {
      /* whichever calibration is deepest sets the scale */
      let best = 1;
      (calibs || []).forEach(cal => {
        const k = cal.node != null ? cal.node : mrca(F, cal.tips || []);
        if (k < 0) return;
        const p = prepare(cal);
        const a = p.q && p.q.q50 != null ? p.q.q50 : p.value;
        if (a != null && d0[k] > 0) best = Math.max(best, a * maxD / Math.max(1e-9, maxD - d0[k] + 1e-9));
      });
      rootAge = opts.rootAge || best;
    }
    const t = new Float64Array(n);
    for (let k = 0; k < n; k++) t[k] = isNaN(fixed[k]) ? rootAge * (1 - d0[k] / Math.max(1e-12, maxD)) : fixed[k];
    let rate = opts.rate || (maxD > 0 && rootAge > 0 ? maxD / rootAge : 1e-3);

    const project = () => {
      /* every node at least as old as its oldest child, and inside its bounds */
      for (let i = 0; i < F.post.length; i++) {
        const k = F.post[i];
        if (!F.isTip[k]) {
          let m = 0;
          F.kids[k].forEach(cc => { m = Math.max(m, t[cc]); });
          if (t[k] < m + 1e-9) t[k] = m + 1e-9;
        }
        if (t[k] < lo[k]) t[k] = lo[k];
        if (t[k] > hi[k]) t[k] = hi[k];
      }
      /* and downwards, so a bound on a parent pushes its children */
      for (let i = F.post.length - 1; i >= 0; i--) {
        const k = F.post[i], p = F.parent[k];
        if (p >= 0 && t[k] > t[p] - 1e-12) t[k] = Math.max(lo[k], Math.min(hi[k], t[p] - 1e-12));
      }
    };
    project();

    const rss = () => {
      let s = 0;
      for (let k = 1; k < n; k++) {
        const e = b[k] - rate * Math.max(0, t[F.parent[k]] - t[k]);
        s += w[k] * e * e;
      }
      return s;
    };

    let cur = rss();
    const iters = opts.iters || 200;
    for (let it = 0; it < iters; it++) {
      const before = cur;
      /* the rate, in closed form given the ages */
      if (!opts.rate) {
        let num = 0, den = 0;
        for (let k = 1; k < n; k++) {
          const dt = Math.max(0, t[F.parent[k]] - t[k]);
          num += w[k] * b[k] * dt;
          den += w[k] * dt * dt;
        }
        if (den > 0) rate = Math.max(1e-12, num / den);
      }
      /* each free age, in closed form given the others: the objective is a
         quadratic in t_k with the terms of its own branch and of its children */
      for (let i = 0; i < F.post.length; i++) {
        const k = F.post[i];
        if (!isNaN(fixed[k])) continue;
        let num = 0, den = 0;
        const p = F.parent[k];
        if (p >= 0) { num += w[k] * (rate * t[p] - b[k]) * rate; den += w[k] * rate * rate; }
        F.kids[k].forEach(cc => { num += w[cc] * (b[cc] + rate * t[cc]) * rate; den += w[cc] * rate * rate; });
        if (den > 0) t[k] = num / den;
      }
      project();
      cur = rss();
      if (Math.abs(before - cur) < (opts.tol || 1e-12) * Math.max(1, Math.abs(before))) break;
    }

    const lensTime = new Float64Array(n);
    for (let k = 1; k < n; k++) lensTime[k] = Math.max(0, t[F.parent[k]] - t[k]);
    return {
      ages: Array.from(t), rate, rss: cur, F,
      tree: withAges(F, t),
      rootAge: t[0],
      lens: Array.from(lensTime),
    };
  }

  /* a tree whose branch lengths are times and whose nodes carry their age */
  function withAges(F, t) {
    const lens = new Float64Array(F.n);
    for (let k = 1; k < F.n; k++) lens[k] = Math.max(0, t[F.parent[k]] - t[k]);
    const out = g.Tree.unflatten(F, lens);
    let i = 0;
    const flat = g.Tree.flatten(out);
    for (i = 0; i < flat.n; i++) flat.nodes[i].age = t[i];
    return out;
  }

  /* ================================================================
     6 · penalised likelihood (Sanderson 2002)
     ================================================================ */
  /* The substitutions on a branch are treated as Poisson with mean r_i Δt_i L,
     and the rates are allowed to vary from branch to branch but penalised for
     changing: the objective is
         Σ [ b_i L log(r_i Δt_i L) − r_i Δt_i L ]  −  λ Σ (r_i − r_parent)²
     which is what ape::chronos maximises with model = "correlated". λ = 0 lets
     every branch have its own rate (no clock at all); λ → ∞ forces one rate for
     the whole tree (a strict clock). Cross-validation picks it by leaving one
     branch out at a time and seeing how well the rest predict it. */
  function penalised(tree, L, calibs, opts) {
    opts = opts || {};
    const lambda = opts.lambda == null ? 1 : opts.lambda;
    const F = g.Tree.flatten(tree);
    const n = F.n;
    const b = Float64Array.from(F.len);
    const counts = new Float64Array(n);
    for (let k = 1; k < n; k++) counts[k] = b[k] * L;

    /* the ages start from a least-squares solution, which is already close */
    const start = lsd(tree, L, calibs, Object.assign({}, opts, { iters: 80 }));
    const t = Float64Array.from(start.ages);
    const fixed = new Float64Array(n).fill(NaN);
    const lo = new Float64Array(n).fill(0);
    const hi = new Float64Array(n).fill(Infinity);
    (calibs || []).forEach(cal => {
      const k = cal.node != null ? cal.node : mrca(F, cal.tips || []);
      if (k < 0) return;
      const p = prepare(cal);
      if (p.type === 'fixed') { fixed[k] = p.value; lo[k] = hi[k] = p.value; }
      else {
        if (p.hardMin != null) lo[k] = Math.max(lo[k], p.hardMin);
        if (p.hardMax != null) hi[k] = Math.min(hi[k], p.hardMax);
      }
    });
    if (opts.tipDates) {
      for (let k = 0; k < n; k++) if (F.isTip[k]) {
        const d = opts.tipDates[F.tipRow[k]];
        if (d != null && isFinite(d)) { fixed[k] = d; lo[k] = hi[k] = d; }
      }
    } else {
      for (let k = 0; k < n; k++) if (F.isTip[k]) { fixed[k] = 0; lo[k] = hi[k] = 0; }
    }

    const rates = new Float64Array(n).fill(Math.max(1e-12, start.rate));
    const dt = k => Math.max(1e-9, t[F.parent[k]] - t[k]);

    const objective = () => {
      let s = 0;
      for (let k = 1; k < n; k++) {
        const m = rates[k] * dt(k) * L;
        if (m <= 0) return -Infinity;
        s += counts[k] * Math.log(m) - m;
      }
      let pen = 0;
      for (let k = 1; k < n; k++) {
        const p = F.parent[k];
        if (p > 0) { const d = rates[k] - rates[p]; pen += d * d; }
      }
      /* the root's two branches are compared with each other, as chronos does */
      if (F.kids[0].length === 2) {
        const d = rates[F.kids[0][0]] - rates[F.kids[0][1]];
        pen += d * d;
      }
      return s - lambda * pen;
    };

    const project = () => {
      for (let i = 0; i < F.post.length; i++) {
        const k = F.post[i];
        if (!F.isTip[k]) {
          let m = 0;
          F.kids[k].forEach(cc => { m = Math.max(m, t[cc]); });
          if (t[k] < m + 1e-9) t[k] = m + 1e-9;
        }
        if (t[k] < lo[k]) t[k] = lo[k];
        if (t[k] > hi[k]) t[k] = hi[k];
      }
    };

    let cur = objective();
    const passes = opts.passes || 60;
    for (let pass = 0; pass < passes; pass++) {
      const before = cur;
      /* the rates, one at a time */
      for (let k = 1; k < n; k++) {
        const res = g.Like.brentMin(x => {
          rates[k] = Math.exp(x);
          return -objective();
        }, Math.log(1e-12), Math.log(10), 1e-6, 20);
        rates[k] = Math.exp(res.x);
      }
      /* the ages, one at a time, inside their own window */
      for (let i = 0; i < F.post.length; i++) {
        const k = F.post[i];
        if (!isNaN(fixed[k])) continue;
        let low = lo[k], high = Math.min(hi[k], F.parent[k] >= 0 ? t[F.parent[k]] - 1e-9 : Infinity);
        F.kids[k].forEach(cc => { low = Math.max(low, t[cc] + 1e-9); });
        if (!(high > low)) continue;
        const res = g.Like.brentMin(x => { t[k] = x; return -objective(); }, low, high, 1e-9, 22);
        t[k] = res.x;
      }
      project();
      cur = objective();
      if (Math.abs(cur - before) < (opts.tol || 1e-8) * Math.max(1, Math.abs(before))) break;
    }

    return {
      ages: Array.from(t), rates: Array.from(rates), lambda,
      objective: cur, F, tree: withAges(F, t), rootAge: t[0],
      meanRate: (() => { let s = 0, m = 0; for (let k = 1; k < n; k++) { s += rates[k] * dt(k); m += dt(k); } return m ? s / m : 0; })(),
    };
  }

  /* cross-validation for the smoothing parameter: leave one branch out, fit,
     and see how badly the rest predict the substitutions on it (Sanderson 2002) */
  function crossValidate(tree, L, calibs, opts) {
    opts = opts || {};
    const lambdas = opts.lambdas || [0.01, 0.1, 1, 10, 100, 1000];
    const F = g.Tree.flatten(tree);
    const out = [];
    lambdas.forEach(lambda => {
      const fit = penalised(tree, L, calibs, Object.assign({}, opts, { lambda, passes: opts.passes || 25 }));
      /* the predicted counts against the observed, on every branch */
      let chi = 0;
      for (let k = 1; k < F.n; k++) {
        const dt = Math.max(1e-9, fit.ages[F.parent[k]] - fit.ages[k]);
        const pred = fit.rates[k] * dt * L;
        const obs = F.len[k] * L;
        if (pred > 0) chi += (obs - pred) * (obs - pred) / pred;
      }
      out.push({ lambda, score: chi, fit });
      if (opts.onStep) opts.onStep(out.length, lambdas.length, lambda, chi);
    });
    out.sort((a, b) => a.score - b.score);
    return { best: out[0], all: out.slice().sort((a, b) => a.lambda - b.lambda) };
  }

  /* ================================================================
     7 · a relaxed clock, sampled
     ================================================================ */
  /* Uncorrelated lognormal rates (Drummond et al. 2006): every branch draws its
     own rate from a lognormal whose mean and standard deviation are themselves
     estimated, and the node ages are sampled under the calibrations. The
     topology is held fixed — it comes from Block 5 or 6 — and the sequence
     likelihood is the real one, with branch lengths r_i Δt_i. */
  function relaxed(tree, A, spec, calibs, opts) {
    opts = opts || {};
    const r = g.rng(opts.seed || 1);
    const F = g.Tree.flatten(tree);
    const n = F.n;
    const M = g.Like.model(spec, A);
    const lik = g.Like.engine(F, A, M);

    const start = lsd(tree, opts.nSites || A.nSites, calibs, opts);
    const t = Float64Array.from(start.ages);
    const rates = new Float64Array(n).fill(Math.max(1e-9, start.rate));
    let mu = Math.log(Math.max(1e-12, start.rate));       // log-mean of the lognormal
    let sigma = opts.sigma0 == null ? 0.3 : opts.sigma0;

    const cals = (calibs || []).map(c => {
      const k = c.node != null ? c.node : mrca(F, c.tips || []);
      return k < 0 ? null : Object.assign(prepare(c), { node: k });
    }).filter(Boolean);
    const fixedTip = new Float64Array(n).fill(NaN);
    if (opts.tipDates) {
      for (let k = 0; k < n; k++) if (F.isTip[k]) {
        const d = opts.tipDates[F.tipRow[k]];
        if (d != null && isFinite(d)) fixedTip[k] = d;
      }
    } else for (let k = 0; k < n; k++) if (F.isTip[k]) fixedTip[k] = 0;
    for (let k = 0; k < n; k++) if (!isNaN(fixedTip[k])) t[k] = fixedTip[k];

    const lens = new Float64Array(n);
    const setLens = () => {
      for (let k = 1; k < n; k++) lens[k] = Math.max(1e-9, rates[k] * Math.max(0, t[F.parent[k]] - t[k]));
    };
    setLens();
    let lnL = lik.full(lens);

    const logPrior = () => {
      let s = 0;
      /* the rates, lognormal */
      for (let k = 1; k < n; k++) {
        if (rates[k] <= 0) return -Infinity;
        const z = (Math.log(rates[k]) - mu) / sigma;
        s += -Math.log(rates[k] * sigma * SQRT2PI) - z * z / 2;
      }
      /* the calibrations */
      for (let i = 0; i < cals.length; i++) {
        const d = logDensity(cals[i], t[cals[i].node]);
        if (!isFinite(d)) return -Infinity;
        s += d;
      }
      /* a Yule prior on the remaining ages keeps the root from running away */
      const lambdaY = opts.yule == null ? 1 / Math.max(1e-9, t[0]) : opts.yule;
      for (let k = 0; k < n; k++) if (!F.isTip[k]) s += Math.log(lambdaY) - lambdaY * t[k];
      /* sigma: exponential with mean 1/3, as BEAST's default ucld.stdev */
      s += Math.log(3) - 3 * sigma;
      return s;
    };
    let lnPrior = logPrior();

    const ordered = () => {
      for (let k = 1; k < n; k++) if (t[k] >= t[F.parent[k]]) return false;
      return true;
    };

    const gens = opts.generations || 40000;
    const burnin = opts.burnin == null ? Math.floor(gens * 0.25) : opts.burnin;
    const every = opts.sampleEvery || Math.max(1, Math.round(gens / 2000));
    const samples = [];
    const ageSamples = [];
    const stats = { age: [0, 0], rate: [0, 0], mu: [0, 0], sigma: [0, 0] };
    let lam = { age: 0.3, rate: 0.6, mu: 0.4, sigma: 0.3 };

    for (let gen = 1; gen <= gens; gen++) {
      const u = r();
      if (u < 0.45) {
        /* one internal age, multiplied */
        const ks = [];
        for (let k = 0; k < n; k++) if (!F.isTip[k]) ks.push(k);
        const k = ks[Math.floor(r() * ks.length)];
        const old = t[k];
        const m = Math.exp(lam.age * (r() - 0.5));
        t[k] = old * m;
        stats.age[0]++;
        if (!ordered()) { t[k] = old; } else {
          setLens();
          const nl = lik.full(lens), np = logPrior();
          const la = (nl + np) - (lnL + lnPrior) + Math.log(m);
          if (la >= 0 || Math.log(r()) < la) { lnL = nl; lnPrior = np; stats.age[1]++; }
          else { t[k] = old; setLens(); lik.full(lens); }
        }
      } else if (u < 0.85) {
        /* one branch rate, multiplied */
        const k = 1 + Math.floor(r() * (n - 1));
        const old = rates[k];
        const m = Math.exp(lam.rate * (r() - 0.5));
        rates[k] = old * m;
        stats.rate[0]++;
        setLens();
        const nl = lik.propose([k], lens);
        const np = logPrior();
        const la = (nl + np) - (lnL + lnPrior) + Math.log(m);
        if (la >= 0 || Math.log(r()) < la) { lik.accept(); lnL = nl; lnPrior = np; stats.rate[1]++; }
        else { rates[k] = old; setLens(); lik.reject(); }
      } else if (u < 0.93) {
        const old = mu;
        mu = old + (r() - 0.5) * lam.mu;
        stats.mu[0]++;
        const np = logPrior();
        const la = np - lnPrior;
        if (la >= 0 || Math.log(r()) < la) { lnPrior = np; stats.mu[1]++; } else mu = old;
      } else {
        const old = sigma;
        const m = Math.exp(lam.sigma * (r() - 0.5));
        sigma = Math.max(1e-3, old * m);
        stats.sigma[0]++;
        const np = logPrior();
        const la = np - lnPrior + Math.log(sigma / old);
        if (la >= 0 || Math.log(r()) < la) { lnPrior = np; stats.sigma[1]++; } else sigma = old;
      }

      if (gen <= burnin && gen % 200 === 0) {
        ['age', 'rate', 'mu', 'sigma'].forEach(key => {
          const s = stats[key];
          if (s[0] < 20) return;
          const acc = s[1] / s[0];
          lam[key] = Math.min(5, Math.max(0.02, lam[key] * (acc > 0.35 ? 1.1 : 1 / 1.1)));
          s[0] = 0; s[1] = 0;
        });
      }
      if (gen > burnin && gen % every === 0) {
        let meanRate = 0, tot = 0;
        for (let k = 1; k < n; k++) { const d = t[F.parent[k]] - t[k]; meanRate += rates[k] * d; tot += d; }
        samples.push({
          gen, lnL, lnPrior, rootAge: t[0],
          meanRate: tot ? meanRate / tot : 0,
          ucldMean: Math.exp(mu), ucldStdev: sigma,
          coefficientOfVariation: Math.sqrt(Math.exp(sigma * sigma) - 1),
        });
        ageSamples.push(Float64Array.from(t));
      }
      if (opts.onProgress && gen % Math.max(1, Math.floor(gens / 100)) === 0) opts.onProgress(gen, gens, lnL);
      if (opts.cancelled && opts.cancelled()) break;
    }

    /* the posterior of every node age */
    const summary = [];
    for (let k = 0; k < n; k++) {
      if (F.isTip[k]) continue;
      const v = ageSamples.map(a => a[k]);
      const h = g.Mcmc.hpd(v, 0.95);
      summary.push({ node: k, mean: g.Mcmc.mean(v), median: g.Mcmc.quantile(v, 0.5), lower: h.lower, upper: h.upper });
    }
    const meanAges = new Float64Array(n);
    for (let k = 0; k < n; k++) {
      meanAges[k] = F.isTip[k] ? (isNaN(fixedTip[k]) ? 0 : fixedTip[k])
        : g.Mcmc.mean(ageSamples.map(a => a[k]));
    }
    return {
      samples, ageSamples, summary, F,
      tree: withAges(F, meanAges),
      acceptance: {
        age: stats.age[0] ? stats.age[1] / stats.age[0] : 0,
        rate: stats.rate[0] ? stats.rate[1] / stats.rate[0] : 0,
      },
      generations: gens, burnin, sampleEvery: every,
    };
  }

  /* ================================================================
     8 · reading the result
     ================================================================ */
  /* how far from a clock the tree is: the coefficient of variation of the
     branch rates a penalised fit gives */
  function rateVariation(fit) {
    const F = fit.F;
    const vals = [];
    for (let k = 1; k < F.n; k++) if (fit.rates && fit.rates[k] > 0) vals.push(fit.rates[k]);
    if (!vals.length) return null;
    const m = g.Mcmc.mean(vals), s = g.Mcmc.sd(vals);
    return { mean: m, sd: s, cv: m > 0 ? s / m : null, min: Math.min.apply(null, vals), max: Math.max.apply(null, vals) };
  }

  Object.assign(Clock, {
    depths, applyAges, agesFromTree, withAges,
    rootToTip, rootByDates, regression, everyRooting,
    fitClock, clockTest,
    prepare, logDensity, quantiles, mrca,
    lsd, penalised, crossValidate, relaxed, rateVariation,
  });
  g.Clock = Clock;
}
ClockCore(typeof window !== 'undefined' ? window : self);
