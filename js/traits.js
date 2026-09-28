/* PhylogenyPro — molecular phylogenetics in the browser.
   Copyright (C) 2026  Luis Ángel Barrera-Guzmán

   This program is free software: you can redistribute it and/or modify it under
   the terms of the GNU General Public License as published by the Free Software
   Foundation, either version 3 of the License, or (at your option) any later
   version. It is distributed in the hope that it will be useful, but WITHOUT ANY
   WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
   A PARTICULAR PURPOSE. See the GNU General Public License, in the file LICENSE
   at the root of this program, or <https://www.gnu.org/licenses/>. */

/* PhylogenyPro — how characters change along a tree.

   Two very different kinds of question share this block, and they share a
   reason: species are not independent observations. A trait measured on two
   sister species carries almost the same history twice, and every method here
   exists to stop that shared history being mistaken for evidence.

   Discrete characters
     · the Mk model of Lewis (2001) with equal rates, symmetric rates or all
       rates different, fitted by maximum likelihood over a general rate matrix
       (so ARD, which is not reversible, works too);
     · ancestral states, both the way ape::ace reports them — the conditional
       likelihoods of the downward pass — and the true marginal, which needs a
       second pass from the root and is what phytools and corHMM report;
     · stochastic character mapping (Huelsenbeck et al. 2003; Nielsen 2002):
       histories drawn along the branches conditioned on the tips, which turn
       "probably state 0 here" into "how long the lineage spent in each state".

   Continuous characters
     · Brownian motion, Ornstein–Uhlenbeck and the early burst, all as a
       multivariate normal whose covariance is a function of the shared path
       from the root, in the parameterisation geiger::fitContinuous uses;
     · Felsenstein's (1985) independent contrasts, which make the observations
       independent again;
     · phylogenetic generalised least squares, whose slope through the origin is
       the same number as the regression on contrasts — an identity worth
       checking, because when it fails one of the two is wrong;
     · Pagel's λ and Blomberg's K, the two measures of phylogenetic signal that
       answer different questions and are often quoted as if they answered the
       same one.

   Wrapped in a named function so a Web Worker can be built from its own source
   text (js/pool.js). Nothing here touches the DOM. */

function TraitsCore(g) {
const Traits = {};

  /* ================================================================
     1 · linear algebra, only as much as is needed
     ================================================================ */
  /* Cholesky of a symmetric positive-definite matrix: returns the lower
     triangle L with A = L Lᵀ, or null when A is not positive definite, which
     is how a badly conditioned covariance announces itself. */
  function chol(A) {
    const n = A.length;
    const L = Array.from({ length: n }, () => new Float64Array(n));
    for (let i = 0; i < n; i++) {
      for (let j = 0; j <= i; j++) {
        let s = A[i][j];
        for (let k = 0; k < j; k++) s -= L[i][k] * L[j][k];
        if (i === j) {
          if (!(s > 0)) return null;
          L[i][i] = Math.sqrt(s);
        } else L[i][j] = s / L[j][j];
      }
    }
    return L;
  }
  /* solve A z = b given the Cholesky factor */
  function cholSolve(L, b) {
    const n = L.length;
    const y = new Float64Array(n), z = new Float64Array(n);
    for (let i = 0; i < n; i++) {
      let s = b[i];
      for (let k = 0; k < i; k++) s -= L[i][k] * y[k];
      y[i] = s / L[i][i];
    }
    for (let i = n - 1; i >= 0; i--) {
      let s = y[i];
      for (let k = i + 1; k < n; k++) s -= L[k][i] * z[k];
      z[i] = s / L[i][i];
    }
    return z;
  }
  const cholLogDet = L => { let s = 0; for (let i = 0; i < L.length; i++) s += 2 * Math.log(L[i][i]); return s; };

  /* the exponential of a general matrix, by scaling and squaring with a Taylor
     series: the rate matrices here are small and this is accurate to the last
     bit for them, and unlike an eigen-decomposition it does not care whether
     the matrix is reversible */
  function expm(Q, t) {
    const n = Q.length;
    let norm = 0;
    for (let i = 0; i < n; i++) { let s = 0; for (let j = 0; j < n; j++) s += Math.abs(Q[i][j]); norm = Math.max(norm, s); }
    norm *= Math.abs(t);
    let sq = 0;
    while (norm > 0.5) { norm /= 2; sq++; }
    const h = t / Math.pow(2, sq);
    let term = eye(n), out = eye(n);
    for (let k = 1; k <= 24; k++) {
      term = scale(matmul(term, Q), h / k);
      out = add(out, term);
      let m = 0;
      for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) m = Math.max(m, Math.abs(term[i][j]));
      if (m < 1e-18) break;
    }
    for (let s = 0; s < sq; s++) out = matmul(out, out);
    return out;
  }
  const eye = n => Array.from({ length: n }, (_, i) => Float64Array.from({ length: n }, (_, j) => (i === j ? 1 : 0)));
  function matmul(A, B) {
    const n = A.length, m = B[0].length, k2 = B.length;
    const C = Array.from({ length: n }, () => new Float64Array(m));
    for (let i = 0; i < n; i++) for (let k = 0; k < k2; k++) {
      const a = A[i][k];
      if (a === 0) continue;
      for (let j = 0; j < m; j++) C[i][j] += a * B[k][j];
    }
    return C;
  }
  const add = (A, B) => A.map((r, i) => Float64Array.from(r, (v, j) => v + B[i][j]));
  const scale = (A, s) => A.map(r => Float64Array.from(r, v => v * s));

  /* ================================================================
     2 · the phylogenetic covariance
     ================================================================ */
  /* C[i][j] is the length of the path the two tips share with the root, which
     is what ape::vcv.phylo returns. Everything continuous in this block is a
     multivariate normal built on it. */
  function vcv(tree, nTaxa) {
    const F = g.Tree.flatten(tree);
    const n = nTaxa || g.Tree.tips(tree).length;
    const depth = new Float64Array(F.n);
    for (let i = F.post.length - 1; i >= 0; i--) {
      const k = F.post[i], p = F.parent[k];
      if (p >= 0) depth[k] = depth[p] + (F.len[k] || 0);
    }
    const C = Array.from({ length: n }, () => new Float64Array(n));
    /* the tips below every node, so that a node contributes its depth to every
       pair that meets there */
    const below = new Array(F.n);
    for (let i = 0; i < F.post.length; i++) {
      const k = F.post[i];
      below[k] = F.isTip[k] ? [F.tipRow[k]] : [].concat.apply([], F.kids[k].map(c => below[c]));
    }
    for (let i = 0; i < F.post.length; i++) {
      const k = F.post[i];
      if (F.isTip[k]) { C[F.tipRow[k]][F.tipRow[k]] = depth[k]; continue; }
      /* every pair whose most recent common ancestor is k shares exactly depth[k] */
      const kids = F.kids[k];
      for (let a = 0; a < kids.length; a++) {
        for (let b = a + 1; b < kids.length; b++) {
          below[kids[a]].forEach(u => below[kids[b]].forEach(v => {
            C[u][v] = depth[k]; C[v][u] = depth[k];
          }));
        }
      }
    }
    return C;
  }

  /* Brownian motion does not care whether the tips are contemporaries, but OU
     and the early burst do: both are written in terms of the time from the root,
     and on a tree whose tips sit at different depths — a phylogram in
     substitutions, or a tree with fossils — they answer a question nobody asked.
     The block warns instead of refusing, because the fit is still computable and
     a reader who knows this may want it. */
  function isUltrametric(C, tol) {
    tol = tol == null ? 1e-4 : tol;
    let lo = Infinity, hi = -Infinity;
    for (let i = 0; i < C.length; i++) { lo = Math.min(lo, C[i][i]); hi = Math.max(hi, C[i][i]); }
    return { ultrametric: hi - lo <= tol * Math.max(1e-12, hi), spread: hi - lo, depth: hi };
  }

  /* the three covariance shapes, from the shared-path matrix */
  function covBM(C, sigma2) { return C.map(r => Float64Array.from(r, v => sigma2 * v)); }
  function covOU(C, sigma2, alpha, T) {
    if (!(alpha > 1e-12)) return covBM(C, sigma2);
    return C.map(r => Float64Array.from(r, s =>
      (sigma2 / (2 * alpha)) * (Math.exp(-2 * alpha * (T - s)) - Math.exp(-2 * alpha * T))));
  }
  /* (exp(a s) − 1) / a tends to s as a tends to zero, and computing it that way
     loses every significant digit long before a reaches zero: exp(x) − 1 for
     tiny x is the textbook cancellation. expm1 is exactly the function that
     does not lose them. */
  function covEB(C, sigma2, a) {
    if (a === 0) return covBM(C, sigma2);
    return C.map(r => Float64Array.from(r, s => sigma2 * Math.expm1(a * s) / a));
  }
  /* Pagel's transformation: the off-diagonals shrink towards a star */
  function covLambda(C, lambda) {
    return C.map((r, i) => Float64Array.from(r, (v, j) => (i === j ? v : lambda * v)));
  }

  /* ================================================================
     3 · the Gaussian likelihood and the three models
     ================================================================ */
  /* With V the covariance and the root state a free parameter, the maximum over
     the root and over the overall scale is available in closed form, so only
     the shape parameters ever need searching. */
  function gaussianFit(V, x) {
    const n = x.length;
    const L = chol(V);
    if (!L) return null;
    const one = Float64Array.from({ length: n }, () => 1);
    const Vi1 = cholSolve(L, one);
    const Vix = cholSolve(L, x);
    let a = 0, b = 0;
    for (let i = 0; i < n; i++) { a += one[i] * Vix[i]; b += one[i] * Vi1[i]; }
    const z0 = a / b;
    const d = Float64Array.from(x, (v, i) => v - z0);
    const Vid = cholSolve(L, d);
    let q = 0;
    for (let i = 0; i < n; i++) q += d[i] * Vid[i];
    /* V was built with the scale already in it, so this is the likelihood at
       that scale, not at its maximum */
    const lnL = -0.5 * (n * Math.log(2 * Math.PI) + cholLogDet(L) + q);
    return { z0, lnL, quad: q, logDet: cholLogDet(L) };
  }
  /* the same, with the scale profiled out: V must be the covariance for σ² = 1 */
  function gaussianProfile(C1, x) {
    const n = x.length;
    const L = chol(C1);
    if (!L) return null;
    const one = Float64Array.from({ length: n }, () => 1);
    const Vi1 = cholSolve(L, one), Vix = cholSolve(L, x);
    let a = 0, b = 0;
    for (let i = 0; i < n; i++) { a += one[i] * Vix[i]; b += one[i] * Vi1[i]; }
    const z0 = a / b;
    const d = Float64Array.from(x, (v, i) => v - z0);
    const Vid = cholSolve(L, d);
    let q = 0;
    for (let i = 0; i < n; i++) q += d[i] * Vid[i];
    const sigma2 = q / n;
    const lnL = -0.5 * (n * Math.log(2 * Math.PI) + n * Math.log(sigma2) + cholLogDet(L) + n);
    return { z0, sigma2, lnL };
  }

  function fitBM(tree, x, opts) {
    opts = opts || {};
    const C = opts.C || vcv(tree, x.length);
    const f = gaussianProfile(C, x);
    if (!f) return null;
    return { model: 'BM', sigma2: f.sigma2, z0: f.z0, lnL: f.lnL, k: 2, C,
      AIC: -2 * f.lnL + 2 * 2 };
  }
  function fitOU(tree, x, opts) {
    opts = opts || {};
    const C = opts.C || vcv(tree, x.length);
    let T = 0;
    for (let i = 0; i < C.length; i++) T = Math.max(T, C[i][i]);
    /* α is a rate, so its useful range is set by the depth of the tree and not
       by any fixed pair of numbers: on a tree 0.08 substitutions deep, an α of
       20 is already a pull the data cannot see, and on one 300 Ma deep an α of
       20 is instantaneous. The bounds are αT ∈ [1e-6, 200], as geiger does. */
    const obj = la => {
      const alpha = Math.exp(la);
      const f = gaussianProfile(covOU(C, 1, alpha, T), x);
      return f ? -f.lnL : 1e100;
    };
    const lo = Math.log(1e-6 / T), hi = Math.log(200 / T);
    const best = g.Like.brentMin(obj, lo, hi, 1e-9, 80);
    const alpha = Math.exp(best.x);
    const f = gaussianProfile(covOU(C, 1, alpha, T), x);
    return { model: 'OU', sigma2: f.sigma2, alpha, z0: f.z0, lnL: f.lnL, k: 3, C, T,
      halfLife: Math.log(2) / alpha,
      /* at the lower bound OU has collapsed into Brownian motion, and the
         half-life it reports is an artefact of where the search stopped. The
         test is relative, because Brent stops near the bound and not on it. */
      atBound: alpha <= 1.01 * Math.exp(lo) || alpha >= 0.99 * Math.exp(hi),
      AIC: -2 * f.lnL + 2 * 3 };
  }
  function fitEB(tree, x, opts) {
    opts = opts || {};
    const C = opts.C || vcv(tree, x.length);
    let T = 0;
    for (let i = 0; i < C.length; i++) T = Math.max(T, C[i][i]);
    /* geiger searches a in [log(1e-5)/T, -1e-6]: the rate may only fall */
    const lo = Math.log(1e-5) / T, hi = -1e-6;
    const obj = a => {
      const f = gaussianProfile(covEB(C, 1, a), x);
      return f ? -f.lnL : 1e100;
    };
    const best = g.Like.brentMin(obj, lo, hi, 1e-10, 60);
    const a = best.x;
    const f = gaussianProfile(covEB(C, 1, a), x);
    return { model: 'EB', sigma2: f.sigma2, a, z0: f.z0, lnL: f.lnL, k: 3, C,
      atBound: Math.abs(a - hi) <= 0.01 * Math.abs(hi) || Math.abs(a - lo) <= 0.01 * Math.abs(lo),
      AIC: -2 * f.lnL + 2 * 3 };
  }
  function compareContinuous(tree, x) {
    const C = vcv(tree, x.length);
    const rows = [fitBM(tree, x, { C }), fitOU(tree, x, { C }), fitEB(tree, x, { C })].filter(Boolean);
    const best = rows.reduce((a, b) => (b.AIC < a.AIC ? b : a), rows[0]);
    let sw = 0;
    rows.forEach(r => { r.dAIC = r.AIC - best.AIC; r.w = Math.exp(-0.5 * r.dAIC); sw += r.w; });
    rows.forEach(r => { r.w /= sw; });
    rows.sort((a, b) => a.AIC - b.AIC);
    return { rows, best: rows[0] };
  }

  /* ================================================================
     4 · independent contrasts
     ================================================================ */
  /* Felsenstein (1985). At every internal node the two daughter values are
     subtracted and the difference divided by the square root of the variance it
     should have; the node then takes a weighted average of its daughters and
     its own branch is lengthened to account for the estimation. The contrasts
     are independent and identically distributed under Brownian motion, which is
     what makes an ordinary regression on them legitimate. */
  function pic(tree, x, opts) {
    opts = opts || {};
    const F = g.Tree.flatten(tree);
    const val = new Float64Array(F.n);
    const bl = Float64Array.from(F.len);
    const contrasts = [], variances = [], nodes = [];
    for (let k = 0; k < F.n; k++) if (F.isTip[k]) val[k] = x[F.tipRow[k]];
    for (let i = 0; i < F.post.length; i++) {
      const k = F.post[i];
      if (F.isTip[k]) continue;
      const kids = F.kids[k];
      if (!kids.length) continue;
      /* a polytomy is walked as a ladder, which is what ape does when it meets
         one: the answer then depends on the order of the daughters, and that is
         a property of the data, not of this code */
      let curV = val[kids[0]], curB = bl[kids[0]];
      for (let j = 1; j < kids.length; j++) {
        const b = kids[j], vb = bl[b];
        contrasts.push((curV - val[b]) / Math.sqrt(curB + vb));
        variances.push(curB + vb);
        nodes.push(k);
        curV = (vb * curV + curB * val[b]) / (curB + vb);
        curB = (curB * vb) / (curB + vb);
      }
      val[k] = curV;
      bl[k] = (F.len[k] || 0) + curB;
    }
    return { contrasts, variances, nodes, rootValue: val[0], rootVariance: bl[0] };
  }

  /* ================================================================
     5 · phylogenetic generalised least squares
     ================================================================ */
  /* The regression species deserve: the residuals are correlated exactly as the
     tree says, so the estimate weights the observations by how much independent
     information each carries. Its slope through the origin is the same number
     as an ordinary regression on the contrasts — the identity of Garland &
     Ives (2000), which is the cheapest check there is that both are right. */
  function pgls(tree, y, X, opts) {
    opts = opts || {};
    const n = y.length;
    const C = opts.C || vcv(tree, n);
    let V = C;
    if (opts.lambda != null) V = covLambda(C, opts.lambda);
    if (opts.alpha != null) {
      let T = 0; for (let i = 0; i < C.length; i++) T = Math.max(T, C[i][i]);
      V = covOU(C, 1, opts.alpha, T);
    }
    /* the design matrix: an intercept plus every predictor */
    const preds = Array.isArray(X[0]) ? X : [X];
    const p = preds.length + 1;
    const D = Array.from({ length: n }, (_, i) => {
      const row = new Float64Array(p);
      row[0] = 1;
      for (let j = 0; j < preds.length; j++) row[j + 1] = preds[j][i];
      return row;
    });
    const L = chol(V);
    if (!L) return null;
    /* (Dᵀ V⁻¹ D) β = Dᵀ V⁻¹ y */
    const ViD = [];
    for (let j = 0; j < p; j++) ViD.push(cholSolve(L, Float64Array.from(D.map(r => r[j]))));
    const Viy = cholSolve(L, Float64Array.from(y));
    const A = Array.from({ length: p }, () => new Float64Array(p));
    const b = new Float64Array(p);
    for (let j = 0; j < p; j++) {
      for (let k = 0; k < p; k++) { let s = 0; for (let i = 0; i < n; i++) s += D[i][j] * ViD[k][i]; A[j][k] = s; }
      let s = 0; for (let i = 0; i < n; i++) s += D[i][j] * Viy[i];
      b[j] = s;
    }
    const LA = chol(A);
    if (!LA) return null;
    const beta = cholSolve(LA, b);
    /* the residual variance, by maximum likelihood as nlme::gls does with
       method = "ML" */
    const resid = Float64Array.from({ length: n }, (_, i) => {
      let f = 0;
      for (let j = 0; j < p; j++) f += D[i][j] * beta[j];
      return y[i] - f;
    });
    const Vir = cholSolve(L, resid);
    let q = 0;
    for (let i = 0; i < n; i++) q += resid[i] * Vir[i];
    const sigma2 = q / n;
    const lnL = -0.5 * (n * Math.log(2 * Math.PI) + n * Math.log(sigma2) + cholLogDet(L) + n);
    /* standard errors from (Dᵀ V⁻¹ D)⁻¹ σ², with the unbiased σ² nlme reports */
    const s2u = q / (n - p);
    const se = new Float64Array(p);
    for (let j = 0; j < p; j++) {
      const e = new Float64Array(p); e[j] = 1;
      const col = cholSolve(LA, e);
      se[j] = Math.sqrt(Math.max(0, col[j] * s2u));
    }
    const tvals = Float64Array.from(beta, (v, j) => (se[j] > 0 ? v / se[j] : 0));
    return {
      beta: Array.from(beta), se: Array.from(se), t: Array.from(tvals),
      p: Array.from(tvals, tv => 2 * (1 - studentCdf(Math.abs(tv), n - p))),
      sigma2, lnL, n, k: p + 1, df: n - p,
      AIC: -2 * lnL + 2 * (p + 1),
    };
  }
  /* Student's t distribution function, through the incomplete beta. The beta
     only sees t², so the sign has to be put back by hand — otherwise the
     function comes out symmetric about zero, which is what a distribution
     function must never be. */
  function studentCdf(t, df) {
    const half = 0.5 * incompleteBeta(df / (df + t * t), df / 2, 0.5);
    return t >= 0 ? 1 - half : half;
  }
  function incompleteBeta(x, a, b) {
    if (x <= 0) return 0;
    if (x >= 1) return 1;
    const lbeta = lgamma(a) + lgamma(b) - lgamma(a + b);
    const front = Math.exp(a * Math.log(x) + b * Math.log(1 - x) - lbeta) / a;
    let f = 1, c = 1, d = 0;
    for (let i = 0; i <= 300; i++) {
      const m = Math.floor(i / 2);
      let num;
      if (i === 0) num = 1;
      else if (i % 2 === 0) num = (m * (b - m) * x) / ((a + 2 * m - 1) * (a + 2 * m));
      else num = -((a + m) * (a + b + m) * x) / ((a + 2 * m) * (a + 2 * m + 1));
      d = 1 + num * d;
      if (Math.abs(d) < 1e-30) d = 1e-30;
      d = 1 / d;
      c = 1 + num / c;
      if (Math.abs(c) < 1e-30) c = 1e-30;
      const cd = c * d;
      f *= cd;
      if (Math.abs(1 - cd) < 1e-12) break;
    }
    const r = front * (f - 1);
    return x < (a + 1) / (a + b + 2) ? r : 1 - incompleteBeta(1 - x, b, a);
  }
  function lgamma(z) {
    const c = [76.18009172947146, -86.50532032941677, 24.01409824083091,
      -1.231739572450155, 0.1208650973866179e-2, -0.5395239384953e-5];
    let y = z, tmp = z + 5.5;
    tmp -= (z + 0.5) * Math.log(tmp);
    let ser = 1.000000000190015;
    for (let j = 0; j < 6; j++) ser += c[j] / ++y;
    return -tmp + Math.log(2.5066282746310005 * ser / z);
  }

  /* ================================================================
     6 · phylogenetic signal
     ================================================================ */
  /* Pagel's λ multiplies the off-diagonals of the covariance: λ = 1 is
     Brownian motion, λ = 0 is a star with no shared history at all, and the
     maximum-likelihood value says how much of the covariance the tree explains. */
  function pagelLambda(tree, x, opts) {
    opts = opts || {};
    const C = opts.C || vcv(tree, x.length);
    const obj = lam => {
      const f = gaussianProfile(covLambda(C, lam), x);
      return f ? -f.lnL : 1e100;
    };
    const best = g.Like.brentMin(obj, 0, 1, 1e-8, 60);
    const lambda = best.x;
    const f = gaussianProfile(covLambda(C, lambda), x);
    const f0 = gaussianProfile(covLambda(C, 0), x);
    const stat = 2 * (f.lnL - f0.lnL);
    return {
      lambda, lnL: f.lnL, lnL0: f0.lnL, sigma2: f.sigma2, z0: f.z0,
      stat, p: 1 - g.Dist.pchisq(Math.max(0, stat), 1),
    };
  }
  /* Blomberg's K compares the variance the trait has around the phylogenetic
     mean with what Brownian motion on this tree would give: K = 1 is exactly
     Brownian, below 1 less signal than the tree implies, above 1 more. */
  function blombergK(tree, x, opts) {
    opts = opts || {};
    const n = x.length;
    const C = opts.C || vcv(tree, n);
    const L = chol(C);
    if (!L) return null;
    const one = Float64Array.from({ length: n }, () => 1);
    const Vi1 = cholSolve(L, one), Vix = cholSolve(L, x);
    let a = 0, b = 0;
    for (let i = 0; i < n; i++) { a += one[i] * Vix[i]; b += one[i] * Vi1[i]; }
    const z0 = a / b;
    const d = Float64Array.from(x, (v, i) => v - z0);
    const Vid = cholSolve(L, d);
    let mse = 0;
    for (let i = 0; i < n; i++) mse += d[i] * Vid[i];
    mse /= (n - 1);
    let mse0 = 0;
    for (let i = 0; i < n; i++) mse0 += d[i] * d[i];
    mse0 /= (n - 1);
    /* the expectation of the ratio under Brownian motion */
    let trC = 0, sumVi = 0;
    for (let i = 0; i < n; i++) { trC += C[i][i]; sumVi += Vi1[i]; }
    const expected = (trC - n / sumVi) / (n - 1);
    const K = (mse0 / mse) / expected;
    return { K, mse, mse0, expected, z0 };
  }
  /* the randomisation test: shuffle the values among the tips and see how often
     a tree this shape produces as much signal by chance */
  function signalTest(tree, x, opts) {
    opts = opts || {};
    const reps = opts.reps || 999;
    const r = g.rng(opts.seed || 3);
    const C = opts.C || vcv(tree, x.length);
    const obs = blombergK(tree, x, { C }).K;
    let ge = 0;
    const nullK = [];
    for (let i = 0; i < reps; i++) {
      const y = Float64Array.from(x);
      for (let j = y.length - 1; j > 0; j--) { const k = Math.floor(r() * (j + 1)); const t = y[j]; y[j] = y[k]; y[k] = t; }
      const kk = blombergK(tree, y, { C });
      if (kk) { nullK.push(kk.K); if (kk.K >= obs) ge++; }
    }
    return { K: obs, p: (ge + 1) / (nullK.length + 1), nullDistribution: nullK, reps: nullK.length };
  }

  /* ================================================================
     7 · ancestral states of a continuous character
     ================================================================ */
  /* The value at every internal node under Brownian motion, with its standard
     error: the conditional expectation given the tips, which for a Gaussian is
     a weighted average of everything below and everything above. */
  function ancestralContinuous(tree, x, opts) {
    opts = opts || {};
    const F = g.Tree.flatten(tree);
    const n = x.length;
    /* the contrast algorithm gives the root; the rest come from rerooting at
       each node in turn, which is the identity Felsenstein used */
    const out = [];
    const fit = fitBM(tree, x, opts);
    for (let k = 0; k < F.n; k++) {
      if (F.isTip[k]) continue;
      const rooted = k === 0 ? g.Tree.clone(tree) : g.Tree.rerootAbove(tree, F.nodes[k], 0);
      const p = pic(rooted, x);
      out.push({ node: k, value: p.rootValue, variance: p.rootVariance * fit.sigma2 });
    }
    return { states: out, sigma2: fit.sigma2, lnL: fit.lnL, root: out.length ? out[0].value : null };
  }

  /* ================================================================
     8 · discrete characters: the Mk model
     ================================================================ */
  /* the rate matrix for a model: 'ER' one rate, 'SYM' one per unordered pair,
     'ARD' one per ordered pair */
  function rateIndex(k, model) {
    const idx = Array.from({ length: k }, () => new Int32Array(k).fill(-1));
    let n = 0;
    if (model === 'ER') {
      for (let i = 0; i < k; i++) for (let j = 0; j < k; j++) if (i !== j) idx[i][j] = 0;
      n = 1;
    } else if (model === 'SYM') {
      for (let i = 0; i < k; i++) for (let j = i + 1; j < k; j++) { idx[i][j] = n; idx[j][i] = n; n++; }
    } else {
      for (let i = 0; i < k; i++) for (let j = 0; j < k; j++) if (i !== j) { idx[i][j] = n; n++; }
    }
    return { idx, n };
  }
  function buildQ(rates, idx, k) {
    const Q = Array.from({ length: k }, () => new Float64Array(k));
    for (let i = 0; i < k; i++) {
      let s = 0;
      for (let j = 0; j < k; j++) if (i !== j) { Q[i][j] = rates[idx[i][j]]; s += Q[i][j]; }
      Q[i][i] = -s;
    }
    return Q;
  }

  /* the downward pass, in ape::ace's convention: every node's partials are
     normalised and the log-likelihood is the sum of the normalising constants,
     which amounts to summing over the root states with equal weight one */
  function mkLogLik(F, tips, Q, k) {
    const partial = new Array(F.n);
    let lnL = 0;
    for (let i = 0; i < F.post.length; i++) {
      const node = F.post[i];
      if (F.isTip[node]) { partial[node] = tips[F.tipRow[node]]; continue; }
      const v = new Float64Array(k).fill(1);
      F.kids[node].forEach(c => {
        const P = expm(Q, F.len[c] || 0);
        const cp = partial[c];
        for (let a = 0; a < k; a++) {
          let s = 0;
          for (let b = 0; b < k; b++) s += P[a][b] * cp[b];
          v[a] *= s;
        }
      });
      let tot = 0;
      for (let a = 0; a < k; a++) tot += v[a];
      if (!(tot > 0)) return { lnL: -Infinity, partial };
      for (let a = 0; a < k; a++) v[a] /= tot;
      partial[node] = v;
      lnL += Math.log(tot);
    }
    return { lnL, partial };
  }

  function fitMk(tree, states, opts) {
    opts = opts || {};
    const model = opts.model || 'ER';
    const levels = opts.levels || [...new Set(states)].sort();
    const k = levels.length;
    const F = g.Tree.flatten(tree);
    const tips = states.map(s => {
      const v = new Float64Array(k);
      const i = levels.indexOf(s);
      if (i >= 0) v[i] = 1; else v.fill(1);       // an unknown state is any state
      return v;
    });
    const { idx, n } = rateIndex(k, model);
    const obj = logRates => {
      const rates = logRates.map(Math.exp);
      const Q = buildQ(rates, idx, k);
      const r = mkLogLik(F, tips, Q, k);
      return isFinite(r.lnL) ? -r.lnL : 1e100;
    };
    const start = new Array(n).fill(Math.log(0.05));
    const best = n === 1
      ? (() => { const b = g.Like.brentMin(v => obj([v]), Math.log(1e-8), Math.log(50), 1e-9, 80); return { x: [b.x], f: b.f }; })()
      : g.Diversify.nelderMead(obj, start, { maxIter: 2000, tol: 1e-12 });
    const rates = best.x.map(Math.exp);
    const Q = buildQ(rates, idx, k);
    const r = mkLogLik(F, tips, Q, k);
    return {
      model, levels, k, rates, Q, lnL: r.lnL, nRates: n,
      AIC: -2 * r.lnL + 2 * n, partial: r.partial, F, tips,
    };
  }
  function compareMk(tree, states, opts) {
    const rows = ['ER', 'SYM', 'ARD'].map(m => fitMk(tree, states, Object.assign({}, opts, { model: m })));
    const uniq = [];
    rows.forEach(r => { if (!uniq.some(u => u.nRates === r.nRates)) uniq.push(r); });
    const use = rows;
    const best = use.reduce((a, b) => (b.AIC < a.AIC ? b : a), use[0]);
    let sw = 0;
    use.forEach(r => { r.dAIC = r.AIC - best.AIC; r.w = Math.exp(-0.5 * r.dAIC); sw += r.w; });
    use.forEach(r => { r.w /= sw; });
    use.sort((a, b) => a.AIC - b.AIC);
    return { rows: use, best: use[0] };
  }

  /* ancestral states.
     'downpass' is what ape::ace returns: the normalised conditional likelihood
     of the subtree below each node, which ignores everything above it.
     'marginal' adds the second pass and is the number phytools and corHMM
     report; it is the one to use, and the two are only equal at the root. */
  function ancestralDiscrete(fit, method) {
    const { F, tips, Q, k, partial } = fit;
    if (method === 'downpass') {
      const out = [];
      for (let i = 0; i < F.n; i++) if (!F.isTip[i]) out.push({ node: i, probs: Array.from(partial[i]) });
      out.sort((a, b) => a.node - b.node);
      return out;
    }
    /* the upward partials: everything except the subtree below the node */
    const up = new Array(F.n);
    up[0] = new Float64Array(k).fill(1);
    for (let i = F.post.length - 1; i >= 0; i--) {
      const node = F.post[i];
      if (F.isTip[node]) continue;
      const kids = F.kids[node];
      kids.forEach(c => {
        const v = Float64Array.from(up[node]);
        kids.forEach(o => {
          if (o === c) return;
          const P = expm(Q, F.len[o] || 0);
          const op = partial[o];
          for (let a = 0; a < k; a++) {
            let s = 0;
            for (let b = 0; b < k; b++) s += P[a][b] * op[b];
            v[a] *= s;
          }
        });
        /* carry it down the branch to the child */
        const Pc = expm(Q, F.len[c] || 0);
        const w = new Float64Array(k);
        for (let b = 0; b < k; b++) {
          let s = 0;
          for (let a = 0; a < k; a++) s += v[a] * Pc[a][b];
          w[b] = s;
        }
        let tot = 0; for (let b = 0; b < k; b++) tot += w[b];
        if (tot > 0) for (let b = 0; b < k; b++) w[b] /= tot;
        up[c] = w;
      });
    }
    const out = [];
    for (let i = 0; i < F.n; i++) {
      if (F.isTip[i]) continue;
      const v = new Float64Array(k);
      let tot = 0;
      for (let a = 0; a < k; a++) { v[a] = partial[i][a] * up[i][a]; tot += v[a]; }
      if (tot > 0) for (let a = 0; a < k; a++) v[a] /= tot;
      out.push({ node: i, probs: Array.from(v) });
    }
    return out;
  }

  /* ================================================================
     9 · stochastic character mapping
     ================================================================ */
  /* Huelsenbeck et al. (2003). A state is drawn for every node from the
     conditional distribution, top down, and then a history is simulated along
     each branch conditioned on its two ends. Doing that many times turns
     "probably state 0 here" into the time each lineage spent in each state,
     which is what a reader can actually use. */
  function simmap(fit, opts) {
    opts = opts || {};
    const reps = opts.reps || 100;
    const r = g.rng(opts.seed || 7);
    const { F, Q, k, partial } = fit;
    const marginal = ancestralDiscrete(fit, 'marginal');
    const byNode = new Map(marginal.map(m => [m.node, m.probs]));
    const timeInState = Array.from({ length: reps }, () => new Float64Array(k));
    const changes = new Float64Array(reps);
    const nodeStates = Array.from({ length: F.n }, () => new Float64Array(k));
    for (let rep = 0; rep < reps; rep++) {
      const st = new Int32Array(F.n).fill(-1);
      /* the root, from its marginal */
      st[0] = draw(byNode.get(0), r);
      nodeStates[0][st[0]]++;
      for (let i = F.post.length - 1; i >= 0; i--) {
        const node = F.post[i];
        if (node === 0) continue;
        const p = F.parent[node];
        if (st[p] < 0) continue;
        const P = expm(Q, F.len[node] || 0);
        /* the state at the child, given the parent and the data below */
        const below = partial[node] || unitVector(k, -1);
        const w = new Float64Array(k);
        let tot = 0;
        for (let b = 0; b < k; b++) { w[b] = P[st[p]][b] * (F.isTip[node] ? fit.tips[F.tipRow[node]][b] : below[b]); tot += w[b]; }
        if (!(tot > 0)) { st[node] = st[p]; continue; }
        for (let b = 0; b < k; b++) w[b] /= tot;
        st[node] = draw(w, r);
        if (!F.isTip[node]) nodeStates[node][st[node]]++;
        /* the history along the branch, by rejection until the ends match */
        const hist = branchHistory(Q, k, F.len[node] || 0, st[p], st[node], r, opts.maxTries || 200);
        hist.time.forEach((t, s) => { timeInState[rep][s] += t; });
        changes[rep] += hist.changes;
      }
    }
    const total = new Float64Array(k);
    timeInState.forEach(t => { for (let a = 0; a < k; a++) total[a] += t[a]; });
    let sum = 0; for (let a = 0; a < k; a++) sum += total[a];
    const nodeProbs = [];
    for (let i = 0; i < F.n; i++) {
      if (F.isTip[i]) continue;
      const v = Array.from(nodeStates[i], c => c / reps);
      nodeProbs.push({ node: i, probs: v });
    }
    return {
      reps, levels: fit.levels,
      timeInState: Array.from(total, v => v / reps),
      proportion: Array.from(total, v => (sum > 0 ? v / sum : 0)),
      changes: { mean: g.Mcmc.mean(Array.from(changes)), sd: g.Mcmc.sd(Array.from(changes)) },
      nodeProbs, perRep: timeInState.map(t => Array.from(t)),
    };
  }
  function unitVector(k, i) {
    const v = new Float64Array(k);
    if (i < 0) v.fill(1); else v[i] = 1;
    return v;
  }
  function draw(probs, r) {
    let u = r(), acc = 0;
    for (let i = 0; i < probs.length; i++) { acc += probs[i]; if (u <= acc) return i; }
    return probs.length - 1;
  }
  /* a continuous-time history along one branch with both ends fixed, by
     simulating forward and rejecting the runs that do not end where they must */
  function branchHistory(Q, k, len, from, to, r, maxTries) {
    for (let attempt = 0; attempt < maxTries; attempt++) {
      const time = new Float64Array(k);
      let t = 0, s = from, changes = 0;
      for (let guard = 0; guard < 10000; guard++) {
        const rate = -Q[s][s];
        if (!(rate > 0)) { time[s] += len - t; t = len; break; }
        const w = -Math.log(Math.max(1e-300, r())) / rate;
        if (t + w >= len) { time[s] += len - t; t = len; break; }
        time[s] += w;
        t += w;
        /* where it jumps to */
        let u = r() * rate, acc = 0, next = s;
        for (let b = 0; b < k; b++) {
          if (b === s) continue;
          acc += Q[s][b];
          if (u <= acc) { next = b; break; }
        }
        s = next;
        changes++;
      }
      if (s === to) return { time, changes };
    }
    /* the rejection failed: put the whole branch in the ancestral state, which
       is the conservative thing and is reported by the caller as a failure */
    const time = new Float64Array(k);
    time[from] = len;
    return { time, changes: from === to ? 0 : 1, failed: true };
  }

  Object.assign(Traits, {
    chol, cholSolve, cholLogDet, expm, vcv,
    covBM, covOU, covEB, covLambda, gaussianFit, gaussianProfile, isUltrametric,
    fitBM, fitOU, fitEB, compareContinuous,
    pic, pgls, studentCdf, incompleteBeta,
    pagelLambda, blombergK, signalTest, ancestralContinuous,
    rateIndex, buildQ, mkLogLik, fitMk, compareMk, ancestralDiscrete, simmap, branchHistory,
  });
  g.Traits = Traits;
}
TraitsCore(typeof window !== 'undefined' ? window : self);
