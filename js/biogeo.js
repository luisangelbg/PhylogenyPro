/* PhylogenyPro — molecular phylogenetics in the browser.
   Copyright (C) 2026  Luis Ángel Barrera-Guzmán

   This program is free software: you can redistribute it and/or modify it under
   the terms of the GNU General Public License as published by the Free Software
   Foundation, either version 3 of the License, or (at your option) any later
   version. It is distributed in the hope that it will be useful, but WITHOUT ANY
   WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
   A PARTICULAR PURPOSE. See the GNU General Public License, in the file LICENSE
   at the root of this program, or <https://www.gnu.org/licenses/>. */

/* PhylogenyPro — where the lineages were.

   Historical biogeography treats a species' geographic range as a character
   whose states are the subsets of a set of areas, and that is the whole reason
   it is harder than every other character in this app: with eight areas there
   are 255 states, and the model has to say not only how a range changes along a
   branch but what happens to it when a lineage splits.

   Along a branch (anagenesis), two things happen, each with one rate:
     d  dispersal — the range gains an area
     e  extinction — the range loses one, and a range of one area that loses it
        is gone, which is why the null range is a state and is absorbing

   At a node (cladogenesis), the ancestral range is divided between the two
   daughters, and which divisions are allowed is what separates the models:
     DEC (Ree & Smith 2008)  narrow sympatry, subset sympatry, and vicariance
                             in which one daughter keeps a single area
     DIVALIKE                vicariance only, with any split into two
                             complementary halves, of any size
     BAYAREALIKE             neither: both daughters inherit the whole range
     +J (Matzke 2013)        adds founder-event jump dispersal, one daughter
                             leaving for a single area outside the range

   The +J models fit better almost always, and that is exactly the problem: Ree
   & Sanmartín (2018) showed the comparison is not a fair one, because j buys
   probability at cladogenesis that d and e have to pay for along branches. The
   block computes them and says so.

   Wrapped in a named function so a Web Worker can be built from its own source
   text (js/pool.js). Nothing here touches the DOM. */

function BiogeoCore(g) {
const Biogeo = {};

  /* ================================================================
     1 · ranges as bit masks
     ================================================================ */
  const popcount = m => { let c = 0; while (m) { m &= m - 1; c++; } return c; };
  const areasOf = (m, n) => { const a = []; for (let i = 0; i < n; i++) if (m & (1 << i)) a.push(i); return a; };
  function rangeName(m, areaNames) {
    if (m === 0) return '—';
    return areasOf(m, areaNames.length).map(i => areaNames[i]).join('');
  }

  /* Every range the model is allowed to be in, smallest first, with the null
     range at the head when it is included. Capping the size is the only thing
     that makes more than six areas tractable, and it is a real assumption: a
     lineage that was once in five of six areas cannot be reconstructed. */
  function makeStates(nAreas, opts) {
    opts = opts || {};
    const maxAreas = Math.min(nAreas, opts.maxAreas || nAreas);
    const includeNull = opts.includeNull !== false;
    const out = [];
    if (includeNull) out.push(0);
    for (let m = 1; m < (1 << nAreas); m++) {
      const k = popcount(m);
      if (k <= maxAreas) out.push(m);
    }
    out.sort((a, b) => (popcount(a) - popcount(b)) || (a - b));
    const index = new Map(out.map((m, i) => [m, i]));
    return { list: out, index, nAreas, maxAreas, includeNull, n: out.length };
  }

  /* ================================================================
     2 · the rate matrix along a branch
     ================================================================ */
  /* m[i][k] is how easily a lineage in area i reaches area k, relative to one:
     that is where a land bridge that opened, or an ocean that did not close,
     enters the model. */
  function buildQ(states, d, e, opts) {
    opts = opts || {};
    const n = states.n, nA = states.nAreas;
    const mult = opts.multipliers || null;
    const allowed = opts.allowedAreas == null ? (1 << nA) - 1 : opts.allowedAreas;
    const Q = Array.from({ length: n }, () => new Float64Array(n));
    for (let s = 0; s < n; s++) {
      const A = states.list[s];
      if (A === 0) continue;                     // the null range is absorbing
      const size = popcount(A);
      let total = 0;
      /* dispersal: one new area at a time */
      if (size < states.maxAreas) {
        for (let k = 0; k < nA; k++) {
          if (A & (1 << k)) continue;
          if (!(allowed & (1 << k))) continue;
          let w = 0;
          for (let i = 0; i < nA; i++) if (A & (1 << i)) w += mult ? mult[i][k] : 1;
          if (w <= 0) continue;
          const t = states.index.get(A | (1 << k));
          if (t == null) continue;
          Q[s][t] += d * w;
          total += d * w;
        }
      }
      /* local extinction: one area at a time, and a single-area range that
         loses it becomes the null range */
      for (let i = 0; i < nA; i++) {
        if (!(A & (1 << i))) continue;
        const B = A & ~(1 << i);
        if (B === 0 && !states.includeNull) continue;
        const t = states.index.get(B);
        if (t == null) continue;
        Q[s][t] += e;
        total += e;
      }
      Q[s][s] = -total;
    }
    return Q;
  }

  /* The exponential of the rate matrix, by scaling and squaring. Q is not
     reversible, so the symmetrisation the substitution models use is not
     available and there is no shortcut.

     This function is called once per branch per likelihood evaluation, which
     over an optimisation is hundreds of thousands of times, so it allocates
     nothing: the three matrices it needs are kept in a scratch pad and reused.
     Writing it the obvious way, with a fresh matrix per Taylor term, spent more
     time in the garbage collector than in the arithmetic. */
  const scratch = new Map();
  function pad(n) {
    let p = scratch.get(n);
    if (!p) {
      p = { A: new Float64Array(n * n), B: new Float64Array(n * n), C: new Float64Array(n * n), n };
      scratch.set(n, p);
    }
    return p;
  }
  function expm(Q, t) {
    const n = Q.length;
    let norm = 0;
    for (let i = 0; i < n; i++) {
      let s = 0;
      const Qi = Q[i];
      for (let j = 0; j < n; j++) s += Math.abs(Qi[j]);
      if (s > norm) norm = s;
    }
    norm *= Math.abs(t);
    let sq = 0;
    while (norm > 0.5) { norm /= 2; sq++; }
    const h = t / Math.pow(2, sq);
    /* Q as one flat array, which is how the inner loop wants it */
    const p = pad(n);
    const q = p.C;
    for (let i = 0; i < n; i++) { const Qi = Q[i]; for (let j = 0; j < n; j++) q[i * n + j] = Qi[j]; }
    let out = p.A, term = p.B;
    out.fill(0); term.fill(0);
    for (let i = 0; i < n; i++) { out[i * n + i] = 1; term[i * n + i] = 1; }
    let tmp = p.D && p.D.length === n * n ? p.D : (p.D = new Float64Array(n * n));
    for (let k = 1; k <= 24; k++) {
      /* term <- term · Q · (h/k) */
      tmp.fill(0);
      for (let i = 0; i < n; i++) {
        const ro = i * n;
        for (let z = 0; z < n; z++) {
          const a = term[ro + z];
          if (a === 0) continue;
          const qo = z * n;
          for (let j = 0; j < n; j++) tmp[ro + j] += a * q[qo + j];
        }
      }
      const f = h / k;
      let mx = 0;
      for (let i = 0; i < n * n; i++) {
        const v = tmp[i] * f;
        term[i] = v;
        out[i] += v;
        const av = v < 0 ? -v : v;
        if (av > mx) mx = av;
      }
      if (mx < 1e-17) break;
    }
    /* the squarings, alternating between the two buffers */
    let src = out, dst = tmp;
    for (let s = 0; s < sq; s++) {
      dst.fill(0);
      for (let i = 0; i < n; i++) {
        const ro = i * n;
        for (let z = 0; z < n; z++) {
          const a = src[ro + z];
          if (a === 0) continue;
          const bo = z * n;
          for (let j = 0; j < n; j++) dst[ro + j] += a * src[bo + j];
        }
      }
      const sw = src; src = dst; dst = sw;
    }
    /* One flat array out, indexed P[a·n + b], with negatives clipped: rounding
       can leave a probability a hair below zero, and a negative probability
       propagates into a negative likelihood and then into the logarithm of a
       negative number, which is a much harder bug to find than this line is to
       write. Flat rather than an array of rows because this result is the one
       thing the function does allocate, and n small arrays per branch per
       likelihood evaluation is a great deal of garbage. */
    const P = new Float64Array(n * n);
    for (let i = 0; i < n * n; i++) { const v = src[i]; P[i] = v > 0 ? v : 0; }
    return P;
  }
  const flatEye = n => {
    const I = new Float64Array(n * n);
    for (let i = 0; i < n; i++) I[i * n + i] = 1;
    return I;
  };
  /* the product of two flat matrices, for a branch that crosses a stratum */
  function matmul(A, B, n) {
    n = n || Math.round(Math.sqrt(A.length));
    const C = new Float64Array(n * n);
    for (let i = 0; i < n; i++) {
      const ro = i * n;
      for (let z = 0; z < n; z++) {
        const a = A[ro + z];
        if (a === 0) continue;
        const bo = z * n;
        for (let j = 0; j < n; j++) C[ro + j] += a * B[bo + j];
      }
    }
    return C;
  }

  /* ================================================================
     3 · what happens at a node
     ================================================================ */
  /* The four models differ only here. The weights follow BioGeoBEARS: every
     event of a type carries the type's weight, the weights of one ancestral
     range are normalised to one, and j is taken out of the others' share so
     that adding jump dispersal does not simply add probability. */
  const MODELS = {
    DEC:          { y: 1, s: 1, v: 1, vAnySize: false, widespreadY: false, jShare: 3 },
    DIVALIKE:     { y: 1, s: 0, v: 1, vAnySize: true,  widespreadY: false, jShare: 2 },
    BAYAREALIKE:  { y: 1, s: 0, v: 0, vAnySize: false, widespreadY: true,  jShare: 1 },
  };

  /* the ordered daughter pairs an ancestral range can leave behind, with their
     probabilities, one list per ancestral state */
  function cladoEvents(states, model, j, opts) {
    opts = opts || {};
    const M = MODELS[model] || MODELS.DEC;
    j = j || 0;
    const share = (M.jShare - j) / M.jShare;
    const wy = M.y * share, ws = M.s * share, wv = M.v * share, wj = j;
    const nA = states.nAreas;
    const mult = opts.multipliers || null;
    const allowed = opts.allowedAreas == null ? (1 << nA) - 1 : opts.allowedAreas;
    const out = new Array(states.n);

    for (let s = 0; s < states.n; s++) {
      const A = states.list[s];
      out[s] = [];
      if (A === 0) continue;                     // a dead lineage does not speciate
      const size = popcount(A);
      const ev = [];

      /* sympatry. In DEC and DIVALIKE only a single-area ancestor copies itself
         into both daughters; in BAYAREALIKE any range does, which is what makes
         that model the one without cladogenetic change at all. */
      if (wy > 0 && (size === 1 || M.widespreadY)) ev.push([s, s, wy]);

      if (size > 1) {
        /* subset sympatry: one daughter keeps everything, the other one area */
        if (ws > 0) {
          for (let i = 0; i < nA; i++) {
            if (!(A & (1 << i))) continue;
            const t = states.index.get(1 << i);
            if (t == null) continue;
            ev.push([s, t, ws]);
            ev.push([t, s, ws]);
          }
        }
        /* vicariance: the range is cut in two complementary pieces. DEC allows
           only a cut that leaves one daughter with a single area; DIVALIKE
           allows any cut, with the sizes of the smaller half sharing the
           weight equally, which is BioGeoBEARS' maximum-entropy setting. */
        if (wv > 0) {
          const half = Math.floor(size / 2);
          const sizeWeight = new Float64Array(half + 1);
          if (M.vAnySize) { for (let k = 1; k <= half; k++) sizeWeight[k] = 1 / half; }
          else sizeWeight[1] = 1;
          const seen = new Set();
          for (let L = A; L > 0; L = (L - 1) & A) {          // every non-empty subset
            const R = A & ~L;
            if (R === 0) continue;
            const kL = popcount(L), kR = popcount(R);
            const smaller = Math.min(kL, kR);
            if (smaller > half || !(sizeWeight[smaller] > 0)) continue;
            const iL = states.index.get(L), iR = states.index.get(R);
            if (iL == null || iR == null) continue;
            const key = iL * states.n + iR;
            if (seen.has(key)) continue;
            seen.add(key);
            ev.push([iL, iR, wv * sizeWeight[smaller]]);
          }
        }
      }
      /* founder-event jump dispersal: one daughter stays, the other leaves for
         a single area the ancestor was never in */
      if (wj > 0) {
        for (let k = 0; k < nA; k++) {
          if (A & (1 << k)) continue;
          if (!(allowed & (1 << k))) continue;
          let rel = 0;
          if (mult) { for (let i = 0; i < nA; i++) if (A & (1 << i)) rel = Math.max(rel, mult[i][k]); }
          else rel = 1;
          if (!(rel > 0)) continue;
          const t = states.index.get(1 << k);
          if (t == null) continue;
          ev.push([s, t, wj * rel]);
          ev.push([t, s, wj * rel]);
        }
      }
      let tot = 0;
      ev.forEach(x => { tot += x[2]; });
      if (tot > 0) ev.forEach(x => { x[2] /= tot; });
      out[s] = ev;
    }
    return out;
  }

  /* ================================================================
     4 · time strata
     ================================================================ */
  /* Geography is not constant: an island rises, an isthmus closes. A stratum is
     a slice of time with its own dispersal multipliers and its own set of areas
     that exist at all, and a branch that crosses a boundary gets the product of
     the two matrices — in the right order, which is the part that is easy to
     get backwards. */
  function normaliseStrata(strata, rootDepth) {
    if (!strata || !strata.length) return null;
    const out = strata.map(s => Object.assign({}, s)).sort((a, b) => a.from - b.from);
    /* from and to are ages before the present, so stratum 0 touches the tips */
    out.forEach((s, i) => { if (s.to == null) s.to = i + 1 < out.length ? out[i + 1].from : rootDepth; });
    return out;
  }
  function strataFor(strata, ageTop, ageBottom) {
    /* the pieces of a branch that runs from ageBottom (older) to ageTop */
    const out = [];
    strata.forEach(s => {
      const lo = Math.max(ageTop, s.from), hi = Math.min(ageBottom, s.to);
      if (hi - lo > 1e-12) out.push({ stratum: s, length: hi - lo });
    });
    return out;
  }

  /* ================================================================
     5 · the likelihood
     ================================================================ */
  /* Felsenstein's pruning with one extra step: at every internal node the
     ancestral range is divided before the daughters start evolving, so the
     partials of the two daughters are combined through the cladogenesis matrix
     instead of being multiplied together. */
  function prepare(tree, tipRanges, states, opts) {
    opts = opts || {};
    const F = g.Tree.flatten(tree);
    const n = states.n;
    /* the age of every node, so the strata know where each branch runs */
    const depth = new Float64Array(F.n);
    for (let i = F.post.length - 1; i >= 0; i--) {
      const k = F.post[i], p = F.parent[k];
      if (p >= 0) depth[k] = depth[p] + (F.len[k] || 0);
    }
    let maxDepth = 0;
    for (let k = 0; k < F.n; k++) if (F.isTip[k]) maxDepth = Math.max(maxDepth, depth[k]);
    const age = Float64Array.from(depth, v => maxDepth - v);

    const tipPartial = new Array(F.n);
    const missing = [];
    for (let k = 0; k < F.n; k++) {
      if (!F.isTip[k]) continue;
      const v = new Float64Array(n);
      const mask = tipRanges[F.tipRow[k]];
      if (mask == null) { v.fill(1); if (states.includeNull) v[0] = 0; missing.push(F.tipRow[k]); }
      else {
        const i = states.index.get(mask);
        if (i == null) { v.fill(1); if (states.includeNull) v[0] = 0; missing.push(F.tipRow[k]); }
        else v[i] = 1;
      }
      tipPartial[k] = v;
    }
    return { F, depth, age, maxDepth, tipPartial, missing };
  }

  /* P(t) for one branch, through however many strata it crosses */
  function branchP(prep, k, Qcache, states, params, strata) {
    const F = prep.F;
    const len = F.len[k] || 0;
    if (len <= 0) return flatEye(states.n);
    if (!strata) return Qcache.get(len);
    const pieces = strataFor(strata, prep.age[k], prep.age[k] + len);
    if (pieces.length <= 1) {
      const st = pieces.length ? pieces[0].stratum : strata[0];
      return Qcache.get(len, st);
    }
    /* the branch runs from the past towards the present, so the matrices
       multiply in that order: the oldest stratum first */
    pieces.sort((a, b) => b.stratum.from - a.stratum.from);
    let P = null;
    pieces.forEach(p => {
      const M = Qcache.get(p.length, p.stratum);
      P = P ? matmul(P, M, states.n) : M;
    });
    return P;
  }
  /* P(t) for every branch of one tree, under one rate matrix.

     Computing each branch's exponential from scratch does the same work over
     and over: the matrix is the same, only the time differs. Uniformization
     (Jensen 1953; Grassmann 1977) does the shared part once. Write

         Q = μ (R − I)     with  μ = max |Q_ii|  and  R = I + Q/μ,

     where R is a proper transition matrix — non-negative, rows summing to one —
     and then

         P(t) = exp(Qt) = Σ_k  Poisson(k ; μt) · R^k.

     The powers of R are built once and shared by every branch; a branch then
     costs a weighted sum instead of twenty matrix products. Every term is
     non-negative, so nothing cancels and the result cannot come out negative,
     which the Taylor series could. On a five-area problem this is some thirty
     times faster than calling the exponential per branch. */
  function makeQCache(states, d, e, opts) {
    const n = states.n;
    const cache = new Map();
    const byStratum = new Map();
    const keyOf = stratum => (stratum ? (stratum.name || String(stratum.from)) : '');

    function unif(stratum) {
      const sk = keyOf(stratum);
      let u = byStratum.get(sk);
      if (u) return u;
      const Q = buildQ(states, d, e, Object.assign({}, opts, stratum ? {
        multipliers: stratum.multipliers || opts.multipliers,
        allowedAreas: stratum.allowedAreas != null ? stratum.allowedAreas : opts.allowedAreas,
      } : {}));
      let mu = 0;
      for (let i = 0; i < n; i++) { const v = -Q[i][i]; if (v > mu) mu = v; }
      const powers = [flatEye(n)];
      let R = null;
      if (mu > 0) {
        R = new Float64Array(n * n);
        for (let i = 0; i < n; i++) {
          const Qi = Q[i];
          for (let j = 0; j < n; j++) R[i * n + j] = (i === j ? 1 : 0) + Qi[j] / mu;
        }
      }
      u = { Q, mu, R, powers };
      byStratum.set(sk, u);
      return u;
    }
    function powerUpTo(u, K) {
      while (u.powers.length <= K) {
        u.powers.push(matmul(u.powers[u.powers.length - 1], u.R, n));
      }
    }
    return {
      get(len, stratum) {
        const key = keyOf(stratum) + '|' + len.toPrecision(12);
        const hit = cache.get(key);
        if (hit) return hit;
        const u = unif(stratum);
        if (!(u.mu > 0) || !(len > 0)) { const I = flatEye(n); cache.set(key, I); return I; }
        const mt = u.mu * len;
        /* enough Poisson terms to leave less than 1e-15 in the tail */
        const K = Math.max(4, Math.ceil(mt + 7 * Math.sqrt(mt) + 12));
        powerUpTo(u, K);
        const P = new Float64Array(n * n);
        let w = Math.exp(-mt), acc = 0;
        for (let k = 0; k <= K; k++) {
          if (w > 1e-18) {
            const Rk = u.powers[k];
            for (let i = 0; i < n * n; i++) P[i] += w * Rk[i];
          }
          acc += w;
          if (acc > 1 - 1e-15 && k > mt) break;
          w *= mt / (k + 1);
        }
        cache.set(key, P);
        return P;
      },
      matrixFor(stratum) { return unif(stratum).Q; },
    };
  }

  /* the log-likelihood, and on the way the partials every later step needs */
  function logLike(prep, states, params, opts) {
    opts = opts || {};
    const F = prep.F, n = states.n;
    const clado = opts.clado || cladoEvents(states, params.model || 'DEC', params.j || 0, opts);
    const Qcache = opts.Qcache || makeQCache(states, params.d, params.e, opts);
    const strata = opts.strata || null;

    const partial = new Array(F.n);      // at the node, after cladogenesis
    const atTop = new Array(F.n);        // at the top of the branch above it
    let lnL = 0;
    for (let i = 0; i < F.post.length; i++) {
      const k = F.post[i];
      if (F.isTip[k]) { partial[k] = prep.tipPartial[k]; continue; }
      const kids = F.kids[k];
      /* carry each daughter up its branch */
      const up = kids.map(c => {
        const P = branchP(prep, c, Qcache, states, params, strata);
        const cp = partial[c];
        const v = new Float64Array(n);
        for (let a = 0; a < n; a++) {
          let s = 0;
          for (let b = 0; b < n; b++) s += P[a * n + b] * cp[b];
          v[a] = s;
        }
        atTop[c] = v;
        return v;
      });
      /* and divide the ancestral range between them */
      const v = new Float64Array(n);
      if (kids.length === 2) {
        const L = up[0], R = up[1];
        for (let a = 0; a < n; a++) {
          const ev = clado[a];
          let s = 0;
          for (let z = 0; z < ev.length; z++) s += ev[z][2] * L[ev[z][0]] * R[ev[z][1]];
          v[a] = s;
        }
      } else {
        /* a polytomy is resolved as a ladder of two-way splits, which is what
           every program does and what the user has to be told */
        let cur = up[0];
        for (let q = 1; q < up.length; q++) {
          const R = up[q], w = new Float64Array(n);
          for (let a = 0; a < n; a++) {
            const ev = clado[a];
            let s = 0;
            for (let z = 0; z < ev.length; z++) s += ev[z][2] * cur[ev[z][0]] * R[ev[z][1]];
            w[a] = s;
          }
          cur = w;
        }
        for (let a = 0; a < n; a++) v[a] = cur[a];
      }
      let tot = 0;
      for (let a = 0; a < n; a++) tot += v[a];
      if (!(tot > 0)) return { lnL: -Infinity, partial, atTop, clado, Qcache };
      for (let a = 0; a < n; a++) v[a] /= tot;
      partial[k] = v;
      lnL += Math.log(tot);
    }
    /* The root. The null range is excluded either way, because a clade cannot
       start out extinct. What is left is a choice of convention:
         'sum'  every allowed range weighted one, which is what LAGRANGE and
                BioGeoBEARS do and therefore what every published DEC
                log-likelihood means;
         'flat' a proper prior, one over the number of ranges.
       They differ by exactly log(number of non-null states), the same constant
       for every model over the same areas, so no model comparison changes — but
       the number a reader compares with a paper does, which is why the app
       reports BioGeoBEARS' convention by default and says which it used. */
    const root = partial[0];
    let s = 0, m = 0;
    for (let a = 0; a < n; a++) { if (states.list[a] === 0) continue; s += root[a]; m++; }
    const prior = opts.rootPrior || 'sum';
    lnL += Math.log(Math.max(1e-300, prior === 'flat' ? s / m : s));
    return { lnL, partial, atTop, clado, Qcache, rootPrior: prior, nRootStates: m };
  }

  /* ================================================================
     6 · fitting
     ================================================================ */
  function fit(tree, tipRanges, states, opts) {
    opts = opts || {};
    const model = opts.model || 'DEC';
    const withJ = !!opts.withJ;
    const prep = opts.prep || prepare(tree, tipRanges, states, opts);
    const jShare = (MODELS[model] || MODELS.DEC).jShare;
    const jMax = jShare - 1e-5;
    const nParams = withJ ? 3 : 2;

    let evals = 0;
    const obj = p => {
      evals++;
      const d = Math.exp(p[0]), e = Math.exp(p[1]);
      let j = 0;
      if (withJ) {
        j = p[2];
        if (!(j >= 0) || j > jMax) return 1e100;
      }
      if (!(d > 0) || !(e >= 0) || d > 20 || e > 20) return 1e100;
      const r = logLike(prep, states, { d, e, j, model }, opts);
      return isFinite(r.lnL) ? -r.lnL : 1e100;
    };
    /* Where to start. d and e are rates per unit of branch length, so a fixed
       number like 0.03 is a good guess on a tree five million years deep and a
       terrible one on a phylogram a tenth of a substitution deep — and a bad
       start is not a wrong answer, it is a hundred extra iterations. Scaling by
       the depth of the tree makes the first guess right to within a factor of a
       few on both. A few starts remain, because the surface with j in it has a
       ridge along which d and j trade against each other. */
    const sc = 1 / Math.max(1e-9, prep.maxDepth);
    const L = v => Math.log(v * sc);
    const starts = withJ
      ? [[L(0.2), L(0.2), 0.01], [L(0.05), L(0.05), 0.2], [L(1), L(1), 0.5], [L(0.02), L(0.02), 0.05]]
      : [[L(0.2), L(0.2)], [L(0.05), L(0.05)], [L(1), L(1)]];
    let best = null, agreed = 0;
    for (let q = 0; q < starts.length; q++) {
      const r = g.Diversify.nelderMead(obj, starts[q], { maxIter: 3000, tol: 1e-10 });
      if (best && Math.abs(r.f - best.f) < 1e-6) agreed++;
      if (!best || r.f < best.f) best = r;
      /* two starts that land on the same likelihood have said what they had to
         say; a third only costs time */
      if (agreed >= 1 && q >= 1) break;
    }
    const d = Math.exp(best.x[0]), e = Math.exp(best.x[1]);
    const j = withJ ? Math.max(0, Math.min(jMax, best.x[2])) : 0;
    const res = logLike(prep, states, { d, e, j, model }, opts);
    const k = nParams;
    const nTips = g.Tree.tips(tree).length;
    return {
      model, withJ, name: model + (withJ ? '+J' : ''),
      d, e, j, lnL: res.lnL, k,
      AIC: -2 * res.lnL + 2 * k,
      AICc: -2 * res.lnL + 2 * k + (2 * k * (k + 1)) / Math.max(1, nTips - k - 1),
      prep, states, partial: res.partial, atTop: res.atTop, clado: res.clado,
      Qcache: res.Qcache, evals, jMax,
    };
  }

  function compare(tree, tipRanges, states, opts) {
    opts = opts || {};
    const which = opts.models || [
      ['DEC', false], ['DEC', true], ['DIVALIKE', false], ['DIVALIKE', true],
      ['BAYAREALIKE', false], ['BAYAREALIKE', true],
    ];
    const prep = prepare(tree, tipRanges, states, opts);
    const rows = which.map(([m, wj]) => fit(tree, tipRanges, states,
      Object.assign({}, opts, { model: m, withJ: wj, prep })));
    const best = rows.reduce((a, b) => (b.AICc < a.AICc ? b : a), rows[0]);
    let sw = 0;
    rows.forEach(r => { r.dAICc = r.AICc - best.AICc; r.w = Math.exp(-0.5 * r.dAICc); sw += r.w; });
    rows.forEach(r => { r.w /= sw; });
    /* the likelihood ratio of each +J against the model it contains: one degree
       of freedom, and j sits on the boundary when it is zero, so the naive
       chi-square is conservative and the block says so */
    rows.forEach(r => {
      if (!r.withJ) return;
      const base = rows.find(q => q.model === r.model && !q.withJ);
      if (!base) return;
      const stat = 2 * (r.lnL - base.lnL);
      r.lrt = { stat, p: stat <= 0 ? 1 : 1 - g.Dist.pchisq(stat, 1), against: base.name };
    });
    const sorted = rows.slice().sort((a, b) => a.AICc - b.AICc);
    return { rows: sorted, best: sorted[0], prep };
  }

  /* ================================================================
     7 · ancestral ranges
     ================================================================ */
  /* The marginal probability of every range at every node, which needs the
     second pass from the root down — and, because a range is divided at the
     node, the probability of each daughter's range immediately after the split.
     Those corners are what a biogeographic figure actually shows. */
  function ancestral(fitted, opts) {
    opts = opts || {};
    const { prep, states, partial, atTop, clado, Qcache } = fitted;
    const F = prep.F, n = states.n;
    const strata = opts.strata || null;
    const params = { d: fitted.d, e: fitted.e, j: fitted.j, model: fitted.model };

    /* the root's marginal is its own partial against the prior */
    const up = new Array(F.n);
    const r0 = new Float64Array(n);
    let tot0 = 0;
    for (let a = 0; a < n; a++) {
      if (states.list[a] === 0) continue;
      r0[a] = 1;
      tot0 += 1;
    }
    for (let a = 0; a < n; a++) r0[a] /= tot0;
    up[0] = r0;

    const marginal = new Array(F.n);
    const corners = new Array(F.n);
    /* preorder: the root first */
    const order = F.post.slice().reverse();
    order.forEach(k => {
      const v = new Float64Array(n);
      let tot = 0;
      for (let a = 0; a < n; a++) { v[a] = up[k][a] * partial[k][a]; tot += v[a]; }
      if (tot > 0) for (let a = 0; a < n; a++) v[a] /= tot;
      marginal[k] = v;
      if (F.isTip[k]) return;
      const kids = F.kids[k];
      if (kids.length !== 2) {
        /* a polytomy has no single division to report */
        kids.forEach(c => { up[c] = Float64Array.from(atTop[c], (x, a) => x); normalise(up[c]); });
        return;
      }
      const [c1, c2] = kids;
      const L = atTop[c1], R = atTop[c2];

      /* What each corner ended up being, given everything: the posterior over
         the divisions. This is for reading, so it uses all the data. */
      const cornerL = new Float64Array(n), cornerR = new Float64Array(n);
      let tt = 0;
      for (let a = 0; a < n; a++) {
        const pa = up[k][a] * partial[k][a];
        if (!(pa > 0)) continue;
        const ev = clado[a];
        let denom = 0;
        for (let z = 0; z < ev.length; z++) denom += ev[z][2] * L[ev[z][0]] * R[ev[z][1]];
        if (!(denom > 0)) continue;
        for (let z = 0; z < ev.length; z++) {
          const w = pa * ev[z][2] * L[ev[z][0]] * R[ev[z][1]] / denom;
          cornerL[ev[z][0]] += w;
          cornerR[ev[z][1]] += w;
          tt += w;
        }
      }
      if (tt > 0) { for (let a = 0; a < n; a++) { cornerL[a] /= tt; cornerR[a] /= tt; } }
      corners[k] = { left: c1, right: c2, leftProbs: Array.from(cornerL), rightProbs: Array.from(cornerR) };

      /* The message a daughter receives is a different quantity, and confusing
         the two is the classic way to get ancestral states that look far more
         certain than they are: what comes down to a daughter must carry
         everything *except* that daughter's own subtree. So it is built from
         the node's own outside message and the sister's likelihood, and the
         daughter's own likelihood is left out. */
      const msgL = new Float64Array(n), msgR = new Float64Array(n);
      for (let a = 0; a < n; a++) {
        const ua = up[k][a];
        if (!(ua > 0)) continue;
        const ev = clado[a];
        for (let z = 0; z < ev.length; z++) {
          msgL[ev[z][0]] += ua * ev[z][2] * R[ev[z][1]];
          msgR[ev[z][1]] += ua * ev[z][2] * L[ev[z][0]];
        }
      }
      normalise(msgL); normalise(msgR);
      [[c1, msgL], [c2, msgR]].forEach(([c, msg]) => {
        const P = branchP(prep, c, Qcache, states, params, strata);
        const w = new Float64Array(n);
        for (let b = 0; b < n; b++) {
          let s = 0;
          for (let a = 0; a < n; a++) s += msg[a] * P[a * n + b];
          w[b] = s;
        }
        normalise(w);
        up[c] = w;
      });
    });
    const nodes = [];
    for (let k = 0; k < F.n; k++) {
      if (F.isTip[k]) continue;
      nodes.push({ node: k, probs: Array.from(marginal[k]), corner: corners[k] || null });
    }
    return { nodes, marginal, corners, states };
  }
  function normalise(v) {
    let t = 0;
    for (let i = 0; i < v.length; i++) t += v[i];
    if (t > 0) for (let i = 0; i < v.length; i++) v[i] /= t;
  }

  /* the readable summary: the most probable range at each node and how sure it is */
  function bestRanges(anc, areaNames, prep) {
    return anc.nodes.map(nd => {
      const i = nd.probs.indexOf(Math.max.apply(null, nd.probs));
      const ordered = nd.probs.map((p, k) => [k, p]).sort((a, b) => b[1] - a[1]);
      return {
        node: nd.node,
        range: rangeName(anc.states.list[i], areaNames),
        probability: nd.probs[i],
        second: ordered[1] ? rangeName(anc.states.list[ordered[1][0]], areaNames) : null,
        secondProbability: ordered[1] ? ordered[1][1] : 0,
        nTips: prep ? countTips(prep.F, nd.node) : null,
      };
    });
  }
  function countTips(F, k) {
    let c = 0;
    (function walk(x) { if (F.isTip[x]) { c++; return; } F.kids[x].forEach(walk); })(k);
    return c;
  }

  /* ================================================================
     8 · biogeographic stochastic mapping
     ================================================================ */
  /* The ancestral probabilities say where a lineage probably was; they do not
     say how many times anything happened. Drawing whole histories does — a
     range for every node and every corner, and then the dispersals and
     extinctions along each branch — and repeating it turns the reconstruction
     into counts of events, which is what a biogeographic story is made of. */
  function stochasticMap(fitted, opts) {
    opts = opts || {};
    const reps = opts.reps || 200;
    const r = g.rng(opts.seed || 31);
    const { prep, states, partial, atTop, clado, Qcache } = fitted;
    const F = prep.F, n = states.n;
    const params = { d: fitted.d, e: fitted.e, j: fitted.j, model: fitted.model };
    const strata = opts.strata || null;
    const anc = opts.ancestral || ancestral(fitted, opts);
    const Q = Qcache.matrixFor(strata ? strata[0] : null);

    const counts = { dispersal: 0, extinction: 0, sympatry: 0, subset: 0, vicariance: 0, jump: 0 };
    const perRep = [];
    const nodeCount = Array.from({ length: F.n }, () => new Float64Array(n));
    let failures = 0;

    for (let rep = 0; rep < reps; rep++) {
      const c = { dispersal: 0, extinction: 0, sympatry: 0, subset: 0, vicariance: 0, jump: 0 };
      const stateAt = new Int32Array(F.n).fill(-1);
      const order = F.post.slice().reverse();
      /* the root, from its marginal */
      stateAt[0] = draw(anc.marginal[0], r);
      nodeCount[0][stateAt[0]]++;
      for (let z = 0; z < order.length; z++) {
        const k = order[z];
        if (F.isTip[k] || stateAt[k] < 0) continue;
        const kids = F.kids[k];
        if (kids.length !== 2) continue;
        const A = stateAt[k];
        const ev = clado[A];
        const L = atTop[kids[0]], R = atTop[kids[1]];
        const w = new Float64Array(ev.length);
        let tot = 0;
        for (let q = 0; q < ev.length; q++) { w[q] = ev[q][2] * L[ev[q][0]] * R[ev[q][1]]; tot += w[q]; }
        if (!(tot > 0)) { failures++; continue; }
        for (let q = 0; q < ev.length; q++) w[q] /= tot;
        const pick = ev[draw(w, r)];
        c[classify(states, A, pick[0], pick[1])]++;
        /* each daughter starts its branch in its corner and ends where the data
           allow: a history along the branch conditioned on both ends */
        [[kids[0], pick[0]], [kids[1], pick[1]]].forEach(([child, startState]) => {
          const P = branchP(prep, child, Qcache, states, params, strata);
          const below = partial[child];
          const wv = new Float64Array(n);
          let t2 = 0;
          for (let b = 0; b < n; b++) { wv[b] = P[startState * n + b] * below[b]; t2 += wv[b]; }
          if (!(t2 > 0)) { stateAt[child] = startState; failures++; return; }
          for (let b = 0; b < n; b++) wv[b] /= t2;
          const endState = draw(wv, r);
          stateAt[child] = endState;
          if (!F.isTip[child]) nodeCount[child][endState]++;
          const h = branchHistory(Q, states, F.len[child] || 0, startState, endState, r, opts.maxTries || 100);
          c.dispersal += h.dispersal;
          c.extinction += h.extinction;
          if (h.failed) failures++;
        });
      }
      Object.keys(c).forEach(key => { counts[key] += c[key]; });
      perRep.push(c);
    }
    const mean = {}, sd = {};
    Object.keys(counts).forEach(key => {
      mean[key] = counts[key] / reps;
      const vals = perRep.map(p => p[key]);
      sd[key] = g.Mcmc.sd(vals);
    });
    const nodeProbs = [];
    for (let k = 0; k < F.n; k++) {
      if (F.isTip[k]) continue;
      nodeProbs.push({ node: k, probs: Array.from(nodeCount[k], v => v / reps) });
    }
    return { reps, mean, sd, perRep, nodeProbs, failures, states };
  }
  /* which of the four cladogenetic events a division was */
  function classify(states, A, L, R) {
    const a = states.list[A], l = states.list[L], rr = states.list[R];
    if (l === a && rr === a) return 'sympatry';
    if ((l === a && (rr & a) === 0) || (rr === a && (l & a) === 0)) return 'jump';
    if (l === a || rr === a) return 'subset';
    return 'vicariance';
  }
  /* A history along one branch with both ends fixed. Plain rejection sampling
     wastes most of its draws when the two ends differ, because the commonest
     forward history is the one where nothing happens at all: Nielsen's (2002)
     modification forces the first jump to fall inside the branch in that case,
     by drawing the waiting time from the exponential truncated at the branch
     length. It is the same distribution, sampled where it is not wasted. */
  function branchHistory(Q, states, len, from, to, r, maxTries) {
    if (len <= 0) return { dispersal: 0, extinction: 0 };
    const mustMove = from !== to;
    for (let attempt = 0; attempt < maxTries; attempt++) {
      let t = 0, s = from, disp = 0, ext = 0;
      for (let guard = 0; guard < 5000; guard++) {
        const rate = -Q[s][s];
        if (!(rate > 0)) { t = len; break; }
        let w;
        if (guard === 0 && mustMove) {
          const u = r();
          const cap = 1 - Math.exp(-rate * len);
          w = -Math.log(Math.max(1e-300, 1 - u * cap)) / rate;
        } else {
          w = -Math.log(Math.max(1e-300, r())) / rate;
        }
        if (t + w >= len) { t = len; break; }
        t += w;
        let u = r() * rate, acc = 0, next = s;
        for (let b = 0; b < states.n; b++) {
          if (b === s) continue;
          acc += Q[s][b];
          if (u <= acc) { next = b; break; }
        }
        if (popcount(states.list[next]) > popcount(states.list[s])) disp++; else ext++;
        s = next;
      }
      if (s === to) return { dispersal: disp, extinction: ext };
    }
    return { dispersal: 0, extinction: 0, failed: true };
  }
  function draw(probs, r) {
    let u = r(), acc = 0;
    for (let i = 0; i < probs.length; i++) { acc += probs[i]; if (u <= acc) return i; }
    return probs.length - 1;
  }

  /* ================================================================
     9 · reading a range table
     ================================================================ */
  /* Two formats, because both are in use: the LAGRANGE/PHYLIP one that
     BioGeoBEARS reads (a header line, then a name and a string of 0 and 1) and
     a plain table of names against area letters. */
  function parseRanges(text, opts) {
    opts = opts || {};
    const lines = String(text).replace(/\r/g, '').split('\n').map(l => l.trim()).filter(l => l && !/^#/.test(l));
    if (!lines.length) return null;
    /* the LAGRANGE header: "19 4 (K O M H)" */
    const head = lines[0].match(/^(\d+)\s+(\d+)\s*(?:\(([^)]*)\))?$/);
    if (head) {
      const nTax = +head[1], nAreas = +head[2];
      let names = head[3] ? head[3].trim().split(/\s+/) : null;
      if (!names || names.length !== nAreas) names = Array.from({ length: nAreas }, (_, i) => String.fromCharCode(65 + i));
      const taxa = [], masks = [];
      lines.slice(1).forEach(l => {
        const m = l.match(/^(\S+)\s+([01?]+)$/);
        if (!m) return;
        taxa.push(m[1]);
        let mask = 0;
        for (let i = 0; i < Math.min(nAreas, m[2].length); i++) if (m[2][i] === '1') mask |= (1 << i);
        masks.push(mask);
      });
      return { format: 'lagrange', areas: names, taxa, masks, expected: nTax };
    }
    /* a table: the first line names the areas */
    const sep = lines[0].indexOf('\t') >= 0 ? '\t' : lines[0].indexOf(';') >= 0 ? ';' : ',';
    const cols = lines[0].split(sep).map(s => s.trim().replace(/^["']|["']$/g, ''));
    const areas = cols.slice(1);
    const taxa = [], masks = [];
    lines.slice(1).forEach(l => {
      const c = l.split(sep).map(s => s.trim().replace(/^["']|["']$/g, ''));
      taxa.push(c[0]);
      let mask = 0;
      for (let i = 0; i < areas.length; i++) {
        const v = c[i + 1];
        if (v === '1' || /^(s[ií]|yes|y|x|true|present)$/i.test(v || '')) mask |= (1 << i);
      }
      masks.push(mask);
    });
    return { format: 'table', areas, taxa, masks };
  }
  function writeLagrange(areas, taxa, masks) {
    const out = [`${taxa.length}\t${areas.length}\t(${areas.join(' ')})`];
    taxa.forEach((t, i) => {
      let s = '';
      for (let k = 0; k < areas.length; k++) s += (masks[i] & (1 << k)) ? '1' : '0';
      out.push(`${t}\t${s}`);
    });
    return out.join('\n');
  }

  Object.assign(Biogeo, {
    popcount, areasOf, rangeName, makeStates, buildQ, expm,
    MODELS, cladoEvents, normaliseStrata, strataFor,
    prepare, logLike, makeQCache, branchP, fit, compare,
    ancestral, bestRanges, stochasticMap, branchHistory, classify,
    parseRanges, writeLagrange,
  });
  g.Biogeo = Biogeo;
}
BiogeoCore(typeof window !== 'undefined' ? window : self);
