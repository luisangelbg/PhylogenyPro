/* PhylogenyPro — molecular phylogenetics in the browser.
   Copyright (C) 2026  Luis Ángel Barrera-Guzmán

   This program is free software: you can redistribute it and/or modify it under
   the terms of the GNU General Public License as published by the Free Software
   Foundation, either version 3 of the License, or (at your option) any later
   version. It is distributed in the hope that it will be useful, but WITHOUT ANY
   WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
   A PARTICULAR PURPOSE. See the GNU General Public License, in the file LICENSE
   at the root of this program, or <https://www.gnu.org/licenses/>. */

/* PhylogenyPro — the likelihood engine.

   One engine for every kind of data: four states for DNA, twenty for proteins,
   k for morphology (the Mk model of Lewis 2001), sixty-one for codons. The
   pieces are the classical ones:

     · site patterns compressed once, with their weights (a 10 000-site
       alignment usually has a few thousand distinct patterns);
     · ambiguities as partial likelihoods — an R contributes to A and to G, a
       gap to everything — which is what lets real data with missing ends be
       analysed without throwing columns away;
     · the rate matrix Q built from a symmetric exchangeability matrix and the
       equilibrium frequencies, normalised to one substitution per site, and
       diagonalised through its symmetrised form (Jacobi), so P(t) = U e^{Λt} U⁻¹
       costs one exponential per eigenvalue and not a matrix exponential;
     · Felsenstein's (1981) pruning with per-node scaling, without which
       anything above ~50 taxa underflows to zero;
     · rate variation among sites with the discrete gamma of Yang (1994) and a
       proportion of invariable sites, with the rates rescaled so the mean rate
       over all sites stays one (the convention of PAML and IQ-TREE);
     · incremental recomputation: changing one branch only recomputes the path
       from that branch to the root, which is what makes NNI, SPR and MCMC
       affordable.

   Validated against phangorn in validation/block3. */

/* Wrapped in a named function so that a Web Worker can be built from its own
   source text: see js/pool.js. Nothing here touches the DOM. */
function LikeCore(g) {
const Like = {};

  const DNA = 'ACGT';
  /* IUPAC ambiguity codes: which bases each letter allows */
  const DNA_AMBIG = {
    A: 'A', C: 'C', G: 'G', T: 'T', U: 'T',
    R: 'AG', Y: 'CT', M: 'AC', K: 'GT', S: 'CG', W: 'AT',
    B: 'CGT', D: 'AGT', H: 'ACT', V: 'ACG', N: 'ACGT', X: 'ACGT', '?': 'ACGT', '-': 'ACGT',
  };
  const AA = () => (g.AA_ORDER || 'ARNDCQEGHILKMFPSTWYV');
  const AA_AMBIG = { B: 'DN', Z: 'EQ', J: 'IL', X: null, '?': null, '-': null, '*': null };

  /* ================================================================
     1 · compressing the alignment into patterns
     ================================================================ */
  /* type: 'dna' | 'aa' | 'morph'. For morphology the alphabet is whatever
     states actually occur, and {01} style polymorphisms are read as partials. */
  function compress(seqs, type, opts) {
    opts = opts || {};
    const n = seqs.length, L = Math.min(...seqs.map(s => s.length));
    let alphabet;
    if (type === 'aa') alphabet = AA().split('');
    else if (type === 'morph') {
      const set = new Set();
      seqs.forEach(s => { for (let i = 0; i < s.length; i++) { const c = s[i]; if (/[0-9A-Za-z]/.test(c)) set.add(c); } });
      alphabet = [...set].sort();
      if (alphabet.length < 2) alphabet = ['0', '1'];
    } else alphabet = DNA.split('');
    const S = alphabet.length;
    const pos = {};
    alphabet.forEach((a, i) => { pos[a] = i; });

    /* which states each character allows, as a bit mask */
    function maskOf(ch) {
      if (type === 'dna') {
        const allowed = DNA_AMBIG[ch];
        if (allowed == null) return (1 << S) - 1;
        let m = 0;
        for (let i = 0; i < allowed.length; i++) m |= 1 << pos[allowed[i]];
        return m;
      }
      if (type === 'aa') {
        if (pos[ch] != null) return 1 << pos[ch];
        const amb = AA_AMBIG[ch];
        if (amb) { let m = 0; for (let i = 0; i < amb.length; i++) if (pos[amb[i]] != null) m |= 1 << pos[amb[i]]; return m; }
        return (1 << S) - 1;
      }
      if (pos[ch] != null) return 1 << pos[ch];
      return (1 << S) - 1;                       // ? and - are "any state"
    }

    /* one key per column, then unique */
    const map = new Map(), patterns = [], weights = [], patternOfSite = new Int32Array(L);
    for (let c = 0; c < L; c++) {
      let key = '';
      for (let i = 0; i < n; i++) key += seqs[i][c];
      let p = map.get(key);
      if (p == null) {
        p = patterns.length;
        map.set(key, p);
        patterns.push(key);
        weights.push(0);
      }
      weights[p]++;
      patternOfSite[c] = p;
    }
    const nPat = patterns.length;

    /* tip partials: for every sequence, nPat × S */
    const tips = [];
    for (let i = 0; i < n; i++) {
      const arr = new Float64Array(nPat * S);
      for (let p = 0; p < nPat; p++) {
        const m = maskOf(patterns[p][i]);
        for (let s = 0; s < S; s++) if (m & (1 << s)) arr[p * S + s] = 1;
      }
      tips.push(arr);
    }

    /* empirical frequencies, and which patterns could be constant (for +I) */
    const counts = new Float64Array(S);
    let tot = 0;
    for (let p = 0; p < nPat; p++) {
      for (let i = 0; i < n; i++) {
        const ch = patterns[p][i];
        const k = pos[ch];
        if (k != null) { counts[k] += weights[p]; tot += weights[p]; }
      }
    }
    const freqs = new Float64Array(S);
    for (let s = 0; s < S; s++) freqs[s] = tot ? counts[s] / tot : 1 / S;
    /* a flat prior smoothing keeps a zero frequency from making Q singular */
    for (let s = 0; s < S; s++) if (freqs[s] <= 0) freqs[s] = 1e-6;
    let fs = 0; for (let s = 0; s < S; s++) fs += freqs[s];
    for (let s = 0; s < S; s++) freqs[s] /= fs;

    const constMask = new Int32Array(nPat);      // bit set = pattern can be constant in that state
    let nConstant = 0, nVariable = 0;
    for (let p = 0; p < nPat; p++) {
      let m = (1 << S) - 1;
      for (let i = 0; i < n; i++) m &= maskOf(patterns[p][i]);
      constMask[p] = m;
      if (m) nConstant += weights[p]; else nVariable += weights[p];
    }
    return {
      nStates: S, alphabet, nPat, nSites: L, nSeq: n,
      weights: Float64Array.from(weights), tips, patterns, patternOfSite,
      freqs, constMask, nConstant, nVariable,
    };
  }

  /* ================================================================
     2 · rate matrices and their eigen-decomposition
     ================================================================ */
  /* DNA rate schemes, as the index of the free rate each of the six
     exchangeabilities uses (order AC, AG, AT, CG, CT, GT) */
  const DNA_SCHEMES = {
    JC: [0, 0, 0, 0, 0, 0], F81: [0, 0, 0, 0, 0, 0],
    K80: [0, 1, 0, 0, 1, 0], HKY: [0, 1, 0, 0, 1, 0],
    TrNef: [0, 1, 0, 0, 2, 0], TrN: [0, 1, 0, 0, 2, 0],
    K81: [0, 1, 2, 2, 1, 0], K81uf: [0, 1, 2, 2, 1, 0],
    TIMef: [0, 1, 2, 2, 3, 0], TIM: [0, 1, 2, 2, 3, 0],
    TVMef: [0, 1, 2, 3, 1, 4], TVM: [0, 1, 2, 3, 1, 4],
    SYM: [0, 1, 2, 3, 4, 5], GTR: [0, 1, 2, 3, 4, 5],
  };
  const DNA_EQUAL_FREQ = { JC: 1, K80: 1, TrNef: 1, K81: 1, TIMef: 1, TVMef: 1, SYM: 1 };
  const PAIRS = [[0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3]];   // AC AG AT CG CT GT

  /* builds the symmetric exchangeability matrix S (S[i][j] = S[j][i]) */
  function exchangeability(spec, S) {
    const n = S;
    const M = Array.from({ length: n }, () => new Float64Array(n));
    if (spec.type === 'aa' && spec.model && g.AAMODELS && g.AAMODELS[spec.model]) {
      const q = g.AAMODELS[spec.model].Q;
      let k = 0;
      for (let i = 1; i < n; i++) for (let j = 0; j < i; j++) { M[i][j] = M[j][i] = q[k++]; }
      return M;
    }
    if (spec.type === 'morph' || spec.type === 'mk') {
      for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) if (i !== j) M[i][j] = 1;
      return M;
    }
    /* DNA */
    const scheme = DNA_SCHEMES[spec.model] || DNA_SCHEMES.GTR;
    const rates = spec.rates || [1, 1, 1, 1, 1, 1];
    PAIRS.forEach((pr, k) => {
      const v = rates[scheme[k]] == null ? 1 : rates[scheme[k]];
      M[pr[0]][pr[1]] = M[pr[1]][pr[0]] = v;
    });
    return M;
  }

  /* Q = S·Π, normalised so that the expected number of substitutions per unit
     of branch length is one */
  function buildQ(Sx, pi) {
    const n = pi.length;
    const Q = Array.from({ length: n }, () => new Float64Array(n));
    for (let i = 0; i < n; i++) {
      let row = 0;
      for (let j = 0; j < n; j++) if (i !== j) { Q[i][j] = Sx[i][j] * pi[j]; row += Q[i][j]; }
      Q[i][i] = -row;
    }
    let scale = 0;
    for (let i = 0; i < n; i++) scale -= pi[i] * Q[i][i];
    if (scale > 0) for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) Q[i][j] /= scale;
    return Q;
  }

  /* Jacobi eigen-decomposition of a real symmetric matrix */
  function jacobi(Ain) {
    const n = Ain.length;
    const a = Ain.map(r => Float64Array.from(r));
    const v = Array.from({ length: n }, (_, i) => { const r = new Float64Array(n); r[i] = 1; return r; });
    for (let sweep = 0; sweep < 100; sweep++) {
      let off = 0;
      for (let i = 0; i < n; i++) for (let j = i + 1; j < n; j++) off += a[i][j] * a[i][j];
      if (off < 1e-26) break;
      for (let p = 0; p < n - 1; p++) for (let q = p + 1; q < n; q++) {
        if (Math.abs(a[p][q]) < 1e-30) continue;
        const theta = (a[q][q] - a[p][p]) / (2 * a[p][q]);
        const t = Math.sign(theta || 1) / (Math.abs(theta) + Math.sqrt(theta * theta + 1));
        const c = 1 / Math.sqrt(t * t + 1), s = t * c;
        for (let k = 0; k < n; k++) {
          const akp = a[k][p], akq = a[k][q];
          a[k][p] = c * akp - s * akq; a[k][q] = s * akp + c * akq;
        }
        for (let k = 0; k < n; k++) {
          const apk = a[p][k], aqk = a[q][k];
          a[p][k] = c * apk - s * aqk; a[q][k] = s * apk + c * aqk;
        }
        for (let k = 0; k < n; k++) {
          const vkp = v[k][p], vkq = v[k][q];
          v[k][p] = c * vkp - s * vkq; v[k][q] = s * vkp + c * vkq;
        }
      }
    }
    return { values: Array.from({ length: n }, (_, i) => a[i][i]), vectors: v };
  }

  /* discrete gamma: the mean of each of k equal-probability categories
     (Yang 1994), which is what PAML, RAxML and IQ-TREE use */
  function discreteGamma(alpha, k) {
    if (!alpha || !isFinite(alpha) || k <= 1) return { rates: [1], weights: [1] };
    const rates = [];
    for (let i = 0; i < k; i++) {
      const lo = i / k, hi = (i + 1) / k;
      const a = i === 0 ? 0 : qgamma(lo, alpha, alpha);
      const b = i === k - 1 ? Infinity : qgamma(hi, alpha, alpha);
      /* mean of the truncated gamma: [P(a+1, x) difference] / (1/k) */
      const m = (pgammaP(alpha + 1, b === Infinity ? Infinity : b * alpha) - pgammaP(alpha + 1, a * alpha)) * k;
      rates.push(m);
    }
    const mean = rates.reduce((x, y) => x + y, 0) / k;
    return { rates: rates.map(r => r / mean), weights: new Array(k).fill(1 / k) };
  }
  function pgammaP(a, x) {
    if (x === Infinity) return 1;
    if (x <= 0) return 0;
    if (x < a + 1) {
      let sum = 1 / a, del = sum, ap = a;
      for (let i = 0; i < 600; i++) { ap++; del *= x / ap; sum += del; if (Math.abs(del) < Math.abs(sum) * 1e-15) break; }
      return sum * Math.exp(-x + a * Math.log(x) - lgamma(a));
    }
    let b = x + 1 - a, c = 1e300, d = 1 / b, h = d;
    for (let i = 1; i < 600; i++) {
      const an = -i * (i - a); b += 2;
      d = an * d + b; if (Math.abs(d) < 1e-300) d = 1e-300;
      c = b + an / c; if (Math.abs(c) < 1e-300) c = 1e-300;
      d = 1 / d; const del = d * c; h *= del;
      if (Math.abs(del - 1) < 1e-15) break;
    }
    return 1 - Math.exp(-x + a * Math.log(x) - lgamma(a)) * h;
  }
  function qgamma(p, shape, rate) {
    let lo = 1e-10, hi = 1;
    while (pgammaP(shape, hi * rate) < p) hi *= 2;
    for (let i = 0; i < 200; i++) { const m = (lo + hi) / 2; if (pgammaP(shape, m * rate) < p) lo = m; else hi = m; }
    return (lo + hi) / 2;
  }
  function lgamma(x) {
    const c = [76.18009172947146, -86.50532032941677, 24.01409824083091, -1.231739572450155, 0.1208650973866179e-2, -0.5395239384953e-5];
    let y = x, tmp = x + 5.5; tmp -= (x + 0.5) * Math.log(tmp);
    let ser = 1.000000000190015;
    for (let j = 0; j < 6; j++) ser += c[j] / ++y;
    return -tmp + Math.log(2.5066282746310005 * ser / x);
  }

  /* a complete model, ready for the engine */
  function model(spec, data) {
    const S = data ? data.nStates : (spec.nStates || 4);
    let pi;
    if (spec.type === 'aa' && spec.model && g.AAMODELS && g.AAMODELS[spec.model] && !spec.empiricalF) {
      pi = Float64Array.from(g.AAMODELS[spec.model].bf);
    } else if (spec.freqs) pi = Float64Array.from(spec.freqs);
    else if (spec.type === 'dna' && DNA_EQUAL_FREQ[spec.model]) pi = Float64Array.from({ length: S }, () => 1 / S);
    else if (spec.type === 'morph') pi = Float64Array.from({ length: S }, () => 1 / S);
    else pi = data ? Float64Array.from(data.freqs) : Float64Array.from({ length: S }, () => 1 / S);
    let ps = 0; for (let i = 0; i < S; i++) ps += pi[i];
    for (let i = 0; i < S; i++) pi[i] /= ps;

    const Sx = exchangeability(spec, S);
    const Q = buildQ(Sx, pi);
    /* symmetrise: B = Π^{1/2} Q Π^{-1/2} is symmetric for a reversible Q */
    const sq = Float64Array.from(pi, v => Math.sqrt(v));
    const B = Array.from({ length: S }, (_, i) => Float64Array.from({ length: S }, (_, j) => Q[i][j] * sq[i] / sq[j]));
    for (let i = 0; i < S; i++) for (let j = i + 1; j < S; j++) { const m = (B[i][j] + B[j][i]) / 2; B[i][j] = B[j][i] = m; }
    const e = jacobi(B);
    const U = Array.from({ length: S }, (_, i) => Float64Array.from({ length: S }, (_, j) => e.vectors[i][j] / sq[i]));
    const Ui = Array.from({ length: S }, (_, i) => Float64Array.from({ length: S }, (_, j) => e.vectors[j][i] * sq[j]));

    const ncat = spec.ncat || (spec.alpha ? 4 : 1);
    const g = spec.alpha ? discreteGamma(spec.alpha, ncat) : { rates: [1], weights: [1] };
    const pInv = spec.pInv || 0;
    /* with invariable sites the variable ones must run faster, so that the mean
       rate over all sites is still one */
    const catRates = g.rates.map(r => r / (1 - pInv));
    return {
      nStates: S, pi, Q, U, Ui, lambda: Float64Array.from(e.values),
      catRates, catWeights: g.weights, ncat: catRates.length, pInv,
      spec: Object.assign({}, spec),
    };
  }

  /* P(t) for one rate category, written into `out` at offset `off` */
  function transition(M, t, out, off) {
    const S = M.nStates, U = M.U, Ui = M.Ui, lam = M.lambda;
    const ex = new Float64Array(S);
    for (let k = 0; k < S; k++) ex[k] = Math.exp(lam[k] * t);
    for (let i = 0; i < S; i++) {
      const Ui_ = U[i];
      for (let j = 0; j < S; j++) {
        let s = 0;
        for (let k = 0; k < S; k++) s += Ui_[k] * ex[k] * Ui[k][j];
        out[off + i * S + j] = s < 0 ? 0 : s;
      }
    }
  }

  /* ================================================================
     3 · the engine
     ================================================================ */
  /* Buffers for an engine, sized once and reused across candidate topologies.
     A tree of n taxa always has the same number of nodes, so a topology search
     can allocate these once instead of a few megabytes per candidate — which is
     the difference between a search that runs and one that stalls the page. */
  function makePool(nNodes, A, M) {
    const S = A.nStates, nPat = A.nPat, C = M.ncat;
    const sz = C * nPat * S;
    const partial = [], scale = [], P = [];
    for (let k = 0; k < nNodes; k++) {
      partial.push([new Float64Array(sz), new Float64Array(sz)]);
      scale.push([new Float64Array(nPat), new Float64Array(nPat)]);
      P.push([new Float64Array(C * S * S), new Float64Array(C * S * S)]);
    }
    return { partial, scale, P, nNodes, nPat, nStates: S, ncat: C };
  }

  /* F: flattened tree (Tree.flatten) · A: compressed data · M: model
     pool: optional buffers from makePool(), reused between candidates */
  function engine(F, A, M, pool) {
    const S = A.nStates, nPat = A.nPat, C = M.ncat;
    const sz = C * nPat * S;
    const reuse = pool && pool.nNodes >= F.n && pool.nPat === nPat && pool.nStates === S && pool.ncat === C;
    const partial = [], scale = [], P = [];
    const lslot = new Uint8Array(F.n), pslot = new Uint8Array(F.n);
    for (let k = 0; k < F.n; k++) {
      if (reuse) {
        partial.push(F.isTip[k] ? null : pool.partial[k]);
        scale.push(F.isTip[k] ? null : pool.scale[k]);
        P.push(pool.P[k]);
        continue;
      }
      partial.push(F.isTip[k] ? null : [new Float64Array(sz), new Float64Array(sz)]);
      scale.push(F.isTip[k] ? null : [new Float64Array(nPat), new Float64Array(nPat)]);
      P.push([new Float64Array(C * S * S), new Float64Array(C * S * S)]);
    }
    const changedL = [], changedP = [];
    let model2 = M;

    function setModel(m) { model2 = m; }
    function writeP(k, t, flip) {
      const slot = flip ? 1 - pslot[k] : pslot[k];
      const out = P[k][slot];
      for (let c = 0; c < C; c++) transition(model2, Math.max(0, t) * model2.catRates[c], out, c * S * S);
      if (flip) { pslot[k] = slot; changedP.push(k); }
    }
    /* partial likelihoods of node k from its children */
    function computeNode(k, flip) {
      const slot = flip ? 1 - lslot[k] : lslot[k];
      const out = partial[k][slot], sc = scale[k][slot];
      const kids = F.kids[k];
      out.fill(0);
      sc.fill(0);
      for (let c = 0; c < C; c++) {
        const base = c * nPat * S;
        for (let p = 0; p < nPat; p++) {
          const o = base + p * S;
          for (let s = 0; s < S; s++) out[o + s] = 1;
        }
        for (let ki = 0; ki < kids.length; ki++) {
          const kid = kids[ki];
          const pm = P[kid][pslot[kid]];
          const off = c * S * S;
          const childPartial = F.isTip[kid] ? A.tips[F.tipRow[kid]] : partial[kid][lslot[kid]];
          const childBase = F.isTip[kid] ? 0 : c * nPat * S;
          if (S === 4) {
            /* DNA is the case that matters for speed — every likelihood in
               Blocks 3, 5 and 6 goes through here — so the four-state kernel is
               written out: no inner loops, no index arithmetic and no test for
               a zero partial. It is the same arithmetic as the general branch
               below and gives the same numbers; it is simply several times
               faster, which is the difference between a Bayesian run of minutes
               and one of an hour. */
            const q0 = pm[off], q1 = pm[off + 1], q2 = pm[off + 2], q3 = pm[off + 3];
            const q4 = pm[off + 4], q5 = pm[off + 5], q6 = pm[off + 6], q7 = pm[off + 7];
            const q8 = pm[off + 8], q9 = pm[off + 9], qa = pm[off + 10], qb = pm[off + 11];
            const qc = pm[off + 12], qd = pm[off + 13], qe = pm[off + 14], qf = pm[off + 15];
            for (let p = 0; p < nPat; p++) {
              const o = base + (p << 2), co = childBase + (p << 2);
              const c0 = childPartial[co], c1 = childPartial[co + 1];
              const c2 = childPartial[co + 2], c3 = childPartial[co + 3];
              out[o] *= q0 * c0 + q1 * c1 + q2 * c2 + q3 * c3;
              out[o + 1] *= q4 * c0 + q5 * c1 + q6 * c2 + q7 * c3;
              out[o + 2] *= q8 * c0 + q9 * c1 + qa * c2 + qb * c3;
              out[o + 3] *= qc * c0 + qd * c1 + qe * c2 + qf * c3;
            }
          } else {
            for (let p = 0; p < nPat; p++) {
              const o = base + p * S, co = childBase + p * S;
              for (let i = 0; i < S; i++) {
                let s = 0;
                const ro = off + i * S;
                for (let j = 0; j < S; j++) {
                  const cv = childPartial[co + j];
                  if (cv !== 0) s += pm[ro + j] * cv;
                }
                out[o + i] *= s;
              }
            }
          }
          if (!F.isTip[kid]) {
            const ks = scale[kid][lslot[kid]];
            for (let p = 0; p < nPat; p++) sc[p] += ks[p];
          }
        }
        /* rescale a pattern whose partials got too small */
        for (let p = 0; p < nPat; p++) {
          const o = base + p * S;
          let mx = 0;
          for (let s = 0; s < S; s++) if (out[o + s] > mx) mx = out[o + s];
          if (mx > 0 && mx < 1e-100) {
            const f = 1 / mx;
            for (let s = 0; s < S; s++) out[o + s] *= f;
            sc[p] += Math.log(mx) / C;    /* the scale is shared across categories */
          }
        }
      }
      if (flip) { lslot[k] = slot; changedL.push(k); }
    }

    function rootLnL() {
      const pi = model2.pi, w = model2.catWeights, root = 0;
      const pr = partial[root][lslot[root]], sc = scale[root][lslot[root]];
      let lnL = 0;
      for (let p = 0; p < nPat; p++) {
        let site = 0;
        for (let c = 0; c < C; c++) {
          const o = c * nPat * S + p * S;
          let s = 0;
          for (let i = 0; i < S; i++) s += pi[i] * pr[o + i];
          site += w[c] * s;
        }
        site *= (1 - model2.pInv);
        if (model2.pInv > 0 && A.constMask[p]) {
          let inv = 0;
          for (let i = 0; i < S; i++) if (A.constMask[p] & (1 << i)) inv += pi[i];
          site += model2.pInv * inv * Math.exp(-sc[p]);
          /* the invariable part is not scaled, so it is brought to the same scale */
        }
        lnL += A.weights[p] * (Math.log(site) + sc[p]);
      }
      return lnL;
    }

    function full(lens) {
      for (let k = 1; k < F.n; k++) writeP(k, lens[k], false);
      for (let i = 0; i < F.post.length; i++) {
        const k = F.post[i];
        if (!F.isTip[k]) computeNode(k, false);
      }
      return rootLnL();
    }

    /* recompute only what changes when the given branches change */
    function propose(edges, lens) {
      changedL.length = 0; changedP.length = 0;
      const dirty = new Set();
      edges.forEach(k => {
        writeP(k, lens[k], true);
        let node2 = F.parent[k];
        while (node2 >= 0) { dirty.add(node2); node2 = F.parent[node2]; }
      });
      for (let i = 0; i < F.post.length; i++) {
        const k = F.post[i];
        if (dirty.has(k)) computeNode(k, true);
      }
      return rootLnL();
    }
    function accept() { changedL.length = 0; changedP.length = 0; }
    function reject() {
      changedL.forEach(k => { lslot[k] = 1 - lslot[k]; });
      changedP.forEach(k => { pslot[k] = 1 - pslot[k]; });
      changedL.length = 0; changedP.length = 0;
    }
    /* per-site log-likelihoods, for the topology tests of Block 5 */
    function siteLnL() {
      const pi = model2.pi, w = model2.catWeights;
      const pr = partial[0][lslot[0]], sc = scale[0][lslot[0]];
      const out = new Float64Array(nPat);
      for (let p = 0; p < nPat; p++) {
        let site = 0;
        for (let c = 0; c < C; c++) {
          const o = c * nPat * S + p * S;
          let s = 0;
          for (let i = 0; i < S; i++) s += pi[i] * pr[o + i];
          site += w[c] * s;
        }
        site *= (1 - model2.pInv);
        if (model2.pInv > 0 && A.constMask[p]) {
          let inv = 0;
          for (let i = 0; i < S; i++) if (A.constMask[p] & (1 << i)) inv += pi[i];
          site += model2.pInv * inv * Math.exp(-sc[p]);
        }
        out[p] = Math.log(site) + sc[p];
      }
      return out;
    }
    return { full, propose, accept, reject, setModel, rootLnL, siteLnL, nPat, nStates: S };
  }

  /* ================================================================
     4 · optimisation
     ================================================================ */
  /* Brent's method without derivatives, minimising f */
  function brentMin(f, a, b, tol, maxIt) {
    const gold = 0.3819660;
    let x = a + gold * (b - a), w = x, v = x;
    let fx = f(x), fw = fx, fv = fx;
    let d = 0, e = 0;
    tol = tol || 1e-6;
    for (let it = 0; it < (maxIt || 60); it++) {
      const m = 0.5 * (a + b);
      const tol1 = tol * Math.abs(x) + 1e-10, tol2 = 2 * tol1;
      if (Math.abs(x - m) <= tol2 - 0.5 * (b - a)) break;
      let u;
      if (Math.abs(e) > tol1) {
        const r = (x - w) * (fx - fv), q0 = (x - v) * (fx - fw);
        let p = (x - v) * q0 - (x - w) * r, q = 2 * (q0 - r);
        if (q > 0) p = -p; else q = -q;
        const etemp = e; e = d;
        if (Math.abs(p) >= Math.abs(0.5 * q * etemp) || p <= q * (a - x) || p >= q * (b - x)) {
          e = (x >= m ? a - x : b - x); d = gold * e;
        } else { d = p / q; u = x + d; if (u - a < tol2 || b - u < tol2) d = (m > x ? tol1 : -tol1); }
      } else { e = (x >= m ? a - x : b - x); d = gold * e; }
      u = Math.abs(d) >= tol1 ? x + d : x + (d > 0 ? tol1 : -tol1);
      const fu = f(u);
      if (fu <= fx) {
        if (u >= x) a = x; else b = x;
        v = w; w = x; x = u; fv = fw; fw = fx; fx = fu;
      } else {
        if (u < x) a = u; else b = u;
        if (fu <= fw || w === x) { v = w; w = u; fv = fw; fw = fu; }
        else if (fu <= fv || v === x || v === w) { v = u; fv = fu; }
      }
    }
    return { x, f: fx };
  }

  /* how many free parameters a model has, for AIC and friends */
  function nParams(spec, nStates) {
    let k = 0;
    if (spec.type === 'dna') {
      const scheme = DNA_SCHEMES[spec.model] || DNA_SCHEMES.GTR;
      k += Math.max(0, Math.max(...scheme) );                 // rates, one fixed to 1
      if (!DNA_EQUAL_FREQ[spec.model]) k += 3;                // base frequencies
    } else if (spec.type === 'aa') {
      if (spec.empiricalF) k += (nStates || 20) - 1;
    } else if (spec.type === 'morph') {
      /* Mk: no free rate parameters; Mkv adds none either */
    }
    if (spec.alpha) k += 1;
    if (spec.pInv) k += 1;
    return k;
  }

  /* fit branch lengths (and the model's free parameters) on a fixed topology */
  /* A model specification carries two arrays, the exchangeability rates and the
     equilibrium frequencies, and the optimiser writes into them. A shallow copy
     would leave the caller's arrays shared, so a second analysis would silently
     start from the first one's fitted values. */
  function cloneSpec(spec) {
    const s = Object.assign({}, spec);
    if (Array.isArray(s.rates) || ArrayBuffer.isView(s.rates)) s.rates = Array.from(s.rates);
    if (Array.isArray(s.freqs) || ArrayBuffer.isView(s.freqs)) s.freqs = Array.from(s.freqs);
    if (s.catRates) s.catRates = Array.from(s.catRates);
    return s;
  }

  function fit(tree, A, spec, opts) {
    opts = opts || {};
    const F = Tree.flatten(tree);
    spec = Object.assign({ type: 'dna', model: 'GTR', ncat: 4 }, cloneSpec(spec));
    if (spec.type === 'dna' && !spec.rates) spec.rates = [1, 1, 1, 1, 1, 1];
    let M = model(spec, A);
    const lik = engine(F, A, M);
    const lens = Float64Array.from(F.len, v => Math.max(1e-6, v || 0.02));
    /* the two branches at the root of a reversible model act as one */
    const rootKids = F.kids[0];
    const pairAtRoot = rootKids.length === 2;
    let lnL = lik.full(lens);

    /* A short sweep over plausible starting values for the model.

       Coordinate descent alternates branch lengths and model parameters, and
       from a bad starting point the two can pull each other into a poor joint
       optimum and stay there. On the Bursera example, starting from α = 0.5 and
       κ = 2 the fit settles at −7112.26 on a topology that reaches −7081.26
       when it starts from α = 0.16 and κ = 10 — thirty-one log-likelihood units
       lost to where the search happened to begin. The Bayesian sampler of Block
       6, which wanders instead of descending, is what found the better basin
       and exposed this. A dozen extra full likelihoods here remove the trap. */
    if (!opts.fixModel && opts.gridStart !== false) {
      const sweep = (get, set, values) => {
        let bestV = get(), bestL = lnL;
        values.forEach(v => {
          set(v);
          M = model(spec, A); lik.setModel(M);
          const l = lik.full(lens);
          if (l > bestL) { bestL = l; bestV = v; }
        });
        set(bestV);
        M = model(spec, A); lik.setModel(M);
        lnL = lik.full(lens);
      };
      if (spec.alpha) sweep(() => spec.alpha, v => { spec.alpha = v; }, [0.05, 0.1, 0.2, 0.5, 1, 2, 5]);
      if (spec.type === 'dna') {
        const sc0 = DNA_SCHEMES[spec.model] || DNA_SCHEMES.GTR;
        const nf0 = Math.max.apply(null, sc0);
        for (let r = 1; r <= nf0; r++) {
          sweep(() => spec.rates[r], v => { spec.rates[r] = v; }, [0.5, 1, 2, 5, 10, 20]);
        }
      }
      if (spec.alpha) sweep(() => spec.alpha, v => { spec.alpha = v; }, [0.05, 0.1, 0.2, 0.5, 1, 2, 5]);
    }

    const maxLen = opts.maxLen || 5;
    /* Coordinate descent converges slowly near the optimum, so a loose stopping
       rule leaves a tenth of a log-likelihood unit on the table: with 10 passes
       and tol 1e-4 this fit reached −7165.4946 where phangorn reaches
       −7165.3960, and with 30 passes and tol 1e-7 it reaches −7165.3959 — a
       shade better. Twenty passes and 1e-6 is where the curve flattens. */
    const passes = opts.passes || 20;
    const tolBranch = opts.tolBranch || 1e-6;
    let it = 0;
    for (let pass = 0; pass < passes; pass++) {
      const before = lnL;
      /* --- branch lengths --- */
      for (let k = 1; k < F.n; k++) {
        if (pairAtRoot && k === rootKids[1]) continue;
        if (pairAtRoot && k === rootKids[0]) {
          const other = rootKids[1];
          const res = brentMin(s => {
            lens[k] = s / 2; lens[other] = s / 2;
            const v = -lik.propose([k, other], lens); lik.accept(); return v;
          }, 1e-9, maxLen * 2, tolBranch, 30);
          lens[k] = res.x / 2; lens[other] = res.x / 2;
        } else {
          const res = brentMin(t => {
            lens[k] = t;
            const v = -lik.propose([k], lens); lik.accept(); return v;
          }, 1e-9, maxLen, tolBranch, 30);
          lens[k] = res.x;
        }
        it++;
      }
      lnL = lik.full(lens);
      /* --- model parameters --- */
      const refit = (lo, hi, get, set) => {
        const res = brentMin(x => {
          set(Math.exp(x));
          M = model(spec, A); lik.setModel(M);
          return -lik.full(lens);
        }, Math.log(lo), Math.log(hi), 1e-4, 25);
        set(Math.exp(res.x));
        M = model(spec, A); lik.setModel(M);
        lnL = lik.full(lens);
      };
      /* opts.fixModel keeps the substitution model exactly as it was given and
         optimises only the branch lengths. That is what a bootstrap replicate
         wants — RAxML and IQ-TREE both hold the model at the values estimated
         from the real data — and it makes a replicate several times cheaper. */
      if (!opts.fixModel && spec.type === 'dna') {
        const scheme = DNA_SCHEMES[spec.model] || DNA_SCHEMES.GTR;
        const nFree = Math.max(...scheme);
        for (let r = 1; r <= nFree; r++) refit(0.01, 300, () => spec.rates[r], v => { spec.rates[r] = v; });
      }
      /* The equilibrium frequencies, when the caller asks for them to be
         ESTIMATED rather than counted.

         By default they are the empirical ones — the proportions of A, C, G and
         T in the alignment — which is what `phangorn::optim.pml` does unless
         optBf is turned on, and what Block 3's model selection assumes. But the
         empirical proportions are not the maximum-likelihood ones: on the
         Bursera example the fitted frequencies come out (0.271, 0.231, 0.249,
         0.249) against the observed (0.266, 0.202, 0.274, 0.258), and the
         difference is worth some thirty log-likelihood units. The Bayesian
         sampler of Block 6, which estimates them as a matter of course, is what
         made the gap visible. */
      if (!opts.fixModel && opts.optFreqs && spec.freqs && spec.freqs.length > 1
          && spec.type !== 'morph' && !DNA_EQUAL_FREQ[spec.model]) {
        const nF = spec.freqs.length;
        for (let i = 0; i < nF; i++) {
          const cur = Array.from(spec.freqs);
          const setF = w => {
            const f = cur.slice();
            f[i] = w;
            let s = 0; f.forEach(z => { s += z; });
            spec.freqs = f.map(z => Math.max(1e-6, z / s));
            M = model(spec, A); lik.setModel(M);
          };
          const res = brentMin(x => { setF(Math.exp(x)); return -lik.full(lens); },
            Math.log(1e-3), Math.log(5), 1e-5, 20);
          setF(Math.exp(res.x));
        }
        lnL = lik.full(lens);
      }
      if (!opts.fixModel && spec.alpha) refit(0.02, 100, () => spec.alpha, v => { spec.alpha = v; });
      if (!opts.fixModel && spec.pInv != null && spec.pInv > 0) {
        /* pInv lives in (0, 1) and cannot exceed the fraction of constant sites */
        const maxInv = Math.min(0.99, A.nConstant / A.nSites - 1e-4);
        if (maxInv > 0.01) {
          const res = brentMin(x => {
            spec.pInv = 1 / (1 + Math.exp(-x));
            M = model(spec, A); lik.setModel(M);
            return -lik.full(lens);
          }, -6, Math.log(maxInv / (1 - maxInv)), 1e-4, 25);
          spec.pInv = 1 / (1 + Math.exp(-res.x));
          M = model(spec, A); lik.setModel(M);
          lnL = lik.full(lens);
        }
      }
      /* --- one global scale of the tree ---
         The shape α, the proportion of invariable sites and the branch lengths
         are strongly correlated: optimising each branch on its own crawls
         towards the optimum one pass at a time. Rescaling the whole tree by a
         single factor costs one parameter and breaks that correlation, which is
         the same idea as the up–down operator of a Bayesian sampler. Without
         this the fit stopped 0.2 log-likelihood units short of phangorn. */
      if (spec.alpha || spec.pInv) {
        const base2 = Float64Array.from(lens);
        const rs = brentMin(x => {
          const s = Math.exp(x);
          for (let k = 1; k < F.n; k++) lens[k] = base2[k] * s;
          return -lik.full(lens);
        }, Math.log(0.5), Math.log(2), 1e-6, 25);
        const s = Math.exp(rs.x);
        for (let k = 1; k < F.n; k++) lens[k] = base2[k] * s;
        lnL = lik.full(lens);
      }
      if (Math.abs(lnL - before) < (opts.tol || 1e-6)) break;
    }

    const k = nParams(spec, A.nStates) + (F.n - 1) - (pairAtRoot ? 1 : 0);
    const nS = A.nSites;
    return {
      lnL, spec: Object.assign({}, spec), k, nSites: nS, nPatterns: A.nPat,
      AIC: -2 * lnL + 2 * k,
      AICc: -2 * lnL + 2 * k + (nS - k - 1 > 0 ? 2 * k * (k + 1) / (nS - k - 1) : Infinity),
      BIC: -2 * lnL + k * Math.log(nS),
      tree: Tree.unflatten(F, lens), lens: Array.from(lens), model: M, iterations: it,
    };
  }
  /* ================================================================
     5 · several partitions, one tree
     ================================================================ */
  /* A concatenated data set is not one alignment: rbcL and ITS do not evolve at
     the same speed, nor with the same base composition, and forcing one model on
     both is the commonest way to get a confident wrong answer. The standard
     answer — the one RAxML calls a proportional model and IQ-TREE calls
     edge-linked — is to share ONE topology and ONE set of relative branch
     lengths across partitions, and give each partition its own substitution
     model plus one rate multiplier that stretches or shrinks the whole tree for
     that partition. The multipliers are normalised so that the average rate,
     weighted by the number of sites, is one: without that constraint the tree
     and the multipliers are not identifiable (any tree twice as long with every
     multiplier halved has the same likelihood).

     parts: [{ A, spec, name }]. Returns the shared tree, the fitted spec of each
     partition, its multiplier and its own log-likelihood. */
  function fitPartitioned(tree, parts, opts) {
    opts = opts || {};
    const F = Tree.flatten(tree);
    const P = parts.map(p => {
      const spec = Object.assign({ type: 'dna', model: 'GTR', ncat: 4 }, cloneSpec(p.spec));
      if (spec.type === 'dna' && !spec.rates) spec.rates = [1, 1, 1, 1, 1, 1];
      if (!spec.freqs) spec.freqs = Array.from(p.A.freqs);
      const M = model(spec, p.A);
      return { A: p.A, name: p.name, spec, M, lik: engine(F, p.A, M), rate: p.rate || 1, lnL: 0 };
    });
    const totalSites = P.reduce((s, p) => s + p.A.nSites, 0);
    const lens = Float64Array.from(F.len, v => Math.max(1e-6, v || 0.02));
    const rootKids = F.kids[0];
    const pairAtRoot = rootKids.length === 2;
    const scaled = new Float64Array(F.n);

    /* every partition sees the same tree, stretched by its own multiplier */
    const put = p => { for (let k = 1; k < F.n; k++) scaled[k] = lens[k] * p.rate; return scaled; };
    const totalLnL = () => {
      let s = 0;
      P.forEach(p => { p.lnL = p.lik.full(put(p)); s += p.lnL; });
      return s;
    };
    /* the multipliers are only defined up to a common factor: fix it by making
       the site-weighted mean rate one, and put the slack into the branches */
    const renormalise = () => {
      let m = 0;
      P.forEach(p => { m += p.rate * p.A.nSites; });
      m /= totalSites;
      if (!(m > 0)) return;
      P.forEach(p => { p.rate /= m; });
      for (let k = 1; k < F.n; k++) lens[k] *= m;
    };

    let lnL = totalLnL();
    const passes = opts.passes || 12;
    const maxLen = opts.maxLen || 5;
    /* A warm start with the models tied.

       Letting every partition find its own κ and α from the very first pass,
       each on a few hundred sites, drops the fit into a poor local optimum: on
       the test data it landed BELOW the single-model fit it strictly contains,
       which is impossible at the true optimum. So the first passes optimise one
       model shared by all partitions — the single-model answer, which is a good
       starting point by construction — and only then is each partition allowed
       its own. */
    /* the shared model, but every partition keeps its own base composition:
       tying the rates and the gamma is a starting point, tying the frequencies
       would be a different (and worse) model */
    const tieModels = () => {
      const src = P[0].spec;
      P.forEach(p => {
        const own = p.spec.freqs;
        p.spec = cloneSpec(src);
        if (own) p.spec.freqs = Array.from(own);
        p.M = model(p.spec, p.A);
        p.lik.setModel(p.M);
      });
    };
    const optimiseTied = () => {
      const s0 = P[0].spec;
      const refit = (lo, hi, set) => {
        const r = brentMin(x => {
          set(Math.exp(x));
          tieModels();
          return -totalLnL();
        }, Math.log(lo), Math.log(hi), 1e-4, 22);
        set(Math.exp(r.x));
        tieModels();
      };
      if (s0.type === 'dna') {
        const scheme = DNA_SCHEMES[s0.model] || DNA_SCHEMES.GTR;
        const nFree = Math.max(...scheme);
        for (let r = 1; r <= nFree; r++) refit(0.01, 300, v => { s0.rates[r] = v; });
      }
      if (s0.alpha) refit(0.02, 100, v => { s0.alpha = v; });
    };

    for (let pass = 0; pass < passes; pass++) {
      const before = lnL;
      const tied = !opts.fixModel && pass < (opts.tiedPasses == null ? 3 : opts.tiedPasses);
      /* --- the shared branch lengths, on the sum over partitions --- */
      for (let k = 1; k < F.n; k++) {
        if (pairAtRoot && k === rootKids[1]) continue;
        const twin = pairAtRoot && k === rootKids[0] ? rootKids[1] : -1;
        const res = brentMin(t => {
          lens[k] = twin >= 0 ? t / 2 : t;
          if (twin >= 0) lens[twin] = t / 2;
          let s = 0;
          P.forEach(p => { const v = put(p); s += p.lik.propose(twin >= 0 ? [k, twin] : [k], v); p.lik.accept(); });
          return -s;
        }, 1e-9, twin >= 0 ? maxLen * 2 : maxLen, opts.tolBranch || 1e-6, 25);
        lens[k] = twin >= 0 ? res.x / 2 : res.x;
        if (twin >= 0) lens[twin] = res.x / 2;
      }
      lnL = totalLnL();
      /* --- each partition's own model and its multiplier --- */
      if (tied) optimiseTied();
      else if (!opts.fixModel) {
        P.forEach(p => {
          const refit = (lo, hi, set) => {
            const r = brentMin(x => {
              set(Math.exp(x));
              p.M = model(p.spec, p.A); p.lik.setModel(p.M);
              return -p.lik.full(put(p));
            }, Math.log(lo), Math.log(hi), 1e-4, 22);
            set(Math.exp(r.x));
            p.M = model(p.spec, p.A); p.lik.setModel(p.M);
          };
          if (p.spec.type === 'dna') {
            const scheme = DNA_SCHEMES[p.spec.model] || DNA_SCHEMES.GTR;
            const nFree = Math.max(...scheme);
            for (let r = 1; r <= nFree; r++) refit(0.01, 300, v => { p.spec.rates[r] = v; });
          }
          if (p.spec.alpha) refit(0.02, 100, v => { p.spec.alpha = v; });
          if (p.spec.pInv != null && p.spec.pInv > 0) {
            const maxInv = Math.min(0.99, p.A.nConstant / p.A.nSites - 1e-4);
            if (maxInv > 0.01) {
              const r = brentMin(x => {
                p.spec.pInv = 1 / (1 + Math.exp(-x));
                p.M = model(p.spec, p.A); p.lik.setModel(p.M);
                return -p.lik.full(put(p));
              }, -6, Math.log(maxInv / (1 - maxInv)), 1e-4, 22);
              p.spec.pInv = 1 / (1 + Math.exp(-r.x));
              p.M = model(p.spec, p.A); p.lik.setModel(p.M);
            }
          }
        });
      }
      /* --- the rate multipliers --- */
      if (P.length > 1 && !opts.fixRates && !tied) {
        P.forEach(p => {
          const r = brentMin(x => {
            p.rate = Math.exp(x);
            return -p.lik.full(put(p));
          }, Math.log(0.02), Math.log(50), 1e-5, 25);
          p.rate = Math.exp(r.x);
        });
        renormalise();
      }
      /* --- one global scale of the shared tree ---
         The same correction the single-partition fit needs, and for the same
         reason: the shape of the gamma and the branch lengths pull against each
         other, and optimising one branch at a time crawls. Without this the
         partitioned fit landed 0.19 log-likelihood units BELOW the single model
         it strictly contains — an impossible result that was purely a failure
         to converge. */
      {
        const base2 = Float64Array.from(lens);
        const rs = brentMin(x => {
          const s = Math.exp(x);
          for (let k = 1; k < F.n; k++) lens[k] = base2[k] * s;
          return -totalLnL();
        }, Math.log(0.5), Math.log(2), 1e-6, 25);
        const s = Math.exp(rs.x);
        for (let k = 1; k < F.n; k++) lens[k] = base2[k] * s;
      }
      lnL = totalLnL();
      if (Math.abs(lnL - before) < (opts.tol || 1e-5)) break;
    }

    /* the parameter count: the shared branches once, plus every partition's own
       model, plus one multiplier per partition beyond the first */
    let k = (F.n - 1) - (pairAtRoot ? 1 : 0) + Math.max(0, P.length - 1);
    P.forEach(p => { k += nParams(p.spec, p.A.nStates); });
    return {
      lnL, k, nSites: totalSites, tree: Tree.unflatten(F, lens), lens: Array.from(lens),
      partitions: P.map(p => ({
        name: p.name, spec: Object.assign({}, p.spec), rate: p.rate, lnL: p.lnL,
        nSites: p.A.nSites, nPatterns: p.A.nPat, model: p.M,
      })),
      AIC: -2 * lnL + 2 * k,
      AICc: -2 * lnL + 2 * k + (totalSites - k - 1 > 0 ? 2 * k * (k + 1) / (totalSites - k - 1) : Infinity),
      BIC: -2 * lnL + k * Math.log(totalSites),
    };
  }

  /* the same likelihood, for a tree that is only being scored (topology search
     over a partitioned model): branch lengths kept, nothing optimised */
  function lnLPartitioned(tree, parts) {
    const F = Tree.flatten(tree);
    let s = 0;
    const scaled = new Float64Array(F.n);
    parts.forEach(p => {
      const spec = Object.assign({ type: 'dna', model: 'GTR', ncat: 4 }, cloneSpec(p.spec));
      if (!spec.freqs) spec.freqs = Array.from(p.A.freqs);
      const M = p.model || model(spec, p.A);
      const lik = engine(F, p.A, M);
      for (let k = 1; k < F.n; k++) scaled[k] = Math.max(1e-9, (F.len[k] || 0.02) * (p.rate || 1));
      s += lik.full(scaled);
    });
    return s;
  }


  /* Fast fit for model selection: the topology AND the relative branch lengths
     are held fixed, and only the model's own parameters plus one global scale
     for the tree are optimised. That is the strategy jModelTest and ModelFinder
     use to compare dozens of models without refitting every branch each time —
     here it is about twenty times faster, and the ranking it produces is the
     same, because every model sees exactly the same tree. The winner is then
     refitted properly with fit(). */
  function fitFast(tree, A, spec, opts) {
    opts = opts || {};
    const F = Tree.flatten(tree);
    spec = Object.assign({ type: 'dna', model: 'GTR', ncat: 4 }, cloneSpec(spec));
    if (spec.type === 'dna' && !spec.rates) spec.rates = [1, 1, 1, 1, 1, 1];
    let M = model(spec, A);
    const lik = engine(F, A, M);
    const base = Float64Array.from(F.len, v => Math.max(1e-6, v || 0.02));
    const lens = Float64Array.from(base);
    let scale = 1;
    const setScale = s => { scale = s; for (let k = 0; k < F.n; k++) lens[k] = base[k] * s; };
    let lnL = lik.full(lens);
    const evalAt = () => { M = model(spec, A); lik.setModel(M); return lik.full(lens); };
    for (let pass = 0; pass < (opts.passes || 6); pass++) {
      const before = lnL;
      /* global scale of the tree */
      const rs = brentMin(x => { setScale(Math.exp(x)); return -lik.full(lens); }, Math.log(0.05), Math.log(20), 1e-5, 30);
      setScale(Math.exp(rs.x));
      lnL = lik.full(lens);
      /* model parameters */
      const refit = (lo, hi, get, set) => {
        const r = brentMin(x => { set(Math.exp(x)); return -evalAt(); }, Math.log(lo), Math.log(hi), 1e-4, 25);
        set(Math.exp(r.x)); lnL = evalAt();
      };
      if (spec.type === 'dna') {
        const scheme = DNA_SCHEMES[spec.model] || DNA_SCHEMES.GTR;
        const nFree = Math.max(...scheme);
        for (let r = 1; r <= nFree; r++) refit(0.01, 300, () => spec.rates[r], v => { spec.rates[r] = v; });
      }
      if (spec.alpha) refit(0.02, 100, () => spec.alpha, v => { spec.alpha = v; });
      if (spec.pInv != null && spec.pInv > 0) {
        const maxInv = Math.min(0.99, A.nConstant / A.nSites - 1e-4);
        if (maxInv > 0.01) {
          const r = brentMin(x => { spec.pInv = 1 / (1 + Math.exp(-x)); return -evalAt(); },
            -6, Math.log(maxInv / (1 - maxInv)), 1e-4, 25);
          spec.pInv = 1 / (1 + Math.exp(-r.x)); lnL = evalAt();
        }
      }
      if (Math.abs(lnL - before) < (opts.tol || 1e-3)) break;
    }
    const k = nParams(spec, A.nStates) + (F.n - 1) - (F.kids[0].length === 2 ? 1 : 0);
    const nS = A.nSites;
    return {
      lnL, spec: Object.assign({}, spec), k, nSites: nS, nPatterns: A.nPat, scale,
      AIC: -2 * lnL + 2 * k,
      AICc: -2 * lnL + 2 * k + (nS - k - 1 > 0 ? 2 * k * (k + 1) / (nS - k - 1) : Infinity),
      BIC: -2 * lnL + k * Math.log(nS),
      lens: Array.from(lens), model: M, fast: true,
    };
  }

  Object.assign(Like, {
    compress, model, engine, makePool, cloneSpec, transition, buildQ, exchangeability, jacobi,
    fitPartitioned, lnLPartitioned,
    discreteGamma, qgamma, pgammaP, lgamma, brentMin, fit, fitFast, nParams,
    DNA_SCHEMES, DNA_EQUAL_FREQ, PAIRS,
  });
  g.Like = Like;
}
LikeCore(typeof window !== 'undefined' ? window : self);
