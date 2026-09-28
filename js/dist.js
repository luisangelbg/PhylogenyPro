/* PhylogenyPro — molecular phylogenetics in the browser.
   Copyright (C) 2026  Luis Ángel Barrera-Guzmán

   This program is free software: you can redistribute it and/or modify it under
   the terms of the GNU General Public License as published by the Free Software
   Foundation, either version 3 of the License, or (at your option) any later
   version. It is distributed in the hope that it will be useful, but WITHOUT ANY
   WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
   A PARTICULAR PURPOSE. See the GNU General Public License, in the file LICENSE
   at the root of this program, or <https://www.gnu.org/licenses/>. */

/* PhylogenyPro — distance matrices corrected by a model.

   Every function returns the number of substitutions per site, or null when the
   correction has no solution: two sequences can be so different that no amount
   of arithmetic recovers how many changes happened, and saying "null" is more
   honest than returning a large number that looks like data. The caller decides
   what to do with it (the matrix builder caps it and warns).

   Missing data are handled by pairwise deletion by default — each pair is
   compared over the sites both of them have — with complete deletion available,
   which is stricter and loses more. */

/* Wrapped in a named function so that a Web Worker can be built from its own
   source text: see js/pool.js. Nothing here touches the DOM. */
function DistCore(g) {
const Dist = {};

  const DNA_OK = { A: 1, C: 1, G: 1, T: 1, U: 1 };
  const purine = { A: 1, G: 1 }, pyrimidine = { C: 1, T: 1, U: 1 };

  /* ---------------- pairwise counts ---------------- */
  function countPair(a, b, type) {
    let n = 0, diff = 0, P = 0, Q = 0, P1 = 0, P2 = 0;
    const L = Math.min(a.length, b.length);
    for (let i = 0; i < L; i++) {
      let x = a[i], y = b[i];
      if (type === 'dna') {
        if (x === 'U') x = 'T';
        if (y === 'U') y = 'T';
        if (!DNA_OK[x] || !DNA_OK[y]) continue;
      } else {
        if (x === '-' || y === '-' || x === '?' || y === '?' || x === 'X' || y === 'X' || x === 'N' || y === 'N') continue;
      }
      n++;
      if (x === y) continue;
      diff++;
      if (type === 'dna') {
        const ts = (purine[x] && purine[y]) || (pyrimidine[x] && pyrimidine[y]);
        if (ts) {
          P++;
          if (purine[x]) P1++; else P2++;       // A<->G and C<->T counted apart, for TN93
        } else Q++;
      }
    }
    return { n, diff, p: n ? diff / n : 0, P: n ? P / n : 0, Q: n ? Q / n : 0, P1: n ? P1 / n : 0, P2: n ? P2 / n : 0 };
  }

  /* base frequencies over a set of sequences */
  function baseFreqs(seqs) {
    const c = { A: 0, C: 0, G: 0, T: 0 };
    let tot = 0;
    seqs.forEach(s => {
      for (let i = 0; i < s.length; i++) {
        let ch = s[i];
        if (ch === 'U') ch = 'T';
        if (DNA_OK[ch]) { c[ch]++; tot++; }
      }
    });
    if (!tot) return { A: 0.25, C: 0.25, G: 0.25, T: 0.25 };
    return { A: c.A / tot, C: c.C / tot, G: c.G / tot, T: c.T / tot };
  }

  /* ---------------- the corrections ---------------- */
  const raw = c => c.p;
  function jc69(c, alpha) {
    const p = c.p;
    if (alpha) return p >= 0.75 ? null : 0.75 * alpha * (Math.pow(1 - 4 * p / 3, -1 / alpha) - 1);
    return p >= 0.75 ? null : -0.75 * Math.log(1 - 4 * p / 3);
  }
  function k80(c, alpha) {
    const P = c.P, Q = c.Q;
    const a = 1 - 2 * P - Q, b = 1 - 2 * Q;
    if (a <= 0 || b <= 0) return null;
    if (alpha) return alpha / 2 * (Math.pow(a, -1 / alpha) + 0.5 * Math.pow(b, -1 / alpha) - 1.5);
    return -0.5 * Math.log(a) - 0.25 * Math.log(b);
  }
  function f81(c, freqs, alpha) {
    const B = 1 - (freqs.A * freqs.A + freqs.C * freqs.C + freqs.G * freqs.G + freqs.T * freqs.T);
    const p = c.p;
    if (p >= B) return null;
    if (alpha) return B * alpha * (Math.pow(1 - p / B, -1 / alpha) - 1);
    return -B * Math.log(1 - p / B);
  }
  /* Tamura & Nei (1993): transitions between purines and between pyrimidines
     get their own rate, and the base frequencies are used */
  function tn93(c, freqs, alpha) {
    const gA = freqs.A, gC = freqs.C, gG = freqs.G, gT = freqs.T;
    const gR = gA + gG, gY = gC + gT;
    if (gR <= 0 || gY <= 0 || gA <= 0 || gG <= 0 || gC <= 0 || gT <= 0) return null;
    const P1 = c.P1, P2 = c.P2, Q = c.Q;
    const a1 = 1 - gR * P1 / (2 * gA * gG) - Q / (2 * gR);
    const a2 = 1 - gY * P2 / (2 * gC * gT) - Q / (2 * gY);
    const b = 1 - Q / (2 * gR * gY);
    if (a1 <= 0 || a2 <= 0 || b <= 0) return null;
    const k1 = 2 * gA * gG / gR, k2 = 2 * gC * gT / gY;
    const k3 = 2 * (gR * gY - gA * gG * gY / gR - gC * gT * gR / gY);
    if (alpha) return k1 * alpha * (Math.pow(a1, -1 / alpha) - 1) + k2 * alpha * (Math.pow(a2, -1 / alpha) - 1) + k3 * alpha * (Math.pow(b, -1 / alpha) - 1);
    return -k1 * Math.log(a1) - k2 * Math.log(a2) - k3 * Math.log(b);
  }
  /* LogDet and paralinear: the two distances that stay right when the base
     composition differs between sequences — and two different formulas that the
     literature and the programs both call "LogDet".

       · Lockhart et al. (1994), which is what ape::dist.dna(model="logdet")
         and MEGA compute:      d = −(1/4)·ln det(F) − ln 4
         (the correction term uses equal frequencies, 1/4 each);
       · Lake's (1994) paralinear distance, which uses the frequencies actually
         observed in each sequence:
                                d = −(1/4)·[ln det(F) − ½·ln(∏f₁·∏f₂)].

     They differ by about 2 % on these data. Both are offered, named after their
     authors, because picking one silently would make the app disagree with
     whichever program the user compares it against. */
  function logdet(a, b, variant) {
    const idx = { A: 0, C: 1, G: 2, T: 3, U: 3 };
    const F = [[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]];
    let n = 0;
    const L = Math.min(a.length, b.length);
    for (let i = 0; i < L; i++) {
      const x = idx[a[i]], y = idx[b[i]];
      if (x == null || y == null) continue;
      F[x][y]++; n++;
    }
    if (!n) return null;
    for (let i = 0; i < 4; i++) for (let j = 0; j < 4; j++) F[i][j] /= n;
    const det = det4(F);
    if (!(det > 0)) return null;
    if (variant === 'paralinear') {
      const f1 = [0, 0, 0, 0], f2 = [0, 0, 0, 0];
      for (let i = 0; i < 4; i++) for (let j = 0; j < 4; j++) { f1[i] += F[i][j]; f2[j] += F[i][j]; }
      const prod1 = f1.reduce((s, v) => s * v, 1), prod2 = f2.reduce((s, v) => s * v, 1);
      if (!(prod1 > 0) || !(prod2 > 0)) return null;
      return -0.25 * (Math.log(det) - 0.5 * Math.log(prod1 * prod2));
    }
    return -0.25 * Math.log(det) - Math.log(4);
  }
  function det4(m) {
    /* expansion by minors, fine for a 4×4 */
    const d = (a, b, c2, d2, e, f, g, h, i) => a * (e * i - f * h) - b * (d2 * i - f * g) + c2 * (d2 * h - e * g);
    return m[0][0] * d(m[1][1], m[1][2], m[1][3], m[2][1], m[2][2], m[2][3], m[3][1], m[3][2], m[3][3])
      - m[0][1] * d(m[1][0], m[1][2], m[1][3], m[2][0], m[2][2], m[2][3], m[3][0], m[3][2], m[3][3])
      + m[0][2] * d(m[1][0], m[1][1], m[1][3], m[2][0], m[2][1], m[2][3], m[3][0], m[3][1], m[3][3])
      - m[0][3] * d(m[1][0], m[1][1], m[1][2], m[2][0], m[2][1], m[2][2], m[3][0], m[3][1], m[3][2]);
  }
  /* proteins: the Poisson correction is Jukes–Cantor with twenty states */
  function poisson(c, alpha) {
    const p = c.p;
    if (p >= 0.95) return null;
    if (alpha) return 0.95 * alpha * (Math.pow(1 - 20 * p / 19, -1 / alpha) - 1);
    return -0.95 * Math.log(1 - 20 * p / 19);
  }
  /* Kimura's (1983) empirical correction for proteins, the one MEGA calls
     "Poisson with a correction for multiple hits" */
  function kimuraProtein(c) {
    const p = c.p;
    const v = 1 - p - 0.2 * p * p;
    if (v <= 0) return null;
    return -Math.log(v);
  }

  /* ---------------- the matrix ---------------- */
  /* model: 'p' | 'jc' | 'k80' | 'f81' | 'tn93' | 'logdet' | 'poisson' | 'kimura' | 'hamming'
     opts: {alpha, complete (complete deletion), type} */
  function matrix(seqs, model, opts) {
    opts = opts || {};
    const type = opts.type || 'dna';
    const n = seqs.length;
    let use = seqs;
    if (opts.complete) use = completeDeletion(seqs, type);
    const freqs = type === 'dna' ? baseFreqs(use) : null;
    const D = Array.from({ length: n }, () => new Float64Array(n));
    let undefinedPairs = 0, minSites = Infinity;
    for (let i = 0; i < n; i++) {
      for (let j = i + 1; j < n; j++) {
        const c = countPair(use[i], use[j], type);
        minSites = Math.min(minSites, c.n);
        let d;
        switch (model) {
          case 'p': case 'hamming': d = raw(c); break;
          case 'jc': d = jc69(c, opts.alpha); break;
          case 'k80': d = k80(c, opts.alpha); break;
          case 'f81': d = f81(c, freqs, opts.alpha); break;
          case 'tn93': d = tn93(c, freqs, opts.alpha); break;
          case 'logdet': d = logdet(use[i], use[j]); break;
          case 'paralinear': d = logdet(use[i], use[j], 'paralinear'); break;
          case 'poisson': d = poisson(c, opts.alpha); break;
          case 'kimura': d = kimuraProtein(c); break;
          default: d = raw(c);
        }
        if (d == null || !isFinite(d)) { undefinedPairs++; d = opts.cap == null ? 3 : opts.cap; }
        D[i][j] = D[j][i] = Math.max(0, d);
      }
    }
    return { D, undefinedPairs, minSites: isFinite(minSites) ? minSites : 0, freqs, model, alpha: opts.alpha || null };
  }

  /* complete deletion: only the columns every sequence has */
  function completeDeletion(seqs, type) {
    const L = Math.min(...seqs.map(s => s.length));
    const keep = [];
    for (let c = 0; c < L; c++) {
      let ok = true;
      for (let i = 0; i < seqs.length; i++) {
        const ch = seqs[i][c];
        if (type === 'dna' ? !DNA_OK[ch === 'U' ? 'T' : ch] : (ch === '-' || ch === '?' || ch === 'X' || ch === 'N')) { ok = false; break; }
      }
      if (ok) keep.push(c);
    }
    return seqs.map(s => keep.map(c => s[c]).join(''));
  }

  /* ---------------- saturation ---------------- */
  /* transitions and transversions against the corrected distance: when the
     curves flatten or fall, the marker is saturated */
  function saturation(seqs, opts) {
    opts = opts || {};
    const n = seqs.length, pts = [];
    const freqs = baseFreqs(seqs);
    for (let i = 0; i < n; i++) for (let j = i + 1; j < n; j++) {
      const c = countPair(seqs[i], seqs[j], 'dna');
      const d = tn93(c, freqs) != null ? tn93(c, freqs) : (jc69(c) != null ? jc69(c) : null);
      pts.push({ i, j, d, ts: c.P, tv: c.Q, p: c.p, ratio: c.Q > 0 ? c.P / c.Q : null, n: c.n });
    }
    /* slope of ts and tv against the corrected distance, through the origin */
    const fit = key => {
      let sxy = 0, sxx = 0;
      pts.forEach(pt => { if (pt.d != null) { sxy += pt.d * pt[key]; sxx += pt.d * pt.d; } });
      return sxx > 0 ? sxy / sxx : 0;
    };
    /* Xia's index of substitution saturation (Xia et al. 2003): the observed
       entropy-based Iss against its critical value is not reproduced here; what
       is computed is the simpler and well-defined part — the proportion of the
       observed differences that the correction has to invent */
    const hidden = pts.filter(pt => pt.d != null && pt.d > 0).map(pt => 1 - pt.p / pt.d);
    return {
      points: pts, slopeTs: fit('ts'), slopeTv: fit('tv'),
      meanHidden: hidden.length ? hidden.reduce((a, b) => a + b, 0) / hidden.length : 0,
      maxD: Math.max(...pts.map(p => p.d == null ? 0 : p.d)),
      undefinedPairs: pts.filter(p => p.d == null).length,
    };
  }

  /* ---------------- composition ---------------- */
  /* Chi-square test of homogeneity of base composition across sequences: the
     assumption every standard model makes, and the one most often broken */
  function compositionTest(seqs, type) {
    const alpha = type === 'aa' ? (g.AA_ORDER || 'ARNDCQEGHILKMFPSTWYV').split('') : ['A', 'C', 'G', 'T'];
    const rows = seqs.map(s => {
      const c = {};
      alpha.forEach(a => { c[a] = 0; });
      let tot = 0;
      for (let i = 0; i < s.length; i++) {
        let ch = s[i];
        if (ch === 'U') ch = 'T';
        if (c[ch] != null) { c[ch]++; tot++; }
      }
      return { counts: c, total: tot, freqs: alpha.map(a => tot ? c[a] / tot : 0) };
    });
    const grand = alpha.map((a, k) => {
      let s = 0, t = 0;
      rows.forEach(r => { s += r.counts[a]; t += r.total; });
      return t ? s / t : 0;
    });
    let chi2 = 0, df = 0;
    rows.forEach(r => {
      alpha.forEach((a, k) => {
        const exp = r.total * grand[k];
        if (exp > 0) chi2 += Math.pow(r.counts[a] - exp, 2) / exp;
      });
    });
    df = (rows.length - 1) * (alpha.filter((a, k) => grand[k] > 0).length - 1);
    return { chi2, df, p: df > 0 ? 1 - pchisq(chi2, df) : 1, rows, grand, alphabet: alpha };
  }
  /* regularised lower incomplete gamma, enough for a chi-square tail */
  function pchisq(x, k) {
    if (x <= 0) return 0;
    const a = k / 2, z = x / 2;
    if (z < a + 1) {
      let sum = 1 / a, del = sum, ap = a;
      for (let i = 0; i < 500; i++) { ap++; del *= z / ap; sum += del; if (Math.abs(del) < Math.abs(sum) * 1e-14) break; }
      return sum * Math.exp(-z + a * Math.log(z) - lgamma(a));
    }
    let b = z + 1 - a, c = 1e300, d = 1 / b, h = d;
    for (let i = 1; i < 500; i++) {
      const an = -i * (i - a); b += 2;
      d = an * d + b; if (Math.abs(d) < 1e-300) d = 1e-300;
      c = b + an / c; if (Math.abs(c) < 1e-300) c = 1e-300;
      d = 1 / d; const del = d * c; h *= del;
      if (Math.abs(del - 1) < 1e-14) break;
    }
    return 1 - Math.exp(-z + a * Math.log(z) - lgamma(a)) * h;
  }
  function lgamma(x) {
    const c = [76.18009172947146, -86.50532032941677, 24.01409824083091, -1.231739572450155, 0.1208650973866179e-2, -0.5395239384953e-5];
    let y = x, tmp = x + 5.5; tmp -= (x + 0.5) * Math.log(tmp);
    let ser = 1.000000000190015;
    for (let j = 0; j < 6; j++) ser += c[j] / ++y;
    return -tmp + Math.log(2.5066282746310005 * ser / x);
  }

  Object.assign(Dist, {
    countPair, baseFreqs, matrix, completeDeletion, saturation, compositionTest,
    jc69, k80, f81, tn93, logdet, poisson, kimuraProtein, pchisq, lgamma,
  });
  g.Dist = Dist;
}
DistCore(typeof window !== 'undefined' ? window : self);
