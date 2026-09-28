/* PhylogenyPro — molecular phylogenetics in the browser.
   Copyright (C) 2026  Luis Ángel Barrera-Guzmán

   This program is free software: you can redistribute it and/or modify it under
   the terms of the GNU General Public License as published by the Free Software
   Foundation, either version 3 of the License, or (at your option) any later
   version. It is distributed in the hope that it will be useful, but WITHOUT ANY
   WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
   A PARTICULAR PURPOSE. See the GNU General Public License, in the file LICENSE
   at the root of this program, or <https://www.gnu.org/licenses/>. */

/* PhylogenyPro — the aligner.

   A progressive aligner of the classical kind, the one Clustal made standard and
   MUSCLE refined:

     1. a quick distance between every pair of sequences from their shared
        k-mers (no alignment needed, O(n² L));
     2. a guide tree by UPGMA on those distances;
     3. profile–profile alignment along the guide tree, in postorder, with
        affine gaps (Gotoh 1982) and gap penalties reduced where the profile
        already has gaps;
     4. optional iterative refinement: every branch of the guide tree splits the
        sequences in two, the two sub-alignments are realigned, and the result is
        kept only if the sum-of-pairs score improves (Berger & Munson's idea, the
        third stage of MUSCLE).

   Coding sequences are aligned by codon: translated, aligned as protein, and the
   gaps mapped back onto the nucleotides three at a time — so a gap of one or two
   nucleotides, which is almost always an alignment artefact, cannot happen.

   The aligner is scored in validation/block2 against the true alignment that
   generated the simulated example data. */

const Align = {};

(function () {

  /* ================================================================
     alphabets and scoring
     ================================================================ */
  const DNA = 'ACGT';
  const AA = 'ARNDCQEGHILKMFPSTWYV';

  /* BLOSUM62 (Henikoff & Henikoff 1992, Proc. Natl. Acad. Sci. USA 89: 10915),
     as distributed by the NCBI at ftp.ncbi.nlm.nih.gov/blast/matrices/BLOSUM62
     — a work of the United States government, in the public domain — reordered
     into the order of AA above. See PROCEDENCIA.md. */
  const BLOSUM62_ROWS = [
    /* A */ [4, -1, -2, -2, 0, -1, -1, 0, -2, -1, -1, -1, -1, -2, -1, 1, 0, -3, -2, 0],
    /* R */ [-1, 5, 0, -2, -3, 1, 0, -2, 0, -3, -2, 2, -1, -3, -2, -1, -1, -3, -2, -3],
    /* N */ [-2, 0, 6, 1, -3, 0, 0, 0, 1, -3, -3, 0, -2, -3, -2, 1, 0, -4, -2, -3],
    /* D */ [-2, -2, 1, 6, -3, 0, 2, -1, -1, -3, -4, -1, -3, -3, -1, 0, -1, -4, -3, -3],
    /* C */ [0, -3, -3, -3, 9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1],
    /* Q */ [-1, 1, 0, 0, -3, 5, 2, -2, 0, -3, -2, 1, 0, -3, -1, 0, -1, -2, -1, -2],
    /* E */ [-1, 0, 0, 2, -4, 2, 5, -2, 0, -3, -3, 1, -2, -3, -1, 0, -1, -3, -2, -2],
    /* G */ [0, -2, 0, -1, -3, -2, -2, 6, -2, -4, -4, -2, -3, -3, -2, 0, -2, -2, -3, -3],
    /* H */ [-2, 0, 1, -1, -3, 0, 0, -2, 8, -3, -3, -1, -2, -1, -2, -1, -2, -2, 2, -3],
    /* I */ [-1, -3, -3, -3, -1, -3, -3, -4, -3, 4, 2, -3, 1, 0, -3, -2, -1, -3, -1, 3],
    /* L */ [-1, -2, -3, -4, -1, -2, -3, -4, -3, 2, 4, -2, 2, 0, -3, -2, -1, -2, -1, 1],
    /* K */ [-1, 2, 0, -1, -3, 1, 1, -2, -1, -3, -2, 5, -1, -3, -1, 0, -1, -3, -2, -2],
    /* M */ [-1, -1, -2, -3, -1, 0, -2, -3, -2, 1, 2, -1, 5, 0, -2, -1, -1, -1, -1, 1],
    /* F */ [-2, -3, -3, -3, -2, -3, -3, -3, -1, 0, 0, -3, 0, 6, -4, -2, -2, 1, 3, -1],
    /* P */ [-1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4, 7, -1, -1, -4, -3, -2],
    /* S */ [1, -1, 1, 0, -1, 0, 0, 0, -1, -2, -2, 0, -1, -2, -1, 4, 1, -3, -2, -2],
    /* T */ [0, -1, 0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1, 1, 5, -2, -2, 0],
    /* W */ [-3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1, 1, -4, -3, -2, 11, 2, -3],
    /* Y */ [-2, -2, -2, -3, -2, -1, -2, -3, 2, -1, -1, -2, -1, 3, -3, -2, -2, 2, 7, -1],
    /* V */ [0, -3, -3, -3, -1, -2, -2, -3, -3, 3, 1, -2, 1, -1, -2, -2, 0, -3, -1, 4],
  ];

  /* DNA: a match, a transition and a transversion are three different things.
     These values are in the same scale as BLOSUM62 so that the gap penalties
     can be shared. */
  function dnaMatrix(opts) {
    const match = opts.match == null ? 5 : opts.match;
    const ts = opts.transition == null ? -2 : opts.transition;
    const tv = opts.transversion == null ? -4 : opts.transversion;
    const M = [];
    for (let i = 0; i < 4; i++) {
      M.push([]);
      for (let j = 0; j < 4; j++) M[i].push(i === j ? match : (isTransition(i, j) ? ts : tv));
    }
    return M;
  }

  function alphabetOf(type) { return type === 'aa' ? AA : DNA; }
  function matrixOf(type, opts) { return type === 'aa' ? BLOSUM62_ROWS : dnaMatrix(opts || {}); }

  /* index of a residue, or -1 for a gap, an ambiguity or anything unknown */
  function indexer(alpha) {
    const map = new Int8Array(128).fill(-1);
    for (let i = 0; i < alpha.length; i++) map[alpha.charCodeAt(i)] = i;
    map['U'.charCodeAt(0)] = alpha === DNA ? 3 : map['U'.charCodeAt(0)];   // RNA
    return map;
  }

  /* ================================================================
     1 · k-mer distances
     ================================================================ */
  /* Counts of every k-mer, compared by the fraction shared. Ambiguities and
     gaps break a k-mer, which is the honest thing to do. */
  function kmerCounts(seq, alpha, k, map) {
    const size = Math.pow(alpha.length, k);
    const counts = new Map();
    let code = 0, run = 0;
    for (let i = 0; i < seq.length; i++) {
      const v = map[seq.charCodeAt(i)];
      if (v < 0) { run = 0; code = 0; continue; }
      code = (code * alpha.length + v) % size;
      run++;
      if (run >= k) counts.set(code, (counts.get(code) || 0) + 1);
    }
    return counts;
  }
  function kmerDistance(seqs, type, k) {
    const alpha = alphabetOf(type);
    k = k || (type === 'aa' ? 3 : 5);
    const map = indexer(alpha);
    const counts = seqs.map(s => kmerCounts(s, alpha, k, map));
    const tot = counts.map(c => { let t = 0; c.forEach(v => { t += v; }); return t; });
    const n = seqs.length;
    const D = Array.from({ length: n }, () => new Float64Array(n));
    for (let i = 0; i < n; i++) {
      for (let j = i + 1; j < n; j++) {
        let shared = 0;
        const a = counts[i], b = counts[j];
        const small = a.size < b.size ? a : b, big = small === a ? b : a;
        small.forEach((v, key) => { const w = big.get(key); if (w) shared += Math.min(v, w); });
        const denom = Math.min(tot[i], tot[j]) || 1;
        const d = 1 - shared / denom;
        D[i][j] = D[j][i] = d;
      }
    }
    return D;
  }

  /* ================================================================
     2 · guide tree (UPGMA)
     ================================================================ */
  function upgma(D, labels) {
    const n = D.length;
    let nodes = [];
    for (let i = 0; i < n; i++) nodes.push({ tip: i, label: labels ? labels[i] : String(i), size: 1, height: 0, children: [] });
    const d = D.map(r => Array.from(r));
    const active = nodes.slice();
    const dist = new Map();
    const key = (a, b) => (a.id < b.id ? a.id + '|' + b.id : b.id + '|' + a.id);
    active.forEach((nd, i) => { nd.id = i; });
    let nextId = n;
    const cur = new Map();
    for (let i = 0; i < n; i++) for (let j = i + 1; j < n; j++) cur.set(i + '|' + j, d[i][j]);
    const getD = (a, b) => cur.get(a.id < b.id ? a.id + '|' + b.id : b.id + '|' + a.id);
    while (active.length > 1) {
      let best = Infinity, bi = 0, bj = 1;
      for (let i = 0; i < active.length; i++) for (let j = i + 1; j < active.length; j++) {
        const v = getD(active[i], active[j]);
        if (v < best) { best = v; bi = i; bj = j; }
      }
      const a = active[bi], b = active[bj];
      const node = { id: nextId++, size: a.size + b.size, height: best / 2, children: [a, b], label: null };
      /* average linkage, weighted by cluster size */
      active.forEach(o => {
        if (o === a || o === b) return;
        const v = (getD(a, o) * a.size + getD(b, o) * b.size) / (a.size + b.size);
        cur.set(node.id < o.id ? node.id + '|' + o.id : o.id + '|' + node.id, v);
      });
      active.splice(bj, 1); active.splice(bi, 1);
      active.push(node);
    }
    return active[0];
  }

  /* ================================================================
     3 · profiles and profile–profile alignment
     ================================================================ */
  /* A profile is a set of aligned sequences: its columns hold the frequency of
     each residue plus the fraction of gaps. */
  function makeProfile(rows, seqs, alpha, map) {
    const L = seqs[0].length, S = alpha.length;
    const freq = new Float32Array(L * S);
    const gap = new Float32Array(L);
    const n = seqs.length;
    for (let s = 0; s < n; s++) {
      const seq = seqs[s];
      for (let c = 0; c < L; c++) {
        const v = map[seq.charCodeAt(c)];
        if (v >= 0) freq[c * S + v] += 1 / n;
        else if (seq[c] === '-') gap[c] += 1 / n;
        /* ambiguities count as neither: they contribute nothing to the score */
      }
    }
    return { rows, seqs, L, S, freq, gap, n };
  }

  /* Score between one column of each profile: the expected substitution score. */
  function columnScore(A, B, i, j, M) {
    const S = A.S, oa = i * S, ob = j * S;
    let s = 0;
    for (let a = 0; a < S; a++) {
      const fa = A.freq[oa + a];
      if (fa === 0) continue;
      const row = M[a];
      for (let b = 0; b < S; b++) {
        const fb = B.freq[ob + b];
        if (fb !== 0) s += fa * fb * row[b];
      }
    }
    return s;
  }

  /* Gotoh's affine-gap alignment of two profiles.
     Three matrices (match, gap in A, gap in B) and a byte of traceback per cell. */
  function alignProfiles(A, B, opts) {
    const M = opts.matrix, GO = opts.gapOpen, GE = opts.gapExt;
    const n = A.L, m = B.L;
    const NEG = -1e30;
    const prevM = new Float64Array(m + 1), curM = new Float64Array(m + 1);
    const prevX = new Float64Array(m + 1), curX = new Float64Array(m + 1);   // gap in A (consume B)
    const prevY = new Float64Array(m + 1), curY = new Float64Array(m + 1);   // gap in B (consume A)
    const tb = new Uint8Array((n + 1) * (m + 1));
    /* a gap opened where the profile is already mostly gaps costs less */
    const openA = new Float64Array(n + 1), openB = new Float64Array(m + 1);
    for (let i = 0; i < n; i++) openA[i] = GO * (1 - 0.85 * A.gap[i]);
    for (let j = 0; j < m; j++) openB[j] = GO * (1 - 0.85 * B.gap[j]);
    const extA = new Float64Array(n + 1), extB = new Float64Array(m + 1);
    for (let i = 0; i < n; i++) extA[i] = GE * (1 - 0.85 * A.gap[i]);
    for (let j = 0; j < m; j++) extB[j] = GE * (1 - 0.85 * B.gap[j]);

    prevM[0] = 0; prevX[0] = NEG; prevY[0] = NEG;
    for (let j = 1; j <= m; j++) {
      prevM[j] = NEG; prevY[j] = NEG;
      prevX[j] = (j === 1 ? openB[0] : prevX[j - 1] + extB[j - 1]);
      tb[j] = 1;                                   // came from the left
    }
    for (let i = 1; i <= n; i++) {
      curM[0] = NEG; curX[0] = NEG;
      curY[0] = (i === 1 ? openA[0] : prevY[0] + extA[i - 1]);
      tb[i * (m + 1)] = 2;                         // came from above
      for (let j = 1; j <= m; j++) {
        const sc = columnScore(A, B, i - 1, j - 1, M);
        /* diagonal */
        let bestD = prevM[j - 1], fromD = 0;
        if (prevX[j - 1] > bestD) { bestD = prevX[j - 1]; fromD = 1; }
        if (prevY[j - 1] > bestD) { bestD = prevY[j - 1]; fromD = 2; }
        curM[j] = bestD + sc;
        /* gap in A: we advance along B */
        const openX = curM[j - 1] + openB[j - 1], extXv = curX[j - 1] + extB[j - 1];
        const isExtX = extXv > openX;
        curX[j] = isExtX ? extXv : openX;
        /* gap in B: we advance along A */
        const openY = prevM[j] + openA[i - 1], extYv = prevY[j] + extA[i - 1];
        const isExtY = extYv > openY;
        curY[j] = isExtY ? extYv : openY;
        tb[i * (m + 1) + j] = fromD | (isExtX ? 4 : 0) | (isExtY ? 8 : 0);
      }
      prevM.set(curM); prevX.set(curX); prevY.set(curY);
    }
    /* traceback from the best of the three states at the corner */
    let i = n, j = m;
    let state = 0, best = prevM[m];
    if (prevX[m] > best) { best = prevX[m]; state = 1; }
    if (prevY[m] > best) { best = prevY[m]; state = 2; }
    const path = [];                                 // 0 = both, 1 = gap in A, 2 = gap in B
    while (i > 0 || j > 0) {
      if (i === 0) { path.push(1); j--; continue; }
      if (j === 0) { path.push(2); i--; continue; }
      const cell = tb[i * (m + 1) + j];
      if (state === 0) {
        path.push(0);
        state = cell & 3;
        i--; j--;
      } else if (state === 1) {
        path.push(1);
        state = (cell & 4) ? 1 : 0;
        j--;
      } else {
        path.push(2);
        state = (cell & 8) ? 2 : 0;
        i--;
      }
    }
    path.reverse();
    return { path, score: best };
  }

  /* apply an alignment path to the two sets of sequences */
  function mergeByPath(A, B, path) {
    const outA = A.seqs.map(() => []), outB = B.seqs.map(() => []);
    let i = 0, j = 0;
    for (let p = 0; p < path.length; p++) {
      const step = path[p];
      if (step === 0) {
        A.seqs.forEach((s, k) => outA[k].push(s[i]));
        B.seqs.forEach((s, k) => outB[k].push(s[j]));
        i++; j++;
      } else if (step === 1) {
        A.seqs.forEach((s, k) => outA[k].push('-'));
        B.seqs.forEach((s, k) => outB[k].push(s[j]));
        j++;
      } else {
        A.seqs.forEach((s, k) => outA[k].push(s[i]));
        B.seqs.forEach((s, k) => outB[k].push('-'));
        i++;
      }
    }
    return { rows: A.rows.concat(B.rows), seqs: outA.map(a => a.join('')).concat(outB.map(a => a.join(''))) };
  }

  /* ================================================================
     4 · the progressive alignment itself
     ================================================================ */
  /* The defaults come from the benchmark in validation/block2, which scores the
     aligner against the true alignment that generated the simulated data:
       rbcL (coding, few indels) ...... 99.98 % of the true residue pairs
       trnL-F (spacer, many indels) ... 98.14 %
       ITS (nuclear, fast) ............ 98.93 %
       protein (MYB family) ........... 99.60 %
       rbcL aligned by codon .......... 100.00 %
     A word on `refine`: iterative refinement always improves the sum-of-pairs
     SCORE, but it optimises the scoring model, not the truth — here it gained
     0.3 points on the indel-rich spacer, changed nothing elsewhere, and costs
     three to four times the running time. One round is the compromise. */
  function defaults(type, opts) {
    opts = opts || {};
    return Object.assign({
      gapOpen: type === 'aa' ? -14 : -20,
      gapExt: -1,
      refine: 1,
      k: type === 'aa' ? 3 : 5,
    }, opts);
  }

  function progressive(seqs, type, opts) {
    const t0 = (typeof performance !== 'undefined' ? performance.now() : Date.now());
    opts = defaults(type, opts);
    const alpha = alphabetOf(type), map = indexer(alpha);
    const M = matrixOf(type, opts);
    const bare = seqs.map(s => String(s).toUpperCase().replace(/[-.\s]/g, ''));
    const n = bare.length;
    if (n < 2) return { aligned: bare.slice(), guide: null, score: 0, ms: 0 };

    const D = kmerDistance(bare, type, opts.k);
    const guide = upgma(D, opts.labels);

    const scoreOpts = { matrix: M, gapOpen: opts.gapOpen, gapExt: opts.gapExt };
    let steps = 0;
    const build = node => {
      if (node.tip != null) return makeProfile([node.tip], [bare[node.tip]], alpha, map);
      const A = build(node.children[0]), B = build(node.children[1]);
      const { path } = alignProfiles(A, B, scoreOpts);
      const merged = mergeByPath(A, B, path);
      steps++;
      if (opts.progress) opts.progress(steps / (n - 1));
      return makeProfile(merged.rows, merged.seqs, alpha, map);
    };
    let prof = build(guide);

    /* put the sequences back in their original order */
    const aligned = new Array(n);
    prof.rows.forEach((r, i) => { aligned[r] = prof.seqs[i]; });

    let score = spScore(aligned, type, opts);
    let refined = 0;
    for (let round = 0; round < (opts.refine || 0); round++) {
      const res = refineOnce(aligned, guide, type, opts, score);
      if (res.improved) { score = res.score; refined++; for (let i = 0; i < n; i++) aligned[i] = res.aligned[i]; }
      else break;
    }
    return {
      aligned, guide, score, refined,
      ms: (typeof performance !== 'undefined' ? performance.now() : Date.now()) - t0,
      length: aligned[0].length,
    };
  }

  /* one round of refinement: every split of the guide tree is tried */
  function refineOnce(aligned, guide, type, opts, curScore) {
    const alpha = alphabetOf(type), map = indexer(alpha);
    const M = matrixOf(type, opts);
    const scoreOpts = { matrix: M, gapOpen: opts.gapOpen, gapExt: opts.gapExt };
    const n = aligned.length;
    const splits = [];
    (function walk(node) {
      if (node.tip != null) return [node.tip];
      const kids = node.children.map(walk);
      const set = kids[0].concat(kids[1]);
      if (kids[0].length && kids[0].length < n) splits.push(kids[0]);
      return set;
    })(guide);

    let best = aligned.slice(), bestScore = curScore, improved = false;
    splits.forEach(group => {
      const inSet = new Set(group);
      const rowsA = [], rowsB = [];
      for (let i = 0; i < n; i++) (inSet.has(i) ? rowsA : rowsB).push(i);
      if (!rowsA.length || !rowsB.length) return;
      const subA = squeeze(rowsA.map(i => best[i]));
      const subB = squeeze(rowsB.map(i => best[i]));
      const A = makeProfile(rowsA, subA, alpha, map);
      const B = makeProfile(rowsB, subB, alpha, map);
      const { path } = alignProfiles(A, B, scoreOpts);
      const merged = mergeByPath(A, B, path);
      const cand = new Array(n);
      merged.rows.forEach((r, i) => { cand[r] = merged.seqs[i]; });
      const sc = spScore(cand, type, opts);
      if (sc > bestScore + 1e-9) { bestScore = sc; best = cand; improved = true; }
    });
    return { aligned: best, score: bestScore, improved };
  }

  /* drop the columns that are all gaps */
  function squeeze(seqs) {
    if (!seqs.length) return seqs;
    const L = seqs[0].length, keep = [];
    for (let c = 0; c < L; c++) {
      let allGap = true;
      for (let i = 0; i < seqs.length; i++) if (seqs[i][c] !== '-') { allGap = false; break; }
      if (!allGap) keep.push(c);
    }
    if (keep.length === L) return seqs.slice();
    return seqs.map(s => keep.map(c => s[c]).join(''));
  }

  /* ================================================================
     5 · scoring an alignment
     ================================================================ */
  /* sum of pairs, with affine gaps counted once per gap run and per pair */
  function spScore(aligned, type, opts) {
    opts = defaults(type, opts);
    const M = matrixOf(type, opts), map = indexer(alphabetOf(type));
    const n = aligned.length, L = aligned[0] ? aligned[0].length : 0;
    let s = 0;
    for (let i = 0; i < n; i++) for (let j = i + 1; j < n; j++) {
      const a = aligned[i], b = aligned[j];
      let inGap = false;
      for (let c = 0; c < L; c++) {
        const ca = a[c], cb = b[c];
        const ga = ca === '-', gb = cb === '-';
        if (ga && gb) continue;
        if (ga || gb) {
          s += inGap ? opts.gapExt : opts.gapOpen;
          inGap = true;
        } else {
          inGap = false;
          const x = map[a.charCodeAt(c)], y = map[b.charCodeAt(c)];
          if (x >= 0 && y >= 0) s += M[x][y];
        }
      }
    }
    return s;
  }

  /* how much of a reference alignment a test alignment recovers:
     SP = fraction of the aligned residue pairs of the reference that the test
     alignment also has; TC = fraction of columns recovered exactly.
     This is the standard way alignment programs are benchmarked. */
  function compareToReference(test, ref) {
    const n = ref.length;
    const pairsOf = aln => {
      /* for every sequence, the index of the residue in each column */
      const idx = aln.map(s => {
        const a = new Int32Array(s.length).fill(-1);
        let k = 0;
        for (let c = 0; c < s.length; c++) a[c] = (s[c] === '-') ? -1 : k++;
        return a;
      });
      const set = new Set();
      for (let i = 0; i < aln.length; i++) for (let j = i + 1; j < aln.length; j++) {
        const L = aln[i].length;
        for (let c = 0; c < L; c++) {
          const x = idx[i][c], y = idx[j][c];
          if (x >= 0 && y >= 0) set.add(i + ',' + j + ',' + x + ',' + y);
        }
      }
      return set;
    };
    const R = pairsOf(ref), T = pairsOf(test);
    let shared = 0;
    R.forEach(p => { if (T.has(p)) shared++; });
    /* columns recovered exactly */
    const colKey = (aln, c) => {
      const parts = [];
      for (let i = 0; i < aln.length; i++) parts.push(aln[i][c] === '-' ? '-' : '');
      return parts.join('');
    };
    const colsOf = aln => {
      const idx = aln.map(s => { const a = []; let k = 0; for (let c = 0; c < s.length; c++) a.push(s[c] === '-' ? -1 : k++); return a; });
      const out = new Set();
      for (let c = 0; c < aln[0].length; c++) {
        const sig = idx.map(a => a[c]).join(',');
        if (/[0-9]/.test(sig)) out.add(sig);
      }
      return out;
    };
    const RC = colsOf(ref), TC = colsOf(test);
    let sharedCols = 0;
    RC.forEach(c => { if (TC.has(c)) sharedCols++; });
    return {
      sp: R.size ? shared / R.size : 1,
      tc: RC.size ? sharedCols / RC.size : 1,
      refPairs: R.size, testPairs: T.size, refCols: RC.size,
    };
  }

  /* ================================================================
     6 · translation and codon-aware alignment
     ================================================================ */
  /* NCBI translation table 1 (the standard code), in the order T C A G for each
     of the three positions. The genetic code is a fact of nature and the NCBI
     tables are a work of the United States government, in the public domain.
     Table 11 (bacterial and plastid) has the same
     amino acids and differs only in the start codons, so plant chloroplast and
     mitochondrial genes use this same table. */
  const CODE1 = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG';
  const ORDER = { T: 0, U: 0, C: 1, A: 2, G: 3 };
  function translateCodon(c) {
    const a = ORDER[c[0]], b = ORDER[c[1]], d = ORDER[c[2]];
    if (a == null || b == null || d == null) return 'X';
    return CODE1[a * 16 + b * 4 + d];
  }
  function translate(seq, frame) {
    frame = (frame || 1) - 1;
    const s = String(seq).toUpperCase().replace(/[\s-]/g, '');
    let out = '';
    for (let i = frame; i + 2 < s.length; i += 3) out += translateCodon(s.slice(i, i + 3));
    return out;
  }
  /* where are the stop codons, ignoring a final one */
  function stopCodons(seq, frame) {
    const aa = translate(seq, frame);
    const out = [];
    for (let i = 0; i < aa.length; i++) if (aa[i] === '*' && i < aa.length - 1) out.push(i + 1);
    return out;
  }
  /* the reading frame with the fewest internal stops */
  function bestFrame(seqs) {
    let best = 1, bestStops = Infinity;
    for (let f = 1; f <= 3; f++) {
      let stops = 0;
      seqs.forEach(s => { stops += stopCodons(s, f).length; });
      if (stops < bestStops) { bestStops = stops; best = f; }
    }
    return { frame: best, stops: bestStops };
  }

  function codonAlign(seqs, opts) {
    opts = opts || {};
    const frame = opts.frame || bestFrame(seqs).frame;
    const bare = seqs.map(s => String(s).toUpperCase().replace(/[-.\s]/g, '').slice(frame - 1));
    const prot = bare.map(s => translate(s, 1));
    const res = progressive(prot, 'aa', Object.assign({}, opts, { refine: opts.refine == null ? 2 : opts.refine }));
    /* map the protein alignment back onto the nucleotides, three at a time */
    const aligned = res.aligned.map((p, i) => {
      const nt = bare[i];
      let k = 0, out = '';
      for (let c = 0; c < p.length; c++) {
        if (p[c] === '-') out += '---';
        else { out += (nt.slice(k * 3, k * 3 + 3) + '---').slice(0, 3); k++; }
      }
      /* a trailing incomplete codon is kept at the end */
      const rest = nt.slice(k * 3);
      if (rest) out += rest;
      return out;
    });
    const L = Math.max(...aligned.map(s => s.length));
    return Object.assign({}, res, {
      aligned: aligned.map(s => s + '-'.repeat(L - s.length)),
      protein: res.aligned, frame, length: L,
    });
  }

  /* ================================================================
     7 · description of an alignment
     ================================================================ */
  function stats(aligned, type) {
    const n = aligned.length, L = aligned[0] ? aligned[0].length : 0;
    let gaps = 0, variable = 0, informative = 0, constant = 0, allGapCols = 0;
    const missing = type === 'aa' ? 'X?' : 'N?';
    for (let c = 0; c < L; c++) {
      const counts = {};
      let nGap = 0, valid = 0;
      for (let i = 0; i < n; i++) {
        const ch = aligned[i][c];
        if (ch === '-') { nGap++; gaps++; continue; }
        if (missing.indexOf(ch) >= 0) continue;
        counts[ch] = (counts[ch] || 0) + 1;
        valid++;
      }
      if (nGap === n) allGapCols++;
      const keys = Object.keys(counts);
      if (keys.length > 1) {
        variable++;
        if (keys.filter(k => counts[k] >= 2).length >= 2) informative++;
      } else if (keys.length === 1 && valid > 0) constant++;
    }
    /* mean pairwise identity over the columns both sequences have */
    let idSum = 0, idPairs = 0;
    for (let i = 0; i < n; i++) for (let j = i + 1; j < n; j++) {
      let same = 0, tot = 0;
      for (let c = 0; c < L; c++) {
        const a = aligned[i][c], b = aligned[j][c];
        if (a === '-' || b === '-') continue;
        tot++; if (a === b) same++;
      }
      if (tot) { idSum += same / tot; idPairs++; }
    }
    return {
      n, length: L, gapFraction: n * L ? gaps / (n * L) : 0,
      variable, informative, constant, allGapCols,
      identity: idPairs ? idSum / idPairs : 1,
    };
  }

  Object.assign(Align, {
    DNA, AA, BLOSUM62_ROWS, dnaMatrix, alphabetOf, matrixOf, indexer,
    kmerDistance, upgma, makeProfile, alignProfiles, mergeByPath, squeeze,
    progressive, refineOnce, spScore, compareToReference,
    translate, translateCodon, stopCodons, bestFrame, codonAlign, stats, defaults,
  });
  window.Align = Align;
})();
