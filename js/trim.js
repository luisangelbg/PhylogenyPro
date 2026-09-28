/* PhylogenyPro — molecular phylogenetics in the browser.
   Copyright (C) 2026  Luis Ángel Barrera-Guzmán

   This program is free software: you can redistribute it and/or modify it under
   the terms of the GNU General Public License as published by the Free Software
   Foundation, either version 3 of the License, or (at your option) any later
   version. It is distributed in the hope that it will be useful, but WITHOUT ANY
   WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
   A PARTICULAR PURPOSE. See the GNU General Public License, in the file LICENSE
   at the root of this program, or <https://www.gnu.org/licenses/>. */

/* PhylogenyPro — trimming of ambiguous blocks and quality control.

   Trimming
     · Gblocks (Castresana 2000): positions are classified as non-conserved,
       conserved or highly conserved by how many sequences share the commonest
       residue; stretches of contiguous non-conserved positions longer than b3
       are dropped, the surviving blocks are flanked by highly conserved
       positions and must be at least b4 long, and gap-rich columns are handled
       by b5. The defaults are the ones the program itself uses.
     · trimAl-style thresholds: drop columns by gap fraction, by mean pairwise
       similarity, or by both.
   Both return a mask over the columns, so the user sees what would be removed
   before anything is removed.

   Quality control
     Everything that quietly ruins a phylogenetic analysis and is invisible in a
     spreadsheet: identical sequences, taxa that are almost all gaps, sequences
     of very different lengths, unexpected characters, internal stop codons,
     columns that are all gaps, and taxa missing from one partition. */

const Trim = {};
const QC = {};

(function () {

  const isGap = c => c === '-' || c === '?' || c === 'N' || c === 'X' || c === '.';

  /* ================================================================
     column statistics, computed once and reused
     ================================================================ */
  function columnStats(seqs) {
    const n = seqs.length, L = seqs[0] ? seqs[0].length : 0;
    const out = [];
    for (let c = 0; c < L; c++) {
      const counts = {};
      let gaps = 0;
      for (let i = 0; i < n; i++) {
        const ch = seqs[i][c];
        if (ch === '-') { gaps++; continue; }
        if (ch === '?' || ch === 'N' || ch === 'X') { gaps++; continue; }
        counts[ch] = (counts[ch] || 0) + 1;
      }
      let best = 0, states = 0;
      for (const k in counts) { states++; if (counts[k] > best) best = counts[k]; }
      out.push({ gaps, best, states, valid: n - gaps, counts });
    }
    return out;
  }

  /* ================================================================
     Gblocks
     ================================================================ */
  function gblocks(seqs, opts) {
    opts = opts || {};
    const n = seqs.length, L = seqs[0] ? seqs[0].length : 0;
    const b1 = opts.b1 || Math.floor(n / 2) + 1;                 // conserved
    const b2 = opts.b2 || Math.ceil(n * 0.85);                   // highly conserved
    const b3 = opts.b3 == null ? 8 : opts.b3;                    // max contiguous non-conserved
    const b4 = opts.b4 == null ? 10 : opts.b4;                   // minimum block length
    const b5 = opts.b5 || 'half';                                // gaps: 'none' | 'half' | 'all'
    const st = opts.stats || columnStats(seqs);

    /* 1 · classify every position */
    const KIND = new Uint8Array(L);          // 0 non-conserved, 1 conserved, 2 highly conserved
    const gapOK = new Uint8Array(L);
    for (let c = 0; c < L; c++) {
      const s = st[c];
      const gapFrac = s.gaps / n;
      gapOK[c] = b5 === 'all' ? 1 : b5 === 'half' ? (gapFrac <= 0.5 ? 1 : 0) : (s.gaps === 0 ? 1 : 0);
      if (!gapOK[c]) { KIND[c] = 0; continue; }
      KIND[c] = s.best >= b2 ? 2 : (s.best >= b1 ? 1 : 0);
    }

    /* 2 · a position is kept only if it is conserved and not inside a long
           stretch of non-conserved positions */
    const keep = new Uint8Array(L);
    for (let c = 0; c < L; c++) keep[c] = KIND[c] > 0 ? 1 : 0;
    let c0 = 0;
    while (c0 < L) {
      if (keep[c0]) { c0++; continue; }
      let c1 = c0;
      while (c1 < L && !keep[c1]) c1++;
      /* a short stretch between two kept regions survives; a long one does not */
      const len = c1 - c0;
      if (len <= b3 && c0 > 0 && c1 < L) for (let c = c0; c < c1; c++) keep[c] = 1;
      c0 = c1;
    }

    /* 3 · the flanks of every block must be highly conserved */
    c0 = 0;
    while (c0 < L) {
      if (!keep[c0]) { c0++; continue; }
      let c1 = c0;
      while (c1 < L && keep[c1]) c1++;
      let a = c0, b = c1 - 1;
      while (a <= b && KIND[a] !== 2) a++;
      while (b >= a && KIND[b] !== 2) b--;
      for (let c = c0; c < a; c++) keep[c] = 0;
      for (let c = b + 1; c < c1; c++) keep[c] = 0;
      c0 = c1;
    }

    /* 4 · blocks shorter than b4 are dropped */
    c0 = 0;
    while (c0 < L) {
      if (!keep[c0]) { c0++; continue; }
      let c1 = c0;
      while (c1 < L && keep[c1]) c1++;
      if (c1 - c0 < b4) for (let c = c0; c < c1; c++) keep[c] = 0;
      c0 = c1;
    }

    return { keep, params: { b1, b2, b3, b4, b5 }, kept: count(keep), L };
  }

  /* ================================================================
     trimAl-style thresholds
     ================================================================ */
  function byGaps(seqs, maxGapFraction, stats) {
    const n = seqs.length, L = seqs[0] ? seqs[0].length : 0;
    const st = stats || columnStats(seqs);
    const keep = new Uint8Array(L);
    for (let c = 0; c < L; c++) keep[c] = (st[c].gaps / n) <= maxGapFraction ? 1 : 0;
    return { keep, kept: count(keep), L };
  }
  /* mean pairwise identity of a column, over the sequences that have a residue */
  function similarity(seqs, stats) {
    const n = seqs.length, L = seqs[0] ? seqs[0].length : 0;
    const st = stats || columnStats(seqs);
    const sim = new Float64Array(L);
    for (let c = 0; c < L; c++) {
      const s = st[c];
      let same = 0, pairs = 0;
      const vals = Object.values(s.counts);
      vals.forEach(v => { same += v * (v - 1) / 2; });
      pairs = s.valid * (s.valid - 1) / 2;
      sim[c] = pairs ? same / pairs : 0;
    }
    return sim;
  }
  function bySimilarity(seqs, minSim, stats) {
    const sim = similarity(seqs, stats);
    const keep = new Uint8Array(sim.length);
    for (let c = 0; c < sim.length; c++) keep[c] = sim[c] >= minSim ? 1 : 0;
    return { keep, sim, kept: count(keep), L: sim.length };
  }
  function count(mask) { let k = 0; for (let i = 0; i < mask.length; i++) if (mask[i]) k++; return k; }
  function intersect(a, b) {
    const out = new Uint8Array(a.length);
    for (let i = 0; i < a.length; i++) out[i] = (a[i] && b[i]) ? 1 : 0;
    return out;
  }
  /* keep whole codons: a codon survives only if its three positions do */
  function codonAware(keep, frame) {
    const out = new Uint8Array(keep.length);
    const off = (frame || 1) - 1;
    for (let c = off; c + 2 < keep.length; c += 3) {
      const all = keep[c] && keep[c + 1] && keep[c + 2];
      out[c] = out[c + 1] = out[c + 2] = all ? 1 : 0;
    }
    return out;
  }
  function apply(seqs, keep) {
    const idx = [];
    for (let c = 0; c < keep.length; c++) if (keep[c]) idx.push(c);
    return seqs.map(s => { let o = ''; for (let i = 0; i < idx.length; i++) o += s[idx[i]]; return o; });
  }
  /* the blocks that survive, as ranges, for the report */
  function blocks(keep) {
    const out = [];
    let c = 0;
    while (c < keep.length) {
      if (!keep[c]) { c++; continue; }
      const start = c;
      while (c < keep.length && keep[c]) c++;
      out.push([start + 1, c]);
    }
    return out;
  }

  Object.assign(Trim, { columnStats, gblocks, byGaps, bySimilarity, similarity, intersect, codonAware, apply, blocks, count });

  /* ================================================================
     QUALITY CONTROL
     ================================================================ */
  /* Each finding is {level, es, en, rows?} — level 'bad' stops the analysis
     being trustworthy, 'warn' deserves a look, 'info' is just information. */
  function check(taxa, seqs, type, opts) {
    opts = opts || {};
    const out = [];
    const n = taxa.length, L = seqs[0] ? seqs[0].length : 0;
    const aligned = seqs.every(s => s.length === L);

    /* --- names --- */
    const seen = {}, dupes = [];
    taxa.forEach(t => { const k = t.trim().toLowerCase(); if (seen[k]) dupes.push(t); else seen[k] = 1; });
    if (dupes.length) out.push({ level: 'bad', rows: dupes,
      es: `Hay ${dupes.length} nombre(s) repetido(s): ${dupes.slice(0, 4).join(', ')}${dupes.length > 4 ? '…' : ''}.`,
      en: `There are ${dupes.length} repeated name(s): ${dupes.slice(0, 4).join(', ')}${dupes.length > 4 ? '…' : ''}.` });
    const odd = taxa.filter(t => /[^\w .\-]/.test(t));
    if (odd.length) out.push({ level: 'warn', rows: odd,
      es: `${odd.length} nombre(s) traen caracteres que algunos programas no aceptan (paréntesis, comas, acentos). Se pueden limpiar aquí.`,
      en: `${odd.length} name(s) contain characters some programs reject (brackets, commas, accents). They can be cleaned here.` });

    /* --- lengths --- */
    const lens = seqs.map(s => s.replace(/-/g, '').length);
    const maxL = Math.max(...lens), minL = Math.min(...lens);
    if (minL < maxL * 0.5) {
      const short = taxa.filter((t, i) => lens[i] < maxL * 0.5);
      out.push({ level: 'warn', rows: short,
        es: `${short.length} secuencia(s) miden menos de la mitad que la más larga (${minL} contra ${maxL} bases). Los datos faltantes en bloque distorsionan las distancias y los soportes.`,
        en: `${short.length} sequence(s) are less than half the length of the longest (${minL} against ${maxL} bases). Block-shaped missing data distort distances and supports.` });
    }
    if (lens.some(l => l === 0)) out.push({ level: 'bad',
      es: 'Hay al menos una secuencia vacía.', en: 'There is at least one empty sequence.' });

    /* --- identical sequences --- */
    if (aligned) {
      const map = new Map();
      seqs.forEach((s, i) => { const k = s; if (!map.has(k)) map.set(k, []); map.get(k).push(taxa[i]); });
      const groups = [...map.values()].filter(g => g.length > 1);
      if (groups.length) out.push({ level: 'info', rows: groups.map(g => g.join(' = ')),
        es: `${groups.length} grupo(s) de secuencias idénticas. No es un error, pero sus ramas quedarán en cero y algunos análisis las colapsan.`,
        en: `${groups.length} group(s) of identical sequences. Not an error, but their branches will be zero and some analyses collapse them.` });
    }

    /* --- unexpected characters --- */
    const allowed = type === 'aa' ? /[ACDEFGHIKLMNPQRSTVWYBZXJUO*\-?.]/i
      : type === 'morph' ? /[0-9A-Z?\-()]/i : /[ACGTUMRWSYKVHDBNX\-?.]/i;
    const badChars = {};
    seqs.forEach((s, i) => {
      for (let c = 0; c < s.length; c++) if (!allowed.test(s[c])) badChars[s[c]] = (badChars[s[c]] || 0) + 1;
    });
    const bk = Object.keys(badChars);
    if (bk.length) out.push({ level: 'warn',
      es: `Caracteres no esperados para estos datos: ${bk.map(k => `«${k}» (${badChars[k]})`).join(', ')}.`,
      en: `Characters not expected for these data: ${bk.map(k => `"${k}" (${badChars[k]})`).join(', ')}.` });

    /* --- gaps and missing data --- */
    if (aligned && L) {
      const perTaxon = seqs.map(s => (s.match(/[-?]/g) || []).length / L);
      const gappy = taxa.filter((t, i) => perTaxon[i] > 0.5);
      if (gappy.length) out.push({ level: 'warn', rows: gappy,
        es: `${gappy.length} taxón(es) tienen más de la mitad de huecos o datos faltantes.`,
        en: `${gappy.length} taxon(a) are more than half gaps or missing data.` });
      const st = columnStats(seqs);
      const allGap = st.filter(s => s.valid === 0).length;
      if (allGap) out.push({ level: 'warn',
        es: `${allGap} columna(s) no tienen ningún residuo; conviene quitarlas.`,
        en: `${allGap} column(s) have no residue at all; they should be removed.` });
      const constant = st.filter(s => s.states === 1).length;
      out.push({ level: 'info',
        es: `${constant} de ${L} columnas son constantes (${fmtPct(constant / L)}); ${st.filter(s => s.states > 1).length} son variables.`,
        en: `${constant} of ${L} columns are constant (${fmtPct(constant / L)}); ${st.filter(s => s.states > 1).length} are variable.` });
    }

    /* --- coding sequences --- */
    if (type === 'dna' && opts.coding) {
      const bf = Align.bestFrame(seqs);
      const stops = seqs.map((s, i) => ({ t: taxa[i], k: Align.stopCodons(s, opts.frame || bf.frame).length })).filter(x => x.k > 0);
      if (stops.length) out.push({ level: 'bad', rows: stops.map(x => `${x.t} (${x.k})`),
        es: `${stops.length} secuencia(s) tienen codones de paro internos en el marco ${opts.frame || bf.frame}. O el marco no es ése, o hay errores de secuenciación, o no todas son codificantes.`,
        en: `${stops.length} sequence(s) have internal stop codons in frame ${opts.frame || bf.frame}. Either that is not the frame, or there are sequencing errors, or not all of them are coding.` });
      else out.push({ level: 'info',
        es: `Ninguna secuencia tiene codones de paro internos en el marco ${opts.frame || bf.frame}.`,
        en: `No sequence has internal stop codons in frame ${opts.frame || bf.frame}.` });
      if (L % 3 !== 0 && aligned) out.push({ level: 'warn',
        es: `La longitud del alineamiento (${L}) no es múltiplo de tres; el último codón queda incompleto.`,
        en: `The alignment length (${L}) is not a multiple of three; the last codon is incomplete.` });
    }

    /* --- composition, which matters for the model in Block 3 --- */
    if (type === 'dna' && aligned) {
      const gc = seqs.map(s => {
        let g = 0, t = 0;
        for (let c = 0; c < s.length; c++) { const ch = s[c]; if (ch === 'G' || ch === 'C') { g++; t++; } else if (ch === 'A' || ch === 'T' || ch === 'U') t++; }
        return t ? g / t : 0;
      });
      const mn = Math.min(...gc), mx = Math.max(...gc);
      if (mx - mn > 0.1) out.push({ level: 'warn',
        es: `El contenido de G+C va de ${fmtPct(mn)} a ${fmtPct(mx)} entre secuencias. Una composición desigual rompe el supuesto de los modelos habituales; el Bloque 3 lo pone a prueba.`,
        en: `G+C content ranges from ${fmtPct(mn)} to ${fmtPct(mx)} between sequences. Unequal composition breaks the assumption of the usual models; Block 3 tests it.` });
    }
    return out;
  }

  /* taxa that are not in every partition */
  function missingTaxa(parts) {
    const all = new Set();
    parts.forEach(p => p.taxa.forEach(t => all.add(t)));
    const rows = [];
    all.forEach(t => {
      const missing = parts.filter(p => p.taxa.indexOf(t) < 0).map(p => p.name);
      if (missing.length) rows.push({ taxon: t, missing });
    });
    return { all: [...all], rows };
  }

  Object.assign(QC, { check, missingTaxa });
  window.Trim = Trim;
  window.QC = QC;
})();
