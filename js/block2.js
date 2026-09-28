/* PhylogenyPro — molecular phylogenetics in the browser.
   Copyright (C) 2026  Luis Ángel Barrera-Guzmán

   This program is free software: you can redistribute it and/or modify it under
   the terms of the GNU General Public License as published by the Free Software
   Foundation, either version 3 of the License, or (at your option) any later
   version. It is distributed in the hope that it will be useful, but WITHOUT ANY
   WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
   A PARTICULAR PURPOSE. See the GNU General Public License, in the file LICENSE
   at the root of this program, or <https://www.gnu.org/licenses/>. */

/* PhylogenyPro — Block 2: data and alignment.

   Holds the partitions, drives the aligner, the viewer, the trimming and the
   quality control, and writes state.data, which every later block reads. */

(function () {

  const B2 = {
    parts: [],          // [{name, genome, type, coding, frame, taxa, seqs, aligned, source}]
    active: -1,
    view: null,
    trim: null,         // current preview mask
    charsets: null,     // when a concatenated matrix is the active partition
  };
  window.B2 = B2;

  const GENOMES = {
    cp: ['cloroplasto', 'chloroplast'], mt: ['mitocondria', 'mitochondrion'],
    nr: ['nuclear ribosomal', 'nuclear ribosomal'], nu: ['nuclear', 'nuclear'], other: ['otro', 'other'],
  };
  const TYPES = { dna: ['ADN', 'DNA'], aa: ['proteína', 'protein'], morph: ['morfología', 'morphology'] };

  /* An <option> cannot hold the two spans that L2 writes: the browser keeps only
     the text and glues the two languages together ("nuclearnuclear"). Options
     carry both texts as attributes, which is what I18N.apply() reads, and start
     out written in the language now on screen. Returns the end of the open tag. */
  function opt(pair) {
    return ` data-es="${esc(pair[0])}" data-en="${esc(pair[1])}">${esc(T(pair[0], pair[1]))}`;
  }

  /* the genome a marker belongs to, guessed from its name — the user can change it */
  function guessGenome(name) {
    const n = String(name).toLowerCase();
    if (/(rbcl|matk|trn|ndh|rpl|rps|atpb|ycf|psb|pet|ccsa|cema|clpp|infa|accd)/.test(n)) return 'cp';
    if (/(its|ets|18s|26s|5\.8s|nrdna)/.test(n)) return 'nr';
    if (/(matr|nad\d|cox\d|atp1|atp9|ccmb|rps3m|mtdna)/.test(n)) return 'mt';
    if (/(waxy|phyc|leafy|g3pdh|gapdh|adh|ncpgs|pgic|rpb2|single.?copy|nuclear)/.test(n)) return 'nu';
    return 'nu';
  }

  /* ================================================================
     partitions
     ================================================================ */
  function addPartition(p) {
    const part = Object.assign({
      genome: guessGenome(p.name), coding: false, frame: 1, derived: false,
    }, p);
    part.aligned = part.seqs.every(s => s.length === part.seqs[0].length);
    B2.parts.push(part);
    return B2.parts.length - 1;
  }

  function renderParts() {
    const host = el('partList');
    if (!host) return;
    host.innerHTML = '';
    B2.parts.forEach((p, i) => {
      const row = mk('div', { class: 'part-row' + (i === B2.active ? ' on' : '') });
      const nSites = p.seqs.length ? Math.max(...p.seqs.map(s => s.length)) : 0;
      row.innerHTML = `
        <div>
          <div class="p-name">${esc(p.name)}
            ${p.aligned ? `<span class="badge ok" style="font-size:.66rem">${L2('alineado', 'aligned')}</span>`
                        : `<span class="badge" style="font-size:.66rem;color:var(--warning)">${L2('sin alinear', 'unaligned')}</span>`}
            ${p.derived ? `<span class="badge" style="font-size:.66rem">${L2('derivada', 'derived')}</span>` : ''}</div>
          <div class="p-meta">${p.taxa.length} ${L2('taxones', 'taxa')} · ${nSites} ${p.type === 'morph' ? L2('caracteres', 'characters') : L2('posiciones', 'positions')}${p.source ? ' · ' + esc(p.source) : ''}</div>
        </div>
        <select data-role="genome" data-i="${i}">${Object.keys(GENOMES).map(g =>
          `<option value="${g}"${p.genome === g ? ' selected' : ''}${opt(GENOMES[g])}</option>`).join('')}</select>
        <select data-role="type" data-i="${i}">${Object.keys(TYPES).map(t =>
          `<option value="${t}"${p.type === t ? ' selected' : ''}${opt(TYPES[t])}</option>`).join('')}</select>
        <button class="btn btn-secondary btn-sm" data-role="open" data-i="${i}">${L2('Ver', 'Open')}</button>
        <button class="btn btn-ghost btn-sm" data-role="drop" data-i="${i}" title="${T('Quitar', 'Remove')}">✕</button>`;
      host.appendChild(row);
    });
    I18N.apply(host);
    host.querySelectorAll('[data-role]').forEach(nd => {
      const i = +nd.dataset.i, role = nd.dataset.role;
      if (role === 'genome' || role === 'type') {
        nd.addEventListener('change', () => {
          B2.parts[i][role] = nd.value;
          if (role === 'type') { B2.parts[i].coding = false; refreshView(); }
          commit();
        });
      } else if (role === 'open') nd.addEventListener('click', () => selectPart(i));
      else if (role === 'drop') nd.addEventListener('click', () => {
        B2.parts.splice(i, 1);
        if (B2.active >= B2.parts.length) B2.active = B2.parts.length - 1;
        renderParts(); refreshView(); commit();
      });
    });
    el('cardParts').style.display = B2.parts.length ? '' : 'none';
    ['cardAln', 'cardTrim', 'cardQC', 'cardExport'].forEach(id => {
      const n = el(id); if (n) n.style.display = B2.parts.length ? '' : 'none';
    });
  }

  function current() { return B2.parts[B2.active] || null; }

  function selectPart(i) {
    B2.active = i;
    B2.trim = null;
    renderParts();
    refreshView();
    runQC();
  }

  /* ================================================================
     the viewer
     ================================================================ */
  function ensureView() {
    if (B2.view) return B2.view;
    B2.view = AlnView.create(el('alnHost'), {
      onSelect: sel => {
        const info = el('alnSel');
        if (!info) return;
        if (sel.col < 0) { info.textContent = ''; return; }
        const p = current();
        const a = Math.min(sel.col, sel.col2), b = Math.max(sel.col, sel.col2);
        const range = a === b ? `${a + 1}` : `${a + 1}–${b + 1} (${b - a + 1})`;
        info.textContent = `${T('columna', 'column')} ${range}` +
          (sel.taxon ? ` · ${sel.taxon}${sel.residue ? ' = ' + sel.residue : ''}` : '');
      },
      onChange: () => {
        const p = current();
        if (!p) return;
        const g = B2.view.get();
        p.taxa = g.taxa; p.seqs = g.seqs;
        p.aligned = p.seqs.every(s => s.length === p.seqs[0].length);
        B2.trim = null;
        updateSummary(); updateTrim(); commit();
      },
    });
    return B2.view;
  }

  function refreshView() {
    const p = current();
    const host = el('alnHost');
    if (!host) return;
    const view = ensureView();
    if (!p) { view.setData([], [], 'dna'); return; }
    view.setData(p.taxa, p.seqs, p.type);
    view.setMode(el('alnMode') ? el('alnMode').value : 'base');
    applyMarks();
    view.setFrame(el('alnFrameChk') && el('alnFrameChk').checked ? (p.frame || 1) : null);
    view.setMask(null);
    updateSummary();
    updateTrim();
    const codon = el('alnCodon');
    if (codon) codon.disabled = p.type !== 'dna';
    const qcCoding = el('qcCoding');
    if (qcCoding) qcCoding.checked = !!p.coding;
  }

  function applyMarks() {
    const p = current();
    if (!p || !B2.view) return;
    if (!el('alnMarks') || !el('alnMarks').checked || !p.aligned) { B2.view.setMarks(null); return; }
    const st = Trim.columnStats(p.seqs);
    const marks = new Uint8Array(st.length);
    st.forEach((s, c) => {
      if (s.states > 1) {
        const informative = Object.values(s.counts).filter(v => v >= 2).length >= 2;
        marks[c] = informative ? 2 : 1;
      }
    });
    B2.view.setMarks(marks);
  }

  function updateSummary() {
    const p = current(), host = el('alnSummary');
    if (!host) return;
    if (!p) { host.textContent = ''; return; }
    if (!p.aligned) {
      const lens = p.seqs.map(s => s.length);
      host.textContent = `${p.taxa.length} ${T('secuencias', 'sequences')} · ${Math.min(...lens)}–${Math.max(...lens)} ${T('bases', 'bases')} · ${T('sin alinear', 'unaligned')}`;
      return;
    }
    const s = Align.stats(p.seqs, p.type);
    host.textContent = `${s.n} × ${s.length} · ${T('identidad', 'identity')} ${fmtPct(s.identity)} · ` +
      `${T('variables', 'variable')} ${s.variable} · ${T('informativos', 'informative')} ${s.informative} · ` +
      `${T('huecos', 'gaps')} ${fmtPct(s.gapFraction)}`;
  }

  /* ================================================================
     loading
     ================================================================ */
  function loadText(text, name) {
    const res = SeqIO.parse(text, name);
    const base = String(name || 'datos').replace(/\.[^.]+$/, '');
    const msgs = el('loadMessages');
    /* a NEXUS or a partition-carrying file may define several partitions at once */
    if (res.charsets && res.charsets.length > 1 && res.aligned) {
      res.charsets.forEach(cs => {
        addPartition({
          name: cs.name, taxa: res.taxa.slice(),
          seqs: res.seqs.map(s => s.slice(cs.from - 1, cs.to)),
          type: res.type, source: name,
        });
      });
      showMessage(msgs, 'success', L2(
        `Se leyeron ${res.taxa.length} taxones y ${res.charsets.length} particiones del archivo ${esc(name)}.`,
        `Read ${res.taxa.length} taxa and ${res.charsets.length} partitions from ${esc(name)}.`));
    } else {
      const i = addPartition({ name: base, taxa: res.taxa, seqs: res.seqs, type: res.type, source: name });
      const p = B2.parts[i];
      showMessage(msgs, 'success', L2(
        `«${esc(base)}»: ${res.taxa.length} secuencias, formato ${res.format.toUpperCase()}, ${p.aligned ? 'alineadas' : 'sin alinear'}.`,
        `"${esc(base)}": ${res.taxa.length} sequences, ${res.format.toUpperCase()} format, ${p.aligned ? 'aligned' : 'unaligned'}.`));
    }
    (res.warnings || []).forEach(w => showMessage(msgs, 'warning', L2(w.es, w.en)));
    if (res.trees && res.trees.length) {
      state.importedTrees = (state.importedTrees || []).concat(res.trees);
      showMessage(msgs, 'info', L2(
        `El archivo traía ${res.trees.length} árbol(es); se guardaron para los bloques que los usan.`,
        `The file carried ${res.trees.length} tree(s); they were kept for the blocks that use them.`));
    }
    return res;
  }

  function loadFiles(files) {
    const msgs = el('loadMessages');
    clearMessages(msgs);
    let pending = files.length;
    [...files].forEach(f => {
      const r = new FileReader();
      r.onload = () => {
        try { loadText(String(r.result), f.name); }
        catch (e) { showMessage(msgs, 'error', `${esc(f.name)}: ${esc(e.message)}`); }
        if (--pending === 0) { if (B2.active < 0) B2.active = 0; renderParts(); refreshView(); runQC(); commit(); }
      };
      r.onerror = () => { showMessage(msgs, 'error', `${esc(f.name)}: ${T('no se pudo leer', 'could not be read')}`); if (--pending === 0) { renderParts(); refreshView(); } };
      r.readAsText(f);
    });
  }

  function loadExample(kind) {
    const msgs = el('loadMessages');
    clearMessages(msgs);
    B2.parts.length = 0;
    const add = (key, name) => {
      const obj = EXAMPLES.seqs[key];
      addPartition({ name, taxa: Object.keys(obj), seqs: Object.values(obj), type: SeqIO.guessType(Object.values(obj)), source: 'ejemplo' });
    };
    if (kind === 'bursera3') {
      add('rbcL', 'rbcL'); add('trnLF', 'trnL-F'); add('ITS', 'ITS');
      B2.parts[0].coding = true;
      showMessage(msgs, 'success', L2(
        'Tres genes de <i>Bursera</i> sin alinear: rbcL (cloroplasto, codificante), trnL-F (cloroplasto, con indeles) e ITS (nuclear ribosomal). Alinéalos uno por uno y compáralos.',
        'Three <i>Bursera</i> genes, unaligned: rbcL (chloroplast, coding), trnL-F (chloroplast, with indels) and ITS (nuclear ribosomal). Align them one by one and compare.'));
    } else if (kind === 'burseraAln') {
      add('concatenado', 'concatenado');
      const parts = SeqIO.parsePartitionText(EXAMPLES.particiones);
      B2.charsets = parts;
      B2.parts[0].derived = true;
      showMessage(msgs, 'success', L2(
        `Matriz concatenada y ya alineada, con su esquema de ${parts.length} particiones (rbcL, trnL-F, ITS).`,
        `A concatenated, already aligned matrix, with its scheme of ${parts.length} partitions (rbcL, trnL-F, ITS).`));
      /* split it straight away, which is what the scheme is for */
      const whole = B2.parts[0];
      parts.forEach(cs => addPartition({
        name: cs.name, taxa: whole.taxa.slice(),
        seqs: whole.seqs.map(s => s.slice(cs.from - 1, cs.to)), type: whole.type, source: 'ejemplo',
      }));
      B2.parts[1].coding = true;
    } else if (kind === 'proteina') {
      add('proteina', 'familia MYB');
      showMessage(msgs, 'success', L2(
        'Veintidós secuencias de aminoácidos de una familia génica, sin alinear.',
        'Twenty-two amino-acid sequences of a gene family, unaligned.'));
    } else if (kind === 'morfologia') {
      const r = SeqIO.parse(EXAMPLES.morfologia, 'morfologia.csv');
      addPartition({ name: 'morfología', taxa: r.taxa, seqs: r.seqs, type: 'morph', source: 'ejemplo' });
      showMessage(msgs, 'success', L2(
        'Matriz morfológica de 25 taxones y 42 caracteres discretos, con algunos datos faltantes.',
        'A morphological matrix of 25 taxa and 42 discrete characters, with some missing entries.'));
    }
    B2.active = 0;
    renderParts(); refreshView(); runQC(); commit();
  }

  /* ================================================================
     aligning
     ================================================================ */
  function runAlign() {
    const p = current();
    if (!p) return;
    const msgs = el('alnMessages');
    clearMessages(msgs);
    if (p.type === 'morph') { showMessage(msgs, 'warning', L2('Una matriz morfológica no se alinea.', 'A morphological matrix is not aligned.')); return; }
    const opts = {
      gapOpen: +el('alnGapOpen').value, gapExt: +el('alnGapExt').value,
      refine: +el('alnRefine').value,
    };
    const codon = el('alnCodon').checked && p.type === 'dna';
    const btn = el('alignBtn');
    btn.disabled = true;
    const old = btn.innerHTML;
    btn.innerHTML = T('Alineando…', 'Aligning…');
    /* let the button repaint before the work starts */
    setTimeout(() => {
      try {
        const res = codon ? Align.codonAlign(p.seqs, Object.assign({ frame: p.frame || 1 }, opts))
          : Align.progressive(p.seqs, p.type, opts);
        p.seqs = res.aligned;
        p.aligned = true;
        if (codon) { p.coding = true; p.frame = res.frame; }
        const s = Align.stats(p.seqs, p.type);
        showMessage(msgs, 'success', L2(
          `Alineado en ${(res.ms / 1000).toFixed(1)} s: ${s.length} columnas, identidad media ${fmtPct(s.identity)}, ${fmtPct(s.gapFraction)} de huecos.` +
            (codon ? ` Marco de lectura ${res.frame}; ningún hueco parte un codón.` : ''),
          `Aligned in ${(res.ms / 1000).toFixed(1)} s: ${s.length} columns, mean identity ${fmtPct(s.identity)}, ${fmtPct(s.gapFraction)} gaps.` +
            (codon ? ` Reading frame ${res.frame}; no gap breaks a codon.` : '')));
        refreshView(); runQC(); commit();
      } catch (e) {
        showMessage(msgs, 'error', esc(e.message));
      }
      btn.disabled = false; btn.innerHTML = old;
      I18N.apply(btn);
    }, 30);
  }

  /* ================================================================
     trimming
     ================================================================ */
  function trimParams() {
    const p = current();
    const n = p ? p.taxa.length : 0;
    return {
      method: el('trimMethod').value,
      gap: +el('trimGap').value, sim: +el('trimSim').value,
      b1: +el('gb1').value, b2: +el('gb2').value, b3: +el('gb3').value, b4: +el('gb4').value, b5: el('gb5').value,
      codon: el('trimCodon').checked, n,
    };
  }
  function updateTrim() {
    const p = current();
    const host = el('trimStats');
    if (!host) return;
    if (!p || !p.aligned || p.type === 'morph') { host.innerHTML = ''; if (B2.view) B2.view.setMask(null); return; }
    const o = trimParams();
    const st = Trim.columnStats(p.seqs);
    let mask;
    if (o.method === 'gblocks') mask = Trim.gblocks(p.seqs, { b1: o.b1, b2: o.b2, b3: o.b3, b4: o.b4, b5: o.b5, stats: st }).keep;
    else if (o.method === 'gaps') mask = Trim.byGaps(p.seqs, o.gap, st).keep;
    else if (o.method === 'sim') mask = Trim.bySimilarity(p.seqs, o.sim, st).keep;
    else mask = Trim.intersect(Trim.byGaps(p.seqs, o.gap, st).keep, Trim.bySimilarity(p.seqs, o.sim, st).keep);
    if (o.codon) mask = Trim.codonAware(mask, p.frame || 1);
    B2.trim = mask;
    if (B2.view) B2.view.setMask(mask);
    const kept = Trim.count(mask), L = mask.length;
    const blocks = Trim.blocks(mask);
    statTiles(host, [
      [T('Columnas conservadas', 'Columns kept'), kept, `${fmtPct(L ? kept / L : 0)} ${T('del alineamiento', 'of the alignment')}`, kept / L > 0.5 ? 'ok' : 'warn'],
      [T('Columnas eliminadas', 'Columns removed'), L - kept, fmtPct(L ? (L - kept) / L : 0)],
      [T('Bloques', 'Blocks'), blocks.length, blocks.length ? `${T('el mayor', 'largest')}: ${Math.max(...blocks.map(b => b[1] - b[0] + 1))}` : ''],
      [T('Informativos que quedan', 'Informative kept'), (() => {
        let k = 0;
        st.forEach((s, c) => { if (mask[c] && s.states > 1 && Object.values(s.counts).filter(v => v >= 2).length >= 2) k++; });
        return k;
      })(), T('sitios', 'sites')],
    ]);
  }
  function applyTrim() {
    const p = current();
    if (!p || !B2.trim) return;
    const kept = Trim.count(B2.trim);
    if (!kept) { showMessage(el('alnMessages'), 'error', L2('El recorte dejaría el alineamiento vacío.', 'That trimming would leave the alignment empty.')); return; }
    p.seqs = Trim.apply(p.seqs, B2.trim);
    p.trimmed = (p.trimmed || 0) + 1;
    B2.trim = null;
    refreshView(); runQC(); commit();
    showMessage(el('alnMessages'), 'success', L2(
      `Recorte aplicado: quedan ${kept} columnas.`, `Trimming applied: ${kept} columns remain.`));
  }

  /* ================================================================
     quality control
     ================================================================ */
  function runQC() {
    const p = current();
    const host = el('qcList');
    if (!host) return;
    host.innerHTML = '';
    if (!p) return;
    const coding = el('qcCoding') && el('qcCoding').checked;
    p.coding = !!coding;
    const items = QC.check(p.taxa, p.seqs, p.type, { coding, frame: p.frame });
    if (B2.parts.length > 1) {
      const m = QC.missingTaxa(B2.parts);
      if (m.rows.length) items.unshift({
        level: 'warn', rows: m.rows.slice(0, 6).map(r => `${r.taxon} → ${T('falta en', 'missing from')} ${r.missing.join(', ')}`),
        es: `${m.rows.length} de ${m.all.length} taxones no están en todas las particiones. Al concatenar, sus huecos serán datos faltantes.`,
        en: `${m.rows.length} of ${m.all.length} taxa are not in every partition. When concatenating, their gaps become missing data.`,
      });
    }
    items.forEach(it => {
      const div = mk('div', { class: 'qc-item ' + it.level });
      div.innerHTML = `<div class="qc-ic">${it.level === 'bad' ? '!' : it.level === 'warn' ? '?' : 'i'}</div>
        <div><div>${L2(it.es, it.en)}</div>${it.rows && it.rows.length ? `<div class="qc-rows">${esc(it.rows.slice(0, 8).join(' · '))}${it.rows.length > 8 ? ' …' : ''}</div>` : ''}</div>`;
      host.appendChild(div);
    });
    I18N.apply(host);
  }

  /* ================================================================
     concatenation, codon split, partition schemes
     ================================================================ */
  function concatenate() {
    const usable = B2.parts.filter(p => !p.derived && p.aligned);
    const msgs = el('partMessages');
    clearMessages(msgs);
    if (usable.length < 2) { showMessage(msgs, 'warning', L2('Hacen falta al menos dos particiones alineadas.', 'At least two aligned partitions are needed.')); return; }
    const types = new Set(usable.map(p => p.type));
    if (types.size > 1) showMessage(msgs, 'warning', L2(
      'Estás concatenando particiones de distinto tipo; el Bloque 3 les pondrá modelos distintos.',
      'You are concatenating partitions of different types; Block 3 will give them different models.'));
    const all = [];
    usable.forEach(p => p.taxa.forEach(t => { if (all.indexOf(t) < 0) all.push(t); }));
    all.sort();
    const seqs = all.map(() => '');
    const charsets = [];
    let at = 1;
    usable.forEach(p => {
      const L = p.seqs[0].length;
      all.forEach((t, i) => {
        const k = p.taxa.indexOf(t);
        seqs[i] += k >= 0 ? p.seqs[k] : '?'.repeat(L);
      });
      charsets.push({ name: p.name, from: at, to: at + L - 1, genome: p.genome, type: p.type });
      at += L;
    });
    B2.charsets = charsets;
    const i = addPartition({ name: T('concatenado', 'concatenated'), taxa: all, seqs, type: usable[0].type, source: T('de %n particiones', 'from %n partitions').replace('%n', usable.length) });
    B2.parts[i].derived = true;
    /* a concatenation of several genomes belongs to none of them */
    const genomes = new Set(usable.map(p => p.genome));
    B2.parts[i].genome = genomes.size === 1 ? [...genomes][0] : 'other';
    B2.active = i;
    renderParts(); refreshView(); runQC(); commit();
    showMessage(msgs, 'success', L2(
      `Matriz concatenada: ${all.length} taxones × ${at - 1} posiciones, en ${charsets.length} particiones.`,
      `Concatenated matrix: ${all.length} taxa × ${at - 1} positions, in ${charsets.length} partitions.`));
  }

  function splitCodon() {
    const p = current();
    const msgs = el('partMessages');
    clearMessages(msgs);
    if (!p || !p.aligned || p.type !== 'dna') { showMessage(msgs, 'warning', L2('Se necesita una partición de ADN alineada.', 'An aligned DNA partition is needed.')); return; }
    const frame = p.frame || 1;
    for (let pos = 0; pos < 3; pos++) {
      const seqs = p.seqs.map(s => {
        let o = '';
        for (let c = (frame - 1) + pos; c < s.length; c += 3) o += s[c];
        return o;
      });
      const i = addPartition({ name: `${p.name}_pos${pos + 1}`, taxa: p.taxa.slice(), seqs, type: 'dna', source: p.name });
      B2.parts[i].genome = p.genome;
      B2.parts[i].derived = true;
    }
    renderParts(); commit();
    showMessage(msgs, 'success', L2(
      `«${esc(p.name)}» se separó en sus tres posiciones de codón. La tercera suele evolucionar mucho más rápido: dale su propio modelo en el Bloque 3.`,
      `"${esc(p.name)}" was split into its three codon positions. The third usually evolves much faster: give it its own model in Block 3.`));
  }

  function importPartitions(text) {
    const p = current();
    const msgs = el('partMessages');
    clearMessages(msgs);
    if (!p || !p.aligned) { showMessage(msgs, 'warning', L2('Primero abre una matriz alineada.', 'Open an aligned matrix first.')); return; }
    let parts = SeqIO.parsePartitionText(text);
    if (!parts.length) {
      try { const nx = SeqIO.parse(text, 'x.nex'); parts = nx.charsets || []; } catch (e) { /* not a NEXUS */ }
    }
    if (!parts.length) { showMessage(msgs, 'error', L2('No se reconoció ningún esquema de particiones.', 'No partition scheme was recognised.')); return; }
    const L = p.seqs[0].length;
    const bad = parts.filter(c => c.to > L);
    if (bad.length) { showMessage(msgs, 'error', L2(
      `El esquema llega hasta la posición ${Math.max(...parts.map(c => c.to))} y la matriz tiene ${L}.`,
      `The scheme reaches position ${Math.max(...parts.map(c => c.to))} and the matrix has ${L}.`)); return; }
    parts.forEach(cs => {
      const i = addPartition({
        name: cs.name, taxa: p.taxa.slice(),
        seqs: p.seqs.map(s => s.slice(cs.from - 1, cs.to)),
        type: p.type, source: p.name,
      });
      B2.parts[i].derived = true;
    });
    B2.charsets = parts;
    renderParts(); commit();
    showMessage(msgs, 'success', L2(`Se crearon ${parts.length} particiones.`, `${parts.length} partitions were created.`));
  }

  /* ================================================================
     export and hand-off
     ================================================================ */
  function exportAs(fmt) {
    const p = current();
    if (!p) return;
    const base = slug(p.name);
    if (fmt === 'fasta') download(SeqIO.writeFasta(p.taxa, p.seqs), base + '.fasta');
    else if (fmt === 'phylip') download(SeqIO.writePhylip(p.taxa, p.seqs, {}), base + '.phy');
    else if (fmt === 'phylip-int') download(SeqIO.writePhylip(p.taxa, p.seqs, { interleaved: true }), base + '_int.phy');
    else if (fmt === 'nexus') download(SeqIO.writeNexus(p.taxa, p.seqs, { type: p.type, charsets: B2.charsets }), base + '.nex');
    else if (fmt === 'parts') {
      const cs = B2.charsets || [{ name: p.name, from: 1, to: p.seqs[0].length, type: p.type }];
      download(SeqIO.writeRaxmlPartitions(cs.map(c => Object.assign({ type: p.type }, c))), base + '_particiones.txt');
    } else if (fmt === 'csv') {
      const L = p.seqs[0].length;
      const header = ['taxon'].concat(Array.from({ length: L }, (_, i) => String(i + 1)));
      download(matrixToCSV(header, p.taxa.map((t, i) => [t].concat(p.seqs[i].split('')))), base + '.csv', 'text/csv');
    }
  }

  /* what the later blocks read */
  function commit() {
    const usable = B2.parts.filter(p => p.aligned);
    if (!usable.length) { state.data = null; enableStep(3, false); return; }
    const all = [];
    usable.forEach(p => p.taxa.forEach(t => { if (all.indexOf(t) < 0) all.push(t); }));
    state.data = {
      taxa: all,
      parts: usable.map(p => ({
        name: p.name, genome: p.genome, type: p.type, coding: !!p.coding, frame: p.frame || 1,
        taxa: p.taxa.slice(), seqs: p.seqs.slice(), length: p.seqs[0].length,
      })),
      charsets: B2.charsets,
      active: Math.max(0, usable.indexOf(current())),
    };
    state.fileName = (current() || {}).source || null;
    /* Block 3 needs at least one aligned partition with variable sites */
    enableStep(3, true);
  }

  /* ================================================================
     wiring
     ================================================================ */
  function wire() {
    const dz = el('dropZone');
    if (!dz) return;
    ['dragenter', 'dragover'].forEach(ev => dz.addEventListener(ev, e => { e.preventDefault(); dz.classList.add('over'); }));
    ['dragleave', 'drop'].forEach(ev => dz.addEventListener(ev, e => { e.preventDefault(); dz.classList.remove('over'); }));
    dz.addEventListener('drop', e => { if (e.dataTransfer && e.dataTransfer.files.length) loadFiles(e.dataTransfer.files); });
    el('pickFiles').addEventListener('click', () => el('fileInput').click());
    el('fileInput').addEventListener('change', e => { if (e.target.files.length) loadFiles(e.target.files); e.target.value = ''; });
    el('pasteBtn').addEventListener('click', () => {
      const a = el('pasteArea');
      a.style.display = a.style.display === 'none' ? '' : 'none';
      if (a.style.display === '') el('pasteText').focus();
    });
    el('pasteCancel').addEventListener('click', () => { el('pasteArea').style.display = 'none'; });
    el('pasteLoad').addEventListener('click', () => {
      const txt = el('pasteText').value.trim();
      const msgs = el('loadMessages');
      clearMessages(msgs);
      if (!txt) return;
      try {
        loadText(txt, T('pegado', 'pasted'));
        if (B2.active < 0) B2.active = 0;
        el('pasteArea').style.display = 'none'; el('pasteText').value = '';
        renderParts(); refreshView(); runQC(); commit();
      } catch (e) { showMessage(msgs, 'error', esc(e.message)); }
    });
    els('[data-example]').forEach(b => b.addEventListener('click', () => loadExample(b.dataset.example)));

    el('alignBtn').addEventListener('click', runAlign);
    el('alnMode').addEventListener('change', () => { if (B2.view) B2.view.setMode(el('alnMode').value); });
    el('alnMarks').addEventListener('change', applyMarks);
    el('alnFrameChk').addEventListener('change', () => {
      const p = current();
      if (B2.view) B2.view.setFrame(el('alnFrameChk').checked && p ? (p.frame || 1) : null);
    });
    el('zoomIn').addEventListener('click', () => B2.view && B2.view.zoom(2));
    el('zoomOut').addEventListener('click', () => B2.view && B2.view.zoom(-2));
    el('gotoCol').addEventListener('change', () => { if (B2.view) B2.view.goTo(+el('gotoCol').value - 1); });

    const ed = {
      edUndo: () => B2.view.undo(), edRedo: () => B2.view.redoLast(),
      edGapOne: () => B2.view.insertGap(false), edGapAll: () => B2.view.insertGap(true),
      edDelGap: () => B2.view.deleteGap(), edDelCol: () => B2.view.deleteColumn(),
      edSqueeze: () => B2.view.squeeze(), edSort: () => B2.view.sortByName(),
      edUp: () => B2.view.moveRow(-1), edDown: () => B2.view.moveRow(1), edDelRow: () => B2.view.removeRow(),
    };
    Object.keys(ed).forEach(id => {
      const n = el(id);
      if (n) n.addEventListener('click', () => { if (!B2.view) return; if (!ed[id]()) flash(n); });
    });

    ['trimMethod', 'gb5'].forEach(id => el(id).addEventListener('change', () => { syncTrimUI(); updateTrim(); }));
    ['trimGap', 'trimSim', 'gb1', 'gb2', 'gb3', 'gb4'].forEach(id => {
      const n = el(id), out = el(id + 'Val');
      n.addEventListener('input', () => {
        if (out) out.textContent = (id === 'trimGap' || id === 'trimSim') ? (+n.value).toFixed(2) : n.value;
        updateTrim();
      });
    });
    el('trimCodon').addEventListener('change', updateTrim);
    el('trimApply').addEventListener('click', applyTrim);
    el('trimClear').addEventListener('click', () => { B2.trim = null; if (B2.view) B2.view.setMask(null); });

    el('qcRun').addEventListener('click', runQC);
    el('qcCoding').addEventListener('change', () => {
      const p = current();
      if (p) p.coding = el('qcCoding').checked;
      if (p && p.coding && p.type === 'dna') {
        const bf = Align.bestFrame(p.seqs);
        p.frame = bf.frame;
        if (el('alnFrameChk').checked && B2.view) B2.view.setFrame(p.frame);
      }
      runQC();
    });
    el('qcCleanNames').addEventListener('click', () => {
      const p = current();
      if (!p) return;
      p.taxa = p.taxa.map(t => t.trim().replace(/\s+/g, '_').replace(/[^\w.\-]/g, ''));
      refreshView(); runQC(); renderParts(); commit();
    });
    el('qcDropEmpty').addEventListener('click', () => { if (B2.view) { B2.view.squeeze(); runQC(); } });

    el('concatBtn').addEventListener('click', concatenate);
    el('splitCodonBtn').addEventListener('click', splitCodon);
    el('importPartsBtn').addEventListener('click', () => el('partsInput').click());
    el('partsInput').addEventListener('change', e => {
      const f = e.target.files[0];
      if (!f) return;
      const r = new FileReader();
      r.onload = () => importPartitions(String(r.result));
      r.readAsText(f);
      e.target.value = '';
    });

    els('[data-export]').forEach(b => b.addEventListener('click', () => exportAs(b.dataset.export)));
    el('toBlock3').addEventListener('click', () => {
      const b = document.querySelector('.step-btn[data-step="3"]');
      if (b && !b.disabled) goStep(3);
      else showMessage(el('alnMessages'), 'info', L2(
        'El Bloque 3 (modelos de sustitución) se construye en la etapa siguiente. Tus datos quedan cargados y listos.',
        'Block 3 (substitution models) is built in the next stage. Your data stay loaded and ready.'));
    });

    document.addEventListener('langchange', () => { renderParts(); runQC(); updateSummary(); updateTrim(); });
    document.addEventListener('stepchange', e => { if (e.detail.step === 2 && B2.view) setTimeout(() => B2.view.redraw(), 30); });
  }

  function flash(node) {
    node.animate([{ opacity: 1 }, { opacity: 0.3 }, { opacity: 1 }], { duration: 260 });
  }

  function syncTrimUI() {
    const m = el('trimMethod').value;
    el('gblocksParams').style.display = m === 'gblocks' ? '' : 'none';
    el('trimGapLab').style.display = (m === 'gaps' || m === 'both') ? '' : 'none';
    el('trimSimLab').style.display = (m === 'sim' || m === 'both') ? '' : 'none';
  }

  /* the Gblocks thresholds are counts of sequences, so they follow the data */
  function syncGblocksRanges() {
    const p = current();
    const n = p ? p.taxa.length : 10;
    const b1 = el('gb1'), b2 = el('gb2');
    b1.max = n; b2.max = n;
    if (+b1.value > n || +b1.value < Math.floor(n / 2) + 1) { b1.value = Math.floor(n / 2) + 1; el('gb1Val').textContent = b1.value; }
    if (+b2.value > n || +b2.value < +b1.value) { b2.value = Math.ceil(n * 0.85); el('gb2Val').textContent = b2.value; }
  }
  const origSelect = selectPart;
  selectPart = function (i) { origSelect(i); syncGblocksRanges(); updateTrim(); };

  function init() {
    if (!el('panel-2')) return;
    wire();
    syncTrimUI();
    renderParts();
  }
  document.addEventListener('DOMContentLoaded', init);

  Object.assign(B2, { addPartition, renderParts, selectPart, loadText, loadExample, runAlign, concatenate, splitCodon, exportAs, commit, current, guessGenome });
})();
