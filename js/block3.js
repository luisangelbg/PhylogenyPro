/* PhylogenyPro — molecular phylogenetics in the browser.
   Copyright (C) 2026  Luis Ángel Barrera-Guzmán

   This program is free software: you can redistribute it and/or modify it under
   the terms of the GNU General Public License as published by the Free Software
   Foundation, either version 3 of the License, or (at your option) any later
   version. It is distributed in the hope that it will be useful, but WITHOUT ANY
   WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
   A PARTICULAR PURPOSE. See the GNU General Public License, in the file LICENSE
   at the root of this program, or <https://www.gnu.org/licenses/>. */

/* PhylogenyPro — Block 3: substitution models.

   Reads the partitions of Block 2, checks the assumptions every model makes,
   fits the candidates and hands the winner to the blocks that follow. The
   fitting runs one model per turn of the event loop, so the page never freezes
   and the run can be cancelled; the Workers arrive with Block 5, where the
   bootstrap makes them worth the trouble. */

(function () {

  const B3 = { part: -1, data: null, tree: null, run: null, result: null, sat: null, comp: null };
  window.B3 = B3;

  function parts() { return (state.data && state.data.parts) || []; }
  function current() { return parts()[B3.part] || null; }

  /* ================================================================
     the data card
     ================================================================ */
  function renderPartSelect() {
    const sel = el('m3Part');
    if (!sel) return;
    const ps = parts();
    sel.innerHTML = ps.map((p, i) => `<option value="${i}">${esc(p.name)} · ${p.taxa.length} × ${p.length}</option>`).join('');
    if (B3.part < 0 || B3.part >= ps.length) B3.part = ps.length ? 0 : -1;
    sel.value = String(Math.max(0, B3.part));
  }

  function prepare() {
    const p = current();
    if (!p) return null;
    const A = Like.compress(p.seqs, p.type);
    B3.data = A;
    return A;
  }

  function renderData() {
    const p = current();
    if (!p) return;
    const A = prepare();
    const st = Align.stats(p.seqs, p.type);
    const gc = p.type === 'dna' ? (() => {
      let g = 0, t = 0;
      p.seqs.forEach(s => { for (let i = 0; i < s.length; i++) { const c = s[i]; if (c === 'G' || c === 'C') { g++; t++; } else if (c === 'A' || c === 'T') t++; } });
      return t ? g / t : 0;
    })() : null;
    statTiles('m3Tiles', [
      [T('Taxones', 'Taxa'), p.taxa.length, esc(p.name)],
      [T('Sitios', 'Sites'), st.length, `${A.nPat} ${T('patrones distintos', 'distinct patterns')}`],
      [T('Variables', 'Variable'), st.variable, fmtPct(st.variable / st.length)],
      [T('Informativos', 'Informative'), st.informative, T('para la parsimonia', 'for parsimony')],
      [T('Constantes', 'Constant'), st.constant, fmtPct(st.constant / st.length)],
      p.type === 'dna' ? [T('Contenido G+C', 'G+C content'), fmtPct(gc), ''] : [T('Estados', 'States'), A.nStates, T('en la matriz', 'in the matrix')],
      [T('Huecos y faltantes', 'Gaps and missing'), fmtPct(st.gapFraction), ''],
      [T('Genoma', 'Genome'), T({ cp: 'cloroplasto', mt: 'mitocondria', nr: 'nuclear rib.', nu: 'nuclear', other: 'otro' }[p.genome] || p.genome,
        { cp: 'chloroplast', mt: 'mitochondrion', nr: 'nuclear rib.', nu: 'nuclear', other: 'other' }[p.genome] || p.genome), ''],
    ]);

    /* composition table */
    if (p.type !== 'morph') {
      const comp = Dist.compositionTest(p.seqs, p.type);
      B3.comp = comp;
      const cols = [{ key: 'taxon', label: T('Taxón', 'Taxon') }].concat(
        comp.alphabet.slice(0, 20).map((a, k) => ({ key: 'f' + k, label: a, num: true, fmt: v => fmtPct(v, 1) })));
      const rows = comp.rows.map((r, i) => {
        const o = { taxon: p.taxa[i] };
        comp.alphabet.forEach((a, k) => { o['f' + k] = r.freqs[k]; });
        return o;
      });
      const all = { taxon: T('todas', 'all'), _class: 'total' };
      comp.alphabet.forEach((a, k) => { all['f' + k] = comp.grand[k]; });
      rows.push(all);
      buildTable('m3CompTable', cols, rows, { limit: 40 });
      el('m3CompFig').innerHTML = Plots3.composition(comp, p.taxa);
    } else {
      el('m3CompTable').innerHTML = '';
      el('m3CompFig').innerHTML = '';
      B3.comp = null;
    }

    /* saturation */
    const msg = el('m3AssumMsg');
    clearMessages(msg);
    if (p.type === 'dna') {
      const sat = Dist.saturation(p.seqs);
      B3.sat = sat;
      el('m3SatFig').innerHTML = Plots3.saturation(sat);
      const ratio = sat.slopeTv > 0 ? sat.slopeTs / sat.slopeTv : null;
      if (sat.undefinedPairs > 0) showMessage(msg, 'warning', L2(
        `En ${sat.undefinedPairs} pares la corrección de distancia no tiene solución: esas secuencias ya no se parecen más que dos al azar. El marcador está saturado para esa comparación.`,
        `In ${sat.undefinedPairs} pairs the distance correction has no solution: those sequences are no more similar than two random ones. The marker is saturated for that comparison.`));
      if (ratio != null && ratio < 1.05) showMessage(msg, 'warning', L2(
        `Las transiciones ya no superan a las transversiones (razón ${ratio.toFixed(2)}): señal de saturación. Considera usar solo transversiones, o un marcador más lento.`,
        `Transitions no longer outnumber transversions (ratio ${ratio.toFixed(2)}): a sign of saturation. Consider using transversions only, or a slower marker.`));
      else if (ratio != null) showMessage(msg, 'success', L2(
        `Las transiciones se acumulan ${ratio.toFixed(2)} veces más rápido que las transversiones y la relación sigue subiendo: el marcador conserva señal. La corrección tiene que recuperar un ${fmtPct(sat.meanHidden)} de cambios ocultos.`,
        `Transitions accumulate ${ratio.toFixed(2)} times faster than transversions and the relation still rises: the marker keeps its signal. The correction has to recover ${fmtPct(sat.meanHidden)} of hidden changes.`));
    } else {
      B3.sat = null;
      el('m3SatFig').innerHTML = '';
    }
    if (B3.comp) {
      const c = B3.comp;
      showMessage(msg, c.p < 0.05 ? 'error' : 'info', L2(
        `Homogeneidad de la composición: χ² = ${fmtFixed(c.chi2, 2)}, ${c.df} grados de libertad, ${pEq(c.p)}. ` +
          (c.p < 0.05 ? 'Se rechaza: las secuencias no tienen la misma composición, y los modelos habituales suponen que sí. Considera la distancia LogDet, o un modelo no homogéneo, y desconfía de los clados que agrupen secuencias por su composición.'
                      : 'No se rechaza: las secuencias comparten composición, que es lo que suponen todos los modelos de esta página.'),
        `Homogeneity of composition: χ² = ${fmtFixed(c.chi2, 2)}, ${c.df} degrees of freedom, ${pEq(c.p)}. ` +
          (c.p < 0.05 ? 'Rejected: the sequences do not share the same composition, and the standard models assume they do. Consider the LogDet distance, or a non-homogeneous model, and be suspicious of clades that group sequences by composition.'
                      : 'Not rejected: the sequences share their composition, which is what every model on this page assumes.')));
    }
  }

  /* ================================================================
     model selection
     ================================================================ */
  function startingTree(p, A) {
    const kind = el('m3Tree').value;
    const model = p.type === 'aa' ? 'poisson' : p.type === 'morph' ? 'hamming' : 'jc';
    const D = Dist.matrix(p.seqs, model, { type: p.type }).D;
    if (kind === 'nj') return Tree.nj(D, p.taxa);
    if (kind === 'upgma') return Tree.upgma(D, p.taxa);
    return Tree.bionj(D, p.taxa);
  }

  function run() {
    const p = current();
    if (!p) return;
    const A = B3.data || prepare();
    const msg = el('m3RunMsg');
    clearMessages(msg);
    const set = el('m3Set').value, ncat = +el('m3Ncat').value;
    const list = Models.buildList(p.type, { set, ncat });
    const btn = el('m3Run'), cancel = el('m3Cancel'), prog = el('m3Progress');
    btn.disabled = true; cancel.style.display = '';
    const t0 = performance.now();

    /* the base tree, optimised once under a middling model — the strategy of
       ModelFinder: every candidate then sees exactly the same tree */
    const bar = window.LABG ? LABG.progressBar(prog, { label: T('Selección de modelo', 'Model selection') }) : null;
    if (bar) bar.update(null, T('preparando el árbol de partida…', 'preparing the starting tree…'));
    else prog.textContent = T('preparando el árbol de partida…', 'preparing the starting tree…');
    setTimeout(() => {
      let base;
      try {
        const t = startingTree(p, A);
        const spec = p.type === 'dna' ? { type: 'dna', model: 'HKY', freqs: Array.from(A.freqs), alpha: 0.5, ncat }
          : p.type === 'aa' ? { type: 'aa', model: 'LG', alpha: 0.5, ncat }
            : { type: 'morph', model: 'Mk', alpha: 0.5, ncat };
        base = Like.fit(t, A, spec, { passes: 6 });
        B3.tree = base.tree;
      } catch (e) {
        showMessage(msg, 'error', esc(e.message));
        if (bar) bar.fail(T('no se pudo ajustar', 'could not fit'));
        btn.disabled = false; cancel.style.display = 'none';
        return;
      }
      const tFit = performance.now();
      const runner = Models.selectModel(B3.tree, A, list, {
        fast: true, fitOpts: { passes: 4 },
        onProgress: (done, total, last) => {
          /* after two models the pace is known, so the wait can be stated
             instead of guessed at */
          const per = (performance.now() - tFit) / done;
          const left = per * (total - done) / 1000;
          const eta = done >= 2 && left > 2
            ? ` · ${T('faltan', 'about')} ${left < 60 ? Math.ceil(left) + ' s' : Math.ceil(left / 60) + ' min'}`
            : '';
          const txt = `${done} / ${total}${eta} · ${last && last.name ? last.name : ''}`;
          if (bar) bar.update(done / total, txt);
          else prog.textContent = txt;
        },
      });
      B3.run = runner;
      runner.promise.then(res => {
        btn.disabled = false; cancel.style.display = 'none';
        if (res.cancelled || !bar) prog.textContent = '';
        if (res.cancelled) { showMessage(msg, 'info', L2('Cancelado.', 'Cancelled.')); return; }
        if (bar) bar.done(T(`${list.length} modelos ajustados`, `${list.length} models fitted`));
        B3.result = res;
        const secs = (performance.now() - t0) / 1000;
        showMessage(msg, 'success', L2(
          `${list.length} modelos ajustados en ${secs.toFixed(1)} s sobre el mismo árbol.`,
          `${list.length} models fitted in ${secs.toFixed(1)} s on the same tree.`));
        refit();
        renderResults();
        commit();
      });
    }, 30);
  }

  /* the winner gets a proper fit, with every branch optimised */
  function refit() {
    const crit = el('m3Crit').value;
    const best = B3.result['best' + crit];
    const A = B3.data;
    const spec = Object.assign({}, best.spec);
    const full = Like.fit(B3.tree, A, spec, { passes: 10 });
    full.name = Models.name(full.spec);
    B3.best = full;
    B3.bestFast = best;
  }

  function renderResults() {
    const res = B3.result, crit = el('m3Crit').value, p = current();
    if (!res) return;
    el('m3Results').style.display = '';
    el('m3Next').style.display = '';
    const best = B3.best, fast = B3.bestFast;
    const second = res.fits.filter(f => f.name !== fast.name).sort((a, b) => a[crit] - b[crit])[0];
    statTiles('m3BestTiles', [
      [T('Modelo elegido', 'Chosen model'), best.name, `${T('por', 'by')} ${crit}`, 'ok'],
      ['lnL', fmtLnL(best.lnL), `${best.k} ${T('parámetros', 'parameters')}`],
      [crit, fmtFixed(best[crit], 2), second ? `${T('el siguiente está a', 'the next one is')} ${fmtFixed(second[crit] - fast[crit], 2)}` : ''],
      [T('Peso de Akaike', 'Akaike weight'), fmtPct(fast['w' + crit] || 0), T('probabilidad relativa', 'relative probability')],
      best.spec.alpha ? [T('Forma de Γ', 'Γ shape') + ' α', fmtFixed(best.spec.alpha, 3), best.spec.alpha < 0.5 ? T('tasas muy desiguales', 'very uneven rates') : T('tasas parejas', 'even rates')] : null,
      best.spec.pInv ? [T('Sitios invariables', 'Invariable sites'), fmtPct(best.spec.pInv), ''] : null,
      best.spec.type === 'dna' && best.spec.rates ? [T('Razón ts/tv', 'ts/tv ratio') + ' κ', fmtFixed(best.spec.rates[1] / (best.spec.rates[0] || 1), 2), ''] : null,
    ].filter(Boolean));

    /* what it means, in words */
    const bits = [];
    bits.push(L2(
      `El modelo que mejor equilibra ajuste y número de parámetros es <b>${best.name}</b>.`,
      `The model that best balances fit against the number of parameters is <b>${best.name}</b>.`));
    if (best.spec.alpha) bits.push(L2(
      `La forma α = ${fmtFixed(best.spec.alpha, 3)} dice que las tasas entre sitios son ${best.spec.alpha < 0.5 ? 'muy desiguales: unos pocos sitios acumulan casi todos los cambios' : 'moderadamente desiguales'}. Ignorar esa variación acorta las ramas largas y, con ellas, las fechas.`,
      `The shape α = ${fmtFixed(best.spec.alpha, 3)} says that rates among sites are ${best.spec.alpha < 0.5 ? 'very uneven: a few sites take almost all the changes' : 'moderately uneven'}. Ignoring that variation shortens long branches and, with them, the dates.`));
    else bits.push(L2(
      'Ningún modelo con variación de tasas entre sitios (+G) ganó, algo poco común: revisa que el alineamiento no sea demasiado corto.',
      'No model with rate variation among sites (+G) won, which is unusual: check that the alignment is not too short.'));
    if (best.spec.pInv) bits.push(L2(
      `Además, un ${fmtPct(best.spec.pInv)} de los sitios se estimó invariable. +I y +G explican el mismo fenómeno de dos maneras y se estorban entre sí: si los dos están, interpreta α con cuidado.`,
      `On top of that, ${fmtPct(best.spec.pInv)} of the sites were estimated invariable. +I and +G explain the same phenomenon in two ways and interfere with each other: if both are in, read α with care.`));
    const wsum = res.fits.filter(f => f['w' + crit] > 0.05).length;
    if (wsum > 1) bits.push(L2(
      `Hay ${wsum} modelos con un peso apreciable: la elección no es rotunda, y en la práctica todos darán árboles casi iguales.`,
      `There are ${wsum} models with an appreciable weight: the choice is not clear-cut, and in practice they will all give nearly the same tree.`));
    el('m3Verdict').innerHTML = `<b>${L2('Qué significa', 'What it means')}</b> ` + bits.join(' ');

    el('m3BarsFig').innerHTML = Plots3.modelBars(res.fits.slice().sort((a, b) => a[crit] - b[crit]), crit, { top: 12 });
    el('m3RatesFig').innerHTML = Plots3.rateCategories(best.model);
    el('m3QFig').innerHTML = Plots3.rateMatrix(best.model, B3.data.alphabet);

    buildTable('m3Table', [
      { key: 'name', label: T('Modelo', 'Model') },
      { key: 'k', label: T('parámetros', 'parameters'), num: true },
      { key: 'lnL', label: 'lnL', num: true, fmt: v => fmtFixed(v, 3) },
      { key: 'AIC', label: 'AIC', num: true, fmt: v => fmtFixed(v, 2) },
      { key: 'AICc', label: 'AICc', num: true, fmt: v => fmtFixed(v, 2) },
      { key: 'BIC', label: 'BIC', num: true, fmt: v => fmtFixed(v, 2) },
      { key: 'd' + crit, label: 'Δ' + crit, num: true, fmt: v => fmtFixed(v, 2) },
      { key: 'w' + crit, label: T('peso', 'weight'), num: true, fmt: v => fmtPct(v, 1) },
    ], res.fits.slice().sort((a, b) => a[crit] - b[crit]).map((f, i) => Object.assign({ _class: i === 0 ? 'total' : '' }, f)), { limit: 60 });

    if (res.lrt && res.lrt.steps.length) {
      buildTable('m3LrtTable', [
        { key: 'step', label: T('Prueba', 'Test'), get: r => `${r.from} → ${r.to}` },
        { key: 'what', label: T('Qué añade', 'What it adds'), get: r => T(r.es, r.en) },
        { key: 'lr', label: '2ΔlnL', num: true, fmt: v => fmtFixed(v, 2) },
        { key: 'df', label: 'df', num: true },
        { key: 'p', label: 'p', num: true, fmt: v => fmtP(v) },
        { key: 'res', label: T('Decisión', 'Decision'), get: r => r.accepted ? T('se acepta el más complejo', 'the more complex one is accepted') : T('se queda el más simple', 'the simpler one stays') },
      ], res.lrt.steps, {});
      const lrtNote = mk('p', { class: 'hint' }, L2(
        `La cadena de pruebas termina en <b>${res.lrt.model}</b>. Si no coincide con el elegido por ${crit}, no es un error: son dos criterios distintos, y el de información no exige que los modelos estén anidados.`,
        `The chain of tests ends at <b>${res.lrt.model}</b>. If that is not what ${crit} chose, it is not a mistake: they are two different criteria, and the information one does not require the models to be nested.`));
      el('m3LrtTable').appendChild(lrtNote);
    }
  }

  /* ================================================================
     partition scheme
     ================================================================ */
  function mergeScheme() {
    const ps = parts().filter(p => p.aligned !== false);
    const msg = el('m3PartMsg');
    clearMessages(msg);
    if (ps.length < 2) { showMessage(msg, 'info', L2('Hace falta más de una partición.', 'More than one partition is needed.')); return; }
    const crit = el('m3Crit').value;
    const btn = el('m3Merge');
    btn.disabled = true;
    const w = window.LABG ? LABG.work({ title: T('Buscando el esquema de particiones', 'Searching for the partition scheme'), delay: 300 }) : null;
    phyAfterPaint(() => {
      try {
        /* every partition needs the same taxa, in the same order */
        const taxa = ps[0].taxa;
        const ok = ps.every(p => p.taxa.length === taxa.length && p.taxa.every((t, i) => t === taxa[i]));
        if (!ok) { showMessage(msg, 'warning', L2('Las particiones no tienen los mismos taxones en el mismo orden; concaténalas en el Bloque 2 antes de buscar el esquema.', 'The partitions do not have the same taxa in the same order; concatenate them in Block 2 first.')); btn.disabled = false; return; }
        const prepared = ps.map(p => {
          const A = Like.compress(p.seqs, p.type);
          const spec = p.type === 'dna' ? { type: 'dna', model: 'GTR', alpha: 0.5, ncat: 4 } : { type: p.type, model: p.type === 'aa' ? 'LG' : 'Mk', alpha: 0.5, ncat: 4 };
          return { name: p.name, A, seqs: p.seqs, type: p.type, fit: Like.fitFast(B3.tree || startingTree(p, A), A, spec, { passes: 4 }) };
        });
        const res = Models.mergeScheme(prepared, B3.tree || startingTree(ps[0], prepared[0].A), { criterion: crit });
        buildTable('m3PartTable', [
          { key: 'name', label: T('Grupo', 'Group') },
          { key: 'n', label: T('particiones', 'partitions'), num: true, get: r => r.ids.length },
          { key: 'lnL', label: 'lnL', num: true, get: r => fmtFixed(r.fit.lnL, 2) },
          { key: 'crit', label: crit, num: true, get: r => fmtFixed(r.fit[crit], 2) },
        ], res.groups, {});
        showMessage(msg, 'success', L2(
          `El esquema que mejor ${crit} da tiene <b>${res.groups.length}</b> grupo(s) de ${ps.length} particiones: ${res.groups.map(g => g.name).join(' | ')}.` +
            (res.groups.length < ps.length ? ' Fusionar particiones que evolucionan parecido gasta menos parámetros y suele dar ramas mejor estimadas.' : ' Ninguna fusión mejoró el criterio: cada partición merece su propio modelo.'),
          `The scheme with the best ${crit} has <b>${res.groups.length}</b> group(s) out of ${ps.length} partitions: ${res.groups.map(g => g.name).join(' | ')}.` +
            (res.groups.length < ps.length ? ' Merging partitions that evolve alike spends fewer parameters and usually gives better branch estimates.' : ' No merge improved the criterion: each partition deserves its own model.')));
        B3.scheme = res;
        commit();
        if (w) w.done();
      } catch (e) {
        showMessage(msg, 'error', esc(e.message));
      }
      btn.disabled = false;
    }, w);
  }

  /* ================================================================
     hand-off
     ================================================================ */
  function commit() {
    if (!B3.best) return;
    state.models = {
      partition: (current() || {}).name,
      best: {
        name: B3.best.name, spec: B3.best.spec, lnL: B3.best.lnL, k: B3.best.k,
        AIC: B3.best.AIC, AICc: B3.best.AICc, BIC: B3.best.BIC,
      },
      criterion: el('m3Crit').value,
      tree: B3.best.tree, lens: B3.best.lens,
      all: B3.result ? B3.result.fits.map(f => ({ name: f.name, lnL: f.lnL, k: f.k, AIC: f.AIC, AICc: f.AICc, BIC: f.BIC })) : [],
      lrt: B3.result ? B3.result.lrt : null,
      composition: B3.comp ? { chi2: B3.comp.chi2, df: B3.comp.df, p: B3.comp.p } : null,
      saturation: B3.sat ? { slopeTs: B3.sat.slopeTs, slopeTv: B3.sat.slopeTv, meanHidden: B3.sat.meanHidden, undefinedPairs: B3.sat.undefinedPairs } : null,
      scheme: B3.scheme ? { groups: B3.scheme.groups.map(g => ({ name: g.name, ids: g.ids })), criterion: B3.scheme.criterion } : null,
    };
    enableStep(4, true);
  }

  function exportModel() {
    if (!B3.best) return;
    const b = B3.best, p = current();
    const lines = [];
    lines.push(`# PhylogenyPro — ${T('modelo elegido', 'chosen model')}`);
    lines.push(`# ${T('partición', 'partition')}: ${p.name} (${p.taxa.length} × ${p.length}, ${p.genome})`);
    lines.push(`model = ${b.name}`);
    lines.push(`lnL = ${b.lnL.toFixed(6)}`);
    lines.push(`parameters = ${b.k}`);
    lines.push(`AIC = ${b.AIC.toFixed(4)}   AICc = ${b.AICc.toFixed(4)}   BIC = ${b.BIC.toFixed(4)}`);
    if (b.spec.rates) lines.push(`rates (AC AG AT CG CT GT) = ${b.spec.rates.map(v => (+v).toFixed(6)).join(' ')}`);
    lines.push(`frequencies = ${Array.from(b.model.pi).map(v => v.toFixed(6)).join(' ')}`);
    if (b.spec.alpha) lines.push(`alpha = ${b.spec.alpha.toFixed(6)}   categories = ${b.spec.ncat}`);
    if (b.spec.pInv) lines.push(`invariable = ${b.spec.pInv.toFixed(6)}`);
    lines.push('');
    lines.push(`# ${T('árbol con las longitudes de este modelo', 'tree with the branch lengths of this model')}`);
    lines.push(Tree.writeNewick(b.tree, { labels: p.taxa }));
    download(lines.join('\n'), slug(p.name) + '_modelo.txt');
  }

  /* ================================================================
     wiring
     ================================================================ */
  function refresh() {
    const has = parts().length > 0;
    el('m3NoData').style.display = has ? 'none' : '';
    ['m3Data', 'm3Assum', 'm3Select'].forEach(id => { el(id).style.display = has ? '' : 'none'; });
    el('m3Parts').style.display = parts().length > 1 ? '' : 'none';
    if (!has) return;
    renderPartSelect();
    renderData();
  }

  function init() {
    if (!el('panel-3')) return;
    el('m3Part').addEventListener('change', () => {
      B3.part = +el('m3Part').value;
      B3.result = null; B3.best = null;
      el('m3Results').style.display = 'none';
      el('m3Next').style.display = 'none';
      renderData();
    });
    el('m3Run').addEventListener('click', run);
    el('m3Cancel').addEventListener('click', () => { if (B3.run) B3.run.cancel(); });
    el('m3Crit').addEventListener('change', () => { if (B3.result) { refit(); renderResults(); commit(); } });
    el('m3Merge').addEventListener('click', mergeScheme);
    el('m3Export').addEventListener('click', exportModel);
    el('m3ToBlock4').addEventListener('click', () => {
      const b = document.querySelector('.step-btn[data-step="4"]');
      if (b && !b.disabled) goStep(4);
      else showMessage(el('m3RunMsg'), 'info', L2(
        'El Bloque 4 (parsimonia y distancias) llega en la etapa siguiente. El modelo elegido queda guardado.',
        'Block 4 (parsimony and distances) arrives in the next stage. The chosen model is saved.'));
    });
    document.addEventListener('stepchange', e => { if (e.detail.step === 3) refresh(); });
    document.addEventListener('langchange', () => { if (parts().length) { renderData(); if (B3.result) renderResults(); } });
    document.addEventListener('themechange', () => { if (parts().length) { renderData(); if (B3.result) renderResults(); } });
  }
  document.addEventListener('DOMContentLoaded', init);

  Object.assign(B3, { refresh, run, renderData, commit, current, mergeScheme });
})();
