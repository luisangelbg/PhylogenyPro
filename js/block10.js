/* PhylogenyPro — molecular phylogenetics in the browser.
   Copyright (C) 2026  Luis Ángel Barrera-Guzmán

   This program is free software: you can redistribute it and/or modify it under
   the terms of the GNU General Public License as published by the Free Software
   Foundation, either version 3 of the License, or (at your option) any later
   version. It is distributed in the hope that it will be useful, but WITHOUT ANY
   WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
   A PARTICULAR PURPOSE. See the GNU General Public License, in the file LICENSE
   at the root of this program, or <https://www.gnu.org/licenses/>. */

/* PhylogenyPro — Block 10: historical biogeography.

   The engine is js/biogeo.js. What this file adds is the reading, and this
   block needs it more than any other, because its output is a picture that
   looks like an answer:

     · choosing the areas is the analysis. Nothing here can check that choice,
       and a reconstruction on badly chosen areas is confident and wrong;
     · the +J models win almost always, and Ree & Sanmartín (2018) showed the
       comparison is not a fair one — j buys at cladogenesis what d and e pay
       for along branches. The block fits them, reports them, and says so;
     · a node reconstructed at 45 % is drawn as 45 %, because a single coloured
       square there would be a claim nobody made;
     · d and e are per unit of time. On a tree in substitutions they are per
       substitution, which is not a rate anybody can compare with a paper. */

(function () {

  const B10 = {
    source: 'dated', tree: null, parsed: null, matched: null,
    states: null, multipliers: null, fits: null, anc: null, bsm: null,
    maxAreas: null, ancModel: null,
  };
  window.B10 = B10;

  /* ================================================================
     the tree
     ================================================================ */
  function sources() {
    const out = [];
    if (state.dated && state.dated.tree) {
      out.push(['dated', T('Bloque 7 · árbol fechado', 'Block 7 · dated tree'), state.dated.tree, true, state.dated.partition]);
    }
    if (state.ml && state.ml.tree) {
      out.push(['ml', T('Bloque 5 · árbol de verosimilitud', 'Block 5 · likelihood tree'), state.ml.tree, false, state.ml.partition]);
    }
    if (state.bayes && state.bayes.consensus) {
      out.push(['bayes', T('Bloque 6 · consenso bayesiano', 'Block 6 · Bayesian consensus'), state.bayes.consensus, false, state.bayes.partition]);
    }
    if (state.quick && state.quick.distance && state.quick.distance.tree) {
      out.push(['distance', T('Bloque 4 · árbol de distancias', 'Block 4 · distance tree'), state.quick.distance.tree, false, state.quick.partition]);
    }
    return out;
  }
  function chosen() { return sources().find(x => x[0] === B10.source) || sources()[0] || null; }
  function chosenTree() { const s = chosen(); return s ? Tree.clone(s[2]) : null; }
  function isDated() { const s = chosen(); return !!(s && s[3]); }
  function unit() { return (state.dated && state.dated.unit) || 'Ma'; }
  function taxaOf() {
    const s = chosen();
    const name = s ? s[4] : null;
    if (name && state.data && state.data.parts) {
      const p = state.data.parts.find(q => q.name === name);
      if (p && p.taxa) return p.taxa;
    }
    return (state.data && state.data.taxa) || [];
  }
  function treeLabels(tree) {
    const taxa = taxaOf();
    const out = [];
    Tree.tips(tree).forEach(t => { out[t.tip] = t.label || taxa[t.tip] || null; });
    return out;
  }

  /* ================================================================
     1 · the areas
     ================================================================ */
  function readAreas() {
    const msg = el('p10DataMsg');
    msg.innerHTML = '';
    const parsed = Biogeo.parseRanges(el('p10Text').value);
    if (!parsed || !parsed.taxa.length || !parsed.areas.length) {
      showMessage(msg, 'error', L2(
        'No se reconoció la tabla. Se espera el formato de LAGRANGE (una línea con el número de taxones, el de áreas y sus nombres entre paréntesis, y luego una fila por especie con ceros y unos) o una tabla separada por comas con una columna por área.',
        'The table was not recognised. Either the LAGRANGE format (a line with the number of taxa, the number of areas and their names in brackets, then one row per species of zeros and ones) or a comma-separated table with one column per area.'));
      return;
    }
    const tree = chosenTree();
    if (!tree) return;
    B10.tree = tree;
    B10.parsed = parsed;
    B10.fits = B10.anc = B10.bsm = null;

    /* line the ranges up with the tree by name */
    const labels = treeLabels(tree);
    const norm = s => String(s).replace(/[\s_]+/g, ' ').trim().toLowerCase();
    const byName = new Map(), byNorm = new Map();
    parsed.taxa.forEach((nm, i) => {
      byName.set(String(nm), i);
      if (!byNorm.has(norm(nm))) byNorm.set(norm(nm), i);
    });
    const ranges = [], missing = [];
    labels.forEach((lab, tipRow) => {
      if (lab == null) return;
      const i = byName.has(lab) ? byName.get(lab) : (byNorm.has(norm(lab)) ? byNorm.get(norm(lab)) : -1);
      if (i < 0) { missing.push(lab); ranges[tipRow] = null; return; }
      ranges[tipRow] = parsed.masks[i];
    });
    const used = new Set(labels.map((lab, i) => (ranges[i] != null ? lab : null)).filter(Boolean));
    const extra = parsed.taxa.filter(nm => !used.has(nm) && !used.has(String(nm)));
    B10.matched = { ranges, labels, missing, extra, matched: ranges.filter(m => m != null).length };

    /* a species with no area at all is not missing data, it is a mistake */
    const empty = [];
    B10.matched.ranges.forEach((m, i) => { if (m === 0) empty.push(labels[i]); });

    B10.multipliers = identity(parsed.areas.length);
    B10.maxAreas = null;
    showAreas();
    if (missing.length) {
      showMessage(msg, missing.length === labels.length ? 'error' : 'warn', L2(
        `Sin áreas en la tabla: ${missing.slice(0, 8).join(', ')}${missing.length > 8 ? '…' : ''}. Todo el bloque necesita un rango por punta, así que hay que completarlos o podar el árbol.`,
        `No areas in the table for: ${missing.slice(0, 8).join(', ')}${missing.length > 8 ? '…' : ''}. The whole block needs a range for every tip, so they have to be filled in or the tree pruned.`));
    }
    if (extra.length) {
      showMessage(msg, 'info', L2(
        `Filas que no corresponden a ninguna punta: ${extra.slice(0, 8).join(', ')}${extra.length > 8 ? '…' : ''}. Se ignoran.`,
        `Rows matching no tip: ${extra.slice(0, 8).join(', ')}${extra.length > 8 ? '…' : ''}. They are ignored.`));
    }
    if (empty.length) {
      showMessage(msg, 'error', L2(
        `Estas especies no aparecen en ninguna área: ${empty.slice(0, 8).join(', ')}. Una fila de ceros es el rango nulo, que en este modelo significa «extinta», y ninguna punta puede estarlo.`,
        `These species are in no area at all: ${empty.slice(0, 8).join(', ')}. A row of zeros is the null range, which in this model means "extinct", and no tip can be.`));
    }
    if (!isDated()) {
      showMessage(msg, 'warn', L2(
        'Este árbol no está fechado, así que <b>d y e quedan por sustitución</b>, no por unidad de tiempo. La reconstrucción de los rangos sigue siendo válida; las tasas no se pueden comparar con las de ningún artículo. El cronograma del Bloque 7 es lo que hace falta.',
        'This tree is not dated, so <b>d and e come out per substitution</b>, not per unit of time. The reconstruction of the ranges is still valid; the rates cannot be compared with any paper\'s. The chronogram of Block 7 is what that needs.'));
    }
    commit();
  }
  const identity = n => Array.from({ length: n }, () => new Float64Array(n).fill(1));

  function showAreas() {
    const p = B10.parsed, m = B10.matched;
    const nA = p.areas.length;
    const counts = p.areas.map((_, i) => m.ranges.filter(r => r != null && (r & (1 << i))).length);
    const sizes = m.ranges.filter(r => r != null && r > 0).map(r => Biogeo.popcount(r));
    statTiles('p10DataTiles', [
      [T('Puntas con rango', 'Tips with a range'), `${m.matched} / ${m.labels.length}`,
        m.missing.length ? T(`${m.missing.length} sin fila`, `${m.missing.length} with no row`) : '',
        m.missing.length ? 'bad' : 'ok'],
      [T('Áreas', 'Areas'), nA, p.areas.join(' ')],
      [T('Especies en más de un área', 'Species in more than one area'),
        sizes.filter(s => s > 1).length, T(`de ${sizes.length}`, `of ${sizes.length}`)],
      [T('Rango observado más ancho', 'Widest observed range'),
        sizes.length ? Math.max.apply(null, sizes) : 0, T('áreas', 'areas')],
    ]);
    buildTable('p10DataTable', [
      { key: 'area', label: T('área', 'area') },
      { key: 'n', label: T('especies', 'species'), num: true },
      { key: 'sole', label: T('sólo ahí', 'only there'), num: true },
      { key: 'who', label: T('quiénes (primeras)', 'which (first few)') },
    ], p.areas.map((a, i) => {
      const who = m.labels.filter((lab, k) => m.ranges[k] != null && (m.ranges[k] & (1 << i)));
      const sole = m.labels.filter((lab, k) => m.ranges[k] === (1 << i));
      return { area: a, n: counts[i], sole: sole.length, who: who.slice(0, 4).join(', ') + (who.length > 4 ? '…' : '') };
    }));
    el('p10DataInfo').textContent = T(`${m.matched} de ${m.labels.length} puntas emparejadas`,
      `${m.matched} of ${m.labels.length} tips matched`);

    /* the state-space controls */
    const sel = el('p10MaxAreas');
    const widest = sizes.length ? Math.max.apply(null, sizes) : 1;
    sel.innerHTML = '';
    for (let k = widest; k <= nA; k++) {
      const nStates = countStates(nA, k) + (el('p10Null').checked ? 1 : 0);
      sel.innerHTML += `<option value="${k}">${k} (${nStates} ${T('estados', 'states')})</option>`;
    }
    if (B10.maxAreas == null || B10.maxAreas < widest || B10.maxAreas > nA) B10.maxAreas = nA;
    sel.value = String(B10.maxAreas);
    ['p10States', 'p10Models', 'p10Export'].forEach(id => { el(id).style.display = ''; });
    el('p10Anc').style.display = 'none';
    el('p10Bsm').style.display = 'none';
    refreshStates();
    refreshMultipliers();
  }
  function countStates(n, maxK) {
    let t = 0;
    for (let k = 1; k <= maxK; k++) t += choose(n, k);
    return t;
  }
  function choose(n, k) { let r = 1; for (let i = 0; i < k; i++) r = r * (n - i) / (i + 1); return Math.round(r); }

  /* ================================================================
     2 · the state space
     ================================================================ */
  function refreshStates() {
    if (!B10.parsed) return;
    const nA = B10.parsed.areas.length;
    B10.maxAreas = +el('p10MaxAreas').value || nA;
    B10.states = Biogeo.makeStates(nA, { maxAreas: B10.maxAreas, includeNull: el('p10Null').checked });
    const n = B10.states.n;
    el('p10StateInfo').textContent = T(`${n} estados`, `${n} states`);

    /* a range that the cap would forbid is a range the data say exists */
    const forbidden = B10.matched.ranges.filter(r => r != null && r > 0 && !B10.states.index.has(r));
    const v = el('p10StateVerdict');
    v.style.display = '';
    const es = [], en = [];
    es.push(`Con ${nA} áreas y un tope de ${B10.maxAreas} hay <b>${n} estados</b>.`);
    en.push(`With ${nA} areas and a cap of ${B10.maxAreas} there are <b>${n} states</b>.`);
    /* the cost grows with the cube of the state count, and saying so before the
       user presses the button is cheaper than saying it afterwards */
    const cost = Math.pow(n / 16, 3);
    if (n > 120) {
      es.push(`Cada rama cuesta la exponencial de una matriz de ${n} × ${n}, unas <b>${cost.toFixed(0)} veces</b> lo que cuesta un problema de cuatro áreas. Ajustar los seis modelos va a tardar minutos, y el panel no se congela pero tampoco corre.`);
      en.push(`Every branch costs the exponential of an ${n} × ${n} matrix, about <b>${cost.toFixed(0)} times</b> what a four-area problem costs. Fitting the six models will take minutes; the panel will not freeze, but it will not be quick either.`);
    } else if (n > 40) {
      es.push(`Eso es unas ${cost.toFixed(1)} veces el coste de un problema de cuatro áreas: del orden de decenas de segundos por modelo.`);
      en.push(`That is about ${cost.toFixed(1)} times a four-area problem: tens of seconds per model.`);
    }
    if (forbidden.length) {
      es.push(`<b>Cuidado:</b> ${forbidden.length} especie(s) ocupan un rango más ancho que el tope, así que el modelo no tiene un estado para lo que se observa. Sube el tope.`);
      en.push(`<b>Careful:</b> ${forbidden.length} species occupy a range wider than the cap, so the model has no state for what is observed. Raise the cap.`);
    }
    if (!el('p10Null').checked) {
      es.push('Sin el rango nulo, un linaje de una sola área no puede extinguirse localmente y e queda parcialmente sin sentido; BioGeoBEARS lo incluye por omisión y la mayoría de los artículos también.');
      en.push('Without the null range, a single-area lineage cannot go locally extinct and e is left partly meaningless; BioGeoBEARS includes it by default and so do most papers.');
    }
    v.innerHTML = `<b>${L2('Cómo se lee', 'How to read it')}</b> ${L2(es.join(' '), en.join(' '))}`;
  }

  function refreshMultipliers() {
    if (!B10.parsed) return;
    const areas = B10.parsed.areas;
    ['p10MultFrom', 'p10MultTo'].forEach((id, k) => {
      const s = el(id);
      s.innerHTML = areas.map((a, i) => `<option value="${i}">${esc(a)}</option>`).join('');
      s.value = String(k === 0 ? 0 : Math.min(1, areas.length - 1));
    });
    el('p10MultFig').innerHTML = Plots10.matrixHeat(B10.multipliers, areas, {
      title: T('de la fila a la columna', 'from the row to the column'),
      caption: T('uno es «sin preferencia»; cero prohíbe el movimiento', 'one is "no preference"; zero forbids the move'),
    });
  }

  /* ================================================================
     3 · the models
     ================================================================ */
  function runModels() {
    const msg = el('p10ModelMsg');
    msg.innerHTML = '';
    if (!B10.states || !B10.matched) return;
    if (B10.matched.missing.length) {
      showMessage(msg, 'error', L2('Faltan rangos para algunas puntas; el ajuste no puede empezar.',
        'Some tips have no range; the fit cannot start.'));
      return;
    }
    const which = el('p10Which').value;
    const list = which === 'dec' ? [['DEC', false], ['DEC', true]]
      : which === 'noj' ? [['DEC', false], ['DIVALIKE', false], ['BAYAREALIKE', false]]
        : [['DEC', false], ['DEC', true], ['DIVALIKE', false], ['DIVALIKE', true],
           ['BAYAREALIKE', false], ['BAYAREALIKE', true]];
    const opts = { multipliers: allOnes(B10.multipliers) ? null : B10.multipliers };
    const bar = window.LABG ? LABG.progressBar(el('p10ModelProgress'), { label: T('Modelos biogeográficos', 'Biogeographic models') }) : null;
    if (bar) bar.update(0, T(`ajustando ${list.length} modelos…`, `fitting ${list.length} models…`));
    else el('p10ModelProgress').innerHTML = L2(`ajustando ${list.length} modelos…`, `fitting ${list.length} models…`);
    el('p10RunModels').disabled = true;
    const t0 = performance.now();
    /* one model per turn of the event loop, so the page keeps breathing and the
       progress note is true rather than decorative */
    const prep = Biogeo.prepare(B10.tree, B10.matched.ranges, B10.states, opts);
    const rows = [];
    let i = 0;
    const step = () => {
      if (i >= list.length) { finish(); return; }
      const [m, wj] = list[i];
      const f = Biogeo.fit(B10.tree, B10.matched.ranges, B10.states,
        Object.assign({}, opts, { model: m, withJ: wj, prep }));
      rows.push(f);
      i++;
      if (bar) bar.update(i / list.length, T(`${i} de ${list.length}…`, `${i} of ${list.length}…`));
      else el('p10ModelProgress').innerHTML = L2(`${i} de ${list.length}…`, `${i} of ${list.length}…`);
      setTimeout(step, 10);
    };
    const finish = () => {
      const best = rows.reduce((a, b) => (b.AICc < a.AICc ? b : a), rows[0]);
      let sw = 0;
      rows.forEach(r => { r.dAICc = r.AICc - best.AICc; r.w = Math.exp(-0.5 * r.dAICc); sw += r.w; });
      rows.forEach(r => { r.w /= sw; });
      rows.forEach(r => {
        if (!r.withJ) return;
        const base = rows.find(q => q.model === r.model && !q.withJ);
        if (!base) return;
        const stat = 2 * (r.lnL - base.lnL);
        r.lrt = { stat, p: stat <= 0 ? 1 : 1 - Dist.pchisq(stat, 1), against: base.name };
      });
      rows.sort((a, b) => a.AICc - b.AICc);
      B10.fits = { rows, best: rows[0], prep, opts };
      if (bar) bar.done(T(
        `${list.length} modelos en ${((performance.now() - t0) / 1000).toFixed(1)} s`,
        `${list.length} models in ${((performance.now() - t0) / 1000).toFixed(1)} s`));
      else el('p10ModelProgress').innerHTML = L2(
        `${list.length} en ${((performance.now() - t0) / 1000).toFixed(1)} s`,
        `${list.length} in ${((performance.now() - t0) / 1000).toFixed(1)} s`);
      el('p10RunModels').disabled = false;
      showModels();
      commit();
    };
    setTimeout(step, 20);
  }
  function allOnes(M) {
    for (let i = 0; i < M.length; i++) for (let j = 0; j < M.length; j++) if (M[i][j] !== 1) return false;
    return true;
  }

  function showModels() {
    const { rows, best } = B10.fits;
    const per = isDated() ? T(`por ${unit()}`, `per ${unit()}`) : T('por sustitución', 'per substitution');
    buildTable('p10ModelTable', [
      { key: 'name', label: T('modelo', 'model') },
      { key: 'what', label: T('qué permite al dividirse', 'what it allows at a split') },
      { key: 'd', label: 'd', num: true },
      { key: 'e', label: 'e', num: true },
      { key: 'j', label: 'j', num: true },
      { key: 'lnL', label: 'lnL', num: true },
      { key: 'k', label: 'k', num: true },
      { key: 'aicc', label: 'AICc', num: true },
      { key: 'd2', label: 'ΔAICc', num: true },
      { key: 'w', label: T('peso', 'weight'), num: true },
    ], rows.map(r => ({
      name: r.name,
      what: r.model === 'DEC' ? T('simpatría, subconjunto y vicarianza estrecha', 'sympatry, subset and narrow vicariance')
        : r.model === 'DIVALIKE' ? T('sólo vicarianza, de cualquier corte', 'vicariance only, any cut')
          : T('nada: las dos hijas heredan todo', 'nothing: both daughters inherit everything'),
      d: fmtNum(r.d), e: fmtNum(r.e), j: r.withJ ? fmtNum(r.j) : '—',
      lnL: fmtLnL(r.lnL), k: r.k, aicc: fmtFixed(r.AICc, 3),
      d2: fmtFixed(r.dAICc, 3), w: fmtFixed(r.w, 3),
      _class: r.name === best.name ? 'row-best' : '',
    })), { caption: T(`d y e están ${per}`, `d and e are ${per}`) });

    el('p10ModelPane').style.display = '';
    el('p10ModelFig').innerHTML = Plots10.modelBars(rows, {});
    registerFigure('models', el('p10ModelFig').innerHTML);

    const v = el('p10ModelVerdict');
    v.style.display = '';
    const es = [], en = [];
    const second = rows[1];
    es.push(`Gana <b>${best.name}</b> con un peso de Akaike de ${fmtFixed(best.w, 3)}${second ? `, ${fmtFixed(second.dAICc, 2)} unidades por delante de ${second.name}` : ''}.`);
    en.push(`<b>${best.name}</b> wins with an Akaike weight of ${fmtFixed(best.w, 3)}${second ? `, ${fmtFixed(second.dAICc, 2)} units ahead of ${second.name}` : ''}.`);
    if (best.e < 1e-6) {
      es.push('La tasa de extinción se estimó en cero. No quiere decir que no haya habido extinción: quiere decir que en este árbol no hay ninguna huella de ella, y con árboles pequeños eso es lo habitual.');
      en.push('The extinction rate was estimated at zero. That does not mean there was no extinction: it means this tree shows no trace of it, and with small trees that is the usual outcome.');
    }
    if (isDated()) {
      const mean = 1 / Math.max(1e-12, best.d);
      es.push(`Con d = ${fmtNum(best.d)}, un linaje en un área tarda en promedio ${fmtNum(mean)} ${unit()} en llegar a cada área nueva.`);
      en.push(`With d = ${fmtNum(best.d)}, a lineage in one area takes on average ${fmtNum(mean)} ${unit()} to reach each new area.`);
    }
    v.innerHTML = `<b>${L2('Cómo se lee', 'How to read it')}</b> ${L2(es.join(' '), en.join(' '))}`;

    /* the warning that has to be there whenever a +J model wins */
    const jRows = rows.filter(r => r.withJ);
    const w = el('p10JWarning');
    if (jRows.length && best.withJ) {
      w.style.display = '';
      const lrt = best.lrt;
      w.innerHTML = `<b>${L2('Sobre el modelo +J', 'About the +J model')}</b> ` + L2(
        `${best.name} gana por ${fmtFixed(2 * (best.lnL - rows.find(r => r.model === best.model && !r.withJ).lnL), 2)} unidades de verosimilitud${lrt ? ` (${pEq(lrt.p)} en la razón de verosimilitud frente a ${lrt.against}, con j en el borde del espacio, así que esa p es conservadora)` : ''}. <b>Esto es lo esperable y no es prueba de nada por sí solo.</b> Ree y Sanmartín (2018) mostraron que la comparación no es limpia: j compra probabilidad en la cladogénesis —un evento instantáneo, sin coste a lo largo de las ramas— mientras que d y e la pagan integrando por todo el tiempo del árbol. En los datos que ellos examinaron, +J ganaba incluso cuando los datos se habían simulado sin ningún salto. Lo que sí se puede decir es <b>cuánto</b> gana y sobre qué modelo, que es lo que está en la tabla.`,
        `${best.name} wins by ${fmtFixed(2 * (best.lnL - rows.find(r => r.model === best.model && !r.withJ).lnL), 2)} likelihood units${lrt ? ` (${pEq(lrt.p)} on the likelihood ratio against ${lrt.against}, with j on the boundary of its space, so that p is conservative)` : ''}. <b>This is the expected outcome and on its own proves nothing.</b> Ree and Sanmartín (2018) showed the comparison is not a clean one: j buys probability at cladogenesis — an instantaneous event, costing nothing along the branches — while d and e pay for theirs by integrating over the whole time of the tree. In the cases they examined, +J won even on data simulated with no jumps at all. What can be said is <b>how much</b> it wins by and over which model, and that is what the table holds.`);
    } else w.style.display = 'none';

    /* the model menu of the next card */
    const sel = el('p10AncModel');
    sel.innerHTML = rows.map(r => `<option value="${esc(r.name)}">${esc(r.name)}</option>`).join('');
    sel.value = best.name;
    B10.ancModel = best.name;
    el('p10Anc').style.display = '';
    el('p10Export').style.display = '';
  }

  /* ================================================================
     4 · ancestral ranges
     ================================================================ */
  function runAnc() {
    const msg = el('p10AncMsg');
    msg.innerHTML = '';
    if (!B10.fits) return;
    const name = el('p10AncModel').value;
    const fit = B10.fits.rows.find(r => r.name === name) || B10.fits.best;
    B10.ancModel = name;
    const anc = Biogeo.ancestral(fit, B10.fits.opts);
    const best = Biogeo.bestRanges(anc, B10.parsed.areas, fit.prep);
    B10.anc = { fit, anc, best };

    const doubtful = best.filter(b => b.probability < 0.8).length;
    const root = best.find(b => b.node === 0);
    statTiles('p10AncTiles', [
      [T('Nodos reconstruidos', 'Nodes reconstructed'), best.length, ''],
      [T('Rango en la raíz', 'Range at the root'), root ? root.range : '—',
        root ? fmtFixed(root.probability, 3) : '', root && root.probability > 0.8 ? 'ok' : 'bad'],
      [T('Nodos por debajo de 0.8', 'Nodes below 0.8'), doubtful,
        T(`de ${best.length}`, `of ${best.length}`), doubtful > best.length / 2 ? 'bad' : ''],
      [T('Modelo', 'Model'), fit.name, `lnL ${fmtLnL(fit.lnL)}`],
    ]);

    const probs = new Map(anc.nodes.map(nd => [nd.node, nd.probs]));
    el('p10AncPane').style.display = '';
    el('p10AncFig').innerHTML = Plots10.rangeTree(B10.tree, {
      states: B10.states, areaNames: B10.parsed.areas, labels: B10.matched.labels,
      tipRanges: B10.matched.ranges, probs,
      corners: el('p10Corners').checked ? anc.corners : null,
      minProbability: Math.max(0, Math.min(0.5, +el('p10MinP').value || 0.05)),
      unit: isDated() ? unit() : T('sust./sitio', 'subst./site'),
      restLabel: T('resto', 'the rest'),
    });
    registerFigure('ranges', el('p10AncFig').innerHTML);

    buildTable('p10AncTable', [
      { key: 'clade', label: T('clado (especies)', 'clade (species)'), num: true },
      { key: 'range', label: T('rango más probable', 'most probable range') },
      { key: 'p', label: T('probabilidad', 'probability'), num: true },
      { key: 'second', label: T('el siguiente', 'the next one') },
      { key: 'p2', label: T('su probabilidad', 'its probability'), num: true },
    ], best.slice().sort((a, b) => b.nTips - a.nTips).map(b => ({
      clade: b.nTips, range: b.range, p: fmtFixed(b.probability, 3),
      second: b.second || '—', p2: fmtFixed(b.secondProbability, 3),
      _class: b.probability < 0.5 ? 'row-flag' : '',
    })), { limit: 40 });

    const v = el('p10AncVerdict');
    v.style.display = '';
    const es = [], en = [];
    if (root) {
      es.push(`El ancestro de todo el grupo sale en <b>${root.range}</b> con probabilidad ${fmtFixed(root.probability, 3)}${root.second ? `, y el siguiente candidato es ${root.second} con ${fmtFixed(root.secondProbability, 3)}` : ''}.`);
      en.push(`The ancestor of the whole group comes out in <b>${root.range}</b> with probability ${fmtFixed(root.probability, 3)}${root.second ? `, the next candidate being ${root.second} at ${fmtFixed(root.secondProbability, 3)}` : ''}.`);
    }
    es.push(doubtful
      ? `${doubtful} de ${best.length} nodos no llegan a 0.8 en ningún rango: son los pasteles repartidos de la figura y no se deben describir como si tuvieran un área.`
      : `Los ${best.length} nodos pasan de 0.8 en algún rango, lo que con este número de áreas es poco común.`);
    en.push(doubtful
      ? `${doubtful} of ${best.length} nodes reach 0.8 in no range: those are the divided pies in the figure, and should not be described as if they had an area.`
      : `All ${best.length} nodes pass 0.8 in some range, which with this many areas is uncommon.`);
    es.push('Y una advertencia que no depende de los números: <b>esta reconstrucción sólo puede poner a los linajes en las áreas que se le dieron.</b> Si el grupo estuvo en una región que no está en la tabla, el modelo repartirá esa historia entre las áreas que sí están, con toda confianza.');
    en.push('And one warning that does not depend on the numbers: <b>this reconstruction can only put lineages in the areas it was given.</b> If the group was once in a region that is not in the table, the model will share that history out among the areas that are, with every appearance of confidence.');
    v.innerHTML = `<b>${L2('Cómo se lee', 'How to read it')}</b> ${L2(es.join(' '), en.join(' '))}`;

    el('p10Bsm').style.display = '';
    el('p10Export').style.display = '';
    commit();
  }

  /* ================================================================
     5 · stochastic mapping
     ================================================================ */
  function runBsm() {
    const msg = el('p10BsmMsg');
    msg.innerHTML = '';
    if (!B10.anc) {
      showMessage(msg, 'error', L2('Primero hay que reconstruir los rangos arriba.', 'The ranges above have to be reconstructed first.'));
      return;
    }
    const reps = Math.max(50, Math.min(20000, +el('p10Reps').value || 500));
    el('p10BsmProgress').innerHTML = L2(`sorteando ${reps} historias…`, `drawing ${reps} histories…`);
    const w = window.LABG ? LABG.work({
      title: T('Mapeo estocástico biogeográfico', 'Biogeographical stochastic mapping'),
      message: T(`sorteando ${reps} historias…`, `drawing ${reps} histories…`), delay: 300,
    }) : null;
    phyAfterPaint(() => {
      const t0 = performance.now();
      const bsm = Biogeo.stochasticMap(B10.anc.fit,
        Object.assign({}, B10.fits.opts, { reps, seed: 41, ancestral: B10.anc.anc, maxTries: 200 }));
      const secs = (performance.now() - t0) / 1000;
      el('p10BsmProgress').innerHTML = L2(`${reps} en ${secs.toFixed(1)} s`, `${reps} in ${secs.toFixed(1)} s`);
      B10.bsm = bsm;

      const labels = {
        dispersal: T('dispersión (a lo largo de una rama)', 'dispersal (along a branch)'),
        extinction: T('extinción local', 'local extinction'),
        sympatry: T('simpatría (sin cambio al dividirse)', 'sympatry (no change at the split)'),
        subset: T('simpatría de subconjunto', 'subset sympatry'),
        vicariance: T('vicarianza', 'vicariance'),
        jump: T('salto fundador', 'founder jump'),
      };
      const nNodes = B10.anc.best.length;
      statTiles('p10BsmTiles', [
        [T('Historias sorteadas', 'Histories drawn'), reps,
          bsm.failures ? T(`${bsm.failures} rechazos agotados`, `${bsm.failures} rejections exhausted`) : '',
          bsm.failures > reps * 0.05 ? 'bad' : ''],
        [T('Dispersiones', 'Dispersals'), fmtFixed(bsm.mean.dispersal, 2), `± ${fmtFixed(bsm.sd.dispersal, 2)}`],
        [T('Extinciones locales', 'Local extinctions'), fmtFixed(bsm.mean.extinction, 2), `± ${fmtFixed(bsm.sd.extinction, 2)}`],
        [T('Saltos fundadores', 'Founder jumps'), fmtFixed(bsm.mean.jump, 2),
          B10.anc.fit.withJ ? `± ${fmtFixed(bsm.sd.jump, 2)}` : T('el modelo no los tiene', 'the model has none')],
      ]);
      el('p10BsmPane').style.display = '';
      el('p10BsmFig').innerHTML = Plots10.eventBars(bsm.mean, bsm.sd, labels, {});
      registerFigure('events', el('p10BsmFig').innerHTML);

      const v = el('p10BsmVerdict');
      v.style.display = '';
      const clado = bsm.mean.sympatry + bsm.mean.subset + bsm.mean.vicariance + bsm.mean.jump;
      const es = [], en = [];
      es.push(`Cada historia necesita ${fmtFixed(bsm.mean.dispersal, 1)} dispersiones y ${fmtFixed(bsm.mean.extinction, 1)} extinciones locales a lo largo de las ramas, y reparte los ${Math.round(clado)} nodos internos entre los cuatro tipos de división.`);
      en.push(`Each history needs ${fmtFixed(bsm.mean.dispersal, 1)} dispersals and ${fmtFixed(bsm.mean.extinction, 1)} local extinctions along the branches, and shares the ${Math.round(clado)} internal nodes out among the four kinds of split.`);
      const dom = ['sympatry', 'subset', 'vicariance', 'jump'].reduce((a, b) => (bsm.mean[b] > bsm.mean[a] ? b : a), 'sympatry');
      es.push(`La división más frecuente es <b>${labels[dom].toLowerCase()}</b>, con ${fmtFixed(bsm.mean[dom], 1)} nodos por historia.`);
      en.push(`The commonest split is <b>${labels[dom].toLowerCase()}</b>, at ${fmtFixed(bsm.mean[dom], 1)} nodes per history.`);
      if (bsm.sd.dispersal > bsm.mean.dispersal * 0.5) {
        es.push('La desviación de las dispersiones es grande comparada con su media: los datos no fijan cuántas hubo, sólo el orden de magnitud.');
        en.push("The spread of the dispersals is large compared with their mean: the data do not pin down how many there were, only the order of magnitude.");
      }
      if (bsm.failures > reps * 0.02) {
        es.push(`En ${bsm.failures} ocasiones el sorteo de la historia de una rama agotó sus intentos y esa rama se contó sin eventos; eso subestima un poco las cuentas de dispersión y extinción. Suele pasar con ramas largas y muchas áreas.`);
        en.push(`On ${bsm.failures} occasions the draw of a branch's history ran out of attempts and that branch was counted with no events; that underestimates the dispersal and extinction counts a little. It happens with long branches and many areas.`);
      }
      v.innerHTML = `<b>${L2('Cómo se lee', 'How to read it')}</b> ${L2(es.join(' '), en.join(' '))}`;
      commit();
      if (w) w.done();
    }, w);
  }

  /* ================================================================
     keeping the result, and taking it away
     ================================================================ */
  function registerFigure(name, svgText) {
    state.figures = state.figures || {};
    state.figures['b10_' + name] = { block: 10, name, svg: svgText };
  }
  function commit() {
    if (!B10.parsed) return;
    state.biogeo = {
      source: B10.source,
      areas: B10.parsed.areas,
      matched: B10.matched ? { matched: B10.matched.matched, nTips: B10.matched.labels.length, missing: B10.matched.missing } : null,
      maxAreas: B10.maxAreas, includeNull: B10.states ? B10.states.includeNull : true,
      nStates: B10.states ? B10.states.n : null,
      multipliers: allOnes(B10.multipliers) ? null : B10.multipliers.map(r => Array.from(r)),
      models: B10.fits ? B10.fits.rows.map(r => ({
        name: r.name, model: r.model, withJ: r.withJ, d: r.d, e: r.e, j: r.j,
        lnL: r.lnL, k: r.k, AICc: r.AICc, dAICc: r.dAICc, w: r.w,
        lrt: r.lrt ? { stat: r.lrt.stat, p: r.lrt.p, against: r.lrt.against } : null,
      })) : null,
      best: B10.fits ? B10.fits.best.name : null,
      ancestral: B10.anc ? { model: B10.anc.fit.name, nodes: B10.anc.best } : null,
      stochastic: B10.bsm ? { reps: B10.bsm.reps, mean: B10.bsm.mean, sd: B10.bsm.sd, failures: B10.bsm.failures } : null,
    };
  }

  function exportAs(kind) {
    const stamp = new Date().toISOString().slice(0, 10);
    const csv = rows => rows.map(r => r.join(',')).join('\n');
    if (kind === 'models' && B10.fits) {
      const rows = [['model', 'd', 'e', 'j', 'lnL', 'k', 'AICc', 'dAICc', 'weight', 'LRT_stat', 'LRT_p', 'LRT_against']];
      B10.fits.rows.forEach(r => rows.push([r.name, r.d, r.e, r.withJ ? r.j : '', r.lnL, r.k, r.AICc, r.dAICc, r.w,
        r.lrt ? r.lrt.stat : '', r.lrt ? r.lrt.p : '', r.lrt ? r.lrt.against : '']));
      download(csv(rows), `phylogenypro_biogeo_modelos_${stamp}.csv`, 'text/csv');
    } else if (kind === 'anc' && B10.anc) {
      const rows = [['model', 'node', 'clade_size', 'range', 'probability', 'next_range', 'next_probability']];
      B10.anc.best.forEach(b => rows.push([B10.anc.fit.name, b.node, b.nTips, b.range, b.probability,
        b.second || '', b.secondProbability]));
      download(csv(rows), `phylogenypro_biogeo_ancestrales_${stamp}.csv`, 'text/csv');
    } else if (kind === 'events' && B10.bsm) {
      const rows = [['event', 'mean_per_history', 'sd']];
      Object.keys(B10.bsm.mean).forEach(k => rows.push([k, B10.bsm.mean[k], B10.bsm.sd[k]]));
      rows.push(['replicates', B10.bsm.reps, '']);
      rows.push(['rejections_exhausted', B10.bsm.failures, '']);
      download(csv(rows), `phylogenypro_biogeo_eventos_${stamp}.csv`, 'text/csv');
    } else if (kind === 'lagrange' && B10.parsed) {
      const labels = B10.matched.labels.filter((l, i) => B10.matched.ranges[i] != null);
      const masks = B10.matched.ranges.filter(m => m != null);
      download(Biogeo.writeLagrange(B10.parsed.areas, labels, masks),
        `phylogenypro_areas_${stamp}.data`, 'text/plain');
    } else if (kind === 'svg') {
      const figs = Object.keys(state.figures || {}).filter(k => k.indexOf('b10_') === 0);
      figs.forEach(k => download('<?xml version="1.0" encoding="UTF-8"?>\n' + state.figures[k].svg,
        `phylogenypro_${k}_${stamp}.svg`, 'image/svg+xml'));
    }
  }

  /* ================================================================
     wiring
     ================================================================ */
  function refresh() {
    const src = sources();
    const has = src.length > 0;
    el('p10NoTree').style.display = has ? 'none' : '';
    el('p10Data').style.display = has ? '' : 'none';
    if (!has) {
      ['p10States', 'p10Models', 'p10Anc', 'p10Bsm', 'p10Export'].forEach(id => { el(id).style.display = 'none'; });
      return;
    }
    const sel = el('p10Source');
    if (!src.some(s => s[0] === B10.source)) B10.source = src[0][0];
    sel.innerHTML = src.map(s => `<option value="${s[0]}">${esc(s[1])}</option>`).join('');
    sel.value = B10.source;
    if (B10.parsed) readAreas();
  }

  function init() {
    if (!el('panel-10')) return;
    el('p10Source').addEventListener('change', () => {
      B10.source = el('p10Source').value;
      if (B10.parsed) readAreas();
    });
    el('p10Upload').addEventListener('click', () => el('p10File').click());
    el('p10File').addEventListener('change', e => {
      const f = e.target.files && e.target.files[0];
      if (!f) return;
      const r = new FileReader();
      r.onload = () => { el('p10Text').value = r.result; readAreas(); };
      r.readAsText(f);
    });
    el('p10Example').addEventListener('click', () => {
      el('p10Text').value = (window.EXAMPLES && EXAMPLES.areas ? EXAMPLES.areas : '').trim();
      readAreas();
    });
    el('p10Read').addEventListener('click', readAreas);
    el('p10MaxAreas').addEventListener('change', refreshStates);
    el('p10Null').addEventListener('change', () => { showAreas(); });
    el('p10MultSet').addEventListener('click', () => {
      const i = +el('p10MultFrom').value, j = +el('p10MultTo').value;
      const v = Math.max(0, Math.min(10, +el('p10MultVal').value));
      if (i === j) return;
      B10.multipliers[i][j] = v;
      B10.multipliers[j][i] = v;
      refreshMultipliers();
    });
    el('p10MultReset').addEventListener('click', () => {
      B10.multipliers = identity(B10.parsed.areas.length);
      refreshMultipliers();
    });
    el('p10RunModels').addEventListener('click', runModels);
    el('p10RunAnc').addEventListener('click', runAnc);
    el('p10AncModel').addEventListener('change', () => { if (B10.anc) runAnc(); });
    el('p10Corners').addEventListener('change', () => { if (B10.anc) runAnc(); });
    el('p10RunBsm').addEventListener('click', runBsm);
    els('[data-p10export]').forEach(b => b.addEventListener('click', () => exportAs(b.dataset.p10export)));
    el('p10ToBlock11').addEventListener('click', () => {
      const b = document.querySelector('.step-btn[data-step="11"]');
      if (b && !b.disabled) goStep(11);
      else showMessage(el('p10DataMsg'), 'info', L2(
        'El Bloque 11 (comparar árboles) llega en la etapa siguiente.',
        'Block 11 (comparing trees) arrives in the next stage.'));
    });
    document.addEventListener('stepchange', e => { if (e.detail.step === 10) refresh(); });
    document.addEventListener('langchange', () => {
      if (!sources().length) return;
      refresh();
      if (B10.fits) showModels();
      if (B10.anc) runAnc();
    });
    document.addEventListener('themechange', () => {
      if (B10.parsed) refreshMultipliers();
      if (B10.anc) runAnc();
    });
  }
  document.addEventListener('DOMContentLoaded', init);

  Object.assign(B10, { refresh, readAreas, refreshStates, runModels, runAnc, runBsm, commit });
})();
