/* PhylogenyPro — molecular phylogenetics in the browser.
   Copyright (C) 2026  Luis Ángel Barrera-Guzmán

   This program is free software: you can redistribute it and/or modify it under
   the terms of the GNU General Public License as published by the Free Software
   Foundation, either version 3 of the License, or (at your option) any later
   version. It is distributed in the hope that it will be useful, but WITHOUT ANY
   WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
   A PARTICULAR PURPOSE. See the GNU General Public License, in the file LICENSE
   at the root of this program, or <https://www.gnu.org/licenses/>. */

/* PhylogenyPro — Block 7: the molecular clock and dating.

   The engine is js/clock.js; this file is the order in which the user is asked
   for things, and it is deliberate: root first, then find out whether the clock
   holds, then write down the calibrations, and only then date. Doing it in any
   other order hides a decision.

   The block refuses to invent a calibration. Without one it will still date the
   tree — relative to a root of age 1 — and says plainly that the numbers are
   proportions, not years. */

(function () {

  const B7 = {
    part: -1, A: null, source: 'ml', tree: null, rooted: null,
    dates: null, regression: null, test: null,
    calibs: [], result: null, method: null, cv: null,
    cancelled: false,
  };
  window.B7 = B7;

  function parts() { return (state.data && state.data.parts) || []; }
  function current() { return parts()[B7.part] || null; }
  function prepare() {
    const p = current();
    if (!p) return null;
    if (!B7.A || B7.A.__part !== p.name) {
      B7.A = Like.compress(p.seqs, p.type);
      B7.A.__part = p.name;
    }
    return B7.A;
  }

  /* the trees this block can start from */
  function sources() {
    const out = [];
    if (state.ml && state.ml.tree) out.push(['ml', T('Bloque 5 · máxima verosimilitud', 'Block 5 · maximum likelihood'), state.ml.tree, state.ml.partition]);
    if (state.bayes && state.bayes.consensus) out.push(['bayes', T('Bloque 6 · consenso bayesiano', 'Block 6 · Bayesian consensus'), state.bayes.consensus, state.bayes.partition]);
    if (state.quick && state.quick.distance && state.quick.distance.tree) out.push(['distance', T('Bloque 4 · distancias', 'Block 4 · distance'), state.quick.distance.tree, state.quick.partition]);
    if (state.quick && state.quick.parsimony && state.quick.parsimony.tree) out.push(['parsimony', T('Bloque 4 · parsimonia', 'Block 4 · parsimony'), state.quick.parsimony.tree, state.quick.partition]);
    return out;
  }
  function chosenTree() {
    const s = sources().find(x => x[0] === B7.source);
    return s ? Tree.clone(s[2]) : null;
  }
  function taxaOf() {
    const p = current();
    return p ? p.taxa : (state.data ? state.data.taxa : []);
  }

  /* ================================================================
     1 · rooting
     ================================================================ */
  function parseDates() {
    const txt = (el('p7Dates').value || '').trim();
    if (!txt) return null;
    const m = new Map();
    txt.split(/\r?\n/).forEach(line => {
      const parts2 = line.split(/[,;\t]/);
      if (parts2.length < 2) return;
      const name = parts2[0].trim();
      const v = parseFloat(parts2[1]);
      if (name && isFinite(v)) m.set(name, v);
    });
    return m.size ? m : null;
  }

  function doRoot() {
    const msg = el('p7RootMsg');
    clearMessages(msg);
    const tree = chosenTree();
    if (!tree) { showMessage(msg, 'warning', L2('No hay árbol.', 'There is no tree.')); return; }
    const taxa = taxaOf();
    const mode = el('p7RootMode').value;
    let rooted = tree, note = '', noteEs = '', noteEn = '';
    B7.regression = null;
    el('p7RegPane').style.display = 'none';

    if (mode === 'outgroup') {
      const sel = [...el('p7Outgroup').selectedOptions].map(o => +o.value);
      if (!sel.length) { showMessage(msg, 'warning', L2('Elige al menos un taxón del grupo externo.', 'Choose at least one outgroup taxon.')); return; }
      const r = Tree.rootByOutgroup(Tree.unroot(tree), sel);
      rooted = r.tree;
      if (!r.monophyletic) showMessage(msg, 'warning', L2(
        'El grupo externo no es monofilético en este árbol, así que la raíz se colocó donde mejor se pudo. Eso ya es un resultado: revisa si esos taxones son de verdad el grupo externo.',
        'The outgroup is not monophyletic on this tree, so the root was placed as best it could be. That is a result in itself: check whether those taxa really are the outgroup.'));
      noteEs = `grupo externo de ${sel.length} taxones`; noteEn = `outgroup of ${sel.length} taxa`;
    } else if (mode === 'midpoint') {
      rooted = Tree.midpointRoot(Tree.unroot(tree));
      noteEs = 'punto medio'; noteEn = 'midpoint';
    } else if (mode === 'dates') {
      const dates = parseDates();
      if (!dates) { showMessage(msg, 'warning', L2('Escribe las fechas de las puntas.', 'Write the tip dates.')); return; }
      B7.dates = dates;
      const best = Clock.rootByDates(Tree.unroot(tree), dates, { labels: taxa });
      if (!best) { showMessage(msg, 'error', L2('No se pudo ajustar la regresión: ¿coinciden los nombres?', 'The regression could not be fitted: do the names match?')); return; }
      rooted = best.tree;
      B7.regression = best.fit;
      noteEs = `regresión, r² ${fmtFixed(best.fit.r2, 3)}`; noteEn = `regression, r² ${fmtFixed(best.fit.r2, 3)}`;
      el('p7RegPane').style.display = '';
      el('p7RegFig').innerHTML = Plots7.regression(best.fit, {});
      showMessage(msg, best.fit.r2 > 0.5 ? 'success' : 'warning', L2(
        `La regresión da r² = ${fmtFixed(best.fit.r2, 4)} y una tasa de ${best.fit.rate.toExponential(3)} sustituciones por sitio y unidad de tiempo. ` +
        (best.fit.direction === 'age'
          ? 'La pendiente es negativa, así que tus fechas son edades hacia atrás.'
          : 'La pendiente es positiva, así que tus fechas van hacia adelante en el tiempo.') +
        (best.fit.r2 < 0.5 ? ' Con una r² tan baja la señal temporal es floja: fechar con estos datos es arriesgado, y conviene decirlo en el artículo.' : ''),
        `The regression gives r² = ${fmtFixed(best.fit.r2, 4)} and a rate of ${best.fit.rate.toExponential(3)} substitutions per site per unit of time. ` +
        (best.fit.direction === 'age'
          ? 'The slope is negative, so your dates are ages before the present.'
          : 'The slope is positive, so your dates run forwards in time.') +
        (best.fit.r2 < 0.5 ? ' With an r² that low the temporal signal is weak: dating these data is risky, and the paper should say so.' : '')));
    } else {
      rooted = tree;
      noteEs = 'sin cambios'; noteEn = 'unchanged';
    }

    B7.rooted = rooted;
    el('p7RootInfo').textContent = T(noteEs, noteEn);
    if (mode !== 'dates') showMessage(msg, 'success', L2(
      `Árbol enraizado (${noteEs}). Ya se puede probar el reloj.`,
      `Tree rooted (${noteEn}). The clock can now be tested.`));
    ['p7Test', 'p7Calib', 'p7Date'].forEach(id => { el(id).style.display = ''; });
    fillCalTaxa();
    renderCalTable();
  }

  /* ================================================================
     2 · the clock test
     ================================================================ */
  function runTest() {
    const msg = el('p7TestMsg');
    clearMessages(msg);
    if (!B7.rooted) { showMessage(msg, 'info', L2('Enraiza primero.', 'Root it first.')); return; }
    const p = current();
    if (!p) return;
    const A = prepare();
    const spec = (state.ml && state.ml.spec && state.ml.partition === p.name)
      ? Like.cloneSpec(state.ml.spec)
      : { type: p.type, model: p.type === 'aa' ? 'LG' : 'HKY', rates: [1, 2, 1, 1, 2, 1], freqs: Array.from(A.freqs), alpha: 0.5, ncat: 4 };
    const btn = el('p7RunTest'), prog = el('p7TestProgress');
    btn.disabled = true;
    prog.textContent = T('ajustando los dos árboles…', 'fitting both trees…');
    const w = window.LABG ? LABG.work({ title: T('Prueba del reloj molecular', 'Molecular clock test'), message: prog.textContent, delay: 300 }) : null;
    phyAfterPaint(() => {
      try {
        const t0 = performance.now();
        const res = Clock.clockTest(B7.rooted, A, spec, { passes: 14 });
        B7.test = res;
        statTiles('p7TestTiles', [
          [T('lnL libre', 'free lnL'), fmtLnL(res.free.lnL), T('una longitud por rama', 'one length per branch')],
          [T('lnL con reloj', 'clock lnL'), fmtLnL(res.clock.lnL), T('árbol ultramétrico', 'ultrametric tree')],
          ['2 ΔlnL', fmtFixed(res.stat, 3), `${res.df} ${T('grados de libertad', 'degrees of freedom')}`],
          ['p', fmtP(res.p), res.p < 0.05 ? T('se rechaza el reloj', 'the clock is rejected') : T('no se rechaza', 'not rejected'),
            res.p < 0.05 ? 'bad' : 'ok'],
          [T('Tiempo', 'Time'), `${((performance.now() - t0) / 1000).toFixed(1)} s`, ''],
        ]);
        const v = el('p7TestVerdict');
        v.style.display = '';
        v.innerHTML = `<b>${L2('Cómo se lee', 'How to read it')}</b> ` + (res.p < 0.05
          ? L2(
            `El reloj estricto se rechaza (${pEq(res.p)}): forzar que todas las puntas queden a la misma distancia de la raíz cuesta ${fmtFixed(res.stat / 2, 1)} unidades de verosimilitud. ` +
            'Esto es lo normal en datos reales y <b>no</b> impide fechar: lo que impide es fechar con un reloj estricto. Usa la verosimilitud penalizada o el reloj relajado de abajo, que dejan que la velocidad cambie entre ramas.',
            `The strict clock is rejected (${pEq(res.p)}): forcing every tip to the same distance from the root costs ${fmtFixed(res.stat / 2, 1)} log-likelihood units. ` +
            'That is the normal outcome with real data and does <b>not</b> stop you dating: it stops you dating with a strict clock. Use the penalised likelihood or the relaxed clock below, which let the rate change between branches.')
          : L2(
            `No hay motivo para rechazar el reloj estricto (${pEq(res.p)}). Eso no lo demuestra —la prueba tiene poca potencia con pocos taxones—, pero permite usar el método más sencillo y explicar por qué.`,
            `There is no reason to reject the strict clock (${pEq(res.p)}). That does not prove it — the test has little power with few taxa — but it lets you use the simplest method and say why.`));
        showMessage(msg, 'success', L2(
          `Prueba terminada. El árbol con reloj es ultramétrico por construcción y se puede usar tal cual si se le da una escala.`,
          `Test finished. The clock tree is ultrametric by construction and can be used as it is once it is given a scale.`));
        if (w) w.done();
      } catch (e) {
        showMessage(msg, 'error', esc(e.message));
        if (w) w.close();
      }
      btn.disabled = false;
      prog.textContent = '';
    }, w);
  }

  /* ================================================================
     3 · calibrations
     ================================================================ */
  function fillCalTaxa() {
    const taxa = taxaOf();
    const sel = el('p7CalTaxa');
    if (sel) sel.innerHTML = taxa.map((t, i) => `<option value="${i}">${esc(t)}</option>`).join('');
    const og = el('p7Outgroup');
    if (og && !og.options.length) og.innerHTML = taxa.map((t, i) => `<option value="${i}">${esc(t)}</option>`).join('');
  }

  /* ================================================================
     3b · el asistente de calibración
     ================================================================ */
  /* Lee lo que hay escrito en los controles y lo devuelve con la forma que
     entienden CalibHelp y Clock. No toca nada. */
  function calibFromControls() {
    const type = el('p7CalType').value;
    const min = +el('p7CalMin').value, max = +el('p7CalMax').value;
    const mean = +el('p7CalMean').value, sd = +el('p7CalSd').value;
    const c = { type };
    if (type === 'fixed') c.value = min;
    else if (type === 'uniform') { c.min = min; c.max = max; }
    else if (type === 'normal') { c.mean = mean; c.sd = sd; }
    else if (type === 'lognormal') { c.min = min; c.M = mean; c.S = sd; }
    else if (type === 'exponential') { c.min = min; c.mean = mean; }
    return c;
  }

  const calUnit = () => (B7.relative ? '' : (el('p7Unit') ? el('p7Unit').value || 'Ma' : 'Ma'));

  function fillSources() {
    const s = el('p7CalSource');
    if (!s || !window.CalibHelp) return;
    s.innerHTML = CalibHelp.SOURCES.map(x =>
      `<option value="${x.id}" data-es="${esc(x.name[0])}" data-en="${esc(x.name[1])}">${esc(T(x.name[0], x.name[1]))}</option>`).join('');
    refreshSourceHint();
  }

  function refreshSourceHint() {
    const why = el('p7CalWhy');
    if (!why || !window.CalibHelp) return;
    const s = CalibHelp.sourceOf(el('p7CalSource').value);
    why.innerHTML = L2(s.why[0], s.why[1]);
    const w = el('p7CalAge2Wrap');
    if (w) w.style.display = s.id === 'geo' ? '' : 'none';
  }

  /* El botón: propone la forma de prior que corresponde a esa procedencia y
     rellena los controles. La decisión sigue siendo del usuario; lo que se
     evita es que tenga que adivinar qué significan M y S. */
  function suggestPrior() {
    if (!window.CalibHelp) return;
    const s = el('p7CalSource').value;
    const sug = CalibHelp.suggest(s, +el('p7CalAge').value, +el('p7CalAge2').value);
    if (!sug) {
      showMessage(el('p7CalMsg'), 'warn', L2('Escribe una edad válida.', 'Enter a valid age.'));
      return;
    }
    el('p7CalType').value = sug.type;
    el('p7CalType').dispatchEvent(new Event('change'));
    if (sug.type === 'fixed') el('p7CalMin').value = sug.value;
    else if (sug.type === 'uniform') { el('p7CalMin').value = sug.min; el('p7CalMax').value = sug.max; }
    else if (sug.type === 'normal') { el('p7CalMean').value = sug.mean; el('p7CalSd').value = sug.sd; }
    else if (sug.type === 'lognormal') { el('p7CalMin').value = sug.min; el('p7CalMean').value = sug.M; el('p7CalSd').value = sug.S; }
    showMessage(el('p7CalMsg'), 'ok', L2(sug.nota[0], sug.nota[1]));
    refreshAdvice();
  }

  /* Lo que el prior dice de verdad, en años, y los avisos. Se recalcula con
     cada tecla: el usuario ve el efecto de mover un número antes de añadir
     nada. */
  function refreshAdvice() {
    const says = el('p7CalSays'), box = el('p7CalAdvice');
    if (!says || !window.CalibHelp) return;
    const c = calibFromControls();
    const d = CalibHelp.describe(c, calUnit());
    says.innerHTML = `<div class="caja-dice">${L2(d[0], d[1])}</div>`;

    box.innerHTML = '';
    const picked = el('p7CalTaxa') ? [...el('p7CalTaxa').selectedOptions].length : 0;
    const avisos = CalibHelp.advise(c, {
      source: el('p7CalSource') ? el('p7CalSource').value : 'fossil',
      others: B7.calibs, isRoot: picked === 0, cladeSize: picked, unit: calUnit(),
    });
    avisos.forEach(a => showMessage(box, a.level === 'bad' ? 'warn' : a.level === 'warn' ? 'warn' : 'info',
      L2(a.es, a.en)));
  }

  function addCalibration() {
    const msg = el('p7CalMsg');
    clearMessages(msg);
    if (!B7.rooted) { showMessage(msg, 'info', L2('Enraiza primero.', 'Root it first.')); return; }
    const picked = [...el('p7CalTaxa').selectedOptions].map(o => +o.value);
    const type = el('p7CalType').value;
    const c = { type, tips: picked.slice() };
    const min = +el('p7CalMin').value, max = +el('p7CalMax').value;
    const mean = +el('p7CalMean').value, sd = +el('p7CalSd').value;
    if (type === 'fixed') c.value = min;
    else if (type === 'uniform') { c.min = min; c.max = max; }
    else if (type === 'normal') { c.mean = mean; c.sd = sd; }
    else if (type === 'lognormal') { c.min = min; c.M = mean; c.S = sd; }
    else if (type === 'exponential') { c.min = min; c.mean = mean; }
    const F = Tree.flatten(B7.rooted);
    c.node = picked.length ? Clock.mrca(F, picked) : 0;
    if (c.node < 0) { showMessage(msg, 'warning', L2('Esos taxones no forman un clado en este árbol.', 'Those taxa do not form a clade on this tree.')); return; }
    /* the split this calibration names, so it survives a change of rooting */
    c.splitKey = c.node === 0 ? 'root' : [...Tree.splits(B7.rooted, taxaOf().length).entries()]
      .filter(([, node]) => node === F.nodes[c.node]).map(([k]) => k)[0] || null;
    B7.calibs.push(c);
    renderCalTable();
    showMessage(msg, 'success', L2(
      `Calibración añadida sobre ${picked.length ? `${picked.length} taxones` : 'la raíz'}.`,
      `Calibration added on ${picked.length ? `${picked.length} taxa` : 'the root'}.`));
  }

  function renderCalTable() {
    const taxa = taxaOf();
    const rows = B7.calibs.map((c, i) => {
      const p = Clock.prepare(c);
      const names = (c.tips || []).map(k => taxa[k]);
      return {
        i: i + 1,
        clade: !names.length ? T('la raíz', 'the root')
          : names.length <= 3 ? names.join(', ') : `${names.slice(0, 2).join(', ')} … (${names.length})`,
        type: c.type,
        summary: p.type === 'fixed' ? fmtFixed(p.value, 3)
          : p.q && p.q.q50 != null
            ? `${fmtFixed(p.q.q025, 2)} – ${fmtFixed(p.q.q50, 2)} – ${fmtFixed(p.q.q975, 2)}`
            : `≥ ${fmtFixed(p.hardMin || 0, 2)}`,
      };
    });
    buildTable('p7CalTable', [
      { key: 'i', label: '#', num: true },
      { key: 'clade', label: T('Clado', 'Clade') },
      { key: 'type', label: T('tipo', 'type') },
      { key: 'summary', label: T('2.5 % – mediana – 97.5 %', '2.5 % – median – 97.5 %') },
    ], rows, {});
  }

  /* ================================================================
     4 · dating
     ================================================================ */
  function runDate() {
    const msg = el('p7DateMsg');
    clearMessages(msg);
    if (!B7.rooted) { showMessage(msg, 'info', L2('Enraiza primero.', 'Root it first.')); return; }
    const p = current();
    const A = prepare();
    const method = el('p7Method').value;
    const btn = el('p7RunDate'), cancel = el('p7CancelDate'), prog = el('p7DateProgress');
    btn.disabled = true; B7.cancelled = false;
    let calibs = B7.calibs.slice();
    let relative = false;
    if (!calibs.length) {
      calibs = [{ node: 0, type: 'fixed', value: 1 }];
      relative = true;
      showMessage(msg, 'warning', L2(
        'No hay ninguna calibración, así que la raíz se fija en 1 y las edades salen <b>relativas</b>: proporciones de la profundidad del árbol, no años. Para obtener años hace falta un fósil, un acontecimiento geológico, una tasa conocida o puntas fechadas.',
        'There is no calibration, so the root is fixed at 1 and the ages come out <b>relative</b>: proportions of the depth of the tree, not years. Getting years needs a fossil, a geological event, a known rate or dated tips.'));
    }
    const nSites = A.nSites;
    prog.textContent = T('fechando…', 'dating…');

    if (method === 'relaxed') {
      cancel.style.display = '';
      const spec = (state.ml && state.ml.spec && state.ml.partition === p.name)
        ? Like.cloneSpec(state.ml.spec)
        : { type: p.type, model: 'HKY', rates: [1, 2, 1, 1, 2, 1], freqs: Array.from(A.freqs), alpha: 0.5, ncat: 4 };
      const gens = Math.max(2000, +el('p7Gens').value || 60000);
      const t0 = performance.now();
      const w = window.LABG ? LABG.work({ title: T('Reloj molecular relajado', 'Relaxed molecular clock'), delay: 300 }) : null;
      phyAfterPaint(() => {
        const res = Clock.relaxed(B7.rooted, A, spec, calibs, {
          generations: gens, burnin: Math.floor(gens * 0.3), sampleEvery: Math.max(1, Math.round(gens / 2000)),
          seed: 17, nSites,
          onProgress: (gen, total) => {
            prog.textContent = `${fmtNum(gen, 0)}/${fmtNum(total, 0)} · ${((performance.now() - t0) / 1000).toFixed(0)} s`;
            if (w) w.update(gen / total, prog.textContent);
          },
          cancelled: () => B7.cancelled,
        });
        B7.result = res; B7.method = 'relaxed'; B7.relative = relative;
        showRelaxed(res, relative, (performance.now() - t0) / 1000);
        btn.disabled = false; cancel.style.display = 'none'; prog.textContent = '';
        if (w) w.done();
      }, w);
      return;
    }

    const w = window.LABG ? LABG.work({
      title: method === 'pl' ? T('Fechando por verosimilitud penalizada', 'Dating by penalised likelihood')
        : T('Fechando por mínimos cuadrados', 'Dating by least squares'),
      delay: 300,
    }) : null;
    phyAfterPaint(() => {
      try {
        const t0 = performance.now();
        let res;
        if (method === 'pl') {
          const lambda = Math.max(1e-3, +el('p7Lambda').value || 1);
          res = Clock.penalised(B7.rooted, nSites, calibs, { lambda, passes: 40 });
        } else {
          res = Clock.lsd(B7.rooted, nSites, calibs, {});
        }
        B7.result = res; B7.method = method; B7.relative = relative;
        showPoint(res, method, relative, (performance.now() - t0) / 1000);
        if (w) w.done();
      } catch (e) {
        showMessage(msg, 'error', esc(e.message));
        if (w) w.close();
      }
      btn.disabled = false; prog.textContent = '';
    }, w);
  }

  function showPoint(res, method, relative, secs) {
    const unit = relative ? '' : ' ' + (el('p7Unit').value || 'Ma');
    const tiles = [
      [T('Edad de la raíz', 'Root age'), fmtFixed(res.rootAge, 3) + unit, relative ? T('relativa', 'relative') : '', 'ok'],
      [T('Método', 'Method'), method === 'pl' ? T('penalizada', 'penalised') : T('mínimos cuadrados', 'least squares'),
        method === 'pl' ? `λ = ${fmtFixed(res.lambda, 3)}` : ''],
      [T('Tiempo', 'Time'), `${secs.toFixed(2)} s`, ''],
    ];
    if (method === 'pl') {
      const rv = Clock.rateVariation(res);
      if (rv) {
        tiles.push([T('Tasa media', 'Mean rate'), res.meanRate.toExponential(3), T('sustituciones/sitio/unidad', 'substitutions/site/unit')]);
        tiles.push([T('Variación de la tasa', 'Rate variation'), fmtFixed(rv.cv, 3), T('coeficiente de variación', 'coefficient of variation'),
          rv.cv > 0.5 ? 'bad' : 'ok']);
      }
    } else {
      tiles.push([T('Tasa', 'Rate'), res.rate.toExponential(3), T('sustituciones/sitio/unidad', 'substitutions/site/unit')]);
      tiles.push(['RSS', fmtFixed(res.rss, 3), T('residuos ponderados', 'weighted residuals')]);
    }
    statTiles('p7DateTiles', tiles);
    const v = el('p7DateVerdict');
    v.style.display = '';
    v.innerHTML = `<b>${L2('Cómo se lee', 'How to read it')}</b> ` + (relative
      ? L2('Sin calibración, estas cifras son proporciones de la profundidad del árbol. Sirven para comparar nodos entre sí, no para escribir «hace X millones de años».',
        'Without a calibration these numbers are proportions of the depth of the tree. They serve to compare nodes with each other, not to write "X million years ago".')
      : L2(
        `La raíz queda en ${fmtFixed(res.rootAge, 2)}${unit}. Un método puntual como este da una cifra por nodo y ninguna incertidumbre: ` +
        'la que importa viene de las calibraciones, y solo el reloj relajado bayesiano la propaga. Si vas a publicar edades, córrelo.',
        `The root comes out at ${fmtFixed(res.rootAge, 2)}${unit}. A point method like this gives one figure per node and no uncertainty: ` +
        'the uncertainty that matters comes from the calibrations, and only the Bayesian relaxed clock propagates it. If you are going to publish ages, run it.'));
    renderAges(res, null);
    el('p7Figure').style.display = '';
    el('p7ExportLog').style.display = 'none';
    redraw();
    commit();
  }

  function showRelaxed(res, relative, secs) {
    const unit = relative ? '' : ' ' + (el('p7Unit').value || 'Ma');
    const roots = res.samples.map(s => s.rootAge);
    const h = Mcmc.hpd(roots, 0.95);
    const cv = Mcmc.mean(res.samples.map(s => s.coefficientOfVariation));
    statTiles('p7DateTiles', [
      [T('Edad de la raíz', 'Root age'), fmtFixed(Mcmc.mean(roots), 3) + unit,
        `HPD 95 % ${fmtFixed(h.lower, 2)} – ${fmtFixed(h.upper, 2)}`, 'ok'],
      [T('Muestras', 'Samples'), fmtNum(res.samples.length, 0), `${fmtNum(res.generations, 0)} ${T('generaciones', 'generations')}`],
      [T('Tasa media', 'Mean rate'), Mcmc.mean(res.samples.map(s => s.meanRate)).toExponential(3), ''],
      [T('Variación de la tasa', 'Rate variation'), fmtFixed(cv, 3), T('coeficiente de variación', 'coefficient of variation'),
        cv > 0.5 ? 'bad' : 'ok'],
      ['ESS ' + T('de la raíz', 'of the root'), fmtFixed(Mcmc.ess(roots), 0), Mcmc.ess(roots) >= 200 ? T('suficiente', 'enough') : T('insuficiente', 'not enough'),
        Mcmc.ess(roots) >= 200 ? 'ok' : 'bad'],
      [T('Tiempo', 'Time'), `${secs.toFixed(1)} s`, ''],
    ]);
    const v = el('p7DateVerdict');
    v.style.display = '';
    v.innerHTML = `<b>${L2('Cómo se lee', 'How to read it')}</b> ` + L2(
      `La raíz tiene una media de ${fmtFixed(Mcmc.mean(roots), 2)}${unit} y un intervalo de máxima densidad del 95 % de ${fmtFixed(h.lower, 2)} a ${fmtFixed(h.upper, 2)}. ` +
      `Ese intervalo <b>es</b> el resultado; la media sola no lo es. El coeficiente de variación de las tasas, ${fmtFixed(cv, 3)}, dice cuánto se aparta el reloj de ser estricto: ` +
      (cv < 0.1 ? 'por debajo de 0.1 el reloj es prácticamente estricto.' : 'por encima de 0.1 las velocidades varían de verdad entre ramas, y usar un reloj estricto habría sesgado las edades.'),
      `The root has a mean of ${fmtFixed(Mcmc.mean(roots), 2)}${unit} and a 95 % highest-density interval from ${fmtFixed(h.lower, 2)} to ${fmtFixed(h.upper, 2)}. ` +
      `That interval <b>is</b> the result; the mean alone is not. The coefficient of variation of the rates, ${fmtFixed(cv, 3)}, says how far the clock is from strict: ` +
      (cv < 0.1 ? 'below 0.1 the clock is practically strict.' : 'above 0.1 the rates really do vary between branches, and a strict clock would have biased the ages.'));
    const intervals = new Map();
    const taxa = taxaOf();
    const splits = Tree.splits(res.tree, taxa.length);
    const F = res.F;
    splits.forEach((node, key) => {
      const k = Tree.flatten(res.tree).nodes.indexOf(node);
      const s = res.summary.find(x => x.node === k);
      if (s) intervals.set(key, { lower: s.lower, upper: s.upper });
    });
    B7.intervals = intervals;
    renderAges(res, res.summary);
    el('p7Figure').style.display = '';
    el('p7ExportLog').style.display = '';
    redraw();
    commit();
  }

  function renderAges(res, summary) {
    const taxa = taxaOf();
    const tree = res.tree;
    const F = Tree.flatten(tree);
    const splits = Tree.splits(tree, taxa.length);
    const rows = [];
    splits.forEach((node, key) => {
      const k = F.nodes.indexOf(node);
      const names = key.split(',').map(Number).map(i => taxa[i]);
      const s = summary ? summary.find(x => x.node === k) : null;
      rows.push({
        clade: names.length <= 3 ? names.join(', ') : `${names.slice(0, 2).join(', ')} … (${names.length})`,
        size: names.length,
        age: node.age != null ? node.age : (res.ages ? res.ages[k] : null),
        lower: s ? s.lower : null,
        upper: s ? s.upper : null,
      });
    });
    rows.sort((a, b) => (b.age || 0) - (a.age || 0));
    buildTable('p7AgeTable', [
      { key: 'clade', label: T('Clado', 'Clade') },
      { key: 'size', label: T('taxones', 'taxa'), num: true },
      { key: 'age', label: T('edad', 'age'), num: true, fmt: v => (v == null ? '—' : fmtFixed(v, 3)) },
      { key: 'lower', label: 'HPD 95 % ' + T('inferior', 'lower'), num: true, fmt: v => (v == null ? '—' : fmtFixed(v, 3)) },
      { key: 'upper', label: 'HPD 95 % ' + T('superior', 'upper'), num: true, fmt: v => (v == null ? '—' : fmtFixed(v, 3)) },
    ], rows, { limit: 60 });
  }

  /* cross-validation for the smoothing parameter */
  function runCV() {
    const msg = el('p7DateMsg');
    clearMessages(msg);
    if (!B7.rooted) return;
    const A = prepare();
    const calibs = B7.calibs.length ? B7.calibs : [{ node: 0, type: 'fixed', value: 1 }];
    const btn = el('p7CV');
    btn.disabled = true;
    el('p7DateProgress').textContent = T('validación cruzada…', 'cross-validation…');
    const w = window.LABG ? LABG.work({ title: T('Validación cruzada de λ', 'Cross-validation of λ'), delay: 300 }) : null;
    phyAfterPaint(() => {
      const cv = Clock.crossValidate(B7.rooted, A.nSites, calibs, { lambdas: [0.01, 0.1, 1, 10, 100, 1000] });
      B7.cv = cv;
      el('p7Lambda').value = cv.best.lambda;
      showMessage(msg, 'success', L2(
        `La validación cruzada elige λ = ${cv.best.lambda}. Valores pequeños dejan que cada rama tenga su velocidad; valores grandes empujan hacia un reloj estricto.`,
        `Cross-validation chooses λ = ${cv.best.lambda}. Small values let each branch have its own rate; large ones push towards a strict clock.`));
      el('p7DateProgress').textContent = '';
      btn.disabled = false;
      if (w) w.done();
    }, w);
  }

  /* ================================================================
     5 · the chronogram
     ================================================================ */
  function redraw() {
    if (!B7.result) return;
    const taxa = taxaOf();
    const tree = Tree.clone(B7.result.tree);
    Tree.ladderize(tree);
    el('p7TreeFig').innerHTML = Plots7.chronogram(tree, {
      labels: taxa, width: 840,
      geo: el('p7Geo').checked && !B7.relative,
      bars: el('p7Bars').checked,
      intervals: B7.intervals || null,
      showAges: el('p7ShowAges').checked,
      unit: B7.relative ? T('profundidad relativa', 'relative depth') : (el('p7Unit').value || 'Ma'),
    });
    /* an age that could not be printed is an age the reader does not get: say so */
    const fm = el('p7FigMsg');
    if (fm) {
      fm.innerHTML = '';
      const n = (Plots7.lastInfo || {}).hiddenAges || 0;
      if (n) {
        showMessage(fm, 'info', L2(
          `Se ocultaron ${n} edad(es) porque se encimaban con otras. Están todas en la tabla y en el paquete; ` +
          'para verlas en la figura, agranda el alto de fila o dibuja menos taxones.',
          `${n} age(s) were hidden because they overlapped. They are all in the table and in the package; ` +
          'to see them on the figure, increase the row height or draw fewer taxa.'));
      }
    }
  }

  function exportAs(kind) {
    const taxa = taxaOf();
    if (!B7.result) return;
    const base = slug((current() || {}).name || 'arbol') + '_fechado';
    const tree = B7.result.tree;
    if (kind === 'newick') download(Tree.writeNewick(tree, { labels: taxa }), base + '.nwk');
    else if (kind === 'nexus') {
      const p = current();
      download(SeqIO.writeNexus(taxa, p ? p.seqs : taxa.map(() => ''), {
        type: p ? p.type : 'dna',
        trees: [{ name: 'cronograma', newick: Tree.writeNewick(tree, { labels: taxa }) }],
      }), base + '.nex');
    } else if (kind === 'svg') download(el('p7TreeFig').innerHTML, base + '.svg', 'image/svg+xml');
    else if (kind === 'ages') {
      const F = Tree.flatten(tree);
      const lines = ['clado,taxones,edad,hpd_inferior,hpd_superior'];
      Tree.splits(tree, taxa.length).forEach((node, key) => {
        const k = F.nodes.indexOf(node);
        const s = B7.result.summary ? B7.result.summary.find(x => x.node === k) : null;
        const names = key.split(',').map(Number).map(i => taxa[i]);
        lines.push(`"${names.join(' ')}",${names.length},${(node.age || 0).toFixed(6)},` +
          `${s ? s.lower.toFixed(6) : ''},${s ? s.upper.toFixed(6) : ''}`);
      });
      download(lines.join('\n') + '\n', base + '_edades.csv', 'text/csv');
    } else if (kind === 'log' && B7.result.samples) {
      download(Mcmc.tracerLog(B7.result.samples), base + '.log', 'text/plain');
    }
  }

  function commit() {
    if (!B7.result) return;
    const taxa = taxaOf();
    state.dated = {
      partition: (current() || {}).name,
      source: B7.source,
      method: B7.method,
      relative: !!B7.relative,
      unit: el('p7Unit').value || 'Ma',
      rootAge: B7.result.rootAge != null ? B7.result.rootAge
        : (B7.result.samples ? Mcmc.mean(B7.result.samples.map(s => s.rootAge)) : null),
      tree: B7.result.tree,
      newick: Tree.writeNewick(B7.result.tree, { labels: taxa }),
      calibrations: B7.calibs.slice(),
      clockTest: B7.test ? { stat: B7.test.stat, df: B7.test.df, p: B7.test.p } : null,
      regression: B7.regression ? { slope: B7.regression.slope, r2: B7.regression.r2, rate: B7.regression.rate } : null,
      intervals: B7.intervals ? [...B7.intervals.entries()] : null,
      samples: B7.result.samples || null,
    };
    enableStep(8, true);
  }

  /* ================================================================
     wiring
     ================================================================ */
  function refresh() {
    const src = sources();
    const has = src.length > 0;
    el('p7NoTree').style.display = has ? 'none' : '';
    ['p7Root'].forEach(id => { el(id).style.display = has ? '' : 'none'; });
    if (!has) { ['p7Test', 'p7Calib', 'p7Date', 'p7Figure'].forEach(id => { el(id).style.display = 'none'; }); return; }
    const sel = el('p7Source');
    sel.innerHTML = src.map(s => `<option value="${s[0]}">${esc(s[1])}</option>`).join('');
    if (!src.some(s => s[0] === B7.source)) B7.source = src[0][0];
    sel.value = B7.source;
    const ps = parts();
    const psel = el('p7Part');
    psel.innerHTML = ps.map((q, i) => `<option value="${i}">${esc(q.name)} · ${q.taxa.length} × ${q.length}</option>`).join('');
    if (B7.part < 0 || B7.part >= ps.length) B7.part = 0;
    psel.value = String(B7.part);
    fillCalTaxa();
  }

  function init() {
    if (!el('panel-7')) return;
    el('p7Source').addEventListener('change', () => {
      B7.source = el('p7Source').value;
      B7.rooted = null; B7.result = null; B7.test = null; B7.intervals = null;
      ['p7Test', 'p7Calib', 'p7Date', 'p7Figure'].forEach(id => { el(id).style.display = 'none'; });
    });
    el('p7Part').addEventListener('change', () => { B7.part = +el('p7Part').value; B7.A = null; });
    el('p7RootMode').addEventListener('change', () => {
      const m = el('p7RootMode').value;
      el('p7OutgroupWrap').style.display = m === 'outgroup' ? '' : 'none';
      el('p7DatesWrap').style.display = m === 'dates' ? '' : 'none';
    });
    el('p7DoRoot').addEventListener('click', doRoot);
    el('p7RunTest').addEventListener('click', runTest);
    el('p7CalType').addEventListener('change', () => {
      const t = el('p7CalType').value;
      el('p7CalMaxWrap').style.display = (t === 'uniform') ? '' : 'none';
      el('p7CalParWrap').style.display = (t === 'normal' || t === 'lognormal' || t === 'exponential') ? '' : 'none';
      el('p7CalMinWrap').style.display = (t === 'normal') ? 'none' : '';
      refreshAdvice();
    });
    el('p7AddCal').addEventListener('click', addCalibration);
    el('p7ClearCal').addEventListener('click', () => { B7.calibs = []; renderCalTable(); refreshAdvice(); });

    /* el asistente */
    if (el('p7CalSource') && window.CalibHelp) {
      fillSources();
      el('p7CalSource').addEventListener('change', () => { refreshSourceHint(); refreshAdvice(); });
      el('p7CalSuggest').addEventListener('click', suggestPrior);
      ['p7CalMin', 'p7CalMax', 'p7CalMean', 'p7CalSd'].forEach(id => {
        const e2 = el(id);
        if (e2) e2.addEventListener('input', refreshAdvice);
      });
      if (el('p7CalTaxa')) el('p7CalTaxa').addEventListener('change', refreshAdvice);
      document.addEventListener('langchange', () => { fillSources(); refreshAdvice(); });
      refreshAdvice();
    }
    el('p7Method').addEventListener('change', () => {
      const m = el('p7Method').value;
      el('p7LambdaWrap').style.display = m === 'pl' ? '' : 'none';
      el('p7GensWrap').style.display = m === 'relaxed' ? '' : 'none';
    });
    el('p7CV').addEventListener('click', runCV);
    el('p7RunDate').addEventListener('click', runDate);
    el('p7CancelDate').addEventListener('click', () => { B7.cancelled = true; });
    ['p7Geo', 'p7Bars', 'p7ShowAges'].forEach(id => el(id).addEventListener('change', redraw));
    el('p7Redraw').addEventListener('click', redraw);
    els('[data-p7export]').forEach(b => b.addEventListener('click', () => exportAs(b.dataset.p7export)));
    el('p7ToBlock8').addEventListener('click', () => {
      const b = document.querySelector('.step-btn[data-step="8"]');
      if (b && !b.disabled) goStep(8);
      else showMessage(el('p7DateMsg'), 'info', L2(
        'El Bloque 8 (diversificación) llega en la etapa siguiente; este cronograma será su punto de partida.',
        'Block 8 (diversification) arrives in the next stage; this chronogram will be its starting point.'));
    });
    document.addEventListener('stepchange', e => { if (e.detail.step === 7) refresh(); });
    document.addEventListener('langchange', () => { refresh(); if (B7.calibs.length) renderCalTable(); if (B7.result) { renderAges(B7.result, B7.result.summary); redraw(); } });
    document.addEventListener('themechange', () => { if (B7.result) redraw(); if (B7.regression) el('p7RegFig').innerHTML = Plots7.regression(B7.regression, {}); });
  }
  document.addEventListener('DOMContentLoaded', init);

  Object.assign(B7, { refresh, doRoot, runTest, addCalibration, runDate, runCV, redraw, commit, current });
})();
