/* PhylogenyPro — molecular phylogenetics in the browser.
   Copyright (C) 2026  Luis Ángel Barrera-Guzmán

   This program is free software: you can redistribute it and/or modify it under
   the terms of the GNU General Public License as published by the Free Software
   Foundation, either version 3 of the License, or (at your option) any later
   version. It is distributed in the hope that it will be useful, but WITHOUT ANY
   WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
   A PARTICULAR PURPOSE. See the GNU General Public License, in the file LICENSE
   at the root of this program, or <https://www.gnu.org/licenses/>. */

/* PhylogenyPro — Block 8: diversification.

   The engine is js/diversify.js. What this file adds is the reading: every
   number here has a way of being wrong that is worth saying out loud, and the
   block says it. A negative γ means a slowdown or missing species, and there is
   no telling which without the Monte Carlo test. A density-dependent model
   winning by AIC does not prove that diversity is bounded. And the estimator
   that needs only an age and a richness needs an extinction rate nobody knows,
   so it is reported at several. */

(function () {

  const B8 = {
    source: 'dated', tree: null, ltt: null, gamma: null,
    mccr: null, models: null, dr: null, ms: null,
    cancelled: false,
  };
  window.B8 = B8;

  function sources() {
    const out = [];
    if (state.dated && state.dated.tree) {
      out.push(['dated', T(`Bloque 7 · ${state.dated.method === 'relaxed' ? 'reloj relajado' : state.dated.method === 'pl' ? 'verosimilitud penalizada' : 'mínimos cuadrados'}`,
        `Block 7 · ${state.dated.method === 'relaxed' ? 'relaxed clock' : state.dated.method === 'pl' ? 'penalised likelihood' : 'least squares'}`), state.dated.tree]);
    }
    return out;
  }
  function chosenTree() {
    const s = sources().find(x => x[0] === B8.source) || sources()[0];
    return s ? Tree.clone(s[2]) : null;
  }
  function taxaOf() {
    if (state.dated && state.dated.partition && state.data) {
      const p = (state.data.parts || []).find(q => q.name === state.dated.partition);
      if (p) return p.taxa;
    }
    return state.data ? state.data.taxa : [];
  }
  function unit() { return (state.dated && state.dated.unit) || 'Ma'; }
  function relative() { return !!(state.dated && state.dated.relative); }

  /* ================================================================
     1 · lineages through time and gamma
     ================================================================ */
  function showLtt() {
    const tree = chosenTree();
    if (!tree) return;
    B8.tree = tree;
    const L = Diversify.ltt(tree);
    const gm = Diversify.gammaStat(tree);
    B8.ltt = L; B8.gamma = gm;
    const n = Tree.tips(tree).length;
    const y = Diversify.yule(tree);

    statTiles('p8LttTiles', [
      [T('Especies en el árbol', 'Species in the tree'), n, ''],
      [T('Edad de la corona', 'Crown age'), fmtFixed(L.times[0], 3) + (relative() ? '' : ' ' + unit()), relative() ? T('relativa', 'relative') : ''],
      ['γ', gm ? fmtFixed(gm.gamma, 4) : '—',
        gm ? (gm.gamma < -1.645 ? T('desaceleración', 'slowdown') : gm.gamma > 1.645 ? T('aceleración', 'speed-up') : T('compatible con tasa constante', 'compatible with a constant rate')) : '',
        gm ? (Math.abs(gm.gamma) > 1.645 ? 'bad' : 'ok') : ''],
      [T('Tasa de nacimiento puro', 'Pure-birth rate'), y.lambda.toExponential(3),
        relative() ? T('por unidad relativa', 'per relative unit') : T('especies por', 'species per') + ' ' + unit()],
    ]);

    const v = el('p8GammaVerdict');
    v.style.display = '';
    const g0 = gm ? gm.gamma : 0;
    v.innerHTML = `<b>${L2('Cómo se lee', 'How to read it')}</b> ` + (!gm
      ? L2('Hacen falta al menos cuatro especies para calcular γ.', 'At least four species are needed to compute γ.')
      : g0 < -1.645
        ? L2(
          `γ = ${fmtFixed(g0, 3)}, por debajo de −1.645: los nodos se acumulan cerca de la raíz y la curva se dobla. Esa es la firma de una desaceleración… <b>y también la que deja el muestreo incompleto</b>, que pierde sobre todo nodos recientes. Antes de llamarlo desaceleración hay que pasar la prueba de abajo con el número verdadero de especies del grupo.`,
          `γ = ${fmtFixed(g0, 3)}, below −1.645: the nodes pile up near the root and the curve bends. That is the signature of a slowdown… <b>and also the one incomplete sampling leaves</b>, because sampling loses recent nodes above all. Before calling it a slowdown, run the test below with the clade's true number of species.`)
        : g0 > 1.645
          ? L2(
            `γ = ${fmtFixed(g0, 3)}, por encima de 1.645: los nodos están más cerca del presente que bajo tasa constante. Esto aparece cuando la diversificación se acelera y también cuando hay extinción, que borra los linajes viejos y deja los recientes.`,
            `γ = ${fmtFixed(g0, 3)}, above 1.645: the nodes sit closer to the present than a constant rate predicts. This appears when diversification speeds up and also when there is extinction, which erases old lineages and leaves recent ones.`)
          : L2(
            `γ = ${fmtFixed(g0, 3)} está dentro de lo que una tasa constante produce. Eso no demuestra que la tasa sea constante: con pocas especies la prueba tiene poca potencia, y con muchas puede haber cambios que se compensan.`,
            `γ = ${fmtFixed(g0, 3)} is within what a constant rate produces. That does not prove the rate was constant: with few species the test has little power, and with many, changes can cancel out.`));

    el('p8LttFig').innerHTML = Plots8.ltt(L, {
      xLabel: relative() ? T('profundidad relativa', 'relative depth') : T(`tiempo antes del presente (${unit()})`, `time before the present (${unit()})`),
    });
    /* sensible defaults for the other cards */
    if (!+el('p8Total').value) el('p8Total').value = n;
    if (!+el('p8N').value) el('p8N').value = n;
    if (!+el('p8Age').value) el('p8Age').value = fmtFixed(L.times[0], 4).replace(/[^0-9.\-]/g, '');
    ['p8Mccr', 'p8Models', 'p8MS', 'p8DR', 'p8Export'].forEach(id => { el(id).style.display = ''; });
  }

  /* ================================================================
     2 · the Monte Carlo constant-rates test
     ================================================================ */
  function runMccr() {
    const msg = el('p8MccrMsg');
    clearMessages(msg);
    if (!B8.tree) return;
    const n = Tree.tips(B8.tree).length;
    const total = Math.max(n, +el('p8Total').value || n);
    const reps = Math.max(50, +el('p8Reps').value || 500);
    const btn = el('p8RunMccr'), cancel = el('p8CancelMccr'), prog = el('p8MccrProgress');
    btn.disabled = true; cancel.style.display = ''; B8.cancelled = false;
    if (total === n) showMessage(msg, 'info', L2(
      'Con el árbol completo la prueba solo confirma lo que dice γ. Su razón de ser es el caso incompleto: escribe arriba cuántas especies tiene el grupo de verdad.',
      'With a complete tree the test only confirms what γ says. Its reason to exist is the incomplete case: write above how many species the clade really has.'));
    const t0 = performance.now();
    setTimeout(() => {
      const lambda = Diversify.yule(B8.tree).lambda;
      const res = Diversify.mccr(B8.tree, {
        total, reps, lambda, seed: 5,
        onProgress: (i, r) => { prog.textContent = `${i}/${r} · ${((performance.now() - t0) / 1000).toFixed(0)} s`; },
        cancelled: () => B8.cancelled,
      });
      B8.mccr = res;
      statTiles('p8MccrTiles', [
        ['γ ' + T('observado', 'observed'), fmtFixed(res.gamma, 4), ''],
        [T('Especies', 'Species'), `${res.nSampled} ${T('de', 'of')} ${res.nTotal}`,
          res.nTotal > res.nSampled ? fmtPct(res.nSampled / res.nTotal, 0) + ' ' + T('muestreado', 'sampled') : T('completo', 'complete')],
        [T('Valor crítico al 5 %', '5 % critical value'), fmtFixed(res.critical, 4), T('de la distribución nula', 'of the null distribution')],
        ['p', fmtP(res.p), res.p < 0.05 ? T('desaceleración significativa', 'significant slowdown') : T('no significativa', 'not significant'),
          res.p < 0.05 ? 'bad' : 'ok'],
        [T('Réplicas', 'Replicates'), fmtNum(res.reps, 0), `${((performance.now() - t0) / 1000).toFixed(1)} s`],
      ]);
      const v = el('p8MccrVerdict');
      v.style.display = '';
      const mean = Mcmc.mean(res.nullDistribution);
      v.innerHTML = `<b>${L2('Cómo se lee', 'How to read it')}</b> ` + L2(
        `Con ${res.nSampled} especies de ${res.nTotal}, una tasa constante produce un γ de ${fmtFixed(mean, 3)} en promedio —negativo, sin que nada se haya desacelerado— y por debajo de ${fmtFixed(res.critical, 3)} solo el 5 % de las veces. ` +
        (res.p < 0.05
          ? `El γ observado, ${fmtFixed(res.gamma, 3)}, cae por debajo: la desaceleración no se explica por las especies que faltan.`
          : `El γ observado, ${fmtFixed(res.gamma, 3)}, no cae por debajo: lo que se ve es compatible con una tasa constante mal muestreada.`),
        `With ${res.nSampled} species of ${res.nTotal}, a constant rate produces a γ of ${fmtFixed(mean, 3)} on average — negative, with nothing having slowed down — and below ${fmtFixed(res.critical, 3)} only 5 % of the time. ` +
        (res.p < 0.05
          ? `The observed γ, ${fmtFixed(res.gamma, 3)}, falls below it: the slowdown is not explained by the missing species.`
          : `The observed γ, ${fmtFixed(res.gamma, 3)}, does not fall below it: what you see is compatible with a constant rate, poorly sampled.`));
      el('p8MccrPane').style.display = '';
      el('p8MccrFig').innerHTML = Plots6.histogram(res.nullDistribution, {
        width: 700, height: 240, bins: 40,
        xLabel: 'γ', yLabel: T('réplicas', 'replicates'),
        marks: [{ x: res.critical, label: T('5 %', '5 %') }, { x: res.gamma, label: T('observado', 'observed') }],
      });
      btn.disabled = false; cancel.style.display = 'none'; prog.textContent = '';
      commit();
    }, 30);
  }

  /* ================================================================
     3 · the models
     ================================================================ */
  function runModels() {
    const msg = el('p8ModelMsg');
    clearMessages(msg);
    if (!B8.tree) return;
    const btn = el('p8RunModels');
    btn.disabled = true;
    el('p8ModelsProgress').textContent = T('ajustando…', 'fitting…');
    setTimeout(() => {
      const cmp = Diversify.compare(B8.tree);
      B8.models = cmp;
      buildTable('p8ModelTable', [
        { key: 'name', label: T('Modelo', 'Model') },
        { key: 'params', label: T('parámetros ajustados', 'fitted parameters') },
        { key: 'loglik', label: 'lnL', num: true, fmt: v => fmtFixed(v, 4) },
        { key: 'k', label: 'k', num: true },
        { key: 'AIC', label: 'AIC', num: true, fmt: v => fmtFixed(v, 3) },
        { key: 'dAIC', label: 'ΔAIC', num: true, fmt: v => fmtFixed(v, 3) },
        { key: 'w', label: T('peso', 'weight'), num: true, fmt: v => fmtFixed(v, 3) },
      ], cmp.rows.map(r => ({
        name: modelName(r.name), params: describeFit(r), loglik: r.loglik, k: r.k,
        AIC: r.AIC, dAIC: r.dAIC, w: r.w,
      })), {});

      const best = cmp.rows[0];
      const second = cmp.rows[1];
      const v = el('p8ModelVerdict');
      v.style.display = '';
      v.innerHTML = `<b>${L2('Cómo se lee', 'How to read it')}</b> ` + L2(
        `Gana <b>${esc(modelName(best.name))}</b> con un peso de Akaike de ${fmtPct(best.w, 0)}. ` +
        (second && second.dAIC < 2
          ? `El siguiente, ${esc(modelName(second.name))}, queda a ${fmtFixed(second.dAIC, 2)} unidades de AIC: por debajo de dos, los datos no distinguen entre ellos y presentar uno como «el modelo» sería exagerar.`
          : `El siguiente queda a ${fmtFixed(second ? second.dAIC : 0, 2)} unidades, así que la diferencia sí es apreciable.`) +
        ' Y una advertencia que vale para todos: que un modelo dependiente de la densidad gane no demuestra que la diversidad tenga un techo; la extinción produce curvas parecidas y estos datos rara vez las separan.',
        `<b>${esc(modelName(best.name))}</b> wins with an Akaike weight of ${fmtPct(best.w, 0)}. ` +
        (second && second.dAIC < 2
          ? `The next, ${esc(modelName(second.name))}, is ${fmtFixed(second.dAIC, 2)} AIC units away: below two, the data cannot tell them apart and presenting one as "the" model would overstate the case.`
          : `The next is ${fmtFixed(second ? second.dAIC : 0, 2)} units away, so the difference is real.`) +
        ' And a warning that holds for all of them: a density-dependent model winning does not prove that diversity has a ceiling; extinction produces similar curves and these data rarely separate them.');

      drawRates(best);
      btn.disabled = false;
      el('p8ModelsProgress').textContent = '';
      commit();
    }, 30);
  }
  function modelName(n) {
    return T({
      'pure birth (Yule)': 'nacimiento puro (Yule)',
      'birth–death': 'nacimiento–muerte',
      'density dependent, logistic': 'dependiente de la densidad, logística',
      'density dependent, exponential': 'dependiente de la densidad, exponencial',
      'two rates, shift estimated': 'dos tasas, con el cambio estimado',
    }[n] || n, n);
  }
  function describeFit(r) {
    const f = r.fit;
    if (r.name.indexOf('pure birth') === 0) return `λ = ${f.r.toExponential(3)}`;
    if (r.name.indexOf('birth–death') === 0) return `λ = ${f.lambda.toExponential(3)} · μ = ${f.mu.toExponential(3)} · ε = ${fmtFixed(f.a, 3)}`;
    if (r.name.indexOf('logistic') >= 0) return `r = ${f.r.toExponential(3)} · K = ${fmtFixed(f.K, 1)}`;
    if (r.name.indexOf('exponential') >= 0) return `r = ${f.r.toExponential(3)} · x = ${fmtFixed(f.x, 3)}`;
    if (r.name.indexOf('two rates') >= 0) {
      return `λ₁ = ${f.r1.toExponential(3)} → λ₂ = ${f.r2.toExponential(3)} ${T('en', 'at')} ${fmtFixed(f.shift, 2)}`;
    }
    return '';
  }
  function drawRates(best) {
    const pane = el('p8RatePane');
    const f = best.fit;
    const n = Tree.tips(B8.tree).length;
    let pts = null, title = '', xLabel = '';
    if (best.name.indexOf('logistic') >= 0) {
      pts = []; for (let i = 2; i <= n; i++) pts.push({ x: i, y: Math.max(0, f.r * (1 - i / f.K)) });
      title = T('La tasa baja conforme se acumulan linajes', 'The rate falls as lineages accumulate');
      xLabel = T('linajes', 'lineages');
    } else if (best.name.indexOf('exponential') >= 0) {
      pts = []; for (let i = 2; i <= n; i++) pts.push({ x: i, y: f.r * Math.pow(i, -f.x) });
      title = T('La tasa baja conforme se acumulan linajes', 'The rate falls as lineages accumulate');
      xLabel = T('linajes', 'lineages');
    } else if (best.name.indexOf('two rates') >= 0) {
      const maxT = B8.ltt.times[0];
      pts = [];
      for (let i = 0; i <= 100; i++) {
        const t = maxT * (1 - i / 100);
        pts.push({ x: t, y: t >= f.shift ? f.r1 : f.r2 });
      }
      title = T('La tasa cambia en un momento estimado', 'The rate changes at an estimated moment');
      xLabel = relative() ? T('profundidad relativa', 'relative depth') : T(`tiempo antes del presente (${unit()})`, `time before the present (${unit()})`);
    }
    if (!pts) { pane.style.display = 'none'; return; }
    pane.style.display = '';
    el('p8RateTitle').textContent = title;
    el('p8RateFig').innerHTML = Plots8.rateCurve(pts, { xLabel });
  }

  /* ================================================================
     4 · age and richness
     ================================================================ */
  function runMS() {
    const msg = el('p8MSMsg');
    clearMessages(msg);
    const age = +el('p8Age').value, n = +el('p8N').value;
    const crown = el('p8Crown').value === 'crown';
    if (!(age > 0) || !(n > 1)) {
      showMessage(msg, 'warning', L2('Hace falta una edad positiva y al menos dos especies.', 'A positive age and at least two species are needed.'));
      return;
    }
    const rows = [0, 0.45, 0.9].map(e => {
      const r = Diversify.magallonSanderson(age, n, e, crown);
      return {
        epsilon: e,
        label: e === 0 ? T('sin extinción', 'no extinction') : e === 0.9 ? T('mucha extinción', 'a lot of extinction') : T('extinción moderada', 'moderate extinction'),
        rate: r,
        doubling: r > 0 ? Math.log(2) / r : null,
      };
    });
    B8.ms = { age, n, crown, rows };
    buildTable('p8MSTable', [
      { key: 'epsilon', label: 'ε = μ/λ', num: true, fmt: v => fmtFixed(v, 2) },
      { key: 'label', label: T('supuesto', 'assumption') },
      { key: 'rate', label: T('tasa neta de diversificación', 'net diversification rate'), num: true, fmt: v => v.toExponential(4) },
      { key: 'doubling', label: T('tiempo de duplicación', 'doubling time'), num: true, fmt: v => (v == null ? '—' : fmtFixed(v, 3)) },
    ], rows, {});
    showMessage(msg, 'success', L2(
      `Con ${n} especies y una edad de ${fmtFixed(age, 2)} ${relative() ? '' : unit()} en ${crown ? 'la corona' : 'el tallo'}, la tasa va de ${rows[0].rate.toExponential(3)} suponiendo que nada se extinguió a ${rows[2].rate.toExponential(3)} suponiendo mucha extinción. ` +
      'Esa horquilla es el resultado: informar solo el primer número es lo que hace que las comparaciones entre clados no se sostengan.',
      `With ${n} species and an age of ${fmtFixed(age, 2)} ${relative() ? '' : unit()} at the ${crown ? 'crown' : 'stem'}, the rate runs from ${rows[0].rate.toExponential(3)} assuming nothing went extinct to ${rows[2].rate.toExponential(3)} assuming a lot did. ` +
      'That range is the result: reporting only the first number is what makes comparisons between clades fall apart.'));
    commit();
  }

  /* ================================================================
     5 · per-species rates
     ================================================================ */
  function runDR() {
    if (!B8.tree) return;
    const taxa = taxaOf();
    const dr = Diversify.drStatistic(B8.tree);
    B8.dr = dr;
    const vals = dr.map(d => d.dr);
    statTiles('p8DRTiles', [
      [T('Mediana', 'Median'), fmtFixed(Mcmc.quantile(vals, 0.5), 4), ''],
      [T('Mínimo', 'Smallest'), fmtFixed(Math.min.apply(null, vals), 4), taxa[dr.reduce((a, b) => (b.dr < a.dr ? b : a)).tip] || ''],
      [T('Máximo', 'Largest'), fmtFixed(Math.max.apply(null, vals), 4), taxa[dr.reduce((a, b) => (b.dr > a.dr ? b : a)).tip] || ''],
      [T('Razón máximo/mínimo', 'Largest over smallest'), fmtFixed(Math.max.apply(null, vals) / Math.min.apply(null, vals), 2), ''],
    ]);
    el('p8DRPane').style.display = '';
    el('p8DRFig').innerHTML = Plots8.ranked(dr.map(d => ({ name: taxa[d.tip] || ('sp' + d.tip), value: d.dr })), {});
    buildTable('p8DRTable', [
      { key: 'name', label: T('Especie', 'Species') },
      { key: 'dr', label: 'DR', num: true, fmt: v => fmtFixed(v, 5) },
      { key: 'es', label: T('longitud de rama equitativa', 'equal-splits length'), num: true, fmt: v => fmtFixed(v, 5) },
    ], dr.map(d => ({ name: taxa[d.tip] || ('sp' + d.tip), dr: d.dr, es: d.es })).sort((a, b) => b.dr - a.dr), { limit: 60 });
    commit();
  }

  /* ================================================================
     export and hand-off
     ================================================================ */
  function exportAs(kind) {
    const taxa = taxaOf();
    const base = 'diversificacion';
    if (kind === 'ltt' && B8.ltt) {
      const lines = ['tiempo,linajes'];
      B8.ltt.times.forEach((t, i) => lines.push(`${t.toFixed(6)},${B8.ltt.lineages[i]}`));
      download(lines.join('\n') + '\n', base + '_ltt.csv', 'text/csv');
    } else if (kind === 'models' && B8.models) {
      const lines = ['modelo,parametros,lnL,k,AIC,dAIC,peso'];
      B8.models.rows.forEach(r => lines.push(
        `"${modelName(r.name)}","${describeFit(r)}",${r.loglik.toFixed(6)},${r.k},${r.AIC.toFixed(6)},${r.dAIC.toFixed(6)},${r.w.toFixed(6)}`));
      download(lines.join('\n') + '\n', base + '_modelos.csv', 'text/csv');
    } else if (kind === 'dr' && B8.dr) {
      const lines = ['especie,DR,longitud_equitativa'];
      B8.dr.forEach(d => lines.push(`"${taxa[d.tip] || ('sp' + d.tip)}",${d.dr.toFixed(8)},${d.es.toFixed(8)}`));
      download(lines.join('\n') + '\n', base + '_dr.csv', 'text/csv');
    } else if (kind === 'svg') {
      download(el('p8LttFig').innerHTML, base + '_ltt.svg', 'image/svg+xml');
    }
  }

  function commit() {
    if (!B8.tree) return;
    /* the model with the lowest AIC, by its English name; the report translates it */
    const best = B8.models && B8.models.rows.length
      ? B8.models.rows.reduce((a, r) => (r.AIC < a.AIC ? r : a)).name : null;
    state.diversification = {
      source: B8.source,
      nTips: Tree.tips(B8.tree).length,
      crownAge: B8.ltt ? B8.ltt.times[0] : null,
      gamma: B8.gamma ? B8.gamma.gamma : null,
      mccr: B8.mccr ? { gamma: B8.mccr.gamma, p: B8.mccr.p, critical: B8.mccr.critical, total: B8.mccr.nTotal, reps: B8.mccr.reps } : null,
      models: B8.models ? B8.models.rows.map(r => ({ name: r.name, loglik: r.loglik, k: r.k, AIC: r.AIC, dAIC: r.dAIC, w: r.w, fit: r.fit })) : null,
      best,
      magallonSanderson: B8.ms || null,
      dr: B8.dr ? B8.dr.map(d => ({ tip: d.tip, dr: d.dr })) : null,
    };
    enableStep(9, true);
  }

  /* ================================================================
     wiring
     ================================================================ */
  function refresh() {
    const src = sources();
    const has = src.length > 0;
    el('p8NoTree').style.display = has ? 'none' : '';
    el('p8Ltt').style.display = has ? '' : 'none';
    if (!has) { ['p8Mccr', 'p8Models', 'p8MS', 'p8DR', 'p8Export'].forEach(id => { el(id).style.display = 'none'; }); return; }
    const sel = el('p8Source');
    sel.innerHTML = src.map(s => `<option value="${s[0]}">${esc(s[1])}</option>`).join('');
    sel.value = B8.source;
    showLtt();
  }

  function init() {
    if (!el('panel-8')) return;
    el('p8Source').addEventListener('change', () => { B8.source = el('p8Source').value; showLtt(); });
    el('p8RunMccr').addEventListener('click', runMccr);
    el('p8CancelMccr').addEventListener('click', () => { B8.cancelled = true; });
    el('p8RunModels').addEventListener('click', runModels);
    el('p8RunMS').addEventListener('click', runMS);
    el('p8RunDR').addEventListener('click', runDR);
    els('[data-p8export]').forEach(b => b.addEventListener('click', () => exportAs(b.dataset.p8export)));
    el('p8ToBlock9').addEventListener('click', () => {
      const b = document.querySelector('.step-btn[data-step="9"]');
      if (b && !b.disabled) goStep(9);
      else showMessage(el('p8LttMsg'), 'info', L2(
        'El Bloque 9 (evolución de caracteres) llega en la etapa siguiente.',
        'Block 9 (trait evolution) arrives in the next stage.'));
    });
    document.addEventListener('stepchange', e => { if (e.detail.step === 8) refresh(); });
    document.addEventListener('langchange', () => { if (sources().length) { refresh(); if (B8.models) runModels(); } });
    document.addEventListener('themechange', () => { if (B8.ltt) showLtt(); });
  }
  document.addEventListener('DOMContentLoaded', init);

  Object.assign(B8, { refresh, showLtt, runMccr, runModels, runMS, runDR, commit });
})();
