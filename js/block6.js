/* PhylogenyPro — molecular phylogenetics in the browser.
   Copyright (C) 2026  Luis Ángel Barrera-Guzmán

   This program is free software: you can redistribute it and/or modify it under
   the terms of the GNU General Public License as published by the Free Software
   Foundation, either version 3 of the License, or (at your option) any later
   version. It is distributed in the hope that it will be useful, but WITHOUT ANY
   WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
   A PARTICULAR PURPOSE. See the GNU General Public License, in the file LICENSE
   at the root of this program, or <https://www.gnu.org/licenses/>. */

/* PhylogenyPro — Block 6: Bayesian inference.

   The sampler lives in js/mcmc.js and the parallelism in js/pool.js; this file
   is the part the user touches: writing down the priors, launching independent
   runs, reading whether they converged before believing anything, summarising
   the posterior as a tree with clade probabilities, and comparing models by
   Bayes factors.

   One decision is deliberate and worth stating: the diagnostics are shown
   BEFORE the tree, and the tree card says plainly when the runs have not
   converged. A posterior probability from a chain that has not mixed is not a
   weak result, it is a wrong one. */

(function () {

  const B6 = {
    part: -1, A: null, runs: null, summary: null, bf: null,
    cancelled: false, running: false, traceKey: 'lnL',
  };
  window.B6 = B6;

  function parts() { return (state.data && state.data.parts) || []; }
  function current() { return parts()[B6.part] || null; }
  function prepare() {
    const p = current();
    if (!p) return null;
    if (!B6.A || B6.A.__part !== p.name) {
      B6.A = Like.compress(p.seqs, p.type);
      B6.A.__part = p.name;
    }
    return B6.A;
  }

  const DNA_MODELS = ['JC', 'F81', 'K80', 'HKY', 'TrN', 'K81uf', 'TIM', 'TVM', 'SYM', 'GTR'];
  const AA_MODELS = ['LG', 'WAG', 'JTT', 'Dayhoff', 'cpREV', 'mtREV24', 'VT', 'Blosum62'];

  function fillModels() {
    const p = current();
    const sel = el('p6Model');
    if (!sel || !p) return;
    const prev = sel.value;
    const list = p.type === 'aa' ? AA_MODELS : p.type === 'morph' ? ['Mk', 'Mkv'] : DNA_MODELS;
    sel.innerHTML = list.map(m => `<option value="${m}">${m}</option>`).join('');
    sel.value = list.indexOf(prev) >= 0 ? prev : (p.type === 'aa' ? 'LG' : p.type === 'morph' ? 'Mkv' : 'GTR');
    const bf = el('p6BFModels');
    if (bf) {
      bf.innerHTML = list.map(m => `<option value="${m}">${m}</option>`).join('');
      [...bf.options].forEach(o => { o.selected = (o.value === 'JC' || o.value === sel.value); });
    }
  }

  function specFromUI(modelName) {
    const p = current();
    const A = prepare();
    const model = modelName || el('p6Model').value;
    const spec = { type: p.type === 'morph' ? 'morph' : p.type, model, ncat: 1 };
    if (p.type === 'dna') spec.rates = [1, 2, 1, 1, 2, 1];
    if (p.type !== 'morph') spec.freqs = Array.from(A.freqs);
    if (el('p6Gamma').checked) { spec.alpha = 0.5; spec.ncat = 4; }
    if (el('p6Inv').checked) spec.pInv = 0.2;
    return spec;
  }
  function priorsFromUI() {
    return {
      branchMean: Math.max(1e-4, +el('p6BranchMean').value || 0.1),
      alphaMean: Math.max(1e-3, +el('p6AlphaMean').value || 1),
      rateMean: Math.max(1e-3, +el('p6RateMean').value || 1),
    };
  }
  function modelLabel(s) { return s.model + (s.alpha ? '+G' : '') + (s.pInv > 0 ? '+I' : ''); }

  function useBlock5() {
    const m = state.ml;
    if (!m || !m.spec) {
      showMessage(el('p6Msg'), 'info', L2('Todavía no hay un modelo del Bloque 5.', 'There is no model from Block 5 yet.'));
      return;
    }
    const sel = el('p6Model');
    if ([...sel.options].some(o => o.value === m.spec.model)) sel.value = m.spec.model;
    el('p6Gamma').checked = !!m.spec.alpha;
    el('p6Inv').checked = !!(m.spec.pInv > 0);
    showMessage(el('p6Msg'), 'success', L2(
      `Modelo del Bloque 5: <b>${esc(modelLabel(m.spec))}</b>.`, `Model from Block 5: <b>${esc(modelLabel(m.spec))}</b>.`));
  }

  /* what the priors mean, spelled out with the numbers the user chose */
  function renderPriorNote() {
    const pr = priorsFromUI();
    const p = current();
    const nB = p ? 2 * p.taxa.length - 3 : 0;
    const box = el('p6PriorNote');
    if (!box) return;
    /* The branch-length prior is the one that bites. On divergent data an
       exponential of mean 0.1 — MrBayes' default — pulls the whole tree in and
       drags the substitution parameters with it (Brown et al. 2010; Marshall
       2010). If Block 5 has measured the tree, the comparison can be made here
       instead of after the run. */
    let warn = '';
    if (state.ml && state.ml.tree && state.ml.partition === (p && p.name)) {
      const mlLen = Tree.totalLength(state.ml.tree);
      const expected = nB * pr.branchMean;
      if (mlLen > 4 * expected || mlLen < expected / 4) {
        warn = ' ' + L2(
          `<b>Ojo:</b> el árbol de máxima verosimilitud del Bloque 5 mide ${fmtFixed(mlLen, 2)}, y este prior espera ${fmtFixed(expected, 2)}. ` +
          `Esa distancia no es un detalle: el prior tirará de las ramas hacia su media y arrastrará consigo a κ y a α. ` +
          `Para estos datos una media de rama cercana a ${fmtFixed(mlLen / nB, 2)} es mucho más razonable.`,
          `<b>Careful:</b> the maximum-likelihood tree of Block 5 measures ${fmtFixed(mlLen, 2)}, and this prior expects ${fmtFixed(expected, 2)}. ` +
          `That gap is not a detail: the prior will pull the branches towards its mean and drag κ and α with them. ` +
          `For these data a mean branch length near ${fmtFixed(mlLen / nB, 2)} is far more reasonable.`);
      }
    }
    box.innerHTML = `<b>${L2('Lo que estás suponiendo', 'What you are assuming')}</b> ` + warn + ' ' + L2(
      `Cada una de las ${nB} ramas tiene de entrada una media de ${fmtFixed(pr.branchMean, 3)} sustituciones por sitio, así que el árbol entero mide en promedio ${fmtFixed(nB * pr.branchMean, 2)} antes de mirar los datos. ` +
      `Las ${Math.max(0, nB)} ramas son independientes entre sí. Todas las topologías arrancan igual de probables. ` +
      `Con datos abundantes el prior de las ramas casi no se nota; con pocos sitios sí, y conviene probar otro valor y ver si el árbol cambia.`,
      `Each of the ${nB} branches starts with a mean of ${fmtFixed(pr.branchMean, 3)} substitutions per site, so the whole tree measures ${fmtFixed(nB * pr.branchMean, 2)} on average before the data are seen. ` +
      `The ${Math.max(0, nB)} branches are independent of each other. Every topology starts equally likely. ` +
      `With plenty of data the branch prior barely shows; with few sites it does, and it is worth trying another value and seeing whether the tree changes.`);
  }

  /* ================================================================
     the starting tree
     ================================================================ */
  function startingTree() {
    const p = current();
    const which = el('p6Start').value;
    if (which === 'ml' && state.ml && state.ml.tree && state.ml.partition === p.name) return Tree.clone(state.ml.tree);
    if (which === 'random') {
      const order = Array.from({ length: p.taxa.length }, (_, i) => i);
      const r = rng((Date.now() & 0xffff) + 1);
      for (let i = order.length - 1; i > 0; i--) { const j = Math.floor(r() * (i + 1)); const t = order[i]; order[i] = order[j]; order[j] = t; }
      let tree = { children: order.slice(0, 3).map(i => ({ node: { tip: i, label: p.taxa[i], children: [] }, len: 0.05 })) };
      for (let k = 3; k < order.length; k++) {
        const edges = Tree.allEdges(tree);
        const at = edges[Math.floor(r() * edges.length)];
        at.edge.node = {
          children: [{ node: at.edge.node, len: (at.edge.len || 0.05) / 2 },
            { node: { tip: order[k], label: p.taxa[order[k]], children: [] }, len: 0.05 }],
        };
      }
      return tree;
    }
    const D = Dist.matrix(p.seqs, p.type === 'aa' ? 'poisson' : p.type === 'morph' ? 'p' : 'jc', { type: p.type }).D;
    return Tree.bionj(D, p.taxa);
  }

  /* ================================================================
     2 · running the chains
     ================================================================ */
  function go() {
    const p = current();
    if (!p) return;
    prepare();
    const msg = el('p6Msg');
    clearMessages(msg);
    const btn = el('p6Go'), cancel = el('p6Cancel'), prog = el('p6Progress');
    const gens = Math.max(1000, +el('p6Gens').value || 200000);
    const nRuns = Math.max(1, +el('p6Runs').value || 2);
    const chains = Math.max(1, +el('p6Chains').value || 4);
    const burnPct = Math.min(90, Math.max(0, +el('p6Burnin').value || 25));
    const every = Math.max(1, +el('p6Every').value || 100);
    const deltaT = Math.max(0.01, +el('p6DeltaT').value || 0.2);
    const spec = specFromUI();
    const start = startingTree();

    btn.disabled = true; cancel.style.display = ''; B6.cancelled = false; B6.running = true;
    const threads = Pool.available();
    showMessage(msg, 'info', L2(
      `${nRuns} corridas de ${fmtNum(gens, 0)} generaciones, ${chains} cadenas cada una, ` +
      `${threads ? `en ${Math.min(nRuns, Pool.cores())} hilos` : 'en la página (este navegador no permite hilos)'}. ` +
      'Las corridas son independientes: es la única manera de comprobar después si convergieron.',
      `${nRuns} runs of ${fmtNum(gens, 0)} generations, ${chains} chains each, ` +
      `${threads ? `across ${Math.min(nRuns, Pool.cores())} threads` : 'on the page (this browser does not allow threads)'}. ` +
      'The runs are independent: that is the only way to check afterwards whether they converged.'));
    const t0 = performance.now();
    const bar = window.LABG ? LABG.progressBar(prog, { label: T('Inferencia bayesiana', 'Bayesian inference') }) : null;
    if (bar) bar.update(null, T('arrancando…', 'starting…'));
    else prog.textContent = T('arrancando…', 'starting…');

    Pool.run({
      task: 'mcmcRun', n: nRuns, seed: 97, minPerWorker: 1, workers: nRuns,
      payload: {
        seqs: p.seqs, labels: p.taxa, type: p.type, spec,
        startNewick: Tree.writeNewick(start, { labels: p.taxa, support: false }),
        opts: {
          generations: gens, burnin: Math.floor(gens * burnPct / 100), sampleEvery: every,
          chains, deltaT, priors: priorsFromUI(),
        },
      },
      onProgress: (done, total) => {
        const txt = `${done}/${total} ${T('corridas', 'runs')} · ${((performance.now() - t0) / 1000).toFixed(0)} s`;
        if (bar) bar.update(done / total, txt);
        else prog.textContent = txt;
      },
      cancelled: () => B6.cancelled,
    }).then(out => {
      B6.running = false;
      btn.disabled = false; cancel.style.display = 'none';
      if (!bar || B6.cancelled) prog.textContent = '';
      if (!out.length) {
        if (bar && !B6.cancelled) bar.fail(T('ninguna corrida', 'no run'));
        showMessage(msg, 'warning', L2('No volvió ninguna corrida.', 'No run came back.'));
        return;
      }
      if (bar && !B6.cancelled) bar.done(T(`${out.length} corridas listas`, `${out.length} runs done`));
      /* the split frequencies travel as pairs; put them back into Maps */
      B6.runs = out.map(r => Object.assign({}, r, { splitFreq: new Map(r.splitFreq) }));
      B6.spec = spec;
      B6.seconds = (performance.now() - t0) / 1000;
      afterRun();
    });
  }

  function afterRun() {
    const p = current();
    const runs = B6.runs;
    const all = [].concat.apply([], runs.map(r => r.samples));
    const lnLs = all.map(s => s.lnL);
    statTiles('p6Tiles', [
      [T('Corridas', 'Runs'), runs.length, `${runs[0].chains} ${T('cadenas cada una', 'chains each')}`],
      [T('Generaciones', 'Generations'), fmtNum(runs[0].generations, 0), `${T('descarte', 'burn-in')} ${fmtNum(runs[0].burnin, 0)}`],
      [T('Muestras', 'Samples'), fmtNum(all.length, 0), T('después del descarte', 'after burn-in')],
      [T('lnL medio', 'mean lnL'), fmtLnL(Mcmc.mean(lnLs)), T('del posterior', 'of the posterior')],
      [T('Intercambios aceptados', 'Swaps accepted'), fmtPct(Mcmc.mean(runs.map(r => r.swapRate)), 1),
        runs[0].chains < 2 ? T('sin cadenas calientes', 'no heated chains')
          : Mcmc.mean(runs.map(r => r.swapRate)) < 0.01 ? T('casi ninguno: baja el calentamiento', 'almost none: lower the heating')
            : T('entre cadenas', 'between chains'),
        runs[0].chains < 2 ? '' : (Mcmc.mean(runs.map(r => r.swapRate)) < 0.01 ? 'bad' : 'ok')],
      [T('Tiempo', 'Time'), `${B6.seconds.toFixed(1)} s`, ''],
    ]);
    showMessage(el('p6Msg'), 'success', L2(
      `${runs.length} corridas terminadas en ${B6.seconds.toFixed(1)} s. Antes de mirar el árbol, mira los diagnósticos.`,
      `${runs.length} runs finished in ${B6.seconds.toFixed(1)} s. Before looking at the tree, look at the diagnostics.`));
    renderDiagnostics();
    summarise();
    ['p6Diag', 'p6Tree', 'p6Bayes', 'p6Export'].forEach(id => { el(id).style.display = ''; });
    commit();
  }

  /* ================================================================
     3 · diagnostics
     ================================================================ */
  function paramKeys() {
    if (!B6.runs || !B6.runs[0].samples.length) return [];
    return Object.keys(B6.runs[0].samples[0]).filter(k => k !== 'gen');
  }

  function renderDiagnostics() {
    const runs = B6.runs;
    const keys = paramKeys();
    const asdsf = runs.length > 1 ? Mcmc.asdsf(runs.map(r => r.splitFreq), 0.1) : null;
    const rows = keys.map(k => {
      const chains = runs.map(r => r.samples.map(s => s[k]));
      const pooled = [].concat.apply([], chains);
      const h = Mcmc.hpd(pooled, 0.95);
      return {
        param: k,
        mean: Mcmc.mean(pooled),
        sd: Mcmc.sd(pooled),
        lower: h.lower, upper: h.upper,
        ess: Math.min.apply(null, chains.map(c => Mcmc.ess(c))),
        psrf: runs.length > 1 ? Mcmc.psrf(chains) : null,
      };
    });
    B6.diag = { rows, asdsf };

    const minESS = Math.min.apply(null, rows.map(r => r.ess));
    const maxPSRF = runs.length > 1 ? Math.max.apply(null, rows.map(r => r.psrf || 1)) : null;
    statTiles('p6DiagTiles', [
      asdsf == null ? null : ['ASDSF', fmtFixed(asdsf, 5), asdsf < 0.01 ? T('por debajo de 0.01', 'below 0.01') : T('todavía alto', 'still high'),
        asdsf < 0.01 ? 'ok' : 'bad'],
      maxPSRF == null ? null : ['PSRF ' + T('máximo', 'largest'), fmtFixed(maxPSRF, 4), maxPSRF < 1.02 ? T('cerca de 1', 'close to 1') : T('lejos de 1', 'far from 1'),
        maxPSRF < 1.02 ? 'ok' : 'bad'],
      [T('ESS mínimo', 'smallest ESS'), fmtFixed(minESS, 0), minESS >= 200 ? T('suficiente', 'enough') : T('insuficiente', 'not enough'),
        minESS >= 200 ? 'ok' : 'bad'],
      [T('Parámetros vigilados', 'Parameters watched'), rows.length, ''],
    ].filter(Boolean));

    buildTable('p6DiagTable', [
      { key: 'param', label: T('Parámetro', 'Parameter') },
      { key: 'mean', label: T('media', 'mean'), num: true, fmt: v => fmtFixed(v, 4) },
      { key: 'sd', label: T('desv. típica', 'sd'), num: true, fmt: v => fmtFixed(v, 4) },
      { key: 'lower', label: 'HPD 95 % ' + T('inferior', 'lower'), num: true, fmt: v => fmtFixed(v, 4) },
      { key: 'upper', label: 'HPD 95 % ' + T('superior', 'upper'), num: true, fmt: v => fmtFixed(v, 4) },
      { key: 'ess', label: 'ESS', num: true, fmt: v => fmtFixed(v, 0) },
      { key: 'psrf', label: 'PSRF', num: true, fmt: v => (v == null ? '—' : fmtFixed(v, 4)) },
    ], rows, {});

    const bad = [];
    if (asdsf != null && asdsf >= 0.01) bad.push(L2('el ASDSF sigue por encima de 0.01', 'the ASDSF is still above 0.01'));
    if (maxPSRF != null && maxPSRF >= 1.02) bad.push(L2('algún PSRF se aleja de 1', 'some PSRF is away from 1'));
    if (minESS < 200) bad.push(L2('algún tamaño efectivo de muestra no llega a 200', 'some effective sample size falls short of 200'));
    const v = el('p6DiagVerdict');
    v.style.display = '';
    v.innerHTML = bad.length
      ? `<b>${L2('Todavía no', 'Not yet')}</b> ` + L2(
        `Hay señales de que las cadenas no han convergido: ${bad.join('; ')}. Alarga la corrida (multiplicar por diez las generaciones es lo normal), sube el número de cadenas calientes o baja el calentamiento. Las probabilidades posteriores de abajo no son de fiar mientras esto siga así.`,
        `There are signs the chains have not converged: ${bad.join('; ')}. Run it longer (ten times the generations is the usual step), raise the number of heated chains, or lower the heating. The posterior probabilities below cannot be trusted while this is the case.`)
      : `<b>${L2('Nada delata falta de convergencia', 'Nothing betrays a lack of convergence')}</b> ` + L2(
        'El ASDSF está por debajo de 0.01, los PSRF cerca de 1 y todos los tamaños efectivos de muestra pasan de 200. Eso no <i>demuestra</i> que la cadena haya visitado todo el posterior —ninguna cifra puede—, pero es lo que se pide para publicar.',
        'The ASDSF is below 0.01, the PSRFs near 1 and every effective sample size is above 200. That does not <i>prove</i> the chain visited the whole posterior — no number can — but it is what is asked for before publishing.');

    /* the traces */
    const sel = el('p6Param');
    const prev = sel.value;
    sel.innerHTML = keys.map(k => `<option value="${k}">${k}</option>`).join('');
    sel.value = keys.indexOf(prev) >= 0 ? prev : (keys.indexOf('alpha') >= 0 ? 'alpha' : keys[0]);
    el('p6TracePane').style.display = '';
    el('p6HistPane').style.display = '';
    drawTrace();
    drawHistogram();
  }

  function drawTrace() {
    const runs = B6.runs;
    if (!runs) return;
    const series = runs.map((r, i) => ({
      name: `${T('corrida', 'run')} ${i + 1}`,
      x: r.samples.map(s => s.gen),
      y: r.samples.map(s => s.lnL),
    }));
    el('p6TraceFig').innerHTML = Plots6.lines(series, {
      width: 760, height: 260,
      xLabel: T('generación', 'generation'), yLabel: 'ln L',
    });
  }
  function drawHistogram() {
    const runs = B6.runs;
    const k = el('p6Param').value;
    if (!runs || !k) return;
    const pooled = [].concat.apply([], runs.map(r => r.samples.map(s => s[k])));
    const h = Mcmc.hpd(pooled, 0.95);
    el('p6HistTitle').textContent = T(`Posterior de ${k}`, `Posterior of ${k}`);
    el('p6HistFig').innerHTML = Plots6.histogram(pooled, {
      width: 760, height: 240, bins: 40,
      xLabel: k, yLabel: T('muestras', 'samples'),
      marks: [{ x: h.lower, label: 'HPD' }, { x: h.upper, label: '' }, { x: Mcmc.mean(pooled), label: T('media', 'mean') }],
    });
  }

  /* ================================================================
     4 · the tree
     ================================================================ */
  function summarise() {
    const p = current();
    const runs = B6.runs;
    /* the posterior of every split, pooled over runs */
    const pooled = new Map();
    let nTot = 0;
    runs.forEach(r => {
      const n = r.trees.length;
      nTot += n;
      r.splitFreq.forEach((v, k) => pooled.set(k, (pooled.get(k) || 0) + v * n));
    });
    const freq = new Map();
    pooled.forEach((v, k) => freq.set(k, v / nTot));
    B6.postSplits = freq;

    /* the majority-rule consensus, built from the splits above one half */
    const chosen = [];
    freq.forEach((v, k) => { if (v >= 0.5) chosen.push({ key: k, set: k.split(',').map(Number), freq: v }); });
    const cons = Consensus.buildFromSplits(chosen, p.taxa.length, p.taxa);
    Pars.applySupport(cons, freq, p.taxa.length, false);
    /* The mean length of every branch, over the samples that contain it. An
       internal branch is named by its split; a terminal branch by its taxon,
       which no split names because a split needs at least two tips. Without the
       terminal ones the consensus would carry the length 1 that the star tree
       starts with, and the scale bar would be a lie. */
    const lenSum = new Map(), lenN = new Map();
    const add = (key, v) => { lenSum.set(key, (lenSum.get(key) || 0) + v); lenN.set(key, (lenN.get(key) || 0) + 1); };
    runs.forEach(r => r.trees.forEach(nw => {
      const t = Tree.parseNewick(nw, p.taxa);
      Tree.splits(t, p.taxa.length).forEach((node, key) => {
        const par = findParent(t, node);
        if (!par) return;
        const e = par.children.find(x => x.node === node);
        if (e) add(key, e.len || 0);
      });
      Tree.nodes(t).forEach(nd => nd.children.forEach(c => {
        if (c.node.tip != null) add('tip:' + c.node.tip, c.len || 0);
      }));
    }));
    const meanOf = key => (lenN.get(key) ? lenSum.get(key) / lenN.get(key) : null);
    Tree.splits(cons, p.taxa.length).forEach((node, key) => {
      const par = findParent(cons, node);
      if (!par) return;
      const e = par.children.find(x => x.node === node);
      const m = meanOf(key);
      if (e && m != null) e.len = m;
    });
    Tree.nodes(cons).forEach(nd => nd.children.forEach(c => {
      if (c.node.tip != null) {
        const m = meanOf('tip:' + c.node.tip);
        c.len = m == null ? 0 : m;
      }
    }));
    B6.consensus = cons;
    B6.consensusHasLengths = true;
    B6.meanLengths = { splits: lenSum, counts: lenN };
    /* the same for the most sampled topology, so the two are comparable */
    B6.tipMean = meanOf;

    /* the topology that was sampled most often */
    const counts = new Map();
    runs.forEach(r => r.trees.forEach(nw => {
      const k = Mcmc.topologyKey(nw, p.taxa, p.taxa.length);
      const e = counts.get(k);
      if (e) e.n++; else counts.set(k, { n: 1, newick: nw });
    }));
    const sorted = [...counts.values()].sort((a, b) => b.n - a.n);
    B6.map = sorted.length ? Tree.parseNewick(sorted[0].newick, p.taxa) : null;
    if (B6.map) Pars.applySupport(B6.map, freq, p.taxa.length, false);
    B6.mapP = sorted.length ? sorted[0].n / nTot : 0;
    B6.nTopologies = sorted.length;
    let acc = 0, credible = 0;
    for (const t of sorted) { acc += t.n / nTot; credible++; if (acc >= 0.95) break; }
    B6.credible = credible;

    renderTree();
    renderSplitTable();
  }
  function findParent(root, target) {
    let found = null;
    (function walk(n) { n.children.forEach(c => { if (c.node === target) found = n; else walk(c.node); }); })(root);
    return found;
  }

  function renderTree() {
    const p = current();
    const which = el('p6Which').value;
    const tree = which === 'map' ? B6.map : B6.consensus;
    if (!tree) return;
    const t = Tree.clone(tree);
    Tree.ladderize(t);
    el('p6TreeFig').innerHTML = TreeView.render(t, {
      labels: p.taxa, width: el('p6Layout').value === 'circular' ? 720 : 780,
      layout: el('p6Layout').value,
      showSupport: true, supportAs: 'node', supportDecimals: 2,
      cladogram: which === 'majority' && !B6.consensusHasLengths,
      scaleLabel: T(' sustituciones/sitio', ' substitutions/site'),
    });
    const msg = el('p6TreeMsg');
    clearMessages(msg);
    showMessage(msg, 'info', L2(
      `El muestreo visitó <b>${fmtNum(B6.nTopologies, 0)}</b> topologías distintas; la más frecuente se llevó ${fmtPct(B6.mapP, 1)} del posterior, y hacen falta <b>${fmtNum(B6.credible, 0)}</b> para juntar el 95 %. ` +
      (B6.nTopologies > 1 ? 'Cuando el conjunto creíble tiene muchas topologías, publicar una sola como «el árbol» es esconder el resultado.' : ''),
      `The sampler visited <b>${fmtNum(B6.nTopologies, 0)}</b> distinct topologies; the commonest took ${fmtPct(B6.mapP, 1)} of the posterior, and <b>${fmtNum(B6.credible, 0)}</b> are needed to make up 95 %. ` +
      (B6.nTopologies > 1 ? 'When the credible set holds many topologies, publishing one as "the" tree hides the result.' : '')));
  }

  function renderSplitTable() {
    const p = current();
    if (!B6.postSplits) return;
    const boot = state.ml && state.ml.support && state.ml.support.boot
      ? new Map(state.ml.support.boot.freq) : null;
    const rows = [];
    B6.postSplits.forEach((v, key) => {
      if (v < 0.05) return;
      const names = key.split(',').map(Number).map(i => p.taxa[i]);
      rows.push({
        clade: names.length <= 4 ? names.join(', ') : `${names.slice(0, 3).join(', ')} … (${names.length})`,
        size: names.length,
        pp: v,
        boot: boot ? (boot.get(key) || 0) * 100 : null,
      });
    });
    rows.sort((a, b) => b.pp - a.pp);
    buildTable('p6SplitTable', [
      { key: 'clade', label: T('Clado', 'Clade') },
      { key: 'size', label: T('taxones', 'taxa'), num: true },
      { key: 'pp', label: T('probabilidad posterior (≥0.95)', 'posterior probability (≥0.95)'), num: true, fmt: v => fmtFixed(v, 3) },
      { key: 'boot', label: T('bootstrap ML % (≥70)', 'ML bootstrap % (≥70)'), num: true, fmt: v => (v == null ? '—' : fmtFixed(v, 0)) },
    ], rows, { limit: 60 });
  }

  /* ================================================================
     5 · Bayes factors
     ================================================================ */
  function runBayesFactors() {
    const p = current();
    if (!p) return;
    const msg = el('p6BFMsg');
    clearMessages(msg);
    const picked = [...el('p6BFModels').selectedOptions].map(o => o.value);
    if (picked.length < 2) {
      showMessage(msg, 'info', L2('Elige al menos dos modelos.', 'Choose at least two models.'));
      return;
    }
    const btn = el('p6RunBF'), prog = el('p6BFProgress');
    btn.disabled = true;
    const steps = Math.max(4, +el('p6Steps').value || 16);
    const stepGens = Math.max(500, +el('p6StepGens').value || 6000);
    const start = state.ml && state.ml.tree && state.ml.partition === p.name ? state.ml.tree : startingTree();
    const startNewick = Tree.writeNewick(start, { labels: p.taxa, support: false });
    const t0 = performance.now();
    showMessage(msg, 'info', L2(
      `Estimando ${picked.length} verosimilitudes marginales con ${steps} escalones de ${fmtNum(stepGens, 0)} generaciones cada uno. ` +
      'Es la parte más cara de todo el bloque: cada modelo es una escalera entera de cadenas.',
      `Estimating ${picked.length} marginal likelihoods with ${steps} steps of ${fmtNum(stepGens, 0)} generations each. ` +
      'It is the most expensive thing in the block: each model is a whole ladder of chains.'));

    const jobs = picked.map(m => ({ name: m, spec: specFromUI(m) }));
    let done = 0;
    const results = [];
    const bar = window.LABG ? LABG.progressBar(prog, { label: T('Factores de Bayes', 'Bayes factors') }) : null;
    const runOne = i => {
      if (i >= jobs.length) {
        finishBF(results, performance.now() - t0); btn.disabled = false;
        if (!bar) prog.textContent = '';
        else if (results.length) bar.done(T(`${results.length} modelos listos`, `${results.length} models done`));
        else bar.fail(T('sin resultados', 'no results'));
        return;
      }
      if (bar) bar.update(done / jobs.length, `${done}/${jobs.length} · ${jobs[i].name}`);
      else prog.textContent = `${done}/${jobs.length} · ${jobs[i].name}`;
      Pool.run({
        task: 'marginalLikelihood', n: 1, seed: 31 + i, minPerWorker: 1, workers: 1,
        payload: {
          seqs: p.seqs, labels: p.taxa, type: p.type, spec: jobs[i].spec, name: jobs[i].name,
          startNewick,
          opts: { steps, stepGenerations: stepGens, priors: priorsFromUI() },
        },
      }).then(out => {
        if (out[0]) results.push(Object.assign({ name: jobs[i].name, spec: jobs[i].spec }, out[0]));
        done++;
        runOne(i + 1);
      });
    };
    runOne(0);
  }

  function finishBF(results, ms) {
    const msg = el('p6BFMsg');
    clearMessages(msg);
    if (!results.length) { showMessage(msg, 'warning', L2('No volvió ningún resultado.', 'Nothing came back.')); return; }
    results.sort((a, b) => b.logMarginalLikelihood - a.logMarginalLikelihood);
    const best = results[0];
    B6.bf = results;
    buildTable('p6BFTable', [
      { key: 'name', label: T('Modelo', 'Model') },
      { key: 'lnML', label: T('log P(datos | modelo)', 'log P(data | model)'), num: true, fmt: v => fmtFixed(v, 3) },
      { key: 'twoLnBF', label: '2 ln BF ' + T('contra el mejor', 'against the best'), num: true, fmt: v => fmtFixed(v, 2) },
      { key: 'verdict', label: T('evidencia', 'evidence') },
    ], results.map(r => {
      const two = 2 * (best.logMarginalLikelihood - r.logMarginalLikelihood);
      return {
        name: modelLabel(r.spec), lnML: r.logMarginalLikelihood, twoLnBF: two,
        verdict: r === best ? T('— el mejor —', '— the best —')
          : two > 10 ? T('muy fuerte en contra', 'very strong against')
            : two > 6 ? T('fuerte en contra', 'strong against')
              : two > 2 ? T('positiva en contra', 'positive against')
                : T('no distingue', 'cannot tell apart'),
      };
    }), {});
    showMessage(msg, 'success', L2(
      `Gana <b>${esc(modelLabel(best.spec))}</b>, en ${(ms / 1000).toFixed(1)} s. ` +
      'La verosimilitud marginal ya penaliza los parámetros de más por sí sola: no hay que añadirle nada, a diferencia del AIC.',
      `<b>${esc(modelLabel(best.spec))}</b> wins, in ${(ms / 1000).toFixed(1)} s. ` +
      'The marginal likelihood already penalises extra parameters on its own: nothing has to be added to it, unlike the AIC.'));
    commit();
  }

  /* ================================================================
     6 · export and hand-off
     ================================================================ */
  function exportAs(kind) {
    const p = current();
    if (!p || !B6.runs) return;
    const base = slug(p.name) + '_bayes';
    if (kind === 'log') {
      const all = [].concat.apply([], B6.runs.map((r, i) => r.samples.map(s => Object.assign({ run: i + 1 }, s))));
      download(Mcmc.tracerLog(all), base + '.log', 'text/plain');
    } else if (kind === 'trees') {
      download([].concat.apply([], B6.runs.map(r => r.trees)).join('\n') + '\n', base + '.t');
    } else if (kind === 'consensus' && B6.consensus) {
      download(Tree.writeNewick(B6.consensus, { labels: p.taxa, supportDecimals: 3 }), base + '_consenso.nwk');
    } else if (kind === 'nexus' && B6.consensus) {
      download(SeqIO.writeNexus(p.taxa, p.seqs, {
        type: p.type,
        trees: [{ name: 'consenso', newick: Tree.writeNewick(B6.consensus, { labels: p.taxa, supportDecimals: 3 }) }],
      }), base + '.nex');
    } else if (kind === 'svg') {
      download(el('p6TreeFig').innerHTML, base + '.svg', 'image/svg+xml');
    } else if (kind === 'splits' && B6.postSplits) {
      const lines = ['clado,taxones,probabilidad_posterior'];
      B6.postSplits.forEach((v, k) => {
        const names = k.split(',').map(Number).map(i => p.taxa[i]);
        lines.push(`"${names.join(' ')}",${names.length},${v.toFixed(6)}`);
      });
      download(lines.join('\n') + '\n', base + '_clados.csv', 'text/csv');
    }
  }

  function commit() {
    const p = current();
    if (!p || !B6.runs) return;
    state.bayes = {
      partition: p.name,
      spec: B6.spec,
      runs: B6.runs.length,
      generations: B6.runs[0].generations,
      burnin: B6.runs[0].burnin,
      chains: B6.runs[0].chains,
      priors: priorsFromUI(),
      asdsf: B6.diag ? B6.diag.asdsf : null,
      diagnostics: B6.diag ? B6.diag.rows : null,
      splits: [...B6.postSplits.entries()],
      consensus: B6.consensus,
      consensusNewick: Tree.writeNewick(B6.consensus, { labels: p.taxa, supportDecimals: 3 }),
      mapNewick: B6.map ? Tree.writeNewick(B6.map, { labels: p.taxa, supportDecimals: 3 }) : null,
      mapProbability: B6.mapP,
      credibleSet: B6.credible,
      nTopologies: B6.nTopologies,
      bayesFactors: B6.bf ? B6.bf.map(r => ({ name: modelLabel(r.spec), lnML: r.logMarginalLikelihood })) : null,
    };
    enableStep(7, true);
  }

  /* ================================================================
     wiring
     ================================================================ */
  function refresh() {
    const has = parts().length > 0;
    el('p6NoData').style.display = has ? 'none' : '';
    ['p6Priors', 'p6Run'].forEach(id => { el(id).style.display = has ? '' : 'none'; });
    if (!has) { ['p6Diag', 'p6Tree', 'p6Bayes', 'p6Export'].forEach(id => { el(id).style.display = 'none'; }); return; }
    const ps = parts();
    const sel = el('p6Part');
    sel.innerHTML = ps.map((q, i) => `<option value="${i}">${esc(q.name)} · ${q.taxa.length} × ${q.length}</option>`).join('');
    if (B6.part < 0 || B6.part >= ps.length) B6.part = 0;
    sel.value = String(B6.part);
    fillModels();
    if (state.ml && state.ml.spec && !B6.runs) useBlock5();
    renderPriorNote();
  }

  function init() {
    if (!el('panel-6')) return;
    el('p6Part').addEventListener('change', () => {
      B6.part = +el('p6Part').value;
      B6.A = null; B6.runs = null; B6.summary = null; B6.bf = null;
      ['p6Diag', 'p6Tree', 'p6Bayes', 'p6Export'].forEach(id => { el(id).style.display = 'none'; });
      refresh();
    });
    el('p6FromBlock5').addEventListener('click', useBlock5);
    ['p6BranchMean', 'p6AlphaMean', 'p6RateMean'].forEach(id => el(id).addEventListener('input', renderPriorNote));
    el('p6Go').addEventListener('click', go);
    el('p6Cancel').addEventListener('click', () => { B6.cancelled = true; });
    el('p6Param').addEventListener('change', drawHistogram);
    el('p6Which').addEventListener('change', renderTree);
    el('p6Layout').addEventListener('change', renderTree);
    el('p6RunBF').addEventListener('click', runBayesFactors);
    els('[data-p6export]').forEach(b => b.addEventListener('click', () => exportAs(b.dataset.p6export)));
    el('p6ToBlock7').addEventListener('click', () => {
      const b = document.querySelector('.step-btn[data-step="7"]');
      if (b && !b.disabled) goStep(7);
      else showMessage(el('p6TreeMsg'), 'info', L2(
        'El Bloque 7 (reloj molecular y tiempos) llega en la etapa siguiente.',
        'Block 7 (molecular clock and dating) arrives in the next stage.'));
    });
    document.addEventListener('stepchange', e => { if (e.detail.step === 6) refresh(); });
    document.addEventListener('langchange', () => { if (parts().length) { refresh(); if (B6.runs) { renderDiagnostics(); renderSplitTable(); renderTree(); } } });
    document.addEventListener('themechange', () => { if (B6.runs) { drawTrace(); drawHistogram(); renderTree(); } });
  }
  document.addEventListener('DOMContentLoaded', init);

  Object.assign(B6, { refresh, go, runBayesFactors, renderTree, renderDiagnostics, commit, current });
})();
