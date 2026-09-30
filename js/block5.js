/* PhylogenyPro — molecular phylogenetics in the browser.
   Copyright (C) 2026  Luis Ángel Barrera-Guzmán

   This program is free software: you can redistribute it and/or modify it under
   the terms of the GNU General Public License as published by the Free Software
   Foundation, either version 3 of the License, or (at your option) any later
   version. It is distributed in the hope that it will be useful, but WITHOUT ANY
   WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
   A PARTICULAR PURPOSE. See the GNU General Public License, in the file LICENSE
   at the root of this program, or <https://www.gnu.org/licenses/>. */

/* PhylogenyPro — Block 5: maximum likelihood.

   The search itself lives in js/mlsearch.js and the likelihood in js/like.js;
   this file is the part of the block the user touches: choosing the model and
   the starting tree, watching the search, reading the four kinds of support
   next to each other, testing a constrained hypothesis, and leaving in
   state.ml everything Blocks 6 to 12 will start from.

   Two things here are deliberate. The bootstrap goes through js/pool.js, so a
   hundred replicates use every core the machine has instead of freezing the
   page; and the model of every replicate is held at the values fitted to the
   real data, as RAxML and IQ-TREE do, because a replicate is there to
   re-estimate the tree, not the model. */

(function () {

  const B5 = {
    part: -1, A: null, res: null, part5: null,
    support: { boot: null, uf: null, alrt: null },
    trees: {}, constrained: null, tests: null,
    cancelled: false, running: false, bootTrees: null,
  };
  window.B5 = B5;

  function parts() { return (state.data && state.data.parts) || []; }
  function current() { return parts()[B5.part] || null; }

  function prepare() {
    const p = current();
    if (!p) return null;
    if (!B5.A || B5.A.__part !== p.name) {
      B5.A = Like.compress(p.seqs, p.type);
      B5.A.__part = p.name;
    }
    return B5.A;
  }

  /* ================================================================
     the model the user asked for
     ================================================================ */
  const DNA_MODELS = ['JC', 'F81', 'K80', 'HKY', 'TrNef', 'TrN', 'K81', 'K81uf', 'TIMef', 'TIM', 'TVMef', 'TVM', 'SYM', 'GTR'];
  const AA_MODELS = ['LG', 'WAG', 'JTT', 'Dayhoff', 'cpREV', 'mtREV24', 'VT', 'Blosum62', 'MtMam', 'RtREV', 'HIVb', 'HIVw'];

  function fillModels() {
    const p = current();
    const sel = el('p5Model');
    if (!sel || !p) return;
    const prev = sel.value;
    let list;
    if (p.type === 'aa') list = AA_MODELS;
    else if (p.type === 'morph') list = ['Mk', 'Mkv'];
    else list = DNA_MODELS;
    sel.innerHTML = list.map(m => `<option value="${m}">${m}</option>`).join('');
    sel.value = list.indexOf(prev) >= 0 ? prev : (p.type === 'aa' ? 'LG' : p.type === 'morph' ? 'Mkv' : 'GTR');
    /* the frequency choice means nothing for the models that fix them */
    el('p5Freq').disabled = p.type === 'morph';
  }

  function specFromUI() {
    const p = current();
    const A = prepare();
    const model = el('p5Model').value;
    const spec = { type: p.type === 'morph' ? 'morph' : p.type, model, ncat: 1 };
    if (p.type === 'dna') spec.rates = [1, 2, 1, 1, 2, 1];
    if (el('p5Freq').value === 'empirical' && p.type !== 'morph') spec.freqs = Array.from(A.freqs);
    if (el('p5Gamma').checked) { spec.alpha = 0.5; spec.ncat = Math.max(2, +el('p5Cats').value || 4); }
    if (el('p5Inv').checked) spec.pInv = 0.2;
    if (model === 'Mkv') spec.ascertainment = 'variable';
    return spec;
  }

  /* the winner of Block 3, when there is one for this partition */
  function block3Spec() {
    const m = state.models;
    if (!m || !m.best || !m.best.spec) return null;
    if (m.partition && current() && m.partition !== current().name) return null;
    return m.best.spec;
  }
  function useBlock3() {
    const s = block3Spec();
    const msg = el('p5Msg');
    clearMessages(msg);
    if (!s) {
      showMessage(msg, 'info', L2(
        'Todavía no hay un modelo elegido en el Bloque 3 para esta partición.',
        'There is no model chosen in Block 3 for this partition yet.'));
      return;
    }
    const sel = el('p5Model');
    if ([...sel.options].some(o => o.value === s.model)) sel.value = s.model;
    el('p5Gamma').checked = !!s.alpha;
    if (s.ncat) el('p5Cats').value = s.ncat;
    el('p5Inv').checked = !!(s.pInv > 0);
    el('p5Freq').value = s.freqs ? 'empirical' : 'equal';
    showMessage(msg, 'success', L2(
      `Modelo del Bloque 3: <b>${esc(modelLabel(s))}</b>.`, `Model from Block 3: <b>${esc(modelLabel(s))}</b>.`));
  }
  function modelLabel(s) {
    return s.model + (s.alpha ? '+G' : '') + (s.pInv > 0 ? '+I' : '') + (s.__optF ? '+F' : '');
  }

  /* ================================================================
     the starting tree
     ================================================================ */
  function startingTree() {
    const p = current();
    const which = el('p5Start').value;
    const q = state.quick || {};
    if (which === 'parsimony' && q.parsimony && q.parsimony.tree) return Tree.clone(q.parsimony.tree);
    if (which === 'distance' && q.distance && q.distance.tree) return Tree.clone(q.distance.tree);
    if (which === 'random') {
      const order = Array.from({ length: p.taxa.length }, (_, i) => i);
      const r = rng(Date.now() & 0xffff);
      for (let i = order.length - 1; i > 0; i--) { const j = Math.floor(r() * (i + 1)); const t = order[i]; order[i] = order[j]; order[j] = t; }
      return randomTree(order, p.taxa);
    }
    const D = Dist.matrix(p.seqs, p.type === 'aa' ? 'poisson' : p.type === 'morph' ? 'p' : 'jc', { type: p.type }).D;
    return Tree.bionj(D, p.taxa);
  }
  /* a random binary tree: the honest starting point when you want to know
     whether the search depends on where it starts */
  function randomTree(order, labels) {
    let tree = { children: order.slice(0, 3).map(i => ({ node: { tip: i, label: labels[i], children: [] }, len: 0.05 })) };
    for (let k = 3; k < order.length; k++) {
      const edges = Tree.allEdges(tree);
      const at = edges[Math.floor(Math.random() * edges.length)];
      const leaf = { tip: order[k], label: labels[order[k]], children: [] };
      at.edge.node = { children: [{ node: at.edge.node, len: (at.edge.len || 0.05) / 2 }, { node: leaf, len: 0.05 }] };
    }
    return tree;
  }

  /* ================================================================
     1 · the search
     ================================================================ */
  function runSearch() {
    const p = current();
    if (!p) return;
    const A = prepare();
    const msg = el('p5Msg');
    clearMessages(msg);
    const btn = el('p5Run'), cancel = el('p5Cancel'), prog = el('p5Progress');
    btn.disabled = true; cancel.style.display = ''; B5.cancelled = false; B5.running = true;
    prog.textContent = T('preparando…', 'preparing…');

    const spec = specFromUI();
    const start = startingTree();
    const t0 = performance.now();
    const w = window.LABG ? LABG.work({ title: T('Buscando el árbol de máxima verosimilitud', 'Searching for the maximum-likelihood tree'), delay: 300 }) : null;
    phyAfterPaint(() => {
      let res;
      try {
        res = ML.search(A, spec, {
          start, labels: p.taxa,
          collect: true, maxCollect: 400,
          spr: el('p5Swap').value === 'spr',
          radius: +el('p5Radius').value || 5,
          optFreqs: el('p5OptF').checked,
          maxRounds: +el('p5Rounds').value || 25,
          onProgress: (round, best, evaluated) => {
            prog.textContent = `${T('ronda', 'round')} ${round} · lnL ${fmtLnL(best)} · ${fmtNum(evaluated, 0)} ${T('topologías', 'topologies')}`;
            if (w) w.message(prog.textContent);
          },
          cancelled: () => B5.cancelled,
        });
      } catch (e) {
        showMessage(msg, 'error', esc(e.message));
        btn.disabled = false; cancel.style.display = 'none'; prog.textContent = ''; B5.running = false;
        if (w) w.close();
        return;
      }
      B5.res = res;
      res.spec.__optF = el('p5OptF').checked;
      B5.trees.ml = res.tree;
      B5.trees.mlName = T(`Máxima verosimilitud · ${modelLabel(res.spec)}`, `Maximum likelihood · ${modelLabel(res.spec)}`);
      B5.support = { boot: null, uf: null, alrt: null };
      B5.bootTrees = null;
      showResult(res, performance.now() - t0);
      btn.disabled = false; cancel.style.display = 'none'; prog.textContent = '';
      B5.running = false;
      commit();
      if (w) w.done();
    }, w);
  }

  function showResult(res, ms) {
    const p = current();
    const s = res.spec;
    const tiles = [
      [T('log-verosimilitud', 'log-likelihood'), fmtLnL(res.lnL), modelLabel(s), 'ok'],
      [T('Parámetros', 'Parameters'), res.k, T('ramas y modelo', 'branches and model')],
      ['AIC', fmtFixed(res.AIC, 2), ''],
      ['BIC', fmtFixed(res.BIC, 2), ''],
      [T('Rondas', 'Rounds'), res.rounds, `${fmtNum(res.evaluated, 0)} ${T('topologías', 'topologies')}`],
      [T('Tiempo', 'Time'), `${(ms / 1000).toFixed(1)} s`, ''],
    ];
    if (s.alpha) tiles.push(['α', fmtFixed(s.alpha, 4), T('forma de la gamma', 'gamma shape')]);
    if (s.pInv > 0) tiles.push(['p(inv)', fmtFixed(s.pInv, 4), T('sitios invariables', 'invariable sites')]);
    if (s.type === 'dna' && s.rates && s.rates[1] != null) {
      tiles.push(['κ / r(A↔G)', fmtFixed(s.rates[1], 3), T('respecto a A↔C', 'relative to A↔C')]);
    }
    statTiles('p5Tiles', tiles);

    /* what the numbers mean, in words */
    const bits = [];
    bits.push(L2(
      `La verosimilitud del árbol es <b>${fmtLnL(res.lnL)}</b>: ese número solo tiene sentido comparado con el de otro árbol o de otro modelo sobre <i>los mismos datos</i>, nunca en solitario.`,
      `The likelihood of the tree is <b>${fmtLnL(res.lnL)}</b>: that number only means something compared with another tree or another model on <i>the same data</i>, never on its own.`));
    if (s.alpha != null) {
      const a = s.alpha;
      bits.push(a < 0.5 ? L2(
        `Una α de ${fmtFixed(a, 3)} es muy baja: la variación de velocidad entre sitios es extrema, con unos pocos sitios cambiando mucho y la mayoría casi congelados. Ignorarla habría acortado las ramas largas.`,
        `An α of ${fmtFixed(a, 3)} is very low: rate variation among sites is extreme, with a few sites changing a lot and most nearly frozen. Ignoring it would have shortened the long branches.`)
        : L2(`Una α de ${fmtFixed(a, 3)} indica una variación de velocidad moderada entre sitios.`,
          `An α of ${fmtFixed(a, 3)} means moderate rate variation among sites.`));
    }
    const start = el('p5Start').value;
    bits.push(L2(
      `La búsqueda se detuvo cuando ningún reordenamiento mejoraba: eso es un óptimo <i>local</i>. Vuelve a correrla desde un árbol al azar (arriba, «árbol inicial») y compara: si llega al mismo sitio, el resultado es sólido.`,
      `The search stopped when no rearrangement improved it: that is a <i>local</i> optimum. Run it again from a random tree ("starting tree" above) and compare: if it lands in the same place, the result is solid.`));
    if (start !== 'random') bits.push(L2(
      'Ahora mismo partió de un árbol razonable, no al azar.', 'Right now it started from a reasonable tree, not a random one.'));
    const v = el('p5Verdict');
    v.style.display = '';
    v.innerHTML = `<b>${L2('Cómo se lee', 'How to read it')}</b> ` + bits.join(' ');

    el('p5Pane').style.display = '';
    el('p5Title').textContent = B5.trees.mlName;
    drawInto('p5Fig', res.tree, p.taxa);
    showMessage(el('p5Msg'), 'success', L2(
      `Árbol encontrado en ${res.rounds} rondas tras evaluar ${fmtNum(res.evaluated, 0)} topologías.`,
      `Tree found in ${res.rounds} rounds after evaluating ${fmtNum(res.evaluated, 0)} topologies.`));
    refreshTreeList();
    fillConstraintList();
  }

  /* ================================================================
     2 · one model per partition
     ================================================================ */
  function runPartitioned() {
    const ps = parts();
    const msg = el('p5PartMsg');
    clearMessages(msg);
    if (ps.length < 2) {
      showMessage(msg, 'info', L2(
        'Solo hay una partición: no hay nada que particionar. Vuelve al Bloque 2 para concatenar varios marcadores o dividir por posición del codón.',
        'There is only one partition: there is nothing to partition. Go back to Block 2 to concatenate several markers or split by codon position.'));
      return;
    }
    const btn = el('p5RunPart'), prog = el('p5PartProgress');
    btn.disabled = true;
    prog.textContent = T('ajustando…', 'fitting…');
    const w = window.LABG ? LABG.work({ title: T('Ajustando un modelo por partición', 'Fitting one model per partition'), delay: 300 }) : null;
    phyAfterPaint(() => {
      try {
        /* the tree they share: the ML tree if there is one, otherwise BIONJ on
           the concatenation */
        let tree = B5.trees.ml ? Tree.clone(B5.trees.ml) : null;
        const taxa = ps[0].taxa;
        if (!tree) {
          const cat = taxa.map((_, i) => ps.map(p => p.seqs[i]).join(''));
          tree = Tree.bionj(Dist.matrix(cat, 'jc', { type: 'dna' }).D, taxa);
        }
        const chosen = ps.map(p => {
          const A = Like.compress(p.seqs, p.type);
          const s = (state.models && state.models.partition === p.name && state.models.best) ? state.models.best.spec : null;
          const spec = s ? Like.cloneSpec(s) : Object.assign(specFromUI(), { freqs: Array.from(A.freqs) });
          return { A, spec, name: p.name };
        });
        const t0 = performance.now();
        const res = Like.fitPartitioned(tree, chosen, { passes: 12 });
        B5.part5 = res;
        B5.trees.partitioned = res.tree;
        B5.trees.partitionedName = T(`Particionado · ${ps.length} particiones`, `Partitioned · ${ps.length} partitions`);

        /* the same tree under one model for everybody, for comparison */
        const catA = Like.compress(taxa.map((_, i) => ps.map(p => p.seqs[i]).join('')), ps[0].type);
        const single = Like.fit(tree, catA, Object.assign(specFromUI(), { freqs: Array.from(catA.freqs) }), { passes: 12 });

        statTiles('p5PartTiles', [
          [T('log-verosimilitud', 'log-likelihood'), fmtLnL(res.lnL), T('particionado', 'partitioned'), 'ok'],
          [T('un solo modelo', 'single model'), fmtLnL(single.lnL), '', ''],
          [T('Parámetros', 'Parameters'), `${res.k} vs ${single.k}`, T('particionado vs único', 'partitioned vs single')],
          ['ΔBIC', fmtFixed(single.BIC - res.BIC, 2), single.BIC - res.BIC > 0 ? T('gana el particionado', 'partitioned wins') : T('gana el modelo único', 'single model wins'),
            single.BIC - res.BIC > 0 ? 'ok' : 'bad'],
          [T('Tiempo', 'Time'), `${((performance.now() - t0) / 1000).toFixed(1)} s`, ''],
        ]);

        buildTable('p5PartTable', [
          { key: 'name', label: T('Partición', 'Partition') },
          { key: 'sites', label: T('sitios', 'sites'), num: true },
          { key: 'model', label: T('modelo', 'model') },
          { key: 'rate', label: T('velocidad relativa', 'relative rate'), num: true, fmt: v => fmtFixed(v, 3) },
          { key: 'alpha', label: 'α', num: true, fmt: v => (v == null ? '—' : fmtFixed(v, 3)) },
          { key: 'lnL', label: 'lnL', num: true, fmt: v => fmtLnL(v) },
        ], res.partitions.map(q => ({
          name: q.name, sites: q.nSites, model: modelLabel(q.spec),
          rate: q.rate, alpha: q.spec.alpha == null ? null : q.spec.alpha, lnL: q.lnL,
        })), {});

        const fastest = res.partitions.reduce((a, b) => (b.rate > a.rate ? b : a));
        const slowest = res.partitions.reduce((a, b) => (b.rate < a.rate ? b : a));
        const vb = el('p5PartVerdict');
        vb.style.display = '';
        vb.innerHTML = `<b>${L2('Cómo se lee', 'How to read it')}</b> ` + L2(
          `Las velocidades relativas están normalizadas a una media de uno: <b>${esc(fastest.name)}</b> evoluciona ${fmtFixed(fastest.rate / slowest.rate, 1)} veces más rápido que <b>${esc(slowest.name)}</b>. ` +
          (single.BIC - res.BIC > 0
            ? `El BIC baja ${fmtFixed(single.BIC - res.BIC, 1)} unidades al dar a cada partición su modelo: el reparto vale lo que cuesta en parámetros.`
            : `El BIC <i>no</i> mejora al particionar: con estos datos, un solo modelo basta y el reparto solo añade parámetros.`),
          `The relative rates are normalised to a mean of one: <b>${esc(fastest.name)}</b> evolves ${fmtFixed(fastest.rate / slowest.rate, 1)} times faster than <b>${esc(slowest.name)}</b>. ` +
          (single.BIC - res.BIC > 0
            ? `The BIC drops by ${fmtFixed(single.BIC - res.BIC, 1)} when each partition gets its own model: the split is worth what it costs in parameters.`
            : `The BIC does <i>not</i> improve with partitioning: for these data one model is enough and the split only adds parameters.`));

        showMessage(msg, 'success', L2(
          `Modelo particionado ajustado sobre ${ps.length} particiones y ${fmtNum(res.nSites, 0)} sitios.`,
          `Partitioned model fitted over ${ps.length} partitions and ${fmtNum(res.nSites, 0)} sites.`));
        refreshTreeList();
        commit();
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
     3 · support
     ================================================================ */
  function runBootstrap() {
    const p = current();
    if (!p || !B5.res) { needTree('p5BootMsg'); return; }
    const msg = el('p5BootMsg');
    clearMessages(msg);
    const btn = el('p5RunBoot'), cancel = el('p5CancelBoot'), prog = el('p5BootProgress');
    const reps = +el('p5Reps').value || 100;
    btn.disabled = true; cancel.style.display = ''; B5.cancelled = false;
    const t0 = performance.now();
    const threads = Pool.available();
    showMessage(msg, 'info', L2(
      `Corriendo ${reps} réplicas ${threads ? `en ${Pool.cores()} hilos` : 'en la página (este navegador no permite hilos)'}. ` +
      'El modelo se mantiene fijo en los valores ajustados a los datos reales, como hacen RAxML e IQ-TREE.',
      `Running ${reps} replicates ${threads ? `across ${Pool.cores()} threads` : 'on the page (this browser does not allow threads)'}. ` +
      'The model is held at the values fitted to the real data, as RAxML and IQ-TREE do.'));
    const bar = window.LABG ? LABG.progressBar(prog, { label: T('Bootstrap de máxima verosimilitud', 'Maximum-likelihood bootstrap') }) : null;
    if (bar) bar.update(0, `0/${reps}`);

    Pool.run({
      task: 'mlBootstrap', n: reps, seed: 7,
      payload: {
        seqs: p.seqs, labels: p.taxa, type: p.type,
        spec: Like.cloneSpec(B5.res.spec),
        search: { maxRounds: 10, spr: false },
      },
      onProgress: (done, total) => {
        /* With W replicates running at once, the first W finish almost together,
           so "elapsed per replicate" is meaningless until a full round is done.
           What matters is how many ROUNDS of W are left. */
        const w = Pool.available() ? Pool.cores() : 1;
        const roundsDone = Math.floor(done / w);
        const left = roundsDone >= 1
          ? (performance.now() - t0) / roundsDone * (Math.ceil(total / w) - roundsDone) / 1000 : null;
        const txt = `${done}/${total}` + (left != null && left > 2
          ? ` · ${T('faltan', 'about')} ${left < 60 ? Math.ceil(left) + ' s' : Math.ceil(left / 60) + ' min'}` : '');
        if (bar) bar.update(done / total, txt);
        else prog.textContent = txt;
      },
      cancelled: () => B5.cancelled,
    }).then(reps2 => {
      const counts = new Map();
      reps2.forEach(r => r.splits.forEach(k => counts.set(k, (counts.get(k) || 0) + 1)));
      const freq = new Map();
      counts.forEach((v, k) => freq.set(k, v / reps2.length));
      B5.support.boot = { freq, reps: reps2.length, ms: performance.now() - t0 };
      B5.bootTrees = reps2.map(r => r.newick);
      renderSupportTable();
      redraw();
      const strong = [...freq.values()].filter(v => v >= 0.7).length;
      const nCl = Tree.splits(B5.trees.ml, p.taxa.length).size;
      clearMessages(msg);
      showMessage(msg, 'success', L2(
        `${reps2.length} réplicas en ${((performance.now() - t0) / 1000).toFixed(1)} s. ${strong} de los ${nCl} clados del árbol llegan al 70 %.`,
        `${reps2.length} replicates in ${((performance.now() - t0) / 1000).toFixed(1)} s. ${strong} of the tree's ${nCl} clades reach 70 %.`));
      btn.disabled = false; cancel.style.display = 'none';
      if (bar && !B5.cancelled) bar.done(T(`${reps2.length} réplicas listas`, `${reps2.length} replicates done`));
      else prog.textContent = '';
      commit();
    });
  }

  function runUFBoot() {
    const p = current();
    if (!p || !B5.res) { needTree('p5BootMsg'); return; }
    const msg = el('p5BootMsg');
    clearMessages(msg);
    if (!B5.res.visited || B5.res.visited.length < 5) {
      showMessage(msg, 'warning', L2(
        'El bootstrap ultrarrápido remuestrea los árboles que la búsqueda visitó, y esta guardó muy pocos. Vuelve a correr la búsqueda.',
        'The ultrafast bootstrap resamples the trees the search visited, and this one kept very few. Run the search again.'));
      return;
    }
    const t0 = performance.now();
    const uf = ML.ufboot(prepare(), B5.res.visited, B5.res.model, { reps: 1000, seed: 11 });
    B5.support.uf = { freq: uf.freq, reps: uf.reps, nTrees: uf.nTrees };
    renderSupportTable();
    redraw();
    showMessage(msg, 'success', L2(
      `${uf.reps} réplicas RELL sobre ${uf.nTrees} topologías en ${((performance.now() - t0) / 1000).toFixed(1)} s. Ojo: estos valores se leen a partir del <b>95 %</b>, no del 70 %.`,
      `${uf.reps} RELL replicates over ${uf.nTrees} topologies in ${((performance.now() - t0) / 1000).toFixed(1)} s. Careful: these are read from <b>95 %</b> up, not 70 %.`));
    commit();
  }

  function runALRT() {
    const p = current();
    if (!p || !B5.res) { needTree('p5BootMsg'); return; }
    const msg = el('p5BootMsg');
    clearMessages(msg);
    const btn = el('p5RunALRT');
    btn.disabled = true;
    const w = window.LABG ? LABG.work({ title: T('Evaluando las ramas (SH-aLRT y aBayes)', 'Testing the branches (SH-aLRT and aBayes)'), delay: 300 }) : null;
    phyAfterPaint(() => {
      const t0 = performance.now();
      const bt = ML.branchTests(B5.trees.ml, prepare(), B5.res.model, {});
      const alrt = new Map(), abayes = new Map();
      bt.tests.forEach(t => { alrt.set(t.split, t.shAlrt); abayes.set(t.split, t.aBayes); });
      B5.support.alrt = { alrt, abayes, ms: performance.now() - t0 };
      renderSupportTable();
      redraw();
      showMessage(msg, 'success', L2(
        `${bt.tests.length} ramas internas evaluadas en ${((performance.now() - t0) / 1000).toFixed(1)} s. SH-aLRT se lee a partir de 80; aBayes, a partir de 0.95.`,
        `${bt.tests.length} internal branches evaluated in ${((performance.now() - t0) / 1000).toFixed(1)} s. SH-aLRT is read from 80 up; aBayes from 0.95 up.`));
      btn.disabled = false;
      commit();
      if (w) w.done();
    }, w);
  }

  function needTree(id) {
    showMessage(el(id), 'info', L2('Primero busca el árbol de máxima verosimilitud.', 'Search for the maximum-likelihood tree first.'));
  }

  function renderSupportTable() {
    const p = current();
    if (!p || !B5.trees.ml) return;
    const n = p.taxa.length;
    const sp = Tree.splits(B5.trees.ml, n);
    const rows = [];
    sp.forEach((node, key) => {
      const names = key.split(',').map(Number).map(i => p.taxa[i]);
      rows.push({
        clade: names.length <= 4 ? names.join(', ') : `${names.slice(0, 3).join(', ')} … (${names.length})`,
        size: names.length,
        boot: B5.support.boot ? (B5.support.boot.freq.get(key) || 0) * 100 : null,
        uf: B5.support.uf ? (B5.support.uf.freq.get(key) || 0) * 100 : null,
        alrt: B5.support.alrt ? B5.support.alrt.alrt.get(key) : null,
        abayes: B5.support.alrt ? B5.support.alrt.abayes.get(key) : null,
      });
    });
    rows.sort((a, b) => (b.boot == null ? (b.uf == null ? -1 : b.uf) : b.boot) - (a.boot == null ? (a.uf == null ? -1 : a.uf) : a.boot));
    buildTable('p5SupportTable', [
      { key: 'clade', label: T('Clado', 'Clade') },
      { key: 'size', label: T('taxones', 'taxa'), num: true },
      { key: 'boot', label: T('bootstrap % (≥70)', 'bootstrap % (≥70)'), num: true, fmt: v => (v == null ? '—' : fmtFixed(v, 0)) },
      { key: 'uf', label: T('ultrarrápido % (≥95)', 'ultrafast % (≥95)'), num: true, fmt: v => (v == null ? '—' : fmtFixed(v, 0)) },
      { key: 'alrt', label: T('SH-aLRT (≥80)', 'SH-aLRT (≥80)'), num: true, fmt: v => (v == null ? '—' : fmtFixed(v, 1)) },
      { key: 'abayes', label: T('aBayes (≥0.95)', 'aBayes (≥0.95)'), num: true, fmt: v => (v == null ? '—' : fmtFixed(v, 3)) },
    ], rows, { limit: 60 });
  }

  /* ================================================================
     4 · constrained trees and topology tests
     ================================================================ */
  function fillConstraintList() {
    const p = current();
    const sel = el('p5Constraint');
    if (!sel || !p) return;
    sel.innerHTML = p.taxa.map((t, i) => `<option value="${i}">${esc(t)}</option>`).join('');
    const og = el('p5Outgroup');
    if (og) og.innerHTML = sel.innerHTML;
  }

  function runConstrained() {
    const p = current();
    if (!p || !B5.res) { needTree('p5TestMsg'); return; }
    const msg = el('p5TestMsg');
    clearMessages(msg);
    const picked = [...el('p5Constraint').selectedOptions].map(o => +o.value);
    if (picked.length < 2 || picked.length >= p.taxa.length - 1) {
      showMessage(msg, 'warning', L2(
        'Elige entre dos taxones y todos menos dos: un grupo de un solo taxón, o de todos, no es una hipótesis que se pueda poner a prueba.',
        'Choose between two taxa and all but two: a group of one taxon, or of all of them, is not a hypothesis that can be tested.'));
      return;
    }
    const btn = el('p5RunConstraint'), prog = el('p5ConProgress');
    btn.disabled = true; prog.textContent = T('buscando…', 'searching…');
    const w = window.LABG ? LABG.work({ title: T('Buscando el árbol restringido', 'Searching for the constrained tree'), delay: 300 }) : null;
    phyAfterPaint(() => {
      try {
        const t0 = performance.now();
        const con = ML.searchConstrained(prepare(), Like.cloneSpec(B5.res.spec), picked, {
          start: B5.trees.ml, labels: p.taxa, seqs: p.seqs, type: p.type, maxRounds: 15,
        });
        B5.constrained = { res: con, taxa: picked.map(i => p.taxa[i]) };
        B5.trees.constrained = con.tree;
        B5.trees.constrainedName = T(`Restringido · ${picked.length} taxones monofiléticos`, `Constrained · ${picked.length} taxa monophyletic`);
        const d = B5.res.lnL - con.lnL;
        showMessage(msg, 'success', L2(
          `El mejor árbol con ese grupo monofilético tiene lnL ${fmtLnL(con.lnL)}, frente a ${fmtLnL(B5.res.lnL)} sin restricción: cuesta <b>${fmtFixed(d, 3)}</b> unidades. ` +
          'Eso todavía no dice si la diferencia es significativa: para eso está la prueba de topología de abajo.',
          `The best tree with that group monophyletic has lnL ${fmtLnL(con.lnL)}, against ${fmtLnL(B5.res.lnL)} unconstrained: it costs <b>${fmtFixed(d, 3)}</b> units. ` +
          'That alone does not say whether the difference is significant: the topology test below is what does.'));
        el('p5ConProgress').textContent = `${((performance.now() - t0) / 1000).toFixed(1)} s`;
        refreshTreeList();
        commit();
        if (w) w.done();
      } catch (e) {
        showMessage(msg, 'error', esc(e.message));
        prog.textContent = '';
        if (w) w.close();
      }
      btn.disabled = false;
    }, w);
  }

  function runTests() {
    const p = current();
    if (!p || !B5.res) { needTree('p5TestMsg'); return; }
    const msg = el('p5TestMsg');
    clearMessages(msg);
    const A = prepare();
    const cands = [];
    const add = (name, tree) => { if (tree) cands.push({ name, tree }); };
    add(B5.trees.mlName, B5.trees.ml);
    if (B5.trees.constrained) add(B5.trees.constrainedName, B5.trees.constrained);
    const q = state.quick || {};
    if (q.parsimony && q.parsimony.tree) add(T('Parsimonia (Bloque 4)', 'Parsimony (Block 4)'), q.parsimony.tree);
    if (q.distance && q.distance.tree) add(T('Distancias (Bloque 4)', 'Distance (Block 4)'), q.distance.tree);
    if (B5.trees.partitioned) add(B5.trees.partitionedName, B5.trees.partitioned);
    if (cands.length < 2) {
      showMessage(msg, 'info', L2(
        'Hace falta más de un árbol para comparar. Busca un árbol restringido arriba, o vuelve al Bloque 4.',
        'More than one tree is needed. Search for a constrained tree above, or go back to Block 4.'));
      return;
    }
    const btn = el('p5RunTests');
    btn.disabled = true;
    const w = window.LABG ? LABG.work({ title: T('Comparando topologías (KH, SH y AU)', 'Comparing topologies (KH, SH and AU)'), delay: 300 }) : null;
    phyAfterPaint(() => {
      /* every candidate gets its branch lengths optimised under the same model,
         because a topology test compares optima, not whatever lengths came with
         the tree */
      const site = cands.map(c => {
        const f = Like.fit(Tree.clone(c.tree), A, Like.cloneSpec(B5.res.spec), { passes: 8, fixModel: true });
        c.lnL = f.lnL;
        return ML.siteLnL(f.tree, A, f.model);
      });
      const tt = ML.topologyTests(site, A.weights, { reps: +el('p5TestReps').value || 1000, seed: 5 });
      B5.tests = { cands: cands.map(c => c.name), tt };
      buildTable('p5TestTable', [
        { key: 'name', label: T('Árbol', 'Tree') },
        { key: 'lnL', label: 'lnL', num: true, fmt: v => fmtLnL(v) },
        { key: 'diff', label: 'Δ lnL', num: true, fmt: v => fmtFixed(v, 3) },
        { key: 'kh', label: 'KH p', num: true, fmt: v => fmtP(v) },
        { key: 'sh', label: 'SH p', num: true, fmt: v => fmtP(v) },
        { key: 'au', label: 'AU p', num: true, fmt: v => (v == null ? '—' : fmtP(v)) },
      ], cands.map((c, i) => ({
        name: c.name, lnL: tt.lnL[i], diff: tt.diff[i], kh: tt.kh[i], sh: tt.sh[i], au: tt.au[i],
      })), {});
      const rejected = cands.filter((c, i) => tt.sh[i] < 0.05);
      showMessage(msg, 'success', L2(
        `El mejor es «${esc(cands[tt.best].name)}». ` +
        (rejected.length
          ? `Con SH se rechazan al 5 %: ${rejected.map(c => '«' + esc(c.name) + '»').join(', ')}. Rechazar quiere decir que los datos <i>sí</i> distinguen, no solo que prefieren otro.`
          : 'Ninguno de los demás se rechaza al 5 %: los datos no distinguen entre ellos, y presentar uno solo como «el» árbol sería exagerar.'),
        `The best is "${esc(cands[tt.best].name)}". ` +
        (rejected.length
          ? `SH rejects at 5 %: ${rejected.map(c => '"' + esc(c.name) + '"').join(', ')}. Rejection means the data <i>do</i> distinguish, not merely that they prefer another.`
          : 'None of the others is rejected at 5 %: the data cannot tell them apart, and presenting one as "the" tree would overstate the case.')));
      btn.disabled = false;
      commit();
      if (w) w.done();
    }, w);
  }

  /* ================================================================
     5 · the figure
     ================================================================ */
  function drawInto(id, tree, labels, extra) {
    const host = el(id);
    if (!host) return;
    host.innerHTML = TreeView.render(tree, Object.assign({
      labels, width: 740, layout: 'rect', showSupport: true, supportAs: 'node',
      scaleLabel: T(' sustituciones/sitio', ' substitutions/site'),
    }, extra || {}));
  }

  function refreshTreeList() {
    const sel = el('p5WhichTree');
    if (!sel) return;
    const prev = sel.value;
    const opts = [];
    if (B5.trees.ml) opts.push(['ml', B5.trees.mlName]);
    if (B5.trees.partitioned) opts.push(['partitioned', B5.trees.partitionedName]);
    if (B5.trees.constrained) opts.push(['constrained', B5.trees.constrainedName]);
    sel.innerHTML = opts.map(([v, n]) => `<option value="${v}">${esc(n)}</option>`).join('');
    if (opts.some(o => o[0] === prev)) sel.value = prev;
    el('p5Figure').style.display = opts.length ? '' : 'none';
    if (opts.length) redraw();
  }

  /* the support the user asked to see, put on the drawn tree */
  function applyChosenSupport(tree, nTaxa) {
    const which = el('p5WhichSupport').value;
    if (which === 'none') { Tree.nodes(tree).forEach(n => { delete n.support; }); return null; }
    if (which === 'boot' && B5.support.boot) { Pars.applySupport(tree, B5.support.boot.freq, nTaxa); return 'bootstrap %'; }
    if (which === 'uf' && B5.support.uf) { Pars.applySupport(tree, B5.support.uf.freq, nTaxa); return 'UFBoot %'; }
    if (which === 'alrt' && B5.support.alrt) {
      const sp = Tree.splits(tree, nTaxa);
      sp.forEach((node, key) => { node.support = B5.support.alrt.alrt.get(key); });
      return 'SH-aLRT';
    }
    if (which === 'abayes' && B5.support.alrt) {
      const sp = Tree.splits(tree, nTaxa);
      sp.forEach((node, key) => { const v = B5.support.alrt.abayes.get(key); node.support = v == null ? null : v * 100; });
      return 'aBayes ×100';
    }
    Tree.nodes(tree).forEach(n => { delete n.support; });
    return null;
  }

  function redraw() {
    const p = current();
    const which = el('p5WhichTree').value;
    if (!p || !B5.trees[which]) return;
    let tree = Tree.clone(B5.trees[which]);
    applyChosenSupport(tree, p.taxa.length);
    const rootMode = el('p5Root').value;
    if (rootMode === 'midpoint') tree = Tree.midpointRoot(tree);
    else if (rootMode === 'outgroup') {
      const sel = [...el('p5Outgroup').selectedOptions].map(o => +o.value);
      if (sel.length) {
        const r = Tree.rootByOutgroup(tree, sel);
        tree = r.tree;
        if (!r.monophyletic) showMessage(el('p5BootMsg'), 'warning', L2(
          'El grupo externo elegido no es monofilético en este árbol. Eso ya es un resultado.',
          'The chosen outgroup is not monophyletic on this tree. That is a result in itself.'));
      }
    }
    Tree.ladderize(tree);
    drawInto('p5MainFig', tree, p.taxa, {
      layout: el('p5Layout').value,
      supportAs: el('p5SupportAs').value,
      showSupport: el('p5WhichSupport').value !== 'none',
      width: el('p5Layout').value === 'circular' ? 720 : 780,
    });
    B5.drawn = tree;
  }

  function exportAs(kind) {
    const p = current();
    if (!p) return;
    const tree = B5.drawn || B5.trees.ml;
    const base = slug(p.name) + '_ml';
    if (kind === 'newick' && tree) download(Tree.writeNewick(tree, { labels: p.taxa }), base + '.nwk');
    else if (kind === 'nexus' && tree) {
      download(SeqIO.writeNexus(p.taxa, p.seqs, {
        type: p.type,
        trees: [{ name: 'ML', newick: Tree.writeNewick(tree, { labels: p.taxa }) }],
      }), base + '.nex');
    } else if (kind === 'svg') download(el('p5MainFig').innerHTML, base + '.svg', 'image/svg+xml');
    else if (kind === 'boottrees' && B5.bootTrees) download(B5.bootTrees.join('\n') + '\n', base + '_bootstrap.nwk');
    else if (kind === 'sitelnl' && B5.res) {
      const A = prepare();
      const s = ML.siteLnL(B5.trees.ml, A, B5.res.model);
      const lines = ['sitio,patron,peso,lnL'];
      for (let i = 0; i < A.nSites; i++) {
        const pat = A.patternOfSite[i];
        lines.push(`${i + 1},${pat + 1},${A.weights[pat]},${s[pat].toFixed(8)}`);
      }
      download(lines.join('\n') + '\n', base + '_lnL_por_sitio.csv', 'text/csv');
    }
  }

  /* ================================================================
     hand-off
     ================================================================ */
  function commit() {
    const p = current();
    if (!p || !B5.res) return;
    state.ml = {
      partition: p.name,
      lnL: B5.res.lnL, spec: Like.cloneSpec(B5.res.spec), k: B5.res.k,
      AIC: B5.res.AIC, AICc: B5.res.AICc, BIC: B5.res.BIC,
      tree: B5.trees.ml,
      newick: Tree.writeNewick(B5.trees.ml, { labels: p.taxa }),
      name: B5.trees.mlName,
      visited: B5.res.visited,
      support: {
        boot: B5.support.boot ? { reps: B5.support.boot.reps, freq: [...B5.support.boot.freq.entries()] } : null,
        uf: B5.support.uf ? { reps: B5.support.uf.reps, freq: [...B5.support.uf.freq.entries()] } : null,
        alrt: B5.support.alrt ? { alrt: [...B5.support.alrt.alrt.entries()], abayes: [...B5.support.alrt.abayes.entries()] } : null,
      },
      partitioned: B5.part5 ? {
        lnL: B5.part5.lnL, k: B5.part5.k, BIC: B5.part5.BIC,
        partitions: B5.part5.partitions.map(q => ({ name: q.name, rate: q.rate, lnL: q.lnL, spec: q.spec })),
        newick: Tree.writeNewick(B5.trees.partitioned, { labels: p.taxa }),
      } : null,
      constrained: B5.constrained ? {
        taxa: B5.constrained.taxa, lnL: B5.constrained.res.lnL,
        newick: Tree.writeNewick(B5.trees.constrained, { labels: p.taxa }),
      } : null,
      tests: B5.tests,
    };
    enableStep(6, true);
  }

  /* ================================================================
     wiring
     ================================================================ */
  function refresh() {
    const has = parts().length > 0;
    el('p5NoData').style.display = has ? 'none' : '';
    ['p5Search', 'p5Partitioned', 'p5Support', 'p5Tests'].forEach(id => { el(id).style.display = has ? '' : 'none'; });
    if (!has) { el('p5Figure').style.display = 'none'; return; }
    const ps = parts();
    const sel = el('p5Part');
    sel.innerHTML = ps.map((q, i) => `<option value="${i}">${esc(q.name)} · ${q.taxa.length} × ${q.length}</option>`).join('');
    if (B5.part < 0 || B5.part >= ps.length) B5.part = 0;
    sel.value = String(B5.part);
    el('p5Partitioned').style.display = ps.length > 1 ? '' : 'none';
    fillModels();
    if (block3Spec() && !B5.res) useBlock3();
    fillConstraintList();
    refreshTreeList();
  }

  function init() {
    if (!el('panel-5')) return;
    el('p5Part').addEventListener('change', () => {
      B5.part = +el('p5Part').value;
      B5.A = null; B5.res = null; B5.trees = {}; B5.support = { boot: null, uf: null, alrt: null };
      B5.constrained = null; B5.tests = null; B5.bootTrees = null;
      el('p5Pane').style.display = 'none';
      el('p5Verdict').style.display = 'none';
      el('p5Figure').style.display = 'none';
      refresh();
    });
    el('p5FromBlock3').addEventListener('click', useBlock3);
    el('p5Run').addEventListener('click', runSearch);
    el('p5Cancel').addEventListener('click', () => { B5.cancelled = true; });
    el('p5RunPart').addEventListener('click', runPartitioned);
    el('p5RunBoot').addEventListener('click', runBootstrap);
    el('p5RunUF').addEventListener('click', runUFBoot);
    el('p5RunALRT').addEventListener('click', runALRT);
    el('p5CancelBoot').addEventListener('click', () => { B5.cancelled = true; });
    el('p5RunConstraint').addEventListener('click', runConstrained);
    el('p5RunTests').addEventListener('click', runTests);
    ['p5WhichTree', 'p5Layout', 'p5WhichSupport', 'p5SupportAs'].forEach(id => el(id).addEventListener('change', redraw));
    el('p5Root').addEventListener('change', () => {
      el('p5OutgroupWrap').style.display = el('p5Root').value === 'outgroup' ? '' : 'none';
      redraw();
    });
    el('p5Outgroup').addEventListener('change', redraw);
    el('p5Redraw').addEventListener('click', redraw);
    els('[data-p5export]').forEach(b => b.addEventListener('click', () => exportAs(b.dataset.p5export)));
    el('p5ToBlock6').addEventListener('click', () => {
      const b = document.querySelector('.step-btn[data-step="6"]');
      if (b && !b.disabled) goStep(6);
      else showMessage(el('p5BootMsg'), 'info', L2(
        'El Bloque 6 (inferencia bayesiana) llega en la etapa siguiente; este árbol y este modelo serán su punto de partida.',
        'Block 6 (Bayesian inference) arrives in the next stage; this tree and this model will be its starting point.'));
    });
    document.addEventListener('stepchange', e => { if (e.detail.step === 5) refresh(); });
    document.addEventListener('langchange', () => { if (parts().length) { refresh(); if (B5.trees.ml) renderSupportTable(); } });
    document.addEventListener('themechange', () => { if (B5.drawn) redraw(); });
  }
  document.addEventListener('DOMContentLoaded', init);

  Object.assign(B5, { refresh, runSearch, runPartitioned, runBootstrap, runUFBoot, runALRT, runConstrained, runTests, redraw, commit, current });
})();
