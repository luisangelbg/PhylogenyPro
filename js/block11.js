/* PhylogenyPro — molecular phylogenetics in the browser.
   Copyright (C) 2026  Luis Ángel Barrera-Guzmán

   This program is free software: you can redistribute it and/or modify it under
   the terms of the GNU General Public License as published by the Free Software
   Foundation, either version 3 of the License, or (at your option) any later
   version. It is distributed in the hope that it will be useful, but WITHOUT ANY
   WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
   A PARTICULAR PURPOSE. See the GNU General Public License, in the file LICENSE
   at the root of this program, or <https://www.gnu.org/licenses/>. */

/* PhylogenyPro — Block 11: comparing trees.

   The engine is js/treecmp.js. What this file adds is the reading, and here it
   mostly consists of refusing to let a number mean more than it does:

     · Robinson–Foulds is a count of splits, not a measure of how wrong a tree
       is; one misplaced tip can send it to its maximum, and the block says so
       beside the number;
     · concordance is not support. A branch with 100 % bootstrap and 35 % gCF
       is well estimated and badly supported at the same time, and that is the
       commonest thing in a phylogenomic paper;
     · a quartet species tree is a heuristic search here, not ASTRAL's exact
       program, so the block reports the score and how it got there;
     · and D ≠ 0 has three explanations, only one of which is introgression. */

(function () {

  const B11 = {
    extra: [], dists: null, tangle: null, conc: null,
    quartet: null, net: null, abba: null, scan: null,
  };
  window.B11 = B11;

  /* ================================================================
     the trees
     ================================================================ */
  function taxa() {
    if (state.data && state.data.parts && state.data.parts.length) {
      const p = state.data.parts[0];
      if (p && p.taxa) return p.taxa;
    }
    return (state.data && state.data.taxa) || [];
  }
  /* every tree the app has made, plus whatever has been pasted in */
  function sources() {
    const out = [];
    const add = (id, name, tree) => { if (tree) out.push({ id, name, tree }); };
    if (state.quick && state.quick.distance) add('distance', T('B4 · distancias', 'B4 · distance'), state.quick.distance.tree);
    if (state.quick && state.quick.parsimony) add('parsimony', T('B4 · parsimonia', 'B4 · parsimony'), state.quick.parsimony.tree);
    if (state.quick && state.quick.consensus) add('consensus4', T('B4 · consenso', 'B4 · consensus'), state.quick.consensus.tree);
    if (state.ml) add('ml', T('B5 · verosimilitud', 'B5 · likelihood'), state.ml.tree);
    if (state.bayes) add('bayes', T('B6 · consenso bayesiano', 'B6 · Bayesian consensus'), state.bayes.consensus);
    if (state.dated) add('dated', T('B7 · fechado', 'B7 · dated'), state.dated.tree);
    B11.extra.forEach((t, i) => add('extra' + i, t.name, t.tree));
    return out;
  }
  const treeById = id => (sources().find(s => s.id === id) || {}).tree;
  function nTaxa() { return taxa().length || (sources()[0] ? Tree.tips(sources()[0].tree).length : 0); }

  /* a pasted tree carries its own tip labels; they have to be turned into the
     same indices the app's trees use, or nothing can be compared */
  function readTrees() {
    const msg = el('p11DataMsg');
    msg.innerHTML = '';
    const text = el('p11Text').value;
    const lines = text.split('\n').map(l => l.trim()).filter(l => l && l.indexOf('(') >= 0);
    if (!lines.length) {
      showMessage(msg, 'error', L2('No hay ningún árbol en el cuadro. Se espera formato Newick, un árbol por línea.',
        'There is no tree in the box. Newick is expected, one tree per line.'));
      return;
    }
    const names = taxa();
    if (!names.length) {
      showMessage(msg, 'error', L2('Primero hace falta un conjunto de datos con nombres de taxones (Bloque 2).',
        'A data set with taxon names is needed first (Block 2).'));
      return;
    }
    const norm = s => String(s).replace(/[\s_]+/g, ' ').trim().toLowerCase();
    const byNorm = new Map(names.map((nm, i) => [norm(nm), i]));
    B11.extra = [];
    const problems = [];
    lines.forEach((l, i) => {
      let t;
      try { t = Tree.parseNewick(l, names); } catch (e) { problems.push(`#${i + 1}: ${e.message}`); return; }
      /* parseNewick assigns an index by name when it can; a label it does not
         know gets a fresh index, which is exactly what must not happen */
      const tips = Tree.tips(t);
      const unknown = [];
      tips.forEach(tp => {
        const lab = tp.label;
        if (lab != null && byNorm.has(norm(lab))) tp.tip = byNorm.get(norm(lab));
        else if (lab != null) unknown.push(lab);
      });
      if (unknown.length) { problems.push(`#${i + 1}: ${unknown.slice(0, 4).join(', ')}`); return; }
      if (tips.length !== names.length) {
        problems.push(T(`#${i + 1}: ${tips.length} puntas de ${names.length}`, `#${i + 1}: ${tips.length} tips of ${names.length}`));
        return;
      }
      B11.extra.push({ name: T(`pegado ${i + 1}`, `pasted ${i + 1}`), tree: t });
    });
    if (problems.length) {
      showMessage(msg, B11.extra.length ? 'warn' : 'error', L2(
        `Estos árboles no se pudieron usar: ${problems.slice(0, 5).join('; ')}. Los nombres de las puntas tienen que ser los del Bloque 2; se perdonan espacios y guiones bajos, nada más.`,
        `These trees could not be used: ${problems.slice(0, 5).join('; ')}. The tip names have to be the ones from Block 2; spaces and underscores are forgiven, nothing else is.`));
    }
    if (B11.extra.length) {
      showMessage(msg, 'success', L2(`${B11.extra.length} árbol(es) añadido(s).`, `${B11.extra.length} tree(s) added.`));
    }
    B11.dists = B11.conc = B11.quartet = B11.net = null;
    refresh();
    commit();
  }

  function showTrees() {
    const src = sources();
    const n = nTaxa();
    buildTable('p11TreeTable', [
      { key: 'name', label: T('árbol', 'tree') },
      { key: 'tips', label: T('puntas', 'tips'), num: true },
      { key: 'internal', label: T('divisiones internas', 'internal splits'), num: true },
      { key: 'len', label: T('longitud total', 'total length'), num: true },
      { key: 'rooted', label: T('raíz', 'root') },
    ], src.map(s => ({
      name: s.name, tips: Tree.tips(s.tree).length,
      internal: Cmp.splitTable(s.tree, n).keys.length,
      len: fmtNum(Tree.totalLength(s.tree)),
      rooted: s.tree.children && s.tree.children.length === 2 ? T('sí', 'yes') : T('no (basal politómica)', 'no (basal polytomy)'),
    })));
    el('p11DataInfo').textContent = T(`${src.length} árboles de ${n} puntas`, `${src.length} trees of ${n} tips`);
    ['p11TangA', 'p11TangB', 'p11RefTree', 'p11QStart', 'p11NetTree'].forEach(id => {
      const sel = el(id);
      if (!sel) return;
      const prev = sel.value;
      sel.innerHTML = src.map(s => `<option value="${s.id}">${esc(s.name)}</option>`).join('');
      if (src.some(s => s.id === prev)) sel.value = prev;
    });
    /* the two sides of a tanglegram must not be the same tree, and they will be
       if the second menu was filled when only one tree existed */
    if (src.length > 1 && el('p11TangB').value === el('p11TangA').value) {
      el('p11TangB').value = src.find(s => s.id !== el('p11TangA').value).id;
    }
    const names = taxa();
    ['p11P1', 'p11P2', 'p11P3', 'p11O'].forEach((id, k) => {
      const sel = el(id);
      const prev = sel.value;
      sel.innerHTML = names.map((nm, i) => `<option value="${i}">${esc(nm)}</option>`).join('');
      if (prev && names[+prev] != null) sel.value = prev;
      else sel.value = String(Math.min(names.length - 1, k === 3 ? names.length - 1 : k));
    });
  }

  /* ================================================================
     2 · distances
     ================================================================ */
  const MEASURES = [
    ['rf', 'Robinson–Foulds', (a, b, n) => Cmp.rf(a, b, n).rf, true],
    ['rfnorm', T('RF normalizada', 'RF normalised'), (a, b, n) => Cmp.rf(a, b, n).normalised, false],
    ['quartet', T('cuartetos', 'quartets'), (a, b, n) => Cmp.quartetDistance(a, b, n).different, true],
    ['wrf', T('RF pesada', 'weighted RF'), (a, b, n) => Cmp.weightedRF(a, b, n), false],
    ['kf', T('puntaje de rama', 'branch score'), (a, b, n) => Cmp.branchScore(a, b, n), false],
    ['path', T('caminos', 'path difference'), (a, b, n) => Cmp.pathDistance(a, b, n, { weighted: false }), false],
  ];
  function runDist() {
    const msg = el('p11DistMsg');
    msg.innerHTML = '';
    const src = sources();
    if (src.length < 2) return;
    const n = nTaxa();
    el('p11DistProgress').innerHTML = L2('calculando…', 'computing…');
    const w = window.LABG ? LABG.work({ title: T('Distancias entre árboles', 'Distances between trees'), delay: 300 }) : null;
    phyAfterPaint(() => {
      const t0 = performance.now();
      const mats = {};
      MEASURES.forEach(([id, , fn]) => {
        const M = Array.from({ length: src.length }, () => new Float64Array(src.length));
        for (let i = 0; i < src.length; i++) {
          for (let j = i + 1; j < src.length; j++) {
            const v = fn(src[i].tree, src[j].tree, n);
            M[i][j] = v; M[j][i] = v;
          }
        }
        mats[id] = M;
      });
      B11.dists = { mats, names: src.map(s => s.name), ids: src.map(s => s.id), n };
      el('p11DistProgress').innerHTML = L2(`${src.length} árboles en ${((performance.now() - t0) / 1000).toFixed(1)} s`,
        `${src.length} trees in ${((performance.now() - t0) / 1000).toFixed(1)} s`);
      showDist();
      commit();
      if (w) w.done();
    }, w);
  }
  function showDist() {
    if (!B11.dists) return;
    const { mats, names, n } = B11.dists;
    const which = el('p11Measure').value;
    const M = mats[which];
    const def = MEASURES.find(m => m[0] === which);
    const k = names.length;
    let sum = 0, cnt = 0, hi = 0, lo = Infinity, hiPair = '', loPair = '';
    for (let i = 0; i < k; i++) for (let j = i + 1; j < k; j++) {
      sum += M[i][j]; cnt++;
      if (M[i][j] > hi) { hi = M[i][j]; hiPair = `${names[i]} / ${names[j]}`; }
      if (M[i][j] < lo) { lo = M[i][j]; loPair = `${names[i]} / ${names[j]}`; }
    }
    const maxRF = 2 * (n - 3);
    const totalQ = (() => { let c = 1; for (let i = 0; i < 4; i++) c = c * (n - i) / (i + 1); return Math.round(c); })();
    statTiles('p11DistTiles', [
      [T('Árboles comparados', 'Trees compared'), k, `${cnt} ` + T('pares', 'pairs')],
      [T('Más parecidos', 'Most alike'), fmtNum(lo), loPair, 'ok'],
      [T('Más distintos', 'Least alike'), fmtNum(hi), hiPair, hi > 0 ? 'bad' : ''],
      [T('Media', 'Mean'), fmtNum(sum / Math.max(1, cnt)),
        which === 'rf' ? T(`de un máximo de ${maxRF}`, `of a maximum of ${maxRF}`)
          : which === 'quartet' ? T(`de ${totalQ} cuartetos`, `of ${totalQ} quartets`) : ''],
    ]);
    el('p11DistPane').style.display = '';
    el('p11DistFig').innerHTML = Plots11.distanceHeat(M, names, {
      integer: !!def[3],
      caption: T(`${def[1]} entre cada par`, `${def[1]} between each pair`),
    });
    registerFigure('distances', el('p11DistFig').innerHTML);

    const rows = [];
    for (let i = 0; i < k; i++) for (let j = i + 1; j < k; j++) {
      const r = { a: names[i], b: names[j] };
      MEASURES.forEach(([id, , , int]) => { r[id] = int ? String(Math.round(mats[id][i][j])) : fmtNum(mats[id][i][j]); });
      rows.push(r);
    }
    buildTable('p11DistTable', [
      { key: 'a', label: T('árbol', 'tree') }, { key: 'b', label: T('contra', 'against') },
      ...MEASURES.map(([id, label]) => ({ key: id, label, num: true })),
    ], rows, { limit: 60 });

    const v = el('p11DistVerdict');
    v.style.display = '';
    const es = [], en = [];
    const rfMean = (() => { let s = 0, c = 0; for (let i = 0; i < k; i++) for (let j = i + 1; j < k; j++) { s += mats.rf[i][j]; c++; } return s / Math.max(1, c); })();
    const qMean = (() => { let s = 0, c = 0; for (let i = 0; i < k; i++) for (let j = i + 1; j < k; j++) { s += mats.quartet[i][j]; c++; } return s / Math.max(1, c); })();
    es.push(`En promedio dos árboles difieren en ${fmtFixed(rfMean, 1)} divisiones de ${maxRF} posibles, y discrepan en ${fmtFixed(qMean, 0)} de ${totalQ} cuartetos (${fmtPct(qMean / totalQ)}).`);
    en.push(`On average two trees differ by ${fmtFixed(rfMean, 1)} splits of a possible ${maxRF}, and disagree on ${fmtFixed(qMean, 0)} of ${totalQ} quartets (${fmtPct(qMean / totalQ)}).`);
    es.push('Las dos cifras dicen cosas distintas a propósito: Robinson–Foulds cuenta divisiones y una sola punta mal colocada puede llevarla al máximo, mientras que la de cuartetos mide cuánta de la estructura cambia realmente. Cuando RF es grande y la de cuartetos pequeña, el desacuerdo está en una o dos puntas inestables, no en la forma del árbol.');
    en.push('The two numbers say different things on purpose: Robinson–Foulds counts splits and one misplaced tip can take it to the maximum, while the quartet distance measures how much of the structure really changes. When RF is large and the quartet distance small, the disagreement is in one or two unstable tips, not in the shape of the tree.');
    if (hi === 0) {
      es.push('Aquí todos los árboles tienen la misma topología: lo que difiere, si algo, son las longitudes de rama.');
      en.push('Here every tree has the same topology: what differs, if anything, are the branch lengths.');
    }
    v.innerHTML = `<b>${L2('Cómo se lee', 'How to read it')}</b> ${L2(es.join(' '), en.join(' '))}`;
    el('p11Export').style.display = '';
  }

  /* ================================================================
     3 · tanglegram
     ================================================================ */
  function runTangle() {
    const a = treeById(el('p11TangA').value), b = treeById(el('p11TangB').value);
    if (!a || !b) return;
    const n = nTaxa(), names = taxa();
    const res = el('p11Untangle').checked ? Cmp.untangle(a, b, {})
      : { tree: Tree.clone(b), crossings: Cmp.crossings(Cmp.leafOrder(a), Cmp.leafOrder(b)),
          before: Cmp.crossings(Cmp.leafOrder(a), Cmp.leafOrder(b)) };
    B11.tangle = res;
    const r = Cmp.rf(a, b, n);
    const q = Cmp.quartetDistance(a, b, n);
    /* the tips that moved are the story, so they are counted and named */
    const o1 = Cmp.leafOrder(a), o2 = Cmp.leafOrder(res.tree);
    const pos = new Map(o1.map((t, i) => [t, i]));
    const moved = o2.map((t, i) => ({ t, d: Math.abs(i - (pos.get(t) || 0)) }))
      .filter(x => x.d > 0).sort((x, y) => y.d - x.d);
    statTiles('p11TangTiles', [
      [T('Cruces', 'Crossings'), res.crossings,
        el('p11Untangle').checked ? T(`${res.before} sin desenredar`, `${res.before} untangled`) : '',
        res.crossings === 0 ? 'ok' : ''],
      ['Robinson–Foulds', r.rf, T(`de ${r.maxRF}`, `of ${r.maxRF}`)],
      [T('Cuartetos discrepantes', 'Disagreeing quartets'), q.different, fmtPct(q.normalised)],
      [T('Puntas que se mueven', 'Tips that move'), moved.length,
        moved.slice(0, 2).map(x => names[x.t] || x.t).join(', ')],
    ]);
    el('p11TanglePane').style.display = '';
    el('p11TangleFig').innerHTML = Plots11.tanglegram(a, res.tree, {
      labels: names, crossings: res.crossings,
      leftLabel: el('p11TangA').selectedOptions[0].textContent,
      rightLabel: el('p11TangB').selectedOptions[0].textContent,
      highlight: moved.slice(0, 3).map(x => x.t),
    });
    registerFigure('tanglegram', el('p11TangleFig').innerHTML);
    const v = el('p11TangVerdict');
    v.style.display = '';
    const es = [], en = [];
    if (r.rf === 0) {
      es.push('Los dos árboles tienen exactamente la misma topología, así que el tanglegrama no puede cruzarse: cualquier cruce que quede es cuestión de cómo se dibujan, no de lo que dicen.');
      en.push('The two trees have exactly the same topology, so the tanglegram cannot cross: any crossing left is a matter of how they are drawn, not of what they say.');
    } else {
      es.push(`Difieren en ${r.rf} divisiones y en ${q.different} cuartetos. Las ${moved.length} puntas resaltadas son las que cambian de sitio${moved.length ? `: ${moved.slice(0, 4).map(x => names[x.t] || x.t).join(', ')}` : ''}.`);
      en.push(`They differ by ${r.rf} splits and ${q.different} quartets. The ${moved.length} highlighted tips are the ones that move${moved.length ? `: ${moved.slice(0, 4).map(x => names[x.t] || x.t).join(', ')}` : ''}.`);
      if (q.normalised < 0.05 && r.rf > 2) {
        es.push('Nótese que la discrepancia de cuartetos es pequeña aunque RF no lo sea: eso es la firma de unas pocas puntas que se mueven mucho, no de dos historias distintas.');
        en.push('Notice that the quartet disagreement is small although RF is not: that is the signature of a few tips that move a lot, not of two different histories.');
      }
      if (res.crossings === 0) {
        es.push('<b>Y una advertencia sobre esta figura en concreto:</b> aquí no queda ningún cruce y los árboles sí difieren. Girando las ramas se consigue que las dos listas de puntas queden en el mismo orden, pero el desacuerdo sigue ahí, en cómo se agrupan por dentro. Un tanglegram sin cruces no demuestra que dos árboles sean iguales; para eso está la RF de al lado.');
        en.push('<b>And a warning about this figure in particular:</b> no crossings are left here and the trees do differ. Rotating branches can put both lists of tips in the same order, but the disagreement is still there, in how they group inside. A tanglegram with no crossings does not prove two trees are the same; the RF beside it is what does.');
      }
    }
    v.innerHTML = `<b>${L2('Cómo se lee', 'How to read it')}</b> ${L2(es.join(' '), en.join(' '))}`;
    el('p11Export').style.display = '';
  }

  /* ================================================================
     4 · concordance
     ================================================================ */
  function alignment() {
    if (!state.data || !state.data.parts || !state.data.parts.length) return null;
    const p = state.data.parts[0];
    if (!p || !p.seqs) return null;
    const names = taxa();
    const byName = new Map(p.taxa.map((t, i) => [t, (p.seqs[i] || '').toUpperCase()]));
    return names.map(nm => byName.get(nm) || null);
  }
  function runConc() {
    const msg = el('p11ConcMsg');
    msg.innerHTML = '';
    const ref = treeById(el('p11RefTree').value);
    if (!ref) return;
    const n = nTaxa(), names = taxa();
    const others = sources().filter(s => s.tree !== ref).map(s => s.tree);
    if (others.length < 2) {
      showMessage(msg, 'error', L2('Hacen falta al menos dos árboles además del de referencia.',
        'At least two trees besides the reference are needed.'));
      return;
    }
    el('p11ConcProgress').innerHTML = L2('calculando…', 'computing…');
    const w = window.LABG ? LABG.work({ title: T('Factores de concordancia (gCF y sCF)', 'Concordance factors (gCF and sCF)'), delay: 300 }) : null;
    phyAfterPaint(() => {
      const gcf = Cmp.gcf(ref, others, n);
      const seqs = alignment();
      const q = Math.max(10, Math.min(2000, +el('p11Quartets').value || 100));
      const scf = seqs && seqs.every(Boolean) ? Cmp.scf(ref, seqs, n, { quartets: q, seed: 99 }) : null;
      B11.conc = { ref, gcf, scf, nGenes: others.length, names };
      el('p11ConcProgress').innerHTML = '';
      showConc();
      if (!scf) {
        showMessage(msg, 'info', L2(
          'Para el sCF hace falta el alineamiento del Bloque 2 con todas las puntas; sin él sólo se calcula el gCF.',
          'The sCF needs the Block 2 alignment with every tip; without it only the gCF is computed.'));
      }
      commit();
      if (w) w.done();
    }, w);
  }
  const cladeName = (clade, names) => clade.map(i => (names[i] || i)).slice(0, 2).join(', ') +
    (clade.length > 2 ? ` +${clade.length - 2}` : '');

  function showConc() {
    const { gcf, scf, nGenes, names } = B11.conc;
    const low = gcf.filter(b => b.gCF < 50).length;
    const veryLow = gcf.filter(b => b.gCF < 33.34).length;
    const meanG = gcf.reduce((a, b) => a + b.gCF, 0) / Math.max(1, gcf.length);
    const meanS = scf ? scf.reduce((a, b) => a + b.sCF, 0) / Math.max(1, scf.length) : null;
    statTiles('p11ConcTiles', [
      [T('Ramas internas', 'Internal branches'), gcf.length, T(`sobre ${nGenes} árboles`, `over ${nGenes} trees`)],
      ['gCF ' + T('medio', 'mean'), fmtFixed(meanG, 1) + ' %', '', meanG > 50 ? 'ok' : 'bad'],
      ...(meanS != null ? [['sCF ' + T('medio', 'mean'), fmtFixed(meanS, 1) + ' %', '', meanS > 50 ? 'ok' : 'bad']] : []),
      [T('Ramas por debajo del azar', 'Branches below chance'), veryLow,
        T(`${low} por debajo del 50 %`, `${low} below 50 %`), veryLow ? 'bad' : ''],
    ]);
    const rows = gcf.slice().sort((a, b) => a.gCF - b.gCF).map(b => ({
      label: cladeName(b.clade, names),
      concordant: b.gCF, alt1: b.gDF1, alt2: b.gDF2, other: b.gDFP,
    }));
    el('p11GcfPane').style.display = '';
    el('p11GcfFig').innerHTML = Plots11.concordance(rows, {
      width: 460, labelWidth: 130,
      legend: [T('con la rama', 'with the branch'), T('alternativa 1', 'alternative 1'),
        T('alternativa 2', 'alternative 2'), T('no decide', 'undecided')],
    });
    registerFigure('gcf', el('p11GcfFig').innerHTML);
    if (scf) {
      const byKey = new Map(scf.map(s => [s.key, s]));
      el('p11ScatterPane').style.display = '';
      el('p11ScatterFig').innerHTML = Plots11.concordanceScatter(
        gcf.filter(b => byKey.has(b.key)).map(b => ({ x: b.gCF, y: byKey.get(b.key).sCF })),
        { width: 400, height: 340 });
      registerFigure('concordance', el('p11ScatterFig').innerHTML);
    } else el('p11ScatterPane').style.display = 'none';

    const byKey = scf ? new Map(scf.map(s => [s.key, s])) : null;
    buildTable('p11ConcTable', [
      { key: 'clade', label: T('clado', 'clade') },
      { key: 'nTips', label: T('puntas', 'tips'), num: true },
      { key: 'len', label: T('longitud', 'length'), num: true },
      { key: 'gcf', label: 'gCF', num: true },
      { key: 'gdf', label: 'gDF1 / gDF2', num: true },
      { key: 'scf', label: 'sCF', num: true },
      { key: 'sdf', label: 'sDF1 / sDF2', num: true },
    ], gcf.slice().sort((a, b) => a.gCF - b.gCF).map(b => {
      const s = byKey ? byKey.get(b.key) : null;
      return {
        clade: cladeName(b.clade, names), nTips: b.nTips, len: fmtNum(b.length),
        gcf: fmtFixed(b.gCF, 1), gdf: `${fmtFixed(b.gDF1, 1)} / ${fmtFixed(b.gDF2, 1)}`,
        scf: s ? fmtFixed(s.sCF, 1) : '—',
        sdf: s ? `${fmtFixed(s.sDF1, 1)} / ${fmtFixed(s.sDF2, 1)}` : '—',
        _class: b.gCF < 33.34 ? 'row-flag' : '',
      };
    }), { limit: 50 });

    const v = el('p11ConcVerdict');
    v.style.display = '';
    const es = [], en = [];
    es.push(`De ${gcf.length} ramas internas, ${low} tienen menos de la mitad de los árboles de genes a favor y ${veryLow} están por debajo del tercio que daría el puro azar.`);
    en.push(`Of ${gcf.length} internal branches, ${low} have less than half the gene trees behind them and ${veryLow} are below the third that pure chance would give.`);
    if (veryLow) {
      es.push('Una rama por debajo del azar no es una rama débil: es una rama que la mayoría de la evidencia contradice, y que probablemente esté mal.');
      en.push('A branch below chance is not a weak branch: it is a branch most of the evidence contradicts, and is probably wrong.');
    }
    if (meanS != null && Math.abs(meanS - meanG) > 15) {
      es.push(`El gCF medio (${fmtFixed(meanG, 1)} %) y el sCF medio (${fmtFixed(meanS, 1)} %) se separan bastante. Cuando el de sitios es el más bajo, suele significar que cada gen tiene poca señal aunque su árbol salga resuelto; cuando el más bajo es el de genes, que el conflicto está entre genes y no dentro de ellos.`);
      en.push(`The mean gCF (${fmtFixed(meanG, 1)} %) and the mean sCF (${fmtFixed(meanS, 1)} %) are some way apart. When the site one is lower it usually means each gene has little signal even though its tree comes out resolved; when the gene one is lower, the conflict is between genes rather than inside them.`);
    }
    es.push('Y lo que no dice ninguna de las dos: <b>la concordancia no es apoyo</b>. Una rama puede tener 100 % de arranque y 35 % de concordancia a la vez, y eso significa que los datos son muchos y están divididos.');
    en.push("And what neither of them says: <b>concordance is not support</b>. A branch can have 100 % bootstrap and 35 % concordance at once, and that means the data are plentiful and divided.");
    v.innerHTML = `<b>${L2('Cómo se lee', 'How to read it')}</b> ${L2(es.join(' '), en.join(' '))}`;
    el('p11Export').style.display = '';
  }

  /* ================================================================
     5 · quartet species tree
     ================================================================ */
  function runQuartet() {
    const msg = el('p11QMsg');
    msg.innerHTML = '';
    const src = sources();
    if (src.length < 3) {
      showMessage(msg, 'error', L2('Hacen falta al menos tres árboles para que contar cuartetos signifique algo.',
        'At least three trees are needed for counting quartets to mean anything.'));
      return;
    }
    const n = nTaxa();
    el('p11QProgress').innerHTML = L2('contando cuartetos…', 'counting quartets…');
    const w = window.LABG ? LABG.work({
      title: T('Árbol de especies a partir de los cuartetos', 'Species tree from the quartets'),
      message: T('contando cuartetos…', 'counting quartets…'), delay: 300,
    }) : null;
    phyAfterPaint(() => {
      const t0 = performance.now();
      const qc = Cmp.quartetCounts(src.map(s => s.tree), n);
      const start = treeById(el('p11QStart').value) || src[0].tree;
      const search = Cmp.quartetSearch(start, qc, n, {});
      const scores = src.map(s => Object.assign({ name: s.name }, Cmp.quartetScore(s.tree, qc, n)));
      B11.quartet = { qc, search, scores, n, secs: (performance.now() - t0) / 1000 };
      el('p11QProgress').innerHTML = L2(`${qc.counts.size} cuartetos en ${B11.quartet.secs.toFixed(1)} s`,
        `${qc.counts.size} quartets in ${B11.quartet.secs.toFixed(1)} s`);
      showQuartet();
      commit();
      if (w) w.done();
    }, w);
  }
  function showQuartet() {
    const { qc, search, scores, n } = B11.quartet;
    const best = scores.reduce((a, b) => (b.score > a.score ? b : a), scores[0]);
    const improved = search.score > best.score + 1e-9;
    statTiles('p11QTiles', [
      [T('Cuartetos', 'Quartets'), qc.counts.size, T(`${qc.nTrees} árboles`, `${qc.nTrees} trees`)],
      [T('Mejor de los árboles dados', 'Best of the given trees'), fmtPct(best.proportion), best.name],
      [T('Árbol hallado por la búsqueda', 'Tree the search found'), fmtPct(search.proportion),
        T(`${search.rounds} rondas`, `${search.rounds} rounds`), improved ? 'ok' : ''],
      [T('Votos recogidos', 'Votes collected'), `${Math.round(search.score)} / ${Math.round(search.total)}`, ''],
    ]);
    buildTable('p11QTable', [
      { key: 'name', label: T('árbol', 'tree') },
      { key: 'score', label: T('cuartetos a favor', 'quartets in favour'), num: true },
      { key: 'prop', label: T('proporción', 'proportion'), num: true },
      { key: 'rf', label: T('RF contra el hallado', 'RF against the one found'), num: true },
    ], scores.slice().sort((a, b) => b.score - a.score).map(s => {
      const src = sources().find(x => x.name === s.name);
      return {
        name: s.name, score: Math.round(s.score), prop: fmtPct(s.proportion),
        rf: src ? Cmp.rf(src.tree, search.tree, n).rf : '—',
        _class: s.name === best.name ? 'row-best' : '',
      };
    }));
    el('p11QPane').style.display = '';
    el('p11QFig').innerHTML = TreeView.render(search.tree, {
      labels: taxa(), width: 720, cladogram: true, showSupport: false,
      title: T('árbol de especies por cuartetos', 'quartet species tree'),
    });
    registerFigure('quartetTree', el('p11QFig').innerHTML);
    const v = el('p11QVerdict');
    v.style.display = '';
    const es = [], en = [];
    es.push(`El árbol que más votos recoge se lleva ${fmtPct(search.proportion)} de los cuartetos decididos. Un valor cercano a uno significa que los árboles casi no discrepan; uno cerca de un tercio, que discrepan tanto como discreparían al azar.`);
    en.push(`The tree that collects the most votes takes ${fmtPct(search.proportion)} of the decided quartets. A value near one means the trees barely disagree; one near a third, that they disagree as much as chance would.`);
    if (improved) {
      es.push(`La búsqueda encontró un árbol mejor que cualquiera de los que se le dieron (${Math.round(search.score)} frente a ${Math.round(best.score)} votos), lo cual es la razón de ser del método: el árbol de especies no tiene por qué ser el árbol de ningún gen.`);
      en.push(`The search found a tree better than any it was given (${Math.round(search.score)} against ${Math.round(best.score)} votes), which is the point of the method: the species tree need not be any gene's tree.`);
    } else {
      es.push('La búsqueda no mejoró ninguno de los árboles dados, así que el mejor de ellos ya estaba en un óptimo local de este criterio.');
      en.push('The search did not improve on any of the given trees, so the best of them was already at a local optimum of this criterion.');
    }
    es.push('Conviene recordar que esto es una búsqueda por intercambios de vecinos desde un punto de partida, no el programa dinámico exacto de ASTRAL: con muchos taxones puede quedarse en un óptimo local.');
    en.push("It is worth remembering that this is a nearest-neighbour search from a starting point, not ASTRAL's exact dynamic program: with many taxa it can settle in a local optimum.");
    v.innerHTML = `<b>${L2('Cómo se lee', 'How to read it')}</b> ${L2(es.join(' '), en.join(' '))}`;
    el('p11Export').style.display = '';
  }

  /* ================================================================
     6 · split network
     ================================================================ */
  function runNet() {
    const msg = el('p11NetMsg');
    msg.innerHTML = '';
    const n = nTaxa(), names = taxa();
    const kind = el('p11NetKind').value;
    if (kind === 'consensus') {
      const src = sources();
      if (src.length < 2) return;
      const th = Math.max(0, Math.min(1, +el('p11Thresh').value || 0.3));
      const cn = Cmp.consensusNetwork(src.map(s => s.tree), n, { threshold: th });
      /* the drawing needs a circular ordering, and any of the trees gives one */
      const order = Cmp.circularOrder(src[0].tree);
      const asCircular = cn.splits.filter(s => isInterval(s.set, order, n))
        .map(s => Object.assign({}, s, intervalOf(s.set, order, n), { weight: s.frequency }));
      const cycle = Cmp.networkOuterCycle({ splits: asCircular, order, nTaxa: n });
      B11.net = { kind, cn, cycle, order, dropped: cn.splits.length - asCircular.length, n };
    } else {
      const t = treeById(el('p11NetTree').value) || sources()[0].tree;
      const v = Cmp.pathVector(t, n, true);
      const D = Array.from({ length: n }, () => new Float64Array(n));
      let k = 0;
      for (let i = 0; i < n; i++) for (let j = i + 1; j < n; j++) { D[i][j] = D[j][i] = v[k++]; }
      const order = Cmp.circularOrder(t);
      const cs = Cmp.circularSplits(D, order, {});
      const cycle = Cmp.networkOuterCycle(cs);
      B11.net = { kind, cs, cycle, order, n };
    }
    showNet();
    commit();
  }
  /* whether a split is a contiguous arc of the circular ordering, which is what
     the drawing needs; the ones that are not cannot be drawn on this cycle */
  function isInterval(set, order, n) {
    const inside = order.map(t => !!set[t]);
    let runs = 0;
    for (let i = 0; i < n; i++) if (inside[i] && !inside[(i - 1 + n) % n]) runs++;
    return runs === 1;
  }
  function intervalOf(set, order, n) {
    const inside = order.map(t => !!set[t]);
    let start = 0;
    for (let i = 0; i < n; i++) if (inside[i] && !inside[(i - 1 + n) % n]) start = i;
    let end = start;
    while (inside[(end + 1) % n]) end = (end + 1) % n;
    return { start, end };
  }
  function showNet() {
    const { kind, cn, cs, cycle, n } = B11.net;
    const names = taxa();
    if (kind === 'consensus') {
      const conflicting = cn.splits.filter(s => s.conflictsWith.length).length;
      statTiles('p11NetTiles', [
        [T('Divisiones', 'Splits'), cn.splits.length, T(`umbral ${cn.threshold}`, `threshold ${cn.threshold}`)],
        [T('En conflicto', 'In conflict'), conflicting, conflicting ? T('la red no es un árbol', 'the network is not a tree') : T('la red es un árbol', 'the network is a tree'),
          conflicting ? 'bad' : 'ok'],
        [T('En todos los árboles', 'In every tree'), cn.splits.filter(s => s.frequency > 0.999).length, ''],
        [T('No dibujables en este orden', 'Not drawable in this order'), B11.net.dropped, ''],
      ]);
      buildTable('p11NetTable', [
        { key: 'clade', label: T('división', 'split') },
        { key: 'size', label: T('puntas', 'tips'), num: true },
        { key: 'freq', label: T('frecuencia', 'frequency'), num: true },
        { key: 'len', label: T('longitud media', 'mean length'), num: true },
        { key: 'conf', label: T('choca con', 'clashes with'), num: true },
      ], cn.splits.map(s => ({
        clade: cladeName(s.members, names), size: s.members.length,
        freq: fmtPct(s.frequency), len: fmtNum(s.weight), conf: s.conflictsWith.length,
        _class: s.conflictsWith.length ? 'row-flag' : '',
      })), { limit: 60 });
    } else {
      statTiles('p11NetTiles', [
        [T('Divisiones con peso', 'Splits with weight'), cs.splits.length, ''],
        [T('Descartadas por peso negativo', 'Dropped for negative weight'), cs.clipped,
          cs.clipped ? T('las distancias no son circulares en este orden', 'the distances are not circular in this order') : '',
          cs.clipped ? 'bad' : 'ok'],
        [T('Peso descartado', 'Weight dropped'), fmtNum(cs.clippedMass), ''],
        [T('Cierre del ciclo', 'Cycle closure'), cycle.closure.toExponential(1), T('debe ser cero', 'should be zero')],
      ]);
      buildTable('p11NetTable', [
        { key: 'clade', label: T('división', 'split') },
        { key: 'size', label: T('puntas', 'tips'), num: true },
        { key: 'w', label: T('peso', 'weight'), num: true },
      ], cs.splits.map(s => ({
        clade: cladeName(s.members, names), size: s.members.length, w: fmtNum(s.weight),
      })), { limit: 60 });
    }
    const conflictIdx = new Set();
    if (kind === 'consensus') {
      cycle.splits.forEach((s, i) => { if (s.conflictsWith && s.conflictsWith.length) conflictIdx.add(i); });
    }
    el('p11NetPane').style.display = '';
    el('p11NetFig').innerHTML = Plots11.network(cycle, {
      labels: names, conflicting: [...conflictIdx],
      caption: kind === 'consensus'
        ? T('ciclo exterior de la red de consenso; en color, las divisiones que chocan con otra',
            'outer cycle of the consensus network; coloured, the splits that clash with another')
        : T('ciclo exterior de la descomposición circular', 'outer cycle of the circular decomposition'),
    });
    registerFigure('network', el('p11NetFig').innerHTML);

    const v = el('p11NetVerdict');
    v.style.display = '';
    const es = [], en = [];
    if (kind === 'consensus') {
      const conflicting = cn.splits.filter(s => s.conflictsWith.length).length;
      es.push(conflicting
        ? `${conflicting} de las ${cn.splits.length} divisiones chocan con alguna otra, así que ningún árbol puede contenerlas todas: eso es exactamente lo que la red dibuja y un consenso tendría que tirar.`
        : `Ninguna de las ${cn.splits.length} divisiones choca con otra a este umbral, de modo que la red es un árbol. Bajar el umbral hará aparecer las divisiones minoritarias y, con ellas, el conflicto.`);
      en.push(conflicting
        ? `${conflicting} of the ${cn.splits.length} splits clash with another, so no tree can contain them all: that is exactly what the network draws and a consensus would have to throw away.`
        : `None of the ${cn.splits.length} splits clashes with another at this threshold, so the network is a tree. Lowering the threshold will bring in the minority splits and, with them, the conflict.`);
      if (B11.net.dropped) {
        es.push(`${B11.net.dropped} división(es) no forman un arco contiguo del orden circular elegido y no se pueden dibujar en este ciclo; siguen en la tabla.`);
        en.push(`${B11.net.dropped} split(s) are not a contiguous arc of the chosen circular ordering and cannot be drawn on this cycle; they are still in the table.`);
      }
    } else {
      es.push(cs.clipped
        ? `${cs.clipped} divisiones salieron con peso negativo, por un total de ${fmtNum(cs.clippedMass)}: las distancias no son exactamente circulares en este orden y esa parte se descarta. Cuanto mayor sea ese total, menos fiel es el dibujo.`
        : 'Ninguna división salió con peso negativo, así que la descomposición reproduce las distancias exactamente y el dibujo no esconde nada.');
      en.push(cs.clipped
        ? `${cs.clipped} splits came out with negative weight, ${fmtNum(cs.clippedMass)} in total: the distances are not exactly circular in this ordering and that part is dropped. The larger that total, the less faithful the drawing.`
        : 'No split came out with negative weight, so the decomposition reproduces the distances exactly and the drawing hides nothing.');
      es.push('Sobre las distancias de un árbol esta descomposición devuelve ese árbol y nada más; las cajas sólo aparecen cuando las distancias no son las de ningún árbol.');
      en.push('On a tree\'s distances this decomposition gives back that tree and nothing else; the boxes only appear when the distances are no tree\'s.');
    }
    v.innerHTML = `<b>${L2('Cómo se lee', 'How to read it')}</b> ${L2(es.join(' '), en.join(' '))}`;
    el('p11Export').style.display = '';
  }

  /* ================================================================
     7 · ABBA-BABA
     ================================================================ */
  function runAbba(scan) {
    const msg = el('p11AbbaMsg');
    msg.innerHTML = '';
    const seqs = alignment();
    const names = taxa();
    if (!seqs || !seqs.every(Boolean)) {
      showMessage(msg, 'error', L2('Hace falta el alineamiento del Bloque 2 con todas las puntas.',
        'The Block 2 alignment with every tip is needed.'));
      return;
    }
    const blocks = Math.max(5, Math.min(200, +el('p11Blocks').value || 20));
    const o = +el('p11O').value;
    el('p11AbbaProgress').innerHTML = L2('contando sitios…', 'counting sites…');
    const w = window.LABG ? LABG.work({
      title: scan ? T('Barrido ABBA-BABA', 'ABBA-BABA scan') : T('Prueba ABBA-BABA', 'ABBA-BABA test'),
      message: T('contando sitios…', 'counting sites…'), delay: 300,
    }) : null;
    phyAfterPaint(() => {
      if (scan) {
        const rows = Cmp.dScan(seqs, names, o, { blocks, minSites: 20 });
        B11.scan = { rows, outgroup: o };
        B11.abba = null;
      } else {
        const p1 = +el('p11P1').value, p2 = +el('p11P2').value, p3 = +el('p11P3').value;
        if (new Set([p1, p2, p3, o]).size < 4) {
          el('p11AbbaProgress').innerHTML = '';
          showMessage(msg, 'error', L2('Los cuatro taxones tienen que ser distintos.', 'The four taxa have to be different.'));
          if (w) w.close();
          return;
        }
        B11.abba = Object.assign({ p1, p2, p3, o }, Cmp.dStatistic(seqs, p1, p2, p3, o, { blocks }));
        B11.scan = null;
      }
      el('p11AbbaProgress').innerHTML = '';
      showAbba();
      commit();
      if (w) w.done();
    }, w);
  }
  function showAbba() {
    const names = taxa();
    const lab = r => `${names[r.p1]} , ${names[r.p2]} | ${names[r.p3]}`;
    if (B11.abba) {
      const d = B11.abba;
      statTiles('p11AbbaTiles', [
        ['ABBA', d.ABBA, T('P2 con P3', 'P2 with P3')],
        ['BABA', d.BABA, T('P1 con P3', 'P1 with P3')],
        ['D', fmtFixed(d.D, 4), `± ${fmtFixed(d.se, 4)}`, Math.abs(d.z) > 1.96 ? 'bad' : 'ok'],
        ['z', fmtFixed(d.z, 3), pEq(d.p), Math.abs(d.z) > 1.96 ? 'bad' : 'ok'],
      ]);
      el('p11AbbaPane').style.display = '';
      el('p11AbbaFig').innerHTML = Plots11.dBars([{ label: lab(d), D: d.D, se: d.se, z: d.z }], { width: 640 });
      el('p11AbbaTable').innerHTML = '';
    } else if (B11.scan) {
      const rows = B11.scan.rows;
      const sig = rows.filter(r => Math.abs(r.z) > 1.96).length;
      statTiles('p11AbbaTiles', [
        [T('Tríos probados', 'Trios tested'), rows.length, T(`exogrupo ${names[B11.scan.outgroup]}`, `outgroup ${names[B11.scan.outgroup]}`)],
        [T('Con |z| > 1.96', 'With |z| > 1.96'), sig, fmtPct(sig / Math.max(1, rows.length)), sig ? 'bad' : 'ok'],
        [T('|D| mayor', 'Largest |D|'), rows.length ? fmtFixed(Math.abs(rows[0].D), 4) : '—',
          rows.length ? lab(rows[0]) : ''],
        [T('Esperados por azar al 5 %', 'Expected by chance at 5 %'), fmtFixed(0.05 * rows.length, 1),
          T('sin corregir por comparaciones múltiples', 'uncorrected for multiple comparisons')],
      ]);
      el('p11AbbaPane').style.display = '';
      el('p11AbbaFig').innerHTML = Plots11.dBars(
        rows.slice(0, 18).map(r => ({ label: lab(r), D: r.D, se: r.se, z: r.z })), { width: 680 });
      buildTable('p11AbbaTable', [
        { key: 'trio', label: 'P1, P2 | P3' },
        { key: 'abba', label: 'ABBA', num: true }, { key: 'baba', label: 'BABA', num: true },
        { key: 'd', label: 'D', num: true }, { key: 'se', label: T('error est.', 'std. error'), num: true },
        { key: 'z', label: 'z', num: true }, { key: 'p', label: 'p', num: true },
      ], rows.map(r => ({
        trio: lab(r), abba: r.ABBA, baba: r.BABA, d: fmtFixed(r.D, 4),
        se: fmtFixed(r.se, 4), z: fmtFixed(r.z, 3), p: pEq(r.p),
        _class: Math.abs(r.z) > 1.96 ? 'row-flag' : '',
      })), { limit: 60 });
    }
    registerFigure('abba', el('p11AbbaFig').innerHTML);

    const v = el('p11AbbaVerdict');
    v.style.display = '';
    const es = [], en = [];
    if (B11.abba) {
      const d = B11.abba;
      const who = d.D > 0 ? names[d.p2] : names[d.p1];
      if (Math.abs(d.z) > 1.96) {
        es.push(`D = ${fmtFixed(d.D, 4)} con z = ${fmtFixed(d.z, 2)} (${pEq(d.p)}): hay ${d.D > 0 ? 'más' : 'menos'} sitios ABBA que BABA de los que la clasificación de linajes produciría, y el exceso apunta a un intercambio entre <b>${who}</b> y ${names[d.p3]}.`);
        en.push(`D = ${fmtFixed(d.D, 4)} with z = ${fmtFixed(d.z, 2)} (${pEq(d.p)}): there are ${d.D > 0 ? 'more' : 'fewer'} ABBA sites than BABA than lineage sorting would produce, and the excess points to exchange between <b>${who}</b> and ${names[d.p3]}.`);
      } else {
        es.push(`D = ${fmtFixed(d.D, 4)} con z = ${fmtFixed(d.z, 2)} (${pEq(d.p)}): la asimetría no se distingue de cero, que es lo que predice la clasificación incompleta de linajes por sí sola.`);
        en.push(`D = ${fmtFixed(d.D, 4)} with z = ${fmtFixed(d.z, 2)} (${pEq(d.p)}): the asymmetry is indistinguishable from zero, which is what incomplete lineage sorting alone predicts.`);
      }
      if (d.ABBA + d.BABA < 100) {
        es.push(`Sólo hay ${d.ABBA + d.BABA} sitios informativos para este trío; con tan pocos, ni un D grande significa gran cosa.`);
        en.push(`There are only ${d.ABBA + d.BABA} informative sites for this trio; with so few, even a large D means little.`);
      }
    } else if (B11.scan) {
      const rows = B11.scan.rows;
      const sig = rows.filter(r => Math.abs(r.z) > 1.96).length;
      es.push(`${sig} de ${rows.length} tríos pasan de |z| = 1.96. Al 5 % y sin corregir se esperarían ${fmtFixed(0.05 * rows.length, 1)} sólo por azar, y además los tríos no son independientes entre sí: comparten taxones y comparten historia.`);
      en.push(`${sig} of ${rows.length} trios exceed |z| = 1.96. At 5 % and uncorrected, ${fmtFixed(0.05 * rows.length, 1)} would be expected by chance alone, and the trios are not independent of each other either: they share taxa and they share history.`);
    }
    es.push('Y la advertencia que acompaña siempre a esta prueba: <b>D ≠ 0 tiene tres explicaciones</b> y sólo una es introgresión. Un exogrupo que no lo es del todo, tasas de sustitución muy desiguales entre P1 y P2, y errores de alineamiento o de ensamblado producen la misma asimetría. Descartarlas es parte del resultado, no un extra.');
    en.push('And the warning that always goes with this test: <b>D ≠ 0 has three explanations</b> and only one is introgression. An outgroup that is not quite one, very unequal substitution rates between P1 and P2, and alignment or assembly errors all produce the same asymmetry. Ruling them out is part of the result, not an extra.');
    v.innerHTML = `<b>${L2('Cómo se lee', 'How to read it')}</b> ${L2(es.join(' '), en.join(' '))}`;
    el('p11Export').style.display = '';
  }

  /* ================================================================
     keeping the result, and taking it away
     ================================================================ */
  function registerFigure(name, svgText) {
    state.figures = state.figures || {};
    state.figures['b11_' + name] = { block: 11, name, svg: svgText };
  }
  function commit() {
    const names = taxa();
    state.compare = {
      trees: sources().map(s => s.id),
      distances: B11.dists ? {
        names: B11.dists.names,
        matrices: Object.fromEntries(Object.keys(B11.dists.mats).map(k =>
          [k, B11.dists.mats[k].map(r => Array.from(r))])),
      } : null,
      concordance: B11.conc ? {
        nGenes: B11.conc.nGenes,
        branches: B11.conc.gcf.map(b => {
          const s = B11.conc.scf ? B11.conc.scf.find(x => x.key === b.key) : null;
          return {
            clade: b.clade.map(i => names[i] || i), nTips: b.nTips, length: b.length,
            gCF: b.gCF, gDF1: b.gDF1, gDF2: b.gDF2, gDFP: b.gDFP, decisive: b.decisive,
            sCF: s ? s.sCF : null, sDF1: s ? s.sDF1 : null, sDF2: s ? s.sDF2 : null,
            sites: s ? s.sites : null,
          };
        }),
      } : null,
      quartet: B11.quartet ? {
        nQuartets: B11.quartet.qc.counts.size,
        score: B11.quartet.search.score, total: B11.quartet.search.total,
        proportion: B11.quartet.search.proportion,
        newick: Tree.writeNewick(B11.quartet.search.tree, { labels: names }),
        scores: B11.quartet.scores.map(s => ({ name: s.name, score: s.score, proportion: s.proportion })),
      } : null,
      network: B11.net ? (B11.net.kind === 'consensus' ? {
        kind: 'consensus', threshold: B11.net.cn.threshold,
        splits: B11.net.cn.splits.map(s => ({ clade: s.members.map(i => names[i] || i),
          frequency: s.frequency, weight: s.weight, conflicts: s.conflictsWith.length })),
      } : {
        kind: 'circular', clipped: B11.net.cs.clipped, clippedMass: B11.net.cs.clippedMass,
        splits: B11.net.cs.splits.map(s => ({ clade: s.members.map(i => names[i] || i), weight: s.weight })),
      }) : null,
      abba: B11.abba ? Object.assign({}, B11.abba, {
        taxa: [names[B11.abba.p1], names[B11.abba.p2], names[B11.abba.p3], names[B11.abba.o]] }) : null,
      abbaScan: B11.scan ? B11.scan.rows.slice(0, 200).map(r => ({
        p1: names[r.p1], p2: names[r.p2], p3: names[r.p3],
        ABBA: r.ABBA, BABA: r.BABA, D: r.D, se: r.se, z: r.z, p: r.p })) : null,
    };
  }

  function exportAs(kind) {
    const stamp = new Date().toISOString().slice(0, 10);
    const csv = rows => rows.map(r => r.join(',')).join('\n');
    const names = taxa();
    if (kind === 'dist' && B11.dists) {
      const rows = [['tree_a', 'tree_b'].concat(MEASURES.map(m => m[0]))];
      const { mats, names: tn } = B11.dists;
      for (let i = 0; i < tn.length; i++) for (let j = i + 1; j < tn.length; j++) {
        rows.push([tn[i], tn[j]].concat(MEASURES.map(m => mats[m[0]][i][j])));
      }
      download(csv(rows), `phylogenypro_distancias_arboles_${stamp}.csv`, 'text/csv');
    } else if (kind === 'conc' && B11.conc) {
      const rows = [['clade', 'n_tips', 'length', 'gCF', 'gDF1', 'gDF2', 'gDFP', 'decisive', 'sCF', 'sDF1', 'sDF2', 'sites']];
      B11.conc.gcf.forEach(b => {
        const s = B11.conc.scf ? B11.conc.scf.find(x => x.key === b.key) : null;
        rows.push(['"' + b.clade.map(i => names[i] || i).join(' ') + '"', b.nTips, b.length,
          b.gCF, b.gDF1, b.gDF2, b.gDFP, b.decisive,
          s ? s.sCF : '', s ? s.sDF1 : '', s ? s.sDF2 : '', s ? s.sites : '']);
      });
      download(csv(rows), `phylogenypro_concordancia_${stamp}.csv`, 'text/csv');
    } else if (kind === 'splits' && B11.net) {
      const rows = [['clade', 'n_tips', 'weight', 'frequency', 'conflicts']];
      const list = B11.net.kind === 'consensus' ? B11.net.cn.splits : B11.net.cs.splits;
      list.forEach(s => rows.push(['"' + s.members.map(i => names[i] || i).join(' ') + '"',
        s.members.length, s.weight, s.frequency == null ? '' : s.frequency,
        s.conflictsWith ? s.conflictsWith.length : '']));
      download(csv(rows), `phylogenypro_divisiones_${stamp}.csv`, 'text/csv');
    } else if (kind === 'abba' && (B11.abba || B11.scan)) {
      const rows = [['P1', 'P2', 'P3', 'outgroup', 'ABBA', 'BABA', 'D', 'se', 'z', 'p', 'sites']];
      const push = r => rows.push([names[r.p1], names[r.p2], names[r.p3],
        names[B11.scan ? B11.scan.outgroup : r.o], r.ABBA, r.BABA, r.D, r.se, r.z, r.p, r.sites]);
      if (B11.abba) push(B11.abba);
      if (B11.scan) B11.scan.rows.forEach(push);
      download(csv(rows), `phylogenypro_abba_baba_${stamp}.csv`, 'text/csv');
    } else if (kind === 'trees') {
      const txt = sources().map(s => Tree.writeNewick(s.tree, { labels: names })).join('\n');
      download(txt, `phylogenypro_arboles_${stamp}.tre`, 'text/plain');
    } else if (kind === 'svg') {
      const figs = Object.keys(state.figures || {}).filter(k => k.indexOf('b11_') === 0);
      figs.forEach(k => download('<?xml version="1.0" encoding="UTF-8"?>\n' + state.figures[k].svg,
        `phylogenypro_${k}_${stamp}.svg`, 'image/svg+xml'));
    }
  }

  /* ================================================================
     wiring
     ================================================================ */
  function refresh() {
    const src = sources();
    const has = src.length >= 2;
    el('p11NoTree').style.display = src.length ? 'none' : '';
    el('p11Data').style.display = '';
    ['p11Dist', 'p11Tangle', 'p11Conc', 'p11Quartet', 'p11Net', 'p11Abba'].forEach(id => {
      el(id).style.display = has ? '' : 'none';
    });
    if (src.length) showTrees();
  }

  function init() {
    if (!el('panel-11')) return;
    el('p11Upload').addEventListener('click', () => el('p11File').click());
    el('p11File').addEventListener('change', e => {
      const f = e.target.files && e.target.files[0];
      if (!f) return;
      const r = new FileReader();
      r.onload = () => { el('p11Text').value = r.result; readTrees(); };
      r.readAsText(f);
    });
    el('p11Example').addEventListener('click', () => {
      el('p11Text').value = (window.EXAMPLES && EXAMPLES.arbolesGenes ? EXAMPLES.arbolesGenes : '').trim();
      readTrees();
    });
    el('p11Read').addEventListener('click', readTrees);
    el('p11RunDist').addEventListener('click', runDist);
    el('p11Measure').addEventListener('change', () => { if (B11.dists) showDist(); });
    el('p11RunTangle').addEventListener('click', runTangle);
    el('p11RunConc').addEventListener('click', runConc);
    el('p11RunQuartet').addEventListener('click', runQuartet);
    el('p11NetKind').addEventListener('change', () => {
      const circ = el('p11NetKind').value === 'circular';
      el('p11ThreshWrap').style.display = circ ? 'none' : '';
      el('p11NetTreeWrap').style.display = circ ? '' : 'none';
    });
    el('p11RunNet').addEventListener('click', runNet);
    el('p11RunAbba').addEventListener('click', () => runAbba(false));
    el('p11RunScan').addEventListener('click', () => runAbba(true));
    els('[data-p11export]').forEach(b => b.addEventListener('click', () => exportAs(b.dataset.p11export)));
    el('p11ToBlock12').addEventListener('click', () => {
      const b = document.querySelector('.step-btn[data-step="12"]');
      if (b && !b.disabled) goStep(12);
      else showMessage(el('p11DataMsg'), 'info', L2(
        'El Bloque 12 (figuras e informe) llega en la etapa siguiente.',
        'Block 12 (figures and report) arrives in the next stage.'));
    });
    document.addEventListener('stepchange', e => { if (e.detail.step === 11) refresh(); });
    document.addEventListener('langchange', () => {
      if (!sources().length) return;
      refresh();
      if (B11.dists) showDist();
      if (B11.conc) showConc();
      if (B11.quartet) showQuartet();
      if (B11.net) showNet();
      if (B11.abba || B11.scan) showAbba();
    });
    document.addEventListener('themechange', () => {
      if (B11.dists) showDist();
      if (B11.conc) showConc();
    });
  }
  document.addEventListener('DOMContentLoaded', init);

  Object.assign(B11, { refresh, readTrees, runDist, runTangle, runConc, runQuartet, runNet, runAbba, commit, sources });
})();
