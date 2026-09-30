/* PhylogenyPro — molecular phylogenetics in the browser.
   Copyright (C) 2026  Luis Ángel Barrera-Guzmán

   This program is free software: you can redistribute it and/or modify it under
   the terms of the GNU General Public License as published by the Free Software
   Foundation, either version 3 of the License, or (at your option) any later
   version. It is distributed in the hope that it will be useful, but WITHOUT ANY
   WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
   A PARTICULAR PURPOSE. See the GNU General Public License, in the file LICENSE
   at the root of this program, or <https://www.gnu.org/licenses/>. */

/* PhylogenyPro — Block 9: trait evolution.

   The engine is js/traits.js. What this file adds is the reading, and in this
   block the reading matters more than usual, because the numbers are easy to
   over-claim:

     · OU beating BM by AIC on twenty species is not evidence of stabilising
       selection; it is what happens when a third parameter is allowed;
     · K and λ answer different questions, so they disagree without either
       being wrong, and the block says so instead of reporting one;
     · a significant PGLS slope is not a cause, and an ancestral state at 62 %
       is not a fact — the pies exist so that it does not look like one;
     · ape's ace and phytools' fitMk quote log-likelihoods that differ by
       log(k), which is only the root prior. Both are shown, named. */

(function () {

  const B9 = {
    source: 'dated', tree: null, table: null, matched: null,
    cont: null, signal: null, two: null, disc: null, sim: null,
    contVar: null, discVar: null,
  };
  window.B9 = B9;

  /* ================================================================
     where the tree comes from
     ================================================================ */
  /* Every tree the app has made, with the partition whose taxon list names its
     tips: the trees carry indices, not names, so the two travel together. */
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
    if (state.quick && state.quick.parsimony && state.quick.parsimony.tree) {
      out.push(['parsimony', T('Bloque 4 · árbol de parsimonia', 'Block 4 · parsimony tree'), state.quick.parsimony.tree, false, state.quick.partition]);
    }
    return out;
  }
  function chosen() {
    const s = sources().find(x => x[0] === B9.source) || sources()[0];
    return s || null;
  }
  function chosenTree() {
    const s = chosen();
    return s ? Tree.clone(s[2]) : null;
  }
  function isDated() { const s = chosen(); return !!(s && s[3]); }
  function unit() { return (state.dated && state.dated.unit) || 'Ma'; }

  /* the taxon names of the partition the chosen tree was built on */
  function taxaOf() {
    const s = chosen();
    const name = s ? s[4] : null;
    if (name && state.data && state.data.parts) {
      const p = state.data.parts.find(q => q.name === name);
      if (p && p.taxa) return p.taxa;
    }
    return (state.data && state.data.taxa) || [];
  }
  /* the tip labels of the chosen tree, in tip-row order. Trees made inside the
     app carry an index into the partition's taxon list; a tree read from a file
     carries its own label, so both are tried. */
  function treeLabels(tree) {
    const taxa = taxaOf();
    const out = [];
    Tree.tips(tree).forEach(t => { out[t.tip] = t.label || taxa[t.tip] || null; });
    return out;
  }

  /* ================================================================
     reading the table
     ================================================================ */
  /* A character table is not an alignment, so it gets its own small reader: any
     of the three usual separators, the first column the name, and a column that
     is numbers all the way down taken as continuous. */
  function parseTable(text) {
    const lines = String(text).replace(/\r/g, '').split('\n').map(l => l.trim()).filter(l => l && !/^#/.test(l));
    if (lines.length < 2) return null;
    const sep = lines[0].indexOf('\t') >= 0 ? '\t' : lines[0].indexOf(';') >= 0 ? ';' : ',';
    const head = lines[0].split(sep).map(s => s.trim().replace(/^["']|["']$/g, ''));
    const rows = lines.slice(1).map(l => l.split(sep).map(s => s.trim().replace(/^["']|["']$/g, '')));
    const missing = v => v === '' || v === '?' || v === 'NA' || v === 'na' || v === '-' || v == null;
    const names = rows.map(r => r[0]);
    const cols = [];
    for (let j = 1; j < head.length; j++) {
      const raw = rows.map(r => (r[j] == null ? '' : r[j]));
      const present = raw.filter(v => !missing(v));
      if (!present.length) continue;
      const numeric = present.every(v => /^[-+]?(\d+\.?\d*|\.\d+)([eE][-+]?\d+)?$/.test(v));
      cols.push({
        name: head[j] || ('col' + j),
        type: numeric ? 'continuous' : 'discrete',
        values: raw.map(v => (missing(v) ? null : (numeric ? +v : v))),
        levels: numeric ? null : [...new Set(present)].sort(),
        nMissing: raw.filter(missing).length,
      });
    }
    return { names, columns: cols, nRows: rows.length };
  }

  /* line the table up with the tree, and say plainly what did not line up */
  function matchToTree(table, tree) {
    const labels = treeLabels(tree);
    const byName = new Map();
    table.names.forEach((nm, i) => byName.set(String(nm), i));
    /* a second try with underscores and spaces made equal, because that is the
       difference that separates most tables from most Newick files */
    const norm = s => String(s).replace(/[\s_]+/g, ' ').trim().toLowerCase();
    const byNorm = new Map();
    table.names.forEach((nm, i) => { if (!byNorm.has(norm(nm))) byNorm.set(norm(nm), i); });

    const rowOf = new Array(labels.length).fill(-1);
    const missing = [];
    labels.forEach((lab, tipRow) => {
      if (lab == null) return;
      let i = byName.has(lab) ? byName.get(lab) : (byNorm.has(norm(lab)) ? byNorm.get(norm(lab)) : -1);
      rowOf[tipRow] = i;
      if (i < 0) missing.push(lab);
    });
    const used = new Set(rowOf.filter(i => i >= 0));
    const extra = table.names.filter((nm, i) => !used.has(i));
    return { rowOf, labels, missing, extra, matched: used.size, nTips: labels.length };
  }

  /* the values of one column, in tip-row order, and the tips that have one */
  function columnFor(col, opts) {
    opts = opts || {};
    const m = B9.matched;
    const vals = [], keep = [];
    m.rowOf.forEach((r, tipRow) => {
      const v = r >= 0 ? col.values[r] : null;
      if (v == null) return;
      keep.push(tipRow);
      vals.push(opts.log ? Math.log(v) : v);
    });
    return { values: vals, tips: keep };
  }
  /* the tree pruned to the tips that have a value, with the tip rows renumbered
     so the trait vector and the tree agree */
  function prunedTree(keep) {
    const tree = chosenTree();
    const set = new Set(keep);
    const idx = new Map(keep.map((t, i) => [t, i]));
    function walk(nd) {
      if (nd.tip != null) return set.has(nd.tip) ? { tip: idx.get(nd.tip), label: nd.label, children: [] } : null;
      const kids = [];
      nd.children.forEach(c => {
        const k = walk(c.node);
        if (k) kids.push({ node: k, len: c.len || 0 });
      });
      if (!kids.length) return null;
      if (kids.length === 1) {
        /* a node with one daughter is not a node: its branch absorbs the edge */
        const only = kids[0];
        only.node.pendingLen = (only.node.pendingLen || 0) + only.len;
        return only.node;
      }
      const nd2 = { label: nd.label, support: nd.support, children: kids };
      kids.forEach(c => { c.len += (c.node.pendingLen || 0); delete c.node.pendingLen; });
      return nd2;
    }
    const r = walk(tree);
    if (r) { delete r.pendingLen; }
    return r;
  }
  function analysisFor(col, opts) {
    const c = columnFor(col, opts);
    return { x: c.values, tips: c.tips, tree: prunedTree(c.tips),
      labels: c.tips.map(t => B9.matched.labels[t]) };
  }

  function showTable() {
    const t = B9.table, m = B9.matched;
    if (!t) return;
    const cols = [
      { key: 'name', label: T('carácter', 'character') },
      { key: 'type', label: T('tipo', 'type') },
      { key: 'detail', label: T('valores', 'values') },
      { key: 'miss', label: T('sin dato', 'missing'), num: true },
    ];
    const rows = t.columns.map(c => {
      const present = c.values.filter(v => v != null);
      let detail;
      if (c.type === 'continuous') {
        const mn = Math.min.apply(null, present), mx = Math.max.apply(null, present);
        detail = `${fmtNum(mn)} – ${fmtNum(mx)}`;
      } else {
        detail = c.levels.join(', ');
      }
      return {
        name: c.name,
        type: c.type === 'continuous' ? T('continuo', 'continuous') : T(`discreto · ${c.levels.length} estados`, `discrete · ${c.levels.length} states`),
        detail, miss: c.nMissing,
      };
    });
    buildTable('p9DataTable', cols, rows);

    const info = [];
    info.push(T(`${m.matched} de ${m.nTips} puntas con datos`, `${m.matched} of ${m.nTips} tips with data`));
    if (m.missing.length) info.push(T(`${m.missing.length} sin fila`, `${m.missing.length} with no row`));
    if (m.extra.length) info.push(T(`${m.extra.length} filas sin punta`, `${m.extra.length} rows with no tip`));
    el('p9DataInfo').textContent = info.join(' · ');

    const msg = el('p9DataMsg');
    msg.innerHTML = '';
    if (m.matched < 4) {
      showMessage(msg, 'error', L2(
        `Solo ${m.matched} puntas quedaron emparejadas. Los nombres de la primera columna tienen que ser los del árbol; se aceptan diferencias de espacios y guiones bajos, pero nada más.`,
        `Only ${m.matched} tips were matched. The names in the first column have to be the tree's; differences in spaces and underscores are forgiven, nothing else is.`));
      return;
    }
    if (m.missing.length) {
      showMessage(msg, 'warn', L2(
        `Sin fila en la tabla: ${m.missing.slice(0, 8).join(', ')}${m.missing.length > 8 ? '…' : ''}. Cada análisis podará el árbol a las puntas que sí tengan el carácter que se esté usando, y lo dirá.`,
        `No row in the table: ${m.missing.slice(0, 8).join(', ')}${m.missing.length > 8 ? '…' : ''}. Each analysis prunes the tree to the tips that do have the character in use, and says so.`));
    }
    if (m.extra.length) {
      showMessage(msg, 'info', L2(
        `Filas que no corresponden a ninguna punta: ${m.extra.slice(0, 8).join(', ')}${m.extra.length > 8 ? '…' : ''}. Se ignoran.`,
        `Rows matching no tip: ${m.extra.slice(0, 8).join(', ')}${m.extra.length > 8 ? '…' : ''}. They are ignored.`));
    }
    if (!isDated()) {
      showMessage(msg, 'warn', L2(
        'Este árbol no está fechado, así que las ramas están en sustituciones por sitio. Todo lo de este bloque sigue funcionando, pero σ² y las tasas quedan <b>por sustitución</b>, no por unidad de tiempo, y no se pueden comparar con las de otro estudio. Para eso hace falta el cronograma del Bloque 7.',
        'This tree is not dated, so its branches are in substitutions per site. Everything in this block still works, but σ² and the rates come out <b>per substitution</b>, not per unit of time, and cannot be compared with another study\'s. The chronogram of Block 7 is what that needs.'));
    }

    /* fill the character menus */
    const cont = t.columns.filter(c => c.type === 'continuous');
    const disc = t.columns.filter(c => c.type === 'discrete');
    fillSelect('p9ContVar', cont);
    fillSelect('p9VarY', cont);
    fillSelect('p9VarX', cont, 1);
    fillSelect('p9DiscVar', disc);
    el('p9Cont').style.display = cont.length ? '' : 'none';
    el('p9Signal').style.display = cont.length ? '' : 'none';
    el('p9Two').style.display = cont.length > 1 ? '' : 'none';
    el('p9Disc').style.display = disc.length ? '' : 'none';
    el('p9Simmap').style.display = 'none';
    el('p9Export').style.display = '';
    if (!cont.length && !disc.length) {
      showMessage(msg, 'error', L2('No se reconoció ninguna columna de caracteres.', 'No character column was recognised.'));
    }
  }
  function fillSelect(id, cols, pick) {
    const s = el(id);
    if (!s) return;
    s.innerHTML = cols.map(c => `<option value="${esc(c.name)}">${esc(c.name)}</option>`).join('');
    if (pick != null && cols[pick]) s.value = cols[pick].name;
  }
  const colByName = nm => (B9.table ? B9.table.columns.find(c => c.name === nm) : null);

  function readTable() {
    const msg = el('p9DataMsg');
    msg.innerHTML = '';
    const t = parseTable(el('p9Text').value);
    if (!t) {
      showMessage(msg, 'error', L2('Hacen falta al menos una línea de encabezado y una de datos.', 'At least a header line and one data line are needed.'));
      return;
    }
    const tree = chosenTree();
    if (!tree) return;
    B9.table = t;
    B9.tree = tree;
    B9.matched = matchToTree(t, tree);
    B9.cont = B9.signal = B9.two = B9.disc = B9.sim = null;
    ['p9ContTiles', 'p9ContTable', 'p9SignalTiles', 'p9TwoTable', 'p9DiscTiles', 'p9DiscTable', 'p9SimTiles'].forEach(id => { if (el(id)) el(id).innerHTML = ''; });
    ['p9ContVerdict', 'p9SignalVerdict', 'p9TwoVerdict', 'p9DiscVerdict', 'p9SimVerdict',
     'p9GramPane', 'p9SignalPane', 'p9ScatterPane', 'p9ContrastPane', 'p9PiePane', 'p9SimPane'].forEach(id => { if (el(id)) el(id).style.display = 'none'; });
    showTable();
    commit();
  }

  /* ================================================================
     2 · one continuous character
     ================================================================ */
  function runCont() {
    const msg = el('p9ContMsg');
    msg.innerHTML = '';
    const col = colByName(el('p9ContVar').value);
    if (!col) return;
    const useLog = el('p9Log').checked;
    if (useLog && col.values.some(v => v != null && v <= 0)) {
      showMessage(msg, 'error', L2('El logaritmo necesita valores positivos, y este carácter tiene ceros o negativos.',
        'The logarithm needs positive values, and this character has zeros or negatives.'));
      return;
    }
    const a = analysisFor(col, { log: useLog });
    if (a.x.length < 5) {
      showMessage(msg, 'error', L2(`Solo ${a.x.length} especies con este carácter: hacen falta al menos cinco.`,
        `Only ${a.x.length} species with this character: at least five are needed.`));
      return;
    }
    const C = Traits.vcv(a.tree, a.x.length);
    const um = Traits.isUltrametric(C);
    const cmp = Traits.compareContinuous(a.tree, a.x);
    const anc = Traits.ancestralContinuous(a.tree, a.x, { C });
    B9.cont = { col: col.name, log: useLog, cmp, anc, a, C, um };
    B9.contVar = col.name;

    const bm = cmp.rows.find(r => r.model === 'BM');
    const ouR = cmp.rows.find(r => r.model === 'OU');
    if (!um.ultrametric) {
      showMessage(msg, 'warn', L2(
        `Las puntas de este árbol no están a la misma distancia de la raíz (varían en ${fmtNum(um.spread)} sobre una profundidad de ${fmtNum(um.depth)}). El movimiento browniano no lo necesita, pero <b>OU y el estallido temprano sí</b>: los dos están escritos en términos del tiempo transcurrido desde la raíz, y sobre un filograma en sustituciones lo que miden no es lo que su nombre dice. Sus filas de la tabla quedan como referencia, no como conclusión.`,
        `The tips of this tree are not all the same distance from the root (they vary by ${fmtNum(um.spread)} over a depth of ${fmtNum(um.depth)}). Brownian motion does not need that, but <b>OU and the early burst do</b>: both are written in terms of the time elapsed since the root, and on a phylogram in substitutions what they measure is not what their names say. Their rows in the table stand as a reference, not as a conclusion.`));
    }
    statTiles('p9ContTiles', [
      [T('Especies usadas', 'Species used'), a.x.length, a.x.length < B9.matched.nTips ? T(`de ${B9.matched.nTips}`, `of ${B9.matched.nTips}`) : ''],
      ['σ²', fmtNum(bm.sigma2), isDated() ? T(`por ${unit()}`, `per ${unit()}`) : T('por sustitución', 'per substitution')],
      [T('Estado en la raíz', 'State at the root'), fmtNum(bm.z0), useLog ? T('en logaritmo', 'on the log scale') : ''],
      [T('Mejor modelo', 'Best model'), cmp.best.model,
        `w = ${fmtFixed(cmp.best.w, 3)}`, cmp.best.w > 0.7 ? 'ok' : ''],
      /* a half-life of six hundred thousand on a tree a tenth of a unit deep is
         not an estimate, it is the lower bound of the search printed out */
      ...(ouR ? [[T('Vida media de OU', 'OU half-life'),
        ouR.atBound ? '—' : fmtNum(ouR.halfLife),
        ouR.atBound ? T('α en el límite: OU colapsó en BM', 'α at its bound: OU collapsed into BM')
          : (isDated() ? unit() : T('sustituciones', 'substitutions')),
        ouR.atBound ? 'bad' : '']] : []),
    ]);

    buildTable('p9ContTable', [
      { key: 'model', label: T('modelo', 'model') },
      { key: 'what', label: T('qué supone', 'what it assumes') },
      { key: 'par', label: T('parámetros', 'parameters') },
      { key: 'lnL', label: 'lnL', num: true },
      { key: 'k', label: 'k', num: true },
      { key: 'aic', label: 'AIC', num: true },
      { key: 'd', label: 'ΔAIC', num: true },
      { key: 'w', label: T('peso', 'weight'), num: true },
    ], cmp.rows.map(r => ({
      model: r.model,
      what: r.model === 'BM' ? T('deriva sin rumbo', 'aimless drift')
        : r.model === 'OU' ? T('atracción a un óptimo', 'pull towards an optimum')
          : T('la tasa decae con el tiempo', 'the rate decays over time'),
      par: r.model === 'BM' ? `σ² = ${fmtNum(r.sigma2)}`
        : r.model === 'OU' ? `σ² = ${fmtNum(r.sigma2)}, α = ${fmtNum(r.alpha)}${r.atBound ? ' *' : ''}`
          : `σ² = ${fmtNum(r.sigma2)}, a = ${fmtNum(r.a)}${r.atBound ? ' *' : ''}`,
      lnL: fmtLnL(r.lnL), k: r.k, aic: fmtFixed(r.AIC, 3),
      d: fmtFixed(r.dAIC, 3), w: fmtFixed(r.w, 3),
    })).map(r => Object.assign(r, { _class: r.model === cmp.best.model ? 'row-best' : '' })));

    const v = el('p9ContVerdict');
    v.style.display = '';
    const dOU = ouR ? ouR.AIC - bm.AIC : 99;
    const ebR = cmp.rows.find(r => r.model === 'EB');
    const partsEs = [], partsEn = [];
    if (cmp.best.model === 'BM' || Math.abs(dOU) < 2) {
      partsEs.push(`Ningún modelo se impone: la diferencia de AIC entre el mejor y el siguiente es ${fmtFixed(cmp.rows[1].dAIC, 2)}, y por debajo de dos unidades no hay con qué elegir. Con ${a.x.length} especies eso es lo normal, y lo honesto es quedarse con el movimiento browniano, que es el más simple.`);
      partsEn.push(`No model wins: the AIC gap between the best and the next is ${fmtFixed(cmp.rows[1].dAIC, 2)}, and below two units there is nothing to choose with. With ${a.x.length} species that is the usual outcome, and the honest thing is to keep Brownian motion, the simplest one.`);
    } else if (cmp.best.model === 'OU') {
      partsEs.push(`OU gana por ${fmtFixed(cmp.rows[1].dAIC, 2)} unidades de AIC. Su vida media es ${fmtNum(ouR.halfLife)} frente a una profundidad del árbol de ${fmtNum(um.depth)}: el carácter olvida la mitad de su historia en ese tiempo. Esto se suele leer como selección estabilizadora, pero un OU también aparece cuando el carácter tiene un límite físico o cuando la muestra no cubre bien el grupo.`);
      partsEn.push(`OU wins by ${fmtFixed(cmp.rows[1].dAIC, 2)} AIC units. Its half-life is ${fmtNum(ouR.halfLife)} against a tree depth of ${fmtNum(um.depth)}: the character forgets half of its history in that time. This is usually read as stabilising selection, but an OU also appears when a character has a physical ceiling or when the sample does not cover the group well.`);
    } else {
      partsEs.push(`El estallido temprano gana por ${fmtFixed(cmp.rows[1].dAIC, 2)} unidades de AIC${ebR && ebR.atBound ? ', pero su parámetro se quedó pegado al borde del intervalo permitido, que es la manera que tiene el ajuste de decir que no hay decaimiento que estimar' : ''}.`);
      partsEn.push(`The early burst wins by ${fmtFixed(cmp.rows[1].dAIC, 2)} AIC units${ebR && ebR.atBound ? ', but its parameter stuck to the edge of the permitted interval, which is how the fit says there is no decay to estimate' : ''}.`);
    }
    if (ebR && ebR.atBound && cmp.best.model !== 'EB') {
      partsEs.push('El parámetro del estallido temprano quedó en el borde (marcado con *): en este árbol no hay desaceleración que medir.');
      partsEn.push('The early burst parameter ended on its bound (marked *): there is no slowdown to measure on this tree.');
    }
    if (ouR && ouR.atBound) {
      partsEs.push('La α de OU también quedó en el borde (*), que es la manera que tiene el ajuste de decir que no ve ninguna atracción: en ese punto OU <i>es</i> el movimiento browniano con un parámetro de más, y por eso su AIC sale dos unidades peor.');
      partsEn.push("OU's α also ended on its bound (*), which is how the fit says it sees no pull at all: at that point OU <i>is</i> Brownian motion with one parameter too many, which is why its AIC comes out two units worse.");
    }
    v.innerHTML = `<b>${L2('Cómo se lee', 'How to read it')}</b> ${L2(partsEs.join(' '), partsEn.join(' '))}`;

    drawTraitgram(a, anc, useLog, col.name);
    el('p9Export').style.display = '';
    commit();
  }

  function drawTraitgram(a, anc, useLog, name) {
    const F = Tree.flatten(a.tree);
    const depth = new Float64Array(F.n);
    for (let i = F.post.length - 1; i >= 0; i--) {
      const k = F.post[i], p = F.parent[k];
      if (p >= 0) depth[k] = depth[p] + (F.len[k] || 0);
    }
    const byNode = new Map(anc.states.map(s => [s.node, s.value]));
    const points = [], edges = [];
    for (let k = 0; k < F.n; k++) {
      const tip = F.isTip[k];
      points.push({
        node: k, tip, tipRow: tip ? F.tipRow[k] : -1, depth: depth[k],
        value: tip ? a.x[F.tipRow[k]] : (byNode.get(k) != null ? byNode.get(k) : 0),
      });
      if (F.parent[k] >= 0) edges.push({ from: F.parent[k], to: k });
    }
    el('p9GramPane').style.display = '';
    el('p9GramFig').innerHTML = Plots9.traitgram(points, edges, {
      labels: a.labels,
      showLabels: a.labels.length <= 40,
      xLabel: isDated() ? T(`tiempo desde la raíz (${unit()})`, `time since the root (${unit()})`)
        : T('sustituciones desde la raíz', 'substitutions since the root'),
      yLabel: useLog ? `log(${name})` : name,
    });
    registerFigure('traitgram', el('p9GramFig').innerHTML);
  }

  /* ================================================================
     3 · phylogenetic signal
     ================================================================ */
  function runSignal() {
    const msg = el('p9SignalMsg');
    msg.innerHTML = '';
    const col = colByName(el('p9ContVar').value);
    if (!col) return;
    const useLog = el('p9Log').checked;
    const a = analysisFor(col, { log: useLog });
    if (a.x.length < 5) {
      showMessage(msg, 'error', L2('Hacen falta al menos cinco especies.', 'At least five species are needed.'));
      return;
    }
    const C = Traits.vcv(a.tree, a.x.length);
    const K = Traits.blombergK(a.tree, a.x, { C });
    const lam = Traits.pagelLambda(a.tree, a.x, { C });
    const reps = Math.max(99, Math.min(9999, +el('p9Reps').value || 999));
    const tst = Traits.signalTest(a.tree, a.x, { C, reps, seed: 17 });
    B9.signal = { col: col.name, log: useLog, K, lam, tst, n: a.x.length };

    statTiles('p9SignalTiles', [
      ['K', fmtFixed(K.K, 4), K.K > 1 ? T('más parecido que browniano', 'more resemblance than Brownian')
        : K.K > 0.6 ? T('cerca de browniano', 'close to Brownian') : T('menos que browniano', 'less than Brownian'),
        tst.p < 0.05 ? 'ok' : 'bad'],
      [T('p de la permutación', 'permutation p'), pEq(tst.p), `${tst.reps} ` + T('permutaciones', 'permutations')],
      ['λ', fmtFixed(lam.lambda, 4), lam.lambda > 0.9 ? T('la estructura del árbol hace falta entera', 'the whole tree structure is needed')
        : lam.lambda > 0.4 ? T('hace falta a medias', 'half of it is needed') : T('casi no hace falta', 'barely needed'),
        lam.p < 0.05 ? 'ok' : 'bad'],
      [T('p de λ (razón de verosimilitud)', 'p for λ (likelihood ratio)'), pEq(lam.p), T('contra λ = 0', 'against λ = 0')],
    ]);

    const v = el('p9SignalVerdict');
    v.style.display = '';
    const sig = tst.p < 0.05 || lam.p < 0.05;
    const es = [], en = [];
    es.push(sig
      ? `Hay señal: repartir los valores al azar entre las especies da una K mayor o igual a ${fmtFixed(K.K, 3)} en ${(tst.p * 100).toFixed(1)} % de las veces, y la prueba de razón de verosimilitud rechaza λ = 0 con ${pEq(lam.p)}. Los parientes se parecen.`
      : `No se detecta señal: ni la permutación de K (${pEq(tst.p)}) ni la razón de verosimilitud de λ (${pEq(lam.p)}) rechazan la hipótesis de que las especies no guarden parecido. Con ${a.x.length} especies eso puede ser falta de potencia tanto como falta de señal.`);
    en.push(sig
      ? `There is signal: shuffling the values among the species gives a K at least ${fmtFixed(K.K, 3)} in ${(tst.p * 100).toFixed(1)} % of the draws, and the likelihood ratio rejects λ = 0 with ${pEq(lam.p)}. Relatives do resemble each other.`
      : `No signal is detected: neither the K permutation (${pEq(tst.p)}) nor the λ likelihood ratio (${pEq(lam.p)}) rejects the hypothesis that species keep no resemblance. With ${a.x.length} species that can be lack of power as much as lack of signal.`);
    if (K.K < 0.8 && lam.lambda > 0.9) {
      es.push(`K = ${fmtFixed(K.K, 3)} y λ = ${fmtFixed(lam.lambda, 3)} parecen contradecirse y no lo hacen: λ dice que la <i>forma</i> del árbol se necesita entera, y K que la <i>cantidad</i> de varianza es menor que la que el browniano predice. Eso es justamente lo que produce un carácter con un óptimo, que no deja que la varianza crezca indefinidamente.`);
      en.push(`K = ${fmtFixed(K.K, 3)} and λ = ${fmtFixed(lam.lambda, 3)} look contradictory and are not: λ says the tree's <i>shape</i> is needed in full, and K that the <i>amount</i> of variance is less than Brownian motion predicts. That is exactly what a character with an optimum produces, since it does not let the variance grow without limit.`);
    }
    if (K.K > 1.2) {
      es.push(`K = ${fmtFixed(K.K, 3)} por encima de uno significa que los parientes se parecen <i>más</i> de lo que el árbol exige. Eso aparece cuando las ramas del árbol están mal estimadas —demasiado largas cerca de las puntas— tanto como cuando hay conservadurismo real.`);
      en.push(`K = ${fmtFixed(K.K, 3)} above one means relatives resemble each other <i>more</i> than the tree requires. That appears when the tree's branches are badly estimated — too long near the tips — as much as when there is real conservatism.`);
    }
    v.innerHTML = `<b>${L2('Cómo se lee', 'How to read it')}</b> ${L2(es.join(' '), en.join(' '))}`;

    el('p9SignalPane').style.display = '';
    el('p9SignalFig').innerHTML = Plots6.histogram(tst.nullDistribution, {
      marks: [{ x: K.K, label: T('K observada', 'observed K') }],
      xLabel: 'K', yLabel: T('permutaciones', 'permutations'),
    });
    registerFigure('signal', el('p9SignalFig').innerHTML);
    commit();
  }

  /* ================================================================
     4 · two characters
     ================================================================ */
  function runTwo() {
    const msg = el('p9TwoMsg');
    msg.innerHTML = '';
    const cy = colByName(el('p9VarY').value), cx = colByName(el('p9VarX').value);
    if (!cy || !cx) return;
    if (cy.name === cx.name) {
      showMessage(msg, 'error', L2('Elige dos caracteres distintos.', 'Pick two different characters.'));
      return;
    }
    const useLog = el('p9LogTwo').checked;
    if (useLog && [cx, cy].some(c => c.values.some(v => v != null && v <= 0))) {
      showMessage(msg, 'error', L2('El logaritmo necesita valores positivos en los dos caracteres.',
        'The logarithm needs positive values in both characters.'));
      return;
    }
    /* only the tips that have both */
    const m = B9.matched;
    const keep = [], xs = [], ys = [];
    m.rowOf.forEach((r, tipRow) => {
      if (r < 0) return;
      const a = cx.values[r], b = cy.values[r];
      if (a == null || b == null) return;
      keep.push(tipRow);
      xs.push(useLog ? Math.log(a) : a);
      ys.push(useLog ? Math.log(b) : b);
    });
    if (keep.length < 6) {
      showMessage(msg, 'error', L2(`Solo ${keep.length} especies tienen los dos caracteres: hacen falta al menos seis.`,
        `Only ${keep.length} species have both characters: at least six are needed.`));
      return;
    }
    const tree = prunedTree(keep);
    const C = Traits.vcv(tree, xs.length);
    const px = Traits.pic(tree, xs), py = Traits.pic(tree, ys);
    let sxy = 0, sxx = 0;
    for (let i = 0; i < px.contrasts.length; i++) { sxy += px.contrasts[i] * py.contrasts[i]; sxx += px.contrasts[i] * px.contrasts[i]; }
    const slopePic = sxy / sxx;
    /* r² of the contrasts regression through the origin */
    let syy = 0;
    py.contrasts.forEach(c => { syy += c * c; });
    const r2pic = (sxy * sxy) / (sxx * syy);
    const dfPic = px.contrasts.length - 1;
    const tPic = slopePic * Math.sqrt(sxx * dfPic / Math.max(1e-300, syy - slopePic * sxy));
    const pPic = 2 * (1 - Traits.studentCdf(Math.abs(tPic), dfPic));

    const gl = Traits.pgls(tree, ys, xs, { C });
    /* and the regression that ignores the tree, for the comparison */
    const n = xs.length;
    const mx = xs.reduce((a, b) => a + b, 0) / n, my = ys.reduce((a, b) => a + b, 0) / n;
    let oxy = 0, oxx = 0, oyy = 0;
    for (let i = 0; i < n; i++) { oxy += (xs[i] - mx) * (ys[i] - my); oxx += (xs[i] - mx) * (xs[i] - mx); oyy += (ys[i] - my) * (ys[i] - my); }
    const slopeOls = oxy / oxx;
    const r2ols = (oxy * oxy) / (oxx * oyy);
    const seOls = Math.sqrt((oyy - slopeOls * oxy) / (n - 2) / oxx);
    const tOls = slopeOls / seOls;
    const pOls = 2 * (1 - Traits.studentCdf(Math.abs(tOls), n - 2));

    B9.two = { x: cx.name, y: cy.name, log: useLog, n, gl, slopePic, r2pic, pPic, slopeOls, r2ols, pOls, px, py, xs, ys, keep };

    buildTable('p9TwoTable', [
      { key: 'meth', label: T('método', 'method') },
      { key: 'slope', label: T('pendiente', 'slope'), num: true },
      { key: 'se', label: T('error estándar', 'standard error'), num: true },
      { key: 't', label: 't', num: true },
      { key: 'p', label: 'p', num: true },
      { key: 'r2', label: 'r²', num: true },
      { key: 'note', label: T('qué supone', 'what it assumes') },
    ], [
      { meth: T('Ordinaria (ignora el árbol)', 'Ordinary (ignores the tree)'), slope: fmtNum(slopeOls), se: fmtNum(seOls),
        t: fmtFixed(tOls, 3), p: pEq(pOls), r2: fmtFixed(r2ols, 4),
        note: T('cada especie un dato independiente — falso', 'every species an independent datum — false') },
      { meth: T('Contrastes independientes', 'Independent contrasts'), slope: fmtNum(slopePic), se: '—',
        t: fmtFixed(tPic, 3), p: pEq(pPic), r2: fmtFixed(r2pic, 4),
        note: T('movimiento browniano, regresión por el origen', 'Brownian motion, regression through the origin') },
      { meth: 'PGLS', slope: fmtNum(gl.beta[1]), se: fmtNum(gl.se[1]),
        t: fmtFixed(gl.t[1], 3), p: pEq(gl.p[1]), r2: '—',
        note: T('lo mismo, con ordenada al origen libre', 'the same, with a free intercept') },
    ]);

    const v = el('p9TwoVerdict');
    v.style.display = '';
    const agree = Math.abs(slopePic - gl.beta[1]) < 1e-8 * Math.max(1, Math.abs(slopePic));
    const shrink = Math.abs(slopeOls) > 1e-12 ? (1 - Math.abs(gl.beta[1] / slopeOls)) * 100 : 0;
    const es = [], en = [];
    es.push(agree
      ? `Las dos correcciones dan exactamente la misma pendiente (${fmtNum(gl.beta[1])}), como tienen que darla: la identidad de Garland e Ives (2000) dice que la regresión por el origen sobre los contrastes y la pendiente de PGLS son el mismo número. Que coincidan hasta el último dígito es la comprobación más barata de que ninguna de las dos está mal calculada.`
      : `Las dos correcciones deberían dar la misma pendiente y difieren en ${fmtNum(Math.abs(slopePic - gl.beta[1]))}; con este árbol y estos datos eso apunta a ramas de longitud cero o a una politomía, que el algoritmo de contrastes resuelve en el orden en que vienen las hijas.`);
    en.push(agree
      ? `Both corrections give exactly the same slope (${fmtNum(gl.beta[1])}), as they must: Garland and Ives's (2000) identity says the through-the-origin regression on contrasts and the PGLS slope are one number. Their agreeing to the last digit is the cheapest check that neither is miscomputed.`
      : `Both corrections should give the same slope and they differ by ${fmtNum(Math.abs(slopePic - gl.beta[1]))}; with this tree and these data that points to zero-length branches or a polytomy, which the contrast algorithm resolves in whatever order the daughters come in.`);
    if (pOls < 0.05 && gl.p[1] >= 0.05) {
      es.push(`La regresión ordinaria es significativa (${pEq(pOls)}) y la corregida no (${pEq(gl.p[1])}). Ese es el caso clásico: la relación aparente la sostenían unos pocos grupos de parientes, y al contar una vez cada historia se desvanece.`);
      en.push(`The ordinary regression is significant (${pEq(pOls)}) and the corrected one is not (${pEq(gl.p[1])}). That is the classic case: the apparent relationship was held up by a few clusters of relatives, and counting each history once makes it vanish.`);
    } else if (pOls >= 0.05 && gl.p[1] < 0.05) {
      es.push(`Aquí pasa lo contrario de lo habitual: la corrección <i>revela</i> la relación (${pEq(gl.p[1])}) en vez de borrarla. Ocurre cuando la estructura filogenética estaba enmascarando la señal, y es un resultado perfectamente legítimo.`);
      en.push(`Here the usual thing is reversed: the correction <i>reveals</i> the relationship (${pEq(gl.p[1])}) instead of erasing it. That happens when the phylogenetic structure was masking the signal, and it is a perfectly legitimate result.`);
    } else if (Math.abs(shrink) > 25) {
      es.push(`La pendiente cambia un ${Math.abs(shrink).toFixed(0)} % al tener en cuenta el árbol, aunque la conclusión sobre la significancia no cambie. Es la magnitud la que hay que informar corregida.`);
      en.push(`The slope changes by ${Math.abs(shrink).toFixed(0)} % once the tree is taken into account, even though the verdict on significance does not change. It is the magnitude that has to be reported corrected.`);
    }
    v.innerHTML = `<b>${L2('Cómo se lee', 'How to read it')}</b> ${L2(es.join(' '), en.join(' '))}`;

    const nx = useLog ? `log(${cx.name})` : cx.name;
    const ny = useLog ? `log(${cy.name})` : cy.name;
    el('p9ScatterPane').style.display = '';
    el('p9ScatterFig').innerHTML = Plots9.scatter(
      xs.map((v2, i) => ({ x: v2, y: ys[i] })),
      [{ slope: slopeOls, intercept: my - slopeOls * mx, label: T('ordinaria', 'ordinary'), colour: 'c10', dashed: true },
       { slope: gl.beta[1], intercept: gl.beta[0], label: 'PGLS', colour: 'c1' }],
      { xLabel: nx, yLabel: ny, width: 480, height: 300 });
    el('p9ContrastPane').style.display = '';
    el('p9ContrastFig').innerHTML = Plots9.scatter(
      px.contrasts.map((v2, i) => ({ x: v2, y: py.contrasts[i] })),
      [{ slope: slopePic, intercept: 0, label: T('por el origen', 'through the origin'), colour: 'c1' }],
      { xLabel: T(`contrastes de ${nx}`, `contrasts of ${nx}`), yLabel: T(`contrastes de ${ny}`, `contrasts of ${ny}`),
        zeroLines: true, width: 480, height: 300 });
    registerFigure('scatter', el('p9ScatterFig').innerHTML);
    registerFigure('contrasts', el('p9ContrastFig').innerHTML);
    el('p9Export').style.display = '';
    commit();
  }

  /* ================================================================
     5 · a discrete character
     ================================================================ */
  function runDisc() {
    const msg = el('p9DiscMsg');
    msg.innerHTML = '';
    const col = colByName(el('p9DiscVar').value);
    if (!col) return;
    const a = analysisFor(col, {});
    if (a.x.length < 5) {
      showMessage(msg, 'error', L2('Hacen falta al menos cinco especies con este carácter.', 'At least five species with this character are needed.'));
      return;
    }
    const levels = col.levels;
    if (levels.length > 8) {
      showMessage(msg, 'error', L2(`${levels.length} estados son demasiados: el modelo ARD necesitaría ${levels.length * (levels.length - 1)} tasas.`,
        `${levels.length} states are too many: the ARD model would need ${levels.length * (levels.length - 1)} rates.`));
      return;
    }
    const cmp = Traits.compareMk(a.tree, a.x, { levels });
    const fit = cmp.best;
    const method = el('p9AncMethod').value;
    const anc = Traits.ancestralDiscrete(fit, method);
    B9.disc = { col: col.name, levels, cmp, fit, anc, method, a };
    B9.discVar = col.name;

    const counts = levels.map(lv => a.x.filter(v => v === lv).length);
    statTiles('p9DiscTiles', [
      [T('Especies usadas', 'Species used'), a.x.length, counts.map((c, i) => `${levels[i]} ${c}`).join(' · ')],
      [T('Mejor modelo', 'Best model'), fit.model, `w = ${fmtFixed(fit.w, 3)}`, fit.w > 0.7 ? 'ok' : ''],
      [T('Tasa', 'Rate'), fit.rates.length === 1 ? fmtNum(fit.rates[0]) : `${fmtNum(Math.min.apply(null, fit.rates))} – ${fmtNum(Math.max.apply(null, fit.rates))}`,
        isDated() ? T(`cambios por ${unit()}`, `changes per ${unit()}`) : T('por sustitución', 'per substitution')],
      ['lnL', fmtLnL(fit.lnL), T('convenio de ape', "ape's convention")],
    ]);

    buildTable('p9DiscTable', [
      { key: 'model', label: T('modelo', 'model') },
      { key: 'what', label: T('qué supone', 'what it assumes') },
      { key: 'r', label: T('tasas', 'rates'), num: true },
      { key: 'lnL', label: 'lnL', num: true },
      { key: 'aic', label: 'AIC', num: true },
      { key: 'd', label: 'ΔAIC', num: true },
      { key: 'w', label: T('peso', 'weight'), num: true },
    ], cmp.rows.map(r => ({
      model: r.model,
      what: r.model === 'ER' ? T('todos los cambios igual de fáciles', 'every change equally easy')
        : r.model === 'SYM' ? T('ida y vuelta iguales, parejas distintas', 'there and back equal, pairs different')
          : T('cada cambio con su propia tasa', 'every change with its own rate'),
      r: r.nRates, lnL: fmtLnL(r.lnL), aic: fmtFixed(r.AIC, 3),
      d: fmtFixed(r.dAIC, 3), w: fmtFixed(r.w, 3),
    })).map(r => Object.assign(r, { _class: r.model === fit.model ? 'row-best' : '' })));

    const v = el('p9DiscVerdict');
    v.style.display = '';
    const root = anc.find(x => x.node === 0);
    const best = root ? root.probs.indexOf(Math.max.apply(null, root.probs)) : -1;
    const doubtful = anc.filter(x => Math.max.apply(null, x.probs) < 0.8).length;
    const es = [], en = [];
    if (root) {
      es.push(`La raíz sale en el estado <b>${levels[best]}</b> con probabilidad ${fmtFixed(root.probs[best], 3)}.`);
      en.push(`The root comes out in state <b>${levels[best]}</b> with probability ${fmtFixed(root.probs[best], 3)}.`);
    }
    es.push(doubtful
      ? `De ${anc.length} nodos internos, ${doubtful} no llegan a 0.8 en ningún estado: son los pasteles repartidos de la figura, y no se deben describir como si tuvieran un estado.`
      : `Los ${anc.length} nodos internos pasan de 0.8 en algún estado, lo que es poco común y suele significar que el carácter cambia despacio comparado con la profundidad del árbol.`);
    en.push(doubtful
      ? `Of ${anc.length} internal nodes, ${doubtful} reach 0.8 in no state: those are the divided pies in the figure, and should not be described as if they had a state.`
      : `All ${anc.length} internal nodes pass 0.8 in some state, which is uncommon and usually means the character changes slowly compared with the depth of the tree.`);
    if (fit.model === 'ARD' && fit.rates.some(r => r < 1e-6)) {
      es.push('Alguna tasa de ARD se estimó en cero. Eso no quiere decir que el cambio sea imposible, sino que en este árbol no se ve ni uno solo, y con una muestra más grande la estimación se movería.');
      en.push('Some ARD rate was estimated at zero. That does not mean the change is impossible, only that not one of them is visible on this tree; with a larger sample the estimate would move.');
    }
    if (method === 'downpass') {
      es.push('Estás viendo el <b>paso descendente</b>: la probabilidad que da el subárbol de abajo, sin lo que dice el resto del árbol. Sirve para entender de dónde sale la reconstrucción, pero lo que se informa es la marginal.');
      en.push('You are seeing the <b>downward pass</b>: the probability the subtree below gives, without what the rest of the tree says. It helps to understand where the reconstruction comes from, but the marginal is what gets reported.');
    }
    v.innerHTML = `<b>${L2('Cómo se lee', 'How to read it')}</b> ${L2(es.join(' '), en.join(' '))}`;

    const probs = new Map(anc.map(x => [x.node, x.probs]));
    el('p9PiePane').style.display = '';
    el('p9PieFig').innerHTML = Plots9.pies(a.tree, probs, {
      labels: a.labels, levels,
      tipState: a.x.map(v => levels.indexOf(v)),
      unit: isDated() ? unit() : T('sust./sitio', 'subst./site'),
      rowHeight: a.labels.length > 45 ? 12 : 16,
    });
    registerFigure('pies', el('p9PieFig').innerHTML);
    el('p9Simmap').style.display = '';
    el('p9Export').style.display = '';
    commit();
  }

  /* ================================================================
     6 · stochastic mapping
     ================================================================ */
  function runSim() {
    const msg = el('p9SimMsg');
    msg.innerHTML = '';
    if (!B9.disc) {
      showMessage(msg, 'error', L2('Primero hay que ajustar el modelo Mk arriba.', 'The Mk model above has to be fitted first.'));
      return;
    }
    const reps = Math.max(50, Math.min(20000, +el('p9SimReps').value || 500));
    el('p9SimProgress').innerHTML = L2(`sorteando ${reps} historias…`, `drawing ${reps} histories…`);
    const w = window.LABG ? LABG.work({
      title: T('Mapeo estocástico de caracteres', 'Stochastic character mapping'),
      message: T(`sorteando ${reps} historias…`, `drawing ${reps} histories…`), delay: 300,
    }) : null;
    phyAfterPaint(() => {
      const t0 = performance.now();
      const sm = Traits.simmap(B9.disc.fit, { reps, seed: 23 });
      const secs = (performance.now() - t0) / 1000;
      el('p9SimProgress').innerHTML = L2(`${reps} en ${secs.toFixed(1)} s`, `${reps} in ${secs.toFixed(1)} s`);
      B9.sim = sm;
      const levels = B9.disc.levels;
      const total = sm.timeInState.reduce((a, b) => a + b, 0);

      statTiles('p9SimTiles', [
        [T('Cambios por historia', 'Changes per history'), fmtFixed(sm.changes.mean, 2), `± ${fmtFixed(sm.changes.sd, 2)}`],
        [T('Longitud del árbol', 'Tree length'), fmtNum(total), isDated() ? unit() : T('sust./sitio', 'subst./site')],
        ...levels.map((lv, i) => [T(`Tiempo en «${lv}»`, `Time in "${lv}"`), fmtNum(sm.timeInState[i]),
          fmtPct(sm.proportion[i])]),
      ]);
      el('p9SimPane').style.display = '';
      el('p9SimFig').innerHTML = Plots9.stateBar(sm.proportion, levels, {
        times: sm.timeInState,
        caption: T(`promedio de ${reps} historias sorteadas`, `average of ${reps} drawn histories`),
      });
      registerFigure('simmap', el('p9SimFig').innerHTML);

      const v = el('p9SimVerdict');
      v.style.display = '';
      const top = sm.proportion.indexOf(Math.max.apply(null, sm.proportion));
      const es = [], en = [];
      es.push(`El árbol pasó el ${fmtPct(sm.proportion[top])} de su longitud en el estado <b>${levels[top]}</b>, y cada historia necesitó ${fmtFixed(sm.changes.mean, 1)} cambios en promedio (desviación ${fmtFixed(sm.changes.sd, 1)}).`);
      en.push(`The tree spent ${fmtPct(sm.proportion[top])} of its length in state <b>${levels[top]}</b>, and each history needed ${fmtFixed(sm.changes.mean, 1)} changes on average (deviation ${fmtFixed(sm.changes.sd, 1)}).`);
      es.push('Estas proporciones son tiempo de rama, no número de especies: un estado puede ocupar poco del árbol y estar en muchas puntas si apareció tarde y varias veces.');
      en.push('These proportions are branch time, not species counts: a state can occupy little of the tree and sit on many tips if it appeared late and more than once.');
      if (sm.changes.sd > sm.changes.mean * 0.5) {
        es.push(`La desviación es grande comparada con la media, lo que quiere decir que los datos no fijan el número de cambios: las historias compatibles con las mismas puntas van de unos pocos a bastantes.`);
        en.push(`The deviation is large compared with the mean, which means the data do not pin the number of changes down: the histories compatible with the same tips range from a few to quite a lot.`);
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
    state.figures['b9_' + name] = { block: 9, name, svg: svgText };
  }
  function commit() {
    state.traits = {
      source: B9.source,
      table: B9.table ? { names: B9.table.names, columns: B9.table.columns.map(c => ({ name: c.name, type: c.type, levels: c.levels })) } : null,
      matched: B9.matched ? { matched: B9.matched.matched, nTips: B9.matched.nTips, missing: B9.matched.missing } : null,
      continuous: B9.cont ? { character: B9.cont.col, log: B9.cont.log, models: B9.cont.cmp.rows.map(r => ({
        model: r.model, lnL: r.lnL, AIC: r.AIC, dAIC: r.dAIC, w: r.w,
        sigma2: r.sigma2, alpha: r.alpha, a: r.a, z0: r.z0, halfLife: r.halfLife })),
        best: B9.cont.cmp.best.model,
        ancestral: B9.cont.anc.states } : null,
      signal: B9.signal ? { character: B9.signal.col, K: B9.signal.K.K, pK: B9.signal.tst.p,
        lambda: B9.signal.lam.lambda, pLambda: B9.signal.lam.p, n: B9.signal.n } : null,
      regression: B9.two ? { x: B9.two.x, y: B9.two.y, log: B9.two.log, n: B9.two.n,
        ols: { slope: B9.two.slopeOls, p: B9.two.pOls, r2: B9.two.r2ols },
        pic: { slope: B9.two.slopePic, p: B9.two.pPic, r2: B9.two.r2pic },
        pgls: { slope: B9.two.gl.beta[1], intercept: B9.two.gl.beta[0], se: B9.two.gl.se[1], p: B9.two.gl.p[1], lnL: B9.two.gl.lnL } } : null,
      discrete: B9.disc ? { character: B9.disc.col, levels: B9.disc.levels,
        models: B9.disc.cmp.rows.map(r => ({ model: r.model, nRates: r.nRates, lnL: r.lnL, AIC: r.AIC, dAIC: r.dAIC, w: r.w, rates: r.rates })),
        best: B9.disc.fit.model, method: B9.disc.method,
        ancestral: B9.disc.anc } : null,
      simmap: B9.sim ? { reps: B9.sim.reps, timeInState: B9.sim.timeInState,
        proportion: B9.sim.proportion, changes: B9.sim.changes } : null,
    };
  }

  function exportAs(kind) {
    const stamp = new Date().toISOString().slice(0, 10);
    if (kind === 'cont' && B9.cont) {
      const rows = [['model', 'lnL', 'k', 'AIC', 'dAIC', 'weight', 'sigma2', 'alpha', 'a', 'z0', 'half_life']];
      B9.cont.cmp.rows.forEach(r => rows.push([r.model, r.lnL, r.k, r.AIC, r.dAIC, r.w,
        r.sigma2, r.alpha == null ? '' : r.alpha, r.a == null ? '' : r.a, r.z0, r.halfLife == null ? '' : r.halfLife]));
      download(rows.map(r => r.join(',')).join('\n'), `phylogenypro_modelos_continuos_${slug(B9.cont.col)}_${stamp}.csv`, 'text/csv');
    } else if (kind === 'anc') {
      const rows = [['character', 'node', 'clade_size', 'value_or_state', 'probability_or_variance']];
      if (B9.cont) {
        const F = Tree.flatten(B9.cont.a.tree);
        B9.cont.anc.states.forEach(s => rows.push([B9.cont.col, s.node, Tree.nTips(F.nodes[s.node]), s.value, Math.sqrt(s.variance)]));
      }
      if (B9.disc) {
        const F = Tree.flatten(B9.disc.a.tree);
        B9.disc.anc.forEach(s => s.probs.forEach((p, i) =>
          rows.push([B9.disc.col, s.node, Tree.nTips(F.nodes[s.node]), B9.disc.levels[i], p])));
      }
      download(rows.map(r => r.join(',')).join('\n'), `phylogenypro_ancestrales_${stamp}.csv`, 'text/csv');
    } else if (kind === 'contrasts' && B9.two) {
      const rows = [[`contrast_${B9.two.x}`, `contrast_${B9.two.y}`, 'variance', 'node']];
      B9.two.px.contrasts.forEach((c, i) => rows.push([c, B9.two.py.contrasts[i], B9.two.px.variances[i], B9.two.px.nodes[i]]));
      download(rows.map(r => r.join(',')).join('\n'), `phylogenypro_contrastes_${stamp}.csv`, 'text/csv');
    } else if (kind === 'mk' && B9.disc) {
      const rows = [['model', 'n_rates', 'lnL_ape', 'lnL_flat_root', 'AIC', 'dAIC', 'weight', 'rates']];
      B9.disc.cmp.rows.forEach(r => rows.push([r.model, r.nRates, r.lnL, r.lnL - Math.log(r.k), r.AIC, r.dAIC, r.w,
        '"' + r.rates.map(v => v.toPrecision(8)).join(' ') + '"']));
      download(rows.map(r => r.join(',')).join('\n'), `phylogenypro_mk_${slug(B9.disc.col)}_${stamp}.csv`, 'text/csv');
    } else if (kind === 'svg') {
      const figs = Object.keys(state.figures || {}).filter(k => k.indexOf('b9_') === 0);
      if (!figs.length) return;
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
    el('p9NoTree').style.display = has ? 'none' : '';
    el('p9Data').style.display = has ? '' : 'none';
    if (!has) {
      ['p9Cont', 'p9Signal', 'p9Two', 'p9Disc', 'p9Simmap', 'p9Export'].forEach(id => { el(id).style.display = 'none'; });
      return;
    }
    const sel = el('p9Source');
    if (!src.some(s => s[0] === B9.source)) B9.source = src[0][0];
    sel.innerHTML = src.map(s => `<option value="${s[0]}">${esc(s[1])}</option>`).join('');
    sel.value = B9.source;
    if (B9.table) { B9.matched = matchToTree(B9.table, chosenTree()); showTable(); }
  }

  function init() {
    if (!el('panel-9')) return;
    el('p9Source').addEventListener('change', () => {
      B9.source = el('p9Source').value;
      if (B9.table) { B9.tree = chosenTree(); B9.matched = matchToTree(B9.table, B9.tree); showTable(); }
    });
    el('p9Upload').addEventListener('click', () => el('p9File').click());
    el('p9File').addEventListener('change', e => {
      const f = e.target.files && e.target.files[0];
      if (!f) return;
      const r = new FileReader();
      r.onload = () => { el('p9Text').value = r.result; readTable(); };
      r.readAsText(f);
    });
    el('p9Example').addEventListener('click', () => {
      el('p9Text').value = (window.EXAMPLES && EXAMPLES.rasgos ? EXAMPLES.rasgos : '').trim();
      readTable();
    });
    el('p9Read').addEventListener('click', readTable);
    el('p9RunCont').addEventListener('click', runCont);
    el('p9RunSignal').addEventListener('click', runSignal);
    el('p9RunTwo').addEventListener('click', runTwo);
    el('p9RunDisc').addEventListener('click', runDisc);
    el('p9AncMethod').addEventListener('change', () => { if (B9.disc) runDisc(); });
    el('p9RunSim').addEventListener('click', runSim);
    els('[data-p9export]').forEach(b => b.addEventListener('click', () => exportAs(b.dataset.p9export)));
    el('p9ToBlock10').addEventListener('click', () => {
      const b = document.querySelector('.step-btn[data-step="10"]');
      if (b && !b.disabled) goStep(10);
      else showMessage(el('p9DataMsg'), 'info', L2(
        'El Bloque 10 (biogeografía histórica) llega en la etapa siguiente.',
        'Block 10 (historical biogeography) arrives in the next stage.'));
    });
    document.addEventListener('stepchange', e => { if (e.detail.step === 9) refresh(); });
    document.addEventListener('langchange', () => {
      if (!sources().length) return;
      refresh();
      if (B9.cont) runCont();
      if (B9.signal) runSignal();
      if (B9.two) runTwo();
      if (B9.disc) runDisc();
    });
    document.addEventListener('themechange', () => { if (B9.cont) runCont(); if (B9.disc) runDisc(); });
  }
  document.addEventListener('DOMContentLoaded', init);

  Object.assign(B9, { refresh, readTable, runCont, runSignal, runTwo, runDisc, runSim, commit, parseTable, matchToTree, prunedTree });
})();
