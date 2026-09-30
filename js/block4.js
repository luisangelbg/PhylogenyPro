/* PhylogenyPro — molecular phylogenetics in the browser.
   Copyright (C) 2026  Luis Ángel Barrera-Guzmán

   This program is free software: you can redistribute it and/or modify it under
   the terms of the GNU General Public License as published by the Free Software
   Foundation, either version 3 of the License, or (at your option) any later
   version. It is distributed in the hope that it will be useful, but WITHOUT ANY
   WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
   A PARTICULAR PURPOSE. See the GNU General Public License, in the file LICENSE
   at the root of this program, or <https://www.gnu.org/licenses/>. */

/* PhylogenyPro — Block 4: parsimony and distance trees.

   Builds distance trees, searches for the most parsimonious ones, measures
   support by resampling and by Bremer, summarises several trees into a
   consensus, and compares everything it has produced. The trees it leaves in
   state.quick are what Block 5 starts from. */

(function () {

  const B4 = {
    part: -1, A: null, trees: {}, parsRes: null, support: null, bremer: null,
    distInfo: null, run: null, cancelled: false,
  };
  window.B4 = B4;

  function parts() { return (state.data && state.data.parts) || []; }
  function current() { return parts()[B4.part] || null; }

  function prepare() {
    const p = current();
    if (!p) return null;
    if (!B4.A || B4.A.__part !== p.name) {
      B4.A = Like.compress(p.seqs, p.type);
      B4.A.__part = p.name;
      Pars.prepare(B4.A);
    }
    return B4.A;
  }

  /* ================================================================
     distance trees
     ================================================================ */
  function runDistance() {
    const p = current();
    if (!p) return;
    const msg = el('p4DistMsg');
    clearMessages(msg);
    const model = el('p4DistModel').value;
    const opts = {
      type: p.type,
      alpha: el('p4Gamma').checked ? +el('p4Alpha').value : null,
      complete: el('p4Complete').checked,
    };
    let res;
    try {
      res = Dist.matrix(p.seqs, model, opts);
    } catch (e) { showMessage(msg, 'error', esc(e.message)); return; }
    B4.distInfo = res;
    const method = el('p4Method').value;
    let tree, extra = '';
    if (method === 'nj') tree = Tree.nj(res.D, p.taxa);
    else if (method === 'upgma') tree = Tree.upgma(res.D, p.taxa);
    else if (method === 'wpgma') tree = Tree.upgma(res.D, p.taxa, true);
    else if (method === 'me') {
      const m = Tree.meMulti(res.D, p.taxa);
      tree = m.tree;
      extra = T(` · criterio de mínima evolución ${fmtFixed(m.length, 6)} desde ${m.starts} arranques`,
        ` · minimum-evolution criterion ${fmtFixed(m.length, 6)} from ${m.starts} starts`);
    } else tree = Tree.bionj(res.D, p.taxa);
    B4.trees.distance = tree;
    B4.trees.distanceName = `${method.toUpperCase()} · ${model}`;

    statTiles('p4DistTiles', [
      [T('Método', 'Method'), method.toUpperCase(), model],
      [T('Longitud total', 'Total length'), fmtFixed(Tree.totalLength(tree), 4), T('sustituciones/sitio', 'substitutions/site')],
      [T('Distancia media', 'Mean distance'), fmtFixed(meanDist(res.D), 4), T('entre pares', 'between pairs')],
      [T('Distancia máxima', 'Largest distance'), fmtFixed(maxDist(res.D), 4), ''],
      [T('Sitios comparados', 'Sites compared'), res.minSites, T('en el peor par', 'in the worst pair')],
      res.undefinedPairs ? [T('Pares sin corrección', 'Pairs without a correction'), res.undefinedPairs, T('saturados', 'saturated'), 'bad'] : null,
    ].filter(Boolean));

    if (res.undefinedPairs) showMessage(msg, 'warning', L2(
      `En ${res.undefinedPairs} pares la corrección no tiene solución y se sustituyó por un tope. Ese árbol no es de fiar en esas partes: usa un marcador más lento o solo transversiones.`,
      `In ${res.undefinedPairs} pairs the correction has no solution and was replaced by a cap. The tree cannot be trusted there: use a slower marker, or transversions only.`));
    showMessage(msg, 'success', L2(
      `Árbol construido con ${p.taxa.length} taxones.${extra}`, `Tree built with ${p.taxa.length} taxa.${extra}`));
    el('p4DistPane').style.display = '';
    el('p4DistTitle').textContent = `${B4.trees.distanceName}`;
    drawInto('p4DistFig', tree, p.taxa);
    refreshTreeList();
    updateRF();
    commit();
  }
  function meanDist(D) {
    let s = 0, k = 0;
    for (let i = 0; i < D.length; i++) for (let j = i + 1; j < D.length; j++) { s += D[i][j]; k++; }
    return k ? s / k : 0;
  }
  function maxDist(D) {
    let m = 0;
    for (let i = 0; i < D.length; i++) for (let j = i + 1; j < D.length; j++) if (D[i][j] > m) m = D[i][j];
    return m;
  }

  /* ================================================================
     parsimony
     ================================================================ */
  function parsOptions() {
    const p = current();
    const A = prepare();
    const kind = el('p4Cost').value;
    const opts = {};
    if (kind === 'ordered') opts.cost = Pars.costMatrix('ordered', A.nStates);
    else if (kind === 'transversion') opts.cost = Pars.costMatrix('transversion', A.nStates, { tsCost: 1, tvCost: +el('p4TvCost').value });
    return opts;
  }

  function runParsimony() {
    const p = current();
    if (!p) return;
    const A = prepare();
    const msg = el('p4ParsMsg');
    clearMessages(msg);
    const btn = el('p4RunPars'), prog = el('p4ParsProgress');
    btn.disabled = true;
    B4.cancelled = false;
    prog.textContent = T('buscando…', 'searching…');
    const w = window.LABG ? LABG.work({ title: T('Buscando el árbol más parsimonioso', 'Searching for the most parsimonious tree'), delay: 300 }) : null;
    phyAfterPaint(() => {
      try {
        const opts = Object.assign(parsOptions(), {
          starts: +el('p4Starts').value,
          swap: el('p4Swap').value,
          maxTrees: +el('p4MaxTrees').value,
          seed: 17,
          onProgress: (done, total, best) => {
            prog.textContent = `${done}/${total} · ${best} ${T('pasos', 'steps')}`;
            if (w) w.update(done / total, prog.textContent);
          },
          cancelled: () => B4.cancelled,
        });
        const res = Pars.search(A, opts);
        B4.parsRes = res;
        B4.trees.parsimony = res.trees[0];
        B4.trees.parsimonyName = T(`Parsimonia · ${res.steps} pasos`, `Parsimony · ${res.steps} steps`);
        const idx = Pars.indices(res.trees[0], A, parsOptions());
        B4.parsIndices = idx;
        statTiles('p4ParsTiles', [
          [T('Pasos', 'Steps'), res.steps, T('en el árbol más corto', 'on the shortest tree'), 'ok'],
          [T('Árboles igual de cortos', 'Equally short trees'), res.nTrees, res.nTrees >= +el('p4MaxTrees').value ? T('se alcanzó el tope', 'the cap was reached') : ''],
          ['CI', fmtFixed(idx.CI, 3), T('consistencia', 'consistency')],
          ['RI', fmtFixed(idx.RI, 3), T('retención', 'retention')],
          ['RC', fmtFixed(idx.RC, 3), T('consistencia reescalada', 'rescaled consistency')],
          ['HI', fmtFixed(idx.HI, 3), T('homoplasia', 'homoplasy')],
          [T('Reordenamientos', 'Rearrangements'), fmtNum(res.rearrangements, 0), `${(res.ms / 1000).toFixed(1)} s`],
        ]);
        /* what the indices mean */
        const bits = [];
        bits.push(L2(
          `El árbol más corto necesita <b>${res.steps}</b> cambios; el mínimo imaginable para estos caracteres sería ${idx.minSteps} y el máximo ${idx.maxSteps}.`,
          `The shortest tree needs <b>${res.steps}</b> changes; the smallest imaginable for these characters would be ${idx.minSteps} and the largest ${idx.maxSteps}.`));
        bits.push(L2(
          `Un CI de ${fmtFixed(idx.CI, 3)} quiere decir que ${fmtPct(1 - idx.CI)} de los cambios son homoplasia: convergencias y reversiones. En moléculas eso es lo normal, y baja al añadir taxones; lo que importa es el RI (${fmtFixed(idx.RI, 3)}), que mide cuánta de la similitud sí tiene estructura de árbol.`,
          `A CI of ${fmtFixed(idx.CI, 3)} means that ${fmtPct(1 - idx.CI)} of the changes are homoplasy: convergence and reversal. In molecules that is normal, and it falls as taxa are added; what matters is the RI (${fmtFixed(idx.RI, 3)}), which measures how much of the similarity does have tree structure.`));
        if (res.nTrees > 1) bits.push(L2(
          `Hay <b>${res.nTrees}</b> árboles igual de cortos: ninguno es «el» resultado. Lo que se publica es el consenso, más abajo.`,
          `There are <b>${res.nTrees}</b> equally short trees: none of them is "the" answer. What gets published is the consensus, below.`));
        const v = el('p4ParsVerdict');
        v.style.display = '';
        v.innerHTML = `<b>${L2('Cómo se lee', 'How to read it')}</b> ` + bits.join(' ');

        showMessage(msg, 'success', L2(
          `Búsqueda con ${res.swap} desde ${res.starts} arranques: ${res.steps} pasos en ${(res.ms / 1000).toFixed(1)} s.`,
          `Search with ${res.swap} from ${res.starts} starts: ${res.steps} steps in ${(res.ms / 1000).toFixed(1)} s.`));
        el('p4ParsPane').style.display = '';
        el('p4ParsTitle').textContent = B4.trees.parsimonyName;
        drawInto('p4ParsFig', res.trees[0], p.taxa);
        renderConsensus();
        refreshTreeList();
        updateRF();
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
     support
     ================================================================ */
  function runResample() {
    const p = current();
    const A = prepare();
    if (!p || !B4.trees.parsimony) {
      showMessage(el('p4BootMsg'), 'info', L2('Primero busca el árbol más corto.', 'Search for the shortest tree first.'));
      return;
    }
    const msg = el('p4BootMsg');
    clearMessages(msg);
    const btn = el('p4RunBoot'), cancel = el('p4CancelBoot'), prog = el('p4BootProgress');
    btn.disabled = true; cancel.style.display = '';
    B4.cancelled = false;
    const reps = +el('p4Reps').value, kind = el('p4Resample').value;
    /* the pace is measured on the first replicates and reported */
    const t0 = performance.now();
    const w = window.LABG ? LABG.work({
      title: kind === 'bootstrap' ? T('Bootstrap de parsimonia', 'Parsimony bootstrap') : T('Jackknife de parsimonia', 'Parsimony jackknife'),
      delay: 300,
    }) : null;
    phyAfterPaint(() => {
      const res = Pars.resample(A, Object.assign(parsOptions(), {
        reps, kind, seed: 99, starts: 1, repSwap: 'NNI',
        onProgress: (done, total) => {
          const per = (performance.now() - t0) / done;
          const left = per * (total - done) / 1000;
          prog.textContent = `${done}/${total}` + (done >= 3 && left > 2 ? ` · ${T('faltan', 'about')} ${left < 60 ? Math.ceil(left) + ' s' : Math.ceil(left / 60) + ' min'}` : '');
          if (w) w.update(done / total, prog.textContent);
        },
        cancelled: () => B4.cancelled,
      }));
      B4.support = res;
      Pars.applySupport(B4.trees.parsimony, res.freq, p.taxa.length);
      if (B4.trees.distance) Pars.applySupport(B4.trees.distance, res.freq, p.taxa.length);
      renderSupportTable();
      drawInto('p4ParsFig', B4.trees.parsimony, p.taxa);
      redraw();
      const vals = [...res.freq.values()].map(v => v * 100);
      const strong = vals.filter(v => v >= 70).length;
      showMessage(msg, 'success', L2(
        `${res.reps} réplicas de ${kind === 'bootstrap' ? 'bootstrap' : 'jackknife'} en ${(res.ms / 1000).toFixed(1)} s. ` +
          `${strong} de ${B4.trees.parsimony ? Tree.splits(B4.trees.parsimony, p.taxa.length).size : 0} clados del árbol tienen 70% o más.`,
        `${res.reps} ${kind} replicates in ${(res.ms / 1000).toFixed(1)} s. ` +
          `${strong} of ${B4.trees.parsimony ? Tree.splits(B4.trees.parsimony, p.taxa.length).size : 0} clades of the tree reach 70% or more.`));
      btn.disabled = false; cancel.style.display = 'none'; prog.textContent = '';
      commit();
      if (w) w.done();
    }, w);
  }

  function runBremer() {
    const p = current();
    const A = prepare();
    if (!p || !B4.trees.parsimony) return;
    const msg = el('p4BootMsg');
    clearMessages(msg);
    const btn = el('p4RunBremer');
    btn.disabled = true;
    const w = window.LABG ? LABG.work({ title: T('Calculando el soporte de Bremer', 'Computing Bremer support'), delay: 300 }) : null;
    phyAfterPaint(() => {
      const br = Pars.bremer(B4.trees.parsimony, A, Object.assign(parsOptions(), { limit: 8000 }));
      B4.bremer = br;
      renderSupportTable();
      showMessage(msg, 'success', L2(
        `Soporte de Bremer calculado examinando ${br.examined} árboles vecinos.`,
        `Bremer support computed by examining ${br.examined} neighbouring trees.`));
      btn.disabled = false;
      commit();
      if (w) w.done();
    }, w);
  }

  function renderSupportTable() {
    const p = current();
    if (!p || !B4.trees.parsimony) return;
    const n = p.taxa.length;
    const sp = Tree.splits(B4.trees.parsimony, n);
    const rows = [];
    sp.forEach((node, key) => {
      const set = key.split(',').map(Number);
      const names = set.map(i => p.taxa[i]);
      const br = B4.bremer ? (B4.bremer.support.find(s => s.split === key) || {}).bremer : null;
      rows.push({
        clade: names.length <= 4 ? names.join(', ') : `${names.slice(0, 3).join(', ')} … (${names.length})`,
        size: names.length,
        boot: B4.support ? (B4.support.freq.get(key) || 0) * 100 : null,
        bremer: br,
      });
    });
    rows.sort((a, b) => (b.boot == null ? -1 : b.boot) - (a.boot == null ? -1 : a.boot));
    buildTable('p4SupportTable', [
      { key: 'clade', label: T('Clado', 'Clade') },
      { key: 'size', label: T('taxones', 'taxa'), num: true },
      { key: 'boot', label: B4.support && B4.support.kind === 'jackknife' ? 'jackknife %' : 'bootstrap %', num: true, fmt: v => fmtFixed(v, 0) },
      { key: 'bremer', label: T('Bremer', 'Bremer'), num: true, fmt: v => fmtFixed(v, 0) },
    ], rows, { limit: 40 });
  }

  /* ================================================================
     consensus and comparison
     ================================================================ */
  function renderConsensus() {
    const p = current();
    if (!p || !B4.parsRes || B4.parsRes.trees.length < 2) {
      el('p4ConsPane').style.display = 'none';
      return;
    }
    const kind = el('p4ConsType').value;
    const fn = kind === 'strict' ? Consensus.strict : kind === 'semistrict' ? Consensus.semistrict
      : kind === 'greedy' ? Consensus.greedy : Consensus.majority;
    const res = fn(B4.parsRes.trees, { nTaxa: p.taxa.length, labels: p.taxa });
    B4.trees.consensus = res.tree;
    B4.trees.consensusName = T(`Consenso ${kind} de ${B4.parsRes.trees.length} árboles`, `${kind} consensus of ${B4.parsRes.trees.length} trees`);
    el('p4ConsPane').style.display = '';
    el('p4ConsTitle').textContent = B4.trees.consensusName;
    drawInto('p4ConsFig', res.tree, p.taxa, { cladogram: true, supportAs: 'node' });
    refreshTreeList();
  }

  function updateRF() {
    const p = current();
    if (!p) return;
    const list = [];
    if (B4.trees.distance) list.push({ name: B4.trees.distanceName, tree: B4.trees.distance });
    if (B4.trees.parsimony) list.push({ name: B4.trees.parsimonyName, tree: B4.trees.parsimony });
    if (B4.trees.consensus) list.push({ name: B4.trees.consensusName, tree: B4.trees.consensus });
    if (state.models && state.models.tree) list.push({ name: T('Bloque 3 · verosimilitud', 'Block 3 · likelihood'), tree: state.models.tree });
    if (list.length < 2) { el('p4RFTable').innerHTML = ''; return; }
    const rows = [];
    for (let i = 0; i < list.length; i++) for (let j = i + 1; j < list.length; j++) {
      const rf = Tree.rfDistance(list[i].tree, list[j].tree, p.taxa.length);
      rows.push({
        a: list[i].name, b: list[j].name, rf: rf.rf, max: rf.maxRF,
        agree: rf.maxRF ? 1 - rf.rf / rf.maxRF : 1,
      });
    }
    buildTable('p4RFTable', [
      { key: 'a', label: T('Árbol', 'Tree') },
      { key: 'b', label: T('contra', 'against') },
      { key: 'rf', label: T('Robinson–Foulds', 'Robinson–Foulds'), num: true },
      { key: 'max', label: T('máximo posible', 'maximum possible'), num: true },
      { key: 'agree', label: T('coincidencia', 'agreement'), num: true, fmt: v => fmtPct(v, 0) },
    ], rows, {});
    const msg = el('p4ConsMsg');
    clearMessages(msg);
    const worst = rows.reduce((a, b) => (b.agree < a.agree ? b : a), rows[0]);
    if (worst.agree < 0.7) showMessage(msg, 'warning', L2(
      `«${esc(worst.a)}» y «${esc(worst.b)}» solo coinciden en ${fmtPct(worst.agree, 0)} de sus clados. Cuando dos criterios difieren tanto, mira las ramas largas: suele ser atracción de ramas largas en la parsimonia, o distancias saturadas.`,
      `"${esc(worst.a)}" and "${esc(worst.b)}" agree on only ${fmtPct(worst.agree, 0)} of their clades. When two criteria disagree this much, look at the long branches: it is usually long-branch attraction in parsimony, or saturated distances.`));
  }

  /* ================================================================
     drawing
     ================================================================ */
  function drawInto(id, tree, labels, extra) {
    const host = el(id);
    if (!host) return;
    const opts = Object.assign({
      labels, width: 720, layout: 'rect',
      showSupport: true, supportAs: 'node',
      scaleLabel: T(' sustituciones/sitio', ' substitutions/site'),
    }, extra || {});
    host.innerHTML = TreeView.render(tree, opts);
  }

  function refreshTreeList() {
    const sel = el('p4WhichTree');
    if (!sel) return;
    const prev = sel.value;
    const opts = [];
    if (B4.trees.distance) opts.push(['distance', B4.trees.distanceName]);
    if (B4.trees.parsimony) opts.push(['parsimony', B4.trees.parsimonyName]);
    if (B4.trees.consensus) opts.push(['consensus', B4.trees.consensusName]);
    sel.innerHTML = opts.map(([v, n]) => `<option value="${v}">${esc(n)}</option>`).join('');
    if (opts.some(o => o[0] === prev)) sel.value = prev;
    el('p4Figure').style.display = opts.length ? '' : 'none';
    /* the outgroup picker */
    const og = el('p4Outgroup');
    const p = current();
    if (og && p) og.innerHTML = p.taxa.map((t, i) => `<option value="${i}">${esc(t)}</option>`).join('');
    if (opts.length) redraw();
  }

  function chosenTree() {
    const which = el('p4WhichTree').value;
    return B4.trees[which] ? Tree.clone(B4.trees[which]) : null;
  }

  function redraw() {
    const p = current();
    let tree = chosenTree();
    if (!p || !tree) return;
    const rootMode = el('p4Root').value;
    if (rootMode === 'midpoint') tree = Tree.midpointRoot(tree);
    else if (rootMode === 'outgroup') {
      const sel = [...el('p4Outgroup').selectedOptions].map(o => +o.value);
      if (sel.length) {
        const r = Tree.rootByOutgroup(tree, sel);
        tree = r.tree;
        if (!r.monophyletic) showMessage(el('p4ConsMsg'), 'warning', L2(
          'El grupo externo elegido no es monofilético en este árbol, así que no se puede enraizar ahí. Eso ya es un resultado: revisa si esos taxones son realmente el grupo externo.',
          'The chosen outgroup is not monophyletic on this tree, so it cannot be rooted there. That is a result in itself: check whether those taxa really are the outgroup.'));
      }
    }
    Tree.ladderize(tree);
    const supportAs = el('p4SupportAs').value;
    drawInto('p4MainFig', tree, p.taxa, {
      layout: el('p4Layout').value,
      cladogram: el('p4Clado').checked,
      showSupport: supportAs !== 'none',
      supportAs,
      width: el('p4Layout').value === 'circular' ? 720 : 760,
    });
    B4.drawn = tree;
  }

  /* ================================================================
     export and hand-off
     ================================================================ */
  function exportAs(kind) {
    const p = current();
    if (!p) return;
    const tree = B4.drawn || chosenTree();
    const base = slug(p.name) + '_' + el('p4WhichTree').value;
    if (kind === 'newick' && tree) download(Tree.writeNewick(tree, { labels: p.taxa }), base + '.nwk');
    else if (kind === 'nexus' && tree) {
      download(SeqIO.writeNexus(p.taxa, p.seqs, {
        type: p.type,
        trees: [{ name: el('p4WhichTree').value, newick: Tree.writeNewick(tree, { labels: p.taxa }) }],
      }), base + '.nex');
    } else if (kind === 'svg') {
      download(el('p4MainFig').innerHTML, base + '.svg', 'image/svg+xml');
    } else if (kind === 'dist' && B4.distInfo) {
      const header = [''].concat(p.taxa);
      const rows = p.taxa.map((t, i) => [t].concat(Array.from(B4.distInfo.D[i]).map(v => v.toFixed(6))));
      download(matrixToCSV(header, rows), slug(p.name) + '_distancias.csv', 'text/csv');
    } else if (kind === 'alltrees' && B4.parsRes) {
      download(B4.parsRes.trees.map(t => Tree.writeNewick(t, { labels: p.taxa })).join('\n') + '\n',
        slug(p.name) + '_arboles_mp.nwk');
    }
  }

  function commit() {
    const p = current();
    if (!p) return;
    state.quick = {
      partition: p.name,
      distance: B4.trees.distance ? { name: B4.trees.distanceName, newick: Tree.writeNewick(B4.trees.distance, { labels: p.taxa }), tree: B4.trees.distance } : null,
      parsimony: B4.parsRes ? {
        name: B4.trees.parsimonyName, steps: B4.parsRes.steps, nTrees: B4.parsRes.nTrees,
        indices: B4.parsIndices, tree: B4.trees.parsimony,
        newick: Tree.writeNewick(B4.trees.parsimony, { labels: p.taxa }),
      } : null,
      consensus: B4.trees.consensus ? { name: B4.trees.consensusName, tree: B4.trees.consensus } : null,
      support: B4.support ? { kind: B4.support.kind, reps: B4.support.reps, freq: [...B4.support.freq.entries()] } : null,
      bremer: B4.bremer ? B4.bremer.support.map(s => ({ split: s.split, bremer: s.bremer })) : null,
    };
    /* Block 5 needs a starting tree, and now it has one */
    enableStep(5, true);
  }

  /* ================================================================
     wiring
     ================================================================ */
  function refresh() {
    const has = parts().length > 0;
    el('p4NoData').style.display = has ? 'none' : '';
    ['p4Dist', 'p4Pars', 'p4Support', 'p4Consensus'].forEach(id => { el(id).style.display = has ? '' : 'none'; });
    if (!has) { el('p4Figure').style.display = 'none'; return; }
    const sel = el('p4Part');
    const ps = parts();
    sel.innerHTML = ps.map((p, i) => `<option value="${i}">${esc(p.name)} · ${p.taxa.length} × ${p.length}</option>`).join('');
    if (B4.part < 0 || B4.part >= ps.length) B4.part = 0;
    sel.value = String(B4.part);
    /* the distance model that fits the data type */
    const p = current();
    if (p && p.type === 'aa') el('p4DistModel').value = 'poisson';
    else if (p && p.type === 'morph') el('p4DistModel').value = 'p';
    /* the gamma shape found in Block 3, if there is one */
    if (state.models && state.models.best && state.models.best.spec && state.models.best.spec.alpha) {
      el('p4Alpha').value = (+state.models.best.spec.alpha).toFixed(3);
    }
    refreshTreeList();
  }

  function init() {
    if (!el('panel-4')) return;
    el('p4Part').addEventListener('change', () => {
      B4.part = +el('p4Part').value;
      B4.A = null; B4.trees = {}; B4.parsRes = null; B4.support = null; B4.bremer = null;
      ['p4DistPane', 'p4ParsPane', 'p4ConsPane'].forEach(id => { el(id).style.display = 'none'; });
      el('p4Figure').style.display = 'none';
      el('p4ParsVerdict').style.display = 'none';
      refresh();
    });
    el('p4RunDist').addEventListener('click', runDistance);
    el('p4RunPars').addEventListener('click', runParsimony);
    el('p4RunBoot').addEventListener('click', runResample);
    el('p4RunBremer').addEventListener('click', runBremer);
    el('p4CancelBoot').addEventListener('click', () => { B4.cancelled = true; });
    el('p4Cost').addEventListener('change', () => {
      el('p4TvWrap').style.display = el('p4Cost').value === 'transversion' ? '' : 'none';
    });
    el('p4ConsType').addEventListener('change', renderConsensus);
    ['p4WhichTree', 'p4Layout', 'p4Clado', 'p4SupportAs'].forEach(id => el(id).addEventListener('change', redraw));
    el('p4Root').addEventListener('change', () => {
      el('p4OutgroupWrap').style.display = el('p4Root').value === 'outgroup' ? '' : 'none';
      redraw();
    });
    el('p4Outgroup').addEventListener('change', redraw);
    el('p4Redraw').addEventListener('click', redraw);
    els('[data-p4export]').forEach(b => b.addEventListener('click', () => exportAs(b.dataset.p4export)));
    el('p4ToBlock5').addEventListener('click', () => {
      const b = document.querySelector('.step-btn[data-step="5"]');
      if (b && !b.disabled) goStep(5);
      else showMessage(el('p4ConsMsg'), 'info', L2(
        'El Bloque 5 (máxima verosimilitud) llega en la etapa siguiente; estos árboles le servirán de punto de partida.',
        'Block 5 (maximum likelihood) arrives in the next stage; these trees will be its starting point.'));
    });
    document.addEventListener('stepchange', e => { if (e.detail.step === 4) refresh(); });
    document.addEventListener('langchange', () => { if (parts().length) { refresh(); if (B4.parsRes) renderSupportTable(); updateRF(); } });
    document.addEventListener('themechange', () => { if (B4.drawn) redraw(); });
  }
  document.addEventListener('DOMContentLoaded', init);

  Object.assign(B4, { refresh, runDistance, runParsimony, runResample, runBremer, redraw, commit, current });
})();
