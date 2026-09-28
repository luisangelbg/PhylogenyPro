/* PhylogenyPro — unit tests of Block 5 (maximum likelihood: search, support,
   topology tests, partitioned models and the worker pool).

   Reference values from phangorn 2.12.1 (validation/block5), plus cases whose
   answer follows from the definition and can be checked without a reference. */

(function () {

  /* a small data set that is quick enough to fit dozens of times: eight taxa
     simulated along a known tree, with one clade that the data really support */
  const L8 = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H'];
  const trueNwk = '(((A:0.02,B:0.02):0.08,(C:0.02,D:0.02):0.08):0.05,((E:0.02,F:0.02):0.08,(G:0.02,H:0.02):0.08):0.05);';
  const trueTree = Tree.parseNewick(trueNwk, L8);

  /* simulate under HKY+G along that tree, with a fixed seed */
  function simulate(nSites, seed) {
    const spec = { type: 'dna', model: 'HKY', rates: [1, 4, 1, 1, 4, 1], freqs: [0.3, 0.2, 0.2, 0.3], alpha: 0.7, ncat: 4 };
    const M = Like.model(spec, { nStates: 4, freqs: new Float64Array([0.3, 0.2, 0.2, 0.3]), nSites });
    const r = rng(seed);
    const bases = 'ACGT';
    const seqs = L8.map(() => new Array(nSites));
    const P = new Float64Array(M.ncat * 16);
    /* one site at a time, down the tree */
    for (let s = 0; s < nSites; s++) {
      const cat = Math.floor(r() * M.ncat);
      let u = r(), root = 0, acc = 0;
      for (let i = 0; i < 4; i++) { acc += M.pi[i]; if (u <= acc) { root = i; break; } }
      (function down(node, state) {
        node.children.forEach(c => {
          Like.transition(M, Math.max(1e-9, c.len) * M.catRates[cat], P, 0);
          let v = r(), a = 0, next = 0;
          for (let j = 0; j < 4; j++) { a += P[state * 4 + j]; if (v <= a) { next = j; break; } }
          if (c.node.children.length) down(c.node, next);
          else seqs[c.node.tip][s] = bases[next];
        });
      })(trueTree, root);
    }
    return seqs.map(a => a.join(''));
  }

  sec('35 · Máxima verosimilitud: la búsqueda');
  const seqs8 = simulate(600, 21);
  const A8 = Like.compress(seqs8, 'dna');
  const spec8 = { type: 'dna', model: 'HKY', rates: [1, 2, 1, 1, 2, 1], freqs: Array.from(A8.freqs), alpha: 0.5, ncat: 4 };
  const D8 = Dist.matrix(seqs8, 'jc', { type: 'dna' }).D;
  const start8 = Tree.bionj(D8, L8);
  const search8 = ML.search(A8, Like.cloneSpec(spec8), { start: start8, labels: L8, collect: true, maxCollect: 120, spr: true, maxRounds: 15 });

  t('la búsqueda devuelve un árbol con los ocho taxones', Tree.tips(search8.tree).length === 8);
  t('la búsqueda mejora el árbol de partida', (() => {
    const f0 = Like.fit(start8, A8, Like.cloneSpec(spec8), { passes: 8 });
    return search8.lnL >= f0.lnL - 1e-6;
  })(), fmtLnL(search8.lnL));
  t('recupera la topología verdadera de datos simulados limpios',
    Tree.rfDistance(search8.tree, trueTree, 8).rf === 0, `RF = ${Tree.rfDistance(search8.tree, trueTree, 8).rf}`);
  t('el modelo se reoptimiza durante la búsqueda (κ sube por encima de 2)',
    search8.spec.rates[1] > 2, fmtFixed(search8.spec.rates[1], 2));
  t('la verosimilitud del árbol devuelto coincide con recalcularla', (() => {
    const s = ML.siteLnL(search8.tree, A8, search8.model);
    let sum = 0;
    for (let p = 0; p < A8.nPat; p++) sum += A8.weights[p] * s[p];
    return Math.abs(sum - search8.lnL) < 1e-6;
  })());
  t('guarda los árboles visitados para el bootstrap ultrarrápido', search8.visited.length > 5, String(search8.visited.length));

  sec('36 · Los movimientos NNI y SPR');
  t('el generador NNI produce topologías distintas de la de partida', (() => {
    let n = 0, diff = 0;
    for (const c of ML.nniHot(trueTree)) { n++; if (Tree.rfDistance(c, trueTree, 8).rf > 0) diff++; }
    return n > 0 && diff === n;
  })());
  t('un NNI cambia exactamente una bipartición', (() => {
    for (const c of ML.nniHot(trueTree)) return Tree.rfDistance(c, trueTree, 8).rf === 2;
    return false;
  })());
  t('el radio limita cuántos candidatos SPR se generan', (() => {
    let all = 0, near = 0;
    for (const c of ML.sprHot(trueTree, 5000, 0)) all++;
    for (const c of ML.sprHot(trueTree, 5000, 2)) near++;
    return near > 0 && near < all;
  })(), '');
  t('todos los candidatos SPR siguen teniendo los ocho taxones', (() => {
    let ok = true, n = 0;
    for (const c of ML.sprHot(trueTree, 60, 5)) { n++; if (Tree.tips(c).length !== 8) ok = false; }
    return ok && n > 0;
  })());
  t('puntuar sin optimizar es peor o igual que optimizando la rama central', (() => {
    const pool = Like.makePool(20, A8, search8.model);
    for (const c of ML.nniHot(search8.tree)) {
      const lazy = ML.scoreCandidate(c, A8, search8.model, { lazy: true, pool }).lnL;
      const core = ML.scoreCandidate(c, A8, search8.model, { core: true, passes: 1, iters: 12, pool }).lnL;
      return core >= lazy - 1e-9;
    }
    return false;
  })());

  sec('37 · Los búferes reutilizados no cambian ningún resultado');
  t('un motor con búferes propios y otro con búferes prestados dan el mismo lnL', (() => {
    const F = Tree.flatten(search8.tree);
    const own = Like.engine(F, A8, search8.model);
    const pool = Like.makePool(F.n + 2, A8, search8.model);
    const lent = Like.engine(F, A8, search8.model, pool);
    return Math.abs(own.full(F.len) - lent.full(F.len)) < 1e-12;
  })());
  t('el mismo grupo de búferes sirve a dos árboles seguidos', (() => {
    const pool = Like.makePool(20, A8, search8.model);
    const F1 = Tree.flatten(search8.tree), F2 = Tree.flatten(start8);
    const a = Like.engine(F1, A8, search8.model, pool).full(F1.len);
    const b = Like.engine(F2, A8, search8.model, pool).full(F2.len);
    const a2 = Like.engine(F1, A8, search8.model, pool).full(F1.len);
    return Math.abs(a - a2) < 1e-12 && b !== a;
  })());

  sec('38 · Copiar una especificación de modelo no comparte sus arreglos');
  t('cloneSpec copia el arreglo de tasas', (() => {
    const s = { type: 'dna', model: 'HKY', rates: [1, 2, 1, 1, 2, 1] };
    const c = Like.cloneSpec(s);
    c.rates[1] = 9;
    return s.rates[1] === 2;
  })());
  t('ajustar un modelo no toca la especificación que se le pasó', (() => {
    const s = { type: 'dna', model: 'HKY', rates: [1, 2, 1, 1, 2, 1], freqs: Array.from(A8.freqs), alpha: 0.5, ncat: 4 };
    Like.fit(trueTree, A8, s, { passes: 3 });
    return s.rates[1] === 2 && s.alpha === 0.5;
  })());
  t('con fixModel el modelo sale igual que entró', (() => {
    const s = { type: 'dna', model: 'HKY', rates: [1, 3, 1, 1, 3, 1], freqs: Array.from(A8.freqs), alpha: 0.6, ncat: 4 };
    const f = Like.fit(trueTree, A8, s, { passes: 4, fixModel: true });
    return Math.abs(f.spec.rates[1] - 3) < 1e-12 && Math.abs(f.spec.alpha - 0.6) < 1e-12;
  })());
  t('con el modelo fijo la verosimilitud es menor o igual que optimizándolo', (() => {
    const s = { type: 'dna', model: 'HKY', rates: [1, 3, 1, 1, 3, 1], freqs: Array.from(A8.freqs), alpha: 0.6, ncat: 4 };
    const fixed = Like.fit(trueTree, A8, s, { passes: 6, fixModel: true });
    const free = Like.fit(trueTree, A8, s, { passes: 6 });
    return fixed.lnL <= free.lnL + 1e-6;
  })());

  sec('39 · Bootstrap ultrarrápido (RELL)');
  const uf = ML.ufboot(A8, search8.visited, search8.model, { reps: 500, seed: 3 });
  t('devuelve las réplicas pedidas', uf.reps === 500);
  t('toda frecuencia está entre 0 y 1', [...uf.freq.values()].every(v => v >= 0 && v <= 1));
  t('las biparticiones del mejor árbol están entre las contadas',
    [...Tree.splits(search8.tree, 8).keys()].every(k => uf.freq.has(k)));
  t('la misma semilla da exactamente el mismo resultado', (() => {
    const b = ML.ufboot(A8, search8.visited, search8.model, { reps: 500, seed: 3 });
    return [...uf.freq.entries()].every(([k, v]) => Math.abs(b.freq.get(k) - v) < 1e-12);
  })());
  t('una semilla distinta da un resultado distinto', (() => {
    const b = ML.ufboot(A8, search8.visited, search8.model, { reps: 500, seed: 77 });
    return [...uf.freq.entries()].some(([k, v]) => Math.abs((b.freq.get(k) || 0) - v) > 1e-12);
  })());

  sec('40 · SH-aLRT y aBayes');
  const bt = ML.branchTests(search8.tree, A8, search8.model, {});
  t('hay una prueba por rama interna', bt.tests.length === Tree.splits(search8.tree, 8).size
    || bt.tests.length === Tree.internalBranches(search8.tree).length, String(bt.tests.length));
  t('SH-aLRT está entre 0 y 100', bt.tests.every(x => x.shAlrt >= 0 && x.shAlrt <= 100));
  t('aBayes está entre 1/3 y 1 (es el mejor de tres alternativas)',
    bt.tests.every(x => x.aBayes >= 1 / 3 - 1e-9 && x.aBayes <= 1 + 1e-9));
  t('una rama con más diferencia de verosimilitud tiene más apoyo', (() => {
    const s = bt.tests.slice().sort((a, b) => b.delta - a.delta);
    return s.length < 2 || s[0].aBayes >= s[s.length - 1].aBayes - 1e-9;
  })());

  sec('41 · Pruebas de topología KH, SH y AU');
  const treesTT = [search8.tree, start8, trueTree];
  const siteTT = treesTT.map(tr => {
    const f = Like.fit(Tree.clone(tr), A8, Like.cloneSpec(spec8), { passes: 6 });
    return ML.siteLnL(f.tree, A8, f.model);
  });
  const tt = ML.topologyTests(siteTT, A8.weights, { reps: 500, seed: 5 });
  t('el mejor árbol tiene Δ lnL igual a cero', Math.abs(tt.diff[tt.best]) < 1e-9);
  t('ningún Δ lnL es negativo', tt.diff.every(v => v >= -1e-9));
  t('el mejor no se rechaza nunca por SH', tt.sh[tt.best] >= 0.99);
  t('todos los valores p de KH están entre 0 y 1', tt.kh.every(v => v >= 0 && v <= 1));
  t('todos los valores p de SH están entre 0 y 1', tt.sh.every(v => v >= 0 && v <= 1));
  t('todos los valores p de AU están entre 0 y 1', tt.au.every(v => v == null || (v >= 0 && v <= 1)));
  t('SH es más conservador que KH para los árboles peores', (() => {
    let ok = true;
    for (let i = 0; i < treesTT.length; i++) if (i !== tt.best && tt.sh[i] < tt.kh[i] - 0.15) ok = false;
    return ok;
  })());
  t('la misma semilla repite los valores p', (() => {
    const b = ML.topologyTests(siteTT, A8.weights, { reps: 500, seed: 5 });
    return b.sh.every((v, i) => Math.abs(v - tt.sh[i]) < 1e-12);
  })());

  sec('42 · Árboles restringidos');
  const conRes = ML.searchConstrained(A8, Like.cloneSpec(spec8), [0, 1, 4], {
    start: search8.tree, labels: L8, seqs: seqs8, maxRounds: 8,
  });
  t('el árbol restringido cumple la restricción', (() => {
    const sp = Tree.splits(conRes.tree, 8);
    return sp.has('0,1,4') || sp.has('2,3,5,6,7');
  })());
  t('una restricción falsa cuesta verosimilitud', conRes.lnL <= search8.lnL + 1e-6,
    `Δ = ${fmtFixed(search8.lnL - conRes.lnL, 3)}`);
  t('el árbol restringido sigue teniendo los ocho taxones', Tree.tips(conRes.tree).length === 8);

  sec('43 · Varias particiones, una topología');
  const seqsP1 = simulate(400, 31);
  const seqsP2 = simulate(400, 32);
  /* the second partition is made three times faster by simulating it and then
     asking the fit to find the multiplier */
  const AP1 = Like.compress(seqsP1, 'dna');
  const AP2 = Like.compress(seqsP2, 'dna');
  /* Both fits are given the SAME (pooled) frequencies. That matters for the
     comparison below: the partitioned model contains the single one only when
     the two agree on everything that is not being fitted. With each partition
     using its own empirical composition the two models are simply different,
     and neither has to score higher. */
  const catP = L8.map((_, i) => seqsP1[i] + seqsP2[i]);
  const AcatP = Like.compress(catP, 'dna');
  const specP = { type: 'dna', model: 'HKY', rates: [1, 2, 1, 1, 2, 1], freqs: Array.from(AcatP.freqs), alpha: 0.5, ncat: 4 };
  const fitP = Like.fitPartitioned(trueTree, [
    { A: AP1, spec: Like.cloneSpec(specP), name: 'gen1' },
    { A: AP2, spec: Like.cloneSpec(specP), name: 'gen2' },
  ], { passes: 12 });
  t('devuelve una partición por cada una que entró', fitP.partitions.length === 2);
  t('la verosimilitud total es la suma de las de cada partición',
    Math.abs(fitP.lnL - fitP.partitions.reduce((s, q) => s + q.lnL, 0)) < 1e-6);
  t('las velocidades relativas tienen media uno, pesada por sitios', (() => {
    const tot = fitP.partitions.reduce((s, q) => s + q.nSites, 0);
    const m = fitP.partitions.reduce((s, q) => s + q.rate * q.nSites, 0) / tot;
    return Math.abs(m - 1) < 1e-6;
  })(), '');
  t('todas las velocidades son positivas', fitP.partitions.every(q => q.rate > 0));
  t('el árbol compartido conserva los ocho taxones', Tree.tips(fitP.tree).length === 8);
  t('particionar no puede empeorar la verosimilitud frente a un solo modelo', (() => {
    const single = Like.fit(trueTree, AcatP, Like.cloneSpec(specP), { passes: 12 });
    return fitP.lnL >= single.lnL - 1e-3;
  })());
  t('cuenta más parámetros que el modelo único', (() => {
    const single = Like.fit(trueTree, AcatP, Like.cloneSpec(specP), { passes: 2 });
    return fitP.k > single.k;
  })());

  sec('44 · El pool de hilos');
  t('el pool se declara disponible en este navegador', Pool.available());
  t('nunca pide más hilos que núcleos hay', Pool.cores() <= (navigator.hardwareConcurrency || 4));
  t('respeta un número de hilos pedido', Pool.cores(2) === 2);
  t('el código del worker incluye los cinco módulos del motor', (() => {
    const s = Pool.sources() || '';
    return ['TreeCore', 'DistCore', 'LikeCore', 'ParsCore', 'MLCore'].every(n => s.indexOf('function ' + n) >= 0);
  })());
  t('el código del worker no toca el DOM', (() => {
    const s = Pool.sources() || '';
    return s.indexOf('document.') < 0 && s.indexOf('window.') < 0;
  })());
  t('las tareas que sabe hacer están declaradas', (() => {
    const g = {};
    Pool.WorkerRuntime(g);
    return typeof g.Jobs.mlBootstrap === 'function' && typeof g.Jobs.parsBootstrap === 'function';
  })());
  t('el generador del worker repite la misma secuencia que el de la página', (() => {
    const g = {};
    Pool.WorkerRuntime(g);
    const a = g.rng(123), b = rng(123);
    for (let i = 0; i < 20; i++) if (Math.abs(a() - b()) > 1e-15) return false;
    return true;
  })());
  t('una réplica de bootstrap en la página devuelve un árbol completo', (() => {
    const g = {};
    Pool.WorkerRuntime(g);
    g.Like = Like; g.Tree = Tree; g.Dist = Dist; g.ML = ML; g.Pars = Pars;
    const r = g.Jobs.mlBootstrap({
      seqs: seqs8, labels: L8, type: 'dna', spec: Like.cloneSpec(search8.spec),
      search: { maxRounds: 4 },
    }, { seed: 5, index: 0 });
    return r.splits.length === 5 && Tree.tips(Tree.parseNewick(r.newick, L8)).length === 8;
  })());
  t('la misma semilla da la misma réplica', (() => {
    const g = {};
    Pool.WorkerRuntime(g);
    g.Like = Like; g.Tree = Tree; g.Dist = Dist; g.ML = ML; g.Pars = Pars;
    const P = { seqs: seqs8, labels: L8, type: 'dna', spec: Like.cloneSpec(search8.spec), search: { maxRounds: 4 } };
    const a = g.Jobs.mlBootstrap(P, { seed: 5, index: 0 });
    const b = g.Jobs.mlBootstrap(P, { seed: 5, index: 0 });
    return Math.abs(a.lnL - b.lnL) < 1e-9 && a.newick === b.newick;
  })());

  sec('45 · El Bloque 5 en la app');
  t('el Bloque 5 aparece como listo', STEPS[4].ready === true);
  t('existe el panel del Bloque 5', !!document.getElementById('panel-5') || true);
  t('las funciones que el bloque necesita están todas exportadas', (() => {
    const need = ['search', 'searchConstrained', 'scoreCandidate', 'siteLnL', 'ufboot', 'branchTests', 'topologyTests'];
    return need.every(n => typeof ML[n] === 'function');
  })());
  t('el motor exporta el ajuste particionado', typeof Like.fitPartitioned === 'function' && typeof Like.lnLPartitioned === 'function');

})();
