/* PhylogenyPro — unit tests of Block 4 (parsimony, minimum evolution, consensus).

   Reference values from phangorn 2.12.1 and ape 5.8.1 (validation/block4), plus
   small cases whose answer can be worked out by hand. */

(function () {

  sec('26 · Parsimonia: Fitch y Sankoff');
  /* five characters on four taxa; the three topologies take 5, 6 and 7 steps
     (phangorn agrees, and it can be checked by eye) */
  const smallSeqs = ['AAAAC', 'AAGAC', 'AGAAT', 'AGGGT'];
  const L4 = ['t1', 't2', 't3', 't4'];
  const As = Like.compress(smallSeqs, 'dna');
  Pars.prepare(As);
  const t12 = Tree.parseNewick('((t1,t2),(t3,t4));', L4);
  const t13 = Tree.parseNewick('((t1,t3),(t2,t4));', L4);
  const t14 = Tree.parseNewick('((t1,t4),(t2,t3));', L4);
  t('((t1,t2),(t3,t4)) necesita 5 pasos (phangorn)', Pars.score(t12, As).steps === 5, String(Pars.score(t12, As).steps));
  t('((t1,t3),(t2,t4)) necesita 6 pasos (phangorn)', Pars.score(t13, As).steps === 6, String(Pars.score(t13, As).steps));
  t('((t1,t4),(t2,t3)) necesita 7 pasos (phangorn)', Pars.score(t14, As).steps === 7, String(Pars.score(t14, As).steps));
  t('un solo carácter constante no cuesta pasos', (() => {
    const B = Like.compress(['A', 'A', 'A', 'A'], 'dna');
    Pars.prepare(B);
    return Pars.score(t12, B).steps === 0;
  })());
  t('un carácter con dos estados repartidos 2+2 cuesta un paso en su topología', (() => {
    const B = Like.compress(['A', 'A', 'G', 'G'], 'dna');
    Pars.prepare(B);
    return Pars.score(t12, B).steps === 1 && Pars.score(t13, B).steps === 2;
  })());
  t('un hueco no cuesta pasos: admite cualquier estado', (() => {
    const B = Like.compress(['A', '-', 'A', 'A'], 'dna');
    Pars.prepare(B);
    return Pars.score(t12, B).steps === 0;
  })());
  t('Sankoff con costes iguales da lo mismo que Fitch', (() => {
    const cost = Pars.costMatrix('equal', As.nStates);
    return Math.abs(Pars.score(t12, As, { cost }).steps - Pars.score(t12, As).steps) < 1e-9;
  })());
  /* the five characters above change only by transitions (A<->G, C<->T), so
     weighting transversions cannot change their cost; a matrix with real
     transversions is needed to see the difference */
  t('Sankoff con transversiones más caras sube el costo', (() => {
    const B = Like.compress(['AAAAC', 'AACAC', 'ACAAT', 'ACCCT'], 'dna');
    Pars.prepare(B);
    const cost = Pars.costMatrix('transversion', B.nStates, { tsCost: 1, tvCost: 5 });
    return Pars.score(t12, B, { cost }).steps > Pars.score(t12, B).steps;
  })());
  t('con transversiones ponderadas, un cambio A<->C cuesta más que uno A<->G', (() => {
    const tv = Like.compress(['A', 'A', 'C', 'C'], 'dna'), ts = Like.compress(['A', 'A', 'G', 'G'], 'dna');
    Pars.prepare(tv); Pars.prepare(ts);
    const cost = Pars.costMatrix('transversion', 4, { tsCost: 1, tvCost: 5 });
    return Pars.score(t12, tv, { cost }).steps === 5 && Pars.score(t12, ts, { cost }).steps === 1;
  })());
  t('la matriz de costes ordenados cuesta |i − j|', (() => {
    const M = Pars.costMatrix('ordered', 4);
    return M[0][3] === 3 && M[1][2] === 1 && M[2][2] === 0;
  })());

  sec('27 · Índices de consistencia y retención');
  const idx4 = Pars.indices(t12, As);
  near('CI de la matriz pequeña (phangorn)', idx4.CI, 0.8, 1e-8);
  near('RI de la matriz pequeña (phangorn)', idx4.RI, 0.66666667, 1e-7);
  t('CI = mínimo/pasos', Math.abs(idx4.CI - idx4.minSteps / idx4.steps) < 1e-12);
  t('RI = (máximo − pasos)/(máximo − mínimo)', Math.abs(idx4.RI - (idx4.maxSteps - idx4.steps) / (idx4.maxSteps - idx4.minSteps)) < 1e-12);
  t('RC = CI × RI', Math.abs(idx4.RC - idx4.CI * idx4.RI) < 1e-12);
  t('HI = 1 − CI', Math.abs(idx4.HI - (1 - idx4.CI)) < 1e-12);
  t('sin homoplasia el CI vale 1', (() => {
    /* two characters that fit ((t1,t2),(t3,t4)) perfectly */
    const B = Like.compress(['AA', 'AA', 'GG', 'GG'], 'dna');
    Pars.prepare(B);
    const ix = Pars.indices(t12, B);
    return Math.abs(ix.CI - 1) < 1e-12;
  })());

  sec('28 · Búsqueda de árboles');
  t('la adición paso a paso construye un árbol con todos los taxones', (() => {
    const tr = Pars.stepwise(As, [0, 1, 2, 3]);
    return Tree.tips(tr).length === 4;
  })());
  t('la búsqueda encuentra la topología más corta del caso pequeño', (() => {
    const res = Pars.search(As, { starts: 2, swap: 'NNI', seed: 1, maxTrees: 5 });
    return res.steps === 5 && Tree.rfDistance(res.trees[0], t12, 4).rf === 0;
  })());
  t('NNI genera los dos reordenamientos de cada rama interna', (() => {
    const tr = Tree.parseNewick('((A,B),(C,D),E);', ['A', 'B', 'C', 'D', 'E']);
    const seen = new Set();
    for (const cand of Tree.nniMoves(tr)) seen.add(Tree.writeNewick(Pars.sortTree(cand), { lengths: false, support: false }));
    return seen.size >= 4;
  })());
  t('SPR genera más vecinos que NNI', (() => {
    const tr = Tree.parseNewick('((A,B),(C,D),E);', ['A', 'B', 'C', 'D', 'E']);
    let nni = 0, spr = 0;
    for (const c of Tree.nniMoves(tr)) nni++;
    for (const c of Tree.sprMoves(tr, false, 200)) spr++;
    return spr > nni;
  })());
  t('todo vecino conserva los mismos taxones', (() => {
    const tr = Tree.parseNewick('((A,B),(C,D),E);', ['A', 'B', 'C', 'D', 'E']);
    for (const c of Tree.sprMoves(tr, false, 60)) {
      const tips = Tree.tips(c).map(x => x.tip).sort().join(',');
      if (tips !== '0,1,2,3,4') return false;
    }
    return true;
  })());
  t('TBR genera más vecinos que SPR', (() => {
    const tr = Tree.parseNewick('(((A,B),(C,D)),(E,F));', ['A', 'B', 'C', 'D', 'E', 'F']);
    let spr = 0, tbr = 0;
    for (const c of Tree.sprMoves(tr, false, 500)) spr++;
    for (const c of Tree.sprMoves(tr, true, 500)) tbr++;
    return tbr >= spr;
  })());

  sec('29 · Mínima evolución balanceada');
  /* an additive matrix: the balanced length must come back exactly */
  t('la longitud BME de un árbol aditivo es la suma de sus ramas', (() => {
    const D = [[0, 0.3, 0.25, 0.35], [0.3, 0, 0.35, 0.45], [0.25, 0.35, 0, 0.3], [0.35, 0.45, 0.3, 0]];
    const tr = Tree.parseNewick('((A:0.1,B:0.2):0.05,(C:0.1,D:0.2));', ['A', 'B', 'C', 'D']);
    return Math.abs(Tree.bmeLength(tr, D) - 0.65) < 1e-9;
  })());
  t('las longitudes balanceadas recuperan las ramas de un árbol aditivo', (() => {
    const D = [[0, 0.3, 0.25, 0.35], [0.3, 0, 0.35, 0.45], [0.25, 0.35, 0, 0.3], [0.35, 0.45, 0.3, 0]];
    const tr = Tree.unroot(Tree.parseNewick('((A:1,B:1):1,(C:1,D:1));', ['A', 'B', 'C', 'D']));
    Tree.setBalancedLengths(tr, D);
    return Math.abs(Tree.totalLength(tr) - 0.65) < 1e-8;
  })());
  t('la búsqueda ME no empeora el árbol de partida', (() => {
    const D = [[0, 0.3, 0.25, 0.35, 0.5], [0.3, 0, 0.35, 0.45, 0.6], [0.25, 0.35, 0, 0.3, 0.5], [0.35, 0.45, 0.3, 0, 0.5], [0.5, 0.6, 0.5, 0.5, 0]];
    const L5 = ['A', 'B', 'C', 'D', 'E'];
    const start = Tree.bionj(D, L5);
    const before = Tree.bmeLength(start, D);
    const after = Tree.me(D, L5, { start }).length;
    return after <= before + 1e-12;
  })());
  t('desenraizar no cambia la longitud total', (() => {
    const tr = Tree.parseNewick('((A:0.1,B:0.2):0.05,(C:0.1,D:0.2):0.05);', ['A', 'B', 'C', 'D']);
    return Math.abs(Tree.totalLength(Tree.unroot(tr)) - Tree.totalLength(tr)) < 1e-12;
  })());
  t('desenraizar deja la raíz con tres ramas', (() => {
    const tr = Tree.unroot(Tree.parseNewick('((A:0.1,B:0.2):0.05,(C:0.1,D:0.2):0.05);', ['A', 'B', 'C', 'D']));
    return tr.children.length === 3;
  })());

  sec('30 · Consensos');
  const L6 = ['A', 'B', 'C', 'D', 'E', 'F'];
  const threeTrees = ['(((A,B),(C,D)),(E,F));', '(((A,B),(C,E)),(D,F));', '(((A,B),(C,D)),(E,F));'].map(s => Tree.parseNewick(s, L6));
  const cs = Consensus.strict(threeTrees, { nTaxa: 6, labels: L6 });
  const cm = Consensus.majority(threeTrees, { nTaxa: 6, labels: L6 });
  t('el consenso estricto sólo conserva lo que está en todos (A+B)', cs.splits.length === 1 && cs.splits[0].set.join(',') === '0,1');
  t('el consenso de mayoría conserva tres divisiones', cm.splits.length === 3);
  t('las frecuencias del consenso son 1, 2/3 y 2/3', (() => {
    const f = cm.splits.map(s => +s.freq.toFixed(4)).sort();
    return f[0] === 0.6667 && f[1] === 0.6667 && f[2] === 1;
  })());
  t('el consenso de árboles idénticos es ese mismo árbol', (() => {
    const t1 = Tree.parseNewick('(((A,B),(C,D)),(E,F));', L6);
    const c = Consensus.strict([t1, Tree.clone(t1), Tree.clone(t1)], { nTaxa: 6, labels: L6 });
    return Tree.rfDistance(c.tree, t1, 6).rf === 0;
  })());
  t('dos divisiones anidadas son compatibles', Consensus.compatible('0,1', '0,1,2', 6));
  t('dos divisiones cruzadas no son compatibles', !Consensus.compatible('0,1', '1,2', 6));
  t('las frecuencias de divisiones suman lo que deben', (() => {
    const f = Consensus.splitFrequencies(threeTrees, 6);
    return Math.abs(f.get('0,1') - 1) < 1e-12 && Math.abs(f.get('2,3') - 2 / 3) < 1e-9;
  })());

  sec('31 · Remuestreo y soporte');
  t('el bootstrap devuelve frecuencias entre 0 y 1', (() => {
    const res = Pars.resample(As, { reps: 10, seed: 3, starts: 1 });
    return [...res.freq.values()].every(v => v >= 0 && v <= 1);
  })());
  t('con datos sin ambigüedad el clado verdadero sale en todas las réplicas', (() => {
    const strong = Like.compress(['AAAAAAAA', 'AAAAAAAA', 'GGGGGGGG', 'GGGGGGGG'], 'dna');
    Pars.prepare(strong);
    const res = Pars.resample(strong, { reps: 10, seed: 5, starts: 1 });
    return (res.freq.get('0,1') || 0) > 0.99;
  })());
  t('el jackknife también corre y da frecuencias', (() => {
    const res = Pars.resample(As, { reps: 10, seed: 3, starts: 1, kind: 'jackknife' });
    return res.kind === 'jackknife' && res.reps === 10;
  })());
  t('los soportes se escriben en los nodos del árbol', (() => {
    const res = Pars.resample(As, { reps: 10, seed: 3, starts: 1 });
    const tr = Pars.applySupport(Tree.clone(t12), res.freq, 4);
    return Tree.nodes(tr).some(n => n.children.length && n.support != null);
  })());
  t('Bremer: perder un clado cuesta cero o más pasos', (() => {
    const br = Pars.bremer(t12, As, { limit: 200 });
    return br.support.every(s => s.bremer == null || s.bremer >= 0);
  })());

  sec('32 · Dibujo de árboles');
  const drawTree = Tree.parseNewick('((A:0.1,B:0.2)90:0.05,(C:0.3,D:0.15)75:0.04,E:0.5);', ['A', 'B', 'C', 'D', 'E']);
  ['rect', 'circular', 'unrooted'].forEach(layout => {
    t(`el dibujo ${layout} produce un SVG válido`, (() => {
      const svg = TreeView.render(drawTree, { layout, labels: ['A', 'B', 'C', 'D', 'E'] });
      return svg.indexOf('<svg') === 0 && svg.length > 200 && !/NaN|undefined/.test(svg);
    })());
  });
  t('el dibujo usa colores del tema', /var\(--/.test(TreeView.render(drawTree, { labels: ['A', 'B', 'C', 'D', 'E'] })));
  t('los soportes aparecen en el dibujo', /90/.test(TreeView.render(drawTree, { labels: ['A', 'B', 'C', 'D', 'E'], showSupport: true })));
  t('el cladograma alinea todas las puntas', (() => {
    const svg = TreeView.render(drawTree, { labels: ['A', 'B', 'C', 'D', 'E'], cladogram: true });
    const xs = [...svg.matchAll(/<circle cx="([\d.]+)"/g)].map(m => +m[1]);
    return xs.length >= 5 && Math.max(...xs) - Math.min(...xs) < 0.6;
  })());
  t('la barra de escala elige un paso redondo', TreeView.niceStep(0.47) === 0.1 && TreeView.niceStep(4.7) === 1);

  sec('33 · El ejemplo de rbcL, contra phangorn');
  const burseraTaxa2 = Object.keys(EXAMPLES.seqs.concatenado);
  const rbcL2 = Object.values(EXAMPLES.seqs.concatenado).map(s => s.slice(0, 621));
  const Ar = Like.compress(rbcL2, 'dna');
  Pars.prepare(Ar);
  const njRef = Tree.parseNewick(window.__refTreeNwk || '((Bursera_grandifolia:0.16,Bursera_linanoe:0.17):0,((((Bursera_cuneata:0.16,Bursera_hindsiana:0.15):0.0098,((Bursera_bipinnata:0.16,Bursera_schlechtendalii:0.136):0.002,Bursera_fagaroides:0.135):0.0205):0.0061,((((Bursera_arborea:0.143,Bursera_excelsa:0.144):0.0067,Commiphora_africana:0.205):0.002,(Bursera_lancifolia:0.162,Boswellia_sacra:0.171):0.014):0.0049,Bursera_copallifera:0.164):0.0137):0.0036,(((Bursera_simaruba:0.141,Bursera_glabrifolia:0.127):0.0098,(Bursera_aptera:0.142,Bursera_morelensis:0.126):0.0216):0.0228,Bursera_palmeri:0.149):0.0153):0.0025);', burseraTaxa2);
  t('el árbol NJ del ejemplo necesita 1940 pasos (phangorn)', Pars.score(njRef, Ar).steps === 1940, String(Pars.score(njRef, Ar).steps));
  t('la búsqueda SPR baja de 1940 pasos', (() => {
    const res = Pars.search(Ar, { starts: 1, swap: 'SPR', seed: 17, maxTrees: 5 });
    return res.steps < 1940;
  })());

})();
