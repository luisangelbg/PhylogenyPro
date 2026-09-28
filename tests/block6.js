/* PhylogenyPro — unit tests of Block 6 (Bayesian inference).

   The sampler is checked against things whose answer is known without running a
   second sampler: its own prior, the analytic moments of the distributions it
   is supposed to reproduce, the symmetry of data that cannot tell two
   topologies apart, and — for the diagnostics — a deterministic series whose
   effective sample size and interval R computes independently with coda
   (validation/block6). */

(function () {

  sec('46 · Priors y densidades');
  const pr = Mcmc.defaultPriors({ type: 'dna', model: 'GTR' });
  t('el prior por omisión de las ramas es exponencial de media 0.1', pr.branchMean === 0.1);
  t('el prior de α es exponencial de media 1', pr.alphaMean === 1);
  t('el log-prior de una rama sigue la fórmula exponencial', (() => {
    const lens = new Float64Array([0, 0.1, 0.2]);
    const mine = Mcmc.logBranchPrior(lens, 3, { branchMean: 0.1 });
    const want = (-Math.log(0.1) - 1) + (-Math.log(0.1) - 2);
    return Math.abs(mine - want) < 1e-12;
  })());
  t('una rama de longitud cero tiene densidad cero', (() => {
    const lens = new Float64Array([0, 0]);
    return Mcmc.logBranchPrior(lens, 2, { branchMean: 0.1 }) === -Infinity;
  })());
  t('el prior del modelo suma las exponenciales de las tasas libres', (() => {
    const spec = { type: 'dna', model: 'HKY', rates: [1, 3, 1, 1, 3, 1] };
    const mine = Mcmc.logModelPrior(spec, { rateMean: 1, alphaMean: 1, pInv: [0, 1], freqsDirichlet: 1 }, Like.DNA_SCHEMES);
    return Math.abs(mine - (-3)) < 1e-12;    // HKY has one free rate: Exp(1) at 3
  })());
  t('una proporción de invariables fuera de rango tiene densidad cero', (() => {
    const spec = { type: 'dna', model: 'JC', pInv: 1.5 };
    return Mcmc.logModelPrior(spec, { rateMean: 1, alphaMean: 1, pInv: [0, 1], freqsDirichlet: 1 }, Like.DNA_SCHEMES) === -Infinity;
  })());

  sec('47 · La densidad de Dirichlet y el sorteo gamma');
  t('la Dirichlet(1,1,1,1) es constante sobre el símplex', (() => {
    const a = [1, 1, 1, 1];
    const d1 = Mcmc.dirichletLogPdf([0.25, 0.25, 0.25, 0.25], a);
    const d2 = Mcmc.dirichletLogPdf([0.7, 0.1, 0.1, 0.1], a);
    return Math.abs(d1 - d2) < 1e-12;
  })());
  t('su constante de normalización es log(3!) = log 6', (() => {
    const d = Mcmc.dirichletLogPdf([0.25, 0.25, 0.25, 0.25], [1, 1, 1, 1]);
    return Math.abs(d - Math.log(6)) < 1e-10;
  })());
  t('lgamma reproduce los factoriales', (() => {
    return Math.abs(Mcmc.lgamma(5) - Math.log(24)) < 1e-9 && Math.abs(Mcmc.lgamma(1) - 0) < 1e-12;
  })());
  t('el sorteo gamma tiene la media y la varianza de su forma', (() => {
    const r = rng(9);
    const n = 40000, k = 3;
    let s = 0, s2 = 0;
    for (let i = 0; i < n; i++) { const v = Mcmc.gammaDraw(r, k); s += v; s2 += v * v; }
    const m = s / n, va = s2 / n - m * m;
    return Math.abs(m - k) < 0.05 && Math.abs(va - k) < 0.12;
  })(), '');
  t('el sorteo normal tiene media cero y varianza uno', (() => {
    const r = rng(4);
    const n = 40000;
    let s = 0, s2 = 0;
    for (let i = 0; i < n; i++) { const v = Mcmc.normalDraw(r); s += v; s2 += v * v; }
    const m = s / n;
    return Math.abs(m) < 0.02 && Math.abs(s2 / n - m * m - 1) < 0.03;
  })());

  sec('48 · El muestreador devuelve su prior (la prueba que delata un Hastings mal puesto)');
  const L4 = ['A', 'B', 'C', 'D'];
  const A4 = Like.compress(['ACGT', 'ACGT', 'ACGT', 'ACGT'], 'dna');
  const start4 = Tree.parseNewick('((A:0.1,B:0.1):0.1,(C:0.1,D:0.1):0.1);', L4);
  const prior4 = Mcmc.run(A4, {
    tree: start4, labels: L4, spec: { type: 'dna', model: 'JC', ncat: 1 }, noLikelihood: true,
    generations: 200000, burnin: 20000, sampleEvery: 40, chains: 1, seed: 7,
    priors: { branchMean: 0.1 },
  });
  const lens4 = [];
  prior4.trees.forEach(nw => {
    const tr = Tree.parseNewick(nw, L4);
    Tree.nodes(tr).forEach(nd => nd.children.forEach(c => lens4.push(c.len)));
  });
  t('la cadena sin verosimilitud da muestras', prior4.samples.length > 1000, String(prior4.samples.length));
  t('todo árbol muestreado no tiene raíz: 2n − 3 = 5 ramas', (() => {
    return prior4.trees.every(nw => {
      const tr = Tree.parseNewick(nw, L4);
      let e = 0;
      Tree.nodes(tr).forEach(nd => { e += nd.children.length; });
      return e === 5;
    });
  })());
  t('la media de las longitudes es la del prior', Math.abs(Mcmc.mean(lens4) - 0.1) < 0.008, fmtFixed(Mcmc.mean(lens4), 5));
  t('su desviación típica también (en una exponencial coincide con la media)',
    Math.abs(Mcmc.sd(lens4) - 0.1) < 0.01, fmtFixed(Mcmc.sd(lens4), 5));
  t('la mediana es media × ln 2', Math.abs(Mcmc.quantile(lens4, 0.5) - 0.0693) < 0.008, fmtFixed(Mcmc.quantile(lens4, 0.5), 5));
  t('la longitud total es una gamma(5, 0.1): media 0.5',
    Math.abs(Mcmc.mean(prior4.samples.map(s => s.treeLength)) - 0.5) < 0.03);
  t('y su desviación típica es 0.1√5 = 0.2236',
    Math.abs(Mcmc.sd(prior4.samples.map(s => s.treeLength)) - 0.2236) < 0.025,
    fmtFixed(Mcmc.sd(prior4.samples.map(s => s.treeLength)), 5));
  t('las tres topologías de cuatro taxones salen igual de probables', (() => {
    const c = new Map();
    prior4.trees.forEach(nw => { const k = Mcmc.topologyKey(nw, L4, 4); c.set(k, (c.get(k) || 0) + 1); });
    if (c.size !== 3) return false;
    const n = prior4.trees.length;
    return [...c.values()].every(v => Math.abs(v / n - 1 / 3) < 0.035);
  })());

  sec('49 · El posterior con datos');
  /* data made symmetric under exchanging C and D: the two topologies that map
     onto each other must get the same posterior probability */
  const rr = rng(5);
  const baseCols = [];
  for (let i = 0; i < 100; i++) baseCols.push([0, 1, 2, 3].map(() => 'ACGT'[Math.floor(rr() * 4)]));
  const symCols = [];
  baseCols.forEach(p => { symCols.push(p); symCols.push([p[0], p[1], p[3], p[2]]); });
  const seqsSym = [0, 1, 2, 3].map(i => symCols.map(c => c[i]).join(''));
  const Asym = Like.compress(seqsSym, 'dna');
  const post = Mcmc.run(Asym, {
    tree: start4, labels: L4, spec: { type: 'dna', model: 'JC', ncat: 1 },
    generations: 150000, burnin: 20000, sampleEvery: 30, chains: 1, seed: 21,
  });
  const cnt = new Map();
  post.trees.forEach(nw => { const k = Mcmc.topologyKey(nw, L4, 4); cnt.set(k, (cnt.get(k) || 0) + 1); });
  const nP = post.trees.length;
  const f = k => (cnt.get(k) || 0) / nP;
  t('con datos simétricos en C y D las dos topologías intercambiables empatan',
    Math.abs(f('0,2') - f('0,3')) < 0.04, `${f('0,2').toFixed(3)} vs ${f('0,3').toFixed(3)}`);
  t('las probabilidades posteriores suman uno', Math.abs(f('0,1') + f('0,2') + f('0,3') - 1) < 1e-9);
  t('ninguna verosimilitud muestreada supera el máximo del mejor árbol', (() => {
    let best = -Infinity;
    ['((A,B),(C,D));', '((A,C),(B,D));', '((A,D),(B,C));'].forEach(s => {
      const fit = Like.fit(Tree.parseNewick(s, L4), Asym, { type: 'dna', model: 'JC', ncat: 1 }, { passes: 12 });
      if (fit.lnL > best) best = fit.lnL;
    });
    return Math.max.apply(null, post.samples.map(s => s.lnL)) <= best + 0.05;
  })());
  t('las frecuencias de bipartición están entre 0 y 1',
    [...post.splitFreq.values()].every(v => v >= 0 && v <= 1));

  sec('50 · Diagnósticos');
  /* the same deterministic AR(1) series R builds in validation/block6 */
  function lcgAR(n, rho, seed) {
    let s = seed >>> 0;
    const u = () => { s = (1103515245 * s + 12345) % 2147483648; return (s + 0.5) / 2147483648; };
    const x = new Float64Array(n);
    let prev = 0;
    for (let i = 0; i < n; i++) {
      const u1 = u(), u2 = u();
      prev = rho * prev + Math.sqrt(-2 * Math.log(u1)) * Math.cos(2 * Math.PI * u2);
      x[i] = prev;
    }
    return x;
  }
  const x1 = lcgAR(20000, 0.9, 12345);
  const x2 = lcgAR(20000, 0.9, 999);
  const x3 = lcgAR(20000, 0.9, 4242);
  t('la serie AR(1) coincide con la que construye R (media)', Math.abs(Mcmc.mean(x1) - 0.095762) < 1e-5, fmtFixed(Mcmc.mean(x1), 6));
  t('y su desviación típica también', Math.abs(Mcmc.sd(x1) - 2.297247) < 1e-5, fmtFixed(Mcmc.sd(x1), 6));
  t('el ESS se acerca al valor teórico n(1−ρ)/(1+ρ) = 1052.6', (() => {
    const e = Mcmc.ess(x1);
    return Math.abs(e - 1052.6) / 1052.6 < 0.25;
  })(), fmtFixed(Mcmc.ess(x1), 1));
  t('el ESS de coda sobre la misma serie es 992.2, del mismo orden',
    Math.abs(Mcmc.ess(x1) - 992.2) / 992.2 < 0.3, fmtFixed(Mcmc.ess(x1), 1));
  t('el HPD coincide con coda::HPDinterval', (() => {
    const h = Mcmc.hpd(x1, 0.95);
    return Math.abs(h.lower - (-4.440432)) < 0.02 && Math.abs(h.upper - 4.511626) < 0.02;
  })(), (() => { const h = Mcmc.hpd(x1, 0.95); return `[${h.lower.toFixed(4)}, ${h.upper.toFixed(4)}]`; })());
  t('el PSRF coincide con coda::gelman.diag (1.000172)',
    Math.abs(Mcmc.psrf([x1, x2, x3]) - 1.000172) < 0.002, fmtFixed(Mcmc.psrf([x1, x2, x3]), 6));
  t('el ESS de ruido sin correlación es casi n', (() => {
    const r = rng(2);
    const w = new Float64Array(5000);
    for (let i = 0; i < w.length; i++) w[i] = Mcmc.normalDraw(r);
    return Mcmc.ess(w) > 0.7 * w.length;
  })());
  t('el ESS nunca pasa del número de muestras', Mcmc.ess(x1) <= x1.length);
  t('el HPD del 95 % es más estrecho que el rango total', (() => {
    const h = Mcmc.hpd(x1, 0.95);
    return (h.upper - h.lower) < (Math.max.apply(null, Array.from(x1)) - Math.min.apply(null, Array.from(x1)));
  })());
  t('el PSRF de una sola cadena no está definido', Mcmc.psrf([x1]) === null);
  t('el PSRF crece cuando las cadenas están desplazadas', (() => {
    const shifted = Array.from(x2).map(v => v + 6);
    return Mcmc.psrf([Array.from(x1), shifted]) > 1.2;
  })());

  sec('51 · ASDSF');
  t('dos corridas idénticas dan ASDSF cero', (() => {
    const a = new Map([['0,1', 0.9], ['2,3', 0.6]]);
    const b = new Map([['0,1', 0.9], ['2,3', 0.6]]);
    return Math.abs(Mcmc.asdsf([a, b], 0.1)) < 1e-12;
  })());
  t('el ASDSF es la media de las desviaciones típicas por bipartición', (() => {
    const a = new Map([['0,1', 1.0]]);
    const b = new Map([['0,1', 0.8]]);
    /* sd of {1.0, 0.8} about their mean, divided by n (population sd) = 0.1 */
    return Math.abs(Mcmc.asdsf([a, b], 0.1) - 0.1) < 1e-12;
  })());
  t('las biparticiones raras se dejan fuera, como hace MrBayes', (() => {
    const a = new Map([['0,1', 0.9], ['4,5', 0.02]]);
    const b = new Map([['0,1', 0.9]]);
    return Math.abs(Mcmc.asdsf([a, b], 0.1)) < 1e-12;
  })());
  t('el ASDSF nunca es negativo', Mcmc.asdsf([new Map([['0,1', 0.2]]), new Map([['0,1', 0.9]])], 0.1) >= 0);

  sec('52 · El registro para Tracer y el resumen del posterior');
  t('el registro tiene una columna por parámetro y una fila por muestra', (() => {
    const log = Mcmc.tracerLog(prior4.samples.slice(0, 5));
    const lines = log.trim().split('\n');
    const head = lines[0].split('\t');
    return lines.length === 6 && head[0] === 'state' && head.indexOf('lnL') > 0 && head.indexOf('treeLength') > 0;
  })());
  t('el resumen encuentra la topología más muestreada y el conjunto creíble', (() => {
    const sum = Mcmc.summarise(post, L4);
    return sum.mapProbability > 0 && sum.mapProbability <= 1
      && sum.credibleSet.length >= 1 && sum.credibleSet.length <= 3
      && sum.nTopologies >= 1;
  })());
  t('las probabilidades del conjunto creíble suman al menos 0.95 o lo agotan', (() => {
    const sum = Mcmc.summarise(post, L4);
    const s = sum.credibleSet.reduce((a, x) => a + x.p, 0);
    return s >= 0.949 || sum.credibleSet.length === sum.nTopologies;
  })());

  sec('53 · Verosimilitud marginal por piedras de paso');
  t('la verosimilitud marginal nunca supera el máximo de la verosimilitud', (() => {
    const ss = Mcmc.steppingStone(Asym, {
      tree: start4, labels: L4, spec: { type: 'dna', model: 'JC', ncat: 1 },
      steps: 8, stepGenerations: 1500, seed: 31,
    });
    let best = -Infinity;
    ['((A,B),(C,D));', '((A,C),(B,D));', '((A,D),(B,C));'].forEach(s => {
      const fit = Like.fit(Tree.parseNewick(s, L4), Asym, { type: 'dna', model: 'JC', ncat: 1 }, { passes: 12 });
      if (fit.lnL > best) best = fit.lnL;
    });
    return ss.logMarginalLikelihood < best;
  })());
  t('la escalera de temperaturas va de 0 a 1', (() => {
    const K = 10;
    const betas = [];
    for (let k = 0; k < K; k++) betas.push(Math.pow(k / (K - 1), 1 / 0.3));
    return Math.abs(betas[0]) < 1e-12 && Math.abs(betas[K - 1] - 1) < 1e-12;
  })());

  sec('54 · El Bloque 6 en la app');
  t('el Bloque 6 aparece como listo', STEPS[5].ready === true);
  t('el muestreador exporta lo que el bloque necesita', (() => {
    const need = ['run', 'steppingStone', 'summarise', 'tracerLog', 'ess', 'hpd', 'psrf', 'asdsf'];
    return need.every(n => typeof Mcmc[n] === 'function');
  })());
  t('las figuras de la traza y del histograma existen',
    typeof Plots6.lines === 'function' && typeof Plots6.histogram === 'function');
  t('el worker sabe correr una cadena y una escalera', (() => {
    const gg = {};
    Pool.WorkerRuntime(gg);
    return typeof gg.Jobs.mcmcRun === 'function' && typeof gg.Jobs.marginalLikelihood === 'function';
  })());
  t('el código del worker incluye el muestreador', (Pool.sources() || '').indexOf('function McmcCore') >= 0);

})();
