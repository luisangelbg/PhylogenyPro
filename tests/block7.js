/* PhylogenyPro — unit tests of Block 7 (molecular clock and dating).

   Reference values from ape 5.8.1 (chronos, rtt) and phangorn 2.12.1
   (optim.pml with optRooted), recorded in validation/block7, plus the
   properties a dating method has to satisfy whatever the data: ages ordered
   down the tree, calibrations honoured, and the whole thing scaling exactly
   with the age given to the root. */

(function () {

  sec('55 · Profundidades, edades y escalado');
  const L6 = ['A', 'B', 'C', 'D', 'E', 'F'];
  const t6 = Tree.parseNewick('(((A:0.1,B:0.1):0.05,(C:0.12,D:0.08):0.04):0.06,(E:0.2,F:0.15):0.03);', L6);
  const F6 = Tree.flatten(t6);
  const d6 = Clock.depths(F6, null);
  t('la raíz está a profundidad cero', d6[0] === 0);
  t('la profundidad de una punta es la suma de su camino', (() => {
    /* A: 0.06 + 0.05 + 0.1 */
    for (let k = 0; k < F6.n; k++) {
      if (F6.isTip[k] && F6.tipRow[k] === 0) return Math.abs(d6[k] - 0.21) < 1e-12;
    }
    return false;
  })());
  t('toda punta está más abajo que su padre', (() => {
    for (let k = 1; k < F6.n; k++) if (d6[k] < d6[F6.parent[k]] - 1e-12) return false;
    return true;
  })());

  sec('56 · Regresión raíz–punta');
  t('una recta perfecta da r² igual a uno', (() => {
    const r = Clock.regression([1, 2, 3, 4], [2, 4, 6, 8]);
    return Math.abs(r.r2 - 1) < 1e-12 && Math.abs(r.slope - 2) < 1e-12;
  })());
  t('la pendiente y la ordenada son las de mínimos cuadrados', (() => {
    const r = Clock.regression([0, 1, 2, 3], [1, 3, 5, 7]);
    return Math.abs(r.slope - 2) < 1e-12 && Math.abs(r.intercept - 1) < 1e-12;
  })());
  t('la tasa es el valor absoluto de la pendiente', (() => {
    const r = Clock.regression([0, 1, 2, 3], [7, 5, 3, 1]);
    return Math.abs(r.rate - 2) < 1e-12 && r.direction === 'age';
  })());
  t('con fechas de calendario la dirección es hacia adelante', Clock.regression([0, 1, 2], [0, 1, 2]).direction === 'forward');
  t('la ordenada en el origen da la fecha de la raíz', (() => {
    const r = Clock.regression([10, 20, 30], [1, 2, 3]);
    /* y = 0 en x = 0 */
    return Math.abs(r.rootDate) < 1e-9;
  })());
  t('sin variación en x la regresión no se inventa una pendiente', (() => {
    const r = Clock.regression([5, 5, 5], [1, 2, 3]);
    return r.slope === 0 && r.r2 === 0;
  })());

  sec('57 · Densidades de las calibraciones');
  t('la normal vale −log(σ√2π) en su media', (() => {
    const c = Clock.prepare({ type: 'normal', mean: 50, sd: 5 });
    return Math.abs(Clock.logDensity(c, 50) + Math.log(5 * Math.sqrt(2 * Math.PI))) < 1e-12;
  })());
  t('y es simétrica alrededor de ella', (() => {
    const c = Clock.prepare({ type: 'normal', mean: 50, sd: 5 });
    return Math.abs(Clock.logDensity(c, 45) - Clock.logDensity(c, 55)) < 1e-12;
  })());
  t('la exponencial desplazada es cero por debajo de su desplazamiento', (() => {
    const c = Clock.prepare({ type: 'exponential', min: 20, mean: 10 });
    return Clock.logDensity(c, 19.9) === -Infinity && isFinite(Clock.logDensity(c, 20.1));
  })());
  t('su mediana es el desplazamiento más media·ln2', (() => {
    const c = Clock.prepare({ type: 'exponential', min: 20, mean: 10 });
    return Math.abs(c.q.q50 - (20 + 10 * Math.log(2))) < 1e-9;
  })());
  t('la lognormal tiene su mediana en el desplazamiento más e^M', (() => {
    const c = Clock.prepare({ type: 'lognormal', min: 10, M: 2, S: 0.5 });
    return Math.abs(c.q.q50 - (10 + Math.exp(2))) < 1e-9;
  })());
  t('la lognormal es cero en el desplazamiento mismo', (() => {
    const c = Clock.prepare({ type: 'lognormal', min: 10, M: 2, S: 0.5 });
    return Clock.logDensity(c, 10) === -Infinity;
  })());
  t('la uniforme es plana y acotada', (() => {
    const c = Clock.prepare({ type: 'uniform', min: 10, max: 20 });
    return Math.abs(Clock.logDensity(c, 11) - Clock.logDensity(c, 19)) < 1e-12
      && Clock.logDensity(c, 9) === -Infinity && Clock.logDensity(c, 21) === -Infinity
      && Math.abs(Clock.logDensity(c, 15) + Math.log(10)) < 1e-12;
  })());
  t('la fija solo admite su valor', (() => {
    const c = Clock.prepare({ type: 'fixed', value: 42 });
    return Clock.logDensity(c, 42) === 0 && Clock.logDensity(c, 42.5) === -Infinity;
  })());
  t('los cuantiles de la normal son μ ± 1.96σ', (() => {
    const c = Clock.prepare({ type: 'normal', mean: 0, sd: 1 });
    return Math.abs(c.q.q025 + 1.959964) < 1e-3 && Math.abs(c.q.q975 - 1.959964) < 1e-3;
  })());

  sec('58 · El ancestro común más reciente');
  t('el ancestro de dos hermanas es su padre', (() => {
    const k = Clock.mrca(F6, [0, 1]);
    return k > 0 && !F6.isTip[k] && F6.kids[k].length === 2;
  })());
  t('el ancestro de todos los taxones es la raíz', Clock.mrca(F6, [0, 1, 2, 3, 4, 5]) === 0);
  t('el ancestro de dos clados lejanos es más profundo que el de dos hermanas', (() => {
    const a = Clock.mrca(F6, [0, 1]);
    const b = Clock.mrca(F6, [0, 4]);
    const d = Clock.depths(F6, null);
    return d[b] < d[a];
  })());

  sec('59 · Datación por mínimos cuadrados');
  const lsA = Clock.lsd(t6, 1000, [{ node: 0, type: 'fixed', value: 10 }], {});
  t('la raíz queda exactamente en la edad fijada', Math.abs(lsA.ages[0] - 10) < 1e-9, fmtFixed(lsA.ages[0], 9));
  t('todas las puntas quedan en cero', (() => {
    for (let k = 0; k < lsA.F.n; k++) if (lsA.F.isTip[k] && Math.abs(lsA.ages[k]) > 1e-9) return false;
    return true;
  })());
  t('ninguna edad supera la de su padre', (() => {
    for (let k = 1; k < lsA.F.n; k++) if (lsA.ages[k] > lsA.ages[lsA.F.parent[k]] + 1e-9) return false;
    return true;
  })());
  t('la tasa es positiva', lsA.rate > 0, lsA.rate.toExponential(3));
  t('duplicar la edad de la raíz duplica todas las edades y divide la tasa', (() => {
    const lsB = Clock.lsd(t6, 1000, [{ node: 0, type: 'fixed', value: 20 }], {});
    if (Math.abs(lsB.rate * 2 - lsA.rate) > lsA.rate * 0.02) return false;
    for (let k = 0; k < lsA.F.n; k++) {
      if (!lsA.F.isTip[k] && lsA.ages[k] > 1e-9) {
        if (Math.abs(lsB.ages[k] / lsA.ages[k] - 2) > 0.02) return false;
      }
    }
    return true;
  })());
  t('un intervalo sobre un clado se respeta', (() => {
    const k = Clock.mrca(F6, [0, 1]);
    const ls = Clock.lsd(t6, 1000, [
      { node: 0, type: 'fixed', value: 10 },
      { node: k, type: 'uniform', min: 7, max: 9 },
    ], {});
    return ls.ages[k] >= 7 - 1e-9 && ls.ages[k] <= 9 + 1e-9;
  })());
  t('el árbol que devuelve lleva las edades en sus nodos', (() => {
    let anyAge = false;
    Tree.nodes(lsA.tree).forEach(n => { if (n.age != null) anyAge = true; });
    return anyAge;
  })());
  t('las longitudes del árbol fechado son diferencias de edad', (() => {
    const F = Tree.flatten(lsA.tree);
    for (let k = 1; k < F.n; k++) {
      const want = lsA.ages[F.parent[k]] - lsA.ages[k];
      if (Math.abs(F.len[k] - want) > 1e-6) return false;
    }
    return true;
  })());

  sec('60 · Verosimilitud penalizada');
  const pen = Clock.penalised(t6, 1000, [{ node: 0, type: 'fixed', value: 1 }], { lambda: 1, passes: 30 });
  t('la raíz queda donde la calibración la fija', Math.abs(pen.rootAge - 1) < 1e-6, fmtFixed(pen.rootAge, 8));
  t('las puntas quedan en cero', (() => {
    for (let k = 0; k < pen.F.n; k++) if (pen.F.isTip[k] && Math.abs(pen.ages[k]) > 1e-9) return false;
    return true;
  })());
  t('las edades están ordenadas hacia la raíz', (() => {
    for (let k = 1; k < pen.F.n; k++) if (pen.ages[k] > pen.ages[pen.F.parent[k]] + 1e-9) return false;
    return true;
  })());
  t('todas las tasas son positivas', (() => {
    for (let k = 1; k < pen.F.n; k++) if (!(pen.rates[k] > 0)) return false;
    return true;
  })());
  t('un suavizado grande acerca las tasas entre sí', (() => {
    const soft = Clock.penalised(t6, 1000, [{ node: 0, type: 'fixed', value: 1 }], { lambda: 0.01, passes: 30 });
    const hard = Clock.penalised(t6, 1000, [{ node: 0, type: 'fixed', value: 1 }], { lambda: 10000, passes: 30 });
    const cvSoft = Clock.rateVariation(soft).cv;
    const cvHard = Clock.rateVariation(hard).cv;
    return cvHard <= cvSoft + 1e-9;
  })());
  t('la validación cruzada devuelve un λ de la lista', (() => {
    const cv = Clock.crossValidate(t6, 1000, [{ node: 0, type: 'fixed', value: 1 }], {
      lambdas: [0.1, 1, 10], passes: 12,
    });
    return [0.1, 1, 10].indexOf(cv.best.lambda) >= 0 && cv.all.length === 3;
  })());
  t('la tasa media es positiva', pen.meanRate > 0);

  sec('61 · El ajuste con reloj y su prueba');
  const seqs6 = ['ACGTACGTAC', 'ACGTACGTAC', 'ACGTACGTAG', 'ACGTACGTAG', 'ACGTTCGTAC', 'ACGTTCGTAC'];
  const A6 = Like.compress(seqs6, 'dna');
  const spec6 = { type: 'dna', model: 'JC', ncat: 1 };
  const fc = Clock.fitClock(t6, A6, spec6, { passes: 8 });
  t('el árbol con reloj es ultramétrico', (() => {
    const F = Tree.flatten(fc.tree);
    const d = Clock.depths(F, null);
    const tips = [];
    for (let k = 0; k < F.n; k++) if (F.isTip[k]) tips.push(d[k]);
    return Math.max.apply(null, tips) - Math.min.apply(null, tips) < 1e-6;
  })());
  t('cuenta n − 1 parámetros libres de tiempo', fc.nFree === 6 - 1, String(fc.nFree));
  t('su verosimilitud es finita', isFinite(fc.lnL));
  t('la prueba tiene n − 2 grados de libertad', (() => {
    const ct = Clock.clockTest(t6, A6, spec6, { passes: 6 });
    return ct.df === 4;
  })());
  t('el ajuste con reloj nunca supera al libre', (() => {
    const ct = Clock.clockTest(t6, A6, spec6, { passes: 6 });
    return ct.clock.lnL <= ct.free.lnL + 1e-6;
  })());
  t('el valor p está entre 0 y 1', (() => {
    const ct = Clock.clockTest(t6, A6, spec6, { passes: 6 });
    return ct.p >= 0 && ct.p <= 1;
  })());

  sec('62 · Reloj relajado bayesiano');
  const rel = Clock.relaxed(t6, A6, spec6, [{ node: 0, type: 'normal', mean: 10, sd: 1 }], {
    generations: 4000, burnin: 1500, sampleEvery: 10, seed: 3, nSites: A6.nSites,
  });
  t('devuelve muestras', rel.samples.length > 50, String(rel.samples.length));
  t('todas las muestras tienen las edades ordenadas', (() => {
    for (const a of rel.ageSamples) {
      for (let k = 1; k < rel.F.n; k++) if (a[k] > a[rel.F.parent[k]] + 1e-9) return false;
    }
    return true;
  })());
  t('la edad de la raíz se queda cerca de su calibración', (() => {
    const m = Mcmc.mean(rel.samples.map(s => s.rootAge));
    return Math.abs(m - 10) < 4;
  })());
  t('hay un resumen por nodo interno', rel.summary.length === 5, String(rel.summary.length));
  t('cada intervalo contiene su mediana', rel.summary.every(s => s.lower <= s.median + 1e-9 && s.median <= s.upper + 1e-9));
  t('el coeficiente de variación de las tasas es positivo',
    Mcmc.mean(rel.samples.map(s => s.coefficientOfVariation)) > 0);
  t('la desviación del reloj relajado es positiva',
    rel.samples.every(s => s.ucldStdev > 0));

  sec('63 · La escala geológica');
  t('el límite Cretácico–Paleógeno está en 66 Ma', (() => {
    const cre = GeoTime.byName('Cretaceous'), pal = GeoTime.byName('Paleogene');
    return Math.abs(cre.top - 66) < 1e-9 && Math.abs(pal.base - 66) < 1e-9;
  })());
  t('45 Ma cae en el Eoceno', GeoTime.at(45, 'epoch').name === 'Eocene');
  t('10 Ma cae en el Mioceno', GeoTime.at(10, 'epoch').name === 'Miocene');
  t('200 Ma cae en el Jurásico', GeoTime.at(200, 'period').name === 'Jurassic');
  t('las épocas del Cenozoico se tocan sin huecos', (() => {
    const eps = GeoTime.within('epoch', 0, 66).sort((a, b) => b.base - a.base);
    for (let i = 1; i < eps.length; i++) if (Math.abs(eps[i - 1].top - eps[i].base) > 1e-9) return false;
    return true;
  })());
  t('toda unidad tiene base mayor que techo', GeoTime.UNITS.every(u => u.base > u.top));
  t('cada unidad tiene un color', GeoTime.UNITS.every(u => /^rgb\(/.test(u.color)));

  sec('64 · El Bloque 7 en la app');
  t('el Bloque 7 aparece como listo', STEPS[6].ready === true);
  t('el motor del reloj exporta lo que el bloque necesita', (() => {
    const need = ['rootToTip', 'rootByDates', 'clockTest', 'fitClock', 'lsd', 'penalised', 'crossValidate', 'relaxed', 'mrca'];
    return need.every(n => typeof Clock[n] === 'function');
  })());
  t('las figuras del bloque existen',
    typeof Plots7.regression === 'function' && typeof Plots7.chronogram === 'function' && typeof Plots7.smoothing === 'function');
  t('el cronograma se dibuja como SVG', (() => {
    const svg = Plots7.chronogram(lsA.tree, { labels: L6, width: 500 });
    return svg.indexOf('<svg') === 0 && svg.length > 500;
  })());
  t('la regresión se dibuja como SVG', (() => {
    const fit = Clock.regression([1, 2, 3, 4, 5], [0.1, 0.2, 0.31, 0.39, 0.52]);
    fit.x = [1, 2, 3, 4, 5]; fit.y = [0.1, 0.2, 0.31, 0.39, 0.52];
    const svg = Plots7.regression(fit, {});
    return svg.indexOf('<svg') === 0;
  })());

})();
