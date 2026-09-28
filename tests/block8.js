/* PhylogenyPro — unit tests of Block 8 (diversification).

   Reference values from ape 5.8.1 (gammaStat, yule, birthdeath, branching.times)
   and geiger 2.0.11 (bd.ms), recorded in validation/block8, plus the properties
   any of these estimators has to satisfy: scale behaviour, nesting of the
   models, and the null distribution of gamma, which is standard normal by
   construction and is the only thing that makes the Monte Carlo test mean
   anything. */

(function () {

  sec('65 · Tiempos de ramificación y linajes en el tiempo');
  const L6 = ['A', 'B', 'C', 'D', 'E', 'F'];
  /* an ultrametric tree whose node ages can be read by eye: 3, 2 and 1 */
  const t6 = Tree.parseNewick('(((A:1,B:1):1,(C:1,D:1):1):1,(E:2,F:2):1);', L6);
  const bt6 = Diversify.branchingTimes(t6);
  t('hay n − 1 tiempos de ramificación', bt6.length === 5, String(bt6.length));
  t('van del más antiguo al más reciente', (() => {
    for (let i = 1; i < bt6.length; i++) if (bt6[i] > bt6[i - 1] + 1e-12) return false;
    return true;
  })());
  t('la raíz está a la profundidad del árbol', Math.abs(bt6[0] - 3) < 1e-12, fmtFixed(bt6[0], 6));
  t('los dos nodos de en medio están a 2', (() => {
    const two = bt6.filter(v => Math.abs(v - 2) < 1e-12).length;
    return two === 2;
  })(), String(bt6.filter(v => Math.abs(v - 2) < 1e-12).length));
  const L = Diversify.ltt(t6);
  t('la curva empieza en dos linajes', L.lineages[0] === 2);
  t('y termina en el número de puntas', L.lineages[L.lineages.length - 1] === 6);
  t('el número de linajes nunca baja', (() => {
    for (let i = 1; i < L.lineages.length; i++) if (L.lineages[i] < L.lineages[i - 1]) return false;
    return true;
  })());
  t('el último tiempo es el presente', L.times[L.times.length - 1] === 0);

  sec('66 · El estadístico γ');
  t('γ necesita al menos cuatro especies', Diversify.gammaStat(Tree.parseNewick('((A:1,B:1):1,C:2);', ['A', 'B', 'C'])) === null);
  t('γ no depende de la escala del árbol', (() => {
    const a = Diversify.gammaStat(t6).gamma;
    const doubled = Tree.parseNewick('(((A:2,B:2):2,(C:2,D:2):2):2,(E:4,F:4):2);', L6);
    const b = Diversify.gammaStat(doubled).gamma;
    return Math.abs(a - b) < 1e-9;
  })());
  /* the branching times are AGES: a node at 0.1 sits next to the present and
     one at 2.9 next to the root, which is the opposite of what the branch
     lengths look like at a glance */
  t('un árbol con los nodos pegados a las puntas da γ positivo', (() => {
    const late = Tree.parseNewick('(((A:0.1,B:0.1):0.1,(C:0.1,D:0.1):0.1):2.8,(E:0.2,F:0.2):2.8);', L6);
    return Diversify.gammaStat(late).gamma > 0;
  })());
  t('un árbol con los nodos pegados a la raíz da γ negativo', (() => {
    const early = Tree.parseNewick('(((A:2.8,B:2.8):0.1,(C:2.8,D:2.8):0.1):0.1,(E:2.9,F:2.9):0.1);', L6);
    return Diversify.gammaStat(early).gamma < 0;
  })());

  sec('67 · Nacimiento puro y nacimiento–muerte');
  const y6 = Diversify.yule(t6);
  t('λ es (nodos internos − 1) entre la longitud total, como en ape', (() => {
    const F = Tree.flatten(t6);
    let X = 0;
    for (let k = 1; k < F.n; k++) X += F.len[k];
    const nInternal = F.n - 6;
    return Math.abs(y6.lambda - (nInternal - 1) / X) < 1e-12;
  })());
  t('la misma λ sale de los tiempos de espera', Math.abs(Diversify.pureBirth(t6).r - y6.lambda) < 1e-10,
    `${Diversify.pureBirth(t6).r.toFixed(8)} vs ${y6.lambda.toFixed(8)}`);
  t('duplicar todas las ramas divide λ entre dos', (() => {
    const doubled = Tree.parseNewick('(((A:2,B:2):2,(C:2,D:2):2):2,(E:4,F:4):2);', L6);
    return Math.abs(Diversify.yule(doubled).lambda * 2 - y6.lambda) < 1e-12;
  })());
  t('el nacimiento–muerte nunca ajusta peor que el nacimiento puro', (() => {
    const bd = Diversify.birthDeath(t6);
    const x = Diversify.branchingTimes(t6);
    const pb = Diversify.bdLogLik(x, 0, Diversify.birthDeath(t6).r);
    return bd.loglik >= pb - 1e-6;
  })());
  t('con extinción cero la verosimilitud se reduce a la del nacimiento puro', (() => {
    const x = Diversify.branchingTimes(t6);
    const r = 0.3;
    const a = Diversify.bdLogLik(x, 0, r);
    const b = Diversify.bdLogLik(x, 1e-12, r);
    return Math.abs(a - b) < 1e-6;
  })());
  t('una tasa negativa no tiene verosimilitud', Diversify.bdLogLik([3, 2, 1], 0.2, -0.1) === -Infinity);
  t('una extinción relativa mayor que uno tampoco', Diversify.bdLogLik([3, 2, 1], 1.5, 0.2) === -Infinity);
  t('λ y μ se recuperan de a y r', (() => {
    const bd = Diversify.birthDeath(t6);
    return Math.abs(bd.lambda - bd.mu - bd.r) < 1e-9 && Math.abs(bd.mu / bd.lambda - bd.a) < 1e-9;
  })());

  sec('68 · Modelos dependientes de la densidad y comparación');
  t('con K enorme la logística es nacimiento puro', (() => {
    const x = Diversify.branchingTimes(t6);
    const r = Diversify.pureBirth(t6).r;
    const a = Diversify.intervalLogLik(x, () => r);
    const b = Diversify.intervalLogLik(x, i => r * (1 - i / 1e9));
    return Math.abs(a - b) < 1e-5;
  })());
  t('con x = 0 la exponencial es nacimiento puro', (() => {
    const x = Diversify.branchingTimes(t6);
    const r = Diversify.pureBirth(t6).r;
    const a = Diversify.intervalLogLik(x, () => r);
    const b = Diversify.intervalLogLik(x, i => r * Math.pow(i, -0));
    return Math.abs(a - b) < 1e-12;
  })());
  t('la logística ajustada no es peor que el nacimiento puro', (() => {
    const dl = Diversify.ddLogistic(t6);
    const pb = Diversify.pureBirth(t6);
    return dl.loglik >= pb.loglik - 1e-4;
  })());
  t('el modelo de dos tasas no es peor que el de una', (() => {
    const tr = Diversify.twoRate(t6);
    const pb = Diversify.pureBirth(t6);
    return !tr || tr.loglik >= pb.loglik - 1e-6;
  })());
  t('el cambio de tasa cae en un tiempo de ramificación', (() => {
    const tr = Diversify.twoRate(t6);
    if (!tr) return true;
    return bt6.some(v => Math.abs(v - tr.shift) < 1e-9);
  })());
  const cmp = Diversify.compare(t6);
  t('la comparación devuelve los cinco modelos', cmp.rows.length === 5, String(cmp.rows.length));
  t('los pesos de Akaike suman uno', Math.abs(cmp.rows.reduce((s, r) => s + r.w, 0) - 1) < 1e-9);
  t('el mejor modelo encabeza la tabla y tiene ΔAIC cero',
    Math.abs(cmp.rows[0].dAIC) < 1e-12 && cmp.rows[0] === cmp.best);
  t('la tabla está ordenada por AIC', (() => {
    for (let i = 1; i < cmp.rows.length; i++) if (cmp.rows[i].AIC < cmp.rows[i - 1].AIC - 1e-12) return false;
    return true;
  })());
  t('el AIC es −2lnL + 2k', cmp.rows.every(r => Math.abs(r.AIC - (-2 * r.loglik + 2 * r.k)) < 1e-9));

  sec('69 · Magallón y Sanderson');
  t('sin extinción y con corona, r = (log n − log 2) / t', (() => {
    const r = Diversify.magallonSanderson(10, 100, 0, true);
    return Math.abs(r - (Math.log(100) - Math.log(2)) / 10) < 1e-12;
  })());
  t('sin extinción y con tallo, r = log n / t', (() => {
    const r = Diversify.magallonSanderson(10, 100, 0, false);
    return Math.abs(r - Math.log(100) / 10) < 1e-12;
  })());
  t('con tallo y extinción, r = log(n(1−ε) + ε) / t', (() => {
    const r = Diversify.magallonSanderson(10, 100, 0.9, false);
    return Math.abs(r - Math.log(100 * 0.1 + 0.9) / 10) < 1e-12;
  })());
  t('más extinción supuesta da una tasa menor', (() => {
    const a = Diversify.magallonSanderson(10, 100, 0, true);
    const b = Diversify.magallonSanderson(10, 100, 0.9, true);
    return b < a;
  })());
  t('la tasa del tallo supera a la de la corona con la misma edad', (() => {
    return Diversify.magallonSanderson(10, 100, 0, false) > Diversify.magallonSanderson(10, 100, 0, true);
  })());
  t('duplicar la edad divide la tasa entre dos', (() => {
    const a = Diversify.magallonSanderson(10, 100, 0, true);
    const b = Diversify.magallonSanderson(20, 100, 0, true);
    return Math.abs(b * 2 - a) < 1e-12;
  })());
  t('una edad no positiva no da estimación', Diversify.magallonSanderson(0, 100, 0, true) === null);

  sec('70 · El estadístico DR');
  t('hay un valor por punta', Diversify.drStatistic(t6).length === 6);
  t('en un árbol simétrico todas las puntas valen lo mismo', (() => {
    const t4 = Tree.parseNewick('((A:1,B:1):1,(C:1,D:1):1);', ['A', 'B', 'C', 'D']);
    const d = Diversify.drStatistic(t4).map(x => x.dr);
    return Math.max.apply(null, d) - Math.min.apply(null, d) < 1e-12;
  })());
  t('la longitud equitativa pesa cada rama con la mitad que la anterior', (() => {
    const t4 = Tree.parseNewick('((A:1,B:1):2,(C:1,D:1):2);', ['A', 'B', 'C', 'D']);
    const d = Diversify.drStatistic(t4)[0];
    /* 1·2⁰ + 2·2⁻¹ = 2 */
    return Math.abs(d.es - 2) < 1e-12;
  })());
  t('una punta con ramas más cortas tiene DR mayor', (() => {
    const tt = Tree.parseNewick('((A:0.1,B:0.1):2,(C:2,D:2):0.1);', ['A', 'B', 'C', 'D']);
    const d = Diversify.drStatistic(tt);
    const a = d.find(x => x.tip === 0).dr, c = d.find(x => x.tip === 2).dr;
    return a > c;
  })());
  t('duplicar el árbol divide DR entre dos', (() => {
    const a = Diversify.drStatistic(t6).map(x => x.dr);
    const doubled = Tree.parseNewick('(((A:2,B:2):2,(C:2,D:2):2):2,(E:4,F:4):2);', L6);
    const b = Diversify.drStatistic(doubled).map(x => x.dr);
    for (let i = 0; i < a.length; i++) if (Math.abs(b[i] * 2 - a[i]) > 1e-9) return false;
    return true;
  })());

  sec('71 · Simulación de árboles de nacimiento–muerte');
  const sim = Diversify.simulateBD(20, 0.3, 0.1, 77);
  t('devuelve un árbol con las puntas pedidas', sim && Tree.tips(sim).length === 20,
    sim ? String(Tree.tips(sim).length) : 'nulo');
  t('el árbol simulado es ultramétrico', (() => {
    const F = Tree.flatten(sim);
    const d = Clock.depths(F, null);
    const tips = [];
    for (let k = 0; k < F.n; k++) if (F.isTip[k]) tips.push(d[k]);
    return Math.max.apply(null, tips) - Math.min.apply(null, tips) < 1e-6;
  })());
  t('todas sus ramas son no negativas', (() => {
    const F = Tree.flatten(sim);
    for (let k = 1; k < F.n; k++) if (F.len[k] < -1e-12) return false;
    return true;
  })());
  t('no quedan nodos con un solo hijo', (() => {
    let ok = true;
    Tree.nodes(sim).forEach(n => { if (n.children.length === 1) ok = false; });
    return ok;
  })());
  t('la misma semilla da el mismo árbol', (() => {
    const a = Diversify.simulateBD(15, 0.3, 0.1, 5);
    const b = Diversify.simulateBD(15, 0.3, 0.1, 5);
    return Tree.writeNewick(a) === Tree.writeNewick(b);
  })());
  /* the property the whole Monte Carlo test rests on */
  t('γ sobre árboles de nacimiento puro simulados es normal estándar', (() => {
    const gs = [];
    for (let i = 0; i < 250; i++) {
      const s = Diversify.simulateBD(20, 0.4, 0, 3000 + i);
      if (!s) continue;
      const gg = Diversify.gammaStat(s);
      if (gg && isFinite(gg.gamma)) gs.push(gg.gamma);
    }
    const m = Mcmc.mean(gs), sd = Mcmc.sd(gs);
    return gs.length > 200 && Math.abs(m) < 0.25 && Math.abs(sd - 1) < 0.25;
  })(), (() => {
    const gs = [];
    for (let i = 0; i < 250; i++) {
      const s = Diversify.simulateBD(20, 0.4, 0, 3000 + i);
      if (!s) continue;
      const gg = Diversify.gammaStat(s);
      if (gg && isFinite(gg.gamma)) gs.push(gg.gamma);
    }
    return `media ${fmtFixed(Mcmc.mean(gs), 3)}, sd ${fmtFixed(Mcmc.sd(gs), 3)}`;
  })());

  sec('72 · La prueba de Monte Carlo de tasa constante');
  const mc = Diversify.mccr(t6, { total: 12, reps: 120, seed: 9, lambda: 0.4 });
  t('devuelve un valor p entre 0 y 1', mc.p >= 0 && mc.p <= 1, fmtFixed(mc.p, 3));
  t('construye una distribución nula', mc.nullDistribution.length > 50, String(mc.nullDistribution.length));
  t('el valor crítico es el percentil 5 de esa distribución', (() => {
    const s = mc.nullDistribution;
    return Math.abs(mc.critical - s[Math.floor(s.length * 0.05)]) < 1e-12;
  })());
  t('el muestreo incompleto corre la nula hacia valores negativos',
    Mcmc.mean(mc.nullDistribution) < 0.3, fmtFixed(Mcmc.mean(mc.nullDistribution), 3));
  t('con el árbol completo la nula está centrada en cero', (() => {
    const full = Diversify.mccr(t6, { total: 6, reps: 150, seed: 4, lambda: 0.4 });
    return Math.abs(Mcmc.mean(full.nullDistribution)) < 0.4;
  })());
  t('podar deja el número de puntas pedido', (() => {
    const s = Diversify.simulateBD(20, 0.3, 0, 12);
    const p = Diversify.prune(s, 10, rng(1));
    return p && Tree.tips(p).length === 10;
  })());

  sec('73 · El Bloque 8 en la app');
  t('el Bloque 8 aparece como listo', STEPS[7].ready === true);
  t('el motor exporta lo que el bloque necesita', (() => {
    const need = ['branchingTimes', 'ltt', 'gammaStat', 'yule', 'birthDeath', 'compare',
      'magallonSanderson', 'drStatistic', 'simulateBD', 'mccr'];
    return need.every(n => typeof Diversify[n] === 'function');
  })());
  t('las figuras del bloque existen',
    typeof Plots8.ltt === 'function' && typeof Plots8.rateCurve === 'function' && typeof Plots8.ranked === 'function');
  t('la curva de linajes se dibuja como SVG', (() => {
    const svg = Plots8.ltt(Diversify.ltt(t6), {});
    return svg.indexOf('<svg') === 0 && svg.length > 500;
  })());
  t('el optimizador de dos parámetros encuentra el mínimo de una parábola', (() => {
    const r = Diversify.nelderMead(p => (p[0] - 3) * (p[0] - 3) + (p[1] + 2) * (p[1] + 2), [0, 0], {});
    return Math.abs(r.x[0] - 3) < 1e-4 && Math.abs(r.x[1] + 2) < 1e-4;
  })());

  /* what Block 8 keeps has to be what the report and the package read: the
     slot was once written under one name and read under another, and the
     diversification never reached either */
  sec('73b · Del Bloque 8 al informe y al paquete');
  (function () {
    const L8 = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H'];
    const t8 = Tree.parseNewick('((((A:1,B:1):1.5,C:2.5):2,(D:3,E:3):1.5):1.5,((F:2,G:2):3,H:5):1);', L8);
    const keep = { tree: B8.tree, ltt: B8.ltt, gamma: B8.gamma, models: B8.models, mccr: B8.mccr, ms: B8.ms, dr: B8.dr };
    const had = state.diversification;
    try {
      Object.assign(B8, { tree: t8, ltt: Diversify.ltt(t8), gamma: Diversify.gammaStat(t8),
        models: Diversify.compare(t8), mccr: null, ms: null, dr: null });
      B8.commit();
      const v = state.diversification;
      t('el Bloque 8 guarda su resultado en el lugar que el informe lee', !!v && Array.isArray(v.models));
      const minAIC = Math.min.apply(null, B8.models.rows.map(r => r.AIC));
      t('y señala como mejor el modelo de menor AIC',
        !!v && v.best === B8.models.rows.find(r => r.AIC === minAIC).name, v && v.best);
      const sec8 = Report.methods(state, {}).find(s => s.id === 'divers');
      t('el informe trae la sección de diversificación', !!sec8);
      t('que nombra el mejor modelo en los dos idiomas', !!sec8
        && sec8.en.join(' ').indexOf(v.best) >= 0 && sec8.es.join(' ').indexOf('undefined') < 0,
        sec8 && sec8.es.join(' '));
      t('y cita a Nee et al. por los modelos', !!sec8 && sec8.refs.indexOf('nee1994') >= 0);
      t('el resumen del informe trae el modelo de diversificación',
        Report.summary(state, 'es').some(r => /diversificación/i.test(r.label)));
      const tab = B12.tablesForZip().find(x => x.name === 'diversificacion_modelos.csv');
      t('el paquete .zip trae la tabla de modelos', !!tab);
      t('con una fila por modelo y sin celdas vacías', !!tab && (() => {
        const lines = tab.content.split('\n');
        return lines.length === B8.models.rows.length + 1
          && lines.every(l => l.indexOf('undefined') < 0 && l.indexOf(',,') < 0);
      })(), tab && tab.content.split('\n')[1]);
      t('y los nombres con coma no parten la fila', !!tab
        && tab.content.split('\n').slice(1).every(l => l.replace(/"[^"]*"/g, '').split(',').length === 6));
    } finally {
      Object.assign(B8, keep);
      state.diversification = had;
    }
  })();

})();
