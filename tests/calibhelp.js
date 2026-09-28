/* PhylogenyPro — unit tests of the calibration assistant.

   The assistant's job is to say, in years, what a prior the user typed actually
   asserts, and to catch the handful of mistakes that ruin a dating study. Both
   halves are testable: the first against the closed-form quantiles of each
   distribution, the second by feeding it the mistake and checking it complains.

   The arithmetic matters more than it looks. An offset lognormal with M = 1.6
   and S = 0.7 offset at 20 has its median at 20 + e^1.6 = 24.95, not at 21.6
   and not at 1.6 — and a user who believes either of those wrong numbers will
   publish a date that is off by millions of years. */

(function () {

  const near = (a, b, tol) => Math.abs(a - b) <= (tol == null ? 1e-6 : tol);
  const Q = CalibHelp.quantiles;

  sec('1 · what the prior actually says');

  {
    const q = Q({ type: 'fixed', value: 20 });
    t('una edad fija no tiene incertidumbre', q.exact === true && q.lo === 20 && q.hi === 20);
  }
  {
    const q = Q({ type: 'uniform', min: 3.1, max: 3.5 });
    t('el uniforme tiene su mediana en el centro', near(q.med, 3.3));
    t('y su 95 % dentro del intervalo', q.lo > 3.1 && q.hi < 3.5);
    t('y recuerda sus dos topes', q.floor === 3.1 && q.ceiling === 3.5);
  }
  {
    const q = Q({ type: 'normal', mean: 30, sd: 5 });
    t('la normal es simétrica', near(q.med, 30) && near(q.hi - q.med, q.med - q.lo, 1e-9));
    t('y su 95 % son ±1.96 desviaciones', near(q.hi, 30 + 1.959963985 * 5, 1e-6));
  }
  {
    /* mediana = desplazamiento + e^M; los extremos, con ±1.96·S en el exponente */
    const q = Q({ type: 'lognormal', min: 20, M: 1.6, S: 0.7 });
    t('la lognormal desplazada pone su mediana en desplazamiento + e^M', near(q.med, 20 + Math.exp(1.6), 1e-9));
    t('su extremo inferior en desplazamiento + e^(M−1.96S)', near(q.lo, 20 + Math.exp(1.6 - 1.959963985 * 0.7), 1e-9));
    t('su extremo superior en desplazamiento + e^(M+1.96S)', near(q.hi, 20 + Math.exp(1.6 + 1.959963985 * 0.7), 1e-9));
    t('y nunca baja del desplazamiento', q.floor === 20 && q.lo > 20);
    t('la cola superior es más larga que la inferior', (q.hi - q.med) > (q.med - q.lo));
  }
  {
    /* mediana = desplazamiento + media·ln2 */
    const q = Q({ type: 'exponential', min: 20, mean: 8 });
    t('la exponencial desplazada pone su mediana en desplazamiento + media·ln2', near(q.med, 20 + 8 * Math.LN2, 1e-9));
    t('su 97.5 % en desplazamiento − media·ln(0.025)', near(q.hi, 20 - 8 * Math.log(0.025), 1e-9));
    t('y tampoco baja del desplazamiento', q.lo > 20);
  }
  {
    t('un prior imposible no devuelve nada',
      Q({ type: 'uniform', min: 10, max: 10 }) === null &&
      Q({ type: 'normal', mean: 10, sd: 0 }) === null &&
      Q({ type: 'lognormal', min: 5, M: 1, S: -1 }) === null);
    t('y una forma desconocida tampoco', Q({ type: 'inventada' }) === null);
  }
  {
    const d = CalibHelp.describe({ type: 'lognormal', min: 20, M: 1.6, S: 0.7 }, 'Ma');
    t('la frase trae los dos idiomas', Array.isArray(d) && d.length === 2 && d[0] !== d[1]);
    t('y dice la mediana en años', d[0].indexOf('24.9') > 0 || d[0].indexOf('25.0') > 0);
    t('y el piso que impone', d[0].indexOf('no puede ser más joven') > 0);
  }

  sec('2 · the mistakes it has to catch');

  const nivel = (c, ctx) => CalibHelp.advise(c, ctx).map(a => a.level);
  const dice = (c, ctx, re) => CalibHelp.advise(c, ctx).some(a => re.test(a.es));

  t('un fósil como edad fija es el error grave',
    nivel({ type: 'fixed', value: 20 }, { source: 'fossil', others: [{}] })[0] === 'bad');
  t('y explica que el fósil da la edad mínima',
    dice({ type: 'fixed', value: 20 }, { source: 'fossil', others: [{}] }, /edad mínima/));
  t('un fósil con prior normal también se avisa',
    dice({ type: 'normal', mean: 20, sd: 3 }, { source: 'fossil', others: [{}] }, /simétrica|más joven/));
  t('una lognormal sin desplazamiento no impone mínimo',
    dice({ type: 'lognormal', min: 0, M: 2, S: 0.6 }, { source: 'fossil', others: [{}] }, /desplazamiento está en cero/));
  t('un intervalo vacío se detecta',
    nivel({ type: 'uniform', min: 10, max: 5 }, { others: [{}] }).indexOf('bad') >= 0);
  t('una normal que llega al futuro se avisa',
    dice({ type: 'normal', mean: 30, sd: 20 }, { source: 'secondary', others: [{}] }, /futuro/));
  /* la anchura se mide por la razón entre los extremos: con la diferencia
     relativa a la mediana, un uniforme no podía pasar de 1.9 y el aviso nunca
     se disparaba en la forma que más se usa para ir a lo vago */
  t('un uniforme enorme se declara poco informativo',
    dice({ type: 'uniform', min: 1, max: 200 }, { others: [{}] }, /apenas informa/));
  t('una lognormal muy vaga también',
    dice({ type: 'lognormal', min: 20, M: 1.6, S: 2 }, { source: 'fossil', others: [{}] }, /apenas informa/));
  t('y un prior razonable no se marca',
    !dice({ type: 'uniform', min: 3.1, max: 3.5 }, { others: [{}] }, /apenas informa/) &&
    !dice({ type: 'lognormal', min: 20, M: 1.6, S: 0.7 }, { source: 'fossil', others: [{}] }, /apenas informa/));
  t('y uno demasiado estrecho, exceso de confianza',
    dice({ type: 'normal', mean: 30, sd: 0.2 }, { source: 'secondary', others: [{}] }, /exceso de confianza/));
  t('avisa de corona contra tallo cuando hay un clado y un fósil',
    dice({ type: 'lognormal', min: 20, M: 1.5, S: 0.7 }, { source: 'fossil', others: [{}], cladeSize: 4 }, /corona/));
  t('y no lo hace cuando se calibra la raíz',
    !dice({ type: 'lognormal', min: 20, M: 1.5, S: 0.7 }, { source: 'fossil', others: [{}], cladeSize: 0 }, /corona/));
  t('avisa cuando va a ser la única calibración',
    dice({ type: 'uniform', min: 3, max: 4 }, { others: [] }, /única/));
  t('y no cuando ya hay otras',
    !dice({ type: 'uniform', min: 3, max: 4 }, { others: [{}, {}] }, /única/));
  t('todos los avisos vienen en los dos idiomas',
    CalibHelp.advise({ type: 'fixed', value: 20 }, { source: 'fossil', others: [] })
      .every(a => a.es && a.en && a.es !== a.en));

  sec('3 · what it proposes');

  {
    const s = CalibHelp.suggest('fossil', 20);
    t('para un fósil propone una lognormal desplazada', s.type === 'lognormal');
    t('desplazada justo en la edad del fósil', s.min === 20);
    const q = Q(s);
    t('y con la mediana por encima del fósil, no en él', q.med > 20 && q.med < 40);
    t('que nunca permite un nodo más joven que el fósil', q.lo > 20);
  }
  {
    const s = CalibHelp.suggest('geo', 3.1, 3.5);
    t('para un acontecimiento geológico propone un intervalo', s.type === 'uniform' && s.min === 3.1 && s.max === 3.5);
  }
  {
    const s = CalibHelp.suggest('secondary', 50);
    t('para una edad heredada propone una normal', s.type === 'normal' && s.mean === 50);
    t('con una desviación que no finge certeza', s.sd > 0 && s.sd / s.mean > 0.05);
  }
  {
    t('para una edad conocida, fija', CalibHelp.suggest('known', 66).type === 'fixed');
    t('y una edad inválida no propone nada', CalibHelp.suggest('fossil', -3) === null && CalibHelp.suggest('fossil', NaN) === null);
  }
  {
    const faltan = CalibHelp.SOURCES.filter(s => !s.name || s.name.length !== 2 || !s.why || s.why.length !== 2 || !s.why[1]);
    t('las cuatro procedencias están explicadas en dos idiomas', CalibHelp.SOURCES.length === 4 && faltan.length === 0);
  }

})();
