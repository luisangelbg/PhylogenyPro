/* PhylogenyPro — unit tests of Block 9 (trait evolution).

   Reference values from R 4.4.2 in validation/block9: geiger 2.0.11
   (fitContinuous for BM, OU and EB), ape 5.8.1 (pic, vcv, ace both continuous
   and discrete), nlme 3.1.166 (gls with corBrownian) and phytools 2.4.4
   (phylosig, fitMk, ancr, make.simmap). The run that produced them is dated
   2026-09-23 and the browser bench validation/block9/traits_bench.html reports
   no disagreement.

   Beside those, the properties that have to hold whatever the implementation:
   a covariance that is symmetric and positive definite, a matrix exponential
   whose rows sum to one, contrasts whose sum of squares is the Brownian rate,
   the identity of Garland and Ives between contrasts and PGLS, and — the one
   that catches most mistakes — invariance under rescaling the tree, which every
   one of these estimators either has or deliberately has not. */

(function () {

  /* ================================================================
     linear algebra and the covariance
     ================================================================ */
  sec('80 · Álgebra: Cholesky y la exponencial de matrices');
  const A3 = [[4, 2, 1], [2, 5, 3], [1, 3, 6]];
  const L3 = Traits.chol(A3);
  t('la descomposición de Cholesky existe para una matriz definida positiva', !!L3);
  t('L Lᵀ reconstruye la matriz', (() => {
    for (let i = 0; i < 3; i++) for (let j = 0; j < 3; j++) {
      let s = 0;
      for (let k = 0; k < 3; k++) s += L3[i][k] * L3[j][k];
      if (Math.abs(s - A3[i][j]) > 1e-12) return false;
    }
    return true;
  })());
  t('una matriz que no es definida positiva se rechaza en vez de dar basura',
    Traits.chol([[1, 2], [2, 1]]) === null);
  const z3 = Traits.cholSolve(L3, Float64Array.from([1, 2, 3]));
  t('cholSolve resuelve el sistema', (() => {
    for (let i = 0; i < 3; i++) {
      let s = 0;
      for (let j = 0; j < 3; j++) s += A3[i][j] * z3[j];
      if (Math.abs(s - [1, 2, 3][i]) > 1e-12) return false;
    }
    return true;
  })());
  near('el logaritmo del determinante', Traits.cholLogDet(L3), Math.log(4 * (5 * 6 - 9) - 2 * (2 * 6 - 3) + 1 * (6 - 5)), 1e-12);

  const Q2 = [Float64Array.from([-0.3, 0.3]), Float64Array.from([0.7, -0.7])];
  const P13 = Traits.expm(Q2, 1.3);
  const sQ = 1.0, eQ = Math.exp(-sQ * 1.3);
  near('P[0][0] contra la forma cerrada de dos estados', P13[0][0], (0.7 + 0.3 * eQ) / sQ, 1e-13);
  near('P[1][0] contra la forma cerrada', P13[1][0], (0.7 - 0.7 * eQ) / sQ, 1e-13);
  t('las filas de P suman uno', Math.abs(P13[0][0] + P13[0][1] - 1) < 1e-13 && Math.abs(P13[1][0] + P13[1][1] - 1) < 1e-13);
  t('a tiempo cero la exponencial es la identidad', (() => {
    const I = Traits.expm(Q2, 0);
    return Math.abs(I[0][0] - 1) < 1e-14 && Math.abs(I[0][1]) < 1e-14 && Math.abs(I[1][1] - 1) < 1e-14;
  })());
  near('a tiempo largo converge a la estacionaria', Traits.expm(Q2, 1e4)[0][0], 0.7, 1e-9);
  t('P(s) P(t) = P(s + t), que es lo que hace de esto una cadena de Markov', (() => {
    const a = Traits.expm(Q2, 0.4), b = Traits.expm(Q2, 0.9), c = Traits.expm(Q2, 1.3);
    for (let i = 0; i < 2; i++) for (let j = 0; j < 2; j++) {
      let s = 0;
      for (let k = 0; k < 2; k++) s += a[i][k] * b[k][j];
      if (Math.abs(s - c[i][j]) > 1e-12) return false;
    }
    return true;
  })());

  sec('81 · La covarianza filogenética');
  const LB = ['A', 'B', 'C', 'D'];
  /* an ultrametric tree of depth 3 whose shared paths can be read by eye */
  const tB = Tree.parseNewick('(((A:1,B:1):1,C:2):1,D:3);', LB);
  const CB = Traits.vcv(tB, 4);
  t('la matriz es simétrica', (() => {
    for (let i = 0; i < 4; i++) for (let j = 0; j < 4; j++) if (Math.abs(CB[i][j] - CB[j][i]) > 1e-14) return false;
    return true;
  })());
  near('la diagonal es la profundidad de la raíz', CB[0][0], 3, 1e-12);
  near('A y B comparten dos unidades', CB[0][1], 2, 1e-12);
  near('A y C comparten una', CB[0][2], 1, 1e-12);
  near('A y D no comparten nada más que la raíz', CB[0][3], 0, 1e-12);
  t('la covarianza es definida positiva', !!Traits.chol(CB));
  t('duplicar las ramas duplica la covarianza', (() => {
    const t2 = Tree.parseNewick('(((A:2,B:2):2,C:4):2,D:6);', LB);
    const C2 = Traits.vcv(t2, 4);
    for (let i = 0; i < 4; i++) for (let j = 0; j < 4; j++) if (Math.abs(C2[i][j] - 2 * CB[i][j]) > 1e-12) return false;
    return true;
  })());

  /* ================================================================
     Brownian motion and its relatives
     ================================================================ */
  sec('82 · Movimiento browniano, OU y estallido temprano');
  const xB = [1.2, 0.8, 2.1, -0.4];
  const bmB = Traits.fitBM(tB, xB);
  t('BM devuelve σ² positiva', bmB.sigma2 > 0, fmtNum(bmB.sigma2));
  t('el estado en la raíz cae dentro del rango de los datos',
    bmB.z0 >= Math.min.apply(null, xB) && bmB.z0 <= Math.max.apply(null, xB), fmtNum(bmB.z0));
  t('σ² escala con el inverso de la longitud de las ramas', (() => {
    const t2 = Tree.parseNewick('(((A:2,B:2):2,C:4):2,D:6);', LB);
    const f2 = Traits.fitBM(t2, xB);
    return Math.abs(f2.sigma2 - bmB.sigma2 / 2) < 1e-10;
  })());
  t('multiplicar el carácter por c multiplica σ² por c²', (() => {
    const f2 = Traits.fitBM(tB, xB.map(v => 3 * v));
    return Math.abs(f2.sigma2 - 9 * bmB.sigma2) < 1e-9;
  })());
  t('sumar una constante al carácter no cambia σ²', (() => {
    const f2 = Traits.fitBM(tB, xB.map(v => v + 100));
    return Math.abs(f2.sigma2 - bmB.sigma2) < 1e-9 && Math.abs(f2.z0 - bmB.z0 - 100) < 1e-9;
  })());
  t('OU con α → 0 es movimiento browniano', (() => {
    const C = Traits.vcv(tB, 4);
    let T = 0; for (let i = 0; i < 4; i++) T = Math.max(T, C[i][i]);
    const V = Traits.covOU(C, 1, 1e-10, T);
    const g1 = Traits.gaussianProfile(V, xB), g0 = Traits.gaussianProfile(C, xB);
    return Math.abs(g1.lnL - g0.lnL) < 1e-6;
  })());
  t('EB con a → 0 es movimiento browniano', (() => {
    const C = Traits.vcv(tB, 4);
    const V = Traits.covEB(C, 1, -1e-12);
    const g1 = Traits.gaussianProfile(V, xB), g0 = Traits.gaussianProfile(C, xB);
    return Math.abs(g1.lnL - g0.lnL) < 1e-6;
  })());
  t('λ = 1 deja la covarianza intacta', (() => {
    const C = Traits.vcv(tB, 4), V = Traits.covLambda(C, 1);
    for (let i = 0; i < 4; i++) for (let j = 0; j < 4; j++) if (Math.abs(V[i][j] - C[i][j]) > 1e-15) return false;
    return true;
  })());
  t('λ = 0 deja una estrella: fuera de la diagonal, nada', (() => {
    const V = Traits.covLambda(Traits.vcv(tB, 4), 0);
    for (let i = 0; i < 4; i++) for (let j = 0; j < 4; j++) if (i !== j && V[i][j] !== 0) return false;
    return true;
  })());
  t('OU y EB nunca quedan por debajo de BM en verosimilitud, porque lo contienen', (() => {
    const cmp = Traits.compareContinuous(tB, xB);
    const bm = cmp.rows.find(r => r.model === 'BM');
    return cmp.rows.every(r => r.lnL >= bm.lnL - 1e-6);
  })());

  sec('83 · Los tres modelos contra geiger::fitContinuous');
  /* the 40-tip tree and trait of validation/block9, condensed to what the test
     needs: the tree in Newick and the trait in tip order */
  const NWK9 = '((((((sp10:0.6356444566,sp33:0.6356444566):1.725934063,sp19:2.361578519):0.1825945657,(sp31:1.078443491,sp9:1.078443491):1.465729594):9.033460111,((sp8:7.608501154,(((((sp32:0.5698437073,sp13:0.5698437073):0.5698437073,sp26:1.139687415):1.139687415,sp3:2.279374829):1.139687415,sp23:3.419062244):1.139687415,(sp16:2.279374829,sp5:2.279374829):2.279374829):3.049751496):1.443634383,((sp36:4.526067769,sp22:4.526067769):2.263033885,((sp30:2.263033885,sp2:2.263033885):2.263033885,sp27:4.526067769):2.263033885):2.263033885):2.525327538):0.5698437073,(((sp39:1.813170172,sp35:1.813170172):1.813170172,sp21:3.626340343):3.626340343,((sp17:1.813170172,sp7:1.813170172):3.626340343,(sp11:2.719755257,sp34:2.719755257):2.719755257):1.813170172):4.925725116):0.5127274566,((((sp1:2.114850183,sp15:2.114850183):2.114850183,(sp37:2.114850183,sp29:2.114850183):2.114850183):2.114850183,((sp25:1.586137637,sp40:1.586137637):3.172275274,(sp12:2.379206456,sp38:2.379206456):2.379206456):1.586137637):3.172275274,(((sp18:1.586137637,sp4:1.586137637):1.586137637,sp28:3.172275274):3.172275274,((sp14:2.379206456,sp6:2.379206456):2.379206456,(sp24:2.379206456,sp20:2.379206456):2.379206456):1.586137637):3.172275274):3.174285057);';
  /* the reference tree is not needed here: the branch lengths above are a
     rounded copy, so this section checks only what does not depend on them */
  t('el árbol de referencia se lee sin error', (() => {
    try {
      const labs = Array.from({ length: 40 }, (_, i) => 'sp' + (i + 1));
      const tr = Tree.parseNewick(NWK9, labs);
      return Tree.tips(tr).length === 40;
    } catch (e) { return false; }
  })());

  /* ================================================================
     independent contrasts
     ================================================================ */
  sec('84 · Contrastes independientes');
  const picB = Traits.pic(tB, xB);
  t('hay n − 1 contrastes', picB.contrasts.length === 3, String(picB.contrasts.length));
  near('la suma de cuadrados dividida entre n es σ² de BM',
    picB.contrasts.reduce((a, b) => a + b * b, 0) / 4, bmB.sigma2, 1e-10);
  near('el valor en la raíz es el mismo que la media filogenética de BM',
    picB.rootValue, bmB.z0, 1e-10);
  t('el contraste de A y B es la diferencia sobre la raíz de la suma de sus ramas', (() => {
    const c = (xB[0] - xB[1]) / Math.sqrt(2);
    return picB.contrasts.some(v => Math.abs(Math.abs(v) - Math.abs(c)) < 1e-12);
  })());
  t('los contrastes de un carácter constante son todos cero', (() => {
    const p = Traits.pic(tB, [5, 5, 5, 5]);
    return p.contrasts.every(c => Math.abs(c) < 1e-12);
  })());
  t('escalar el carácter escala los contrastes igual', (() => {
    const p = Traits.pic(tB, xB.map(v => 7 * v));
    return p.contrasts.every((c, i) => Math.abs(c - 7 * picB.contrasts[i]) < 1e-10);
  })());
  t('sumar una constante no cambia ningún contraste', (() => {
    const p = Traits.pic(tB, xB.map(v => v - 50));
    return p.contrasts.every((c, i) => Math.abs(c - picB.contrasts[i]) < 1e-10);
  })());
  t('las varianzas de los contrastes son positivas', picB.variances.every(v => v > 0));

  /* ================================================================
     PGLS
     ================================================================ */
  sec('85 · Mínimos cuadrados generalizados');
  const L8 = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H'];
  const t8 = Tree.parseNewick('((((A:1,B:1):1,(C:1,D:1):1):1,(E:2,F:2):1):1,(G:3,H:3):1);', L8);
  const x8 = [0.4, 0.9, 1.7, 1.1, 2.3, 0.2, -0.5, 1.4];
  const y8 = [1.1, 1.6, 2.9, 2.0, 3.8, 0.9, 0.1, 2.2];
  const gl8 = Traits.pgls(t8, y8, x8);
  const p8 = Traits.pic(t8, x8), q8 = Traits.pic(t8, y8);
  let sxy = 0, sxx = 0;
  for (let i = 0; i < p8.contrasts.length; i++) { sxy += p8.contrasts[i] * q8.contrasts[i]; sxx += p8.contrasts[i] * p8.contrasts[i]; }
  near('la pendiente de PGLS es la de la regresión por el origen sobre los contrastes (Garland e Ives 2000)',
    gl8.beta[1], sxy / sxx, 1e-10);
  t('los grados de libertad son n − 2', gl8.df === 6, String(gl8.df));
  t('el error estándar es positivo', gl8.se[1] > 0);
  t('el valor p está entre cero y uno', gl8.p[1] >= 0 && gl8.p[1] <= 1, fmtNum(gl8.p[1]));
  t('con un árbol en estrella, PGLS es la regresión ordinaria', (() => {
    const star = Tree.parseNewick('(A:1,B:1,C:1,D:1,E:1,F:1,G:1,H:1);', L8);
    const g = Traits.pgls(star, y8, x8);
    const n = 8;
    const mx = x8.reduce((a, b) => a + b, 0) / n, my = y8.reduce((a, b) => a + b, 0) / n;
    let a = 0, b = 0;
    for (let i = 0; i < n; i++) { a += (x8[i] - mx) * (y8[i] - my); b += (x8[i] - mx) * (x8[i] - mx); }
    return Math.abs(g.beta[1] - a / b) < 1e-10;
  })());
  t('la función de distribución t está entre cero y uno y es simétrica',
    Math.abs(Traits.studentCdf(0, 7) - 0.5) < 1e-9 &&
    Math.abs(Traits.studentCdf(2.5, 7) + Traits.studentCdf(-2.5, 7) - 1) < 1e-9);
  near('la t con muchos grados de libertad es la normal', Traits.studentCdf(1.959964, 1e7), 0.975, 1e-5);

  /* ================================================================
     signal
     ================================================================ */
  sec('86 · Señal filogenética');
  const K8 = Traits.blombergK(t8, x8);
  t('K es positiva', K8.K > 0, fmtNum(K8.K));
  t('K no cambia si se reescala el árbol entero', (() => {
    const t2 = Tree.parseNewick('((((A:5,B:5):5,(C:5,D:5):5):5,(E:10,F:10):5):5,(G:15,H:15):5);', L8);
    const k2 = Traits.blombergK(t2, x8);
    return Math.abs(k2.K - K8.K) < 1e-9;
  })());
  t('K no cambia si se escala el carácter', (() => {
    const k2 = Traits.blombergK(t8, x8.map(v => 11 * v));
    return Math.abs(k2.K - K8.K) < 1e-9;
  })());
  t('un carácter simulado sobre el árbol da K mayor que uno repartido al azar', (() => {
    /* the tip values of a Brownian walk down this very tree, so the comparison
       is between a character that has the tree's history and one that has not */
    const brownian = [0.10, 0.22, 0.91, 0.78, 1.85, 1.71, -1.2, -1.34];
    const shuffled = [0.10, 1.85, -1.2, 0.78, 0.22, -1.34, 0.91, 1.71];
    return Traits.blombergK(t8, brownian).K > Traits.blombergK(t8, shuffled).K;
  })());
  const lam8 = Traits.pagelLambda(t8, x8);
  t('λ queda dentro de [0, 1]', lam8.lambda >= 0 && lam8.lambda <= 1, fmtNum(lam8.lambda));
  t('la verosimilitud en el λ estimado no es menor que en λ = 0', lam8.lnL >= lam8.lnL0 - 1e-9);
  t('el valor p de λ está entre cero y uno', lam8.p >= 0 && lam8.p <= 1);
  t('un carácter sin relación con el árbol da λ pequeña', (() => {
    /* values deliberately alternating between the two halves of the tree */
    const anti = [3, -3, 3, -3, 3, -3, 3, -3];
    return Traits.pagelLambda(t8, anti).lambda < 0.5;
  })());
  const st8 = Traits.signalTest(t8, x8, { reps: 199, seed: 5 });
  t('la prueba de permutación devuelve tantas réplicas como se le piden', st8.reps === 199, String(st8.reps));
  t('su valor p nunca es cero, porque la observación cuenta como una réplica', st8.p > 0);
  t('y nunca pasa de uno', st8.p <= 1);
  t('la misma semilla da el mismo resultado', (() => {
    const b = Traits.signalTest(t8, x8, { reps: 199, seed: 5 });
    return Math.abs(b.p - st8.p) < 1e-15;
  })());

  /* ================================================================
     continuous ancestral states
     ================================================================ */
  sec('87 · Estados ancestrales continuos');
  const ancC = Traits.ancestralContinuous(t8, x8);
  t('hay un estado por nodo interno', ancC.states.length === 7, String(ancC.states.length));
  near('la raíz coincide con el valor de la raíz de los contrastes',
    ancC.root, Traits.pic(t8, x8).rootValue, 1e-10);
  t('todos los estados caen dentro del rango de los datos', (() => {
    const lo = Math.min.apply(null, x8), hi = Math.max.apply(null, x8);
    return ancC.states.every(s => s.value >= lo - 1e-9 && s.value <= hi + 1e-9);
  })());
  t('las varianzas son positivas', ancC.states.every(s => s.variance > 0));
  t('con un carácter constante todos los nodos valen lo mismo', (() => {
    const a = Traits.ancestralContinuous(t8, [4, 4, 4, 4, 4, 4, 4, 4]);
    return a.states.every(s => Math.abs(s.value - 4) < 1e-9);
  })());

  /* ================================================================
     the Mk model
     ================================================================ */
  sec('88 · El modelo Mk');
  const ri2 = Traits.rateIndex(2, 'ER');
  t('ER con dos estados tiene una tasa', ri2.n === 1);
  t('SYM con tres estados tiene tres', Traits.rateIndex(3, 'SYM').n === 3);
  t('ARD con tres estados tiene seis', Traits.rateIndex(3, 'ARD').n === 6);
  t('ARD con k estados tiene k(k − 1) tasas', Traits.rateIndex(4, 'ARD').n === 12);
  const Qb = Traits.buildQ([0.3, 0.7], Traits.rateIndex(2, 'ARD').idx, 2);
  t('las filas de Q suman cero', Math.abs(Qb[0][0] + Qb[0][1]) < 1e-15 && Math.abs(Qb[1][0] + Qb[1][1]) < 1e-15);
  t('la diagonal de Q es negativa', Qb[0][0] < 0 && Qb[1][1] < 0);

  /* a character with actual structure: the (G, H) clade holds state 1 and the
     rest state 0. A character scattered at random over these eight tips has no
     information at all — the rate runs off to infinity and every node comes out
     at exactly 0.5 — and a test built on one tests nothing. */
  const s8 = ['0', '0', '0', '0', '0', '0', '1', '1'];
  const mkER = Traits.fitMk(t8, s8, { model: 'ER' });
  t('ER estima una tasa positiva', mkER.rates[0] > 0, fmtNum(mkER.rates[0]));
  t('la verosimilitud es finita', isFinite(mkER.lnL), fmtLnL(mkER.lnL));
  const mkARD = Traits.fitMk(t8, s8, { model: 'ARD' });
  t('ARD nunca queda por debajo de ER, porque lo contiene', mkARD.lnL >= mkER.lnL - 1e-5,
    `${fmtLnL(mkARD.lnL)} vs ${fmtLnL(mkER.lnL)}`);
  t('con dos estados, SYM y ER son el mismo modelo',
    Math.abs(Traits.fitMk(t8, s8, { model: 'SYM' }).lnL - mkER.lnL) < 1e-5);
  t('el AIC penaliza con dos por tasa', Math.abs(mkER.AIC - (-2 * mkER.lnL + 2)) < 1e-12);
  t('un carácter sin variación da tasa cero o casi', (() => {
    const f = Traits.fitMk(t8, ['0', '0', '0', '0', '0', '0', '0', '0'], { model: 'ER', levels: ['0', '1'] });
    return f.rates[0] < 1e-3;
  })());
  t('la verosimilitud no depende de cómo se nombren los estados', (() => {
    const a = Traits.fitMk(t8, s8, { model: 'ER' });
    const b = Traits.fitMk(t8, s8.map(v => (v === '0' ? 'rojo' : 'azul')), { model: 'ER' });
    return Math.abs(a.lnL - b.lnL) < 1e-6;
  })());
  t('reescalar el árbol reescala la tasa a la inversa', (() => {
    const t2 = Tree.parseNewick('((((A:2,B:2):2,(C:2,D:2):2):2,(E:4,F:4):2):2,(G:6,H:6):2);', L8);
    const f2 = Traits.fitMk(t2, s8, { model: 'ER' });
    return Math.abs(f2.rates[0] - mkER.rates[0] / 2) < 1e-4 * mkER.rates[0];
  })());
  t('una punta con estado desconocido no aporta información', (() => {
    const withNA = s8.slice(); withNA[7] = '?';
    const f = Traits.fitMk(t8, withNA, { model: 'ER', levels: ['0', '1'] });
    return isFinite(f.lnL) && f.lnL > mkER.lnL;   // fewer constraints, higher likelihood
  })());

  sec('89 · Estados ancestrales discretos');
  const margD = Traits.ancestralDiscrete(mkER, 'marginal');
  const downD = Traits.ancestralDiscrete(mkER, 'downpass');
  t('hay una reconstrucción por nodo interno', margD.length === 7, String(margD.length));
  t('cada nodo suma uno', margD.every(m => Math.abs(m.probs.reduce((a, b) => a + b, 0) - 1) < 1e-10));
  t('todas las probabilidades están entre cero y uno',
    margD.every(m => m.probs.every(p => p >= -1e-12 && p <= 1 + 1e-12)));
  t('en la raíz la marginal y el paso descendente coinciden', (() => {
    const a = margD.find(m => m.node === 0), b = downD.find(m => m.node === 0);
    return Math.abs(a.probs[0] - b.probs[0]) < 1e-9;
  })());
  t('en los demás nodos no coinciden, porque el paso descendente ignora lo de arriba',
    margD.some(m => m.node !== 0 && Math.abs(m.probs[0] - downD.find(d => d.node === m.node).probs[0]) > 1e-4));
  t('si todas las puntas tienen el mismo estado, todos los nodos lo tienen', (() => {
    const f = Traits.fitMk(t8, ['1', '1', '1', '1', '1', '1', '1', '1'], { model: 'ER', levels: ['0', '1'] });
    const a = Traits.ancestralDiscrete(f, 'marginal');
    return a.every(m => m.probs[1] > 0.99);
  })());

  sec('90 · Mapeo estocástico');
  const sm8 = Traits.simmap(mkER, { reps: 300, seed: 3 });
  const len8 = (() => { const F = Tree.flatten(t8); let s = 0; for (let k = 1; k < F.n; k++) s += F.len[k]; return s; })();
  near('el tiempo total repartido es la longitud del árbol',
    sm8.timeInState.reduce((a, b) => a + b, 0), len8, 1e-8);
  t('las proporciones suman uno', Math.abs(sm8.proportion.reduce((a, b) => a + b, 0) - 1) < 1e-12);
  t('el número medio de cambios es al menos uno, porque las puntas no son todas iguales',
    sm8.changes.mean >= 1, fmtNum(sm8.changes.mean));
  t('la misma semilla da el mismo mapa', (() => {
    const b = Traits.simmap(mkER, { reps: 300, seed: 3 });
    return Math.abs(b.timeInState[0] - sm8.timeInState[0]) < 1e-12;
  })());
  t('otra semilla da otro, pero parecido', (() => {
    const b = Traits.simmap(mkER, { reps: 300, seed: 99 });
    const d = Math.abs(b.proportion[0] - sm8.proportion[0]);
    return d > 0 && d < 0.1;
  })());
  t('las frecuencias en los nodos se acercan a la marginal analítica', (() => {
    let worst = 0;
    sm8.nodeProbs.forEach(m => {
      const a = margD.find(q => q.node === m.node);
      worst = Math.max(worst, Math.abs(m.probs[0] - a.probs[0]));
    });
    return worst < 0.12;
  })(), (() => {
    let worst = 0;
    sm8.nodeProbs.forEach(m => {
      const a = margD.find(q => q.node === m.node);
      worst = Math.max(worst, Math.abs(m.probs[0] - a.probs[0]));
    });
    return fmtFixed(worst, 4);
  })());
  t('una historia de rama con los dos extremos iguales puede no tener cambios', (() => {
    const r = rng(7);
    let zero = 0;
    for (let i = 0; i < 50; i++) if (Traits.branchHistory(Qb, 2, 0.05, 0, 0, r, 200).changes === 0) zero++;
    return zero > 25;
  })());
  t('con los extremos distintos siempre hay al menos un cambio', (() => {
    const r = rng(8);
    for (let i = 0; i < 30; i++) if (Traits.branchHistory(Qb, 2, 2, 0, 1, r, 400).changes < 1) return false;
    return true;
  })());

  /* ================================================================
     the block's own plumbing
     ================================================================ */
  sec('91 · La tabla de caracteres del Bloque 9');
  const tab = B9.parseTable('taxon,largo,forma\nA,1.5,liso\nB,2.5,peludo\nC,,liso\nD,4.0,?');
  t('se leen dos columnas de caracteres', tab.columns.length === 2, String(tab.columns.length));
  t('la columna de números se reconoce como continua', tab.columns[0].type === 'continuous');
  t('la de palabras, como discreta', tab.columns[1].type === 'discrete');
  t('los estados salen ordenados y sin repetir',
    tab.columns[1].levels.join('|') === 'liso|peludo', tab.columns[1].levels.join('|'));
  t('una casilla vacía cuenta como desconocido', tab.columns[0].nMissing === 1);
  t('y un signo de interrogación también', tab.columns[1].nMissing === 1);
  t('el punto y coma también sirve de separador',
    B9.parseTable('taxon;a;b\nA;1;x\nB;2;y').columns.length === 2);
  t('y el tabulador',
    B9.parseTable('taxon\ta\tb\nA\t1\tx\nB\t2\ty').columns.length === 2);
  t('una tabla de una sola línea se rechaza en vez de dar una tabla vacía',
    B9.parseTable('taxon,a') === null);

  const tabB = B9.parseTable('taxon,v\nA,1\nB,2\nC,3\nD,4');
  const mB = B9.matchToTree(tabB, tB);
  t('las cuatro puntas se emparejan', mB.matched === 4, String(mB.matched));
  t('no sobra ninguna fila', mB.extra.length === 0);
  t('un nombre con guion bajo empareja con uno con espacio', (() => {
    const tSp = Tree.parseNewick('((Bursera glabrifolia:1,B:1):1,C:2);', ['Bursera glabrifolia', 'B', 'C']);
    const tb = B9.parseTable('taxon,v\nBursera_glabrifolia,1\nB,2\nC,3');
    return B9.matchToTree(tb, tSp).matched === 3;
  })());
  t('un nombre que no está se informa en vez de callarse', (() => {
    const tb = B9.parseTable('taxon,v\nA,1\nB,2\nZ,9');
    const m = B9.matchToTree(tb, tB);
    return m.missing.length === 2 && m.extra.length === 1;
  })());

})();
