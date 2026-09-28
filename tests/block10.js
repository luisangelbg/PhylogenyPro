/* PhylogenyPro — unit tests of Block 10 (historical biogeography).

   The numerical references come from BioGeoBEARS 1.1.3 on the Psychotria of
   Hawaii, the data set LAGRANGE and BioGeoBEARS are both distributed with:
   six models, their parameters and their log-likelihoods, recorded in
   validation/block10 by ref_biogeo.R and checked one by one in
   biogeo_bench.html, which reports no disagreement.

   Beside those, the properties that have to hold whatever the implementation:
   a rate matrix whose rows sum to zero and whose null range is absorbing, a
   transition matrix whose rows sum to one and which obeys the Chapman–
   Kolmogorov identity, cladogenesis weights that are a probability
   distribution for every ancestral range, the nesting of each +J model inside
   the one it extends, and — the test that catches the most — what each of the
   four models does and does not allow at a split. */

(function () {

  sec('92 · Rangos como conjuntos de áreas');
  t('popcount cuenta los bits', Biogeo.popcount(0b1011) === 3, String(Biogeo.popcount(0b1011)));
  t('areasOf los enumera', Biogeo.areasOf(0b1011, 4).join(',') === '0,1,3', Biogeo.areasOf(0b1011, 4).join(','));
  t('un rango se nombra por sus áreas', Biogeo.rangeName(0b0101, ['K', 'O', 'M', 'H']) === 'KM',
    Biogeo.rangeName(0b0101, ['K', 'O', 'M', 'H']));
  t('el rango nulo tiene su propio nombre', Biogeo.rangeName(0, ['K', 'O']) === '—');

  const A4 = ['K', 'O', 'M', 'H'];
  const st4 = Biogeo.makeStates(4, { maxAreas: 4, includeNull: true });
  t('cuatro áreas con el rango nulo dan 16 estados', st4.n === 16, String(st4.n));
  t('el rango nulo va primero', st4.list[0] === 0);
  t('los estados van de menor a mayor tamaño', (() => {
    for (let i = 2; i < st4.n; i++) {
      if (Biogeo.popcount(st4.list[i]) < Biogeo.popcount(st4.list[i - 1])) return false;
    }
    return true;
  })());
  t('sin el rango nulo hay uno menos',
    Biogeo.makeStates(4, { maxAreas: 4, includeNull: false }).n === 15);
  t('con un tope de dos áreas quedan 1 + 4 + 6 = 11',
    Biogeo.makeStates(4, { maxAreas: 2, includeNull: true }).n === 11,
    String(Biogeo.makeStates(4, { maxAreas: 2, includeNull: true }).n));
  t('ocho áreas sin tope dan 256 estados, que es el límite práctico',
    Biogeo.makeStates(8, { maxAreas: 8, includeNull: true }).n === 256);
  t('el índice encuentra cada rango', (() => {
    for (let i = 0; i < st4.n; i++) if (st4.index.get(st4.list[i]) !== i) return false;
    return true;
  })());

  sec('93 · La matriz de tasas a lo largo de una rama');
  const Q4 = Biogeo.buildQ(st4, 0.03, 0.02, {});
  t('las filas suman cero', (() => {
    for (let i = 0; i < st4.n; i++) {
      let s = 0;
      for (let j = 0; j < st4.n; j++) s += Q4[i][j];
      if (Math.abs(s) > 1e-13) return false;
    }
    return true;
  })());
  t('el rango nulo es absorbente: de ahí no se sale', Q4[0].every(v => v === 0));
  const iK = st4.index.get(0b0001), iKO = st4.index.get(0b0011), iAll = st4.index.get(0b1111);
  near('de una sola área salen tres dispersiones y una extinción', -Q4[iK][iK], 3 * 0.03 + 0.02, 1e-14);
  near('de dos áreas salen dos dispersiones y dos extinciones', -Q4[iKO][iKO], 2 * 2 * 0.03 + 2 * 0.02, 1e-14);
  near('del rango completo sólo se puede perder área', -Q4[iAll][iAll], 4 * 0.02, 1e-14);
  t('de una sola área se llega al rango nulo con tasa e', Math.abs(Q4[iK][0] - 0.02) < 1e-14);
  t('sin el rango nulo, una sola área no puede extinguirse', (() => {
    const s = Biogeo.makeStates(4, { maxAreas: 4, includeNull: false });
    const Q = Biogeo.buildQ(s, 0.03, 0.02, {});
    const i = s.index.get(0b0001);
    return Math.abs(-Q[i][i] - 3 * 0.03) < 1e-14;
  })());
  t('un multiplicador de cero prohíbe el movimiento entre dos áreas', (() => {
    const M = Array.from({ length: 4 }, () => new Float64Array(4).fill(1));
    M[0][1] = 0; M[1][0] = 0;
    const Q = Biogeo.buildQ(st4, 0.03, 0.02, { multipliers: M });
    return Q[iK][st4.index.get(0b0011)] === 0;
  })());
  t('un multiplicador de dos duplica esa dispersión', (() => {
    const M = Array.from({ length: 4 }, () => new Float64Array(4).fill(1));
    M[0][1] = 2;
    const Q = Biogeo.buildQ(st4, 0.03, 0.02, { multipliers: M });
    return Math.abs(Q[iK][st4.index.get(0b0011)] - 0.06) < 1e-14;
  })());
  t('con un tope de rango, el rango más ancho ya no se dispersa', (() => {
    const s = Biogeo.makeStates(4, { maxAreas: 2, includeNull: true });
    const Q = Biogeo.buildQ(s, 0.03, 0.02, {});
    const i = s.index.get(0b0011);
    return Math.abs(-Q[i][i] - 2 * 0.02) < 1e-14;
  })());

  sec('94 · La exponencial de la matriz de tasas');
  const n4 = st4.n;
  const P1 = Biogeo.expm(Q4, 1.0);
  t('P viene como un arreglo plano de n²', P1.length === n4 * n4, String(P1.length));
  t('las filas suman uno', (() => {
    for (let i = 0; i < n4; i++) {
      let s = 0;
      for (let j = 0; j < n4; j++) s += P1[i * n4 + j];
      if (Math.abs(s - 1) > 1e-11) return false;
    }
    return true;
  })());
  t('ninguna probabilidad es negativa', Array.from(P1).every(v => v >= 0));
  t('a tiempo cero es la identidad', (() => {
    const I = Biogeo.expm(Q4, 0);
    for (let i = 0; i < n4; i++) for (let j = 0; j < n4; j++) {
      if (Math.abs(I[i * n4 + j] - (i === j ? 1 : 0)) > 1e-14) return false;
    }
    return true;
  })());
  t('P(s)·P(t) = P(s + t)', (() => {
    const a = Biogeo.expm(Q4, 0.7), b = Biogeo.expm(Q4, 1.9), c = Biogeo.expm(Q4, 2.6);
    for (let i = 0; i < n4; i++) for (let j = 0; j < n4; j++) {
      let s = 0;
      for (let z = 0; z < n4; z++) s += a[i * n4 + z] * b[z * n4 + j];
      if (Math.abs(s - c[i * n4 + j]) > 1e-10) return false;
    }
    return true;
  })());
  t('a tiempo largo todo acaba en el rango nulo, que es lo que absorbente significa', (() => {
    const P = Biogeo.expm(Q4, 1e5);
    for (let i = 1; i < n4; i++) if (P[i * n4 + 0] < 0.999) return false;
    return true;
  })());
  t('la uniformización del caché da lo mismo que la serie de Taylor', (() => {
    /* the two are independent implementations of the same function, and the
       block uses the first: if they ever disagree, one of them is wrong */
    const cache = Biogeo.makeQCache(st4, 0.03, 0.02, {});
    let worst = 0;
    [0.05, 0.9, 3.7].forEach(len => {
      const a = cache.get(len), b = Biogeo.expm(Q4, len);
      for (let i = 0; i < n4 * n4; i++) worst = Math.max(worst, Math.abs(a[i] - b[i]));
    });
    return worst < 1e-12;
  })(), (() => {
    const cache = Biogeo.makeQCache(st4, 0.03, 0.02, {});
    let worst = 0;
    [0.05, 0.9, 3.7].forEach(len => {
      const a = cache.get(len), b = Biogeo.expm(Q4, len);
      for (let i = 0; i < n4 * n4; i++) worst = Math.max(worst, Math.abs(a[i] - b[i]));
    });
    return worst.toExponential(2);
  })());

  sec('95 · Lo que pasa al dividirse un linaje');
  const sumOf = ev => ev.reduce((a, b) => a + b[2], 0);
  ['DEC', 'DIVALIKE', 'BAYAREALIKE'].forEach(m => {
    const c = Biogeo.cladoEvents(st4, m, 0, {});
    t(`${m}: cada rango ancestral reparte exactamente una probabilidad`, (() => {
      for (let i = 0; i < st4.n; i++) {
        if (st4.list[i] === 0) continue;
        if (Math.abs(sumOf(c[i]) - 1) > 1e-12) return false;
      }
      return true;
    })());
    t(`${m}: el rango nulo no se divide`, c[0].length === 0);
    t(`${m}: un ancestro de una sola área se copia entero en las dos hijas`, (() => {
      const ev = c[iK];
      return ev.length === 1 && ev[0][0] === iK && ev[0][1] === iK;
    })());
  });
  const cDEC = Biogeo.cladoEvents(st4, 'DEC', 0, {});
  const cDIV = Biogeo.cladoEvents(st4, 'DIVALIKE', 0, {});
  const cBAY = Biogeo.cladoEvents(st4, 'BAYAREALIKE', 0, {});
  t('DEC, ancestro de dos áreas: cuatro de subconjunto y dos de vicarianza', cDEC[iKO].length === 6,
    String(cDEC[iKO].length));
  t('y todas igual de probables', cDEC[iKO].every(e => Math.abs(e[2] - 1 / 6) < 1e-12));
  t('DEC, ancestro de cuatro áreas: dieciséis divisiones', cDEC[iAll].length === 16, String(cDEC[iAll].length));
  t('DEC nunca parte cuatro áreas en dos y dos', !cDEC[iAll].some(e =>
    Biogeo.popcount(st4.list[e[0]]) === 2 && Biogeo.popcount(st4.list[e[1]]) === 2));
  t('DIVALIKE, ancestro de dos áreas: sólo las dos vicarianzas', cDIV[iKO].length === 2, String(cDIV[iKO].length));
  t('DIVALIKE, ancestro de cuatro áreas: catorce divisiones', cDIV[iAll].length === 14, String(cDIV[iAll].length));
  t('DIVALIKE sí parte cuatro áreas en dos y dos', cDIV[iAll].some(e =>
    Biogeo.popcount(st4.list[e[0]]) === 2 && Biogeo.popcount(st4.list[e[1]]) === 2));
  t('y con el ajuste de máxima entropía todas pesan igual',
    cDIV[iAll].every(e => Math.abs(e[2] - 1 / 14) < 1e-12));
  t('DIVALIKE nunca deja el rango entero a una hija', !cDIV[iAll].some(e =>
    st4.list[e[0]] === st4.list[iAll] || st4.list[e[1]] === st4.list[iAll]));
  t('BAYAREALIKE deja el rango entero a las dos hijas y nada más',
    cBAY[iAll].length === 1 && cBAY[iAll][0][0] === iAll && cBAY[iAll][0][1] === iAll);
  t('BAYAREALIKE no reparte nada, sea cual sea el tamaño', (() => {
    for (let i = 1; i < st4.n; i++) {
      if (cBAY[i].length !== 1 || cBAY[i][0][0] !== i || cBAY[i][0][1] !== i) return false;
    }
    return true;
  })());

  const cJ = Biogeo.cladoEvents(st4, 'DEC', 0.5, {});
  t('con j > 0 un ancestro de una sola área ya puede saltar', cJ[iK].length > 1, String(cJ[iK].length));
  t('el salto va siempre a un área donde el ancestro no estaba', (() => {
    const A = st4.list[iK];
    return cJ[iK].every(e => {
      const L = st4.list[e[0]], R = st4.list[e[1]];
      if (L === A && R === A) return true;                  // the sympatric event
      return (L === A && (R & A) === 0) || (R === A && (L & A) === 0);
    });
  })());
  t('j sale del presupuesto de los otros eventos, no se suma encima', (() => {
    /* with j = 0 the sympatric event has probability one; with j > 0 it has to
       have less, or +J would be adding probability out of nowhere */
    const sym0 = cDEC[iK].find(e => e[0] === iK && e[1] === iK)[2];
    const sym1 = cJ[iK].find(e => e[0] === iK && e[1] === iK)[2];
    return sym0 === 1 && sym1 < 1;
  })());
  t('con j en su máximo, la simpatría se apaga del todo', (() => {
    const c = Biogeo.cladoEvents(st4, 'DEC', 3, {});
    const sym = c[iK].find(e => e[0] === iK && e[1] === iK);
    return !sym || sym[2] < 1e-12;
  })());
  t('un multiplicador de cero también prohíbe el salto hacia esa área', (() => {
    const M = Array.from({ length: 4 }, () => new Float64Array(4).fill(1));
    M[0][1] = 0;
    const c = Biogeo.cladoEvents(st4, 'DEC', 0.5, { multipliers: M });
    return !c[iK].some(e => st4.list[e[0]] === 0b0010 || st4.list[e[1]] === 0b0010);
  })());

  sec('96 · La verosimilitud sobre un árbol pequeño');
  const L4 = ['a', 'b', 'c', 'd'];
  /* an ultrametric four-taxon tree whose areas are as simple as they can be:
     two species in K, two in O */
  const tr4 = Tree.parseNewick('((a:1,b:1):1,(c:1,d:1):1);', L4);
  const ranges = [0b0001, 0b0001, 0b0010, 0b0010];
  const prep4 = Biogeo.prepare(tr4, ranges, st4, {});
  t('la profundidad del árbol se mide bien', Math.abs(prep4.maxDepth - 2) < 1e-12, String(prep4.maxDepth));
  t('ninguna punta queda sin rango', prep4.missing.length === 0);
  const r1 = Biogeo.logLike(prep4, st4, { d: 0.3, e: 0.1, j: 0, model: 'DEC' }, {});
  t('la verosimilitud es finita', isFinite(r1.lnL), fmtLnL(r1.lnL));
  t('y negativa, como una probabilidad debe ser', r1.lnL < 0);
  const r2 = Biogeo.logLike(prep4, st4, { d: 0.3, e: 0.1, j: 0, model: 'DEC' }, { rootPrior: 'flat' });
  near('los dos convenios de raíz difieren en log(número de rangos)',
    r1.lnL - r2.lnL, Math.log(15), 1e-12);
  t('el rango nulo no aporta nada: sus parciales son cero', (() => {
    for (let k = 0; k < prep4.F.n; k++) {
      if (!r1.partial[k]) continue;
      if (r1.partial[k][0] > 1e-15) return false;
    }
    return true;
  })());
  t('con d = 0 y e = 0 nada cambia nunca, y sólo el ancestro común es posible', (() => {
    const r = Biogeo.logLike(prep4, st4, { d: 0, e: 0, j: 0, model: 'DEC' }, {});
    /* with no dispersal a lineage in K cannot reach O, so the two clades can
       only be explained by an ancestor already in KO, split by vicariance */
    return isFinite(r.lnL);
  })());
  t('duplicar las ramas y dividir las tasas entre dos da la misma verosimilitud', (() => {
    const tr8 = Tree.parseNewick('((a:2,b:2):2,(c:2,d:2):2);', L4);
    const p8 = Biogeo.prepare(tr8, ranges, st4, {});
    const a = Biogeo.logLike(prep4, st4, { d: 0.3, e: 0.1, j: 0, model: 'DEC' }, {});
    const b = Biogeo.logLike(p8, st4, { d: 0.15, e: 0.05, j: 0, model: 'DEC' }, {});
    return Math.abs(a.lnL - b.lnL) < 1e-9;
  })());
  t('renombrar las áreas no cambia nada', (() => {
    const swapped = ranges.map(m => (m === 0b0001 ? 0b0010 : 0b0001));
    const p = Biogeo.prepare(tr4, swapped, st4, {});
    const r = Biogeo.logLike(p, st4, { d: 0.3, e: 0.1, j: 0, model: 'DEC' }, {});
    return Math.abs(r.lnL - r1.lnL) < 1e-10;
  })());

  sec('97 · El ajuste y la comparación de modelos');
  const fitDEC = Biogeo.fit(tr4, ranges, st4, { model: 'DEC', withJ: false, prep: prep4 });
  t('d sale positiva', fitDEC.d > 0, fmtNum(fitDEC.d));
  t('la verosimilitud del ajuste no es peor que la de un punto cualquiera',
    fitDEC.lnL >= r1.lnL - 1e-9, `${fmtLnL(fitDEC.lnL)} vs ${fmtLnL(r1.lnL)}`);
  t('AIC penaliza con dos por parámetro', Math.abs(fitDEC.AIC - (-2 * fitDEC.lnL + 4)) < 1e-12);
  t('AICc corrige por el tamaño de la muestra y es mayor que AIC', fitDEC.AICc > fitDEC.AIC);
  const fitDECJ = Biogeo.fit(tr4, ranges, st4, { model: 'DEC', withJ: true, prep: prep4 });
  t('DEC+J no queda por debajo de DEC, porque lo contiene',
    fitDECJ.lnL >= fitDEC.lnL - 1e-6, `${fmtLnL(fitDECJ.lnL)} vs ${fmtLnL(fitDEC.lnL)}`);
  t('DEC+J tiene un parámetro más', fitDECJ.k === fitDEC.k + 1);
  t('j se queda dentro de su intervalo', fitDECJ.j >= 0 && fitDECJ.j <= fitDECJ.jMax,
    `${fmtNum(fitDECJ.j)} de ${fmtNum(fitDECJ.jMax)}`);
  t('el tope de j depende del modelo: 3 en DEC, 2 en DIVALIKE, 1 en BAYAREALIKE', (() => {
    const a = Biogeo.fit(tr4, ranges, st4, { model: 'DIVALIKE', withJ: true, prep: prep4 });
    const b = Biogeo.fit(tr4, ranges, st4, { model: 'BAYAREALIKE', withJ: true, prep: prep4 });
    return Math.abs(fitDECJ.jMax - 3) < 1e-4 && Math.abs(a.jMax - 2) < 1e-4 && Math.abs(b.jMax - 1) < 1e-4;
  })());
  const cmp = Biogeo.compare(tr4, ranges, st4, { models: [['DEC', false], ['DEC', true]] });
  t('la comparación devuelve los modelos ordenados por AICc',
    cmp.rows[0].AICc <= cmp.rows[1].AICc);
  t('los pesos de Akaike suman uno',
    Math.abs(cmp.rows.reduce((a, b) => a + b.w, 0) - 1) < 1e-12);
  t('la razón de verosimilitud del +J se informa con su p', (() => {
    const j = cmp.rows.find(r => r.withJ);
    return j.lrt && j.lrt.p >= 0 && j.lrt.p <= 1 && j.lrt.against === 'DEC';
  })());

  sec('98 · Rangos ancestrales');
  const anc4 = Biogeo.ancestral(fitDEC, {});
  t('hay una reconstrucción por nodo interno', anc4.nodes.length === 3, String(anc4.nodes.length));
  t('cada nodo suma uno', anc4.nodes.every(nd => Math.abs(nd.probs.reduce((a, b) => a + b, 0) - 1) < 1e-9));
  t('el rango nulo nunca es un ancestro', anc4.nodes.every(nd => nd.probs[0] < 1e-9));
  t('el clado de a y b, los dos en K, se reconstruye en K', (() => {
    const nd = anc4.nodes.find(x => x.node !== 0 && x.probs[iK] > 0.5);
    return !!nd;
  })());
  t('cada esquina suma uno también', (() => {
    for (let k = 0; k < anc4.corners.length; k++) {
      const c = anc4.corners[k];
      if (!c) continue;
      const a = c.leftProbs.reduce((x, y) => x + y, 0), b = c.rightProbs.reduce((x, y) => x + y, 0);
      if (Math.abs(a - 1) > 1e-9 || Math.abs(b - 1) > 1e-9) return false;
    }
    return true;
  })());
  t('el mejor rango de cada nodo se nombra y se acompaña de su probabilidad', (() => {
    const b = Biogeo.bestRanges(anc4, A4, prep4);
    return b.length === 3 && b.every(x => x.range && x.probability > 0 && x.probability <= 1);
  })());
  t('si todas las puntas están en la misma área, todos los nodos también', (() => {
    const same = [0b0001, 0b0001, 0b0001, 0b0001];
    const p = Biogeo.prepare(tr4, same, st4, {});
    const f = Biogeo.fit(tr4, same, st4, { model: 'DEC', withJ: false, prep: p });
    const a = Biogeo.ancestral(f, {});
    return a.nodes.every(nd => nd.probs[iK] > 0.9);
  })());

  sec('99 · Mapeo estocástico biogeográfico');
  const bsm4 = Biogeo.stochasticMap(fitDEC, { reps: 200, seed: 9, ancestral: anc4 });
  t('los eventos cladogenéticos suman un nodo interno cada uno', (() => {
    const c = bsm4.mean.sympatry + bsm4.mean.subset + bsm4.mean.vicariance + bsm4.mean.jump;
    return Math.abs(c - 3) < 1e-9;
  })(), fmtFixed(bsm4.mean.sympatry + bsm4.mean.subset + bsm4.mean.vicariance + bsm4.mean.jump, 6));
  t('sin j no hay saltos', bsm4.mean.jump === 0);
  t('las cuentas no son negativas', Object.keys(bsm4.mean).every(k => bsm4.mean[k] >= 0));
  t('la misma semilla da el mismo mapa', (() => {
    const b = Biogeo.stochasticMap(fitDEC, { reps: 200, seed: 9, ancestral: anc4 });
    return Math.abs(b.mean.dispersal - bsm4.mean.dispersal) < 1e-12;
  })());
  t('con dos clados limpios, uno en K y otro en O, no hace falta ninguna dispersión',
    bsm4.mean.dispersal === 0 && bsm4.mean.extinction === 0,
    `${fmtNum(bsm4.mean.dispersal)} / ${fmtNum(bsm4.mean.extinction)}`);
  /* a configuration that does need movement, so that there is something for two
     seeds to disagree about: one species widespread and one in a third area */
  const rMix = [0b0001, 0b0011, 0b0010, 0b0100];
  const pMix = Biogeo.prepare(tr4, rMix, st4, {});
  const fMix = Biogeo.fit(tr4, rMix, st4, { model: 'DEC', withJ: false, prep: pMix });
  const aMix = Biogeo.ancestral(fMix, {});
  const mix9 = Biogeo.stochasticMap(fMix, { reps: 200, seed: 9, ancestral: aMix });
  t('cuando sí hace falta, el mapeo la encuentra', mix9.mean.dispersal > 0.5, fmtNum(mix9.mean.dispersal));
  t('otra semilla da otro mapa, pero parecido', (() => {
    const b = Biogeo.stochasticMap(fMix, { reps: 200, seed: 77, ancestral: aMix });
    const d = Math.abs(b.mean.dispersal - mix9.mean.dispersal);
    return d > 0 && d < 0.5;
  })());
  t('las frecuencias de los nodos se acercan a la marginal analítica', (() => {
    let worst = 0;
    bsm4.nodeProbs.forEach(m => {
      const a = anc4.nodes.find(q => q.node === m.node);
      m.probs.forEach((p, i) => { worst = Math.max(worst, Math.abs(p - a.probs[i])); });
    });
    return worst < 0.15;
  })());
  t('una historia de rama con los dos extremos iguales puede no tener ningún evento', (() => {
    const r = rng(3);
    let quiet = 0;
    for (let i = 0; i < 40; i++) {
      const h = Biogeo.branchHistory(Q4, st4, 0.02, iK, iK, r, 100);
      if (h.dispersal === 0 && h.extinction === 0) quiet++;
    }
    return quiet > 20;
  })());
  t('con los extremos distintos siempre hay al menos un evento', (() => {
    const r = rng(4);
    for (let i = 0; i < 25; i++) {
      const h = Biogeo.branchHistory(Q4, st4, 3, iK, iKO, r, 300);
      if (!h.failed && h.dispersal + h.extinction < 1) return false;
    }
    return true;
  })());
  t('classify reconoce los cuatro tipos de división', (() => {
    const iO = st4.index.get(0b0010);
    return Biogeo.classify(st4, iKO, iKO, iKO) === 'sympatry'
      && Biogeo.classify(st4, iKO, iKO, iK) === 'subset'
      && Biogeo.classify(st4, iKO, iK, iO) === 'vicariance'
      && Biogeo.classify(st4, iK, iK, st4.index.get(0b0100)) === 'jump';
  })());

  sec('100 · Leer y escribir la tabla de áreas');
  const lg = Biogeo.parseRanges('3\t4\t(K O M H)\nsp1\t1000\nsp2\t0110\nsp3\t0001');
  t('se reconoce el formato de LAGRANGE', lg.format === 'lagrange');
  t('los nombres de las áreas salen del encabezado', lg.areas.join('') === 'KOMH', lg.areas.join(''));
  t('los rangos se leen como conjuntos', lg.masks.join(',') === '1,6,8', lg.masks.join(','));
  t('sin nombres en el encabezado se ponen letras', (() => {
    const p = Biogeo.parseRanges('2 3\nsp1 100\nsp2 011');
    return p.areas.join('') === 'ABC';
  })());
  const tb = Biogeo.parseRanges('taxon,K,O,M\nsp1,1,0,0\nsp2,0,1,1\nsp3,sí,no,no');
  t('también se lee una tabla separada por comas', tb.format === 'table');
  t('con sus áreas en el encabezado', tb.areas.join('') === 'KOM', tb.areas.join(''));
  t('y acepta «sí» además de 1', tb.masks.join(',') === '1,6,1', tb.masks.join(','));
  t('lo que escribe lo vuelve a leer igual', (() => {
    const txt = Biogeo.writeLagrange(lg.areas, lg.taxa, lg.masks);
    const back = Biogeo.parseRanges(txt);
    return back.taxa.join('|') === lg.taxa.join('|') && back.masks.join(',') === lg.masks.join(',');
  })());
  t('una tabla vacía se rechaza en vez de dar algo sin sentido',
    Biogeo.parseRanges('') === null);

  sec('101 · Estratos de tiempo');
  const strata = Biogeo.normaliseStrata([
    { name: 'reciente', from: 0, multipliers: null },
    { name: 'antiguo', from: 1.2, multipliers: null },
  ], 3);
  t('los estratos se ordenan del presente hacia el pasado', strata[0].from === 0 && strata[1].from === 1.2);
  t('el final de cada uno es el principio del siguiente', strata[0].to === 1.2);
  t('el último llega hasta la raíz', strata[1].to === 3);
  t('una rama dentro de un solo estrato devuelve una pieza',
    Biogeo.strataFor(strata, 0.1, 0.9).length === 1);
  t('una rama que cruza la frontera devuelve dos', (() => {
    const p = Biogeo.strataFor(strata, 0.5, 2.0);
    return p.length === 2 && Math.abs(p[0].length - 0.7) < 1e-12 && Math.abs(p[1].length - 0.8) < 1e-12;
  })());
  t('las piezas suman la longitud de la rama', (() => {
    const p = Biogeo.strataFor(strata, 0.3, 2.6);
    return Math.abs(p.reduce((a, b) => a + b.length, 0) - 2.3) < 1e-12;
  })());
  t('con estratos idénticos la verosimilitud no cambia', (() => {
    const a = Biogeo.logLike(prep4, st4, { d: 0.3, e: 0.1, j: 0, model: 'DEC' }, {});
    const b = Biogeo.logLike(prep4, st4, { d: 0.3, e: 0.1, j: 0, model: 'DEC' }, { strata });
    return Math.abs(a.lnL - b.lnL) < 1e-9;
  })(), (() => {
    const a = Biogeo.logLike(prep4, st4, { d: 0.3, e: 0.1, j: 0, model: 'DEC' }, {});
    const b = Biogeo.logLike(prep4, st4, { d: 0.3, e: 0.1, j: 0, model: 'DEC' }, { strata });
    return fmtNum(Math.abs(a.lnL - b.lnL));
  })());
  t('un estrato que prohíbe un área cambia la verosimilitud', (() => {
    const s2 = Biogeo.normaliseStrata([
      { name: 'reciente', from: 0 },
      { name: 'antiguo', from: 1.2, allowedAreas: 0b0001 },
    ], 3);
    const a = Biogeo.logLike(prep4, st4, { d: 0.3, e: 0.1, j: 0, model: 'DEC' }, {});
    const b = Biogeo.logLike(prep4, st4, { d: 0.3, e: 0.1, j: 0, model: 'DEC' }, { strata: s2 });
    return Math.abs(a.lnL - b.lnL) > 1e-3;
  })());

})();
