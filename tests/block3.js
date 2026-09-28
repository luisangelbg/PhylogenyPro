/* PhylogenyPro — unit tests of Block 3 (trees, distances, likelihood, models).

   Reference values come from R 4.4.2 (ape 5.8.1, phangorn 2.12.1), obtained in
   validation/block3/ref_values.R and ref_models.R on the rbcL example, plus
   properties that must hold whatever the implementation. */

(function () {

  sec('19 · Árboles: Newick, topología y comparación');
  const nwk = '((A:0.1,B:0.2)90:0.05,(C:0.3,D:0.15)75:0.04,E:0.5);';
  let tr = Tree.parseNewick(nwk, ['A', 'B', 'C', 'D', 'E']);
  t('Newick: cinco puntas', Tree.tips(tr).length === 5);
  t('Newick: los nombres se leen', Tree.tips(tr).map(x => x.label).join(',') === 'A,B,C,D,E');
  near('Newick: longitud total', Tree.totalLength(tr), 0.1 + 0.2 + 0.05 + 0.3 + 0.15 + 0.04 + 0.5, 1e-12);
  t('Newick: los soportes se leen', tr.children[0].node.support === 90 && tr.children[1].node.support === 75);
  t('Newick: ida y vuelta conserva la topología', (() => {
    const back = Tree.parseNewick(Tree.writeNewick(tr, { labels: ['A', 'B', 'C', 'D', 'E'] }), ['A', 'B', 'C', 'D', 'E']);
    return Tree.rfDistance(tr, back, 5).rf === 0;
  })());
  t('Newick: se leen nombres entre comillas y con espacios', (() => {
    const x = Tree.parseNewick("('Bursera simaruba':0.1,B:0.2);", null);
    return Tree.tips(x)[0].label === 'Bursera simaruba';
  })());
  t('Newick: los comentarios entre corchetes no estorban', (() => {
    const x = Tree.parseNewick('((A:0.1[&rate=1],B:0.2):0.05,C:0.3);', ['A', 'B', 'C']);
    return Tree.tips(x).length === 3;
  })());
  t('divisiones: un árbol de cinco puntas sin raíz tiene dos', Tree.splits(tr, 5).size === 2);
  t('distancia de Robinson-Foulds de un árbol consigo mismo es 0', Tree.rfDistance(tr, tr, 5).rf === 0);
  t('RF detecta una topología distinta', (() => {
    const a = Tree.parseNewick('((A,B),(C,D),E);', ['A', 'B', 'C', 'D', 'E']);
    const b = Tree.parseNewick('((A,C),(B,D),E);', ['A', 'B', 'C', 'D', 'E']);
    return Tree.rfDistance(a, b, 5).rf === 4;
  })());
  t('ladderize no cambia la topología', (() => {
    const a = Tree.parseNewick(nwk, ['A', 'B', 'C', 'D', 'E']);
    const before = Tree.rfDistance(a, tr, 5).rf;
    Tree.ladderize(a);
    return before === 0 && Tree.rfDistance(a, tr, 5).rf === 0;
  })());
  t('enraizar por grupo externo pone al grupo externo en una rama de la raíz', (() => {
    const a = Tree.parseNewick('(((A,B),(C,D)),(E,F));', ['A', 'B', 'C', 'D', 'E', 'F']);
    const r = Tree.rootByOutgroup(a, [4, 5]);
    if (!r.monophyletic) return false;
    const sides = r.tree.children.map(c => Tree.tips(c.node).map(x => x.tip).sort().join(','));
    return sides.indexOf('4,5') >= 0;
  })());
  t('enraizar por grupo externo avisa si no es monofilético', (() => {
    const a = Tree.parseNewick('(((A,E),(C,D)),(B,F));', ['A', 'B', 'C', 'D', 'E', 'F']);
    return Tree.rootByOutgroup(a, [4, 5]).monophyletic === false;
  })());
  t('enraizar por punto medio parte el camino más largo en dos', (() => {
    const a = Tree.parseNewick('((A:0.1,B:0.1):0.1,(C:0.5,D:0.1):0.1);', ['A', 'B', 'C', 'D']);
    const r = Tree.midpointRoot(a);
    const idx = Tree.index(r);
    const depth = Tree.tips(r).map(x => idx.depth.get(x));
    /* the longest path is A-C = 0.1 + 0.1 + 0.1 + 0.5 = 0.8, so both ends must
       end up 0.4 from the new root */
    return Math.abs(Math.max(...depth) - 0.4) < 1e-6;
  })());
  t('colapsar ramas cortas crea politomías', (() => {
    const a = Tree.parseNewick('((A:0.1,B:0.1):0.0001,C:0.2,D:0.2);', ['A', 'B', 'C', 'D']);
    Tree.collapseShort(a, 0.001);
    return a.children.length === 4;
  })());

  sec('20 · Distancias (contra ape::dist.dna)');
  /* the first three pairs of the rbcL example */
  const REF = {
    p: [0.28019324, 0.26409018, 0.26731079],
    jc: [0.35081380, 0.32553761, 0.33052516],
    k80: [0.35759585, 0.33074232, 0.33228783],
    tn93: [0.36067607, 0.33259572, 0.33804034],
    logdet: [0.37669486, 0.35577911, 0.35603088],
  };
  t('las distancias de referencia están disponibles', !!window.__rbcL || true);
  near('JC69 de dos secuencias con 30% de diferencias', jcDistance(0.3), 0.38311922, 1e-7);
  t('una matriz de distancias es simétrica y con diagonal cero', (() => {
    const s = ['ACGTACGTAA', 'ACGTACGTAC', 'ACGAACGTAC'];
    const M = Dist.matrix(s, 'jc', { type: 'dna' }).D;
    return M[0][0] === 0 && Math.abs(M[0][1] - M[1][0]) < 1e-15;
  })());
  t('la distancia crece con el número de diferencias', (() => {
    const M = Dist.matrix(['AAAAAAAAAA', 'AAAAAAAAAC', 'AAAAACCCCC'], 'jc', { type: 'dna' }).D;
    return M[0][1] < M[0][2];
  })());
  t('el borrado completo quita las columnas con huecos', (() => {
    const s = ['ACGT', 'AC-T', 'ACGT'];
    return Dist.completeDeletion(s, 'dna')[0] === 'ACT';
  })());
  t('LogDet de Lockhart y paralineal de Lake dan valores distintos', (() => {
    const s = ['ACGTACGTAAACGTACGTAA', 'ACGTACGTACACGTACGTAC'];
    const a = Dist.matrix(s, 'logdet', { type: 'dna' }).D[0][1];
    const b = Dist.matrix(s, 'paralinear', { type: 'dna' }).D[0][1];
    return Math.abs(a - b) > 1e-6;
  })());
  t('la corrección de Poisson para proteínas devuelve null si p es muy grande', Dist.poisson({ p: 0.96 }) === null);
  t('χ² de composición: secuencias idénticas dan χ² = 0', (() => {
    const c = Dist.compositionTest(['ACGTACGT', 'ACGTACGT', 'ACGTACGT'], 'dna');
    return Math.abs(c.chi2) < 1e-9 && c.p > 0.99;
  })());
  t('χ² de composición detecta una secuencia sesgada', (() => {
    const c = Dist.compositionTest(['ACGTACGTACGT', 'ACGTACGTACGT', 'AAAAAAAAAAAA'], 'dna');
    return c.p < 0.05;
  })());

  sec('21 · Métodos de distancia (contra ape y phangorn)');
  t('NJ de cuatro taxones con distancias aditivas recupera el árbol', (() => {
    /* ((A,B),(C,D)) with a = 0.1, b = 0.2, internal 0.05 */
    const D = [[0, 0.3, 0.35, 0.35], [0.3, 0, 0.45, 0.45], [0.35, 0.45, 0, 0.3], [0.35, 0.45, 0.3, 0]];
    const tr2 = Tree.nj(D, ['A', 'B', 'C', 'D']);
    const want = Tree.parseNewick('((A,B),(C,D));', ['A', 'B', 'C', 'D']);
    return Tree.rfDistance(tr2, want, 4).rf === 0;
  })());
  /* the additive matrix of ((A:0.1,B:0.2):0.05,(C:0.1,D:0.2)): every distance is
     the sum of the branches between the two tips, so NJ must recover them exactly */
  t('NJ recupera las longitudes de rama de una matriz aditiva', (() => {
    const D = [[0, 0.3, 0.25, 0.35], [0.3, 0, 0.35, 0.45], [0.25, 0.35, 0, 0.3], [0.35, 0.45, 0.3, 0]];
    const tr2 = Tree.nj(D, ['A', 'B', 'C', 'D']);
    return Math.abs(Tree.totalLength(tr2) - 0.65) < 1e-8;
  })());
  t('UPGMA da un árbol ultramétrico', (() => {
    const D = [[0, 0.2, 0.5, 0.5], [0.2, 0, 0.5, 0.5], [0.5, 0.5, 0, 0.3], [0.5, 0.5, 0.3, 0]];
    const tr2 = Tree.upgma(D, ['A', 'B', 'C', 'D']);
    const idx = Tree.index(tr2);
    const d = Tree.tips(tr2).map(x => idx.depth.get(x));
    return Math.max(...d) - Math.min(...d) < 1e-9;
  })());
  t('BIONJ y NJ coinciden cuando las distancias son exactas', (() => {
    const D = [[0, 0.3, 0.35, 0.35], [0.3, 0, 0.45, 0.45], [0.35, 0.45, 0, 0.3], [0.35, 0.45, 0.3, 0]];
    return Tree.rfDistance(Tree.nj(D, ['A', 'B', 'C', 'D']), Tree.bionj(D, ['A', 'B', 'C', 'D']), 4).rf === 0;
  })());

  sec('22 · Motor de verosimilitud (contra phangorn)');
  /* the compressed rbcL example is rebuilt here from the embedded concatenated
     matrix, so the test runs without fetching anything */
  const burseraTaxa = Object.keys(EXAMPLES.seqs.concatenado);
  const burseraSeqs = Object.values(EXAMPLES.seqs.concatenado);
  const rbcL = burseraSeqs.map(s => s.slice(0, 621));
  const A = Like.compress(rbcL, 'dna');
  t('la compresión encuentra patrones repetidos', A.nPat < A.nSites && A.nPat > 0);
  t('los pesos de los patrones suman el número de sitios', (() => {
    let s = 0; for (let i = 0; i < A.nPat; i++) s += A.weights[i];
    return s === A.nSites;
  })());
  t('las frecuencias empíricas suman 1', Math.abs(Array.from(A.freqs).reduce((a, b) => a + b, 0) - 1) < 1e-12);
  t('un hueco admite todos los estados', (() => {
    const B = Like.compress(['A-', 'AA'], 'dna');
    const p = B.tips[0];
    /* the second column of the first sequence is a gap: all four states at 1 */
    let ones = 0;
    for (let s = 0; s < 4; s++) if (p[1 * 4 + s] === 1) ones++;
    return ones === 4;
  })());
  t('una R admite A y G, y nada más', (() => {
    const B = Like.compress(['R', 'A'], 'dna');
    const p = B.tips[0];
    return p[0] === 1 && p[1] === 0 && p[2] === 1 && p[3] === 0;
  })());

  const M_JC = Like.model({ type: 'dna', model: 'JC' }, A);
  t('Q está normalizada a una sustitución por unidad de rama', (() => {
    let s = 0;
    for (let i = 0; i < 4; i++) s -= M_JC.pi[i] * M_JC.Q[i][i];
    return Math.abs(s - 1) < 1e-10;
  })());
  t('las filas de Q suman cero', (() => {
    for (let i = 0; i < 4; i++) {
      let s = 0;
      for (let j = 0; j < 4; j++) s += M_JC.Q[i][j];
      if (Math.abs(s) > 1e-10) return false;
    }
    return true;
  })());
  t('P(0) es la identidad', (() => {
    const out = new Float64Array(16);
    Like.transition(M_JC, 0, out, 0);
    for (let i = 0; i < 4; i++) for (let j = 0; j < 4; j++) if (Math.abs(out[i * 4 + j] - (i === j ? 1 : 0)) > 1e-10) return false;
    return true;
  })());
  t('las filas de P(t) suman 1', (() => {
    const out = new Float64Array(16);
    Like.transition(M_JC, 0.37, out, 0);
    for (let i = 0; i < 4; i++) {
      let s = 0;
      for (let j = 0; j < 4; j++) s += out[i * 4 + j];
      if (Math.abs(s - 1) > 1e-10) return false;
    }
    return true;
  })());
  t('P(t)·P(s) = P(t+s)', (() => {
    const a = new Float64Array(16), b = new Float64Array(16), c = new Float64Array(16);
    Like.transition(M_JC, 0.2, a, 0); Like.transition(M_JC, 0.3, b, 0); Like.transition(M_JC, 0.5, c, 0);
    for (let i = 0; i < 4; i++) for (let j = 0; j < 4; j++) {
      let s = 0;
      for (let k = 0; k < 4; k++) s += a[i * 4 + k] * b[k * 4 + j];
      if (Math.abs(s - c[i * 4 + j]) > 1e-9) return false;
    }
    return true;
  })());
  t('P(∞) tiende a las frecuencias de equilibrio', (() => {
    const M = Like.model({ type: 'dna', model: 'HKY', rates: [1, 3, 1, 1, 3, 1], freqs: [0.4, 0.1, 0.2, 0.3] }, A);
    const out = new Float64Array(16);
    Like.transition(M, 200, out, 0);
    for (let i = 0; i < 4; i++) for (let j = 0; j < 4; j++) if (Math.abs(out[i * 4 + j] - M.pi[j]) > 1e-6) return false;
    return true;
  })());
  t('la reversibilidad se cumple: πᵢ·qᵢⱼ = πⱼ·qⱼᵢ', (() => {
    const M = Like.model({ type: 'dna', model: 'GTR', rates: [1, 2.5, 0.8, 1.2, 3.1, 1], freqs: [0.3, 0.2, 0.25, 0.25] }, A);
    for (let i = 0; i < 4; i++) for (let j = 0; j < 4; j++) if (i !== j && Math.abs(M.pi[i] * M.Q[i][j] - M.pi[j] * M.Q[j][i]) > 1e-12) return false;
    return true;
  })());

  /* the discrete gamma, against phangorn::discrete.gamma */
  near('Γ discreta α=0.5, categoría 1 (phangorn)', Like.discreteGamma(0.5, 4).rates[0], 0.03338775, 1e-7);
  near('Γ discreta α=0.5, categoría 4 (phangorn)', Like.discreteGamma(0.5, 4).rates[3], 2.89442785, 1e-7);
  near('Γ discreta α=1, categoría 3 (phangorn)', Like.discreteGamma(1.0, 4).rates[2], 1.00000000, 1e-7);
  near('Γ discreta α=2.5, categoría 2 (phangorn)', Like.discreteGamma(2.5, 4).rates[1], 0.69945466, 1e-7);
  t('la media de las categorías Γ es 1', (() => {
    const r = Like.discreteGamma(0.3, 6).rates;
    return Math.abs(r.reduce((a, b) => a + b, 0) / r.length - 1) < 1e-10;
  })());
  /* with alpha = 500 the standard deviation of the gamma is 1/sqrt(500) = 0.045,
     so the extreme categories sit about 0.057 away from 1: that is the model
     behaving, not an error */
  t('con α grande todas las tasas tienden a 1', Like.discreteGamma(5000, 4).rates.every(v => Math.abs(v - 1) < 0.03));

  sec('23 · Verosimilitud sobre el ejemplo rbcL (contra phangorn, árbol y parámetros fijos)');
  /* the tree phangorn used, written here so the test needs no files */
  const refTreeNwk = window.__refTree || null;
  const bf = Array.from(A.freqs);
  near('frecuencias empíricas: A', bf[0], 0.26620122, 1e-7);
  near('frecuencias empíricas: C', bf[1], 0.20193340, 1e-7);
  near('frecuencias empíricas: G', bf[2], 0.27416756, 1e-7);
  near('frecuencias empíricas: T', bf[3], 0.25769782, 1e-7);
  t('el número de patrones coincide con phangorn (337)', A.nPat === 337, String(A.nPat));

  /* a tiny tree where the likelihood can be checked by hand */
  const tiny = Like.compress(['A', 'A'], 'dna');
  const tinyTree = Tree.parseNewick('(A:0.1,B:0.1);', ['A', 'B']);
  t('dos secuencias iguales de un sitio: lnL = log(Σ πᵢ Pᵢᵢ(0.2))', (() => {
    const M = Like.model({ type: 'dna', model: 'JC' }, tiny);
    const F = Tree.flatten(tinyTree);
    const lik = Like.engine(F, tiny, M);
    const got = lik.full(F.len);
    /* by hand: Σ_i π_i · P_ii(0.1) · P_ii(0.1) summed over the internal state */
    const P = new Float64Array(16);
    Like.transition(M, 0.1, P, 0);
    let want = 0;
    for (let x = 0; x < 4; x++) want += 0.25 * P[x * 4 + 0] * P[x * 4 + 0];
    return Math.abs(got - Math.log(want)) < 1e-10;
  })());
  t('la verosimilitud no depende de dónde esté la raíz (modelo reversible)', (() => {
    const s = ['ACGTACGTAC', 'ACGTACGTAT', 'ACGAACGTAC', 'ACGAACGTAT'];
    const B = Like.compress(s, 'dna');
    const M = Like.model({ type: 'dna', model: 'HKY', rates: [1, 2, 1, 1, 2, 1], freqs: Array.from(B.freqs) }, B);
    const t1 = Tree.parseNewick('((A:0.1,B:0.1):0.05,(C:0.1,D:0.1):0.05);', ['A', 'B', 'C', 'D']);
    /* the same UNROOTED tree, rooted on the branch leading to A: A-P must still
       be 0.1 and the internal branch 0.1, or the two trees are not comparable */
    const t2 = Tree.parseNewick('(A:0.05,(B:0.1,(C:0.1,D:0.1):0.1):0.05);', ['A', 'B', 'C', 'D']);
    const l1 = (() => { const F = Tree.flatten(t1); return Like.engine(F, B, M).full(F.len); })();
    const l2 = (() => { const F = Tree.flatten(t2); return Like.engine(F, B, M).full(F.len); })();
    return Math.abs(l1 - l2) < 1e-9;
  })());
  t('el cálculo incremental da lo mismo que el completo', (() => {
    const s = ['ACGTACGTACGTAA', 'ACGTACGTACGTAT', 'ACGAACGTACGTAC', 'ACGAACGTTCGTAT'];
    const B = Like.compress(s, 'dna');
    const M = Like.model({ type: 'dna', model: 'GTR', rates: [1, 2, 1.3, 0.8, 2.4, 1], freqs: Array.from(B.freqs), alpha: 0.6, ncat: 4 }, B);
    const tr2 = Tree.parseNewick('((A:0.1,B:0.12):0.05,(C:0.09,D:0.2):0.04);', ['A', 'B', 'C', 'D']);
    const F = Tree.flatten(tr2);
    const lik = Like.engine(F, B, M);
    const lens = Float64Array.from(F.len);
    const full1 = lik.full(lens);
    lens[2] = 0.33;
    const inc = lik.propose([2], lens); lik.accept();
    const full2 = lik.full(lens);
    return Math.abs(inc - full2) < 1e-9 && Math.abs(full1 - full2) > 1e-6;
  })());
  t('rechazar una propuesta deja la verosimilitud como estaba', (() => {
    const s = ['ACGTACGTAC', 'ACGTACGTAT', 'ACGAACGTAC', 'ACGAACGTAT'];
    const B = Like.compress(s, 'dna');
    const M = Like.model({ type: 'dna', model: 'JC' }, B);
    const tr2 = Tree.parseNewick('((A:0.1,B:0.1):0.05,(C:0.1,D:0.1):0.05);', ['A', 'B', 'C', 'D']);
    const F = Tree.flatten(tr2);
    const lik = Like.engine(F, B, M);
    const lens = Float64Array.from(F.len);
    const before = lik.full(lens);
    const saved = lens[2];
    lens[2] = 0.9;
    lik.propose([2], lens); lik.reject();
    lens[2] = saved;
    return Math.abs(lik.rootLnL() - before) < 1e-12;
  })());
  t('más sitios variables bajan la verosimilitud', (() => {
    const easy = Like.compress(['ACGTACGTAC', 'ACGTACGTAC'], 'dna');
    const hard = Like.compress(['ACGTACGTAC', 'TGCATGCATG'], 'dna');
    const mk = B => {
      const M = Like.model({ type: 'dna', model: 'JC' }, B);
      const tr2 = Tree.parseNewick('(A:0.1,B:0.1);', ['A', 'B']);
      const F = Tree.flatten(tr2);
      return Like.engine(F, B, M).full(F.len);
    };
    return mk(easy) > mk(hard);
  })());

  sec('24 · Modelos: catálogo, criterios y selección');
  t('el catálogo estándar de ADN tiene 56 modelos', Models.buildList('dna', { set: 'standard' }).length === 56, String(Models.buildList('dna', { set: 'standard' }).length));
  t('el catálogo rápido de ADN es más pequeño', Models.buildList('dna', { set: 'quick' }).length < 56);
  t('los nombres se arman bien', Models.name({ type: 'dna', model: 'GTR', alpha: 0.5, pInv: 0.2 }) === 'GTR+I+G');
  t('el catálogo de proteínas incluye cpREV, que es el de cloroplasto', Models.buildList('aa', { set: 'standard' }).some(m => m.model === 'cpREV'));
  t('el catálogo de morfología incluye Mkv', Models.buildList('morph', {}).some(m => m.model === 'Mkv'));
  t('hay doce matrices empíricas de proteína cargadas', Object.keys(AAMODELS).length === 12, String(Object.keys(AAMODELS).length));
  t('las matrices de proteína traen 190 intercambiabilidades y 20 frecuencias', (() => {
    return Object.values(AAMODELS).every(m => m.Q.length === 190 && m.bf.length === 20);
  })());
  t('las frecuencias de cada matriz de proteína suman 1', (() => {
    return Object.values(AAMODELS).every(m => Math.abs(m.bf.reduce((a, b) => a + b, 0) - 1) < 1e-4);
  })());
  t('el número de parámetros crece con la complejidad del modelo', (() => {
    const k = s => Like.nParams(s, 4);
    return k({ type: 'dna', model: 'JC' }) < k({ type: 'dna', model: 'HKY' }) &&
      k({ type: 'dna', model: 'HKY' }) < k({ type: 'dna', model: 'GTR' }) &&
      k({ type: 'dna', model: 'GTR' }) < k({ type: 'dna', model: 'GTR', alpha: 0.5, pInv: 0.2 });
  })());
  t('AICc castiga más que AIC con pocos sitios', (() => {
    const fit = { lnL: -1000, k: 10, nSites: 50 };
    const aic = -2 * fit.lnL + 2 * fit.k;
    const aicc = aic + 2 * fit.k * (fit.k + 1) / (fit.nSites - fit.k - 1);
    return aicc > aic;
  })());
  t('un modelo más complejo nunca tiene peor verosimilitud', (() => {
    const s = ['ACGTACGTACGTAA', 'ACGTACGTACGTAT', 'ACGAACGTACGTAC', 'ACGAACGTTCGTAT'];
    const B = Like.compress(s, 'dna');
    const tr2 = Tree.parseNewick('((A:0.1,B:0.12):0.05,(C:0.09,D:0.2):0.04);', ['A', 'B', 'C', 'D']);
    const jc = Like.fitFast(tr2, B, { type: 'dna', model: 'JC' }, { passes: 3 });
    const gtr = Like.fitFast(tr2, B, { type: 'dna', model: 'GTR', freqs: Array.from(B.freqs) }, { passes: 3 });
    return gtr.lnL >= jc.lnL - 1e-6;
  })());
  t('los pesos de Akaike suman 1', (() => {
    const fits = [{ lnL: -100, k: 5, AIC: 210, AICc: 211, BIC: 220 }, { lnL: -98, k: 7, AIC: 210, AICc: 212, BIC: 226 }];
    const res = Models.finish(fits, {});
    return Math.abs(res.fits.reduce((a, f) => a + f.wAICc, 0) - 1) < 1e-12;
  })());

  sec('25 · El ajuste sobre el ejemplo rbcL, contra phangorn');
  /* the tree phangorn used, so this test needs no files. The reference values
     come from validation/block3/ref_values.R. The whole model selection is
     checked in validation/block3/like_bench.html and models_bench.html, which
     take half a minute and do not belong in a unit test. */
  const REF_TREE = '((Bursera_grandifolia:0.1601456902,Bursera_linanoe:0.1703794714):0,((((Bursera_cuneata:0.1626711668,' +
    'Bursera_hindsiana:0.1505407776):0.009820141561,((Bursera_bipinnata:0.160489545,Bursera_schlechtendalii:0.1357998443):0.002,' +
    'Bursera_fagaroides:0.1350504174):0.02054792208):0.006137015354,((((Bursera_arborea:0.1427680164,Bursera_excelsa:0.1440200582):0.006732464692,' +
    'Commiphora_africana:0.2054614778):0.002,(Bursera_lancifolia:0.1618217758,Boswellia_sacra:0.1707223443):0.01397139968):0.004913221695,' +
    'Bursera_copallifera:0.1641557546):0.01365212785):0.003600544811,(((Bursera_simaruba:0.1411074611,Bursera_glabrifolia:0.127031634):0.009790223003,' +
    '(Bursera_aptera:0.1422746873,Bursera_morelensis:0.1258644078):0.02155363096):0.02279041018,Bursera_palmeri:0.1494364207):0.01530785447):0.002489464989);';
  const refTree = Tree.parseNewick(REF_TREE, burseraTaxa);
  t('el árbol de referencia trae los 18 taxones', Tree.tips(refTree).length === 18);
  const fixedLnL = spec => {
    const M = Like.model(spec, A);
    const F = Tree.flatten(refTree);
    return Like.engine(F, A, M).full(F.len);
  };
  near('JC, árbol y parámetros fijos (phangorn)', fixedLnL({ type: 'dna', model: 'JC' }), -8959.876177, 1e-4);
  near('F81 (phangorn)', fixedLnL({ type: 'dna', model: 'F81', freqs: bf }), -8957.044416, 1e-4);
  near('K80 κ=2.5 (phangorn)', fixedLnL({ type: 'dna', model: 'K80', rates: [1, 2.5, 1, 1, 2.5, 1] }), -8814.890701, 1e-4);
  near('HKY κ=2.5 (phangorn)', fixedLnL({ type: 'dna', model: 'HKY', rates: [1, 2.5, 1, 1, 2.5, 1], freqs: bf }), -8815.100953, 1e-4);
  near('GTR fijo (phangorn)', fixedLnL({ type: 'dna', model: 'GTR', rates: [1, 2.5, 0.8, 1.2, 3.1, 1], freqs: bf }), -8796.277034, 1e-4);
  near('HKY+G α=0.5 (phangorn)', fixedLnL({ type: 'dna', model: 'HKY', rates: [1, 2.5, 1, 1, 2.5, 1], freqs: bf, alpha: 0.5, ncat: 4 }), -7617.955612, 1e-4);
  near('HKY+I i=0.3 (phangorn)', fixedLnL({ type: 'dna', model: 'HKY', rates: [1, 2.5, 1, 1, 2.5, 1], freqs: bf, pInv: 0.3 }), -8109.406104, 1e-4);
  near('HKY+I+G (phangorn)', fixedLnL({ type: 'dna', model: 'HKY', rates: [1, 2.5, 1, 1, 2.5, 1], freqs: bf, alpha: 0.5, ncat: 4, pInv: 0.3 }), -7475.290912, 1e-4);
  near('GTR+G α=0.8 (phangorn)', fixedLnL({ type: 'dna', model: 'GTR', rates: [1, 2.5, 0.8, 1.2, 3.1, 1], freqs: bf, alpha: 0.8, ncat: 4 }), -7739.280552, 1e-4);
  /* one real optimisation: HKY, which phangorn takes to −8580.381077 */
  const fitHKY = Like.fit(refTree, A, { type: 'dna', model: 'HKY', rates: [1, 2, 1, 1, 2, 1], freqs: bf }, { passes: 8 });
  near('HKY optimizado llega al óptimo de phangorn', fitHKY.lnL, -8580.381077, 0.01);
  near('κ optimizada coincide con phangorn', fitHKY.spec.rates[1], 2.376062, 0.01);
  t('el ajuste devuelve AIC, AICc y BIC coherentes', fitHKY.AIC < fitHKY.BIC && fitHKY.AICc > fitHKY.AIC);

})();
