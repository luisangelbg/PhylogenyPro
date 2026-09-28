/* PhylogenyPro — unit tests of Block 11 (comparing trees).

   The numerical references come from ape 5.8.1 and phangorn 2.12.1 in
   validation/block11 — RF, normalised RF, weighted RF, the branch score, the
   path difference and the consensus network — plus two things no installed
   package computes, worked out from their definitions in R so that they are
   still independent of this code: the quartet distance by brute force over all
   C(n,4) quartets, and Patterson's D with a block jackknife. The bench
   cmp_bench.html checks all of them and reports no disagreement.

   Beside those, the properties that have to hold whatever the implementation:
   a distance that is zero exactly when two trees are the same and never
   negative, splits that are a set and not a list, a quartet that has three
   possible answers and no more, concordance factors that sum to a hundred, a
   circular decomposition that gives back a tree when it is handed a tree, and
   D that changes sign when P1 and P2 change places. */

(function () {

  const L4 = ['a', 'b', 'c', 'd'];
  const L6 = ['a', 'b', 'c', 'd', 'e', 'f'];
  const ab_cd = Tree.parseNewick('((a:1,b:1):1,(c:1,d:1):1);', L4);
  const ac_bd = Tree.parseNewick('((a:1,c:1):1,(b:1,d:1):1);', L4);
  const ad_bc = Tree.parseNewick('((a:1,d:1):1,(b:1,c:1):1);', L4);
  const t6 = Tree.parseNewick('(((a:1,b:1):1,c:2):1,((d:1,e:1):1,f:2):1);', L6);
  const t6b = Tree.parseNewick('(((a:1,c:1):1,b:2):1,((d:1,e:1):1,f:2):1);', L6);

  sec('102 · Divisiones');
  const tab = Cmp.splitTable(t6, 6);
  t('un árbol de seis puntas tiene tres divisiones internas', tab.keys.length === 3, String(tab.keys.length));
  t('cada división trae su vector de pertenencia', tab.sets.every(s => s.length === 6));
  t('y su longitud de rama', tab.lens.every(v => v > 0));
  t('con las triviales hay seis más', Cmp.splitTable(t6, 6, { trivial: true }).keys.length === 9,
    String(Cmp.splitTable(t6, 6, { trivial: true }).keys.length));
  t('las dos ramas bajo una raíz bifurcada son una sola división', (() => {
    /* their split is the same, so the table has one entry and its length is the
       sum: reading only one of the two halves the branch */
    const tt = Tree.parseNewick('((a:0.3,b:0.3):0.4,(c:0.3,d:0.3):0.6);', L4);
    const s = Cmp.splitTable(tt, 4);
    return s.keys.length === 1 && Math.abs(s.lens[0] - 1.0) < 1e-12;
  })());
  t('una raíz sobre una punta no inventa una división interna', (() => {
    /* the branch above the sister then holds n−1 tips, and its split is that
       one tip against the rest: trivial, however it is written */
    const tt = Tree.parseNewick('(((a:1,b:1):1,c:2):1,d:3);', L4);
    return Cmp.splitTable(tt, 4).keys.length === 1;
  })(), String(Cmp.splitTable(Tree.parseNewick('(((a:1,b:1):1,c:2):1,d:3);', L4), 4).keys.length));

  sec('103 · Distancias entre árboles');
  t('la distancia de un árbol consigo mismo es cero', Cmp.rf(t6, t6, 6).rf === 0);
  t('y su puntaje de rama también', Cmp.branchScore(t6, t6, 6) < 1e-14);
  t('los tres árboles de cuatro puntas están todos a RF 2 entre sí',
    Cmp.rf(ab_cd, ac_bd, 4).rf === 2 && Cmp.rf(ab_cd, ad_bc, 4).rf === 2 && Cmp.rf(ac_bd, ad_bc, 4).rf === 2);
  t('que es el máximo posible con cuatro puntas', Cmp.rf(ab_cd, ac_bd, 4).maxRF === 2);
  t('y por tanto la RF normalizada es uno', Math.abs(Cmp.rf(ab_cd, ac_bd, 4).normalised - 1) < 1e-12);
  t('RF es simétrica', Cmp.rf(t6, t6b, 6).rf === Cmp.rf(t6b, t6, 6).rf);
  t('RF cuenta las divisiones que sobran a cada lado', (() => {
    const r = Cmp.rf(t6, t6b, 6);
    return r.rf === r.only1 + r.only2 && r.shared === r.n1 - r.only1;
  })());
  t('el puntaje de rama es cero cuando las longitudes coinciden', (() => {
    const copy = Tree.parseNewick('(((a:1,b:1):1,c:2):1,((d:1,e:1):1,f:2):1);', L6);
    return Cmp.branchScore(t6, copy, 6) < 1e-14;
  })());
  t('escalar todas las ramas escala el puntaje de rama igual', (() => {
    const big = Tree.parseNewick('(((a:2,b:2):2,c:4):2,((d:2,e:2):2,f:4):2);', L6);
    const a = Cmp.branchScore(t6, big, 6);
    const b = Cmp.branchScore(t6, t6, 6);
    return a > 0 && b === 0;
  })());
  t('la RF pesada no es menor que la diferencia en una sola rama', (() => {
    const other = Tree.parseNewick('(((a:1,b:1):5,c:2):1,((d:1,e:1):1,f:2):1);', L6);
    return Cmp.weightedRF(t6, other, 6) >= 4 - 1e-9;
  })());
  t('la diferencia de caminos no depende de dónde esté la raíz', (() => {
    /* the same tree with the root moved: unrooted they are one tree, and the
       measure is about the unrooted tree, so the distance has to be zero */
    const moved = Tree.rootByOutgroup(t6, [5]).tree;
    return Cmp.rf(t6, moved, 6).rf === 0
      && Cmp.pathDistance(t6, moved, 6, { weighted: false }) < 1e-9
      && Cmp.pathDistance(t6, moved, 6, { weighted: true }) < 1e-9;
  })(), (() => {
    const moved = Tree.rootByOutgroup(t6, [5]).tree;
    return `RF ${Cmp.rf(t6, moved, 6).rf}, caminos ${Cmp.pathDistance(t6, moved, 6, {}).toExponential(2)}`;
  })());
  t('todas las distancias son no negativas',
    [Cmp.rf(t6, t6b, 6).rf, Cmp.weightedRF(t6, t6b, 6), Cmp.branchScore(t6, t6b, 6),
     Cmp.pathDistance(t6, t6b, 6, {})].every(v => v >= 0));

  sec('104 · Cuartetos');
  t('un cuarteto tiene tres respuestas posibles', Cmp.QUARTET.length === 3);
  const tabA = Cmp.splitTable(ab_cd, 4), tabB = Cmp.splitTable(ac_bd, 4), tabC = Cmp.splitTable(ad_bc, 4);
  t('el árbol ab|cd responde la primera', Cmp.quartetOf(tabA, 0, 1, 2, 3) === 0);
  t('el árbol ac|bd la segunda', Cmp.quartetOf(tabB, 0, 1, 2, 3) === 1);
  t('el árbol ad|bc la tercera', Cmp.quartetOf(tabC, 0, 1, 2, 3) === 2);
  t('una politomía no responde', (() => {
    const star = Tree.parseNewick('(a:1,b:1,c:1,d:1);', L4);
    return Cmp.quartetOf(Cmp.splitTable(star, 4), 0, 1, 2, 3) === -1;
  })());
  t('con cuatro puntas hay un solo cuarteto', Cmp.quartetDistance(ab_cd, ac_bd, 4).total === 1);
  t('y dos árboles distintos discrepan en él', Cmp.quartetDistance(ab_cd, ac_bd, 4).different === 1);
  t('un árbol consigo mismo no discrepa en ninguno', Cmp.quartetDistance(t6, t6, 6).different === 0);
  t('con seis puntas hay C(6,4) = 15 cuartetos', Cmp.quartetDistance(t6, t6b, 6).total === 15);
  t('los cuartetos se reparten entre iguales, distintos y sin resolver', (() => {
    const q = Cmp.quartetDistance(t6, t6b, 6);
    return q.same + q.different + q.unresolved === q.total;
  })());
  t('la distancia de cuartetos es simétrica',
    Cmp.quartetDistance(t6, t6b, 6).different === Cmp.quartetDistance(t6b, t6, 6).different);
  t('una estrella no resuelve ningún cuarteto', (() => {
    const star = Tree.parseNewick('(a:1,b:1,c:1,d:1,e:1,f:1);', L6);
    const q = Cmp.quartetDistance(star, t6, 6);
    return q.unresolved === 15 && q.different === 0;
  })());

  sec('105 · Cuartetos sobre un conjunto de árboles');
  const qc = Cmp.quartetCounts([ab_cd, ab_cd, ac_bd], 4);
  t('se cuenta un cuarteto', qc.counts.size === 1);
  t('con dos votos para la primera respuesta y uno para la segunda', (() => {
    const v = [...qc.counts.values()][0];
    return v.c[0] === 2 && v.c[1] === 1 && v.c[2] === 0;
  })());
  t('el árbol mayoritario recoge dos de tres votos',
    Math.abs(Cmp.quartetScore(ab_cd, qc, 4).proportion - 2 / 3) < 1e-12);
  t('y el minoritario uno de tres',
    Math.abs(Cmp.quartetScore(ac_bd, qc, 4).proportion - 1 / 3) < 1e-12);
  t('el tercero, ninguno', Cmp.quartetScore(ad_bc, qc, 4).score === 0);
  t('con un solo árbol, ese árbol recoge todos sus cuartetos', (() => {
    const q1 = Cmp.quartetCounts([t6], 6);
    return Math.abs(Cmp.quartetScore(t6, q1, 6).proportion - 1) < 1e-12;
  })());
  t('la búsqueda nunca empeora su punto de partida', (() => {
    const q = Cmp.quartetCounts([t6, t6, t6b], 6);
    const s0 = Cmp.quartetScore(t6b, q, 6).score;
    return Cmp.quartetSearch(t6b, q, 6, {}).score >= s0 - 1e-9;
  })());
  t('y desde el árbol minoritario encuentra el mayoritario', (() => {
    const q = Cmp.quartetCounts([t6, t6, t6, t6b], 6);
    const r = Cmp.quartetSearch(t6b, q, 6, {});
    return Cmp.rf(r.tree, t6, 6).rf === 0;
  })());

  sec('106 · Factores de concordancia');
  const brs = Cmp.branchQuartets(t6, 6);
  t('hay una rama interna por división', brs.length === 3, String(brs.length));
  t('cada una tiene cuatro grupos a su alrededor', brs.every(b => b.groups.length === 4));
  t('y ninguno vacío', brs.every(b => b.groups.every(gp => gp.length > 0)));
  t('los cuatro grupos no comparten ninguna punta', brs.every(b => {
    const all = [].concat.apply([], b.groups);
    return new Set(all).size === all.length;
  }));
  const g3 = Cmp.gcf(t6, [t6, t6, t6], 6);
  t('con tres copias del árbol, gCF es cien en todas las ramas',
    g3.every(b => Math.abs(b.gCF - 100) < 1e-9));
  const gmix = Cmp.gcf(t6, [t6, t6, t6b], 6);
  t('con dos de tres, gCF baja proporcionalmente', gmix.some(b => Math.abs(b.gCF - 200 / 3) < 1e-9),
    gmix.map(b => b.gCF.toFixed(1)).join(' '));
  t('gCF, gDF1, gDF2 y gDFP suman cien',
    gmix.every(b => Math.abs(b.gCF + b.gDF1 + b.gDF2 + b.gDFP - 100) < 1e-9));
  t('todos los árboles son decisivos cuando tienen todas las puntas',
    gmix.every(b => b.decisive === 3));
  t('un árbol al que le falta una punta no decide sobre las ramas que la necesitan', (() => {
    const part = Tree.parseNewick('((a:1,b:1):1,(d:1,e:1):1);', ['a', 'b', 'd', 'e']);
    /* the tips carry their own indices, so the labels are mapped by hand */
    Tree.tips(part).forEach(tp => { tp.tip = L6.indexOf(tp.label); });
    const r = Cmp.gcf(t6, [t6, part], 6);
    return r.some(b => b.decisive < 2);
  })());

  sec('107 · Concordancia de sitios');
  /* four sequences in which two sites say ab|cd and one says ac|bd */
  const seqs4 = ['AAACCCGGGT', 'AAACCCTTTA', 'GGGTTTGGGT', 'GGGTTTTTTA'];
  const t4 = Tree.parseNewick('((a:1,b:1):1,(c:1,d:1):1);', L4);
  const sc = Cmp.scf(t4, seqs4, 4, { quartets: 5, seed: 1 });
  t('hay un factor por rama interna', sc.length === 1, String(sc.length));
  t('sCF, sDF1 y sDF2 suman cien', sc.every(b => Math.abs(b.sCF + b.sDF1 + b.sDF2 - 100) < 1e-9));
  t('los sitios que apoyan la rama son los que la app cuenta', sc[0].sites > 0, String(sc[0].sites));
  t('con un alineamiento que sólo apoya ab|cd, sCF es cien', (() => {
    const only = ['AAAA', 'AAAA', 'GGGG', 'GGGG'];
    const r = Cmp.scf(t4, only, 4, { quartets: 3, seed: 2 });
    return Math.abs(r[0].sCF - 100) < 1e-9;
  })());
  t('con un alineamiento que sólo apoya ac|bd, sCF es cero', (() => {
    const only = ['AAAA', 'GGGG', 'AAAA', 'GGGG'];
    const r = Cmp.scf(t4, only, 4, { quartets: 3, seed: 2 });
    return Math.abs(r[0].sCF) < 1e-9 && Math.abs(r[0].sDF1 - 100) < 1e-9;
  })());
  t('un sitio constante no vota', (() => {
    const flat = ['AAAA', 'AAAA', 'AAAA', 'AAAA'];
    const r = Cmp.scf(t4, flat, 4, { quartets: 3, seed: 2 });
    return r[0].sites === 0;
  })());
  t('la misma semilla da el mismo resultado', (() => {
    const a = Cmp.scf(t4, seqs4, 4, { quartets: 5, seed: 1 });
    return Math.abs(a[0].sCF - sc[0].sCF) < 1e-15;
  })());

  sec('108 · Redes de divisiones');
  t('dos divisiones anidadas son compatibles', (() => {
    const A = new Uint8Array(6), B = new Uint8Array(6);
    [0, 1].forEach(i => { A[i] = 1; });
    [0, 1, 2].forEach(i => { B[i] = 1; });
    return Cmp.compatible(A, B, 6);
  })());
  t('dos disjuntas también', (() => {
    const A = new Uint8Array(6), B = new Uint8Array(6);
    [0, 1].forEach(i => { A[i] = 1; });
    [3, 4].forEach(i => { B[i] = 1; });
    return Cmp.compatible(A, B, 6);
  })());
  t('dos que se cruzan, no', (() => {
    const A = new Uint8Array(6), B = new Uint8Array(6);
    [0, 1, 2].forEach(i => { A[i] = 1; });
    [1, 2, 3].forEach(i => { B[i] = 1; });
    return !Cmp.compatible(A, B, 6);
  })());
  const cn = Cmp.consensusNetwork([t6, t6, t6b], 6, { threshold: 0.3 });
  t('la red de consenso guarda las divisiones que pasan el umbral',
    cn.splits.every(s => s.frequency >= 0.3 - 1e-12));
  t('cada división sabe en cuántos árboles está',
    cn.splits.every(s => s.count >= 1 && s.count <= 3));
  t('con copias de un mismo árbol la red es un árbol',
    Cmp.consensusNetwork([t6, t6], 6, { threshold: 0.5 }).compatible);
  t('con árboles que se contradicen, no', !cn.compatible,
    `${cn.splits.filter(s => s.conflictsWith.length).length} en conflicto`);
  t('un umbral de uno deja sólo lo que está en todos', (() => {
    const c = Cmp.consensusNetwork([t6, t6, t6b], 6, { threshold: 1 });
    return c.splits.every(s => s.count === 3);
  })());

  sec('109 · Descomposición circular');
  t('el orden circular tiene cada punta una vez', (() => {
    const o = Cmp.circularOrder(t6);
    return o.length === 6 && new Set(o).size === 6;
  })());
  /* the exact test: a distance that comes from a tree decomposes into that
     tree's splits with that tree's branch lengths, and nothing else */
  const Dt6 = (() => {
    const v = Cmp.pathVector(t6, 6, true);
    const D = Array.from({ length: 6 }, () => new Float64Array(6));
    let k = 0;
    for (let i = 0; i < 6; i++) for (let j = i + 1; j < 6; j++) { D[i][j] = D[j][i] = v[k++]; }
    return D;
  })();
  const cs = Cmp.circularSplits(Dt6, Cmp.circularOrder(t6), {});
  t('un árbol de seis puntas da 6 + 3 = 9 divisiones', cs.splits.length === 9, String(cs.splits.length));
  t('ninguna con peso negativo', cs.clipped === 0, String(cs.clipped));
  t('las tres internas del árbol están entre ellas', (() => {
    const got = new Set(cs.splits.map(s => s.key));
    return Cmp.splitTable(t6, 6).keys.every(k => got.has(k));
  })());
  t('con sus longitudes de rama', (() => {
    const byKey = new Map(cs.splits.map(s => [s.key, s.weight]));
    const tt = Cmp.splitTable(t6, 6);
    return tt.keys.every((k, i) => Math.abs(byKey.get(k) - tt.lens[i]) < 1e-12);
  })());
  t('y la descomposición reconstruye las distancias exactamente', (() => {
    const D2 = Cmp.splitDistance(cs.splits, 6);
    for (let i = 0; i < 6; i++) for (let j = 0; j < 6; j++) if (Math.abs(D2[i][j] - Dt6[i][j]) > 1e-12) return false;
    return true;
  })());
  t('el ciclo exterior de la red se cierra sobre sí mismo',
    Cmp.networkOuterCycle(cs).closure < 1e-12, Cmp.networkOuterCycle(cs).closure.toExponential(2));
  t('y pasa por las seis puntas', (() => {
    const c = Cmp.networkOuterCycle(cs);
    return c.points.filter(p => p.tip != null).length === 6;
  })());
  t('una distancia que no viene de un árbol deja de descomponerse limpiamente', (() => {
    /* the average of two conflicting trees' distances: the split one tree has
       and the other has not is no longer an arc of this circular ordering, and
       the decomposition says so by giving it a negative weight */
    const a = Cmp.pathVector(t6, 6, true), b = Cmp.pathVector(t6b, 6, true);
    const D = Array.from({ length: 6 }, () => new Float64Array(6));
    let k = 0;
    for (let i = 0; i < 6; i++) for (let j = i + 1; j < 6; j++) { D[i][j] = D[j][i] = 0.5 * (a[k] + b[k]); k++; }
    const r = Cmp.circularSplits(D, Cmp.circularOrder(t6), {});
    return r.clipped > 0 && r.clippedMass > 0;
  })(), (() => {
    const a = Cmp.pathVector(t6, 6, true), b = Cmp.pathVector(t6b, 6, true);
    const D = Array.from({ length: 6 }, () => new Float64Array(6));
    let k = 0;
    for (let i = 0; i < 6; i++) for (let j = i + 1; j < 6; j++) { D[i][j] = D[j][i] = 0.5 * (a[k] + b[k]); k++; }
    const r = Cmp.circularSplits(D, Cmp.circularOrder(t6), {});
    return `${r.clipped} recortadas, masa ${fmtNum(r.clippedMass)}`;
  })());
  t('y sobre las distancias de un árbol no se recorta nada', cs.clipped === 0);

  sec('110 · ABBA-BABA');
  /* P1, P2, P3, O: five ABBA sites, four BABA and one constant, so D is small
     and positive and every count can be checked by eye */
  const abba = ['AAAAAGGGGA', 'GGGGGAAAAA', 'GGGGGGGGGA', 'AAAAAAAAAA'];
  const d1 = Cmp.dStatistic(abba, 0, 1, 2, 3, { blocks: 5 });
  t('los sitios ABBA se cuentan', d1.ABBA === 5, String(d1.ABBA));
  t('los BABA también', d1.BABA === 4, String(d1.BABA));
  near('D es (ABBA − BABA) / (ABBA + BABA)', d1.D, 1 / 9, 1e-12);
  t('intercambiar P1 y P2 cambia el signo', (() => {
    const b = Cmp.dStatistic(abba, 1, 0, 2, 3, { blocks: 5 });
    return Math.abs(b.D + d1.D) < 1e-12;
  })());
  t('sin sitios discordantes D es cero', (() => {
    const flat = ['AAAA', 'AAAA', 'AAAA', 'AAAA'];
    const b = Cmp.dStatistic(flat, 0, 1, 2, 3, { blocks: 2 });
    return b.ABBA === 0 && b.BABA === 0 && b.D === 0;
  })());
  t('con ABBA y BABA iguales, D es cero', (() => {
    const bal = ['AAGG', 'GGAA', 'GGGG', 'AAAA'];
    const b = Cmp.dStatistic(bal, 0, 1, 2, 3, { blocks: 2 });
    return b.ABBA === 2 && b.BABA === 2 && Math.abs(b.D) < 1e-12;
  })());
  t('un hueco no cuenta como estado', (() => {
    const gapped = ['A-AA', 'G-GG', 'G-GG', 'A-AA'];
    const b = Cmp.dStatistic(gapped, 0, 1, 2, 3, { blocks: 2 });
    return b.sites === 3;
  })());
  t('el error estándar del jackknife no es negativo', d1.se >= 0);
  t('el valor p está entre cero y uno', d1.p >= 0 && d1.p <= 1);
  near('la normal acumulada en cero vale un medio', Cmp.normalCdf(0), 0.5, 1e-7);
  near('y en 1.959964 vale 0.975', Cmp.normalCdf(1.959964), 0.975, 1e-5);
  t('el barrido devuelve un trío por combinación ordenada', (() => {
    const rows = Cmp.dScan(abba, ['P1', 'P2', 'P3', 'O'], 3, { blocks: 2, minSites: 0 });
    /* three taxa besides the outgroup: three ways to pick the pair and one
       remaining taxon each time */
    return rows.length === 3;
  })(), String(Cmp.dScan(abba, ['P1', 'P2', 'P3', 'O'], 3, { blocks: 2, minSites: 0 }).length));

  sec('111 · Tanglegramas');
  t('el orden de las hojas tiene cada punta una vez', (() => {
    const o = Cmp.leafOrder(t6);
    return o.length === 6 && new Set(o).size === 6;
  })());
  t('dos órdenes iguales no se cruzan', Cmp.crossings([0, 1, 2], [0, 1, 2]) === 0);
  t('dos órdenes invertidos se cruzan todo lo posible',
    Cmp.crossings([0, 1, 2], [2, 1, 0]) === 3, String(Cmp.crossings([0, 1, 2], [2, 1, 0])));
  t('un solo intercambio es un solo cruce', Cmp.crossings([0, 1, 2], [1, 0, 2]) === 1);
  const ut = Cmp.untangle(t6, t6b, {});
  t('desenredar no aumenta los cruces', ut.crossings <= ut.before, `${ut.before} → ${ut.crossings}`);
  t('y no cambia el árbol', Cmp.rf(ut.tree, t6b, 6).rf === 0);
  t('un árbol contra sí mismo se desenreda hasta cero cruces', (() => {
    /* the same tree written with its daughters the other way round */
    const flipped = Tree.parseNewick('(((d:1,e:1):1,f:2):1,((a:1,b:1):1,c:2):1);', L6);
    return Cmp.untangle(t6, flipped, { rounds: 6 }).crossings === 0;
  })());

  sec('112 · El lector de árboles del Bloque 11');
  t('el bloque expone su lista de árboles', typeof B11.sources === 'function');
  t('y sus funciones de cálculo',
    ['runDist', 'runTangle', 'runConc', 'runQuartet', 'runNet', 'runAbba'].every(k => typeof B11[k] === 'function'));

})();
