/* PhylogenyPro — unit tests of the tree studio: groups, labels that keep
   apart, and the pictures of the OTUs.

   Two of these tests exist because the code failed them first. Cutting a tree
   into k groups by its longest branches has two traps, and the engine fell into
   both: the two branches under a bifurcating root are ONE branch of the
   unrooted tree, so cutting them both buys one group and pays for two; and
   cutting both branches of an internal node leaves that node holding no tips at
   all, a group that exists in the graph and not in the data. Asking for k = 3
   returned 2, and asking for 5 returned 4. */

(function () {

  const near = (a, b, tol) => Math.abs(a - b) <= (tol == null ? 1e-9 : tol);

  /* a tree with a bifurcating root, nested cuts, and ties in branch length */
  const NWK = "(((A:0.1,B:0.1):0.4,(C:0.1,D:0.1):0.4):0.5,((E:0.1,F:0.1):0.4,G:0.5):0.5);";
  const LAB = ['A', 'B', 'C', 'D', 'E', 'F', 'G'];
  const mk = () => Tree.parseNewick(NWK, LAB);

  /* ================================================================ */
  sec('1 · k groups by the longest branches');

  for (let k = 1; k <= 7; k++) {
    const g = Groups.byLongestBranches(mk(), k);
    t(`k = ${k} da exactamente ${k} grupo(s)`, g.k === k);
  }
  {
    const g = Groups.byLongestBranches(mk(), 3);
    t('con k = 3 los tamaños suman las 7 puntas', Groups.sizes(g).reduce((a, b) => a + b, 0) === 7);
    t('ningún grupo queda vacío', Groups.sizes(g).every(n => n > 0));
    t('todos los grupos son clados', Groups.monophyly(mk(), g).every(m => m.monophyletic));
  }
  {
    /* no se pueden pedir más grupos que puntas, y hay que decirlo */
    const g = Groups.byLongestBranches(mk(), 12);
    t('pedir más grupos que puntas devuelve una punta por grupo', g.k === 7);
    t('y lo declara con exact = false', g.exact === false);
  }
  {
    /* el corte más largo es la rama de la raíz: debe partir el árbol en las dos
       mitades verdaderas, no dejar a nadie fuera */
    const g = Groups.byLongestBranches(mk(), 2);
    const of = g.of;
    const izq = [0, 1, 2, 3].map(i => of.get(i));
    const der = [4, 5, 6].map(i => of.get(i));
    t('con k = 2 las dos mitades de la raíz quedan separadas',
      izq.every(v => v === izq[0]) && der.every(v => v === der[0]) && izq[0] !== der[0]);
  }
  {
    /* una politomía en la raíz: no hay par que fusionar y el conteo no cambia */
    const p = Tree.parseNewick('(A:0.5,B:0.5,(C:0.1,D:0.1):0.5);', ['A', 'B', 'C', 'D']);
    const g = Groups.byLongestBranches(p, 3);
    t('con la raíz politómica k = 3 sigue dando 3', g.k === 3);
  }

  /* ================================================================ */
  sec('2 · k groups by depth');

  const ULT = Tree.parseNewick('(((A:1,B:1):1,(C:1,D:1):1):1,((E:1,F:1):1,(G:1,H:1):1):1);',
    ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']);
  t('reconoce un árbol ultramétrico', Groups.isUltrametric(ULT, 1e-6) === true);
  /* NWK, pese a las ramas de largos distintos, deja todas las puntas a 1.0:
     sí es ultramétrico, y hay que comprobarlo con uno que no lo sea */
  t('y el de puntas parejas también lo es aunque las ramas difieran',
    Groups.isUltrametric(mk(), 1e-6) === true);
  t('reconoce que un filograma de tasas desiguales no es ultramétrico',
    Groups.isUltrametric(Tree.parseNewick('((A:0.9,B:0.1):0.2,C:0.3);', ['A', 'B', 'C']), 1e-6) === false);
  {
    const g = Groups.byDepth(Tree.clone(ULT), 4);
    t('corta el ultramétrico en 4 grupos', g.k === 4);
    t('de dos puntas cada uno', Groups.sizes(g).every(n => n === 2));
    t('todos son clados', Groups.monophyly(Tree.clone(ULT), g).every(m => m.monophyletic));
    t('y dice a qué altura cortó', typeof g.height === 'number' && g.height > 0);
  }
  {
    const g = Groups.byDepth(Tree.clone(ULT), 2);
    t('con k = 2 deja las dos mitades', g.k === 2 && Groups.sizes(g).every(n => n === 4));
  }

  /* ================================================================ */
  sec('3 · groups by hand and from clades');

  {
    const g = Groups.fromClades(mk(), [[0, 1], [4, 5]]);
    t('lo no asignado cae en un grupo "el resto"', g.k === 3);
    t('y el resto es el grupo 0', g.of.get(2) === 0 && g.of.get(6) === 0);
    t('los clados dados conservan su orden', g.of.get(0) === 1 && g.of.get(4) === 2);
  }
  {
    const g0 = Groups.empty(mk());
    const g1 = Groups.assign(g0, 3, 2);
    t('asignar a mano no toca a los demás', g1.of.get(0) === 0 && g1.of.get(3) === 2);
    t('y el original queda intacto', g0.of.get(3) === 0);
  }
  {
    /* un grupo parafilético: el aviso del panel depende de esto */
    const g = Groups.fromClades(mk(), [[0, 6]]);
    const m = Groups.monophyly(mk(), g);
    t('detecta que {A, G} no es un clado', m.find(x => x.group === 1).monophyletic === false);
  }
  t('nombra los grupos A, B, C…', Groups.nameOf({ k: 3 }, 0) === 'A' && Groups.nameOf({ k: 3 }, 2) === 'C');
  t('y respeta los nombres que le pongan', Groups.nameOf({ k: 2, names: ['Norte', 'Sur'] }, 1) === 'Sur');
  t('la paleta da tantos colores como grupos', Groups.palette(14).length === 14);
  t('y no repite ninguno', new Set(Groups.palette(14)).size === 14);

  /* ================================================================ */
  sec('4 · labels that keep apart');

  {
    const want = [10, 12, 13, 40, 41, 80];
    const r = TreeView.spreadLabels(want, 10, 0, 200);
    const sorted = r.pos.slice().sort((a, b) => a - b);
    let minGap = Infinity;
    for (let i = 1; i < sorted.length; i++) minGap = Math.min(minGap, sorted[i] - sorted[i - 1]);
    t('ningún par queda más cerca que el hueco pedido', minGap >= 10 - 1e-9);
    t('todo cae dentro de los límites', r.pos.every(p => p >= -1e-9 && p <= 200 + 1e-9));
    t('no se declara apretado cuando sí cabía', r.crowded === false);
    t('el orden relativo se conserva',
      r.pos.map((p, i) => [p, want[i]]).sort((a, b) => a[1] - b[1]).every((v, i, arr) => i === 0 || v[0] >= arr[i - 1][0]));
  }
  {
    /* sin espacio: se reparte parejo y se avisa, en vez de amontonar */
    const r = TreeView.spreadLabels([1, 2, 3, 4, 5], 30, 0, 50);
    t('cuando no cabe, lo dice', r.crowded === true);
    const sorted = r.pos.slice().sort((a, b) => a - b);
    t('y reparte parejo', near(sorted[1] - sorted[0], sorted[2] - sorted[1], 1e-9));
  }
  {
    const r = TreeView.spreadLabels([], 10, 0, 100);
    t('una lista vacía no rompe nada', r.pos.length === 0);
  }
  {
    /* el que se sale por abajo se jala hacia arriba */
    const r = TreeView.spreadLabels([95, 96, 97], 10, 0, 100);
    t('lo que se pasaba del borde se recoge', Math.max.apply(null, r.pos) <= 100 + 1e-9);
  }

  /* ================================================================ */
  sec('5 · labels that cannot move: dropping the ones that collide');

  {
    const boxes = [
      { x: 0, y: 0, w: 20, h: 10, rank: 0 },
      { x: 5, y: 2, w: 20, h: 10, rank: 1 },     // encima del primero
      { x: 60, y: 0, w: 20, h: 10, rank: 2 },    // lejos
    ];
    const keep = TreeView.dropColliding(boxes);
    t('conserva el de mayor prioridad', keep[0] === true);
    t('tira el que se le encima', keep[1] === false);
    t('y conserva el que no estorba', keep[2] === true);
  }
  {
    const boxes = [];
    for (let i = 0; i < 40; i++) boxes.push({ x: i, y: 0, w: 10, h: 10, rank: i });
    const keep = TreeView.dropColliding(boxes);
    const kept = boxes.filter((b, i) => keep[i]);
    let choca = false;
    for (let i = 0; i < kept.length; i++) for (let j = i + 1; j < kept.length; j++) {
      const a = kept[i], b = kept[j];
      if (!(b.x + b.w < a.x || a.x + a.w < b.x || b.y + b.h < a.y || a.y + a.h < b.y)) choca = true;
    }
    t('ningún par de los que quedan se encima', choca === false);
    t('y quedan los que caben', kept.length === 4);
  }

  /* ================================================================ */
  sec('6 · the figure with groups and pictures');

  {
    const g = Groups.byLongestBranches(mk(), 3);
    const svg = TreeView.render(mk(), {
      labels: LAB, groups: g, groupColours: Groups.palette(3),
      groupStrip: true, legend: true, alignTips: true, width: 600,
    });
    t('el SVG sale bien formado', /^<svg[\s\S]*<\/svg>$/.test(svg.trim()));
    t('lleva la tira de colores', (svg.match(/<rect/g) || []).length >= 3);
    t('y anota cuántos grupos dibujó', TreeView.lastInfo.k === 3);
    t('los colores de la paleta están en el dibujo',
      Groups.palette(3).every(c => svg.indexOf(c) >= 0));
  }
  {
    /* sin imágenes confirmadas no debe dibujarse ninguna: una foto sin
       confirmar es una afirmación que nadie hizo */
    const svg = TreeView.render(mk(), { labels: LAB, images: 'tips', imgSize: 40, width: 600 });
    t('sin confirmar, no se dibuja ninguna imagen', (svg.match(/<image/g) || []).length === 0);
    t('y el informe de la figura lo dice', TreeView.lastInfo.images === 0);
  }
  {
    const svg = TreeView.render(mk(), { labels: LAB, images: 'none', width: 600 });
    t('sin imágenes no se define ningún recorte', svg.indexOf('clipPath') < 0);
  }
  t('el recorte circular es un círculo', /circle/.test(TreeView.imageDefs('circle', 50, 'x')));
  t('el de esquinas suaves es un rectángulo con radio', /rect[^>]*rx="9\.0"/.test(TreeView.imageDefs('round', 50, 'x')));
  t('«sin recorte» no define nada', TreeView.imageDefs('none', 50, 'x') === '');

  /* ================================================================ */
  sec('7 · pictures: names, and the promise not to draw an unconfirmed one');

  t('empareja el archivo con su taxón',
    OTUImg.matchName('Bursera_simaruba.jpg', ['Bursera simaruba', 'Bursera aptera']) === 'Bursera simaruba');
  t('no se inventa un emparejamiento ambiguo',
    OTUImg.matchName('Bursera.jpg', ['Bursera simaruba', 'Bursera aptera']) === null);
  t('ignora acentos y signos',
    OTUImg.matchName('Piñón-mexicano.png', ['Piñon mexicano', 'Otra cosa']) === 'Piñon mexicano');
  t('un nombre sin registro no tiene imagen', OTUImg.url('no existe este taxón') === null);
  t('los ajustes de fábrica no cambian nada',
    OTUImg.ADJ0.zoom === 1 && OTUImg.ADJ0.dx === 0 && OTUImg.ADJ0.removeBg === false);

  /* ================================================================ */
  sec('8 · circular and unrooted: the same options, the same promises');

  {
    const g = Groups.byLongestBranches(mk(), 3);
    const base = {
      labels: LAB, groups: g, groupColours: Groups.palette(3),
      groupStrip: true, legend: true, alignTips: true, width: 600,
    };
    ['circular', 'unrooted'].forEach(lay => {
      const svg = TreeView.render(mk(), Object.assign({}, base, { layout: lay }));
      t(`${lay}: el SVG sale bien formado`, /^<svg[\s\S]*<\/svg>$/.test(svg.trim()));
      t(`${lay}: pinta los tres colores de grupo`,
        Groups.palette(3).every(c => svg.indexOf(c) >= 0));
      t(`${lay}: dibuja la leyenda`, (svg.match(/<rect/g) || []).length >= 3);
      t(`${lay}: informa cuántos grupos`, TreeView.lastInfo.k === 3);
    });
  }
  {
    /* Sin imágenes confirmadas no se dibuja ninguna, en ningún trazado. */
    ['rect', 'circular', 'unrooted'].forEach(lay => {
      const svg = TreeView.render(mk(), { labels: LAB, images: 'tips', imgSize: 40, width: 600, layout: lay });
      t(`${lay}: sin confirmar no dibuja imágenes`, (svg.match(/<image/g) || []).length === 0);
    });
  }
  {
    /* Pedir imágenes debe AGRANDAR el lienzo, no encoger el árbol: reservar el
       espacio quitándoselo al radio dejaba el árbol en siete píxeles en cuanto
       los nombres eran largos. */
    const opts = { labels: LAB, width: 420, tipFont: 12, layout: 'circular' };
    const sinImg = TreeView.render(mk(), opts);
    TreeView.render(mk(), Object.assign({}, opts, { images: 'tips', imgSize: 48 }));
    const conImg = TreeView.lastInfo;
    const anchoSin = +sinImg.match(/viewBox="0 0 ([\d.]+)/)[1];
    t('circular: el lienzo crece cuando se piden imágenes', conImg.width > anchoSin);
    t('circular: y crece justo lo que ocupan', conImg.width - anchoSin >= 96);
  }
  {
    const opts = { labels: LAB, width: 420, layout: 'unrooted' };
    const sinImg = TreeView.render(mk(), opts);
    const altoSin = +sinImg.match(/viewBox="0 0 [\d.]+ ([\d.]+)/)[1];
    TreeView.render(mk(), Object.assign({}, opts, { images: 'tips', imgSize: 48 }));
    t('no enraizado: el lienzo también crece', TreeView.lastInfo.height > altoSin);
  }
  {
    /* el no enraizado no puede repartir por un eje: oculta y lo cuenta */
    const muchos = [];
    for (let i = 0; i < 40; i++) muchos.push('Especie con nombre largo ' + i);
    const nk = '(' + muchos.map((m, i) => `'${m}':0.1`).join(',') + ');';
    const t40 = Tree.parseNewick(nk, muchos);
    TreeView.render(t40, { labels: muchos, layout: 'unrooted', width: 300, tipFont: 12 });
    t('no enraizado: oculta nombres encimados y los cuenta', TreeView.lastInfo.hiddenLabels > 0);
    TreeView.render(t40, { labels: muchos, layout: 'unrooted', width: 300, tipFont: 12, labelDeclutter: false });
    t('y se puede apagar si el usuario los quiere todos', TreeView.lastInfo.hiddenLabels === 0);
  }

  /* ================================================================ */
  sec('9 · the other figures: where labels really did collide');

  /* Cuántos pares de etiquetas se encimarían, contado sobre el SVG. */
  function choques(svg, re, dx, dy) {
    const a = [...svg.matchAll(re)].map(m => [+m[1], +m[2]]);
    let c = 0;
    for (let i = 0; i < a.length; i++) for (let j = i + 1; j < a.length; j++) {
      if (Math.abs(a[i][1] - a[j][1]) < dy && Math.abs(a[i][0] - a[j][0]) < dx) c++;
    }
    return { n: a.length, c };
  }

  {
    /* El cronograma denso: 64 puntas y las edades encima unas de otras. */
    const names = [];
    for (let i = 0; i < 64; i++) names.push('sp' + i);
    const bal = (ns) => (ns.length === 1 ? `${ns[0]}:0.05`
      : `(${bal(ns.slice(0, Math.ceil(ns.length / 2)))},${bal(ns.slice(Math.ceil(ns.length / 2)))}):0.05`);
    const tb = Tree.parseNewick(bal(names) + ';', names);
    const re = /<text x="([\d.]+)" y="([\d.]+)"[^>]*text-anchor="end"/g;
    const o = { labels: names, showAges: true, width: 520, rowHeight: 7 };

    const sin = choques(Plots7.chronogram(tb, Object.assign({ ageDeclutter: false }, o)), re, 26, 8);
    t('sin limpiar, el cronograma denso encima decenas de edades', sin.c > 20);
    const con = choques(Plots7.chronogram(tb, o), re, 26, 8);
    t('con la limpieza no queda ni un par encimado', con.c === 0);
    t('y dice cuántas edades ocultó', Plots7.lastInfo.hiddenAges === sin.n - con.n);
    t('conservando más de la mitad', con.n > sin.n * 0.4);
  }
  {
    /* El diagrama gCF–sCF con etiquetas: es opcional, y cuando se enciende
       sobre ramas amontonadas era una mancha. */
    const rows = [];
    for (let i = 0; i < 60; i++) {
      rows.push({ x: 40 + (i % 12) * 0.7, y: 42 + Math.floor(i / 12) * 0.9, label: 'rama ' + i });
    }
    const re = /<text x="([\d.]+)" y="([\d.]+)"[^>]*font-size="8"/g;
    const sin = choques(Plots11.concordanceScatter(rows, { labelAll: true, labelDeclutter: false, width: 640, height: 420 }), re, 40, 9);
    t('sin limpiar, 60 etiquetas en un rincón son una mancha', sin.c > 500);
    const con = choques(Plots11.concordanceScatter(rows, { labelAll: true, width: 640, height: 420 }), re, 40, 9);
    t('con la limpieza no se encima ninguna', con.c === 0);
    t('y lo dice', Plots11.lastInfo.hiddenLabels === 60 - con.n);
  }
  {
    /* El traitgram del Bloque 9, el peor caso del principio. */
    const pts = [], edges = [], labels = [];
    for (let i = 0; i < 25; i++) {
      pts.push({ node: i, tip: true, tipRow: i, depth: 1, value: 1 + i * 0.004 });
      labels.push('Bursera especie ' + (i + 1));
      edges.push({ from: 99, to: i });
    }
    pts.push({ node: 99, tip: false, depth: 0, value: 1.05 });
    const svg = Plots9.traitgram(pts, edges, { labels, width: 740, height: 400 });
    const re = /<text x="([\d.]+)" y="([\d.]+)"[^>]*font-size="8.5"/g;
    const r = choques(svg, re, 90, 9);
    t('el traitgram dibuja las 25 puntas', r.n === 25);
    t('sin encimar ninguna', r.c === 0);
    t('y sin tener que ocultar nada', Plots9.lastInfo.hiddenLabels === 0);
  }
  {
    /* Las marcas de los ejes salen de un generador que apunta a 5–7 por eje:
       no chocan, y conviene dejarlo comprobado para que siga siendo cierto. */
    const pts = [], edges = [];
    for (let i = 0; i < 12; i++) { pts.push({ node: i, tip: true, tipRow: i, depth: 1, value: i * 0.001 }); edges.push({ from: 99, to: i }); }
    pts.push({ node: 99, tip: false, depth: 0, value: 0 });
    const svg = Plots9.traitgram(pts, edges, { width: 740, height: 400 });
    const re = /<text x="([\d.]+)" y="([\d.]+)"[^>]*font-size="9"/g;
    const r = choques(svg, re, 34, 10);
    t('las marcas de los ejes no se enciman entre sí', r.c === 0);
  }

})();
