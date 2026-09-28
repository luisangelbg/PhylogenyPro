/* PhylogenyPro — unit tests of the figure studio.

   The studio restyles every figure by writing one CSS rule instead of touching
   any plotting module, and that trick has exactly one place it can fail: the
   export. An SVG loaded into an <img> carries no stylesheet, so if the chosen
   colours and the chosen type are not written *into* the file before it leaves
   the page, the figure on screen and the figure in the PNG are different
   pictures. Most of what follows tests that they are the same one. */

(function () {

  const LAB = ['A', 'B', 'C', 'D'];
  const tree = () => Tree.parseNewick('((A:0.1,B:0.1):0.2,(C:0.1,D:0.1):0.2);', LAB);
  const before = FigStyle.get();

  sec('1 · the catalogue of looks');

  t('hay varias paletas', FigStyle.PALETTES.length >= 6);
  t('cada paleta trae diez colores', FigStyle.PALETTES.every(p => p.c.length === 10));
  t('todos los colores son hexadecimales válidos',
    FigStyle.PALETTES.every(p => p.c.every(c => /^#[0-9a-f]{6}$/i.test(c))));
  t('ninguna paleta repite un color',
    FigStyle.PALETTES.every(p => new Set(p.c.map(c => c.toLowerCase())).size === 10));
  t('cada paleta se explica en los dos idiomas',
    FigStyle.PALETTES.every(p => p.name.length === 2 && p.note && p.note.length === 2 && p.note[0] && p.note[1]));
  t('hay una declarada segura al daltonismo', FigStyle.PALETTES.some(p => p.id === 'okabe'));
  t('y una en escala de grises, para revistas que cobran el color',
    FigStyle.PALETTES.some(p => p.id === 'grey'));
  t('hay varias familias de letra', FigStyle.FONTS.length >= 5);
  t('y ninguna se baja de la red',
    FigStyle.FONTS.every(f => !/https?:|url\(|@import/.test(f.css)));
  t('cada familia declara una alternativa genérica',
    FigStyle.FONTS.every(f => /(sans-serif|serif|monospace)\s*$/.test(f.css)));

  sec('2 · the state, and the rule it writes');

  {
    FigStyle.reset();
    const r0 = FigStyle.cssText();
    t('la regla se aplica a svg y sólo a svg', /^svg\{/.test(r0));
    t('define los diez colores', [1, 2, 3, 4, 5, 6, 7, 8, 9, 10].every(i => r0.indexOf('--c' + i + ':') > 0));
    t('y la variable de la letra', r0.indexOf('--fig-font:') > 0);
  }
  {
    FigStyle.set({ palette: 'okabe' });
    t('cambiar de paleta cambia los colores', FigStyle.colours()[0] === '#0072B2');
    FigStyle.set({ custom: { c1: '#123456' } });
    t('un color a mano gana a la paleta', FigStyle.colours()[0] === '#123456');
    t('y no toca a los demás', FigStyle.colours()[1] === '#E69F00');
    FigStyle.set({ palette: 'phylo' });
    t('el color a mano sobrevive al cambio de paleta', FigStyle.colours()[0] === '#123456');
    FigStyle.reset();
    t('volver al original lo borra', FigStyle.colours()[0] === '#2a4a94');
  }
  {
    FigStyle.set({ fontScale: 1.4, lineScale: 2 });
    const r = FigStyle.cssText();
    t('el tamaño de letra entra en la regla', /font-size:140\.0%/.test(r));
    t('y el grosor del trazo también', /stroke-width:2px/.test(r));
    FigStyle.reset();
    t('con los valores de fábrica no se escribe ninguna de las dos',
      FigStyle.cssText().indexOf('stroke-width') < 0 && FigStyle.cssText().indexOf('font-size:1') < 0);
  }
  t('get() devuelve una copia, no el estado vivo', (() => {
    const g = FigStyle.get(); g.palette = 'destrozado';
    return FigStyle.get().palette !== 'destrozado';
  })());

  sec('3 · what the export carries');

  {
    FigStyle.set({ palette: 'okabe', font: 'serif' });
    /* un árbol liso dibuja sus ramas con el color del texto y no toca la
       paleta: para ver la paleta en el archivo hace falta algo que la use */
    const g = Groups.byLongestBranches(tree(), 2);
    const svg = TreeView.render(tree(), {
      labels: LAB, groups: g, groupColours: FigStyle.colours().slice(0, 2), width: 400,
    });
    t('la figura sin resolver sigue usando variables', svg.indexOf('var(--') > 0);
    const res = Export.inlineVars(svg, Export.readVars());
    t('al resolverla no queda ninguna variable', res.indexOf('var(--') < 0);
    t('lleva el color de la paleta elegida', res.indexOf('#0072B2') >= 0);
    t('y la letra elegida', res.indexOf('Georgia') >= 0);
    t('y ya no la del sistema', res.indexOf('font-family:system-ui,sans-serif') < 0);
  }
  {
    /* .art-ink lo usaba el cronograma y no estaba definida en la hoja de
       estilo: sus nombres salían negros e invisibles en el tema oscuro */
    const res = Export.inlineVars('<svg><text class="art-ink">x</text></svg>', Export.readVars());
    t('el estilo exportado define .art-ink', res.indexOf('.art-ink{') >= 0);
  }
  {
    FigStyle.set({ background: 'white' });
    const res = Export.inlineVars(TreeView.render(tree(), { labels: LAB, width: 300 }), Export.readVars());
    t('un fondo elegido se escribe dentro del archivo', res.indexOf('data-figbg') >= 0);
    t('y va primero, debajo de todo', res.indexOf('data-figbg') < res.indexOf('<line'));
    FigStyle.set({ background: 'none' });
    const res2 = Export.inlineVars(TreeView.render(tree(), { labels: LAB, width: 300 }), Export.readVars());
    t('sin fondo no se añade nada', res2.indexOf('data-figbg') < 0);
  }
  {
    FigStyle.reset();
    const v = Export.readVars();
    t('readVars sigue trayendo el tema completo', v.text && v.border && v.c1);
    t('y los colores del estudio ganan',
      (FigStyle.set({ custom: { c5: '#abcdef' } }), Export.readVars().c5 === '#abcdef'));
    FigStyle.reset();
  }

  sec('4 · the groups of the tree follow the palette');

  {
    FigStyle.set({ palette: 'contrast' });
    const g = Groups.byLongestBranches(tree(), 2);
    const cols = FigStyle.colours().slice(0, 2);
    const svg = TreeView.render(tree(), { labels: LAB, groups: g, groupColours: cols, width: 400 });
    t('el árbol se pinta con los colores del estudio', cols.every(c => svg.indexOf(c) >= 0));
    FigStyle.reset();
  }

  /* dejar el estudio como estaba, para no torcer las demás pruebas */
  FigStyle.reset();
  FigStyle.set(before);
  FigStyle.reset();

})();
