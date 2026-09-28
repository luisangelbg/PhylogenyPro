/* PhylogenyPro — unit tests of Block 12 (figures, formats and the report).

   This block hands the work over, so what has to be true of it is not a number
   but a property: that another program can read what it writes. The reference
   is therefore R itself — validation/block12/validate_block12.R opens the
   package the app produced and reads every file in it with ape, xml2 and
   unzip, and reports 38 agreements and no mismatches.

   Here the checks are the ones a unit test can make: that every writer round
   trips through its own reader, that the arithmetic of dots per inch is the
   arithmetic of pixels, that the CRC32 of "123456789" is the value the ZIP
   specification names, and — the one that matters most — that the report says
   only what ran. */

(function () {

  const L4 = ['a', 'b', 'c', 'd'];
  const L6 = ['a', 'b', 'c', 'd', 'e', 'f'];
  const t4 = Tree.parseNewick('((a:0.1,b:0.2):0.3,(c:0.15,d:0.25):0.35);', L4);
  const t6 = Tree.parseNewick('(((a:1,b:1):1,c:2):1,((d:1,e:1):1,f:2):1);', L6);
  /* a tree carrying what only phyloXML can hold: support, age, and both */
  const rich = Tree.clone(t6);
  (function mark(nd, d) {
    if (nd.tip != null) return;
    nd.support = 70 + d * 7;
    nd.age = 3 - d * 0.5;
    nd.children.forEach(c => mark(c.node, d + 1));
  })(rich, 0);

  sec('113 · Newick');
  const nk = Export.newick(t4, { labels: L4, decimals: 8 });
  t('termina en punto y coma', /;$/.test(nk.trim()));
  t('lleva los cuatro nombres', L4.every(l => nk.indexOf(l) >= 0));
  t('lo que escribe lo vuelve a leer igual', Cmp.rf(Tree.parseNewick(nk, L4), t4, 4).rf === 0);
  t('con las mismas longitudes de rama', Cmp.branchScore(Tree.parseNewick(nk, L4), t4, 4) < 1e-8,
    Cmp.branchScore(Tree.parseNewick(nk, L4), t4, 4).toExponential(2));
  t('sin longitudes no escribe ningún dos puntos',
    Export.newick(t4, { labels: L4, lengths: false }).indexOf(':') < 0);
  t('los decimales se respetan', (() => {
    const two = Export.newick(t4, { labels: L4, decimals: 2 });
    return /0\.10[,)]/.test(two) && two.indexOf('0.100000') < 0;
  })());
  t('un nombre con espacios se entrecomilla',
    Export.newick(t4, { labels: ['a b', 'c', 'd', 'e'] }).indexOf("'a b'") >= 0);
  t('y uno con coma también',
    Export.newick(t4, { labels: ['x,y', 'c', 'd', 'e'] }).indexOf("'x,y'") >= 0);

  sec('114 · NEXUS');
  const nex = Export.nexusTrees([{ name: 'uno', tree: t6 }, { name: 'dos', tree: rich }],
    { labels: L6, taxa: L6, figtree: true });
  t('empieza por #NEXUS', /^#NEXUS/.test(nex));
  t('declara los seis taxones', /NTAX=6/.test(nex));
  t('los nombra en TAXLABELS', L6.every(l => new RegExp('^\\s*' + l + '\\s*$', 'm').test(nex)));
  t('trae el bloque TRANSLATE', /TRANSLATE/.test(nex));
  t('y los árboles llaman a los nombres por número', /TREE uno = \[&R\] \(\(\(1/.test(nex));
  t('trae los dos árboles', (nex.match(/^\s*TREE /gm) || []).length === 2,
    String((nex.match(/^\s*TREE /gm) || []).length));
  t('cierra todos sus bloques', (nex.match(/END;/g) || []).length >= 2);
  t('trae el bloque de FigTree cuando se pide', /begin figtree;/.test(nex));
  t('y no lo trae cuando no', !/begin figtree;/.test(
    Export.nexusTrees([{ name: 'uno', tree: t6 }], { labels: L6, taxa: L6 })));
  t('sin TRANSLATE escribe los nombres enteros',
    /TREE uno = \[&R\] \(\(\(a/.test(Export.nexusTrees([{ name: 'uno', tree: t6 }],
      { labels: L6, taxa: L6, translate: false })));
  t('quoteNexus entrecomilla lo que hay que entrecomillar',
    Export.quoteNexus('a b') === "'a b'" && Export.quoteNexus('ab') === 'ab'
    && Export.quoteNexus("d'or") === "'d''or'");

  sec('115 · phyloXML');
  const xml = Export.phyloXML([{ name: 'con todo', tree: rich, supportType: 'bootstrap', unit: 'Ma' }],
    { labels: L6 });
  t('es XML bien formado', (() => {
    const doc = new DOMParser().parseFromString(xml, 'application/xml');
    return doc.getElementsByTagName('parsererror').length === 0;
  })());
  t('declara el espacio de nombres', /xmlns="http:\/\/www\.phyloxml\.org"/.test(xml));
  t('lleva los apoyos con su tipo', /<confidence type="bootstrap">/.test(xml));
  t('y las edades con su unidad', /<date unit="Ma">/.test(xml));
  t('nombra las seis puntas', L6.every(l => xml.indexOf(`<name>${l}</name>`) >= 0));
  const back = Export.readPhyloXML(xml, L6);
  t('se vuelve a leer', !!back && back.length === 1);
  t('con el mismo árbol', back && Cmp.rf(back[0].tree, rich, 6).rf === 0);
  t('con las mismas longitudes', back && Cmp.branchScore(back[0].tree, rich, 6) < 1e-8,
    back ? Cmp.branchScore(back[0].tree, rich, 6).toExponential(2) : '');
  t('y con los apoyos intactos', (() => {
    if (!back) return false;
    const gather = nd => {
      const o = [];
      (function w(x) { if (x.tip != null) return; if (x.support != null) o.push(+x.support); x.children.forEach(c => w(c.node)); })(nd);
      return o;
    };
    const a = gather(rich), b = gather(back[0].tree);
    return a.length === b.length && a.every((v, i) => Math.abs(v - b[i]) < 1e-9);
  })());
  t('un & en un nombre no rompe el XML', (() => {
    const s = Export.phyloXML([{ tree: t4 }], { labels: ['A & B', '<c>', 'd', 'e'] });
    const doc = new DOMParser().parseFromString(s, 'application/xml');
    return doc.getElementsByTagName('parsererror').length === 0 && s.indexOf('A &amp; B') >= 0;
  })());
  t('un árbol sin raíz se declara como tal',
    /rooted="false"/.test(Export.phyloXML([{ tree: t4, rooted: false }], { labels: L4 })));
  t('un XML roto se rechaza en vez de dar un árbol inventado',
    Export.readPhyloXML('<phyloxml><phylogeny>', L4) === null);

  sec('116 · Los píxeles de una resolución');
  const svg = '<svg viewBox="0 0 800 600"></svg>';
  const s96 = Export.pngSize(svg, { dpi: 96 });
  t('a 96 ppp el PNG mide lo que el dibujo', s96.width === 800 && s96.height === 600);
  t('a 300 ppp mide 3.125 veces más', (() => {
    const s = Export.pngSize(svg, { dpi: 300 });
    return s.width === 2500 && s.height === 1875;
  })());
  t('y a 900 ppp, 9.375 veces', (() => {
    const s = Export.pngSize(svg, { dpi: 900 });
    return s.width === 7500 && s.height === 5625;
  })());
  t('el número de píxeles es el producto', (() => {
    const s = Export.pngSize(svg, { dpi: 600 });
    return s.pixels === s.width * s.height;
  })());
  t('una figura que no cabe en un lienzo se marca', (() => {
    const huge = '<svg viewBox="0 0 4000 4000"></svg>';
    return Export.pngSize(huge, { dpi: 900 }).tooBig;
  })());
  t('y una normal, no', !Export.pngSize(svg, { dpi: 900 }).tooBig);
  t('sin viewBox se supone un tamaño y no se falla', (() => {
    const s = Export.pngSize('<svg></svg>', { dpi: 96 });
    return s.width > 0 && s.height > 0;
  })());
  t('las variables de color se resuelven', (() => {
    const r = Export.inlineVars('<svg viewBox="0 0 9 9"><line stroke="var(--accent)"/></svg>', { accent: '#b76a17' });
    return r.indexOf('var(--') < 0 && r.indexOf('#b76a17') >= 0;
  })());
  t('una variable que no existe cae en un color, no en nada', (() => {
    const r = Export.inlineVars('<svg viewBox="0 0 9 9"><line stroke="var(--noexiste)"/></svg>', {});
    return r.indexOf('var(--') < 0 && /stroke="#/.test(r);
  })());
  t('y se añade el estilo de los textos, que un <img> no hereda',
    Export.inlineVars('<svg viewBox="0 0 9 9"></svg>', { text: '#111' }).indexOf('.art-txt') >= 0);

  sec('117 · CRC32 y el ZIP');
  t('CRC32 de "123456789" es 0xCBF43926',
    Export.crc32(new TextEncoder().encode('123456789')) === 0xCBF43926,
    '0x' + Export.crc32(new TextEncoder().encode('123456789')).toString(16).toUpperCase());
  t('CRC32 de la cadena vacía es cero', Export.crc32(new Uint8Array(0)) === 0);
  t('CRC32 de "a" es 0xE8B7BE43',
    Export.crc32(new TextEncoder().encode('a')) === 0xE8B7BE43,
    '0x' + Export.crc32(new TextEncoder().encode('a')).toString(16).toUpperCase());
  t('nunca sale negativo', Export.crc32(new TextEncoder().encode('ñ á é í ó ú')) >= 0);

  sec('118 · El informe dice lo que se corrió');
  const soloDatos = {
    data: { taxa: L6, parts: [{ name: 'gen1', taxa: L6, seqs: L6.map(() => 'ACGT'.repeat(50)) }] },
  };
  const s1 = Report.methods(soloDatos, {});
  t('con sólo datos hay dos secciones', s1.length === 2, s1.map(s => s.id).join(','));
  t('los datos y el programa', s1.map(s => s.id).join(',') === 'data,software');
  const conPars = Object.assign({}, soloDatos, {
    quick: { parsimony: { steps: 412, nTrees: 3, indices: { ci: 0.61, ri: 0.78 } } },
  });
  const s2 = Report.methods(conPars, {});
  t('añadir parsimonia añade una sección', s2.length === 3);
  t('con los números de la corrida, no de una plantilla', (() => {
    const txt = Report.methodsText(s2, 'es');
    return txt.indexOf('412') >= 0 && txt.indexOf('0.61') >= 0 && txt.indexOf('3 árboles') >= 0;
  })());
  t('no inventa una sección de verosimilitud', s2.every(s => s.id !== 'ml'));
  t('ni de reloj molecular', s2.every(s => s.id !== 'dated'));
  const r2 = Report.references(s2);
  t('cita a Fitch y a Sankoff, que es lo que usó',
    r2.some(r => r.id === 'fitch1971') && r2.some(r => r.id === 'sankoff1975'));
  t('y no cita a Sanderson, que no usó', !r2.some(r => r.id === 'sanderson2002'));
  t('ninguna cita queda sin texto', r2.every(r => r.text && r.text.length > 40));
  t('la bibliografía va en orden alfabético', (() => {
    const texts = r2.map(r => r.text);
    return texts.slice().sort((a, b) => a.localeCompare(b, 'es')).join('|') === texts.join('|');
  })());
  t('un apoyo que no se calculó no se menciona', (() => {
    const soloBoot = Object.assign({}, conPars, { ml: { lnL: -100, support: { boot: { reps: 250 } } } });
    const txt = Report.methodsText(Report.methods(soloBoot, {}), 'es');
    return txt.indexOf('250 réplicas') >= 0 && txt.indexOf('ultrarrápido') < 0;
  })());
  t('y el número de réplicas viene de la corrida', (() => {
    const dos = Object.assign({}, conPars, { ml: { lnL: -100, support: { uf: { reps: 1234 } } } });
    return Report.methodsText(Report.methods(dos, {}), 'es').indexOf('1234') >= 0;
  })());
  t('el aviso sobre +J aparece cuando gana un modelo +J', (() => {
    const st = Object.assign({}, conPars, {
      biogeo: { areas: ['K', 'O'], maxAreas: 2, nStates: 4, includeNull: true, best: 'DEC+J', models: [1, 2] },
    });
    const s = Report.methods(st, {}).find(x => x.id === 'biogeo');
    return s.es.join(' ').indexOf('Ree y Sanmartín') >= 0 && s.en.join(' ').indexOf('Ree & Sanmartín') >= 0;
  })());
  t('y no aparece cuando gana uno sin +J', (() => {
    const st = Object.assign({}, conPars, {
      biogeo: { areas: ['K', 'O'], maxAreas: 2, nStates: 4, includeNull: true, best: 'DEC', models: [1, 2] },
    });
    const s = Report.methods(st, {}).find(x => x.id === 'biogeo');
    return s.es.join(' ').indexOf('Ree y Sanmartín') < 0;
  })());

  sec('119 · Los dos idiomas');
  const bi = Report.methods({
    data: { taxa: L6, parts: [{ name: 'g', taxa: L6, seqs: L6.map(() => 'ACGT') }] },
    ml: { lnL: -100, support: { boot: { reps: 100 }, uf: { reps: 1000 }, alrt: {} } },
    dated: { method: 'pl', rootMethod: 'outgroup' },
  }, {});
  t('cada sección trae su título en los dos idiomas',
    bi.every(s => s.titleEs && s.titleEn && s.titleEs !== s.titleEn));
  t('y su texto en los dos', bi.every(s => s.es.length && s.en.length));
  t('el párrafo en inglés no lleva palabras en español', (() => {
    const en = bi.map(s => s.en.join(' ')).join(' ');
    return en.indexOf('arranque') < 0 && en.indexOf('réplicas') < 0
      && en.indexOf('grupo externo') < 0 && en.indexOf('penalizada') < 0;
  })());
  t('ni el párrafo en español palabras en inglés', (() => {
    const es = bi.map(s => s.es.join(' ')).join(' ');
    return es.indexOf('bootstrap (') < 0 && es.indexOf('by outgroup') < 0;
  })());
  t('el texto plano usa el título del idioma que se le pide',
    Report.methodsText(bi, 'en').indexOf('Maximum likelihood') >= 0
    && Report.methodsText(bi, 'es').indexOf('Máxima verosimilitud') >= 0);
  t('y no deja etiquetas HTML en el texto plano',
    Report.methodsText(bi, 'es').indexOf('<b>') < 0);

  sec('120 · El informe como documento');
  const full = {
    data: { taxa: L6, parts: [{ name: 'g', taxa: L6, seqs: L6.map(() => 'ACGT'.repeat(50)) }] },
    models: { best: { name: 'GTR+G', lnL: -1234.5 }, criterion: 'BIC' },
    ml: { lnL: -1200.1, support: { boot: { reps: 100 } } },
  };
  const secs = Report.methods(full, {});
  const html = Report.html(full, { lang: 'es', sections: secs, figures: [{ svg, caption: 'un árbol' }] });
  t('es un documento completo', /^<!DOCTYPE html>/.test(html) && /<\/html>\s*$/.test(html));
  t('no pide ningún archivo externo', html.indexOf('<link') < 0 && html.indexOf('<script') < 0);
  t('ni ninguna imagen por URL', html.indexOf('<img ') < 0);
  t('lleva su propio estilo dentro', html.indexOf('<style>') >= 0);
  t('trae la figura', html.indexOf('<figure>') >= 0 && html.indexOf('un árbol') >= 0);
  t('trae la tabla de resumen', html.indexOf('<table>') >= 0);
  t('y la bibliografía', html.indexOf('Referencias') >= 0);
  t('la versión en inglés cambia los encabezados', (() => {
    const en = Report.html(full, { lang: 'en', sections: secs });
    return en.indexOf('References') >= 0 && en.indexOf('Referencias') < 0
      && en.indexOf('Methods') >= 0;
  })());
  t('un título con < no rompe el documento', (() => {
    const h = Report.html(full, { lang: 'es', sections: secs, title: 'a < b' });
    return h.indexOf('a &lt; b') >= 0;
  })());
  t('el resumen sólo trae las filas que existen', (() => {
    const rows = Report.summary({ data: full.data }, 'es');
    return rows.every(r => r.value != null && r.value !== '')
      && !rows.some(r => /verosimilitud/i.test(r.label));
  })());
  t('y crece cuando hay más', Report.summary(full, 'es').length > Report.summary({ data: full.data }, 'es').length);

  sec('121 · Lo que el Bloque 12 expone');
  t('el bloque ofrece sus funciones',
    ['drawStudio', 'refreshFormats', 'buildReport', 'buildZip', 'trees'].every(k => typeof B12[k] === 'function'));
  t('el módulo de exportación ofrece las suyas',
    ['newick', 'nexusTrees', 'phyloXML', 'readPhyloXML', 'pngSize', 'svgToPng', 'crc32', 'zip']
      .every(k => typeof Export[k] === 'function'));
  t('y el del informe las suyas',
    ['methods', 'references', 'methodsText', 'html', 'summary'].every(k => typeof Report[k] === 'function'));
  t('la bibliografía que el programa conoce tiene más de cuarenta obras',
    Object.keys(Report.REFS).length > 40, String(Object.keys(Report.REFS).length));
  t('y ninguna entrada vacía',
    Object.keys(Report.REFS).every(k => Report.REFS[k] && Report.REFS[k].length > 40));

})();
