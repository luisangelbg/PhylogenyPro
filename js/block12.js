/* PhylogenyPro — molecular phylogenetics in the browser.
   Copyright (C) 2026  Luis Ángel Barrera-Guzmán

   This program is free software: you can redistribute it and/or modify it under
   the terms of the GNU General Public License as published by the Free Software
   Foundation, either version 3 of the License, or (at your option) any later
   version. It is distributed in the hope that it will be useful, but WITHOUT ANY
   WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
   A PARTICULAR PURPOSE. See the GNU General Public License, in the file LICENSE
   at the root of this program, or <https://www.gnu.org/licenses/>. */

/* PhylogenyPro — Block 12: figures, formats and the report.

   The engines are js/exportfmt.js and js/report.js. What this file adds is the
   studio — every decision a tree figure contains, made visible and reversible —
   and the honesty of the export:

     · a PNG at 900 dpi is a specific number of pixels, and that number is shown
       before anything is generated, because a browser canvas has a ceiling and
       a silently blank image is worse than a refusal;
     · the report says what ran. A block that was never opened leaves no
       paragraph, and the bibliography holds only what was cited;
     · the ZIP is written here, byte by byte, so the package works from a
       double-clicked file with no network and no library. */

(function () {

  const B12 = {
    opts: null, tree: null, colours: null, report: null, zip: null, figs: null,
  };
  window.B12 = B12;

  /* ================================================================
     what there is to work with
     ================================================================ */
  function taxa() {
    if (state.data && state.data.parts && state.data.parts.length) {
      const p = state.data.parts[0];
      if (p && p.taxa) return p.taxa;
    }
    return (state.data && state.data.taxa) || [];
  }
  function trees() {
    const out = [];
    const add = (id, name, tree, extra) => { if (tree) out.push(Object.assign({ id, name, tree }, extra || {})); };
    if (state.quick && state.quick.distance) add('distance', T('Bloque 4 · distancias', 'Block 4 · distance'), state.quick.distance.tree);
    if (state.quick && state.quick.parsimony) add('parsimony', T('Bloque 4 · parsimonia', 'Block 4 · parsimony'), state.quick.parsimony.tree);
    if (state.quick && state.quick.consensus) add('consensus4', T('Bloque 4 · consenso', 'Block 4 · consensus'), state.quick.consensus.tree);
    if (state.ml) add('ml', T('Bloque 5 · máxima verosimilitud', 'Block 5 · maximum likelihood'), state.ml.tree,
      { supportType: 'bootstrap' });
    if (state.bayes) add('bayes', T('Bloque 6 · consenso bayesiano', 'Block 6 · Bayesian consensus'), state.bayes.consensus,
      { supportType: 'posterior probability' });
    if (state.dated) add('dated', T('Bloque 7 · árbol fechado', 'Block 7 · dated tree'), state.dated.tree,
      { unit: state.dated.unit || 'Ma' });
    if (state.compare && state.compare.quartet && state.compare.quartet.newick) {
      try {
        add('quartet', T('Bloque 11 · árbol por cuartetos', 'Block 11 · quartet species tree'),
          Tree.parseNewick(state.compare.quartet.newick, taxa()));
      } catch (e) { /* a tree that will not parse is simply not offered */ }
    }
    return out;
  }
  const treeById = id => trees().find(t => t.id === id) || trees()[0] || null;
  function figures() {
    const f = state.figures || {};
    return Object.keys(f).sort().map(k => Object.assign({ key: k }, f[k]));
  }

  /* ================================================================
     1 · the studio
     ================================================================ */
  function studioOptions() {
    const o = {
      layout: el('p12Layout').value,
      cladogram: el('p12Clado').checked,
      showSupport: el('p12SupportAs').value !== 'none',
      supportAs: el('p12SupportAs').value === 'none' ? 'node' : el('p12SupportAs').value,
      supportMin: +el('p12SupportMin').value || 0,
      width: Math.max(300, Math.min(2000, +el('p12Width').value || 760)),
      rowHeight: Math.max(8, Math.min(40, +el('p12RowHeight').value || 18)),
      tipFont: Math.max(5, Math.min(24, +el('p12TipFont').value || 12)),
      lineWidth: Math.max(0.5, Math.min(6, +el('p12LineWidth').value || 1.8)),
      italicTips: el('p12Italic').checked,
      showScale: el('p12Scale').checked,
      labels: taxa(),
      title: el('p12Title').value || null,
      colours: B12.colours,
      alignTips: el('p12AlignTips') ? el('p12AlignTips').checked : false,
    };
    /* groups: the cut is recomputed from the tree as the studio has shaped it,
       so rerooting or ladderising never leaves a stale colouring behind */
    if (B12.groups) {
      const show = el('p12GroupShow') ? el('p12GroupShow').value : 'branches';
      o.groups = B12.groups;
      o.groupColours = window.FigStyle ? FigStyle.colours().slice(0, Math.max(1, B12.groups.k))
        : Groups.palette(B12.groups.k);
      o.colourBranches = true;
      o.colourLabels = true;
      o.groupStrip = show === 'strip' || show === 'bands';
      o.groupBands = show === 'bands';
      o.legend = el('p12Legend') ? el('p12Legend').checked : true;
    }
    if (el('p12Images')) {
      o.images = el('p12Images').value;
      o.imgSize = Math.max(18, Math.min(180, +el('p12ImgSize').value || 54));
      o.imgShape = el('p12ImgShape').value;
      o.imgFrame = el('p12ImgFrame').checked;
      o.cladeImages = B12.cladeImages || null;
    }
    return o;
  }

  /* ================================================================
     1b · the figure studio
     ================================================================ */
  /* The controls only write into FigStyle; FigStyle writes one CSS rule and
     every figure of the twelve blocks follows it. Nothing here knows how a
     tree or a curve is drawn. */
  function fillStyleControls() {
    if (!el('p12Palette') || !window.FigStyle) return;
    const st = FigStyle.get();

    el('p12Palette').innerHTML = FigStyle.PALETTES.map(p =>
      `<option value="${p.id}"${p.id === st.palette ? ' selected' : ''} data-es="${esc(p.name[0])}" data-en="${esc(p.name[1])}">${esc(T(p.name[0], p.name[1]))}</option>`).join('');
    el('p12Font').innerHTML = FigStyle.FONTS.map(f =>
      `<option value="${f.id}"${f.id === st.font ? ' selected' : ''} data-es="${esc(f.name[0])}" data-en="${esc(f.name[1])}" style="font-family:${f.css}">${esc(T(f.name[0], f.name[1]))}</option>`).join('');
    el('p12Background').innerHTML = FigStyle.BACKGROUNDS.map(b =>
      `<option value="${b.id}"${b.id === st.background ? ' selected' : ''} data-es="${esc(b.name[0])}" data-en="${esc(b.name[1])}">${esc(T(b.name[0], b.name[1]))}</option>`).join('');

    renderSwatches();
    renderPaletteNote();
  }

  function renderPaletteNote() {
    const n = el('p12PaletteNote');
    if (!n || !window.FigStyle) return;
    const p = FigStyle.paletteOf();
    n.innerHTML = `<span class="pal-strip">${FigStyle.colours().map(c => `<i style="background:${esc(c)}"></i>`).join('')}</span> ` +
      L2(p.note[0], p.note[1]);
  }

  function renderSwatches() {
    const box = el('p12Swatches');
    if (!box || !window.FigStyle) return;
    const st = FigStyle.get();
    const cols = FigStyle.colours();
    box.innerHTML = '';
    cols.forEach((c, i) => {
      const key = 'c' + (i + 1);
      const w = mk('label', { class: 'swatch' + (st.custom[key] ? ' edited' : '') });
      w.innerHTML = `<input type="color" value="${esc(toHex(c))}"><span>${i + 1}</span>`;
      w.querySelector('input').addEventListener('input', e => {
        FigStyle.set({ custom: { [key]: e.target.value } });
      });
      box.appendChild(w);
    });

    const tb = el('p12TokenSwatches');
    if (!tb) return;
    const NAMES = {
      text: ['texto', 'text'], 'text-muted': ['texto tenue', 'muted text'], ink: ['curvas', 'curves'],
      border: ['líneas finas', 'thin lines'], 'border-strong': ['ejes', 'axes'],
      accent: ['acento', 'accent'], leaf: ['apoyo alto', 'high support'],
      gold: ['apoyo medio', 'mid support'], rose: ['apoyo bajo', 'low support'],
    };
    tb.innerHTML = '';
    Object.keys(NAMES).forEach(k => {
      const cur = st.tokens[k] || readToken(k);
      const w = mk('label', { class: 'swatch' + (st.tokens[k] ? ' edited' : '') });
      w.innerHTML = `<input type="color" value="${esc(toHex(cur))}"><span>${esc(T(NAMES[k][0], NAMES[k][1]))}</span>`;
      w.querySelector('input').addEventListener('input', e => {
        FigStyle.set({ tokens: { [k]: e.target.value } });
      });
      tb.appendChild(w);
    });
  }

  const readToken = k => getComputedStyle(document.documentElement).getPropertyValue('--' + k).trim() || '#888888';

  /* an <input type="color"> only understands #rrggbb, and the theme is written
     in every notation CSS allows */
  function toHex(c) {
    c = String(c || '').trim();
    if (/^#[0-9a-f]{6}$/i.test(c)) return c.toLowerCase();
    if (/^#[0-9a-f]{3}$/i.test(c)) return '#' + c[1] + c[1] + c[2] + c[2] + c[3] + c[3];
    try {
      const d = document.createElement('span');
      d.style.color = c;
      document.body.appendChild(d);
      const rgb = getComputedStyle(d).color.match(/\d+/g);
      d.remove();
      if (rgb) return '#' + rgb.slice(0, 3).map(v => (+v).toString(16).padStart(2, '0')).join('');
    } catch (e) { /* fall through */ }
    return '#888888';
  }

  function initStyleControls() {
    if (!el('p12Palette') || !window.FigStyle) return;
    fillStyleControls();
    el('p12Palette').addEventListener('change', e => { FigStyle.set({ palette: e.target.value }); renderSwatches(); renderPaletteNote(); });
    el('p12Font').addEventListener('change', e => FigStyle.set({ font: e.target.value }));
    el('p12Background').addEventListener('change', e => FigStyle.set({ background: e.target.value }));
    const scale = (id, out, key, fmt) => {
      const r = el(id);
      if (!r) return;
      r.addEventListener('input', () => {
        el(out).textContent = fmt(+r.value);
        FigStyle.set({ [key]: +r.value });
      });
      el(out).textContent = fmt(+r.value);
    };
    scale('p12FontScale', 'p12FontScaleV', 'fontScale', v => v.toFixed(2) + '×');
    scale('p12LineScale', 'p12LineScaleV', 'lineScale', v => v.toFixed(1) + '×');
    el('p12StyleReset').addEventListener('click', () => {
      FigStyle.reset();
      el('p12FontScale').value = 1; el('p12FontScaleV').textContent = '1.00×';
      el('p12LineScale').value = 1; el('p12LineScaleV').textContent = '1.0×';
      fillStyleControls();
      drawStudio();
    });
    /* the tree has to be redrawn because its group colours are read at draw
       time; everything else follows the CSS rule on its own */
    FigStyle.onChange(() => { renderPaletteNote(); if (el('p12TreeFig')) drawStudio(); });
    document.addEventListener('langchange', fillStyleControls);
  }

  /* ---------------- groups ---------------- */
  /* Recomputed on the shaped tree, and reported: how many groups actually came
     out (a cut cannot always give exactly the k asked for), and whether each
     one is a clade. */
  function computeGroups(tree) {
    const msg = el('p12GroupMsg');
    if (msg) msg.innerHTML = '';
    const how = el('p12GroupMethod') ? el('p12GroupMethod').value : 'none';
    if (how === 'none' || !tree) { B12.groups = null; return; }
    const k = Math.max(1, Math.min(24, +el('p12K').value || 3));
    let g;
    if (how === 'depth') {
      if (!Groups.isUltrametric(tree, 1e-6)) {
        showMessage(msg, 'warn', L2(
          'Cortar por profundidad sólo tiene sentido en un árbol ultramétrico —un cronograma o un UPGMA—, y éste no lo es: ' +
          'las puntas no llegan todas a la misma altura. Se cortó por las ramas más largas.',
          'Cutting by depth only makes sense on an ultrametric tree — a chronogram or a UPGMA — and this one is not: ' +
          'the tips do not all reach the same height. It was cut by the longest branches instead.'));
        g = Groups.byLongestBranches(tree, k);
      } else {
        g = Groups.byDepth(tree, k);
        if (g.exact === false) {
          showMessage(msg, 'info', L2(
            `Ningún corte por profundidad deja exactamente ${k} grupos en este árbol; el más cercano deja ${g.k}.`,
            `No cut by depth leaves exactly ${k} groups on this tree; the nearest one leaves ${g.k}.`));
        }
      }
    } else {
      g = Groups.byLongestBranches(tree, k);
    }
    if (B12.groupNames) g.names = B12.groupNames;
    B12.groups = g;

    const sz = Groups.sizes(g);
    const mono = Groups.monophyly(tree, g);
    const notMono = mono.filter(m => !m.monophyletic && m.size > 0);
    const parts = sz.map((n, i) => `${Groups.nameOf(g, i)} (${n})`).join(', ');
    showMessage(msg, notMono.length ? 'warn' : 'ok', L2(
      `${g.k} grupo(s): ${parts}.` + (notMono.length
        ? ` Ojo: ${notMono.map(m => Groups.nameOf(g, m.group)).join(', ')} no ${notMono.length === 1 ? 'es' : 'son'} monofilético(s), así que el color no marca un clado.`
        : ' Todos son clados.'),
      `${g.k} group(s): ${parts}.` + (notMono.length
        ? ` Note: ${notMono.map(m => Groups.nameOf(g, m.group)).join(', ')} ${notMono.length === 1 ? 'is' : 'are'} not monophyletic, so the colour does not mark a clade.`
        : ' Every one of them is a clade.')));
  }
  /* the tree as the studio's rooting and ladderising leave it */
  function shapedTree() {
    const src = treeById(el('p12Tree').value);
    if (!src) return null;
    let t = Tree.clone(src.tree);
    const root = el('p12Root').value;
    if (root === 'midpoint') t = Tree.midpointRoot(t);
    else if (root === 'outgroup') {
      const picked = [...el('p12Outgroup').selectedOptions].map(o => +o.value);
      if (picked.length) {
        const r = Tree.rootByOutgroup(t, picked);
        if (r && r.tree) {
          t = r.tree;
          if (!r.monophyletic) {
            showMessage(el('p12StudioMsg'), 'warn', L2(
              'El grupo externo elegido no es monofilético en este árbol, así que no se puede enraizar sobre él; el árbol queda como estaba.',
              'The chosen outgroup is not monophyletic on this tree, so it cannot be rooted on it; the tree is left as it was.'));
          }
        }
      }
    }
    const lad = el('p12Ladder').value;
    if (lad !== 'none') Tree.ladderize(t, lad === 'up');
    return Object.assign({}, src, { tree: t });
  }
  function drawStudio() {
    el('p12StudioMsg').innerHTML = '';
    const src = shapedTree();
    if (!src) return;
    B12.tree = src;
    computeGroups(src.tree);
    B12.opts = studioOptions();
    const svg = TreeView.render(src.tree, B12.opts);
    el('p12TreeFig').innerHTML = svg;
    reportFigure();
    state.figures = state.figures || {};
    state.figures.b12_tree = { block: 12, name: 'tree', svg, caption: src.name };
    refreshFormats();
    refreshFigures();
  }

  /* What the figure could not show, said out loud instead of silently dropped. */
  function reportFigure() {
    const msg = el('p12ImgMsg');
    if (!msg) return;
    msg.innerHTML = '';
    const info = TreeView.lastInfo || {};
    const notes = [];
    if (info.hiddenLabels) {
      notes.push(L2(
        `Se ocultaron ${info.hiddenLabels} nombre(s) de punta porque se encimaban con otros. En un árbol no enraizado ` +
        'nada impide que dos puntas caigan juntas: agranda la figura o baja el tamaño de la letra.',
        `${info.hiddenLabels} tip name(s) were hidden because they overlapped. Nothing keeps two tips of an unrooted ` +
        'tree from landing beside each other: make the figure larger or the font smaller.'));
    }
    if (info.hiddenSupport) {
      notes.push(L2(
        `Se ocultaron ${info.hiddenSupport} valor(es) de apoyo porque se encimaban. Sube el umbral, agranda el alto de fila, ` +
        'o cambia el apoyo a puntos.',
        `${info.hiddenSupport} support value(s) were hidden because they overlapped. Raise the threshold, increase the row ` +
        'height, or show support as dots.'));
    }
    const want = el('p12Images') ? el('p12Images').value : 'none';
    if (want === 'tips' && window.OTUImg) {
      const names = taxa() || [];
      const missing = names.filter(n => n && !OTUImg.url(n)).length;
      if (missing) {
        notes.push(L2(
          `${missing} de ${names.length} puntas no tienen imagen confirmada todavía.`,
          `${missing} of ${names.length} tips have no confirmed picture yet.`));
      }
    }
    const cr = window.OTUImg ? OTUImg.credits() : [];
    if (cr.length && want !== 'none') {
      notes.push(L2(`Créditos de las imágenes: ${cr.join(' · ')}.`, `Picture credits: ${cr.join(' · ')}.`));
    }
    notes.forEach((t, i) => showMessage(msg, i === 0 && (info.hiddenSupport || info.hiddenLabels) ? 'warn' : 'info', t));
  }

  function refreshTrees() {
    const list = trees();
    const sel = el('p12Tree');
    const prev = sel.value;
    sel.innerHTML = list.map(t => `<option value="${t.id}">${esc(t.name)}</option>`).join('');
    if (list.some(t => t.id === prev)) sel.value = prev;
    const names = taxa();
    ['p12Outgroup', 'p12ColourClade'].forEach(id => {
      const s = el(id);
      s.innerHTML = names.map((nm, i) => `<option value="${i}">${esc(nm)}</option>`).join('');
    });
  }

  /* ================================================================
     2 · tree formats
     ================================================================ */
  function formatText() {
    const list = el('p12AllTrees').checked ? trees().map(t => {
      const shaped = t.id === el('p12Tree').value ? shapedTree() : t;
      return shaped || t;
    }) : [shapedTree()].filter(Boolean);
    if (!list.length) return '';
    const names = taxa();
    const dec = Math.max(0, Math.min(12, +el('p12Decimals').value));
    const nwkOpts = {
      labels: names, decimals: dec,
      lengths: el('p12WithLengths').checked,
      support: el('p12WithSupport').checked,
    };
    const fmt = el('p12Format').value;
    if (fmt === 'newick') return list.map(t => Export.newick(t.tree, nwkOpts)).join('\n') + '\n';
    if (fmt === 'nexus') {
      return Export.nexusTrees(list.map(t => ({ name: t.name, tree: t.tree })),
        { labels: names, taxa: names, newick: nwkOpts, figtree: true });
    }
    return Export.phyloXML(list.map(t => ({
      name: t.name, tree: t.tree, supportType: t.supportType, unit: t.unit,
    })), { labels: names });
  }
  function refreshFormats() {
    const txt = formatText();
    el('p12FormatText').value = txt;
    const lines = txt ? txt.split('\n').length : 0;
    el('p12FormatInfo').textContent = T(`${(txt.length / 1024).toFixed(1)} kB, ${lines} líneas`,
      `${(txt.length / 1024).toFixed(1)} kB, ${lines} lines`);
  }
  function downloadFormat() {
    const txt = formatText();
    if (!txt) return;
    const stamp = new Date().toISOString().slice(0, 10);
    const ext = el('p12Format').value === 'newick' ? 'tre' : el('p12Format').value === 'nexus' ? 'nex' : 'xml';
    const mime = ext === 'xml' ? 'application/xml' : 'text/plain';
    download(txt, `phylogenypro_arbol_${stamp}.${ext}`, mime);
  }

  /* ================================================================
     3 · figures
     ================================================================ */
  function refreshFigures() {
    const figs = figures();
    B12.figs = figs;
    const dpi = +el('p12Dpi').value || 300;
    buildTable('p12FigTable', [
      { key: 'name', label: T('figura', 'figure') },
      { key: 'block', label: T('bloque', 'block'), num: true },
      { key: 'size', label: T('tamaño del dibujo', 'drawing size'), num: true },
      { key: 'px', label: T(`píxeles a ${dpi} ppp`, `pixels at ${dpi} dpi`), num: true },
      { key: 'act', label: T('descargar', 'download'), html: true },
    ], figs.map(f => {
      const s = Export.pngSize(f.svg, { dpi });
      return {
        name: f.caption || f.name, block: f.block,
        size: `${Math.round(s.baseWidth)} × ${Math.round(s.baseHeight)}`,
        px: `${s.width} × ${s.height}`,
        act: `<button class="btn btn-ghost btn-sm" data-p12svg="${esc(f.key)}">SVG</button> ` +
          `<button class="btn btn-ghost btn-sm" data-p12png="${esc(f.key)}"${s.tooBig ? ' disabled' : ''}>PNG</button>`,
        _class: s.tooBig ? 'row-flag' : '',
      };
    }), { limit: 40 });
    els('[data-p12svg]').forEach(b => b.addEventListener('click', () => saveSvg(b.dataset.p12svg)));
    els('[data-p12png]').forEach(b => b.addEventListener('click', () => savePng(b.dataset.p12png)));
    const big = figs.filter(f => Export.pngSize(f.svg, { dpi }).tooBig);
    const msg = el('p12FigMsg');
    msg.innerHTML = '';
    if (big.length) {
      showMessage(msg, 'warn', L2(
        `A ${dpi} ppp, ${big.length} figura(s) pasarían de los ${(Export.MAX_PIXELS / 1e6).toFixed(0)} millones de píxeles que un lienzo de navegador admite, así que su PNG queda desactivado. El SVG no tiene ese límite y es lo que conviene mandar; si la revista exige PNG, basta con bajar la resolución o el ancho del dibujo.`,
        `At ${dpi} dpi, ${big.length} figure(s) would exceed the ${(Export.MAX_PIXELS / 1e6).toFixed(0)} million pixels a browser canvas allows, so their PNG is disabled. The SVG has no such limit and is what to send; if the journal insists on PNG, lowering the resolution or the drawing width is enough.`));
    }
  }
  function saveSvg(key) {
    const f = (state.figures || {})[key];
    if (!f) return;
    download('<?xml version="1.0" encoding="UTF-8"?>\n' + Export.inlineVars(f.svg, Export.readVars()),
      `phylogenypro_${key}.svg`, 'image/svg+xml');
  }
  function savePng(key) {
    const f = (state.figures || {})[key];
    if (!f) return;
    const dpi = +el('p12Dpi').value || 300;
    const msg = el('p12FigMsg');
    Export.svgToPng(f.svg, { dpi }).then(r => {
      const url = URL.createObjectURL(r.blob);
      const a = document.createElement('a');
      a.href = url; a.download = `phylogenypro_${key}_${dpi}dpi.png`;
      document.body.appendChild(a); a.click(); a.remove();
      setTimeout(() => URL.revokeObjectURL(url), 2000);
    }).catch(e => {
      showMessage(msg, 'error', L2(
        `No se pudo generar el PNG (${e.message}). El SVG sí se puede descargar y no pierde nada al imprimirse.`,
        `The PNG could not be generated (${e.message}). The SVG can be downloaded and loses nothing in print.`));
    });
  }

  /* ================================================================
     4 · the report
     ================================================================ */
  function buildReport() {
    const lang = el('p12Lang').value;
    const sections = Report.methods(state, {});
    const refs = Report.references(sections);
    const figs = el('p12RepFigures').checked
      ? figures().map(f => ({ svg: Export.inlineVars(f.svg, Export.readVars()), caption: f.caption || f.name, name: f.name }))
      : [];
    const html = Report.html(state, {
      lang, sections, figures: figs,
      title: el('p12RepTitle').value || undefined,
      author: el('p12RepAuthor').value || undefined,
      tables: reportTables(lang),
    });
    B12.report = { html, sections, refs, lang };
    statTiles('p12RepTiles', [
      [T('Secciones de métodos', 'Methods sections'), sections.length, ''],
      [T('Referencias citadas', 'References cited'), refs.length,
        T(`de ${Object.keys(Report.REFS).length} que el programa conoce`, `of ${Object.keys(Report.REFS).length} the program knows`)],
      [T('Figuras incluidas', 'Figures included'), figs.length, ''],
      [T('Tamaño', 'Size'), `${(html.length / 1024).toFixed(0)} kB`, T('un solo archivo', 'a single file')],
    ]);
    el('p12RepPane').style.display = '';
    const frame = el('p12RepFrame');
    frame.srcdoc = html;
    el('p12Zip').style.display = '';
    const msg = el('p12RepMsg');
    msg.innerHTML = '';
    const missing = [];
    if (!state.ml && !state.bayes) missing.push(T('una inferencia (Bloques 5 o 6)', 'an inference (Blocks 5 or 6)'));
    if (!state.models) missing.push(T('la selección de modelo (Bloque 3)', 'the model selection (Block 3)'));
    if (missing.length) {
      showMessage(msg, 'info', L2(
        `El informe sólo describe lo que se corrió, así que todavía no menciona ${missing.join(' ni ')}. No es un error: si esos bloques no se usaron, no deben aparecer en los métodos.`,
        `The report describes only what was run, so it does not yet mention ${missing.join(' or ')}. That is not an error: if those blocks were not used, they have no business in the methods.`));
    }
  }
  /* the tables worth carrying into a supplementary file */
  function reportTables(lang) {
    const en = lang === 'en';
    const out = [];
    if (state.models && state.models.fits && state.models.fits.length) {
      const top = state.models.fits.slice(0, 10);
      out.push({
        caption: en ? 'Substitution models, best ten by the chosen criterion'
          : 'Modelos de sustitución, los diez mejores por el criterio elegido',
        columns: [en ? 'model' : 'modelo', 'lnL', 'k', 'AIC', 'BIC'],
        rows: top.map(f => [f.name || f.model, fmtLnL(f.lnL), f.k, fmtFixed(f.AIC, 3), fmtFixed(f.BIC, 3)]),
      });
    }
    if (state.biogeo && state.biogeo.models) {
      out.push({
        caption: en ? 'Biogeographic models' : 'Modelos biogeográficos',
        columns: [en ? 'model' : 'modelo', 'd', 'e', 'j', 'lnL', 'AICc', 'ΔAICc', en ? 'weight' : 'peso'],
        rows: state.biogeo.models.map(m => [m.name, fmtNum(m.d), fmtNum(m.e), m.withJ ? fmtNum(m.j) : '—',
          fmtLnL(m.lnL), fmtFixed(m.AICc, 3), fmtFixed(m.dAICc, 3), fmtFixed(m.w, 3)]),
      });
    }
    if (state.compare && state.compare.concordance) {
      const br = state.compare.concordance.branches.slice().sort((a, b) => a.gCF - b.gCF).slice(0, 15);
      out.push({
        caption: en ? 'Concordance factors, the fifteen least concordant branches'
          : 'Factores de concordancia, las quince ramas menos concordantes',
        columns: [en ? 'clade' : 'clado', en ? 'tips' : 'puntas', 'gCF', 'gDF1', 'gDF2', 'sCF'],
        rows: br.map(b => [b.clade.slice(0, 3).join(', ') + (b.clade.length > 3 ? '…' : ''),
          b.nTips, fmtFixed(b.gCF, 1), fmtFixed(b.gDF1, 1), fmtFixed(b.gDF2, 1),
          b.sCF == null ? '—' : fmtFixed(b.sCF, 1)]),
      });
    }
    if (state.traits && state.traits.continuous) {
      out.push({
        caption: en ? 'Models of continuous character evolution' : 'Modelos de evolución del carácter continuo',
        columns: [en ? 'model' : 'modelo', 'lnL', 'AIC', 'ΔAIC', en ? 'weight' : 'peso'],
        rows: state.traits.continuous.models.map(m => [m.model, fmtLnL(m.lnL),
          fmtFixed(m.AIC, 3), fmtFixed(m.dAIC, 3), fmtFixed(m.w, 3)]),
      });
    }
    return out;
  }

  /* ================================================================
     5 · the package
     ================================================================ */
  async function buildZip() {
    const msg = el('p12ZipMsg');
    msg.innerHTML = '';
    const bar = window.LABG ? LABG.progressBar(el('p12ZipProgress'), { label: T('Paquete reproducible', 'Reproducible package') }) : null;
    if (bar) bar.update(null, T('armando…', 'building…'));
    else el('p12ZipProgress').innerHTML = L2('armando…', 'building…');
    el('p12BuildZip').disabled = true;
    try {
      /* let the bar show before the files are gathered, which takes the page */
      if (bar) await LABG.nextPaint();
      const stamp = new Date().toISOString().slice(0, 10);
      const names = taxa();
      const files = [];
      const manifest = [];
      const put = (name, content, what) => { files.push({ name, content }); manifest.push({ name, what, size: content.length }); };

      if (el('p12ZipData').checked && state.data && state.data.parts) {
        state.data.parts.forEach(p => {
          if (!p.seqs || !p.seqs.length) return;
          put(`datos/${slug(p.name)}.fasta`, SeqIO.writeFasta(p.taxa, p.seqs),
            T(`alineamiento de la partición ${p.name}`, `alignment of partition ${p.name}`));
        });
        if (state.data.parts.length > 1) {
          put('datos/particiones.txt',
            state.data.parts.map(p => `${p.type === 'aa' ? 'PROT' : 'DNA'}, ${p.name} = ?`).join('\n'),
            T('esquema de particiones', 'partition scheme'));
        }
      }
      if (el('p12ZipTrees').checked) {
        const list = trees();
        if (list.length) {
          put('arboles/arboles.tre', list.map(t => Tree.writeNewick(t.tree, { labels: names })).join('\n') + '\n',
            T('todos los árboles, Newick', 'every tree, Newick'));
          put('arboles/arboles.nex', Export.nexusTrees(list.map(t => ({ name: t.name, tree: t.tree })),
            { labels: names, taxa: names, figtree: true }),
            T('todos los árboles, NEXUS con bloque de FigTree', 'every tree, NEXUS with a FigTree block'));
          put('arboles/arboles.xml', Export.phyloXML(list.map(t => ({
            name: t.name, tree: t.tree, supportType: t.supportType, unit: t.unit })), { labels: names }),
            T('todos los árboles, phyloXML con apoyos y edades', 'every tree, phyloXML with supports and ages'));
        }
      }
      if (el('p12ZipFigures').checked) {
        const vars = Export.readVars();
        figures().forEach(f => {
          put(`figuras/${f.key}.svg`, '<?xml version="1.0" encoding="UTF-8"?>\n' + Export.inlineVars(f.svg, vars),
            f.caption || f.name);
        });
      }
      if (el('p12ZipTables').checked) {
        tablesForZip().forEach(t => put(`tablas/${t.name}`, t.content, t.what));
      }
      if (el('p12ZipReport').checked) {
        if (!B12.report) buildReport();
        put('informe.html', B12.report.html, T('informe con los métodos redactados', 'report with the methods written'));
        put('metodos.txt', Report.methodsText(B12.report.sections, B12.report.lang) + '\n\n' +
          (B12.report.lang === 'en' ? 'REFERENCES\n' : 'REFERENCIAS\n') +
          B12.report.refs.map(r => r.text).join('\n'),
          T('los métodos en texto plano, para pegar', 'the methods as plain text, to paste'));
      }
      /* the file that says what everything is: a package nobody can read is a
         backup, not a reproducible analysis */
      const readme = buildReadme(manifest, stamp);
      files.unshift({ name: 'LEEME.txt', content: readme });
      manifest.unshift({ name: 'LEEME.txt', what: T('este índice', 'this index'), size: readme.length });

      const blob = await Export.zip(files);
      B12.zip = { blob, manifest, files: files.length };
      const url = URL.createObjectURL(blob);
      const a = document.createElement('a');
      a.href = url; a.download = `phylogenypro_paquete_${stamp}.zip`;
      document.body.appendChild(a); a.click(); a.remove();
      setTimeout(() => URL.revokeObjectURL(url), 4000);

      statTiles('p12ZipTiles', [
        [T('Archivos', 'Files'), files.length, ''],
        [T('Tamaño del paquete', 'Package size'), `${(blob.size / 1024).toFixed(0)} kB`, ''],
        [T('Sin comprimir', 'Uncompressed'), `${(manifest.reduce((a, m) => a + m.size, 0) / 1024).toFixed(0)} kB`, ''],
        [T('Escrito', 'Written'), T('por la propia app', 'by the app itself'), T('sin red ni biblioteca', 'no network, no library')],
      ]);
      buildTable('p12ZipTable', [
        { key: 'name', label: T('archivo', 'file') },
        { key: 'what', label: T('qué es', 'what it is') },
        { key: 'size', label: T('bytes', 'bytes'), num: true },
      ], manifest.map(m => ({ name: m.name, what: m.what, size: m.size })), { limit: 60 });
      if (bar) bar.done(T(`paquete listo · ${files.length} archivos`, `package ready · ${files.length} files`));
      else el('p12ZipProgress').innerHTML = L2('listo', 'done');
    } catch (e) {
      showMessage(msg, 'error', L2(`No se pudo armar el paquete: ${e.message}`,
        `The package could not be built: ${e.message}`));
      if (bar) bar.fail(T('no se pudo armar', 'could not be built'));
      else el('p12ZipProgress').innerHTML = '';
    }
    el('p12BuildZip').disabled = false;
  }
  function tablesForZip() {
    const out = [];
    const csv = rows => rows.map(r => r.join(',')).join('\n');
    if (state.biogeo && state.biogeo.models) {
      out.push({ name: 'biogeografia_modelos.csv', what: T('modelos biogeográficos', 'biogeographic models'),
        content: csv([['model', 'd', 'e', 'j', 'lnL', 'k', 'AICc', 'dAICc', 'weight']].concat(
          state.biogeo.models.map(m => [m.name, m.d, m.e, m.j, m.lnL, m.k, m.AICc, m.dAICc, m.w]))) });
    }
    if (state.compare && state.compare.concordance) {
      out.push({ name: 'concordancia.csv', what: T('factores de concordancia', 'concordance factors'),
        content: csv([['clade', 'n_tips', 'gCF', 'gDF1', 'gDF2', 'sCF']].concat(
          state.compare.concordance.branches.map(b => ['"' + b.clade.join(' ') + '"', b.nTips, b.gCF, b.gDF1, b.gDF2, b.sCF]))) });
    }
    if (state.traits && state.traits.continuous) {
      out.push({ name: 'caracteres_modelos.csv', what: T('modelos de carácter continuo', 'continuous character models'),
        content: csv([['model', 'lnL', 'k', 'AIC', 'dAIC', 'weight', 'sigma2']].concat(
          state.traits.continuous.models.map(m => [m.model, m.lnL, m.k, m.AIC, m.dAIC, m.w, m.sigma2]))) });
    }
    if (state.diversification && state.diversification.models) {
      /* the model names carry commas ("density dependent, logistic"), so they go quoted */
      out.push({ name: 'diversificacion_modelos.csv', what: T('modelos de diversificación', 'diversification models'),
        content: csv([['model', 'lnL', 'k', 'AIC', 'dAIC', 'weight']].concat(
          state.diversification.models.map(m => ['"' + m.name + '"', m.loglik, m.k, m.AIC, m.dAIC, m.w]))) });
    }
    return out;
  }
  function buildReadme(manifest, stamp) {
    const es = L2on();
    const L = [];
    L.push(es ? 'PAQUETE DE ANÁLISIS — PhylogenyPro' : 'ANALYSIS PACKAGE — PhylogenyPro');
    L.push('='.repeat(48));
    L.push('');
    L.push((es ? 'Generado el ' : 'Produced on ') + stamp);
    L.push('');
    L.push(es
      ? 'Todo lo que hay aquí lo escribió el programa a partir de la corrida, no de\nuna plantilla. Los árboles vienen en los tres formatos porque cada programa\nlee uno distinto, y el informe describe únicamente los análisis que se\nhicieron: si un bloque no se usó, no aparece.'
      : 'Everything here was written by the program from the run, not from a\ntemplate. The trees come in all three formats because every program reads a\ndifferent one, and the report describes only the analyses that were carried\nout: if a block was not used, it does not appear.');
    L.push('');
    L.push(es ? 'CONTENIDO' : 'CONTENTS');
    L.push('-'.repeat(48));
    manifest.forEach(m => { L.push(`${m.name}\n    ${m.what}`); });
    L.push('');
    L.push(es ? 'CÓMO VOLVER A ABRIRLO' : 'HOW TO OPEN IT AGAIN');
    L.push('-'.repeat(48));
    L.push(es
      ? 'Los archivos .fasta y .tre se leen en PhylogenyPro (Bloque 2 y Bloque 11),\nen R con ape (read.dna, read.tree, read.nexus) y en cualquier visor de\nárboles. El informe .html se abre con un navegador y se imprime a PDF desde\nahí. Los .svg se abren en Inkscape o Illustrator y no pierden calidad a\nninguna resolución.'
      : 'The .fasta and .tre files are read by PhylogenyPro (Block 2 and Block 11),\nby R with ape (read.dna, read.tree, read.nexus) and by any tree viewer. The\n.html report opens in a browser and prints to PDF from there. The .svg files\nopen in Inkscape or Illustrator and lose nothing at any resolution.');
    return L.join('\n') + '\n';
  }
  const L2on = () => (typeof Prefs !== 'undefined' && Prefs && Prefs.lang ? Prefs.lang !== 'en' : true);

  /* ================================================================
     wiring
     ================================================================ */
  function refresh() {
    const list = trees();
    const has = list.length > 0;
    el('p12NoTree').style.display = has ? 'none' : '';
    ['p12Studio', 'p12Formats', 'p12Figures', 'p12Report'].forEach(id => { el(id).style.display = has ? '' : 'none'; });
    el('p12Zip').style.display = has ? '' : 'none';
    if (!has) return;
    refreshTrees();
    if (!el('p12RepTitle').value) {
      el('p12RepTitle').value = T('Análisis filogenético', 'Phylogenetic analysis');
    }
    drawStudio();
  }

  function init() {
    if (!el('panel-12')) return;
    const redraw = () => drawStudio();
    ['p12Tree', 'p12Layout', 'p12Clado', 'p12SupportAs', 'p12SupportMin', 'p12Width',
     'p12RowHeight', 'p12TipFont', 'p12LineWidth', 'p12Italic', 'p12Scale', 'p12Ladder',
     'p12Title', 'p12AlignTips',
     'p12GroupMethod', 'p12K', 'p12GroupShow', 'p12Legend',
     'p12Images', 'p12ImgSize', 'p12ImgShape', 'p12ImgFrame'].forEach(id => {
      const e2 = el(id);
      if (e2) e2.addEventListener(e2.tagName === 'INPUT' && e2.type === 'text' ? 'input' : 'change', redraw);
    });
    el('p12Root').addEventListener('change', () => {
      el('p12OutWrap').style.display = el('p12Root').value === 'outgroup' ? '' : 'none';
      drawStudio();
    });

    /* pictures: the manager works on whatever the figure is about to draw —
       the tips, or the groups — so the list is never a guess */
    if (el('p12ImgManage')) {
      el('p12ImgManage').addEventListener('click', () => {
        const want = el('p12Images').value;
        let names;
        if (want === 'groups' && B12.groups) {
          names = [];
          for (let g = 0; g < B12.groups.k; g++) names.push(Groups.nameOf(B12.groups, g));
        } else {
          names = (taxa() || []).filter(Boolean);
        }
        OTUImg.openManager({ names, onDone: drawStudio });
      });
    }
    if (el('p12GroupNames')) {
      el('p12GroupNames').addEventListener('click', () => {
        if (!B12.groups) {
          showMessage(el('p12GroupMsg'), 'info',
            L2('Primero corta el árbol en grupos.', 'Cut the tree into groups first.'));
          return;
        }
        const cur = [];
        for (let g = 0; g < B12.groups.k; g++) cur.push(Groups.nameOf(B12.groups, g));
        const ans = prompt(T(
          'Nombres de los grupos, separados por comas:', 'Group names, separated by commas:'), cur.join(', '));
        if (ans == null) return;
        const parts = ans.split(',').map(s => s.trim());
        B12.groupNames = cur.map((d, i) => parts[i] || d);
        drawStudio();
      });
    }
    initStyleControls();
    if (window.OTUImg) OTUImg.onChange(() => { if (el('p12Images') && el('p12Images').value !== 'none') drawStudio(); });
    el('p12Outgroup').addEventListener('change', redraw);
    el('p12AddColour').addEventListener('click', () => {
      const tip = +el('p12ColourClade').value;
      const n = taxa().length;
      if (!B12.colours) B12.colours = new Array(n).fill(0);
      /* colour the clade that tip belongs to, one group at a time */
      const src = shapedTree();
      if (!src) return;
      const F = Tree.flatten(src.tree);
      const below = new Array(F.n);
      for (let i = 0; i < F.post.length; i++) {
        const k = F.post[i];
        below[k] = F.isTip[k] ? [F.tipRow[k]] : [].concat.apply([], F.kids[k].map(c => below[c]));
      }
      /* the smallest clade of two or more tips that contains it */
      let best = null;
      for (let k = 0; k < F.n; k++) {
        if (F.isTip[k]) continue;
        if (below[k].indexOf(tip) < 0) continue;
        if (below[k].length === n) continue;
        if (!best || below[k].length < below[best].length) best = k;
      }
      const group = (Math.max.apply(null, B12.colours) || 0) + 1;
      (best != null ? below[best] : [tip]).forEach(t => { B12.colours[t] = group; });
      drawStudio();
    });
    el('p12ClearColour').addEventListener('click', () => { B12.colours = null; drawStudio(); });

    ['p12Format', 'p12WithLengths', 'p12WithSupport', 'p12AllTrees', 'p12Decimals'].forEach(id => {
      el(id).addEventListener('change', refreshFormats);
    });
    el('p12Download').addEventListener('click', downloadFormat);
    el('p12Copy').addEventListener('click', () => {
      const t = el('p12FormatText');
      t.select();
      try {
        document.execCommand('copy');
        showMessage(el('p12FormatMsg'), 'success', L2('Copiado al portapapeles.', 'Copied to the clipboard.'));
      } catch (e) {
        showMessage(el('p12FormatMsg'), 'info', L2('Selecciona el texto y cópialo con Ctrl+C.',
          'Select the text and copy it with Ctrl+C.'));
      }
    });
    el('p12Dpi').addEventListener('change', refreshFigures);
    el('p12BuildReport').addEventListener('click', buildReport);
    el('p12DownloadReport').addEventListener('click', () => {
      if (!B12.report) buildReport();
      download(B12.report.html, `phylogenypro_informe_${new Date().toISOString().slice(0, 10)}.html`, 'text/html');
    });
    el('p12PrintReport').addEventListener('click', () => {
      if (!B12.report) buildReport();
      const w = window.open('', '_blank');
      if (!w) {
        showMessage(el('p12RepMsg'), 'warn', L2(
          'El navegador bloqueó la ventana. Descarga el informe y ábrelo para imprimirlo.',
          'The browser blocked the window. Download the report and open it to print.'));
        return;
      }
      w.document.write(B12.report.html);
      w.document.close();
      setTimeout(() => w.print(), 600);
    });
    el('p12BuildZip').addEventListener('click', buildZip);
    document.addEventListener('stepchange', e => { if (e.detail.step === 12) refresh(); });
    document.addEventListener('langchange', () => { if (trees().length) refresh(); });
    document.addEventListener('themechange', () => { if (trees().length) drawStudio(); });
  }
  document.addEventListener('DOMContentLoaded', init);

  Object.assign(B12, { refresh, drawStudio, refreshFormats, refreshFigures,
    buildReport, buildZip, formatText, trees, shapedTree, reportTables, tablesForZip });
})();
