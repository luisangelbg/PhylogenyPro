/* PhylogenyPro — molecular phylogenetics in the browser.
   Copyright (C) 2026  Luis Ángel Barrera-Guzmán

   This program is free software: you can redistribute it and/or modify it under
   the terms of the GNU General Public License as published by the Free Software
   Foundation, either version 3 of the License, or (at your option) any later
   version. It is distributed in the hope that it will be useful, but WITHOUT ANY
   WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
   A PARTICULAR PURPOSE. See the GNU General Public License, in the file LICENSE
   at the root of this program, or <https://www.gnu.org/licenses/>.

   ---------------------------------------------------------------------------

   The look of every figure, in one place.

   Every plot in this program draws its colours through CSS custom properties
   and its text through two or three classes — never a literal colour. That was
   done so the figures would follow the light and the dark theme, and it has a
   second consequence worth taking: overriding those properties restyles the
   whole app's figures at once, without a single plotting module knowing.

   So the studio does not edit figures one by one. It writes one rule,

       svg { --c1: …; --text: …; --fig-font: … }

   and every tree, every curve and every heat map in the twelve blocks follows.
   The rule is scoped to `svg` on purpose: the coloured chips of the interface
   use the same tokens, and a figure palette has no business recolouring the
   menus.

   Export is the part that is easy to get wrong. An SVG loaded into an <img>
   has no stylesheet, so the colours have to be written into the file before it
   leaves the page. `js/exportfmt.js` already did that for the theme; here its
   two functions are wrapped so that they resolve to what the studio chose and
   not to what the stylesheet says. What you see is what the PNG carries. */

const FigStyle = {};

(function () {

  /* =====================================================================
     1 · palettes
     ===================================================================== */
  /* Ten colours each, in the order the figures use them. They are chosen to
     stay apart in hue AND in lightness, so a figure still reads in greyscale
     and in the two commonest kinds of colour blindness. */
  const PALETTES = [
    { id: 'phylo', name: ['PhylogenyPro', 'PhylogenyPro'],
      note: ['La del programa: índigo, ámbar y verde de linaje.', 'The program\'s own: indigo, amber and lineage green.'],
      c: ['#2a4a94', '#b76a17', '#227c60', '#8e3b6c', '#2d84b8', '#9a7b12', '#6b4bab', '#1f6f7a', '#a8432f', '#4f6b23'] },
    { id: 'okabe', name: ['Okabe–Ito (segura al daltonismo)', 'Okabe–Ito (colour-blind safe)'],
      note: ['Diseñada para distinguirse con deuteranopía y protanopía. La opción prudente para publicar.',
        'Designed to stay distinguishable under deuteranopia and protanopia. The prudent choice for publishing.'],
      c: ['#0072B2', '#E69F00', '#009E73', '#CC79A7', '#56B4E9', '#F0E442', '#D55E00', '#000000', '#7F7F7F', '#B2DF8A'] },
    { id: 'viridis', name: ['Viridis', 'Viridis'],
      note: ['Perceptualmente uniforme y monótona en luminancia: ordena bien y sobrevive al blanco y negro.',
        'Perceptually uniform and monotonic in lightness: it orders well and survives greyscale.'],
      c: ['#440154', '#472D7B', '#3B518B', '#2C718E', '#21908C', '#27AD81', '#5CC863', '#AADC32', '#D8E219', '#FDE725'] },
    { id: 'earth', name: ['Tierra', 'Earth'],
      note: ['Ocres, verdes y arcillas: para figuras de biogeografía y de suelos.',
        'Ochres, greens and clays: for biogeography and soil figures.'],
      c: ['#6B4423', '#A0713B', '#C9A227', '#7D8C3C', '#4A6B3A', '#2E5A50', '#8C5A3C', '#B5793C', '#5C6B47', '#3F4A3C'] },
    { id: 'cool', name: ['Fríos', 'Cool'],
      note: ['Azules, verdes y violetas. Discreta, buena para muchos grupos.',
        'Blues, greens and violets. Quiet, good for many groups.'],
      c: ['#1F4E79', '#2E75B6', '#41A5C4', '#3FA37A', '#5FB56B', '#7A6CC4', '#4B3F8F', '#2C8A8A', '#6E8FB5', '#39566B'] },
    { id: 'warm', name: ['Cálidos', 'Warm'],
      note: ['Rojos, naranjas y dorados. Llama la atención; úsala cuando la figura sea el argumento.',
        'Reds, oranges and golds. It shouts; use it when the figure is the argument.'],
      c: ['#8C2F1F', '#C1441E', '#E07B24', '#D9A521', '#A85C2E', '#7A3B52', '#B5563F', '#E0A96D', '#96421F', '#C97B2E'] },
    { id: 'pastel', name: ['Pasteles', 'Pastels'],
      note: ['Claros y suaves, para figuras con mucha superficie coloreada y texto encima.',
        'Light and soft, for figures with a lot of coloured area and text on top.'],
      c: ['#9DB4D8', '#EBC08C', '#9ED0B8', '#D9A8C4', '#A8CFE3', '#DCCF95', '#BCAEDC', '#9CC8C8', '#E0A99B', '#BCCB9C'] },
    { id: 'contrast', name: ['Alto contraste', 'High contrast'],
      note: ['Máxima separación entre colores contiguos. Para proyectar o imprimir en papel pobre.',
        'Maximum separation between neighbouring colours. For projecting or printing on poor paper.'],
      c: ['#000000', '#E6194B', '#3CB44B', '#4363D8', '#F58231', '#911EB4', '#008080', '#9A6324', '#800000', '#808000'] },
    { id: 'grey', name: ['Escala de grises', 'Greyscale'],
      note: ['Para revistas que cobran el color. Comprueba aquí si tu figura se entiende sin él.',
        'For journals that charge for colour. Check here whether your figure survives without it.'],
      c: ['#111111', '#3D3D3D', '#5C5C5C', '#777777', '#919191', '#A8A8A8', '#BDBDBD', '#D0D0D0', '#E0E0E0', '#EFEFEF'] },
  ];

  /* =====================================================================
     2 · fonts
     ===================================================================== */
  /* Only families that exist on a normal computer, because a figure that falls
     back to something else on the reviewer's machine is a figure you did not
     design. Nothing is fetched from the network. */
  const FONTS = [
    { id: 'system', name: ['Del sistema', 'System'], css: 'system-ui, sans-serif' },
    { id: 'humanist', name: ['Sans humanista', 'Humanist sans'], css: '"Segoe UI", "Helvetica Neue", Arial, sans-serif' },
    { id: 'grotesque', name: ['Sans neogrotesca', 'Neo-grotesque sans'], css: 'Helvetica, Arial, sans-serif' },
    { id: 'serif', name: ['Serif', 'Serif'], css: 'Georgia, "Times New Roman", Times, serif' },
    { id: 'slab', name: ['Serif de trazo grueso', 'Slab serif'], css: '"Bookman Old Style", "Palatino Linotype", Palatino, Georgia, serif' },
    { id: 'condensed', name: ['Estrecha', 'Condensed'], css: '"Arial Narrow", "Liberation Sans Narrow", "Segoe UI", sans-serif' },
    { id: 'mono', name: ['Monoespaciada', 'Monospaced'], css: 'ui-monospace, "Cascadia Code", Consolas, "Courier New", monospace' },
  ];

  const BACKGROUNDS = [
    { id: 'none', name: ['Sin fondo (transparente)', 'No background (transparent)'], css: 'transparent' },
    { id: 'card', name: ['El del tema', 'The theme\'s'], css: 'var(--card-bg)' },
    { id: 'white', name: ['Blanco', 'White'], css: '#ffffff' },
    { id: 'paper', name: ['Papel crema', 'Cream paper'], css: '#fbf8f1' },
  ];

  /* the tokens a figure may recolour; the rest of the theme is left alone */
  const TOKENS = ['text', 'text-muted', 'ink', 'border', 'border-strong', 'bg-soft', 'card-bg',
    'accent', 'primary', 'gold', 'leaf', 'sky', 'rose'];

  /* =====================================================================
     3 · the state and the rule it writes
     ===================================================================== */
  const DEF = {
    palette: 'phylo',
    custom: {},          // c1…c10 the user overrode by hand
    tokens: {},          // any of TOKENS the user overrode by hand
    font: 'system',
    fontScale: 1,        // multiplies every font-size in the figure
    lineScale: 1,        // multiplies every stroke-width
    background: 'none',
    grid: true,
  };
  let S = Object.assign({}, DEF, { custom: {}, tokens: {} });
  const listeners = new Set();

  function paletteOf(id) { return PALETTES.find(p => p.id === (id || S.palette)) || PALETTES[0]; }
  function fontOf(id) { return FONTS.find(f => f.id === (id || S.font)) || FONTS[0]; }
  function backgroundOf(id) { return BACKGROUNDS.find(b => b.id === (id || S.background)) || BACKGROUNDS[0]; }

  /* the ten colours as they stand, palette plus whatever was edited by hand */
  function colours() {
    const base = paletteOf().c.slice();
    for (let i = 1; i <= 10; i++) if (S.custom['c' + i]) base[i - 1] = S.custom['c' + i];
    return base;
  }

  /* everything the figures should resolve to: what `Export.inlineVars` needs */
  function vars() {
    const out = {};
    colours().forEach((c, i) => { out['c' + (i + 1)] = c; });
    Object.keys(S.tokens).forEach(k => { if (S.tokens[k]) out[k] = S.tokens[k]; });
    return out;
  }

  const STYLE_ID = 'figstyle-rule';

  function cssText() {
    const v = vars();
    const decl = Object.keys(v).map(k => `--${k}:${v[k]}`).join(';');
    const bg = backgroundOf().css;
    return `svg{${decl};--fig-font:${fontOf().css};` +
      `font-size:${S.fontScale === 1 ? '' : (100 * S.fontScale).toFixed(1) + '%'};` +
      (bg === 'transparent' ? '' : `background:${bg};`) + '}' +
      (S.fontScale === 1 ? '' : `svg text{font-size:inherit}`) +
      (S.lineScale === 1 ? '' : `svg line,svg path,svg circle,svg rect,svg polyline,svg polygon{stroke-width:${S.lineScale}px}`);
  }

  function apply() {
    if (typeof document === 'undefined') return;
    let st = document.getElementById(STYLE_ID);
    if (!st) {
      st = document.createElement('style');
      st.id = STYLE_ID;
      document.head.appendChild(st);
    }
    st.textContent = cssText();
    listeners.forEach(fn => { try { fn(); } catch (e) { console.error(e); } });
  }

  function set(patch) {
    if (!patch) return;
    if (patch.custom) { S.custom = Object.assign({}, S.custom, patch.custom); delete patch.custom; }
    if (patch.tokens) { S.tokens = Object.assign({}, S.tokens, patch.tokens); delete patch.tokens; }
    S = Object.assign(S, patch);
    apply();
  }
  function reset() { S = Object.assign({}, DEF, { custom: {}, tokens: {} }); apply(); }
  function get() { return JSON.parse(JSON.stringify(S)); }
  function onChange(fn) { listeners.add(fn); return () => listeners.delete(fn); }

  /* =====================================================================
     4 · making the export carry the same look
     ===================================================================== */
  /* `Export.readVars` reads the computed theme off the document root, which
     knows nothing of a rule scoped to `svg`. Wrapping it is enough: the studio's
     colours win, everything it did not touch keeps coming from the theme.
     And the font, which is a class and not a variable, is written into the
     <style> that `inlineVars` already injects into every exported SVG. */
  function patchExport() {
    if (typeof window === 'undefined' || !window.Export) return false;
    if (Export.readVars && !Export.readVars.__styled) {
      const orig = Export.readVars;
      const wrapped = function () { return Object.assign(orig.apply(this, arguments), vars()); };
      wrapped.__styled = true;
      Export.readVars = wrapped;
    }
    if (Export.inlineVars && !Export.inlineVars.__styled) {
      const orig = Export.inlineVars;
      const wrapped = function (svgText, v) {
        let out = orig.call(this, svgText, Object.assign({}, v, vars()));
        const fam = fontOf().css;
        /* the injected block hardcodes system-ui; make it say what was chosen,
           and give .art-ink the fill it never had in the stylesheet */
        out = out.split('font-family:system-ui,sans-serif').join('font-family:' + fam);
        if (out.indexOf('.art-ink{') < 0) {
          const ink = (v && v.ink) || '#161b2a';
          out = out.replace('</style>', `.art-ink{fill:${ink};font-family:${fam}}</style>`);
        }
        const bg = backgroundOf();
        if (bg.id !== 'none' && out.indexOf('data-figbg') < 0) {
          const css = bg.css.indexOf('var(') === 0 ? ((v && v['card-bg']) || '#ffffff') : bg.css;
          out = out.replace(/(<svg[^>]*>)/, `$1<rect data-figbg="1" x="0" y="0" width="100%" height="100%" fill="${css}"/>`);
        }
        return out;
      };
      wrapped.__styled = true;
      Export.inlineVars = wrapped;
    }
    return true;
  }

  Object.assign(FigStyle, {
    PALETTES, FONTS, BACKGROUNDS, TOKENS, DEF,
    paletteOf, fontOf, backgroundOf, colours, vars, cssText,
    apply, set, reset, get, onChange, patchExport,
  });
  if (typeof window !== 'undefined') window.FigStyle = FigStyle;

  /* Patch as soon as Export is there, which is while the page is still being
     parsed if its script came first, and again when the document is ready in
     case it came later. Waiting only for DOMContentLoaded would leave anything
     that exports during parsing — the test page does exactly that — resolving
     its colours against the stylesheet instead of against the studio. */
  if (typeof document !== 'undefined') {
    apply();
    if (!patchExport()) {
      document.addEventListener('DOMContentLoaded', () => { apply(); patchExport(); });
    } else if (document.readyState === 'loading') {
      document.addEventListener('DOMContentLoaded', () => { apply(); patchExport(); });
    }
  }

})();
