/* PhylogenyPro — molecular phylogenetics in the browser.
   Copyright (C) 2026  Luis Ángel Barrera-Guzmán

   This program is free software: you can redistribute it and/or modify it under
   the terms of the GNU General Public License as published by the Free Software
   Foundation, either version 3 of the License, or (at your option) any later
   version. It is distributed in the hope that it will be useful, but WITHOUT ANY
   WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
   A PARTICULAR PURPOSE. See the GNU General Public License, in the file LICENSE
   at the root of this program, or <https://www.gnu.org/licenses/>. */

/* PhylogenyPro — figures of Block 3: composition, saturation and model fit.
   Every figure is a real SVG built here, with theme colours, so it exports at
   any resolution and follows light and dark mode. */

const Plots3 = {};

(function () {

  const V = n => `var(--${n})`;
  const f1 = v => (+v).toFixed(1);
  const svg = (vb, inner, plot) => `<svg viewBox="${vb}"${plot ? ` data-plot="${plot}"` : ''} xmlns="http://www.w3.org/2000/svg">${inner}</svg>`;
  /* marks for the figure studio: the area of the data (data-plot) and the
     legend, in one group (data-role="legend") with its entries numbered (data-li) */
  const plotArea = (x, y, w, h) => [x, y, w, h].map(v => +(+v).toFixed(2)).join(' ');
  const line = (x1, y1, x2, y2, st, w, ex) => `<line x1="${f1(x1)}" y1="${f1(y1)}" x2="${f1(x2)}" y2="${f1(y2)}" stroke="${st}" stroke-width="${w || 1}" ${ex || ''}/>`;
  const circ = (cx, cy, r, fill, ex) => `<circle cx="${f1(cx)}" cy="${f1(cy)}" r="${f1(r)}" fill="${fill}" ${ex || ''}/>`;
  const rect = (x, y, w, h, fill, ex) => `<rect x="${f1(x)}" y="${f1(y)}" width="${f1(Math.max(0, w))}" height="${f1(Math.max(0, h))}" fill="${fill}" ${ex || ''}/>`;
  const pathEl = (d, st, w, ex) => `<path d="${d}" stroke="${st}" stroke-width="${w || 1}" fill="none" ${ex || ''}/>`;
  const poly = pts => pts.map((p, i) => (i ? 'L' : 'M') + f1(p[0]) + ' ' + f1(p[1])).join(' ');
  function txt(x, y, s, cls, size, anchor, ex) {
    return `<text x="${f1(x)}" y="${f1(y)}" class="${cls || 'art-mut'}" font-size="${size || 9}" text-anchor="${anchor || 'start'}" ${ex || ''}>${esc(s)}</text>`;
  }

  /* ================================================================
     saturation: transitions and transversions against the corrected distance
     ================================================================ */
  function saturation(sat, opts) {
    opts = opts || {};
    const W = 620, H = 380, mx = 62, my = 52;
    const pts = sat.points.filter(p => p.d != null);
    const maxD = Math.max(0.05, ...pts.map(p => p.d)) * 1.05;
    const maxY = Math.max(0.05, ...pts.map(p => Math.max(p.ts, p.tv))) * 1.1;
    const X = d => mx + (W - mx - 20) * d / maxD;
    const Y = v => (H - my) - (H - my - 26) * v / maxY;
    let s = '';
    /* axes */
    s += line(mx, H - my, W - 20, H - my, V('border-strong'), 1.4);
    s += line(mx, 20, mx, H - my, V('border-strong'), 1.4);
    for (let i = 0; i <= 5; i++) {
      const d = maxD * i / 5;
      s += line(X(d), H - my, X(d), H - my + 5, V('border-strong'), 1.2);
      s += txt(X(d), H - my + 18, d.toFixed(2), 'art-mut', 9, 'middle');
      const v = maxY * i / 5;
      s += line(mx - 5, Y(v), mx, Y(v), V('border-strong'), 1.2);
      s += txt(mx - 8, Y(v) + 3, v.toFixed(2), 'art-mut', 9, 'end');
    }
    /* the fitted lines through the origin */
    s += line(X(0), Y(0), X(maxD), Y(Math.min(maxY, sat.slopeTs * maxD)), V('accent'), 2, 'opacity="0.75"');
    s += line(X(0), Y(0), X(maxD), Y(Math.min(maxY, sat.slopeTv * maxD)), V('sky'), 2, 'opacity="0.75"');
    pts.forEach(p => {
      s += circ(X(p.d), Y(p.ts), 2.6, V('accent'), 'opacity="0.55"');
      s += circ(X(p.d), Y(p.tv), 2.6, V('sky'), 'opacity="0.55"');
    });
    s += txt(W / 2, H - 16, T('distancia corregida (sustituciones por sitio)', 'corrected distance (substitutions per site)'), 'art-mut', 10, 'middle');
    s += txt(-H / 2, 16, T('proporción observada', 'observed proportion'), 'art-mut', 10, 'middle', 'transform="rotate(-90)"');
    s += '<g data-role="legend">';
    s += circ(W - 150, 28, 4, V('accent'), 'data-li="0"'); s += txt(W - 142, 31, T('transiciones', 'transitions'), 'art-mut', 9.5, 'start', 'data-li="0"');
    s += circ(W - 150, 44, 4, V('sky'), 'data-li="1"'); s += txt(W - 142, 47, T('transversiones', 'transversions'), 'art-mut', 9.5, 'start', 'data-li="1"');
    s += '</g>';
    return svg(`0 0 ${W} ${H}`, s, plotArea(mx, 20, W - mx - 20, H - my - 20));
  }

  /* ================================================================
     base composition per sequence
     ================================================================ */
  function composition(comp, taxa, opts) {
    opts = opts || {};
    const n = comp.rows.length;
    const alpha = comp.alphabet;
    const W = 620, H = Math.max(220, 34 + n * 18 + 46), mx = 150, my = 34;
    const barW = W - mx - 30;
    const colours = alpha.length === 4 ? ['base-a', 'base-c', 'base-g', 'base-t'] : alpha.map((a, i) => 'c' + (1 + i % 10));
    let s = '';
    comp.rows.forEach((r, i) => {
      const y = my + i * 18;
      s += txt(mx - 8, y + 11, taxa[i], 'art-mut', 9.5, 'end');
      let x = mx;
      alpha.forEach((a, k) => {
        const w = barW * r.freqs[k];
        s += rect(x, y + 2, w, 13, V(colours[k]), 'opacity="0.85"');
        x += w;
      });
    });
    /* the overall frequencies as a reference line */
    let x = mx;
    comp.grand.forEach((g, k) => {
      x += barW * g;
      if (k < alpha.length - 1) s += line(x, my - 4, x, my + n * 18 + 2, V('text'), 1, 'stroke-dasharray="3 3" opacity="0.55"');
    });
    s += txt(mx, 18, T('cada barra es una secuencia; la línea punteada marca la composición media', 'each bar is one sequence; the dashed line marks the overall composition'), 'art-mut', 9.5);
    let lx = mx;
    s += '<g data-role="legend">';
    alpha.slice(0, 8).forEach((a, k) => {
      s += rect(lx, H - 22, 11, 11, V(colours[k]), `opacity="0.85" data-li="${k}"`);
      s += txt(lx + 15, H - 13, a, 'art-mut', 9.5, 'start', `data-li="${k}"`);
      lx += 34;
    });
    s += '</g>';
    return svg(`0 0 ${W} ${H}`, s, plotArea(mx, my, barW, n * 18));
  }

  /* ================================================================
     model comparison: ΔAICc and the Akaike weights
     ================================================================ */
  function modelBars(fits, key, opts) {
    opts = opts || {};
    const top = fits.slice(0, opts.top || 12);
    const W = 620, H = Math.max(200, 40 + top.length * 22 + 30), mx = 110;
    const maxD = Math.max(1, ...top.map(f => f['d' + key]));
    let s = txt(mx, 20, T(`Δ${key} de los mejores modelos (0 = el mejor)`, `Δ${key} of the best models (0 = the best)`), 'art-mut', 10);
    top.forEach((f, i) => {
      const y = 32 + i * 22;
      const w = (W - mx - 90) * Math.min(1, f['d' + key] / maxD);
      s += txt(mx - 8, y + 12, f.name, 'art-txt', 10, 'end', i === 0 ? 'font-weight="700"' : '');
      s += rect(mx, y + 2, Math.max(2, w), 14, i === 0 ? V('leaf') : V('primary'), 'rx="3" opacity="0.85"');
      s += txt(mx + Math.max(2, w) + 6, y + 13, f['d' + key].toFixed(2) + (f['w' + key] != null ? `   ${(f['w' + key] * 100).toFixed(1)}%` : ''), 'art-mut', 9);
    });
    /* the conventional threshold of 2 units */
    const x2 = mx + (W - mx - 90) * Math.min(1, 2 / maxD);
    s += line(x2, 28, x2, 32 + top.length * 22, V('accent'), 1.2, 'stroke-dasharray="4 3"');
    s += txt(x2 + 3, 26, 'Δ = 2', 'art-mut', 8.5, 'start', `fill="${V('accent')}"`);
    return svg(`0 0 ${W} ${H}`, s);
  }

  /* ================================================================
     the fitted rate matrix, drawn as a grid
     ================================================================ */
  function rateMatrix(M, alphabet, opts) {
    opts = opts || {};
    const S = alphabet.length;
    const cell = S <= 4 ? 52 : Math.max(14, Math.min(26, 420 / S));
    const W = Math.max(260, cell * S + 120), H = cell * S + 96;
    const x0 = 70, y0 = 46;
    let maxOff = 0;
    for (let i = 0; i < S; i++) for (let j = 0; j < S; j++) if (i !== j && M.Q[i][j] > maxOff) maxOff = M.Q[i][j];
    let s = txt(x0, 24, T('matriz de tasas Q (filas: desde · columnas: hacia)', 'rate matrix Q (rows: from · columns: to)'), 'art-mut', 10);
    for (let i = 0; i < S; i++) {
      s += txt(x0 - 8, y0 + i * cell + cell * 0.62, alphabet[i], 'art-txt', Math.min(12, cell * 0.4), 'end');
      s += txt(x0 + i * cell + cell / 2, y0 - 6, alphabet[i], 'art-txt', Math.min(12, cell * 0.4), 'middle');
      for (let j = 0; j < S; j++) {
        const v = i === j ? 0 : M.Q[i][j] / (maxOff || 1);
        s += rect(x0 + j * cell, y0 + i * cell, cell - 1.5, cell - 1.5, i === j ? V('bg-soft') : V('primary'),
          `opacity="${i === j ? 1 : (0.12 + 0.8 * v).toFixed(2)}" rx="2"`);
        if (cell >= 34 && i !== j) s += txt(x0 + j * cell + cell / 2, y0 + i * cell + cell * 0.6, M.Q[i][j].toFixed(2), 'art-txt', 10, 'middle',
          `fill="${v > 0.45 ? 'white' : V('text')}"`);
      }
    }
    /* equilibrium frequencies */
    s += txt(x0, y0 + S * cell + 24, T('frecuencias de equilibrio: ', 'equilibrium frequencies: ') +
      alphabet.map((a, i) => `${a} ${M.pi[i].toFixed(3)}`).join('   '), 'art-mut', 9.5);
    return svg(`0 0 ${W} ${H}`, s);
  }

  /* ================================================================
     the distribution of rates among sites, under +G and +I
     ================================================================ */
  function rateCategories(M, opts) {
    const W = 620, H = 240, mx = 56, my = 44;
    const cats = M.catRates, w = M.catWeights;
    const maxR = Math.max(1.5, ...cats) * 1.1;
    const X = r => mx + (W - mx - 24) * r / maxR;
    let s = line(mx, H - my, W - 24, H - my, V('border-strong'), 1.4);
    s += txt(mx, 22, T('tasa relativa de cada categoría de sitios', 'relative rate of each category of sites'), 'art-mut', 10);
    const barH = 30;
    if (M.pInv > 0) {
      s += rect(mx, H - my - barH - 2, Math.max(3, X(0) - mx + 3), barH, V('border-strong'), 'rx="3"');
      s += txt(mx + 6, H - my - barH - 8, `${T('invariables', 'invariable')} ${(M.pInv * 100).toFixed(1)}%`, 'art-mut', 9);
    }
    cats.forEach((r, i) => {
      const x = X(r);
      s += line(x, H - my, x, H - my - 8 - 120 * w[i] * cats.length / 2, V('primary'), 8, 'stroke-linecap="round" opacity="0.85"');
      s += txt(x, H - my + 16, r.toFixed(3), 'art-mut', 9, 'middle');
    });
    s += txt(W / 2, H - 10, T('tasa relativa (1 = la media)', 'relative rate (1 = the mean)'), 'art-mut', 9.5, 'middle');
    s += line(X(1), H - my, X(1), 40, V('accent'), 1.2, 'stroke-dasharray="4 3"');
    s += txt(X(1) + 4, 46, T('media', 'mean'), 'art-mut', 9, 'start', `fill="${V('accent')}"`);
    return svg(`0 0 ${W} ${H}`, s);
  }

  Object.assign(Plots3, { saturation, composition, modelBars, rateMatrix, rateCategories });
  window.Plots3 = Plots3;
})();
