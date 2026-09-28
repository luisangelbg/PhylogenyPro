/* PhylogenyPro — molecular phylogenetics in the browser.
   Copyright (C) 2026  Luis Ángel Barrera-Guzmán

   This program is free software: you can redistribute it and/or modify it under
   the terms of the GNU General Public License as published by the Free Software
   Foundation, either version 3 of the License, or (at your option) any later
   version. It is distributed in the hope that it will be useful, but WITHOUT ANY
   WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
   A PARTICULAR PURPOSE. See the GNU General Public License, in the file LICENSE
   at the root of this program, or <https://www.gnu.org/licenses/>. */

/* PhylogenyPro — figures of Block 7: the regression and the chronogram.

   The chronogram is the figure the whole app has been walking towards, and it
   is not a phylogram with different numbers: the horizontal axis is time, it
   runs backwards from the present, the geological periods are drawn behind the
   tree because a date without its epoch means little to a reader, and every
   node carries the interval its estimate really has. A chronogram without those
   bars claims a precision that no dating method has. */

const Plots7 = {};

(function () {

  const V = n => `var(--${n})`;
  const f1 = v => (+v).toFixed(1);
  const svg = (vb, inner) => `<svg viewBox="${vb}" xmlns="http://www.w3.org/2000/svg">${inner}</svg>`;
  const line = (x1, y1, x2, y2, st, w, ex) => `<line x1="${f1(x1)}" y1="${f1(y1)}" x2="${f1(x2)}" y2="${f1(y2)}" stroke="${st}" stroke-width="${w || 1}" ${ex || ''}/>`;
  const rect = (x, y, w, h, fill, ex) => `<rect x="${f1(x)}" y="${f1(y)}" width="${f1(Math.max(0, w))}" height="${f1(Math.max(0, h))}" fill="${fill}" ${ex || ''}/>`;
  const circ = (cx, cy, r, fill, ex) => `<circle cx="${f1(cx)}" cy="${f1(cy)}" r="${f1(r)}" fill="${fill}" ${ex || ''}/>`;
  function txt(x, y, s, cls, size, anchor, ex) {
    return `<text x="${f1(x)}" y="${f1(y)}" class="${cls || 'art-mut'}" font-size="${size || 9}" text-anchor="${anchor || 'start'}" ${ex || ''}>${esc(s)}</text>`;
  }
  function niceStep(span, want) {
    const raw = span / Math.max(1, want);
    const mag = Math.pow(10, Math.floor(Math.log10(Math.max(1e-12, raw))));
    const n = raw / mag;
    return (n <= 1 ? 1 : n <= 2 ? 2 : n <= 5 ? 5 : 10) * mag;
  }
  const fmtTick = (v, step) => (+v).toFixed(Math.max(0, Math.ceil(-Math.log10(step))));

  /* ================================================================
     root-to-tip regression
     ================================================================ */
  function regression(fit, opts) {
    opts = opts || {};
    const W = opts.width || 660, H = opts.height || 340;
    const mx = 72, my = 52, top = 18, right = 20;
    if (!fit || !fit.x || fit.x.length < 2) return '';
    const xs = fit.x, ys = fit.y;
    let x0 = Math.min.apply(null, xs), x1 = Math.max.apply(null, xs);
    let y0 = Math.min.apply(null, ys), y1 = Math.max.apply(null, ys);
    const padX = (x1 - x0) * 0.08 || 1, padY = (y1 - y0) * 0.12 || 1e-3;
    x0 -= padX; x1 += padX; y0 -= padY; y1 += padY;
    const X = v => mx + (W - mx - right) * (v - x0) / (x1 - x0);
    const Y = v => (H - my) - (H - my - top) * (v - y0) / (y1 - y0);

    let s = '';
    s += line(mx, H - my, W - right, H - my, V('border-strong'), 1.4);
    s += line(mx, top, mx, H - my, V('border-strong'), 1.4);
    const sx = niceStep(x1 - x0, 6), sy = niceStep(y1 - y0, 5);
    for (let v = Math.ceil(x0 / sx) * sx; v <= x1; v += sx) {
      s += line(X(v), H - my, X(v), H - my + 5, V('border-strong'), 1.1);
      s += txt(X(v), H - my + 17, fmtTick(v, sx), 'art-mut', 9, 'middle');
    }
    for (let v = Math.ceil(y0 / sy) * sy; v <= y1; v += sy) {
      s += line(mx - 5, Y(v), W - right, Y(v), V('border'), 0.7, 'opacity="0.4"');
      s += txt(mx - 8, Y(v) + 3, fmtTick(v, sy), 'art-mut', 9, 'end');
    }
    /* the fitted line */
    s += line(X(x0), Y(fit.intercept + fit.slope * x0), X(x1), Y(fit.intercept + fit.slope * x1),
      V('accent'), 2, 'opacity="0.85"');
    xs.forEach((v, i) => {
      s += circ(X(v), Y(ys[i]), 3.2, V('sky'), 'opacity="0.8"');
    });
    s += txt(mx + (W - mx - right) / 2, H - 12, opts.xLabel || T('fecha de la punta', 'tip date'), 'art-mut', 10, 'middle');
    s += txt(-(top + (H - my - top) / 2), 15,
      opts.yLabel || T('distancia raíz–punta (sustituciones/sitio)', 'root-to-tip distance (substitutions/site)'),
      'art-mut', 10, 'middle', 'transform="rotate(-90)"');
    const lines = [
      `r² = ${fit.r2.toFixed(4)}`,
      `${T('tasa', 'rate')} = ${fit.rate.toExponential(3)}`,
      fit.rootDate != null ? `${T('raíz en', 'root at')} ${fit.rootDate.toFixed(2)}` : '',
    ].filter(Boolean);
    lines.forEach((l, i) => { s += txt(W - right - 8, top + 12 + i * 13, l, 'art-mut', 9.5, 'end'); });
    return svg(`0 0 ${W} ${H}`, s);
  }

  /* ================================================================
     the chronogram
     ================================================================ */
  /* tree: nodes carry `age`; opts.intervals: Map from a split key to
     {lower, upper} for the uncertainty bars. */
  function chronogram(tree, opts) {
    opts = opts || {};
    const labels = opts.labels || [];
    const tips = Tree.tips(tree);
    const nT = tips.length;
    const rowH = opts.rowHeight || Math.max(13, Math.min(26, 520 / Math.max(1, nT)));
    const W = opts.width || 820;
    const geoH = opts.geo ? 34 : 0;
    const top = 16, bottom = 46 + geoH, left = 16;
    const labelW = opts.labelWidth || Math.min(260, 9 + Math.max.apply(null, tips.map(t =>
      String(labels[t.tip] != null ? labels[t.tip] : (t.label || '')).length)) * 6.2);
    const H = top + nT * rowH + bottom;
    const plotW = W - left - labelW - 14;

    /* every node's age; tips at their own age (0 unless dated) */
    const F = Tree.flatten(tree);
    const age = new Float64Array(F.n);
    for (let k = 0; k < F.n; k++) age[k] = F.nodes[k].age != null ? F.nodes[k].age : 0;
    let maxAge = 0;
    for (let k = 0; k < F.n; k++) maxAge = Math.max(maxAge, age[k]);
    if (opts.intervals) opts.intervals.forEach(v => { if (v && v.upper != null) maxAge = Math.max(maxAge, v.upper); });
    if (!(maxAge > 0)) maxAge = 1;
    const X = a => left + labelW * 0 + plotW * (1 - a / maxAge) + 0;   // present at the right
    const y = new Map();
    let row = 0;
    (function assign(node) {
      if (!node.children.length) { y.set(node, top + rowH * (row + 0.5)); row++; return; }
      node.children.forEach(c => assign(c.node));
      const ys = node.children.map(c => y.get(c.node));
      y.set(node, (Math.min.apply(null, ys) + Math.max.apply(null, ys)) / 2);
    })(tree);

    let s = '';
    /* the geological band */
    if (opts.geo && typeof GeoTime !== 'undefined') {
      const bandY = H - bottom + 6;
      ['period', 'epoch'].forEach((rank, ri) => {
        const units = GeoTime.within(rank, 0, maxAge);
        units.forEach(u => {
          const xa = X(Math.min(maxAge, u.base)), xb = X(Math.max(0, u.top));
          const w = xb - xa;
          if (w < 1) return;
          s += rect(xa, bandY + ri * 16, w, 15, u.color, 'opacity="0.55" stroke="var(--border)" stroke-width="0.4"');
          if (w > 26) s += txt((xa + xb) / 2, bandY + ri * 16 + 11, w > 64 ? u.name : u.name.slice(0, 2), 'art-mut', 8, 'middle');
        });
      });
    }
    /* the time axis */
    const axisY = H - bottom - 2;
    s += line(X(maxAge), axisY, X(0), axisY, V('border-strong'), 1.3);
    const st = niceStep(maxAge, 7);
    for (let a = 0; a <= maxAge + 1e-9; a += st) {
      s += line(X(a), axisY, X(a), axisY + 5, V('border-strong'), 1.1);
      s += txt(X(a), axisY + 17, fmtTick(a, st), 'art-mut', 9, 'middle');
      s += line(X(a), top, X(a), axisY, V('border'), 0.6, 'opacity="0.35"');
    }
    s += txt((X(maxAge) + X(0)) / 2, H - 6 - geoH, opts.unit || 'Ma', 'art-mut', 10, 'middle');

    /* the branches */
    (function draw(node) {
      const k = F.nodes.indexOf(node);
      const xn = X(age[k]);
      node.children.forEach(c => {
        const kc = F.nodes.indexOf(c.node);
        const xc = X(age[kc]);
        s += line(xn, y.get(c.node), xc, y.get(c.node), V('ink'), 1.4);
        s += line(xn, y.get(node), xn, y.get(c.node), V('ink'), 1.4);
        draw(c.node);
      });
    })(tree);

    /* the uncertainty bars */
    if (opts.intervals && opts.bars !== false) {
      const splits = Tree.splits(tree, nT);
      splits.forEach((node, key) => {
        const iv = opts.intervals.get(key);
        if (!iv || iv.lower == null) return;
        const yy = y.get(node);
        s += line(X(iv.upper), yy, X(iv.lower), yy, V('accent'), 5.5, 'opacity="0.35" stroke-linecap="round"');
      });
    }
    /* The ages. Two nodes of similar age and similar height print their numbers
       on top of each other, and a chronogram of any size has plenty of those.
       The most inclusive clades keep their age — those are the dates a reader
       is looking for — and the rest are dropped and counted, rather than left
       to smear into an unreadable blur. */
    let hiddenAges = 0;
    if (opts.showAges) {
      const fs = Math.max(7.5, rowH * 0.42);
      const marks = [];
      Tree.nodes(tree).forEach(node => {
        if (!node.children.length) return;
        const k = F.nodes.indexOf(node);
        const lab = (+age[k]).toFixed(opts.ageDecimals == null ? 1 : opts.ageDecimals);
        marks.push({ node, lab, x: X(age[k]) - 3, y: y.get(node) - 4, n: Tree.nTips(node) });
      });
      const boxes = marks.map(m => ({
        x: m.x - m.lab.length * fs * 0.56, y: m.y - fs,
        w: m.lab.length * fs * 0.56 + 3, h: fs + 2, rank: -m.n,
      }));
      const keep = (window.TreeView && opts.ageDeclutter !== false)
        ? TreeView.dropColliding(boxes) : boxes.map(() => true);
      marks.forEach((m, i) => {
        if (!keep[i]) { hiddenAges++; return; }
        s += txt(m.x, m.y, m.lab, 'art-mut', fs, 'end');
      });
    }
    Plots7.lastInfo = { hiddenAges };
    /* the tips */
    tips.forEach(t => {
      const k = F.nodes.indexOf(t);
      const name = labels[t.tip] != null ? labels[t.tip] : (t.label || '');
      s += circ(X(age[k]), y.get(t), 2.2, V('accent'));
      s += txt(X(age[k]) + 6, y.get(t) + rowH * 0.3, name, 'art-ink', Math.max(8, rowH * 0.6));
    });
    return svg(`0 0 ${W} ${H}`, s);
  }

  /* the cross-validation curve for the smoothing parameter */
  function smoothing(points, opts) {
    opts = opts || {};
    const W = opts.width || 520, H = opts.height || 220;
    const mx = 62, my = 44, top = 16, right = 16;
    if (!points || points.length < 2) return '';
    const lx = points.map(p => Math.log10(p.lambda));
    const ly = points.map(p => Math.log10(Math.max(1e-12, p.score)));
    const x0 = Math.min.apply(null, lx), x1 = Math.max.apply(null, lx);
    const y0 = Math.min.apply(null, ly), y1 = Math.max.apply(null, ly);
    const X = v => mx + (W - mx - right) * (v - x0) / Math.max(1e-9, x1 - x0);
    const Y = v => (H - my) - (H - my - top) * (v - y0) / Math.max(1e-9, y1 - y0);
    let s = '';
    s += line(mx, H - my, W - right, H - my, V('border-strong'), 1.3);
    s += line(mx, top, mx, H - my, V('border-strong'), 1.3);
    let d = '';
    lx.forEach((v, i) => { d += (i ? 'L' : 'M') + f1(X(v)) + ' ' + f1(Y(ly[i])); });
    s += `<path d="${d}" stroke="${V('accent')}" stroke-width="1.8" fill="none"/>`;
    let bi = 0;
    points.forEach((p, i) => { if (p.score < points[bi].score) bi = i; });
    points.forEach((p, i) => { s += circ(X(lx[i]), Y(ly[i]), i === bi ? 4.2 : 2.6, i === bi ? V('danger') : V('sky')); });
    lx.forEach((v, i) => { if (i % Math.ceil(lx.length / 6) === 0) s += txt(X(v), H - my + 16, points[i].lambda.toString(), 'art-mut', 8.5, 'middle'); });
    s += txt(mx + (W - mx - right) / 2, H - 8, T('suavizado λ (escala logarítmica)', 'smoothing λ (log scale)'), 'art-mut', 10, 'middle');
    s += txt(-(top + (H - my - top) / 2), 14, T('error de validación cruzada', 'cross-validation error'), 'art-mut', 10, 'middle', 'transform="rotate(-90)"');
    return svg(`0 0 ${W} ${H}`, s);
  }

  Object.assign(Plots7, { regression, chronogram, smoothing });
  window.Plots7 = Plots7;
})();
