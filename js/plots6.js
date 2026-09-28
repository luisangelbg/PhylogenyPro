/* PhylogenyPro — molecular phylogenetics in the browser.
   Copyright (C) 2026  Luis Ángel Barrera-Guzmán

   This program is free software: you can redistribute it and/or modify it under
   the terms of the GNU General Public License as published by the Free Software
   Foundation, either version 3 of the License, or (at your option) any later
   version. It is distributed in the hope that it will be useful, but WITHOUT ANY
   WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
   A PARTICULAR PURPOSE. See the GNU General Public License, in the file LICENSE
   at the root of this program, or <https://www.gnu.org/licenses/>. */

/* PhylogenyPro — figures of Block 6: the pictures a Bayesian run is judged by.

   A trace and a histogram are not decoration here. The trace is how a reader
   sees that the chain stopped climbing and started wandering — the caterpillar
   MrBayes users look for — and that two independent runs wander over the same
   band. The histogram is the posterior itself, with its highest-density
   interval marked, which is the estimate. Both are real SVG built here, with
   the theme's colours, so they export at any resolution. */

const Plots6 = {};

(function () {

  const V = n => `var(--${n})`;
  const f1 = v => (+v).toFixed(1);
  const svg = (vb, inner) => `<svg viewBox="${vb}" xmlns="http://www.w3.org/2000/svg">${inner}</svg>`;
  const line = (x1, y1, x2, y2, st, w, ex) => `<line x1="${f1(x1)}" y1="${f1(y1)}" x2="${f1(x2)}" y2="${f1(y2)}" stroke="${st}" stroke-width="${w || 1}" ${ex || ''}/>`;
  const rect = (x, y, w, h, fill, ex) => `<rect x="${f1(x)}" y="${f1(y)}" width="${f1(Math.max(0, w))}" height="${f1(Math.max(0, h))}" fill="${fill}" ${ex || ''}/>`;
  const pathEl = (d, st, w, ex) => `<path d="${d}" stroke="${st}" stroke-width="${w || 1}" fill="none" ${ex || ''}/>`;
  const poly = pts => pts.map((p, i) => (i ? 'L' : 'M') + f1(p[0]) + ' ' + f1(p[1])).join(' ');
  function txt(x, y, s, cls, size, anchor, ex) {
    return `<text x="${f1(x)}" y="${f1(y)}" class="${cls || 'art-mut'}" font-size="${size || 9}" text-anchor="${anchor || 'start'}" ${ex || ''}>${esc(s)}</text>`;
  }
  const COLOURS = ['accent', 'sky', 'success', 'warning', 'danger', 'violet', 'teal', 'amber'];

  /* a tick step that lands on 1, 2 or 5 times a power of ten */
  function niceStep(span, want) {
    const raw = span / Math.max(1, want);
    const mag = Math.pow(10, Math.floor(Math.log10(Math.max(1e-12, raw))));
    const norm = raw / mag;
    return (norm <= 1 ? 1 : norm <= 2 ? 2 : norm <= 5 ? 5 : 10) * mag;
  }
  function fmtTick(v, step) {
    const dec = Math.max(0, Math.min(6, Math.ceil(-Math.log10(step)) + (step < 1 ? 0 : 0)));
    if (Math.abs(v) >= 10000) return (v / 1000).toFixed(0) + 'k';
    return (+v).toFixed(dec);
  }

  /* ================================================================
     traces: one line per run, over the generations
     ================================================================ */
  function lines(series, opts) {
    opts = opts || {};
    const W = opts.width || 760, H = opts.height || 260;
    const mx = 68, my = 42, top = 18, right = 16;
    let xMin = Infinity, xMax = -Infinity, yMin = Infinity, yMax = -Infinity;
    series.forEach(s => {
      s.x.forEach(v => { if (v < xMin) xMin = v; if (v > xMax) xMax = v; });
      s.y.forEach(v => { if (v < yMin) yMin = v; if (v > yMax) yMax = v; });
    });
    if (!isFinite(xMin)) return '';
    if (yMax - yMin < 1e-12) { yMax += 0.5; yMin -= 0.5; }
    const pad = (yMax - yMin) * 0.06;
    yMin -= pad; yMax += pad;
    const X = v => mx + (W - mx - right) * (v - xMin) / Math.max(1e-12, xMax - xMin);
    const Y = v => (H - my) - (H - my - top) * (v - yMin) / Math.max(1e-12, yMax - yMin);

    let s = '';
    s += line(mx, H - my, W - right, H - my, V('border-strong'), 1.4);
    s += line(mx, top, mx, H - my, V('border-strong'), 1.4);
    const sx = niceStep(xMax - xMin, 6), sy = niceStep(yMax - yMin, 5);
    for (let v = Math.ceil(xMin / sx) * sx; v <= xMax; v += sx) {
      s += line(X(v), H - my, X(v), H - my + 5, V('border-strong'), 1.1);
      s += txt(X(v), H - my + 17, fmtTick(v, sx), 'art-mut', 9, 'middle');
    }
    for (let v = Math.ceil(yMin / sy) * sy; v <= yMax; v += sy) {
      s += line(mx - 5, Y(v), W - right, Y(v), V('border'), 0.7, 'opacity="0.45"');
      s += txt(mx - 8, Y(v) + 3, fmtTick(v, sy), 'art-mut', 9, 'end');
    }
    series.forEach((ser, i) => {
      const col = V(COLOURS[i % COLOURS.length]);
      const pts = ser.x.map((v, k) => [X(v), Y(ser.y[k])]);
      s += pathEl(poly(pts), col, 1.1, 'opacity="0.85"');
    });
    if (opts.xLabel) s += txt(mx + (W - mx - right) / 2, H - 10, opts.xLabel, 'art-mut', 10, 'middle');
    if (opts.yLabel) s += txt(-(top + (H - my - top) / 2), 14, opts.yLabel, 'art-mut', 10, 'middle', 'transform="rotate(-90)"');
    series.forEach((ser, i) => {
      const x = W - right - 150 + (i % 2) * 74, y = top + 4 + Math.floor(i / 2) * 14;
      s += rect(x, y - 6, 9, 3, V(COLOURS[i % COLOURS.length]));
      s += txt(x + 13, y - 1, ser.name, 'art-mut', 9);
    });
    return svg(`0 0 ${W} ${H}`, s);
  }

  /* ================================================================
     the posterior of one parameter, with its highest-density interval
     ================================================================ */
  function histogram(values, opts) {
    opts = opts || {};
    const W = opts.width || 760, H = opts.height || 240;
    const mx = 62, my = 44, top = 16, right = 16;
    const n = values.length;
    if (!n) return '';
    let lo = Infinity, hi = -Infinity;
    for (let i = 0; i < n; i++) { const v = values[i]; if (v < lo) lo = v; if (v > hi) hi = v; }
    if (hi - lo < 1e-12) { hi = lo + 1e-6; }
    const bins = Math.max(8, Math.min(120, opts.bins || 40));
    const w = (hi - lo) / bins;
    const count = new Float64Array(bins);
    for (let i = 0; i < n; i++) {
      let b = Math.floor((values[i] - lo) / w);
      if (b >= bins) b = bins - 1;
      if (b < 0) b = 0;
      count[b]++;
    }
    let maxC = 0;
    for (let b = 0; b < bins; b++) if (count[b] > maxC) maxC = count[b];
    const X = v => mx + (W - mx - right) * (v - lo) / (hi - lo);
    const Y = c => (H - my) - (H - my - top) * c / Math.max(1, maxC);

    let s = '';
    s += line(mx, H - my, W - right, H - my, V('border-strong'), 1.4);
    s += line(mx, top, mx, H - my, V('border-strong'), 1.4);
    for (let b = 0; b < bins; b++) {
      const x0 = X(lo + b * w), x1 = X(lo + (b + 1) * w);
      s += rect(x0, Y(count[b]), x1 - x0 - 0.6, (H - my) - Y(count[b]), V('accent'), 'opacity="0.68"');
    }
    const sx = niceStep(hi - lo, 6);
    for (let v = Math.ceil(lo / sx) * sx; v <= hi; v += sx) {
      s += line(X(v), H - my, X(v), H - my + 5, V('border-strong'), 1.1);
      s += txt(X(v), H - my + 17, fmtTick(v, sx), 'art-mut', 9, 'middle');
    }
    const sy = niceStep(maxC, 4);
    for (let c = sy; c <= maxC; c += sy) {
      s += line(mx - 5, Y(c), mx, Y(c), V('border-strong'), 1.1);
      s += txt(mx - 8, Y(c) + 3, fmtTick(c, sy), 'art-mut', 9, 'end');
    }
    (opts.marks || []).forEach((m, i) => {
      if (m.x == null || m.x < lo || m.x > hi) return;
      s += line(X(m.x), top, X(m.x), H - my, V(i === 2 ? 'danger' : 'sky'), 1.6, 'stroke-dasharray="4 3"');
      if (m.label) s += txt(X(m.x) + 3, top + 10 + i * 11, m.label, 'art-mut', 9);
    });
    if (opts.xLabel) s += txt(mx + (W - mx - right) / 2, H - 10, opts.xLabel, 'art-mut', 10, 'middle');
    if (opts.yLabel) s += txt(-(top + (H - my - top) / 2), 13, opts.yLabel, 'art-mut', 10, 'middle', 'transform="rotate(-90)"');
    return svg(`0 0 ${W} ${H}`, s);
  }

  Object.assign(Plots6, { lines, histogram });
  window.Plots6 = Plots6;
})();
