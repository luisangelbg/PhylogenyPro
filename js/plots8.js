/* PhylogenyPro — molecular phylogenetics in the browser.
   Copyright (C) 2026  Luis Ángel Barrera-Guzmán

   This program is free software: you can redistribute it and/or modify it under
   the terms of the GNU General Public License as published by the Free Software
   Foundation, either version 3 of the License, or (at your option) any later
   version. It is distributed in the hope that it will be useful, but WITHOUT ANY
   WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
   A PARTICULAR PURPOSE. See the GNU General Public License, in the file LICENSE
   at the root of this program, or <https://www.gnu.org/licenses/>. */

/* PhylogenyPro — figures of Block 8: the shape of a radiation.

   The lineages-through-time plot is the one figure of this block that everybody
   reads, and it is read on a logarithmic axis: under a constant rate the curve
   is a straight line there, so a bend is the whole message. The straight line a
   constant rate would give is drawn behind it, because a curve without that
   reference says nothing. */

const Plots8 = {};

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
     lineages through time
     ================================================================ */
  /* curve: {times, lineages}; band (optional): {times, median, lower, upper} */
  function ltt(curve, opts) {
    opts = opts || {};
    const W = opts.width || 720, H = opts.height || 360;
    const mx = 70, my = 52, top = 18, right = 20;
    const maxT = curve.times[0];
    const maxN = curve.lineages[curve.lineages.length - 1];
    const X = t => mx + (W - mx - right) * (1 - t / maxT);
    const logN = n => Math.log(Math.max(1, n));
    const Y = n => (H - my) - (H - my - top) * logN(n) / logN(maxN);

    let s = '';
    s += line(mx, H - my, W - right, H - my, V('border-strong'), 1.4);
    s += line(mx, top, mx, H - my, V('border-strong'), 1.4);
    const st = niceStep(maxT, 6);
    for (let t = 0; t <= maxT + 1e-9; t += st) {
      s += line(X(t), H - my, X(t), H - my + 5, V('border-strong'), 1.1);
      s += txt(X(t), H - my + 17, fmtTick(t, st), 'art-mut', 9, 'middle');
    }
    /* a logarithmic vertical axis, labelled at the round numbers it contains */
    const ticks = [];
    for (let e = 0; Math.pow(10, e) <= maxN * 1.5; e++) {
      [1, 2, 5].forEach(m => { const v = m * Math.pow(10, e); if (v >= 1 && v <= maxN) ticks.push(v); });
    }
    if (ticks.indexOf(maxN) < 0) ticks.push(maxN);
    ticks.forEach(v => {
      s += line(mx - 5, Y(v), W - right, Y(v), V('border'), 0.7, 'opacity="0.35"');
      s += txt(mx - 8, Y(v) + 3, String(v), 'art-mut', 9, 'end');
    });

    /* the band a sample of trees gives */
    if (opts.band) {
      const b = opts.band;
      let d = '';
      b.times.forEach((t, i) => { d += (i ? 'L' : 'M') + f1(X(t)) + ' ' + f1(Y(b.upper[i])); });
      for (let i = b.times.length - 1; i >= 0; i--) d += 'L' + f1(X(b.times[i])) + ' ' + f1(Y(b.lower[i]));
      d += 'Z';
      s += `<path d="${d}" fill="${V('accent')}" opacity="0.16"/>`;
    }
    /* the straight line a constant rate would give */
    if (opts.expected !== false) {
      s += line(X(maxT), Y(2), X(0), Y(maxN), V('sky'), 1.6, 'stroke-dasharray="5 4" opacity="0.8"');
    }
    /* the curve itself, as the steps it really is */
    let d = 'M' + f1(X(curve.times[0])) + ' ' + f1(Y(curve.lineages[0]));
    for (let i = 1; i < curve.times.length; i++) {
      d += 'L' + f1(X(curve.times[i])) + ' ' + f1(Y(curve.lineages[i - 1]));
      d += 'L' + f1(X(curve.times[i])) + ' ' + f1(Y(curve.lineages[i]));
    }
    s += `<path d="${d}" stroke="${V('ink')}" stroke-width="1.8" fill="none"/>`;

    s += txt(mx + (W - mx - right) / 2, H - 14, opts.xLabel || T('tiempo antes del presente', 'time before the present'), 'art-mut', 10, 'middle');
    s += txt(-(top + (H - my - top) / 2), 15, T('linajes (escala logarítmica)', 'lineages (log scale)'), 'art-mut', 10, 'middle', 'transform="rotate(-90)"');
    s += line(W - right - 150, top + 6, W - right - 132, top + 6, V('ink'), 1.8);
    s += txt(W - right - 128, top + 9, T('observado', 'observed'), 'art-mut', 9);
    s += line(W - right - 150, top + 20, W - right - 132, top + 20, V('sky'), 1.6, 'stroke-dasharray="5 4"');
    s += txt(W - right - 128, top + 23, T('tasa constante', 'constant rate'), 'art-mut', 9);
    return svg(`0 0 ${W} ${H}`, s);
  }

  /* ================================================================
     the speciation rate a model implies
     ================================================================ */
  /* For the models whose rate depends on how many lineages there are, the
     horizontal axis is the number of lineages, not time: that is what the model
     actually says. */
  function rateCurve(points, opts) {
    opts = opts || {};
    const W = opts.width || 560, H = opts.height || 240;
    const mx = 68, my = 46, top = 16, right = 16;
    if (!points.length) return '';
    const xs = points.map(p => p.x), ys = points.map(p => p.y);
    const x0 = Math.min.apply(null, xs), x1 = Math.max.apply(null, xs);
    let y0 = Math.min.apply(null, ys), y1 = Math.max.apply(null, ys);
    if (y1 - y0 < 1e-12) { y0 = Math.max(0, y0 * 0.5); y1 = y1 * 1.5 + 1e-6; }
    const X = v => mx + (W - mx - right) * (v - x0) / Math.max(1e-9, x1 - x0);
    const Y = v => (H - my) - (H - my - top) * (v - y0) / Math.max(1e-12, y1 - y0);
    let s = '';
    s += line(mx, H - my, W - right, H - my, V('border-strong'), 1.3);
    s += line(mx, top, mx, H - my, V('border-strong'), 1.3);
    const sx = niceStep(x1 - x0, 6), sy = niceStep(y1 - y0, 4);
    for (let v = Math.ceil(x0 / sx) * sx; v <= x1; v += sx) {
      s += line(X(v), H - my, X(v), H - my + 5, V('border-strong'), 1.1);
      s += txt(X(v), H - my + 16, fmtTick(v, sx), 'art-mut', 9, 'middle');
    }
    for (let v = Math.ceil(y0 / sy) * sy; v <= y1; v += sy) {
      s += line(mx - 5, Y(v), W - right, Y(v), V('border'), 0.7, 'opacity="0.35"');
      s += txt(mx - 8, Y(v) + 3, fmtTick(v, sy), 'art-mut', 9, 'end');
    }
    let d = '';
    points.forEach((p, i) => { d += (i ? 'L' : 'M') + f1(X(p.x)) + ' ' + f1(Y(p.y)); });
    s += `<path d="${d}" stroke="${V('accent')}" stroke-width="2" fill="none"/>`;
    s += txt(mx + (W - mx - right) / 2, H - 10, opts.xLabel || '', 'art-mut', 10, 'middle');
    s += txt(-(top + (H - my - top) / 2), 14, opts.yLabel || T('tasa de especiación', 'speciation rate'), 'art-mut', 10, 'middle', 'transform="rotate(-90)"');
    return svg(`0 0 ${W} ${H}`, s);
  }

  /* ================================================================
     a rate for every species, ranked
     ================================================================ */
  function ranked(values, opts) {
    opts = opts || {};
    const W = opts.width || 620, H = opts.height || 240;
    const mx = 66, my = 46, top = 16, right = 16;
    const v = values.slice().sort((a, b) => b.value - a.value);
    if (!v.length) return '';
    const maxV = Math.max.apply(null, v.map(x => x.value));
    const bw = (W - mx - right) / v.length;
    let s = '';
    s += line(mx, H - my, W - right, H - my, V('border-strong'), 1.3);
    s += line(mx, top, mx, H - my, V('border-strong'), 1.3);
    const sy = niceStep(maxV, 4);
    for (let y = 0; y <= maxV; y += sy) {
      const yy = (H - my) - (H - my - top) * y / maxV;
      s += line(mx - 5, yy, W - right, yy, V('border'), 0.7, 'opacity="0.35"');
      s += txt(mx - 8, yy + 3, fmtTick(y, sy), 'art-mut', 9, 'end');
    }
    v.forEach((x, i) => {
      const h = (H - my - top) * x.value / maxV;
      s += rect(mx + i * bw, (H - my) - h, Math.max(1, bw - 0.7), h, V('accent'), 'opacity="0.75"');
    });
    s += txt(mx + (W - mx - right) / 2, H - 10, opts.xLabel || T('especies, de mayor a menor', 'species, largest first'), 'art-mut', 10, 'middle');
    s += txt(-(top + (H - my - top) / 2), 14, opts.yLabel || 'DR', 'art-mut', 10, 'middle', 'transform="rotate(-90)"');
    /* the three fastest, named */
    v.slice(0, 3).forEach((x, i) => {
      s += txt(W - right - 6, top + 12 + i * 12, `${x.name} ${x.value.toFixed(3)}`, 'art-mut', 9, 'end');
    });
    return svg(`0 0 ${W} ${H}`, s);
  }

  Object.assign(Plots8, { ltt, rateCurve, ranked });
  window.Plots8 = Plots8;
})();
