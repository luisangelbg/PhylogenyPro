/* PhylogenyPro — molecular phylogenetics in the browser.
   Copyright (C) 2026  Luis Ángel Barrera-Guzmán

   This program is free software: you can redistribute it and/or modify it under
   the terms of the GNU General Public License as published by the Free Software
   Foundation, either version 3 of the License, or (at your option) any later
   version. It is distributed in the hope that it will be useful, but WITHOUT ANY
   WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
   A PARTICULAR PURPOSE. See the GNU General Public License, in the file LICENSE
   at the root of this program, or <https://www.gnu.org/licenses/>. */

/* PhylogenyPro — figures of Block 9: what a character did on the tree.

   Four pictures, each answering a question the numbers alone cannot:

   traitgram   Felsenstein's phenogram. Time runs across, the trait value runs
               up, and the reconstructed value at every node is placed where the
               model says. It is the one figure that shows a continuous trait and
               the tree in the same plane, and where a group of species has ended
               up far from its relatives it is visible at once.
   pies        the tree with a pie at every node. Its virtue is honesty: a node
               reconstructed at 55 % looks half-and-half, which a coloured branch
               never would.
   scatter     two traits against each other, with the contrasts and the
               generalised fit drawn over the raw points, so the difference
               between the two regressions can be seen rather than asserted.
   stateBar    how much of the tree each state occupied, from the stochastic map.

   Colours come from the theme tokens, so everything follows light and dark. */

const Plots9 = {};

(function () {

  const V = n => `var(--${n})`;
  const f1 = v => (+v).toFixed(1);
  const esc2 = s => String(s == null ? '' : s).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
  const svg = (vb, inner, plot) => `<svg viewBox="${vb}"${plot ? ` data-plot="${plot}"` : ''} xmlns="http://www.w3.org/2000/svg">${inner}</svg>`;
  /* marks for the figure studio: the area of the data (data-plot) and the
     legend, in one group (data-role="legend") with its entries numbered (data-li) */
  const plotArea = (x, y, w, h) => [x, y, w, h].map(v => +(+v).toFixed(2)).join(' ');
  const line = (x1, y1, x2, y2, st, w, ex) => `<line x1="${f1(x1)}" y1="${f1(y1)}" x2="${f1(x2)}" y2="${f1(y2)}" stroke="${st}" stroke-width="${w || 1}" ${ex || ''}/>`;
  const circ = (cx, cy, r, fill, ex) => `<circle cx="${f1(cx)}" cy="${f1(cy)}" r="${f1(r)}" fill="${fill}" ${ex || ''}/>`;
  const rect = (x, y, w, h, fill, ex) => `<rect x="${f1(x)}" y="${f1(y)}" width="${f1(Math.max(0, w))}" height="${f1(Math.max(0, h))}" fill="${fill}" ${ex || ''}/>`;
  function txt(x, y, s, cls, size, anchor, ex) {
    return `<text x="${f1(x)}" y="${f1(y)}" class="${cls || 'art-mut'}" font-size="${size || 9}" text-anchor="${anchor || 'start'}" ${ex || ''}>${esc2(s)}</text>`;
  }
  function niceStep(span, want) {
    const raw = span / Math.max(1, want);
    const mag = Math.pow(10, Math.floor(Math.log10(Math.max(1e-12, raw))));
    const n = raw / mag;
    return (n <= 1 ? 1 : n <= 2 ? 2 : n <= 5 ? 5 : 10) * mag;
  }
  const fmtTick = (v, step) => (+v).toFixed(Math.max(0, Math.min(6, Math.ceil(-Math.log10(step)))));

  /* the palette for discrete states: the app's own accents, in an order that
     stays distinguishable when printed in grey */
  const STATE_COLOURS = ['c2', 'c1', 'c3', 'c4', 'c5', 'c7', 'c8', 'c6', 'c9', 'c10'];
  const stateColour = i => V(STATE_COLOURS[i % STATE_COLOURS.length]);

  /* ================================================================
     the traitgram
     ================================================================ */
  /* nodes: [{node, value, depth, parent}] with the tips included */
  function traitgram(points, edges, opts) {
    opts = opts || {};
    const W = opts.width || 740, H = opts.height || 400;
    const mx = 62, my = 48, top = 18, right = opts.labels ? 120 : 20;
    const ds = points.map(p => p.depth), vs = points.map(p => p.value);
    const d0 = 0, d1 = Math.max.apply(null, ds);
    let v0 = Math.min.apply(null, vs), v1 = Math.max.apply(null, vs);
    const pad = (v1 - v0) * 0.08 || 0.5;
    v0 -= pad; v1 += pad;
    const X = d => mx + (W - mx - right) * (d - d0) / Math.max(1e-9, d1 - d0);
    const Y = v => (H - my) - (H - my - top) * (v - v0) / Math.max(1e-12, v1 - v0);

    let s = '';
    s += line(mx, H - my, W - right, H - my, V('border-strong'), 1.4);
    s += line(mx, top, mx, H - my, V('border-strong'), 1.4);
    const sx = niceStep(d1 - d0, 6), sy = niceStep(v1 - v0, 5);
    for (let d = 0; d <= d1 + 1e-9; d += sx) {
      s += line(X(d), H - my, X(d), H - my + 5, V('border-strong'), 1.1);
      s += txt(X(d), H - my + 17, fmtTick(opts.reverseTime ? d1 - d : d, sx), 'art-mut', 9, 'middle');
    }
    for (let v = Math.ceil(v0 / sy) * sy; v <= v1; v += sy) {
      s += line(mx - 5, Y(v), W - right, Y(v), V('border'), 0.7, 'opacity="0.3"');
      s += txt(mx - 8, Y(v) + 3, fmtTick(v, sy), 'art-mut', 9, 'end');
    }
    const byId = new Map(points.map(p => [p.node, p]));
    edges.forEach(e => {
      const a = byId.get(e.from), b = byId.get(e.to);
      if (!a || !b) return;
      s += line(X(a.depth), Y(a.value), X(b.depth), Y(b.value), V('ink'), 1.1, 'opacity="0.55"');
    });
    /* The tips all sit at the right edge, at whatever value they reached, so
       their names pile up wherever the trait is crowded. They are nudged apart
       vertically first — with a thin leader back to the point, so no name ends
       up claiming a value it does not have — and whatever still collides after
       that is dropped and counted, because two names printed on top of each
       other say less than one name and an honest note. */
    let hiddenLabels = 0;
    if (opts.labels && opts.showLabels !== false) {
      const tipsP = points.filter(p => p.tip);
      const fs = 8.5;
      const spread = window.TreeView
        ? TreeView.spreadLabels(tipsP.map(p => Y(p.value)), fs + 1.5, top + fs, H - my - 2)
        : { pos: tipsP.map(p => Y(p.value)) };
      const boxes = tipsP.map((p, i) => {
        const nm = opts.labels[p.tipRow] || '';
        return { x: X(p.depth) + 5, y: spread.pos[i] - fs, w: nm.length * fs * 0.52, h: fs + 1, rank: i };
      });
      const keep = window.TreeView ? TreeView.dropColliding(boxes) : boxes.map(() => true);
      tipsP.forEach((p, i) => {
        if (!keep[i]) { hiddenLabels++; return; }
        const y0 = Y(p.value), y1 = spread.pos[i];
        if (Math.abs(y1 - y0) > 1.5) {
          s += line(X(p.depth) + 2, y0, X(p.depth) + 4.5, y1, V('border'), 0.7, 'opacity="0.7"');
        }
        s += txt(X(p.depth) + 5, y1 + 3, opts.labels[p.tipRow] || '', 'art-txt', fs, 'start', 'font-style="italic"');
      });
    }
    Plots9.lastInfo = { hiddenLabels };

    points.forEach(p => {
      if (p.tip) {
        s += circ(X(p.depth), Y(p.value), 2.6, V('accent'), 'opacity="0.9"');
      } else {
        s += circ(X(p.depth), Y(p.value), 1.9, V('sky'), 'opacity="0.7"');
      }
    });
    s += txt(mx + (W - mx - right) / 2, H - 13, opts.xLabel || '', 'art-mut', 10, 'middle');
    s += txt(-(top + (H - my - top) / 2), 14, opts.yLabel || '', 'art-mut', 10, 'middle', 'transform="rotate(-90)"');
    return svg(`0 0 ${W} ${H}`, s);
  }

  /* ================================================================
     the tree with a pie at every node
     ================================================================ */
  /* tree: the usual node format; probs: Map from the flattened node index to an
     array of state probabilities; tipState: array of state indices per tip row */
  function pies(tree, probs, opts) {
    opts = opts || {};
    const labels = opts.labels;
    const tips = Tree.tips(tree);
    const n = tips.length;
    const rowH = opts.rowHeight || 16;
    const pad = 14, tipFont = opts.tipFont || 10.5;
    const longest = tips.reduce((a, b) => {
      const la = String((labels ? labels[a.tip] : a.label) || ''), lb = String((labels ? labels[b.tip] : b.label) || '');
      return la.length > lb.length ? a : b;
    }, tips[0]);
    const labW = Math.min(220, String((labels ? labels[longest.tip] : longest.label) || '').length * tipFont * 0.56 + 26);
    const W = opts.width || 760;
    const legendH = opts.levels ? 22 : 0;
    const H = pad * 2 + n * rowH + 30 + legendH;
    const x0 = pad + 4, x1 = W - pad - labW - 8;
    const y0 = pad + legendH, y1 = H - pad - 30;

    const F = Tree.flatten(tree);
    const depth = new Float64Array(F.n);
    for (let i = F.post.length - 1; i >= 0; i--) {
      const k = F.post[i], p = F.parent[k];
      if (p >= 0) depth[k] = depth[p] + (opts.cladogram ? 1 : (F.len[k] || 0));
    }
    let maxD = 0;
    for (let k = 0; k < F.n; k++) if (F.isTip[k]) maxD = Math.max(maxD, depth[k]);
    const X = d => x0 + (x1 - x0) * (maxD > 0 ? d / maxD : 0);
    /* every tip gets a row, every internal node the mean of its daughters */
    const yOf = new Float64Array(F.n);
    let row = 0;
    for (let i = 0; i < F.post.length; i++) {
      const k = F.post[i];
      if (F.isTip[k]) { yOf[k] = y0 + rowH / 2 + (row++) * rowH; continue; }
      const ys = F.kids[k].map(c => yOf[c]);
      yOf[k] = (Math.min.apply(null, ys) + Math.max.apply(null, ys)) / 2;
    }

    let s = '';
    for (let k = 0; k < F.n; k++) {
      const p = F.parent[k];
      if (p < 0) continue;
      s += line(X(depth[p]), yOf[k], X(depth[k]), yOf[k], V('border-strong'), 1.5);
    }
    for (let k = 0; k < F.n; k++) {
      if (F.isTip[k]) continue;
      const ys = F.kids[k].map(c => yOf[c]);
      s += line(X(depth[k]), Math.min.apply(null, ys), X(depth[k]), Math.max.apply(null, ys), V('border-strong'), 1.5);
    }
    /* the tips, as a filled square in their own state */
    for (let k = 0; k < F.n; k++) {
      if (!F.isTip[k]) continue;
      const st = opts.tipState ? opts.tipState[F.tipRow[k]] : -1;
      if (st >= 0) s += rect(X(depth[k]) + 2, yOf[k] - 3.4, 6.8, 6.8, stateColour(st), 'rx="1.4"');
      if (labels !== false) {
        const lab = labels ? labels[F.tipRow[k]] : (F.nodes[k].label || '');
        s += txt(X(depth[k]) + 12, yOf[k] + 3.4, lab, 'art-txt', tipFont, 'start', 'font-style="italic"');
      }
    }
    /* the pies */
    const r = opts.pieRadius || 5.2;
    for (let k = 0; k < F.n; k++) {
      if (F.isTip[k]) continue;
      const pr = probs.get ? probs.get(k) : probs[k];
      if (!pr) continue;
      const cx = X(depth[k]), cy = yOf[k];
      let a0 = -Math.PI / 2;
      pr.forEach((p, i) => {
        if (p <= 1e-6) return;
        const a1 = a0 + 2 * Math.PI * p;
        if (p > 0.9999) { s += circ(cx, cy, r, stateColour(i)); a0 = a1; return; }
        const large = (a1 - a0) > Math.PI ? 1 : 0;
        const px0 = cx + r * Math.cos(a0), py0 = cy + r * Math.sin(a0);
        const px1 = cx + r * Math.cos(a1), py1 = cy + r * Math.sin(a1);
        s += `<path d="M${f1(cx)} ${f1(cy)} L${f1(px0)} ${f1(py0)} A${f1(r)} ${f1(r)} 0 ${large} 1 ${f1(px1)} ${f1(py1)} Z" fill="${stateColour(i)}"/>`;
        a0 = a1;
      });
      s += circ(cx, cy, r, 'none', `stroke="${V('bg')}" stroke-width="0.8"`);
    }
    /* the scale bar */
    if (!opts.cladogram && maxD > 0) {
      const step = niceStep(maxD, 5);
      const bx = x0, by = H - pad - 12;
      s += line(bx, by, X(step), by, V('border-strong'), 1.4);
      s += line(bx, by - 4, bx, by + 4, V('border-strong'), 1.2);
      s += line(X(step), by - 4, X(step), by + 4, V('border-strong'), 1.2);
      s += txt((bx + X(step)) / 2, by + 15, fmtTick(step, step) + (opts.unit ? ' ' + opts.unit : ''), 'art-mut', 9, 'middle');
    }
    /* the legend */
    if (opts.levels) {
      let lx = pad + 4;
      s += '<g data-role="legend">';
      opts.levels.forEach((lv, i) => {
        s += rect(lx, pad - 2, 9, 9, stateColour(i), `rx="1.6" data-li="${i}"`);
        s += txt(lx + 13, pad + 6, String(lv), 'art-txt', 10, 'start', `data-li="${i}"`);
        lx += 22 + String(lv).length * 6;
      });
      s += '</g>';
    }
    return svg(`0 0 ${W} ${H}`, s, plotArea(x0, y0, x1 - x0, y1 - y0));
  }

  /* ================================================================
     two traits, three regressions
     ================================================================ */
  /* raw: [{x,y,label}]; fits: [{slope, intercept, label, colour, dashed}] */
  function scatter(raw, fits, opts) {
    opts = opts || {};
    const W = opts.width || 620, H = opts.height || 330;
    const mx = 62, my = 50, top = 18, right = 18;
    const xs = raw.map(p => p.x), ys = raw.map(p => p.y);
    let x0 = Math.min.apply(null, xs), x1 = Math.max.apply(null, xs);
    let y0 = Math.min.apply(null, ys), y1 = Math.max.apply(null, ys);
    const px = (x1 - x0) * 0.06 || 0.5, py = (y1 - y0) * 0.08 || 0.5;
    x0 -= px; x1 += px; y0 -= py; y1 += py;
    const X = v => mx + (W - mx - right) * (v - x0) / Math.max(1e-12, x1 - x0);
    const Y = v => (H - my) - (H - my - top) * (v - y0) / Math.max(1e-12, y1 - y0);
    let s = '';
    s += line(mx, H - my, W - right, H - my, V('border-strong'), 1.3);
    s += line(mx, top, mx, H - my, V('border-strong'), 1.3);
    const sx = niceStep(x1 - x0, 6), sy = niceStep(y1 - y0, 5);
    for (let v = Math.ceil(x0 / sx) * sx; v <= x1; v += sx) {
      s += line(X(v), H - my, X(v), H - my + 5, V('border-strong'), 1.1);
      s += txt(X(v), H - my + 16, fmtTick(v, sx), 'art-mut', 9, 'middle');
    }
    for (let v = Math.ceil(y0 / sy) * sy; v <= y1; v += sy) {
      s += line(mx - 5, Y(v), W - right, Y(v), V('border'), 0.7, 'opacity="0.3"');
      s += txt(mx - 8, Y(v) + 3, fmtTick(v, sy), 'art-mut', 9, 'end');
    }
    if (opts.zeroLines) {
      if (x0 < 0 && x1 > 0) s += line(X(0), top, X(0), H - my, V('border-strong'), 0.9, 'stroke-dasharray="3 3"');
      if (y0 < 0 && y1 > 0) s += line(mx, Y(0), W - right, Y(0), V('border-strong'), 0.9, 'stroke-dasharray="3 3"');
    }
    raw.forEach(p => { s += circ(X(p.x), Y(p.y), 3, V('accent'), 'opacity="0.7"'); });
    /* the legend entries are gathered with their fits and drawn after them, in one group */
    let lg = '';
    (fits || []).forEach((f, i) => {
      const col = f.colour ? V(f.colour) : V('sky');
      s += line(X(x0), Y(f.intercept + f.slope * x0), X(x1), Y(f.intercept + f.slope * x1), col, 1.8,
        f.dashed ? 'stroke-dasharray="5 4"' : '');
      lg += line(W - right - 152, top + 6 + i * 13, W - right - 134, top + 6 + i * 13, col, 1.8, (f.dashed ? 'stroke-dasharray="5 4" ' : '') + `data-li="${i}"`);
      lg += txt(W - right - 130, top + 9 + i * 13, f.label || '', 'art-mut', 9, 'start', `data-li="${i}"`);
    });
    s += `<g data-role="legend">${lg}</g>`;
    s += txt(mx + (W - mx - right) / 2, H - 12, opts.xLabel || '', 'art-mut', 10, 'middle');
    s += txt(-(top + (H - my - top) / 2), 14, opts.yLabel || '', 'art-mut', 10, 'middle', 'transform="rotate(-90)"');
    return svg(`0 0 ${W} ${H}`, s, plotArea(mx, top, W - mx - right, H - my - top));
  }

  /* ================================================================
     how much of the tree each state held
     ================================================================ */
  function stateBar(proportions, levels, opts) {
    opts = opts || {};
    const W = opts.width || 620, H = opts.height || 96;
    const mx = 16, top = 26;
    const barH = 30;
    let s = '';
    let x = mx;
    const inner = W - mx * 2;
    proportions.forEach((p, i) => {
      const w = inner * p;
      s += rect(x, top, w, barH, stateColour(i), 'rx="2"');
      if (w > 34) s += txt(x + w / 2, top + barH / 2 + 4, (p * 100).toFixed(1) + ' %', 'art-txt', 10.5, 'middle', `fill="${V('bg')}"`);
      x += w;
    });
    let lx = mx;
    s += '<g data-role="legend">';
    levels.forEach((lv, i) => {
      s += rect(lx, 6, 9, 9, stateColour(i), `rx="1.6" data-li="${i}"`);
      s += txt(lx + 13, 14, String(lv) + (opts.times ? `  ${opts.times[i].toFixed(2)}` : ''), 'art-txt', 10, 'start', `data-li="${i}"`);
      lx += 30 + (String(lv).length + (opts.times ? 6 : 0)) * 6;
    });
    s += '</g>';
    if (opts.caption) s += txt(mx, top + barH + 18, opts.caption, 'art-mut', 9.5);
    return svg(`0 0 ${W} ${H}`, s, plotArea(mx, top, inner, barH));
  }

  Object.assign(Plots9, { traitgram, pies, scatter, stateBar, stateColour, STATE_COLOURS });
  window.Plots9 = Plots9;
})();
