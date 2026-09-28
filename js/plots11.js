/* PhylogenyPro — molecular phylogenetics in the browser.
   Copyright (C) 2026  Luis Ángel Barrera-Guzmán

   This program is free software: you can redistribute it and/or modify it under
   the terms of the GNU General Public License as published by the Free Software
   Foundation, either version 3 of the License, or (at your option) any later
   version. It is distributed in the hope that it will be useful, but WITHOUT ANY
   WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
   A PARTICULAR PURPOSE. See the GNU General Public License, in the file LICENSE
   at the root of this program, or <https://www.gnu.org/licenses/>. */

/* PhylogenyPro — figures of Block 11: the disagreement, drawn.

   tanglegram   two trees facing each other with a line per tip. The lines cross
                exactly where the trees disagree, and the crossings are the
                figure: a tanglegram with none of them says the two trees are
                the same, which no table says as fast.
   concordance  every branch of a reference tree as a stacked bar — how much of
                the evidence is with it, how much is with each of the two
                alternatives, and how much cannot say. A branch whose bar is
                mostly not its own colour is a branch to be careful about, even
                at 100 % bootstrap.
   network      the outer cycle of a split network. Where the splits agree the
                cycle is thin and tree-like; where they conflict it opens into
                boxes, and the width of a box is how much support the conflict
                has.
   distanceHeat the pairwise distances between a set of trees, so that "these
                two disagree" becomes "these two, and not those".
   dBars        Patterson's D for a set of trios, with the interval the block
                jackknife gives and the line at zero that the test is about.

   Colours come from the theme tokens, so everything follows light and dark. */

const Plots11 = {};

(function () {

  const V = n => `var(--${n})`;
  const f1 = v => (+v).toFixed(1);
  const esc2 = s => String(s == null ? '' : s).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
  const svg = (vb, inner) => `<svg viewBox="${vb}" xmlns="http://www.w3.org/2000/svg">${inner}</svg>`;
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

  /* the layout every tree figure here needs: a y for every node and an x that
     is either the accumulated branch length or the depth in edges */
  function layout(tree, opts) {
    const F = Tree.flatten(tree);
    const depth = new Float64Array(F.n);
    for (let i = F.post.length - 1; i >= 0; i--) {
      const k = F.post[i], p = F.parent[k];
      if (p >= 0) depth[k] = depth[p] + (opts.cladogram ? 1 : (F.len[k] || 0));
    }
    let maxD = 0;
    for (let k = 0; k < F.n; k++) if (F.isTip[k]) maxD = Math.max(maxD, depth[k]);
    const y = new Float64Array(F.n);
    let row = 0;
    for (let i = 0; i < F.post.length; i++) {
      const k = F.post[i];
      if (F.isTip[k]) { y[k] = row++; continue; }
      const ys = F.kids[k].map(c => y[c]);
      y[k] = (Math.min.apply(null, ys) + Math.max.apply(null, ys)) / 2;
    }
    return { F, depth, y, maxD, nTips: row };
  }

  /* ================================================================
     tanglegram
     ================================================================ */
  function tanglegram(left, right, opts) {
    opts = opts || {};
    const labels = opts.labels || [];
    const L = layout(left, opts), R = layout(right, opts);
    const n = L.nTips;
    const rowH = opts.rowHeight || (n > 40 ? 12 : 17);
    const pad = 14;
    const gap = opts.gap || 190;                      // room for the labels and links
    const treeW = opts.treeWidth || 200;
    const W = opts.width || (pad * 2 + treeW * 2 + gap);
    const H = pad * 2 + n * rowH + 26;
    const y0 = pad + rowH / 2;
    const Y = r => y0 + r * rowH;
    const XL = d => pad + treeW * (L.maxD > 0 ? d / L.maxD : 0);
    const XR = d => W - pad - treeW * (R.maxD > 0 ? d / R.maxD : 0);

    let s = '';
    const drawTree = (lay, X, flip) => {
      let t = '';
      for (let k = 0; k < lay.F.n; k++) {
        const p = lay.F.parent[k];
        if (p < 0) continue;
        t += line(X(lay.depth[p]), Y(lay.y[k]), X(lay.depth[k]), Y(lay.y[k]), V('border-strong'), 1.4);
      }
      for (let k = 0; k < lay.F.n; k++) {
        if (lay.F.isTip[k]) continue;
        const ys = lay.F.kids[k].map(c => lay.y[c]);
        t += line(X(lay.depth[k]), Y(Math.min.apply(null, ys)), X(lay.depth[k]), Y(Math.max.apply(null, ys)), V('border-strong'), 1.4);
      }
      void flip;
      return t;
    };
    s += drawTree(L, XL, false);
    s += drawTree(R, XR, true);

    /* the links, and their crossings, which are what the figure is for */
    const rowOfLeft = new Map(), rowOfRight = new Map();
    for (let k = 0; k < L.F.n; k++) if (L.F.isTip[k]) rowOfLeft.set(L.F.tipRow[k], L.y[k]);
    for (let k = 0; k < R.F.n; k++) if (R.F.isTip[k]) rowOfRight.set(R.F.tipRow[k], R.y[k]);
    const xa = pad + treeW + 6, xb = W - pad - treeW - 6;
    const highlight = new Set(opts.highlight || []);
    rowOfLeft.forEach((ry, tip) => {
      if (!rowOfRight.has(tip)) return;
      const y1 = Y(ry), y2 = Y(rowOfRight.get(tip));
      const moved = Math.abs(ry - rowOfRight.get(tip)) > 0.5;
      const col = highlight.has(tip) ? V('c2') : moved ? V('c4') : V('border');
      s += `<path d="M${f1(xa)} ${f1(y1)} C${f1(xa + 40)} ${f1(y1)}, ${f1(xb - 40)} ${f1(y2)}, ${f1(xb)} ${f1(y2)}" fill="none" stroke="${col}" stroke-width="${moved ? 1.6 : 1}" opacity="${moved ? 0.9 : 0.5}"/>`;
      const lab = labels[tip] != null ? labels[tip] : String(tip);
      s += txt(xa + 4, y1 - 2.5, lab, moved ? 'art-txt' : 'art-mut', opts.tipFont || 9.5, 'start', 'font-style="italic"');
    });
    s += txt(pad, H - 8, opts.leftLabel || '', 'art-mut', 10);
    s += txt(W - pad, H - 8, opts.rightLabel || '', 'art-mut', 10, 'end');
    if (opts.crossings != null) {
      s += txt(W / 2, pad - 2, T(`${opts.crossings} cruces`, `${opts.crossings} crossings`), 'art-mut', 10, 'middle');
    }
    return svg(`0 0 ${W} ${H}`, s);
  }

  /* ================================================================
     concordance factors, branch by branch
     ================================================================ */
  /* rows: [{label, concordant, alt1, alt2, other}] as percentages */
  function concordance(rows, opts) {
    opts = opts || {};
    const W = opts.width || 720;
    const rowH = opts.rowHeight || 20;
    const mx = opts.labelWidth || 150, right = 70, top = 34;
    const H = top + rows.length * rowH + 26;
    const bw = W - mx - right;
    let s = '';
    const cols = [V('c3'), V('c2'), V('c6'), V('border-strong')];
    const names = opts.legend || ['concordante', 'alternativa 1', 'alternativa 2', 'no decide'];
    let lx = mx;
    names.forEach((nm, i) => {
      s += rect(lx, 8, 9, 9, cols[i], 'rx="1.6"');
      s += txt(lx + 13, 16, nm, 'art-txt', 9.5);
      lx += 26 + String(nm).length * 5.6;
    });
    rows.forEach((r, i) => {
      const y = top + i * rowH;
      let x = mx;
      [r.concordant, r.alt1, r.alt2, r.other].forEach((v, k) => {
        const w = bw * Math.max(0, v || 0) / 100;
        if (w > 0.4) s += rect(x, y, w, rowH - 6, cols[k], `rx="2" opacity="${k === 3 ? 0.4 : 0.85}"`);
        x += w;
      });
      s += txt(mx - 8, y + rowH - 9, r.label, 'art-txt', 10, 'end');
      s += txt(mx + bw + 6, y + rowH - 9, `${(r.concordant || 0).toFixed(1)} %`, 'art-mut', 9.5);
    });
    for (let p = 0; p <= 100; p += 25) {
      const x = mx + bw * p / 100;
      s += line(x, top - 4, x, top + rows.length * rowH - 4, V('border'), 0.7, 'opacity="0.35"');
      s += txt(x, H - 8, p + ' %', 'art-mut', 9, 'middle');
    }
    /* a third of the evidence is what pure chance gives between three
       resolutions, and the line says so */
    const x33 = mx + bw / 3;
    s += line(x33, top - 4, x33, top + rows.length * rowH - 4, V('c4'), 1.2, 'stroke-dasharray="4 3"');
    s += txt(x33, top - 8, opts.chanceLabel || '⅓', 'art-mut', 9, 'middle');
    return svg(`0 0 ${W} ${H}`, s);
  }

  /* gCF against sCF, which is the plot that shows they are different questions */
  function concordanceScatter(rows, opts) {
    opts = opts || {};
    const W = opts.width || 420, H = opts.height || 360;
    const mx = 54, my = 46, top = 16, right = 16;
    const X = v => mx + (W - mx - right) * v / 100;
    const Y = v => (H - my) - (H - my - top) * v / 100;
    let s = '';
    s += line(mx, H - my, W - right, H - my, V('border-strong'), 1.3);
    s += line(mx, top, mx, H - my, V('border-strong'), 1.3);
    for (let p = 0; p <= 100; p += 25) {
      s += line(X(p), H - my, X(p), H - my + 5, V('border-strong'), 1.1);
      s += txt(X(p), H - my + 16, String(p), 'art-mut', 9, 'middle');
      s += line(mx - 5, Y(p), W - right, Y(p), V('border'), 0.7, 'opacity="0.3"');
      s += txt(mx - 8, Y(p) + 3, String(p), 'art-mut', 9, 'end');
    }
    s += line(X(0), Y(0), X(100), Y(100), V('border-strong'), 1, 'stroke-dasharray="4 3"');
    s += line(X(33.33), top, X(33.33), H - my, V('c4'), 1, 'stroke-dasharray="3 3"');
    s += line(mx, Y(33.33), W - right, Y(33.33), V('c4'), 1, 'stroke-dasharray="3 3"');
    rows.forEach(r => { s += circ(X(r.x), Y(r.y), 3.4, V('c1'), 'opacity="0.78"'); });
    /* Branch labels on a scatter land wherever the data land, and concordance
       factors cluster: with fifty branches the names become one black stain.
       The ones furthest from the diagonal are kept first — a branch whose gene
       and site concordance disagree is the one worth naming — and whatever
       still collides is dropped and counted. */
    let hiddenLabels = 0;
    if (opts.labelAll) {
      const fs = 8;
      const marked = rows.filter(r => r.label);
      const boxes = marked.map(r => ({
        x: X(r.x) + 5, y: Y(r.y) + 3 - fs, w: String(r.label).length * fs * 0.54, h: fs + 1,
        rank: -Math.abs((+r.x || 0) - (+r.y || 0)),
      }));
      const keep = (window.TreeView && opts.labelDeclutter !== false)
        ? TreeView.dropColliding(boxes) : boxes.map(() => true);
      marked.forEach((r, i) => {
        if (!keep[i]) { hiddenLabels++; return; }
        s += txt(X(r.x) + 5, Y(r.y) + 3, r.label, 'art-mut', fs);
      });
    }
    Plots11.lastInfo = { hiddenLabels };
    s += txt(mx + (W - mx - right) / 2, H - 10, opts.xLabel || 'gCF (%)', 'art-mut', 10, 'middle');
    s += txt(-(top + (H - my - top) / 2), 14, opts.yLabel || 'sCF (%)', 'art-mut', 10, 'middle', 'transform="rotate(-90)"');
    return svg(`0 0 ${W} ${H}`, s);
  }

  /* ================================================================
     the split network
     ================================================================ */
  function network(cycle, opts) {
    opts = opts || {};
    const labels = opts.labels || [];
    const pts = cycle.points;
    let x0 = Infinity, x1 = -Infinity, y0 = Infinity, y1 = -Infinity;
    pts.forEach(p => {
      if (p.x < x0) x0 = p.x; if (p.x > x1) x1 = p.x;
      if (p.y < y0) y0 = p.y; if (p.y > y1) y1 = p.y;
    });
    const pad = 70;
    const W = opts.width || 680, H = opts.height || 520;
    const sx = (W - 2 * pad) / Math.max(1e-9, x1 - x0);
    const sy = (H - 2 * pad) / Math.max(1e-9, y1 - y0);
    const k = Math.min(sx, sy);
    const X = v => pad + (v - x0) * k + (W - 2 * pad - (x1 - x0) * k) / 2;
    const Y = v => pad + (v - y0) * k + (H - 2 * pad - (y1 - y0) * k) / 2;

    let s = '';
    /* the cycle itself, each segment coloured by whether its split is one that
       some other split contradicts */
    const conflict = new Set(opts.conflicting || []);
    for (let i = 0; i < pts.length; i++) {
      const a = pts[i], b = pts[(i + 1) % pts.length];
      const sp = b.split != null ? b.split : a.split;
      const isCon = sp != null && conflict.has(sp);
      s += line(X(a.x), Y(a.y), X(b.x), Y(b.y), isCon ? V('c2') : V('border-strong'), isCon ? 2.2 : 1.6);
    }
    /* the tips */
    pts.forEach(p => {
      if (p.tip == null) return;
      const cx = X(p.x), cy = Y(p.y);
      s += circ(cx, cy, 3, V('c1'));
      const lab = labels[p.tip] != null ? labels[p.tip] : String(p.tip);
      /* push the label outwards from the centre of the figure */
      const dx = cx - W / 2, dy = cy - H / 2;
      const m = Math.hypot(dx, dy) || 1;
      s += txt(cx + 7 * dx / m, cy + 7 * dy / m + 3, lab, 'art-txt', opts.tipFont || 9.5,
        dx > 0 ? 'start' : 'end', 'font-style="italic"');
    });
    if (opts.caption) s += txt(12, H - 10, opts.caption, 'art-mut', 9.5);
    return svg(`0 0 ${W} ${H}`, s);
  }

  /* ================================================================
     the distances between a set of trees
     ================================================================ */
  function distanceHeat(M, names, opts) {
    opts = opts || {};
    const n = names.length;
    const cell = opts.cell || Math.max(22, Math.min(40, 420 / n));
    const left = opts.labelWidth || 92, top = opts.labelHeight || 76;
    const W = left + n * cell + 20, H = top + n * cell + 30;
    let hi = 0;
    for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) if (M[i][j] > hi) hi = M[i][j];
    if (hi <= 0) hi = 1;
    let s = '';
    for (let i = 0; i < n; i++) {
      s += txt(left - 8, top + i * cell + cell / 2 + 4, names[i], 'art-txt', 10, 'end');
      s += txt(left + i * cell + cell / 2, top - 8, names[i], 'art-txt', 10, 'start',
        `transform="rotate(-55 ${f1(left + i * cell + cell / 2)} ${f1(top - 8)})"`);
      for (let j = 0; j < n; j++) {
        const v = M[i][j], t = v / hi;
        s += rect(left + j * cell, top + i * cell, cell - 1.5, cell - 1.5,
          i === j ? V('bg-soft') : V('c1'), `rx="3" opacity="${i === j ? 0.45 : (0.1 + 0.85 * t).toFixed(3)}"`);
        if (i !== j && cell > 26) {
          s += txt(left + j * cell + cell / 2 - 0.75, top + i * cell + cell / 2 + 4,
            opts.integer ? String(Math.round(v)) : (+v).toFixed(2),
            t > 0.55 ? 'art-txt' : 'art-mut', 9, 'middle', t > 0.55 ? `fill="${V('bg')}"` : '');
        }
      }
    }
    if (opts.caption) s += txt(left, H - 8, opts.caption, 'art-mut', 9.5);
    return svg(`0 0 ${W} ${H}`, s);
  }

  /* ================================================================
     Patterson's D, with its interval
     ================================================================ */
  function dBars(rows, opts) {
    opts = opts || {};
    const W = opts.width || 680;
    const rowH = 22, top = 20, mx = opts.labelWidth || 190, right = 80;
    const H = top + rows.length * rowH + 34;
    let hi = 0.05;
    rows.forEach(r => { hi = Math.max(hi, Math.abs(r.D) + 1.96 * (r.se || 0)); });
    hi = Math.ceil(hi * 20) / 20;
    const bw = W - mx - right;
    const X = v => mx + bw * (v + hi) / (2 * hi);
    let s = '';
    s += line(X(0), top - 6, X(0), top + rows.length * rowH - 4, V('border-strong'), 1.3);
    rows.forEach((r, i) => {
      const y = top + i * rowH + rowH / 2 - 3;
      const lo = r.D - 1.96 * (r.se || 0), up = r.D + 1.96 * (r.se || 0);
      const sig = Math.abs(r.z) > 1.96;
      s += line(X(lo), y, X(up), y, V('border-strong'), 1.3);
      s += line(X(lo), y - 4, X(lo), y + 4, V('border-strong'), 1.1);
      s += line(X(up), y - 4, X(up), y + 4, V('border-strong'), 1.1);
      s += circ(X(r.D), y, 4, sig ? V('c2') : V('c10'), sig ? '' : 'opacity="0.7"');
      s += txt(mx - 8, y + 4, r.label, 'art-txt', 10, 'end');
      s += txt(mx + bw + 6, y + 4, `${r.D >= 0 ? '' : '−'}${Math.abs(r.D).toFixed(3)}  z ${r.z.toFixed(2)}`,
        sig ? 'art-txt' : 'art-mut', 9.5);
    });
    const step = niceStep(2 * hi, 6);
    for (let v = -hi; v <= hi + 1e-9; v += step) {
      s += line(X(v), top + rows.length * rowH - 4, X(v), top + rows.length * rowH, V('border-strong'), 1);
      s += txt(X(v), top + rows.length * rowH + 13, fmtTick(v, step), 'art-mut', 9, 'middle');
    }
    s += txt(mx + bw / 2, H - 4, opts.xLabel || 'D', 'art-mut', 10, 'middle');
    return svg(`0 0 ${W} ${H}`, s);
  }

  Object.assign(Plots11, { tanglegram, concordance, concordanceScatter, network, distanceHeat, dBars, layout });
  window.Plots11 = Plots11;
})();
