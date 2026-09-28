/* PhylogenyPro — molecular phylogenetics in the browser.
   Copyright (C) 2026  Luis Ángel Barrera-Guzmán

   This program is free software: you can redistribute it and/or modify it under
   the terms of the GNU General Public License as published by the Free Software
   Foundation, either version 3 of the License, or (at your option) any later
   version. It is distributed in the hope that it will be useful, but WITHOUT ANY
   WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
   A PARTICULAR PURPOSE. See the GNU General Public License, in the file LICENSE
   at the root of this program, or <https://www.gnu.org/licenses/>. */

/* PhylogenyPro — hand-drawn SVG illustrations.
   Every picture is generated here with CSS-variable colours, so the whole app
   follows the light/dark theme and nothing depends on external images.
   Each function returns an SVG string.

   Most pictures are built from a handful of primitives — a tree, an alignment,
   a curve, a grid, a map — parameterised through the SPECS table at the end.
   That keeps forty method cards honest: each one really draws its own thing. */

(function () {

  const V = n => `var(--${n})`;
  const f1 = v => (+v).toFixed(1);
  const R = seed => rng(seed || 7);
  function wrap(vb, inner, extra) { return `<svg viewBox="${vb}" xmlns="http://www.w3.org/2000/svg" ${extra || ''}>${inner}</svg>`; }
  function txt(x, y, s, cls, size, anchor, extra) {
    return `<text x="${f1(x)}" y="${f1(y)}" class="${cls || 'art-mut'}" font-size="${size || 8}" text-anchor="${anchor || 'start'}" ${extra || ''}>${s}</text>`;
  }
  const line = (x1, y1, x2, y2, stroke, w, extra) => `<line x1="${f1(x1)}" y1="${f1(y1)}" x2="${f1(x2)}" y2="${f1(y2)}" stroke="${stroke}" stroke-width="${w || 1}" stroke-linecap="round" ${extra || ''}/>`;
  const circ = (cx, cy, r, fill, extra) => `<circle cx="${f1(cx)}" cy="${f1(cy)}" r="${f1(r)}" fill="${fill}" ${extra || ''}/>`;
  const rect = (x, y, w, h, fill, extra) => `<rect x="${f1(x)}" y="${f1(y)}" width="${f1(Math.max(0, w))}" height="${f1(Math.max(0, h))}" fill="${fill}" ${extra || ''}/>`;
  const path = (d, stroke, w, extra) => `<path d="${d}" stroke="${stroke}" stroke-width="${w || 1}" fill="none" stroke-linecap="round" stroke-linejoin="round" ${extra || ''}/>`;
  const poly = pts => pts.map((p, i) => (i ? 'L' : 'M') + f1(p[0]) + ' ' + f1(p[1])).join(' ');
  const BASE_V = ['base-a', 'base-c', 'base-g', 'base-t'];

  /* ============================================================
     TREES — the primitive everything else leans on
     ============================================================ */

  /* A random rooted tree of n tips: every step joins two random lineages,
     which is the coalescent read backwards and gives balanced-looking trees. */
  function randomTree(n, r, spread) {
    spread = spread == null ? 0.55 : spread;
    let pool = Array.from({ length: n }, (_, i) => ({ tip: i, depth: 0, children: [] }));
    while (pool.length > 1) {
      const i = Math.floor(r() * pool.length);
      let j = Math.floor(r() * (pool.length - 1)); if (j >= i) j++;
      const a = pool[i], b = pool[j];
      const d = Math.max(a.depth, b.depth) + 0.18 + spread * r();
      const node = { depth: d, children: [a, b] };
      pool = pool.filter((_, k) => k !== i && k !== j);
      pool.push(node);
    }
    return pool[0];
  }
  /* Ultrametric version: internal nodes keep their order in time and every tip
     ends at the same depth — what a dated tree looks like. */
  function ultrametric(root, r) {
    const md = maxDepth(root) || 1;
    (function walk(n, parentDepth) {
      if (!n.children.length) { n.depth = md; return; }
      if (n !== root) n.depth = Math.min(md * 0.92, parentDepth + (md - parentDepth) * (0.25 + 0.45 * r()));
      n.children.forEach(c => walk(c, n.depth));
    })(root, 0);
    return root;
  }
  /* tips in drawing order, each with its y */
  function layout(root, y0, y1) {
    const tips = [];
    (function collect(n) { if (!n.children.length) tips.push(n); else n.children.forEach(collect); })(root);
    tips.forEach((t, i) => { t.y = tips.length === 1 ? (y0 + y1) / 2 : y0 + (y1 - y0) * i / (tips.length - 1); });
    (function setY(n) {
      if (!n.children.length) return n.y;
      const ys = n.children.map(setY);
      n.y = (Math.min(...ys) + Math.max(...ys)) / 2;
      return n.y;
    })(root);
    return tips;
  }
  function maxDepth(n) { return n.children.length ? Math.max(...n.children.map(maxDepth)) : n.depth; }

  /* Rectangular phylogram. opts:
     n tips · seed · ultra (ultrametric) · support (dots on nodes) · clades (colour groups)
     axis ('time' draws a time scale) · hpd (bars on nodes) · outgroup (marks the first tip)
     tipDots · labels (array) · w/h */
  function tree(opts) {
    opts = opts || {};
    const W = opts.W || 200, H = opts.H || 124;
    const n = opts.n || 8, r = R(opts.seed);
    const root = randomTree(n, r, opts.spread);
    if (opts.ultra) ultrametric(root, r);
    const x0 = opts.x0 == null ? 8 : opts.x0, x1 = W - (opts.x1 == null ? 22 : opts.x1);
    const md = maxDepth(root) || 1;
    const X = d => x0 + (x1 - x0) * d / md;
    const tips = layout(root, opts.y0 == null ? 12 : opts.y0, opts.y1 == null ? H - 16 : opts.y1);
    /* tip index, so a clade colour can be looked up, and the clade of a node
       when every tip below it agrees */
    tips.forEach((t, i) => { t.ti = i; });
    const cladeOf = nd => {
      if (!opts.clades) return null;
      const seen = new Set();
      (function walk(x) { if (!x.children.length) seen.add(opts.clades[x.ti]); else x.children.forEach(walk); })(nd);
      return seen.size === 1 ? [...seen][0] : null;
    };
    let s = '';
    if (opts.axis === 'time') {
      for (let i = 0; i < 4; i++) {
        const xa = X(md * i / 4), xb = X(md * (i + 1) / 4);
        if (i % 2 === 0) s += rect(xa, 6, xb - xa, H - 20, V('time-band'));
      }
      s += line(x0, H - 12, x1, H - 12, V('border-strong'), 1);
      for (let i = 0; i <= 4; i++) s += line(X(md * i / 4), H - 12, X(md * i / 4), H - 9, V('border-strong'), 1);
    }
    const base = opts.stroke || V('primary');
    const colour = nd => { const g = cladeOf(nd); return g == null ? base : V('c' + (1 + (g % 10))); };
    (function draw(nd, parentX) {
      const x = X(nd.depth), col = colour(nd);
      s += line(parentX, nd.y, x, nd.y, col, opts.lw || 1.6);
      if (nd.children.length) {
        const ys = nd.children.map(c => c.y);
        s += line(x, Math.min(...ys), x, Math.max(...ys), col, opts.lw || 1.6);
        nd.children.forEach(c => draw(c, x));
      } else {
        if (opts.tipDots !== false) s += circ(x, nd.y, opts.tipR || 2.2, col);
        if (opts.labels && opts.labels[nd.ti]) s += txt(x + 4, nd.y + 2.6, opts.labels[nd.ti], 'art-mut', opts.fs || 7);
      }
    })(root, x0);
    if (opts.support) {
      (function dots(nd) {
        if (!nd.children.length) return;
        const x = X(nd.depth);
        if (nd !== root) {
          const v = r();
          s += circ(x, nd.y, 2.6, v > 0.45 ? V('leaf') : v > 0.2 ? V('gold') : V('rose'), 'opacity="0.95"');
        }
        nd.children.forEach(dots);
      })(root);
    }
    if (opts.hpd) {
      (function bars(nd) {
        if (!nd.children.length) return;
        const x = X(nd.depth), w = (x1 - x0) * 0.07 * (0.6 + r());
        s += rect(x - w, nd.y - 2.2, 2 * w, 4.4, V('accent'), 'opacity="0.3" rx="2"');
        nd.children.forEach(bars);
      })(root);
    }
    if (opts.outgroup && tips.length) {
      const t = tips[tips.length - 1];
      s += line(X(t.depth) + 5, t.y - 5, X(t.depth) + 5, t.y + 5, V('accent'), 2);
    }
    if (opts.fossil) {
      const nodes = [];
      (function all(nd) { if (nd.children.length) nodes.push(nd); nd.children.forEach(all); })(root);
      const pick = nodes[Math.floor(nodes.length / 2)] || root;
      s += `<path d="M${f1(X(pick.depth))} ${f1(pick.y - 4)} l3.4 4 -3.4 4 -3.4 -4 Z" fill="${V('accent')}"/>`;
    }
    return opts.raw ? s : wrap(`0 0 ${W} ${H}`, s);
  }

  /* Circular tree */
  function treeCircular(opts) {
    opts = opts || {};
    const W = opts.W || 190, H = opts.H || 124, cx = W / 2, cy = H / 2;
    const n = opts.n || 14, r = R(opts.seed), root = randomTree(n, r, 0.5);
    const md = maxDepth(root) || 1, R0 = 8, R1 = Math.min(W, H) / 2 - 12;
    const open = opts.open || 0.35;
    const tips = [];
    (function collect(nd) { if (!nd.children.length) tips.push(nd); else nd.children.forEach(collect); })(root);
    tips.forEach((t, i) => { t.a = -Math.PI / 2 + (2 * Math.PI - open) * i / (tips.length - 1 || 1) + open / 2; });
    (function setA(nd) {
      if (!nd.children.length) return nd.a;
      const as = nd.children.map(setA);
      nd.a = (Math.min(...as) + Math.max(...as)) / 2;
      return nd.a;
    })(root);
    const rad = d => R0 + (R1 - R0) * d / md;
    const P = (rr, a) => [cx + rr * Math.cos(a), cy + rr * Math.sin(a)];
    let s = '';
    let ti = 0;
    (function draw(nd, pr) {
      const rr = rad(nd.depth);
      const col = opts.clades ? V('c' + (1 + (opts.clades[Math.min(ti, opts.clades.length - 1)] % 10))) : V('primary');
      const [ax, ay] = P(pr, nd.a), [bx, by] = P(rr, nd.a);
      s += line(ax, ay, bx, by, col, 1.4);
      if (nd.children.length) {
        const as = nd.children.map(c => c.a);
        const [sx, sy] = P(rr, Math.min(...as)), [ex, ey] = P(rr, Math.max(...as));
        s += `<path d="M${f1(sx)} ${f1(sy)} A ${f1(rr)} ${f1(rr)} 0 ${Math.abs(Math.max(...as) - Math.min(...as)) > Math.PI ? 1 : 0} 1 ${f1(ex)} ${f1(ey)}" stroke="${col}" stroke-width="1.4" fill="none"/>`;
        nd.children.forEach(c => draw(c, rr));
      } else { s += circ(bx, by, 1.8, col); ti++; }
    })(root, R0 * 0.4);
    if (opts.ring) {
      tips.forEach((t, i) => {
        const [x1, y1] = P(R1 + 3, t.a - 0.05), [x2, y2] = P(R1 + 7, t.a - 0.05);
        s += line(x1, y1, x2, y2, V('c' + (1 + (i % 10))), 3.4);
      });
    }
    return wrap(`0 0 ${W} ${H}`, s);
  }

  /* Unrooted tree, drawn by the equal-angle algorithm */
  function unrooted(opts) {
    opts = opts || {};
    const W = opts.W || 190, H = opts.H || 124;
    const n = opts.n || 10, r = R(opts.seed), root = randomTree(n, r, 0.5);
    let s = '';
    const nTips = n;
    let used = 0;
    (function draw(nd, x, y, a0, a1) {
      const leaves = countTips(nd);
      const a = (a0 + a1) / 2;
      const len = (nd.children.length ? 12 : 16) * (0.6 + r() * 0.9);
      const x2 = x + len * Math.cos(a), y2 = y + len * Math.sin(a);
      s += line(x, y, x2, y2, opts.stroke || V('primary'), 1.5);
      if (!nd.children.length) { s += circ(x2, y2, 2, V('c' + (1 + (used++ % 10)))); return; }
      let acc = a0;
      nd.children.forEach(c => {
        const share = (a1 - a0) * countTips(c) / leaves;
        draw(c, x2, y2, acc, acc + share);
        acc += share;
      });
    })(root, W / 2, H / 2, 0, Math.PI * 2);
    return wrap(`0 0 ${W} ${H}`, s);
  }
  function countTips(n) { return n.children.length ? n.children.reduce((a, c) => a + countTips(c), 0) : 1; }

  /* ============================================================
     ALIGNMENT
     ============================================================ */
  /* opts: rows, cols, gaps, highlight ('informative'|'variable'|null),
     protein (grey-blue tones), trimmed (shows the blocks that would be cut) */
  function align(opts) {
    opts = opts || {};
    const W = opts.W || 200, H = opts.H || 124;
    const rows = opts.rows || 8, cols = opts.cols || 26, r = R(opts.seed);
    const x0 = opts.labels === false ? 4 : 18, y0 = 10;
    const cw = (W - x0 - 6) / cols, ch = Math.min(9, (H - y0 - 12) / rows);
    let s = '';
    /* a column is either conserved (one base) or variable */
    const colType = [];
    for (let c = 0; c < cols; c++) colType.push(r() < (opts.varFrac || 0.3) ? (r() < 0.55 ? 2 : 1) : 0);
    for (let c = 0; c < cols; c++) {
      const base0 = Math.floor(r() * 4), alt = (base0 + 1 + Math.floor(r() * 3)) % 4;
      for (let i = 0; i < rows; i++) {
        const gap = opts.gaps && r() < 0.05;
        let b = base0;
        if (colType[c] === 1 && i === Math.floor(r() * rows)) b = alt;          // singleton
        else if (colType[c] === 2 && r() < 0.45) b = alt;                        // informative
        const fill = gap ? V('base-gap') : (opts.protein ? V('c' + (1 + (b * 2 % 10))) : V(BASE_V[b]));
        s += rect(x0 + c * cw, y0 + i * ch, cw - 0.6, ch - 0.6, fill, `opacity="${gap ? 0.45 : 0.9}" rx="0.8"`);
      }
      if (opts.highlight === 'informative' && colType[c] === 2) s += rect(x0 + c * cw - 0.4, y0 - 4, cw, 3, V('accent'), 'rx="1"');
      if (opts.highlight === 'variable' && colType[c] > 0) s += rect(x0 + c * cw - 0.4, y0 - 4, cw, 3, V('primary'), 'rx="1"');
    }
    if (opts.labels !== false) for (let i = 0; i < rows; i++) s += rect(2, y0 + i * ch + 1.4, 13, ch - 3.4, V('border'), 'rx="1.4"');
    if (opts.trimmed) {
      const a = Math.floor(cols * 0.18), b = Math.floor(cols * 0.32);
      s += rect(x0 + a * cw, y0 - 5, (b - a) * cw, rows * ch + 8, V('danger'), 'opacity="0.16" rx="2"');
      const c2 = Math.floor(cols * 0.72), d2 = Math.floor(cols * 0.85);
      s += rect(x0 + c2 * cw, y0 - 5, (d2 - c2) * cw, rows * ch + 8, V('danger'), 'opacity="0.16" rx="2"');
      s += txt(W / 2, H - 2, opts.trimLabel || '', 'art-mut', 7, 'middle');
    }
    if (opts.partitions) {
      const cuts = [0, Math.floor(cols * 0.38), Math.floor(cols * 0.68), cols];
      const names = opts.partNames || ['', '', ''];
      for (let k = 0; k < 3; k++) {
        s += rect(x0 + cuts[k] * cw, H - 10, (cuts[k + 1] - cuts[k]) * cw - 1.5, 5, V('c' + (k * 2 + 1)), 'rx="1.6"');
        if (names[k]) s += txt(x0 + (cuts[k] + cuts[k + 1]) / 2 * cw, H - 12, names[k], 'art-mut', 6.4, 'middle');
      }
    }
    return wrap(`0 0 ${W} ${H}`, s);
  }

  /* ============================================================
     CURVES, BARS, GRIDS
     ============================================================ */
  /* opts.series: [{f(t)->y in 0..1, colour, dash}] · marks · shade */
  function curves(opts) {
    opts = opts || {};
    const W = opts.W || 200, H = opts.H || 124;
    const mx = 20, my = 14;
    let s = '';
    s += line(mx, H - my, W - 6, H - my, V('border-strong'), 1);
    s += line(mx, 6, mx, H - my, V('border-strong'), 1);
    const px = t => mx + (W - mx - 8) * t, py = v => (H - my) - (H - my - 8) * Math.max(0, Math.min(1, v));
    if (opts.shade) {
      const pts = [];
      for (let t = 0; t <= 1.001; t += 0.02) pts.push([px(t), py(opts.shade(t))]);
      s += `<path d="${poly(pts)} L${f1(px(1))} ${f1(py(0))} L${f1(px(0))} ${f1(py(0))} Z" fill="${V('primary')}" opacity="0.12"/>`;
    }
    (opts.series || []).forEach(se => {
      const pts = [];
      for (let t = 0; t <= 1.001; t += 0.02) pts.push([px(t), py(se.f(t))]);
      s += path(poly(pts), se.c || V('primary'), se.w || 1.8, se.dash ? `stroke-dasharray="${se.dash}"` : '');
    });
    (opts.points || []).forEach(p => { s += circ(px(p[0]), py(p[1]), p[2] || 1.8, p[3] || V('accent')); });
    if (opts.xlab) s += txt(W / 2, H - 3, opts.xlab, 'art-mut', 7, 'middle');
    if (opts.ylab) s += txt(-H / 2, 8, opts.ylab, 'art-mut', 7, 'middle', 'transform="rotate(-90)"');
    return wrap(`0 0 ${W} ${H}`, s);
  }

  function bars(opts) {
    opts = opts || {};
    const W = opts.W || 200, H = opts.H || 124;
    const vals = opts.values || [0.9, 0.7, 0.55, 0.45, 0.3, 0.2];
    const mx = 16, my = 14, bw = (W - mx - 8) / vals.length;
    let s = line(mx, H - my, W - 6, H - my, V('border-strong'), 1);
    vals.forEach((v, i) => {
      const h = (H - my - 10) * v;
      const col = opts.colours ? V(opts.colours[i % opts.colours.length]) : (i === (opts.best == null ? 0 : opts.best) ? V('accent') : V('primary'));
      s += rect(mx + i * bw + bw * 0.16, H - my - h, bw * 0.68, h, col, 'rx="1.6" opacity="0.92"');
    });
    if (opts.line != null) s += line(mx, H - my - (H - my - 10) * opts.line, W - 6, H - my - (H - my - 10) * opts.line, V('danger'), 1.2, 'stroke-dasharray="3 2"');
    if (opts.xlab) s += txt(W / 2, H - 3, opts.xlab, 'art-mut', 7, 'middle');
    return wrap(`0 0 ${W} ${H}`, s);
  }

  /* a k×k matrix: substitution matrix Q, distance matrix, relationship matrix */
  function grid(opts) {
    opts = opts || {};
    const W = opts.W || 190, H = opts.H || 124, k = opts.k || 4, r = R(opts.seed);
    const cell = Math.min((W - 40) / k, (H - 30) / k);
    const x0 = (W - cell * k) / 2, y0 = (H - cell * k) / 2 + 4;
    let s = '';
    for (let i = 0; i < k; i++) for (let j = 0; j < k; j++) {
      const diag = i === j;
      const v = diag ? 1 : (opts.f ? opts.f(i, j) : r());
      const fill = diag && opts.diagDark ? V('primary') : V('primary');
      s += rect(x0 + j * cell, y0 + i * cell, cell - 1, cell - 1, fill, `opacity="${(0.12 + 0.75 * v).toFixed(2)}" rx="1.4"`);
    }
    if (opts.labels) opts.labels.forEach((L, i) => {
      s += txt(x0 - 3, y0 + i * cell + cell * 0.68, L, 'art-mut', 7, 'end');
      s += txt(x0 + i * cell + cell / 2, y0 - 3, L, 'art-mut', 7, 'middle');
    });
    return wrap(`0 0 ${W} ${H}`, s);
  }

  /* ============================================================
     MAPS AND NETWORKS
     ============================================================ */
  /* a schematic map of areas with dispersal arrows — historical biogeography */
  function mapAreas(opts) {
    opts = opts || {};
    const W = opts.W || 200, H = opts.H || 124, r = R(opts.seed || 3);
    const areas = [
      { x: 44, y: 38, rx: 26, ry: 18, c: 'c1' }, { x: 118, y: 30, rx: 22, ry: 14, c: 'c3' },
      { x: 150, y: 76, rx: 24, ry: 16, c: 'c2' }, { x: 62, y: 88, rx: 28, ry: 15, c: 'c5' },
    ];
    let s = rect(0, 0, W, H, V('bg-soft'), 'rx="6"');
    areas.forEach((a, i) => {
      s += `<ellipse cx="${a.x}" cy="${a.y}" rx="${a.rx}" ry="${a.ry}" fill="${V(a.c)}" opacity="0.75"/>`;
      if (opts.letters) s += txt(a.x, a.y + 3, 'ABCD'[i], 'art-txt', 9, 'middle', `fill="${V('card-bg')}" font-weight="700"`);
    });
    if (opts.arrows !== false) {
      const arrow = (a, b, col) => {
        const mx = (a.x + b.x) / 2, my = (a.y + b.y) / 2 - 16;
        return path(`M${a.x} ${a.y} Q ${mx} ${my} ${b.x} ${b.y}`, col, 1.6, 'stroke-dasharray="4 3" marker-end="url(#phyArrow)"');
      };
      s = `<defs><marker id="phyArrow" markerWidth="6" markerHeight="6" refX="5" refY="3" orient="auto"><path d="M0 0 L6 3 L0 6 z" fill="${V('text-muted')}"/></marker></defs>` + s;
      s += arrow(areas[0], areas[1], V('text-muted'));
      s += arrow(areas[1], areas[2], V('text-muted'));
    }
    if (opts.pie) areas.forEach(a => {
      s += circ(a.x, a.y, 7, V('card-bg'), 'opacity="0.9"');
      s += `<path d="M${a.x} ${a.y} L${a.x} ${a.y - 7} A 7 7 0 0 1 ${f1(a.x + 7 * Math.sin(2.2))} ${f1(a.y - 7 * Math.cos(2.2))} Z" fill="${V('accent')}"/>`;
    });
    return wrap(`0 0 ${W} ${H}`, s);
  }

  /* split network / haplotype-style network */
  function network(opts) {
    opts = opts || {};
    const W = opts.W || 190, H = opts.H || 124, r = R(opts.seed || 11);
    const nodes = [];
    const k = opts.n || 9;
    for (let i = 0; i < k; i++) {
      const a = 2 * Math.PI * i / k + r() * 0.3;
      const rr = 26 + r() * 18;
      nodes.push([W / 2 + rr * Math.cos(a), H / 2 + rr * Math.sin(a) * 0.82]);
    }
    let s = '';
    for (let i = 0; i < k; i++) {
      s += line(W / 2, H / 2, nodes[i][0], nodes[i][1], V('border-strong'), 1.2);
      if (opts.boxes && i < k - 1) s += line(nodes[i][0], nodes[i][1], nodes[i + 1][0], nodes[i + 1][1], V('border-strong'), 1);
    }
    s += circ(W / 2, H / 2, 6, V('primary'));
    nodes.forEach((n2, i) => s += circ(n2[0], n2[1], 3 + r() * 3, V('c' + (1 + i % 10)), 'opacity="0.9"'));
    return wrap(`0 0 ${W} ${H}`, s);
  }

  /* two trees facing each other with links: tanglegram */
  function tanglegram(opts) {
    opts = opts || {};
    const W = 200, H = 124, n = 7;
    const r = R(opts.seed || 5);
    const left = tree({ n, seed: opts.seed || 5, W: 92, H, raw: true, x0: 4, x1: 12, tipDots: true, stroke: V('primary') });
    let s = `<g>${left}</g>`;
    s += `<g transform="translate(${W},0) scale(-1,1)">${tree({ n, seed: (opts.seed || 5) + 9, W: 92, H, raw: true, x0: 4, x1: 12, tipDots: true, stroke: V('leaf') })}</g>`;
    const ys = [];
    for (let i = 0; i < n; i++) ys.push(12 + (H - 28) * i / (n - 1));
    const perm = ys.slice(); shuffle(perm, r);
    ys.forEach((y, i) => { s += line(84, y, W - 84, perm[i], V('text-muted'), 0.9, 'opacity="0.6"'); });
    return wrap(`0 0 ${W} ${H}`, s);
  }

  /* ============================================================
     HERO — an alignment turns into a tree, and the tree into time
     ============================================================ */
  function hero() {
    const W = 540, H = 430;
    let s = `<defs>
      <linearGradient id="phFade" x1="0" y1="0" x2="1" y2="0">
        <stop offset="0" stop-color="${V('primary')}" stop-opacity="0.16"/>
        <stop offset="1" stop-color="${V('primary')}" stop-opacity="0"/>
      </linearGradient>
      <linearGradient id="phAmber" x1="0" y1="0" x2="0" y2="1">
        <stop offset="0" stop-color="${V('accent')}" stop-opacity="0.35"/>
        <stop offset="1" stop-color="${V('accent')}" stop-opacity="0.05"/>
      </linearGradient>
    </defs>`;

    /* --- 1. the alignment, at the left --- */
    const r = R(20260922);
    const rows = 9, cols = 13, cw = 11, ch = 13, ax = 14, ay = 44;
    for (let c = 0; c < cols; c++) {
      const b0 = Math.floor(r() * 4), alt = (b0 + 1 + Math.floor(r() * 3)) % 4;
      const informative = r() < 0.32;
      for (let i = 0; i < rows; i++) {
        const b = informative && r() < 0.45 ? alt : b0;
        s += rect(ax + c * cw, ay + i * ch, cw - 1.4, ch - 1.4, V(BASE_V[b]),
          `rx="1.4" opacity="0.88" class="${informative && c % 3 === 0 ? 'site-pulse' : ''}"`);
      }
      if (informative) s += rect(ax + c * cw, ay - 7, cw - 1.4, 3.4, V('accent'), 'rx="1.4"');
    }
    s += txt(ax, ay - 14, '', 'art-mut', 9);
    for (let i = 0; i < rows; i++) s += rect(2, ay + i * ch + 2, 9, ch - 5, V('border'), 'rx="1.2"');
    s += rect(ax - 6, ay - 12, cols * cw + 10, rows * ch + 18, 'none', `stroke="${V('border')}" rx="6"`);

    /* --- 2. the tree, in the middle --- */
    const tx0 = 182, tx1 = 392, ty0 = 30, ty1 = H - 96;
    const root = randomTree(9, R(77), 0.5);
    const md = maxDepth(root) || 1;
    const X = d => tx0 + (tx1 - tx0) * d / md;
    const tips = layout(root, ty0, ty1);
    const clade = i => V('c' + (1 + (i < 3 ? 0 : i < 6 ? 2 : 4)));
    let ti = 0;
    (function draw(nd, px) {
      const x = X(nd.depth);
      const col = nd.children.length ? V('primary') : clade(ti);
      s += line(px, nd.y, x, nd.y, col, nd.children.length ? 2.6 : 2.6);
      if (nd.children.length) {
        const ys = nd.children.map(c2 => c2.y);
        s += line(x, Math.min(...ys), x, Math.max(...ys), V('primary'), 2.6);
        s += circ(x, nd.y, 3.4, V('card-bg'), `stroke="${V('primary')}" stroke-width="1.6"`);
        nd.children.forEach(c2 => draw(c2, x));
      } else {
        s += circ(x, nd.y, 4.2, clade(ti));
        s += rect(x + 8, nd.y - 3, 26 + 14 * ((ti * 7) % 3), 6, V('border'), 'rx="3" opacity="0.8"');
        ti++;
      }
    })(root, tx0 - 16, null);
    /* the dashed path from the alignment into the tree */
    s += path(`M${ax + cols * cw + 6} ${ay + rows * ch / 2} C 150 ${ay + rows * ch / 2}, 150 ${(ty0 + ty1) / 2}, ${tx0 - 18} ${(ty0 + ty1) / 2}`,
      V('accent'), 2, 'stroke-dasharray="5 4" class="branch-grow" opacity="0.85"');

    /* --- 3. the time scale below the tree --- */
    const gy = H - 74, gh = 20;
    const bands = [['', 0, 0.28, 'c5'], ['', 0.28, 0.55, 'c3'], ['', 0.55, 0.78, 'c2'], ['', 0.78, 1, 'c6']];
    bands.forEach(b => {
      s += rect(tx0 + (tx1 - tx0) * b[1], gy, (tx1 - tx0) * (b[2] - b[1]) - 1.5, gh, V(b[3]), 'rx="3" opacity="0.55"');
    });
    s += line(tx0, gy - 6, tx1, gy - 6, V('border-strong'), 1.4);
    for (let i = 0; i <= 4; i++) s += line(tx0 + (tx1 - tx0) * i / 4, gy - 6, tx0 + (tx1 - tx0) * i / 4, gy - 11, V('border-strong'), 1.4);
    /* a fossil calibration sitting on a node */
    const fx = tx0 + (tx1 - tx0) * 0.42, fy = gy - 24;
    s += rect(fx - 26, fy - 14, 52, 20, V('card-bg'), `rx="5" stroke="${V('accent')}"`);
    s += `<path d="M${fx - 14} ${fy - 4} l4.6 5.4 -4.6 5.4 -4.6 -5.4 Z" fill="${V('accent')}"/>`;
    s += rect(fx - 6, fy - 7, 26, 3.4, V('accent'), 'rx="1.7" opacity="0.5"');
    s += rect(fx - 6, fy - 1, 18, 3.4, V('accent'), 'rx="1.7" opacity="0.35"');

    /* --- 4. posterior density card, top right --- */
    const px0 = 412, py0 = 26, pw = 116, ph = 74;
    s += rect(px0, py0, pw, ph, V('card-bg'), `rx="8" stroke="${V('border')}"`);
    const dens = [];
    for (let t = 0; t <= 1.0001; t += 0.02) {
      const z = (t - 0.46) / 0.17;
      dens.push([px0 + 10 + (pw - 20) * t, py0 + ph - 16 - (ph - 30) * Math.exp(-0.5 * z * z)]);
    }
    s += `<path d="${poly(dens)} L${f1(px0 + pw - 10)} ${f1(py0 + ph - 16)} L${f1(px0 + 10)} ${f1(py0 + ph - 16)} Z" fill="url(#phAmber)"/>`;
    s += path(poly(dens), V('accent'), 2);
    s += line(px0 + 10, py0 + ph - 16, px0 + pw - 10, py0 + ph - 16, V('border-strong'), 1);
    s += rect(px0 + 10 + (pw - 20) * 0.24, py0 + ph - 13, (pw - 20) * 0.45, 3.4, V('accent'), 'rx="1.7"');

    /* --- 5. support card, bottom right --- */
    const sx = 412, sy = 116, sw = 116, sh = 92;
    s += rect(sx, sy, sw, sh, V('card-bg'), `rx="8" stroke="${V('border')}"`);
    [[0.98, 'leaf'], [0.86, 'leaf'], [0.71, 'gold'], [0.55, 'gold'], [0.34, 'rose']].forEach((b, i) => {
      s += rect(sx + 12, sy + 14 + i * 15, (sw - 34) * b[0], 8, V(b[1]), 'rx="4" opacity="0.85"');
      s += circ(sx + sw - 14, sy + 18 + i * 15, 3, V(b[1]));
    });

    /* --- 6. ancestral areas strip, bottom --- */
    const my0 = H - 44;
    for (let i = 0; i < 4; i++) {
      s += `<ellipse cx="${f1(200 + i * 46)}" cy="${f1(my0 + 16)}" rx="20" ry="12" fill="${V('c' + (1 + i * 2))}" opacity="0.55"/>`;
      s += txt(200 + i * 46, my0 + 19, 'ABCD'[i], 'art-txt', 9, 'middle', `fill="${V('card-bg')}" font-weight="700"`);
    }
    for (let i = 0; i < 3; i++) s += path(`M${f1(220 + i * 46)} ${f1(my0 + 10)} Q ${f1(223 + i * 46)} ${f1(my0)} ${f1(226 + i * 46)} ${f1(my0 + 10)}`, V('text-muted'), 1.4, 'stroke-dasharray="3 2"');

    return wrap(`0 0 ${W} ${H}`, s, 'style="overflow:visible"');
  }

  /* ============================================================
     THEORY FIGURES (they carry translated labels, so they are redrawn
     whenever the language changes)
     ============================================================ */
  function theoryAnatomy() {
    const W = 320, H = 200;
    const L = (es, en) => (window.I18N && I18N.lang === 'en' ? en : es);
    let s = '';
    const tips = [30, 62, 94, 130, 166];
    const X = { root: 24, n1: 70, n2: 120, n3: 96, tip: 236 };
    s += line(10, 98, X.root, 98, V('text-muted'), 2.2);
    s += line(X.root, 30, X.root, 166, V('primary'), 2.2);
    s += line(X.root, 98, X.root, 98, V('primary'), 2.2);
    /* upper clade */
    s += line(X.root, 30, X.n1, 30, V('primary'), 2.2);
    s += line(X.n1, 16, X.n1, 46, V('c1'), 2.2);
    s += line(X.n1, 16, X.tip, 16, V('c1'), 2.2);
    s += line(X.n1, 46, X.tip, 46, V('c1'), 2.2);
    /* middle clade */
    s += line(X.root, 98, X.n2, 98, V('primary'), 2.2);
    s += line(X.n2, 80, X.n2, 118, V('c3'), 2.2);
    s += line(X.n2, 80, X.tip, 80, V('c3'), 2.2);
    s += line(X.n2, 118, X.tip, 118, V('c3'), 2.2);
    /* outgroup */
    s += line(X.root, 166, X.tip, 166, V('accent'), 2.2);
    [16, 46, 80, 118, 166].forEach((y, i) => s += circ(X.tip, y, 3.4, i === 4 ? V('accent') : (y < 60 ? V('c1') : V('c3'))));
    s += circ(X.root, 98, 4, V('card-bg'), `stroke="${V('primary')}" stroke-width="2"`);
    s += circ(X.n1, 30, 3.4, V('card-bg'), `stroke="${V('c1')}" stroke-width="1.8"`);
    s += circ(X.n2, 98, 3.4, V('card-bg'), `stroke="${V('c3')}" stroke-width="1.8"`);
    s += txt(X.root - 6, 92, L('raíz', 'root'), 'art-mut', 9, 'end');
    s += txt(X.n1 + 4, 26, L('nodo', 'node'), 'art-mut', 9);
    s += txt(X.tip + 8, 19, L('punta (taxón)', 'tip (taxon)'), 'art-mut', 9);
    s += txt(X.tip + 8, 169, L('grupo externo', 'outgroup'), 'art-mut', 9, 'start', `fill="${V('accent')}"`);
    s += txt(150, 60, L('rama', 'branch'), 'art-mut', 9, 'middle');
    s += path('M150 64 L150 78', V('text-muted'), 1, 'stroke-dasharray="2 2"');
    s += `<path d="M${X.n1 - 4} 12 L${X.n1 - 4} 50 L${X.tip + 20} 50 L${X.tip + 20} 12 Z" fill="none" stroke="${V('c1')}" stroke-dasharray="3 3" opacity="0.6" rx="4"/>`;
    s += txt(X.tip + 24, 34, L('clado', 'clade'), 'art-mut', 9, 'start', `fill="${V('c1')}"`);
    return wrap(`0 0 ${W} ${H}`, s);
  }

  function theorySaturation() {
    const W = 320, H = 200;
    const L = (es, en) => (window.I18N && I18N.lang === 'en' ? en : es);
    const mx = 40, my = 30;
    let s = line(mx, H - my, W - 14, H - my, V('border-strong'), 1.4) + line(mx, 14, mx, H - my, V('border-strong'), 1.4);
    const px = t => mx + (W - mx - 18) * t, py = v => (H - my) - (H - my - 20) * v;
    /* true distance = identity line; observed = saturating curve */
    const obs = [], tru = [];
    for (let t = 0; t <= 1.0001; t += 0.02) {
      tru.push([px(t), py(t)]);
      obs.push([px(t), py(0.75 * (1 - Math.exp(-4 * t / 3 * 1.6)) / 0.75 * 0.74)]);
    }
    s += path(poly(tru), V('text-muted'), 1.6, 'stroke-dasharray="4 3"');
    s += path(poly(obs), V('primary'), 2.4);
    s += line(px(0), py(0.74), px(1), py(0.74), V('danger'), 1.2, 'stroke-dasharray="3 3"');
    s += txt(px(1) - 2, py(0.74) - 5, '3/4', 'art-mut', 8, 'end', `fill="${V('danger')}"`);
    s += txt(px(0.62), py(0.86), L('esperado sin corregir', 'expected if uncorrected'), 'art-mut', 8, 'middle');
    s += txt(W / 2, H - 8, L('distancia real (sustituciones por sitio)', 'true distance (substitutions per site)'), 'art-mut', 8.5, 'middle');
    s += txt(-H / 2, 12, L('diferencias observadas', 'observed differences'), 'art-mut', 8.5, 'middle', 'transform="rotate(-90)"');
    return wrap(`0 0 ${W} ${H}`, s);
  }

  function theoryLongBranch() {
    const W = 320, H = 200;
    const L = (es, en) => (window.I18N && I18N.lang === 'en' ? en : es);
    let s = '';
    const draw = (ox, oy, labels, colours, title) => {
      let t = '';
      t += line(ox + 40, oy + 40, ox + 64, oy + 40, V('text-muted'), 2);       // internal branch
      t += line(ox + 40, oy + 12, ox + 40, oy + 68, V('text-muted'), 2);
      t += line(ox + 64, oy + 12, ox + 64, oy + 68, V('text-muted'), 2);
      t += line(ox + 40, oy + 12, ox + 6, oy + 12, colours[0], 3);             // long
      t += line(ox + 40, oy + 68, ox + 28, oy + 68, colours[1], 3);            // short
      t += line(ox + 64, oy + 12, ox + 98, oy + 12, colours[2], 3);            // long
      t += line(ox + 64, oy + 68, ox + 76, oy + 68, colours[3], 3);            // short
      t += txt(ox + 2, oy + 9, labels[0], 'art-mut', 9, 'end');
      t += txt(ox + 24, oy + 71, labels[1], 'art-mut', 9, 'end');
      t += txt(ox + 102, oy + 15, labels[2], 'art-mut', 9);
      t += txt(ox + 80, oy + 71, labels[3], 'art-mut', 9);
      t += txt(ox + 52, oy + 92, title, 'art-mut', 9, 'middle');
      return t;
    };
    s += draw(16, 20, ['A', 'B', 'C', 'D'], [V('danger'), V('primary'), V('danger'), V('primary')], L('árbol verdadero', 'true tree'));
    s += draw(186, 20, ['A', 'B', 'C', 'D'], [V('danger'), V('primary'), V('primary'), V('danger')], L('lo que recupera la parsimonia', 'what parsimony recovers'));
    s += txt(W / 2, 140, L('Las dos ramas largas (rojas) acumulan tantos cambios que coinciden por azar', 'The two long branches (red) pile up so many changes that they match by chance'), 'art-mut', 8.5, 'middle');
    s += txt(W / 2, 152, L('y la parsimonia las junta: atracción de ramas largas.', 'and parsimony joins them: long-branch attraction.'), 'art-mut', 8.5, 'middle');
    s += txt(W / 2, 172, L('La verosimilitud, que sí cuenta los cambios invisibles, resiste.', 'Likelihood, which does count the invisible changes, resists.'), 'art-mut', 8.5, 'middle', `fill="${V('leaf')}"`);
    return wrap(`0 0 ${W} ${H}`, s);
  }

  function theoryClock() {
    const W = 320, H = 190;
    const L = (es, en) => (window.I18N && I18N.lang === 'en' ? en : es);
    let s = '';
    for (let i = 0; i < 4; i++) s += rect(24 + i * 68, 14, 66, 120, V('time-band'), i % 2 ? 'opacity="0"' : '');
    const nodes = [[24, 74], [92, 42], [160, 100], [206, 26]];
    s += line(24, 26, 24, 122, V('primary'), 2.2);
    s += line(24, 26, 92, 26, V('primary'), 2.2);
    s += line(24, 122, 160, 122, V('primary'), 2.2);
    s += line(92, 14, 92, 42, V('primary'), 2.2);
    s += line(92, 14, 288, 14, V('c1'), 2.2);
    s += line(92, 42, 288, 42, V('c1'), 2.2);
    s += line(160, 100, 160, 134, V('primary'), 2.2);
    s += line(160, 100, 288, 100, V('c3'), 2.2);
    s += line(160, 134, 288, 134, V('c3'), 2.2);
    [14, 42, 100, 134].forEach((y, i) => s += circ(288, y, 3.4, y < 60 ? V('c1') : V('c3')));
    s += line(24, 150, 288, 150, V('border-strong'), 1.4);
    [0, 1, 2, 3, 4].forEach(i => {
      const x = 24 + (288 - 24) * i / 4;
      s += line(x, 150, x, 155, V('border-strong'), 1.4);
      s += txt(x, 165, String(40 - i * 10), 'art-mut', 8, 'middle');
    });
    s += txt(156, 180, L('millones de años antes del presente', 'million years before present'), 'art-mut', 8.5, 'middle');
    s += `<path d="M92 34 l5 6 -5 6 -5 -6 Z" fill="${V('accent')}"/>`;
    s += rect(66, 52, 52, 4, V('accent'), 'rx="2" opacity="0.45"');
    s += txt(92, 66, L('fósil: edad mínima', 'fossil: minimum age'), 'art-mut', 8, 'middle', `fill="${V('accent')}"`);
    return wrap(`0 0 ${W} ${H}`, s);
  }

  function theoryGenomes() {
    const W = 320, H = 180;
    const L = (es, en) => (window.I18N && I18N.lang === 'en' ? en : es);
    let s = '';
    const cell = (x, y, c, title, sub) => {
      let t = rect(x, y, 96, 62, V('card-bg'), `rx="8" stroke="${V('border')}"`);
      t += circ(x + 22, y + 30, 13, V(c), 'opacity="0.8"');
      t += txt(x + 42, y + 26, title, 'art-txt', 9);
      t += txt(x + 42, y + 38, sub, 'art-mut', 7.6);
      return t;
    };
    s += cell(6, 12, 'leaf', L('cloroplasto', 'chloroplast'), L('materno, sin recombinar', 'maternal, no recombination'));
    s += cell(112, 12, 'rose', L('mitocondria', 'mitochondrion'), L('materno, muy lento', 'maternal, very slow'));
    s += cell(218, 12, 'accent', L('ITS / ETS', 'ITS / ETS'), L('nuclear, rápido', 'nuclear, fast'));
    s += cell(60, 92, 'primary', L('genes nucleares', 'nuclear genes'), L('biparental, recombina', 'biparental, recombines'));
    s += cell(166, 92, 'c8', L('codones', 'codons'), L('selección dN/dS', 'dN/dS selection'));
    s += txt(W / 2, 172, L('cada partición declara su genoma y la app lo respeta en todo el análisis', 'each partition declares its genome and the app honours it throughout'), 'art-mut', 8, 'middle');
    return wrap(`0 0 ${W} ${H}`, s);
  }

  /* ============================================================
     SPECS — every named illustration of the app
     ============================================================ */
  const SPECS = {
    /* --- blocks --- */
    blkData: () => align({ rows: 9, cols: 30, gaps: true, highlight: 'variable', seed: 3 }),
    blkModels: () => bars({ values: [1, 0.93, 0.88, 0.7, 0.62, 0.5, 0.34], best: 0, seed: 4 }),
    blkParsimony: () => tree({ n: 9, seed: 5, support: true, spread: 0.7 }),
    blkLikelihood: () => tree({ n: 10, seed: 6, support: true, clades: [0, 0, 0, 2, 2, 2, 4, 4, 6, 6] }),
    blkBayes: () => curves({
      series: [{ f: t => 0.55 + 0.32 * Math.sin(t * 22) * Math.exp(-t * 0.6), c: V('primary'), w: 1.2 },
               { f: t => 0.5 + 0.3 * Math.sin(t * 19 + 1.2) * Math.exp(-t * 0.5), c: V('accent'), w: 1.2 }],
    }),
    blkDating: () => tree({ n: 9, seed: 8, ultra: true, axis: 'time', hpd: true, fossil: true }),
    blkDiversification: () => curves({
      series: [{ f: t => Math.pow(t, 0.55) * 0.92, c: V('primary'), w: 2.2 },
               { f: t => t * 0.9, c: V('text-muted'), w: 1.2, dash: '4 3' }],
      shade: t => Math.pow(t, 0.55) * 0.92,
    }),
    blkTraits: () => tree({ n: 9, seed: 12, clades: [0, 0, 2, 2, 2, 4, 4, 6, 6], support: true }),
    blkBiogeo: () => mapAreas({ letters: true, seed: 9 }),
    blkCompare: () => tanglegram({ seed: 4 }),
    blkReport: () => {
      const W = 200, H = 124;
      let s = rect(20, 10, 112, 104, V('card-bg'), `rx="7" stroke="${V('border')}"`);
      s += rect(30, 22, 60, 5, V('primary'), 'rx="2.5"');
      for (let i = 0; i < 5; i++) s += rect(30, 34 + i * 9, 92 - (i % 3) * 18, 4, V('border-strong'), 'rx="2" opacity="0.7"');
      s += tree({ n: 5, seed: 3, W: 60, H: 42, raw: true, x0: 92, x1: 6, y0: 80, y1: 106, tipR: 1.6, lw: 1.2 });
      s += rect(140, 26, 50, 36, V('bg-soft'), `rx="5" stroke="${V('border')}"`);
      s += rect(140, 70, 50, 36, V('bg-soft'), `rx="5" stroke="${V('border')}"`);
      s += txt(165, 48, 'SVG', 'art-mut', 9, 'middle');
      s += txt(165, 92, 'ZIP', 'art-mut', 9, 'middle');
      return wrap(`0 0 ${W} ${H}`, s);
    },

    /* --- materials (what can be analysed) --- */
    matCp: () => align({ rows: 5, cols: 20, labels: false, H: 64, W: 170, varFrac: 0.2, seed: 21 }),
    matNuclear: () => align({ rows: 5, cols: 20, labels: false, H: 64, W: 170, varFrac: 0.45, gaps: true, seed: 22 }),
    matProtein: () => align({ rows: 5, cols: 20, labels: false, H: 64, W: 170, protein: true, varFrac: 0.35, seed: 23 }),
    matMorph: () => grid({ k: 5, W: 170, H: 64, seed: 24, f: (i, j) => ((i * 3 + j * 5) % 4) / 3 }),
    matCodon: () => {
      let s = '';
      for (let c = 0; c < 18; c++) for (let i = 0; i < 3; i++) {
        s += rect(6 + c * 9, 12 + i * 13, 7.6, 11, V(BASE_V[(c * 3 + i * 5) % 4]), 'rx="1.2" opacity="0.9"');
      }
      for (let c = 0; c < 6; c++) s += rect(6 + c * 27, 8, 25, 2.6, V('accent'), 'rx="1.3"');
      return wrap('0 0 170 64', s);
    },
    matDated: () => tree({ n: 6, W: 170, H: 64, ultra: true, axis: 'time', seed: 26, tipR: 1.6, lw: 1.2 }),

    /* --- theory figures --- */
    theoryAnatomy, theorySaturation, theoryLongBranch, theoryClock, theoryGenomes,

    /* --- brand mark --- */
    brand: () => `<svg viewBox="0 0 32 32" xmlns="http://www.w3.org/2000/svg" aria-hidden="true">
      <path d="M4 16 H10 M10 6 V26 M10 6 H18 M10 26 H16 M18 2 V10 M18 2 H27 M18 10 H27 M16 20 V31 M16 20 H27 M16 31 H27"
        stroke="${V('primary')}" stroke-width="2.4" fill="none" stroke-linecap="round" stroke-linejoin="round"/>
      <circle cx="27" cy="2" r="2.6" fill="${V('c1')}"/><circle cx="27" cy="10" r="2.6" fill="${V('c3')}"/>
      <circle cx="27" cy="20" r="2.6" fill="${V('accent')}"/><circle cx="27" cy="31" r="2.6" fill="${V('c5')}"/></svg>`,
  };

  /* --- method cards: one entry per method of the gallery --- */
  const METHOD_ART = {
    mFasta: () => align({ rows: 7, cols: 24, gaps: true, seed: 31 }),
    mAlign: () => align({ rows: 7, cols: 24, gaps: true, highlight: 'variable', seed: 32 }),
    mTrim: () => align({ rows: 7, cols: 24, gaps: true, trimmed: true, seed: 33 }),
    mPartition: () => align({ rows: 6, cols: 24, partitions: true, partNames: ['cp', 'ITS', 'nuc'], seed: 34 }),
    mCodonAlign: () => align({ rows: 6, cols: 24, partitions: true, partNames: ['1', '2', '3'], seed: 35 }),
    mQmatrix: () => grid({ k: 4, labels: ['A', 'C', 'G', 'T'], seed: 36, f: (i, j) => (isTransition(i, j) ? 0.85 : 0.3) }),
    mModelTest: () => bars({ values: [1, 0.96, 0.9, 0.84, 0.66, 0.52, 0.38, 0.24], best: 0 }),
    mGamma: () => curves({ series: [{ f: t => Math.exp(-3.2 * t) * 0.95 + 0.03, c: V('accent'), w: 2.2 }, { f: t => 0.42 * Math.exp(-Math.pow((t - 0.42) / 0.26, 2)), c: V('primary'), w: 1.8 }] }),
    mInvariant: () => align({ rows: 7, cols: 24, varFrac: 0.12, highlight: 'variable', seed: 37 }),
    mSaturation: () => curves({ series: [{ f: t => 0.74 * (1 - Math.exp(-2.6 * t)), c: V('primary'), w: 2.2 }, { f: t => t * 0.85, c: V('text-muted'), w: 1.2, dash: '4 3' }] }),
    mComposition: () => bars({ values: [0.62, 0.38, 0.4, 0.6, 0.58, 0.42], colours: ['base-a', 'base-c', 'base-g', 'base-t', 'base-a', 'base-c'] }),
    mProteinModel: () => grid({ k: 6, seed: 38 }),
    mMk: () => grid({ k: 3, labels: ['0', '1', '2'], seed: 39 }),
    mDistance: () => grid({ k: 6, seed: 40, f: (i, j) => Math.abs(i - j) / 6 }),
    mNJ: () => tree({ n: 9, seed: 41 }),
    mBionj: () => tree({ n: 9, seed: 42, spread: 0.8 }),
    mUPGMA: () => tree({ n: 9, seed: 43, ultra: true }),
    mME: () => unrooted({ n: 10, seed: 44 }),
    mParsimony: () => tree({ n: 9, seed: 45, support: true }),
    mSankoff: () => grid({ k: 4, labels: ['A', 'C', 'G', 'T'], seed: 46, f: (i, j) => (i === j ? 0 : Math.abs(i - j) / 3) }),
    mBremer: () => tree({ n: 8, seed: 47, support: true, spread: 0.4 }),
    mConsensus: () => tree({ n: 9, seed: 48, support: true, spread: 0.2 }),
    mML: () => tree({ n: 10, seed: 49, support: true }),
    mNNI: () => {
      let s = tree({ n: 6, seed: 50, W: 96, H: 110, raw: true, x1: 8 });
      s += `<g transform="translate(100,0)">${tree({ n: 6, seed: 57, W: 96, H: 110, raw: true, x1: 8 })}</g>`;
      s += path('M92 60 L108 60', V('accent'), 2, 'marker-end="url(#phyArrow2)"');
      s = `<defs><marker id="phyArrow2" markerWidth="6" markerHeight="6" refX="5" refY="3" orient="auto"><path d="M0 0 L6 3 L0 6 z" fill="${V('accent')}"/></marker></defs>` + s;
      return wrap('0 0 200 124', s);
    },
    mBootstrap: () => {
      let s = '';
      for (let i = 0; i < 3; i++) s += `<g transform="translate(${i * 8},${i * 6}) scale(0.82)" opacity="${0.35 + i * 0.3}">${tree({ n: 7, seed: 51 + i, W: 200, H: 124, raw: true })}</g>`;
      return wrap('0 0 200 124', s);
    },
    mUFBoot: () => bars({ values: [0.98, 0.95, 0.92, 0.84, 0.72, 0.58, 0.41], line: 0.95 }),
    mSHaLRT: () => tree({ n: 9, seed: 52, support: true }),
    mTopoTest: () => bars({ values: [1, 0.72, 0.55, 0.3], best: 0, colours: ['primary', 'c5', 'c5', 'c5'] }),
    mConstraint: () => tree({ n: 9, seed: 53, clades: [0, 0, 0, 0, 2, 2, 4, 4, 4] }),
    mMCMC: () => curves({ series: [{ f: t => 0.3 + 0.55 * (1 - Math.exp(-6 * t)) + 0.05 * Math.sin(t * 40), c: V('primary'), w: 1.2 }, { f: t => 0.3 + 0.55 * (1 - Math.exp(-5 * t)) + 0.05 * Math.sin(t * 37 + 2), c: V('accent'), w: 1.2 }] }),
    mPosterior: () => curves({ series: [{ f: t => Math.exp(-Math.pow((t - 0.5) / 0.16, 2)) * 0.92, c: V('accent'), w: 2.2 }], shade: t => Math.exp(-Math.pow((t - 0.5) / 0.16, 2)) * 0.92 }),
    mConvergence: () => curves({ series: [{ f: t => 0.75 * Math.exp(-4 * t) + 0.03, c: V('primary'), w: 2.2 }], points: [[0.75, 0.06], [0.9, 0.04]] }),
    mMC3: () => bars({ values: [1, 0.78, 0.6, 0.45], colours: ['danger', 'accent', 'c5', 'primary'] }),
    mClockTest: () => tree({ n: 8, seed: 54, ultra: true }),
    mRooting: () => tree({ n: 9, seed: 55, outgroup: true }),
    mLSD: () => tree({ n: 9, seed: 56, ultra: true, axis: 'time' }),
    mFossil: () => tree({ n: 8, seed: 58, ultra: true, axis: 'time', fossil: true, hpd: true }),
    mRelaxed: () => tree({ n: 9, seed: 59, ultra: true, axis: 'time', clades: [0, 0, 2, 2, 4, 4, 6, 6, 8] }),
    mTipDating: () => curves({ series: [{ f: t => 0.15 + 0.7 * t, c: V('primary'), w: 2 }], points: [[0.1, 0.22], [0.3, 0.35], [0.5, 0.52], [0.7, 0.62], [0.9, 0.78]] }),
    mLTT: () => curves({ series: [{ f: t => Math.pow(t, 0.6) * 0.92, c: V('primary'), w: 2.2 }], shade: t => Math.pow(t, 0.6) * 0.92 }),
    mBirthDeath: () => curves({ series: [{ f: t => 0.2 + 0.7 * t * t, c: V('leaf'), w: 2.2 }, { f: t => 0.18 + 0.25 * t, c: V('rose'), w: 2.2 }] }),
    mRateShift: () => tree({ n: 10, seed: 60, clades: [0, 0, 0, 0, 0, 2, 2, 4, 4, 4], support: true }),
    mAncestral: () => tree({ n: 9, seed: 61, clades: [0, 0, 2, 0, 2, 2, 4, 4, 2], support: true }),
    mSimmap: () => tree({ n: 9, seed: 62, clades: [2, 0, 0, 2, 2, 4, 2, 4, 0] }),
    mBMOU: () => curves({ series: [{ f: t => 0.5 + 0.28 * Math.sin(t * 13) * Math.sqrt(t), c: V('primary'), w: 1.6 }, { f: t => 0.5 + 0.2 * Math.sin(t * 17 + 1) * Math.sqrt(t), c: V('c5'), w: 1.6 }, { f: t => 0.5 + 0.24 * Math.sin(t * 11 + 2) * Math.sqrt(t), c: V('leaf'), w: 1.6 }] }),
    mSignal: () => curves({ series: [{ f: t => 0.12 + 0.78 * t, c: V('primary'), w: 2 }], points: [[0.15, 0.2], [0.35, 0.42], [0.55, 0.5], [0.75, 0.66], [0.92, 0.85]] }),
    mPGLS: () => curves({ series: [{ f: t => 0.2 + 0.62 * t, c: V('accent'), w: 2 }], points: [[0.12, 0.3], [0.3, 0.36], [0.48, 0.5], [0.66, 0.55], [0.85, 0.74]] }),
    mDEC: () => mapAreas({ letters: true, seed: 63 }),
    mDECJ: () => mapAreas({ letters: true, seed: 64 }),
    mStratified: () => mapAreas({ letters: true, seed: 65, arrows: false, pie: true }),
    mBSM: () => tree({ n: 9, seed: 66, clades: [0, 0, 2, 2, 4, 4, 6, 0, 2] }),
    mRF: () => tanglegram({ seed: 67 }),
    mConcordance: () => bars({ values: [0.92, 0.8, 0.62, 0.5, 0.38, 0.2], colours: ['leaf', 'leaf', 'gold', 'gold', 'rose', 'rose'] }),
    mSpeciesTree: () => {
      let s = tree({ n: 7, seed: 68, W: 200, H: 124, raw: true, lw: 4, stroke: V('border-strong') });
      s += tree({ n: 7, seed: 69, W: 200, H: 124, raw: true, lw: 1.2, tipR: 1.4 });
      return wrap('0 0 200 124', s);
    },
    mNetwork: () => network({ n: 9, boxes: true, seed: 70 }),
    mABBA: () => bars({ values: [0.55, 0.85, 0.3, 0.6], colours: ['c5', 'accent', 'c5', 'c5'] }),
    mTreeStudio: () => treeCircular({ n: 16, seed: 71, ring: true }),
    mExport: () => SPECS.blkReport(),
  };

  window.Art = Object.assign({}, SPECS, METHOD_ART, {
    tree, treeCircular, unrooted, align, curves, bars, grid, mapAreas, network, tanglegram, hero,
  });
})();
