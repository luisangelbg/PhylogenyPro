/* PhylogenyPro — molecular phylogenetics in the browser.
   Copyright (C) 2026  Luis Ángel Barrera-Guzmán

   This program is free software: you can redistribute it and/or modify it under
   the terms of the GNU General Public License as published by the Free Software
   Foundation, either version 3 of the License, or (at your option) any later
   version. It is distributed in the hope that it will be useful, but WITHOUT ANY
   WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
   A PARTICULAR PURPOSE. See the GNU General Public License, in the file LICENSE
   at the root of this program, or <https://www.gnu.org/licenses/>. */

/* PhylogenyPro — drawing trees.

   One renderer for the whole app: Block 4 shows parsimony and distance trees
   with it, Block 5 the likelihood tree, Block 7 the dated tree, and Block 12
   will wrap it in the tree studio. It returns a real SVG string, so the same
   picture exports at any resolution and follows the light and dark theme.

   Layouts
     rect       phylogram or cladogram, the usual shape
     circular   the same, wrapped around, for many taxa
     unrooted   equal-angle, for a tree whose root is arbitrary

   Everything optional is really optional: support values, branch lengths, a
   scale bar, coloured clades, highlighted tips. */

const TreeView = {};

(function () {

  const V = n => `var(--${n})`;
  const f1 = v => (+v).toFixed(1);
  const esc2 = s => String(s == null ? '' : s).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
  const line = (x1, y1, x2, y2, st, w, ex) => `<line x1="${f1(x1)}" y1="${f1(y1)}" x2="${f1(x2)}" y2="${f1(y2)}" stroke="${st}" stroke-width="${w}" stroke-linecap="round" ${ex || ''}/>`;
  const circ = (cx, cy, r, fill, ex) => `<circle cx="${f1(cx)}" cy="${f1(cy)}" r="${f1(r)}" fill="${fill}" ${ex || ''}/>`;
  const txt = (x, y, s, size, anchor, cls, ex) =>
    `<text x="${f1(x)}" y="${f1(y)}" font-size="${size}" text-anchor="${anchor || 'start'}" class="${cls || 'art-txt'}" ${ex || ''}>${esc2(s)}</text>`;
  /* marks for the figure studio: the area of the drawing (data-plot) and the
     legend, in one group (data-role="legend") with its entries numbered (data-li) */
  const plotArea = (x, y, w, h) => [x, y, w, h].map(v => +(+v).toFixed(2)).join(' ');

  function defaults(opts) {
    return Object.assign({
      layout: 'rect', width: 760, cladogram: false,
      showSupport: true, supportMin: 0, supportAs: 'node', supportLabel: null,
      showScale: true, tipFont: 12, lineWidth: 1.8, italicTips: true,
      padding: 14, labels: null, colours: null, highlight: null, title: null,
      rowHeight: 18, maxLabel: 260,
    }, opts || {});
  }

  /* how wide the longest label will be, roughly */
  function labelWidth(labels, font) {
    if (!labels || !labels.length) return 0;
    const longest = labels.reduce((a, b) => (String(a).length > String(b).length ? a : b), '');
    return Math.min(260, String(longest).length * font * 0.56 + 10);
  }

  /* ================================================================
     keeping labels apart
     ================================================================ */
  /* Nudge a set of one-dimensional positions until no two are closer than gap,
     staying as near as possible to where each one wanted to be and inside
     [lo, hi]. Two passes of the classic isotonic push: forward, then backward,
     which is enough because the wanted positions are already sorted.

     Items: [{pos, gap?}] in any order; returns the new positions in the same
     order. When there is simply no room — more labels than the space allows —
     it spreads them evenly and reports `crowded`, so the caller can hide some
     instead of printing a black smear. */
  function spreadLabels(wanted, gap, lo, hi) {
    const n = wanted.length;
    if (!n) return { pos: [], crowded: false, shown: [] };
    const idx = wanted.map((p, i) => i).sort((a, b) => wanted[a] - wanted[b]);
    const p = idx.map(i => Math.min(hi, Math.max(lo, wanted[i])));
    const crowded = (n - 1) * gap > (hi - lo) + 1e-9;

    if (crowded) {
      const step = n > 1 ? (hi - lo) / (n - 1) : 0;
      const out = new Array(n);
      idx.forEach((orig, k) => { out[orig] = lo + k * step; });
      return { pos: out, crowded: true, shown: null };
    }
    /* forward: push each one down until it clears the previous */
    for (let k = 1; k < n; k++) if (p[k] - p[k - 1] < gap) p[k] = p[k - 1] + gap;
    /* backward: if the last one ran past the end, pull the tail back up */
    if (p[n - 1] > hi) {
      p[n - 1] = hi;
      for (let k = n - 2; k >= 0; k--) if (p[k + 1] - p[k] < gap) p[k] = p[k + 1] - gap;
    }
    const out = new Array(n);
    idx.forEach((orig, k) => { out[orig] = p[k]; });
    return { pos: out, crowded: false, shown: null };
  }

  /* When even spreading will not do — a scatter, where a label may move in two
     directions — drop the ones that still collide, keeping the most important.
     `boxes` are {x, y, w, h, rank}; a lower rank is kept first. */
  function dropColliding(boxes) {
    const order = boxes.map((b, i) => i).sort((a, b) => (boxes[a].rank || 0) - (boxes[b].rank || 0));
    const kept = [];
    const keep = new Array(boxes.length).fill(false);
    for (const i of order) {
      const b = boxes[i];
      const hits = kept.some(j => {
        const a = boxes[j];
        return !(b.x + b.w < a.x || a.x + a.w < b.x || b.y + b.h < a.y || a.y + a.h < b.y);
      });
      if (!hits) { keep[i] = true; kept.push(i); }
    }
    return keep;
  }

  /* ================================================================
     pictures of the OTUs and of the clades
     ================================================================ */
  /* A clip path per shape, defined once and referred to by every picture. */
  function imageDefs(shape, size, id) {
    const h = size / 2, r = shape === 'circle' ? h : shape === 'round' ? size * 0.18 : 0;
    if (shape === 'none') return '';
    const body = shape === 'circle'
      ? `<circle cx="${f1(h)}" cy="${f1(h)}" r="${f1(h)}"/>`
      : `<rect x="0" y="0" width="${f1(size)}" height="${f1(size)}" rx="${f1(r)}" ry="${f1(r)}"/>`;
    return `<clipPath id="${id}" clipPathUnits="userSpaceOnUse">${body}</clipPath>`;
  }

  /* one picture, centred on (x, y), with its frame and the line that ties it to
     the branch it belongs to */
  function imageNode(item, x, y, size, shape, frame) {
    const h = size / 2;
    const id = item.clip;
    let s = '';
    if (item.linkTo != null) {
      /* The leader stops at the edge of the picture on the side the branch is
         actually on. Ending it at a fixed point to the left works in the
         rectangular layout and is nonsense in a radial one, where a picture at
         the top of the circle would be pointed at from its side. */
      const dx = item.linkTo[0] - x, dy = item.linkTo[1] - y;
      const d = Math.hypot(dx, dy) || 1;
      const ex = x + (dx / d) * (h + 2), ey = y + (dy / d) * (h + 2);
      s += `<line x1="${f1(item.linkTo[0])}" y1="${f1(item.linkTo[1])}" x2="${f1(ex)}" y2="${f1(ey)}" ` +
           `stroke="${item.colour || V('border')}" stroke-width="0.9" stroke-dasharray="2 3" opacity="0.8"/>`;
    }
    const clip = shape === 'none' ? '' : ` clip-path="url(#${id})"`;
    s += `<g transform="translate(${f1(x - h)},${f1(y - h)})">` +
      `<image x="0" y="0" width="${f1(size)}" height="${f1(size)}" href="${item.url}" ` +
      `preserveAspectRatio="${shape === 'none' ? 'xMidYMid meet' : 'xMidYMid slice'}"${clip}/>`;
    if (frame !== false) {
      const r = shape === 'circle' ? h : shape === 'round' ? size * 0.18 : 0;
      s += shape === 'circle'
        ? `<circle cx="${f1(h)}" cy="${f1(h)}" r="${f1(h - 0.6)}" fill="none" stroke="${item.colour || V('border')}" stroke-width="1.6"/>`
        : `<rect x="0.8" y="0.8" width="${f1(size - 1.6)}" height="${f1(size - 1.6)}" rx="${f1(r)}" ry="${f1(r)}" fill="none" stroke="${item.colour || V('border')}" stroke-width="1.6"/>`;
    }
    s += '</g>';
    return s;
  }

  /* Collect what has to be drawn: a picture per tip, per named group, or per
     clade the user marked. Nothing is invented — a name with no confirmed
     picture simply has none. */
  function collectImages(cfg, ctx) {
    const src = (typeof window !== 'undefined' && window.OTUImg) ? window.OTUImg : null;
    if (!src || !cfg.images || cfg.images === 'none') return [];
    const out = [];
    if (cfg.images === 'tips') {
      ctx.tips.forEach(t => {
        const nm = ctx.labelOf(t);
        const u = src.url(nm);
        if (u) out.push({ name: nm, url: u, at: ctx.yOf(t), range: [ctx.yOf(t), ctx.yOf(t)], colour: ctx.colourOf(t) });
      });
    } else if (cfg.images === 'groups' && cfg.groups) {
      const g = cfg.groups;
      for (let k = 0; k < g.k; k++) {
        const nm = window.Groups ? Groups.nameOf(g, k) : String(k);
        const u = src.url(nm);
        if (!u) continue;
        const ys = ctx.tips.filter(t => g.of.get(t.tip) === k).map(t => ctx.yOf(t));
        if (!ys.length) continue;
        out.push({ name: nm, url: u, at: (Math.min.apply(null, ys) + Math.max.apply(null, ys)) / 2,
                   range: [Math.min.apply(null, ys), Math.max.apply(null, ys)], colour: ctx.groupColour(k) });
      }
    } else if (cfg.images === 'clades' && cfg.cladeImages) {
      cfg.cladeImages.forEach(cl => {
        const u = src.url(cl.name);
        if (!u) return;
        const ys = ctx.tips.filter(t => cl.tips.indexOf(t.tip) >= 0).map(t => ctx.yOf(t));
        if (!ys.length) return;
        out.push({ name: cl.name, url: u, at: (Math.min.apply(null, ys) + Math.max.apply(null, ys)) / 2,
                   range: [Math.min.apply(null, ys), Math.max.apply(null, ys)], colour: cl.colour || V('text') });
      });
    }
    out.sort((a, b) => a.at - b.at);
    out.forEach((it, i) => { it.clip = 'pp-clip-' + i; });
    return out;
  }

  /* ================================================================
     colour by group
     ================================================================ */
  function groupPalette(cfg) {
    if (cfg.groupColours && cfg.groupColours.length) return cfg.groupColours;
    if (typeof window !== 'undefined' && window.Groups && cfg.groups) return Groups.palette(cfg.groups.k);
    return [];
  }

  /* ================================================================
     rectangular
     ================================================================ */
  function rect(tree, opts) {
    opts = defaults(opts);
    const labels = opts.labels;
    const tips = Tree.tips(tree);
    const n = tips.length;
    const pal = groupPalette(opts);
    const gOf = t => (opts.groups ? opts.groups.of.get(t) : null);
    const groupColour = g => (g == null || !pal.length ? V('text') : pal[g % pal.length]);

    /* pictures need a column of their own on the right, and a row tall enough
       to hold one; the tree simply gets shorter, never overlapped */
    const wantImg = opts.images && opts.images !== 'none';
    const imgSize = wantImg ? Math.max(18, +opts.imgSize || 54) : 0;
    const imgGap = wantImg ? 12 : 0;
    const rowH = Math.max(opts.rowHeight, wantImg && opts.images === 'tips' ? imgSize + 4 : 0);

    const legendH = (opts.legend && opts.groups && opts.groups.k > 1) ? 22 : 0;
    const H = Math.max(120, opts.padding * 2 + n * rowH + (opts.showScale ? 34 : 0) + (opts.title ? 22 : 0) + legendH);
    const labW = labelWidth(tips.map(t => (labels ? labels[t.tip] : t.label) || ''), opts.tipFont);
    const W = opts.width;
    const stripW = opts.groupStrip && opts.groups && opts.groups.k > 1 ? 9 : 0;
    const x0 = opts.padding + 4;
    const x1 = W - opts.padding - labW - 8 - stripW - (wantImg ? imgSize + imgGap : 0);
    const y0 = opts.padding + (opts.title ? 22 : 0);
    const y1 = H - opts.padding - (opts.showScale ? 34 : 0) - legendH;

    /* depth of every node: additive lengths, or topological for a cladogram */
    const depth = new Map();
    let maxDepth = 0;
    (function walk(nd, d) {
      depth.set(nd, d);
      if (d > maxDepth) maxDepth = d;
      nd.children.forEach(c => walk(c.node, d + (opts.cladogram ? 1 : (c.len || 0))));
    })(tree, 0);
    /* a cladogram puts every tip at the right edge */
    if (opts.cladogram) {
      const tipDepth = Math.max(...tips.map(t => depth.get(t)));
      tips.forEach(t => depth.set(t, tipDepth));
      maxDepth = tipDepth;
    }
    const X = d => x0 + (maxDepth > 0 ? (x1 - x0) * d / maxDepth : 0);

    /* y of every node */
    const Y = new Map();
    tips.forEach((t, i) => Y.set(t, n === 1 ? (y0 + y1) / 2 : y0 + (y1 - y0) * i / (n - 1)));
    (function setY(nd) {
      if (!nd.children.length) return Y.get(nd);
      const ys = nd.children.map(c => setY(c.node));
      Y.set(nd, (Math.min(...ys) + Math.max(...ys)) / 2);
      return Y.get(nd);
    })(tree);

    /* Colour, in order of precedence: a clade the user painted by hand, then the
       group a branch belongs to when all of its tips agree, then plain text
       colour. A branch whose descendants are of two groups stays neutral — the
       one thing it must not do is pick one of them and mislead. */
    const colourOf = nd => {
      if (opts.colours) {
        const g = cladeGroup(nd, opts.colours);
        if (g != null) return V('c' + (1 + (g % 10)));
      }
      if (opts.groups && opts.colourBranches !== false && pal.length) {
        const g = window.Groups ? Groups.groupOfNode(nd, opts.groups) : null;
        if (g != null) return groupColour(g);
      }
      return V('text');
    };

    let s = '';
    let defs = '';
    if (opts.title) s += txt(opts.padding, 16, opts.title, 13, 'start', 'art-txt', 'font-weight="600"');

    /* A soft band behind each group, drawn first so everything sits on top.
       Only for a group whose tips are consecutive in the drawing: a band that
       spanned other groups' tips would claim they belong to it. A paraphyletic
       "rest of the tree" therefore gets no band, which is the honest outcome. */
    if (opts.groupBands && opts.groups && opts.groups.k > 1) {
      for (let g = 0; g < opts.groups.k; g++) {
        const rows = tips.map((t, i) => (gOf(t.tip) === g ? i : -1)).filter(i => i >= 0);
        if (rows.length < 2) continue;
        if (rows[rows.length - 1] - rows[0] !== rows.length - 1) continue;   // not consecutive
        const ys = rows.map(i => Y.get(tips[i]));
        const a = Math.min.apply(null, ys) - rowH * 0.45, b = Math.max.apply(null, ys) + rowH * 0.45;
        s += `<rect x="${f1(x0 - 2)}" y="${f1(a)}" width="${f1(W - opts.padding - x0 + 2)}" height="${f1(b - a)}" ` +
             `fill="${groupColour(g)}" opacity="0.07"/>`;
      }
    }

    (function draw(nd, parentX) {
      const x = X(depth.get(nd)), y = Y.get(nd);
      const col = colourOf(nd);
      s += line(parentX, y, x, y, col, opts.lineWidth);
      if (nd.children.length) {
        const ys = nd.children.map(c => Y.get(c.node));
        s += line(x, Math.min(...ys), x, Math.max(...ys), col, opts.lineWidth);
        nd.children.forEach(c => draw(c.node, x));
      }
    })(tree, X(0));

    /* tips */
    const labelX = x1 + 6;
    tips.forEach(t => {
      const x = X(depth.get(t)), y = Y.get(t);
      const name = (labels ? labels[t.tip] : t.label) || '';
      const hl = opts.highlight && opts.highlight.indexOf(t.tip) >= 0;
      const g = gOf(t.tip);
      s += circ(x, y, 2.4, colourOf(t));
      /* aligned names read as a column and let the pictures line up with them;
         the dotted leader keeps the eye on the right branch */
      const tx = opts.alignTips ? labelX : x + 6;
      if (opts.alignTips && labelX - x > 10) {
        s += line(x + 4, y, labelX - 3, y, V('border'), 0.8, 'stroke-dasharray="1.5 3" opacity="0.7"');
      }
      const fill = hl ? V('accent') : (opts.colourLabels !== false && g != null && pal.length ? groupColour(g) : null);
      s += txt(tx, y + opts.tipFont * 0.36, name, opts.tipFont, 'start', 'art-txt',
        (opts.italicTips ? 'font-style="italic" ' : '') +
        (hl ? 'font-weight="700" ' : '') + (fill ? `fill="${fill}"` : ''));
    });

    /* a strip of group colour, right of the names: the grouping stated plainly
       instead of left to be guessed from the branch colours */
    if (stripW) {
      const sx = W - opts.padding - (wantImg ? imgSize + imgGap : 0) - stripW;
      for (let g = 0; g < opts.groups.k; g++) {
        const ys = tips.filter(t => gOf(t.tip) === g).map(t => Y.get(t));
        if (!ys.length) continue;
        const a = Math.min.apply(null, ys) - rowH * 0.42, b = Math.max.apply(null, ys) + rowH * 0.42;
        s += `<rect x="${f1(sx)}" y="${f1(a)}" width="${f1(stripW - 3)}" height="${f1(b - a)}" rx="2" ` +
             `fill="${groupColour(g)}" opacity="0.85"/>`;
        if (b - a > 26 && window.Groups) {
          s += txt(sx + (stripW - 3) / 2, (a + b) / 2, Groups.nameOf(opts.groups, g), 9, 'middle', 'art-txt',
            `fill="#fff" font-weight="700" transform="rotate(-90 ${f1(sx + (stripW - 3) / 2)} ${f1((a + b) / 2)})"`);
        }
      }
    }

    /* support values */
    let hiddenSupport = 0;
    if (opts.showSupport) {
      const marks = [];
      Tree.nodes(tree).forEach(nd => {
        if (!nd.children.length || nd.support == null) return;
        if (nd.support < opts.supportMin) return;
        marks.push({ nd, x: X(depth.get(nd)), y: Y.get(nd) });
      });
      if (opts.supportAs === 'dot') {
        marks.forEach(m => {
          const f = m.nd.support > 1 ? m.nd.support / 100 : m.nd.support;
          s += circ(m.x, m.y, 3.4, f >= 0.95 ? V('leaf') : f >= 0.7 ? V('gold') : V('rose'));
        });
      } else {
        const fs = Math.max(8, opts.tipFont - 3);
        /* A bootstrap arrives as a percentage and a posterior probability as a
           fraction: print the first whole and the second with two decimals,
           unless the caller asks for a precision. */
        marks.forEach(m => {
          const v = opts.supportAs === 'percent' ? Math.round(m.nd.support) : m.nd.support;
          m.lab = opts.supportDecimals != null ? (+v).toFixed(opts.supportDecimals)
            : (v <= 1 ? (+v).toFixed(2) : String(Math.round(v)));
        });
        /* Two numbers printed on top of each other are worse than one number
           and an honest count of what was left out. The most inclusive clades
           keep their value; the caller is told how many were dropped. */
        const boxes = marks.map(m => ({
          x: m.x - 3 - m.lab.length * fs * 0.56, y: m.y - 4 - fs,
          w: m.lab.length * fs * 0.56 + 3, h: fs + 2,
          rank: -Tree.nTips(m.nd),
        }));
        const keep = opts.supportDeclutter === false ? boxes.map(() => true) : dropColliding(boxes);
        marks.forEach((m, i) => {
          if (!keep[i]) { hiddenSupport++; return; }
          s += txt(m.x - 3, m.y - 4, m.lab, fs, 'end', 'art-mut');
        });
      }
    }

    /* pictures of the OTUs, of the groups or of the clades the user marked */
    const imgs = collectImages(opts, {
      tips, yOf: t => Y.get(t), labelOf: t => (labels ? labels[t.tip] : t.label) || '',
      colourOf: t => colourOf(t), groupColour,
    });
    if (imgs.length) {
      const shape = opts.imgShape || 'round';
      const cx = W - opts.padding - imgSize / 2;
      /* one column, each picture as near as possible to what it illustrates and
         never on top of the next one */
      const lo = y0 + imgSize / 2, hi = y1 - imgSize / 2;
      const placed = spreadLabels(imgs.map(it => it.at), imgSize + 4, lo, hi);
      imgs.forEach((it, i) => {
        defs += imageDefs(shape, imgSize, it.clip);
        it.linkTo = [W - opts.padding - imgSize - imgGap + 2, it.at];
        s += imageNode(it, cx, placed.pos[i], imgSize, shape, opts.imgFrame);
      });
    }

    /* the legend of the groups */
    if (legendH) {
      const y = H - opts.padding - 4;
      let lx = opts.padding + 2;
      s += '<g data-role="legend">';
      for (let g = 0; g < opts.groups.k; g++) {
        const nm = (window.Groups ? Groups.nameOf(opts.groups, g) : String(g + 1));
        s += `<rect x="${f1(lx)}" y="${f1(y - 9)}" width="10" height="10" rx="2" fill="${groupColour(g)}" data-li="${g}"/>`;
        s += txt(lx + 14, y, nm, 10, 'start', 'art-mut', `data-li="${g}"`);
        lx += 14 + nm.length * 5.6 + 14;
      }
      s += '</g>';
    }

    /* scale bar */
    if (opts.showScale && !opts.cladogram && maxDepth > 0) {
      const nice = niceStep(maxDepth);
      const px = (x1 - x0) * nice / maxDepth;
      const y = H - opts.padding - 12 - legendH;
      s += line(x0, y, x0 + px, y, V('text'), 1.6);
      s += line(x0, y - 4, x0, y + 4, V('text'), 1.6);
      s += line(x0 + px, y - 4, x0 + px, y + 4, V('text'), 1.6);
      s += txt(x0 + px + 8, y + 4, nice + (opts.scaleLabel || ''), 10, 'start', 'art-mut');
    }
    TreeView.lastInfo = { hiddenSupport, images: imgs.length, width: W, height: H, k: opts.groups ? opts.groups.k : 0 };
    return `<svg viewBox="0 0 ${W} ${H}" data-plot="${plotArea(x0, y0, x1 - x0, y1 - y0)}" xmlns="http://www.w3.org/2000/svg">` +
      (defs ? `<defs>${defs}</defs>` : '') + s + '</svg>';
  }
  function niceStep(max) {
    const raw = max / 5;
    const mag = Math.pow(10, Math.floor(Math.log10(raw)));
    const n = raw / mag;
    const step = n < 1.5 ? 1 : n < 3 ? 2 : n < 7 ? 5 : 10;
    return +(step * mag).toPrecision(2);
  }
  /* the colour group of a node, when every tip below it belongs to the same one */
  function cladeGroup(nd, groups) {
    const seen = new Set();
    (function walk(x) { if (!x.children.length) seen.add(groups[x.tip]); else x.children.forEach(c => walk(c.node)); })(nd);
    return seen.size === 1 ? [...seen][0] : null;
  }

  /* ================================================================
     circular
     ================================================================ */
  function circular(tree, opts) {
    opts = defaults(opts);
    const labels = opts.labels;
    const tips = Tree.tips(tree);
    const n = tips.length;
    const labW = labelWidth(tips.map(t => (labels ? labels[t.tip] : t.label) || ''), opts.tipFont);

    const pal = groupPalette(opts);
    const gOf = t => (opts.groups ? opts.groups.of.get(t) : null);
    const groupColour = g => (g == null || !pal.length ? V('text') : pal[g % pal.length]);
    const wantImg = opts.images && opts.images !== 'none';
    const imgSize = wantImg ? Math.max(18, +opts.imgSize || 54) : 0;
    const ringW = opts.groupStrip && opts.groups && opts.groups.k > 1 ? 8 : 0;

    /* The ring and the pictures go OUTSIDE, on a wider canvas. Taking their
       room out of the radius instead would shrink the tree towards nothing as
       soon as the names were long: with a 430-wide figure and names like
       "Bursera copallifera" the tree collapsed to seven pixels. The width the
       caller asks for is the width of the tree and its names; whatever is drawn
       around them is added to the figure. */
    const R1 = opts.width / 2 - labW - 14;
    const extra = ringW + (wantImg ? imgSize + 12 : 0);
    const W = opts.width + 2 * extra, H = W;
    const cx = W / 2, cy = H / 2;
    const open = opts.openAngle == null ? 0.25 : opts.openAngle;

    const depth = new Map();
    let maxDepth = 0;
    (function walk(nd, d) {
      depth.set(nd, d);
      if (d > maxDepth) maxDepth = d;
      nd.children.forEach(c => walk(c.node, d + (opts.cladogram ? 1 : (c.len || 0))));
    })(tree, 0);
    if (opts.cladogram) {
      const td = Math.max(...tips.map(t => depth.get(t)));
      tips.forEach(t => depth.set(t, td));
      maxDepth = td;
    }
    const rad = d => 10 + (R1 - 10) * (maxDepth > 0 ? d / maxDepth : 0);
    const ang = new Map();
    tips.forEach((t, i) => ang.set(t, -Math.PI / 2 + open / 2 + (2 * Math.PI - open) * i / (n - 1 || 1)));
    (function setA(nd) {
      if (!nd.children.length) return ang.get(nd);
      const as = nd.children.map(c => setA(c.node));
      ang.set(nd, (Math.min(...as) + Math.max(...as)) / 2);
      return ang.get(nd);
    })(tree);
    const P = (r, a) => [cx + r * Math.cos(a), cy + r * Math.sin(a)];

    let s = '';
    let defs = '';
    if (opts.title) s += txt(12, 18, opts.title, 13, 'start', 'art-txt', 'font-weight="600"');
    const colourOf = nd => {
      if (opts.colours) {
        const g = cladeGroup(nd, opts.colours);
        if (g != null) return V('c' + (1 + (g % 10)));
      }
      if (opts.groups && opts.colourBranches !== false && pal.length) {
        const g = window.Groups ? Groups.groupOfNode(nd, opts.groups) : null;
        if (g != null) return groupColour(g);
      }
      return V('text');
    };

    /* an arc of group colour outside the names — the circular answer to the
       strip of the rectangular layout; only for a group whose tips are
       consecutive, so an arc never spans taxa that are not in it */
    const arcOf = (r, a0, a1, col, w) => {
      const [sx, sy] = P(r, a0), [ex, ey] = P(r, a1);
      return `<path d="M${f1(sx)} ${f1(sy)} A ${f1(r)} ${f1(r)} 0 ${a1 - a0 > Math.PI ? 1 : 0} 1 ${f1(ex)} ${f1(ey)}" ` +
        `stroke="${col}" stroke-width="${w}" fill="none" stroke-linecap="butt" opacity="0.85"/>`;
    };
    (function draw(nd, parentR) {
      const r = rad(depth.get(nd)), a = ang.get(nd), col = colourOf(nd);
      const [x1, y1] = P(parentR, a), [x2, y2] = P(r, a);
      s += line(x1, y1, x2, y2, col, opts.lineWidth);
      if (nd.children.length) {
        const as = nd.children.map(c => ang.get(c.node));
        const a0 = Math.min(...as), a1 = Math.max(...as);
        const [sx, sy] = P(r, a0), [ex, ey] = P(r, a1);
        s += `<path d="M${f1(sx)} ${f1(sy)} A ${f1(r)} ${f1(r)} 0 ${a1 - a0 > Math.PI ? 1 : 0} 1 ${f1(ex)} ${f1(ey)}" stroke="${col}" stroke-width="${opts.lineWidth}" fill="none"/>`;
        nd.children.forEach(c => draw(c.node, r));
      } else {
        s += circ(x2, y2, 2.2, col);
        const name = (labels ? labels[nd.tip] : nd.label) || '';
        const deg = a * 180 / Math.PI;
        const flip = Math.cos(a) < 0;
        /* aligned names put every picture at the same radius, which is what
           makes a ring of photographs read as a ring */
        const lr = opts.alignTips ? R1 + 5 : r + 5;
        if (opts.alignTips && lr - r > 8) {
          const [dx0, dy0] = P(r + 3, a), [dx1, dy1] = P(lr - 2, a);
          s += line(dx0, dy0, dx1, dy1, V('border'), 0.8, 'stroke-dasharray="1.5 3" opacity="0.7"');
        }
        const [lx, ly] = P(lr, a);
        const g = gOf(nd.tip);
        const fill = opts.colourLabels !== false && g != null && pal.length ? ` fill="${groupColour(g)}"` : '';
        s += `<text x="${f1(lx)}" y="${f1(ly)}" font-size="${opts.tipFont}" class="art-txt"${fill} ${opts.italicTips ? 'font-style="italic"' : ''}
          text-anchor="${flip ? 'end' : 'start'}"
          transform="rotate(${flip ? deg + 180 : deg} ${f1(lx)} ${f1(ly)})" dominant-baseline="middle">${esc2(name)}</text>`;
      }
    })(tree, 4);
    if (opts.showSupport) {
      Tree.nodes(tree).forEach(nd => {
        if (!nd.children.length || nd.support == null || nd.support < opts.supportMin) return;
        const [x, y] = P(rad(depth.get(nd)), ang.get(nd));
        const f = nd.support > 1 ? nd.support / 100 : nd.support;
        s += circ(x, y, 3, f >= 0.95 ? V('leaf') : f >= 0.7 ? V('gold') : V('rose'));
      });
    }

    /* the ring of groups, outside the names */
    const ringR = R1 + labW + 6;
    if (ringW) {
      for (let g = 0; g < opts.groups.k; g++) {
        const rows = tips.map((t, i) => (gOf(t.tip) === g ? i : -1)).filter(i => i >= 0);
        if (!rows.length) continue;
        if (rows[rows.length - 1] - rows[0] !== rows.length - 1) continue;   // not consecutive
        const as = rows.map(i => ang.get(tips[i]));
        const step = n > 1 ? (2 * Math.PI - open) / (n - 1) : 0;
        s += arcOf(ringR, Math.min.apply(null, as) - step * 0.42,
          Math.max.apply(null, as) + step * 0.42, groupColour(g), ringW - 2);
      }
    }

    /* the pictures, spread around the circle so none sits on the next */
    const imgs = collectImages(opts, {
      tips, yOf: t => ang.get(t), labelOf: t => (labels ? labels[t.tip] : t.label) || '',
      colourOf: t => colourOf(t), groupColour,
    });
    if (imgs.length) {
      const shape = opts.imgShape || 'round';
      const Rimg = ringR + ringW + imgSize / 2 + 2;
      const a0 = -Math.PI / 2 + open / 2, a1 = -Math.PI / 2 + open / 2 + (2 * Math.PI - open);
      /* the gap is an angle: how much of the circle one picture takes up at
         this radius. Spread on the angle, then place. */
      const gapA = (imgSize + 5) / Math.max(1, Rimg);
      const placed = spreadLabels(imgs.map(it => it.at), gapA, a0, a1);
      imgs.forEach((it, i) => {
        defs += imageDefs(shape, imgSize, it.clip);
        const [px, py] = P(Rimg, placed.pos[i]);
        const [ax, ay] = P(ringR + ringW, it.at);
        it.linkTo = [ax, ay];
        s += imageNode(it, px, py, imgSize, shape, opts.imgFrame);
      });
    }

    if (opts.legend && opts.groups && opts.groups.k > 1) {
      let lx = 12, ly = H - 12;
      s += '<g data-role="legend">';
      for (let g = 0; g < opts.groups.k; g++) {
        const nm = (window.Groups ? Groups.nameOf(opts.groups, g) : String(g + 1));
        s += `<rect x="${f1(lx)}" y="${f1(ly - 9)}" width="10" height="10" rx="2" fill="${groupColour(g)}" data-li="${g}"/>`;
        s += txt(lx + 14, ly, nm, 10, 'start', 'art-mut', `data-li="${g}"`);
        lx += 14 + nm.length * 5.6 + 14;
      }
      s += '</g>';
    }

    TreeView.lastInfo = { hiddenSupport: 0, images: imgs.length, width: W, height: H, k: opts.groups ? opts.groups.k : 0 };
    return `<svg viewBox="0 0 ${W} ${H}" data-plot="${plotArea(0, 0, W, H)}" xmlns="http://www.w3.org/2000/svg">` +
      (defs ? `<defs>${defs}</defs>` : '') + s + '</svg>';
  }

  /* ================================================================
     unrooted (equal angle)
     ================================================================ */
  function unrooted(tree, opts) {
    opts = defaults(opts);
    const W = opts.width, H = Math.round(opts.width * 0.8);
    const labels = opts.labels;
    const t = Tree.unroot(Tree.clone(tree));
    const total = Tree.totalLength(t) || 1;

    const pal = groupPalette(opts);
    const gOf = tp => (opts.groups ? opts.groups.of.get(tp) : null);
    const groupColour = g => (g == null || !pal.length ? V('text') : pal[g % pal.length]);
    const wantImg = opts.images && opts.images !== 'none';
    const imgSize = wantImg ? Math.max(18, +opts.imgSize || 54) : 0;

    /* Same rule as the circular layout: the pictures go on a ring outside, on a
       wider canvas, instead of eating the space the tree needs. */
    const extra = wantImg ? imgSize + 14 : 0;
    const scale = Math.min(W, H) * 0.34 / (total / Math.max(2, Tree.tips(t).length) * 3);
    const WW = W + 2 * extra, HH = H + 2 * extra;
    const cx = WW / 2, cy = HH / 2;

    let s = '';
    let defs = '';
    if (opts.title) s += txt(12, 18, opts.title, 13, 'start', 'art-txt', 'font-weight="600"');

    const colourOf = nd => {
      if (opts.colours) {
        const g = cladeGroup(nd, opts.colours);
        if (g != null) return V('c' + (1 + (g % 10)));
      }
      if (opts.groups && opts.colourBranches !== false && pal.length) {
        const g = window.Groups ? Groups.groupOfNode(nd, opts.groups) : null;
        if (g != null) return groupColour(g);
      }
      return V('text');
    };

    const pts = [];
    (function draw(nd, x, y, a0, a1) {
      const leaves = Tree.nTips(nd);
      let acc = a0;
      nd.children.forEach(c => {
        const share = (a1 - a0) * Tree.nTips(c.node) / leaves;
        const mid = acc + share / 2;
        const len = Math.max(2, (c.len || 0) * scale);
        const x2 = x + len * Math.cos(mid), y2 = y + len * Math.sin(mid);
        s += line(x, y, x2, y2, colourOf(c.node), opts.lineWidth);
        if (!c.node.children.length) {
          const g = gOf(c.node.tip);
          s += circ(x2, y2, 2.2, g != null && pal.length ? groupColour(g) : V('primary'));
          const name = (labels ? labels[c.node.tip] : c.node.label) || '';
          pts.push({ x: x2, y: y2, name, a: mid, tip: c.node.tip, node: c.node, g });
        } else draw(c.node, x2, y2, acc, acc + share);
        acc += share;
      });
    })(t, cx, cy, 0, Math.PI * 2);

    /* Nothing keeps two tips of an unrooted tree from landing beside each
       other, so their names land on each other too. A name may move in two
       directions here, so there is no one axis to spread along: the ones that
       still collide are dropped, the most peripheral kept first — those are
       the ones a reader can actually follow back to a branch. */
    let hiddenLabels = 0;
    const fs = opts.tipFont;
    const boxes = pts.map(p => {
      const flip = Math.cos(p.a) < 0;
      const w = p.name.length * fs * 0.52;
      const x = p.x + Math.cos(p.a) * 5 - (flip ? w : 0);
      return { x, y: p.y + Math.sin(p.a) * 5 - fs * 0.6, w, h: fs + 1,
        rank: -Math.hypot(p.x - cx, p.y - cy) };
    });
    const keep = opts.labelDeclutter === false ? boxes.map(() => true) : dropColliding(boxes);
    pts.forEach((p, i) => {
      if (!keep[i]) { hiddenLabels++; return; }
      const flip = Math.cos(p.a) < 0;
      const fill = opts.colourLabels !== false && p.g != null && pal.length ? ` fill="${groupColour(p.g)}"` : '';
      s += `<text x="${f1(p.x + Math.cos(p.a) * 5)}" y="${f1(p.y + Math.sin(p.a) * 5)}" font-size="${fs}" class="art-txt"${fill}
        ${opts.italicTips ? 'font-style="italic"' : ''} text-anchor="${flip ? 'end' : 'start'}" dominant-baseline="middle">${esc2(p.name)}</text>`;
    });

    /* the pictures, on a ring around the whole drawing */
    const imgs = collectImages(opts, {
      tips: pts.map(p => ({ tip: p.tip, label: p.name })),
      yOf: tp => { const p = pts.find(q => q.tip === tp.tip); return p ? Math.atan2(p.y - cy, p.x - cx) : 0; },
      labelOf: tp => tp.label,
      colourOf: tp => { const p = pts.find(q => q.tip === tp.tip); return p && p.g != null ? groupColour(p.g) : V('text'); },
      groupColour,
    });
    if (imgs.length) {
      const shape = opts.imgShape || 'round';
      let far = 0;
      pts.forEach(p => { far = Math.max(far, Math.hypot(p.x - cx, p.y - cy)); });
      const Rimg = Math.min(Math.min(WW, HH) / 2 - imgSize / 2 - 4, far + imgSize / 2 + 18);
      /* angles wrap, so they are spread over one full turn starting at the
         first picture; the gap is how much of the circle one picture occupies */
      const base = imgs[0].at;
      const rel = imgs.map(it => { let d = it.at - base; while (d < 0) d += 2 * Math.PI; return d; });
      const gapA = (imgSize + 5) / Math.max(1, Rimg);
      const placed = spreadLabels(rel, gapA, 0, 2 * Math.PI - gapA);
      imgs.forEach((it, i) => {
        defs += imageDefs(shape, imgSize, it.clip);
        const a = base + placed.pos[i];
        const px = cx + Rimg * Math.cos(a), py = cy + Rimg * Math.sin(a);
        const src = pts.find(q => q.name === it.name);
        if (src) it.linkTo = [src.x, src.y];
        s += imageNode(it, px, py, imgSize, shape, opts.imgFrame);
      });
    }

    if (opts.legend && opts.groups && opts.groups.k > 1) {
      let lx = 12; const ly = HH - 12;
      s += '<g data-role="legend">';
      for (let g = 0; g < opts.groups.k; g++) {
        const nm = (window.Groups ? Groups.nameOf(opts.groups, g) : String(g + 1));
        s += `<rect x="${f1(lx)}" y="${f1(ly - 9)}" width="10" height="10" rx="2" fill="${groupColour(g)}" data-li="${g}"/>`;
        s += txt(lx + 14, ly, nm, 10, 'start', 'art-mut', `data-li="${g}"`);
        lx += 14 + nm.length * 5.6 + 14;
      }
      s += '</g>';
    }

    TreeView.lastInfo = { hiddenSupport: 0, hiddenLabels, images: imgs.length,
      width: WW, height: HH, k: opts.groups ? opts.groups.k : 0 };
    return `<svg viewBox="0 0 ${WW} ${HH}" data-plot="${plotArea(0, 0, WW, HH)}" xmlns="http://www.w3.org/2000/svg">` +
      (defs ? `<defs>${defs}</defs>` : '') + s + '</svg>';
  }

  function render(tree, opts) {
    opts = opts || {};
    if (opts.layout === 'circular') return circular(tree, opts);
    if (opts.layout === 'unrooted') return unrooted(tree, opts);
    return rect(tree, opts);
  }

  Object.assign(TreeView, { render, rect, circular, unrooted, niceStep, cladeGroup,
    spreadLabels, dropColliding, collectImages, groupPalette, imageDefs, imageNode });
  window.TreeView = TreeView;
})();
