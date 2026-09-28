/* PhylogenyPro — molecular phylogenetics in the browser.
   Copyright (C) 2026  Luis Ángel Barrera-Guzmán

   This program is free software: you can redistribute it and/or modify it under
   the terms of the GNU General Public License as published by the Free Software
   Foundation, either version 3 of the License, or (at your option) any later
   version. It is distributed in the hope that it will be useful, but WITHOUT ANY
   WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
   A PARTICULAR PURPOSE. See the GNU General Public License, in the file LICENSE
   at the root of this program, or <https://www.gnu.org/licenses/>. */

/* PhylogenyPro — figures of Block 10: where the lineages were.

   A biogeographic figure has one job that no table does: putting the
   reconstruction next to the tree so the reader can see that the confident
   nodes and the doubtful ones are not the same nodes. Four pictures:

   rangeTree   the tree with a pie at every node over the ranges that have any
               probability, and the two corner boxes that say how the range was
               divided. Ranges are coloured by composition, so a widespread
               range is visibly a mixture of its areas rather than an
               arbitrary new colour.
   eventBars   the counts of dispersal, extinction and the four cladogenetic
               events, from the stochastic map, with their spread.
   matrixHeat  the dispersal multipliers between areas — the part of the model
               the user writes, and therefore the part worth showing back.
   modelBars   the Akaike weights of the six models, side by side.

   Colours come from the theme tokens, so everything follows light and dark. */

const Plots10 = {};

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

  /* one colour per area; a range takes the colour of its first area and a
     hatch of the rest, so widespread ranges read as mixtures */
  const AREA_COLOURS = ['c1', 'c2', 'c3', 'c4', 'c5', 'c6', 'c7', 'c8', 'c9', 'c10'];
  const areaColour = i => V(AREA_COLOURS[i % AREA_COLOURS.length]);

  /* a wedge of a pie, filled with the areas of its range in equal stripes */
  function wedge(cx, cy, r, a0, a1, areas, id) {
    if (!areas.length) return '';
    const large = (a1 - a0) > Math.PI ? 1 : 0;
    const p0 = [cx + r * Math.cos(a0), cy + r * Math.sin(a0)];
    const p1 = [cx + r * Math.cos(a1), cy + r * Math.sin(a1)];
    const d = `M${f1(cx)} ${f1(cy)} L${f1(p0[0])} ${f1(p0[1])} A${f1(r)} ${f1(r)} 0 ${large} 1 ${f1(p1[0])} ${f1(p1[1])} Z`;
    if (areas.length === 1) return `<path d="${d}" fill="${areaColour(areas[0])}"/>`;
    /* more than one area: a striped pattern, declared once per combination */
    return `<path d="${d}" fill="url(#${id})"/>`;
  }
  function stripePattern(id, areas) {
    const w = 4 * areas.length;
    let s = `<pattern id="${id}" width="${w}" height="${w}" patternUnits="userSpaceOnUse" patternTransform="rotate(45)">`;
    areas.forEach((a, i) => { s += `<rect x="${i * 4}" y="0" width="4" height="${w}" fill="${areaColour(a)}"/>`; });
    return s + '</pattern>';
  }

  /* ================================================================
     the tree with its ranges
     ================================================================ */
  /* nodes: [{node, probs}]; states: the Biogeo state table; corners optional */
  function rangeTree(tree, opts) {
    opts = opts || {};
    const states = opts.states, areaNames = opts.areaNames || [];
    const nA = areaNames.length;
    const labels = opts.labels;
    const probsOf = opts.probs;              // Map node -> array over states
    const tipRanges = opts.tipRanges || [];
    const corners = opts.corners || null;
    const minP = opts.minProbability == null ? 0.05 : opts.minProbability;

    const tips = Tree.tips(tree);
    const n = tips.length;
    const rowH = opts.rowHeight || (n > 45 ? 13 : 17);
    const pad = 14, tipFont = opts.tipFont || 10.5;
    const longest = tips.reduce((a, b) => {
      const la = String((labels ? labels[a.tip] : a.label) || '');
      const lb = String((labels ? labels[b.tip] : b.label) || '');
      return la.length >= lb.length ? a : b;
    }, tips[0]);
    const labW = Math.min(240, String((labels ? labels[longest.tip] : longest.label) || '').length * tipFont * 0.56 + 42);
    const W = opts.width || 780;
    const legendH = 24;
    const H = pad * 2 + n * rowH + 34 + legendH;
    const x0 = pad + 4, x1 = W - pad - labW - 8;
    const y0 = pad + legendH, y1 = H - pad - 34;

    const F = Tree.flatten(tree);
    const depth = new Float64Array(F.n);
    for (let i = F.post.length - 1; i >= 0; i--) {
      const k = F.post[i], p = F.parent[k];
      if (p >= 0) depth[k] = depth[p] + (opts.cladogram ? 1 : (F.len[k] || 0));
    }
    let maxD = 0;
    for (let k = 0; k < F.n; k++) if (F.isTip[k]) maxD = Math.max(maxD, depth[k]);
    const X = d => x0 + (x1 - x0) * (maxD > 0 ? d / maxD : 0);
    const yOf = new Float64Array(F.n);
    let row = 0;
    for (let i = 0; i < F.post.length; i++) {
      const k = F.post[i];
      if (F.isTip[k]) { yOf[k] = y0 + rowH / 2 + (row++) * rowH; continue; }
      const ys = F.kids[k].map(c => yOf[c]);
      yOf[k] = (Math.min.apply(null, ys) + Math.max.apply(null, ys)) / 2;
    }

    /* every range that will actually be drawn, so the patterns are declared once */
    const used = new Set();
    const noteRange = mask => { if (mask > 0) used.add(mask); };
    tipRanges.forEach(noteRange);
    if (probsOf) {
      for (let k = 0; k < F.n; k++) {
        const p = probsOf.get ? probsOf.get(k) : probsOf[k];
        if (!p) continue;
        p.forEach((v, i) => { if (v >= minP) noteRange(states.list[i]); });
      }
    }
    const patternId = mask => 'p10r' + mask;
    let defs = '<defs>';
    used.forEach(mask => {
      const areas = Biogeo.areasOf(mask, nA);
      if (areas.length > 1) defs += stripePattern(patternId(mask), areas);
    });
    defs += '</defs>';

    let s = defs;
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
    /* the tips, with their observed range */
    for (let k = 0; k < F.n; k++) {
      if (!F.isTip[k]) continue;
      const mask = tipRanges[F.tipRow[k]];
      const areas = mask ? Biogeo.areasOf(mask, nA) : [];
      areas.forEach((a, i) => {
        s += rect(X(depth[k]) + 3 + i * 7, yOf[k] - 3.6, 6, 7.2, areaColour(a), 'rx="1.3"');
      });
      const lab = labels ? labels[F.tipRow[k]] : (F.nodes[k].label || '');
      s += txt(X(depth[k]) + 6 + Math.max(1, areas.length) * 7, yOf[k] + 3.4, lab, 'art-txt', tipFont, 'start', 'font-style="italic"');
    }
    /* the corners: small boxes on each daughter branch, just after the split */
    if (corners && opts.showCorners !== false) {
      for (let k = 0; k < F.n; k++) {
        const c = corners[k];
        if (!c) continue;
        [[c.left, c.leftProbs], [c.right, c.rightProbs]].forEach(([child, pr]) => {
          const i = pr.indexOf(Math.max.apply(null, pr));
          if (pr[i] < minP) return;
          const mask = states.list[i];
          const cx = X(depth[k]) + 5, cy = yOf[child] + (yOf[child] < yOf[k] ? -6 : 6);
          const areas = Biogeo.areasOf(mask, nA);
          s += rect(cx, cy - 3, Math.max(5, areas.length * 4.6), 6,
            areas.length > 1 ? `url(#${patternId(mask)})` : areaColour(areas[0] || 0),
            `rx="1.2" opacity="${(0.35 + 0.65 * pr[i]).toFixed(2)}"`);
        });
      }
    }
    /* the pies */
    const r = opts.pieRadius || 6;
    if (probsOf) {
      for (let k = 0; k < F.n; k++) {
        if (F.isTip[k]) continue;
        const pr = probsOf.get ? probsOf.get(k) : probsOf[k];
        if (!pr) continue;
        const shown = pr.map((v, i) => [i, v]).filter(x => x[1] >= minP).sort((a, b) => b[1] - a[1]);
        const shownSum = shown.reduce((a, b) => a + b[1], 0);
        const rest = 1 - shownSum;
        const cx = X(depth[k]), cy = yOf[k];
        let a0 = -Math.PI / 2;
        shown.forEach(([i, v]) => {
          const a1 = a0 + 2 * Math.PI * v;
          const areas = Biogeo.areasOf(states.list[i], nA);
          if (v > 0.9995) s += circ(cx, cy, r, areas.length > 1 ? `url(#${patternId(states.list[i])})` : areaColour(areas[0] || 0));
          else s += wedge(cx, cy, r, a0, a1, areas, patternId(states.list[i]));
          a0 = a1;
        });
        if (rest > 1e-6) {
          const a1 = a0 + 2 * Math.PI * rest;
          const large = (a1 - a0) > Math.PI ? 1 : 0;
          const p0 = [cx + r * Math.cos(a0), cy + r * Math.sin(a0)];
          const p1 = [cx + r * Math.cos(a1), cy + r * Math.sin(a1)];
          s += `<path d="M${f1(cx)} ${f1(cy)} L${f1(p0[0])} ${f1(p0[1])} A${f1(r)} ${f1(r)} 0 ${large} 1 ${f1(p1[0])} ${f1(p1[1])} Z" fill="${V('border-strong')}" opacity="0.45"/>`;
        }
        s += circ(cx, cy, r, 'none', `stroke="${V('bg')}" stroke-width="0.9"`);
      }
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
    /* the legend of areas */
    let lx = pad + 4;
    areaNames.forEach((nm, i) => {
      s += rect(lx, pad - 2, 9, 9, areaColour(i), 'rx="1.6"');
      s += txt(lx + 13, pad + 6, nm, 'art-txt', 10);
      lx += 24 + String(nm).length * 6.2;
    });
    if (lx < W - 120) {
      s += rect(lx, pad - 2, 9, 9, V('border-strong'), 'rx="1.6" opacity="0.45"');
      s += txt(lx + 13, pad + 6, opts.restLabel || 'resto', 'art-mut', 10);
    }
    return svg(`0 0 ${W} ${H}`, s);
  }

  /* ================================================================
     the events a stochastic map counted
     ================================================================ */
  function eventBars(mean, sd, labels, opts) {
    opts = opts || {};
    const keys = Object.keys(mean);
    const W = opts.width || 620, H = opts.height || 44 + keys.length * 26;
    const mx = opts.labelWidth || 140, right = 60, top = 16;
    const maxV = Math.max.apply(null, keys.map(k => mean[k] + (sd ? sd[k] : 0))) || 1;
    const X = v => mx + (W - mx - right) * v / maxV;
    let s = '';
    s += line(mx, top - 6, mx, top + keys.length * 26 - 4, V('border-strong'), 1.2);
    keys.forEach((k, i) => {
      const y = top + i * 26;
      s += rect(mx, y, X(mean[k]) - mx, 14, V(i < 2 ? 'c5' : 'c2'), 'rx="2" opacity="0.8"');
      if (sd && sd[k] > 0) {
        s += line(X(mean[k] - sd[k]), y + 7, X(mean[k] + sd[k]), y + 7, V('border-strong'), 1.3);
        s += line(X(mean[k] - sd[k]), y + 3, X(mean[k] - sd[k]), y + 11, V('border-strong'), 1.1);
        s += line(X(mean[k] + sd[k]), y + 3, X(mean[k] + sd[k]), y + 11, V('border-strong'), 1.1);
      }
      s += txt(mx - 8, y + 11, (labels && labels[k]) || k, 'art-txt', 10.5, 'end');
      s += txt(X(mean[k]) + 6, y + 11, mean[k].toFixed(2) + (sd ? ` ± ${sd[k].toFixed(2)}` : ''), 'art-mut', 9.5);
    });
    const step = niceStep(maxV, 5);
    for (let v = 0; v <= maxV; v += step) {
      const y = top + keys.length * 26 - 2;
      s += line(X(v), y, X(v), y + 4, V('border-strong'), 1);
      s += txt(X(v), y + 15, fmtTick(v, step), 'art-mut', 9, 'middle');
    }
    return svg(`0 0 ${W} ${H + 10}`, s);
  }

  /* ================================================================
     the dispersal multipliers
     ================================================================ */
  function matrixHeat(M, areaNames, opts) {
    opts = opts || {};
    const n = areaNames.length;
    const cell = opts.cell || 34;
    const left = 74, top = 40;
    const W = left + n * cell + 20, H = top + n * cell + 26;
    let lo = Infinity, hi = -Infinity;
    for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) { lo = Math.min(lo, M[i][j]); hi = Math.max(hi, M[i][j]); }
    if (hi - lo < 1e-12) { lo = 0; hi = Math.max(1, hi); }
    let s = '';
    for (let i = 0; i < n; i++) {
      s += txt(left - 8, top + i * cell + cell / 2 + 4, areaNames[i], 'art-txt', 10.5, 'end');
      s += txt(left + i * cell + cell / 2, top - 10, areaNames[i], 'art-txt', 10.5, 'middle');
      for (let j = 0; j < n; j++) {
        const v = M[i][j];
        const t = (v - lo) / (hi - lo);
        s += rect(left + j * cell, top + i * cell, cell - 1.5, cell - 1.5,
          i === j ? V('bg-soft') : V('primary'), `rx="3" opacity="${i === j ? 0.5 : (0.12 + 0.8 * t).toFixed(3)}"`);
        if (i !== j) s += txt(left + j * cell + cell / 2 - 0.75, top + i * cell + cell / 2 + 4,
          v === 1 ? '1' : (+v).toFixed(2).replace(/0+$/, '').replace(/\.$/, ''),
          t > 0.55 ? 'art-txt' : 'art-mut', 9.5, 'middle', t > 0.55 ? `fill="${V('bg')}"` : '');
      }
    }
    s += txt(left, 16, opts.title || '', 'art-mut', 10);
    s += txt(left, H - 6, opts.caption || '', 'art-mut', 9);
    return svg(`0 0 ${W} ${H}`, s);
  }

  /* ================================================================
     the models, by Akaike weight
     ================================================================ */
  function modelBars(rows, opts) {
    opts = opts || {};
    const W = opts.width || 620, H = 30 + rows.length * 24;
    const mx = 120, right = 120, top = 12;
    let s = '';
    rows.forEach((r, i) => {
      const y = top + i * 24;
      const w = (W - mx - right) * r.w;
      s += rect(mx, y, Math.max(1, w), 13, V(r.withJ ? 'c2' : 'c1'), 'rx="2" opacity="0.82"');
      s += txt(mx - 8, y + 11, r.name, 'art-txt', 10.5, 'end');
      s += txt(mx + Math.max(1, w) + 6, y + 11,
        `w ${r.w.toFixed(3)}   ΔAICc ${r.dAICc.toFixed(2)}`, 'art-mut', 9.5);
    });
    s += line(mx, top - 4, mx, top + rows.length * 24 - 6, V('border-strong'), 1.2);
    return svg(`0 0 ${W} ${H}`, s);
  }

  Object.assign(Plots10, { rangeTree, eventBars, matrixHeat, modelBars, areaColour, AREA_COLOURS });
  window.Plots10 = Plots10;
})();
