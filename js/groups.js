/* PhylogenyPro — molecular phylogenetics in the browser.
   Copyright (C) 2026  Luis Ángel Barrera-Guzmán

   This program is free software: you can redistribute it and/or modify it under
   the terms of the GNU General Public License as published by the Free Software
   Foundation, either version 3 of the License, or (at your option) any later
   version. It is distributed in the hope that it will be useful, but WITHOUT ANY
   WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
   A PARTICULAR PURPOSE. See the GNU General Public License, in the file LICENSE
   at the root of this program, or <https://www.gnu.org/licenses/>. */

/* PhylogenyPro — cutting a tree into k groups, and colouring by them.

   Three ways of getting k groups, because they answer different questions:

     · longest branches — remove the k−1 longest internal edges and take the
       pieces that fall apart. Works on any tree, including one with very
       unequal rates, and is what the eye does when it says "these go together".
     · by depth — the threshold on node age that leaves exactly k clades. Only
       honest on an ultrametric tree (a chronogram, or UPGMA), so it refuses
       elsewhere and says why.
     · by hand — the user names the groups and drags taxa into them.

   A group is never imposed on the figure: the cut is a hypothesis drawn over
   the tree, and the panel says which of the three produced it.

   Groups are plain data — {of: Map(tipIndex → g), k, names, method} — so the
   drawing code and the report can both read them without knowing how they
   were obtained. */

const Groups = {};

(function () {

  /* Ten colours that stay apart both in hue and in lightness, so a group is
     still a group in greyscale and for the most common colour blindness.
     Beyond ten the hues are walked round the circle. */
  const BASE = [
    '#2a4a94', '#b76a17', '#227c60', '#8e3b6c', '#2d84b8',
    '#9a7b12', '#6b4bab', '#1f6f7a', '#a8432f', '#4f6b23',
  ];

  function palette(k) {
    const out = [];
    for (let i = 0; i < k; i++) {
      if (i < BASE.length) { out.push(BASE[i]); continue; }
      const t = (i - BASE.length) / Math.max(1, k - BASE.length);
      out.push(`hsl(${Math.round(360 * ((t * 0.618) % 1))} 52% ${i % 2 ? 38 : 56}%)`);
    }
    return out;
  }

  /* ================================================================
     1 · the pieces a set of cut edges leaves behind
     ================================================================ */
  /* Every node gets the group of the nearest ancestor edge that was NOT cut;
     crossing a cut edge starts a new group. Tips carry the answer. */
  function componentsAfterCutting(root, cutSet) {
    const of = new Map();
    let next = 0;
    (function walk(nd, g) {
      if (Tree.isTip(nd)) { of.set(nd.tip, g); return; }
      nd.children.forEach(c => walk(c.node, cutSet.has(c.node) ? ++next : g));
    })(root, 0);
    return of;
  }

  /* Renumber so that group 0 is the one whose first tip comes first in the
     drawing order — otherwise the legend jumps about between redraws. */
  function renumber(of, order) {
    const seen = new Map();
    let n = 0;
    order.forEach(t => { const g = of.get(t); if (g != null && !seen.has(g)) seen.set(g, n++); });
    const out = new Map();
    of.forEach((g, t) => out.set(t, seen.has(g) ? seen.get(g) : 0));
    return { of: out, k: n };
  }

  /* ================================================================
     2 · by the longest internal branches
     ================================================================ */
  function byLongestBranches(root, k) {
    const order = Tree.tips(root).map(t => t.tip);
    if (k <= 1) return { of: new Map(order.map(t => [t, 0])), k: 1, method: 'longest' };

    /* The two branches hanging off a bifurcating root are ONE branch of the
       unrooted tree: cutting either of them splits the tree into the same two
       halves, so counting them separately spends two cuts to buy one group.
       They are merged into a single edge whose length is their sum — which is
       also the length the unrooted tree actually has there. */
    const edges = [];
    const rootPair = root.children.length === 2;
    (function walk(nd) {
      nd.children.forEach(c => {
        if (!(nd === root && rootPair)) edges.push({ node: c.node, len: c.len || 0 });
        walk(c.node);
      });
    })(root);
    if (rootPair) {
      edges.push({
        node: root.children[0].node,
        len: (root.children[0].len || 0) + (root.children[1].len || 0),
      });
    }

    edges.sort((a, b) => b.len - a.len);

    /* Cutting both branches of an internal node leaves that node with no tips
       at all: a component that exists in the graph and not in the data. So the
       cuts are taken one at a time, longest first, and one that does not
       actually add a group is skipped rather than wasted. */
    const cut = new Set();
    let have = 1;
    for (const e of edges) {
      if (have >= k) break;
      cut.add(e.node);
      const got = renumber(componentsAfterCutting(root, cut), order).k;
      if (got > have) have = got; else cut.delete(e.node);
    }
    const r = renumber(componentsAfterCutting(root, cut), order);
    return {
      of: r.of, k: r.k, method: 'longest',
      exact: r.k === k,
      cuts: [...cut].length,
    };
  }

  /* ================================================================
     3 · by depth (only on an ultrametric tree)
     ================================================================ */
  /* depth from the root to every node */
  function depths(root) {
    const d = new Map();
    (function walk(nd, x) { d.set(nd, x); nd.children.forEach(c => walk(c.node, x + (c.len || 0))); })(root, 0);
    return d;
  }

  function isUltrametric(root, tol) {
    const d = depths(root);
    const t = Tree.tips(root).map(n => d.get(n));
    if (!t.length) return false;
    const mx = Math.max.apply(null, t), mn = Math.min.apply(null, t);
    return mx - mn <= (tol == null ? 1e-6 : tol) * Math.max(1e-12, mx);
  }

  /* The threshold that leaves exactly k groups: cut every edge that crosses it.
     Walking the internal depths in order gives every attainable k. */
  function byDepth(root, k) {
    const order = Tree.tips(root).map(t => t.tip);
    if (k <= 1) return { of: new Map(order.map(t => [t, 0])), k: 1, method: 'depth' };
    const d = depths(root);
    const internal = Tree.nodes(root).filter(n => !Tree.isTip(n) && n !== root)
      .map(n => d.get(n)).sort((a, b) => a - b);
    /* candidate thresholds: just above each internal node's depth */
    const cand = [0].concat(internal);
    let best = null;
    for (const h of cand) {
      const cut = new Set();
      (function walk(nd) {
        nd.children.forEach(c => {
          if (d.get(nd) < h + 1e-12 && d.get(c.node) >= h - 1e-12 && d.get(c.node) > d.get(nd)) cut.add(c.node);
          else walk(c.node);
        });
      })(root);
      const r = renumber(componentsAfterCutting(root, cut), order);
      if (best === null || Math.abs(r.k - k) < Math.abs(best.k - k)) best = { of: r.of, k: r.k, height: h };
      if (r.k === k) break;
    }
    return { of: best.of, k: best.k, method: 'depth', height: best.height, exact: best.k === k };
  }

  /* ================================================================
     4 · from clades the user picked, and by hand
     ================================================================ */
  /* clades: an array of arrays of tip indices. Anything left over is group 0
     ("the rest"), which is what an outgroup usually is. */
  function fromClades(root, clades) {
    const order = Tree.tips(root).map(t => t.tip);
    const of = new Map(order.map(t => [t, -1]));
    clades.forEach((set, i) => set.forEach(t => { if (of.has(t)) of.set(t, i); }));
    let rest = false;
    of.forEach(v => { if (v < 0) rest = true; });
    const shift = rest ? 1 : 0;
    const out = new Map();
    of.forEach((v, t) => out.set(t, v < 0 ? 0 : v + shift));
    return { of: out, k: clades.length + shift, method: 'clades' };
  }

  function empty(root) {
    const order = Tree.tips(root).map(t => t.tip);
    return { of: new Map(order.map(t => [t, 0])), k: 1, method: 'manual' };
  }

  function assign(groups, tip, g) {
    const of = new Map(groups.of);
    of.set(tip, g);
    let k = 0; of.forEach(v => { if (v + 1 > k) k = v + 1; });
    return { of, k: Math.max(k, groups.k), method: 'manual', names: groups.names };
  }

  /* ================================================================
     5 · reading a grouping
     ================================================================ */
  function sizes(groups) {
    const s = new Array(groups.k).fill(0);
    groups.of.forEach(g => { if (g >= 0 && g < groups.k) s[g]++; });
    return s;
  }

  function nameOf(groups, g) {
    if (groups.names && groups.names[g]) return groups.names[g];
    return String.fromCharCode(65 + (g % 26)) + (g >= 26 ? Math.floor(g / 26) : '');
  }

  /* Is every group a clade? A cut by branches always is; a grouping by hand
     may not be, and a figure that colours a non-monophyletic group must say so
     rather than let the reader assume it. */
  function monophyly(root, groups) {
    const out = [];
    for (let g = 0; g < groups.k; g++) {
      const want = new Set();
      groups.of.forEach((v, t) => { if (v === g) want.add(t); });
      if (!want.size) { out.push({ group: g, size: 0, monophyletic: true }); continue; }
      let found = false;
      Tree.nodes(root).forEach(nd => {
        if (found) return;
        const here = Tree.tips(nd).map(t => t.tip);
        if (here.length === want.size && here.every(t => want.has(t))) found = true;
      });
      out.push({ group: g, size: want.size, monophyletic: found });
    }
    return out;
  }

  /* the group of an internal node, when all of its tips agree; otherwise null */
  function groupOfNode(nd, groups) {
    const t = Tree.tips(nd).map(x => x.tip);
    if (!t.length) return null;
    const g = groups.of.get(t[0]);
    return t.every(x => groups.of.get(x) === g) ? g : null;
  }

  Object.assign(Groups, {
    BASE, palette,
    byLongestBranches, byDepth, fromClades, empty, assign,
    sizes, nameOf, monophyly, groupOfNode, isUltrametric, depths,
    componentsAfterCutting,
  });
  if (typeof window !== 'undefined') window.Groups = Groups;

})();
