/* PhylogenyPro — molecular phylogenetics in the browser.
   Copyright (C) 2026  Luis Ángel Barrera-Guzmán

   This program is free software: you can redistribute it and/or modify it under
   the terms of the GNU General Public License as published by the Free Software
   Foundation, either version 3 of the License, or (at your option) any later
   version. It is distributed in the hope that it will be useful, but WITHOUT ANY
   WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
   A PARTICULAR PURPOSE. See the GNU General Public License, in the file LICENSE
   at the root of this program, or <https://www.gnu.org/licenses/>. */

/* PhylogenyPro — trees: structure, Newick, distance methods and the operations
   every block needs.

   Node format, shared with the whole app (and with PopGeneticsPro, so trees can
   travel between them):

     { tip: index | undefined, label, children: [{ node, len }], support, age }

   Parents are never stored, so a tree can be cloned and posted to a Worker;
   they are rebuilt on demand with Tree.index(). */

/* Wrapped in a named function so that a Web Worker can be built from its own
   source text: see js/pool.js. Nothing here touches the DOM. */
function TreeCore(g) {
const Tree = {};

  const isTip = n => n.tip != null;

  function tips(node, out) { out = out || []; if (isTip(node)) out.push(node); node.children.forEach(c => tips(c.node, out)); return out; }
  function nodes(node, out) { out = out || []; out.push(node); node.children.forEach(c => nodes(c.node, out)); return out; }
  function nTips(node) { return isTip(node) ? 1 : node.children.reduce((a, c) => a + nTips(c.node), 0); }

  function index(root) {
    const parent = new Map(), lenUp = new Map(), depth = new Map(), order = [];
    (function walk(n, d) {
      depth.set(n, d); order.push(n);
      n.children.forEach(c => { parent.set(c.node, n); lenUp.set(c.node, c.len); walk(c.node, d + (c.len || 0)); });
    })(root, 0);
    return { parent, lenUp, depth, order, postorder: order.slice().reverse() };
  }

  function clone(node) {
    const out = Object.assign({}, node, { children: node.children.map(c => ({ node: clone(c.node), len: c.len })) });
    return out;
  }

  /* ================================================================
     Newick
     ================================================================ */
  function parseNewick(text, labels) {
    const s = String(text).trim().replace(/;\s*$/, '');
    let i = 0;
    const byName = labels ? new Map(labels.map((l, k) => [l, k])) : null;
    let nextTip = 0;

    function skipSpace() { while (i < s.length && /\s/.test(s[i])) i++; }
    function readLabel() {
      skipSpace();
      if (s[i] === "'" || s[i] === '"') {
        const q = s[i++];
        let out = '';
        while (i < s.length && s[i] !== q) out += s[i++];
        i++;
        return out;
      }
      let out = '';
      while (i < s.length && !/[(),:;\[\]]/.test(s[i])) out += s[i++];
      return out.trim().replace(/_/g, ' ').trim() === '' ? out.trim() : out.trim();
    }
    function readNumber() {
      skipSpace();
      const m = s.slice(i).match(/^[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?/);
      if (!m) return null;
      i += m[0].length;
      return parseFloat(m[0]);
    }
    function skipComment() {
      if (s[i] === '[') { let d = 0; do { if (s[i] === '[') d++; else if (s[i] === ']') d--; i++; } while (i < s.length && d > 0); }
    }
    function node() {
      skipSpace(); skipComment(); skipSpace();
      const n = { children: [] };
      if (s[i] === '(') {
        i++;
        for (;;) {
          const child = node();
          skipSpace(); skipComment(); skipSpace();
          let len = 0;
          if (s[i] === ':') { i++; skipComment(); skipSpace(); len = readNumber(); }
          /* an annotation may follow the length, as FigTree and BEAST write it:
             (A:0.1[&rate=1],B:0.2) — without skipping it here the comment would
             be read as another taxon */
          skipSpace(); skipComment(); skipSpace();
          n.children.push({ node: child, len: len == null ? 0 : len });
          if (s[i] === ',') { i++; continue; }
          if (s[i] === ')') { i++; break; }
          if (i >= s.length) break;
          i++;                       // tolerate stray characters
        }
        skipComment();
        const lab = readLabel();
        if (lab !== '') {
          const v = parseFloat(lab);
          if (isFinite(v)) n.support = v; else n.label = lab;
        }
      } else {
        const lab = readLabel();
        n.label = lab;
        n.tip = byName && byName.has(lab) ? byName.get(lab) : nextTip;
        if (!byName || !byName.has(lab)) nextTip++;
      }
      skipComment();
      return n;
    }
    const root = node();
    /* a child that carries its own branch length after the root */
    skipSpace();
    if (s[i] === ':') { i++; readNumber(); }
    return root;
  }

  function writeNewick(root, opts) {
    opts = opts || {};
    const dec = opts.decimals == null ? 6 : opts.decimals;
    const name = n => {
      let l = opts.labels && n.tip != null ? opts.labels[n.tip] : (n.label || '');
      l = String(l);
      return /[\s(),:;\[\]']/.test(l) ? "'" + l.replace(/'/g, "''") + "'" : l;
    };
    function rec(n) {
      if (isTip(n)) return name(n);
      const inner = n.children.map(c => rec(c.node) + (opts.lengths === false ? '' : ':' + (+c.len || 0).toFixed(dec))).join(',');
      let lab = '';
      if (opts.support !== false && n.support != null) lab = String(+(+n.support).toFixed(opts.supportDecimals == null ? 2 : opts.supportDecimals));
      else if (n.label) lab = name(n);
      return '(' + inner + ')' + lab;
    }
    return rec(root) + ';';
  }

  /* ================================================================
     splits (bipartitions) and tree comparison
     ================================================================ */
  /* every internal branch of an unrooted tree defines a split of the tips;
     a split is stored as a sorted string of the smaller side */
  function splits(root, nTaxa) {
    const out = new Map();
    const n = nTaxa || tips(root).length;
    (function walk(node2) {
      if (isTip(node2)) return [node2.tip];
      let set = [];
      node2.children.forEach(c => { set = set.concat(walk(c.node)); });
      if (set.length > 1 && set.length < n) {
        const a = set.slice().sort((x, y) => x - y);
        const b = [];
        for (let k = 0; k < n; k++) if (a.indexOf(k) < 0) b.push(k);
        const key = (a.length < b.length || (a.length === b.length && a[0] < b[0]) ? a : b).join(',');
        out.set(key, node2);
      }
      return set;
    })(root);
    return out;
  }
  /* Robinson–Foulds distance: splits in one tree and not in the other */
  function rfDistance(t1, t2, nTaxa) {
    const a = splits(t1, nTaxa), b = splits(t2, nTaxa);
    let diff = 0;
    a.forEach((v, k) => { if (!b.has(k)) diff++; });
    b.forEach((v, k) => { if (!a.has(k)) diff++; });
    return { rf: diff, maxRF: a.size + b.size, shared: a.size - [...a.keys()].filter(k => !b.has(k)).length, n1: a.size, n2: b.size };
  }

  /* ================================================================
     rooting and shaping
     ================================================================ */
  function totalLength(root) { return nodes(root).reduce((s, n) => s + n.children.reduce((a, c) => a + (c.len || 0), 0), 0); }

  /* Reroot on the branch above `target`, splitting its length at `pos` (0..1).

     Everything below the target keeps its shape; the rest of the tree is
     rebuilt walking upwards, turning each parent edge around. A node that ends
     up with a single child had degree two in the unrooted tree — the old root
     is the usual case — and is spliced out, its two branch lengths added, so
     that no artificial node of degree two survives. */
  function rerootAbove(root, target, pos) {
    pos = pos == null ? 0.5 : pos;
    const idx = index(root);
    const parent = idx.parent.get(target);
    if (!parent) return clone(root);
    const lenT = idx.lenUp.get(target) || 0;

    /* the tree as seen from `node`, without the branch towards `from` */
    function up(node, from) {
      const out = { children: [], support: node.support, label: node.label };
      node.children.forEach(c => { if (c.node !== from) out.children.push({ node: clone(c.node), len: c.len || 0 }); });
      const p = idx.parent.get(node);
      if (p) {
        const r = up(p, node);
        out.children.push({ node: r.node, len: (idx.lenUp.get(node) || 0) + (r.extra || 0) });
      }
      if (out.children.length === 1) return { node: out.children[0].node, extra: out.children[0].len };
      return { node: out, extra: 0 };
    }
    const rest = up(parent, target);
    return {
      children: [
        { node: clone(target), len: lenT * pos },
        { node: rest.node, len: lenT * (1 - pos) + (rest.extra || 0) },
      ],
    };
  }

  /* root with an outgroup: the branch above the MRCA of the outgroup tips */
  function rootByOutgroup(root, outgroupTips) {
    const want = new Set(outgroupTips);
    if (!want.size) return clone(root);
    let best = null;
    (function walk(n) {
      const set = isTip(n) ? [n.tip] : n.children.reduce((a, c) => a.concat(walk(c.node)), []);
      const inSet = set.filter(t => want.has(t)).length;
      if (inSet === want.size && set.length === want.size && !best) best = n;
      return set;
    })(root);
    if (!best) return { tree: clone(root), monophyletic: false };
    return { tree: rerootAbove(root, best, 0.5), monophyletic: true };
  }

  /* midpoint rooting: the middle of the longest path between two tips */
  function midpointRoot(root) {
    const idx = index(root);
    const all = tips(root);
    /* distance from every node to every tip below it, by post-order */
    const down = new Map();   // node -> [{tip, dist}]
    idx.postorder.forEach(n => {
      if (isTip(n)) { down.set(n, [{ tip: n, dist: 0 }]); return; }
      let list = [];
      n.children.forEach(c => { down.get(c.node).forEach(x => list.push({ tip: x.tip, dist: x.dist + (c.len || 0) })); });
      down.set(n, list);
    });
    /* the longest path passes through some node; find it */
    let bestPair = null, bestD = -1;
    nodes(root).forEach(n => {
      if (isTip(n)) return;
      const lists = n.children.map(c => down.get(c.node).map(x => ({ tip: x.tip, dist: x.dist + (c.len || 0) })));
      for (let a = 0; a < lists.length; a++) for (let b = a + 1; b < lists.length; b++) {
        lists[a].forEach(x => lists[b].forEach(y => {
          if (x.dist + y.dist > bestD) { bestD = x.dist + y.dist; bestPair = [x, y]; }
        }));
      }
    });
    if (!bestPair) return clone(root);
    /* walk from the deeper tip towards the root until half the path is behind */
    const [x, y] = bestPair;
    const half = bestD / 2;
    let cur = (x.dist >= y.dist ? x.tip : y.tip), acc = 0;
    while (idx.parent.has(cur)) {
      const l = idx.lenUp.get(cur) || 0;
      if (acc + l >= half) {
        /* the part of this branch that stays on the side of the deeper tip is
           half − acc, and rerootAbove measures `pos` from that same side */
        const pos = l > 0 ? (half - acc) / l : 0.5;
        return rerootAbove(root, cur, pos);
      }
      acc += l;
      cur = idx.parent.get(cur);
    }
    return clone(root);
  }

  /* ladderise: children sorted by how many tips they carry */
  function ladderize(root, ascending) {
    (function walk(n) {
      n.children.forEach(c => walk(c.node));
      n.children.sort((a, b) => (ascending ? 1 : -1) * (nTips(a.node) - nTips(b.node)));
    })(root);
    return root;
  }
  /* branches shorter than eps become polytomies */
  function collapseShort(root, eps) {
    (function walk(n) {
      const kids = [];
      n.children.forEach(c => {
        walk(c.node);
        if (!isTip(c.node) && (c.len || 0) <= eps) c.node.children.forEach(g => kids.push({ node: g.node, len: g.len }));
        else kids.push(c);
      });
      n.children = kids;
    })(root);
    return root;
  }
  /* support values below a threshold become polytomies */
  function collapseSupport(root, minSupport) {
    (function walk(n) {
      const kids = [];
      n.children.forEach(c => {
        walk(c.node);
        if (!isTip(c.node) && c.node.support != null && c.node.support < minSupport) c.node.children.forEach(g => kids.push({ node: g.node, len: (g.len || 0) + (c.len || 0) }));
        else kids.push(c);
      });
      n.children = kids;
    })(root);
    return root;
  }

  /* ================================================================
     distance methods
     ================================================================ */
  /* Neighbour joining (Saitou & Nei 1987), with the O(n³) formulation of
     Studier & Keppler (1988). D is a square array of arrays. */
  function nj(D0, labels) {
    const n0 = D0.length;
    if (n0 < 2) return { tip: 0, label: labels ? labels[0] : '0', children: [] };
    const D = D0.map(r => Array.from(r));
    const active = [];
    for (let i = 0; i < n0; i++) active.push({ node: { tip: i, label: labels ? labels[i] : String(i), children: [] }, idx: i });
    const alive = new Set(active.map(a => a.idx));
    const nodeOf = new Map(active.map(a => [a.idx, a.node]));
    let next = n0;
    while (alive.size > 2) {
      const ids = [...alive];
      const r = ids.length;
      const S = new Map();
      ids.forEach(i => { let s = 0; ids.forEach(j => { if (i !== j) s += D[i][j]; }); S.set(i, s); });
      let bi = ids[0], bj = ids[1], best = Infinity;
      for (let a = 0; a < ids.length; a++) for (let b = a + 1; b < ids.length; b++) {
        const i = ids[a], j = ids[b];
        const q = (r - 2) * D[i][j] - S.get(i) - S.get(j);
        if (q < best) { best = q; bi = i; bj = j; }
      }
      const dij = D[bi][bj];
      let li = 0.5 * dij + (S.get(bi) - S.get(bj)) / (2 * (r - 2));
      let lj = dij - li;
      /* negative branches are set to zero and the difference passed to the sister,
         which is what every implementation does */
      if (li < 0) { lj -= li; li = 0; }
      if (lj < 0) { li -= lj; lj = 0; }
      const parent = { children: [{ node: nodeOf.get(bi), len: li }, { node: nodeOf.get(bj), len: lj }] };
      const k = next++;
      nodeOf.set(k, parent);
      /* distances from the new node */
      D[k] = [];
      [...alive].forEach(m => {
        if (m === bi || m === bj) return;
        const d = 0.5 * (D[bi][m] + D[bj][m] - dij);
        D[k][m] = d;
        if (!D[m]) D[m] = [];
        D[m][k] = d;
      });
      D[k][k] = 0;
      alive.delete(bi); alive.delete(bj); alive.add(k);
    }
    const [a, b] = [...alive];
    return { children: [{ node: nodeOf.get(a), len: Math.max(0, D[a][b] / 2) }, { node: nodeOf.get(b), len: Math.max(0, D[a][b] / 2) }] };
  }

  /* BIONJ (Gascuel 1997): like NJ, but the new distances are a weighted average
     whose weights come from the variances of the estimates */
  function bionj(D0, labels) {
    const n0 = D0.length;
    if (n0 < 3) return nj(D0, labels);
    const D = D0.map(r => Array.from(r));
    const V = D0.map(r => Array.from(r));         // variance matrix, starts equal to D
    const alive = new Set();
    const nodeOf = new Map();
    for (let i = 0; i < n0; i++) { alive.add(i); nodeOf.set(i, { tip: i, label: labels ? labels[i] : String(i), children: [] }); }
    let next = n0;
    while (alive.size > 2) {
      const ids = [...alive], r = ids.length;
      const S = new Map();
      ids.forEach(i => { let s = 0; ids.forEach(j => { if (i !== j) s += D[i][j]; }); S.set(i, s); });
      let bi = ids[0], bj = ids[1], best = Infinity;
      for (let a = 0; a < ids.length; a++) for (let b = a + 1; b < ids.length; b++) {
        const i = ids[a], j = ids[b];
        const q = (r - 2) * D[i][j] - S.get(i) - S.get(j);
        if (q < best) { best = q; bi = i; bj = j; }
      }
      const dij = D[bi][bj];
      let li = 0.5 * dij + (S.get(bi) - S.get(bj)) / (2 * (r - 2));
      let lj = dij - li;
      if (li < 0) { lj -= li; li = 0; }
      if (lj < 0) { li -= lj; lj = 0; }
      /* lambda minimises the variance of the new distances */
      let num = 0, den = 0;
      ids.forEach(m => {
        if (m === bi || m === bj) return;
        num += V[bj][m] - V[bi][m];
        den += 1;
      });
      let lambda = 0.5 + (den > 0 ? num / (2 * den * Math.max(V[bi][bj], 1e-12)) : 0);
      lambda = Math.max(0, Math.min(1, lambda));
      const k = next++;
      nodeOf.set(k, { children: [{ node: nodeOf.get(bi), len: li }, { node: nodeOf.get(bj), len: lj }] });
      D[k] = []; V[k] = [];
      ids.forEach(m => {
        if (m === bi || m === bj) return;
        const d = lambda * (D[bi][m] - li) + (1 - lambda) * (D[bj][m] - lj);
        D[k][m] = d; if (!D[m]) D[m] = []; D[m][k] = d;
        const v = lambda * V[bi][m] + (1 - lambda) * V[bj][m] - lambda * (1 - lambda) * V[bi][bj];
        V[k][m] = v; if (!V[m]) V[m] = []; V[m][k] = v;
      });
      D[k][k] = 0; V[k][k] = 0;
      alive.delete(bi); alive.delete(bj); alive.add(k);
    }
    const [a, b] = [...alive];
    return { children: [{ node: nodeOf.get(a), len: Math.max(0, D[a][b] / 2) }, { node: nodeOf.get(b), len: Math.max(0, D[a][b] / 2) }] };
  }

  /* UPGMA and WPGMA: rooted, ultrametric, and only right under a clock */
  function upgma(D0, labels, weighted) {
    const n0 = D0.length;
    const D = D0.map(r => Array.from(r));
    const clusters = [];
    for (let i = 0; i < n0; i++) clusters.push({ id: i, size: 1, height: 0, node: { tip: i, label: labels ? labels[i] : String(i), children: [] } });
    const alive = new Map(clusters.map(c => [c.id, c]));
    let next = n0;
    while (alive.size > 1) {
      const ids = [...alive.keys()];
      let bi = ids[0], bj = ids[1], best = Infinity;
      for (let a = 0; a < ids.length; a++) for (let b = a + 1; b < ids.length; b++) {
        const v = D[ids[a]][ids[b]];
        if (v < best) { best = v; bi = ids[a]; bj = ids[b]; }
      }
      const ci = alive.get(bi), cj = alive.get(bj);
      const h = best / 2;
      const node = { children: [{ node: ci.node, len: Math.max(0, h - ci.height) }, { node: cj.node, len: Math.max(0, h - cj.height) }] };
      const k = next++;
      const merged = { id: k, size: ci.size + cj.size, height: h, node };
      D[k] = [];
      alive.forEach((c, id) => {
        if (id === bi || id === bj) return;
        const d = weighted ? (D[bi][id] + D[bj][id]) / 2
          : (D[bi][id] * ci.size + D[bj][id] * cj.size) / (ci.size + cj.size);
        D[k][id] = d; if (!D[id]) D[id] = []; D[id][k] = d;
      });
      D[k][k] = 0;
      alive.delete(bi); alive.delete(bj); alive.set(k, merged);
    }
    return [...alive.values()][0].node;
  }

  /* Remove an artificial root of degree two, so that the tree is what it really
     is: unrooted, with every internal node of degree three. Distance methods
     produce such a root, and anything that counts neighbours — balanced
     averages, minimum evolution, quartets — must not see it. */
  function unroot(tree) {
    const t = clone(tree);
    if (t.children.length !== 2) return t;
    const a = t.children[0], b = t.children[1];
    const host = a.node.children.length ? a : (b.node.children.length ? b : null);
    if (!host) return t;                       // two tips: nothing to unroot
    const other = host === a ? b : a;
    host.node.children.push({ node: other.node, len: (a.len || 0) + (b.len || 0) });
    return host.node;
  }

  /* ================================================================
     rearrangements: NNI, SPR and TBR
     ================================================================
     Shared by parsimony, minimum evolution and (from Block 5) maximum
     likelihood, so that every criterion searches the same space the same way. */

  /* Every branch whose far end is an internal node — but counted the way an
     unrooted tree counts them. When the root is bifurcating and both its
     children are internal, its two edges are the two halves of ONE unrooted
     branch: they carry the same split, an NNI on either produces the same set
     of topologies (and four of them are no-ops), and a branch test on both
     reports the same clade twice. So only the first of the pair is returned. */
  function internalBranches(tree) {
    const out = [];
    const rootPair = tree.children.length === 2
      && tree.children[0].node.children.length && tree.children[1].node.children.length;
    (function walk(node) {
      node.children.forEach((c, i) => {
        if (c.node.children.length && !(node === tree && rootPair && i === 1)) out.push({ parent: node, edge: c });
        walk(c.node);
      });
    })(tree);
    return out;
  }

  /* The slots on the far side of an internal branch that may be swapped with
     one of its children.

     Normally they are the branch's siblings. The exception is the branch that
     hangs from a bifurcating root: there the root is a node of degree two that
     does not exist in the unrooted tree, so the neighbours on that side are not
     the other root child but ITS children. Swapping the whole other side
     instead would only re-root the tree — a move that changes nothing and
     wastes a likelihood evaluation, which is exactly what it used to do. */
  function nniPartners(tree, br) {
    if (br.parent === tree && tree.children.length === 2) {
      const other = tree.children.find(c => c !== br.edge);
      if (other && other.node.children.length) return other.node.children;
      return [];
    }
    return br.parent.children.filter(c => c !== br.edge);
  }

  /* every partner of an internal branch swapped against every child of it */
  function* nniMoves(tree) {
    const branches = internalBranches(tree);
    for (let b = 0; b < branches.length; b++) {
      const probe = branches[b];
      if (!probe) continue;
      const nSibs = nniPartners(tree, probe).length;
      const nGrand = probe.edge.node.children.length;
      for (let si = 0; si < nSibs; si++) {
        for (let gi = 0; gi < nGrand; gi++) {
          const cand = clone(tree);
          const cb = internalBranches(cand)[b];
          if (!cb) continue;
          const sibs = nniPartners(cand, cb);
          if (si >= sibs.length || gi >= cb.edge.node.children.length) continue;
          const sib = sibs[si], grand = cb.edge.node.children[gi];
          const tmp = sib.node; sib.node = grand.node; grand.node = tmp;
          yield cand;
        }
      }
    }
  }

  function allEdges(tree) {
    const out = [];
    (function walk(node) { node.children.forEach(c => { out.push({ parent: node, edge: c }); walk(c.node); }); })(tree);
    return out;
  }
  function findParentOf(root, target) {
    let found = null;
    (function walk(node) { node.children.forEach(c => { if (c.node === target) found = node; else walk(c.node); }); })(root);
    return found;
  }
  /* every way of rerooting a pruned subtree, capped so that TBR stays finite */
  function tbrRootings(sub, cap) {
    const out = [sub];
    nodes(sub).forEach(nd => {
      if (nd === sub) return;
      try { out.push(rerootAbove(sub, nd, 0.5)); } catch (e) { /* skip */ }
    });
    return out.slice(0, cap || 12);
  }

  /* Subtree pruning and regrafting. With reroot = true it becomes TBR, which
     also tries every rooting of the pruned subtree. */
  function* sprMoves(tree, reroot, maxMoves) {
    const nEdges = allEdges(tree).length;
    let done = 0;
    for (let i = 0; i < nEdges; i++) {
      const base = clone(tree);
      const cut = allEdges(base)[i];
      if (!cut) continue;
      const pruned = cut.edge.node;
      const parent = cut.parent;
      parent.children = parent.children.filter(c => c !== cut.edge);
      let rest = base;
      if (parent.children.length === 1 && parent !== base) {
        const gp = findParentOf(base, parent);
        if (gp) {
          const e = gp.children.find(c => c.node === parent);
          e.node = parent.children[0].node;
          e.len = (e.len || 0) + (parent.children[0].len || 0);
        }
      } else if (base.children.length === 1) {
        rest = base.children[0].node;
      }
      if (tips(rest).length < 2) continue;
      const rootings = reroot && pruned.children.length ? tbrRootings(pruned) : [pruned];
      for (let r = 0; r < rootings.length; r++) {
        const nRest = allEdges(rest).length;
        for (let j = 0; j < nRest; j++) {
          const cand = clone(rest);
          const at = allEdges(cand)[j];
          if (!at) continue;
          const attach = clone(rootings[r]);
          at.edge.node = { children: [{ node: at.edge.node, len: at.edge.len || 1 }, { node: attach, len: cut.edge.len || 1 }] };
          yield cand;
          if (maxMoves && ++done >= maxMoves) return;
        }
      }
    }
  }

  /* ================================================================
     balanced minimum evolution (Desper & Gascuel 2002)
     ================================================================ */
  /* Pauplin's (2000) formula: the length of a tree under balanced averaging is
     Σ d(i,j) / 2^p(i,j), where p is the number of internal nodes between the
     two tips. No branch lengths are needed to score a topology, which is what
     makes the NNI search affordable. */
  function bmeLength(tree0, D) {
    const tree = tree0.children.length === 2 ? unroot(tree0) : tree0;
    const tipsList = tips(tree);
    const n = tipsList.length;
    const idx = index(tree);
    /* number of internal nodes on the path between every pair */
    const depth = new Map();
    (function walk(nd, d) { depth.set(nd, d); nd.children.forEach(c => walk(c.node, d + 1)); })(tree, 0);
    let total = 0;
    for (let a = 0; a < n; a++) {
      for (let b = a + 1; b < n; b++) {
        const p = internalNodesBetween(tipsList[a], tipsList[b], idx, depth);
        total += D[tipsList[a].tip][tipsList[b].tip] / Math.pow(2, p);
      }
    }
    return total;
  }
  function internalNodesBetween(x, y, idx, depth) {
    /* walk both up to their common ancestor, counting internal nodes */
    let a = x, b = y, count = 0;
    while (a !== b) {
      if (depth.get(a) >= depth.get(b)) { a = idx.parent.get(a); count++; }
      else { b = idx.parent.get(b); count++; }
    }
    /* the two tips themselves were counted as steps, so the internal nodes on
       the path are count − 1 when the ancestor is counted once */
    return count - 1;
  }

  /* balanced average distance between the two sides of every edge, used for the
     branch lengths of a minimum-evolution tree */
  function balancedAverages(tree, D) {
    const cache = new Map();
    const key = (a, b) => (a.__id < b.__id ? a.__id + '|' + b.__id : b.__id + '|' + a.__id);
    let id = 0;
    nodes(tree).forEach(n => { n.__id = id++; });
    const idx = index(tree);
    /* the set of tips on the far side of a node, seen from a neighbour */
    function avg(A, fromA, B, fromB) {
      const k = A.__id + '>' + (fromA ? fromA.__id : -1) + '|' + B.__id + '>' + (fromB ? fromB.__id : -1);
      if (cache.has(k)) return cache.get(k);
      let v;
      const kidsA = neighbours(A, fromA, idx);
      const kidsB = neighbours(B, fromB, idx);
      if (!kidsA.length && !kidsB.length) v = D[A.tip][B.tip];
      else if (kidsA.length) v = (avg(kidsA[0], A, B, fromB) + avg(kidsA[1], A, B, fromB)) / 2;
      else v = (avg(A, fromA, kidsB[0], B) + avg(A, fromA, kidsB[1], B)) / 2;
      cache.set(k, v);
      return v;
    }
    return avg;
  }
  /* the neighbours of a node other than `from`, in the unrooted sense */
  function neighbours(node, from, idx) {
    const out = [];
    node.children.forEach(c => { if (c.node !== from) out.push(c.node); });
    const p = idx.parent.get(node);
    if (p && p !== from) out.push(p);
    return out;
  }

  /* Minimum evolution from several starting trees, keeping the shortest.
     Balanced minimum evolution is a heuristic search like any other, and where
     it lands depends on where it starts: on the rbcL example, starting from
     BIONJ gives 2.903640, from UPGMA 2.899076 — shorter than what ape's
     fastme.bal reaches (2.901969) from its own greedy start. Trying the three
     obvious starts costs milliseconds and removes the lottery. */
  function meMulti(D, labels, opts) {
    opts = opts || {};
    const starts = opts.start ? [opts.start] : [bionj(D, labels), nj(D, labels), upgma(D, labels)];
    let best = null;
    starts.forEach(s => {
      const r = me(D, labels, Object.assign({}, opts, { start: s }));
      if (!best || r.length < best.length - 1e-12) best = r;
    });
    best.starts = starts.length;
    return best;
  }

  /* one search, from one starting tree */
  function me(D, labels, opts) {
    opts = opts || {};
    /* unrooted from the start: the balanced averages count neighbours, and an
       artificial root of degree two would be counted as a real node */
    let tree = unroot(opts.start || bionj(D, labels));
    let best = bmeLength(tree, D);
    let rounds = 0;
    let improved = true;
    while (improved && rounds < (opts.maxRounds || 50)) {
      improved = false;
      rounds++;
      /* NNI first, then SPR — which is what ape's fastme.bal does, and the
         reason it reaches a shorter tree than an NNI-only search */
      for (const cand of nniMoves(tree)) {
        const L = bmeLength(cand, D);
        if (L < best - 1e-12) { best = L; tree = cand; improved = true; break; }
      }
      if (!improved && opts.spr !== false) {
        for (const cand of sprMoves(tree, false, opts.maxMoves || 6000)) {
          const L = bmeLength(cand, D);
          if (L < best - 1e-12) { best = L; tree = cand; improved = true; break; }
        }
      }
    }
    setBalancedLengths(tree, D);
    return { tree, length: best, rounds };
  }

  /* branch lengths under balanced minimum evolution (Desper & Gascuel 2002) */
  function setBalancedLengths(tree, D) {
    const idx = index(tree);
    const avg = balancedAverages(tree, D);
    nodes(tree).forEach(nd => {
      nd.children.forEach(c => {
        const child = c.node, parent = nd;
        const childSides = neighbours(child, parent, idx);
        const parentSides = neighbours(parent, child, idx);
        if (!childSides.length && parentSides.length === 2) {
          /* external branch: ½ (Δ(i|B) + Δ(i|C) − Δ(B|C)) */
          const [B, C] = parentSides;
          c.len = Math.max(0, 0.5 * (avg(child, parent, B, parent) + avg(child, parent, C, parent) - avg(B, parent, C, parent)));
        } else if (childSides.length === 2 && parentSides.length === 2) {
          const [A, B] = childSides, [C, E] = parentSides;
          const v = 0.5 * ((avg(A, child, C, parent) + avg(B, child, E, parent)) / 2
            + (avg(A, child, E, parent) + avg(B, child, C, parent)) / 2
            - avg(A, child, B, child) - avg(C, parent, E, parent));
          c.len = Math.max(0, v);
        } else if (!childSides.length && parentSides.length) {
          /* a tip hanging from a node of degree three at the root */
          const others = parentSides;
          let s = 0;
          others.forEach(o => { s += avg(child, parent, o, parent); });
          c.len = Math.max(0, s / others.length);
        }
      });
    });
    nodes(tree).forEach(n => { delete n.__id; });
    return tree;
  }

  /* ================================================================
     flat form for the likelihood engine
     ================================================================ */
  /* Turns the tree into arrays: parent, branch length, post-order, and which
     row of the alignment each tip is. Node 0 is the root. */
  function flatten(root) {
    const list = [], parent = [], len = [], tipRow = [], kids = [];
    (function walk(n, p, l) {
      const k = list.length;
      list.push(n); parent.push(p); len.push(l); kids.push([]);
      tipRow.push(isTip(n) ? n.tip : -1);
      if (p >= 0) kids[p].push(k);
      n.children.forEach(c => walk(c.node, k, c.len || 0));
    })(root, -1, 0);
    const post = [];
    (function order(k) { kids[k].forEach(order); post.push(k); })(0);
    return { n: list.length, nodes: list, parent, len: Float64Array.from(len), kids, tipRow, post, isTip: tipRow.map(r => r >= 0) };
  }
  function unflatten(F, lens) {
    const build = k => {
      const n = { label: F.nodes[k].label, support: F.nodes[k].support, children: [] };
      if (F.isTip[k]) n.tip = F.tipRow[k];
      F.kids[k].forEach(c => n.children.push({ node: build(c), len: lens ? lens[c] : F.len[c] }));
      return n;
    };
    return build(0);
  }

  Object.assign(Tree, {
    isTip, tips, nodes, nTips, index, clone,
    parseNewick, writeNewick, splits, rfDistance, totalLength,
    rerootAbove, rootByOutgroup, midpointRoot, ladderize, collapseShort, collapseSupport,
    nj, bionj, upgma, unroot, me, meMulti, bmeLength, internalBranches, nniMoves, nniPartners, sprMoves, tbrRootings, allEdges, setBalancedLengths, balancedAverages, neighbours, flatten, unflatten,
  });
  g.Tree = Tree;
}
TreeCore(typeof window !== 'undefined' ? window : self);
