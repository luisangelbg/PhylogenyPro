/* PhylogenyPro — molecular phylogenetics in the browser.
   Copyright (C) 2026  Luis Ángel Barrera-Guzmán

   This program is free software: you can redistribute it and/or modify it under
   the terms of the GNU General Public License as published by the Free Software
   Foundation, either version 3 of the License, or (at your option) any later
   version. It is distributed in the hope that it will be useful, but WITHOUT ANY
   WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
   A PARTICULAR PURPOSE. See the GNU General Public License, in the file LICENSE
   at the root of this program, or <https://www.gnu.org/licenses/>. */

/* PhylogenyPro — Block 1: the two interactive labs on the home page.

   Nothing here is an animation. Both labs compute real phylogenetics:

   LAB 1 · The Felsenstein zone
     A four-taxon tree ((A,B),(C,D)) with two long branches that are not
     sisters. Under Jukes–Cantor the probability of every site pattern has a
     closed form, and by symmetry the 256 patterns collapse into the 15 set
     partitions of {A,B,C,D}, so everything is exact and fast:
       · infinite data  — the expected frequency of each informative pattern
         decides whether parsimony is consistent, and the maximum-likelihood
         score is maximised over branch lengths on each of the three
         topologies (this is the Kullback–Leibler argument of Felsenstein 1978);
       · finite data    — sites are drawn from that same distribution, and
         parsimony, neighbour joining (four-point condition on Jukes–Cantor
         distances) and maximum likelihood are asked to recover the tree.
     The boundary of the zone (where the pattern supporting the wrong tree
     becomes as frequent as the one supporting the true tree) is found by
     bisection and drawn in the (q, p) plane.

   LAB 2 · Saturation and distance correction
     Two sequences diverging under K80 with gamma-distributed rates. The
     observed proportion of differences, the Jukes–Cantor and the Kimura
     two-parameter distances are compared with the true number of
     substitutions per site, and the transition/transversion ratio is followed
     until it collapses — which is what saturation looks like in real data. */

(function () {

  const V = n => `var(--${n})`;
  const f1 = v => (+v).toFixed(1);
  function svg(vb, inner) { return `<svg viewBox="${vb}" xmlns="http://www.w3.org/2000/svg">${inner}</svg>`; }
  const line = (x1, y1, x2, y2, st, w, ex) => `<line x1="${f1(x1)}" y1="${f1(y1)}" x2="${f1(x2)}" y2="${f1(y2)}" stroke="${st}" stroke-width="${w || 1}" stroke-linecap="round" ${ex || ''}/>`;
  const circ = (cx, cy, r, fill, ex) => `<circle cx="${f1(cx)}" cy="${f1(cy)}" r="${f1(r)}" fill="${fill}" ${ex || ''}/>`;
  const rect = (x, y, w, h, fill, ex) => `<rect x="${f1(x)}" y="${f1(y)}" width="${f1(Math.max(0, w))}" height="${f1(Math.max(0, h))}" fill="${fill}" ${ex || ''}/>`;
  const pathEl = (d, st, w, ex) => `<path d="${d}" stroke="${st}" stroke-width="${w || 1}" fill="none" stroke-linecap="round" ${ex || ''}/>`;
  const poly = pts => pts.map((p, i) => (i ? 'L' : 'M') + f1(p[0]) + ' ' + f1(p[1])).join(' ');
  function txt(x, y, s, cls, size, anchor, ex) {
    return `<text x="${f1(x)}" y="${f1(y)}" class="${cls || 'art-mut'}" font-size="${size || 8}" text-anchor="${anchor || 'start'}" ${ex || ''}>${s}</text>`;
  }

  /* ================================================================
     1 · THE 15 PATTERN CLASSES OF FOUR TAXA
     Under any model with equal base frequencies and a single exchangeability
     (Jukes–Cantor), the probability of a site pattern depends only on which
     taxa share a state, not on which states. The set partitions of
     {A, B, C, D} are therefore the natural units: 15 of them, with the
     multiplicity of the concrete patterns each one stands for.
     ================================================================ */
  const CLASSES = (function () {
    const out = [];
    /* every assignment of the four taxa to blocks, in canonical (restricted
       growth) form: g[0] = 0 and g[i] ≤ max(g[0..i-1]) + 1 */
    const g = [0, 0, 0, 0];
    const rec = i => {
      if (i === 4) {
        const k = Math.max(...g) + 1;
        /* multiplicity: ordered choices of k distinct bases out of 4 */
        let mult = 1;
        for (let j = 0; j < k; j++) mult *= (4 - j);
        out.push({ g: g.slice(), k, mult });
        return;
      }
      const lim = Math.max(...g.slice(0, i)) + 1;
      for (let v = 0; v <= lim; v++) { g[i] = v; rec(i + 1); }
    };
    rec(1);
    return out.map(c => {
      /* which pairs differ, for the pairwise distances */
      const diff = {};
      for (let i = 0; i < 4; i++) for (let j = i + 1; j < 4; j++) diff[i + '' + j] = c.g[i] !== c.g[j] ? 1 : 0;
      /* informative patterns: exactly two blocks of two */
      const sizes = {};
      c.g.forEach(b => sizes[b] = (sizes[b] || 0) + 1);
      const vals = Object.values(sizes).sort();
      let informs = null;
      if (c.k === 2 && vals[0] === 2 && vals[1] === 2) {
        if (c.g[0] === c.g[1]) informs = 0;          // AB | CD  → topology 1
        else if (c.g[0] === c.g[2]) informs = 1;     // AC | BD  → topology 2
        else informs = 2;                            // AD | BC  → topology 3
      }
      return Object.assign({}, c, { diff, informs });
    });
  })();

  /* Jukes–Cantor transition probability of staying (same) or changing (diff) */
  function jcP(t) {
    const e = Math.exp(-4 * t / 3);
    return { same: 0.25 + 0.75 * e, diff: 0.25 - 0.25 * e };
  }
  /* Probability of each of the 15 classes on the topology given by `pairing`
     (0 = ((A,B),(C,D)), 1 = ((A,C),(B,D)), 2 = ((A,D),(B,C))) with branch
     lengths [tA, tB, tC, tD, tI] written in the order of that pairing.
     Sum over the two internal nodes: 4 × 4 = 16 terms. */
  const PAIRINGS = [[0, 1, 2, 3], [0, 2, 1, 3], [0, 3, 1, 2]];
  function classProbs(pairing, lens) {
    const ord = PAIRINGS[pairing];
    const P = lens.map(jcP);
    const probs = new Array(CLASSES.length).fill(0);
    CLASSES.forEach((cl, ci) => {
      /* one representative concrete pattern of the class: block index = state */
      const st = [0, 0, 0, 0];
      for (let i = 0; i < 4; i++) st[i] = cl.g[i];
      let p = 0;
      for (let x = 0; x < 4; x++) for (let y = 0; y < 4; y++) {
        /* x is the internal node joining ord[0], ord[1]; y joins ord[2], ord[3] */
        let term = 0.25;
        term *= (st[ord[0]] === x ? P[0].same : P[0].diff);
        term *= (st[ord[1]] === x ? P[1].same : P[1].diff);
        term *= (x === y ? P[4].same : P[4].diff);
        term *= (st[ord[2]] === y ? P[2].same : P[2].diff);
        term *= (st[ord[3]] === y ? P[3].same : P[3].diff);
        p += term;
      }
      probs[ci] = p * cl.mult;
    });
    return probs;
  }

  /* log-likelihood of observed class counts under a topology and branch lengths */
  function lnLik(counts, pairing, lens) {
    const p = classProbs(pairing, lens);
    let s = 0;
    for (let i = 0; i < counts.length; i++) if (counts[i] > 0) s += counts[i] * Math.log(Math.max(p[i] / CLASSES[i].mult, 1e-300));
    return s;
  }
  /* golden-section search of the maximum along one branch */
  function optimiseBranch(f, lo, hi) {
    const gr = (Math.sqrt(5) - 1) / 2;
    let a = lo, b = hi, c = b - gr * (b - a), d = a + gr * (b - a);
    let fc = f(c), fd = f(d);
    for (let it = 0; it < 60 && b - a > 1e-9; it++) {
      if (fc > fd) { b = d; d = c; fd = fc; c = b - gr * (b - a); fc = f(c); }
      else { a = c; c = d; fc = fd; d = a + gr * (b - a); fd = f(d); }
    }
    return (a + b) / 2;
  }
  /* maximum likelihood on a fixed topology: cyclic optimisation of the five
     branches (the tree is unrooted, so five free lengths is exactly right).
     The lower bound has to be a true zero, not a small positive number: on a
     wrong topology the internal branch collapses, and holding it at 1e-5 costs
     about 0.001 log-likelihood units — enough to disagree with phangorn, as the
     validation of this block showed. */
  function fitTopology(counts, pairing, passes) {
    const lens = [0.1, 0.1, 0.1, 0.1, 0.05];
    let best = lnLik(counts, pairing, lens);
    for (let pass = 0; pass < (passes || 4); pass++) {
      const before = best;
      for (let k = 0; k < 5; k++) {
        const x = optimiseBranch(t => { const L = lens.slice(); L[k] = t; return lnLik(counts, pairing, L); }, 1e-9, 5);
        lens[k] = x;
        best = lnLik(counts, pairing, lens);
      }
      if (Math.abs(best - before) < 1e-9) break;
    }
    return { lnL: best, lens };
  }

  /* ================================================================
     2 · LAB 1 — THE FELSENSTEIN ZONE
     ================================================================ */
  const F = {
    p: 0.55,      // long branches (A and C)
    q: 0.05,      // short branches (B, D) and the internal branch
    L: 500,       // sites
    reps: 100,    // replicates
    seed: 1,
  };

  /* true tree: topology 1, lengths [A=p, B=q, C=p, D=q, internal=q] */
  function trueLens() { return [F.p, F.q, F.p, F.q, F.q]; }

  /* expected frequencies of the three informative patterns */
  function expectedInformative(p, q) {
    const probs = classProbs(0, [p, q, p, q, q]);
    const f = [0, 0, 0];
    CLASSES.forEach((c, i) => { if (c.informs != null) f[c.informs] += probs[i]; });
    return f;
  }
  /* boundary of the zone: for a given q, the value of p at which the pattern
     supporting the wrong topology becomes as frequent as the true one */
  function zoneBoundary(q) {
    const g = p => { const f = expectedInformative(p, q); return f[1] - f[0]; };
    let lo = q, hi = 3;
    if (g(hi) < 0) return null;               // parsimony stays consistent
    if (g(lo) > 0) return lo;
    for (let i = 0; i < 40; i++) { const m = (lo + hi) / 2; if (g(m) > 0) hi = m; else lo = m; }
    return (lo + hi) / 2;
  }

  /* draw L sites from the class distribution of the true tree */
  function simulate(probs, L, r) {
    const cum = [];
    let acc = 0;
    probs.forEach(p => { acc += p; cum.push(acc); });
    const counts = new Array(probs.length).fill(0);
    for (let s = 0; s < L; s++) {
      const u = r() * acc;
      let lo = 0, hi = cum.length - 1;
      while (lo < hi) { const m = (lo + hi) >> 1; if (u <= cum[m]) hi = m; else lo = m + 1; }
      counts[lo]++;
    }
    return counts;
  }
  /* pairwise Jukes–Cantor distances from class counts */
  function pairDistances(counts, L) {
    const d = {};
    ['01', '02', '03', '12', '13', '23'].forEach(key => {
      let diff = 0;
      CLASSES.forEach((c, i) => { if (c.diff[key]) diff += counts[i]; });
      const p = diff / L;
      const jc = jcDistance(p);
      d[key] = jc == null ? 3 : jc;      // saturated pairs are capped, as in any program
    });
    return d;
  }
  /* which topology the four-point condition picks (this is what neighbour
     joining does with four taxa) */
  function njChoice(d) {
    const s = [d['01'] + d['23'], d['02'] + d['13'], d['03'] + d['12']];
    const min = Math.min(...s);
    const ties = s.filter(v => v === min).length;
    return ties > 1 ? -1 : s.indexOf(min);
  }
  function mpChoice(counts) {
    const f = [0, 0, 0];
    CLASSES.forEach((c, i) => { if (c.informs != null) f[c.informs] += counts[i]; });
    const max = Math.max(...f);
    return f.filter(v => v === max).length > 1 ? -1 : f.indexOf(max);
  }
  function mlChoice(counts) {
    const fits = [0, 1, 2].map(t => fitTopology(counts, t, 3));
    const best = Math.max(...fits.map(f2 => f2.lnL));
    const win = fits.findIndex(f2 => f2.lnL === best);
    return { win, fits };
  }

  let felsenResult = null;
  function runFelsen() {
    const r = rng(F.seed);
    const probs = classProbs(0, trueLens());
    const res = { mp: 0, nj: 0, ml: 0, mpTie: 0, njTie: 0, reps: F.reps };
    for (let rep = 0; rep < F.reps; rep++) {
      const counts = simulate(probs, F.L, r);
      const mp = mpChoice(counts);
      const nj = njChoice(pairDistances(counts, F.L));
      const ml = mlChoice(counts).win;
      if (mp === 0) res.mp++; else if (mp === -1) res.mpTie++;
      if (nj === 0) res.nj++; else if (nj === -1) res.njTie++;
      if (ml === 0) res.ml++;
    }
    /* infinite data: the expected pattern frequencies decide */
    const fInf = expectedInformative(F.p, F.q);
    const bigCounts = probs.map(p => p * 1e6);
    const mlInf = mlChoice(bigCounts);
    res.expected = fInf;
    res.mpConsistent = fInf[0] > fInf[1];
    res.mlInfinite = mlInf.win === 0;
    res.mlInfiniteGap = mlInf.fits[0].lnL - Math.max(mlInf.fits[1].lnL, mlInf.fits[2].lnL);
    felsenResult = res;
    return res;
  }

  /* ---------------- drawings of lab 1 ---------------- */
  function drawFourTaxon(host) {
    const W = 300, H = 220;
    const L = (es, en) => T(es, en);
    const scale = 150;                       // pixels per substitution/site, clipped
    const px = t => Math.min(120, t * scale);
    /* the internal branch is drawn horizontally between the two cherries */
    const mid = 150, half = Math.max(6, px(F.q)) / 2;
    const xL = mid - half, xR = mid + half;
    let s = '';
    s += line(xL, 110, xR, 110, V('text-muted'), 3);
    s += line(xL, 70, xL, 150, V('text-muted'), 3);
    s += line(xR, 70, xR, 150, V('text-muted'), 3);
    /* A long, B short (left cherry); C long, D short (right cherry) */
    const A = xL - px(F.p), B = xL - px(F.q), C = xR + px(F.p), D = xR + px(F.q);
    s += line(xL, 70, A, 70, V('danger'), 3.4);
    s += line(xL, 150, B, 150, V('primary'), 3.4);
    s += line(xR, 70, C, 70, V('danger'), 3.4);
    s += line(xR, 150, D, 150, V('primary'), 3.4);
    s += circ(A, 70, 4, V('danger')) + circ(B, 150, 4, V('primary'));
    s += circ(C, 70, 4, V('danger')) + circ(D, 150, 4, V('primary'));
    s += txt(A - 7, 73, 'A', 'art-txt', 11, 'end', `fill="${V('danger')}" font-weight="700"`);
    s += txt(B - 7, 153, 'B', 'art-txt', 11, 'end', `fill="${V('primary')}" font-weight="700"`);
    s += txt(C + 7, 73, 'C', 'art-txt', 11, 'start', `fill="${V('danger')}" font-weight="700"`);
    s += txt(D + 7, 153, 'D', 'art-txt', 11, 'start', `fill="${V('primary')}" font-weight="700"`);
    s += txt(150, 26, L('Árbol verdadero: A y C tienen ramas largas y NO son hermanas',
      'True tree: A and C have long branches and are NOT sisters'), 'art-mut', 8.4, 'middle');
    /* each branch labelled beside itself, never on top of the internal branch */
    s += txt((A + xL) / 2, 62, 'p = ' + F.p.toFixed(2), 'art-mut', 8, 'middle', `fill="${V('danger')}"`);
    s += txt((C + xR) / 2, 62, 'p = ' + F.p.toFixed(2), 'art-mut', 8, 'middle', `fill="${V('danger')}"`);
    s += txt(B - 4, 166, 'q = ' + F.q.toFixed(3), 'art-mut', 8, 'end', `fill="${V('primary')}"`);
    s += txt(D + 4, 166, 'q = ' + F.q.toFixed(3), 'art-mut', 8, 'start', `fill="${V('primary')}"`);
    s += txt(mid, 126, L('rama interna = q', 'internal branch = q'), 'art-mut', 7.6, 'middle');
    s += line(20, 196, 20 + px(0.1), 196, V('text-muted'), 2);
    s += txt(20, 208, L('0.1 sustituciones/sitio', '0.1 substitutions/site'), 'art-mut', 7.6);
    host.innerHTML = svg(`0 0 ${W} ${H}`, s);
  }

  function drawZone(host) {
    const W = 300, H = 220, mx = 42, my = 34;
    const L = (es, en) => T(es, en);
    const qMax = 0.3, pMax = 1.2;
    const X = q => mx + (W - mx - 14) * q / qMax;
    const Y = p => (H - my) - (H - my - 16) * p / pMax;
    let s = '';
    /* the zone: above the boundary parsimony is inconsistent */
    const pts = [];
    for (let q = 0.002; q <= qMax + 1e-9; q += qMax / 60) {
      const b = zoneBoundary(q);
      pts.push([X(q), Y(b == null ? pMax * 1.2 : Math.min(pMax * 1.2, b))]);
    }
    s += `<path d="${poly(pts)} L${f1(X(qMax))} ${f1(Y(pMax * 1.2))} L${f1(X(0.002))} ${f1(Y(pMax * 1.2))} Z" fill="${V('danger')}" opacity="0.14"/>`;
    s += pathEl(poly(pts), V('danger'), 2);
    s += line(mx, H - my, W - 14, H - my, V('border-strong'), 1.4);
    s += line(mx, 12, mx, H - my, V('border-strong'), 1.4);
    for (let i = 0; i <= 3; i++) {
      const q = qMax * i / 3;
      s += line(X(q), H - my, X(q), H - my + 4, V('border-strong'), 1.2);
      s += txt(X(q), H - my + 14, q.toFixed(2), 'art-mut', 7.6, 'middle');
    }
    for (let i = 0; i <= 4; i++) {
      const p = pMax * i / 4;
      s += line(mx - 4, Y(p), mx, Y(p), V('border-strong'), 1.2);
      s += txt(mx - 6, Y(p) + 3, p.toFixed(1), 'art-mut', 7.6, 'end');
    }
    s += txt(W / 2, H - 6, L('ramas cortas  q', 'short branches  q'), 'art-mut', 8.4, 'middle');
    s += txt(-H / 2, 11, L('ramas largas  p', 'long branches  p'), 'art-mut', 8.4, 'middle', 'transform="rotate(-90)"');
    s += txt(X(qMax * 0.52), Y(pMax * 0.86), L('la parsimonia falla aquí', 'parsimony fails here'), 'art-mut', 8.4, 'middle', `fill="${V('danger')}"`);
    s += txt(X(qMax * 0.5), Y(pMax * 0.18), L('zona segura', 'safe zone'), 'art-mut', 8.4, 'middle', `fill="${V('leaf')}"`);
    /* the current point */
    const cx = X(Math.min(F.q, qMax)), cy = Y(Math.min(F.p, pMax));
    s += circ(cx, cy, 6.5, V('card-bg'), `stroke="${V('text')}" stroke-width="2"`);
    s += circ(cx, cy, 3, V('accent'));
    host.innerHTML = svg(`0 0 ${W} ${H}`, s);
  }

  function drawPatterns(host) {
    const W = 300, H = 220;
    const L = (es, en) => T(es, en);
    const f = felsenResult ? felsenResult.expected : expectedInformative(F.p, F.q);
    const tot = f[0] + f[1] + f[2] || 1;
    const names = ['AB|CD', 'AC|BD', 'AD|BC'];
    const sub = [L('apoya el árbol verdadero', 'supports the true tree'),
                 L('junta las dos ramas largas', 'joins the two long branches'),
                 L('la tercera topología', 'the third topology')];
    const cols = [V('leaf'), V('danger'), V('text-muted')];
    let s = txt(W / 2, 16, L('Frecuencia esperada de los sitios informativos', 'Expected frequency of the informative sites'), 'art-mut', 8.6, 'middle');
    const maxv = Math.max(...f) || 1;
    f.forEach((v, i) => {
      const y = 40 + i * 52;
      s += txt(14, y, names[i], 'art-txt', 10, 'start', `fill="${cols[i]}" font-weight="700"`);
      s += txt(14, y + 12, sub[i], 'art-mut', 7.8);
      s += rect(14, y + 18, W - 28, 12, V('bg-soft'), 'rx="6"');
      s += rect(14, y + 18, (W - 28) * v / maxv, 12, cols[i], 'rx="6" opacity="0.9"');
      s += txt(W - 14, y + 2, (v * 100).toFixed(2) + '%', 'art-mut', 9, 'end');
    });
    const win = f[1] > f[0];
    s += txt(W / 2, 205, win
      ? L('El patrón equivocado es el más frecuente: con datos infinitos, la parsimonia se equivoca.',
          'The wrong pattern is the most frequent: with infinite data, parsimony gets it wrong.')
      : L('El patrón verdadero es el más frecuente: la parsimonia es consistente aquí.',
          'The true pattern is the most frequent: parsimony is consistent here.'),
      'art-mut', 7.8, 'middle', `fill="${win ? V('danger') : V('leaf')}"`);
    host.innerHTML = svg(`0 0 ${W} ${H}`, s);
  }

  function renderFelsenBars() {
    const host = el('fzBars');
    if (!host || !felsenResult) return;
    const r = felsenResult;
    const rows = [
      ['ml', T('Máxima verosimilitud', 'Maximum likelihood'), r.ml / r.reps, 'leaf'],
      ['nj', T('Distancias + NJ', 'Distances + NJ'), r.nj / r.reps, 'primary'],
      ['mp', T('Máxima parsimonia', 'Maximum parsimony'), r.mp / r.reps, 'accent'],
    ];
    const best = Math.max(...rows.map(x => x[2]));
    host.innerHTML = rows.map(x => `<div class="win-row${x[2] === best ? ' best' : ''}">
        <span class="wl">${x[1]}</span>
        <span class="wb"><i style="width:${(x[2] * 100).toFixed(1)}%;background:var(--${x[3]})"></i></span>
        <span class="wv">${(x[2] * 100).toFixed(0)}%</span></div>`).join('');
  }

  function felsenStatus() {
    const host = el('fzStatus');
    if (!host || !felsenResult) return;
    const r = felsenResult;
    const inZone = !r.mpConsistent;
    const parts = [];
    parts.push(T(
      `Con <b>${r.reps}</b> réplicas de <b>${F.L}</b> sitios, la verosimilitud recuperó el árbol verdadero en <b>${(r.ml / r.reps * 100).toFixed(0)}%</b> de los casos, el vecino más cercano en <b>${(r.nj / r.reps * 100).toFixed(0)}%</b> y la parsimonia en <b>${(r.mp / r.reps * 100).toFixed(0)}%</b>.`,
      `With <b>${r.reps}</b> replicates of <b>${F.L}</b> sites, likelihood recovered the true tree in <b>${(r.ml / r.reps * 100).toFixed(0)}%</b> of the runs, neighbour joining in <b>${(r.nj / r.reps * 100).toFixed(0)}%</b> and parsimony in <b>${(r.mp / r.reps * 100).toFixed(0)}%</b>.`));
    if (inZone) {
      parts.push(T(
        'Estás <b>dentro de la zona de Felsenstein</b>: por más sitios que añadas, la parsimonia converge con seguridad creciente hacia el árbol equivocado. Eso es inconsistencia estadística, no mala suerte.',
        'You are <b>inside the Felsenstein zone</b>: however many sites you add, parsimony converges with growing confidence on the wrong tree. That is statistical inconsistency, not bad luck.'));
      parts.push(r.mlInfinite
        ? T('Con datos infinitos la verosimilitud sigue prefiriendo el árbol verdadero: por eso se dice que es consistente bajo el modelo correcto.',
            'With infinite data likelihood still prefers the true tree: that is what being consistent under the correct model means.')
        : T('Ni siquiera la verosimilitud lo resuelve con estas longitudes: las ramas largas ya borraron la señal.',
            'Not even likelihood resolves it with these lengths: the long branches have erased the signal.'));
    } else {
      parts.push(T('Estás <b>fuera de la zona</b>: aquí los tres métodos son consistentes y la diferencia entre ellos es solo de eficiencia (cuántos sitios necesitan).',
        'You are <b>outside the zone</b>: here all three methods are consistent and the difference between them is only efficiency (how many sites they need).'));
    }
    host.innerHTML = parts.join(' ');
  }

  function updateFelsen(run) {
    if (run !== false) runFelsen();
    const t = el('fzTree'), z = el('fzZone'), p = el('fzPatterns');
    if (t) drawFourTaxon(t);
    if (z) drawZone(z);
    if (p) drawPatterns(p);
    renderFelsenBars();
    felsenStatus();
  }

  /* ================================================================
     3 · LAB 2 — SATURATION AND DISTANCE CORRECTION
     ================================================================ */
  const S = { kappa: 4, alpha: 0.5, gamma: true, L: 600, seed: 2, tMax: 1.5 };

  /* K80 transition probabilities for a branch of length t and rate ratio kappa.
     Returns the probability of a transition (P) and of a transversion (Q). */
  function k80Probs(t, kappa) {
    const b = 1 / (kappa + 2);              // rate of each transversion type
    const a = kappa * b;                    // rate of the transition
    const e1 = Math.exp(-4 * b * t);
    const e2 = Math.exp(-2 * (a + b) * t);
    const P = 0.25 + 0.25 * e1 - 0.5 * e2;  // transition
    const Q = 0.5 - 0.5 * e1;               // both transversions together
    return { P, Q };
  }
  /* averaged over gamma-distributed rates, by a four-category discrete
     approximation (Yang 1994) using the MEDIAN of each equal-probability
     category rescaled to mean 1. The analysis blocks will offer the mean
     method as well; for this lab the medians are enough and cheaper. */
  function gammaRates(alpha, k) {
    if (!alpha || alpha <= 0) return [1];
    const out = [];
    /* medians of the k equal-probability categories, rescaled to mean 1 */
    for (let i = 0; i < k; i++) {
      const p = (i + 0.5) / k;
      out.push(invGamma(p, alpha));
    }
    const m = out.reduce((a, b) => a + b, 0) / k;
    return out.map(v => v / m);
  }
  /* quantile of Gamma(alpha, alpha) by bisection on a series for the CDF */
  function invGamma(p, a) {
    const cdf = x => {
      if (x <= 0) return 0;
      /* regularised lower incomplete gamma P(a, a·x) by series/continued fraction */
      const z = a * x;
      if (z < a + 1) {
        let sum = 1 / a, del = sum, ap = a;
        for (let n = 0; n < 300; n++) { ap++; del *= z / ap; sum += del; if (Math.abs(del) < Math.abs(sum) * 1e-13) break; }
        return sum * Math.exp(-z + a * Math.log(z) - lgamma(a));
      }
      let b = z + 1 - a, c = 1e300, d = 1 / b, h = d;
      for (let i = 1; i < 300; i++) {
        const an = -i * (i - a); b += 2;
        d = an * d + b; if (Math.abs(d) < 1e-300) d = 1e-300;
        c = b + an / c; if (Math.abs(c) < 1e-300) c = 1e-300;
        d = 1 / d; const del = d * c; h *= del;
        if (Math.abs(del - 1) < 1e-13) break;
      }
      return 1 - Math.exp(-z + a * Math.log(z) - lgamma(a)) * h;
    };
    let lo = 1e-6, hi = 50;
    for (let i = 0; i < 200; i++) { const m = (lo + hi) / 2; if (cdf(m) < p) lo = m; else hi = m; }
    return (lo + hi) / 2;
  }
  function lgamma(x) {
    const c = [76.18009172947146, -86.50532032941677, 24.01409824083091, -1.231739572450155, 0.1208650973866179e-2, -0.5395239384953e-5];
    let y = x, tmp = x + 5.5; tmp -= (x + 0.5) * Math.log(tmp);
    let ser = 1.000000000190015;
    for (let j = 0; j < 6; j++) ser += c[j] / ++y;
    return -tmp + Math.log(2.5066282746310005 * ser / x);
  }

  /* expected proportions of transitions and transversions at true distance t */
  function expectedPQ(t) {
    const rates = S.gamma ? gammaRates(S.alpha, 4) : [1];
    let P = 0, Q = 0;
    rates.forEach(r => { const v = k80Probs(t * r, S.kappa); P += v.P / rates.length; Q += v.Q / rates.length; });
    return { P, Q };
  }

  let satPoints = null;
  function runSaturation() {
    const r = rng(S.seed);
    const pts = [];
    const nT = 26;
    for (let i = 1; i <= nT; i++) {
      const t = S.tMax * i / nT;
      const e = expectedPQ(t);
      /* one realisation with L sites: each site is a transition, a transversion
         or unchanged, drawn from the expected probabilities */
      let nP = 0, nQ = 0;
      for (let s = 0; s < S.L; s++) {
        const u = r();
        if (u < e.P) nP++; else if (u < e.P + e.Q) nQ++;
      }
      const P = nP / S.L, Q = nQ / S.L;
      pts.push({ t, P, Q, p: P + Q, jc: jcDistance(P + Q), k2p: k2pDistance(P, Q), eP: e.P, eQ: e.Q });
    }
    satPoints = pts;
    return pts;
  }

  function drawSaturation(host) {
    const W = 320, H = 230, mx = 44, my = 34;
    const L = (es, en) => T(es, en);
    const dMax = Math.max(S.tMax, 1) * 1.05;
    const X = t => mx + (W - mx - 12) * t / dMax;
    const Y = d => (H - my) - (H - my - 16) * Math.min(d, dMax) / dMax;
    let s = '';
    s += line(mx, H - my, W - 12, H - my, V('border-strong'), 1.4);
    s += line(mx, 12, mx, H - my, V('border-strong'), 1.4);
    /* the identity line: a perfect estimator */
    s += line(X(0), Y(0), X(dMax), Y(dMax), V('text-muted'), 1.2, 'stroke-dasharray="4 3"');
    const series = [
      ['p', V('rose'), p => p.p],
      ['jc', V('primary'), p => p.jc],
      ['k2p', V('leaf'), p => p.k2p],
    ];
    series.forEach(([, col, get]) => {
      const pts = satPoints.filter(p => get(p) != null).map(p => [X(p.t), Y(get(p))]);
      if (pts.length > 1) s += pathEl(poly(pts), col, 2);
      pts.forEach(p => { s += circ(p[0], p[1], 1.8, col); });
    });
    /* points where the correction breaks down */
    const broke = satPoints.find(p => p.k2p == null || p.jc == null);
    if (broke) {
      s += line(X(broke.t), 12, X(broke.t), H - my, V('danger'), 1.2, 'stroke-dasharray="3 3"');
      s += txt(X(broke.t) + 3, 22, L('la corrección se rompe', 'the correction breaks'), 'art-mut', 7.6, 'start', `fill="${V('danger')}"`);
    }
    for (let i = 0; i <= 3; i++) {
      s += line(X(dMax * i / 3), H - my, X(dMax * i / 3), H - my + 4, V('border-strong'), 1.2);
      s += txt(X(dMax * i / 3), H - my + 14, (dMax * i / 3).toFixed(1), 'art-mut', 7.6, 'middle');
      s += line(mx - 4, Y(dMax * i / 3), mx, Y(dMax * i / 3), V('border-strong'), 1.2);
      s += txt(mx - 6, Y(dMax * i / 3) + 3, (dMax * i / 3).toFixed(1), 'art-mut', 7.6, 'end');
    }
    s += txt(W / 2, H - 6, L('distancia real (sustituciones por sitio)', 'true distance (substitutions per site)'), 'art-mut', 8.4, 'middle');
    s += txt(-H / 2, 11, L('distancia estimada', 'estimated distance'), 'art-mut', 8.4, 'middle', 'transform="rotate(-90)"');
    host.innerHTML = svg(`0 0 ${W} ${H}`, s);
  }

  function drawTsTv(host) {
    const W = 320, H = 230, mx = 40, my = 34;
    const L = (es, en) => T(es, en);
    const dMax = Math.max(S.tMax, 1) * 1.05;
    const X = t => mx + (W - mx - 12) * t / dMax;
    const yMax = 0.6;
    const Y = v => (H - my) - (H - my - 16) * Math.min(v, yMax) / yMax;
    let s = '';
    s += line(mx, H - my, W - 12, H - my, V('border-strong'), 1.4);
    s += line(mx, 12, mx, H - my, V('border-strong'), 1.4);
    const ts = satPoints.map(p => [X(p.t), Y(p.eP)]);
    const tv = satPoints.map(p => [X(p.t), Y(p.eQ)]);
    s += pathEl(poly(ts), V('accent'), 2.2);
    s += pathEl(poly(tv), V('sky'), 2.2);
    satPoints.forEach(p => { s += circ(X(p.t), Y(p.P), 1.7, V('accent')); s += circ(X(p.t), Y(p.Q), 1.7, V('sky')); });
    /* the peak of the transition curve is where saturation starts to bite */
    let peak = satPoints[0];
    satPoints.forEach(p => { if (p.eP > peak.eP) peak = p; });
    if (peak.t < S.tMax * 0.95) {
      s += line(X(peak.t), Y(peak.eP) - 6, X(peak.t), H - my, V('accent'), 1, 'stroke-dasharray="3 3"');
      s += txt(X(peak.t) + 3, Y(peak.eP) - 9, L('aquí empieza la saturación', 'saturation starts here'), 'art-mut', 7.6, 'start', `fill="${V('accent')}"`);
    }
    for (let i = 0; i <= 3; i++) {
      s += line(X(dMax * i / 3), H - my, X(dMax * i / 3), H - my + 4, V('border-strong'), 1.2);
      s += txt(X(dMax * i / 3), H - my + 14, (dMax * i / 3).toFixed(1), 'art-mut', 7.6, 'middle');
    }
    for (let i = 0; i <= 3; i++) {
      s += line(mx - 4, Y(yMax * i / 3), mx, Y(yMax * i / 3), V('border-strong'), 1.2);
      s += txt(mx - 6, Y(yMax * i / 3) + 3, (yMax * i / 3).toFixed(1), 'art-mut', 7.6, 'end');
    }
    s += txt(W / 2, H - 6, L('distancia real', 'true distance'), 'art-mut', 8.4, 'middle');
    s += txt(mx + 8, 20, L('— transiciones', '— transitions'), 'art-mut', 8, 'start', `fill="${V('accent')}"`);
    s += txt(mx + 8, 31, L('— transversiones', '— transversions'), 'art-mut', 8, 'start', `fill="${V('sky')}"`);
    host.innerHTML = svg(`0 0 ${W} ${H}`, s);
  }

  function satStatus() {
    const host = el('satStatus');
    if (!host || !satPoints) return;
    const last = satPoints[satPoints.length - 1];
    const broke = satPoints.find(p => p.k2p == null);
    const bias = last.k2p != null ? (1 - last.p / last.t) : null;
    const parts = [];
    parts.push(T(
      `A una distancia real de <b>${last.t.toFixed(2)}</b> sustituciones por sitio, solo se ven <b>${(last.p * 100).toFixed(1)}%</b> de posiciones diferentes: ${bias != null ? `se pierde el <b>${(bias * 100).toFixed(0)}%</b> de los cambios` : 'la mayoría de los cambios ya no se ven'}, porque un mismo sitio cambió varias veces.`,
      `At a true distance of <b>${last.t.toFixed(2)}</b> substitutions per site, only <b>${(last.p * 100).toFixed(1)}%</b> of the positions look different: ${bias != null ? `<b>${(bias * 100).toFixed(0)}%</b> of the changes are lost` : 'most of the changes can no longer be seen'}, because the same site changed more than once.`));
    parts.push(T(
      'La corrección de Jukes y Cantor y la de Kimura devuelven la cuenta de cambios ocultos; por eso una filogenia no se hace con el porcentaje de diferencias crudo.',
      'The Jukes–Cantor and Kimura corrections give back the hidden changes; that is why a phylogeny is not built on the raw percentage of differences.'));
    if (broke) parts.push(T(
      `Pasando <b>${broke.t.toFixed(2)}</b> sustituciones por sitio la corrección ya no tiene solución: las secuencias se parecen tanto como dos al azar. Ahí ese marcador dejó de servir y hay que buscar uno más lento —del cloroplasto, por ejemplo— o usar solo las transversiones.`,
      `Past <b>${broke.t.toFixed(2)}</b> substitutions per site the correction has no solution: the sequences are as similar as two random ones. That marker has stopped working, and a slower one — from the chloroplast, say — or transversions only, is needed.`));
    host.innerHTML = parts.join(' ');
  }

  function updateSaturation(run) {
    if (run !== false) runSaturation();
    const a = el('satPlot'), b = el('satTsTv');
    if (a) drawSaturation(a);
    if (b) drawTsTv(b);
    satStatus();
  }

  /* ================================================================
     4 · WIRING
     ================================================================ */
  function bindRange(id, valId, set, fmt, after) {
    const inp = el(id), out = el(valId);
    if (!inp) return;
    const show = () => { if (out) out.textContent = fmt ? fmt(+inp.value) : inp.value; };
    show();
    inp.addEventListener('input', () => { show(); set(+inp.value); if (after) after(); });
  }

  function init() {
    if (!el('fzTree')) return;

    /* lab 1 */
    bindRange('fzP', 'fzPVal', v => { F.p = v; }, v => v.toFixed(2), () => updateFelsen());
    bindRange('fzQ', 'fzQVal', v => { F.q = v; }, v => v.toFixed(3), () => updateFelsen());
    bindRange('fzL', 'fzLVal', v => { F.L = v; }, v => String(v), () => updateFelsen());
    bindRange('fzReps', 'fzRepsVal', v => { F.reps = v; }, v => String(v), () => updateFelsen());
    const preset = el('fzPreset');
    if (preset) preset.addEventListener('change', () => {
      const v = preset.value;
      const set = (p, q) => {
        F.p = p; F.q = q;
        const a = el('fzP'), b = el('fzQ');
        if (a) { a.value = p; el('fzPVal').textContent = p.toFixed(2); }
        if (b) { b.value = q; el('fzQVal').textContent = q.toFixed(3); }
      };
      if (v === 'zone') set(0.65, 0.04);
      else if (v === 'safe') set(0.12, 0.09);
      else if (v === 'edge') set(0.38, 0.05);
      else if (v === 'hopeless') set(1.1, 0.02);
      updateFelsen();
    });
    const again = el('fzAgain');
    if (again) again.addEventListener('click', () => { F.seed = (F.seed + 1) % 100000; updateFelsen(); });

    /* lab 2 */
    bindRange('satKappa', 'satKappaVal', v => { S.kappa = v; }, v => v.toFixed(1), () => updateSaturation());
    bindRange('satAlpha', 'satAlphaVal', v => { S.alpha = v; }, v => v.toFixed(2), () => updateSaturation());
    bindRange('satT', 'satTVal', v => { S.tMax = v; }, v => v.toFixed(1), () => updateSaturation());
    bindRange('satL', 'satLVal', v => { S.L = v; }, v => String(v), () => updateSaturation());
    const gam = el('satGamma');
    if (gam) gam.addEventListener('change', () => { S.gamma = gam.checked; updateSaturation(); });
    const satAgain = el('satAgain');
    if (satAgain) satAgain.addEventListener('click', () => { S.seed = (S.seed + 1) % 100000; updateSaturation(); });

    /* lab tabs */
    els('.lab-tab').forEach(b => b.addEventListener('click', () => {
      els('.lab-tab').forEach(x => x.classList.toggle('on', x === b));
      els('.lab').forEach(x => x.classList.toggle('on', x.id === b.dataset.lab));
    }));

    updateFelsen();
    updateSaturation();
    document.addEventListener('langchange', () => { updateFelsen(false); updateSaturation(false); });
    document.addEventListener('themechange', () => { updateFelsen(false); updateSaturation(false); });
  }

  document.addEventListener('DOMContentLoaded', init);

  /* exported for the test bench */
  window.Playground = {
    CLASSES, classProbs, lnLik, fitTopology, expectedInformative, zoneBoundary,
    simulate, pairDistances, njChoice, mpChoice, mlChoice, jcP,
    k80Probs, gammaRates, invGamma, lgamma, expectedPQ,
    F, S, runFelsen, runSaturation,
  };
})();
