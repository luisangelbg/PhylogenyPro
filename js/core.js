/* PhylogenyPro — molecular phylogenetics in the browser.
   Copyright (C) 2026  Luis Ángel Barrera-Guzmán

   This program is free software: you can redistribute it and/or modify it under
   the terms of the GNU General Public License as published by the Free Software
   Foundation, either version 3 of the License, or (at your option) any later
   version. It is distributed in the hope that it will be useful, but WITHOUT ANY
   WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
   A PARTICULAR PURPOSE. See the GNU General Public License, in the file LICENSE
   at the root of this program, or <https://www.gnu.org/licenses/>. */

/* PhylogenyPro — global state and shared utilities.
   No ES modules: everything hangs from window so the app also works when
   index.html is opened with a double click (file://). */

const state = {
  /* --- Block 2: data --- */
  fileName: null,
  data: null,        // alignment(s), partitions and taxa, read by every other block
  /* --- results, one slot per block --- */
  models: null,      // Block 3: model selection, partitioning scheme, saturation
  quick: null,       // Block 4: distances, NJ/BIONJ/ME, parsimony
  ml: null,          // Block 5: maximum likelihood tree, supports, topology tests
  bayes: null,       // Block 6: posterior sample, consensus, convergence
  dated: null,       // Block 7: rooting, clock tests, divergence times
  diversification: null, // Block 8: LTT, birth-death, rate shifts
  traits: null,      // Block 9: ancestral states, BM/OU, phylogenetic signal
  biogeo: null,      // Block 10: DEC, DEC+J, DIVALIKE, BAYAREALIKE
  compare: null,     // Block 11: tree distances, concordance, species tree, networks
  figures: {},       // every figure registered for the report and the ZIP
};
window.state = state;

/* the twelve blocks of the app, in navigation order */
const STEPS = [
  { n: 1, es: 'Inicio', en: 'Home', ready: true },
  { n: 2, es: 'Datos', en: 'Data', ready: true },
  { n: 3, es: 'Modelos', en: 'Models', ready: true },
  { n: 4, es: 'Parsimonia', en: 'Parsimony', ready: true },
  { n: 5, es: 'Verosimilitud', en: 'Likelihood', ready: true },
  { n: 6, es: 'Bayesiano', en: 'Bayesian', ready: true },
  { n: 7, es: 'Tiempos', en: 'Dating', ready: true },
  { n: 8, es: 'Diversificación', en: 'Diversification', ready: true },
  { n: 9, es: 'Caracteres', en: 'Traits', ready: true },
  { n: 10, es: 'Biogeografía', en: 'Biogeography', ready: true },
  { n: 11, es: 'Comparar árboles', en: 'Compare trees', ready: true },
  { n: 12, es: 'Informe', en: 'Report', ready: true },
];

/* ---------------- DOM ---------------- */
function el(id) { return document.getElementById(id); }
function els(sel, root) { return [...(root || document).querySelectorAll(sel)]; }
function mk(tag, attrs, html) {
  const n = document.createElement(tag);
  if (attrs) for (const k in attrs) {
    if (k === 'class') n.className = attrs[k];
    else if (k === 'style') n.setAttribute('style', attrs[k]);
    else if (k.startsWith('on') && typeof attrs[k] === 'function') n.addEventListener(k.slice(2), attrs[k]);
    else if (attrs[k] != null) n.setAttribute(k, attrs[k]);
  }
  if (html != null) n.innerHTML = html;
  return n;
}
function esc(s) {
  return String(s ?? '').replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;').replace(/"/g, '&quot;');
}
/* bilingual inline HTML: both spans are written, CSS shows the active one */
function L2(es, en) { return `<span data-l="es">${es}</span><span data-l="en">${en}</span>`; }
/* Labels shown in capitals would turn π into Π and Γ into something else:
   Greek letters, symbols with a sub- or superscript and bracketed parameters
   keep their own case inside an uppercase label. */
function keepGreek(s) {
  return String(s).replace(/([Ͱ-Ͽ][²³₀-₉]*[A-Za-z]{0,3}|[A-Za-z][²³₀-₉]+|\[[a-z]\])/g, '<span class="nc">$1</span>');
}
function svgEl(tag, attrs, text) {
  const n = document.createElementNS('http://www.w3.org/2000/svg', tag);
  if (attrs) for (const k in attrs) if (attrs[k] != null) n.setAttribute(k, attrs[k]);
  if (text != null) n.textContent = text;
  return n;
}
function showMessage(container, type, text) {
  if (typeof container === 'string') container = el(container);
  if (!container) return null;
  const div = mk('div', { class: 'msg msg-' + type }, text);
  /* errors and warnings are announced to screen readers */
  if (window.LABG) LABG.messageRole(div, type);
  container.appendChild(div);
  return div;
}
function clearMessages(container) {
  if (typeof container === 'string') container = el(container);
  if (container) container.innerHTML = '';
}

/* ---------------- shared result components ----------------
   Every block reports with the same two pieces: a row of tiles for the handful
   of numbers that matter, and a table for the rest. */
function statTiles(container, tiles) {
  if (typeof container === 'string') container = el(container);
  if (!container) return;
  container.innerHTML = '';
  container.classList.add('results-summary');
  tiles.forEach(t => {
    const [label, value, sub, level] = Array.isArray(t) ? t : [t.label, t.value, t.sub, t.level];
    const d = mk('div', { class: 'stat-tile' + (level ? ' ' + level : '') });
    d.innerHTML = `<div class="stat-label">${keepGreek(label)}</div><div class="stat-value">${value}</div>` +
      (sub ? `<div class="stat-sub">${sub}</div>` : '');
    container.appendChild(d);
  });
}

/* columns: [{key, label, get?, fmt?, num?, html?}] */
function buildTable(container, columns, rows, opts) {
  opts = opts || {};
  if (typeof container === 'string') container = el(container);
  if (!container) return null;
  container.innerHTML = '';
  const table = mk('table');
  if (opts.className) table.className = opts.className;
  if (opts.caption) table.appendChild(mk('caption', null, opts.caption));
  const thead = mk('thead'), trh = mk('tr');
  columns.forEach(c => {
    const th = mk('th', { class: c.num ? 'num' : null });
    th.innerHTML = c.label != null ? c.label : c.key;
    trh.appendChild(th);
  });
  thead.appendChild(trh); table.appendChild(thead);
  const tbody = mk('tbody');
  const shown = opts.limit ? rows.slice(0, opts.limit) : rows;
  shown.forEach(r => {
    const tr = mk('tr');
    if (r && r._class) tr.className = r._class;
    columns.forEach(c => {
      const td = mk('td', { class: c.num ? 'num' : null });
      let v = c.get ? c.get(r) : r[c.key];
      if (c.html) { td.innerHTML = v == null ? '—' : v; tr.appendChild(td); return; }
      if (c.fmt && v != null && v !== '') v = c.fmt(v);
      td.textContent = (v === null || v === undefined || v === '' ||
        (typeof v === 'number' && !isFinite(v))) ? '—' : v;
      tr.appendChild(td);
    });
    tbody.appendChild(tr);
  });
  table.appendChild(tbody);
  container.appendChild(table);
  if (opts.limit && rows.length > opts.limit) {
    const p = mk('p', { class: 'hint', style: 'padding:6px 12px;margin:0' });
    p.innerHTML = L2(`Se muestran ${opts.limit} de ${rows.length} filas.`, `Showing ${opts.limit} of ${rows.length} rows.`);
    container.appendChild(p);
  }
  return table;
}
function tableToCSV(table) {
  const rows = [...table.querySelectorAll('tr')].map(tr =>
    [...tr.children].map(td => csvEscape(td.textContent.trim())).join(','));
  return '﻿' + rows.join('\r\n');
}

/* ---------------- numbers ----------------
   Numbers use the decimal point in both languages (the convention of the
   scientific literature in Mexico and in English). */
function fmtNum(v, d) {
  if (v === null || v === undefined || v === '' || (typeof v === 'number' && !isFinite(v))) return '—';
  const n = Number(v);
  if (!isFinite(n)) return String(v);
  if (n === 0) return '0';
  const abs = Math.abs(n);
  if (abs < 1e-4 || abs >= 1e7) { const e = n.toExponential(d != null ? d : 2); return e.startsWith('-') ? '−' + e.slice(1) : e; }
  const s = n.toLocaleString('en-US', { maximumFractionDigits: d != null ? d : 3 });
  if (/^-0(\.0*)?$/.test(s)) return s.slice(1);   /* a negative value that rounds to zero loses its sign */
  return s.startsWith('-') ? '−' + s.slice(1) : s;
}
function fmtFixed(v, d) {
  if (v === Infinity) return '∞';
  if (v === -Infinity) return '−∞';
  if (v == null || !isFinite(v)) return '—';
  const s = Number(v).toFixed(d == null ? 3 : d);
  if (/^-0(\.0*)?$/.test(s)) return s.slice(1);
  return s.startsWith('-') ? '−' + s.slice(1) : s;
}
function fmtP(p) {
  if (p == null || !isFinite(p)) return '—';
  if (p < 0.0001) return '< 0.0001';
  return Number(p).toFixed(4);
}
/* "p = 0.0123" or "p < 0.0001": the sign that goes with the value */
function pEq(p) { const s = fmtP(p); return s.startsWith('<') ? 'p ' + s : 'p = ' + s; }
function stars(p) {
  if (p == null || !isFinite(p)) return '';
  if (p < 0.001) return '***';
  if (p < 0.01) return '**';
  if (p < 0.05) return '*';
  return 'ns';
}
function fmtPct(x, d) {
  if (x == null || !isFinite(x)) return '—';
  return (x * 100).toLocaleString('en-US', { maximumFractionDigits: d == null ? 1 : d }) + '%';
}
/* a log-likelihood always with the same number of decimals, and a real minus sign */
function fmtLnL(v) { return v == null || !isFinite(v) ? '—' : (v < 0 ? '−' : '') + Math.abs(v).toFixed(3); }
/* "1 taxón", "2 taxones" */
function plural(n, one, many) { return `${n} ${n === 1 ? one : many}`; }

/* ---------------- CSV / downloads ---------------- */
function csvEscape(v) {
  const s = String(v ?? '');
  return /[",\n;]/.test(s) ? '"' + s.replace(/"/g, '""') + '"' : s;
}
function matrixToCSV(header, rows) {
  const head = header.map(csvEscape).join(',');
  const body = rows.map(r => r.map(csvEscape).join(','));
  return '﻿' + [head, ...body].join('\r\n');
}
function download(content, filename, mime) {
  const blob = content instanceof Blob ? content : new Blob([content], { type: mime || 'text/plain;charset=utf-8' });
  const url = URL.createObjectURL(blob);
  const a = mk('a', { href: url, download: filename });
  document.body.appendChild(a); a.click(); document.body.removeChild(a);
  setTimeout(() => URL.revokeObjectURL(url), 2000);
}
function slug(s) {
  return String(s || 'phylogenypro').replace(/\.[^.]+$/, '')
    .normalize('NFD').replace(new RegExp('[\\u0300-\\u036f]', 'g'), '')
    .replace(/[^\w\-]+/g, '_').replace(/^_+|_+$/g, '').slice(0, 60) || 'phylogenypro';
}

/* ---------------- step navigation ---------------- */
function goStep(n) {
  n = String(n);
  els('.step-panel').forEach(p => p.classList.toggle('active', p.id === 'panel-' + n));
  els('.step-btn').forEach(b => b.classList.toggle('active', b.dataset.step === n));
  document.body.classList.toggle('on-home', n === '1');
  window.scrollTo({ top: 0, behavior: 'smooth' });
  const btn = stepBtn(n);
  if (window.LABG) {
    LABG.setCurrentStep(n);
    if (btn) LABG.announce(T('Bloque ', 'Block ') + stepName(n));
  } else if (btn && btn.scrollIntoView) btn.scrollIntoView({ block: 'nearest', inline: 'nearest' });
  refreshStepMarks();
  refreshStepFooters();
  document.dispatchEvent(new CustomEvent('stepchange', { detail: { step: Number(n) } }));
}
function enableStep(n, on) {
  const b = stepBtn(n);
  if (b) b.disabled = (on === false);
  refreshStepMarks();
  refreshStepFooters();
}

/* ---------------- common bar of the LABG Suite ----------------
   Every call to labg-core.js is guarded: the test page loads this file
   without it. */
const stepBtn = n => document.querySelector('.step-btn[data-step="' + n + '"]');
const stepOn = n => { const b = stepBtn(n); return !!b && !b.disabled; };
function stepName(n) {
  const s = STEPS.find(x => String(x.n) === String(n));
  return s ? s.n + ' · ' + T(s.es, s.en) : String(n);
}

/* In this app every block can be opened from the start, so a block counts as
   done when the result it hands on to the others is in `state` (Blocks 1 and
   12 hand nothing on). Recomputed on every enableStep and every change of
   block, so loading new data clears the marks that no longer hold. */
const STEP_RESULT = {
  2: 'data', 3: 'models', 4: 'quick', 5: 'ml', 6: 'bayes', 7: 'dated',
  8: 'diversification', 9: 'traits', 10: 'biogeo', 11: 'compare',
};
function refreshStepMarks() {
  if (!window.LABG) return;
  STEPS.forEach(s => {
    const key = STEP_RESULT[s.n];
    if (!key) return;
    LABG.markStep(s.n, stepOn(s.n) && state[key] ? 'done' : null);
  });
}

/* Footer of every block: Previous / Next, with the name of the block. Both
   languages are written side by side (L2), so a change of language needs no
   redraw. */
function refreshStepFooters() {
  if (!window.LABG) return;
  els('.step-panel').forEach(p => {
    const n = Number(p.id.replace('panel-', ''));
    const i = STEPS.findIndex(s => s.n === n);
    if (i < 0) return;
    let f = p.querySelector(':scope > .step-footer');
    if (!f) {
      f = mk('nav', { class: 'step-footer no-print' });
      f.innerHTML = '<button type="button" class="btn btn-secondary prev"></button><button type="button" class="btn btn-primary next"></button>';
      f.addEventListener('click', e => { const b = e.target.closest('button[data-go]'); if (b && !b.disabled) goStep(b.dataset.go); });
      p.appendChild(f);
    }
    f.setAttribute('aria-label', T('Bloques', 'Blocks'));
    const prev = STEPS.slice(0, i).reverse().find(s => stepOn(s.n));
    const next = STEPS.slice(i + 1).find(s => stepBtn(s.n));
    const label = s => L2(s.n + ' · ' + s.es, s.n + ' · ' + s.en);
    const bp = f.querySelector('.prev'), bn = f.querySelector('.next');
    bp.hidden = !prev;
    if (prev) { bp.dataset.go = prev.n; bp.innerHTML = `← <span><small>${L2('Anterior', 'Previous')}</small>${label(prev)}</span>`; }
    bn.hidden = !next;
    if (next) {
      bn.dataset.go = next.n; bn.disabled = !stepOn(next.n);
      bn.innerHTML = `<span><small>${L2('Siguiente', 'Next')}</small>${label(next)}</span> →`;
    }
  });
}

/* Theme button and block bar: labels that depend on the state and the language */
function paintCommonBar() {
  if (!window.LABG) return;
  LABG.theme.key = 'phylogenypro:theme';
  LABG.theme.paint();
  const nav = el('stepper');
  if (nav) nav.setAttribute('aria-label', T('Bloques', 'Blocks'));
}

/* Wired after every DOMContentLoaded handler has run: home.js builds the
   block bar in its own. */
document.addEventListener('DOMContentLoaded', () => setTimeout(() => {
  if (!window.LABG) return;
  const hb = el('helpBtn');
  if (hb) hb.addEventListener('click', () => LABG.showShortcuts());
  LABG.shortcuts([]);
  LABG.bindStepKeys(goStep);
  LABG.guardUnload(() => !!state.data);
  LABG.setCurrentStep((document.querySelector('.step-btn.active') || {}).dataset?.step || '1');
  document.addEventListener('themechange', paintCommonBar);
  document.addEventListener('langchange', () => { paintCommonBar(); refreshStepFooters(); });
  paintCommonBar();
  refreshStepMarks();
  refreshStepFooters();
}, 0));

/* Persisted preferences (figure style, last settings) */
const Prefs = {
  get(k, d) { try { const v = localStorage.getItem('phylogenypro:' + k); return v == null ? d : JSON.parse(v); } catch (e) { return d; } },
  set(k, v) { try { localStorage.setItem('phylogenypro:' + k, JSON.stringify(v)); } catch (e) { /* ignore */ } },
};

/* ---------------- random numbers ----------------
   Everything stochastic (the simulators of this page, bootstrap resampling,
   the MCMC, stochastic mapping) draws from a seeded generator, so a run is
   reproducible and its seed can be reported in the methods section. sfc32
   passes the usual statistical batteries, unlike a 32-bit linear congruential
   generator, which matters when a chain draws hundreds of millions of values. */
function rng(seed) {
  let a = 0x9e3779b9, b = 0x243f6a88, c = 0xb7e15162, d = (seed >>> 0) ^ 0xdeadbeef;
  const next = () => {
    a >>>= 0; b >>>= 0; c >>>= 0; d >>>= 0;
    let t = (a + b) | 0;
    a = b ^ (b >>> 9);
    b = (c + (c << 3)) | 0;
    c = (c << 21) | (c >>> 11);
    d = (d + 1) | 0;
    t = (t + d) | 0;
    c = (c + t) | 0;
    return (t >>> 0) / 4294967296;
  };
  for (let i = 0; i < 15; i++) next();
  return next;
}
function randn(r) { let u = 0, v = 0; while (u === 0) u = r(); while (v === 0) v = r(); return Math.sqrt(-2 * Math.log(u)) * Math.cos(2 * Math.PI * v); }
/* exponential waiting time with rate lam — the engine of every simulation of
   substitutions along a branch */
function rexp(r, lam) { return -Math.log(1 - r()) / lam; }
function shuffle(arr, r) {
  for (let i = arr.length - 1; i > 0; i--) {
    const j = Math.floor(r() * (i + 1));
    const t = arr[i]; arr[i] = arr[j]; arr[j] = t;
  }
  return arr;
}
function cssVar(name, fallback) {
  const v = getComputedStyle(document.documentElement).getPropertyValue(name).trim();
  return v || fallback || '#2a4a94';
}

/* ---------------- sequences ----------------
   Shared by the illustrations, the simulators and, later, by every analysis
   block: the four bases in a fixed order, so that index 0 is always A. */
const BASES = ['A', 'C', 'G', 'T'];
const BASE_IDX = { A: 0, C: 1, G: 2, T: 3 };
/* a transition keeps the chemical class (purine A↔G, pyrimidine C↔T) */
function isTransition(a, b) {
  return (a === 0 && b === 2) || (a === 2 && b === 0) || (a === 1 && b === 3) || (a === 3 && b === 1);
}
/* Jukes–Cantor correction of an observed proportion of differences.
   Returns null when p ≥ 3/4, where the correction is undefined: at that point
   the sequences are no more similar than two random ones. */
function jcDistance(p) {
  if (p == null || !isFinite(p)) return null;
  if (p <= 0) return 0;
  if (p >= 0.75) return null;
  return -0.75 * Math.log(1 - 4 * p / 3);
}
/* Kimura's two-parameter distance from the proportions of transitions (P)
   and transversions (Q) */
function k2pDistance(P, Q) {
  const a = 1 - 2 * P - Q, b = 1 - 2 * Q;
  if (a <= 0 || b <= 0) return null;
  return -0.5 * Math.log(a) - 0.25 * Math.log(b);
}

Object.assign(window, {
  STEPS, el, els, mk, esc, L2, keepGreek, svgEl, showMessage, clearMessages, statTiles, buildTable, tableToCSV,
  fmtNum, fmtFixed, fmtP, pEq, fmtPct, fmtLnL, plural, stars, csvEscape, matrixToCSV, download, slug,
  goStep, enableStep, Prefs, rng, randn, rexp, shuffle, cssVar,
  BASES, BASE_IDX, isTransition, jcDistance, k2pDistance,
});
