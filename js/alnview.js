/* PhylogenyPro — molecular phylogenetics in the browser.
   Copyright (C) 2026  Luis Ángel Barrera-Guzmán

   This program is free software: you can redistribute it and/or modify it under
   the terms of the GNU General Public License as published by the Free Software
   Foundation, either version 3 of the License, or (at your option) any later
   version. It is distributed in the hope that it will be useful, but WITHOUT ANY
   WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
   A PARTICULAR PURPOSE. See the GNU General Public License, in the file LICENSE
   at the root of this program, or <https://www.gnu.org/licenses/>. */

/* PhylogenyPro — the alignment viewer and editor.

   Drawn on a canvas, not with DOM nodes: an alignment of 200 taxa × 5000
   columns is a million cells, and only the visible window is ever painted, so
   scrolling stays smooth. The scrollbars are a real, native scrolling div with
   a spacer of the full size, which is what keeps the keyboard and the trackpad
   behaving the way people expect.

   What it can do
     · colour by residue, by agreement with the consensus, or not at all;
     · mark the variable and the parsimony-informative columns;
     · show which columns a trimming would remove, before removing them;
     · show codon positions when the partition is coding;
     · select a column, a residue or a range;
     · edit by hand: insert or delete a gap, delete a column, move a sequence,
       with unlimited undo — because an alignment is a hypothesis and sometimes
       the eye of whoever knows the group beats the algorithm. */

const AlnView = {};

(function () {

  const AA_GROUPS = {
    /* Clustal-like colouring by physicochemical class */
    A: 'hydro', V: 'hydro', L: 'hydro', I: 'hydro', M: 'hydro', F: 'aroma', W: 'aroma', Y: 'aroma',
    K: 'basic', R: 'basic', H: 'basic', D: 'acid', E: 'acid',
    S: 'polar', T: 'polar', N: 'polar', Q: 'polar', C: 'cys', G: 'gly', P: 'pro',
  };

  function create(host, opts) {
    opts = opts || {};
    const state = {
      taxa: [], seqs: [], type: 'dna',
      cw: 10, ch: 16, nameW: 170,
      mode: 'base', showMask: null, marks: null, frame: null,
      sel: { row: -1, col: -1, col2: -1 },
      undo: [], redo: [],
      onChange: opts.onChange || null,
      onSelect: opts.onSelect || null,
    };

    /* --- DOM --- */
    host.innerHTML = '';
    host.classList.add('aln-host');
    const scroller = mk('div', { class: 'aln-scroll' });
    const spacer = mk('div', { class: 'aln-spacer' });
    const canvas = mk('canvas', { class: 'aln-canvas' });
    scroller.appendChild(spacer);
    host.appendChild(canvas);
    host.appendChild(scroller);
    const ctx = canvas.getContext('2d');

    function colours() {
      const css = getComputedStyle(document.documentElement);
      const v = n => css.getPropertyValue(n).trim();
      return {
        A: v('--base-a'), C: v('--base-c'), G: v('--base-g'), T: v('--base-t'), U: v('--base-t'),
        gap: v('--base-gap'),
        hydro: v('--c5'), aroma: v('--c7'), basic: v('--c1'), acid: v('--c4'),
        polar: v('--c3'), cys: v('--c2'), gly: v('--c6'), pro: v('--c8'),
        text: v('--text'), muted: v('--text-muted'), bg: v('--card-bg'), soft: v('--bg-soft'),
        border: v('--border'), primary: v('--primary'), accent: v('--accent'), danger: v('--danger'),
      };
    }

    function layout() {
      /* the box follows the data: a small alignment does not leave a blank half */
      if (!host.dataset.fixed) {
        const want = 28 + state.taxa.length * state.ch + 26;
        host.style.height = Math.max(120, Math.min(440, want)) + 'px';
      }
      const rect = host.getBoundingClientRect();
      const dpr = window.devicePixelRatio || 1;
      const w = Math.max(200, rect.width), h = Math.max(120, rect.height);
      canvas.width = Math.round(w * dpr); canvas.height = Math.round(h * dpr);
      canvas.style.width = w + 'px'; canvas.style.height = h + 'px';
      ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
      const L = state.seqs[0] ? state.seqs[0].length : 0;
      spacer.style.width = (state.nameW + L * state.cw + 20) + 'px';
      spacer.style.height = (28 + state.taxa.length * state.ch + 20) + 'px';
      return { w, h };
    }

    function consensus() {
      const L = state.seqs[0] ? state.seqs[0].length : 0;
      const out = new Array(L);
      for (let c = 0; c < L; c++) {
        const counts = {};
        let best = null, bestN = 0;
        for (let i = 0; i < state.seqs.length; i++) {
          const ch = state.seqs[i][c];
          if (ch === '-' || ch === '?') continue;
          counts[ch] = (counts[ch] || 0) + 1;
          if (counts[ch] > bestN) { bestN = counts[ch]; best = ch; }
        }
        out[c] = best;
      }
      return out;
    }
    let cons = null;

    function draw() {
      const { w, h } = layout();
      const C = colours();
      const sx = scroller.scrollLeft, sy = scroller.scrollTop;
      const L = state.seqs[0] ? state.seqs[0].length : 0;
      const n = state.taxa.length;
      ctx.clearRect(0, 0, w, h);
      ctx.fillStyle = C.bg; ctx.fillRect(0, 0, w, h);
      if (!n) { ctx.fillStyle = C.muted; ctx.font = '13px system-ui'; ctx.fillText('—', 12, 24); return; }

      const headH = 28;
      const firstCol = Math.max(0, Math.floor(sx / state.cw));
      const lastCol = Math.min(L - 1, Math.ceil((sx + w - state.nameW) / state.cw));
      const firstRow = Math.max(0, Math.floor(sy / state.ch));
      const lastRow = Math.min(n - 1, Math.ceil((sy + h - headH) / state.ch));
      const X = c => state.nameW + c * state.cw - sx;
      const Y = r => headH + r * state.ch - sy;

      /* --- ruler and marks --- */
      ctx.fillStyle = C.soft; ctx.fillRect(0, 0, w, headH);
      ctx.font = '10px ui-monospace, monospace';
      ctx.textBaseline = 'middle';
      for (let c = firstCol; c <= lastCol; c++) {
        const x = X(c);
        if (state.showMask && !state.showMask[c]) { ctx.fillStyle = C.danger; ctx.globalAlpha = 0.16; ctx.fillRect(x, headH, state.cw, h - headH); ctx.globalAlpha = 1; }
        if (state.marks) {
          const m = state.marks[c];
          if (m) { ctx.fillStyle = m === 2 ? C.accent : C.primary; ctx.fillRect(x, headH - 5, Math.max(1, state.cw - 1), 3); }
        }
        if (state.frame && state.cw >= 7) {
          const pos = ((c - (state.frame - 1)) % 3 + 3) % 3;
          ctx.fillStyle = pos === 2 ? C.accent : C.border;
          ctx.globalAlpha = pos === 2 ? 0.5 : 0.35;
          ctx.fillRect(x, headH - 9, Math.max(1, state.cw - 1), 3);
          ctx.globalAlpha = 1;
        }
        const step = state.cw >= 9 ? 10 : state.cw >= 5 ? 20 : 50;
        if ((c + 1) % step === 0) {
          ctx.fillStyle = C.muted;
          ctx.fillText(String(c + 1), x, 10);
          ctx.strokeStyle = C.border; ctx.beginPath(); ctx.moveTo(x, 16); ctx.lineTo(x, 20); ctx.stroke();
        }
      }

      /* --- residues --- */
      const small = state.cw < 7;
      ctx.font = Math.min(state.ch - 3, state.cw + 2) + 'px ui-monospace, monospace';
      ctx.textAlign = 'center';
      for (let r = firstRow; r <= lastRow; r++) {
        const seq = state.seqs[r], y = Y(r);
        for (let c = firstCol; c <= lastCol; c++) {
          const ch = seq[c] || '-';
          const x = X(c);
          let fill = null;
          if (ch === '-' || ch === '?') fill = null;
          else if (state.mode === 'base') {
            fill = state.type === 'aa' ? C[AA_GROUPS[ch] || 'polar'] : C[ch];
          } else if (state.mode === 'identity') {
            if (!cons) cons = consensus();
            fill = ch === cons[c] ? C.border : C.accent;
          }
          if (fill) {
            ctx.fillStyle = fill;
            ctx.globalAlpha = state.mode === 'identity' ? (ch === cons[c] ? 0.35 : 0.85) : 0.82;
            ctx.fillRect(x, y, Math.max(1, state.cw - 0.6), state.ch - 1);
            ctx.globalAlpha = 1;
          } else if (ch === '-') {
            ctx.strokeStyle = C.border;
            ctx.beginPath(); ctx.moveTo(x + 1, y + state.ch / 2); ctx.lineTo(x + state.cw - 1.5, y + state.ch / 2); ctx.stroke();
          }
          if (!small) {
            ctx.fillStyle = fill ? '#fff' : C.muted;
            if (state.mode === 'identity' && fill === C.border) ctx.fillStyle = C.muted;
            if (ch !== '-') ctx.fillText(ch, x + state.cw / 2, y + state.ch / 2);
          }
        }
      }

      /* --- selection --- */
      const s = state.sel;
      if (s.col >= 0) {
        const a = Math.min(s.col, s.col2 < 0 ? s.col : s.col2), b = Math.max(s.col, s.col2 < 0 ? s.col : s.col2);
        ctx.strokeStyle = C.primary; ctx.lineWidth = 2;
        ctx.strokeRect(X(a), headH, (b - a + 1) * state.cw, Math.min(h - headH, n * state.ch - sy));
        ctx.lineWidth = 1;
      }
      if (s.row >= 0) {
        ctx.fillStyle = C.primary; ctx.globalAlpha = 0.08;
        ctx.fillRect(state.nameW, Y(s.row), w - state.nameW, state.ch);
        ctx.globalAlpha = 1;
      }

      /* --- names, drawn last so they sit on top --- */
      ctx.fillStyle = C.bg; ctx.fillRect(0, headH, state.nameW, h - headH);
      ctx.strokeStyle = C.border; ctx.beginPath(); ctx.moveTo(state.nameW - 0.5, 0); ctx.lineTo(state.nameW - 0.5, h); ctx.stroke();
      ctx.textAlign = 'left';
      ctx.font = Math.min(12, state.ch - 3) + 'px system-ui';
      for (let r = firstRow; r <= lastRow; r++) {
        const y = Y(r);
        if (r === state.sel.row) { ctx.fillStyle = C.primary; ctx.globalAlpha = 0.12; ctx.fillRect(0, y, state.nameW, state.ch); ctx.globalAlpha = 1; }
        ctx.fillStyle = r === state.sel.row ? C.primary : C.text;
        let nm = state.taxa[r];
        const maxW = state.nameW - 12;
        while (ctx.measureText(nm).width > maxW && nm.length > 4) nm = nm.slice(0, -2) + '…';
        ctx.fillText(nm, 6, y + state.ch / 2);
      }
      ctx.fillStyle = C.soft; ctx.fillRect(0, 0, state.nameW, headH);
      ctx.fillStyle = C.muted; ctx.font = '11px system-ui';
      ctx.fillText(`${n} × ${L}`, 6, 14);
    }

    /* --- interaction --- */
    function cellAt(ev) {
      const rect = host.getBoundingClientRect();
      const x = ev.clientX - rect.left + scroller.scrollLeft - state.nameW;
      const y = ev.clientY - rect.top + scroller.scrollTop - 28;
      return { col: Math.floor(x / state.cw), row: Math.floor(y / state.ch) };
    }
    let dragging = false;
    scroller.addEventListener('scroll', draw);
    scroller.addEventListener('mousedown', ev => {
      const { col, row } = cellAt(ev);
      const L = state.seqs[0] ? state.seqs[0].length : 0;
      if (col < 0 || col >= L) { state.sel = { row: row >= 0 && row < state.taxa.length ? row : -1, col: -1, col2: -1 }; draw(); notifySel(); return; }
      state.sel = { row: (row >= 0 && row < state.taxa.length) ? row : -1, col, col2: col };
      dragging = true;
      draw(); notifySel();
    });
    scroller.addEventListener('mousemove', ev => {
      if (!dragging) return;
      const { col } = cellAt(ev);
      const L = state.seqs[0] ? state.seqs[0].length : 0;
      state.sel.col2 = Math.max(0, Math.min(L - 1, col));
      draw(); notifySel();
    });
    window.addEventListener('mouseup', () => { dragging = false; });
    scroller.addEventListener('wheel', ev => {
      if (!ev.ctrlKey) return;
      ev.preventDefault();
      zoom(ev.deltaY < 0 ? 1 : -1);
    }, { passive: false });

    function notifySel() {
      if (!state.onSelect) return;
      const s = state.sel;
      state.onSelect({
        row: s.row, col: s.col, col2: s.col2,
        taxon: s.row >= 0 ? state.taxa[s.row] : null,
        residue: (s.row >= 0 && s.col >= 0) ? state.seqs[s.row][s.col] : null,
      });
    }

    /* --- editing, with undo --- */
    function snapshot() {
      state.undo.push({ taxa: state.taxa.slice(), seqs: state.seqs.slice() });
      if (state.undo.length > 200) state.undo.shift();
      state.redo.length = 0;
    }
    function changed() { cons = null; draw(); if (state.onChange) state.onChange(); }

    const api = {
      setData(taxa, seqs, type) {
        state.taxa = taxa.slice(); state.seqs = seqs.slice(); state.type = type || 'dna';
        state.sel = { row: -1, col: -1, col2: -1 };
        state.undo.length = 0; state.redo.length = 0;
        cons = null;
        scroller.scrollLeft = 0; scroller.scrollTop = 0;
        draw();
      },
      get() { return { taxa: state.taxa.slice(), seqs: state.seqs.slice() }; },
      setMode(m) { state.mode = m; cons = null; draw(); },
      setMask(mask) { state.showMask = mask; draw(); },
      setMarks(marks) { state.marks = marks; draw(); },
      setFrame(f) { state.frame = f; draw(); },
      setNameWidth(w) { state.nameW = w; draw(); },
      zoom(d) { zoom(d); },
      redraw: draw,
      selection() { return Object.assign({}, state.sel); },
      canUndo() { return state.undo.length > 0; },
      undo() {
        if (!state.undo.length) return false;
        state.redo.push({ taxa: state.taxa.slice(), seqs: state.seqs.slice() });
        const s = state.undo.pop();
        state.taxa = s.taxa; state.seqs = s.seqs;
        changed(); return true;
      },
      redoLast() {
        if (!state.redo.length) return false;
        state.undo.push({ taxa: state.taxa.slice(), seqs: state.seqs.slice() });
        const s = state.redo.pop();
        state.taxa = s.taxa; state.seqs = s.seqs;
        changed(); return true;
      },
      /* insert a gap in one sequence (or in all of them) at the selected column */
      insertGap(all) {
        const s = state.sel;
        if (s.col < 0) return false;
        snapshot();
        state.seqs = state.seqs.map((seq, i) => (all || i === s.row) ? seq.slice(0, s.col) + '-' + seq.slice(s.col) : seq);
        if (!all) padRight();
        changed(); return true;
      },
      /* delete the selected column, or the gap under the cursor in one sequence */
      deleteColumn() {
        const s = state.sel;
        if (s.col < 0) return false;
        const a = Math.min(s.col, s.col2 < 0 ? s.col : s.col2), b = Math.max(s.col, s.col2 < 0 ? s.col : s.col2);
        snapshot();
        state.seqs = state.seqs.map(seq => seq.slice(0, a) + seq.slice(b + 1));
        state.sel = { row: s.row, col: Math.min(a, (state.seqs[0] || '').length - 1), col2: -1 };
        changed(); return true;
      },
      deleteGap() {
        const s = state.sel;
        if (s.col < 0 || s.row < 0) return false;
        if (state.seqs[s.row][s.col] !== '-') return false;
        snapshot();
        state.seqs[s.row] = state.seqs[s.row].slice(0, s.col) + state.seqs[s.row].slice(s.col + 1) + '-';
        changed(); return true;
      },
      /* drop every column that is all gaps */
      squeeze() {
        snapshot();
        state.seqs = Align.squeeze(state.seqs);
        changed(); return true;
      },
      moveRow(dir) {
        const r = state.sel.row;
        if (r < 0) return false;
        const j = r + dir;
        if (j < 0 || j >= state.taxa.length) return false;
        snapshot();
        const t = state.taxa[r]; state.taxa[r] = state.taxa[j]; state.taxa[j] = t;
        const s2 = state.seqs[r]; state.seqs[r] = state.seqs[j]; state.seqs[j] = s2;
        state.sel.row = j;
        changed(); return true;
      },
      removeRow() {
        const r = state.sel.row;
        if (r < 0) return false;
        snapshot();
        state.taxa.splice(r, 1); state.seqs.splice(r, 1);
        state.sel.row = Math.min(r, state.taxa.length - 1);
        changed(); return true;
      },
      sortByName() {
        snapshot();
        const idx = state.taxa.map((t, i) => i).sort((a, b) => state.taxa[a].localeCompare(state.taxa[b]));
        state.taxa = idx.map(i => state.taxa[i]);
        state.seqs = idx.map(i => state.seqs[i]);
        changed(); return true;
      },
      goTo(col) {
        const L = state.seqs[0] ? state.seqs[0].length : 0;
        col = Math.max(0, Math.min(L - 1, col));
        state.sel = { row: state.sel.row, col, col2: col };
        scroller.scrollLeft = Math.max(0, col * state.cw - 200);
        draw(); notifySel();
      },
    };
    function padRight() {
      const L = Math.max(...state.seqs.map(s => s.length));
      state.seqs = state.seqs.map(s => s + '-'.repeat(L - s.length));
    }
    function zoom(d) {
      state.cw = Math.max(2, Math.min(22, state.cw + d));
      state.ch = Math.max(6, Math.min(26, Math.round(state.cw * 1.6)));
      draw();
    }

    document.addEventListener('themechange', draw);
    window.addEventListener('resize', draw);
    draw();
    return api;
  }

  Object.assign(AlnView, { create });
  window.AlnView = AlnView;
})();
