/* PhylogenyPro — molecular phylogenetics in the browser.
   Copyright (C) 2026  Luis Ángel Barrera-Guzmán

   This program is free software: you can redistribute it and/or modify it under
   the terms of the GNU General Public License as published by the Free Software
   Foundation, either version 3 of the License, or (at your option) any later
   version. It is distributed in the hope that it will be useful, but WITHOUT ANY
   WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
   A PARTICULAR PURPOSE. See the GNU General Public License, in the file LICENSE
   at the root of this program, or <https://www.gnu.org/licenses/>. */

/* PhylogenyPro — images of the OTUs and of the clades, drawn beside a tree.

   Every image belongs to one name: a tip of the tree (a species, a variety, an
   accession) or a group. It is drawn only once the user has confirmed that it
   shows exactly that taxon — a photograph pinned to the wrong branch is a
   claim, not a decoration.

   Adjustments are applied on a canvas and cached as a PNG data URL, so the
   figure stays self-contained: it exports to SVG, PNG and TIFF identically, and
   the SVG carries its pictures inside it rather than pointing at files that
   will not be there tomorrow.

   The images live in the browser (IndexedDB), keyed by the data set they belong
   to, and never leave the computer — the same promise the rest of the app makes.

   Sister module of PopGeneticsPro's image library, by the same author. */

const OTUImg = {};

(function () {

  const DB_NAME = 'phylogenypro', STORE = 'otuImages';
  const OUT = 640;               // side of the processed square image, px
  const MAX_ORIGINAL = 1400;     // originals are downscaled to this on import

  const recs = new Map();        // name → record
  const listeners = new Set();
  let scope = 'default';
  let dbp = null;

  const ADJ0 = {
    zoom: 1, dx: 0, dy: 0, rotate: 0,
    brightness: 0, contrast: 0, saturation: 0, warmth: 0, sharpness: 0,
    gray: false, removeBg: false, bgTol: 28,
  };

  /* ================================================================
     1 · persistence
     ================================================================ */
  function db() {
    if (dbp) return dbp;
    dbp = new Promise(resolve => {
      try {
        const req = indexedDB.open(DB_NAME, 1);
        req.onupgradeneeded = () => {
          if (!req.result.objectStoreNames.contains(STORE)) req.result.createObjectStore(STORE);
        };
        req.onsuccess = () => resolve(req.result);
        req.onerror = () => resolve(null);
      } catch (e) { resolve(null); }
    });
    return dbp;
  }
  const keyOf = name => scope + '::' + name;

  async function persist(rec) {
    const d = await db(); if (!d) return;
    try {
      const tx = d.transaction(STORE, 'readwrite');
      tx.objectStore(STORE).put(Object.assign({}, rec), keyOf(rec.name));
    } catch (e) { /* storage full or blocked: the record stays in memory */ }
  }
  async function unpersist(name) {
    const d = await db(); if (!d) return;
    try { d.transaction(STORE, 'readwrite').objectStore(STORE).delete(keyOf(name)); } catch (e) { /* ignore */ }
  }

  /* Images belong to a data set, so opening another one does not drag the
     previous pictures along. */
  async function useScope(s) {
    scope = String(s || 'default');
    recs.clear();
    const d = await db();
    if (d) {
      await new Promise(resolve => {
        try {
          const st = d.transaction(STORE, 'readonly').objectStore(STORE);
          const req = st.openCursor();
          req.onsuccess = () => {
            const cur = req.result;
            if (!cur) return resolve();
            if (String(cur.key).startsWith(scope + '::')) recs.set(cur.value.name, cur.value);
            cur.continue();
          };
          req.onerror = () => resolve();
        } catch (e) { resolve(); }
      });
    }
    emit();
  }

  function emit() { listeners.forEach(fn => { try { fn(); } catch (e) { console.error(e); } }); }
  function onChange(fn) { listeners.add(fn); return () => listeners.delete(fn); }

  /* ================================================================
     2 · import and processing
     ================================================================ */
  function loadImage(src) {
    return new Promise((resolve, reject) => {
      const im = new Image();
      im.onload = () => resolve(im);
      im.onerror = () => reject(new Error('unreadable'));
      im.src = src;
    });
  }
  function fileToDataURL(file) {
    return new Promise((resolve, reject) => {
      const fr = new FileReader();
      fr.onload = () => resolve(fr.result);
      fr.onerror = () => reject(new Error('unreadable'));
      fr.readAsDataURL(file);
    });
  }

  async function importFile(file) {
    if (!/^image\//.test(file.type) && !/\.(png|jpe?g|webp|gif|bmp|svg)$/i.test(file.name))
      throw new Error('not_an_image');
    const src = await fileToDataURL(file);
    const im = await loadImage(src);
    const k = Math.min(1, MAX_ORIGINAL / Math.max(im.naturalWidth, im.naturalHeight));
    const c = document.createElement('canvas');
    c.width = Math.max(1, Math.round(im.naturalWidth * k));
    c.height = Math.max(1, Math.round(im.naturalHeight * k));
    c.getContext('2d').drawImage(im, 0, 0, c.width, c.height);
    const keepAlpha = /png|webp|gif|svg/i.test(file.type) || /\.(png|webp|gif|svg)$/i.test(file.name);
    return {
      original: keepAlpha ? c.toDataURL('image/png') : c.toDataURL('image/jpeg', 0.92),
      w: c.width, h: c.height, fileName: file.name,
    };
  }

  const clamp = v => (v < 0 ? 0 : v > 255 ? 255 : v);

  /* A flood fill from the border: a plain background — paper, a sheet, the sky —
     becomes transparent, while light areas *inside* the object are kept, which is
     what a threshold on brightness alone would destroy. */
  function removeBackground(px, W, H, tol) {
    const n = W * H;
    const idx = new Uint8Array(n);                 // 1 = background
    const border = [];
    for (let x = 0; x < W; x++) border.push(x, (H - 1) * W + x);
    for (let y = 0; y < H; y++) border.push(y * W, y * W + W - 1);
    const comp = ch => {
      const v = border.map(p => px[p * 4 + ch]).sort((a, b) => a - b);
      return v[v.length >> 1];                     // median of the border
    };
    const ref = [comp(0), comp(1), comp(2)];
    const t2 = Math.pow(tol * 4.42, 2);            // tol 0–100 → distance in RGB
    const near = p => {
      if (px[p * 4 + 3] < 10) return true;
      const dr = px[p * 4] - ref[0], dg = px[p * 4 + 1] - ref[1], db2 = px[p * 4 + 2] - ref[2];
      return dr * dr + dg * dg + db2 * db2 <= t2;
    };
    const stack = [];
    border.forEach(p => { if (!idx[p] && near(p)) { idx[p] = 1; stack.push(p); } });
    while (stack.length) {
      const p = stack.pop(), x = p % W, y = (p / W) | 0;
      const nb = [x > 0 ? p - 1 : -1, x < W - 1 ? p + 1 : -1, y > 0 ? p - W : -1, y < H - 1 ? p + W : -1];
      for (const q of nb) if (q >= 0 && !idx[q] && near(q)) { idx[q] = 1; stack.push(q); }
    }
    for (let p = 0; p < n; p++) {
      if (idx[p]) { px[p * 4 + 3] = 0; continue; }
      /* soften the outline: a pixel touching the background goes partly transparent */
      const x = p % W, y = (p / W) | 0;
      let touch = 0;
      if (x > 0 && idx[p - 1]) touch++;
      if (x < W - 1 && idx[p + 1]) touch++;
      if (y > 0 && idx[p - W]) touch++;
      if (y < H - 1 && idx[p + W]) touch++;
      if (touch) px[p * 4 + 3] = Math.round(px[p * 4 + 3] * (1 - 0.18 * touch));
    }
  }

  /* separable 1-2-1 blur, the reference for the unsharp mask */
  function blur3(src, W, H) {
    const out = new Float32Array(src.length);
    const tmp = new Float32Array(src.length);
    const k = [1, 2, 1];
    for (let y = 0; y < H; y++) for (let x = 0; x < W; x++) for (let ch = 0; ch < 3; ch++) {
      let s = 0, w = 0;
      for (let d = -1; d <= 1; d++) {
        const xx = Math.min(W - 1, Math.max(0, x + d));
        s += src[(y * W + xx) * 4 + ch] * k[d + 1]; w += k[d + 1];
      }
      tmp[(y * W + x) * 4 + ch] = s / w;
    }
    for (let y = 0; y < H; y++) for (let x = 0; x < W; x++) for (let ch = 0; ch < 3; ch++) {
      let s = 0, w = 0;
      for (let d = -1; d <= 1; d++) {
        const yy = Math.min(H - 1, Math.max(0, y + d));
        s += tmp[(yy * W + x) * 4 + ch] * k[d + 1]; w += k[d + 1];
      }
      out[(y * W + x) * 4 + ch] = s / w;
    }
    return out;
  }

  async function processRecord(rec, size) {
    size = size || OUT;
    const a = Object.assign({}, ADJ0, rec.adj || {});
    const im = await loadImage(rec.original);
    const c = document.createElement('canvas');
    c.width = size; c.height = size;
    const ctx = c.getContext('2d');
    const W = im.naturalWidth, H = im.naturalHeight;
    const base = Math.max(size / W, size / H) * a.zoom;   // cover the square
    const dw = W * base, dh = H * base;
    ctx.save();
    ctx.translate(size / 2, size / 2);
    ctx.rotate((a.rotate || 0) * Math.PI / 180);
    /* dx, dy in −1…1 move the picture by up to half of its overflow, or half of
       the frame when the picture is smaller than the frame */
    const ox = Math.max((dw - size) / 2, size / 2) * a.dx;
    const oy = Math.max((dh - size) / 2, size / 2) * a.dy;
    ctx.imageSmoothingQuality = 'high';
    ctx.drawImage(im, -dw / 2 - ox, -dh / 2 - oy, dw, dh);
    ctx.restore();

    const data = ctx.getImageData(0, 0, size, size);
    const px = data.data;
    if (a.removeBg) removeBackground(px, size, size, a.bgTol);

    const b = a.brightness * 2.55;
    const cc = a.contrast * 2.55, cf = (259 * (cc + 255)) / (255 * (259 - cc));
    const sat = 1 + a.saturation / 100;
    const warm = a.warmth * 0.6;
    for (let p = 0; p < px.length; p += 4) {
      if (px[p + 3] === 0) continue;
      let r = px[p], g = px[p + 1], bl = px[p + 2];
      r = cf * (r + b - 128) + 128; g = cf * (g + b - 128) + 128; bl = cf * (bl + b - 128) + 128;
      const L = 0.2126 * r + 0.7152 * g + 0.0722 * bl;
      if (a.gray) { r = g = bl = L; }
      else { r = L + (r - L) * sat; g = L + (g - L) * sat; bl = L + (bl - L) * sat; }
      r += warm; bl -= warm;
      px[p] = clamp(r); px[p + 1] = clamp(g); px[p + 2] = clamp(bl);
    }

    /* sharpness: an unsharp mask above zero, a gentle blur below */
    if (a.sharpness) {
      const blr = blur3(px, size, size);
      const amt = a.sharpness / 50;
      for (let p = 0; p < px.length; p += 4) {
        if (px[p + 3] === 0) continue;
        for (let ch = 0; ch < 3; ch++) {
          const v = px[p + ch], m = blr[p + ch];
          px[p + ch] = clamp(amt > 0 ? v + amt * (v - m) : v + (-amt / 2) * (m - v));
        }
      }
    }
    ctx.putImageData(data, 0, 0);
    return c.toDataURL('image/png');
  }

  /* ================================================================
     3 · the records
     ================================================================ */
  function get(name) { return recs.get(name) || null; }

  /* the picture to draw for a name — only when it has been confirmed */
  function url(name) {
    const r = recs.get(name);
    return r && r.confirmed && r.processed ? r.processed : null;
  }
  function names() { return [...recs.keys()]; }
  function count() { let n = 0; recs.forEach(r => { if (r.confirmed && r.processed) n++; }); return n; }

  async function setImage(name, file) {
    const imp = await importFile(file);
    const old = recs.get(name);
    const rec = {
      name, original: imp.original, w: imp.w, h: imp.h, fileName: imp.fileName,
      adj: Object.assign({}, ADJ0),
      credit: old ? old.credit : '', license: old ? old.license : '',
      confirmed: false, processed: null, updated: Date.now(),
    };
    rec.processed = await processRecord(rec);
    recs.set(name, rec); persist(rec); emit();
    return rec;
  }

  async function update(name, patch) {
    const rec = recs.get(name); if (!rec) return null;
    if (patch.adj) rec.adj = Object.assign({}, rec.adj, patch.adj);
    ['credit', 'license', 'confirmed'].forEach(k => { if (patch[k] !== undefined) rec[k] = patch[k]; });
    if (patch.adj) rec.processed = await processRecord(rec);
    rec.updated = Date.now();
    persist(rec); emit();
    return rec;
  }

  function remove(name) { recs.delete(name); unpersist(name); emit(); }

  async function removeAll() {
    const all = [...recs.keys()];
    recs.clear();
    for (const n of all) await unpersist(n);
    emit();
  }

  /* file names → OTU names: "Quercus_rugosa.jpg" finds "Quercus rugosa" */
  const norm = s => String(s).toLowerCase().normalize('NFD')
    .replace(/[̀-ͯ]/g, '')
    .replace(/\.[a-z0-9]{2,4}$/, '')
    .replace(/[^a-z0-9]+/g, '');

  function matchName(fileName, candidates) {
    const f = norm(fileName);
    if (!f) return null;
    let hit = candidates.find(c => norm(c) === f);
    if (hit) return hit;
    hit = candidates.filter(c => norm(c).length >= 3 && (f.startsWith(norm(c)) || f.includes(norm(c))));
    return hit.length === 1 ? hit[0] : null;
  }

  /* the credit line a figure has to carry, for the ones that have one */
  function credits(list) {
    const out = [];
    (list || names()).forEach(n => {
      const r = recs.get(n);
      if (r && r.confirmed && (r.credit || r.license)) {
        out.push(`${n}: ${[r.credit, r.license].filter(Boolean).join(', ')}`);
      }
    });
    return out;
  }

  /* ================================================================
     4 · the manager
     ================================================================ */
  let modal = null, list = [], current = null, onDone = null;

  const TT = (es, en) => (typeof T === 'function' ? T(es, en) : es);

  function buildModal() {
    modal = document.createElement('div');
    modal.className = 'oi-back';
    modal.innerHTML = `
      <div class="oi-modal" role="dialog" aria-modal="true">
        <div class="oi-head">
          <h3 class="oi-title"></h3>
          <button class="btn btn-ghost btn-sm oi-close">✕</button>
        </div>
        <div class="oi-body">
          <div class="oi-list">
            <div class="oi-drop">
              <b class="oi-droptitle"></b>
              <p class="oi-dropsub"></p>
              <button class="btn btn-secondary btn-sm oi-pick"></button>
              <input type="file" accept="image/*" multiple hidden class="oi-files">
            </div>
            <div class="oi-names"></div>
          </div>
          <div class="oi-editor"></div>
        </div>
      </div>`;
    document.body.appendChild(modal);
    modal.querySelector('.oi-close').addEventListener('click', close);
    modal.addEventListener('click', e => { if (e.target === modal) close(); });
    document.addEventListener('keydown', e => {
      if (e.key === 'Escape' && modal && modal.style.display !== 'none') close();
    });
    const files = modal.querySelector('.oi-files');
    modal.querySelector('.oi-pick').addEventListener('click', () => files.click());
    files.addEventListener('change', () => { takeFiles([...files.files]); files.value = ''; });
    const drop = modal.querySelector('.oi-drop');
    ['dragenter', 'dragover'].forEach(ev => drop.addEventListener(ev, e => {
      e.preventDefault(); drop.classList.add('over');
    }));
    ['dragleave', 'drop'].forEach(ev => drop.addEventListener(ev, e => {
      e.preventDefault(); drop.classList.remove('over');
      if (ev === 'drop') takeFiles([...(e.dataTransfer.files || [])]);
    }));
  }

  /* Several files at once: each one goes to the taxon its file name names.
     Anything ambiguous is left for the user rather than guessed. */
  async function takeFiles(files) {
    if (!files.length) return;
    const unmatched = [];
    for (const f of files) {
      let target = files.length === 1 && current ? current : matchName(f.name, list);
      if (!target) { unmatched.push(f.name); continue; }
      try { await setImage(target, f); current = target; } catch (e) { unmatched.push(f.name); }
    }
    render();
    if (unmatched.length) {
      const box = modal.querySelector('.oi-names');
      const p = document.createElement('p');
      p.className = 'oi-warn';
      p.textContent = TT(
        `No supe a qué taxón corresponde: ${unmatched.join(', ')}. Elige el nombre en la lista y vuelve a soltar el archivo.`,
        `I could not tell which taxon these belong to: ${unmatched.join(', ')}. Pick the name in the list and drop the file again.`);
      box.insertBefore(p, box.firstChild);
    }
  }

  function stateOf(nm) {
    const r = recs.get(nm);
    if (!r) return { cls: 'none', txt: TT('sin imagen', 'no picture') };
    if (!r.confirmed) return { cls: 'warn', txt: TT('sin confirmar', 'unconfirmed') };
    return { cls: 'ok', txt: TT('lista', 'ready') };
  }

  function render() {
    if (!modal) return;
    modal.querySelector('.oi-title').textContent = TT('Imágenes de los taxones', 'Pictures of the taxa');
    modal.querySelector('.oi-droptitle').textContent = TT('Arrastra aquí tus fotografías', 'Drag your photographs here');
    modal.querySelector('.oi-dropsub').textContent = TT(
      'Si el archivo se llama como el taxón (Bursera_simaruba.jpg), va solo a su lugar.',
      'A file named after its taxon (Bursera_simaruba.jpg) finds its own place.');
    modal.querySelector('.oi-pick').textContent = TT('Elegir archivos…', 'Choose files…');

    const box = modal.querySelector('.oi-names');
    box.innerHTML = '';
    list.forEach(nm => {
      const st = stateOf(nm);
      const r = recs.get(nm);
      const row = document.createElement('button');
      row.className = 'oi-row' + (nm === current ? ' on' : '');
      row.innerHTML =
        `<span class="oi-thumb">${r && r.processed ? `<img src="${r.processed}" alt="">` : ''}</span>` +
        `<span class="oi-nm"></span><span class="oi-state ${st.cls}"></span>`;
      row.querySelector('.oi-nm').textContent = nm;
      row.querySelector('.oi-state').textContent = st.txt;
      row.addEventListener('click', () => { current = nm; render(); });
      box.appendChild(row);
    });
    renderEditor();
  }

  const SLIDERS = [
    ['zoom', 'acercamiento', 'zoom', 1, 4, 0.01],
    ['dx', 'mover ←→', 'move ←→', -1, 1, 0.01],
    ['dy', 'mover ↑↓', 'move ↑↓', -1, 1, 0.01],
    ['rotate', 'girar', 'rotate', -180, 180, 1],
    ['brightness', 'brillo', 'brightness', -60, 60, 1],
    ['contrast', 'contraste', 'contrast', -60, 60, 1],
    ['saturation', 'saturación', 'saturation', -100, 100, 1],
    ['warmth', 'calidez', 'warmth', -60, 60, 1],
    ['sharpness', 'nitidez', 'sharpness', -50, 100, 1],
  ];

  function renderEditor() {
    const host = modal.querySelector('.oi-editor');
    host.innerHTML = '';
    if (!current) {
      host.innerHTML = `<p class="oi-empty">${esc(TT('Elige un taxón de la lista.', 'Pick a taxon from the list.'))}</p>`;
      return;
    }
    const rec = recs.get(current);
    const h = document.createElement('div');
    h.innerHTML = `<h4 class="oi-eh"></h4>`;
    h.querySelector('.oi-eh').textContent = current;
    host.appendChild(h);

    if (!rec) {
      const p = document.createElement('p');
      p.className = 'oi-empty';
      p.textContent = TT('Todavía no tiene imagen. Suelta un archivo o elígelo arriba.',
                         'It has no picture yet. Drop a file or choose one above.');
      host.appendChild(p);
      return;
    }

    const prev = document.createElement('div');
    prev.className = 'oi-prev';
    prev.innerHTML = `<img src="${rec.processed}" alt="">`;
    host.appendChild(prev);

    const adj = Object.assign({}, ADJ0, rec.adj);
    const ctr = document.createElement('div');
    ctr.className = 'oi-ctrls';
    SLIDERS.forEach(([k, es, en, lo, hi, step]) => {
      const row = document.createElement('label');
      row.className = 'oi-sl';
      row.innerHTML = `<span></span><input type="range" min="${lo}" max="${hi}" step="${step}" value="${adj[k]}"><b></b>`;
      row.querySelector('span').textContent = TT(es, en);
      row.querySelector('b').textContent = adj[k];
      const inp = row.querySelector('input');
      inp.addEventListener('input', () => { row.querySelector('b').textContent = inp.value; });
      inp.addEventListener('change', async () => {
        await update(current, { adj: { [k]: +inp.value } });
        render();
      });
      ctr.appendChild(row);
    });
    host.appendChild(ctr);

    const tog = document.createElement('div');
    tog.className = 'oi-togs';
    const mkTog = (label, checked, fn) => {
      const l = document.createElement('label');
      l.className = 'inline-label';
      l.innerHTML = `<input type="checkbox"${checked ? ' checked' : ''}><span></span>`;
      l.querySelector('span').textContent = label;
      l.querySelector('input').addEventListener('change', e => fn(e.target.checked));
      tog.appendChild(l);
      return l;
    };
    mkTog(TT('blanco y negro', 'greyscale'), adj.gray, v => update(current, { adj: { gray: v } }).then(render));
    mkTog(TT('quitar el fondo', 'remove the background'), adj.removeBg, v => update(current, { adj: { removeBg: v } }).then(render));
    if (adj.removeBg) {
      const l = document.createElement('label');
      l.className = 'oi-sl';
      l.innerHTML = `<span></span><input type="range" min="4" max="90" step="1" value="${adj.bgTol}"><b>${adj.bgTol}</b>`;
      l.querySelector('span').textContent = TT('tolerancia del fondo', 'background tolerance');
      l.querySelector('input').addEventListener('change', e => update(current, { adj: { bgTol: +e.target.value } }).then(render));
      tog.appendChild(l);
    }
    host.appendChild(tog);

    const meta = document.createElement('div');
    meta.className = 'oi-meta';
    meta.innerHTML =
      `<label class="inline-label"><span></span><input type="text" class="oi-credit" value="${esc(rec.credit || '')}"></label>` +
      `<label class="inline-label"><span></span><input type="text" class="oi-lic" value="${esc(rec.license || '')}"></label>`;
    meta.querySelectorAll('span')[0].textContent = TT('crédito', 'credit');
    meta.querySelectorAll('span')[1].textContent = TT('licencia', 'licence');
    meta.querySelector('.oi-credit').addEventListener('change', e => update(current, { credit: e.target.value }));
    meta.querySelector('.oi-lic').addEventListener('change', e => update(current, { license: e.target.value }));
    host.appendChild(meta);

    /* the confirmation: nothing is drawn until this is ticked */
    const conf = document.createElement('label');
    conf.className = 'oi-confirm' + (rec.confirmed ? ' on' : '');
    conf.innerHTML = `<input type="checkbox"${rec.confirmed ? ' checked' : ''}><span></span>`;
    conf.querySelector('span').textContent = TT(
      `Confirmo que esta imagen es de ${current}.`, `I confirm this picture is of ${current}.`);
    conf.querySelector('input').addEventListener('change', e => update(current, { confirmed: e.target.checked }).then(render));
    host.appendChild(conf);

    const rm = document.createElement('button');
    rm.className = 'btn btn-ghost btn-sm';
    rm.textContent = TT('Quitar esta imagen', 'Remove this picture');
    rm.addEventListener('click', () => { remove(current); render(); });
    host.appendChild(rm);
  }

  function openManager(opts) {
    opts = opts || {};
    list = (opts.names || []).filter(Boolean);
    current = list.indexOf(current) >= 0 ? current : (list[0] || null);
    onDone = opts.onDone || null;
    if (!modal) buildModal();
    modal.style.display = 'flex';
    document.body.style.overflow = 'hidden';
    render();
  }

  function close() {
    if (modal) modal.style.display = 'none';
    document.body.style.overflow = '';
    if (onDone) { try { onDone(); } catch (e) { console.error(e); } }
  }

  Object.assign(OTUImg, {
    ADJ0, OUT, openManager, close,
    useScope, onChange, get, url, names, count,
    setImage, update, remove, removeAll,
    processRecord, importFile, matchName, credits,
    /* exposed for the tests */
    _removeBackground: removeBackground, _blur3: blur3, _norm: norm,
  });
  if (typeof window !== 'undefined') window.OTUImg = OTUImg;

})();
