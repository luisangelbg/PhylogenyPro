/* PhylogenyPro — molecular phylogenetics in the browser.
   Copyright (C) 2026  Luis Ángel Barrera-Guzmán

   This program is free software: you can redistribute it and/or modify it under
   the terms of the GNU General Public License as published by the Free Software
   Foundation, either version 3 of the License, or (at your option) any later
   version. It is distributed in the hope that it will be useful, but WITHOUT ANY
   WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
   A PARTICULAR PURPOSE. See the GNU General Public License, in the file LICENSE
   at the root of this program, or <https://www.gnu.org/licenses/>. */

/* PhylogenyPro — getting the work out of the browser.

   A result that cannot leave the program is not a result. This module writes
   the four things a phylogenetic study has to hand over:

     trees      Newick, NEXUS with a TRANSLATE block and a FigTree block, and
                phyloXML, which is the only one of the three that can carry a
                node's support, its age and its confidence interval at once
     figures    the SVG the app drew, and a PNG rendered from it at whatever
                resolution a journal asks for — 900 dpi included, which is a
                real number of pixels and is refused politely when it would be
                larger than the browser can hold
     data       the alignment, in the formats Block 2 reads back
     everything ZIP, written here rather than fetched, so the app keeps working
                with no network and from a double-clicked file

   Wrapped in a named function so a Web Worker can be built from its own source
   text (js/pool.js). The only DOM it touches is the canvas the PNG needs, and
   that is behind a check. */

function ExportCore(g) {
const Export = {};

  const esc = s => String(s == null ? '' : s)
    .replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;')
    .replace(/"/g, '&quot;');

  /* ================================================================
     1 · trees
     ================================================================ */
  function newick(tree, opts) {
    return g.Tree.writeNewick(tree, opts || {});
  }

  /* NEXUS with a TRANSLATE block: the names are written once and the trees
     refer to them by number, which is what keeps a file of a thousand trees
     from being mostly names. FigTree reads the extra block and nothing else
     minds it. */
  function nexusTrees(entries, opts) {
    opts = opts || {};
    const labels = opts.labels || [];
    const taxa = opts.taxa || labels;
    let out = '#NEXUS\n\n[ Written by PhylogenyPro ]\n\n';
    out += 'BEGIN TAXA;\n';
    out += `  DIMENSIONS NTAX=${taxa.length};\n`;
    out += '  TAXLABELS\n';
    taxa.forEach(t => { out += `    ${quoteNexus(t)}\n`; });
    out += '  ;\nEND;\n\nBEGIN TREES;\n';
    if (opts.translate !== false) {
      out += '  TRANSLATE\n';
      taxa.forEach((t, i) => {
        out += `    ${i + 1} ${quoteNexus(t)}${i < taxa.length - 1 ? ',' : ''}\n`;
      });
      out += '  ;\n';
    }
    const numberLabels = taxa.map((t, i) => String(i + 1));
    entries.forEach((e, i) => {
      const nm = (e.name || ('tree' + (i + 1))).replace(/[\s(),:;[\]']/g, '_');
      const nwk = g.Tree.writeNewick(e.tree, Object.assign({}, opts.newick || {},
        { labels: opts.translate === false ? labels : numberLabels }));
      out += `  TREE ${nm} = ${e.rooted === false ? '[&U] ' : '[&R] '}${nwk}\n`;
    });
    out += 'END;\n';
    if (opts.figtree) {
      out += '\nbegin figtree;\n';
      out += '  set appearance.branchLineWidth=1.5;\n';
      out += '  set nodeLabels.isShown=true;\n';
      out += '  set nodeLabels.displayAttribute="label";\n';
      out += '  set tipLabels.fontSize=11;\n';
      out += '  set scaleBar.isShown=true;\n';
      out += 'end;\n';
    }
    return out;
  }
  const quoteNexus = s => (/[\s(),:;[\]'=]/.test(String(s)) ? "'" + String(s).replace(/'/g, "''") + "'" : String(s));

  /* phyloXML (Han & Zmasek 2009). Verbose, and the only common format that can
     carry the support, the age and the confidence interval of a node together
     with its name, which is exactly what Blocks 5, 6 and 7 produce and what
     Newick has to smuggle in comments. */
  function phyloXML(entries, opts) {
    opts = opts || {};
    const labels = opts.labels || [];
    let out = '<?xml version="1.0" encoding="UTF-8"?>\n';
    out += '<phyloxml xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"';
    out += ' xsi:schemaLocation="http://www.phyloxml.org http://www.phyloxml.org/1.10/phyloxml.xsd"';
    out += ' xmlns="http://www.phyloxml.org">\n';
    entries.forEach(e => {
      out += `  <phylogeny rooted="${e.rooted === false ? 'false' : 'true'}">\n`;
      if (e.name) out += `    <name>${esc(e.name)}</name>\n`;
      if (e.description) out += `    <description>${esc(e.description)}</description>\n`;
      out += node(e.tree, 2, e);
      out += '  </phylogeny>\n';
    });
    out += '</phyloxml>\n';
    return out;

    function node(nd, depth, entry) {
      const pad = '  '.repeat(depth + 1);
      let s = `${pad}<clade>\n`;
      if (nd.len != null) s += `${pad}  <branch_length>${(+nd.len).toPrecision(10)}</branch_length>\n`;
      const nm = nd.tip != null ? (labels[nd.tip] != null ? labels[nd.tip] : nd.label) : nd.label;
      if (nm) s += `${pad}  <name>${esc(nm)}</name>\n`;
      if (nd.support != null) {
        const type = (entry && entry.supportType) || (nd.support > 1 ? 'bootstrap' : 'posterior probability');
        s += `${pad}  <confidence type="${esc(type)}">${+nd.support}</confidence>\n`;
      }
      if (nd.age != null) {
        s += `${pad}  <date unit="${esc((entry && entry.unit) || 'Ma')}">\n`;
        s += `${pad}    <value>${+nd.age}</value>\n`;
        if (nd.ageLow != null) s += `${pad}    <minimum>${+nd.ageLow}</minimum>\n`;
        if (nd.ageHigh != null) s += `${pad}    <maximum>${+nd.ageHigh}</maximum>\n`;
        s += `${pad}  </date>\n`;
      }
      nd.children.forEach(c => {
        const child = Object.assign({}, c.node, { len: c.len });
        s += node(child, depth + 1, entry);
      });
      s += `${pad}</clade>\n`;
      return s;
    }
  }
  /* reading it back, which is how the writer gets tested */
  function readPhyloXML(text, labels) {
    if (typeof DOMParser === 'undefined') return null;
    const doc = new DOMParser().parseFromString(text, 'application/xml');
    if (doc.getElementsByTagName('parsererror').length) return null;
    const byName = labels ? new Map(labels.map((l, i) => [l, i])) : null;
    const out = [];
    const phys = doc.getElementsByTagName('phylogeny');
    for (let i = 0; i < phys.length; i++) {
      const top = firstChild(phys[i], 'clade');
      if (!top) continue;
      out.push({
        name: textOf(phys[i], 'name'),
        rooted: phys[i].getAttribute('rooted') !== 'false',
        tree: build(top),
      });
    }
    return out;

    function firstChild(el, tag) {
      for (let i = 0; i < el.children.length; i++) if (el.children[i].localName === tag) return el.children[i];
      return null;
    }
    function childrenNamed(el, tag) {
      const r = [];
      for (let i = 0; i < el.children.length; i++) if (el.children[i].localName === tag) r.push(el.children[i]);
      return r;
    }
    function textOf(el, tag) {
      const c = firstChild(el, tag);
      return c ? c.textContent.trim() : null;
    }
    function build(cl) {
      const kids = childrenNamed(cl, 'clade');
      const nm = textOf(cl, 'name');
      const nd = { label: nm, children: [] };
      const conf = firstChild(cl, 'confidence');
      if (conf) nd.support = +conf.textContent;
      const date = firstChild(cl, 'date');
      if (date) {
        const v = textOf(date, 'value');
        if (v != null) nd.age = +v;
      }
      if (!kids.length) {
        nd.tip = byName && byName.has(nm) ? byName.get(nm) : out.length * 0;
        if (byName && !byName.has(nm)) nd.tip = -1;
        return nd;
      }
      kids.forEach(k => {
        const child = build(k);
        const bl = textOf(k, 'branch_length');
        nd.children.push({ node: child, len: bl == null ? 0 : +bl });
      });
      return nd;
    }
  }

  /* ================================================================
     2 · figures
     ================================================================ */
  /* An SVG is resolution-free; a journal asking for 900 dpi is asking for a
     particular number of pixels. This turns one into the other, refuses the
     sizes a browser canvas cannot hold instead of returning a blank image, and
     paints the background rather than leaving it transparent, because a
     transparent PNG printed on paper is a white rectangle with black text only
     by luck. */
  const MAX_PIXELS = 268435456;                 // 16384², which browsers allow
  function pngSize(svgText, opts) {
    opts = opts || {};
    const vb = /viewBox\s*=\s*"([^"]+)"/.exec(svgText);
    let w = 800, h = 600;
    if (vb) {
      const p = vb[1].trim().split(/[\s,]+/).map(Number);
      if (p.length === 4 && p[2] > 0 && p[3] > 0) { w = p[2]; h = p[3]; }
    }
    /* the SVG's own units are treated as points at 96 dpi, which is what a
       browser means by a pixel */
    const dpi = opts.dpi || 300;
    const scale = opts.scale || (dpi / 96);
    const W = Math.round(w * scale), H = Math.round(h * scale);
    return { width: W, height: H, baseWidth: w, baseHeight: h, scale, dpi,
      pixels: W * H, tooBig: W * H > MAX_PIXELS, maxPixels: MAX_PIXELS };
  }
  function svgToPng(svgText, opts) {
    opts = opts || {};
    const size = pngSize(svgText, opts);
    if (size.tooBig) {
      return Promise.reject(new Error(`too_big:${size.width}x${size.height}`));
    }
    if (typeof document === 'undefined') return Promise.reject(new Error('no_dom'));
    /* the colours are CSS variables, so they have to be resolved before the
       image leaves the page: an SVG loaded into an <img> has no page to ask */
    const resolved = inlineVars(svgText, opts.vars || readVars());
    const withSize = resolved.replace(/<svg/, `<svg width="${size.width}" height="${size.height}"`);
    const blob = new Blob([withSize], { type: 'image/svg+xml;charset=utf-8' });
    const url = URL.createObjectURL(blob);
    return new Promise((resolve, reject) => {
      const img = new Image();
      img.onload = () => {
        try {
          const canvas = document.createElement('canvas');
          canvas.width = size.width; canvas.height = size.height;
          const ctx = canvas.getContext('2d');
          ctx.fillStyle = opts.background || '#ffffff';
          ctx.fillRect(0, 0, size.width, size.height);
          ctx.drawImage(img, 0, 0, size.width, size.height);
          canvas.toBlob(b => {
            URL.revokeObjectURL(url);
            b ? resolve({ blob: b, size }) : reject(new Error('canvas_failed'));
          }, 'image/png');
        } catch (e) { URL.revokeObjectURL(url); reject(e); }
      };
      img.onerror = () => { URL.revokeObjectURL(url); reject(new Error('svg_load_failed')); };
      img.src = url;
    });
  }
  /* every CSS custom property the figures use, read once from the live page */
  function readVars() {
    const out = {};
    if (typeof document === 'undefined') return out;
    const cs = getComputedStyle(document.documentElement);
    const names = ['bg', 'bg-soft', 'card-bg', 'text', 'ink', 'text-muted', 'border', 'border-strong',
      'primary', 'primary-dark', 'accent', 'gold', 'leaf', 'sky', 'rose', 'success', 'danger', 'warning',
      'heat-lo', 'heat-mid', 'heat-hi', 'base-a', 'base-c', 'base-g', 'base-t', 'base-gap'];
    for (let i = 1; i <= 10; i++) names.push('c' + i);
    names.forEach(n => {
      const v = cs.getPropertyValue('--' + n).trim();
      if (v) out[n] = v;
    });
    return out;
  }
  function inlineVars(svgText, vars) {
    let out = svgText.replace(/var\(--([a-z0-9-]+)\)/gi, (m, name) => (vars[name] || '#444444'));
    /* the two classes the figures use for text come from the stylesheet, and an
       <img> has no stylesheet */
    if (out.indexOf('<style') < 0) {
      const style = `<style>.art-txt{fill:${vars.text || '#161b2a'};font-family:system-ui,sans-serif}` +
        `.art-mut{fill:${vars['text-muted'] || '#5a6479'};font-family:system-ui,sans-serif}</style>`;
      out = out.replace(/(<svg[^>]*>)/, '$1' + style);
    }
    return out;
  }

  /* ================================================================
     3 · ZIP
     ================================================================ */
  const CRC_TABLE = (() => {
    const t = new Uint32Array(256);
    for (let n = 0; n < 256; n++) {
      let c = n;
      for (let k = 0; k < 8; k++) c = (c & 1) ? (0xEDB88320 ^ (c >>> 1)) : (c >>> 1);
      t[n] = c >>> 0;
    }
    return t;
  })();
  function crc32(bytes) {
    let c = 0xFFFFFFFF;
    for (let i = 0; i < bytes.length; i++) c = CRC_TABLE[(c ^ bytes[i]) & 0xFF] ^ (c >>> 8);
    return (c ^ 0xFFFFFFFF) >>> 0;
  }
  const utf8 = s => new TextEncoder().encode(s);

  /* A ZIP written by hand: no library, no network, works from a file:// page.
     Entries are stored uncompressed unless the browser offers a deflate stream,
     which most do; a reader cannot tell the difference. */
  async function zip(files, opts) {
    opts = opts || {};
    const canDeflate = opts.compress !== false && typeof CompressionStream === 'function';
    const now = opts.date || new Date();
    const dosTime = ((now.getHours() & 31) << 11) | ((now.getMinutes() & 63) << 5) | ((now.getSeconds() / 2) & 31);
    const dosDate = (((now.getFullYear() - 1980) & 127) << 9) | (((now.getMonth() + 1) & 15) << 5) | (now.getDate() & 31);

    const parts = [], central = [];
    let offset = 0;
    for (let i = 0; i < files.length; i++) {
      const f = files[i];
      const nameBytes = utf8(f.name);
      const raw = typeof f.content === 'string' ? utf8(f.content)
        : (f.content instanceof Uint8Array ? f.content : new Uint8Array(await f.content.arrayBuffer()));
      const crc = crc32(raw);
      let body = raw, method = 0;
      if (canDeflate && raw.length > 64) {
        try {
          const cs = new CompressionStream('deflate-raw');
          const stream = new Blob([raw]).stream().pipeThrough(cs);
          const packed = new Uint8Array(await new Response(stream).arrayBuffer());
          if (packed.length < raw.length) { body = packed; method = 8; }
        } catch (e) { /* stored, which is always valid */ }
      }
      const local = new Uint8Array(30 + nameBytes.length);
      const dv = new DataView(local.buffer);
      dv.setUint32(0, 0x04034b50, true);
      dv.setUint16(4, 20, true);
      dv.setUint16(6, 0x0800, true);              // UTF-8 names
      dv.setUint16(8, method, true);
      dv.setUint16(10, dosTime, true);
      dv.setUint16(12, dosDate, true);
      dv.setUint32(14, crc, true);
      dv.setUint32(18, body.length, true);
      dv.setUint32(22, raw.length, true);
      dv.setUint16(26, nameBytes.length, true);
      dv.setUint16(28, 0, true);
      local.set(nameBytes, 30);
      parts.push(local, body);

      const cen = new Uint8Array(46 + nameBytes.length);
      const cv = new DataView(cen.buffer);
      cv.setUint32(0, 0x02014b50, true);
      cv.setUint16(4, 20, true);
      cv.setUint16(6, 20, true);
      cv.setUint16(8, 0x0800, true);
      cv.setUint16(10, method, true);
      cv.setUint16(12, dosTime, true);
      cv.setUint16(14, dosDate, true);
      cv.setUint32(16, crc, true);
      cv.setUint32(20, body.length, true);
      cv.setUint32(24, raw.length, true);
      cv.setUint16(28, nameBytes.length, true);
      cv.setUint32(42, offset, true);
      cen.set(nameBytes, 46);
      central.push(cen);
      offset += local.length + body.length;
    }
    let centralSize = 0;
    central.forEach(c => { centralSize += c.length; });
    const end = new Uint8Array(22);
    const ev = new DataView(end.buffer);
    ev.setUint32(0, 0x06054b50, true);
    ev.setUint16(8, files.length, true);
    ev.setUint16(10, files.length, true);
    ev.setUint32(12, centralSize, true);
    ev.setUint32(16, offset, true);
    return new Blob(parts.concat(central, [end]), { type: 'application/zip' });
  }

  Object.assign(Export, {
    newick, nexusTrees, phyloXML, readPhyloXML, quoteNexus,
    pngSize, svgToPng, readVars, inlineVars, MAX_PIXELS,
    crc32, zip,
  });
  g.Export = Export;
}
ExportCore(typeof window !== 'undefined' ? window : self);
