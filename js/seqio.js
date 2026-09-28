/* PhylogenyPro — molecular phylogenetics in the browser.
   Copyright (C) 2026  Luis Ángel Barrera-Guzmán

   This program is free software: you can redistribute it and/or modify it under
   the terms of the GNU General Public License as published by the Free Software
   Foundation, either version 3 of the License, or (at your option) any later
   version. It is distributed in the hope that it will be useful, but WITHOUT ANY
   WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
   A PARTICULAR PURPOSE. See the GNU General Public License, in the file LICENSE
   at the root of this program, or <https://www.gnu.org/licenses/>. */

/* PhylogenyPro — readers and writers of sequence formats.

   One entry point, SeqIO.parse(text, name), which recognises the format and
   always returns the same shape:

     { format, taxa:[names], seqs:[strings], type:'dna'|'aa'|'morph',
       charsets:[{name, from, to, every?}], trees:[{name, newick}],
       codes:{missing, gap, matchchar}, warnings:[] }

   Sequences come back upper-cased, with the gap written '-' and missing data
   written '?' (or 'N'/'X' where the file used them), and every sequence of an
   aligned file padded to the same length. Nothing here assumes the file was
   written by a particular program: the formats in the wild are far looser than
   their specifications, so each reader is deliberately forgiving and reports
   what it had to guess in `warnings`. */

const SeqIO = {};

(function () {

  const DNA_CHARS = /^[ACGTUMRWSYKVHDBNX\-?.]+$/i;
  const AA_CHARS = /^[ACDEFGHIKLMNPQRSTVWYBZXJUO*\-?.]+$/i;

  /* ---------------- helpers ---------------- */
  const clean = s => String(s == null ? '' : s).replace(/\r\n?/g, '\n');
  const isBlank = l => /^\s*$/.test(l);
  /* a NEXUS/PHYLIP token: quoted, or a run of non-space characters */
  function splitNameSeq(line) {
    const m = line.match(/^\s*(?:'([^']*)'|"([^"]*)"|(\S+))\s+(.*)$/);
    if (!m) return null;
    return { name: (m[1] != null ? m[1] : m[2] != null ? m[2] : m[3]), rest: m[4] };
  }
  const stripSpaces = s => s.replace(/[\s\d]/g, '');

  /* ---------------- type guessing ----------------
     Counting letters is not enough: a protein alignment of mostly alanine and
     glycine would look like DNA. The rule used here is the practical one: if
     the residues outside the nucleotide alphabet are more than 1% of the
     unambiguous characters, it is a protein. */
  function guessType(seqs) {
    let dnaLike = 0, aaOnly = 0, digits = 0, total = 0;
    seqs.forEach(s => {
      for (let i = 0; i < s.length; i++) {
        const c = s[i].toUpperCase();
        if (c === '-' || c === '?' || c === '.' || c === ' ') continue;
        total++;
        if ('ACGTUN'.indexOf(c) >= 0) dnaLike++;
        else if ('EFILPQZJO*'.indexOf(c) >= 0) aaOnly++;
        else if (c >= '0' && c <= '9') digits++;
      }
    });
    if (!total) return 'dna';
    if (digits / total > 0.5) return 'morph';
    if (aaOnly / total > 0.01) return 'aa';
    if (dnaLike / total > 0.85) return 'dna';
    return 'aa';
  }

  /* ---------------- format detection ---------------- */
  function detect(text, name) {
    /* leading blank lines and a byte-order mark are common in files that went
       through a spreadsheet or a copy-paste, and they must not hide the format */
    const t = clean(text).replace(/^﻿/, '').replace(/^\s*\n/, '').replace(/^[ \t]+/, '');
    const head = t.slice(0, 4000);
    const ext = (String(name || '').match(/\.([a-z0-9]+)$/i) || [, ''])[1].toLowerCase();
    if (/^\s*#NEXUS/i.test(head)) return 'nexus';
    if (/^\s*CLUSTAL/i.test(head)) return 'clustal';
    if (/^\s*#MEGA/i.test(head)) return 'mega';
    if (/^\s*LOCUS\s/i.test(head) || /^\s*ORIGIN\s*$/im.test(head) && /^\s*\/\//m.test(head)) return 'genbank';
    if (/^\s*>/.test(head)) return 'fasta';
    if (/^\s*\d+\s+\d+\s*$/m.test(head.split('\n')[0] || '')) return 'phylip';
    if (/^\s*\(.*\)\s*;\s*$/s.test(t.trim())) return 'newick';
    const lines = t.split('\n').filter(l => !isBlank(l)).slice(0, 5);
    const first = lines[0] || '';
    if ((ext === 'csv' || ext === 'tsv' || ext === 'txt') && /[,;\t]/.test(first)) return 'table';
    /* a table with no helpful extension: the same number of separators on every
       line, and more than one column */
    if (lines.length > 1) {
      const sepCount = l => Math.max((l.match(/,/g) || []).length, (l.match(/\t/g) || []).length, (l.match(/;/g) || []).length);
      const k = sepCount(first);
      if (k >= 2 && lines.every(l => sepCount(l) === k)) return 'table';
    }
    return 'unknown';
  }

  /* ================================================================
     FASTA
     ================================================================ */
  function parseFasta(text) {
    const out = { taxa: [], seqs: [], warnings: [] };
    let cur = null, buf = [];
    clean(text).split('\n').forEach(line => {
      if (line[0] === '>' || line[0] === ';') {
        if (line[0] === ';') return;                         // old comment style
        if (cur != null) { out.taxa.push(cur); out.seqs.push(buf.join('')); }
        cur = line.slice(1).trim();
        buf = [];
      } else if (cur != null) {
        buf.push(line.replace(/\s/g, ''));
      }
    });
    if (cur != null) { out.taxa.push(cur); out.seqs.push(buf.join('')); }
    return out;
  }

  /* ================================================================
     PHYLIP — sequential and interleaved, strict and relaxed
     ================================================================ */
  function parsePhylip(text) {
    const lines = clean(text).split('\n');
    let i = 0;
    while (i < lines.length && isBlank(lines[i])) i++;
    const dim = (lines[i] || '').match(/^\s*(\d+)\s+(\d+)/);
    if (!dim) throw new Error('PHYLIP: the first line must give the number of taxa and of sites.');
    const ntax = +dim[1], nchar = +dim[2];
    i++;
    const warnings = [];
    const body = lines.slice(i);

    /* PHYLIP comes in two flavours and files rarely say which: strict puts the
       name in the first ten columns, relaxed separates it with white space. The
       file is read both ways and the reading that matches the announced number
       of sites wins; if neither matches exactly, the relaxed one is kept,
       because that is what every modern program writes. */
    function readWith(strict) {
      const taxa = [], seqs = [];
      let seen = 0, k = 0, j = 0;
      while (j < body.length && seen < ntax) {
        const line = body[j++];
        if (isBlank(line)) continue;
        let name, rest;
        if (strict) { name = line.slice(0, 10).trim(); rest = line.slice(10); }
        else {
          const sp = splitNameSeq(line);
          if (sp) { name = sp.name; rest = sp.rest; }
          else { name = line.trim(); rest = ''; }
        }
        taxa.push(name); seqs.push(stripSpaces(rest)); seen++;
      }
      while (j < body.length) {                       // interleaved blocks: sequence only
        const line = body[j++];
        if (isBlank(line)) { k = 0; continue; }
        if (k >= ntax) k = 0;
        seqs[k] = (seqs[k] || '') + stripSpaces(line);
        k++;
      }
      return { taxa, seqs };
    }
    const relaxed = readWith(false), strict = readWith(true);
    const fits = r => r.seqs.length === ntax && r.seqs.every(s => s.length === nchar);
    const chosen = fits(relaxed) ? relaxed : (fits(strict) ? strict : relaxed);
    const taxa = chosen.taxa, seqs = chosen.seqs;
    if (!fits(relaxed) && fits(strict)) warnings.push({
      es: 'El archivo está en PHYLIP estricto: los nombres ocupan las diez primeras columnas.',
      en: 'The file is in strict PHYLIP: names occupy the first ten columns.' });
    if (seqs.some(s => s.length !== nchar)) {
      warnings.push({ es: `El encabezado anuncia ${nchar} sitios y no todas las secuencias los tienen; se usó lo que traía el archivo.`,
        en: `The header announces ${nchar} sites and not every sequence has them; what the file actually contained was used.` });
    }
    return { taxa, seqs, warnings };
  }

  /* ================================================================
     NEXUS — DATA/CHARACTERS, SETS (charset) and TREES
     ================================================================ */
  function stripNexusComments(t) {
    /* comments are in square brackets and can nest; anything inside is dropped,
       except the [&...] annotations inside a tree string, which are kept */
    let out = '', depth = 0;
    for (let i = 0; i < t.length; i++) {
      const c = t[i];
      if (c === '[') { depth++; continue; }
      if (c === ']') { if (depth > 0) depth--; continue; }
      if (!depth) out += c;
    }
    return out;
  }
  function parseNexus(text) {
    const raw = clean(text);
    const t = stripNexusComments(raw);
    const warnings = [], charsets = [], trees = [];
    let taxa = [], seqs = [], type = null;
    const codes = { missing: '?', gap: '-', matchchar: null };
    let ntax = null, nchar = null, interleave = false;

    /* --- blocks --- */
    const blockRe = /BEGIN\s+(\w+)\s*;([\s\S]*?)(?:^|\n)\s*END(?:BLOCK)?\s*;/gi;
    let m;
    while ((m = blockRe.exec(t)) !== null) {
      const kind = m[1].toUpperCase(), body = m[2];
      if (kind === 'DATA' || kind === 'CHARACTERS') {
        const dim = body.match(/DIMENSIONS([^;]*);/i);
        if (dim) {
          const nt = dim[1].match(/NTAX\s*=\s*(\d+)/i), nc = dim[1].match(/NCHAR\s*=\s*(\d+)/i);
          if (nt) ntax = +nt[1];
          if (nc) nchar = +nc[1];
        }
        const fmt = body.match(/FORMAT([^;]*);/i);
        if (fmt) {
          const f = fmt[1];
          const dt = f.match(/DATATYPE\s*=\s*(\w+)/i);
          if (dt) {
            const d = dt[1].toLowerCase();
            type = (d === 'protein' || d === 'aa') ? 'aa' : (d === 'standard' || d === 'restriction') ? 'morph' : 'dna';
          }
          const mi = f.match(/MISSING\s*=\s*(\S)/i); if (mi) codes.missing = mi[1];
          const ga = f.match(/GAP\s*=\s*(\S)/i); if (ga) codes.gap = ga[1];
          const mc = f.match(/MATCHCHAR\s*=\s*(\S)/i); if (mc) codes.matchchar = mc[1];
          if (/INTERLEAVE\s*(=\s*YES)?/i.test(f) && !/INTERLEAVE\s*=\s*NO/i.test(f)) interleave = true;
        }
        const mx = body.match(/MATRIX([\s\S]*?);\s*$/i) || body.match(/MATRIX([\s\S]*)/i);
        if (mx) {
          const rows = new Map(); const order = [];
          mx[1].split('\n').forEach(line => {
            if (isBlank(line)) return;
            const sp = splitNameSeq(line.replace(/;\s*$/, ''));
            if (!sp) return;
            const seq = stripSpaces(sp.rest);
            if (!seq) return;
            if (!rows.has(sp.name)) { rows.set(sp.name, ''); order.push(sp.name); }
            rows.set(sp.name, rows.get(sp.name) + seq);
          });
          taxa = order; seqs = order.map(k2 => rows.get(k2));
        }
      } else if (kind === 'SETS') {
        const re = /CHARSET\s+([^\s=]+)\s*=\s*([^;]+);/gi;
        let c;
        while ((c = re.exec(body)) !== null) {
          const name = c[1].replace(/['"]/g, '');
          c[2].trim().split(/\s+/).forEach(rangeTxt => {
            const r = rangeTxt.match(/^(\d+)\s*-\s*(\d+|\.)(?:\\(\d+))?$/) || rangeTxt.match(/^(\d+)$/);
            if (!r) return;
            const from = +r[1];
            const to = r[2] == null ? from : (r[2] === '.' ? (nchar || from) : +r[2]);
            charsets.push({ name, from, to, every: r[3] ? +r[3] : null });
          });
        }
      } else if (kind === 'TREES') {
        const tr = /TREE\s+\*?\s*([^\s=]+)\s*=\s*(?:\[[^\]]*\]\s*)?([^;]+);/gi;
        let c;
        const transl = {};
        const tb = body.match(/TRANSLATE([\s\S]*?);/i);
        if (tb) tb[1].split(',').forEach(p => {
          const q = p.trim().match(/^(\S+)\s+(.+)$/);
          if (q) transl[q[1]] = q[2].replace(/['"]/g, '');
        });
        while ((c = tr.exec(body)) !== null) {
          let nwk = c[2].trim();
          if (Object.keys(transl).length) nwk = nwk.replace(/([(,])\s*(\d+)/g, (s2, pre, num) => transl[num] ? pre + transl[num] : s2);
          trees.push({ name: c[1], newick: nwk + ';' });
        }
      }
    }
    if (!taxa.length) throw new Error('NEXUS: no MATRIX was found in a DATA or CHARACTERS block.');
    if (ntax && taxa.length !== ntax) warnings.push({
      es: `El bloque anuncia ${ntax} taxones y se leyeron ${taxa.length}.`,
      en: `The block announces ${ntax} taxa and ${taxa.length} were read.` });
    if (interleave) { /* already handled: rows were concatenated by name */ }
    return { taxa, seqs, type, charsets, trees, codes, warnings };
  }

  /* ================================================================
     CLUSTAL
     ================================================================ */
  function parseClustal(text) {
    const lines = clean(text).split('\n');
    const rows = new Map(), order = [];
    lines.slice(1).forEach(line => {
      if (isBlank(line) || /^\s/.test(line)) return;         // conservation line starts with spaces
      const sp = splitNameSeq(line);
      if (!sp) return;
      const seq = stripSpaces(sp.rest);
      if (!seq) return;
      if (!rows.has(sp.name)) { rows.set(sp.name, ''); order.push(sp.name); }
      rows.set(sp.name, rows.get(sp.name) + seq);
    });
    return { taxa: order, seqs: order.map(k => rows.get(k)), warnings: [] };
  }

  /* ================================================================
     MEGA (.meg)
     ================================================================ */
  function parseMega(text) {
    const lines = clean(text).split('\n');
    const rows = new Map(), order = [];
    let type = null, cur = null;
    lines.forEach(line => {
      const l = line.trim();
      if (!l || /^#MEGA/i.test(l)) return;
      if (/^!/.test(l)) {
        const dt = l.match(/DataType\s*=\s*(\w+)/i);
        if (dt) type = /prot/i.test(dt[1]) ? 'aa' : 'dna';
        return;
      }
      if (l[0] === '#') {
        const m = l.slice(1).match(/^(\S+)\s*(.*)$/);
        cur = m ? m[1] : l.slice(1).trim();
        if (!rows.has(cur)) { rows.set(cur, ''); order.push(cur); }
        if (m && m[2]) rows.set(cur, rows.get(cur) + stripSpaces(m[2]));
        return;
      }
      if (cur) rows.set(cur, rows.get(cur) + stripSpaces(l));
    });
    return { taxa: order, seqs: order.map(k => rows.get(k)), type, warnings: [] };
  }

  /* ================================================================
     GenBank flat file — several records in one file
     ================================================================ */
  function parseGenbank(text) {
    const recs = clean(text).split(/^\/\/\s*$/m);
    const taxa = [], seqs = [], warnings = [];
    recs.forEach(rec => {
      if (!/LOCUS/.test(rec)) return;
      const acc = (rec.match(/^ACCESSION\s+(\S+)/m) || [])[1] || (rec.match(/^LOCUS\s+(\S+)/m) || [])[1] || 'seq';
      const org = (rec.match(/^\s+ORGANISM\s+(.+)$/m) || [])[1];
      const gene = (rec.match(/\/gene="([^"]+)"/) || [])[1];
      const ori = rec.split(/^ORIGIN\s*$/m)[1];
      if (!ori) return;
      const seq = ori.replace(/[\s\d\/]/g, '');
      let name = org ? org.trim().replace(/\s+/g, '_') : acc;
      if (gene) name += '_' + gene;
      if (taxa.indexOf(name) >= 0) name = name + '_' + acc;
      taxa.push(name); seqs.push(seq);
    });
    if (!taxa.length) throw new Error('GenBank: no ORIGIN block with sequence was found.');
    warnings.push({ es: 'Los nombres se tomaron del campo ORGANISM (y del gen, si estaba); revísalos antes de seguir.',
      en: 'Names were taken from the ORGANISM field (and the gene, when present); check them before going on.' });
    return { taxa, seqs, warnings };
  }

  /* ================================================================
     Table: a morphological matrix, one row per taxon
     ================================================================ */
  function parseTable(text) {
    const lines = clean(text).split('\n').filter(l => !isBlank(l));
    if (!lines.length) throw new Error('The table is empty.');
    const sep = (lines[0].match(/\t/g) || []).length >= (lines[0].match(/[,;]/g) || []).length ? '\t' : (lines[0].indexOf(';') >= 0 && lines[0].indexOf(',') < 0 ? ';' : ',');
    const cells = lines.map(l => l.split(sep).map(c => c.trim().replace(/^"|"$/g, '')));
    const header = cells[0];
    /* a header row is one whose cells are not single-character states */
    const looksHeader = header.slice(1).some(c => c.length > 1 || /[A-Za-z]{2}/.test(c));
    const body = looksHeader ? cells.slice(1) : cells;
    const charNames = looksHeader ? header.slice(1) : header.slice(1).map((_, i) => 'car' + (i + 1));
    const taxa = body.map(r => r[0]);
    const seqs = body.map(r => r.slice(1).map(v => (v === '' || v === 'NA' || v === '-9') ? '?' : v).join(''));
    const warnings = [];
    if (body.some(r => r.length !== body[0].length)) warnings.push({
      es: 'No todas las filas tienen el mismo número de columnas.', en: 'Not every row has the same number of columns.' });
    if (seqs.some(s => s.length !== charNames.length)) warnings.push({
      es: 'Algún estado ocupa más de un carácter; los estados con dos dígitos deben ir entre paréntesis.',
      en: 'Some state takes more than one character; two-digit states must be written in brackets.' });
    return { taxa, seqs, type: 'morph', charNames, warnings };
  }

  /* ================================================================
     the entry point
     ================================================================ */
  function parse(text, name) {
    const format = detect(text, name);
    let r;
    switch (format) {
      case 'fasta': r = parseFasta(text); break;
      case 'phylip': r = parsePhylip(text); break;
      case 'nexus': r = parseNexus(text); break;
      case 'clustal': r = parseClustal(text); break;
      case 'mega': r = parseMega(text); break;
      case 'genbank': r = parseGenbank(text); break;
      case 'table': r = parseTable(text); break;
      case 'newick': return { format, taxa: [], seqs: [], trees: [{ name: 'tree', newick: clean(text).trim() }], warnings: [] };
      default: throw new Error('No se reconoció el formato del archivo. / The file format was not recognised.');
    }
    const out = Object.assign({ format, charsets: [], trees: [], codes: { missing: '?', gap: '-', matchchar: null }, warnings: [] }, r);
    if (!out.taxa.length) throw new Error('No sequences were found in the file.');

    /* normalise: upper case, one gap symbol, one missing symbol */
    const gap = out.codes.gap || '-', miss = out.codes.missing || '?';
    out.seqs = out.seqs.map(s => {
      let x = String(s || '').toUpperCase();
      if (gap !== '-') x = x.split(gap).join('-');
      if (miss !== '?' && miss !== 'N' && miss !== 'X') x = x.split(miss.toUpperCase()).join('?');
      return x;
    });
    /* MATCHCHAR: a dot means "the same as the first sequence" */
    const mc = out.codes.matchchar;
    if (mc && out.seqs.length) {
      const ref = out.seqs[0];
      out.seqs = out.seqs.map((s, i) => i === 0 ? s : s.split('').map((c, j) => c === mc ? (ref[j] || '?') : c).join(''));
    } else if (out.format !== 'table' && out.seqs.some(s => s.indexOf('.') >= 0)) {
      const ref = out.seqs[0];
      out.seqs = out.seqs.map((s, i) => i === 0 ? s : s.split('').map((c, j) => c === '.' ? (ref[j] || '?') : c).join(''));
      out.warnings.push({ es: 'Los puntos se leyeron como "igual que la primera secuencia".',
        en: 'Dots were read as "same as the first sequence".' });
    }
    if (!out.type) out.type = guessType(out.seqs);

    /* duplicated names get a suffix, or two different sequences would be merged */
    const seen = {};
    out.taxa = out.taxa.map(nm => {
      let base = String(nm).trim() || 'sin_nombre';
      if (seen[base] == null) { seen[base] = 0; return base; }
      seen[base]++;
      out.warnings.push({ es: `El nombre «${base}» venía repetido; el segundo se renombró.`,
        en: `The name "${base}" was repeated; the second one was renamed.` });
      return base + '_' + (seen[base] + 1);
    });

    /* aligned? every sequence of the same length */
    const lens = out.seqs.map(s => s.length);
    out.aligned = lens.every(l => l === lens[0]);
    out.nchar = Math.max(...lens);
    if (!out.aligned && (format === 'phylip' || format === 'nexus' || format === 'clustal')) {
      out.warnings.push({ es: 'El archivo dice estar alineado pero las secuencias no miden lo mismo; se rellenaron con huecos al final.',
        en: 'The file claims to be aligned but the sequences differ in length; they were padded with gaps at the end.' });
      out.seqs = out.seqs.map(s => s + '-'.repeat(out.nchar - s.length));
      out.aligned = true;
    }
    return out;
  }

  /* ================================================================
     writers
     ================================================================ */
  function writeFasta(taxa, seqs, width) {
    width = width || 60;
    return taxa.map((t, i) => {
      const s = seqs[i] || '';
      const lines = [];
      for (let k = 0; k < s.length; k += width) lines.push(s.slice(k, k + width));
      return '>' + t + '\n' + (lines.join('\n') || '');
    }).join('\n') + '\n';
  }
  /* PHYLIP: relaxed by default (names are not cut at ten characters, which is
     what every modern program reads); strict on request */
  function writePhylip(taxa, seqs, opts) {
    opts = opts || {};
    const n = taxa.length, L = seqs[0] ? seqs[0].length : 0;
    const names = opts.strict
      ? taxa.map(t => (t.slice(0, 10) + '          ').slice(0, 10))
      : (() => { const w = Math.max(...taxa.map(t => t.length)) + 2; return taxa.map(t => t + ' '.repeat(w - t.length)); })();
    let out = ' ' + n + ' ' + L + '\n';
    if (!opts.interleaved) {
      taxa.forEach((t, i) => { out += names[i] + seqs[i] + '\n'; });
      return out;
    }
    const block = opts.block || 60;
    for (let start = 0; start < L; start += block) {
      taxa.forEach((t, i) => {
        out += (start === 0 ? names[i] : '') + seqs[i].slice(start, start + block) + '\n';
      });
      out += '\n';
    }
    return out;
  }
  function writeNexus(taxa, seqs, opts) {
    opts = opts || {};
    const type = opts.type === 'aa' ? 'PROTEIN' : opts.type === 'morph' ? 'STANDARD' : 'DNA';
    const L = seqs[0] ? seqs[0].length : 0;
    const w = Math.max(...taxa.map(t => t.length)) + 2;
    let out = '#NEXUS\n\n[ Written by PhylogenyPro ]\n\nBEGIN DATA;\n';
    out += `  DIMENSIONS NTAX=${taxa.length} NCHAR=${L};\n`;
    out += `  FORMAT DATATYPE=${type} MISSING=? GAP=-${type === 'STANDARD' ? ' SYMBOLS="0123456789"' : ''} INTERLEAVE=NO;\n  MATRIX\n`;
    taxa.forEach((t, i) => { out += '    ' + t + ' '.repeat(w - t.length) + seqs[i] + '\n'; });
    out += '  ;\nEND;\n';
    if (opts.charsets && opts.charsets.length) {
      out += '\nBEGIN SETS;\n';
      opts.charsets.forEach(c => { out += `  CHARSET ${c.name.replace(/\s+/g, '_')} = ${c.from}-${c.to}${c.every ? '\\' + c.every : ''};\n`; });
      if (opts.charsets.length > 1) {
        out += '  CHARPARTITION genes = ' + opts.charsets.map((c, i) => `${c.name.replace(/\s+/g, '_')}:${i + 1}`).join(', ') + ';\n';
      }
      out += 'END;\n';
    }
    if (opts.trees && opts.trees.length) {
      out += '\nBEGIN TREES;\n';
      opts.trees.forEach(t => { out += `  TREE ${t.name || 'tree'} = ${t.newick.replace(/;\s*$/, '')};\n`; });
      out += 'END;\n';
    }
    return out;
  }
  /* the partition file RAxML and IQ-TREE read */
  function writeRaxmlPartitions(parts) {
    return parts.map(p => `${p.model || (p.type === 'aa' ? 'WAG' : 'DNA')}, ${p.name.replace(/\s+/g, '_')} = ${p.from}-${p.to}${p.every ? '\\' + p.every : ''}`).join('\n') + '\n';
  }
  /* reading that same file back */
  function parsePartitionText(text) {
    const out = [];
    clean(text).split('\n').forEach(line => {
      const l = line.trim();
      if (!l || l[0] === '#') return;
      const m = l.match(/^(?:([^,]+)\s*,\s*)?([^=]+)=\s*(.+)$/);
      if (!m) return;
      const model = (m[1] || '').trim();
      const name = m[2].trim();
      m[3].split(',').forEach(rangeTxt => {
        const r = rangeTxt.trim().match(/^(\d+)\s*-\s*(\d+)(?:\\(\d+))?$/) || rangeTxt.trim().match(/^(\d+)$/);
        if (!r) return;
        out.push({ name, model, from: +r[1], to: r[2] ? +r[2] : +r[1], every: r[3] ? +r[3] : null });
      });
    });
    return out;
  }

  Object.assign(SeqIO, {
    detect, parse, guessType,
    parseFasta, parsePhylip, parseNexus, parseClustal, parseMega, parseGenbank, parseTable,
    writeFasta, writePhylip, writeNexus, writeRaxmlPartitions, parsePartitionText,
  });
  window.SeqIO = SeqIO;
})();
