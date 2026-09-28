/* PhylogenyPro — unit tests of Block 2 (formats, aligner, trimming, quality).

   The reference values come from:
     · small files written by hand, where the right answer is obvious;
     · the files written by ape 5.8.1 in data/ejemplo_bursera/ (PHYLIP
       sequential and interleaved, NEXUS), which must read back identical to
       the FASTA of the same matrix;
     · the true alignment that generated the simulated data, which the aligner
       is scored against (validation/block2/align_bench.html);
     · properties that must hold whatever the implementation. */

(function () {

  sec('11 · Lectura de formatos');

  const FASTA = '>Taxon_A\nACGT-ACGT\nAAAA\n>Taxon_B\nACGTTACGT\nAAAA\n>Taxon_C\nACGTTACGTAAAA\n';
  let r = SeqIO.parse(FASTA, 'x.fasta');
  t('FASTA: tres secuencias', r.taxa.length === 3);
  t('FASTA: nombres', r.taxa.join(',') === 'Taxon_A,Taxon_B,Taxon_C');
  t('FASTA: líneas unidas', r.seqs[0] === 'ACGT-ACGTAAAA');
  t('FASTA: detectado como alineado', r.aligned === true);
  t('FASTA: tipo ADN', r.type === 'dna');

  const PHY_RELAXED = ' 3 13\nTaxon_A   ACGT-ACGTAAAA\nTaxon_B   ACGTTACGTAAAA\nTaxon_C   ACGTTACGTAAAA\n';
  r = SeqIO.parse(PHY_RELAXED, 'x.phy');
  t('PHYLIP relajado: nombres completos', r.taxa.join(',') === 'Taxon_A,Taxon_B,Taxon_C');
  t('PHYLIP relajado: secuencias', r.seqs[0] === 'ACGT-ACGTAAAA');

  const PHY_STRICT = ' 3 13\nTaxonA    ACGT-ACGTAAAA\nTaxonB    ACGTTACGTAAAA\nTaxonC    ACGTTACGTAAAA\n';
  r = SeqIO.parse(PHY_STRICT, 'x.phy');
  t('PHYLIP estricto: se lee igual', r.seqs[0] === 'ACGT-ACGTAAAA' && r.taxa[0] === 'TaxonA');

  const PHY_INT = ' 3 13\nTaxon_A   ACGT-ACGT\nTaxon_B   ACGTTACGT\nTaxon_C   ACGTTACGT\n\nAAAA\nAAAA\nAAAA\n';
  r = SeqIO.parse(PHY_INT, 'x.phy');
  t('PHYLIP intercalado: bloques unidos', r.seqs[0] === 'ACGT-ACGTAAAA' && r.seqs.length === 3);

  const NEX = `#NEXUS
[ a comment that must be ignored ]
BEGIN DATA;
  DIMENSIONS NTAX=3 NCHAR=13;
  FORMAT DATATYPE=DNA MISSING=? GAP=- INTERLEAVE=NO;
  MATRIX
    Taxon_A ACGT-ACGTAAAA
    Taxon_B ACGTTACGTAAAA
    Taxon_C ACGTTACGTAAAA
  ;
END;
BEGIN SETS;
  CHARSET gen1 = 1-9;
  CHARSET gen2 = 10-13;
END;
BEGIN TREES;
  TREE uno = ((Taxon_A,Taxon_B),Taxon_C);
END;`;
  r = SeqIO.parse(NEX, 'x.nex');
  t('NEXUS: matriz leída', r.seqs[0] === 'ACGT-ACGTAAAA');
  t('NEXUS: comentarios ignorados', r.taxa.length === 3);
  t('NEXUS: dos charsets', r.charsets.length === 2 && r.charsets[1].from === 10 && r.charsets[1].to === 13);
  t('NEXUS: un árbol', r.trees.length === 1 && /Taxon_C/.test(r.trees[0].newick));

  const NEX_MATCH = `#NEXUS
BEGIN DATA;
  DIMENSIONS NTAX=2 NCHAR=5;
  FORMAT DATATYPE=DNA MATCHCHAR=. GAP=-;
  MATRIX
    A ACGTA
    B ..T..
  ;
END;`;
  r = SeqIO.parse(NEX_MATCH, 'y.nex');
  t('NEXUS: MATCHCHAR resuelto contra la primera fila', r.seqs[1] === 'ACTTA');

  const CLU = 'CLUSTAL W (1.83) multiple sequence alignment\n\nTaxon_A   ACGT-ACGT\nTaxon_B   ACGTTACGT\n          **** ****\n\nTaxon_A   AAAA\nTaxon_B   AAAA\n          ****\n';
  r = SeqIO.parse(CLU, 'x.aln');
  t('Clustal: dos secuencias, línea de conservación ignorada', r.taxa.length === 2 && r.seqs[0] === 'ACGT-ACGTAAAA');

  const MEG = '#MEGA\n!Title Prueba;\n!Format DataType=DNA;\n\n#Taxon_A\nACGT-ACGT\nAAAA\n#Taxon_B\nACGTTACGTAAAA\n';
  r = SeqIO.parse(MEG, 'x.meg');
  t('MEGA: leído', r.taxa.join(',') === 'Taxon_A,Taxon_B' && r.seqs[0] === 'ACGT-ACGTAAAA');

  const GB = `LOCUS       AB000001    13 bp    DNA     linear   PLN 01-JAN-2020
ACCESSION   AB000001
  ORGANISM  Bursera simaruba
            Eukaryota; Viridiplantae.
FEATURES             Location/Qualifiers
     gene            1..13
                     /gene="rbcL"
ORIGIN
        1 acgtaacgta aaa
//
`;
  r = SeqIO.parse(GB, 'x.gb');
  t('GenBank: nombre del organismo y del gen', r.taxa[0] === 'Bursera_simaruba_rbcL');
  t('GenBank: secuencia sin números ni espacios', r.seqs[0] === 'ACGTAACGTAAAA');

  const CSV = 'taxon,car1,car2,car3\nT1,0,1,?\nT2,1,1,0\nT3,0,0,1\n';
  r = SeqIO.parse(CSV, 'x.csv');
  t('Tabla morfológica: tipo morph', r.type === 'morph');
  t('Tabla morfológica: estados', r.seqs[0] === '01?' && r.taxa[2] === 'T3');

  r = SeqIO.parse('>a\nACGT\n>a\nTTTT\n', 'd.fasta');
  t('nombres repetidos: se renombra el segundo y se avisa', r.taxa[1] !== r.taxa[0] && r.warnings.length > 0);

  sec('12 · Escritura de formatos (ida y vuelta)');
  const taxa = ['Bursera_simaruba', 'Bursera_fagaroides', 'Commiphora_africana'];
  const seqs = ['ACGT-ACGTAAAA', 'ACGTTACGTAAAA', 'ACGTTACGTAAAT'];
  const back = (txt, name) => SeqIO.parse(txt, name);
  let b = back(SeqIO.writeFasta(taxa, seqs), 'x.fasta');
  t('FASTA ida y vuelta', b.taxa.join() === taxa.join() && b.seqs.join() === seqs.join());
  b = back(SeqIO.writePhylip(taxa, seqs, {}), 'x.phy');
  t('PHYLIP secuencial ida y vuelta', b.taxa.join() === taxa.join() && b.seqs.join() === seqs.join());
  b = back(SeqIO.writePhylip(taxa, seqs, { interleaved: true, block: 5 }), 'x.phy');
  t('PHYLIP intercalado ida y vuelta', b.taxa.join() === taxa.join() && b.seqs.join() === seqs.join());
  b = back(SeqIO.writeNexus(taxa, seqs, { type: 'dna', charsets: [{ name: 'g1', from: 1, to: 9 }, { name: 'g2', from: 10, to: 13 }] }), 'x.nex');
  t('NEXUS ida y vuelta, con charsets', b.taxa.join() === taxa.join() && b.seqs.join() === seqs.join() && b.charsets.length === 2);
  const partTxt = SeqIO.writeRaxmlPartitions([{ name: 'rbcL', type: 'dna', from: 1, to: 620 }, { name: 'ITS', type: 'dna', from: 621, to: 1200 }]);
  const parts = SeqIO.parsePartitionText(partTxt);
  t('esquema de particiones ida y vuelta', parts.length === 2 && parts[1].from === 621 && parts[1].to === 1200);
  t('esquema de particiones: se lee el formato de RAxML con modelo', SeqIO.parsePartitionText('DNA, gene1 = 1-100\nWAG, gene2 = 101-200\n').length === 2);

  sec('13 · Alineador');
  /* two sequences that differ by one insertion: the gap must fall where the
     insertion is, and nowhere else */
  let res = Align.progressive(['ACGTACGTACGT', 'ACGTAAACGTACGT'], 'dna', { refine: 0 });
  t('un hueco de dos bases en la secuencia corta', res.aligned[0].length === 14 && (res.aligned[0].match(/-/g) || []).length === 2);
  t('la secuencia larga queda sin huecos', res.aligned[1].indexOf('-') < 0);
  res = Align.progressive(['ACGTACGT', 'ACGTACGT', 'ACGTACGT'], 'dna', { refine: 0 });
  t('secuencias idénticas: ningún hueco', res.aligned.every(s => s === 'ACGTACGT'));
  t('el alineamiento no pierde residuos', (() => {
    const inp = ['ACGTTTACGT', 'ACGTACGT', 'ACGTTTTTACGTA'];
    const r2 = Align.progressive(inp, 'dna', { refine: 1 });
    return r2.aligned.every((s, i) => s.replace(/-/g, '') === inp[i]);
  })());
  t('todas las filas quedan del mismo largo', (() => {
    const r2 = Align.progressive(['ACGTACGTAA', 'ACGTAA', 'ACGTTTACGTAA', 'AC'], 'dna', { refine: 1 });
    return r2.aligned.every(s => s.length === r2.aligned[0].length);
  })());
  t('k-meros: distancia de una secuencia consigo misma es 0', (() => {
    const D = Align.kmerDistance(['ACGTACGTACGTAC', 'ACGTACGTACGTAC'], 'dna', 4);
    return Math.abs(D[0][1]) < 1e-12;
  })());
  t('k-meros: dos secuencias sin nada en común dan distancia 1', (() => {
    const D = Align.kmerDistance(['AAAAAAAAAA', 'CCCCCCCCCC'], 'dna', 4);
    return Math.abs(D[0][1] - 1) < 1e-12;
  })());
  t('UPGMA: junta primero a los más parecidos', (() => {
    const D = [[0, 0.1, 0.9], [0.1, 0, 0.9], [0.9, 0.9, 0]];
    const tree = Align.upgma(D, ['a', 'b', 'c']);
    const first = tree.children.find(c => c.children.length);
    return first && first.size === 2;
  })());
  t('BLOSUM62 es simétrica', (() => {
    const M = Align.BLOSUM62_ROWS;
    for (let i = 0; i < 20; i++) for (let j = 0; j < 20; j++) if (M[i][j] !== M[j][i]) return false;
    return true;
  })());
  t('BLOSUM62: valores conocidos (W-W = 11, C-C = 9, A-A = 4)', (() => {
    const A = Align.AA, M = Align.BLOSUM62_ROWS;
    return M[A.indexOf('W')][A.indexOf('W')] === 11 && M[A.indexOf('C')][A.indexOf('C')] === 9 && M[A.indexOf('A')][A.indexOf('A')] === 4;
  })());
  t('BLOSUM62: la diagonal es positiva en los veinte', Align.BLOSUM62_ROWS.every((row, i) => row[i] > 0));
  t('la puntuación mejora al alinear bien', (() => {
    const good = ['ACGTACGT', 'ACGTACGT'];
    const bad = ['ACGTACGT', 'TGCATGCA'];
    return Align.spScore(good, 'dna') > Align.spScore(bad, 'dna');
  })());
  t('comparar con una referencia: idéntica da SP = 1', (() => {
    const a = ['AC-GT', 'ACGGT'];
    return Math.abs(Align.compareToReference(a, a).sp - 1) < 1e-12;
  })());
  t('comparar con una referencia detecta un desfase', (() => {
    const ref = ['AC-GT', 'ACGGT'], test = ['ACG-T', 'ACGGT'];
    return Align.compareToReference(test, ref).sp < 1;
  })());

  sec('14 · Traducción y codones');
  t('código estándar: ATG = M, TGG = W, TAA = paro', Align.translateCodon('ATG') === 'M' && Align.translateCodon('TGG') === 'W' && Align.translateCodon('TAA') === '*');
  t('código estándar: los seis codones de leucina', ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'].every(c => Align.translateCodon(c) === 'L'));
  t('código estándar: los seis de arginina', ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'].every(c => Align.translateCodon(c) === 'R'));
  t('código estándar: los tres de paro', ['TAA', 'TAG', 'TGA'].every(c => Align.translateCodon(c) === '*'));
  t('traducción completa', Align.translate('ATGGCCTGGTAA') === 'MAW*');
  t('marco 2 cambia la traducción', Align.translate('AATGGCCTGGTAA', 2) === 'MAW*');
  t('codones de paro internos detectados, el final no cuenta', Align.stopCodons('ATGTAAATGTAA').length === 1);
  /* TAA GCC TGG … in frame 1 (an internal stop), AAG CCT GGG … in frame 2 (clean) */
  t('el mejor marco es el que no tiene paros', (() => {
    const s = 'TAAGCCTGGGCCTGGTTT';
    return Align.stopCodons(s, 1).length === 1 && Align.stopCodons(s, 2).length === 0 && Align.bestFrame([s]).frame === 2;
  })());
  t('alineamiento por codón: ningún hueco parte un codón', (() => {
    const seqs = ['ATGGCCTGGGCCTGG', 'ATGGCCGCCTGG'];
    const r2 = Align.codonAlign(seqs, { frame: 1, refine: 0 });
    return r2.aligned.every(s => {
      for (let i = 0; i < s.length; i += 3) {
        const cod = s.slice(i, i + 3);
        if (cod.indexOf('-') >= 0 && cod !== '---') return false;
      }
      return true;
    });
  })());
  t('alineamiento por codón: longitud múltiplo de tres', (() => {
    const r2 = Align.codonAlign(['ATGGCCTGGGCCTGG', 'ATGGCCGCCTGG'], { frame: 1, refine: 0 });
    return r2.length % 3 === 0;
  })());

  sec('15 · Recorte y estadísticas de columnas');
  const ALN = ['ACGTACGTAC', 'ACGTACGTAC', 'ACGTACGTAC', 'ACGAAC--AC'];
  let st = Trim.columnStats(ALN);
  t('estadísticas: columna constante', st[0].states === 1 && st[0].gaps === 0);
  t('estadísticas: columna con hueco', st[6].gaps === 1);
  t('estadísticas: columna variable', st[3].states === 2);
  let g = Trim.byGaps(ALN, 0.2);
  t('recorte por huecos: se quitan las columnas con hueco', g.kept === 8);
  g = Trim.byGaps(ALN, 0.5);
  t('recorte por huecos con umbral alto: no se quita nada', g.kept === 10);
  const sim = Trim.similarity(ALN);
  t('similitud: columna idéntica = 1', Math.abs(sim[0] - 1) < 1e-12);
  t('similitud: columna con tres iguales y uno distinto = 0.5', Math.abs(sim[3] - 0.5) < 1e-12);
  t('Gblocks conserva el bloque conservado de un alineamiento limpio', (() => {
    const clean = Array.from({ length: 8 }, () => 'ACGTACGTACGTACGTACGTACGT');
    const gb = Trim.gblocks(clean, {});
    return gb.kept === 24;
  })());
  t('Gblocks quita una región no conservada larga', (() => {
    const n = 8, L = 60;
    const seqs2 = [];
    for (let i = 0; i < n; i++) {
      let s = '';
      for (let c = 0; c < L; c++) s += (c >= 20 && c < 40) ? 'ACGT'[(i + c) % 4] : 'A';
      seqs2.push(s);
    }
    const gb = Trim.gblocks(seqs2, {});
    let insideKept = 0;
    for (let c = 20; c < 40; c++) if (gb.keep[c]) insideKept++;
    return insideKept === 0 && gb.kept === 40;
  })());
  t('recorte por codón: se conservan tríos completos', (() => {
    const keep = new Uint8Array([1, 1, 0, 1, 1, 1]);
    const out = Trim.codonAware(keep, 1);
    return out[0] === 0 && out[1] === 0 && out[2] === 0 && out[3] === 1 && out[4] === 1 && out[5] === 1;
  })());
  t('aplicar el recorte deja las columnas elegidas', Trim.apply(['ABCDE', 'FGHIJ'], new Uint8Array([1, 0, 1, 0, 1]))[0] === 'ACE');
  t('bloques: se reportan como rangos de base 1', (() => {
    const bl = Trim.blocks(new Uint8Array([0, 1, 1, 0, 1]));
    return bl.length === 2 && bl[0][0] === 2 && bl[0][1] === 3 && bl[1][0] === 5;
  })());

  sec('16 · Control de calidad');
  let q = QC.check(['a', 'a', 'b'], ['ACGT', 'ACGT', 'ACGT'], 'dna');
  t('detecta nombres repetidos', q.some(x => x.level === 'bad' && /repetido|repeated/.test(x.es + x.en)));
  q = QC.check(['a', 'b'], ['ACGTACGTAC', 'AC'], 'dna');
  t('detecta una secuencia mucho más corta', q.some(x => /mitad|half/.test(x.es)));
  q = QC.check(['a', 'b'], ['ACGTZZ', 'ACGTAA'], 'dna');
  t('detecta caracteres no esperados', q.some(x => /no esperados|not expected/.test(x.es + x.en)));
  q = QC.check(['a', 'b'], ['ATGTAAATGA', 'ATGGCCATGA'], 'dna', { coding: true, frame: 1 });
  t('detecta codones de paro internos', q.some(x => x.level === 'bad' && /paro|stop/.test(x.es + x.en)));
  q = QC.check(['a', 'b'], ['ATGGCCTGA', 'ATGGCCTGA'], 'dna', { coding: true, frame: 1 });
  t('no inventa paros donde no los hay', !q.some(x => x.level === 'bad'));
  t('taxones faltantes entre particiones', (() => {
    const m = QC.missingTaxa([{ name: 'g1', taxa: ['a', 'b'] }, { name: 'g2', taxa: ['a'] }]);
    return m.rows.length === 1 && m.rows[0].taxon === 'b';
  })());

  sec('17 · Genoma de cada marcador');
  t('rbcL, matK, trnL-F y ycf1 son de cloroplasto', ['rbcL', 'matK', 'trnL-F', 'ycf1', 'ndhF'].every(n => B2.guessGenome(n) === 'cp'));
  t('ITS y ETS son nucleares ribosomales', ['ITS', 'ITS2', 'ETS'].every(n => B2.guessGenome(n) === 'nr'));
  t('matR, nad1 y cox1 son mitocondriales', ['matR', 'nad1', 'cox1', 'atp1'].every(n => B2.guessGenome(n) === 'mt'));
  t('waxy, PHYC y LEAFY son nucleares', ['waxy', 'PHYC', 'LEAFY', 'G3pdh'].every(n => B2.guessGenome(n) === 'nu'));

  sec('18 · Datos de ejemplo');
  t('los ejemplos están incrustados', !!window.EXAMPLES && Object.keys(EXAMPLES.seqs).length >= 5);
  t('Bursera: 18 taxones en los tres genes', ['rbcL', 'trnLF', 'ITS'].every(k => Object.keys(EXAMPLES.seqs[k]).length === 18));
  t('el concatenado está alineado', (() => {
    const s = Object.values(EXAMPLES.seqs.concatenado);
    return s.every(x => x.length === s[0].length);
  })());
  t('los genes sueltos vienen sin alinear', (() => {
    const s = Object.values(EXAMPLES.seqs.rbcL);
    return s.some(x => x.length !== s[0].length) || s.every(x => x.indexOf('-') < 0);
  })());
  t('el esquema de particiones del ejemplo cubre todo el concatenado', (() => {
    const parts = SeqIO.parsePartitionText(EXAMPLES.particiones);
    const L = Object.values(EXAMPLES.seqs.concatenado)[0].length;
    return parts.length === 3 && parts[0].from === 1 && parts[2].to === L;
  })());
  t('el árbol verdadero de Bursera trae los 18 taxones', (EXAMPLES.trees.bursera.match(/Bursera_|Commiphora_|Boswellia_/g) || []).length === 18);
  t('rbcL del ejemplo no tiene codones de paro internos en el marco 1', (() => {
    const seqs2 = Object.values(EXAMPLES.seqs.rbcL);
    return seqs2.every(s => Align.stopCodons(s, 1).length === 0);
  })());
  t('la matriz morfológica se lee', (() => {
    const m = SeqIO.parse(EXAMPLES.morfologia, 'm.csv');
    return m.type === 'morph' && m.taxa.length === 25 && m.seqs[0].length === 42;
  })());

})();
