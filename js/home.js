/* PhylogenyPro — molecular phylogenetics in the browser.
   Copyright (C) 2026  Luis Ángel Barrera-Guzmán

   This program is free software: you can redistribute it and/or modify it under
   the terms of the GNU General Public License as published by the Free Software
   Foundation, either version 3 of the License, or (at your option) any later
   version. It is distributed in the hope that it will be useful, but WITHOUT ANY
   WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
   A PARTICULAR PURPOSE. See the GNU General Public License, in the file LICENSE
   at the root of this program, or <https://www.gnu.org/licenses/>. */

/* PhylogenyPro — Block 1: the home page.
   Builds the stepper, the block cards, the data strip, the marker list, the
   method gallery, the comparison of inference criteria, the table of sizes the
   app can really handle, and the reference list. The content lives here as
   bilingual data, so the map of the app is one editable list. */

(function () {

  /* ---------------- the eleven blocks that follow the home page ---------------- */
  const BLOCKS = [
    { n: 2, art: 'blkData', tag: ['datos', 'data'],
      t: ['Datos y alineamiento', 'Data and alignment'],
      d: ['Lee FASTA, PHYLIP, NEXUS, Clustal, MEGA y GenBank; alinea tus secuencias con un alineador progresivo propio (o por codones), deja editar el alineamiento a mano, recorta los bloques ambiguos, concatena genes y arma el mapa de particiones declarando el genoma de cada una.',
        'Reads FASTA, PHYLIP, NEXUS, Clustal, MEGA and GenBank; aligns your sequences with its own progressive aligner (or by codons), lets you edit the alignment by hand, trims the ambiguous blocks, concatenates genes and builds the partition map, declaring the genome each one belongs to.'] },
    { n: 3, art: 'blkModels', tag: ['modelos', 'models'],
      t: ['Modelos de sustitución', 'Substitution models'],
      d: ['Composición de bases, prueba de saturación y selección del modelo entre los de ADN (JC a GTR, con +F, +I, +G y +R), los empíricos de proteína (LG, WAG, JTT, cpREV…), los de codón y los de morfología (Mk, Mkv), por AIC, AICc y BIC, con el esquema de particiones elegido por fusión.',
        'Base composition, test of saturation and model selection among the DNA models (JC to GTR, with +F, +I, +G and +R), the empirical protein matrices (LG, WAG, JTT, cpREV…), codon models and morphology (Mk, Mkv), by AIC, AICc and BIC, with the partitioning scheme chosen by merging.'] },
    { n: 4, art: 'blkParsimony', tag: ['árboles rápidos', 'quick trees'],
      t: ['Parsimonia y distancias', 'Parsimony and distances'],
      d: ['Distancias corregidas por modelo, NJ, BIONJ, UPGMA y mínima evolución; máxima parsimonia con Fitch y Sankoff, búsqueda heurística con TBR, índices de consistencia y retención, soporte de Bremer, bootstrap, jackknife y árboles de consenso.',
        'Model-corrected distances, NJ, BIONJ, UPGMA and minimum evolution; maximum parsimony with Fitch and Sankoff, heuristic search with TBR, consistency and retention indices, Bremer support, bootstrap, jackknife and consensus trees.'] },
    { n: 5, art: 'blkLikelihood', tag: ['verosimilitud', 'likelihood'],
      t: ['Máxima verosimilitud', 'Maximum likelihood'],
      d: ['Búsqueda de topología con NNI y SPR desde varios arranques, longitudes y parámetros optimizados a la vez, particiones con tasas ligadas o libres; soporte por bootstrap, bootstrap ultrarrápido, SH-aLRT y aBayes; pruebas de topología KH, SH y AU, y árboles restringidos para poner a prueba una hipótesis de monofilia.',
        'Topology search with NNI and SPR from several starts, branch lengths and parameters optimised together, partitions with linked or free rates; support by bootstrap, ultrafast bootstrap, SH-aLRT and aBayes; KH, SH and AU topology tests, and constrained trees to put a hypothesis of monophyly to the test.'] },
    { n: 6, art: 'blkBayes', tag: ['bayesiano', 'Bayesian'],
      t: ['Inferencia bayesiana', 'Bayesian inference'],
      d: ['MCMC sobre el espacio de árboles con cadenas acopladas (MC³) repartidas en hilos del navegador, distribuciones previas explícitas, diagnóstico de convergencia en vivo (traza, ESS, ASDSF entre corridas), consenso de mayoría con probabilidades posteriores y bitácora exportable a Tracer.',
        'MCMC over tree space with coupled chains (MC³) spread across browser threads, explicit priors, live convergence diagnostics (trace, ESS, ASDSF between runs), majority-rule consensus with posterior probabilities and a log file that opens in Tracer.'] },
    { n: 7, art: 'blkDating', tag: ['tiempo', 'time'],
      t: ['Reloj molecular y tiempos', 'Molecular clock and dating'],
      d: ['Enraizamiento por grupo externo, punto medio o mínima varianza; pruebas de reloj y regresión raíz–punta para muestras fechadas; datación por mínimos cuadrados y por MCMC con reloj estricto o relajado, calibraciones fósiles con su distribución, intervalos de credibilidad y escala geológica.',
        'Rooting by outgroup, midpoint or minimum variance; clock tests and root-to-tip regression for dated samples; least-squares and MCMC dating under a strict or relaxed clock, fossil calibrations with their own distribution, credibility intervals and a geological time scale.'] },
    { n: 8, art: 'blkDiversification', tag: ['diversificación', 'diversification'],
      t: ['Diversificación', 'Diversification'],
      d: ['Curvas de linajes en el tiempo, estadístico γ, modelos de nacimiento y muerte con muestreo incompleto, modelos episódicos y detección de clados que se diversificaron más rápido que el resto.',
        'Lineages-through-time curves, the γ statistic, birth–death models with incomplete sampling, episodic models and detection of the clades that diversified faster than the rest.'] },
    { n: 9, art: 'blkTraits', tag: ['caracteres', 'traits'],
      t: ['Evolución de caracteres', 'Trait evolution'],
      d: ['Estados ancestrales por parsimonia y por el modelo Mk (ER, SYM, ARD) con mapeo estocástico; caracteres continuos con movimiento browniano, Ornstein–Uhlenbeck y explosión temprana; contrastes independientes, regresión filogenética y señal filogenética (K de Blomberg, λ de Pagel).',
        'Ancestral states by parsimony and under the Mk model (ER, SYM, ARD) with stochastic mapping; continuous traits under Brownian motion, Ornstein–Uhlenbeck and early burst; independent contrasts, phylogenetic regression and phylogenetic signal (Blomberg\'s K, Pagel\'s λ).'] },
    { n: 10, art: 'blkBiogeo', tag: ['biogeografía', 'biogeography'],
      t: ['Biogeografía histórica', 'Historical biogeography'],
      d: ['Rangos ancestrales con DEC, DEC+J, DIVALIKE y BAYAREALIKE comparados por AICc, con restricciones de dispersión por intervalos de tiempo, áreas adyacentes y tamaño máximo de rango; mapeo estocástico biogeográfico y el mapa de las áreas junto al árbol fechado.',
        'Ancestral ranges under DEC, DEC+J, DIVALIKE and BAYAREALIKE compared by AICc, with time-stratified dispersal constraints, adjacency and a maximum range size; biogeographical stochastic mapping and the map of the areas beside the dated tree.'] },
    { n: 11, art: 'blkCompare', tag: ['comparar', 'compare'],
      t: ['Comparación de árboles', 'Comparing trees'],
      d: ['Distancias entre árboles (Robinson–Foulds, cuartetos, camino), tanglegramas, consensos, concordancia de genes y de sitios, árbol de especies por cuartetos, redes de divisiones y la prueba ABBA-BABA: qué significa que el árbol del cloroplasto y el nuclear no coincidan.',
        'Tree distances (Robinson–Foulds, quartets, path), tanglegrams, consensus trees, gene and site concordance, a quartet species tree, split networks and the ABBA-BABA test: what it means when the chloroplast tree and the nuclear tree disagree.'] },
    { n: 12, art: 'blkReport', tag: ['publicar', 'publish'],
      t: ['Figuras e informe', 'Figures and report'],
      d: ['Estudio de árboles con presentación rectangular, circular o no enraizada, clados coloreados y colapsados, soportes múltiples, barras de credibilidad, escala geológica, anillos e imágenes de las puntas; exportación en Newick, NEXUS y PhyloXML, figuras hasta 900 ppp, informe con la metodología redactada y ZIP reproducible.',
        'A tree studio with rectangular, circular or unrooted layouts, coloured and collapsed clades, several support values, credibility bars, a geological scale, rings and images at the tips; export as Newick, NEXUS and PhyloXML, figures at up to 900 dpi, a report with the methods already written and a reproducible ZIP.'] },
  ];

  /* ---------------- what can be analysed ---------------- */
  const MATERIALS = [
    { art: 'matCp', k: ['cloroplasto', 'chloroplast'], kc: 'l', t: ['ADN de cloroplasto', 'Chloroplast DNA'], s: ['rbcL, matK, trnL-F, ndhF, ycf1: herencia materna, sin recombinación, evolución lenta.', 'rbcL, matK, trnL-F, ndhF, ycf1: maternal inheritance, no recombination, slow evolution.'] },
    { art: 'matNuclear', k: ['nuclear', 'nuclear'], kc: '', t: ['ADN nuclear e ITS', 'Nuclear DNA and ITS'], s: ['ITS, ETS y genes de copia baja: biparentales, rápidos y con recombinación.', 'ITS, ETS and low-copy genes: biparental, fast and recombining.'] },
    { art: 'matProtein', k: ['proteínas', 'proteins'], kc: 'g', t: ['Secuencias de aminoácidos', 'Amino-acid sequences'], s: ['Familias génicas y filogenias profundas con LG, WAG, JTT o cpREV.', 'Gene families and deep phylogenies with LG, WAG, JTT or cpREV.'] },
    { art: 'matCodon', k: ['codones', 'codons'], kc: 'g', t: ['Alineamientos por codón', 'Codon alignments'], s: ['Particiones por posición y modelos de codón para medir selección (dN/dS).', 'Partitions by position and codon models to measure selection (dN/dS).'] },
    { art: 'matMorph', k: ['morfología', 'morphology'], kc: '', t: ['Caracteres morfológicos', 'Morphological characters'], s: ['Matrices discretas con el modelo Mk, solas o combinadas con moléculas.', 'Discrete matrices under the Mk model, alone or combined with molecules.'] },
    { art: 'matDated', k: ['fechas', 'dates'], kc: 'l', t: ['Fósiles y fechas de muestreo', 'Fossils and sampling dates'], s: ['Calibraciones con su distribución, o puntas fechadas para patógenos y poblaciones.', 'Calibrations with their own distribution, or dated tips for pathogens and populations.'] },
  ];

  /* ---------------- markers the app recognises by name ---------------- */
  const MARKERS = [
    { g: 'cp', n: 'rbcL', d: ['El caballo de batalla: lento, fácil de amplificar, bueno de familia hacia arriba.', 'The workhorse: slow, easy to amplify, good from family level upwards.'] },
    { g: 'cp', n: 'matK', d: ['Más variable que rbcL; el otro marcador del código de barras de plantas.', 'More variable than rbcL; the other plant barcode marker.'] },
    { g: 'cp', n: 'trnL-F', d: ['Espaciador no codificante: resuelve entre especies cercanas.', 'Non-coding spacer: resolves between close species.'] },
    { g: 'cp', n: 'trnH-psbA', d: ['El espaciador más variable del cloroplasto; con muchos indeles.', 'The most variable chloroplast spacer; full of indels.'] },
    { g: 'cp', n: 'ndhF · atpB · rpl16 · rps16 · ycf1', d: ['El resto del repertorio habitual en sistemática de plantas.', 'The rest of the usual repertoire in plant systematics.'] },
    { g: 'nr', n: 'ITS · ITS2 · ETS', d: ['Nuclear ribosomal: rápido y muy usado, pero con copias divergentes y evolución concertada incompleta.', 'Nuclear ribosomal: fast and widely used, but with divergent copies and incomplete concerted evolution.'] },
    { g: 'nu', n: 'waxy · PHYC · LEAFY · G3pdh', d: ['Nucleares de copia baja: biparentales, detectan hibridación.', 'Low-copy nuclear genes: biparental, they detect hybridisation.'] },
    { g: 'mt', n: 'matR · nad1 · cox1 · atp1', d: ['Mitocondriales de planta: los más lentos, útiles en filogenias profundas.', 'Plant mitochondrial: the slowest, useful in deep phylogenies.'] },
  ];
  const GENOME_NAMES = {
    cp: ['cloroplasto', 'chloroplast'], mt: ['mitocondria', 'mitochondrion'],
    nr: ['nuclear ribosomal', 'nuclear ribosomal'], nu: ['nuclear', 'nuclear'],
  };

  /* ---------------- method gallery ---------------- */
  const FAMS = {
    dat: ['datos', 'data'], mod: ['modelos', 'models'], tre: ['árboles rápidos', 'quick trees'],
    inf: ['inferencia', 'inference'], tim: ['tiempo', 'time'], mac: ['macroevolución', 'macroevolution'],
    bio: ['biogeografía', 'biogeography'], cmp: ['comparar y publicar', 'compare and publish'],
  };
  const METHODS = [
    { art: 'mFasta', fam: 'dat', n: ['Lectura de formatos', 'Format readers'], s: ['FASTA, PHYLIP, NEXUS, Clustal, MEGA', 'FASTA, PHYLIP, NEXUS, Clustal, MEGA'] },
    { art: 'mAlign', fam: 'dat', n: ['Alineamiento progresivo', 'Progressive alignment'], s: ['k-meros, árbol guía y perfil–perfil', 'k-mers, guide tree and profile–profile'] },
    { art: 'mCodonAlign', fam: 'dat', n: ['Alineamiento por codón', 'Codon-aware alignment'], s: ['traduce, alinea y regresa', 'translate, align and back-translate'] },
    { art: 'mTrim', fam: 'dat', n: ['Recorte de bloques', 'Block trimming'], s: ['columnas ambiguas y huecos', 'ambiguous columns and gaps'] },
    { art: 'mPartition', fam: 'dat', n: ['Particiones y genomas', 'Partitions and genomes'], s: ['cloroplasto, mitocondria, núcleo', 'chloroplast, mitochondrion, nucleus'] },
    { art: 'mQmatrix', fam: 'mod', n: ['Matriz de sustitución', 'Substitution matrix'], s: ['de JC69 a GTR', 'from JC69 to GTR'] },
    { art: 'mModelTest', fam: 'mod', n: ['Selección de modelo', 'Model selection'], s: ['AIC, AICc, BIC y pesos', 'AIC, AICc, BIC and weights'] },
    { art: 'mGamma', fam: 'mod', n: ['Tasas variables +Γ', 'Rate variation +Γ'], s: ['gamma discreta y categorías libres', 'discrete gamma and free rates'] },
    { art: 'mInvariant', fam: 'mod', n: ['Sitios invariables +I', 'Invariable sites +I'], s: ['la fracción que nunca cambia', 'the fraction that never changes'] },
    { art: 'mSaturation', fam: 'mod', n: ['Saturación', 'Saturation'], s: ['índice de Xia y gráfica ts/tv', 'Xia\'s index and the ts/tv plot'] },
    { art: 'mComposition', fam: 'mod', n: ['Composición de bases', 'Base composition'], s: ['prueba de homogeneidad', 'homogeneity test'] },
    { art: 'mProteinModel', fam: 'mod', n: ['Modelos de proteína', 'Protein models'], s: ['LG, WAG, JTT, cpREV, mtREV', 'LG, WAG, JTT, cpREV, mtREV'] },
    { art: 'mMk', fam: 'mod', n: ['Modelo Mk de morfología', 'Mk model for morphology'], s: ['Lewis 2001, con corrección Mkv', 'Lewis 2001, with the Mkv correction'] },
    { art: 'mDistance', fam: 'tre', n: ['Matrices de distancia', 'Distance matrices'], s: ['p, JC, K2P, TN93, LogDet', 'p, JC, K2P, TN93, LogDet'] },
    { art: 'mNJ', fam: 'tre', n: ['Vecino más cercano', 'Neighbour joining'], s: ['Saitou y Nei 1987', 'Saitou and Nei 1987'] },
    { art: 'mBionj', fam: 'tre', n: ['BIONJ', 'BIONJ'], s: ['NJ con varianzas', 'NJ with variances'] },
    { art: 'mUPGMA', fam: 'tre', n: ['UPGMA y WPGMA', 'UPGMA and WPGMA'], s: ['solo si hay reloj', 'only under a clock'] },
    { art: 'mME', fam: 'tre', n: ['Mínima evolución', 'Minimum evolution'], s: ['equilibrada, con NNI', 'balanced, with NNI'] },
    { art: 'mParsimony', fam: 'tre', n: ['Máxima parsimonia', 'Maximum parsimony'], s: ['Fitch, búsqueda con TBR', 'Fitch, search with TBR'] },
    { art: 'mSankoff', fam: 'tre', n: ['Parsimonia con pesos', 'Weighted parsimony'], s: ['matriz de costos de Sankoff', 'Sankoff cost matrix'] },
    { art: 'mBremer', fam: 'tre', n: ['Soporte de Bremer', 'Bremer support'], s: ['pasos extra para perder un clado', 'extra steps to lose a clade'] },
    { art: 'mConsensus', fam: 'tre', n: ['Árboles de consenso', 'Consensus trees'], s: ['estricto, mayoría, semiestricto', 'strict, majority-rule, semistrict'] },
    { art: 'mML', fam: 'inf', n: ['Árbol de máxima verosimilitud', 'Maximum-likelihood tree'], s: ['Felsenstein 1981', 'Felsenstein 1981'] },
    { art: 'mNNI', fam: 'inf', n: ['Búsqueda NNI y SPR', 'NNI and SPR search'], s: ['reordenamientos con varios arranques', 'rearrangements from several starts'] },
    { art: 'mBootstrap', fam: 'inf', n: ['Bootstrap', 'Bootstrap'], s: ['réplicas repartidas en hilos', 'replicates spread across threads'] },
    { art: 'mUFBoot', fam: 'inf', n: ['Bootstrap ultrarrápido', 'Ultrafast bootstrap'], s: ['aproximación RELL', 'RELL approximation'] },
    { art: 'mSHaLRT', fam: 'inf', n: ['SH-aLRT y aBayes', 'SH-aLRT and aBayes'], s: ['pruebas por rama', 'branch tests'] },
    { art: 'mTopoTest', fam: 'inf', n: ['Pruebas de topología', 'Topology tests'], s: ['KH, SH y AU', 'KH, SH and AU'] },
    { art: 'mConstraint', fam: 'inf', n: ['Árboles restringidos', 'Constrained trees'], s: ['poner a prueba una monofilia', 'testing a monophyly'] },
    { art: 'mMCMC', fam: 'inf', n: ['MCMC bayesiano', 'Bayesian MCMC'], s: ['topología, ramas y parámetros', 'topology, branches and parameters'] },
    { art: 'mMC3', fam: 'inf', n: ['Cadenas acopladas MC³', 'Coupled chains MC³'], s: ['una cadena caliente por hilo', 'one heated chain per thread'] },
    { art: 'mConvergence', fam: 'inf', n: ['Convergencia', 'Convergence'], s: ['ASDSF, ESS y traza', 'ASDSF, ESS and trace'] },
    { art: 'mPosterior', fam: 'inf', n: ['Probabilidades posteriores', 'Posterior probabilities'], s: ['consenso de mayoría al 50%', '50% majority-rule consensus'] },
    { art: 'mRooting', fam: 'tim', n: ['Enraizamiento', 'Rooting'], s: ['grupo externo, punto medio, varianza', 'outgroup, midpoint, variance'] },
    { art: 'mClockTest', fam: 'tim', n: ['Prueba de reloj molecular', 'Molecular clock test'], s: ['razón de verosimilitudes y tasas relativas', 'likelihood ratio and relative rates'] },
    { art: 'mLSD', fam: 'tim', n: ['Datación por mínimos cuadrados', 'Least-squares dating'], s: ['rápida, para árboles grandes', 'fast, for large trees'] },
    { art: 'mFossil', fam: 'tim', n: ['Calibraciones fósiles', 'Fossil calibrations'], s: ['uniforme, exponencial, lognormal', 'uniform, exponential, lognormal'] },
    { art: 'mRelaxed', fam: 'tim', n: ['Reloj relajado', 'Relaxed clock'], s: ['lognormal no correlacionado', 'uncorrelated lognormal'] },
    { art: 'mTipDating', fam: 'tim', n: ['Puntas fechadas', 'Dated tips'], s: ['regresión raíz–punta', 'root-to-tip regression'] },
    { art: 'mLTT', fam: 'mac', n: ['Linajes en el tiempo', 'Lineages through time'], s: ['curva LTT y estadístico γ', 'LTT curve and the γ statistic'] },
    { art: 'mBirthDeath', fam: 'mac', n: ['Nacimiento y muerte', 'Birth–death'], s: ['con muestreo incompleto', 'with incomplete sampling'] },
    { art: 'mRateShift', fam: 'mac', n: ['Cambios de tasa', 'Rate shifts'], s: ['qué clado se aceleró', 'which clade sped up'] },
    { art: 'mAncestral', fam: 'mac', n: ['Estados ancestrales', 'Ancestral states'], s: ['Mk con ER, SYM y ARD', 'Mk with ER, SYM and ARD'] },
    { art: 'mSimmap', fam: 'mac', n: ['Mapeo estocástico', 'Stochastic mapping'], s: ['historias muestreadas del carácter', 'sampled histories of the trait'] },
    { art: 'mBMOU', fam: 'mac', n: ['Browniano, OU y explosión', 'Brownian, OU and early burst'], s: ['caracteres continuos', 'continuous traits'] },
    { art: 'mSignal', fam: 'mac', n: ['Señal filogenética', 'Phylogenetic signal'], s: ['K de Blomberg y λ de Pagel', 'Blomberg\'s K and Pagel\'s λ'] },
    { art: 'mPGLS', fam: 'mac', n: ['Regresión filogenética', 'Phylogenetic regression'], s: ['contrastes independientes y PGLS', 'independent contrasts and PGLS'] },
    { art: 'mDEC', fam: 'bio', n: ['DEC', 'DEC'], s: ['dispersión, extinción, cladogénesis', 'dispersal, extinction, cladogenesis'] },
    { art: 'mDECJ', fam: 'bio', n: ['DEC+J y DIVALIKE', 'DEC+J and DIVALIKE'], s: ['con especiación fundadora', 'with founder-event speciation'] },
    { art: 'mStratified', fam: 'bio', n: ['Restricciones por tiempo', 'Time-stratified constraints'], s: ['qué áreas estaban conectadas', 'which areas were connected'] },
    { art: 'mBSM', fam: 'bio', n: ['Mapeo biogeográfico', 'Biogeographical mapping'], s: ['cuántos eventos y de qué tipo', 'how many events and of what kind'] },
    { art: 'mRF', fam: 'cmp', n: ['Distancias entre árboles', 'Tree distances'], s: ['Robinson–Foulds, cuartetos', 'Robinson–Foulds, quartets'] },
    { art: 'mConcordance', fam: 'cmp', n: ['Concordancia de genes', 'Gene concordance'], s: ['gCF y sCF por rama', 'gCF and sCF per branch'] },
    { art: 'mSpeciesTree', fam: 'cmp', n: ['Árbol de especies', 'Species tree'], s: ['resumen por cuartetos', 'quartet summary'] },
    { art: 'mNetwork', fam: 'cmp', n: ['Redes de divisiones', 'Split networks'], s: ['NeighborNet y descomposición', 'NeighborNet and decomposition'] },
    { art: 'mABBA', fam: 'cmp', n: ['ABBA-BABA', 'ABBA-BABA'], s: ['introgresión y captura de cloroplasto', 'introgression and chloroplast capture'] },
    { art: 'mTreeStudio', fam: 'cmp', n: ['Estudio de árboles', 'Tree studio'], s: ['la figura final, editable', 'the final figure, editable'] },
    { art: 'mExport', fam: 'cmp', n: ['Informe y exportación', 'Report and export'], s: ['Newick, NEXUS, PhyloXML, ZIP', 'Newick, NEXUS, PhyloXML, ZIP'] },
  ];

  /* ---------------- what it brings together ---------------- */
  const ICONS = {
    seq: '<svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round"><path d="M7 3c0 6 10 6 10 12s-10 6-10 6M17 3c0 6-10 6-10 12"/><path d="M8 7h8M8 17h8"/></svg>',
    tree: '<svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round"><path d="M3 12h4M7 5v14M7 5h6M7 19h5M13 2v6M13 2h7M13 8h7M12 15v6M12 15h8M12 21h8"/></svg>',
    scale: '<svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round"><circle cx="12" cy="12" r="9"/><path d="M12 7v5l3 2"/></svg>',
    globe: '<svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round"><circle cx="12" cy="12" r="9"/><path d="M3 12h18M12 3c3 3 3 15 0 18M12 3c-3 3-3 15 0 18"/></svg>',
    check: '<svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round"><path d="M12 3l8 4v5c0 5-3.5 8-8 9-4.5-1-8-4-8-9V7z"/><path d="M8.5 12l2.5 2.5 4.5-5"/></svg>',
    chat: '<svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round"><path d="M4 5h16v11H9l-5 4z"/><path d="M8 9h8M8 12h5"/></svg>',
    fig: '<svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round"><path d="M4 20V4M4 20h16"/><rect x="7" y="11" width="3" height="6"/><rect x="12" y="7" width="3" height="10"/><rect x="17" y="13" width="3" height="4"/></svg>',
    doc: '<svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round"><path d="M6 3h9l4 4v14H6z"/><path d="M14 3v5h5M9 13h7M9 17h5"/></svg>',
    lock: '<svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round"><rect x="5" y="11" width="14" height="10" rx="2"/><path d="M8 11V7a4 4 0 0 1 8 0v4"/></svg>',
    cpu: '<svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round"><rect x="7" y="7" width="10" height="10" rx="1.5"/><path d="M4 10h3M4 14h3M17 10h3M17 14h3M10 4v3M14 4v3M10 17v3M14 17v3"/></svg>',
  };
  const BRING = [
    ['seq', ['De las secuencias al artículo, sin cambiar de programa', 'From the sequences to the paper, without switching programs'],
      ['Alinear, elegir modelo, inferir, fechar y dibujar dejó de ser una cadena de siete programas y tres formatos.', 'Aligning, choosing a model, inferring, dating and drawing is no longer a chain of seven programs and three file formats.']],
    ['tree', ['Los cuatro criterios, con los mismos datos', 'All four criteria, on the same data'],
      ['Distancias, parsimonia, verosimilitud y bayesiano lado a lado, y la comparación explícita de en qué difieren.', 'Distances, parsimony, likelihood and Bayesian side by side, with an explicit comparison of where they differ.']],
    ['scale', ['El tiempo, no solo la forma', 'Time, not only shape'],
      ['Pruebas de reloj, calibraciones fósiles con su distribución e intervalos de credibilidad sobre una escala geológica.', 'Clock tests, fossil calibrations with their own distribution and credibility intervals on a geological scale.']],
    ['globe', ['Dónde estuvo cada linaje', 'Where each lineage was'],
      ['DEC y su familia sobre el árbol fechado: áreas ancestrales, dispersión, vicarianza y extinción.', 'DEC and its family on the dated tree: ancestral areas, dispersal, vicariance and extinction.']],
    ['check', ['Soporte que se puede defender', 'Support you can defend'],
      ['Bootstrap, bootstrap ultrarrápido, SH-aLRT, probabilidades posteriores y pruebas de topología, con su interpretación.', 'Bootstrap, ultrafast bootstrap, SH-aLRT, posterior probabilities and topology tests, each with its interpretation.']],
    ['chat', ['Interpretación en palabras', 'Interpretation in words'],
      ['Qué significa un soporte bajo, una rama larga sospechosa o un conflicto entre el cloroplasto y el núcleo.', 'What a low support, a suspicious long branch or a conflict between chloroplast and nucleus actually means.']],
    ['cpu', ['Aprovecha todos los núcleos de tu computadora', 'It uses every core of your computer'],
      ['El bootstrap y las cadenas calientes se reparten en hilos del navegador, y cada análisis largo dice cuánto va a tardar antes de empezar.', 'Bootstrap replicates and heated chains are spread across browser threads, and every long analysis says how long it will take before it starts.']],
    ['fig', ['Figuras listas para publicar', 'Publication-ready figures'],
      ['Estudio de árboles editable en colores, fuentes y etiquetas; hasta 900 ppp en PNG, TIFF o SVG.', 'A tree studio with editable colours, fonts and labels; up to 900 dpi as PNG, TIFF or SVG.']],
    ['doc', ['Informe reproducible', 'Reproducible report'],
      ['Metodología redactada con los ajustes reales, las citas de cada método y la semilla de cada corrida.', 'Methods written with the actual settings, the citation of every method and the seed of every run.']],
    ['lock', ['Tus datos no salen de tu computadora', 'Your data never leave your computer'],
      ['Todo se calcula en el navegador. Sin servidor, sin cola, sin conexión y sin subir nada.', 'Everything is computed in the browser. No server, no queue, no connection and nothing uploaded.']],
  ];

  /* ---------------- comparison of the inference criteria ---------------- */
  const Y = (es, en) => ({ c: 'yes', es, en }), N = (es, en) => ({ c: 'no', es, en }), O = (es, en) => ({ c: 'opt', es, en }), P = (es, en) => ({ c: '', es, en });
  const COMPARE_COLS = [
    ['Distancias (NJ)', 'Distances (NJ)'],
    ['Parsimonia', 'Parsimony'],
    ['Máxima verosimilitud', 'Maximum likelihood'],
    ['Bayesiano', 'Bayesian'],
  ];
  const COMPARE = [
    [['Qué optimiza', 'What it optimises'],
      P('Ajuste a una matriz de distancias', 'Fit to a distance matrix'),
      P('El menor número de cambios', 'The fewest changes'),
      P('La probabilidad de los datos dado el árbol', 'The probability of the data given the tree'),
      P('La probabilidad del árbol dados los datos', 'The probability of the tree given the data')],
    [['Usa un modelo de sustitución', 'Uses a substitution model'],
      O('Sí, al corregir las distancias', 'Yes, when correcting the distances'),
      N('No (o solo pesos)', 'No (or only weights)'),
      Y('Sí, explícito', 'Yes, explicit'),
      Y('Sí, explícito y con previas', 'Yes, explicit and with priors')],
    [['Consistente con ramas largas desiguales', 'Consistent with unequal long branches'],
      O('Depende de la corrección', 'Depends on the correction'),
      N('No: zona de Felsenstein', 'No: the Felsenstein zone'),
      Y('Sí, bajo el modelo correcto', 'Yes, under the correct model'),
      Y('Sí, bajo el modelo correcto', 'Yes, under the correct model')],
    [['Medida de incertidumbre', 'Measure of uncertainty'],
      P('Bootstrap', 'Bootstrap'),
      P('Bootstrap, jackknife, Bremer', 'Bootstrap, jackknife, Bremer'),
      P('Bootstrap, UFBoot, SH-aLRT', 'Bootstrap, UFBoot, SH-aLRT'),
      P('Probabilidad posterior', 'Posterior probability')],
    [['Estima longitudes de rama en sustituciones', 'Estimates branch lengths in substitutions'],
      Y('Sí', 'Yes'), O('En pasos, no en sustituciones', 'In steps, not substitutions'), Y('Sí', 'Yes'), Y('Sí, con intervalo', 'Yes, with an interval')],
    [['Permite fechar el árbol', 'Allows the tree to be dated'],
      O('Solo con métodos aproximados', 'Only with approximate methods'), N('No', 'No'), Y('Sí (reloj, LSD)', 'Yes (clock, LSD)'), Y('Sí, con intervalos de credibilidad', 'Yes, with credibility intervals')],
    [['Costo de cómputo', 'Computational cost'],
      P('Segundos', 'Seconds'), P('Segundos a minutos', 'Seconds to minutes'), P('Minutos', 'Minutes'), P('De minutos a horas', 'Minutes to hours')],
    [['Para qué sirve mejor', 'What it is best for'],
      P('Un primer vistazo y el árbol de arranque', 'A first look and the starting tree'),
      P('Morfología y datos con pocos cambios', 'Morphology and data with few changes'),
      P('El árbol que se publica', 'The tree that gets published'),
      P('Incertidumbre, fechas y modelos complejos', 'Uncertainty, dates and complex models')],
  ];

  /* ---------------- sizes the app really handles (measured, not guessed) ---------------- */
  const LIMITS = [
    [['Distancias, NJ, parsimonia', 'Distances, NJ, parsimony'], '≤ 300', '≤ 1 000', ['segundos', 'seconds']],
    [['Verosimilitud + bootstrap ultrarrápido', 'Likelihood + ultrafast bootstrap'], '≤ 120', '≤ 300', ['minutos', 'minutes']],
    [['Verosimilitud + bootstrap clásico', 'Likelihood + classical bootstrap'], '≤ 60', '≤ 120', ['minutos con varios hilos', 'minutes with several threads']],
    [['Bayesiano (MC³)', 'Bayesian (MC³)'], '≤ 50', '≤ 100', ['de minutos a horas', 'minutes to hours']],
    [['Datación por mínimos cuadrados', 'Least-squares dating'], '≤ 500', '≤ 1 000', ['segundos', 'seconds']],
    [['Datación bayesiana', 'Bayesian dating'], '≤ 60', '≤ 120', ['de minutos a horas', 'minutes to hours']],
    [['Biogeografía DEC', 'DEC biogeography'], ['≤ 8 áreas', '≤ 8 areas'], ['≤ 10 áreas', '≤ 10 areas'], ['minutos', 'minutes']],
  ];

  /* ---------------- references ---------------- */
  const REFS = [
    ['dat', 'Needleman, S.B. & Wunsch, C.D. (1970)', 'A general method applicable to the search for similarities in the amino acid sequence of two proteins. <i>Journal of Molecular Biology</i> 48: 443–453.'],
    ['dat', 'Gotoh, O. (1982)', 'An improved algorithm for matching biological sequences. <i>Journal of Molecular Biology</i> 162: 705–708.'],
    ['dat', 'Thompson, J.D., Higgins, D.G. & Gibson, T.J. (1994)', 'CLUSTAL W: improving the sensitivity of progressive multiple sequence alignment. <i>Nucleic Acids Research</i> 22: 4673–4680.'],
    ['dat', 'Edgar, R.C. (2004)', 'MUSCLE: multiple sequence alignment with high accuracy and high throughput. <i>Nucleic Acids Research</i> 32: 1792–1797.'],
    ['dat', 'Katoh, K. & Standley, D.M. (2013)', 'MAFFT multiple sequence alignment software version 7. <i>Molecular Biology and Evolution</i> 30: 772–780.'],
    ['dat', 'Castresana, J. (2000)', 'Selection of conserved blocks from multiple alignments for their use in phylogenetic analysis. <i>Molecular Biology and Evolution</i> 17: 540–552.'],
    ['dat', 'Capella-Gutiérrez, S., Silla-Martínez, J.M. & Gabaldón, T. (2009)', 'trimAl: a tool for automated alignment trimming. <i>Bioinformatics</i> 25: 1972–1973.'],
    ['mod', 'Jukes, T.H. & Cantor, C.R. (1969)', 'Evolution of protein molecules. In: Munro, H.N. (ed.) <i>Mammalian Protein Metabolism</i>. Academic Press, New York, pp. 21–132.'],
    ['mod', 'Kimura, M. (1980)', 'A simple method for estimating evolutionary rates of base substitutions through comparative studies of nucleotide sequences. <i>Journal of Molecular Evolution</i> 16: 111–120.'],
    ['mod', 'Felsenstein, J. (1981)', 'Evolutionary trees from DNA sequences: a maximum likelihood approach. <i>Journal of Molecular Evolution</i> 17: 368–376.'],
    ['mod', 'Hasegawa, M., Kishino, H. & Yano, T. (1985)', 'Dating of the human–ape splitting by a molecular clock of mitochondrial DNA. <i>Journal of Molecular Evolution</i> 22: 160–174.'],
    ['mod', 'Tavaré, S. (1986)', 'Some probabilistic and statistical problems in the analysis of DNA sequences. <i>Lectures on Mathematics in the Life Sciences</i> 17: 57–86.'],
    ['mod', 'Tamura, K. & Nei, M. (1993)', 'Estimation of the number of nucleotide substitutions in the control region of mitochondrial DNA in humans and chimpanzees. <i>Molecular Biology and Evolution</i> 10: 512–526.'],
    ['mod', 'Yang, Z. (1994)', 'Maximum likelihood phylogenetic estimation from DNA sequences with variable rates over sites: approximate methods. <i>Journal of Molecular Evolution</i> 39: 306–314.'],
    ['mod', 'Yang, Z. (1996)', 'Among-site rate variation and its impact on phylogenetic analyses. <i>Trends in Ecology & Evolution</i> 11: 367–372.'],
    ['mod', 'Lewis, P.O. (2001)', 'A likelihood approach to estimating phylogeny from discrete morphological character data. <i>Systematic Biology</i> 50: 913–925.'],
    ['mod', 'Le, S.Q. & Gascuel, O. (2008)', 'An improved general amino acid replacement matrix. <i>Molecular Biology and Evolution</i> 25: 1307–1320.'],
    ['mod', 'Whelan, S. & Goldman, N. (2001)', 'A general empirical model of protein evolution derived from multiple protein families using a maximum-likelihood approach. <i>Molecular Biology and Evolution</i> 18: 691–699.'],
    ['mod', 'Goldman, N. & Yang, Z. (1994)', 'A codon-based model of nucleotide substitution for protein-coding DNA sequences. <i>Molecular Biology and Evolution</i> 11: 725–736.'],
    ['mod', 'Posada, D. & Crandall, K.A. (1998)', 'MODELTEST: testing the model of DNA substitution. <i>Bioinformatics</i> 14: 817–818.'],
    ['mod', 'Kalyaanamoorthy, S., Minh, B.Q., Wong, T.K.F., von Haeseler, A. & Jermiin, L.S. (2017)', 'ModelFinder: fast model selection for accurate phylogenetic estimates. <i>Nature Methods</i> 14: 587–589.'],
    ['mod', 'Lanfear, R., Frandsen, P.B., Wright, A.M., Senfeld, T. & Calcott, B. (2017)', 'PartitionFinder 2: new methods for selecting partitioned models of evolution. <i>Molecular Biology and Evolution</i> 34: 772–773.'],
    ['mod', 'Xia, X., Xie, Z., Salemi, M., Chen, L. & Wang, Y. (2003)', 'An index of substitution saturation and its application. <i>Molecular Phylogenetics and Evolution</i> 26: 1–7.'],
    ['tre', 'Fitch, W.M. (1971)', 'Toward defining the course of evolution: minimum change for a specific tree topology. <i>Systematic Zoology</i> 20: 406–416.'],
    ['tre', 'Sankoff, D. (1975)', 'Minimal mutation trees of sequences. <i>SIAM Journal on Applied Mathematics</i> 28: 35–42.'],
    ['tre', 'Kluge, A.G. & Farris, J.S. (1969)', 'Quantitative phyletics and the evolution of anurans. <i>Systematic Zoology</i> 18: 1–32.'],
    ['tre', 'Farris, J.S. (1989)', 'The retention index and the rescaled consistency index. <i>Cladistics</i> 5: 417–419.'],
    ['tre', 'Bremer, K. (1994)', 'Branch support and tree stability. <i>Cladistics</i> 10: 295–304.'],
    ['tre', 'Saitou, N. & Nei, M. (1987)', 'The neighbor-joining method: a new method for reconstructing phylogenetic trees. <i>Molecular Biology and Evolution</i> 4: 406–425.'],
    ['tre', 'Studier, J.A. & Keppler, K.J. (1988)', 'A note on the neighbor-joining algorithm of Saitou and Nei. <i>Molecular Biology and Evolution</i> 5: 729–731.'],
    ['tre', 'Gascuel, O. (1997)', 'BIONJ: an improved version of the NJ algorithm based on a simple model of sequence data. <i>Molecular Biology and Evolution</i> 14: 685–695.'],
    ['tre', 'Desper, R. & Gascuel, O. (2002)', 'Fast and accurate phylogeny reconstruction algorithms based on the minimum-evolution principle. <i>Journal of Computational Biology</i> 9: 687–705.'],
    ['tre', 'Goloboff, P.A., Farris, J.S. & Nixon, K.C. (2008)', 'TNT, a free program for phylogenetic analysis. <i>Cladistics</i> 24: 774–786.'],
    ['inf', 'Felsenstein, J. (1978)', 'Cases in which parsimony or compatibility methods will be positively misleading. <i>Systematic Zoology</i> 27: 401–410.'],
    ['inf', 'Felsenstein, J. (1985)', 'Confidence limits on phylogenies: an approach using the bootstrap. <i>Evolution</i> 39: 783–791.'],
    ['inf', 'Hillis, D.M. & Bull, J.J. (1993)', 'An empirical test of bootstrapping as a method for assessing confidence in phylogenetic analysis. <i>Systematic Biology</i> 42: 182–192.'],
    ['inf', 'Kishino, H. & Hasegawa, M. (1989)', 'Evaluation of the maximum likelihood estimate of the evolutionary tree topologies from DNA sequence data. <i>Journal of Molecular Evolution</i> 29: 170–179.'],
    ['inf', 'Shimodaira, H. & Hasegawa, M. (1999)', 'Multiple comparisons of log-likelihoods with applications to phylogenetic inference. <i>Molecular Biology and Evolution</i> 16: 1114–1116.'],
    ['inf', 'Shimodaira, H. (2002)', 'An approximately unbiased test of phylogenetic tree selection. <i>Systematic Biology</i> 51: 492–508.'],
    ['inf', 'Guindon, S. & Gascuel, O. (2003)', 'A simple, fast, and accurate algorithm to estimate large phylogenies by maximum likelihood. <i>Systematic Biology</i> 52: 696–704.'],
    ['inf', 'Anisimova, M. & Gascuel, O. (2006)', 'Approximate likelihood-ratio test for branches: a fast, accurate, and powerful alternative. <i>Systematic Biology</i> 55: 539–552.'],
    ['inf', 'Guindon, S., Dufayard, J.-F., Lefort, V., Anisimova, M., Hordijk, W. & Gascuel, O. (2010)', 'New algorithms and methods to estimate maximum-likelihood phylogenies: assessing the performance of PhyML 3.0. <i>Systematic Biology</i> 59: 307–321.'],
    ['inf', 'Stamatakis, A. (2014)', 'RAxML version 8: a tool for phylogenetic analysis and post-analysis of large phylogenies. <i>Bioinformatics</i> 30: 1312–1313.'],
    ['inf', 'Nguyen, L.-T., Schmidt, H.A., von Haeseler, A. & Minh, B.Q. (2015)', 'IQ-TREE: a fast and effective stochastic algorithm for estimating maximum-likelihood phylogenies. <i>Molecular Biology and Evolution</i> 32: 268–274.'],
    ['inf', 'Hoang, D.T., Chernomor, O., von Haeseler, A., Minh, B.Q. & Vinh, L.S. (2018)', 'UFBoot2: improving the ultrafast bootstrap approximation. <i>Molecular Biology and Evolution</i> 35: 518–522.'],
    ['inf', 'Minh, B.Q., Schmidt, H.A., Chernomor, O., Schrempf, D., Woodhams, M.D., von Haeseler, A. & Lanfear, R. (2020)', 'IQ-TREE 2: new models and efficient methods for phylogenetic inference in the genomic era. <i>Molecular Biology and Evolution</i> 37: 1530–1534.'],
    ['inf', 'Kozlov, A.M., Darriba, D., Flouri, T., Morel, B. & Stamatakis, A. (2019)', 'RAxML-NG: a fast, scalable and user-friendly tool for maximum likelihood phylogenetic inference. <i>Bioinformatics</i> 35: 4453–4455.'],
    ['inf', 'Rannala, B. & Yang, Z. (1996)', 'Probability distribution of molecular evolutionary trees: a new method of phylogenetic inference. <i>Journal of Molecular Evolution</i> 43: 304–311.'],
    ['inf', 'Yang, Z. & Rannala, B. (1997)', 'Bayesian phylogenetic inference using DNA sequences: a Markov chain Monte Carlo method. <i>Molecular Biology and Evolution</i> 14: 717–724.'],
    ['inf', 'Mau, B., Newton, M.A. & Larget, B. (1999)', 'Bayesian phylogenetic inference via Markov chain Monte Carlo methods. <i>Biometrics</i> 55: 1–12.'],
    ['inf', 'Huelsenbeck, J.P. & Ronquist, F. (2001)', 'MRBAYES: Bayesian inference of phylogenetic trees. <i>Bioinformatics</i> 17: 754–755.'],
    ['inf', 'Ronquist, F., Teslenko, M., van der Mark, P., Ayres, D.L., Darling, A., Höhna, S., Larget, B., Liu, L., Suchard, M.A. & Huelsenbeck, J.P. (2012)', 'MrBayes 3.2: efficient Bayesian phylogenetic inference and model choice across a large model space. <i>Systematic Biology</i> 61: 539–542.'],
    ['inf', 'Gelman, A. & Rubin, D.B. (1992)', 'Inference from iterative simulation using multiple sequences. <i>Statistical Science</i> 7: 457–472.'],
    ['inf', 'Rambaut, A., Drummond, A.J., Xie, D., Baele, G. & Suchard, M.A. (2018)', 'Posterior summarization in Bayesian phylogenetics using Tracer 1.7. <i>Systematic Biology</i> 67: 901–904.'],
    ['tim', 'Zuckerkandl, E. & Pauling, L. (1965)', 'Evolutionary divergence and convergence in proteins. In: Bryson, V. & Vogel, H.J. (eds.) <i>Evolving Genes and Proteins</i>. Academic Press, New York, pp. 97–166.'],
    ['tim', 'Tajima, F. (1993)', 'Simple methods for testing the molecular evolutionary clock hypothesis. <i>Genetics</i> 135: 599–607.'],
    ['tim', 'Thorne, J.L., Kishino, H. & Painter, I.S. (1998)', 'Estimating the rate of evolution of the rate of molecular evolution. <i>Molecular Biology and Evolution</i> 15: 1647–1657.'],
    ['tim', 'Sanderson, M.J. (2002)', 'Estimating absolute rates of molecular evolution and divergence times: a penalized likelihood approach. <i>Molecular Biology and Evolution</i> 19: 101–109.'],
    ['tim', 'Drummond, A.J., Ho, S.Y.W., Phillips, M.J. & Rambaut, A. (2006)', 'Relaxed phylogenetics and dating with confidence. <i>PLoS Biology</i> 4: e88.'],
    ['tim', 'Ho, S.Y.W. & Phillips, M.J. (2009)', 'Accounting for calibration uncertainty in phylogenetic estimation of evolutionary divergence times. <i>Systematic Biology</i> 58: 367–380.'],
    ['tim', 'Tamura, K., Battistuzzi, F.U., Billing-Ross, P., Murillo, O., Filipski, A. & Kumar, S. (2012)', 'Estimating divergence times in large molecular phylogenies. <i>Proceedings of the National Academy of Sciences USA</i> 109: 19333–19338.'],
    ['tim', 'To, T.-H., Jung, M., Lycett, S. & Gascuel, O. (2016)', 'Fast dating using least-squares criteria and algorithms. <i>Systematic Biology</i> 65: 82–97.'],
    ['tim', 'Rambaut, A., Lam, T.T., Max Carvalho, L. & Pybus, O.G. (2016)', 'Exploring the temporal structure of heterochronous sequences using TempEst. <i>Virus Evolution</i> 2: vew007.'],
    ['tim', 'Bouckaert, R. et al. (2019)', 'BEAST 2.5: an advanced software platform for Bayesian evolutionary analysis. <i>PLoS Computational Biology</i> 15: e1006650.'],
    ['mac', 'Nee, S., May, R.M. & Harvey, P.H. (1994)', 'The reconstructed evolutionary process. <i>Philosophical Transactions of the Royal Society B</i> 344: 305–311.'],
    ['mac', 'Pybus, O.G. & Harvey, P.H. (2000)', 'Testing macro-evolutionary models using incomplete molecular phylogenies. <i>Proceedings of the Royal Society B</i> 267: 2267–2272.'],
    ['mac', 'Alfaro, M.E., Santini, F., Brock, C., Alamillo, H., Dornburg, A., Rabosky, D.L., Carnevale, G. & Harmon, L.J. (2009)', 'Nine exceptional radiations plus high turnover explain species diversity in jawed vertebrates. <i>Proceedings of the National Academy of Sciences USA</i> 106: 13410–13414.'],
    ['mac', 'Felsenstein, J. (1985)', 'Phylogenies and the comparative method. <i>The American Naturalist</i> 125: 1–15.'],
    ['mac', 'Grafen, A. (1989)', 'The phylogenetic regression. <i>Philosophical Transactions of the Royal Society B</i> 326: 119–157.'],
    ['mac', 'Pagel, M. (1994)', 'Detecting correlated evolution on phylogenies: a general method for the comparative analysis of discrete characters. <i>Proceedings of the Royal Society B</i> 255: 37–45.'],
    ['mac', 'Pagel, M. (1999)', 'Inferring the historical patterns of biological evolution. <i>Nature</i> 401: 877–884.'],
    ['mac', 'Hansen, T.F. (1997)', 'Stabilizing selection and the comparative analysis of adaptation. <i>Evolution</i> 51: 1341–1351.'],
    ['mac', 'Butler, M.A. & King, A.A. (2004)', 'Phylogenetic comparative analysis: a modeling approach for adaptive evolution. <i>The American Naturalist</i> 164: 683–695.'],
    ['mac', 'Blomberg, S.P., Garland, T. & Ives, A.R. (2003)', 'Testing for phylogenetic signal in comparative data: behavioral traits are more labile. <i>Evolution</i> 57: 717–745.'],
    ['mac', 'Huelsenbeck, J.P., Nielsen, R. & Bollback, J.P. (2003)', 'Stochastic mapping of morphological characters. <i>Systematic Biology</i> 52: 131–158.'],
    ['mac', 'Harmon, L.J., Weir, J.T., Brock, C.D., Glor, R.E. & Challenger, W. (2008)', 'GEIGER: investigating evolutionary radiations. <i>Bioinformatics</i> 24: 129–131.'],
    ['mac', 'Revell, L.J. (2012)', 'phytools: an R package for phylogenetic comparative biology (and other things). <i>Methods in Ecology and Evolution</i> 3: 217–223.'],
    ['bio', 'Ronquist, F. (1997)', 'Dispersal–vicariance analysis: a new approach to the quantification of historical biogeography. <i>Systematic Biology</i> 46: 195–203.'],
    ['bio', 'Ree, R.H. & Smith, S.A. (2008)', 'Maximum likelihood inference of geographic range evolution by dispersal, local extinction, and cladogenesis. <i>Systematic Biology</i> 57: 4–14.'],
    ['bio', 'Landis, M.J., Matzke, N.J., Moore, B.R. & Huelsenbeck, J.P. (2013)', 'Bayesian analysis of biogeography when the number of areas is large. <i>Systematic Biology</i> 62: 789–804.'],
    ['bio', 'Matzke, N.J. (2014)', 'Model selection in historical biogeography reveals that founder-event speciation is a crucial process in island clades. <i>Systematic Biology</i> 63: 951–970.'],
    ['bio', 'Ree, R.H. & Sanmartín, I. (2018)', 'Conceptual and statistical problems with the DEC+J model of founder-event speciation and its comparison with DEC via model selection. <i>Journal of Biogeography</i> 45: 741–749.'],
    ['cmp', 'Robinson, D.F. & Foulds, L.R. (1981)', 'Comparison of phylogenetic trees. <i>Mathematical Biosciences</i> 53: 131–147.'],
    ['cmp', 'Bryant, D. & Moulton, V. (2004)', 'Neighbor-Net: an agglomerative method for the construction of phylogenetic networks. <i>Molecular Biology and Evolution</i> 21: 255–265.'],
    ['cmp', 'Huson, D.H. & Bryant, D. (2006)', 'Application of phylogenetic networks in evolutionary studies. <i>Molecular Biology and Evolution</i> 23: 254–267.'],
    ['cmp', 'Durand, E.Y., Patterson, N., Reich, D. & Slatkin, M. (2011)', 'Testing for ancient admixture between closely related populations. <i>Molecular Biology and Evolution</i> 28: 2239–2252.'],
    ['cmp', 'Mirarab, S. & Warnow, T. (2015)', 'ASTRAL-II: coalescent-based species tree estimation with many hundreds of taxa and thousands of genes. <i>Bioinformatics</i> 31: i44–i52.'],
    ['cmp', 'Minh, B.Q., Hahn, M.W. & Lanfear, R. (2020)', 'New methods to calculate concordance factors for phylogenomic datasets. <i>Molecular Biology and Evolution</i> 37: 2727–2733.'],
    ['cmp', 'Paradis, E. & Schliep, K. (2019)', 'ape 5.0: an environment for modern phylogenetics and evolutionary analyses in R. <i>Bioinformatics</i> 35: 526–528.'],
    ['cmp', 'Schliep, K.P. (2011)', 'phangorn: phylogenetic analysis in R. <i>Bioinformatics</i> 27: 592–593.'],
    ['cmp', 'Tamura, K., Stecher, G. & Kumar, S. (2021)', 'MEGA11: Molecular Evolutionary Genetics Analysis version 11. <i>Molecular Biology and Evolution</i> 38: 3022–3027.'],
    ['cmp', 'Felsenstein, J. (2004)', '<i>Inferring Phylogenies</i>. Sinauer Associates, Sunderland.'],
    ['cmp', 'Yang, Z. (2014)', '<i>Molecular Evolution: A Statistical Approach</i>. Oxford University Press, Oxford.'],
    ['cmp', 'Nei, M. & Kumar, S. (2000)', '<i>Molecular Evolution and Phylogenetics</i>. Oxford University Press, New York.'],
  ];

  /* ---------------- renderers ---------------- */
  const two = pair => L2(pair[0], pair[1]);
  const tw = v => Array.isArray(v) ? L2(v[0], v[1]) : v;

  function renderStepper() {
    const nav = el('stepper');
    if (!nav) return;
    nav.innerHTML = STEPS.map(s => `<button class="step-btn${s.n === 1 ? ' active' : ''}" data-step="${s.n}"${s.ready ? '' : ' disabled'}><span class="step-num">${s.n}</span>${L2(s.es, s.en)}</button>`).join('');
  }

  function renderFeatures() {
    const g = el('featureGrid');
    if (!g) return;
    g.innerHTML = '';
    BLOCKS.forEach(b => {
      const card = mk('div', { class: 'feature', tabindex: '0', role: 'button' });
      const btn = document.querySelector(`.step-btn[data-step="${b.n}"]`);
      const ready = btn && !btn.disabled;
      card.innerHTML = `<div class="f-num">${b.n}</div>` +
        `<div class="f-art">${Art[b.art] ? Art[b.art]() : ''}</div>` +
        `<div class="f-tag">${two(b.tag)}${ready ? '' : ` <span class="f-soon">${L2('en construcción', 'coming next')}</span>`}</div>` +
        `<h3>${two(b.t)}</h3><p>${two(b.d)}</p>`;
      const open = () => {
        const b2 = document.querySelector('.step-btn[data-step="' + b.n + '"]');
        if (b2 && !b2.disabled) goStep(b.n);
        else {
          const m = el('homeMessages');
          if (m) {
            clearMessages(m);
            showMessage(m, 'info', L2(`El bloque ${b.n} se construye en una etapa posterior; la portada ya muestra lo que hará.`,
              `Block ${b.n} is built in a later stage; the home page already shows what it will do.`));
            m.scrollIntoView({ behavior: 'smooth', block: 'center' });
          }
        }
      };
      card.addEventListener('click', open);
      card.addEventListener('keydown', e => { if (e.key === 'Enter' || e.key === ' ') { e.preventDefault(); open(); } });
      g.appendChild(card);
    });
  }

  function renderMaterials() {
    const g = el('materialStrip');
    if (!g) return;
    g.innerHTML = MATERIALS.map(m => `<div class="material">${Art[m.art]()}<div class="mt-t">${two(m.t)}</div><div class="mt-s">${two(m.s)}</div><span class="mt-k ${m.kc}">${two(m.k)}</span></div>`).join('');
  }

  function renderMarkers() {
    const g = el('markerGrid');
    if (!g) return;
    g.innerHTML = MARKERS.map(m => `<div class="marker">
      <span class="genome-tag ${m.g}">${two(GENOME_NAMES[m.g])}</span>
      <div class="mk-n${/[ ·]/.test(m.n) && m.g === 'nu' ? '' : ''}">${m.n}</div>
      <div class="mk-d">${two(m.d)}</div></div>`).join('');
  }

  let famFilter = 'all';
  function renderMethods() {
    const g = el('methodGallery'), f = el('methodFilter');
    if (!g) return;
    if (f && !f.dataset.built) {
      f.dataset.built = '1';
      f.innerHTML = `<button class="chip on" data-fam="all">${L2('todos', 'all')} · ${METHODS.length}</button>` +
        Object.keys(FAMS).map(k => `<button class="chip" data-fam="${k}">${two(FAMS[k])} · ${METHODS.filter(m => m.fam === k).length}</button>`).join('');
      f.addEventListener('click', e => {
        const b = e.target.closest('.chip'); if (!b) return;
        famFilter = b.dataset.fam;
        els('.chip', f).forEach(c => c.classList.toggle('on', c === b));
        renderMethods();
      });
    }
    g.innerHTML = METHODS.filter(m => famFilter === 'all' || m.fam === famFilter).map(m =>
      `<div class="method-card"><span class="m-fam ${m.fam}">${two(FAMS[m.fam])}</span>${Art[m.art] ? Art[m.art]() : ''}<div class="m-name">${two(m.n)}</div><div class="m-sub">${two(m.s)}</div></div>`).join('');
  }

  function renderBring() {
    const g = el('bringGrid');
    if (!g) return;
    g.innerHTML = BRING.map(b => `<div class="bring"><div class="b-ic">${ICONS[b[0]]}</div><div><b>${two(b[1])}</b><span>${two(b[2])}</span></div></div>`).join('');
  }

  function renderCompare() {
    const wrap = el('methodCompare');
    if (!wrap) return;
    let html = '<table><thead><tr><th></th>' + COMPARE_COLS.map(c => `<th>${two(c)}</th>`).join('') + '</tr></thead><tbody>';
    COMPARE.forEach(row => {
      html += `<tr><td>${two(row[0])}</td>` + row.slice(1).map(c =>
        `<td>${c.c ? `<span class="${c.c}">${c.c === 'yes' ? '●' : c.c === 'no' ? '○' : '◐'}</span> ` : ''}${L2(c.es, c.en)}</td>`).join('') + '</tr>';
    });
    wrap.innerHTML = html + '</tbody></table>';
  }

  function renderLimits() {
    const wrap = el('limitTable');
    if (!wrap) return;
    let html = `<table class="limit-table"><thead><tr>
      <th>${L2('Análisis', 'Analysis')}</th>
      <th class="num">${L2('Cómodo (taxones)', 'Comfortable (taxa)')}</th>
      <th class="num">${L2('Con aviso', 'With a warning')}</th>
      <th>${L2('Tiempo típico', 'Typical time')}</th></tr></thead><tbody>`;
    LIMITS.forEach(r => {
      html += `<tr><td>${tw(r[0])}</td><td class="num">${tw(r[1])}</td><td class="num">${tw(r[2])}</td><td>${tw(r[3])}</td></tr>`;
    });
    wrap.innerHTML = html + '</tbody></table>';
  }

  let refFilter = 'all';
  function renderRefs() {
    const g = el('refList'), f = el('refFilter');
    if (!g) return;
    if (f && !f.dataset.built) {
      f.dataset.built = '1';
      f.innerHTML = `<button class="chip on" data-fam="all">${L2('todas', 'all')} · ${REFS.length}</button>` +
        Object.keys(FAMS).map(k => `<button class="chip" data-fam="${k}">${two(FAMS[k])}</button>`).join('');
      f.addEventListener('click', e => {
        const b = e.target.closest('.chip'); if (!b) return;
        refFilter = b.dataset.fam;
        els('.chip', f).forEach(c => c.classList.toggle('on', c === b));
        renderRefs();
      });
    }
    const list = REFS.filter(r => refFilter === 'all' || r[0] === refFilter);
    g.innerHTML = list.length
      ? list.map(r => `<li><b>${r[1]}</b> ${r[2]}</li>`).join('')
      : `<li>${L2('Las referencias de esta familia llegan con su bloque.', 'The references of this family arrive with its block.')}</li>`;
  }

  /* Blocks 2 to 12 exist as panels from the start, so the stepper, the block
     cards and the report all have somewhere to point. Each one shows what it
     will contain until it is built. */
  function renderPlaceholders() {
    const host = el('placeholders');
    if (!host) return;
    host.innerHTML = BLOCKS.map(b => {
      const step = STEPS.find(s => s.n === b.n) || {};
      if (step.ready) return '';
      return `<section class="step-panel" id="panel-${b.n}">
        <div class="panel-title"><div class="pt-num">${b.n}</div>
          <div><h2>${two(b.t)}</h2><p>${two(b.tag)}</p></div></div>
        <div class="card"><div class="soon-box">
          ${Art[b.art] ? Art[b.art]() : ''}
          <h3>${L2('Este bloque se está construyendo', 'This block is being built')}</h3>
          <p style="max-width:720px;margin:6px auto 0">${two(b.d)}</p>
          <div class="btn-row" style="justify-content:center">
            <button class="btn btn-secondary btn-sm" data-goto="1">${L2('← Volver al inicio', '← Back to the home page')}</button>
          </div>
        </div></div></section>`;
    }).join('');
    host.addEventListener('click', e => {
      const b = e.target.closest('[data-goto]');
      if (b) goStep(b.dataset.goto);
    });
  }

  /* illustrations that carry translated labels are redrawn when the language changes */
  function renderArt() {
    const h = el('heroArt'); if (h) h.innerHTML = Art.hero();
    const figs = {
      theoryAnatomyFig: 'theoryAnatomy', theorySatFig: 'theorySaturation', theoryLbaFig: 'theoryLongBranch',
      theoryClockFig: 'theoryClock', theoryGenomeFig: 'theoryGenomes',
    };
    for (const id in figs) {
      const n = el(id);
      if (n) { const cap = n.querySelector('.cap'); n.innerHTML = Art[figs[id]](); if (cap) n.appendChild(cap); }
    }
    const b = el('brandLogo');
    if (b) b.innerHTML = Art.brand();
  }

  /* ---------------- navigation ---------------- */
  function wire() {
    const nav = el('stepper');
    if (nav) nav.addEventListener('click', e => { const b = e.target.closest('.step-btn'); if (b && !b.disabled) goStep(b.dataset.step); });
    const brand = el('brand');
    if (brand) brand.addEventListener('click', () => goStep(1));
    const scrollTo = id => { const n = el(id); if (n) n.scrollIntoView({ behavior: 'smooth', block: 'start' }); };
    const on = (id, fn) => { const n = el(id); if (n) n.addEventListener('click', fn); };
    on('startBtn', () => {
      const b = document.querySelector('.step-btn[data-step="2"]');
      if (b && !b.disabled) goStep(2);
      else {
        const m = el('homeMessages');
        if (m) {
          clearMessages(m);
          showMessage(m, 'info', L2('La carga de datos y el alineador llegan con el Bloque 2. Mientras tanto, prueba los laboratorios y la teoría de esta página.',
            'Data import and the aligner arrive with Block 2. Meanwhile, try the labs and the theory on this page.'));
          m.scrollIntoView({ behavior: 'smooth', block: 'center' });
        }
      }
    });
    on('simBtn', () => scrollTo('labs'));
    on('theoryBtn', () => { scrollTo('theory'); const first = document.querySelector('#theory .acc'); if (first) first.open = true; });
    on('citeBtn', () => scrollTo('cite'));
    on('copyCite', () => {
      const t = el('citeText');
      if (!t || !navigator.clipboard) return;
      const txt = [...t.querySelectorAll('[data-l="' + I18N.lang + '"]')].map(n => n.textContent).join('') || t.textContent;
      navigator.clipboard.writeText(txt.trim()).then(() => {
        const b = el('copyCite');
        if (b) { b.textContent = T('✓ Copiada', '✓ Copied'); setTimeout(() => I18N.apply(b.parentNode), 1800); }
      });
    });
    document.addEventListener('langchange', () => { renderArt(); renderMethods(); renderRefs(); });
    document.addEventListener('themechange', renderArt);
  }

  function init() {
    renderStepper();
    renderPlaceholders();
    renderArt();
    renderFeatures();
    renderMaterials();
    renderMarkers();
    renderMethods();
    renderBring();
    renderCompare();
    renderLimits();
    renderRefs();
    wire();
    I18N.apply();
  }

  document.addEventListener('DOMContentLoaded', init);
  window.Home = { BLOCKS, METHODS, REFS, MATERIALS, MARKERS, FAMS, LIMITS, COMPARE };
})();
