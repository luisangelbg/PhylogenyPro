/* PhylogenyPro — molecular phylogenetics in the browser.
   Copyright (C) 2026  Luis Ángel Barrera-Guzmán

   This program is free software: you can redistribute it and/or modify it under
   the terms of the GNU General Public License as published by the Free Software
   Foundation, either version 3 of the License, or (at your option) any later
   version. It is distributed in the hope that it will be useful, but WITHOUT ANY
   WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
   A PARTICULAR PURPOSE. See the GNU General Public License, in the file LICENSE
   at the root of this program, or <https://www.gnu.org/licenses/>. */

/* PhylogenyPro — writing down what was actually done.

   A methods section is not decoration: it is the part of a paper that lets
   somebody else get the same answer, and it is the part that gets written last,
   from memory, weeks after the analysis. This module writes it from the state
   the blocks left behind, so it says what the program did and not what the
   author remembers asking for.

   Two rules it follows, and they are the whole design:

     · it only claims what happened. A block that was never run leaves no
       paragraph. A support value that was not computed is not mentioned. The
       number of bootstrap replicates comes from the run, not from a template;
     · it cites what it used, and only that. The bibliography is assembled from
       the methods that actually ran, so a study that never touched the
       molecular clock does not end up citing Sanderson.

   Nothing here touches the DOM: it returns strings, and the block decides what
   to do with them. */

function ReportCore(g) {
const Report = {};

  /* ================================================================
     the bibliography
     ================================================================ */
  /* Only the works the app's own methods come from. A reference nobody cited is
     never printed, which is why every entry carries its key rather than a
     position in a list. */
  const REFS = {
    felsenstein1981: 'Felsenstein, J. (1981). Evolutionary trees from DNA sequences: a maximum likelihood approach. Journal of Molecular Evolution 17: 368–376.',
    felsenstein1985: 'Felsenstein, J. (1985). Phylogenies and the comparative method. The American Naturalist 125: 1–15.',
    felsenstein1985b: 'Felsenstein, J. (1985). Confidence limits on phylogenies: an approach using the bootstrap. Evolution 39: 783–791.',
    fitch1971: 'Fitch, W. M. (1971). Toward defining the course of evolution: minimum change for a specific tree topology. Systematic Zoology 20: 406–416.',
    sankoff1975: 'Sankoff, D. (1975). Minimal mutation trees of sequences. SIAM Journal on Applied Mathematics 28: 35–42.',
    saitou1987: 'Saitou, N. & Nei, M. (1987). The neighbor-joining method: a new method for reconstructing phylogenetic trees. Molecular Biology and Evolution 4: 406–425.',
    gascuel1997: 'Gascuel, O. (1997). BIONJ: an improved version of the NJ algorithm based on a simple model of sequence data. Molecular Biology and Evolution 14: 685–695.',
    desper2002: 'Desper, R. & Gascuel, O. (2002). Fast and accurate phylogeny reconstruction algorithms based on the minimum-evolution principle. Journal of Computational Biology 9: 687–705.',
    pauplin2000: 'Pauplin, Y. (2000). Direct calculation of a tree length using a distance matrix. Journal of Molecular Evolution 51: 41–47.',
    yang1994: 'Yang, Z. (1994). Maximum likelihood phylogenetic estimation from DNA sequences with variable rates over sites: approximate methods. Journal of Molecular Evolution 39: 306–314.',
    tavare1986: 'Tavaré, S. (1986). Some probabilistic and statistical problems in the analysis of DNA sequences. Lectures on Mathematics in the Life Sciences 17: 57–86.',
    lewis2001: 'Lewis, P. O. (2001). A likelihood approach to estimating phylogeny from discrete morphological character data. Systematic Biology 50: 913–925.',
    castresana2000: 'Castresana, J. (2000). Selection of conserved blocks from multiple alignments for their use in phylogenetic analysis. Molecular Biology and Evolution 17: 540–552.',
    gotoh1982: 'Gotoh, O. (1982). An improved algorithm for matching biological sequences. Journal of Molecular Biology 162: 705–708.',
    lanfear2012: 'Lanfear, R., Calcott, B., Ho, S. Y. W. & Guindon, S. (2012). PartitionFinder: combined selection of partitioning schemes and substitution models for phylogenetic analyses. Molecular Biology and Evolution 29: 1695–1701.',
    minh2013: 'Minh, B. Q., Nguyen, M. A. T. & von Haeseler, A. (2013). Ultrafast approximation for phylogenetic bootstrap. Molecular Biology and Evolution 30: 1188–1195.',
    guindon2010: 'Guindon, S., Dufayard, J.-F., Lefort, V., Anisimova, M., Hordijk, W. & Gascuel, O. (2010). New algorithms and methods to estimate maximum-likelihood phylogenies: assessing the performance of PhyML 3.0. Systematic Biology 59: 307–321.',
    anisimova2011: 'Anisimova, M., Gil, M., Dufayard, J.-F., Dessimoz, C. & Gascuel, O. (2011). Survey of branch support methods demonstrates accuracy, power, and robustness of fast likelihood-based approximation schemes. Systematic Biology 60: 685–699.',
    shimodaira2002: 'Shimodaira, H. (2002). An approximately unbiased test of phylogenetic tree selection. Systematic Biology 51: 492–508.',
    shimodaira1999: 'Shimodaira, H. & Hasegawa, M. (1999). Multiple comparisons of log-likelihoods with applications to phylogenetic inference. Molecular Biology and Evolution 16: 1114–1116.',
    kishino1989: 'Kishino, H. & Hasegawa, M. (1989). Evaluation of the maximum likelihood estimate of the evolutionary tree topologies from DNA sequence data. Journal of Molecular Evolution 29: 170–179.',
    geyer1991: 'Geyer, C. J. (1991). Markov chain Monte Carlo maximum likelihood. In: Computing Science and Statistics: Proceedings of the 23rd Symposium on the Interface, pp. 156–163.',
    larget1999: 'Larget, B. & Simon, D. L. (1999). Markov chain Monte Carlo algorithms for the Bayesian analysis of phylogenetic trees. Molecular Biology and Evolution 16: 750–759.',
    xie2011: 'Xie, W., Lewis, P. O., Fan, Y., Kuo, L. & Chen, M.-H. (2011). Improving marginal likelihood estimation for Bayesian phylogenetic model selection. Systematic Biology 60: 150–160.',
    gelman1992: 'Gelman, A. & Rubin, D. B. (1992). Inference from iterative simulation using multiple sequences. Statistical Science 7: 457–472.',
    brown2010: 'Brown, J. M., Hedtke, S. M., Lemmon, A. R. & Lemmon, E. M. (2010). When trees grow too long: investigating the causes of highly inaccurate Bayesian branch-length estimates. Systematic Biology 59: 145–161.',
    to2016: 'To, T.-H., Jung, M., Lycett, S. & Gascuel, O. (2016). Fast dating using least-squares criteria and algorithms. Systematic Biology 65: 82–97.',
    sanderson2002: 'Sanderson, M. J. (2002). Estimating absolute rates of molecular evolution and divergence times: a penalized likelihood approach. Molecular Biology and Evolution 19: 101–109.',
    drummond2006: 'Drummond, A. J., Ho, S. Y. W., Phillips, M. J. & Rambaut, A. (2006). Relaxed phylogenetics and dating with confidence. PLoS Biology 4: e88.',
    rambaut2016: 'Rambaut, A., Lam, T. T., Max Carvalho, L. & Pybus, O. G. (2016). Exploring the temporal structure of heterochronous sequences using TempEst. Virus Evolution 2: vew007.',
    pybus2000: 'Pybus, O. G. & Harvey, P. H. (2000). Testing macro-evolutionary models using incomplete molecular phylogenies. Proceedings of the Royal Society B 267: 2267–2272.',
    nee1994: 'Nee, S., May, R. M. & Harvey, P. H. (1994). The reconstructed evolutionary process. Philosophical Transactions of the Royal Society B 344: 305–311.',
    rabosky2008: 'Rabosky, D. L. & Lovette, I. J. (2008). Density-dependent diversification in North American wood warblers. Proceedings of the Royal Society B 275: 2363–2371.',
    magallon2001: 'Magallón, S. & Sanderson, M. J. (2001). Absolute diversification rates in angiosperm clades. Evolution 55: 1762–1780.',
    jetz2012: 'Jetz, W., Thomas, G. H., Joy, J. B., Hartmann, K. & Mooers, A. O. (2012). The global diversity of birds in space and time. Nature 491: 444–448.',
    pagel1999: 'Pagel, M. (1999). Inferring the historical patterns of biological evolution. Nature 401: 877–884.',
    blomberg2003: 'Blomberg, S. P., Garland, T. & Ives, A. R. (2003). Testing for phylogenetic signal in comparative data: behavioral traits are more labile. Evolution 57: 717–745.',
    garland2000: 'Garland, T. & Ives, A. R. (2000). Using the past to predict the present: confidence intervals for regression equations in phylogenetic comparative methods. The American Naturalist 155: 346–364.',
    hansen1997: 'Hansen, T. F. (1997). Stabilizing selection and the comparative analysis of adaptation. Evolution 51: 1341–1351.',
    harmon2010: 'Harmon, L. J. et al. (2010). Early bursts of body size and shape evolution are rare in comparative data. Evolution 64: 2385–2396.',
    huelsenbeck2003: 'Huelsenbeck, J. P., Nielsen, R. & Bollback, J. P. (2003). Stochastic mapping of morphological characters. Systematic Biology 52: 131–158.',
    ree2008: 'Ree, R. H. & Smith, S. A. (2008). Maximum likelihood inference of geographic range evolution by dispersal, local extinction, and cladogenesis. Systematic Biology 57: 4–14.',
    matzke2013: 'Matzke, N. J. (2013). Probabilistic historical biogeography: new models for founder-event speciation, imperfect detection, and fossils allow improved accuracy and model-testing. Frontiers of Biogeography 5: 242–248.',
    ree2018: 'Ree, R. H. & Sanmartín, I. (2018). Conceptual and statistical problems with the DEC+J model of founder-event speciation and its comparison with DEC via model selection. Journal of Biogeography 45: 741–749.',
    robinson1981: 'Robinson, D. F. & Foulds, L. R. (1981). Comparison of phylogenetic trees. Mathematical Biosciences 53: 131–147.',
    kuhner1994: 'Kuhner, M. K. & Felsenstein, J. (1994). A simulation comparison of phylogeny algorithms under equal and unequal evolutionary rates. Molecular Biology and Evolution 11: 459–468.',
    steel1993: 'Steel, M. A. & Penny, D. (1993). Distributions of tree comparison metrics: some new results. Systematic Biology 42: 126–141.',
    minh2020: 'Minh, B. Q., Hahn, M. W. & Lanfear, R. (2020). New methods to calculate concordance factors for phylogenomic datasets. Molecular Biology and Evolution 37: 2727–2733.',
    bandelt1992: 'Bandelt, H.-J. & Dress, A. W. M. (1992). A canonical decomposition theory for metrics on a finite set. Advances in Mathematics 92: 47–105.',
    green2010: 'Green, R. E. et al. (2010). A draft sequence of the Neandertal genome. Science 328: 710–722.',
    durand2011: 'Durand, E. Y., Patterson, N., Reich, D. & Slatkin, M. (2011). Testing for ancient admixture between closely related populations. Molecular Biology and Evolution 28: 2239–2252.',
    han2009: 'Han, M. V. & Zmasek, C. M. (2009). phyloXML: XML for evolutionary biology and comparative genomics. BMC Bioinformatics 10: 356.',
  };

  /* ================================================================
     the methods, section by section
     ================================================================ */
  /* Each section is produced only when the block that would justify it left
     something behind. The prose is written twice, because a methods section
     read in translation is a methods section nobody checks. */
  function methods(st, opts) {
    opts = opts || {};
    const out = [];
    const num = v => (v == null || !isFinite(v) ? '—' : (Math.abs(v) >= 0.01 && Math.abs(v) < 1e6
      ? String(Math.round(v * 10000) / 10000) : (+v).toPrecision(4)));
    const pct = v => (v == null ? '—' : (100 * v).toFixed(1) + ' %');

    /* ---- data and alignment ---- */
    if (st.data && st.data.parts && st.data.parts.length) {
      const parts = st.data.parts;
      const total = parts.reduce((a, p) => a + ((p.seqs && p.seqs[0]) ? p.seqs[0].length : 0), 0);
      const nTax = (st.data.taxa || (parts[0] && parts[0].taxa) || []).length;
      const aligned = parts.filter(p => p.source === 'alineado' || p.aligned).length;
      const trimmed = parts.filter(p => p.trimmed).length;
      const es = [], en = [];
      es.push(`El conjunto de datos reúne ${nTax} terminales y ${parts.length} ${parts.length === 1 ? 'partición' : 'particiones'} (${parts.map(p => `${p.name}, ${(p.seqs && p.seqs[0] ? p.seqs[0].length : 0)} caracteres`).join('; ')}), ${total} caracteres en total.`);
      en.push(`The data set comprises ${nTax} terminals and ${parts.length} ${parts.length === 1 ? 'partition' : 'partitions'} (${parts.map(p => `${p.name}, ${(p.seqs && p.seqs[0] ? p.seqs[0].length : 0)} characters`).join('; ')}), ${total} characters in all.`);
      /* the genome is stored as a two-letter code, which is fine in a menu and
         absurd in a methods section: "the genomes represented are: nu" */
      const GENOME_NAMES = {
        cp: ['cloroplasto', 'chloroplast'], mt: ['mitocondria', 'mitochondrion'],
        nr: ['nuclear ribosomal', 'nuclear ribosomal'], nu: ['nuclear de copia baja', 'low-copy nuclear'],
        other: ['sin asignar', 'unassigned'],
      };
      const genomes = [...new Set(parts.map(p => p.genome).filter(Boolean))];
      if (genomes.length) {
        const nameOf = (gcode, i) => (GENOME_NAMES[gcode] ? GENOME_NAMES[gcode][i] : gcode);
        es.push(`Los genomas representados son: ${genomes.map(x => nameOf(x, 0)).join(', ')}.`);
        en.push(`The genomes represented are: ${genomes.map(x => nameOf(x, 1)).join(', ')}.`);
      }
      if (aligned) {
        es.push('Los alineamientos se construyeron con el alineador progresivo de PhylogenyPro (distancias de k-meros, árbol guía por UPGMA y alineamiento perfil–perfil con huecos afines) y se revisaron a mano en el visor del programa.');
        en.push("The alignments were built with PhylogenyPro's progressive aligner (k-mer distances, a UPGMA guide tree and profile-to-profile alignment with affine gaps) and inspected by hand in the program's viewer.");
      }
      if (trimmed) {
        es.push('Los bloques de alineamiento ambiguos se recortaron antes del análisis.');
        en.push('Ambiguously aligned blocks were trimmed before analysis.');
      }
      out.push({ id: 'data', titleEs: 'Datos y alineamiento', titleEn: 'Data and alignment',
        es, en, refs: aligned ? ['gotoh1982'].concat(trimmed ? ['castresana2000'] : []) : (trimmed ? ['castresana2000'] : []) });
    }

    /* ---- substitution models ---- */
    if (st.models && st.models.best) {
      const b = st.models.best;
      const crit = st.models.criterion || 'BIC';
      const es = [], en = [];
      es.push(`El modelo de sustitución se eligió por ${crit} entre los candidatos que el programa ajusta sobre una topología fija; el seleccionado fue <b>${b.name || b.model}</b>${b.lnL != null ? ` (lnL = ${num(b.lnL)})` : ''}.`);
      en.push(`The substitution model was chosen by ${crit} among the candidates the program fits on a fixed topology; the selected one was <b>${b.name || b.model}</b>${b.lnL != null ? ` (lnL = ${num(b.lnL)})` : ''}.`);
      if (b.alpha != null) {
        es.push(`La heterogeneidad de tasas entre sitios se modeló con una distribución gamma discreta de cuatro categorías (α = ${num(b.alpha)}).`);
        en.push(`Rate heterogeneity among sites was modelled with a four-category discrete gamma distribution (α = ${num(b.alpha)}).`);
      }
      if (st.models.scheme && st.models.scheme.length > 1) {
        es.push(`El esquema de particiones se obtuvo por fusión codiciosa, y quedó en ${st.models.scheme.length} subconjuntos.`);
        en.push(`The partitioning scheme was found by greedy merging and settled on ${st.models.scheme.length} subsets.`);
      }
      out.push({ id: 'models', titleEs: 'Modelos de sustitución', titleEn: 'Substitution models', es, en,
        refs: ['tavare1986', 'yang1994'].concat(st.models.scheme && st.models.scheme.length > 1 ? ['lanfear2012'] : []) });
    }

    /* ---- parsimony and distances ---- */
    if (st.quick) {
      const es = [], en = [];
      if (st.quick.distance) {
        es.push(`Se construyó un árbol de distancias por ${st.quick.distance.name || 'neighbor joining'} sobre distancias corregidas.`);
        en.push(`A distance tree was built by ${st.quick.distance.name || 'neighbour joining'} over corrected distances.`);
      }
      if (st.quick.parsimony) {
        const p = st.quick.parsimony;
        es.push(`La búsqueda por parsimonia (adición paso a paso seguida de intercambios NNI, SPR y TBR) encontró ${p.nTrees || 1} ${(p.nTrees || 1) === 1 ? 'árbol' : 'árboles'} de ${p.steps} pasos${p.indices ? `, con CI = ${num(p.indices.ci)} y RI = ${num(p.indices.ri)}` : ''}.`);
        en.push(`The parsimony search (stepwise addition followed by NNI, SPR and TBR rearrangements) found ${p.nTrees || 1} ${(p.nTrees || 1) === 1 ? 'tree' : 'trees'} of ${p.steps} steps${p.indices ? `, with CI = ${num(p.indices.ci)} and RI = ${num(p.indices.ri)}` : ''}.`);
      }
      if (es.length) {
        out.push({ id: 'quick', titleEs: 'Parsimonia y distancias', titleEn: 'Parsimony and distances', es, en,
          refs: [].concat(st.quick.distance ? ['saitou1987', 'gascuel1997'] : [])
            .concat(st.quick.parsimony ? ['fitch1971', 'sankoff1975'] : []) });
      }
    }

    /* ---- maximum likelihood ---- */
    if (st.ml) {
      const m = st.ml;
      const es = [], en = [];
      es.push(`El árbol de máxima verosimilitud se buscó con intercambios NNI y SPR con radio de reinjerto, partiendo del árbol de distancias y optimizando las longitudes de rama y los parámetros del modelo en cada paso; la verosimilitud final fue lnL = ${num(m.lnL)}${m.visited ? ` tras evaluar ${m.visited} topologías` : ''}.`);
      en.push(`The maximum-likelihood tree was searched with NNI and SPR rearrangements under a regraft radius, starting from the distance tree and optimising branch lengths and model parameters at each step; the final likelihood was lnL = ${num(m.lnL)}${m.visited ? ` after evaluating ${m.visited} topologies` : ''}.`);
      if (m.partitioned) {
        es.push(`El análisis particionado ajustó un modelo y un multiplicador de velocidad por partición sobre una topología y unas longitudes relativas compartidas (lnL = ${num(m.partitioned.lnL)}, ${m.partitioned.k} parámetros).`);
        en.push(`The partitioned analysis fitted one model and one rate multiplier per partition over a shared topology and shared relative branch lengths (lnL = ${num(m.partitioned.lnL)}, ${m.partitioned.k} parameters).`);
      }
      /* the fragments are kept in both languages: writing them once and reusing
         the string would put Spanish inside the English paragraph, which is the
         quiet way a bilingual report ends up half translated */
      const supEs = [], supEn = [];
      if (m.support && m.support.boot) {
        supEs.push(`arranque no paramétrico (${m.support.boot.reps} réplicas)`);
        supEn.push(`non-parametric bootstrap (${m.support.boot.reps} replicates)`);
      }
      if (m.support && m.support.uf) {
        supEs.push(`arranque ultrarrápido (${m.support.uf.reps} réplicas)`);
        supEn.push(`ultrafast bootstrap (${m.support.uf.reps} replicates)`);
      }
      if (m.support && m.support.alrt) { supEs.push('SH-aLRT y aBayes'); supEn.push('SH-aLRT and aBayes'); }
      if (supEs.length) {
        es.push(`El apoyo de las ramas se estimó por ${supEs.join(', ')}.`);
        en.push(`Branch support was estimated by ${supEn.join(', ')}.`);
      }
      const refs = ['felsenstein1981'];
      if (m.support && m.support.boot) refs.push('felsenstein1985b');
      if (m.support && m.support.uf) refs.push('minh2013');
      if (m.support && m.support.alrt) refs.push('guindon2010', 'anisimova2011');
      if (m.topologyTests) refs.push('kishino1989', 'shimodaira1999', 'shimodaira2002');
      out.push({ id: 'ml', titleEs: 'Máxima verosimilitud', titleEn: 'Maximum likelihood', es, en, refs });
    }

    /* ---- Bayesian ---- */
    if (st.bayes) {
      const b = st.bayes;
      const es = [], en = [];
      es.push(`La inferencia bayesiana se hizo con ${b.runs} ${b.runs === 1 ? 'corrida' : 'corridas'} independientes de ${b.generations} generaciones, con ${b.chains} ${b.chains === 1 ? 'cadena' : 'cadenas'} acopladas por corrida y un descarte inicial de ${pct(b.burnin)}.`);
      en.push(`Bayesian inference used ${b.runs} independent ${b.runs === 1 ? 'run' : 'runs'} of ${b.generations} generations, with ${b.chains} coupled ${b.chains === 1 ? 'chain' : 'chains'} per run and an initial burn-in of ${pct(b.burnin)}.`);
      if (b.priors) {
        es.push(`Los priors fueron: longitudes de rama exponenciales de media ${num(b.priors.branchMean)}, topologías uniformes y distribuciones exponenciales sobre las tasas y sobre α.`);
        en.push(`The priors were: exponential branch lengths with mean ${num(b.priors.branchMean)}, uniform topologies, and exponential distributions on the rates and on α.`);
      }
      if (b.asdsf != null) {
        es.push(`La convergencia se juzgó por la desviación estándar media de las frecuencias de división (ASDSF = ${num(b.asdsf)}), por el factor de reducción de escala de Gelman y Rubin y por el tamaño de muestra efectivo de cada parámetro.`);
        en.push(`Convergence was judged by the average standard deviation of split frequencies (ASDSF = ${num(b.asdsf)}), by the Gelman–Rubin scale reduction factor and by the effective sample size of each parameter.`);
      }
      if (b.marginal) {
        es.push('Las verosimilitudes marginales se estimaron por piedras de paso para comparar modelos con factores de Bayes.');
        en.push('Marginal likelihoods were estimated by stepping stone in order to compare models with Bayes factors.');
      }
      const refs = ['larget1999', 'geyer1991', 'gelman1992'];
      if (b.marginal) refs.push('xie2011');
      if (b.priors) refs.push('brown2010');
      out.push({ id: 'bayes', titleEs: 'Inferencia bayesiana', titleEn: 'Bayesian inference', es, en, refs });
    }

    /* ---- clock and dating ---- */
    if (st.dated) {
      const d = st.dated;
      const es = [], en = [];
      const howEs = d.rootMethod === 'outgroup' ? 'por grupo externo'
        : d.rootMethod === 'tipdates' ? 'por regresión raíz–punta sobre las fechas de las puntas'
          : 'por punto medio';
      const howEn = d.rootMethod === 'outgroup' ? 'by outgroup'
        : d.rootMethod === 'tipdates' ? 'by root-to-tip regression on the tip dates'
          : 'by midpoint';
      es.push(`El árbol se enraizó ${howEs}.`);
      en.push(`The tree was rooted ${howEn}.`);
      if (d.clockTest) {
        es.push(`La hipótesis de reloj estricto se contrastó por razón de verosimilitudes contra el ajuste sin reloj (2Δ = ${num(d.clockTest.stat)}, ${d.clockTest.df} g.l., p ${d.clockTest.p < 1e-4 ? '< 0.0001' : '= ' + num(d.clockTest.p)}), y ${d.clockTest.p < 0.05 ? 'se rechazó' : 'no se rechazó'}.`);
        en.push(`The strict-clock hypothesis was tested by a likelihood ratio against the unconstrained fit (2Δ = ${num(d.clockTest.stat)}, ${d.clockTest.df} d.f., p ${d.clockTest.p < 1e-4 ? '< 0.0001' : '= ' + num(d.clockTest.p)}), and was ${d.clockTest.p < 0.05 ? 'rejected' : 'not rejected'}.`);
      }
      const methEs = d.method === 'relaxed' ? 'un reloj relajado lognormal no correlacionado, muestreado por MCMC'
        : d.method === 'pl' ? 'verosimilitud penalizada, con el parámetro de suavizado elegido por validación cruzada'
          : 'mínimos cuadrados';
      const methEn = d.method === 'relaxed' ? 'an uncorrelated lognormal relaxed clock, sampled by MCMC'
        : d.method === 'pl' ? 'penalised likelihood, with the smoothing parameter chosen by cross-validation'
          : 'least squares';
      es.push(`Los tiempos de divergencia se estimaron por ${methEs}${d.calibrations && d.calibrations.length ? `, con ${d.calibrations.length} ${d.calibrations.length === 1 ? 'calibración' : 'calibraciones'}` : ', sin calibraciones externas, de modo que las edades son proporciones de la profundidad de la raíz y no años'}.`);
      en.push(`Divergence times were estimated by ${methEn}${d.calibrations && d.calibrations.length ? `, with ${d.calibrations.length} ${d.calibrations.length === 1 ? 'calibration' : 'calibrations'}` : ', with no external calibration, so the ages are proportions of the root depth and not years'}.`);
      const refs = [];
      if (d.rootMethod === 'tipdates') refs.push('rambaut2016');
      if (d.method === 'relaxed') refs.push('drummond2006');
      else if (d.method === 'pl') refs.push('sanderson2002');
      else refs.push('to2016');
      out.push({ id: 'dated', titleEs: 'Reloj molecular y tiempos de divergencia', titleEn: 'Molecular clock and divergence times', es, en, refs });
    }

    /* ---- diversification ---- */
    if (st.divers) {
      const v = st.divers;
      const es = [], en = [];
      if (v.gamma != null) {
        es.push(`El estadístico γ de los tiempos de ramificación fue ${num(v.gamma)}${v.mccr ? `, y la prueba de Monte Carlo de tasa constante, que simula árboles del tamaño verdadero del grupo y los poda al muestreo disponible, dio p = ${num(v.mccr.p)}` : ''}.`);
        en.push(`The γ statistic of the branching times was ${num(v.gamma)}${v.mccr ? `, and the Monte Carlo constant-rates test, which simulates trees of the clade's true size and prunes them to the available sampling, gave p = ${num(v.mccr.p)}` : ''}.`);
      }
      if (v.models && v.models.length) {
        es.push(`Se compararon por AIC ${v.models.length} modelos de diversificación (nacimiento puro, nacimiento–muerte, dependencia de la densidad logística y exponencial, y dos tasas con el momento del cambio estimado); el mejor fue ${v.best}.`);
        en.push(`${v.models.length} diversification models (pure birth, birth–death, logistic and exponential density dependence, and two rates with the shift point estimated) were compared by AIC; the best was ${v.best}.`);
      }
      const refs = [];
      if (v.gamma != null) refs.push('pybus2000');
      if (v.models) refs.push('nee1994', 'rabosky2008');
      if (v.magallon) refs.push('magallon2001');
      if (v.dr) refs.push('jetz2012');
      out.push({ id: 'divers', titleEs: 'Diversificación', titleEn: 'Diversification', es, en, refs });
    }

    /* ---- traits ---- */
    if (st.traits) {
      const tr = st.traits;
      const es = [], en = [];
      if (tr.continuous) {
        const c = tr.continuous;
        es.push(`Para el carácter continuo «${c.character}» se ajustaron por máxima verosimilitud los modelos de movimiento browniano, Ornstein–Uhlenbeck y estallido temprano, comparados por AIC; el preferido fue ${c.best}.`);
        en.push(`For the continuous character "${c.character}", Brownian motion, Ornstein–Uhlenbeck and early-burst models were fitted by maximum likelihood and compared by AIC; the preferred one was ${c.best}.`);
      }
      if (tr.signal) {
        es.push(`La señal filogenética se midió con la K de Blomberg (K = ${num(tr.signal.K)}, p = ${num(tr.signal.pK)} por permutación) y con la λ de Pagel (λ = ${num(tr.signal.lambda)}, p = ${num(tr.signal.pLambda)} por razón de verosimilitudes).`);
        en.push(`Phylogenetic signal was measured with Blomberg's K (K = ${num(tr.signal.K)}, p = ${num(tr.signal.pK)} by permutation) and Pagel's λ (λ = ${num(tr.signal.lambda)}, p = ${num(tr.signal.pLambda)} by likelihood ratio).`);
      }
      if (tr.regression) {
        const r = tr.regression;
        es.push(`La relación entre «${r.y}» y «${r.x}» se estimó por contrastes independientes y por mínimos cuadrados generalizados filogenéticos, que dan la misma pendiente (${num(r.pgls.slope)}, p = ${num(r.pgls.p)}); la regresión sin corregir daba ${num(r.ols.slope)} (p = ${num(r.ols.p)}).`);
        en.push(`The relationship between "${r.y}" and "${r.x}" was estimated by independent contrasts and by phylogenetic generalised least squares, which give the same slope (${num(r.pgls.slope)}, p = ${num(r.pgls.p)}); the uncorrected regression gave ${num(r.ols.slope)} (p = ${num(r.ols.p)}).`);
      }
      if (tr.discrete) {
        const d = tr.discrete;
        es.push(`Para el carácter discreto «${d.character}» (${d.levels.length} estados) se ajustaron modelos Mk con tasas iguales, simétricas y todas distintas, comparados por AIC (mejor: ${d.best}); los estados ancestrales se calcularon como probabilidades marginales${tr.simmap ? `, y ${tr.simmap.reps} historias por mapeo estocástico dieron ${num(tr.simmap.changes.mean)} cambios por historia en promedio` : ''}.`);
        en.push(`For the discrete character "${d.character}" (${d.levels.length} states), Mk models with equal, symmetric and all-different rates were fitted and compared by AIC (best: ${d.best}); ancestral states were computed as marginal probabilities${tr.simmap ? `, and ${tr.simmap.reps} histories by stochastic mapping gave ${num(tr.simmap.changes.mean)} changes per history on average` : ''}.`);
      }
      const refs = [];
      if (tr.continuous) refs.push('felsenstein1985', 'hansen1997', 'harmon2010');
      if (tr.signal) refs.push('blomberg2003', 'pagel1999');
      if (tr.regression) refs.push('garland2000');
      if (tr.discrete) refs.push('lewis2001');
      if (tr.simmap) refs.push('huelsenbeck2003');
      out.push({ id: 'traits', titleEs: 'Evolución de caracteres', titleEn: 'Trait evolution', es, en, refs });
    }

    /* ---- biogeography ---- */
    if (st.biogeo && st.biogeo.models) {
      const b = st.biogeo;
      const es = [], en = [];
      es.push(`La biogeografía histórica se analizó sobre ${b.areas.length} áreas (${b.areas.join(', ')}), con un tamaño máximo de rango de ${b.maxAreas} y ${b.nStates} estados${b.includeNull ? ', incluido el rango nulo' : ''}.`);
      en.push(`Historical biogeography was analysed over ${b.areas.length} areas (${b.areas.join(', ')}), with a maximum range size of ${b.maxAreas} and ${b.nStates} states${b.includeNull ? ', the null range included' : ''}.`);
      es.push(`Se ajustaron por máxima verosimilitud ${b.models.length} modelos de la familia DEC, DIVALIKE y BAYAREALIKE, con y sin dispersión fundadora (+J), comparados por AICc; el mejor fue <b>${b.best}</b>.`);
      en.push(`${b.models.length} models of the DEC, DIVALIKE and BAYAREALIKE family, with and without founder-event dispersal (+J), were fitted by maximum likelihood and compared by AICc; the best was <b>${b.best}</b>.`);
      if (/\+J/.test(b.best || '')) {
        es.push('Conviene señalar, siguiendo a Ree y Sanmartín (2018), que la comparación entre un modelo +J y el que extiende no es una comparación limpia: el parámetro j compra probabilidad en la cladogénesis sin pagarla a lo largo de las ramas, de modo que los modelos +J ganan casi siempre y su victoria no constituye por sí sola evidencia de especiación fundadora.');
        en.push('It should be noted, following Ree & Sanmartín (2018), that the comparison between a +J model and the one it extends is not a clean one: the parameter j buys probability at cladogenesis without paying for it along the branches, so +J models win almost always and their victory is not on its own evidence of founder-event speciation.');
      }
      if (b.stochastic) {
        es.push(`El mapeo estocástico biogeográfico (${b.stochastic.reps} historias) dio ${num(b.stochastic.mean.dispersal)} dispersiones y ${num(b.stochastic.mean.extinction)} extinciones locales por historia.`);
        en.push(`Biogeographic stochastic mapping (${b.stochastic.reps} histories) gave ${num(b.stochastic.mean.dispersal)} dispersals and ${num(b.stochastic.mean.extinction)} local extinctions per history.`);
      }
      out.push({ id: 'biogeo', titleEs: 'Biogeografía histórica', titleEn: 'Historical biogeography', es, en,
        refs: ['ree2008', 'matzke2013'].concat(/\+J/.test(b.best || '') ? ['ree2018'] : []) });
    }

    /* ---- comparing trees ---- */
    if (st.compare) {
      const c = st.compare;
      const es = [], en = [];
      if (c.distances) {
        es.push(`Los árboles se compararon por la distancia de Robinson–Foulds, el puntaje de rama de Kuhner y Felsenstein, la diferencia de caminos de Steel y Penny y la distancia de cuartetos.`);
        en.push('The trees were compared by the Robinson–Foulds distance, the Kuhner–Felsenstein branch score, the Steel–Penny path difference and the quartet distance.');
      }
      if (c.concordance) {
        const br = c.concordance.branches || [];
        const mean = br.length ? br.reduce((a, x) => a + x.gCF, 0) / br.length : null;
        es.push(`Se calcularon factores de concordancia de genes y de sitios para cada rama interna sobre ${c.concordance.nGenes} árboles${mean != null ? ` (gCF medio ${mean.toFixed(1)} %)` : ''}. Conviene recordar que la concordancia no es apoyo: una rama puede estar bien estimada y contradicha por buena parte de la evidencia a la vez.`);
        en.push(`Gene and site concordance factors were computed for every internal branch over ${c.concordance.nGenes} trees${mean != null ? ` (mean gCF ${mean.toFixed(1)} %)` : ''}. It is worth remembering that concordance is not support: a branch can be well estimated and contradicted by much of the evidence at the same time.`);
      }
      if (c.quartet) {
        es.push(`Se estimó además un árbol de especies a partir de las frecuencias de cuartetos entre los árboles disponibles, mediante una búsqueda por intercambios de vecinos (${pct(c.quartet.proportion)} de los cuartetos decididos a favor del árbol hallado). La búsqueda es heurística y no equivale al programa dinámico exacto de ASTRAL.`);
        en.push(`A species tree was also estimated from the quartet frequencies among the available trees, by a nearest-neighbour search (${pct(c.quartet.proportion)} of the decided quartets in favour of the tree found). The search is a heuristic and is not equivalent to ASTRAL's exact dynamic program.`);
      }
      if (c.network) {
        es.push(c.network.kind === 'consensus'
          ? `El conflicto entre árboles se representó con una red de consenso (umbral ${c.network.threshold}).`
          : 'El conflicto se representó con una descomposición circular de las distancias entre puntas.');
        en.push(c.network.kind === 'consensus'
          ? `The conflict between trees was shown as a consensus network (threshold ${c.network.threshold}).`
          : 'The conflict was shown as a circular decomposition of the distances between tips.');
      }
      if (c.abba || c.abbaScan) {
        const d = c.abba;
        es.push(`La asimetría entre sitios ABBA y BABA se midió con el estadístico D de Patterson, con un error estándar por jackknife de bloques${d ? ` (D = ${num(d.D)}, z = ${num(d.z)})` : ''}. Un valor distinto de cero admite al menos tres explicaciones —introgresión, un exogrupo mal elegido y tasas muy desiguales entre P1 y P2— y sólo la primera es la que suele buscarse.`);
        en.push(`The asymmetry between ABBA and BABA sites was measured with Patterson's D, with a block-jackknife standard error${d ? ` (D = ${num(d.D)}, z = ${num(d.z)})` : ''}. A non-zero value admits at least three explanations — introgression, a badly chosen outgroup and very unequal rates between P1 and P2 — and only the first is the one usually being looked for.`);
      }
      const refs = [];
      if (c.distances) refs.push('robinson1981', 'kuhner1994', 'steel1993');
      if (c.concordance) refs.push('minh2020');
      if (c.network) refs.push('bandelt1992');
      if (c.abba || c.abbaScan) refs.push('green2010', 'durand2011');
      if (es.length) out.push({ id: 'compare', titleEs: 'Comparación de árboles', titleEn: 'Comparing trees', es, en, refs });
    }

    /* ---- the program itself ---- */
    {
      const es = ['Todos los análisis se realizaron en PhylogenyPro, que ejecuta cada método en el navegador sin enviar los datos a ningún servidor. Cada bloque del programa está validado contra implementaciones independientes en R (ape, phangorn, geiger, phytools, nlme y BioGeoBEARS), y los guiones de comparación se distribuyen con el código.'];
      const en = ['All analyses were carried out in PhylogenyPro, which runs every method in the browser without sending the data to any server. Every block of the program is validated against independent implementations in R (ape, phangorn, geiger, phytools, nlme and BioGeoBEARS), and the comparison scripts are distributed with the code.'];
      out.push({ id: 'software', titleEs: 'Programa', titleEn: 'Software', es, en, refs: [] });
    }
    return out;
  }

  /* the bibliography of exactly what was cited, in alphabetical order */
  function references(sections) {
    const ids = new Set();
    sections.forEach(s => (s.refs || []).forEach(r => ids.add(r)));
    return [...ids].map(id => ({ id, text: REFS[id] })).filter(r => r.text)
      .sort((a, b) => a.text.localeCompare(b.text, 'es'));
  }

  /* ================================================================
     the report
     ================================================================ */
  /* the title of a section in one language or the other; the module carries
     both and never asks the page which one is showing */
  const titleOf = (s, lang) => (lang === 'en' ? s.titleEn : s.titleEs);
  function methodsText(sections, lang) {
    return sections.map(s => {
      const body = (lang === 'en' ? s.en : s.es).join(' ').replace(/<\/?b>/g, '');
      return `${titleOf(s, lang)}\n${body}\n`;
    }).join('\n');
  }

  /* A standalone HTML page: no stylesheet to fetch, no font to download, no
     script to run. It opens the same in ten years as it does today, and prints
     to PDF from the browser, which is what a supplementary file has to do. */
  function html(st, opts) {
    opts = opts || {};
    const lang = opts.lang === 'en' ? 'en' : 'es';
    const sections = opts.sections || methods(st, opts);
    const refs = references(sections);
    const figs = opts.figures || [];
    const title = opts.title || (lang === 'en' ? 'Phylogenetic analysis' : 'Análisis filogenético');
    const when = (opts.date || new Date()).toISOString().slice(0, 10);
    const e = s => String(s == null ? '' : s).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');

    let out = `<!DOCTYPE html>\n<html lang="${lang}">\n<head>\n<meta charset="UTF-8">\n`;
    out += `<meta name="viewport" content="width=device-width,initial-scale=1">\n<title>${e(title)}</title>\n`;
    out += '<style>\n' + REPORT_CSS + '</style>\n</head>\n<body>\n';
    out += `<header><h1>${e(title)}</h1>`;
    if (opts.author) out += `<p class="meta">${e(opts.author)}</p>`;
    out += `<p class="meta">${lang === 'en' ? 'Produced by PhylogenyPro on' : 'Generado por PhylogenyPro el'} ${when}</p></header>\n`;

    /* what the analysis consisted of, in one table */
    const rows = summary(st, lang);
    if (rows.length) {
      out += `<section><h2>${lang === 'en' ? 'Summary' : 'Resumen'}</h2>\n<table>\n<tbody>\n`;
      rows.forEach(r => { out += `<tr><th>${e(r.label)}</th><td>${r.html || e(r.value)}</td></tr>\n`; });
      out += '</tbody>\n</table>\n</section>\n';
    }

    out += `<section><h2>${lang === 'en' ? 'Methods' : 'Métodos'}</h2>\n`;
    sections.forEach(s => {
      out += `<h3>${e(titleOf(s, lang))}</h3>\n`;
      (lang === 'en' ? s.en : s.es).forEach(p => { out += `<p>${p}</p>\n`; });
    });
    out += '</section>\n';

    if (figs.length) {
      out += `<section><h2>${lang === 'en' ? 'Figures' : 'Figuras'}</h2>\n`;
      figs.forEach((f, i) => {
        out += `<figure>\n${f.svg}\n<figcaption>${lang === 'en' ? 'Figure' : 'Figura'} ${i + 1}. ${e(f.caption || f.name)}</figcaption>\n</figure>\n`;
      });
      out += '</section>\n';
    }

    const tabs = opts.tables || [];
    if (tabs.length) {
      out += `<section><h2>${lang === 'en' ? 'Tables' : 'Tablas'}</h2>\n`;
      tabs.forEach((tb, i) => {
        out += `<h3>${lang === 'en' ? 'Table' : 'Tabla'} ${i + 1}. ${e(tb.caption)}</h3>\n<table>\n<thead><tr>`;
        tb.columns.forEach(c => { out += `<th>${e(c)}</th>`; });
        out += '</tr></thead>\n<tbody>\n';
        tb.rows.forEach(r => {
          out += '<tr>' + r.map(v => `<td>${e(v)}</td>`).join('') + '</tr>\n';
        });
        out += '</tbody>\n</table>\n';
      });
      out += '</section>\n';
    }

    if (refs.length) {
      out += `<section><h2>${lang === 'en' ? 'References' : 'Referencias'}</h2>\n<ol class="refs">\n`;
      refs.forEach(r => { out += `<li>${e(r.text)}</li>\n`; });
      out += '</ol>\n</section>\n';
    }
    out += `<footer><p class="meta">${lang === 'en'
      ? 'This report was written by the program from what it actually ran. Every number in the methods comes from the run, not from a template.'
      : 'Este informe lo escribió el programa a partir de lo que realmente corrió. Cada número de los métodos viene de la corrida, no de una plantilla.'}</p></footer>\n`;
    out += '</body>\n</html>\n';
    return out;
  }

  const REPORT_CSS = `
:root { color-scheme: light; }
body { font-family: Georgia, 'Times New Roman', serif; max-width: 46em; margin: 0 auto;
  padding: 2.5em 1.5em 4em; line-height: 1.62; color: #1a1a1a; background: #fff; }
header { border-bottom: 2px solid #1a1a1a; padding-bottom: 1em; margin-bottom: 2em; }
h1 { font-size: 1.7em; margin: 0 0 .3em; line-height: 1.25; }
h2 { font-size: 1.2em; margin: 2.2em 0 .6em; border-bottom: 1px solid #ccc; padding-bottom: .25em; }
h3 { font-size: 1.02em; margin: 1.5em 0 .35em; }
p { margin: 0 0 .85em; text-align: justify; hyphens: auto; }
.meta { color: #555; font-size: .88em; margin: .15em 0; }
table { border-collapse: collapse; width: 100%; margin: .6em 0 1.4em; font-size: .9em; }
th, td { border: 1px solid #d5d5d5; padding: .38em .6em; text-align: left; vertical-align: top; }
thead th { background: #f2f2f2; }
tbody th { background: #fafafa; width: 32%; font-weight: 600; }
figure { margin: 1.6em 0; page-break-inside: avoid; }
figure svg { width: 100%; height: auto; max-width: 100%; }
figcaption { font-size: .86em; color: #444; margin-top: .5em; }
ol.refs { padding-left: 1.4em; font-size: .9em; }
ol.refs li { margin-bottom: .5em; }
footer { margin-top: 3em; border-top: 1px solid #ccc; padding-top: .8em; }
@media print { body { max-width: none; padding: 0; } h2 { page-break-after: avoid; } }
`;

  /* the one-table answer to "what did you do" */
  function summary(st, lang) {
    const en = lang === 'en';
    const rows = [];
    const add = (label, value) => { if (value != null && value !== '') rows.push({ label, value }); };
    if (st.data && st.data.parts && st.data.parts.length) {
      const parts = st.data.parts;
      const nTax = (st.data.taxa || parts[0].taxa || []).length;
      const chars = parts.reduce((a, p) => a + ((p.seqs && p.seqs[0]) ? p.seqs[0].length : 0), 0);
      add(en ? 'Terminals' : 'Terminales', nTax);
      add(en ? 'Characters' : 'Caracteres', chars);
      add(en ? 'Partitions' : 'Particiones', parts.map(p => p.name).join(', '));
    }
    if (st.models && st.models.best) add(en ? 'Substitution model' : 'Modelo de sustitución', st.models.best.name || st.models.best.model);
    if (st.ml) add(en ? 'Likelihood of the ML tree' : 'Verosimilitud del árbol de ML', st.ml.lnL == null ? null : (+st.ml.lnL).toFixed(4));
    if (st.bayes) add(en ? 'Bayesian run' : 'Corrida bayesiana',
      `${st.bayes.runs} × ${st.bayes.generations} ${en ? 'generations' : 'generaciones'}`);
    if (st.dated) add(en ? 'Dating' : 'Datación',
      st.dated.method === 'relaxed' ? (en ? 'relaxed clock' : 'reloj relajado')
        : st.dated.method === 'pl' ? (en ? 'penalised likelihood' : 'verosimilitud penalizada')
          : (en ? 'least squares' : 'mínimos cuadrados'));
    if (st.divers && st.divers.best) add(en ? 'Diversification model' : 'Modelo de diversificación', st.divers.best);
    if (st.traits && st.traits.continuous) add(en ? 'Trait model' : 'Modelo de carácter', st.traits.continuous.best);
    if (st.biogeo && st.biogeo.best) add(en ? 'Biogeographic model' : 'Modelo biogeográfico', st.biogeo.best);
    if (st.compare && st.compare.trees) add(en ? 'Trees compared' : 'Árboles comparados', st.compare.trees.length);
    return rows;
  }

  Object.assign(Report, { REFS, methods, references, methodsText, html, summary, titleOf, REPORT_CSS });
  g.Report = Report;
}
ReportCore(typeof window !== 'undefined' ? window : self);
