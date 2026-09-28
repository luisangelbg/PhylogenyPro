# PhylogenyPro

**Molecular phylogenetics in the browser, from the alignment to the dated tree.** All twelve blocks are complete.

A web platform (HTML + JavaScript, no installation; it also runs offline from a local copy) that takes a set of
sequences — chloroplast, mitochondrial, nuclear ribosomal or low-copy nuclear DNA, amino acids, codons or
morphological characters — through the whole of a real phylogenetic study: alignment, substitution-model selection,
inference by parsimony, maximum likelihood and Bayesian MCMC, branch support and topology tests, molecular clock and
fossil calibration, diversification, ancestral traits and ancestral areas, ending in publication-ready figures and a
reproducible report.

It fills the gap left by [PopGeneticsPro](https://github.com/luisangelbg/PopGeneticsPro), which covers population
genetics (within species) but does not search tree space, select models or reconstruct history above the species
level.

## How to open it

1. Right-click **`server.ps1`** → *Run with PowerShell* (or double-click `Open PhylogenyPro.bat`).
   The browser opens at `http://localhost:9800`.
   If the port is busy: `powershell -ExecutionPolicy Bypass -File server.ps1 -Port 9801`
2. Double-clicking `index.html` also works — nothing is fetched from disk at run time.
3. To use it from a tablet on the same Wi-Fi network, run the script *as administrator*; it prints the address.

## Status of the blocks

| Block | Content | Status |
|---|---|---|
| 1 | Home: theory, method gallery, the honest table of sizes it can handle, and two live labs — the Felsenstein zone and substitution saturation | ✅ done |
| 2 | Data and alignment: FASTA, PHYLIP, NEXUS, Clustal, MEGA and GenBank readers; progressive aligner (also codon-aware); alignment viewer and manual editing; trimming of ambiguous blocks; concatenation and partitions, each declaring its genome | ✅ done |
| 3 | Substitution models: base composition, saturation, model selection (DNA, protein, morphology) by AIC/AICc/BIC, hierarchical LRT, partition scheme by greedy merging | ✅ done |
| 4 | Parsimony and distances: corrected distances, NJ, BIONJ, UPGMA, balanced minimum evolution, Fitch and Sankoff parsimony with NNI/SPR/TBR, CI/RI/RC/HI, Bremer, bootstrap and jackknife, consensus trees, tree drawing | ✅ done |
| 5 | Maximum likelihood: NNI and SPR search with a regraft radius, partitioned models, bootstrap in Web Workers, ultrafast bootstrap, SH-aLRT, aBayes, KH/SH/AU topology tests, constrained trees | ✅ done |
| 6 | Bayesian inference: Metropolis-coupled MCMC in Web Workers, priors written down, ASDSF, PSRF, ESS, HPD, traces, majority-rule consensus with posterior probabilities, marginal likelihoods by stepping stone, Tracer-compatible log | ✅ done |
| 7 | Molecular clock and dating: rooting by outgroup, midpoint or tip dates, the strict-clock test, fossil calibrations, least squares, penalised likelihood with cross-validation, relaxed-clock MCMC, and a chronogram on the geological scale | ✅ done |
| 8 | Diversification: lineages through time, γ with the Monte Carlo test for incomplete sampling, five models compared by AIC, Magallón–Sanderson from age and richness, and the DR statistic | ✅ done |
| 9 | Trait evolution: a character table of its own, Mk with ER/SYM/ARD, marginal ancestral states drawn as pies, stochastic mapping, BM/OU/EB, independent contrasts, PGLS, Blomberg's K and Pagel's λ | ✅ done |
| 10 | Historical biogeography: DEC, DEC+J, DIVALIKE, BAYAREALIKE compared by AICc, dispersal multipliers, time strata, ancestral ranges with their cladogenetic corners, biogeographic stochastic mapping | ✅ done |
| 11 | Comparing trees: five distances including quartets, tanglegrams that untangle, gene and site concordance factors, a quartet species tree, consensus and circular split networks, ABBA-BABA with a block jackknife | ✅ done |
| 12 | Figures and report: tree studio, Newick/NEXUS/phyloXML export, figures at any resolution with the pixel count shown first, a methods section written from the run, and a reproducible ZIP the app writes itself | ✅ done |

## What Block 1 already computes

The home page is not a brochure. Its two labs run real phylogenetics in the page:

- **The Felsenstein zone.** For a four-taxon tree, the 256 site patterns collapse by symmetry into the 15 set
  partitions of {A, B, C, D}, so the probability of every pattern under Jukes–Cantor is exact. The lab shows the
  infinite-data answer (which pattern wins, and therefore whether parsimony is consistent), finds the boundary of
  the zone by bisection, and simulates finite alignments that are analysed by parsimony, by neighbour joining over
  Jukes–Cantor distances and by maximum likelihood with the branch lengths optimised on each of the three topologies.
- **Saturation and distance correction.** Two sequences diverging under K80 with gamma-distributed rates, with the
  raw proportion of differences, the Jukes–Cantor and Kimura corrections and the transition/transversion curves,
  including the point where the correction stops having a solution.

## What Block 2 already does

- **Reads** FASTA, PHYLIP (strict and relaxed, sequential and interleaved), NEXUS (DATA, CHARACTERS, SETS with
  charsets, TREES, MATCHCHAR), Clustal, MEGA, GenBank flat files and spreadsheets for morphological matrices,
  guessing the format and the data type and reporting whatever it had to guess.
- **Aligns** with its own progressive aligner: k-mer distances → UPGMA guide tree → profile–profile alignment with
  affine gaps (Gotoh), optional iterative refinement, and **codon-aware alignment** that no gap can break.
- **Shows** the alignment on a canvas that handles thousands of columns, coloured by residue or by difference from
  the consensus, marking variable and parsimony-informative columns and the codon positions.
- **Lets you edit it by hand**: insert or remove gaps, delete columns, move and drop sequences, with unlimited undo.
- **Trims** ambiguous blocks with Gblocks (Castresana 2000) or trimAl-style thresholds, showing what would be
  removed before removing it, and optionally keeping whole codons.
- **Checks the data**: repeated names, identical sequences, sequences of very different length, unexpected
  characters, gap-heavy taxa and columns, internal stop codons, unequal G+C, and taxa missing from a partition.
- **Builds the partition map**, with the genome of each partition (chloroplast, mitochondrion, nuclear ribosomal,
  low-copy nuclear) carried into every later block, concatenation with missing data handled, splitting by codon
  position, and import and export of partition schemes.

## What Block 3 already does

- **Checks the assumptions before choosing anything**: base composition per sequence with a χ² test of homogeneity,
  and saturation, with transitions and transversions plotted against the corrected distance and the proportion of
  changes the correction has to recover.
- **Fits every candidate model** on one fixed tree — 56 DNA models (seven rate schemes × equal or empirical
  frequencies × nothing, +I, +G, +I+G), twelve empirical protein matrices including cpREV for chloroplast genes,
  and Mk and Mkv for morphology — and ranks them by AIC, AICc and BIC with Akaike weights.
- **Refits the winner properly**, with every branch optimised, and reports its rate matrix, its equilibrium
  frequencies, the shape of the gamma and the proportion of invariable sites, each with what it means.
- **Runs the hierarchical likelihood ratio tests** as well, for reviewers who ask for them, saying plainly why the
  information criteria are the default.
- **Finds a partition scheme** by greedy merging in the style of PartitionFinder, so that partitions that evolve
  alike can share a model instead of spending parameters on each.
- One model per turn of the event loop: the page never freezes, the wait is estimated out loud after two models,
  and the run can be cancelled.

## What Block 4 already does

- **Distance trees**: NJ, BIONJ, UPGMA, WPGMA and balanced minimum evolution (Pauplin's criterion with NNI and SPR
  swapping and Desper–Gascuel branch lengths), over any of the corrected distances, with or without a gamma
  correction and with pairwise or complete deletion.
- **Parsimony**: Fitch for unordered characters and Sankoff for a general cost matrix (ordered characters,
  transversion weighting), a heuristic search by stepwise addition from several random orders followed by NNI, SPR
  or TBR swapping, and every equally short tree kept.
- **The indices** — consistency, retention, rescaled consistency and homoplasy — with what they mean in words:
  a low CI is normal in molecules and falls as taxa are added; what matters is the RI.
- **Support**: bootstrap and jackknife with a live estimate of the time left, and Bremer support (how many extra
  steps it costs to lose each clade).
- **Consensus** trees: strict, majority-rule, semistrict and greedy, with the split frequencies.
- **A comparison table** of every tree the block has produced (and the one from Block 3) by Robinson–Foulds
  distance, with a warning when two criteria disagree badly — the usual cause being long-branch attraction.
- **Tree drawing** (`js/treeview.js`), shared with the blocks that follow: rectangular, circular and unrooted
  layouts, cladograms, support as numbers or coloured dots, rooting by midpoint or outgroup, and SVG export.

## What Block 5 already does

- **Searches tree space** by NNI and then SPR, scoring each candidate the way PhyML and RAxML do: only the branches
  the rearrangement touched are reoptimised, over an incremental likelihood that recomputes just the path to the
  root. SPR is capped by a regraft radius and runs lazily first — every candidate scored with the lengths it
  inherits, only the best dozen refined — which is the difference between a search that finishes and a frozen page.
- **Measures support four ways at once**, each with its own threshold printed in the table header because they do
  not mean the same thing: the classical bootstrap (≥ 70), the ultrafast bootstrap by RELL resampling of the trees
  the search visited (≥ 95), SH-aLRT (≥ 80) and aBayes (≥ 0.95).
- **Runs the bootstrap on every core but one**, in Web Workers built from the engine's own source text, so it works
  the same when the app is opened by double-click; if workers are blocked the same jobs run on the page in chunks.
  The model is held at the values fitted to the real data, as RAxML and IQ-TREE do.
- **Tests hypotheses properly**: force a group to be monophyletic, find the best tree that obeys the constraint, and
  let KH, SH and AU say whether the data actually reject it — rather than concluding anything from its absence.
- **Fits a model per partition** over one shared topology and one set of relative branch lengths, with a rate
  multiplier each, normalised to a site-weighted mean of one, and says whether the BIC pays for the extra
  parameters.

## What Block 6 already does

- **Samples the posterior** with Metropolis–Hastings over topology, branch lengths and model parameters, coupled
  across several heated chains (MC³) so the sampler can leave one island of trees for another, and across several
  independent runs, because without them nothing can be said about convergence.
- **Makes the assumptions explicit.** The priors are written down and the panel says what the chosen numbers mean:
  with thirty-three branches and a mean of 0.1, the tree measures 3.3 substitutions per site before the data are
  seen. Other methods have assumptions too; this one has them in a box you can change.
- **Shows the diagnostics before the tree**: ASDSF between runs, PSRF, effective sample size, HPD intervals, the
  trace of the log-likelihood run by run, and the posterior of any parameter as a histogram with its interval
  marked. If any of them says the chains have not converged, the block says so in words instead of drawing a
  confident tree.
- **Summarises the posterior** as a majority-rule consensus with posterior probabilities and mean branch lengths,
  next to the most sampled topology, the size of the 95 % credible set and — when Block 5 has run — the bootstrap
  of the same clades side by side, because a posterior probability and a bootstrap proportion are not the same
  quantity and are not read against the same threshold.
- **Compares models by Bayes factors**, with marginal likelihoods estimated by stepping stone, on Kass and
  Raftery's scale.

## What Block 7 already does

- **Asks for things in the order that matters**: a root, then whether the clock holds, then the calibrations, then
  the dates. Any other order hides a decision.
- **Roots by tip dates** with a root-to-tip regression that searches every branch and the position along it, and
  says from the sign of the slope whether the dates given are calendar years or ages before the present.
- **Tests the strict clock** against the ultrametric fit and explains that rejection calls for a relaxed clock, not
  for giving up on dating.
- **Takes calibrations in the shapes a fossil justifies**: a fossil is a minimum age and almost never a maximum, so
  the offset exponential and offset lognormal are there next to fixed ages and uniform intervals.
- **Dates three ways**: least squares, penalised likelihood with the smoothing chosen by cross-validation, and a
  relaxed-clock MCMC that gives a 95 % interval for every node — the only one of the three that propagates the
  uncertainty of the calibrations.
- **Draws a chronogram** with time running backwards, the geological periods and epochs behind the tree, and the
  interval of every node as a bar.
- **Refuses to invent a calibration**: without one it dates the tree with the root at 1 and says the numbers are
  proportions, not years.

## What Block 8 already does

- **Draws the lineages through time** on a logarithmic axis with the constant-rate line behind them, because the
  curve alone says nothing.
- **Computes γ and says what it can and cannot mean**: a negative value is a slowdown or a sampling gap, and the
  number by itself cannot tell them apart.
- **Settles that with the Monte Carlo test**, simulating trees of the clade's real size and pruning them down to the
  species in the tree. On the example, eighteen species of a hundred and twenty give a mean γ of −2.2 with nothing
  having slowed down at all — which is exactly why the raw statistic is not enough.
- **Fits five models to the same branching times** and compares them by AIC with Akaike weights, warning that a
  density-dependent model winning does not prove diversity has a ceiling.
- **Estimates a rate from an age and a richness alone** when that is all there is, at three assumed extinction
  fractions, because the assumption moves the answer by half.
- **Gives every species its own rate** with the DR statistic.

## What Block 9 already does

- **Reads a character table of its own** — commas, semicolons or tabs, the type guessed per column, unknowns as
  `?`, `NA`, `-` or empty — and lines it up with the tree, forgiving the difference between `Bursera_palmeri` and
  `Bursera palmeri` and listing by name whatever still does not match. Each analysis prunes the tree to the tips
  that have the character it needs, and says how many that left.
- **Fits Brownian motion, Ornstein–Uhlenbeck and the early burst** to a continuous character, compares them by AIC
  with Akaike weights, and refuses to dress up a two-unit gap as a result. It also **checks that the tips are
  contemporaries**, because OU and the early burst assume it and a phylogram in substitutions is not.
- **Reconstructs ancestral values** and draws the traitgram, the one figure that puts the character and the tree in
  the same plane.
- **Measures phylogenetic signal both ways**: Blomberg's K with a permutation test and Pagel's λ with a likelihood
  ratio, and explains why K < 1 together with λ ≈ 1 is a description of a character with an optimum, not a
  contradiction.
- **Regresses two characters three ways at once** — ordinary, contrasts and PGLS — and checks in the interface that
  the last two give the same slope, which is Garland & Ives's (2000) identity and the cheapest test there is that
  neither is miscomputed.
- **Fits the Mk model** with ER, SYM and ARD over a general rate matrix, and draws the **marginal** ancestral states
  as pies, so a node reconstructed at 55 % looks doubtful instead of looking decided. The downward pass alone is
  available for anyone who wants to see where the reconstruction comes from.
- **Maps characters stochastically**: five hundred whole histories in under a second, turning node probabilities
  into how much of the tree each state held and how many changes that took.

## What Block 10 already does

- **Reads a range table** in the LAGRANGE format BioGeoBEARS also uses, or as a plain comma-separated table, lines
  it up with the tree by name, and says by name whatever does not match — including the species that end up in no
  area at all, which in this model means "extinct" and no tip can be.
- **Builds the state space out loud**: subsets of areas, the null range as an absorbing state, and a cap on the
  range size whose cost it states before the button is pressed, because the work grows with the cube of the number
  of states.
- **Fits the four models and their +J versions** by maximum likelihood, compares them by AICc with Akaike weights,
  and adds a likelihood ratio of each +J against the model it contains.
- **Says what a +J victory is worth.** It is the expected outcome, not evidence: Ree & Sanmartín (2018) showed that
  j buys probability at cladogenesis, where it costs nothing along the branches, while d and e pay for theirs over
  the whole time of the tree. The block reports how much it wins by and over which model, and leaves the claim
  there.
- **Reconstructs ancestral ranges** as marginal probabilities, and the **corners** with them: what each daughter
  started with immediately after the split. Both are drawn on the tree as pies, with a grey slice for everything
  below the threshold, so a node at 45 % looks like 45 %.
- **Maps histories stochastically**, turning those probabilities into counts of dispersals, local extinctions,
  sympatry, subset sympatry, vicariance and founder jumps, with their spread.
- **Lets geography change**: dispersal multipliers between areas, edited in the panel, and time strata in which
  both the multipliers and the set of areas that exist at all can differ.

## What Block 11 already does

- **Gathers every tree the app has made** and lets you paste more in Newick, checking the names and saying by name
  whatever does not match.
- **Measures the distance between each pair five ways** — Robinson–Foulds raw and normalised, weighted, the branch
  score, the path difference and the quartet distance — and explains beside them why RF and quartets can disagree:
  one misplaced tip sends RF to its maximum and barely moves the quartet count, and when that happens the
  disagreement is a tip, not a history.
- **Draws tanglegrams and untangles them first**, rotating branches until the lines cross as little as possible —
  and then says out loud that a tanglegram with no crossings does **not** mean the trees are the same.
- **Puts a concordance factor on every branch**: gCF from the gene trees, sCF from the alignment, against the third
  that pure chance gives, and keeps repeating that concordance is not support.
- **Estimates a species tree from quartet frequencies**, reporting the score of each candidate and naming its search
  as a heuristic rather than ASTRAL's exact program.
- **Builds split networks two ways** — a consensus network with the incompatible splits found and coloured, and the
  circular decomposition of a distance matrix, which returns a tree exactly when it is handed a tree's distances.
- **Runs the ABBA-BABA test**, for one trio or for every trio against an outgroup, with a block jackknife and the
  warning that D ≠ 0 has three explanations and only one is introgression.

## What Block 12 already does

- **A tree studio** in which every decision a tree figure contains is visible and reversible: layout, phylogram or
  cladogram, which support is shown and above what threshold, where the root goes, how the branches ladder, the
  fonts and widths, coloured clades — redrawn at every change.
- **Writes the tree in the three formats that matter**, each for what it can carry: Newick, NEXUS with a TRANSLATE
  block and the block FigTree reads, and phyloXML, the only one that holds a node's support, its age and its
  interval at once.
- **Shows the pixel count before making a PNG.** At 900 dpi a half-page figure is several thousand pixels a side;
  the block prints the number, and when it goes beyond what a browser canvas can hold it says so instead of
  returning a blank image. The SVG has no such limit.
- **Writes the methods section from the run**: no block that was not opened, no support that was not computed, no
  replicate count from a template, and a bibliography holding only what was cited. One self-contained HTML file
  that prints to PDF from the browser.
- **Packs everything into a ZIP it writes itself** — data, trees, figures, tables, report and an index explaining
  each file — with no library and no network, so it works from a double-clicked page.

## Validation

Everything numerical is checked against independent implementations in R 4.4.2 (ape 5.8.1, phangorn 2.12.1):

- `validation/block7/validate_block7.R`, with `clock_bench.html`: rooting by tip dates **reproducing `ape::rtt`**
  (slope −0.00158106 against −0.00158106, a difference of 3·10⁻⁹; r² 0.783919 against 0.783919); penalised
  likelihood against `ape::chronos` on the same tree, same λ and same calibration (15 clades, largest difference in
  age 0.0382 of a root fixed at 1, correlation r² 0.9906); the strict-clock test against
  `phangorn::optim.pml(optRooted = TRUE)`, where the app's ultrametric fit is **better by 11.8 log-likelihood
  units** and both reject the clock; and the exact identities of least-squares dating. **17 agreements, no
  mismatches.**

- `validation/block8/validate_block8.R`, with `div_bench.html`: γ identical to `ape::gammaStat` (0.44734655); pure
  birth identical to `ape::yule` in λ, its standard error and its log-likelihood; birth–death against
  `ape::birthdeath` (a 0.24684678 against 0.24684645, r 0.20638381 against 0.20638489, the same log-likelihood);
  Magallón–Sanderson identical to `geiger::bd.ms` at three extinction fractions; branching times against
  `ape::branching.times` to 1.2·10⁻⁹; and the null distribution of γ from the app's own simulator with mean −0.066
  and standard deviation 0.9999, which is the standard normal the Monte Carlo test rests on. **16 agreements, no
  mismatches.**

- `validation/block9/validate_block9.R`, with `traits_bench.html`: Brownian motion, OU and the early burst against
  `geiger::fitContinuous` (σ² 0.06350382, z₀ 0.25281475, lnL −26.628333; α 0.06985066 against 0.06985065); the 39
  contrasts of two characters against `ape::pic` to 1.7·10⁻¹⁴; PGLS against `nlme::gls` with `corBrownian` in every
  coefficient, standard error, t, p and likelihood, with the slope agreeing with the contrasts regression to
  4·10⁻¹⁷; K 0.90386290 and λ 0.97380521 against `phytools::phylosig`; the Mk rate 0.02920872 and lnL −11.307863 in
  `ape::ace`'s convention, which is `phytools::fitMk`'s plus log(k) — the whole difference is the root prior, and
  both are reported; marginal ancestral states matching `ape::ace` over all 39 nodes to 5.7·10⁻⁸ and
  `phytools::ancr` to 3.5·10⁻¹⁰; and stochastic mapping against `phytools::make.simmap` in branch time per state
  and in changes per history, from independent draws. **42 agreements, no mismatches.**

- `validation/block10/validate_block10.R`, with `biogeo_bench.html`: the six biogeographic models against
  **BioGeoBEARS 1.1.3** on the *Psychotria* of Hawaii — DEC −34.54195625 against −34.54195758, DEC+J −20.94758871
  against −20.94758855, DIVALIKE −33.14967582 against −33.14967538, DIVALIKE+J −21.08620942 against −21.08620912,
  BAYAREALIKE −40.33442556 against −40.33442795 and BAYAREALIKE+J −21.55262993 against −21.55262941, with d and e
  agreeing to 1e-5; the eighteen ancestral range reconstructions agreeing to **1.75·10⁻⁸** at the same parameters,
  the root in KO at 0.51655966 against 0.51655968; and the cladogenesis weight tables of the three models checked
  event by event against `cladoRcpp`. **41 agreements, no mismatches.**

- `validation/block11/validate_block11.R`, with `cmp_bench.html`: the five tree distances over 36 pairs each — 216
  numbers — against `phangorn::RF.dist`, `wRF.dist`, `KF.dist` and `path.dist`, RF exact and the rest to 10⁻¹⁰;
  the quartet distance against a brute-force count over all 1820 quartets of sixteen taxa (1767/53, 1691/129 and
  1820/0, exact); gCF against the same definition recomputed in R over all thirteen internal branches, to zero
  difference; the consensus network split for split against `phangorn::consensusNet`; and Patterson's D on a
  60 000-site alignment with 15 % introgression — ABBA 1089, BABA 753, D 0.18241042, jackknife standard error
  0.07884390, z 2.31356419, every digit. The circular decomposition has no package to check against, so it is
  checked against an identity: fed a tree's own distances it returns that tree's 13 internal and 16 terminal
  splits with its branch lengths to 3·10⁻¹⁶, reconstructing the distances to 9·10⁻¹⁶ and clipping nothing.
  **25 agreements, no mismatches.**

- `validation/block12/validate_block12.R`, with `export_bench.html`: this block hands the work over, so it is
  validated by making R do what a reader would do. R's own `unzip()` lists and extracts the package the app wrote,
  verifying every CRC32 on the way; `ape::read.tree` and `ape::read.nexus` both read the trees, with a topological
  distance of 0 between the two files and branch lengths identical to 10⁻⁶; `xml2` parses the phyloXML, whose tip
  names are the tree's and whose branch lengths sum to what the Newick sums to (2.446035 against 2.446036);
  `ape::read.dna` reads the alignment at 18 × 1711; the report is one file with no `<link>`, no `<script>` and no
  `<img>`, and — built from a distance tree alone — cites no clock, no biogeography and no MCMC, because none of
  them ran; and the figures are well-formed XML with no `var(--…)` left unresolved. **38 agreements, no
  mismatches.**

### The totals

Across the twelve folders: **910 unit tests**, and **258 numerical agreements with no mismatch** in the eleven
blocks that report a count, against ape 5.8.1, phangorn 2.12.1, geiger 2.0.11, phytools 2.4.4, nlme 3.1.166,
BioGeoBEARS 1.1.3 and xml2 1.3.8. Block 6 is not counted that way because it is checked in three different ways
instead: with the likelihood switched off its sampler reproduces its own prior exactly, its HPD and PSRF are
identical to `coda`'s, and an audit compares the likelihood the chain carries with a calculation from scratch at
every move, finding a mismatch of zero.

Every block also carries a browser bench under `validation/blockN/` that re-runs the comparison against the same
reference files, so the claim can be checked without R by opening one page.

- `validation/block1/validate_block1.R`: pattern probabilities recomputed from the eigen-decomposition of Q,
  log-likelihoods of the three topologies against `phangorn::optim.pml`, the boundary of the Felsenstein zone
  against `uniroot`, distance corrections against `ape::dist.dna`, gamma rate categories against `qgamma`.
  **17 agreements, no mismatches.**
- `validation/block2/validate_block2.R`: every format the app writes is read back identical by ape; the trimming
  rules and the column statistics are re-derived from scratch in R; the genetic code is checked against `ape::trans`.
  **11 agreements, no mismatches.**
- `validation/block2/align_bench.html`: the aligner scored against the true alignment that generated the simulated
  data — **99.98 % of the true residue pairs on rbcL, 98.14 % on the indel-rich trnL-F spacer, 98.93 % on ITS,
  99.60 % on the protein family, and 100 % when rbcL is aligned by codon.**
- `validation/block3/validate_block3.R`, with `like_bench.html` and `models_bench.html`: pairwise distances against
  `ape::dist.dna`; NJ, BIONJ and UPGMA identical in topology (RF = 0) and in total length; **the likelihood itself
  against `phangorn::pml` to 1e-7 for JC, F81, K80, HKY, GTR and their +I, +G and +I+G variants**; the discrete
  gamma against `phangorn::discrete.gamma`; the optimiser reaching the same optimum as `optim.pml`; and the model
  ranking agreeing with `phangorn::modelTest`. **29 agreements, no mismatches.**

- `validation/block5/validate_block5.R`, with `ml_bench.html` and `pool_bench.html`: the log-likelihood of given trees
  against `phangorn::optim.pml` (the app reaches **−7165.395908** where phangorn reaches −7165.396045); the
  **topology search reaching −7108.800752**, better than phangorn's NNI search (−7126.989550) and a shade better
  than its stochastic search (−7108.802104); the SH test agreeing with `phangorn::SH.test` (NJ tree p = 0.029
  against 0.026, ML tree 0.666 against 0.651); the bootstrap against `bootstrap.pml`; and the partitioned model
  against `phangorn::pmlPart`. **9 agreements, no mismatches.**

- `validation/block6/validate_block6.R`, with `prior_bench.html` and `post_bench.html`: **the sampler run with the
  likelihood switched off reproduces its prior** — topologies 0.3413 / 0.3311 / 0.3277 against 1/3, branch lengths
  with mean 0.10039 and sd 0.10038 against 0.1, tree length 0.50195 ± 0.22252 against Gamma(5, 0.1), and all 105
  unrooted topologies of six taxa visited at 0.00952 against 1/105; the HPD interval and the PSRF matching
  `coda::HPDinterval` and `coda::gelman.diag` on a series R and the browser build identically; and the 95 % HPD of
  κ and α containing the maximum-likelihood estimates from `phangorn::optim.pml`.

- `validation/block4/validate_block4.R`, with `pars_bench.html`: Fitch parsimony against `phangorn::parsimony` on a
  hand-checkable matrix and on rbcL; the search reaching **1916 steps, the same optimum as `pratchet`** and thirteen
  steps shorter than `optim.parsimony` with SPR; CI and RI against `phangorn::CI` and `RI`; the balanced
  minimum-evolution length reproducing `ape::fastme.bal` exactly on ape's own tree; consensus trees against
  `ape::consensus`. **13 agreements, no mismatches.**

The twelve empirical amino-acid matrices in `js/aamodels.js` are the ones published in their own articles, each
cited beside it in the file. `data/verificar_matrices_aa.R` checks every value in it against five independent
canonical distributions — PAML, IQ-TREE, MrBayes, RAxML, and the WAG file the authors themselves put on the EBI —
**5806 comparisons, not one of the 2280 exchangeabilities disagreeing**, and prints, one by one, the eight
equilibrium frequencies that do. BLOSUM62 comes from the matrix the NCBI distributes, which is in the public
domain. [`PROCEDENCIA.md`](PROCEDENCIA.md) sets out where every non-original element of the program comes from.

The unit tests live in `tests/index.html`; open that page and every line must be green (**910 tests**).

The measurements that set the size limits stated on the home page are in `research/02_viabilidad_computacional.md`,
with the benchmarks that produced them in `research/bench/`.

## What runs where

Everything is computed in the browser with plain JavaScript: no server, no upload, no installation. Long analyses
run in Web Workers built by the page itself, so the interface never freezes, every one of them estimates its own
running time before starting, and all of them can be cancelled. Everything stochastic draws from a seeded generator,
so a run can be reproduced from the seed reported with the results.

How fast, honestly: on the *Bursera* rbcL example (18 taxa × 621 sites) the maximum-likelihood search with NNI and
SPR takes about eight seconds, a hundred bootstrap replicates a few minutes across sixteen threads, and a hundred
thousand generations of coupled Bayesian chains a few minutes per run. This is
JavaScript rather than C, and the page says what it is about to cost before it starts. A worker that is handed a
single job pays its whole warm-up — a fresh JavaScript context starts interpreted — so the pool never gives out
fewer than three replicates per thread.

## Licence

Copyright © 2026 Luis Ángel Barrera-Guzmán.

PhylogenyPro is free software, released under the **GNU General Public License, version 3 or later**. The full text
is in [`LICENSE`](LICENSE), verbatim as the Free Software Foundation publishes it; the notice that applies it to
this program is at the top of `index.html`. You may use, study, share and modify it, provided that anything you
distribute built from it carries the same freedoms and the same licence.

The licence is a choice, not an inheritance: the program contains no third-party code, so nothing obliged it to be
copyleft. [`PROCEDENCIA.md`](PROCEDENCIA.md) sets out where every element of it that is not original comes from —
twelve published substitution matrices, the NCBI's BLOSUM62 and genetic code, and the ICS time scale — with each
citation and its standing.

Developed as a teaching and research tool for botany, agronomy, forestry and evolutionary biology.

## How to cite

If it helps with a thesis or a paper, please cite the software **and the original method papers it implements** —
each one is named at the point where it is used, listed on the home page, and written into the methods section the
program generates.

> Barrera-Guzmán, L.Á. (2026). *PhylogenyPro: a browser-based platform for molecular phylogenetics, from the
> alignment to the dated tree* (Version 1.4.0) [Computer software].

The machine-readable metadata is in [`CITATION.cff`](CITATION.cff) and [`codemeta.json`](codemeta.json). Once the
first release is deposited in Zenodo, the concept DOI goes into both of them and into the line above.
