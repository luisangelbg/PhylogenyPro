# Changelog

## [0.1.0] — 2026-09-22 · Block 1

First block: the home page.

- Bilingual (Spanish/English) interface with a light and a dark theme, and an original palette.
- Hero illustration, workflow strip, the twelve block cards, the data strip, the marker list by genome
  (chloroplast, mitochondrion, nuclear ribosomal, low-copy nuclear), a gallery of 58 methods, a plain comparison of
  the four inference criteria, and the measured table of sizes the app can handle.
- Two live labs that compute real phylogenetics in the page: **the Felsenstein zone** (exact pattern probabilities
  over the 15 set partitions of four taxa, the boundary of the zone by bisection, and finite-data simulation
  analysed by parsimony, neighbour joining and maximum likelihood) and **saturation and distance correction**
  (K80 with gamma rates, raw differences against the Jukes–Cantor and Kimura corrections).
- Nine sections of theory, the citation, and 94 references.
- Blocks 2 to 12 exist as panels showing what each one will do.
- `tests/index.html`: 89 unit tests, all passing.
- `validation/block1/validate_block1.R`: 17 checks against ape 5.8.1 and phangorn 2.12.1, no mismatches.
- Fixed during validation: the branch-length optimiser had a lower bound of 1e-5, which cost about 0.001
  log-likelihood units on a topology whose internal branch collapses; the bound is now a true zero.

## [0.2.0] — 2026-09-23 · Block 2

Data and alignment.

- **Readers** for FASTA, PHYLIP (strict and relaxed, sequential and interleaved), NEXUS (DATA, CHARACTERS, SETS,
  TREES, MATCHCHAR, comments), Clustal, MEGA, GenBank flat files and spreadsheets, with format and data-type
  detection and a report of everything that had to be guessed. **Writers** for FASTA, PHYLIP, NEXUS with charsets,
  RAxML partition schemes and CSV.
- **Progressive aligner**: k-mer distances, UPGMA guide tree, profile–profile alignment with affine gaps (Gotoh),
  optional iterative refinement, and codon-aware alignment. BLOSUM62 taken from the matrix distributed by the NCBI.
- **Alignment viewer and editor** on a canvas, with colouring, marks for variable and informative columns, codon
  positions, selection, and hand editing with unlimited undo.
- **Trimming**: Gblocks (Castresana 2000) and trimAl-style thresholds, with a preview of what would be removed.
- **Quality control**: repeated names, identical sequences, length imbalance, unexpected characters, gap-heavy taxa
  and columns, internal stop codons, unequal G+C, taxa missing from a partition.
- **Partitions**, each declaring its genome, with concatenation, splitting by codon position, and import and export
  of schemes.
- Example data simulated on a known tree (`data/generar_ejemplos.R`): three *Bursera* genes with three different
  genomes, a protein family and a morphological matrix, each with its true tree and true alignment.
- 170 unit tests; validation against R in `validation/block2/`.
- Fixed during validation: the PHYLIP reader preferred the strict interpretation and cut long names; the format
  detector was blind to a leading blank line, which hid a CSV.

## [0.3.0] — 2026-09-23 · Block 3

Substitution models, and the likelihood engine the rest of the app will run on.

- **Likelihood engine** (`js/like.js`): pattern compression, ambiguities as partial likelihoods, the rate matrix
  built from exchangeabilities and frequencies and diagonalised through its symmetrised form, Felsenstein pruning
  with per-node scaling, discrete gamma (Yang 1994) and invariable sites with the rates rescaled as PAML does, and
  incremental recomputation along the path to the root. **Matches phangorn to 1e-7 on nine models.**
- **Trees** (`js/tree.js`): node format shared with PopGeneticsPro, Newick reader and writer (quoted names,
  comments, support values), splits, Robinson–Foulds, rerooting by outgroup, midpoint and arbitrary branch,
  ladderising, collapsing by length or support, NJ, BIONJ and UPGMA/WPGMA.
- **Distances** (`js/dist.js`): p, JC69, K80, F81, TN93, LogDet (Lockhart) and paralinear (Lake), Poisson and
  Kimura for proteins, gamma-corrected variants, pairwise or complete deletion, saturation summary and the χ² test
  of compositional homogeneity.
- **Model catalogue and selection** (`js/models.js`): 56 DNA models, twelve empirical protein matrices, Mk and Mkv
  with the ascertainment correction of Lewis (2001), AIC/AICc/BIC with Akaike weights, hierarchical LRT, and a
  greedy partition-merging scheme.
- **Block 3 panel**: composition and saturation figures, model table, rate matrix, rate categories, the verdict in
  words, and the model handed to the blocks that follow.
- 251 unit tests; validation against R in `validation/block3/`.
- Found and fixed while validating: the Newick reader turned a `[&rate=1]` annotation into a taxon; `rerootAbove`
  left nodes of degree two behind, and `midpointRoot` placed the root on the wrong side of the branch; the
  optimiser stopped 0.2 log-likelihood units short of phangorn until a global tree rescaling was added to break the
  correlation between α and the branch lengths.
- Named honestly: what ape calls `logdet` is Lockhart et al. (1994); the formula with the observed frequencies is
  Lake's paralinear distance. The app offers both under their own names.

## [0.4.0] — 2026-09-23 · Block 4

Parsimony, distance trees and the first real tree drawings.

- **Parsimony** (`js/pars.js`): Fitch (1971) on compressed patterns with weights, Sankoff (1975) with a general cost
  matrix (ordered characters, transversion weighting), stepwise addition from several random orders, swapping by
  NNI, SPR or TBR, every equally short tree kept, CI/RI/RC/HI, Bremer support, bootstrap and jackknife.
- **Consensus** (`Consensus` in the same file): strict, majority-rule, semistrict and greedy, all built from split
  frequencies with a compatibility check.
- **Rearrangements moved into `js/tree.js`** (`nniMoves`, `sprMoves`, TBR rootings), so parsimony, minimum evolution
  and — from Block 5 — likelihood all search the same space the same way.
- **Balanced minimum evolution**: Pauplin's (2000) criterion, Desper & Gascuel (2002) branch lengths from balanced
  averages, NNI and SPR swapping, and several starting trees.
- **`Tree.unroot`**: removes the artificial root of degree two that every distance method leaves behind.
- **Tree drawing** (`js/treeview.js`): rectangular, circular and unrooted layouts, cladograms, support as numbers or
  coloured dots, scale bar, coloured clades, SVG export — the base of the tree studio of Block 12.
- **Block 4 panel**: distance trees, parsimony search with live progress, support, consensus, a Robinson–Foulds
  comparison of everything produced, figure options and exports.
- 300 unit tests; validation against R in `validation/block4/`.
- Found while validating: the NNI generator produced only two of the rearrangements per branch when a node had
  degree three, which left the minimum-evolution search in a local optimum; and `me()` searched with NNI only,
  while `ape::fastme.bal` also uses SPR.
- Reported honestly: the balanced minimum-evolution search from three starts reaches a **shorter** tree
  (2.899076) than `ape::fastme.bal` (2.901969) on the example. Both are local optima of the same criterion — the
  formula agrees to 5e-9 when applied to ape's own tree — and where a heuristic lands depends on where it starts.

## [0.5.0] — 2026-09-23 · Block 5

Maximum likelihood: searching tree space, and measuring how much the data support what is found.

- **Topology search** (`js/mlsearch.js`): NNI and SPR, with each candidate scored the way PhyML and RAxML score
  one — reoptimising only the branches the move touched, over the incremental engine of `js/like.js` — and the
  winner refitted in full. SPR is capped by a **regraft radius** and runs in two passes: every candidate within the
  radius is first scored with the lengths it inherits (the "lazy SPR" of RAxML), and only the best dozen get their
  neighbourhood reoptimised.
- **Branch support, four ways**: classical bootstrap (read from 70 % up), ultrafast bootstrap by RELL resampling of
  the per-site likelihoods of the trees the search visited (Minh et al. 2013, read from 95 % up), SH-aLRT
  (Guindon et al. 2010) and aBayes (Anisimova et al. 2011) — side by side in one table, each with its own threshold
  written in the header, because they do not mean the same thing.
- **Topology tests**: KH, SH and AU, all on RELL replicates of the per-site log-likelihoods, with AU following the
  multiscale bootstrap of Shimodaira (2002).
- **Constrained search**: force a group to be monophyletic, find the best tree that obeys the constraint, and let the
  topology test say whether the data actually reject it — instead of concluding from its mere absence.
- **Partitioned models** (`Like.fitPartitioned`): one topology and one set of relative branch lengths shared by every
  partition, each with its own substitution model and its own rate multiplier, normalised to a site-weighted mean of
  one. The block reports the multipliers and whether the BIC pays for the extra parameters.
- **Web Workers** (`js/pool.js`): the bootstrap runs on every core but one. A worker is built from the engine's own
  source text — every module of the engine is now a named function that installs itself onto whatever global it is
  handed — so it works from `file://` too, with no second implementation to keep in step and an automatic fallback
  to chunked execution on the page when workers are blocked.
- **Block 5 panel**: model and starting tree, the search with live progress and a cancel button, the partitioned fit,
  the four supports, the constrained hypothesis and the tests, the figure with any of the supports drawn on it, and
  exports (Newick, NEXUS, SVG, the bootstrap trees, the per-site likelihoods).
- 357 unit tests; validation against R in `validation/block5/`.

Found and fixed while validating:

- **The search stalled the page.** Every candidate topology built a fresh likelihood engine, allocating megabytes
  each time. `Like.makePool` now hands out buffers that candidates borrow, and the search creates them once.
- **The optimiser stopped too early.** With ten passes and a tolerance of 1e-4 the fit landed 0.0985 log-likelihood
  units below `phangorn::optim.pml` on the same tree — which looked like a numerical disagreement and was simply a
  stopping rule. With twenty passes and 1e-6 (the new defaults) it reaches −7165.395908 where phangorn reaches
  −7165.396045: a shade better.
- **A model specification was shared, not copied.** `Object.assign({}, spec)` leaves the rate and frequency arrays
  aliased, so a second analysis started from the first one's fitted values — which is exactly what made a bootstrap
  replicate computed on the page disagree with the same replicate computed in a worker. `Like.cloneSpec` fixes it,
  and the two now agree to the last digit.
- **Two root edges were counted as two branches.** In a rooted tree the two edges below a bifurcating root are the
  two halves of one unrooted branch: `Tree.internalBranches` returned both, so four of every twelve NNI candidates
  were no-ops and the branch tests reported one clade twice. `Tree.nniPartners` now swaps across that edge properly.
- **The partitioned fit landed below the single model it contains** — impossible at the optimum, and caused by
  letting every partition chase its own κ and α from the first pass. The first passes now optimise one model shared
  by all partitions (each keeping its own base composition) before releasing them, and the fit also gained the
  global tree rescaling that the single-partition fit already had.

Reported honestly: the search reaches **−7108.800752** on the *Bursera* rbcL example, better than
`phangorn::optim.pml` with NNI (−7126.989550) and a shade better than its stochastic search (−7108.802104), in 52 s
in the browser. A hundred bootstrap replicates take a few minutes across sixteen threads; this is JavaScript, not C,
and the page says so before it starts.

## [0.6.0] — 2026-09-23 · Block 6

Bayesian inference: not "the best tree" but the probability of every tree.

- **The sampler** (`js/mcmc.js`): Metropolis–Hastings over a state of (topology, branch lengths, model parameters),
  with **Metropolis coupling** (MC³, Geyer 1991) — several chains at different temperatures, only the cold one
  sampled, swapping states so the sampler can cross the valleys between islands of trees.
- **Priors written down, not assumed**: independent exponentials on the branch lengths, a uniform distribution over
  labelled topologies, exponentials on the relative rates and on the shape of the gamma, a Dirichlet on the base
  frequencies, a uniform on the proportion of invariable sites. The panel says out loud what the chosen numbers
  imply — how long the tree is *before* the data are seen.
- **Proposals with their Hastings ratios derived, not guessed**: the multiplier of Larget & Simon (1999), a
  stochastic NNI, an SPR whose Jacobian is T/(a+b) because it merges two branches and splits one, a Dirichlet
  proposal on the frequencies, and a whole-tree rescaling with its (n−1)·log m.
- **Diagnostics before the tree**: ASDSF between independent runs, the Gelman–Rubin PSRF, effective sample sizes by
  the initial-positive-sequence estimator, HPD intervals, traces and posterior histograms. The block refuses to
  present a tree as settled while any of them says otherwise.
- **Marginal likelihoods by stepping stone** (Xie et al. 2011) and Bayes factors on Kass and Raftery's scale.
- **Export**: a `.log` Tracer reads, the sampled trees as `.t`, the consensus with posterior probabilities, the
  clade table as CSV.
- Runs go through `js/pool.js`, one Web Worker per independent run.
- 400 unit tests; validation in `validation/block6/`.

The validation is of a different kind here, because two correct samplers do not give the same numbers:

- **The sampler against its own prior** (`prior_bench.html`). With the likelihood switched off the samples must
  reproduce the prior exactly, and they do: the three topologies of four taxa come out at 0.3413 / 0.3311 / 0.3277
  against 1/3; the branch lengths have mean 0.10039 and standard deviation 0.10038 against 0.1 and 0.1, median
  0.06968 against 0.06931, and the 99th percentile 0.45472 against 0.46052; the tree length has mean 0.50195 and
  standard deviation 0.22252 against the Gamma(5, 0.1) values 0.5 and 0.22361; α and κ reproduce their
  exponentials; and with six taxa **all 105 unrooted topologies are visited**, with a mean frequency of 0.00952
  against 1/105 = 0.009524. A single wrong Hastings ratio fails this test, and one did.
- **The diagnostics against `coda`** (`validate_block6.R`), on a deterministic AR(1) series that R and the browser
  build identically: the HPD interval and the potential scale reduction factor match digit for digit, and the
  effective sample size (992 in coda, theory 1053) agrees to within the difference between the two standard
  estimators.
- **The posterior against maximum likelihood**: no sampled log-likelihood above the maximum, and the 95 % HPD of κ
  and of α containing the maximum-likelihood estimates.
- **Symmetry**: on data built to be invariant under exchanging two taxa, the two topologies that map onto each
  other get the same posterior probability — 0.1751 and 0.1732.

And one finding that is not a bug but a lesson, reproduced here on purpose: with the exponential prior of mean 0.1
that MrBayes uses by default, the posterior on these very divergent simulated data comes out at half the
maximum-likelihood tree length, and κ and α go with it. That is the branch-length prior problem of Brown et al.
(2010) and Marshall (2010). The panel now compares the prior's expected tree length with the one Block 5 measured
and says so before the run rather than after.

Found and fixed while validating:

- **The incremental NNI broke the order of computation.** Scoring a rearrangement without rebuilding the tree means
  swapping two subtrees inside the flattened arrays — and the postorder that tells the engine which node to compute
  first survives that swap only when the moved subtree already happened to come before its new parent. Half the
  time it did not, and a node was combined from a child that had not been recomputed yet. Rebuilding the postorder
  after the swap costs a walk over thirty nodes.
- **Every chain computed its likelihood in another chain's memory.** The pooled buffers that made Block 5's search
  affordable are one engine's working memory, and Metropolis coupling keeps several chains alive at once: they were
  all handed the same pool. The symptom was loud once looked at — the sampled log-likelihoods came out four
  thousand units ABOVE the maximum, which cannot happen, and the substitution parameters never moved at all — and
  invisible with a single chain, which is how it survived the prior test. One pool per chain.
- **Topology moves reset the branch lengths.** The chain's lengths live in a flat array the incremental engine
  works on, and the tree object's own lengths are stale between moves. Rebuilding the tree from the node objects
  therefore threw the state away on every accepted NNI. The symptom was subtle and would have been easy to publish:
  the sampled branch lengths had exactly the right mean and far too little spread — the standard deviation of the
  tree length came out 0.161 where the prior says 0.2236.
- **A rooted tree has one branch too many.** The state is now an unrooted tree with a trifurcating root, so it has
  2n − 3 branches. With a bifurcating root the two basal branches are the two halves of one real branch, and the
  prior would have put 2n − 2 exponentials on 2n − 3 branches.

And one gap in Block 5 that only the Bayesian sampler could have exposed. The chain kept reporting log-likelihoods
above the "maximum" — thirty units above — which is impossible, and the first explanation tried here (a trapped
optimiser) turned out to be wrong. The real reason is simpler: the sampler ESTIMATES the base frequencies, while
`Like.fit` holds them at the empirical proportions, as `phangorn::optim.pml` does unless optBf is turned on. On the
*Bursera* example the fitted frequencies are (0.271, 0.231, 0.249, 0.249) against the observed (0.266, 0.202,
0.274, 0.258), and that difference is worth those thirty units. Block 5 now offers **+F** — estimate them — as a
checkbox, off by default so that the agreement with phangorn documented in Block 3 still holds.

`Like.fit` also gained a short sweep over starting values for the gamma shape and the exchangeability rates before
the first branch pass. That was written while chasing the wrong explanation, but it earns its place: without it the
fit of one topology could land in different places depending on where it started, and with it two very different
starting models converge to the same optimum.

Speed, honestly: the four-state kernel of the likelihood was written out by hand for this block — DNA is the case
that matters — and the NNI is now scored incrementally instead of rebuilding the tree. Block 5's own benchmark,
rerun unchanged afterwards, went from **77.8 s to 7.9 s** for the same search reaching the same −7108.800700, and
the agreement with phangorn is untouched (−7165.395908 against −7165.396045). A hundred thousand generations of
coupled chains on 18 taxa and 621 sites takes minutes, not the hour it took before. It is still JavaScript and not
C, and the panel says what it is about to cost before it starts.

## [0.7.0] — 2026-09-23 · Block 7

The molecular clock, and turning substitutions into years.

- **The order is the method** (`js/block7.js`): root first, then find out whether the clock holds, then write the
  calibrations down, and only then date. Any other order hides a decision.
- **Rooting** (`js/clock.js`): by outgroup, by midpoint, and — when the tips carry dates — by root-to-tip
  regression, searching every branch *and* the position along it, which is what `ape::rtt` does and what makes the
  difference between reproducing it and missing by half the r².
- **The strict-clock test**: the free tree against the ultrametric one, with n − 2 degrees of freedom. The panel
  says what rejection means, which is that a relaxed clock is needed, not that dating is impossible.
- **Calibrations** as fixed ages, uniform intervals, and normal, offset lognormal and offset exponential densities —
  the shapes a fossil actually justifies, since a fossil gives a minimum and almost never a maximum.
- **Three ways to date**: least squares in the manner of To et al. (2016), penalised likelihood (Sanderson 2002,
  the method behind `ape::chronos`) with the smoothing chosen by cross-validation, and a **relaxed-clock MCMC** with
  uncorrelated lognormal rates (Drummond et al. 2006) that gives a posterior and a 95 % interval for every node.
- **The chronogram** (`js/plots7.js`): time running backwards, the geological periods and epochs drawn behind the
  tree from the ICS chart (`js/geotime.js`), and the uncertainty of every node as a bar. A chronogram without those
  bars claims a precision no dating method has.
- **Without a calibration** the block still dates the tree, with the root at 1, and says in as many words that the
  numbers are proportions of the depth of the tree and not years.
- 466 unit tests; validation in `validation/block7/`.

Validated against ape and phangorn:

- **Rooting by tip dates reproduces `ape::rtt`**: slope −0.00158106 against −0.00158106 (a difference of 3·10⁻⁹),
  intercept 0.19218168 against 0.19218230, r² 0.783919 against 0.783919.
- **Penalised likelihood against `ape::chronos`** on the same rooted tree with the same λ and the same calibration:
  15 clades in common, largest difference in age 0.0382 of a root fixed at 1, mean 0.0125, correlation r² 0.9906.
- **The clock test against `phangorn::optim.pml(optRooted = TRUE)`**: the free fits agree to 0.005; the ultrametric
  fits do not, and the app's is **better by 11.8 log-likelihood units**, so its statistic is smaller (49.12 against
  72.81). Both reject the strict clock on these data at any threshold, which is the conclusion that matters.
- **The identities least squares has to satisfy**, checked exactly: the root lands on its fixed age to nine
  decimals, every tip at zero, no node older than its parent, doubling the root age doubles every age and halves
  the rate, and an interval calibration is honoured to the boundary.
- **The calibration densities** against their closed forms, and the geological chart against the ICS boundaries.

Found while validating: the search over rootings was filtering candidates by the *sign* of the regression slope,
on the assumption that root-to-tip distance grows with the date. It does when the dates are calendar years and
falls when they are ages before the present, so the filter threw away every correct rooting and kept noise — r²
0.006 where ape reaches 0.784. The sign now only reports which convention the dates follow; the rate is the
magnitude of the slope.

## [0.8.0] — 2026-09-23 · Block 8

Diversification: reading the shape of a dated tree.

- **Lineages through time** on a logarithmic axis, with the straight line a constant rate would give drawn behind
  it, because the curve says nothing without that reference.
- **The γ statistic** of Pybus & Harvey (2000), in exactly the form `ape::gammaStat` computes it, with the reading
  spelled out: a negative γ is a slowdown *or* missing species, and nothing in the number itself tells which.
- **The Monte Carlo constant-rates test**, which is what settles that: the null distribution is built by simulating
  trees of the clade's *true* size under a constant rate and pruning them at random down to the species in hand.
  On the example, 18 species of 120 give a mean γ of −2.21 with nothing having slowed down at all.
- **Five models on the same branching times** — pure birth, birth–death, logistic and exponential density
  dependence (Rabosky & Lovette 2008), and two rates with the moment of the shift estimated — compared by AIC with
  Akaike weights, and a warning that a density-dependent model winning does not prove diversity is bounded.
- **Magallón & Sanderson (2001)** for the case where only a clade's age and richness are known, reported at three
  assumed extinction fractions, because that assumption moves the answer by half.
- **The DR statistic** of Jetz et al. (2012), a rate for every species.
- **Simulation** of reconstructed birth–death trees, which the Monte Carlo test runs on.
- 525 unit tests; validation in `validation/block8/`.

Validated against ape and geiger, digit for digit:

- **γ**: 0.44734655 against `ape::gammaStat`'s 0.44734655.
- **Pure birth**: λ 0.24064051, standard error 0.03903705 and log-likelihood 14.502617, all identical to `ape::yule`.
  The same λ comes out of the waiting-time formulation the density-dependent models use, which ties the two
  parameterisations together.
- **Birth–death**: a = 0.24684678 against `ape::birthdeath`'s 0.24684645, r = 0.20638381 against 0.20638489, and the
  same log-likelihood 14.620461. On a tree simulated without extinction both implementations find exactly zero
  extinction and a net rate equal to the Yule λ.
- **Magallón & Sanderson**: 0.23604983, 0.21510392 and 0.12113234 for crown ages at ε = 0, 0.5 and 0.9, identical to
  `geiger::bd.ms` to 1e-7, and the stem cases too.
- **Branching times** against `ape::branching.times` to 1.2·10⁻⁹.

Found while validating: the simulator kept the extinct lineages in the tree, so it returned more tips than it was
asked for and trees that were not ultrametric — and then it stopped the clock at the instant the last lineage
appeared. That second mistake was the dangerous one: the present is somewhere *after* that split, and putting it
exactly on it squeezed every node towards the tips, so γ over simulated constant-rate trees averaged **+0.32**
instead of 0. Every Monte Carlo test built on that null would have been wrong in the same direction. With the
present drawn properly, γ over 300 simulated pure-birth trees has mean −0.066 and standard deviation 0.9999 — the
standard normal the theory says, and the same as `ape::rphylo` gives (0.073, 0.952).

## [0.9.0] — 2026-09-23 · Block 9

Trait evolution: characters read on a tree, so that shared history is not mistaken for evidence.

- **A character table of its own**: one row per species and one column per character, in commas, semicolons or
  tabs, with the type guessed per column (all numbers means continuous), unknowns as `?`/`NA`/`-`/empty, and
  names matched to the tree forgiving the difference between `Bursera_palmeri` and `Bursera palmeri`. Whatever
  does not line up is listed, not swallowed; each analysis prunes the tree to the tips that have the character
  it needs and says how many that left.
- **Three continuous models** — Brownian motion, Ornstein–Uhlenbeck and the early burst — in geiger's
  parameterisation, compared by AIC with Akaike weights, plus **ancestral states** and the **traitgram** that
  puts the character and the tree in one plane.
- **Phylogenetic signal** both ways: Blomberg's K with a permutation test and Pagel's λ with a likelihood ratio,
  and the block explains why K < 1 with λ ≈ 1 is not a contradiction.
- **Two characters**: independent contrasts, PGLS and the uncorrected regression side by side, with the identity
  of Garland & Ives (2000) checked in the interface — the contrasts slope through the origin and the PGLS slope
  have to be the same number, and when they are not, something is wrong.
- **The Mk model** of Lewis (2001) with ER, SYM and ARD over a general rate matrix, **marginal ancestral states**
  drawn as pies so a node at 55 % looks like one, and the downward pass available separately for anyone who
  wants to see where the reconstruction comes from.
- **Stochastic character mapping** (Huelsenbeck et al. 2003): 500 histories in under a second, turning node
  probabilities into how much of the tree each state actually held and how many changes that took.
- 628 unit tests; validation in `validation/block9/`.

Validated against geiger, ape, nlme and phytools, digit for digit — 42 checks in `validate_block9.R` and every
comparison in `traits_bench.html`, with no mismatch:

- **BM** σ² 0.06350382, z₀ 0.25281475, lnL −26.628333; **OU** σ² 0.08197347 and α 0.06985066 against geiger's
  0.08197349 and 0.06985065; **EB** at its bound, exactly where `fitContinuous` leaves it.
- **The 39 contrasts** of both characters against `ape::pic` to 1.7·10⁻¹⁴, and their sum of squares over n equal
  to the Brownian rate.
- **PGLS** against `nlme::gls` with `corBrownian`: intercept 0.11628292, slope 0.00796916, standard error
  0.12963002, t 0.06147623, p 0.95130215 and lnL −17.657220 — and the slope agrees with the contrasts
  regression to 4·10⁻¹⁷.
- **Signal**: K 0.90386290 and λ 0.97380521 against `phytools::phylosig`'s 0.90386290 and 0.97380668.
- **Mk**: rate 0.02920872 and lnL −11.307863 in `ape::ace`'s convention, which is **`phytools::fitMk`'s plus
  log(k)** — the whole difference between the two is the flat prior at the root, and the block reports both.
- **Ancestral states**: the marginal matches `ape::ace`'s `lik.anc` over all 39 nodes to 5.7·10⁻⁸ and
  `phytools::ancr` to 3.5·10⁻¹⁰; the continuous ones match `ace(method = "ML")` to 2.8·10⁻⁷.
- **Stochastic mapping**: 129.43 and 28.48 units of branch time in the two states against `make.simmap`'s 129.73
  and 28.19, and 4.60 changes per history against 4.65, from independent draws.

Found while building and validating:

- `(exp(a·s) − 1)/a`, the early burst's covariance, is the textbook cancellation: for small |a| it loses every
  significant digit, and a model that is supposed to become Brownian motion as a → 0 instead drifted away from
  it. `expm1` is exactly the function that does not lose them.
- The Student t distribution function was built on the incomplete beta at `df/(df + t²)`, which only sees t², so
  it came out symmetric about zero — a distribution function that is its own reflection. The p-values were right
  because they were always called with |t|, but the function was not; the sign is now put back by hand.
- **α in the OU fit was bounded between fixed numbers**, 10⁻⁶ and 20. α is a rate, so its useful range is set by
  the depth of the tree: on a phylogram 0.08 substitutions deep an α of 20 is a pull nothing can see, and the fit
  reported a half-life of 691 179 substitutions on a tree a tenth of a unit long. The bounds are now αT ∈ [10⁻⁶,
  200], as geiger does, and when α lands on one the tile says so instead of printing the number.
- **OU and the early burst assume the tips are contemporaries** and nothing said so. On the app's own distance
  tree — tips varying by 0.29 over a depth of 0.392 — the early burst wins with weight 0.92 and means nothing.
  The block now measures the spread and warns before the table is read.
- `--ink`, the token Blocks 7 and 8 draw their curves with, **was never defined in the stylesheet**. An undefined
  custom property makes the stroke invalid, so those curves fell back to black and vanished into the dark theme.
  It is now declared in both themes.

## [0.10.0] — 2026-09-23 · Block 10

Historical biogeography: where the lineages were, and what happened to a range when a lineage split.

- **Ranges as a character whose states are subsets of areas**, with the null range as an absorbing state, a cap on
  the range size (the only thing that makes more than six areas tractable, and a real assumption, said out loud)
  and a matrix of **dispersal multipliers** the user edits: that is where a land bridge that opened, or an ocean
  that did not close, enters the model.
- **The four models and their +J versions** — DEC (Ree & Smith 2008), DIVALIKE, BAYAREALIKE, each with and without
  Matzke's (2013) founder-event jump dispersal — fitted by maximum likelihood and compared by AICc with Akaike
  weights and a likelihood ratio against the model each +J contains.
- **The warning that has to be there**: whenever a +J model wins, the block says why that is the expected outcome
  and not evidence on its own. Ree & Sanmartín (2018) showed the comparison is not clean — j buys probability at
  cladogenesis, an instantaneous event costing nothing along the branches, while d and e pay for theirs by
  integrating over the whole time of the tree.
- **Ancestral ranges**, marginal, with the **corners**: the range each daughter starts with immediately after the
  split, which is the half of the story a node-only reconstruction leaves out. Drawn as pies over the ranges that
  have any probability, with a grey slice for everything below the threshold.
- **Biogeographic stochastic mapping**: whole histories drawn conditioned on the tips, turning probabilities into
  counts of dispersals, local extinctions and the four kinds of split.
- **Time strata**: a branch that crosses a boundary gets the product of the matrices of the periods it crosses,
  each with its own multipliers and its own set of areas that exist at all.
- Reads and writes the LAGRANGE range format that BioGeoBEARS also uses, and a plain comma-separated table.
- 736 unit tests; validation in `validation/block10/`.

Validated against **BioGeoBEARS 1.1.3** on the *Psychotria* of Hawaii — the data set LAGRANGE and BioGeoBEARS are
both distributed with, and the one every paper about DEC+J re-analyses. 41 checks in `validate_block10.R` and every
comparison in `biogeo_bench.html`, with no mismatch:

| model | PhylogenyPro | BioGeoBEARS |
|---|---|---|
| DEC | −34.54195625 | −34.54195758 |
| DEC+J | −20.94758871 | −20.94758855 |
| DIVALIKE | −33.14967582 | −33.14967538 |
| DIVALIKE+J | −21.08620942 | −21.08620912 |
| BAYAREALIKE | −40.33442556 | −40.33442795 |
| BAYAREALIKE+J | −21.55262993 | −21.55262941 |

- **d and e** agree to about 1e-5 (DEC: 0.03505 and 0.02831 against 0.03505 and 0.02836).
- **The ancestral ranges** agree to **1.75·10⁻⁸** over all eighteen internal nodes when both are evaluated at the
  same parameters; the root comes out in KO at 0.51655966 against 0.51655968. At PhylogenyPro's own optimum the
  largest difference rises to 4.3·10⁻⁴, which is what a difference of 1e-5 in d costs — worth knowing, and the
  bench reports both.

Found while building and validating:

- **The reference was wrong before the code was.** DIVALIKE disagreed by 0.033, and the cause was in the R script:
  BioGeoBEARS reads the parameter `mx01v`, not `maxent01v`. Setting the second one does nothing at all and leaves
  DIVALIKE without its 2–2 vicariance splits, which is not DIVALIKE. The cladogenesis weight tables were then
  extracted from `cladoRcpp` directly (`validation/block10/pesos_cladogenesis.R`) and match event for event.
- **The upward message in the ancestral reconstruction included the daughter's own likelihood**, counting it
  twice. The symptom is the one to remember: reconstructions far *too confident* and biased towards widespread
  ranges — one node at 0.96 for KOMH where the truth was 0.11. The message that comes down to a daughter has to
  carry everything *except* that daughter's subtree; the corner posterior, which is for reading, uses all the data.
  They are two different quantities and the code now computes both.
- **The root convention** is a choice, not a fact: LAGRANGE and BioGeoBEARS weight every allowed range by one,
  which is what every published DEC log-likelihood means, and a proper flat prior differs by exactly log(number of
  ranges) — 2.708050 here. The same constant for every model over the same areas, so no comparison changes, but
  the number a reader checks against a paper does. The app reports BioGeoBEARS' convention and says so.
- **Speed, three times.** Fitting one model took 7.4 s. Reusing the buffers inside the matrix exponential instead
  of allocating a matrix per Taylor term took it to 0.3 s — the garbage collector had been doing most of the work.
  Tying the optimiser's starting values to the depth of the tree, instead of the fixed 0.03 that only suits a tree
  a few million years deep, removed most of the remaining iterations. And replacing the per-branch exponential
  with **uniformization** — P(t) = Σ Poisson(k; μt)·Rᵏ, where the powers of R are built once and shared by every
  branch — took it to **0.1 s, some seventy times the original speed**, with every digit unchanged. Every term of
  that sum is non-negative, so it cannot produce the slightly negative probabilities the Taylor series could.
- The rejection sampler for a branch's history wasted most of its draws when the two ends differed, because the
  commonest forward history is the one where nothing happens. Nielsen's (2002) modification — draw the first
  waiting time from the exponential truncated at the branch length when the ends differ — cut the exhausted
  rejections from 5 % to under 1 %.

## [0.11.0] — 2026-09-23 · Block 11

Comparing trees, and admitting they disagree.

- **Five ways of saying how far apart two trees are**, in one table and one heat map: Robinson–Foulds raw and
  normalised, the weighted version, the branch score of Kuhner & Felsenstein, the path difference of Steel & Penny
  and the **quartet distance**. The block says beside them why the last one degrades gracefully and the first does
  not: one misplaced tip can send RF to its maximum while the quartet distance barely moves, and when that happens
  the disagreement is a tip, not a history.
- **Tanglegrams**, with the untangling step that makes them readable — rotating branches, which changes neither
  tree, until the lines cross as little as possible — and a warning nobody else prints: **a tanglegram with no
  crossings does not mean the trees are the same**, because a rotation can align two different topologies.
- **Concordance factors** (Minh et al. 2020): gCF from a set of gene trees and sCF from the alignment, branch by
  branch, as stacked bars against the third that pure chance would give, and as a gCF-against-sCF scatter. The
  block keeps repeating the point of them: **concordance is not support**, and a branch at 100 % bootstrap and
  35 % concordance is well estimated and badly supported at once.
- **A species tree from the quartets**: the frequency of each of the three resolutions of every quartet across the
  gene trees, the score of each candidate tree, and a nearest-neighbour search — a heuristic, named as one, not
  ASTRAL's exact program.
- **Split networks**, two ways: a **consensus network** of the splits enough trees have, with the incompatible
  pairs found and drawn in colour, and the **circular decomposition** of a distance matrix, which on a tree's own
  distances gives that tree back exactly and on anything else shows the difference.
- **ABBA-BABA**: Patterson's D with a block jackknife, for one trio or for every trio against an outgroup, with
  the warning that **D ≠ 0 has three explanations** and only one is introgression.
- 829 unit tests; validation in `validation/block11/`.

Validated against ape 5.8.1 and phangorn 2.12.1 — 25 checks in `validate_block11.R` and every comparison in
`cmp_bench.html`, with no mismatch:

- **The five distances over 36 pairs each**, 216 numbers: RF and normalised RF exact, weighted RF to 4.5·10⁻¹⁰,
  the branch score to 1.5·10⁻¹⁰, the path difference to 5·10⁻¹⁴ and its weighted form to 8·10⁻¹⁰.
- **The quartet distance** against a brute-force count over all 1820 quartets of sixteen taxa, resolving each one
  by pruning both trees to it: 1767/53, 1691/129 and 1820/0, exact.
- **gCF** against the same definition recomputed in R, over all thirteen internal branches, to zero difference.
- **The consensus network** split for split against `phangorn::consensusNet`, the thirteen non-trivial splits
  identical.
- **Patterson's D** against the definition computed in R on a 60 000-site alignment with 15 % introgression:
  ABBA 1089, BABA 753, D 0.18241042, jackknife standard error 0.07884390 and z 2.31356419, every digit.
- **The circular decomposition** against an identity rather than a package, because none is installed for it: fed
  the path distances of a sixteen-taxon tree it returns exactly that tree's 13 internal and 16 terminal splits,
  with its branch lengths to 3·10⁻¹⁶ and the distances reconstructed to 9·10⁻¹⁶, nothing clipped; and the outer
  cycle of the network closes on itself to 3·10⁻¹⁶.

Found while building and validating:

- **The weighted measures were reading half of one branch.** A rooted tree has two edges under its root and they
  are *one* edge of the unrooted tree; taking the length from whichever node the split map happened to keep
  halved it. Every weighted comparison inherited the error — weighted RF came out 2.54 where phangorn says 3.78.
  The split table now accumulates lengths by split instead of reading them off a node.
- **The path difference counted the root.** In edges, the path between two tips depends on where the root sits, so
  two identical trees rooted differently came out 7.94 apart instead of 0. The measure is about the unrooted tree,
  and now the root goes first.
- **A tree rooted on a single tip grew an extra internal branch.** The branch above its sister holds n−1 tips, and
  the test for "trivial" only looked at the side below the node. One of the app's own example trees showed 16
  internal splits where 15 exist, which inflates every comparison with a tree rooted elsewhere by one. Found by
  running the block on real output, not by the bench, whose reference trees never rooted that way.
- **The four groups around a branch were three under the root.** A branch directly under a bifurcating root has
  nothing "above" it — the root is not a node of the unrooted tree — and the division happens at the sister's own
  daughters instead. Missing that case dropped exactly one branch, always the deepest, from every concordance
  table: twelve branches where thirteen exist.
- **The circular decomposition had its four terms in the wrong pairs.** The weight of the split {xᵢ … xⱼ} is a
  four-point difference across the two gaps that bound it, and pairing the distances by gap instead of by endpoint
  gives a number that looks like a weight and is not one: 91 splits where a sixteen-taxon tree has 29, and
  distances it could not reproduce. The corrected formula returns the tree exactly.
- `opts.minSites || 20` turned a deliberate zero into twenty, so asking the ABBA-BABA scan for every trio returned
  none.

## [1.0.0] — 2026-09-23 · Block 12 · the twelve blocks complete

Figures, formats and the report: handing the work over.

- **A tree studio** where every decision a tree figure contains is visible and reversible — layout, cladogram or
  phylogram, which support is shown and above what threshold, where the root goes, how the branches ladder, the
  fonts, the widths, coloured clades — with the drawing redrawn at every change.
- **Three tree formats**, each for what it can carry: **Newick**, universal and narrow; **NEXUS** with a TRANSLATE
  block and the block FigTree reads; and **phyloXML**, the only one of the three that can hold a node's support,
  its age and its confidence interval at once, which is exactly what Blocks 5, 6 and 7 produce.
- **Figures at any resolution**, with the arithmetic shown first: a PNG at 900 dpi is a particular number of
  pixels, and the block prints that number before generating anything and refuses — with a reason — the sizes a
  browser canvas cannot hold, instead of returning a blank image. The SVG has no such limit and is what to send.
- **A methods section written from the run.** A block that was never opened leaves no paragraph; a support that
  was not computed is not mentioned; the number of replicates comes from the run and not from a template. The
  bibliography is assembled from what was cited and holds nothing else — 52 works the program knows, and only the
  ones it used appear. The report is one self-contained HTML file: no stylesheet to fetch, no font to download, no
  script to run, and it prints to PDF from the browser.
- **A reproducible package**: the data as they came in, the trees in all three formats, the figures as SVG with
  their colours resolved, the tables as CSV, the report, and an index saying what each file is and how to open it
  again. The ZIP is written here, byte by byte, with no library and no network, so it works from a double-clicked
  page.
- 910 unit tests; validation in `validation/block12/`.

Validated by making R do what a reader would do — 38 checks in `validate_block12.R` and every check in
`export_bench.html`, with no mismatch:

- **The ZIP**: R's own `unzip()` lists and extracts every file, which verifies each CRC32 on the way, and the
  extracted sizes are the declared ones.
- **Newick and NEXUS**: both read by `ape::read.tree` and `ape::read.nexus`, the same 18 tips, topological
  distance 0 between them and branch lengths identical to 10⁻⁶.
- **phyloXML**: parsed by `xml2`, its tip names are the tree's, and its branch lengths sum to what the Newick sums
  to (2.446035 against 2.446036).
- **The alignment**: read by `ape::read.dna` at 18 × 1711, with the tree's taxa.
- **The report**: one file with no `<link>`, no `<script>` and no `<img>`, carrying its figures as inline SVG —
  and, in a package built from a distance tree alone, citing no clock, no biogeography and no MCMC, because none
  of them ran.
- **The figures**: well-formed XML with no `var(--…)` left unresolved, which would be a black line in Illustrator.
- **CRC32** against the value the ZIP specification names for "123456789": 0xCBF43926.

Found while building and validating:

- **The report module depended on the page's translator**, which broke it outside the app and, worse, produced a
  half-translated document: the fragments it assembled — "non-parametric bootstrap (100 replicates)", "by
  outgroup", "penalised likelihood" — were resolved once, in whatever language the interface happened to be in,
  and then pasted into *both* the Spanish and the English paragraph. Every such fragment is now kept as a pair,
  and the module no longer knows what language the page is showing.
- The genome of a partition is stored as a two-letter code, which is fine in a menu and absurd in a methods
  section: "the genomes represented are: nu".

---

With Block 12 the twelve blocks are complete. The app takes a set of sequences from a FASTA file to a dated tree,
its ancestral traits, its ancestral areas and a reproducible package, entirely in the browser, and every block is
checked against an independent implementation in R. The totals, across the twelve validation folders: **910 unit
tests**, and **258 numerical agreements with no mismatch** in the eleven blocks that report a count, against ape,
phangorn, geiger, phytools, nlme, BioGeoBEARS and xml2. Block 6 is not counted the same way because it is checked
in three other ways instead: its sampler reproduces its own prior exactly with the likelihood switched off, its HPD
and PSRF are identical to `coda`'s, and an audit compares the likelihood the chain carries with a calculation from
scratch at every move, finding a mismatch of zero.

## [1.0.1] — 2026-09-24 · Provenance

No change to what the program computes. What changed is what it says about where its constants come from, before
the work is filed with INDAUTOR.

- **The twelve empirical amino-acid matrices now cite their own articles.** Until now the header of
  `js/aamodels.js` said the values had been pulled out of phangorn's internals, which is what
  `data/extraer_matrices_aa.R` had done. That was a true statement about how the file was built and a false one
  about where the numbers come from: they are the matrices published by Le & Gascuel, Whelan & Goldman, Jones,
  Taylor & Thornton, Dayhoff, Adachi, Müller & Vingron, Dimmic, Henikoff & Henikoff, Dang, Nickle and Yang.
  Each of the twelve now carries its own citation beside it in the file.
- **`data/extraer_matrices_aa.R` is gone**, replaced by `data/verificar_matrices_aa.R`, which extracts nothing and
  verifies everything: it reads the canonical distributions — PAML's `dat/`, IQ-TREE's `modelprotein.cpp`,
  MrBayes' `model.c`, RAxML's `models.c`, and the WAG file Whelan and Goldman themselves put on the EBI server —
  transposes each lower triangle into the order this program uses, and compares every value.
  **28 contrasts, 5806 comparisons, 5242 identical to the bit, and not one of the 2280 exchangeabilities
  disagreeing.** The eight differences it does find are all equilibrium frequencies, and it prints them one by one
  instead of hiding them in a tolerance.
- Three of those eight are worth recording, because in each the program is right and the secondary distribution is
  not: IQ-TREE's Dayhoff normalises H and K from a copy whose last published digit differs from PAML's, and
  MrBayes writes 0.168 for leucine in mtREV24 where PAML and Adachi & Hasegawa write 0.169 — with 0.168 its twenty
  frequencies sum to 0.999. The other five are FLU frequencies differing from IQ-TREE's in the seventh decimal,
  below the precision at which FLU was published.
- `data/matrices_fuente/descargar_fuentes.R` fetches those witnesses on demand, and `LEEME.md` says what each one
  is and that none of them is part of the program.
- **New `PROCEDENCIA.md`**: every element of the program that is not original — the twelve matrices, BLOSUM62, the
  NCBI genetic code, the ICS time scale — with its source, its citation and its standing, plus the reasons the
  work is primigenia and what belongs in the deposited copy.
- `js/align.js` and `js/geotime.js` state that the NCBI material is a work of the US government in the public
  domain and that the ICS chart is published for free use with the citation the program already gives.
- The 910 unit tests still pass and the app loads the twelve matrices unchanged.
- The GPL-3.0 the README had always declared now has the file that makes it real: `LICENSE`, the Free Software
  Foundation's text verbatim, and the copyright line in the README. The licence is a choice and not an
  inheritance — the program contains no third-party code — and `PROCEDENCIA.md` says so, along with the fact that
  licensing is not assigning authorship: moral rights are inalienable under Mexican law, and the registration and
  the licence answer two different questions.
- **The notice is now at the head of every source file of the program**: the 46 files in `js/`, `css/style.css` and
  `index.html`, 48 in all, each carrying the copyright line and the licence paragraph above the descriptive header
  it already had. `data/hacer_examples.sh` emits it too, so rebuilding `js/examples.js` no longer drops it.
  The Web Workers are unaffected: `js/pool.js` builds them from the modules' own `toString()`, which is function
  source and not file source, so `Pool.sources()` still begins at `function TreeCore(g) {`. Checked by running six
  real parsimony bootstrap replicates through the pool, and the 910 unit tests still pass.
- **The author's surname was wrong throughout, and is now right.** The citation the home page has shown since Block
  1 read "Barrera-García"; the ten sibling applications and the ORCID record all read **Barrera-Guzmán**. The
  mistake had been copied into the copyright line of all 48 source files when the licence notices went in. Every
  occurrence in the program, the README, PROCEDENCIA, the changelog and `data/hacer_examples.sh` — 51 in 50
  files — now reads Barrera-Guzmán.
- **`CITATION.cff` and `codemeta.json`**, in the form the other LABG applications use: version 1.0.1, GPL-3.0-or-
  later, ORCID 0000-0001-8057-2583, twelve keywords, and an abstract that says what the twelve blocks do. Both
  parse cleanly and carry every required field. There is no DOI yet: the CFF has a commented line marking where the
  Zenodo concept DOI goes once the first release is deposited, and the README says the same. The repository and
  Pages URLs follow the convention of the sibling apps and go live when the repository is created.
- A **How to cite** section in the README, which it did not have.
- **Fixed: the partition menus showed both languages glued together** — "nuclearnuclear", "cloroplastochloroplast",
  "ADNDNA". The genome and data-type `<option>`s were built with `L2()`, which writes one `<span>` per language and
  lets CSS hide the inactive one. That works anywhere else in the page, but an `<option>` may not contain elements:
  the browser keeps the text of both spans and concatenates it. Those two menus now carry the two texts as
  `data-es` / `data-en` attributes, which is what `I18N.apply()` already reads, and are written in the language on
  screen when they are rendered. They were the only two places in the program where `L2()` was used inside an
  `<option>`.

## [1.1.0] — 2026-09-24 · The tree studio grows up

Pictures of the taxa, groups, and labels that stop sitting on top of each other.

- **`js/otuimg.js`** — pictures of the OTUs and of the clades. Import by drag and drop (a file named after its taxon
  finds its own place), a square preview, zoom, pan, rotation, brightness, contrast, saturation, warmth, sharpness,
  greyscale, and a **background remover** that floods in from the border so light areas *inside* the object survive,
  which a threshold on brightness would destroy. Credit and licence per picture. Images live in IndexedDB, keyed by
  data set, and never leave the computer. **Nothing is drawn until the user ticks "I confirm this picture is of X"**:
  a photograph pinned to the wrong branch is a claim, not a decoration. Sister module of PopGeneticsPro's, by the
  same author.
- **`js/groups.js`** — cutting a tree into k groups, three ways: by the longest branches (any tree), by depth (only
  on an ultrametric tree — it says so and falls back rather than lying), and from clades or by hand. Reports the
  sizes and **whether each group is actually a clade**, because colouring a paraphyletic group and saying nothing
  is how a figure misleads.
- **Tree studio**: colour by group on branches, names, a colour strip and soft bands; a legend; group names the user
  can set; aligned tip names with a dotted leader; pictures per tip, per group or per marked clade, at any size, as
  circles, soft-cornered or square, framed or not. The pictures travel **inside** the SVG as data URLs, so the
  exported figure depends on no loose files — checked by exporting to PNG and finding the photograph's pixels in it.
- **Labels that keep apart.** Two new helpers in `js/treeview.js`: `spreadLabels` nudges positions along one axis
  until none is closer than a given gap, staying inside the frame and reporting when there is simply no room;
  `dropColliding` keeps the most important of a set of boxes and drops the rest. Support values on a crowded tree
  no longer print on top of each other, and **the panel says how many were hidden and what to change** instead of
  silently losing them. The **traitgram of Block 9** was the worst offender — every tip labelled at whatever value
  it reached — and now spreads its names vertically with a thin leader back to the point: 25 tips crammed into a
  range of 0.1 draw all 25 names with no overlap at all.
- Two bugs found and fixed while testing the cut into k groups, both of them mine, and both with the same shape —
  a branch counted twice: **the two branches under a bifurcating root are one branch of the unrooted tree**, so
  cutting both bought one group and paid for two (asking for 3 returned 2); and **cutting both branches of an
  internal node leaves it holding no tips**, a group that exists in the graph and not in the data (asking for 5
  returned 4). Cuts are now taken one at a time, longest first, and one that does not actually add a group is
  skipped rather than wasted. Asking for k now gives k, for every k up to the number of tips, and says so plainly
  when it cannot.
- Also fixed: the legend was drawn over the scale bar, and a band for a paraphyletic group spanned other groups'
  tips — it would have claimed they belonged to it. A band is now drawn only for a group whose tips are
  consecutive in the figure.
- `tests/studio.js`: **60 new unit tests**, total **970**, all passing.
- **Circular and unrooted now do everything the rectangular layout does**: branches and names coloured by group, a
  ring of group colour (the circular answer to the strip), a legend, aligned names with a dotted leader, and the
  pictures. In the unrooted layout nothing keeps two tips from landing beside each other, so the names that still
  collide after placement are dropped, the most peripheral kept, and the count reported — `labelDeclutter: false`
  turns that off for anyone who wants them all.
- Found while drawing them: reserving room for the pictures by **taking it out of the radius** shrank the tree
  towards nothing as soon as the names were long — at 430 px wide with names like *Bursera copallifera* the whole
  tree collapsed into a seven-pixel knot at the centre. Both radial layouts now put the ring and the pictures on a
  **wider canvas** instead: the width asked for is the width of the tree and its names, and whatever is drawn
  around them is added to the figure.
- Also found: the dotted leader tying a picture to its branch always ended at a point to the left of the picture.
  That is right in the rectangular layout and nonsense in a radial one, where a picture at the top of the circle
  was pointed at from its side. It now stops at the edge of the picture on the side the branch is really on.
- `tests/studio.js` grew to **76 tests**; total **986**, all passing.
- **Audit of every other figure for overlapping labels.** Each of the seven plot modules was gone through, and every
  label whose position comes from the data — rather than from a row or an axis tick — was checked with crowded
  input. Two more real cases, and a good deal of false alarm:
  - **The chronogram of Block 7.** A balanced tree of 64 tips at a row height of 7 printed **62 colliding pairs of
    node ages**. The most inclusive clades now keep their date, the rest are dropped, and the panel says how many
    and where to find them — they are all in the table and in the package regardless.
  - **The gCF–sCF scatter of Block 11**, when branch labels are switched on: sixty branches clustered near the
    diagonal gave **1319 colliding pairs**. The branches furthest from the diagonal are kept first, since a branch
    whose gene and site concordance disagree is the one worth naming. This one is a latent fault rather than
    something anyone was seeing: nothing in the app turns `labelAll` on today.
  - **Everything else was fine, and for a reason worth recording**: the tip names of the pies tree, of the
    biogeography tree and of the chronogram are placed by row, so they cannot collide; the axis ticks come from a
    generator that aims for five to seven per axis with the decimals tied to the step, so they cannot either. A
    test now pins that last one down so it stays true.
- `tests/studio.js`: **87 tests**; total **997**, all passing. Several of them measure the collisions on the SVG
  itself, before and after, so the numbers above are the test's own output and not an impression.

## [1.2.0] — 2026-09-24 · Reading the numbers

A number a reader cannot place is not a result. **`js/help.js`**: 46 entries covering the twelve blocks, each with
what the thing measures, how it is read, the scale that says whether a value is low or high, and the mistake most
commonly made with it — in both languages, with the colours taken from the stylesheet so both themes follow.

- **The scale is applied to the number, not only explained.** `Help.band(key, value)` says which band a computed
  value falls into and `Help.tag` prints it as a coloured chip, so the reader does not have to carry the thresholds
  in their head. A scale printed in a manual is documentation; a scale applied to what is on screen is a decision
  rule, which is what was actually needed.
- **Where a threshold is a habit, the entry says so.** Half the cut-offs in phylogenetics are reading conventions
  rather than statistical laws, and each one is labelled as such.
- The entries make a point of the distinctions that are most often got wrong: that **70 % bootstrap and 95 % UFBoot
  mean roughly the same thing** and confusing the two scales is the commonest error with ultrafast bootstrap; that a
  posterior probability and a bootstrap value are not on the same scale; that **a branch can have 100 % bootstrap
  and 35 % gCF**, which is the signature of real conflict between loci; that missing species imitate a
  diversification slowdown exactly; that comparing DEC with DEC+J by likelihood is biased (Ree & Sanmartín 2018);
  and that μ and ε from a tree of living species alone are not to be trusted.
- **It installs itself.** A collapsible interpretation guide is appended to each of the twelve block panels, and the
  two functions every block already uses to print a table or a row of tiles — `statTiles` and `buildTable` — are
  wrapped once, so a label the app prints anywhere gets its badge without any panel having to remember to ask.
  The wrapper swallows its own errors: the help must never break what it is decorating.
- `tests/help.js`: **41 tests**, total **1038**, all passing. They check that every entry is complete and bilingual,
  and — the one a decision rule cannot get wrong — that **every scale is a real partition**: ordered, with no gaps
  and no overlaps, so no value can fall between two bands and come back without a verdict.
- **The missing entries, found by counting rather than guessing.** Every label the twelve blocks print in a tile or a
  table header was extracted from the source and compared against the registry. **25 more entries**, total **71**,
  and the labels the help recognises went from 33 to 253: trimming of ambiguous blocks; constant, variable and
  informative sites; G+C content; the ts/tv ratio; the hierarchical likelihood ratio test; corrected distances;
  equally parsimonious trees; tree search by NNI, SPR and TBR; partition rate multipliers; the marginal likelihood
  as distinct from the Bayes factor; the substitution rate and its variation among branches; crown against stem age;
  the MCCR test; doubling time; independent contrasts; stochastic mapping; d and e of DEC as distinct from j;
  biogeographic stochastic mapping; splits; quartets; the discordance factors gDF and sDF; the tanglegram and its
  crossings; split networks; and the reproducible report.
- What is deliberately left without an entry is structural — "Taxon", "Method", "file", "bytes", "Time" — and the
  bare **p, t, z, D and Δ**, which in this app belong to several different tests: sending all of them to one entry
  would be worse than leaving them alone.

## [1.3.0] — 2026-09-25 · The look of every figure

**`js/figstyle.js`**: colours, type, stroke width and background for **every figure in the twelve blocks**, changed
from one place.

- It does not edit figures one at a time. Every plot in this program already draws its colours through CSS custom
  properties and its text through two or three classes — never a literal colour, so that the figures would follow
  the light and the dark theme. That has a second consequence worth taking: **overriding those properties restyles
  the whole app at once**, and the studio writes a single rule, `svg { --c1: …; --fig-font: … }`. Not one plotting
  module knows it exists. The rule is scoped to `svg` on purpose — the coloured chips of the interface use the same
  tokens, and a figure palette has no business recolouring the menus.
- **Nine palettes**, each explained: the program's own, **Okabe–Ito** (colour-blind safe), viridis, earth, cool,
  warm, pastels, high contrast, and **greyscale — for checking whether your figure survives a journal that charges
  for colour**. Plus the ten colours editable one by one, and nine more for text, axes and the support colours.
- **Seven families of type**, all of them present on a normal computer: nothing is fetched from the network, because
  a figure that falls back to something else on the reviewer's machine is a figure you did not design. Type size and
  stroke width scale from 0.7× to 1.6× and 0.5× to 3×. Four backgrounds, including transparent and white.
- **The export carries what you see.** An SVG loaded into an `<img>` has no stylesheet, so the two functions of
  `js/exportfmt.js` that resolve the theme are wrapped: the chosen colours win over the stylesheet, the chosen
  family replaces the hardcoded `system-ui`, and a chosen background is written in as a rectangle beneath
  everything. Checked by exporting a tree and finding the Okabe blue and the white ground in the PNG's own pixels.
- The tree's **group colours now follow the palette**, so the k groups of Block 12 are recoloured with everything
  else instead of being the one thing that stays fixed.
- **Fixed: `.art-ink` was used by the chronogram for its tip names and was never defined in the stylesheet.** Those
  names fell back to black, which is invisible in the dark theme — the same fault as the `--ink` token found while
  building Block 9, in the one place that had been missed. It is defined now, and written into exported figures too.
- **Fixed: the studio patched the export on `DOMContentLoaded`, and `js/figstyle.js` loaded before
  `js/exportfmt.js`** — so there was nothing to patch when it looked, and anything exporting before the document was
  ready resolved its colours against the stylesheet instead of against the studio. The test page does exactly that,
  which is how it surfaced. The script now loads after the one it wraps.
- `tests/figstyle.js`: **34 tests**, total **1072**, all passing. Most of them check the one thing this design can
  get wrong: that the figure on screen and the figure in the exported file are the same picture.

## [1.4.0] — 2026-09-25 · The calibration assistant

**`js/calibhelp.js`**: help for the one step of a dating study that goes wrong most often, without taking the
decision away from the user.

- **The program still ships no fossil database, and never will.** Assigning a fossil to a node is a taxonomic
  judgement about that fossil's characters and about which clade it is sister to — the part reviewers argue about.
  A lookup table would invite people to publish a date nobody examined. What the assistant does is everything
  *around* that judgement.
- **It says what the prior actually asserts, in years.** "Offset lognormal, M = 1.6, S = 0.7, offset 20" means
  nothing to most people; "the node's age is between 21.3 and 39.7 Ma with a median of 25.0, and cannot be younger
  than 20" means everything. Every shape — fixed, uniform, normal, offset lognormal, offset exponential — is reduced
  to the same 2.5 %, median and 97.5 %, recomputed on every keystroke, before anything is added.
- **It suggests the shape that matches where the age came from.** A fossil is a minimum, so an offset lognormal; a
  geological event is a window, so a uniform; an age borrowed from another study carries that study's uncertainty,
  so a normal, and one worth widening rather than narrowing.
- **It catches the mistakes that ruin a dating study**: a fossil entered as a fixed age — which asserts the lineage
  arose exactly when that individual died and squeezes the whole tree towards the present; a fossil given a
  symmetric normal, which allows the node to be younger than the fossil; an offset of zero, which imposes no minimum
  at all; a normal whose lower tail reaches into the future; a prior so wide it barely informs or so narrow it
  feigns certainty; **the crown-versus-stem confusion**, with the reminder that the selected taxa calibrate the
  crown and that a sister-group fossil calibrates the older stem; and a chronogram hanging from a single date.
- Found while testing: the "this prior barely informs" rule measured width as the 95 % range **relative to the
  median**, which for a uniform interval cannot exceed 1.9 by construction — so the warning could never fire on the
  very shape people reach for when they are being vague. It now uses the **ratio between the two extremes**, which
  depends on neither the scale nor the shape.
- `tests/calibhelp.js`: **44 tests**, total **1116**, all passing. The quantiles are checked against the closed form
  of each distribution, because a user who believes the median of that lognormal is 21.6 instead of 24.95 will
  publish a date that is wrong by millions of years.
