# PhylogenyPro — generation of the example data sets (Block 2)
#
# The sequences are SIMULATED on a known tree under a known model, and the
# indels are inserted by lineage, so every example carries its own truth:
#   · the true tree                      -> *_true_tree.nwk
#   · the true alignment                 -> *_true_alignment.fasta
#   · the unaligned sequences the user reads into the app -> *.fasta
# That lets the aligner of Block 2 be scored against the alignment that really
# generated the data, which is how alignment benchmarks are built.
#
# Species names are real names of Mexican plants; the sequences are not.
#
# Run:  "C:\Program Files\R\R-4.4.2\bin\Rscript.exe" generar_ejemplos.R

suppressMessages({library(ape); library(phangorn)})
set.seed(20260923)

taxa <- c("Bursera_simaruba","Bursera_fagaroides","Bursera_morelensis","Bursera_copallifera",
          "Bursera_bipinnata","Bursera_glabrifolia","Bursera_grandifolia","Bursera_linanoe",
          "Bursera_schlechtendalii","Bursera_aptera","Bursera_lancifolia","Bursera_cuneata",
          "Bursera_palmeri","Bursera_arborea","Bursera_excelsa","Bursera_hindsiana",
          "Commiphora_africana","Boswellia_sacra")            # last two = outgroup

n <- length(taxa)
# the ingroup is drawn at random, the outgroup is glued on, so that the outgroup
# is monophyletic by construction
ing <- rtree(16, tip.label = taxa[1:16], br = function(k) rexp(k, 60))
ing_txt <- sub(";$", "", write.tree(ing))
tree <- read.tree(text = paste0("(", ing_txt, ":0.018,(Commiphora_africana:0.030,Boswellia_sacra:0.028):0.020);"))
tree$edge.length <- pmax(tree$edge.length, 0.0008)

## ---------------------------------------------------------------- indels by lineage
# A deletion turns a run of columns into gaps in every descendant of one node;
# an insertion adds columns that only that clade has.
clade_tips <- function(tr, node) if (node <= Ntip(tr)) tr$tip.label[node] else extract.clade(tr, node)$tip.label

add_indels <- function(aln, tr, n_del, n_ins, maxlen = 12) {
  nodes <- c(seq_len(Ntip(tr)), (Ntip(tr) + 2):(Ntip(tr) + Nnode(tr)))
  for (i in seq_len(n_del)) {
    who <- clade_tips(tr, sample(nodes, 1))
    if (length(who) >= nrow(aln) - 1) next
    len <- 1 + rgeom(1, 0.25); len <- min(len, maxlen)
    pos <- sample(seq_len(ncol(aln) - len), 1)
    aln[who, pos:(pos + len - 1)] <- "-"
  }
  for (i in seq_len(n_ins)) {
    who <- clade_tips(tr, sample(nodes, 1))
    if (length(who) >= nrow(aln) - 1) next
    len <- 1 + rgeom(1, 0.3); len <- min(len, maxlen)
    pos <- sample(seq_len(ncol(aln) - 1), 1)
    block <- matrix("-", nrow(aln), len, dimnames = list(rownames(aln), NULL))
    block[who, ] <- sample(c("a","c","g","t"), length(who) * len, replace = TRUE)
    aln <- cbind(aln[, seq_len(pos), drop = FALSE], block, aln[, (pos + 1):ncol(aln), drop = FALSE])
  }
  aln
}

write_fasta <- function(mat, file, degap = FALSE, width = 60) {
  con <- file(file, "w", encoding = "UTF-8")
  for (i in seq_len(nrow(mat))) {
    s <- toupper(paste(mat[i, ], collapse = ""))
    if (degap) s <- gsub("-", "", s, fixed = TRUE)
    cat(">", rownames(mat)[i], "\n", sep = "", file = con)
    for (k in seq(1, nchar(s), by = width)) cat(substr(s, k, k + width - 1), "\n", sep = "", file = con)
  }
  close(con)
  invisible(nchar(paste(mat[1, ], collapse = "")))
}

sim_gene <- function(tr, len, rate, model, kappa, alpha, n_del, n_ins, tag) {
  t2 <- tr; t2$edge.length <- t2$edge.length * rate
  d <- simSeq(t2, l = len, type = "DNA",
              bf = c(0.30, 0.19, 0.21, 0.30), Q = c(1, kappa, 1, 1, kappa, 1),
              rate = 1, ancestral = FALSE)
  aln <- as.character(as.DNAbin(d))
  aln <- aln[tr$tip.label, ]
  aln <- add_indels(aln, tr, n_del, n_ins)
  cat(sprintf("%-10s %4d columns in the true alignment\n", tag, ncol(aln)))
  aln
}

## ---------------------------------------------------------------- 1. three genes, three genomes
# rbcL is a CODING gene, so it is simulated as a protein and back-translated with
# synonymous codons: no internal stop codons, a real reading frame, and third
# positions far more variable than first and second — which is what makes it
# useful for codon-aware alignment and for partitioning by codon position.
codon_table <- list(
  F = c("TTT","TTC"), L = c("TTA","TTG","CTT","CTC","CTA","CTG"), I = c("ATT","ATC","ATA"),
  M = "ATG", V = c("GTT","GTC","GTA","GTG"), S = c("TCT","TCC","TCA","TCG","AGT","AGC"),
  P = c("CCT","CCC","CCA","CCG"), T = c("ACT","ACC","ACA","ACG"), A = c("GCT","GCC","GCA","GCG"),
  Y = c("TAT","TAC"), H = c("CAT","CAC"), Q = c("CAA","CAG"), N = c("AAT","AAC"),
  K = c("AAA","AAG"), D = c("GAT","GAC"), E = c("GAA","GAG"), C = c("TGT","TGC"),
  W = "TGG", R = c("CGT","CGC","CGA","CGG","AGA","AGG"), G = c("GGT","GGC","GGA","GGG"))

sim_coding <- function(tr, aa_len, rate, tag) {
  t2 <- tr; t2$edge.length <- t2$edge.length * rate
  prot <- as.character(simSeq(t2, l = aa_len, type = "AA", model = "LG"))
  prot <- toupper(prot[tr$tip.label, ])
  nt <- t(apply(prot, 1, function(row) {
    unlist(lapply(row, function(a) {
      cods <- codon_table[[a]]
      if (is.null(cods)) cods <- "NNN"
      strsplit(cods[sample.int(length(cods), 1)], "")[[1]]
    }))
  }))
  rownames(nt) <- rownames(prot)
  nt <- tolower(nt)
  # one in-frame deletion of a whole codon, as happens in real coding genes
  who <- clade_tips(tr, sample((Ntip(tr) + 2):(Ntip(tr) + Nnode(tr)), 1))
  pos <- 3 * sample.int(aa_len - 4, 1) + 1
  nt[who, pos:(pos + 2)] <- "-"
  cat(sprintf("%-10s %4d columns in the true alignment (%d codons, coding)\n", tag, ncol(nt), aa_len))
  nt
}
rbcL <- sim_coding(tree, 207, 1.1, "rbcL")
# trnL-F: chloroplast spacer, faster and full of indels
trnLF <- sim_gene(tree, 430, 2.6, "HKY", 2.2, NULL, 9, 7, "trnL-F")
# ITS: nuclear ribosomal, the fastest
ITS   <- sim_gene(tree, 610, 4.2, "HKY", 2.0, NULL, 7, 5, "ITS")

dir.create("ejemplo_bursera", showWarnings = FALSE)
write_fasta(rbcL,  "ejemplo_bursera/rbcL_true_alignment.fasta")
write_fasta(trnLF, "ejemplo_bursera/trnLF_true_alignment.fasta")
write_fasta(ITS,   "ejemplo_bursera/ITS_true_alignment.fasta")
write_fasta(rbcL,  "ejemplo_bursera/rbcL.fasta",  degap = TRUE)
write_fasta(trnLF, "ejemplo_bursera/trnLF.fasta", degap = TRUE)
write_fasta(ITS,   "ejemplo_bursera/ITS.fasta",   degap = TRUE)
write.tree(tree, "ejemplo_bursera/true_tree.nwk")

# the same three genes concatenated, in the formats the app must read
conc <- cbind(rbcL, trnLF, ITS)
write_fasta(conc, "ejemplo_bursera/concatenado.fasta")
write.dna(as.DNAbin(conc), "ejemplo_bursera/concatenado_interleaved.phy", format = "interleaved", nbcol = -1, colsep = "")
write.dna(as.DNAbin(conc), "ejemplo_bursera/concatenado_sequential.phy", format = "sequential", nbcol = -1, colsep = "")
write.nexus.data(as.list(as.data.frame(t(conc), stringsAsFactors = FALSE)),
                 "ejemplo_bursera/concatenado.nex", format = "dna", interleaved = FALSE)
cat(sprintf("DNA, rbcL = 1-%d\n", ncol(rbcL)),
    sprintf("DNA, trnLF = %d-%d\n", ncol(rbcL) + 1, ncol(rbcL) + ncol(trnLF)),
    sprintf("DNA, ITS = %d-%d\n", ncol(rbcL) + ncol(trnLF) + 1, ncol(conc)),
    sep = "", file = "ejemplo_bursera/particiones.txt")

## ---------------------------------------------------------------- 2. a protein family
ptree <- rtree(22, tip.label = paste0("MYB", sprintf("%02d", 1:22), c(rep("_Zm", 8), rep("_At", 7), rep("_Os", 7))),
               br = function(k) rexp(k, 18))
paln <- as.character(simSeq(ptree, l = 185, type = "AA", model = "LG"))
paln <- add_indels(paln, ptree, 6, 4)
write_fasta(paln, "ejemplo_proteina_true_alignment.fasta")
write_fasta(paln, "ejemplo_proteina.fasta", degap = TRUE)
write.tree(ptree, "ejemplo_proteina_true_tree.nwk")
cat(sprintf("protein   %4d columns, %d sequences\n", ncol(paln), nrow(paln)))

## ---------------------------------------------------------------- 3. a morphological matrix
mtree <- rtree(25, tip.label = paste0("Taxon_", sprintf("%02d", 1:25)), br = function(k) rexp(k, 12))
nchar_m <- 42
mat <- sapply(seq_len(nchar_m), function(j) {
  k <- sample(2:3, 1, prob = c(0.72, 0.28))
  rTraitDisc(mtree, model = "ER", k = k, rate = runif(1, 0.6, 3), states = as.character(0:(k - 1)))
})
rownames(mat) <- mtree$tip.label
colnames(mat) <- paste0("car", sprintf("%02d", seq_len(nchar_m)))
mat[sample(length(mat), 18)] <- "?"          # a few missing entries, as in any real matrix
write.csv(cbind(taxon = rownames(mat), as.data.frame(mat)), "ejemplo_morfologia.csv", row.names = FALSE, quote = FALSE)
write.tree(mtree, "ejemplo_morfologia_true_tree.nwk")
cat(sprintf("morphology  %d characters, %d taxa\n", nchar_m, nrow(mat)))

cat("\nDone.\n")
