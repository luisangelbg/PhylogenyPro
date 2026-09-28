## The example trait table of Block 9. The species are the real ones of the
## app's Bursera example; the measurements are SIMULATED on that tree, so they
## carry genuine phylogenetic signal and nobody can mistake them for data.
suppressPackageStartupMessages({library(ape); library(phytools)})
nwk <- "((Bursera_arborea:0.01682243586,((((Bursera_cuneata:0.01990515465,Bursera_palmeri:0.01705012692):0.005192182685,((Bursera_simaruba:0.03419069402,Bursera_aptera:0.01049894283):0.01653436683,(Bursera_morelensis:0.004776271518,Bursera_glabrifolia:0.006036682047):0.01030173116):0.05872499753):0.00845428959,((Bursera_linanoe:0.02372140996,Bursera_grandifolia:0.01266004867):0.01400102644,(Bursera_hindsiana:0.01949399542,(Bursera_bipinnata:0.05386133442,(Bursera_schlechtendalii:0.01147089096,Bursera_fagaroides:0.004736136458):0.01180161975):0.02184093824):0.009034191833):0.0107184034):0.02175372587,(Bursera_copallifera:0.02265991548,(Bursera_excelsa:0.02026025054,Bursera_lancifolia:0.02962202253):0.0008):0.002379033105):0.005138137357):0.018,(Commiphora_africana:0.03,Boswellia_sacra:0.028):0.02);"
tr <- read.tree(text = nwk)
tr <- chronos(tr, lambda = 1, quiet = TRUE)   # a dated tree, as Block 7 leaves it
class(tr) <- "phylo"; tr$edge.length <- tr$edge.length * 40  # 40 Ma of crown age
set.seed(4021)

## leaf length: log-normal, Brownian on the log
lh <- fastBM(tr, a = log(9), sig2 = 0.007)
hoja <- round(exp(lh), 2)

## seed mass: correlated with the leaf through a shared Brownian component
sm <- 0.55 * lh + fastBM(tr, a = log(24) - 0.55 * log(9), sig2 = 0.004)
semilla <- round(exp(sm), 1)

## habit: three ordered states, from a slow Mk chain
Q <- matrix(c(-0.05, 0.05, 0, 0.03, -0.06, 0.03, 0, 0.05, -0.05), 3, 3, byrow = TRUE)
dimnames(Q) <- list(c("arbusto","arbolito","arbol"), c("arbusto","arbolito","arbol"))
h <- sim.history(tr, Q, anc = "arbolito", message = FALSE)
habito <- as.character(h$states[tr$tip.label])

d <- data.frame(taxon = tr$tip.label, hoja_cm = hoja[tr$tip.label],
                semilla_mg = semilla[tr$tip.label], habito = habito,
                stringsAsFactors = FALSE)
write.csv(d, "data/rasgos_bursera.csv", row.names = FALSE, quote = FALSE)
cat(paste(apply(d, 1, paste, collapse = ","), collapse = "\n"), "\n")
cat("\nestados:", table(habito), "\n")
cat("K de la hoja:", sprintf("%.4f", phylosig(tr, setNames(lh, tr$tip.label), method = "K")), "\n")
