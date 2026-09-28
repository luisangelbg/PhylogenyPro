## Four "gene trees" for the Bursera example of Block 11, SIMULATED by moving
## the app's own example tree a little: each is the species tree with one or two
## nearest-neighbour interchanges and branch lengths jittered. They are what a
## real set of gene trees looks like — mostly agreeing, disagreeing in different
## places — without pretending to be anybody's data.
suppressMessages({library(ape); library(phangorn)})
set.seed(2718)
nwk <- "((Bursera_arborea:0.01682243586,((((Bursera_cuneata:0.01990515465,Bursera_palmeri:0.01705012692):0.005192182685,((Bursera_simaruba:0.03419069402,Bursera_aptera:0.01049894283):0.01653436683,(Bursera_morelensis:0.004776271518,Bursera_glabrifolia:0.006036682047):0.01030173116):0.05872499753):0.00845428959,((Bursera_linanoe:0.02372140996,Bursera_grandifolia:0.01266004867):0.01400102644,(Bursera_hindsiana:0.01949399542,(Bursera_bipinnata:0.05386133442,(Bursera_schlechtendalii:0.01147089096,Bursera_fagaroides:0.004736136458):0.01180161975):0.02184093824):0.009034191833):0.0107184034):0.02175372587,(Bursera_copallifera:0.02265991548,(Bursera_excelsa:0.02026025054,Bursera_lancifolia:0.02962202253):0.0008):0.002379033105):0.005138137357):0.018,(Commiphora_africana:0.03,Boswellia_sacra:0.028):0.02);"
tr <- read.tree(text = nwk)
out <- lapply(1:4, function(i) {
  g <- rNNI(tr, moves = sample(1:2, 1))
  g$edge.length <- pmax(1e-4, g$edge.length * rlnorm(length(g$edge.length), 0, 0.3))
  g
})
class(out) <- "multiPhylo"
for (i in 1:4) cat(write.tree(out[[i]]), "\n")
cat("\nRF contra el arbol base:", sapply(out, function(g) RF.dist(tr, g)), "\n")
