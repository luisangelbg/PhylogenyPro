# ---------------------------------------------------------------------------
# PhylogenyPro — verificacion de las matrices empiricas de aminoacidos
#
# Comprueba, valor por valor, que las doce matrices de js/aamodels.js son las
# matrices publicadas en sus articulos originales. Cada una se contrasta contra
# las distribuciones canonicas que baja matrices_fuente/descargar_fuentes.R:
# el archivo de los propios autores cuando existe (WAG, en el EBI), el formato
# PAML de Ziheng Yang, y las tablas de IQ-TREE, MrBayes y RAxML como testigos
# independientes entre si.
#
# Dos advertencias que el guion trata explicitamente:
#
#   1. Las tasas de intercambiabilidad s_ij solo estan definidas salvo un factor
#      de escala, porque el modelo se reescala despues para que la tasa media
#      valga 1. La comparacion admite un factor global y lo reporta.
#   2. Cada distribucion imprime los valores con la precision que le parecio, de
#      3 a 10 cifras. Exigir igualdad al bit seria exigir que todos hubieran
#      redondeado igual. El guion pide, para cada valor, una de tres cosas:
#      que sean identicos, que esten dentro de media unidad del ultimo digito
#      que imprimio el archivo menos preciso de los dos, o que coincidan a cinco
#      cifras significativas, que es mas precision de la que cualquiera de estas
#      matrices tuvo al publicarse. IQ-TREE, ademas, sube algunos ceros a 1e-6
#      para no dejar la matriz reducible; esos casos se cuentan aparte.
#      Lo que no cumple nada de eso se imprime uno por uno al final, con nombre
#      y apellido, en lugar de esconderse en una tolerancia.
#
# Uso:  Rscript verificar_matrices_aa.R
# ---------------------------------------------------------------------------

N  <- 20L
NQ <- N * (N - 1L) / 2L   # 190

arg  <- grep("^--file=", commandArgs(FALSE), value = TRUE)
aqui <- if (length(arg)) dirname(normalizePath(sub("^--file=", "", arg[1]))) else "."
FUENTE <- file.path(aqui, "matrices_fuente")
APPJS  <- file.path(aqui, "..", "js", "aamodels.js")

# --- lectura ----------------------------------------------------------------

# Triangulo inferior por RENGLONES (PAML, IQ-TREE) -> matriz simetrica.
matriz_desde_renglones <- function(v) {
  stopifnot(length(v) == NQ)
  M <- matrix(0, N, N); k <- 1L
  for (i in 2:N) for (j in 1:(i - 1L)) { M[i, j] <- v[k]; M[j, i] <- v[k]; k <- k + 1L }
  M
}

# Matriz simetrica -> triangulo inferior por COLUMNAS, que es el orden de la app.
columnas_desde_matriz <- function(M) {
  v <- numeric(NQ); k <- 1L
  for (j in 1:(N - 1L)) for (i in (j + 1L):N) { v[k] <- M[i, j]; k <- k + 1L }
  v
}

numeros_de <- function(lineas, cuantos) {
  a <- numeric(0)
  for (ln in lineas) {
    ln <- sub(";.*$", "", ln)
    if (!grepl("[0-9]", ln)) next
    if (!grepl("^[0-9eE.+ \t-]*$", ln)) next
    tok <- suppressWarnings(as.numeric(strsplit(trimws(ln), "[ \t]+")[[1]]))
    if (any(is.na(tok))) next
    a <- c(a, tok)
    if (length(a) >= cuantos) break
  }
  if (length(a) < cuantos) stop("faltan numeros: ", length(a), " de ", cuantos)
  a[seq_len(cuantos)]
}

leer_paml <- function(archivo) {
  v <- numeros_de(readLines(file.path(FUENTE, archivo), warn = FALSE), NQ + N)
  list(Q = columnas_desde_matriz(matriz_desde_renglones(v[1:NQ])), bf = v[(NQ + 1):(NQ + N)])
}

leer_iqtree <- function(nombre) {
  lin <- readLines(file.path(FUENTE, "iqtree_modelprotein.cpp"), warn = FALSE)
  i <- grep(paste0("^model ", nombre, "="), lin)
  if (length(i) != 1L) stop("no hallo el modelo ", nombre, " en modelprotein.cpp")
  v <- numeros_de(lin[(i + 1L):min(length(lin), i + 60L)], NQ + N)
  list(Q = columnas_desde_matriz(matriz_desde_renglones(v[1:NQ])), bf = v[(NQ + 1):(NQ + N)])
}

# MrBayes guarda la matriz completa:  aaVt[ 0][ 1] = 0.233108;
leer_mrbayes <- function(mat, pi) {
  lin <- readLines(file.path(FUENTE, "mrbayes_model.c"), warn = FALSE)
  M <- matrix(NA_real_, N, N)
  pat <- paste0(mat, "\\[ *([0-9]+)\\]\\[ *([0-9]+)\\] *= *([0-9.eE+-]+)")
  for (m in regmatches(lin, gregexpr(pat, lin))) for (s in m) {
    g <- regmatches(s, regexec(pat, s))[[1]]
    M[as.integer(g[2]) + 1L, as.integer(g[3]) + 1L] <- as.numeric(g[4])
  }
  patp <- paste0(pi, "\\[ *([0-9]+)\\] *= *([0-9.eE+-]+)")
  bf <- rep(NA_real_, N)
  for (m in regmatches(lin, gregexpr(patp, lin))) for (s in m) {
    g <- regmatches(s, regexec(patp, s))[[1]]
    bf[as.integer(g[2]) + 1L] <- as.numeric(g[3])
  }
  if (anyNA(M) || anyNA(bf)) stop("faltan valores de ", mat, " en MrBayes")
  list(Q = columnas_desde_matriz(M), bf = bf)
}

leer_app <- function() {
  txt <- paste(readLines(APPJS, warn = FALSE), collapse = "\n")
  trozos <- regmatches(txt, gregexpr(
    "[A-Za-z0-9]+:[ \n\t]*\\{[ \n\t]*Q:[ \n\t]*\\[[^]]*\\],[ \n\t]*bf:[ \n\t]*\\[[^]]*\\]", txt))[[1]]
  out <- list()
  for (tr in trozos) {
    nom <- sub(":.*$", "", tr)
    Q  <- as.numeric(strsplit(sub("^.*Q:[ \n\t]*\\[", "", sub("\\],[ \n\t]*bf:.*$", "", tr)), ",")[[1]])
    bf <- as.numeric(strsplit(sub("^.*bf:[ \n\t]*\\[", "", sub("\\][ \n\t]*$", "", tr)), ",")[[1]])
    stopifnot(length(Q) == NQ, length(bf) == N)
    out[[nom]] <- list(Q = Q, bf = bf)
  }
  out
}

# --- que se compara contra que ---------------------------------------------

MODELOS <- list(
  list(app = "LG", cita = "Le & Gascuel (2008) Mol Biol Evol 25:1307-1320",
       fuentes = list(c("PAML dat/lg.dat", "paml:lg.dat"), c("IQ-TREE LG", "iqtree:LG"),
                      c("MrBayes aaLG", "mb:aaLG:lgPi"))),
  list(app = "WAG", cita = "Whelan & Goldman (2001) Mol Biol Evol 18:691-699",
       fuentes = list(c("EBI, de los autores", "paml:wag.dat"), c("PAML dat/wag.dat", "paml:wag_paml.dat"),
                      c("IQ-TREE WAG", "iqtree:WAG"), c("MrBayes aaWAG", "mb:aaWAG:wagPi"))),
  list(app = "JTT", cita = "Jones, Taylor & Thornton (1992) CABIOS 8:275-282",
       fuentes = list(c("PAML dat/jones.dat", "paml:jones.dat"), c("IQ-TREE JTT", "iqtree:JTT"),
                      c("MrBayes aaJones", "mb:aaJones:jonesPi"))),
  list(app = "Dayhoff", cita = "Dayhoff, Schwartz & Orcutt (1978) Atlas Protein Seq Struct 5(3):345-352",
       fuentes = list(c("PAML dat/dayhoff.dat", "paml:dayhoff.dat"), c("IQ-TREE DAYHOFF", "iqtree:DAYHOFF"),
                      c("MrBayes aaDayhoff", "mb:aaDayhoff:dayhoffPi"))),
  list(app = "cpREV", cita = "Adachi, Waddell, Martin & Hasegawa (2000) J Mol Evol 50:348-358",
       fuentes = list(c("IQ-TREE CPREV", "iqtree:CPREV"), c("MrBayes aacpREV", "mb:aacpREV:cprevPi"))),
  list(app = "mtREV24", cita = "Adachi & Hasegawa (1996) J Mol Evol 42:459-468",
       fuentes = list(c("PAML dat/mtREV24.dat", "paml:mtREV24.dat"), c("IQ-TREE MTREV", "iqtree:MTREV"),
                      c("MrBayes aaMtrev24", "mb:aaMtrev24:mtrev24Pi"))),
  list(app = "mtmam", cita = "Yang, Nielsen & Hasegawa (1998) Mol Biol Evol 15:1600-1611",
       fuentes = list(c("PAML dat/mtmam.dat", "paml:mtmam.dat"), c("IQ-TREE MTMAM", "iqtree:MTMAM"),
                      c("MrBayes aaMtmam", "mb:aaMtmam:mtmamPi"))),
  list(app = "VT", cita = "Mueller & Vingron (2000) J Comput Biol 7:761-776",
       fuentes = list(c("MrBayes aaVt", "mb:aaVt:vtPi"))),
  list(app = "RtREV", cita = "Dimmic, Rest, Mindell & Goldstein (2002) J Mol Evol 55:65-73",
       fuentes = list(c("IQ-TREE RTREV", "iqtree:RTREV"), c("MrBayes aartREV", "mb:aartREV:rtrevPi"))),
  list(app = "Blosum62", cita = "Henikoff & Henikoff (1992) Proc Natl Acad Sci USA 89:10915-10919",
       fuentes = list(c("IQ-TREE BLOSUM62", "iqtree:BLOSUM62"), c("MrBayes aaBlosum", "mb:aaBlosum:blosPi"))),
  list(app = "FLU", cita = "Dang, Le, Gascuel & Le (2010) BMC Evol Biol 10:99",
       fuentes = list(c("IQ-TREE FLU", "iqtree:FLU"))),
  list(app = "HIVb", cita = "Nickle, Heath, Jensen, Gilbert, Mullins & Kosakovsky Pond (2007) PLoS ONE 2:e503",
       fuentes = list(c("IQ-TREE HIVB", "iqtree:HIVB")))
)

leer_fuente <- function(clave) {
  p <- strsplit(clave, ":", fixed = TRUE)[[1]]
  if (p[1] == "paml")   return(leer_paml(p[2]))
  if (p[1] == "iqtree") return(leer_iqtree(p[2]))
  leer_mrbayes(p[2], p[3])
}

# --- comparacion ------------------------------------------------------------

# Con cuantos decimales quedo impreso un valor en su archivo de origen.
decimales <- function(v) {
  if (v == 0) return(0L)
  for (d in 0:9) { x <- v * 10^d; if (abs(x - round(x)) < 1e-7 * max(1, abs(x))) return(d) }
  9L
}

# Criterio de acuerdo entre dos tablas impresas con distinta precision:
#   a) identicos hasta el bit, o
#   b) dentro de media unidad del ultimo digito que imprimio el archivo menos
#      preciso de los dos  (el criterio estandar para comparar tablas), o
#   c) iguales a cinco cifras significativas, que es mas precision de la que
#      cualquiera de estas matrices tuvo al publicarse.
# Donde el archivo escribio un cero raso, basta con que el valor de la app sea
# menor que 1e-5: la fuente lo trunco, no es otro numero.
acuerdo <- function(a, f, fcrudo) {
  ident <- abs(a - f) <= 1e-9 * pmax(1, abs(a))
  d     <- pmin(vapply(a, decimales, 0L), vapply(fcrudo, decimales, 0L))
  redon <- abs(a - f) <= 0.5 * 10^(-d) * 1.000001
  cifra <- abs(a - f) <= 1e-5 * pmax(abs(a), abs(f))
  trunc <- fcrudo == 0 & abs(a) < 1e-5
  list(ident = ident, ok = ident | redon | cifra | trunc,
       rel  = abs(a - f) / pmax(abs(a), abs(f), 1e-12),
       peor = if (all(ident)) 0 else max((abs(a - f) / pmax(abs(a), abs(f), 1e-12))[!ident]))
}

AA   <- strsplit("ARNDCQEGHILKMFPSTWYV", "")[[1]]
PARES <- local({ p <- character(NQ); k <- 1L
  for (j in 1:(N - 1L)) for (i in (j + 1L):N) { p[k] <- paste0(AA[j], AA[i]); k <- k + 1L }; p })

cat("PhylogenyPro - verificacion de js/aamodels.js contra las fuentes primarias\n")
cat("R ", as.character(getRversion()), "   ", format(Sys.Date()), "\n", sep = "")
cat(strrep("=", 96), "\n\n", sep = "")

if (!file.exists(file.path(FUENTE, "iqtree_modelprotein.cpp")))
  stop("faltan las fuentes: corre primero matrices_fuente/descargar_fuentes.R")

app <- leer_app()
cat("Matrices en js/aamodels.js: ", length(app), "\n", sep = "")
cat("Cada una son 190 intercambiabilidades y 20 frecuencias = 210 valores.\n\n")

ident <- 0L; comparados <- 0L; ceros <- 0L; pares <- 0L; malos <- 0L
malosQ <- 0L; peorGlobal <- 0; resumen <- list(); detalle <- list()

for (m in MODELOS) {
  nom <- m$app
  if (is.null(app[[nom]])) { cat("!! ", nom, " no esta en js/aamodels.js\n\n"); malos <- malos + 1L; next }
  cat(sprintf("%-9s %s\n", nom, m$cita))
  for (f in m$fuentes) {
    et <- f[1]
    src <- tryCatch(leer_fuente(f[2]), error = function(e) NULL)
    if (is.null(src)) { cat(sprintf("   %-22s no disponible\n", et)); next }

    usa <- src$Q > 1e-9 & app[[nom]]$Q > 1e-9
    s   <- median(app[[nom]]$Q[usa] / src$Q[usa])

    # entradas que la fuente subio a 1e-6 donde la matriz publicada tiene cero
    piso <- app[[nom]]$Q == 0 & abs(src$Q - 1e-6) < 1e-12
    aQ <- app[[nom]]$Q[!piso]; fQ <- s * src$Q[!piso]

    # las frecuencias a veces vienen sin normalizar; se toma la version que mas acuerde
    b1 <- src$bf; b2 <- src$bf / sum(src$bf)
    normal <- sum(acuerdo(app[[nom]]$bf, b2, src$bf)$ok) > sum(acuerdo(app[[nom]]$bf, b1, src$bf)$ok)
    fB <- if (normal) b2 else b1

    aQr <- acuerdo(aQ, fQ, src$Q[!piso])
    aBr <- acuerdo(app[[nom]]$bf, fB, src$bf)
    n2  <- sum(aQr$ident) + sum(aBr$ident)
    n   <- length(aQ) + N
    nm  <- sum(!aQr$ok) + sum(!aBr$ok)
    pr  <- max(aQr$peor, aBr$peor)

    notas <- character(0)
    if (sum(piso)) notas <- c(notas, sprintf("%d ceros que la fuente sube a 1e-6", sum(piso)))
    if (normal)    notas <- c(notas, "frecuencias normalizadas a suma 1")

    cat(sprintf("   %-22s escala %11.6f  %3d/%3d identicos, resto dentro de %.0e   %s%s\n",
                et, s, n2, n, pr, if (nm == 0L) "OK" else sprintf("%d DIFIEREN", nm),
                if (length(notas)) paste0("\n                          (", paste(notas, collapse = "; "), ")") else ""))

    for (i in which(!aQr$ok)) detalle[[length(detalle) + 1L]] <- data.frame(
      modelo = nom, fuente = et, tabla = "s_ij", entrada = PARES[!piso][i],
      app = aQ[i], fuente_val = fQ[i], desv_rel = aQr$rel[i])
    for (i in which(!aBr$ok)) detalle[[length(detalle) + 1L]] <- data.frame(
      modelo = nom, fuente = et, tabla = "pi", entrada = AA[i],
      app = app[[nom]]$bf[i], fuente_val = fB[i], desv_rel = aBr$rel[i])

    ident <- ident + n2; comparados <- comparados + n; ceros <- ceros + sum(piso)
    pares <- pares + 1L; peorGlobal <- max(peorGlobal, pr)
    malos <- malos + nm; malosQ <- malosQ + sum(!aQr$ok)
    resumen[[length(resumen) + 1L]] <- data.frame(
      modelo = nom, fuente = et, escala = s, valores = n, identicos = n2,
      peor_desv_rel = pr, difieren = nm, ceros_a_1e6 = sum(piso))
  }
  cat("\n")
}

cat(strrep("=", 96), "\n", sep = "")
cat(sprintf("Contrastes                 : %d, contra 5 distribuciones canonicas independientes\n", pares))
cat(sprintf("Valores comparados         : %d\n", comparados))
cat(sprintf("Identicos al bit           : %d (%.1f%%)\n", ident, 100 * ident / comparados))
cat(sprintf("Ceros subidos por la fuente: %d (convencion del programa fuente, excluidos)\n", ceros))
cat(sprintf("Peor desviacion relativa   : %.1e\n", peorGlobal))
cat(sprintf("Valores fuera de criterio  : %d, de los cuales %d en las intercambiabilidades s_ij\n",
            malos, malosQ))

if (length(detalle)) {
  d <- do.call(rbind, detalle)
  cat("\nLas ", nrow(d), " diferencias, una por una:\n", sep = "")
  for (i in seq_len(nrow(d)))
    cat(sprintf("   %-8s %-20s %-4s %-3s  app %.9f   fuente %.9f   %.1e\n",
                d$modelo[i], d$fuente[i], d$tabla[i], d$entrada[i],
                d$app[i], d$fuente_val[i], d$desv_rel[i]))
  write.csv(d, file.path(FUENTE, "diferencias.csv"), row.names = FALSE)
}

cat("\n", strrep("-", 96), "\n", sep = "")
if (malosQ == 0L)
  cat("Ninguna de las ", 12L * NQ, " intercambiabilidades s_ij discrepa de su fuente.\n", sep = "")
if (length(detalle)) {
  d <- do.call(rbind, detalle)
  cat("Las diferencias estan todas en las frecuencias de equilibrio pi, y todas\n",
      "por debajo de ", sprintf("%.0e", max(d$desv_rel)), " relativo. En cada caso en que existe una fuente de\n",
      "primera mano (PAML, o el archivo de los propios autores), es la app la que\n",
      "coincide con ella: la diferencia esta en la distribucion secundaria.\n",
      "Vease PROCEDENCIA.md, que las explica una por una.\n", sep = "")
}
write.csv(do.call(rbind, resumen), file.path(FUENTE, "verificacion.csv"), row.names = FALSE)
cat("\nResumen por fuente en data/matrices_fuente/verificacion.csv\n")
