# Procedencia de PhylogenyPro

Este documento dice, elemento por elemento, qué partes de PhylogenyPro son obra
propia y cuáles no lo son, de dónde salen estas últimas y en qué situación
quedan. Se escribió para acompañar el registro de la obra ante el INDAUTOR y
para que cualquiera pueda comprobar lo que aquí se afirma sin tener que creerlo.

Última revisión: 24 de septiembre de 2026.

---

## 1. En una frase

PhylogenyPro no contiene código de nadie más. Contiene cuatro tablas de
constantes científicas publicadas —las doce matrices empíricas de sustitución de
aminoácidos, la matriz BLOSUM62, el código genético y la escala de tiempo
geológico—, todas citadas a su publicación original, todas verificables con los
guiones que vienen incluidos.

---

## 2. Lo que el programa no tiene

Conviene decirlo primero porque es lo que más suele ensuciar un programa así:

- **No usa ninguna biblioteca de terceros.** Ni jQuery, ni D3, ni Plotly, ni
  SheetJS, ni un motor de gráficas, ni nada. Las figuras son SVG escrito a mano,
  las tablas son HTML generado por `js/core.js`, el ZIP del Bloque 12 lo arma
  `js/exportfmt.js` byte por byte con el `CompressionStream` del propio
  navegador, y el álgebra lineal (Cholesky, exponencial de matrices,
  descomposición espectral) está escrita en `js/like.js`, `js/traits.js` y
  `js/biogeo.js`.
- **No carga nada de la red.** `index.html` no tiene una sola etiqueta que
  apunte fuera de la carpeta: el ícono es un SVG incrustado en el propio
  atributo, la hoja de estilo y los 47 archivos de JavaScript son locales. Por
  eso el programa abre con doble clic, sin servidor y sin conexión.
- **No hay código minificado, empaquetado ni pegado.** Los renglones más largos
  del proyecto son los arreglos de las matrices y de las secuencias de ejemplo.
- **No hay carpeta `vendor/`, `node_modules/`, ni dependencia de compilación.**

---

## 3. Las doce matrices empíricas de aminoácidos

**Dónde:** `js/aamodels.js`. **Qué son:** LG, WAG, JTT, Dayhoff, cpREV, mtREV24,
VT, RtREV, Blosum62, FLU, HIVb y mtMam: los parámetros de doce modelos de
sustitución de proteínas, cada uno con 190 tasas de intercambiabilidad y 20
frecuencias de equilibrio.

**De dónde salen:** de los artículos en que se publicaron. Cada matriz lleva su
cita al lado en el propio archivo:

| Modelo | Publicación original |
|---|---|
| LG | Le, S.Q. & Gascuel, O. (2008). *Mol. Biol. Evol.* 25: 1307–1320. |
| WAG | Whelan, S. & Goldman, N. (2001). *Mol. Biol. Evol.* 18: 691–699. |
| JTT | Jones, D.T., Taylor, W.R. & Thornton, J.M. (1992). *Comput. Appl. Biosci.* 8: 275–282. |
| Dayhoff | Dayhoff, M.O., Schwartz, R.M. & Orcutt, B.C. (1978). *Atlas of Protein Sequence and Structure* 5(3): 345–352. |
| cpREV | Adachi, J., Waddell, P.J., Martin, W. & Hasegawa, M. (2000). *J. Mol. Evol.* 50: 348–358. |
| mtREV24 | Adachi, J. & Hasegawa, M. (1996). *J. Mol. Evol.* 42: 459–468. |
| VT | Müller, T. & Vingron, M. (2000). *J. Comput. Biol.* 7: 761–776. |
| RtREV | Dimmic, M.W., Rest, J.S., Mindell, D.P. & Goldstein, R.A. (2002). *J. Mol. Evol.* 55: 65–73. |
| Blosum62 | Henikoff, S. & Henikoff, J.G. (1992). *Proc. Natl. Acad. Sci. USA* 89: 10915–10919. |
| FLU | Dang, C.C., Le, Q.S., Gascuel, O. & Le, V.S. (2010). *BMC Evol. Biol.* 10: 99. |
| HIVb | Nickle, D.C. *et al.* (2007). *PLoS ONE* 2: e503. |
| mtMam | Yang, Z., Nielsen, R. & Hasegawa, M. (1998). *Mol. Biol. Evol.* 15: 1600–1611. |

**Por qué se pueden incluir:** son las constantes numéricas de un modelo
científico publicado, no la expresión creativa de nadie. Se distribuyen, con los
mismos valores, en media docena de programas independientes entre sí, cada uno en
un orden y con una precisión distintos: PAML las escribe en formato de triángulo
inferior por renglones, IQ-TREE las repite en un archivo de C++, MrBayes las
guarda como matrices completas de 20×20, RAxML como asignaciones sueltas, y
PhylogenyPro las lee por columnas. Son los mismos números en cuatro
disposiciones distintas, que es exactamente lo que se espera de un dato y no de
una obra.

**Cómo se comprobó:** con `data/verificar_matrices_aa.R`, que contrasta cada
valor contra cinco distribuciones canónicas —entre ellas el archivo de los
propios autores de WAG, en el servidor del grupo de Nick Goldman en el EBI—.
Resultado de la corrida del 24 de septiembre de 2026:

```
Contrastes                 : 28, contra 5 distribuciones canonicas independientes
Valores comparados         : 5806
Identicos al bit           : 5242 (90.3%)
Valores fuera de criterio  : 8, de los cuales 0 en las intercambiabilidades s_ij
```

**Ninguna de las 2 280 intercambiabilidades discrepa.** Las ocho diferencias
están todas en frecuencias de equilibrio, y todas se explican:

| Modelo | Fuente | Letra | PhylogenyPro | La fuente | Explicación |
|---|---|---|---|---|---|
| Dayhoff | IQ-TREE | H | 0.033617966 | 0.033618970 | el archivo de PAML, que es la fuente de primera mano, dice 0.033618, y con ese coincide PhylogenyPro; IQ-TREE partió de una copia con el último dígito distinto |
| Dayhoff | IQ-TREE | K | 0.080481920 | 0.080480920 | igual que el anterior |
| mtREV24 | MrBayes | L | 0.169 | 0.168 | PAML y el artículo de Adachi y Hasegawa dicen 0.169; con 0.168 las veinte frecuencias de MrBayes suman 0.999 |
| FLU | IQ-TREE | C, Q, H, M, W | ≈1.3e-5 de diferencia relativa | | ambos renormalizaron a suma 1 partiendo de copias que difieren en el último dígito publicado; la diferencia está en el séptimo decimal de una frecuencia |

En los tres primeros casos existe fuente de primera mano y es PhylogenyPro el
que coincide con ella. En el cuarto no la hay, y la diferencia es menor que la
precisión con la que se publicó FLU.

**Un apunte sobre la versión anterior de este archivo.** Hasta el 23 de
septiembre de 2026 el encabezado de `js/aamodels.js` decía que los valores se
habían sacado del interior del paquete phangorn, que es lo que había hecho el
guion `data/extraer_matrices_aa.R`. Ese guion se eliminó y lo sustituyó
`data/verificar_matrices_aa.R`, que no extrae nada: verifica contra las fuentes
primarias. Los números no cambiaron —son los mismos que publicaron sus autores—,
pero ahora la procedencia que declara el archivo es la correcta, y se puede
comprobar.

---

## 4. La matriz BLOSUM62 del alineador

**Dónde:** `js/align.js`. Se usa para puntuar el alineamiento progresivo de
proteínas, que es un uso distinto del modelo de sustitución del mismo nombre.

**De dónde sale:** Henikoff, S. & Henikoff, J.G. (1992), *Proc. Natl. Acad. Sci.
USA* 89: 10915–10919, en la forma en que la distribuye el NCBI en
`ftp.ncbi.nlm.nih.gov/blast/matrices/BLOSUM62`.

**Situación:** obra del gobierno de los Estados Unidos, de dominio público. Se
reordenó al orden de aminoácidos que usa el programa.

---

## 5. El código genético

**Dónde:** `js/align.js`, tabla de traducción 1 del NCBI, usada por el
alineamiento consciente de codones y por la partición por posición.

**Situación:** el código genético es un hecho de la naturaleza y las tablas del
NCBI son obra del gobierno de los Estados Unidos, de dominio público.

---

## 6. La escala de tiempo geológico

**Dónde:** `js/geotime.js`: 66 unidades del Fanerozoico con su edad de base y de
techo en millones de años y su color.

**De dónde sale:** la Carta Cronoestratigráfica Internacional de la Comisión
Internacional de Estratigrafía, versión 2020/03, citada como Cohen, K.M., Finney,
S.C., Gibbard, P.L. & Fan, J.-X. (2013), *Episodes* 36: 199–204.

**Situación:** la ICS publica la carta para su uso libre a condición de citarla,
y se cita tanto en el encabezado del archivo como en la portada del programa. Las
edades son mediciones y los colores son los estándar de la cartografía geológica.

---

## 7. Los datos de ejemplo

**Dónde:** `js/examples.js` y `data/`.

Las secuencias **están simuladas**. `data/generar_ejemplos.R` las genera sobre un
árbol conocido y bajo un modelo conocido, e inserta los huecos por linaje, de
modo que cada ejemplo trae su propia verdad: el árbol verdadero y el alineamiento
verdadero. Eso es lo que permite calificar al alineador del Bloque 2 contra el
alineamiento que realmente generó los datos.

Los nombres de especie sí son reales —son plantas mexicanas, sobre todo del
género *Bursera*—, pero un nombre científico no es materia de derecho de autor y
las secuencias que llevan esos nombres no son las de ningún organismo.

No hay una sola secuencia descargada de GenBank ni de ninguna otra base en el
programa.

---

## 8. Los programas que se mencionan en los comentarios

En los comentarios del código y en los guiones de `validation/` aparecen los
nombres de ape, phangorn, geiger, phytools, nlme, BioGeoBEARS, IQ-TREE, RAxML,
MrBayes y PAML. Aparecen siempre en el mismo papel: el de **implementación de
referencia contra la cual se comprobó un resultado**. Frases como «validado
contra phangorn en validation/block3» o «con esta tolerancia el ajuste llegó a
−7165.4946 donde phangorn llega a −7165.4948» son constancias de una
comparación numérica, no señales de que se haya tomado algo de ellos.

Citar el programa contra el que uno comprobó su trabajo es práctica científica
normal y es lo contrario de ocultar una procedencia.

---

## 9. El material de verificación, que no es parte de la obra

Tres carpetas del proyecto no son el programa y conviene distinguirlas:

- **`validation/`** — los guiones de R que reproducen cada bloque en las
  implementaciones de referencia, con sus resultados. Incluye, en
  `validation/block10/`, el conjunto de datos de *Psychotria* de Hawái que
  distribuye BioGeoBEARS y que se remonta a Ree y Smith (2008): está ahí porque
  es el ejemplo con el que se comprueba que los seis modelos biogeográficos dan
  los mismos números que BioGeoBEARS, y por ningún otro motivo.
- **`data/matrices_fuente/`** — los archivos de PAML, IQ-TREE, MrBayes y RAxML
  que sirven de testigo a la verificación del punto 3. No se distribuyen: los
  baja `descargar_fuentes.R` cuando hacen falta. Véase el
  [LEEME](data/matrices_fuente/LEEME.md) de esa carpeta.
- **`research/`** — las notas de alcance y viabilidad previas al desarrollo.

**Recomendación para el depósito ante el INDAUTOR:** depositar el programa, esto
es `index.html`, `js/`, `css/`, `data/` sin `matrices_fuente/`, el `README.md` y
este documento, más el manual de usuario cuando esté. Las tres carpetas
anteriores son material de trabajo: no forman parte de la obra que se registra y
meterlas solo invita a preguntas sobre datos que no son del autor.

---

## 10. Sobre el carácter de la obra

PhylogenyPro es una **obra primigenia**. No es adaptación, traducción ni arreglo
de ningún programa anterior: los 47 archivos de JavaScript están escritos desde
cero, y los cuatro conjuntos de constantes de los puntos 3 a 6 son datos
científicos publicados, citados, y de libre uso por su naturaleza o por la
licencia de quien los publica.

Que un programa reproduzca los resultados de otro no lo convierte en obra
derivada de él, igual que dos implementaciones de la transformada rápida de
Fourier no son derivadas una de la otra por dar el mismo número.

Esta es la lectura del autor y el motivo por el que se documenta con este detalle;
la calificación del trámite corresponde al INDAUTOR.

---

## 11. La licencia

PhylogenyPro se distribuye bajo la **Licencia Pública General de GNU, versión 3 o
posterior**. El texto íntegro está en el archivo `LICENSE`, tal como lo publica la
Free Software Foundation, y el aviso que la aplica a este programa está en la
cabecera de `index.html`.

Dos cosas conviene dejar dichas, porque suelen confundirse:

- **La licencia es una decisión del autor, no una herencia.** El programa no
  contiene código de nadie más (punto 2), así que nada lo obligaba a ser
  copyleft. Se eligió la GPL-3.0 porque es la que corresponde a una herramienta
  de enseñanza e investigación que se quiere libre, y para que quien la
  modifique tenga que conservar esas mismas libertades.
- **Licenciar no es ceder la autoría.** Publicar bajo GPL no afecta al registro
  de la obra ni a los derechos morales, que en México son inalienables e
  imprescriptibles (Ley Federal del Derecho de Autor, artículos 18 a 21). El
  registro acredita quién es el autor; la licencia dice qué puede hacer el
  público con la obra. Son dos cosas distintas y compatibles.

---

## 12. Cómo comprobar todo esto

```
cd data/matrices_fuente
Rscript descargar_fuentes.R
cd ..
Rscript verificar_matrices_aa.R
```

Imprime los 28 contrastes, el factor de escala de cada uno, cuántos valores
resultaron idénticos y, una por una, cualquier diferencia que encuentre. Deja el
resumen en `data/matrices_fuente/verificacion.csv` y el detalle en
`data/matrices_fuente/diferencias.csv`.
