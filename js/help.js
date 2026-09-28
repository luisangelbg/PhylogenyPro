/* PhylogenyPro — molecular phylogenetics in the browser.
   Copyright (C) 2026  Luis Ángel Barrera-Guzmán

   This program is free software: you can redistribute it and/or modify it under
   the terms of the GNU General Public License as published by the Free Software
   Foundation, either version 3 of the License, or (at your option) any later
   version. It is distributed in the hope that it will be useful, but WITHOUT ANY
   WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
   A PARTICULAR PURPOSE. See the GNU General Public License, in the file LICENSE
   at the root of this program, or <https://www.gnu.org/licenses/>. */

/* PhylogenyPro — interpretation help.

   A number a reader cannot place is not a result. This is a registry of short
   explanations, one per concept or statistic the app prints, each with what it
   measures, how it is read, the scale that says whether a value is low or high,
   and the mistake most commonly made with it.

     Help.badge(key)             the circled question mark to put next to a label
     Help.panel(keys, title)     a collapsible guide listing several entries
     Help.markTable(box, map)    decorates the headers of a table already built
     Help.markTiles(box, map)    the same for the summary tiles
     Help.hydrate(root)          turns every <span data-help="key"> into a badge
     Help.band(key, value)       which band of the scale a computed value falls in
     Help.tag(key, value)        that band as a small coloured chip, ready to print

   The last two are the point. A scale printed in a manual is documentation; a
   scale applied to the number on screen is a decision rule, which is what the
   reader actually needed. Where a scale is a convention rather than a law, the
   entry says so — half of the thresholds in phylogenetics are habits, and a
   habit stated as a habit is honest.

   Everything is written in both languages with L2(), so switching the language
   needs no redraw, and every colour is a token of the stylesheet, so both
   themes follow on their own. Popovers are built the first time they are
   opened, never at load. */

const Help = {};

(function () {

  const two = p => (Array.isArray(p) ? L2(p[0], p[1]) : String(p || ''));
  const num = v => (Math.abs(v) >= 1000 || (v !== 0 && Math.abs(v) < 0.01) ? String(v) : String(+(+v).toFixed(4)));

  /* a band of a scale: from/to may be null for an open end */
  const S = (from, to, label, tone, txt) => ({ from, to, label, tone, txt });
  const HELP = {};
  const E = (key, def) => { HELP[key] = def; return key; };

  /* =====================================================================
     1 · the catalogue
     ===================================================================== */

  /* ---------- Block 1 · the ideas the whole app rests on ---------- */

  E('longbranch', {
    t: ['Atracción de ramas largas', 'Long-branch attraction'],
    what: ['Dos linajes que acumulan muchos cambios por su cuenta terminan pareciéndose <b>por azar</b>: con cuatro estados posibles, dos ramas largas coinciden en un sitio más veces de lo que su parentesco justifica. El método los junta en el árbol aunque no sean parientes.',
      'Two lineages that accumulate many changes on their own end up resembling each other <b>by chance</b>: with four possible states, two long branches match at a site more often than their relatedness warrants. The method then joins them on the tree although they are not relatives.'],
    read: ['Sospecha cuando dos ramas muy largas salgan juntas y el resto del árbol sea de ramas cortas, o cuando el grupo externo se pegue al taxón de evolución más rápida. La parsimonia es la más vulnerable; la verosimilitud con un modelo adecuado resiste mucho más.',
      'Suspect it when two very long branches come out together and the rest of the tree has short ones, or when the outgroup sticks to the fastest-evolving taxon. Parsimony is the most vulnerable; likelihood with an adequate model resists far better.'],
    care: ['No se arregla con más sitios: con datos infinitos la parsimonia converge <b>con más seguridad</b> al árbol equivocado. Se arregla partiendo las ramas largas con más taxones, quitando el taxón problemático a ver si el resto cambia, o usando un modelo mejor.',
      'More sites do not fix it: with infinite data parsimony converges <b>more confidently</b> on the wrong tree. It is fixed by breaking long branches with more taxa, by removing the offending taxon to see whether the rest changes, or by using a better model.'],
    ref: 'Felsenstein 1978; Bergsten 2005',
  });

  E('saturation', {
    t: ['Saturación', 'Saturation'],
    what: ['Un sitio que ha cambiado varias veces solo enseña el último cambio. Al aumentar la distancia, las diferencias observadas dejan de crecer y se estancan: los datos <b>ya no distinguen</b> parientes lejanos de muy lejanos.',
      'A site that has changed several times shows only the last change. As distance grows, the observed differences stop growing and plateau: the data <b>no longer distinguish</b> distant relatives from very distant ones.'],
    read: ['En la gráfica de diferencias contra distancia corregida, mientras la nube siga una recta el marcador informa. Cuando se dobla y se aplana, está saturado. Las transiciones se saturan antes que las transversiones, y la tercera posición del codón antes que la primera y la segunda.',
      'On the plot of differences against corrected distance, as long as the cloud follows a straight line the marker is informative. When it bends and flattens, it is saturated. Transitions saturate before transversions, and the third codon position before the first and the second.'],
    care: ['Saturado no es inservible: sigue resolviendo nodos recientes. Lo que no se puede es usarlo para los nodos profundos y creerle a su apoyo.',
      'Saturated is not useless: it still resolves recent nodes. What cannot be done is to use it for the deep nodes and believe its support.'],
    ref: 'Philippe et al. 2011',
  });

  /* ---------- Block 2 · data and alignment ---------- */

  E('informative', {
    t: ['Sitios informativos para parsimonia', 'Parsimony-informative sites'],
    what: ['Un sitio es informativo si tiene al menos <b>dos estados distintos, cada uno en al menos dos taxones</b>. Solo esos pueden preferir una topología sobre otra: un sitio con un único taxón distinto cuesta un paso en cualquier árbol.',
      'A site is informative if it has at least <b>two different states, each in at least two taxa</b>. Only those can prefer one topology over another: a site where a single taxon differs costs one step on any tree.'],
    read: ['Es la materia prima real del análisis. Cuenta los informativos, no las posiciones totales: mil posiciones con treinta informativas resuelven poco, y la cifra que suele citarse en los artículos es la primera.',
      'This is the real raw material of the analysis. Count the informative ones, not the total positions: a thousand positions with thirty informative ones resolve little, and the figure usually quoted in papers is the first.'],
    scaleTitle: ['Informativos por taxón', 'Informative sites per taxon'],
    scale: [
      S(null, 2, ['muy pocos', 'very few'], 'bad'),
      S(2, 10, ['justos', 'tight'], 'warn'),
      S(10, 40, ['razonables', 'reasonable'], 'ok'),
      S(40, null, ['holgados', 'ample'], 'good'),
    ],
    conv: true,
    care: ['Muchos sitios informativos no salvan un alineamiento mal hecho: columnas mal alineadas son sitios informativos <b>falsos</b>, y son los que más torcerán el árbol.',
      'Many informative sites do not rescue a bad alignment: misaligned columns are <b>false</b> informative sites, and they are the ones that will bend the tree most.'],
  });

  E('gaps', {
    t: ['Huecos y datos faltantes', 'Gaps and missing data'],
    what: ['La proporción de la matriz que no es un nucleótido o un aminoácido observado: inserciones, deleciones, extremos sin secuenciar y taxones ausentes de una partición.',
      'The proportion of the matrix that is not an observed nucleotide or amino acid: insertions, deletions, unsequenced ends and taxa absent from a partition.'],
    read: ['Importa menos el total que <b>cómo se reparte</b>. Un 40 % repartido parejo estorba poco; un 40 % concentrado en tres taxones los deja colgando de casi nada y su posición en el árbol será inestable.',
      'What matters is less the total than <b>how it is distributed</b>. 40 % spread evenly is little trouble; 40 % concentrated in three taxa leaves them hanging from almost nothing and their position on the tree will be unstable.'],
    scaleTitle: ['Huecos en la matriz', 'Gaps in the matrix'],
    scale: [
      S(null, 0.05, ['muy limpia', 'very clean'], 'good'),
      S(0.05, 0.2, ['normal', 'normal'], 'ok'),
      S(0.2, 0.5, ['mucha ausencia', 'a lot of absence'], 'warn'),
      S(0.5, null, ['más hueco que dato', 'more gap than data'], 'bad'),
    ],
    conv: true,
    care: ['La mayoría de los modelos tratan el hueco como «dato desconocido», no como un estado: una deleción compartida, que es evidencia filogenética real, se desperdicia. Por eso conviene mirarla aparte.',
      'Most models treat a gap as "unknown", not as a state: a shared deletion, which is real phylogenetic evidence, is thrown away. That is why it is worth looking at separately.'],
  });

  E('sp', {
    t: ['Puntaje SP del alineamiento', 'SP score of the alignment'],
    what: ['De todos los pares de residuos que el alineamiento verdadero pone en la misma columna, la fracción que el alineamiento estimado también junta. Solo se puede calcular cuando se conoce la verdad, es decir sobre datos simulados.',
      'Of all the pairs of residues the true alignment puts in the same column, the fraction the estimated alignment also puts together. It can only be computed when the truth is known, that is, on simulated data.'],
    read: ['Es la nota del alineador contra la respuesta correcta. Sirve para comparar opciones —por codón contra por nucleótido, con y sin refinamiento— en los ejemplos que traen su verdad.',
      'It is the aligner\'s mark against the right answer. It is useful for comparing options — by codon against by nucleotide, with and without refinement — on the examples that carry their own truth.'],
    scaleTitle: ['SP', 'SP'],
    scale: [
      S(null, 0.7, ['mal alineado', 'badly aligned'], 'bad'),
      S(0.7, 0.9, ['aceptable', 'acceptable'], 'warn'),
      S(0.9, 0.98, ['bueno', 'good'], 'ok'),
      S(0.98, null, ['casi perfecto', 'nearly perfect'], 'good'),
    ],
    conv: true,
    care: ['Con datos reales no hay verdad contra la cual medir. Un SP alto en el ejemplo no garantiza que tus datos queden igual de bien.',
      'With real data there is no truth to measure against. A high SP on the example does not guarantee your own data will come out as well.'],
  });

  E('chicomp', {
    t: ['Prueba χ² de composición', 'χ² test of composition'],
    what: ['Compara la composición de bases (o de aminoácidos) de cada secuencia contra la media de todas. Pregunta si alguna secuencia tiene una composición <b>distinta de las demás</b>.',
      'Compares the base (or amino-acid) composition of each sequence against the average of all. It asks whether any sequence has a composition <b>different from the rest</b>.'],
    read: ['Un valor p pequeño dice que la composición no es homogénea, y eso es un problema: casi todos los modelos suponen que sí lo es, y dos secuencias con composición parecida se atraen en el árbol aunque no sean parientes.',
      'A small p-value says the composition is not homogeneous, and that is a problem: nearly every model assumes it is, and two sequences of similar composition attract each other on the tree although they are not relatives.'],
    scaleTitle: ['Valor p', 'p-value'],
    scale: [
      S(null, 0.01, ['heterogénea, cuidado', 'heterogeneous, careful'], 'bad'),
      S(0.01, 0.05, ['dudosa', 'doubtful'], 'warn'),
      S(0.05, null, ['compatible con homogeneidad', 'compatible with homogeneity'], 'good'),
    ],
    care: ['La prueba compara cada secuencia con la media <b>de todas, incluida ella</b>, así que pierde potencia cuando hay pocos taxones. Y no rechazar no demuestra homogeneidad: demuestra que no había con qué rechazarla.',
      'The test compares each sequence against the mean <b>of all of them, itself included</b>, so it loses power when taxa are few. And failing to reject does not prove homogeneity: it proves there was nothing to reject it with.'],
    ref: 'Foster 2004',
  });

  /* ---------- Block 3 · substitution models ---------- */

  E('lnl', {
    t: ['Log-verosimilitud (lnL)', 'Log-likelihood (lnL)'],
    what: ['El logaritmo de la probabilidad de haber observado <b>estos datos</b> bajo <b>este árbol y este modelo</b>. Siempre negativo; menos negativo es mejor ajuste.',
      'The logarithm of the probability of having observed <b>these data</b> under <b>this tree and this model</b>. Always negative; less negative is a better fit.'],
    read: ['Su valor absoluto no significa nada por sí solo: depende del número de sitios y de taxones. Solo tiene sentido <b>comparar</b> lnL calculadas sobre exactamente los mismos datos.',
      'Its absolute value means nothing on its own: it depends on the number of sites and taxa. It only makes sense to <b>compare</b> lnL values computed on exactly the same data.'],
    formula: 'lnL = Σ_sitios ln P(sitio | árbol, modelo)',
    care: ['Un modelo con más parámetros <b>siempre</b> da mejor lnL, aunque los parámetros de más no sirvan para nada. Por eso se compara con AIC o BIC y no con la lnL a secas.',
      'A model with more parameters <b>always</b> gives a better lnL, even when the extra parameters are useless. That is why models are compared with AIC or BIC and not with lnL alone.'],
  });

  E('aic', {
    t: ['AIC, AICc y BIC', 'AIC, AICc and BIC'],
    what: ['Tres formas de premiar el ajuste y castigar los parámetros. AIC = −2·lnL + 2k; AICc le añade una corrección para muestras pequeñas; BIC castiga más fuerte, con ln(n) por parámetro. En los tres, <b>menor es mejor</b>.',
      'Three ways of rewarding fit and penalising parameters. AIC = −2·lnL + 2k; AICc adds a small-sample correction; BIC penalises harder, with ln(n) per parameter. In all three, <b>lower is better</b>.'],
    read: ['Lo que se lee es la <b>diferencia</b> con el mejor modelo, ΔAIC. Un modelo con Δ menor que 2 es prácticamente indistinguible del ganador; a partir de 10 está descartado.',
      'What is read is the <b>difference</b> from the best model, ΔAIC. A model with Δ below 2 is practically indistinguishable from the winner; from 10 up it is ruled out.'],
    scaleTitle: ['ΔAIC frente al mejor', 'ΔAIC from the best'],
    scale: [
      S(0, 2, ['igual de bueno', 'just as good'], 'good'),
      S(2, 4, ['casi igual', 'nearly as good'], 'ok'),
      S(4, 10, ['claramente peor', 'clearly worse'], 'warn'),
      S(10, null, ['descartado', 'ruled out'], 'bad'),
    ],
    conv: true,
    care: ['BIC tiende a elegir modelos más simples que AIC, y con miles de sitios los dos suelen elegir el más complejo disponible. Elegir el modelo no es el análisis: un GTR+I+G mal alineado sigue estando mal.',
      'BIC tends to pick simpler models than AIC, and with thousands of sites both usually pick the most complex one available. Choosing the model is not the analysis: a badly aligned GTR+I+G is still wrong.'],
    ref: 'Burnham & Anderson 2002; Posada & Buckley 2004',
  });

  E('akaike', {
    t: ['Peso de Akaike', 'Akaike weight'],
    what: ['Convierte las ΔAIC en probabilidades que suman uno: el peso de un modelo es la probabilidad de que <b>sea el mejor del conjunto comparado</b>.',
      'Turns the ΔAIC values into probabilities that sum to one: the weight of a model is the probability that it <b>is the best of the set compared</b>.'],
    read: ['Un ganador con peso 0.9 se lleva la decisión. Varios modelos repartiéndose el peso significa que los datos no distinguen entre ellos, y entonces conviene comprobar que el árbol no cambia según cuál se use.',
      'A winner with a weight of 0.9 settles it. Several models sharing the weight means the data cannot tell them apart, and then it is worth checking that the tree does not change with the choice.'],
    scaleTitle: ['Peso del mejor modelo', 'Weight of the best model'],
    scale: [
      S(null, 0.5, ['sin ganador claro', 'no clear winner'], 'warn'),
      S(0.5, 0.9, ['preferido', 'preferred'], 'ok'),
      S(0.9, null, ['ganador neto', 'clear winner'], 'good'),
    ],
    conv: true,
    care: ['«El mejor del conjunto» no es «bueno». Si todos los modelos candidatos son malos para estos datos, el peso 0.99 se lo lleva el menos malo.',
      '"The best of the set" is not "good". If every candidate model is bad for these data, the 0.99 weight goes to the least bad one.'],
  });

  E('alpha', {
    t: ['Forma de la gamma (α)', 'Gamma shape (α)'],
    what: ['No todos los sitios cambian a la misma velocidad. α describe cuán desigual es ese reparto: α pequeña significa unos pocos sitios muy rápidos y muchos casi congelados; α grande, velocidades parecidas.',
      'Not every site changes at the same speed. α describes how uneven that spread is: a small α means a few very fast sites and many almost frozen; a large α, similar speeds.'],
    read: ['Es el parámetro que más suele importar después de la topología. Por debajo de 1 la heterogeneidad es fuerte y no modelarla deforma las ramas largas.',
      'It is the parameter that usually matters most after the topology. Below 1 the heterogeneity is strong, and not modelling it distorts the long branches.'],
    scaleTitle: ['α', 'α'],
    scale: [
      S(null, 0.5, ['muy desigual', 'very uneven'], 'bad'),
      S(0.5, 1, ['desigual', 'uneven'], 'warn'),
      S(1, 3, ['moderada', 'moderate'], 'ok'),
      S(3, null, ['casi homogénea', 'almost homogeneous'], 'good'),
    ],
    conv: true,
    care: ['α y las longitudes de rama están correlacionadas: un optimizador que se detiene pronto deja las dos mal. Y +I con +G compiten por explicar lo mismo, así que sus valores por separado se interpretan con cuidado.',
      'α and the branch lengths are correlated: an optimiser that stops early leaves both wrong. And +I and +G compete to explain the same thing, so their separate values are read with care.'],
    ref: 'Yang 1994',
  });

  E('pinv', {
    t: ['Proporción de sitios invariables (I)', 'Proportion of invariable sites (I)'],
    what: ['La fracción de sitios que el modelo considera <b>incapaces</b> de cambiar, no simplemente sitios que no cambiaron por azar.',
      'The fraction of sites the model treats as <b>unable</b> to change, not merely sites that happened not to change.'],
    read: ['Valores altos son normales en genes con mucha restricción funcional. Lo que se compara es el ajuste con y sin +I, no el valor en sí.',
      'High values are normal in genes under strong functional constraint. What is compared is the fit with and without +I, not the value itself.'],
    care: ['+I y +G se solapan: juntos suelen dar valores extremos de los dos (I alto y α baja) que no significan lo que parecen. Muchos autores recomiendan usar +G solo.',
      '+I and +G overlap: together they often give extreme values of both (a high I and a low α) that do not mean what they appear to. Many authors recommend using +G alone.'],
    ref: 'Yang 2006, §4.3',
  });

  /* ---------- Block 4 · parsimony and distances ---------- */

  E('steps', {
    t: ['Pasos (longitud del árbol)', 'Steps (tree length)'],
    what: ['El número mínimo de cambios de estado que el árbol necesita para explicar la matriz. Menos pasos, árbol preferido por parsimonia.',
      'The smallest number of state changes the tree needs to explain the matrix. Fewer steps, the tree parsimony prefers.'],
    read: ['Solo se comparan longitudes sobre la <b>misma</b> matriz. Una diferencia de un paso entre dos topologías no es evidencia de nada: por eso existen el soporte de Bremer y el remuestreo.',
      'Lengths are only compared on the <b>same</b> matrix. A one-step difference between two topologies is evidence of nothing: that is what Bremer support and resampling are for.'],
  });

  E('ci', {
    t: ['Índice de consistencia (CI)', 'Consistency index (CI)'],
    what: ['El mínimo de cambios posible dividido entre los que el árbol necesita: CI = m/s. Vale 1 cuando ningún carácter cambia dos veces —ninguna homoplasia— y baja conforme aparece convergencia o reversión.',
      'The minimum possible number of changes divided by the number the tree needs: CI = m/s. It is 1 when no character changes twice — no homoplasy — and falls as convergence or reversal appear.'],
    read: ['Da idea de cuánta homoplasia hay, pero <b>baja sistemáticamente con el número de taxones</b>, así que no se comparan CI de estudios con tamaños distintos.',
      'It gives a sense of how much homoplasy there is, but it <b>falls systematically with the number of taxa</b>, so CI values from studies of different size are not comparable.'],
    scaleTitle: ['CI en matrices de tamaño medio', 'CI on mid-sized matrices'],
    scale: [
      S(null, 0.3, ['mucha homoplasia', 'much homoplasy'], 'bad'),
      S(0.3, 0.5, ['bastante', 'considerable'], 'warn'),
      S(0.5, 0.75, ['moderada', 'moderate'], 'ok'),
      S(0.75, null, ['poca', 'little'], 'good'),
    ],
    conv: true,
    care: ['Un CI alto con quince taxones y uno bajo con doscientos pueden significar lo mismo. El RI corrige parte de este problema.',
      'A high CI with fifteen taxa and a low one with two hundred can mean the same thing. RI corrects part of this problem.'],
    ref: 'Kluge & Farris 1969; Sanderson & Donoghue 1989',
  });

  E('ri', {
    t: ['Índice de retención (RI)', 'Retention index (RI)'],
    what: ['Cuánta de la similitud que <b>podría</b> ser sinapomorfía lo es de verdad: RI = (g−s)/(g−m), donde g es la longitud en el peor árbol posible.',
      'How much of the similarity that <b>could</b> be synapomorphy actually is: RI = (g−s)/(g−m), where g is the length on the worst possible tree.'],
    read: ['Menos sensible al número de taxones que el CI, así que es el que conviene reportar. Alto significa que los caracteres, aunque tengan homoplasia, siguen agrupando de forma congruente.',
      'Less sensitive to the number of taxa than CI, so it is the one worth reporting. High means that the characters, homoplasy and all, still group congruently.'],
    scaleTitle: ['RI', 'RI'],
    scale: [
      S(null, 0.4, ['poca estructura', 'little structure'], 'bad'),
      S(0.4, 0.6, ['algo', 'some'], 'warn'),
      S(0.6, 0.8, ['buena', 'good'], 'ok'),
      S(0.8, null, ['muy buena', 'very good'], 'good'),
    ],
    conv: true,
    ref: 'Farris 1989',
  });

  E('bootstrap', {
    t: ['Bootstrap', 'Bootstrap'],
    what: ['Se remuestrean las columnas del alineamiento con reemplazo, se rehace el árbol cientos de veces y se cuenta en qué porcentaje de esas réplicas aparece cada rama.',
      'The alignment columns are resampled with replacement, the tree is rebuilt hundreds of times, and the percentage of those replicates in which each branch appears is counted.'],
    read: ['Mide la <b>consistencia de la señal en los datos</b>, no la probabilidad de que la rama sea correcta. La costumbre es tomar 70 % como umbral de «bien apoyada» y 95 % como fuerte.',
      'It measures the <b>consistency of the signal in the data</b>, not the probability that the branch is right. The habit is to take 70 % as the threshold for "well supported" and 95 % as strong.'],
    scaleTitle: ['Bootstrap (%)', 'Bootstrap (%)'],
    scale: [
      S(null, 50, ['sin apoyo', 'unsupported'], 'bad'),
      S(50, 70, ['débil', 'weak'], 'warn'),
      S(70, 95, ['apoyada', 'supported'], 'ok'),
      S(95, null, ['fuerte', 'strong'], 'good'),
    ],
    conv: true,
    care: ['Es <b>conservador</b>: 70 % ya indica bastante apoyo. Pero no protege contra el error sistemático: si el modelo está mal, el bootstrap apoyará con firmeza la rama equivocada, y con más datos la apoyará más.',
      'It is <b>conservative</b>: 70 % already indicates fair support. But it does not protect against systematic error: if the model is wrong, the bootstrap will firmly support the wrong branch, and with more data it will support it more.'],
    ref: 'Felsenstein 1985; Hillis & Bull 1993',
  });

  E('bremer', {
    t: ['Soporte de Bremer (índice de decaimiento)', 'Bremer support (decay index)'],
    what: ['Cuántos pasos de más hay que aceptar antes de encontrar un árbol donde <b>esa rama ya no aparece</b>.',
      'How many extra steps have to be accepted before finding a tree in which <b>that branch no longer appears</b>.'],
    read: ['Se lee en pasos, no en porcentaje, y depende del tamaño de la matriz: en matrices grandes un Bremer de 5 puede ser poco. Un valor de 1 significa que basta un paso para deshacer la rama.',
      'It is read in steps, not in per cent, and depends on the size of the matrix: on large matrices a Bremer of 5 may be little. A value of 1 means one step is enough to undo the branch.'],
    care: ['No es comparable entre estudios sin normalizar por la longitud del árbol.',
      'It is not comparable between studies without normalising by the tree length.'],
    ref: 'Bremer 1994',
  });

  E('rf', {
    t: ['Distancia de Robinson–Foulds', 'Robinson–Foulds distance'],
    what: ['Cuenta las divisiones (bipartitions) que están en un árbol y no en el otro, en las dos direcciones. Cero significa la misma topología.',
      'Counts the splits (bipartitions) present in one tree and absent from the other, in both directions. Zero means the same topology.'],
    read: ['El máximo posible es 2(n−3) para árboles sin raíz de n puntas, así que conviene leer la <b>versión normalizada</b>, entre 0 y 1. Ignora las longitudes de rama por completo.',
      'The largest possible value is 2(n−3) for unrooted trees of n tips, so it is worth reading the <b>normalised version</b>, between 0 and 1. It ignores branch lengths entirely.'],
    scaleTitle: ['RF normalizada', 'Normalised RF'],
    scale: [
      S(0, 0.1, ['casi el mismo árbol', 'nearly the same tree'], 'good'),
      S(0.1, 0.3, ['diferencias menores', 'minor differences'], 'ok'),
      S(0.3, 0.6, ['discrepan', 'they disagree'], 'warn'),
      S(0.6, null, ['topologías distintas', 'different topologies'], 'bad'),
    ],
    conv: true,
    care: ['Es brusca: mover <b>una sola punta</b> de un extremo al otro del árbol puede llevar la RF casi al máximo. Dos árboles con RF grande pueden contar casi la misma historia.',
      'It is abrupt: moving <b>a single tip</b> from one end of the tree to the other can take RF nearly to its maximum. Two trees with a large RF can tell almost the same story.'],
    ref: 'Robinson & Foulds 1981',
  });

  /* ---------- Block 5 · maximum likelihood ---------- */

  E('ufboot', {
    t: ['Bootstrap ultrarrápido (UFBoot)', 'Ultrafast bootstrap (UFBoot)'],
    what: ['Una aproximación al bootstrap que reutiliza las verosimilitudes por sitio en lugar de reoptimizar cada réplica desde cero. Cientos de veces más rápida.',
      'An approximation to the bootstrap that reuses the per-site likelihoods instead of reoptimising every replicate from scratch. Hundreds of times faster.'],
    read: ['Su escala <b>no es la del bootstrap clásico</b>: está construida para ser casi insesgada, así que el umbral de «bien apoyada» es <b>95 %</b>, no 70 %. Confundir las dos escalas es el error más común con esta medida.',
      'Its scale is <b>not that of the classical bootstrap</b>: it is built to be nearly unbiased, so the threshold for "well supported" is <b>95 %</b>, not 70 %. Confusing the two scales is the commonest mistake with this measure.'],
    scaleTitle: ['UFBoot (%)', 'UFBoot (%)'],
    scale: [
      S(null, 85, ['sin apoyo', 'unsupported'], 'bad'),
      S(85, 95, ['dudosa', 'doubtful'], 'warn'),
      S(95, null, ['apoyada', 'supported'], 'good'),
    ],
    care: ['Al ser casi insesgada, es <b>menos</b> conservadora que el bootstrap clásico: un 90 % de UFBoot no equivale a un 90 % de bootstrap, sino a bastante menos apoyo.',
      'Being nearly unbiased, it is <b>less</b> conservative than the classical bootstrap: 90 % UFBoot is not the same as 90 % bootstrap, but considerably less support.'],
    ref: 'Minh, Nguyen & von Haeseler 2013',
  });

  E('shalrt', {
    t: ['SH-aLRT', 'SH-aLRT'],
    what: ['Una prueba rápida sobre cada rama interna: compara la verosimilitud del árbol con la de las dos reorganizaciones NNI alternativas de esa rama.',
      'A quick test on each internal branch: it compares the likelihood of the tree against that of the two alternative NNI rearrangements of that branch.'],
    read: ['Se lee en porcentaje, con <b>80 %</b> como umbral habitual. Es información distinta de la del bootstrap: mide si la rama mejora la verosimilitud, no si la señal se sostiene al remuestrear.',
      'It is read as a percentage, with <b>80 %</b> as the usual threshold. It is different information from the bootstrap: it measures whether the branch improves the likelihood, not whether the signal holds up under resampling.'],
    scaleTitle: ['SH-aLRT (%)', 'SH-aLRT (%)'],
    scale: [
      S(null, 80, ['sin apoyo', 'unsupported'], 'bad'),
      S(80, null, ['apoyada', 'supported'], 'good'),
    ],
    conv: true,
    care: ['Solo mira <b>una rama a la vez</b> y suponiendo el resto del árbol correcto. La costumbre es exigir las dos cosas: SH-aLRT ≥ 80 y UFBoot ≥ 95.',
      'It looks at <b>one branch at a time</b> and assumes the rest of the tree is right. The habit is to require both: SH-aLRT ≥ 80 and UFBoot ≥ 95.'],
    ref: 'Guindon et al. 2010',
  });

  E('autest', {
    t: ['Pruebas de topología (KH, SH, AU)', 'Topology tests (KH, SH, AU)'],
    what: ['Preguntan si un árbol alternativo —el de tu hipótesis, o el publicado— es <b>significativamente peor</b> que el mejor encontrado, o si la diferencia cabe dentro del ruido.',
      'They ask whether an alternative tree — your hypothesis, or a published one — is <b>significantly worse</b> than the best one found, or whether the difference fits within the noise.'],
    read: ['Un valor p por debajo de 0.05 rechaza el árbol alternativo. La prueba AU es la que corrige mejor el sesgo de selección y es la que conviene reportar.',
      'A p-value below 0.05 rejects the alternative tree. The AU test corrects the selection bias best and is the one worth reporting.'],
    scaleTitle: ['Valor p (AU)', 'p-value (AU)'],
    scale: [
      S(null, 0.05, ['se rechaza', 'rejected'], 'bad'),
      S(0.05, 0.1, ['al borde', 'borderline'], 'warn'),
      S(0.1, null, ['no se rechaza', 'not rejected'], 'good'),
    ],
    care: ['«No rechazado» no es «igual de bueno»: con pocos datos no se rechaza casi nada. Y el conjunto de árboles a comparar hay que fijarlo <b>antes</b> de mirar los resultados.',
      '"Not rejected" is not "just as good": with little data almost nothing is rejected. And the set of trees to compare has to be fixed <b>before</b> looking at the results.'],
    ref: 'Shimodaira 2002',
  });

  /* ---------- Block 6 · Bayesian inference ---------- */

  E('pp', {
    t: ['Probabilidad posterior', 'Posterior probability'],
    what: ['La probabilidad de que una rama exista, <b>dado el modelo, los priors y los datos</b>. A diferencia del bootstrap, es directamente una probabilidad.',
      'The probability that a branch exists, <b>given the model, the priors and the data</b>. Unlike the bootstrap, it is directly a probability.'],
    read: ['Se toma <b>0.95</b> como umbral de apoyo. Nunca se compara con el bootstrap en la misma escala: para una misma rama la posterior suele ser bastante más alta.',
      '<b>0.95</b> is taken as the support threshold. It is never compared with the bootstrap on the same scale: for the same branch the posterior is usually considerably higher.'],
    scaleTitle: ['Probabilidad posterior', 'Posterior probability'],
    scale: [
      S(null, 0.5, ['sin apoyo', 'unsupported'], 'bad'),
      S(0.5, 0.95, ['dudosa', 'doubtful'], 'warn'),
      S(0.95, null, ['apoyada', 'supported'], 'good'),
    ],
    conv: true,
    care: ['Es conocida por ser <b>demasiado segura</b> cuando el modelo está mal especificado: puede dar 1.00 a una rama que el bootstrap apoya al 60 %. Un 1.00 no es una garantía, es una afirmación condicionada al modelo.',
      'It is known to be <b>overconfident</b> when the model is misspecified: it can give 1.00 to a branch the bootstrap supports at 60 %. A 1.00 is not a guarantee, it is a statement conditional on the model.'],
    ref: 'Huelsenbeck & Rannala 2004; Douady et al. 2003',
  });

  E('ess', {
    t: ['Tamaño de muestra efectivo (ESS)', 'Effective sample size (ESS)'],
    what: ['Cuántas muestras <b>independientes</b> equivalen a las que la cadena tomó. Las muestras consecutivas de un MCMC están correlacionadas, así que diez mil muestras pueden valer por treinta.',
      'How many <b>independent</b> samples the ones the chain took are worth. Consecutive MCMC samples are correlated, so ten thousand samples can be worth thirty.'],
    read: ['La regla habitual es <b>ESS ≥ 200</b> para cada parámetro que se vaya a reportar. Por debajo, la estimación y su intervalo no son de fiar, por larga que se vea la corrida.',
      'The usual rule is <b>ESS ≥ 200</b> for every parameter to be reported. Below that, the estimate and its interval are not trustworthy, however long the run may look.'],
    scaleTitle: ['ESS', 'ESS'],
    scale: [
      S(null, 100, ['insuficiente', 'insufficient'], 'bad'),
      S(100, 200, ['justo', 'marginal'], 'warn'),
      S(200, null, ['suficiente', 'sufficient'], 'good'),
    ],
    conv: true,
    care: ['Un ESS alto dice que la cadena <b>mezcló bien donde estuvo</b>, no que haya encontrado todo el espacio. Una cadena atrapada en un óptimo local puede tener ESS excelente.',
      'A high ESS says the chain <b>mixed well where it was</b>, not that it found the whole space. A chain stuck in a local optimum can have an excellent ESS.'],
  });

  E('psrf', {
    t: ['PSRF (R̂, factor de reducción de escala)', 'PSRF (R̂, potential scale reduction)'],
    what: ['Compara la varianza <b>dentro</b> de cada cadena con la varianza <b>entre</b> cadenas independientes. Si todas exploran la misma distribución, las dos se parecen y el cociente tiende a 1.',
      'Compares the variance <b>within</b> each chain against the variance <b>between</b> independent chains. If they are all exploring the same distribution, the two are similar and the ratio tends to 1.'],
    read: ['Se exige <b>PSRF < 1.01</b> (algunos aceptan 1.05). Por encima, las cadenas no han convergido a lo mismo y no se pueden juntar.',
      '<b>PSRF < 1.01</b> is required (some accept 1.05). Above that, the chains have not converged on the same thing and cannot be pooled.'],
    scaleTitle: ['PSRF', 'PSRF'],
    scale: [
      S(null, 1.01, ['convergieron', 'converged'], 'good'),
      S(1.01, 1.05, ['al borde', 'borderline'], 'warn'),
      S(1.05, null, ['no convergieron', 'not converged'], 'bad'),
    ],
    conv: true,
    care: ['Necesita cadenas que empiecen en sitios <b>distintos</b>: arrancadas todas del mismo punto, pueden dar 1.00 y estar todas igual de equivocadas.',
      'It needs chains started from <b>different</b> places: all started from the same point, they can give 1.00 and all be equally wrong.'],
    ref: 'Gelman & Rubin 1992',
  });

  E('asdsf', {
    t: ['ASDSF', 'ASDSF'],
    what: ['La desviación estándar media de las frecuencias de cada división entre corridas independientes: cuánto se contradicen dos corridas sobre cuánto apoya cada rama.',
      'The average standard deviation of each split frequency across independent runs: how much two runs contradict each other about how much each branch is supported.'],
    read: ['Se busca <b>por debajo de 0.01</b>. Entre 0.01 y 0.05 la corrida es utilizable para explorar pero no para publicar.',
      'One looks for <b>below 0.01</b>. Between 0.01 and 0.05 the run is usable for exploring but not for publishing.'],
    scaleTitle: ['ASDSF', 'ASDSF'],
    scale: [
      S(null, 0.01, ['convergieron', 'converged'], 'good'),
      S(0.01, 0.05, ['aún no', 'not yet'], 'warn'),
      S(0.05, null, ['lejos', 'far from it'], 'bad'),
    ],
    conv: true,
    ref: 'Ronquist et al. 2012',
  });

  E('bayesfactor', {
    t: ['Factor de Bayes', 'Bayes factor'],
    what: ['El cociente de las verosimilitudes marginales de dos modelos: cuántas veces más probable hace un modelo a los datos observados que el otro. Se suele reportar como 2·ln(BF).',
      'The ratio of the marginal likelihoods of two models: how many times more probable one model makes the observed data than the other. Usually reported as 2·ln(BF).'],
    read: ['La escala de Kass y Raftery sobre 2·ln(BF): menos de 2, nada que decir; 2 a 6, positivo; 6 a 10, fuerte; más de 10, muy fuerte.',
      'The Kass and Raftery scale on 2·ln(BF): below 2, nothing to say; 2 to 6, positive; 6 to 10, strong; above 10, very strong.'],
    scaleTitle: ['2·ln(BF)', '2·ln(BF)'],
    scale: [
      S(null, 2, ['sin diferencia', 'no difference'], 'neutral'),
      S(2, 6, ['positivo', 'positive'], 'ok'),
      S(6, 10, ['fuerte', 'strong'], 'warn'),
      S(10, null, ['muy fuerte', 'very strong'], 'good'),
    ],
    conv: true,
    care: ['La verosimilitud marginal es <b>muy sensible al prior</b>, mucho más que la posterior. Dos personas con el mismo modelo y priors distintos pueden obtener factores de Bayes opuestos.',
      'The marginal likelihood is <b>very sensitive to the prior</b>, far more than the posterior. Two people with the same model and different priors can obtain opposite Bayes factors.'],
    ref: 'Kass & Raftery 1995; Xie et al. 2011',
  });

  E('burnin', {
    t: ['Quema (burn-in)', 'Burn-in'],
    what: ['Las primeras generaciones de la cadena, cuando todavía va subiendo desde donde la soltaron y no está muestreando la distribución posterior. Se descartan.',
      'The first generations of the chain, while it is still climbing from wherever it was released and is not yet sampling the posterior. They are discarded.'],
    read: ['Se mira la traza: donde deja de subir y empieza a oscilar alrededor de un valor, ahí empieza lo aprovechable. Descartar el 25 % es la costumbre, pero <b>mirar la traza es mejor que la costumbre</b>.',
      'One looks at the trace: where it stops climbing and begins to oscillate around a value, the usable part starts. Discarding 25 % is the habit, but <b>looking at the trace beats the habit</b>.'],
    care: ['Una traza plana no prueba convergencia: puede ser una cadena atascada. Por eso hacen falta varias corridas y el PSRF.',
      'A flat trace does not prove convergence: it can be a stuck chain. That is why several runs and the PSRF are needed.'],
  });

  /* ---------- Block 7 · the molecular clock and dating ---------- */

  E('rtt', {
    t: ['Regresión raíz–punta (r²)', 'Root-to-tip regression (r²)'],
    what: ['Se mide la distancia desde la raíz hasta cada punta y se compara con la fecha de esa punta. Si hay reloj, las dos suben juntas en línea recta y la pendiente <b>es la tasa</b>.',
      'The distance from the root to each tip is measured and compared with that tip\'s date. If there is a clock, the two rise together in a straight line and the slope <b>is the rate</b>.'],
    read: ['r² dice cuánta de la variación en distancia explica el tiempo. Solo sirve cuando las puntas tienen fechas distintas —virus, ADN antiguo, fósiles datados—: con puntas todas del presente no hay nada que regresar.',
      'r² says how much of the variation in distance is explained by time. It only works when the tips have different dates — viruses, ancient DNA, dated fossils: with tips all from the present there is nothing to regress.'],
    scaleTitle: ['r²', 'r²'],
    scale: [
      S(null, 0.2, ['sin señal temporal', 'no temporal signal'], 'bad'),
      S(0.2, 0.5, ['señal débil', 'weak signal'], 'warn'),
      S(0.5, 0.8, ['señal clara', 'clear signal'], 'ok'),
      S(0.8, null, ['reloj muy marcado', 'very clock-like'], 'good'),
    ],
    conv: true,
    care: ['La distancia raíz–punta de dos taxones hermanos <b>no es independiente</b>: comparten casi todo el camino. El r² de esta regresión es optimista y no debe usarse como una prueba formal.',
      'The root-to-tip distances of two sister taxa are <b>not independent</b>: they share almost the whole path. The r² of this regression is optimistic and should not be used as a formal test.'],
    ref: 'Rambaut et al. 2016',
  });

  E('clocktest', {
    t: ['Prueba del reloj estricto', 'Strict-clock test'],
    what: ['Compara la verosimilitud del árbol con las ramas libres contra la del mismo árbol forzado a ser ultramétrico. La diferencia se contrasta con una χ² de n−2 grados de libertad.',
      'Compares the likelihood of the tree with free branches against the same tree forced to be ultrametric. The difference is tested against a χ² with n−2 degrees of freedom.'],
    read: ['Un p pequeño rechaza el reloj estricto y manda usar un reloj relajado. No rechazarlo <b>no demuestra</b> que haya reloj: la prueba tiene poca potencia con pocos taxones.',
      'A small p rejects the strict clock and calls for a relaxed one. Failing to reject it <b>does not prove</b> there is a clock: the test has little power with few taxa.'],
    scaleTitle: ['Valor p', 'p-value'],
    scale: [
      S(null, 0.01, ['se rechaza el reloj', 'clock rejected'], 'bad'),
      S(0.01, 0.05, ['se rechaza', 'rejected'], 'warn'),
      S(0.05, null, ['no se rechaza', 'not rejected'], 'good'),
    ],
    care: ['Con muchos sitios la prueba rechaza casi siempre, porque detecta desviaciones diminutas y sin importancia práctica. Y la χ² es una aproximación asintótica.',
      'With many sites the test almost always rejects, because it detects tiny and practically unimportant deviations. And the χ² is an asymptotic approximation.'],
    ref: 'Felsenstein 1981',
  });

  E('hpd', {
    t: ['Intervalo de credibilidad (HPD 95 %)', 'Credible interval (95 % HPD)'],
    what: ['El intervalo más corto que contiene el 95 % de la distribución posterior de una edad. A diferencia de un intervalo de confianza, sí se lee como «hay un 95 % de probabilidad de que la edad esté aquí dentro».',
      'The shortest interval containing 95 % of the posterior distribution of an age. Unlike a confidence interval, it does read as "there is a 95 % probability the age lies in here".'],
    read: ['Lo que importa es su <b>anchura relativa</b> a la edad. Un nodo fechado en 30 Ma con un HPD de 28–32 está bien determinado; el mismo 30 Ma con 8–70 no dice casi nada, aunque la media se vea igual de redonda.',
      'What matters is its <b>width relative</b> to the age. A node dated at 30 Ma with an HPD of 28–32 is well determined; the same 30 Ma with 8–70 says almost nothing, although the mean looks just as tidy.'],
    care: ['El HPD recoge la incertidumbre <b>dentro del modelo y de las calibraciones que le diste</b>. No incluye el error de haber escogido mal una calibración, que suele ser el más grande de todos. Una edad nunca se reporta sin su intervalo.',
      'The HPD captures the uncertainty <b>within the model and the calibrations you gave it</b>. It does not include the error of having chosen a calibration badly, which is usually the largest of all. An age is never reported without its interval.'],
  });

  E('calibration', {
    t: ['Calibraciones', 'Calibrations'],
    what: ['Lo único que convierte cambios en años. Un fósil da la edad <b>mínima</b> de un nodo; una tasa publicada convierte directamente; una edad tomada de otro estudio arrastra la incertidumbre de aquel estudio entero.',
      'The only thing that turns changes into years. A fossil gives a node\'s <b>minimum</b> age; a published rate converts directly; an age taken from another study drags in that whole study\'s uncertainty.'],
    read: ['Un fósil casi nunca es la edad del nodo: es una cota inferior, y por eso se usan priors asimétricos (lognormal o exponencial desplazada) que permiten que el nodo sea mucho más viejo pero no más joven.',
      'A fossil is almost never the node\'s age: it is a lower bound, which is why asymmetric priors (lognormal or offset exponential) are used, allowing the node to be much older but not younger.'],
    care: ['Las fechas resultantes son tan buenas como la peor calibración, y dependen más de la <b>asignación filogenética del fósil</b> que de cualquier decisión estadística. Una calibración secundaria es una opinión heredada.',
      'The resulting dates are only as good as the worst calibration, and depend more on the <b>phylogenetic placement of the fossil</b> than on any statistical decision. A secondary calibration is an inherited opinion.'],
    ref: 'Parham et al. 2012',
  });

  E('smoothing', {
    t: ['Suavizado (λ) de la verosimilitud penalizada', 'Smoothing (λ) of penalised likelihood'],
    what: ['Controla cuánto se permite que la tasa cambie de una rama a su hija. λ grande se acerca a un reloj estricto; λ pequeña deja que cada rama tenga su propia tasa.',
      'Controls how much the rate is allowed to change from a branch to its daughter. A large λ approaches a strict clock; a small λ lets every branch have its own rate.'],
    read: ['No se elige a ojo: se escoge por <b>validación cruzada</b>, quitando una rama a la vez y viendo qué λ predice mejor. La curva de validación cruzada del panel es la que manda.',
      'It is not chosen by eye: it is chosen by <b>cross-validation</b>, removing one branch at a time and seeing which λ predicts best. The cross-validation curve in the panel is what decides.'],
    care: ['Si el mínimo de la curva cae en un extremo del rango probado, el rango era demasiado estrecho y hay que ampliarlo.',
      'If the minimum of the curve falls at one end of the range tried, the range was too narrow and has to be widened.'],
    ref: 'Sanderson 2002',
  });

  /* ---------- Block 8 · diversification ---------- */

  E('gamma', {
    t: ['Estadístico γ de Pybus y Harvey', 'Pybus & Harvey γ statistic'],
    what: ['Compara dónde caen los nodos de un árbol respecto de donde caerían con una tasa constante. Negativo: los nodos se acumulan cerca de la raíz, señal de <b>desaceleración</b>. Positivo: se acumulan cerca del presente.',
      'Compares where a tree\'s nodes fall against where they would fall under a constant rate. Negative: the nodes pile up near the root, a sign of <b>slowdown</b>. Positive: they pile up near the present.'],
    read: ['Bajo tasa constante γ se distribuye normal estándar, así que <b>γ < −1.645</b> rechaza la tasa constante a favor de una desaceleración, con una cola.',
      'Under a constant rate γ has a standard normal distribution, so <b>γ < −1.645</b> rejects the constant rate in favour of a slowdown, one-tailed.'],
    scaleTitle: ['γ', 'γ'],
    scale: [
      S(null, -1.645, ['desaceleración', 'slowdown'], 'warn'),
      S(-1.645, 1.645, ['compatible con tasa constante', 'compatible with a constant rate'], 'good'),
      S(1.645, null, ['aceleración', 'speed-up'], 'warn'),
    ],
    care: ['<b>Faltar especies imita exactamente una desaceleración</b>: si no muestreaste todo el grupo, γ saldrá negativo aunque la tasa sea constante. Por eso existe la prueba MCCR, que es la que hay que usar cuando el muestreo es incompleto.',
      '<b>Missing species imitate a slowdown exactly</b>: if you did not sample the whole group, γ will come out negative even under a constant rate. That is what the MCCR test is for, and it is the one to use when sampling is incomplete.'],
    ref: 'Pybus & Harvey 2000',
  });

  E('birthdeath', {
    t: ['Nacimiento–muerte: λ, μ, r y ε', 'Birth–death: λ, μ, r and ε'],
    what: ['λ es la tasa de especiación, μ la de extinción, r = λ−μ la <b>diversificación neta</b> y ε = μ/λ el <b>recambio</b>. El árbol reconstruido solo contiene linajes que sobrevivieron, y de ahí sale casi toda la dificultad.',
      'λ is the speciation rate, μ the extinction rate, r = λ−μ the <b>net diversification</b> and ε = μ/λ the <b>turnover</b>. The reconstructed tree contains only lineages that survived, and nearly all the difficulty comes from that.'],
    read: ['r se estima razonablemente bien; <b>μ y ε son notoriamente poco fiables</b> con árboles solo de especies vivas. Reportar r con confianza y ε con reservas es lo honesto.',
      'r is estimated reasonably well; <b>μ and ε are notoriously unreliable</b> from trees of living species alone. Reporting r with confidence and ε with reservations is the honest course.'],
    scaleTitle: ['Recambio ε = μ/λ', 'Turnover ε = μ/λ'],
    scale: [
      S(0, 0.2, ['poca extinción', 'little extinction'], 'good'),
      S(0.2, 0.6, ['moderada', 'moderate'], 'ok'),
      S(0.6, 0.9, ['alta', 'high'], 'warn'),
      S(0.9, 1, ['casi todo se extingue', 'almost everything goes extinct'], 'bad'),
    ],
    conv: true,
    care: ['Una ε estimada en 0 no significa que no hubo extinción: significa que el árbol de los supervivientes no puede verla. Con fósiles la estimación cambia por completo.',
      'An ε estimated at 0 does not mean there was no extinction: it means the tree of the survivors cannot see it. With fossils the estimate changes completely.'],
    ref: 'Nee, May & Harvey 1994; Rabosky 2010',
  });

  E('dr', {
    t: ['Estadístico DR', 'DR statistic'],
    what: ['Reparte una tasa entre las puntas: para cada especie pesa las ramas del camino hasta la raíz, dando más peso a las recientes. Sirve para comparar especies, no para estimar una tasa del grupo.',
      'Splits a rate among the tips: for each species it weights the branches on the path to the root, giving more weight to the recent ones. It is for comparing species, not for estimating a rate for the group.'],
    read: ['Se lee <b>relativo</b>: qué puntas tienen DR alto frente a las demás del mismo árbol. Los valores absolutos dependen de la escala del árbol.',
      'It is read <b>relatively</b>: which tips have a high DR compared with the rest of the same tree. Absolute values depend on the tree\'s scale.'],
    care: ['Se ha demostrado que DR estima mejor la <b>especiación</b> que la diversificación neta, y que se comporta mal cuando la extinción es alta. No es una tasa de diversificación por especie, aunque se use así a menudo.',
      'DR has been shown to estimate <b>speciation</b> better than net diversification, and to behave badly when extinction is high. It is not a per-species diversification rate, although it is often used as one.'],
    ref: 'Jetz et al. 2012; Title & Rabosky 2019',
  });

  /* ---------- Block 9 · trait evolution ---------- */

  E('pagel', {
    t: ['λ de Pagel', 'Pagel\'s λ'],
    what: ['Multiplica las covarianzas entre especies —lo que comparten por descendencia— sin tocar las varianzas. λ = 1 es movimiento browniano puro; λ = 0 es un carácter repartido como si la filogenia no existiera.',
      'Multiplies the covariances between species — what they share by descent — without touching the variances. λ = 1 is pure Brownian motion; λ = 0 is a character distributed as though the phylogeny did not exist.'],
    read: ['Es una medida de <b>señal filogenética</b>: cuánto se parecen los parientes más de lo esperado por azar. Se prueba contra 0 (¿hay señal?) y contra 1 (¿es browniano?).',
      'It is a measure of <b>phylogenetic signal</b>: how much relatives resemble each other beyond chance. It is tested against 0 (is there signal?) and against 1 (is it Brownian?).'],
    scaleTitle: ['λ', 'λ'],
    scale: [
      S(0, 0.2, ['sin señal', 'no signal'], 'bad'),
      S(0.2, 0.6, ['señal parcial', 'partial signal'], 'warn'),
      S(0.6, 0.9, ['buena señal', 'good signal'], 'ok'),
      S(0.9, 1, ['browniano', 'Brownian'], 'good'),
    ],
    conv: true,
    care: ['λ mide el <b>patrón</b>, no el proceso: λ = 1 no demuestra deriva neutra, y λ bajo puede deberse a error de medición y no a la biología. Con menos de veinte especies su estimación es muy imprecisa.',
      'λ measures the <b>pattern</b>, not the process: λ = 1 does not prove neutral drift, and a low λ may come from measurement error rather than biology. With fewer than twenty species its estimate is very imprecise.'],
    ref: 'Pagel 1999; Freckleton et al. 2002',
  });

  E('blombergk', {
    t: ['K de Blomberg', 'Blomberg\'s K'],
    what: ['Compara la varianza observada entre especies con la que el movimiento browniano predice sobre ese árbol. K = 1 es exactamente browniano; K > 1 significa que los parientes se parecen <b>más</b> de lo que el browniano espera.',
      'Compares the observed variance among species against what Brownian motion predicts on that tree. K = 1 is exactly Brownian; K > 1 means relatives resemble each other <b>more</b> than Brownian motion expects.'],
    read: ['Se acompaña siempre de una prueba de permutación: K sin su valor p no dice si la señal es distinguible del azar.',
      'It always comes with a permutation test: K without its p-value does not say whether the signal is distinguishable from chance.'],
    scaleTitle: ['K', 'K'],
    scale: [
      S(0, 0.5, ['menos señal que el browniano', 'less signal than Brownian'], 'warn'),
      S(0.5, 1, ['algo menos', 'somewhat less'], 'ok'),
      S(1, null, ['más que el browniano', 'more than Brownian'], 'good'),
    ],
    conv: true,
    care: ['K y λ pueden parecer contradictorios y no serlo: <b>λ mira la forma del árbol y K la magnitud de la varianza</b>. Son preguntas distintas sobre los mismos datos.',
      'K and λ can look contradictory without being so: <b>λ looks at the shape of the tree and K at the magnitude of the variance</b>. They are different questions about the same data.'],
    ref: 'Blomberg, Garland & Ives 2003',
  });

  E('ou', {
    t: ['Ornstein–Uhlenbeck: α y vida media', 'Ornstein–Uhlenbeck: α and half-life'],
    what: ['Un modelo con una fuerza que jala el carácter hacia un óptimo. α es esa fuerza; la <b>vida media</b> = ln(2)/α es el tiempo que tarda en recorrer la mitad del camino hacia el óptimo, y se lee mucho mejor que α.',
      'A model with a force pulling the character towards an optimum. α is that force; the <b>half-life</b> = ln(2)/α is the time it takes to cover half the way to the optimum, and it reads far better than α.'],
    read: ['Compara la vida media con la <b>profundidad del árbol</b>. Mucho menor: selección estabilizadora fuerte, el árbol casi no importa. Mayor que la profundidad: el OU es indistinguible de un browniano.',
      'Compare the half-life against the <b>depth of the tree</b>. Much smaller: strong stabilising selection, the tree hardly matters. Larger than the depth: the OU is indistinguishable from Brownian motion.'],
    care: ['El OU supone <b>puntas contemporáneas</b>: en un filograma no ultramétrico sus parámetros no significan lo que dicen. Y con menos de cincuenta especies el OU se prefiere a menudo por razones estadísticas y no biológicas: es un modelo muy hábil para ajustar ruido.',
      'The OU assumes <b>contemporaneous tips</b>: on a non-ultrametric phylogram its parameters do not mean what they say. And with fewer than fifty species the OU is often preferred for statistical rather than biological reasons: it is very good at fitting noise.'],
    ref: 'Hansen 1997; Cooper et al. 2016',
  });

  E('pgls', {
    t: ['PGLS', 'PGLS'],
    what: ['Una regresión que corrige el hecho de que las especies <b>no son observaciones independientes</b>: dos especies hermanas comparten casi toda su historia, y tratarlas como puntos sueltos infla la significancia.',
      'A regression that corrects for the fact that species are <b>not independent observations</b>: two sister species share nearly all their history, and treating them as separate points inflates the significance.'],
    read: ['Se lee como cualquier regresión —pendiente, error estándar, p— pero el p ya está corregido por el parentesco. Comparar el p de PGLS con el de una regresión ordinaria enseña cuánto se estaba exagerando.',
      'It reads like any regression — slope, standard error, p — but the p is already corrected for relatedness. Comparing the PGLS p with that of an ordinary regression shows how much was being exaggerated.'],
    care: ['Corrige la <b>no independencia</b>, no la causalidad ni la variable omitida. Y si la señal filogenética del residuo es nula, PGLS y mínimos cuadrados dan lo mismo: la corrección no siempre cambia la conclusión.',
      'It corrects <b>non-independence</b>, not causality nor an omitted variable. And if the phylogenetic signal of the residual is nil, PGLS and ordinary least squares agree: the correction does not always change the conclusion.'],
    ref: 'Grafen 1989; Freckleton 2009',
  });

  E('mk', {
    t: ['Modelo Mk (ER, SYM, ARD)', 'Mk model (ER, SYM, ARD)'],
    what: ['El modelo de un carácter discreto sobre el árbol. ER: una sola tasa para todos los cambios. SYM: la tasa de i a j igual que la de j a i. ARD: todas distintas.',
      'The model of a discrete character on the tree. ER: a single rate for every change. SYM: the rate from i to j equals that from j to i. ARD: all different.'],
    read: ['Se comparan por AIC como cualquier otro modelo. ARD tiene k(k−1) parámetros: con cuatro estados son doce, y raras veces hay datos para estimarlos.',
      'They are compared by AIC like any other model. ARD has k(k−1) parameters: with four states that is twelve, and there is rarely data enough to estimate them.'],
    care: ['Con pocas especies casi siempre gana ER, y no porque la biología sea simple. Una tasa que se va a infinito es el síntoma de un carácter <b>saturado</b>: los nodos salen todos al 0.5 y la reconstrucción no dice nada.',
      'With few species ER almost always wins, and not because the biology is simple. A rate running off to infinity is the symptom of a <b>saturated</b> character: every node comes out at 0.5 and the reconstruction says nothing.'],
    ref: 'Lewis 2001; Pagel 1994',
  });

  E('ancstate', {
    t: ['Estados ancestrales', 'Ancestral states'],
    what: ['La probabilidad de cada estado en cada nodo interno, dibujada como un pastel. La versión <b>marginal</b> —la que se reporta— integra sobre todo lo demás del árbol.',
      'The probability of each state at each internal node, drawn as a pie. The <b>marginal</b> version — the one reported — integrates over everything else on the tree.'],
    read: ['Un pastel casi de un solo color es una reconstrucción firme; uno repartido significa que <b>los datos no saben</b>, y ese reparto es el resultado, no un defecto del dibujo.',
      'A pie of almost one colour is a firm reconstruction; a divided one means <b>the data do not know</b>, and that division is the result, not a defect of the drawing.'],
    care: ['La incertidumbre crece hacia la raíz, y la raíz suele ser el nodo que a uno le interesa. Las reconstrucciones profundas rara vez son firmes, y presentarlas sin el pastel es esconder el resultado.',
      'Uncertainty grows towards the root, and the root is usually the node one cares about. Deep reconstructions are rarely firm, and presenting them without the pie hides the result.'],
  });

  /* ---------- Block 10 · historical biogeography ---------- */

  E('decj', {
    t: ['DEC, DIVALIKE, BAYAREALIKE y el parámetro +J', 'DEC, DIVALIKE, BAYAREALIKE and the +J parameter'],
    what: ['Tres modelos que difieren en <b>qué permiten al dividirse un linaje</b>: DEC admite simpatría estrecha, de subconjunto y vicarianza; DIVALIKE favorece la vicarianza; BAYAREALIKE deja el rango idéntico. El +J añade la <b>especiación fundadora</b>: un salto a un área nueva en el momento de la división.',
      'Three models differing in <b>what they allow when a lineage splits</b>: DEC admits narrow sympatry, subset sympatry and vicariance; DIVALIKE favours vicariance; BAYAREALIKE leaves the range identical. The +J adds <b>founder-event speciation</b>: a jump to a new area at the moment of the split.'],
    read: ['Se comparan por AICc con pesos de Akaike. En islas el +J casi siempre gana, y con mucha diferencia.',
      'They are compared by AICc with Akaike weights. On islands the +J almost always wins, and by a wide margin.'],
    care: ['Ree y Sanmartín (2018) mostraron que <b>comparar DEC con DEC+J por verosimilitud está sesgado</b>: el +J puede ganar sin que el proceso exista, porque los dos modelos no son estrictamente anidados de la forma que la comparación supone. Reporta el +J, pero no lo presentes como demostración de dispersión fundadora.',
      'Ree and Sanmartín (2018) showed that <b>comparing DEC with DEC+J by likelihood is biased</b>: the +J can win without the process existing, because the two models are not nested in the way the comparison assumes. Report the +J, but do not present it as a demonstration of founder dispersal.'],
    ref: 'Ree & Smith 2008; Matzke 2013; Ree & Sanmartín 2018',
  });

  E('ancrange', {
    t: ['Rangos ancestrales', 'Ancestral ranges'],
    what: ['La probabilidad de cada combinación de áreas en cada nodo. El número de estados crece como 2 elevado al número de áreas, así que con ocho áreas ya hay 256 rangos posibles.',
      'The probability of each combination of areas at each node. The number of states grows as 2 to the power of the number of areas, so with eight areas there are already 256 possible ranges.'],
    read: ['Lo mismo que con los estados ancestrales: un nodo repartido entre muchos rangos es un nodo del que no se sabe. Conviene fijarse en el rango <b>más probable</b> y en cuánta probabilidad se lleva, no solo en el nombre del ganador.',
      'The same as with ancestral states: a node split among many ranges is a node about which nothing is known. It is worth looking at the <b>most probable</b> range and at how much probability it takes, not only at the winner\'s name.'],
    care: ['El resultado depende muchísimo de <b>cómo definiste las áreas</b>, que es una decisión tuya y no del programa. Cambiar el número de áreas cambia las conclusiones más que cambiar de modelo.',
      'The result depends enormously on <b>how you defined the areas</b>, which is your decision and not the program\'s. Changing the number of areas changes the conclusions more than changing the model.'],
  });

  /* ---------- Block 11 · comparing trees ---------- */

  E('gcf', {
    t: ['Factor de concordancia génica (gCF)', 'Gene concordance factor (gCF)'],
    what: ['El porcentaje de árboles de genes que contienen <b>esa misma rama</b>. Es un recuento directo, no una prueba estadística.',
      'The percentage of gene trees that contain <b>that same branch</b>. It is a direct count, not a statistical test.'],
    read: ['Responde una pregunta distinta de la del bootstrap. Una rama puede tener bootstrap del 100 % y gCF del 35 %: la concatenación está segura, y aun así <b>dos tercios de los genes dicen otra cosa</b>. Esa combinación es la señal clásica de conflicto real entre loci.',
      'It answers a different question from the bootstrap. A branch can have 100 % bootstrap and 35 % gCF: the concatenation is certain, and yet <b>two thirds of the genes say otherwise</b>. That combination is the classic signature of real conflict between loci.'],
    scaleTitle: ['gCF (%)', 'gCF (%)'],
    scale: [
      S(null, 33, ['la mayoría discrepa', 'most genes disagree'], 'bad'),
      S(33, 50, ['conflicto fuerte', 'strong conflict'], 'warn'),
      S(50, 75, ['mayoría a favor', 'majority in favour'], 'ok'),
      S(75, null, ['acuerdo amplio', 'broad agreement'], 'good'),
    ],
    conv: true,
    care: ['Un gCF bajo puede venir de conflicto biológico real —clasificación incompleta de linajes, hibridación— o simplemente de que cada gen tiene poca señal. Mirar también el sCF ayuda a distinguirlos.',
      'A low gCF can come from real biological conflict — incomplete lineage sorting, hybridisation — or simply from each gene having little signal. Looking at the sCF as well helps tell them apart.'],
    ref: 'Minh, Hahn & Lanfear 2020',
  });

  E('scf', {
    t: ['Factor de concordancia de sitios (sCF)', 'Site concordance factor (sCF)'],
    what: ['El porcentaje de sitios <b>informativos</b> que apoyan esa rama frente a las dos alternativas. Se puede calcular aunque haya un solo gen.',
      'The percentage of <b>informative</b> sites supporting that branch against the two alternatives. It can be computed even with a single gene.'],
    read: ['Su mínimo práctico ronda el <b>33 %</b>, porque hay tres resoluciones posibles: un sCF cerca de 33 significa que los sitios están repartidos al azar entre las tres, es decir que no hay señal.',
      'Its practical floor is about <b>33 %</b>, because there are three possible resolutions: an sCF near 33 means the sites are split at random among the three, that is, there is no signal.'],
    scaleTitle: ['sCF (%)', 'sCF (%)'],
    scale: [
      S(null, 36, ['sin señal', 'no signal'], 'bad'),
      S(36, 45, ['muy repartida', 'very split'], 'warn'),
      S(45, 60, ['mayoría', 'majority'], 'ok'),
      S(60, null, ['clara', 'clear'], 'good'),
    ],
    conv: true,
    ref: 'Minh, Hahn & Lanfear 2020',
  });

  E('dstat', {
    t: ['Estadístico D de Patterson (ABBA-BABA)', 'Patterson\'s D (ABBA-BABA)'],
    what: ['Con cuatro taxones ((P1,P2),P3),O), cuenta los sitios ABBA y BABA. Bajo clasificación incompleta de linajes sola, los dos patrones deberían ser <b>igual de frecuentes</b>. D mide el desequilibrio.',
      'With four taxa ((P1,P2),P3),O), it counts ABBA and BABA sites. Under incomplete lineage sorting alone, the two patterns should be <b>equally frequent</b>. D measures the imbalance.'],
    read: ['D distinto de cero indica flujo génico: positivo entre P2 y P3, negativo entre P1 y P3. Lo que decide es el <b>z</b> del jackknife por bloques: |z| > 3 es la convención.',
      'A D different from zero indicates gene flow: positive between P2 and P3, negative between P1 and P3. What decides is the block-jackknife <b>z</b>: |z| > 3 is the convention.'],
    scaleTitle: ['|z| del jackknife', 'Jackknife |z|'],
    scale: [
      S(null, 2, ['sin evidencia', 'no evidence'], 'good'),
      S(2, 3, ['sugerente', 'suggestive'], 'warn'),
      S(3, null, ['flujo génico', 'gene flow'], 'bad'),
    ],
    conv: true,
    care: ['El jackknife tiene que ser <b>por bloques</b>, porque los sitios cercanos están ligados y contarlos como independientes infla el z. Y un D significativo dice que hubo flujo <b>en algún lugar de ese cuarteto</b>, no necesariamente entre los dos taxones que nombra.',
      'The jackknife has to be <b>by blocks</b>, because nearby sites are linked and counting them as independent inflates z. And a significant D says there was flow <b>somewhere in that quartet</b>, not necessarily between the two taxa it names.'],
    ref: 'Green et al. 2010; Durand et al. 2011',
  });

  /* ---------- Block 12 · the figure and the report ---------- */

  E('dpi', {
    t: ['Resolución (ppp) y píxeles', 'Resolution (dpi) and pixels'],
    what: ['Los puntos por pulgada con que se convierte la figura vectorial en imagen. El número de píxeles crece con el <b>cuadrado</b> de la resolución: duplicar los ppp cuadruplica el archivo.',
      'The dots per inch at which the vector figure is turned into an image. The pixel count grows with the <b>square</b> of the resolution: doubling the dpi quadruples the file.'],
    read: ['300 ppp es lo que piden casi todas las revistas para una figura de línea; 600 para una con detalle fino. El panel calcula los píxeles <b>antes</b> de generar nada y avisa si no caben en un lienzo de navegador.',
      '300 dpi is what nearly every journal asks for a line figure; 600 for one with fine detail. The panel computes the pixels <b>before</b> generating anything and warns when they do not fit in a browser canvas.'],
    care: ['Subir los ppp no añade información a una figura vectorial: el SVG ya tiene resolución infinita. Si la revista acepta vectores, mándale el SVG o el PDF y olvídate de los ppp.',
      'Raising the dpi adds no information to a vector figure: the SVG already has infinite resolution. If the journal accepts vectors, send the SVG or the PDF and forget about dpi.'],
  });

  /* ---------- a second pass: what the panels print and had no entry ---------- */

  E('trimming', {
    t: ['Recorte de bloques ambiguos', 'Trimming of ambiguous blocks'],
    what: ['Quitar las columnas donde el alineamiento no es de fiar: zonas con muchos huecos, sin conservación, o rodeadas de otras igual de dudosas. Gblocks las escoge por reglas de vecindad; los umbrales tipo trimAl, por columna.',
      'Removing the columns where the alignment is not trustworthy: stretches with many gaps, without conservation, or surrounded by equally doubtful ones. Gblocks picks them by neighbourhood rules; trimAl-style thresholds, column by column.'],
    read: ['Se mira <b>cuánto se quita</b>. Perder un 10–30 % es normal. Perder más de la mitad significa que el alineamiento era malo, y el recorte no lo arregla: lo esconde.',
      'What is looked at is <b>how much goes</b>. Losing 10–30 % is normal. Losing more than half means the alignment was bad, and trimming does not fix that: it hides it.'],
    scaleTitle: ['Columnas eliminadas', 'Columns removed'],
    scale: [
      S(0, 0.1, ['casi nada', 'almost nothing'], 'good'),
      S(0.1, 0.3, ['normal', 'normal'], 'ok'),
      S(0.3, 0.5, ['mucho', 'a lot'], 'warn'),
      S(0.5, 1, ['el alineamiento era el problema', 'the alignment was the problem'], 'bad'),
    ],
    conv: true,
    care: ['Recortar no siempre mejora el árbol: se han publicado casos en que quitar columnas <b>empeora</b> la inferencia, porque se van también sitios informativos. Conviene comparar el árbol con y sin recorte.',
      'Trimming does not always improve the tree: cases have been published where removing columns <b>worsens</b> the inference, because informative sites go too. It is worth comparing the tree with and without.'],
    ref: 'Castresana 2000; Tan et al. 2015',
  });

  E('variable', {
    t: ['Sitios constantes, variables e informativos', 'Constant, variable and informative sites'],
    what: ['Constante: todos los taxones tienen lo mismo. Variable: hay al menos dos estados. Informativo: hay dos estados <b>cada uno en al menos dos taxones</b>. Cada categoría contiene a la siguiente.',
      'Constant: every taxon has the same thing. Variable: there are at least two states. Informative: there are two states <b>each in at least two taxa</b>. Each category contains the next.'],
    read: ['La cifra que decide es la de <b>informativos</b>: los variables no informativos —un solo taxón distinto— cuestan un paso en cualquier árbol y no eligen topología. Los constantes tampoco, pero sí informan sobre las tasas en verosimilitud.',
      'The figure that decides is the <b>informative</b> one: variable but uninformative sites — a single differing taxon — cost one step on any tree and choose no topology. Constant sites do not either, but they do inform the rates under likelihood.'],
    care: ['Quitar los sitios constantes «porque no aportan» sesga las longitudes de rama hacia arriba. Solo se quitan cuando el modelo lo sabe y lo corrige, que es lo que hace Mkv en morfología.',
      'Removing constant sites "because they add nothing" biases branch lengths upwards. They are only removed when the model knows and corrects for it, which is what Mkv does in morphology.'],
  });

  E('gc', {
    t: ['Contenido G+C', 'G+C content'],
    what: ['La proporción de guaninas y citosinas de una secuencia. Varía entre linajes, entre genomas y entre posiciones del codón.',
      'The proportion of guanines and cytosines in a sequence. It varies between lineages, between genomes and between codon positions.'],
    read: ['Lo que importa no es el valor sino <b>cuánto difieren las secuencias entre sí</b>. Dos taxones no emparentados con G+C parecido se atraen en el árbol, porque comparten estados por composición y no por ascendencia.',
      'What matters is not the value but <b>how much the sequences differ from each other</b>. Two unrelated taxa with a similar G+C attract each other on the tree, because they share states through composition and not through descent.'],
    care: ['Casi todos los modelos suponen una composición <b>constante en todo el árbol</b>. Cuando la prueba χ² la rechaza, el remedio no es cambiar de modelo de sustitución sino quitar la tercera posición, usar aminoácidos o recodificar.',
      'Nearly every model assumes a composition <b>constant across the whole tree</b>. When the χ² test rejects it, the remedy is not a different substitution model but dropping the third position, switching to amino acids, or recoding.'],
  });

  E('titv', {
    t: ['Razón transiciones/transversiones (ts/tv, κ)', 'Transition/transversion ratio (ts/tv, κ)'],
    what: ['Las transiciones (A↔G, C↔T) cambian una purina por otra purina, o una pirimidina por otra; las transversiones cruzan los dos grupos. Las primeras ocurren mucho más a menudo aunque haya el doble de transversiones posibles.',
      'Transitions (A↔G, C↔T) swap a purine for a purine, or a pyrimidine for a pyrimidine; transversions cross the two groups. The former happen far more often although twice as many transversions are possible.'],
    read: ['En ADN nuclear y de cloroplasto lo normal ronda 2; en mitocondrial animal puede pasar de 10. Un valor cercano a 1 en un gen que debería tener más es señal de <b>saturación</b>: las transiciones, que se saturan primero, ya se borraron.',
      'In nuclear and chloroplast DNA the usual value is around 2; in animal mitochondrial DNA it can exceed 10. A value near 1 in a gene that should have more is a sign of <b>saturation</b>: the transitions, which saturate first, have already been erased.'],
    scaleTitle: ['κ', 'κ'],
    scale: [
      S(null, 1.2, ['sospechosamente baja', 'suspiciously low'], 'warn'),
      S(1.2, 4, ['normal', 'normal'], 'good'),
      S(4, 10, ['alta', 'high'], 'ok'),
      S(10, null, ['muy alta', 'very high'], 'ok'),
    ],
    conv: true,
    ref: 'Hasegawa, Kishino & Yano 1985',
  });

  E('lrt', {
    t: ['Prueba de razón de verosimilitudes (LRT)', 'Likelihood ratio test (LRT)'],
    what: ['Compara dos modelos <b>anidados</b> —uno es caso particular del otro— con el estadístico 2ΔlnL, contrastado contra una χ² con tantos grados de libertad como parámetros de diferencia.',
      'Compares two <b>nested</b> models — one a special case of the other — with the statistic 2ΔlnL, tested against a χ² with as many degrees of freedom as the difference in parameters.'],
    read: ['Un p pequeño dice que el modelo más complejo <b>vale lo que cuesta</b>. La cadena jerárquica de pruebas termina en el modelo más simple que no se rechaza.',
      'A small p says the more complex model <b>is worth its cost</b>. The hierarchical chain of tests ends at the simplest model that is not rejected.'],
    scaleTitle: ['Valor p', 'p-value'],
    scale: [
      S(null, 0.01, ['el complejo gana claro', 'the complex one clearly wins'], 'good'),
      S(0.01, 0.05, ['el complejo gana', 'the complex one wins'], 'ok'),
      S(0.05, null, ['basta el simple', 'the simple one suffices'], 'neutral'),
    ],
    care: ['El resultado <b>depende del orden</b> en que se encadenan las pruebas, y eso es una debilidad conocida: por eso AIC y BIC son hoy la vía preferida. Además la χ² no vale cuando el parámetro queda en el borde de su rango, que es el caso de +I.',
      'The result <b>depends on the order</b> in which the tests are chained, and that is a known weakness: which is why AIC and BIC are the preferred route today. Besides, the χ² does not hold when a parameter sits on the boundary of its range, which is the case for +I.'],
    ref: 'Posada & Crandall 2001',
  });

  E('distance', {
    t: ['Distancias corregidas', 'Corrected distances'],
    what: ['La proporción de posiciones distintas (distancia p) subestima los cambios reales, porque un sitio que cambió dos veces se ve como uno o como ninguno. Los modelos —JC69, K80, TN93, LogDet— estiman cuántos cambios hubo <b>de verdad</b>.',
      'The proportion of differing positions (the p-distance) underestimates the real changes, because a site that changed twice looks like one change or none. The models — JC69, K80, TN93, LogDet — estimate how many changes there really were.'],
    read: ['La corrección crece rápido: una p de 0.5 puede corresponder a más de un cambio por sitio. Compara la distancia media con la máxima: si la máxima se dispara, ese par está en la zona donde la corrección deja de tener solución.',
      'The correction grows fast: a p of 0.5 can correspond to more than one change per site. Compare the mean distance with the maximum: if the maximum runs away, that pair is in the zone where the correction stops having a solution.'],
    care: ['Cada corrección supone su modelo. Con secuencias muy divergentes la fórmula puede no tener solución y devolver infinito: eso no es un error del programa, es el dato diciendo que ya no sabe.',
      'Each correction assumes its model. With very divergent sequences the formula may have no solution and return infinity: that is not a bug, it is the data saying it no longer knows.'],
  });

  E('equaltrees', {
    t: ['Árboles igual de cortos', 'Equally parsimonious trees'],
    what: ['Muchas matrices no tienen un solo árbol más corto sino decenas o miles, todos con el mismo número de pasos. La parsimonia no prefiere ninguno.',
      'Many matrices have not one shortest tree but tens or thousands, all with the same number of steps. Parsimony prefers none of them.'],
    read: ['Cuando hay varios, lo que se reporta es su <b>consenso</b>, y las politomías del consenso son partes del árbol que los datos no resuelven. Un solo árbol más corto es poco común y conviene desconfiar si la búsqueda fue corta.',
      'When there are several, what is reported is their <b>consensus</b>, and the polytomies of the consensus are parts of the tree the data do not resolve. A single shortest tree is uncommon, and worth doubting if the search was short.'],
    care: ['Encontrar «un solo árbol» puede significar que la búsqueda no exploró bastante. Y un consenso estricto de miles de árboles puede quedar casi todo colapsado, lo cual es el resultado honesto y no un fracaso.',
      'Finding "a single tree" can mean the search did not explore enough. And a strict consensus of thousands of trees can end up almost entirely collapsed, which is the honest result and not a failure.'],
  });

  E('search', {
    t: ['Búsqueda en el espacio de árboles (NNI, SPR, TBR)', 'Tree search (NNI, SPR, TBR)'],
    what: ['El número de topologías crece más rápido que factorialmente: con 20 taxones hay más de 10²⁰ árboles sin raíz. No se pueden evaluar todos, así que se parte de uno y se prueban reordenamientos. NNI mueve poco, SPR arranca una rama y la reinjerta, TBR corta el árbol en dos y lo vuelve a unir de todas las formas.',
      'The number of topologies grows faster than factorially: with 20 taxa there are more than 10²⁰ unrooted trees. They cannot all be evaluated, so one starts from a tree and tries rearrangements. NNI moves little, SPR tears off a branch and regrafts it, TBR cuts the tree in two and rejoins it every way.'],
    read: ['A más alcance del reordenamiento, menos riesgo de quedarse en un óptimo local y más tiempo. Varios puntos de partida distintos que lleguen al mismo árbol es la mejor señal de que la búsqueda fue suficiente.',
      'The wider the rearrangement, the smaller the risk of stopping at a local optimum and the longer it takes. Several different starting points arriving at the same tree is the best sign the search was enough.'],
    care: ['Ninguna búsqueda heurística <b>garantiza</b> el óptimo. Que dos programas den árboles distintos suele ser esto y no un desacuerdo de método: son dos óptimos locales del mismo criterio.',
      'No heuristic search <b>guarantees</b> the optimum. Two programs giving different trees is usually this and not a disagreement of method: they are two local optima of the same criterion.'],
  });

  E('partrate', {
    t: ['Velocidad relativa de una partición', 'Relative rate of a partition'],
    what: ['En un análisis particionado, todas las particiones comparten la topología y las proporciones de las ramas, pero cada una tiene un multiplicador que dice si va más rápido o más despacio que la media.',
      'In a partitioned analysis every partition shares the topology and the branch proportions, but each has a multiplier saying whether it runs faster or slower than the average.'],
    read: ['Los multiplicadores se normalizan a media uno pesada por sitios. Una partición con 0.3 evoluciona a un tercio de la velocidad media; con 2.5, dos veces y media más rápido.',
      'The multipliers are normalised to a site-weighted mean of one. A partition at 0.3 evolves at a third of the average speed; at 2.5, two and a half times faster.'],
    care: ['Que cada partición tenga su velocidad no la libera de compartir el árbol. Si dos genes tienen historias distintas de verdad, el análisis particionado los promedia en vez de destaparlo: eso se ve en el Bloque 11.',
      'Giving each partition its own speed does not free it from sharing the tree. If two genes really have different histories, the partitioned analysis averages them instead of revealing it: that is what Block 11 is for.'],
  });

  E('marginal', {
    t: ['Verosimilitud marginal', 'Marginal likelihood'],
    what: ['La probabilidad de los datos bajo un modelo, <b>promediada sobre todos los valores posibles de sus parámetros</b> según el prior. No es el máximo, es un promedio, y por eso castiga sola la complejidad inútil.',
      'The probability of the data under a model, <b>averaged over all the possible values of its parameters</b> according to the prior. It is not the maximum, it is an average, which is why it penalises useless complexity by itself.'],
    read: ['Se estima con piedras de paso o con muestreo por caminos, nunca con la media armónica, que es notoriamente mala. Solo se usa para comparar, a través del factor de Bayes.',
      'It is estimated by stepping stones or path sampling, never by the harmonic mean, which is notoriously bad. It is only used for comparing, through the Bayes factor.'],
    care: ['Su valor depende mucho del prior, así que dos corridas con priors distintos no se pueden comparar. Y la estimación tiene su propio error: repetirla con otra semilla dice cuánto.',
      'Its value depends strongly on the prior, so two runs with different priors cannot be compared. And the estimate has its own error: repeating it with another seed says how much.'],
    ref: 'Xie et al. 2011; Baele et al. 2012',
  });

  E('rate', {
    t: ['Tasa de sustitución', 'Substitution rate'],
    what: ['Sustituciones por sitio y por unidad de tiempo. Es el factor que convierte longitudes de rama en años, y sale de la regresión raíz–punta, de una calibración, o de un valor publicado para el mismo marcador.',
      'Substitutions per site per unit of time. It is the factor that turns branch lengths into years, and comes from the root-to-tip regression, from a calibration, or from a published value for the same marker.'],
    read: ['Los órdenes de magnitud son muy distintos según el marcador: un virus de ARN ronda 10⁻³ por sitio y año; el rbcL de plantas, del orden de 10⁻⁹ a 10⁻¹⁰. Un resultado fuera del orden esperado es casi siempre un error de unidades.',
      'The orders of magnitude differ wildly by marker: an RNA virus is around 10⁻³ per site per year; plant rbcL, of the order of 10⁻⁹ to 10⁻¹⁰. A result outside the expected order is almost always a units mistake.'],
    care: ['Una tasa tomada de otro estudio arrastra el árbol, el modelo y las calibraciones de aquel estudio. Es la fuente de error más grande de una datación y la que menos se reporta.',
      'A rate taken from another study drags in that study\'s tree, model and calibrations. It is the largest source of error in a dating analysis and the least often reported.'],
  });

  E('ratevar', {
    t: ['Variación de la tasa entre ramas', 'Rate variation among branches'],
    what: ['Cuánto se aparta el árbol de un reloj estricto. En un reloj relajado lognormal se resume con el <b>coeficiente de variación</b> de las tasas de las ramas.',
      'How far the tree departs from a strict clock. Under a lognormal relaxed clock it is summarised by the <b>coefficient of variation</b> of the branch rates.'],
    read: ['Un coeficiente cercano a 0 dice que había reloj y que un reloj estricto habría bastado. Por encima de 0.5 las tasas son muy desiguales y el reloj estricto habría dado fechas equivocadas.',
      'A coefficient near 0 says there was a clock and a strict one would have sufficed. Above 0.5 the rates are very uneven and a strict clock would have given wrong dates.'],
    scaleTitle: ['Coeficiente de variación', 'Coefficient of variation'],
    scale: [
      S(0, 0.1, ['prácticamente reloj', 'practically clock-like'], 'good'),
      S(0.1, 0.3, ['moderada', 'moderate'], 'ok'),
      S(0.3, 0.6, ['alta', 'high'], 'warn'),
      S(0.6, null, ['muy alta', 'very high'], 'bad'),
    ],
    conv: true,
    care: ['Si el intervalo del coeficiente incluye el cero, los datos no distinguen el reloj relajado del estricto y conviene reportar los dos.',
      'If the interval of the coefficient includes zero, the data cannot tell the relaxed clock from the strict one, and both are worth reporting.'],
    ref: 'Drummond et al. 2006',
  });

  E('crown', {
    t: ['Edad de corona y edad de tallo', 'Crown age and stem age'],
    what: ['La <b>corona</b> es el ancestro común más reciente de las especies vivas del grupo. El <b>tallo</b> es la divergencia del grupo respecto de su hermano, y siempre es más viejo. Entre los dos hay todo el tiempo en que el linaje existió sin diversificarse o con linajes hoy extintos.',
      'The <b>crown</b> is the most recent common ancestor of the living species of the group. The <b>stem</b> is the group\'s divergence from its sister, and is always older. Between the two lies all the time the lineage existed without diversifying, or with lineages now extinct.'],
    read: ['Las tasas de diversificación se calculan con una u otra y <b>dan resultados distintos</b>: con el tallo la tasa sale menor, porque se reparte el mismo número de especies en más tiempo. Hay que decir cuál se usó.',
      'Diversification rates are computed from one or the other and <b>give different results</b>: with the stem the rate comes out lower, because the same number of species is spread over more time. Which one was used has to be stated.'],
    care: ['Confundirlas es un error frecuente al comparar con la literatura: una «edad del grupo» sin apellido no se puede reutilizar.',
      'Confusing them is a frequent mistake when comparing with the literature: an "age of the group" without a qualifier cannot be reused.'],
    ref: 'Magallón & Sanderson 2001',
  });

  E('mccr', {
    t: ['Prueba MCCR', 'MCCR test'],
    what: ['La prueba de tasa constante de Monte Carlo. Simula árboles de tasa constante <b>con el mismo muestreo incompleto que los datos</b>, mide su γ, y compara el γ observado contra esa distribución en vez de contra la normal.',
      'The Monte Carlo constant-rates test. It simulates constant-rate trees <b>with the same incomplete sampling as the data</b>, measures their γ, and compares the observed γ against that distribution instead of against the normal.'],
    read: ['Es la forma correcta de usar γ cuando no se muestreó todo el grupo. Su valor crítico al 5 % es <b>más negativo</b> que −1.645, y cuánto más lo sea dice cuánta desaceleración explicaba el muestreo por sí solo.',
      'It is the right way to use γ when the whole group was not sampled. Its 5 % critical value is <b>more negative</b> than −1.645, and how much more says how much of the slowdown the sampling alone explained.'],
    care: ['Hay que darle el número real de especies del grupo, no el del árbol. Si esa cifra está mal, la prueba corrige de menos o de más.',
      'It has to be given the real number of species in the group, not the number in the tree. If that figure is wrong, the test corrects too little or too much.'],
    ref: 'Pybus & Harvey 2000',
  });

  E('doubling', {
    t: ['Tiempo de duplicación', 'Doubling time'],
    what: ['ln(2)/r: cuánto tarda el grupo en duplicar su número de linajes a la tasa neta estimada. Es la misma información que r, dicha en unidades que se entienden.',
      'ln(2)/r: how long the group takes to double its number of lineages at the estimated net rate. It is the same information as r, said in units one can picture.'],
    read: ['Un tiempo de duplicación de 5 Ma en un grupo de 60 Ma significa que hubo margen de sobra para la diversidad actual; uno de 50 Ma en el mismo grupo dice que la diversidad no se explica por la tasa sola.',
      'A doubling time of 5 Ma in a 60 Ma old group means there was ample room for today\'s diversity; one of 50 Ma in the same group says the diversity is not explained by the rate alone.'],
    care: ['Hereda toda la incertidumbre de r, que a su vez depende de la edad usada y de si es de corona o de tallo.',
      'It inherits all the uncertainty of r, which in turn depends on the age used and on whether it is crown or stem.'],
  });

  E('contrasts', {
    t: ['Contrastes independientes', 'Independent contrasts'],
    what: ['Convierten los valores de n especies en n−1 diferencias entre nodos hermanas, cada una escalada por la longitud de sus ramas. Esas diferencias <b>sí</b> son independientes, que es lo que las especies no eran.',
      'They turn the values of n species into n−1 differences between sister nodes, each scaled by its branch lengths. Those differences <b>are</b> independent, which is what the species were not.'],
    read: ['La regresión de contrastes se ajusta <b>forzada por el origen</b>: el signo de un contraste es arbitrario, así que una ordenada al origen no significaría nada. Su pendiente es idéntica a la de PGLS con λ = 1.',
      'The contrast regression is fitted <b>through the origin</b>: the sign of a contrast is arbitrary, so an intercept would mean nothing. Its slope is identical to that of PGLS with λ = 1.'],
    care: ['Suponen movimiento browniano y longitudes de rama correctas. Si los contrastes estandarizados muestran tendencia contra su desviación esperada, el supuesto falla y hay que transformar las ramas.',
      'They assume Brownian motion and correct branch lengths. If the standardised contrasts show a trend against their expected deviation, the assumption fails and the branches need transforming.'],
    ref: 'Felsenstein 1985; Garland et al. 1992',
  });

  E('simmap', {
    t: ['Mapeo estocástico', 'Stochastic mapping'],
    what: ['En vez de dar una probabilidad por nodo, sortea <b>historias completas</b> del carácter sobre el árbol —dónde y cuándo ocurrió cada cambio— coherentes con las puntas y con el modelo. Cientos de historias resumen la incertidumbre.',
      'Instead of giving a probability per node, it draws <b>complete histories</b> of the character over the tree — where and when each change happened — consistent with the tips and the model. Hundreds of histories summarise the uncertainty.'],
    read: ['Lo que se reporta es el <b>promedio sobre las historias</b>: número esperado de cambios, tiempo pasado en cada estado, y con qué frecuencia cada rama contiene un cambio. Es más informativo que un solo estado por nodo.',
      'What is reported is the <b>average over the histories</b>: the expected number of changes, the time spent in each state, and how often each branch contains a change. It is more informative than a single state per node.'],
    care: ['Las historias son tan buenas como el modelo Mk que las genera: si las tasas están mal estimadas, todas las historias lo estarán en la misma dirección y el promedio no lo revela.',
      'The histories are only as good as the Mk model that generates them: if the rates are badly estimated, every history will be wrong in the same direction and the average will not reveal it.'],
    ref: 'Huelsenbeck, Nielsen & Bollback 2003',
  });

  E('de', {
    t: ['d y e: dispersión y extinción local', 'd and e: dispersal and local extinction'],
    what: ['Los dos parámetros que gobiernan el cambio de rango <b>a lo largo</b> de una rama: d es la tasa a la que un linaje añade un área, e la tasa a la que pierde una. El j, cuando está, actúa solo <b>en el momento</b> de una división.',
      'The two parameters governing range change <b>along</b> a branch: d is the rate at which a lineage adds an area, e the rate at which it loses one. The j, when present, acts only <b>at the moment</b> of a split.'],
    read: ['Se leen en relación con la profundidad del árbol: d·T dice cuántas áreas gana un linaje en toda la historia del grupo. Una e mucho mayor que d produce rangos estrechos; al revés, rangos anchos.',
      'They are read relative to the depth of the tree: d·T says how many areas a lineage gains over the whole history of the group. An e much larger than d yields narrow ranges; the reverse, wide ones.'],
    care: ['Una e estimada en cero es frecuente y no significa que nunca hubo extinción local: el árbol de los rangos actuales suele no poder verla, igual que pasa con μ en diversificación.',
      'An e estimated at zero is frequent and does not mean local extinction never happened: the tree of present-day ranges usually cannot see it, just as with μ in diversification.'],
    ref: 'Ree & Smith 2008',
  });

  E('bsm', {
    t: ['Mapeo estocástico biogeográfico', 'Biogeographic stochastic mapping'],
    what: ['El equivalente del mapeo estocástico para los rangos: sortea historias completas y cuenta <b>cuántos eventos de cada tipo</b> hubo — dispersiones, extinciones locales, vicarianzas, saltos fundadores.',
      'The equivalent of stochastic mapping for ranges: it draws complete histories and counts <b>how many events of each kind</b> there were — dispersals, local extinctions, vicariance events, founder jumps.'],
    read: ['Los conteos se dan como media sobre las historias, con su dispersión. Esa dispersión es el resultado tanto como la media: 12 ± 9 dispersiones no es lo mismo que 12 ± 1.',
      'The counts are given as an average over the histories, with their spread. That spread is as much the result as the mean: 12 ± 9 dispersals is not the same as 12 ± 1.'],
    care: ['Los conteos dependen del modelo que los generó. Con +J, los saltos fundadores aparecerán por construcción; contarlos no demuestra que ocurrieran.',
      'The counts depend on the model that generated them. With +J, founder jumps will appear by construction; counting them does not prove they happened.'],
    ref: 'Dupin et al. 2017',
  });

  E('splits', {
    t: ['Divisiones (bipartitions)', 'Splits (bipartitions)'],
    what: ['Cada rama interna de un árbol sin raíz parte las puntas en dos conjuntos. Ese par es una <b>división</b>, y un árbol no es más que un conjunto de divisiones compatibles entre sí.',
      'Every internal branch of an unrooted tree splits the tips into two sets. That pair is a <b>split</b>, and a tree is nothing more than a set of mutually compatible splits.'],
    read: ['Casi todo lo que compara árboles se hace en divisiones: la RF las cuenta, el bootstrap mide su frecuencia, el consenso se arma con las que pasan un umbral. Un árbol de n puntas sin raíz tiene n−3 divisiones internas.',
      'Almost everything that compares trees works on splits: RF counts them, the bootstrap measures their frequency, the consensus is built from those passing a threshold. An unrooted tree of n tips has n−3 internal splits.'],
    care: ['Dos divisiones son <b>incompatibles</b> cuando no pueden estar en el mismo árbol. Un conjunto de divisiones con incompatibilidades no es un árbol: es una red, y dibujarlo como árbol esconde el conflicto.',
      'Two splits are <b>incompatible</b> when they cannot be on the same tree. A set of splits with incompatibilities is not a tree: it is a network, and drawing it as a tree hides the conflict.'],
  });

  E('quartet', {
    t: ['Cuartetos', 'Quartets'],
    what: ['Cuatro puntas cualesquiera se pueden resolver de tres maneras. Un árbol grande queda determinado por cómo resuelve todos sus cuartetos, y contarlos es una forma de comparar árboles que no depende de la raíz.',
      'Any four tips can be resolved in three ways. A large tree is determined by how it resolves all of its quartets, and counting them is a way of comparing trees that does not depend on the root.'],
    read: ['La distancia de cuartetos es la fracción resuelta de forma distinta. A diferencia de la RF, <b>no salta al máximo por mover una sola punta</b>: es una medida gradual, y por eso suele reflejar mejor cuánto discrepan dos árboles.',
      'The quartet distance is the fraction resolved differently. Unlike RF, it <b>does not jump to its maximum when a single tip moves</b>: it is a gradual measure, and usually reflects better how much two trees disagree.'],
    care: ['El número de cuartetos crece como n⁴, así que sobre árboles grandes se estima por muestreo. La búsqueda de árbol de especies por cuartetos de esta app es una heurística, no el algoritmo exacto de ASTRAL.',
      'The number of quartets grows as n⁴, so on large trees it is estimated by sampling. This app\'s quartet species-tree search is a heuristic, not ASTRAL\'s exact algorithm.'],
  });

  E('dcf', {
    t: ['Factores de discordancia (gDF1, gDF2, sDF1, sDF2)', 'Discordance factors (gDF1, gDF2, sDF1, sDF2)'],
    what: ['Lo que no está de acuerdo, repartido entre las <b>dos</b> resoluciones alternativas de la rama. gDF para genes, sDF para sitios.',
      'What disagrees, split between the <b>two</b> alternative resolutions of the branch. gDF for genes, sDF for sites.'],
    read: ['Ahí está la pregunta interesante. Si las dos alternativas se reparten el desacuerdo <b>por partes iguales</b>, eso es lo que la clasificación incompleta de linajes predice. Si una alternativa domina claramente a la otra, hay algo más: flujo génico, hibridación o un error sistemático.',
      'This is where the interesting question lies. If the two alternatives share the disagreement <b>equally</b>, that is what incomplete lineage sorting predicts. If one alternative clearly dominates the other, something else is going on: gene flow, hybridisation or a systematic error.'],
    care: ['El desequilibrio entre gDF1 y gDF2 es una <b>pista</b>, no una prueba. La prueba es el estadístico D del mismo bloque.',
      'The imbalance between gDF1 and gDF2 is a <b>hint</b>, not a test. The test is the D statistic in the same block.'],
    ref: 'Minh, Hahn & Lanfear 2020',
  });

  E('tanglegram', {
    t: ['Tanglegrama y cruces', 'Tanglegram and crossings'],
    what: ['Dos árboles enfrentados con las mismas puntas unidas por líneas. Cada cruce de líneas es un desacuerdo en el orden de las puntas.',
      'Two trees face to face with the same tips joined by lines. Each crossing of lines is a disagreement in the order of the tips.'],
    read: ['El número de cruces depende del <b>orden en que se dibujan las puntas</b>, no solo de los árboles: por eso se desenreda primero, girando los nodos, que no cambia ninguno de los dos árboles. Lo que queda después de desenredar es el conflicto real.',
      'The number of crossings depends on the <b>order in which the tips are drawn</b>, not only on the trees: which is why it is untangled first, by rotating nodes, which changes neither tree. What remains after untangling is the real conflict.'],
    care: ['Un tanglegram impresiona pero <b>no es una medida</b>: el mismo par de árboles puede verse enmarañado o limpio según el dibujo. Para cuantificar están las cinco distancias del bloque.',
      'A tanglegram is striking but <b>is not a measurement</b>: the same pair of trees can look tangled or clean depending on the drawing. The five distances of the block are what quantifies it.'],
  });

  E('network', {
    t: ['Redes de divisiones', 'Split networks'],
    what: ['Cuando las divisiones que apoyan los datos son incompatibles entre sí, no hay un árbol que las contenga a todas. Una red las dibuja todas: cada división es un haz de líneas paralelas, y los rectángulos son el conflicto.',
      'When the splits the data support are mutually incompatible, no tree contains them all. A network draws them all: each split is a bundle of parallel lines, and the boxes are the conflict.'],
    read: ['Ancho del rectángulo = fuerza del conflicto. Una red casi arbórea dice que los datos son coherentes; una llena de cajas dice que empeñarse en un árbol es forzar los datos.',
      'The width of a box is the strength of the conflict. A nearly tree-like network says the data are coherent; one full of boxes says insisting on a tree forces the data.'],
    care: ['Una red <b>no es una filogenia con hibridación</b>: las cajas pueden venir de hibridación, de clasificación incompleta de linajes, de recombinación o simplemente de ruido. La red muestra el conflicto, no lo explica. Y la descomposición circular solo puede dibujar las divisiones que caben en el orden circular elegido; las demás se reportan aparte.',
      'A network <b>is not a phylogeny with hybridisation</b>: the boxes can come from hybridisation, from incomplete lineage sorting, from recombination or simply from noise. The network shows the conflict, it does not explain it. And the circular decomposition can only draw the splits that fit the chosen circular order; the rest are reported separately.'],
    ref: 'Bandelt & Dress 1992; Huson & Bryant 2006',
  });

  E('reproducible', {
    t: ['Informe y paquete reproducible', 'Reproducible report and package'],
    what: ['Los métodos redactados <b>a partir de lo que realmente se corrió</b> —no de una plantilla— con sus números y sus referencias, más un paquete con los datos, los árboles, las figuras y un LEEME que explica cada archivo.',
      'The methods written <b>from what was actually run</b> — not from a template — with their numbers and their references, plus a package with the data, the trees, the figures and a README explaining every file.'],
    read: ['Si un bloque no se corrió, no deja párrafo: el informe no menciona reloj, biogeografía ni MCMC si no los hubo. Las réplicas y los parámetros que aparecen son los de tu corrida.',
      'If a block was not run, it leaves no paragraph: the report mentions no clock, no biogeography and no MCMC if there were none. The replicates and parameters that appear are those of your run.'],
    care: ['El informe describe lo que el programa hizo, no justifica que estuviera bien hecho. Las decisiones —el alineamiento, las calibraciones, las áreas— siguen siendo tuyas y hay que defenderlas.',
      'The report describes what the program did, it does not justify that it was the right thing to do. The decisions — the alignment, the calibrations, the areas — remain yours and have to be defended.'],
  });

  /* =====================================================================
     2 · reading a scale, and applying it to a number
     ===================================================================== */

  function bandRange(b) {
    if (b.txt) return esc(b.txt);
    if (b.from == null && b.to == null) return '';
    if (b.from == null) return '&lt; ' + num(b.to);
    if (b.to == null) return '&gt; ' + num(b.from);
    if (b.from === b.to) return '= ' + num(b.from);
    return num(b.from) + ' – ' + num(b.to);
  }

  /* Widths follow each band's span, open ends get an average share, and nothing
     is allowed to become too thin to read. */
  function bandWidths(scale) {
    const spans = scale.map(b => (b.from != null && b.to != null && b.to > b.from) ? (b.to - b.from) : null);
    const finite = spans.filter(v => v != null);
    const mean = finite.length ? finite.reduce((s, v) => s + v, 0) / finite.length : 1;
    const raw = spans.map(v => (v == null || v === 0) ? mean : v);
    const m = raw.reduce((s, v) => s + v, 0) / raw.length || 1;
    return raw.map(v => Math.max(0.55, Math.min(2.2, v / m)));
  }

  /* Which band a computed value falls into. This is the part that turns a
     scale from documentation into a decision rule: the reader does not have to
     carry the thresholds in their head and compare by eye. */
  function band(key, value) {
    const d = HELP[key];
    if (!d || !d.scale || value == null || !isFinite(value)) return null;
    for (const b of d.scale) {
      const okLo = b.from == null || value >= b.from;
      const okHi = b.to == null || value < b.to;
      if (okLo && okHi) return b;
    }
    /* a value exactly at the top of the last closed band belongs to it */
    const last = d.scale[d.scale.length - 1];
    return (last && last.to != null && value === last.to) ? last : null;
  }

  /* the band as a small chip, ready to print beside the number itself */
  function tag(key, value, opts) {
    const b = band(key, value);
    if (!b) return '';
    const cls = 'help-tag help-t-' + (b.tone || 'neutral');
    const ttl = (opts && opts.title === false) ? ''
      : ` title="${esc(T('Escala de ', 'Scale of ') + (HELP[key].t[I18N.lang === 'en' ? 1 : 0]))}"`;
    return `<span class="${cls}"${ttl}>${two(b.label)}</span>`;
  }

  /* the plain words, for a sentence written by a block */
  function verdict(key, value) {
    const b = band(key, value);
    return b ? T(b.label[0], b.label[1]) : '';
  }

  function scaleHTML(d) {
    if (!d.scale || !d.scale.length) return '';
    const w = bandWidths(d.scale);
    const title = d.scaleTitle ? `<div class="help-lead">${two(d.scaleTitle)}</div>` : '';
    const bands = d.scale.map((b, i) =>
      `<div class="help-band help-t-${b.tone || 'neutral'}" style="flex:${w[i].toFixed(2)} 1 0">` +
      `<b>${bandRange(b)}</b><span>${two(b.label)}</span></div>`).join('');
    return `<div class="help-scale">${title}<div class="help-bands">${bands}</div></div>`;
  }

  const CONV = ['Esta escala es una convención de lectura de uso común, no una ley estadística.',
    'This scale is a widely used reading convention, not a statistical law.'];

  /* the body of one entry, used both inside a popover and inside the guide */
  function entryBody(d) {
    return `<p>${two(d.what)}</p>` +
      (d.formula ? `<div class="help-formula">${esc(d.formula)}</div>` : '') +
      `<div class="help-lead">${L2('Cómo se lee', 'How to read it')}</div><p>${two(d.read)}</p>` +
      scaleHTML(d) +
      (d.conv ? `<p class="help-conv">${two(CONV)}</p>` : '') +
      (d.care ? `<div class="help-lead">${L2('Cuidado con', 'Watch out for')}</div>` +
        `<p class="help-care">${two(d.care)}</p>` : '') +
      (d.ref ? `<p class="help-ref">${L2('Referencia:', 'Reference:')} ${esc(d.ref)}</p>` : '');
  }

  /* =====================================================================
     3 · the badge and its popover
     ===================================================================== */

  const badgeTitle = d => T('Qué significa: ' + d.t[0], 'What it means: ' + d.t[1]);

  function badge(key) {
    const d = HELP[key];
    if (!d) { console.warn('Help: unknown key', key); return ''; }
    const ttl = esc(badgeTitle(d));
    return `<button type="button" class="help-badge" data-help-key="${esc(key)}" aria-expanded="false" ` +
      `aria-haspopup="dialog" title="${ttl}" aria-label="${ttl}">?</button>`;
  }

  let pop = null, openBtn = null, posTimer = null;

  function ensurePop() {
    if (pop) return pop;
    pop = document.createElement('div');
    pop.className = 'help-pop';
    pop.setAttribute('role', 'dialog');
    pop.setAttribute('tabindex', '-1');
    pop.style.display = 'none';
    document.body.appendChild(pop);
    return pop;
  }

  function fillPop(key) {
    const d = HELP[key];
    ensurePop().innerHTML =
      `<button type="button" class="help-close" data-help-close="1" ` +
      `title="${esc(T('Cerrar', 'Close'))}" aria-label="${esc(T('Cerrar', 'Close'))}">✕</button>` +
      `<h4>${two(d.t)}</h4>` + entryBody(d);
  }

  /* placed beside the badge when there is room, so it never covers the number
     it is explaining */
  function placePop() {
    if (!pop || !openBtn || !openBtn.isConnected) return;
    const r = openBtn.getBoundingClientRect(), vw = innerWidth, vh = innerHeight;
    if (vw <= 560) { pop.style.top = ''; pop.style.left = ''; return; }  // the stylesheet makes it a sheet
    pop.style.top = '0px'; pop.style.left = '0px';
    const w = pop.offsetWidth, h = pop.offsetHeight;
    let left, top;
    if (r.right + 10 + w <= vw - 8) left = r.right + 10;
    else if (r.left - 10 - w >= 8) left = r.left - 10 - w;
    else left = Math.max(8, Math.min(vw - w - 8, r.left - w / 2));
    if (left === r.right + 10 || left === r.left - 10 - w) top = r.top - 6;
    else top = (r.bottom + 8 + h <= vh - 8) ? r.bottom + 8 : r.top - 8 - h;
    top = Math.max(8, Math.min(vh - h - 8, top));
    pop.style.left = left + 'px';
    pop.style.top = top + 'px';
  }
  const schedulePlace = () => { cancelAnimationFrame(posTimer); posTimer = requestAnimationFrame(placePop); };

  function openPop(btn) {
    const key = btn.dataset.helpKey;
    if (!HELP[key]) return;
    if (openBtn === btn) { closePop(true); return; }
    closePop(false);
    fillPop(key);
    openBtn = btn;
    btn.setAttribute('aria-expanded', 'true');
    pop.style.display = 'block';
    placePop();
    pop.focus({ preventScroll: true });
  }

  function closePop(refocus) {
    if (!pop || pop.style.display === 'none') { openBtn = null; return; }
    pop.style.display = 'none';
    if (openBtn) {
      openBtn.setAttribute('aria-expanded', 'false');
      if (refocus && openBtn.isConnected) openBtn.focus({ preventScroll: true });
    }
    openBtn = null;
  }

  if (typeof document !== 'undefined') {
    document.addEventListener('click', e => {
      const close = e.target.closest && e.target.closest('[data-help-close]');
      if (close) { e.preventDefault(); closePop(true); return; }
      const b = e.target.closest && e.target.closest('.help-badge');
      if (b) { e.preventDefault(); e.stopPropagation(); openPop(b); return; }
      if (pop && pop.style.display !== 'none' && !(e.target.closest && e.target.closest('.help-pop'))) closePop(false);
    }, true);
    document.addEventListener('keydown', e => {
      if (e.key === 'Escape' && pop && pop.style.display !== 'none') { e.stopPropagation(); closePop(true); }
    });
    addEventListener('resize', schedulePlace);
    addEventListener('scroll', schedulePlace, true);
    document.addEventListener('langchange', () => { if (openBtn) fillPop(openBtn.dataset.helpKey); });
  }

  /* =====================================================================
     4 · decorating what the blocks already build
     ===================================================================== */

  /* the texts a header can be matched against: the whole text and each language */
  function candidates(node) {
    const out = [(node.textContent || '').trim()];
    node.querySelectorAll('[data-l]').forEach(s => out.push((s.textContent || '').trim()));
    return out.filter(Boolean).map(s => s.toLowerCase());
  }

  function decorate(nodes, map, prefix) {
    const pairs = Object.keys(map).map(k => [k.trim().toLowerCase(), map[k]]);
    nodes.forEach(n => {
      if (n.dataset.helpDone) return;
      const cand = candidates(n);
      const hit = prefix ? pairs.find(([txt]) => cand.some(c => c.startsWith(txt)))
        : pairs.find(([txt]) => cand.includes(txt));
      n.dataset.helpDone = '1';
      if (!hit || !HELP[hit[1]]) return;
      n.insertAdjacentHTML('beforeend', badge(hit[1]));
    });
  }

  function markTable(container, map, prefix) {
    if (typeof container === 'string') container = el(container);
    if (!container) return;
    decorate([...container.querySelectorAll('thead th')], map, prefix);
  }
  function markTiles(container, map, prefix) {
    if (typeof container === 'string') container = el(container);
    if (!container) return;
    decorate([...container.querySelectorAll('.stat-label')], map, prefix);
  }
  /* every <span data-help="key"> written in the page becomes a badge */
  function hydrate(root) {
    (root || document).querySelectorAll('[data-help]:not([data-help-ready])').forEach(n => {
      const key = n.getAttribute('data-help');
      n.setAttribute('data-help-ready', '1');
      if (HELP[key]) n.innerHTML = badge(key);
    });
  }

  /* =====================================================================
     5 · the collapsible guide of a block
     ===================================================================== */

  const GUIDES = {};

  function panel(keys, title) {
    const id = 'helpGuide' + (panel._n = (panel._n || 0) + 1);
    GUIDES[id] = keys.filter(k => HELP[k]);
    return `<details class="acc help-guide" id="${id}">` +
      `<summary>${title ? two(title) : L2('📖 Guía de interpretación de este bloque', '📖 Interpretation guide for this block')}</summary>` +
      `<div class="acc-body"><p class="hint">${L2(
        'Qué mide cada número, en qué escala se lee y cuál es el error más común con él. Las escalas marcadas como convención son costumbres de lectura, no leyes.',
        'What each number measures, on what scale it is read and the commonest mistake made with it. The scales marked as a convention are reading habits, not laws.')}</p><div class="help-guide-grid"></div></div></details>`;
  }

  /* a guide is written the first time it is opened, never at load */
  function fillGuide(det) {
    const box = det.querySelector('.help-guide-grid');
    if (!box || box.dataset.ready) return;
    box.dataset.ready = '1';
    box.innerHTML = (GUIDES[det.id] || []).map(k =>
      `<section class="help-entry"><h5>${two(HELP[k].t)}</h5>${entryBody(HELP[k])}</section>`).join('');
  }
  if (typeof document !== 'undefined') {
    document.addEventListener('toggle', e => {
      if (e.target.classList && e.target.classList.contains('help-guide') && e.target.open) fillGuide(e.target);
    }, true);
  }

  /* =====================================================================
     6 · which entries belong to which block, and installing itself
     ===================================================================== */

  const BLOCK_KEYS = {
    1: ['longbranch', 'saturation', 'splits'],
    2: ['informative', 'variable', 'gaps', 'trimming', 'sp', 'gc', 'chicomp'],
    3: ['lnl', 'aic', 'akaike', 'lrt', 'alpha', 'pinv', 'titv', 'gc', 'saturation', 'chicomp'],
    4: ['distance', 'steps', 'equaltrees', 'ci', 'ri', 'bremer', 'bootstrap', 'search', 'rf'],
    5: ['lnl', 'search', 'partrate', 'bootstrap', 'ufboot', 'shalrt', 'autest'],
    6: ['pp', 'ess', 'psrf', 'asdsf', 'burnin', 'marginal', 'bayesfactor'],
    7: ['rtt', 'clocktest', 'rate', 'ratevar', 'calibration', 'smoothing', 'hpd', 'crown'],
    8: ['gamma', 'mccr', 'birthdeath', 'doubling', 'crown', 'dr'],
    9: ['pagel', 'blombergk', 'ou', 'contrasts', 'pgls', 'mk', 'ancstate', 'simmap'],
    10: ['decj', 'de', 'ancrange', 'bsm'],
    11: ['splits', 'rf', 'quartet', 'gcf', 'scf', 'dcf', 'dstat', 'tanglegram', 'network'],
    12: ['dpi', 'reproducible'],
  };

  /* Labels the app already prints, in either language, and the entry each one
     belongs to. Matching is by prefix and case-insensitive, so "Bootstrap (%)"
     and "bootstrap" both find the same entry. */
  const LABELS = {
    'lnl': 'lnl', 'log-verosimilitud': 'lnl', 'log-likelihood': 'lnl', 'verosimilitud': 'lnl',
    'aic': 'aic', 'aicc': 'aic', 'bic': 'aic', 'δaic': 'aic', 'Δaic': 'aic',
    'peso': 'akaike', 'weight': 'akaike', 'peso de akaike': 'akaike', 'akaike weight': 'akaike',
    'α': 'alpha', 'alpha': 'alpha', 'forma de la gamma': 'alpha', 'gamma shape': 'alpha',
    'i': 'pinv', 'invariables': 'pinv', 'invariable': 'pinv',
    'sitios informativos': 'informative', 'informative sites': 'informative', 'informativos': 'informative',
    'huecos': 'gaps', 'gaps': 'gaps',
    'pasos': 'steps', 'steps': 'steps', 'longitud': 'steps', 'length': 'steps',
    'ci': 'ci', 'ri': 'ri', 'bremer': 'bremer',
    'bootstrap': 'bootstrap', 'jackknife': 'bootstrap',
    'rf': 'rf', 'robinson': 'rf',
    'ufboot': 'ufboot', 'sh-alrt': 'shalrt', 'abayes': 'shalrt',
    'au': 'autest', 'kh': 'autest', 'sh': 'autest',
    'pp': 'pp', 'probabilidad posterior': 'pp', 'posterior': 'pp',
    'ess': 'ess', 'psrf': 'psrf', 'asdsf': 'asdsf',
    'factor de bayes': 'bayesfactor', 'bayes factor': 'bayesfactor',
    'r²': 'rtt', 'r2': 'rtt',
    'hpd': 'hpd', 'λ (suavizado)': 'smoothing',
    'γ': 'gamma', 'gamma': 'gamma',
    'ε': 'birthdeath', 'recambio': 'birthdeath', 'turnover': 'birthdeath',
    'λ': 'birthdeath', 'μ': 'birthdeath', 'r (neta)': 'birthdeath',
    'dr': 'dr',
    'λ de pagel': 'pagel', "pagel's λ": 'pagel', 'lambda': 'pagel',
    'k de blomberg': 'blombergk', "blomberg's k": 'blombergk', 'k': 'blombergk',
    'vida media': 'ou', 'half-life': 'ou', 'α (ou)': 'ou',
    'gcf': 'gcf', 'scf': 'scf', 'd de patterson': 'dstat', "patterson's d": 'dstat',
    'ppp': 'dpi', 'dpi': 'dpi', 'resolución': 'dpi', 'resolution': 'dpi',

    /* the second pass: labels the panels print that had no entry before */
    'columnas conservadas': 'trimming', 'columnas eliminadas': 'trimming',
    'columns kept': 'trimming', 'columns removed': 'trimming', 'bloques': 'trimming',
    'variables': 'variable', 'constantes': 'variable', 'constant': 'variable', 'variable': 'variable',
    'sitios': 'variable', 'sites': 'variable',
    'contenido g+c': 'gc', 'g+c content': 'gc', 'g+c': 'gc',
    'razón ts/tv': 'titv', 'ts/tv': 'titv', 'κ': 'titv', 'kappa': 'titv',
    '2δlnl': 'lrt', '2dlnl': 'lrt', 'df': 'lrt', 'prueba': 'lrt', 'decisión': 'lrt', 'decision': 'lrt',
    'forma de γ': 'alpha', 'forma de la gamma': 'alpha', 'shape of γ': 'alpha',
    'sitios invariables': 'pinv', 'invariable sites': 'pinv',
    'distancia media': 'distance', 'distancia máxima': 'distance', 'mean distance': 'distance',
    'pares sin corrección': 'distance', 'uncorrected pairs': 'distance',
    'árboles igual de cortos': 'equaltrees', 'equally short trees': 'equaltrees',
    'reordenamientos': 'search', 'rearrangements': 'search', 'rondas': 'search', 'rounds': 'search',
    'velocidad relativa': 'partrate', 'relative rate': 'partrate',
    'log p(datos': 'marginal', 'log p(data': 'marginal', 'verosimilitud marginal': 'marginal',
    'marginal likelihood': 'marginal',
    '2 ln bf': 'bayesfactor', 'evidencia': 'bayesfactor', 'evidence': 'bayesfactor',
    'tasa': 'rate', 'tasa media': 'rate', 'rate': 'rate', 'mean rate': 'rate',
    'variación de la tasa': 'ratevar', 'rate variation': 'ratevar',
    'edad de la corona': 'crown', 'crown age': 'crown', 'edad del tallo': 'crown', 'stem age': 'crown',
    'edad de la raíz': 'crown', 'root age': 'crown',
    'valor crítico al 5': 'mccr', 'critical value': 'mccr', 'réplicas': 'mccr', 'replicates': 'mccr',
    'tiempo de duplicación': 'doubling', 'doubling time': 'doubling',
    'tasa neta de diversificación': 'birthdeath', 'net diversification': 'birthdeath',
    'tasa de nacimiento puro': 'birthdeath', 'pure birth': 'birthdeath',
    'por el origen': 'contrasts', 'through the origin': 'contrasts', 'contrastes': 'contrasts',
    'contrasts': 'contrasts', 'pendiente': 'pgls', 'slope': 'pgls',
    'p de la permutación': 'blombergk', 'permutation p': 'blombergk',
    'p de λ': 'pagel', 'cambios por historia': 'simmap', 'changes per history': 'simmap',
    'estado en la raíz': 'ancstate', 'root state': 'ancstate', 'historias sorteadas': 'simmap',
    'd': 'de', 'e': 'de', 'j': 'decj',
    'dispersiones': 'bsm', 'dispersals': 'bsm', 'extinciones locales': 'bsm',
    'local extinctions': 'bsm', 'saltos fundadores': 'bsm', 'founder jumps': 'bsm',
    'rango en la raíz': 'ancrange', 'rango más probable': 'ancrange', 'nodos por debajo': 'ancrange',
    'divisiones': 'splits', 'splits': 'splits', 'división': 'splits', 'split': 'splits',
    'cuartetos': 'quartet', 'quartets': 'quartet',
    'gdf1': 'dcf', 'gdf': 'dcf', 'sdf1': 'dcf', 'sdf': 'dcf',
    'cruces': 'tanglegram', 'crossings': 'tanglegram', 'puntas que se mueven': 'tanglegram',
    'descartadas por peso negativo': 'network', 'no dibujables': 'network', 'red': 'network',
    'abba': 'dstat', 'baba': 'dstat', 'z': 'dstat', 'tríos probados': 'dstat',
    'referencias citadas': 'reproducible', 'secciones de métodos': 'reproducible',
    'tamaño del paquete': 'reproducible', 'figuras incluidas': 'reproducible',
    'tamaño del dibujo': 'dpi',

    /* a third pass: labels that needed no new entry, only pointing at one.
       Deliberately absent: the bare "p", "t", "z", "D" and "Δ", which appear in
       several different tests in this app — sending all of them to one entry
       would be worse than leaving them alone. */
    'ultrarrápido': 'ufboot', 'ultrafast': 'ufboot',
    'δ lnl': 'lnl', 'd lnl': 'lnl', 'Δ lnl': 'lnl',
    'mejor modelo': 'aic', 'modelo elegido': 'aic', 'best model': 'aic', 'chosen model': 'aic',
    'parámetros vigilados': 'ess', 'watched parameters': 'ess',
    '2.5 %': 'hpd', 'mediana – 97.5': 'hpd',
    'ramas internas': 'splits', 'internal branches': 'splits',
    'ramas por debajo del azar': 'scf', 'branches below chance': 'scf',
    'razón máximo/mínimo': 'dr', 'max/min ratio': 'dr',
    'nodos reconstruidos': 'ancrange', 'puntas con rango': 'ancrange',
    'rango observado más ancho': 'ancrange', 'quiénes': 'ancrange', 'su probabilidad': 'ancrange',
    'qué permite al dividirse': 'decj', 'what it allows at a split': 'decj',
    'árboles comparados': 'quartet', 'mejor de los árboles dados': 'quartet',
    'votos recogidos': 'quartet', 'cuartetos a favor': 'quartet', 'cuartetos discrepantes': 'quartet',
    'más parecidos': 'rf', 'más distintos': 'rf', 'coincidencia': 'rf', 'máximo posible': 'rf',
    'pgls': 'pgls', 'ordinaria': 'pgls', 'error estándar': 'pgls', 'error est.': 'dstat',
    'un solo modelo': 'partrate', 'p1, p2': 'dstat', 'p3': 'dstat',
    'esperados por azar': 'dstat', 'con |z|': 'dstat',
  };

  /* Rather than asking twelve panels to remember to decorate themselves, the
     two functions every block already uses to print a table or a row of tiles
     are wrapped once. A label the app prints anywhere gets its badge. */
  function wrapBuilders() {
    if (typeof window === 'undefined') return;
    ['statTiles', 'buildTable'].forEach(name => {
      const orig = window[name];
      if (typeof orig !== 'function' || orig.__helped) return;
      const wrapped = function () {
        const out = orig.apply(this, arguments);
        try {
          const c = typeof arguments[0] === 'string' ? el(arguments[0]) : arguments[0];
          if (c) (name === 'statTiles' ? markTiles : markTable)(c, LABELS, true);
        } catch (e) { /* the help must never break what it decorates */ }
        return out;
      };
      wrapped.__helped = true;
      window[name] = wrapped;
    });
  }

  /* The guide is appended to every block panel on its own, so no panel has to
     remember to ask for it and none can be forgotten. */
  function install(root) {
    const scope = root || document;
    Object.keys(BLOCK_KEYS).forEach(n => {
      const p = scope.querySelector ? scope.querySelector('#panel-' + n) : null;
      if (!p || p.querySelector(':scope > .help-guide')) return;
      const box = document.createElement('div');
      box.innerHTML = panel(BLOCK_KEYS[n]);
      const node = box.firstElementChild;
      if (node) p.appendChild(node);
    });
    wrapBuilders();
    hydrate(scope);
  }

  if (typeof document !== 'undefined') {
    if (document.readyState === 'loading') document.addEventListener('DOMContentLoaded', () => install());
    else install();
  }

  Object.assign(Help, {
    HELP, BLOCK_KEYS, LABELS, install, badge, panel, markTable, markTiles, hydrate,
    band, tag, verdict, entryBody, scaleHTML,
    keys: () => Object.keys(HELP),
    _S: S, _E: E, _two: two, _num: num, _bandWidths: bandWidths, _bandRange: bandRange,
  });
  if (typeof window !== 'undefined') window.Help = Help;

})();
