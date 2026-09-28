# Manual de usuario de PhylogenyPro

Manual en español, en HTML, pensado para imprimirse a PDF desde el navegador.
Sigue la misma forma que los de las demás aplicaciones LABG: hoja tamaño carta
sin márgenes añadidos, Cormorant para los títulos y Jost para los rótulos, y una
portada que es **una sola ilustración vectorial original**, dibujada por código
con semilla fija para que salga idéntica cada vez e imprima nítida a cualquier
tamaño.

## Cómo está partido

Un archivo por parte, en `es/`, y cada capítulo explica el bloque del mismo
número — la introducción no lleva número:

| Archivo | Contenido | Estado |
|---|---|---|
| `00a-portada.html` | Portada | ✅ |
| `00b-introduccion.html` | Créditos, índice e introducción (secciones I.1 a I.10) | ✅ |
| `01-bloque1.html` | Capítulo 1 · Bloque 1: la portada, los dos laboratorios y los límites | ✅ |
| `02-bloque2.html` | Capítulo 2 · Bloque 2: datos, alineamiento, recorte y control de calidad | ✅ |
| `03-bloque3.html` | Capítulo 3 · Bloque 3: supuestos, 56 modelos de ADN, matrices de proteína y particiones | ✅ |
| `04-bloque4.html` | Capítulo 4 · Bloque 4: distancias, parsimonia, índices, soporte y consensos | ✅ |
| `05-bloque5.html` | Capítulo 5 · Bloque 5: búsqueda ML, particiones, los cuatro apoyos y las pruebas de topología | ✅ |
| `06-bloque6.html` | Capítulo 6 · Bloque 6: priors, MCMC acoplado, convergencia y factores de Bayes | ✅ |
| `07-bloque7.html` | Capítulo 7 · Bloque 7: raíz, reloj, calibraciones, el asistente y los tres métodos de datación | ✅ |
| `08-bloque8.html` | Capítulo 8 · Bloque 8: linajes en el tiempo, γ, MCCR, cinco modelos y DR | ✅ |
| `09-bloque9.html` | Capítulo 9 · Bloque 9: BM/OU/EB, señal filogenética, contrastes, PGLS, Mk y mapeo estocástico | ✅ |
| `10-bloque10.html` | Capítulo 10 · Bloque 10: áreas, seis modelos DEC/DIVALIKE/BAYAREALIKE ±J, esquinas y mapeo biogeográfico | ✅ |
| `11-bloque11.html` | Capítulo 11 · Bloque 11: distancias, tanglegrama, gCF/sCF, cuartetos, redes y ABBA-BABA | ✅ |
| `12-bloque12.html` | Capítulo 12 · Bloque 12: estudio de árboles y de figuras, tres formatos, resolución, informe y paquete | ✅ |
| `13-apendices.html` | Apéndices A–G: archivos, catálogo de modelos, reglas de decisión, glosario, fallos comunes, referencias y licencia | ✅ |
| `manual-completo.html` | Las partes unidas, que es lo que se imprime — lo genera `herramientas/unir-manual.ps1` | ✅ |

La portada usa sólo `manual.css`, porque es una hoja entera de dibujo. Todo lo
demás usa `interior.css` y `paginar.js`, que reparten el texto en hojas con su
encabezado corrido, su número de página y la pestaña del capítulo en el margen.
Los números de página del índice se llenan solos al abrir el archivo; mientras
falten capítulos, sus renglones quedan en blanco, que es lo esperado hasta que
se unan las partes.

Al partir una tabla, una lista o una lista de definiciones entre dos hojas, `paginar.js`
deja **al menos dos filas o dos incisos de cada lado**; si no puede, mueve el elemento
entero a la hoja siguiente. Así no aparece un inciso solitario al pie de una hoja con
el resto de su lista al dorso. La numeración continúa sola, y una tabla partida repite
su encabezado y añade «(continuación)» a su leyenda.

## Cómo se ve mientras se escribe

Abre cualquier parte directamente en el navegador. Para verlas todas juntas hace
falta un servidor local, porque el navegador no deja que un archivo `file://`
cargue a sus hermanos:

```bash
powershell -ExecutionPolicy Bypass -File server.ps1 -Port 9801
```

y luego `http://localhost:9801/manual/es/00-portada.html`.

## Cómo se imprime

**Sólo al final, nunca tras cada parte.** Se unen las partes en
`manual-completo.html` y se imprime desde el navegador con:

- tamaño de papel **Carta**,
- márgenes **ninguno**,
- **gráficos de fondo activados** (si no, las franjas de tiempo geológico y los
  colores de las figuras salen en blanco).

## La portada

Cuenta el recorrido de la aplicación en una sola imagen, de izquierda a derecha:
un **alineamiento** cuyas hebras se trenzan hacia la raíz de un **árbol fechado**
que asciende por las franjas de la escala de tiempo geológico —con los **colores
reales de la Carta Cronoestratigráfica Internacional**, los mismos que el
programa usa en el cronograma del Bloque 7, de modo que la portada no es un
adorno sino una figura correcta—, con sus **probabilidades posteriores** en los
nodos, un **fósil** que calibra una edad con su barra de incertidumbre, dos
**pasteles de estado ancestral**, y las puntas ilustradas con dibujos originales
de lo que los ejemplos del programa usan: la hoja compuesta y el fruto de un
copal, una gota de resina con un insecto dentro, una flor, un cono y un helecho.
Al pie, dos viñetas: un **árbol circular con k grupos** y una **red de
divisiones** donde los datos no caben en un árbol.

Nada de la portada se descarga: las tipografías vienen de Google Fonts y todo lo
demás está dibujado en el propio archivo.

## Cómo se une el manual

Desde la carpeta `manual/`:

```
powershell -ExecutionPolicy Bypass -File herramientas\unir-manual.ps1 es
```

Escribe `es/manual-completo.html` con la portada como primera hoja y las catorce
partes en orden, con un solo `paginar.js`, de modo que **la numeración de páginas
sea continua y el índice general encuentre las páginas de todos los capítulos**.
Los estilos propios de cada parte se funden; si dos partes escriben el mismo
selector con valores distintos, gana el más frecuente y el otro se acota a las
hojas de su capítulo con `.hoja[data-pestana="N"]`.

**El archivo generado no se edita a mano**: se corrigen las partes y se vuelve a
correr. Y **el script tiene que guardarse con marca de orden de bytes (BOM)**:
Windows PowerShell 5.1 lee como ANSI un `.ps1` sin BOM, y el `·` del título sale
como `Â·` en el documento unido.

Estado al 25 de septiembre de 2026: **112 hojas** (portada + 111), la última
numerada **109**, sin desbordes, con el índice completo.

## Cómo se obtiene el PDF

El PDF se arma **al final**, cuando todas las partes están escritas y revisadas. Desde
la carpeta `manual/`, con el servidor de la app corriendo:

```
powershell -ExecutionPolicy Bypass -File herramientas\imprimir-pdf.ps1
```

El guion maneja Chrome por el protocolo de DevTools. **No se usa `--print-to-pdf`**:
ese imprime en cuanto la página carga, es decir, antes de que `paginar.js` reparta las
hojas y antes de que terminen de llegar las tipografías. El guion espera a que las tres
familias estén cargadas —comprobando Cormorant en el peso **500**, que es el que pide
la hoja de estilo; con el 400 por omisión la comprobación falla siempre— y a que el
número de hojas deje de cambiar, y sólo entonces llama a `Page.printToPDF` con
`printBackground` activado y hoja carta sin márgenes.

También se puede imprimir a mano desde el navegador (destino **Guardar como PDF**,
márgenes **Ninguno**, **Gráficos de fondo** activado), que es como se hicieron los
manuales de las demás apps. El manual se imprime **de una sola vez**: unir PDF sueltos
pierde los enlaces del índice y reinicia la numeración.

El resultado se guarda como `manual/PhylogenyPro_Manual_de_usuario_ES.pdf`, igual que en las
demás apps. Al 25 de septiembre de 2026: **112 páginas** tamaño carta, 3.8 MB.
