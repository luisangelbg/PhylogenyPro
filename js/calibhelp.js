/* PhylogenyPro — molecular phylogenetics in the browser.
   Copyright (C) 2026  Luis Ángel Barrera-Guzmán

   This program is free software: you can redistribute it and/or modify it under
   the terms of the GNU General Public License as published by the Free Software
   Foundation, either version 3 of the License, or (at your option) any later
   version. It is distributed in the hope that it will be useful, but WITHOUT ANY
   WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
   A PARTICULAR PURPOSE. See the GNU General Public License, in the file LICENSE
   at the root of this program, or <https://www.gnu.org/licenses/>.

   ---------------------------------------------------------------------------

   The calibration assistant.

   The program does not ship a database of fossils and never will: assigning a
   fossil to a node is a taxonomic judgement about that fossil's characters and
   about which clade it is sister to, and it is the part of a dating study that
   reviewers argue about. A lookup table would invite people to publish a date
   nobody examined.

   What it can do is everything around that judgement, which is where the
   avoidable mistakes live:

     · suggest the shape of prior that matches the SOURCE of the age — a fossil
       is a minimum, a geological event is a window, an age borrowed from
       another study is a normal with somebody else's uncertainty inside it;
     · say out loud what the prior the user just typed actually asserts, in
       years, because "lognormal, M = 1.5, S = 0.6, offset 20" means nothing to
       most people and "median 24.5, 95 % between 21.4 and 34.5" means
       everything;
     · and warn about the handful of errors that turn up again and again: a
       fossil entered as a fixed age, a crown calibrated with a stem fossil, a
       prior so wide it says nothing, a calibration younger than its own
       descendant, and a whole chronogram hanging from a single date.

   Nothing here decides anything. It tells the user what they just said. */

const CalibHelp = {};

(function () {

  const ln = Math.log, exp = Math.exp;
  const Z = 1.959963985;                       /* the 95 % of a normal */

  /* =====================================================================
     1 · where the age comes from
     ===================================================================== */
  const SOURCES = [
    {
      id: 'fossil',
      name: ['Un fósil asignado a este clado', 'A fossil assigned to this clade'],
      prior: 'lognormal',
      why: ['Un fósil da la <b>edad mínima</b> del nodo, nunca su edad: el linaje existía ya cuando ese individuo murió, y pudo existir mucho antes. Por eso se usa un prior asimétrico y desplazado, que permite que el nodo sea mucho más viejo pero no más joven que el fósil.',
        'A fossil gives the node\'s <b>minimum age</b>, never its age: the lineage already existed when that individual died, and may have existed long before. Hence an asymmetric, offset prior, which lets the node be much older but not younger than the fossil.'],
    },
    {
      id: 'geo',
      name: ['Un acontecimiento geológico', 'A geological event'],
      prior: 'uniform',
      why: ['El levantamiento de una isla, el cierre de un istmo o la apertura de un mar dan una <b>ventana</b> de tiempo con principio y fin, no un punto. El intervalo uniforme dice justo eso y nada más.',
        'An island emerging, an isthmus closing or a seaway opening give a <b>window</b> with a beginning and an end, not a point. The uniform interval says exactly that and nothing more.'],
    },
    {
      id: 'secondary',
      name: ['Una edad tomada de otro estudio', 'An age taken from another study'],
      prior: 'normal',
      why: ['Una calibración secundaria es una <b>opinión heredada</b>: arrastra el árbol, el modelo y las calibraciones de aquel trabajo. Se representa con una normal centrada en su media y con su desviación, y conviene <b>ensancharla</b>, no estrecharla: la incertidumbre real es siempre mayor que la publicada.',
        'A secondary calibration is an <b>inherited opinion</b>: it drags in that study\'s tree, model and calibrations. It is represented by a normal at its mean with its standard deviation, and it is worth <b>widening</b> it, not narrowing it: the real uncertainty is always larger than the published one.'],
    },
    {
      id: 'known',
      name: ['Una edad que se conoce con certeza', 'An age known for certain'],
      prior: 'fixed',
      why: ['Casi nunca existe en biología. Se usa para fijar la raíz en una edad de referencia cuando el objetivo es comparar nodos entre sí, no publicar fechas. Si el dato viene de un fósil, <b>ésta no es la opción</b>.',
        'It almost never exists in biology. It is used to pin the root at a reference age when the aim is to compare nodes with each other, not to publish dates. If the figure comes from a fossil, <b>this is not the option</b>.'],
    },
  ];
  const sourceOf = id => SOURCES.find(s => s.id === id) || SOURCES[0];

  /* =====================================================================
     2 · what the prior the user typed actually says
     ===================================================================== */
  /* Every shape is reduced to the same three numbers — the 2.5 %, the median
     and the 97.5 % — because that is what a reader of the paper will want and
     what the author almost never checks. */
  function quantiles(c) {
    if (!c) return null;
    const t = c.type;
    if (t === 'fixed') {
      const v = +c.value;
      return isFinite(v) ? { lo: v, med: v, hi: v, exact: true } : null;
    }
    if (t === 'uniform') {
      const a = +c.min, b = +c.max;
      if (!isFinite(a) || !isFinite(b) || b <= a) return null;
      return { lo: a + 0.025 * (b - a), med: (a + b) / 2, hi: a + 0.975 * (b - a), floor: a, ceiling: b };
    }
    if (t === 'normal') {
      const m = +c.mean, s = +c.sd;
      if (!isFinite(m) || !isFinite(s) || s <= 0) return null;
      return { lo: m - Z * s, med: m, hi: m + Z * s };
    }
    if (t === 'lognormal') {
      /* edad = desplazamiento + lognormal(M, S), con M y S en escala logarítmica */
      const off = +c.min, M = +c.M, S = +c.S;
      if (!isFinite(off) || !isFinite(M) || !isFinite(S) || S <= 0) return null;
      return { lo: off + exp(M - Z * S), med: off + exp(M), hi: off + exp(M + Z * S), floor: off };
    }
    if (t === 'exponential') {
      /* edad = desplazamiento + exponencial de media dada */
      const off = +c.min, mu = +c.mean;
      if (!isFinite(off) || !isFinite(mu) || mu <= 0) return null;
      return { lo: off - mu * ln(0.975), med: off + mu * ln(2), hi: off - mu * ln(0.025), floor: off };
    }
    return null;
  }

  const f2 = v => (Math.abs(v) >= 100 ? v.toFixed(0) : Math.abs(v) >= 10 ? v.toFixed(1) : v.toFixed(2));

  /* one sentence saying what the prior asserts */
  function describe(c, unit) {
    const q = quantiles(c);
    const u = unit ? ' ' + unit : '';
    if (!q) return ['Faltan datos o no son válidos para esta forma de prior.',
      'Data are missing or invalid for this prior shape.'];
    if (q.exact) return [`El nodo queda fijo en ${f2(q.med)}${u}, sin incertidumbre.`,
      `The node is pinned at ${f2(q.med)}${u}, with no uncertainty.`];
    const piso = q.floor != null
      ? [`, y no puede ser más joven que ${f2(q.floor)}${u}`, `, and cannot be younger than ${f2(q.floor)}${u}`]
      : ['', ''];
    return [
      `Este prior dice que la edad del nodo está, con 95 % de probabilidad, entre <b>${f2(q.lo)}</b> y <b>${f2(q.hi)}</b>${u}, con mediana <b>${f2(q.med)}</b>${u}${piso[0]}.`,
      `This prior says the node's age lies, with 95 % probability, between <b>${f2(q.lo)}</b> and <b>${f2(q.hi)}</b>${u}, with a median of <b>${f2(q.med)}</b>${u}${piso[1]}.`,
    ];
  }

  /* =====================================================================
     3 · the mistakes that turn up again and again
     ===================================================================== */
  /* ctx: { source, others: [calibs], isRoot, nTips, cladeSize, unit } */
  function advise(c, ctx) {
    ctx = ctx || {};
    const out = [];
    const add = (level, es, en) => out.push({ level, es, en });
    const q = quantiles(c);
    const u = ctx.unit ? ' ' + ctx.unit : '';

    /* — la confusión que más caro sale — */
    if (ctx.source === 'fossil' && c.type === 'fixed') {
      add('bad',
        'Un fósil <b>no da la edad de un nodo</b>: da su edad mínima. Fijar el nodo en la edad del fósil afirma que el linaje nació justo cuando murió ese individuo, lo que es casi con seguridad falso y comprime todo el árbol hacia el presente. Usa una lognormal o una exponencial desplazadas, con el desplazamiento en la edad del fósil.',
        'A fossil <b>does not give a node\'s age</b>: it gives its minimum age. Pinning the node at the fossil\'s age asserts that the lineage arose exactly when that individual died, which is almost certainly false and squeezes the whole tree towards the present. Use an offset lognormal or exponential, with the offset at the fossil\'s age.');
    }
    if (ctx.source === 'fossil' && c.type === 'normal') {
      add('warn',
        'Una normal es simétrica: permite que el nodo sea <b>más joven</b> que el fósil, lo que es imposible. Para un fósil, una lognormal o una exponencial desplazadas.',
        'A normal is symmetric: it allows the node to be <b>younger</b> than the fossil, which is impossible. For a fossil, use an offset lognormal or exponential.');
    }
    if ((c.type === 'lognormal' || c.type === 'exponential') && +c.min === 0) {
      add('warn',
        'El desplazamiento está en cero, así que este prior <b>no impone ninguna edad mínima</b>: permite que el nodo sea de ayer. Si viene de un fósil, pon su edad en el desplazamiento.',
        'The offset is zero, so this prior <b>imposes no minimum age</b>: it allows the node to be from yesterday. If it comes from a fossil, put its age in the offset.');
    }
    if (c.type === 'uniform' && isFinite(+c.min) && isFinite(+c.max) && +c.max <= +c.min) {
      add('bad',
        'El máximo no es mayor que el mínimo: el intervalo está vacío.',
        'The maximum is not greater than the minimum: the interval is empty.');
    }
    if (c.type === 'normal' && q && q.lo < 0) {
      add('warn',
        `Con esa media y esa desviación, el 2.5 % inferior del prior cae en <b>${f2(q.lo)}${u}</b>, es decir en el futuro. Reduce la desviación o usa una lognormal.`,
        `With that mean and standard deviation, the lower 2.5 % of the prior falls at <b>${f2(q.lo)}${u}</b>, that is, in the future. Reduce the standard deviation or use a lognormal.`);
    }

    /* — priors que no dicen nada —
       La anchura se mide como la RAZÓN entre los dos extremos, no como su
       diferencia relativa a la mediana: con un intervalo uniforme esa
       diferencia no puede pasar de 1.9 por construcción, de modo que el aviso
       nunca se habría disparado justo en la forma que se usa para ir a lo
       vago. La razón no depende de la escala ni de la forma. */
    if (q && !q.exact && q.lo > 0 && q.hi / q.lo > 4) {
      add('warn',
        'El intervalo del 95 % es más de tres veces la mediana: este prior es tan ancho que apenas informa. No está mal —a veces es lo honesto—, pero no esperes que fije la escala del árbol.',
        'The 95 % interval is more than three times the median: this prior is so wide that it barely informs. That is not wrong — sometimes it is the honest thing — but do not expect it to set the tree\'s scale.');
    }
    if (q && !q.exact && q.med > 0 && (q.hi - q.lo) / q.med < 0.05) {
      add('warn',
        'El intervalo del 95 % es menos del 5 % de la mediana: es un prior <b>muy seguro</b>. Para una calibración secundaria eso suele ser un exceso de confianza, porque la incertidumbre real incluye la del estudio de donde salió.',
        'The 95 % interval is less than 5 % of the median: this is a <b>very confident</b> prior. For a secondary calibration that is usually overconfidence, because the real uncertainty includes that of the study it came from.');
    }

    /* — corona contra tallo — */
    if (ctx.cladeSize != null && ctx.cladeSize >= 2 && ctx.source === 'fossil') {
      add('info',
        'Los taxones que elegiste calibran el <b>nodo de corona</b> de ese grupo: su ancestro común más reciente. Si tu fósil es <b>hermano</b> del grupo y no un miembro suyo, lo que calibra es el <b>tallo</b>, que es más viejo: para eso hay que seleccionar también al grupo hermano. Confundir los dos desplaza la fecha, y es de los errores más frecuentes al datar.',
        'The taxa you chose calibrate the <b>crown node</b> of that group: its most recent common ancestor. If your fossil is <b>sister</b> to the group rather than a member of it, what it calibrates is the <b>stem</b>, which is older: for that you have to select the sister group as well. Confusing the two shifts the date, and it is one of the commonest mistakes in dating.');
    }

    /* — el conjunto de calibraciones — */
    const otras = (ctx.others || []).filter(Boolean);
    if (otras.length === 0) {
      add('info',
        'Ésta será la <b>única</b> calibración: toda la escala temporal del árbol dependerá de ella, y su error se propagará a cada nodo sin que nada lo contrapese. Dos o tres calibraciones repartidas por el árbol son mucho mejores que una sola, por buena que sea.',
        'This will be the <b>only</b> calibration: the whole time scale of the tree will depend on it, and its error will propagate to every node with nothing to balance it. Two or three calibrations spread over the tree are far better than one, however good.');
    }
    if (ctx.isRoot && otras.length === 0) {
      add('info',
        'Calibrar sólo la raíz fija la escala pero no comprueba nada: cualquier error queda repartido por todo el árbol sin que se note. Con una calibración interna adicional, las dos se contrastan entre sí.',
        'Calibrating only the root sets the scale but checks nothing: any error is spread over the whole tree unnoticed. With one additional internal calibration, the two check each other.');
    }

    return out;
  }

  /* =====================================================================
     4 · suggesting parameters from a source and an age
     ===================================================================== */
  /* Given the source and the age the user has in hand, propose a prior that a
     reviewer would accept, and say why. The numbers are a starting point, not
     a recommendation: the offset is the fossil's age, and the spread is the
     conventional one — it is the user who knows how confident to be. */
  function suggest(sourceId, age, age2) {
    const s = sourceOf(sourceId);
    const a = +age;
    if (!isFinite(a) || a < 0) return null;
    if (s.id === 'fossil') {
      /* lognormal desplazada: mediana un 25 % por encima del fósil, cola larga */
      const M = ln(Math.max(0.5, a * 0.25));
      return { type: 'lognormal', min: a, M: +M.toFixed(3), S: 0.7,
        nota: ['Desplazada en la edad del fósil, con la mediana algo por encima y una cola larga hacia atrás.',
          'Offset at the fossil\'s age, with the median a little above it and a long tail backwards.'] };
    }
    if (s.id === 'geo') {
      const b = isFinite(+age2) && +age2 > a ? +age2 : a * 1.2;
      return { type: 'uniform', min: a, max: +b.toFixed(3),
        nota: ['El intervalo del acontecimiento, sin preferir ningún punto dentro de él.',
          'The window of the event, with no point inside it preferred.'] };
    }
    if (s.id === 'secondary') {
      return { type: 'normal', mean: a, sd: +(a * 0.12).toFixed(3),
        nota: ['Centrada en la edad publicada, con una desviación del 12 %. Ensánchala si el estudio de origen tenía pocas calibraciones.',
          'Centred on the published age, with a 12 % standard deviation. Widen it if the source study had few calibrations.'] };
    }
    return { type: 'fixed', value: a,
      nota: ['Sin incertidumbre: sólo para fijar una escala de referencia.',
        'With no uncertainty: only to pin a reference scale.'] };
  }

  Object.assign(CalibHelp, { SOURCES, sourceOf, quantiles, describe, advise, suggest, _f2: f2, Z });
  if (typeof window !== 'undefined') window.CalibHelp = CalibHelp;

})();
