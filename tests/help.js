/* PhylogenyPro — unit tests of the interpretation help.

   Two things have to be true of a registry like this, and neither is obvious
   from reading it: that every entry is complete and bilingual, and that the
   scales are real partitions — ordered, without gaps and without overlaps — so
   that every possible value falls into exactly one band. A scale with a hole in
   it would silently return "no verdict" for the values in the hole, which is
   the one failure a decision rule must not have. */

(function () {

  const H = Help.HELP;
  const KEYS = Object.keys(H);

  sec('1 · the catalogue is complete');

  t('hay fichas para todos los bloques', Object.keys(Help.BLOCK_KEYS).length === 12);
  t('y suficientes en total', KEYS.length >= 40);

  {
    const faltan = KEYS.filter(k => !H[k].t || !H[k].what || !H[k].read);
    t('cada ficha dice qué es y cómo se lee', faltan.length === 0);
  }
  {
    /* toda la prosa va en pares [es, en]: si alguno quedara suelto, la app
       mostraría español en la versión en inglés sin avisar */
    const malos = [];
    KEYS.forEach(k => {
      ['t', 'what', 'read', 'care', 'scaleTitle'].forEach(f => {
        const v = H[k][f];
        if (v == null) return;
        if (!Array.isArray(v) || v.length !== 2 || !v[0] || !v[1]) malos.push(k + '.' + f);
      });
    });
    t('todo el texto está en los dos idiomas', malos.length === 0);
  }
  {
    const sinCuidado = KEYS.filter(k => !H[k].care);
    t('casi todas advierten del error más común', sinCuidado.length <= KEYS.length * 0.35);
  }
  {
    const todas = new Set();
    Object.keys(Help.BLOCK_KEYS).forEach(n => Help.BLOCK_KEYS[n].forEach(k => todas.add(k)));
    const huerfanas = [...todas].filter(k => !H[k]);
    t('ningún bloque apunta a una ficha que no existe', huerfanas.length === 0);
    const sueltas = KEYS.filter(k => !todas.has(k));
    t('y ninguna ficha se queda sin bloque', sueltas.length === 0);
  }
  {
    const malas = Object.keys(Help.LABELS).filter(l => !H[Help.LABELS[l]]);
    t('las etiquetas que se decoran apuntan a fichas reales', malas.length === 0);
  }

  sec('2 · the scales are real partitions');

  {
    const rotas = [];
    KEYS.forEach(k => {
      const sc = H[k].scale;
      if (!sc) return;
      for (let i = 0; i < sc.length; i++) {
        const b = sc[i];
        if (b.from != null && b.to != null && !(b.to > b.from)) rotas.push(k + ' banda ' + i + ': vacía');
        if (i > 0 && sc[i - 1].to !== b.from) rotas.push(k + ' banda ' + i + ': no empalma con la anterior');
      }
      if (sc[0].from != null && sc[0].from !== 0) rotas.push(k + ': la primera banda no abre en 0 ni al infinito');
    });
    t('cada escala está ordenada y sin huecos ni traslapes', rotas.length === 0);
  }
  {
    const sinTono = [];
    KEYS.forEach(k => (H[k].scale || []).forEach((b, i) => {
      if (!['good', 'ok', 'warn', 'bad', 'neutral'].includes(b.tone || 'neutral')) sinTono.push(k + '/' + i);
    }));
    t('todos los tonos son de los cinco definidos', sinTono.length === 0);
  }
  {
    const sinTitulo = KEYS.filter(k => H[k].scale && !H[k].scaleTitle);
    t('cada escala dice de qué cantidad es', sinTitulo.length === 0);
  }
  {
    /* Una escala que no es una ley estadística tiene que declararse como
       convención: la mitad de los umbrales de esta disciplina son costumbres. */
    const conEscala = KEYS.filter(k => H[k].scale);
    const marcadas = conEscala.filter(k => H[k].conv);
    t('la mayoría de las escalas se declaran convención', marcadas.length >= conEscala.length * 0.5);
  }

  sec('3 · the decision rule, applied to a number');

  t('bootstrap 45 no tiene apoyo', Help.verdict('bootstrap', 45) === 'sin apoyo');
  t('bootstrap 72 está apoyada', Help.verdict('bootstrap', 72) === 'apoyada');
  t('bootstrap 96 es fuerte', Help.verdict('bootstrap', 96) === 'fuerte');

  /* el error clásico: creer que las dos escalas de apoyo son la misma */
  t('UFBoot 90, en cambio, es dudosa', Help.verdict('ufboot', 90) === 'dudosa');
  t('y UFBoot 96 sí está apoyada', Help.verdict('ufboot', 96) === 'apoyada');
  t('las dos escalas no coinciden en 90',
    Help.verdict('bootstrap', 90) !== Help.verdict('ufboot', 90));

  t('ESS 150 es justo', Help.verdict('ess', 150) === 'justo');
  t('ESS 240 es suficiente', Help.verdict('ess', 240) === 'suficiente');
  t('PSRF 1.002 convergió', Help.verdict('psrf', 1.002) === 'convergieron');
  t('PSRF 1.2 no', Help.verdict('psrf', 1.2) === 'no convergieron');
  t('γ = −2 es desaceleración', Help.verdict('gamma', -2) === 'desaceleración');
  t('γ = 0 es compatible con tasa constante', Help.verdict('gamma', 0).indexOf('constante') > 0);
  t('gCF 35 es conflicto fuerte', Help.verdict('gcf', 35) === 'conflicto fuerte');
  t('sCF 34 no tiene señal', Help.verdict('scf', 34) === 'sin señal');

  {
    /* el límite pertenece a la banda de arriba, y hay que fijarlo */
    t('el valor del límite cae en la banda superior', Help.verdict('bootstrap', 70) === 'apoyada');
    t('y justo por debajo, en la inferior', Help.verdict('bootstrap', 69.999) === 'débil');
  }
  {
    /* ningún valor razonable debe quedarse sin veredicto */
    const huecos = [];
    KEYS.filter(k => H[k].scale).forEach(k => {
      const sc = H[k].scale;
      const lo = sc[0].from == null ? -5 : sc[0].from;
      const hi = sc[sc.length - 1].to == null ? (lo + 100) : sc[sc.length - 1].to;
      for (let i = 0; i <= 40; i++) {
        const v = lo + (hi - lo) * i / 40;
        if (!Help.band(k, v)) huecos.push(k + ' @ ' + v);
      }
    });
    t('ningún valor del rango se queda sin banda', huecos.length === 0);
  }
  t('un valor que no es número no inventa veredicto',
    Help.band('bootstrap', null) === null && Help.band('bootstrap', NaN) === null);
  t('una ficha sin escala tampoco', Help.band('hpd', 1) === null);
  t('y una clave inexistente no rompe nada', Help.band('no existe', 1) === null);

  sec('4 · what it prints');

  t('la insignia lleva su clave', Help.badge('ci').indexOf('data-help-key="ci"') > 0);
  t('y no se imprime para una clave desconocida', Help.badge('inventada') === '');
  t('el chip trae el tono de la banda', Help.tag('bootstrap', 96).indexOf('help-t-good') > 0);
  t('y no aparece si no hay veredicto', Help.tag('bootstrap', null) === '');
  {
    const cuerpo = Help.entryBody(H.bootstrap);
    t('la ficha imprime su escala', cuerpo.indexOf('help-bands') > 0);
    t('y su advertencia', cuerpo.indexOf('help-care') > 0);
    t('y dice que el umbral es una convención', cuerpo.indexOf('help-conv') > 0);
  }
  {
    const g = Help.panel(['ci', 'ri']);
    t('la guía de un bloque es un desplegable', /^<details/.test(g));
    t('y no se llena hasta abrirla', g.indexOf('help-entry') < 0);
  }

})();
