/* PhylogenyPro — molecular phylogenetics in the browser.
   Copyright (C) 2026  Luis Ángel Barrera-Guzmán

   This program is free software: you can redistribute it and/or modify it under
   the terms of the GNU General Public License as published by the Free Software
   Foundation, either version 3 of the License, or (at your option) any later
   version. It is distributed in the hope that it will be useful, but WITHOUT ANY
   WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
   A PARTICULAR PURPOSE. See the GNU General Public License, in the file LICENSE
   at the root of this program, or <https://www.gnu.org/licenses/>. */

/* PhylogenyPro — the catalogue of substitution models and how one is chosen.

   Choosing a model is not a matter of taste: it is a comparison of fit against
   cost. Every candidate is fitted on the same fixed topology and the same
   branch lengths, and then ranked by AIC, corrected AIC and BIC, with Akaike
   weights so that "how much better" has a number. The hierarchical likelihood
   ratio test is offered too, because it is what many reviewers still expect,
   with the warning that it only compares nested models and depends on the order
   the tests are done in.

   The candidate lists follow jModelTest and ModelFinder: seven rate schemes for
   DNA (JC, K80, TrN, K81, TIM, TVM, SYM) with equal or empirical frequencies,
   each with nothing, +I, +G or +I+G; the usual empirical matrices for proteins;
   and Mk with or without the correction for only-variable characters (Mkv,
   Lewis 2001) for morphology. */

const Models = {};

(function () {

  const DNA_SCHEMES = [
    { ef: 'JC', uf: 'F81', label: 'JC / F81', rates: 0 },
    { ef: 'K80', uf: 'HKY', label: 'K80 / HKY', rates: 1 },
    { ef: 'TrNef', uf: 'TrN', label: 'TrN', rates: 2 },
    { ef: 'K81', uf: 'K81uf', label: 'K81', rates: 2 },
    { ef: 'TIMef', uf: 'TIM', label: 'TIM', rates: 3 },
    { ef: 'TVMef', uf: 'TVM', label: 'TVM', rates: 4 },
    { ef: 'SYM', uf: 'GTR', label: 'SYM / GTR', rates: 5 },
  ];
  const QUICK_DNA = ['JC', 'F81', 'K80', 'HKY', 'TrN', 'TIM', 'TVM', 'SYM', 'GTR'];
  const AA_QUICK = ['LG', 'WAG', 'JTT', 'cpREV', 'Dayhoff', 'VT'];
  const AA_ALL = ['LG', 'WAG', 'JTT', 'Dayhoff', 'cpREV', 'mtREV24', 'VT', 'RtREV', 'Blosum62', 'FLU', 'HIVb', 'mtmam'];

  /* the four rate-variation variants every model comes in */
  const VARIANTS = [
    { tag: '', alpha: null, pInv: 0 },
    { tag: '+I', alpha: null, pInv: 0.2 },
    { tag: '+G', alpha: 0.5, pInv: 0 },
    { tag: '+I+G', alpha: 0.5, pInv: 0.2 },
  ];

  function name(spec) {
    let n = spec.model;
    if (spec.type === 'aa' && spec.empiricalF) n += '+F';
    if (spec.pInv) n += '+I';
    if (spec.alpha) n += '+G';
    return n;
  }

  /* build the list of candidates */
  function buildList(type, opts) {
    opts = opts || {};
    const set = opts.set || 'standard';
    const out = [];
    if (type === 'dna') {
      DNA_SCHEMES.forEach(s => {
        [s.ef, s.uf].forEach(m => {
          if (set === 'quick' && QUICK_DNA.indexOf(m) < 0) return;
          VARIANTS.forEach(v => {
            if (set === 'quick' && v.tag === '+I') return;
            out.push({ type: 'dna', model: m, rates: [1, 1, 1, 1, 1, 1], alpha: v.alpha, pInv: v.pInv, ncat: opts.ncat || 4 });
          });
        });
      });
    } else if (type === 'aa') {
      const list = set === 'quick' ? AA_QUICK : AA_ALL;
      list.forEach(m => {
        [false, true].forEach(empF => {
          if (set === 'quick' && empF) return;
          VARIANTS.forEach(v => {
            if (set === 'quick' && v.tag === '+I') return;
            out.push({ type: 'aa', model: m, empiricalF: empF, alpha: v.alpha, pInv: v.pInv, ncat: opts.ncat || 4 });
          });
        });
      });
    } else {
      /* morphology: Mk with equal rates, with or without the correction for
         having scored only variable characters */
      [false, true].forEach(asc => {
        [null, 0.5].forEach(alpha => {
          out.push({ type: 'morph', model: asc ? 'Mkv' : 'Mk', ascBias: asc, alpha, pInv: 0, ncat: opts.ncat || 4 });
        });
      });
    }
    return out;
  }

  /* ================================================================
     the correction of Lewis (2001) for morphological matrices
     ================================================================ */
  /* When only variable characters are scored — which is what a morphologist
     does — the likelihood must be conditioned on the character being variable:
     lnL corrected = lnL − N · log(1 − P(constant)). P(constant) is computed by
     running the engine on the S artificial constant patterns. */
  function ascertainmentCorrection(tree, A, M) {
    const S = A.nStates;
    const fake = {
      nStates: S, nPat: S, nSites: S, nSeq: A.nSeq,
      weights: Float64Array.from({ length: S }, () => 1),
      tips: A.tips.map(() => {
        const arr = new Float64Array(S * S);
        for (let p = 0; p < S; p++) arr[p * S + p] = 1;
        return arr;
      }),
      constMask: Int32Array.from({ length: S }, (_, p) => 1 << p),
      freqs: A.freqs, nConstant: S, nVariable: 0,
    };
    const F = Tree.flatten(tree);
    const lik = Like.engine(F, fake, M);
    lik.full(F.len);
    const site = lik.siteLnL();
    let pConst = 0;
    for (let p = 0; p < S; p++) pConst += Math.exp(site[p]);
    return { pConst, correction: -Math.log(Math.max(1e-12, 1 - pConst)) };
  }

  /* ================================================================
     fitting every candidate
     ================================================================ */
  /* Runs one model per turn of the event loop, so the interface never freezes
     and the run can be cancelled. onProgress(done, total, lastFit) is called
     after each one. */
  function selectModel(tree, A, list, opts) {
    opts = opts || {};
    const fits = [];
    let cancelled = false;
    const state2 = { cancel() { cancelled = true; } };
    const promise = new Promise((resolve, reject) => {
      let i = 0;
      const step = () => {
        if (cancelled) { resolve({ fits, cancelled: true }); return; }
        if (i >= list.length) { resolve(finish(fits, opts)); return; }
        const spec = Object.assign({}, list[i]);
        try {
          const fitter = opts.fast === false ? Like.fit : Like.fitFast;
          const f = fitter(tree, A, spec, opts.fitOpts || {});
          if (spec.ascBias) {
            const corr = ascertainmentCorrection(tree, A, f.model);
            f.lnL -= A.nSites * Math.log(Math.max(1e-12, 1 - corr.pConst));
            f.pConst = corr.pConst;
            f.AIC = -2 * f.lnL + 2 * f.k;
            f.AICc = -2 * f.lnL + 2 * f.k + (A.nSites - f.k - 1 > 0 ? 2 * f.k * (f.k + 1) / (A.nSites - f.k - 1) : Infinity);
            f.BIC = -2 * f.lnL + f.k * Math.log(A.nSites);
          }
          f.name = name(f.spec);
          fits.push(f);
        } catch (e) {
          fits.push({ name: name(spec), spec, error: e.message, lnL: -Infinity, AIC: Infinity, AICc: Infinity, BIC: Infinity, k: 0 });
        }
        i++;
        if (opts.onProgress) opts.onProgress(i, list.length, fits[fits.length - 1]);
        setTimeout(step, 0);
      };
      setTimeout(step, 0);
    });
    return { promise, cancel: state2.cancel };
  }

  function finish(fits, opts) {
    const ok = fits.filter(f => isFinite(f.lnL));
    const best = key => ok.reduce((a, b) => (b[key] < a[key] ? b : a), ok[0]);
    const bAIC = best('AIC'), bAICc = best('AICc'), bBIC = best('BIC');
    ['AIC', 'AICc', 'BIC'].forEach(key => {
      const b = key === 'AIC' ? bAIC : key === 'AICc' ? bAICc : bBIC;
      let sum = 0;
      ok.forEach(f => { f['d' + key] = f[key] - b[key]; sum += Math.exp(-f['d' + key] / 2); });
      ok.forEach(f => { f['w' + key] = Math.exp(-f['d' + key] / 2) / sum; });
    });
    ok.sort((a, b) => a.AICc - b.AICc);
    return {
      fits: ok, all: fits,
      bestAIC: bAIC, bestAICc: bAICc, bestBIC: bBIC,
      lrt: hLRT(fits),
    };
  }

  /* ================================================================
     hierarchical likelihood ratio tests
     ================================================================ */
  /* Each step compares two nested models; the chain is the classical one of
     Modeltest: equal frequencies → unequal, one rate → transitions apart, and
     so on, then +I and +G on top of the winner. */
  const CHAIN = [
    ['JC', 'F81', 'frecuencias de base desiguales', 'unequal base frequencies', 3],
    ['F81', 'HKY', 'transiciones distintas de transversiones', 'transitions differ from transversions', 1],
    ['HKY', 'TrN', 'las dos transiciones con tasas distintas', 'the two transitions with different rates', 1],
    ['TrN', 'TIM', 'transversiones con tasas distintas', 'transversions with different rates', 1],
    ['TIM', 'GTR', 'las seis tasas libres', 'all six rates free', 2],
  ];
  function hLRT(fits) {
    const byName = new Map(fits.map(f => [f.name, f]));
    const out = [];
    let current = 'JC';
    CHAIN.forEach(([a, b, es, en, df]) => {
      const fa = byName.get(a), fb = byName.get(b);
      if (!fa || !fb || !isFinite(fa.lnL) || !isFinite(fb.lnL)) return;
      const lr = 2 * (fb.lnL - fa.lnL);
      const p = 1 - Dist.pchisq(Math.max(0, lr), df);
      const accepted = p < 0.01;
      out.push({ from: a, to: b, es, en, lr, df, p, accepted });
      if (accepted) current = b;
    });
    /* +G and +I on the winner of the chain */
    [['+G', 1, 'variación de tasas entre sitios', 'rate variation among sites'],
     ['+I', 1, 'una fracción de sitios invariables', 'a fraction of invariable sites']].forEach(([suffix, df, es, en]) => {
      const fa = byName.get(current), fb = byName.get(current + suffix);
      if (!fa || !fb) return;
      const lr = 2 * (fb.lnL - fa.lnL);
      const p = 1 - Dist.pchisq(Math.max(0, lr), df);
      const accepted = p < 0.01;
      out.push({ from: current, to: current + suffix, es, en, lr, df, p, accepted });
      if (accepted) current = current + suffix;
    });
    return { steps: out, model: current };
  }

  /* ================================================================
     partition schemes
     ================================================================ */
  /* Greedy merging in the style of PartitionFinder: start with one model per
     partition, then try merging the pair whose merged model loses the least,
     and keep merging while the score improves. The score is BIC by default,
     because AIC merges too little on long alignments. */
  function mergeScheme(parts, tree, opts) {
    opts = opts || {};
    const key = opts.criterion || 'BIC';
    const list = parts.map((p, i) => ({ ids: [i], name: p.name, A: p.A, fit: p.fit }));
    const scoreOf = group => group.fit[key];
    let total = list.reduce((s, g) => s + scoreOf(g), 0);
    const steps = [{ scheme: list.map(g => g.name), score: total, k: list.length }];
    const fitGroup = ids => {
      /* concatenating the compressed data is not valid, so the sequences are
         concatenated and compressed again */
      const seqs = parts[0].seqs.map((_, r) => ids.map(i => parts[i].seqs[r]).join(''));
      const A2 = Like.compress(seqs, parts[ids[0]].type);
      const spec = Object.assign({}, opts.spec || { type: parts[ids[0]].type, model: 'GTR', alpha: 0.5, ncat: 4 });
      return { A: A2, fit: Like.fitFast(tree, A2, spec, { passes: 4 }) };
    };
    let changed = true;
    while (changed && list.length > 1) {
      changed = false;
      let best = null;
      for (let a = 0; a < list.length; a++) for (let b = a + 1; b < list.length; b++) {
        const ids = list[a].ids.concat(list[b].ids);
        const merged = fitGroup(ids);
        const delta = merged.fit[key] - (scoreOf(list[a]) + scoreOf(list[b]));
        if (!best || delta < best.delta) best = { a, b, ids, merged, delta };
      }
      if (best && best.delta < 0) {
        const a = list[best.a], b = list[best.b];
        const merged = { ids: best.ids, name: a.name + '+' + b.name, A: best.merged.A, fit: best.merged.fit };
        list.splice(best.b, 1); list.splice(best.a, 1);
        list.push(merged);
        total = list.reduce((s, g) => s + scoreOf(g), 0);
        steps.push({ scheme: list.map(g => g.name), score: total, k: list.length, merged: merged.name, delta: best.delta });
        changed = true;
      }
    }
    return { groups: list, score: total, steps, criterion: key };
  }

  Object.assign(Models, {
    DNA_SCHEMES, VARIANTS, AA_QUICK, AA_ALL, QUICK_DNA,
    name, buildList, selectModel, hLRT, mergeScheme, ascertainmentCorrection, finish,
  });
  window.Models = Models;
})();
