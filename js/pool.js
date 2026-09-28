/* PhylogenyPro — molecular phylogenetics in the browser.
   Copyright (C) 2026  Luis Ángel Barrera-Guzmán

   This program is free software: you can redistribute it and/or modify it under
   the terms of the GNU General Public License as published by the Free Software
   Foundation, either version 3 of the License, or (at your option) any later
   version. It is distributed in the hope that it will be useful, but WITHOUT ANY
   WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
   A PARTICULAR PURPOSE. See the GNU General Public License, in the file LICENSE
   at the root of this program, or <https://www.gnu.org/licenses/>. */

/* PhylogenyPro — running the slow things off the main thread.

   A bootstrap of maximum likelihood is a hundred independent searches, and a
   Bayesian analysis is several independent chains: both are embarrassingly
   parallel, and both take minutes. Doing them on the page freezes the browser;
   doing them in Web Workers uses every core the machine has and leaves the
   interface alive.

   The awkward part is that a Worker cannot `importScripts` a sibling file when
   the app is opened by double-click (a file:// page cannot fetch its own
   neighbours). The way round it is the one the other LABG apps use: every
   module of the engine is a *named function* that installs itself onto whatever
   global it is handed, so the worker's source can be assembled from those
   functions' own `toString()` and turned into a Blob. Nothing is fetched, and
   the same code runs on the page and in the worker — there is no second
   implementation to keep in step.

   If Workers are unavailable or silent (some privacy settings block Blob
   workers), the same jobs run on the page in small chunks with setTimeout
   between them, slower but never frozen. Everything is seeded, so a run can be
   reproduced whichever path it took. */

const Pool = {};

(function () {

  /* ================================================================
     1 · the code a worker runs
     ================================================================ */
  /* The seeded generator, repeated here because core.js is full of DOM helpers
     that have no business inside a worker. Same algorithm (sfc32), so a seed
     gives the same stream on the page and off it. */
  function WorkerRuntime(g) {
    g.rng = function rng(seed) {
      let a = 0x9e3779b9, b = 0x243f6a88, c = 0xb7e15162, d = (seed >>> 0) ^ 0xdeadbeef;
      const next = () => {
        a >>>= 0; b >>>= 0; c >>>= 0; d >>>= 0;
        let t = (a + b) | 0;
        a = b ^ (b >>> 9);
        b = (c + (c << 3)) | 0;
        c = (c << 21) | (c >>> 11);
        d = (d + 1) | 0;
        t = (t + d) | 0;
        c = (c + t) | 0;
        return (t >>> 0) / 4294967296;
      };
      for (let i = 0; i < 15; i++) next();
      return next;
    };

    /* ---- the jobs a worker knows how to do ---------------------- */
    g.Jobs = {
      /* one bootstrap replicate: resample the columns with replacement, then
         search. Returns the tree as Newick, which is all the caller needs. */
      mlBootstrap(P, job) {
        const r = g.rng(job.seed);
        const L = P.seqs[0].length;
        const pick = new Int32Array(L);
        for (let i = 0; i < L; i++) pick[i] = Math.floor(r() * L);
        const seqs = P.seqs.map(s => {
          const b = new Array(L);
          for (let i = 0; i < L; i++) b[i] = s[pick[i]];
          return b.join('');
        });
        const A = g.Like.compress(seqs, P.type || 'dna');
        const spec = Object.assign({}, P.spec);
        if (spec.freqs === 'empirical' || !spec.freqs) spec.freqs = Array.from(A.freqs);
        const D = g.Dist.matrix(seqs, P.distance || 'jc', { type: P.type || 'dna' }).D;
        const start = g.Tree.bionj(D, P.labels);
        /* the model is held at the values fitted to the real data, as RAxML and
           IQ-TREE do: a replicate re-estimates the tree, not the model */
        const res = g.ML.search(A, spec, Object.assign({
          start, labels: P.labels, collect: false, spr: false, maxRounds: 12,
          /* a replicate is there for its TOPOLOGY: its branch lengths are never
             published, so they are optimised roughly and no further */
          initPasses: 3, refitPasses: 2, finalPasses: 3,
          roughTol: 0.01, finalTol: 0.01,
          fixModel: P.fixModel !== false,
        }, P.search || {}));
        /* the splits are what a bootstrap is counted on; the Newick is kept so the
           replicate trees can be exported or drawn */
        return {
          newick: g.Tree.writeNewick(res.tree, { support: false }),
          splits: [...g.Tree.splits(res.tree, P.labels.length).keys()],
          lnL: res.lnL,
        };
      },

      /* one independent MCMC run, with as many heated chains as asked for.
         Two or four of these in parallel is what makes the convergence
         diagnostics possible: ASDSF and PSRF need runs that never talked to
         each other. */
      mcmcRun(P, job) {
        const A = g.Like.compress(P.seqs, P.type || 'dna');
        const spec = g.Like.cloneSpec(P.spec);
        if (!spec.freqs) spec.freqs = Array.from(A.freqs);
        const res = g.Mcmc.run(A, Object.assign({}, P.opts, {
          tree: g.Tree.parseNewick(P.startNewick, P.labels),
          labels: P.labels, spec, seed: job.seed,
        }));
        /* Maps do not survive postMessage: send the split frequencies as pairs */
        return {
          samples: res.samples, trees: res.trees,
          splitFreq: [...res.splitFreq.entries()],
          acceptance: res.acceptance, swapRate: res.swapRate,
          generations: res.generations, burnin: res.burnin,
          sampleEvery: res.sampleEvery, chains: res.chains, ms: res.ms,
        };
      },

      /* one stepping-stone ladder: the marginal likelihood of one model */
      marginalLikelihood(P, job) {
        const A = g.Like.compress(P.seqs, P.type || 'dna');
        const spec = g.Like.cloneSpec(P.spec);
        if (!spec.freqs) spec.freqs = Array.from(A.freqs);
        const res = g.Mcmc.steppingStone(A, Object.assign({}, P.opts, {
          tree: g.Tree.parseNewick(P.startNewick, P.labels),
          labels: P.labels, spec, seed: job.seed,
        }));
        return { logMarginalLikelihood: res.logMarginalLikelihood, steps: res.steps, name: P.name };
      },

      /* one parsimony bootstrap or jackknife replicate, for Block 4's support
         when the user asks for it from here */
      parsBootstrap(P, job) {
        const r = g.rng(job.seed);
        const L = P.seqs[0].length;
        const keep = [];
        if (P.jackknife) { for (let i = 0; i < L; i++) if (r() > (P.deleteFraction || 0.37)) keep.push(i); }
        else for (let i = 0; i < L; i++) keep.push(Math.floor(r() * L));
        const seqs = P.seqs.map(s => keep.map(i => s[i]).join(''));
        const A = g.Like.compress(seqs, P.type || 'dna');
        const res = g.Pars.search(A, { labels: P.labels, swap: P.swap || 'NNI', starts: P.starts || 1, seed: job.seed });
        return {
          newick: g.Tree.writeNewick(res.trees[0], { support: false }),
          splits: [...g.Tree.splits(res.trees[0], P.labels.length).keys()],
          steps: res.steps,
        };
      },
    };

    g.onmessage = ev => {
      const m = ev.data;
      if (m.type !== 'run') return;
      if (m.AAMODELS) { g.AAMODELS = m.AAMODELS; g.AA_ORDER = m.AA_ORDER; }
      const fn = g.Jobs[m.task];
      if (!fn) { g.postMessage({ type: 'error', message: 'unknown task ' + m.task }); return; }
      for (let i = 0; i < m.jobs.length; i++) {
        try {
          const out = fn(m.payload, m.jobs[i]);
          g.postMessage({ type: 'result', index: m.jobs[i].index, result: out });
        } catch (e) {
          g.postMessage({ type: 'error', index: m.jobs[i].index, message: e && e.message });
        }
      }
      g.postMessage({ type: 'done' });
    };
  }

  /* the modules a worker needs, in dependency order */
  function sources() {
    const mods = [];
    [['TreeCore', typeof TreeCore], ['DistCore', typeof DistCore], ['LikeCore', typeof LikeCore],
     ['ParsCore', typeof ParsCore], ['MLCore', typeof MLCore], ['McmcCore', typeof McmcCore]].forEach(([name, t]) => {
      if (t === 'function') mods.push(name);
    });
    if (mods.length < 6) return null;                    // something is not loaded: no worker
    const src = mods.map(n => window[n] ? window[n].toString() : eval(n).toString()).join('\n\n');
    return src + '\n\n' + WorkerRuntime.toString() + '\n' +
      mods.map(n => n + '(self);').join('\n') + '\nWorkerRuntime(self);\n';
  }

  let cachedURL = null;
  function workerURL() {
    if (cachedURL) return cachedURL;
    const src = sources();
    if (!src) return null;
    try {
      cachedURL = URL.createObjectURL(new Blob([src], { type: 'text/javascript' }));
      return cachedURL;
    } catch (e) { return null; }
  }
  function available() {
    return typeof Worker === 'function' && typeof Blob === 'function' && !!sources();
  }
  /* How many workers to start: every core but one, so the interface keeps a
     thread to itself, and never more than sixteen — beyond that the replicates
     start competing for memory bandwidth instead of for cores. */
  function cores(want) {
    const hw = (typeof navigator !== 'undefined' && navigator.hardwareConcurrency) || 4;
    return Math.max(1, Math.min(want || (hw - 1), 16, hw - 1 || 1));
  }

  /* ================================================================
     2 · running a batch of jobs
     ================================================================ */
  /* opts: { task, payload, n, workers, onProgress(done, n), onResult(i, r),
             cancelled() }
     Resolves with an array of results in job order. Falls back to the page. */
  function run(opts) {
    const n = opts.n;
    const results = new Array(n);
    let done = 0, failed = 0, cancelled = false;
    const report = () => { if (opts.onProgress) opts.onProgress(done, n, failed); };

    const url = opts.inline ? null : workerURL();
    if (!url) return runInline();

    return new Promise(resolve => {
      /* Never more workers than jobs — an idle worker that reports "done"
         immediately would otherwise end the batch before the busy ones answer —
         and never so many that each gets a single job. A fresh worker starts in
         the interpreter and only reaches optimised code after a few seconds of
         the same loops, so one job per worker pays the whole warm-up and throws
         it away: measured here, sixteen replicates on sixteen workers took
         longer than the same sixteen on six. */
      const nw = Math.max(1, Math.min(cores(opts.workers), Math.ceil(n / (opts.minPerWorker || 3))));
      const workers = [];
      let heard = false, pending = nw, gaveUp = false;
      /* jobs are dealt round-robin so every worker gets a similar share */
      const lists = Array.from({ length: nw }, () => []);
      for (let i = 0; i < n; i++) lists[i % nw].push({ index: i, seed: (opts.seed || 1) * 7919 + i * 104729 });

      const giveUp = why => {
        if (gaveUp) return;
        gaveUp = true;
        console.warn('PhylogenyPro: ' + why + ' — running on the page instead');
        workers.forEach(w => w.terminate());
        runInline().then(resolve);
      };
      /* a worker that is silently blocked must not leave the page waiting */
      const silence = setTimeout(() => { if (!heard) giveUp('the workers never answered'); }, 10000);

      const AAM = (typeof AAMODELS !== 'undefined' && opts.payload && opts.payload.type === 'aa') ? AAMODELS : null;
      for (let k = 0; k < nw; k++) {
        if (!lists[k].length) { pending--; continue; }
        let w;
        try { w = new Worker(url); } catch (e) { clearTimeout(silence); giveUp('workers are not allowed'); return; }
        workers.push(w);
        w.onmessage = e => {
          heard = true;
          const m = e.data;
          if (m.type === 'result') {
            results[m.index] = m.result;
            done++;
            if (opts.onResult) opts.onResult(m.index, m.result);
            report();
            if (opts.cancelled && opts.cancelled() && !cancelled) {
              cancelled = true;
              workers.forEach(x => x.terminate());
              clearTimeout(silence);
              resolve(results.filter(x => x !== undefined));
            }
          } else if (m.type === 'error') { failed++; done++; report(); }
          else if (m.type === 'done') {
            if (--pending <= 0) {
              clearTimeout(silence);
              workers.forEach(x => x.terminate());
              resolve(results.filter(x => x !== undefined));
            }
          }
        };
        w.onerror = err => { if (!heard) { clearTimeout(silence); giveUp('a worker failed: ' + (err && err.message)); } };
        w.postMessage({
          type: 'run', task: opts.task, payload: opts.payload, jobs: lists[k],
          AAMODELS: AAM, AA_ORDER: typeof AA_ORDER !== 'undefined' ? AA_ORDER : null,
        });
      }
      if (!workers.length) { clearTimeout(silence); resolve([]); }
    });

    /* the same jobs, on the page, one at a time, with the loop handed back
       between them so the interface keeps breathing */
    function runInline() {
      const g = window;
      if (!g.Jobs) WorkerRuntime(g);          // installs Jobs and rng onto the page
      return new Promise(resolve => {
        let i = 0;
        const step = () => {
          if (opts.cancelled && opts.cancelled()) { resolve(results.filter(x => x !== undefined)); return; }
          if (i >= n) { resolve(results.filter(x => x !== undefined)); return; }
          const job = { index: i, seed: (opts.seed || 1) * 7919 + i * 104729 };
          try {
            results[i] = g.Jobs[opts.task](opts.payload, job);
            if (opts.onResult) opts.onResult(i, results[i]);
          } catch (e) { failed++; console.warn('job ' + i + ': ' + e.message); }
          i++; done++;
          report();
          setTimeout(step, 0);
        };
        setTimeout(step, 0);
      });
    }
  }

  /* an estimate of how long a batch will take, from the first few jobs: the
     app promises a time before it starts, and this is how it keeps that promise */
  function estimate(msPerJob, n, workers) {
    const w = cores(workers);
    return msPerJob * Math.ceil(n / w);
  }

  Object.assign(Pool, { run, available, cores, estimate, WorkerRuntime, sources });
  window.Pool = Pool;
})();
