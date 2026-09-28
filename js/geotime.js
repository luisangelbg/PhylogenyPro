/* PhylogenyPro — molecular phylogenetics in the browser.
   Copyright (C) 2026  Luis Ángel Barrera-Guzmán

   This program is free software: you can redistribute it and/or modify it under
   the terms of the GNU General Public License as published by the Free Software
   Foundation, either version 3 of the License, or (at your option) any later
   version. It is distributed in the hope that it will be useful, but WITHOUT ANY
   WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
   A PARTICULAR PURPOSE. See the GNU General Public License, in the file LICENSE
   at the root of this program, or <https://www.gnu.org/licenses/>. */

/* PhylogenyPro — geological time scale for dated trees and fossil calibrations.

   Ages (Ma, base → top) and colours follow the International Chronostratigraphic
   Chart of the International Commission on Stratigraphy, v2020/03 (Cohen et al.
   2013, updated). The ICS publishes the chart for free use with that citation,
   which is given here and on the home page; the ages are measurements and the
   colours the standard stratigraphic ones. See PROCEDENCIA.md.
   Stages are listed for the Cenozoic and the Cretaceous, the
   interval in which almost every angiosperm fossil calibration falls. Always
   check a boundary against the current chart before publishing a calibration. */

(function () {

  const CHART = 'ICS International Chronostratigraphic Chart v2020/03';

  /* [name, rank, base, top, [r,g,b], parent] */
  const UNITS = [
    /* eras */
    ['Cenozoic', 'era', 66.0, 0, [242, 249, 29], null],
    ['Mesozoic', 'era', 251.902, 66.0, [103, 197, 202], null],
    ['Paleozoic', 'era', 541.0, 251.902, [153, 192, 141], null],
    /* periods */
    ['Quaternary', 'period', 2.58, 0, [249, 249, 127], 'Cenozoic'],
    ['Neogene', 'period', 23.03, 2.58, [255, 230, 25], 'Cenozoic'],
    ['Paleogene', 'period', 66.0, 23.03, [253, 154, 82], 'Cenozoic'],
    ['Cretaceous', 'period', 145.0, 66.0, [127, 198, 78], 'Mesozoic'],
    ['Jurassic', 'period', 201.3, 145.0, [52, 178, 201], 'Mesozoic'],
    ['Triassic', 'period', 251.902, 201.3, [129, 43, 146], 'Mesozoic'],
    ['Permian', 'period', 298.9, 251.902, [240, 64, 40], 'Paleozoic'],
    ['Carboniferous', 'period', 358.9, 298.9, [103, 165, 153], 'Paleozoic'],
    ['Devonian', 'period', 419.2, 358.9, [203, 140, 55], 'Paleozoic'],
    ['Silurian', 'period', 443.8, 419.2, [179, 225, 182], 'Paleozoic'],
    ['Ordovician', 'period', 485.4, 443.8, [0, 146, 112], 'Paleozoic'],
    ['Cambrian', 'period', 541.0, 485.4, [127, 160, 86], 'Paleozoic'],
    /* epochs */
    ['Holocene', 'epoch', 0.0117, 0, [254, 242, 224], 'Quaternary'],
    ['Pleistocene', 'epoch', 2.58, 0.0117, [255, 242, 174], 'Quaternary'],
    ['Pliocene', 'epoch', 5.333, 2.58, [255, 255, 153], 'Neogene'],
    ['Miocene', 'epoch', 23.03, 5.333, [255, 255, 0], 'Neogene'],
    ['Oligocene', 'epoch', 33.9, 23.03, [254, 192, 122], 'Paleogene'],
    ['Eocene', 'epoch', 56.0, 33.9, [253, 180, 108], 'Paleogene'],
    ['Paleocene', 'epoch', 66.0, 56.0, [253, 167, 95], 'Paleogene'],
    ['Late Cretaceous', 'epoch', 100.5, 66.0, [166, 216, 74], 'Cretaceous'],
    ['Early Cretaceous', 'epoch', 145.0, 100.5, [140, 205, 87], 'Cretaceous'],
    ['Late Jurassic', 'epoch', 163.5, 145.0, [179, 227, 238], 'Jurassic'],
    ['Middle Jurassic', 'epoch', 174.1, 163.5, [128, 207, 216], 'Jurassic'],
    ['Early Jurassic', 'epoch', 201.3, 174.1, [66, 174, 208], 'Jurassic'],
    ['Late Triassic', 'epoch', 237.0, 201.3, [189, 140, 195], 'Triassic'],
    ['Middle Triassic', 'epoch', 247.2, 237.0, [177, 104, 177], 'Triassic'],
    ['Early Triassic', 'epoch', 251.902, 247.2, [152, 57, 153], 'Triassic'],
    /* stages — Cenozoic */
    ['Meghalayan', 'stage', 0.0042, 0, [253, 237, 236], 'Holocene'],
    ['Northgrippian', 'stage', 0.0082, 0.0042, [253, 236, 228], 'Holocene'],
    ['Greenlandian', 'stage', 0.0117, 0.0082, [254, 236, 219], 'Holocene'],
    ['Upper Pleistocene', 'stage', 0.129, 0.0117, [255, 242, 211], 'Pleistocene'],
    ['Chibanian', 'stage', 0.774, 0.129, [255, 242, 199], 'Pleistocene'],
    ['Calabrian', 'stage', 1.80, 0.774, [255, 242, 186], 'Pleistocene'],
    ['Gelasian', 'stage', 2.58, 1.80, [255, 237, 179], 'Pleistocene'],
    ['Piacenzian', 'stage', 3.600, 2.58, [255, 255, 191], 'Pliocene'],
    ['Zanclean', 'stage', 5.333, 3.600, [255, 255, 179], 'Pliocene'],
    ['Messinian', 'stage', 7.246, 5.333, [255, 255, 115], 'Miocene'],
    ['Tortonian', 'stage', 11.63, 7.246, [255, 255, 102], 'Miocene'],
    ['Serravallian', 'stage', 13.82, 11.63, [255, 255, 89], 'Miocene'],
    ['Langhian', 'stage', 15.97, 13.82, [255, 255, 77], 'Miocene'],
    ['Burdigalian', 'stage', 20.44, 15.97, [255, 255, 65], 'Miocene'],
    ['Aquitanian', 'stage', 23.03, 20.44, [255, 255, 51], 'Miocene'],
    ['Chattian', 'stage', 27.82, 23.03, [254, 230, 170], 'Oligocene'],
    ['Rupelian', 'stage', 33.9, 27.82, [254, 217, 154], 'Oligocene'],
    ['Priabonian', 'stage', 37.71, 33.9, [253, 205, 161], 'Eocene'],
    ['Bartonian', 'stage', 41.2, 37.71, [253, 192, 145], 'Eocene'],
    ['Lutetian', 'stage', 47.8, 41.2, [252, 180, 130], 'Eocene'],
    ['Ypresian', 'stage', 56.0, 47.8, [252, 167, 115], 'Eocene'],
    ['Thanetian', 'stage', 59.2, 56.0, [253, 191, 111], 'Paleocene'],
    ['Selandian', 'stage', 61.6, 59.2, [254, 191, 101], 'Paleocene'],
    ['Danian', 'stage', 66.0, 61.6, [253, 180, 98], 'Paleocene'],
    /* stages — Cretaceous */
    ['Maastrichtian', 'stage', 72.1, 66.0, [242, 250, 140], 'Late Cretaceous'],
    ['Campanian', 'stage', 83.6, 72.1, [230, 244, 127], 'Late Cretaceous'],
    ['Santonian', 'stage', 86.3, 83.6, [217, 239, 116], 'Late Cretaceous'],
    ['Coniacian', 'stage', 89.8, 86.3, [204, 233, 104], 'Late Cretaceous'],
    ['Turonian', 'stage', 93.9, 89.8, [191, 227, 93], 'Late Cretaceous'],
    ['Cenomanian', 'stage', 100.5, 93.9, [179, 222, 83], 'Late Cretaceous'],
    ['Albian', 'stage', 113.0, 100.5, [204, 234, 151], 'Early Cretaceous'],
    ['Aptian', 'stage', 125.0, 113.0, [191, 228, 138], 'Early Cretaceous'],
    ['Barremian', 'stage', 129.4, 125.0, [179, 223, 127], 'Early Cretaceous'],
    ['Hauterivian', 'stage', 132.6, 129.4, [166, 217, 117], 'Early Cretaceous'],
    ['Valanginian', 'stage', 139.8, 132.6, [153, 211, 106], 'Early Cretaceous'],
    ['Berriasian', 'stage', 145.0, 139.8, [140, 205, 96], 'Early Cretaceous'],
  ].map(([name, rank, base, top, c, parent]) => ({ name, rank, base, top, color: `rgb(${c[0]},${c[1]},${c[2]})`, parent }));

  const RANKS = ['era', 'period', 'epoch', 'stage'];

  function byName(name) {
    const k = String(name || '').trim().toLowerCase().replace(/\s+/g, ' ');
    return UNITS.find(u => u.name.toLowerCase() === k) || null;
  }

  /* units of a rank overlapping [young, old] Ma */
  function within(rank, young, old) {
    return UNITS.filter(u => u.rank === rank && u.base > young && u.top < old);
  }

  /* the finest unit containing an age */
  function at(age, rank) {
    const ranks = rank ? [rank] : RANKS.slice().reverse();
    for (const r of ranks) {
      const u = UNITS.find(x => x.rank === r && age <= x.base && age >= x.top);
      if (u) return u;
    }
    return null;
  }

  /* "Late Miocene" style descriptions of an age, for texts */
  function describe(age) {
    const e = at(age, 'epoch'), p = at(age, 'period');
    if (!e && !p) return age.toFixed(1) + ' Ma';
    return e ? e.name + (e.name.includes(p ? p.name : '#') ? '' : '') : p.name;
  }

  window.GeoTime = { CHART, UNITS, RANKS, byName, within, at, describe };
})();
