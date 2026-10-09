// Checks js/core.js against angles that were used to BUILD the test ensemble.
// Run from the tests folder:  python3 make_test_ensemble.py && node test_core.js
const fs = require('fs');
const path = require('path');
const C = require(path.join(__dirname, '..', 'js', 'core.js'));

const pdb = fs.readFileSync(path.join(__dirname, 'test_ensemble.pdb'), 'utf8');
const expected = fs.readFileSync(path.join(__dirname, 'expected_angles.csv'), 'utf8').trim().split('\n').slice(1).map(l => l.split(','));

const models = C.parseModels(pdb).map(C.modelTorsions);
let worst = 0, checked = 0;
for (const [m, res, phi, psi] of expected) {
  const t = models[+m - 1]['A|' + res + '|'];
  for (const [want, got] of [[phi, t.phi], [psi, t.psi]]) {
    if (want === '') { if (got !== null) throw new Error('expected null at model ' + m + ' res ' + res); continue; }
    worst = Math.max(worst, C.angDiff(+want, got)); checked++;
  }
}
console.log('angles checked:', checked, ' worst error (deg):', worst.toFixed(3));
if (worst > 0.15) throw new Error('angle error too large');

const rows = C.aggregate(models);
const two = rows.filter(r => r.twoState).map(r => r.resSeq);
console.log('two-state residues found:', two.join(', ') || 'none');
if (two.join() !== '5') throw new Error('expected residue 5 only to be two-state');

const flex = Object.fromEntries(rows.map(r => [r.resSeq, r.flex]));
console.log('flex rigid res 3:', flex[3], ' tail res 12:', flex[12]);
if (!(flex[3] < 3 && flex[12] > 40)) throw new Error('flexibility profile wrong');

// circular stats across the +/-180 seam: 170 and -170 must average to 180, SD small
const cs = C.circStats([170, -170]);
console.log('circular mean of 170 & -170:', cs.mean, ' SD:', cs.sd);
if (Math.abs(Math.abs(cs.mean) - 180) > 0.1 || cs.sd > 15) throw new Error('circular stats wrong');
if (C.circStats([42]).sd !== null) throw new Error('single value must give null SD');

const s = C.spearman([1, 2, 3, 4, 5, 6, 7, 8], [8, 7, 6, 5, 4, 3, 2, 1]);
console.log('spearman perfect negative:', s.rho, ' p:', s.p.toFixed(4));
if (s.rho !== -1) throw new Error('spearman wrong');
console.log('ALL TESTS PASSED');
