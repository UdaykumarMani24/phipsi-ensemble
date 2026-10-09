/* =========================================================================
   Phi/Psi Ensemble Analyzer - core calculations
   Developed by Dr. Udayakumar Mani and Dr. Senthilkumar Rathnasamy
   Pure functions only (no page code), so they can be tested with Node.
   ========================================================================= */
(function (root, factory) {
  if (typeof module === 'object' && module.exports) module.exports = factory();
  else root.PhiPsiCore = factory();
})(typeof self !== 'undefined' ? self : this, function () {
  'use strict';

  const DEG = 180 / Math.PI;

  /* ---------- vector maths ---------- */
  const sub = (a, b) => [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
  const dot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
  const cross = (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
  const norm = v => Math.sqrt(dot(v, v));

  /** Signed dihedral angle a-b-c-d in degrees (-180..180), IUPAC sign convention. */
  function dihedral(a, b, c, d) {
    const b1 = sub(b, a), b2 = sub(c, b), b3 = sub(d, c);
    const n1 = cross(b1, b2), n2 = cross(b2, b3);
    return Math.atan2(norm(b2) * dot(b1, n2), dot(n1, n2)) * DEG;
  }

  /** Parse every MODEL of a PDB-format text. Keeps backbone N, CA, C only. */
  function parseModels(text, chainFilter = '', maxModels = 500) {
    const models = [];
    let cur = [], lastKey = null;
    for (const line of text.split(/\r\n|\r|\n/)) {
      const rec = line.slice(0, 6);
      if (rec === 'ENDMDL') {
        if (cur.length) models.push(cur);
        cur = []; lastKey = null;
        if (models.length >= maxModels) return models;
        continue;
      }
      if (rec !== 'ATOM  ' || line.length < 54) continue;
      const name = line.slice(12, 16).trim();
      if (name !== 'N' && name !== 'CA' && name !== 'C') continue;
      const alt = line[16];
      if (alt !== ' ' && alt !== 'A') continue;
      const chain = line[21];
      if (chainFilter && chain !== chainFilter) continue;
      const resSeq = parseInt(line.slice(22, 26), 10);
      const icode = (line[26] || ' ').trim();
      const key = chain + '|' + resSeq + '|' + icode;
      if (key !== lastKey) {
        cur.push({ chain, resSeq, icode, name: line.slice(17, 20).trim(), atoms: {} });
        lastKey = key;
      }
      cur[cur.length - 1].atoms[name] = [
        parseFloat(line.slice(30, 38)), parseFloat(line.slice(38, 46)), parseFloat(line.slice(46, 54)),
      ];
    }
    if (cur.length) models.push(cur);
    return models;
  }

  const r1 = v => Math.round(v * 10) / 10;

  /** phi/psi of every residue in one model. Chain breaks (C-N > 2 A) give null. */
  function modelTorsions(res) {
    const out = {};
    for (let i = 0; i < res.length; i++) {
      const r = res[i], A = r.atoms;
      if (!A.N || !A.CA || !A.C) continue;
      let phi = null, psi = null;
      const p = res[i - 1], n = res[i + 1];
      if (p && p.chain === r.chain && p.atoms.C && norm(sub(p.atoms.C, A.N)) < 2.0)
        phi = r1(dihedral(p.atoms.C, A.N, A.CA, A.C));
      if (n && n.chain === r.chain && n.atoms.N && norm(sub(A.C, n.atoms.N)) < 2.0)
        psi = r1(dihedral(A.N, A.CA, A.C, n.atoms.N));
      out[r.chain + '|' + r.resSeq + '|' + r.icode] =
        { chain: r.chain, resSeq: r.resSeq, icode: r.icode, name: r.name, phi, psi };
    }
    return out;
  }

  /** Circular mean and circular SD (degrees). SD is null when fewer than 2 values. */
  function circStats(vals) {
    let s = 0, c = 0, n = 0;
    for (const v of vals) {
      if (v === null || v === undefined) continue;
      s += Math.sin(v / DEG); c += Math.cos(v / DEG); n++;
    }
    if (!n) return { n: 0, mean: null, sd: null };
    s /= n; c /= n;
    const R = Math.min(1, Math.hypot(s, c));
    const mean = r1(Math.atan2(s, c) * DEG);
    if (n < 2) return { n, mean, sd: null };
    const sd = R > 1e-9 ? Math.min(180, Math.sqrt(Math.max(0, -2 * Math.log(R))) * DEG) : 180;
    return { n, mean, sd: r1(sd) };
  }

  const angDiff = (a, b) => { const d = Math.abs(a - b) % 360; return d > 180 ? 360 - d : d; };

  /**
   * Flags residues whose (phi, psi) points form two separate clusters
   * (for example a peptide flip). Simple 2-means on the circle.
   * Returns null, or {sizes:[a,b], centres:[[phi,psi],[phi,psi]], sep}.
   */
  function twoStates(phi, psi) {
    const pts = [];
    for (let i = 0; i < phi.length; i++) if (phi[i] !== null && psi[i] !== null) pts.push([phi[i], psi[i]]);
    if (pts.length < 6) return null;
    const dist = (a, b) => Math.hypot(angDiff(a[0], b[0]), angDiff(a[1], b[1]));
    // start from the two points farthest apart
    let A = pts[0], B = pts[0], best = -1;
    for (const p of pts) { const d = dist(pts[0], p); if (d > best) { best = d; B = p; } }
    best = -1;
    for (const p of pts) { const d = dist(B, p); if (d > best) { best = d; A = p; } }
    let ca = A.slice(), cb = B.slice(), lab = [];
    const centre = g => {
      const f = k => circStats(g.map(p => p[k])).mean;
      return [f(0), f(1)];
    };
    for (let it = 0; it < 20; it++) {
      lab = pts.map(p => (dist(p, ca) <= dist(p, cb) ? 0 : 1));
      const g0 = pts.filter((_, i) => lab[i] === 0), g1 = pts.filter((_, i) => lab[i] === 1);
      if (!g0.length || !g1.length) return null;
      ca = centre(g0); cb = centre(g1);
    }
    const g0 = pts.filter((_, i) => lab[i] === 0), g1 = pts.filter((_, i) => lab[i] === 1);
    const minFrac = Math.min(g0.length, g1.length) / pts.length;
    const sep = dist(ca, cb);
    const spread = g => Math.sqrt(g.reduce((s, p) => s + dist(p, centre(g)) ** 2, 0) / g.length);
    const within = Math.max(spread(g0), spread(g1));
    if (minFrac >= 0.15 && sep >= 60 && within < sep / 2.5)
      return { sizes: [g0.length, g1.length], centres: [ca, cb], sep: r1(sep) };
    return null;
  }

  /** Combine per-model torsions into one row per residue. */
  function aggregate(allModels, plddtMap = {}, offset = 0) {
    const meta = new Map();
    for (const t of allModels) for (const k in t) if (!meta.has(k)) meta.set(k, t[k]);
    const rows = [];
    for (const [k, m] of meta) {
      const phi = allModels.map(t => (t[k] ? t[k].phi : null));
      const psi = allModels.map(t => (t[k] ? t[k].psi : null));
      const P = circStats(phi), Q = circStats(psi);
      const sds = [P.sd, Q.sd].filter(v => v !== null);
      const pl = m.icode === '' ? plddtMap[m.resSeq + offset] : undefined;
      rows.push({
        chain: m.chain, resSeq: m.resSeq, icode: m.icode, name: m.name, phi, psi,
        n: Math.max(P.n, Q.n),
        phiMean: P.mean, phiSD: P.sd, psiMean: Q.mean, psiSD: Q.sd,
        flex: sds.length ? r1(sds.reduce((a, b) => a + b, 0) / sds.length) : null,
        plddt: pl === undefined ? null : r1(pl),
        twoState: twoStates(phi, psi),
      });
    }
    return rows;
  }

  /** AlphaFold pLDDT lives in the B-factor column of CA atoms. */
  function parsePlddt(text) {
    const map = {};
    for (const line of text.split(/\r\n|\r|\n/)) {
      if (!line.startsWith('ATOM  ') || line.length < 66) continue;
      if (line.slice(12, 16).trim() !== 'CA') continue;
      map[parseInt(line.slice(22, 26), 10)] = parseFloat(line.slice(60, 66));
    }
    return map;
  }

  /* ---------- Spearman correlation with permutation p-value ---------- */
  function ranks(v) {
    const idx = v.map((_, i) => i).sort((a, b) => v[a] - v[b]);
    const r = new Array(v.length);
    for (let i = 0; i < idx.length;) {
      let j = i;
      while (j + 1 < idx.length && v[idx[j + 1]] === v[idx[i]]) j++;
      for (let k = i; k <= j; k++) r[idx[k]] = (i + j) / 2 + 1;
      i = j + 1;
    }
    return r;
  }
  function pearson(x, y) {
    const n = x.length, mx = x.reduce((a, b) => a + b) / n, my = y.reduce((a, b) => a + b) / n;
    let sxy = 0, sxx = 0, syy = 0;
    for (let i = 0; i < n; i++) { const dx = x[i] - mx, dy = y[i] - my; sxy += dx * dy; sxx += dx * dx; syy += dy * dy; }
    const d = Math.sqrt(sxx * syy);
    return d > 0 ? sxy / d : null;
  }
  function seededRandom(seed) {   // mulberry32, so p-values are reproducible
    return function () {
      seed |= 0; seed = (seed + 0x6D2B79F5) | 0;
      let t = Math.imul(seed ^ (seed >>> 15), 1 | seed);
      t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
      return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
    };
  }
  function spearman(x, y, perms = 5000) {
    const n = x.length;
    if (n < 5) return null;
    const rx = ranks(x), ry = ranks(y).slice();
    const rho = pearson(rx, ry);
    if (rho === null) return null;
    const rand = seededRandom(12345);
    let ge = 0;
    for (let p = 0; p < perms; p++) {
      for (let i = n - 1; i > 0; i--) { const j = Math.floor(rand() * (i + 1)); [ry[i], ry[j]] = [ry[j], ry[i]]; }
      if (Math.abs(pearson(rx, ry)) >= Math.abs(rho) - 1e-12) ge++;
    }
    return { rho: Math.round(rho * 1000) / 1000, p: (ge + 1) / (perms + 1), n };
  }

  return { dihedral, parseModels, modelTorsions, circStats, twoStates, aggregate, parsePlddt, spearman, angDiff };
});
