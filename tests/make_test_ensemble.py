"""Builds a small synthetic NMR-style ensemble with KNOWN phi/psi angles.

Used to check the analyzer: every angle written here must be recovered.
  - residues 2-8 : rigid (tiny noise)
  - residue 5    : two states (half the models helical psi, half extended psi)
  - residues 11-14: floppy tail (random angles)
Run:  python3 make_test_ensemble.py  ->  writes test_ensemble.pdb and expected_angles.csv
"""
import csv
import math
import random

import numpy as np

N_RES, N_MODELS = 14, 10
BOND = {"N-CA": 1.458, "CA-C": 1.525, "C-N": 1.329}
ANGLE = {"N-CA-C": 111.2, "CA-C-N": 116.2, "C-N-CA": 121.7}
OMEGA = 180.0


def place(a, b, c, bond, angle, torsion):
    """NeRF: position of atom d given a, b, c, bond length c-d, angle b-c-d, torsion a-b-c-d (degrees)."""
    angle, torsion = math.radians(angle), math.radians(torsion)
    bc = (c - b) / np.linalg.norm(c - b)
    n = np.cross(b - a, bc)
    n /= np.linalg.norm(n)
    m = np.cross(n, bc)
    d2 = np.array([-bond * math.cos(angle), bond * math.sin(angle) * math.cos(torsion), bond * math.sin(angle) * math.sin(torsion)])
    return c + d2[0] * bc + d2[1] * m + d2[2] * n


def build(phis, psis):
    atoms = [np.array([0.0, 0.0, 0.0]), np.array([BOND["N-CA"], 0.0, 0.0])]
    a = ANGLE["N-CA-C"]
    atoms.append(atoms[1] + BOND["CA-C"] * np.array([-math.cos(math.radians(a)), math.sin(math.radians(a)), 0.0]))
    for i in range(1, N_RES):
        n_, ca, c = atoms[-3], atoms[-2], atoms[-1]
        N = place(n_, ca, c, BOND["C-N"], ANGLE["CA-C-N"], psis[i - 1])
        CA = place(ca, c, N, BOND["N-CA"], ANGLE["C-N-CA"], OMEGA)
        C = place(c, N, CA, BOND["CA-C"], ANGLE["N-CA-C"], phis[i])
        atoms += [N, CA, C]
    return atoms


def main():
    rng = random.Random(7)
    rows, lines = [], []
    for m in range(N_MODELS):
        phis, psis = [], []
        for r in range(N_RES):
            resno = r + 1
            if resno == 5:
                phi, psi = (-65, -40) if m % 2 == 0 else (-80, 135)
            elif resno >= 11:
                phi, psi = rng.uniform(-180, 180), rng.uniform(-180, 180)
            elif resno % 2:
                phi, psi = -63, -42
            else:
                phi, psi = -120, 130
            phis.append(phi + rng.gauss(0, 1.0))
            psis.append(psi + rng.gauss(0, 1.0))
        xyz = build(phis, psis)
        lines.append(f"MODEL     {m + 1:4d}")
        serial = 1
        for r in range(N_RES):
            for k, name in enumerate(["N", "CA", "C"]):
                x, y, z = xyz[3 * r + k]
                lines.append(f"ATOM  {serial:5d}  {name:<3s} ALA A{r + 1:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           {name[0]}")
                serial += 1
        lines.append("ENDMDL")
        for r in range(N_RES):
            rows.append([m + 1, r + 1, "" if r == 0 else round(phis[r], 3), "" if r == N_RES - 1 else round(psis[r], 3)])
    lines.append("END")
    with open("test_ensemble.pdb", "w") as f:
        f.write("\n".join(lines) + "\n")
    with open("expected_angles.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["model", "resSeq", "phi", "psi"])
        w.writerows(rows)


if __name__ == "__main__":
    main()
