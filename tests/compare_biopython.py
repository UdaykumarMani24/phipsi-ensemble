"""Independent check: phi/psi from Biopython vs the analyzer's per-model CSV.

Usage:
  python3 compare_biopython.py structure.pdb phi_psi_per_model.csv
(the CSV is the "Per-model torsions CSV" downloaded from the app for the same structure)
"""
import csv
import math
import sys

from Bio.PDB import PDBParser, PPBuilder


def main(pdb_path, csv_path):
    structure = PDBParser(QUIET=True).get_structure("s", pdb_path)
    ours = {}
    with open(csv_path) as f:
        for row in csv.DictReader(f):
            model_no = int(row["model"].split("#")[1]) if "#" in row["model"] else 1
            ours[(model_no, row["chain"], int(row["resSeq"]))] = (row["phi"], row["psi"])
    worst, n = 0.0, 0
    for mi, model in enumerate(structure, start=1):
        for pp in PPBuilder().build_peptides(model):
            for res, (phi, psi) in zip(pp, pp.get_phi_psi_list()):
                key = (mi, res.get_parent().id, res.id[1])
                if key not in ours:
                    continue
                for bio, mine in ((phi, ours[key][0]), (psi, ours[key][1])):
                    if bio is None or mine == "":
                        continue
                    d = abs(math.degrees(bio) - float(mine)) % 360
                    worst = max(worst, min(d, 360 - d))
                    n += 1
    print(f"compared {n} angles; largest difference = {worst:.3f} degrees")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
