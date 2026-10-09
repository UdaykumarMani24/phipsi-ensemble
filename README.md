# Phi/Psi Ensemble Analyzer

**Developed by Dr. Udayakumar Mani and Dr. Senthilkumar Rathnasamy**

A free web and Android app that measures **per-residue backbone flexibility** from NMR ensembles or several crystal structures of the same protein, and compares it with **AlphaFold pLDDT** confidence.

- Computes backbone torsion angles φ and ψ for every residue in every model
- Measures their spread with **circular statistics** (angles wrap at ±180°)
- Flags **two-state residues** (for example peptide flips) where a single average is misleading
- Compares flexibility with AlphaFold pLDDT (Spearman ρ with a permutation test)
- Interactive flexibility profile, 3D view, Ramachandran plot, scatter plot and CSV export
- Runs entirely on your device: no server, works offline once installed

See **[How it works](how-it-works.html)** for the method, and **[DEPLOY.md](DEPLOY.md)** to publish on GitHub and Android.

## Use it

Open the app at `https://<your-github-username>.github.io/phipsi-ensemble/` (after setup below), or open `index.html` locally.

Example: PDB ID `1D3Z`, chain `A`, UniProt `P62987`, offset `0` (ubiquitin, 10 NMR models).

| Input | What it is |
|---|---|
| PDB ID(s) | One NMR entry (all models used), or several IDs of the same protein separated by commas |
| Chain | Optional; needed when comparing several entries |
| UniProt accession | Optional; loads the AlphaFold model for comparison |
| Residue offset | UniProt position = PDB residue number + offset |

## Files

```
index.html            the app
how-it-works.html     method explained
js/core.js            all calculations (no page code; testable with Node)
js/app.js             page logic and plots
css/style.css         styles
vendor/3Dmol-min.js   3D viewer (3Dmol.js, BSD-3-Clause)
manifest.webmanifest, sw.js, icons/   make it installable as an app
tests/                validation scripts
```

## Validation

```bash
cd tests
python3 make_test_ensemble.py      # builds a test ensemble from KNOWN angles
node test_core.js                  # checks the app recovers them (largest error 0.13°)
pip install biopython
python3 compare_biopython.py structure.pdb phi_psi_per_model.csv   # independent check on real data
```

## Cite

If you use this tool, please cite it (see `CITATION.cff`):

> Mani U., Rathnasamy S. *Phi/Psi Ensemble Analyzer*. Software, 2026.

## License

MIT License (see `LICENSE`). Uses 3Dmol.js (BSD-3-Clause). Structure data: RCSB PDB. Predictions: AlphaFold Protein Structure Database (CC-BY 4.0).
