# Machine-Learning

## Peroxisome protein localization (`Peroxisome_protein_localization/`)

Early work toward a support vector machine (SVM) that predicts whether a plant protein goes to the
peroxisome, using its C-terminal sequence (where the PTS1 targeting signal sits). Each protein is turned
into a PSI-BLAST position-specific scoring matrix (PSSM) to use as features.

**Status:** work in progress. Several scripts are drafts and do not run as they are (see the notes).

| File | What it does |
|---|---|
| `form_neg_peroxisome_fasta.ipynb` | Builds the training sets. Removes known peroxisome and glyoxysome proteins (UniProt and GO lists) from reviewed plant UniProt proteins to get negatives, trims sequences to the C-terminus (last 50 aa), and removes negatives whose last 15 aa match a positive (Reumann et al. supplemental datasets) |
| `pssm-object.py` | `Generate_pssm` class: writes each FASTA record to a file, runs `psiblast` (3 iterations, E-value 0.001) against a local NR database, and reads the ASCII PSSM |
| `gpssm.py` | PSSM helpers: parses PSSM files, normalises scores with the logistic function 1 / (1 + e<sup>−x</sup>), and keeps track of sequences with no hits |
| `dbconn.py` | Saves protein name, last 15 aa, PSSM table and class label in a SQLite database (`pssmdb.sqlite`) |
| `model.py` | Start of the SVM model: reads the Reumann dataset (C-terminal 14 aa and class) for scikit-learn `SVC` / `LinearSVC` |
| `Peroxisome_svm_model/object-jupyter.ipynb` | Notebook version of the PSSM generator |
| `gpssm_deneme.py` | Scratch test of the class (*deneme* = "test") |

Requirements: Python 3, Biopython, NumPy, pandas, scikit-learn, BLAST+ (`psiblast`) and a protein database
(NR or UniProt).

### Notes

- `Bio.Alphabet` was removed in Biopython 1.78. Use an older Biopython or delete those imports.
- `model.py` has a syntax error on its last line, and `gpssm_deneme.py` uses `self` outside a method.
- Database and output paths (`/media/yunus/…`, `/home/yunus/…`) are local and must be changed.
- `a` is an empty placeholder file.
