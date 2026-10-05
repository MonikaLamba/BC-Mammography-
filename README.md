# Leakage-free 2-Stage GeS + HWNB: analysis code

Reproduces every number and figure in the conference paper.

## Order of execution
1. `export.R` – loads the Bioconductor packages breastCancerTRANSBIG, breastCancerUPP,
   breastCancerMAINZ and breastCancerUNT, maps probes to genes (highest-mean probe) and
   writes `<cohort>_expr.csv` and `<cohort>_pheno.csv`. Edit the paths at the top of
   `pipeline.py` (`DATA`) to point to these files.
2. `run_main.py <cohort>` – 10 x 5-fold leakage-free CV (MDL, SMOTE-N, CFS-BFS,
   CONSISTENCY-BFS, 11 classifiers), grouped permutation importance, KernelSHAP for the first
   repetition, and the leaky protocol of the IJCDS paper. Resumable: set `BUDGET=<seconds>`.
   Run for GSE7390, UPP, MAINZ and UNT.
3. `run_lodo.py <held-out cohort>` – leave-one-dataset-out (per-cohort z-scored genes).
4. `run_ablation.py` – selection-stage ablation on GSE7390 (25 folds).
5. `aggregate.py` – tables, Friedman test, stability, robust rank aggregation, SHAP and
   permutation summaries -> `summary.pkl`.
6. `survival.py` – genomic-grade score, published four genes, grade 2 re-stratification.
7. `xai_fig.py` – proliferation-module enrichment and figure data.

## Requirements
Python 3.12: numpy, pandas, scipy, numba, scikit-learn, imbalanced-learn, shap, lifelines,
scikit-survival, skrebate, wittgenstein. R with Biobase for the export step.

## Adding the five other cohorts of the IJCDS paper
Export each GEO series matrix to the same `<name>_expr.csv` / `<name>_pheno.csv` format
(columns: grade, er, node, size, age, t.rfs, e.rfs, t.dmfs, e.dmfs), add the cohort to
`COHORTS` in `pipeline.py`, and rerun steps 2-7.

## Refined model (Section 5.7)
8. `run_binary.py <cohort>` – grade 1 vs grade 3: 10 x 5 leakage-free CV of 2-Stage GeS + HWNB
   and the label-free proliferation score (AURKA module, defined inside each training fold);
   final model applied to the cohort's grade 2 tumours. Resumable with `BUDGET`.
9. `survival_refined.py` – LODO proliferation score, grade 2 re-stratification by all four
   splits (Table 16) and the continuous score as a prognostic factor.

## Jupyter notebooks
`notebooks/` holds one notebook per script, in the same run order (Steps 2–9). Put the CSV files
from `export.R` in `notebooks/data/` and `common_genes.csv` next to the notebooks. `core.ipynb` and
`pipeline.ipynb` are loaded by the others with `%run`, so keep all notebooks in one folder.
Notebooks that run per cohort have a **Parameters** cell at the top (`cohort` or `held_out`).
