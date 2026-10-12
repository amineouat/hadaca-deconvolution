# HADACA – Tumor Microenvironment Deconvolution

Solution to the [HADACA 3 health data challenge](https://www.codabench.org/competitions/11973) on Codabench.

The goal is to estimate the cellular composition of **pancreatic cancer (PDAC)** samples from bulk
**RNA-seq** and **DNA methylation** profiles. For each sample, the method predicts the proportions of
5 cell populations: endothelial cells, cancer-associated fibroblasts, immune cells, and two cancer
subtypes (classic and basal PDAC).

**Result: median score of 0.90** across the 9 evaluation datasets.

## Data

The challenge provides 9 heterogeneous datasets (234 samples in total: in silico, in vitro and in vivo
mixtures) and a reference file with the molecular profile of each pure cell type.

The data is not included in this repository. Download the starting kit from the
[challenge page](https://www.codabench.org/competitions/11973) and place its `data/` folder at the root
of this repository:

```
data/
├── reference_pdac.h5
├── mixes_SBN5.h5
├── ...
└── mixes_VIVO.h5
```

## Method

The deconvolution relies on **Non-Negative Least Squares (NNLS)** against reference profiles: each
mixture is modelled as a non-negative combination of the pure cell-type profiles, and the coefficients
are rescaled to sum to 1.

1. **Feature selection.** Genes (RNA) or CpG sites (methylation) are ranked with a composite score
   (40% variance, 40% fold-change, 20% coefficient of variation). Cell-type-specific markers are added,
   then a correlation filter removes redundant features.
2. **Dataset-adaptive modality.** DNA methylation is the primary signal for most datasets. RNA-seq is used
   for the in vivo dataset (VIVO) and combined with methylation (90% RNA / 10% methylation) for the
   in vitro dataset (VITR). The number of markers is adapted for some datasets (SBN5, SDEL).
3. **Post-processing.**
   - Per-dataset bias correction of the methylation-based estimates.
   - Mutual exclusivity between the two cancer subtypes (classic / basal) on VITR and SBN5.
   - SDN4 contains no basal cells, so its deconvolution is re-run without the basal reference.

The bias corrections and marker counts were tuned on the challenge datasets themselves, so the score
reflects performance on these datasets rather than on unseen data.

## Results

The score is a weighted geometric mean of Pearson/Spearman correlations (1/3), MAE and RMSE (1/3),
and the Aitchison distance (1/3). The final ranking uses the median over the 9 datasets.

| Dataset | Score |
|---------|:-----:|
| SDN4 | 0.957 |
| SDN5 | 0.922 |
| SDN6 | 0.916 |
| SDE5 | 0.912 |
| **SDC5 (median)** | **0.900** |
| SDEL | 0.898 |
| SBN5 | 0.862 |
| VIVO | 0.810 |
| VITR | 0.796 |

## Usage

```bash
pip install -r requirements.txt
python deconvolution.py
```

Predictions are written to `predictions/prediction.h5`.

## Repository structure

```
├── deconvolution.py        # Deconvolution method and entry point
├── data_processing.py      # HDF5 utilities provided by the challenge organizers
├── notebooks/
│   └── data_exploration.ipynb
└── requirements.txt
```

`data_processing.py` comes from the challenge starting kit
([bcm-uga/hadaca3](https://github.com/bcm-uga/hadaca3)).
