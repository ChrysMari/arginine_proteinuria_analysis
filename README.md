# Arginine proteinuria analysis

Reproducible analysis code accompanying the manuscript:

**“Dietary isotope-labeling of Arginine reveals interorgan nitrogen shifts in proteinuria.”**

This repository contains the custom Python analyses used downstream of the deposited mass-spectrometry data. Spectronaut was used for proteomics channel assignment; the Python scripts provided here perform downstream data aggregation, statistical testing, multiple-testing correction, visualization, and supplementary-table generation.

## Repository structure

```text
arginine_proteinuria_analysis/
├── proteomics/
│   ├── systemwide_protein_channel_analysis.py
│   └── supplementary_table_4_protein_abundance.py
├── metabolomics/
│   ├── isotopologue_heatmap_analysis.py
│   └── unlabeled_metabolite_abundance_analysis.py
├── data/
│   ├── proteomics/
│   │   └── README.md
│   └── metabolomics/
│       └── README.md
├── results/
│   └── .gitkeep
├── requirements.txt
└── .gitignore
```

## Scripts

### `proteomics/systemwide_protein_channel_analysis.py`

Processes Spectronaut protein-group exports across tissues and summarizes protein abundance by MS2 channel and experimental group. The script produces absolute and percentage channel-abundance summaries and the corresponding figures.

Channel interpretation:

- `PG.MS2Channel1`: unlabeled
- `PG.MS2Channel2`: Arg10 / Pro6
- `PG.MS2Channel3`: Arg7 / Pro6

For manuscript significance testing, protein abundance is summed per biological replicate. Control and Proteinuric groups are compared using a two-sided Welch's t-test separately within ND and AD conditions. Benjamini–Hochberg correction is applied across tissues and both labeled channels separately within each dietary state.

### `proteomics/supplementary_table_4_protein_abundance.py`

Reproduces Supplementary Table 4 using the same per-mouse aggregation and statistical procedure as the manuscript proteomics analysis:

- two-sided Welch's t-test;
- Control vs Proteinuric comparisons separately within ND and AD;
- Benjamini–Hochberg correction across all tissue × labeled-channel comparisons separately within ND and AD.

The output is an Excel workbook containing the protein-abundance summaries and associated statistics.

### `metabolomics/isotopologue_heatmap_analysis.py`

Analyzes isotope/light ratios for labeled metabolite isotopologues across tissues in the 2-week heavy-arginine diet experiment.

For each isotopologue × tissue comparison, the script calculates:

- mean isotope/light ratios;
- Proteinuric/Control ratio;
- log2 Proteinuric/Control effect size;
- bootstrap 95% confidence interval;
- two-sided Mann–Whitney U test;
- global Benjamini–Hochberg FDR;
- Benjamini–Hochberg FDR within each parent metabolite.

The script exports significance-display versions based on:

- nominal `P < 0.05`;
- global `BH-FDR < 0.20`;
- within-metabolite `BH-FDR < 0.05`.

The manuscript figure uses the **global BH-FDR < 0.20** significance criterion.

### `metabolomics/unlabeled_metabolite_abundance_analysis.py`

Analyzes unlabeled metabolite abundance (`IsotopeLabelType == 0`) in control and proteinuric mice across tissues.

For each metabolite × tissue comparison, the script calculates:

- group mean abundance;
- Proteinuric/Control ratio;
- log10 Proteinuric/Control effect size;
- bootstrap 95% confidence interval;
- two-sided Mann–Whitney U test;
- nominal P value;
- global Benjamini–Hochberg FDR for reference;
- Benjamini–Hochberg FDR within each metabolite.

Two significance-display versions are exported:

- nominal `P < 0.05`;
- within-metabolite `BH-FDR < 0.20`.

The manuscript figure uses the **within-metabolite BH-FDR < 0.20** significance criterion.

## Data

The mass-spectrometry data and processed input files used by these scripts are deposited in public repositories and are not duplicated in this GitHub repository.

Before running the code, download the relevant processed input files and place them under:

```text
data/proteomics/
data/metabolomics/
```

The exact filenames required by the scripts are listed in the README file within each data directory.

### Public data repositories

**Proteomics**

PRIDE accession: `PXD075940`

**Metabolomics**

MassIVE accessions:

- `MSV000101202`
- `MSV000101204`

## Installation

Python 3.9 or later is recommended.

Create a virtual environment:

```bash
python -m venv .venv
```

On Windows:

```bash
.venv\Scripts\activate
pip install -r requirements.txt
```

## Running the analyses

Run all commands from the repository root.

### Proteomics

```bash
python proteomics/systemwide_protein_channel_analysis.py
python proteomics/supplementary_table_4_protein_abundance.py
```

### Metabolomics

```bash
python metabolomics/isotopologue_heatmap_analysis.py
python metabolomics/unlabeled_metabolite_abundance_analysis.py
```

Outputs are written under the `results/` directory.

## Reproducibility notes

The metabolomics scripts use fixed random seeds for bootstrap confidence intervals. Input filenames are explicit so that each analysis can be matched directly to the publicly deposited source files.

The scripts in this repository correspond to the analyses used for manuscript revision and provide a versioned record of the derivation of the reported quantities, figures, and supplementary tables.

## License

This code is released under the MIT License. See `LICENSE` for details.
