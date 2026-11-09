# Metabolomics input data

The processed targeted-metabolomics CSV reports used by the metabolomics analyses are available from MassIVE:

**MassIVE accessions:**

- `MSV000101202`
- `MSV000101204`

Download the required files and place them in this directory:

```text
data/metabolomics/
```

Expected files:

- `AVRQ_ArgDiet_2weeks_Ctx_Med_Liver_Heart_report.csv`
- `AVRQ_ArgDiet_2weeks_Lung_Pancreas_Spleen_Muscle_report.csv`
- `AVRQ_ArgDiet_2weeks_Brain_WAT_Duodenum_Ileum_report.csv`
- `AVRQ_ArgDiet_2weeks_Jejenum_Colon_Erythrocytes_report.csv`
- `AVRQ_ArgDiet_2weeks_TCA_report.csv`

These files are used by:

- `metabolomics/unlabeled_metabolite_abundance_analysis.py` — uses the first four files.
- `metabolomics/isotopologue_heatmap_analysis.py` — uses all five files.

The deposited data files themselves are not included in this GitHub repository.
