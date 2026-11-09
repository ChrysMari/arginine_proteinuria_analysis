# Proteomics input data

The processed Spectronaut protein-group export files used by the proteomics analysis are available from PRIDE:

**PRIDE accession:** `PXD075940`

Download the required files and place them in this directory:

```text
data/proteomics/
```

The proteomics scripts expect files with the following stems:

- `20240714_ArgDiet_2w_Liver_Double_Labelled_Channels_report`
- `20240718_ArgDiet_2w_Lung_Double_Labelled_Channels_report`
- `20250430_ArgDiet_2w_Brain_Double_Labelled_Channels_report`
- `20250430_ArgDiet_2w_Colon_Double_Labelled_Channels_report`
- `20250430_ArgDiet_2w_Cortex_Double_Labelled_Channels_report`
- `20250430_ArgDiet_2w_Duodenum_Double_Labelled_Channels_report`
- `20250430_ArgDiet_2w_Erythrocytes_Double_Labelled_Channels_report`
- `20250430_ArgDiet_2w_Heart_repet_Double_Labelled_Channels_report`
- `20250430_ArgDiet_2w_Ileum_Double_Labelled_Channels_report`
- `20250430_ArgDiet_2w_Jejunum_Double_Labelled_Channels_report`
- `20250430_ArgDiet_2w_Medulla_Double_Labelled_Channels_report`
- `20250430_ArgDiet_2w_Muscle_Double_Labelled_Channels_report`
- `20250430_ArgDiet_2w_Plasma_Double_Labelled_Channels_report`
- `20250430_ArgDiet_2w_Spleen_Double_Labelled_Channels_report`
- `20250430_ArgDiet_2w_Urine_Double_Labelled_Channels_report`

The scripts accept the deposited `.csv` and `.tsv` report files corresponding to these stems.

These input files are used by:

- `proteomics/systemwide_protein_channel_analysis.py`
- `proteomics/supplementary_table_4_protein_abundance.py`

The deposited data files themselves are not included in this GitHub repository.
