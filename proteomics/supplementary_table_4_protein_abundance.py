"""
Generate Supplementary Table 4: system-wide protein abundance and isotope-channel statistics.

This script processes the 15 deposited Spectronaut protein-group exports used
for Supplementary Table 4. Files may be CSV or TSV. Sample names are mapped to
Control/Proteinuric and ND/AD groups using the study-specific SAMPLE_MAP below.

For each tissue, the script:
1. Retains the reported Spectronaut protein-group metrics.
2. Renames biological samples using the study-specific sample map.
3. Pivots protein-level metrics into the wide format used for Supplementary Table 4.
4. Sums PG.MS2Channel2 and PG.MS2Channel3 across protein groups per mouse.
5. Compares Control vs Proteinuric separately within ND and AD using two-sided
   Welch's t-tests.
6. Applies Benjamini-Hochberg FDR correction across all tissue x labeled-channel
   tests together, separately within ND and AD.
7. Writes one Excel workbook containing one sheet per tissue plus a combined
   Stats_AllOrgans sheet.

Default input directory: data/proteomics/

Default output: results/supplementary_table_4/Supplementary_Table_4_protein_abundance.xlsx
"""

from pathlib import Path
import argparse
import re
import numpy as np
import pandas as pd
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests

def parse_arguments():
    parser = argparse.ArgumentParser(description="Generate Supplementary Table 4 from deposited Spectronaut exports.")
    parser.add_argument("--input", type=Path, default=Path("data/proteomics"),
                        help="Directory containing the deposited Spectronaut CSV/TSV files.")
    parser.add_argument("--output", type=Path, default=Path("results/supplementary_table_4"),
                        help="Directory for the Supplementary Table 4 Excel workbook.")
    return parser.parse_args()

EXPECTED_FILES = {'Liver': '20250430_ArgDiet_2w_Liver_Double_Labelled_Channels_BGS_Mari_plusPGchannelsReport_GWTC',
                  'Lung': '20250430_ArgDiet_2w_Lung_Double_Labelled_Channels_BGS_Mari_plusPGchannelsReport_GWTC',
                  'Brain': '20250430_ArgDiet_2w_Brain_Double_Labelled_Channels_report',
                  'Colon': '20250430_ArgDiet_2w_Colon_Double_Labelled_Channels_report',
                  'Cortex': '20250430_ArgDiet_2w_Cortex_Double_Labelled_Channels_report',
                  'Duodenum': '20250430_ArgDiet_2w_Duodenum_Double_Labelled_Channels_report',
                  'Erythrocytes': '20250430_ArgDiet_2w_Erythrocytes_Double_Labelled_Channels_report',
                  'Heart': '20250430_ArgDiet_2w_Heart_repet_Double_Labelled_Channels_report',
                  'Ileum': '20250430_ArgDiet_2w_Ileum_Double_Labelled_Channels_report',
                  'Jejunum': '20250430_ArgDiet_2w_Jejunum_Double_Labelled_Channels_report',
                  'Medulla': '20250430_ArgDiet_2w_Medulla_Double_Labelled_Channels_report',
                  'Muscle': '20250430_ArgDiet_2w_Muscle_Double_Labelled_Channels_report',
                  'Plasma': '20250430_ArgDiet_2w_Plasma_Double_Labelled_Channels_report',
                  'Spleen': '20250430_ArgDiet_2w_Spleen_Double_Labelled_Channels_report',
                  'Urine': '20250430_ArgDiet_2w_Urine_Double_Labelled_Channels_report'}

ORGAN_ORDER = ['Cortex', 'Medulla', 'Liver', 'Heart', 'Muscle', 'Lung', 'Brain', 'Spleen', 'Duodenum', 'Jejunum',
               'Ileum', 'Colon', 'Erythrocytes', 'Plasma', 'Urine']

# Study-specific mapping from Spectronaut R.FileName to manuscript sample groups.
SAMPLE_MAP = {'Cortex': {'22_AVRQ_ArgDiet_Study_2w_Cortex_repet_No147': 'Control_ND_1',
                         '23_AVRQ_ArgDiet_Study_2w_Cortex_repet_No157': 'Proteinuric_AD_1',
                         '24_AVRQ_ArgDiet_Study_2w_Cortex_repet_No146': 'Control_ND_2',
                         '25_AVRQ_ArgDiet_Study_2w_Cortex_repet_No166': 'Proteinuric_ND_1',
                         '26_AVRQ_ArgDiet_Study_2w_Cortex_repet_No148': 'Proteinuric_AD_2',
                         '27_AVRQ_ArgDiet_Study_2w_Cortex_repet_No141': 'Control_AD_1',
                         '28_AVRQ_ArgDiet_Study_2w_Cortex_repet_No164': 'Control_ND_3',
                         '29_AVRQ_ArgDiet_Study_2w_Cortex_repet_No154': 'Control_AD_2',
                         '30_AVRQ_ArgDiet_Study_2w_Cortex_repet_No145': 'Control_ND_4',
                         '31_AVRQ_ArgDiet_Study_2w_Cortex_repet_No167': 'Proteinuric_ND_2',
                         '32_AVRQ_ArgDiet_Study_2w_Cortex_repet_No150': 'Proteinuric_AD_3',
                         '33_AVRQ_ArgDiet_Study_2w_Cortex_repet_No144': 'Control_ND_5',
                         '34_AVRQ_ArgDiet_Study_2w_Cortex_repet_No143': 'Control_AD_3',
                         '35_AVRQ_ArgDiet_Study_2w_Cortex_repet_No168': 'Proteinuric_ND_3',
                         '36_AVRQ_ArgDiet_Study_2w_Cortex_repet_No140': 'Control_AD_4',
                         '37_AVRQ_ArgDiet_Study_2w_Cortex_repet_No153': 'Proteinuric_AD_4',
                         '38_AVRQ_ArgDiet_Study_2w_Cortex_repet_No160': 'Proteinuric_AD_5',
                         '40_AVRQ_ArgDiet_Study_2w_Cortex_repet_No142': 'Proteinuric_ND_4',
                         '42_AVRQ_ArgDiet_Study_2w_Cortex_repet_No138': 'Control_AD_5'},
 'Medulla': {'22_AVRQ_ArgDiet_2w_Med_#147_20240226': 'Control_ND_1',
             '23_AVRQ_ArgDiet_2w_Med_#157_20240226': 'Proteinuric_AD_1',
             '24_AVRQ_ArgDiet_2w_Med_#146_20240226': 'Control_ND_2',
             '25_AVRQ_ArgDiet_2w_Med_#166_20240226': 'Proteinuric_ND_1',
             '26_AVRQ_ArgDiet_2w_Med_#148_20240226': 'Proteinuric_AD_2',
             '27_AVRQ_ArgDiet_2w_Med_#141_20240226': 'Control_AD_1',
             '28_AVRQ_ArgDiet_2w_Med_#164_20240226': 'Control_ND_3',
             '29_AVRQ_ArgDiet_2w_Med_#154_20240226': 'Control_AD_2',
             '30_AVRQ_ArgDiet_2w_Med_#145_20240226': 'Control_ND_4',
             '31_AVRQ_ArgDiet_2w_Med_#167_20240226': 'Proteinuric_ND_2',
             '32_AVRQ_ArgDiet_2w_Med_#150_20240226': 'Proteinuric_AD_3',
             '33_AVRQ_ArgDiet_2w_Med_#144_20240226': 'Control_ND_5',
             '34_AVRQ_ArgDiet_2w_Med_#143_20240226': 'Control_AD_3',
             '35_AVRQ_ArgDiet_2w_Med_#168_20240226': 'Proteinuric_ND_3',
             '36_AVRQ_ArgDiet_2w_Med_#140_20240226': 'Control_AD_4',
             '37_AVRQ_ArgDiet_2w_Med_#153_20240226': 'Control_ND_6',
             '38_AVRQ_ArgDiet_2w_Med_#160_20240226': 'Proteinuric_AD_4',
             '39_AVRQ_ArgDiet_2w_Med_#139_20240226': 'Proteinuric_AD_5',
             '40_AVRQ_ArgDiet_2w_Med_#142_20240226': 'Control_AD_5',
             '41_AVRQ_ArgDiet_2w_Med_#161_20240226': 'Proteinuric_ND_4',
             '42_AVRQ_ArgDiet_2w_Med_#138_20240226': 'Control_AD_6'},
 'Liver': {'43_AVRQ_ArgDiet_2w_Liver_#147_20240226': 'Proteinuric_AD_1',
           '44_AVRQ_ArgDiet_2w_Liver_#157_20240226': 'Proteinuric_AD_2',
           '45_AVRQ_ArgDiet_2w_Liver_#146_20240226': 'Proteinuric_AD_3',
           '46_AVRQ_ArgDiet_2w_Liver_#166_20240226': 'Control_AD_1',
           '47_AVRQ_ArgDiet_2w_Liver_#148_20240226': 'Control_AD_2',
           '48_AVRQ_ArgDiet_2w_Liver_#141_20240226': 'Control_AD_3',
           '49_AVRQ_ArgDiet_2w_Liver_#164_20240226': 'Control_AD_4',
           '50_AVRQ_ArgDiet_2w_Liver_#154_20240226': 'Control_AD_5',
           '51_AVRQ_ArgDiet_2w_Liver_#145_20240226': 'Control_AD_6',
           '52_AVRQ_ArgDiet_2w_Liver_#167_20240226': 'Proteinuric_AD_4',
           '53_AVRQ_ArgDiet_2w_Liver_#150_20240226': 'Proteinuric_AD_5',
           '54_AVRQ_ArgDiet_2w_Liver_#144_20240226': 'Proteinuric_ND_1',
           '55_AVRQ_ArgDiet_2w_Liver_#143_20240226': 'Proteinuric_ND_2',
           '56_AVRQ_ArgDiet_2w_Liver_#168_20240226': 'Proteinuric_ND_3',
           '57_AVRQ_ArgDiet_2w_Liver_#140_20240226': 'Proteinuric_ND_4',
           '58_AVRQ_ArgDiet_2w_Liver_#153_20240226': 'Control_ND_1',
           '59_AVRQ_ArgDiet_2w_Liver_#160_20240226': 'Control_ND_2',
           '60_AVRQ_ArgDiet_2w_Liver_#139_20240226': 'Control_ND_3',
           '61_AVRQ_ArgDiet_2w_Liver_#142_20240226': 'Control_ND_4',
           '62_AVRQ_ArgDiet_2w_Liver_#161_20240226': 'Control_ND_5',
           '63_AVRQ_ArgDiet_2w_Liver_#138_20240226': 'Control_ND_6'},
 'Heart': {'1_AVRQ_ArgDiet_Study_2w_Heart_repet_No150': 'Proteinuric_AD_1',
           '2_AVRQ_ArgDiet_Study_2w_Heart_repet_No139': 'Proteinuric_AD_2',
           '3_AVRQ_ArgDiet_Study_2w_Heart_repet_No157': 'Proteinuric_AD_3',
           '4_AVRQ_ArgDiet_Study_2w_Heart_repet_No142': 'Control_AD_1',
           '5_AVRQ_ArgDiet_Study_2w_Heart_repete_No143': 'Control_AD_2',
           '6_AVRQ_ArgDiet_Study_2w_Heart_repet_No154': 'Control_AD_3',
           '7_AVRQ_ArgDiet_Study_2w_Heart_repet_No138': 'Control_AD_4',
           '8_AVRQ_ArgDiet_Study_2w_Heart_repett_No140': 'Control_AD_5',
           '9_AVRQ_ArgDiet_Study_2w_Heart_repet_No141': 'Control_AD_6',
           '10_AVRQ_ArgDiet_Study_2w_Heart_repet_No148': 'Proteinuric_AD_4',
           '11_AVRQ_ArgDiet_Study_2w_Heart_repet_No160': 'Proteinuric_AD_5',
           '12_AVRQ_ArgDiet_Study_2w_Heart_repet_No166': 'Proteinuric_ND_1',
           '13_AVRQ_ArgDiet_Study_2w_Heart_repet_No167': 'Proteinuric_ND_2',
           '14_AVRQ_ArgDiet_Study_2w_Heart_repet_No168': 'Proteinuric_ND_3',
           '15_AVRQ_ArgDiet_Study_2w_Heart_repet_No161': 'Proteinuric_ND_4',
           '16_AVRQ_ArgDiet_Study_2w_Heart_repet_No147': 'Control_ND_1',
           '17_AVRQ_ArgDiet_Study_2w_Heart_repet_No153': 'Control_ND_2',
           '18_AVRQ_ArgDiet_Study_2w_Heart_repet_No164': 'Control_ND_3',
           '19_AVRQ_ArgDiet_Study_2w_Heart_repet_No144': 'Control_ND_4',
           '20_AVRQ_ArgDiet_Study_2w_Heart_repet_No145': 'Control_ND_5',
           '21_AVRQ_ArgDiet_Study_2w_Heart_repet_No146': 'Control_ND_6'},
 'Muscle': {'22_AVRQ_2kDiet_Study_Muscle_No150': 'Proteinuric_AD_1',
            '23_AVRQ_2kDiet_Study_Muscle_No139': 'Proteinuric_AD_2',
            '24_AVRQ_2kDiet_Study_Muscle_No157': 'Proteinuric_AD_3',
            '25_AVRQ_2kDiet_Study_Muscle_No142': 'Control_AD_1',
            '26_AVRQ_2kDiet_Study_Muscle_No143': 'Control_AD_2',
            '27_AVRQ_2kDiet_Study_Muscle_No154': 'Control_AD_3',
            '28_AVRQ_2kDiet_Study_Muscle_No138': 'Control_AD_4',
            '29_AVRQ_2kDiet_Study_Muscle_No140': 'Control_AD_5',
            '30_AVRQ_2kDiet_Study_Muscle_No141': 'Control_AD_6',
            '31_AVRQ_2kDiet_Study_Muscle_No148': 'Proteinuric_AD_4',
            '32_AVRQ_2kDiet_Study_Muscle_No160': 'Proteinuric_AD_5',
            '33_AVRQ_2kDiet_Study_Muscle_No166': 'Proteinuric_ND_1',
            '34_AVRQ_2kDiet_Study_Muscle_No167': 'Proteinuric_ND_2',
            '35_AVRQ_2kDiet_Study_Muscle_No168': 'Proteinuric_ND_3',
            '36_AVRQ_2kDiet_Study_Muscle_No161': 'Proteinuric_ND_4',
            '37_AVRQ_2kDiet_Study_Muscle_No147': 'Control_ND_1',
            '38_AVRQ_2kDiet_Study_Muscle_No153': 'Control_ND_2',
            '39_AVRQ_2kDiet_Study_Muscle_No164': 'Control_ND_3',
            '40_AVRQ_2kDiet_Study_Muscle_No144': 'Control_ND_4',
            '41_AVRQ_2kDiet_Study_Muscle_No145': 'Control_ND_5',
            '42_AVRQ_2kDiet_Study_Muscle_No146': 'Control_ND_6'},
 'Brain': {'1_AVRQ_ArgDiet_2w_Brain_#147_rep_20240604': 'Proteinuric_AD_1',
           '2_AVRQ_ArgDiet_2w_Brain_#157_rep_20240604': 'Proteinuric_AD_2',
           '3_AVRQ_ArgDiet_2w_Brain_#146_rep_20240604': 'Proteinuric_AD_3',
           '4_AVRQ_ArgDiet_2w_Brain_#166_rep_20240604': 'Control_AD_1',
           '5_AVRQ_ArgDiet_2w_Brain_#148_rep_20240604': 'Control_AD_2',
           '6_AVRQ_ArgDiet_2w_Brain_#141_rep_20240604': 'Control_AD_3',
           '7_AVRQ_ArgDiet_2w_Brain_#164_rep_20240604': 'Control_AD_4',
           '8_AVRQ_ArgDiet_2w_Brain_#154_rep_20240604': 'Control_AD_5',
           '9_AVRQ_ArgDiet_2w_Brain_#145_rep_20240604': 'Control_AD_6',
           '10_AVRQ_ArgDiet_2w_Brain_#167_rep_20240604': 'Proteinuric_AD_4',
           '11_AVRQ_ArgDiet_2w_Brain_#150_rep_20240604': 'Proteinuric_AD_5',
           '12_AVRQ_ArgDiet_2w_Brain_#144_rep_20240604': 'Proteinuric_ND_1',
           '13_AVRQ_ArgDiet_2w_Brain_#143_rep_20240604': 'Proteinuric_ND_2',
           '14_AVRQ_ArgDiet_2w_Brain_#168_rep_20240604': 'Proteinuric_ND_3',
           '15_AVRQ_ArgDiet_2w_Brain_#140_rep_20240604': 'Proteinuric_ND_4',
           '16_AVRQ_ArgDiet_2w_Brain_#153_rep_20240604': 'Control_ND_1',
           '17_AVRQ_ArgDiet_2w_Brain_#160_rep_20240604': 'Control_ND_2',
           '18_AVRQ_ArgDiet_2w_Brain_#139_rep_20240604': 'Control_ND_3',
           '19_AVRQ_ArgDiet_2w_Brain_#142_rep_20240604': 'Control_ND_4',
           '20_AVRQ_ArgDiet_2w_Brain_#161_rep_20240604': 'Control_ND_5',
           '21_AVRQ_ArgDiet_2w_Brain_#138_rep_20240604': 'Control_ND_6'},
 'Lung': {'64_AVRQ_ArgDiet_2w_Lung_#147_20240226': 'Proteinuric_AD_1',
          '65_AVRQ_ArgDiet_2w_Lung_#157_20240226': 'Proteinuric_AD_2',
          '66_AVRQ_ArgDiet_2w_Lung_#146_20240226': 'Proteinuric_AD_3',
          '67_AVRQ_ArgDiet_2w_Lung_#166_20240226': 'Control_AD_1',
          '68_AVRQ_ArgDiet_2w_Lung_#148_20240226': 'Control_AD_2',
          '69_AVRQ_ArgDiet_2w_Lung_#141_20240226': 'Control_AD_3',
          '70_AVRQ_ArgDiet_2w_Lung_#164_20240226': 'Control_AD_4',
          '71_AVRQ_ArgDiet_2w_Lung_#154_20240226': 'Control_AD_5',
          '72_AVRQ_ArgDiet_2w_Lung_#145_20240226': 'Control_AD_6',
          '73_AVRQ_ArgDiet_2w_Lung_#167_20240226': 'Proteinuric_AD_4',
          '74_AVRQ_ArgDiet_2w_Lung_#150_20240226_20240306151522': 'Proteinuric_AD_5',
          '75_AVRQ_ArgDiet_2w_Lung_#144_20240226': 'Proteinuric_ND_1',
          '76_AVRQ_ArgDiet_2w_Lung_#143_20240226': 'Proteinuric_ND_2',
          '77_AVRQ_ArgDiet_2w_Lung_#168_20240226': 'Proteinuric_ND_3',
          '78_AVRQ_ArgDiet_2w_Lung_#140_20240226': 'Proteinuric_ND_4',
          '79_AVRQ_ArgDiet_2w_Lung_#153_20240226': 'Control_ND_1',
          '80_AVRQ_ArgDiet_2w_Lung_#160_20240226': 'Control_ND_2',
          '81_AVRQ_ArgDiet_2w_Lung_#139_20240226': 'Control_ND_3',
          '82_AVRQ_ArgDiet_2w_Lung_#142_20240226': 'Control_ND_4',
          '83_AVRQ_ArgDiet_2w_Lung_#161_20240226': 'Control_ND_5',
          '84_AVRQ_ArgDiet_2w_Lung_#138_20240226': 'Control_ND_6'},
 'Ileum': {'AVRQ_2kDiet_Study_Ileum_No147': 'Control_ND_1',
           'AVRQ_2kDiet_Study_Ileum_No153': 'Control_ND_2',
           'AVRQ_2kDiet_Study_Ileum_No164': 'Control_ND_3',
           'AVRQ_2kDiet_Study_Ileum_No144': 'Control_ND_4',
           'AVRQ_2kDiet_Study_Ileum_No145': 'Control_ND_5',
           'AVRQ_2kDiet_Study_Ileum_No146': 'Control_ND_6',
           'AVRQ_2kDiet_Study_Ileum_No166': 'Proteinuric_ND_1',
           'AVRQ_2kDiet_Study_Ileum_No167': 'Proteinuric_ND_2',
           'AVRQ_2kDiet_Study_Ileum_No168': 'Proteinuric_ND_3',
           'AVRQ_2kDiet_Study_Ileum_No161': 'Proteinuric_ND_4',
           'AVRQ_2kDiet_Study_Ileum_No142': 'Control_AD_1',
           'AVRQ_2kDiet_Study_Ileum_No143': 'Control_AD_2',
           'AVRQ_2kDiet_Study_Ileum_No154': 'Control_AD_3',
           'AVRQ_2kDiet_Study_Ileum_No138': 'Control_AD_4',
           'AVRQ_2kDiet_Study_Ileum_No140': 'Control_AD_5',
           'AVRQ_2kDiet_Study_Ileum_No141': 'Control_AD_6',
           'AVRQ_2kDiet_Study_Ileum_No150': 'Proteinuric_AD_1',
           'AVRQ_2kDiet_Study_Ileum_No139': 'Proteinuric_AD_2',
           'AVRQ_2kDiet_Study_Ileum_No157': 'Proteinuric_AD_3',
           'AVRQ_2kDiet_Study_Ileum_No148': 'Proteinuric_AD_4',
           'AVRQ_2kDiet_Study_Ileum_No160': 'Proteinuric_AD_5'},
 'Duodenum': {'1_AVRQ_ArgDiet_2w_Duodenum_#150': 'Proteinuric_AD_1',
              '2_AVRQ_ArgDiet_2w_Duodenum_#139': 'Proteinuric_AD_2',
              '3_AVRQ_ArgDiet_2w_Duodenum_#157': 'Proteinuric_AD_3',
              '4_AVRQ_ArgDiet_2w_Duodenum_#142': 'Control_AD_1',
              '5_AVRQ_ArgDiet_2w_Duodenum_#143': 'Control_AD_2',
              '6_AVRQ_ArgDiet_2w_Duodenum_#154': 'Control_AD_3',
              '7_AVRQ_ArgDiet_2w_Duodenum_#138': 'Control_AD_4',
              '8_AVRQ_ArgDiet_2w_Duodenum_#140': 'Control_AD_5',
              '9_AVRQ_ArgDiet_2w_Duodenum_#141': 'Control_AD_6',
              '10_AVRQ_ArgDiet_2w_Duodenum_#148': 'Proteinuric_AD_4',
              '11_AVRQ_ArgDiet_2w_Duodenum_#160': 'Proteinuric_AD_5',
              '12_AVRQ_ArgDiet_2w_Duodenum_#166': 'Proteinuric_ND_1',
              '13_AVRQ_ArgDiet_2w_Duodenum_#167': 'Proteinuric_ND_2',
              '14_AVRQ_ArgDiet_2w_Duodenum_#168': 'Proteinuric_ND_3',
              '15_AVRQ_ArgDiet_2w_Duodenum_#161': 'Proteinuric_ND_4',
              '16_AVRQ_ArgDiet_2w_Duodenum_#147': 'Control_ND_1',
              '17_AVRQ_ArgDiet_2w_Duodenum_#153': 'Control_ND_2',
              '18_AVRQ_ArgDiet_2w_Duodenum_#164': 'Control_ND_3',
              '19_AVRQ_ArgDiet_2w_Duodenum_#144': 'Control_ND_4',
              '20_AVRQ_ArgDiet_2w_Duodenum_#145': 'Control_ND_5',
              '21_AVRQ_ArgDiet_2w_Duodenum_#146': 'Control_ND_6'},
 'Plasma': {'1_AVRQ_ArgDiet_2w_Plasma_#150_20240510': 'Proteinuric_AD_1',
            '2_AVRQ_ArgDiet_2w_Plasma_#139_20240510': 'Proteinuric_AD_2',
            '3_AVRQ_ArgDiet_2w_Plasma_#157_20240510': 'Proteinuric_AD_3',
            '4_AVRQ_ArgDiet_2w_Plasma_#142_20240510': 'Control_AD_1',
            '5_AVRQ_ArgDiet_2w_Plasma_#143_20240510': 'Control_AD_2',
            '6_AVRQ_ArgDiet_2w_Plasma_#154_20240510': 'Control_AD_3',
            '7_AVRQ_ArgDiet_2w_Plasma_#138_20240510': 'Control_AD_4',
            '8_AVRQ_ArgDiet_2w_Plasma_#140_20240510': 'Control_AD_5',
            '9_AVRQ_ArgDiet_2w_Plasma_#141_20240510': 'Control_AD_6',
            '10_AVRQ_ArgDiet_2w_Plasma_#148_20240510': 'Proteinuric_AD_4',
            '11_AVRQ_ArgDiet_2w_Plasma_#160_20240510': 'Proteinuric_AD_5',
            '12_AVRQ_ArgDiet_2w_Plasma_#166_20240510': 'Proteinuric_ND_1',
            '13_AVRQ_ArgDiet_2w_Plasma_#167_20240510': 'Proteinuric_ND_2',
            '14_AVRQ_ArgDiet_2w_Plasma_#168_20240510': 'Proteinuric_ND_3',
            '15_AVRQ_ArgDiet_2w_Plasma_#161_20240510': 'Proteinuric_ND_4',
            '16_AVRQ_ArgDiet_2w_Plasma_#147_20240510': 'Control_ND_1',
            '17_AVRQ_ArgDiet_2w_Plasma_#153_20240510': 'Control_ND_2',
            '18_AVRQ_ArgDiet_2w_Plasma_#164_20240510': 'Control_ND_3',
            '19_AVRQ_ArgDiet_2w_Plasma_#144_20240510': 'Control_ND_4',
            '20_AVRQ_ArgDiet_2w_Plasma_#145_20240510': 'Control_ND_5',
            '21_AVRQ_ArgDiet_2w_Plasma_#146_20240510': 'Control_ND_6'},
 'Urine': {'1_AVRQ_ArgDiet_2w_Urine_#150_20240513': 'Proteinuric_AD_1',
           '2_AVRQ_ArgDiet_2w_Urine_#168_20240513': 'Proteinuric_AD_2',
           '3_AVRQ_ArgDiet_2w_Urine_#139_20240513': 'Proteinuric_AD_3',
           '4_AVRQ_ArgDiet_2w_Urine_#153_20240513': 'Proteinuric_AD_4',
           '5_AVRQ_ArgDiet_2w_Urine_#142_20240513': 'Control_AD_1',
           '6_AVRQ_ArgDiet_2w_Urine_#164_20240513': 'Control_AD_2',
           '7_AVRQ_ArgDiet_2w_Urine_#143_20240513': 'Control_AD_3',
           '8_AVRQ_ArgDiet_2w_Urine_#144_20240513': 'Control_AD_4',
           '9_AVRQ_ArgDiet_2w_Urine_#154_20240513': 'Control_AD_5',
           '10_AVRQ_ArgDiet_2w_Urine_#146_20240513': 'Control_AD_6',
           '11_AVRQ_ArgDiet_2w_Urine_#138_20240513': 'Control_ND_1',
           '12_AVRQ_ArgDiet_2w_Urine_#140_20240513': 'Control_ND_2',
           '13_AVRQ_ArgDiet_2w_Urine_#141_20240513': 'Control_ND_3',
           '14_AVRQ_ArgDiet_2w_Urine_#148_20240513': 'Control_ND_4',
           '15_AVRQ_ArgDiet_2w_Urine_#160_20240513': 'Proteinuric_ND_1',
           '16_AVRQ_ArgDiet_2w_Urine_#166_20240513': 'Proteinuric_ND_2',
           '17_AVRQ_ArgDiet_2w_Urine_#167_20240513': 'Proteinuric_ND_3'},
 'Erythrocytes': {'1_AVRQ_ArgDiet_Study_2w_Erythrocytes_No150': 'Proteinuric_AD_1',
                  '2_AVRQ_ArgDiet_Study_2w_Erythrocytes_No139': 'Proteinuric_AD_2',
                  '3_AVRQ_ArgDiet_Study_2w_Erythrocytes_No157': 'Proteinuric_AD_3',
                  '4_AVRQ_ArgDiet_Study_2w_Erythrocytes_No142': 'Control_AD_1',
                  '5_AVRQ_ArgDiet_Study_2w_Erythrocytes_No143': 'Control_AD_2',
                  '6_AVRQ_ArgDiet_Study_2w_Erythrocytes_No154': 'Control_AD_3',
                  '7_AVRQ_ArgDiet_Study_2w_Erythrocytes_No138': 'Control_AD_4',
                  '8_AVRQ_ArgDiet_Study_2w_Erythrocytes_No140': 'Control_AD_5',
                  '9_AVRQ_ArgDiet_Study_2w_Erythrocytes_No141': 'Control_AD_6',
                  '10_AVRQ_ArgDiet_Study_2w_Erythrocytes_No148': 'Proteinuric_AD_4',
                  '11_AVRQ_ArgDiet_Study_2w_Erythrocytes_No160': 'Proteinuric_AD_5',
                  '12_AVRQ_ArgDiet_Study_2w_Erythrocytes_No166': 'Proteinuric_ND_1',
                  '13_AVRQ_ArgDiet_Study_2w_Erythrocytes_No167': 'Proteinuric_ND_2',
                  '14_AVRQ_ArgDiet_Study_2w_Erythrocytes_No168': 'Proteinuric_ND_3',
                  '15_AVRQ_ArgDiet_Study_2w_Erythrocytes_No161': 'Proteinuric_ND_4',
                  '16_AVRQ_ArgDiet_Study_2w_Erythrocytes_No147': 'Control_ND_1',
                  '17_AVRQ_ArgDiet_Study_2w_Erythrocytes_No153': 'Control_ND_2',
                  '18_AVRQ_ArgDiet_Study_2w_Erythrocytes_No164': 'Control_ND_3',
                  '19_AVRQ_ArgDiet_Study_2w_Erythrocytes_No144': 'Control_ND_4',
                  '20_AVRQ_ArgDiet_Study_2w_Erythrocytes_No145': 'Control_ND_5',
                  '21_AVRQ_ArgDiet_Study_2w_Erythrocytes_No146': 'Control_ND_6'},
 'Colon': {'1_AVRQ_ArgDiet_Study_2w_Colon_repet_No150': 'Proteinuric_AD_1',
           '2_AVRQ_ArgDiet_Study_2w_Colon_repet_No139': 'Proteinuric_AD_2',
           '3_AVRQ_ArgDiet_Study_2w_Colon_repet_No157': 'Proteinuric_AD_3',
           '4_AVRQ_ArgDiet_Study_2w_Colon_repet_No142': 'Control_AD_1',
           '5_AVRQ_ArgDiet_Study_2w_Colon_repet_No143': 'Control_AD_2',
           '6_AVRQ_ArgDiet_Study_2w_Colon_repet_No154': 'Control_AD_3',
           '7_AVRQ_ArgDiet_Study_2w_Colon_repet_No138': 'Control_AD_4',
           '8_AVRQ_ArgDiet_Study_2w_Colon_repet_No140': 'Control_AD_5',
           '9_AVRQ_ArgDiet_Study_2w_Colon_repet_No141': 'Control_AD_6',
           '10_AVRQ_ArgDiet_Study_2w_Colon_repet_No148': 'Proteinuric_AD_4',
           '11_AVRQ_ArgDiet_Study_2w_Colon_repet_No160': 'Proteinuric_AD_5',
           '12_AVRQ_ArgDiet_Study_2w_Colon_repet_No166': 'Proteinuric_ND_1',
           '13_AVRQ_ArgDiet_Study_2w_Colon_repet_No167': 'Proteinuric_ND_2',
           '14_AVRQ_ArgDiet_Study_2w_Colon_repet_No168': 'Proteinuric_ND_3',
           '15_AVRQ_ArgDiet_Study_2w_Colon_repet_No161': 'Proteinuric_ND_4',
           '16_AVRQ_ArgDiet_Study_2w_Colon_repet_No147': 'Control_ND_1',
           '17_AVRQ_ArgDiet_Study_2w_Colon_repet_No153': 'Control_ND_2',
           '18_AVRQ_ArgDiet_Study_2w_Colon_repet_No164': 'Control_ND_3',
           '19_AVRQ_ArgDiet_Study_2w_Colon_repet_No144': 'Control_ND_4',
           '20_AVRQ_ArgDiet_Study_2w_Colon_repet_No145': 'Control_ND_5',
           '21_AVRQ_ArgDiet_Study_2w_Colon_repet_No146': 'Control_ND_6'},
 'Spleen': {'22_AVRQ_ArgDiet_Study_2w_Spleen_repet_No150': 'Proteinuric_AD_1',
            '23_AVRQ_ArgDiet_Study_2w_Spleen_repet_No139': 'Proteinuric_AD_2',
            '24_AVRQ_ArgDiet_Study_2w_Spleen_repet_No157': 'Proteinuric_AD_3',
            '25_AVRQ_ArgDiet_Study_2w_Spleen_repet_No142': 'Control_AD_1',
            '26_AVRQ_ArgDiet_Study_2w_Spleen_repet_No143': 'Control_AD_2',
            '27_AVRQ_ArgDiet_Study_2w_Spleen_repet_No154': 'Control_AD_3',
            '28_AVRQ_ArgDiet_Study_2w_Spleen_repet_No138': 'Control_AD_4',
            '29_AVRQ_ArgDiet_Study_2w_Spleen_repet_No140': 'Control_AD_5',
            '30_AVRQ_ArgDiet_Study_2w_Spleen_repet_No141': 'Control_AD_6',
            '31_AVRQ_ArgDiet_Study_2w_Spleen_repet_No148': 'Proteinuric_AD_4',
            '32_AVRQ_ArgDiet_Study_2w_Spleen_repet_No160': 'Proteinuric_AD_5',
            '33_AVRQ_ArgDiet_Study_2w_Spleen_repet_No166': 'Proteinuric_ND_1',
            '34_AVRQ_ArgDiet_Study_2w_Spleen_repet_No167': 'Proteinuric_ND_2',
            '35_AVRQ_ArgDiet_Study_2w_Spleen_repet_No168': 'Proteinuric_ND_3',
            '36_AVRQ_ArgDiet_Study_2w_Spleen_repet_No161': 'Proteinuric_ND_4',
            '37_AVRQ_ArgDiet_Study_2w_Spleen_repet_No147': 'Control_ND_1',
            '38_AVRQ_ArgDiet_Study_2w_Spleen_repet_No153': 'Control_ND_2',
            '39_AVRQ_ArgDiet_Study_2w_Spleen_repet_No164': 'Control_ND_3',
            '40_AVRQ_ArgDiet_Study_2w_Spleen_repet_No144': 'Control_ND_4',
            '41_AVRQ_ArgDiet_Study_2w_Spleen_repet_No145': 'Control_ND_5',
            '42_AVRQ_ArgDiet_Study_2w_Spleen_repet_No146': 'Control_ND_6'},
 'Jejunum': {'1_AVRQ_ArgDiet_2w_Jejenum_#150': 'Proteinuric_AD_1',
             '2_AVRQ_ArgDiet_2w_Jejenum_#139': 'Proteinuric_AD_2',
             '3_AVRQ_ArgDiet_2w_Jejenum_#157': 'Proteinuric_AD_3',
             '4_AVRQ_ArgDiet_2w_Jejenum_#142': 'Control_AD_1',
             '5_AVRQ_ArgDiet_2w_Jejenum_#143': 'Control_AD_2',
             '6_AVRQ_ArgDiet_2w_Jejenum_#154': 'Control_AD_3',
             '7_AVRQ_ArgDiet_2w_Jejenum_#138': 'Control_AD_4',
             '8_AVRQ_ArgDiet_2w_Jejenum_#140': 'Control_AD_5',
             '9_AVRQ_ArgDiet_2w_Jejenum_#141': 'Control_AD_6',
             '10_AVRQ_ArgDiet_2w_Jejenum_#148': 'Proteinuric_AD_4',
             '11_AVRQ_ArgDiet_2w_Jejenum_#160': 'Proteinuric_AD_5',
             '12_AVRQ_ArgDiet_2w_Jejenum_#166': 'Proteinuric_ND_1',
             '13_AVRQ_ArgDiet_2w_Jejenum_#167': 'Proteinuric_ND_2',
             '14_AVRQ_ArgDiet_2w_Jejenum_#168': 'Proteinuric_ND_3',
             '15_AVRQ_ArgDiet_2w_Jejenum_#161': 'Proteinuric_ND_4',
             '16_AVRQ_ArgDiet_2w_Jejenum_#147': 'Control_ND_1',
             '17_AVRQ_ArgDiet_2w_Jejenum_#153': 'Control_ND_2',
             '18_AVRQ_ArgDiet_2w_Jejenum_#164': 'Control_ND_3',
             '19_AVRQ_ArgDiet_2w_Jejenum_#144': 'Control_ND_4',
             '20_AVRQ_ArgDiet_2w_Jejenum_#145': 'Control_ND_5',
             '21_AVRQ_ArgDiet_2w_Jejenum_#146': 'Control_ND_6'}}

FAIL_ON_UNMAPPED = True
BASE_COLUMNS = ["PG.Organisms", "PG.FastaFiles", "R.FileName", "R.Condition", "PG.Genes", "PG.ProteinNames",
                "PG.UniProtIds", "PG.NrOfStrippedSequencesUsedForQuantification", "PG.Coverage", "PG.Pvalue",
                "PG.Qvalue", "PG.Log2Quantity", "PG.IBAQ", "PG.MS1Quantity", "PG.MS2Quantity", "PG.Quantity",
                "PG.MS2Channel1", "PG.MS2Channel2", "PG.MS2Channel3"]
METRICS = ["PG.Log2Quantity", "PG.IBAQ", "PG.MS1Quantity", "PG.MS2Quantity", "PG.Quantity",
           "PG.MS2Channel1", "PG.MS2Channel2", "PG.MS2Channel3"]
COMPARISON_PAIRS = [("ND", "Control_ND", "Proteinuric_ND"), ("AD", "Control_AD", "Proteinuric_AD")]
STAT_CHANNELS = ["PG.MS2Channel2", "PG.MS2Channel3"]


def find_input_file(input_dir: Path, file_stem: str) -> Path:
    """Locate one deposited file by exact stem, allowing CSV or TSV."""
    candidates = [input_dir / f"{file_stem}.csv", input_dir / f"{file_stem}.tsv"]
    existing = [path for path in candidates if path.exists()]

    if len(existing) == 0:
        raise FileNotFoundError(f"Could not find '{file_stem}.csv' or '{file_stem}.tsv' in {input_dir.resolve()}")
    if len(existing) > 1:
        raise RuntimeError(f"Both CSV and TSV versions were found for '{file_stem}'. "
                           "Keep only the deposited file used for this analysis.")
    return existing[0]


def read_spectronaut_file(file_path: Path) -> pd.DataFrame:
    """Read one Spectronaut protein-group export."""
    if file_path.suffix.lower() == ".tsv":
        return pd.read_csv(file_path, sep="\t", low_memory=False)
    if file_path.suffix.lower() == ".csv":
        return pd.read_csv(file_path, low_memory=False)
    raise ValueError(f"Unsupported file type: {file_path}")


def load_one_export(file_path: Path, organ: str) -> pd.DataFrame:
    df = read_spectronaut_file(file_path)
    keep = [column for column in BASE_COLUMNS if column in df.columns]
    df = df[keep].copy()
    df["Organ"] = organ

    mapping = SAMPLE_MAP.get(organ, {})
    if mapping:
        df["R.FileName"] = df["R.FileName"].map(mapping).fillna(df["R.FileName"])

    return df

def check_unmapped_samples(df_all: pd.DataFrame) -> None:
    problems = []

    for organ, df_organ in df_all.groupby("Organ", sort=False):
        mapping = SAMPLE_MAP.get(organ, {})
        if not mapping:
            continue

        mapped_names = set(mapping.values())
        unmapped = (df_organ.loc[~df_organ["R.FileName"].isin(mapped_names), ["Organ", "R.FileName"]].drop_duplicates()
                    .sort_values("R.FileName"))

        if not unmapped.empty:
            problems.append(unmapped)

    if problems:
        problems_df = pd.concat(problems, ignore_index=True)
        print("\nUnmapped samples detected:")
        print(problems_df.to_string(index=False))

        if FAIL_ON_UNMAPPED:
            raise SystemExit("\nStopping because FAIL_ON_UNMAPPED=True. "
                             "The deposited R.FileName values must match SAMPLE_MAP.")


def make_wide(df: pd.DataFrame) -> pd.DataFrame:
    present_metrics = [metric for metric in METRICS if metric in df.columns]
    if not present_metrics:
        raise ValueError("None of the requested metrics exist in the dataframe.")

    id_cols = ["Organ", "PG.Genes", "PG.ProteinNames", "PG.UniProtIds"]
    id_cols = [column for column in id_cols if column in df.columns]

    long = df.melt(id_vars=id_cols + ["R.FileName"], value_vars=present_metrics, var_name="Metric", value_name="Value")

    long["Sample_Metric"] = (long["R.FileName"].astype(str) + "__" + long["Metric"].astype(str))

    wide = (long.pivot_table(index=id_cols, columns="Sample_Metric", values="Value", aggfunc="first").reset_index())
    wide.columns.name = None

    sample_metric_cols = [column for column in wide.columns if column not in id_cols and "__" in column]
    samples = sorted({column.split("__", 1)[0] for column in sample_metric_cols})

    ordered = []
    for metric in present_metrics:
        for sample in samples:
            column = f"{sample}__{metric}"
            if column in wide.columns:
                ordered.append(column)

    return wide[id_cols + ordered]

def infer_sample_group(short_name: str) -> str:
    if pd.isna(short_name):
        return np.nan

    match = re.match(r"^(Control|Proteinuric)_(ND|AD)_\d+$", str(short_name))
    if match:
        return f"{match.group(1)}_{match.group(2)}"
    return np.nan

def make_stats_table(df_organ: pd.DataFrame) -> pd.DataFrame:
    """
    Compare Control vs Proteinuric within ND and AD separately for Channel 2
    and Channel 3. Raw p-values are calculated here; BH-FDR is applied later
    across all tissues and both labeled channels within each state.
    """
    rows = []

    needed = {"R.FileName", "PG.MS2Channel2", "PG.MS2Channel3"}
    if not needed.issubset(df_organ.columns):
        return pd.DataFrame()

    per_mouse = (df_organ.groupby("R.FileName")[STAT_CHANNELS].sum().reset_index())
    per_mouse["SampleGroup"] = per_mouse["R.FileName"].apply(infer_sample_group)

    for state_label, control_label, proteinuric_label in COMPARISON_PAIRS:
        subset = per_mouse[per_mouse["SampleGroup"].isin([control_label, proteinuric_label])].copy()

        if subset.empty:
            continue

        for channel in STAT_CHANNELS:
            control = subset.loc[subset["SampleGroup"] == control_label, channel].astype(float).dropna()
            proteinuric = subset.loc[subset["SampleGroup"] == proteinuric_label, channel].astype(float).dropna()
            n_control = len(control)
            n_proteinuric = len(proteinuric)

            if n_control == 0 or n_proteinuric == 0:
                continue

            mean_control = control.mean()
            mean_proteinuric = proteinuric.mean()
            sd_control = control.std(ddof=1) if n_control > 1 else np.nan
            sd_proteinuric = proteinuric.std(ddof=1) if n_proteinuric > 1 else np.nan
            sem_control = sd_control / np.sqrt(n_control) if n_control > 1 else np.nan
            sem_proteinuric = sd_proteinuric / np.sqrt(n_proteinuric) if n_proteinuric > 1 else np.nan

            log2fc = (np.log2(mean_proteinuric / mean_control) if mean_control > 0 and mean_proteinuric > 0 else np.nan)

            if n_control > 1 and n_proteinuric > 1:
                t_statistic, p_value = ttest_ind(control, proteinuric, equal_var=False, nan_policy="omit")
            else:
                t_statistic, p_value = np.nan, np.nan

            rows.append({"State": state_label, "Channel": channel, "ControlGroup": control_label,
                         "ProteinuricGroup": proteinuric_label, "N_control": n_control, "N_proteinuric": n_proteinuric,
                         "Mean_control": mean_control, "Mean_proteinuric": mean_proteinuric, "SD_control": sd_control,
                         "SD_proteinuric": sd_proteinuric, "SEM_control": sem_control,
                         "SEM_proteinuric": sem_proteinuric, "log2FC_proteinuric_vs_control": log2fc,
                         "t_statistic": t_statistic, "p_value": p_value})
    return pd.DataFrame(rows)


def main():
    args = parse_arguments()
    input_dir = args.input
    output_dir = args.output

    if not input_dir.exists():
        raise FileNotFoundError(f"Input directory does not exist: {input_dir.resolve()}")

    output_dir.mkdir(parents=True, exist_ok=True)
    output_excel = output_dir / "Supplementary_Table_4_protein_abundance.xlsx"
    dataframes = []
    files_used = []

    for organ in ORGAN_ORDER:
        file_path = find_input_file(input_dir, EXPECTED_FILES[organ])
        df = load_one_export(file_path, organ)
        dataframes.append(df)
        files_used.append((organ, file_path.name))
    df_all = pd.concat(dataframes, ignore_index=True)

    print("Input files used:")
    for organ, filename in files_used:
        print(f"  {organ}: {filename}")

    check_unmapped_samples(df_all)
    wide_tables = {}
    all_stats = []

    for organ in ORGAN_ORDER:
        df_organ = df_all.loc[df_all["Organ"] == organ].copy()
        if df_organ.empty:
            continue

        wide_tables[organ] = make_wide(df_organ)

        stats_organ = make_stats_table(df_organ)
        if not stats_organ.empty:
            stats_organ = stats_organ.copy()
            stats_organ.insert(0, "Organ", organ)
            all_stats.append(stats_organ)

    if all_stats:
        stats_all_df = pd.concat(all_stats, ignore_index=True)
        stats_all_df["p_adj_BH"] = np.nan

        for state in ["ND", "AD"]:
            mask = ((stats_all_df["State"] == state) & stats_all_df["p_value"].notna())

            if mask.any():
                stats_all_df.loc[mask, "p_adj_BH"] = multipletests(stats_all_df.loc[mask, "p_value"], method="fdr_bh")[1]

        stats_all_df["FDR_significant_0.05"] = (stats_all_df["p_adj_BH"] <= 0.05)
    else:
        stats_all_df = pd.DataFrame()

    with pd.ExcelWriter(output_excel, engine="xlsxwriter") as writer:
        for organ in ORGAN_ORDER:
            if organ not in wide_tables:
                continue

            wide_organ = wide_tables[organ]
            sheet_name = str(organ)[:31]

            wide_organ.to_excel(writer, sheet_name=sheet_name, index=False, startrow=0)
            worksheet = writer.sheets[sheet_name]

            if not stats_all_df.empty:
                stats_organ = (stats_all_df.loc[stats_all_df["Organ"] == organ].drop(columns=["Organ"])
                               .reset_index(drop=True))
            else:
                stats_organ = pd.DataFrame()

            startrow_stats = len(wide_organ) + 3
            worksheet.write(startrow_stats - 1, 0, "Control vs Proteinuric within ND and AD; Welch's t-test; "
                                                   "BH-FDR corrected across all tissues and labeled channels "
                                                   "separately within each state")

            stats_organ.to_excel(writer, sheet_name=sheet_name, index=False, startrow=startrow_stats)

        if not stats_all_df.empty:
            stats_all_df.to_excel(writer, sheet_name="Stats_AllOrgans", index=False)

    print(f"\nSaved: {output_excel.resolve()}")

    if not stats_all_df.empty:
        for state in ["ND", "AD"]:
            n_tests = ((stats_all_df["State"] == state) & stats_all_df["p_value"].notna()).sum()
            print(f"{state}: BH-FDR correction applied across "
                  f"{n_tests} valid tissue x channel tests.")

if __name__ == "__main__":
    main()
