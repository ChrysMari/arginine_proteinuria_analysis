"""
System-wide protein isotope-channel analysis.

This script processes the 15 deposited Spectronaut protein-group exports used
for the system-wide protein isotope-channel analysis.

Channel assignment is performed in Spectronaut:
    PG.MS2Channel1 = Unlabeled
    PG.MS2Channel2 = Arg10, Pro6
    PG.MS2Channel3 = Arg7, Pro6

For each tissue and experimental condition, the script:
1. Sums protein-group channel intensities within each biological replicate.
2. Calculates mean replicate-level channel intensities within each condition.
3. Calculates the percentage contribution of each isotope channel.
4. Reproduces the statistical analysis used for Supplementary Table:
   - sample identities are mapped to Control/Proteinuric and ND/AD groups;
   - Channel 2 and Channel 3 are summed per mouse;
   - Control vs Proteinuric are compared separately within ND and AD using two-sided Welch's t-tests;
   - Benjamini-Hochberg FDR correction is applied across all tissue x labeled-
     channel tests together, separately within ND and AD.
5. Generates stacked-bar plots. FDR significance stars are shown only on the absolute-abundance plots and use the same
statistics as Supplementary Table.

Significance display:
    Orange * : Channel 2 (Arg10, Pro6), BH-FDR <= 0.05
    Blue   * : Channel 3 (Arg7, Pro6), BH-FDR <= 0.05

Input: Place the deposited Spectronaut files in: data/proteomics/

Files may be .csv or .tsv - Exact deposited file stems are listed in EXPECTED_FILES below.

Output: Derived tables, statistics, and figures are written to: results/channel_intensities/"""

from pathlib import Path
import argparse
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests

def parse_arguments():
    parser = argparse.ArgumentParser(description="Analyze Spectronaut protein isotope-channel intensities.")
    parser.add_argument("--input", type=Path, default=Path("data/proteomics"),
                        help="Directory containing the deposited Spectronaut CSV/TSV files.")
    parser.add_argument("--output", type=Path, default=Path("results/channel_intensities"),
                        help="Directory for derived tables, statistics, and figures.")
    return parser.parse_args()

# Exact deposited file stems. Files may have either .csv or .tsv extension.
EXPECTED_FILES = {"Liver": "20250430_ArgDiet_2w_Liver_Double_Labelled_Channels_BGS_Mari_plusPGchannelsReport_GWTC",
                  "Lung": "20250430_ArgDiet_2w_Lung_Double_Labelled_Channels_BGS_Mari_plusPGchannelsReport_GWTC",
                  "Brain": "20250430_ArgDiet_2w_Brain_Double_Labelled_Channels_report",
                  "Colon": "20250430_ArgDiet_2w_Colon_Double_Labelled_Channels_report",
                  "Cortex": "20250430_ArgDiet_2w_Cortex_Double_Labelled_Channels_report",
                  "Duodenum": "20250430_ArgDiet_2w_Duodenum_Double_Labelled_Channels_report",
                  "Erythrocytes": "20250430_ArgDiet_2w_Erythrocytes_Double_Labelled_Channels_report",
                  "Heart": "20250430_ArgDiet_2w_Heart_repet_Double_Labelled_Channels_report",
                  "Ileum": "20250430_ArgDiet_2w_Ileum_Double_Labelled_Channels_report",
                  "Jejunum": "20250430_ArgDiet_2w_Jejunum_Double_Labelled_Channels_report",
                  "Medulla": "20250430_ArgDiet_2w_Medulla_Double_Labelled_Channels_report",
                  "Muscle": "20250430_ArgDiet_2w_Muscle_Double_Labelled_Channels_report",
                  "Plasma": "20250430_ArgDiet_2w_Plasma_Double_Labelled_Channels_report",
                  "Spleen": "20250430_ArgDiet_2w_Spleen_Double_Labelled_Channels_report",
                  "Urine": "20250430_ArgDiet_2w_Urine_Double_Labelled_Channels_report"}

CHANNEL_COLUMNS = ["PG.MS2Channel1", "PG.MS2Channel2", "PG.MS2Channel3"]
STAT_CHANNELS = ["PG.MS2Channel2", "PG.MS2Channel3"]
CHANNEL_LABELS = {"PG.MS2Channel1": "Unlabeled", "PG.MS2Channel2": "Arg10, Pro6", "PG.MS2Channel3": "Arg7, Pro6"}
EXPERIMENTAL_GROUPS = {"Normal": ["Normal_Control", "Normal_Diseased"], "Heavy": ["Heavy_Control", "Heavy_Diseased"]}

# The Spectronaut plot-group names correspond to the ND/AD labels used in the Supplementary Table.
PLOT_GROUP_TO_STATE = {"Normal": "ND", "Heavy": "AD"}
COMPARISON_PAIRS = [("ND", "Control_ND", "Proteinuric_ND"), ("AD", "Control_AD", "Proteinuric_AD")]
ORGAN_ORDER = ["Cortex", "Medulla", "Liver", "Heart", "Muscle", "Lung", "Brain", "Spleen", "Duodenum", "Jejunum",
               "Ileum", "Colon", "Erythrocytes", "Plasma", "Urine"]

COLORS = [
    "#d3d3d3",  # Unlabeled
    "#FFB450",   # Arg10, Pro6
    "#BCCADB",   # Arg7, Pro6
]

CHANNEL_STAR_COLORS = {"PG.MS2Channel2": "#FFB450", "PG.MS2Channel3": "#BCCADB"}

SAMPLE_MAP = {
    "Cortex": {"22_AVRQ_ArgDiet_Study_2w_Cortex_repet_No147": "Control_ND_1",
               "23_AVRQ_ArgDiet_Study_2w_Cortex_repet_No157": "Proteinuric_AD_1",
               "24_AVRQ_ArgDiet_Study_2w_Cortex_repet_No146": "Control_ND_2",
               "25_AVRQ_ArgDiet_Study_2w_Cortex_repet_No166": "Proteinuric_ND_1",
               "26_AVRQ_ArgDiet_Study_2w_Cortex_repet_No148": "Proteinuric_AD_2",
               "27_AVRQ_ArgDiet_Study_2w_Cortex_repet_No141": "Control_AD_1",
               "28_AVRQ_ArgDiet_Study_2w_Cortex_repet_No164": "Control_ND_3",
               "29_AVRQ_ArgDiet_Study_2w_Cortex_repet_No154": "Control_AD_2",
               "30_AVRQ_ArgDiet_Study_2w_Cortex_repet_No145": "Control_ND_4",
               "31_AVRQ_ArgDiet_Study_2w_Cortex_repet_No167": "Proteinuric_ND_2",
               "32_AVRQ_ArgDiet_Study_2w_Cortex_repet_No150": "Proteinuric_AD_3",
               "33_AVRQ_ArgDiet_Study_2w_Cortex_repet_No144": "Control_ND_5",
               "34_AVRQ_ArgDiet_Study_2w_Cortex_repet_No143": "Control_AD_3",
               "35_AVRQ_ArgDiet_Study_2w_Cortex_repet_No168": "Proteinuric_ND_3",
               "36_AVRQ_ArgDiet_Study_2w_Cortex_repet_No140": "Control_AD_4",
               "37_AVRQ_ArgDiet_Study_2w_Cortex_repet_No153": "Proteinuric_AD_4",
               "38_AVRQ_ArgDiet_Study_2w_Cortex_repet_No160": "Proteinuric_AD_5",
               "40_AVRQ_ArgDiet_Study_2w_Cortex_repet_No142": "Proteinuric_ND_4",
               "42_AVRQ_ArgDiet_Study_2w_Cortex_repet_No138": "Control_AD_5"},
    "Medulla": {"22_AVRQ_ArgDiet_2w_Med_#147_20240226": "Control_ND_1",
                "23_AVRQ_ArgDiet_2w_Med_#157_20240226": "Proteinuric_AD_1",
                "24_AVRQ_ArgDiet_2w_Med_#146_20240226": "Control_ND_2",
                "25_AVRQ_ArgDiet_2w_Med_#166_20240226": "Proteinuric_ND_1",
                "26_AVRQ_ArgDiet_2w_Med_#148_20240226": "Proteinuric_AD_2",
                "27_AVRQ_ArgDiet_2w_Med_#141_20240226": "Control_AD_1",
                "28_AVRQ_ArgDiet_2w_Med_#164_20240226": "Control_ND_3",
                "29_AVRQ_ArgDiet_2w_Med_#154_20240226": "Control_AD_2",
                "30_AVRQ_ArgDiet_2w_Med_#145_20240226": "Control_ND_4",
                "31_AVRQ_ArgDiet_2w_Med_#167_20240226": "Proteinuric_ND_2",
                "32_AVRQ_ArgDiet_2w_Med_#150_20240226": "Proteinuric_AD_3",
                "33_AVRQ_ArgDiet_2w_Med_#144_20240226": "Control_ND_5",
                "34_AVRQ_ArgDiet_2w_Med_#143_20240226": "Control_AD_3",
                "35_AVRQ_ArgDiet_2w_Med_#168_20240226": "Proteinuric_ND_3",
                "36_AVRQ_ArgDiet_2w_Med_#140_20240226": "Control_AD_4",
                "37_AVRQ_ArgDiet_2w_Med_#153_20240226": "Control_ND_6",
                "38_AVRQ_ArgDiet_2w_Med_#160_20240226": "Proteinuric_AD_4",
                "39_AVRQ_ArgDiet_2w_Med_#139_20240226": "Proteinuric_AD_5",
                "40_AVRQ_ArgDiet_2w_Med_#142_20240226": "Control_AD_5",
                "41_AVRQ_ArgDiet_2w_Med_#161_20240226": "Proteinuric_ND_4",
                "42_AVRQ_ArgDiet_2w_Med_#138_20240226": "Control_AD_6"},
    "Liver": {"43_AVRQ_ArgDiet_2w_Liver_#147_20240226": "Proteinuric_AD_1",
              "44_AVRQ_ArgDiet_2w_Liver_#157_20240226": "Proteinuric_AD_2",
              "45_AVRQ_ArgDiet_2w_Liver_#146_20240226": "Proteinuric_AD_3",
              "46_AVRQ_ArgDiet_2w_Liver_#166_20240226": "Control_AD_1",
              "47_AVRQ_ArgDiet_2w_Liver_#148_20240226": "Control_AD_2",
              "48_AVRQ_ArgDiet_2w_Liver_#141_20240226": "Control_AD_3",
              "49_AVRQ_ArgDiet_2w_Liver_#164_20240226": "Control_AD_4",
              "50_AVRQ_ArgDiet_2w_Liver_#154_20240226": "Control_AD_5",
              "51_AVRQ_ArgDiet_2w_Liver_#145_20240226": "Control_AD_6",
              "52_AVRQ_ArgDiet_2w_Liver_#167_20240226": "Proteinuric_AD_4",
              "53_AVRQ_ArgDiet_2w_Liver_#150_20240226": "Proteinuric_AD_5",
              "54_AVRQ_ArgDiet_2w_Liver_#144_20240226": "Proteinuric_ND_1",
              "55_AVRQ_ArgDiet_2w_Liver_#143_20240226": "Proteinuric_ND_2",
              "56_AVRQ_ArgDiet_2w_Liver_#168_20240226": "Proteinuric_ND_3",
              "57_AVRQ_ArgDiet_2w_Liver_#140_20240226": "Proteinuric_ND_4",
              "58_AVRQ_ArgDiet_2w_Liver_#153_20240226": "Control_ND_1",
              "59_AVRQ_ArgDiet_2w_Liver_#160_20240226": "Control_ND_2",
              "60_AVRQ_ArgDiet_2w_Liver_#139_20240226": "Control_ND_3",
              "61_AVRQ_ArgDiet_2w_Liver_#142_20240226": "Control_ND_4",
              "62_AVRQ_ArgDiet_2w_Liver_#161_20240226": "Control_ND_5",
              "63_AVRQ_ArgDiet_2w_Liver_#138_20240226": "Control_ND_6"},
    "Heart": {"1_AVRQ_ArgDiet_Study_2w_Heart_repet_No150": "Proteinuric_AD_1",
              "2_AVRQ_ArgDiet_Study_2w_Heart_repet_No139": "Proteinuric_AD_2",
              "3_AVRQ_ArgDiet_Study_2w_Heart_repet_No157": "Proteinuric_AD_3",
              "4_AVRQ_ArgDiet_Study_2w_Heart_repet_No142": "Control_AD_1",
              "5_AVRQ_ArgDiet_Study_2w_Heart_repete_No143": "Control_AD_2",
              "6_AVRQ_ArgDiet_Study_2w_Heart_repet_No154": "Control_AD_3",
              "7_AVRQ_ArgDiet_Study_2w_Heart_repet_No138": "Control_AD_4",
              "8_AVRQ_ArgDiet_Study_2w_Heart_repett_No140": "Control_AD_5",
              "9_AVRQ_ArgDiet_Study_2w_Heart_repet_No141": "Control_AD_6",
              "10_AVRQ_ArgDiet_Study_2w_Heart_repet_No148": "Proteinuric_AD_4",
              "11_AVRQ_ArgDiet_Study_2w_Heart_repet_No160": "Proteinuric_AD_5",
              "12_AVRQ_ArgDiet_Study_2w_Heart_repet_No166": "Proteinuric_ND_1",
              "13_AVRQ_ArgDiet_Study_2w_Heart_repet_No167": "Proteinuric_ND_2",
              "14_AVRQ_ArgDiet_Study_2w_Heart_repet_No168": "Proteinuric_ND_3",
              "15_AVRQ_ArgDiet_Study_2w_Heart_repet_No161": "Proteinuric_ND_4",
              "16_AVRQ_ArgDiet_Study_2w_Heart_repet_No147": "Control_ND_1",
              "17_AVRQ_ArgDiet_Study_2w_Heart_repet_No153": "Control_ND_2",
              "18_AVRQ_ArgDiet_Study_2w_Heart_repet_No164": "Control_ND_3",
              "19_AVRQ_ArgDiet_Study_2w_Heart_repet_No144": "Control_ND_4",
              "20_AVRQ_ArgDiet_Study_2w_Heart_repet_No145": "Control_ND_5",
              "21_AVRQ_ArgDiet_Study_2w_Heart_repet_No146": "Control_ND_6"},
    "Muscle": {"22_AVRQ_2kDiet_Study_Muscle_No150": "Proteinuric_AD_1",
               "23_AVRQ_2kDiet_Study_Muscle_No139": "Proteinuric_AD_2",
               "24_AVRQ_2kDiet_Study_Muscle_No157": "Proteinuric_AD_3",
               "25_AVRQ_2kDiet_Study_Muscle_No142": "Control_AD_1",
               "26_AVRQ_2kDiet_Study_Muscle_No143": "Control_AD_2",
               "27_AVRQ_2kDiet_Study_Muscle_No154": "Control_AD_3",
               "28_AVRQ_2kDiet_Study_Muscle_No138": "Control_AD_4",
               "29_AVRQ_2kDiet_Study_Muscle_No140": "Control_AD_5",
               "30_AVRQ_2kDiet_Study_Muscle_No141": "Control_AD_6",
               "31_AVRQ_2kDiet_Study_Muscle_No148": "Proteinuric_AD_4",
               "32_AVRQ_2kDiet_Study_Muscle_No160": "Proteinuric_AD_5",
               "33_AVRQ_2kDiet_Study_Muscle_No166": "Proteinuric_ND_1",
               "34_AVRQ_2kDiet_Study_Muscle_No167": "Proteinuric_ND_2",
               "35_AVRQ_2kDiet_Study_Muscle_No168": "Proteinuric_ND_3",
               "36_AVRQ_2kDiet_Study_Muscle_No161": "Proteinuric_ND_4",
               "37_AVRQ_2kDiet_Study_Muscle_No147": "Control_ND_1",
               "38_AVRQ_2kDiet_Study_Muscle_No153": "Control_ND_2",
               "39_AVRQ_2kDiet_Study_Muscle_No164": "Control_ND_3",
               "40_AVRQ_2kDiet_Study_Muscle_No144": "Control_ND_4",
               "41_AVRQ_2kDiet_Study_Muscle_No145": "Control_ND_5",
               "42_AVRQ_2kDiet_Study_Muscle_No146": "Control_ND_6"},
    "Brain": {"1_AVRQ_ArgDiet_2w_Brain_#147_rep_20240604": "Proteinuric_AD_1",
              "2_AVRQ_ArgDiet_2w_Brain_#157_rep_20240604": "Proteinuric_AD_2",
              "3_AVRQ_ArgDiet_2w_Brain_#146_rep_20240604": "Proteinuric_AD_3",
              "4_AVRQ_ArgDiet_2w_Brain_#166_rep_20240604": "Control_AD_1",
              "5_AVRQ_ArgDiet_2w_Brain_#148_rep_20240604": "Control_AD_2",
              "6_AVRQ_ArgDiet_2w_Brain_#141_rep_20240604": "Control_AD_3",
              "7_AVRQ_ArgDiet_2w_Brain_#164_rep_20240604": "Control_AD_4",
              "8_AVRQ_ArgDiet_2w_Brain_#154_rep_20240604": "Control_AD_5",
              "9_AVRQ_ArgDiet_2w_Brain_#145_rep_20240604": "Control_AD_6",
              "10_AVRQ_ArgDiet_2w_Brain_#167_rep_20240604": "Proteinuric_AD_4",
              "11_AVRQ_ArgDiet_2w_Brain_#150_rep_20240604": "Proteinuric_AD_5",
              "12_AVRQ_ArgDiet_2w_Brain_#144_rep_20240604": "Proteinuric_ND_1",
              "13_AVRQ_ArgDiet_2w_Brain_#143_rep_20240604": "Proteinuric_ND_2",
              "14_AVRQ_ArgDiet_2w_Brain_#168_rep_20240604": "Proteinuric_ND_3",
              "15_AVRQ_ArgDiet_2w_Brain_#140_rep_20240604": "Proteinuric_ND_4",
              "16_AVRQ_ArgDiet_2w_Brain_#153_rep_20240604": "Control_ND_1",
              "17_AVRQ_ArgDiet_2w_Brain_#160_rep_20240604": "Control_ND_2",
              "18_AVRQ_ArgDiet_2w_Brain_#139_rep_20240604": "Control_ND_3",
              "19_AVRQ_ArgDiet_2w_Brain_#142_rep_20240604": "Control_ND_4",
              "20_AVRQ_ArgDiet_2w_Brain_#161_rep_20240604": "Control_ND_5",
              "21_AVRQ_ArgDiet_2w_Brain_#138_rep_20240604": "Control_ND_6"},
    "Lung": {"64_AVRQ_ArgDiet_2w_Lung_#147_20240226": "Proteinuric_AD_1",
             "65_AVRQ_ArgDiet_2w_Lung_#157_20240226": "Proteinuric_AD_2",
             "66_AVRQ_ArgDiet_2w_Lung_#146_20240226": "Proteinuric_AD_3",
             "67_AVRQ_ArgDiet_2w_Lung_#166_20240226": "Control_AD_1",
             "68_AVRQ_ArgDiet_2w_Lung_#148_20240226": "Control_AD_2",
             "69_AVRQ_ArgDiet_2w_Lung_#141_20240226": "Control_AD_3",
             "70_AVRQ_ArgDiet_2w_Lung_#164_20240226": "Control_AD_4",
             "71_AVRQ_ArgDiet_2w_Lung_#154_20240226": "Control_AD_5",
             "72_AVRQ_ArgDiet_2w_Lung_#145_20240226": "Control_AD_6",
             "73_AVRQ_ArgDiet_2w_Lung_#167_20240226": "Proteinuric_AD_4",
             "74_AVRQ_ArgDiet_2w_Lung_#150_20240226_20240306151522": "Proteinuric_AD_5",
             "75_AVRQ_ArgDiet_2w_Lung_#144_20240226": "Proteinuric_ND_1",
             "76_AVRQ_ArgDiet_2w_Lung_#143_20240226": "Proteinuric_ND_2",
             "77_AVRQ_ArgDiet_2w_Lung_#168_20240226": "Proteinuric_ND_3",
             "78_AVRQ_ArgDiet_2w_Lung_#140_20240226": "Proteinuric_ND_4",
             "79_AVRQ_ArgDiet_2w_Lung_#153_20240226": "Control_ND_1",
             "80_AVRQ_ArgDiet_2w_Lung_#160_20240226": "Control_ND_2",
             "81_AVRQ_ArgDiet_2w_Lung_#139_20240226": "Control_ND_3",
             "82_AVRQ_ArgDiet_2w_Lung_#142_20240226": "Control_ND_4",
             "83_AVRQ_ArgDiet_2w_Lung_#161_20240226": "Control_ND_5",
             "84_AVRQ_ArgDiet_2w_Lung_#138_20240226": "Control_ND_6"},
    "Ileum": {"AVRQ_2kDiet_Study_Ileum_No147": "Control_ND_1",
              "AVRQ_2kDiet_Study_Ileum_No153": "Control_ND_2",
              "AVRQ_2kDiet_Study_Ileum_No164": "Control_ND_3",
              "AVRQ_2kDiet_Study_Ileum_No144": "Control_ND_4",
              "AVRQ_2kDiet_Study_Ileum_No145": "Control_ND_5",
              "AVRQ_2kDiet_Study_Ileum_No146": "Control_ND_6",
              "AVRQ_2kDiet_Study_Ileum_No166": "Proteinuric_ND_1",
              "AVRQ_2kDiet_Study_Ileum_No167": "Proteinuric_ND_2",
              "AVRQ_2kDiet_Study_Ileum_No168": "Proteinuric_ND_3",
              "AVRQ_2kDiet_Study_Ileum_No161": "Proteinuric_ND_4",
              "AVRQ_2kDiet_Study_Ileum_No142": "Control_AD_1",
              "AVRQ_2kDiet_Study_Ileum_No143": "Control_AD_2",
              "AVRQ_2kDiet_Study_Ileum_No154": "Control_AD_3",
              "AVRQ_2kDiet_Study_Ileum_No138": "Control_AD_4",
              "AVRQ_2kDiet_Study_Ileum_No140": "Control_AD_5",
              "AVRQ_2kDiet_Study_Ileum_No141": "Control_AD_6",
              "AVRQ_2kDiet_Study_Ileum_No150": "Proteinuric_AD_1",
              "AVRQ_2kDiet_Study_Ileum_No139": "Proteinuric_AD_2",
              "AVRQ_2kDiet_Study_Ileum_No157": "Proteinuric_AD_3",
              "AVRQ_2kDiet_Study_Ileum_No148": "Proteinuric_AD_4",
              "AVRQ_2kDiet_Study_Ileum_No160": "Proteinuric_AD_5"},
    "Jejunum": {"1_AVRQ_ArgDiet_2w_Jejenum_#150": "Proteinuric_AD_1",
                "2_AVRQ_ArgDiet_2w_Jejenum_#139": "Proteinuric_AD_2",
                "3_AVRQ_ArgDiet_2w_Jejenum_#157": "Proteinuric_AD_3",
                "4_AVRQ_ArgDiet_2w_Jejenum_#142": "Control_AD_1",
                "5_AVRQ_ArgDiet_2w_Jejenum_#143": "Control_AD_2",
                "6_AVRQ_ArgDiet_2w_Jejenum_#154": "Control_AD_3",
                "7_AVRQ_ArgDiet_2w_Jejenum_#138": "Control_AD_4",
                "8_AVRQ_ArgDiet_2w_Jejenum_#140": "Control_AD_5",
                "9_AVRQ_ArgDiet_2w_Jejenum_#141": "Control_AD_6",
                "10_AVRQ_ArgDiet_2w_Jejenum_#148": "Proteinuric_AD_4",
                "11_AVRQ_ArgDiet_2w_Jejenum_#160": "Proteinuric_AD_5",
                "12_AVRQ_ArgDiet_2w_Jejenum_#166": "Proteinuric_ND_1",
                "13_AVRQ_ArgDiet_2w_Jejenum_#167": "Proteinuric_ND_2",
                "14_AVRQ_ArgDiet_2w_Jejenum_#168": "Proteinuric_ND_3",
                "15_AVRQ_ArgDiet_2w_Jejenum_#161": "Proteinuric_ND_4",
                "16_AVRQ_ArgDiet_2w_Jejenum_#147": "Control_ND_1",
                "17_AVRQ_ArgDiet_2w_Jejenum_#153": "Control_ND_2",
                "18_AVRQ_ArgDiet_2w_Jejenum_#164": "Control_ND_3",
                "19_AVRQ_ArgDiet_2w_Jejenum_#144": "Control_ND_4",
                "20_AVRQ_ArgDiet_2w_Jejenum_#145": "Control_ND_5",
                "21_AVRQ_ArgDiet_2w_Jejenum_#146": "Control_ND_6"},
    "Duodenum": {"1_AVRQ_ArgDiet_2w_Duodenum_#150": "Proteinuric_AD_1",
                 "2_AVRQ_ArgDiet_2w_Duodenum_#139": "Proteinuric_AD_2",
                 "3_AVRQ_ArgDiet_2w_Duodenum_#157": "Proteinuric_AD_3",
                 "4_AVRQ_ArgDiet_2w_Duodenum_#142": "Control_AD_1",
                 "5_AVRQ_ArgDiet_2w_Duodenum_#143": "Control_AD_2",
                 "6_AVRQ_ArgDiet_2w_Duodenum_#154": "Control_AD_3",
                 "7_AVRQ_ArgDiet_2w_Duodenum_#138": "Control_AD_4",
                 "8_AVRQ_ArgDiet_2w_Duodenum_#140": "Control_AD_5",
                 "9_AVRQ_ArgDiet_2w_Duodenum_#141": "Control_AD_6",
                 "10_AVRQ_ArgDiet_2w_Duodenum_#148": "Proteinuric_AD_4",
                 "11_AVRQ_ArgDiet_2w_Duodenum_#160": "Proteinuric_AD_5",
                 "12_AVRQ_ArgDiet_2w_Duodenum_#166": "Proteinuric_ND_1",
                 "13_AVRQ_ArgDiet_2w_Duodenum_#167": "Proteinuric_ND_2",
                 "14_AVRQ_ArgDiet_2w_Duodenum_#168": "Proteinuric_ND_3",
                 "15_AVRQ_ArgDiet_2w_Duodenum_#161": "Proteinuric_ND_4",
                 "16_AVRQ_ArgDiet_2w_Duodenum_#147": "Control_ND_1",
                 "17_AVRQ_ArgDiet_2w_Duodenum_#153": "Control_ND_2",
                 "18_AVRQ_ArgDiet_2w_Duodenum_#164": "Control_ND_3",
                 "19_AVRQ_ArgDiet_2w_Duodenum_#144": "Control_ND_4",
                 "20_AVRQ_ArgDiet_2w_Duodenum_#145": "Control_ND_5",
                 "21_AVRQ_ArgDiet_2w_Duodenum_#146": "Control_ND_6"},
    "Plasma": {"1_AVRQ_ArgDiet_2w_Plasma_#150_20240510": "Proteinuric_AD_1",
               "2_AVRQ_ArgDiet_2w_Plasma_#139_20240510": "Proteinuric_AD_2",
               "3_AVRQ_ArgDiet_2w_Plasma_#157_20240510": "Proteinuric_AD_3",
               "4_AVRQ_ArgDiet_2w_Plasma_#142_20240510": "Control_AD_1",
               "5_AVRQ_ArgDiet_2w_Plasma_#143_20240510": "Control_AD_2",
               "6_AVRQ_ArgDiet_2w_Plasma_#154_20240510": "Control_AD_3",
               "7_AVRQ_ArgDiet_2w_Plasma_#138_20240510": "Control_AD_4",
               "8_AVRQ_ArgDiet_2w_Plasma_#140_20240510": "Control_AD_5",
               "9_AVRQ_ArgDiet_2w_Plasma_#141_20240510": "Control_AD_6",
               "10_AVRQ_ArgDiet_2w_Plasma_#148_20240510": "Proteinuric_AD_4",
               "11_AVRQ_ArgDiet_2w_Plasma_#160_20240510": "Proteinuric_AD_5",
               "12_AVRQ_ArgDiet_2w_Plasma_#166_20240510": "Proteinuric_ND_1",
               "13_AVRQ_ArgDiet_2w_Plasma_#167_20240510": "Proteinuric_ND_2",
               "14_AVRQ_ArgDiet_2w_Plasma_#168_20240510": "Proteinuric_ND_3",
               "15_AVRQ_ArgDiet_2w_Plasma_#161_20240510": "Proteinuric_ND_4",
               "16_AVRQ_ArgDiet_2w_Plasma_#147_20240510": "Control_ND_1",
               "17_AVRQ_ArgDiet_2w_Plasma_#153_20240510": "Control_ND_2",
               "18_AVRQ_ArgDiet_2w_Plasma_#164_20240510": "Control_ND_3",
               "19_AVRQ_ArgDiet_2w_Plasma_#144_20240510": "Control_ND_4",
               "20_AVRQ_ArgDiet_2w_Plasma_#145_20240510": "Control_ND_5",
               "21_AVRQ_ArgDiet_2w_Plasma_#146_20240510": "Control_ND_6"},
    "Urine": {"1_AVRQ_ArgDiet_2w_Urine_#150_20240513": "Proteinuric_AD_1",
              "2_AVRQ_ArgDiet_2w_Urine_#168_20240513": "Proteinuric_AD_2",
              "3_AVRQ_ArgDiet_2w_Urine_#139_20240513": "Proteinuric_AD_3",
              "4_AVRQ_ArgDiet_2w_Urine_#153_20240513": "Proteinuric_AD_4",
              "5_AVRQ_ArgDiet_2w_Urine_#142_20240513": "Control_AD_1",
              "6_AVRQ_ArgDiet_2w_Urine_#164_20240513": "Control_AD_2",
              "7_AVRQ_ArgDiet_2w_Urine_#143_20240513": "Control_AD_3",
              "8_AVRQ_ArgDiet_2w_Urine_#144_20240513": "Control_AD_4",
              "9_AVRQ_ArgDiet_2w_Urine_#154_20240513": "Control_AD_5",
              "10_AVRQ_ArgDiet_2w_Urine_#146_20240513": "Control_AD_6",
              "11_AVRQ_ArgDiet_2w_Urine_#138_20240513": "Control_ND_1",
              "12_AVRQ_ArgDiet_2w_Urine_#140_20240513": "Control_ND_2",
              "13_AVRQ_ArgDiet_2w_Urine_#141_20240513": "Control_ND_3",
              "14_AVRQ_ArgDiet_2w_Urine_#148_20240513": "Control_ND_4",
              "15_AVRQ_ArgDiet_2w_Urine_#160_20240513": "Proteinuric_ND_1",
              "16_AVRQ_ArgDiet_2w_Urine_#166_20240513": "Proteinuric_ND_2",
              "17_AVRQ_ArgDiet_2w_Urine_#167_20240513": "Proteinuric_ND_3"},
    "Erythrocytes": {"1_AVRQ_ArgDiet_Study_2w_Erythrocytes_No150": "Proteinuric_AD_1",
                     "2_AVRQ_ArgDiet_Study_2w_Erythrocytes_No139": "Proteinuric_AD_2",
                     "3_AVRQ_ArgDiet_Study_2w_Erythrocytes_No157": "Proteinuric_AD_3",
                     "4_AVRQ_ArgDiet_Study_2w_Erythrocytes_No142": "Control_AD_1",
                     "5_AVRQ_ArgDiet_Study_2w_Erythrocytes_No143": "Control_AD_2",
                     "6_AVRQ_ArgDiet_Study_2w_Erythrocytes_No154": "Control_AD_3",
                     "7_AVRQ_ArgDiet_Study_2w_Erythrocytes_No138": "Control_AD_4",
                     "8_AVRQ_ArgDiet_Study_2w_Erythrocytes_No140": "Control_AD_5",
                     "9_AVRQ_ArgDiet_Study_2w_Erythrocytes_No141": "Control_AD_6",
                     "10_AVRQ_ArgDiet_Study_2w_Erythrocytes_No148": "Proteinuric_AD_4",
                     "11_AVRQ_ArgDiet_Study_2w_Erythrocytes_No160": "Proteinuric_AD_5",
                     "12_AVRQ_ArgDiet_Study_2w_Erythrocytes_No166": "Proteinuric_ND_1",
                     "13_AVRQ_ArgDiet_Study_2w_Erythrocytes_No167": "Proteinuric_ND_2",
                     "14_AVRQ_ArgDiet_Study_2w_Erythrocytes_No168": "Proteinuric_ND_3",
                     "15_AVRQ_ArgDiet_Study_2w_Erythrocytes_No161": "Proteinuric_ND_4",
                     "16_AVRQ_ArgDiet_Study_2w_Erythrocytes_No147": "Control_ND_1",
                     "17_AVRQ_ArgDiet_Study_2w_Erythrocytes_No153": "Control_ND_2",
                     "18_AVRQ_ArgDiet_Study_2w_Erythrocytes_No164": "Control_ND_3",
                     "19_AVRQ_ArgDiet_Study_2w_Erythrocytes_No144": "Control_ND_4",
                     "20_AVRQ_ArgDiet_Study_2w_Erythrocytes_No145": "Control_ND_5",
                     "21_AVRQ_ArgDiet_Study_2w_Erythrocytes_No146": "Control_ND_6"},
    "Colon": {"1_AVRQ_ArgDiet_Study_2w_Colon_repet_No150": "Proteinuric_AD_1",
              "2_AVRQ_ArgDiet_Study_2w_Colon_repet_No139": "Proteinuric_AD_2",
              "3_AVRQ_ArgDiet_Study_2w_Colon_repet_No157": "Proteinuric_AD_3",
              "4_AVRQ_ArgDiet_Study_2w_Colon_repet_No142": "Control_AD_1",
              "5_AVRQ_ArgDiet_Study_2w_Colon_repet_No143": "Control_AD_2",
              "6_AVRQ_ArgDiet_Study_2w_Colon_repet_No154": "Control_AD_3",
              "7_AVRQ_ArgDiet_Study_2w_Colon_repet_No138": "Control_AD_4",
              "8_AVRQ_ArgDiet_Study_2w_Colon_repet_No140": "Control_AD_5",
              "9_AVRQ_ArgDiet_Study_2w_Colon_repet_No141": "Control_AD_6",
              "10_AVRQ_ArgDiet_Study_2w_Colon_repet_No148": "Proteinuric_AD_4",
              "11_AVRQ_ArgDiet_Study_2w_Colon_repet_No160": "Proteinuric_AD_5",
              "12_AVRQ_ArgDiet_Study_2w_Colon_repet_No166": "Proteinuric_ND_1",
              "13_AVRQ_ArgDiet_Study_2w_Colon_repet_No167": "Proteinuric_ND_2",
              "14_AVRQ_ArgDiet_Study_2w_Colon_repet_No168": "Proteinuric_ND_3",
              "15_AVRQ_ArgDiet_Study_2w_Colon_repet_No161": "Proteinuric_ND_4",
              "16_AVRQ_ArgDiet_Study_2w_Colon_repet_No147": "Control_ND_1",
              "17_AVRQ_ArgDiet_Study_2w_Colon_repet_No153": "Control_ND_2",
              "18_AVRQ_ArgDiet_Study_2w_Colon_repet_No164": "Control_ND_3",
              "19_AVRQ_ArgDiet_Study_2w_Colon_repet_No144": "Control_ND_4",
              "20_AVRQ_ArgDiet_Study_2w_Colon_repet_No145": "Control_ND_5",
              "21_AVRQ_ArgDiet_Study_2w_Colon_repet_No146": "Control_ND_6"},
    "Spleen": {"22_AVRQ_ArgDiet_Study_2w_Spleen_repet_No150": "Proteinuric_AD_1",
               "23_AVRQ_ArgDiet_Study_2w_Spleen_repet_No139": "Proteinuric_AD_2",
               "24_AVRQ_ArgDiet_Study_2w_Spleen_repet_No157": "Proteinuric_AD_3",
               "25_AVRQ_ArgDiet_Study_2w_Spleen_repet_No142": "Control_AD_1",
               "26_AVRQ_ArgDiet_Study_2w_Spleen_repet_No143": "Control_AD_2",
               "27_AVRQ_ArgDiet_Study_2w_Spleen_repet_No154": "Control_AD_3",
               "28_AVRQ_ArgDiet_Study_2w_Spleen_repet_No138": "Control_AD_4",
               "29_AVRQ_ArgDiet_Study_2w_Spleen_repet_No140": "Control_AD_5",
               "30_AVRQ_ArgDiet_Study_2w_Spleen_repet_No141": "Control_AD_6",
               "31_AVRQ_ArgDiet_Study_2w_Spleen_repet_No148": "Proteinuric_AD_4",
               "32_AVRQ_ArgDiet_Study_2w_Spleen_repet_No160": "Proteinuric_AD_5",
               "33_AVRQ_ArgDiet_Study_2w_Spleen_repet_No166": "Proteinuric_ND_1",
               "34_AVRQ_ArgDiet_Study_2w_Spleen_repet_No167": "Proteinuric_ND_2",
               "35_AVRQ_ArgDiet_Study_2w_Spleen_repet_No168": "Proteinuric_ND_3",
               "36_AVRQ_ArgDiet_Study_2w_Spleen_repet_No161": "Proteinuric_ND_4",
               "37_AVRQ_ArgDiet_Study_2w_Spleen_repet_No147": "Control_ND_1",
               "38_AVRQ_ArgDiet_Study_2w_Spleen_repet_No153": "Control_ND_2",
               "39_AVRQ_ArgDiet_Study_2w_Spleen_repet_No164": "Control_ND_3",
               "40_AVRQ_ArgDiet_Study_2w_Spleen_repet_No144": "Control_ND_4",
               "41_AVRQ_ArgDiet_Study_2w_Spleen_repet_No145": "Control_ND_5",
               "42_AVRQ_ArgDiet_Study_2w_Spleen_repet_No146": "Control_ND_6"}}

def find_input_file(input_dir, file_stem):
    """Locate one deposited file by exact stem, accepting CSV or TSV."""
    candidates = [input_dir / f"{file_stem}.csv", input_dir / f"{file_stem}.tsv"]
    existing = [path for path in candidates if path.exists()]

    if len(existing) == 0:
        raise FileNotFoundError(f"Could not find '{file_stem}.csv' or '{file_stem}.tsv' "
                                f"in {input_dir.resolve()}")
    if len(existing) > 1:
        raise RuntimeError(f"Both CSV and TSV versions were found for '{file_stem}'. "
                           "Keep only the deposited file that should be used for this analysis.")
    return existing[0]

def read_spectronaut_file(file_path):
    """Read a Spectronaut CSV or TSV export."""
    suffix = file_path.suffix.lower()

    if suffix == ".tsv":
        return pd.read_csv(file_path, sep="\t", low_memory=False)
    if suffix == ".csv":
        return pd.read_csv(file_path, low_memory=False)
    raise ValueError(f"Unsupported file type: {file_path}")

def infer_sample_group(short_name):
    """Convert mapped sample names to Control_ND, Proteinuric_ND, etc."""
    if pd.isna(short_name):
        return np.nan

    match = re.match(
        r"^(Control|Proteinuric)_(ND|AD)_\d+$", str(short_name))

    if match:
        return f"{match.group(1)}_{match.group(2)}"

    return np.nan

def load_spectronaut_exports():
    """Load the 15 deposited exports and apply the manuscript sample map."""
    required_columns = ["R.FileName", "R.Condition", *CHANNEL_COLUMNS]
    dataframes = []
    unmapped_rows = []
    print("\nInput files used:")

    for organ in ORGAN_ORDER:
        file_path = find_input_file(INPUT_DIR, EXPECTED_FILES[organ])
        temp_df = read_spectronaut_file(file_path)

        missing_columns = [column for column in required_columns if column not in temp_df.columns]

        if missing_columns:
            raise ValueError(f"{file_path.name} is missing required columns: {missing_columns}")

        temp_df = temp_df[required_columns].copy()
        temp_df["Organ"] = organ
        mapping = SAMPLE_MAP.get(organ, {})
        temp_df["MappedSample"] = temp_df["R.FileName"].map(mapping)
        missing_map = (temp_df.loc[temp_df["MappedSample"].isna(), ["Organ", "R.FileName"]].drop_duplicates())

        if not missing_map.empty:
            unmapped_rows.append(missing_map)

        temp_df["SampleGroup"] = temp_df["MappedSample"].apply(infer_sample_group)
        dataframes.append(temp_df)
        print(f"  {organ}: {file_path.name}")

    if unmapped_rows:
        unmapped = pd.concat(unmapped_rows, ignore_index=True)
        raise ValueError("Unmapped R.FileName values were detected. The statistical analysis "
                         "requires the exact Supplementary Table 4 sample mapping.\n\n"
                         + unmapped.to_string(index=False))

    combined = pd.concat(dataframes, ignore_index=True)

    if combined["SampleGroup"].isna().any():
        bad = (combined.loc[combined["SampleGroup"].isna(), ["Organ", "R.FileName", "MappedSample"]].drop_duplicates())
        raise ValueError("Some mapped sample names could not be assigned to ND/AD groups:\n\n"
                         + bad.to_string(index=False))

    return combined


def calculate_plot_results(protein_data):
    """Calculate the values used for the stacked-bar plots."""
    results = {}

    for group_name, conditions in EXPERIMENTAL_GROUPS.items():
        intensity_rows = []
        percentage_rows = []

        for organ in ORGAN_ORDER:
            organ_data = protein_data.loc[(protein_data["Organ"] == organ) &
                                          (protein_data["R.Condition"].isin(conditions))].copy()

            if organ_data.empty:
                print(f"Warning: no rows found for {group_name} | {organ}")
                continue

            # Same calculation as the original plotting script:
            # sum protein-group channel intensities per mouse, then average mice.
            per_mouse_sum = (organ_data.groupby(["R.Condition", "R.FileName"])[CHANNEL_COLUMNS].sum())
            mean_channel_intensity = (per_mouse_sum.groupby(level="R.Condition").mean().reindex(conditions))
            sample_counts = (per_mouse_sum.groupby(level="R.Condition").size().reindex(conditions))
            totals = mean_channel_intensity.sum(axis=1)
            channel_percentage = mean_channel_intensity.div(totals, axis=0) * 100

            for condition in conditions:
                if (condition not in mean_channel_intensity.index
                        or mean_channel_intensity.loc[condition].isna().all()):
                    continue

                intensity_row = {"Organ": organ, "Group": group_name, "Condition": condition,
                                 "N": sample_counts.loc[condition]}
                percentage_row = {"Organ": organ, "Group": group_name, "Condition": condition,
                                  "N": sample_counts.loc[condition]}

                for channel in CHANNEL_COLUMNS:
                    intensity_row[channel] = mean_channel_intensity.loc[condition, channel]
                    percentage_row[channel] = channel_percentage.loc[condition, channel]
                intensity_rows.append(intensity_row)
                percentage_rows.append(percentage_row)
        results[group_name] = {"intensity": pd.DataFrame(intensity_rows), "percentage": pd.DataFrame(percentage_rows)}

    return results


def calculate_statistics(protein_data):
    """Reproduce Supplementary Table statistics.

    For each organ, channel intensities are summed per mapped mouse. Control and Proteinuric groups are compared
    separately in ND and AD using two-sided Welch's t-tests.
    BH-FDR is then applied across ALL organ x Channel 2/3 tests together, separately within ND and AD."""
    rows = []
    per_mouse_rows = []

    for organ in ORGAN_ORDER:
        organ_data = protein_data.loc[protein_data["Organ"] == organ].copy()

        per_mouse = (organ_data.groupby("MappedSample")[STAT_CHANNELS].sum().reset_index())
        per_mouse["SampleGroup"] = per_mouse["MappedSample"].apply(infer_sample_group)
        per_mouse.insert(0, "Organ", organ)
        per_mouse_rows.append(per_mouse.copy())

        for state, control_label, proteinuric_label in COMPARISON_PAIRS:
            sub = per_mouse.loc[per_mouse["SampleGroup"].isin([control_label, proteinuric_label])].copy()

            for channel in STAT_CHANNELS:
                control_values = (sub.loc[sub["SampleGroup"] == control_label, channel].astype(float).dropna())
                proteinuric_values = (sub.loc[sub["SampleGroup"] == proteinuric_label, channel].astype(float).dropna())
                n_control = len(control_values)
                n_proteinuric = len(proteinuric_values)

                if n_control == 0 or n_proteinuric == 0:
                    continue

                mean_control = control_values.mean()
                mean_proteinuric = proteinuric_values.mean()
                sd_control = (control_values.std(ddof=1) if n_control > 1 else np.nan)
                sd_proteinuric = (proteinuric_values.std(ddof=1) if n_proteinuric > 1 else np.nan)
                sem_control = (sd_control / np.sqrt(n_control) if n_control > 1 else np.nan)
                sem_proteinuric = (sd_proteinuric / np.sqrt(n_proteinuric) if n_proteinuric > 1 else np.nan)
                log2fc = (np.log2(mean_proteinuric / mean_control) if mean_control > 0 and mean_proteinuric > 0
                          else np.nan)

                if n_control > 1 and n_proteinuric > 1:
                    t_statistic, p_value = ttest_ind(control_values, proteinuric_values, equal_var=False,
                                                     nan_policy="omit")
                else:
                    t_statistic, p_value = np.nan, np.nan

                rows.append({"Organ": organ, "State": state, "Channel": channel,
                             "Channel_label": CHANNEL_LABELS[channel], "ControlGroup": control_label,
                             "ProteinuricGroup": proteinuric_label, "N_control": n_control,
                             "N_proteinuric": n_proteinuric, "Mean_control": mean_control,
                             "Mean_proteinuric": mean_proteinuric, "SD_control": sd_control,
                             "SD_proteinuric": sd_proteinuric, "SEM_control": sem_control,
                             "SEM_proteinuric": sem_proteinuric, "log2FC_proteinuric_vs_control": log2fc,
                             "t_statistic": t_statistic, "p_value": p_value})

    stats_df = pd.DataFrame(rows)
    stats_df["p_adj_BH"] = np.nan

    # Exact FDR family used in Supplementary Table:
    # all organ x Channel 2/3 tests together, separately in ND and AD.
    for state in ["ND", "AD"]:
        mask = ((stats_df["State"] == state) & stats_df["p_value"].notna())

        if mask.any():
            stats_df.loc[mask, "p_adj_BH"] = multipletests(stats_df.loc[mask, "p_value"], method="fdr_bh")[1]

    stats_df["FDR_significant_0.05"] = (stats_df["p_adj_BH"] <= 0.05)
    per_mouse_df = pd.concat(per_mouse_rows, ignore_index=True)

    return stats_df, per_mouse_df

def export_tables(results, stats_df, per_mouse_df):
    """Save derived numerical results as CSV files."""
    stats_df.to_csv(OUTPUT_DIR / "channel2_channel3_statistics_FDR.csv", index=False)

    per_mouse_df.to_csv(OUTPUT_DIR / "per_mouse_channel_intensities_for_statistics.csv", index=False)

    for group_name, group_results in results.items():
        intensity_table = group_results["intensity"].rename(columns=CHANNEL_LABELS)
        percentage_table = group_results["percentage"].rename(columns=CHANNEL_LABELS)

        intensity_table.to_csv(OUTPUT_DIR / f"{group_name}_mean_channel_intensities.csv", index=False)
        percentage_table.to_csv(OUTPUT_DIR / f"{group_name}_channel_percentages.csv", index=False)


def plot_results(results, stats_df):
    """Generate manuscript-style stacked-bar plots."""
    global_max_y = 0

    for group_results in results.values():
        intensity_table = group_results["intensity"]
        if intensity_table.empty:
            continue

        group_max = (intensity_table[CHANNEL_COLUMNS].sum(axis=1).max())
        global_max_y = max(global_max_y, group_max)

    for group_name, group_results in results.items():
        datasets = [(group_results["intensity"], "Intensity", "Protein abundance", "Protein abundance"),
                    (group_results["percentage"], "Percentage", "Protein abundance (%)", "Protein abundance (%)")]

        for table, plot_type, y_label, title_prefix in datasets:
            if table.empty:
                continue

            plotting_rows = []
            hatched_indices = []
            organ_centers = []
            organ_labels = []
            organ_bar_positions = {}
            bar_index = 0

            for organ in ORGAN_ORDER:
                organ_table = table.loc[table["Organ"] == organ]
                if organ_table.empty:
                    continue

                pair_positions = []

                for condition_index, condition in enumerate(EXPERIMENTAL_GROUPS[group_name]):
                    row = organ_table.loc[organ_table["Condition"] == condition]
                    if row.empty:
                        continue

                    row = row.iloc[0]
                    plotting_rows.append(row[CHANNEL_COLUMNS].astype(float).values)
                    pair_positions.append(bar_index)

                    if condition_index == 1:
                        hatched_indices.append(bar_index)
                    bar_index += 1

                if len(pair_positions) == 2:
                    center = (pair_positions[0] + pair_positions[1]) / 2
                    organ_centers.append(center)
                    organ_labels.append(organ)
                    organ_bar_positions[organ] = pair_positions

                # Spacer between tissue pairs.
                plotting_rows.append([0.0] * len(CHANNEL_COLUMNS))
                bar_index += 1

            plot_df = pd.DataFrame(plotting_rows, columns=CHANNEL_COLUMNS).rename(columns=CHANNEL_LABELS)
            fig, ax = plt.subplots(figsize=(12, 4.8))

            plot_df.plot(kind="bar", stacked=True, color=COLORS, width=1.0, ax=ax)

            for container in ax.containers:
                for index in hatched_indices:
                    if index < len(container.patches):
                        container.patches[index].set_hatch("///")

            if plot_type == "Intensity":
                ax.set_ylim(0, global_max_y * 1.30)
            else:
                ax.set_ylim(0, 100)

            ax.set_xticks(organ_centers)
            ax.set_xticklabels(organ_labels, rotation=90)
            ax.set_ylabel(y_label)
            ax.set_xlabel("")
            ax.set_title(f"{title_prefix} - {group_name}")
            channel_legend = ax.legend(title="Label", bbox_to_anchor=(1.02, 1.00), loc="upper left")
            ax.add_artist(channel_legend)
            condition_handles = [Patch(facecolor="white", edgecolor="black", label="Control"),
                                 Patch(facecolor="white", edgecolor="black", hatch="///", label="Proteinuric")]
            ax.legend(handles=condition_handles, title="Condition", bbox_to_anchor=(1.02, 0.55), loc="upper left")

            # Significance stars come from the exact Supplementary Table
            # FDR analysis and are shown only on absolute-abundance plots.
            if plot_type == "Intensity":
                state = PLOT_GROUP_TO_STATE[group_name]
                plot_stats = stats_df.loc[stats_df["State"] == state].copy()

                for organ in organ_labels:
                    if organ not in organ_bar_positions:
                        continue

                    x1, x2 = organ_bar_positions[organ]
                    star_x = (x1 + x2) / 2
                    significant_channels = []

                    for channel in STAT_CHANNELS:
                        stat_row = plot_stats.loc[(plot_stats["Organ"] == organ) & (plot_stats["Channel"] == channel)]

                        if stat_row.empty:
                            continue

                        q_value = stat_row.iloc[0]["p_adj_BH"]
                        if np.isfinite(q_value) and q_value <= 0.05:
                            significant_channels.append(channel)

                    if not significant_channels:
                        continue

                    two_rows = table.loc[table["Organ"] == organ]
                    bar_max = (two_rows[CHANNEL_COLUMNS].sum(axis=1).max())
                    first_star_y = bar_max + global_max_y * 0.035
                    star_spacing = global_max_y * 0.045

                    for star_number, channel in enumerate(significant_channels):
                        ax.text(star_x, first_star_y + star_number * star_spacing, "*", ha="center", va="bottom",
                                fontsize=14, fontweight="bold", color=CHANNEL_STAR_COLORS[channel], clip_on=False)

            plt.tight_layout()
            output_stem = (OUTPUT_DIR / f"{group_name}_{plot_type.lower()}_protein_channel_abundance")
            fig.savefig(f"{output_stem}.pdf", dpi=300, bbox_inches="tight")
            fig.savefig(f"{output_stem}.svg", dpi=300, bbox_inches="tight")
            plt.close(fig)


def main():
    args = parse_arguments()

    global INPUT_DIR, OUTPUT_DIR
    INPUT_DIR = args.input
    OUTPUT_DIR = args.output

    if not INPUT_DIR.exists():
        raise FileNotFoundError(f"Input directory does not exist: {INPUT_DIR.resolve()}")

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    print(f"Input directory: {INPUT_DIR.resolve()}")
    print(f"Output directory: {OUTPUT_DIR.resolve()}")

    protein_data = load_spectronaut_exports()
    results = calculate_plot_results(protein_data)
    stats_df, per_mouse_df = calculate_statistics(protein_data)

    export_tables(results, stats_df, per_mouse_df)

    plot_results(results, stats_df)

    print("\nFDR-significant tests (BH-FDR <= 0.05):")
    significant = stats_df.loc[stats_df["FDR_significant_0.05"],
                               ["Organ", "State", "Channel_label", "p_value", "p_adj_BH"]]

    if significant.empty:
        print("  None")
    else:
        print(significant.to_string(index=False))

    for state in ["ND", "AD"]:
        n_tests = ((stats_df["State"] == state) & stats_df["p_value"].notna()).sum()
        print(f"{state}: BH-FDR correction applied across "
              f"{n_tests} valid tissue x labeled-channel tests.")

    print("\nAnalysis complete.")

if __name__ == "__main__":
    main()