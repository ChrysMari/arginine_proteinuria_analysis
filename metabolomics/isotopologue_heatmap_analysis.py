"""
Isotopologue heatmap analysis for the 2-week heavy-arginine diet experiment.

The script reads the five deposited targeted-metabolomics CSV reports and
compares isotope/light ratios between proteinuric (AVRQ_Heavy) and control
(wtRQ_Heavy) mice across tissues.

For each isotopologue x tissue comparison it:
- calculates group means and the Proteinuric/Control log2 ratio;
- calculates a bootstrap 95% confidence interval for the log2 ratio;
- performs a two-sided Mann-Whitney U test;
- calculates global BH-FDR across all valid comparisons;
- calculates BH-FDR within each parent metabolite;
- generates three otherwise identical dot plots using:
    1. nominal P < 0.05,
    2. global BH-FDR < 0.20,
    3. within-metabolite BH-FDR < 0.05;
- exports a supplementary Excel workbook containing all statistics.
"""

from pathlib import Path
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.collections import PatchCollection
from matplotlib import colors
from matplotlib.lines import Line2D
from scipy.stats import mannwhitneyu

# SETTINGS AND COMMAND-LINE ARGUMENTS

# Significance thresholds
ALPHA = 0.05
GLOBAL_FDR_CUTOFF = 0.20
WITHIN_METABOLITE_FDR_CUTOFF = 0.05
N_BOOTSTRAP = 10000
RANDOM_SEED = 12345

def parse_arguments():
    parser = argparse.ArgumentParser(
        description=("Generate isotopologue dot plots and supplementary statistics "
                     "for control versus proteinuric mice."))

    parser.add_argument("--input", type=Path, default=Path("data/metabolomics"),
                        help="Directory containing the five deposited metabolomics CSV files.")
    parser.add_argument("--output", type=Path, default=Path("results/metabolomics/isotopologue_heatmap"),
                        help="Directory for figures and supplementary tables.")
    return parser.parse_args()

ARGS = parse_arguments()
INPUT_DIR = ARGS.input
OUTPUT_DIR = ARGS.output
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

OUT_XLSX = (OUTPUT_DIR  / "isotopologue_stats_three_methods.xlsx")

# RENAME DICTIONARIES
peptide_rename = {"ADMA": "ADMA", "Agmatine": "Agm", "Arginine": "Arg", "Argininosuccinic acid": "ArgSuc",
                  "Aspartic acid": "Asp", "Citrulline": "Citr", "GABA": "GABA", "Glutamic acid": "Glu",
                  "Glutamine": "Gln", "L-Homoarginine": "HArg", "MMA": "MMA", "N-a-Acetyl-L-arginine": "AcArg",
                  "Na-Acetyl-L-ornithine": "AcOrn", "N-Acetyl-L-glutamic acid": "AcGlu", "Ornithine": "Orn",
                  "Proline": "Pro", "Putrescine": "Putr", "Spermine": "Spn", "trans-4-Hydroxy-L-proline": "4OHPro",
                  "Urea": "Urea", "Citric acid": "CA", "Malic acid": "MA", "Glucose": "Gluc"}

organ_rename = {"Cortex": "Ct", "Medulla": "Md", "Liver": "Lv", "Heart": "He", "Muscle": "Ms", "Brain": "Br",
                "Lung": "Lg", "Spleen": "Sp", "Duodenum": "Dd", "Jejenum": "Jj", "Jejunum": "Jj", "Ileum": "Il",
                "Colon": "Cl", "Erythrocytes": "Er", "Plasma": "Pl", "Urine": "Ur"}

peptides_to_exclude = ["Guanidineacetic acid"]
peptides_to_include = ["Arg_13C10", "Orn_13C7", "Urea_13C3", "Citr_13C7", "ArgSuc_13C7", "Arg_13C7", "ADMA_13C7",
                       "Citr_13C9", "Arg_13C9", "MMA_13C9", "Pro_13C6", "Gln_13C6", "Glu_13C6", "Glu_13C5", "Asp_13C4",
                       "Asp_13C5", "Asp_13C1", "Orn_13C6", "Citr_13C6", "ArgSuc_13C6", "Arg_13C6", "Urea_13C1",
                       "ArgSuc_13C10", "ADMA_13C6", "MMA_13C6", "AcArg_13C6", "Arg_13C8", "MMA_13C8", "Agm_13C9",
                       "AcArg_13C10", "MMA_13C10", "ADMA_13C10", "CA_13C4", "MA_13C4", "Gluc_13C3"]
organs = ["Ct", "Md", "Lv", "He", "Ms", "Br", "Lg", "Sp", "Dd", "Jj", "Il", "Cl", "Er", "Pl", "Ur"]

# DEPOSITED INPUT FILES
CSV_FILES = [
    "AVRQ_ArgDiet_2weeks_Ctx_Med_Liver_Heart_report.csv",
    "AVRQ_ArgDiet_2weeks_Lung_Pancreas_Spleen_Muscle_report.csv",
    "AVRQ_ArgDiet_2weeks_Brain_WAT_Duodenum_Ileum_report.csv",
    "AVRQ_ArgDiet_2weeks_Jejenum_Colon_Erythrocytes_report.csv",
    "AVRQ_ArgDiet_2weeks_TCA_report.csv",
]

def validate_input_files():
    missing = [ filename for filename in CSV_FILES if not (INPUT_DIR / filename).exists()]

    if missing:
        missing_text = "\n".join(f"  - {filename}" for filename in missing)
        raise FileNotFoundError(
            "The following deposited metabolomics files were not found in "
            f"{INPUT_DIR.resolve()}:\n{missing_text}"
        )
validate_input_files()

# HELPER FUNCTIONS
def bh_fdr(pvals):
    """Benjamini-Hochberg false discovery rate correction."""

    p = np.asarray(pvals, dtype=float)
    out = np.full_like(p, np.nan, dtype=float)
    mask = np.isfinite(p)

    if mask.sum() == 0:
        return out

    pv = p[mask]
    n = len(pv)
    order = np.argsort(pv)
    ranked = pv[order]
    q = ranked * n / np.arange(1, n + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.clip(q, 0, 1)
    inv = np.empty_like(order)
    inv[order] = np.arange(n)
    out[mask] = q[inv]

    return out

def bootstrap_log2_ratio_ci(
    control,
    proteinuric,
    n_bootstrap=10000,
    confidence_level=0.95,
    rng=None
):
    """Bootstrap 95% confidence interval for:
        log2(mean isotope/light ratio in Proteinuric mice /  mean isotope/light ratio in Control mice)
    Bootstrap resampling is performed independently within the two biological groups."""

    control = np.asarray(control, dtype=float)
    proteinuric = np.asarray(proteinuric, dtype=float)
    control = control[np.isfinite(control)]
    proteinuric = proteinuric[np.isfinite(proteinuric)]

    if (len(control) == 0 or len(proteinuric) == 0):
        return np.nan, np.nan

    if rng is None:
        rng = np.random.default_rng()

    bootstrap_effects = []

    for _ in range(n_bootstrap):

        ctrl_boot = rng.choice(control, size=len(control), replace=True)
        prot_boot = rng.choice(proteinuric, size=len(proteinuric), replace=True)
        ctrl_mean = np.mean(ctrl_boot)
        prot_mean = np.mean(prot_boot)

        if (ctrl_mean > 0 and prot_mean > 0):

            effect = np.log2(prot_mean / ctrl_mean)
            bootstrap_effects.append(effect)

    if len(bootstrap_effects) == 0:
        return np.nan, np.nan

    bootstrap_effects = np.asarray(bootstrap_effects)
    alpha_ci = (1 - confidence_level)
    lower = np.percentile(bootstrap_effects, 100 * alpha_ci / 2)
    upper = np.percentile(bootstrap_effects, 100 * (1 - alpha_ci / 2))

    return float(lower), float(upper)

def format_effect_ci(effect, lower, upper):

    if (np.isfinite(effect) and np.isfinite(lower) and np.isfinite(upper)):

        return (
            f"{effect:.2f} "
            f"({lower:.2f} to {upper:.2f})"
        )

    return ""

# Reproducible bootstrap
rng = np.random.default_rng(RANDOM_SEED)

# LOAD AND COMBINE FILES
dfs = []
print("Input directory:", INPUT_DIR.resolve())
print("Output directory:", OUTPUT_DIR.resolve())
print("\nInput files used:")

for filename in CSV_FILES:

    file_path = INPUT_DIR / filename
    temp_df = pd.read_csv(file_path, sep=",", low_memory=False)
    temp_df["SourceFile"] = filename
    dfs.append(temp_df)
    print(f"  - {filename}")

df = pd.concat(dfs, ignore_index=True)

# Remove exact duplicate rows
df = df.drop_duplicates()

# CLEAN DATA
df = df[~df["Peptide"].isin(peptides_to_exclude)].copy()
df["Peptide_short"] = (df["Peptide"].map(peptide_rename).fillna(df["Peptide"]))
df["Organ"] = (df["Organ"].map(organ_rename).fillna(df["Organ"]))

# CLEAN ISOTOPE LABEL TYPE
iso = pd.to_numeric(df["IsotopeLabelType"], errors="coerce")
df = df[iso.notna()].copy()
df["IsotopeLabelType"] = (iso[iso.notna()].astype(int))

# CREATE METABOLITE + ISOTOPOLOGUE IDENTIFIER
df["Peptide_Ratio"] = (df["Peptide_short"] + "_13C" + df["IsotopeLabelType"].astype(str))

# FILTER INCLUDED ISOTOPOLOGUES AND ORGANS
df = df[df["Peptide_Ratio"].isin(peptides_to_include)].copy()
df = df[df["Organ"].isin(organs)].copy()

# CONDITIONS
sub_df_proteinuric = df[df["Condition"] == "AVRQ_Heavy"].copy()
sub_df_control = df[df["Condition"] == "wtRQ_Heavy"].copy()

# PREPARE MATRICES
N = len(peptides_to_include)
M = len(organs)
effect_matrix = np.full((N, M), np.nan, dtype=float)
area_matrix = np.full((N, M), np.nan, dtype=float)
p_matrix = np.full((N, M), np.nan, dtype=float)
supp_rows = []

# LOOP THROUGH ISOTOPOLOGUES × ORGANS
for i, peptide in enumerate(peptides_to_include):

    for j, organ in enumerate(organs):

        # Select corresponding samples
        proteinuric_df = sub_df_proteinuric[(sub_df_proteinuric["Peptide_Ratio"] == peptide) &
                                            (sub_df_proteinuric["Organ"] == organ)].copy()

        control_df = sub_df_control[(sub_df_control["Peptide_Ratio"] == peptide) &
                                    (sub_df_control["Organ"] == organ)].copy()

        # Determine isotopologue-specific ratio column
        carbon_num = peptide.split("13C", 1)[1]
        ratio_col = ("Ratio" + carbon_num + "To0")

        # If this ratio column does not exist, skip comparison
        if (ratio_col not in proteinuric_df.columns or ratio_col not in control_df.columns):
            continue

        # Biological-replicate isotope/light ratios
        proteinuric_values = pd.to_numeric(proteinuric_df[ratio_col], errors="coerce").dropna()
        control_values = pd.to_numeric(control_df[ratio_col], errors="coerce").dropna()
        n_proteinuric_mice = int(len(proteinuric_values))
        n_control_mice = int(len(control_values))

        # Group means
        mean_ratio_proteinuric = (proteinuric_values.mean() if len(proteinuric_values) > 0 else np.nan)
        mean_ratio_control = (control_values.mean() if len(control_values) > 0 else np.nan)

        # Proteinuric / Control ratio
        if (np.isfinite(mean_ratio_proteinuric) and np.isfinite(mean_ratio_control) and mean_ratio_control > 0):
            ratio_proteinuric_vs_control = (mean_ratio_proteinuric / mean_ratio_control)
        else:
            ratio_proteinuric_vs_control = np.nan


        # EFFECT SIZE SHOWN BY CIRCLE COLOR
        if (np.isfinite(ratio_proteinuric_vs_control) and ratio_proteinuric_vs_control > 0):
            log2_effect = np.log2(ratio_proteinuric_vs_control)
        else:
            log2_effect = np.nan
        effect_matrix[i, j] = log2_effect

        # DOT SIZE
        proteinuric_area = pd.to_numeric(proteinuric_df["Area"], errors="coerce").dropna()
        mean_area_proteinuric = (proteinuric_area.mean() if len(proteinuric_area) > 0 else np.nan)

        if (np.isfinite(mean_area_proteinuric) and mean_area_proteinuric > 0):
            log10_mean_area_proteinuric = (np.log10(mean_area_proteinuric))
        else:
            log10_mean_area_proteinuric = np.nan
        area_matrix[i, j] = log10_mean_area_proteinuric

        # BOOTSTRAP 95% CI
        if (len(proteinuric_values) > 0 and len(control_values) > 0):

            ci95_lower, ci95_upper = (bootstrap_log2_ratio_ci(control=control_values.to_numpy(),
                                                              proteinuric=proteinuric_values.to_numpy(),
                                                              n_bootstrap=N_BOOTSTRAP, confidence_level=0.95, rng=rng))
        else:
            ci95_lower = np.nan
            ci95_upper = np.nan

        # TWO-SIDED MANN-WHITNEY U TEST
        if (len(proteinuric_values) > 1 and len(control_values) > 1):

            try:
                u_stat, p_value = mannwhitneyu(proteinuric_values, control_values, alternative="two-sided")
                u_stat = float(u_stat)
                p_value = float(p_value)

            except ValueError:
                u_stat = np.nan
                p_value = np.nan

        else:
            u_stat = np.nan
            p_value = np.nan

        p_matrix[i, j] = p_value

        # Nominal significance retained for reporting
        significant_nominal = (bool(p_value < ALPHA) if np.isfinite(p_value) else False)

        # SUPPLEMENTARY TABLE
        supp_rows.append({ "Isotopologue": peptide, "Organ": organ, "Ratio_column": ratio_col,
                           "Control_condition": "wtRQ_Heavy", "Proteinuric_condition": "AVRQ_Heavy",
                           "n_Control_mice": n_control_mice, "n_Proteinuric_mice": n_proteinuric_mice,
                           "mean_isotope_to_light_ratio_Control":
                               (float(mean_ratio_control) if np.isfinite(mean_ratio_control) else np.nan),
                           "mean_isotope_to_light_ratio_Proteinuric":
                               (float(mean_ratio_proteinuric) if np.isfinite(mean_ratio_proteinuric) else np.nan),
                           "ratio_Proteinuric_vs_Control":
                               (float(ratio_proteinuric_vs_control) if np.isfinite(ratio_proteinuric_vs_control) else np.nan),
                           "log2_ratio_Proteinuric_vs_Control":
                               (float(log2_effect) if np.isfinite(log2_effect) else np.nan),
                           "log2_ratio_CI95_lower": ci95_lower, "log2_ratio_CI95_upper": ci95_upper,
                           "mean_area_Proteinuric":
                               (float(mean_area_proteinuric) if np.isfinite(mean_area_proteinuric) else np.nan),
                           "log10_mean_area_Proteinuric":
                               (float(log10_mean_area_proteinuric) if np.isfinite(log10_mean_area_proteinuric) else np.nan),
                           "U_statistic": u_stat, "P_value_nominal": p_value, "nominal_P_lt_0.05": significant_nominal})


# CREATE SUPPLEMENTARY DATAFRAME
supp = pd.DataFrame(supp_rows)

# MULTIPLE-TESTING CORRECTIONS
# Three significance definitions are retained and plotted:
# 1) Nominal P < 0.05
# 2) Global BH-FDR < 0.20
#    BH correction is performed across ALL available isotopologue x organ comparisons.
# 3) Within-metabolite BH-FDR < 0.05
#    A metabolite is defined by the part of the isotopologue identifier before "_13C" (for example Arg_13C10 -> Arg).
#    BH correction is performed separately for each metabolite across all of its available isotopologue x organ tests.

# Parent metabolite identifier used for within-metabolite correction.
supp["Metabolite"] = (supp["Isotopologue"].astype(str).str.split("_13C", n=1).str[0])

# 1) NOMINAL P < 0.05
supp["Nominal_P_lt_0.05"] = (supp["P_value_nominal"].lt(ALPHA) & np.isfinite(supp["P_value_nominal"]))

# 2) GLOBAL BH-FDR < 0.20
supp["P_value_BH_FDR_global"] = bh_fdr(supp["P_value_nominal"].to_numpy(dtype=float))
supp["Global_FDR_lt_0.20"] = (supp["P_value_BH_FDR_global"].lt(GLOBAL_FDR_CUTOFF) &
                              np.isfinite(supp["P_value_BH_FDR_global"]))

# 3) WITHIN-METABOLITE BH-FDR < 0.05
supp["P_value_BH_FDR_within_metabolite"] = np.nan

for metabolite, metabolite_index in supp.groupby("Metabolite", sort=False).groups.items():

    metabolite_index = list(metabolite_index)
    supp.loc[metabolite_index, "P_value_BH_FDR_within_metabolite"] = \
        bh_fdr(supp.loc[metabolite_index, "P_value_nominal"].to_numpy(dtype=float))

supp["Within_metabolite_FDR_lt_0.05"] = (supp["P_value_BH_FDR_within_metabolite"].lt(WITHIN_METABOLITE_FDR_CUTOFF) &
                                         np.isfinite(supp["P_value_BH_FDR_within_metabolite"]))

# SIGNIFICANCE MATRICES FOR THE THREE PLOTS
def significance_matrix(value_column):
    return (
        supp.pivot(index="Isotopologue", columns="Organ", values=value_column)
            .reindex(index=peptides_to_include, columns=organs).fillna(False).astype(bool).to_numpy())

sig_matrix_nominal = significance_matrix("Nominal_P_lt_0.05")
sig_matrix_global_fdr020 = significance_matrix("Global_FDR_lt_0.20")
sig_matrix_within_metabolite_fdr005 = significance_matrix("Within_metabolite_FDR_lt_0.05")

# EFFECT SIZE + 95% CI
supp["log2_ratio_Proteinuric_vs_Control_95CI"] = [format_effect_ci(effect, lower, upper)
    for effect, lower, upper
    in zip(supp["log2_ratio_Proteinuric_vs_Control"], supp["log2_ratio_CI95_lower"], supp["log2_ratio_CI95_upper"])]

# SORT RESULTS
supp["Isotopologue"] = pd.Categorical(supp["Isotopologue"], categories=peptides_to_include, ordered=True)
supp["Organ"] = pd.Categorical(supp["Organ"], categories=organs, ordered=True)
supp = (supp.sort_values(["Isotopologue", "Organ"]).reset_index(drop=True))

# SUPPLEMENTARY TABLE
supp_r = supp[["Isotopologue", "Metabolite", "Organ", "n_Control_mice", "n_Proteinuric_mice",
                      "mean_isotope_to_light_ratio_Control", "mean_isotope_to_light_ratio_Proteinuric",
                      "ratio_Proteinuric_vs_Control", "log2_ratio_Proteinuric_vs_Control_95CI", "U_statistic",
                      "P_value_nominal", "Nominal_P_lt_0.05", "P_value_BH_FDR_global", "Global_FDR_lt_0.20",
                      "P_value_BH_FDR_within_metabolite", "Within_metabolite_FDR_lt_0.05"]].copy()

supp_r = supp_r.rename(
    columns={"n_Control_mice": "n Control mice",
             "n_Proteinuric_mice": "n Proteinuric mice",
             "mean_isotope_to_light_ratio_Control": "Mean isotope/light ratio Control",
             "mean_isotope_to_light_ratio_Proteinuric": "Mean isotope/light ratio Proteinuric",
             "ratio_Proteinuric_vs_Control": "Proteinuric/Control ratio",
             "log2_ratio_Proteinuric_vs_Control_95CI": "log2 ratio Proteinuric/Control (95% CI)",
             "U_statistic": "Mann-Whitney U",
             "P_value_nominal": "Nominal P value",
             "Nominal_P_lt_0.05": "Nominal P < 0.05",
             "P_value_BH_FDR_global": "Global BH-adjusted P value",
             "Global_FDR_lt_0.20": "Global BH-FDR < 0.20",
             "P_value_BH_FDR_within_metabolite": "Within-metabolite BH-adjusted P value",
             "Within_metabolite_FDR_lt_0.05": "Within-metabolite BH-FDR < 0.05"})

# ROUND NUMERICAL COLUMNS
columns_to_round_3 = ["Mean isotope/light ratio Control", "Mean isotope/light ratio Proteinuric",
                      "Proteinuric/Control ratio", "Mann-Whitney U"]

for col_name in columns_to_round_3:
    supp_r[col_name] = pd.to_numeric(supp_r[col_name], errors="coerce").round(3)

for col_name in ["Nominal P value", "Global BH-adjusted P value", "Within-metabolite BH-adjusted P value"]:
    supp_r[col_name] = pd.to_numeric(supp_r[col_name], errors="coerce").round(6)

# PRINT BIOLOGICAL REPLICATE SUMMARY
print("\nBiological replicate numbers")
print("----------------------------")
print("Control mice, minimum n:", supp["n_Control_mice"].min())
print("Control mice, maximum n:", supp["n_Control_mice"].max())
print("Proteinuric mice, minimum n:", supp["n_Proteinuric_mice"].min())
print("Proteinuric mice, maximum n:", supp["n_Proteinuric_mice"].max())

# PRINT STATISTICAL SUMMARY
n_tests = int(supp["P_value_nominal"].notna().sum())
n_nominal = int(supp["Nominal_P_lt_0.05"].sum())
n_global_fdr020 = int(supp["Global_FDR_lt_0.20"].sum())
n_within_metabolite_fdr005 = int(supp["Within_metabolite_FDR_lt_0.05"].sum())

print("\nStatistical summary")
print("-------------------")
print("Number of valid comparisons:", n_tests)
print("Nominal P < 0.05:", n_nominal)
print("Global BH-FDR < 0.20:", n_global_fdr020)
print("Within-metabolite BH-FDR < 0.05:", n_within_metabolite_fdr005)

# DOT SIZE:
# Dot size represents log10 mean peak area of the corresponding labeled isotopologue in proteinuric mice.

valid_areas = area_matrix[np.isfinite(area_matrix) & (area_matrix > 0)]

if len(valid_areas) > 0:
    max_area = np.max(valid_areas)
else:
    max_area = np.nan

if (np.isfinite(max_area) and max_area > 0):
    # A circle cannot have a negative radius. Values <= 0 have no circle.
    positive_areas = np.where(np.isfinite(area_matrix) & (area_matrix > 0), area_matrix, 0.0)
    R = positive_areas / max_area / 2.3
else:
    R = np.zeros_like(area_matrix)

# PLOT FUNCTION
def make_dot_plot(significance_matrix, significance_label, output_name):
    """Create one isotopologue dot plot.

    The only difference between the three exported versions is
    the significance criterion used for the brown cell shading."""

    fig, ax = plt.subplots(figsize=(10, 9.8))
    ax.set_aspect("equal", adjustable="box")
    cs = ["#022c53", "#044685", "white", "#a90704", "#770502"]
    cmap = (colors.LinearSegmentedColormap.from_list("cmap_name", cs))

    # Significance shading
    for i in range(N):
        for j in range(M):
            if significance_matrix[i, j]:

                rect = mpatches.Rectangle((j - 0.5, i - 0.5), 1, 1, facecolor="#b19a6a", alpha=0.4, edgecolor="black",
                                          linewidth=1, zorder=0.5,)
                ax.add_patch(rect)

    # Circles
    x, y = np.meshgrid(np.arange(M), np.arange(N))
    circles = [plt.Circle((j, i), radius=r, zorder=2) for r, j, i in zip(R.flat, x.flat, y.flat)]
    plot_effects = np.ma.masked_invalid(effect_matrix)
    col = PatchCollection(circles, array=plot_effects.flatten(), cmap=cmap, norm=colors.CenteredNorm(halfrange=3))
    ax.add_collection(col)

    # Grid and labels
    ax.set_xticks(np.arange(M), labels=organs)
    ax.set_yticks(np.arange(N), labels=peptides_to_include)
    ax.set_xticks(np.arange(M + 1) - 0.5, minor=True)
    ax.set_yticks(np.arange(N + 1) - 0.5, minor=True)
    ax.tick_params(which="minor", size=0)
    ax.grid(which="minor", linestyle="-", linewidth=0.5, color="gray")
    ax.set_xlim(-0.5, M - 0.5)
    ax.set_ylim(N - 0.5, -0.5)

    # Colorbar
    cbar = fig.colorbar(col, ax=ax, shrink=0.3)
    cbar.set_label("log2(mean isotope/light ratio Proteinuric / Control)")

    # Dot-size legend: reference values use the same radius scale as plot dots.
    reference_values = [float(np.max(valid_areas)) if len(valid_areas) else np.nan,
                        float(np.median(valid_areas)) if len(valid_areas) else np.nan,
                        float(np.min(valid_areas)) if len(valid_areas) else np.nan]
    reference_radii = [value / max_area / 2.3 if np.isfinite(max_area) and max_area > 0 else 0.0
                       for value in reference_values]
    legend_handles = [Line2D([], [], marker="o", linestyle="None", color="gray", markersize=1,
                             label=f"{value:.1f}") for value in reference_values]
    legend2 = ax.legend(handles=legend_handles, loc="upper left", bbox_to_anchor=(1.02, 0.8),
                        title="log10 mean area", labelspacing=1.5)
    ax.add_artist(legend2)

    # Significance legend
    legend_patches = [mpatches.Patch(color="#b19a6a", alpha=0.4, label=significance_label)]
    legend1 = ax.legend(handles=legend_patches, loc="upper left", bbox_to_anchor=(1.02, 1), title="Mann-Whitney U")
    ax.add_artist(legend1)
    plt.tight_layout(pad=5.0)

    # Match legend diameters to actual plot-dot diameters after final layout.
    fig.canvas.draw()
    pixels_per_data_unit = abs(ax.transData.transform((1, 0))[0] - ax.transData.transform((0, 0))[0])
    for handle, radius in zip(legend2.get_lines(), reference_radii):
        handle.set_markersize(2 * radius * pixels_per_data_unit * 72 / fig.dpi)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    output_stem = (OUTPUT_DIR / output_name)

    for extension in ["svg", "pdf", "png"]:
        fig.savefig(f"{output_stem}.{extension}", dpi=300, transparent=True)

    plt.close(fig)

# EXPORT THREE SIGNIFICANCE VERSIONS
make_dot_plot(significance_matrix=sig_matrix_nominal, significance_label="P < 0.05", output_name="nominal_P005")
make_dot_plot(significance_matrix=sig_matrix_global_fdr020, significance_label="BH-FDR < 0.20",
              output_name="global_FDR020")
make_dot_plot(significance_matrix=sig_matrix_within_metabolite_fdr005,
              significance_label="BH-FDR < 0.05", output_name="within_metabolite_FDR005")

# SUPPLEMENTARY MATRIX OUTPUTS
mat_ratio = (supp.pivot(index="Isotopologue", columns="Organ", values="ratio_Proteinuric_vs_Control")
             .reindex(index=peptides_to_include, columns=organs))
mat_effect = (supp.pivot(index="Isotopologue", columns="Organ", values="log2_ratio_Proteinuric_vs_Control")
              .reindex(index=peptides_to_include, columns=organs))
mat_ci_lower = (supp.pivot(index="Isotopologue", columns="Organ", values="log2_ratio_CI95_lower")
                .reindex(index=peptides_to_include, columns=organs))
mat_ci_upper = (supp.pivot(index="Isotopologue", columns="Organ", values="log2_ratio_CI95_upper")
                .reindex(index=peptides_to_include, columns=organs))
mat_p = (supp.pivot(index="Isotopologue", columns="Organ", values="P_value_nominal")
         .reindex(index=peptides_to_include, columns=organs))
mat_global_fdr = (supp.pivot(index="Isotopologue", columns="Organ", values="P_value_BH_FDR_global")
                  .reindex(index=peptides_to_include, columns=organs))
mat_within_metabolite_fdr = (supp.pivot(index="Isotopologue", columns="Organ", values="P_value_BH_FDR_within_metabolite")
                             .reindex(index=peptides_to_include, columns=organs))
mat_nominal_sig = (supp.pivot(index="Isotopologue", columns="Organ", values="Nominal_P_lt_0.05")
                   .reindex(index=peptides_to_include, columns=organs))
mat_global_fdr_sig = (supp.pivot(index="Isotopologue", columns="Organ", values="Global_FDR_lt_0.20")
                      .reindex(index=peptides_to_include, columns=organs))
mat_within_metabolite_fdr_sig = (supp.pivot(index="Isotopologue", columns="Organ", values="Within_metabolite_FDR_lt_0.05")
                                 .reindex(index=peptides_to_include, columns=organs))

# WRITE SUPPLEMENTARY EXCEL FILE
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

with pd.ExcelWriter(OUT_XLSX, engine="openpyxl") as writer:

    supp_r.to_excel(writer, sheet_name="Results_long", index=False)
    mat_ratio.to_excel(writer, sheet_name="Mean_ratio")
    mat_effect.to_excel(writer, sheet_name="Log2_effect")
    mat_ci_lower.to_excel(writer, sheet_name="CI95_lower")
    mat_ci_upper.to_excel(writer, sheet_name="CI95_upper")
    mat_p.to_excel(writer, sheet_name="Nominal_P")
    mat_global_fdr.to_excel(writer, sheet_name="Global_BH_FDR")
    mat_within_metabolite_fdr.to_excel(writer, sheet_name="WithinMetab_BH_FDR")
    mat_nominal_sig.to_excel(writer, sheet_name="Nominal_P_lt_005")
    mat_global_fdr_sig.to_excel(writer, sheet_name="Global_FDR_lt_020")
    mat_within_metabolite_fdr_sig.to_excel(writer, sheet_name="WithinMet_FDR_lt_005")

print("\nWrote supplementary table:", OUT_XLSX)
print("\nAll outputs saved in:", OUTPUT_DIR.resolve())