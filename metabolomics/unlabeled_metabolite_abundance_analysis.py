"""
Unlabeled metabolite abundance analysis for the 2-week arginine-diet study.

This script reads the deposited targeted-metabolomics CSV reports, filters to
unlabeled metabolites (IsotopeLabelType == 0), compares proteinuric
(AVRQ_Light) versus control (wtRQ_Light) mice across tissues, and exports:

1. a dot plot shaded by nominal P < 0.05;
2. a dot plot shaded by within-metabolite BH-FDR < 0.20;
3. an Excel workbook containing effect sizes, bootstrap 95% confidence
   intervals, nominal P values, global BH-FDR values (reference only), and
   within-metabolite BH-FDR values.

Expected repository layout:
    data/metabolomics/
    metabolomics/
    results/metabolomics/unlabeled_metabolites/

Run from the repository root:
    python metabolomics/unlabeled_metabolite_abundance_analysis.py
"""

import os
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib import colors
import matplotlib.patches as patches
from scipy.stats import mannwhitneyu

# SETTINGS
INPUT_DIR = Path("data/metabolomics")
OUTPUT_DIR = Path("results/metabolomics/unlabeled_metabolites")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
OUT_XLSX = ("SupplementaryTable_Fig2B_effectsize_CI_"
            "NominalP_and_WithinMetaboliteFDR020.xlsx")
OUT_PREFIX_NOMINAL = ("dotBlot_Heatmap_areasTotal_Diff_NominalP005")
OUT_PREFIX_FDR = ("dotBlot_Heatmap_areasTotal_Diff_WithinMetaboliteFDR020")

# Nominal P-value threshold retained for reporting
ALPHA = 0.05

# Within-metabolite BH-FDR threshold used for one figure version
FDR_CUTOFF = 0.20

# Bootstrap settings for 95% confidence intervals
N_BOOTSTRAP = 10000
RANDOM_SEED = 12345

organs = ["Ct", "Md", "Lv", "He", "Ms", "Br", "Lg", "Sp", "Dd", "Jj", "Il", "Cl", "Er", "Pl", "Ur"]
peptide_to_include = ["Arg", "Orn", "Urea", "Citr", "ArgSuc", "Glu", "Gln", "Asp", "Pro", "4OHPro", "Agm", "ADMA",
                      "MMA", "AcArg"]

csv_files = ["AVRQ_ArgDiet_2weeks_Ctx_Med_Liver_Heart_report.csv",
             "AVRQ_ArgDiet_2weeks_Lung_Pancreas_Spleen_Muscle_report.csv",
             "AVRQ_ArgDiet_2weeks_Brain_WAT_Duodenum_Ileum_report.csv",
             "AVRQ_ArgDiet_2weeks_Jejenum_Colon_Erythrocytes_report.csv"]

BASE_DIR = INPUT_DIR

print("Input directory used for files:", BASE_DIR)
print("Output directory:", OUTPUT_DIR)
print("\nChecking expected deposited CSV files:")

# RENAME DICTIONARIES
peptide_rename = {"ADMA": "ADMA", "Agmatine": "Agm", "Arginine": "Arg", "Argininosuccinic acid": "ArgSuc",
                  "Aspartic acid": "Asp", "Citrulline": "Citr", "GABA": "GABA", "Glutamic acid": "Glu",
                  "Glutamine": "Gln", "L-Homoarginine": "HArg", "MMA": "MMA", "N-a-Acetyl-L-arginine": "AcArg",
                  "Na-Acetyl-L-ornithine": "AcOrn", "N-Acetyl-L-glutamic acid": "AcGlu", "Ornithine": "Orn",
                  "Proline": "Pro", "Putrescine": "Putr", "Spermine": "Spn", "trans-4-Hydroxy-L-proline": "4OHPro",
                  "Urea": "Urea"}

organ_rename = {"Cortex": "Ct", "Medulla": "Md", "Liver": "Lv", "Heart": "He", "Muscle": "Ms", "Brain": "Br",
                "Lung": "Lg", "Spleen": "Sp", "Duodenum": "Dd", "Jejenum": "Jj", "Jejunum": "Jj", "Ileum": "Il",
                "Colon": "Cl", "Erythrocytes": "Er", "Plasma": "Pl", "Urine": "Ur"}

peptides_to_exclude = ["Guanidineacetic acid"]

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
    q = (ranked * n / np.arange(1, n + 1))
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.clip(q, 0, 1)
    inv = np.empty_like(order)
    inv[order] = np.arange(n)
    out[mask] = q[inv]

    return out

def bootstrap_log10_ratio_ci(control, proteinuric, n_bootstrap=10000, confidence_level=0.95, rng=None):
    """ Bootstrap confidence interval for:
        log10(mean abundance Proteinuric / mean abundance Control)
    Biological replicates are resampled independently within each experimental group.
    """

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
            effect = np.log10(prot_mean / ctrl_mean)
            bootstrap_effects.append(effect)

    if len(bootstrap_effects) == 0:
        return np.nan, np.nan

    bootstrap_effects = np.asarray(bootstrap_effects)
    alpha_ci = (1 - confidence_level)
    lower = np.percentile(bootstrap_effects, 100 * alpha_ci / 2)
    upper = np.percentile(bootstrap_effects, 100 * (1 - alpha_ci / 2 ))

    return (float(lower), float(upper))


def format_effect_ci(effect, lower, upper):
    """ Representation of effect size and 95% confidence interval. """

    if (np.isfinite(effect) and np.isfinite(lower) and np.isfinite(upper)):

        return (f"{effect:.2f} "
                f"({lower:.2f} to {upper:.2f})")
    return ""

def win_long_path(pathlike):
    p = os.fspath(pathlike)
    p = os.path.abspath(p)

    if os.name != "nt":
        return p
    if p.startswith("\\\\?\\"):
        return p
    if p.startswith("\\\\"):
        return ("\\\\?\\UNC\\" + p.lstrip("\\"))
    return ("\\\\?\\" + p)

# REPRODUCIBLE BOOTSTRAP
rng = np.random.default_rng(RANDOM_SEED)

# LOAD AND COMBINE SOURCE FILES
dfs = []

for file in csv_files:
    file_path = (BASE_DIR / file)
    file_path_long = win_long_path(file_path)

    print(f"  {file_path.name} -> "
          f"len={len(str(file_path))}")

    if not os.path.exists(file_path_long):
        print("\nCSV files visible "
              "in working directory:")

        for f in sorted(BASE_DIR.glob("*.csv")):
            print(" ", repr(f.name), "len=", len(str(BASE_DIR / f.name)))

        raise FileNotFoundError(f"Missing file: "
                                f"{file_path}")

    temp_df = pd.read_csv(file_path_long, sep=",")
    temp_df["SourceFile"] = file
    dfs.append(temp_df)

df = pd.concat(dfs, ignore_index=True)

# Remove exact duplicates
df = df.drop_duplicates()

# CLEAN / RENAME
df = df[~df["Peptide"].isin(peptides_to_exclude)].copy()
df["Peptide"] = (df["Peptide"].map(peptide_rename).fillna(df["Peptide"]))
df["Organ"] = (df["Organ"].map(organ_rename).fillna(df["Organ"]))
iso = pd.to_numeric(df["IsotopeLabelType"], errors="coerce")
df = df[iso.notna()].copy()
df["IsotopeLabelType"] = (iso[iso.notna()].astype(int))

# FILTER DATA FOR THIS FIGURE
df = df[df["Organ"].isin(organs)].copy()
df = df[df["Peptide"].isin(peptide_to_include)].copy()

# Unlabeled metabolites only
df = df[df["IsotopeLabelType"] == 0].copy()
df = df[df["Condition"].isin(["wtRQ_Light", "AVRQ_Light"])].copy()
sub_df_ctrl = df[df["Condition"] == "wtRQ_Light"].copy()
sub_df_dis = df[df["Condition"] == "AVRQ_Light"].copy()

# PREPARE MATRICES
N = len(peptide_to_include)
M = len(organs)
effect_matrix = np.full((N, M), np.nan, dtype=float)
p_matrix = np.full((N, M), np.nan, dtype=float)
u_matrix = np.full((N, M), np.nan, dtype=float)
supp_rows = []

# LOOP THROUGH METABOLITES × ORGANS

for i, peptide in enumerate(peptide_to_include):
    for j, organ in enumerate(organs):

        # Biological replicate measurements
        proteinuric = (sub_df_dis[(sub_df_dis["Peptide"] == peptide) & (sub_df_dis["Organ"] == organ)]["Area"].dropna())
        control = (sub_df_ctrl[(sub_df_ctrl["Peptide"] == peptide) & (sub_df_ctrl["Organ"] == organ)]["Area"].dropna())

        # Biological replicate numbers
        n_control_mice = int(len(control))
        n_proteinuric_mice = int(len(proteinuric))

        # Mean abundance
        mean_control = (control.mean() if len(control) > 0 else np.nan)
        mean_proteinuric = (proteinuric.mean() if len(proteinuric) > 0 else np.nan)

        # Proteinuric / Control mean ratio
        if (np.isfinite(mean_control) and np.isfinite(mean_proteinuric) and mean_control > 0):
            mean_ratio = (mean_proteinuric / mean_control)

        else:
            mean_ratio = np.nan

        # EFFECT SIZE REPRESENTED BY CIRCLE COLOR - log10(mean Proteinuric / mean Control)
        if (np.isfinite(mean_ratio) and mean_ratio > 0):
            log10_ratio = np.log10(mean_ratio)
        else:
            log10_ratio = np.nan

        # BOOTSTRAP 95% CI
        if (len(control) > 0 and len(proteinuric) > 0):
            (ci95_lower, ci95_upper) = bootstrap_log10_ratio_ci(control=control.to_numpy(),
                                                                proteinuric=proteinuric.to_numpy(),
                                                                n_bootstrap=N_BOOTSTRAP, confidence_level=0.95,
                                                                rng=rng)
        else:
            ci95_lower = np.nan
            ci95_upper = np.nan

        # TWO-SIDED MANN-WHITNEY U TEST
        if (len(control) > 0 and len(proteinuric) > 0):

            try:
                (u_stat, p_value) = mannwhitneyu(proteinuric, control, alternative="two-sided")
                u_stat = float(u_stat)
                p_value = float(p_value)

            except ValueError:
                u_stat = np.nan
                p_value = np.nan

        else:
            u_stat = np.nan
            p_value = np.nan

        # NOMINAL SIGNIFICANCE
        significant_nominal = (bool(p_value < ALPHA) if np.isfinite(p_value) else False)

        # STORE MATRICES
        effect_matrix[i, j] = log10_ratio
        p_matrix[i, j] = p_value
        u_matrix[i, j] = u_stat

        # SUPPLEMENTARY TABLE ROW
        supp_rows.append({"Metabolite": peptide, "Organ": organ, "Control_condition": "wtRQ_Light",
                          "Proteinuric_condition": "AVRQ_Light", "IsotopeLabelType": 0,
                          "n_Control_mice": n_control_mice, "n_Proteinuric_mice": n_proteinuric_mice,
                          "mean_abundance_Control": (float(mean_control) if np.isfinite(mean_control) else np.nan),
                          "mean_abundance_Proteinuric": (float(mean_proteinuric) if np.isfinite(mean_proteinuric)
                                                         else np.nan),
                          "mean_ratio_Proteinuric_vs_Control": (float(mean_ratio) if np.isfinite(mean_ratio)
                                                                else np.nan),
                          "log10_ratio_Proteinuric_vs_Control": (float(log10_ratio) if np.isfinite(log10_ratio)
                                                                 else np.nan),
                          "log10_ratio_CI95_lower": ci95_lower, "log10_ratio_CI95_upper": ci95_upper,
                          "U_statistic": u_stat, "P_value_nominal": p_value, "nominal_P_lt_0.05": significant_nominal})

# CREATE SUPPLEMENTARY DATAFRAME
supp = pd.DataFrame(supp_rows)

# BENJAMINI-HOCHBERG FDR CORRECTIONS

# Global correction retained for reference only.
supp["P_value_BH_FDR_global"] = bh_fdr(supp["P_value_nominal"].to_numpy(dtype=float))

# Within-metabolite correction:
# for each metabolite, BH correction is performed across organs.
supp["P_value_BH_FDR_within_metabolite"] = np.nan

for metabolite, indices in (supp.groupby("Metabolite", sort=False).groups.items()):

    pvals_this_metabolite = (supp.loc[indices, "P_value_nominal"].to_numpy(dtype=float))

    corrected = bh_fdr(pvals_this_metabolite)

    supp.loc[indices, "P_value_BH_FDR_within_metabolite"] = corrected

supp["Within_metabolite_FDR_lt_0.20"] = ((supp["P_value_BH_FDR_within_metabolite"] < FDR_CUTOFF) &
                                         np.isfinite(supp["P_value_BH_FDR_within_metabolite"]))

# SIGNIFICANCE MATRICES FOR THE TWO FIGURE VERSIONS
sig_matrix_nominal = (supp.pivot(index="Metabolite", columns="Organ", values="nominal_P_lt_0.05")
                      .reindex(index=peptide_to_include, columns=organs).fillna(False).astype(bool).to_numpy())

sig_matrix_fdr = (supp.pivot(index="Metabolite", columns="Organ", values="Within_metabolite_FDR_lt_0.20")
                  .reindex(index=peptide_to_include, columns=organs).fillna(False).astype(bool).to_numpy())

# EFFECT SIZE + 95% CI
supp["log10_ratio_Proteinuric_vs_Control_95CI"] = [format_effect_ci(effect, lower, upper)
                                                   for effect, lower, upper
                                                   in zip(supp["log10_ratio_Proteinuric_vs_Control"],
                                                          supp["log10_ratio_CI95_lower"],
                                                          supp["log10_ratio_CI95_upper"])]

# SORT RESULTS
supp["Metabolite"] = pd.Categorical(supp["Metabolite"], categories=peptide_to_include, ordered=True)
supp["Organ"] = pd.Categorical(supp["Organ"], categories=organs, ordered=True)
supp = (supp.sort_values(["Metabolite", "Organ"]).reset_index(drop=True))

# CLEAN SUPPLEMENTARY TABLE
supp_r = supp[["Metabolite", "Organ", "n_Control_mice", "n_Proteinuric_mice", "mean_abundance_Control",
                      "mean_abundance_Proteinuric", "mean_ratio_Proteinuric_vs_Control",
                      "log10_ratio_Proteinuric_vs_Control_95CI", "U_statistic", "P_value_nominal",
                      "P_value_BH_FDR_global", "P_value_BH_FDR_within_metabolite",
                      "Within_metabolite_FDR_lt_0.20"]].copy()

# RENAME COLUMNS
supp_r = supp_r.rename(columns={"n_Control_mice": "n Control mice",
                                              "n_Proteinuric_mice": "n Proteinuric mice",
                                              "mean_abundance_Control": "Mean abundance Control",
                                              "mean_abundance_Proteinuric": "Mean abundance Proteinuric",
                                              "mean_ratio_Proteinuric_vs_Control": "Proteinuric/Control ratio",
                                              "log10_ratio_Proteinuric_vs_Control_95CI":
                                                  "log10 ratio Proteinuric/Control (95% CI)",
                                              "U_statistic": "Mann-Whitney U",
                                              "P_value_nominal": "Nominal P value",
                                              "P_value_BH_FDR_global": "Global BH-adjusted P value",
                                              "P_value_BH_FDR_within_metabolite":
                                                  "Within-metabolite BH-adjusted P value",
                                              "Within_metabolite_FDR_lt_0.20": "Within-metabolite BH-FDR < 0.20"})

# ROUND NUMERICAL VALUES
supp_r["Proteinuric/Control ratio"] = pd.to_numeric(supp_r["Proteinuric/Control ratio"], errors="coerce")\
    .round(3)
supp_r["Mann-Whitney U"] = pd.to_numeric(supp_r["Mann-Whitney U"], errors="coerce").round(3)

# Peak areas may be large, so retain useful numerical precision
supp_r["Mean abundance Control"] = pd.to_numeric(supp_r["Mean abundance Control"], errors="coerce")
supp_r["Mean abundance Proteinuric"] = pd.to_numeric(supp_r["Mean abundance Proteinuric"], errors="coerce")

# Retain more precision for P values
supp_r["Nominal P value"] = pd.to_numeric(supp_r["Nominal P value"], errors="coerce").round(6)

for col_name in ["Global BH-adjusted P value", "Within-metabolite BH-adjusted P value"]:

    supp_r[col_name] = pd.to_numeric(supp_r[col_name], errors="coerce").round(6)

# PRINT BIOLOGICAL REPLICATE SUMMARY
print("\nBiological replicate numbers")
print("----------------------------")
print("Control mice, minimum n:", supp["n_Control_mice"].min())
print("Control mice, maximum n:", supp["n_Control_mice"].max())
print("Proteinuric mice, minimum n:", supp["n_Proteinuric_mice"].min())
print("Proteinuric mice, maximum n:", supp["n_Proteinuric_mice"].max())

# PRINT STATISTICAL SUMMARY
n_tests = (supp["P_value_nominal"].notna().sum())
n_nominal = int(supp["nominal_P_lt_0.05"].sum())
n_fdr = int(supp["Within_metabolite_FDR_lt_0.20"].sum())

print("\nStatistical summary")
print("-------------------")
print("Number of valid comparisons:", n_tests)
print("Nominal P < 0.05:", n_nominal)
print("Within-metabolite BH-FDR < 0.20:", n_fdr)

# PLOT
valid_values = effect_matrix[np.isfinite(effect_matrix)]

if len(valid_values) > 0:
    print("\nMinimum log10 ratio:", np.min(valid_values))
    print("Maximum log10 ratio:", np.max(valid_values))

plot_values = np.nan_to_num(effect_matrix, nan=0.0)
cs = ["#044685", "white", "#a90704"]
cmap = (colors.LinearSegmentedColormap.from_list("cmap_name", cs))
x, y = np.meshgrid(np.arange(M), np.arange(N))
R = (np.ones((N, M), dtype=float) / 4)

def make_plot(significance_matrix, significance_label, output_prefix):

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.set_aspect("equal", adjustable="box")

    # SIGNIFICANCE SHADING
    for i in range(N):
        for j in range(M):
            if significance_matrix[i, j]:
                rect = patches.Rectangle((j - 0.5, i - 0.5), 1, 1, linewidth=0, edgecolor="yellow", facecolor="#b19a6a",
                                         alpha=0.4, zorder=0)
                ax.add_patch(rect)

    # ADD CIRCLES
    circles = [plt.Circle((j, i), radius=r) for r, j, i in zip(R.flat, x.flat, y.flat)]
    col = PatchCollection(circles, array=plot_values.flatten(), cmap=cmap, norm=colors.CenteredNorm(halfrange=0.4))
    ax.add_collection(col)

    # LABELS / GRID
    ax.set(xticks=np.arange(M), yticks=np.arange(N), xticklabels=organs, yticklabels=peptide_to_include)
    ax.set_xticks(np.arange(M + 1) - 0.5, minor=True)
    ax.set_yticks(np.arange(N + 1) - 0.5, minor=True)
    ax.grid(which="minor")
    ax.set_xlim(-0.5, M - 0.5)
    ax.set_ylim(N - 0.5, -0.5)

    # COLORBAR
    cbar = fig.colorbar(col, ax=ax)
    cbar.set_label("log10(mean Proteinuric / mean Control)")

    # SIGNIFICANCE LEGEND
    sig_patch = patches.Patch(facecolor="#b19a6a", alpha=0.4, label=significance_label)
    ax.legend(handles=[sig_patch], loc="upper left", bbox_to_anchor=(1.02, 1), frameon=False)
    plt.tight_layout()

    for extension in ["svg", "pdf", "png"]:
        plt.savefig(win_long_path(OUTPUT_DIR / f"{output_prefix}.{extension}"), dpi=300, transparent=True)
    plt.close(fig)

make_plot(significance_matrix=sig_matrix_nominal, significance_label="Nominal P < 0.05",
          output_prefix=OUT_PREFIX_NOMINAL)
make_plot(significance_matrix=sig_matrix_fdr, significance_label="Within-metabolite BH-FDR < 0.20",
          output_prefix=OUT_PREFIX_FDR)

# MATRIX OUTPUTS
mat_effect = (supp.pivot(index="Metabolite", columns="Organ", values="log10_ratio_Proteinuric_vs_Control")
              .reindex(index=peptide_to_include, columns=organs))
mat_ratio = (supp.pivot(index="Metabolite", columns= "Organ", values="mean_ratio_Proteinuric_vs_Control")
             .reindex(index=peptide_to_include, columns=organs))
mat_ci_lower = (supp.pivot(index="Metabolite", columns="Organ", values="log10_ratio_CI95_lower")
                .reindex(index=peptide_to_include, columns=organs))
mat_ci_upper = (supp.pivot(index="Metabolite", columns="Organ", values="log10_ratio_CI95_upper")
                .reindex(index=peptide_to_include, columns=organs))
mat_p = (supp.pivot(index="Metabolite", columns="Organ", values="P_value_nominal")
         .reindex(index=peptide_to_include, columns=organs))
mat_global_fdr = (supp.pivot(index="Metabolite", columns="Organ", values="P_value_BH_FDR_global")
                  .reindex(index=peptide_to_include, columns=organs))
mat_within_fdr = (supp.pivot(index="Metabolite", columns="Organ", values="P_value_BH_FDR_within_metabolite")
                  .reindex(index=peptide_to_include, columns=organs))
mat_nominal_sig = (supp.pivot(index="Metabolite", columns="Organ", values="nominal_P_lt_0.05")
                   .reindex(index=peptide_to_include, columns=organs))
mat_within_fdr_sig = (supp.pivot(index="Metabolite", columns="Organ", values="Within_metabolite_FDR_lt_0.20")
                      .reindex(index=peptide_to_include, columns=organs))

# WRITE SUPPLEMENTARY EXCEL FILE
out_xlsx_path = win_long_path(OUTPUT_DIR / OUT_XLSX)
with pd.ExcelWriter(out_xlsx_path, engine="openpyxl") as writer:

    # Main table
    supp_r.to_excel(writer, sheet_name= "Results_long", index=False)

    # Numerical matrices
    mat_ratio.to_excel(writer, sheet_name="Mean_ratio")
    mat_effect.to_excel(writer, sheet_name="Log10_ratio")
    mat_ci_lower.to_excel(writer, sheet_name="CI95_lower")
    mat_ci_upper.to_excel(writer, sheet_name="CI95_upper")
    mat_p.to_excel(writer, sheet_name="Nominal_P")
    mat_global_fdr.to_excel(writer, sheet_name="Global_BH_FDR")
    mat_within_fdr.to_excel(writer, sheet_name="WithinMetabolite_FDR")
    mat_nominal_sig.to_excel(writer, sheet_name="Nominal_P_sig")
    mat_within_fdr_sig.to_excel(writer, sheet_name="WithinMetab_FDR020_sig")

print("\nWrote supplementary table:", OUTPUT_DIR / OUT_XLSX)
print("\nFinished. All outputs saved in:", OUTPUT_DIR)