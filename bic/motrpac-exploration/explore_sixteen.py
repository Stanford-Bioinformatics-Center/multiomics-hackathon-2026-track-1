"""Reproducible significance audit of the 16 mapped plasma candidates.

Uses existing consortium estimates; does not fit participant-level models.
Outputs a Markdown report, exhaustive JSON results, and a scientific figure.
"""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.patches import Patch

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "sixteen_candidates"
OUT.mkdir(exist_ok=True)
all_rows = pd.read_csv(ROOT / "table3/all_candidate_results.csv")
original_plasma = pd.read_csv(ROOT / "table3/plasma_results.csv")
plasma = original_plasma[original_plasma.identity_resolved].copy()
candidate_ids = list(plasma.candidate_id.unique())
data = all_rows[all_rows.candidate_id.isin(candidate_ids)].copy()
direct = data[(data.interpretation_scope == "circulating_abundance") &
              (data.contrast_type == "Endur_vs_Resist")].copy()
context = data[(data.interpretation_scope != "circulating_abundance") &
               (data.contrast_type == "exercise_with_controls")].copy()

def bh(values):
    p = np.asarray(values, dtype=float)
    order = np.argsort(p, kind="stable")
    ranked = p[order] * len(p) / np.arange(1, len(p) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1].clip(0, 1)
    result = np.empty(len(p))
    result[order] = ranked
    return result

assert len(candidate_ids) == 16
assert len(plasma) == 230 and len(direct) == 92
assert not data.duplicated(["candidate_id", "source_object", "feature_id", "contrast_type", "contrast_category", "Timepoint"]).any()
assert data[["logFC", "p_value", "adj_p_value", "CI.L_calculated", "CI.R_calculated"]].notna().all().all()
assert data.p_value.between(0, 1).all() and data.adj_p_value.between(0, 1).all()
assert (data["CI.L_calculated"] <= data.logFC).all() and (data.logFC <= data["CI.R_calculated"]).all()
plasma["candidate_family_bh"] = bh(plasma.p_value)
assert np.allclose(plasma.candidate_family_bh, plasma.table3_screen_bh, rtol=1e-12, atol=1e-15)
assert set(plasma.loc[plasma.adj_p_value < .05, "candidate_name"]) == {"Fractalkine", "Lactate"}
assert set(plasma.loc[plasma.candidate_family_bh < .05, "candidate_name"]) == {"Fractalkine", "Lactate", "FGF21", "IL-15"}
assert set(direct.loc[direct.adj_p_value < .05, "candidate_name"]) == {"Lactate"}

times = ["during_20_min", "during_40_min", "post_10_min", "post_15_30_45_min", "post_3.5_4_hr", "post_24_hr"]
short = {"Angiopoietin-like protein 4": "ANGPTL4", "Angiopoietin 1": "Angiopoietin-1", "Myostatin (GDF8)": "Myostatin", "Fractalkine": "Fractalkine / CX3CL1", "VEGF": "VEGF / VEGFA", "HSP72": "HSP72 / HSPA1A"}

def time_label(row):
    t, tissue = row["Timepoint"], row["tissue"]
    mapping = {"during_20_min": "during 20 min", "during_40_min": "during 40 min", "post_10_min": "post 10 min", "post_24_hr": "post 24 h"}
    if t == "post_15_30_45_min":
        return "post " + {"blood": "30", "muscle": "15", "adipose": "45"}[tissue] + " min"
    if t == "post_3.5_4_hr":
        return "post " + ("4" if tissue == "adipose" else "3.5") + " h"
    return mapping[t]

def fmt(x):
    return f"{x:.3g}"

def direction(x):
    return "increase" if x > 0 else "decrease" if x < 0 else "zero"

def describe(row):
    return f'{row["contrast_category"]}, {time_label(row)}'

def table(headers, rows):
    return "\n".join(["| " + " | ".join(headers) + " |", "|" + "|".join(["---"] * len(headers)) + "|"] +
                     ["| " + " | ".join(map(str, row)) + " |" for row in rows])

best = plasma.loc[plasma.groupby("candidate_id").adj_p_value.idxmin()].sort_values("adj_p_value")
summaries = []
for _, b in best.iterrows():
    q = plasma[plasma.candidate_id == b.candidate_id]
    mode = direct[direct.candidate_id == b.candidate_id]
    tissue = context[context.candidate_id == b.candidate_id]
    summaries.append({
        "candidate_id": b.candidate_id, "candidate_name": b.candidate_name,
        "n_plasma_features": len(q[["source_object", "feature_id"]].drop_duplicates()),
        "n_plasma_tests": len(q), "n_nominal_hits": int((q.p_value < .05).sum()),
        "n_source_fdr_hits": int((q.adj_p_value < .05).sum()),
        "n_candidate_family_bh_hits": int((q.candidate_family_bh < .05).sum()),
        "best_source_fdr_row": b.to_dict(),
        "minimum_candidate_family_bh": float(q.candidate_family_bh.min()),
        "minimum_ee_source_fdr": float(q.loc[q.contrast_category == "EE-CON", "adj_p_value"].min()),
        "minimum_re_source_fdr": float(q.loc[q.contrast_category == "RE-CON", "adj_p_value"].min()),
        "minimum_direct_mode_source_fdr": float(mode.adj_p_value.min()),
        "n_direct_mode_source_fdr_hits": int((mode.adj_p_value < .05).sum()),
        "n_context_source_fdr_hits": int((tissue.adj_p_value < .05).sum()),
    })

def records(frame):
    return json.loads(frame.to_json(orient="records", double_precision=15))

payload = {
    "collection": "human-precovid-sed-adu c2.0",
    "package_version": "2.0.8",
    "commit": "535b4044e7417413de471104c619120337602b77",
    "analysis_date": "2026-09-26",
    "primary_threshold": "Original MoTrPAC adj_p_value < 0.05; per tissue/assay/platform/contrast",
    "secondary_family": "Exploratory BH across all 230 resolved plasma exercise-versus-control tests; excludes the 92 direct mode tests",
    "summaries": summaries, "plasma_exercise_control": records(plasma),
    "plasma_direct_mode": records(direct), "tissue_context": records(context),
    "all_candidate_rows": records(data),
    "input_sha256": {str(f.relative_to(ROOT)): hashlib.sha256(f.read_bytes()).hexdigest() for f in
                     [ROOT / "table3/all_candidate_results.csv", ROOT / "table3/plasma_results.csv"]}
}
# Replace NaN in optional gene annotations with JSON null.
payload = json.loads(pd.io.json.dumps(payload)) if hasattr(pd.io.json, "dumps") else payload
def clean(v):
    if isinstance(v, dict): return {k: clean(x) for k, x in v.items()}
    if isinstance(v, list): return [clean(x) for x in v]
    if isinstance(v, (float, np.floating)) and not np.isfinite(v): return None
    return v
(OUT / "results.json").write_text(json.dumps(clean(payload), indent=2, allow_nan=False))

report = ["# Significance of 16 Table 3 candidates in human acute exercise", "",
          "Analysis date: 2026-09-26. MoTrPAC human pre-suspension sedentary adult collection c2.0; public Analysis package 2.0.8, commit `535b4044e7417413de471104c619120337602b77`.", "",
          "## Findings", "",
          "- **Fractalkine/CX3CL1 and lactate** have circulating responses passing the original MoTrPAC FDR threshold.",
          "- **FGF21 and IL-15** have decreasing plasma responses passing a separate exploratory BH correction over these 230 candidate tests, but not the original source correction.",
          "- **Lactate alone** has a source-FDR-significant direct endurance-versus-resistance plasma contrast among these 16 candidates.",
          "- Several candidates with nonsignificant plasma results have tissue RNA responses. These are separate observations, not evidence that the tissue released the protein into circulation.", "",
          "## What was tested", "",
          "Fifteen protein candidates comprise 21 Olink assay IDs: IL-6 and IL-8 each have four IDs, retained separately. Lactate has a research metabolomics measurement and a clinical chemistry measurement. The resulting 23 assay features each have six endurance-versus-control and four resistance-versus-control tests: 230 tests. The same 23 features each have four direct endurance-versus-resistance tests: 92 tests. There are no during-resistance samples.", "",
          "All effects are consortium difference-in-changes estimates: change from pre-exercise in an exercise arm minus the time-matched change in controls. EE = endurance; RE = resistance; CON = non-exercise control. We did not refit individual-level models.", "",
          "## Strongest source-adjusted plasma result for every candidate", "",
          "Each row selects the smallest original adjusted p-value across that candidate's tested assay/time/mode conditions. This is a descriptive minimum, not a candidate-level omnibus p-value. Directions for nonsignificant results describe estimates only. Effect and confidence interval use the native source model scale; do not compare effect magnitudes across the clinical lactate and other assays.", ""]
rows = []
for b in best.to_dict("records"):
    status = "Source FDR < 0.05" if b["adj_p_value"] < .05 else "Secondary family only" if b["candidate_family_bh"] < .05 else "No source-FDR hit"
    rows.append([b["candidate_name"], describe(b), f'{b["logFC"]:+.3f} [{b["CI.L_calculated"]:.3f}, {b["CI.R_calculated"]:.3f}]', fmt(b["p_value"]), fmt(b["adj_p_value"]), fmt(b["candidate_family_bh"]), status])
report += [table(["Candidate", "Condition", "Effect [pointwise 95% CI]", "Raw p", "Source adjusted p", "Candidate-family BH", "Interpretation"], rows), "",
           "The lactate row above is the clinical assay. The research metabolomics assay also passes source FDR, with its smallest adjusted p = " + fmt(plasma.loc[(plasma.candidate_name == "Lactate") & (plasma.platform == "metab-u-ionpneg"), "adj_p_value"].min()) + ". Agreement of two assays in this cohort is not independent-cohort replication.", "",
           "## Every source-FDR-significant plasma result", ""]
hits = plasma[plasma.adj_p_value < .05].sort_values(["candidate_name", "platform", "contrast_category", "Timepoint"])
report += [table(["Candidate", "Platform", "Condition", "Effect", "Source adjusted p"],
                 [[r.candidate_name, r.platform, describe(r), f"{r.logFC:+.4f}", fmt(r.adj_p_value)] for _, r in hits.iterrows()]), "",
           "## Exercise-mode comparisons", "",
           "Only direct EE-RE contrasts support a difference between exercise protocols. Significance in one exercise-versus-control arm alone is insufficient. A negative EE-RE effect means the endurance change is lower than the resistance change. Timing refers to time after the end of each protocol, not equal elapsed time after starting exercise.", ""]
dh = direct[direct.adj_p_value < .05]
report += [table(["Candidate", "Platform", "Time", "EE-RE effect", "Source adjusted p"],
                 [[r.candidate_name, r.platform, time_label(r), f"{r.logFC:+.4f}", fmt(r.adj_p_value)] for _, r in dh.iterrows()]), "",
           "CX3CL1 has no significant direct post-exercise EE-RE contrast (minimum source adjusted p = " + fmt(direct.loc[direct.candidate_name == "Fractalkine", "adj_p_value"].min()) + "). Its during-endurance response cannot be compared with during-resistance, which was not sampled.", "",
           "## Tissue and blood-cell context", "",
           "Shown below is the strongest source-FDR result per candidate/tissue/assay layer with at least one hit. The full JSON retains all significant and nonsignificant rows. Blood RNA is a cellular expression measurement, distinct from plasma protein. No new cross-tissue omnibus correction was applied. Phosphosite effects do not directly establish protein activity.", ""]
sig_context = context[context.adj_p_value < .05]
cb = sig_context.loc[sig_context.groupby(["candidate_id", "tissue", "assay"]).adj_p_value.idxmin()]
report += [table(["Candidate", "Tissue", "Layer / matched gene", "Condition", "Effect", "Source adjusted p"],
                 [[r.candidate_name, r.tissue, r.assay + " / " + (str(r.matched_gene) if pd.notna(r.matched_gene) else r.feature_id), describe(r), f"{r.logFC:+.3f}", fmt(r.adj_p_value)] for _, r in cb.iterrows()]), "",
           "HSP72 plasma is mapped to HSPA1A. Tissue context also includes HSPA1B as an explicitly labeled related gene; its RNA response is not a measurement of plasma HSPA1A.", "",
           "## Multiple testing and uncertainty", "",
           "The primary threshold is the original `adj_p_value < 0.05`, which uses BH within tissue, assay, platform, and contrast. It does not control an additional omnibus family spanning every timepoint and assay in this report. The secondary candidate-family calculation applies BH to all 230 raw plasma exercise-versus-control p-values. This is the same family used in the earlier Table 3 screen: the 12 unmatched/ambiguous candidates supplied no included tests. It is exploratory because these data have already been examined. FGF21 and IL-15 should be described as secondary findings requiring confirmation, not upgraded to primary discoveries.", "",
           "A 95% CI excluding zero or a raw p below 0.05 can coexist with a nonsignificant adjusted p because the displayed CI is pointwise, not multiplicity-adjusted. Selecting a candidate's smallest p-value can exaggerate its apparent strength. All assay IDs, timepoints, and nonsignificant results are retained in `results.json` to make that selection visible. A nonsignificant result is not an equivalence test. Neither plasma abundance nor tissue RNA establishes secretion, disease benefit, or causal mediation.", "",
           "## Interpretation for follow-up", "",
           "1. **CX3CL1** is the clearest protein candidate for a disease-evidence overlay: circulating endurance-associated increases and strong muscle RNA responses, with source attribution and disease direction still unresolved.",
           "2. **Lactate** provides a strong acute metabolic response and a directly supported difference between exercise protocols. Its known abundance response does not by itself establish a novel disease mechanism.",
           "3. **FGF21 and IL-15** merit targeted validation of the observed post-resistance decreases, retaining their exploratory status and specific sampling times.",
           "4. **MSTN, ANGPTL4, VEGFA, and SDC4** are candidates for tissue-response hypotheses; their RNA responses should not be presented as significant circulating exerkine responses.", "",
           "## Reproducibility", "",
           "The source-object verification independently checks extracted effects, confidence intervals, raw p-values, and adjusted p-values against the pinned R data objects. The Python analysis checks unique keys, valid p-value ranges, confidence-interval bounds, candidate counts, contrast coverage, and reproduces the 230-test BH correction. No restricted participant data were used.", "",
           "Sources: [public MoTrPAC package](https://github.com/MoTrPAC/MotrpacHumanPreSuspensionAnalysis/tree/535b4044e7417413de471104c619120337602b77), [contrast and FDR definitions](https://github.com/MoTrPAC/MotrpacHumanPreSuspensionAnalysis/blob/535b4044e7417413de471104c619120337602b77/R/load_differential_analysis.R), [portal collection](https://motrpac-data.org/data-download/file-browser/analysis/human-precovid-sed-adu/c2.0).", ""]
(OUT / "REPORT.md").write_text("\n".join(report))

# A categorical significance map retains every physical assay and sampled contrast.
# The sign marks the model estimate; the color marks significance, not effect size.
feature_rows = plasma[["candidate_id", "candidate_name", "source_object", "platform", "feature_id"]].drop_duplicates().sort_values(["candidate_id", "platform", "feature_id"])
columns = [("EE-CON", t) for t in times] + [("RE-CON", t) for t in times[2:]]
matrix = np.full((len(feature_rows), len(columns)), np.nan)
signs = np.empty(matrix.shape, dtype=object)
labels = []
for y, f in enumerate(feature_rows.to_dict("records")):
    q = plasma[(plasma.candidate_id == f["candidate_id"]) & (plasma.source_object == f["source_object"]) & (plasma.feature_id == f["feature_id"])]
    label = short.get(f["candidate_name"], f["candidate_name"])
    label += " / clinical" if f["platform"] == "metab-t-clinical" else " / metabolomics" if f["candidate_name"] == "Lactate" else " / " + f["feature_id"]
    labels.append(label)
    for x, (arm, t) in enumerate(columns):
        row = q[(q.contrast_category == arm) & (q.Timepoint == t)]
        assert len(row) == 1
        r = row.iloc[0]
        matrix[y, x] = (2 if r.logFC > 0 else 1) if r.adj_p_value < .05 else 3 if r.candidate_family_bh < .05 else 0
        signs[y, x] = "+" if r.logFC > 0 else "−" if r.logFC < 0 else "0"

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10})
fig, ax = plt.subplots(figsize=(13, 11))
colors = ["#edf0f3", "#2564a4", "#b64b35", "#f3cb75"]
ax.imshow(matrix, cmap=ListedColormap(colors), norm=BoundaryNorm(np.arange(-.5, 4.5), 4), aspect="auto")
for y in range(matrix.shape[0]):
    for x in range(matrix.shape[1]):
        ax.text(x, y, signs[y, x], ha="center", va="center", color="white" if matrix[y, x] in (1, 2) else "#465464", fontsize=11)
ax.set_yticks(np.arange(len(labels)), labels, fontsize=9.5)
time_labels = ["During\n20 min", "During\n40 min", "Post\n10 min", "Post\n30 min", "Post\n3.5 h", "Post\n24 h"]
ax.set_xticks(range(10), time_labels + time_labels[2:], fontsize=9)
ax.xaxis.tick_top()
ax.tick_params(axis="both", length=0, pad=7)
ax.set_xticks(np.arange(-.5, 10), minor=True)
ax.set_yticks(np.arange(-.5, len(labels)), minor=True)
ax.grid(which="minor", color="white", linewidth=1.5)
ax.tick_params(which="minor", length=0)
ax.axvline(5.5, color="#667788", linewidth=2)
ax.text(2.5, -2.1, "ENDURANCE − CONTROL", ha="center", fontweight="bold", color="#30475c")
ax.text(7.5, -2.1, "RESISTANCE − CONTROL", ha="center", fontweight="bold", color="#30475c")
for edge in ax.spines.values(): edge.set_visible(False)
fig.suptitle("Plasma significance across all 16 candidates", x=.025, y=.985, ha="left", fontsize=17, fontweight="bold")
fig.text(.025, .952, "MoTrPAC human acute exercise · c2.0 · 23 assay features, 230 exercise-versus-control tests", fontsize=11, color="#465464")
legend = [Patch(facecolor=colors[2], label="Increase: source FDR < 0.05"), Patch(facecolor=colors[1], label="Decrease: source FDR < 0.05"), Patch(facecolor=colors[3], label="Secondary candidate-family BH only"), Patch(facecolor=colors[0], label="Neither threshold passed")]
fig.legend(handles=legend, loc="lower left", bbox_to_anchor=(.02, .045), ncol=2, frameon=False, fontsize=10)
fig.text(.025, .025, "+ / − shows estimate direction, including nonsignificant estimates. Orange findings are exploratory.\nIL-6 and IL-8 retain four assay IDs each. No during-resistance samples. Colors do not encode effect size.", fontsize=9, color="#465464")
fig.subplots_adjust(left=.36, right=.985, top=.84, bottom=.15)
fig.savefig(OUT / "plasma_significance.png", dpi=180, facecolor="white")
fig.savefig(OUT / "plasma_significance.svg", facecolor="white")
plt.close(fig)

print(json.dumps({"candidates": len(candidate_ids), "all_verified_rows_expected": len(data), "plasma_tests": len(plasma), "direct_mode_tests": len(direct), "context_tests": len(context), "source_fdr_hits": len(hits), "secondary_bh_hits": int((plasma.candidate_family_bh < .05).sum()), "direct_mode_source_hits": len(dh), "context_source_hits": len(sig_context)}, indent=2))
print(best[["candidate_name", "contrast_category", "Timepoint", "adj_p_value", "candidate_family_bh"]].to_string(index=False))
print("Wrote", OUT)
