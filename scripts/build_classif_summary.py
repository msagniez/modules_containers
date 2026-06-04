#!/usr/bin/env python3
"""
Rebuilds MPXXXX_Classification-summary.csv from classifier subdirectories.

Expected directory layout (relative to --classifiers-dir):
  ALLcatchR/results/ALLcatchR_lineage_results.csv
  ALLcatchR/results/ALLcatchR_ball_results.csv   (used for B-ALL and Unclassified samples)
  ALLcatchR/results/ALLcatchR_tall_results.csv   (used for T-ALL samples)
  MD-ALL/results/MD-ALL_predictions.csv
  MnM/results/MnM_lineage-predictions.csv
  MnM/results/MnM_subtype-predictions.csv
  SIGNATURE/results/SIGNATURE_B_predictions.csv
  SIGNATURE/results/SIGNATURE_T_predictions.csv
  Maylis/results/Maylis_results.csv
  AMLmapR/results/AMLmapR_predictions.csv
  AttentionAML/results/AttentionAML_predictions.csv

Output columns: Sample, classifier, subtype, score
"""

import argparse
import os
import pandas as pd


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def normalise_sample_name(name: str) -> str:
    """Strip leading 'X' and replace dots used instead of dashes (R convention)."""
    name = str(name).strip()
    if name.startswith("X"):
        name = name[1:]
    # R replaces '-' with '.' in row names for some files
    name = name.replace(".", "-")
    return name


# ---------------------------------------------------------------------------
# Per-classifier parsers
# ---------------------------------------------------------------------------

def parse_allcatchr_lineage(results_dir: str) -> pd.DataFrame:
    """
    ALLcatchR_lineage_results.csv
    Columns: Sample, sample, B-ALL, T-ALL, prediction
    Score = max(B-ALL, T-ALL)
    """
    path = os.path.join(results_dir, "ALLcatchR_lineage_results.csv")
    df = pd.read_csv(path)
    rows = []
    for _, r in df.iterrows():
        sample = str(r["Sample"]).strip()
        subtype = str(r["prediction"]).strip()
        score = max(float(r["B-ALL"]), float(r["T-ALL"]))
        # Mirror the summary: score is the winning class probability
        score = float(r["B-ALL"]) if subtype == "B-ALL" else float(r["T-ALL"])
        # Unclassified: use the higher of the two
        if subtype == "Unclassified":
            score = max(float(r["B-ALL"]), float(r["T-ALL"]))
        rows.append({"Sample": sample, "classifier": "ALLcatchR_lineage",
                     "subtype": subtype, "score": round(score, 9)})
    return pd.DataFrame(rows)


def parse_allcatchr_subtype(results_dir: str) -> pd.DataFrame:
    """
    Routes each sample to the correct subtype file based on ALLcatchR lineage:
      - B-ALL (or Unclassified) → ALLcatchR_ball_results.csv
        Columns of interest: Sample, Score, Prediction
      - T-ALL                   → ALLcatchR_tall_results.csv
        Columns of interest: Sample, BC_pred  (the main-cluster prediction string,
        e.g. "C8 (BCL11B)"), plus the per-cluster probability columns to derive
        the winning score.
    """
    # --- Load lineage predictions to build a per-sample routing map -----------
    lineage_path = os.path.join(results_dir, "ALLcatchR_lineage_results.csv")
    lineage_df = pd.read_csv(lineage_path)
    lineage_map = {
        str(r["Sample"]).strip(): str(r["prediction"]).strip()
        for _, r in lineage_df.iterrows()
    }

    rows = []

    # --- B-ALL subtype (also used for Unclassified) ---------------------------
    ball_path = os.path.join(results_dir, "ALLcatchR_ball_results.csv")
    ball_df = pd.read_csv(ball_path)
    for _, r in ball_df.iterrows():
        sample = str(r["Sample"]).strip()
        lineage = lineage_map.get(sample, "Unclassified")
        if lineage == "T-ALL":
            continue  # will be handled by the T-ALL block below
        rows.append({"Sample": sample, "classifier": "ALLcatchR_subtype",
                     "subtype": str(r["Prediction"]).strip(),
                     "score": round(float(r["Score"]), 9)})

    # --- T-ALL subtype --------------------------------------------------------
    # BC_pred contains the winning cluster label (or empty string if none).
    # The per-cluster probability columns are named like "C1 (TAL1 ...)", "C2 ...", etc.
    # We derive the score as the probability of the winning cluster.
    tall_path = os.path.join(results_dir, "ALLcatchR_tall_results.csv")
    tall_df = pd.read_csv(tall_path)

    # Identify the numeric cluster-probability columns (all columns between
    # "C1 (...)" and the last "C..." column, i.e. everything before the
    # non-numeric annotation columns that follow).
    cluster_cols = [c for c in tall_df.columns
                    if c.startswith("C") and c not in ("BC_pred",)]

    for _, r in tall_df.iterrows():
        sample = str(r["Sample"]).strip()
        lineage = lineage_map.get(sample, "Unclassified")
        if lineage != "T-ALL":
            continue  # only process samples called T-ALL by the lineage step

        bc_pred = str(r.get("BC_pred", "")).strip()
        subtype = bc_pred if bc_pred else "Unclassified"

        # Score = probability of the predicted cluster, or 0 if not found
        score = 0.0
        if bc_pred and cluster_cols:
            # BC_pred may contain multiple entries separated by ";"
            # e.g. "C13 (NUP214, MLLT10, KMT2A);C8 (BCL11B)"
            # Take the first one as the primary call
            primary = bc_pred.split(";")[0].strip()
            # Find the matching column (exact match)
            if primary in tall_df.columns:
                val = r[primary]
                score = float(val) if pd.notna(val) else 0.0
            else:
                # Fall back to the highest probability across all cluster cols
                vals = pd.to_numeric(r[cluster_cols], errors="coerce")
                score = float(vals.max()) if not vals.isna().all() else 0.0

        rows.append({"Sample": sample, "classifier": "ALLcatchR_subtype",
                     "subtype": subtype,
                     "score": round(score, 9)})

    return pd.DataFrame(rows)


def parse_md_all(results_dir: str) -> pd.DataFrame:
    """
    MD-ALL_predictions.csv
    Columns: id, PhenoGraph_pred, PhenoGraph_predScore, PhenoGraph_predLabel,
             svm_pred, svm_predScore, svm_predLabel
    Emits two rows per sample: MD-ALL_phenograph and MD-ALL_svm.
    """
    path = os.path.join(results_dir, "MD-ALL_predictions.csv")
    df = pd.read_csv(path)
    rows = []
    for _, r in df.iterrows():
        sample = str(r["id"]).strip()
        rows.append({"Sample": sample, "classifier": "MD-ALL_phenograph",
                     "subtype": str(r["PhenoGraph_pred"]).strip(),
                     "score": round(float(r["PhenoGraph_predScore"]), 9)})
        rows.append({"Sample": sample, "classifier": "MD-ALL_svm",
                     "subtype": str(r["svm_pred"]).strip(),
                     "score": round(float(r["svm_predScore"]), 9)})
    return pd.DataFrame(rows)


def parse_mnm_lineage(results_dir: str) -> pd.DataFrame:
    """
    MnM_lineage-predictions.csv
    Columns: (row name), predict, probability1, predict2, probability2, ...
    """
    path = os.path.join(results_dir, "MnM_lineage-predictions.csv")
    df = pd.read_csv(path, index_col=0)
    rows = []
    for idx, r in df.iterrows():
        sample = normalise_sample_name(idx)
        rows.append({"Sample": sample, "classifier": "MnM_lineage",
                     "subtype": str(r["predict"]).strip(),
                     "score": round(float(r["probability1"]), 9)})
    return pd.DataFrame(rows)


def parse_mnm_subtype(results_dir: str) -> pd.DataFrame:
    """
    MnM_subtype-predictions.csv — same layout as lineage file.
    """
    path = os.path.join(results_dir, "MnM_subtype-predictions.csv")
    df = pd.read_csv(path, index_col=0)
    rows = []
    for idx, r in df.iterrows():
        sample = normalise_sample_name(idx)
        rows.append({"Sample": sample, "classifier": "MnM_subtype",
                     "subtype": str(r["predict"]).strip(),
                     "score": round(float(r["probability1"]), 9)})
    return pd.DataFrame(rows)


def parse_signature(results_dir: str) -> pd.DataFrame:
    """
    SIGNATURE_B_predictions.csv  and  SIGNATURE_T_predictions.csv
    Format: rows = signature names, columns = sample names.
    For each sample, pick the signature with the highest score.
    B and T files cover different (possibly overlapping) samples; merge both.
    """
    rows = []
    for fname in ("SIGNATURE_B_predictions.csv", "SIGNATURE_T_predictions.csv"):
        path = os.path.join(results_dir, fname)
        if not os.path.exists(path):
            continue
        df = pd.read_csv(path, index_col=0)
        for col in df.columns:
            sample = normalise_sample_name(col)
            best_idx = df[col].astype(float).idxmax()
            best_score = float(df[col][best_idx])
            rows.append({"Sample": sample, "classifier": "SIGNATURE",
                         "subtype": str(best_idx).strip(),
                         "score": round(best_score, 9)})

    # If a sample appeared in both B and T files, keep the one with the higher score
    result = pd.DataFrame(rows)
    if result.empty:
        return result
    result = (result
              .sort_values("score", ascending=False)
              .drop_duplicates(subset="Sample")
              .sort_index())
    return result.reset_index(drop=True)


def parse_maylis(results_dir: str) -> pd.DataFrame:
    """
    Maylis_results.csv
    Columns: Sample, Subgroup_def, Subgroup_fr, Probability   (Probability is 0-100)
    """
    path = os.path.join(results_dir, "Maylis_results.csv")
    df = pd.read_csv(path)
    rows = []
    for _, r in df.iterrows():
        sample = str(r["Sample"]).strip()
        subtype = str(r["Subgroup_fr"]).strip()
        score = round(float(r["Probability"]) / 100.0, 9)
        rows.append({"Sample": sample, "classifier": "Maylis",
                     "subtype": subtype, "score": score})
    return pd.DataFrame(rows)


def parse_amlmapr(results_dir: str) -> pd.DataFrame:
    """
    AMLmapR_predictions.csv
    Columns: AML.MRC.1., ..., prediction, pass_cutoff, sample_id
    Score = the numeric value of the predicted subtype column.
    Subtype label is taken from the 'prediction' column as-is.
    """
    path = os.path.join(results_dir, "AMLmapR_predictions.csv")
    df = pd.read_csv(path)
    # Numeric score columns are everything except the last three meta columns
    score_cols = [c for c in df.columns if c not in ("prediction", "pass_cutoff", "sample_id")]
    rows = []
    for _, r in df.iterrows():
        sample = str(r["sample_id"]).strip()
        subtype = str(r["prediction"]).strip()
        # The score column name uses dots instead of special chars (R convention)
        # Try to find the matching score column by exact name first, then fuzzy
        score = float("nan")
        if subtype in score_cols:
            score = float(r[subtype])
        else:
            # Map subtype back to column name: replace non-alphanumeric with '.'
            candidate = re.sub(r"[^A-Za-z0-9]", ".", subtype)
            # Try exact and trailing-dot variants
            for col in score_cols:
                if col == candidate or col == candidate + ".":
                    score = float(r[col])
                    break
        rows.append({"Sample": sample, "classifier": "AMLmapR",
                     "subtype": subtype, "score": round(score, 9)})
    return pd.DataFrame(rows)


def parse_attentionaml(results_dir: str) -> pd.DataFrame:
    """
    AttentionAML_predictions.csv
    Columns: Samples, Subtype_pred, Probabolity  (typo in header is intentional)
    """
    path = os.path.join(results_dir, "AttentionAML_predictions.csv")
    df = pd.read_csv(path)
    rows = []
    for _, r in df.iterrows():
        sample = str(r["Samples"]).strip()
        subtype = str(r["Subtype_pred"]).strip()
        # Handle the typo gracefully — check both spellings
        if "Probabolity" in df.columns:
            score = float(r["Probabolity"])
        else:
            score = float(r["Probability"])
        rows.append({"Sample": sample, "classifier": "AttentionAML",
                     "subtype": subtype, "score": round(score, 9)})
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

PARSERS = [
    ("ALLcatchR/results", parse_allcatchr_lineage),
    ("ALLcatchR/results", parse_allcatchr_subtype),
    ("MD-ALL/results",    parse_md_all),
    ("MnM/results",       parse_mnm_lineage),
    ("MnM/results",       parse_mnm_subtype),
    ("SIGNATURE/results", parse_signature),
    ("Maylis/results",    parse_maylis),
    ("AMLmapR",     "results", parse_amlmapr),
    ("AttentionAML","results", parse_attentionaml),
]


def build_summary(classifiers_dir: str) -> pd.DataFrame:
    all_parts = []
    for classifier_dir, results_subdir, parser_fn in PARSERS:
        classifier_path = os.path.join(classifiers_dir, classifier_dir)
        results_dir     = os.path.join(classifier_path, results_subdir)

        # Skip entire classifier if its top-level directory is absent
        if not os.path.isdir(classifier_path):
            print(f"  [SKIP] classifier not found: {classifier_dir}/")
            continue

        # Skip if the results sub-directory is missing (ran but produced nothing)
        if not os.path.isdir(results_dir):
            print(f"  [SKIP] results dir missing:  {classifier_dir}/{results_subdir}/")
            continue

        try:
            part = parser_fn(results_dir)
            print(f"  [OK]   {parser_fn.__name__:35s}  → {len(part)} rows")
            all_parts.append(part)
        except Exception as exc:
            print(f"  [ERR]  {parser_fn.__name__}: {exc}")

    if not all_parts:
        raise RuntimeError("No data parsed — check your classifiers directory.")

    summary = pd.concat(all_parts, ignore_index=True)
    summary["score"] = summary["score"].round(9)
    return summary


def main():
    parser = argparse.ArgumentParser(
        description="Build MPXXXX_Classification-summary.csv from classifier subdirectories."
    )
    parser.add_argument(
        "--classifiers-dir", "-d",
        default=".",
        help="Path to the classifiers directory (default: current directory)"
    )
    parser.add_argument(
        "--output", "-o",
        default="MPXXXX_Classification-summary.csv",
        help="Output CSV path (default: MPXXXX_Classification-summary.csv)"
    )
    args = parser.parse_args()

    classifiers_dir = os.path.abspath(args.classifiers_dir)
    print(f"Classifiers directory : {classifiers_dir}")
    print(f"Output file           : {args.output}\n")

    summary = build_summary(classifiers_dir)

    summary.to_csv(args.output, index=False)
    print(f"\nWrote {len(summary)} rows to {args.output}")


if __name__ == "__main__":
    main()
