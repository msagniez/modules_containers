#!/usr/bin/env python3
"""
Builds MPXXXX_Classification-summary.csv from classifier subdirectories.

Expected directory layout (relative to --classifiers-dir):
  ALLcatchR/results/ALLcatchR_lineage_results.csv
  ALLcatchR/results/ALLcatchR_ball_results.csv   (used for B-ALL and Unclassified samples)
  ALLcatchR/results/ALLcatchR_tall_results.csv   (used for T-ALL samples)
  MD-ALL/results/MD-ALL_predictions.tsv
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
import csv as _csv
import os
import re
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


# ALLcatchR filename conventions
#   new: ALLcatchR_lineage_results.csv / ALLcatchR_ball_results.csv / ALLcatchR_tall_results.csv
#   old: gene_counts_lineage-prediction.tsv / gene_counts_ball-prediction.tsv / gene_counts_tall-prediction.tsv
_MD_ALL_NAMES = [
    "MD-ALL_predictions.csv",
    "MD-ALL_predictions.tsv",
]

def _md_all_path(results_dir: str) -> tuple[str, str]:
    """Return (path, sep) for the MD-ALL predictions file."""
    for fname in _MD_ALL_NAMES:
        path = os.path.join(results_dir, fname)
        if os.path.exists(path):
            sep = "\t" if fname.endswith(".tsv") else ","
            return path, sep
    checked = ", ".join(_MD_ALL_NAMES)
    raise FileNotFoundError(
        f"MD-ALL predictions file not found in {results_dir}. Tried: {checked}"
    )


_ALLCATCHR_NAMES = {
    "lineage": [
        "ALLcatchR_lineage_results.csv",
        "gene_counts_lineage-prediction.tsv",
    ],
    "ball": [
        "ALLcatchR_ball_results.csv",
        "gene_counts_ball-prediction.tsv",
    ],
    "tall": [
        "ALLcatchR_tall_results.csv",
        "gene_counts_tall-prediction.tsv",
    ],
}

def _allcatchr_path(results_dir: str, key: str) -> tuple[str, str]:
    """
    Return (path, sep) for the requested ALLcatchR file (key: lineage/ball/tall).
    Tries each known filename in order and returns the first one that exists.
    Raises FileNotFoundError listing all candidates if none is found.
    """
    candidates = _ALLCATCHR_NAMES[key]
    for fname in candidates:
        path = os.path.join(results_dir, fname)
        if os.path.exists(path):
            sep = "\t" if fname.endswith(".tsv") else ","
            return path, sep
    checked = ", ".join(candidates)
    raise FileNotFoundError(
        f"ALLcatchR {key} file not found in {results_dir}. Tried: {checked}"
    )


# ---------------------------------------------------------------------------
# Per-classifier parsers
# ---------------------------------------------------------------------------

def parse_allcatchr_lineage(results_dir: str) -> pd.DataFrame:
    """
    ALLcatchR lineage file — two naming conventions supported:
      - ALLcatchR_lineage_results.csv  (newer runs)
      - gene_counts_lineage-prediction.tsv  (older runs)
    Columns: Sample, sample, B-ALL, T-ALL, prediction
    Score = winning class probability; for Unclassified, the higher of the two.
    """
    path, sep = _allcatchr_path(results_dir, "lineage")
    df = pd.read_csv(path, sep=sep)
    df.columns = [c.strip() for c in df.columns]
    df = df.rename(columns={"sample": "Sample"})  # normalise case
    rows = []
    for _, r in df.iterrows():
        sample  = str(r["Sample"]).strip()
        subtype = str(r["prediction"]).strip()
        if subtype == "B-ALL":
            score = float(r["B-ALL"])
        elif subtype == "T-ALL":
            score = float(r["T-ALL"])
        else:
            score = max(float(r["B-ALL"]), float(r["T-ALL"]))
        rows.append({"Sample": sample, "classifier": "ALLcatchR_lineage",
                     "subtype": subtype, "score": round(score, 9)})
    return pd.DataFrame(rows)


def parse_allcatchr_subtype(results_dir: str) -> pd.DataFrame:
    """
    Routes each sample to the correct subtype file based on ALLcatchR lineage.
    Two naming conventions are supported for each file:
      - ALLcatchR_ball_results.csv / gene_counts_ball-prediction.tsv  (B-ALL)
      - ALLcatchR_tall_results.csv / gene_counts_tall-prediction.tsv  (T-ALL)
    """
    lineage_path, lineage_sep = _allcatchr_path(results_dir, "lineage")
    lineage_df   = pd.read_csv(lineage_path, sep=lineage_sep)
    lineage_df = lineage_df.rename(columns={"sample": "Sample"})  # normalise case
    lineage_map  = {
        str(r["Sample"]).strip(): str(r["prediction"]).strip()
        for _, r in lineage_df.iterrows()
    }

    rows = []

    # B-ALL / Unclassified
    ball_path, ball_sep = _allcatchr_path(results_dir, "ball")
    ball_df   = pd.read_csv(ball_path, sep=ball_sep)
    ball_df.columns = [c.strip() for c in ball_df.columns]
    ball_df = ball_df.rename(columns={"sample": "Sample", "prediction": "Prediction", "score": "Score"})
    for _, r in ball_df.iterrows():
        sample  = str(r["Sample"]).strip()
        lineage = lineage_map.get(sample, "Unclassified")
        if lineage == "T-ALL":
            continue
        rows.append({"Sample": sample, "classifier": "ALLcatchR_subtype",
                     "subtype": str(r["Prediction"]).strip(),
                     "score":   round(float(r["Score"]), 9)})

    # T-ALL
    tall_path, tall_sep = _allcatchr_path(results_dir, "tall")
    tall_df = pd.read_csv(tall_path, sep=tall_sep)
    tall_df.columns = [c.strip() for c in tall_df.columns]
    tall_df = tall_df.rename(columns={"sample": "Sample"})

    # Cluster probability columns: positions 2–23 (1-indexed) = indices 1–22
    score_cols = list(tall_df.columns[1:23])

    for _, r in tall_df.iterrows():
        sample  = str(r["Sample"]).strip()
        lineage = lineage_map.get(sample, "Unclassified")
        if lineage != "T-ALL":
            continue

        vals    = pd.to_numeric(r[score_cols], errors="coerce")
        subtype = vals.idxmax() if not vals.isna().all() else "Unclassified"
        score   = float(vals.max()) if not vals.isna().all() else 0.0

        # Strip cluster number, keep only the description in parentheses
        # e.g. "C13 (NUP214, MLLT10, KMT2A)" -> "NUP214, MLLT10, KMT2A"
        subtype = re.sub(r"^C\d+[\.\d]*\s*\((.+)\)$", r"\1", subtype)

        rows.append({"Sample": sample, "classifier": "ALLcatchR_subtype",
                     "subtype": subtype, "score": round(score, 9)})

    return pd.DataFrame(rows)


def parse_md_all(results_dir: str) -> pd.DataFrame:
    """
    MD-ALL_predictions.csv
    Nominal columns: id, PhenoGraph_pred, PhenoGraph_predScore, PhenoGraph_predLabel,
                     svm_pred, svm_predScore, svm_predLabel

    The *predLabel columns contain embedded commas
    (e.g. "PAX5alt,9(0.82)|PAX5::ETV6,2(0.18)") which breaks standard CSV
    parsing and shifts all subsequent column indices.

    Strategy: read line-by-line with csv.reader, anchor on fixed left-side
    positions (id=0, phenograph_pred=1, phenograph_score=2) and scan
    right-to-left for the first parseable float to locate svm_predScore,
    with svm_pred one field to its left.
    """
    path, sep = _md_all_path(results_dir)
    # TSV files don't have the embedded-comma problem — use pandas directly
    if sep == "\t":
        df = pd.read_csv(path, sep="\t")
        rows = []
        for _, r in df.iterrows():
            sample = str(r["id"]).strip()
            rows.append({"Sample": sample, "classifier": "MD-ALL_phenograph",
                         "subtype": str(r["PhenoGraph_pred"]).strip(),
                         "score":   round(float(r["PhenoGraph_predScore"]), 9)})
            rows.append({"Sample": sample, "classifier": "MD-ALL_svm",
                         "subtype": str(r["svm_pred"]).strip(),
                         "score":   round(float(r["svm_predScore"]), 9)})
        return pd.DataFrame(rows)
    # CSV: label columns contain embedded commas — parse line by line
    rows = []
    with open(path, newline="") as fh:
        reader = _csv.reader(fh)
        next(reader)  # skip header
        for fields in reader:
            if len(fields) < 7:
                continue
            sample           = fields[0].strip()
            phenograph_pred  = fields[1].strip()
            phenograph_score = float(fields[2])

            # Scan right-to-left for svm_predScore
            svm_pred  = None
            svm_score = None
            for i in range(len(fields) - 1, 3, -1):
                try:
                    svm_score = float(fields[i])
                    svm_pred  = fields[i - 1].strip()
                    break
                except ValueError:
                    continue

            rows.append({"Sample": sample, "classifier": "MD-ALL_phenograph",
                         "subtype": phenograph_pred,
                         "score":   round(phenograph_score, 9)})
            if svm_pred is not None:
                rows.append({"Sample": sample, "classifier": "MD-ALL_svm",
                             "subtype": svm_pred,
                             "score":   round(svm_score, 9)})
    return pd.DataFrame(rows)


def parse_mnm_lineage(results_dir: str) -> pd.DataFrame:
    """
    MnM_lineage-predictions.csv
    Columns: (row name), predict, probability1, predict2, probability2, ...
    """
    path = os.path.join(results_dir, "MnM_lineage-predictions.csv")
    df   = pd.read_csv(path, index_col=0)
    rows = []
    for idx, r in df.iterrows():
        sample = normalise_sample_name(idx)
        rows.append({"Sample": sample, "classifier": "MnM_lineage",
                     "subtype": str(r["predict"]).strip(),
                     "score":   round(float(r["probability1"]), 9)})
    return pd.DataFrame(rows)


def parse_mnm_subtype(results_dir: str) -> pd.DataFrame:
    """
    MnM_subtype-predictions.csv -- same layout as lineage file.
    """
    path = os.path.join(results_dir, "MnM_subtype-predictions.csv")
    df   = pd.read_csv(path, index_col=0)
    rows = []
    for idx, r in df.iterrows():
        sample = normalise_sample_name(idx)
        rows.append({"Sample": sample, "classifier": "MnM_subtype",
                     "subtype": str(r["predict"]).strip(),
                     "score":   round(float(r["probability1"]), 9)})
    return pd.DataFrame(rows)


def parse_signature(results_dir: str) -> pd.DataFrame:
    """
    SIGNATURE_B_predictions.csv and SIGNATURE_T_predictions.csv
    Format: rows = signature names, columns = sample names.
    For each sample, pick the signature with the highest score.
    If a sample appears in both files, keep the higher-scoring entry.
    """
    rows = []
    for fname in ("SIGNATURE_B_predictions.csv", "SIGNATURE_T_predictions.csv"):
        path = os.path.join(results_dir, fname)
        if not os.path.exists(path):
            continue
        df = pd.read_csv(path, index_col=0)
        for col in df.columns:
            sample     = normalise_sample_name(col)
            col_vals = df[col].astype(float)
            if col_vals.isna().all():
                print(f"  [WARN] SIGNATURE: all-NaN scores for sample "
                      f"'{sample}' in {fname}; skipping")
                continue
            best_idx   = col_vals.idxmax(skipna=True)
            best_score = float(df[col][best_idx])
            rows.append({"Sample": sample, "classifier": "SIGNATURE",
                         "subtype": str(best_idx).strip(),
                         "score":   round(best_score, 9)})

    result = pd.DataFrame(rows)
    if result.empty:
        return result
    result = (result
              .sort_values("score", ascending=False)
              .drop_duplicates(subset="Sample")
              .sort_index()
              .reset_index(drop=True))
    return result


def parse_amlmapr(results_dir: str) -> pd.DataFrame:
    """
    AMLmapR_predictions.csv
    Columns: AML.MRC.1., ..., prediction, pass_cutoff, sample_id
    Score = the numeric value of the predicted subtype column.
    """
    path       = os.path.join(results_dir, "AMLmapR_predictions.csv")
    df         = pd.read_csv(path)
    score_cols = [c for c in df.columns
                  if c not in ("prediction", "pass_cutoff", "sample_id")]
    rows = []
    for _, r in df.iterrows():
        sample  = str(r["sample_id"]).strip()
        subtype = str(r["prediction"]).strip()
        score   = float("nan")
        if subtype in score_cols:
            score = float(r[subtype])
        else:
            candidate = re.sub(r"[^A-Za-z0-9]", ".", subtype)
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
    Columns: Samples, Subtype_pred, Probabolity  (typo in source header)
    """
    path = os.path.join(results_dir, "AttentionAML_predictions.csv")
    df   = pd.read_csv(path)
    rows = []
    for _, r in df.iterrows():
        sample  = str(r["Samples"]).strip()
        subtype = str(r["Subtype_pred"]).strip()
        score   = float(r["Probabolity"] if "Probabolity" in df.columns
                        else r["Probability"])
        rows.append({"Sample": sample, "classifier": "AttentionAML",
                     "subtype": subtype, "score": round(score, 9)})
    return pd.DataFrame(rows)


def parse_maylis(results_dir: str) -> pd.DataFrame:
    """
    Maylis_results.csv
    Columns: Sample, Subgroup_def, Subgroup_fr, Probability  (0-100 scale)
    """
    path = os.path.join(results_dir, "Maylis_results.csv")
    df   = pd.read_csv(path)
    rows = []
    for _, r in df.iterrows():
        sample  = str(r["Sample"]).strip()
        subtype = str(r["Subgroup_fr"]).strip()
        score   = round(float(r["Probability"]) / 100.0, 9)
        rows.append({"Sample": sample, "classifier": "Maylis",
                     "subtype": subtype, "score": score})
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Parser registry  --  (classifier_subdir, results_subdir, parser_function)
# Entries are skipped silently when the classifier directory is absent,
# so the script works regardless of which classifiers were actually run.
# ---------------------------------------------------------------------------

PARSERS = [
    ("ALLcatchR",    "results", parse_allcatchr_lineage),
    ("ALLcatchR",    "results", parse_allcatchr_subtype),
    ("MD-ALL",       "results", parse_md_all),
    ("MnM",          "results", parse_mnm_lineage),
    ("MnM",          "results", parse_mnm_subtype),
    ("SIGNATURE",    "results", parse_signature),
    ("Maylis",       "results", parse_maylis),
    ("AMLmapR",      "results", parse_amlmapr),
    ("AttentionAML", "results", parse_attentionaml),
]


def build_summary(classifiers_dir: str) -> pd.DataFrame:
    all_parts = []
    for classifier_dir, results_subdir, parser_fn in PARSERS:
        classifier_path = os.path.join(classifiers_dir, classifier_dir)
        results_dir     = os.path.join(classifier_path, results_subdir)

        if not os.path.isdir(classifier_path):
            print(f"  [SKIP] classifier not found: {classifier_dir}/")
            continue

        if not os.path.isdir(results_dir):
            print(f"  [SKIP] results dir missing:  {classifier_dir}/{results_subdir}/")
            continue

        try:
            part = parser_fn(results_dir)
            print(f"  [OK]   {parser_fn.__name__:35s}  -> {len(part)} rows")
            all_parts.append(part)
        except Exception as exc:
            print(f"  [ERR]  {parser_fn.__name__}: {exc}")

    if not all_parts:
        raise RuntimeError("No data parsed -- check your classifiers directory.")

    summary = pd.concat(all_parts, ignore_index=True)
    summary["score"] = summary["score"].round(9)
    return summary


def main():
    parser = argparse.ArgumentParser(
        description="Rebuild Classification-summary.csv from classifier subdirectories."
    )
    parser.add_argument(
        "--classifiers-dir", "-d",
        default=".",
        help="Path to the classifiers directory (default: current directory)"
    )
    parser.add_argument(
        "--output", "-o",
        default="Classification-summary.csv",
        help="Output CSV path"
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