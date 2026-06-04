#!/usr/bin/env Rscript
# plot_classification_summary.R
# Reads a Classification-summary CSV and produces a per-sample barplot PDF.
#
# Usage:
#   Rscript plot_classification_summary.R \
#     --input  /data/MPXXXX_Classification-summary.csv \
#     --output /data/MPXXXX_Classification-summary.pdf

suppressPackageStartupMessages(library(ggplot2))

# ---------------------------------------------------------------------------
# CLI arguments
# ---------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

parse_args <- function(args) {
  input  <- NULL
  output <- NULL
  i <- 1
  while (i <= length(args)) {
    if (args[i] %in% c("--input",  "-i")) { input  <- args[i + 1]; i <- i + 2 }
    else if (args[i] %in% c("--output", "-o")) { output <- args[i + 1]; i <- i + 2 }
    else i <- i + 1
  }
  list(input = input, output = output)
}

opt <- parse_args(args)
if (is.null(opt$input))  stop("--input <csv> is required")
if (is.null(opt$output)) {
  opt$output <- sub("\\.csv$", ".pdf", opt$input)
  message("No --output specified; writing to: ", opt$output)
}

# ---------------------------------------------------------------------------
# Colour maps
# ---------------------------------------------------------------------------
useColsSubtype <- c(
  # Others / Unclassified
  "Others"                   = "#999999",
  "Unclassified"             = "#999999",
  "Unassigned"               = "#999999",
  "diverse"                  = "#999999",
  # Control
  "BMF"                      = "#999999",
  "Not malignant bone marrow"= "#999999",
  "Not malignant lymph node" = "#999999",
  "Not malignant lung"       = "#999999",
  "Not malignant subcutis"   = "#999999",
  "Not malignant (Hemato)"   = "#999999",

  # B-ALL
  "ALL"                      = "#9966CC",
  "B-ALL"                    = "#9966CC",
  "NOS B-ALL"                = "#9966CC",
  "BCL2/MYC"                 = "#B284BE",
  "CDX2/UBTF"                = "#B284BE",
  "CEBP"                     = "#B284BE",
  "ZEB2/CEBP"                = "#B284BE",
  "BCR-ABL1_or_like"         = "#C9A0DC",
  "Ph-like"                  = "#C9A0DC",
  "Ph"                       = "#C9A0DC",
  "Ph-pos"                   = "#C9A0DC",
  "BCR::ABL1 (like) B-ALL"   = "#C9A0DC",
  "DUX4-rearranged"          = "#db97ca",
  "DUX4"                     = "#db97ca",
  "DUX4 rearranged B-ALL"    = "#db97ca",
  "ETV6-RUNX1_or_like"       = "#b897bd",
  "ETV6::RUNX1"              = "#b897bd",
  "ETV6::RUNX1 (like) B-ALL" = "#b897bd",
  "ETV6::RUNX1-like"         = "#b897bd",
  "Hyperdiploidy"            = "#b162bd",
  "Hyperdiploid"             = "#b162bd",
  "Low hyperdiploid"         = "#b162bd",
  "High hyperdiploid B-ALL"  = "#b162bd",
  "Hypodiploidy"             = "#7B2FBE",
  "Low hypodiploid"          = "#7B2FBE",
  "Low hypodiploid B-ALL"    = "#7B2FBE",
  "Near haploid"             = "#7B2FBE",
  "Near haploid B-ALL"       = "#7B2FBE",
  "iAMP21"                   = "#984ea3",
  "iAMP21 B-ALL"             = "#984ea3",
  "IKZF1_N159Y"              = "#2D1B69",
  "IKZF1 N159Y"              = "#2D1B69",
  "KMT2A"                    = "#6C6C9D",
  "KMT2A_or_like"            = "#6C6C9D",
  "KMT2A rearranged (like) B-ALL" = "#6C6C9D",
  "MEF2D_rearrangement"      = "#5C4A72",
  "MEF2D"                    = "#5C4A72",
  "NUTM1"                    = "#DEC4F0",
  "NUTM1_rearrangement"      = "#DEC4F0",
  "PAX5_alteration"          = "#4B0082",
  "PAX5alt"                  = "#4B0082",
  "PAX5 altered B-ALL"       = "#4B0082",
  "PAX5::ETV6"               = "#4B0082",
  "PAX5_P80R"                = "#4B0082",
  "PAX5 P80R"                = "#4B0082",
  "TCF3-PBX1"                = "#BC80BD",
  "TCF3::PBX1"               = "#BC80BD",
  "TCF3::PBX1 positive B-ALL"= "#BC80BD",
  "ZNF384_or_like"           = "#86608E",
  "ZNF384 Group"             = "#86608E",
  "ZNF384"                   = "#86608E",
  "ZNF384 rearranged B-ALL"  = "#86608E",
  "HLF"                      = "#b162bd",

  # T-ALL
  "T-cell"                   = "#15a2a2",
  "T-ALL"                    = "#15a2a2",
  "HOXA"                     = "#43B3AE",
  "HOXA9/10 TCR"             = "#43B3AE",
  "LMO1/2"                   = "#87decd",
  "LMO2_LYL1"                = "#87decd",
  "NKX2_1"                   = "#0A9396",
  "TAL1"                     = "#00A693",
  "TAL1 \u03b1\u03b2-like"   = "#00A693",
  "TAL-LMO subtype T-ALL"    = "#00A693",
  "TAL1 DP-like"             = "#00A693",
  "TAL2"                     = "#94D2BD",
  "TLX1"                     = "#39A78D",
  "TLX3"                     = "#00CCCC",
  "TLX3 subtype T-ALL"       = "#00CCCC",
  "BCL11B"                   = "#81D8D0",
  "DNMT3A/IDH"               = "#81D8D0",
  "HOXA9"                    = "#43B3AE",
  "BCL11B-ARID1B"            = "#81D8D0",
  "MEF2C subtype T-ALL"      = "#81D8D0",

  # AML
  "AML"                      = "#80B1D3",
  "NOS AML"                  = "#80B1D3",
  "AML-MRC (myelodysplasia-related changes)" = "#80B1D3",
  "AML.MRC.5."               = "#80B1D3",
  "RUNX1 rearranged AML"     = "#6CB4EE",
  "RUNX1::RUNX1T1"           = "#6CB4EE",
  "RUNX1-RUNX1T1"            = "#6CB4EE",
  "RUNX1-RUNX1T1-like"       = "#6CB4EE",
  "RUNX1::RUNX1T1-like"      = "#6CB4EE",
  "t(8;21)(q22;q22.1)/RUNX1\u2212RUNX1T1" = "#6CB4EE",
  "t(8;21)(q22;q22.1)/RUNX1-RUNX1T1" = "#6CB4EE",
  "CBFB::MYH11"              = "#00B9E8",
  "CBFB::MYH11 positive AML" = "#00B9E8",
  "CBFB-MYH11"               = "#00B9E8",
  "inv(16)(p13.1q22)/CBFB\u2212MYH11" = "#00B9E8",
  "CBFB..MYH11"              = "#00B9E8",
  "CBFB-GDXY"                = "#00B9E8",
  "CBFA2T3-GLIS2"            = "#00B9E8",
  "FLT3-ITD, RUNX1 altered"  = "#6CB4EE",
  "UBTF"                     = "#6082B6",
  "UBTF, FLT3 mutated"       = "#6082B6",
  "UBTF ITD AML"             = "#6082B6",
  "SET-NUP214"               = "#9ECAE1",
  "DEK-NUP214"               = "#9ECAE1",
  "DEK::NUP214"              = "#9ECAE1",
  "NUP98r"                   = "#6BAED6",
  "NUP98-NSD1"               = "#6BAED6",
  "NUP98::NSD1 positive AML" = "#6BAED6",
  "NUP98-KDM5A"              = "#4292C6",
  "MDS, NUP98-KDM5A"         = "#4292C6",
  "AMKL_mixed"               = "#1164B4",
  "PML::RARA"                = "#2072AF",
  "PML-RARA"                 = "#2072AF",
  "GATA1 mutated AML/TMD"    = "#1164B4",
  "GATA1"                    = "#1164B4",
  "CEBPA"                    = "#7BAFD4",
  "CEBPA mutant AML"         = "#7BAFD4",
  "CEBPA bZIP indel"         = "#7BAFD4",
  "FUS_r"                    = "#468FEA",
  "ETV6::MNX1"               = "#B9D9EB",
  "ETV6r"                    = "#B9D9EB",
  "HNRNPH1::ERG"             = "#1E90FF",
  "CBFA2T3::GLIS2"           = "#6082B6",
  "KMT2Ar"                   = "#2171B5",
  "KMT2A rearranged (like) AML" = "#2171B5",
  "KMT2A other"              = "#2171B5",
  "KMT2A-PTD"                = "#2171B5",
  "HOXr"                     = "#0072BB",
  "HOX Grp 1"                = "#0072BB",
  "HOX Grp 2"                = "#0072BB",
  "HOX Grp 3"                = "#0072BB",
  "HOX Grp 4"                = "#0072BB",
  "HOX Grp 6"                = "#80B1D3",
  "HOX Grp 7"                = "#80B1D3",
  "HOX Grp 8"                = "#80B1D3",
  "HOX Grp 9"                = "#80B1D3",
  "IDH1_2"                   = "#3f96d4",
  "IDH2 mutation"            = "#3f96d4",
  "TP53_Aneuploidy"          = "#3f96d4",
  "Chromatin_Splicosome"     = "#3f96d4",
  "HOX Grp 5"                = "#80B1D3",
  "NPM1 altered AML"         = "#AFDFEF",
  "NPM1 mutation"            = "#AFDFEF",
  "NPM1"                     = "#AFDFEF",
  "ASXL1, BCOR, EZH2, RUNX1, SF3B1, SRSF2, STAG2, U2AF1, ZRSR2 mutation" = "#6082B6",
  "Normal or intermediate karyotype (<= 2 cytogenetic anomalies)"          = "#6082B6",
  "Normal or intermediate karyotype"                                        = "#6082B6",
  "Complex karyotype (TP53 mutation and/or >=3 cytogenetic anomalies)"      = "#6082B6",
  "Complex karyotype and/or TP53 mutations"                                 = "#6082B6",
  "del5q and/or del7q"       = "#6082B6",
  "ETS family (ERG, ETV6, FLI1, FEV, ETS2, MEF)" = "#6082B6",
  "ETS family"               = "#6082B6",
  "GLISr"                    = "#6082B6",
  "KAT6Ar"                   = "#6082B6",
  "MECOMr"                   = "#6082B6",
  "MECOM"                    = "#6082B6",
  "MNX1"                     = "#6082B6",
  "PICALM-MLLT10"            = "#6082B6",
  "RBM15-MKL1"               = "#6082B6",

  # Other cancers – haematological
  "APL"                                   = "#D9AC7A",
  "PML::RARA positive APL"                = "#D9AC7A",
  "JMML"                                  = "#D9AC7A",
  "Langerhans cell histiocytosis"         = "#D9AC7A",
  "Mature B-cell lymphoma"                = "#E88200",
  "Diffuse large B-NHL"                   = "#E88200",
  "Diffuse large B-cell lymphoma"         = "#E88200",
  "Primary mediastinal large B-NHL"       = "#E88200",
  "Immature B-cell lymphoma"              = "#E88200",
  "B-LBL"                                 = "#E88200",
  "Immature T-cell lymphoma"              = "#BC8B5A",
  "T-LBL"                                 = "#BC8B5A",
  "Mature T-cell lymphoma"                = "#BC8B5A",
  "ALK positive anaplastic large cell TL" = "#BC8B5A",
  "Hodgkin lymphoma"                      = "#DEA563",
  "Classical mixed cellularity HL"        = "#DEA563",
  "Classical nodular sclerosing hodgkin lymphoma" = "#DEA563",
  "Nodular lymphocyte predominant HL"     = "#DEA563",
  "Burkitt lymphoma"                      = "#C06605",
  "MYC rearranged BL"                     = "#C06605",
  "MYC rearranged Burkitt lymphoma"       = "#C06605",
  "CML"                                   = "#DBAA74",
  "BCR::ABL1 positive CML"               = "#DBAA74",

  # Solid tumors
  "Hepatoblastoma"                        = "#EFD4B5",
  "Epithelial-mesenchymal HBL"            = "#EFD4B5",
  "Epithelial HBL"                        = "#EFD4B5",
  "Lipoblastoma"                          = "#BC8B5A",
  "LB"                                    = "#BC8B5A",
  "Pheochromocytoma"                      = "#BC8B5A",
  "PHEO"                                  = "#BC8B5A",
  "Peripheral neuroblastic tumor"         = "#D9AC7A",
  "Neuroblastoma"                         = "#D9AC7A",
  "Pretreated necrotic poorly differentiated NB" = "#D9AC7A",
  "NB"                                    = "#D9AC7A",
  "Ganglio NB"                            = "#D9AC7A",
  "GN"                                    = "#D9AC7A",
  "Neurofibroma"                          = "#D19653",
  "NF"                                    = "#D19653",
  "Schwannoma"                            = "#D19653",
  "SCHW"                                  = "#D19653",
  "Cystic nephroma"                       = "#CC7D23",
  "CN"                                    = "#CC7D23",
  "Wilms tumor"                           = "#D79D5B",
  "Mixed WT"                              = "#D79D5B",
  "Blastemic WT"                          = "#D79D5B",
  "Regressive WT"                         = "#D79D5B",
  "Epithelial WT"                         = "#D79D5B",
  "Stromal WT"                            = "#D79D5B",
  "Dysgerminoma/Germinoma"                = "#E4C6A4",
  "GE"                                    = "#E4C6A4",
  "Teratoma"                              = "#C2792A",
  "Mature TRT"                            = "#C2792A",
  "Immature TRT"                          = "#C2792A",
  "Yolk sac tumor"                        = "#C2792A",
  "YST"                                   = "#C2792A",
  "Ovarian tumor"                         = "#AA6011",
  "Mucinous adenoma of the ovary"         = "#AA6011",
  "Sertoli (leydig) cell tumor"           = "#F8B162",
  "Poorly differentiated SLCT"            = "#F8B162",
  "Granulosa cell tumor"                  = "#F8B162",
  "Juvenile granulosa cell tumor"         = "#F8B162",
  "Juvenile GrCT"                         = "#F8B162",
  "Adrenal cortical carcinoma"            = "#E6C399",
  "ACC"                                   = "#E6C399",
  "positive"                              = "#F8B162",
  "PTC"                                   = "#F8B162",
  "Lymphoepithelial carcinoma"            = "#F8B162",
  "LEC"                                   = "#F8B162",
  "Malignant rhabdoid tumor"              = "#F8B162",
  "MRT"                                   = "#F8B162",
  "Small cell carcinoma of the ovary, hypercalcemic type" = "#E6C399",
  "SCCOHT"                                = "#E6C399",
  "Arteriovenous malformation"            = "#E6C399",
  "AVM"                                   = "#E6C399",
  "BCOR sarcoma"                          = "#E6C399",
  "BCOR(S)"                               = "#E6C399",
  "Ewing sarcoma"                         = "#D08938",
  "EWSR1 rearranged EWS"                  = "#D08938",
  "Osteosarcoma"                          = "#BA6C15",
  "Conventional type high grade osteosarcoma" = "#BA6C15",
  "Conventional type high grade OS"       = "#BA6C15",
  "Spindle cell sarcoma"                  = "#F69A3A",
  "NTRK sarcoma"                          = "#F69A3A",
  "Infantile fibrosarcoma"                = "#F69A3A",
  "Synovio sarcoma"                       = "#FF921E",
  "SS"                                    = "#FF921E",
  "Rhabdomyosarcoma"                      = "#BB6400",
  "Embryonal RMS"                         = "#BB6400",
  "Embryonal rhabdomyosarcoma"            = "#BB6400",
  "MYOD1 mutant RMS"                      = "#BB6400",
  "MYOD1 mutant rhabdomyosarcoma"         = "#BB6400",
  "Alveolar RMS"                          = "#BB6400",
  "Alveolar rhabdomyosarcoma"             = "#BB6400",
  "Pleuropulmonary blastoma"              = "#DF8F34",
  "PPB"                                   = "#DF8F34",

  # Neuro
  "Atypical teratoid/rhabdoid tumor"      = "#E4943F",
  "ATRT subclass SHH"                     = "#E4943F",
  "ATRT subclass TYR"                     = "#E4943F",
  "Medulloblastoma"                       = "#B3660F",
  "MB group 3"                            = "#B3660F",
  "MB group 4"                            = "#B3660F",
  "MB, WNT-activated"                     = "#B3660F",
  "MB, SSH-activated"                     = "#B3660F",
  "Medulloblastoma, SHH-activated"        = "#B3660F",
  "Meningioma"                            = "#E4943F",
  "MEN"                                   = "#E4943F",
  "Infant-type hemispheric glioma"        = "#B97022",
  "IHG"                                   = "#B97022",
  "Paediatric-type diffuse high-grade glioma" = "#CE8231",
  "Diffuse pHGG,  H3-WT and IDH-WT"       = "#CE8231",
  "Paediatric-type diffuse low-grade glioma" = "#CE8231",
  "Diffuse pLGG, MYB- or MYBL1-altered"   = "#CE8231",
  "Ganglioglioma"                         = "#CE8231",
  "GG"                                    = "#CE8231",
  "Craniopharyngioma"                     = "#CE8231",
  "Adamantinomatous CP"                   = "#CE8231",
  "Ependymoma"                            = "#DF9C53",
  "Myxopapillary EP"                      = "#DF9C53",
  "Posterior fossa A EP"                  = "#DF9C53",
  "Choroid plexus papilloma"              = "#B97022",
  "CPP"                                   = "#B97022",
  "Xanthoastrocytoma"                     = "#B97022",
  "PXA"                                   = "#B97022",
  "Pilocytic astrocytoma"                 = "#B97022",
  "PA"                                    = "#B97022",
  "Embryonal tumor multilayered rosettes" = "#B97022",
  "ETMR"                                  = "#B97022"
)

thresholds <- data.frame(
  classifier = c("ALLcatchR_lineage", "ALLcatchR_subtype", "SIGNATURE",
                 "MD-ALL_phenograph", "MD-ALL_svm", "ALLSorts",
                 "Maylis", "AMLmapR", "AttentionAML",
                 "MnM_lineage", "MnM_subtype"),
  threshold  = c(0.7, 0.5, 0.8, 0.5, 0.5, 0.5, 0.5, 0, 0.8, 0.82, 0.72)
)

classifier_order <- thresholds$classifier

# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------
message("Reading: ", opt$input)
df <- read.csv(opt$input, header = TRUE, sep = ",")
df$comb <- paste(df$classifier, df$subtype)

# Keep only classifiers present in the data (handles optional classifiers)
classifier_order <- classifier_order[classifier_order %in% unique(df$classifier)]

df_list <- split(df, df$Sample)

# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------
make_plot <- function(s, df_sub) {
  df_sub <- merge(df_sub, thresholds, by = "classifier", all.x = TRUE)

  df_sub$classifier <- factor(df_sub$classifier, levels = classifier_order)
  df_sub <- df_sub[order(df_sub$classifier), ]

  comb_levels <- unique(df_sub$comb)
  df_sub$comb <- factor(df_sub$comb, levels = comb_levels)

  threshold_segments <- do.call(rbind, lapply(unique(df_sub$classifier), function(cl) {
    combs_for_classifier <- df_sub$comb[df_sub$classifier == cl]
    x_positions <- which(comb_levels %in% combs_for_classifier)
    data.frame(
      classifier = cl,
      threshold  = unique(df_sub$threshold[df_sub$classifier == cl]),
      xmin       = min(x_positions) - 0.5,
      xmax       = max(x_positions) + 0.5
    )
  }))

  # Any subtype not in the colour map gets a fallback grey
  missing_cols <- setdiff(unique(df_sub$subtype), names(useColsSubtype))
  extra_cols   <- setNames(rep("#CCCCCC", length(missing_cols)), missing_cols)
  all_cols     <- c(useColsSubtype, extra_cols)

  ggplot(df_sub, aes(x = comb, y = score, fill = subtype)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = all_cols) +
    geom_text(aes(label = round(score, 2)), vjust = -0.5, size = 3.5) +
    geom_segment(data = threshold_segments,
                 aes(x = xmin, xend = xmax,
                     y = threshold, yend = threshold),
                 linetype = "dashed", linewidth = 0.5,
                 inherit.aes = FALSE) +
    ylim(0, 1) +
    xlab("") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    ggtitle(s)
}

message("Generating plots for ", length(df_list), " sample(s)...")
plots <- mapply(make_plot, names(df_list), df_list, SIMPLIFY = FALSE)

# ---------------------------------------------------------------------------
# Write PDF
# ---------------------------------------------------------------------------
message("Writing: ", opt$output)
pdf(opt$output, width = 10)
for (p in plots) print(p)
invisible(dev.off())
message("Done.")