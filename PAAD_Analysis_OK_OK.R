# ============================================================
# TCGA-PAAD survival analysis (download via TCGAbiolinks)
# Genes: BAZ2A, SCARB2, ARPC1A, GSPT1
#
# Goals
#   (1) Single-gene survival association (HIGH/LOW by median; log2(TPM+1); 0–60 months)
#   (2) Assess whether gene–survival association persists after clinical adjustment
#       - Clinical feature selection via elastic net Cox (lambda.min / lambda.1se)
#       - Adjusted Cox per gene (Route B simple + 2b-plus using EN lambda.min terms)
#   (3) Joint 4-gene signature via GSVA (custom gene set) + survival (median HIGH/LOW)
#
# Outputs
#   OUTDIR/
#     01_data/
#     02_qc/
#     03_gene_survival/
#     04_clinical_selection/
#     05_adjusted_models/
#     06_gsva_signature/
#
# Notes
#   - This script does NOT use any local paths. It runs in the current working directory.
#   - TCGA downloads can take time and require internet access.
#   - Survival is restricted to MAX_FU_MONTHS (default 60 months) to reduce tail effects.
#   - All survival plots: HIGH = red, LOW = blue, with larger fonts and log-rank p-value.
#
# Author: (your lab)
# ============================================================

# -------------------------
# 1) Libraries
# -------------------------
suppressPackageStartupMessages({
  library(tidyverse)
  library(openxlsx)
  library(survival)
  library(survminer)
  library(broom)
  library(glmnet)
  library(GSVA)
  library(matrixStats)
  library(TCGAbiolinks)
  library(SummarizedExperiment)
})

# -------------------------
# 2) User parameters
# -------------------------
GENES <- c("BAZ2A", "SCARB2", "ARPC1A", "GSPT1")
MAX_FU_MONTHS <- 60
SEED <- 123

# Output root + subfolders
OUTDIR <- "PAAD_TCGA_pipeline_median_60m"
DIR_01_DATA   <- file.path(OUTDIR, "01_data")
DIR_02_QC     <- file.path(OUTDIR, "02_qc")
DIR_03_GENE   <- file.path(OUTDIR, "03_gene_survival")
DIR_04_CLIN   <- file.path(OUTDIR, "04_clinical_selection")
DIR_05_ADJ    <- file.path(OUTDIR, "05_adjusted_models")
DIR_06_GSVA   <- file.path(OUTDIR, "06_gsva_signature")

dir.create(DIR_01_DATA, showWarnings = FALSE, recursive = TRUE)
dir.create(DIR_02_QC,   showWarnings = FALSE, recursive = TRUE)
dir.create(DIR_03_GENE, showWarnings = FALSE, recursive = TRUE)
dir.create(DIR_04_CLIN, showWarnings = FALSE, recursive = TRUE)
dir.create(DIR_05_ADJ,  showWarnings = FALSE, recursive = TRUE)
dir.create(DIR_06_GSVA, showWarnings = FALSE, recursive = TRUE)

save_xlsx <- function(x, path) openxlsx::write.xlsx(x, path, overwrite = TRUE)
save_rds  <- function(x, path) saveRDS(x, path)

# -------------------------
# 3) Helper functions
# -------------------------

# Safe numeric coercion (avoids breaking pipelines)
as_num <- function(x) suppressWarnings(as.numeric(as.character(x)))

# TCGA sample barcode -> patient barcode (first 12 chars)
sample_to_patient <- function(x) substr(x, 1, 12)

# Restrict follow-up (months)
add_restricted_fu <- function(df, time_col, event_col, max_m = MAX_FU_MONTHS) {
  df %>%
    dplyr::mutate(
      time_restricted  = pmin(.data[[time_col]], max_m),
      event_restricted = ifelse(.data[[time_col]] <= max_m, .data[[event_col]], 0)
    )
}

# HIGH/LOW by median (returns factor LOW/HIGH)
median_group <- function(x) {
  med <- median(x, na.rm = TRUE)
  factor(ifelse(x > med, "HIGH", "LOW"), levels = c("LOW", "HIGH"))
}

# Build OS in months from TCGA-like clinical columns (days_to_death / days_to_last_followup)
# Returns OS_days, OS_event (0/1), OS_months; filters invalid rows.
build_os_from_tcga_clinical <- function(df,
                                        vital_col = "vital_status",
                                        death_col = "days_to_death",
                                        fu_col    = "days_to_last_followup") {
  df <- df %>%
    dplyr::mutate(
      vital_status = as.character(.data[[vital_col]]),
      days_to_death = as_num(.data[[death_col]]),
      days_to_last_followup = as_num(.data[[fu_col]]),
      
      OS_days = dplyr::if_else(!is.na(days_to_death), days_to_death, days_to_last_followup),
      OS_event = dplyr::case_when(
        toupper(vital_status) %in% c("DEAD", "DECEASED") ~ 1,
        toupper(vital_status) %in% c("ALIVE", "LIVING") ~ 0,
        TRUE ~ NA_real_
      ),
      OS_months = OS_days / 30.44
    ) %>%
    dplyr::filter(!is.na(OS_months), OS_months > 0, !is.na(OS_event))
  
  df
}

# Collapse duplicated gene symbols by median across samples (counts rounded)
dedup_genes_by_median <- function(mat, round_counts = FALSE) {
  stopifnot(is.matrix(mat))
  df <- as.data.frame(mat, check.names = FALSE)
  df$gene <- rownames(mat)
  
  out <- df %>%
    dplyr::group_by(gene) %>%
    dplyr::summarise(dplyr::across(where(is.numeric), median, na.rm = TRUE), .groups = "drop")
  
  rn <- out$gene
  out <- out %>% dplyr::select(-gene)
  out <- as.data.frame(out, check.names = FALSE)
  
  if (round_counts) out[] <- lapply(out, function(x) round(x))
  
  m <- as.matrix(out)
  rownames(m) <- rn
  m
}

# Collapse duplicated patient columns by row median
collapse_cols_by_patient_median <- function(mat, patient_ids) {
  stopifnot(ncol(mat) == length(patient_ids))
  u <- unique(patient_ids)
  
  if (length(u) == length(patient_ids)) {
    colnames(mat) <- patient_ids
    return(mat)
  }
  
  out <- matrix(NA_real_, nrow = nrow(mat), ncol = length(u))
  rownames(out) <- rownames(mat)
  colnames(out) <- u
  
  idx <- split(seq_along(patient_ids), patient_ids)
  for (pid in names(idx)) {
    ii <- idx[[pid]]
    if (length(ii) == 1) out[, pid] <- mat[, ii]
    else out[, pid] <- matrixStats::rowMedians(mat[, ii, drop = FALSE], na.rm = TRUE)
  }
  out
}

# Simple imputation for model.matrix:
# - character -> factor
# - factor NA -> "Missing"
# - numeric NA -> median + missing indicator
# - drop constant columns after imputation
impute_simple <- function(df) {
  df <- df %>% dplyr::mutate(dplyr::across(where(is.character), as.factor))
  
  fac_cols <- names(df)[vapply(df, is.factor, logical(1))]
  for (v in fac_cols) {
    x <- droplevels(df[[v]])
    if (any(is.na(x))) {
      x <- addNA(x)
      levels(x)[is.na(levels(x))] <- "Missing"
      x[is.na(x)] <- "Missing"
      x <- droplevels(x)
    }
    df[[v]] <- x
  }
  
  num_cols <- names(df)[vapply(df, is.numeric, logical(1))]
  for (v in num_cols) {
    miss <- is.na(df[[v]])
    if (any(miss)) {
      df[[paste0(v, "_missing")]] <- as.integer(miss)
      med <- median(df[[v]], na.rm = TRUE)
      if (!is.finite(med)) med <- 0
      df[[v]][miss] <- med
    }
  }
  
  # drop constants
  all_cols <- names(df)
  keep <- vapply(all_cols, function(v) {
    x <- df[[v]]
    if (is.factor(x)) return(nlevels(droplevels(x)) >= 2)
    if (is.numeric(x)) {
      x2 <- x[is.finite(x)]
      if (length(x2) <= 1) return(FALSE)
      return(stats::var(x2, na.rm = TRUE) > 0)
    }
    TRUE
  }, logical(1))
  
  df[, keep, drop = FALSE]
}

# Uniform KM plot style: HIGH red, LOW blue; larger fonts; always show log-rank p-value
km_plot_HL <- function(df, time_col, event_col, group_col, title, out_png, break_by = 12) {
  d <- df %>%
    dplyr::filter(!is.na(.data[[time_col]]), !is.na(.data[[event_col]]), !is.na(.data[[group_col]])) %>%
    dplyr::mutate(group_tmp = factor(.data[[group_col]], levels = c("LOW", "HIGH")))
  
  if (nrow(d) < 20 || nlevels(d$group_tmp) < 2) return(NULL)
  
  fml <- as.formula(paste0("Surv(", time_col, ", ", event_col, ") ~ group_tmp"))
  fit <- survival::survfit(fml, data = d)
  sd  <- survival::survdiff(fml, data = d)
  p_lr <- 1 - pchisq(sd$chisq, df = length(sd$n) - 1)
  p_txt <- paste0("Log-rank p = ", format.pval(p_lr, digits = 3, eps = 1e-3))
  
  p <- survminer::ggsurvplot(
    fit, data = d,
    title = title,
    pval = p_txt,
    xlab = "Months",
    break.time.by = break_by,
    risk.table = "abs_pct",
    conf.int = TRUE,
    palette = c("LOW" = "blue", "HIGH" = "red"),
    ggtheme = theme_light(base_size = 16),
    risk.table.height = 0.25,
    risk.table.fontsize = 5,
    fontsize = 5 # p-value size (survminer)
  )
  
  ggsave(out_png, p$plot, width = 9, height = 6.5, dpi = 300)
  list(
    plot = p,
    logrank_p = p_lr,
    n = nrow(d),
    events = sum(d[[event_col]] == 1, na.rm = TRUE)
  )
}

# -------------------------
# 4) Download TCGA-PAAD RNA-seq (STAR - Counts) + clinical (TCGAbiolinks)
# -------------------------

message("Downloading / preparing TCGA-PAAD RNA-seq (STAR - Counts) ...")

query_exp <- TCGAbiolinks::GDCquery(
  project = "TCGA-PAAD",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  experimental.strategy = "RNA-Seq",
  workflow.type = "STAR - Counts"
)

TCGAbiolinks::GDCdownload(query_exp)
PAAD_rnaseq <- TCGAbiolinks::GDCprepare(query_exp)

# Save raw SummarizedExperiment
save_rds(PAAD_rnaseq, file.path(DIR_01_DATA, "PAAD_RNAseq_STARcounts_SummarizedExperiment.rds"))

# Extract TPM and counts
tpm_all    <- SummarizedExperiment::assay(PAAD_rnaseq, "tpm_unstrand")
counts_all <- SummarizedExperiment::assay(PAAD_rnaseq, "unstranded")

# Convert ENSEMBL -> gene symbols when needed
if (identical(SummarizedExperiment::rowData(PAAD_rnaseq)$gene_id, rownames(tpm_all))) {
  rownames(tpm_all) <- SummarizedExperiment::rowData(PAAD_rnaseq)$gene_name
}
if (identical(SummarizedExperiment::rowData(PAAD_rnaseq)$gene_id, rownames(counts_all))) {
  rownames(counts_all) <- SummarizedExperiment::rowData(PAAD_rnaseq)$gene_name
}

# Restrict to primary tumor samples (01A)
col_prim <- grep("^.{13}01A", colnames(counts_all))
counts_primary <- counts_all[, col_prim, drop = FALSE]
tpm_primary    <- tpm_all[, col_prim, drop = FALSE]

# Deduplicate gene symbols (median across duplicated symbols)
counts_dedup <- dedup_genes_by_median(counts_primary, round_counts = TRUE)
tpm_dedup    <- dedup_genes_by_median(tpm_primary,    round_counts = FALSE)

# Collapse duplicate patients (if multiple aliquots) by median across columns after converting sample->patient
pat_ids <- sample_to_patient(colnames(tpm_dedup))
tpm_pat    <- collapse_cols_by_patient_median(tpm_dedup, pat_ids)
counts_pat <- collapse_cols_by_patient_median(counts_dedup, pat_ids)

# Save matrices
save_rds(tpm_pat,    file.path(DIR_01_DATA, "TPM_PAAD_primary_dedup_patient.rds"))
save_rds(counts_pat, file.path(DIR_01_DATA, "Counts_PAAD_primary_dedup_patient.rds"))
save_xlsx(as.data.frame(tpm_pat) %>% rownames_to_column("gene"), file.path(DIR_01_DATA, "TPM_PAAD_OK.xlsx"))
save_xlsx(as.data.frame(counts_pat) %>% rownames_to_column("gene"), file.path(DIR_01_DATA, "Counts_PAAD_OK.xlsx"))

# Clinical (patient)
message("Downloading / preparing TCGA-PAAD clinical (patient) ...")

query_clin <- TCGAbiolinks::GDCquery(
  project = "TCGA-PAAD",
  data.category = "Clinical",
  data.format = "bcr xml"
)
TCGAbiolinks::GDCdownload(query_clin)
clinical <- TCGAbiolinks::GDCprepare_clinic(query_clin, clinical.info = "patient")

save_xlsx(clinical, file.path(DIR_01_DATA, "clinical_PAAD_patient.xlsx"))
save_rds(clinical,  file.path(DIR_01_DATA, "clinical_PAAD_patient.rds"))

# -------------------------
# 5) Build analysis cohort: merge clinical + gene TPM
# -------------------------

# Create patient barcode column consistent with expression matrices
# Most clinical patient identifiers are already TCGA-XX-XXXX (12 chars); keep as-is.
# If your clinical uses submitter_id, it is typically the patient barcode.
if ("bcr_patient_barcode" %in% colnames(clinical)) {
  clinical$patient_id <- as.character(clinical$bcr_patient_barcode)
} else if ("submitter_id" %in% colnames(clinical)) {
  clinical$patient_id <- as.character(clinical$submitter_id)
} else {
  stop("Clinical table does not contain bcr_patient_barcode nor submitter_id. Please inspect clinical columns.")
}

# Build OS (months) + restricted FU
clin0 <- clinical %>%
  dplyr::mutate(patient_id = as.character(patient_id)) %>%
  build_os_from_tcga_clinical(
    vital_col = "vital_status",
    death_col = "days_to_death",
    fu_col    = "days_to_last_followup"
  ) %>%
  add_restricted_fu(time_col = "OS_months", event_col = "OS_event", max_m = MAX_FU_MONTHS)

# Merge 4 genes TPM into clinical table (patient-level)
common_patients <- intersect(clin0$patient_id, colnames(tpm_pat))
message("Patients with BOTH clinical OS and primary-tumor TPM: ", length(common_patients))

clin <- clin0 %>%
  dplyr::filter(patient_id %in% common_patients) %>%
  dplyr::arrange(patient_id)

# Add TPM (raw) for the 4 genes
missing_genes <- setdiff(GENES, rownames(tpm_pat))
if (length(missing_genes) > 0) {
  stop("Genes not found in TPM matrix rownames: ", paste(missing_genes, collapse = ", "))
}

for (g in GENES) {
  clin[[paste0(g, "_TPM")]]  <- as.numeric(tpm_pat[g, clin$patient_id])
  clin[[paste0(g, "_log2")]] <- log2(clin[[paste0(g, "_TPM")]] + 1)
  clin[[paste0(g, "_HL")]]   <- median_group(clin[[paste0(g, "_log2")]])
}

# Save survival-ready cohort
save_xlsx(clin, file.path(DIR_01_DATA, "clinical_survival_ready_plus_4genes.xlsx"))
save_rds(clin,  file.path(DIR_01_DATA, "clinical_survival_ready_plus_4genes.rds"))

# -------------------------
# 6) QC
# -------------------------
qc_tbl <- tibble::tibble(
  n_total = nrow(clin),
  events_total = sum(clin$OS_event == 1, na.rm = TRUE),
  events_restricted_60m = sum(clin$event_restricted == 1, na.rm = TRUE),
  median_followup_months = median(clin$OS_months, na.rm = TRUE)
)
save_xlsx(qc_tbl, file.path(DIR_02_QC, "QC_cohort_summary.xlsx"))

miss_vars <- c("OS_months","OS_event","time_restricted","event_restricted", paste0(GENES, "_log2"))
miss_tbl <- tibble::tibble(
  variable = miss_vars,
  n_missing = sapply(miss_vars, function(v) sum(is.na(clin[[v]]))),
  pct_missing = round(100 * sapply(miss_vars, function(v) mean(is.na(clin[[v]]))), 2)
)
save_xlsx(miss_tbl, file.path(DIR_02_QC, "QC_missingness.xlsx"))

# -------------------------
# 7) Objective 1: single-gene survival (HIGH/LOW median; 0–60m)
# -------------------------
gene_surv_results <- list()

for (g in GENES) {
  grp  <- paste0(g, "_HL")
  xlog <- paste0(g, "_log2")
  
  d <- clin %>%
    dplyr::filter(!is.na(.data[[grp]]), !is.na(.data[[xlog]])) %>%
    dplyr::mutate(group_tmp = factor(.data[[grp]], levels = c("LOW","HIGH")))
  
  # KM with standard colors + larger fonts + log-rank
  km_out <- km_plot_HL(
    df = d,
    time_col = "time_restricted",
    event_col = "event_restricted",
    group_col = grp,
    title = paste0("TCGA-PAAD OS (0–", MAX_FU_MONTHS, "m): ", g, " (median HIGH/LOW; log2(TPM+1))"),
    out_png = file.path(DIR_03_GENE, paste0("KM_", g, "_median_HL.png"))
  )
  
  # Cox univariate: HIGH vs LOW
  m1 <- survival::coxph(Surv(time_restricted, event_restricted) ~ group_tmp, data = d)
  tt1 <- broom::tidy(m1, exponentiate = TRUE, conf.int = TRUE) %>%
    dplyr::mutate(gene = g, model = "Cox_univ_median_HL",
                  n = nrow(d), events = sum(d$event_restricted == 1, na.rm = TRUE)) %>%
    dplyr::relocate(gene, model, n, events)
  
  # Cox univariate: continuous log2(TPM+1)
  m2 <- survival::coxph(Surv(time_restricted, event_restricted) ~ .data[[xlog]], data = d)
  tt2 <- broom::tidy(m2, exponentiate = TRUE, conf.int = TRUE) %>%
    dplyr::mutate(gene = g, model = "Cox_univ_continuous_log2",
                  n = nrow(d), events = sum(d$event_restricted == 1, na.rm = TRUE)) %>%
    dplyr::relocate(gene, model, n, events)
  
  gene_surv_results[[g]] <- dplyr::bind_rows(tt1, tt2)
}

gene_surv_tbl <- dplyr::bind_rows(gene_surv_results) %>% dplyr::arrange(p.value)
save_xlsx(gene_surv_tbl, file.path(DIR_03_GENE, "Gene_survival_unadjusted_results.xlsx"))

# -------------------------
# 8) Objective 2: clinical conditioning
#   8.1 clinical variable selection via elastic net Cox
#   8.2 adjusted models (Route B simple + 2b-plus)
# -------------------------

# 8.1 Candidate clinical covariates: exclude IDs, survival, gene columns
exclude_cols <- c(
  "patient_id",
  "OS_days","OS_months","OS_event","time_restricted","event_restricted",
  paste0(GENES,"_TPM"), paste0(GENES,"_log2"), paste0(GENES,"_HL")
)
candidate_cols <- setdiff(colnames(clin), exclude_cols)

clin_for_sel <- clin %>%
  dplyr::select(dplyr::all_of(candidate_cols), time_restricted, event_restricted) %>%
  dplyr::filter(complete.cases(time_restricted, event_restricted))

clin_X <- clin_for_sel %>%
  dplyr::select(-time_restricted, -event_restricted) %>%
  impute_simple()

drop_report <- tibble::tibble(
  n_rows = nrow(clin_X),
  n_cols_after_impute = ncol(clin_X)
)
save_xlsx(drop_report, file.path(DIR_04_CLIN, "Clinical_design_precheck.xlsx"))
save_xlsx(tibble::tibble(final_candidate_columns = colnames(clin_X)),
          file.path(DIR_04_CLIN, "Clinical_vars_kept_before_modelmatrix.xlsx"))

if (ncol(clin_X) == 0) {
  warning("No clinical covariates left after cleaning/imputation. Elastic net will be skipped.")
  X <- matrix(nrow = nrow(clin_for_sel), ncol = 0)
  colnames(X) <- character(0)
} else {
  X <- model.matrix(~ ., data = clin_X, na.action = na.pass)
  X <- X[, colnames(X) != "(Intercept)", drop = FALSE]
}

y <- survival::Surv(clin_for_sel$time_restricted, clin_for_sel$event_restricted)
stopifnot(nrow(X) == length(y))

save_xlsx(tibble::tibble(n_rows = nrow(X), n_cols = ncol(X)),
          file.path(DIR_04_CLIN, "Design_matrix_dimensions.xlsx"))

if (ncol(X) == 0) {
  SELECTED_TERMS <- character(0)
  save_xlsx(tibble::tibble(
    selected_lambda_1se = "",
    selected_lambda_min = "",
    n_selected_1se = 0,
    n_selected_min = 0
  ), file.path(DIR_04_CLIN, "Clinical_variable_selection_elasticnet.xlsx"))
  
  save_xlsx(tibble::tibble(term = character(), beta_lambda_min = numeric(), beta_lambda_1se = numeric()),
            file.path(DIR_04_CLIN, "ElasticNet_coefficients.xlsx"))
} else {
  set.seed(SEED)
  cv <- glmnet::cv.glmnet(x = X, y = y, family = "cox", alpha = 0.5, nfolds = 5)
  
  coef_1se <- as.matrix(coef(cv, s = "lambda.1se"))
  coef_min <- as.matrix(coef(cv, s = "lambda.min"))
  
  sel_1se <- rownames(coef_1se)[abs(coef_1se[,1]) > 0]
  sel_min <- rownames(coef_min)[abs(coef_min[,1]) > 0]
  
  save_xlsx(
    tibble::tibble(
      selected_lambda_1se = paste(sel_1se, collapse = "; "),
      selected_lambda_min = paste(sel_min, collapse = "; "),
      n_selected_1se = length(sel_1se),
      n_selected_min = length(sel_min)
    ),
    file.path(DIR_04_CLIN, "Clinical_variable_selection_elasticnet.xlsx")
  )
  
  coef_tbl <- tibble::tibble(
    term = rownames(coef_min),
    beta_lambda_min = as.numeric(coef_min[,1]),
    beta_lambda_1se = as.numeric(coef_1se[rownames(coef_min),1])
  ) %>% dplyr::arrange(dplyr::desc(abs(beta_lambda_min)))
  
  save_xlsx(coef_tbl, file.path(DIR_04_CLIN, "ElasticNet_coefficients.xlsx"))
  
  # Selection rule: prefer 1se if it selects >=2 terms, else lambda.min
  SELECTED_TERMS <- if (length(sel_1se) >= 2) sel_1se else sel_min
}

# 8.2 Route B SIMPLE: adjusted Cox per gene using a minimal required set (age + sex if present)
REQ_VARS <- c("age_at_initial_pathologic_diagnosis", "age", "gender", "sex")
REQ_AVAILABLE <- intersect(REQ_VARS, colnames(clin))
save_xlsx(tibble::tibble(required_covariates_used = REQ_AVAILABLE),
          file.path(DIR_05_ADJ, "Required_covariates_used.xlsx"))

# Helper to pick base vars consistently (age/sex variants in TCGA clinical)
pick_base_vars <- function(df) {
  # age
  age_var <- intersect(c("age_at_initial_pathologic_diagnosis", "age"), colnames(df))
  age_var <- if (length(age_var) > 0) age_var[1] else character(0)
  
  # sex
  sex_var <- intersect(c("gender", "sex"), colnames(df))
  sex_var <- if (length(sex_var) > 0) sex_var[1] else character(0)
  
  c(age_var, sex_var)
}

# Route B SIMPLE
adj_results <- list()

for (g in GENES) {
  grp <- paste0(g, "_HL")
  
  d0 <- clin %>%
    dplyr::select(time_restricted, event_restricted, dplyr::all_of(grp), dplyr::all_of(candidate_cols)) %>%
    dplyr::filter(!is.na(.data[[grp]])) %>%
    dplyr::filter(complete.cases(time_restricted, event_restricted)) %>%
    dplyr::mutate(group_tmp = factor(.data[[grp]], levels = c("LOW","HIGH"))) %>%
    dplyr::select(-dplyr::all_of(grp))
  
  base_vars <- pick_base_vars(d0)
  base_vars <- intersect(base_vars, colnames(d0))
  
  d_cox <- d0 %>%
    dplyr::select(time_restricted, event_restricted, group_tmp, dplyr::all_of(base_vars))
  
  rhs <- c("group_tmp", base_vars)
  fml <- reformulate(rhs, response = "Surv(time_restricted, event_restricted)")
  
  m <- tryCatch(survival::coxph(fml, data = d_cox), error = function(e) NULL)
  
  if (is.null(m)) {
    adj_results[[g]] <- tibble::tibble(
      gene = g, model = "Cox_adjusted_RouteB_SIMPLE_FAILED",
      n = nrow(d_cox), events = sum(d_cox$event_restricted == 1, na.rm = TRUE),
      term = "group_tmpHIGH", estimate = NA_real_, conf.low = NA_real_, conf.high = NA_real_, p.value = NA_real_
    )
    next
  }
  
  tt_all <- broom::tidy(m, exponentiate = TRUE, conf.int = TRUE)
  save_xlsx(tt_all %>% dplyr::mutate(gene = g), file.path(DIR_05_ADJ, paste0("Cox_full_terms_", g, ".xlsx")))
  
  ph_obj <- survival::cox.zph(m)
  ph_df <- as.data.frame(ph_obj$table) %>%
    tibble::rownames_to_column("term") %>%
    dplyr::mutate(gene = g) %>%
    dplyr::relocate(gene, term)
  save_xlsx(ph_df, file.path(DIR_05_ADJ, paste0("PHtest_", g, ".xlsx")))
  
  tt_grp <- tt_all %>% dplyr::filter(term == "group_tmpHIGH")
  if (nrow(tt_grp) == 0) {
    tt_grp <- tibble::tibble(term = "group_tmpHIGH", estimate = NA_real_, conf.low = NA_real_, conf.high = NA_real_, p.value = NA_real_)
  }
  
  adj_results[[g]] <- tt_grp %>%
    dplyr::mutate(
      gene = g,
      model = "Cox_adjusted_RouteB_SIMPLE_base_only",
      n = nrow(d_cox),
      events = sum(d_cox$event_restricted == 1, na.rm = TRUE),
      n_covars_used = length(base_vars)
    ) %>%
    dplyr::relocate(gene, model, n, events, n_covars_used)
}

adj_tbl <- dplyr::bind_rows(adj_results)
if ("p.value" %in% colnames(adj_tbl)) adj_tbl <- adj_tbl %>% dplyr::arrange(p.value)
save_xlsx(adj_tbl, file.path(DIR_05_ADJ, "Gene_survival_adjusted_results_RouteB_SIMPLE.xlsx"))

# 8.2PLUS (2b-plus): base_vars + EN lambda.min dummies from ElasticNet_coefficients.xlsx
en_coef_path <- file.path(DIR_04_CLIN, "ElasticNet_coefficients.xlsx")
stopifnot(file.exists(en_coef_path))
en_coef <- openxlsx::read.xlsx(en_coef_path)

EN_TERMS_ORIG <- en_coef %>%
  dplyr::filter(!is.na(beta_lambda_min), abs(beta_lambda_min) > 0) %>%
  dplyr::pull(term) %>%
  unique()

save_xlsx(tibble::tibble(EN_terms_lambda_min = EN_TERMS_ORIG),
          file.path(DIR_05_ADJ, "EN_terms_lambda_min_used.xlsx"))

adjplus_results <- list()

for (g in GENES) {
  
  grp <- paste0(g, "_HL")
  
  d0 <- clin %>%
    dplyr::select(time_restricted, event_restricted, dplyr::all_of(grp), dplyr::all_of(candidate_cols)) %>%
    dplyr::filter(!is.na(.data[[grp]])) %>%
    dplyr::filter(complete.cases(time_restricted, event_restricted)) %>%
    dplyr::mutate(group_tmp = factor(.data[[grp]], levels = c("LOW","HIGH"))) %>%
    dplyr::select(-dplyr::all_of(grp))
  
  base_vars <- pick_base_vars(d0)
  base_vars <- intersect(base_vars, colnames(d0))
  
  d_cox <- d0 %>%
    dplyr::select(time_restricted, event_restricted, group_tmp, dplyr::all_of(base_vars))
  
  # Build EN dummy matrix from candidate_cols
  cov_df <- d0 %>% dplyr::select(dplyr::all_of(candidate_cols)) %>% impute_simple()
  
  Xorig <- model.matrix(~ ., data = cov_df, na.action = na.pass)
  Xorig <- Xorig[, colnames(Xorig) != "(Intercept)", drop = FALSE]
  
  orig_names <- colnames(Xorig)
  safe_names <- make.names(orig_names, unique = TRUE)
  colnames(Xorig) <- safe_names
  
  name_map <- tibble::tibble(original = orig_names, safe = safe_names)
  
  EN_TERMS_SAFE <- name_map$safe[match(EN_TERMS_ORIG, name_map$original)]
  EN_TERMS_SAFE <- EN_TERMS_SAFE[!is.na(EN_TERMS_SAFE)]
  EN_TERMS_SAFE <- setdiff(EN_TERMS_SAFE, c(base_vars, "group_tmp"))
  EN_TERMS_SAFE <- intersect(colnames(Xorig), EN_TERMS_SAFE)
  
  X_en <- if (length(EN_TERMS_SAFE) > 0) Xorig[, EN_TERMS_SAFE, drop = FALSE] else NULL
  
  if (!is.null(X_en) && ncol(X_en) > 0) {
    X_en <- X_en[, setdiff(colnames(X_en), colnames(d_cox)), drop = FALSE]
    if (ncol(X_en) == 0) X_en <- NULL
  }
  
  if (!is.null(X_en) && ncol(X_en) > 0) {
    d_cox <- dplyr::bind_cols(d_cox, as.data.frame(X_en))
    rhs <- c("group_tmp", base_vars, colnames(X_en))
    model_tag <- "Cox_adjusted_2bPLUS_base+EN_lambda_min"
  } else {
    rhs <- c("group_tmp", base_vars)
    model_tag <- "Cox_adjusted_2bPLUS_base_only_EN_missing_or_dropped"
  }
  
  fml <- reformulate(rhs, response = "Surv(time_restricted, event_restricted)")
  
  m <- tryCatch(survival::coxph(fml, data = d_cox), error = function(e) e)
  
  if (inherits(m, "error")) {
    adjplus_results[[g]] <- tibble::tibble(
      gene = g, model = paste0(model_tag, "_FAILED"),
      n = nrow(d_cox),
      events = sum(d_cox$event_restricted == 1, na.rm = TRUE),
      n_covars_base = length(base_vars),
      n_EN_terms = ifelse(is.null(X_en), 0, ncol(X_en)),
      estimate = NA_real_, conf.low = NA_real_, conf.high = NA_real_, p.value = NA_real_,
      error_msg = conditionMessage(m)
    )
    next
  }
  
  tt_all <- broom::tidy(m, exponentiate = TRUE, conf.int = TRUE)
  save_xlsx(tt_all %>% dplyr::mutate(gene = g, model = model_tag),
            file.path(DIR_05_ADJ, paste0("CoxPLUS_full_terms_", g, ".xlsx")))
  
  ph_obj <- survival::cox.zph(m)
  ph_df <- as.data.frame(ph_obj$table) %>%
    tibble::rownames_to_column("term") %>%
    dplyr::mutate(gene = g, model = model_tag) %>%
    dplyr::relocate(gene, model, term)
  save_xlsx(ph_df, file.path(DIR_05_ADJ, paste0("PHtestPLUS_", g, ".xlsx")))
  
  tt_grp <- tt_all %>% dplyr::filter(term == "group_tmpHIGH")
  if (nrow(tt_grp) == 0) {
    tt_grp <- tibble::tibble(term = "group_tmpHIGH",
                             estimate = NA_real_, conf.low = NA_real_, conf.high = NA_real_, p.value = NA_real_)
  }
  
  adjplus_results[[g]] <- tt_grp %>%
    dplyr::mutate(
      gene = g,
      model = model_tag,
      n = nrow(d_cox),
      events = sum(d_cox$event_restricted == 1, na.rm = TRUE),
      n_covars_base = length(base_vars),
      n_EN_terms = ifelse(is.null(X_en), 0, ncol(X_en))
    ) %>%
    dplyr::relocate(gene, model, n, events, n_covars_base, n_EN_terms)
}

adjplus_tbl <- dplyr::bind_rows(adjplus_results)
if ("p.value" %in% colnames(adjplus_tbl)) adjplus_tbl <- adjplus_tbl %>% dplyr::arrange(p.value)
save_xlsx(adjplus_tbl, file.path(DIR_05_ADJ, "Gene_survival_adjusted_results_2bPLUS.xlsx"))

# -------------------------
# 9) Objective 3: GSVA joint signature (custom 4-gene set) + survival (median HL)
# -------------------------

# GSVA requires a gene x sample matrix. We use log2(TPM+1) of the full TPM patient matrix.
common_ids_gsva <- intersect(colnames(tpm_pat), clin$patient_id)
message("Patients available for GSVA + clinical: ", length(common_ids_gsva))

if (length(common_ids_gsva) < 30) {
  warning("Too few patients for GSVA analysis. Skipping GSVA.")
} else {
  expr_log2 <- log2(tpm_pat[, common_ids_gsva, drop = FALSE] + 1)
  
  gene_set <- list(GSVA_4GENES = GENES)
  gsva_scores <- GSVA::gsva(
    expr = expr_log2,
    gset.idx.list = gene_set,
    method = "gsva",
    kcdf = "Gaussian",
    min.sz = 1,
    max.sz = 500,
    parallel.sz = 1
  )
  
  sig_score <- as.numeric(gsva_scores["GSVA_4GENES", ])
  names(sig_score) <- colnames(gsva_scores)
  
  save_xlsx(
    tibble::tibble(patient_id = names(sig_score), GSVA_4GENES_score = sig_score),
    file.path(DIR_06_GSVA, "GSVA_4GENES_scores_per_patient.xlsx")
  )
  save_rds(gsva_scores, file.path(DIR_06_GSVA, "GSVA_4GENES_scores_matrix.rds"))
  
  clin_gsva <- clin %>%
    dplyr::filter(patient_id %in% names(sig_score)) %>%
    dplyr::mutate(
      GSVA_4GENES_score = sig_score[patient_id],
      GSVA_4GENES_HL = median_group(GSVA_4GENES_score)
    )
  
  save_xlsx(clin_gsva, file.path(DIR_06_GSVA, "clinical_plus_GSVA4genes.xlsx"))
  save_rds(clin_gsva,  file.path(DIR_06_GSVA, "clinical_plus_GSVA4genes.rds"))
  
  # KM (HIGH/LOW) with consistent colors
  km_plot_HL(
    df = clin_gsva,
    time_col = "time_restricted",
    event_col = "event_restricted",
    group_col = "GSVA_4GENES_HL",
    title = paste0("TCGA-PAAD OS (0–", MAX_FU_MONTHS, "m): GSVA 4-gene signature (median HIGH/LOW)"),
    out_png = file.path(DIR_06_GSVA, "KM_GSVA_4GENES_median_HL.png")
  )
  
  # Cox univariate (HL)
  m_g1 <- survival::coxph(Surv(time_restricted, event_restricted) ~ GSVA_4GENES_HL, data = clin_gsva)
  tt_g1 <- broom::tidy(m_g1, exponentiate = TRUE, conf.int = TRUE) %>%
    dplyr::mutate(model = "Cox_univ_GSVA_HL", n = nrow(clin_gsva),
                  events = sum(clin_gsva$event_restricted == 1, na.rm = TRUE)) %>%
    dplyr::relocate(model, n, events)
  
  # Cox univariate (continuous score)
  m_g1c <- survival::coxph(Surv(time_restricted, event_restricted) ~ GSVA_4GENES_score, data = clin_gsva)
  tt_g1c <- broom::tidy(m_g1c, exponentiate = TRUE, conf.int = TRUE) %>%
    dplyr::mutate(model = "Cox_univ_GSVA_continuous", n = nrow(clin_gsva),
                  events = sum(clin_gsva$event_restricted == 1, na.rm = TRUE)) %>%
    dplyr::relocate(model, n, events)
  
  # Cox adjusted SIMPLE (base vars only)
  base_vars_g <- pick_base_vars(clin_gsva)
  base_vars_g <- intersect(base_vars_g, colnames(clin_gsva))
  if (length(base_vars_g) > 0) {
    fml_simple <- reformulate(c("GSVA_4GENES_HL", base_vars_g),
                              response = "Surv(time_restricted, event_restricted)")
    m_g2s <- survival::coxph(fml_simple, data = clin_gsva)
    tt_g2s <- broom::tidy(m_g2s, exponentiate = TRUE, conf.int = TRUE) %>%
      dplyr::mutate(model = "Cox_adj_GSVA_HL_base_only", n = nrow(clin_gsva),
                    events = sum(clin_gsva$event_restricted == 1, na.rm = TRUE)) %>%
      dplyr::relocate(model, n, events)
  } else {
    tt_g2s <- tibble::tibble(
      model = "Cox_adj_GSVA_HL_base_only_SKIPPED_no_covars",
      n = nrow(clin_gsva),
      events = sum(clin_gsva$event_restricted == 1, na.rm = TRUE)
    )
  }
  
  # Cox adjusted 2b-plus (base + EN dummies derived from candidate_cols in clin_gsva)
  # Reuse EN_TERMS_ORIG from section 8.2PLUS.
  cov_df_g <- clin_gsva %>% dplyr::select(dplyr::all_of(candidate_cols)) %>% impute_simple()
  Xg <- model.matrix(~ ., data = cov_df_g, na.action = na.pass)
  Xg <- Xg[, colnames(Xg) != "(Intercept)", drop = FALSE]
  
  orig_names_g <- colnames(Xg)
  safe_names_g <- make.names(orig_names_g, unique = TRUE)
  colnames(Xg) <- safe_names_g
  name_map_g <- tibble::tibble(original = orig_names_g, safe = safe_names_g)
  
  EN_TERMS_SAFE_g <- name_map_g$safe[match(EN_TERMS_ORIG, name_map_g$original)]
  EN_TERMS_SAFE_g <- EN_TERMS_SAFE_g[!is.na(EN_TERMS_SAFE_g)]
  EN_TERMS_SAFE_g <- setdiff(EN_TERMS_SAFE_g, c(base_vars_g, "GSVA_4GENES_HL", "GSVA_4GENES_score"))
  EN_TERMS_SAFE_g <- intersect(colnames(Xg), EN_TERMS_SAFE_g)
  
  X_en_g <- if (length(EN_TERMS_SAFE_g) > 0) Xg[, EN_TERMS_SAFE_g, drop = FALSE] else NULL
  
  d_cox_g <- clin_gsva %>%
    dplyr::select(time_restricted, event_restricted, GSVA_4GENES_HL, dplyr::all_of(base_vars_g))
  
  if (!is.null(X_en_g) && ncol(X_en_g) > 0) {
    X_en_g <- X_en_g[, setdiff(colnames(X_en_g), colnames(d_cox_g)), drop = FALSE]
    if (ncol(X_en_g) == 0) X_en_g <- NULL
  }
  
  if (!is.null(X_en_g) && ncol(X_en_g) > 0) {
    d_cox_g <- dplyr::bind_cols(d_cox_g, as.data.frame(X_en_g))
    rhs_g <- c("GSVA_4GENES_HL", base_vars_g, colnames(X_en_g))
    tag_g <- "Cox_adj_GSVA_HL_2bPLUS_base+EN_lambda_min"
  } else {
    rhs_g <- c("GSVA_4GENES_HL", base_vars_g)
    tag_g <- "Cox_adj_GSVA_HL_2bPLUS_base_only_EN_missing_or_dropped"
  }
  
  fml_plus_g <- reformulate(rhs_g, response = "Surv(time_restricted, event_restricted)")
  m_g2p <- tryCatch(survival::coxph(fml_plus_g, data = d_cox_g), error = function(e) e)
  
  if (inherits(m_g2p, "error")) {
    tt_g2p <- tibble::tibble(
      model = paste0(tag_g, "_FAILED"),
      n = nrow(d_cox_g),
      events = sum(d_cox_g$event_restricted == 1, na.rm = TRUE),
      error_msg = conditionMessage(m_g2p)
    )
  } else {
    tt_g2p <- broom::tidy(m_g2p, exponentiate = TRUE, conf.int = TRUE) %>%
      dplyr::mutate(model = tag_g, n = nrow(d_cox_g),
                    events = sum(d_cox_g$event_restricted == 1, na.rm = TRUE)) %>%
      dplyr::relocate(model, n, events)
  }
  
  save_xlsx(dplyr::bind_rows(tt_g1, tt_g1c, tt_g2s, tt_g2p),
            file.path(DIR_06_GSVA, "GSVA_4genes_survival_results.xlsx"))
}

# -------------------------
# 10) Run summary
# -------------------------
message("============================================================")
message("Pipeline completed.")
message("Outputs saved in: ", normalizePath(OUTDIR))
message("Key outputs:")
message("  - 01_data/: clinical + TPM/counts matrices (primary tumor, dedup, patient-collapsed)")
message("  - 03_gene_survival/: KM plots + unadjusted Cox results per gene (median HL)")
message("  - 04_clinical_selection/: elastic net selection artifacts (coefficients, selected terms)")
message("  - 05_adjusted_models/: adjusted Cox per gene (Route B simple + 2b-plus) + PH tests")
message("  - 06_gsva_signature/: GSVA 4-gene score + survival models + KM plot")
message("============================================================")