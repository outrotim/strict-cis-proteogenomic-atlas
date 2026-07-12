#!/usr/bin/env Rscript

# Minimal reusable implementation of the principal Study 10 methods.
# Source summary statistics are not distributed with this repository.

suppressPackageStartupMessages({
  library(data.table)
  library(TwoSampleMR)
  library(coloc)
})

canonical_outcomes <- data.table(
  outcome = c(
    "MDD_PGC", "AD_Bellenguez", "Anxiety", "Dementia_FinnGen",
    "CogPerformance", "HeartFailure", "T2D"
  ),
  outcome_id = c(
    "ieu-b-102", "ebi-a-GCST90027158", "finn-b-KRA_PSY_ANXIETY",
    "finn-b-F5_DEMENTIA", "ebi-a-GCST006572", "ebi-a-GCST009541",
    "ebi-a-GCST90018926"
  ),
  outcome_type = c(
    "binary", "binary", "binary", "binary", "continuous", "binary", "binary"
  ),
  outcome_stratum = c(rep("neuropsychiatric", 5L), rep("cardiometabolic", 2L))
)

analysis_parameters <- list(
  cis_window_bp = 500000L,
  pqtl_p_threshold = 5e-8,
  f_threshold = 10,
  clump_r2 = 0.001,
  clump_kb = 10000L,
  coloc_p1 = 1e-4,
  coloc_p2 = 1e-4,
  coloc_p12 = c(1e-6, 1e-5, 1e-4),
  main_p12 = 1e-5,
  strong_pph4 = 0.80,
  pairwise_min_overlap = 10L,
  hyprcoloc_min_overlap = 20L
)

define_cis_regions <- function(protein_catalog, gene_coordinates,
                               window_bp = analysis_parameters$cis_window_bp) {
  catalog <- as.data.table(copy(protein_catalog))
  coordinates <- as.data.table(copy(gene_coordinates))
  required_catalog <- c("prot_id", "gene_symbol")
  required_coordinates <- c("gene_symbol", "gene_chr", "gene_start", "gene_end")
  stopifnot(all(required_catalog %chin% names(catalog)))
  stopifnot(all(required_coordinates %chin% names(coordinates)))

  coordinates <- coordinates[
    gene_chr %in% 1:22 & is.finite(gene_start) & is.finite(gene_end)
  ]
  coordinates[, gene_length := gene_end - gene_start]
  coordinates <- coordinates[, .SD[which.max(gene_length)], by = gene_symbol]
  regions <- merge(catalog, coordinates, by = "gene_symbol", all.x = TRUE)
  regions[, `:=`(
    region_start = pmax(1L, as.integer(gene_start) - window_bp),
    region_end = as.integer(gene_end) + window_bp,
    region_definition = "coding_gene_GRCh37_plus_minus_500kb"
  )]
  unique(regions, by = "prot_id")
}

filter_strict_cis <- function(pqtl, region,
                              p_threshold = analysis_parameters$pqtl_p_threshold,
                              f_threshold = analysis_parameters$f_threshold) {
  dat <- as.data.table(copy(pqtl))
  required <- c("rsid", "chr", "position", "beta", "se", "p")
  stopifnot(all(required %chin% names(dat)), nrow(region) == 1L)
  dat <- dat[
    as.integer(chr) == as.integer(region$gene_chr) &
      position >= region$region_start & position <= region$region_end &
      is.finite(beta) & is.finite(se) & se > 0 & is.finite(p) & p <= p_threshold
  ]
  dat[, F_statistic := (beta / se)^2]
  dat <- dat[F_statistic >= f_threshold]
  setorder(dat, p)
  unique(dat, by = "rsid")
}

clump_local <- function(pqtl, bfile = Sys.getenv("LD_BFILE"),
                        plink = Sys.getenv("PLINK_BIN"),
                        p_threshold = analysis_parameters$pqtl_p_threshold,
                        r2 = analysis_parameters$clump_r2,
                        kb = analysis_parameters$clump_kb) {
  dat <- as.data.table(copy(pqtl))
  if (nrow(dat) <= 1L) return(dat)
  if (!nzchar(bfile) || !nzchar(plink)) {
    stop("Set LD_BFILE and PLINK_BIN before local clumping.", call. = FALSE)
  }
  required_reference <- paste0(bfile, c(".bed", ".bim", ".fam"))
  if (any(!file.exists(c(required_reference, plink)))) {
    stop("PLINK executable or LD reference files are unavailable.", call. = FALSE)
  }

  prefix <- tempfile("strict_cis_clump_")
  on.exit(unlink(paste0(prefix, "*")), add = TRUE)
  fwrite(dat[, .(SNP = rsid, P = p)], prefix, sep = " ")
  fwrite(dat[, .(rsid)], paste0(prefix, ".extract"), col.names = FALSE)
  status <- system2(plink, c(
    "--bfile", shQuote(bfile), "--extract", shQuote(paste0(prefix, ".extract")),
    "--clump", shQuote(prefix), "--clump-p1", format(p_threshold, scientific = TRUE),
    "--clump-r2", as.character(r2), "--clump-kb", as.character(kb),
    "--out", shQuote(prefix)
  ), stdout = FALSE, stderr = FALSE)
  clumped <- paste0(prefix, ".clumped")
  if (!identical(status, 0L) || !file.exists(clumped)) {
    stop("Local PLINK clumping failed.", call. = FALSE)
  }
  retained <- fread(clumped, select = "SNP")$SNP
  dat[rsid %chin% retained]
}

fetch_exact_outcome <- function(snps, outcome_id) {
  if (!nzchar(Sys.getenv("OPENGWAS_JWT"))) {
    stop("Set OPENGWAS_JWT in the process environment.", call. = FALSE)
  }
  as.data.table(TwoSampleMR::extract_outcome_data(
    snps = unique(snps), outcomes = outcome_id, proxies = FALSE,
    splitsize = 10000L, proxy_splitsize = 500L
  ))
}

estimate_mr <- function(harmonised) {
  dat <- as.data.table(copy(harmonised))
  required <- c("beta.exposure", "se.exposure", "beta.outcome", "se.outcome")
  stopifnot(all(required %chin% names(dat)))
  dat <- dat[
    is.finite(beta.exposure) & is.finite(se.exposure) & se.exposure > 0 &
      is.finite(beta.outcome) & is.finite(se.outcome) & se.outcome > 0
  ]
  if (nrow(dat) == 0L) return(NULL)
  if (nrow(dat) == 1L) {
    fit <- TwoSampleMR::mr_wald_ratio(
      dat$beta.exposure, dat$beta.outcome, dat$se.exposure, dat$se.outcome
    )
    method <- "Wald ratio"
  } else {
    fit <- TwoSampleMR::mr_ivw(
      dat$beta.exposure, dat$beta.outcome, dat$se.exposure, dat$se.outcome
    )
    method <- "Inverse variance weighted"
  }
  data.table(
    method = method, nsnp = as.integer(fit$nsnp), beta = fit$b, se = fit$se,
    p = fit$pval, beta_low = fit$b - 1.96 * fit$se,
    beta_high = fit$b + 1.96 * fit$se
  )
}

run_pairwise_coloc <- function(pqtl, outcome, outcome_type,
                               outcome_n, outcome_case_fraction = NA_real_,
                               p12_grid = analysis_parameters$coloc_p12) {
  p <- as.data.table(copy(pqtl))
  o <- as.data.table(copy(outcome))
  shared <- merge(
    p[, .(snp = rsid, beta_p = beta, varbeta_p = se^2, n_p = n)],
    o[, .(snp = rsid, beta_o = beta, varbeta_o = se^2, maf)],
    by = "snp"
  )
  if (nrow(shared) < analysis_parameters$pairwise_min_overlap) {
    return(data.table(p12 = p12_grid, PP.H4 = NA_real_, status = "insufficient_overlap"))
  }
  dataset1 <- list(
    beta = shared$beta_p, varbeta = shared$varbeta_p, snp = shared$snp,
    N = median(shared$n_p, na.rm = TRUE), type = "quant", sdY = 1
  )
  dataset2 <- list(
    beta = shared$beta_o, varbeta = shared$varbeta_o, snp = shared$snp,
    N = outcome_n, MAF = shared$maf,
    type = if (outcome_type == "binary") "cc" else "quant"
  )
  if (outcome_type == "binary") dataset2$s <- outcome_case_fraction
  if (outcome_type != "binary") dataset2$sdY <- 1

  rbindlist(lapply(p12_grid, function(p12) {
    fit <- coloc::coloc.abf(
      dataset1, dataset2, p1 = analysis_parameters$coloc_p1,
      p2 = analysis_parameters$coloc_p2, p12 = p12
    )
    data.table(p12 = p12, PP.H4 = unname(fit$summary[["PP.H4.abf"]]))
  }))
}

score_candidates <- function(candidate_table) {
  dat <- as.data.table(copy(candidate_table))
  components <- c(
    "score_stage1_fdr", "score_global_fdr", "score_bonferroni",
    "score_multi_outcome", "score_regional_main_strong",
    "score_regional_prior_stable", "score_hypr_pqtl_strong",
    "score_eqtlgen", "score_tissue_eqtl", "score_smr_heidi",
    "score_ukbppp", "score_steiger", "score_pleiotropy",
    "score_mediation", "score_druggable"
  )
  stopifnot(all(components %chin% names(dat)))
  dat[, evidence_score := rowSums(.SD), .SDcols = components]
  dat[, evidence_score_max := 16L]
  setorder(dat, -evidence_score, gene_symbol)
  dat[, priority_rank := seq_len(.N)]
  dat[, priority_tier := fcase(
    evidence_score >= 10L, "Tier 1",
    evidence_score >= 7L, "Tier 2",
    default = "Tier 3"
  )]
  dat[]
}

self_check <- function() {
  mr <- estimate_mr(data.table(
    beta.exposure = 0.20, se.exposure = 0.04,
    beta.outcome = -0.10, se.outcome = 0.03
  ))
  stopifnot(nrow(mr) == 1L, mr$method == "Wald ratio")
  candidate <- data.table(gene_symbol = "TEST")
  for (name in c(
    "score_stage1_fdr", "score_global_fdr", "score_bonferroni",
    "score_multi_outcome", "score_regional_main_strong",
    "score_regional_prior_stable", "score_hypr_pqtl_strong",
    "score_eqtlgen", "score_tissue_eqtl", "score_smr_heidi",
    "score_ukbppp", "score_steiger", "score_pleiotropy",
    "score_mediation", "score_druggable"
  )) candidate[, (name) := 0L]
  candidate[, `:=`(score_stage1_fdr = 1L, score_multi_outcome = 2L)]
  scored <- score_candidates(candidate)
  stopifnot(scored$evidence_score == 3L, scored$priority_tier == "Tier 3")
  message("Self-check passed: MR estimator and 16-point ranking rubric.")
}

args <- commandArgs(trailingOnly = TRUE)
if ("--self-check" %chin% args) self_check()
if ("--print-outcomes" %chin% args) print(canonical_outcomes)
if (!length(args)) {
  message("Source this file to reuse the functions, or run with --self-check.")
}
