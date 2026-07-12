#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(grid)
  library(gridExtra)
})

script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_file <- if (length(script_arg)) sub("^--file=", "", script_arg[[1]]) else "R/02_rebuild_main_figures.R"
root <- normalizePath(file.path(dirname(script_file), ".."), mustWork = TRUE)
input_file <- file.path(root, "data", "main_figure_inputs.tsv")
output_dir <- Sys.getenv("OUTPUT_DIR", file.path(root, "output"))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(input_file)) stop("Missing figure input: ", input_file)

bundle <- fread(input_file, na.strings = "NA", check.names = FALSE)
read_set <- function(name) copy(bundle[dataset == name])
to_num <- function(dat, columns) {
  for (column in intersect(columns, names(dat))) set(dat, j = column, value = as.numeric(dat[[column]]))
  dat
}
to_logical <- function(dat, columns) {
  for (column in intersect(columns, names(dat))) {
    value <- tolower(as.character(dat[[column]]))
    set(dat, j = column, value = fifelse(value == "true", TRUE,
      fifelse(value == "false", FALSE, NA)))
  }
  dat
}

C <- list(
  blue = "#0072B2", sky = "#56B4E9", teal = "#009E73",
  amber = "#E69F00", coral = "#D55E00", ink = "#222222",
  mid = "#66727A", line = "#C9D0D4", pale = "#F4F6F7"
)
stratum_colors <- c(cardiometabolic = C$amber, neuropsychiatric = C$blue)
outcome_labels <- c(
  AD_Bellenguez = "Alzheimer disease", Anxiety = "Anxiety",
  CogPerformance = "Cognitive performance", Dementia_FinnGen = "All-cause dementia",
  HeartFailure = "Heart failure", MDD_PGC = "Major depression", T2D = "Type 2 diabetes"
)
theme_journal <- function(base = 7.5) {
  theme_classic(base_size = base, base_family = "sans") +
    theme(
      plot.title = element_text(face = "bold", colour = C$ink),
      plot.subtitle = element_text(colour = C$mid),
      strip.background = element_rect(fill = "#EEF2F4", colour = NA),
      strip.text = element_text(face = "bold"),
      legend.title = element_text(face = "bold"),
      panel.spacing = unit(0.35, "lines"),
      plot.margin = margin(7, 8, 7, 8)
    )
}
save_figure <- function(plot, stem, width, height) {
  ggsave(file.path(output_dir, paste0(stem, ".pdf")), plot,
    width = width, height = height, device = cairo_pdf, bg = "white")
  ggsave(file.path(output_dir, paste0(stem, ".png")), plot,
    width = width, height = height, dpi = 600, bg = "white")
}

flow <- to_num(read_set("flow"), "value")
stage1 <- read_set("stage1")
regional <- to_num(read_set("regional"), c("PP.H4", "min_PP.H4"))
regional <- to_logical(regional, c("strong_all_priors", "hypr_pqtl_strong_support"))
score_columns <- c(
  "priority_rank", "evidence_score", "evidence_score_max", "evidence_evaluable_max",
  "score_stage1_fdr", "score_global_fdr", "score_bonferroni", "score_multi_outcome",
  "score_regional_main_strong", "score_regional_prior_stable", "score_hypr_pqtl_strong",
  "score_eqtlgen", "score_tissue_eqtl", "score_smr_heidi", "score_ukbppp",
  "score_steiger", "score_pleiotropy", "score_mediation", "score_druggable"
)
score <- to_num(read_set("score"), score_columns)
score <- to_logical(score, c(
  "is_druggable", "eqtlgen_evaluable", "eqtlgen_validated", "gtex_evaluable",
  "brain_eqtl_validated", "heart_eqtl_validated", "smr_evaluable", "heidi_evaluable",
  "smr_heidi_joint_pass", "ukbppp_evaluable", "ukbppp_replicated",
  "mediation_evaluable", "mediation_supported"
))
triple <- to_num(read_set("triple_prior"), c("prior_c", "best_posterior_prob", "best_regional_prob"))
separation <- to_num(read_set("separation"), c(
  "binary_exact_p", "null_monte_carlo_p_lower_tail",
  "observed_signed_correlation_contrast", "bootstrap_ci_low", "bootstrap_ci_high"
))
external_summary <- to_num(read_set("external_summary"), c(
  "n_exact_computable", "n_exact_instrument_pass", "n_exact_joint_mr_coloc_pass"
))
external_pairs <- to_num(read_set("external_pairs"), c(
  "wald_p", "wald_p_fdr_exact", "external_coloc_main_pph4"
))
external_pairs <- to_logical(external_pairs, c(
  "exact_phenotype", "external_computable", "direction_concordant", "replication_pass",
  "external_coloc_main_strong", "external_coloc_prior_stable", "joint_mr_coloc_pass"
))
grouped <- to_num(read_set("grouped_ablation"), c(
  "baseline_rank", "baseline_score", "removed_points", "revised_score", "revised_rank", "rank_change"
))

# Figure 1: analysis universe, outcome composition, and interpretation boundary.
flow1 <- flow[stage %chin% c("Assay inventory", "Gene-mapped assays", "Proteins with instruments", "Evaluable tests")]
flow1[, x := seq_len(.N)]
p1a <- ggplot(flow1, aes(x, 1)) +
  geom_tile(width = .76, height = .46, fill = "white", colour = C$line, linewidth = .45) +
  geom_segment(data = flow1[x < max(x)], aes(x = x + .39, xend = x + .61, y = 1, yend = 1),
    arrow = arrow(length = unit(.06, "inches")), colour = C$mid, linewidth = .35) +
  geom_text(aes(label = format(value, big.mark = ",")), y = 1.05, size = 4.0, fontface = "bold") +
  geom_text(aes(label = stage), y = .89, size = 2.1, colour = C$mid) +
  scale_x_continuous(limits = c(.55, 4.45)) +
  scale_y_continuous(limits = c(.72, 1.30)) +
  labs(title = "A  Strict-cis analysis universe") + theme_void() +
  theme(plot.title = element_text(face = "bold", size = 8.5, margin = margin(b = 4)))

outcome_counts <- stage1[, .N, by = .(outcome, anchor_class)]
outcome_counts[, label := factor(outcome_labels[outcome], levels = rev(outcome_labels))]
p1b <- ggplot(outcome_counts, aes(N, label, fill = anchor_class)) +
  geom_col(width = .56) + geom_text(aes(label = N), hjust = -.25, size = 2.5) +
  scale_fill_manual(values = stratum_colors) +
  scale_x_continuous(limits = c(0, 21), expand = c(0, 0)) +
  labs(title = "B  Outcome composition", subtitle = "47 primary pairs across 45 proteins",
    x = "Protein-outcome pairs", y = NULL) + theme_journal() +
  theme(legend.position = "none", axis.text.y = element_text(size = 6))

taxonomy <- score[, .N, by = anchor_class]
taxonomy[, anchor_class := factor(anchor_class, levels = c("cardiometabolic", "neuropsychiatric"))]
p1c <- ggplot(taxonomy, aes(N, "Primary atlas", fill = anchor_class)) +
  geom_col(width = .38) +
  geom_text(aes(label = N), position = position_stack(vjust = .5), colour = "white", fontface = "bold") +
  scale_fill_manual(values = stratum_colors, labels = c("Cardiometabolic outcome", "Neuropsychiatric outcome")) +
  scale_x_continuous(limits = c(0, 45), expand = c(0, 0)) +
  labs(title = "C  Outcome-defined taxonomy", subtitle = "12 / 33; observed overlap 0",
    x = "Proteins", y = NULL, fill = NULL) + theme_journal() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), legend.position = "bottom")

p1d <- ggplot() +
  annotate("segment", x = separation$bootstrap_ci_low, xend = separation$bootstrap_ci_high,
    y = 1, yend = 1, colour = C$blue, linewidth = 1) +
  annotate("point", x = separation$observed_signed_correlation_contrast, y = 1,
    shape = 21, fill = "white", colour = C$blue, size = 2.6) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = C$mid) +
  annotate("text", x = .10, y = 1.20, hjust = 1, size = 2.5,
    label = sprintf("Null P=%.3f | Partition P=%.3f",
      separation$null_monte_carlo_p_lower_tail, separation$binary_exact_p)) +
  annotate("text", x = .10, y = .80, hjust = 1, size = 2.5, fontface = "bold", colour = C$coral,
    label = "Biological separation not supported") +
  scale_x_continuous(limits = c(-.06, .105), breaks = c(-.05, 0, .05, .10)) +
  scale_y_continuous(limits = c(.65, 1.32)) +
  labs(title = "D  Claim calibration", x = "Signed-correlation contrast", y = NULL) +
  theme_journal() + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

fig1 <- p1a / (p1b + p1c + plot_layout(widths = c(1.35, .85))) / p1d +
  plot_layout(heights = c(.72, 1.28, .80))
save_figure(fig1, "Figure1", 7.2, 6.7)

# Figure 2: evidence contraction, prior sensitivity, and evidence states.
attrition <- flow[stage %chin% c("Primary pairs", "Main-prior regional", "Discovery-stable")]
attrition[, x := seq_len(.N)]
p2a <- ggplot(attrition, aes(x, 1)) +
  geom_tile(width = .75, height = .48, fill = "white", colour = C$line, linewidth = .45) +
  geom_segment(data = attrition[x < max(x)], aes(x = x + .38, xend = x + .62, y = 1, yend = 1),
    arrow = arrow(length = unit(.06, "inches")), colour = C$mid) +
  geom_text(aes(label = value), y = 1.05, fontface = "bold", size = 4.2) +
  geom_text(aes(label = stage), y = .89, size = 2.2, colour = C$mid) +
  scale_x_continuous(limits = c(.55, 3.45)) + scale_y_continuous(limits = c(.70, 1.30)) +
  labs(title = "A  Evidence contraction") + theme_void() +
  theme(plot.title = element_text(face = "bold", size = 8.5))

strong <- regional[PP.H4 > .80]
strong[, pair_label := sprintf("%s | %s", gene_symbol, outcome_labels[outcome])]
setorder(strong, min_PP.H4, PP.H4)
strong[, pair_label := factor(pair_label, levels = pair_label)]
p2b <- ggplot(strong, aes(y = pair_label)) +
  geom_vline(xintercept = .8, linetype = "dashed", colour = C$mid) +
  geom_segment(aes(x = min_PP.H4, xend = PP.H4, yend = pair_label), colour = C$sky, linewidth = .8) +
  geom_point(aes(x = min_PP.H4), shape = 21, fill = "white", colour = C$blue, size = 1.8) +
  geom_point(aes(x = PP.H4), colour = C$blue, size = 1.8) +
  scale_x_continuous(limits = c(0, 1.02), breaks = c(0, .4, .8, 1)) +
  labs(title = "B  Regional posterior sensitivity", x = "PP.H4", y = NULL) +
  theme_journal(6.5) + theme(axis.text.y = element_text(size = 4.8))

pair_ev <- merge(
  strong[, .(prot_id, pair_id = paste(gene_symbol, outcome, sep = "__"), pair_label,
    strong_all_priors, hypr_pqtl_strong_support)],
  score[, .(prot_id, eqtlgen_evaluable, eqtlgen_validated, gtex_evaluable,
    tissue = brain_eqtl_validated | heart_eqtl_validated,
    smr_eval = smr_evaluable & heidi_evaluable, smr_heidi_joint_pass,
    ukbppp_evaluable, ukbppp_replicated, mediation_evaluable, mediation_supported)],
  by = "prot_id", all.x = TRUE
)
pair_ev <- merge(pair_ev,
  external_pairs[exact_phenotype == TRUE, .(
    pair_id, independent_evaluable = external_computable,
    independent_joint = joint_mr_coloc_pass
  )], by = "pair_id", all.x = TRUE)
state <- function(evaluable, supported) fifelse(
  is.na(evaluable) | evaluable != TRUE, "NA", fifelse(supported == TRUE, "Yes", "No")
)
matrix <- rbindlist(list(
  pair_ev[, .(pair_label, group = "Regional", layer = "Prior-stable", value = fifelse(strong_all_priors, "Yes", "No"))],
  pair_ev[, .(pair_label, group = "Regional", layer = "pQTL HyPr", value = fifelse(hypr_pqtl_strong_support, "Yes", "No"))],
  pair_ev[, .(pair_label, group = "Context", layer = "eQTLGen", value = state(eqtlgen_evaluable, eqtlgen_validated))],
  pair_ev[, .(pair_label, group = "Context", layer = "Tissue", value = state(gtex_evaluable, tissue))],
  pair_ev[, .(pair_label, group = "Context", layer = "SMR+HEIDI", value = state(smr_eval, smr_heidi_joint_pass))],
  pair_ev[, .(pair_label, group = "Context", layer = "Mediation", value = state(mediation_evaluable, mediation_supported))],
  pair_ev[, .(pair_label, group = "Platform", layer = "UKB-PPP", value = state(ukbppp_evaluable, ukbppp_replicated))],
  pair_ev[, .(pair_label, group = "External", layer = "Independent", value = state(independent_evaluable, independent_joint))]
))
matrix[, pair_label := factor(pair_label, levels = levels(strong$pair_label))]
p2c <- ggplot(matrix, aes(layer, pair_label)) +
  geom_tile(fill = "white", colour = "#E6EAED") +
  geom_point(data = matrix[value == "Yes"], colour = C$teal, size = 1.8) +
  geom_point(data = matrix[value == "No"], shape = 4, colour = C$coral, size = 1.8) +
  geom_text(data = matrix[value == "NA"], aes(label = "-"), colour = "#9AA3A8") +
  facet_grid(. ~ group, scales = "free_x", space = "free_x") +
  labs(title = "C  Evidence states remain non-interchangeable", x = NULL, y = NULL) +
  theme_journal(6.2) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 4.5),
    axis.text.y = element_blank(), axis.ticks.y = element_blank())

fig2 <- p2a / (p2b + p2c + plot_layout(widths = c(1, 1.35))) +
  plot_layout(heights = c(.68, 1.32))
save_figure(fig2, "Figure2", 7.2, 6.1)

# Figure 3: falsification, GRN convergence, and independent disease assessment.
tests <- data.table(
  test = factor(c("Signed-correlation contrast", "Exact partition", "Null simulation"),
    levels = rev(c("Signed-correlation contrast", "Exact partition", "Null simulation"))),
  estimate = c(separation$observed_signed_correlation_contrast, NA, NA),
  low = c(separation$bootstrap_ci_low, NA, NA), high = c(separation$bootstrap_ci_high, NA, NA),
  label = c(
    sprintf("%.3f (95%% CI %.3f to %.3f)", separation$observed_signed_correlation_contrast,
      separation$bootstrap_ci_low, separation$bootstrap_ci_high),
    sprintf("P=%.3f", separation$binary_exact_p),
    sprintf("P=%.3f", separation$null_monte_carlo_p_lower_tail)
  )
)
p3a <- ggplot(tests, aes(y = test)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = C$mid) +
  geom_segment(data = tests[!is.na(estimate)], aes(x = low, xend = high, yend = test), colour = C$blue) +
  geom_point(data = tests[!is.na(estimate)], aes(x = estimate), shape = 21, fill = "white", colour = C$blue) +
  geom_text(data = tests[is.na(estimate)], aes(x = .103, label = label), hjust = 1, size = 2.0) +
  geom_text(data = tests[!is.na(estimate)], aes(x = .103, label = label), hjust = 1,
    size = 2.0, position = position_nudge(y = .25)) +
  scale_x_continuous(limits = c(-.055, .105)) +
  labs(title = "A  Outcome-stratum falsification", x = "Signed-correlation contrast", y = NULL) +
  theme_journal(6.5)

grn <- triple[target_gene == "GRN"]
grn[, disease := fifelse(outcome == "AD_Bellenguez", "Alzheimer disease", "All-cause dementia")]
grn[, prior := factor(sprintf("%.3g", prior_c), levels = c("0.002", "0.02", "0.2"))]
p3b <- ggplot(grn, aes(prior, disease)) +
  geom_tile(fill = "#E6F2EF", colour = "white") +
  geom_text(aes(label = sprintf("%.3f", best_posterior_prob)), fontface = "bold", size = 2.3) +
  labs(title = "B  GRN prior stability", subtitle = "2/10 chains; shared-signal support only",
    x = "HyPrColoc prior.c", y = NULL) + theme_journal(6.5)

external_exact <- external_pairs[exact_phenotype == TRUE]
external_exact[, pair_label := fifelse(gene_symbol == "GRN", "GRN | dementia", "BAG3 | heart failure")]
external_matrix <- rbindlist(list(
  external_exact[, .(pair_label, layer = "Instrument", support = replication_pass,
    detail = sprintf("FDR P=%.3f", wald_p_fdr_exact))],
  external_exact[, .(pair_label, layer = "Regional", support = external_coloc_main_strong,
    detail = sprintf("PP.H4=%.3f", external_coloc_main_pph4))],
  external_exact[, .(pair_label, layer = "Joint MR-coloc", support = joint_mr_coloc_pass,
    detail = "Joint decision")]
))
p3c <- ggplot(external_matrix, aes(layer, pair_label)) +
  geom_tile(fill = "white", colour = "#E1E6E8") +
  geom_point(aes(shape = support, colour = support), size = 2.2) +
  geom_text(aes(label = detail), position = position_nudge(y = -.20), size = 1.8) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 4), na.value = 4) +
  scale_colour_manual(values = c(`TRUE` = C$teal, `FALSE` = C$coral), na.value = C$coral) +
  labs(title = "C  Independent disease assessment",
    subtitle = sprintf("Instrument support %d/%d; joint MR-coloc %d/%d",
      external_summary$n_exact_instrument_pass, external_summary$n_exact_computable,
      external_summary$n_exact_joint_mr_coloc_pass, external_summary$n_exact_computable),
    x = NULL, y = NULL, shape = NULL, colour = NULL) +
  theme_journal(6.5) + theme(legend.position = "none")

fig3 <- (p3a + p3b) / p3c + plot_layout(heights = c(1, .9))
save_figure(fig3, "Figure3", 7.2, 4.8)

# Figure 4: target prioritization and grouped evidence-family ablation.
genetic_components <- c(
  "score_stage1_fdr", "score_global_fdr", "score_bonferroni", "score_multi_outcome",
  "score_regional_main_strong", "score_regional_prior_stable", "score_hypr_pqtl_strong",
  "score_eqtlgen", "score_tissue_eqtl", "score_smr_heidi", "score_ukbppp",
  "score_steiger", "score_pleiotropy"
)
score[, genetic_support := rowSums(.SD, na.rm = TRUE), .SDcols = genetic_components]
score[, drug_axis := fifelse(is_druggable == TRUE, druggable_tier_label, "No annotation")]
drug_levels <- c("No annotation", "Tier 3B", "Tier 3A", "Tier 2", "Tier 1")
score[, drug_axis := factor(drug_axis, levels = drug_levels)]
score[, plot_x := genetic_support + ((priority_rank %% 5) - 2) * .035]
score[, plot_y := as.numeric(drug_axis) + ((priority_rank %% 3) - 1) * .035]
top_targets <- c("GRN", "RMDN1", "SIRPA", "SERPING1", "BAG3")
candidate_label_spec <- data.table(
  gene_symbol = top_targets,
  label_x = c(9.98, 8.95, 7.62, 7.95, 6.95),
  label_y = c(2.82, .82, 2.82, 3.26, 1.22),
  label_hjust = c(0, .5, .5, .5, .5)
)
candidate_labels <- merge(
  score[gene_symbol %chin% top_targets, .(gene_symbol, plot_x, plot_y)],
  candidate_label_spec, by = "gene_symbol", sort = FALSE
)
p4a <- ggplot(score, aes(plot_x, plot_y, colour = anchor_class, shape = is_druggable)) +
  geom_point(alpha = .55, size = 1.6) +
  geom_point(data = score[gene_symbol %chin% top_targets], shape = 21, fill = NA,
    colour = C$ink, size = 2.8) +
  geom_segment(data = candidate_labels,
    aes(x = plot_x, y = plot_y, xend = label_x, yend = label_y),
    inherit.aes = FALSE, colour = C$mid, linewidth = .22) +
  geom_text(data = candidate_labels,
    aes(x = label_x, y = label_y, label = gene_symbol, hjust = label_hjust),
    inherit.aes = FALSE, colour = C$ink, fontface = "italic", size = 2.0) +
  scale_colour_manual(values = stratum_colors) +
  scale_shape_manual(values = c(`FALSE` = 1, `TRUE` = 16)) +
  scale_x_continuous(breaks = 1:10, limits = c(.5, 10.5)) +
  scale_y_continuous(breaks = seq_along(drug_levels), labels = drug_levels,
    limits = c(.65, length(drug_levels) + .35)) +
  labs(title = "A  Genetic support and druggability annotation",
    x = "Genetic-support components", y = NULL, colour = "Outcome stratum", shape = "Positive annotation") +
  theme_journal(6.5) + theme(legend.position = "bottom")

family_labels <- c(
  discovery_significance = "Discovery", regional_shared_signal = "Regional",
  molecular_context = "Context", platform_and_bias = "Platform-bias",
  translation_annotation = "Translation"
)
ablation <- grouped[gene_symbol %chin% top_targets]
ablation[, family := factor(family_labels[removed_family], levels = family_labels)]
ablation[, row_label := sprintf("%s (baseline %d)", gene_symbol, baseline_rank)]
ablation[, cell_label := sprintf("%d\n(%+d)", revised_rank, rank_change)]
limit <- max(abs(ablation$rank_change), na.rm = TRUE)
p4b <- ggplot(ablation, aes(family, row_label, fill = rank_change)) +
  geom_tile(colour = "white") + geom_text(aes(label = cell_label), size = 2.0, fontface = "bold") +
  scale_fill_gradient2(low = C$blue, mid = "#F1F3F4", high = C$coral,
    midpoint = 0, limits = c(-limit, limit)) +
  labs(title = "B  Rank sensitivity to evidence-family removal",
    subtitle = "Cell: revised rank (change); positive change = worse rank",
    x = NULL, y = NULL, fill = "Rank change") + theme_journal(6.2) +
  theme(axis.text.x = element_text(size = 5), axis.text.y = element_text(size = 5),
    legend.position = "bottom")

fig4 <- p4a + p4b + plot_layout(widths = c(.92, 1.18))
save_figure(fig4, "Figure4", 7.2, 4.9)

expected <- file.path(output_dir, c(
  paste0("Figure", 1:4, ".pdf"), paste0("Figure", 1:4, ".png")
))
if (any(!file.exists(expected))) stop("One or more figure outputs were not created.")
message("Created Figure 1-4 as PDF and 600-dpi PNG in: ", output_dir)
