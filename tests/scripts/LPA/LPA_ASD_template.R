# ============================================================
# LPA_ASD_template.R
# Requires: sex_label, result_dir set before sourcing
#           config.R sourced before sourcing this file
# ============================================================

library(tidyverse)
library(tidyLPA)
library(ggplot2)
library(ggpubr)
library(rstatix)
library(FSA)
library(patchwork)

dir.create(result_dir, showWarnings = FALSE, recursive = TRUE)

# --- Check input file exists before attempting to load ---
asd_input_path <- file.path(PATHS[[sex_label]], LPA_CONFIG$asd_input_file)

if (!file.exists(asd_input_path)) {
  stop(paste0(
    "\n================================================================\n",
    "ERROR: ASD input file not found for ", sex_label, ":\n",
    asd_input_path, "\n",
    "----------------------------------------------------------------\n",
    "Have you run 03_subset_neuromeasures.R for ", sex_label, "?\n",
    "Check that NA_SAMPLE_IDS$", sex_label, " is set in config.R\n",
    "then rerun script 03 before running this script.\n",
    "================================================================\n"
  ))
}

# --- Load inputs ---
ASD_lpa_input <- readRDS(asd_input_path)
ASD_mullen_vineland_full <- readRDS(
  file.path(PATHS[[sex_label]], "ASD_mullen_vineland_full.rds"))

cat("Running ASD LPA for", sex_label, "| n =", nrow(ASD_lpa_input), "\n")
cat("Input file:", LPA_CONFIG$asd_input_file, "\n")
cat("Variables used:\n")
print(colnames(ASD_lpa_input))

# --- Pull sex-specific settings from nested config ---
best_K               <- LPA_CONFIG$best_K_ASD[[sex_label]]
best_model_num       <- LPA_CONFIG$best_model_ASD[[sex_label]]
k_limit_m1           <- LPA_CONFIG$max_K_ASD_m1[[sex_label]]
k_limit_m2           <- LPA_CONFIG$max_K_ASD_m2[[sex_label]]
asd_classed_filename <- LPA_CONFIG$asd_classed_file[[sex_label]]

cat("Settings for", sex_label, ":\n")
cat("  Model 1 max K:", k_limit_m1, "\n")
cat("  Model 2 max K:", k_limit_m2, "\n")
cat("  Best K:       ", best_K, "\n")
cat("  Best model:   ", best_model_num, "\n")

# --- Model comparison --- run separately to avoid convergence failures ---
cat("Running Model 1 (equal variances) K = 1 to", k_limit_m1, "\n")
LPA_m1 <- ASD_lpa_input %>%
  estimate_profiles(1:k_limit_m1, models = 1)

cat("Running Model 2 (varying variances) K = 1 to", k_limit_m2, "\n")
LPA_m2 <- ASD_lpa_input %>%
  estimate_profiles(1:k_limit_m2, models = 2)

# --- Combine fit tables for reporting ---
fit_m1 <- get_fit(LPA_m1)
fit_m2 <- get_fit(LPA_m2)
fit_table <- bind_rows(fit_m1, fit_m2) %>%
  arrange(Model, Classes)

print(fit_table)
write.csv(fit_table,
          file.path(result_dir, "ASD_model_comparison_fit_table.csv"),
          row.names = FALSE)

# --- Extract best model safely by name ---
best_model_name <- paste0("model_", best_model_num, "_class_", best_K)
cat("Extracting solution:", best_model_name, "\n")

LPA_best_obj <- if (best_model_num == 1) LPA_m1 else LPA_m2

available_models <- names(LPA_best_obj)
cat("Available solutions in selected model object:\n")
print(available_models)

if (!best_model_name %in% available_models) {
  stop(paste0(
    "\n================================================================\n",
    "ERROR: Requested solution '", best_model_name, "' not found.\n",
    "----------------------------------------------------------------\n",
    "Available solutions:\n",
    paste(" -", available_models, collapse = "\n"), "\n",
    "----------------------------------------------------------------\n",
    "Update best_K_ASD$", sex_label, " and/or best_model_ASD$",
    sex_label, " in config.R\n",
    "================================================================\n"
  ))
}

best_model_data <- get_data(LPA_best_obj[[best_model_name]])
ASD_lpa_input$Class <- as.factor(best_model_data$Class)

cat("Class distribution (K =", best_K, ", Model =", best_model_num, "):\n")
print(table(ASD_lpa_input$Class))

neuro_vars    <- colnames(ASD_lpa_input)[colnames(ASD_lpa_input) != "Class"]
n_classes_ASD <- nlevels(ASD_lpa_input$Class)

cat("Number of classes:", n_classes_ASD, "\n")
cat("Statistical test:",
    if (n_classes_ASD == 2) "Wilcoxon rank-sum (non-parametric, 2 groups)" else
    "Kruskal-Wallis + Dunn (non-parametric, 3+ groups)", "\n")

# --- Profile plot ---
p_profiles <- plot_profiles(LPA_best_obj[[best_model_name]]) +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(result_dir,
                 paste0("ASD_profile_plot_k", best_K, ".pdf")),
       p_profiles, width = 14, height = 6)

# --- Line plot ---
long_lpa <- ASD_lpa_input %>%
  pivot_longer(cols = -Class, names_to = "Variable", values_to = "Score")

p_line <- ggplot(long_lpa,
                 aes(x = Variable, y = Score, group = Class, color = Class)) +
  stat_summary(fun = mean, geom = "line", linewidth = 1) +
  stat_summary(fun = mean, geom = "point") +
  theme_minimal() +
  labs(title = paste("ASD Latent Profile Plot —", sex_label),
       x = "Neurodevelopmental Measure", y = "Mean Score") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(result_dir,
                 paste0("ASD_line_plot_k", best_K, ".pdf")),
       p_line, width = 15, height = 6)

# --- Summary table ---
summary_table <- ASD_lpa_input %>%
  group_by(Class) %>%
  summarise(across(all_of(neuro_vars),
                   list(mean = ~mean(.x, na.rm = TRUE),
                        sd   = ~sd(.x,   na.rm = TRUE)),
                   .names = "{.col}_{.fn}"),
            .groups = "drop")
write.csv(summary_table,
          file.path(result_dir,
                    paste0("ASD_summary_by_class_k", best_K, ".csv")),
          row.names = FALSE)

# --- PCA plot ---
pca_result <- prcomp(ASD_lpa_input[, neuro_vars],
                     center = TRUE, scale. = TRUE)
var_exp    <- round(summary(pca_result)$importance[2, 1:2] * 100, 1)
pca_scores <- as.data.frame(pca_result$x[, 1:2])
colnames(pca_scores) <- c("PC1", "PC2")
pca_scores$Class <- ASD_lpa_input$Class

p_pca <- ggplot(pca_scores, aes(x = PC1, y = PC2, color = Class)) +
  geom_point(size = 3, alpha = 0.8) +
  # Only add ellipse if more than 1 class and enough points
  { if (n_classes_ASD > 1)
      stat_ellipse(type = "norm", linewidth = 1) } +
  labs(title = paste("PCA by LPA Class — ASD", sex_label),
       x = paste0("PC1 (", var_exp[1], "%)"),
       y = paste0("PC2 (", var_exp[2], "%)")) +
  theme_minimal()
ggsave(file.path(result_dir,
                 paste0("ASD_PCA_k", best_K, ".pdf")),
       p_pca, width = 8, height = 6)

# --- Statistical tests + boxplots ---
# Adapts automatically based on number of classes:
#   K = 2: Wilcoxon rank-sum (non-parametric equivalent of t-test)
#   K > 2: Kruskal-Wallis + Dunn pairwise (non-parametric multi-group)

plot_list <- list()

if (n_classes_ASD == 2) {

  cat("K=2 — using Wilcoxon rank-sum test\n")

  wilcox_results <- map_dfr(neuro_vars, function(v) {
    wt <- ASD_lpa_input %>%
      wilcox_test(as.formula(paste(v, "~ Class"))) %>%
      adjust_pvalue(method = "BH") %>%
      mutate(p.adj.signif = case_when(
        p.adj <= 0.001 ~ "***",
        p.adj <= 0.01  ~ "**",
        p.adj <= 0.05  ~ "*",
        p.adj <= 0.1   ~ ".",
        TRUE           ~ "ns"
      ))
    wt$Variable <- v
    wt
  })
  write.csv(wilcox_results,
            file.path(result_dir, "ASD_wilcoxon_results.csv"),
            row.names = FALSE)

  for (v in neuro_vars) {
    wilcox_tbl <- ASD_lpa_input %>%
      wilcox_test(as.formula(paste(v, "~ Class"))) %>%
      adjust_pvalue(method = "BH") %>%
      mutate(p.adj.signif = case_when(
        p.adj <= 0.001 ~ "***",
        p.adj <= 0.01  ~ "**",
        p.adj <= 0.05  ~ "*",
        p.adj <= 0.1   ~ ".",
        TRUE           ~ "ns"
      )) %>%
      add_xy_position(x = "Class")

    wilcox_p <- wilcox_results$p.adj[wilcox_results$Variable == v]

    p_box <- ggboxplot(ASD_lpa_input, x = "Class", y = v,
                       color = "Class", add = "jitter") +
      stat_pvalue_manual(wilcox_tbl, label = "p.adj.signif",
                         tip.length = 0.01, hide.ns = TRUE) +
      labs(title    = v,
           subtitle = paste0("Wilcoxon BH-adjusted p = ",
                             signif(wilcox_p, 3))) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))

    plot_list[[v]] <- p_box
  }

} else {

  cat("K>2 — using Kruskal-Wallis + Dunn pairwise tests\n")

  kw_results <- map_dfr(neuro_vars, function(v) {
    kt <- kruskal_test(ASD_lpa_input, as.formula(paste(v, "~ Class")))
    kt$Variable <- v
    kt
  })
  write.csv(kw_results,
            file.path(result_dir, "ASD_kruskal_results.csv"),
            row.names = FALSE)

  dunn_results <- map_dfr(neuro_vars, function(v) {
    dt <- ASD_lpa_input %>%
      dunn_test(as.formula(paste(v, "~ Class")), p.adjust.method = "BH")
    dt$Variable <- v
    dt
  })
  write.csv(dunn_results,
            file.path(result_dir, "ASD_dunn_results.csv"),
            row.names = FALSE)

  for (v in neuro_vars) {
    pairwise_tbl <- ASD_lpa_input %>%
      dunn_test(as.formula(paste(v, "~ Class")), p.adjust.method = "BH") %>%
      add_xy_position(x = "Class")
    kw_p <- kw_results$p[kw_results$Variable == v]

    p_box <- ggboxplot(ASD_lpa_input, x = "Class", y = v,
                       color = "Class", add = "jitter") +
      stat_pvalue_manual(pairwise_tbl,
                         label      = "p.adj.signif",
                         tip.length = 0.01,
                         hide.ns    = TRUE) +
      labs(title    = v,
           subtitle = paste0("Kruskal-Wallis p = ", signif(kw_p, 3))) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))

    plot_list[[v]] <- p_box
  }
}

# --- Save boxplots --- runs for both branches ---
wrap_plots(plot_list, ncol = 2) %>%
  ggsave(file.path(result_dir,
                   paste0("ASD_boxplots_k", best_K, ".pdf")),
         plot   = .,
         width  = 12,
         height = ceiling(length(plot_list) / 2) * 5,
         bg     = "white")

# --- Save classed data ---
saveRDS(ASD_lpa_input,
        file.path(PATHS[[sex_label]], asd_classed_filename))

cat("Done: ASD LPA for", sex_label, "\n")
cat("Classed data saved to:", asd_classed_filename, "\n")