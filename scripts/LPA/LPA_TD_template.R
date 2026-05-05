# ============================================================
# LPA_TD_template.R
# Requires: sex_label, result_dir set before sourcing
#           config.R sourced before sourcing this file
# ============================================================

library(tidyverse)
library(tidyLPA)
library(ggplot2)
library(ggpubr)
library(rstatix)
library(patchwork)

dir.create(result_dir, showWarnings = FALSE, recursive = TRUE)

# --- Load input ---
td_input_path <- file.path(PATHS[[sex_label]], "TD_lpa_input_Mullen_Vine_SCQ.rds")

if (!file.exists(td_input_path)) {
  stop(paste0(
    "\n================================================================\n",
    "ERROR: TD input file not found for ", sex_label, ":\n",
    td_input_path, "\n",
    "----------------------------------------------------------------\n",
    "Have you run 03_subset_neuromeasures.R for ", sex_label, "?\n",
    "================================================================\n"
  ))
}

TD_lpa_input <- readRDS(td_input_path)

cat("Running TD LPA for", sex_label, "| n =", nrow(TD_lpa_input), "\n")

# --- Pull sex-specific settings from nested config ---
best_K_TD      <- LPA_CONFIG$best_K_TD[[sex_label]]
best_model_num <- LPA_CONFIG$best_model_TD[[sex_label]]
k_limit_TD_m1  <- LPA_CONFIG$max_K_TD_m1[[sex_label]]
k_limit_TD_m2  <- LPA_CONFIG$max_K_TD_m2[[sex_label]]

cat("Settings for", sex_label, ":\n")
cat("  TD Model 1 max K:", k_limit_TD_m1, "\n")
cat("  TD Model 2 max K:", k_limit_TD_m2, "\n")
cat("  Best K:          ", best_K_TD, "\n")
cat("  Best model:      ", best_model_num, "\n")

# --- Model comparison --- run separately to avoid convergence failures ---
cat("Running Model 1 (equal variances) K = 1 to", k_limit_TD_m1, "\n")
LPA_TD_m1 <- TD_lpa_input %>%
  estimate_profiles(1:k_limit_TD_m1, models = 1)

cat("Running Model 2 (varying variances) K = 1 to", k_limit_TD_m2, "\n")
LPA_TD_m2 <- TD_lpa_input %>%
  estimate_profiles(1:k_limit_TD_m2, models = 2)

# --- Combine fit tables for reporting ---
fit_m1_TD <- get_fit(LPA_TD_m1)
fit_m2_TD <- get_fit(LPA_TD_m2)
fit_table_TD <- bind_rows(fit_m1_TD, fit_m2_TD) %>%
  arrange(Model, Classes)

print(fit_table_TD)
write.csv(fit_table_TD,
          file.path(result_dir, "TD_model_comparison_fit_table.csv"),
          row.names = FALSE)

# --- Extract best model safely ---
best_model_name_TD <- paste0("model_", best_model_num, "_class_", best_K_TD)
cat("Extracting solution:", best_model_name_TD, "\n")

LPA_TD_best_obj <- if (best_model_num == 1) LPA_TD_m1 else LPA_TD_m2

available_models_TD <- names(LPA_TD_best_obj)
cat("Available solutions in selected model object:\n")
print(available_models_TD)

if (!best_model_name_TD %in% available_models_TD) {
  stop(paste0(
    "\n================================================================\n",
    "ERROR: Requested solution '", best_model_name_TD, "' not found.\n",
    "----------------------------------------------------------------\n",
    "Available solutions:\n",
    paste(" -", available_models_TD, collapse = "\n"), "\n",
    "----------------------------------------------------------------\n",
    "Update best_K_TD$", sex_label, " and/or best_model_TD$",
    sex_label, " in config.R\n",
    "================================================================\n"
  ))
}

best_model_data_TD <- get_data(LPA_TD_best_obj[[best_model_name_TD]])
TD_lpa_input$Class <- as.factor(best_model_data_TD$Class)

cat("Class distribution (K =", best_K_TD, ", Model =", best_model_num, "):\n")
print(table(TD_lpa_input$Class))

neuro_vars_TD <- colnames(TD_lpa_input)[colnames(TD_lpa_input) != "Class"]
n_classes_TD  <- nlevels(TD_lpa_input$Class)

cat("Number of classes:", n_classes_TD, "\n")
cat("Statistical test:",
    if (n_classes_TD == 1) "None (K=1)" else
    if (n_classes_TD == 2) "Wilcoxon rank-sum (non-parametric, 2 groups)" else
    "Kruskal-Wallis + Dunn (non-parametric, 3+ groups)", "\n")

# --- Profile plot ---
p_profiles <- plot_profiles(LPA_TD_best_obj[[best_model_name_TD]]) +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(result_dir,
                 paste0("TD_profile_plot_k", best_K_TD, ".pdf")),
       p_profiles, width = 14, height = 6)

# --- Line plot ---
long_lpa_TD <- TD_lpa_input %>%
  pivot_longer(cols = -Class, names_to = "Variable", values_to = "Score")

p_line <- ggplot(long_lpa_TD,
                 aes(x = Variable, y = Score, group = Class, color = Class)) +
  stat_summary(fun = mean, geom = "line", linewidth = 1) +
  stat_summary(fun = mean, geom = "point") +
  theme_minimal() +
  labs(title = paste("TD Latent Profile Plot —", sex_label),
       x = "Neurodevelopmental Measure", y = "Mean Score") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(result_dir,
                 paste0("TD_line_plot_k", best_K_TD, ".pdf")),
       p_line, width = 15, height = 6)

# --- Summary table ---
summary_table_TD <- TD_lpa_input %>%
  group_by(Class) %>%
  summarise(across(all_of(neuro_vars_TD),
                   list(mean = ~mean(.x, na.rm = TRUE),
                        sd   = ~sd(.x,   na.rm = TRUE)),
                   .names = "{.col}_{.fn}"),
            .groups = "drop")
write.csv(summary_table_TD,
          file.path(result_dir,
                    paste0("TD_summary_by_class_k", best_K_TD, ".csv")),
          row.names = FALSE)

# --- Identify most quintessential TD class ---
# Select largest class as refined TD comparison group
class_sizes   <- table(TD_lpa_input$Class)
largest_class <- names(which.max(class_sizes))

cat("Largest TD class (most quintessential TD):", largest_class, "\n")
cat("Class sizes:\n"); print(class_sizes)

TD_refined <- TD_lpa_input[TD_lpa_input$Class == largest_class, ]
TD_refined$Freezerworks_ID <- rownames(TD_refined)
cat("Refined TD group n =", nrow(TD_refined), "\n")

# Save class selection rationale
class_summary <- as.data.frame(class_sizes)
colnames(class_summary) <- c("Class", "n")
class_summary$selected <- class_summary$Class == largest_class
write.csv(class_summary,
          file.path(result_dir,
                    paste0("TD_class_sizes_k", best_K_TD, ".csv")),
          row.names = FALSE)

# --- PCA plot ---
pca_result_TD <- prcomp(TD_lpa_input[, neuro_vars_TD],
                        center = TRUE, scale. = TRUE)
var_exp_TD    <- round(summary(pca_result_TD)$importance[2, 1:2] * 100, 1)
pca_scores_TD <- as.data.frame(pca_result_TD$x[, 1:2])
colnames(pca_scores_TD) <- c("PC1", "PC2")
pca_scores_TD$Class <- TD_lpa_input$Class

p_pca <- ggplot(pca_scores_TD, aes(x = PC1, y = PC2, color = Class)) +
  geom_point(size = 3, alpha = 0.8) +
  # Only add ellipse if more than 1 class and enough points
  { if (n_classes_TD > 1)
      stat_ellipse(type = "norm", linewidth = 1) } +
  labs(title = paste("PCA by LPA Class — TD", sex_label),
       x = paste0("PC1 (", var_exp_TD[1], "%)"),
       y = paste0("PC2 (", var_exp_TD[2], "%)")) +
  theme_minimal()
ggsave(file.path(result_dir,
                 paste0("TD_PCA_k", best_K_TD, ".pdf")),
       p_pca, width = 8, height = 6)

# --- Statistical tests + boxplots ---
# Adapts automatically based on number of classes:
#   K = 1: descriptive summary only
#   K = 2: Wilcoxon rank-sum (non-parametric equivalent of t-test)
#   K > 2: Kruskal-Wallis + Dunn pairwise (non-parametric multi-group)

if (n_classes_TD == 1) {

  cat("K=1 — skipping tests and boxplots (no groups to compare)\n")
  cat("Saving descriptive summary for single TD group\n")

  desc_summary <- TD_lpa_input %>%
    dplyr::select(all_of(neuro_vars_TD)) %>%
    summarise(across(everything(),
                     list(mean   = ~mean(.x,   na.rm = TRUE),
                          sd     = ~sd(.x,     na.rm = TRUE),
                          median = ~median(.x, na.rm = TRUE),
                          min    = ~min(.x,    na.rm = TRUE),
                          max    = ~max(.x,    na.rm = TRUE)),
                     .names = "{.col}_{.fn}"))
  write.csv(desc_summary,
            file.path(result_dir, "TD_descriptive_summary_k1.csv"),
            row.names = FALSE)

} else if (n_classes_TD == 2) {

  cat("K=2 — using Wilcoxon rank-sum test\n")

  wilcox_results <- map_dfr(neuro_vars_TD, function(v) {
    wt <- TD_lpa_input %>%
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
            file.path(result_dir, "TD_wilcoxon_results.csv"),
            row.names = FALSE)

  plot_list_TD <- list()
  for (v in neuro_vars_TD) {
    wilcox_tbl <- TD_lpa_input %>%
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

    p_box <- ggboxplot(TD_lpa_input, x = "Class", y = v,
                       color = "Class", add = "jitter") +
      stat_pvalue_manual(wilcox_tbl, label = "p.adj.signif",
                         tip.length = 0.01, hide.ns = TRUE) +
      labs(title    = v,
           subtitle = paste0("Wilcoxon BH-adjusted p = ",
                             signif(wilcox_p, 3))) +
      theme_classic(base_size = 14) +
      theme(panel.grid = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1))

    plot_list_TD[[v]] <- p_box
    ggsave(file.path(result_dir, paste0("TD_boxplot_", v, ".pdf")),
           p_box, width = 6, height = 5, bg = "white")
  }

  wrap_plots(plot_list_TD, ncol = 2) %>%
    ggsave(file.path(result_dir,
                     paste0("TD_all_boxplots_k", best_K_TD, ".pdf")),
           plot   = .,
           width  = 12,
           height = ceiling(length(plot_list_TD) / 2) * 5,
           bg     = "white")

} else {

  cat("K >2 — using Kruskal-Wallis + Dunn pairwise tests\n")

  kw_results_TD <- map_dfr(neuro_vars_TD, function(v) {
    kt <- kruskal_test(TD_lpa_input, as.formula(paste(v, "~ Class")))
    kt$Variable <- v
    kt
  })
  write.csv(kw_results_TD,
            file.path(result_dir, "TD_kruskal_results.csv"),
            row.names = FALSE)

  dunn_results_TD <- map_dfr(neuro_vars_TD, function(v) {
    dt <- TD_lpa_input %>%
      dunn_test(as.formula(paste(v, "~ Class")), p.adjust.method = "BH")
    dt$Variable <- v
    dt
  })
  write.csv(dunn_results_TD,
            file.path(result_dir, "TD_dunn_results.csv"),
            row.names = FALSE)

  plot_list_TD <- list()
  for (v in neuro_vars_TD) {
    pairwise_tbl <- TD_lpa_input %>%
      dunn_test(as.formula(paste(v, "~ Class")), p.adjust.method = "BH") %>%
      add_xy_position(x = "Class")
    kw_p <- kw_results_TD$p[kw_results_TD$Variable == v]

    p_box <- ggboxplot(TD_lpa_input, x = "Class", y = v,
                       color = "Class", add = "jitter") +
      stat_pvalue_manual(pairwise_tbl, label = "p.adj.signif",
                         tip.length = 0.01, hide.ns = TRUE) +
      labs(title    = v,
           subtitle = paste0("Kruskal-Wallis p = ", signif(kw_p, 3))) +
      theme_classic(base_size = 14) +
      theme(panel.grid = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1))

    plot_list_TD[[v]] <- p_box
    ggsave(file.path(result_dir, paste0("TD_boxplot_", v, ".pdf")),
           p_box, width = 6, height = 5, bg = "white")
  }

  wrap_plots(plot_list_TD, ncol = 2) %>%
    ggsave(file.path(result_dir,
                     paste0("TD_all_boxplots_k", best_K_TD, ".pdf")),
           plot   = .,
           width  = 12,
           height = ceiling(length(plot_list_TD) / 2) * 5,
           bg     = "white")
}

# --- Save outputs ---
saveRDS(TD_lpa_input,
        file.path(PATHS[[sex_label]],
                  paste0("TD_lpa_classed_k", best_K_TD, ".rds")))
saveRDS(TD_refined,
        file.path(PATHS[[sex_label]], "TD_refined_comparison_group.rds"))

cat("Done: TD LPA for", sex_label, "\n")