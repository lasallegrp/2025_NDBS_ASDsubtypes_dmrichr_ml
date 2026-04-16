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

# --- Load inputs ---
ASD_lpa_input <- readRDS(
  file.path(PATHS[[sex_label]], "ASD_lpa_input_ADOS_ADIR.rds"))
ASD_mullen_vineland_full <- readRDS(
  file.path(PATHS[[sex_label]], "ASD_mullen_vineland_full.rds"))

cat("Running ASD LPA for", sex_label, "| n =", nrow(ASD_lpa_input), "\n")

# --- Model comparison ---
LPA_model_comparison <- ASD_lpa_input %>%
  estimate_profiles(1:LPA_CONFIG$max_K_ASD,
                    models = LPA_CONFIG$models)

fit_table <- get_fit(LPA_model_comparison)
print(fit_table)
write.csv(fit_table,
          file.path(result_dir, "ASD_model_comparison_fit_table.csv"),
          row.names = FALSE)

# --- Extract best model ---
best_K <- LPA_CONFIG$best_K_ASD
best_model_data <- get_data(LPA_model_comparison[[best_K]])
ASD_lpa_input$Class <- as.factor(best_model_data$Class)

cat("Class distribution (K =", best_K, "):\n")
print(table(ASD_lpa_input$Class))

neuro_vars <- colnames(ASD_lpa_input)[colnames(ASD_lpa_input) != "Class"]

# --- Profile plot ---
p_profiles <- plot_profiles(LPA_model_comparison[[best_K]]) +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(result_dir, paste0("ASD_profile_plot_k", best_K, ".pdf")),
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
ggsave(file.path(result_dir, paste0("ASD_line_plot_k", best_K, ".pdf")),
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
          file.path(result_dir, paste0("ASD_summary_by_class_k", best_K, ".csv")),
          row.names = FALSE)

# --- PCA plot ---
pca_result <- prcomp(ASD_lpa_input[, neuro_vars], center = TRUE, scale. = TRUE)
var_exp    <- round(summary(pca_result)$importance[2, 1:2] * 100, 1)
pca_scores <- as.data.frame(pca_result$x[, 1:2])
colnames(pca_scores) <- c("PC1", "PC2")
pca_scores$Class <- ASD_lpa_input$Class

p_pca <- ggplot(pca_scores, aes(x = PC1, y = PC2, color = Class)) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(type = "norm", linewidth = 1) +
  labs(title = paste("PCA by LPA Class — ASD", sex_label),
       x = paste0("PC1 (", var_exp[1], "%)"),
       y = paste0("PC2 (", var_exp[2], "%)")) +
  theme_minimal()
ggsave(file.path(result_dir, paste0("ASD_PCA_k", best_K, ".pdf")),
       p_pca, width = 8, height = 6)

# --- Kruskal-Wallis + Dunn ---
kw_results <- map_dfr(neuro_vars, function(v) {
  kt <- kruskal_test(ASD_lpa_input, as.formula(paste(v, "~ Class")))
  kt$Variable <- v; kt
})
write.csv(kw_results,
          file.path(result_dir, "ASD_kruskal_results.csv"), row.names = FALSE)

dunn_results <- map_dfr(neuro_vars, function(v) {
  dt <- ASD_lpa_input %>%
    dunn_test(as.formula(paste(v, "~ Class")), p.adjust.method = "BH")
  dt$Variable <- v; dt
})
write.csv(dunn_results,
          file.path(result_dir, "ASD_dunn_results.csv"), row.names = FALSE)

# --- Boxplots ---
plot_list <- list()
for (v in neuro_vars) {
  pairwise_tbl <- ASD_lpa_input %>%
    dunn_test(as.formula(paste(v, "~ Class")), p.adjust.method = "BH") %>%
    add_xy_position(x = "Class")
  kw_p <- kw_results$p[kw_results$Variable == v]

  p_box <- ggboxplot(ASD_lpa_input, x = "Class", y = v,
                     color = "Class", add = "jitter") +
    stat_pvalue_manual(pairwise_tbl, label = "p.adj.signif",
                       tip.length = 0.01, hide.ns = TRUE) +
    labs(title = v,
         subtitle = paste0("Kruskal-Wallis p = ", signif(kw_p, 3))) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  plot_list[[v]] <- p_box
}

wrap_plots(plot_list, ncol = 2) %>%
  ggsave(file.path(result_dir, paste0("ASD_boxplots_k", best_K, ".pdf")),
         plot = ., width = 12,
         height = ceiling(length(plot_list) / 2) * 5, bg = "white")

# --- Save classed data ---
saveRDS(ASD_lpa_input,
        file.path(PATHS[[sex_label]],
                  paste0("ASD_lpa_classed_k", best_K, ".rds")))

cat("Done: ASD LPA for", sex_label, "\n")