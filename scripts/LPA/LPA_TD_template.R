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
TD_lpa_input <- readRDS(
  file.path(PATHS[[sex_label]], "TD_lpa_input_Mullen_Vine_SCQ.rds"))

cat("Running TD LPA for", sex_label, "| n =", nrow(TD_lpa_input), "\n")

# --- Model comparison ---
LPA_model_comparison_TD <- TD_lpa_input %>%
  estimate_profiles(1:LPA_CONFIG$max_K_TD,
                    models = LPA_CONFIG$models)

fit_table_TD <- get_fit(LPA_model_comparison_TD)
print(fit_table_TD)
write.csv(fit_table_TD,
          file.path(result_dir, "TD_model_comparison_fit_table.csv"),
          row.names = FALSE)

# --- Extract best model ---
# Review fit_table_TD output and update LPA_CONFIG$best_K_TD
# in config.R if the female solution differs from males
best_K_TD <- LPA_CONFIG$best_K_TD

best_model_data_TD <- get_data(LPA_model_comparison_TD[[best_K_TD]])
TD_lpa_input$Class <- as.factor(best_model_data_TD$Class)

cat("Class distribution (K =", best_K_TD, "):\n")
print(table(TD_lpa_input$Class))

neuro_vars_TD <- colnames(TD_lpa_input)[colnames(TD_lpa_input) != "Class"]

# --- Profile plot ---
p_profiles <- plot_profiles(LPA_model_comparison_TD[[best_K_TD]]) +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(result_dir, paste0("TD_profile_plot_k", best_K_TD, ".pdf")),
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
ggsave(file.path(result_dir, paste0("TD_line_plot_k", best_K_TD, ".pdf")),
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
  stat_ellipse(type = "norm", linewidth = 1) +
  labs(title = paste("PCA by LPA Class — TD", sex_label),
       x = paste0("PC1 (", var_exp_TD[1], "%)"),
       y = paste0("PC2 (", var_exp_TD[2], "%)")) +
  theme_minimal()
ggsave(file.path(result_dir, paste0("TD_PCA_k", best_K_TD, ".pdf")),
       p_pca, width = 8, height = 6)

# --- T-tests + boxplots ---
# T-test used here because TD solution is 2-class
# If best_K_TD > 2, switch to Kruskal + Dunn (same pattern as ASD template)
plot_list_TD <- list()

for (v in neuro_vars_TD) {
  ttest_result <- TD_lpa_input %>%
    t_test(as.formula(paste(v, "~ Class"))) %>%
    adjust_pvalue(method = "BH") %>%
    mutate(p.adj.signif = case_when(
      p.adj <= 0.001 ~ "***",
      p.adj <= 0.01  ~ "**",
      p.adj <= 0.05  ~ "*",
      p.adj <= 0.1   ~ ".",
      TRUE           ~ "ns"
    )) %>%
    add_xy_position(x = "Class")

  p_box <- ggboxplot(TD_lpa_input, x = "Class", y = v,
                     color = "Class", add = "jitter") +
    stat_pvalue_manual(ttest_result, label = "p.adj.signif",
                       tip.length = 0.01) +
    labs(title = v,
         subtitle = paste0("BH-adjusted p = ",
                           signif(ttest_result$p.adj, 3))) +
    theme_classic(base_size = 14) +
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1))

  plot_list_TD[[v]] <- p_box
  ggsave(file.path(result_dir, paste0("TD_boxplot_", v, ".pdf")),
         p_box, width = 6, height = 5, bg = "white")
}

wrap_plots(plot_list_TD, ncol = 2) %>%
  ggsave(file.path(result_dir, paste0("TD_all_boxplots_k", best_K_TD, ".pdf")),
         plot = .,
         width  = 12,
         height = ceiling(length(plot_list_TD) / 2) * 5,
         bg     = "white")

# --- Save outputs ---
saveRDS(TD_lpa_input,
        file.path(PATHS[[sex_label]],
                  paste0("TD_lpa_classed_k", best_K_TD, ".rds")))
saveRDS(TD_refined,
        file.path(PATHS[[sex_label]], "TD_refined_comparison_group.rds"))

cat("Done: TD LPA for", sex_label, "\n")