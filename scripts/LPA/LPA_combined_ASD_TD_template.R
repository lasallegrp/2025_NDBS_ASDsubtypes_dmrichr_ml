# ============================================================
# LPA_combined_ASD_TD_template.R
# Requires: sex_label, result_dir set before sourcing
#           config.R sourced before sourcing this file
# ============================================================

library(tidyverse)
library(ggplot2)
library(ggpubr)
library(rstatix)
library(FSA)
library(patchwork)

dir.create(result_dir, showWarnings = FALSE, recursive = TRUE)

# --- Load classed ASD and TD data ---
# Uses LPA_CONFIG$asd_classed_file so filename stays in sync with config
ASD_classed <- readRDS(
  file.path(PATHS[[sex_label]], LPA_CONFIG$asd_classed_file))

TD_classed <- readRDS(
  file.path(PATHS[[sex_label]],
            paste0("TD_lpa_classed_k", LPA_CONFIG$best_K_TD, ".rds")))

ASD_mullen_vineland_full <- readRDS(
  file.path(PATHS[[sex_label]], "ASD_mullen_vineland_full.rds"))

cat("Combining ASD (n =", nrow(ASD_classed),
    ") and TD (n =", nrow(TD_classed), ") for", sex_label, "\n")

# --- Recode TD class labels to avoid overlap with ASD classes ---
# ASD classes: 1, 2, 3  →  TD classes become 4, 5
# This makes all groups distinguishable in combined plots
n_asd_classes <- nlevels(ASD_classed$Class)

TD_classed$Class <- factor(
  TD_classed$Class,
  levels = levels(TD_classed$Class),
  labels = as.character(
    seq(n_asd_classes + 1,
        n_asd_classes + nlevels(TD_classed$Class))
  )
)

cat("ASD class distribution:\n");  print(table(ASD_classed$Class))
cat("TD class distribution (relabelled):\n"); print(table(TD_classed$Class))

# --- Identify shared neuromeasure columns automatically ---
# intersect() finds all common columns across both datasets
# We then exclude Class since it is metadata not a neuromeasure
common_cols       <- intersect(colnames(ASD_classed), colnames(TD_classed))
shared_neuro_vars <- common_cols[common_cols != "Class"]

cat("ASD columns:  ", paste(colnames(ASD_classed),  collapse = ", "), "\n")
cat("TD  columns:  ", paste(colnames(TD_classed),   collapse = ", "), "\n")
cat("Shared neuromeasure columns detected:\n")
print(shared_neuro_vars)

# --- Validate the overlap is meaningful ---
if (length(shared_neuro_vars) == 0) {
  stop(paste0(
    "\n================================================================\n",
    "ERROR: No shared neuromeasure columns found between ASD and TD.\n",
    "----------------------------------------------------------------\n",
    "ASD columns: ", paste(colnames(ASD_classed), collapse = ", "), "\n",
    "TD  columns: ", paste(colnames(TD_classed),  collapse = ", "), "\n",
    "----------------------------------------------------------------\n",
    "Check that both classed .rds files contain overlapping measures.\n",
    "ASD classed file: ", LPA_CONFIG$asd_classed_file, "\n",
    "TD  classed file: TD_lpa_classed_k", LPA_CONFIG$best_K_TD, ".rds\n",
    "================================================================\n"
  ))
}

if (length(shared_neuro_vars) < 2) {
  warning(paste0(
    "\nWARNING: Only ", length(shared_neuro_vars),
    " shared neuromeasure(s) found: ",
    paste(shared_neuro_vars, collapse = ", "), "\n",
    "Combined plots will have limited interpretability.\n",
    "Consider whether ASD and TD input files have sufficient overlap.\n"
  ))
}

cat("Proceeding with", length(shared_neuro_vars),
    "shared neuromeasure(s) for combined analysis.\n")

# --- Combine using shared columns ---
ASD_subset_combined <- ASD_classed[, c(shared_neuro_vars, "Class")]
TD_subset_combined  <- TD_classed[,  c(shared_neuro_vars, "Class")]

combined_data <- rbind(ASD_subset_combined, TD_subset_combined)
combined_data$Group <- c(
  rep("ASD", nrow(ASD_classed)),
  rep("TD",  nrow(TD_classed))
)

cat("Combined dataset: n =", nrow(combined_data), "\n")

plot_vars <- shared_neuro_vars

# --- Combined profile line plot ---
long_combined <- combined_data %>%
  pivot_longer(cols      = all_of(plot_vars),
               names_to  = "Variable",
               values_to = "Score")

p_combined <- ggplot(long_combined,
                     aes(x     = Variable,
                         y     = Score,
                         group = Class,
                         color = Class)) +
  stat_summary(fun = mean, geom = "line",  linewidth = 1) +
  stat_summary(fun = mean, geom = "point", size = 2) +
  theme_minimal(base_size = 14) +
  labs(title = paste("ASD + TD Combined Latent Profiles —", sex_label),
       x     = "Neurodevelopmental Measure",
       y     = "Mean Score") +
  theme(axis.text.x = element_text(angle = 50, hjust = 1),
        plot.margin = margin(t = 10, r = 10, b = 10, l = 20))

ggsave(file.path(result_dir, "combined_ASD_TD_line_plot.pdf"),
       p_combined, width = 18, height = 6, dpi = 300)

# --- Kruskal-Wallis across all classes ---
kw_results_combined <- map_dfr(plot_vars, function(v) {
  kt <- kruskal_test(combined_data, as.formula(paste(v, "~ Class")))
  kt$Variable <- v
  kt
})
write.csv(kw_results_combined,
          file.path(result_dir, "combined_kruskal_results.csv"),
          row.names = FALSE)

# --- Dunn pairwise tests ---
dunn_results_combined <- map_dfr(plot_vars, function(v) {
  dt <- combined_data %>%
    dunn_test(as.formula(paste(v, "~ Class")), p.adjust.method = "BH")
  dt$Variable <- v
  dt
})
write.csv(dunn_results_combined,
          file.path(result_dir, "combined_dunn_results.csv"),
          row.names = FALSE)

# --- Boxplots ---
plot_list_combined <- list()

for (v in plot_vars) {
  kw_p <- kw_results_combined$p[kw_results_combined$Variable == v]

  pairwise_tbl <- combined_data %>%
    dunn_test(as.formula(paste(v, "~ Class")), p.adjust.method = "BH") %>%
    add_xy_position(x = "Class")

  p_box <- ggboxplot(combined_data, x = "Class", y = v,
                     color = "Class", add = "jitter") +
    stat_pvalue_manual(pairwise_tbl,
                       label      = "p.adj.signif",
                       tip.length = 0.01,
                       hide.ns    = TRUE) +
    labs(title    = v,
         subtitle = paste0("Kruskal-Wallis p = ", signif(kw_p, 3))) +
    theme_classic(base_size = 13) +
    theme(panel.grid  = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1))

  plot_list_combined[[v]] <- p_box
  ggsave(file.path(result_dir, paste0("boxplot_", v, ".pdf")),
         p_box, width = 6, height = 5, bg = "white")
}

wrap_plots(plot_list_combined, ncol = 2) %>%
  ggsave(file.path(result_dir, "all_combined_boxplots.pdf"),
         plot   = .,
         width  = 12,
         height = ceiling(length(plot_list_combined) / 2) * 5,
         bg     = "white")

# --- Save combined classed data ---
saveRDS(combined_data,
        file.path(PATHS[[sex_label]], "ASD_TD_combined_classed.rds"))

# --- Save group membership summary ---
group_summary <- combined_data %>%
  group_by(Group, Class) %>%
  summarise(n = n(), .groups = "drop")

write.csv(group_summary,
          file.path(result_dir, "combined_group_class_summary.csv"),
          row.names = FALSE)

cat("Done: Combined ASD + TD LPA for", sex_label, "\n")