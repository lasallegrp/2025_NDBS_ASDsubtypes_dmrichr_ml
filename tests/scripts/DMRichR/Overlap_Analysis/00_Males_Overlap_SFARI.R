#This analysis for the overlap first ASD vs controls and next ASD subtypes

# ================================
# Set up environment and libraries
# ================================
# library(GeneOverlap)
library(openxlsx)

# library(DMRichR)
# library(GenomicRanges)
# library(GeneOverlap)
library(openxlsx)
# library(regioneR)
# library(TxDb.Hsapiens.UCSC.hg38.knownGene)
# library(org.Hs.eg.db)
library(ggvenn)
library(ggplot2)

# # ================================
# # Load annotated batch1 DMRs
# # ================================
# setwd("/quobyte/lasallegrp/George/ECHO_NBDS/DM_cyto/cytosine_reports_combined/DMR_pairwise_batches_ASD_TD/Idiopathic_ASD_subtype1_vs_Control_subtype1/V1/DMRs")
# getwd()
# batch1 <- read.xlsx("/quobyte/lasallegrp/George/ECHO_NBDS/DM_cyto/cytosine_reports_combined/DMR_pairwise_batches_ASD_TD/Idiopathic_ASD_subtype1_vs_Control_subtype1/V1/DMRs/DMRs_annotated.xlsx")
# dim(batch1)
# ================================
# Load and clean SFARI gene list
# ================================
SFARI_data <- read.csv("/quobyte/lasallegrp/projects/CHARGE/WGBS/2025_NDBS_ASDsubtypes_dmrichr_ml_George/data/metadata/raw/SFARI-Gene_genes_04-03-2025release_07-02-2025export.csv")
dim(SFARI_data)
head(SFARI_data)
# Ensure gene symbols are a character vector
SFARI_gene_symbols <- unique(as.character(SFARI_data$gene.symbol))
dim(SFARI_gene_symbols)
# # ================================
# # Load annotated universe for batch1
# # ================================
# univ1 <- read.xlsx("/quobyte/lasallegrp/George/ECHO_NBDS/DM_cyto/cytosine_reports_combined/DMR_pairwise_batches_ASD_TD/Idiopathic_ASD_subtype1_vs_Control_subtype1/V1/DMRs/background_annotated.xlsx")
# dim(univ1)

# # ================================
# # Perform GeneOverlap comparison
# # ================================
# overlap_sfari <- newGeneOverlap(
#   listA = batch1$geneSymbol,
#   listB = SFARI_gene_symbols,
#   genome.size = length(unique(univ1$geneSymbol))
# )

# overlap_sfari <- testGeneOverlap(overlap_sfari)
# overlap_sfari

# # ================================
# # Save results to Excel
# # ================================
# sfari_overlap_df <- data.frame(
#   Metric = c("Jaccard Index", "Odds Ratio", "P Value"),
#   SFARI_vs_Batch1 = c(overlap_sfari@Jaccard,
#                       overlap_sfari@odds.ratio,
#                       overlap_sfari@pval)
# )

# write.xlsx(sfari_overlap_df,
#            file = "Control1vsASD1_vs_SFARI_GeneName_Overlap_Stats.xlsx",
#            rowNames = FALSE)


# # =====================================================================================
# # Code chunk 4 - Venn Diagram
# #=====================================================================================
# #
#  # Ensure SFARI gene list is still in memory
# SFARI_gene_symbols <- unique(as.character(SFARI_gene_symbols$gene.symbol))

# # Prepare list
# venn_list_sfari <- list(
#   ASD_subtype1_vs_TD_subtype1 = unique(batch1$geneSymbol),
#   SFARI = SFARI_gene_symbols
# )

# # Plot and save
# p_sfari <- ggvenn(
#   venn_list_sfari,
#   fill_color = c("skyblue", "orchid"),
#   stroke_size = 0.5,
#   text_size = 5,
#   set_name_size = 5
# )

# ggsave(
#   filename = "Control1vsASD1_vs_SFARI_Venn.pdf",
#   plot = p_sfari,
#   width = 6,
#   height = 6
# )

# # =====================================================================================
# # Code chunk 3 - Save
# # #=====================================================================================

# # # Extract gene sets
# # genes_batch1 <- unique(batch1$geneSymbol)
# # sfari_genes <- SFARI_gene_symbols

# # # Compute overlap and exclusives
# # only_batch1_vs_sfari <- setdiff(genes_batch1, sfari_genes)
# # only_sfari_vs_batch1 <- setdiff(sfari_genes, genes_batch1)
# # overlap_sfari_batch1 <- intersect(genes_batch1, sfari_genes)

# # # Create list for export
# # sfari_venn_df <- list(
# #   Only_in_Batch1 = only_batch1_vs_sfari,
# #   Only_in_SFARI = only_sfari_vs_batch1,
# #   Overlap = overlap_sfari_batch1
# # )

# # # Normalize and save
# # max_len_sfari <- max(lengths(sfari_venn_df))
# # sfari_venn_df_equal <- lapply(sfari_venn_df, function(x) { length(x) <- max_len_sfari; return(x) })
# # sfari_venn_df_final <- as.data.frame(sfari_venn_df_equal)

# # write.xlsx(sfari_venn_df_final,
# #   file = "Control1vsASD1_vs_SFARI_GeneName_Overlap_genes.xlsx",
# #   rowNames = FALSE
# # )


# =====================================================================================
# Males
#=====================================================================================

# =====================================================================================
# Code chunk 4 - Venn Diagram  Subtype 1 vs TD
#=====================================================================================
#
setwd("/quobyte/lasallegrp/projects/CHARGE/WGBS/2025_NDBS_ASDsubtypes_dmrichr_ml_George/results/DMRichR/males/ASD_neurotypical_comparisons/Subtype1_ASD_vs_Neurotypical/DMRs")
getwd()
batch1 <- read.xlsx("/quobyte/lasallegrp/projects/CHARGE/WGBS/2025_NDBS_ASDsubtypes_dmrichr_ml_George/results/DMRichR/males/ASD_neurotypical_comparisons/Subtype1_ASD_vs_Neurotypical/DMRs/DMRs_annotated.xlsx")
dim(batch1)
univ1 <- read.xlsx("/quobyte/lasallegrp/projects/CHARGE/WGBS/2025_NDBS_ASDsubtypes_dmrichr_ml_George/results/DMRichR/males/ASD_neurotypical_comparisons/Subtype1_ASD_vs_Neurotypical/DMRs/background_annotated.xlsx")
dim(univ1)

 # Ensure SFARI gene list is still in memory
SFARI_gene_symbols <- unique(as.character(SFARI_data$gene.symbol))

# Prepare list
venn_list_sfari <- list(
  Subtype1_ASD_vs_Neurotypical = unique(batch1$geneSymbol),
  SFARI = SFARI_gene_symbols
)

# Plot and save
p_sfari <- ggvenn(
  venn_list_sfari,
  fill_color = c("skyblue", "orchid"),
  stroke_size = 0.5,
  text_size = 5,
  set_name_size = 5
)

ggsave(
  filename = "Subtype1_ASD_vs_Neurotypical_vs_SFARI_Venn.pdf",
  plot = p_sfari,
  width = 6,
  height = 6
)

# =====================================================================================
# Code chunk 5 - Venn Diagram Subtype 2 vs TD
#=====================================================================================
#
setwd("/quobyte/lasallegrp/projects/CHARGE/WGBS/2025_NDBS_ASDsubtypes_dmrichr_ml_George/results/DMRichR/males/ASD_neurotypical_comparisons/Subtype2_ASD_vs_Neurotypical/DMRs")
getwd()
batch1 <- read.xlsx("/quobyte/lasallegrp/projects/CHARGE/WGBS/2025_NDBS_ASDsubtypes_dmrichr_ml_George/results/DMRichR/males/ASD_neurotypical_comparisons/Subtype2_ASD_vs_Neurotypical/DMRs/DMRs_annotated.xlsx")
dim(batch1)
univ1 <- read.xlsx("/quobyte/lasallegrp/projects/CHARGE/WGBS/2025_NDBS_ASDsubtypes_dmrichr_ml_George/results/DMRichR/males/ASD_neurotypical_comparisons/Subtype2_ASD_vs_Neurotypical/DMRs/background_annotated.xlsx")
dim(univ1)

 # Ensure SFARI gene list is still in memory
SFARI_gene_symbols <- unique(as.character(SFARI_data$gene.symbol))

# Prepare list
venn_list_sfari <- list(
  Subtype2_ASD_vs_Neurotypical = unique(batch1$geneSymbol),
  SFARI = SFARI_gene_symbols
)

# Plot and save
p_sfari <- ggvenn(
  venn_list_sfari,
  fill_color = c("skyblue", "orchid"),
  stroke_size = 0.5,
  text_size = 5,
  set_name_size = 5
)

ggsave(
  filename = "Subtype2_ASD_vs_Neurotypical_vs_SFARI_Venn.pdf",
  plot = p_sfari,
  width = 6,
  height = 6
)

# =====================================================================================
# Code chunk 6 - Venn Diagram Subtype 3 vs TD
#=====================================================================================
#
setwd("/quobyte/lasallegrp/projects/CHARGE/WGBS/2025_NDBS_ASDsubtypes_dmrichr_ml_George/results/DMRichR/males/ASD_neurotypical_comparisons/Subtype3_ASD_vs_Neurotypical/DMRs")
getwd()
batch1 <- read.xlsx("/quobyte/lasallegrp/projects/CHARGE/WGBS/2025_NDBS_ASDsubtypes_dmrichr_ml_George/results/DMRichR/males/ASD_neurotypical_comparisons/Subtype3_ASD_vs_Neurotypical/DMRs/DMRs_annotated.xlsx")
dim(batch1)
univ1 <- read.xlsx("/quobyte/lasallegrp/projects/CHARGE/WGBS/2025_NDBS_ASDsubtypes_dmrichr_ml_George/results/DMRichR/males/ASD_neurotypical_comparisons/Subtype3_ASD_vs_Neurotypical/DMRs/background_annotated.xlsx")
dim(univ1)

 # Ensure SFARI gene list is still in memory
SFARI_gene_symbols <- unique(as.character(SFARI_data$gene.symbol))

# Prepare list
venn_list_sfari <- list(
  Subtype3_ASD_vs_Neurotypical = unique(batch1$geneSymbol),
  SFARI = SFARI_gene_symbols
)

# Plot and save
p_sfari <- ggvenn(
  venn_list_sfari,
  fill_color = c("skyblue", "orchid"),
  stroke_size = 0.5,
  text_size = 5,
  set_name_size = 5
)

ggsave(
  filename = "Subtype3_ASD_vs_Neurotypical_vs_SFARI_Venn.pdf",
  plot = p_sfari,
  width = 6,
  height = 6
)



# =====================================================================================
# Code chunk 7 - Venn Diagram  Subtype 1 vs Subtype 2
#=====================================================================================
#
setwd("/quobyte/lasallegrp/projects/CHARGE/WGBS/2025_NDBS_ASDsubtypes_dmrichr_ml_George/results/DMRichR/males/ASD_comparisons/Subtype1_ASD_vs_Subtype2_ASD/DMRs")
getwd()
batch1 <- read.xlsx("/quobyte/lasallegrp/projects/CHARGE/WGBS/2025_NDBS_ASDsubtypes_dmrichr_ml_George/results/DMRichR/males/ASD_comparisons/Subtype1_ASD_vs_Subtype2_ASD/DMRs/DMRs_annotated.xlsx")
dim(batch1)
univ1 <- read.xlsx("/quobyte/lasallegrp/projects/CHARGE/WGBS/2025_NDBS_ASDsubtypes_dmrichr_ml_George/results/DMRichR/males/ASD_comparisons/Subtype1_ASD_vs_Subtype2_ASD/DMRs/background_annotated.xlsx")
dim(univ1)

 # Ensure SFARI gene list is still in memory
SFARI_gene_symbols <- unique(as.character(SFARI_data$gene.symbol))

# Prepare list
venn_list_sfari <- list(
  Subtype1_ASD_vs_Subtype2_ASD = unique(batch1$geneSymbol),
  SFARI = SFARI_gene_symbols
)

# Plot and save
p_sfari <- ggvenn(
  venn_list_sfari,
  fill_color = c("skyblue", "orchid"),
  stroke_size = 0.5,
  text_size = 5,
  set_name_size = 5
)

ggsave(
  filename = "Subtype1_ASD_vs_Subtype2_ASD_vs_SFARI_Venn.pdf",
  plot = p_sfari,
  width = 6,
  height = 6
)

# =====================================================================================
# Code chunk 8 - Venn Diagram Subtype 1 vs Subtype 3
#=====================================================================================
#
setwd("/quobyte/lasallegrp/projects/CHARGE/WGBS/2025_NDBS_ASDsubtypes_dmrichr_ml_George/results/DMRichR/males/ASD_comparisons/Subtype1_ASD_vs_Subtype3_ASD/DMRs")
getwd()
batch1 <- read.xlsx("/quobyte/lasallegrp/projects/CHARGE/WGBS/2025_NDBS_ASDsubtypes_dmrichr_ml_George/results/DMRichR/males/ASD_comparisons/Subtype1_ASD_vs_Subtype3_ASD/DMRs/DMRs_annotated.xlsx")
dim(batch1)
univ1 <- read.xlsx("/quobyte/lasallegrp/projects/CHARGE/WGBS/2025_NDBS_ASDsubtypes_dmrichr_ml_George/results/DMRichR/males/ASD_comparisons/Subtype1_ASD_vs_Subtype3_ASD/DMRs/background_annotated.xlsx")
dim(univ1)

 # Ensure SFARI gene list is still in memory
SFARI_gene_symbols <- unique(as.character(SFARI_data$gene.symbol))

# Prepare list
venn_list_sfari <- list(
  Subtype1_ASD_vs_Subtype3_ASD = unique(batch1$geneSymbol),
  SFARI = SFARI_gene_symbols
)

# Plot and save
p_sfari <- ggvenn(
  venn_list_sfari,
  fill_color = c("skyblue", "orchid"),
  stroke_size = 0.5,
  text_size = 5,
  set_name_size = 5
)

ggsave(
  filename = "Subtype1_ASD_vs_Subtype3_ASD_vs_SFARI_Venn.pdf",
  plot = p_sfari,
  width = 6,
  height = 6
)

# =====================================================================================
# Code chunk 9 - Venn Diagram Subtype 2 vs Subtype 3
#=====================================================================================
#
setwd("/quobyte/lasallegrp/projects/CHARGE/WGBS/2025_NDBS_ASDsubtypes_dmrichr_ml_George/results/DMRichR/males/ASD_comparisons/Subtype2_ASD_vs_Subtype3_ASD/DMRs")
getwd()
batch1 <- read.xlsx("/quobyte/lasallegrp/projects/CHARGE/WGBS/2025_NDBS_ASDsubtypes_dmrichr_ml_George/results/DMRichR/males/ASD_comparisons/Subtype2_ASD_vs_Subtype3_ASD/DMRs/DMRs_annotated.xlsx")
dim(batch1)
univ1 <- read.xlsx("/quobyte/lasallegrp/projects/CHARGE/WGBS/2025_NDBS_ASDsubtypes_dmrichr_ml_George/results/DMRichR/males/ASD_comparisons/Subtype2_ASD_vs_Subtype3_ASD/DMRs/background_annotated.xlsx")
dim(univ1)

 # Ensure SFARI gene list is still in memory
SFARI_gene_symbols <- unique(as.character(SFARI_data$gene.symbol))

# Prepare list
venn_list_sfari <- list(
  Subtype2_ASD_vs_Subtype3_ASD = unique(batch1$geneSymbol),
  SFARI = SFARI_gene_symbols
)

# Plot and save
p_sfari <- ggvenn(
  venn_list_sfari,
  fill_color = c("skyblue", "orchid"),
  stroke_size = 0.5,
  text_size = 5,
  set_name_size = 5
)

ggsave(
  filename = "Subtype2_ASD_vs_Subtype3_ASD_vs_SFARI_Venn.pdf",
  plot = p_sfari,
  width = 6,
  height = 6
)

# =====================================================================================
# Females
#=====================================================================================

# =====================================================================================
# Code chunk 10 - Venn Diagram  Subtype 1 vs TD
#=====================================================================================
#
setwd("/quobyte/lasallegrp/projects/CHARGE/WGBS/2025_NDBS_ASDsubtypes_dmrichr_ml_George/results/DMRichR/females/ASD_neurotypical_comparisons/Subtype1_ASD_vs_Neurotypical/DMRs")
getwd()
batch1 <- read.xlsx("/quobyte/lasallegrp/projects/CHARGE/WGBS/2025_NDBS_ASDsubtypes_dmrichr_ml_George/results/DMRichR/females/ASD_neurotypical_comparisons/Subtype1_ASD_vs_Neurotypical/DMRs/DMRs_annotated.xlsx")
dim(batch1)
univ1 <- read.xlsx("/quobyte/lasallegrp/projects/CHARGE/WGBS/2025_NDBS_ASDsubtypes_dmrichr_ml_George/results/DMRichR/females/ASD_neurotypical_comparisons/Subtype1_ASD_vs_Neurotypical/DMRs/background_annotated.xlsx")
dim(univ1)

 # Ensure SFARI gene list is still in memory
SFARI_gene_symbols <- unique(as.character(SFARI_data$gene.symbol))

# Prepare list
venn_list_sfari <- list(
  Subtype1_ASD_vs_Neurotypical = unique(batch1$geneSymbol),
  SFARI = SFARI_gene_symbols
)

# Plot and save
p_sfari <- ggvenn(
  venn_list_sfari,
  fill_color = c("skyblue", "orchid"),
  stroke_size = 0.5,
  text_size = 5,
  set_name_size = 5
)

ggsave(
  filename = "Subtype1_ASD_vs_Neurotypical_vs_SFARI_Venn.pdf",
  plot = p_sfari,
  width = 6,
  height = 6
)

# =====================================================================================
# Code chunk 11 - Venn Diagram Subtype 2 vs TD
#=====================================================================================
#
setwd("/quobyte/lasallegrp/projects/CHARGE/WGBS/2025_NDBS_ASDsubtypes_dmrichr_ml_George/results/DMRichR/females/ASD_neurotypical_comparisons/Subtype2_ASD_vs_Neurotypical/DMRs")
getwd()
batch1 <- read.xlsx("/quobyte/lasallegrp/projects/CHARGE/WGBS/2025_NDBS_ASDsubtypes_dmrichr_ml_George/results/DMRichR/females/ASD_neurotypical_comparisons/Subtype2_ASD_vs_Neurotypical/DMRs/DMRs_annotated.xlsx")
dim(batch1)
univ1 <- read.xlsx("/quobyte/lasallegrp/projects/CHARGE/WGBS/2025_NDBS_ASDsubtypes_dmrichr_ml_George/results/DMRichR/females/ASD_neurotypical_comparisons/Subtype2_ASD_vs_Neurotypical/DMRs/background_annotated.xlsx")
dim(univ1)

 # Ensure SFARI gene list is still in memory
SFARI_gene_symbols <- unique(as.character(SFARI_data$gene.symbol))

# Prepare list
venn_list_sfari <- list(
  Subtype2_ASD_vs_Neurotypical = unique(batch1$geneSymbol),
  SFARI = SFARI_gene_symbols
)

# Plot and save
p_sfari <- ggvenn(
  venn_list_sfari,
  fill_color = c("skyblue", "orchid"),
  stroke_size = 0.5,
  text_size = 5,
  set_name_size = 5
)

ggsave(
  filename = "Subtype2_ASD_vs_Neurotypical_vs_SFARI_Venn.pdf",
  plot = p_sfari,
  width = 6,
  height = 6
)

# =====================================================================================
# Code chunk 12 - Venn Diagram  Subtype 1 vs Subtype 2
#=====================================================================================
#
setwd("/quobyte/lasallegrp/projects/CHARGE/WGBS/2025_NDBS_ASDsubtypes_dmrichr_ml_George/results/DMRichR/females/ASD_comparisons/Subtype1_ASD_vs_Subtype2_ASD/DMRs")
getwd()
batch1 <- read.xlsx("/quobyte/lasallegrp/projects/CHARGE/WGBS/2025_NDBS_ASDsubtypes_dmrichr_ml_George/results/DMRichR/females/ASD_comparisons/Subtype1_ASD_vs_Subtype2_ASD/DMRs/DMRs_annotated.xlsx")
dim(batch1)
univ1 <- read.xlsx("/quobyte/lasallegrp/projects/CHARGE/WGBS/2025_NDBS_ASDsubtypes_dmrichr_ml_George/results/DMRichR/females/ASD_comparisons/Subtype1_ASD_vs_Subtype2_ASD/DMRs/background_annotated.xlsx")
dim(univ1)

 # Ensure SFARI gene list is still in memory
SFARI_gene_symbols <- unique(as.character(SFARI_data$gene.symbol))

# Prepare list
venn_list_sfari <- list(
  Subtype1_ASD_vs_Subtype2_ASD = unique(batch1$geneSymbol),
  SFARI = SFARI_gene_symbols
)

# Plot and save
p_sfari <- ggvenn(
  venn_list_sfari,
  fill_color = c("skyblue", "orchid"),
  stroke_size = 0.5,
  text_size = 5,
  set_name_size = 5
)

ggsave(
  filename = "Subtype1_ASD_vs_Subtype2_ASD_vs_SFARI_Venn.pdf",
  plot = p_sfari,
  width = 6,
  height = 6
)
# #
# setwd("/quobyte/lasallegrp/George/ECHO_NBDS/DM_cyto/cytosine_reports_combined/DMR_pairwise_batches_ASD_subtypes/Idiopathic_ASD_subtype3_vs_Idiopathic_ASD_subtype2/V3/DMRs")
# getwd()
# batch1 <- read.xlsx("/quobyte/lasallegrp/George/ECHO_NBDS/DM_cyto/cytosine_reports_combined/DMR_pairwise_batches_ASD_subtypes/Idiopathic_ASD_subtype3_vs_Idiopathic_ASD_subtype2/V3/DMRs/DMRs_annotated.xlsx")
# dim(batch1)
# univ1 <- read.xlsx("/quobyte/lasallegrp/George/ECHO_NBDS/DM_cyto/cytosine_reports_combined/DMR_pairwise_batches_ASD_subtypes/Idiopathic_ASD_subtype3_vs_Idiopathic_ASD_subtype2/V3/DMRs/background_annotated.xlsx")
# dim(univ1)

#  # Ensure SFARI gene list is still in memory
# SFARI_gene_symbols <- unique(as.character(SFARI_data$gene.symbol))

# # Prepare list
# venn_list_sfari <- list(
#   ASD_subtype3_vs_ASD_subtype2 = unique(batch1$geneSymbol),
#   SFARI = SFARI_gene_symbols
# )

# # Plot and save
# p_sfari <- ggvenn(
#   venn_list_sfari,
#   fill_color = c("skyblue", "orchid"),
#   stroke_size = 0.5,
#   text_size = 5,
#   set_name_size = 5
# )

# ggsave(
#   filename = "ASD_subtype3_vs_ASD_subtype2_vs_SFARI_Venn.pdf",
#   plot = p_sfari,
#   width = 6,
#   height = 6
# )


# #
# setwd("/quobyte/lasallegrp/George/ECHO_NBDS/DM_cyto/cytosine_reports_combined/DMR_pairwise_batches_ASD_subtypes/Idiopathic_ASD_subtype4_vs_Idiopathic_ASD_subtype1/V3/DMRs")
# getwd()
# batch1 <- read.xlsx("/quobyte/lasallegrp/George/ECHO_NBDS/DM_cyto/cytosine_reports_combined/DMR_pairwise_batches_ASD_subtypes/Idiopathic_ASD_subtype4_vs_Idiopathic_ASD_subtype1/V3/DMRs/DMRs_annotated.xlsx")
# dim(batch1)
# univ1 <- read.xlsx("/quobyte/lasallegrp/George/ECHO_NBDS/DM_cyto/cytosine_reports_combined/DMR_pairwise_batches_ASD_subtypes/Idiopathic_ASD_subtype4_vs_Idiopathic_ASD_subtype1/V3/DMRs/background_annotated.xlsx")
# dim(univ1)

#  # Ensure SFARI gene list is still in memory
# SFARI_gene_symbols <- unique(as.character(SFARI_data$gene.symbol))

# # Prepare list
# venn_list_sfari <- list(
#   ASD_subtype4_vs_ASD_subtype1 = unique(batch1$geneSymbol),
#   SFARI = SFARI_gene_symbols
# )

# # Plot and save
# p_sfari <- ggvenn(
#   venn_list_sfari,
#   fill_color = c("skyblue", "orchid"),
#   stroke_size = 0.5,
#   text_size = 5,
#   set_name_size = 5
# )

# ggsave(
#   filename = "ASD_subtype4_vs_ASD_subtype1_vs_SFARI_Venn.pdf",
#   plot = p_sfari,
#   width = 6,
#   height = 6
# )

# #
# setwd("/quobyte/lasallegrp/George/ECHO_NBDS/DM_cyto/cytosine_reports_combined/DMR_pairwise_batches_ASD_subtypes/Idiopathic_ASD_subtype4_vs_Idiopathic_ASD_subtype2/V3/DMRs")
# getwd()
# batch1 <- read.xlsx("/quobyte/lasallegrp/George/ECHO_NBDS/DM_cyto/cytosine_reports_combined/DMR_pairwise_batches_ASD_subtypes/Idiopathic_ASD_subtype4_vs_Idiopathic_ASD_subtype2/V3/DMRs/DMRs_annotated.xlsx")
# dim(batch1)
# univ1 <- read.xlsx("/quobyte/lasallegrp/George/ECHO_NBDS/DM_cyto/cytosine_reports_combined/DMR_pairwise_batches_ASD_subtypes/Idiopathic_ASD_subtype4_vs_Idiopathic_ASD_subtype2/V3/DMRs/background_annotated.xlsx")
# dim(univ1)

#  # Ensure SFARI gene list is still in memory
# SFARI_gene_symbols <- unique(as.character(SFARI_data$gene.symbol))

# # Prepare list
# venn_list_sfari <- list(
#   ASD_subtype4_vs_ASD_subtype2 = unique(batch1$geneSymbol),
#   SFARI = SFARI_gene_symbols
# )

# # Plot and save
# p_sfari <- ggvenn(
#   venn_list_sfari,
#   fill_color = c("skyblue", "orchid"),
#   stroke_size = 0.5,
#   text_size = 5,
#   set_name_size = 5
# )

# ggsave(
#   filename = "ASD_subtype4_vs_ASD_subtype2_vs_SFARI_Venn.pdf",
#   plot = p_sfari,
#   width = 6,
#   height = 6
# )

