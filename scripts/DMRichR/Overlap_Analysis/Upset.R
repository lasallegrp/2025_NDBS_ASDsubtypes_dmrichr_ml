srun --job-name=meth      --account=genome-center-grp      --partition=high      --cpus-per-task=10      --mem=200G      --time=8:00:00      --pty bash

source ~/.profile
conda activate /quobyte/lasallegrp/programs/.conda/envs/DMRichR4.2_Updated/

#you can also just run the R script directly which is called upset_runner.R
#Rscript upset_runner.R

# activate your conda env
source ~/.profile
conda activate /quobyte/lasallegrp/programs/.conda/envs/DMRichR4.2_Updated/

# change into a comparison folder, then run: #e.g
cd /quobyte/lasallegrp/George/ECHO_NBDS/.../Idiopathic_ASD_subtype2_vs_Control_subtype1

# Option A: look for all V1/V2/V3... under the parent folder
Rscript upset_runner.R             # root=.. (default)

# Option B: only look under the current folder
# old version Rscript /quobyte/lasallegrp/George/ECHO_NBDS/DM_cyto/cytosine_reports_combined/Overlap_Analysis/Males_Overlap/upset_runner.R --root . --out upset_outputs --min_size 1 --width_ratio auto --wrap 36 --dpi 300

Rscript /quobyte/lasallegrp/George/ECHO_NBDS/DM_cyto/cytosine_reports_combined/Overlap_Analysis/Males_Overlap/upset_runner.R --root . --use_background_universe true --min_size 1 --width_ratio auto --wrap 36 --dpi 300


setwd("/quobyte/lasallegrp/George/ECHO_NBDS/DM_cyto/cytosine_reports_combined/Overlap_Analysis/Males_Overlap")

# ---- Setup ----
library(readxl)
library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)
library(ComplexUpset)
library(readr)
# ---- File paths & labels ----
files <- c(
  Idiopathic_ASD_subtype1_vs_Control_subtype1_71 = "/quobyte/lasallegrp/George/ECHO_NBDS/DM_cyto/cytosine_reports_combined/DMR_pairwise_batches_ASD_TD/Idiopathic_ASD_subtype1_vs_Control_subtype1/V3_cutoff_003_maxPerms_71/DMRs/DMRs_annotated.xlsx", 
  Idiopathic_ASD_subtype2_vs_Control_subtype1_64   = "/quobyte/lasallegrp/George/ECHO_NBDS/DM_cyto/cytosine_reports_combined/DMR_pairwise_batches_ASD_TD/Idiopathic_ASD_subtype2_vs_Control_subtype1/V3_cutoff_003_maxPerms_64/DMRs/DMRs_annotated.xlsx",
  Idiopathic_ASD_subtype3_vs_Control_subtype1_25 = "/quobyte/lasallegrp/George/ECHO_NBDS/DM_cyto/cytosine_reports_combined/DMR_pairwise_batches_ASD_TD/Idiopathic_ASD_subtype3_vs_Control_subtype1/V2_cutoff_003_maxPerms_25/DMRs/DMRs_annotated.xlsx",
  Idiopathic_ASD_subtype4_vs_Control_subtype1_82 = "/quobyte/lasallegrp/George/ECHO_NBDS/DM_cyto/cytosine_reports_combined/DMR_pairwise_batches_ASD_TD/Idiopathic_ASD_subtype4_vs_Control_subtype1/V4_cutoff_003_maxPerms_82_perGroup_1/DMRs/DMRs_annotated.xlsx",

  Idiopathic_ASD_subtype1_vs_Control_subtype2_25   = "/quobyte/lasallegrp/George/ECHO_NBDS/DM_cyto/cytosine_reports_combined/Males/DMR_pairwise_batches_ASD_TD2/Idiopathic_ASD_subtype1_vs_Control_subtype2/V2_cutoff_003_maxPerms_25/DMRs/DMRs_annotated.xlsx",
  Idiopathic_ASD_subtype2_vs_Control_subtype2_25   = "/quobyte/lasallegrp/George/ECHO_NBDS/DM_cyto/cytosine_reports_combined/Males/DMR_pairwise_batches_ASD_TD2/Idiopathic_ASD_subtype2_vs_Control_subtype2/V2_cutoff_003_maxPerms_25/DMRs/DMRs_annotated.xlsx",
  Idiopathic_ASD_subtype3_vs_Control_subtype2_25  = "/quobyte/lasallegrp/George/ECHO_NBDS/DM_cyto/cytosine_reports_combined/Males/DMR_pairwise_batches_ASD_TD2/Idiopathic_ASD_subtype3_vs_Control_subtype2/V2_cutoff_003_maxPerms_25/DMRs/DMRs_annotated.xlsx",
  Idiopathic_ASD_subtype4_vs_Control_subtype2_25  = "/quobyte/lasallegrp/George/ECHO_NBDS/DM_cyto/cytosine_reports_combined/Males/DMR_pairwise_batches_ASD_TD2/Idiopathic_ASD_subtype4_vs_Control_subtype2/V2_cutoff_003_maxPerms_25/DMRs/DMRs_annotated.xlsx"
  
)

dir.create("upset_outputs", showWarnings = FALSE)

# ---- Helper: read one file and return gene symbols + direction ----
read_dmrs_genes <- function(path) {
  x <- read_excel(path, guess_max = 1e6)

  # extract genes
  genes <- x$geneSymbol
  genes <- genes[!is.na(genes)]
  genes <- unlist(strsplit(genes, "[;,/\\|]"))  # split multi-gene cells
  genes <- str_trim(genes)
  genes <- genes[nzchar(genes)]
  genes <- unique(genes)

  # direction
  dir_vec <- NULL
  if ("direction" %in% names(x)) {
    df <- x %>%
      filter(!is.na(geneSymbol)) %>%
      mutate(gene = strsplit(geneSymbol, "[;,/\\|]")) %>%
      tidyr::unnest_longer(gene) %>%
      mutate(gene = str_trim(gene)) %>%
      filter(gene != "") %>%
      select(gene, direction)

    dir_map <- df %>%
      group_by(gene) %>%
      summarize(direction = {
        d <- unique(na.omit(direction))
        if (!length(d)) NA_character_
        else if (length(d) == 1) d
        else "Mixed"
      }, .groups = "drop")

    dir_vec <- setNames(dir_map$direction, dir_map$gene)
  }

  tibble(gene = genes,
         direction = if (is.null(dir_vec)) NA_character_ else unname(dir_vec[genes]))
}

# ---- Load all versions ----
gene_list <- lapply(files, read_dmrs_genes)

# ---- Build membership ----
all_genes <- sort(unique(unlist(lapply(gene_list, `[[`, "gene"))))
membership <- tibble(gene = all_genes)
for (lab in names(gene_list)) {
  membership[[lab]] <- membership$gene %in% gene_list[[lab]]$gene
}
write.csv(membership, "upset_outputs/GENE_membership_all.csv", row.names = FALSE)

# --- 1) Read & clean SFARI gene symbols ---
sfari_path <- "/quobyte/lasallegrp/George/ECHO_NBDS/DM_cyto/cytosine_reports_combined/SFARI-Gene_genes_04-03-2025release_07-02-2025export.csv"

SFARI_Genes_raw <- read_csv(sfari_path, show_col_types = FALSE)

# Be robust to column naming variants
possible_cols <- c("gene.symbol","gene-symbol","Symbol","symbol","Gene Symbol","GeneSymbol","Gene")
gcol <- possible_cols[possible_cols %in% names(SFARI_Genes_raw)][1]
if (is.na(gcol)) stop("Couldn't find a gene symbol column in the SFARI CSV. Found cols: ", paste(names(SFARI_Genes_raw), collapse=", "))

SFARI_gene_symbols <- SFARI_Genes_raw[[gcol]] %>%
  as.character() %>%
  .[!is.na(.)] %>%
  str_split("[;,/\\|]") %>%                     # split if any multi-gene cells exist
  unlist() %>%
  str_trim() %>%
  .[nzchar(.)] %>%
  unique()

length(SFARI_gene_symbols)  # quick sanity check

# --- 2) Append to your membership (gene-based) ---
# membership must already exist from earlier (first column = gene)
if (!exists("membership") || !"gene" %in% names(membership)) {
  stop("membership (gene-based) not found. Run the earlier code that builds `membership` first.")
}

# Add SFARI column
membership <- membership %>%
  mutate(SFARI = gene %in% SFARI_gene_symbols)

# Optional: save updated matrix
write.csv(membership, "upset_outputs/GENE_membership_all_with_SFARI.csv", row.names = FALSE)

# --- 3) Plot with SFARI included ---
# Use membership columns (all set columns, including SFARI), not names(files)

# sets = all logical columns except the 'gene' column
sets <- setdiff(names(membership), "gene")

# # Use namespaced layer constructors (so we don't hit your helper)
# ba <- list(
#   "Intersection size" = ComplexUpset::intersection_size(),
#   "Set size"          = ComplexUpset::upset_set_size()
# )

p_all <- ComplexUpset::upset(
  membership,
  intersect = sets,
  name = "Intersection size (genes)",
  width_ratio = 0.33,
  sort_intersections = "descending",
  min_size = 3
  # base_annotations = ba
)

ggsave("upset_outputs/UpSet_all_GENES_with_SFARI.pdf", p_all, width = 30, height = 15)
# ggsave("upset_outputs/UpSet_all_GENES_with_SFARI.png", p_all, width = 20, height = 6, dpi = 300)



#This worked,I like it more because it shows the correct setsize for the SFARI.1230
setwd("/quobyte/lasallegrp/George/ECHO_NBDS/DM_cyto/cytosine_reports_combined/Overlap_Analysis/Males_Overlap")

# ---- Setup ----
library(readxl)
library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)
library(ComplexUpset)
library(readr)
# ---- File paths & labels ----
files <- c(
  Idiopathic_ASD_subtype1_vs_Control_subtype1_71 = "/quobyte/lasallegrp/George/ECHO_NBDS/DM_cyto/cytosine_reports_combined/DMR_pairwise_batches_ASD_TD/Idiopathic_ASD_subtype1_vs_Control_subtype1/V3_cutoff_003_maxPerms_71/DMRs/DMRs_annotated.xlsx", 
  Idiopathic_ASD_subtype2_vs_Control_subtype1_64   = "/quobyte/lasallegrp/George/ECHO_NBDS/DM_cyto/cytosine_reports_combined/DMR_pairwise_batches_ASD_TD/Idiopathic_ASD_subtype2_vs_Control_subtype1/V3_cutoff_003_maxPerms_64/DMRs/DMRs_annotated.xlsx",
  Idiopathic_ASD_subtype3_vs_Control_subtype1_25 = "/quobyte/lasallegrp/George/ECHO_NBDS/DM_cyto/cytosine_reports_combined/DMR_pairwise_batches_ASD_TD/Idiopathic_ASD_subtype3_vs_Control_subtype1/V2_cutoff_003_maxPerms_25/DMRs/DMRs_annotated.xlsx",
  Idiopathic_ASD_subtype4_vs_Control_subtype1_82 = "/quobyte/lasallegrp/George/ECHO_NBDS/DM_cyto/cytosine_reports_combined/DMR_pairwise_batches_ASD_TD/Idiopathic_ASD_subtype4_vs_Control_subtype1/V4_cutoff_003_maxPerms_82_perGroup_1/DMRs/DMRs_annotated.xlsx",

  Idiopathic_ASD_subtype1_vs_Control_subtype2_25   = "/quobyte/lasallegrp/George/ECHO_NBDS/DM_cyto/cytosine_reports_combined/Males/DMR_pairwise_batches_ASD_TD2/Idiopathic_ASD_subtype1_vs_Control_subtype2/V2_cutoff_003_maxPerms_25/DMRs/DMRs_annotated.xlsx",
  Idiopathic_ASD_subtype2_vs_Control_subtype2_25   = "/quobyte/lasallegrp/George/ECHO_NBDS/DM_cyto/cytosine_reports_combined/Males/DMR_pairwise_batches_ASD_TD2/Idiopathic_ASD_subtype2_vs_Control_subtype2/V2_cutoff_003_maxPerms_25/DMRs/DMRs_annotated.xlsx",
  Idiopathic_ASD_subtype3_vs_Control_subtype2_25  = "/quobyte/lasallegrp/George/ECHO_NBDS/DM_cyto/cytosine_reports_combined/Males/DMR_pairwise_batches_ASD_TD2/Idiopathic_ASD_subtype3_vs_Control_subtype2/V2_cutoff_003_maxPerms_25/DMRs/DMRs_annotated.xlsx",
  Idiopathic_ASD_subtype4_vs_Control_subtype2_25  = "/quobyte/lasallegrp/George/ECHO_NBDS/DM_cyto/cytosine_reports_combined/Males/DMR_pairwise_batches_ASD_TD2/Idiopathic_ASD_subtype4_vs_Control_subtype2/V2_cutoff_003_maxPerms_25/DMRs/DMRs_annotated.xlsx"
  
)

dir.create("upset_outputs", showWarnings = FALSE)

# ---- Helper: read one file and return gene symbols + direction ----
read_dmrs_genes <- function(path) {
  x <- read_excel(path, guess_max = 1e6)

  # extract genes
  genes <- x$geneSymbol
  genes <- genes[!is.na(genes)]
  genes <- unlist(strsplit(genes, "[;,/\\|]"))  # split multi-gene cells
  genes <- str_trim(genes)
  genes <- genes[nzchar(genes)]
  genes <- unique(genes)

  # direction
  dir_vec <- NULL
  if ("direction" %in% names(x)) {
    df <- x %>%
      filter(!is.na(geneSymbol)) %>%
      mutate(gene = strsplit(geneSymbol, "[;,/\\|]")) %>%
      tidyr::unnest_longer(gene) %>%
      mutate(gene = str_trim(gene)) %>%
      filter(gene != "") %>%
      select(gene, direction)

    dir_map <- df %>%
      group_by(gene) %>%
      summarize(direction = {
        d <- unique(na.omit(direction))
        if (!length(d)) NA_character_
        else if (length(d) == 1) d
        else "Mixed"
      }, .groups = "drop")

    dir_vec <- setNames(dir_map$direction, dir_map$gene)
  }

  tibble(gene = genes,
         direction = if (is.null(dir_vec)) NA_character_ else unname(dir_vec[genes]))
}

# ---- Load all versions ----
gene_list <- lapply(files, read_dmrs_genes)

# ---- Build membership ----
all_genes <- sort(unique(unlist(lapply(gene_list, `[[`, "gene"))))
membership <- tibble(gene = all_genes)
for (lab in names(gene_list)) {
  membership[[lab]] <- membership$gene %in% gene_list[[lab]]$gene
}
write.csv(membership, "upset_outputs/GENE_membership_all.csv", row.names = FALSE)

# --- 1) Read & clean SFARI gene symbols (robust to header variants) ---
sfari_path <- "/quobyte/lasallegrp/George/ECHO_NBDS/DM_cyto/cytosine_reports_combined/SFARI-Gene_genes_04-03-2025release_07-02-2025export.csv"
# Re-read the file to reset the mangled column names
SFARI_Genes_raw <- readr::read_csv(sfari_path, show_col_types = FALSE)

# SAFER normalization: use POSIX classes (no perl), or set perl=TRUE if you prefer \s
orig_names <- names(SFARI_Genes_raw)
# POSIX class for whitespace
norm_names <- tolower(gsub("[[:space:].-]+", "_", orig_names))

names(SFARI_Genes_raw) <- norm_names

# Now pick the gene symbol column robustly
gene_col_candidates <- c("gene_symbol", "gene_name", "symbol", "gene", "genes")
gcol <- intersect(gene_col_candidates, names(SFARI_Genes_raw))[1]
if (is.na(gcol)) {
  stop("Couldn't find a gene symbol column. Columns now are: ",
       paste(names(SFARI_Genes_raw), collapse = ", "))
}

# Extract/clean symbols
SFARI_gene_symbols <- SFARI_Genes_raw[[gcol]] |>
  as.character() |>
  (\(x) x[!is.na(x)])() |>
  stringr::str_split("[;,/\\|]") |>
  unlist() |>
  stringr::str_trim() |>
  (\(x) x[nzchar(x)])() |>
  unique()

cat("SFARI genes detected:", length(SFARI_gene_symbols), "\n")

# --- 2) Build FULL membership from raw set genes (DMR sets + SFARI) ---
# You already have gene_list: a list of tibbles each with a column `gene` from your Excel files.

# Collect each set’s full gene vector from the original gene_list (not the trimmed membership)
set_genes <- lapply(gene_list, function(df) unique(as.character(df$gene)))
# add SFARI as its own set
set_genes$SFARI <- unique(as.character(SFARI_gene_symbols))

# Universe = union of ALL genes across ALL sets
all_genes_expanded <- sort(unique(unlist(set_genes)))

# Build expanded membership: one logical column per set
membership_expanded <- tibble::tibble(gene = all_genes_expanded)
for (s in names(set_genes)) {
  membership_expanded[[s]] <- membership_expanded$gene %in% set_genes[[s]]
}

# Doublecheck: each set size (Observed) should match raw vector length (Expected)
expected_sizes <- vapply(set_genes, length, integer(1))
observed_sizes <- vapply(names(set_genes), function(s) sum(membership_expanded[[s]]), integer(1))
sizes_check <- rbind(Expected = expected_sizes, Observed = observed_sizes)
print(sizes_check)

if (!all(expected_sizes == observed_sizes)) {
  warning("One or more set sizes do not match expected counts. Check upstream gene cleaning.")
}

# Save expanded matrix
dir.create("upset_outputs", showWarnings = FALSE)
write.csv(membership_expanded, "upset_outputs/GENE_membership_ALLSETS_expanded.csv", row.names = FALSE)

# --- 3) Plot with expanded membership (all bars now reflect total per source) ---
sets <- setdiff(names(membership_expanded), "gene")

# (Avoid any helper-name collisions by namespacing)
p_all <- ComplexUpset::upset(
  membership_expanded,
  intersect = sets,
  name = "Intersection size (genes)",
  width_ratio = 0.33,        # tweak 0.2–0.5 for readability
  sort_intersections = "descending",
  min_size = 3
)

ggsave("upset_outputs/UpSet_ALLSETS_with_SFARI_expanded.pdf", p_all, width = 30, height = 15)
# ggsave("upset_outputs/UpSet_ALLSETS_with_SFARI_expanded.png", p_all, width = 30, height = 15, dpi = 300)

# # ---- UpSet: ALL ----
# p_all <- upset(
#   membership,
#   intersect = names(files),
#   name = "Intersection size (genes)",
#   width_ratio = 0.33,
#   sort_intersections = "descending",
#   min_size = 1
# )
# ggsave("upset_outputs/UpSet_all_GENES.pdf", p_all, width = 9, height = 6)
# ggsave("upset_outputs/UpSet_all_GENES.png", p_all, width = 9, height = 6, dpi = 300)





# ---- Split by direction ----
long_dir <- bind_rows(lapply(names(gene_list), function(lab) {
  gene_list[[lab]] %>% mutate(set = lab) %>% select(gene, set, direction)
}))

make_dir_membership_gene <- function(target_dir = c("Hyper","Hypo")) {
  target_dir <- match.arg(target_dir)
  by_set <- lapply(names(files), function(lab) {
    ids <- long_dir %>%
      filter(set == lab,
             !is.na(direction),
             direction == target_dir) %>%
      pull(gene)
    tibble(set = lab, gene = unique(ids))
  }) %>% bind_rows()

  ids <- sort(unique(by_set$gene))
  if (!length(ids)) return(NULL)

  out <- tibble(gene = ids)
  for (lab in names(files)) {
    out[[lab]] <- out$gene %in% (by_set %>% filter(set == lab) %>% pull(gene))
  }
  out
}

# Hyper genes
df_hyper_g <- make_dir_membership_gene("Hyper")
if (!is.null(df_hyper_g)) {
  p_hyper <- upset(
    df_hyper_g,
    intersect = names(files),
    name = "Hyper genes",
    width_ratio = 0.33,
    sort_intersections = "descending",
    min_size = 1
  )
  ggsave("upset_outputs/UpSet_hyper_GENES.pdf", p_hyper, width = 9, height = 6)
}

# Hypo genes
df_hypo_g <- make_dir_membership_gene("Hypo")
if (!is.null(df_hypo_g)) {
  p_hypo <- upset(
    df_hypo_g,
    intersect = names(files),
    name = "Hypo genes",
    width_ratio = 0.33,
    sort_intersections = "descending",
    min_size = 1
  )
  ggsave("upset_outputs/UpSet_hypo_GENES.pdf", p_hypo, width = 9, height = 6)
}

# ---- Export intersections ----
extract_intersection <- function(df_bool, sets_keep) {
  others <- setdiff(names(files), sets_keep)
  keep <- df_bool
  for (s in sets_keep) keep <- keep %>% filter(.data[[s]] == TRUE)
  if (length(others)) for (s in others) keep <- keep %>% filter(.data[[s]] == FALSE)
  keep$gene
}

# Example: genes unique to V1 ∩ V2
genes_v1v2 <- extract_intersection(membership, c("V1_cutoff_003","V2_cutoff_003_maxPerms_25"))
 write_csv(genes_v1v2, "upset_outputs/intersection_GENES_V1V2_only.csv")
















#The code below does the upset by the genomic regions. But it failed when i need to plot it. 
# # ---- File paths & labels ----
# files <- c(
#   V1_cutoff_003 = "/quobyte/lasallegrp/George/ECHO_NBDS/DM_cyto/cytosine_reports_combined/DMR_pairwise_batches_ASD_TD/Idiopathic_ASD_subtype1_vs_Control_subtype1/V1_cutoff_003/DMRs/DMRs_annotated.xlsx",
#   V2_cutoff_003_maxPerms_25 = "/quobyte/lasallegrp/George/ECHO_NBDS/DM_cyto/cytosine_reports_combined/DMR_pairwise_batches_ASD_TD/Idiopathic_ASD_subtype1_vs_Control_subtype1/V2_cutoff_003_maxPerms_25/DMRs/DMRs_annotated.xlsx",
#   V3_cutoff_003_maxPerms_71 = "/quobyte/lasallegrp/George/ECHO_NBDS/DM_cyto/cytosine_reports_combined/DMR_pairwise_batches_ASD_TD/Idiopathic_ASD_subtype1_vs_Control_subtype1/V3_cutoff_003_maxPerms_71/DMRs/DMRs_annotated.xlsx",
#   V4_cutoof_003_maxPerms_71_perGroup_1 = "/quobyte/lasallegrp/George/ECHO_NBDS/DM_cyto/cytosine_reports_combined/DMR_pairwise_batches_ASD_TD/Idiopathic_ASD_subtype1_vs_Control_subtype1/V4_cutoof_003_maxPerms_71_perGroup_1/DMRs/DMRs_annotated.xlsx"
# )

# dir.create("upset_outputs", showWarnings = FALSE)

# # ---- Helper: read a DMR file & standardize columns ----
# read_dmrs <- function(path) {
#   x <- read_excel(path, guess_max = 1e6)
#   nm <- tolower(names(x))
#   names(x) <- nm

#   # Build genomic ID
#   if (all(c("chr","start","end") %in% nm)) {
#     id <- paste0(x$chr, ":", x$start, "-", x$end)
#   } else if ("region" %in% nm) {
#     parts <- tidyr::separate_wider_regex(
#       x$region,
#       patterns = c("chr" = "chr[^:]+", ":", "start" = "[0-9]+", "-", "end" = "[0-9]+"),
#       cols_remove = FALSE
#     )
#     id <- paste0(parts$chr, ":", parts$start, "-", parts$end)
#   } else {
#     stop("Could not find chr/start/end or region column in: ", path)
#   }

#   # Direction (prefer explicit; else infer from diff/stat column)
#   dir_col <- NULL
#   if ("direction" %in% nm) dir_col <- x$direction
#   if (is.null(dir_col)) {
#     cand <- c("meandiff","methdiff","avgmethdiff","mean_diff","diff","stat","areastat")
#     cand <- cand[cand %in% nm]
#     if (length(cand)) {
#       v <- suppressWarnings(as.numeric(x[[cand[1]]]))
#       if (any(!is.na(v))) {
#         dir_col <- ifelse(v > 0, "Hyper", ifelse(v < 0, "Hypo", NA))
#       }
#     }
#   }

#   tibble(id = id, direction = if (!is.null(dir_col)) as.character(dir_col) else NA_character_)
# }

# # ---- Load all versions ----
# dmr_list <- lapply(files, read_dmrs)

# # ---- Build membership table (ALL) ----
# all_ids <- sort(unique(unlist(lapply(dmr_list, `[[`, "id"))))
# membership <- tibble(id = all_ids)
# for (lab in names(dmr_list)) {
#   membership[[lab]] <- membership$id %in% dmr_list[[lab]]$id
# }
# write.csv(membership, "upset_outputs/DMR_membership_all.csv", row.names = FALSE)

# # ---- UpSet: ALL DMRs ----
# p_all <- upset(
#   membership,
#   intersect = names(files),
#   name = "Intersection size",
#   width_ratio = 0.33,
#   sort_intersections = "descending",
#   min_size = 1,
#   base_annotations = list(
#     "Intersection size" = intersection_size(text = list(size = 3)),
#     "Set size" = upset_set_size()
#   )
# )
# ggsave("upset_outputs/UpSet_all_DMRs.pdf", p_all, width = 9, height = 6)
# ggsave("upset_outputs/UpSet_all_DMRs.png", p_all, width = 9, height = 6, dpi = 300)

# # ---- Split by direction (Hyper / Hypo) ----
# long_dir <- bind_rows(lapply(names(dmr_list), function(lab) {
#   dmr_list[[lab]] %>% mutate(set = lab) %>% select(id, set, direction)
# }))

# make_dir_membership <- function(target_dir = c("Hyper","Hypo")) {
#   target_dir <- match.arg(target_dir)
#   by_set <- lapply(names(files), function(lab) {
#     ids <- long_dir %>%
#       filter(set == lab, !is.na(direction), grepl(target_dir, direction, ignore.case = TRUE)) %>%
#       pull(id)
#     tibble(set = lab, id = ids)
#   }) %>% bind_rows()

#   ids <- sort(unique(by_set$id))
#   if (!length(ids)) return(NULL)

#   out <- tibble(id = ids)
#   for (lab in names(files)) {
#     out[[lab]] <- out$id %in% (by_set %>% filter(set == lab) %>% pull(id))
#   }
#   out
# }

# # Hyper
# df_hyper <- make_dir_membership("Hyper")
# if (!is.null(df_hyper)) {
#   p_hyper <- upset(
#     df_hyper,
#     intersect = names(files),
#     name = "Hyper DMRs",
#     width_ratio = 0.33,
#     sort_intersections = "descending",
#     min_size = 1,
#     base_annotations = list(
#       "Intersection size" = intersection_size(text = list(size = 3)),
#       "Set size" = upset_set_size()
#     )
#   )
#   ggsave("upset_outputs/UpSet_hyper_DMRs.pdf", p_hyper, width = 9, height = 6)
#   ggsave("upset_outputs/UpSet_hyper_DMRs.png", p_hyper, width = 9, height = 6, dpi = 300)
# }

# # Hypo
# df_hypo <- make_dir_membership("Hypo")
# if (!is.null(df_hypo)) {
#   p_hypo <- upset(
#     df_hypo,
#     intersect = names(files),
#     name = "Hypo DMRs",
#     width_ratio = 0.33,
#     sort_intersections = "descending",
#     min_size = 1,
#     base_annotations = list(
#       "Intersection size" = intersection_size(text = list(size = 3)),
#       "Set size" = upset_set_size()
#     )
#   )
#   ggsave("upset_outputs/UpSet_hypo_DMRs.pdf", p_hypo, width = 9, height = 6)
#   ggsave("upset_outputs/UpSet_hypo_DMRs.png", p_hypo, width = 9, height = 6, dpi = 300)
# }

# # ---- Export intersection lists (for downstream ML/features) ----
# extract_intersection <- function(df_bool, sets_keep) {
#   # df_bool: id + boolean columns
#   # sets_keep: character vector of set names, e.g. c("V1_cutoff_003","V2_cutoff_003_maxPerms_25")
#   others <- setdiff(names(files), sets_keep)
#   keep <- df_bool
#   for (s in sets_keep) keep <- keep %>% filter(.data[[s]] == TRUE)
#   if (length(others)) {
#     for (s in others) keep <- keep %>% filter(.data[[s]] == FALSE)
#   }
#   keep$id
# }

# # Example: IDs unique to V1_cutoff_003 ∩ V2_cutoff_003_maxPerms_25 only
# ids_v1v2 <- extract_intersection(membership, c("V1_cutoff_003","V2_cutoff_003_maxPerms_25"))
# writeLines(ids_v1v2, "upset_outputs/intersection_V1V2_only.txt")

# message("Done. Outputs in upset_outputs/ folder.")
