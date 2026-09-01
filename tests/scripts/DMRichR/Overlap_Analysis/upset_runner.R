# #!/usr/bin/env Rscript


# # source ~/.profile
# # conda activate /quobyte/lasallegrp/programs/.conda/envs/DMRichR4.2_Updated/

# # #you can also just run the R script directly which is called upset_runner.R
# # #Rscript upset_runner.R

# # # activate your conda env
# # source ~/.profile
# # conda activate /quobyte/lasallegrp/programs/.conda/envs/DMRichR4.2_Updated/

# # # change into a comparison folder, then run: #e.g
# # cd /quobyte/lasallegrp/George/ECHO_NBDS/.../Idiopathic_ASD_subtype2_vs_Control_subtype1

# # # Option A: look for all V1/V2/V3... under the parent folder
# # Rscript upset_runner.R             # root=.. (default)

# # # Option B: only look under the current folder
#you are currently running ooption B.
# Rscript /quobyte/lasallegrp/George/ECHO_NBDS/DM_cyto/cytosine_reports_combined/Overlap_Analysis/scripts/upset_runner.R --root . --use_background_universe true --min_size 1 --width_ratio auto --wrap 36 --dpi 300


## this code works, I want to play around with font size. 
## so i added the theme to alter font size
## the theme option actually didnt work
## so what works is changing the height and width of your ggsave output currently using 9 by 6 
suppressPackageStartupMessages({
  library(readxl)
  library(readr)
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(ggplot2)
  library(ComplexUpset)
  library(purrr)
})

# ========================== CLI ARGS ==========================
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  hit <- which(args == flag)
  if (length(hit) == 0 || hit == length(args)) return(default)
  args[hit + 1]
}

root        <- get_arg("--root", "..")          # default: parent of working dir
out_dir     <- get_arg("--out",  "upset_outputs")
min_size    <- as.integer(get_arg("--min_size", "1"))
width_ratio <- get_arg("--width_ratio", "auto") # numeric in (0,1] or "auto"
wrap_chars  <- as.integer(get_arg("--wrap", "36"))
png_dpi     <- as.integer(get_arg("--dpi", "300"))
max_k       <- as.integer(get_arg("--max_k", "2"))  # export intersections up to k sets (2=uniques+pairwise)
use_bg_universe <- tolower(get_arg("--use_background_universe", "false")) %in% c("1","true","t","yes","y")

# --- SFARI gene file (fixed path) ---
sfari_path <- "/quobyte/lasallegrp/George/ECHO_NBDS/DM_cyto/cytosine_reports_combined/SFARI-Gene_genes_04-03-2025release_07-02-2025export.csv"

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(out_dir, "intersections"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(out_dir, "intersections_dmr_rows"), showWarnings = FALSE, recursive = TRUE)

# ========================== LOG HELPERS ==========================
msg  <- function(...) cat("[upset]", ..., "\n")
warn <- function(...) cat("[upset][warn]", ..., "\n")

# ========================== FILE DISCOVERY ==========================
find_dmr_files <- function(root) {
  f <- list.files(path = root, pattern = "^DMRs_annotated\\.xlsx$", recursive = TRUE, full.names = TRUE)
  f <- f[!startsWith(basename(f), "._")]
  info <- file.info(f)
  f[is.finite(info$size) & info$size > 1024]
}

label_from_path <- function(p) {
  dmr_dir  <- dirname(p)                             # .../DMRs
  ver_dir  <- basename(dirname(dmr_dir))             # version folder
  comp_dir <- basename(dirname(dirname(dmr_dir)))    # comparison folder
  paste0(comp_dir, "_", ver_dir)
}

# the version folder containing both DMRs/ and Extra/
version_dir_from_dmr <- function(p) dirname(dirname(p))  # .../<comparison>/<version>

safe_combo_name <- function(sets) {
  nm <- paste(sets, collapse = "__AND__")
  gsub("[^A-Za-z0-9._-]+", "_", nm)
}

# ========================== CLEANING ==========================
clean_genes <- function(v) {
  v <- as.character(v)
  v <- v[!is.na(v)]
  v <- unlist(strsplit(v, "[;,/\\|]"))
  v <- stringr::str_trim(v)
  unique(v[nzchar(v)])
}

# ========================== READERS ==========================
read_dmrs_genes <- function(path) {
  x <- readxl::read_excel(path, guess_max = 1e6)
  if (!"geneSymbol" %in% names(x)) stop("geneSymbol column not found in: ", path)
  clean_genes(x$geneSymbol)
}

# normalize "Hypermethylated"/"Hypomethylated" -> "Hyper"/"Hypo"
read_dmrs_gene_dir <- function(path) {
  x <- readxl::read_excel(path, guess_max = 1e6)
  if (!"geneSymbol" %in% names(x)) stop("geneSymbol column not found in: ", path)
  dir_col_candidates <- c("direction", "Direction")
  dir_col <- dir_col_candidates[dir_col_candidates %in% names(x)][1]
  if (!is.na(dir_col)) {
    dir_vec <- tolower(as.character(x[[dir_col]]))
    dir_vec[grepl("^hyper", dir_vec)] <- "Hyper"
    dir_vec[grepl("^hypo",  dir_vec)] <- "Hypo"
    dir_vec[!dir_vec %in% c("Hyper","Hypo")] <- NA_character_
  } else dir_vec <- NULL

  x %>%
    dplyr::mutate(.dir = dir_vec) %>%
    dplyr::filter(!is.na(geneSymbol)) %>%
    dplyr::mutate(gene = strsplit(geneSymbol, "[;,/\\|]")) %>%
    tidyr::unnest_longer(gene) %>%
    dplyr::mutate(gene = stringr::str_trim(gene)) %>%
    dplyr::filter(gene != "") %>%
    dplyr::group_by(gene) %>%
    dplyr::summarize(
      direction = {
        d <- unique(na.omit(.dir))
        if (!length(d)) NA_character_
        else if (length(d) == 1) d
        else "Mixed"
      },
      .groups = "drop"
    )
}

read_dmrs_rows_exploded <- function(path) {
  x <- readxl::read_excel(path, guess_max = 1e6)
  if (!"geneSymbol" %in% names(x)) stop("geneSymbol column not found in: ", path)
  x %>%
    dplyr::mutate(gene = strsplit(geneSymbol, "[;,/\\|]")) %>%
    tidyr::unnest_longer(gene) %>%
    dplyr::mutate(gene = stringr::str_trim(gene)) %>%
    dplyr::filter(!is.na(gene), gene != "")
}

read_sfari_genes <- function(path) {
  x <- readr::read_csv(path, show_col_types = FALSE)
  orig_names <- names(x)
  names(x) <- tolower(gsub("[[:space:].-]+", "_", orig_names))
  cand <- c("gene_symbol", "gene_name", "symbol", "gene", "genes")
  gcol <- cand[cand %in% names(x)][1]
  if (is.na(gcol)) stop("Could not find a gene symbol column in SFARI CSV. Columns: ",
                        paste(orig_names, collapse = ", "))
  clean_genes(x[[gcol]])
}

# ========================== BACKGROUND UNIVERSE ==========================
annotate_background_for_version <- function(ver_dir) {
  bg_csv  <- file.path(ver_dir, "Extra", "bsseq_background.csv")
  bg_xlsx <- file.path(ver_dir, "Extra", "bsseq_background_annotated.xlsx")

  if (file.exists(bg_xlsx)) {
    msg("Using existing annotated background:", bg_xlsx)
    x <- readxl::read_excel(bg_xlsx, guess_max = 1e6)
    if (!"geneSymbol" %in% names(x)) {
      stop("Background XLSX exists but has no geneSymbol column: ", bg_xlsx)
    }
    return(clean_genes(x$geneSymbol))
  }

  if (!file.exists(bg_csv)) {
    warn("No background files found in ", file.path(ver_dir, "Extra"), " — skipping this version for universe.")
    return(character(0))
  }

  msg("Annotating background CSV → XLSX for:", ver_dir)
  suppressPackageStartupMessages({
    library(GenomicRanges)
    library(ChIPseeker)
    library(TxDb.Hsapiens.UCSC.hg38.knownGene)
    library(org.Hs.eg.db)
    library(openxlsx)
  })

  bg <- readr::read_csv(bg_csv, show_col_types = FALSE)
  nms <- tolower(names(bg))
  chr_ix   <- which(nms %in% c("chr","chrom","chromosome"))[1]
  start_ix <- which(nms %in% c("start","beg","begin","startpos","start_pos"))[1]
  end_ix   <- which(nms %in% c("end","stop","endpos","end_pos"))[1]
  if (any(is.na(c(chr_ix,start_ix,end_ix)))) {
    warn("Cannot find chr/start/end columns in ", bg_csv, " — skipping.")
    return(character(0))
  }

  chr <- bg[[chr_ix]]; start <- bg[[start_ix]]; end <- bg[[end_ix]]
  gr <- GenomicRanges::GRanges(seqnames = chr, ranges = IRanges::IRanges(start, end))
  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
  ann <- ChIPseeker::annotatePeak(gr, TxDb = txdb, annoDb = "org.Hs.eg.db", verbose = FALSE)
  out <- as.data.frame(ann) %>%
    dplyr::transmute(seqnames, start, end, geneSymbol = SYMBOL)

  openxlsx::write.xlsx(out, bg_xlsx, overwrite = TRUE)
  msg("Wrote annotated background:", bg_xlsx, " (", nrow(out), " rows )")
  clean_genes(out$geneSymbol)
}

# ========================== MEMBERSHIP & UTILS ==========================
build_membership <- function(set_genes) {
  all_genes <- sort(unique(unlist(set_genes)))
  mem <- tibble::tibble(gene = all_genes)
  for (s in names(set_genes)) mem[[s]] <- mem$gene %in% set_genes[[s]]
  mem
}

auto_width_ratio <- function(n_sets) {
  if (n_sets >= 16) return(0.18)
  if (n_sets >= 12) return(0.20)
  if (n_sets >= 10) return(0.25)
  if (n_sets >= 8)  return(0.30)
  if (n_sets >= 6)  return(0.33)
  0.40
}

extract_intersection_genes <- function(df_bool, sets_keep) {
  stopifnot("gene" %in% names(df_bool))
  all_sets <- setdiff(names(df_bool), "gene")
  others <- setdiff(all_sets, sets_keep)
  keep <- df_bool
  for (s in sets_keep) keep <- dplyr::filter(keep, .data[[s]] == TRUE)
  if (length(others)) for (s in others) keep <- dplyr::filter(keep, .data[[s]] == FALSE)
  keep$gene
}

make_combos <- function(sets, max_k = 2) {
  combos <- list()
  for (k in 1:min(max_k, length(sets))) {
    combos <- c(combos, combn(sets, k, simplify = FALSE))
  }
  combos
}

# ========================== BIG-TEXT THEME + VERSION-SAFE LABELS ==========================
style_upset_big <- function(p) {
  p + theme(
    text        = element_text(size = 30),
    axis.text.x = element_text(size = 30, angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 30),
    strip.text  = element_text(size = 30, face = "bold"),
    plot.title  = element_text(size = 35, face = "bold", hjust = 0.5),
    legend.title = element_text(size = 30),
    legend.text  = element_text(size = 30)
  )
}

# Try to place numbers on bars when supported; fall back silently if not.
upset_with_labels <- function(data, sets, name, wr, min_size,
                              inter_size = 12, set_size = 12) {

  base_call <- function() {
    ComplexUpset::upset(
      data,
      intersect = sets,
      name = name,
      wrap = as.integer(wrap_chars),
      width_ratio = wr,
      sort_intersections = "descending",
      min_size = min_size
    )
  }

  p <- tryCatch({
    ComplexUpset::upset(
      data,
      intersect = sets,
      name = name,
      wrap = as.integer(wrap_chars),
      width_ratio = wr,
      sort_intersections = "descending",
      min_size = min_size,
      base_annotations = list(
        "Intersection size" = ComplexUpset::intersection_size(
          counts = TRUE,
          text   = list(size = inter_size, vjust = -0.5, hjust = 0.5)
        ),
        "Set size" = ComplexUpset::upset_set_size(
          counts = TRUE,
          text   = list(size = set_size, hjust = -0.1)
        )
      )
    )
  }, error = function(e1) {
    tryCatch({
      ComplexUpset::upset(
        data,
        intersect = sets,
        name = name,
        wrap = as.integer(wrap_chars),
        width_ratio = wr,
        sort_intersections = "descending",
        min_size = min_size,
        base_annotations = list(
          "Intersection size" = ComplexUpset::intersection_size(
            text = list(size = inter_size, vjust = -0.5, hjust = 0.5)
          ),
          "Set size" = ComplexUpset::upset_set_size(
            text = list(size = set_size, hjust = -0.1)
          )
        )
      )
    }, error = function(e2) {
      msg("ComplexUpset is older; drawing without numeric labels on bars. To enable labels, install ComplexUpset >= 1.3.0.")
      base_call()
    })
  })

  style_upset_big(p)
}

# ========================== MAIN ==========================
msg("Searching for DMR files under:", normalizePath(root))
dmr_files <- find_dmr_files(root)
if (length(dmr_files) == 0) stop("No DMRs_annotated.xlsx found under: ", root)

labels <- vapply(dmr_files, label_from_path, character(1))
if (anyDuplicated(labels)) {
  t <- table(labels); dups <- names(t[t > 1])
  for (d in dups) {
    idx <- which(labels == d)
    labels[idx] <- paste0(d, "_", seq_along(idx))
  }
}
names(dmr_files) <- labels

msg("Found", length(dmr_files), "DMR files (sets):")
msg(paste0(" - ", labels, collapse = "\n"))

# Read sets + direction
set_genes <- vector("list", length(labels)); names(set_genes) <- labels
gene_dir_list <- vector("list", length(labels)); names(gene_dir_list) <- labels
for (nm in labels) {
  set_genes[[nm]]     <- read_dmrs_genes(dmr_files[[nm]])
  gene_dir_list[[nm]] <- read_dmrs_gene_dir(dmr_files[[nm]])
}

# SFARI
msg("Reading SFARI genes from:", sfari_path)
set_genes$SFARI <- read_sfari_genes(sfari_path)
msg("SFARI genes:", length(set_genes$SFARI))

# -------------------------- UNIVERSE FROM BACKGROUND --------------------------
U_genes <- NULL
if (use_bg_universe) {
  ver_dirs <- unique(vapply(dmr_files, version_dir_from_dmr, character(1)))
  msg("Preparing background-based universe from ", length(ver_dirs), " version folder(s)…")
  msg(paste0("  • ", ver_dirs, collapse = "\n"))

  invisible(lapply(ver_dirs, function(vd) {
    xlsx <- file.path(vd, "Extra", "bsseq_background_annotated.xlsx")
    csv  <- file.path(vd, "Extra", "bsseq_background.csv")
    if (file.exists(xlsx)) msg("  ✓ Found annotated background: ", xlsx)
    else if (file.exists(csv)) msg("  … Will annotate: ", csv)
    else warn("  ✗ No background in: ", file.path(vd, "Extra"))
  }))

  bg_gene_lists <- lapply(ver_dirs, annotate_background_for_version)
  bg_gene_lists <- bg_gene_lists[lengths(bg_gene_lists) > 0]
  if (length(bg_gene_lists) == 0) {
    warn("Requested --use_background_universe but found no usable backgrounds; using observed-union universe.")
  } else {
    U_genes <- sort(unique(unlist(bg_gene_lists)))
    msg("Universe size (background union): ", length(U_genes))
    for (nm in names(set_genes)) set_genes[[nm]] <- intersect(set_genes[[nm]], U_genes)
  }
} else {
  msg("Using observed-union universe (genes seen in any set, incl. SFARI).")
}

# membership matrix
membership <- build_membership(set_genes)
readr::write_csv(membership, file.path(out_dir, "GENE_membership_ALLSETS_expanded.csv"))

# If background universe not built, use all genes seen
if (is.null(U_genes)) U_genes <- membership$gene

# sanity check
expected_sizes <- vapply(set_genes, length, integer(1))
observed_sizes <- vapply(names(set_genes), function(s) sum(membership[[s]]), integer(1))
sizes_check <- rbind(Expected = expected_sizes, Observed = observed_sizes)
readr::write_csv(as.data.frame(sizes_check), file.path(out_dir, "set_size_check.csv"))
if (!all(expected_sizes == observed_sizes)) {
  warning("One or more set sizes differ from source totals. See set_size_check.csv")
}

# -------------------------- PLOT: MAIN UPSET (BIG TEXT + LABELS IF SUPPORTED) --------------------------
sets <- setdiff(names(membership), "gene")
wr <- if (identical(width_ratio, "auto")) auto_width_ratio(length(sets)) else as.numeric(width_ratio)
msg("Plotting UpSet…  sets:", length(sets), "  min_size:", min_size, "  width_ratio:", wr)

p_main <- upset_with_labels(
  data = membership,
  sets = sets,
  name = "Intersection size (genes)",
  wr = wr,
  min_size = min_size,
  inter_size = 30,  # numbers on intersection bars
  set_size  =  30   # numbers on set-size bars
)

pdf_path <- file.path(out_dir, "UpSet_ALLSETS_expanded.pdf")
png_path <- file.path(out_dir, "UpSet_ALLSETS_expanded.png")
ggsave(pdf_path, p_main, width = 21, height = 4) #16 by 4 for combined plots
ggsave(png_path, p_main, width = 9, height = 6, dpi = png_dpi)

# -------------------------- PLOTS: DIRECTION SPLIT --------------------------
long_dir <- dplyr::bind_rows(lapply(names(gene_dir_list), function(lab) {
  gene_dir_list[[lab]] %>% dplyr::mutate(set = lab) %>% dplyr::select(gene, set, direction)
}))
dir_counts <- long_dir %>%
  dplyr::mutate(direction = dplyr::coalesce(direction, "NA")) %>%
  dplyr::count(set, direction, name = "n") %>%
  dplyr::arrange(set, desc(n))
print(dir_counts)

make_dir_membership_gene <- function(target_dir = c("Hyper","Hypo")) {
  target_dir <- match.arg(target_dir)
  by_set <- lapply(names(gene_dir_list), function(lab) {
    ids <- long_dir %>%
      dplyr::filter(set == lab, !is.na(direction), direction == target_dir) %>%
      dplyr::pull(gene) %>%
      unique()
    tibble::tibble(set = lab, gene = ids)
  }) %>% dplyr::bind_rows()

  ids <- sort(unique(by_set$gene))
  if (!length(ids)) return(NULL)

  out <- tibble::tibble(gene = ids)
  for (lab in names(gene_dir_list)) {
    out[[lab]] <- out$gene %in% (by_set %>% dplyr::filter(set == lab) %>% dplyr::pull(gene))
  }
  out
}

df_hyper_g <- make_dir_membership_gene("Hyper")
if (!is.null(df_hyper_g)) {
  p_hyper <- upset_with_labels(
    data = df_hyper_g,
    sets = setdiff(names(df_hyper_g), "gene"),
    name = "Hyper genes",
    wr   = wr,
    min_size = max(1, min_size),
    inter_size = 14,
    set_size  = 14
  )
  ggsave(file.path(out_dir, "UpSet_HYPER_GENES.pdf"), p_hyper, width = 21, height = 4)
  ggsave(file.path(out_dir, "UpSet_HYPER_GENES.png"), p_hyper, width = 9, height = 6, dpi = png_dpi)
}

df_hypo_g <- make_dir_membership_gene("Hypo")
if (!is.null(df_hypo_g)) {
  p_hypo <- upset_with_labels(
    data = df_hypo_g,
    sets = setdiff(names(df_hypo_g), "gene"),
    name = "Hypo genes",
    wr   = wr,
    min_size = max(1, min_size),
    inter_size = 14,
    set_size  = 14
  )
  ggsave(file.path(out_dir, "UpSet_HYPO_GENES.pdf"), p_hypo, width = 21, height = 4)
  ggsave(file.path(out_dir, "UpSet_HYPO_GENES.png"), p_hypo, width = 9, height = 6, dpi = png_dpi)
}

# -------------------------- EXPORT: PER-INTERSECTION GENES & DMR ROWS --------------------------
msg("Exporting per-intersection gene lists and DMR rows (up to k =", max_k, ")…")

dmr_rows_exploded <- list()
for (lab in names(dmr_files)) {
  dmr_rows_exploded[[lab]] <- read_dmrs_rows_exploded(dmr_files[[lab]]) %>%
    mutate(.set = lab)
}

combos <- make_combos(sets, max_k = max_k)

intersection_index <- tibble::tibble(
  combo_label = character(0),
  k = integer(0),
  n_genes = integer(0),
  n_dmr_rows = integer(0),
  file_genes = character(0),
  file_dmr_rows = character(0)
)

for (cmb in combos) {
  genes_cmb <- extract_intersection_genes(membership, cmb)
  if (length(genes_cmb) < min_size) next

  combo_nm   <- safe_combo_name(cmb)
  genes_path <- file.path(out_dir, "intersections", paste0("genes__", combo_nm, "__only.txt"))
  writeLines(sort(unique(genes_cmb)), genes_path)

  rows_cmb <- bind_rows(lapply(cmb, function(lab) {
    if (lab %in% names(dmr_rows_exploded)) {
      dmr_rows_exploded[[lab]] %>% filter(gene %in% genes_cmb) %>% mutate(.source_set = lab)
    } else NULL
  }))

  dmr_path <- file.path(out_dir, "intersections_dmr_rows", paste0("dmr_rows__", combo_nm, "__only.csv"))
  if (nrow(rows_cmb) > 0) {
    rows_cmb <- distinct(rows_cmb)
    write_csv(rows_cmb, dmr_path)
  } else {
    write_csv(tibble(), dmr_path)
  }

  intersection_index <- add_row(
    intersection_index,
    combo_label = paste(cmb, collapse = " ∩ "),
    k = length(cmb),
    n_genes = length(genes_cmb),
    n_dmr_rows = nrow(rows_cmb),
    file_genes = genes_path,
    file_dmr_rows = dmr_path
  )
}

write_csv(intersection_index, file.path(out_dir, "intersections_index.csv"))
msg("Per-intersection exports complete. Wrote index with", nrow(intersection_index), "entries.")

# -------------------------- STATS: PAIRWISE FISHER (ENRICHMENT-ONLY) --------------------------
sets_pw <- setdiff(names(membership), "gene")
U_genes <- unique(U_genes)
N <- length(U_genes)

gene.list <- lapply(sets_pw, function(s) intersect(membership$gene[membership[[s]]], U_genes))
names(gene.list) <- sets_pw

fisher_enrich_greater <- function(A, B, U) {
  A <- intersect(A, U); B <- intersect(B, U)
  n1 <- length(A); n2 <- length(B); k <- length(intersect(A,B)); N <- length(U)
  a <- k; b <- n1 - k; c <- n2 - k; d <- N - n1 - n2 + k
  ft <- fisher.test(matrix(c(a,b,c,d), 2, byrow = TRUE), alternative = "greater")
  expected <- (n1 * n2) / N
  tibble::tibble(
    k = k, n1 = n1, n2 = n2, N = N,
    expected   = expected,
    enrichment = ifelse(expected > 0, k/expected, NA_real_),
    jaccard    = ifelse(n1 + n2 - k > 0, k/(n1 + n2 - k), NA_real_),
    odds_ratio = unname(ifelse(is.finite(ft$estimate), ft$estimate, NA_real_)),
    p_value    = ft$p.value
  )
}

if (length(sets_pw) >= 2) {
  pairs <- combn(sets_pw, 2, simplify = FALSE)

  pair_stats <- purrr::map_df(pairs, function(pr) {
    s1 <- pr[[1]]; s2 <- pr[[2]]
    fisher_enrich_greater(gene.list[[s1]], gene.list[[s2]], U_genes) %>%
      dplyr::mutate(set1 = s1, set2 = s2, .before = 1)
  }) %>%
    dplyr::mutate(p_adj_BH = p.adjust(p_value, method = "BH")) %>%
    dplyr::arrange(p_value)

  readr::write_csv(pair_stats, file.path(out_dir, "pairwise_overlap_fisher_enrichment_only.csv"))

  if ("SFARI" %in% sets_pw) {
    sfari_pairs <- setdiff(sets_pw, "SFARI") %>%
      purrr::map_df(function(s) {
        fisher_enrich_greater(gene.list[[s]], gene.list[["SFARI"]], U_genes) %>%
          dplyr::mutate(set1 = s, set2 = "SFARI", .before = 1)
      }) %>%
      dplyr::mutate(p_adj_BH = p.adjust(p_value, method = "BH")) %>%
      dplyr::arrange(p_value)
    readr::write_csv(sfari_pairs, file.path(out_dir, "pairwise_overlap_vs_SFARI_enrichment_only.csv"))
  }
} else {
  msg("Skipping pairwise stats: need at least 2 sets (found ", length(sets_pw), ").")
}

# -------------------------- FINAL MESSAGES --------------------------
msg("Done.")
msg("Wrote:")
msg(" - ", normalizePath(file.path(out_dir, "GENE_membership_ALLSETS_expanded.csv")))
msg(" - ", normalizePath(file.path(out_dir, "set_size_check.csv")))
msg(" - ", normalizePath(file.path(out_dir, "UpSet_ALLSETS_expanded.pdf")))
msg(" - ", normalizePath(file.path(out_dir, "UpSet_ALLSETS_expanded.png")))
msg(" - ", normalizePath(file.path(out_dir, "intersections_index.csv")))
msg("   + per-intersection genes in ./intersections/")
msg("   + per-intersection DMR rows in ./intersections_dmr_rows/")
if (file.exists(file.path(out_dir, "UpSet_HYPER_GENES.pdf"))) msg(" - ", normalizePath(file.path(out_dir, "UpSet_HYPER_GENES.pdf")))
if (file.exists(file.path(out_dir, "UpSet_HYPER_GENES.png"))) msg(" - ", normalizePath(file.path(out_dir, "UpSet_HYPER_GENES.png")))
if (file.exists(file.path(out_dir, "UpSet_HYPO_GENES.pdf")))  msg(" - ", normalizePath(file.path(out_dir, "UpSet_HYPO_GENES.pdf")))
if (file.exists(file.path(out_dir, "UpSet_HYPO_GENES.png")))  msg(" - ", normalizePath(file.path(out_dir, "UpSet_HYPO_GENES.png")))
if (file.exists(file.path(out_dir, "pairwise_overlap_fisher_enrichment_only.csv"))) msg(" - ", normalizePath(file.path(out_dir, "pairwise_overlap_fisher_enrichment_only.csv")))
if (file.exists(file.path(out_dir, "pairwise_overlap_vs_SFARI_enrichment_only.csv"))) msg(" - ", normalizePath(file.path(out_dir, "pairwise_overlap_vs_SFARI_enrichment_only.csv")))

## this code works, I just commented it out because I want to play around with font size. 
## so i added the theme to alter font size
# suppressPackageStartupMessages({
#   library(readxl)
#   library(readr)
#   library(dplyr)
#   library(stringr)
#   library(tidyr)
#   library(ggplot2)
#   library(ComplexUpset)
#   library(purrr)
# })

# # ========================== CLI ARGS ==========================
# args <- commandArgs(trailingOnly = TRUE)
# get_arg <- function(flag, default = NULL) {
#   hit <- which(args == flag)
#   if (length(hit) == 0 || hit == length(args)) return(default)
#   args[hit + 1]
# }

# root        <- get_arg("--root", "..")          # default: parent of working dir
# out_dir     <- get_arg("--out",  "upset_outputs")
# min_size    <- as.integer(get_arg("--min_size", "1"))
# width_ratio <- get_arg("--width_ratio", "auto") # numeric in (0,1] or "auto"
# wrap_chars  <- as.integer(get_arg("--wrap", "36"))
# png_dpi     <- as.integer(get_arg("--dpi", "300"))
# max_k       <- as.integer(get_arg("--max_k", "2"))  # export intersections up to k sets (2=uniques+pairwise)
# use_bg_universe <- tolower(get_arg("--use_background_universe", "false")) %in% c("1","true","t","yes","y")

# # --- SFARI gene file (fixed path) ---
# sfari_path <- "/quobyte/lasallegrp/George/ECHO_NBDS/DM_cyto/cytosine_reports_combined/SFARI-Gene_genes_04-03-2025release_07-02-2025export.csv"

# dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
# dir.create(file.path(out_dir, "intersections"), showWarnings = FALSE, recursive = TRUE)
# dir.create(file.path(out_dir, "intersections_dmr_rows"), showWarnings = FALSE, recursive = TRUE)

# # ========================== HELPERS ==========================
# msg  <- function(...) cat("[upset]", ..., "\n")
# warn <- function(...) cat("[upset][warn]", ..., "\n")

# # Find DMR excel tables
# find_dmr_files <- function(root) {
#   f <- list.files(path = root, pattern = "^DMRs_annotated\\.xlsx$", recursive = TRUE, full.names = TRUE)
#   f <- f[!startsWith(basename(f), "._")]
#   info <- file.info(f)
#   f[is.finite(info$size) & info$size > 1024]
# }

# # Version label & folder
# label_from_path <- function(p) {
#   dmr_dir  <- dirname(p)                             # .../DMRs
#   ver_dir  <- basename(dirname(dmr_dir))             # version folder
#   comp_dir <- basename(dirname(dirname(dmr_dir)))    # comparison folder
#   paste0(comp_dir, "_", ver_dir)
# }

# # *** FIXED: the version folder containing both DMRs/ and Extra/ ***
# version_dir_from_dmr <- function(p) dirname(dirname(p))  # .../<comparison>/<version>

# safe_combo_name <- function(sets) {
#   nm <- paste(sets, collapse = "__AND__")
#   gsub("[^A-Za-z0-9._-]+", "_", nm)
# }

# clean_genes <- function(v) {
#   v <- as.character(v)
#   v <- v[!is.na(v)]
#   v <- unlist(strsplit(v, "[;,/\\|]"))
#   v <- stringr::str_trim(v)
#   unique(v[nzchar(v)])
# }

# # ----- Read DMRs -----
# read_dmrs_genes <- function(path) {
#   x <- readxl::read_excel(path, guess_max = 1e6)
#   if (!"geneSymbol" %in% names(x)) stop("geneSymbol column not found in: ", path)
#   clean_genes(x$geneSymbol)
# }

# read_dmrs_gene_dir <- function(path) {
#   x <- readxl::read_excel(path, guess_max = 1e6)
#   if (!"geneSymbol" %in% names(x)) stop("geneSymbol column not found in: ", path)
#   dir_col_candidates <- c("direction", "Direction")
#   dir_col <- dir_col_candidates[dir_col_candidates %in% names(x)][1]
#   if (!is.na(dir_col)) {
#     dir_vec <- tolower(as.character(x[[dir_col]]))
#     dir_vec[grepl("^hyper", dir_vec)] <- "Hyper"  # "Hypermethylated" etc.
#     dir_vec[grepl("^hypo",  dir_vec)] <- "Hypo"   # "Hypomethylated" etc.
#     dir_vec[!dir_vec %in% c("Hyper","Hypo")] <- NA_character_
#   } else dir_vec <- NULL

#   x %>%
#     dplyr::mutate(.dir = dir_vec) %>%
#     dplyr::filter(!is.na(geneSymbol)) %>%
#     dplyr::mutate(gene = strsplit(geneSymbol, "[;,/\\|]")) %>%
#     tidyr::unnest_longer(gene) %>%
#     dplyr::mutate(gene = stringr::str_trim(gene)) %>%
#     dplyr::filter(gene != "") %>%
#     dplyr::group_by(gene) %>%
#     dplyr::summarize(
#       direction = {
#         d <- unique(na.omit(.dir))
#         if (!length(d)) NA_character_
#         else if (length(d) == 1) d
#         else "Mixed"
#       },
#       .groups = "drop"
#     )
# }

# read_dmrs_rows_exploded <- function(path) {
#   x <- readxl::read_excel(path, guess_max = 1e6)
#   if (!"geneSymbol" %in% names(x)) stop("geneSymbol column not found in: ", path)
#   x %>%
#     dplyr::mutate(gene = strsplit(geneSymbol, "[;,/\\|]")) %>%
#     tidyr::unnest_longer(gene) %>%
#     dplyr::mutate(gene = stringr::str_trim(gene)) %>%
#     dplyr::filter(!is.na(gene), gene != "")
# }

# # ----- Read SFARI -----
# read_sfari_genes <- function(path) {
#   x <- readr::read_csv(path, show_col_types = FALSE)
#   orig_names <- names(x)
#   names(x) <- tolower(gsub("[[:space:].-]+", "_", orig_names))
#   cand <- c("gene_symbol", "gene_name", "symbol", "gene", "genes")
#   gcol <- cand[cand %in% names(x)][1]
#   if (is.na(gcol)) stop("Could not find a gene symbol column in SFARI CSV. Columns: ",
#                         paste(orig_names, collapse = ", "))
#   clean_genes(x[[gcol]])
# }

# # ========================== BACKGROUND ANNOTATION ==========================
# annotate_background_for_version <- function(ver_dir) {
#   bg_csv  <- file.path(ver_dir, "Extra", "bsseq_background.csv")
#   bg_xlsx <- file.path(ver_dir, "Extra", "bsseq_background_annotated.xlsx")

#   if (file.exists(bg_xlsx)) {
#     msg("Using existing annotated background:", bg_xlsx)
#     x <- readxl::read_excel(bg_xlsx, guess_max = 1e6)
#     if (!"geneSymbol" %in% names(x)) {
#       stop("Background XLSX exists but has no geneSymbol column: ", bg_xlsx)
#     }
#     return(clean_genes(x$geneSymbol))
#   }

#   if (!file.exists(bg_csv)) {
#     warn("No background files found in ", file.path(ver_dir, "Extra"), " — skipping this version for universe.")
#     return(character(0))
#   }

#   msg("Annotating background CSV → XLSX for:", ver_dir)
#   suppressPackageStartupMessages({
#     library(GenomicRanges)
#     library(ChIPseeker)
#     library(TxDb.Hsapiens.UCSC.hg38.knownGene)
#     library(org.Hs.eg.db)
#     library(openxlsx)
#   })

#   bg <- readr::read_csv(bg_csv, show_col_types = FALSE)
#   nms <- tolower(names(bg))
#   chr_ix   <- which(nms %in% c("chr","chrom","chromosome"))[1]
#   start_ix <- which(nms %in% c("start","beg","begin","startpos","start_pos"))[1]
#   end_ix   <- which(nms %in% c("end","stop","endpos","end_pos"))[1]
#   if (any(is.na(c(chr_ix,start_ix,end_ix)))) {
#     warn("Cannot find chr/start/end columns in ", bg_csv, " — skipping.")
#     return(character(0))
#   }

#   chr <- bg[[chr_ix]]; start <- bg[[start_ix]]; end <- bg[[end_ix]]
#   gr <- GenomicRanges::GRanges(seqnames = chr, ranges = IRanges::IRanges(start, end))
#   txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
#   ann <- ChIPseeker::annotatePeak(gr, TxDb = txdb, annoDb = "org.Hs.eg.db", verbose = FALSE)
#   out <- as.data.frame(ann) %>%
#     dplyr::transmute(seqnames, start, end, geneSymbol = SYMBOL)

#   openxlsx::write.xlsx(out, bg_xlsx, overwrite = TRUE)
#   msg("Wrote annotated background:", bg_xlsx, " (", nrow(out), " rows )")
#   clean_genes(out$geneSymbol)
# }

# # ========================== OTHER HELPERS ==========================
# build_membership <- function(set_genes) {
#   all_genes <- sort(unique(unlist(set_genes)))
#   mem <- tibble::tibble(gene = all_genes)
#   for (s in names(set_genes)) mem[[s]] <- mem$gene %in% set_genes[[s]]
#   mem
# }

# auto_width_ratio <- function(n_sets) {
#   if (n_sets >= 16) return(0.18)
#   if (n_sets >= 12) return(0.20)
#   if (n_sets >= 10) return(0.25)
#   if (n_sets >= 8)  return(0.30)
#   if (n_sets >= 6)  return(0.33)
#   0.40
# }

# extract_intersection_genes <- function(df_bool, sets_keep) {
#   stopifnot("gene" %in% names(df_bool))
#   all_sets <- setdiff(names(df_bool), "gene")
#   others <- setdiff(all_sets, sets_keep)
#   keep <- df_bool
#   for (s in sets_keep) keep <- dplyr::filter(keep, .data[[s]] == TRUE)
#   if (length(others)) for (s in others) keep <- dplyr::filter(keep, .data[[s]] == FALSE)
#   keep$gene
# }

# make_combos <- function(sets, max_k = 2) {
#   combos <- list()
#   for (k in 1:min(max_k, length(sets))) {
#     combos <- c(combos, combn(sets, k, simplify = FALSE))
#   }
#   combos
# }

# # ========================== MAIN ==========================
# msg("Searching for DMR files under:", normalizePath(root))
# dmr_files <- find_dmr_files(root)
# if (length(dmr_files) == 0) stop("No DMRs_annotated.xlsx found under: ", root)

# labels <- vapply(dmr_files, label_from_path, character(1))
# if (anyDuplicated(labels)) {
#   t <- table(labels); dups <- names(t[t > 1])
#   for (d in dups) {
#     idx <- which(labels == d)
#     labels[idx] <- paste0(d, "_", seq_along(idx))
#   }
# }
# names(dmr_files) <- labels

# msg("Found", length(dmr_files), "DMR files (sets):")
# msg(paste0(" - ", labels, collapse = "\n"))

# # read sets + direction
# set_genes <- vector("list", length(labels)); names(set_genes) <- labels
# gene_dir_list <- vector("list", length(labels)); names(gene_dir_list) <- labels
# for (nm in labels) {
#   set_genes[[nm]]     <- read_dmrs_genes(dmr_files[[nm]])
#   gene_dir_list[[nm]] <- read_dmrs_gene_dir(dmr_files[[nm]])
# }

# # SFARI
# msg("Reading SFARI genes from:", sfari_path)
# set_genes$SFARI <- read_sfari_genes(sfari_path)
# msg("SFARI genes:", length(set_genes$SFARI))

# # -------------------------- UNIVERSE FROM BACKGROUND --------------------------
# U_genes <- NULL
# if (use_bg_universe) {
#   ver_dirs <- unique(vapply(dmr_files, version_dir_from_dmr, character(1)))
#   msg("Preparing background-based universe from ", length(ver_dirs), " version folder(s)…")
#   msg(paste0("  • ", ver_dirs, collapse = "\n"))

#   invisible(lapply(ver_dirs, function(vd) {
#     xlsx <- file.path(vd, "Extra", "bsseq_background_annotated.xlsx")
#     csv  <- file.path(vd, "Extra", "bsseq_background.csv")
#     if (file.exists(xlsx)) msg("  ✓ Found annotated background: ", xlsx)
#     else if (file.exists(csv)) msg("  … Will annotate: ", csv)
#     else warn("  ✗ No background in: ", file.path(vd, "Extra"))
#   }))

#   bg_gene_lists <- lapply(ver_dirs, annotate_background_for_version)
#   bg_gene_lists <- bg_gene_lists[lengths(bg_gene_lists) > 0]
#   if (length(bg_gene_lists) == 0) {
#     warn("Requested --use_background_universe but found no usable backgrounds; using observed-union universe.")
#   } else {
#     U_genes <- sort(unique(unlist(bg_gene_lists)))
#     msg("Universe size (background union): ", length(U_genes))
#     # Restrict each set to universe genes so set sizes reflect the universe
#     for (nm in names(set_genes)) set_genes[[nm]] <- intersect(set_genes[[nm]], U_genes)
#   }
# } else {
#   msg("Using observed-union universe (genes seen in any set, incl. SFARI).")
# }

# # membership matrix
# membership <- build_membership(set_genes)
# readr::write_csv(membership, file.path(out_dir, "GENE_membership_ALLSETS_expanded.csv"))

# # if background universe not built, use all genes seen
# if (is.null(U_genes)) U_genes <- membership$gene

# # sanity check
# expected_sizes <- vapply(set_genes, length, integer(1))
# observed_sizes <- vapply(names(set_genes), function(s) sum(membership[[s]]), integer(1))
# sizes_check <- rbind(Expected = expected_sizes, Observed = observed_sizes)
# readr::write_csv(as.data.frame(sizes_check), file.path(out_dir, "set_size_check.csv"))
# if (!all(expected_sizes == observed_sizes)) {
#   warning("One or more set sizes differ from source totals. See set_size_check.csv")
# }

# # -------------------------- PLOT: MAIN UPSET --------------------------
# sets <- setdiff(names(membership), "gene")
# wr <- if (identical(width_ratio, "auto")) auto_width_ratio(length(sets)) else as.numeric(width_ratio)

# msg("Plotting UpSet…  sets:", length(sets), "  min_size:", min_size, "  width_ratio:", wr)
# p <- ComplexUpset::upset(
#   membership,
#   intersect = sets,
#   name = "Intersection size (genes)",
#   wrap = as.integer(wrap_chars),
#   width_ratio = wr,
#   sort_intersections = "descending",
#   min_size = min_size
# )
# pdf_path <- file.path(out_dir, "UpSet_ALLSETS_expanded.pdf")
# png_path <- file.path(out_dir, "UpSet_ALLSETS_expanded.png")
# ggsave(pdf_path, p, width = 30, height = 12)
# ggsave(png_path, p, width = 24, height = 12, dpi = png_dpi)

# # -------------------------- PLOTS: DIRECTION SPLIT --------------------------
# long_dir <- dplyr::bind_rows(lapply(names(gene_dir_list), function(lab) {
#   gene_dir_list[[lab]] %>% dplyr::mutate(set = lab) %>% dplyr::select(gene, set, direction)
# }))
# dir_counts <- long_dir %>%
#   dplyr::mutate(direction = dplyr::coalesce(direction, "NA")) %>%
#   dplyr::count(set, direction, name = "n") %>%
#   dplyr::arrange(set, desc(n))
# print(dir_counts)

# make_dir_membership_gene <- function(target_dir = c("Hyper","Hypo")) {
#   target_dir <- match.arg(target_dir)
#   by_set <- lapply(names(gene_dir_list), function(lab) {
#     ids <- long_dir %>%
#       dplyr::filter(set == lab, !is.na(direction), direction == target_dir) %>%
#       dplyr::pull(gene) %>%
#       unique()
#     tibble::tibble(set = lab, gene = ids)
#   }) %>% dplyr::bind_rows()

#   ids <- sort(unique(by_set$gene))
#   if (!length(ids)) return(NULL)

#   out <- tibble::tibble(gene = ids)
#   for (lab in names(gene_dir_list)) {
#     out[[lab]] <- out$gene %in% (by_set %>% dplyr::filter(set == lab) %>% dplyr::pull(gene))
#   }
#   out
# }

# df_hyper_g <- make_dir_membership_gene("Hyper")
# if (!is.null(df_hyper_g)) {
#   p_hyper <- ComplexUpset::upset(
#     df_hyper_g,
#     intersect = setdiff(names(df_hyper_g), "gene"),
#     name = "Hyper genes",
#     width_ratio = wr,
#     sort_intersections = "descending",
#     min_size = max(1, min_size)
#   )
#   ggsave(file.path(out_dir, "UpSet_HYPER_GENES.pdf"), p_hyper, width = 24, height = 10)
#   ggsave(file.path(out_dir, "UpSet_HYPER_GENES.png"), p_hyper, width = 24, height = 10, dpi = png_dpi)
# }

# df_hypo_g <- make_dir_membership_gene("Hypo")
# if (!is.null(df_hypo_g)) {
#   p_hypo <- ComplexUpset::upset(
#     df_hypo_g,
#     intersect = setdiff(names(df_hypo_g), "gene"),
#     name = "Hypo genes",
#     width_ratio = wr,
#     sort_intersections = "descending",
#     min_size = max(1, min_size)
#   )
#   ggsave(file.path(out_dir, "UpSet_HYPO_GENES.pdf"), p_hypo, width = 24, height = 10)
#   ggsave(file.path(out_dir, "UpSet_HYPO_GENES.png"), p_hypo, width = 24, height = 10, dpi = png_dpi)
# }

# # -------------------------- EXPORT: PER-INTERSECTION GENES & DMR ROWS --------------------------
# msg("Exporting per-intersection gene lists and DMR rows (up to k =", max_k, ")…")

# dmr_rows_exploded <- list()
# for (lab in names(dmr_files)) {
#   dmr_rows_exploded[[lab]] <- read_dmrs_rows_exploded(dmr_files[[lab]]) %>%
#     mutate(.set = lab)
# }

# combos <- make_combos(sets, max_k = max_k)

# intersection_index <- tibble::tibble(
#   combo_label = character(0),
#   k = integer(0),
#   n_genes = integer(0),
#   n_dmr_rows = integer(0),
#   file_genes = character(0),
#   file_dmr_rows = character(0)
# )

# for (cmb in combos) {
#   genes_cmb <- extract_intersection_genes(membership, cmb)
#   if (length(genes_cmb) < min_size) next

#   combo_nm   <- safe_combo_name(cmb)
#   genes_path <- file.path(out_dir, "intersections", paste0("genes__", combo_nm, "__only.txt"))
#   writeLines(sort(unique(genes_cmb)), genes_path)

#   rows_cmb <- bind_rows(lapply(cmb, function(lab) {
#     if (lab %in% names(dmr_rows_exploded)) {
#       dmr_rows_exploded[[lab]] %>% filter(gene %in% genes_cmb) %>% mutate(.source_set = lab)
#     } else NULL
#   }))

#   dmr_path <- file.path(out_dir, "intersections_dmr_rows", paste0("dmr_rows__", combo_nm, "__only.csv"))
#   if (nrow(rows_cmb) > 0) {
#     rows_cmb <- distinct(rows_cmb)
#     write_csv(rows_cmb, dmr_path)
#   } else {
#     write_csv(tibble(), dmr_path)
#   }

#   intersection_index <- add_row(
#     intersection_index,
#     combo_label = paste(cmb, collapse = " ∩ "),
#     k = length(cmb),
#     n_genes = length(genes_cmb),
#     n_dmr_rows = nrow(rows_cmb),
#     file_genes = genes_path,
#     file_dmr_rows = dmr_path
#   )
# }

# write_csv(intersection_index, file.path(out_dir, "intersections_index.csv"))
# msg("Per-intersection exports complete. Wrote index with", nrow(intersection_index), "entries.")

# # -------------------------- STATS: PAIRWISE FISHER (ENRICHMENT-ONLY) --------------------------
# sets_pw <- setdiff(names(membership), "gene")
# U_genes <- unique(U_genes)
# N <- length(U_genes)

# gene.list <- lapply(sets_pw, function(s) intersect(membership$gene[membership[[s]]], U_genes))
# names(gene.list) <- sets_pw

# fisher_enrich_greater <- function(A, B, U) {
#   A <- intersect(A, U); B <- intersect(B, U)
#   n1 <- length(A); n2 <- length(B); k <- length(intersect(A,B)); N <- length(U)
#   a <- k; b <- n1 - k; c <- n2 - k; d <- N - n1 - n2 + k
#   ft <- fisher.test(matrix(c(a,b,c,d), 2, byrow = TRUE), alternative = "greater")
#   expected <- (n1 * n2) / N
#   tibble::tibble(
#     k = k, n1 = n1, n2 = n2, N = N,
#     expected   = expected,
#     enrichment = ifelse(expected > 0, k/expected, NA_real_),
#     jaccard    = ifelse(n1 + n2 - k > 0, k/(n1 + n2 - k), NA_real_),
#     odds_ratio = unname(ifelse(is.finite(ft$estimate), ft$estimate, NA_real_)),
#     p_value    = ft$p.value
#   )
# }

# if (length(sets_pw) >= 2) {
#   pairs <- combn(sets_pw, 2, simplify = FALSE)

#   pair_stats <- purrr::map_df(pairs, function(pr) {
#     s1 <- pr[[1]]; s2 <- pr[[2]]
#     fisher_enrich_greater(gene.list[[s1]], gene.list[[s2]], U_genes) %>%
#       dplyr::mutate(set1 = s1, set2 = s2, .before = 1)
#   }) %>%
#     dplyr::mutate(p_adj_BH = p.adjust(p_value, method = "BH")) %>%
#     dplyr::arrange(p_value)

#   readr::write_csv(pair_stats, file.path(out_dir, "pairwise_overlap_fisher_enrichment_only.csv"))

#   if ("SFARI" %in% sets_pw) {
#     sfari_pairs <- setdiff(sets_pw, "SFARI") %>%
#       purrr::map_df(function(s) {
#         fisher_enrich_greater(gene.list[[s]], gene.list[["SFARI"]], U_genes) %>%
#           dplyr::mutate(set1 = s, set2 = "SFARI", .before = 1)
#       }) %>%
#       dplyr::mutate(p_adj_BH = p.adjust(p_value, method = "BH")) %>%
#       dplyr::arrange(p_value)
#     readr::write_csv(sfari_pairs, file.path(out_dir, "pairwise_overlap_vs_SFARI_enrichment_only.csv"))
#   }
# } else {
#   msg("Skipping pairwise stats: need at least 2 sets (found ", length(sets_pw), ").")
# }

# # -------------------------- FINAL MESSAGES --------------------------
# msg("Done.")
# msg("Wrote:")
# msg(" - ", normalizePath(file.path(out_dir, "GENE_membership_ALLSETS_expanded.csv")))
# msg(" - ", normalizePath(file.path(out_dir, "set_size_check.csv")))
# msg(" - ", normalizePath(file.path(out_dir, "UpSet_ALLSETS_expanded.pdf")))
# msg(" - ", normalizePath(file.path(out_dir, "UpSet_ALLSETS_expanded.png")))
# msg(" - ", normalizePath(file.path(out_dir, "intersections_index.csv")))
# msg("   + per-intersection genes in ./intersections/")
# msg("   + per-intersection DMR rows in ./intersections_dmr_rows/")
# if (file.exists(file.path(out_dir, "UpSet_HYPER_GENES.pdf"))) msg(" - ", normalizePath(file.path(out_dir, "UpSet_HYPER_GENES.pdf")))
# if (file.exists(file.path(out_dir, "UpSet_HYPER_GENES.png"))) msg(" - ", normalizePath(file.path(out_dir, "UpSet_HYPER_GENES.png")))
# if (file.exists(file.path(out_dir, "UpSet_HYPO_GENES.pdf")))  msg(" - ", normalizePath(file.path(out_dir, "UpSet_HYPO_GENES.pdf")))
# if (file.exists(file.path(out_dir, "UpSet_HYPO_GENES.png")))  msg(" - ", normalizePath(file.path(out_dir, "UpSet_HYPO_GENES.png")))
# if (file.exists(file.path(out_dir, "pairwise_overlap_fisher_enrichment_only.csv"))) msg(" - ", normalizePath(file.path(out_dir, "pairwise_overlap_fisher_enrichment_only.csv")))
# if (file.exists(file.path(out_dir, "pairwise_overlap_vs_SFARI_enrichment_only.csv"))) msg(" - ", normalizePath(file.path(out_dir, "pairwise_overlap_vs_SFARI_enrichment_only.csv")))

#I dont remember why I commented this one out
# suppressPackageStartupMessages({
#   library(readxl)
#   library(readr)
#   library(dplyr)
#   library(stringr)
#   library(tidyr)
#   library(ggplot2)
#   library(ComplexUpset)
# })

# # ========================== CLI ARGS ==========================
# args <- commandArgs(trailingOnly = TRUE)
# get_arg <- function(flag, default = NULL) {
#   hit <- which(args == flag)
#   if (length(hit) == 0 || hit == length(args)) return(default)
#   args[hit + 1]
# }

# root        <- get_arg("--root", "..")          # search base (default: parent of wd)
# out_dir     <- get_arg("--out",  "upset_outputs")
# min_size    <- as.integer(get_arg("--min_size", "1"))
# width_ratio <- get_arg("--width_ratio", "auto") # numeric in (0,1] or "auto"
# wrap_chars  <- as.integer(get_arg("--wrap", "36"))
# png_dpi     <- as.integer(get_arg("--dpi", "300"))

# # --- SFARI gene file (fixed path as requested) ---
# sfari_path <- "/quobyte/lasallegrp/George/ECHO_NBDS/DM_cyto/cytosine_reports_combined/SFARI-Gene_genes_04-03-2025release_07-02-2025export.csv"

# dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# # ========================== HELPERS ==========================
# msg <- function(...) cat("[upset]", ..., "\n")

# # Find DMR files recursively
# find_dmr_files <- function(root) {
#   list.files(
#     path = root,
#     pattern = "^DMRs_annotated\\.xlsx$",
#     recursive = TRUE,
#     full.names = TRUE
#   )
# }

# # Label from path: <comparison>/<version>/DMRs/DMRs_annotated.xlsx -> "<comparison>_<version>"
# label_from_path <- function(p) {
#   dmr_dir  <- dirname(p)                             # .../DMRs
#   ver_dir  <- basename(dirname(dmr_dir))             # version folder
#   comp_dir <- basename(dirname(dirname(dmr_dir)))    # comparison folder
#   paste0(comp_dir, "_", ver_dir)
# }

# # Read DMR Excel -> vector of unique gene symbols (expects 'geneSymbol')
# read_dmrs_genes <- function(path) {
#   x <- read_excel(path, guess_max = 1e6)
#   if (!"geneSymbol" %in% names(x)) {
#     stop("geneSymbol column not found in: ", path)
#   }
#   genes <- x$geneSymbol
#   genes <- genes[!is.na(genes)]
#   genes <- unlist(strsplit(genes, "[;,/\\|]"))
#   genes <- stringr::str_trim(genes)
#   unique(genes[nzchar(genes)])
# }

# # Read SFARI CSV robustly (header variants)
# read_sfari_genes <- function(path) {
#   x <- read_csv(path, show_col_types = FALSE)
#   orig_names <- names(x)
#   # safer normalization (no perl needed): spaces, dots, dashes -> underscore
#   names(x) <- tolower(gsub("[[:space:].-]+", "_", orig_names))
#   cand <- c("gene_symbol", "gene_name", "symbol", "gene", "genes")
#   gcol <- cand[cand %in% names(x)][1]
#   if (is.na(gcol)) {
#     stop("Could not find a gene symbol column in SFARI CSV. Columns: ",
#          paste(orig_names, collapse = ", "))
#   }
#   genes <- x[[gcol]] |> as.character() |> (\(v) v[!is.na(v)])()
#   genes <- unlist(strsplit(genes, "[;,/\\|]")) |> stringr::str_trim()
#   unique(genes[nzchar(genes)])
# }

# # Build expanded membership from named list: set -> gene vector
# build_membership <- function(set_genes) {
#   all_genes <- sort(unique(unlist(set_genes)))
#   mem <- tibble(gene = all_genes)
#   for (s in names(set_genes)) {
#     mem[[s]] <- mem$gene %in% set_genes[[s]]
#   }
#   mem
# }

# # Auto width ratio based on number of sets
# auto_width_ratio <- function(n_sets) {
#   if (n_sets >= 16) return(0.18)
#   if (n_sets >= 12) return(0.20)
#   if (n_sets >= 10) return(0.25)
#   if (n_sets >= 8)  return(0.30)
#   if (n_sets >= 6)  return(0.33)
#   0.40
# }

# # # Version-safe annotations; add labels on bars if supported
# # make_base_annotations <- function() {
# #   ann <- list(
# #     "Intersection size" = ComplexUpset::intersection_size(),
# #     "Set size"          = ComplexUpset::upset_set_size()
# #   )
# #   # If your ComplexUpset supports 'text=', use it to show numbers on bars
# #   try({
# #     if (packageVersion("ComplexUpset") >= "1.3.0") {
# #       ann[["Intersection size"]] <- ComplexUpset::intersection_size(text = list(size = 3))
# #       ann[["Set size"]]          <- ComplexUpset::upset_set_size(text = list(size = 3, hjust = -0.1))
# #     }
# #   }, silent = TRUE)
# #   ann
# # }

# # ========================== MAIN ==========================
# msg("Searching for DMR files under:", normalizePath(root))
# dmr_files <- find_dmr_files(root)
# if (length(dmr_files) == 0) stop("No DMRs_annotated.xlsx found under: ", root)

# # Build labels; ensure uniqueness if collisions
# labels <- vapply(dmr_files, label_from_path, character(1))
# if (anyDuplicated(labels)) {
#   t <- table(labels)
#   dups <- names(t[t > 1])
#   for (d in dups) {
#     idx <- which(labels == d)
#     labels[idx] <- paste0(d, "_", seq_along(idx))
#   }
# }
# names(dmr_files) <- labels

# msg("Found", length(dmr_files), "DMR files (sets):")
# msg(paste0(" - ", labels, collapse = "\n"))

# # Read all DMR sets
# set_genes <- vector("list", length(labels))
# names(set_genes) <- labels
# for (nm in labels) {
#   set_genes[[nm]] <- read_dmrs_genes(dmr_files[[nm]])
# }

# # Always add SFARI
# msg("Reading SFARI genes from:", sfari_path)
# set_genes$SFARI <- read_sfari_genes(sfari_path)
# msg("SFARI genes:", length(set_genes$SFARI))

# # Expanded membership (each set bar equals its source total)
# membership <- build_membership(set_genes)
# readr::write_csv(membership, file.path(out_dir, "GENE_membership_ALLSETS_expanded.csv"))

# # Sanity: expected vs observed set sizes
# expected_sizes <- vapply(set_genes, length, integer(1))
# observed_sizes <- vapply(names(set_genes), function(s) sum(membership[[s]]), integer(1))
# sizes_check <- rbind(Expected = expected_sizes, Observed = observed_sizes)
# readr::write_csv(as.data.frame(sizes_check), file.path(out_dir, "set_size_check.csv"))
# if (!all(expected_sizes == observed_sizes)) {
#   warning("One or more set sizes differ from source totals. See set_size_check.csv")
# }

# # Plot
# sets <- setdiff(names(membership), "gene")
# wr <- if (identical(width_ratio, "auto")) auto_width_ratio(length(sets)) else as.numeric(width_ratio)
# # ba <- make_base_annotations()

# msg("Plotting UpSet…  sets:", length(sets), "  min_size:", min_size, "  width_ratio:", wr)
# p <- ComplexUpset::upset(
#   membership,
#   intersect = sets,
#   name = "Intersection size (genes)",
#   wrap = as.integer(wrap_chars),
#   width_ratio = wr,
#   sort_intersections = "descending",
#   min_size = min_size
# #   base_annotations = ba
# )

# pdf_path <- file.path(out_dir, "UpSet_ALLSETS_expanded.pdf")
# png_path <- file.path(out_dir, "UpSet_ALLSETS_expanded.png")
# ggsave(pdf_path, p, width = 30, height = 10)
# ggsave(png_path, p, width = 20, height = 10, dpi = png_dpi)

# msg("Done.")
# msg("Wrote:")
# msg(" - ", normalizePath(file.path(out_dir, "GENE_membership_ALLSETS_expanded.csv")))
# msg(" - ", normalizePath(file.path(out_dir, "set_size_check.csv")))
# msg(" - ", normalizePath(pdf_path))
# msg(" - ", normalizePath(png_path))

# # =========================
# # Pairwise overlap statistics
# # =========================
# suppressPackageStartupMessages({
#   library(dplyr)
#   library(readr)
#   library(tidyr)
#   library(purrr)
# })

# # 1) Choose sets (all logical columns except 'gene')
# sets <- setdiff(names(membership_expanded), "gene")
# N <- nrow(membership_expanded)

# # Small helper to run Fisher’s exact test (enrichment)
# fisher_enrich <- function(x_in, y_in, N) {
#   k  <- sum(x_in & y_in)        # overlap
#   n1 <- sum(x_in)               # size set1
#   n2 <- sum(y_in)               # size set2
#   # contingency table
#   mat <- matrix(c(k, n1 - k,
#                   n2 - k, N - n1 - n2 + k),
#                 nrow = 2, byrow = TRUE)
#   ft <- fisher.test(mat, alternative = "greater")
#   tibble(
#     k = k, n1 = n1, n2 = n2, N = N,
#     expected = n1 * n2 / N,
#     enrichment = ifelse(expected > 0, k / expected, NA_real_),
#     jaccard = ifelse(n1 + n2 - k > 0, k / (n1 + n2 - k), NA_real_),
#     odds_ratio = unname(ifelse(is.finite(ft$estimate), ft$estimate, NA_real_)),
#     p_value = ft$p.value
#   )
# }

# # 2) Compute all pairwise overlaps
# pairs <- combn(sets, 2, simplify = FALSE)

# pair_stats <- map_df(pairs, function(pr) {
#   s1 <- pr[[1]]; s2 <- pr[[2]]
#   res <- fisher_enrich(membership_expanded[[s1]], membership_expanded[[s2]], N)
#   res |>
#     mutate(set1 = s1, set2 = s2, .before = 1)
# })

# # 3) Multiple testing correction (BH/FDR)
# pair_stats <- pair_stats |>
#   mutate(p_adj_BH = p.adjust(p_value, method = "BH")) |>
#   arrange(p_adj_BH, p_value)

# # 4) Save results
# dir.create("upset_outputs", showWarnings = FALSE)
# write_csv(pair_stats, "upset_outputs/pairwise_overlap_fisher.csv")

# # 5) (Optional) quick pretty print of the top hits
# print(
#   pair_stats |>
#     transmute(set1, set2, k, n1, n2, N,
#               expected = round(expected, 2),
#               enrichment = round(enrichment, 2),
#               jaccard = round(jaccard, 3),
#               odds_ratio = round(odds_ratio, 2),
#               p_value = signif(p_value, 3),
#               p_adj_BH = signif(p_adj_BH, 3)) |>
#     head(20)
# )






















#same here I dont remember why I commented this one out tto
# #!/usr/bin/env Rscript

# suppressPackageStartupMessages({
#   library(readxl)
#   library(readr)
#   library(dplyr)
#   library(stringr)
#   library(tidyr)
#   library(ggplot2)
#   library(ComplexUpset)
#   library(purrr)
# })

# # ========================== CLI ARGS ==========================
# args <- commandArgs(trailingOnly = TRUE)
# get_arg <- function(flag, default = NULL) {
#   hit <- which(args == flag)
#   if (length(hit) == 0 || hit == length(args)) return(default)
#   args[hit + 1]
# }

# root        <- get_arg("--root", "..")          # default: parent of working dir
# out_dir     <- get_arg("--out",  "upset_outputs")
# min_size    <- as.integer(get_arg("--min_size", "1"))
# width_ratio <- get_arg("--width_ratio", "auto") # numeric in (0,1] or "auto"
# wrap_chars  <- as.integer(get_arg("--wrap", "36"))
# png_dpi     <- as.integer(get_arg("--dpi", "300"))
# max_k       <- as.integer(get_arg("--max_k", "2"))  # export intersections up to k sets (2=uniques+pairwise)

# # --- SFARI gene file (fixed path) ---
# sfari_path <- "/quobyte/lasallegrp/George/ECHO_NBDS/DM_cyto/cytosine_reports_combined/SFARI-Gene_genes_04-03-2025release_07-02-2025export.csv"

# dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
# dir.create(file.path(out_dir, "intersections"), showWarnings = FALSE, recursive = TRUE)
# dir.create(file.path(out_dir, "intersections_dmr_rows"), showWarnings = FALSE, recursive = TRUE)

# # ========================== HELPERS ==========================
# msg <- function(...) cat("[upset]", ..., "\n")

# # Recursively find real DMR files and drop macOS resource forks / empty files
# find_dmr_files <- function(root) {
#   f <- list.files(path = root, pattern = "^DMRs_annotated\\.xlsx$", recursive = TRUE, full.names = TRUE)
#   f <- f[!startsWith(basename(f), "._")]
#   info <- file.info(f)
#   f <- f[is.finite(info$size) & info$size > 1024]  # > 1KB
#   unique(f)
# }

# # Path -> label: <comparison>/<version>/DMRs/DMRs_annotated.xlsx -> "<comparison>_<version>"
# label_from_path <- function(p) {
#   dmr_dir  <- dirname(p)                             # .../DMRs
#   ver_dir  <- basename(dirname(dmr_dir))             # version folder
#   comp_dir <- basename(dirname(dirname(dmr_dir)))    # comparison folder
#   paste0(comp_dir, "_", ver_dir)
# }

# # Safe filename from vector of set names (order matters)
# safe_combo_name <- function(sets) {
#   nm <- paste(sets, collapse = "__AND__")
#   gsub("[^A-Za-z0-9._-]+", "_", nm)
# }

# # Read DMR Excel -> vector of unique gene symbols (expects 'geneSymbol')
# read_dmrs_genes <- function(path) {
#   x <- readxl::read_excel(path, guess_max = 1e6)
#   if (!"geneSymbol" %in% names(x)) stop("geneSymbol column not found in: ", path)
#   genes <- x$geneSymbol
#   genes <- genes[!is.na(genes)]
#   genes <- unlist(strsplit(genes, "[;,/\\|]"))
#   genes <- stringr::str_trim(genes)
#   unique(genes[nzchar(genes)])
# }

# # Read DMR Excel -> tibble(gene, direction) with normalization:
# # "Hypermethylated"/"Hypomethylated" -> "Hyper"/"Hypo"
# read_dmrs_gene_dir <- function(path) {
#   x <- readxl::read_excel(path, guess_max = 1e6)
#   if (!"geneSymbol" %in% names(x)) stop("geneSymbol column not found in: ", path)

#   # normalize textual direction if present
#   dir_col_candidates <- c("direction", "Direction")
#   dir_col <- dir_col_candidates[dir_col_candidates %in% names(x)][1]
#   if (!is.na(dir_col)) {
#     dir_vec <- tolower(as.character(x[[dir_col]]))
#     dir_vec[grepl("^hyper", dir_vec)] <- "Hyper"   # catches "Hypermethylated"
#     dir_vec[grepl("^hypo",  dir_vec)] <- "Hypo"    # catches "Hypomethylated"
#     dir_vec[!dir_vec %in% c("Hyper","Hypo")] <- NA_character_
#   } else {
#     dir_vec <- NULL
#   }

#   # explode genes and collapse to per-gene direction (resolve conflicts as "Mixed")
#   df <- x %>%
#     dplyr::mutate(.dir = dir_vec) %>%
#     dplyr::filter(!is.na(geneSymbol)) %>%
#     dplyr::mutate(gene = strsplit(geneSymbol, "[;,/\\|]")) %>%
#     tidyr::unnest_longer(gene) %>%
#     dplyr::mutate(gene = stringr::str_trim(gene)) %>%
#     dplyr::filter(gene != "")

#   df %>%
#     dplyr::group_by(gene) %>%
#     dplyr::summarize(
#       direction = {
#         d <- unique(na.omit(.dir))
#         if (!length(d)) NA_character_
#         else if (length(d) == 1) d
#         else "Mixed"
#       },
#       .groups = "drop"
#     )
# }

# # Read full DMR table and explode to long gene level (returns original cols + 'gene')
# read_dmrs_rows_exploded <- function(path) {
#   x <- readxl::read_excel(path, guess_max = 1e6)
#   if (!"geneSymbol" %in% names(x)) stop("geneSymbol column not found in: ", path)
#   x %>%
#     dplyr::mutate(gene = strsplit(geneSymbol, "[;,/\\|]")) %>%
#     tidyr::unnest_longer(gene) %>%
#     dplyr::mutate(gene = stringr::str_trim(gene)) %>%
#     dplyr::filter(!is.na(gene), gene != "")
# }

# # Read SFARI CSV robustly (header variants)
# read_sfari_genes <- function(path) {
#   x <- readr::read_csv(path, show_col_types = FALSE)
#   orig_names <- names(x)
#   names(x) <- tolower(gsub("[[:space:].-]+", "_", orig_names))  # safe: no perl
#   cand <- c("gene_symbol", "gene_name", "symbol", "gene", "genes")
#   gcol <- cand[cand %in% names(x)][1]
#   if (is.na(gcol)) stop("Could not find a gene symbol column in SFARI CSV. Columns: ",
#                         paste(orig_names, collapse = ", "))
#   genes <- x[[gcol]] |> as.character() |> (\(v) v[!is.na(v)])()
#   genes <- unlist(strsplit(genes, "[;,/\\|]")) |> stringr::str_trim()
#   unique(genes[nzchar(genes)])
# }

# # Build expanded membership from named list: set -> gene vector
# build_membership <- function(set_genes) {
#   all_genes <- sort(unique(unlist(set_genes)))
#   mem <- tibble::tibble(gene = all_genes)
#   for (s in names(set_genes)) mem[[s]] <- mem$gene %in% set_genes[[s]]
#   mem
# }

# # Auto width ratio based on number of sets
# auto_width_ratio <- function(n_sets) {
#   if (n_sets >= 16) return(0.18)
#   if (n_sets >= 12) return(0.20)
#   if (n_sets >= 10) return(0.25)
#   if (n_sets >= 8)  return(0.30)
#   if (n_sets >= 6)  return(0.33)
#   0.40
# }

# # Extract exact intersection (present in sets_keep, absent in others)
# extract_intersection_genes <- function(df_bool, sets_keep) {
#   stopifnot("gene" %in% names(df_bool))
#   all_sets <- setdiff(names(df_bool), "gene")
#   others <- setdiff(all_sets, sets_keep)
#   keep <- df_bool
#   for (s in sets_keep) keep <- dplyr::filter(keep, .data[[s]] == TRUE)
#   if (length(others)) for (s in others) keep <- dplyr::filter(keep, .data[[s]] == FALSE)
#   keep$gene
# }

# # Create combinations up to size K (1=uniques, 2=pairwise, etc.)
# make_combos <- function(sets, max_k = 2) {
#   combos <- list()
#   for (k in 1:min(max_k, length(sets))) {
#     combos <- c(combos, combn(sets, k, simplify = FALSE))
#   }
#   combos
# }

# # ========================== MAIN ==========================
# msg("Searching for DMR files under:", normalizePath(root))
# dmr_files <- find_dmr_files(root)
# if (length(dmr_files) == 0) stop("No DMRs_annotated.xlsx found under: ", root)

# # Build labels; ensure uniqueness if collisions
# labels <- vapply(dmr_files, label_from_path, character(1))
# if (anyDuplicated(labels)) {
#   t <- table(labels); dups <- names(t[t > 1])
#   for (d in dups) {
#     idx <- which(labels == d)
#     labels[idx] <- paste0(d, "_", seq_along(idx))
#   }
# }
# names(dmr_files) <- labels

# msg("Found", length(dmr_files), "DMR files (sets):")
# msg(paste0(" - ", labels, collapse = "\n"))

# # Read all DMR sets (genes only)
# set_genes <- vector("list", length(labels)); names(set_genes) <- labels
# for (nm in labels) set_genes[[nm]] <- read_dmrs_genes(dmr_files[[nm]])

# # Also read gene + direction per set (for direction-split)
# gene_dir_list <- vector("list", length(labels)); names(gene_dir_list) <- labels
# for (nm in labels) gene_dir_list[[nm]] <- read_dmrs_gene_dir(dmr_files[[nm]])

# # Always add SFARI
# msg("Reading SFARI genes from:", sfari_path)
# set_genes$SFARI <- read_sfari_genes(sfari_path)
# msg("SFARI genes:", length(set_genes$SFARI))

# # Expanded membership (each set bar equals its source total)
# membership <- build_membership(set_genes)
# readr::write_csv(membership, file.path(out_dir, "GENE_membership_ALLSETS_expanded.csv"))

# # Sanity: expected vs observed set sizes
# expected_sizes <- vapply(set_genes, length, integer(1))
# observed_sizes <- vapply(names(set_genes), function(s) sum(membership[[s]]), integer(1))
# sizes_check <- rbind(Expected = expected_sizes, Observed = observed_sizes)
# readr::write_csv(as.data.frame(sizes_check), file.path(out_dir, "set_size_check.csv"))
# if (!all(expected_sizes == observed_sizes)) {
#   warning("One or more set sizes differ from source totals. See set_size_check.csv")
# }

# # -------------------------- PLOT: MAIN UPSET --------------------------
# sets <- setdiff(names(membership), "gene")
# wr <- if (identical(width_ratio, "auto")) auto_width_ratio(length(sets)) else as.numeric(width_ratio)

# msg("Plotting UpSet…  sets:", length(sets), "  min_size:", min_size, "  width_ratio:", wr)
# p <- ComplexUpset::upset(
#   membership,
#   intersect = sets,
#   name = "Intersection size (genes)",
#   wrap = as.integer(wrap_chars),
#   width_ratio = wr,
#   sort_intersections = "descending",
#   min_size = min_size
# )
# pdf_path <- file.path(out_dir, "UpSet_ALLSETS_expanded.pdf")
# png_path <- file.path(out_dir, "UpSet_ALLSETS_expanded.png")
# ggsave(pdf_path, p, width = 30, height = 12)
# ggsave(png_path, p, width = 24, height = 12, dpi = png_dpi)

# # -------------------------- PLOTS: DIRECTION SPLIT --------------------------
# # Build long table: gene, set, direction (SFARI excluded—no direction)
# long_dir <- dplyr::bind_rows(lapply(names(gene_dir_list), function(lab) {
#   gene_dir_list[[lab]] %>% dplyr::mutate(set = lab) %>% dplyr::select(gene, set, direction)
# }))

# # Debug: counts per set to verify direction pickup
# dir_counts <- long_dir %>%
#   dplyr::mutate(direction = dplyr::coalesce(direction, "NA")) %>%
#   dplyr::count(set, direction, name = "n") %>%
#   dplyr::arrange(set, desc(n))
# print(dir_counts)

# make_dir_membership_gene <- function(target_dir = c("Hyper","Hypo")) {
#   target_dir <- match.arg(target_dir)
#   by_set <- lapply(names(gene_dir_list), function(lab) {
#     ids <- long_dir %>%
#       dplyr::filter(set == lab, !is.na(direction), direction == target_dir) %>%
#       dplyr::pull(gene) %>%
#       unique()
#     tibble::tibble(set = lab, gene = ids)
#   }) %>% dplyr::bind_rows()

#   ids <- sort(unique(by_set$gene))
#   if (!length(ids)) return(NULL)

#   out <- tibble::tibble(gene = ids)
#   for (lab in names(gene_dir_list)) {
#     out[[lab]] <- out$gene %in% (by_set %>% dplyr::filter(set == lab) %>% dplyr::pull(gene))
#   }
#   out
# }

# df_hyper_g <- make_dir_membership_gene("Hyper")
# if (!is.null(df_hyper_g)) {
#   p_hyper <- ComplexUpset::upset(
#     df_hyper_g,
#     intersect = setdiff(names(df_hyper_g), "gene"),
#     name = "Hyper genes",
#     width_ratio = wr,
#     sort_intersections = "descending",
#     min_size = max(1, min_size)
#   )
#   ggsave(file.path(out_dir, "UpSet_HYPER_GENES.pdf"), p_hyper, width = 24, height = 10)
#   ggsave(file.path(out_dir, "UpSet_HYPER_GENES.png"), p_hyper, width = 24, height = 10, dpi = png_dpi)
# }

# df_hypo_g <- make_dir_membership_gene("Hypo")
# if (!is.null(df_hypo_g)) {
#   p_hypo <- ComplexUpset::upset(
#     df_hypo_g,
#     intersect = setdiff(names(df_hypo_g), "gene"),
#     name = "Hypo genes",
#     width_ratio = wr,
#     sort_intersections = "descending",
#     min_size = max(1, min_size)
#   )
#   ggsave(file.path(out_dir, "UpSet_HYPO_GENES.pdf"), p_hypo, width = 24, height = 10)
#   ggsave(file.path(out_dir, "UpSet_HYPO_GENES.png"), p_hypo, width = 24, height = 10, dpi = png_dpi)
# }

# # -------------------------- EXPORT: PER-INTERSECTION GENES & DMR ROWS --------------------------
# msg("Exporting per-intersection gene lists and DMR rows (up to k =", max_k, ")…")

# combos <- make_combos(sets, max_k = max_k)

# # Preload exploded DMR rows for each non-SFARI set to speed up subsetting
# dmr_rows_exploded <- list()
# for (lab in names(dmr_files)) {
#   dmr_rows_exploded[[lab]] <- read_dmrs_rows_exploded(dmr_files[[lab]]) %>%
#     mutate(.set = lab)
# }

# intersection_index <- tibble::tibble(
#   combo_label = character(0),
#   k = integer(0),
#   n_genes = integer(0),
#   n_dmr_rows = integer(0),
#   file_genes = character(0),
#   file_dmr_rows = character(0)
# )

# for (cmb in combos) {
#   genes_cmb <- extract_intersection_genes(membership, cmb)
#   if (length(genes_cmb) < min_size) next

#   combo_nm  <- safe_combo_name(cmb)
#   genes_path <- file.path(out_dir, "intersections", paste0("genes__", combo_nm, "__only.txt"))
#   writeLines(sort(unique(genes_cmb)), genes_path)

#   rows_cmb <- bind_rows(lapply(cmb, function(lab) {
#     if (lab %in% names(dmr_rows_exploded)) {
#       dmr_rows_exploded[[lab]] %>% filter(gene %in% genes_cmb) %>% mutate(.source_set = lab)
#     } else {
#       NULL
#     }
#   }))

#   dmr_path <- file.path(out_dir, "intersections_dmr_rows", paste0("dmr_rows__", combo_nm, "__only.csv"))
#   if (nrow(rows_cmb) > 0) {
#     rows_cmb <- distinct(rows_cmb)
#     write_csv(rows_cmb, dmr_path)
#   } else {
#     write_csv(tibble(), dmr_path)
#   }

#   intersection_index <- add_row(
#     intersection_index,
#     combo_label = paste(cmb, collapse = " ∩ "),
#     k = length(cmb),
#     n_genes = length(genes_cmb),
#     n_dmr_rows = nrow(rows_cmb),
#     file_genes = genes_path,
#     file_dmr_rows = dmr_path
#   )
# }

# write_csv(intersection_index, file.path(out_dir, "intersections_index.csv"))
# msg("Per-intersection exports complete. Wrote index with", nrow(intersection_index), "entries.")

# # -------------------------- STATS: PAIRWISE FISHER (two-sided) --------------------------
# sets_pw <- setdiff(names(membership), "gene")
# N <- nrow(membership)

# if (length(sets_pw) >= 2) {
#   fisher_enrich <- function(x_in, y_in, N) {
#     k  <- sum(x_in & y_in)
#     n1 <- sum(x_in)
#     n2 <- sum(y_in)
#     mat <- matrix(c(k, n1 - k,
#                     n2 - k, N - n1 - n2 + k),
#                   nrow = 2, byrow = TRUE)
#     ft <- fisher.test(mat, alternative = "two.sided")
#     tibble::tibble(
#       k = k, n1 = n1, n2 = n2, N = N,
#       expected   = n1 * n2 / N,
#       enrichment = ifelse(expected > 0, k / expected, NA_real_), # >1 enrichment, <1 depletion
#       jaccard    = ifelse(n1 + n2 - k > 0, k / (n1 + n2 - k), NA_real_),
#       odds_ratio = unname(ifelse(is.finite(ft$estimate), ft$estimate, NA_real_)),
#       p_value    = ft$p.value
#     )
#   }

#   pairs <- combn(sets_pw, 2, simplify = FALSE)

#   pair_stats <- purrr::map_df(pairs, function(pr) {
#     s1 <- pr[[1]]; s2 <- pr[[2]]
#     res <- fisher_enrich(membership[[s1]], membership[[s2]], N)
#     dplyr::mutate(res, set1 = s1, set2 = s2, .before = 1)
#   }) %>%
#     dplyr::mutate(p_adj_BH = p.adjust(p_value, method = "BH")) %>%
#     dplyr::arrange(p_adj_BH, p_value)

#   readr::write_csv(pair_stats, file.path(out_dir, "pairwise_overlap_fisher.csv"))

#   if ("SFARI" %in% sets_pw) {
#     sfari_pairs <- setdiff(sets_pw, "SFARI") %>%
#       purrr::map_df(function(s) {
#         fisher_enrich(membership[[s]], membership$SFARI, N) %>%
#           dplyr::mutate(set1 = s, set2 = "SFARI", .before = 1)
#       }) %>%
#       dplyr::mutate(p_adj_BH = p.adjust(p_value, method = "BH")) %>%
#       dplyr::arrange(p_adj_BH, p_value)
#     readr::write_csv(sfari_pairs, file.path(out_dir, "pairwise_overlap_vs_SFARI.csv"))
#   }
# } else {
#   msg("Skipping pairwise stats: need at least 2 sets (found ", length(sets_pw), ").")
# }

# # -------------------------- FINAL MESSAGES --------------------------
# msg("Done.")
# msg("Wrote:")
# msg(" - ", normalizePath(file.path(out_dir, "GENE_membership_ALLSETS_expanded.csv")))
# msg(" - ", normalizePath(file.path(out_dir, "set_size_check.csv")))
# msg(" - ", normalizePath(file.path(out_dir, "UpSet_ALLSETS_expanded.pdf")))
# msg(" - ", normalizePath(file.path(out_dir, "UpSet_ALLSETS_expanded.png")))
# msg(" - ", normalizePath(file.path(out_dir, "intersections_index.csv")))
# msg("   + per-intersection genes in ./intersections/")
# msg("   + per-intersection DMR rows in ./intersections_dmr_rows/")
# if (file.exists(file.path(out_dir, "UpSet_HYPER_GENES.pdf"))) msg(" - ", normalizePath(file.path(out_dir, "UpSet_HYPER_GENES.pdf")))
# if (file.exists(file.path(out_dir, "UpSet_HYPER_GENES.png"))) msg(" - ", normalizePath(file.path(out_dir, "UpSet_HYPER_GENES.png")))
# if (file.exists(file.path(out_dir, "UpSet_HYPO_GENES.pdf")))  msg(" - ", normalizePath(file.path(out_dir, "UpSet_HYPO_GENES.pdf")))
# if (file.exists(file.path(out_dir, "UpSet_HYPO_GENES.png")))  msg(" - ", normalizePath(file.path(out_dir, "UpSet_HYPO_GENES.png")))
# if (file.exists(file.path(out_dir, "pairwise_overlap_fisher.csv"))) msg(" - ", normalizePath(file.path(out_dir, "pairwise_overlap_fisher.csv")))
# if (file.exists(file.path(out_dir, "pairwise_overlap_vs_SFARI.csv"))) msg(" - ", normalizePath(file.path(out_dir, "pairwise_overlap_vs_SFARI.csv")))

