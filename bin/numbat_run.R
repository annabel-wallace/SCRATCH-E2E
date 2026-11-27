#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Run Numbat on one GROUP (one individual, all of their libraries) and emit
# per-cell malignancy calls keyed by the pipeline's own cell ids.
#
# Numbat combines expression-derived CNV signal with allele-specific signal from
# phased heterozygous SNPs. The allelic channel is what expression-only callers
# (inferCNV, SCEVAN, CopyKAT) structurally cannot see: a copy-neutral LOH event
# changes allele ratios without changing dosage, so it produces no expression
# signal at all.
#
# WHY PER INDIVIDUAL AND NOT PER LIBRARY. numbat:::mark_tumor_lineage() assigns
# the tumour compartment with
#     tumor_root = mut_nodes[which.max(mut_burdens)]
#     compartment = ifelse(seq > 0, "tumor", "normal")
# There is no significance test and no null: whatever segments survive min_LLR
# go on a tree, the heaviest clade becomes "tumor". Numbat cannot answer "there
# is no tumour in this sample". Run per library, a PBMC library therefore comes
# back with its own noise partitioned into tumour and normal rather than with
# nothing. Pooled with the individual's tumour libraries, the real aneuploid
# clones dominate the mutation burden, tumor_root lands on the real tumour, and
# the normal cells stay normal.
#
# The pileup is also joint: pileup_and_phase.R was given every library for this
# individual at once (--label <group> --samples lib1,lib2,...), so the
# heterozygous sites are genotyped and phased across all of them. It writes ONE
# allele table per library, with that library's raw barcodes; this script
# namespaces each to "<library>_<barcode>" so they can be concatenated and still
# line up with the pooled count matrix.
#
# Output: <group>_numbat_calls.tsv with cell_id, numbat_call, p_cnv, clone.
# STRATIFY reads cell_id and numbat_call; the rest is for the report.
#
# Usage: numbat_run.R <counts.rds> <cellmap.tsv> <refcells.tsv> <group> <outdir> \
#                     <ncores> <min_cells> <t> <max_iter> <min_LLR> <min_ref_cells> \
#                     <max_entropy> <alleles...>
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(numbat)
  library(Matrix)
  library(data.table)
})

a <- commandArgs(trailingOnly = TRUE)
if (length(a) < 13L)
  stop("usage: numbat_run.R <counts.rds> <cellmap.tsv> <refcells.tsv> <group> <outdir> ",
       "<ncores> <min_cells> <t> <max_iter> <min_LLR> <min_ref_cells> <max_entropy> ",
       "<allele_counts.tsv.gz...>")

counts_path <- a[[1]]; map_path <- a[[2]]; refcells_path <- a[[3]]
group_id    <- a[[4]]; outdir   <- a[[5]]
ncores      <- as.integer(a[[6]])
min_cells   <- as.integer(a[[7]])
t_param     <- as.numeric(a[[8]])
max_iter    <- as.integer(a[[9]])
# Numbat's own failure message says "Consider reducing min_LLR", but the value
# was not reachable from the pipeline. Exposed so that advice is actionable.
min_llr       <- as.numeric(a[[10]])
min_ref_cells <- as.integer(a[[11]])
# Single-cell filtering threshold. Numbat's maintainer names this and min_LLR as
# the two levers for "a tumor proportion much higher than orthogonal
# measurements suggest" -- i.e. CNV false positives in normal cells
# (kharchenkolab/numbat#224). Exposed at Numbat's own default of 0.5, so passing
# it changes nothing until it is deliberately lowered.
max_entropy   <- as.numeric(a[[12]])
if (is.na(max_entropy)) max_entropy <- 0.5
allele_paths  <- a[13:length(a)]

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

count_mat <- readRDS(counts_path)
cellmap   <- fread(map_path)

# A group that yields no calls is a RESULT, not a pipeline failure. Emitting an
# empty-but-valid table keeps the output contract, lets STRATIFY see that this
# group contributed nothing, and -- the reason this exists -- stops one group
# from destroying the whole run. A 7h42m run died here with the OTHER sample's
# calls already complete, and threw them away with it.
emit_empty <- function(outdir, group_id, why) {
  message(sprintf("[%s] NO CALLS: %s", group_id, why))
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  empty <- data.table::data.table(
    cell = character(), barcode = character(), cell_id = character(),
    sample = character(), group = character(),
    numbat_call = character(), p_cnv = numeric(),
    clone_opt = integer(), compartment_opt = character())
  data.table::fwrite(empty,
    file.path(outdir, paste0(group_id, "_numbat_calls.tsv")), sep = "\t")
  # Machine-readable reason, so a report can say WHY, not merely "0 rows".
  writeLines(why, file.path(outdir, paste0(group_id, "_numbat_nocalls.txt")))
  quit(save = "no", status = 0)
}

# ---------------------------------------------------------------------------
# Concatenate the per-library allele tables, namespacing each cell.
#
# pileup_and_phase.R writes {outdir}/{sample}_allele_counts.tsv.gz per library,
# and the `cell` column holds that library's RAW cellranger barcodes. Two
# libraries of one patient both contain AACGT...-1, so concatenating without a
# prefix silently merges distinct cells. The prefix is the library id recovered
# from the filename, which is the same string numbat_prep.R used to build the
# count matrix column names.
# ---------------------------------------------------------------------------
libs_seen <- character(0)
allele_list <- lapply(allele_paths, function(p) {
  lb <- sub("_allele_counts\\.tsv\\.gz$", "", basename(p))
  d  <- fread(p)
  if (!nrow(d)) {
    message(sprintf("[%s] %s: allele table is empty", group_id, lb))
    return(NULL)
  }
  if (!"cell" %in% names(d))
    stop(sprintf("[%s] %s: no 'cell' column in %s", group_id, lb, basename(p)))
  d[, cell := paste0(lb, "_", as.character(cell))]
  libs_seen <<- c(libs_seen, lb)
  d
})
allele_list <- Filter(Negate(is.null), allele_list)
if (!length(allele_list))
  emit_empty(outdir, group_id, sprintf(
    "All %d allele table(s) for this group were empty; the pileup produced no phased sites.",
    length(allele_paths)))

df_allele <- rbindlist(allele_list, use.names = TRUE, fill = TRUE)
rm(allele_list); invisible(gc())

message(sprintf("[%s] %d genes x %d cells | %d allele records from %d librar%s (%s)",
                group_id, nrow(count_mat), ncol(count_mat), nrow(df_allele),
                length(libs_seen), if (length(libs_seen) == 1L) "y" else "ies",
                paste(libs_seen, collapse = ", ")))

# The pooled count matrix and the concatenated allele table must address the
# same cells. A total mismatch means the prefixes disagree, which is silent and
# fatal to the result rather than to the run: Numbat would model expression with
# no allelic evidence at all and fall back to behaving like inferCNV.
shared <- length(intersect(colnames(count_mat), unique(df_allele$cell)))
message(sprintf("[%s] %d cell(s) carry both expression and allele data", group_id, shared))
if (shared == 0L)
  emit_empty(outdir, group_id, paste(
    "No cell appears in both the count matrix and the allele tables, so the",
    "allelic channel is empty. The two are namespaced '<library>_<barcode>' from",
    "the library id; check that the BAM directory names match the library column."))
if (shared < 0.5 * ncol(count_mat))
  message(sprintf(paste("[%s] WARNING: only %.1f%% of cells have allele data.",
                        "Numbat will model the rest from expression alone."),
                  group_id, 100 * shared / ncol(count_mat)))

if (ncol(count_mat) < min_cells)
  emit_empty(outdir, group_id,
             sprintf("%d cells is below --numbat_min_cells=%d; group skipped.",
                     ncol(count_mat), min_cells))

# ---------------------------------------------------------------------------
# Expression reference. Prefer the one numbat_prep.R built from THIS
# individual's own immune and stromal cells: numbat::ref_hca is 12 generic
# profiles from a different study, and the gap between it and the real cohort is
# read as expression CNV signal. The internal reference is why running per
# individual helps twice -- it is also what makes a normal library useful rather
# than merely harmless.
# ---------------------------------------------------------------------------
# numbat_prep.R SELECTS the cells (it has the annotation and the shared label
# matcher); the aggregation happens here because aggregate_counts() needs the
# numbat package, and prep runs in image_cnv, which has Seurat and no numbat.
# The two images are mirrors of each other -- verified, not assumed.
#
# The `group` column handed to aggregate_counts is the lambdas_ref COLUMN name.
# prep writes it as <celltype>__c<cluster>, one column per cluster, because
# Numbat picks the best-fitting column per pseudobulk and the author's own
# reference construction uses raw clusters (vignettes/mouse.Rmd section 4).
# `ref_group` is read in preference to `celltype`; falling back to `celltype`
# keeps an older refcells.tsv working.
#
# count_mat here is RAW counts and aggregate_counts does the normalisation
# (normalized = TRUE). That ordering is a requirement, not a preference: the
# maintainer's instruction is "created using numbat::aggregate_counts() from raw
# counts ... Don't use Seurat normalization or other methods"
# (kharchenkolab/numbat#224). numbat_prep.R reads layer = "counts" for the same
# reason.
ref <- NULL
ref_basis <- "ref_hca"
refcells <- tryCatch(fread(refcells_path), error = function(e) NULL)
grp_col <- if (!is.null(refcells) && "ref_group" %in% names(refcells)) "ref_group" else "celltype"
if (!is.null(refcells) && nrow(refcells) >= min_ref_cells &&
    all(c("cell", grp_col) %in% names(refcells))) {
  refcells <- refcells[cell %in% colnames(count_mat)]
  if (nrow(refcells) >= min_ref_cells) {
    annot <- data.frame(cell = refcells$cell,
                        group = as.character(refcells[[grp_col]]),
                        stringsAsFactors = FALSE)
    ref <- tryCatch(
      numbat::aggregate_counts(count_mat[, annot$cell, drop = FALSE], annot,
                               normalized = TRUE, verbose = FALSE),
      error = function(e) {
        message(sprintf("[%s] aggregate_counts failed (%s); falling back to ref_hca",
                        group_id, conditionMessage(e))); NULL })
    if (!is.null(ref)) {
      ref_basis <- "internal"
      message(sprintf(paste("[%s] internal reference from %d of this individual's own cells",
                            "over %d column(s), grouped by '%s'"),
                      group_id, nrow(annot), length(unique(annot$group)), grp_col))
      if (grp_col == "celltype")
        message(sprintf(paste("[%s] [warn] refcells.tsv has no ref_group column, so the",
                              "reference is ONE COLUMN PER CELL TYPE. That averages several",
                              "clusters into each column; set --numbat_cluster_col so prep",
                              "can write cluster-granular columns."), group_id))
    }
  }
}
if (is.null(ref)) {
  ref <- tryCatch(numbat::ref_hca, error = function(e) NULL)
  message(sprintf(paste("[%s] no internal reference (need %d reference-lineage cells);",
                        "using numbat::ref_hca"), group_id, min_ref_cells))
}
if (is.null(ref))
  stop(sprintf("[%s] no expression reference available (internal missing and numbat::ref_hca unreadable)",
               group_id))
message(sprintf("[%s] expression reference: %s (%d genes x %d profile(s): %s)",
                group_id, ref_basis, nrow(ref), ncol(ref),
                paste(colnames(ref), collapse = ", ")))

out_sub <- file.path(outdir, group_id)
dir.create(out_sub, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Scaling to large individuals (set by NUMBAT_RUN through the environment, so
# the positional arguments above are untouched).
#
# Numbat's initial clustering (numbat:::exp_hclust) smooths every cell's
# expression along the genome, takes Euclidean distances between EVERY pair of
# cells with parallelDist (one pair at a time over all genes), then builds a
# Ward tree. On DFCI1 that was 6 h 25 m of an 8 h run at 35k cells, ~30 h at
# 60k. The tree only seeds Numbat (cutree at init_k = 3).
#
# NUMBAT_INIT_DIST:
#   exact  (default) the SAME Euclidean distances, computed as
#          ||x||^2 + ||y||^2 - 2 x.y with blocked matrix multiplication (BLAS,
#          all task CPUs) and written straight into a 'dist' object. Same
#          numbers to floating-point rounding, so the same tree; verified
#          against Numbat's own hc.rds on DFCI1 (see gram_dist).
#   pca    distances on the top NUMBAT_INIT_PCS principal components
#          (randomized SVD). Faster again, but NOT the same tree: on DFCI1 it
#          moved 3-5% of cells between the init_k = 3 starting groups (30-200
#          PCs). For individuals too large for exact (dist memory is cells^2 x
#          4 bytes; Ward holds about two copies).
#   numbat Numbat's own function, unchanged.
#
# NUMBAT_SKIP_NJ (auto | true | false), NUMBAT_NJ_MAX_CELLS: after calling
#   cells Numbat builds a UPGMA tree AND a neighbour-joining tree over every
#   cell and seeds NNI with the better one. ape::nj is O(cells^3); UPGMA is
#   O(cells^2). auto skips NJ above NUMBAT_NJ_MAX_CELLS cells; NNI still refines
#   the UPGMA seed (run_numbat's own skip_nj).
# ---------------------------------------------------------------------------
gram_dist <- function(X, block = 1024L) {
  # Exact Euclidean distances between the ROWS of X as a 'dist' object.
  # G = X %*% t(X[I, ]) per block of rows I (n x |I|, contiguous per column);
  # column i of the dist vector holds rows i+1..n, starting at
  # n*(i-1) - i*(i-1)/2 + 1. Memory: the dist vector plus one n x block block.
  n <- nrow(X)
  sq <- numeric(n)
  for (b in split(seq_len(n), ceiling(seq_len(n) / 5000))) sq[b] <- rowSums(X[b, , drop = FALSE]^2)
  out <- numeric(n * (n - 1) / 2)
  for (i0 in seq(1L, n - 1L, by = block)) {
    I <- i0:min(i0 + block - 1L, n - 1L)
    G <- tcrossprod(X, X[I, , drop = FALSE])
    for (r in seq_along(I)) {
      i <- I[r]; jj <- (i + 1L):n
      st <- n * (i - 1) - i * (i - 1) / 2
      out[(st + 1):(st + n - i)] <- sqrt(pmax(sq[i] + sq[jj] - 2 * G[jj, r], 0))
    }
  }
  structure(out, Size = n, Labels = rownames(X), Diag = FALSE, Upper = FALSE,
            method = "euclidean", class = "dist", call = quote(gram_dist(X)))
}
# ---- end gram_dist

# R's stats::hclust refuses more than 65,536 observations ("size cannot be NA
# nor exceed 65536"). Numbat builds trees over CELLS twice -- the initial Ward
# tree (exp_hclust) and the UPGMA phylogeny seed (numbat:::upgma) -- so an
# individual above 65,536 cells cannot finish. GBM DFCI7 (79,687 cells) failed
# exactly there, after ~40 h of distance computation. fastcluster::hclust is the
# drop-in replacement: identical merges and heights to stats::hclust (ward.D2
# and average checked, heights within 1e-14) and no size limit (70k cells in
# 3.8 min, emulated). It is in the Numbat image from tag procps-fastcluster.
has_fastcluster <- requireNamespace("fastcluster", quietly = TRUE)
hclust_fun <- if (has_fastcluster) fastcluster::hclust else stats::hclust
if (ncol(count_mat) > 65536L && !has_fastcluster)
  stop(sprintf(paste(
    "[%s] %d cells: above stats::hclust's 65,536-observation limit, and the",
    "fastcluster package is not in this image. Numbat would fail at its first tree",
    "after computing every pairwise distance. Use the procps-fastcluster Numbat image."),
    group_id, ncol(count_mat)))
if (has_fastcluster) {
  .upgma <- function(D, method = "average", ...) {
    # numbat:::upgma (1.4.2) verbatim, with the unlimited hclust.
    DD <- as.dist(D)
    hc <- hclust_fun(DD, method = method, ...)
    result <- ape::as.phylo(hc)
    result <- reorder(result, "postorder")
    result
  }
  .uenv <- new.env(parent = asNamespace("numbat")); .uenv$hclust_fun <- hclust_fun
  environment(.upgma) <- .uenv
  utils::assignInNamespace("upgma", .upgma, ns = "numbat")
}

init_dist <- tolower(Sys.getenv("NUMBAT_INIT_DIST", "numbat"))
init_pcs  <- suppressWarnings(as.integer(Sys.getenv("NUMBAT_INIT_PCS", "50")))
# Numbat's own exp_hclust calls stats::hclust, so above the limit it cannot be
# used even when asked for; the exact path computes the same distances.
if (init_dist == "numbat" && ncol(count_mat) > 65536L) {
  message(sprintf("[%s] %d cells > 65,536: Numbat's own initial clustering cannot run; using 'exact' (same distances).",
                  group_id, ncol(count_mat)))
  init_dist <- "exact"
}
if (init_dist %in% c("exact", "pca")) {
  .env <- new.env(parent = asNamespace("numbat"))
  .env$init_dist <- init_dist; .env$init_pcs <- init_pcs; .env$gram_dist <- gram_dist
  .env$hclust_fun <- hclust_fun
  fast_exp_hclust <- function(count_mat, lambdas_ref, gtf, sc_refs = NULL, window = 101,
                              ncores = 1, verbose = TRUE) {
    # Same first three steps as numbat:::exp_hclust (1.4.2), verbatim.
    count_mat <- check_matrix(count_mat)
    if (is.null(sc_refs)) sc_refs <- choose_ref_cor(count_mat, lambdas_ref, gtf)
    lambdas_bar <- get_lambdas_bar(lambdas_ref, sc_refs, verbose = FALSE)
    gexp_roll_wide <- t(smooth_expression(count_mat, lambdas_bar, gtf,
                                          window = window, verbose = verbose))
    X <- gexp_roll_wide
    # parDist propagates an NA into the distance and exp_hclust then sets that
    # distance to 0; both fast paths need finite input, so say so and use 0.
    if (anyNA(X)) { log_warn("NAs in smoothed expression, set to 0 before the distances."); X[is.na(X)] <- 0 }
    if (init_dist == "pca") {
      k <- min(init_pcs, ncol(X), nrow(X) - 1L)
      # Top-k PCs by randomized SVD (Halko, Martinsson & Tropp 2011; 10
      # oversampling columns, 2 power iterations), centring applied implicitly
      # so the dense matrix is never copied. irlba/RSpectra are not in the
      # Numbat image. Seeded, so a rerun gives the same tree.
      mu <- colMeans(X); n <- nrow(X); l <- min(k + 10L, ncol(X))
      xc_mul  <- function(M) X %*% M - matrix(drop(crossprod(mu, M)), n, ncol(M), byrow = TRUE)
      xct_mul <- function(Q) crossprod(X, Q) - tcrossprod(mu, colSums(Q))
      set.seed(20261117)
      Q <- qr.Q(qr(xc_mul(matrix(stats::rnorm(ncol(X) * l), ncol(X), l))))
      for (it in 1:2) { Q <- qr.Q(qr(xct_mul(Q))); Q <- qr.Q(qr(xc_mul(Q))) }
      sv <- svd(t(xct_mul(Q)), nu = 0, nv = k)
      V  <- sv$v[, seq_len(k), drop = FALSE]
      ss <- 0
      for (b in split(seq_len(n), ceiling(seq_len(n) / 5000))) ss <- ss + sum(X[b, , drop = FALSE]^2)
      kept <- 100 * sum(sv$d[seq_len(k)]^2) / (ss - n * sum(mu^2))
      # Distances are translation-invariant, so the projection needs no centring.
      X <- X %*% V; rownames(X) <- rownames(gexp_roll_wide)
      log_message(sprintf("initial distances on %d principal components (%.1f%% of smoothed-expression variance)",
                          k, kept), verbose = verbose)
    } else {
      log_message(sprintf("initial distances: exact Euclidean by blocked BLAS (%d cells x %d genes)",
                          nrow(X), ncol(X)), verbose = verbose)
    }
    dist_mat <- gram_dist(X); rm(X); gc()
    log_message("running hclust...")
    hc <- hclust_fun(dist_mat, method = "ward.D2")
    list(gexp_roll_wide = gexp_roll_wide, hc = hc)
  }
  environment(fast_exp_hclust) <- .env
  utils::assignInNamespace("exp_hclust", fast_exp_hclust, ns = "numbat")
}
message(sprintf("[%s] tree building: %s", group_id,
                if (has_fastcluster) "fastcluster::hclust (no cell limit; same result as stats::hclust)"
                else "stats::hclust (limit 65,536 cells)"))
message(sprintf("[%s] initial clustering distances: %s (%d cells)", group_id,
                switch(init_dist, exact = "exact, blocked BLAS",
                       pca = sprintf("top %d principal components", init_pcs),
                       "Numbat's own parallelDist"), ncol(count_mat)))

# init_k: number of starting groups cut from the initial Ward tree (Numbat's
# own default is 3). Set by --numbat_init_k via NUMBAT_INIT_K.
init_k <- suppressWarnings(as.integer(Sys.getenv("NUMBAT_INIT_K", "3")))
if (is.na(init_k) || init_k < 1L) init_k <- 3L
message(sprintf("[%s] init_k = %d starting cell groups", group_id, init_k))

skip_nj_mode <- tolower(Sys.getenv("NUMBAT_SKIP_NJ", "false"))
nj_max_cells <- suppressWarnings(as.integer(Sys.getenv("NUMBAT_NJ_MAX_CELLS", "10000")))
skip_nj <- switch(skip_nj_mode,
  "true" = TRUE, "yes" = TRUE, "1" = TRUE,
  "auto" = ncol(count_mat) > nj_max_cells,
  FALSE)
message(sprintf("[%s] phylogeny seed: %s", group_id,
                if (skip_nj) sprintf("UPGMA only (skip_nj; %d cells%s)", ncol(count_mat),
                                     if (skip_nj_mode == "auto") sprintf(" > %d", nj_max_cells) else "")
                else "better of UPGMA and neighbour-joining"))

# ---- worker count from memory, not CPUs ------------------------------------
# run_numbat forks `ncores` R workers (parallel::mclapply) for the per-cell
# evaluation and NNI, and each fork ends up holding its own copy of what it
# touches. Cirro 2026-11-17 (792f4c3d): with ncores = 48 both DFCI1 (59.5k
# cells) and DFCI7 (79.7k) were OOM-killed at 128 GB in "Evaluating CNV per
# cell", while the main process reported 11.5 / 15.6 GB. The previous run
# (3f4cc38b), 16 workers at 128 GB, got DFCI1 through that step -- the one
# measurement this estimate rests on: about 8 GB per worker at 60k cells,
# scaled linearly by cell count, after reserving a base for the parent.
# The BLAS distance step uses threads (OPENBLAS_NUM_THREADS), not forks, so it
# keeps every CPU. NUMBAT_NCORES (a number) overrides; 'auto' computes.
mem_gb <- suppressWarnings(as.numeric(Sys.getenv("NUMBAT_MEM_GB", "")))
nc_env <- Sys.getenv("NUMBAT_NCORES", "auto")
if (nc_env != "auto" && nzchar(nc_env)) {
  workers <- max(1L, min(ncores, as.integer(nc_env)))
  why <- "set by --numbat_ncores"
} else if (is.finite(mem_gb) && mem_gb > 0) {
  per_worker <- max(2, 8 * ncol(count_mat) / 60000)
  # parent peaked at 11.5 GB (59.5k cells) / 15.6 GB (79.7k): keep ~2x that
  base_gb    <- max(24, 0.4 * ncol(count_mat) / 1000)
  workers <- as.integer(max(1, min(ncores, floor((mem_gb - base_gb) / per_worker))))
  why <- sprintf("%.0f GB task, ~%.1f GB per worker at %d cells, %.0f GB kept for the parent",
                 mem_gb, per_worker, ncol(count_mat), base_gb)
} else {
  workers <- ncores
  why <- "task memory unknown, using every CPU"
}
message(sprintf("[%s] numbat workers: %d of %d CPUs (%s)", group_id, workers, ncores, why))

res <- try(
  numbat::run_numbat(
    count_mat = count_mat,
    lambdas_ref = ref,
    df_allele = df_allele,
    genome = Sys.getenv("NUMBAT_GENOME", "hg38"),
    t = t_param,
    ncores = workers,
    ncores_nni = workers,
    max_iter = max_iter,
    min_LLR = min_llr,
    max_entropy = max_entropy,
    init_k = init_k,
    skip_nj = skip_nj,
    plot = FALSE,
    out_dir = out_sub
  ),
  silent = TRUE
)

if (inherits(res, "try-error")) {
  msg <- as.character(res)
  # Numbat drops cells with no usable SNP coverage BEFORE it checks the two
  # inputs against each other, so a pileup that produced nothing usable surfaces
  # as "No matching cell names between count_mat and df_allele" -- which reads
  # like the namespacing is broken when it is not. The explicit overlap check
  # above already proved the names line up, so reaching here means the coverage
  # is the problem. That is a RESULT for this individual, not a broken run, and
  # killing the task would take the whole cohort's other individuals with it.
  if (grepl("No matching cell names|0 coverage", msg)) {
    emit_empty(outdir, group_id, sprintf(paste(
      "%d cell(s) carried both expression and allele data, but Numbat found no",
      "usable SNP coverage in any of them and dropped them all. The pileup ran",
      "and produced %d allele record(s) across %d librar%s; none fell on a",
      "heterozygous site this individual's reads cover deeply enough. Check the",
      "phasing panel matches --numbat_genome (%s) and that the BAMs are the ones",
      "these barcodes came from."),
      shared, nrow(df_allele), length(libs_seen),
      if (length(libs_seen) == 1L) "y" else "ies",
      Sys.getenv("NUMBAT_GENOME", "hg38")))
  }
  stop(sprintf("[%s] run_numbat failed: %s", group_id, msg))
}

# run_numbat writes its per-cell posterior to clone_post_2.tsv. Later Numbat
# versions may emit a different iteration index, so take the highest available.
posts <- list.files(out_sub, pattern = "^clone_post_\\d+\\.tsv$", full.names = TRUE)
if (!length(posts)) {
  emit_empty(outdir, group_id, paste(
    "Numbat completed but no CNV survived filtering by log-likelihood ratio,",
    "so no clone posterior was produced. Either this individual carries no",
    "detectable clonal copy-number change, or the threshold is too strict for",
    "their coverage -- lower --numbat_min_llr to test the second possibility."))
}
post_path <- posts[order(as.integer(sub(".*clone_post_(\\d+)\\.tsv$", "\\1", posts)))][length(posts)]
message("[", group_id, "] reading ", basename(post_path))

post <- fread(post_path)

cell_col <- intersect(c("cell", "cell_id", "barcode"), names(post))[1]
if (is.na(cell_col)) stop(sprintf("[%s] no cell column in %s", group_id, basename(post_path)))
setnames(post, cell_col, "cell")
post[, cell := as.character(cell)]

# `compartment_opt` is Numbat's own tumour/normal assignment; p_cnv is the
# posterior probability of carrying any CNV. Prefer the explicit assignment and
# fall back to thresholding the posterior.
if ("compartment_opt" %in% names(post)) {
  post[, numbat_call := fifelse(grepl("tumou?r", tolower(compartment_opt)), "malignant",
                        fifelse(grepl("normal", tolower(compartment_opt)), "non_malignant",
                                NA_character_))]
} else if ("p_cnv" %in% names(post)) {
  post[, numbat_call := fifelse(is.na(p_cnv), NA_character_,
                        fifelse(p_cnv > 0.5, "malignant", "non_malignant"))]
} else {
  stop(sprintf("[%s] neither compartment_opt nor p_cnv present", group_id))
}

keep <- intersect(c("cell", "numbat_call", "p_cnv", "clone_opt", "compartment_opt"), names(post))
post <- post[, ..keep]

# Join back onto the pipeline's cell ids via the explicit map. `cell` is the
# namespaced "<library>_<barcode>" that both sides were built from.
# in_reference_lineage comes from prep and is the SAME predicate
# stratify_compartments.qmd gates on. Kept optional so an older cellmap still
# joins; the gated column is then simply absent rather than wrong.
map_cols <- c("cell", "cell_id", "barcode", "library", "group")
has_gate <- "in_reference_lineage" %in% names(cellmap)
if (has_gate) map_cols <- c(map_cols, "in_reference_lineage")
# The compartment the library came from. Optional for the same reason the gate
# is: an older cellmap simply has no such column, and the compartment split is
# then absent rather than wrong. A column present but all-NA counts as absent,
# because a table split into one group called "NA" is worse than no split.
has_stype <- "sample_type" %in% names(cellmap) &&
             any(!is.na(cellmap$sample_type) & nzchar(as.character(cellmap$sample_type)))
if (has_stype) map_cols <- c(map_cols, "sample_type")
res_dt <- merge(cellmap[, ..map_cols], post, by = "cell", all.x = FALSE)
setnames(res_dt, "library", "sample")
if (!nrow(res_dt)) stop(sprintf("[%s] no Numbat cell matched the cell map", group_id))

fwrite(res_dt, file.path(outdir, paste0(group_id, "_numbat_calls.tsv")), sep = "\t")

# ---------------------------------------------------------------------------
# Per-library breakdown. This is the number the per-library configuration got
# wrong and the number a reviewer will look at first: a library that should be
# ~0% neoplastic (PBMC, matched normal) saying otherwise is visible here,
# in .command.log and in a published table, rather than only after someone
# plots the final object.
# ---------------------------------------------------------------------------
by_lib <- res_dt[, .(
  n_cells    = .N,
  malignant  = sum(numbat_call == "malignant", na.rm = TRUE),
  pct_malignant = round(100 * mean(numbat_call == "malignant", na.rm = TRUE), 1)
), by = .(sample)][order(sample)]

# The RAW percentage above is Numbat's compartment_opt, ungated, and for a
# library with no malignant cells it is structurally guaranteed to be non-zero:
# mark_tumor_lineage() always labels the heaviest clade tumour. The percentage
# that survives the lineage gate is what the final compartments actually use.
# Reporting only the first is how a PBMC library comes to look neoplastic in a
# table while the stratified object is fine -- two different numbers, one of
# which nobody had reason to doubt.
if (has_gate) {
  gated <- res_dt[, .(
    n_reference_lineage = sum(in_reference_lineage %in% c(TRUE, "TRUE", "True")),
    malignant_gated = sum(numbat_call == "malignant" &
                          !(in_reference_lineage %in% c(TRUE, "TRUE", "True")), na.rm = TRUE)
  ), by = .(sample)]
  by_lib <- merge(by_lib, gated, by = "sample", all.x = TRUE)
  by_lib[, pct_malignant_gated := round(100 * malignant_gated / pmax(n_cells, 1L), 1)]
  by_lib[, gate_overrode := malignant - malignant_gated]
  setnames(by_lib, "pct_malignant", "pct_malignant_raw")
  setcolorder(by_lib, c("sample", "n_cells", "n_reference_lineage", "malignant",
                        "pct_malignant_raw", "malignant_gated",
                        "pct_malignant_gated", "gate_overrode"))
}
if (has_stype) {
  st_by_lib <- res_dt[, .(sample_type = sample_type[1]), by = .(sample)]
  by_lib <- merge(by_lib, st_by_lib, by = "sample", all.x = TRUE)
  setcolorder(by_lib, c("sample", "sample_type"))
}
by_lib[, reference_basis := ref_basis]
by_lib[, n_reference_columns := if (ref_basis == "internal") ncol(ref) else NA_integer_]
fwrite(by_lib, file.path(outdir, paste0(group_id, "_numbat_by_library.tsv")), sep = "\t")

# ---------------------------------------------------------------------------
# Per-COMPARTMENT breakdown: the reviewer's check, computed.
#
# "In theory there should be almost zero neoplastic cells in PBMCs" is a
# statement about a compartment, not about a library, and until the compartment
# travelled with the cells it could only be answered by joining this table to
# the metadata by hand. Pooling the libraries of one compartment also makes the
# number stabler than any single library's: a 358-cell library and a 10,517-cell
# library carry very different weight, and reading them as peers is its own
# source of alarm.
#
# Both percentages appear here for the same reason they do per library. The raw
# one cannot be ~0 for any compartment -- mark_tumor_lineage() always labels its
# heaviest clade tumour -- so judging a blood compartment by it will always look
# like an overcall. The gated one is what the stratified object uses.
# ---------------------------------------------------------------------------
if (has_stype) {
  by_st <- res_dt[, .(
    n_libraries = uniqueN(sample),
    n_cells     = .N,
    malignant   = sum(numbat_call == "malignant", na.rm = TRUE),
    pct_malignant_raw = round(100 * mean(numbat_call == "malignant", na.rm = TRUE), 1)
  ), by = .(sample_type)][order(sample_type)]
  if (has_gate) {
    gated_st <- res_dt[, .(
      n_reference_lineage = sum(in_reference_lineage %in% c(TRUE, "TRUE", "True")),
      malignant_gated = sum(numbat_call == "malignant" &
                            !(in_reference_lineage %in% c(TRUE, "TRUE", "True")), na.rm = TRUE)
    ), by = .(sample_type)]
    by_st <- merge(by_st, gated_st, by = "sample_type", all.x = TRUE)
    by_st[, pct_malignant_gated := round(100 * malignant_gated / pmax(n_cells, 1L), 1)]
    by_st[, gate_overrode := malignant - malignant_gated]
    setcolorder(by_st, c("sample_type", "n_libraries", "n_cells",
                         "n_reference_lineage", "malignant", "pct_malignant_raw",
                         "malignant_gated", "pct_malignant_gated", "gate_overrode"))
  }
  by_st[, reference_basis := ref_basis]
  fwrite(by_st, file.path(outdir, paste0(group_id, "_numbat_by_sample_type.tsv")), sep = "\t")
}

message(sprintf("[%s] %d cells called | per library:", group_id, nrow(res_dt)))
pct_raw <- if (has_gate) by_lib$pct_malignant_raw else by_lib$pct_malignant
for (i in seq_len(nrow(by_lib))) {
  if (has_gate) {
    message(sprintf("  %-40s %6d cells  %6d malignant raw (%.1f%%)  %6d after lineage gate (%.1f%%)",
                    by_lib$sample[i], by_lib$n_cells[i],
                    by_lib$malignant[i], pct_raw[i],
                    by_lib$malignant_gated[i], by_lib$pct_malignant_gated[i]))
  } else {
    message(sprintf("  %-40s %6d cells  %6d malignant (%.1f%%)",
                    by_lib$sample[i], by_lib$n_cells[i],
                    by_lib$malignant[i], pct_raw[i]))
  }
}
if (has_gate)
  message(sprintf(paste("[%s] 'raw' is Numbat's ungated compartment_opt and is non-zero even",
                        "for a library with no malignant cells; 'after lineage gate' is what",
                        "the stratified compartments use."), group_id))

if (has_stype) {
  message(sprintf("[%s] per sample type:", group_id))
  for (i in seq_len(nrow(by_st))) {
    if (has_gate) {
      message(sprintf("  %-16s %2d librar%s %7d cells  %6d malignant raw (%.1f%%)  %6d after lineage gate (%.1f%%)",
                      by_st$sample_type[i], by_st$n_libraries[i],
                      if (by_st$n_libraries[i] == 1L) "y " else "ies",
                      by_st$n_cells[i], by_st$malignant[i], by_st$pct_malignant_raw[i],
                      by_st$malignant_gated[i], by_st$pct_malignant_gated[i]))
    } else {
      message(sprintf("  %-16s %2d librar%s %7d cells  %6d malignant (%.1f%%)",
                      by_st$sample_type[i], by_st$n_libraries[i],
                      if (by_st$n_libraries[i] == 1L) "y " else "ies",
                      by_st$n_cells[i], by_st$malignant[i], by_st$pct_malignant_raw[i]))
    }
  }
  # Which compartment is supposed to be tumour-free is a biological fact this
  # script has no way to know -- sample_type is a free-text column -- so the
  # comparison is stated rather than judged. Naming the span is enough: a
  # compartment whose gated fraction is far below the others is the expected
  # shape, and one where they are all alike is the finding.
  if (has_gate && nrow(by_st) > 1L) {
    lo <- by_st[which.min(pct_malignant_gated)]
    hi <- by_st[which.max(pct_malignant_gated)]
    message(sprintf(paste("[%s] gated malignant fraction spans %.1f%% (%s) to %.1f%% (%s).",
                          "A compartment expected to be tumour-free should sit at the low",
                          "end; if it does not, raise --numbat_min_llr or lower",
                          "--numbat_max_entropy before trusting the calls."),
                    group_id, lo$pct_malignant_gated, lo$sample_type,
                    hi$pct_malignant_gated, hi$sample_type))
  }
}
