# ============================================================================
# Isotope identification for mass spectrometry imaging
# ----------------------------------------------------------------------------
# Method (after Guo et al.): two m/z features i, j are isotopic if
#     |dm - k * iso_spacing| <= tol(ppm)   AND   image_correlation(i, j) > cutoff
# where dm = mz_j - mz_i and k is an integer isotope order (1, 2, ...).
# Transitivity is resolved by taking connected components of the pass/fail
# graph, so an isotope "group" is a full envelope (M+0, M+1, M+2, ...).
#
# The correlation cutoff is chosen by a MetaSpace-style target-decoy FDR
# (Palmer et al. 2017, Nat. Methods; https://doi.org/10.1038/nmeth.4072):
# real candidates are scored against decoy pairs built from implausible mass
# offsets, and the cutoff is the smallest correlation whose estimated FDR meets
# a target.
#
# The pipeline is exposed as small, composable steps so the process is legible
# from a caller (see 2-Filtering-Aggregation.Rmd):
#   find_isotope_candidates()  - mass-delta candidate pairs (no images).
#   score_pairs()              - per-run, aggregated image correlation per pair.
#   make_decoy_pairs()         - MetaSpace-style implausible decoy pairs.
#   estimate_isotope_fdr()     - target-decoy FDR curve + cutoff.
#   assign_isotope_groups()    - connected components -> per-feature groups.
#
# detect_isotope_relationships() composes those steps (plus optional adducts)
# and annotate_isotopes() wraps a Cardinal MSImagingExperiment. Base R only.
# ============================================================================


# ---- small union-find (disjoint set) --------------------------------------
.uf_new  <- function(n) seq_len(n)
.uf_find <- function(parent, i) {
  i <- as.integer(i)
  root <- i
  while (parent[root] != root) root <- parent[root]
  root
}
.uf_union <- function(parent, a, b) {
  ra <- .uf_find(parent, a); rb <- .uf_find(parent, b)
  if (ra != rb) parent[ra] <- rb
  parent
}
.uf_flatten <- function(parent) vapply(seq_along(parent), function(i) .uf_find(parent, i), integer(1))


# ---- step 1: mass-delta isotope candidates --------------------------------
#' Candidate isotope pairs by mass alone (no ion images).
#'
#' Windowed scan over sorted m/z: for each feature, any heavier feature whose
#' mass gap matches k * iso_spacing within ppm (for some k in iso_orders) is a
#' candidate. Returned in the ORIGINAL feature indexing (i, j = row numbers).
#'
#' @param mz          numeric vector of feature m/z (need not be sorted).
#' @param ppm         mass tolerance in ppm, applied at the anchor (lower) m/z.
#' @param iso_spacing 13C-12C spacing in Da. Default 1.0033548.
#' @param iso_orders  integer isotope orders to test. Default 1:2.
#' @return data.frame(i, j, mz_i, mz_j, delta_m, order), one row per candidate.
find_isotope_candidates <- function(mz, ppm = 20, iso_spacing = 1.0033548,
                                    iso_orders = 1:2) {
  F <- length(mz)
  empty <- data.frame(i = integer(), j = integer(), mz_i = numeric(),
                      mz_j = numeric(), delta_m = numeric(), order = integer())
  if (F < 2) return(empty)

  ord    <- order(mz)
  mzs    <- mz[ord]
  maxwin <- max(iso_orders) * iso_spacing
  out    <- list()
  for (a in seq_len(F - 1L)) {
    mz_a  <- mzs[a]
    tol_a <- ppm * mz_a / 1e6
    b <- a + 1L
    while (b <= F && (mzs[b] - mz_a) <= (maxwin + tol_a)) {
      dm <- mzs[b] - mz_a
      for (k in iso_orders) {
        if (abs(dm - k * iso_spacing) <= tol_a) {
          out[[length(out) + 1L]] <- c(i = ord[a], j = ord[b], delta_m = dm, order = k)
          break
        }
      }
      b <- b + 1L
    }
  }
  if (!length(out)) return(empty)

  cand <- as.data.frame(do.call(rbind, out))
  cand$i     <- as.integer(cand$i)
  cand$j     <- as.integer(cand$j)
  cand$order <- as.integer(cand$order)
  cand$mz_i  <- mz[cand$i]
  cand$mz_j  <- mz[cand$j]
  cand[, c("i", "j", "mz_i", "mz_j", "delta_m", "order")]
}


# ---- step 2: ion-image correlation for a set of pairs ---------------------
#' Aggregated spatial correlation for feature pairs.
#'
#' One correlation per row of `pairs`: correlate the two ion images, optionally
#' within each run then aggregate across runs. This is the single scorer used
#' for both real candidates and decoys.
#'
#' @param pairs     data.frame with integer columns i, j (feature rows).
#' @param intensity F x P matrix (features in ROWS, pixels in COLUMNS).
#' @param run       optional length-P run/sample vector. NULL = one run.
#' @param cor_method   "pearson" (default) or "spearman".
#' @param pixel_mask   "all" (default), "anchor_signal", or "either_signal".
#' @param signal_threshold intensity above which a pixel counts as signal.
#' @param by_run    correlate within each run then aggregate (default TRUE).
#' @param run_agg   combine per-run correlations: "median" (default), "min", "mean".
#' @param min_run_pixels minimum masked pixels for a run to count. Default 20.
#' @param min_runs  minimum evaluable runs required. Default 1.
#' @return numeric vector length nrow(pairs); NA where not evaluable.
score_pairs <- function(pairs, intensity, run = NULL,
                        cor_method = c("pearson", "spearman"),
                        pixel_mask = c("all", "anchor_signal", "either_signal"),
                        signal_threshold = 0,
                        by_run = TRUE,
                        run_agg = c("median", "min", "mean"),
                        min_run_pixels = 20L, min_runs = 1L) {
  cor_method <- match.arg(cor_method)
  pixel_mask <- match.arg(pixel_mask)
  run_agg    <- match.arg(run_agg)

  intensity <- as.matrix(intensity)
  P <- ncol(intensity)
  if (is.null(run)) run <- rep("run1", P)
  run <- as.factor(run)
  if (length(run) != P)
    stop("run must have one entry per pixel: length(run) != ncol(intensity).")
  run_idx <- if (by_run) split(seq_len(P), run) else list(all = seq_len(P))
  agg_fun <- switch(run_agg, min = min, median = stats::median, mean = mean)
  npair   <- nrow(pairs)
  if (!npair) return(numeric(0))

  # fast path: Pearson over all pixels == dot product of per-run z-scored rows.
  # Standardise each feature row ONCE per run, then score every pair vectorised.
  if (pixel_mask == "all" && cor_method == "pearson" && all(is.finite(intensity))) {
    ii <- pairs$i; jj <- pairs$j
    contrib <- matrix(NA_real_, npair, length(run_idx))
    r <- 0L
    for (idx in run_idx) {
      r <- r + 1L
      n <- length(idx)
      if (n < min_run_pixels) next
      M   <- intensity[, idx, drop = FALSE]
      mu  <- rowMeans(M)
      sdv <- sqrt(rowSums((M - mu)^2) / (n - 1))
      Z   <- (M - mu) / sdv                      # rows with sd 0 -> NaN
      cc  <- rowSums(Z[ii, , drop = FALSE] * Z[jj, , drop = FALSE]) / (n - 1)
      cc[!is.finite(sdv[ii]) | sdv[ii] == 0 | !is.finite(sdv[jj]) | sdv[jj] == 0] <- NA_real_
      contrib[, r] <- cc
    }
    return(vapply(seq_len(npair), function(p) {
      v <- contrib[p, ]; v <- v[is.finite(v)]
      if (length(v) < min_runs) NA_real_ else agg_fun(v)
    }, numeric(1)))
  }

  cor1 <- function(a, b) {
    va <- intensity[a, ]; vb <- intensity[b, ]
    cors <- numeric(0)
    for (idx in run_idx) {
      ai <- va[idx]; bi <- vb[idx]
      keep <- switch(pixel_mask,
                     all           = rep(TRUE, length(ai)),
                     anchor_signal = ai > signal_threshold,
                     either_signal = (ai > signal_threshold) | (bi > signal_threshold))
      keep[is.na(keep)] <- FALSE
      ai <- ai[keep]; bi <- bi[keep]
      ok <- is.finite(ai) & is.finite(bi)
      ai <- ai[ok]; bi <- bi[ok]
      if (length(ai) < min_run_pixels) next
      if (stats::sd(ai) == 0 || stats::sd(bi) == 0) next
      cv <- suppressWarnings(stats::cor(ai, bi, method = cor_method))
      if (is.finite(cv)) cors <- c(cors, cv)
    }
    if (length(cors) < min_runs) return(NA_real_)
    agg_fun(cors)
  }

  vapply(seq_len(nrow(pairs)),
         function(r) cor1(pairs$i[r], pairs$j[r]), numeric(1))
}


# ---- step 3: MetaSpace-style decoy pairs ----------------------------------
#' Implausible decoy pairs for target-decoy FDR.
#'
#' Mirrors MetaSpace's implausible-adduct decoys: each candidate anchor is
#' re-paired with a random real feature at an implausible mass offset
#' (|dm| > min_mz_gap, which cannot be an isotope). Repeated `n_rep` times so a
#' median FDR can be taken downstream. One decoy per candidate per repetition.
#'
#' @param candidates data.frame from find_isotope_candidates() (uses column i).
#' @param mz         numeric vector of all feature m/z.
#' @param n_rep      number of decoy repetitions (default 20, MetaSpace-like).
#' @param min_mz_gap minimum |m/z| offset for a decoy partner. Default 3 Da
#'                    (must exceed the isotope window so decoys aren't isotopes).
#' @param seed       RNG seed for reproducibility.
#' @return data.frame(rep, i, j, mz_i, mz_j), one row per candidate per rep.
make_decoy_pairs <- function(candidates, mz, n_rep = 20L, min_mz_gap = 3,
                             seed = 1L) {
  empty <- data.frame(rep = integer(), i = integer(), j = integer(),
                      mz_i = numeric(), mz_j = numeric())
  if (!nrow(candidates)) return(empty)

  set.seed(seed)
  anchors <- candidates$i
  rows <- vector("list", n_rep)
  for (r in seq_len(n_rep)) {
    j <- vapply(anchors, function(i) {
      far <- which(abs(mz - mz[i]) > min_mz_gap)
      if (!length(far)) NA_integer_ else far[sample.int(length(far), 1L)]
    }, integer(1))
    rows[[r]] <- data.frame(rep = r, i = anchors, j = j,
                            mz_i = mz[anchors], mz_j = mz[j])
  }
  do.call(rbind, rows)
}


# ---- step 4: target-decoy FDR curve + cutoff ------------------------------
#' Correlation cutoff from a MetaSpace-style target-decoy FDR.
#'
#' At each candidate cutoff t, FDR(t) = median-over-repetitions(#decoy >= t)
#' divided by (#target >= t). The curve is monotonised (non-increasing in t) and
#' the returned cutoff is the smallest t whose FDR meets `target_fdr`.
#'
#' @param real_corr  target correlations (one per candidate).
#' @param decoy_corr decoy correlations (from score_pairs() on decoy pairs).
#' @param rep        repetition index aligned to decoy_corr (from make_decoy_pairs).
#' @param target_fdr FDR to control. Default 0.10 (MetaSpace default).
#' @param grid       cutoffs to evaluate. Default seq(-1, 1, 0.01).
#' @return list(cutoff, table = data.frame(thr, n_real, n_decoy, fdr)).
estimate_isotope_fdr <- function(real_corr, decoy_corr, rep,
                                 target_fdr = 0.10,
                                 grid = seq(-1, 1, 0.01)) {
  real_corr <- real_corr[is.finite(real_corr)]
  d <- split(decoy_corr, rep)

  n_real  <- vapply(grid, function(t) sum(real_corr >= t), numeric(1))
  n_decoy <- vapply(grid, function(t)
    stats::median(vapply(d, function(dc) sum(dc >= t, na.rm = TRUE), numeric(1))),
    numeric(1))

  fdr      <- ifelse(n_real > 0, n_decoy / n_real, 0)
  fdr_mono <- rev(cummax(rev(fdr)))                 # worst FDR at/above each cutoff
  pass     <- which(fdr_mono <= target_fdr)
  cutoff   <- if (length(pass)) grid[min(pass)] else max(grid)

  list(cutoff = cutoff,
       table  = data.frame(thr = grid, n_real = n_real,
                           n_decoy = n_decoy, fdr = fdr_mono))
}


# ---- step 5: groups from surviving edges ----------------------------------
#' Assign isotope groups from candidate pairs that passed the cutoff.
#'
#' Connected components of the passing pairs; every feature gets a group
#' (unmatched features are singletons). Groups are relabelled 1..G by ascending
#' minimum m/z for stable output.
#'
#' @param edges       passing candidate rows (columns i, j).
#' @param mz          numeric vector of all feature m/z.
#' @param iso_spacing 13C-12C spacing, used to label iso_order. Default 1.0033548.
#' @return data.frame in feature order: feature_id, mz, group_id, iso_order,
#'         is_monoisotopic, group_size, monoisotopic_mz.
assign_isotope_groups <- function(edges, mz, iso_spacing = 1.0033548) {
  F <- length(mz)
  parent <- .uf_new(F)
  for (r in seq_len(nrow(edges)))
    parent <- .uf_union(parent, edges$i[r], edges$j[r])
  comp <- .uf_flatten(parent)

  grp_min <- tapply(mz, comp, min)
  ord     <- order(grp_min)
  relabel <- integer(length(grp_min)); relabel[ord] <- seq_along(ord)
  names(relabel) <- names(grp_min)
  group_id <- relabel[as.character(comp)]

  group_size <- as.integer(table(group_id)[as.character(group_id)])
  mono_mz    <- tapply(mz, group_id, min)[as.character(group_id)]

  data.frame(
    feature_id      = seq_len(F),
    mz              = mz,
    group_id        = as.integer(group_id),
    iso_order       = as.integer(round((mz - mono_mz) / iso_spacing)),
    is_monoisotopic = mz == mono_mz,
    group_size      = group_size,
    monoisotopic_mz = as.numeric(mono_mz),
    row.names       = NULL,
    stringsAsFactors = FALSE
  )
}


# ---- convenience: full pipeline (isotopes + optional adducts) -------------
#' Detect isotope (and optionally adduct) relationships in one call.
#'
#' Composes the steps above with a FIXED correlation cutoff (`cor_threshold`).
#' For a data-driven cutoff, call the steps directly and use estimate_isotope_fdr().
#'
#' @param mz,intensity,run see score_pairs().
#' @param ppm,iso_spacing,iso_orders see find_isotope_candidates().
#' @param cor_threshold correlation cutoff for isotopes. Default 0.7.
#' @param adducts named numeric vector of positive Da offsets between
#'        MONOISOTOPIC peaks, e.g. c("Na/H"=21.981944). NULL disables adducts.
#' @param adduct_ppm,adduct_cor_threshold adduct tolerances (default to isotope ones).
#' @param cor_method,pixel_mask,signal_threshold,by_run,run_agg,min_run_pixels,min_runs
#'        see score_pairs().
#' @param verbose print a summary. Default TRUE.
#' @return list(features, isotope_edges, adduct_edges, params).
detect_isotope_relationships <- function(
    mz, intensity, run = NULL,
    ppm = 10,
    iso_spacing = 1.0033548,
    iso_orders = 1:2,
    cor_threshold = 0.7,
    adducts = c("Na/H" = 21.981944, "K/H" = 37.955882),
    adduct_ppm = ppm,
    adduct_cor_threshold = cor_threshold,
    cor_method = c("pearson", "spearman"),
    pixel_mask = c("all", "anchor_signal", "either_signal"),
    signal_threshold = 0,
    by_run = TRUE,
    run_agg = c("median", "min", "mean"),
    min_run_pixels = 20L,
    min_runs = 1L,
    verbose = TRUE
) {
  cor_method <- match.arg(cor_method)
  pixel_mask <- match.arg(pixel_mask)
  run_agg    <- match.arg(run_agg)

  intensity <- as.matrix(intensity)
  F <- length(mz)
  if (nrow(intensity) != F)
    stop("intensity must have one row per feature: nrow(intensity) != length(mz).")

  score <- function(pairs, method_thr)
    score_pairs(pairs, intensity, run, cor_method = cor_method,
                pixel_mask = pixel_mask, signal_threshold = signal_threshold,
                by_run = by_run, run_agg = run_agg,
                min_run_pixels = min_run_pixels, min_runs = min_runs)

  # 1-2. candidates scored by image correlation
  cand <- find_isotope_candidates(mz, ppm, iso_spacing, iso_orders)
  cand$corr <- if (nrow(cand)) score(cand) else numeric(0)
  if (verbose) message(sprintf("Isotopes: scored %d mass-candidate pair(s).", nrow(cand)))

  iso_edges <- cand[is.finite(cand$corr) & cand$corr > cor_threshold, , drop = FALSE]
  iso_edges <- iso_edges[, c("i", "j", "mz_i", "mz_j", "delta_m", "order", "corr")]

  # 3. groups from passing isotope edges
  grp <- assign_isotope_groups(iso_edges, mz, iso_spacing)

  # 4. adducts between monoisotopic peaks (optional)
  adduct_edges <- data.frame(i = integer(), j = integer(), mz_i = numeric(),
                             mz_j = numeric(), delta_m = numeric(),
                             adduct = character(), corr = numeric())
  if (length(adducts)) {
    mono_idx <- which(grp$is_monoisotopic)
    if (length(mono_idx) >= 2) {
      mo   <- mono_idx[order(mz[mono_idx])]
      mz_m <- mz[mo]
      apairs <- list()
      for (ai in seq_along(mo)) {
        base_mz  <- mz_m[ai]
        tol_base <- adduct_ppm * base_mz / 1e6
        for (nm in names(adducts)) {
          d    <- adducts[[nm]]
          hits <- which(mz_m >= base_mz + d - tol_base & mz_m <= base_mz + d + tol_base)
          for (h in hits)
            apairs[[length(apairs) + 1L]] <-
              data.frame(i = mo[ai], j = mo[h], adduct = nm, stringsAsFactors = FALSE)
        }
      }
      if (length(apairs)) {
        ap <- do.call(rbind, apairs)
        ap$corr <- score(ap)
        ap <- ap[is.finite(ap$corr) & ap$corr > adduct_cor_threshold, , drop = FALSE]
        if (nrow(ap))
          adduct_edges <- data.frame(i = ap$i, j = ap$j, mz_i = mz[ap$i],
                                     mz_j = mz[ap$j], delta_m = mz[ap$j] - mz[ap$i],
                                     adduct = ap$adduct, corr = ap$corr)
      }
    }
  }

  # 5. molecule groups = isotope edges + adduct edges
  parent_mol <- .uf_new(F)
  for (r in seq_len(nrow(iso_edges)))
    parent_mol <- .uf_union(parent_mol, iso_edges$i[r], iso_edges$j[r])
  for (r in seq_len(nrow(adduct_edges)))
    parent_mol <- .uf_union(parent_mol, adduct_edges$i[r], adduct_edges$j[r])
  comp_mol <- .uf_flatten(parent_mol)
  molm     <- tapply(mz, comp_mol, min)
  mrelab   <- integer(length(molm)); mrelab[order(molm)] <- seq_along(molm)
  names(mrelab) <- names(molm)
  molecule_group <- mrelab[as.character(comp_mol)]

  adduct_label <- rep(NA_character_, F)
  for (r in seq_len(nrow(adduct_edges))) {
    j <- adduct_edges$j[r]; lab <- adduct_edges$adduct[r]
    adduct_label[j] <- if (is.na(adduct_label[j])) lab else paste(adduct_label[j], lab, sep = ";")
  }

  features <- data.frame(
    feature         = grp$feature_id,
    mz              = grp$mz,
    isotope_group   = grp$group_id,
    iso_order       = grp$iso_order,
    is_monoisotopic = grp$is_monoisotopic,
    group_size      = grp$group_size,
    monoisotopic_mz = grp$monoisotopic_mz,
    molecule_group  = as.integer(molecule_group),
    adduct          = adduct_label,
    stringsAsFactors = FALSE
  )

  if (verbose) {
    n_multi <- sum(tapply(features$group_size, features$isotope_group, `[`, 1) > 1)
    message(sprintf(
      "Done. %d feature(s) -> %d isotope group(s) (%d multi-peak), %d isotope edge(s), %d adduct link(s).",
      F, length(unique(features$isotope_group)), n_multi,
      nrow(iso_edges), nrow(adduct_edges)))
  }

  list(features = features, isotope_edges = iso_edges,
       adduct_edges = adduct_edges,
       params = list(ppm = ppm, iso_spacing = iso_spacing, iso_orders = iso_orders,
                     cor_threshold = cor_threshold, adducts = adducts,
                     cor_method = cor_method, pixel_mask = pixel_mask,
                     by_run = by_run, run_agg = run_agg))
}


# ---- Cardinal wrapper ------------------------------------------------------
#' Annotate isotope/adduct groups on a Cardinal MSImagingExperiment.
#'
#' Extracts m/z, the intensity matrix and the run factor, runs
#' detect_isotope_relationships(), and writes the per-feature results back into
#' featureData(). Assumes peak-picked / centroided data with a shared m/z axis.
#'
#' @param x          an MSImagingExperiment (Cardinal).
#' @param ...        passed to detect_isotope_relationships() (ppm, cor_threshold, ...).
#' @param materialize coerce spectra(x) into an in-memory matrix (default TRUE).
#' @param return     "object" (default) or "both" (list(object, result)).
annotate_isotopes <- function(x, ..., materialize = TRUE, return = c("object", "both")) {
  return <- match.arg(return)
  if (!methods::is(x, "MSImagingExperiment"))
    stop("x must be a Cardinal MSImagingExperiment.")

  cen <- tryCatch(Cardinal::centroided(x), error = function(e) NA)
  if (isFALSE(cen))
    warning("centroided(x) is FALSE: this method expects peak-picked features. ",
            "Run peakProcess()/peakPick()+peakAlign() first, or results will be meaningless.")

  mzv  <- Cardinal::mz(x)
  runf <- tryCatch(Cardinal::run(x), error = function(e) Cardinal::pData(x)$run)
  runf <- as.factor(runf)

  X <- Cardinal::spectra(x)                 # features x pixels
  if (materialize) X <- as.matrix(X)

  res <- detect_isotope_relationships(mz = mzv, intensity = X, run = runf, ...)

  f <- res$features
  Cardinal::fData(x)$isotope_group      <- f$isotope_group
  Cardinal::fData(x)$iso_order          <- f$iso_order
  Cardinal::fData(x)$is_monoisotopic    <- f$is_monoisotopic
  Cardinal::fData(x)$isotope_group_size <- f$group_size
  Cardinal::fData(x)$monoisotopic_mz    <- f$monoisotopic_mz
  Cardinal::fData(x)$molecule_group     <- f$molecule_group
  Cardinal::fData(x)$adduct             <- f$adduct

  if (return == "both") list(object = x, result = res) else x
}
