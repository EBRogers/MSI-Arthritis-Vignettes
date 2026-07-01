connected_components <- function(nb) {
	n <- length(nb)
	comp <- integer(n)
	cid <- 0L
	for (i in seq_len(n)) {
		if (comp[[i]] != 0L) next
		cid <- cid + 1L
		stack <- integer(0L)
		stack[[1L]] <- i
		sp <- 1L
		comp[[i]] <- cid
		while (sp > 0L) {
			v <- stack[[sp]]
			sp <- sp - 1L
			nei <- nb[[v]]
			if (length(nei) == 0L) next
			for (u in nei) {
				if (comp[[u]] == 0L) {
					comp[[u]] <- cid
					sp <- sp + 1L
					stack[[sp]] <- u
				}
			}
		}
	}
	comp
}

pixel_islands_mask_xy <- function(x, y = NULL, n = 1L, group = NULL, nb = NULL, include_diagonal = TRUE) {
	if (is.data.frame(x)) {
		if (!all(c("x", "y") %in% colnames(x))) stop("x data.frame must contain columns named 'x' and 'y'")
		if (is.null(y)) y <- x[["y"]]
		x <- x[["x"]]
	}
	if (is.null(y)) stop("y must be provided")
	if (length(x) != length(y)) stop("x and y must have the same length")

	n <- as.integer(n)
	if (length(n) != 1L || is.na(n) || n < 0L) stop("n must be a single non-negative integer")

	npix <- length(x)
	if (npix < 1L) return(logical(0L))

	x <- as.numeric(x)
	y <- as.numeric(y)
	if (anyNA(x) || anyNA(y)) stop("x and y must not contain NA")

	if (is.null(group)) {
		group <- rep.int("1", npix)
	} else {
		if (length(group) != npix) stop("group must have length equal to length(x)")
		if (anyNA(group)) stop("group contains NA")
		group <- as.character(group)
	}

	if (is.null(nb)) {
		tol <- sqrt(.Machine$double.eps)
		if (any(abs(x - round(x)) > tol) || any(abs(y - round(y)) > tol)) {
			stop("x and y must be integer-valued pixel coordinates when nb is NULL")
		}
		x <- as.integer(round(x))
		y <- as.integer(round(y))
	} else {
		if (!is.list(nb)) stop("nb must be a list when provided")
		if (length(nb) != npix) stop("length(nb) must equal number of pixels")
		for (i in seq_len(npix)) {
			nei <- nb[[i]]
			if (length(nei) == 0L) next
			if (anyNA(nei)) stop("nb contains NA at index ", i)
			if (any(nei < 1L | nei > npix)) stop("nb contains out-of-range indices at index ", i)
		}
	}

	if (!isTRUE(include_diagonal) && !identical(include_diagonal, FALSE)) {
		stop("include_diagonal must be TRUE or FALSE")
	}

	off <- if (isTRUE(include_diagonal)) {
		expand.grid(dx = -1L:1L, dy = -1L:1L)
	} else {
		data.frame(dx = c(-1L, 1L, 0L, 0L), dy = c(0L, 0L, -1L, 1L))
	}
	if (nrow(off) > 0L) off <- off[!(off$dx == 0L & off$dy == 0L), , drop = FALSE]

	out <- logical(npix)
	idx_by_group <- split(seq_len(npix), group)
	for (gn in names(idx_by_group)) {
		idx <- idx_by_group[[gn]]
		if (length(idx) == 0L) next

		nb_local <- vector("list", length(idx))
		if (is.null(nb)) {
			xi <- x[idx]
			yi <- y[idx]
			keys <- paste(xi, yi, sep = ":")
			for (j in seq_along(idx)) {
				nb_keys <- paste(xi[[j]] + off$dx, yi[[j]] + off$dy, sep = ":")
				loc <- match(nb_keys, keys, nomatch = 0L)
				nb_local[[j]] <- as.integer(loc[loc > 0L])
			}
		} else {
			local_map <- integer(npix)
			local_map[idx] <- seq_along(idx)
			for (j in seq_along(idx)) {
				gidx <- idx[[j]]
				nei <- as.integer(nb[[gidx]])
				if (length(nei) == 0L) {
					nb_local[[j]] <- integer(0L)
					next
				}
				nei <- nei[nei != gidx]
				loc <- local_map[nei]
				loc <- loc[loc > 0L]
				nb_local[[j]] <- unique(as.integer(loc))
			}
		}

		comp <- connected_components(nb_local)
		sz <- table(comp)
		out[idx] <- as.integer(sz[as.character(comp)]) > n
	}

	nms <- names(x)
	if (!is.null(nms) && length(nms) == length(out)) names(out) <- nms
	out
}

drop_pixel_islands_mask_mse <- function(mse, n = 1L, r = sqrt(2), metric = "maximum", group = Cardinal::run(mse)) {
	if (!is(mse, "SpectralImagingData")) stop("mse must inherit from SpectralImagingData")
	n <- as.integer(n)
	if (length(n) != 1L || is.na(n) || n < 0L) stop("n must be a single non-negative integer")
	r <- as.numeric(r)
	if (length(r) != 1L || !is.finite(r) || r < 0) stop("r must be a single non-negative number")

	npix <- length(mse)
	if (npix < 1L) return(logical(0L))

	if (length(group) != npix) stop("group must have length equal to number of pixels in mse")
	if (anyNA(group)) stop("group contains NA")

	nb <- Cardinal::findNeighbors(mse, r = r, metric = metric)
	coord_df <- Cardinal::coord(mse)
	if (!is.data.frame(coord_df)) coord_df <- as.data.frame(coord_df)
	if (!all(c("x", "y") %in% colnames(coord_df))) stop("coord(mse) must contain numeric columns named 'x' and 'y'")

	keep <- pixel_islands_mask_xy(
		x = coord_df[["x"]],
		y = coord_df[["y"]],
		n = n,
		group = as.character(group),
		nb = nb,
		include_diagonal = TRUE
	)

	nms <- try(Cardinal::pixelNames(mse), silent = TRUE)
	if (!inherits(nms, "try-error") && length(nms) == length(keep)) names(keep) <- nms
	keep
}

near_edge_xy <- function(x, y = NULL, radius = 1L, group = NULL, contiguous = FALSE) {
	if (is.data.frame(x)) {
		if (!all(c("x", "y") %in% colnames(x))) stop("x data.frame must contain columns named 'x' and 'y'")
		y <- x[["y"]]
		x <- x[["x"]]
	}
	if (is.null(y)) stop("y must be provided")
	if (length(x) != length(y)) stop("x and y must have the same length")
	if (length(x) < 1L) return(factor(logical(0L), levels = c(FALSE, TRUE)))

	r <- as.numeric(radius)
	if (length(r) != 1L || !is.finite(r) || r < 0) stop("radius must be a single non-negative number")
	if (abs(r - round(r)) > 0) stop("radius must be an integer in pixel-coordinate units")
	r <- as.integer(round(r))

	x <- as.numeric(x)
	y <- as.numeric(y)
	if (anyNA(x) || anyNA(y)) stop("x and y must not contain NA")

	if (is.null(group)) {
		group <- rep.int("1", length(x))
	} else {
		if (length(group) != length(x)) stop("group must have length equal to length(x)")
		if (anyNA(group)) stop("group contains NA")
		group <- as.character(group)
	}

	out <- logical(length(x))
	idx_by_group <- split(seq_along(x), group)
	for (gn in names(idx_by_group)) {
		idx <- idx_by_group[[gn]]
		if (length(idx) == 0L) next
		xi <- x[idx]
		yi <- y[idx]
		keys_all <- paste(xi, yi, sep = ":")

		minx <- ave(xi, yi, FUN = min)
		maxx <- ave(xi, yi, FUN = max)
		miny <- ave(yi, xi, FUN = min)
		maxy <- ave(yi, xi, FUN = max)

		ed <- (xi == minx) | (xi == maxx) | (yi == miny) | (yi == maxy)
		if (isTRUE(contiguous) && any(ed)) {
			gi <- which(ed)
			ex <- xi[gi]
			ey <- yi[gi]
			outer <- (ex == min(xi)) | (ex == max(xi)) | (ey == min(yi)) | (ey == max(yi))
			if (any(outer)) {
				ek <- paste(ex, ey, sep = ":")
				off <- expand.grid(dx = -1:1, dy = -1:1)
				off <- off[!(off$dx == 0 & off$dy == 0), , drop = FALSE]
				nb <- vector("list", length(ek))
				if (length(ek) > 0L) {
					from_all <- integer(0L)
					to_all <- integer(0L)
					for (k in seq_len(nrow(off))) {
						k2 <- paste(ex + off$dx[[k]], ey + off$dy[[k]], sep = ":")
						m <- match(k2, ek, nomatch = NA_integer_)
						ok <- which(!is.na(m))
						if (length(ok) > 0L) {
							from_all <- c(from_all, ok)
							to_all <- c(to_all, m[ok])
						}
					}
					if (length(from_all) > 0L) {
						nb2 <- split(to_all, from_all)
						idx2 <- as.integer(names(nb2))
						for (ii in seq_along(idx2)) nb[[idx2[[ii]]]] <- unname(nb2[[ii]])
					}
				}
				cc <- connected_components(nb)
				keep_cc <- unique(cc[outer])
				keep <- cc %in% keep_cc
				ed[] <- FALSE
				ed[gi[keep]] <- TRUE
			}
		}

		if (!any(ed)) {
			out[idx] <- FALSE
			next
		}
		if (r == 0L) {
			near <- ed
		} else {
			xe <- xi[ed]
			ye <- yi[ed]
			off <- expand.grid(dx = (-r):r, dy = (-r):r)
			near_keys <- character(0L)
			need <- length(xe) * nrow(off)
			if (need > 0L) {
				near_keys <- character(need)
				pos <- 1L
				for (k in seq_len(nrow(off))) {
					k2 <- paste(xe + off$dx[[k]], ye + off$dy[[k]], sep = ":")
					end <- pos + length(k2) - 1L
					near_keys[pos:end] <- k2
					pos <- end + 1L
				}
				near_keys <- unique(near_keys)
			}
			near <- match(keys_all, near_keys, nomatch = 0L) != 0L
		}
		out[idx] <- near
	}

	factor(out, levels = c(FALSE, TRUE))
}

add_edge_factor_mse <- function(mse, radius = 1L, contiguous = FALSE) {
	if (!is(mse, "SpectralImagingData")) stop("mse must inherit from SpectralImagingData")
	if (length(mse) < 1L) {
		mse$edge <- factor(logical(0L), levels = c(FALSE, TRUE))
		return(mse)
	}
	cn <- Cardinal::coordNames(mse)
	if (!all(c("x", "y") %in% cn)) stop("coord(mse) must contain numeric columns named 'x' and 'y'")
	coord_df <- Cardinal::coord(mse)
	if (!is.data.frame(coord_df)) coord_df <- as.data.frame(coord_df)
	x <- as.numeric(coord_df[["x"]])
	y <- as.numeric(coord_df[["y"]])
	if (anyNA(x) || anyNA(y)) stop("coord(mse)$x and coord(mse)$y must not contain NA")
	runv <- Cardinal::run(mse)
	if (anyNA(runv)) stop("run(mse) contains NA")
	edge <- near_edge_xy(x = x, y = y, radius = radius, group = runv, contiguous = contiguous)
	mse$edge <- edge
	mse
}