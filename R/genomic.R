# funcrop: Genomic and pedigree relationship matrix utilities
#
# Provides functions for constructing and manipulating relationship matrices
# (genomic G, pedigree A, blended H) for use with mixed model fitting in
# multi-environment trial analyses.
#
# Author: Max Moldovan
# Licence: GPL (>= 3)


# ---- make_genomic_matrix -----------------------------------------------------

#' Construct a Genomic Relationship Matrix
#'
#' Computes the genomic relationship matrix (G) from marker data using
#' the VanRaden (2008) Method 1 or Method 2. This is a convenience function;
#' users may also supply pre-computed matrices to fitting functions directly.
#'
#' @param markers Numeric matrix of dimension n_varieties x n_markers.
#'   Coded as 0, 1, 2 (minor allele counts) or -1, 0, 1 (centred). Row names
#'   should identify varieties; column names (if present) identify markers.
#' @param method Character: `"vanraden1"` (default) or `"vanraden2"`.
#'   - **VanRaden Method 1**: \eqn{G = ZZ' / (2 \sum p_j(1-p_j))}, where
#'     \eqn{Z_{ij} = M_{ij} - 2p_j} and \eqn{p_j} is the allele frequency
#'     at marker \eqn{j}.
#'   - **VanRaden Method 2**: \eqn{G = ZDZ' / n_m}, where
#'     \eqn{D = \mathrm{diag}(1/(2 p_j (1-p_j)))}.
#' @param min_maf Numeric: minimum minor allele frequency filter. Markers with
#'   MAF below this threshold are removed. Default 0.01.
#' @param scale Logical: if `TRUE` (default), scale G so that the mean
#'   diagonal equals 1.
#'
#' @return A symmetric positive (semi-)definite matrix of dimension
#'   n_varieties x n_varieties, with row/colnames from `rownames(markers)`.
#'
#' @details
#' Missing marker values (`NA`) are imputed to the column mean (i.e., twice
#' the allele frequency) prior to computation. If many markers are missing,
#' consider using a dedicated imputation tool first.
#'
#' @references
#' VanRaden, P.M. (2008). Efficient methods to compute genomic predictions.
#' *Journal of Dairy Science*, 91(11), 4414--4423.
#' \doi{10.3168/jds.2007-0980}
#'
#' @seealso [make_pedigree_matrix()], [make_H_matrix()],
#'   [check_relationship_matrix()]
#'
#' @examples
#' set.seed(123)
#' # Simulate 20 varieties, 500 markers (coded 0/1/2)
#' M <- matrix(sample(0:2, 20 * 500, replace = TRUE,
#'                     prob = c(0.25, 0.50, 0.25)),
#'             nrow = 20, ncol = 500)
#' rownames(M) <- paste0("var_", seq_len(20))
#'
#' G <- make_genomic_matrix(M, method = "vanraden1")
#' dim(G)
#' isSymmetric(G)
#'
#' @export
make_genomic_matrix <- function(markers,
                                method = c("vanraden1", "vanraden2"),
                                min_maf = 0.01,
                                scale = TRUE) {

  # --- Input validation -------------------------------------------------------
  method <- match.arg(method)

  if (!is.matrix(markers) || !is.numeric(markers)) {
    stop("`markers` must be a numeric matrix.", call. = FALSE)
  }
  n <- nrow(markers)
  m <- ncol(markers)

  if (n < 2L) {
    stop("`markers` must have at least 2 rows (varieties).", call. = FALSE)
  }
  if (m < 1L) {
    stop("`markers` must have at least 1 column (marker).", call. = FALSE)
  }

  if (!is.numeric(min_maf) || length(min_maf) != 1L ||
      min_maf < 0 || min_maf > 0.5) {
    stop("`min_maf` must be a single numeric value in [0, 0.5].", call. = FALSE)
  }

  # --- Impute missing values to column means ----------------------------------
  na_count <- sum(is.na(markers))
  if (na_count > 0L) {
    .msg(sprintf("Imputing %d missing marker values to column means.", na_count))
    col_means <- colMeans(markers, na.rm = TRUE)
    for (j in seq_len(m)) {
      na_idx <- which(is.na(markers[, j]))
      if (length(na_idx) > 0L) {
        markers[na_idx, j] <- col_means[j]
      }
    }
  }

  # --- Compute allele frequencies and filter by MAF ---------------------------
  # Detect coding: if min is negative, assume centred (-1, 0, 1)
  if (min(markers) < -0.5) {
    # Centred coding: p = (mean + 1) / 2
    p <- (colMeans(markers) + 1) / 2
  } else {
    # Dosage coding (0, 1, 2): p = mean / 2
    p <- colMeans(markers) / 2
  }

  # Minor allele frequency
  maf <- pmin(p, 1 - p)
  keep <- which(maf >= min_maf & maf > 0)

  if (length(keep) == 0L) {
    stop("No markers remaining after MAF filter (min_maf = ", min_maf, ").",
         call. = FALSE)
  }

  n_removed <- m - length(keep)
  if (n_removed > 0L) {
    .msg(sprintf("Removed %d markers with MAF < %.4f. %d markers retained.",
                 n_removed, min_maf, length(keep)))
  }

  markers <- markers[, keep, drop = FALSE]
  p <- p[keep]
  m <- length(keep)

  # --- Centre marker matrix ---------------------------------------------------
  # Z_{ij} = M_{ij} - 2p_j  (regardless of original coding)
  if (min(markers) < -0.5) {
    # Already centred around 0; shift to dosage then centre
    Z <- markers - rep(2 * p - 1, each = n)
  } else {
    Z <- markers - rep(2 * p, each = n)
  }

  # --- Compute G --------------------------------------------------------------
  if (method == "vanraden1") {
    # G = ZZ' / (2 * sum(p * (1 - p)))
    denom <- 2 * sum(p * (1 - p))
    if (denom < .Machine$double.eps) {
      stop("Denominator is zero: all markers are monomorphic.", call. = FALSE)
    }
    G <- tcrossprod(Z) / denom

  } else {
    # VanRaden Method 2: G = Z D Z' / m
    # D = diag(1 / (2 * p * (1 - p)))
    d_vals <- 1 / (2 * p * (1 - p))
    # Z_scaled: multiply each column j of Z by sqrt(d_vals[j])
    Z_scaled <- Z * rep(sqrt(d_vals), each = n)
    G <- tcrossprod(Z_scaled) / m
  }

  # --- Scale to mean diagonal = 1 --------------------------------------------
  if (isTRUE(scale)) {
    mean_diag <- mean(diag(G))
    if (mean_diag > .Machine$double.eps) {
      G <- G / mean_diag
    }
  }

  # --- Set names --------------------------------------------------------------
  if (!is.null(rownames(markers))) {
    # markers may have been subsetted; use original row names
    # (rownames persist through column subsetting)
  }
  rn <- rownames(markers)
  if (is.null(rn)) {
    rn <- paste0("V", seq_len(n))
  }
  rownames(G) <- rn
  colnames(G) <- rn

  G
}


# ---- make_pedigree_matrix ----------------------------------------------------

#' Construct a Pedigree-Based Relationship Matrix (A Matrix)
#'
#' Computes the numerator relationship matrix (A) from a pedigree data.frame
#' using Henderson's tabular method. The A matrix encodes the expected
#' proportion of alleles shared identical by descent between all pairs of
#' individuals.
#'
#' @param pedigree A data.frame (or data.table) with exactly three columns
#'   in order: `id`, `sire`, `dam`. Unknown parents should be coded as `NA`,
#'   `"0"`, or `0`. The pedigree must be sorted such that parents appear
#'   before their offspring (chronological order).
#'
#' @return A symmetric positive definite matrix of dimension n x n, where
#'   n is the number of individuals. Row and column names are set to the
#'   individual IDs.
#'
#' @details
#' Implements Henderson's (1976) tabular method without external dependencies.
#' Time complexity is O(n^2) in the number of individuals, which is adequate
#' for typical crop trial pedigrees (up to ~5000 individuals). For larger
#' pedigrees, consider using the `pedigreemm` or `nadiv` packages.
#'
#' Inbreeding coefficients are computed on the diagonal: \eqn{A_{ii} = 1 + F_i},
#' where \eqn{F_i} is the inbreeding coefficient of individual \eqn{i}.
#'
#' @references
#' Henderson, C.R. (1976). A simple method for computing the inverse of a
#' numerator relationship matrix used in prediction of breeding values.
#' *Biometrics*, 32(1), 69--83. \doi{10.2307/2529339}
#'
#' @seealso [make_genomic_matrix()], [make_H_matrix()],
#'   [check_relationship_matrix()]
#'
#' @examples
#' ped <- data.frame(
#'   id   = c("A", "B", "C", "D", "E"),
#'   sire = c(NA,  NA,  "A", "A", "C"),
#'   dam  = c(NA,  NA,  "B", "B", "D")
#' )
#' A <- make_pedigree_matrix(ped)
#' print(round(A, 3))
#'
#' @export
make_pedigree_matrix <- function(pedigree) {

  # --- Input validation -------------------------------------------------------
  if (!is.data.frame(pedigree)) {
    stop("`pedigree` must be a data.frame or data.table.", call. = FALSE)
  }

  if (ncol(pedigree) < 3L) {
    stop("`pedigree` must have at least 3 columns: id, sire, dam.",
         call. = FALSE)
  }

  # Use first three columns
  ped <- data.table::as.data.table(pedigree)[, 1:3]
  data.table::setnames(ped, c("id", "sire", "dam"))

  # Coerce to character
  ped[, id   := as.character(id)]
  ped[, sire := as.character(sire)]
  ped[, dam  := as.character(dam)]

  # Recode unknown parents
  ped[sire %in% c("0", "NA", ""), sire := NA_character_]
  ped[dam  %in% c("0", "NA", ""), dam  := NA_character_]
  ped[is.na(sire), sire := NA_character_]
  ped[is.na(dam),  dam  := NA_character_]

  # Check for duplicates
  if (anyDuplicated(ped$id)) {
    stop("Duplicate IDs found in pedigree.", call. = FALSE)
  }

  n <- nrow(ped)
  ids <- ped$id

  # Map parents to integer indices (NA if unknown)
  sire_idx <- match(ped$sire, ids)
  dam_idx  <- match(ped$dam, ids)

  # Validate parent ordering: parents must appear before offspring
  for (i in seq_len(n)) {
    if (!is.na(sire_idx[i]) && sire_idx[i] >= i) {
      stop(sprintf(
        "Pedigree not in chronological order: sire '%s' of '%s' must appear earlier.",
        ped$sire[i], ids[i]
      ), call. = FALSE)
    }
    if (!is.na(dam_idx[i]) && dam_idx[i] >= i) {
      stop(sprintf(
        "Pedigree not in chronological order: dam '%s' of '%s' must appear earlier.",
        ped$dam[i], ids[i]
      ), call. = FALSE)
    }
  }

  # --- Henderson's tabular method ---------------------------------------------
  A <- matrix(0, nrow = n, ncol = n)

  for (i in seq_len(n)) {
    s <- sire_idx[i]
    d <- dam_idx[i]

    # Diagonal: A[i,i] = 1 + F_i, where F_i = 0.5 * A[s,d] if both known
    if (!is.na(s) && !is.na(d)) {
      A[i, i] <- 1 + 0.5 * A[s, d]
    } else {
      A[i, i] <- 1
    }

    # Off-diagonals: fill row i and column i for individuals j > i
    # (done in subsequent iterations when j uses i as parent)
    # For individuals j < i: A[i,j] = 0.5*(A[s,j] + A[d,j])
    if (i > 1L) {
      for (j in seq_len(i - 1L)) {
        val <- 0
        if (!is.na(s)) val <- val + 0.5 * A[s, j]
        if (!is.na(d)) val <- val + 0.5 * A[d, j]
        A[i, j] <- val
        A[j, i] <- val
      }
    }
  }

  rownames(A) <- ids
  colnames(A) <- ids

  A
}


# ---- make_H_matrix -----------------------------------------------------------

#' Blend Genomic and Pedigree Relationship Matrices (H Matrix)
#'
#' Creates the H-matrix for single-step GBLUP by blending a genomic
#' relationship matrix (G) with a pedigree-based numerator relationship
#' matrix (A), following the Legarra et al. (2009) formulation.
#'
#' @param G Genomic relationship matrix (square, symmetric) for genotyped
#'   individuals. Must have row/colnames identifying genotyped varieties.
#' @param A Pedigree relationship matrix (square, symmetric) for all
#'   individuals (genotyped and non-genotyped). Must have row/colnames
#'   identifying all varieties.
#' @param genotyped Character vector of genotyped individual IDs. Must be a
#'   subset of `rownames(A)` and must match `rownames(G)`.
#' @param tau Numeric: weight for \eqn{A_{22}^{-1}} in blending. Default 1.
#' @param omega Numeric: weight for \eqn{G^{-1}} in blending. Default 1.
#'
#' @return The H matrix of dimension n_all x n_all, with row/colnames from A.
#'
#' @details
#' The inverse of H is constructed as:
#' \deqn{H^{-1} = A^{-1} + \begin{pmatrix} 0 & 0 \\ 0 &
#'   \omega G^{-1} - \tau A_{22}^{-1} \end{pmatrix}}
#'
#' where \eqn{A_{22}} is the submatrix of A corresponding to genotyped
#' individuals. H itself is then obtained by inverting \eqn{H^{-1}}.
#'
#' Setting \eqn{\tau = \omega = 1} gives the standard Legarra formulation.
#' Values different from 1 can be used for regularisation or when G and A
#' are on different scales.
#'
#' @references
#' Legarra, A., Aguilar, I. and Misztal, I. (2009). A relationship matrix
#' including full pedigree and genomic information. *Journal of Dairy Science*,
#' 92(9), 4656--4663. \doi{10.3168/jds.2009-2061}
#'
#' Aguilar, I., Misztal, I., Johnson, D.L., Legarra, A., Tsuruta, S. and
#' Lawlor, T.J. (2010). Hot topic: A unified approach to utilize phenotypic,
#' full pedigree, and genomic information for genetic evaluation of Holstein
#' final score. *Journal of Dairy Science*, 93(2), 743--752.
#'
#' @seealso [make_genomic_matrix()], [make_pedigree_matrix()],
#'   [check_relationship_matrix()]
#'
#' @examples
#' \dontrun{
#' # Suppose A is a 50x50 pedigree matrix and G is 20x20 genomic matrix
#' # for the first 20 individuals
#' H <- make_H_matrix(G, A, genotyped = rownames(G))
#' dim(H)  # 50 x 50
#' }
#'
#' @export
make_H_matrix <- function(G, A, genotyped, tau = 1, omega = 1) {

  # --- Input validation -------------------------------------------------------
  if (!is.matrix(G) || !is.numeric(G)) {
    stop("`G` must be a numeric matrix.", call. = FALSE)
  }
  if (!is.matrix(A) || !is.numeric(A)) {
    stop("`A` must be a numeric matrix.", call. = FALSE)
  }
  if (nrow(G) != ncol(G)) {
    stop("`G` must be square.", call. = FALSE)
  }
  if (nrow(A) != ncol(A)) {
    stop("`A` must be square.", call. = FALSE)
  }
  if (!is.character(genotyped) || length(genotyped) < 1L) {
    stop("`genotyped` must be a character vector of genotyped IDs.",
         call. = FALSE)
  }

  # Validate names
  if (is.null(rownames(G)) || is.null(rownames(A))) {
    stop("Both `G` and `A` must have row/colnames identifying individuals.",
         call. = FALSE)
  }

  # Check genotyped IDs match G
  if (!setequal(genotyped, rownames(G))) {
    stop("`genotyped` must match `rownames(G)` exactly.", call. = FALSE)
  }

  # Check genotyped IDs are subset of A
  missing_in_a <- setdiff(genotyped, rownames(A))
  if (length(missing_in_a) > 0L) {
    stop(sprintf(
      "%d genotyped IDs not found in `A`: %s",
      length(missing_in_a),
      paste(head(missing_in_a, 5L), collapse = ", ")
    ), call. = FALSE)
  }

  if (!is.numeric(tau) || length(tau) != 1L || tau <= 0) {
    stop("`tau` must be a positive numeric scalar.", call. = FALSE)
  }
  if (!is.numeric(omega) || length(omega) != 1L || omega <= 0) {
    stop("`omega` must be a positive numeric scalar.", call. = FALSE)
  }

  all_ids <- rownames(A)
  n_all <- length(all_ids)
  n_g   <- length(genotyped)

  # --- Reorder G to match A's ordering of genotyped IDs -----------------------
  g_order <- match(genotyped, all_ids)
  g_idx_in_a <- sort(g_order)  # positions in A of genotyped individuals
  genotyped_ordered <- all_ids[g_idx_in_a]

  G <- G[genotyped_ordered, genotyped_ordered, drop = FALSE]

  # --- Extract A22 (submatrix for genotyped individuals) ----------------------
  A22 <- A[g_idx_in_a, g_idx_in_a, drop = FALSE]

  # --- Compute inverses -------------------------------------------------------
  A_inv <- tryCatch(
    solve(A),
    error = function(e) {
      stop("Failed to invert A matrix: ", conditionMessage(e),
           "\nConsider using check_relationship_matrix() first.", call. = FALSE)
    }
  )

  G_inv <- tryCatch(
    solve(G),
    error = function(e) {
      stop("Failed to invert G matrix: ", conditionMessage(e),
           "\nConsider using check_relationship_matrix(G) first.", call. = FALSE)
    }
  )

  A22_inv <- tryCatch(
    solve(A22),
    error = function(e) {
      stop("Failed to invert A22 submatrix: ", conditionMessage(e),
           call. = FALSE)
    }
  )

  # --- Construct H^{-1} = A^{-1} + [0, 0; 0, omega*G^{-1} - tau*A22^{-1}] --
  H_inv <- A_inv
  delta <- omega * G_inv - tau * A22_inv
  H_inv[g_idx_in_a, g_idx_in_a] <- H_inv[g_idx_in_a, g_idx_in_a] + delta

  # --- Invert to get H --------------------------------------------------------
  H <- tryCatch(
    solve(H_inv),
    error = function(e) {
      stop("Failed to invert H^{-1}: ", conditionMessage(e), call. = FALSE)
    }
  )

  # Ensure symmetry (numerical)
  H <- 0.5 * (H + t(H))

  rownames(H) <- all_ids
  colnames(H) <- all_ids

  H
}


# ---- check_relationship_matrix -----------------------------------------------

#' Check and Repair a Relationship Matrix
#'
#' Ensures a relationship matrix is symmetric, positive definite (bending if
#' necessary), and properly scaled. Useful as a pre-processing step before
#' passing relationship matrices to mixed model fitting functions.
#'
#' @param K Numeric square matrix (relationship matrix).
#' @param tol Numeric: eigenvalue tolerance. Eigenvalues below this value are
#'   raised to `tol` when bending. Default `1e-6`.
#' @param bend Logical: if `TRUE` (default), bend the matrix to ensure positive
#'   definiteness by replacing negative/near-zero eigenvalues with `tol`.
#'   If `FALSE`, only checks and reports issues without modifying.
#'
#' @return If `bend = TRUE`, returns the repaired symmetric positive definite
#'   matrix (with row/colnames preserved). If `bend = FALSE`, returns the
#'   symmetrised matrix (or raises an error if not PD).
#'
#' @details
#' The bending procedure uses spectral decomposition:
#' 1. Decompose \eqn{K = U \Lambda U'}.
#' 2. Replace any \eqn{\lambda_i < \mathrm{tol}} with \eqn{\mathrm{tol}}.
#' 3. Reconstruct \eqn{K^* = U \Lambda^* U'}.
#'
#' This is the simplest form of matrix bending (Jorjani et al., 2003). For
#' more sophisticated approaches (e.g., preserving diagonal elements), see
#' the `Matrix` or `nadiv` packages.
#'
#' @references
#' Jorjani, H., Klei, L. and Emanuelson, U. (2003). A simple method for
#' weighted bending of genetic (co)variance matrices. *Journal of Dairy
#' Science*, 86(2), 677--679. \doi{10.3168/jds.S0022-0302(03)73646-7}
#'
#' @seealso [make_genomic_matrix()], [make_pedigree_matrix()],
#'   [make_H_matrix()]
#'
#' @examples
#' # Create a near-singular matrix
#' set.seed(42)
#' K <- crossprod(matrix(rnorm(30), 6, 5))
#' K[1, 1] <- 0  # Make it singular
#' K_fixed <- check_relationship_matrix(K)
#' eigen(K_fixed)$values  # All positive
#'
#' @export
check_relationship_matrix <- function(K, tol = 1e-6, bend = TRUE) {

  # --- Input validation -------------------------------------------------------
  if (!is.matrix(K) || !is.numeric(K)) {
    stop("`K` must be a numeric matrix.", call. = FALSE)
  }
  if (nrow(K) != ncol(K)) {
    stop("`K` must be a square matrix.", call. = FALSE)
  }
  if (!is.numeric(tol) || length(tol) != 1L || tol <= 0) {
    stop("`tol` must be a positive numeric scalar.", call. = FALSE)
  }

  n <- nrow(K)
  rn <- rownames(K)
  cn <- colnames(K)

  # --- Enforce symmetry -------------------------------------------------------
  max_asym <- max(abs(K - t(K)))
  if (max_asym > .Machine$double.eps * 100) {
    .msg(sprintf("Matrix asymmetry detected (max |K - K'| = %.2e). Symmetrising.",
                 max_asym))
  }
  K <- 0.5 * (K + t(K))

  # --- Check positive definiteness --------------------------------------------
  eig <- eigen(K, symmetric = TRUE)
  min_eval <- min(eig$values)
  n_negative <- sum(eig$values < tol)

  if (n_negative > 0L) {
    .msg(sprintf(
      "Matrix has %d eigenvalue(s) below tolerance (min = %.4e).",
      n_negative, min_eval
    ))

    if (!isTRUE(bend)) {
      stop(sprintf(
        "Relationship matrix is not positive definite (min eigenvalue = %.4e). ",
        min_eval,
        "Set `bend = TRUE` to repair."
      ), call. = FALSE)
    }

    # Bend: replace small eigenvalues with tol
    eig$values[eig$values < tol] <- tol
    K <- eig$vectors %*% diag(eig$values) %*% t(eig$vectors)
    K <- 0.5 * (K + t(K))  # re-symmetrise after reconstruction

    .msg("Matrix bent to positive definiteness.")
  } else {
    .msg(sprintf("Matrix is positive definite (min eigenvalue = %.4e).",
                 min_eval))
  }

  # --- Restore names ----------------------------------------------------------
  rownames(K) <- rn
  colnames(K) <- cn

  K
}
