#' @keywords internal
randomized.svd <- function(A, K, q, method = 'rsvd', mkl.seed = -1) {
    out <- stats::setNames(vector("list", 3), c("u", "d", "v"))
    if (method == 'rsvd') {
        out <- rsvd::rsvd(A, K, q = q)
    } else if (method == 'rsvd-mkl') {
        rlang::check_installed("fastRPCA")
        fastPCAOut <- fastRPCA::fastPCA(
            A,
            k = K,
            its = q,
            l = (K + 10),
            seed = mkl.seed
        )
        out$u <- fastPCAOut$U
        out$v <- fastPCAOut$V
        out$d <- diag(fastPCAOut$S)
    } else {
        stop('Method not recognized')
    }
    out
}


#' Simple convenience function to library and log normalize a matrix
#'
#'
#' @param A A matrix of gene expression counts (genes x cells)
#' @param verbose Logical indicating whether to print progress messages
#'
#' @export
#'
normalize_data <- function(A, verbose = TRUE) {
    totalUMIPerCell <- SigBridgeRUtils::rowSums3(A)
    if (any(totalUMIPerCell == 0)) {
        toRemove <- which(totalUMIPerCell == 0)
        A <- A[-toRemove, ]
        totalUMIPerCell <- totalUMIPerCell[-toRemove]
        if (verbose) {
            ts_cli$cli_alert_info(sprintf(
                "Removed {.val %d} cells which did not express any genes",
                length(toRemove)
            ))
        }
    }

    # A_norm <- sweep(A, 1, totalUMIPerCell, '/')
    # A_norm <- A_norm * 10E3
    # A_norm <- log(A_norm + 1)
    A_norm <- log((Matrix::Diagonal(x = 1 / totalUMIPerCell) %*% A) * 1e4 + 1L)
}

#' Heuristic for choosing rank k for the low rank approximation based on
#' statistics of the spacings between consecutive singular values. Finds
#' the smallest singular value \eqn{\sigma_i} such that \eqn{\sigma_i - \sigma_{i-1}}
#' is significantly different than spacings in the tail of the singular values.
#'
#' @param A_norm The log-transformed expression matrix of cells (rows) vs. genes (columns)
#' @param K Number of singular values to compute. Must be less than the
#' smallest dimension of the matrix.
#' @param thresh Number of standard deviations away from the ``noise'' singular
#' values which you consider to be signal
#' @param noise_start Index for which all smaller singular values are
#' considered noise
#' @param q Number of additional power iterations
#' @param use.mkl Use the Intel MKL based implementation of SVD. Needs to be
#' installed from https://github.com/KlugerLab/rpca-mkl.
#' @param mkl.seed Only relevant if use.mkl=T. Set the seed for the random
#' generator for the Intel MKL implementation of SVD. Any number <0 will
#' use the current timestamp. If use.mkl=F, set the seed using
#' set.seed() function as usual.
#'
#' @returns
#' A list with three items
#' 1) Chosen k
#' 2) P values of each possible k
#' 3) Singular values of the matrix A_norm
#' @export
choose_k <- function(
    A_norm,
    K = 100,
    thresh = 6,
    noise_start = 80,
    q = 2,
    use.mkl = F,
    mkl.seed = -1
) {
    if (K > min(dim(A_norm))) {
        cli::cli_abort(
            "For an m by n matrix, K must be smaller than the min(m,n)."
        )
    }
    if (noise_start > K - 5) {
        cli::cli_abort(
            "There need to be at least 5 singular values considered noise."
        )
    }
    noise_svals <- noise_start:K

    if (!use.mkl) {
        rsvd_out <- randomized.svd(A = A_norm, K = K, q = q)
    } else {
        rsvd_out <- randomized.svd(
            A = A_norm,
            K = K,
            q = q,
            method = 'rsvd-mkl',
            mkl.seed = mkl.seed
        )
    }

    diffs <- rsvd_out$d[seq_len(length(rsvd_out$d) - 1)] -
        rsvd_out$d[2:length(rsvd_out$d)]

    mu <- mean(diffs[noise_svals - 1])

    sigma <- stats::sd(diffs[noise_svals - 1])

    num_of_sds <- (diffs - mu) / sigma

    k <- max(which(num_of_sds > thresh))

    list(k = k, num_of_sds = num_of_sds, d = rsvd_out$d)
}

#' Computes the k-rank approximation to A_norm and adjusts it according to the
#' error distribution learned from the negative values.
#'
#' @param A_norm The log-transformed expression matrix of cells (rows) vs. genes (columns)
#' @param k the rank of the rank-k approximation. Set to 0 for automated choice of k.
#' @param q the number of additional power iterations in randomized SVD
#' @param quantile.prob quantile probaility, default is 0.001
#' @param use.mkl Use the Intel MKL based implementation of SVD. Needs to be
#' installed from https://github.com/KlugerLab/rpca-mkl.
#' @param mkl.seed Only relevant if use.mkl=T. Set the seed for the random
#' generator for the Intel MKL implementation of SVD. Any number <0 will
#' use the current timestamp. If use.mkl=F, set the seed using
#' set.seed() function as usual.
#' @param ... additional parameters like verbose
#'
#' @returns
#' A list with three items
#' 1) The rank k approximation of A_norm.
#' 2) The rank k approximation of A_norm, adaptively thresholded
#' 3) The rank k approximation of A_norm, adaptively thresholded and
#' with the first two moments of the non-zero values matched to the
#' first two moments of the non-zeros of A_norm. This is the completed
#' matrix most people will want to work with
#'
#' @examples
#' A_norm <- normalize_data(b_nk_example)
#' result.completed <- alra(A_norm,15)
#' # The low rank approximation for reference purposes...not suggested for matrix completion
#' A_norm_rank15 <- result.completed[[1]]
#'
#' # The actual adjusted, completed matrix
#' A_norm_rank15_cor <- result.completed[[3]]
#' @export
alra <- function(
    A_norm,
    k = 0L,
    q = 10L,
    quantile.prob = 0.001,
    use.mkl = FALSE,
    mkl.seed = -1L,
    ...
) {
    dots <- rlang::list2(...)
    verbose <- dots$verbose %||%
        SigBridgeRUtils::getFuncOption("verbose") %||%
        TRUE

    if (verbose) {
        ts_cli$cli_alert_info(sprintf(
            "Read matrix with {.val {nrow(A_norm)}} cells and {.val {ncol(A_norm)}} genes"
        ))
    }

    if (!inherits(A_norm, "matrix") && !inherits(A_norm, "Matrix")) {
        cli::cli_abort(sprintf(
            "{.arg A_norm} is of class {.type %s}, but it should be of class matrix. Did you forget to run as.matrix()?",
            class(A_norm)
        ))
    }

    if (k == 0) {
        k_choice <- choose_k(A_norm)
        k <- k_choice$k
        if (verbose) {
            ts_cli$cli_alert_info("Chose k = {.val {k}}")
        }
    }

    if (verbose) {
        ts_cli$cli_alert_info("Getting nonzeros")
    }
    originally_nonzero <- A_norm > 0

    if (verbose) {
        ts_cli$cli_alert_info("Randomized SVD")
    }

    if (!use.mkl) {
        fastDecomp_noc <- randomized.svd(A = A_norm, K = k, q = q)
    } else {
        fastDecomp_noc <- randomized.svd(
            A = A_norm,
            K = k,
            q = q,
            method = 'rsvd-mkl',
            mkl.seed = mkl.seed
        )
    }

    A_norm_rank_k <- fastDecomp_noc$u[, seq_len(k)] %*%
        diag(fastDecomp_noc$d[seq_len(k)]) %*%
        t(fastDecomp_noc$v[, seq_len(k)])

    if (verbose) {
        ts_cli$cli_alert_info(
            "Find the {.val {quantile.prob}} quantile of each gene"
        )
    }
    # A_norm_rank_k_mins <- abs(apply(A_norm_rank_k,2,min))
    A_norm_rank_k_mins <- abs(
        # convert to vector
        drop(SigBridgeRUtils::colQuantiles3(
            A_norm_rank_k,
            probs = quantile.prob
        ))
    )

    if (verbose) {
        ts_cli$cli_alert_info("Sweeping")
    }

    A_norm_rank_k_cor <- replace(
        A_norm_rank_k,
        A_norm_rank_k <= A_norm_rank_k_mins[col(A_norm_rank_k)],
        0
    )

    sd_nonzero <- function(x) stats::sd(x[!x == 0])
    sigma_1 <- apply(A_norm_rank_k_cor, 2, sd_nonzero)
    sigma_2 <- apply(A_norm, 2, sd_nonzero)

    colSums3 <- SigBridgeRUtils::colSums3

    mu_1 <- colSums3(A_norm_rank_k_cor) / colSums3(!!A_norm_rank_k_cor)
    mu_2 <- colSums3(A_norm) / colSums3(!!A_norm)

    toscale <- !is.na(sigma_1) &
        !is.na(sigma_2) &
        !(sigma_1 == 0 & sigma_2 == 0) &
        !(sigma_1 == 0)

    if (verbose) {
        ts_cli$cli_alert_info(
            "Scaling all except for {.val {sum(!toscale)}} columns"
        )
    }

    sigma_1_2 <- sigma_2 / sigma_1
    toadd <- -1 * mu_1 * sigma_2 / sigma_1 + mu_2

    A_norm_rank_k_temp <- A_norm_rank_k_cor[, toscale]
    A_norm_rank_k_temp <- sweep(
        A_norm_rank_k_temp,
        2,
        sigma_1_2[toscale],
        FUN = "*"
    )
    A_norm_rank_k_temp <- sweep(
        A_norm_rank_k_temp,
        2,
        toadd[toscale],
        FUN = "+"
    )

    A_norm_rank_k_cor_sc <- A_norm_rank_k_cor
    A_norm_rank_k_cor_sc[, toscale] <- A_norm_rank_k_temp
    A_norm_rank_k_cor_sc[A_norm_rank_k_cor == 0] <- 0

    lt0 <- A_norm_rank_k_cor_sc < 0
    A_norm_rank_k_cor_sc[lt0] <- 0

    if (verbose) {
        ts_cli$cli_alert_info(sprintf(
            "{.val {100 * sum(lt0) / (nrow(A_norm) * ncol(A_norm))}}% of the values became negative in the scaling process and were set to zero"
        ))
    }

    # A_norm_rank_k_cor_sc[
    #     originally_nonzero & A_norm_rank_k_cor_sc == 0
    # ] <- A_norm[originally_nonzero & A_norm_rank_k_cor_sc == 0]

    nnz_orig <- which(as.matrix(originally_nonzero), arr.ind = TRUE)

    vals_at_orig_nnz <- A_norm_rank_k_cor_sc[nnz_orig]
    target_idx <- nnz_orig[vals_at_orig_nnz == 0, , drop = FALSE]

    if (nrow(target_idx) > 0) {
        A_norm_rank_k_cor_sc[target_idx] <- A_norm[target_idx]
    }

    col_names <- colnames(A_norm)
    colnames(A_norm_rank_k_cor) <- col_names
    colnames(A_norm_rank_k_cor_sc) <- col_names
    colnames(A_norm_rank_k) <- col_names
    row_names <- rownames(A_norm)
    rownames(A_norm_rank_k_cor) <- row_names
    rownames(A_norm_rank_k_cor_sc) <- row_names
    rownames(A_norm_rank_k) <- row_names

    original_nz <- sum(A_norm > 0) / (nrow(A_norm) * ncol(A_norm))
    completed_nz <- sum(A_norm_rank_k_cor_sc > 0) /
        (nrow(A_norm) * ncol(A_norm))

    if (verbose) {
        ts_cli$cli_alert_info(
            "The matrix went from {.val {round(100 * original_nz, 3)}}% nonzero to {.val {round(100 * completed_nz, 3)}}% nonzero"
        )
    }

    return_Matrix <- inherits(A_norm, "Matrix")
    list(
        A_norm_rank_k = if (return_Matrix) {
            Matrix::Matrix(A_norm_rank_k)
        } else {
            A_norm_rank_k
        },
        A_norm_rank_k_cor = if (return_Matrix) {
            Matrix::Matrix(A_norm_rank_k_cor)
        } else {
            A_norm_rank_k_cor
        },
        A_norm_rank_k_cor_sc = if (return_Matrix) {
            Matrix::Matrix(A_norm_rank_k_cor_sc)
        } else {
            A_norm_rank_k_cor_sc
        }
    )
}
