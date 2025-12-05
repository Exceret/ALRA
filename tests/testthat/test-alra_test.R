data("b_nk_example")
data("labels_example")

A <- b_nk_example
labels <- labels_example

test_that("alra works", {
    # Library and log normalize the data
    A_norm <- normalize_data(A)

    # Choose k.
    k_choice <- choose_k(A_norm)

    # run ALPRA
    A_norm_completed <- alra(A_norm, k = k_choice$k)[[3]]
    # ℹ [2025/12/05 20:33:56] Read matrix with 16427 cells and 12776 genes
    # ℹ [2025/12/05 20:33:56] Getting nonzeros
    # ℹ [2025/12/05 20:33:56] Randomized SVD
    # ℹ [2025/12/05 20:33:58] Find the 0.001 quantile of each gene
    # ℹ [2025/12/05 20:34:03] Sweeping
    # ℹ [2025/12/05 20:34:26] Scaling all except for 351 columns
    # ℹ [2025/12/05 20:34:37] 0% of the values became negative in the scaling process and were set to zero
    # ℹ [2025/12/05 20:34:41] The matrix went from 5.107% nonzero to 73.457% nonzero

    # type check
    expect_equal(class(k_choice), "list")
    expect_equal(class(k_choice$k), "integer")
    expect_equal(class(k_choice$num_of_sds), "numeric")
    expect_equal(class(k_choice$d), "numeric")
    expect_s4_class(A_norm, "dgeMatrix")
    expect_s4_class(A_norm_completed, "dgeMatrix")

    # quantity check
    expect_true(k_choice$k > 1)
    expect_true(length(k_choice$num_of_sds) > 0)
    expect_true(length(k_choice$d) > 0)

    # NCAM1 should be expressed in all NK cells, but is only expressed in 4% of of NK cells in the original data.
    expect_lte(
        Matrix::mean(A_norm[labels == "cd56_nk", "NCAM1", drop = FALSE] > 0),
        .1
    )
    expect_gte(
        Matrix::mean(
            A_norm_completed[labels == "cd56_nk", "NCAM1", drop = FALSE] > 0
        ),
        .5
    )

    # CR2 should be expressed in all B cells, but is only expressed in 1% of of B cells in the original data.
    expect_lte(
        Matrix::mean(A_norm[labels == "b_cells", "CR2", drop = FALSE] > 0),
        .1
    )
    expect_gte(
        Matrix::mean(
            A_norm_completed[labels == "b_cells", "CR2", drop = FALSE] > 0
        ),
        .5
    )
})

test_that("alra works with S4 Matrix", {
    A2 <- Matrix::Matrix(A)

    # Library and log normalize the data
    A_norm2 <- normalize_data(A2)

    # Choose k.
    k_choice2 <- choose_k(A_norm2)

    # run ALPRA
    A_norm_completed2 <- alra(A_norm2, k = k_choice2$k)[[3]]
    # ℹ [2025/12/04 19:54:09] Read matrix with 16427 cells and 12776 genes
    # ℹ [2025/12/04 19:54:09] Getting nonzeros
    # ℹ [2025/12/04 19:54:11] Randomized SVD
    # ℹ [2025/12/04 19:54:12] Find the 0.001 quantile of each gene
    # ℹ [2025/12/04 19:54:16] Sweeping
    # ℹ [2025/12/04 19:54:36] Scaling all except for 351 columns
    # ℹ [2025/12/04 19:54:48] 0.00% of the values became negative in the scaling process and were set to zero
    # ℹ [2025/12/04 19:54:52] The matrix went from 5.107% nonzero to 75.206% nonzero

    # NCAM1 should be expressed in all NK cells, but is only expressed in 4% of of NK cells in the original data.
    # NCAM1 should be expressed in all NK cells, but is only expressed in 4% of of NK cells in the original data.
    expect_lte(
        Matrix::mean(A_norm[labels == "cd56_nk", "NCAM1", drop = FALSE] > 0),
        .1
    )
    expect_gte(
        Matrix::mean(
            A_norm_completed[labels == "cd56_nk", "NCAM1", drop = FALSE] > 0
        ),
        .5
    )

    # CR2 should be expressed in all B cells, but is only expressed in 1% of of B cells in the original data.
    expect_lte(
        Matrix::mean(A_norm[labels == "b_cells", "CR2", drop = FALSE] > 0),
        .1
    )
    expect_gte(
        Matrix::mean(
            A_norm_completed[labels == "b_cells", "CR2", drop = FALSE] > 0
        ),
        .5
    )
})
