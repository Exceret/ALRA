A <- b_nk_example
labels <- labels_example

test_that("alra works", {
    # Library and log normalize the data
    A_norm <- normalize_data(A)

    # Choose k.
    k_choice <- choose_k(A_norm)

    # run ALPRA
    A_norm_completed <- alra(A_norm, k = k_choice$k)[[3]]

    # NCAM1 should be expressed in all NK cells, but is only expressed in 4% of of NK cells in the original data.
    expect_lte(mean(A_norm[labels == "cd56_nk", "NCAM1", drop = FALSE] > 0), .1)
    expect_gte(
        mean(A_norm_completed[labels == "cd56_nk", "NCAM1", drop = FALSE] > 0),
        .5
    )

    # CR2 should be expressed in all B cells, but is only expressed in 1% of of B cells in the original data.
    expect_lte(mean(A_norm[labels == "b_cells", "CR2", drop = FALSE] > 0), .1)
    expect_gte(
        mean(A_norm_completed[labels == "b_cells", "CR2", drop = FALSE] > 0),
        .5
    )
})

test_that("alra works with S4 Matrix", {
    A <- Matrix::Matrix(A)

    # Library and log normalize the data
    A_norm <- normalize_data(A)

    # Choose k.
    k_choice <- choose_k(A_norm)

    # run ALPRA
    A_norm_completed <- alra(A_norm, k = k_choice$k)[[3]]
    # ℹ [2025/12/04 19:54:09] Read matrix with 16427 cells and 12776 genes
    # ℹ [2025/12/04 19:54:09] Getting nonzeros
    # ℹ [2025/12/04 19:54:11] Randomized SVD
    # ℹ [2025/12/04 19:54:12] Find the 0.001 quantile of each gene
    # ℹ [2025/12/04 19:54:16] Sweeping
    # ℹ [2025/12/04 19:54:36] Scaling all except for 351 columns
    # ℹ [2025/12/04 19:54:48] 0.00% of the values became negative in the scaling process and were set to zero
    # ℹ [2025/12/04 19:54:52] The matrix went from 5.107% nonzero to 75.206% nonzero

    # NCAM1 should be expressed in all NK cells, but is only expressed in 4% of of NK cells in the original data.
    expect_lte(mean(A_norm[labels == "cd56_nk", "NCAM1"] > 0), .1)
    expect_gte(mean(A_norm_completed[labels == "cd56_nk", "NCAM1"] > 0), .5)

    # CR2 should be expressed in all B cells, but is only expressed in 1% of of B cells in the original data.
    expect_lte(mean(A_norm[labels == "b_cells", "CR2"]), .1)
    expect_gte(mean(A_norm_completed[labels == "b_cells", "CR2"] > 0), .5)
})
