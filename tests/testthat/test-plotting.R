test_that("Plotting Function - Test plotPP", {
  tmp <- plotPP(TNBC_pheno[
    TNBC_pheno$Person == 1,
    c("cellx", "celly", "Phenotype")
  ])

  expect_s3_class(tmp, "ggplot")
})


test_that("Plotting Function - Test plot_K_function", {
  tmp <- getKFunction(TNBC_pheno[TNBC_pheno$Class == 0, -1],
    agents = c("Tumor", "B"), unit = "Person",
    rCheckVals = seq(0, 50, 1)
  )
  tmp1 <- getKFunction(TNBC_pheno[TNBC_pheno$Class == 1, -1],
    agents = c("Tumor", "B"), unit = "Person",
    rCheckVals = seq(0, 50, 1)
  )

  tmp_1 <- tidyr::pivot_longer(data = tmp, cols = K1:K18)
  tmp1_1 <- tidyr::pivot_longer(data = tmp1, cols = K1:K15)

  data_plot <- rbind(
    data.frame(
      "r" = tmp_1$r,
      "K" = tmp_1$value,
      "unit" = tmp_1$name,
      "outcome" = "0"
    ),
    data.frame(
      "r" = tmp1_1$r,
      "K" = tmp1_1$value,
      "unit" = paste0(tmp1_1$name, "_1"),
      "outcome" = "1"
    )
  )

  tmp <- plot_K_functions(data_plot)

  expect_s3_class(tmp, "ggplot")
})
