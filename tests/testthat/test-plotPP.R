test_that("plotPP Function - Test plot", {
  tmp <- plotPP(TNBC_pheno[TNBC_pheno$Person == 1, c("cellx", "celly", "Phenotype")])

  expect_s3_class(tmp, "ggplot")
})
