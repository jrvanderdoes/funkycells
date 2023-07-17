test_that("PCA Function - All Classes, people, pca present", {
  dataPCA_pheno <- getKsPCAData(
    data = TNBC_pheno, unit = "Person",
    agents_df = data.frame("B", "Tumor"),
    nPCs = 3,
    rCheckVals = seq(0, 50, 1), silent = TRUE
  )

  # All classes present
  expect_setequal(c(0, 1), unique(dataPCA_pheno$Class))
  # All People computed
  expect_equal(c(1:14, 16:18, 20:21, 23, 27:29, 31:37, 39:41),
               dataPCA_pheno$Person)
  # All pca's considered
  expect_equal(5, ncol(dataPCA_pheno))
})


test_that("PCA Function - Give NA when non-existant Interaction", {
  dataPCA_pheno <- getKsPCAData(
    data = TNBC_pheno, unit = "Person",
    agents_df = data.frame("B", "FAKE"),
    nPCs = 3,
    rCheckVals = seq(0, 50, 1), silent = TRUE
  )

  # Only one PCA
  expect_equal(3, ncol(dataPCA_pheno))
  # All NA in PC
  expect_equal(nrow(dataPCA_pheno), sum(is.na(dataPCA_pheno$B_FAKE_PC)))
})


test_that("PCA Function - Check values", {
  dataPCA_pheno1 <- getKsPCAData(
    data = TNBC_pheno, unit = "Person",
    agents_df = data.frame("B", "B"),
    nPCs = 1,
    rCheckVals = seq(0, 50, 1), silent = TRUE
  )
  dataPCA_pheno3 <- getKsPCAData(
    data = TNBC_pheno, unit = "Person",
    agents_df = data.frame("B", "B"),
    nPCs = 3,
    rCheckVals = seq(0, 50, 1), silent = TRUE
  )

  # Only one PCA
  expect_equal(
    c(
      -239796.36, -384747.62, -279741.44, -305658.81, -384747.62,
      187159.99, NA, -384747.62, -83035.85, -292083.60, NA,
      -171701.38, -314443.39, NA, 717477.43, -243126.94,
      -384747.62, -73120.55, NA, -384747.62, NA, -310412.18,
      NA, NA, 58436.28, -384747.62, 3822328.01, -314793.06,
      788627.16, -278852.58, 91650.73, -243926.94, -206500.81
    ),
    dataPCA_pheno1$B_B_PC1
  )
  expect_equal(
    c(
      -239319.79, -384168.52, -279233.51, -305145.18, -384168.52,
      187116.97, NA, -384168.52, -82864.59, -291194.25, NA,
      -171323.79, -313902.44, NA, 717236.33, -242664.79, -384168.52,
      -72793.58, NA, -384168.52, NA, -309887.79, NA, NA, 58775.92,
      -384168.52, 3817344.68, -314264.52, 783198.30, -278377.71,
      91828.64, -243453.86, -206063.92
    ),
    dataPCA_pheno3$B_B_PC1
  )
})
