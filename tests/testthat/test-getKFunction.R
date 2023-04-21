test_that("K Function - Check existant K", {
  KFunction <- getKFunction(
    agents = c("B", "Tumour"), unit = "Person",
    data = TNBC_pheno[TNBC_pheno$Person == 1, -1],
    rCheckVals = seq(0, 50, 1),
    edgeCorrection = "isotropic"
  )

  # Check Radi
  expect_equal(0:50, KFunction$r)
  # Check K
  expect_equal(c(
    0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
    0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
    1.381006, 4.143018, 12.573413, 18.097437, 44.336550, 89.242898,
    119.552954, 166.507156, 223.304416, 289.194499, 361.562910,
    455.437013, 526.922843, 601.497163, 701.528950, 767.466969,
    863.567377, 939.746958, 1037.002729, 1119.274950, 1209.951945,
    1304.491041, 1387.524038, 1483.379199, 1563.344451,
    1660.082556, 1768.439650, 1867.575888, 1981.296872,
    2085.793752, 2193.149356, 2316.333468, 2422.376019,
    2536.355854, 2666.937927, 2790.817201, 2938.092315,
    3058.739663, 3199.854683
  ), KFunction$K1)
})


test_that("K Function - Check non-existant K", {
  KFunction <- getKFunction(
    agents = c("B", "FAKE"), unit = "Person",
    data = TNBC_pheno[TNBC_pheno$Person == 1, -1],
    rCheckVals = seq(0, 50, 1),
    edgeCorrection = "isotropic"
  )

  # Only one PCA
  expect_equal(0, nrow(KFunction))
})
