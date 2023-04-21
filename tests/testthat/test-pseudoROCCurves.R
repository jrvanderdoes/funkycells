test_that("PseudoROC - Ensure figure produced", {
  skip_if_not_installed("pROC")

  data_pp_roc <- simulatePP(
    cellVarData =
      data.frame(
        "stage" = c(0, 1),
        "A" = c(0, 0),
        "B" = c(1 / 50, 1 / 50)
      ),
    cellKappaData = data.frame(
      "cell" = c("A", "B"),
      "clusterCell" = c(NA, "A"),
      "kappa" = c(20, 5)
    ),
    peoplePerStage = 25,
    imagesPerPerson = 1,
    silent = TRUE
  )
  pcaData_roc <- getPCAData(data_pp_roc,
    repeatedUniqueId = "Image",
    xRange = c(0, 1), yRange = c(0, 1), silent = TRUE
  )
  RF_roc <- funkyForest(data = pcaData_roc[-2], nTrees = 5)
  pred_roc <- predict_funkyForest(
    model = RF_roc$model,
    data_pred = pcaData_roc[-2],
    data = pcaData_roc[-2]
  )
  tmp <- computePseudoROCCurves(
    trueOutcomes = pcaData_roc$Stage,
    modelPercents = pred_roc$PredPerc[-1]
  )

  expect_s3_class(tmp, "ggplot")
})
