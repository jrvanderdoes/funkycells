## All important
dat <- simulatePP(
  cellVarData =
    data.frame("stage" = c(0, 1),
               "A" = c(0, 1/100),
               "B" = c(1 / 50, 1 / 200)),
  cellKappaData = data.frame(
    "cell" = c("A", "B"),
    "clusterCell" = c(NA, "A"),
    "kappa" = c(20, 5)
  ),
  peoplePerStage = 20,
  imagesPerPerson = 1,
  silent = FALSE
)
pcaData <- getPCAData(dat,
                      repeatedUniqueId = "Image",
                      xRange = c(0, 1), yRange = c(0, 1), silent = FALSE
)
pcaMeta <- simulateMeta(pcaData,
                        metaInfo = data.frame(
                          "var" = c("corrNorm"),
                          "rdist" = c("rnorm"),
                          "Stage_0" = c("1"),
                          "Stage_1" = c("2")
                        )
)
saveRDS(list('data'=pcaMeta,'dat'=dat),'C:/Users/jerem/Downloads/allVal.rds')
rfcv <- funkyModel(
  data = pcaMeta, outcome = "Stage", unit = "Person",
  metaNames = c("corrNorm"),
  subsetPlotSize = 25,nTrees = 1000)


## No important
dat <- simulatePP(
  cellVarData =
    data.frame("stage" = c(0, 1),
               "A" = c(0, 0),
               "B" = c(1 / 50, 1 / 50)),
  cellKappaData = data.frame(
    "cell" = c("A", "B"),
    "clusterCell" = c(NA, "A"),
    "kappa" = c(20, 5)
  ),
  peoplePerStage = 20,
  imagesPerPerson = 1,
  silent = FALSE
)
pcaData <- getPCAData(dat,
                      repeatedUniqueId = "Image",
                      xRange = c(0, 1), yRange = c(0, 1), silent = FALSE
)
pcaMeta <- simulateMeta(pcaData,
                        metaInfo = data.frame(
                          "var" = c("corrNorm"),
                          "rdist" = c("rnorm"),
                          "Stage_0" = c("1"),
                          "Stage_1" = c("1")
                        )
)
saveRDS(list('data'=pcaMeta,'dat'=dat),'C:/Users/jerem/Downloads/noval.rds')
rfcv <- funkyModel(
  data = pcaMeta, outcome = "Stage", unit = "Person",
  metaNames = c("corrNorm"),
  subsetPlotSize = 25)
