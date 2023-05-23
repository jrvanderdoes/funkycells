## All important
dat <- simulatePP(
  cellVarData =
    data.frame("stage" = c(0, 1),
               "A" = c(0, 1/100),
               "B" = c(1 / 50, 1 / 200),
               "C" = c(1 / 100, 1 / 20),
               "D" = c(1 / 5, 1 / 100),
               "E" = c(1 / 20, 1 / 70),
               "F" = c(1 / 50, 1 / 150)),
  cellKappaData = data.frame(
    "cell" = c("A", "B","C", "D", "E", "F"),
    "clusterCell" = c(NA, "A", "A", "B", NA, "E"),
    "kappa" = c(20, 5, 3, 2, 20, 5)
  ),
  peoplePerStage = 20,
  imagesPerPerson = 1,
  silent = FALSE
)
interactions <- data.frame(c('A',"B","C","D","E","F",
                             'A','A','D','F'),
                           c('A',"B","C","D","E","F",
                             'B','C','B','E'))

# cells <- unique(dat$cellType)
# interactions <- rbind(data.frame(t(combn(cells,2))),
#                       data.frame('X1'=cells,'X2'=cells))


pcaData <- getPCAData(dat,agents_df = interactions,
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
saveRDS(list('data'=pcaMeta,'dat'=dat),'C:/Users/jerem/Downloads/allVal_big.rds')
rfcv <- funkyModel(
  data = pcaMeta, outcome = "Stage", unit = "Person",
  metaNames = c("corrNorm"),
  subsetPlotSize = 25)


## No important
dat <- simulatePP(
  cellVarData =
    data.frame("stage" = c(0, 1),
               "A" = c(1/100, 1/100),
               "B" = c(1 / 200, 1 / 200),
               "C" = c(1 / 20, 1 / 20),
               "D" = c(1 / 100, 1 / 100),
               "E" = c(1 / 70, 1 / 70),
               "F" = c(1 / 150, 1 / 150)),
  cellKappaData = data.frame(
    "cell" = c("A", "B","C", "D", "E", "F"),
    "clusterCell" = c(NA, "A", "A", "B", NA, "E"),
    "kappa" = c(20, 5, 3, 2, 20, 5)
  ),
  peoplePerStage = 20,
  imagesPerPerson = 1,
  silent = FALSE
)
interactions <- data.frame(c('A',"B","C","D","E","F",
                             'A','A','D','F'),
                           c('A',"B","C","D","E","F",
                             'B','C','B','E'))

# cells <- unique(dat$cellType)
# interactions <- rbind(data.frame(t(combn(cells,2))),
#                       data.frame('X1'=cells,'X2'=cells))


pcaData <- getPCAData(dat,agents_df = interactions,
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
saveRDS(list('data'=pcaMeta,'dat'=dat),'C:/Users/jerem/Downloads/noval_big.rds')
rfcv <- funkyModel(
  data = pcaMeta, outcome = "Stage", unit = "Person",
  metaNames = c("corrNorm"),
  subsetPlotSize = 25)
