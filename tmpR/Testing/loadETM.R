
tmp <- readRDS('C:/Users/jerem/Downloads/allVal_extreme.rds')
pcaMeta <- tmp$data
rfcv <- funkyModel(
  data = pcaMeta, outcome = "Stage", unit = "Person",
  metaNames = c("corrNorm"),
  subsetPlotSize = 25)
rfcv1 <- funkyModel1(
  data = pcaMeta, outcome = "Stage", unit = "Person",
  metaNames = c("corrNorm"),
  subsetPlotSize = 25)
rfcv2 <- funkyModel2(
  data = pcaMeta, outcome = "Stage", unit = "Person",
  metaNames = c("corrNorm"),
  subsetPlotSize = 25)
rfcv3 <- funkyModel3(
  data = pcaMeta, outcome = "Stage", unit = "Person",
  metaNames = c("corrNorm"),
  subsetPlotSize = 25)
rfcv4 <- funkyModel4(
  data = pcaMeta, outcome = "Stage", unit = "Person",
  metaNames = c("corrNorm"),
  subsetPlotSize = 25)

rfcv$viPlot
rfcv1$viPlot
rfcv2$viPlot
rfcv3$viPlot
rfcv4$viPlot
