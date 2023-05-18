
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

rfcv$viPlot # Bad
rfcv1$viPlot #Bad
rfcv2$viPlot # Oaky
rfcv3$viPlot # Okay
rfcv4$viPlot # Huge intervals

####################

tmp <- readRDS('C:/Users/jerem/Downloads/noval_extreme.rds')
pcaMeta <- tmp$data
rfcv_no <- funkyModel(
  data = pcaMeta, outcome = "Stage", unit = "Person",
  metaNames = c("corrNorm"),
  subsetPlotSize = 25)
rfcv1_no <- funkyModel1(
  data = pcaMeta, outcome = "Stage", unit = "Person",
  metaNames = c("corrNorm"),
  subsetPlotSize = 25)
rfcv2_no <- funkyModel2(
  data = pcaMeta, outcome = "Stage", unit = "Person",
  metaNames = c("corrNorm"),
  subsetPlotSize = 25)
rfcv3_no <- funkyModel3(
  data = pcaMeta, outcome = "Stage", unit = "Person",
  metaNames = c("corrNorm"),
  subsetPlotSize = 25)
rfcv4_no <- funkyModel4(
  data = pcaMeta, outcome = "Stage", unit = "Person",
  metaNames = c("corrNorm"),
  subsetPlotSize = 25)


rfcv_no$viPlot # bad
rfcv1_no$viPlot # Good
rfcv2_no$viPlot # okay
rfcv3_no$viPlot # okay
rfcv4_no$viPlot # good


####################

tmp <- readRDS('C:/Users/jerem/Downloads/noval.rds')
pcaMeta <- tmp$data
rfcv_no_sm <- funkyModel(
  data = pcaMeta, outcome = "Stage", unit = "Person",
  metaNames = c("corrNorm"),
  subsetPlotSize = 25)
rfcv1_no_sm <- funkyModel1(
  data = pcaMeta, outcome = "Stage", unit = "Person",
  metaNames = c("corrNorm"),
  subsetPlotSize = 25)
rfcv2_no_sm <- funkyModel2(
  data = pcaMeta, outcome = "Stage", unit = "Person",
  metaNames = c("corrNorm"),
  subsetPlotSize = 25)
rfcv3_no_sm <- funkyModel3(
  data = pcaMeta, outcome = "Stage", unit = "Person",
  metaNames = c("corrNorm"),
  subsetPlotSize = 25)
rfcv4_no_sm <- funkyModel4(
  data = pcaMeta, outcome = "Stage", unit = "Person",
  metaNames = c("corrNorm"),
  subsetPlotSize = 25)


rfcv_no_sm$viPlot # Good
rfcv1_no_sm$viPlot #  Good
rfcv2_no_sm$viPlot # Okay
rfcv3_no_sm$viPlot # Good
rfcv4_no_sm$viPlot # Good

####################

tmp <- readRDS('C:/Users/jerem/Downloads/allval.rds')
pcaMeta <- tmp$data
rfcv_all_sm <- funkyModel(
  data = pcaMeta, outcome = "Stage", unit = "Person",
  metaNames = c("corrNorm"),
  subsetPlotSize = 25)
rfcv1_all_sm <- funkyModel1(
  data = pcaMeta, outcome = "Stage", unit = "Person",
  metaNames = c("corrNorm"),
  subsetPlotSize = 25)
rfcv2_all_sm <- funkyModel2(
  data = pcaMeta, outcome = "Stage", unit = "Person",
  metaNames = c("corrNorm"),
  subsetPlotSize = 25)
rfcv3_all_sm <- funkyModel3(
  data = pcaMeta, outcome = "Stage", unit = "Person",
  metaNames = c("corrNorm"),
  subsetPlotSize = 25)
rfcv4_all_sm <- funkyModel4(
  data = pcaMeta, outcome = "Stage", unit = "Person",
  metaNames = c("corrNorm"),
  subsetPlotSize = 25)


rfcv_all_sm$viPlot # okay
rfcv1_all_sm$viPlot # Bad
rfcv2_all_sm$viPlot # Good
rfcv3_all_sm$viPlot # Good
rfcv4_all_sm$viPlot # Good

####################

tmp <- readRDS('C:/Users/jerem/Downloads/noval_big.rds')
pcaMeta <- tmp$data
rfcv_no_bg <- funkyModel(
  data = pcaMeta, outcome = "Stage", unit = "Person",
  metaNames = c("corrNorm"),
  subsetPlotSize = 25)
rfcv1_no_bg <- funkyModel1(
  data = pcaMeta, outcome = "Stage", unit = "Person",
  metaNames = c("corrNorm"),
  subsetPlotSize = 25)
rfcv2_no_bg <- funkyModel2(
  data = pcaMeta, outcome = "Stage", unit = "Person",
  metaNames = c("corrNorm"),
  subsetPlotSize = 25)
rfcv3_no_bg <- funkyModel3(
  data = pcaMeta, outcome = "Stage", unit = "Person",
  metaNames = c("corrNorm"),
  subsetPlotSize = 25)
rfcv4_no_bg <- funkyModel4(
  data = pcaMeta, outcome = "Stage", unit = "Person",
  metaNames = c("corrNorm"),
  subsetPlotSize = 25)


rfcv_no_bg$viPlot # Bad
rfcv1_no_bg$viPlot #  Good
rfcv2_no_bg$viPlot # Bad
rfcv3_no_bg$viPlot # Bad
rfcv4_no_bg$viPlot # Bad

####################

tmp <- readRDS('C:/Users/jerem/Downloads/allval_big.rds')
pcaMeta <- tmp$data
rfcv_all_bg <- funkyModel(
  data = pcaMeta, outcome = "Stage", unit = "Person",
  metaNames = c("corrNorm"),
  subsetPlotSize = 25)
rfcv1_all_bg <- funkyModel1(
  data = pcaMeta, outcome = "Stage", unit = "Person",
  metaNames = c("corrNorm"),
  subsetPlotSize = 25)
rfcv2_all_bg <- funkyModel2(
  data = pcaMeta, outcome = "Stage", unit = "Person",
  metaNames = c("corrNorm"),
  subsetPlotSize = 25)
rfcv3_all_bg <- funkyModel3(
  data = pcaMeta, outcome = "Stage", unit = "Person",
  metaNames = c("corrNorm"),
  subsetPlotSize = 25)
rfcv4_all_bg <- funkyModel4(
  data = pcaMeta, outcome = "Stage", unit = "Person",
  metaNames = c("corrNorm"),
  subsetPlotSize = 25)


rfcv_all_bg$viPlot # Bad
rfcv1_all_bg$viPlot #  Bad
rfcv2_all_bg$viPlot # Okay
rfcv3_all_bg$viPlot # Okay
rfcv4_all_bg$viPlot # Okay


##############################

rfcv$viPlot # Bad
rfcv1$viPlot #Bad
rfcv2$viPlot # Oaky
rfcv3$viPlot # Okay
rfcv4$viPlot # Huge intervals

rfcv_no$viPlot # bad
rfcv1_no$viPlot # Good
rfcv2_no$viPlot # okay
rfcv3_no$viPlot # okay
rfcv4_no$viPlot # good

rfcv_no_sm$viPlot # Good
rfcv1_no_sm$viPlot #  Good
rfcv2_no_sm$viPlot # Okay
rfcv3_no_sm$viPlot # Good
rfcv4_no_sm$viPlot # Good

rfcv_all_sm$viPlot # okay
rfcv1_all_sm$viPlot # Bad
rfcv2_all_sm$viPlot # Good
rfcv3_all_sm$viPlot # Good
rfcv4_all_sm$viPlot # Good

rfcv_no_bg$viPlot # Bad
rfcv1_no_bg$viPlot #  Good
rfcv2_no_bg$viPlot # Bad
rfcv3_no_bg$viPlot # Bad
rfcv4_no_bg$viPlot # Bad

rfcv_all_bg$viPlot # Bad
rfcv1_all_bg$viPlot #  Bad
rfcv2_all_bg$viPlot # Okay
rfcv3_all_bg$viPlot # Okay
rfcv4_all_bg$viPlot # Okay

