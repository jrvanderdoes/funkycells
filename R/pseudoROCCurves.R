#' Compute Pseudo-ROC Curves
#'
#' True ROC can only be computed for two-options, but we can consider each one
#'     as right or wrong and overlay the plots. Note this means the lines may
#'     cover each other
#'
#' Upcoming: Suggests for pROC and RCOlorBrewer then check. See the following:
#'     - https://stackoverflow.com/questions/6895852/load-a-package-only-when-needed-in-r-packagesee
#'     - https://r-pkgs.org/metadata.html#sec-description
#'
#' @param trueOutcomes Vector of the true results
#' @param modelPercents Data.frame with columns named after the true outcomes,
#'     giving the percent of selecting that outcome. This is what is returned
#'     .predict.RandomForest_PC when type='all' and look at PredPerc[-1] (first
#'     column is the predictions).
#'
#' @return NULL. See plot of ROC curves.
#' @export
#'
#' @examples
#' data_pp_roc <- simulatePP(cellVarData=
#'                             data.frame('stage'=c(0,1),
#'                                        'A'=c(0,0),
#'                                        'B'=c(1/50,1/50)),
#'                           cellKappaData=data.frame(
#'                                        'cell'=c('A','B'),
#'                                        'clusterCell'=c(NA,'A'),
#'                                        'kappa'=c(20,5)),
#'                           peoplePerStage=100,
#'                           imagesPerPerson=1,
#'                           reduceEdge=0.025,
#'                           silent=F )
#' pcaData_roc <- getPCAData(data_pp_roc,repeatedUniqueId='Image',
#'                           xRange = c(0,1),  yRange = c(0,1), silent=F)
#' RF_roc <- .computeRandomForest_PC(data=pcaData_roc[-2], nTrees = 5)
#' pred_roc <- .predict.RandomForest_PC(model = RF_roc[[1]],
#'                                      data_pred = pcaData_roc[-2],
#'                                      data=pcaData_roc[-2])
#' computePseudoROCCurves(pcaData_roc$Stage,pred_roc$PredPerc[-1])
computePseudoROCCurves <- function(trueOutcomes, modelPercents){
  options <- colnames(modelPercents)
  cols <- suppressWarnings(RColorBrewer::brewer.pal(length(options),'Set1'))
  # 1 ROC curve, mock vs non mock
  roc.curve <- pROC::roc(ifelse(trueOutcomes==options[1],
                                options[1], paste0('Not-',options[1])),
                         as.numeric(modelPercents[,options[1]]),
                         levels=c(options[1], paste0('Not-',options[1])),
                         direction='>'
  )
  plot(roc.curve, col = cols[1])

  for(i in 2:length(options)){
    roc.curve <- pROC::roc(ifelse(trueOutcomes==options[i],
                                  options[i], paste0('Not-',options[i])),
                           as.numeric(modelPercents[,options[i]]),
                           levels=c(options[i], paste0('Not-',options[i])),
                           direction='>')
    lines(roc.curve, col = cols[i])
  }
}
