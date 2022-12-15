#' Title
#'
#' @param dat
#' @param nSims
#' @param outcome
#' @param unit
#' @param repeatedId
#' @param noiseMap
#' @param KFunctions
#' @param metaNames
#' @param syntheticMetaNames
#' @param alpha
#' @param silent
#' @param alignmentMethod
#'
#' @return
#' @export
#'
#' @examples
.generateCurveNoise <- function(dat, nSims,
                    outcome, unit, repeatedId, noiseMap,
                    KFunctions, metaNames, syntheticMetaNames,
                    alpha=0.05, silent=F, alignmentMethod='Add'){
  if(nSims<=0) return(NULL)
  if(nSims==1){
    warning('Warning: curvedSigSims must be greater than 1 to work. Setting to 0')
    return(NULL)
  }
  # Organize Variables
  underlyingDataAlignedFunctions <- .getUnderlyingVariable(colnames(dat),
                                                           returnUnique = F)
  underlyingVars <- unique(underlyingDataAlignedFunctions)
  underlyingVars <- underlyingVars[!(underlyingVars %in%
                                       c(outcome,unit,repeatedId))]

  # Setup Variables
  #simData <- RF <- list()
  simData <- list()
  avgVI <- avgGini <- data.frame('var'=underlyingVars)

  # Bootstrap data and simulate
  for(sim in 1:nSims){
    if(!silent) cat(sim,'/',nSims,'\n')
    tmpDF <- dat[,colnames(dat) %in% c(outcome,unit,repeatedId)]

    for(varIdx in 1:length(underlyingVars)){
      tmpDF <- cbind(tmpDF,
                     dat[sample(1:nrow(dat)),
                         underlyingDataAlignedFunctions==underlyingVars[varIdx],
                         drop=F])
    }
    # simData[[sim]] <- tmpDF

    # Get RF and VI
    #RF[[sim]] <- .computeRandomForest_PC(data=tmpDF,#simData[[sim]],
    RF <- .computeRandomForest_PC(data=tmpDF,#simData[[sim]],
                                  outcome = outcome,
                                  unit=unit, repeatedId=repeatedId,
                                  varImpPlot = F,
                                  metaNames=c(metaNames,syntheticMetaNames))

    #data_merge <- RF[[sim]][[2]][,c('var','avgGini')]
    data_merge <- RF[[2]][,c('var','avgGini')]
    colnames(data_merge) <- c('var',paste0('avgGiniK',sim))
    avgGini <- merge(avgGini, data_merge, by='var')

    #data_merge <- RF[[sim]][[2]][,c('var','avgVI')]
    data_merge <- RF[[2]][,c('var','avgVI')]
    colnames(data_merge) <- c('var',paste0('avgVIK',sim))
    avgVI <- merge(avgVI, data_merge, by='var')
  }

  ## Standardize results
  if(length(alignmentMethod)==2){
    # Return both methods
    stdCurveList_add <- .getStdCurveData(avgGini, avgVI, alpha,
                                         KFunctions, metaNames, noiseMap,
                                         'Add')

    stdCurveList_mult <- .getStdCurveData(avgGini, avgVI, alpha,
                                          KFunctions, metaNames, noiseMap,
                                          'Mult')

    ## Get Curve Values
    dataReturn <-
      data.frame(
        'gini_add'=.getCurveValues(stdCurveList_add$gini[-1],alpha),
        'vi_add'=.getCurveValues(stdCurveList_add$vi[-1],alpha),
        'gini_mult'=.getCurveValues(stdCurveList_mult$gini[-1],alpha),
        'vi_mult'=.getCurveValues(stdCurveList_mult$vi[-1],alpha) )
  } else{
    # Return the desired method
    stdCurveList <- .getStdCurveData(avgGini, avgVI, alpha,
                                     KFunctions, metaNames, noiseMap,
                                     alignmentMethod)

    ## Get Curve Values
    dataReturn <-
      data.frame(
        'gini'=.getCurveValues(stdCurveList$gini[-1],alpha),
        'vi'=.getCurveValues(stdCurveList$vi[-1],alpha) )
  }

  dataReturn
}



#' Title
#'
#' @param avgGini
#' @param avgVI
#' @param alpha
#' @param KFunctions
#' @param metaNames
#' @param noiseMap
#' @param alignmentMethod
#'
#' @return
#' @export
#'
#' @examples
.getStdCurveData <- function(avgGini, avgVI, alpha,
                             KFunctions, metaNames, noiseMap,
                             alignmentMethod){

  # VarImpList only for noise organization
  varImpList <- .getVariableImportanceMetrics(avgGini, avgVI, NULL, NULL, alpha)

  cutoffs<- .getIndividualCutoffs(noiseMap = noiseMap,
                                       giniData = varImpList$giniData,
                                       viData = varImpList$viData,
                                       alignmentMethod = alignmentMethod,
                                  alpha=alpha)

  # Return Std Curve Data
  list('gini'=.stdCurveData(avgGini, cutoffs$giniCutoff_df,
                                KFunctions, metaNames,alignmentMethod),
       'vi'=.stdCurveData(avgVI, cutoffs$viCutoff_df,
                              KFunctions, metaNames,alignmentMethod))
}


#' Title
#'
#' @param datAvg
#' @param cutoff_df
#' @param KFunctions
#' @param metaNames
#' @param alignmentMethod
#'
#' @return
#' @export
#'
#' @examples
.stdCurveData <- function(datAvg, cutoff_df, KFunctions, metaNames,
                          alignmentMethod){

  tmp <- apply(X=datAvg[datAvg$var%in%c(KFunctions, metaNames),],1,
               function(rowData,cutoff,KFunctions,metaNames,alignmentMethod){
                 if(rowData[1] %in% cutoff$dataVar) {
                   adjVal <- cutoff[cutoff$dataVar==rowData[1],'adj']
                 }else if(rowData[1] %in% KFunctions){
                   adjVal <- cutoff[cutoff$dataVar=='syntheticK','adj']
                 }else if(trowData[1] %in% metaNames){
                   adjVal <- cutoff[cutoff$dataVar=='syntheticMeta','adj']
                 }else{
                   warning(paste0(rowData[1],' was not standardized.\n'))
                   if(alignmentMethod=='Add'){
                     adjVal <- 0
                   }else if(alignmentMethod=='Mult'){
                     adjVal <- 1
                   }
                 }

                 if(alignmentMethod=='Add'){
                   retDat <- c(rowData[1], as.numeric(rowData[-1])+adjVal)
                 }else if(alignmentMethod=='Mult'){
                   retDat <- c(rowData[1], as.numeric(rowData[-1])*adjVal)
                 }

                 retDat
               }, cutoff=cutoff_df,
               KFunctions=KFunctions,metaNames=metaNames,
               alignmentMethod=alignmentMethod)

  stdData <- as.data.frame(t(tmp))
  colnames(stdData) <- colnames(datAvg)
  for(i in 2:ncol(stdData)){
    stdData[,i] <- as.numeric(stdData[,i])
  }

  stdData
}


#' Title
#'
#' @param stdData
#' @param alpha
#'
#' @return
#' @export
#'
#' @examples
.getCurveValues <- function(stdData,alpha){
  apply(apply(stdData[-1], MARGIN=2, FUN=sort, decreasing = T),
        MARGIN = 1, FUN = quantile, probs = 1-alpha)

}
