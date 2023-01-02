#' Generate Noise Curve
#'
#' This (internal) function generates the noise course of data. It indicates
#'     the variable importance expected for a group of comparable noise
#'     variables. This can combat mistakenly thinking the variable with the
#'     greatest variable importance is significant.
#'
#' @param dat Data.frame used to fit random forest. See data given to
#'     computeRandomForest_CVPC, but note synthetics will also be attached.
#' @param nSims Numeric indicating the number of simulations to bootstrap the
#'     data .
#' @param outcome String indicating the column name with the outcome in dat.
#' @param unit String indicating the column name with the unit in dat.
#' @param repeatedId String indicating the column name with the unique id from
#'     repeated measures in dat.
#' @param noiseMap Data.frame for matching noise variables to noise groups. The
#'     columns are noiseVar (the synthetic interactions of which there are
#'     syntheticKs), dataVar (name of noise group), and synType ('K', which will
#'     compare to 'Meta' from other functions).
#' @param KFunctions Vector indicating the K functions to be considered. These
#'     string names should be in the column names of dat.
#' @param metaNames Vector indicating the meta-variables to be considered. These
#'     string names should be in the column names of dat.
#' @param syntheticMetaNames Vector indicating the synthetic meta-variables to
#'     be considered. These string names should be in the column names of dat.
#' @param alpha (Optional) Numeric in (0,1) indicating the significance used
#'     throughout the analysis. Default is 0.05.
#' @param silent (Optional) Boolean indicating if output should be suppressed
#'     when the function is running. Default is FALSE.
#' @param alignmentMethod (Optional) String indicating the method of aligning
#'     the variables. The options are 'Add', 'Mult', or c('Add','Mult'). The
#'     default is 'Add'.
#'
#' @return Data.frame with 2 - 4 columns. If only one alignmentMethod is given,
#'     a 2 column data.frame is returned with the first column for gini and the
#'     second column for variable importance (both standardized). If two
#'     alignmentMethod values are given, the first two columns related to gini
#'     and variable importance for the additional standardization while columns
#'     3-4  related to the multiplication standardization.
#' @export
#'
#' @examples
#' # See code for computeRandomForest_CVPC. This is not an outward
#' #     function so won't be viewable.
.generateNoiseCurve <- function(dat, nSims,
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
  simData <- list()
  avgVI <- avgGini <- data.frame('var'=underlyingVars)

  # Bootstrap data and simulate
  if(!silent) cat('Curved Sims (',nSims,'): ',sep='')
  for(sim in 1:nSims){
    if(!silent) cat(sim,', ',sep='')
    tmpDF <- dat[,colnames(dat) %in% c(outcome,unit,repeatedId)]

    for(varIdx in 1:length(underlyingVars)){
      tmpDF <- cbind(tmpDF,
                     dat[sample(1:nrow(dat)),
                         underlyingDataAlignedFunctions==underlyingVars[varIdx],
                         drop=F])
    }
    # simData[[sim]] <- tmpDF

    # Get RF and VI
    RF <- .computeRandomForest_PC(data=tmpDF,#simData[[sim]],
                                  outcome = outcome,
                                  unit=unit, repeatedId=repeatedId,
                                  varImpPlot = F,
                                  metaNames=c(metaNames,syntheticMetaNames))

    data_merge <- RF[[2]][,c('var','avgGini')]
    colnames(data_merge) <- c('var',paste0('avgGiniK',sim))
    avgGini <- merge(avgGini, data_merge, by='var')

    data_merge <- RF[[2]][,c('var','avgVI')]
    colnames(data_merge) <- c('var',paste0('avgVIK',sim))
    avgVI <- merge(avgVI, data_merge, by='var')
  }
  if(!silent) cat('\n')

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



#' Get standardized Curve Data
#'
#' This (internal) function gets variable importance metrics and standardizes
#'     the data used in the curved cutoff.
#'
#' @param avgGini Data.frame with gini index values from a random forest in
#'     columns, each column representing a trial and a row for each variable,
#'     with the first column being a variable name column.
#' @param avgVI Data.frame with variable importance values from a random forest
#'     in columns, each column representing a trial and a row for each variable,
#'     with the first column being a variable name column.
#' @param alpha Numeric in (0,1) indicating the significance used
#'     throughout the analysis. Default is 0.05.
#' @param KFunctions Vector indicating the K functions of interest.
#' @param metaNames Vector indicating the meta-variables of interest.
#' @param noiseMap Data.frame for matching noise variables to noise groups.
#'                 The columns are noiseVar (the synthetic interactions of which
#'                 there are syntheticKs), dataVar (name of noise group), and
#'                 synType ('K', which will compare to 'Meta' from other
#'                 functions).
#' @param alignmentMethod String indicating the method of aligning the
#'     variables. The options are 'Add' or 'Mult'.
#'
#' @return List containing two data.frame elements of the same structure. The
#'     first column is var and indicates the variable and subsequent columns are
#'     the standardized values for each sim. The list elements are:
#'     1. gini: Data.frame built using the gini index.
#'     2. vi: Data.frame built using the variable importance metric.
#' @export
#'
#' @examples
#' # See code for .generateNoiseCurve. This is not an outward
#' #     function so won't be viewable.
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


#' Standardize Curve Data
#'
#' This (internal) function standardizes the data in datAvg using cutoffs in
#'     cutoff_df. It returns the standardization for variables of interest,
#'     indicated in KFunctions and metaNames. Standardization is completed
#'     according to alignmentMethod.
#'
#' @param datAvg Data.frame with first column, var, indicating variable (inc.
#'     noise variables) and subsequent columns indicating the values for
#'     importance of each variable.
#' @param cutoff_df Data.frame with 4 columns. The column dataVar indicates the
#'     variables which are related, synType indicates if the variable is K or
#'     meta-variable, adj indicates how much the related variables will need to
#'     be adjusted, and cutoff indicates the current cutoffs for the variables.
#' @param KFunctions Vector indicating the K functions of interest.
#' @param metaNames Vector indicating the meta-variables of interest.
#' @param alignmentMethod String indicating the method of aligning the
#'     variables. The options are 'Add' or 'Mult'.
#'
#' @return Data.frame with the first column (var) indicating variables of
#'     interest and the rest being the standardized values from datAvg.
#' @export
#'
#' @examples
#' # See code for .getStdCurveData. This is not an outward
#' #     function so won't be viewable.
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


#' Get Curve Values
#'
#' The (internal) function gets the 1-alpha quantile of each relative location
#'     in data. This accounts for the fact that variables may change relative
#'     locations on different runs.
#'
#' @param stdData Data.frame with the variables as rows and the point estimates
#'     for each run, given in each column.
#' @param alpha Numeric in (0,1) for finding the 1-alpha quantile.
#'
#' @return This function returns a vector of numerics indicating the expected
#'     noise values for an equivalent number of random variables.
#' @export
#'
#' @examples
.getCurveValues <- function(stdData,alpha){
  apply(apply(stdData[-1], MARGIN=2, FUN=sort, decreasing = T),
        MARGIN = 1, FUN = quantile, probs = 1-alpha)

}
