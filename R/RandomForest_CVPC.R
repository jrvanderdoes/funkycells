#' Fit a Random Forest model with PC data (Using CV for Improvements)
#'
#' The function fits a random forest model to the data along with using cross-
#'     validation to quantify variable importance.
#'
#' Upcoming: Consider multiplying for the noise alignment. This will be values
#'     near 0 won't seem suddenly significant.
#'     Non-noise full model
#'
#' @param data Data.frame of outcome and predictors (PCs and meta-variables).
#'     Note, currently Unit or repeated measures should not be included.
#'     Generally use the results from getPCAData, potentially with meta-
#'     variables attached.
#' @param K (Optional) Numeric indicating the number of folds to use in K-fold
#'     CV. The default is 10.
#' @param outcome (Optional) String indicating the column name with the outcome
#'     in both data and cellData. Default is the first column of data.
#' @param unit (Optional) String indicating the column name with the unit in
#'     both data and cellData. Default is the second column of data.
#' @param repeatedId (Optional) String indicating the column name with the
#'     unique id from repeated measures in both data and cellData. Default is
#'     NULL.
#' @param metaNames (Optional) Vector indicating the meta-variables to be
#'     considered. Default is NULL.
#' @param cellData (Optional) Data.frame indicating the cells used to create
#'     data. It contains outcome, unit, (possibly) repeatedId, agentType.
#' @param syntheticKs (Optional) Numeric indicating the number of variables
#'     in the K noise groups. If the value is 0, no K noise variables are
#'     generated. Default is 100.
#' @param syntheticMetas (Optional) Numeric indicating the number of variables
#'     in the meta noise groups. If the value is 0, no meta noise variables are
#'     generated. Default is 100.
#' @param generalSyntheticK (Optional) Boolean indicating if a general K noise
#'     group should be used or specialized K noise groups should be used.
#'     Default is TRUE.
#' @param alpha (Optional) Numeric in (0,1) indicating the significance used
#'     throughout the analysis. Default is 0.05.
#'
#' This is used in selection of noise variables and in building the CV
#'     intervals.
#'
#' @param silent (Optional) Boolean indicating if output should be suppressed
#'     when the function is running. Default is FALSE.
#' @param ... (Optional) These extra parameters are passed into the
#'     .generateCSRData function.
#'
#' @return List with the following items:
#'     1. Gini: Data.frame with the results of gini indices from the models and
#'                  CV. The columns are var, avg, sd, lower, and upper. The
#'                  columns lower and upper are made with significance alpha.
#'     2. VI: Data.frame with the results of variable importance indices from
#'                the models and CV. The columns are var, avg, sd, lower, and
#'                upper. The columns lower and upper are made with significance
#'                alpha.
#'     3. Accuracy: Data.frame with results of cross validation. The columns are
#'                      avg, sd, lower and upper.
#'     4. varImpPlot: ggplot2 object for a plot of both gini and vi plots. See
#'                        descriptions below.
#'     5. viPlot: ggplot2 object for a plot of vi plot. This will display
#'                    ordered underlying functions and meta-variables with point
#'                    estimates, intervals, and the red (standardized) noise
#'                    cutoff. Values are based on variable importance values.
#'     6. giniPlot: ggplot2 object for a plot of gini plot. This will display
#'                    ordered underlying functions and meta-variables with point
#'                    estimates, intervals, and the red (standardized) noise
#'                    cutoff. Values are based on gini index values.
#'
#' @export
#'
#' @examples
#' data <- simulatePP(cellVarData=
#'                        data.frame('stage'=c(0,1),
#'                                   'A'=c(0,0),
#'                                   'B'=c(1/50,1/50)),
#'                    cellKappaData=data.frame(
#'                                   'cell'=c('A','B'),
#'                                   'clusterCell'=c(NA,'A'),
#'                                   'kappa'=c(20,5)),
#'                    peoplePerStage=100,
#'                    imagesPerPerson=1,
#'                    reduceEdge=0.025,
#'                    silent=F )
#' pcaData <- getPCAData(data,repeatedUniqueId='Image',
#'                           xRange = c(0,1),  yRange = c(0,1), silent=F)
#' pcaMeta <- simulateMeta(pcaData,
#'                 metaInfo = data.frame(
#'                       'var'=c('randUnif','randBin','corrNorm'),
#'                       'rdist'=c('runif','rbinom','rnorm'),
#'                       'Stage_0'=c('0.5','0.5','1'),
#'                       'Stage_1'=c('0.5','0.5','2')))
#' rfcv <- computeRandomForest_CVPC(data=pcaMeta,repeatedId='Image',
#'                                  metaNames=c('randUnif','randBin','corrNorm'),
#'                                  cellData=data)
computeRandomForest_CVPC <- function(data, K=10,
                    outcome=colnames(data)[1],
                    unit=colnames(data)[2],
                    repeatedId=NULL,#colnames(data)
                    metaNames=NULL,
                    cellData=NULL,
                    syntheticKs=100, syntheticMetas=100,
                    generalSyntheticK=T,
                    alpha=0.05, silent=F, ...){

  KFunctions <- .getUnderlyingVariable(
    colnames(data)[!(colnames(data)%in%c(outcome,unit,repeatedId,metaNames))])
  ## Generate Synthetics And Connect
  SyntheticK <- .generateSyntheticKs(syntheticKs = syntheticKs,
                                     generalSyntheticK = generalSyntheticK,
                                     data = data,cellData = cellData,
                                     outcome = outcome, unit = unit,
                                     repeatedId = repeatedId,
                                     metaNames = metaNames, ...)
  SyntheticMeta <- .generateSyntheticMetas(metaNames = metaNames,
                                           syntheticMetas = syntheticMetas,
                                           data = data,
                                           outcome = outcome, unit = unit,
                                           repeatedId = repeatedId)

  data <- merge(merge(data, SyntheticK$pcaData,
              by=colnames(data)[colnames(data)%in%c(outcome,unit,repeatedId)]),
            SyntheticMeta$syntheticMetaData,
              by=colnames(data)[colnames(data)%in%c(outcome,unit,repeatedId)])
  noiseMap <- rbind(SyntheticK$noiseMap, SyntheticMeta$noiseMap)
  # Organize and Prepare
  syntheticMetaNames <- colnames(SyntheticMeta$syntheticMetaData)[
    !(colnames(SyntheticMeta$syntheticMetaData)%in%c(outcome,unit,repeatedId))]
  syntheticKFunctions <- .getUnderlyingVariable(
    colnames(data)[!(colnames(data)%in%c(outcome,unit,repeatedId,
                                         metaNames,syntheticMetaNames))])
  syntheticKFunctions <- syntheticKFunctions[!(syntheticKFunctions %in% KFunctions)]
  allKFunctionsMetas <- c(KFunctions, syntheticKFunctions,
                          metaNames, syntheticMetaNames)

  avgVI <- avgGini <- data.frame('var'=allKFunctionsMetas)
  oobAcc <- rep(NA,K)
  groups <- .getFolds(1:nrow(data), K)

  ## Fit Models
  # TODO:: Update this without noise
  # if(!silent)
  #   cat(paste0('Full Model\n'))
  # model <- .computeRandomForest_PC(data=data,
  #                                  outcome = outcome,
  #                                  unit=unit, repeatedId=repeatedId,
  #                                  varImpPlot = F,
  #                                  metaNames=c(metaNames, syntheticMetaNames))

  # K-Fold Cross-Validation
  for(i in 1:K){
    if(!silent)
      cat(paste0('CV Trial: ',i,'/',K,'\n'))
    RF <- .computeRandomForest_PC(data=data[-groups[[i]],],
                                  outcome = outcome,
                                  unit=unit, repeatedId=repeatedId,
                                  varImpPlot = F,
                                  metaNames=c(metaNames, syntheticMetaNames))

    data_merge <- RF[[2]][,c('var','avgGini')]
    colnames(data_merge) <- c('var',paste0('avgGiniK',i))
    avgGini <- merge(avgGini, data_merge, by='var')

    data_merge <- RF[[2]][,c('var','avgVI')]
    colnames(data_merge) <- c('var',paste0('avgVIK',i))
    avgVI <- merge(avgVI, data_merge, by='var')

    oobAcc[i] <- sum(data[groups[[i]],'Stage']==
                       .predict.RandomForest_PC(model = RF[[1]],
                                data_pred = data[groups[[i]],],
                                type = 'pred', data = data)) /
                  nrow(data[groups[[i]],])
  }

  # Organize and return
  varImpList <- .getVariableImportanceMetrics(avgGini, avgVI, oobAcc, alpha)

  append(
  append(
    list('Gini'=varImpList$giniData,
         'VI'=varImpList$viData,
         'Accuracy'=varImpList$oobAccData),
    .generateCVVariableImportancePlot(giniData=varImpList$giniData,
                                    viData=varImpList$viData,
                                    oobAccData=varImpList$oobAccData,
                                    noiseMap=noiseMap,
                                    KFunctions=KFunctions,metaNames=metaNames,
                                    alpha=alpha)
    ),
  list('a'=.generateCVVariableImportancePlot_mult(giniData=varImpList$giniData,
                                         viData=varImpList$viData,
                                         oobAccData=varImpList$oobAccData,
                                         noiseMap=noiseMap,
                                         KFunctions=KFunctions,metaNames=metaNames,
                                         alpha=alpha))
  )
}


#' Generate Synthetic K functions
#'
#' This (internal) function generates the synthetic K functions (in PCs).
#'
#' Upcoming: See following experiments:
#'     - Test if Kappa should be related to data or not
#'     - Determine if window size c(1,1) is reasonable for all situation
#'
#' @param syntheticKs Numeric indicating the number of variables
#'     in the K noise groups. If the value is 0, no K noise variables are
#'     generated.
#' @param data Data.frame of outcome and predictors (PCs and meta-variables).
#'     Note, currently Unit or repeated measures should not be included.
#'     Generally use the results from getPCAData, potentially with meta-
#'     variables attached.
#' @param cellData Data.frame indicating the cells used to create data. It
#'     contains outcome, unit, (possibly) repeatedId, agentType.
#' @param outcome String indicating the column name with the outcome in both
#'     data and cellData.
#' @param unit String indicating the column name with the unit in both data and
#'     cellData.
#' @param repeatedId String indicating the column name with the unique id from
#'     repeated measures in both data and cellData.
#' @param metaNames Vector indicating the meta-variables to be considered.
#' @param generalSyntheticK Boolean indicating if a general K noise group
#'     should be used or specialized K noise groups should be used.
#' @param ... (Optional) Parameters sent to .generateKNoiseGroup. These extra
#'     parameters are passed into the .generateCSRData function.
#'
#' @return List with two elements:
#'         1. noiseMap: Data.frame for matching noise variables to noise groups.
#'                          The columns are noiseVar (the synthetic interactions
#'                          of which there are syntheticKs), dataVar (name of
#'                          noise group), and synType ('K', which will compare
#'                          to 'Meta' from other functions).
#'         2. pcaData: Data.frame with nPCs equal to that in data. This contains
#'                         the pcaData for the noise variables. See values
#'                         returns from getPCAData, but outcome, person,
#'                         repeatedId, and PCs for synthetic Ks.
#'
#' @export
#'
#' @examples
#' # See code for computeRandomForest_CVPC. This is not an outward function so
#' #     won't be viewable.
.generateSyntheticKs <- function(syntheticKs, data, cellData,
                                 outcome, unit, repeatedId, metaNames,
                                 generalSyntheticK, ...){
  # Get Underlying K function names and number of PCs
  NamesPCs <- .getUnderlyingNamesAndPCsCount(
                colnames(data)[!(colnames(data)%in%
                                 c(outcome, unit, repeatedId, metaNames))])

  # Return if no synthetics need to be generated
  if(syntheticKs<1 || is.null(NamesPCs))
    return()

  noiseMap <- NULL

  if(generalSyntheticK){
    if(!is.null(cellData[['cellType']])){
      kappa <- round(mean(table(cellData[,'cellType'])) /
                       nrow(unique(cellData[,c(outcome, unit, repeatedId)])))
    } else{
      kappa <- (sum(cellData[
                      unique(stringr::str_remove_all(NamesPCs$Vars,'_.*$'))]) /
                  length(NamesPCs$Vars) ) /
                  nrow(unique(cellData[,c(outcome, unit, repeatedId)]))
    }

    noiseData <- .generateKNoiseGroup(
                    syntheticKs = syntheticKs, cellData = cellData,
                    agentName='syntheticCell',
                    kappaL = kappa, kappaR=kappa,
                    outcome = outcome, unit = unit,
                    repeatedId = repeatedId, ...)
    # Interactions to check
    agent_df <- data.frame('c1'=paste0('syntheticCell',1:syntheticKs,'L'),
                           'c2'=paste0('syntheticCell',1:syntheticKs,'R'))
    noiseMap <- data.frame(
              'noiseVar'=paste(agent_df[,1],agent_df[,2],sep='_'),
              'dataVar'='syntheticK',
              'synType'='K')
  }else{
    # Setup agent counts
    agentCounts <- as.data.frame(table(cellData[,'cellType']))
    nrowUnique <- nrow(unique(cellData[,c(outcome, unit, repeatedId)]))
    agentCounts$measureCounts <- agentCounts$Freq/nrowUnique

    # Other prep
    noiseData <- cell_df <- NULL
    agentL <- agentR <- rep(NA, length(NamesPCs$Vars))
    for(i in 1:length(NamesPCs$Vars)){
      # Determine mean of each agent
      PCagents <- unlist(strsplit(NamesPCs$Vars[1],'_'))
      kapL <- agentCounts[agentCounts$Var1==PCagents[1],'measureCounts']
      kapR <- agentCounts[agentCounts$Var1==PCagents[2],'measureCounts']

      noiseData <- rbind(noiseData,
                      .generateKNoiseGroup(
                         syntheticKs = syntheticKs, cellData = cellData,
                         agentName=paste0('Synthetic',NamesPCs$Vars[i]),
                         kappaL = kapL, kappaR = kapR,
                         outcome = outcome, unit = unit,
                         repeatedId = repeatedId, ...))

      # Interactions to check
      agent_df <- rbind(agent_df,
                data.frame(
                  'c1'=paste0('Synthetic',NamesPCs$Vars[i],1:syntheticKs,'L'),
                  'c2'=paste0('Synthetic',NamesPCs$Vars[i],1:syntheticKs,'R')))
      noiseMap <- rbind(noiseMap,
                        data.frame(
                          'noiseVar'=paste(agent_df[,1],agent_df[,2],sep='_'),
                          'dataVar'=NamesPCs$Vars[i],
                          'synType'='K'))
    }
  }

  list('noiseMap'=noiseMap,
       'pcaData'= getPCAData(data=noiseData, nPCs = NamesPCs$PCs,
                    outcome = outcome,unit=unit,repeatedUniqueId = repeatedId,
                    xRange = c(0,1),yRange = c(0,1), agents_df = agent_df, ...)
  )
}


#' Generate K Noise Group
#'
#' This (internal) function generates noise groups of syntheticKs patterns.
#'
#' @param syntheticKs Numeric indicating the number of variables
#'     in the K noise groups. If the value is 0, no K noise variables are
#'     generated.
#' @param cellData Data.frame indicating the cells used to create data. It
#'     contains outcome, unit, (possibly) repeatedId, and agentType.
#' @param kappaL Numeric value for intensity of left agent (A in traditional A_A
#'     writing of interaction terms).
#' @param kappaR Numeric value for intensity of right agent (B in traditional
#'     A_B writing of interaction terms).
#' @param agentName String indicating the K function it is matched with. This
#'     may be a general term.
#' @param outcome String indicating the column name with the outcome in
#'     cellData.
#' @param unit String indicating the column name with the unit in cellData.
#' @param repeatedId String indicating the column name with the unique id from
#'     repeated measures in cellData.
#' @param ... (Optional) Parameters sent to .generateKNoiseGroup. These extra
#'     parameters are passed into the .generateNoiseVariable function.
#'
#' @return Data.frame for all the cells in the noise group. This will contain
#'     columns of x, y, and cellType.
#' @export
#'
#' @examples
#' # See code for .generateSyntheticKs. This is not an outward function so won't
#' #     be viewable.
.generateKNoiseGroup <- function(syntheticKs,cellData,kappaL,kappaR,agentName,
                                 outcome, unit, repeatedId, ...){
  allSituations <- unique(cellData[,c(outcome, unit, repeatedId)])
  syntheticData <- NULL

  # Generate Cells
  for(i in 1:nrow(allSituations)){
    syntheticData <- rbind(syntheticData,
                           .generateNoiseVariable(
                                  syntheticKs = syntheticKs,
                                  kappaL = kappaL, kappaR = kappaR,
                                  currentInfo = allSituations[i,],
                                  agentName = agentName,
                                  outcome = outcome, unit = unit,
                                  repeatedId = repeatedId, ...))
  }

  syntheticData
}


#' Generate Noise Variable
#'
#' This (internal) function generates the cells for the syntheticKs K functions
#'     --two patterns for each side.
#'
#' @param syntheticKs Numeric indicating the number of variables
#'     in the K noise groups. If the value is 0, no K noise variables are
#'     generated.
#' @param kappaL Numeric value for intensity of left agent (A in traditional A_A
#'     writing of interaction terms).
#' @param kappaR Numeric value for intensity of right agent (B in traditional
#'     A_B writing of interaction terms).
#' @param currentInfo Data.frame with the information for the patterns. Includes
#'     outcome, unit, and (potentially) repeatedId values.
#' @param agentName String indicating the K function it is matched with. This
#'     may be a general term.
#' @param outcome String indicating the column name with the outcome in
#'     cellData.
#' @param unit String indicating the column name with the unit in cellData.
#' @param repeatedId String indicating the column name with the unique id from
#'     repeated measures in cellData.
#' @param ... (Optional) These extra parameters are passed into the
#'     .generateCSRData function.
#'
#' @return Data.frame with cells for noise variables. The columns are x, y, and
#'     cellType.
#' @export
#'
#' @examples
#' # See code for .generateKNoiseGroup. This is not an outward function so won't
#' #     be viewable.
.generateNoiseVariable <- function(syntheticKs, kappaL, kappaR, currentInfo,
                                   agentName, outcome, unit, repeatedId, ...){
  syntheticCells <- NULL

  for(p in 1:syntheticKs){
    # Generate syntheticCellsL
    syntheticCellsL <- .generateCSRData(xRange = c(0,1), yRange = c(0,1),
                                    kappa = kappaL,
                                    cellType = paste0(agentName,p,'L'),
                                    ...)
    syntheticCellsL[[outcome]] <- currentInfo[[outcome]]
    syntheticCellsL[[unit]] <- currentInfo[[unit]]
    if(!is.null(repeatedId))
      syntheticCellsL[[repeatedId]] <- currentInfo[[repeatedId]]

    # Generate syntheticCellsR
    syntheticCellsR <- .generateCSRData(xRange = c(0,1), yRange = c(0,1),
                                    kappa = kappaR,
                                    cellType = paste0(agentName,p,'R'),
                                    ...)
    syntheticCellsR[[outcome]] <- currentInfo[[outcome]]
    syntheticCellsR[[unit]] <- currentInfo[[unit]]
    if(!is.null(repeatedId))
      syntheticCellsR[[repeatedId]] <- currentInfo[[repeatedId]]

    syntheticCells <- rbind(syntheticCells, syntheticCellsL, syntheticCellsR)
  }

  syntheticCells
}


#' Get Underlying Names and Principal Components Counts
#'
#' This (internal) function get the underlying functions and the total principal
#'     components from the variable function names.
#'
#' @param namesOfInterest Vector of variable names, names with _PC# with them
#'
#' @return List with two elements:
#'     1. Vars: Underlying function names
#'     2. PCs: Numeric indicating the number of principal components
#' @export
#'
#' @examples
#' # See code for .generateSyntheticKs. This is not an outward function so won't
#' #     be viewable.
.getUnderlyingNamesAndPCsCount <- function(namesOfInterest){
  splitNames <- stringr::str_split(namesOfInterest, "_PC(?=[0-9]$)")
  tryCatch({
    return(list('Vars'=unique(sapply(splitNames,'[[',1)),
                'PCs'=max(as.numeric(sapply(splitNames,'[[',2)))))
  },error=function(e){
    return()
  }, warning=function(w){
    return()
  })
}


#' Generate Synthetic Metas
#'
#' This (internal) function generates the synthetic meta-variables and also
#'     returns the noiseMap between noise meta variables and underlying meta
#'     variables.
#'
#' @param metaNames Vector of meta variable names as columns in data.
#' @param syntheticMetas Numeric indicating the number of synthetic meta
#'     variables in each noise group.
#' @param data Data.frame of outcome and predictors (PCs and meta-variables).
#'     Note, currently Unit or repeated measures should not be included.
#'     Generally use the results from getPCAData, potentially with meta-
#'     variables attached.
#' @param outcome String indicating the column name with the outcome in
#'     cellData.
#' @param unit String indicating the column name with the unit in cellData.
#' @param repeatedId String indicating the column name with the unique id from
#'     repeated measures in cellData.
#'
#' @return List with two elements:
#'     1. noiseMap: Data.frame for matching noise variables to noise groups. The
#'                  columns are noiseVar (the synthetic interactions of which
#'                  there are syntheticMetas), dataVar (name of noise group),
#'                  and synType ('Meta', which will compare to 'K' from other
#'                  functions).
#'     2. syntheticMetaData: Data.frame with results and synthetic meta
#'                           variables. This contains the outcome, unit, and
#'                           (potentially) repeatedId along with synthetic meta
#'                           variables.
#' @export
#'
#' @examples
#' # See code for computeRandomForest_CVPC. This is not an outward function so
#' #     won't be viewable.
.generateSyntheticMetas <- function(metaNames,syntheticMetas,data,
                                    outcome,unit,repeatedId){

  syntheticMetaData <- data[colnames(data) %in% c(outcome,unit,repeatedId)]
  noiseMap <- NULL

  if(is.null(metaNames) || syntheticMetas<1) {
    return(list('noiseMap'=noiseMap,
                'syntheticMetaData'=syntheticMetaData))
  }

  # Bootstrap each meta-variable
  for(i in 1:length(metaNames)){
    syntheticMetaData <- cbind(syntheticMetaData,
              .generateMetaNoiseGroup(data, metaNames[i],syntheticMetas))

    noiseMap <- rbind(noiseMap,
                      data.frame(
                        'noiseVar'=paste0('SyntheticMeta',1:syntheticMetas,
                                          '_',metaNames[i]),
                        'dataVar'=metaNames[i],
                        'synType'='Meta'))
  }

  list('noiseMap'=noiseMap, 'syntheticMetaData'=syntheticMetaData)
}


#' Generate Meta Noise Group
#'
#' This (internal) function that will bootstrap meta variables to create
#'     synthetic variables.
#'
#' @param data Data.frame of outcome and predictors (PCs and meta-variables).
#'     Note, currently Unit or repeated measures should not be included.
#'     Generally use the results from getPCAData, potentially with meta-
#'     variables attached. Really only column with meta column needed
#' @param meta String of the column name indicating the meta variable to be
#'     generated.
#' @param syntheticMetas Numeric indicating the number of synthetic meta
#'     variables in each noise group.
#'
#' @return Data.frame with syntheticMetas columns and observations (rows) of
#'     data. Each column is a bootstrapped sample of the meta variable.
#' @export
#'
#' @examples
#' # See code for .generateSyntheticMetas. This is not an outward function so
#' #     won't be viewable.
.generateMetaNoiseGroup <- function(data, meta, syntheticMetas){
  noiseData <- NULL
  for(m in 1:syntheticMetas){
    noiseData <- cbind(noiseData, sample(data[[meta]],replace = T))
  }
  colnames(noiseData) <- paste0('SyntheticMeta',1:syntheticMetas,'_',meta)

  noiseData
}


#' Generate Variable Importance Plots
#'
#' This (internal) function generates the CV variable importance plots for the
#'     RandomForest_CVPC function. It add CV errors and manages error alignment.
#'
#' Upcoming: Remove Gini.
#'           Remove upper1/lower1
#'
#' @param giniData Data.frame
#' @param viData Data.frame
#' @param oobAccData Data.frame
#' @param noiseMap Data.frame
#' @param KFunctions Vector of the underlying K functions of interest. These
#'     will be the ones plotted.
#' @param metaNames Vector of the underlying Meta-variables of interest. These
#'     will be the ones plotted.
#' @param alpha Numeric in (0,1) indicating the significance for the cutoff.
#'
#' @return List containing three plots.
#'     1. varImpPlot: giniPlot and viPlot together with cross-validated error.
#'     2. viPlot: plot with each K function and Meta-variable aligned with
#'                    noise variable cutoff using variable importance metric.
#'     3. giniPlot: plot with each K function and Meta-variable aligned with
#'                    noise variable cutoff using gini metric.
#' @export
#'
#' @examples
#' # See code for computeRandomForest_CVPC. This is not an outward function so
#' #     won't be viewable.
.generateCVVariableImportancePlot <- function(giniData, viData, oobAccData,
                                              noiseMap, KFunctions, metaNames,
                                              alpha ){
  # Organize Noise
  giniCutoff_df <- viCutoff_df <- unique(noiseMap[,c('dataVar','synType')])
  if(nrow(giniCutoff_df)){
    giniCutoff_df$cutoff <- viCutoff_df$cutoff <- NA

    # Get all alpha confidence of the related noise
    for(i in 1:nrow(giniCutoff_df)){

      giniCutoff_df$cutoff[i] <-
        quantile(giniData[giniData$var %in%
                            noiseMap[noiseMap$dataVar==giniCutoff_df$dataVar[i],
                                     'noiseVar'],
                          'avg'],1-alpha)[[1]]
      viCutoff_df$cutoff[i] <-
        quantile(viData[viData$var %in%
                          noiseMap[noiseMap$dataVar==viCutoff_df$dataVar[i],
                                   'noiseVar'],
                        'avg'],1-alpha)[[1]]
    }

    giniCutoff <- max(giniCutoff_df$cutoff)
    viCutoff <- max(viCutoff_df$cutoff)
    giniCutoff_df$adj <- giniCutoff-giniCutoff_df$cutoff
    viCutoff_df$adj <- viCutoff-viCutoff_df$cutoff
  }

  # Standardize
  underlyingGini <- underlyingVI <- data.frame('var'=c(KFunctions, metaNames),
                                               'avgVal'=NA,
                                               'lower'=NA, 'upper'=NA,
                                               'lower1'=NA, 'upper1'=NA)
  for(uVar in c(KFunctions, metaNames)){
    if(uVar %in% giniCutoff_df$dataVar) {
      giniAdjVal <- giniCutoff_df[giniCutoff_df$dataVar==uVar,'adj']
      viAdjVal <- viCutoff_df[viCutoff_df$dataVar==uVar,'adj']
    }else if(uVar %in% KFunctions){
      giniAdjVal <- giniCutoff_df[giniCutoff_df$dataVar=='syntheticK','adj']
      viAdjVal <- viCutoff_df[viCutoff_df$dataVar=='syntheticK','adj']
    }else if(uVar %in% metaNames){
      giniAdjVal <- giniCutoff_df[giniCutoff_df$dataVar=='syntheticMeta','adj']
      viAdjVal <- viCutoff_df[viCutoff_df$dataVar=='syntheticMeta','adj']
    }

    underlyingGini[underlyingGini$var==uVar,2:4] <- giniAdjVal +
      giniData[giniData$var==uVar,c('avg','lower','upper')]
    underlyingVI[underlyingVI$var==uVar,2:4] <- viAdjVal +
      viData[viData$var==uVar,c('avg','lower','upper')]

  }


  # Generate figures
  giniPlot <- .plotCVVariableImportance(underlyingGini,'Gini',giniCutoff)
  viPlot <- .plotCVVariableImportance(underlyingVI,'Variable Importance',viCutoff)

  varImpPlot <- gridExtra::arrangeGrob(giniPlot,viPlot,
                                       layout_matrix = rbind(c(1,2)),
                                       bottom = paste0('Variable Importance - OOB-Acc ',
                                                       min(1,max(0,.specify_decimal(oobAccData$avg, 2))),
                                                       ' (',max(0,.specify_decimal(oobAccData$lower,2)),'-',
                                                       min(1,.specify_decimal(oobAccData$upper,2)),')' ))
  giniPlot <- giniPlot+ggplot2::ylab(paste0('Gini - OOB-Acc ',
                                            min(1,max(0,.specify_decimal(oobAccData$avg, 2))),
                                            ' (', max(0,.specify_decimal(oobAccData$lower,2)),'-',
                                            min(1,.specify_decimal(oobAccData$upper,2)),')' ))
  viPlot <- viPlot+ggplot2::ylab(paste0('Variable Importance - OOB-Acc ',
                                        min(1,max(0,.specify_decimal(oobAccData$avg, 2))),
                                        ' (', max(0,.specify_decimal(oobAccData$lower,2)),'-',
                                        min(1,.specify_decimal(oobAccData$upper,2)),')' ))

  list('varImpPlot'=varImpPlot,'viPlot'=viPlot,'giniPlot'=giniPlot)
}


.generateCVVariableImportancePlot_mult <- function(giniData, viData, oobAccData,
                                              noiseMap, KFunctions, metaNames,
                                              alpha ){
  # Organize Noise
  giniCutoff_df <- viCutoff_df <- unique(noiseMap[,c('dataVar','synType')])
  if(nrow(giniCutoff_df)){
    giniCutoff_df$cutoff <- viCutoff_df$cutoff <- NA

    # Get all alpha confidence of the related noise
    for(i in 1:nrow(giniCutoff_df)){

      giniCutoff_df$cutoff[i] <-
        quantile(giniData[giniData$var %in%
                            noiseMap[noiseMap$dataVar==giniCutoff_df$dataVar[i],
                                     'noiseVar'],
                          'avg'],1-alpha)[[1]]
      viCutoff_df$cutoff[i] <-
        quantile(viData[viData$var %in%
                          noiseMap[noiseMap$dataVar==viCutoff_df$dataVar[i],
                                   'noiseVar'],
                        'avg'],1-alpha)[[1]]
    }

    giniCutoff <- max(giniCutoff_df$cutoff)
    viCutoff <- max(viCutoff_df$cutoff)
    giniCutoff_df$adj <- giniCutoff/giniCutoff_df$cutoff
    viCutoff_df$adj <- viCutoff/viCutoff_df$cutoff
  }

  # Standardize
  underlyingGini <- underlyingVI <- data.frame('var'=c(KFunctions, metaNames),
                                               'avgVal'=NA,
                                               'lower'=NA, 'upper'=NA)
  for(uVar in c(KFunctions, metaNames)){
    if(uVar %in% giniCutoff_df$dataVar) {
      giniAdjVal <- giniCutoff_df[giniCutoff_df$dataVar==uVar,'adj']
      viAdjVal <- viCutoff_df[viCutoff_df$dataVar==uVar,'adj']
    }else if(uVar %in% KFunctions){
      giniAdjVal <- giniCutoff_df[giniCutoff_df$dataVar=='syntheticK','adj']
      viAdjVal <- viCutoff_df[viCutoff_df$dataVar=='syntheticK','adj']
    }else if(uVar %in% metaNames){
      giniAdjVal <- giniCutoff_df[giniCutoff_df$dataVar=='syntheticMeta','adj']
      viAdjVal <- viCutoff_df[viCutoff_df$dataVar=='syntheticMeta','adj']
    }

    underlyingGini[underlyingGini$var==uVar,2:4] <-
      giniData[giniData$var==uVar,c('avg','lower','upper')]*giniAdjVal
    underlyingVI[underlyingVI$var==uVar,2:4] <-
      viData[viData$var==uVar,c('avg','lower','upper')]*viAdjVal

  }


  # Generate figures
  giniPlot <- .plotCVVariableImportance(underlyingGini,'Gini',giniCutoff)
  viPlot <- .plotCVVariableImportance(underlyingVI,'Variable Importance',viCutoff)

  varImpPlot <- gridExtra::arrangeGrob(giniPlot,viPlot,
                                       layout_matrix = rbind(c(1,2)),
                                       bottom = paste0('Variable Importance - OOB-Acc ',
                                                       min(1,max(0,.specify_decimal(oobAccData$avg, 2))),
                                                       ' (',max(0,.specify_decimal(oobAccData$lower,2)),'-',
                                                       min(1,.specify_decimal(oobAccData$upper,2)),')' ))
  giniPlot <- giniPlot+ggplot2::ylab(paste0('Gini - OOB-Acc ',
                                            min(1,max(0,.specify_decimal(oobAccData$avg, 2))),
                                            ' (', max(0,.specify_decimal(oobAccData$lower,2)),'-',
                                            min(1,.specify_decimal(oobAccData$upper,2)),')' ))
  viPlot <- viPlot+ggplot2::ylab(paste0('Variable Importance - OOB-Acc ',
                                        min(1,max(0,.specify_decimal(oobAccData$avg, 2))),
                                        ' (', max(0,.specify_decimal(oobAccData$lower,2)),'-',
                                        min(1,.specify_decimal(oobAccData$upper,2)),')' ))

  list('varImpPlot'=varImpPlot,'viPlot'=viPlot,'giniPlot'=giniPlot)
}


#' Get Variable Importance Metrics
#'
#' @param avgGini Data.frame obtained from CV sample of random forest. Columns
#'     of var (K functions and metas, including synthetics) and average gini
#'     index value from each CV random forest.
#'     variable
#' @param avgVI Data.frame obtained from CV sample of random forest. Columns
#'     of var (K functions and metas, including synthetics) and average variable
#'     importance value from each CV random forest.
#' @param oobAcc Vector of numerics for out-of-bag accuracy from each fold.
#'     Numerics are thus in (0,1).
#' @param alpha (Optional) Numeric in (0,1) indicating the significance for the
#'     upper and lower values (quantiles).
#'
#' @return List with three elements:
#'     1. giniData: Data.frame organizing the gini index metrics of the
#'                  variables. Columns var, avg, sd, lower, and upper.
#'     2. viData: Data.frame organizing the variable importance index metrics of
#'                the variables. Columns var, avg, sd, lower, and upper.
#'     3. oobAccData: Data.frame organizing the out-of-bag accuracy of the
#'                    models (only 1 row). Columns avg, sd, lower, and upper.
#' @export
#'
#' @examples
#' # See code for computeRandomForest_CVPC. This is not an outward function so
#' #     won't be viewable.
.getVariableImportanceMetrics <- function(avgGini, avgVI, oobAcc, alpha=0.05){
  # Variable Importance values (Empirical close to normal)
  giniData <- data.frame('var'=avgGini$var,
                         'avg'=rowMeans(avgGini[-1]),
                         'sd'=apply(avgGini[-1],MARGIN=1,FUN=function(x){sd(x)}),
                         'lower'= apply(avgGini[-1], 1, quantile, probs = alpha/2),
                         'upper'= apply(avgGini[-1], 1, quantile, probs = 1-alpha/2))
  viData <- data.frame('var'=avgVI$var,
                       'avg'=rowMeans(avgVI[-1]),
                       'sd'=apply(avgVI[-1],MARGIN=1,FUN=function(x){sd(x)}),
                       'lower'= apply(avgVI[-1], 1, quantile, probs = alpha/2),
                       'upper'= apply(avgVI[-1], 1, quantile, probs = 1-alpha/2))
  # OOB (Length K so use normal approx rather than empirical percent)
  oobAccMeans <- mean(oobAcc)
  oobAccSD <- sd(oobAcc)
  oobAccData <- data.frame('avg'=oobAccMeans,
                           'sd'=oobAccSD,
                           'lower'= oobAccMeans +qnorm(alpha) * oobAccSD/sqrt(length(oobAcc)),
                           'upper'= oobAccMeans -qnorm(alpha) * oobAccSD/sqrt(length(oobAcc)))

  list('giniData'=giniData,'viData'=viData,'oobAccData'=oobAccData)
}


#' Plot CV Variable Importance
#'
#' This (internal) function produces the CV variable importance plots with the
#'     noise cutoff and CV intervals.
#'
#' Upcoming: Rename ylabel (as we flip cooridinates its actually x-axis we see)
#'
#' @param underlyingData Data.frame
#' @param ylabel String for the y-axis label
#' @param cutoff Numeric in [0,1] for the red noise line
#'
#' @return ggplot2 plot for the CV variable importance plot
#' @export
#'
#' @examples
#' # See code for .generateCVVariableImportancePlot. This is not an outward
#' #     function so won't be viewable.
.plotCVVariableImportance <- function(underlyingData, ylabel,cutoff){
  maxVal <- max(underlyingData$avgVal)
  ggplot2::ggplot(data=underlyingData,
         mapping=ggplot2::aes(x=reorder(var,avgVal),
                     y=ifelse(avgVal/maxVal>1,1,
                              ifelse(avgVal/maxVal<0,0,
                                     avgVal/maxVal)))) +
    ggplot2::geom_point(color='black') +
    ggplot2::geom_errorbar(ggplot2::aes(ymin=ifelse(lower/maxVal<0,0,lower/maxVal),
                      ymax = ifelse(upper/maxVal>1,1,upper/maxVal)),
                  color='black', width=0.2) +
    ggplot2::geom_hline(ggplot2::aes(yintercept=cutoff/maxVal),
               color='red', linetype='dashed') +
    ggplot2::coord_flip() +
    ggplot2::xlab(NULL) +
    ggplot2::ylim(c(0,1)) +
    ggplot2::ylab(ylabel) +
    ggplot2::theme_bw()
}
