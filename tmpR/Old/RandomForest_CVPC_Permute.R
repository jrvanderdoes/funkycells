#' Fit a Random Forest model with PC data (Using CV for Improvements)
#'
#' The function fits a random forest model to the data along with using cross-
#'     validation to quantify variable importance.
#'
#' Warning: If there are no sythentics, this may break (will fix it eventually)
#'
#' @param data Data.frame of outcome and predictors (PCs and meta-variables).
#'     Generally use the results from getKsPCAData, potentially with meta-
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
#' @param curvedSigSims (Optional) Numeric indicating the number of simulations
#'     used to create the curved lines. Default is 100.
#' @param alpha (Optional) Numeric in (0,1) indicating the significance used
#'     throughout the analysis. Default is 0.05.
#'
#' This is used in selection of noise variables and in building the CV
#'     intervals.
#'
#' @param silent (Optional) Boolean indicating if output should be suppressed
#'     when the function is running. Default is FALSE.
#' @param rGuessSims (Optional) Numeric value indicating the number of
#'     simulations used for guessing and creating the guess estimate on the
#'     plot. Default is 500.
#' @param alignmentMethod (Optional) String indicating the method of aligning
#'     variables via noise variables. The options are 'Add', 'Mult', or c('Add',
#'     'Mult'). The default value is c('Add','Mult').
#' @param subsetPlotSize (Optional) Numeric indicating the number of top
#'     variables to include in a subset graph (note if there are less variables)
#'     than this value indicates then no subset graph will be produced. Default
#'     is 25.
#' @param nTrees (Optional) Numeric indicating the number of trees in each
#'     forest. The default is 500.
#'
#' @return List with the following items:
#'     \enumerate{
#'         \item Gini: Data.frame with the results of gini indices from the
#'                     models and CV. The columns are var, avg, sd, lower, and
#'                     upper. The columns lower and upper are made with
#'                     significance alpha.
#'         \item VI: Data.frame with the results of variable importance indices
#'                   from the models and CV. The columns are var, avg, sd,
#'                   lower, and upper. The columns lower and upper are made with
#'                   significance alpha.
#'         \item Accuracy: Data.frame with results of cross validation. The
#'                         columns are avg, sd, lower and upper.
#'         \item NoiseCurve: (Optional) Contains columns for noise curve (orange)
#'         \item varImpPlot: ggplot2 object (may be in list or seperate) for a
#'                           plot of both gini and vi plots. See following
#'                           descriptions.
#'               \enumerate{
#'                   \item viPlot: ggplot2 object for a plot of vi plot. This
#'                                 will display ordered underlying functions and
#'                                 meta-variables with point estimates,
#'                                 intervals, and the red (standardized) noise
#'                                 cutoff. Values are based on variable
#'                                 importance values.
#'                   \item giniPlot: ggplot2 object for a plot of gini plot.
#'                                   This will display ordered underlying
#'                                   functions and meta-variables with point
#'                                   estimates, intervals, and the red
#'                                   (standardized) noise cutoff. Values are
#'                                   based on gini index values.
#'               }
#'         \item Subset (Optional) These next parts are similar to the previous
#'               part, but subset figures. Only the top subsetPlotSize number of
#'               variables are displayed in the graph - for better seeing
#'               interesting patterns.
#'            }
#' @export
#'
#' @examples
#' dat <- simulatePP(cellVarData=
#'             data.frame('stage'=c(0,1), 'A'=c(0,0), 'B'=c(1/50,1/50)),
#'             cellKappaData=data.frame(
#'                   'cell'=c('A','B'),
#'                   'clusterCell'=c(NA,'A'),
#'                   'kappa'=c(20,5)),
#'             peoplePerStage=100,
#'             imagesPerPerson=1,
#'             silent=F )
#' pcaData <- getKsPCAData(dat,repeatedUniqueId='Image',
#'                       xRange = c(0,1),  yRange = c(0,1), silent=F)
#' pcaMeta <- simulateMeta(pcaData,
#'                         metaInfo = data.frame(
#'                            'var'=c('randUnif','randBin','corrNorm'),
#'                            'rdist'=c('runif','rbinom','rnorm'),
#'                            'Stage_0'=c('0.5','0.5','1'),
#'                            'Stage_1'=c('0.5','0.5','2')))
#' rfcv <- computeRandomForest_CVPC_Permute(data=pcaMeta,repeatedId='Image',
#'                             metaNames=c('randUnif','randBin','corrNorm'),
#'                             cellData=dat)
#'
#' dat <- simulatePP(cellVarData=
#'                      data.frame('stage'=c(0,1), 'A'=c(0,0), 'B'=c(1/50,1/100)),
#'                  cellKappaData=data.frame(
#'                      'cell'=c('A','B'),
#'                      'clusterCell'=c(NA,'A'),
#'                      'kappa'=c(20,5)),
#'                  peoplePerStage=100,
#'                  imagesPerPerson=1,
#'                  silent=F )
#' pcaData <- getKsPCAData(dat,repeatedUniqueId='Image',
#'                       xRange = c(0,1),  yRange = c(0,1), silent=F)
#' pcaMeta <- simulateMeta(pcaData,
#'                         metaInfo = data.frame(
#'                             'var'=c('randUnif','randBin','corrNorm'),
#'                             'rdist'=c('runif','rbinom','rnorm'),
#'                             'Stage_0'=c('0.5','0.5','1'),
#'                             'Stage_1'=c('0.5','0.5','2')))
#' rfcv <- computeRandomForest_CVPC_Permute(data=pcaMeta,repeatedId='Image',
#'                     metaNames=c('randUnif','randBin','corrNorm'),
#'                     cellData=dat)
#' #funkyRandomForest
computeRandomForest_CVPC_Permute <- function(data, K=10,
                                             outcome=colnames(data)[1],
                                             unit=colnames(data)[2],
                                             metaNames=NULL,
                                             synthetics=100,
                                             alpha=0.05,
                                             silent=F,
                                             rGuessSims=500,
                                             subsetPlotSize=25, nTrees=500){
  ## From moving code over, will remove
  repeatedId=NULL

  ## Error checking
  #.checkData(alignmentMethod)

  ## Generate Synthetics And Connect
  components <- colnames(data)[!(colnames(data)%in%c(outcome,unit,repeatedId,metaNames))]
  KFunctions <- .getUnderlyingVariable(components)
  # Get Var and data column alignmnet
  underlyingDataAlignedFunctions <- .getUnderlyingVariable(colnames(data),
                                                           returnUnique = F)
  underlyingVars <- unique(underlyingDataAlignedFunctions)
  underlyingVars <- underlyingVars[!(underlyingVars %in%
                                       c(outcome,unit,repeatedId))]

  avgVI <- data.frame('var'=underlyingVars)
  oobAcc <- rep(NA,K)
  groups <- .getFolds(1:nrow(data), K)

  # K-Fold Cross-Validation (True Data)
  if(!silent) cat('CV Trials (',K,'): ',sep='')
  for(i in 1:K){
    if(!silent) cat(i,', ',sep='')
    RF <- computeRandomForest_PC(data=data[-groups[[i]],],
                                 outcome = outcome,
                                 unit=unit, repeatedId=repeatedId,
                                 varImpPlot = F,
                                 metaNames=c(metaNames),
                                 nTrees=nTrees)

    data_merge <- RF$varImportanceData[,c('var','avgVI')]
    colnames(data_merge) <- c('var',paste0('avgVIK',i))
    avgVI <- merge(avgVI, data_merge, by='var')

    oobAcc[i] <- sum(data[groups[[i]],outcome]==
                       predict.RandomForest_PC(model = RF$model,
                                               data_pred = data[groups[[i]],],
                                               type = 'pred', data = data)) /
      nrow(data[groups[[i]],])
  }
  if(!silent) cat('\n')

  # Average on Full for estimate
  avgVI_old <- data.frame('var'=underlyingVars)
  for(i in 1:K){

    RF_old <- computeRandomForest_PC(data=data,
                                     outcome = outcome,
                                     unit=unit, repeatedId=repeatedId,
                                     varImpPlot = F,
                                     metaNames=c(metaNames),
                                     nTrees=nTrees)

    data_merge <- RF_old$varImportanceData[,c('var','avgVI')]
    colnames(data_merge) <- c('var',paste0('avgVIK',i))
    avgVI_old <- merge(avgVI_old, data_merge, by='var')
  }
  data_summ <- .summData(avgVI_old)[,c('var','est')]
  tmp <- .summData(avgVI)
  data_summ <- merge(data_summ,
                     tmp[,c('var','sd')])


  ## Permutation
  avgVI_perm <- data.frame('var'=underlyingVars)

  if(!silent) cat('Permutation Trials (',synthetics,'): ',sep='')
  for(sim in 1:synthetics){
    if(!silent) cat(sim,', ',sep='')

    # Permute but do the functional components together
    data_permute <- data[colnames(data) %in% c(outcome,unit,repeatedId)]
    data_permute <- cbind(data_permute,
                   as.data.frame(sapply(1:length(underlyingVars),
                                        function(idx, DF) {
                                          DF[sample.int(nrow(DF)),
                                             underlyingDataAlignedFunctions==underlyingVars[idx],
                                             drop=FALSE]
                                        }, simplify = 'array',
                                        DF=data)))
    ## Ensure 3 are permuted together
    #   Set PC to 1

    # Get RF and VI
    #numDrop <- as.integer(nrow(data)/K) +
    #              rbinom(1,1,nrow(data_permute)/K-as.integer(nrow(data)/K))
    #drop3Idx <- sample(1:nrow(data_permute),numDrop)
    RF <- computeRandomForest_PC(data=data_permute,#[-drop3Idx,],
                                 outcome = outcome,
                                 unit=unit, repeatedId=repeatedId,
                                 varImpPlot = F,
                                 metaNames=metaNames,
                                 nTrees=nTrees,keepModels=F)

    data_merge <- RF$varImportanceData[,c('var','avgVI')]
    colnames(data_merge) <- c('var',paste0('avgVIK',sim))
    avgVI_perm <- merge(avgVI_perm, data_merge, by='var')
  }
  if(!silent) cat('\n')

  # Check rescaling
  #   No meta-vars and no rescaling

  # Summarize Data
  #data_summ <- .summData(avgVI)
  data_perm_summ <- .summData(avgVI_perm)
  tmp <- as.data.frame(t(sapply(X = data_perm_summ$var,
                  FUN = function(var, alpha, data_vi){
                    c(var,
                      unlist(stats::quantile(as.numeric(data_vi[data_vi$var==var,-1]),
                                           probs=1-alpha)))
                  }, alpha=alpha, data_vi=avgVI_perm, USE.NAMES = FALSE)))
  colnames(tmp) <- c('var','upAlpha')
  tmp$upAlpha <- as.numeric(tmp$upAlpha)
  data_perm_summ <- merge(data_perm_summ,tmp)

  ## Standardize Data
  # Get Noise cutoff
  cutoffs_noise <- data.frame('Type'=c('KFunction',metaNames),
                        'Cutoff'=c(max(data_perm_summ[data_perm_summ$var %in% KFunctions,'upAlpha']),
                                   rep(NA,length(metaNames))),
                        'Adj'=NA)

  for(meta in metaNames){
    cutoffs_noise[cutoffs_noise$Type==meta,'Cutoff'] <-
            data_perm_summ[data_perm_summ$var==meta,'upAlpha']
  }
  cutoffs_noise$Adj <- max(cutoffs_noise$Cutoff)/cutoffs_noise$Cutoff

  noiseCO <- max(cutoffs_noise$Cutoff)

  data_perm_summ[data_perm_summ$var%in%KFunctions,'AdjAmt'] <-
    cutoffs_noise[cutoffs_noise$Type=='KFunction','Adj']
  if(length(metaNames)>0){
    for(var in metaNames){
      data_perm_summ[data_perm_summ$var==var,'AdjAmt'] <-
        cutoffs_noise[cutoffs_noise$Type==var,'Adj']
    }
  }

  data_perm_summ$adj <- data_perm_summ$est*data_perm_summ$AdjAmt
  data_perm_summ$adj <- ifelse(is.infinite(data_perm_summ$adj),0,data_perm_summ$adj)

  data_summ_std <- data_summ
  data_summ_std$est <- data_summ_std$est * data_perm_summ$AdjAmt
  data_summ_std$sd <- data_summ_std$sd * data_perm_summ$AdjAmt

  avgVI_perm_std <- avgVI_perm
  avgVI_perm_std[-1] <- avgVI_perm_std[-1] * data_perm_summ$AdjAmt

  # Get Interpolation cutoff
  interpolationCO <- apply(apply(avgVI_perm_std[-1], MARGIN=2, FUN=sort, decreasing = T),
                           MARGIN = 1, FUN = quantile, probs = 1-alpha)

  # Model accuracy estimates
  data_modelAcc <- .computeModelAccuracy(oobAcc=oobAcc, outcomes=data[,outcome])


  ## Get plots and organize results
  append(list('VariableImportance'=data_summ_std,
              'AccuracyEstimate'=data_modelAcc,
              'NoiseCutoff'=noiseCO,
              'InterpolationCutoff'=interpolationCO,
              'AdditionalParams'=list('alpha'=alpha,
                                      'subsetPlotSize'=subsetPlotSize)),
        .generateVIPlot(
          viData=data_summ_std,
          accData=data_modelAcc,
          NoiseCutoff=noiseCO,
          InterpolationCutoff=interpolationCO,
          subsetPlotSize=subsetPlotSize)
  )
}

#' Title
#'
#' @param dat
#'
#' @return
#'
#' @examples
.summData <- function(dat){
  data.frame('var'=dat$var,
             'est'=rowMeans(dat[-1]),
             'sd'=apply(dat[-1],MARGIN=1,FUN=function(x){sd(x)}))
}

#' Title
#'
#' @param oobAcc
#' @param outcomes
#' @param rGuessSims
#'
#' @return
#'
#' @examples
.computeModelAccuracy <- function(oobAcc, outcomes, rGuessSims=500){

  optVals <- as.numeric(table(outcomes))
  n <- length(outcomes)

  # OOB
  if(!is.null(oobAcc)){
    accData <- data.frame('OOB'=mean(oobAcc))
  } else{
    accData <- NULL
  }

  if(!is.null(outcomes)){
    # Bias guess - Pick most popular group
    accData$bias <- max(optVals)/sum(optVals)

    # Random guessing based on total sample
    acc <- rep(NA, rGuessSims)
    for(i in 1:rGuessSims){
      # Columns relate to outcomes
      guesses <- rmultinom(length(outcomes),size=1,optVals/sum(optVals))
      acc[i] <- sum(which(guesses==1,arr.ind = T)[,1] ==
                      as.integer(as.factor(outcomes)))/n
    }
    accData$guess <- mean(acc)
  }

  accData
}

#' Title
#'
#' @param viData
#' @param accData
#' @param NoiseCutoff
#' @param InterpolationCutoff
#' @param subsetPlotSize
#'
#' @return
#'
#' @examples
.generateVIPlot <- function(viData, accData, NoiseCutoff,
                            InterpolationCutoff, subsetPlotSize){
  viPlot <- .plotVI(viData, accData, NoiseCutoff, InterpolationCutoff)

  data_return <- list('viPlot'=viPlot)

  if(nrow(viData)>subsetPlotSize){
    # Order Data to take top
    tmpStdVI <- viData[order(-viData$est),]

    viPlot <- .plotVI(tmpStdVI[1:subsetPlotSize,], accData, NoiseCutoff,
                      InterpolationCutoff[1:subsetPlotSize])

    data_return <- append(
      data_return,
      list('subset_viPlot'=viPlot)
    )
  }

  data_return
}


#' Title
#'
#' @param viData
#' @param accData
#' @param NoiseCutoff
#' @param InterpolationCutoff
#'
#' @return
#'
#' @examples
.plotVI <- function(viData, accData, NoiseCutoff,
                    InterpolationCutoff){
  maxVal <- max(InterpolationCutoff,NoiseCutoff,viData$est)
  returnPlot <- ggplot2::ggplot(data=viData,
                                mapping=ggplot2::aes(x=factor(reorder(var,est)),
                                                     y=ifelse(est/maxVal>1,1,
                                                              ifelse(est/maxVal<0,0,
                                                                     est/maxVal)),
                                                     group=1)) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin=ifelse((est-sd)/maxVal<0,0,(est-sd)/maxVal),
                                        ymax = ifelse((est+sd)/maxVal>1,1,(est+sd)/maxVal)),
                           color='black', width=0.2) +
    ggplot2::geom_point(color='black') +
    ggplot2::geom_hline(ggplot2::aes(yintercept=max(0,min(1,NoiseCutoff/maxVal))),
                        color='red', linetype='dotted',linewidth=1) +
    ggplot2::coord_flip(ylim = c(0,1)) +
    ggplot2::xlab(NULL) +
    ggplot2::ylim(c(0,1)) +
    ggplot2::ylab(NULL) +
    ggplot2::theme_bw() +
    ggplot2::geom_line(ggplot2::aes(
      x=ordered(viData[order(-est),'var']),
      y=InterpolationCutoff/maxVal),color='orange', linetype='dashed',size=1) +
    ggplot2::ylab(paste0('Variable Importance - OOB (',
                         .specify_decimal(accData$OOB,2),
                         '), Guess (',
                         .specify_decimal(accData$guess,2),
                         '), Bias (',
                         .specify_decimal(accData$bias,2),')'))
}
