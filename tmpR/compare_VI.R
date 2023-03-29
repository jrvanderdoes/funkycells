computeRandomForest_CVPC_Permute <- function(data, K=10,
                                             outcome=colnames(data)[1],
                                             unit=colnames(data)[2],
                                             repeatedId=NULL,
                                             metaNames=NULL,
                                             cellData=NULL,
                                             synthetics=100,
                                             alpha=0.05,
                                             silent=F,
                                             rGuessSims=500,
                                             subsetPlotSize=25, nTrees=500){
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

  avgVI <- avgVI1 <- data.frame('var'=underlyingVars)
  oobAcc <- rep(NA,K)
  groups <- .getFolds(1:nrow(data), K)

  # K-Fold Cross-Validation (True Data)
  data1 <- data
  data1[[outcome]] <- as.numeric(data1[[outcome]])
  if(!silent) cat('CV Trials (',K,'): ',sep='')
  for(i in 1:K){
    if(!silent) cat(i,', ',sep='')
    RF <- computeRandomForest_PC(data=data[-groups[[i]],],
                                 outcome = outcome,
                                 unit=unit, repeatedId=repeatedId,
                                 varImpPlot = F,
                                 metaNames=c(metaNames),
                                 nTrees=nTrees)
    RF1 <- party::cforest(Stage ~ ., data=data1[-groups[[i]],!(colnames(data1)%in%c(unit,repeatedId))])

    data_merge <- RF$varImportanceData[,c('var','avgVI')]
    colnames(data_merge) <- c('var',paste0('avgVIK',i))
    avgVI <- merge(avgVI, data_merge, by='var')

    tmp <- as.data.frame(varimp(RF1))
    tmp1 <- data.frame('var'=rownames(tmp),
                             'avgVI'=tmp[,1])
    tmp1$underlying <- .getUnderlyingVariable(tmp1$var,returnUnique = F)
    data_merge <- data.frame('var'=unique(tmp1$underlying),
                             'avgVI'=NA)
    for(j in 1:nrow(data_merge)){
      data_merge[j,'avgVI'] <- sum(tmp1[tmp1$underlying==data_merge$var[j],'avgVI'])
    }

    colnames(data_merge) <- c('var',paste0('avgVIK',i))
    avgVI1 <- merge(avgVI1, data_merge, by='var')

    oobAcc[i] <- sum(data[groups[[i]],outcome]==
                       predict.RandomForest_PC(model = RF$model,
                                               data_pred = data[groups[[i]],],
                                               type = 'pred', data = data)) /
      nrow(data[groups[[i]],])
  }
  if(!silent) cat('\n')


  ## Permutation
  avgVI_perm <- avgVI_perm1 <- data.frame('var'=underlyingVars)

  if(!silent) cat('Permutation Trials (',synthetics,'): ',sep='')
  for(sim in 1:synthetics){
    if(!silent) cat(sim,', ',sep='')

    # Permute but do the functional components together
    data_permute <- data[,colnames(data) %in% c(outcome,unit,repeatedId)]
    data_permute <- cbind(data_permute,
                          as.data.frame(sapply(1:length(underlyingVars),
                                               function(idx, DF) {
                                                 DF[sample.int(nrow(DF)),
                                                    underlyingDataAlignedFunctions==underlyingVars[idx],
                                                    drop=FALSE]
                                               }, simplify = 'array',
                                               DF=data)))

    # Get RF and VI
    RF <- computeRandomForest_PC(data=data_permute,
                                 outcome = outcome,
                                 unit=unit, repeatedId=repeatedId,
                                 varImpPlot = F,
                                 metaNames=metaNames,
                                 nTrees=nTrees,keepModels=F)


    data_merge <- RF$varImportanceData[,c('var','avgVI')]
    colnames(data_merge) <- c('var',paste0('avgVIK',sim))
    avgVI_perm <- merge(avgVI_perm, data_merge, by='var')

    # Get RF and VI
    data_permute1 <- data_permute
    data_permute1[[outcome]] <- as.numeric(data_permute1[[outcome]])
    RF1 <- party::cforest(Stage ~ ., data=data_permute1[,!(colnames(data_permute1)%in%c(unit,repeatedId))])

    tmp <- as.data.frame(varimp(RF1))
    tmp1 <- data.frame('var'=rownames(tmp),
                       'avgVI'=tmp[,1])
    tmp1$underlying <- .getUnderlyingVariable(tmp1$var,returnUnique = F)
    data_merge <- data.frame('var'=unique(tmp1$underlying),
                             'avgVI'=NA)
    for(j in 1:nrow(data_merge)){
      data_merge[j,'avgVI'] <- sum(tmp1[tmp1$underlying==data_merge$var[j],'avgVI'])
    }

    colnames(data_merge) <- c('var',paste0('avgVIK',sim))
    avgVI_perm1 <- merge(avgVI_perm1, data_merge, by='var')
  }
  if(!silent) cat('\n')

  # Summarize Data
  data_summ <- .summData(avgVI)
  data_summ1 <- .summData(avgVI1)
  data_perm_summ <- .summData(avgVI_perm)
  data_perm_summ1 <- .summData(avgVI_perm1)
  tmp <- as.data.frame(t(sapply(X = data_perm_summ$var,
                                FUN = function(var, alpha, data_vi){
                                  c(var,
                                    unlist(stats::quantile(as.numeric(data_vi[data_vi$var==var,-1]),
                                                           probs=1-alpha)))
                                }, alpha=alpha, data_vi=avgVI_perm, USE.NAMES = FALSE)))
  colnames(tmp) <- c('var','upAlpha')
  tmp$upAlpha <- as.numeric(tmp$upAlpha)
  data_perm_summ <- merge(data_perm_summ,tmp)

  tmp <- as.data.frame(t(sapply(X = data_perm_summ1$var,
                                FUN = function(var, alpha, data_vi){
                                  c(var,
                                    unlist(stats::quantile(as.numeric(data_vi[data_vi$var==var,-1]),
                                                           probs=1-alpha)))
                                }, alpha=alpha, data_vi=avgVI_perm1, USE.NAMES = FALSE)))
  colnames(tmp) <- c('var','upAlpha')
  tmp$upAlpha <- as.numeric(tmp$upAlpha)
  data_perm_summ1 <- merge(data_perm_summ1,tmp)

  # Get Noise cutoff
  noiseCO <- max(data_perm_summ$upAlpha)
  noiseCO1 <- max(data_perm_summ1$upAlpha)

  # Standardize Data
  data_perm_summ$adj <- noiseCO/data_perm_summ$upAlpha
  data_perm_summ$adj <- ifelse(is.infinite(data_perm_summ$adj),0,data_perm_summ$adj)
  data_perm_summ1$adj <- noiseCO1/data_perm_summ1$upAlpha
  data_perm_summ1$adj <- ifelse(is.infinite(data_perm_summ1$adj),0,data_perm_summ1$adj)

  data_summ_std <- data_summ
  data_summ_std$est <- data_summ_std$est * data_perm_summ$adj
  data_summ_std$sd <- data_summ_std$sd * data_perm_summ$adj
  data_summ_std1 <- data_summ1
  data_summ_std1$est <- data_summ_std1$est * data_perm_summ1$adj
  data_summ_std1$sd <- data_summ_std1$sd * data_perm_summ1$adj

  avgVI_perm_std <- avgVI_perm
  avgVI_perm_std[-1] <- avgVI_perm_std[-1] * data_perm_summ$adj
  avgVI_perm_std1 <- avgVI_perm1
  avgVI_perm_std1[-1] <- avgVI_perm_std1[-1] * data_perm_summ1$adj

  # Get Interpolation cutoff
  interpolationCO <- apply(apply(avgVI_perm_std[-1], MARGIN=2, FUN=sort, decreasing = T),
                           MARGIN = 1, FUN = quantile, probs = 1-alpha)
  interpolationCO1 <- apply(apply(avgVI_perm_std1[-1], MARGIN=2, FUN=sort, decreasing = T),
                           MARGIN = 1, FUN = quantile, probs = 1-alpha)

  # Model accuracy estimates
  data_modelAcc <- .computeModelAccuracy(oobAcc=oobAcc, outcomes=data[,outcome])


  ## Get plots and organize results
  base <- append(list('VariableImportance'=data_summ_std,
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
  base1 <- append(list('VariableImportance'=data_summ_std1,
                      'AccuracyEstimate'=data_modelAcc,
                      'NoiseCutoff'=noiseCO1,
                      'InterpolationCutoff'=interpolationCO1,
                      'AdditionalParams'=list('alpha'=alpha,
                                              'subsetPlotSize'=subsetPlotSize)),
                 .generateVIPlot(
                   viData=data_summ_std1,
                   accData=data_modelAcc,
                   NoiseCutoff=noiseCO1,
                   InterpolationCutoff=interpolationCO1,
                   subsetPlotSize=subsetPlotSize)
  )
}
