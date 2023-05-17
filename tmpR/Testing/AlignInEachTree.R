# FK4
#   CV intervals based on only data
#   OOB based on only data
#   Point estimates (real data + 100vars per var)
#   Curve is permuted data + 100synt per tree to std

funkyModel4 <- function(data, K = 10,
                        outcome = colnames(data)[1],
                        unit = colnames(data)[2],
                        metaNames = NULL,
                        synthetics = 100,
                        alpha = 0.05,
                        silent = FALSE,
                        rGuessSims = 500,
                        subsetPlotSize = 25, nTrees = 500,
                        method = "class") {
  ## From moving code over, will remove
  repeatedId <- NULL

  ## Error checking
  # .checkData(alignmentMethod) ## TODO:: Add something in

  ## Generate Synthetics And Connect
  components <- colnames(data)[!(colnames(data) %in%
                                   c(outcome, unit, repeatedId, metaNames))]
  KFunctions <- .getUnderlyingVariable(components)
  # Get Var and data column alignment
  underlyingDataAlignedFunctions <- .getUnderlyingVariable(colnames(data),
                                                           returnUnique = FALSE
  )
  underlyingVars <- unique(underlyingDataAlignedFunctions)
  underlyingVars <- underlyingVars[!(underlyingVars %in%
                                       c(outcome, unit, repeatedId))]

  avgVI <- avgVI_full <-
    data.frame("var" = c(underlyingVars,
                         paste0('permuteInternal',1:synthetics,'K_K'),
                         as.vector(sapply(metaNames,function(m){
                           paste0(paste0('permuteInternal',1:synthetics),m)}))))
  oobAcc <- rep(NA, K)
  groups <- .getFolds(1:nrow(data), K)

  # K-Fold Cross-Validation (True Data)
  if (!silent) cat("CV Trials (", K, "): ", sep = "")
  for (i in 1:K) {
    if (!silent) cat(i, ", ", sep = "")

    # Do CV without noise for Sd and OOB
    RF <- funkyForest(
      data = data[-groups[[i]], ],
      outcome = outcome,
      unit = unit, repeatedId = repeatedId,
      varImpPlot = FALSE,
      metaNames = c(metaNames),
      nTrees = nTrees, method = method
    )

    data_merge <- RF$varImportanceData[, c("var", "avgVI")]
    colnames(data_merge) <- c("var", paste0("avgVIK", i))
    avgVI <- merge(avgVI, data_merge, by = "var")

    oobAcc[i] <- sum(data[groups[[i]], outcome] ==
                       predict_funkyForest(
                         model = RF$model,
                         data_pred = data[groups[[i]], ],
                         type = "pred", data = data
                       )) /
      nrow(data[groups[[i]], ])

    ## Run on all for VI estimate

    # Permute Noise but do the functional components together
    data_full <- data
    synths <- c('K',metaNames)
    for(j in 1:synthetics) {
      ## TODO FIX
      data_permute <-
        as.data.frame(sapply(1:length(synths),
                             function(idx, DF) {
                               if(idx==1){
                                 matchedVal <- sample(KFunctions,1)
                               }else{
                                 matchedVal <- synths[idx]
                               }
                               DF[sample.int(nrow(DF)),
                                  underlyingDataAlignedFunctions == matchedVal,
                                  drop = FALSE
                               ]
                             },
                             simplify = "array",
                             DF = data
        ))
      colnames(data_permute)[!(colnames(data_permute) %in% metaNames)] <-
        paste0('permuteInternal',j,'K_K_PC',1:3)
        colnames(data_permute)[colnames(data_permute) %in% metaNames] <-
        paste0('permuteInternal',j,metaNames)
      data_full <- cbind(data_full, data_permute)
    }

    RF_full <- funkyForest(
      data = data_full,
      outcome = outcome,
      unit = unit, repeatedId = repeatedId,
      varImpPlot = FALSE,
      metaNames = c(metaNames,
                    as.vector(sapply(metaNames,function(m){
                      paste0(paste0('permuteInternal',1:synthetics),m)}))),
      nTrees = nTrees,keepModels = FALSE
    )

    data_merge <- RF_full$varImportanceData[, c("var", "avgVI")]
    colnames(data_merge) <- c("var", paste0("avgVIK", i))
    avgVI_full <- merge(avgVI_full, data_merge, by = "var")
  }
  if (!silent) cat("\n")

  ## Permutation
  avgVI_perm <-
    data.frame("var" = c(underlyingVars,
                         paste0('permuteInternal',1:synthetics,'K_K'),
                         as.vector(sapply(metaNames,function(m){
                           paste0(paste0('permuteInternal',1:synthetics),m)}))))

  if (!silent) cat("Permutation Trials (", synthetics, "): ", sep = "")
  for (sim in 1:synthetics) {
    if (!silent) cat(sim, ", ", sep = "")

    # Permute but do the functional components together
    data_permute1 <-
      as.data.frame(sapply(1:length(underlyingVars),
                           function(idx, DF) {
                             DF[sample.int(nrow(DF)),
                                underlyingDataAlignedFunctions == underlyingVars[idx],
                                drop = FALSE
                             ]
                           },
                           simplify = "array",
                           DF = data
      ))
    data_permute2 <- data[!(colnames(data) %in%
                         c(outcome, unit, repeatedId))]
    synths <- c('K',metaNames)
    for(j in 1:synthetics) {
      data_tmp <-
        as.data.frame(sapply(1:length(synths),
                             function(idx, DF) {
                               if(idx==1){
                                 matchedVal <- sample(KFunctions,1)
                               }else{
                                 matchedVal <- synths[idx]
                               }
                               DF[sample.int(nrow(DF)),
                                  underlyingDataAlignedFunctions == matchedVal,
                                  drop = FALSE
                               ]
                             },
                             simplify = "array",
                             DF = data
        ))
      colnames(data_tmp)[!(colnames(data_tmp) %in% metaNames)] <-
        paste0('permuteInternal',j,'K_K_PC',1:3)
      colnames(data_tmp)[colnames(data_tmp) %in% metaNames] <-
        paste0('permuteInternal',j,metaNames)
      data_permute2 <- cbind(data_permute2, data_tmp)
    }

    data_permute <- cbind(data[colnames(data) %in%
                                 c(outcome, unit, repeatedId)],
                          data_permute1,
                          data_permute2)

    # Get RF and VI
    RF <- funkyForest(
      data = data_permute,
      outcome = outcome,
      unit = unit, repeatedId = repeatedId,
      varImpPlot = FALSE,
      metaNames = c(metaNames,
                    as.vector(sapply(metaNames,function(m){
                      paste0(paste0('permuteInternal',1:synthetics),m)}))),
      nTrees = nTrees, keepModels = FALSE
    )

    data_merge <- RF$varImportanceData[, c("var", "avgVI")]
    colnames(data_merge) <- c("var", paste0("avgVIK", sim))
    avgVI_perm <- merge(avgVI_perm, data_merge, by = "var")
  }
  if (!silent) cat("\n")


  # Summarize Data
  #data_summ <- .summData(avgVI_full)[, c("var", "est","sd")]

  tmp <- as.data.frame(t(sapply(
    X = underlyingVars,
    FUN = function(var1, alpha, data_vi, KFunctions) {
      if(var1 %in% KFunctions){
        tmp  <- data.frame(
          var1,
          t(as.numeric(apply(X = data_vi[data_vi$var %in% paste0('permuteInternal',1:100,'K_K'), -1],
                             MARGIN = 2,FUN = quantile, probs=1-alpha)))
          )
      }else{
        tmp  <- data.frame(
          var1,
          t(as.numeric(apply(X = data_vi[data_vi$var %in% paste0('permuteInternal',1:100,var1), -1],
                MARGIN = 2,FUN = quantile, probs=1-alpha)))
        )
      }
      tmp
    },
    alpha = alpha, data_vi = avgVI_full, KFunctions=KFunctions,
    USE.NAMES = FALSE,simplify = 'matrix'
  )))
  for(i in 2:ncol(tmp)){
    tmp[i]<-unlist(tmp[i])
  }

  noiseCO <- max(tmp[-1])

  # If value is 0, will jsut take rather than scale
  tmp1 <- noiseCO/tmp[,-1]
  tmp1[sapply(tmp1, simplify = 'matrix', is.infinite)] <- 1

  avgVI_std <- cbind('var'=avgVI_full[avgVI_full$var %in% underlyingVars,1],
                     avgVI_full[avgVI_full$var %in% underlyingVars,-1] * tmp1)
  data_summ <- .summData(avgVI_std)

  # Add CV SDs (which do not include any noise Vars)
  tmp <- .summData(avgVI)
  colnames(tmp)[3]<- 'cvSD'
  data_summ <- merge(
    data_summ,
    tmp[, c("var", "cvSD")],
    all.x = TRUE
  )

  ## Standardize Data
  # Get Noise cutoff

  tmp <- as.data.frame(t(sapply(
    X = underlyingVars,
    FUN = function(var, alpha, data_vi, KFunctions) {
      if(var %in% KFunctions){
        tmp  <- data.frame(
          var,
          t(as.numeric(apply(X = data_vi[data_vi$var %in% paste0('permuteInternal',1:100,'K_K'), -1],
                             MARGIN = 2,FUN = quantile, probs=1-alpha)))
        )
      }else{
        tmp  <- data.frame(
          var,
          t(as.numeric(apply(X = data_vi[data_vi$var %in% paste0('permuteInternal',1:100,var), -1],
                             MARGIN = 2,FUN = quantile, probs=1-alpha)))
        )
      }
      tmp
    },
    alpha = alpha, data_vi = avgVI_perm, KFunctions=KFunctions,
    USE.NAMES = FALSE,simplify = 'matrix'
  )))
  for(i in 2:ncol(tmp)){
    tmp[i]<-unlist(tmp[i])
  }

  # If value is 0, will jsut take rather than scale
  tmp1 <- noiseCO/tmp[,-1]
  tmp1[sapply(tmp1, simplify = 'matrix', is.infinite)] <- 1

  avgVI_perm_std <- cbind('var'=avgVI_perm[avgVI_perm$var %in% underlyingVars,1],
                          avgVI_perm[avgVI_perm$var %in% underlyingVars,-1] * tmp1)
  data_perm_summ <- .summData(avgVI_perm_std)


  # Get Interpolation cutoff
  interpolationCO <- apply(
    apply(avgVI_perm_std[-1],
          MARGIN = 2, FUN = sort,
          decreasing = TRUE
    ),
    MARGIN = 1, FUN = stats::quantile, probs = 1 - alpha
  )

  # Model accuracy estimates
  data_modelAcc <- .computeModelAccuracy(oobAcc = oobAcc,
                                         outcomes = data[, outcome])


  ## Get plots and organize results
  append(
    list(
      "VariableImportance" = data_summ,
      "AccuracyEstimate" = data_modelAcc,
      "NoiseCutoff" = noiseCO,
      "InterpolationCutoff" = interpolationCO,
      "AdditionalParams" = list(
        "alpha" = alpha,
        "subsetPlotSize" = subsetPlotSize
      )
    ),
    .generateVIPlot(
      viData = data_summ,
      accData = data_modelAcc,
      NoiseCutoff = noiseCO,
      InterpolationCutoff = interpolationCO,
      subsetPlotSize = subsetPlotSize
    )
  )
}
