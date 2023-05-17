# FK3
#   CV intervals based on only data
#   OOB based on only data
#   Point estimates (real data + 1 permute per var)
#   Curve is permuted data + 1 permute to std
#   Diff std method


funkyModel3 <- function(data, K = 10,
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

  avgVI <- avgVI_full <- data.frame("var" = c(underlyingVars,
                                              paste0('permuteInternal',
                                                     underlyingVars)))
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
    data_permute <-
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
    colnames(data_permute) <- paste0('permuteInternal',colnames(data_permute))
    data_full <- cbind(data, data_permute)

    RF_full <- funkyForest(
      data = data_full,
      outcome = outcome,
      unit = unit, repeatedId = repeatedId,
      varImpPlot = FALSE,
      metaNames = c(metaNames,paste0('permuteInternal',metaNames)),
      nTrees = nTrees
    )

    data_merge <- RF_full$varImportanceData[, c("var", "avgVI")]
    colnames(data_merge) <- c("var", paste0("avgVIK", i))
    avgVI_full <- merge(avgVI_full, data_merge, by = "var")
  }
  if (!silent) cat("\n")

  ## Permutation
  avgVI_perm <- data.frame("var" = c(underlyingVars,
                                     paste0('permuteInternal',underlyingVars)))

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
    data_permute2 <-
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
    colnames(data_permute2) <- paste0('permuteInternal',colnames(data_permute2))
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
      metaNames = c(metaNames,paste0('permuteInternal',metaNames)),
      nTrees = nTrees, keepModels = FALSE
    )

    data_merge <- RF$varImportanceData[, c("var", "avgVI")]
    colnames(data_merge) <- c("var", paste0("avgVIK", sim))
    avgVI_perm <- merge(avgVI_perm, data_merge, by = "var")
  }
  if (!silent) cat("\n")


  # Summarize Data
  data_summ <- .summData(avgVI_full)[, c("var", "est","sd")]

  # Add upAlpha for upper \alpha percentile to standardize
  tmp <- as.data.frame(t(sapply(
    X = paste0('permuteInternal',underlyingVars),#data_summ$var,
    FUN = function(var, alpha, data_vi) {
      c(
        var,
        unlist(stats::quantile(as.numeric(data_vi[data_vi$var == var, -1]),
                               probs = 1 - alpha
        ))
      )
    }, alpha = alpha, data_vi = avgVI_full, USE.NAMES = FALSE
  )))
  colnames(tmp) <- c("var", "upAlpha")
  tmp$upAlpha <- as.numeric(tmp$upAlpha)
  data_summ <- merge(data_summ, tmp, all.x = TRUE)

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
  noiseCO <- max(data_summ[data_summ$var %in%
                             paste0('permuteInternal',underlyingVars),
                           "upAlpha"])
  # noiseCO <- max(data_perm_summ[data_perm_summ$var %in%
  #                                 paste0('permuteInternal',underlyingVars),
  #                               "upAlpha"])

  data_summ[["adjAmt"]] <- NA
  data_summ[data_summ$var %in% paste0('permuteInternal',underlyingVars), "adjAmt"] <-
    noiseCO / data_summ[data_summ$var %in% paste0('permuteInternal',underlyingVars),
                        "upAlpha"]
  data_summ$adjAmt <-
    ifelse(is.infinite(data_summ$adjAmt), 1, data_summ$adjAmt)

  tmp <- c(NA)
  tmp_vals <- c(NA)
  for(uVar in paste0('permuteInternal',underlyingVars)){
    if(!(uVar %in% metaNames)){
      tmp <- c(tmp,uVar)
      tmp_vals <- c(tmp_vals,data_summ[data_summ$var==uVar,"upAlpha"])
    }
  }
  data_summ[data_summ$var %in% tmp,"upAlpha"] <- mean(tmp_vals)

  data_summ[['stdEst']] <- data_summ[['stdSD']] <- NA
  # data_perm_summ[['stdEst']] <- data_perm_summ[['stdSD']] <- NA
  for(uVar in underlyingVars){
    data_summ[data_summ$var==uVar,c('stdEst','stdSD')] <-
      data_summ[data_summ$var==uVar,c('est','cvSD')] *
      data_summ[data_summ$var==paste0('permuteInternal',uVar),'adjAmt']
    # data_perm_summ[data_perm_summ$var==uVar,c('stdEst','stdSD')] <-
    #   data_perm_summ[data_perm_summ$var==uVar,c('est','sd')] *
    #   data_perm_summ[data_perm_summ$var==paste0('permuteInternal',uVar),'adjAmt']
  }

  ## For use in current code
  data_summ_std <- data.frame('var'=data_summ[data_summ$var %in% underlyingVars,'var'],
                              'est'=data_summ[data_summ$var %in% underlyingVars,'stdEst'],
                              'sd'=data_summ[data_summ$var %in% underlyingVars,'stdSD'])


  # data_summ[,c('est','sd')] <- data_summ[,c('est','sd')] * data_summ$adjAmt
  # data_perm_summ[,c('est','sd')] <- data_perm_summ[,c('est','sd')] * data_perm_summ$adjAmt



  # data_perm_summ$adjEst <- data_perm_summ$est * data_perm_summ$adjAmt
  # data_perm_summ$adjEst <-
  #   ifelse(is.infinite(data_perm_summ$adjEst), 0, data_perm_summ$adjEst)
  # data_perm_summ$adjSd <- data_perm_summ$sd * data_perm_summ$adjAmt

  ### CURVE DATA

  data_perm_summ <- .summData(avgVI_perm)
  # Add upAlpha to std.
  tmp <- as.data.frame(t(sapply(
    X = paste0('permuteInternal',underlyingVars),#data_perm_summ$var,
    FUN = function(var, alpha, data_vi) {
      c(
        var,
        unlist(stats::quantile(as.numeric(data_vi[data_vi$var == var, -1]),
                               probs = 1 - alpha
        ))
      )
    }, alpha = alpha, data_vi = avgVI_perm, USE.NAMES = FALSE
  )))
  colnames(tmp) <- c("var", "upAlpha")
  tmp$upAlpha <- as.numeric(tmp$upAlpha)
  data_perm_summ <- merge(data_perm_summ, tmp, all.x = TRUE)

  noiseCO_p <- max(data_perm_summ[data_perm_summ$var %in%
                                    paste0('permuteInternal',underlyingVars),
                                  "upAlpha"])
  data_perm_summ[["adjAmt"]] <- NA
  data_perm_summ[data_perm_summ$var %in% paste0('permuteInternal',underlyingVars), "adjAmt"] <-
    noiseCO_p / data_perm_summ[data_perm_summ$var %in% paste0('permuteInternal',underlyingVars),
                               "upAlpha"]
  data_perm_summ$adjAmt <-
    ifelse(is.infinite(data_perm_summ$adjAmt), 1, data_perm_summ$adjAmt)

  for(uVar in paste0('permuteInternal',underlyingVars)){
    if(!(uVar %in% metaNames)){
      tmp <- c(tmp,uVar)
      tmp_vals <- c(tmp_vals,data_perm_summ[data_perm_summ$var==uVar,"upAlpha"])
    }
  }
  data_perm_summ[data_perm_summ$var %in% tmp,"upAlpha"] <- mean(tmp_vals)

  avgVI_perm_std <- avgVI_perm
  for(uVar in underlyingVars){
    avgVI_perm_std[avgVI_perm_std$var==uVar,-1] <-
      avgVI_perm_std[avgVI_perm_std$var==uVar,-1] *
      data_perm_summ[data_perm_summ$var==paste0('permuteInternal',uVar),'adjAmt']
    # data_perm_summ[data_perm_summ$var==uVar,c('stdEst','stdSD')] <-
    #   data_perm_summ[data_perm_summ$var==uVar,c('est','sd')] *
    #   data_perm_summ[data_perm_summ$var==paste0('permuteInternal',uVar),'adjAmt']
  }

  # Get Interpolation cutoff
  interpolationCO <- apply(
    #apply(avgVI_perm_std[avgVI_perm_std$var %in% underlyingVars,-1],
    apply(avgVI_perm[avgVI_perm$var %in% underlyingVars,-1],
          MARGIN = 2, FUN = sort,
          decreasing = TRUE
    ),
    MARGIN = 1, FUN = stats::quantile, probs = 1 - alpha
  )

  avgVI_perm_std <- avgVI_perm
  for(uVar in underlyingVars){
    avgVI_perm_std[avgVI_perm_std$var==uVar,-1] <-
      avgVI_perm_std[avgVI_perm_std$var==uVar,-1] *
      data_perm_summ[data_perm_summ$var==paste0('permuteInternal',uVar),'adjAmt']
    # data_perm_summ[data_perm_summ$var==uVar,c('stdEst','stdSD')] <-
    #   data_perm_summ[data_perm_summ$var==uVar,c('est','sd')] *
    #   data_perm_summ[data_perm_summ$var==paste0('permuteInternal',uVar),'adjAmt']
  }

  # avgVI_perm_std <- avgVI_perm
  # avgVI_perm_std[-1] <- avgVI_perm_std[-1] * data_perm_summ$adjAmt

  # cutoffs_noise <- data.frame(
  #   "Type" = c("KFunction", metaNames),
  #   "Cutoff" = c(
  #     max(data_perm_summ[data_perm_summ$var %in% KFunctions, "upAlpha"]),
  #     rep(NA, length(metaNames))
  #   ),
  #   "Adj" = NA
  # )

  # for (meta in metaNames) {
  #   cutoffs_noise[cutoffs_noise$Type == meta, "Cutoff"] <-
  #     data_perm_summ[data_perm_summ$var == meta, "upAlpha"]
  # }
  # cutoffs_noise$Adj <- max(cutoffs_noise$Cutoff) / cutoffs_noise$Cutoff

  # noiseCO <- max(cutoffs_noise$Cutoff)

  # data_perm_summ[data_perm_summ$var %in% KFunctions, "AdjAmt"] <-
  #   cutoffs_noise[cutoffs_noise$Type == "KFunction", "Adj"]
  # if (length(metaNames) > 0) {
  #   for (var in metaNames) {
  #     data_perm_summ[data_perm_summ$var == var, "AdjAmt"] <-
  #       cutoffs_noise[cutoffs_noise$Type == var, "Adj"]
  #   }
  # }
  #
  # data_perm_summ$adj <- data_perm_summ$est * data_perm_summ$AdjAmt
  # data_perm_summ$adj <-
  #   ifelse(is.infinite(data_perm_summ$adj), 0, data_perm_summ$adj)

  # data_summ_std <- data_summ
  # data_summ_std$est <- data_summ_std$est * data_perm_summ$AdjAmt
  # data_summ_std$sd <- data_summ_std$sd * data_perm_summ$AdjAmt

  # avgVI_perm_std <- avgVI_perm
  # avgVI_perm_std[-1] <- avgVI_perm_std[-1] * data_perm_summ$AdjAmt
  #
  # # Get Interpolation cutoff
  # interpolationCO <- apply(
  #   apply(avgVI_perm_std[-1],
  #         MARGIN = 2, FUN = sort,
  #         decreasing = TRUE
  #   ),
  #   MARGIN = 1, FUN = stats::quantile, probs = 1 - alpha
  # )

  # Model accuracy estimates
  data_modelAcc <- .computeModelAccuracy(oobAcc = oobAcc,
                                         outcomes = data[, outcome])


  ## Get plots and organize results
  append(
    list(
      "VariableImportance" = data_summ_std,
      "AccuracyEstimate" = data_modelAcc,
      "NoiseCutoff" = noiseCO,
      "InterpolationCutoff" = interpolationCO,
      "AdditionalParams" = list(
        "alpha" = alpha,
        "subsetPlotSize" = subsetPlotSize
      )
    ),
    .generateVIPlot(
      viData = data_summ_std,
      accData = data_modelAcc,
      NoiseCutoff = noiseCO,
      InterpolationCutoff = interpolationCO,
      subsetPlotSize = subsetPlotSize
    )
  )
}
