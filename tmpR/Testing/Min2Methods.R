
funkyModel1 <- function(data, K = 10,
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
  avgVI_full_old <-data.frame("var"=c(underlyingVars))
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
    RF_full_old <- funkyForest(
      data = data,
      outcome = outcome,
      unit = unit, repeatedId = repeatedId,
      varImpPlot = FALSE,
      metaNames = metaNames,
      nTrees = nTrees,keepModels = FALSE
    )

    data_merge <- RF_full$varImportanceData[, c("var", "avgVI")]
    colnames(data_merge) <- c("var", paste0("avgVIK", i))
    avgVI_full <- merge(avgVI_full, data_merge, by = "var")

    data_merge <- RF_full_old$varImportanceData[, c("var", "avgVI")]
    colnames(data_merge) <- c("var", paste0("avgVIK", i))
    avgVI_full_old <- merge(avgVI_full_old, data_merge, by = "var")
  }
  if (!silent) cat("\n")

  ## Permutation
  avgVI_perm <-
    data.frame("var" = c(underlyingVars,
                         paste0('permuteInternal',1:synthetics,'K_K'),
                         as.vector(sapply(metaNames,function(m){
                           paste0(paste0('permuteInternal',1:synthetics),m)}))))
  avgVI_perm_old <-
    data.frame("var" = c(underlyingVars))

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
    RF_old <- funkyForest(
      data = cbind(data[colnames(data) %in%
                          c(outcome, unit, repeatedId)],
                   data_permute1),
      outcome = outcome,
      unit = unit, repeatedId = repeatedId,
      varImpPlot = FALSE,
      metaNames = c(metaNames),
      nTrees = nTrees, keepModels = FALSE
    )

    data_merge <- RF$varImportanceData[, c("var", "avgVI")]
    colnames(data_merge) <- c("var", paste0("avgVIK", sim))
    avgVI_perm <- merge(avgVI_perm, data_merge, by = "var")

    data_merge <- RF_old$varImportanceData[, c("var", "avgVI")]
    colnames(data_merge) <- c("var", paste0("avgVIK", sim))
    avgVI_perm_old <- merge(avgVI_perm_old, data_merge, by = "var")
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

  ## PERMUTATION
  getOld_intCO <- function(avgVI_full_old, avgVI,
                           avgVI_perm_old,
                           KFunctions, metaNames){
    # Summarize Data
    data_summ_old <- .summData(avgVI_full_old)[, c("var", "est")]
    tmp <- .summData(avgVI)
    data_summ_old <- merge(
      data_summ_old,
      tmp[, c("var", "sd")]
    )

    data_perm_summ_old <- .summData(avgVI_perm_old)
    tmp <- as.data.frame(t(sapply(
      X = data_perm_summ_old$var,
      FUN = function(var, alpha, data_vi) {
        c(
          var,
          unlist(stats::quantile(as.numeric(data_vi[data_vi$var == var, -1]),
                                 probs = 1 - alpha
          ))
        )
      }, alpha = alpha, data_vi = avgVI_perm_old, USE.NAMES = FALSE
    )))
    colnames(tmp) <- c("var", "upAlpha")
    tmp$upAlpha <- as.numeric(tmp$upAlpha)
    data_perm_summ_old <- merge(data_perm_summ_old, tmp)

    ## Standardize Data
    # Get Noise cutoff
    cutoffs_noise <- data.frame(
      "Type" = c("KFunction", metaNames),
      "Cutoff" = c(
        max(data_perm_summ_old[data_perm_summ_old$var %in% KFunctions, "upAlpha"]),
        rep(NA, length(metaNames))
      ),
      "Adj" = NA
    )

    for (meta in metaNames) {
      cutoffs_noise[cutoffs_noise$Type == meta, "Cutoff"] <-
        data_perm_summ_old[data_perm_summ_old$var == meta, "upAlpha"]
    }
    cutoffs_noise$Adj <- max(cutoffs_noise$Cutoff) / cutoffs_noise$Cutoff

    noiseCO_old <- max(cutoffs_noise$Cutoff)

    data_perm_summ[data_perm_summ$var %in% KFunctions, "AdjAmt"] <-
      cutoffs_noise[cutoffs_noise$Type == "KFunction", "Adj"]
    if (length(metaNames) > 0) {
      for (var in metaNames) {
        data_perm_summ[data_perm_summ$var == var, "AdjAmt"] <-
          cutoffs_noise[cutoffs_noise$Type == var, "Adj"]
      }
    }

    data_perm_summ$adj <- data_perm_summ$est * data_perm_summ$AdjAmt
    data_perm_summ$adj <-
      ifelse(is.infinite(data_perm_summ$adj), 0, data_perm_summ$adj)

    data_summ_std <- data_summ_old
    data_summ_std$est <- data_summ_std$est * data_perm_summ$AdjAmt
    data_summ_std$sd <- data_summ_std$sd * data_perm_summ$AdjAmt

    avgVI_perm_std <- avgVI_perm_old
    avgVI_perm_std[-1] <- avgVI_perm_std[-1] * data_perm_summ$AdjAmt

    # Get Interpolation cutoff
    interpolationCO <- apply(
      apply(avgVI_perm_std[-1],
            MARGIN = 2, FUN = sort,
            decreasing = TRUE
      ),
      MARGIN = 1, FUN = stats::quantile, probs = 1 - alpha
    )

    interpolationCO
  }

  interpolation2 <- getOld_intCO(avgVI_full_old, avgVI,
                                 avgVI_perm_old,
                                 KFunctions, metaNames)

  ## FIX
  interpolationCO <- pmin(interpolationCO,interpolation2)

  # Model accuracy estimates
  data_modelAcc <- .computeModelAccuracy(oobAcc = oobAcc,
                                         outcomes = data[, outcome])


  ## Get plots and organize results
  append(
    list(
      "model"= funkyForest(
        data = data,
        outcome = outcome,
        unit = unit,
        repeatedId = repeatedId,
        varImpPlot = FALSE,
        metaNames = c(metaNames),
        nTrees = nTrees,
        method = method
      )$model,
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


## 20 Patients
set.seed(123)
cell_data_20 <- simulatePP(
  cellVarData =
    data.frame(
      "stage" = c(0, 1),
      "A" = c(0, 0),
      "B" = c(1 / 100, 1 / 500),
      "C" = c(1 / 500, 1 / 250),
      "D" = c(1 / 100, 1 / 100)
    ),
  cellKappaData =
    data.frame(
      "cell" = c("A", "B", "C", "D"),
      "clusterCell" = c(NA, "A", NA, NA),
      "kappa" = c(20, 5, 15, 15)
    ),
  peoplePerStage = 15,
  imagesPerPerson = 1,
  silent = FALSE
)

cells <- unique(cell_data_20$cellType)
cells_interactions <- rbind(data.frame(t(combn(cells,2))),
                            data.frame('X1'=cells,'X2'=cells))

pca_data_20 <- getKsPCAData(
  data = cell_data_20,
  outcome = 'Stage',
  unit = "Person",
  repeatedUniqueId = 'Image',
  rCheckVals = seq(0,0.25,0.01),
  agents_df = cells_interactions,
  xRange = c(0,1),  yRange = c(0,1),
  nPCs = 3)


set.seed(123)
pcaMeta_20 <- simulateMeta(pca_data_20,
                           outcome = 'Stage',
                           metaInfo = data.frame(
                             "var" = c("gender", "age"),
                             "rdist" = c("rbinom", "rnorm"),
                             "Stage_0" = c("0.5", "30"),
                             "Stage_1" = c("0.5", "31"))
)


set.seed(123)
model_fc_20 <- funkyModel(data=pcaMeta_20,
                          outcome = 'Stage',
                          unit = 'Person',
                          metaNames=c('gender','age'))

set.seed(123)
model_fc_new_20 <- funkyModel1(data=pcaMeta_20,
                               outcome = 'Stage',
                               unit = 'Person',
                               metaNames=c('gender','age'))
