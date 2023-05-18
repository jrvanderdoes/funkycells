#' Fit a Random Forest model with PC data (Using CV for Improvements)
#'
#' The function fits a random forest model to the data along with using cross-
#'     validation to quantify variable importance. Warning, if there are no
#'     synthetics, this may break (will fix it eventually).
#'
#' @param data Data.frame of outcome and predictors (PCs and meta-variables).
#'     Generally use the results from getPCAData, potentially with meta-
#'     variables attached.
#' @param K (Optional) Numeric indicating the number of folds to use in K-fold
#'     CV. The default is 10.
#' @param metaNames (Optional) Vector indicating the meta-variables to be
#'     considered. Default is NULL.
#' @param synthetics (Optional)
#' @param alpha (Optional) Numeric in (0,1) indicating the significance used
#'     throughout the analysis. Default is 0.05.#'
#' @param silent (Optional) Boolean indicating if output should be suppressed
#'     when the function is running. Default is FALSE.
#' @param rGuessSims (Optional) Numeric value indicating the number of
#'     simulations used for guessing and creating the guess estimate on the
#'     plot. Default is 500.
#' @param subsetPlotSize (Optional) Numeric indicating the number of top
#'     variables to include in a subset graph (note if there are less variables)
#'     than this value indicates then no subset graph will be produced. Default
#'     is 25.
#' @inheritParams funkyForest
#'
#' @return List with the following items:
#'     \enumerate{
#'         \item VariableImportance: Data.frame with the results of variable
#'                   importance indices from the models and CV. The columns are
#'                   var, est, and sd. The columns lower and upper are made with
#'                   significance alpha.
#'         \item AccuracyEstimate: Data.frame with model accuracy estimates:
#'                   out-of-bag accuracy (OOB), biased estimate (bias), and
#'                   random guess (guess). The columns are OOB, bias, and guess.
#'         \item NoiseCutoff: Numeric indicating noise cutoff (vertical line)
#'         \item InterpolationCutoff: Vector of numerics indicating the
#'                   interpolation cutoff (curved line)
#'         \item AdditionalParams: List of additional params for reference:
#'                   Alpha and subsetPlotSize.
#'         \item viPlot: ggplot2 object for vi plot with standardized results.
#'                   It displays ordered underlying functions and meta-variables
#'                   with point estimates, sd, noise cutoff, and interpolation
#'                   cutoff all based on variable importance values
#'         \item subset_viPlot: ggplot2 object for vi plot with standardized
#'                   results and only top subsetPlotSize variables. It displays
#'                   ordered underlying functions and meta-variables with point
#'                   estimates, sd, noise cutoff, and interpolation cutoff all
#'                   based on variable importance values
#'         }
#' @export
#'
#' @examples
#' \dontrun{
#' dat <- simulatePP(
#'   cellVarData =
#'     data.frame("stage" = c(0, 1), "A" = c(0, 0), "B" = c(1 / 50, 1 / 50)),
#'   cellKappaData = data.frame(
#'     "cell" = c("A", "B"),
#'     "clusterCell" = c(NA, "A"),
#'     "kappa" = c(20, 5)
#'   ),
#'   peoplePerStage = 50,
#'   imagesPerPerson = 1,
#'   silent = FALSE
#' )
#' pcaData <- getPCAData(dat,
#'   repeatedUniqueId = "Image",
#'   xRange = c(0, 1), yRange = c(0, 1), silent = FALSE
#' )
#' pcaMeta <- simulateMeta(pcaData,
#'   metaInfo = data.frame(
#'     "var" = c("randUnif", "randBin", "corrNorm"),
#'     "rdist" = c("runif", "rbinom", "rnorm"),
#'     "Stage_0" = c("0.5", "0.5", "1"),
#'     "Stage_1" = c("0.5", "0.5", "2")
#'   )
#' )
#' rfcv <- funkyModel(
#'   data = pcaMeta, outcome = "Stage", unit = "Person",
#'   metaNames = c("randUnif", "randBin", "corrNorm")
#' )
#'
#' dat <- simulatePP(
#'   cellVarData =
#'     data.frame("stage" = c(0, 1), "A" = c(0, 0), "B" = c(1 / 50, 1 / 100)),
#'   cellKappaData = data.frame(
#'     "cell" = c("A", "B"),
#'     "clusterCell" = c(NA, "A"),
#'     "kappa" = c(20, 5)
#'   ),
#'   peoplePerStage = 20,
#'   imagesPerPerson = 3,
#'   silent = FALSE
#' )
#' pcaData <- getPCAData(dat,
#'   repeatedUniqueId = "Image",
#'   xRange = c(0, 1), yRange = c(0, 1), silent = FALSE
#' )
#' pcaMeta <- simulateMeta(pcaData,
#'   metaInfo = data.frame(
#'     "var" = c("randUnif", "randBin", "corrNorm"),
#'     "rdist" = c("runif", "rbinom", "rnorm"),
#'     "Stage_0" = c("0.5", "0.5", "1"),
#'     "Stage_1" = c("0.5", "0.5", "2")
#'   )
#' )
#' rfcv <- funkyModel(
#'   data = pcaMeta, outcome = "Stage", unit = "Person",
#'   metaNames = c("randUnif", "randBin", "corrNorm"),
#'   subsetPlotSize = 2
#' )
#' }
funkyModel <- function(data, K = 10,
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


# funkyModel <- function(data, K = 10,
#                        outcome = colnames(data)[1],
#                        unit = colnames(data)[2],
#                        metaNames = NULL,
#                        synthetics = 100,
#                        alpha = 0.05,
#                        silent = FALSE,
#                        rGuessSims = 500,
#                        subsetPlotSize = 25, nTrees = 500,
#                        method = "class") {
#   ## From moving code over, will remove
#   repeatedId <- NULL
#
#   ## Error checking
#   # .checkData(alignmentMethod) ## TODO:: Add something in
#
#   ## Generate Synthetics And Connect
#   components <- colnames(data)[!(colnames(data) %in%
#                                    c(outcome, unit, repeatedId, metaNames))]
#   KFunctions <- .getUnderlyingVariable(components)
#   # Get Var and data column alignment
#   underlyingDataAlignedFunctions <- .getUnderlyingVariable(colnames(data),
#     returnUnique = FALSE
#   )
#   underlyingVars <- unique(underlyingDataAlignedFunctions)
#   underlyingVars <- underlyingVars[!(underlyingVars %in%
#     c(outcome, unit, repeatedId))]
#
#   avgVI <- avgVI_full <- data.frame("var" = underlyingVars)
#   oobAcc <- rep(NA, K)
#   groups <- .getFolds(1:nrow(data), K)
#
#   # K-Fold Cross-Validation (True Data)
#   if (!silent) cat("CV Trials (", K, "): ", sep = "")
#   for (i in 1:K) {
#     if (!silent) cat(i, ", ", sep = "")
#     RF <- funkyForest(
#       data = data[-groups[[i]], ],
#       outcome = outcome,
#       unit = unit, repeatedId = repeatedId,
#       varImpPlot = FALSE,
#       metaNames = c(metaNames),
#       nTrees = nTrees, method = method
#     )
#
#     data_merge <- RF$varImportanceData[, c("var", "avgVI")]
#     colnames(data_merge) <- c("var", paste0("avgVIK", i))
#     avgVI <- merge(avgVI, data_merge, by = "var")
#
#     oobAcc[i] <- sum(data[groups[[i]], outcome] ==
#       predict_funkyForest(
#         model = RF$model,
#         data_pred = data[groups[[i]], ],
#         type = "pred", data = data
#       )) /
#       nrow(data[groups[[i]], ])
#
#     # Run on all for VI estimate
#     RF_full <- funkyForest(
#       data = data,
#       outcome = outcome,
#       unit = unit, repeatedId = repeatedId,
#       varImpPlot = FALSE,
#       metaNames = c(metaNames),
#       nTrees = nTrees
#     )
#
#     data_merge <- RF_full$varImportanceData[, c("var", "avgVI")]
#     colnames(data_merge) <- c("var", paste0("avgVIK", i))
#     avgVI_full <- merge(avgVI_full, data_merge, by = "var")
#   }
#   if (!silent) cat("\n")
#
#   ## Permutation
#   avgVI_perm <- data.frame("var" = underlyingVars)
#
#   if (!silent) cat("Permutation Trials (", synthetics, "): ", sep = "")
#   for (sim in 1:synthetics) {
#     if (!silent) cat(sim, ", ", sep = "")
#
#     # Permute but do the functional components together
#     data_permute <- data[colnames(data) %in% c(outcome, unit, repeatedId)]
#     data_permute <- cbind(
#       data_permute,
#       as.data.frame(sapply(1:length(underlyingVars),
#         function(idx, DF) {
#           DF[sample.int(nrow(DF)),
#             underlyingDataAlignedFunctions == underlyingVars[idx],
#             drop = FALSE
#           ]
#         },
#         simplify = "array",
#         DF = data
#       ))
#     )
#
#     # Get RF and VI
#     RF <- funkyForest(
#       data = data_permute,
#       outcome = outcome,
#       unit = unit, repeatedId = repeatedId,
#       varImpPlot = FALSE,
#       metaNames = metaNames,
#       nTrees = nTrees, keepModels = FALSE
#     )
#
#     data_merge <- RF$varImportanceData[, c("var", "avgVI")]
#     colnames(data_merge) <- c("var", paste0("avgVIK", sim))
#     avgVI_perm <- merge(avgVI_perm, data_merge, by = "var")
#   }
#   if (!silent) cat("\n")
#
#   # Summarize Data
#   data_summ <- .summData(avgVI_full)[, c("var", "est")]
#   tmp <- .summData(avgVI)
#   data_summ <- merge(
#     data_summ,
#     tmp[, c("var", "sd")]
#   )
#
#   data_perm_summ <- .summData(avgVI_perm)
#   tmp <- as.data.frame(t(sapply(
#     X = data_perm_summ$var,
#     FUN = function(var, alpha, data_vi) {
#       c(
#         var,
#         unlist(stats::quantile(as.numeric(data_vi[data_vi$var == var, -1]),
#           probs = 1 - alpha
#         ))
#       )
#     }, alpha = alpha, data_vi = avgVI_perm, USE.NAMES = FALSE
#   )))
#   colnames(tmp) <- c("var", "upAlpha")
#   tmp$upAlpha <- as.numeric(tmp$upAlpha)
#   data_perm_summ <- merge(data_perm_summ, tmp)
#
#   ## Standardize Data
#   # Get Noise cutoff
#   cutoffs_noise <- data.frame(
#     "Type" = c("KFunction", metaNames),
#     "Cutoff" = c(
#       max(data_perm_summ[data_perm_summ$var %in% KFunctions, "upAlpha"]),
#       rep(NA, length(metaNames))
#     ),
#     "Adj" = NA
#   )
#
#   for (meta in metaNames) {
#     cutoffs_noise[cutoffs_noise$Type == meta, "Cutoff"] <-
#       data_perm_summ[data_perm_summ$var == meta, "upAlpha"]
#   }
#   cutoffs_noise$Adj <- max(cutoffs_noise$Cutoff) / cutoffs_noise$Cutoff
#
#   noiseCO <- max(cutoffs_noise$Cutoff)
#
#   data_perm_summ[data_perm_summ$var %in% KFunctions, "AdjAmt"] <-
#     cutoffs_noise[cutoffs_noise$Type == "KFunction", "Adj"]
#   if (length(metaNames) > 0) {
#     for (var in metaNames) {
#       data_perm_summ[data_perm_summ$var == var, "AdjAmt"] <-
#         cutoffs_noise[cutoffs_noise$Type == var, "Adj"]
#     }
#   }
#
#   data_perm_summ$adj <- data_perm_summ$est * data_perm_summ$AdjAmt
#   data_perm_summ$adj <-
#     ifelse(is.infinite(data_perm_summ$adj), 0, data_perm_summ$adj)
#
#   data_summ_std <- data_summ
#   data_summ_std$est <- data_summ_std$est * data_perm_summ$AdjAmt
#   data_summ_std$sd <- data_summ_std$sd * data_perm_summ$AdjAmt
#
#   avgVI_perm_std <- avgVI_perm
#   avgVI_perm_std[-1] <- avgVI_perm_std[-1] * data_perm_summ$AdjAmt
#
#   # Get Interpolation cutoff
#   interpolationCO <- apply(
#     apply(avgVI_perm_std[-1],
#       MARGIN = 2, FUN = sort,
#       decreasing = TRUE
#     ),
#     MARGIN = 1, FUN = stats::quantile, probs = 1 - alpha
#   )
#
#   # Model accuracy estimates
#   data_modelAcc <- .computeModelAccuracy(oobAcc = oobAcc,
#                                          outcomes = data[, outcome])
#
#
#   ## Get plots and organize results
#   append(
#     list(
#       "VariableImportance" = data_summ_std,
#       "AccuracyEstimate" = data_modelAcc,
#       "NoiseCutoff" = noiseCO,
#       "InterpolationCutoff" = interpolationCO,
#       "AdditionalParams" = list(
#         "alpha" = alpha,
#         "subsetPlotSize" = subsetPlotSize
#       )
#     ),
#     .generateVIPlot(
#       viData = data_summ_std,
#       accData = data_modelAcc,
#       NoiseCutoff = noiseCO,
#       InterpolationCutoff = interpolationCO,
#       subsetPlotSize = subsetPlotSize
#     )
#   )
# }

#' Get Data Summary
#'
#' This is an internal functional for summarizing data from different iterations.
#'
#' @param dat Data.frame with columns for different iterations, rows are the
#'     variables.
#'
#' @return Data.frame with var, est, and sd as columns
#' @noRd
.summData <- function(dat) {
  data.frame(
    "var" = dat$var,
    "est" = rowMeans(dat[-1]),
    "sd" = apply(dat[-1], MARGIN = 1, FUN = function(x) {
      stats::sd(x)
    })
  )
}

#' Compute Model Accuracy
#'
#' This (internal) function compute model accuracy based on results from CV.
#'
#' @param oobAcc Vector of OOB for each iteration
#' @param outcomes Vector of outcomes for the data modeled
#' @param rGuessSims Integer indicating the number of iterations used in guess
#'     method
#'
#' @return Data.frame with OOB, guess, and bias accuracy estimates
#' @noRd
.computeModelAccuracy <- function(oobAcc, outcomes, rGuessSims = 500) {
  optVals <- as.numeric(table(outcomes))
  n <- length(outcomes)

  # OOB
  if (!is.null(oobAcc)) {
    accData <- data.frame("OOB" = mean(oobAcc))
  } else {
    accData <- NULL
  }

  if (!is.null(outcomes)) {
    # Bias guess - Pick most popular group
    accData$bias <- max(optVals) / sum(optVals)

    # Random guessing based on total sample
    acc <- rep(NA, rGuessSims)
    for (i in 1:rGuessSims) {
      # Columns relate to outcomes
      guesses <- stats::rmultinom(length(outcomes),
                                  size = 1,
                                  optVals / sum(optVals))
      acc[i] <- sum(which(guesses == 1, arr.ind = TRUE)[, 1] ==
        as.integer(as.factor(outcomes))) / n
    }
    accData$guess <- mean(acc)
  }

  accData
}

#' Generate VI Plot
#'
#' This (internal) function is used in creation of the VI plot.
#'
#' @param viData Data.frame with var, est, and sd
#' @param accData Data.frame with oob, bias, and guess estimates of model
#'     accuracy
#' @param NoiseCutoff Numeric indicating the vertical noise cutoff
#' @param InterpolationCutoff Vector of numerics for the curved cutoff
#' @param subsetPlotSize Integer indicating number of interactions to have on
#'     the subset plot. Note if there are fewer in the base, no subset plot will
#'     be created
#'
#' @return List with VI figures (one or two)
#' @noRd
.generateVIPlot <- function(viData, accData, NoiseCutoff,
                            InterpolationCutoff, subsetPlotSize) {
  viPlot <- .plotVI(viData, accData, NoiseCutoff, InterpolationCutoff)

  data_return <- list("viPlot" = viPlot)

  if (nrow(viData) > subsetPlotSize) {
    # Order Data to take top
    tmpStdVI <- viData[order(-viData$est), ]

    viPlot <- .plotVI(
      tmpStdVI[1:subsetPlotSize, ], accData, NoiseCutoff,
      InterpolationCutoff[1:subsetPlotSize]
    )

    data_return <- append(
      data_return,
      list("subset_viPlot" = viPlot)
    )
  }

  data_return
}


#' Plot VI
#'
#' This (internal) function standardizing and plots the variable importance
#'     values.
#'
#' @param viData Data.frame with var, est, and sd
#' @param accData Data.frame with oob, bias, and guess estimates of model
#'     accuracy
#' @param NoiseCutoff Numeric indicating the vertical noise cutoff
#' @param InterpolationCutoff Vector of numerics for the curved cutoff
#'
#' @return ggplot2 figure object for the VI plot
#' @noRd
.plotVI <- function(viData, accData, NoiseCutoff,
                    InterpolationCutoff) {
  # Add this to stop NOTEs in building package
  est <- sd <- var <- NULL

  maxVal <- max(InterpolationCutoff, NoiseCutoff, viData$est)

  ggplot2::ggplot(
    data = viData,
    mapping = ggplot2::aes(
      x = factor(stats::reorder(var, est)),
      y = ifelse(est / maxVal > 1, 1,
        ifelse(est / maxVal < 0, 0,
          est / maxVal
        )
      ),
      group = 1
    )
  ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(
        ymin = ifelse((est - sd) / maxVal < 0, 0, (est - sd) / maxVal),
        ymax = ifelse((est + sd) / maxVal > 1, 1, (est + sd) / maxVal)
      ),
      color = "black", width = 0.2
    ) +
    ggplot2::geom_point(color = "black") +
    ggplot2::geom_hline(
      ggplot2::aes(yintercept = max(0, min(1, NoiseCutoff / maxVal))),
      color = "red", linetype = "dotted", linewidth = 1
    ) +
    ggplot2::coord_flip(ylim = c(0, 1)) +
    ggplot2::xlab(NULL) +
    ggplot2::ylim(c(0, 1)) +
    ggplot2::ylab(NULL) +
    ggplot2::theme_bw() +
    ggplot2::geom_line(ggplot2::aes(
      x = ordered(viData[order(-est), "var"]),
      y = InterpolationCutoff / maxVal
    ), color = "orange", linetype = "dashed", size = 1) +
    ggplot2::ylab(paste0(
      "Variable Importance - OOB (",
      .specify_decimal(accData$OOB, 2),
      "), Guess (",
      .specify_decimal(accData$guess, 2),
      "), Bias (",
      .specify_decimal(accData$bias, 2), ")"
    ))
}
