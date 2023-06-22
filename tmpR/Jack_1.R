#' Fit a Random Forest model with PC data (Using CV for Improvements)
#'
#' The function fits a random forest model to the data along with using cross-
#'     validation to quantify variable importance. Warning, if there are no
#'     synthetics, this may break (will fix it eventually).
#'
#' @param data Data.frame of outcome and predictors (PCs and meta-variables).
#'     Generally use the results from getKsPCAData, potentially with meta-
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
#'
#' # Short Example
#' #   Parameters are reduced beyond recommended levels for speed
#' set.seed(1234)
#' rfcv <- funkyModel(
#'   data = TNBC[, c(1:8, ncol(TNBC))],
#'   outcome = "Class", unit = "Person",
#'   metaNames = c("Age"),
#'   nTrees = 50, synthetics = 10,
#'   silent = TRUE
#' )
#'
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
#' pcaData <- getKsPCAData(dat,
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
#' pcaData <- getKsPCAData(dat,
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
#'
#' \dontrun{
#' set.seed(1234567)
#' dat <- simulatePP(
#'   cellVarData =
#'     data.frame(
#'       "stage" = c(0, 1),
#'       "A" = c(0, 0),
#'       "B" = c(1 / 50, 1 / 100)
#'     ),
#'   cellKappaData = data.frame(
#'     "cell" = c("A", "B"),
#'     "clusterCell" = c(NA, "A"),
#'     "kappa" = c(20, 5)
#'   ),
#'   peoplePerStage = 20,
#'   imagesPerPerson = 3,
#'   silent = FALSE
#' )
#' pcaData <- getKsPCAData(dat,
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
#'
#' rfcv1 <- funkyModel(
#'   data = pcaMeta, outcome = "Stage", unit = "Person",
#'   metaNames = c("randUnif", "randBin", "corrNorm"),
#'   nTrees = 10, synthetics = 25
#' )
#'
#' rfcv2 <- funkyModel(
#'   data = pcaData, outcome = "Stage", unit = "Person",
#'   nTrees = 10, synthetics = 25
#' )
#'
#' onlyMeta <- pcaMeta[, c("Stage", "Person", "randUnif", "randBin", "corrNorm")]
#' rfcv3 <- funkyModel(
#'   data = onlyMeta, outcome = "Stage", unit = "Person",
#'   metaNames = c("randUnif", "randBin", "corrNorm"),
#'   nTrees = 10, synthetics = 25
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
  nPCs <- ifelse(length(components) == 0,
                 0,
                 as.numeric(max(sub(".*_PC", "", components)))
  )

  KFunctions <- .getUnderlyingVariable(components)
  # Get Var and data column alignment
  underlyingDataAlignedFunctions <- .getUnderlyingVariable(colnames(data),
                                                           returnUnique = FALSE
  )
  underlyingVars <- unique(underlyingDataAlignedFunctions)
  underlyingVars <- underlyingVars[!(underlyingVars %in%
                                       c(outcome, unit, repeatedId))]

  underlyingNoiseVars <- c(underlyingVars)
  if (length(KFunctions) > 0) {
    underlyingNoiseVars <- c(
      underlyingNoiseVars,
      paste0("permuteInternal", 1:synthetics, "K_K")
    )
  }
  if (length(metaNames) > 0) {
    underlyingNoiseVars <- c(
      underlyingNoiseVars,
      as.vector(sapply(metaNames, function(m) {
        paste0(paste0("permuteInternal", 1:synthetics), m)
      }))
    )
  }

  avgVI <- data.frame("var" = c(underlyingVars))
  avgVI_full <- data.frame("var" = c(underlyingNoiseVars))

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
    data_full <- .permuteData(
      data_base = data, outcome = outcome, unit = unit,
      synthetics = synthetics,
      KFunctions = KFunctions, metaNames = metaNames,
      underlyingDataAlignedFunctions = underlyingDataAlignedFunctions,
      nPCs = nPCs, attach.data = TRUE, permute.data = FALSE
    )

    RF_full <- funkyForest(
      data = data_full$data,
      outcome = outcome,
      unit = unit, repeatedId = repeatedId,
      varImpPlot = FALSE,
      metaNames = data_full$metaNames,
      nTrees = nTrees, keepModels = FALSE
    )

    data_merge <- RF_full$varImportanceData[, c("var", "avgVI")]
    colnames(data_merge) <- c("var", paste0("avgVIK", i))
    avgVI_full <- merge(avgVI_full, data_merge, by = "var")
  }
  if (!silent) cat("\n")

  ## Permutation
  avgVI_perm <- permuteComputeIntCurve(data, outcome, unit,
                                       underlyingNoiseVars,
                                       underlyingDataAlignedFunctions,
                                       KFunctions, metaNames, nPCs, synthetics)


  # Summarize Data
  # data_summ <- .summData(avgVI_full)[, c("var", "est","sd")]

  tmp <- as.data.frame(t(sapply(
    X = underlyingVars,
    FUN = function(var1, alpha, data_vi, KFunctions) {
      if (var1 %in% KFunctions) {
        tmp <- data.frame(
          var1,
          t(as.numeric(apply(
            X = data_vi[data_vi$var %in% paste0("permuteInternal", 1:100, "K_K"), -1],
            MARGIN = 2, FUN = stats::quantile, probs = 1 - alpha
          )))
        )
      } else {
        tmp <- data.frame(
          var1,
          t(as.numeric(apply(
            X = data_vi[data_vi$var %in% paste0("permuteInternal", 1:100, var1), -1],
            MARGIN = 2, FUN = stats::quantile, probs = 1 - alpha
          )))
        )
      }
      tmp
    },
    alpha = alpha, data_vi = avgVI_full, KFunctions = KFunctions,
    USE.NAMES = FALSE, simplify = "matrix"
  )))
  for (i in 2:ncol(tmp)) {
    tmp[i] <- unlist(tmp[i])
  }

  noiseCO <- max(tmp[-1])

  # If value is 0, will jsut take rather than scale
  tmp1 <- noiseCO / tmp[, -1]
  tmp1[sapply(tmp1, simplify = "matrix", is.infinite)] <- 1

  avgVI_std <- cbind(
    "var" = avgVI_full[avgVI_full$var %in% underlyingVars, 1],
    avgVI_full[avgVI_full$var %in% underlyingVars, -1] * tmp1
  )
  data_summ <- .summData(avgVI_std)

  # Add CV SDs (which do not include any noise Vars)
  tmp <- .summData(avgVI)
  colnames(tmp)[3] <- "cvSD"
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
      if (var %in% KFunctions) {
        tmp <- data.frame(
          var,
          t(as.numeric(apply(
            X = data_vi[data_vi$var %in% paste0("permuteInternal", 1:100, "K_K"), -1],
            MARGIN = 2, FUN = stats::quantile, probs = 1 - alpha
          )))
        )
      } else {
        tmp <- data.frame(
          var,
          t(as.numeric(apply(
            X = data_vi[data_vi$var %in% paste0("permuteInternal", 1:100, var), -1],
            MARGIN = 2, FUN = stats::quantile, probs = 1 - alpha
          )))
        )
      }
      tmp
    },
    alpha = alpha, data_vi = avgVI_perm, KFunctions = KFunctions,
    USE.NAMES = FALSE, simplify = "matrix"
  )))
  for (i in 2:ncol(tmp)) {
    tmp[i] <- unlist(tmp[i])
  }

  # If value is 0, will just take rather than scale
  tmp1 <- noiseCO / tmp[, -1]
  tmp1[sapply(tmp1, simplify = "matrix", is.infinite)] <- 1

  avgVI_perm_std <- cbind(
    "var" = avgVI_perm[avgVI_perm$var %in% underlyingVars, 1],
    avgVI_perm[avgVI_perm$var %in% underlyingVars, -1] * tmp1
  )
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
  data_modelAcc <- .computeModelAccuracy(
    oobAcc = oobAcc,
    outcomes = data[, outcome]
  )


  ## Get plots and organize results
  append(
    list(
      "model" = funkyForest(
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
