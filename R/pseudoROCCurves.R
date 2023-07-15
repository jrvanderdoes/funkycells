#' Compute Pseudo-ROC Curves
#'
#' An receiver operating characteristic (ROC) curve is a curve showing the
#'  performance of a classification model at all classification thresholds.
#'  True ROC can only be computed for two-options, but we can consider each
#'  classification, i.e. prediction, correct or incorrect and overlay the
#'  curves. Note this means the lines may cover each other and be difficult to
#'  see.
#'
#'  This function requires the package 'pROC' to be installed.
#'
#' @param trueOutcomes Vector of the true results
#' @param modelPercents Data.frame with columns named after the true outcomes,
#'     giving the percent of selecting that outcome. This is what is returned
#'     predict.RandomForest_PC with type='all' in object `PredPerc[-1]` (first
#'     column is the predictions).
#'
#' @return ggplot object containing the ROC curves.
#' @export
#'
#' @examples
#' set.seed(123)
#' data_pp_roc <- simulatePP(
#'   cellVarData =
#'     data.frame(
#'       "stage" = c(0, 1),
#'       "A" = c(0, 0),
#'       "B" = c(1 / 50, 1 / 50)
#'     ),
#'   cellKappaData = data.frame(
#'     "cell" = c("A", "B"),
#'     "clusterCell" = c(NA, "A"),
#'     "kappa" = c(20, 5)
#'   ),
#'   peoplePerStage = 10,
#'   imagesPerPerson = 1,
#'   silent = FALSE
#' )
#' # Caution, in general use more than 5 nTrees (Default is 500)
#' pcaData_roc <- getKsPCAData(data_pp_roc,
#'   repeatedUniqueId = "Image",
#'   xRange = c(0, 1), yRange = c(0, 1), silent = FALSE
#' )
#' RF_roc <- funkyForest(data = pcaData_roc[-2], nTrees = 5)
#' pred_roc <- predict_funkyForest(
#'   model = RF_roc$model,
#'   data_pred = pcaData_roc[-2],
#'   data = pcaData_roc[-2]
#' )
#' computePseudoROCCurves(pcaData_roc$Stage, pred_roc$PredPerc[-1])
computePseudoROCCurves <- function(trueOutcomes, modelPercents) {
  if (!requireNamespace("pROC", quietly = TRUE)) {
    stop(
      "Package \"pROC\" must be installed to use this function.",
      call. = FALSE
    )
  }

  if (ncol(modelPercents) < 2) {
    stop("Error: Only one outcome (col) in modelPercents")
  }

  options_roc <- colnames(modelPercents)
  data_auc <- rep(NA, length(options_roc))
  data_plot <- data.frame(matrix(nrow = 0, ncol = 3))
  colnames(data_plot) <- c("Type", "Spec", "Sens")

  for (i in 1:length(options_roc)) {
    roc.curve <- pROC::roc(
      ifelse(trueOutcomes == options_roc[i],
        options_roc[i], paste0("Not-", options_roc[i])
      ),
      as.numeric(modelPercents[, options_roc[i]]),
      levels = c(options_roc[i], paste0("Not-", options_roc[i])),
      direction = ">"
    )

    data_plot <- rbind(
      data_plot,
      data.frame(
        "Type" = options_roc[i],
        "Spec" = roc.curve$specificities,
        "Sens" = roc.curve$sensitivities
      )
    )
    data_auc[i] <- roc.curve$auc
  }

  # Add this to stop NOTEs in building package
  Spec <- Sens <- Type <- NULL

  ggplot2::ggplot(
    data_plot,
    ggplot2::aes(x = Spec, y = Sens, color = Type)
  ) +
    ggplot2::geom_path(linewidth = 0.75) +
    ggplot2::geom_abline(slope = 1, intercept = 1, color = "gray") +
    ggplot2::xlab("Specificity") +
    ggplot2::ylab("Sensitivity") +
    ggplot2::xlim(c(1, 0)) +
    ggplot2::ylim(c(0, 1)) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid = ggplot2::element_blank()) +
    ggplot2::scale_color_discrete(
      labels = paste0(options_roc, " (AUC: ", data_auc, ")"),
      name = "Outcome"
    )
}
