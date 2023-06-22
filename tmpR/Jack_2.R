permuteComputeIntCurve <- function(data, outcome, unit,
                                   underlyingNoiseVars,
                                   underlyingDataAlignedFunctions,
                                   KFunctions, metaNames, nPCs, synthetics){

  ## Permutation
  avgVI_perm <- list()

  if (!silent) cat("Permutation Trials (", synthetics, "): ", sep = "")
  for (sim in 1:synthetics) {
    if (!silent) cat(sim, ", ", sep = "")
    avgVI_perm[[sim]] <- data.frame("var" = c(underlyingNoiseVars))

    for(topVal in 1:length(underlyingNoiseVars)){
      removeVars <- sample(underlyingNoiseVars,topVal-1)
      data_minusData <- data[,!(underlyingDataAlignedFunctions %in%
                                  removeVars)]
      metaNames_still <- metaNames[metaNames %in% colnames(data_minusData)]

      # Permute but do the functional components together
      data_permute <- .permuteData(data_minusData, outcome, unit,
                                   synthetics, KFunctions, metaNames_still,
                                   underlyingDataAlignedFunctions, nPCs,
                                   attach.data = TRUE, permute.data = TRUE
      )

      # Get RF and VI
      RF <- funkyForest(
        data = data_permute$data,
        outcome = outcome,
        unit = unit, repeatedId = repeatedId,
        varImpPlot = FALSE,
        metaNames = data_permute$metaNames,
        nTrees = nTrees, keepModels = FALSE
      )

      data_merge <- RF$varImportanceData[, c("var", "avgVI")]
      colnames(data_merge) <- c("var", paste0("avgVIK", sim))
      avgVI_perm[[sim]] <- merge(avgVI_perm[[sim]], data_merge, by = "var")
    }

  }
  if (!silent) cat("\n")

  return(avgVI_perm)
}
