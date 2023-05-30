#' Get Cell Count Data
#'
#' This (underwork) function gets the average cell counts per image. It is
#'   classified as underwork because it requires the 'cellType' column.
#'
#' @param cell_data Data.frame of cell data information. Must have cellType
#'  column. Future extensions will remove this restriction.
#' @param data_append (Optional) Data.frame with outcome, patient that the
#'  results can be appended to if desired. Default is NULL.
#' @inheritParams funkyForest
#'
#' @return List with two elements:
#'  1. dat: Data.frame with outcome, unit, any data_append and the count
#'    data. Columns of the count data are named after the cellType and are given
#'    in the next list entry.
#'  2. cells: Vector of cellTypes, i.e. the column names for the new count
#'    data. This is treated as meta data for funkyForest.
#' @export
#'
#' @examples
#' dat <- simulatePP(
#'   cellVarData =
#'     data.frame(
#'       "stage" = c(0, 1),
#'       "A" = c(0, 0),
#'       "B" = c(1 / 50, 1 / 100),
#'       "C" = c(1 / 100, 1 / 100)
#'     ),
#'   cellKappaData = data.frame(
#'     "cell" = c("A", "B", "C"),
#'     "clusterCell" = c(NA, "A", NA),
#'     "kappa" = c(20, 5, 10)
#'   ),
#'   peoplePerStage = 10,
#'   imagesPerPerson = 2,
#'   silent = FALSE
#' )
#'
#' # Use Case 1
#' data_ct <- getCountData(dat, "Stage", "Person", "Image")
#'
#' \dontrun{
#' dat_mod <- funkyModel(
#'   data = data_ct$dat,
#'   outcome = "Stage", unit = "Person",
#'   metaNames = data_ct$cells,
#'   synthetics = 25, nTrees = 50
#' )
#' }
#'
#'
#' # Use Case 2
#' \dontrun{
#' dat <- simulatePP(
#'   cellVarData =
#'     data.frame(
#'       "stage" = c(0, 1),
#'       "A" = c(0, 0),
#'       "B" = c(1 / 50, 1 / 100),
#'       "C" = c(1 / 100, 1 / 100),
#'       "D" = c(1 / 50, 1 / 100)
#'     ),
#'   cellKappaData = data.frame(
#'     "cell" = c("A", "B", "C", "D"),
#'     "clusterCell" = c(NA, "A", NA, "B"),
#'     "kappa" = c(20, 5, 10, 4)
#'   ),
#'   peoplePerStage = 15,
#'   imagesPerPerson = 2,
#'   silent = FALSE
#' )
#'
#' pcaData <- getKsPCAData(dat,
#'   repeatedUniqueId = "Image",
#'   xRange = c(0, 1), yRange = c(0, 1), silent = FALSE
#' )
#' data_ct1 <- getCountData(dat, "Stage", "Person", "Image", data_append = pcaData)
#'
#' dat_mod <- funkyModel(
#'   data = data_ct1$dat,
#'   outcome = "Stage", unit = "Person",
#'   metaNames = data_ct1$cells,
#'   synthetics = 25, nTrees = 50
#' )
#' }
getCountData <- function(cell_data, outcome, unit, repeatedId = NULL,
                         data_append = NULL) {
  # Setup Data
  results <- unique(cell_data[, c(outcome, unit)])

  units <- unique(cell_data[[unit]])
  cellTypes <- unique(cell_data[["cellType"]])

  results <- cbind(
    results,
    matrix(ncol = length(cellTypes))
  )
  colnames(results) <- c(outcome, unit, cellTypes)

  # Get Average Counts
  for (u in units) {
    for (ct in cellTypes) {
      results[results[[unit]] == u, ct] <-
        nrow(cell_data[cell_data[[unit]] == u & cell_data[["cellType"]] == ct, ])
    }

    if (!is.null(repeatedId)) {
      results[results[[unit]] == u, -c(1:2)] <- results[results[[unit]] == u, -c(1:2)] /
        length(unique(cell_data[cell_data[[unit]] == u, repeatedId]))
    }
  }
  rownames(results) <- NULL

  # Return Data
  if (!is.null(data_append)) {
    return(list(
      "dat" = merge(data_append, results),
      "cells" = cellTypes
    ))
  }

  list("dat" = results, "cells" = cellTypes)
}
