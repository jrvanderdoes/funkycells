#' Get Cell Count Data
#'
#' This (underwork) function gets the average percent cell counts per image, if
#'  there are repeated images (i.e. repeatedId is not NULL), then the cell
#'  percents are calculated for each image and then then these percentages are
#'  averaged. The method is classified as underwork because it requires the
#'  'cellType' column rather than having multiple options.
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
  results_tmp <- unique(cell_data[, c(outcome, unit, repeatedId)])

  cellTypes <- unique(cell_data[["cellType"]])

  results <- cbind(
    results,
    matrix(ncol = length(cellTypes))
  )
  colnames(results) <- c(outcome, unit, cellTypes)

  # Get Average Counts
  for (i in 1:nrow(results_tmp)) {
    for (ct in cellTypes) {
      if(!is.null(repeatedId)){
        results_tmp[i, ct] <-
          nrow(cell_data[cell_data[,unit] == results_tmp[i,unit] &
                           cell_data[,repeatedId] == results_tmp[i,repeatedId] &
                           cell_data[,"cellType"] == ct, ])
      }else{
        results_tmp[i, ct] <-
          nrow(cell_data[cell_data[,unit] == results_tmp[i,unit] &
                           cell_data[,"cellType"] == ct, ])

      }
    }
  }
  # Make Percents
  results_tmp[!(colnames(results_tmp) %in% c(outcome, unit, repeatedId))] <-
    results_tmp[!(colnames(results_tmp) %in% c(outcome, unit, repeatedId))] /
    rowSums(results_tmp[!(colnames(results_tmp) %in% c(outcome, unit, repeatedId))])

  for(i in 1:nrow(results)){
    results[i,-c(1:2)] <- colMeans(
      results_tmp[results_tmp[[outcome]]==results[i,outcome] &
                    results_tmp[[unit]]==results[i,unit],
                  !(colnames(results_tmp) %in%
                      c(outcome,unit,repeatedId))])
  }
  rownames(results) <- NULL

  # Return Data
  if (!is.null(data_append)) {
    return(list(
      "dat" = merge(data_append, results, keep.x=TRUE),
      "cells" = cellTypes
    ))
  }

  list("dat" = results, "cells" = cellTypes)
}
