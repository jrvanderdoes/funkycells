#' Get Agent Count Data
#'
#' This function gets the average percent agent counts per replicate, if there
#'  are replicates (i.e. replicate is not NULL), then the agent percents are
#'  calculated for each replicate and these percentages are averaged.
#'
#' @param agent_data Data.frame of agent data information, with columns as
#'  defined in subsequent parameters
#' @param outcome String of the column name in data indicating the outcome or
#'  response.
#' @param unit String of the column name in data indicating a unit or
#'     base thing. Note this unit may have replicates.
#' @param replicate (Optional) String of the column name in data indicating the
#'  replicate id. Default is NULL.
#' @param type (Optional) String of the column name in data indicating the
#'  type. Default is type.
#' @param data_append (Optional) Data.frame with outcome, patient that the
#'  results can be appended to if desired. Default is NULL.
#'
#' @return List with two elements:
#'     \itemize{
#'         \item dat: Data.frame with outcome, unit, data_append, and the count
#'    data. Columns of the count data are named after the type and are given
#'    in the next list entry.
#'         \item agents: Vector of the the types, i.e. the column names for the
#'    new count data. This can be treated as meta data for funkyForest.
#'     }
#' @export
#'
#' @examples
#' data_ct <- getCountData(TNBC_pheno[TNBC_pheno$Phenotype %in% c('Tumor','B'),],
#'                         outcome="Class", unit="Person",type="Phenotype")
getCountData <- function(agent_data, outcome, unit, replicate = NULL,
                         type = "type", data_append = NULL) {
  # Setup Data
  results <- unique(agent_data[, c(outcome, unit)])
  results_tmp <- unique(agent_data[, c(outcome, unit, replicate)])

  types <- unique(agent_data[[type]])

  results <- cbind(
    results,
    matrix(ncol = length(types))
  )
  colnames(results) <- c(outcome, unit, types)

  # Get Average Counts
  for (i in 1:nrow(results_tmp)) {
    for (ct in types) {
      if (!is.null(replicate)) {
        results_tmp[i, ct] <-
          nrow(agent_data[agent_data[, unit] == results_tmp[i, unit] &
            agent_data[, replicate] == results_tmp[i, replicate] &
            agent_data[, type] == ct, ])
      } else {
        results_tmp[i, ct] <-
          nrow(agent_data[agent_data[, unit] == results_tmp[i, unit] &
            agent_data[, type] == ct, ])
      }
    }
  }
  # Make Percents
  results_tmp[!(colnames(results_tmp) %in% c(outcome, unit, replicate))] <-
    results_tmp[!(colnames(results_tmp) %in% c(outcome, unit, replicate))] /
      rowSums(results_tmp[!(colnames(results_tmp) %in% c(outcome, unit, replicate))])

  for (i in 1:nrow(results)) {
    results[i, -c(1:2)] <- colMeans(
      results_tmp[
        results_tmp[[outcome]] == results[i, outcome] &
          results_tmp[[unit]] == results[i, unit],
        !(colnames(results_tmp) %in%
          c(outcome, unit, replicate))
      ]
    )
  }
  rownames(results) <- NULL

  # Return Data
  if (!is.null(data_append)) {
    return(list(
      "dat" = merge(data_append, results, keep.x = TRUE),
      "types" = types
    ))
  }

  list("dat" = results, "types" = types)
}
