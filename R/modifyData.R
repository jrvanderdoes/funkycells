#' Get Agent Count Data
#'
#' This (in development) function gets the average percent agent counts per
#'  replicate, if there are replicates (i.e. replicate is not NULL), then the
#'  agent percents are calculated for each replicate and then then these
#'  percentages are averaged. The method is classified as in development.
#'
#' @param agent_data Data.frame of agent data information. Must have type
#'  column. Future extensions will remove this restriction.#'
#' @param outcome String of the column name in data indicating the outcome or
#'  response.
#' @param unit String of the column name in data indicating a unit or
#'     base thing. Note this unit may have replicates.
#' @param replicate (Optional) String of the column name in data indicating the
#'  replicate id. Default is NULL.
#' @param data_append (Optional) Data.frame with outcome, patient that the
#'  results can be appended to if desired. Default is NULL.
#' @inheritParams funkyForest
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
#' dat <- simulatePP(
#'   agentVarData =
#'     data.frame(
#'       "outcome" = c(0, 1),
#'       "A" = c(0, 0),
#'       "B" = c(1 / 50, 1 / 100),
#'       "C" = c(1 / 100, 1 / 100)
#'     ),
#'   agentKappaData = data.frame(
#'     "agent" = c("A", "B", "C"),
#'     "clusterAgent" = c(NA, "A", NA),
#'     "kappa" = c(20, 5, 10)
#'   ),
#'   unitsPerOutcome = 10,
#'   replicatesPerUnit = 2,
#'   silent = FALSE
#' )
#'
#' # Use Case 1
#' data_ct <- getCountData(dat, "outcome", "unit", "replicate")
#'
#' \dontrun{
#' dat_mod <- funkyModel(
#'   data = data_ct$dat,
#'   outcome = "outcome", unit = "unit",
#'   metaNames = data_ct$agents,
#'   synthetics = 25, nTrees = 50
#' )
#' }
#'
#'
#' # Use Case 2
#' \dontrun{
#' dat <- simulatePP(
#'   agentVarData =
#'     data.frame(
#'       "outcome" = c(0, 1),
#'       "A" = c(0, 0),
#'       "B" = c(1 / 50, 1 / 100),
#'       "C" = c(1 / 100, 1 / 100),
#'       "D" = c(1 / 50, 1 / 100)
#'     ),
#'   agentKappaData = data.frame(
#'     "agent" = c("A", "B", "C", "D"),
#'     "clusterAgent" = c(NA, "A", NA, "B"),
#'     "kappa" = c(20, 5, 10, 4)
#'   ),
#'   unitsPerOutcome = 15,
#'   replicatesPerUnit = 2,
#'   silent = FALSE
#' )
#'
#' pcaData <- getKsPCAData(dat,
#'   replicate = "replicate",
#'   xRange = c(0, 1), yRange = c(0, 1), silent = FALSE
#' )
#' data_ct1 <- getCountData(dat, "outcome", "unit", "replicate", data_append = pcaData)
#'
#' dat_mod <- funkyModel(
#'   data = data_ct1$dat,
#'   outcome = "outcome", unit = "unit",
#'   metaNames = data_ct1$agents,
#'   synthetics = 25, nTrees = 50
#' )
#' }
getCountData <- function(agent_data, outcome, unit, replicate = NULL,
                         data_append = NULL) {
  # Setup Data
  results <- unique(agent_data[, c(outcome, unit)])
  results_tmp <- unique(agent_data[, c(outcome, unit, replicate)])

  types <- unique(agent_data[["type"]])

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
            agent_data[, "type"] == ct, ])
      } else {
        results_tmp[i, ct] <-
          nrow(agent_data[agent_data[, unit] == results_tmp[i, unit] &
            agent_data[, "type"] == ct, ])
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
