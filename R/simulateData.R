#' Simulate a Point Processes
#'
#' This function simulates a point pattern with optional clustering ( visible
#'     and invisible). Multiple outcomes, people, and repeated measure are
#'     possible.
#'
#' Upcoming: cellVarData unnamed. Rename to general names (agent, unit, etc.)
#'
#' @param cellVarData (Optional) Data.frame describing variances with each agent
#'     type.
#'
#' The data.frame has a stage column and a named column for each agent type.
#'     Currently, these names are mandatory.
#' @param cellKappaData (Optional) Data.frame describing cell interactions.
#'
#' The data.frame has a cell column giving agent names (matching cellVarData),
#'     a clusterCell column indicating which agent the agent clusters (put NA
#'     if the agent doesn't cluster or clusters a hidden agent / self-clusters),
#'     and a kappa column directing the number of agents of per image.
#' @param peoplePerStage (Optional) Numeric indicating the number of units per
#'     outcome.
#' @param imagesPerPerson (Optional) Numeric indicating the number of repeated
#'     measures.
#' @param silent (Optional) Boolean indicating if progress output should be
#'     printed.
#'
#' @return Data.frame containing each point the defined patterns.
#'
#' The data.frame has columns for outcome, x coordinate, y coordinate, agent
#'     type, unit, and unique repeated measure id.
#' @export
#'
#' @examples
#' \dontrun{
#' data <- simulatePP()
#' }
#' data <- simulatePP(
#'   cellVarData = data.frame(
#'     "stage" = c(0, 1),
#'     "A" = c(0, 0),
#'     "B" = c(1 / 100, 1 / 500),
#'     "C" = c(1 / 500, 1 / 250),
#'     "D" = c(1 / 100, 1 / 100),
#'     "E" = c(1 / 500, 1 / 500)
#'   ),
#'   cellKappaData = data.frame(
#'     "cell" = c("A", "B", "C", "D", "E"),
#'     "clusterCell" = c(NA, "A", "B", "C", NA),
#'     "kappa" = c(10, 3, 2, 1, 8)
#'   ),
#'   peoplePerStage = 4,
#'   imagesPerPerson = 1
#' )
simulatePP <- function(cellVarData =
                         data.frame(
                           "stage" = c(0, 1, 2),
                           "A" = c(0, 0, 0),
                           "B" = c(1 / 100, 1 / 500, 1 / 500),
                           "C" = c(1 / 500, 1 / 250, 1 / 100),
                           "D" = c(1 / 100, 1 / 100, 1 / 100),
                           "E" = c(1 / 500, 1 / 500, 1 / 500),
                           "F" = c(1 / 250, 1 / 250, 1 / 250)
                         ),
                       cellKappaData = data.frame(
                         "cell" = c("A", "B", "C", "D", "E", "F"),
                         "clusterCell" = c(NA, "A", "B", "C", NA, "A"),
                         "kappa" = c(20, 5, 4, 2, 15, 5)
                       ),
                       peoplePerStage = 20,
                       imagesPerPerson = 5,
                       silent = FALSE) {
  ## Setup
  data_stages <- list()
  data_stages1 <- list()
  # Go through each stage
  for (stageIdx in 1:nrow(cellVarData)) {
    # General Vars
    stage <- cellVarData$stage[stageIdx] # Current Stage
    imageAdj <- (stageIdx - 1) * (peoplePerStage * imagesPerPerson) # Adjustment for images due to stage
    personAdj <- (stageIdx - 1) * (peoplePerStage) # Adjustment for person due to stage

    if (!silent) {
      cat(paste0("Stage: ", stage, " (", stageIdx, "/", nrow(cellVarData), ")\n"))
    }

    ## Do all non-clustering or inv-clustering first
    clusterCellsNA_cKD_Idx <- which(is.na(cellKappaData$clusterCell))
    clusterCellsNA_Names <- cellKappaData$cell[clusterCellsNA_cKD_Idx]
    clusterCellsNA_Vars <-
      cellVarData[cellVarData[, "stage"] == stage, clusterCellsNA_Names]

    nonClusterCells_cKD_Idx <- clusterCellsNA_cKD_Idx[clusterCellsNA_Vars == 0]
    invClusterCells_cKD_Idx <- clusterCellsNA_cKD_Idx[clusterCellsNA_Vars > 0]

    ## Non-clustering
    if (length(nonClusterCells_cKD_Idx) > 0) {
      nonClusterCells_data <-
        data.frame(
          "cell" = cellKappaData[nonClusterCells_cKD_Idx, "cell"],
          "kappa" = cellKappaData[nonClusterCells_cKD_Idx, "kappa"]
        )

      data_stages[[stageIdx]] <-
        .generateCSRPatterns(
          stageName = stage,
          peoplePerStage = peoplePerStage,
          imagesPerPerson = imagesPerPerson,
          kappas = nonClusterCells_data$kappa,
          cellTypes = nonClusterCells_data$cell,
          kappaSep = TRUE, imageAdj = imageAdj, personAdj = personAdj
        )
    }
    ## Inv-clustering
    if (length(invClusterCells_cKD_Idx) > 0) {
      invClusterCells_data <-
        data.frame(
          "cell" = cellKappaData[invClusterCells_cKD_Idx, "cell"],
          "kappa" = cellKappaData[invClusterCells_cKD_Idx, "kappa"],
          "var" = NA
        )
      invClusterCells_data$var <-
        as.numeric(cellVarData[cellVarData[, "stage"] == stage, invClusterCells_data$cell])

      data_stages[[stageIdx]] <- rbind(
        data_stages[[stageIdx]],
        .generateInvClusterPatterns(
          stageName = stage,
          peoplePerStage = peoplePerStage,
          imagesPerPerson = imagesPerPerson,
          kappas = invClusterCells_data$kappa,
          cellTypes = invClusterCells_data$cell,
          cellVars = invClusterCells_data$var,
          kappaSep = TRUE, imageAdj = imageAdj,
          personAdj = personAdj
        )
      )
    }

    ## Recursively plot clusters
    completeCells <- clusterCellsNA_Names # Record completed cell generation
    while (length(completeCells) != nrow(cellKappaData)) {
      # See all cells that cluster around newly added (not previously added)
      nextCell_cKD_Idx <- which(cellKappaData$clusterCell %in% completeCells &
        !(cellKappaData$cell %in% completeCells))
      nextCell_cKD <- cellKappaData[nextCell_cKD_Idx, ]
      nextCell_Vars <-
        cellVarData[cellVarData[, "stage"] == stage, nextCell_cKD$cell]

      if (nrow(nextCell_cKD) == 0) {
        stop("Error: There is an impossibility in cell placement")
      }

      data_stages[[stageIdx]] <- rbind(
        data_stages[[stageIdx]],
        .clusterAroundCells(
          clusterData = data_stages[[stageIdx]][
            data_stages[[stageIdx]]$cellType %in% unique(nextCell_cKD$clusterCell),
          ],
          cellVarData = as.numeric(nextCell_Vars),
          stageName = stage,
          cells = nextCell_cKD$cell,
          clusterCells = nextCell_cKD$clusterCell,
          kappas = nextCell_cKD$kappa,
          minPts = 1
        )
      )

      # Record generation
      completeCells <- c(completeCells, nextCell_cKD$cell)
    }
  }

  ## Organize and return
  data_ret <- .convertList2Dataframe(data_stages, typeBind = "row")[, c(6, 1:3, 5, 4)]
  data_ret$Stage <- as.character(data_ret$Stage)
  data_ret
}


#' Simulate Meta Variables
#'
#' This function simulates meta-variables to append to pca data.
#'
#' Notes: runif may induce useless information so don't recommend correlating it
#'
#' @param pcaData Data.frame with the outcome, unit and principle components of
#'     computed K functions.
#' @param outcome (Optional) Column title for the outcome in the pcaData.
#' @param metaInfo (Optional) Data.frame indicating the metavariables to
#'     include.
#'
#' The data.frame has a var column, rdist column, and columns for each outcome.
#'     The var column names the meta-variables, rdist indicates the distribution
#'     (the options as runif, rbinom, and rnorm), and each outcome. column gives
#'     the expected value (must be positive) of the random variables for each
#'     outcome.
#'
#' In order to allow designation of the expected values, the following rules are
#'     imposed on each distribution:
#'     \itemize{
#'         \item runif: a=0, so b is modified,
#'         \item rbinom: n=1, so this defines the probability
#'         \item runif: the standard deviation is set to 1
#'     }
#'
#' @return Data.frame with the outcome, unit, principle components of
#'     computed K functions, and the meta-variables. pcaData with appended
#'     appended columns at the end.
#' @export
#'
#' @examples
#' data <- simulatePP(
#'   cellVarData = data.frame(
#'     "stage" = c(0, 1, 2),
#'     "A" = c(0, 0, 0),
#'     "B" = c(1 / 100, 1 / 500, 1 / 1000)
#'   ),
#'   cellKappaData = data.frame(
#'     "cell" = c("A", "B"),
#'     "clusterCell" = c(NA, "A"),
#'     "kappa" = c(10, 3)
#'   ),
#'   peoplePerStage = 5,
#'   imagesPerPerson = 1
#' )
#' pcaData <- getPCAData(
#'   data = data, repeatedUniqueId = "Image",
#'   xRange = c(0, 1), yRange = c(0, 1)
#' )
#' pcaMeta <- simulateMeta(pcaData)
simulateMeta <- function(pcaData,
                         outcome = colnames(pcaData)[1],
                         metaInfo = data.frame(
                           "var" = c(
                             "randUnif", "randBin", "rNorm",
                             "corrUnif", "corrBin", "corrNorm"
                           ),
                           "rdist" = c(
                             "runif", "rbinom", "rnorm",
                             "runif", "rbinom", "rnorm"
                           ),
                           "Stage_0" = c(
                             "0.5", "0.5", "1",
                             "0.5", "0.6", "1"
                           ),
                           "Stage_1" = c(
                             "0.5", "0.5", "1",
                             "0.75", "0.65", "1.5"
                           ),
                           "Stage_2" = c(
                             "0.5", "0.5", "1",
                             "0.95", "0.75", "1.5"
                           )
                         )) {
  # Setup Outcome_df
  outcomes_df <- data.frame(
    "outcome" = unique(pcaData[[outcome]]),
    "MetaInfoCol" = NA
  )
  for (i in 1:nrow(outcomes_df)) {
    outcomes_df$MetaInfoCol[i] <-
      which(colnames(metaInfo) == paste0("Stage_", outcomes_df[i, "outcome"]))
  }


  for (i in 1:nrow(metaInfo)) {
    # Setup
    pcaData[[metaInfo[i, 1]]] <- NA
    typeDist <- metaInfo[i, 2]

    # Stages
    for (j in 1:nrow(outcomes_df)) {
      EX <- as.numeric(metaInfo[i, outcomes_df$MetaInfoCol[j]])
      n <- nrow(pcaData[pcaData[, outcome] == outcomes_df[j, "outcome"], ])
      if (EX <= 0) stop("Error: EX must be positive")
      if (typeDist == "runif") {
        # 1/2(b-a)=EX -> b = 2*EX
        parseString <- paste0(
          typeDist, "(", n,
          ", min=0, max=", 2 * EX, ")"
        )
      } else if (typeDist == "rnorm") {
        # mu = EX
        parseString <- paste0(typeDist, "(", n, ", mean=", EX, ")")
      } else if (typeDist == "rbinom") {
        # EX=p
        parseString <- paste0(
          typeDist, "(", n,
          ", size=1, prob=", EX, ")"
        )
      } else {
        stop("Error: Only runif, rnorm, and rbinom accepted.")
      }

      pcaData[pcaData[, outcome] == outcomes_df[j, "outcome"], metaInfo[i, 1]] <-
        eval(parse(text = parseString))
    }
  }

  pcaData
}


#' Generate CSR Point Patterns
#'
#' This (internal) function is used to generate CSR point patterns. It calls
#'     .generateCSRData while also adding necessary data and allowing vectors to
#'     be entered into the function.
#'
#' See usage in simulatePP.
#'
#' Upcoming: Try to eliminate or indicate when kappaSep is used.
#'
#' @param stageName String inidicating the outcome name that should be given.
#' @param peoplePerStage Numeric indicating the number of units per
#'     outcome.
#' @param imagesPerPerson Numeric indicating the number of repeated
#'     measures.
#' @param kappas Numeric vector directing the number of agents per type per
#'     image.
#' @param cellTypes Vector with the agent types of interest being generated.
#'     Corresponds to the kappas.
#' @param kappaSep (Optional) Boolean. Current use understood but not used. May
#'     remove.
#'
#' @return Data.frame containing each point in the defined CSR patterns.
#'
#' The data.frame has columns for outcome, x coordinate, y coordinate, agent
#'     type, unit, and unique repeated measure id.
#' @noRd
.generateCSRPatterns <- function(stageName,
                                 peoplePerStage,
                                 imagesPerPerson,
                                 kappas,
                                 cellTypes,
                                 kappaSep = TRUE,
                                 imageAdj = 0,
                                 personAdj = 0) {
  data <- NULL

  for (personCt in 1:peoplePerStage) {
    for (imageCt in 1:imagesPerPerson) {
      data_tmp <- list()
      for (cell in 1:length(cellTypes)) {
        if (cell == 1 || kappaSep) {
          kapVal <- kappas[cell]
        } else {
          ## Used so that there are more cells when there is no clustering
          #   TODO: See if ever used
          kapVal <- nrow(data_tmp[[cell - 1]]) * kappas[cell]
        }
        data_tmp[[cell]] <- .generateCSRData(
          xRange = c(0, 1), yRange = c(0, 1),
          kappa = kapVal,
          cellType = cellTypes[cell]
        )
      }
      # Clean data (with correct info)
      data_tmp_df <- .convertList2Dataframe(data_tmp, typeBind = "row")
      data_tmp_df$Image <- imageCt + (personCt - 1) * imagesPerPerson + imageAdj
      data_tmp_df$Person <- paste0("p", personCt + personAdj)
      data_tmp_df$Stage <- stageName

      # Save to full data
      data <- rbind(data, data_tmp_df)
    }
  }

  data
}


#' Generate Self/Invisible Clustering Patterns
#'
#' This (internal) function creates self or invisible clustering patterns. It
#'     relies on .generateCSRPatterns and .clusterAroundCells, dropping
#'     unnecessary data.
#'
#' See usage in simulatePP.
#'
#' Upcoming: Determine (generally) how many clusters and how big of clusters to
#'     build.
#'
#' @param stageName String inidicating the outcome name that should be given.
#' @param peoplePerStage Numeric indicating the number of units per
#'     outcome.
#' @param imagesPerPerson Numeric indicating the number of repeated
#'     measures.
#' @param kappas Numeric vector directing the number of agents per type per
#'     image.
#' @param cellTypes Vector with the agent types of interest being generated.
#'     Corresponds to the kappas.
#' @param cellVars Vector with the variables for each agent type being
#'     generated. Corresponds to the kappas and cellTypes.
#' @param kappaSep (Optional) Boolean. Current use understood but not used. May
#'     remove.
#' @param imageAdj (Optional) Numeric to add to image count. Used when images
#'     have already been generated in the resulting dataset and this will be
#'     appended.
#' @return Data.frame containing each point in the defined clustered patterns.
#'     Note, this does not include the data that it clustered around, hence
#'     creating self/invisible clustering.
#'
#' The data.frame has columns for outcome, x coordinate, y coordinate, agent
#'     type, unit, and unique repeated measure id.
#' @noRd
.generateInvClusterPatterns <- function(stageName,
                                        peoplePerStage,
                                        imagesPerPerson,
                                        kappas,
                                        cellTypes,
                                        cellVars,
                                        kappaSep = FALSE,
                                        imageAdj = 0,
                                        personAdj = 0) {
  data <- NULL
  # Consider ways to decide how many clusters
  clusterKappas <- rep(3, length(kappas))
  dataKappas <- kappas / clusterKappas

  # Generate data to cluster around
  cluster_data <- .generateCSRPatterns(stageName,
    peoplePerStage,
    imagesPerPerson,
    clusterKappas,
    paste0(cellTypes, "_Cluster"),
    kappaSep = TRUE,
    imageAdj = imageAdj,
    personAdj = personAdj
  )
  # Cluster data
  data <- .clusterAroundCells(
    clusterData = cluster_data,
    cellVarData = cellVars,
    stageName = stageName,
    cells = cellTypes,
    clusterCells = paste0(cellTypes, "_Cluster"),
    kappas = dataKappas,
    minPts = 1
  )
  data
}


#' Cluster Point Around Cells
#'
#' This (internal) function generates cells that cluster around the given
#'     points.
#'
#' See usage in simulatePP.
#'
#' Upcoming: Convert from cells to agents. Allow different clusterData column
#'     names.
#'
#' @param clusterData Data.frame with columns x, y, cellType, Image, Person, and
#'     Stage. These names are required.
#' @param stageName String inidicating the outcome name that should be given.
#' @param cells Vector indicating the agents that are clustering.
#' @param clusterCells Vector indicating the agents in which each agent in cells
#'     will cluster around. These should be found in clusterData column
#'     cellType. Corresponds to cells.
#' @param kappas Numeric vector directing the number of agents per type per
#'     image. Corresponds to cells.
#' @param minPts (Optional) Numeric indicating the minimum number of points.
#'     Although potentially distributive to the process, the default is 1 to
#'     ensure the cells are placed and any cells that cluster around these can
#'     be placed.
#'
#' @return Data.frame containing each point in the defined clustered patterns.
#'     Note, this does not include the data that it clustered around.
#'
#' The data.frame has columns for outcome, x coordinate, y coordinate, agent
#'     type, unit, and unique repeated measure id.
#' @noRd
.clusterAroundCells <- function(clusterData, cellVarData,
                                stageName,
                                cells, clusterCells, kappas,
                                minPts = 1) {
  newData <- data.frame()
  # Go through each cell
  for (i in 1:length(cells)) {
    clusterCellData <- clusterData[clusterData$cellType == clusterCells[i], ]

    # Go through each cell and develop clusters
    for (j in 1:nrow(clusterCellData)) {
      data_pts <- .placeClusteredPts(
        currXY = as.numeric(clusterCellData[j, c("x", "y")]),
        cell = cells[i],
        numPts = stats::rpois(1, kappas[i]),
        varValue = cellVarData[i]
      )
      if (!is.null(data_pts)) {
        data_pts$Image <- clusterCellData[j, "Image"]
        data_pts$Person <- clusterCellData[j, "Person"]
        data_pts$Stage <- clusterCellData[j, "Stage"]

        newData <- rbind(newData, data_pts)
      }
    }

    # Require minPts
    #     Note, this can change your distribution if not thought about!
    while (nrow(newData[newData$cellType == cells[i], ]) < minPts) {
      # Select a point from previous iteration
      preItrPt <- sample(nrow(clusterCellData), 1)
      # Fill with enough pts
      numPts <- stats::rpois(1, kappas[i])
      if (numPts + nrow(newData[newData$cellType == cells[i], ] < minPts)) {
        numPts <- minPts - nrow(newData[newData$cellType == cells[i], ])
      }

      data_pts <- .placeClusteredPts(
        currXY = as.numeric(clusterCellData[j, c("x", "y")]),
        cell = cells[i],
        numPts = stats::rpois(1, kappas[i]),
        varValue = cellVarData[i]
      )
      data_pts$Image <- clusterCellData[j, "Image"]
      data_pts$Person <- clusterCellData[j, "Person"]
      data_pts$Stage <- clusterCellData[j, "Stage"]

      newData <- rbind(newData, data_pts)
    }
  }

  newData
}


#' Place Clustered Points
#'
#' This (internal) function places points following a normal distribution around
#'     some given corridinates.
#'
#' See usage in simulatePP.
#'
#' Upcoming: Generalize naming
#'
#' @param currXY Vector of two numerics, the first relating the x coordinate
#'     and the second to the y coordinate. These turn into the mean of the
#'     normal distribution
#' @param cell String indicating the agent type of the currently placed agents.
#' @param numPts Numeric indicating the number of agents to place.
#' @param varValue Numeric giving the variance of the normal distribution for
#'     placement of the agents.
#' @param xRange (Optional) Vector of two values indicating the x range of the
#'     region. Default is c(0,1).
#' @param yRange (Optional) Vector of two values indicating the x range of the
#'     region. Default is c(0,1).
#'
#' @return Data.frame with placed cells.
#'
#' The data.frame has 3 columns, x, y, and cellType.
#' @noRd
.placeClusteredPts <- function(currXY, cell, numPts, varValue,
                               xRange = c(0, 1), yRange = c(0, 1)) {
  if (numPts <= 0) {
    return()
  }

  data_ret <- data.frame("x" = rep(NA, numPts), "y" = NA, "cellType" = cell)

  compPts <- 0
  while (compPts < numPts) {
    data_ret[compPts + 1, "x"] <-
      stats::rnorm(1, mean = currXY[1], sd = sqrt(varValue))
    data_ret[compPts + 1, "y"] <-
      stats::rnorm(1, mean = currXY[2], sd = sqrt(varValue))

    # Ensure its in boundaries
    if (data_ret[compPts + 1, "x"] >= xRange[1] &
      data_ret[compPts + 1, "x"] <= xRange[2] &
      data_ret[compPts + 1, "y"] >= yRange[1] &
      data_ret[compPts + 1, "y"] <= yRange[2]) {
      compPts <- compPts + 1
    }
  }

  data_ret
}
