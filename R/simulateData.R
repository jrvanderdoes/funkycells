#' Simulate a Point Process
#'
#' This function simulates a point pattern with optional clustering (visible
#'     and invisible). Multiple outcomes, units, and replicates are
#'     possible, e.g. a 3 stage disease (outcomes) over 20 people (units) with
#'     3 images each (replicates).
#'
#' @param agentVarData (Optional) Data.frame describing variances with each
#'  agent type.
#'
#' The data.frame has a outcome column and a named column for each agent type.
#'     Currently, these names are mandatory.
#' @param agentKappaData (Optional) Data.frame describing agent interactions.
#'
#' The data.frame has a agent column giving agent names (matching agentVarData),
#'     a clusterAgent column indicating which agent the agent clusters (put NA
#'     if the agent doesn't cluster or clusters a hidden agent / self-clusters),
#'     and a kappa column directing the number of agents of per replicate.
#' @param unitsPerOutcome (Optional) Numeric indicating the number of units per
#'     outcome.
#' @param replicatesPerUnit (Optional) Numeric indicating the number of
#'  replicates, or repeated measures, per unit.
#' @param silent (Optional) Boolean indicating if progress output should be
#'     printed.
#'
#' @return Data.frame containing each point the defined patterns.
#'
#' The data.frame has columns for outcome, x coordinate, y coordinate, agent
#'     type, unit, and replicate id.
#' @export
#'
#' @examples
#' data <- simulatePP(
#'   agentVarData = data.frame(
#'     "outcome" = c(0, 1),
#'     "A" = c(0, 0),
#'     "B" = c(1 / 100, 1 / 500),
#'     "C" = c(1 / 500, 1 / 250),
#'     "D" = c(1 / 100, 1 / 100),
#'     "E" = c(1 / 500, 1 / 500)
#'   ),
#'   agentKappaData = data.frame(
#'     "agent" = c("A", "B", "C", "D", "E"),
#'     "clusterAgent" = c(NA, "A", "B", "C", NA),
#'     "kappa" = c(10, 3, 2, 1, 8)
#'   ),
#'   unitsPerOutcome = 4,
#'   replicatesPerUnit = 1
#' )
simulatePP <- function(agentVarData =
                         data.frame(
                           "outcome" = c(0, 1, 2),
                           "A" = c(0, 0, 0),
                           "B" = c(1 / 100, 1 / 500, 1 / 500),
                           "C" = c(1 / 500, 1 / 250, 1 / 100),
                           "D" = c(1 / 100, 1 / 100, 1 / 100),
                           "E" = c(1 / 500, 1 / 500, 1 / 500),
                           "F" = c(1 / 250, 1 / 250, 1 / 250)
                         ),
                       agentKappaData = data.frame(
                         "agent" = c("A", "B", "C", "D", "E", "F"),
                         "clusterAgent" = c(NA, "A", "B", "C", NA, "A"),
                         "kappa" = c(20, 5, 4, 2, 15, 5)
                       ),
                       unitsPerOutcome = 20,
                       replicatesPerUnit = 5,
                       silent = FALSE) {
  ## Setup
  data_outcomes <- list()
  data_outcomes1 <- list()
  # Go through each outcome
  for (outcomeIdx in 1:nrow(agentVarData)) {
    # General Vars
    outcome <- agentVarData$outcome[outcomeIdx] # Current Stage
    replicateAdj <- (outcomeIdx - 1) * (unitsPerOutcome * replicatesPerUnit) # Adjustment for replicates due to stage
    unitAdj <- (outcomeIdx - 1) * (unitsPerOutcome) # Adjustment for unit due to stage
    data_outcomes[[outcomeIdx]] <- data.frame() # Make this so no errors
    if (!silent) {
      cat(paste0(
        "Outcome: ", outcome, " (", outcomeIdx, "/",
        nrow(agentVarData), ")\n"
      ))
    }

    ## Do all non-clustering or inv-clustering first
    clusterAgentsNA_cKD_Idx <- which(is.na(agentKappaData$clusterAgent))
    clusterAgentsNA_Names <- agentKappaData$agent[clusterAgentsNA_cKD_Idx]
    clusterAgentsNA_Vars <-
      agentVarData[agentVarData[, "outcome"] == outcome, clusterAgentsNA_Names]

    nonClusterAgents_cKD_Idx <- clusterAgentsNA_cKD_Idx[clusterAgentsNA_Vars == 0]
    invClusterAgents_cKD_Idx <- clusterAgentsNA_cKD_Idx[clusterAgentsNA_Vars > 0]

    ## Non-clustering
    if (length(nonClusterAgents_cKD_Idx) > 0) {
      nonClusterAgents_data <-
        data.frame(
          "agent" = agentKappaData[nonClusterAgents_cKD_Idx, "agent"],
          "kappa" = agentKappaData[nonClusterAgents_cKD_Idx, "kappa"]
        )

      data_outcomes[[outcomeIdx]] <-
        .generateCSRPatterns(
          outcomeName = outcome,
          unitsPerOutcome = unitsPerOutcome,
          replicatesPerUnit = replicatesPerUnit,
          kappas = nonClusterAgents_data$kappa,
          types = nonClusterAgents_data$agent,
          replicateAdj = replicateAdj,
          unitAdj = unitAdj
        )
    }
    ## Inv-clustering
    if (length(invClusterAgents_cKD_Idx) > 0) {
      invClusterAgents_data <-
        data.frame(
          "agent" = agentKappaData[invClusterAgents_cKD_Idx, "agent"],
          "kappa" = agentKappaData[invClusterAgents_cKD_Idx, "kappa"],
          "var" = NA
        )
      invClusterAgents_data$var <-
        as.numeric(agentVarData[
          agentVarData[, "outcome"] == outcome,
          invClusterAgents_data$agent
        ])

      data_outcomes[[outcomeIdx]] <- rbind(
        data_outcomes[[outcomeIdx]],
        .generateInvClusterPatterns(
          outcomeName = outcome,
          unitsPerOutcome = unitsPerOutcome,
          replicatesPerUnit = replicatesPerUnit,
          kappas = invClusterAgents_data$kappa,
          types = invClusterAgents_data$agent,
          agentVars = invClusterAgents_data$var,
          replicateAdj = replicateAdj,
          unitAdj = unitAdj
        )
      )
    }

    ## Recursively plot clusters
    completeAgents <- clusterAgentsNA_Names # Record completed agent generation
    while (length(completeAgents) != nrow(agentKappaData)) {
      # See all types that cluster around newly added (not previously added)
      nextAgent_cKD_Idx <- which(agentKappaData$clusterAgent %in% completeAgents &
        !(agentKappaData$agent %in% completeAgents))
      nextAgent_cKD <- agentKappaData[nextAgent_cKD_Idx, ]
      nextAgent_Vars <-
        agentVarData[agentVarData[, "outcome"] == outcome, nextAgent_cKD$agent]

      if (nrow(nextAgent_cKD) == 0) {
        stop("Error: There is an impossibility in agent placement")
      }

      data_outcomes[[outcomeIdx]] <- rbind(
        data_outcomes[[outcomeIdx]],
        .clusterAroundAgents(
          clusterData = data_outcomes[[outcomeIdx]][
            data_outcomes[[outcomeIdx]]$type %in% unique(nextAgent_cKD$clusterAgent),
          ],
          agentVarData = as.numeric(nextAgent_Vars),
          outcomeName = outcome,
          types = nextAgent_cKD$agent,
          clusterAgents = nextAgent_cKD$clusterAgent,
          kappas = nextAgent_cKD$kappa,
          minPts = 1
        )
      )

      # Record generation
      completeAgents <- c(completeAgents, nextAgent_cKD$agent)
    }
  }

  ## Organize and return
  data_ret <- do.call("rbind", data_outcomes)[, c(6, 1:3, 5, 4)]
  data_ret$outcome <- as.character(data_ret$outcome)

  data_ret
}


#' Simulate Meta Variables
#'
#' This function simulates meta-variables with varying distributions to append
#'  to some data.
#'
#' Notes: runif may induce useless information so don't recommend correlating it
#'
#' @param data Data.frame with the outcome and unit. Typically this also
#'  includes PCA data as it is run after computing the principle components (see
#'  examples).
#' @param outcome (Optional) String for column title of the data's outcome.
#'  Default is the first column.
#' @param metaInfo (Optional) Data.frame indicating the meta-variables (and
#'  properties) to generate. Default has some examples of possible options.
#'
#' The data.frame has a var column, rdist column, and columns for each outcome.
#'     The var column names the meta-variables, rdist indicates the distribution
#'     (options are runif, rbinom, and rnorm), and the outcome columns indicate
#'     mean of the variable for that outcome.
#'
#' In order to allow designation of the expected values, the following rules are
#'     imposed on each distribution:
#'     \itemize{
#'         \item runif: a=0, so b is modified,
#'         \item rbinom: n=1, so this defines the probability
#'         \item runif: variance is set to 1
#'     }
#'
#' @return Data.frame of the original data with meta-variables appended (as
#' columns) at the end.
#' @export
#'
#' @examples
#' data <- simulatePP(
#'   agentVarData = data.frame(
#'     "outcome" = c(0, 1, 2),
#'     "A" = c(0, 0, 0),
#'     "B" = c(1 / 100, 1 / 500, 1 / 1000)
#'   ),
#'   agentKappaData = data.frame(
#'     "agent" = c("A", "B"),
#'     "clusterAgent" = c(NA, "A"),
#'     "kappa" = c(10, 3)
#'   ),
#'   unitsPerOutcome = 5,
#'   replicatesPerUnit = 1
#' )
#' pcaData <- getKsPCAData(
#'   data = data, replicate = "replicate",
#'   xRange = c(0, 1), yRange = c(0, 1)
#' )
#' pcaMeta <- simulateMeta(pcaData)
#'
#' ## Another simple example
#' data <- simulateMeta(
#'   data.frame("outcome" = c(0, 0, 0, 1, 1, 1), "unit" = 1:6)
#' )
simulateMeta <- function(data,
                         outcome = colnames(data)[1],
                         metaInfo = data.frame(
                           "var" = c(
                             "randUnif", "randBin", "rNorm",
                             "corrUnif", "corrBin", "corrNorm"
                           ),
                           "rdist" = c(
                             "runif", "rbinom", "rnorm",
                             "runif", "rbinom", "rnorm"
                           ),
                           "outcome_0" = c(
                             "0.5", "0.5", "1",
                             "0.5", "0.6", "1"
                           ),
                           "outcome_1" = c(
                             "0.5", "0.5", "1",
                             "0.75", "0.65", "1.5"
                           ),
                           "outcome_2" = c(
                             "0.5", "0.5", "1",
                             "0.95", "0.75", "1.5"
                           )
                         )) {
  # Setup Outcome_df
  outcomes_df <- data.frame(
    "outcome" = unique(data[[outcome]]),
    "MetaInfoCol" = NA
  )
  for (i in 1:nrow(outcomes_df)) {
    outcomes_df$MetaInfoCol[i] <-
      which(colnames(metaInfo) == paste0("outcome_", outcomes_df[i, "outcome"]))
  }


  for (i in 1:nrow(metaInfo)) {
    # Setup
    data[[metaInfo[i, 1]]] <- NA
    typeDist <- metaInfo[i, 2]

    # Outcomes
    for (j in 1:nrow(outcomes_df)) {
      EX <- as.numeric(metaInfo[i, outcomes_df$MetaInfoCol[j]])
      signEX <- sign(EX)
      EX <- abs(EX)
      n <- nrow(data[data[, outcome] == outcomes_df[j, "outcome"], ])
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

      data[data[, outcome] == outcomes_df[j, "outcome"], metaInfo[i, 1]] <-
        eval(parse(text = parseString)) * signEX
    }
  }

  data
}


#' Generate CSR Point Patterns
#'
#' This (internal) function is used to generate CSR point patterns. It calls
#'     .generateCSRData while also adding necessary data and allowing vectors to
#'     be entered into the function.
#'
#' See usage in simulatePP.
#'
#' @param outcomeName String indicating the outcome name that should be given.
#' @param unitsPerOutcome Numeric indicating the number of units per
#'     outcome.
#' @param replicatesPerUnit Numeric indicating the number of replicates, or
#'  repeated measures, per unit.
#' @param kappas Numeric vector directing the number of agents per type per
#'     replicate.
#' @param types Character vector with the agent types of interest being
#'  generated. Corresponds to the kappas.
#' @param replicateAdj (Optional) Numeric that increments replicate name. Used
#'  when replicates have already been generated in the resulting data set and
#'  this will be appended.
#' @param unitAdj (Optional) Numeric that increments unit name. Used when replicates
#'     have already been generated in the resulting data set and this will be
#'     appended.
#'
#' @return Data.frame containing each point in the defined CSR patterns.
#'
#' The data.frame has columns for outcome, x coordinate, y coordinate, agent
#'     type, unit, and replicate.
#' @noRd
.generateCSRPatterns <- function(outcomeName,
                                 unitsPerOutcome,
                                 replicatesPerUnit,
                                 kappas,
                                 types,
                                 replicateAdj = 0,
                                 unitAdj = 0) {
  data <- NULL

  for (unitCt in 1:unitsPerOutcome) {
    for (replicateCt in 1:replicatesPerUnit) {
      data_tmp <- list()
      for (agent in 1:length(types)) {
        kapVal <- kappas[agent]
        data_tmp[[agent]] <- .generateCSRData(
          xRange = c(0, 1), yRange = c(0, 1),
          kappa = kapVal,
          type = types[agent]
        )
      }
      # Clean data (with correct info)
      data_tmp_df <- do.call("rbind", data_tmp)
      data_tmp_df$replicate <- replicateCt + (unitCt - 1) *
        replicatesPerUnit + replicateAdj
      data_tmp_df$unit <- paste0("u", unitCt + unitAdj)
      data_tmp_df$outcome <- outcomeName

      # Add to full data
      data <- rbind(data, data_tmp_df)
    }
  }

  data
}


#' Generate Self/Invisible Clustering Patterns
#'
#' This (internal) function creates self or invisible clustering patterns. It
#'     relies on .generateCSRPatterns and .clusterAroundAgents, dropping
#'     unnecessary data.
#'
#' See usage in simulatePP.
#'
#' Upcoming: Determine (generally) how many clusters and how big of clusters to
#'     build.
#' @inheritParams .generateCSRPatterns
#' @param agentVars Vector with the variables for each agent type being
#'     generated. Corresponds to the kappas and types.
#'
#' @return Data.frame containing each point in the defined clustered patterns.
#'     Note, this does not include the data that it clustered around, hence
#'     creating self/invisible clustering.
#'
#' The data.frame has columns for outcome, x coordinate, y coordinate, agent
#'     type, unit, and replicate.
#' @noRd
.generateInvClusterPatterns <- function(outcomeName,
                                        unitsPerOutcome,
                                        replicatesPerUnit,
                                        kappas,
                                        types,
                                        agentVars,
                                        replicateAdj = 0,
                                        unitAdj = 0) {
  data <- NULL
  # Consider ways to decide how many clusters
  clusterKappas <- rep(3, length(kappas))
  dataKappas <- kappas / clusterKappas

  # Generate data to cluster around
  cluster_data <- .generateCSRPatterns(outcomeName,
    unitsPerOutcome,
    replicatesPerUnit,
    clusterKappas,
    paste0(types, "_Cluster"),
    replicateAdj = replicateAdj,
    unitAdj = unitAdj
  )
  # Cluster data
  data <- .clusterAroundAgents(
    clusterData = cluster_data,
    agentVarData = agentVars,
    outcomeName = outcomeName,
    types = types,
    clusterAgents = paste0(types, "_Cluster"),
    kappas = dataKappas,
    minPts = 1
  )
  data
}


#' Cluster Point Around Types
#'
#' This (internal) function generates types that cluster around the given type.
#'
#' See usage in simulatePP.
#'
#' @inheritParams .generateCSRPatterns
#' @param clusterData Data.frame with columns x, y, type, replicate, unit, and
#'     outcome. These names are required (could be extended, but internal
#'     function, so no rush).
#' @param clusterAgents Vector indicating the agents in which each agent in types
#'     will cluster around. These should be found in clusterData column
#'     type. Corresponds to types.
#' @param minPts (Optional) Numeric indicating the minimum number of points.
#'     Although potentially distributive to the process, the default is 1 to
#'     ensure the types are placed and any types that cluster around these can
#'     be placed.
#'
#' @return Data.frame containing each point in the defined clustered patterns.
#'     Note, this does not include the data that it clustered around.
#'
#' The data.frame has columns for outcome, x coordinate, y coordinate, agent
#'     type, unit, and replicate.
#' @noRd
.clusterAroundAgents <- function(clusterData, agentVarData,
                                 outcomeName,
                                 types, clusterAgents, kappas,
                                 minPts = 1) {
  newData <- data.frame()
  # Go through each agent
  for (i in 1:length(types)) {
    clusterAgentData <- clusterData[clusterData$type == clusterAgents[i], ]

    # Go through each agent and develop clusters
    for (j in 1:nrow(clusterAgentData)) {
      data_pts <- .placeClusteredPts(
        currXY = as.numeric(clusterAgentData[j, c("x", "y")]),
        agent = types[i],
        numPts = stats::rpois(1, kappas[i]),
        varValue = agentVarData[i]
      )
      if (!is.null(data_pts)) {
        data_pts$replicate <- clusterAgentData[j, "replicate"]
        data_pts$unit <- clusterAgentData[j, "unit"]
        data_pts$outcome <- clusterAgentData[j, "outcome"]

        newData <- rbind(newData, data_pts)
      }
    }

    # Require minPts
    #     Note, this can change your distribution if not thought about!
    while (nrow(newData[newData$type == types[i], ]) < minPts) {
      # Select a point from previous iteration
      preItrPt <- sample(nrow(clusterAgentData), 1)
      # Fill with enough pts
      numPts <- stats::rpois(1, kappas[i])
      if (numPts + nrow(newData[newData$type == types[i], ]) < minPts) {
        numPts <- minPts - nrow(newData[newData$type == types[i], ])
      }

      data_pts <- .placeClusteredPts(
        currXY = as.numeric(clusterAgentData[j, c("x", "y")]),
        agent = types[i],
        numPts = stats::rpois(1, kappas[i]),
        varValue = agentVarData[i]
      )
      data_pts$replicate <- clusterAgentData[j, "replicate"]
      data_pts$unit <- clusterAgentData[j, "unit"]
      data_pts$outcome <- clusterAgentData[j, "outcome"]

      newData <- rbind(newData, data_pts)
    }
  }

  newData
}


#' Place Clustered Points
#'
#' This (internal) function places points following a normal distribution around
#'     some given coordinates.
#'
#' See usage in simulatePP.
#'
#' @param currXY Vector of two numerics, the first relating the x coordinate
#'     and the second to the y coordinate. These turn into the mean of the
#'     normal distribution
#' @param agent String indicating the agent type of the currently placed agents.
#' @param numPts Numeric indicating the number of agents to place.
#' @param varValue Numeric giving the variance of the normal distribution for
#'     placement of the agents.
#' @param xRange (Optional) Vector of two values indicating the x range of the
#'     region. Default is c(0,1).
#' @param yRange (Optional) Vector of two values indicating the x range of the
#'     region. Default is c(0,1).
#'
#' @return Data.frame with placed types.
#'
#' The data.frame has 3 columns, x, y, and type.
#' @noRd
.placeClusteredPts <- function(currXY, agent, numPts, varValue,
                               xRange = c(0, 1), yRange = c(0, 1)) {
  if (numPts <= 0) {
    return()
  }

  data_ret <- data.frame("x" = rep(NA, numPts), "y" = NA, "type" = agent)

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
