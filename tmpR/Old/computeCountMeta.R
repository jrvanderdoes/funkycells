

# Data should have outcome, unit, and column for agents or a column for each
#     agent

#' Compute Count Meta
#'
#' This function compute the average count of cellType per repeatedUniqueId.
#'
#' @param data Data.frame with outcome, unit, and potentially any other
#'     variables to use in model. Often this would be run after getPCAData.
#' @param outcome String of the column name in data indicating the outcome or
#'     response.
#' @param unit String of the column name in data indicating a unit or base thing.
#'     Note this unit may have repeated measures. Default is the 4th column.
#' @param agents Vector with strings indicating the possible agents.
#' @param cellData Data.frame with outcome, unit, and column for agents or a
#'     column for each agent.
#' @param repeatedUniqueId (optional) String of the column name in data
#'     indicating the unique ID when using repeated measures. Default is NULL
#'
#' @return Data.frame with outcome, unit, any other variables in data and counts
#'     for agents (per repeatedUniqueId if non-null).
#' @export
#'
#' @examples
#' TNBC <- TNBC
#' TNBC_Meta <- TNBC_Meta
#' agents_df <- expand.grid(c("FoxP3","Lag3","CD4","CD16","CD56","PD1","PD.L1",
#'                            "Ki67","CD138","CD68","CD8","CD3","IDO","CD45RO",
#'                            "CD20","p53","MPO","tumorYN"),c("FoxP3","Lag3",
#'                            "CD4","CD16","CD56","PD1","PD.L1","Ki67","CD138",
#'                            "CD68","CD8","CD3","IDO","CD45RO","CD20","p53",
#'                            "MPO","tumorYN"))
#' dataPCA <- getPCAData(data = TNBC[,-3], unit='Person', agents_df=agents_df,
#'                       rCheckVals = seq(0,50,1))
#' dataCt <- computeCountMeta(data = dataPCA, outcome = 'Class', unit = 'Person',
#'                        agents = c("FoxP3","Lag3","CD4","CD16","CD56","PD1",
#'                                   "PD.L1","Ki67","CD138","CD68","CD8","CD3",
#'                                   "IDO","CD45RO","CD20","p53","MPO","tumorYN"),
#'                        cellData = TNBC[,-c(3:5)])
computeCountMeta <- function(data, outcome, unit, agents, cellData,
                         repeatedUniqueId=NULL){
  # Verify way agents are defined:
  #   1. One column naming the agent
  #   2. A column for each agent, values would be intensities or T/F
  if(ncol(cellData[!(colnames(cellData)%in%c(outcome,unit,repeatedUniqueId))])==1){
    metaNames <- c()
    retData <- data
    # Count repeatedImages (default of 1 will change nothing later)
    imgCts <- rep(1,nrow(data))
    if(!is.null(repeatedUniqueId)){
      repData <- unique(cellData[colnames(cellData) %in%
                                   c(outcome,unit,repeatedUniqueId)])
      divData <- aggregate(as.formula(paste('. ~ ',unit,'+',outcome)),
                           repData, FUN=length)
      imgCts <- merge(retData, divData,
                      by=c(outcome,unit),sort=F)[[repeatedUniqueId]]
    }

    # Go through agents
    for(i in 1:length(agents)){
      ctData <- as.data.frame(table(cellData[cellData[,!(colnames(cellData)%in%
                        c(outcome,unit,repeatedUniqueId))] == agents[i],unit]))
      # Standardize to per image
      ctData$Freq <- ctData$Freq / imgCts
      names(ctData)[names(ctData) == 'Freq'] <- paste0(agents[i],'_Ct')
      metaNames <- c(metaNames, paste0(agents[i],'_Ct'))

      retData <- merge(retData, ctData, by.x=unit, by.y='Var1', sort=F)
    }

  }else{
    aggData <- aggregate(as.formula(paste('. ~ ',unit,'+',outcome)),
                         cellData, FUN=sum)
    # Divide to number per image
    if(!is.null(repeatedUniqueId)){
      repData <- unique(cellData[colnames(cellData) %in%
                                   c(outcome,unit,repeatedUniqueId)])
      divData <- aggregate(as.formula(paste('. ~ ',unit,'+',outcome)),
                           repData, FUN=length)

      # TODO::Do this vector-ly
      for(i in 1:nrow(aggData)){
        aggData[i,!(colnames(aggData) %in% c(outcome,unit))] <-
          aggData[i,!(colnames(aggData) %in% c(outcome,unit))] /
          divData[divData[[outcome]]==aggData[i,outcome] &
                    divData[[unit]]==aggData[i,unit],
                  !(colnames(divData) %in% c(outcome,unit)) ]
      }
    }

    metaNames <-
      paste0(colnames(aggData[!(colnames(aggData)%in%c(outcome,unit))]), '_Ct')
    colnames(aggData) <- c(colnames(aggData)[1:2], metaNames)
    retData <- merge(data, aggData, by=c(outcome, unit), sort=F)
  }

  list('data'=retData, 'metaNames'= metaNames)
}
