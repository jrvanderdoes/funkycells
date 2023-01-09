
simulatePP_1 <- function(cellVarData=
                         data.frame('stage'=c('Stage_0','Stage_1'),
                                    'A'=c(0,0),
                                    'B'=c(1/100,1/100),
                                    'C'=c(1/500,1/250)),
                         clusterKappaData=data.frame(
                           'cell'=c('A','B','C'),
                           'Stage_0'=c(NA,NA,NA),
                           'Stage_1'=c(NA,'A',NA),
                           'kappa'=c(NA,5,NA)),
                         indKappaData=data.frame(
                           'cell'=c('A','B','C'),
                           'Stage_0'=c(NA,NA,NA),
                           'Stage_1'=c(NA,'A',NA),
                           'kappa'=c(20,5,4)),
                       peoplePerStage=20,
                       imagesPerPerson=5,
                       reduceEdge=0.025,
                       silent=F ){
  ## Setup
  data_stages <- list()
  data_stages1 <- list()
  # Go through each stage
  for(stageIdx in 1:nrow(cellVarData)){
    # General Vars
    stage <- cellVarData$stage[stageIdx] # Current Stage
    imageAdj <- (stageIdx-1)*(peoplePerStage*imagesPerPerson) # Adjustment for images due to stage
    personAdj <- (stageIdx-1)*(peoplePerStage) # Adjustment for person due to stage

    if(!silent)
      cat(paste0('Stage: ',stage, ' (',stageIdx,'/',nrow(cellVarData),')\n'))

    ## Do all non-clustering or inv-clustering first
    clusterCellsNA_cKD_Idx <- which(is.na(cellKappaData$clusterCell))
    clusterCellsNA_Names <- cellKappaData$cell[clusterCellsNA_cKD_Idx]
    clusterCellsNA_Vars <-
      cellVarData[cellVarData[,'stage']==stage,clusterCellsNA_Names]

    nonClusterCells_cKD_Idx <- clusterCellsNA_cKD_Idx[clusterCellsNA_Vars==0]
    invClusterCells_cKD_Idx <- clusterCellsNA_cKD_Idx[clusterCellsNA_Vars>0]

    ## Non-clustering
    if(length(nonClusterCells_cKD_Idx)>0){
      nonClusterCells_data <-
        data.frame('cell'=cellKappaData[nonClusterCells_cKD_Idx,'cell'],
                   'kappa'=cellKappaData[nonClusterCells_cKD_Idx,'kappa'])

      data_stages[[stageIdx]] <-
        .generateCSRPatterns(stageName=stage,
                             reduceEdge=reduceEdge,
                             peoplePerStage=peoplePerStage,
                             imagesPerPerson=imagesPerPerson,
                             kappas=nonClusterCells_data$kappa,
                             cellTypes=nonClusterCells_data$cell,
                             kappaSep=T, imageAdj=imageAdj, personAdj=personAdj)

    }

    ## Inv-clustering
    if(length(invClusterCells_cKD_Idx)>0){
      invClusterCells_data <-
        data.frame('cell'=cellKappaData[invClusterCells_cKD_Idx,'cell'],
                   'kappa'=cellKappaData[invClusterCells_cKD_Idx,'kappa'],
                   'var'=NA)
      invClusterCells_data$var <-
        cellVarData[cellVarData[,'stage']==stage,invClusterCells_data$cell]

      data_stages[[stageIdx]] <- rbind(data_stages[[stageIdx]],
                                       .generateInvClusterPatterns(
                                         stageName=stage,
                                         reduceEdge=reduceEdge,
                                         peoplePerStage=peoplePerStage,
                                         imagesPerPerson=imagesPerPerson,
                                         kappas=invClusterCells_data$kappa,
                                         cellTypes=invClusterCells_data$cell,
                                         cellVars=invClusterCells_data$var,
                                         kappaSep=T, imageAdj=imageAdj,
                                         personAdj=personAdj)
      )
    }

    ## Recursively plot clusters
    completeCells <- clusterCellsNA_Names # Record completed cell generation
    while(length(completeCells) != nrow(cellKappaData)){

      # See all cells that cluster around newly added (not previously added)
      nextCell_cKD_Idx <- which(cellKappaData$clusterCell %in% completeCells &
                                  !(cellKappaData$cell %in% completeCells))
      nextCell_cKD <- cellKappaData[nextCell_cKD_Idx,]
      nextCell_Vars <-
        cellVarData[cellVarData[,'stage']==stage,nextCell_cKD$cell]

      if(nrow(nextCell_cKD)==0)
        stop('Error: There is an impossibility in cell placement')

      data_stages[[stageIdx]] <- rbind(data_stages[[stageIdx]],
                                       .clusterAroundCells(
                                         clusterData=data_stages[[stageIdx]][
                                           data_stages[[stageIdx]]$cellType %in% unique(nextCell_cKD$clusterCell),],
                                         cellVarData=as.numeric(nextCell_Vars),
                                         stageName=stage, reduceEdge=reduceEdge,
                                         cells=nextCell_cKD$cell,
                                         clusterCells=nextCell_cKD$clusterCell,
                                         kappas=nextCell_cKD$kappa,
                                         minPts=1)
      )

      # Add Seperate Clustering
      if(makeInteractionUniDirectional){

        uniDirection_data <-
          data.frame('cell'=cellKappaData[nextCell_cKD_Idx,'cell'],
                     'kappa'=cellKappaData[nextCell_cKD_Idx,'kappa'],
                     'var'=NA)

        uniDirection_data$var <-
          as.numeric(cellVarData[cellVarData[,'stage']==stage,
                                 uniDirection_data$cell])
        data_stages[[stageIdx]] <- rbind(data_stages[[stageIdx]],
                                         .generateInvClusterPatterns(
                                           stageName=stage,
                                           reduceEdge=reduceEdge,
                                           peoplePerStage=peoplePerStage,
                                           imagesPerPerson=imagesPerPerson,
                                           kappas=uniDirection_data$kappa,
                                           cellTypes=uniDirection_data$cell,
                                           cellVars=uniDirection_data$var,
                                           kappaSep=T, imageAdj=imageAdj,
                                           personAdj=personAdj)
        )
      }

      # Record generation
      completeCells <- c(completeCells, nextCell_cKD$cell)
    }

  }

  ## Organize and return
  data_ret <- .convertList2Dataframe(data_stages, typeBind = 'row')[,c(6,1:3,5,4)]
  data_ret$Stage <- as.character(data_ret$Stage)
  data_ret
}

dat <- simulatePP_1(cellVarData=
                      data.frame('stage'=c(0,1),
                                 'A'=c(0,0),
                                 'B'=c(1/100,1/500)),
                    cellKappaData=data.frame(
                      'cell'=c('A','B'),
                      'clusterCell'=c(NA,'A'),
                      'kappa'=c(20,5)),
                    makeInteractionUniDirectional=T)
pcaDat <- getPCAData(dat,repeatedUniqueId='Image',outcome = 'Stage',unit = 'Person',
                          xRange = c(0,1),  yRange = c(0,1), silent=F)
rfcv <- computeRandomForest_CVPC(data=pcaDat,repeatedId='Image',
                                 outcome = "Stage",unit = 'Person',
                                 cellData=dat)
