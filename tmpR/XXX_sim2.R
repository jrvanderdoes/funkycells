#############################

simulatePP1 <- function(cellVarData=
                         data.frame('stage'=c(0,1,2),
                                    'A'=c(0,0,0),
                                    'B'=c(1/100,1/500,1/500),
                                    'C'=c(1/500,1/250,1/100),
                                    'D'=c(1/100,1/100,1/100),
                                    'E'=c(1/500,1/500,1/500),
                                    'F'=c(1/250,1/250,1/250)),
                       cellKappaData=data.frame(
                         'cell'=c('A','B','C','D','E','F'),
                         'clusterCell'=c(NA,'A','B','C',NA,'A'),
                         'kappa'=c(20,5,4,2,15,5)),
                       peoplePerStage=20,
                       imagesPerPerson=5,
                       reduceEdge=0.025,
                       silent=F){
  ## Setup
  data_stages <- list()
  data_stages1 <- list()
  # Go through each stage
  for(stageIdx in 1:nrow(cellVarData)){
    # General Vars
    stage <- cellVarData$stage[stageIdx] # Current Stage
    imageAdj <- (stageIdx-1)*(peoplePerStage*imagesPerPerson) # Adjustment for images due to stage
    personAdj <- (stageIdx-1)*(peoplePerStage) # Adjustment for person due to stage

    cellKappaData_stage <- cellKappaData[cellKappaData$stage==stage,]

    if(!silent)
      cat(paste0('Stage: ',stage, ' (',stageIdx,'/',nrow(cellVarData),')\n'))

    ## Do all non-clustering or inv-clustering first
    clusterCellsNA_cKD_Idx <-
      which(is.na(cellKappaData_stage$clusterCell))
    clusterCellsNA_Names <- cellKappaData_stage$cell[clusterCellsNA_cKD_Idx]
    clusterCellsNA_Vars <-
      cellVarData[cellVarData[,'stage']==stage,clusterCellsNA_Names]

    nonClusterCells_cKD_Idx <- clusterCellsNA_cKD_Idx[clusterCellsNA_Vars==0]
    invClusterCells_cKD_Idx <- clusterCellsNA_cKD_Idx[clusterCellsNA_Vars>0]

    ## Non-clustering
    if(length(nonClusterCells_cKD_Idx)>0){
      nonClusterCells_data <-
        data.frame('cell'=cellKappaData_stage[nonClusterCells_cKD_Idx,'cell'],
                   'kappa'=cellKappaData_stage[nonClusterCells_cKD_Idx,'kappa'])

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
        data.frame('cell'=cellKappaData_stage[invClusterCells_cKD_Idx,'cell'],
                   'kappa'=cellKappaData_stage[invClusterCells_cKD_Idx,'kappa'],
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
    while(length(completeCells) !=
          nrow(cellKappaData_stage[cellKappaData_stage$stage==stage,])){

      # See all cells that cluster around newly added (not previously added)
      nextCell_cKD_Idx <-
            which(cellKappaData_stage$clusterCell %in% completeCells &
                                !(cellKappaData_stage$cell %in% completeCells))
      nextCell_cKD <- cellKappaData_stage[nextCell_cKD_Idx,]
      nextCell_Vars <-
        cellVarData[cellVarData[,'stage']==stage,nextCell_cKD$cell]

      if(nrow(nextCell_cKD)==0)
        stop('Error: There is an impossibility in cell placement')

      data_stages[[stageIdx]] <-
          rbind(data_stages[[stageIdx]],
                .clusterAroundXCells(
                   usageRate=nextCell_cKD$clusterAmount,
                   clusterData=data_stages[[stageIdx]][
                     data_stages[[stageIdx]]$cellType %in% unique(nextCell_cKD$clusterCell),],
                   cellVarData=as.numeric(nextCell_Vars),
                   stageName=stage, reduceEdge=reduceEdge,
                   cells=nextCell_cKD$cell,
                   clusterCells=nextCell_cKD$clusterCell,
                   kappas=nextCell_cKD$kappa,
                   minPts=1)
      )

      # Record generation
      completeCells <- c(completeCells, nextCell_cKD$cell)
    }

  }

  ## Organize and return
  data_ret <- .convertList2Dataframe(data_stages, typeBind = 'row')[,c(6,1:3,5,4)]
  data_ret$Stage <- as.character(data_ret$Stage)
  data_ret
}



.clusterAroundXCells <- function(usageRate,clusterData, cellVarData,
                                 stageName, reduceEdge,
                                 cells, clusterCells, kappas,
                                 minPts=1){
  newData <- data.frame()
  # Go through each cell
  for(i in 1:length(cells)){
    clusterCellData <- clusterData[clusterData$cellType==clusterCells[i],]

    # Go through each cell and develop clusters
    for(j in 1:nrow(clusterCellData)){
      if(rbinom(1,1,usageRate)){
        currXY <- as.numeric(clusterCellData[j, c('x','y')])
      } else{
        currXY <- runif(2)
      }
      data_pts <- .placeClusteredPts(currXY=currXY,
                                     cell=cells[i],
                                     numPts=rpois(1,kappas[i]),
                                     varValue=cellVarData[i],
                                     reduceEdge=reduceEdge)
      if(!is.null(data_pts)){
        data_pts$Image <- clusterCellData[j,'Image']
        data_pts$Person <- clusterCellData[j,'Person']
        data_pts$Stage <- clusterCellData[j,'Stage']

        newData <- rbind(newData, data_pts)
      }
    }

    # Require minPts
    #     Note, this can change your distribution if not thought about!
    while(nrow(newData[newData$cellType==cells[i],])<minPts){
      # Select a point from previous iteration
      preItrPt <- sample(nrow(clusterCellData),1)
      # Fill with enough pts
      numPts <- rpois(1,kappas[i])
      if(numPts+nrow(newData[newData$cellType==cells[i],]<minPts))
        numPts <- minPts - nrow(newData[newData$cellType==cells[i],])


      if(rbinom(1,1,usageRate)){
        currXY <- as.numeric(clusterCellData[j, c('x','y')])
      } else{
        currXY <- runif(2)
      }
      data_pts <- .placeClusteredPts(currXY=currXY,
                                     cell=cells[i],
                                     numPts=rpois(1,kappas[i]),
                                     varValue=cellVarData[i],
                                     reduceEdge=reduceEdge)
      data_pts$Image <- clusterCellData[j,'Image']
      data_pts$Person <- clusterCellData[j,'Person']
      data_pts$Stage <- clusterCellData[j,'Stage']

      newData <- rbind(newData, data_pts)
    }

  }

  newData
}

#############################
dat <- simulatePP1(cellVarData=
             data.frame('stage'=c(0,1),
                        'A'=c(0,0),
                        'B'=c(1/1250,1/1250)),
           cellKappaData=data.frame(
             'cell'=c('A','B','A','B'),
             'stage'=c(0,0,1,1),
             'clusterCell'=c(NA,NA,NA,'A'),
             'clusterAmount'=c(NA,NA,NA,0.5),
             'kappa'=c(7,15,7,7)),
           peoplePerStage=20,
           imagesPerPerson=5,
           reduceEdge=0.025,
           silent=F)

plotPP(dat[dat$Image==1,c('x','y','cellType')],xlim = c(0,1),ylim = c(0,1))
plotPP(dat[dat$Image==25,c('x','y','cellType')],xlim = c(0,1),ylim = c(0,1))
plotPP(dat[dat$Image==50,c('x','y','cellType')],xlim = c(0,1),ylim = c(0,1))
plotPP(dat[dat$Image==100,c('x','y','cellType')],xlim = c(0,1),ylim = c(0,1))

plotPP(dat[dat$Image==101,c('x','y','cellType')],xlim = c(0,1),ylim = c(0,1))
plotPP(dat[dat$Image==125,c('x','y','cellType')],xlim = c(0,1),ylim = c(0,1))
plotPP(dat[dat$Image==150,c('x','y','cellType')],xlim = c(0,1),ylim = c(0,1))
plotPP(dat[dat$Image==200,c('x','y','cellType')],xlim = c(0,1),ylim = c(0,1))

pcaData <- getPCAData(dat,repeatedUniqueId='Image',
                      xRange = c(0,1),  yRange = c(0,1), silent=F)
rfcv <- computeRandomForest_CVPC(data=pcaData,repeatedId='Image',
                                 cellData=dat)
