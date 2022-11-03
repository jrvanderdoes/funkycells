#' Fit a Random Forest model with PC data (Using CV for Improvements)
#' 
#' The function fits a random forest model to the data along with using cross-
#'     validation to quantify variable importance.
#'
#' @param data Data.frame of outcome and predictors (PCs and meta-variables). 
#'     Note, currently Unit or repeated measures should not be included. 
#'     Generally use the results from getPCAData, potentially with meta-
#'     variables attached.
#' @param K 
#' @param plotVarImp 
#' @param metaNames 
#' @param cellData 
#' @param maxR 
#' @param precR 
#' @param dropCols 
#' @param FakePairs 
#' @param FakeMetas 
#' @param generalFakeCell 
#' @param reduceEdge 
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
computeRandomForest_CVPC <- function(data, K=nrow(data), plotVarImp=T,
                    outcome=colnames(data)[1],
                    unit=colnames(data)[2],
                    metaNames=NULL, 
                    cellData=NULL,
                    maxR=0.5, precR=1/100,
                    FakePairs=100, FakeMetas=100,
                    generalFakeCell=T,reduceEdge=0.05,
                    ...){
  ## Warnings
  if(!is.null(metaNames) && FakeMetas<1){
    warning('Warning: If PCA>1, Metas are not comparable to (PC) K functions')
  }
  
  ## Generate Fakes
  if(FakePairs>0) {
    ## Determine KFunctions of interest
    tmp <- colnames(data)[!(colnames(data)%in%c(outcome, unit, metaNames))]
    cellsSetup <- str_extract_all(tmp, "[^_]*_[^_]*")
    KFunctionsNames <- 
            unique(unlist(cellsSetup)[-seq(2,2*length(cellsSetup),2)])
    
    if(length(KFunctionsNames)<1){
      FakePairs <- 0
    }else{
      data <- .generateFakeKFunctions()
    }
  }
  
  if(FakeMetas>0){
    if(is.null(metaNames)){
      # No point of FakeMeta without Meta
      FakeMetas <- 0
    }else{
      data <- .generateFakeMetas()
      
      ## Add FakeMeta to metaNames
      metaNames <-  unique(c(metaNames, 
            paste0('FakeMeta',
                   rep(1:FakeMetas, each = length(metaNames)), 
                   metaNames)))
    }
  }
  
  
  if(!is.null(dropCols))
    data <- data[,!(colnames(data) %in% dropCols)]
  
  # Get columns to build
  vars <- colnames(data)
  for(i in 2:length(vars)){
    vars[i] <- substr(vars[i],1,
                      tail(unlist(gregexpr('_', vars[i])), n=1)-1)
  }
  ## Only keep data I'll use (this can throw off predictions)
  data <- data[vars != '' | colnames(data) %in% metaNames]
  
  # Drop Stage, simplify vars so any of the PC count to one variable
  vars <- unique(vars[-1])
  ## Drop blanks
  vars <- vars[vars!=""]
  
  ## Add in meta
  vars <- c(vars, metaNames)
  
  ## Prep Data
  avgGini <- data.frame(matrix(ncol=K, nrow=length(vars)))
  avgVI <- avgGini
  oobAcc <- rep(NA,K)
  
  # Setup groups for K-fold CV
  n <- nrow(data)
  groups <- .getChunks(1:n, K)
  
  
  ## Run CV
  for(i in 1:K){
    cat(paste0('Trial ',i,'/',K,'\n'))
    RF <- .computeRandomForest_PC(data=data[-groups[[i]],], 
                                  outcome=outcome,
                                  varImpPlot = F, 
                                  metaNames=metaNames, ...) 
    avgGini[,i] <- RF[[2]]$avgGini
    avgVI[,i] <- RF[[2]]$avgVI
    oobAcc[i] <- sum(data[groups[[i]],'Stage']==
                       .predict.RandomForest_PC(model = RF[[1]], 
                                data_pred = data[groups[[i]],], 
                                type = 'pred', data = data)) / 
                  nrow(data[groups[[i]],])
  }
  
  # Organize results
  gMeans <- rowMeans(avgGini)
  gSD <- apply(avgGini,MARGIN=1,FUN=function(x){sd(x)})
  giniData <- data.frame('var'=vars,
                         'avgGini'=gMeans,
                         'sdGini'=gSD,
                         'l95'= apply(avgGini, 1, quantile, probs = c(0.025)), #gMeans - 1.96 * gSD/sqrt(nrow(data)),
                         'u95'= apply(avgGini, 1, quantile, probs = c(0.975)),
                         'justUp95'= apply(avgGini, 1, quantile, probs = c(0.95)))
  viMeans <- rowMeans(avgVI)
  viSD <- apply(avgVI,MARGIN=1,FUN=function(x){sd(x)})
  viData <- data.frame('var'=vars,
                       'avgVI'=viMeans,
                       'sdVI'=viSD,
                       'l95'= apply(avgVI, 1, quantile, probs = c(0.025)),
                       'u95'= apply(avgVI, 1, quantile, probs = c(0.975)),
                       'justUp95'= apply(avgVI, 1, quantile, probs = c(0.95)))
  oobAccMeans <- mean(oobAcc)
  oobAccSD <- sd(oobAcc)
  ## TODO: Explain why these are different from others
  oobAccData <- data.frame('avgOOB'=oobAccMeans,
                           'sdOOB'=oobAccSD,
                           'l95q'= quantile(oobAcc, probs = c(0.025))[[1]],
                           'l95'= oobAccMeans - 1.96 * oobAccSD/sqrt(nrow(data)),
                           'u95q'= quantile(oobAcc, probs = c(0.975))[[1]],
                           'u95'= oobAccMeans + 1.96 * oobAccSD/sqrt(nrow(data)))
  
  # Plot Variable importance with confidence bands
  if(plotVarImp){
    # Setup
    # Check if there is a fake cell (If so, find biggest, mark and drop others)
    
    giniCutoff_df <- data.frame()
    viCutoff_df <- giniCutoff_df
    
    if(FakePairs>0){
      giniCutoff_df <- rbind(giniCutoff_df, 
                             data.frame('var'=KFunctions,'mean95th'=NA))
      viCutoff_df <- rbind(viCutoff_df, 
                           data.frame('var'=KFunctions,'mean95th'=NA))
      if(generalFakeCell){
        ## Gini
        FakeCellsLoc <- grep("FakeCell[1-9]*L_FakeCell[1-9]*R", giniData$var)
        giniCutoff_df[giniCutoff_df$var%in%KFunctions,'mean95th'] <- 
          quantile(giniData[FakeCellsLoc,'avgGini'],0.95)[[1]]
        giniData <- giniData[-FakeCellsLoc,]
        
        ## VI
        FakeCellsLoc <- grep("FakeCell[1-9]*L_FakeCell[1-9]*R", viData$var)
        viCutoff_df[viCutoff_df$var%in%KFunctions,'mean95th'] <- 
          quantile(viData[FakeCellsLoc,'avgVI'],0.95)[[1]]
        viData <- viData[-FakeCellsLoc,]
        
      }else{
        
        for(k in 1:length(KFunctions)){
          pat <- paste0('FakeCell',cellL[k],'[0-9]*L_FakeCell',cellR[k],'[0-9]*R')
          ## Gini
          FakeCellsLoc <- grep(pat, giniData$var)
          giniCutoff_df[giniCutoff_df$var==KFunctions[k],'mean95th'] <- 
            quantile(giniData[FakeCellsLoc,'avgGini'],0.95)[[1]]
          giniData <- giniData[-FakeCellsLoc,]
          
          ## VI
          FakeCellsLoc <- grep(pat, viData$var)
          viCutoff_df[viCutoff_df$var==KFunctions[k],'mean95th'] <- 
            quantile(viData[FakeCellsLoc,'avgVI'],0.95)[[1]]
          viData <- viData[-FakeCellsLoc,]
        }
      }
      
    } 
    
    if(FakeMetas>0){
      giniCutoff_df <- rbind(giniCutoff_df, 
                             data.frame('var'=metaNames[-grep('Fake',metaNames)],
                                        'mean95th'=NA))
      viCutoff_df <- rbind(viCutoff_df, 
                           data.frame('var'=metaNames[-grep('Fake',metaNames)],
                                      'mean95th'=NA))
      
      for(i in 1:length(metaNames[-grep('Fake',metaNames)])){
        FakeMetaLoc <- grep(paste0("FakeMeta.*",metaNames[i]), giniData$var)
        giniCutoff_df[giniCutoff_df$var==metaNames[i],'mean95th'] <- 
          quantile(giniData[FakeMetaLoc,'avgGini'],0.95)[[1]]
        giniData <- giniData[-FakeMetaLoc,]
        
        FakeMetaLoc <- grep(paste0("FakeMeta.*",metaNames[i]), viData$var)
        viCutoff_df[viCutoff_df$var==metaNames[i],'mean95th'] <- 
          quantile(viData[FakeMetaLoc,'avgVI'],0.95)[[1]]
        viData <- viData[-FakeMetaLoc,]
      }
    }
    
    if((FakePairs+FakeMetas)>0){
      ## Gini
      # Get Cutoff to be biggest number
      giniCellMean <- max(giniCutoff_df$mean95th)
      # Standardize all results
      giniData <- merge(giniData,giniCutoff_df)
      maxInd <- which(giniData$mean95th==giniCellMean)
      giniData[-maxInd,c(2,4,5,6)] <-
        giniData[-maxInd,c(2,4,5,6)] +
        (giniCellMean-giniData[-maxInd,'mean95th'])
      giniData <- giniData[,-7]
      
      ## VI
      # Get Cutoff to be biggest number
      viCellMean <- max(viCutoff_df$mean95th)
      # Standardize all results
      viData <- merge(viData,viCutoff_df)
      maxInd <- which(viData$mean95th==viCellMean)
      viData[-maxInd,c(2,4,5,6)] <-
        viData[-maxInd,c(2,4,5,6)] +
        (giniCellMean-viData[-maxInd,'mean95th'])
      viData <- viData[,-7]
    }
    
    ## Make plots
    maxValGini <- max(giniData$avgGini, giniCellMean)
    avgGiniPlot <- ggplot(data=giniData,
                          mapping=aes(x=reorder(var,avgGini), 
                                      y=ifelse(avgGini/maxValGini>1,1,
                                               ifelse(avgGini/maxValGini<0,0,
                                                      avgGini/maxValGini)))) + 
      geom_point(color='black') + 
      geom_errorbar(aes(ymin = ifelse(l95/maxValGini<0,0,l95/maxValGini),
                        ymax = ifelse(u95/maxValGini>1,1,u95/maxValGini)),
                    color='black', width=0.2) +
      geom_hline(aes(yintercept=giniCellMean/maxValGini), color='red', linetype='dashed') +
      coord_flip() +
      xlab(NULL) +
      ylim(c(0,1)) +
      ylab('Gini') +
      theme_bw()
    
    maxValFakeVi <- max(viData$avgVI, viCellMean)
    avgVIPlot <- ggplot(data=viData,
                        mapping=aes(x=reorder(var,avgVI), 
                                    y=ifelse(avgVI/maxValFakeVi>1,1,
                                             ifelse(avgVI/maxValFakeVi<0,0,
                                                    avgVI/maxValFakeVi)))) + 
      geom_point(color='black') + 
      geom_errorbar(aes(ymin = ifelse(l95/maxValFakeVi<0,0,l95/maxValFakeVi),
                        ymax = ifelse(u95/maxValFakeVi>1,1,u95/maxValFakeVi)),
                    color='black', width=0.2) +
      geom_hline(aes(yintercept=viCellMean/maxValFakeVi), color='red', linetype='dashed') +
      coord_flip() +
      xlab(NULL) +
      ylim(c(0,1)) +
      ylab('Variable Importance') +
      theme_bw()
    
    specify_decimal <- function(x, k) trimws(format(round(x, k), nsmall=k))
    
    varImportance <- arrangeGrob(avgGiniPlot,avgVIPlot,
                                 layout_matrix = rbind(c(1,2)),
                                 bottom = paste0('Variable Importance - OOB-Acc ',
                                                 min(1,max(0,specify_decimal(oobAccData$avgOOB, 2))),' (',
                                                 max(0,specify_decimal(oobAccData$l95,2)),'-',
                                                 min(1,specify_decimal(oobAccData$u95,2)),') : q(',
                                                 max(0,specify_decimal(oobAccData$l95q,2)),'-',
                                                 min(1,specify_decimal(oobAccData$u95q,2)),')'
                                 ))
    return(list('Gini'=giniData, 'VI'=viData,
                'Accuracy'=oobAccData, 
                'plotGrob'=varImportance,
                'plots'=list(avgGiniPlot,avgVIPlot)))
  }
  
  list('Gini'=giniData, 'VI'=viData, 'Accuracy'=oobAccData)
}


#' Title
#'
#' @return
#' @export
#'
#' @examples
.generateFakeKFunctions <- function(data, cellData,
                                    maxR=NULL, precR=NULL){
  # Ensure requirements met
  if(sum(grepl("(Fake)*[_](Fake*)", colnames(data)))>0)
    stop('Error: Colname with Fake*_Fake* must be renamed')
  if(is.null(maxR) || is.null(precR))
    stop('Error: Need maxR and precR to add fake cells')
  if(sum(c('Person','Stage') %in% colnames(cellData))!=2)
    stop('Error: Need Person and Stage in PcA data to match to fake')
  
  ## Get Number of principle components (PCs)
  nPCs <- max(suppressWarnings(as.numeric(
                  sub(".*_PC", "", colnames(data)) )), na.rm = T)

  if(generalFakeCell){
    ## This makes one fake and compares all K-functions to it
    #     Generally effective, but can consider specific if
    #       cell-count is an issue (particular for K and crossK)
    data_unique <- unique(cellData[,c('Image','Person','Stage')])
    fake_data <- NULL
    # Get the avg number of cells per type
    kap_val <- round(mean(table(cellData[,'cellType']))/nrow(data_unique))
    
    # Generate Cells
    for(i in 1:nrow(data_unique)){
      for(p in 1:FakePairs){
        # Generate FakeCellL
        fake_cellsL <- generateRandomPP(xRange = c(0,1), yRange = c(0,1),
                                        kappa = kap_val,
                                        percentAwayFromEdge = reduceEdge,
                                        cellType = paste0('FakeCell',p,'L'))
        fake_cellsL$Image <- data_unique[i,'Image']
        fake_cellsL$Person <- data_unique[i,'Person']
        fake_cellsL$Stage <- data_unique[i,'Stage']
        
        
        # Generate FakeCellR
        fake_cellsR <- generateRandomPP(xRange = c(0,1), yRange = c(0,1),
                                        kappa = kap_val,
                                        percentAwayFromEdge = reduceEdge,
                                        cellType = paste0('FakeCell',p,'R'))
        fake_cellsR$Image <- data_unique[i,'Image']
        fake_cellsR$Person <- data_unique[i,'Person']
        fake_cellsR$Stage <- data_unique[i,'Stage']
        
        fake_data <- rbind(fake_data, fake_cellsL,fake_cellsR)
      }
    }
    
    # Setup for PCA
    cell_df <- data.frame('c1'=paste0('FakeCell',1:FakePairs,'L'),
                          'c2'=paste0('FakeCell',1:FakePairs,'R'))
    
  }else{
    cellL <- rep(NA,length(KFunctionsNames))
    cellR <- rep(NA,length(KFunctionsNames))
    data_unique <- unique(cellData[,c('Image','Person','Stage')])
    fake_data <- NULL
    
    # Go over each unique K-Function 
    for(k in 1:length(KFunctionsNames)){
      tmp <- unlist(strsplit(KFunctionsNames[k],'_'))
      cellL[k] <- tmp[1]
      cellR[k] <- tmp[2]
      
      tab_tmp <- as.data.frame(table(cellData[,'cellType']))
      tab_tmp$Freq <- tab_tmp$Freq/nrow(data_unique)
      kapL <- tab_tmp[tab_tmp$Var1==cellL[k],'Freq']
      kapR <- tab_tmp[tab_tmp$Var1==cellR[k],'Freq']
      
      # Go through each image-person pair
      for(i in 1:nrow(data_unique)){
        # Get the number of cells for each type
        # tab_tmp <- as.data.frame(table(
        #   cellData[cellData$Image==data_unique$Image[i] &
        #              cellData$Person==data_unique$Person[i],'cellType']))
        # kapL <- tab_tmp[tab_tmp$Var1==cellL[k],'Freq']
        # kapR <- tab_tmp[tab_tmp$Var1==cellR[k],'Freq']
        
        for(p in 1:FakePairs){
          # Generate FakeCellL
          fake_cellsL <- generateRandomPP(xRange = c(0,1), yRange = c(0,1), 
                                          kappa = kapL, 
                                          percentAwayFromEdge = reduceEdge,
                                          cellType = paste0('FakeCell',cellL[k],p,'L'))
          fake_cellsL$Image <- data_unique[i,'Image']
          fake_cellsL$Person <- data_unique[i,'Person']
          fake_cellsL$Stage <- data_unique[i,'Stage']
          
          # Generate FakeCellR
          fake_cellsR <- generateRandomPP(xRange = c(0,1), yRange = c(0,1), 
                                          kappa = kapR, 
                                          percentAwayFromEdge = reduceEdge,
                                          cellType = paste0('FakeCell',cellR[k],p,'R'))
          fake_cellsR$Image <- data_unique[i,'Image']
          fake_cellsR$Person <- data_unique[i,'Person']
          fake_cellsR$Stage <- data_unique[i,'Stage']
          
          fake_data <- rbind(fake_data, fake_cellsL,fake_cellsR)
        }
      }
    }
    
    # Setup for PCA
    cell_df <- data.frame()
    for(k in 1:length(KFunctionsNames)){
      cell_df <- rbind(cell_df,
                       data.frame('c1'=paste0('FakeCell',cellL[k],1:FakePairs,'L'),
                                  'c2'=paste0('FakeCell',cellR[k],1:FakePairs,'R')))
    }
  }
  
  fakeCellPCA <- getPCAData(data=fake_data,
                            nPCs = nPCs, maxR = maxR,precR = precR,
                            xRange = c(0,1),yRange = c(0,1),
                            cell_df = cell_df)
  
  data <- merge(data, fakeCellPCA,by = c('Person','Stage'))
  data <- data[order(as.numeric(substr(data$Person,2,100))),]
  
}


#' Generate Fake Metavariables
#' 
#' This (internal) function bootstraps metavariables to create FakeMetas number
#'     of fake meta-variables for each meta-variable.
#' 
#' @param data Data.frame with (minimum) of columns names in metaNames. These 
#'     can be of any data type.
#' @param metaNames Vector with the column names that correspond to 
#'     meta-variables.
#' @param FakeMetas Numeric indicating the number of fake meta-variables per
#'     meta-variable.
#' 
#' @return
#' @export Data.frame of data with fakeMetas appended.
#'
#' @examples
#' # See code for computeRandomForest_CVPC. This is not an outward function so 
#' #     won't be viewable.
.generateFakeMetas <- function(data,metaNames,FakeMetas){
  # Ensure requirements met
  if(sum(grepl("FakeMeta", colnames(data)))>0)
    stop('Error: Colname with FakeMeta must be renamed')
  
  # Bootstrap each meta-variable
  for(i in 1:length(metaNames)){
    for(m in 1:FakeMetas){
      data[[paste0('FakeMeta',m,metaNames[i])]] <- 
                  sample(data[[metaNames[i]]],replace = T)
    }
  }
  
  data
}

#' Create K folds 
#' 
#' This (internal) function creates chunksN of scrambled x.
#'
#' @param x Vector of data (or often indices)
#' @param chunksN Numeric for the number of chunks, or folds, desired
#'
#' @return a list with length chunksN, with each containing approximately the 
#'     same number of points
#' @export
#'
#' @examples
#' .getChunks(1:10,3)
.getChunks <- function(x,chunksN) {
  if(chunksN>length(x))
    warning(paste0('Warning: Only ',length(x),' chunks are used due to data'))
  if(chunksN<2)
    return(sample(x))
  split(x, sample(cut(x, chunksN, labels = FALSE)))
}
