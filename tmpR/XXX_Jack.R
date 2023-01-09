myRF_CV_multi_singleimage <- function(data,
                                      K=nrow(data),
                                      plotVarImp=T,
                                      metaNames=NULL,
                                      MarkNames=NULL,
                                      cellData=NULL,
                                      dropCols=c('Person'),
                                      FakePairs=100, FakeMetas=100,
                                      generalFakeCell=T,
                                      ...){
  ## This generates fake data and fits and performs CV on a random forest

  ## Internal Function
  getChunks <- function(x,chunksN) split(x, sample(cut(x, chunksN, labels = FALSE)))

  ## Code

  ## Warnings
  if(!is.null(metaNames) && FakeMetas<1){
    warning('Warning: If PCA>1, Metas are not comparable to others')
  }

  ## Generate Fakes
  if(FakePairs>0){
    # Ensure requirements met
    if(sum(grepl("(Fake)*[_](Fake*)", colnames(data)))>0)
      stop('Error: Colname with Fake*_Fake* must be renamed')
    if(sum(c('Person','Stage') %in% colnames(cellData))!=2)
      stop('Error: Need Person and Stage in PcA data to match to fake')

    ## Get Number of principle components (PCs)
    nPCs <- max(suppressWarnings(as.numeric(
      sub(".*_PC", "", colnames(data)) )), na.rm = T)
    ## Determine KFunctions of interest
    tmp <- colnames(data)[!(colnames(data)%in%c(metaNames,dropCols,'Stage'))]
    cellsSetup <- str_extract_all(tmp, "[^_]*_[^_]*")
    KFunctions <- unique(unlist(cellsSetup)[-seq(2,2*length(cellsSetup),2)])

    if(length(KFunctions)<1){
      FakePairs <- 0
    }else{

      ## This makes one fake and compares all K-functions to it
      #     Generally effective, but can consider specific if
      #       cell-count is an issue (particular for K and crossK)
      data_unique <- unique(cellData[,c('Person','Stage')])
      fake_data <- NULL
      # Get the avg number of cells per type
      # kap_val <- round(mean(colSums(cellData[,MarkNames]==TRUE))/nrow(data_unique))

      # JH: Questions:
      # 1. does the kap_val we use need to be related to the average number of cells per image in the data?
      # Can we just set it to something fixed?
      # 2. Along these lines, is it ok to just use a 1x1 window?

      # For simplicity I have changed kap_val to a fixed value.

      kap_val <- 100

      # Generate Cells

      # I have changed the below just to use rpoispp to simplify things.
      # Maybe change kap_val to lambda to avoid confusion?

      # Also, for simplicity to avoid confusing myself I have just deleted the part that
      # deals with generalFakeCell=FALSE

      for(i in 1:nrow(data_unique)){
        print(c(i,data_unique[i,'Person']))
        for(p in 1:FakePairs){
          # Generate FakeCellL
          fake_cellsL.pp = rpoispp(lambda = kap_val,win=owin(c(0,1),c(0,1)))
          fake_cellsL = data.frame('x' = fake_cellsL.pp$x,
                                   'y' = fake_cellsL.pp$y,
                                   'cellType' = paste0('FakeCell',p,'L'))
          fake_cellsL$Person <- data_unique[i,'Person']
          fake_cellsL$Stage <- data_unique[i,'Stage']


          # Generate FakeCellR
          fake_cellsR.pp = rpoispp(lambda = kap_val,win=owin(c(0,1),c(0,1)))
          fake_cellsR = data.frame('x' = fake_cellsR.pp$x,
                                   'y' = fake_cellsR.pp$y,
                                   'cellType' = paste0('FakeCell',p,'R'))
          fake_cellsR$Person <- data_unique[i,'Person']
          fake_cellsR$Stage <- data_unique[i,'Stage']

          fake_data <- rbind(fake_data, fake_cellsL,fake_cellsR)
        }
      }

      # Setup for PCA
      Cells_df <- data.frame('c1'=paste0('FakeCell',1:FakePairs,'L'),
                             'c2'=paste0('FakeCell',1:FakePairs,'R'))


      fake_data[,'cellType'] = as.factor(fake_data[,'cellType'])

      # We use "_cross" here as the fake cell data consists of two cell types, rather than marked cells:

      fakeCellPCA <- getPCAData_cross(data=fake_data,
                                      outcome = "Stage",
                                      window = as.owin(c(0,1,0,1)), # defines a window in spatstat.
                                      rmax = NULL,
                                      nPCs = nPCs,
                                      Cells_df = Cells_df,
                                      edgeCorrection="isotropic",
                                      nbasis=21)

      data <- merge(data, fakeCellPCA,by = c('Person','Stage'))
      data <- data[order(as.numeric(substr(data$Person,2,100))),]
    }

  }

  if(FakeMetas>0 && is.null(metaNames)){
    # No point of FakeMeta without Meta
    FakeMetas <- 0
  }

  if(FakeMetas>0){
    # Ensure requirements met
    if(sum(grepl("FakeMeta", colnames(data)))>0)
      stop('Error: Colname with FakeMeta must be renamed')

    for(i in 1:length(metaNames)){
      for(m in 1:FakeMetas){
        data[[paste0('FakeMeta',m,metaNames[i])]] <-
          sample(data[[metaNames[i]]],replace = T)
      }

      ## Add FakeMeta to metaNames
      metaNames <-
        unique(c(metaNames, paste0('FakeMeta',1:FakeMetas,metaNames[i])))
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
  groups <- getChunks(1:n, K)

  ## Run CV
  for(i in 1:K){
    cat(paste0('Trial ',i,'/',K,'\n'))
    RF <- myRandomForest_multi(data[-groups[[i]],], varImpPlot = F,
                               metaNames=metaNames)

    avgGini[,i] <- RF[[2]]$avgGini
    avgVI[,i] <- RF[[2]]$avgVI
    oobAcc[i] <- sum(data[groups[[i]],'Stage']==
                       predict.myRandomForest(model = RF[[1]],
                                              data_pred = data[groups[[i]],],
                                              type = 'pred',
                                              data = data)) /
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
                                                 min(1,specify_decimal(oobAccData$u95,2)),')'
                                 ))
    ## Added to immediately display the desired plot
    avgVIGrob <- grid.arrange(avgVIPlot,
                              layout_matrix = rbind(c(1,1)),
                              bottom = paste0('Variable Importance - OOB-Acc ',
                                              min(1,max(0,specify_decimal(oobAccData$avgOOB, 2))),' (',
                                              max(0,specify_decimal(oobAccData$l95,2)),'-',
                                              min(1,specify_decimal(oobAccData$u95,2)),')'
                              ))

    return(list('Gini'=giniData, 'VI'=viData,
                'Accuracy'=oobAccData,
                'plotGrob'=avgVIGrob,
                'plotGrobOld'=varImportance,
                'plots'=list(avgGiniPlot,avgVIPlot)))
  }

  list('Gini'=giniData, 'VI'=viData, 'Accuracy'=oobAccData)
}
