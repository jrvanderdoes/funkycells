########
#
## Overview
#
########
# Overview
#   Run this entire file to get all the functions in the environment
#   Each function should have a comment near the top as to its purpose
#
# TODO
#   Stop standardizing images
#
# Old Notes
#   Best Info: https://bioconductor.org/packages/release/bioc/vignettes/spicyR/inst/doc/spicy.html
#   Note: Traditionally first iteration of Thomas process is dropped ('A')
#   See also: BiocManager::install("SPIAT")

########
#
## Libraries
#
########
library(spicyR)
library(fda)
library(rpart)
library(randomForest)
library(ggplot2)
library(gridExtra)
library(pROC)
library(RColorBrewer)
library(refund)
library(tidyverse)
library(spatstat)

########
#
## Functions
#
########
#####
# Get Data
#####
## Wrote File for
getDiabetesData <- function(){
  # This function get the diabetes data and organizes it from spicyR
  
  #####
  # Get Data
  data("diabetesData")
  
  tmp <- as.data.frame(cellSummary(diabetesData))
  tmp1 <- as.data.frame(imagePheno(diabetesData))
  
  #####
  # Organize Data
  data <- data.frame('Person'=rep(NA,length(tmp$imageID)),
                     'Image'=NA,
                     'x'=NA, 'y'=NA,
                     'cellType'=NA,
                     'Stage'=NA)
  
  data$Image <- tmp$imageID
  data$x <- tmp$x
  data$y <- tmp$y
  data$cellType <- as.character(tmp$cellType)
  
  for(i in 1:length(data$Person)){
    data$Person[i] <-
      tmp1[as.character(tmp1$imageID)==as.character(data$Image[i]), 'case']
    data$Stage[i] <- 
      as.character(tmp1[as.character(tmp1$imageID)==as.character(data$Image[i]), 'stage'])
  }
  
  data
}

## Wrote File for
getPCAData <- function(data, nPCs=1, maxR=NA, precR=1/10,
                       xRange = c(0,max(data$x)), 
                       yRange = c(0,max(data$y)),
                       people=unique(data$Person),
                       cellTypes=NULL, cell_df = NULL,
                       edgeCorrection=T, nbasis=21){
  ## This function get the principal components based on some point process data
  
  ## Functions
  getPCA <- function(KFunctions, cells,  nPCs, evalPts, nbasis=21){
    # Replicate is each column
    # Check if any are missing
    dropIdx <- which(colSums(is.na(KFunctions))!=0)
    
    if(length(dropIdx)>0){ 
      # Any Replicate doesn't have K-function
      data_tmp <- KFunctions[,-dropIdx]
      
      if(ncol(data_tmp)<=1){ 
        # No multiple people have K-function so no PCA
        K_pca_scores <- data.frame('res'=NA)
      } else{ 
        # Drop the missing person and replace with NA after
        K_func <- Data2fd(argvals = evalPts,
                          y=as.matrix(data_tmp),
                          basisobj = 
                            create.bspline.basis(rangeval = range(evalPts),
                                                 nbasis = nbasis))
        K_pca <- pca.fd(K_func, nharm = nPCs)
        K_pca_scores <- insertMissingRows(K_pca$scores, dropIdx, nPCs)
      }
    }else{
      ## K-Functions all present
      K_func <- Data2fd(argvals = evalPts,
                        y=as.matrix(KFunctions),
                        basisobj = 
                          create.bspline.basis(rangeval = range(evalPts),
                                               nbasis = nbasis))
      K_pca <- pca.fd(K_func, nharm = nPCs)
      K_pca_scores <- K_pca$scores
    }
    colnames(K_pca_scores) <- paste0(cells[1],'_',cells[2],'_PC',1:nPCs)
    
    K_pca_scores
  }
  
  ## Start Code
  
  # Define pcaData
  pcaData <- unique(data[,c('Person','Stage')])
  
  # Define facts for K-functions
  if(is.na(maxR))
    maxR <- sqrt(((xRange[2]-xRange[1])/2)^2 + ((yRange[2]-yRange[1])/2)^2)
  
  if(precR!='E'){
    rCheckVals <- seq(0, maxR, precR)
  } else{
    rCheckVals <- maxR
  }
  
  ## Compute PCA for each cell-cell K-function
  
  if(is.null(cell_df)){
    if(is.null(cellTypes))
      stop('Error: cell_df or cellTypes must be non-null')
    if(length(cellTypes)>10)
      stop("Error: cellTypes too long, expansion may take a long time")
    
    cell_df <- expand.grid('c1'=cellTypes,
                           'c2'=cellTypes)
  }
  
  pcaData_list <- apply(as.data.frame(cell_df),
                        MARGIN=1,
                        FUN=function(cells, nPCs, rCheckVals, people,
                                     edgeCorrection, nbasis, ...){
                          cat(cells,'\n')
                          
                          ## Compute K-Function for each person
                          evaled_fd_K <- getKFunction('Person',
                                                      xRange=xRange, yRange=yRange,
                                                      people=people, cells=cells,
                                                      data=data, rCheckVals=rCheckVals,
                                                      edgeCorrection=edgeCorrection,
                                                      combineImage=F)
                          
                          ## Get PCA Scores
                          K_pca_scores <- getPCA(
                            KFunctions=evaled_fd_K, 
                            cells=cells, nPCs=nPCs, evalPts=rCheckVals,
                            nbasis=nbasis)
                          
                          as.data.frame(K_pca_scores)
                        }, 
                        nPCs=nPCs, rCheckVals=rCheckVals, people=people,
                        data=data,
                        edgeCorrection=edgeCorrection,
                        nbasis=nbasis)
  
  cbind(pcaData, convertList2Dataframe(pcaData_list))
}

## Try to integrate and delete
getRFData <- function(pcaData, dropNACol=F){
  ## This function reorganizes pcaData
  
  data_rf <- pcaData
  for(i in 3:ncol(data_rf)){
    data_rf[,i]<-as.numeric(data_rf[,i])
  }
  if(dropNACol)
    data_rf[ , colSums(is.na(data_rf))==0]
  
  data_rf
}

## Test when needed
getLMData <- function(data, byVar = c('Person','Image'),
                      maxR=NA, precR=1/10,
                      xRange = c(0,max(data$x)), 
                      yRange = c(0,max(data$y)),
                      people=unique(data$Person),
                      cells=c('ductal','ductal'),
                      edgeCorrection=T,
                      nbasis=21){
  ## This function does work on linear model approach
  
  ## Functions
  getLMData_Images <- function(data, rCheckVals,
                               xRange, yRange,
                               people=unique(data$Person),
                               cells=c('ductal','ductal'),
                               edgeCorrection=F,nbasis=21){
    
    ## Get K-Functions for each image
    eval_K <- sapply(people, function(person, data, rCheckVals,
                                      xRange, yRange,
                                      cells, edgeCorrection,
                                      nbasis){
      imgK <- getKFunction(type='Image', 
                           images = unique(data[data$Person==person,'Image']),
                           data=data[data$Person==person,],
                           cell1=cells[1], cell2=cells[2],
                           xRange=xRange, yRange=yRange,
                           rCheckVals = rCheckVals,
                           edgeCorrection=edgeCorrection,
                           combineImage = T,
                           nbasis=nbasis)
      colnames(imgK) <- unique(data[data$Person==person,'Image'])
      #imgK <- data.frame(rep(NA,length(rCheckVals)))
      #imgK <- convertList2Dataframe(imgKList)
      #for(i in 1:length(imgKList)){
      #  imgK[,i] <- convertKMatrix2Points(K=imgKList[[i]], 
      #                                    evalPts = evalPts, maxR=maxR)
      #}
      
      #as.matrix(imgK)
      #convertKMatrix2Points(K=imgK, evalPts=rCheckVals, nbasis=21)
      imgK
    },
    data=data,rCheckVals=rCheckVals,
    xRange=xRange, yRange=yRange,
    cells=cells,
    edgeCorrection=edgeCorrection,
    nbasis=nbasis,
    simplify = F)
    
    # Included for real data
    names(eval_K) <- people
    
    # Clean-up data for use in pffr
    data_return <- data.frame(eval_K[[1]])
    colnames(data_return) <- colnames(eval_K[[1]])
    People <- rep(names(eval_K)[1],ncol(eval_K[[1]]))
    Stages <- rep(unique(data[data$Person==names(eval_K)[1],'Stage']),ncol(eval_K[[1]]))
    Images <- colnames(eval_K[[1]])
    for(i in 2:length(eval_K)){
      data_return <- cbind(data_return,eval_K[[i]])
      People <- c(People, rep(names(eval_K)[i],ncol(eval_K[[i]])))
      Stages <- c(Stages,
                  rep(unique(data[data$Person==names(eval_K)[i],'Stage']),ncol(eval_K[[i]])))
      Images <- c(Images, colnames(eval_K[[i]]))
    }
    
    list('Y'=t(data_return),
         'Person'=as.factor(People),
         'Stage'=as.factor(Stages),
         'Image'=as.factor(Images), 
         'EvalPts'=rCheckVals)
  }
  
  ## Code
  if(length(byVar) >1)
    stop('Error: Please Select Person or Image for byVar')
  
  # Define facts for K-functions
  #xLen <- (max(data$x)-min(data$x))
  #yLen <- (max(data$y)-min(data$y))
  #maxR <- ceiling(sqrt((xLen/2)^2 + (yLen/2)^2))
  #area <- xLen * yLen
  if(is.na(maxR))
    maxR <- sqrt(((xRange[2]-xRange[1])/2)^2 + ((yRange[2]-yRange[1])/2)^2)
  
  if(precR!='E'){
    rCheckVals <- seq(0, maxR, precR)
  } else{
    rCheckVals <- maxR
  }
  
  if(byVar == 'Person'){
    eval_K <- getKFunction(type='Person', people=people,
                           cells=cells, data=data, 
                           xRange = xRange, yRange=yRange,
                           rCheckVals = rCheckVals,
                           edgeCorrection=edgeCorrection,
                           nbasis=nbasis)
    # Included for real data
    colnames(eval_K) <- people
    
    LMData <- list('Y'=t(eval_K),
                   'Stage'=as.factor(unique(data[,c('Person','Stage')])$Stage),
                   'Person'=as.factor(unique(data[,c('Person','Stage')])$Person),
                   'EvalPts'=rCheckVals)
    
  }else if(byVar == 'Image'){
    LMData <- getLMData_Images(data=data, rCheckVals=rCheckVals, 
                               xRange=xRange, yRange=yRange,
                               people=people,
                               cells=cells,
                               edgeCorrection=edgeCorrection,
                               nbasis=nbasis)
  } else{
    stop('Error: Only Person and Image available for byVar')
  }
  
  LMData
}

#####
# Random PP
#####
## Wrote file for
generateRandomPP <- function(xRange = c(0,1), yRange = c(0,1),
                             kappa=25, percentAwayFromEdge=0.05,
                             requireOne=T, cellType='A'){
  # This generates a simple random point process
  
  area <- (xRange[2]-xRange[1]) * (yRange[2]-yRange[1])
  intensity <- kappa*area
  numPts <- rpois(1, intensity)
  
  xReduceEdge <- (xRange[2]-xRange[1])*percentAwayFromEdge
  yReduceEdge <- (yRange[2]-yRange[1])*percentAwayFromEdge
  
  pointPattern <- data.frame(
    'x'=runif(numPts, min=xRange[1]+xReduceEdge, max=xRange[2]-xReduceEdge),
    'y'=runif(numPts, min=yRange[1]+yReduceEdge, max=yRange[2]-yReduceEdge),
    'cellType'=cellType)
  
  if(nrow(unique(pointPattern)) != nrow(pointPattern)){
    warning('Points placed on top of each other, so dropped (not replaced)')
    pointPattern <- unique(pointPattern) 
  }
  
  if(requireOne && nrow(pointPattern)==0)
    pointPattern <- generateRandomPP(xRange=xRange, yRange=yRange,
                                     kappa=kappa,percentAwayFromEdge=percentAwayFromEdge,
                                     requireOne=requireOne)
  
  pointPattern
}

## See if can delete
generateClusterAroundPP <- function(data,
                                    ptDist='Poisson', ptPars=list(list('Lam'=5)),
                                    placeDist='Normal', placePars=list(list('Var'=1/10000)),
                                    minPt=1, maxTries=3, xRange=c(0,1), yRange=c(0,1)){
  ## This clusters cells (Based on Thomas process originally)
  
  ## Functions
  getPtCount <- function(dist, pars){
    numPts <- 0
    if(dist=='Poisson'){
      numPts <- rpois(1,pars$Lam)
    } else{
      stop('Only Poisson allowed for ptDist currently')
    }
    
    numPts
  }
  
  placePts <- function(currXY, data, colVal, numPts, dist, pars, xRange, yRange){
    if(numPts==0)
      return()
    
    data_tmp <- data.frame('x'=rep(NA,numPts), 'y'=NA, 'cellType'=NA)
    done <- FALSE
    
    while(!done){
      compPts <- 0
      
      if(dist=='Normal'){
        while(compPts < numPts){
          data_tmp[compPts+1,'x'] <- rnorm(1,mean=currXY[1],sd=sqrt(pars$Var))
          data_tmp[compPts+1,'y'] <- rnorm(1,mean=currXY[2],sd=sqrt(pars$Var))
          
          # Ensure its in boundaries
          if(data_tmp[compPts+1,'x'] >= xRange[1] &
             data_tmp[compPts+1,'x'] <= xRange[2] & 
             data_tmp[compPts+1,'y'] >= yRange[1] & 
             data_tmp[compPts+1,'y'] <= yRange[2]){
            compPts <- compPts + 1
          }
        }
      } else if(dist=='Donut'){
        while(compPts < numPts){
          radius <- sqrt(runif(1,pars$minR2,pars$maxR2))
          theta <- 2*pi*runif(1)
          
          data_tmp$x[compPts+1] <- radius * cos(theta) + currXY[1]
          data_tmp$y[compPts+1] <- radius * sin(theta) + currXY[2]
          
          # Ensure its in boundaries
          if(data_tmp[compPts+1,'x'] >= xRange[1] &
             data_tmp[compPts+1,'x'] <= xRange[2] & 
             data_tmp[compPts+1,'y'] >= yRange[1] & 
             data_tmp[compPts+1,'y'] <= yRange[2]){
            compPts <- compPts + 1
          }
        }
      } else{
        stop('Sorry, only Normal or Donut placeDist currently allowed')
      }
      
      anyOverlap <- (sum(apply(data_tmp[,c('x','y')], MARGIN=1, 
                               function(row, dat){ 
                                 nrow(dat[dat$x==row[1] & dat$y==row[2],])>0
                               }, dat = data[,c('x','y')])) > 0 )
      
      if(!anyOverlap){
        data_tmp['cellType'] <- colVal
        #data_tmp[setdiff(names(data), names(data_tmp))] <- 0
        #data <- rbind(data,data_tmp)
        
        done <- TRUE
      }
      
    }
    
    data_tmp
  }
  
  # Get Columns of points
  iterations <- length(ptDist)
  cols <- LETTERS[1:(1+iterations)] 
  
  for(i in 1:iterations){
    # Make new points storage
    # Get points from previous iteration
    data_tmp <- data[data$cellType==cols[i],] 
    
    for(j in 1:nrow(data_tmp)){
      numPts <- getPtCount(ptDist[i], ptPars[[i]])
      data_new <- placePts(currXY=as.numeric(data_tmp[j, c('x','y')]), 
                           data=data, colVal=cols[i+1], numPts=numPts, 
                           dist=placeDist[i], pars=placePars[[i]],
                           xRange=xRange, yRange=yRange)
      
      data <- rbind(data, data_new)
    }
    
    ## Verify at least minPt was placed
    #     Note, this can change your distribution if not thought about!
    ct <- 0
    while(nrow(data[data$cellType==cols[i+1],])<minPt){
      ct <- ct + 1
      # Select a point from previous iteration
      preItrPt <- sample(nrow(data_tmp),1)
      
      numPts <- getPtCount(ptDist[i], ptPars[[i]])
      if(ct >= maxTries) # If we try to many times, just complete
        numPts <- minPt
      
      data_new <- placePts(currXY=as.numeric(data_tmp[preItrPt, c('x','y')]),  
                           data=data, colVal=cols[i+1],
                           numPts=numPts, dist=placeDist[i], pars=placePars[[i]],
                           xRange=xRange, yRange=yRange)
      
      data <- rbind(data, data_new)
    }
  }
  
  data
}

## See if can delete
generateReplicatedRandomPP <- function(images=100, 
                                       xRange = c(0,1), yRange = c(0,1),
                                       kappa=25, percentAwayFromEdge = 0.05,
                                       ...){
  ## This just create replicated simulated data
  
  
  pp_full <- data.frame()
  
  for(i in 1:images){
    # Sometimes with low lam never get a point, which breaks clustering
    #   but low levels may be wanted to ensure viewable data. So now just re-run
    pp_start <- data.frame()
    while(nrow(pp_start)==0){
      pp_start <- generateRandomPP(xRange=xRange, yRange=yRange,
                                   kappa=kappa, percentAwayFromEdge=percentAwayFromEdge)
    }
    pp_image <- generateClusterAroundPP(data=pp_start, 
                                        xRange=xRange, yRange=yRange,
                                        ...)
    pp_image$Image <- i
    
    pp_full <- rbind(pp_full, pp_image)
  }
  
  pp_full
}

## See if can delete
simulate3Stages <- function(){
  ## This simulates a process with 3 stages
  
  stageVar <- c(1,1/100,1/10000)
  finalData <- NULL
  for(i in 1:length(stageVar)){
    # This group has the change Var
    pp1 <- generateReplicatedRandomPP(images = 100,
                                      kappa=10, percentAwayFromEdge = 0.05, 
                                      ptDist = c('Poisson','Poisson'), 
                                      ptPars = list(list('Lam'=3),list('Lam'=5)),
                                      placeDist=c('Normal','Normal'), 
                                      placePars = list(list('Var'=1/500),list('Var'=stageVar[i])))
    # This group will have no relation
    pp2 <- generateReplicatedRandomPP(images = 100,
                                      kappa=10, percentAwayFromEdge = 0.05, 
                                      ptDist = c('Poisson'), 
                                      ptPars = list(list('Lam'=3)),
                                      placeDist=c('Normal','Normal'), 
                                      placePars = list(list('Var'=1/500)))
    ## Combine them
    pp2$cellType <- ifelse(pp2$cellType=='B','D',pp2$cellType)
    ppComb <- rbind(pp1, pp2)
    ppComb <- ppComb[ppComb$cellType!='A',]
    ppComb$Stage <-  as.character(i)
    ppComb$Image <- ppComb$Image + (i-1)*100
    ppComb$Person <- NA
    for(j in (1:100 + (i-1)*100)){
      ppComb[ppComb$Image==j,'Person'] <- paste0('p',((j-1) %% 10) + (i-1)*10+1)
    }
    
    if(is.null(finalData)){
      finalData <- ppComb
    }else{
      finalData <- rbind(finalData, ppComb)
    }
  }
  
  finalData
}

## Wrote file for
plotPP <- function(data, colorGuide = NULL,ptSize=1,
                   xlim=c(min(data$x),max(data$x)),
                   ylim=c(min(data$y),max(data$y)),
                   dropAxes=F,
                   colors=NULL){
  ## This is used to plot point process data
  
  retPlot <- ggplot() +
    geom_point(mapping=aes(x=x,y=y, col=cellType), data=data,
               size=ptSize) +
    theme_bw() +
    xlim(xlim) + 
    ylim(ylim) +
    guides(color=colorGuide) 
  if(dropAxes)
    retPlot <- retPlot +
      theme(axis.title = element_blank(),
            axis.text = element_blank())
  if(!is.null(colors))
    retPlot <- retPlot +
      scale_color_manual(values = colors)
  
  retPlot
}

#####
# Example and Simulations (Unused)
#####
## See if can delete
compare_K_functions_examples <- function(example=c(1,2,3)){
  # This just plays around with some K examples
  
  if(length(example)!=1)
    stop('Error: Select one example')
  
  if(example==1){
    ## All Works
    kap <- 15 
    sig2 <- 100
    mu <- 10
    xRange <- c(0,1)
    yRange <- c(0,1)
    
  } else if(example==2){
    ## Shows My Improvement
    kap <- 25 
    sig2 <- 100
    mu <- 1#10
    xRange <- c(0,1)
    yRange <- c(0,1)
    
  } else if(example==3){
    ## Clustering
    kap <- 25 #25 
    sig2 <- 1/1000 #100
    mu <- 10
    xRange <- c(0,1)
    yRange <- c(0,1)
    
  }else{
    stop('Error: Select an example 1, 2, or 3')
  }
  
  ## Generate Data
  pp <- generateRandomPP(kappa=kap, xRange = xRange, yRange = yRange)
  pp2 <- generateClusterAroundPP(
    data=pp,
    ptDist='Poisson',
    ptPars=list(list('Lam'=mu)),
    placeDist='Normal',
    placePars=list(list('Var'=sig2)),
    minPt=1,
    maxTries=3,
    xRange = xRange, 
    yRange = yRange)
  
  plot_pp <- plotPP(pp2)
  
  ## Spatstat
  tmp<-as.ppp(pp2[pp2$cellType=='B',c('x','y')], square(1)) ## Update with win
  K2 <- thomas.estK(tmp)
  
  ## My Code
  K3 <- getKFunction(type='All',cell1='B',cell2='B',data = pp2[pp2$cellType=='B',], 
                     precR=1/10000, maxR = 0.25,
                     xRange = xRange, yRange = yRange,
                     edgeCorrection=T, nbasis=21)
  r1 <- seq(0, 0.25, 1/10000) ## update maxR with ratio
  
  ## Plot and Compare
  plot_K <- ggplot() +
    geom_line(aes(x=r, y=fit, col='Emp'), data=K2$fit) +
    geom_line(aes(x=r, y=pi*r^2 + (1-exp(-r^2/(4*sig2)))/kap, col='Theo'),
              linetype='dashed', data=K2$fit) + 
    geom_line(aes(x=r1, y=K3, col='Mine')) +
    theme_bw()
  ## Failure with low intensities
  
  ## Convert to L
  L2 <- Lest(tmp)
  L3 <- convertK2L(K3, r1)
  plot_L <- ggplot() +
    geom_line(aes(x=r, y=iso, col='Emp'),data=L2) +
    geom_abline(mapping=aes(col='45D',intercept = 0, slope = 1), linetype='dashed') +
    geom_line(aes(x=r1, y=L3, col='Mine')) +
    theme_bw()
  
  ## Convert K to L-line
  L4 <- convertK2L(K3, r1, remove.time=T)
  plot_Lline <- ggplot() +
    geom_line(aes(x=r, y=iso-r, col='Emp'),data=L2) +
    geom_hline(aes(yintercept=0, col='No effect'), linetype='dashed') +
    geom_line(aes(x=r1, y=L4, col='Mine')) +
    theme_bw()
  
  list(plot_pp, plot_K, plot_L, plot_Lline)
}

## See if can delete
noEffect_CI_sim <- function(){
  ## This plays around with CI
  
  numPeople <- 10
  imgVector <- c(1,5,25,125)
  results <- list()
  for(iter in 1:length(imgVector)){
    numImgPerson <- imgVector[iter]
    pp <- generateReplicatedRandomPP(images = numImgPerson*numPeople,
                                     kappa=10, percentAwayFromEdge = 0.05,
                                     ptDist = c('Poisson','Poisson'), 
                                     ptPars = list(list('Lam'=3),list('Lam'=5)),
                                     placeDist=c('Normal','Normal'), 
                                     placePars = list(list('Var'=1/500),list('Var'=1)))
    pp$Stage <- numImgPerson
    pp$Person <- NA
    for(i in 1:(numImgPerson*numPeople)){
      pp[pp$Image==i,'Person'] <- paste0('p',((i-1) %% 10) + 1)
    }
    
    data_lm <- getLMData(data = pp,byVar = 'Image',maxR = 0.2,precR=1/100,
                         xRange=c(0,1),yRange=c(0,1),
                         people=unique(pp$Person),cells = c('B','C'),
                         edgeCorrection = T, nbasis=21)
    data_lm_l <- list(Y=convertK2L(data_lm$Y,data_lm$EvalPts),
                      Person=data_lm$Person,
                      Stage=data_lm$Stage,
                      Image=data_lm$Image,
                      EvalPts=data_lm$EvalPts)
    model <- pffr(Y ~ -1 + Stage + c(s(Person,bs='re')) + c(s(Image,bs='re')),
                  yind=data_lm_l$EvalPts,
                  data=data_lm_l)
    plots <- plotLMCoefs(coef(model), reCt=2)
    
    results[[iter]] <- list('pp'=pp, 'data_lm'=data_lm, 'data_lm_l'=data_lm_l, 
                            'model'=model, 'plots'=plots, 
                            'numImgPerson'=numImgPerson)
  }
  
  results
}

## See if can delete
varImp_sim <- function(){
  ## This plays around with VarImp
  
  sim_vi <- simulate3Stages()
  sim_vi_pca <- getPCAData(sim_vi, nPCs=3, maxR=0.4, precR=1/100,
                           xRange = c(0,1), yRange = c(0,1),
                           cellTypes=c('B', 'C','D'),
                           #cell_df = data.frame('cell1'=c('B','B'),
                           #                      'cell2'=c('C','D'))
                           fakeCell = 'Fake',
                           edgeCorrection=T, nbasis=21)
  dsim_vi_rf <- getRFData(sim_vi_pca, dropNACol=F) 
  
  use_data <- cbind(dsim_vi_rf, data.frame('Person'=sim_vi_pca$Person))
  use_data$Age <- c(rep(20,5), rep(25,10), rep(30,10), rep(35,5))
  use_data$FakeMeta <- rbinom(30,1,0.5)
  
  cat("-----------\n")
  
  # B-C (and related), Age matters, Others not
  cvResult <- myRF_CV(data=use_data,fakeCell = 'Fake',
                      metaNames = c('Age','FakeMeta'), fakeMeta = 'FakeMeta')
  
  list(sim_vi, sim_vi_pca, cvResult)
}

#####
# Conversion Functions
#####
## See if can delete
convertK2L <- function(Kmatrix, evalPts, remove.time=F){
  ## This converts a K function into an L Function
  
  # Make Lmatrix (entries not negative)
  Lmatrix <- sqrt(pmax(Kmatrix,0)/pi)
  
  if(remove.time){
    Lmatrix <- mapply(function(LRow,time){
      LRow - time
    }, LRow=Lmatrix, time=evalPts)
  }
  
  Lmatrix
}

## Wrote function for
convertList2Dataframe <- function(data_list, na.omit=F){
  ## This converts a list into a DF
  
  if(dim(data_list[[1]])[[1]]>1){
    data_df <- data_list[[1]]
  }else{
    data_df <- data.frame('V1'=t(data_list[[1]]))
  }
  
  if(length(data_list)!=1){
    for(ii in 2:length(data_list)){
      if(dim(data_list[[ii]])[[1]]>1){
        data_df <- cbind(data_df,data_list[[ii]])
      }else{
        data_df <- cbind(data_df, t(data_list[[ii]]))
      }
    }
  }
  if(na.omit)
    data_df <- na.omit(data_df)
  
  data_df
}

## See if can delete
convertKMatrices2DF <- function(K, rValsKMatrix, firstColRVals = F){
  ## This converts the individual K matrices into a cohesive dataframe
  
  # data_ret will have K functions down cols, rownames will be evaled pts
  if(length(rValsKMatrix)==1){
    # Get all r vals
    rVals <- c()
    for(i in 1:length(K)){
      rVals <- c(rVals, K[[i]][,1])
    }
    rVals <- unique(0, rVals, rValsKMatrix)
    # Rename rValsKMatrix to use
    rValsKMatrix <- rVals[order(rVals)]
  }
  
  # Get K function at each r
  data_ret <- data.frame(matrix(NA, ncol=length(K), nrow=length(rValsKMatrix)))
  tryCatch({
    data_ret <- sapply(K, function(Kmat,rValsKMatrix){
      # Setup data for this iter
      data_ret <- rep(0, length(Kmat$r))
      
      # Take the max K for the r that is <= the given r
      for(j in 1:length(Kmat$r)){
        data_ret <- ifelse(Kmat$r[j] >= rValsKMatrix,
                           data_ret,
                           Kmat$K[j])
      }
      data_ret
    }, rValsKMatrix=rValsKMatrix)
  }, error=function(e){
    # View(K)
    # cat(length(K))
    # cat(rep(NA,length(K)))
    # View(as.data.frame(rValsKMatrix))
    # readline('p')
  })
  
  rownames(data_ret) <- rValsKMatrix
  
  if(firstColRVals)
    return(data.frame('r'=rValsKMatrix,data_ret))
  else
    return(data_ret)
}

## See if can delete
convertKMatrix2Points <- function(K, evalPts, nbasis=21){
  ## This converts the K matrix to specific evaluation points
  
  tmp = Data2fd(argvals = evalPts,#/max(evalPts),
                y=as.matrix(K),
                basisobj = 
                  create.bspline.basis(rangeval = range(evalPts),#c(0,1),
                                       nbasis = nbasis))
  
  eval.fd(evalPts,mean.fd(tmp))
}

#####
# LM Functions
#####

## See if can delete
plotLMCoefs <- function(coefs, reCt=0, intercept=F){
  ## Playing with LM model
  
  lenSM <- length(coefs$smterms)
  result <- list()
  
  se_re <- 0
  if(reCt>0){
    for(i in 0:(reCt-1)){
      se_re <- se_re + coefs$smterms[[lenSM-i]]$se[1]
    }
  }
  
  for(i in (1+intercept):(lenSM-reCt)){
    result[[i]] <- 
      ggplot(data=coefs$smterms[[i]]$coef) +
      geom_abline(intercept = coefs$smterms[[i]]$coef$value[1], 
                  slope = 1, col='red', linetype='dashed') +
      geom_line(aes(x=coefs$smterms[[i]]$coef[,1], y=value)) +
      geom_ribbon(aes(x=coefs$smterms[[i]]$coef[,1], 
                      ymin=value-1.96*(se + se_re), 
                      ymax=value+1.96*(se + se_re)),
                  alpha=0.5) +
      ylim(c(min(coefs$smterms[[i]]$coef$value - 
                   1.96*(coefs$smterms[[i]]$coef$se+se_re)) -
               min((coefs$smterms[[i]]$coef$se+se_re)),
             max(coefs$smterms[[i]]$coef$value + 
                   1.96*(coefs$smterms[[i]]$coef$se+se_re))+
               max((coefs$smterms[[i]]$coef$se+se_re)))) +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5),
            axis.text = element_text(size=14),
            axis.title = element_text(size=20)) +
      xlab("Standardized Distance") +
      ggtitle(coefs$smterms[[i]]$main[[1]])
  }
  
  result
}

#####
# Utility Functions
#####
## Wrote function for
insertMissingRows <- function(data_add, insertRows,npca=3){
  ## Not sure if this is used anymore. To fill missing.
  
  data_return <- matrix(ncol=npca, 
                        nrow=(length(data_add[,1])+length(insertRows)))
  currentAdd <- 1 # Current row I take data from
  
  for(i in 1:length(data_return[,1])){
    if(!(i %in% insertRows)){ # If this row shouldn't be NA
      data_return[i,] <- as.numeric(data_add[currentAdd,])
      currentAdd <- currentAdd + 1
    }
  }
  
  data_return
}

## Wrote function for
getKFunction <- function(type, data=data, reps=NULL, ...){
  ## Gets a K function based on some type (calls itself as needed)
  
  ## Functions
  getKFunction_Images <- function(images, data, 
                                  cell1, cell2,
                                  xRange, yRange, rCheckVals,
                                  edgeCorrection,
                                  combineImage = F,
                                  nbasis=21){
    ## This returns Kmatrix for an image
    tmp <- sapply(images, 
                  function(img, data, cell1, cell2,
                           xRange, yRange, rCheckVals,
                           edgeCorrection, combineImage,
                           nbasis){
                    imageKFunctions <- getKMatrices(
                      data = data[data$Image==img,], 
                      rCheckVals=rCheckVals, 
                      cell1 = cell1, 
                      cell2 = cell2,
                      xRange = xRange, 
                      yRange = yRange,
                      edgeCorrection = F)
                    
                    if(length(imageKFunctions)==1 && 
                       is.na(sum(imageKFunctions[[1]])))
                      return()
                    
                    tmp <- convertKMatrices2DF(imageKFunctions, rCheckVals, 
                                               firstColRVals = T)
                    
                    imageKFunctions_std <- tmp[,-1]
                    imageKFunctions_std <- imageKFunctions_std[colSums(!is.na(imageKFunctions_std))>0]
                    if(length(rCheckVals)==1)
                      rCheckVals <- tmp[,1]
                    
                    if(combineImage){
                      imageKFunctions_std <- convertKMatrix2Points(
                        K = imageKFunctions_std, 
                        evalPts = rCheckVals, 
                        nbasis=nbasis)
                    }
                    
                    imageKFunctions_std
                  }, 
                  data=data,
                  cell1=cell1, cell2=cell2,
                  xRange=xRange, yRange=yRange, 
                  rCheckVals=rCheckVals,
                  edgeCorrection=edgeCorrection,
                  combineImage=combineImage,
                  nbasis=nbasis,
                  simplify = F)
    
    convertList2Dataframe(tmp)
  }
  
  getKFunction_Person <- function(people, cells, data, 
                                  xRange, yRange,
                                  rCheckVals, edgeCorrection,
                                  combineImage = F, nbasis=21){
    
    ## Compute K-Function for each person
    evaled_fd_K <- sapply(people, function(person, cells,
                                           data,
                                           xRange, yRange,
                                           rCheckVals,
                                           edgeCorrection,
                                           combineImage,
                                           nbasis){
      ## CODE BEGIN
      ## Compute K-Functions for person (do NOT combine image K functions into 1)
      singlePersonKFunctions <- getKFunction('Image',
                                             images=unique(data[data$Person==person,"Image"]),
                                             data = data[data$Person==person, ], 
                                             cell1=cells[[1]], cell2=cells[[2]],
                                             xRange=xRange, yRange=yRange,
                                             rCheckVals=rCheckVals,
                                             edgeCorrection=edgeCorrection,
                                             combineImage=combineImage,
                                             nbasis=nbasis)
      
      if(is.list(singlePersonKFunctions) && !is.data.frame(singlePersonKFunctions)){
        if(length(singlePersonKFunctions)==0)
          evaled_fd_K <- rep(NA,length(rCheckVals))
        else
          K_mat_person <- convertList2Dataframe(singlePersonKFunctions, na.omit=T)
      } else{
        K_mat_person <- singlePersonKFunctions 
      }
      
      ## Get K Function (avg) if exists
      if(sum(!is.na(K_mat_person))>1){ # Not all NA
        ## Get K function for cell-cell on this person
        evaled_fd_K <- convertKMatrix2Points(
          K=K_mat_person, evalPts=rCheckVals, nbasis=nbasis)
      } else{
        evaled_fd_K <- rep(NA,length(rCheckVals))
      }
      
      evaled_fd_K
    },
    cells=cells, data=data, 
    xRange=xRange, yRange=yRange,
    rCheckVals=rCheckVals,
    edgeCorrection=edgeCorrection,
    combineImage=combineImage, nbasis=nbasis )
    
    evaled_fd_K
  }
  
  getKFunction_All <- function(data, cell1, cell2, 
                               maxR=NA, xRange=c(0,1), yRange=c(0,1),
                               precR = 1/100, edgeCorrection=F, nbasis=21){
    
    if(is.na(maxR))
      maxR <- sqrt(((xRange[2]-xRange[1])/2)^2 + ((yRange[2]-yRange[1])/2)^2)
    
    if(precR!='E'){
      rCheckVals <- seq(0, maxR, precR)
    } else{
      rCheckVals <- maxR
    }
    
    # Get list of K matrices
    repKFunctions <- getKMatrices(
      data = data,
      cell1 = cell1, 
      cell2 = cell2,
      xRange = xRange, 
      yRange = yRange,
      rCheckVals = rCheckVals,
      edgeCorrection=edgeCorrection)
    
    if(length(repKFunctions)==1 && is.na(repKFunctions))
      return()
    
    tmp <- convertKMatrices2DF(repKFunctions, rCheckVals, 
                               firstColRVals = T)
    
    repKFunctions_std <- tmp[,-1]
    repKFunctions_std <- repKFunctions_std[colSums(!is.na(repKFunctions_std))>0]
    
    repKFunctions_std <- convertKMatrix2Points(
      K = repKFunctions_std[colSums(!is.na(repKFunctions_std))>0], 
      evalPts = rCheckVals, 
      nbasis=nbasis)
    
    repKFunctions_std
  }
  
  ## Code
  if(type=='Person'){
    retVal <- getKFunction_Person(data=data, ...)
  }else if(type=='Image'){
    retVal <- getKFunction_Images(data=data, ...)
  }else if(type=='All'){
    retVal <- getKFunction_All(data=data, ...)
  }
  
  retVal
}

## Can delete
getKMatrices <- function(data, rCheckVals, 
                         cell1, cell2=NA, 
                         xRange = c(0,1), yRange = c(0,1),
                         edgeCorrection = F,
                         ...){
  ## Gets the individual K matrices
  
  ## Functions
  ## TODO:: Some roundoff possibility
  dist <- function(c2,c1x,c1y){
    tmp <- max(abs(min(c2[2]-c1y, c2[1]-c1x)), 1)
    abs(tmp) * sqrt( ((c2[2]-c1y)/tmp)^2 + ((c2[1]-c1x)/tmp)^2)
  }
  
  ## Data
  if(is.na(cell2))
    cell2 <- cell1
  
  cell1Data <- data[data$cellType==cell1,]
  cell2Data <- data[data$cellType==cell2,]
  
  if(nrow(cell1Data)==0 || nrow(cell2Data)==0){
    return(rep(NA, length(rCheckVals)))
  }
  
  # Vars
  maxR <- max(rCheckVals)
  numAgents <- nrow(cell2Data)#+nrow(cell1Data) # TODO:: Verify this choice
  lambdaR <- numAgents/( (xRange[2]-xRange[1]) * (yRange[2]-yRange[1]) )#/(pi*maxR^2)
  
  if(edgeCorrection){
    
    cell2Data_u <- cell2Data
    cell2Data_u$y <- cell2Data_u$y + yRange[2]
    cell2Data_d <- cell2Data
    cell2Data_d$y <- cell2Data_d$y - yRange[2]
    cell2Data_r <- cell2Data
    cell2Data_r$x <- cell2Data_r$x + xRange[2]
    cell2Data_l <- cell2Data
    cell2Data_l$x <- cell2Data_l$x - xRange[2]
    
    cell2Data_ur <- cell2Data
    cell2Data_ur$x <- cell2Data_ur$x + xRange[2]
    cell2Data_ur$y <- cell2Data_ur$y + yRange[2]
    cell2Data_ul <- cell2Data
    cell2Data_ul$x <- cell2Data_ul$x - xRange[2]
    cell2Data_ul$y <- cell2Data_ul$y + yRange[2]
    cell2Data_dr <- cell2Data
    cell2Data_dr$x <- cell2Data_dr$x + xRange[2]
    cell2Data_dr$y <- cell2Data_dr$y - yRange[2]
    cell2Data_dl <- cell2Data
    cell2Data_dl$x <- cell2Data_dl$x - xRange[2]
    cell2Data_dl$y <- cell2Data_dl$y - yRange[2]
    
    cell2Data <- rbind(cell2Data,
                       cell2Data_u, cell2Data_d,
                       cell2Data_r, cell2Data_l,
                       cell2Data_ur, cell2Data_ul,
                       cell2Data_dr, cell2Data_dl)
    
    ## Remove all those "too far out"
    #     That would be outside maxR distance from edges
    cell2Data <- cell2Data[cell2Data$x<= max(xRange)+maxR &
                             cell2Data$x>= min(xRange)-maxR,]
    cell2Data <- cell2Data[cell2Data$y<= max(yRange)+maxR &
                             cell2Data$y>= min(yRange)-maxR,]
  }
  
  # Prep vars
  cell1s <- nrow(cell1Data)
  cell2s <- nrow(cell2Data)
  
  if(cell2s<=20000){
    # This method is faster, but may break with large data
    
    # Get c1 (col) to c2 (row) distances 
    dists <- mapply(function(c1x,c1y,c2){ dist(c2,c1x=c1x,c1y=c1y)[[1]] }, 
                    c1x = cell1Data[,'x'],
                    c1y = cell1Data[,'y'], 
                    c2=list(cell2Data[,c('x','y')]))
    # Remove 0's (same cell), if any
    dists <- ifelse(dists==0,NA,dists)
    
    # Find K values
    K <- sapply(as.data.frame(dists),function(dat, maxR, lambdaR){
      dist_tmp <- na.omit(dat[dat<maxR])
      dist_tmp <- dist_tmp[order(dist_tmp)]
      if(length(dist_tmp)>0){
        retVal <- data.frame('r'=dist_tmp,
                             'K'=(1:length(dist_tmp))/lambdaR)
      } else{
        retVal <- data.frame('r'=NA,
                             'K'=NA)
      }
      retVal
    }, maxR=maxR, lambdaR=lambdaR, simplify = F)
    
  } else{
    K <- apply(cell1Data[(colnames(cell1Data) %in% c('x','y'))], 
               MARGIN=1, FUN=function(c1, c2s, maxR, lambdaR){
                 
                 dists <- apply(c2s,MARGIN=1, FUN=function(c2,c1x, c1y){
                   val <- dist(c2,c1x=c1x,c1y=c1y)[[1]]
                   ifelse(val<=10^-8,NA,val)
                 }, c1x=as.numeric(c1[1]),c1y=as.numeric(c1[2]))
                 
                 
                 if(length(dists)>0){
                   K_tmp <- sapply(as.data.frame(dists),function(dat, maxR, lambdaR){
                     dist_tmp <- na.omit(dat[dat<maxR])
                     dist_tmp <- dist_tmp[order(dist_tmp)]
                     if(length(dist_tmp)>0){
                       retVal <- data.frame('r'=dist_tmp,
                                            'K'=(1:length(dist_tmp))/lambdaR)
                     } else{
                       retVal <- data.frame('r'=NA,
                                            'K'=NA)
                     }
                   }, maxR=maxR, lambdaR=lambdaR, simplify = F)
                 } else{
                   K_tmp <- list(NA)
                 }
                 
                 K_tmp[[1]]
               }, 
               c2s=cell2Data[,c('x','y')], maxR=maxR, 
               lambdaR=lambdaR,simplify = F)
    
  }
  
  K
}

#####
# myRandomForest Code
#####
# Wrote function for
myRandomForest <- function(data, trees = 1000, varImpPlot=T, metaNames=NULL){
  ## My implementation of a random forest
  
  # Setup
  # columns: all columns in data
  # cols_preds: all predictors name
  columns <- colnames(data)
  for(i in 2:length(columns)){
    columns[i] <- substr(columns[i],1,
                         tail(unlist(gregexpr('_', columns[i])), n=1)-1)
  }
  columns <- columns[columns!=""]
  columns <- c(columns,metaNames)
  
  cols_preds <- unique(columns[-1])
  
  data_result <- data.frame('var'=cols_preds,
                            'splits'=0,
                            'giniDec'=0,
                            'varImp'=0,
                            'varImpCt'=0,
                            'vi'=0)
  RF <- list()
  if(length(unique(data$Stage))==1)
    stop('Error: Only 1 stage in data, cannot do RF')
  # To do CART
  for(i in 1:trees){
    # Subset data
    data_rf <- data[sample(1:nrow(data),nrow(data), replace = T),]
    # Stupid check. If only one-class, rebuild
    while(length(unique(data_rf$Stage))==1){
      data_rf <- data[sample(1:nrow(data),nrow(data), replace = T),]
    }
    
    # Subset predictors
    preds <- unique(sample(cols_preds, 0.8*length(cols_preds), 
                           replace = F))
    data_rf <- data_rf[columns %in% c('Stage', preds)]
    data_rf <- cbind(data_rf[1],data_rf[,sample(2:ncol(data_rf))])
    
    # Fit CART
    model <- rpart(Stage ~ . , data=data_rf, method="class",
                   control =rpart.control(minsplit =1,minbucket=1, cp=0))
    
    # Find Error at nodes (Gini is 1-p_1^2-p_2^2)
    frame <- model$frame
    frame[['gini']] = 1 - (frame[['dev']] / frame[['n']])^2 - 
      (1 - frame[['dev']] / frame[['n']])^2
    #frame[,c('var','n','dev','gini')]
    
    # Get Decrease for each split
    frame[['improve']] = NA
    for (j in 1:nrow(frame)) {
      if (frame[j,'var'] == '<leaf>') next
      
      ind = which(rownames(frame) %in% (as.numeric(rownames(frame)[j])*2+c(0,1)))
      frame[j,'improve'] = frame[j,'n']*frame[j,'gini'] - frame[ind[1],'n']*frame[ind[1],'gini'] - frame[ind[2],'n']*frame[ind[2],'gini']
      
      # Get Gini improvement
      name <- frame[j,'var']
      
      # Get name of bigger group
      if(sum(gregexpr('_', name)[[1]]!=-1)>0){
        # If cell interaction
        name_sub <- substring(name,1,tail(unlist(gregexpr('_', name)), n=1)-1)
      }else{
        # If meta data
        name_sub <- name
      }
      
      data_result[data_result$var==name_sub,'giniDec'] <- 
        data_result[data_result$var==name_sub,'giniDec'] +
        frame[j,'improve']
      # Count split by variable
      data_result[data_result$var==name_sub,'splits'] <- 
        data_result[data_result$var==name_sub,'splits'] + 1
    }
    
    # Get and count Var Import
    for(j in 1:length(model$variable.importance)){
      name <- names(model$variable.importance[j])
      
      # Get name of bigger group
      if(sum(gregexpr('_', name)[[1]]!=-1)>0){
        # If cell interaction
        name_sub <- substring(name,1,tail(unlist(gregexpr('_', name)), n=1)-1)
      }else{
        # If meta data
        name_sub <- name
      }
      
      data_result[data_result$var==name_sub,'varImp'] <-
        data_result[data_result$var==name_sub,'varImp'] +
        model$variable.importance[[j]]
      data_result[data_result$var==name_sub,'varImpCt'] <-
        data_result[data_result$var==name_sub,'varImpCt'] +
        1
    }
    
    ## Save model so I can reuse this RF
    RF[[i]] <- model
  }
  
  # Get mean gini decrease for nodes
  data_result$avgGini <- data_result$giniDec / trees#data_result$splits
  data_result$avgVI <- data_result$varImp / trees#data_result$varImpCt
  
  
  # Get Variable importance
  if(varImpPlot){
    avgGini <- ggplot() + 
      geom_point(aes(x=reorder(var,avgGini), y=avgGini/max(avgGini)), 
                 data=data_result) + 
      coord_flip() +
      xlab(NULL) +
      ylim(c(0,1)) +
      ylab('Gini') +
      theme_bw()
    avgVI <- ggplot() + 
      geom_point(aes(x=reorder(var,avgVI), y=avgVI/max(avgVI)), 
                 data=data_result) + 
      coord_flip() +
      xlab(NULL) +
      ylim(c(0,1)) +
      ylab('VarImp') +
      theme_bw()
    
    varImportance <- grid.arrange(avgGini,avgVI,
                                  layout_matrix = rbind(c(1,2)),
                                  bottom = 'Variable Importance - Percent of Max')
    
    return(list(RF, data_result, varImportance))
  }
  
  list(RF, data_result)
}

# Wrote function for
predict.myRandomForest <- function(model, data_pred, type='all', data=NULL){
  ## My implementation of prediction by the RF
  
  ## Function
  getMode <- function(x) {
    uniqx <- unique(na.omit(x))
    uniqx[which.max(tabulate(match(x, uniqx)))]
  }
  
  ## Prediction is the majority vote
  predictions <- data.frame(matrix(nrow=nrow(data_pred),
                                   ncol=length(model)))
  
  if(!is.null(data)){
    checkcols <- c()
    for(i in 1:ncol(data)){
      checkcols <- c(checkcols,
                     ifelse(is.character(data_pred[,i]) || 
                              is.factor(data_pred[,i]),i,NA))
    }
    checkcols <- checkcols[!is.na(checkcols)]
    if(length(checkcols)==0)
      data <- NULL # No need for data if nothing to check
  }
  data_pred_use <- data_pred
  drop_pred_rows <- rep(0, nrow(data_pred))
  
  for(i in 1:length(model)){
    # If original data given, this will allow verification
    if(!is.null(data)){
      # Reset
      data_pred_use <- data_pred
      
      # Get rows of data used
      data_row <- names(unlist(model[[i]]$where))
      unique_rows <- unique(gsub("\\..*","", data_row ))
      
      # Check each variable to see if contained in fit data
      #   Only check chars and factors (nums won't have this problem)
      drop_spot <- data.frame(matrix(0, 
                                     ncol=length(checkcols),
                                     nrow=nrow(data_pred_use))) # 0 keep, 1 drop
      for(j in checkcols){ 
        # If missing data, see if for single entry or all
        if(sum(!(data_pred_use[,j] %in% data[unique_rows,j]))>0){
          for(k in 1:nrow(data_pred_use)){
            #drop_spot[k,which(checkcols==j)] <- 
            #  !(data_pred_use[k,j] %in% data[unique_rows,j])
            if(!(data_pred_use[k,j] %in% data[unique_rows,j]))
              data_pred_use[k,j] <- NA
          }
        }
      }
      ## TODO:: Figure out how to drop only issue parts
      # Keep all pred_data with no issues
      #drop_pred_rows <- rowSums(drop_spot)>0
      #data_pred_use_drop <- data_pred_use[drop_pred_rows,]
      #data_pred_use_keep <- data_pred_use[!drop_pred_rows,]
      # Don't forget to shift predictions!
      
      # If all pred_data has issues, skip iteration
      if(nrow(data_pred_use)==0)
        next
      
    }
    
    pred_mat <- predict(model[[i]], data_pred_use)
    #pred_mat <- predict(model[[i]], data_pred_use_drop)
    
    # See if anything was drop
    if(nrow(pred_mat)!= length(drop_pred_rows)){
      stop('Error: This should not happen')
      # tmp <- data.frame(matrix(NA, ncol = ncol(pred_mat), nrow = nrow(data_pred)))
      # colnames(tmp) <- colnames(pred_mat)
      # rownames(tmp) <- rownames(data_pred)
      # 
      # jidx <- 0
      # for(j in 1:length(drop_pred_rows)){
      #   if(!drop_pred_rows[j]){# Kept
      #     jidx <- jidx + 1
      #     tmp[j,] <- pred_mat[jidx,]
      #   }
      # }
      # pred_mat <- tmp
    }
    
    sols <- colnames(pred_mat)
    for(j in 1:nrow(data_pred)){
      max_pred <- which.max(pred_mat[j,])
      if(length(max_pred)!=0)
        predictions[j,i] <- sols[max_pred]
    }
  }
  # Vote
  if(type=='pred'){
    modelPredictions <- apply(predictions, MARGIN = 1, function(x){ getMode(x) })
    
    return(modelPredictions) 
  }else if(type=='all'){
    modelPredictions <- apply(predictions, MARGIN = 1, function(x){ getMode(x) })
    modelPercents <- as.data.frame(t(apply(predictions, MARGIN = 1, 
                                           function(x,options){
                                             sapply(options, function(option,x){
                                               sum(x==option)/length(x)},
                                               x=x)
                                           }, options=unique(data_pred[,1]))))
    
    return(list('PredPerc'=cbind(modelPredictions,modelPercents),
                'Acc'=sum(modelPredictions==data_pred[,1])/nrow(data_pred),
                'ROC'=multiclass.roc(data_pred[,1], modelPercents)))
  }else if(type=='ROC'){
    
    modelPredictions <- as.data.frame(t(apply(predictions, MARGIN = 1, 
                                              function(x,options){
                                                sapply(options, function(option,x){
                                                  sum(x==option)/length(x)},
                                                  x=x)
                                              }, options=unique(data_pred[,1]))))
    
    return(list('Percents'=modelPredictions,
                'ROC'=multiclass.roc(data_pred[,1], modelPredictions),
                'Acc'=sum(modelPredictions[,1]==data_pred[,1])/nrow(data_pred)))
  }
}

# Inluded in RF to test
cheat.roc <- function(modelPredictions, modelPercents){
  ## My implementation of a fake ROC curve 
  
  options <- unique(modelPredictions)
  cols <- brewer.pal(length(options),'Set1')
  # 1 ROC curve, mock vs non mock
  roc.curve <- pROC::roc(ifelse(modelPredictions==options[1], 
                                options[1], paste0('Not-',options[1])),
                         as.numeric(modelPercents[,options[1]]),
                         levels=c(options[1], paste0('Not-',options[1])),
                         direction='>'
  )
  plot(roc.curve, col = cols[1])
  
  for(i in 2:length(options)){
    roc.curve <- pROC::roc(ifelse(modelPredictions==options[i], 
                                  options[i], paste0('Not-',options[i])), 
                           as.numeric(modelPercents[,options[i]]),
                           levels=c(options[i], paste0('Not-',options[i])),
                           direction='>')
    lines(roc.curve, col = cols[i])
  }
}

myRF_CV <- function(data, K=nrow(data), plotVarImp=T, 
                    #incFakeCell=T, incFakeMeta=T,
                    metaNames=NULL, 
                    cellData=NULL,
                    maxR=0.5, precR=1/100,
                    dropCols=c('Person'),
                    FakePairs=100, FakeMetas=100,
                    generalFakeCell=T,reduceEdge=0.05,
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
    if(is.null(maxR) || is.null(precR))
      stop('Error: Need maxR and precR to add fake cells')
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
        cellL <- rep(NA,length(KFunctions))
        cellR <- rep(NA,length(KFunctions))
        data_unique <- unique(cellData[,c('Image','Person','Stage')])
        fake_data <- NULL
        
        # Go over each unique K-Function 
        for(k in 1:length(KFunctions)){
          tmp <- unlist(strsplit(KFunctions[k],'_'))
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
        for(k in 1:length(KFunctions)){
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
    RF <- myRandomForest(data[-groups[[i]],], varImpPlot = F, 
                         metaNames=metaNames, ...) 
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

# Unused. May delete
importanceSplits <- function(fit){
  ## This is to determine variable importance
  
  # Based on https://github.com/cran/rpart/blob/master/R/importance.R
  
  ff <- fit$frame
  fpri <- which(ff$var != "<leaf>")  # points to primary splits in ff
  spri <- 1 + cumsum(c(0, 1 + ff$ncompete[fpri] + ff$nsurrogate[fpri]))
  spri <- spri[seq_along(fpri)] # points to primaries in the splits matrix
  nsurr <- ff$nsurrogate[fpri]  # number of surrogates each has
  
  sname <- vector("list", length(fpri))
  sval <- sname
  
  ## The importance for primary splits needs to be scaled
  ## It was a printout choice for the anova method to list % improvement in
  ##  the sum of squares, an importance calculation needs the total SS.
  ## All the other methods report an unscaled change.
  scaled.imp <- if (fit$method == "anova")
    fit$splits[spri, "improve"] * ff$dev[fpri]
  else fit$splits[spri, "improve"]
  
  sdim <- rownames(fit$splits)
  for (i in seq_along(fpri)) {
    ## points to surrogates
    if (nsurr[i] > 0L) {
      indx <- spri[i] + ff$ncompete[fpri[i]] + seq_len(nsurr[i])
      sname[[i]] <- sdim[indx]
      sval[[i]] <- scaled.imp[i] * fit$splits[indx, "adj"]
    }
  }
  
  import <- tapply(c(scaled.imp, unlist(sval)),
                   c(as.character(ff$var[fpri]), unlist(sname)),
                   sum)
  sort(c(import), decreasing = TRUE) # a named vector
}

#####
# Data Simulation
#####
# Likely unneeded
simulatePP_PCA_Meta <- function(stages=2,
                                peoplePerStage=40,
                                imagesPerPerson=1,
                                AB_Effect=TRUE,
                                age_Effect=TRUE,
                                laterStageVars=c(1/10000),
                                kappaA=20, kappaB=5,
                                nPCs=3, reduceEdge=0.05,
                                maxR=0.5, precR=1/100){
  ## Simulate a point process with some meta info
  
  ## Functions
  generateRandomAB <- function(stageName, reduceEdge,
                               peoplePerStage,
                               imagesPerPerson,
                               kappaA, kappaB,
                               cellA='A', cellB='B'){
    data <- NULL
    
    for(personCt in 1:peoplePerStage){
      for(imageCt in 1:imagesPerPerson){
        
        data_tmp <- generateRandomPP(xRange = c(0,1), yRange = c(0,1), 
                                     kappa = kappaA, percentAwayFromEdge = reduceEdge,
                                     cellType=cellA)
        
        data_Bpart <- generateRandomPP(xRange = c(0,1), yRange = c(0,1), 
                                       kappa = nrow(data_tmp)*kappaB, 
                                       percentAwayFromEdge = reduceEdge,
                                       cellType=cellB)
        
        data_tmp <- rbind(data_tmp, data_Bpart)
        data_tmp$Image <- imageCt+(personCt-1)*imagesPerPerson
        data_tmp$Person <- paste0('p',personCt)
        data_tmp$Stage <- stageName
        
        data <- rbind(data, data_tmp)
      }
    }
    
    data
  }
  
  ## Code
  
  ## Stage 0
  #   This stage randomly places A and B
  stage0 <- generateRandomAB(stageName = '0',
                             reduceEdge=reduceEdge,
                             peoplePerStage=peoplePerStage,
                             imagesPerPerson=imagesPerPerson,
                             kappaA=kappaA, kappaB=kappaB)
  
  ## Other Stages
  data_stages <- list()
  for(stage in 1:(stages-1)){
    if(!AB_Effect){
      # No AB Effect
      data_stages[[stage]] <- 
        generateRandomAB(stageName = as.character(stage),
                         reduceEdge=reduceEdge,
                         peoplePerStage=peoplePerStage,
                         imagesPerPerson=imagesPerPerson,
                         kappaA=kappaA, kappaB=kappaB)
    }else{
      # Some AB Effect
      data_stages[[stage]] <- generateReplicatedRandomPP(
        images = peoplePerStage*imagesPerPerson,
        kappa=kappaA, percentAwayFromEdge = reduceEdge, 
        ptDist = c('Poisson'), 
        ptPars = list(list('Lam'=kappaB)),
        placeDist=c('Normal'), 
        placePars = list(list('Var'=laterStageVars[stage])))
    }
    
    ## Align image, person, etc
    data_stages[[stage]]$Image <- 
      data_stages[[stage]]$Image + stage*(peoplePerStage*imagesPerPerson)
    ## This is only needed for the else part
    data_stages[[stage]]$Stage <- 
      as.character(stage)
    
    data_stages[[stage]]$Person <- NA
    for(j in 1:(peoplePerStage*imagesPerPerson)){
      data_stages[[stage]][data_stages[[stage]]$Image==j + stage*(peoplePerStage*imagesPerPerson),'Person'] <- 
        paste0('p',((j-1) %% peoplePerStage) + stage*peoplePerStage + 1)
    }
  }
  
  ## Combine Data
  data <- stage0
  for(stage in 1:(stages-1)){
    data <- rbind(data,data_stages[[stage]])
  }
  
  ## Get PCA
  data_pca <- getPCAData(data, nPCs=nPCs, maxR=maxR, precR=precR,
                         xRange = c(0,1), yRange = c(0,1),
                         cell_df = data.frame(
                           'cell1'=c('A','A','B'),
                           'cell2'=c('A','B','B')),
                         edgeCorrection=T, nbasis=21)
  data_rf <- getRFData(data_pca, dropNACol=F) 
  
  ## Add Meta
  data_rf$RandUnif <- runif(nrow(data_rf))
  data_rf$RandBin <- rbinom(nrow(data_rf),1,0.5)
  data_rf$Age <- rnorm(nrow(data_rf),40,5)
  
  if(age_Effect){
    for(stage in 0:(stages-1)){
      data_rf[data_rf$Stage==stage,'Age'] <-
        (stage*5+20) + 
        rbinom(nrow(data_rf[data_rf$Stage==stage,]),1,0.5)*5
    }
  }
  
  # Return
  list(data_rf,data)
}

# Likely unneeded
simulatePP_Mult_PCA_Meta <- function(peoplePerStage=40,
                                     imagesPerPerson=1,
                                     clusterEffect=TRUE,
                                     ageEffect=TRUE,
                                     laterStageVars=c(1/1000,1/10000),
                                     kappas=c(20,5,5),
                                     wasteKappas=c(10,15),
                                     nPCs=3, reduceEdge=0.05,
                                     maxR=0.5, precR=1/100){
  ## Simulate a pp with multiple outcomes
  
  ## Functions
  generateRandomPatterns <- function(stageName, reduceEdge,
                                     peoplePerStage,
                                     imagesPerPerson,
                                     kappas,
                                     cellTypes){
    
    ## Code
    data <- NULL
    
    for(personCt in 1:peoplePerStage){
      for(imageCt in 1:imagesPerPerson){
        data_tmp <- list()
        for(cell in 1:length(cellTypes)){
          if(cell==1){
            kapVal <- kappas[cell]
          }else{
            kapVal <- nrow(data_tmp[[cell-1]])*kappas[cell]
          }
          
          data_tmp[[cell]] <- generateRandomPP(xRange = c(0,1), yRange = c(0,1), 
                                               kappa = kapVal, 
                                               percentAwayFromEdge = reduceEdge,
                                               cellType=cellTypes[cell])
        }
        
        data_tmp_df <- convertList2DF(data_tmp)
        
        data_tmp_df$Image <- imageCt+(personCt-1)*imagesPerPerson
        data_tmp_df$Person <- paste0('p',personCt)
        data_tmp_df$Stage <- stageName
        
        data <- rbind(data, data_tmp_df)
      }
    }
    
    data
  }
  
  convertList2DF <- function(listData){
    data <- data.frame()
    for(i in 1:length(listData)){
      data <- rbind(data,listData[[i]])
    }
    
    data
  }
  
  ## Code
  if(length(kappas)!= (length(laterStageVars)+1))
    stop('Error: Make sure len(kappas)=len(Vars)+1')
  
  cellTypes <- LETTERS[1:length(kappas)]
  wasteCellTypes <- LETTERS[length(kappas) + 1:length(wasteKappas)]
  
  
  ## Stage 0
  #   This stage randomly places cells
  stage0 <- generateRandomPatterns(stageName = '0',
                                   reduceEdge=reduceEdge,
                                   peoplePerStage=peoplePerStage,
                                   imagesPerPerson=imagesPerPerson,
                                   kappas=kappas,
                                   cellTypes=cellTypes)
  
  stage0 <- rbind(stage0,
                  generateRandomPatterns(
                    stageName='0',
                    reduceEdge=reduceEdge,
                    peoplePerStage=peoplePerStage,
                    imagesPerPerson=imagesPerPerson,
                    kappas=wasteKappas,
                    cellTypes=wasteCellTypes))
  
  ## Other Stages
  data_stages <- list()
  stages <- 2 # Only two stage currently works
  for(stage in 1:(stages-1)){
    if(!clusterEffect){
      # No AB Effect
      data_stages[[stage]] <- 
        generateRandomPatterns(stageName = as.character(stage),
                               reduceEdge = reduceEdge,
                               peoplePerStage = peoplePerStage,
                               imagesPerPerson = imagesPerPerson,
                               kappas = kappas,
                               cellTypes = cellTypes)
    }else{
      # Some Clustering Effect
      ptDist <- rep('Poisson',length(laterStageVars))
      placeDist <- rep('Normal',length(laterStageVars))
      ptPars <- list()
      placePars <- list()
      for(i in 1:length(laterStageVars)){
        ptPars[[i]] <- list('Lam'=kappas[i+1])
        placePars[[i]] <- list('Var'=laterStageVars[i])
      }
      
      data_stages[[stage]] <- 
        generateReplicatedRandomPP(
          images = peoplePerStage*imagesPerPerson,
          kappa = kappas[1], percentAwayFromEdge = reduceEdge, 
          ptDist = ptDist, 
          ptPars = ptPars,
          placeDist = placeDist, 
          placePars = placePars)
    }
    
    data_stages[[stage]] <- rbind(data.frame(data_stages[[stage]],
                                             'Person'=NA,'Stage'=NA),
                                  generateRandomPatterns(
                                    stageName=as.character(stage),
                                    reduceEdge=reduceEdge,
                                    peoplePerStage=peoplePerStage,
                                    imagesPerPerson=imagesPerPerson,
                                    kappas=wasteKappas,
                                    cellTypes=wasteCellTypes))
    
    ## Align image, person, etc
    data_stages[[stage]]$Image <- 
      data_stages[[stage]]$Image + stage*(peoplePerStage*imagesPerPerson)
    ## This is only needed for the else part
    data_stages[[stage]]$Stage <- 
      as.character(stage)
    
    data_stages[[stage]]$Person <- NA
    for(j in 1:(peoplePerStage*imagesPerPerson)){
      data_stages[[stage]][data_stages[[stage]]$Image==j + stage*(peoplePerStage*imagesPerPerson),'Person'] <- 
        paste0('p',((j-1) %% peoplePerStage) + stage*peoplePerStage + 1)
    }
  }
  
  ## Combine Data
  data <- rbind(stage0, convertList2DF(data_stages))
  
  ## Get PCA
  allCellTypes <- c(cellTypes, wasteCellTypes)
  cell_df <- data.frame('cell1'=rep(NA,choose(length(allCellTypes),2)), 
                        'cell2'=NA)
  idx <- 1
  for(i in 1:length(allCellTypes)){
    for(j in i:length(allCellTypes)){
      cell_df[idx,] <- c(allCellTypes[i],
                         allCellTypes[j]) 
      idx <- idx +1
    }
  }
  
  data_pca <- getPCAData(data, nPCs=nPCs, maxR=maxR, precR=precR,
                         xRange = c(0,1), yRange = c(0,1),
                         cell_df = cell_df,
                         edgeCorrection=T, nbasis=21)
  data_rf <- getRFData(data_pca, dropNACol=F) 
  
  ## Add Meta
  data_rf$RandUnif <- runif(nrow(data_rf))
  data_rf$RandBin <- rbinom(nrow(data_rf),1,0.5)
  data_rf$Age <- rnorm(nrow(data_rf),40,5)
  
  if(ageEffect){
    for(stage in 0:(stages-1)){
      data_rf[data_rf$Stage==stage,'Age'] <-
        (stage*5+20) + 
        rbinom(nrow(data_rf[data_rf$Stage==stage,]),1,0.5)*5
    }
  }
  
  # Return
  list(data_rf,data)
}


simulatePP_Mult_PCA_Meta <- function(cellVarData=
                                       data.frame('stage'=c(0,1,2),
                                                  'A'=c(0,0,0),
                                                  'B'=c(0,1/1000,1/1000),
                                                  'C'=c(1/10000,1/1000,1/500),
                                                  'D'=c(1/100,1/100,1/100),
                                                  'E'=c(0,0,0),
                                                  'F'=c(1/250,1/250,1/250)),
                                     cellKappaData=data.frame(
                                       'cell'=c('A','B','C','D','E','F'),
                                       'clusterCell'=c(NA,'A','B','C',NA,'A'),
                                       'kappa'=c(20,5,4,2,15,5)),
                                     peoplePerStage=40,
                                     imagesPerPerson=1,
                                     ageEffect=TRUE,
                                     nPCs=3, reduceEdge=0.05,
                                     maxR=0.5, precR=1/100){
  ## A very comprehensive PP generation method
  
  ## Functions
  generateRandomPatterns <- function(stageName, reduceEdge,
                                     peoplePerStage,
                                     imagesPerPerson,
                                     kappas,
                                     cellTypes,
                                     kappaSep=F){
    
    ## Code
    data <- NULL
    
    for(personCt in 1:peoplePerStage){
      for(imageCt in 1:imagesPerPerson){
        data_tmp <- list()
        for(cell in 1:length(cellTypes)){
          if(cell==1 || kappaSep){
            kapVal <- kappas[cell]
          }else{
            kapVal <- nrow(data_tmp[[cell-1]])*kappas[cell]
          }
          
          data_tmp[[cell]] <- generateRandomPP(xRange = c(0,1), yRange = c(0,1), 
                                               kappa = kapVal, 
                                               percentAwayFromEdge = reduceEdge,
                                               cellType=cellTypes[cell])
        }
        
        data_tmp_df <- convertList2DF(data_tmp)
        
        data_tmp_df$Image <- imageCt+(personCt-1)*imagesPerPerson
        data_tmp_df$Person <- paste0('p',personCt)
        data_tmp_df$Stage <- stageName
        
        data <- rbind(data, data_tmp_df)
      }
    }
    
    data
  }
  
  convertList2DF <- function(listData){
    data <- data.frame()
    for(i in 1:length(listData)){
      data <- rbind(data,listData[[i]])
    }
    
    data
  }
  
  clusterAroundCells <- function(clusterData, stageName, reduceEdge,
                                 cells, clusterCells, kappas,
                                 minPts=1){
    ## Function
    placePts <- function(currXY, cell, numPts, varValue, reduceEdge){
      if(numPts==0)
        return()
      xRange <- c(0,1)
      yRange <- xRange
      
      data_ret <- data.frame('x'=rep(NA,numPts), 'y'=NA, 'cellType'=cell)
      
      compPts <- 0
      while(compPts < numPts){
        data_ret[compPts+1,'x'] <- rnorm(1,mean=currXY[1],sd=sqrt(varValue))
        data_ret[compPts+1,'y'] <- rnorm(1,mean=currXY[2],sd=sqrt(varValue))
        
        # Ensure its in boundaries
        if(data_ret[compPts+1,'x'] >= xRange[1]+reduceEdge &
           data_ret[compPts+1,'x'] <= xRange[2]-reduceEdge & 
           data_ret[compPts+1,'y'] >= yRange[1]+reduceEdge & 
           data_ret[compPts+1,'y'] <= yRange[2]-reduceEdge){
          compPts <- compPts + 1
        }
      }
      
      data_ret
    }
    
    ## Code
    newData <- data.frame()
    # Go through each cell
    for(i in 1:length(cells)){
      clusterCellData <- clusterData[clusterData$cellType==clusterCells[i],]
      
      # Go through each cell and develop clusters
      for(j in 1:nrow(clusterCellData)){
        data_pts <- placePts(currXY=as.numeric(clusterCellData[j, c('x','y')]), 
                             cell=cells[i], 
                             numPts=rpois(1,kappas[i]), 
                             varValue=cellVarData[cellVarData$stage==stageName,
                                                  cells[i]],
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
        
        data_pts <- placePts(currXY=as.numeric(clusterCellData[j, c('x','y')]), 
                             cell=cells[i], 
                             numPts=rpois(1,kappas[i]), 
                             varValue=cellVarData[cellVarData$stage==stageName,
                                                  cells[i]],
                             reduceEdge=reduceEdge)
        data_pts$Image <- clusterCellData[j,'Image']
        data_pts$Person <- clusterCellData[j,'Person']
        data_pts$Stage <- clusterCellData[j,'Stage']
        
        newData <- rbind(newData, data_pts)
      }
      
    }
    
    newData
  }
  
  ## Code
  
  ## Cell Data
  data_stages <- list()
  for(stageIdx in 1:nrow(cellVarData)){
    
    stage <- cellVarData$stage[stageIdx]
    
    ## Do all non-clustering first
    nextCellIdx <- which(is.na(cellKappaData$clusterCell))
    
    data_stages[[stageIdx]] <- 
      generateRandomPatterns(stageName=stage, 
                             reduceEdge=reduceEdge,
                             peoplePerStage=peoplePerStage,
                             imagesPerPerson=imagesPerPerson,
                             kappas=cellKappaData[nextCellIdx,'kappa'],
                             cellTypes=cellKappaData[nextCellIdx,'cell'],
                             kappaSep=T)
    # ## Align image, person, etc
    data_stages[[stageIdx]]$Image <-
      data_stages[[stageIdx]]$Image + (stageIdx-1)*(peoplePerStage*imagesPerPerson)
    ## This is only needed for the clustered part
    #data_stages[[stageIdx]]$Stage <- as.character(stage)
    
    data_stages[[stageIdx]]$Person <- NA
    for(j in 1:(peoplePerStage*imagesPerPerson)){
      data_stages[[stageIdx]][data_stages[[stageIdx]]$Image==j + (stageIdx-1)*(peoplePerStage*imagesPerPerson),'Person'] <-
        paste0('p',((j-1) %% peoplePerStage) + (stageIdx-1)*peoplePerStage + 1)
    }
    
    # Record generation
    completeCells <- c(cellKappaData[nextCellIdx,'cell'])
    
    
    
    ## Recursively plot clusters
    while(length(completeCells)!= nrow(cellKappaData)){
      
      # See all cells that cluster around newly added (not previously added)
      nextCellIdx <- which(cellKappaData$clusterCell %in% completeCells &
                             !(cellKappaData$cell %in% completeCells))
      data_stages[[stageIdx]] <- rbind(data_stages[[stageIdx]],
                                       clusterAroundCells(
                                         clusterData=data_stages[[stageIdx]][
                                           data_stages[[stageIdx]]$cellType %in% 
                                             cellKappaData[nextCellIdx,'clusterCell'],],
                                         stageName=cellVarData[stageIdx,'stage'], 
                                         reduceEdge=reduceEdge,
                                         #peoplePerStage=peoplePerStage,
                                         #imagesPerPerson=imagesPerPerson,
                                         cells=cellKappaData[nextCellIdx,'cell'],
                                         clusterCells=cellKappaData[nextCellIdx,'clusterCell'],
                                         kappas=cellKappaData[nextCellIdx,'kappa']))
      
      # Record generation
      completeCells <- c(completeCells, cellKappaData[nextCellIdx,'cell'])
    }
    
  }
  
  ## Combine Data
  data <- convertList2DF(data_stages)
  
  ## Get PCA
  cell_df <- data.frame('cell1'=rep(NA,choose(nrow(cellKappaData),2)), 
                        'cell2'=NA)
  idx <- 1
  for(i in 1:nrow(cellKappaData)){
    for(j in i:nrow(cellKappaData)){
      cell_df[idx,] <- c(cellKappaData$cell[i],
                         cellKappaData$cell[j]) 
      idx <- idx +1
    }
  }
  
  data_pca <- getPCAData(data, nPCs=nPCs, maxR=maxR, precR=precR,
                         xRange = c(0,1), yRange = c(0,1),
                         cell_df = cell_df,
                         edgeCorrection=T, nbasis=21)
  data_rf <- getRFData(data_pca, dropNACol=F) 
  
  ## Add Meta
  data_rf$RandUnif <- runif(nrow(data_rf))
  data_rf$RandBin <- rbinom(nrow(data_rf),1,0.5)
  data_rf$Age <- rnorm(nrow(data_rf),40,5)
  
  if(ageEffect){
    for(stageIdx in 1:length(cellVarData$stage)){
      stage <- cellVarData$stage[stageIdx]
      data_rf[data_rf$Stage==stage,'Age'] <-
        (stageIdx*5+20) + 
        rbinom(nrow(data_rf[data_rf$Stage==stage,]),1,0.5)*5
    }
  }
  
  # Return
  list(data_rf,data)
}
