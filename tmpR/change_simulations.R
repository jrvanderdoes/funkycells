nSims <- 100
changes <- c(1/25,1/50,1/100)
results <- list()

set.seed(12345)
for(c in 1:length(changes)){
  cat(paste0('\n- Change: ',changes[c],'\n-- (',nSims,'): '))

  VarVI <- VarLVI <- data.frame('var'=c('age','c1_c1','c1_c2','c1_c3',
                                        'c1_c4','c1_c5','c2_c1','c2_c2',
                                        'c2_c3','c2_c4','c2_c5','c3_c1',
                                        'c3_c2','c3_c3','c3_c4','c3_c5',
                                        'c4_c1','c4_c2','c4_c3','c4_c4',
                                        'c4_c5','c5_c1','c5_c2','c5_c3',
                                        'c5_c4','c5_c5','gender'))
  RedVI1 <- RedLVI1 <- OrangeVI1 <- OrangeLVI1 <-
    RedVI2 <- RedLVI2 <- OrangeVI2 <- OrangeLVI2 <- rep(NA,nSims)
  for(i in 1:nSims){
    cat(paste0(i,', '))
    # Generate
    dat <- simulatePP(cellVarData=
                        data.frame('stage'=c(0,1),
                                   'c1'=c(0,0),
                                   'c2'=c(1/25,changes[c]), 'c3'=c(1/50,1/50),
                                   'c4'=c(0,0)),
                      cellKappaData=data.frame(
                        'cell'=paste0('c',1:4),
                        'clusterCell'=c(NA,'c1','c1', NA),
                        'kappa'=c(rbinom(1,50,0.5),
                                  rbinom(2,20,0.5),
                                  rbinom(1,50,0.5))),
                      peoplePerStage=17,
                      imagesPerPerson=1,
                      silent=T)
    pcaData <- getPCAData(data = dat,repeatedUniqueId='Image',
                          xRange = c(0,1), yRange = c(0,1),
                          silent=T)
    pcaMeta <- simulateMeta(pcaData,
                            metaInfo = data.frame(
                              'var'=c('gender','age'),
                              'rdist'=c('rbinom','rnorm'),
                              'Stage_0'=c('0.5','25'),
                              'Stage_1'=c('0.5','25')))
    # Fit RF
    rfcv <- computeRandomForest_CVPC_Permute(data=pcaMeta,
                                             outcome = 'Stage',
                                             unit = 'Person',
                                             metaNames=c('gender','age'),
                                             silent = T)

    # Org Data
    tmp <- rfcv$VariableImportance[,c('var','est','sd')]
    tmp <- tmp[order(-tmp$est),]
    tmp$lower1sd <- tmp$est-tmp$sd

    # Orange/Red for one c1_c2
    idx <- which(tmp$var=='c1_c2')
    RedVI1[i] <- tmp[idx,'est']> rfcv$NoiseCutoff
    RedLVI1[i] <- tmp[idx,'lower1sd']> rfcv$NoiseCutoff
    OrangeVI1[i] <- tmp[idx,'est']> rfcv$InterpolationCutoff[idx]
    OrangeLVI1[i] <- tmp[idx,'lower1sd']> rfcv$InterpolationCutoff[idx]

    # Orange/Red for one c2_c1
    idx <- which(tmp$var=='c2_c1')
    RedVI2[i] <- tmp[idx,'est']> rfcv$NoiseCutoff
    RedLVI2[i] <- tmp[idx,'lower1sd']> rfcv$NoiseCutoff
    OrangeVI2[i] <- tmp[idx,'est']> rfcv$InterpolationCutoff[idx]
    OrangeLVI2[i] <- tmp[idx,'lower1sd']> rfcv$InterpolationCutoff[idx]

    # Above both lines
    tmp$IndCurve <- (tmp$est > rfcv$InterpolationCutoff)
    tmp$IndNoise <- (tmp$est > rfcv$NoiseCutoff)
    tmp[[paste0('iter',i)]] <- tmp$IndCurve * tmp$IndNoise
    VarVI <- merge(VarVI,tmp[,c('var',paste0('iter',i))])
    tmp[[paste0('iter',i)]] <- tmp$IndNoise
    VarVI_r <- merge(VarVI,tmp[,c('var',paste0('iter',i))])
    tmp[[paste0('iter',i)]] <- tmp$IndCurve
    VarVI_o <- merge(VarVI,tmp[,c('var',paste0('iter',i))])

    # 1 SD above both lines
    tmp$SDIndCurve <- (tmp$lower1sd > rfcv$InterpolationCutoff)
    tmp$SDIndNoise <- (tmp$lower1sd > rfcv$NoiseCutoff)
    tmp[[paste0('iter',i)]] <- tmp$SDIndCurve * tmp$SDIndNoise
    VarLVI <- merge(VarLVI,tmp[,c('var',paste0('iter',i))])
    tmp[[paste0('iter',i)]] <- tmp$SDIndNoise
    VarLVI_r <- merge(VarLVI,tmp[,c('var',paste0('iter',i))])
    tmp[[paste0('iter',i)]] <- tmp$SDIndCurve
    VarLVI_o <- merge(VarLVI,tmp[,c('var',paste0('iter',i))])
  }

  results[[paste0('Change_',changes[c])]] <-
            list('VarVI'=VarVI, 'VarLVI'=VarLVI,
                 'VarVI_r'=VarVI_r, 'VarLVI_r'=VarLVI_r,
                 'VarVI_o'=VarVI_o, 'VarLVI_o'=VarLVI_o,
                 'RedVI1'=RedVI1, 'RedLVI1'=RedLVI1,
                    'OrangeVI1'=OrangeVI1, 'OrangeLVI1'=OrangeLVI1,
                 'RedVI2'=RedVI2, 'RedLVI2'=RedLVI2,
                    'OrangeVI2'=OrangeVI2, 'OrangeLVI2'=OrangeLVI2)
}

## Save RDS
saveRDS(results,
        paste0('C:/Users/jerem/OneDrive/Documents/School/Waterloo/',
               'Research/RPackages/funkycells/tmpR/small_simulation.rds'))

# Look at c1_c2 / c2_c1
for(i in 1:length(results)){
  # See how red/orange does for c1_c2
  results[[i]]$RedVI1 / length(results[[i]]$RedVI1)
  results[[i]]$RedLVI1 / length(results[[i]]$RedLVI1)
  results[[i]]$OrangeVI1 / length(results[[i]]$OrangeVI1)
  results[[i]]$OrangeLVI1 / length(results[[i]]$OrangeLVI1)

  # See how red/orange does for c2_c1
  results[[i]]$RedVI2 / length(results[[i]]$RedVI2)
  results[[i]]$RedLVI2 / length(results[[i]]$RedLVI2)
  results[[i]]$OrangeVI2 / length(results[[i]]$OrangeVI2)
  results[[i]]$OrangeLVI2 / length(results[[i]]$OrangeLVI2)

  # Red
  sum(results[[i]]$VarVI_r[results[[i]]$VarVI_r$var=='c1_c2',-1]) /
    length(results[[i]]$VarVI_r[results[[i]]$VarVI_r$var=='c1_c2',-1])
  sum(results[[i]]$VarLVI_r[results[[i]]$VarLVI_r$var=='c1_c2',-1]) /
    length(results[[i]]$VarLVI_r[results[[i]]$VarLVI_r$var=='c1_c2',-1])

  sum(results[[i]]$VarVI_r[results[[i]]$VarVI_r$var=='c2_c1',-1]) /
    length(results[[i]]$VarVI_r[results[[i]]$VarVI_r$var=='c2_c1',-1])
  sum(results[[i]]$VarLVI_r[results[[i]]$VarLVI_r$var=='c2_c1',-1]) /
    length(results[[i]]$VarLVI_r[results[[i]]$VarLVI_r$var=='c2_c1',-1])


  # Orange
  sum(results[[i]]$VarVI_o[results[[i]]$VarVI_o$var=='c1_c2',-1]) /
    length(results[[i]]$VarVI_o[results[[i]]$VarVI_o$var=='c1_c2',-1])
  sum(results[[i]]$VarLVI_o[results[[i]]$VarLVI_o$var=='c1_c2',-1]) /
    length(results[[i]]$VarLVI_o[results[[i]]$VarLVI_o$var=='c1_c2',-1])

  sum(results[[i]]$VarVI_o[results[[i]]$VarVI_o$var=='c2_c1',-1]) /
    length(results[[i]]$VarVI_o[results[[i]]$VarVI_o$var=='c2_c1',-1])
  sum(results[[i]]$VarLVI_o[results[[i]]$VarLVI_o$var=='c2_c1',-1]) /
    length(results[[i]]$VarLVI_o[results[[i]]$VarLVI_o$var=='c2_c1',-1])

  # Either
  sum(results[[i]]$VarVI[results[[i]]$VarVI$var=='c1_c2',-1]) /
    length(results[[i]]$VarVI[results[[i]]$VarVI$var=='c1_c2',-1])
  sum(results[[i]]$VarLVI[results[[i]]$VarLVI$var=='c1_c2',-1]) /
    length(results[[i]]$VarLVI[results[[i]]$VarLVI$var=='c1_c2',-1])

  sum(results[[i]]$VarVI[results[[i]]$VarVI$var=='c2_c1',-1]) /
    length(results[[i]]$VarVI[results[[i]]$VarVI$var=='c2_c1',-1])
  sum(results[[i]]$VarLVI[results[[i]]$VarLVI$var=='c2_c1',-1]) /
    length(results[[i]]$VarLVI[results[[i]]$VarLVI$var=='c2_c1',-1])
}

# Look Null red rejection
sum(results$Change_0$VarVI_r[-1]) /
  (ncol(results$Change_0$VarVI_r[-1])*nrow(results$Change_0$VarVI_r[-1]))
sum(results$Change_0$VarLVI_r[-1]) /
  (ncol(results$Change_0$VarLVI_r[-1])*nrow(results$Change_0$VarLVI_r[-1]))

# Look Null orange rejection
sum(results$Change_0$VarVI_o[-1]) /
  (ncol(results$Change_0$VarVI_o[-1])*nrow(results$Change_0$VarVI_o[-1]))
sum(results$Change_0$VarLVI_o[-1]) /
  (ncol(results$Change_0$VarLVI_o[-1])*nrow(results$Change_0$VarLVI_o[-1]))

# Look Null both rejection
sum(results$Change_0$VarVI[-1]) /
  (ncol(results$Change_0$VarVI[-1])*nrow(results$Change_0$VarVI[-1]))
sum(results$Change_0$VarLVI[-1]) /
  (ncol(results$Change_0$VarLVI[-1])*nrow(results$Change_0$VarLVI[-1]))
