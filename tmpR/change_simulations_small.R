nSims <- 100
changes <- c(1.5/10, 2/10, 4/10, 8/10, 16/10,
             1/10,
             1/15, 1/20, 1/40, 1/80, 1/160)
baseline <- changes[6]
results <- list()

cells <- paste0('c',1:4)
cells_interactions <- rbind(data.frame(t(combn(cells,2))),
                            data.frame('X1'=cells,'X2'=cells))


for(c in 1:length(changes)){
  cat(paste0('\n- Change: ',changes[c],'\n-- (',nSims,'): '))

  VarVI <- VarLVI <-
    VarVI_r <- VarLVI_r <-
    VarVI_o <- VarLVI_o <-
    VarVI_max <- VarVI_o_max <-
    VarLVI_max <- VarLVI_o_max <-
    data.frame('var'=c('age',
                       paste(cells_interactions$X1,cells_interactions$X2,sep='_'),
                       'gender'))

  set.seed(c*12345)

  for(i in 1:nSims){
    cat(paste0(i,', '))
    # Generate
    dat <- simulatePP(cellVarData=
                        data.frame('stage'=c(0,1),
                                   'c1'=c(0,0),
                                   'c2'=c(baseline,changes[c]),
                                   'c3'=c(1/50,1/50),
                                   'c4'=c(1/10,1/10)),
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
                          agents_df = cells_interactions,
                          silent=TRUE)
    pcaMeta <- simulateMeta(pcaData,
                            metaInfo = data.frame(
                              'var'=c('gender','age'),
                              'rdist'=c('rbinom','rnorm'),
                              'Stage_0'=c('0.5','25'),
                              'Stage_1'=c('0.5','25')))
    # Fit RF
    rfcv <- funkyModel(data=pcaMeta,
                       outcome = 'Stage', unit = 'Person',
                       metaNames=c('gender','age'), silent = TRUE)

    # Org Data
    tmp <- rfcv$VariableImportance[,c('var','est','sd')]
    tmp <- tmp[order(-tmp$est),]
    tmp$lower1sd <- tmp$est-tmp$sd

    ## Orange/Red Information

    # Above both lines
    tmp$IndCurve <- (tmp$est > rfcv$InterpolationCutoff)
    tmp$IndNoise <- (tmp$est > rfcv$NoiseCutoff)
    tmp[[paste0('iter',i)]] <- tmp$IndCurve * tmp$IndNoise
    VarVI <- merge(VarVI,tmp[,c('var',paste0('iter',i))])
    tmp[[paste0('iter',i)]] <- tmp$IndNoise
    VarVI_r <- merge(VarVI_r,tmp[,c('var',paste0('iter',i))])
    tmp[[paste0('iter',i)]] <- tmp$IndCurve
    VarVI_o <- merge(VarVI_o,tmp[,c('var',paste0('iter',i))])

    # max Orange
    tmp$IndCurve <- (tmp$est > max(rfcv$InterpolationCutoff))
    tmp[[paste0('iter',i)]] <- tmp$IndCurve * tmp$IndNoise
    VarVI_max <- merge(VarVI_max,tmp[,c('var',paste0('iter',i))])
    tmp[[paste0('iter',i)]] <- tmp$IndCurve
    VarVI_o_max <- merge(VarVI_o_max,tmp[,c('var',paste0('iter',i))])



    # 1 SD above both lines
    tmp$SDIndCurve <- (tmp$lower1sd > rfcv$InterpolationCutoff)
    tmp$SDIndNoise <- (tmp$lower1sd > rfcv$NoiseCutoff)
    tmp[[paste0('iter',i)]] <- tmp$SDIndCurve * tmp$SDIndNoise
    VarLVI <- merge(VarLVI,tmp[,c('var',paste0('iter',i))])
    tmp[[paste0('iter',i)]] <- tmp$SDIndNoise
    VarLVI_r <- merge(VarLVI_r,tmp[,c('var',paste0('iter',i))])
    tmp[[paste0('iter',i)]] <- tmp$SDIndCurve
    VarLVI_o <- merge(VarLVI_o,tmp[,c('var',paste0('iter',i))])

    # max Orange
    tmp$IndCurve <- (tmp$lower1sd > max(rfcv$InterpolationCutoff))
    tmp[[paste0('iter',i)]] <- tmp$IndCurve * tmp$IndNoise
    VarLVI_max <- merge(VarLVI_max,tmp[,c('var',paste0('iter',i))])
    tmp[[paste0('iter',i)]] <- tmp$IndCurve
    VarLVI_o_max <- merge(VarLVI_o_max,tmp[,c('var',paste0('iter',i))])
  }

  results[[paste0('Change_',changes[c])]] <-
    list('VarVI'=VarVI, 'VarLVI'=VarLVI,
         'VarVI_r'=VarVI_r, 'VarLVI_r'=VarLVI_r,
         'VarVI_o'=VarVI_o, 'VarLVI_o'=VarLVI_o,
         'VarVI_max'=VarVI_max, 'VarVI_o_max'=VarVI_o_max,
         'VarLVI_max'=VarLVI_max, 'VarLVI_o_max'=VarLVI_o_max)


  ## Save RDS
  saveRDS(results, paste0('C:/Users/j53vande/Downloads/change_sim_small.rds'))
}


####################################
# Look at c1_c2
data_summary <- data.frame('var'=changes,
                           'vert'=NA, 'lvert'=NA,
                           'curve'=NA, 'lcurve'=NA,
                           'mcurve'=NA, 'lmcurve'=NA,
                           'both'=NA, 'lboth'=NA,
                           'mboth'=NA, 'lmboth'=NA)

for(i in 1:length(results)){
  # (vert)
  data_summary[i,'vert'] <-
    sum(results[[i]]$VarVI_r[results[[i]]$VarVI_r$var=='c1_c2',-1]) /
    length(results[[i]]$VarVI_r[results[[i]]$VarVI_r$var=='c1_c2',-1])
  data_summary[i,'lvert'] <-
    sum(results[[i]]$VarLVI_r[results[[i]]$VarLVI_r$var=='c1_c2',-1]) /
    length(results[[i]]$VarLVI_r[results[[i]]$VarLVI_r$var=='c1_c2',-1])

  # (curve)
  data_summary[i,'curve'] <-
    sum(results[[i]]$VarVI_o[results[[i]]$VarVI_o$var=='c1_c2',-1]) /
    length(results[[i]]$VarVI_o[results[[i]]$VarVI_o$var=='c1_c2',-1])
  data_summary[i,'lcurve'] <-
    sum(results[[i]]$VarLVI_o[results[[i]]$VarLVI_o$var=='c1_c2',-1]) /
    length(results[[i]]$VarLVI_o[results[[i]]$VarLVI_o$var=='c1_c2',-1])

  # (max curve)
  data_summary[i,'mcurve'] <-
    sum(results[[i]]$VarVI_o_max[results[[i]]$VarVI_o_max$var=='c1_c2',-1]) /
    length(results[[i]]$VarVI_o_max[results[[i]]$VarVI_o_max$var=='c1_c2',-1])
  data_summary[i,'lmcurve'] <-
    sum(results[[i]]$VarLVI_o_max[results[[i]]$VarLVI_o_max$var=='c1_c2',-1]) /
    length(results[[i]]$VarLVI_o_max[results[[i]]$VarLVI_o_max$var=='c1_c2',-1])

  # (either)
  data_summary[i,'both'] <-
    sum(results[[i]]$VarVI[results[[i]]$VarVI$var=='c1_c2',-1]) /
    length(results[[i]]$VarVI[results[[i]]$VarVI$var=='c1_c2',-1])
  data_summary[i,'lboth'] <-
    sum(results[[i]]$VarLVI[results[[i]]$VarLVI$var=='c1_c2',-1]) /
    length(results[[i]]$VarLVI[results[[i]]$VarLVI$var=='c1_c2',-1])

  # (max either)
  data_summary[i,'mboth'] <-
    sum(results[[i]]$VarVI_max[results[[i]]$VarVI_max$var=='c1_c2',-1]) /
    length(results[[i]]$VarVI_max[results[[i]]$VarVI_max$var=='c1_c2',-1])
  data_summary[i,'lmboth'] <-
    sum(results[[i]]$VarLVI_max[results[[i]]$VarLVI_max$var=='c1_c2',-1]) /
    length(results[[i]]$VarLVI_max[results[[i]]$VarLVI_max$var=='c1_c2',-1])


  # # See how red/orange does in general
  # # (vert)
  # sum(results[[i]]$VarVI_r[-1]) /
  #   (ncol(results[[i]]$VarVI_r[-1])*nrow(results[[i]]$VarVI_r[-1]))
  # sum(results[[i]]$VarLVI_r[-1]) /
  #   (ncol(results[[i]]$VarLVI_r[-1])*nrow(results[[i]]$VarLVI_r[-1]))
  #
  # # (curve)
  # sum(results[[i]]$VarVI_o[-1]) /
  #   (ncol(results[[i]]$VarVI_o[-1])*nrow(results[[i]]$VarVI_o[-1]))
  # sum(results[[i]]$VarLVI_o[-1]) /
  #   (ncol(results[[i]]$VarLVI_o[-1])*nrow(results[[i]]$VarLVI_o[-1]))
  #
  # # (either)
  # sum(results[[i]]$VarVI[-1]) /
  #   (ncol(results[[i]]$VarVI[-1])*nrow(results[[i]]$VarVI[-1]))
  # sum(results[[i]]$VarLVI[-1]) /
  #   (ncol(results[[i]]$VarLVI[-1])*nrow(results[[i]]$VarLVI[-1]))
}


data_plot <-
  tidyr::pivot_longer(data_summary,cols=vert:lmboth)

ggplot2::ggplot(data_plot[data_plot$name %in% c('vert','curve','mcurve','both'),],
                ggplot2::aes(x=var, y=value, color=name, group=name)) +
  ggplot2::geom_line(linewidth=1.25) +
  ggplot2::geom_point(size=3) +
  ggplot2::geom_vline(ggplot2::aes(xintercept=baseline),
                      color='black', linetype='dashed', linewidth=1.25) +
  ggplot2::geom_hline(ggplot2::aes(yintercept=0.05), linetype='dotted') +
  ggplot2::theme_bw() +
  ggplot2::theme(axis.text=ggplot2::element_text(size=18),
                 axis.title = ggplot2::element_text(size=22),
                 legend.position = "none",
                 legend.title = ggplot2::element_text(size=22),
                 legend.text = ggplot2::element_text(size=18)) +
  ggplot2::scale_color_discrete(name='Cutoff',
                                labels=c('Both','Interp',
                                         'Max Interp','Noise')) +
  ggplot2::xlab('Variance') +
  ggplot2::ylab(NULL) +
  ggplot2::scale_x_log10()
