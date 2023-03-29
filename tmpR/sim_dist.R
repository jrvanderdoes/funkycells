set.seed(12345)
n <- 200
trials <- 1000
evalPts <- seq(0,1,0.05)

results_p <- results_loc <-
  data.frame('MST5_OL2_ori'=1:trials, 'MST5_OL2_wei'=NA,
             'MST5_OL2_max'=NA, 'MST5_OL2_gen'=NA,
             'MST5_OL1_ori'=NA, 'MST5_OL1_wei'=NA,
             'MST5_OL1_max'=NA, 'MST5_OL1_gen'=NA,
             'MST5_OMax_ori'=NA, 'MST5_OMax_wei'=NA,
             'MST5_OMax_max'=NA, 'MST5_OMax_gen'=NA,

             'MST5_CEL2_ori'=1:trials, 'MST5_CEL2_wei'=NA,
             'MST5_CEL2_max'=NA, 'MST5_CEL2_gen'=NA,
             'MST5_CEL1_ori'=NA, 'MST5_CEL1_wei'=NA,
             'MST5_CEL1_max'=NA, 'MST5_CEL1_gen'=NA,
             'MST5_CEMax_ori'=NA, 'MST5_CEMax_wei'=NA,
             'MST5_CEMax_max'=NA, 'MST5_CEMax_gen'=NA,

             'meanChange'=NA
  )
data_cp <- list()
for(i in 1:trials){
  cat(paste0(i,'\n'))
  ## Generate Data
  data_cp[[i]] <- generate_data_fd(
    ns = c(n),
    eigsList = list(c(3,2,1,0.5)),
    basesList = list(create.bspline.basis(nbasis=4, norder=4)),
    meansList = c(0),
    distsArray = c('Normal'),
    evals = evalPts,
    kappasArray = c(0),
    silent = T)
  data_dist0 <- calculateDistanceMatrix(data=data_cp[[i]], silent=T,
                                       errType='L2', dataUse='Orig')
  data_dist1 <- calculateDistanceMatrix(data=data_cp[[i]], silent=T,
                                        errType='L1', dataUse='Orig')
  data_dist2 <- calculateDistanceMatrix(data=data_cp[[i]], silent=T,
                                        errType='Max', dataUse='Orig')

  data_dist3 <- calculateDistanceMatrix(data=data_cp[[i]], silent=T,
                                       errType='L2', dataUse='CE')
  data_dist4 <- calculateDistanceMatrix(data=data_cp[[i]], silent=T,
                                        errType='L1', dataUse='CE')
  data_dist5 <- calculateDistanceMatrix(data=data_cp[[i]], silent=T,
                                        errType='Max', dataUse='CE')


  MST5_0 <- mstree(as.dist(data_dist0), ngmax=5)
  gSeg_MST5_0 <- gseg1(n, MST5_0, pval.perm = T, pval.appr = F)
  MST5_1 <- mstree(as.dist(data_dist1), ngmax=5)
  gSeg_MST5_1 <- gseg1(n, MST5_1, pval.perm = T, pval.appr = F)
  MST5_2 <- mstree(as.dist(data_dist2), ngmax=5)
  gSeg_MST5_2 <- gseg1(n, MST5_2, pval.perm = T, pval.appr = F)

  MST5_3 <- mstree(as.dist(data_dist3), ngmax=5)
  gSeg_MST5_3 <- gseg1(n, MST5_3, pval.perm = T, pval.appr = F)
  MST5_4 <- mstree(as.dist(data_dist4), ngmax=5)
  gSeg_MST5_4 <- gseg1(n, MST5_4, pval.perm = T, pval.appr = F)
  MST5_5 <- mstree(as.dist(data_dist5), ngmax=5)
  gSeg_MST5_5 <- gseg1(n, MST5_5, pval.perm = T, pval.appr = F)

  # Comparision
  mean_power <- mean_change(data_cp[[i]])

  results_p[i, ] <- c(
    # MST 5
    gSeg_MST5_0$pval.perm$ori$pval,
    gSeg_MST5_0$pval.perm$weighted$pval,
    gSeg_MST5_0$pval.perm$max.type$pval,
    gSeg_MST5_0$pval.perm$generalized$pval,
    # MST 5
    gSeg_MST5_1$pval.perm$ori$pval,
    gSeg_MST5_1$pval.perm$weighted$pval,
    gSeg_MST5_1$pval.perm$max.type$pval,
    gSeg_MST5_1$pval.perm$generalized$pval,
    # MST 5
    gSeg_MST5_2$pval.perm$ori$pval,
    gSeg_MST5_2$pval.perm$weighted$pval,
    gSeg_MST5_2$pval.perm$max.type$pval,
    gSeg_MST5_2$pval.perm$generalized$pval,
    # MST 5
    gSeg_MST5_3$pval.perm$ori$pval,
    gSeg_MST5_3$pval.perm$weighted$pval,
    gSeg_MST5_3$pval.perm$max.type$pval,
    gSeg_MST5_3$pval.perm$generalized$pval,
    # MST 5
    gSeg_MST5_4$pval.perm$ori$pval,
    gSeg_MST5_4$pval.perm$weighted$pval,
    gSeg_MST5_4$pval.perm$max.type$pval,
    gSeg_MST5_4$pval.perm$generalized$pval,
    # MST 5
    gSeg_MST5_5$pval.perm$ori$pval,
    gSeg_MST5_5$pval.perm$weighted$pval,
    gSeg_MST5_5$pval.perm$max.type$pval,
    gSeg_MST5_5$pval.perm$generalized$pval,

    # Mean
    ifelse(is.na(mean_power),1,0)
  )
}

## Null
colSums(results_p<=0.05)/nrow(results_p)
saveRDS(results_p,paste0('C:/Users/jerem/OneDrive/Documents/School/',
                         'Waterloo/Research/GraphCP/Data/Dist_pvals_200.rds'))

## Examined
#   MDT (1,3,5), ori, weighted, max.type, generalized
#   NNL (1,3,5), ori, weighted, max.type, generalized
#   MST (1,3,5), ori, weighted, max.type, generalized
## L2 Error
### N = 200, break at half
### 1.75 eig change
#   MDT -
#
#
#   NNL -
#
#
#   MST -
#
#
