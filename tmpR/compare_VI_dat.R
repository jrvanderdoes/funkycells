# dataPCA <- getPCAData(data = TNBC[,-3],
#                       unit='Person', nPCs = 2,
#                       agents_df=expand.grid(proteins,proteins),
#                       rCheckVals = seq(0,50,1))
# 0.98

proteins <- c("FoxP3","Lag3","CD4","CD16","CD56","PD1",
              "PD.L1","Ki67","CD138","CD68","CD8","CD3",
              "IDO","CD45RO","CD20","p53","MPO","tumorYN")


set.seed(12345)
dat0 <- simulatePP(cellVarData=
                     data.frame('stage'=c(0,1),
                                'c1'=c(0,0),
                                'c2'=c(1/25,1/100), 'c3'=c(1/50,1/10),
                                'c4'=c(0,0),
                                'c5'=c(0,0),'c6'=c(0,0),
                                'c7'=c(0,0),'c8'=c(1/100,1/100),
                                'c9'=c(1/20,1/20), 'c10'=c(1/250,1/250),
                                'c11'=c(1/100,1/100),'c12'=c(1/80,1/80),
                                'c13'=c(0,0),'c14'=c(0,0),
                                'c15'=c(0,0),'c16'=c(0,0),
                                'c17'=c(0,0),'c18'=c(0,0)),
                   cellKappaData=data.frame(
                     'cell'=paste0('c',1:18),
                     'clusterCell'=c(NA,'c1','c1', rep(NA,15)),
                     'kappa'=c(rbinom(1,100,0.5),
                               rbinom(2,70,0.5),
                               rbinom(1,60,0.5),
                               rbinom(2,300,0.5),
                               rbinom(2,120,0.5),
                               rbinom(4,600,0.5),
                               rbinom(2,120,0.5),
                               rbinom(2,50,0.5),
                               rbinom(2,20,0.5))),
                   peoplePerStage=17,
                   imagesPerPerson=1,
                   silent=F)
dat1 <- simulatePP(cellVarData=
                     data.frame('stage'=c(0,1),
                                'c1'=c(0,0),
                                'c2'=c(1/25,1/25), 'c3'=c(1/50,1/50),
                                'c4'=c(0,0),
                                'c5'=c(0,0),'c6'=c(0,0),
                                'c7'=c(0,0),'c8'=c(1/100,1/100),
                                'c9'=c(1/20,1/20), 'c10'=c(1/250,1/250),
                                'c11'=c(1/100,1/100),'c12'=c(1/80,1/80),
                                'c13'=c(0,0),'c14'=c(0,0),
                                'c15'=c(0,0),'c16'=c(0,0),
                                'c17'=c(0,0),'c18'=c(0,0)),
                   cellKappaData=data.frame(
                     'cell'=paste0('c',1:18),
                     'clusterCell'=c(NA,'c1','c1', rep(NA,15)),
                     'kappa'=c(rbinom(1,100,0.5),
                               rbinom(2,70,0.5),
                               rbinom(1,60,0.5),
                               rbinom(2,300,0.5),
                               rbinom(2,120,0.5),
                               rbinom(4,600,0.5),
                               rbinom(2,120,0.5),
                               rbinom(2,50,0.5),
                               rbinom(2,20,0.5))),
                   peoplePerStage=1,
                   imagesPerPerson=1,
                   silent=F)
dat1$Person <- ifelse(dat1$Person=='p1','p35',
                      ifelse(dat1$Person=='p2','p36', NA))
dat1$Image <- ifelse(dat1$Image=='1','35',
                     ifelse(dat1$Image=='2','36',NA))
dat <- rbind(dat0,dat1)
set.seed(12345)
pcaData <- getPCAData(data = dat,repeatedUniqueId='Image',
                      xRange = c(0,1),  yRange = c(0,1), nPCs = 1,
                      silent=F)
set.seed(12345)
pcaMeta <- simulateMeta(pcaData,
                        metaInfo = data.frame(
                          'var'=c('gender','age'),
                          'rdist'=c('rbinom','rnorm'),
                          'Stage_0'=c('0.5','25'),
                          'Stage_1'=c('0.5','27.5')))
####
data=pcaMeta
K=10
outcome=colnames(data)[1]
unit=colnames(data)[2]
repeatedId='Image'
metaNames=c('age','gender')
cellData=dat
synthetics=100
alpha=0.05
silent=F
rGuessSims=500
subsetPlotSize=25
nTrees=500


#
# tmp1 <- pcaMeta[,!(colnames(pcaMeta) %in% c('Person'))]
# tmp = party::ctree(Stage ~ ., data=tmp1)
#
# tmp = party::cforest(Stage ~ ., data=tmp1)
# tmpval = party::varimp(tmp)
