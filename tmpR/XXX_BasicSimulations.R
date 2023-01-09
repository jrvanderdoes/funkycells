
set.seed(123456)
dat <- simulatePP(cellVarData=
                    data.frame('stage'=c(0,1),
                               'A'=c(0,0),
                               'B'=c(1/50,1/50)),
                  cellKappaData=data.frame(
                    'cell'=c('A','B'),
                    'clusterCell'=c(NA,'A'),
                    'kappa'=c(15,5)),
                  peoplePerStage=100,
                  imagesPerPerson=1,
                  silent=F )
pcaData <- getPCAData(dat,repeatedUniqueId='Image',
                      xRange = c(0,1),  yRange = c(0,1), silent=F)
pcaMeta <- simulateMeta(pcaData,
                        metaInfo = data.frame(
                          'var'=c('randUnif','randBin','corrNorm'),
                          'rdist'=c('runif','rbinom','rnorm'),
                          'Stage_0'=c('0.5','0.5','1'),
                          'Stage_1'=c('0.5','0.5','2')))

rfcv <- computeRandomForest_CVPC(data=pcaMeta,repeatedId='Image',
                                 metaNames=c('randUnif','randBin','corrNorm'),
                                 cellData=dat,
                                 curvedSigSims=100)


################
set.seed(123456)
dat1 <- simulatePP(cellVarData=
                    data.frame('stage'=c(0,1),
                               'A'=c(0,0),
                               'B'=c(1/50,1/250)),
                  cellKappaData=data.frame(
                    'cell'=c('A','B'),
                    'clusterCell'=c(NA,'A'),
                    'kappa'=c(15,5)),
                  peoplePerStage=100,
                  imagesPerPerson=1,
                  silent=F )
pcaData1 <- getPCAData(dat1,repeatedUniqueId='Image',
                      xRange = c(0,1),  yRange = c(0,1), silent=F)
pcaMeta1 <- simulateMeta(pcaData1,
                        metaInfo = data.frame(
                          'var'=c('randUnif','randBin','corrNorm'),
                          'rdist'=c('runif','rbinom','rnorm'),
                          'Stage_0'=c('0.5','0.5','1'),
                          'Stage_1'=c('0.5','0.5','2')))

rfcv1 <- computeRandomForest_CVPC(data=pcaMeta1,repeatedId='Image',
                                 metaNames=c('randUnif','randBin','corrNorm'),
                                 cellData=dat1,curvedSigSims=100)
