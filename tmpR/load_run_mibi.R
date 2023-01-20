library(data.table)
data1 <- fread("C:/Users/j53vande/Downloads/1/NSCLCcohort_final_cleaned.csv")

data_simple <- data1[,-c(6:2285)]
colnames(data_simple) <- c('Image','Class','X','Y','InTumour','CancerType',
                           'Survival','SmokingStatus','Patient')
outcome <- 'CancerType'
data_simple <- as.data.frame(data_simple)[colnames(data_simple)!='Survival']
unit <- 'Patient'
repeatedId <- 'Image'

data_forK <- as.data.frame(data_simple)[,c(outcome,unit,repeatedId,'X','Y','Class')]
classes <- unique(data_simple$Class)
classes <- classes[classes!='Other']
classes_df <- expand.grid(classes,classes)


pcaData <- getPCAData(data=data_forK,
                      outcome = outcome,
                      unit = unit,
                      repeatedUniqueId = repeatedId,
                      agents_df = classes_df)
pcaMeta <- merge(x = pcaData$pcaData,
                 y = unique(data_simple[,c('Patient','SmokingStatus')]),
                 by='Patient')
tmp <- as.data.frame(data_simple)
colnames(tmp)[2] <- 'cellType'
rfcv <- computeRandomForest_CVPC(data = pcaMeta, outcome = 'CancerType',
                                 unit = 'Patient',
                                 repeatedId = repeatedId,
                                 metaNames=c(pcaData$metaNames,'SmokingStatus'),
                                 cellData = tmp,curvedSigSims=100,
                                 subsetPlotSize=25)

rfcv$multIntVI$viPlot


############
dropTypes <- c('Stromal cells', 'Epithelial cells')
classes_df_sub <- expand.grid(classes[!(classes%in%dropTypes)],
                              classes[!(classes%in%dropTypes)])


pcaData1 <- getPCAData(data=data_forK,outcome = outcome,unit = unit,
                       repeatedUniqueId = repeatedId,
                       agents_df = classes_df_sub)
pcaMeta1 <- merge(x = pcaData1,
                  y = unique(data_simple[,c('Patient','SmokingStatus')]),
                  by='Patient')
tmp1 <- as.data.frame(data_simple)
colnames(tmp1)[2] <- 'cellType'
rfcv1 <- computeRandomForest_CVPC(data = pcaMeta1, outcome = 'CancerType',
                                  unit = 'Patient',
                                  repeatedId = repeatedId,
                                  metaNames=c('SmokingStatus'),
                                  cellData = tmp1,curvedSigSims=100,
                                  subsetPlotSize=25)

rfcv1$multIntVI$viPlot
