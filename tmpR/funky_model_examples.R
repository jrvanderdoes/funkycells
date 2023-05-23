agents_df <- data.frame(Var1 = c("Tumour", "Tumour", "Tumour",
                                 "Tumour", "Tumour", "Tumour",
                                 "Tumour", "Tumour", "Tumour",
                                 "Tumour", "Tumour"),
                        Var2 = c("Tumour", "CD3T", "CD4T",
                                 "CD8T", "B", "DC",
                                 "DCMono", "Macrophage", "MonoNeu",
                                 "NK", "Treg"))
dataPCA <- getPCAData(data = TNBC_pheno, unit = 'Person',
                      agents_df = agents_df, rCheckVals = seq(0,50,1))

dataPCA_Age  <- merge(dataPCA,TNBC_Meta)

set.seed(12345)
rfcv <- funkyModel(data=dataPCA_Age, K=10, outcome='Class',
                   unit='Person',metaNames=c('Age'),
                   synthetics=100, alpha=0.05, silent=FALSE,
                   subsetPlotSize = 25, nTrees=500)
rfcv$viPlot

# K function 1
tmp <- getKFunction(TNBC_pheno[TNBC_pheno$Class==0,-1],
                    agents = c('Tumour','Tumour'),unit = 'Person',
                    rCheckVals = seq(0,50,1))
tmp1 <- getKFunction(TNBC_pheno[TNBC_pheno$Class==1,-1],
                    agents = c('Tumour','Tumour'),unit = 'Person',
                    rCheckVals = seq(0,50,1))
data_plot <- rbind(data.frame('r'=seq(0,50,1),
                        'K'=rowMeans(tmp[,-1],na.rm = TRUE),
                        'Class'="0"),
                   data.frame('r'=seq(0,50,1),
                              'K'=rowMeans(tmp1[,-1],na.rm = TRUE),
                              'Class'='1'))

ggplot2::ggplot(data=data_plot,
                ggplot2::aes(x=r,y=K,group=Class,color=Class))+
  ggplot2::geom_line() +
  ggplot2::theme_bw()

# K function 2
tmp <- getKFunction(TNBC_pheno[TNBC_pheno$Class==0,-1],
                    agents = c('Tumour','MonoNeu'),unit = 'Person',
                    rCheckVals = seq(0,50,1))
tmp1 <- getKFunction(TNBC_pheno[TNBC_pheno$Class==1,-1],
                     agents = c('Tumour','MonoNeu'),unit = 'Person',
                     rCheckVals = seq(0,50,1))
data_plot <- rbind(data.frame('r'=seq(0,50,1),
                              'K'=rowMeans(tmp[,-1],na.rm = TRUE),
                              'Class'="0"),
                   data.frame('r'=seq(0,50,1),
                              'K'=rowMeans(tmp1[,-1],na.rm = TRUE),
                              'Class'='1'))

ggplot2::ggplot(data=data_plot,
                ggplot2::aes(x=r,y=K,group=Class,color=Class))+
  ggplot2::geom_line() +
  ggplot2::theme_bw()

# K function 3 (Fake)
tmp <- getKFunction(TNBC_pheno[TNBC_pheno$Class==0,-1],
                    agents = c('Tumour','Treg'),unit = 'Person',
                    rCheckVals = seq(0,50,1))
tmp1 <- getKFunction(TNBC_pheno[TNBC_pheno$Class==1,-1],
                     agents = c('Tumour','Treg'),unit = 'Person',
                     rCheckVals = seq(0,50,1))
data_plot <- rbind(data.frame('r'=seq(0,50,1),
                              'K'=rowMeans(tmp[,-1],na.rm = TRUE),
                              'Class'="0"),
                   data.frame('r'=seq(0,50,1),
                              'K'=rowMeans(tmp1[,-1],na.rm = TRUE),
                              'Class'='1'))

ggplot2::ggplot(data=data_plot,
                ggplot2::aes(x=r,y=K,group=Class,color=Class))+
  ggplot2::geom_line() +
  ggplot2::theme_bw()
################################################

phenos <- unique(TNBC_pheno$Phenotype)
dataPCA_noInv <- getPCAData(data = TNBC_pheno, unit='Person',
                                  agents_df=rbind(data.frame(t(combn(phenos,2))),
                                                  data.frame('X1'=phenos,'X2'=phenos)),
                                  rCheckVals = seq(0,50,1))

dataPCA_Age_NoInv  <- merge(dataPCA_noInv,TNBC_Meta)

set.seed(12345)
rfcv1 <- funkyModel(data=dataPCA_Age_NoInv, K=10, outcome='Class',
                   unit='Person',metaNames=c('Age'),
                   synthetics=100, alpha=0.05, silent=FALSE,
                   subsetPlotSize = 25, nTrees=500)
rfcv1$subset_viPlot

################################################

proteins <- c("FoxP3","Lag3","CD4","CD16","CD56","PD1","PD.L1",
              "Ki67","CD138","CD68","CD8","CD3","IDO","CD45RO",
              "CD20","p53","MPO","tumorYN")
dataPCA_int_noInv <- getPCAData(data = TNBC[,-3], unit='Person',
                            agents_df=rbind(data.frame(t(combn(proteins,2))),
                                            data.frame('X1'=proteins,'X2'=proteins)),
                            rCheckVals = seq(0,50,1))

dataPCA_int_Age_noInv  <- merge(dataPCA_int_noInv,TNBC_Meta)

set.seed(12345)
rfcv2 <- funkyModel(data=dataPCA_int_Age_noInv, K=10, outcome='Class',
                    unit='Person',metaNames=c('Age'),
                    synthetics=100, alpha=0.05, silent=FALSE,
                    subsetPlotSize = 25, nTrees=500)
rfcv2$subset_viPlot

################################################
phenos <- unique(TNBC_pheno$Phenotype)
pheno_interactions <- rbind(data.frame(t(combn(phenos,2))),
                            data.frame('X1'=phenos,'X2'=phenos))

phenos_subset <- c("Tumour","CD3T", "CD4T", "CD8T", "B", "DC",
                   "DCMono", "Macrophage", "MonoNeu", "NK", "Treg")
pheno_interactions_subset <- data.frame(
  Var1 = rep("Tumour",11),
  Var2 = c("Tumour", "CD3T", "CD4T", "CD8T", "B", "DC",
           "DCMono", "Macrophage", "MonoNeu", "NK", "Treg"))

dataPCA_pheno<- getPCAData(data = TNBC_pheno, unit='Person',
                           agents_df=pheno_interactions_subset,
                           rCheckVals = seq(0,50,1))

dataPCAAge_pheno  <- merge(dataPCA_pheno,TNBC_Meta)

set.seed(12345)
rfcv <- funkyModel(data=dataPCAAge_pheno, K=10,
                   outcome='Class',unit='Person',
                   metaNames=c('Age'), synthetics=500,
                   alpha=0.05, silent=FALSE,
                   subsetPlotSize = 25,
                   rGuessSims = 500, nTrees = 500,method = 'class')
