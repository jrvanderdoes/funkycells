# library(fda)
# library(survival)
# library(survminer)
# library(dplyr)
# library(tidyverse)
# library(janitor)
# library(skimr)
# library(data.table)
# library(pivottabler)
# library(ggpubr)
# library(zoo)
# library(lubridate)
# library(reReg)
# library(spatstat)
# library(mgcv)
# library(refund)
# library(ggplot2)
# library(gridExtra)
# library(rainbow)
## Confirm
library(readxl)
library(gtools)
library(raster)

#####

# Various lines below are commented out after initially reading the Tiff data, etc.
# Cell position data saved, and read each time I want to start analysis.

#####

# Load source functions

#setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/IHC_spatial/TNBC_shareCellData")
#source("~/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/IHC_spatial/SPF-main/ref_functions.R")

# read in TNBC MIBI data

CDir =  "C:/Users/jerem/Downloads/TNBC/"

# Read the class data:

tnbc_class = read.csv(paste0(CDir, "downloaded/patient_class.csv"), header = F)#[,2]
colnames(tnbc_class) = c("Person","Class")

# The classes are associated with the following:
#0: (mixed) tumor mixed with immune cells
#1: (compartmentalized) tumor cells and immune cells spatially separated
#2: (cold): no infiltrate

# Read clinical data:
# This file is found in the supplementary file of the main TNBC MIBI Angelo paper.

tnbc_clinical = read_excel(paste0(CDir,"downloaded/mmc2.xlsx"),
                           skip = 1, n_max = 40)


### Uncomment to load tif images. This takes a long time.

# PREPROCESS DATA
# read in multiple tiff files

# There's no subject ID 30 in the csv file.
# setwd("C:/Users/jerem/Downloads/TNBC/downloaded/")
# tnbc = read.csv(paste0(CDir,"downloaded/cellData.csv"))
# temp_file = list.files(pattern = "*.tiff")
# temp_file = mixedsort(temp_file) # reverse the name ordering
# all_file = lapply(temp_file, raster) #2048x2048 images
# pixel_xy  = lapply(all_file, function(x) coordinates(x))
# pixel_val = lapply(all_file, function(x) as.vector(x))
# N = length(all_file)
# all_idx = lapply(seq(1, N, 1), function(x) which(tnbc$SampleID == x))
# # extract locations & pixel values
# center_x = list()
# center_y = list()
# for (i in 1:N){
#   print(i)
#   temp_idx = all_idx[[i]]
#   temp_val = pixel_val[[i]]
#   temp_xy = pixel_xy[[i]]
#   temp_val = lapply(tnbc$cellLabelInImage[temp_idx], function(x) which(temp_val == x))
#   center_x[[i]] = lapply(temp_val, function(x) mean(temp_xy[x,1]))
#   center_y[[i]] = lapply(temp_val, function(x) mean(temp_xy[x,2]))
# }


# Store the unlisted x and y data:
# write.csv(unlist(center_x), file = paste0(CDir,"generated/center_x.csv"))
# write.csv(unlist(center_y), file = paste0(CDir,"generated/center_y.csv"))

# Need to double check the details, but there are differences in the patient numbers for
# MIBI data file, Tiff files, clinical data. Continue considering the "41" patient data set.

# Store data frame of x-y data and marks:

# tnbc41 = tnbc[tnbc$SampleID %in% 1:41,]
# tnbc41_xy = cbind(tnbc41,unlist(center_x),unlist(center_y))
# colnames(tnbc41_xy)[58:59] = c("cellx","celly")

# write.csv(tnbc41_xy, file = paste0(CDir,"generated/cellData_xy.csv"))


# Read the data generated above for analysis on subsequent runs:

tnbc41_xy_read = read.csv(paste0(CDir,"generated/cellData_xy.csv"))


# Continue by making marks T/F factors based on the 0.5 cutoff as used in the Angelo paper:

tnbc41_xy_read_fact = tnbc41_xy_read
for(i in colnames(tnbc41_xy_read)[5:54]){
  tnbc41_xy_read_fact[,i] = as.factor(tnbc41_xy_read[,i]>0.5)
}
colnames(tnbc41_xy_read_fact)[2] = "Person"

# Plot survival data as hist.
# hist((tnbc_clinical$`Survival_days_capped*`)/365,breaks = 20)


# Add class data:
tnbc41_xy_read_fact_class = data.frame(tnbc41_xy_read_fact,
                                       Class = NA)
for(i in unique(tnbc41_xy_read_fact_class$Person)){
  tnbc41_xy_read_fact_class[tnbc41_xy_read_fact_class$Person==i,]$Class = tnbc_class[tnbc_class$Person==i,]$Class
}

# write.csv(tnbc41_xy_read_fact_class, file =  paste0(CDir,"generated/cellData_xy_Class.csv"))


# Similarly to SPF paper, just look at classes 0 and 1.
tnbc_spatial_01 = tnbc41_xy_read_fact_class[-(which(tnbc41_xy_read_fact_class$Class==2)),]

persons_clinical01 = intersect(tnbc_clinical$InternalId,tnbc_spatial_01$Person)
tnbc_clinical01 = tnbc_clinical[tnbc_clinical$InternalId %in% persons_clinical01,]
tnbc_spatial_01 = tnbc_spatial_01[tnbc_spatial_01$Person %in% tnbc_clinical01$InternalId,]

# x-y data to front:
tnbc_class01_final = tnbc_spatial_01[,c(59,60,61,1:58)]


### Plot
data = tnbc_class01_final


# Following lines are copied from the bottom README file.
# Details the entries of the Group and ImmuneGroup columns.

# 55. Group - One of six primary cell classifications
#     (1- Unidentified, 2- Immune, 3- Endothelial, 4- Mesenchymal-like, 5- Tumor, 6- Keratin-positive tumor
# 57. ImmuneGroup - Breakdown of the immune cells into categories
#     (1- Tregs, 2- CD4 T, 3- CD8 T, 4- CD3 T, 5- NK, 6- B, 7- Neutrophils, 8- Macrophages, 9- DC,
#      10- DC/Mono, 11- Mono/Neu, 12- Other immune)

# Following just makes new column called Phenotype - following fills column with the cell phenotypes:

data$Phenotype = data$Group
data$Phenotype[which(data$Phenotype==1)] = "Other"
data$Phenotype[which(data$Phenotype==3)] = "Endothelial"
data$Phenotype[which(data$Phenotype==4)] = "Mesenchymal"
data$Phenotype[which(data$Phenotype==5)] = "Tumour"
data$Phenotype[which(data$Phenotype==6)] = "Tumour"

data_immune = data$immuneGroup
data_immune[which(data_immune==1)] = "Treg"
data_immune[which(data_immune==2)] = "CD4 T"
data_immune[which(data_immune==3)] = "CD8 T"
data_immune[which(data_immune==4)] = "CD3 T"
data_immune[which(data_immune==5)] = "NK"
data_immune[which(data_immune==6)] = "B"
data_immune[which(data_immune==7)] = "Neutrophil"
data_immune[which(data_immune==8)] = "Macrophage"
data_immune[which(data_immune==9)] = "DC"
data_immune[which(data_immune==10)] = "DC/Mono"
data_immune[which(data_immune==11)] = "Mono/Neu"
data_immune[which(data_immune==12)] = "Other immune"


data$Phenotype[which(data_immune!=0)] = data_immune[which(data_immune!=0)]


# Separate data based on Classes 0 and 1.

data0 = data[which(data$Class==0),]
data1 = data[which(data$Class==1),]


# Can use the following code to plot individual plots of cell types:

datap2 =data0[which(data0$Person==unique(data0$Person)[5]),]

# I have tested out and used the SPIAT package for plotting:

phenotypes <- datap2$Phenotype
coord_x <- datap2$cellx
coord_y <- datap2$celly
imagep2 <- format_image_to_spe(format = "general",
                               phenotypes = phenotypes,
                               coord_x = coord_x,coord_y = coord_y)
my_colors <- c("red", "blue", "darkcyan", "darkgreen")

# pdf(paste("0_2", ".pdf", sep=""))
plot_cell_categories(spe_object = imagep2,
                     categories_of_interest =
                       c("Tumour", "CD4 T", "CD8 T", "B"),
                     colour_vector = my_colors, feature_colname = "Phenotype")


kfun <- getKFunction(data = datap2[,!(colnames(datap2)%in%c('Class'))],
             agents = c('CD4','PD.L1'),
             unit = 'Person',rCheckVals = seq(0,50,1))
ggplot(data=kfun) +
  geom_line(aes(x=r,y=K1),size=1.5) +
  geom_line(aes(x=r,y=pi*r^2),size=1.5, linetype='dashed', color='red') +
  theme_bw() +
  theme(axis.text.y = element_blank(),
        axis.text.x = element_text(size=18),
        axis.title = element_blank())






## source various functions
# source("~/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/IHC_spatial/PCAdata_cross.R")
# source("~/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/IHC_spatial/PCAdata_multi_K.R")
# source("~/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/IHC_spatial/Kfunctions_Kcross.R")
# source("~/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/IHC_spatial/Kfunctions_Kmulti.R")
# source("~/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/IHC_spatial/GeneralFunctions.R")
# source("~/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/IHC_spatial/Randomforests_multi.R")
# source("~/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/IHC_spatial/myRF_CV_multi.R")


# Following allows you to check the min and max x-y values that cells have,
# and also the number of cells per mark per image (need to change number of cols
# depending on number of marks you are looking at)

# rRange = matrix(0,nrow = length(unique(tnbc_class01_final$Person)),ncol = 4)
# cellsm = matrix(0,nrow = length(unique(tnbc_class01_final$Person)),ncol = 8)
#
# for(i in 1:length(unique(tnbc_class01_final$Person))){
#   data_person = tnbc_class01_final[tnbc_class01_final$Person==unique(tnbc_class01_final$Person)[i],]
#   rRange[i,] = c(min(data_person$cellx),max(data_person$cellx),min(data_person$celly),max(data_person$celly))
#   cellsm[i,] = colSums(data_person[,c('tumorYN','CD4','CD8','CD31','CD68','PD.L1','PD1','FoxP3')]==TRUE)
# }
# rmeans = colMeans(rRange)
# rvol = (rmeans[2] - rmeans[1])*(rmeans[4] - rmeans[3])
# cmeans = colMeans(cellsm)


# From the above the window size is quite large. Setting the maximum r value to the default
# of (minimum side length of the window)/4 is probably too large. Ranges of interactions are probably
# smaller than this. Arbitrarily continue by setting rmax = 50. Need to 1) discuss with ML regarding
# ranges of r values we should consider, or 2) develop some approach to optimise r/check multiple r values.

# I have adapted the getPCAdata function to account for the marked nature of cells, rather than dealing
# with distinct cell types.

# In addition, this function also spits out the K functions and the r values so that you can plot these.

# Analysis for two groups of Class response variable, non-exhaustive or optimised collection of mark interactions:


cellTypes <- c("FoxP3", "Lag3", "CD4", "CD16",
               "CD56", "PD1", "PD.L1", "Ki67",
               "CD138", "CD68", "CD8", "CD3",
               "IDO", "CD45RO", "CD20", "p53",
               "MPO", "tumorYN")
cellDat <- tnbc_class01_final[!(colnames(tnbc_class01_final) %in%
                                  c('X','cellSize','cellLabelInImage','tumorCluster',
                                    'Group','immuneCluster','immuneGroup'))]
colnames(cellDat)[3] <- 'Stage'
for(i in 5:ncol(cellDat)){
  cellDat[,i] <- as.logical(cellDat[,i])
}

agents_df_tmp <- as.data.frame(expand.grid(cellTypes,cellTypes,
                                           stringsAsFactors=F))
data_pca <- getPCAData(cellDat,
           outcome = "Stage", unit='Person', repeatedUniqueId=NULL,
           rCheckVals=seq(0,50,0.05), nPCs=3,
           agents_df = agents_df_tmp,
           xRange=NULL, yRange=NULL,
           edgeCorrection="isotropic", nbasis=21,
           silent=F)


# Add age meta variable (though on subsequent analysis this does not seem to be a significant
# variable)

data_pca_age  = merge(data_pca,tnbc_clinical01[c('InternalId','AGE_AT_DX')],
                      by.x='Person',by.y='InternalId')
colnames(data_pca_age)[dim(data_pca_age)[2]] = "Age"

# I am continuing by just changing the names of response variables to "Stage", whatever
# their actual name is.

tnbc_class01_final_stage = tnbc_class01_final
colnames(tnbc_class01_final_stage)[which(colnames(tnbc_class01_final_stage)=="Class")] = "Stage"
colnames(data_pca_age)[which(colnames(data_pca_age)=="Class")] = "Stage"

# Following gives an idea about what might be important variables on subsequent analysis.

cellDat <- cellDat[c('Stage','Person',cellTypes)]
set.seed(12345)
RFm = computeRandomForest_CVPC(data=data_pca_age, K=10,
                               outcome='Stage', unit='Person',
                               repeatedId=NULL, metaNames="Age",
                               cellData=cellDat,
                               synthetics=100,
                               generalSyntheticK=T,
                               curvedSigSims=100, alpha=0.05,
                               silent=F,
                               rGuessSims=500, alignmentMethod=c('Add','Mult'),
                               subsetPlotSize=25, nTrees=500)
RFm$multIntVI$subset_viPlot

# # Arbitrarily continue by analysing interactions with varImp >=0.25.
# #
# # newVars = RFm[[2]][((RFm[[2]]$varImp)/max(RFm[[2]]$varImp))>=0.25,]
# # newVars = newVars[order(newVars$varImp),]
# #
# # Marks_df = data.frame('Mark1'=NA,
# #                       'Mark2'=NA)
# # for(i in 1:(dim(newVars)[[1]])){
# #   Marks_df[i,] = str_split(newVars$var[i],"_",simplify=T)
# # }
# #
# # # > Marks_df
# # # Mark1   Mark2
# # # 1      CD3 tumorYN
# # # 2  tumorYN  CD45RO
# # # 3      CD4   PD.L1
# # # 4    PD.L1     CD4
# # # 5   CD45RO tumorYN
# # # 6  tumorYN     CD3
# # # 7    PD.L1  CD45RO
# # # 8   CD45RO   PD.L1
# # # 9  tumorYN     CD4
# # # 10   PD.L1     CD3
# # # 11 tumorYN   PD.L1
# # # 12     CD3   PD.L1
# # # 13     CD4 tumorYN
# # # 14   PD.L1 tumorYN
# # # 15     PD1   PD.L1
# # # 16   PD.L1     PD1
# # # 17    CD16   PD.L1
# # # 18   PD.L1    CD16
# # # 19    CD68   PD.L1
# # # 20   PD.L1    CD68
# # # 21 tumorYN tumorYN
# #
# # data_pca2 <- getKPCAData_multi(tnbc_class01_final,
# #                                outcome = "Class",
# #                                nPCs=3,
# #                                Marks_df = Marks_df,
# #                                rmax = 50,
# #                                edgeCorrection="isotropic",
# #                                nbasis=21)
# #
# #
# # data_pca_final = data_pca2[[1]]
# # colnames(data_pca_final)[which(colnames(data_pca_final)=="Class")] = "Stage"
# #
# # RFm2 = myRF_CV_multi_singleimage(data_pca_final,
# #                                  plotVarImp=T,
# #                                  metaNames=NULL,
# #                                  MarkNames=NULL,
# #                                  cellData=tnbc_class01_final_stage,
# #                                  dropCols=c('Person'),
# #                                  FakePairs=100,
# #                                  FakeMetas=100,
# #                                  generalFakeCell=T)
#
#
#
#
#
#
#
# # Following is a very basic approach to plotting the K functions for various interactions.
# # "InterMark" specifies the interaction of interest, and is just the number in which the pair
# # of marks occurs in the Marks_df variable in getKPCAData_multi.
#
# Kdata = data_pca[[2]]
# InterMark = 2
#
# plot(Kdata[[InterMark]][[2]][,1],
#      Kdata[[InterMark]][[1]][,1],
#      type="l",
#      col=1+data_pca_age$Stage[1],
#      xlim=range(Kdata[[InterMark]][[2]],na.rm = T),
#      ylim=range(Kdata[[InterMark]][[1]],na.rm = T))
# for(i in 2:dim(Kdata[[InterMark]][[2]])[2]){
#   lines(Kdata[[InterMark]][[2]][,i],
#         Kdata[[InterMark]][[1]][,i],
#         col=1+data_pca_age$Stage[i])
# }
#
#
#
#
#
#
# tnbc_spatial_0 = tnbc41_xy_read_fact_class[(which(tnbc41_xy_read_fact_class$Class==0)),]
#
# persons_clinical0 = intersect(tnbc_clinical$InternalId,tnbc_spatial_0$Person)
# tnbc_clinical0 = tnbc_clinical[tnbc_clinical$InternalId %in% persons_clinical0,]
# tnbc_spatial_0 = tnbc_spatial_0[tnbc_spatial_0$Person %in% tnbc_clinical0$InternalId,]
#
# # x-y data to front:
# tnbc_class0_final = tnbc_spatial_0[,c(59,60,61,1:58)]
#
#
#
# Surv_group = as.integer(tnbc_clinical0$`Survival_days_capped*`>1825)
#
# tnbc_surv0 = data.frame(tnbc_class0_final,
#                         Surv_group = NA)
# ptsurvclass = data.frame(Person = NA,Class = NA,Surv_group=NA)
# for(i in 1:length(unique(tnbc_surv0$Person))){
#   tnbc_surv0[tnbc_surv0$Person==(unique(tnbc_surv0$Person)[i]),]$Surv_group = Surv_group[i]
#   cls = unique(tnbc_surv0[tnbc_surv0$Person==(unique(tnbc_surv0$Person)[i]),]$Class)
#   ptsurvclass[i,] = c(unique(tnbc_surv0$Person)[i],cls,Surv_group[i])
# }
#
# tnbc_surv0_final = tnbc_surv0
#
#
# # Define additional meta variables based on variables that might influence outcomes;
# # These include p53 (~marker of higher rate of mutations),, Ki67 (marker of mitotic activity),
# # size of tumour cells (probably not important. Size of the nucleus, rather than the entire
# # cell probably more useful - as typically would mean higher grade, but worth considering)
# p53prop = c()
# Ki67prop = c()
# tumsize = c()
# for(i in 1:length(unique(tnbc_surv0_final$Person))){
#   personi = tnbc_surv0_final[tnbc_surv0_final$Person == unique(tnbc_surv0_final$Person)[i],]
#   p53prop[i] = mean(personi[personi$tumorYN == TRUE,"p53"]==TRUE)
#   Ki67prop[i] = mean(personi[personi$tumorYN == TRUE,"Ki67"]==TRUE)
#   tumsize[i] = mean(personi[personi$tumorYN == TRUE,"cellSize"])
# }
#
# data_pca0 <- getKPCAData_multi(tnbc_surv0_final,
#                                outcome = "Surv_group",
#                                nPCs=3,
#                                MarksTypes = c("FoxP3",
#                                               "Lag3",
#                                               "CD4",
#                                               "CD16",
#                                               "CD56",
#                                               "PD1",
#                                               "PD.L1",
#                                               "Ki67",
#                                               "CD138",
#                                               "CD68",
#                                               "CD8",
#                                               "CD3",
#                                               "IDO",
#                                               "CD45RO",
#                                               "CD20",
#                                               "p53",
#                                               "MPO",
#                                               "tumorYN"),
#                                rmax = 50,
#                                edgeCorrection="isotropic",
#                                nbasis=21)
#
#
# # Add age meta variable (though on subsequent analysis this does not seem to be a significant
# # variable)
#
# data_pca_age  = cbind(data_pca0[[1]],tnbc_clinical0$AGE_AT_DX)
# colnames(data_pca_age )[dim(data_pca_age )[2]] = "Age"
#
# # I am continuing by just changing the names of response variables to "Stage", whatever
# # their actual name is.
#
# tnbc_surv0_final_stage = tnbc_surv0_final
# colnames(tnbc_surv0_final_stage)[which(colnames(tnbc_surv0_final_stage)=="Surv_group")] = "Stage"
# colnames(data_pca_age)[which(colnames(data_pca_age)=="Surv_group")] = "Stage"
#
# # Following gives an idea about what might be important variables on subsequent analysis.
#
# data_pca_add  = cbind(data_pca_age,p53prop,Ki67prop,tumsize)
# colnames(data_pca_add)[which(colnames(data_pca_add)=="Surv_group")] = "Stage"
#
# RFm = myRandomForest_multi(data=data_pca_add, varImpPlot = T,
#                            metaNames=c("p53prop","Ki67prop"))
#
#
# # Arbitrarily continue by analysing interactions with varImp >=0.25.
#
# newVars = RFm[[2]][((RFm[[2]]$varImp)/max(RFm[[2]]$varImp))>=0.20,]
# newVars = newVars[order(newVars$varImp),]
#
# # Probably no meta variables in it but check:
# print(newVars)
#
# Marks_df = data.frame('Mark1'=NA,
#                       'Mark2'=NA)
# for(i in 1:(dim(newVars)[[1]])){
#   Marks_df[i,] = str_split(newVars$var[i],"_",simplify=T)
# }
#
# # > Marks_df
# # Mark1  Mark2
# # 1   CD45RO   Ki67
# # 2   CD45RO    CD8
# # 3      CD8 CD45RO
# # 4      CD8   Ki67
# # 5     Ki67    CD8
# # 6  tumorYN    CD8
# # 7      CD3   Ki67
# # 8     Ki67    CD3
# # 9     Ki67    CD4
# # 10     CD4   Ki67
#
# data_pca <- getKPCAData_multi(tnbc_surv0_final,
#                               outcome = "Surv_group",
#                               nPCs=3,
#                               Marks_df = Marks_df,
#                               rmax = 50,
#                               edgeCorrection="isotropic",
#                               nbasis=21)
#
#
# data_pca_final = data_pca[[1]]
# colnames(data_pca_final)[which(colnames(data_pca_final)=="Surv_group")] = "Stage"
#
# RFm2 = myRF_CV_multi_singleimage(data_pca_final,
#                                  plotVarImp=T,
#                                  metaNames=NULL,
#                                  MarkNames=NULL,
#                                  cellData=tnbc_surv0_final_stage,
#                                  dropCols=c('Person'),
#                                  FakePairs=100,
#                                  FakeMetas=100,
#                                  generalFakeCell=T)
#


data$Phenotype <- gsub('[ &//]','', data$Phenotype)
data_pca <- getPCAData(data[,c('Class','Person','cellx','celly','Phenotype')],
                       outcome = "Class", unit='Person', repeatedUniqueId=NULL,
                       rCheckVals=seq(0,50,0.05), nPCs=2,
                       agents_df = expand.grid(unique(data$Phenotype),unique(data$Phenotype)),
                       xRange=NULL, yRange=NULL,
                       edgeCorrection="isotropic", nbasis=21,
                       silent=F)

data_pca <- merge(data_pca,tnbc_clinical01[c('InternalId','AGE_AT_DX')],
                  by.x='Person',by.y='InternalId')
colnames(data_pca) <- c(colnames(data_pca)[-length(colnames(data_pca))],'Age')


data$cellType <- data$Phenotype
rf <- funkyRandomForest(data = data_pca,
                        K = 10,outcome = 'Class',unit = 'Person',
                        metaNames = c('Age'),synthetics = 100)
rf$subset_viPlot



