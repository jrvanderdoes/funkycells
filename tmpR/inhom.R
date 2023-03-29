# library(fda)
# library(survival)
# library(survminer)
# library(dplyr)
# library(tidyverse)
# library(readxl)
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
# library(gtools)
# library(raster)
# library(rainbow)
## Confirm
library(readxl)
library(gtools)
library(raster)
require(SPIAT)

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("SPIAT")
# devtools::install_github("TrigosTeam/SPIAT")


# The classes are associated with the following:
#0: (mixed) tumor mixed with immune cells
#1: (compartmentalized) tumor cells and immune cells spatially separated
#2: (cold): no infiltrate

#####
# Load source functions
# setwd("C:/Users/jerem/Downloads/NBC")

## Load Data
CDir =  "C:/Users/jerem/Downloads/TNBC_MyCode/"
tnbc = read.csv(paste0(CDir, "downloaded/cellData.csv"))

tnbc_class = read.csv(paste0(CDir, "downloaded/patient_class.csv"), header = F)#[,2]
colnames(tnbc_class) = c("Person","Class")

tnbc_clinical = read_excel(paste0(CDir,"downloaded/mmc2.xlsx"), skip = 1, n_max = 40)

preprocess <- function(savePath, loadpath){
  orig_wd <- getwd(); on.exit(setwd(orig_wd))
  setwd(loadpath)

  temp_file = list.files(pattern = "*.tiff")
  temp_file = mixedsort(temp_file) # reverse the name ordering
  all_file = lapply(temp_file, raster) #2048x2048 images
  pixel_xy  = lapply(all_file, function(x) coordinates(x))
  pixel_val = lapply(all_file, function(x) as.vector(x))
  N = length(all_file)
  all_idx = lapply(seq(1, N, 1), function(x) which(tnbc$SampleID == x))
  # extract locations & pixel values
  center_x = list()
  center_y = list()
  for (i in 1:N){
    print(i)
    temp_idx = all_idx[[i]]
    temp_val = pixel_val[[i]]
    temp_xy = pixel_xy[[i]]
    temp_val = lapply(tnbc$cellLabelInImage[temp_idx], function(x) which(temp_val == x))
    center_x[[i]] = lapply(temp_val, function(x) mean(temp_xy[x,1]))
    center_y[[i]] = lapply(temp_val, function(x) mean(temp_xy[x,2]))
  }

  # Store the unlisted x and y data:
  write.csv(unlist(center_x), file = paste0(savePath, "center_x.csv"))
  write.csv(unlist(center_y), file = paste0(savePath, "center_y.csv"))

  # Store data frame of x-y data and marks:
  tnbc41 = tnbc[tnbc$SampleID %in% 1:41,]
  tnbc41_xy = cbind(tnbc41,unlist(center_x),unlist(center_y))
  colnames(tnbc41_xy)[58:59] = c("cellx","celly")

  write.csv(tnbc41_xy, file = paste0(savePath, "cellData_xy.csv"))
}
# preprocess(paste0(CDir,'generated/'),paste0(CDir,'downloaded/'))
tnbc41_xy_read = read.csv(paste0(CDir,'generated/',"cellData_xy.csv"))
colnames(tnbc41_xy_read)[2] = "Person"

tnbc41_xy_read_class = data.frame(tnbc41_xy_read, Class = NA)
for(i in unique(tnbc41_xy_read_class$Person)){
  tnbc41_xy_read_class[tnbc41_xy_read_class$Person==i,]$Class <-
    tnbc_class[tnbc_class$Person==i,]$Class
}
write.csv(tnbc41_xy_read_class, file = paste0(CDir,"generated/cellData_xy_Class.csv"))

# Similarly to SPF paper, just look at classes 0 and 1.
tnbc_spatial_01 = tnbc41_xy_read_class[-(which(tnbc41_xy_read_class$Class==2)),]

persons_clinical01 = intersect(tnbc_clinical$InternalId,tnbc_spatial_01$Person)
tnbc_clinical01 = tnbc_clinical[tnbc_clinical$InternalId %in% persons_clinical01,]
tnbc_spatial_01 = tnbc_spatial_01[tnbc_spatial_01$Person %in% tnbc_clinical01$InternalId,]

#tnbc_class01_final = tnbc_spatial_01[,c(59,60,61,1:58)]
data = tnbc_spatial_01[,c(59,60,61,1:58)]


# Following lines are copied from the bottom README file.
# Details the entries of the Group and ImmuneGroup columns.


# 55. Group - One of six primary cell classifications (1- Unidentified, 2- Immune, 3- Endothelial, 4- Mesenchymal-like, 5- Tumor, 6- Keratin-positive tumor
# 57. ImmuneGroup - Breakdown of the immune cells into categories (1- Tregs, 2- CD4 T, 3- CD8 T, 4- CD3 T, 5- NK, 6- B, 7- Neutrophils, 8- Macrophages, 9- DC, 10- DC/Mono, 11- Mono/Neu, 12- Other immune)

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

data_plot <- data[which(data$Person==unique(data$Person)[11]),] # nums 2, 11, 17
data_plot$cellType <-
  ifelse(data_plot$Phenotype%in%c("Tumour", "CD4 T", "CD8 T", "B"),
         data_plot$Phenotype,'Other')
data_plot <- data_plot[,c('cellx','celly','cellType')]

ggplot2::ggplot() +
  ggplot2::geom_point(mapping=ggplot2::aes(x=data_plot[,1],y=data_plot[,2],
                                           col=data_plot[,3]),
                      size=2) +
  ggplot2::theme_bw() +
  ggplot2::xlim(c(min(data_plot[,1]),max(data_plot[,1]))) +
  ggplot2::ylim(c(min(data_plot[,2]),max(data_plot[,2]))) +
  ggplot2::guides(color='none') + #guide_legend(title="Phenotype")) +#
  ggplot2::theme(axis.title = ggplot2::element_blank(),
                 axis.text = ggplot2::element_blank(),
                 legend.text = ggplot2::element_text(size=24),
                 legend.title = ggplot2::element_text(size=28)) +
  ggplot2::scale_color_manual(values = c( "red", "blue", "darkcyan", "darkgreen", 'gray'),
                              breaks=c("Tumour", "CD4 T", "CD8 T", "B", "Other"))

# phenotypes <- data_plot$Phenotype
# coord_x <- data_plot$cellx
# coord_y <- data_plot$celly
# imagep2 <- format_image_to_spe(format = "general",
#                                phenotypes = phenotypes,
#                                coord_x = coord_x,coord_y = coord_y)
# my_colors <- c("red", "blue", "darkcyan", "darkgreen")
#
# # pdf(paste("0_2", ".pdf", sep=""))
# plot_cell_categories(spe_object = imagep2,
#                      categories_of_interest =
#                        c("Tumour", "CD4 T", "CD8 T", "B"),
#                      colour_vector = my_colors, feature_colname = "Phenotype")
# dev.off()


# Example of plotting a "homogeneous" tumour data.

datap0 =data0[which(data0$Person==unique(data0$Person)[5]),]
datap0$Phenotype = as.factor(datap0$Phenotype)

window = as.owin(c(c(30,2018),c(30,2018)))

ppp.data_person0 = as.ppp(datap0[,c('cellx','celly','Phenotype')],
                          W=window)

# Can produce K cross function for above data for cell types in Phenotype.

K0 = Kcross(ppp.data_person0,"CD4 T", "CD8 T", r=seq(0,200,length.out=513))
plot(K0)


# Analysis of inhomogeneous tumour.

datap1 =data1[which(data1$Person==unique(data1$Person)[3]),]

datap1$Phenotype = as.factor(datap1$Phenotype)
window = as.owin(c(c(30,2018),c(30,2018)))
ppp.data_person1 = as.ppp(datap1[,c('cellx','celly','Phenotype')],
                          W=window)

# Homogeneous K cross function

K = Kcross(ppp.data_person1,"CD4 T","CD8 T",r=seq(0,100,length.out=513))
plot(K)

# Following appears to produce a point pattern made up of cells that
# are not tumour cells.

ppp.other = ppp.data_person1[ppp.data_person1$marks != "Tumour"]

# Make a density plot of the non-tumour cells.
# See "density.ppp" for details.

dd2 = density(ppp.other, sigma=50)

# Following gives the proportion of the cell of interest, some non-tumour cells,
# w.r.t. the total non-tumour cells.

fact1 = sum(ppp.data_person1$marks == "CD4 T")/sum(ppp.data_person1$marks != "Tumour")
fact2 = sum(ppp.data_person1$marks == "CD8 T")/sum(ppp.data_person1$marks != "Tumour")

# We adjust the smoothed density by these factors.

dd1f = dd2 * fact1
dd2f = dd2 * fact2

# We plot one of the above (they will look the same, and just differ by the scale bar).
# We add the non-tumour points to see how smoothed density matches up with the points.

plot(dd2f)
points(ppp.other)

# Note that we have not smoothed each specific cell type, e.g. CD4 and CD8 T cells.
# The rationale behind using ALL the non-tumour cells for the smoothing part,
# and then adjusting by the proportions of cell type is the following:
# Visually from the plots of cell positions, it looks like the density of non-tumour
# cells is pretty uniform within the stroma. We want to compare the distribution
# of cell phenotypes of interest, e.g. CD4 and CD8 T cells, against the situation
# in which these cell types are randomly distributed amongst the stromal cells

# It might conceivably be better to define distinct boundaries of the tumour/stroma interface
# and then estimate a constant density of non-tummour cells within the stroma, and setting
# the tumour/stroma interface to be a boundary, with associted boundary correction for cells
# bordering this.

Kin = Kcross.inhom(ppp.data_person1,"CD4 T","CD8 T",dd1f,dd2f,r=seq(0,100,length.out=513))
plot(Kin)

# In addition, it might even be better to avoid usng the inhomogeneous K function in its entirety.
# One could compare the homogeneous K function against a simulated the homogeneous K function obtained
# by randomly selecting cells from the stromal cells and constructing K functions - the average would
# conceivably give the K functions associated with CSR for the given distribution of stromal cells.
# selecting non-tumour cells

