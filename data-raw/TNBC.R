#' Title
#'
#' This function organizes the TNBC data.
#'
#' @return
#' @export
#'
#' @source <https://www.angelolab.com/mibi-data>
#'
#' @examples
organizeTNBC <- function(CDir){
  # Read in TNBC MIBI data
  tnbc <- read.csv(paste0(CDir, "downloaded/cellData.csv"))

  # Read the class data
  #   Classes:
  #       0: (mixed) tumor mixed with immune cells
  #       1: (compartmentalized) tumor cells and immune cells spatially separated
  #       2: (cold): no infiltrate
  tnbc_class <- read.csv(paste0(CDir, "/downloaded/patient_class.csv"), header=F)
  colnames(tnbc_class) <- c("Person","Class")

  # Read clinical data:
  tnbc_clinical = readxl::read_excel(paste0(CDir,"/downloaded/mmc2.xlsx"), skip = 1, n_max = 40)

  # Preprocss tiff files
  temp_file <- list.files(pattern = "*.tiff",path = paste0(CDir,'downloaded/'))
  temp_file <- gtools::mixedsort(temp_file) # reverse the name ordering
  all_file <- lapply(paste0(CDir,'downloaded/',temp_file), raster::raster)
  pixel_xy  <- lapply(all_file, function(x) raster::coordinates(x))
  pixel_val <- lapply(all_file, function(x) as.vector(x))
  N <- length(all_file)
  all_idx <- lapply(seq(1, N, 1), function(x) which(tnbc$SampleID == x))

  center_x <- center_y <- list()
  for (i in 1:N){
    temp_idx <- all_idx[[i]]
    temp_val <- pixel_val[[i]]
    temp_xy <- pixel_xy[[i]]
    temp_val <- lapply(tnbc$cellLabelInImage[temp_idx],
                       function(x) which(temp_val == x))
    center_x[[i]] <- lapply(temp_val, function(x) mean(temp_xy[x,1]))
    center_y[[i]] <- lapply(temp_val, function(x) mean(temp_xy[x,2]))
  }

  # Store the unlisted x and y data:
  write.csv(unlist(center_x), file = paste0(CDir,'generated/',"center_x.csv"))
  write.csv(unlist(center_y), file = paste0(CDir,'generated/',"center_y.csv"))

  # Store data frame of x-y data and marks:
  tnbc41 <- tnbc[tnbc$SampleID %in% 1:41,]
  tnbc41_xy <- cbind(tnbc41,unlist(center_x),unlist(center_y))
  colnames(tnbc41_xy)[58:59] <- c("cellx","celly")

  write.csv(tnbc41_xy, file = paste0(CDir,'generated/',"cellData_xy.csv"))

  # Read the data generated above for analysis on subsequent runs:
  tnbc41_xy_read = read.csv(paste0(CDir,'generated/',"cellData_xy.csv"))

  # Marking T/F factors based on the 0.5 cutoff as used in the Angelo paper:
  tnbc41_xy_read_fact = tnbc41_xy_read
  for(i in colnames(tnbc41_xy_read)[5:54]){
    tnbc41_xy_read_fact[,i] <- as.factor(tnbc41_xy_read[,i]>0.5)
  }
  colnames(tnbc41_xy_read_fact)[2] <- "Person"

  # Add class data:
  tnbc41_xy_read_fact_class = data.frame(tnbc41_xy_read_fact, Class = NA)
  for(i in unique(tnbc41_xy_read_fact_class$Person)){
    tnbc41_xy_read_fact_class[tnbc41_xy_read_fact_class$Person==i,]$Class <-
      tnbc_class[tnbc_class$Person==i,]$Class
  }

  write.csv(tnbc41_xy_read_fact_class,
            file = paste0(CDir,'generated/',"cellData_xy_Class.csv"))
}

generateTNBC <- function(CDir){
  tnbc_clinical <- readxl::read_excel(paste0(CDir,"/downloaded/mmc2.xlsx"),
                                      skip = 1, n_max = 40)

  tnbc41_xy_read_fact_class <- read.csv(
    paste0(CDir,'generated/',"cellData_xy_Class.csv"))[,-1]
  # Similarly to SPF paper, just look at classes 0 and 1.
  tnbc_spatial_01 = tnbc41_xy_read_fact_class[-(which(tnbc41_xy_read_fact_class$Class==2)),]

  persons_clinical01 = intersect(tnbc_clinical$InternalId,tnbc_spatial_01$Person)
  tnbc_clinical01 = tnbc_clinical[tnbc_clinical$InternalId %in% persons_clinical01,]
  tnbc_spatial_01 = tnbc_spatial_01[tnbc_spatial_01$Person %in% tnbc_clinical01$InternalId,]

  write.csv(tnbc_clinical01,
            file = paste0(CDir,'generated/',"tnbc_clinical01.csv"))

  # x-y data to front:
  tnbc_class01_final <- tnbc_spatial_01[,c(59,60,61,1:58)]


  write.csv(tnbc_class01_final,
            file = paste0(CDir,'generated/',"tnbc_class01_final.csv"))

  tnbc_class01_final
}

cleanTNBC <- function(TNBC){
  TNBC[,c("Class","Person",'cellLabelInImage',"cellx","celly",
          "C","Na","Si","P","Ca","Fe","dsDNA","Vimentin","SMA","Background",
          "B7H3","FoxP3","Lag3","CD4","CD16","CD56","OX40","PD1","CD31","PD.L1",
          "EGFR","Ki67","CD209","CD11c","CD138","CD163","CD68","CSF.1R","CD8",
          "CD3","IDO","Keratin17","CD63","CD45RO","CD20","p53","Beta.catenin",
          "HLA.DR","CD11b","CD45","H3K9ac","Pan.Keratin","H3K27me3","phospho.S6",
          "MPO","Keratin6","HLA_Class_1","Ta","Au","tumorYN")]
}

# Save CSVs
organizeTNBC(CDir=
               "~/School/Waterloo/Research/RPackages/funkycells/data-raw/TNBC/")
# Organize TNBC Cells
TNBC_basic <- generateTNBC(
          CDir="~/School/Waterloo/Research/RPackages/funkycells/data-raw/TNBC/")
TNBC <- cleanTNBC(TNBC_basic)

# Get PCA
agents_df <- expand.grid(c("FoxP3","Lag3","CD4","CD16", "CD56","PD1","PD.L1",
                           "Ki67","CD138","CD68","CD8","CD3","IDO","CD45RO",
                           "CD20","p53","MPO","tumorYN"),
                         c("FoxP3","Lag3","CD4","CD16", "CD56","PD1","PD.L1",
                           "Ki67","CD138","CD68","CD8","CD3","IDO","CD45RO",
                           "CD20","p53","MPO","tumorYN"))
dataPCA <- getPCAData(data = TNBC[,-3], unit='Person', agents_df=agents_df,
                      rCheckVals = seq(0,50,1))

# Setup PCA Data
tnbc_clinical01 <- read.csv(
  paste0("~/School/Waterloo/Research/RPackages/funkycells/data-raw/TNBC/",
         'generated/',"tnbc_clinical01.csv"))[-1]
data_pca_age  = merge(dataPCA,tnbc_clinical01[,c('InternalId','AGE_AT_DX')],
                      by.x='Person',by.y='InternalId')
colnames(data_pca_age )[dim(data_pca_age )[2]] = "Age"

# Save Data
usethis::use_data(TNBC, overwrite = TRUE)

## Delete me
data_pca_age

rfcv <- computeRandomForest_CVPC(data=data_pca_age, K=10,
                                 outcome='Class',unit='Person',repeatedId=NULL,
                                 metaNames=c('Age'),cellData= TNBC[,-3],
                                 syntheticKs=100, syntheticMetas=100,
                                 generalSyntheticK=T,alpha=0.05, silent=F)
