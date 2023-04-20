#' Organize TNBC
#'
#' This function organizes the TNBC data.
#'
#' @param CDir The directory to save/read files
#'
#' @source <https://www.angelolab.com/mibi-data>
organizeTNBCs <- function(CDir){
  # Read in TNBC MIBI data
  tnbc <- read.csv(paste0(CDir, "downloaded/cellData.csv"))

  # Read the class data
  #   Classes:
  #       0: (mixed) tumor mixed with immune cells
  #       1: (compartmentalized) tumor cells and immune cells spatially separated
  #       2: (cold): no infiltrate
  tnbc_class <- read.csv(paste0(CDir, "/downloaded/patient_class.csv"),
                         header=FALSE)
  colnames(tnbc_class) <- c("Person","Class")

  # Read clinical data:
  tnbc_clinical = readxl::read_excel(paste0(CDir,"/downloaded/mmc2.xlsx"),
                                     skip = 1, n_max = 40)

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
  #   Eventually could remove
  tnbc41_xy = read.csv(paste0(CDir,'generated/',"cellData_xy.csv"))

  # Marking T/F factors based on the 0.5 cutoff as used in the Angelo paper:
  tnbc41_xy_fact = tnbc41_xy
  for(i in colnames(tnbc41_xy)[5:54]){
    tnbc41_xy_fact[,i] <- as.factor(tnbc41_xy[,i]>0.5)
  }
  colnames(tnbc41_xy_fact)[2] <- "Person"
  colnames(tnbc41_xy)[2] <- "Person"

  # Add class data:
  tnbc41_xy_class <- data.frame(tnbc41_xy, Class = NA)
  for(i in unique(tnbc41_xy_class$Person)){
    tnbc41_xy_class[tnbc41_xy_class$Person==i,'Class'] <-
      tnbc_class[tnbc_class$Person==i,'Class']
  }

  tnbc41_xy_fact_class <- data.frame(tnbc41_xy_fact, Class = NA)
  for(i in unique(tnbc41_xy_fact_class$Person)){
    tnbc41_xy_fact_class[tnbc41_xy_fact_class$Person==i,'Class'] <-
      tnbc_class[tnbc_class$Person==i,'Class']
  }

  write.csv(tnbc41_xy_class,
            file = paste0(CDir,'generated/',"cellData_xy_Class_int.csv"))
  write.csv(tnbc41_xy_fact_class,
            file = paste0(CDir,'generated/',"cellData_xy_Class.csv"))
}


#' Generate TNBC Data
#'
#' This function generates the TNBC data.
#'
#' @param CDir The directory to save/read files
#'
#' @source <https://www.angelolab.com/mibi-data>
generateTNBCs <- function(CDir){
  tnbc_clinical <- readxl::read_excel(paste0(CDir,"/downloaded/mmc2.xlsx"),
                                      skip = 1, n_max = 40)

  tnbc41_xy_read_fact_class <- read.csv(
    paste0(CDir,'generated/',"cellData_xy_Class.csv"))[,-1]
  tnbc41_xy_class <- read.csv(
    paste0(CDir,'generated/',"cellData_xy_Class_int.csv"))[,-1]

  # Similarly to SPF paper, just look at classes 0 and 1.
  tnbc_spatial_01 = tnbc41_xy_read_fact_class[-(which(tnbc41_xy_read_fact_class$Class==2)),]
  tnbc_spatial_01_int = tnbc41_xy_class[-(which(tnbc41_xy_class$Class==2)),]

  persons_clinical01 = intersect(tnbc_clinical$InternalId,tnbc_spatial_01$Person)
  tnbc_clinical01 = tnbc_clinical[tnbc_clinical$InternalId %in% persons_clinical01,]
  tnbc_spatial_01 = tnbc_spatial_01[tnbc_spatial_01$Person %in% tnbc_clinical01$InternalId,]

  persons_clinical01_int = intersect(tnbc_clinical$InternalId,tnbc_spatial_01_int$Person)
  tnbc_clinical01_int = tnbc_clinical[tnbc_clinical$InternalId %in% persons_clinical01_int,]
  tnbc_spatial_01_int = tnbc_spatial_01_int[tnbc_spatial_01_int$Person %in% tnbc_clinical01_int$InternalId,]

  write.csv(tnbc_clinical01,
            file = paste0(CDir,'generated/',"tnbc_clinical01.csv"))
  write.csv(tnbc_clinical01_int,
            file = paste0(CDir,'generated/',"tnbc_clinical01_int.csv"))

  # x-y data to front:
  tnbc_class01_final <- tnbc_spatial_01[,c(59,60,61,1:58)]
  tnbc_class01_final_int <- tnbc_spatial_01_int[,c(59,60,61,1:58)]


  write.csv(tnbc_class01_final,
            file = paste0(CDir,'generated/',"tnbc_class01_final.csv"))
  write.csv(tnbc_class01_final_int,
            file = paste0(CDir,'generated/',"tnbc_class01_final_int.csv"))

  list('fact'=tnbc_class01_final,
       'int'=tnbc_class01_final_int)
}


#' Generate Phenotypes for TNBC
#'
#' This function add phenotypes to the TNBC data.
#'
#' @param CDir The directory to save/read files
#'
#' @source <https://www.angelolab.com/mibi-data>
generatePheno <- function(CDir){
  data <- read.csv(paste0(CDir,'generated/',"tnbc_class01_final.csv"))

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

  data$Phenotype <- gsub('[ &//]','', data$Phenotype)

  data[,c('Class','Person','cellx','celly','Phenotype')]
}


#' Organize TNBC
#'
#' This function cleans the TNBC data.
#'
#' @param CDir The directory to save/read files
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
organizeTNBCs(CDir=
               "~/School/Waterloo/Research/RPackages/funkycells/data-raw/TNBC/")
# Organize TNBC Cells
TNBC_basic <- generateTNBCs(
          CDir="~/School/Waterloo/Research/RPackages/funkycells/data-raw/TNBC/")
TNBC_pheno <- generatePheno(
          CDir="~/School/Waterloo/Research/RPackages/funkycells/data-raw/TNBC/")
# TNBC_int <- cleanTNBC(TNBC_basic$int)
TNBC <- cleanTNBC(TNBC_basic$fact)

# Get Meta-Info
tnbc_clinical01 <- read.csv(
  paste0("~/School/Waterloo/Research/RPackages/funkycells/data-raw/TNBC/",
         'generated/',"tnbc_clinical01.csv"))[-1]
TNBC_Meta <- tnbc_clinical01[,c('InternalId','AGE_AT_DX')]
colnames(TNBC_Meta) <- c('Person','Age')

# Save files
usethis::use_data(TNBC, overwrite = TRUE,compress = "xz")
# usethis::use_data(TNBC_int, overwrite = TRUE,compress = "xz")
usethis::use_data(TNBC_pheno, overwrite = TRUE,compress = "xz")
usethis::use_data(TNBC_Meta, overwrite = TRUE,compress = "xz")

