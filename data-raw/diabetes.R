## code to prepare `diabetes` dataset goes here
library(spicyR)

getDiabetesData <- function(){
  # This function get the diabetes data and organizes it from spicyR

  #####
  # Get Data
  data("diabetesData")

  tmp <- as.data.frame(cellSummary(diabetesData))
  tmp1 <- as.data.frame(imagePheno(diabetesData))

  #####
  # Organize Data
  data <- data.frame('Person'=rep(NA,length(tmp$imageID)),
                     'Image'=NA,
                     'x'=NA, 'y'=NA,
                     'cellType'=NA,
                     'Stage'=NA)

  data$Image <- tmp$imageID
  data$x <- tmp$x
  data$y <- tmp$y
  data$cellType <- as.character(tmp$cellType)

  for(i in 1:length(data$Person)){
    data$Person[i] <-
      tmp1[as.character(tmp1$imageID)==as.character(data$Image[i]), 'case']
    data$Stage[i] <-
      as.character(tmp1[as.character(tmp1$imageID)==as.character(data$Image[i]), 'stage'])
  }

  data
}
diabetes <- getDiabetesData()

usethis::use_data(diabetes, overwrite = TRUE,compress = 'bzip2')

