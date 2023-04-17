#' Triple Negative Breast Cancer Data - Indicators
#'
#' A subset of triple negative breast cancer data from patients indicating
#'     cell information for each person. The proteins are coded as T/F values.
#'
#' @format ## `TNBC`
#' A data frame with 170,171 rows and 60 columns:
#' \describe{
#'   \item{Class}{Outcome of each patient}
#'   \item{Person}{Person for each image}
#'   \item{cellLabelInImage}{Label of cell}
#'   \item{cellx, celly}{x and y location for each cell}
#'   \item{C,Na,Si,P,Ca,Fe,dsDNA,Vimentin,SMA,Background,B7H3,FoxP3,Lag3,CD4,
#'         CD16,CD56,OX40,PD1,CD31,PD.L1,EGFR,Ki67,CD209,CD11c,CD138,CD163,CD68,
#'         CSF.1R,CD8,CD3,IDO,Keratin17,CD63,CD4RO,CD20,p53,Beta.catenin,HLA.DR,
#'         CD11b,CD45,H3K9ac,Pan.Keratin,H3K27me3,phospho.S6,MPO,Keratin6,
#'         HLA_Class_1,Ta,Au,tumerYN}{Categories for marked point patterns}
#'   ...
#' }
#' @source <https://www.angelolab.com/mibi-data>
"TNBC"


#' Triple Negative Breast Cancer Data - Intensities
#'
#' A subset of triple negative breast cancer data from patients indicating
#'     cell information for each person. The proteins are coded as intensity
#'     values.
#'
#' @format ## `TNBC`
#' A data frame with 170,171 rows and 60 columns:
#' \describe{
#'   \item{Class}{Outcome of each patient}
#'   \item{Person}{Person for each image}
#'   \item{cellLabelInImage}{Label of cell}
#'   \item{cellx, celly}{x and y location for each cell}
#'   \item{C,Na,Si,P,Ca,Fe,dsDNA,Vimentin,SMA,Background,B7H3,FoxP3,Lag3,CD4,
#'         CD16,CD56,OX40,PD1,CD31,PD.L1,EGFR,Ki67,CD209,CD11c,CD138,CD163,CD68,
#'         CSF.1R,CD8,CD3,IDO,Keratin17,CD63,CD4RO,CD20,p53,Beta.catenin,HLA.DR,
#'         CD11b,CD45,H3K9ac,Pan.Keratin,H3K27me3,phospho.S6,MPO,Keratin6,
#'         HLA_Class_1,Ta,Au,tumerYN}{Categories for marked point patterns}
#'   ...
#' }
#' @source <https://www.angelolab.com/mibi-data>
"TNBC_int"


#' Triple Negative Breast Cancer Data - Phenotypes
#'
#' A subset of triple negative breast cancer data from patients indicating
#'     cell information for each person. The phenotypes are given.
#'
#' @format ## `TNBC`
#' A data frame with 170,171 rows and 5 columns:
#' \describe{
#'   \item{Class}{Outcome of each patient}
#'   \item{Person}{Person for each image}
#'   \item{cellx, celly}{x and y location for each cell}
#'   \item{Phenotype}{Categories for cell phenotype}
#'   ...
#' }
#' @source <https://www.angelolab.com/mibi-data>
"TNBC_pheno"


#' Triple Negative Breast Cancer Data
#'
#' A subset of triple negative breast cancer data containing the meta-data
#'     for each person.
#'
#' @format ## `TNBC_Meta`
#' A data frame with 33 rows and 2 columns:
#' \describe{
#'   \item{Person}{Person which the meta-data is related to}
#'   \item{Age}{Age for the related person}
#'   ...
#' }
#' @source <https://www.angelolab.com/mibi-data>
"TNBC_Meta"
