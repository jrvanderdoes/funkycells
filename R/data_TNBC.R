#' Triple Negative Breast Cancer Data
#'
#' A funky model ready set of principle components from K functions based on
#'   triple negative breast cancer data from patients. The original data was
#'   proteins are coded as T/F values. Additionally, the age meta-variable was
#'   added.
#'
#' @format ## `TNBC`
#' A data frame with 33 rows and 1398 columns:
#' \describe{
#'   \item{Class}{Outcome of each patient}
#'   \item{Person}{Person for each image}
#'   \item{NA_Si_PC1 through tumerYN_tumerYN_PC3}{Principle components of the K
#'     functions for the named interactions}
#'   \item{age}{Meta-variable for patient age}
#'   ...
#' }
#' @source <https://www.angelolab.com/mibi-data>
"TNBC"

#' Triple Negative Breast Cancer Phenotypes
#'
#' Data of triple negative breast cancer biopsies from patients, with. Additionally, the age meta-variable was
#'   added.
#'
#' @format ## `TNBC_pheno`
#' A data frame with 170,171 rows and 5 columns:
#' \describe{
#'   \item{Class}{Outcome of each patient}
#'   \item{Person}{Person for which each cell is related}
#'   \item{cellx, celly}{The x-y coordinates of the cell}
#'   \item{Phenotype}{The classified phenotype for the cecll}
#'   ...
#' }
#' @source <https://www.angelolab.com/mibi-data>
"TNBC_pheno"


#' Triple Negative Breast Cancer Meta Data
#'
#' Data of triple negative breast cancer biopsies from patients, with. Additionally, the age meta-variable was
#'   added.
#'
#' @format ## `TNBC_Meta`
#' A data frame with 33 rows and 2 columns:
#' \describe{
#'   \item{Person}{Person for each image}
#'   \item{Age}{Meta-variable for patient age}
#'   ...
#' }
#' @source <https://www.angelolab.com/mibi-data>
"TNBC_Meta"
