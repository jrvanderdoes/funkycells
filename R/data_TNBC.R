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
#'   \item{Person}{Person for each image}
#'   \item{Class}{Outcome of each patient}
#'   \item{NA_Si_PC1 through tumerYN_tumerYN_PC3}{Principle components of the K
#'     functions for the named interactions}
#'   \item{age}{Meta-variable for patient age}
#'   ...
#' }
#' @source <https://www.angelolab.com/mibi-data>
"TNBC"

