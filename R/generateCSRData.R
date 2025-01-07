#' Generate Completely Spatially Random Point Process Data
#'
#' This (internal) function generates a completely spatially random (CSR)
#'     process on some region. Some Poisson number of points are placed
#'     uniformly in the region by both x and y.
#'
#' See use in simulatePP_Mult_PCA_Meta.
#'
#' @param xRange,yRange (Optional) Two value numeric vector indicating the size
#'  of the region in the x/y-direction. Default is c(0,1).
#' @param kappa (Optional) Numeric influencing the number of points. Kappa times
#'     area defines the mean of a Poisson random variable which is randomly
#'     drawn to indicate the number of agents in the area. Default is 25.
#' @param requireOne (Optional) Boolean that indicating if at least one point
#'     must be placed. This can affect the Poisson distribution, but is valuable
#'     in ensuring a non-empty region. Default is TRUE.
#' @param type (Optional) Value (often character) indicating the label the
#'     points in the process should be given. Default is A.
#'
#' @return Data.frame with x, y, and type specified. Each row is a new
#'     point.
#'
#' @keywords internal
#' @noRd
.generateCSRData <- function(xRange = c(0, 1), yRange = c(0, 1),
                             kappa = 25, requireOne = TRUE, type = "A") {
  area <- (xRange[2] - xRange[1]) * (yRange[2] - yRange[1])
  intensityValue <- kappa * area
  numPts <- stats::rpois(1, intensityValue)
  if (numPts == 0) {
    if (!requireOne) {
      return() # TEST: Not sure how this will work
    } else if (intensityValue == 0) {
      numPts <- 1 # 0 intensityValue always gives 0 so just give 1
    } else {
      idx <- 0
      while (numPts == 0) {
        numPts <- stats::rpois(1, intensityValue)
        idx <- idx + 1
        if (idx > 5) numPts <- 1 # Avoid long loops if intensityValue is small
      }
    }
  }

  pointPattern <- data.frame(
    "x" = stats::runif(numPts, min = xRange[1], max = xRange[2]),
    "y" = stats::runif(numPts, min = yRange[1], max = yRange[2]),
    "type" = type
  )

  if (nrow(unique(pointPattern)) != nrow(pointPattern)) {
    warning("Points placed on top of each other, so dropped (not replaced)")
    pointPattern <- unique(pointPattern)
  }

  pointPattern
}
