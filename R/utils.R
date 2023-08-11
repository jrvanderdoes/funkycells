#' Get Mode of vector
#'
#' This ignores NA. See below internal for other things which may be
#'  worth trying later.
#'
#' Internal note: Few lines temporarily saved:
#'           uniqx <- unique(na.omit(x))
#'           uniqx[which.max(tabulate(match(x, uniqx)))]
#'
#' @param x Vector of some items
#'
#' @return The mode(s) of the vector
#'
#' @examples
#' m1 <- .getMode(c(1, 2, 3, 1))
#' m2 <- .getMode(c("A", "B", "C", "A", "B"))
#' m3 <- .getMode(c(1, 2, 3, "A", "A"))
#' m4 <- .getMode(c(1, 2, 3, "A", "A", 2, 2, 2))
#' @noRd
.getMode <- function(x) {
  a <- table(x)
  mode <- names(a)[a == max(a)]
  if (is.numeric(x)) {
    mode <- as.numeric(mode)
  }
  mode
}


#' Merge Lists into Base Data.frame
#'
#' This (internal) function appends items in list to a data.frame by the given
#'     column names.
#'
#' See usage in getKsPCAData.
#'
#' @param df Data.frame with a column named by baseCol
#' @param lists Lists of data.frames with each containing listCol column
#'              matching the values in df$baseCol.
#' @param dfCol String with column name for matching column in df
#' @param listsDFCol String with column name for matching column in lists
#'
#' @return Data.frame df with values from lists appended
#' @noRd
.mergeListsToDF <- function(df, lists, dfCol, listsDFCol) {
  for (i in 1:length(lists)) {
    df <- merge(
      x = df, by.x = dfCol,
      y = lists[[i]], by.y = listsDFCol,
      sort = FALSE
    )[, union(names(df), names(lists[[i]])[!(names(lists[[i]]) %in% listsDFCol)])]
  }

  df
}

#' Create K folds
#'
#' This (internal) function creates K folds of x. Note these folds are scrambled
#'     data, but then ordered in each fold.
#'
#' @param x Vector of data (or often indices)
#' @param chunksN Numeric for the number of chunks, or folds, desired
#'
#' @return a list with length chunksN, with each containing approximately the
#'     same number of points
#'
#' @examples
#' gf <- .getFolds(1:10, 3)
#' @noRd
.getFolds <- function(x, K) {
  if (K > length(x)) {
    warning(paste0("Warning: Only ", length(x), " chunks are used due to data"))
  }
  if (K < 2) {
    return(sample(x))
  }
  split(x, sample(cut(x, K, labels = FALSE)))
}


#' Specify Decimals
#'
#' This (internal) function ensures a numeric displays a given number of
#'     decimals. This is different from "round" in that it will add zeros as
#'     necessary.
#'
#' @param x Numeric to display k decimals of.
#' @param k Numeric indicating the number of decimals to display
#'
#' @return String with numeric x to k decimal places
#'
#' @examples
#' .specify_decimal(10.123, 1)
#' .specify_decimal(10.123, 3)
#' .specify_decimal(10.123, 5)
#' @noRd
.specify_decimal <- function(x, k) {
  trimws(format(round(x, k), nsmall = k))
}
