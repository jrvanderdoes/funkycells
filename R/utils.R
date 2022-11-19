#' Get Mode of vector
#'
#' @param x Vector of some items
#'
#' @return The mode(s) of the vector
#' @export
#'
#' @examples
#' .getMode(c(1,2,3,1))
#' .getMode(c('A','B','C','A','B'))
#' .getMode(c(1,2,3,'A','A'))
#' .getMode(c(1,2,3,'A','A',2,2,2))
.getMode <- function(x) {
  a <- table(x)
  mode <- names(a)[a == max(a)]
  if(is.numeric(x))
    mode <- as.numeric(mode)
  mode
  #uniqx <- unique(na.omit(x))
  #uniqx[which.max(tabulate(match(x, uniqx)))]
}


#' Convert a list to data.frame
#'
#' This (internal) function converts a list into a data.frame, either by row or
#'     column binding. Some error handling exists, but care should be taken.
#'
#' @param data_list List of data.frames to be combined
#' @param typeBind String of row or col indicating if the lists should be
#'     combined using the rows or columns.
#' @param na.omit (Optional) Boolean that drops with any NAs if TRUE. The
#'     default is FALSE.
#'
#' @return Data.frame with the data from the list.
#' @export
#'
#' @examples
#' # See code for .getPCs or simulatePP. This is not an outward function so
#' #     won't be viewable.
.convertList2Dataframe <- function(data_list,
                                   typeBind=c('row','col'),
                                   na.omit=F){
  if(length(typeBind)!=1 || !(typeBind%in%c('row','col')))
    stop('Error: Select row or col for typeBind')

  data_df <- data.frame()

  if(typeBind=='row'){
    # No error checking. See simulatePP or .generateCSRPatterns
    for(i in 1:length(data_list)){
      data_df <- rbind(data_df,data_list[[i]])
    }
  }else if(typeBind=='col'){
    # See if first column is DF or not. Can occur based on getPCAData
    if(dim(data_list[[1]])[[1]]>1){
      data_df <- data_list[[1]]
    }else{
      data_df <- data.frame('V1'=t(data_list[[1]]))
    }

    if(length(data_list)!=1){
      for(ii in 2:length(data_list)){
        if(dim(data_list[[ii]])[[1]]>1){
          data_df <- cbind(data_df,data_list[[ii]])
        }else{
          data_df <- cbind(data_df, t(data_list[[ii]]))
        }
      }
    }
  }

  if(na.omit)
    data_df <- na.omit(data_df)

  data_df
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
#' @export
#'
#' @examples
#' .getFolds(1:10,3)
.getFolds <- function(x,K) {
  if(K>length(x))
    warning(paste0('Warning: Only ',length(x),' chunks are used due to data'))
  if(K<2)
    return(sample(x))
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
#' @export
#'
#' @examples
#' .specify_decimal(10.123,1)
#' .specify_decimal(10.123,3)
#' .specify_decimal(10.123,5)
.specify_decimal <- function(x, k) {trimws(format(round(x, k), nsmall=k))}
