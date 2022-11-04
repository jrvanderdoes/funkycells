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
