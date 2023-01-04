#' Get K function
#'
#' This function computes the K function between the two agents for each unit,
#'     potentially averaging over repeated measures.
#'
#' @param data Dataframe with column titles for at least x, y, agents,
#'     and unit. For consistency (and avoiding errors), use that order.
#'     Additionally, repeatedUniqueId can be added.
#' @param agents Two value vector indicating the two agents to use for the K
#'     function, the first to the second. These should be in the unit column.
#' @param unit String of the column name in data indicating a unit or base
#'     thing. Note this unit may have repeated measures.
#' @param repeatedUniqueId (optional) String of the column name in data
#'     indicating the unique ID when using repeated measures.
#' @param rCheckVals (optional) numeric vector indicating the radius to check.
#'    Note, if note specified, this could take a lot of memory, particularly
#'    with many units and repeated measures.
#' @param xRange (optional) two value numeic vector indicating the min and max
#'     x values. Note this is re-used for all images. The default just takes
#'     the min and max x from each image. This allows different sized images,
#'     but note that the edges are defined by some cell.
#' @param yRange (optional) two value numeic vector indicating the min and max
#'     y values. Note this is re-used for all images. The default just takes
#'     the min and max y from each image. This allows different sized images,
#'     but note that the edges are defined by some cell.
#' @param edgeCorrection (optional) String indicating type of edgeCorrection to
#'     use in spatStat for computation of K functions.
#'
#' @return data.frame with the first column being the checked radius and the
#'     rest relating to the K function for each unit at those points. If a K
#'     function could not be computed, perhaps due to lack of data, an NA is
#'     returned for the K function.
#' @export
#'
#' @examples
#' data1 <- simulatePP(cellVarData=
#'                       data.frame('stage'=c(0,1,2),
#'                       'A'=c(0,0,0),
#'                       'B'=c(1/100,1/500,1/500)),
#'                    cellKappaData=data.frame(
#'                       'cell'=c('A','B'),
#'                       'clusterCell'=c(NA,'A'),
#'                       'kappa'=c(20,5)))
#' # Select the data
#' KData <- data1[data1$Person=='p1',colnames(data1)!='Stage']
#' KFunction <- getKFunction(agents=c('A','B'), unit='Person',
#'              repeatedUniqueId='Image',
#'              data=KData,
#'              rCheckVals=seq(0,0.25,0.01),
#'              xRange=c(0,1),yRange=c(0,1),
#'              edgeCorrection="isotropic")
#' plot(KFunction, type='l')
getKFunction <- function(data, agents, unit,
                         repeatedUniqueId=NULL,
                         rCheckVals=NULL,
                         xRange=NULL, yRange=NULL,
                         edgeCorrection="isotropic"){

  ## Must save as a list to ensure multiple sizes are handled (rCheckVals=NULL)
  matrix_r <- matrix_K <- list()
  units <- unique(data[,unit])

  for(i in 1:length(units)){
    # Unit Data
    data_unit <- data[data[,unit]==units[i],]
    data_unit <- data_unit[,colnames(data_unit)!=unit]

    ## Handle repeated measures (including if none)
    repeats <- NA
    if(!is.null(repeatedUniqueId))
      repeats <- unique(data_unit[,repeatedUniqueId])
    # Set up var
    matrix_K_tmp <- matrix_r_tmp <- NULL
    agent1Counts <- rep(NA, length(repeats))

    for(j in 1:length(repeats)){
      # Unit's Repeat Data
      if(length(repeats)>1){
        data_repeat <- data_unit[data_unit[,repeatedUniqueId]==repeats[j],]
      }else{
        data_repeat <- data_unit
      }
      if(!is.null(repeatedUniqueId))
        data_repeat <- data_repeat[,colnames(data_unit)!=repeatedUniqueId]

      # Make marks all factors
      for(k in 3:ncol(data_repeat)){
        if(is.character(data_repeat[[k]]))
          data_repeat[[k]] <- as.factor(data_repeat[[k]])
      }

      # Define xRange and yRange for this image. Take given if exists
      if(is.null(xRange)){
        xRange_rm <- c(min(data_repeat[,1]),max(data_repeat[,1]))
      } else{
        xRange_rm <- xRange
      }
      if(is.null(yRange)) {
        yRange_rm <- c(min(data_repeat[,2]),max(data_repeat[,2]))
      } else{
        yRange_rm <- yRange
      }

      # Define the marked point pattern using the data set for an individual person.
      # For this data the person data has had the marks listed as True/False, as factors.
      # as.ppp interprets the first two columns as x-y coordinates and remaining colums as marks
      data_ppp <- spatstat.geom::as.ppp(data_repeat,
                                        W=spatstat.geom::as.owin(c(xRange_rm,yRange_rm)))

      ## Marked since additional properties
      if(ncol(data_repeat)>3){
        # Apply Kmulti, which generates Kcross functions for subsets of cells that we can
        # define by functions of the cell marks.
        # Since we have Marks as T/F, we set the functions f1 and f2 to return the
        # Subsets of two given marks that are true. The marks can be the same.

        # Because we may be defining different window sizes for each person
        # we return the vector r as well, as it will be different for each patient.
        f1 = function(X){which(spatstat.geom::marks(X)[,agents[1]]==TRUE)}
        f2 = function(X){which(spatstat.geom::marks(X)[,agents[2]]==TRUE)}
        K <- tryCatch({spatstat.core::Kmulti(data_ppp,f1,f2,
                             correction = edgeCorrection,
                             r=rCheckVals)
        }, error=function(e){
          list(NA,NA,NA)
        })

        # Save weights
        agent1Counts[j] <- nrow(data_repeat[data_repeat[,agents[1]],])
      } else{
        K <- tryCatch({
          spatstat.core::Kcross(data_ppp,
                           agents[1], agents[2],
                           correction = edgeCorrection,
                           r=rCheckVals)
        }, error=function(e){
          list(NA,NA,NA)
        })

        # Save weights
        agent1Counts[j] <- nrow(data_repeat[data_repeat[,3]==agents[1],])
      }

      ## AlignK and r
      result <- .alignKr(K = matrix_K_tmp, newK = K[[3]],
                         r = matrix_r_tmp, newr = K[[1]])
      matrix_K_tmp <- cbind(result[[1]],result[[2]])
      matrix_r_tmp <- result[[3]]

    }
    # Drop any NA
    dropIdxs <- which(colSums(is.na(matrix_K_tmp))!=0)
    if(length(dropIdxs)>0){
      matrix_K_tmp <- matrix_K_tmp[,-dropIdxs]
      agent1Counts <- agent1Counts[-dropIdxs]
    }

    if(length(matrix_K_tmp)!=0){
      ## Take weighted average
      matrix_K[[i]] <- (as.matrix(matrix_K_tmp) %*% agent1Counts)/sum(agent1Counts)
      matrix_r[[i]] <- data.frame(matrix_r_tmp)
    } else{
      matrix_K[[i]] <- data.frame(NA)
      matrix_r[[i]] <- data.frame(NA)
    }
  }

  # Convert into a single df
  data_return <- .rK2DF(K_list = matrix_K, r_list = matrix_r)

  data_return
}


#' Align K and r matrices
#'
#' This (internal) function makes a master r vector and ensures all K functions
#'     are defined at each r.
#'
#' @param K Data.frame of (potentially) columns of K functions for each previous
#'     unit (i.e. the previous K matrix). These are evaluated at the values in r.
#' @param newKData.frame of column of K functions for current unit, evaluated
#'     at the values in newr.
#' @param r Data.frame of one column for the evaluated r values in K (i.e. the
#'     previous r matrix).
#' @param newr  Data.frame of one column for the evaluated r values in newK.
#'
#' @return list of three data.frames
#'     1. K which is K evaluated at the (potentially) new r
#'     2. newK which is newK evaluated at the (potentially) new r
#'     3. r which is the new master vector of evaluated radius r values
#'
#' @examples
#' # See code for getKFunction This is not an outward function so won't be
#' #     viewable.
.alignKr <- function(K, newK, r, newr){
  # On the first loop (i.e. no existing data)
  if(is.null(K) && is.null(r)){
    # If failure in K function computation
    if(length(newK)==1 && length(newr)==1 &&
       is.na(newK) && is.na(newr)){
          ## Check, this may not work
          return(list('K'=NA,
                      'newK'=NULL,
                      'r'=NA))
    }else{
          return(list('K'=NULL,
                      'newK'=newK,
                      'r'=newr))
    }
  }

  # Failure in K function computation
  if(length(newK)==1 && length(newr)==1 &&
     is.na(newK) && is.na(newr)){
        return(list('K'=K,
                    'newK'=rep(NA, length(r)),
                    'r'=r))
  }

  bigr <- na.omit(unique(c(r, newr)))
  bigr <- bigr[order(bigr)]

  K_ret <- merge(data.frame('r'=bigr),
               data.frame('r'=r,K),
               by=c("r"), all.x=TRUE)
  K_ret <- tidyr::fill(K_ret,-r,.direction='down')#K_ret %>% fill (-r,.direction='down')

  newK_ret <- merge(data.frame('r'=bigr),
                data.frame('r'=newr,'K'=newK),
                by=c("r"), all.x=TRUE)
  newK_ret <- tidyr::fill(newK_ret, K,.direction='down')#newK_ret %>% fill(K,.direction='down')

  list('K'=K_ret[-1],
       'newK'=newK_ret[-1],
       'r'=bigr)
}


#' Convert r and K to a Data.frame
#'
#' This converts the r and K matrices from different units into a single
#'     data.frame.
#'
#' @param K_list List of data.frames from each unit. Each data.frame has the
#'     K function for each unit evaluated at the corresponding r in r_list.
#' @param r_list List of data.frames from each unit. Each data.frame has the
#'     evaluted r values for each unit used in the corresponding K in K_list.
#'
#' @return Data.frame with the first column being the evaluated r and the rest
#'     being the evaluted K functions for each unit.
#'
#' @examples
#' # See code for getKFunction This is not an outward function so won't be
#' #     viewable.
.rK2DF <- function(K_list, r_list){
  # Define r
  r <- c()
  for(i in 1:length(r_list)){
    r <- c(r, r_list[[i]][,1])
  }
  r <- unique(na.omit(r))

  # Define K and set up
  data_ret <- data.frame('r'=r[order(r)])
  for(i in 1:length(K_list)){
    tmp <- data.frame('r'=r_list[[i]][[1]],
                      'K'=K_list[[i]][,1])
    colnames(tmp) <- c('r',paste0('K',i))
    data_ret <- merge(data_ret, tmp,
                      by=c("r"), all.x=TRUE)
  }
  data_ret <- tidyr::fill(data_ret, -r,.direction='down')#data_ret %>% fill (-r,.direction='down')

  data_ret
}
