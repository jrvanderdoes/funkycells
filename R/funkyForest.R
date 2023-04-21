#' Compute Random Forest for Data with multiple PC (Along With Meta-Variables)
#'
#' This function creates a random forest model for data with PCs and
#'     meta-variables. This includes proper combination for variable importance.
#'     Recommend to users to use randomForest_CVPC in general and perhaps just
#'     this for a final model
#'
#' @param data Data.frame of outcome and predictors (PCs and meta-variables).
#'     Note, currently Unit or repeated measures should not be included.
#'     Generally use the results from getPCAData, potentially with meta-
#'     variables attached.
#' @param outcome (Optional) String indicating the outcome column name in data.
#'   Default is the first column of data.
#' @param unit (Optional) String indicating the unit column name in data.
#'   Default is the second column of data.
#' @param repeatedId (Optional) String indicating the repeated measure column
#'   name in data (if present). Default is NULL indicating no repeated measures.
#' @param nTrees (Optional) Numeric indicating the number of trees to use in the
#'     random forest model. Default is 500.
#' @param varImpPlot (Optional) Boolean indicating if variable importance plots
#'     should also be returned with the model. Default is TRUE.
#' @param metaNames (Optional) Vector with the column names of data that
#'     correspond to metavariables. Default is NULL.
#' @param keepModels (Optional) Boolean indicating if the individual models
#'     should be kept. Can get large in size. Default is TRUE.
#' @param varSelPercent Numeric in (0,1) indicating (approx) percentage of
#'     features to keep for each tree.
#' @param method (Optional) Method for rpart tree to build random forest. Default
#'   is "class".
#'
#' @return A list with  entries
#'     \enumerate{
#'         \item varImportanceData: Data.frame for variable importance
#'                                  information.
#'         \item (Optional) model: List of CART that builds the random forest model.
#'         \item (Optional) varImportancePlot: Variable importance plots.
#'     }
#' @export
#'
#' @examples
#' data <- simulatePP(cellVarData=
#'                     data.frame('stage'=c(0,1),
#'                                'A'=c(0,0),
#'                                'B'=c(1/100,1/500),
#'                                'C'=c(1/500,1/250),
#'                                'D'=c(1/100,1/100),
#'                                'E'=c(1/500,1/500)),
#'                    cellKappaData=data.frame(
#'                                'cell'=c('A','B','C','D','E'),
#'                                'clusterCell'=c(NA,'A','B','C',NA),
#'                                'kappa'=c(10,3,2,1,8)),
#'                    peoplePerStage=4,
#'                    imagesPerPerson=1)
#' pcaData <- getPCAData(data=data, repeatedUniqueId='Image',
#'                       xRange = c(0,1),  yRange = c(0,1))
#' RF1 <- funkyForest(data=pcaData[-2])
#'
#' \dontrun{
#' pcaMeta <- simulateMeta(pcaData)
#' RF2 <- funkyForest(pcaMeta[-2],
#'                                metaNames=c("randUnif","randBin","rNorm",
#'                                            "corrUnif","corrBin","corrNorm"))
#' }
funkyForest  <- function(data, outcome=colnames(data)[1],
                            unit=colnames(data)[2], repeatedId=NULL,
                            nTrees=500, varImpPlot=TRUE, metaNames=NULL,
                            keepModels=TRUE, varSelPercent=0.8, method="class"){
  # Ensure this is worthwhile
  if(length(unique(data[,outcome]))==1)
    stop('Error: Only 1 outcome in data, cannot do RF')

  # Setup
  underlyingDataAlignedFunctions <- colnames(data)[!(colnames(data) %in%
                                                       c(outcome,unit,repeatedId,metaNames))]
  for(i in 1:length(underlyingDataAlignedFunctions)){
    underlyingDataAlignedFunctions[i] <-
      substr(underlyingDataAlignedFunctions[i],1,
             utils::tail(unlist(gregexpr('_', underlyingDataAlignedFunctions[i])), n=1)-1)
  }
  underlyingDataAlignedFunctions <- underlyingDataAlignedFunctions[underlyingDataAlignedFunctions!=""]
  underlyingDataAlignedFunctions <- c(colnames(data)[colnames(data)%in%c(outcome,unit,repeatedId)],
                                      underlyingDataAlignedFunctions,metaNames)

  underlyingVars <- unique(
    underlyingDataAlignedFunctions[!(underlyingDataAlignedFunctions%in%
                                       c(outcome,unit,repeatedId))])

  data_result <- data.frame('var'=underlyingVars,
                            # xx 'splits'=0, 'giniDec'=0,
                            # xx 'varImp'=0, 'varImpCt'=0,
                            'VI'=0)
  if(keepModels) RF <- list()
  # To do CART
  for(i in 1:nTrees){
    # Bootstrap observations
    data_rf <- data[sample(1:nrow(data), nrow(data), replace = TRUE),]
    # Ensure there are multiple responses in bootstrapped sample.
    while(length(unique(data_rf[,outcome]))==1){
      data_rf <- data[sample(1:nrow(data), nrow(data), replace = TRUE),]
    }

    # Variable selection
    underlyingVarsSelected <- unique(sample(underlyingVars,
                                            varSelPercent*length(underlyingVars),
                                            replace = FALSE))
    # Build RF data. Note estimators are scrambled to avoid computational bias
    data_rf_vars <-
      data_rf[underlyingDataAlignedFunctions%in%underlyingVarsSelected]
    data_rf <- cbind(data_rf[outcome],data_rf_vars[,sample(1:ncol(data_rf_vars))])

    # Fit CART
    model <- rpart::rpart(paste0(outcome,' ~ .-',outcome), data=data_rf, method=method,
                          control=rpart::rpart.control(minsplit =1,minbucket=1, cp=0))

    # Get Var importance
    data_result <- .computeTotalVarImportance(model, data_result)

    if(keepModels) RF[[i]] <- model # Save model
  }

  # Get mean results for variable importance
  #     TEST: Use average over all or just when appears in model
  # data_result$avgGini <- data_result$giniDec / nTrees#data_result$splits
  data_result$avgVI <- data_result$VI / nTrees#data_result$varImpCt

  returnResults <- list('varImportanceData'=data_result)

  if(keepModels){
    returnResults <- append(returnResults, list('model'=RF))
  }

  # Get Variable importance
  if(varImpPlot){
    returnResults <- append(returnResults,
                            list('varImportancePlot'=.plotVariableImportance(varImportanceData=data_result)))
  }

  returnResults
}


#' Compute Total Variable Importances
#'
#' This (internal) function computes the variable importance for a tree model
#'     and combines it with totVarImportance data.frame.
#'
#' See use in funkyForest.
#'
#' @param tree Model from rpart for fitting a CART tree.
#' @param totVarImportance Data.frame indicating the total (cumulative) variable
#'     importance.
#'
#' @return Data.frame which is the updated totVarImportance adding in the new
#'     tree information.
#' @noRd
.computeTotalVarImportance <- function(tree, totVarImportance=NULL){
  tmp1 <- as.data.frame(tree$variable.importance)
  tmp2 <- data.frame('var'=rownames(tmp1),
                     'VI'=tmp1[,1])
  tmp2$var <- .getUnderlyingVariable(rownames(tmp1),returnUnique = FALSE)

  # Combine Data
  if(is.null(totVarImportance))
    totVarImportance <- data.frame('var' = unique(tmp2$var),'VI'=0)

  for(v in unique(tmp2$var)){
    totVarImportance[totVarImportance$var==v,"VI"] <-
      totVarImportance[totVarImportance$var==v,"VI"] + sum(tmp2[tmp2$var==v,'VI'])
  }

  totVarImportance
}


#' Plot Variable Importance (Ensure Data)
#'
#' This (internal) function organizes the data to plot the variable importance
#'     for a funkyForest model. One of the inputs is necessary.
#'
#' See use in funkyForest.
#'
#' Upcoming: See TEST comment
#'
#' @param varImportanceData (Optional) Data.frame for the variable importance
#'     information.
#' @param model (Optional) Random forest model from funkyForest.
#'
#' @return grid.arrange containing two ggplots
#' @noRd
.plotVariableImportance <- function(varImportanceData=NULL, model=NULL){
  if((is.null(varImportanceData) && is.null(model))||
     (!is.null(varImportanceData) && !is.null(model))){
    stop('Error: Give one varImportanceData or model')
  }

  if(is.null(varImportanceData)){
    stop('Give Var Importance')
    # cols_preds never defined
    # varImportanceData <- data.frame('var'=cols_preds,
    #                                 # 'splits'=0, 'giniDec'=0,
    #                                 # 'varImp'=0, 'varImpCt'=0,
    #                                 'VI'=0)
    #
    # for(i in 1:length(model)){
    #   varImportanceData <- .computeTotalVarImportance(model[[i]], varImportanceData)
    # }

    # Get mean results for variable importance
    #     TEST: Use average over all or just when appears in model
    # varImportanceData$avgVI <- varImportanceData$VI / nTrees#data_result$varImpCt

  }

  .plotImportance(varImportanceData)
}


#' Plot Variable Importance (Plot Data)
#'
#' This (internal) function plots the variable importance for a
#' funkyForest model.
#'
#' See use in .plotVariableImportance.
#'
#' @param varImportanceData Data.frame for the variable importance information.
#'
#' @return grid.arrange containing ggplot
#' @noRd
.plotImportance <- function(varImportanceData){
  varImportance <- ggplot2::ggplot() +
    ggplot2::geom_point(ggplot2::aes(x=stats::reorder(`var`,`avgVI`),
                                     y=avgVI/max(avgVI)),
                        data=varImportanceData) +
    ggplot2::coord_flip() +
    ggplot2::xlab(NULL) +
    ggplot2::ylim(c(0,1)) +
    ggplot2::ylab('VarImp') +
    ggplot2::theme_bw()

  varImportance
}


#' Predict Using funkyForest
#'
#' This function gets the predicted value from a funkyForest
#'     model.
#'
#' @param model funkyForest model. See funkyForest. A list of
#'     CART models from rpart.
#' @param data_pred data.frame of the data to be predicted.
#' @param type (Optional) String indicating type of analysis. Options are pred
#'     or all. The choice changes the return to best fit intended use.
#' @param data (Optional) data.frame of full data. The data used to fit the
#'     model will be extracted (by row name) to determine data the model knows.
#'
#' @return The returned data depends on type:
#'     \itemize{
#'         \item type='pred': returns a vector of the predictions
#'         \item type='all': returns a vector of the predictions
#'     }
#' @export
#'
#' @examples
#' data_pp <- simulatePP(cellVarData=
#'                       data.frame('stage'=c(0,1),
#'                                  'A'=c(0,0),
#'                                  'B'=c(1/50,1/50)),
#'                     cellKappaData=data.frame(
#'                                  'cell'=c('A','B'),
#'                                  'clusterCell'=c(NA,'A'),
#'                                  'kappa'=c(20,5)),
#'                     peoplePerStage=50,
#'                     imagesPerPerson=1,
#'                     silent=FALSE )
#' pcaData <- getPCAData(data_pp,repeatedUniqueId='Image',
#'                       xRange = c(0,1),  yRange = c(0,1), silent=FALSE)
#' RF <- funkyForest(data=pcaData[-2], nTrees = 5)
#' pred <- predict_funkyForest(model = RF$model, type='all',
#'                                  data_pred = pcaData[-2],
#'                                  data=pcaData[-2])
predict_funkyForest <- function(model, data_pred, type='all', data=NULL){
  # Verification
  #   This checks the data to see if there are any columns that are characters
  #   or factors. This will be used later to manage missing data. We obviously
  #   can't predict a group we've never seen.
  if(!is.null(data)){
    checkcols <- c()
    for(i in 1:ncol(data)){
      checkcols <- c(checkcols,
                     ifelse(is.character(data_pred[,i]) ||
                              is.factor(data_pred[,i]),i,NA))
    }
    checkcols <- checkcols[!is.na(checkcols)]
    if(length(checkcols)==0)
      data <- NULL # No need for data if nothing to check
  }

  # Setup
  #   predictions: each row is an observation prediction and columns are models
  predictions <- data.frame(matrix(nrow=nrow(data_pred),
                                   ncol=length(model)))
  data_pred_use <- data_pred # This will be used to manage missing data
  drop_pred_rows <- rep(0, nrow(data_pred)) # Rows of data_pred_use

  # Get prediction for each tree
  for(i in 1:length(model)){
    # If original data given, this will allow verification
    #   Will ensure that all characters/factors are present in the fit model or
    #   set the data to missing. Missing data is managed by model, but not
    #   previously unseen data.
    if(!is.null(data)){
      # Reset
      data_pred_use <- data_pred

      # Get rows of data used
      data_row <- names(unlist(model[[i]]$where))
      unique_rows <- unique(gsub("\\..*","", data_row )) # Since bootstrap

      # Check each variable to see if contained in fit data
      #   Only check chars and factors (nums won't have this problem)
      drop_spot <- data.frame(matrix(0,
                                     ncol=length(checkcols),
                                     nrow=nrow(data_pred_use))) # 0 keep, 1 drop
      for(j in checkcols){
        # This flags any NA or data factors/characters the model hasn't seen
        #   If any are seen, set their values to NA. This will allow the model
        #   to use secondary splits.
        if(sum(!(data_pred_use[,j] %in% data[unique_rows,j]))>0){
          for(k in 1:nrow(data_pred_use)){
            if(!(data_pred_use[k,j] %in% data[unique_rows,j]))
              data_pred_use[k,j] <- NA
          }
        }
      }
    }

    ## Predict the outcome using the model
    pred_mat <- stats::predict(model[[i]], data_pred_use)

    sols <- colnames(pred_mat)
    for(j in 1:nrow(data_pred)){
      max_pred <- which.max(pred_mat[j,])
      if(length(max_pred)!=0)
        predictions[j,i] <- sols[max_pred]
    }
  }


  # Combine predictions
  if(is.numeric(predictions)){
    warning(paste('Sorry: I haven\'t decided how to deal with modelPercents.',
                  'Only Preds will be returned'))
    type = 'pred'
    # Take Mean
    modelPredictions <- rowMeans(predictions)
    # Perhaps return L2 error instead?
  }else{
    # Vote (Randomly selects on ties)
    modelPredictions <- apply(predictions, MARGIN = 1,
                              function(x){ tmp<-.getMode(x)
                              tmp[sample(1:length(tmp),1)]})
    # Percent correctly classified
    modelPercents <- as.data.frame(t(apply(predictions, MARGIN = 1,
                                           function(x,options){
                                             sapply(options, function(option,x){
                                               sum(x==option)/length(x)},
                                               x=x)
                                           }, options=unique(data_pred[,1]))))
    colnames(modelPercents) <- unique(data_pred[,1])
  }

  # Return Results
  if(type=='pred'){
    # This returns only the predictions
    return(modelPredictions)
  }else if(type=='all'){
    return(list('PredPerc'=cbind(modelPredictions,modelPercents),
                'Acc'=sum(modelPredictions==data_pred[,1])/nrow(data_pred)))#,
    #"ROC"=multiclass.roc(data_pred[,1], modelPercents)))
  }
}


#' Get Underlying Function
#'
#' This (internal) function is used to get the underlying functino (i.e. no _PC#
#'     in it). This is used to centralize called from funkyForest.R and
#'     RandomForest_CVPC.R. If there is no PC (i.e. in meta-varaiables), the
#'     meta-variable is returned.
#'
#' See use in funkyForest.
#'
#' Note: The following line is saved temporarily.
#'     substr(test, 1, max(1,regexpr("_PC[0-9]+$", test)-1))
#'
#' @param names Vector of names to get underlying function on each
#'
#' @return Vector of underlying functions for each in name.
#' @noRd
.getUnderlyingVariable <- function(names, returnUnique=TRUE){
  if(returnUnique)
    return(unique(stringr::str_remove(names,'(_PC)[0-9]+$')))
  return(stringr::str_remove(names,'(_PC)[0-9]+$'))
}
