#' Compute Random Forest for Data with multiple PC (Along With Meta-Variables)
#'
#' This (internal) function creates a random forest model for data with PCs and
#'     meta-variables. This includes proper combination for variable importance.
#'     This is internal because we want users to use randomForest_CVPC.
#'
#' Upcoming: See TEST comment
#'           Drop parts to make smaller memory impact.
#'               https://rstudio-pubs-static.s3.amazonaws.com/381812_7b5343e0f1a7403fba31ff7d45cf6d88.html
#'
#' @param data Data.frame of outcome and predictors (PCs and meta-variables).
#'     Note, currently Unit or repeated measures should not be included.
#'     Generally use the results from getPCAData, potentially with meta-
#'     variables attached.
#' @param outcome (Optional) String indicating the outcome column name in data
#' @param nTrees (Optional) Numeric indicating the number of trees to use in the
#'     random forest model. Default is 1000.
#' @param varImpPlot (Optional) Boolean indicating if variable importance plots
#'     should also be returned with the model. Default is TRUE.
#' @param metaNames (Optional) Vector with the column names of data that
#'     correspond to metavariables. Default is NULL.
#'
#' @return A list with two (three) entries
#'     1. model: List of CART that builds the andom forest model
#'     2. varImportanceData: Data.frame for variable importance information
#'     3. (Optional) varImportancePlot: Variable importance plots.
#' @export
#'
#' @examples
#' data <- simulatePP()
#' pcaData <- getPCAData(data=data, repeatedUniqueId='Image',
#'                       xRange = c(0,1),  yRange = c(0,1), silent=F)
#' RF1 <- .computeRandomForest_PC(data=pcaData[-2])
#'
#' pcaMeta <- simulateMeta(pcaData)
#' RF2 <- .computeRandomForest_PC(pcaMeta[-2],
#'                                metaNames=c("randUnif","randBin","rNorm",
#'                                            "corrUnif","corrBin","corrNorm"))
.computeRandomForest_PC  <- function(data, outcome=colnames(data)[1],
                              unit=colnames(data)[2], repeatedId=NULL,
                              nTrees=1000, varImpPlot=T, metaNames=NULL){
  # Ensure this is worthwhile
  if(length(unique(data[,outcome]))==1)
    stop('Error: Only 1 outcome in data, cannot do RF')

  # Setup
  underlyingDataAlignedFunctions <- colnames(data)[!(colnames(data) %in%
                                                      c(outcome,unit,repeatedId,metaNames))]
  for(i in 1:length(underlyingDataAlignedFunctions)){
    underlyingDataAlignedFunctions[i] <-
      substr(underlyingDataAlignedFunctions[i],1,
        tail(unlist(gregexpr('_', underlyingDataAlignedFunctions[i])), n=1)-1)
  }
  underlyingDataAlignedFunctions <- underlyingDataAlignedFunctions[underlyingDataAlignedFunctions!=""]
  underlyingDataAlignedFunctions <- c(colnames(data)[colnames(data)%in%c(outcome,unit,repeatedId)],
                                      underlyingDataAlignedFunctions,metaNames)

  underlyingVars <- unique(
    underlyingDataAlignedFunctions[!(underlyingDataAlignedFunctions%in%
                                     c(outcome,unit,repeatedId))])

  data_result <- data.frame('var'=underlyingVars,
                            'splits'=0,
                            'giniDec'=0,
                            'varImp'=0,
                            'varImpCt'=0,
                            'vi'=0)
  RF <- list()
  # To do CART
  varSelPercent <- 0.8 # Percent of variable selection
  for(i in 1:nTrees){
    # Bootstrap observations
    data_rf <- data[sample(1:nrow(data), nrow(data), replace = T),]
    # Ensure there are multiple responses in bootstrapped sample.
    while(length(unique(data_rf[,outcome]))==1){
      data_rf <- data[sample(1:nrow(data), nrow(data), replace = T),]
    }

    # Variable selection
    underlyingVarsSelected <- unique(sample(underlyingVars,
                    varSelPercent*length(underlyingVars), replace = F))
    # Build RF data. Note estimators are scrambled to avoid computational bias
    data_rf_vars <-
          data_rf[underlyingDataAlignedFunctions%in%underlyingVarsSelected]
    data_rf <- cbind(data_rf[outcome],data_rf_vars[,sample(1:ncol(data_rf_vars))])

    # Fit CART
    model <- rpart::rpart(paste0(outcome,' ~ .'), data=data_rf, method="class",
                   control=rpart::rpart.control(minsplit =1,minbucket=1, cp=0))
    data_result <- .computeTotalVarImportance(model, data_result)

    RF[[i]] <- model # Save model
  }

  # Get mean results for variable importance
  #     TEST: Use average over all or just when appears in model
  data_result$avgGini <- data_result$giniDec / nTrees#data_result$splits
  data_result$avgVI <- data_result$varImp / nTrees#data_result$varImpCt


  # Get Variable importance
  if(varImpPlot){
    return(list('model'=RF, 'varImportanceData'=data_result,
                'varImportancePlot'=.plotVariableImportance(varImportanceData=data_result)))
  }

  list('model'=RF, 'varImportanceData'=data_result)
}


#' Compute Total Variable Importances
#'
#' This (internal) function computes the variable importance for a tree model
#'     and combines it with totVarImportance data.frame.
#'
#' @param tree Model from rpart for fitting a CART tree.
#' @param totVarImportance Data.frame indicating the total (cumulative) variable
#'     importance.
#'
#' @return Data.frame which is the updated totVarImportance adding in the new
#'     tree information.
#' @export
#'
#' @examples
#' # See code for .computeRandomForest_PC. This is not an outward function so won't be
#' #     viewable.
.computeTotalVarImportance <- function(tree, totVarImportance){

  # Find Error at nodes (Gini is 1-p_1^2-p_2^2)
  frame <- tree$frame
  frame[['gini']] = 1 - (frame[['dev']] / frame[['n']])^2 -
    (1 - frame[['dev']] / frame[['n']])^2
  #frame[,c('var','n','dev','gini')]

  # Get Decrease for each split
  frame[['improve']] = NA
  for (j in 1:nrow(frame)) {
    name <- frame[j,'var']
    if (name == '<leaf>') next

    # Get name of bigger group (if exists)
    name_sub <- .getUnderlyingVariable(name)

    # Get Gini improvement
    ind = which(rownames(frame) %in%
                  (as.numeric(rownames(frame)[j])*2+c(0,1)))
    frame[j,'improve'] = frame[j,'n']*frame[j,'gini'] - frame[ind[1],'n']*
      frame[ind[1],'gini'] - frame[ind[2],'n']*frame[ind[2],'gini']

    totVarImportance[totVarImportance$var==name_sub,'giniDec'] <-
      totVarImportance[totVarImportance$var==name_sub,'giniDec'] +
      frame[j,'improve']
    # Count split by variable
    totVarImportance[totVarImportance$var==name_sub,'splits'] <-
      totVarImportance[totVarImportance$var==name_sub,'splits'] + 1
  }

  # Get and count Var Import
  for(j in 1:length(tree$variable.importance)){
    name <- names(tree$variable.importance[j])

    # Get name of bigger group (if exists)
    name_sub <- .getUnderlyingVariable(name)

    totVarImportance[totVarImportance$var==name_sub,'varImp'] <-
      totVarImportance[totVarImportance$var==name_sub,'varImp'] +
      tree$variable.importance[[j]]
    totVarImportance[totVarImportance$var==name_sub,'varImpCt'] <-
      totVarImportance[totVarImportance$var==name_sub,'varImpCt'] +
      1
  }

  totVarImportance
}


#' Plot Variable Importance (Ensure Data)
#'
#' This (internal) function organizes the data to plot the variable importance
#'     for a .computeRandomForest_PC model. One of the inputs is necessary.
#'
#' Upcoming: See TEST comment
#'
#' @param varImportanceData (Optional) Data.frame for the variable importance
#'     information.
#' @param model (Optional) Random forest model from .computeRandomForest_PC.
#'
#' @return grid.arrange containing two ggplots
#' @export
#'
#' @examples
#' # See code for .computeRandomForest_PC. This is not an outward function so won't be
#' #     viewable.
.plotVariableImportance <- function(varImportanceData=NULL, model=NULL){
  if((is.null(varImportanceData) && is.null(model))||
     (!is.null(varImportanceData) && !is.null(model))){
    stop('Error: Give one varImportanceData or model')
  }
  if(is.null(varImportanceData)){

    varImportanceData <- data.frame('var'=cols_preds,
                              'splits'=0,
                              'giniDec'=0,
                              'varImp'=0,
                              'varImpCt'=0,
                              'vi'=0)

    for(i in 1:length(model)){
      varImportanceData <- .computeTotalVarImportance(model[[i]], varImportanceData)
    }

    # Get mean results for variable importance
    #     TEST: Use average over all or just when appears in model
    varImportanceData$avgGini <- varImportanceData$giniDec / nTrees#data_result$splits
    varImportanceData$avgVI <- varImportanceData$varImp / nTrees#data_result$varImpCt

  }

  .plotImportance(varImportanceData)
}


#' Plot Variable Importance (Plot Data)
#'
#' This (internal) function plots the variable importance for a
#' .computeRandomForest_PC model.
#'
#' @param varImportanceData Data.frame for the variable importance information.
#'
#' @return grid.arrange containing two ggplots
#' @export
#'
#' @examples
#' # See code for .plotVariableImportance. This is not an outward function so won't be
#' #     viewable.
.plotImportance <- function(varImportanceData){

  avgGini <- ggplot2::ggplot() +
    ggplot2::geom_point(ggplot2::aes(x=reorder(var,avgGini), y=avgGini/max(avgGini)),
                        data=varImportanceData) +
    ggplot2::coord_flip() +
    ggplot2::xlab(NULL) +
    ggplot2::ylim(c(0,1)) +
    ggplot2::ylab('Gini') +
    ggplot2::theme_bw()
  avgVI <- ggplot2::ggplot() +
    ggplot2::geom_point(ggplot2::aes(x=reorder(var,avgVI), y=avgVI/max(avgVI)),
                        data=varImportanceData) +
    ggplot2::coord_flip() +
    ggplot2::xlab(NULL) +
    ggplot2::ylim(c(0,1)) +
    ggplot2::ylab('VarImp') +
    ggplot2::theme_bw()

  varImportance <- gridExtra::arrangeGrob(avgGini,avgVI,
                                layout_matrix = rbind(c(1,2)),
                                bottom = 'Variable Importance - Percent of Max')

  varImportance
}


#' Predict Using RandomForest_PC
#'
#' This (internal) function gets the predicted value from a RandomForest_PC
#'     model.
#'
#' @param model RandomForest_PC model. See .computeRandomForest_PC. A list of
#'     CART models from rpart.
#' @param data_pred data.frame of the data to be predicted.
#' @param type (Optional) String indicating type of analysis. Options are pred
#'     or all. The choice changes the return to best fit intended use.
#' @param data (Optional) data.frame of full data. The data used to fit the
#'     model will be extracted (by row name) to determine data the model knows.
#'
#' @return The returned data depends on type:
#'     1. type='pred': returns a vector of the predictions
#'     2. type='all': returns a vector of the predictions
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
#'                     peoplePerStage=100,
#'                     imagesPerPerson=1,
#'                     reduceEdge=0.025,
#'                     silent=F )
#' pcaData <- getPCAData(data_pp,repeatedUniqueId='Image',
#'                       xRange = c(0,1),  yRange = c(0,1), silent=F)
#' RF <- .computeRandomForest_PC(data=pcaData[-2], nTrees = 5)
#' pred <- .predict.RandomForest_PC(model = RF[[1]], type='all',
#'                                  data_pred = pcaData[-2],
#'                                  data=pcaData[-2])
.predict.RandomForest_PC <- function(model, data_pred, type='all', data=NULL){
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
                #'ROC'=multiclass.roc(data_pred[,1], modelPercents)))
  }
}


#' Get Underlying Function
#'
#' This (internal) function is used to get the underlying functino (i.e. no _PC#
#'     in it). This is used to centralize called from RandomForest_PC.R and
#'     RandomForest_CVPC.R. If there is no PC (i.e. in meta-varaiables), the
#'     meta-variable is returned.
#'
#' @param names Vector of names to get underlying function on each
#'
#' @return Vector of underlying functions for each in name.
#' @export
#'
#' @examples
#' # See code for .computeRandomForest_PC. This is not an outward function so
#' #     won't be viewable.
.getUnderlyingVariable <- function(names, returnUnique=T){
  if(returnUnique)
    return(unique(stringr::str_remove(names,'(_PC)[0-9]+$')))
  return(stringr::str_remove(names,'(_PC)[0-9]+$'))
}
# substr(test, 1, max(1,regexpr("_PC[0-9]+$", test)-1))
