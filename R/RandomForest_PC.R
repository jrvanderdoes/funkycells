
#' Compute Random Forest for Data with multiple PC (Along With Meta-Variables)
#' 
#' This (internal) function creates a random forest model for data with PCs and 
#'     meta-variables. This includes proper combination for variable importance. 
#'     This is internal because we want users to use randomForest_CVPC.
#'
#' Upcoming: See TEST comment
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
#' @return A list with two entries
#'     1. model: List of CART that builds the andom forest model 
#'     2. varImportanceData: Data.frame for variable importance information
#' @export
#' 
#' @examples
#' # See code for randomForest_CVPC. This is not an outward function so won't be 
#' #     viewable.
.computeRandomForest_PC  <- function(data, outcome=colnames(data)[1],
                              nTrees = 1000, varImpPlot=T, metaNames=NULL){
  # Ensure this is worthwhile
  if(length(unique(data[,outcome]))==1)
    stop('Error: Only 1 outcome in data, cannot do RF')
  
  # Setup
  # columns: all columns in data
  # cols_preds: all predictors name
  columns <- colnames(data)
  for(i in 2:length(columns)){
    columns[i] <- substr(columns[i],1,
                         tail(unlist(gregexpr('_', columns[i])), n=1)-1)
  }
  columns <- columns[columns!=""] # Drop all blanks, which will be meta-vars
  columns <- c(columns,metaNames)
  
  cols_preds <- unique(columns!=outcome) # Get overall vars names
  
  data_result <- data.frame('var'=cols_preds,
                            'splits'=0,
                            'giniDec'=0,
                            'varImp'=0,
                            'varImpCt'=0,
                            'vi'=0)
  RF <- list()
  # To do CART
  for(i in 1:nTrees){
    # Bootstrap observations
    data_rf <- data[sample(1:nrow(data),nrow(data), replace = T),]
    # Ensure there are multiple responses in bootstrapped sample.
    while(length(unique(data_rf[,outcome]))==1){
      data_rf <- data[sample(1:nrow(data),nrow(data), replace = T),]
    }
    
    # Variable selection. Currently 80%
    preds <- unique(sample(cols_preds, 0.8*length(cols_preds), replace = F))
    # Build RF data. Note estimators are scrambled to avoid computational bias
    data_rf <- data_rf[columns %in% c(outcome, preds)]
    data_rf <- cbind(data_rf[outcome],data_rf[,sample(2:ncol(data_rf))])
    
    # Fit CART
    model <- rpart(paste0(outcome,' ~ .'), data=data_rf, method="class",
                   control=rpart.control(minsplit =1,minbucket=1, cp=0))
    data_result <- .computeTotalVarImportance(model, data_result)
    
    RF[[i]] <- model # Save model 
  }
  
  # Get mean results for variable importance
  #     TEST: Use average over all or just when appears in model
  data_result$avgGini <- data_result$giniDec / nTrees#data_result$splits
  data_result$avgVI <- data_result$varImp / nTrees#data_result$varImpCt
  

  # Get Variable importance
  if(varImpPlot){
    return(list(RF, data_result, 
                .plotVariableImportance(varImportanceData=data_result)))
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
  frame <- model$frame
  frame[['gini']] = 1 - (frame[['dev']] / frame[['n']])^2 - 
    (1 - frame[['dev']] / frame[['n']])^2
  #frame[,c('var','n','dev','gini')]
  
  # Get Decrease for each split
  frame[['improve']] = NA
  for (j in 1:nrow(frame)) {
    if (frame[j,'var'] == '<leaf>') next
    
    ind = which(rownames(frame) %in% 
                  (as.numeric(rownames(frame)[j])*2+c(0,1)))
    frame[j,'improve'] = frame[j,'n']*frame[j,'gini'] - frame[ind[1],'n']*
      frame[ind[1],'gini'] - frame[ind[2],'n']*frame[ind[2],'gini']
    
    # Get Gini improvement
    name <- frame[j,'var']
    
    # Get name of bigger group
    if(sum(gregexpr('_', name)[[1]]!=-1)>0){
      # If cell interaction
      name_sub <- substring(name,1,tail(unlist(gregexpr('_', name)), n=1)-1)
    }else{
      # If meta data
      name_sub <- name
    }
    
    totVarImportance[totVarImportance$var==name_sub,'giniDec'] <- 
      totVarImportance[totVarImportance$var==name_sub,'giniDec'] +
      frame[j,'improve']
    # Count split by variable
    totVarImportance[totVarImportance$var==name_sub,'splits'] <- 
      totVarImportance[totVarImportance$var==name_sub,'splits'] + 1
  }
  
  # Get and count Var Import
  for(j in 1:length(model$variable.importance)){
    name <- names(model$variable.importance[j])
    
    # Get name of bigger group
    if(sum(gregexpr('_', name)[[1]]!=-1)>0){
      # If cell interaction
      name_sub <- substring(name,1,tail(unlist(gregexpr('_', name)), n=1)-1)
    }else{
      # If meta data
      name_sub <- name
    }
    
    totVarImportance[totVarImportance$var==name_sub,'varImp'] <-
      totVarImportance[totVarImportance$var==name_sub,'varImp'] +
      model$variable.importance[[j]]
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
  if((is.null(varImportanceData) && !is.null(model))||
     (!is.null(varImportanceData) && is.null(model))){
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
      varImportanceData <- .computeTotalVarImportance(model, varImportanceData)
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
  
  avgGini <- ggplot() + 
    geom_point(aes(x=reorder(var,avgGini), y=avgGini/max(avgGini)), 
               data=varImportanceData) + 
    coord_flip() +
    xlab(NULL) +
    ylim(c(0,1)) +
    ylab('Gini') +
    theme_bw()
  avgVI <- ggplot() + 
    geom_point(aes(x=reorder(var,avgVI), y=avgVI/max(avgVI)), 
               data=varImportanceData) + 
    coord_flip() +
    xlab(NULL) +
    ylim(c(0,1)) +
    ylab('VarImp') +
    theme_bw()
  
  varImportance <- grid.arrange(avgGini,avgVI,
                                layout_matrix = rbind(c(1,2)),
                                bottom = 'Variable Importance - Percent of Max')
  
  varImportance
}

#' Predict Using RandomForest_PC
#' 
#' This (internal) function gets the predicted value from a RandomForest_PC 
#'     model.
#' 
#' Upcoming: Solve TODO and document this better
#'          Test this whole thing
#'
#' @param model RandomForest_PC model. See .computeRandomForest_PC. A list of 
#'     CART models from rpart.
#' @param data_pred 
#' @param type (Optional) String indicating type of analysis. Options are pred, 
#'     all or ROC. The choice changes the return to best fit intended use.
#' @param data (Optional) 
#'
#' @return The returned data depends on type/
#'     1. type='pred': returns a (??Vector??) of the predictions
#' @export
#'
#' @examples
.predict.RandomForest_PC <- function(model, data_pred, type='all', data=NULL){
  # Verification - If data 
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
  predictions <- data.frame(matrix(nrow=nrow(data_pred),
                                   ncol=length(model)))
  data_pred_use <- data_pred
  drop_pred_rows <- rep(0, nrow(data_pred))
  
  # Get prediction for each tree
  for(i in 1:length(model)){
    # If original data given, this will allow verification
    if(!is.null(data)){
      # Reset
      data_pred_use <- data_pred
      
      # Get rows of data used
      data_row <- names(unlist(model[[i]]$where))
      unique_rows <- unique(gsub("\\..*","", data_row ))
      
      # Check each variable to see if contained in fit data
      #   Only check chars and factors (nums won't have this problem)
      drop_spot <- data.frame(matrix(0, 
                                     ncol=length(checkcols),
                                     nrow=nrow(data_pred_use))) # 0 keep, 1 drop
      for(j in checkcols){ 
        # If missing data, see if for single entry or all
        if(sum(!(data_pred_use[,j] %in% data[unique_rows,j]))>0){
          for(k in 1:nrow(data_pred_use)){
            if(!(data_pred_use[k,j] %in% data[unique_rows,j]))
              data_pred_use[k,j] <- NA
          }
        }
      }
      # If all pred_data has issues, skip iteration
      if(nrow(data_pred_use)==0)
        next
    }
    
    ## Predict the outcome using the model
    pred_mat <- predict(model[[i]], data_pred_use)
    
    # See if anything was dropped
    if(nrow(pred_mat)!= length(drop_pred_rows)){
      stop('Error: Contact us. This should not happen.')
    }
    
    sols <- colnames(pred_mat)
    for(j in 1:nrow(data_pred)){
      max_pred <- which.max(pred_mat[j,])
      if(length(max_pred)!=0)
        predictions[j,i] <- sols[max_pred]
    }
  }
  
  # Combine predictions
  if(is.numeric(predictions)){
    stop('Error: Sorry still need to figure out numeric')
  }else{
    # Vote (Randomly selects on ties)
    modelPredictions <- apply(predictions, MARGIN = 1, 
                              function(x){ tmp<-getMode(x)
                              tmp[sample(1:length(tmp),1)]}) 
    # see percent of each (TODO: Review this)
    modelPercents <- as.data.frame(t(apply(predictions, MARGIN = 1, 
                                           function(x,options){
                                             sapply(options, function(option,x){
                                               sum(x==option)/length(x)},
                                               x=x)
                                           }, options=unique(data_pred[,1]))))
  }
  
  # Return Results
  if(type=='pred'){
    # This returns only the predictions
    return(modelPredictions) 
  }else if(type=='all'){
    
    return(list('PredPerc'=cbind(modelPredictions,modelPercents),
                'Acc'=sum(modelPredictions==data_pred[,1])/nrow(data_pred),
                'ROC'=multiclass.roc(data_pred[,1], modelPercents)))
  }else if(type=='ROC'){
    ## TODO:: Look into this
    modelPredictions <- as.data.frame(t(apply(predictions, MARGIN = 1, 
                                              function(x,options){
                                                sapply(options, function(option,x){
                                                  sum(x==option)/length(x)},
                                                  x=x)
                                              }, options=unique(data_pred[,1]))))
    
    return(list('Percents'=modelPredictions,
                'ROC'=multiclass.roc(data_pred[,1], modelPredictions),
                'Acc'=sum(modelPredictions[,1]==data_pred[,1])/nrow(data_pred)))
  }
}


#' Title
#' 
#' This may be removed entirely. I will play around and see first
#'
#' @param modelPredictions 
#' @param modelPercents 
#'
#' @return
#' @export
#'
#' @examples
.computePseudoROCCurve <- function(modelPredictions, modelPercents){
  ## My implementation of a fake ROC curve 
  
  options <- unique(modelPredictions)
  cols <- brewer.pal(length(options),'Set1')
  # 1 ROC curve, mock vs non mock
  roc.curve <- pROC::roc(ifelse(modelPredictions==options[1], 
                                options[1], paste0('Not-',options[1])),
                         as.numeric(modelPercents[,options[1]]),
                         levels=c(options[1], paste0('Not-',options[1])),
                         direction='>'
  )
  plot(roc.curve, col = cols[1])
  
  for(i in 2:length(options)){
    roc.curve <- pROC::roc(ifelse(modelPredictions==options[i], 
                                  options[i], paste0('Not-',options[i])), 
                           as.numeric(modelPercents[,options[i]]),
                           levels=c(options[i], paste0('Not-',options[i])),
                           direction='>')
    lines(roc.curve, col = cols[i])
  }
}