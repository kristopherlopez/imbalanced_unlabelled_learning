#############################################
### Positive unlabelled ensemble learning ###
#############################################

library(caret)
library(randomForest)
library(glmnet)
library(xgboost)
library(e1071)

####################################
##  Notes from ensemblePrediction ##
####################################

#Parameters
    #substrates: vector of labelled substrates
    #phospho.feature.data: the full dataset containing all substrates and features
    #ensemble.size: number of base classifiers in the model
    #size.negative: sample size of unlabelled instances used as 'negative cases' for training base classifiers
    #kernelType: underlying svm package from e1071 library. Choices: 'radial', 'polynomial', 'linear'.

#Underlying procedure
    # 1. Define set of positive labelled instances and unlabelled instances as negative
    # 2. Build an ensemble of r = ensemble.size base classifiers by
        # 3. Sampling without replacement from the unlabelled 'negative' pool, i.e. 'downweighting'
        # 4. Training an kernelised svm as a 'base classifier' on balanced training set, using 3-fold cross validation
        # 5. On the validation set for each base classifier, estimate the value of 'c', for adjustment as part of the
             #positive unlabelled ensemble learning procedure. This is achieved by assessing probabilities outputted
             #by the prediction model on known positive cases
    # 6. Now apply the ensemble of svms with paired correction factors 'c' to the full dataset. Then return:
        # $prediction: The average adjusted probability of 'positive class' across all ensembles, 
        #  normalised by the highest probability prediction
        # $svm_models: the unadjusted svm models for all base classifiers

#################################################################        
##  Apply the ensemblePrediction to the feature enhanced data  ##
#################################################################

#####################################
##  Data preparation and cleansing ##
#####################################


# Convert various features to numeric vectors
# ===========================================

#dat_modified$site <- as.numeric(dat_modified$site)
#dat_modified$`{sw_1=F,sw_2=R}` <- as.numeric(dat_modified$`{sw_1=F,sw_2=R}`)
#dat_modified$`{sw_1=F,sw_4=R}` <- as.numeric(dat_modified$`{sw_1=F,sw_4=R}`)
#dat_modified$`{sw_1=F,sw_7=S}` <- as.numeric(dat_modified$`{sw_1=F,sw_7=S}`)
#dat_modified$`{sw_2=R,sw_4=R}` <- as.numeric(dat_modified$`{sw_2=R,sw_4=R}`)
#dat_modified$`{sw_2=R,sw_7=S}` <- as.numeric(dat_modified$`{sw_2=R,sw_7=S}`)
#dat_modified$`{sw_4=R,sw_7=S}` <- as.numeric(dat_modified$`{sw_4=R,sw_7=S}`)
#dat_modified$`{sw_7=S,sw_8=P}` <- as.numeric(dat_modified$`{sw_7=S,sw_8=P}`)

# Get idx for positive and negative labels
# ========================================

get_negative_idx <- function(df, target_attr, attr_values) {
  return(rownames(df[!(df[[target_attr]] %in% attr_values), ]))
}

get_positive_idx <- function(df, target_attr, attr_values) {
  return(as.numeric(rownames(df[(df[[target_attr]] %in% attr_values), ])))
}

############################################################        
##  Build ensemble Predictors: SVM, GLM, RF, XGBoost, kNN  ##
############################################################

get_prediction <- function(
                    kinase, 
                    kinases, 
                    full_dataset, 
                    test_idx,
                    ensemble_size,
                    negative_size,
                    model_type,
                    svm_kernel,
                    svm_cost,
                    glm_family,
                    glm_alpha,
                    xgb_rounds,
                    xgb_depth,
                    rf_ntree,
                    rf_nodesize,
                    estimate_c = FALSE,
                    other_ksr = FALSE
                    ){
  
  positive_idx <- setdiff(kinases[[kinase]][['positive-idx']], test_idx)
  positive_train <- full_dataset[positive_idx, ]
  test_data <- full_dataset[test_idx, ]

  # Building the positive training set
  negative_unlabelled <- full_dataset[kinases[[kinase]][['unlabelled-idx']], ]
  
  # Sampling from the pool of 'negative' unlabelled
  pred_models <- list();
  model_estimated_c <- c();
  
  for (r in 1:ensemble_size){
    set.seed(r)
    
    idx <- sample(1:nrow(negative_unlabelled), size = negative_size, replace = F)
    negative_train <- negative_unlabelled[idx, ]

    # Creating the training set with positive and negative samples
    train_df <- rbind(positive_train, negative_train)
    rownames(train_df) <- NULL;
    
    # estimating c
    k <- 3
    cls <- as.factor(rep(c(1, 2), times=c(nrow(positive_train), negative_size)))
    fold <- createFolds(cls, k);
    
    # label 1 correspond to positive labelled examples
    p <- which(cls[fold$Fold1] == 1)
    
    if(model_type == "svm"){
      
      train_df <- train_df[, !(names(train_df) %in% features$nonnumeric)]
      train_df <- train_df[, !(names(train_df) %in% c('substrate_type'))]
      if(other_ksr == TRUE){ train_df <- train_df[, !(names(train_df) %in% kinase)]}
 
      pred_model <- svm(train_df, cls, kernel=svm_kernel, probability=TRUE, scale = FALSE)
      c_model <- svm(train_df[-fold$Fold1,], cls[-fold$Fold1], kernel=svm_kernel, probability=TRUE, scale = FALSE)
      c_pred <- predict(c_model, train_df[fold$Fold1,][p,], decision.values=F, probability=T);
      estimated_c <- sum(attr(c_pred, "probabilities")[,1]) / nrow(attr(c_pred, "probabilities"))
      
    } else if (model_type == 'glm'){
      
      if(ncol(train_df) > nrow(train_df)) { train_df = train_df[,1:nrow(train_df)-1]}
      train_df$substrate_type <- ifelse(train_df$substrate_type == "<unlabeled>", 0, 1)
      if(other_ksr == TRUE){ train_df <- train_df[, !(names(train_df) %in% kinase)]}
      
      train_x <- as.matrix(train_df[, !(names(train_df) %in% features$nonnumeric)])
      train_y <- as.matrix(train_df$substrate_type)
      
      cdf <- train_df[-fold$Fold1,]
      cx <- as.matrix(cdf[, !(colnames(train_df) %in% features$nonnumeric)])
      cy <- as.matrix(cdf$substrate_type)
      tx = as.matrix(cdf[fold$Fold1, !(colnames(train_df) %in% c('substrate_type'))][p,])
      
      pred_model <- suppressWarnings(cv.glmnet(x=train_x, y=train_y, type.measure='mse', family=glm_family, alpha=glm_alpha))
      c_model <- suppressWarnings(cv.glmnet(x=cx, y=cy, type.measure='mse', family=glm_family, alpha=glm_alpha))
      c_pred <- predict(c_model, s=c_model$lambda.1se, newx=tx, type='response');
      estimated_c <- sum(c_pred) / length(c_pred)
      
      
    } else if (model_type == 'rf'){
      
      if(other_ksr == TRUE){ 
        
        train_df <- train_df[, !(names(train_df) %in% 'substrate_type')]
        train_df[[kinase]] <- factor(train_df[[kinase]])

        pred_model <- randomForest(y=train_df[[kinase]], x=train_df[, !(names(train_df) %in% kinase)], norm.votes=TRUE, proximity=TRUE, ntree=rf_ntree, maxnodes=rf_nodesize, keep.forest = TRUE)

        c_model <- randomForest(y=train_df[-fold$Fold1,kinase], x=train_df[-fold$Fold1, !(names(train_df) %in% kinase)], norm.votes=TRUE, proximity=TRUE, ntree=rf_ntree, maxnodes=rf_nodesize, keep.forest = TRUE)
        c_pred <- predict(c_model, train_df[fold$Fold1,][p,], type='prob');
        
      } else {
        
        train_df$substrate_type <- factor(train_df$substrate_type)
        pred_model <- randomForest(y=train_df$substrate_type, x=train_df[, !(names(train_df) %in% c('substrate_type')) ], norm.votes=TRUE, proximity=TRUE, ntree=rf_ntree, maxnodes=rf_nodesize, keep.forest = TRUE)
        c_model <- randomForest(y=train_df[-fold$Fold1,'substrate_type'], x=train_df[-fold$Fold1, !(names(train_df) %in% c('substrate_type'))], norm.votes=TRUE, proximity=TRUE, ntree=rf_ntree, maxnodes=rf_nodesize, keep.forest = TRUE)
        c_pred <- predict(c_model, train_df[fold$Fold1,][p,], type='prob');
        
      }
      
      estimated_c <- sum(c_pred[,2]) / nrow(c_pred)
      
    } else if (model_type == 'xgb'){
      
      if(other_ksr == TRUE){ train_df <- train_df[, !(names(train_df) %in% kinase)]}
      train_df <- train_df[, !(names(train_df) %in% features$nonnumeric[-31])]
      
      param <- list(max.depth = xgb_depth, eta = 0.01,  objective="binary:logistic", subsample=0.9)
      train_df$substrate_type <- ifelse(train_df$substrate_type == "<unlabeled>", 0, 1)
      pred_model <- xgboost(param, label=data.matrix(train_df$substrate_type), data=data.matrix(train_df[, !(names(train_df) %in% c('substrate_type')) ]), objective='binary:logistic', nrounds=xgb_rounds, verbose=0)
      c_model <- xgboost(param, label=data.matrix(train_df[-fold$Fold1,'substrate_type']), data=data.matrix(train_df[-fold$Fold1, !(names(train_df) %in% c('substrate_type'))]), objective='binary:logistic', nrounds=xgb_rounds, verbose=0)
      c_pred <- predict(c_model, data.matrix(train_df[fold$Fold1,][p,]), type='prob');
      estimated_c <- sum(c_pred) / length(c_pred)
      
    } else if (model_type == 'knn'){
      
    } else if (model_type == 'nn'){
      
    }
    
    model_estimated_c <- c(model_estimated_c, estimated_c)
    
    # training base classifiers
    pred_models[[r]] <- pred_model
    
  }
    
  # an ensemble approach for prediction
  
  predict_df <- test_data
  predict_full <- full_dataset
  
  if (model_type == 'xgb') {
    
    predict_df <- predict_df[, !(names(predict_df) %in% features$nonnumeric[-31])]
    predict_df$substrate_type <- ifelse(predict_df$substrate_type == "<unlabeled>", 0, 1)
    
    predict_full <- predict_full[, !(names(predict_full) %in% features$nonnumeric[-31])]
    predict_full$substrate_type <- ifelse(predict_full$substrate_type == "<unlabeled>", 0, 1)
    
    
  } else if (model_type == 'rf') {
    
    predict_df$substrate_type = factor(predict_df$substrate_type)
    predict_full$substrate_type = factor(predict_full$substrate_type)
    
  } else if (model_type == 'glm') {
    
    predict_df <- predict_df[, !(names(predict_df) %in% features$nonnumeric)]
    predict_df <- as.matrix(predict_df)
    
    predict_full <- predict_full[, !(names(predict_full) %in% features$nonnumeric)]
    predict_full <- as.matrix(predict_full)

  }
  
  # correcting base classifiers with estimated c values
  pred_corrected <- 0
  pred_corrected_full <- 0
  
  for(i in 1:length(pred_models)) {
    if(model_type == 'svm'){
      
      predict_df <- predict_df[, !(names(predict_df) %in% features$nonnumeric)]
      predict_df <- predict_df[, !(names(predict_df) %in% c('substrate_type'))]
      if(other_ksr == TRUE){ predict_df <- predict_df[, !(names(predict_df) %in% kinase)]}

      pred <- predict(pred_models[[i]], predict_df, decision.values=T, probability=T)
      pred <- attr(pred,"probabilities")[,1]
      
      if(estimate_c == TRUE) { pred <- pred / model_estimated_c[i]}
      pred_corrected <- pred_corrected + pred
      
      predict_full <- predict_full[, !(names(predict_full) %in% features$nonnumeric)]
      predict_full <- predict_full[, !(names(predict_full) %in% c('substrate_type'))]
      if(other_ksr == TRUE){ predict_full <- predict_full[, !(names(predict_full) %in% kinase)]}
      
      full_pred <- predict(pred_models[[i]], predict_full, decision.values=T, probability=T)
      full_pred <- attr(full_pred,"probabilities")[,1]
      
      if(estimate_c == TRUE) { full_pred <- full_pred / model_estimated_c[i]}
      pred_corrected_full <- pred_corrected_full + full_pred
      
    } else if (model_type == 'glm') {
      
      if(other_ksr == TRUE){ predict_df <- predict_df[, !(names(predict_df) %in% kinase)]}

      pred <- predict(pred_models[[i]], s=pred_models[[i]]$lambda.1se, newx=predict_df, type='response');
      if(estimate_c == TRUE) { pred <- pred / model_estimated_c[i]}
      pred_corrected <- pred_corrected + pred
      
      if(other_ksr == TRUE){ predict_full <- predict_full[, !(names(predict_full) %in% kinase)]}
      
      full_pred <- predict(pred_models[[i]], s=pred_models[[i]]$lambda.1se, newx=predict_full, type='response');
      if(estimate_c == TRUE) { full_pred <- full_pred / model_estimated_c[i]}
      pred_corrected_full <- pred_corrected_full + full_pred
      
    } else if (model_type == 'rf') {
      
      if(other_ksr == TRUE){ predict_df <- predict_df[, !(names(predict_df) %in% kinase)]}
      
      pred <- predict(pred_models[[i]], predict_df, type='prob')[, 2]
      if(estimate_c == TRUE) { pred <- pred / model_estimated_c[i]}
      pred_corrected <- pred_corrected + pred
      
      if(other_ksr == TRUE){ predict_full <- predict_full[, !(names(predict_full) %in% kinase)]}
      
      full_pred <- predict(pred_models[[i]], predict_full, type='prob')[, 2]
      if(estimate_c == TRUE) { full_pred <- full_pred / model_estimated_c[i]}
      pred_corrected_full <- pred_corrected_full + full_pred

      
    } else if (model_type == 'xgb') {
      
      if(other_ksr == TRUE){ predict_df <- predict_df[, !(names(predict_df) %in% kinase)]}
      
      pred <- predict(pred_models[[i]], data.matrix(predict_df))
      if(estimate_c == TRUE) { pred <- pred / model_estimated_c[i]}
      pred_corrected <- pred_corrected + pred
      
      if(other_ksr == TRUE){ predict_full <- predict_full[, !(names(predict_full) %in% kinase)]}
      
      full_pred <- predict(pred_models[[i]], data.matrix(predict_full))
      if(estimate_c == TRUE) { full_pred <- full_pred / model_estimated_c[i]}
      pred_corrected_full <- pred_corrected_full + full_pred
      
    }
    
  }
  
  # return prediction results and base classifiers
  results <- list()
  predicts <- (pred_corrected / ensemble_size)
  predicts_full <- (pred_corrected_full / ensemble_size)
  results$prediction <- predicts
  results$prediction <- predicts / max(predicts)
  results$prediction_full <- predicts_full
  results$prediction_full <- predicts_full / max(predicts_full)
  results$pred_models <- pred_models;
  results$estimated_c <- estimated_c
  return(results);

}



