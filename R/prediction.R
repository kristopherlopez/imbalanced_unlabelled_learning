###########################
### Prediction function ###
###########################

create_prediction <- function(x_tr, x_ts, y_tr, y_ts, model, k=3, kernel='radial'){
  
  if(model == 'knn'){
    
    y_hat <- knn(x_tr, x_ts, y_tr, k=k)
    
  } else if(model == 'glm'){
    
    data_tr <- cbind(x_tr, y_tr)
    data_ts <- cbind(x_ts, y_ts)
    
    names(data_tr) <- c(names(x_tr), 'y_cls')
    
    y_model <- glm(y_cls~., family=binomial(link='logit'), data=data_tr, control = list(maxit = 50))
    y_prob  <- predict(y_model, newdata=data_ts, type='response')
    y_hat   <- ifelse(y_prob > 0.5, 1, 0)
    
  } else if(model == 'svm'){
    
    y_model <- svm(x_tr, y_tr, kernel=kernel, type="C-classification")
    y_hat  <- predict(y_model, newdata=x_ts, type='response')
    
  } else if(model == 'lda'){
    
    data_tr <- cbind(x_tr, y_tr)
    data_ts <- cbind(x_ts, y_ts)
    names(data_tr) <- c(names(x_tr), 'y_cls')
    
    y_model <- lda(y_cls~., data = data_tr)
    y_prob  <- predict(y_model, data_ts)$posterior[, "1"]
    y_hat   <- ifelse(y_prob > 0.5, 1, 0)
    
  } else if(model == 'lasso'){
    
    grid <- 10^seq(10,-2, length=100)
    
    y_model <- glmnet(as.matrix(x_tr), as.matrix(y_tr), alpha=1, lambda=grid)
    y_lambd <- cv.glmnet(as.matrix(x_tr), as.matrix(y_tr), alpha=1)$lambda.min
    y_prob   <- predict(y_model, s=y_lambd, newx=as.matrix(x_ts), type='response')
    y_hat   <- ifelse(y_prob > 0.5, 1, 0)
    
  } else if(model == 'ridge'){
    
    grid <- 10^seq(10,-2, length=100)
    
    y_model <- glmnet(as.matrix(x_tr), as.matrix(y_tr), alpha=0, lambda=grid)
    y_lambd <- cv.glmnet(as.matrix(x_tr), as.matrix(y_tr), alpha=0)$lambda.min
    y_prob  <- predict(y_model, s=y_lambd, newx=as.matrix(x_ts), type='response')
    y_hat   <- ifelse(y_prob > 0.5, 1, 0)
    
  } else if(model == 'rf'){
    
    y_tr    <- as.factor(y_tr)
    y_model <- randomForest(y=y_tr, x=x_tr, proximity=TRUE, ntree=100)
    y_hat   <- predict(y_model, x_ts);
    
  } else if(model == 'gbm'){
    
    data_tr <- cbind(x_tr, y_tr)
    data_ts <- cbind(x_ts, y_ts)
    
    names(data_tr) <- c(names(x_tr), 'y_cls')
    
    y_model <- gbm(y_cls~., data=data_tr, distribution="bernoulli", n.trees=5000)
    y_prob   <- predict(y_model, newdata=data_ts, n.trees = 1500, type="response")
    y_hat   <- ifelse(y_prob > 0.5, 1, 0)
    
  }
  
  return(y_hat)
  
} 


#######################################
### Models used in final prediction ###
#######################################

######################################
### Models used for Akt prediction ###
######################################

generate_akt_prob <- function(){ 

  kinase_count <- floor(length(kinases[['Akt']][['positive-idx']]))

  # Model 1: XGB
  model1 <- get_prediction(
    kinase = 'Akt',
    kinases = kinases,
    full_dataset = full_data[, features$temporal],
    test_idx = 1,
    ensemble_size = 50,
    negative_size = kinase_count * 1,
    model_type = 'xgb',
    xgb_rounds = 10,
    xgb_depth = 2,
    estimate_c = TRUE
  )$prediction_full
 
  # Model 2: Random Forest
  model2 <- get_prediction(
    kinase = 'Akt',
    kinases = kinases,
    full_dataset = full_data[, features$final],
    test_idx = 1,
    ensemble_size = 100,
    negative_size = kinase_count * 1,
    model_type = 'rf',
    rf_ntree = 500,
    rf_nodesize = 4,
    estimate_c = TRUE
  )$prediction_full

  # Model 3: SVM (Radial Kernel)
  model3 <- get_prediction(
    kinase = 'Akt',
    kinases = kinases,
    full_dataset = full_data[, features$relative],
    test_idx = 1,
    ensemble_size = 100,
    negative_size = kinase_count * 1,
    model_type = 'svm',
    svm_kernel =  'radial',
    svm_cost = 1,
    estimate_c = TRUE
  )$prediction_full

  # Model 4: Generalised Linear Model
  model4 <- get_prediction(
    kinase = 'Akt',
    kinases = kinases,
    full_dataset = full_data[, features$temporal],
    test_idx = 1,
    ensemble_size = 75,
    negative_size = kinase_count * 1.1,
    model_type = 'glm',
    glm_family = 'gaussian',
    glm_alpha = 1,
    estimate_c = FALSE
  )$prediction_full

  # Model 5: Random Forest (Motif classifier)
  model5 <- get_prediction(
    kinase = 'Akt',
    kinases = kinases,
    full_dataset = full_data[, features$sequence],
    test_idx = 1,
    ensemble_size = 25,
    negative_size = kinase_count * 1,
    model_type = 'rf',
    rf_ntree = 750,
    rf_nodesize = 4,
    estimate_c = TRUE
  )$prediction_full
  
  # Create ensemble of individual Akt models
  
  all_predictions <- cbind(model1, model2, model3, model4, model5)
  
  prediction_threshold <- 0.8
  
  # Three voting options: majority vote, all vote and average vote by prediction threshold
  
  ensemble_probability <- rowSums(all_predictions) / ncol(all_predictions)
  ensemble_yes_vote <- apply(all_predictions, 1, function(x) sum(x > prediction_threshold))
  ensemble_no_vote <- apply(all_predictions, 1, function(x) sum(x <= prediction_threshold))
  probability_decision <- ifelse(ensemble_probability > prediction_threshold, 1, 0)
  majority_decision <- ifelse(ensemble_yes_vote > ensemble_no_vote, 1, 0)
  all_decision <- ifelse(ensemble_yes_vote == ncol(all_predictions), 1, 0)
  
  akt_prob <- cbind(all_predictions, ensemble_yes_vote, ensemble_no_vote, ensemble_probability, probability_decision, majority_decision, all_decision)

  return(akt_prob)
}

#######################################
### Models used for mTOR prediction ###
#######################################

generate_mtor_prob <- function(){

  kinase_count <- floor(length(kinases[['mTOR']][['positive-idx']]))
  
  # Model 1: XGB
  model1 <- get_prediction(
    kinase = 'mTOR',
    kinases = kinases,
    full_dataset = full_data[, features$relative],
    test_idx = 1,
    ensemble_size = 50,
    negative_size = kinase_count * 1.2,
    model_type = 'xgb',
    xgb_rounds = 50,
    xgb_depth = 2,
    estimate_c = FALSE
  )$prediction_full
  
  # Model 2: Random Forest
  model2 <- get_prediction(
    kinase = 'mTOR',
    kinases = kinases,
    full_dataset = full_data[, features$final],
    test_idx = 1,
    ensemble_size = 25,
    negative_size = kinase_count * 1,
    model_type = 'rf',
    rf_ntree = 750,
    rf_nodesize = 4,
    estimate_c = TRUE
  )$prediction_full
  
  # Model 3: Support Vector Machine
  model3 <- get_prediction(
    kinase = 'mTOR',
    kinases = kinases,
    full_dataset = full_data[, features$relative],
    test_idx = 1,
    ensemble_size = 100,
    negative_size = kinase_count * 1,
    model_type = 'svm',
    svm_kernel =  'radial',
    svm_cost = 1,
    estimate_c = TRUE
  )$prediction_full
  
  # Model 4: Generalised Linear Model (Lasso)
  model4 <- get_prediction(
    kinase = 'mTOR',
    kinases = kinases,
    full_dataset = full_data[, features$temporal],
    test_idx = 1,
    ensemble_size = 100,
    negative_size = kinase_count * 1.1,
    model_type = 'glm',
    glm_family = 'binomial',
    glm_alpha = 1,
    estimate_c = FALSE
  )$prediction_full
  
  # Model 5: Random Forest (Motif classifier)
  model5 <- get_prediction(
    kinase = 'mTOR',
    kinases = kinases,
    full_dataset = full_data[, features$sequence],
    test_idx = 1,
    ensemble_size = 25,
    negative_size = kinase_count * 1,
    model_type = 'rf',
    rf_ntree = 500,
    rf_nodesize = 4,
    estimate_c = TRUE
  )$prediction_full
  
  # Create ensemble of individual mTOR models
  
  all_predictions <- cbind(model1, model2, model3, model4, model5)
  
  prediction_threshold <- 0.7
  
  # Three voting options: majority vote, all vote and average vote by prediction threshold
  
  ensemble_probability <- rowSums(all_predictions) / ncol(all_predictions)
  ensemble_yes_vote <- apply(all_predictions, 1, function(x) sum(x > prediction_threshold))
  ensemble_no_vote <- apply(all_predictions, 1, function(x) sum(x <= prediction_threshold))
  probability_decision <- ifelse(ensemble_probability > prediction_threshold, 1, 0)
  majority_decision <- ifelse(ensemble_yes_vote > ensemble_no_vote, 1, 0)
  all_decision <- ifelse(ensemble_yes_vote == ncol(all_predictions), 1, 0)
  
  mtor_prob <- cbind(all_predictions, ensemble_yes_vote, ensemble_no_vote, ensemble_probability, probability_decision, majority_decision, all_decision)

  return(mtor_prob)
  
}

# Create single probability file

generate_all_prob <- function(akt_prob, mtor_prob){
  
  all_prob <- cbind(full_data[, c('Identifier', 'Seq.Window', 'substrate_type')], akt_prob[,'ensemble_probability'], mtor_prob[,'ensemble_probability'])
  
  colnames(all_prob) <- c('Identifier', 'Seq.Window', 'substrate_type', 'akt$prediction_full', 'mtor$prediction_full')
  
  return(all_prob)
  
}

