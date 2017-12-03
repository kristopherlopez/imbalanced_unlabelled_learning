###########################
### Prediction function ###
###########################

create_prediction <- function(x_tr, x_ts, y_tr, y_ts, model, params){
  
  prediction <- list()
  y_hat   <- NULL
  y_model <- NULL

  numeric_columns <- sapply(x_tr, is.numeric)
  nx_tr <- x_tr[, numeric_columns]
  nx_ts <- x_ts[, numeric_columns]
    
  if(model == 'knn'){
    
   y_hat <- knn(nx_tr, nx_ts, y_tr, k=params[['knn_k']])

  } else if(model == 'glm'){
    
    data_tr <- cbind(nx_tr, y_tr)
    data_ts <- cbind(nx_ts, y_ts)
    
    names(data_tr) <- c(names(nx_tr), 'y_cls')
    
    y_model <- glm(y_cls~., family=binomial(link='logit'), data=data_tr, control = list(maxit = params[['glm_maxit']]))
    y_prob  <- predict(y_model, newdata=data_ts, type='response')
    y_hat   <- ifelse(y_prob > 0.5, 1, 0)
    
  } else if(model == 'svm'){
    
    y_model <- svm(nx_tr, y_tr, kernel=params[['svm_kernel']], type="C-classification")
    y_hat  <- predict(y_model, newdata=nx_ts, type='response')
    
  } else if(model == 'lda'){
    
    data_tr <- cbind(nx_tr, y_tr)
    data_ts <- cbind(nx_ts, y_ts)
    names(data_tr) <- c(names(nx_tr), 'y_cls')
    
    y_model <- lda(y_cls~., data = data_tr)
    y_prob  <- predict(y_model, data_ts)$posterior[, "1"]
    y_hat   <- ifelse(y_prob > 0.5, 1, 0)
    
  } else if(model == 'lasso'){
    
    grid <- 10^seq(10,-2, length=100)
    
    y_model <- glmnet(as.matrix(nx_tr), as.matrix(y_tr), alpha=1, lambda=grid)
    y_lambd <- cv.glmnet(as.matrix(nx_tr), as.matrix(y_tr), alpha=1)$lambda.min
    y_prob   <- predict(y_model, s=y_lambd, newx=as.matrix(nx_ts), type='response')
    y_hat   <- ifelse(y_prob > 0.5, 1, 0)
    
  } else if(model == 'ridge'){
    
    grid <- 10^seq(10,-2, length=100)
    
    y_model <- glmnet(as.matrix(nx_tr), as.matrix(y_tr), alpha=0, lambda=grid)
    y_lambd <- cv.glmnet(as.matrix(nx_tr), as.matrix(y_tr), alpha=0)$lambda.min
    y_prob  <- predict(y_model, s=y_lambd, newx=as.matrix(nx_ts), type='response')
    y_hat   <- ifelse(y_prob > 0.5, 1, 0)
    
  } else if(model == 'rf'){
    
    y_tr    <- as.factor(y_tr)
    y_model <- randomForest(y=y_tr, x=x_tr, 
                            norm.votes=TRUE,
                            proximity=TRUE, 
                            ntree=params[['rf_ntree']],
                            maxnodes=params[['rf_nodesize']],
                            keep.forest = TRUE
                            )
    y_hat   <- predict(y_model, x_ts);
    
  } else if(model == 'gbm'){
    
    data_tr <- cbind(nx_tr, y_tr)
    data_ts <- cbind(nx_ts, y_ts)
    
    names(data_tr) <- c(names(nx_tr), 'y_cls')
    
    y_model <- gbm(y_cls~., 
                  data=data_tr, 
                  distribution=params[['gbm_distribution']], 
                  n.trees=params[['gbm_ntrees']]
                  )
    y_prob  <- predict(
                  y_model, 
                  newdata=data_ts, 
                  n.trees = params[['gbm_ttrees']],
                  type="response"
                  )
    y_hat   <- ifelse(y_prob > 0.5, 1, 0)
    
  } else if(model == 'xgb'){
    
    
    
  }
  
  prediction[['y_hat']]   <- y_hat
  prediction[['y_model']] <- y_model
  
  return(prediction)
  
}