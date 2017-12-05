
downsampled_indices <- function(y, majority_class, balance=1){
  
  y <- as.data.frame(y)
  
  minority_size <- length(which(y != majority_class))
  minority_ix   <- which(y != majority_class)
  majority_ix   <- which(y == majority_class)
  
  y_pool <- y[majority_ix, ]

  downsampled_ix <- sample(1:length(y_pool), size = minority_size * balance, replace = F)
  
  new_ix <- c(minority_ix, downsampled_ix)

  return(new_ix)
}

loocv_unlabelled_preparation <- function(x, y, unlabelled_class){
  
  l_ix  <- which(y != unlablelled_class)
  ul_ix <- which(y == unlabelled_class)
  
  loocv_ix <- sample(1:length(l_ix), size=1)
  loocv_df <- 
  
}

cross_validation <- function(x, y, folds, model, params, labelled_ix=NULL){
  
  results_predict <- init_results()
  fold <- createFolds(y, folds)
  
  for(i in 1:folds){
    
    x_train <- as.data.frame(x[-fold[[i]], ])
    x_test  <- as.data.frame(x[fold[[i]], ])
    
    y_train <- y[-fold[[i]]]
    y_test  <- y[fold[[i]]]

    y_hat <- create_prediction(x_train, x_test, y_train, y_test, model, params)[['y_hat']]
    results_predict <- evaluate_predict(results_predict, y_test, y_hat)
  }
  
  results <- evaluate_results(
    results_predict[['TN']], 
    results_predict[['FP']],
    results_predict[['TP']], 
    results_predict[['FN']]
  )
  
  return(results)
  
}

hold_out <- function(x, y, part, model){
  
  partition <- createDataPartition(y, p = part)[[1]]
  
  x_train <- x[partition, ]
  x_test  <- x[-partition, ]
  
  y_train <- y[partition]
  y_test  <- y[-partition]
  
  y_hat <- create_predict(x_train, x_test, y_train, y_test, model)
  
  results_predict <- init_results()
  results_predict <- evaluate_predict(results_predict, y_test, y_hat)
  
  results <- evaluate_results(
    results_predict[['TN']], 
    results_predict[['FP']],
    results_predict[['TP']], 
    results_predict[['FN']]
  )
  
  return(results)
}