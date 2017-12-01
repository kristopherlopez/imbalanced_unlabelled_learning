cross_validation <- function(x, y, folds, model){
  
  results_predict <- init_results()
  fold <- createFolds(y, folds)
  
  for(i in 1:folds){
    
    x_train <- as.data.frame(x[-fold[[i]], ])
    x_test  <- as.data.frame(x[fold[[i]], ])
    
    y_train <- y[-fold[[i]]]
    y_test  <- y[fold[[i]]]
    
    y_hat <- create_predict(x_train, x_test, y_train, y_test, model)
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