# Forward stepwise

forwardSelect <- function(d_x, d_y, model, validation) {
  
  best_features <- c()
  current_best_accuracy <- -Inf
  
  for(i in 1:ncol(d_x)) {
    
    selected_f <- NULL
    
    for(i in 1:ncol(d_x)) {
      
      current_f <- colnames(d_x)[i]
      
      if(!current_f %in% best_features) {
        
        df_x <- d_x[, names(d_x) %in% c(best_features, current_f)]
        df_x <- as.data.frame(df_x)
        
        if(validation == 'cv'){ v_res <- cross_validation(df_x, d_y, 10, model) }
        if(validation == 'ho'){ v_res <- hold_out(df_x, d_y, 0.7, model) } 
        
        if(v_res[['Ac']] > current_best_accuracy) {
          
          current_best_accuracy <- v_res[['Ac']]
          selected_f <- colnames(d_x)[i]
        }
      }
    }
    
    if(is.null(selected_f)) { break } else { best_features <- c(best_features, selected_f) }
    
  }
  
  return(best_features)
}

# Backward stepwise

backwardSelect <- function(d_x, d_y, model, validation) {
  
  worst_features <- c()
  current_best_accuracy <- -Inf
  
  if(validation == 'cv'){ current_best_accuracy <- cross_validation(d_x, d_y, 10, model)[['Ac']] }
  if(validation == 'ho'){ current_best_accuracy <- hold_out(d_x, d_y, 0.7, model)[['Ac']] } 
  
  for(i in 1:ncol(d_x)) {
    
    selected_f <- NULL
    
    for(i in 1:ncol(d_x)) {
      
      current_f <- colnames(d_x)[i]
      
      if(!current_f %in% worst_features) {
        
        df_x <- d_x[, setdiff(names(d_x), c(worst_features, current_f))]
        df_x <- as.data.frame(df_x)
        
        if(validation == 'cv'){ v_res <- cross_validation(df_x, d_y, 10, model) }
        if(validation == 'ho'){ v_res <- hold_out(df_x, d_y, 0.7, model) } 
        
        if(v_res[['Ac']] > current_best_accuracy) {
          current_best_accuracy <- v_res[['Ac']]
          selected_f <- colnames(d_x)[i]
        }
      }
    }
    
    if(is.null(selected_f)) { break } else { worst_features <- c(worst_features, selected_f) }
    
  }
  
  return(setdiff(names(d_x), worst_features))
}

# Best subset

bestSubset <- function(d_x, d_y, model, validation, max_f=2) {
  
  best_subset <- NULL
  current_best_accuracy <- -Inf
  
  for(i in 1:ncol(d_x)) {
    
    if(i <= max_f){
      
      combinations <- combn(names(d_x), i)
      
      for(j in 1:ncol(combinations)){
        
        current_f <- combinations[, j]
        df_x <- d_x[, names(d_x) %in% current_f]
        df_x <- as.data.frame(df_x)
        
        if(validation == 'cv'){ v_res <- cross_validation(df_x, d_y, 10, model) }
        if(validation == 'ho'){ v_res <- hold_out(df_x, d_y, 0.7, model) } 
        
        if(v_res[['Ac']] > current_best_accuracy) {
          current_best_accuracy <- v_res[['Ac']]
          best_subset <- current_f
        }
      }
    }
  }
  
  return(best_subset)
}