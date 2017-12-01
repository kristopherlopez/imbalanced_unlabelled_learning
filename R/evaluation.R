init_results <- function(){
  
  r <- list()
  r[['TP']] <- c()
  r[['TN']] <- c()
  r[['FP']] <- c()
  r[['FN']] <- c()
  
  return(r)
}

# populate list with results
evaluate_predict <- function(results, truth, preds){
  
  results$TP <- c(results$TP, sum((truth == preds)[truth == 1]))
  results$TN <- c(results$TN, sum((truth == preds)[truth == 0]))
  results$FP <- c(results$FP, sum((truth != preds)[truth == 0]))
  results$FN <- c(results$FN, sum((truth != preds)[truth == 1]))
  
  return(results)
}

# accuracy
Ac <- function(mat){
  apply(mat, 2, function(x) {
    TN <- x[1]
    FP <- x[2]
    TP <- x[3]
    FN <- x[4]
    (TP+TN)/(TP+TN+FP+FN)
  })  
}

# sensitivity
Se <- function(mat) {
  apply(mat, 2, function(x) {
    TN <- x[1]
    FP <- x[2]
    TP <- x[3]
    FN <- x[4]
    TP/(TP+FN)
  })
}

# specificity
Sp <- function(mat) {
  apply(mat, 2, function(x) {
    TN <- x[1]
    FP <- x[2]
    TP <- x[3]
    FN <- x[4]
    TN/(FP+TN)
  })
}

# F1 score
F1 <- function(mat) {
  apply(mat, 2, function(x){
    TN <- x[1]
    FP <- x[2]
    TP <- x[3]
    FN <- x[4]
    2*TP/(2*TP+FP+FN)
  })
}

# geometric mean
GMean <- function(mat) {
  apply(mat, 2, function(x){
    TN <- x[1]
    FP <- x[2]
    TP <- x[3]
    FN <- x[4]
    sqrt((TP/(TP+FN))*(TP/(TP+FP)))
  })
}   

evaluate_results <- function(TN, FP, TP, FN) {
  mat <- rbind(TN, FP, TP, FN)
  new_results <- list()
  new_results[['Ac']] <- round(mean(Ac(mat)), digits=3)
  new_results[['Se']] <- round(mean(Se(mat)), digits=3)
  new_results[['Sp']] <- round(mean(Sp(mat)), digits=3)
  new_results[['F1']] <- round(mean(F1(mat)), digits=3)
  new_results[['GM']] <- round(mean(GMean(mat)), digits=3)
  
  return(new_results)
}

# wrapper function for evaluation
evaluate_print <- function(TN, FP, TP, FN) {
  mat <- rbind(TN, FP, TP, FN)
  
  cat(c("Accuracy:", round(mean(Ac(mat)), digits=3)))
  cat(" ")
  
  cat(c("Sensitivity:", round(mean(Se(mat)), digits=3)))
  cat(" ")
  
  cat(c("Specificity:", round(mean(Sp(mat)), digits=3)))
  cat(" ")
  
  cat(c("F1 score:", round(mean(F1(mat)), digits=3)))
  cat(" ")
  
  cat(c("Geometric Mean:", round(mean(GMean(mat)), digits=3)))
  cat(" ")
  
}