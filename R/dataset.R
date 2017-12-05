



# Parameters: class types, feature sets and model options

params <- list()
params[['substrate_type']] <- c('Akt', 'mTOR')
params[['feature_sets']] <- c('temporal', 'relative', 'retained', 'final')
params[['ensemble_size']] <- seq(0, 100, 25)[-1]
params[['negative_size']] <- c(1, 1.1, 1.2)
params[['bmm_threshold']] <- seq(0, 1, 0.01)
params[['model_type']] <- c('svm', 'glm', 'xgb', 'rf')
params[['knn_k']] <- 5
params[['svm_kernel']] <- c('radial')
params[['svm_cost']] <- c(0, 0.5, 1, 2)
params[['glm_family']] <- c('binomial', 'gaussian')
params[['glm_alpha']] <- c(0, 1)
params[['glm_lambda']] <- 10^seq(10, -2, length = 100)
params[['glm_maxit']] <- 50
params[['gbm_ntrees']] <- c(2500, 3500, 5000, 7500, 10000)
params[['gbm_ttrees']] <- 1500
params[['gbm_distribution']] <- c('binomial', 'adaboost')
params[['xgb_nrounds']] <- c(10, 15, 25, 50)
params[['xgb_depth']] <- c(2, 3, 4)
params[['rf_ntree']] <- c(10, 20, 50, 100)
params[['rf_nodesize']] <- c(4, 6, 8)
params[['estimate_c']] <- c(TRUE, FALSE)
