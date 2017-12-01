#####################################
### Common functions and settings ###
#####################################

###############
##  Settings ##
###############

# Paths
# =====

paths <- list()
paths[["root"]] <- dirname(rstudioapi::getSourceEditorContext()$path)
paths[["data"]] <- file.path(paths[["root"]], "data")
paths[["output"]] <- file.path(paths[["root"]], "output")