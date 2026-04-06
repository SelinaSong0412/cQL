#' cQL: Clustered Q-Learning for cSMART Data
#'
#' `cQL` implements the manuscript's clustered Q-learning workflow for
#' two-stage clustered SMART data using the M-out-of-N cluster bootstrap.
#' The main entry point is [cQL()], which fits the stage-2 and stage-1
#' Q-functions on one user-supplied cSMART dataset and returns bootstrap-based
#' confidence intervals and p-values for the regression coefficients.
#'
#' @docType package
#' @name cQL
#' @importFrom geepack geeglm
#' @importFrom sandwich vcovCL
#' @importFrom stats coef formula lm model.matrix plogis predict quantile qt sd terms
"_PACKAGE"
