#' Tabular maximum likelihood classifier
#'
#' Maximum likelihood is a common classifier used for land use classification.
#' It calculates the likelihood of an object to belong to each class based on
#' an expected distribution and a metric of distance. 
#' 
#' The most common implementation, like in this package, will assume normal
#' distributed variables within classes, and calculate the distance, based
#' on Mahalanobis distance.
#' 
#' Imports
#' @useDynLib tabularMLC, .registration = TRUE
#' @import Rcpp
#' @importFrom Rcpp sourceCpp
"_PACKAGE"