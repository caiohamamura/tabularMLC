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
#' @references 
#' Mather, P. M. (1985). Remote sensing letters: A computationally efficient 
#' maximum-likelihood classifier employing prior probabilities for remotely-sensed 
#' data. International Journal of Remote Sensing, 6(2), 369–376. 
#' \doi{10.1080/01431168508948456}
#' 
#' Imports
#' @useDynLib tabularMLC, .registration = TRUE
#' @import Rcpp
#' @importFrom Rcpp sourceCpp
"_PACKAGE"