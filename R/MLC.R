#' Maximum likelihood model class
#' 
#' @slot k the constant fraction to be used in model \eqn{\frac{1}{(2 \pi)^{\frac{L}{2}}  \sqrt{\left | \Sigma_i  \right |}}}
#' @slot mu mean (\eqn{\mu_i}) list for each variable and class
#' @slot invCovs inverted covariance matrix (\eqn{\Sigma_i}) for each class
#' @slot groups the classification levels
#' @slot vars the variables used for trainning the model
#' 
#' @export
setClass(
  "MLC.model",
    slots = list(
      k="numeric",
      mu="list",
      invCovs="list",
      groups="character",
      vars="character"
    )
)


#' MLC
#' 
#' Function to create classifier class with dataset
#' 
#' @param x feature vector
#' @param y response vector
#' 
#' @return An object of class \code{\link[tabularMLC]{MLC.model-class}} parameters used for the model  
#' 
#' @export
MLC = function(x, y, ridge=1e-9) {
  colNames = colnames(x)
  rowNames = rownames(x)
  if (class(x)[1] != "matrix") 
    x = as.matrix(x)
  if (class(y) != "factor")
    stop("Response vector is not a factor")
  
  xClasses = sapply(1:nrow(x), function(i) class(x[[i]]))
  if (!all(xClasses %in% c("numeric", "integer"))) {
    stop("Only continuous features are accepted!")
  }
  
  res = cpp_MLC(x, y, levels(y))
  model.class = new("MLC.model",
                    k = res[["k"]],
                    invCovs = res[["invCovs"]],
                    groups = res[["groups"]],
                    mu = res[["mu"]],
                    vars = colNames)
  
  return(model.class)
}


#' predict.MLC
#' 
#' Function to create classifier class with dataset
#' 
#' @param model \code{\link[tabularMLC]{MLC.model-class}} model class to use for prediction
#' @param x feature vector to predict
#' 
#' @return a factor vector with the predicted value
#' 
#' @export
predict.MLC = function(model, x, likelihood=FALSE) {
  nVars = length(model@vars)
  xVars = colnames(x)
  commonVars = intersect(model@vars, xVars)
  if (length(commonVars) != nVars) {
    stop(sprintf("It looks like the x data.frame is missing some variables found in model:
                 %s\n", paste0(setdiff(model@vars, xVars), collapse="\n")))
  }
  if (class(x)[1] != "matrix") {
    x = as.matrix(as.data.frame(x)[, model@vars])
  }
  
  likelihoods = cpp_predict(model, x)
  predict = factor(max.col(likelihoods), levels=1:length(model@groups))
  levels(predict) = model@groups
  if (likelihood==FALSE) 
    return(predict)
  
  colnames(likelihoods) = model@groups
  relative.like = likelihoods/rowSums(likelihoods)
  list(predict=predict, likelihoods=likelihoods, relative.like=relative.like)
}


#' Predict function for class MLC.model
#' 
#' @export
setMethod("predict", signature = c("MLC.model"), function(object, x, likelihood=FALSE) {
  predict.MLC(object, x, likelihood)
})
