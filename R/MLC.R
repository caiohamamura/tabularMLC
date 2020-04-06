#' Maximum likelihood model class
#' 
#' @slot 
#' 
#' @export
setClass(
  "MLC.model",
    slots = list(
      k="numeric",
      invCovs="list",
      groups="character"
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
MLC = function(x, y) {
  colNames = colnames(x)
  rowNames = rownames(x)
  if (class(x) != "matrix") 
    x = as.matrix(x)
  if (class(y) != "factor")
    stop("Response vector is not a factor")
  
  res = cpp_MLC(x, y, levels(y))
  model.class = new("MLC.model",
                    k = res[["k"]],
                    invCovs = res[["invCovs"]],
                    groups = res[["groups"]],
                    mu = res[["mu"]],
                    vars = colNames)
  
  return(model.class)
}
