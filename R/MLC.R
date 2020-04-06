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
                    k=res[[1]],
                    invCovs=res[[2]],
                    groups=res[[3]])
  
  return(model.class)
}
