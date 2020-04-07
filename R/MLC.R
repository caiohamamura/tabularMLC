#' Maximum likelihood model class
#' 
#' @slot k the constant fraction to be used in model \eqn{\frac{1}{(2 \pi)^{\frac{L}{2}}  \sqrt{\left | \Sigma_i  \right |}}}
#' @slot mu mean (\eqn{\mu_i}) list for each variable and class
#' @slot inverseCovarianceMatrices inverted covariance matrix (\eqn{\Sigma_i}) for each class
#' @slot groups the classification levels
#' @slot vars the variables used for trainning the model
#' 
#' @seealso \code{\link[tabularMLC:MLC]{MLC}} which creates this class
#' 
#' @export
setClass(
  "MLC.model",
    slots = list(
      k="numeric",
      mu="list",
      inverseCovarianceMatrices="list",
      groups="character",
      vars="character"
    )
)


#' Maximum Likelihood Classifier
#' 
#' Function to create the classifier class from the training set
#' 
#' @aliases MLC.default MLC.formula
#' 
#' @param x feature vector for the training set
#' @param y factor vector with the training set labels
#' @param ... for other signatures
#' 
#' @return An object of class \code{\link[tabularMLC:MLC.model-class]{MLC.model}} parameters used for the model  
#' 
#' @examples 
#' data(iris)
#' 
#' x = iris[, -5]
#' y = iris$Species
#' 
#' # Default x y interface
#' mlcModel1 = MLC(x, y)
#' 
#' # Formula interface
#' mlcModel2 = MLC(Species ~ Petal.Length + Petal.Width, iris)
#' 
#' # Formula except one column
#' mlcModel3 = MLC(Species ~ . - Sepal.Length, iris)
#' 
#' @export
MLC = function(x, ...) {
  UseMethod("MLC")
}

#' 
#' @param formula \code{formula}. The formula for defining the model.
#' @param data the dataset
#' 
#' @rdname MLC
#' @export
MLC.formula = function(formula, data=NULL, ...) {
  if (is.null(data)) {
    stop("You need to provide a data argument!")
  }
  results = parseFormula(formula, data)
  MLC.default(results[[1]], results[[2]])
}

#' @rdname MLC
#' @importFrom methods new
#' @export
MLC.default = function(x, y=NULL, ...) {
  if (is.null(y)) {
    stop("You need to provide a y argument!")
  }
  colNames = colnames(x)
  rowNames = rownames(x)
  if (class(x)[1] != "matrix") 
    x = as.matrix(x)
  if (class(y) != "factor")
    if (class(y[,1]) != "factor") {
      stop("Response vector is not a factor")
    } else {
      y = y[,1]
    }
  
  xClasses = sapply(1:nrow(x), function(i) class(x[[i]]))
  if (!all(xClasses %in% c("numeric", "integer"))) {
    stop("Only continuous features are accepted!")
  }
  
  res = cpp_MLC(x, y, levels(y))
  model.class = new("MLC.model",
                    k = res[["k"]],
                    inverseCovarianceMatrices = res[["inverseCovarianceMatrices"]],
                    groups = res[["groups"]],
                    mu = res[["mu"]],
                    vars = colNames)
  
  return(model.class)
}

#' Predict function for MLC.model-class
#' 
#' \code{predict} is inherited from the generic function for predictions from the results. 
#' 
#' @param object \code{\link[tabularMLC]{MLC.model-class}} model class to use for prediction
#' @param x data.frame. The feature vector to predict
#' @param likelihood logical. Whether to return or not the likelihood values, default FALSE.
#' @param ... inherited from generic function (not in use)
#' 
#' @return a factor vector with the predicted value. If
#' likelihood is TRUE, then it will also return the calculated likelihoods.
#' 
#' @examples 
#' data(iris)
#' 
#' n = length(iris$Species)
#' 
#' # Split training by sample
#' training = sample(1:n, size=n*0.7)
#' validation = (1:n)[-training]
#' 
#' # Train model with training dataset
#' mlcModel = MLC(Species ~ ., iris[training,])
#' 
#' # Predict using validation dataset
#' predict = predict(mlcModel, iris[validation,])
#' 
#' # Print confusion matrix
#' confusionMatrix = table(predicted=predict, observed=iris$Species[validation])
#' print(confusionMatrix)
#'
#' @rdname predict
#' @method predict MLC.model
#' @export
predict.MLC.model = function(object, x=NULL, likelihood=FALSE, ...) {
  if (is.null(x)) {
    stop("You must provide an x vector to predict!")
  }
  model = object
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

#' @import stats
parseFormula = function(formula, data) {
  ### code gratefully stolen from randomForest.formula (package randomForest).
  ###
  if (!inherits(formula, "formula"))
    stop("method is only for formula objects")
  m <- match.call(expand.dots = FALSE)
  ## Catch xtest and ytest in arguments.
  if (any(c("xtest", "ytest") %in% names(m)))
    stop("xtest/ytest not supported through the formula interface")
  names(m)[2] <- "formula"
  if (is.matrix(eval(m$data, parent.frame())))
    m$data <- as.data.frame(data)
  m$... <- NULL
  m[[1]] <- as.name("model.frame")
  m <- eval(m, parent.frame())
  
  y <- stats::model.response(m)
  Terms <- attr(m, "terms")
  attr(Terms, "intercept") <- 0
  ## Drop any "negative" terms in the formula.
  m <- stats::model.frame(stats::terms(stats::reformulate(attributes(Terms)$term.labels)),
                   data.frame(m))
  ## if (!is.null(y)) m <- m[, -1, drop=FALSE]
  for (i in seq(along=m)) {
    if (is.ordered(m[[i]])) m[[i]] <- as.numeric(m[[i]])
  }
  return (list(m, y))
}