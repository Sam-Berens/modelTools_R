#############################
# getXEUL.R
#############################
getXEUL <- function(mdl, ...) {
  # mdl is assumed to be a list with elements:
  #   Coefficients (a data.frame with column Estimate),
  #   CoefficientCovariance (a matrix),
  #   DFE (degrees of freedom for error),
  #   CoefficientNames (a character vector),
  #   and optionally Link (with an element $Inverse)
  predictors <- list(...)
  tryCatch(b <- fixef(mdl), error = function(e) b <- coef(mdl))
  d <- length(b)
  dIn <- length(predictors)
  coefName <- names(b)
  numSimple <- sum(!grepl(":", coefName))
  if ((dIn < d) && (dIn != numSimple)) {
    stop("You must provide either all the model terms or only the simple effects.")
  }
  n <- max(sapply(predictors, length))
  X <- matrix(0, nrow = n, ncol = d)
  for (ii in 1:dIn) {
    vec <- predictors[[ii]]
    if (length(vec) == 1) {
      X[, ii] <- rep(vec, n)
    } else {
      if (length(vec) == n) {
        X[, ii] <- vec
      } else {
        stop("Incompatible input dimensions.")
      }
    }
  }
  # Add any remaining interaction terms.
  xTermStruct <- getXtermStruct(coefName)
  X <- addXterms(X, xTermStruct)
  
  C <- vcov(mdl)
  df <- df.residual(mdl)
  est <- X %*% b
  se <- sqrt(rowSums((X %*% C) * X))
  tval_low <- qt(0.025, df)
  tval_upp <- qt(0.975, df)
  Ci95 <- cbind(est + se * tval_low, est + se * tval_upp)
  
  tryCatch(invLink <- mdl@resp$family$linkinv, error = function(e) invLink <- function(x) x)
  est <- invLink(est)
  Ci95 <- invLink(Ci95)
  
  return(list(X = X,
              est = as.numeric(est),
              low = as.numeric(Ci95[, 1]),
              upp = as.numeric(Ci95[, 2])))
}

if (!exists("getXtermStruct")) {
  source("./getXtermStruct.R")
}

if (!exists("addXterms")) {
  source("./addXterms.R")
}