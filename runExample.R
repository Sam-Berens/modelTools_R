#############################
# runExample.R
#############################
runExample <- function() {
  # Obtain an example model and factor structure.
  model_and_factors <- getModel()
  mdl <- model_and_factors$mdl
  factorStruct <- model_and_factors$factorStruct
  coefName <- mdl$CoefficientNames
  
  # Compute ANOVA contrast vectors and condition weights.
  treatRes <- getTreatWeights(coefName, factorStruct)
  anovaHs <- treatRes$anovaHs
  condXs <- treatRes$condXs
  linXs <- treatRes$linXs
  
  # Run hypothesis contrasts and compute condition estimates.
  anovaHs <- runHCons(mdl, anovaHs)
  condEsts <- getCondEsts(mdl, condXs)
  
  # Plot the average effect of Complexity as a linear predictor.
  domain <- seq(-4, 4, by = 0.1)
  # For demonstration, we create a design matrix X with constant values
  # (in practice you would compute this from your linXs$Complexity data).
  w <- kronecker(rep(1/6, 6), c(0, 1))
  meanX <- w %*% Re(do.call(rbind,linXs$Complexity$X))
  meanX <- matrix(meanX[1:6],1,6)
  X <- meanX[rep(1:nrow(meanX), times = length(domain)), ]
  # Replace any NA values in X with the corresponding domain value.
  for (i in 1:nrow(X)) {
    X[i, is.na(X[i, ])] <- domain[i]
  }
  predictors <- split(X, col(X))
  xeulRes <- do.call(getXEUL, c(list(mdl), predictors))
  est <- xeulRes$est
  low <- xeulRes$low
  upp <- xeulRes$upp
  
  # Create the plot.
  plot(domain, est, type = "l", lwd = 2, col = "black",
       xlab = "Complexity", ylab = "y", ylim = c(-1, 1))
  lines(domain, low, lwd = 2, col = "red")
  lines(domain, upp, lwd = 2, col = "blue")
  legend("topright", legend = c("est", "low", "upp"),
         col = c("black", "red", "blue"), lwd = 2)
  
  # Report results.
  cat("Model Summary:\n")
  print(summary(mdl))
  cat("\nANOVA Contrasts:\n")
  print(anovaHs)
  cat("\nCondition Estimates:\n")
  print(condEsts)
}

# Helper: getModel (used by runExample)
getModel <- function() {
  # Define the model formula.
  formula <- y ~ Session * ActiveDrug * (Complex + Complexity)
  factorStruct <- list(
    list(name = "Session", coding = "bin", coefs = c(2, 3)),
    list(name = "Treatment", coding = "bin", coefs = 4),
    list(name = "Trialtype", coding = "bin", coefs = 5),
    list(name = "Complexity", coding = "lin", coefs = 6, centre = 0)
  )
  n <- 12^3
  Session <- sample(0:2, n, replace = TRUE)
  Session01 <- as.numeric(Session == 1)
  Session02 <- as.numeric(Session == 2)
  Session_cat <- factor(Session)
  ActiveDrug <- sample(0:1, n, replace = TRUE)
  Complex <- sample(0:1, n, replace = TRUE)
  Complexity <- rep(NA, n)
  for (i in 1:n) {
    if (Complex[i] == 1) {
      Complexity[i] <- rnorm(1)
    }
  }
  nanzscore <- function(v) (v - mean(v, na.rm = TRUE)) / sd(v, na.rm = TRUE)
  Complexity <- nanzscore(Complexity)
  Complexity[is.na(Complexity)] <- 0
  y <- rnorm(n)
  DataTable <- data.frame(Session = Session_cat,
                          Session01 = Session01,
                          Session02 = Session02,
                          ActiveDrug = ActiveDrug,
                          Complex = Complex,
                          Complexity = Complexity,
                          y = y)
  mdl <- lm(formula, data = DataTable)
  # Add additional fields to mimic the MATLAB model structure.
  mdl$CoefficientNames <- names(coef(mdl))
  mdl$Coefficients <- data.frame(Estimate = coef(mdl))
  mdl$CoefficientCovariance <- vcov(mdl)
  mdl$DFE <- df.residual(mdl)
  mdl$Link <- NULL  # or provide a list with an Inverse function if needed
  return(list(mdl = mdl, factorStruct = factorStruct))
}

if (!exists("getTreatWeights")) {
  source("./getTreatWeights.R")
}

if (!exists("runHCons")) {
  source("./runHCons.R")
}

if (!exists("getCondEsts")) {
  source("./getCondEsts.R")
}

if (!exists("getXEUL")) {
  source("./getXEUL.R")
}