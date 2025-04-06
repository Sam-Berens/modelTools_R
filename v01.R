#############################
# 1. runHCons.R
#############################
runHCons <- function(mdl, conTable) {
  # conTable is assumed to be a data.frame whose row names indicate effect names
  # and which contains a single list-column "H" with contrast matrices.
  conName <- rownames(conTable)
  nCon <- length(conName)
  testType <- character(nCon)
  testStat <- rep(NA, nCon)
  pValue <- rep(NA, nCon)
  
  # Check that the contrast column is a list
  if (!is.list(conTable$H)) {
    stop("The conTable input must have a single variable corresponding to the desired contrasts.")
  }
  
  for (iCon in seq_len(nCon)) {
    res <- runFCon(conTable$H[[iCon]], mdl)
    pValue[iCon] <- res$pValue
    fStat <- res$fStat
    df1 <- res$df1
    df2 <- res$df2
    rowSign <- res$rowSign
    if (df1 == 1) {
      testType[iCon] <- sprintf("t(%d)", df2)
      testStat[iCon] <- sqrt(fStat) * rowSign
    } else {
      testType[iCon] <- sprintf("F(%d,%d)", df1, df2)
      testStat[iCon] <- fStat
    }
  }
  conTable$testType <- testType
  conTable$testStat <- testStat
  conTable$pValue <- pValue
  return(conTable)
}

#############################
# 2. addXterms.R
#############################
addXterms <- function(X, xTermStruct) {
  # xTermStruct is a list where each element is a vector of column indices
  # that should be multiplied together to yield a (possibly interaction) term.
  k <- length(xTermStruct)
  n <- nrow(X)
  # Pad X with additional columns if needed
  if (ncol(X) < k) {
    X <- cbind(X, matrix(NA, n, k - ncol(X)))
  }
  for (iK in seq_len(k)) {
    idx <- xTermStruct[[iK]]
    # Compute the product of the specified columns row-wise.
    X[, iK] <- apply(X[, idx, drop = FALSE], 1, prod)
  }
  return(X)
}

#############################
# 3. getCondEsts.R
#############################
getCondEsts <- function(mdl, condXs) {
  # condXs is a data.frame with row names (the condition labels)
  # and a column "X" containing condition weight vectors (each as a numeric vector).
  condName <- rownames(condXs)
  n <- nrow(condXs)
  est <- rep(NA, n)
  low <- rep(NA, n)
  upp <- rep(NA, n)
  
  # If condXs$X is not already a list, convert each row to a vector.
  if (!is.list(condXs$X)) {
    Xs <- split(as.matrix(condXs$X), seq_len(n))
  } else {
    Xs <- condXs$X
  }
  
  for (i in seq_len(n)) {
    # Pass the condition vector as separate arguments via do.call.
    res <- do.call(getXEUL, c(list(mdl), as.list(Xs[[i]])))
    est[i] <- res$est
    low[i] <- res$low
    upp[i] <- res$upp
  }
  condXs <- data.frame(est = est, low = low, upp = upp, row.names = condName)
  return(condXs)
}

#############################
# 4. getTreatWeights.R
#############################
getTreatWeights <- function(coefName, factorStruct) {
  # If no inputs are provided then run in test mode.
  if (missing(coefName)) {
    # Test mode (example values)
    coefName <- c("(Intercept)", "Session01", "Session02", "Activedrug",
                  "Complex", "Complexity", "Session01:Activedrug",
                  "Session02:Activedrug", "Session01:Complex",
                  "Session02:Complex", "Activedrug:Complex",
                  "Session01:Complexity", "Session02:Complexity",
                  "Activedrug:Complexity", "Session01:Activedrug:Complex",
                  "Session02:Activedrug:Complex", "Session01:Activedrug:Complexity",
                  "Session02:Activedrug:Complexity")
    factorStruct <- list(
      list(name = "Session", coding = "bin", coefs = c(2, 3)),
      list(name = "Treatment", coding = "bin", coefs = 4),
      list(name = "Trialtype", coding = "bin", coefs = 5),
      list(name = "Complexity", coding = "lin", coefs = 6, centre = 0)
    )
  }
  if (missing(factorStruct))
    stop("Two inputs are required: coefName and a valid factorStruct.")
  
  # Validate factorStruct: each factor's coding must be 'bin' or 'lin'.
  for (fs in factorStruct) {
    if (!(tolower(fs$coding) %in% c("bin", "lin"))) {
      stop("Unrecognised coding for one of the factors.")
    }
    if (tolower(fs$coding) == "lin" && length(fs$coefs) > 1)
      stop("This function only supports binary coding for factors with more than 2 levels.")
    if (tolower(fs$coding) == "lin" && is.null(fs$centre))
      stop("Mean values (centre) must be specified for linearly coded factors.")
  }
  
  k <- length(coefName)
  f <- length(factorStruct)
  l <- sapply(factorStruct, function(fs) length(fs$coefs) + 1)
  xTermStruct <- getXtermStruct(coefName)
  Xidx <- getXidx(l, f)
  X <- getX(factorStruct, k, xTermStruct, Xidx)
  G_main <- getG_main(factorStruct, X)
  
  linFacs <- which(tolower(sapply(factorStruct, function(fs) fs$coding)) == "lin")
  iH <- 0
  Hname <- vector("list", 2^f - 1)
  Hmat <- vector("list", 2^f - 1)
  
  for (iCol in seq_len(f)) {
    # Get all combinations of factors of size iCol.
    combMat <- combn(seq_len(f), iCol, simplify = FALSE)
    for (t in combMat) {
      fxName <- paste(sapply(t, function(i) factorStruct[[i]]$name), collapse = ":")
      # Build contrast from G_main matrices via successive cross-products.
      for (iT in seq_along(t)) {
        if (iT == 1) {
          gg <- G_main[[t[iT]]]
        } else {
          gg <- crossGs(gg, G_main[[t[iT]]])
        }
      }
      # Determine which linear factors (if any) to scale.
      facsToScale <- linFacs[!linFacs %in% t]
      if (!rlang::is_empty(facsToScale)) {
        scaleRes <- getCoefsToScale(coefName, factorStruct, facsToScale)
        coefsToScale <- scaleRes$coefsToScale
        centres <- scaleRes$means
      } else {
        coefsToScale <- list()
        centres <- list()
      }
      SX <- X
      if (length(coefsToScale) > 0) {
        # Multiply each designated column by its centre value.
        for (j in seq_along(coefsToScale)) {
          SX[, coefsToScale[j]] <- SX[, coefsToScale[j]] * centres[j]
        }
      }
      H <- gg %*% SX
      if (sum(H^2) > 0) {
        iH <- iH + 1
        # Scale H by its smallest nonzero absolute value
        nonzero <- abs(H)[abs(H) > 0]
        if (length(nonzero) > 0)
          H <- H / min(nonzero)
        Hname[[iH]] <- fxName
        Hmat[[iH]] <- H
      }
    }
  }
  # Remove unused (NULL) entries
  Hname <- unlist(Hname[1:iH])
  Hmat <- Hmat[1:iH]
  anovaHs <- data.frame(H = I(Hmat), row.names = Hname)
  
  # --- Construct condition weight tables ---
  # Select rows in Xidx that have level 1 for all linear factors.
  rowsToSel <- which(apply(Xidx[, linFacs, drop = FALSE], 1, function(x) all(x == 1)))
  # Only select columns corresponding to main effects (where xTermStruct elements are of length 1)
  mainEff <- which(sapply(xTermStruct, function(x) length(x) == 1))
  X_reduced <- X[rowsToSel, mainEff, drop = FALSE]
  
  # Generate condition names by pasting each non-linear factor's name and level.
  nonLinFacs <- setdiff(seq_len(f), linFacs)
  Xnames <- apply(Xidx[rowsToSel, nonLinFacs, drop = FALSE], 1, function(r) {
    paste(mapply(function(i, level) {
      paste(factorStruct[[i]]$name, level, sep = "")
    }, nonLinFacs, r), collapse = ":")
  })
  
  # For each linear factor, replace the designated columns with their centre value.
  for (i in seq_along(linFacs)) {
    iCol <- factorStruct[[linFacs[i]]]$coefs
    centre <- factorStruct[[linFacs[i]]]$centre
    X_reduced[, iCol] <- centre
  }
  
  # Prepare linXs: for each linear factor, copy X_reduced and replace the corresponding
  # column with 1i (complex 1) so that later arithmetic works as in MATLAB.
  linXs <- list()
  for (i in seq_along(linFacs)) {
    iCol <- factorStruct[[linFacs[i]]]$coefs
    linX <- X_reduced
    linX[, iCol] <- 1i
    linXs[[ factorStruct[[linFacs[i]]]$name ]] <- linX
  }
  
  # Add interaction terms.
  X_final <- addXterms(X_reduced, xTermStruct)
  for (i in seq_along(linXs)) {
    linXs[[i]] <- addXterms(linXs[[i]], xTermStruct)
    # Replace any complex 1i values with NA.
    linXs[[i]][linXs[[i]] == 1i] <- NA
    # Convert each row into a list element if desired.
    linXs[[i]] <- data.frame(X = I(split(linXs[[i]], seq_len(nrow(linXs[[i]])))),
                             row.names = Xnames)
  }
  
  condXs <- data.frame(X = I(split(X_final, seq_len(nrow(X_final)))),
                       row.names = Xnames)
  
  # Return a list with the three objects.
  return(list(anovaHs = anovaHs, condXs = condXs, linXs = linXs))
}

#############################
# Helper functions for getTreatWeights
#############################

getXidx <- function(l, f) {
  m <- prod(l)
  Xid <- matrix(NA, nrow = m, ncol = f)
  scale <- rev(cumprod(rev(c(l[-1], 1))))
  for (iM in 1:m) {
    c_val <- iM - 1
    for (ig in 1:f) {
      n_val <- floor(c_val / scale[ig])
      Xid[iM, ig] <- n_val
      c_val <- c_val - n_val * scale[ig]
    }
  }
  return(Xid)
}

getX <- function(factorStruct, k, xTermStruct, Xidx) {
  f <- length(factorStruct)
  m <- nrow(Xidx)
  X <- matrix(NA, nrow = m, ncol = k)
  X[, 1] <- 1
  for (iX in 1:m) {
    for (iFac in 1:f) {
      level <- Xidx[iX, iFac]
      ffx <- factorStruct[[iFac]]$coefs
      for (iffx in seq_along(ffx)) {
        if (level == 0) {
          X[iX, ffx[iffx]] <- 0
        } else if (level == iffx) {
          X[iX, ffx[iffx]] <- 1
        } else {
          X[iX, ffx[iffx]] <- 0
        }
      }
    }
  }
  X <- addXterms(X, xTermStruct)
  return(X)
}

getG_main <- function(factorStruct, X) {
  f <- length(factorStruct)
  G_main <- vector("list", f)
  for (iFac in 1:f) {
    coefs <- factorStruct[[iFac]]$coefs
    nRows <- length(coefs)
    G_current <- matrix(NA, nrow = nRows, ncol = nrow(X))
    for (iRow in 1:nRows) {
      for (iX in 1:nrow(X)) {
        row_val <- X[iX, ]
        if (sum(row_val[coefs]) == 0) {
          G_current[iRow, iX] <- -1
        } else if (row_val[coefs[iRow]] == 1) {
          G_current[iRow, iX] <- 1
        } else {
          G_current[iRow, iX] <- 0
        }
      }
    }
    G_main[[iFac]] <- G_current
  }
  return(G_main)
}

getCoefsToScale <- function(coefName, factorStruct, facsToScale) {
  simpleCoeffsToScale <- unlist(lapply(facsToScale, function(i) factorStruct[[i]]$coefs))
  if (length(simpleCoeffsToScale) > 0) {
    simpleMeans <- unlist(lapply(facsToScale, function(i) factorStruct[[i]]$centre))
  } else {
    simpleMeans <- numeric(0)
  }
  coefsToScale <- c()
  means <- c()
  for (ii in seq_along(simpleCoeffsToScale)) {
    term <- coefName[simpleCoeffsToScale[ii]]
    terms <- grep(term, coefName)
    coefsToScale <- c(coefsToScale, terms)
    means <- c(means, rep(simpleMeans[ii], length(terms)))
  }
  ord <- order(coefsToScale)
  coefsToScale <- coefsToScale[ord]
  means <- means[ord]
  return(list(coefsToScale = coefsToScale, means = means))
}

crossGs <- function(A, B) {
  # For each column j, compute the outer product and flatten.
  n <- nrow(A)
  k <- ncol(A)
  m <- nrow(B)
  C <- matrix(NA, nrow = n * m, ncol = k)
  for (j in 1:k) {
    outer_prod <- outer(A[, j], B[, j])
    C[, j] <- as.vector(outer_prod)
  }
  return(C)
}

#############################
# 5. getXEUL.R
#############################
getXEUL <- function(mdl, ...) {
  # mdl is assumed to be a list with elements:
  #   Coefficients (a data.frame with column Estimate),
  #   CoefficientCovariance (a matrix),
  #   DFE (degrees of freedom for error),
  #   CoefficientNames (a character vector),
  #   and optionally Link (with an element $Inverse)
  predictors <- list(...)
  d <- length(mdl$Coefficients$Estimate)
  dIn <- length(predictors)
  coefName <- mdl$CoefficientNames
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
  
  b <- mdl$Coefficients$Estimate
  C <- mdl$CoefficientCovariance
  df <- mdl$DFE
  est <- X %*% b
  se <- sqrt(rowSums((X %*% C) * X))
  tval_low <- qt(0.025, df)
  tval_upp <- qt(0.975, df)
  Ci95 <- cbind(est + se * tval_low, est + se * tval_upp)
  
  if (!is.null(mdl$Link)) {
    invLink <- mdl$Link$Inverse
  } else {
    invLink <- function(x) x
  }
  est <- invLink(est)
  Ci95 <- invLink(Ci95)
  
  return(list(X = X,
              est = as.numeric(est),
              low = as.numeric(Ci95[, 1]),
              upp = as.numeric(Ci95[, 2])))
}

#############################
# 6. getXtermStruct.R
#############################
getXtermStruct <- function(coefName) {
  isXterm <- grepl(":", coefName)
  xTermStruct <- vector("list", length(coefName))
  for (i in seq_along(coefName)) {
    if (isXterm[i]) {
      parts <- unlist(strsplit(coefName[i], ":"))
      # Find indices of the main effects (i.e. those without ':')
      mainEffectIndices <- which(!grepl(":", coefName))
      # For each part, find its match among the main effects.
      indices <- sapply(parts, function(p) match(p, coefName[mainEffectIndices]))
      xTermStruct[[i]] <- indices
    } else {
      xTermStruct[[i]] <- i
    }
  }
  return(xTermStruct)
}

#############################
# 7. runExample.R
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
       xlab = "Complexity", ylab = "y")
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

#############################
# 8. runFCon.R
#############################
runFCon <- function(H, glme) {
  b <- glme$Coefficients$Estimate
  C <- glme$CoefficientCovariance
  df1 <- qr(H)$rank
  df2 <- glme$DFE
  Hb <- as.vector(H %*% b)
  HCHt <- H %*% C %*% t(H)
  fStat <- as.numeric((t(Hb) %*% solve(HCHt) %*% Hb) / df1)
  pValue <- 1 - pf(fStat, df1, df2)
  rowSign <- sign(Hb)
  return(list(pValue = pValue, fStat = fStat, df1 = df1, df2 = df2, rowSign = rowSign))
}
