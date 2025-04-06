#############################
# getTreatWeights.R
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

if (!exists("getXtermStruct")) {
  source("./getXtermStruct.R")
}

if (!exists("addXterms")) {
  source("./addXterms.R")
}