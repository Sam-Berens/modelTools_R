#############################
# runHCons.R
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

if (!exists("runFCon")) {
  source("./runFCon.R")
}