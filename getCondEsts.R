#############################
# getCondEsts.R
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

if (!exists("getXEUL")) {
  source("./getXEUL.R")
}