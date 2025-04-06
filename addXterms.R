#############################
# addXterms.R
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