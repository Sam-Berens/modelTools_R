#############################
# getXtermStruct.R
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