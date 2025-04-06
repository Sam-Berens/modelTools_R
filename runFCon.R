#############################
# runFCon.R
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