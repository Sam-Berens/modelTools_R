#############################
# runFCon.R
#############################
runFCon <- function(H, mdl) {
  tryCatch(b <- fixef(mdl), error = function(e) b <- coef(mdl))
  C <- vcov(mdl)
  df1 <- qr(H)$rank
  df2 <- df.residual(mdl)
  Hb <- as.vector(H %*% b)
  HCHt <- H %*% C %*% t(H)
  fStat <- as.numeric((t(Hb) %*% solve(HCHt) %*% Hb) / df1)
  pValue <- 1 - pf(fStat, df1, df2)
  rowSign <- sign(Hb)
  return(list(pValue = pValue, fStat = fStat, df1 = df1, df2 = df2, rowSign = rowSign))
}