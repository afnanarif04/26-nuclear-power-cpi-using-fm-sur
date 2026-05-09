###############################################################################
# 06_fmsur_main.R
# FM-SUR primary estimator (Wagner, Grabarczyk & Hong 2020) — Table 6
# Implements: long-run covariance (Bartlett kernel), FM bias correction,
# country-specific long-run coefficients, system poolability Wald test.
###############################################################################
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(dplyr); library(Matrix)
dir.create("output", showWarnings=FALSE)

panel    <- readRDS("data/panel_clean.rds")
COUNTRIES <- c("FRA","KOR","SVK","HUN","CZE","FIN","CHE","SWE")
regs     <- c("log_N","log_N_sq","log_G","log_M","log_Y","log_K","L")

# ── Bartlett kernel long-run covariance ────────────────────────────────────
lr_cov <- function(mat, bw=NULL) {
  T_ <- nrow(mat); k <- ncol(mat)
  if (is.null(bw)) bw <- max(1L, floor(4*(T_/100)^(2/9)))
  G0 <- (1/T_) * t(mat) %*% mat
  L  <- matrix(0,k,k)
  for (j in seq_len(bw)) {
    w   <- 1 - j/(bw+1)
    Gj  <- (1/T_) * t(mat[1:(T_-j),,drop=FALSE]) %*% mat[(j+1):T_,,drop=FALSE]
    L   <- L + w*Gj
  }
  list(G0=G0, Lambda=L, Omega=G0+L+t(L))
}

# ── Per-country FM-OLS ─────────────────────────────────────────────────────
fmsur <- list()
for (co in COUNTRIES) {
  d  <- panel[panel$country_iso==co,]; d <- d[order(d$year),]
  T_ <- nrow(d)
  Y  <- d$log_E
  X  <- cbind(1, as.matrix(d[,regs]))
  k_ <- ncol(X)
  
  # OLS residuals
  beta_ols <- solve(t(X)%*%X) %*% (t(X)%*%Y)
  u_ols    <- Y - X%*%beta_ols
  
  # Long-run covariance of (residual, Delta_X)
  dX   <- apply(X[,-1,drop=FALSE], 2, diff)
  jmat <- cbind(u_ols[-1], dX)
  lrc  <- lr_cov(jmat)
  Om   <- lrc$Omega; Lam <- lrc$Lambda
  
  om_uu <- Om[1,1]; om_uX <- Om[1,-1]
  Om_XX <- Om[-1,-1]
  Om_XX_inv <- tryCatch(solve(Om_XX), error=function(e) MASS::ginv(Om_XX))
  
  # FM-corrected dependent variable and bias correction
  Y_star     <- Y[-1] - dX %*% Om_XX_inv %*% om_uX
  delta_plus <- Lam[1,-1] - om_uX %*% Om_XX_inv %*% Lam[-1,-1]
  bias_corr  <- c(0, delta_plus)
  
  X_post     <- X[-1,]
  XtX_inv    <- solve(t(X_post)%*%X_post)
  beta_fm    <- XtX_inv %*% (t(X_post)%*%Y_star - T_*bias_corr)
  
  # FM standard errors
  sigma2     <- max(1e-10, as.numeric(om_uu - om_uX%*%Om_XX_inv%*%om_uX))
  V_fm       <- sigma2 * XtX_inv
  se_fm      <- sqrt(pmax(0, diag(V_fm)))
  
  R2_        <- 1 - sum((Y-X%*%beta_ols)^2)/sum((Y-mean(Y))^2)
  
  fmsur[[co]] <- list(beta=as.numeric(beta_fm), se=se_fm, R2=round(R2_,4))
}

# ── Cross-country average & turning point ─────────────────────────────────
beta_mat <- do.call(rbind, lapply(fmsur, function(x) x$beta[-1]))
colnames(beta_mat) <- regs
avg_b1   <- mean(beta_mat[,"log_N"])
avg_b2   <- mean(beta_mat[,"log_N_sq"])
tp_logN  <- -avg_b1/(2*avg_b2)
tp_pct   <- exp(tp_logN)*100

cat("\n=== FM-SUR Long-Run Estimates ===\n")
cat(sprintf("%-14s %8s %8s %8s\n","Country","beta1_N","beta2_Nsq","R2"))
for (co in COUNTRIES) {
  b <- fmsur[[co]]$beta
  cat(sprintf("%-14s %8.4f %9.4f %8.4f\n", co, b[2], b[3], fmsur[[co]]$R2))
}
cat(sprintf("%-14s %8.4f %9.4f\n","[Average]",avg_b1,avg_b2))
cat(sprintf("Panel-average turning point: N* = %.2f%%\n", tp_pct))

# ── Poolability Wald test ─────────────────────────────────────────────────
# Compare country coefficients to cross-country average
W_pool <- sum(sapply(seq_along(COUNTRIES), function(i) {
  diff_ <- beta_mat[i,] - colMeans(beta_mat)
  se_   <- fmsur[[COUNTRIES[i]]]$se[-1]
  sum((diff_/se_)^2, na.rm=TRUE)
}))
df_W  <- length(COUNTRIES)*length(regs)
p_W   <- pchisq(W_pool, df_W, lower.tail=FALSE)
cat(sprintf("Poolability Wald: W=%.3f df=%d p=%.4f\n", W_pool, df_W, p_W))

# ── Save Table 6 ──────────────────────────────────────────────────────────
table6 <- data.frame(
  Country  = c(COUNTRIES,"Average"),
  beta1_N  = round(c(beta_mat[,"log_N"],   avg_b1), 4),
  beta2_Nsq= round(c(beta_mat[,"log_N_sq"],avg_b2), 4),
  gamma_G  = round(c(beta_mat[,"log_G"],   mean(beta_mat[,"log_G"])), 4),
  gamma_M  = round(c(beta_mat[,"log_M"],   mean(beta_mat[,"log_M"])), 4),
  gamma_Y  = round(c(beta_mat[,"log_Y"],   mean(beta_mat[,"log_Y"])), 4),
  gamma_K  = round(c(beta_mat[,"log_K"],   mean(beta_mat[,"log_K"])), 4),
  gamma_L  = round(c(beta_mat[,"L"],       mean(beta_mat[,"L"])),     4),
  R2       = c(sapply(fmsur,"[[","R2"), NA),
  tp_pct   = round(c(sapply(COUNTRIES,function(co){
                       b<-fmsur[[co]]$beta; exp(-b[2]/(2*b[3]))*100}), tp_pct),2)
)
print(table6); write.csv(table6,"output/table_6_fmsur.csv",row.names=FALSE)
cat("Table 6 -> output/table_6_fmsur.csv\n")

# Save for figure scripts
saveRDS(fmsur, "data/fmsur_coefs.rds")
cat("FM-SUR coefficients saved to data/fmsur_coefs.rds\n")
cat("NOTE: Use beta1_N and beta2_Nsq from table6 to update figure_3_fmsur_results.R\n")
