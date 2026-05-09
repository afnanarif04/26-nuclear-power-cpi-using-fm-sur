###############################################################################
# 09_fmsur_ardl.R
# Robustness 3: Bounds-test ARDL (Pesaran, Shin & Smith 2001) — Table 7 row
# Country-specific ARDL(1,1,...) conditional error-correction; averages across
# countries for Table 7. Accommodates mixed I(0)/I(1) regressors.
###############################################################################
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(dplyr); library(car)
dir.create("output", showWarnings=FALSE)

panel    <- readRDS("data/panel_clean.rds")
COUNTRIES <- c("FRA","KOR","SVK","HUN","CZE","FIN","CHE","SWE")
regs     <- c("log_N","log_N_sq","log_G","log_M","log_Y","log_K","L")

ardl_lr <- matrix(NA, length(COUNTRIES), length(regs),
                   dimnames=list(COUNTRIES,regs))
ardl_F  <- numeric(length(COUNTRIES))

for (i in seq_along(COUNTRIES)) {
  co <- COUNTRIES[i]
  d  <- panel[panel$country_iso==co,]; d<-d[order(d$year),]; T_<-nrow(d)
  
  # ECM representation: Delta_y = rho*y_lag + theta*X_lag + gamma*Delta_X + eps
  dY   <- c(NA,diff(d$log_E));    Yl <- c(NA,head(d$log_E,-1))
  dmat <- apply(as.matrix(d[,regs]),2,diff)
  dmat_full <- rbind(NA,dmat)
  Xlags <- apply(as.matrix(d[,regs]),2,function(x) c(NA,head(x,-1)))
  colnames(Xlags)     <- paste0(regs,"_l")
  colnames(dmat_full) <- paste0(regs,"_d")
  df_ <- data.frame(dY=dY, Yl=Yl, dmat_full, Xlags) |> na.omit()
  
  rhs_l <- paste(paste0(regs,"_l"), collapse="+")
  rhs_d <- paste(paste0(regs,"_d"), collapse="+")
  fml   <- as.formula(paste("dY ~ Yl +", rhs_l, "+", rhs_d))
  
  m <- tryCatch(lm(fml, data=df_), error=function(e) NULL)
  if (is.null(m)) next
  
  rho <- coef(m)["Yl"]
  if (!is.na(rho) && abs(rho) > 1e-8) {
    theta <- coef(m)[paste0(regs,"_l")]
    ardl_lr[i,] <- -theta / rho
  }
  
  # Bounds F-test: joint significance of Yl and all X_lags
  lv_terms <- c("Yl", paste0(regs,"_l"))
  ardl_F[i] <- tryCatch({
    lhyp <- paste0(lv_terms," = 0")
    lh   <- linearHypothesis(m, lhyp, test="F", white.adjust=FALSE)
    as.numeric(lh$F[2])
  }, error=function(e) NA_real_)
}

avg_lr  <- colMeans(ardl_lr, na.rm=TRUE)
avg_b1  <- avg_lr["log_N"]; avg_b2 <- avg_lr["log_N_sq"]
avg_tp  <- exp(-avg_b1/(2*avg_b2))*100
avg_F   <- mean(ardl_F, na.rm=TRUE)

cat(sprintf("Bounds-test avg: beta1=%.4f beta2=%.4f N*=%.2f%% avg_F=%.3f\n",
            avg_b1,avg_b2,avg_tp,avg_F))
cat("Country F-stats:\n"); print(setNames(round(ardl_F,3),COUNTRIES))

row7c <- data.frame(Estimator="Bounds-test SUR (Pesaran et al. 2001)",
  beta1=round(avg_b1,4), se1=NA,
  beta2=round(avg_b2,4), se2=NA,
  turning_point_pct=round(avg_tp,2), avg_bounds_F=round(avg_F,3))
write.csv(row7c,"output/table_7_row_bounds_ardl.csv",row.names=FALSE)
cat("-> output/table_7_row_bounds_ardl.csv\n")

# Combined Table 7
t7a <- tryCatch(read.csv("output/table_7_row_groupmean_fmols.csv"),error=function(e)NULL)
t7b <- tryCatch(read.csv("output/table_7_row_pooled_cpr.csv"),error=function(e)NULL)
fm6 <- tryCatch(read.csv("output/table_6_fmsur.csv"),error=function(e)NULL)
if (!is.null(fm6)) {
  avg_row <- fm6[fm6$Country=="Average",]
  t7_fmsur <- data.frame(Estimator="FM-SUR (Wagner et al. 2020)",
    beta1=avg_row$beta1_N, se1=NA, beta2=avg_row$beta2_Nsq, se2=NA,
    turning_point_pct=avg_row$tp_pct)
  all_rows <- rbind(t7_fmsur,t7a,t7b,row7c[,names(t7_fmsur)])
  write.csv(all_rows,"output/table_7_robustness.csv",row.names=FALSE)
  cat("Combined Table 7 -> output/table_7_robustness.csv\n")
  print(all_rows)
}
