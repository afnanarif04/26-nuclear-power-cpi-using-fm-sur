###############################################################################
# 07_groupmean_fmols.R
# Robustness 1: Group-mean FMOLS (Wagner & Reichold 2023) — Table 7 row
###############################################################################
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(dplyr)
dir.create("output", showWarnings=FALSE)

panel    <- readRDS("data/panel_clean.rds")
COUNTRIES <- c("FRA","KOR","SVK","HUN","CZE","FIN","CHE","SWE")
regs     <- c("log_N","log_N_sq","log_G","log_M","log_Y","log_K","L")
N_       <- length(COUNTRIES)

lr_cov1 <- function(mat, bw=NULL) {
  T_<-nrow(mat); k<-ncol(mat)
  if(is.null(bw)) bw<-max(1L,floor(4*(T_/100)^(2/9)))
  G0<-(1/T_)*t(mat)%*%mat; L<-matrix(0,k,k)
  for(j in seq_len(bw)){ w<-1-j/(bw+1)
    L<-L+w*(1/T_)*t(mat[1:(T_-j),,drop=FALSE])%*%mat[(j+1):T_,,drop=FALSE] }
  list(G0=G0,Lambda=L,Omega=G0+L+t(L))
}

beta_list <- list()
for (co in COUNTRIES) {
  d  <- panel[panel$country_iso==co,]; d<-d[order(d$year),]; T_<-nrow(d)
  Y  <- d$log_E; X <- cbind(1,as.matrix(d[,regs])); k_<-ncol(X)
  bo <- solve(t(X)%*%X)%*%(t(X)%*%Y); uo <- Y-X%*%bo
  dX <- apply(X[,-1,drop=FALSE],2,diff)
  lrc<- lr_cov1(cbind(uo[-1],dX))
  Om <-lrc$Omega; Lam<-lrc$Lambda
  om_uu<-Om[1,1]; om_uX<-Om[1,-1]; OmXX<-Om[-1,-1]
  OmXXi<-tryCatch(solve(OmXX),error=function(e) diag(1/diag(OmXX)))
  Ys    <- Y[-1]-dX%*%OmXXi%*%om_uX
  bias  <- c(0, Lam[1,-1]-om_uX%*%OmXXi%*%Lam[-1,-1])
  Xp    <- X[-1,]
  bf    <- solve(t(Xp)%*%Xp)%*%(t(Xp)%*%Ys - T_*bias)
  beta_list[[co]] <- as.numeric(bf)[-1]
}
bm      <- do.call(rbind, beta_list); colnames(bm) <- regs
gm_avg  <- colMeans(bm)
gm_se   <- apply(bm,2,sd)/sqrt(N_)
gm_b1   <- gm_avg["log_N"]; gm_b2 <- gm_avg["log_N_sq"]
gm_tp   <- exp(-gm_b1/(2*gm_b2))*100

cat(sprintf("Group-mean FMOLS: beta1=%.4f (se=%.4f) beta2=%.4f (se=%.4f) N*=%.2f%%\n",
            gm_b1,gm_se["log_N"],gm_b2,gm_se["log_N_sq"],gm_tp))

row7a <- data.frame(Estimator="Group-mean FMOLS (Wagner & Reichold 2023)",
  beta1=round(gm_b1,4),se1=round(gm_se["log_N"],4),
  beta2=round(gm_b2,4),se2=round(gm_se["log_N_sq"],4),
  turning_point_pct=round(gm_tp,2))
write.csv(row7a,"output/table_7_row_groupmean_fmols.csv",row.names=FALSE)
cat("-> output/table_7_row_groupmean_fmols.csv\n")
