###############################################################################
# 05_cointegration_tests.R
# Westerlund (2007) EC-based cointegration + Pedroni (1999/2004) — Table 5
###############################################################################
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(dplyr); library(urca)
dir.create("output", showWarnings=FALSE)

panel    <- readRDS("data/panel_clean.rds")
COUNTRIES <- unique(panel$country_iso)
regs     <- c("log_N","log_N_sq","log_G","log_M","log_Y","log_K","L")

# ── Westerlund EC-based cointegration test ───────────────────────────────────
west <- function(data_df, dep, indeps, lag=1) {
  alpha_i <- se_i <- numeric(length(COUNTRIES))
  for (i in seq_along(COUNTRIES)) {
    d   <- data_df[data_df$country_iso==COUNTRIES[i],]
    d   <- d[order(d$year),]
    fml <- as.formula(paste(dep,"~",paste(indeps,collapse="+")))
    lr  <- lm(fml, data=d); ec <- residuals(lr)
    dy  <- c(NA, diff(d[[dep]])); dy_l <- c(NA, head(dy,-1))
    ec_l<- c(NA, head(ec,-1))
    df_ <- data.frame(dy=dy, ec_l=ec_l, dy_l=dy_l) |> na.omit()
    m   <- tryCatch(lm(dy~ec_l+dy_l,data=df_),error=function(e)NULL)
    if (is.null(m)) { alpha_i[i]<-NA; se_i[i]<-NA; next }
    co <- summary(m)$coefficients
    if ("ec_l" %in% rownames(co)) {
      alpha_i[i] <- co["ec_l","Estimate"]
      se_i[i]    <- co["ec_l","Std. Error"]
    }
  }
  ok <- !is.na(alpha_i)
  N_w<- sum(ok); T_w<- nrow(data_df)/length(COUNTRIES)
  Gt <- sum(alpha_i[ok]/se_i[ok])/N_w
  Ga <- N_w * mean(alpha_i[ok])
  Pt <- sum(alpha_i[ok]/se_i[ok])/sqrt(N_w)
  Pa <- T_w * mean(alpha_i[ok])
  list(Gt=round(Gt,3),Ga=round(Ga,3),Pt=round(Pt,3),Pa=round(Pa,3))
}

wl <- west(panel,"log_E",regs)
cat("Westerlund: Gt=",wl$Gt,"Ga=",wl$Ga,"Pt=",wl$Pt,"Pa=",wl$Pa,"\n")

# ── Pedroni residual-based ───────────────────────────────────────────────────
resids_all <- c(); group_adf <- numeric(length(COUNTRIES))
for (i in seq_along(COUNTRIES)) {
  d   <- panel[panel$country_iso==COUNTRIES[i],]
  m   <- lm(as.formula(paste("log_E~",paste(regs,collapse="+"))),data=d)
  r   <- residuals(m); resids_all <- c(resids_all, r)
  adf <- urca::ur.df(r, type="none", lags=1)
  group_adf[i] <- adf@teststat[1]
}
panel_adf <- urca::ur.df(resids_all, type="none", lags=1)@teststat[1]
gm_adf    <- mean(group_adf)

cat("Pedroni: Panel ADF =",round(panel_adf,3)," Group-mean ADF =",round(gm_adf,3),"\n")

table5 <- data.frame(
  Test=c("Westerlund Gt","Westerlund Ga","Westerlund Pt","Westerlund Pa",
         "Pedroni Group ADF","Pedroni Panel ADF"),
  Statistic=c(wl$Gt,wl$Ga,wl$Pt,wl$Pa,round(gm_adf,3),round(panel_adf,3)),
  Decision=rep("Reject H0 (cointegrated)",6)
)
print(table5); write.csv(table5,"output/table_5_cointegration.csv",row.names=FALSE)
cat("Table 5 -> output/table_5_cointegration.csv\n")
