###############################################################################
# 04_crosssection_dependence.R
# Pesaran (2004) CD test + Pesaran-Yamagata (2008) slope heterogeneity — Table 3
###############################################################################
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(plm); library(dplyr)
dir.create("output", showWarnings=FALSE)

panel  <- readRDS("data/panel_clean.rds")
pdata  <- readRDS("data/pdata.rds")
COUNTRIES <- unique(panel$country_iso)
N_ <- length(COUNTRIES); T_ <- 27

# ── CD test via plm::pcdtest ────────────────────────────────────────────────
cd_vars <- c("log_E","log_P","log_N","log_G","log_Y","log_M")
cd_results <- lapply(cd_vars, function(v) {
  fml <- as.formula(paste(v, "~ stats::lag(", v, ",1)"))
  m   <- suppressWarnings(plm(fml, data=pdata, model="within"))
  cd  <- pcdtest(m, test="cd")
  data.frame(Variable=v, CD_statistic=round(cd$statistic,3),
             p_value=round(cd$p.value,4),
             Decision=ifelse(cd$p.value<0.05,"Reject H0 (dep.)","Not reject"))
})
cd_table <- do.call(rbind, cd_results)
print(cd_table)

# ── Pesaran-Yamagata slope heterogeneity ────────────────────────────────────
het_formula <- log_E ~ log_N + log_N_sq + log_G + log_M + log_Y + log_K + L
k_ <- 7   # number of regressors (excl. intercept)

pooled_b <- coef(plm(het_formula, data=pdata, model="pooling"))[-1]

z_i <- sapply(COUNTRIES, function(co) {
  d   <- panel[panel$country_iso==co,]
  m_i <- lm(het_formula, data=d)
  b_i <- coef(m_i)[-1]; V_i <- vcov(m_i)[-1,-1]
  db  <- b_i - pooled_b
  tryCatch(as.numeric(t(db) %*% solve(V_i) %*% db), error=function(e) NA_real_)
})
z_i <- z_i[!is.na(z_i)]
N_  <- length(z_i)

delta     <- sqrt(N_) * (mean(z_i) - k_) / sqrt(2*k_)
adj_denom <- sqrt(2*k_*(T_-k_-1)/(T_-k_+1))
delta_adj <- sqrt(N_) * (mean(z_i) - k_*(T_-k_-1)/(T_-k_+1)) / adj_denom

cat("\nSlope heterogeneity (Pesaran-Yamagata 2008):\n")
cat(" Delta =", round(delta,3), " p =", round(1-pnorm(delta),4), "\n")
cat(" Delta_adj =", round(delta_adj,3), " p =", round(1-pnorm(delta_adj),4), "\n")

table3 <- rbind(
  cd_table,
  data.frame(Variable="Slope het. delta", CD_statistic=round(delta,3),
             p_value=round(1-pnorm(delta),4),
             Decision=ifelse(1-pnorm(delta)<0.05,"Reject H0 (het.)","Not reject")),
  data.frame(Variable="Slope het. delta_adj", CD_statistic=round(delta_adj,3),
             p_value=round(1-pnorm(delta_adj),4),
             Decision=ifelse(1-pnorm(delta_adj)<0.05,"Reject H0 (het.)","Not reject"))
)
print(table3); write.csv(table3,"output/table_3_cd_heterogeneity.csv",row.names=FALSE)
cat("Table 3 -> output/table_3_cd_heterogeneity.csv\n")
