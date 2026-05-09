###############################################################################
# 08_pooled_cpr.R
# Robustness 2: Pooled CPR with CCE (de Jong & Wagner 2022) — Table 7 row
###############################################################################
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(plm); library(dplyr); library(sandwich)
dir.create("output", showWarnings=FALSE)

panel    <- readRDS("data/panel_clean.rds")
pdata    <- readRDS("data/pdata.rds")

# Cross-section means (CCE augmentation)
cs <- panel %>% group_by(year) %>%
  summarise(across(c(log_E,log_N,log_N_sq,log_G,log_M,log_Y,log_K,L),
                   ~mean(.x,na.rm=TRUE),.names="cs_{.col}"),.groups="drop")
pa <- left_join(panel, cs, by="year")
pd2 <- pdata.frame(pa, index=c("country_iso","year"))

fml <- log_E ~ log_N + log_N_sq + log_G + log_M + log_Y + log_K + L +
  cs_log_E + cs_log_N + cs_log_N_sq + cs_log_G +
  cs_log_M + cs_log_Y + cs_log_K + cs_L

m   <- plm(fml, data=pd2, model="pooling")
vcv <- vcovHC(m, type="HC3", cluster="group", method="arellano")
se_ <- sqrt(diag(vcv))
b   <- coef(m)
b1  <- b["log_N"]; b2 <- b["log_N_sq"]
tp  <- exp(-b1/(2*b2))*100

cat(sprintf("Pooled CPR (CCE): beta1=%.4f (se=%.4f) beta2=%.4f (se=%.4f) N*=%.2f%%\n",
            b1,se_["log_N"],b2,se_["log_N_sq"],tp))

row7b <- data.frame(Estimator="Pooled CPR/CCE (de Jong & Wagner 2022)",
  beta1=round(b1,4), se1=round(se_["log_N"],4),
  beta2=round(b2,4), se2=round(se_["log_N_sq"],4),
  turning_point_pct=round(tp,2))
write.csv(row7b,"output/table_7_row_pooled_cpr.csv",row.names=FALSE)
cat("-> output/table_7_row_pooled_cpr.csv\n")
