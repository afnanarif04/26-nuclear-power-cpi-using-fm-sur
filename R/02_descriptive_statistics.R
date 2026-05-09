###############################################################################
# 02_descriptive_statistics.R
# Descriptive statistics — Table 2
###############################################################################
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(dplyr); library(tidyr); library(ggplot2)
dir.create("output", showWarnings=FALSE)

panel <- readRDS("data/panel_clean.rds")

desc_vars <- c("log_E","log_P","log_N","log_N_sq","log_G","log_M","log_Y","log_K","L")
labels    <- c("log E","log P","log N","(log N)\u00B2","log G","log M","log Y","log K","L")

sk <- function(x) { x<-x[!is.na(x)]; n<-length(x); m<-mean(x)
  (sum((x-m)^3)/n) / (sum((x-m)^2)/n)^1.5 }
kt <- function(x) { x<-x[!is.na(x)]; n<-length(x); m<-mean(x)
  (sum((x-m)^4)/n) / (sum((x-m)^2)/n)^2 }

table2 <- do.call(rbind, lapply(seq_along(desc_vars), function(i) {
  x <- panel[[desc_vars[i]]]
  data.frame(Variable=labels[i], Mean=round(mean(x,na.rm=TRUE),4),
             SD=round(sd(x,na.rm=TRUE),4), Min=round(min(x,na.rm=TRUE),4),
             Max=round(max(x,na.rm=TRUE),4), Skewness=round(sk(x),4),
             Kurtosis=round(kt(x),4))
}))
print(table2); write.csv(table2,"output/table_2_descriptive.csv",row.names=FALSE)

country_sum <- panel %>% group_by(country_iso) %>%
  summarise(nuc_mean=round(mean(nuclear_share,na.rm=TRUE),2),
            nuc_min=round(min(nuclear_share,na.rm=TRUE),2),
            nuc_max=round(max(nuclear_share,na.rm=TRUE),2),.groups="drop")
print(country_sum); write.csv(country_sum,"output/country_summary.csv",row.names=FALSE)
cat("Table 2 -> output/table_2_descriptive.csv\n")
