###############################################################################
# 10_causality.R
# Dumitrescu-Hurlin (2012) panel non-causality test — Table 8
###############################################################################
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(plm); library(dplyr)
dir.create("output", showWarnings=FALSE)

panel  <- readRDS("data/panel_clean.rds")
pdata  <- readRDS("data/pdata.rds")
COUNTRIES <- unique(panel$country_iso)
N_    <- length(COUNTRIES); lag <- 1

# ── Dumitrescu-Hurlin test ─────────────────────────────────────────────────
dh_test <- function(pdata, y_var, x_var, lag_order=1) {
  W_i <- sapply(COUNTRIES, function(co) {
    d   <- panel[panel$country_iso==co,]; d<-d[order(d$year),]
    y   <- d[[y_var]]; x <- d[[x_var]]; T_<-length(y)
    yl  <- c(rep(NA,lag_order), head(y,-lag_order))
    xl  <- c(rep(NA,lag_order), head(x,-lag_order))
    df_ <- data.frame(y=y,yl=yl,xl=xl) |> na.omit()
    if (nrow(df_) < 5) return(NA_real_)
    mr <- lm(y~yl,      data=df_)  # restricted
    mu <- lm(y~yl+xl,   data=df_)  # unrestricted
    Rr <- sum(residuals(mr)^2); Ru <- sum(residuals(mu)^2)
    df_num <- lag_order; df_den <- max(1, nrow(df_)-2*lag_order-1)
    ((Rr-Ru)/df_num)/(Ru/df_den) * df_num  # W_i = lag * F_i
  })
  W_i    <- W_i[!is.na(W_i)]; N_ok <- length(W_i)
  W_bar  <- mean(W_i)
  Z_bar  <- sqrt(N_ok/(2*lag_order)) * (W_bar - lag_order)
  p_val  <- 2*(1-pnorm(abs(Z_bar)))
  data.frame(W_bar=round(W_bar,3), Z_bar=round(Z_bar,3),
             p_value=round(p_val,4),
             Decision=ifelse(p_val<0.05,"Reject H0","Not reject"))
}

pairs_ <- list(c("log_N","log_E"), c("log_E","log_N"),
               c("log_G","log_E"), c("log_Y","log_E"))
labels_ <- c("log N -> log E","log E -> log N",
             "log G -> log E","log Y -> log E")

rows <- lapply(seq_along(pairs_), function(i)
  cbind(data.frame(Direction=labels_[i]),
        dh_test(pdata, pairs_[[i]][2], pairs_[[i]][1])))
table8 <- do.call(rbind, rows)
print(table8); write.csv(table8,"output/table_8_causality.csv",row.names=FALSE)
cat("Table 8 -> output/table_8_causality.csv\n")
