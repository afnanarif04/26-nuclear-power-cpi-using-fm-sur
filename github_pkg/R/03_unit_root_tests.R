###############################################################################
# 03_unit_root_tests.R
# CIPS panel unit root tests (Pesaran 2007) — Table 4
# Implementation: cross-sectionally augmented ADF (CADF) per country,
# then average -> CIPS statistic. No external package needed beyond urca.
###############################################################################
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(dplyr); library(urca)
dir.create("output", showWarnings=FALSE)

panel <- readRDS("data/panel_clean.rds")
COUNTRIES <- unique(panel$country_iso)
YEARS     <- sort(unique(panel$year))
T_        <- length(YEARS)

# ── CIPS: CADF test per country, then average ────────────────────────────────
cips_stat <- function(data_df, var, type="trend") {
  # Cross-section means at each t
  cs_mean <- tapply(data_df[[var]], data_df$year, mean, na.rm=TRUE)
  
  stats <- sapply(COUNTRIES, function(co) {
    d   <- data_df[data_df$country_iso==co, ]; d <- d[order(d$year),]
    y   <- d[[var]]; T_i <- length(y)
    if (any(is.na(y)) || T_i < 10) return(NA_real_)
    ym  <- cs_mean[as.character(d$year)]
    dy  <- c(NA, diff(y)); yl <- c(NA, head(y,-1))
    yml <- c(NA, head(ym,-1)); dym <- c(NA, diff(ym))
    df_ <- data.frame(dy=dy, yl=yl, ym=ym, yml=yml, dym=dym,
                      tr=seq_len(T_i)) |> na.omit()
    if (nrow(df_) < 6) return(NA_real_)
    fml <- if (type=="trend") dy~yl+ym+yml+dym+tr else dy~yl+ym+yml+dym
    m   <- tryCatch(lm(fml, data=df_), error=function(e) NULL)
    if (is.null(m)) return(NA_real_)
    co2 <- summary(m)$coefficients
    if ("yl" %in% rownames(co2)) co2["yl","t value"] else NA_real_
  })
  CIPS <- mean(stats, na.rm=TRUE)
  cv5  <- if (type=="trend") -2.42 else -2.07   # Pesaran (2007) 5% CV, T~27
  list(CIPS=round(CIPS,3), reject=(CIPS < cv5))
}

# Difference helper
first_diff <- function(df, var) {
  df %>% arrange(country_iso, year) %>% group_by(country_iso) %>%
    mutate(across(all_of(var), ~c(NA, diff(.x)), .names="d_{.col}")) %>%
    ungroup() %>% filter(!is.na(.data[[paste0("d_",var)]]))
}

ur_vars <- c("log_E","log_P","log_N","log_N_sq","log_G","log_M","log_Y","log_K","L")
ul      <- c("log E","log P","log N","(log N)\u00B2","log G","log M","log Y","log K","L")

rows <- lapply(seq_along(ur_vars), function(i) {
  v    <- ur_vars[i]
  lev  <- cips_stat(panel, v, "trend")
  diff_data <- first_diff(panel, v)
  dv   <- paste0("d_",v)
  # rename column for cips_stat
  diff_data[[v]] <- diff_data[[dv]]
  dif  <- cips_stat(diff_data, v, "drift")
  ord  <- if (!lev$reject & dif$reject) "I(1)" else
          if (lev$reject & dif$reject) "Mixed/I(0)" else "Higher"
  data.frame(Variable=ul[i], CIPS_level=lev$CIPS, Reject_level=lev$reject,
             CIPS_diff=dif$CIPS, Reject_diff=dif$reject, Order=ord)
})
table4 <- do.call(rbind, rows)
print(table4); write.csv(table4,"output/table_4_unit_root.csv",row.names=FALSE)
cat("Table 4 -> output/table_4_unit_root.csv\n")
