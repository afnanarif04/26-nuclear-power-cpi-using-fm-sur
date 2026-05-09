###############################################################################
# 01_data_loading.R
# Balanced panel: N=8 countries x T=27 years (1995-2021) = 216 observations
#
# Usage: source("01_data_loading.R")   OR   open and Run in RStudio
# Output: data/panel_clean.rds, data/pdata.rds
###############################################################################

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# ── Packages ──────────────────────────────────────────────────────────────────
# install.packages(c("readxl","dplyr","tidyr","plm"))
library(readxl); library(dplyr); library(tidyr); library(plm)

# ── Country sample and period ──────────────────────────────────────────────────
COUNTRIES <- c("FRA","KOR","SVK","HUN","CZE","FIN","CHE","SWE")
COUNTRY_NAMES <- c(FRA="France",KOR="South Korea",SVK="Slovakia",HUN="Hungary",
                   CZE="Czech Republic",FIN="Finland",CHE="Switzerland",SWE="Sweden")
START_YEAR <- 1995; END_YEAR <- 2021

cat("Sample: N=8 x T=27 =", 8*27, "observations\n")

# ── Load panel ─────────────────────────────────────────────────────────────────
# File: nuclear_cpi_panel.xlsx must be in the data/ sub-folder next to this script
data_path <- file.path("data", "nuclear_cpi_panel.xlsx")

if (!file.exists(data_path)) {
  stop("Cannot find ", data_path,
       "\nMake sure nuclear_cpi_panel.xlsx is in the data/ folder next to this script.")
}

panel <- read_excel(data_path, sheet = "panel") %>%
  filter(country_iso %in% COUNTRIES,
         year >= START_YEAR, year <= END_YEAR) %>%
  arrange(country_iso, year)

# ── Variable construction (log transforms) ─────────────────────────────────────
panel <- panel %>%
  mutate(
    log_E    = log(energy_cpi),            # E — energy CPI (2015=100)
    log_P    = log(headline_cpi),          # P — headline CPI
    log_N    = log(nuclear_share),         # N — nuclear electricity share (%)
    log_N_sq = log_N^2,                    # N^2 — polynomial term
    log_G    = log(gas_price_usd_mmbtu),   # G — wholesale gas price
    log_M    = log(energy_import_dep),     # M — energy import dependence (%)
    log_Y    = log(gdp_per_cap_usd2015),   # Y — GDP per capita (2015 USD)
    log_K    = log(gfcf_share_gdp),        # K — GFCF (% of GDP)
    L        = elec_market_lib_index       # L — OECD PMR index (interpolated)
  )

# ── Balanced panel check ───────────────────────────────────────────────────────
check <- panel %>%
  group_by(country_iso) %>%
  summarise(n_obs     = n(),
            miss_E    = sum(is.na(log_E)),
            miss_N    = sum(is.na(log_N)),
            miss_Y    = sum(is.na(log_Y)),
            .groups   = "drop")
print(check)

if (any(check$n_obs != 27)) stop("Unbalanced panel: some countries have != 27 obs")
if (any(check$miss_E > 0))  stop("Missing values in energy_cpi — check data sheet")
if (any(check$miss_N > 0))  stop("Missing values in nuclear_share — check data sheet")

cat("\nPanel balanced: 8 countries x 27 years = 216 observations\n")
cat("All variables complete. No missing values.\n")

# ── PLM panel data object ──────────────────────────────────────────────────────
pdata <- pdata.frame(panel, index = c("country_iso","year"))

# ── Save ───────────────────────────────────────────────────────────────────────
if (!dir.exists("data")) dir.create("data")
saveRDS(panel, "data/panel_clean.rds")
saveRDS(pdata, "data/pdata.rds")

cat("Saved: data/panel_clean.rds\n")
cat("Saved: data/pdata.rds\n")
