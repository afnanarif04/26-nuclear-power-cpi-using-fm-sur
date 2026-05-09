# Data Sources

The panel dataset (`nuclear_cpi_panel.xlsx`) covers **8 countries × 27 years (1995–2021)**.

| Variable | Column name | Source | Direct link |
|---|---|---|---|
| Energy consumer price index (2015 = 100) | `energy_cpi` | Eurostat HICP (CP04 — Housing, water, electricity, gas); ILO ILOSTAT for Korea and Switzerland | https://ec.europa.eu/eurostat/databrowser/view/PRC_HICP_AIND/default/table · https://ilostat.ilo.org/data/ |
| Headline CPI (2010 = 100) | `headline_cpi` | World Bank World Development Indicators (series FP.CPI.TOTL) | https://data.worldbank.org/indicator/FP.CPI.TOTL |
| Nuclear electricity share (% of total generation) | `nuclear_share` | Our World in Data — Share of electricity production from nuclear (sourced from Ember) | https://ourworldindata.org/grapher/share-electricity-nuclear |
| Wholesale natural gas price (USD per MMBtu) | `gas_price_usd_mmbtu` | World Bank Commodity Markets (Pink Sheet) — European natural gas (NGAS_EUR); LNG Japan (NGAS_JP) for Korea | https://www.worldbank.org/en/research/commodity-markets |
| Energy import dependence (% of total energy use) | `energy_import_dep` | World Bank World Development Indicators (series EG.IMP.CONS.ZS) | https://data.worldbank.org/indicator/EG.IMP.CONS.ZS |
| GDP per capita, constant 2015 USD | `gdp_per_cap_usd2015` | World Bank World Development Indicators (series NY.GDP.PCAP.KD) | https://data.worldbank.org/indicator/NY.GDP.PCAP.KD |
| Gross fixed capital formation (% of GDP) | `gfcf_share_gdp` | World Bank World Development Indicators (series NE.GDI.FTOT.ZS) | https://data.worldbank.org/indicator/NE.GDI.FTOT.ZS |
| Electricity market liberalisation index | `elec_market_lib_index` | OECD Product Market Regulation database — energy sector sub-index (benchmarks 1998, 2003, 2008, 2013, 2018, 2023; linearly interpolated to annual) | https://stats.oecd.org/index.aspx?queryid=96764 |

## Notes

- Energy CPI for Hungary 1995–2000 is chain-linked using OECD/FRED growth rates (series CPGREN01HUA659N) to extend coverage to 1995. Cells affected are marked in the `notes` sheet of the workbook.
- The electricity market liberalisation index is available at five-year benchmark years only. Annual values are obtained by linear interpolation between benchmarks; variation within each five-year interval is therefore mechanical.
- All continuous variables are transformed to natural logarithms in the R scripts. The nuclear share polynomial term `log_N_sq` = (log `nuclear_share`)² is constructed in `01_data_loading.R`.

## Country sample

| ISO code | Country |
|---|---|
| FRA | France |
| KOR | South Korea |
| SVK | Slovakia |
| HUN | Hungary |
| CZE | Czech Republic |
| FIN | Finland |
| CHE | Switzerland |
| SWE | Sweden |
