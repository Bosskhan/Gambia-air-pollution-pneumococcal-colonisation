# Personal environmental exposure and *Streptococcus pneumoniae* colonisation in children 5–14 years in rural Gambia

R analysis code for the manuscript:

> Khan D, Kwambana-Adams B, Tientcheu P-E, Thomas S-A, Ndow PS, Kinteh M-A, Worwui A, Mohammed NI, Nicol MP, Antonio M.
> *Longitudinal personal environmental exposure and *Streptococcus pneumoniae* colonisation in a cohort of children 5–14 years old in rural Gambia.*

## Important: code only, data not included

**This repository contains the R analysis code only.** The study data (minute-level personal monitor data, REDCap metadata and microbiology results) are **not** included because they contain participant-level information from a cohort study.

**To request access to the data, please email Dam Khan: [Dam.khan@lshtm.ac.uk](mailto:Dam.khan@lshtm.ac.uk).**

Without the data, the script cannot be run end to end. It is shared for transparency so that the methods (quality control, exposure summaries, mixed-effects models, figures and tables) can be inspected and reproduced by approved data users.

## Study summary

- 134 children aged 5–14 years, Foni Bintang-Karanai District, The Gambia, followed over eight visits (2019–2022).
- 24-hour personal exposure to PM2.5, CO, temperature and relative humidity measured with the PATS+ monitor.
- *S. pneumoniae* colonisation determined by culture of combined nasopharyngeal/oropharyngeal swabs.
- Mixed-effects logistic regression (conventional and within-between decomposition models).

## Repository contents

| File | Description |
|---|---|
| `script.R` | Full analysis: QC, exposure summaries, tables, figures, models |
| `README.md` | This file |

## Input data required (not in this repository)

Place the following files in a `data/` folder. File names and formats are as assumed in the script; edit section 1 of `analysis_clean.R` if yours differ.

| File | Content |
|---|---|
| `Merged_PATS_V1_V4.xlsx` | Raw minute-level PATS+ data, visits 1–4 |
| `Merged_PATS_V5_V8.xlsx` | Raw minute-level PATS+ data, visits 5–8 |
| `Metadataairpollution.xlsx` | REDCap export (one row per child per visit): season, age, sex, household variables, *S. pneumoniae* result |

Figure 1 (map) also downloads the GADM Gambia administrative boundaries from the internet on first run.

## Software

- R (analyses were run in RStudio 2024.12.0.467)
- Packages: `tidyverse`, `lubridate`, `janitor`, `readxl`, `lme4`, `broom.mixed`, `car`, `patchwork`, `cowplot`, `scales`, `ggcorrplot`, `flextable`, `officer`, `sf`, `rnaturalearth`, `rnaturalearthdata`, `ggspatial`

Install with:

```r
install.packages(c("tidyverse", "lubridate", "janitor", "readxl", "lme4",
                   "broom.mixed", "car", "patchwork", "cowplot", "scales",
                   "ggcorrplot", "flextable", "officer", "sf",
                   "rnaturalearth", "rnaturalearthdata", "ggspatial"))
```

## How to run

1. Obtain the data (see above) and place the files in `data/`.
2. Open `analysis_clean.R` in R or RStudio and set the working directory to the repository root.
3. Run the script from top to bottom. Outputs are written to `output/`.

## Script structure

| Section | Content |
|---|---|
| 1 | Load raw data |
| 2 | PATS+ quality control (audit table saved to `output/qc_audit.csv`) |
| 3 | Child-visit 24-hour exposure summaries |
| 4 | Join metadata and define season and colonisation outcome |
| 5 | Baseline (visit 1) household variables, age, analysis dataset |
| 6 | Numbers quoted in the Results text |
| 7 | Tables 1–2 and Table S1 |
| 8 | Figures 1–3 and S1–S4 |
| 9 | Aim 3 models (Tables 3, S2–S5), Figures 4, 5 and S5 |

## Outputs

| Manuscript item | Output file |
|---|---|
| Figure 1 | `Figure1A_AfricaMap.png`, `Figure1B_GambiaMap.png` |
| Figure 2 | `Figure2_Monthly.png` |
| Figure 3 | `Figure3_Diurnal.png` |
| Figure 4 | `Figure4_Quartiles.png` |
| Figure 5 | `Figure5_WithinBetween.png` |
| Figures S1–S5 | `FigureS1_…` to `FigureS5_…` |
| Table 1 | `Table1_Baseline.docx` |
| Table 2 | `Table2_*.csv` |
| Tables S1–S5 | `TableS1_…` to `TableS5_…` (CSV) |
| Table 3 | `Table3_WithinBetween.csv` |

## Key methodological notes

- PM2.5 minute values at the lower detection limit (10 µg/m³) were substituted with LDL/√2 (7.07 µg/m³).
- A constant of 0.01 ppm was added to CO values before log transformation.
- Child-visits with less than 50% valid wear time (720 minutes) were excluded.
- Exposures were standardised for modelling; odds ratios are reported per IQR increase (PM2.5 and CO on the log scale).
- Within-between models decompose each exposure into a child-level mean (between) and the visit deviation from that mean (within).

## Contact

Dam Khan, MRC Unit The Gambia at LSHTM and London School of Hygiene & Tropical Medicine
Email: [Dam.khan@lshtm.ac.uk](mailto:Dam.khan@lshtm.ac.uk)

