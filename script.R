
##Air pollution and pneumococcal colonisation##########


library(tidyverse); library(lubridate); library(janitor); library(readxl)
library(lme4); library(broom.mixed); library(car)
library(patchwork); library(cowplot); library(scales); library(ggcorrplot)
library(flextable); library(officer)

LDL     <- 10                 # PATS+ PM lower detection limit (ug/m3)
LDL_SUB <- LDL / sqrt(2)      # 7.07 substitution value
WHO_PM  <- 15                 # WHO 2021 24 h PM2.5
WHO_CO  <- 3.5                # WHO 24 h CO (ppm)

metadata_raw <- read_excel("input/Metadataairpollution.xlsx")
pats_v1_v4 <- read_excel("input/Merged_PATS_V1V4.xlsx")
pats_v5_v8 <- read_excel("input/Merged_PATS_V5V8.xlsx")


# 1. Combine, clean the column names, and rename specific columns

pats_raw <- bind_rows(pats_v1_v4, pats_v5_v8) |>
  clean_names() |>
  rename(
    temp = deg_c_air, 
    rh   = percent_rh_air, 
    co   = co_ppm, 
    pm   = pm_estimate
  )


# 2. PATS+ QUALITY CONTROL  ->  pats_final (minute level)

audit <- tibble()
log_audit <- function(label, df) {
  audit <<- bind_rows(audit, tibble(
    step = label, n_rows = nrow(df),
    n_visits = n_distinct(paste(df$rfchildid, df$visitmp))))
}
log_audit("0. Raw", pats_raw)

# duplicate timestamps
qc1 <- pats_raw |> arrange(rfchildid, visitmp, date_time) |>
  distinct(rfchildid, visitmp, date_time, .keep_all = TRUE)
log_audit("1. De-duplicated", qc1)

# 2. device 'na' flag missing
qc2 <- qc1 |> filter(!is.na(na))
log_audit("2. na-flag removed", qc2)

# 3. visits with < 16 h clock duration
ok_dur <- qc2 |> group_by(rfchildid, visitmp) |>
  summarise(h = as.numeric(difftime(max(date_time), min(date_time), units = "hours")),
            .groups = "drop") |>
  filter(h >= 16) |> select(rfchildid, visitmp)
qc3 <- semi_join(qc2, ok_dur, by = c("rfchildid", "visitmp"))
log_audit("3. >=16 h duration", qc3)

# 4. trim to 25 h from first reading
qc4 <- qc3 |> group_by(rfchildid, visitmp) |>
  filter(as.numeric(difftime(date_time, min(date_time), units = "hours")) <= 25) |>
  ungroup()
log_audit("4. Trimmed to 25 h", qc4)

# 5a. device error flags (status != 0)
qc5a <- qc4 |> filter(status == 0)
log_audit("5a. status == 0", qc5a)

# 5b. manufacturer operating ranges
qc5b <- qc5a |> filter(pm >= LDL, pm <= 50000, co >= 0, co <= 500,
                       temp >= 0, temp <= 50, rh >= 0, rh <= 95)
log_audit("5b. In operating range", qc5b)

# 6. drop visits where >80% of rows (pre status filter) were pump error (status 8)
pump_bad <- qc4 |> group_by(rfchildid, visitmp) |>
  summarise(pct8 = mean(status == 8) * 100, .groups = "drop") |>
  filter(pct8 > 80) |> select(rfchildid, visitmp)
qc6 <- anti_join(qc5b, pump_bad, by = c("rfchildid", "visitmp"))
log_audit("6. Pump-error visits removed", qc6)

# 7. LDL substitution; fix doubled visit; require >= 720 valid minutes (50% of 24 h)
pats_final <- qc6 |>
  mutate(pm_adj      = if_else(pm == LDL, LDL_SUB, pm),
         hour_of_day = hour(date_time)) |>
  arrange(rfchildid, visitmp, date_time) |>
  group_by(rfchildid, visitmp) |>
  mutate(row_in_visit = row_number()) |>
  ungroup() |>
  # Si0140H1C1 V4 contains two recordings: keep first 1500 rows only
  filter(!(rfchildid == "Si0140H1C1" & visitmp == "V4" & row_in_visit > 1500)) |>
  select(-row_in_visit) |>
  add_count(rfchildid, visitmp, name = "n_valid") |>
  filter(n_valid >= 720) |> select(-n_valid)
log_audit("7. >=720 valid minutes (final)", pats_final)

print(audit)
write_csv(audit, "output/qc_audit.csv")


# 3. CHILD-VISIT 24 h SUMMARIES  ->  cv_summary

cv_summary <- pats_final |>
  group_by(rfchildid, visitmp) |>
  summarise(n_minutes = n(),
            pm_mean = mean(pm_adj), pm_max = max(pm_adj),
            co_mean = mean(co),     co_max = max(co),
            temp_mean = mean(temp), temp_max = max(temp),
            rh_mean = mean(rh),     rh_max = max(rh),
            .groups = "drop") |>
  mutate(pm_exceeds_WHO2021 = pm_mean > WHO_PM,
         co_exceeds_WHO     = co_mean > WHO_CO)
cv_summary


# 4. JOIN METADATA  ->  cv_joined

meta <- metadata_raw |>
  select(rfchildid = q3_record_id, visitmp = redcap_event_name,
         date_sampled, seasonsample_collected,
         sex = q12_sex, age_at_visit = Age_at_sample_collection,
         fuel_charcoal = q1_flsrc___1, fuel_wood = q1_flsrc___2,
         dist_road = q6_dstrd,
         burn_herbs = q10_brntitms___1, burn_incense = q10_brntitms___2,
         floor_type = q9_flrg, hh_income = q11_hhldinc,
         electricity = q5_elechom,
         smoker_compound = q_smoker, waste_disposal = q_waste_disposal,
         spn = Spn) |>
  filter(!is.na(visitmp))

View(cv_joined)

cv_joined <- cv_summary |>
  left_join(meta, by = c("rfchildid", "visitmp")) |>
  mutate(
    season = factor(case_when(seasonsample_collected == "Dry season" ~ "Dry",
                              seasonsample_collected == "Wet season" ~ "Wet"),
                    levels = c("Dry", "Wet")),
    spn_bin = case_when(spn == "Pos" ~ 1L, spn == "Neg" ~ 0L))


# season for minute-level data (used in diurnal summaries)
pats_final <- pats_final |>
  left_join(cv_joined |> select(rfchildid, visitmp, season),
            by = c("rfchildid", "visitmp"))

pats_final


# 5. BASELINE (V1) VARIABLES

baseline_final <- metadata_raw |>
  filter(redcap_event_name == "V1", q3_record_id %in% cv_joined$rfchildid) |>
  transmute(
    rfchildid = q3_record_id,
    fuel_type = factor(ifelse(q1_flsrc___1 == 1, "Charcoal", "Wood"),
                       levels = c("Wood", "Charcoal")),
    floor_cat = factor(case_when(
      q9_flrg == "Earth/Sand" ~ "Earthen",
      q9_flrg %in% c("Cement", "Vinyl or asphalt strips", "Other") ~ "Basic hard",
      q9_flrg %in% c("Ceramic tile", "Carpet") ~ "Tiled/finished"),
      levels = c("Earthen", "Basic hard", "Tiled/finished")),
    dist_road_cat = factor(case_when(
      q6_dstrd == 1 ~ "<50 m", q6_dstrd %in% 2:3 ~ "50-200 m",
      q6_dstrd %in% 4:5 ~ ">200 m"),
      levels = c("<50 m", "50-200 m", ">200 m")),
    hh_income_num = as.numeric(q11_hhldinc),
    hh_income_num = ifelse(hh_income_num < 1000, NA, hh_income_num),
    electricity_bin = factor(ifelse(q5_elechom == 1, "Yes", "No"), levels = c("No", "Yes")),
    sex_label = factor(q12_sex, levels = c("Male", "Female")),
    age_v1 = as.numeric(as.character(Age_at_sample_collection))) |>
  mutate(income_tertile = factor(ntile(hh_income_num, 3), 1:3,
                                 c("Low", "Medium", "High")))
stopifnot(nrow(baseline_final) == 132)


v1_dates <- metadata_raw |> filter(redcap_event_name == "V1") |>
  transmute(rfchildid = q3_record_id, date_v1 = as.Date(date_sampled))

aim2_data <- cv_joined |>
  left_join(baseline_final, by = "rfchildid") |>
  left_join(v1_dates, by = "rfchildid") |>
  mutate(
    age_recorded = suppressWarnings(as.numeric(age_at_visit)),
    age_calc     = age_v1 + as.numeric(difftime(as.Date(date_sampled), date_v1,
                                                units = "days")) / 365.25,
    age_at_visit_final = coalesce(age_recorded, age_calc),
    log_pm = log(pm_mean),
    log_co = log(co_mean + 0.01))         

aim3_data <- aim2_data |> filter(!is.na(spn_bin))
cat("Aim 3: ", nrow(aim3_data), "child-visits;",
    sum(aim3_data$spn_bin), "positive;",
    round(100 * mean(aim3_data$spn_bin), 1), "% prevalence\n")  


# 6. NUMBERS QUOTED IN THE RESULTS TEXT

geo_ci <- function(x) {                     
  lx <- log(x[x > 0]); se <- sd(lx) / sqrt(length(lx))
  c(GM = exp(mean(lx)), lo = exp(mean(lx) - 1.96 * se), hi = exp(mean(lx) + 1.96 * se))
}
cat("\n-- Overall exposure --\n")
print(round(geo_ci(cv_summary$pm_mean), 2))                 
print(round(geo_ci(cv_summary$co_mean + 0.01), 2))          
cat("PM AM (SD):", round(mean(cv_summary$pm_mean), 1), round(sd(cv_summary$pm_mean), 1), "\n")  
cat("CO AM (SD):", round(mean(cv_summary$co_mean), 1), round(sd(cv_summary$co_mean), 1), "\n")  
cat("PM <=10:", round(100 * mean(cv_summary$pm_mean <= 10), 1), "%\n")              
cat("CO <0.1 :", round(100 * mean(cv_summary$co_mean < 0.1), 1), "%\n")             
cat("CO ==0  :", round(100 * mean(cv_summary$co_mean < 0.01), 1), "%\n")            
cat("PM p95, max:", round(quantile(cv_summary$pm_mean, .95), 1), round(max(cv_summary$pm_mean), 1), "\n")  
cat("CO p95, max:", round(quantile(cv_summary$co_mean, .95), 1), round(max(cv_summary$co_mean), 1), "\n")  
cat("WHO exceedance PM / CO (%):", round(100 * mean(cv_summary$pm_exceeds_WHO2021), 1),
    round(100 * mean(cv_summary$co_exceeds_WHO), 1), "\n")                           
cat("\nWHO exceedance by season:\n")
cv_joined |> filter(!is.na(season)) |> group_by(season) |>
  summarise(n = n(), pm_pct = 100 * mean(pm_exceeds_WHO2021),
            co_pct = 100 * mean(co_exceeds_WHO)) |> print()                         

cat("\nVisits per child (AP data):\n")
vpc <- cv_joined |> count(rfchildid, name = "n_visits")
print(quantile(vpc$n_visits, c(.25, .5, .75))); print(range(vpc$n_visits))           
cat("Children with >=6 visits:", sum(vpc$n_visits >= 6),
    round(100 * mean(vpc$n_visits >= 6), 1), "%\n")                                  

month_order <- c("Nov","Dec","Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct")
cv_month <- cv_joined |>
  filter(!is.na(date_sampled)) |>
  mutate(date_sampled = as.Date(date_sampled),
         season = if_else(month(date_sampled) == 11, factor("Dry", levels = levels(season)), season),
         month_name = factor(month(date_sampled, label = TRUE, abbr = TRUE), levels = month_order)) |>
  filter(!is.na(month_name))
cv_month |> group_by(month_name) |>
  summarise(n = n(), pm = mean(pm_mean), co = mean(co_mean),
            temp = mean(temp_mean), rh = mean(rh_mean)) |> print(n = 12)

# Hourly means (text)
diurnal <- pats_final |> filter(!is.na(season)) |> group_by(hour_of_day) |>
  summarise(pm = mean(pm_adj), co = mean(co), temp = mean(temp), rh = mean(rh), .groups = "drop")
print(diurnal |> mutate(across(-hour_of_day, ~round(., 2))), n = 24)
diurnal


# 7. TABLES

# ---- Table 1: baseline characteristics (denominator = 132) ------------------
make_rows <- function(label, x, levels, labels = levels, show_missing = FALSE) {
  n <- length(x)
  f <- function(k) sprintf("%d (%.1f)", k, 100 * k / n)
  rows <- tibble(Characteristic = paste0("  ", labels),
                 Value = map_chr(levels, ~f(sum(x == .x, na.rm = TRUE))))
  if (show_missing)
    rows <- bind_rows(rows, tibble(Characteristic = "  Missing", Value = f(sum(is.na(x)))))
  bind_rows(tibble(Characteristic = label, Value = ""), rows)
}
age_q <- quantile(baseline_final$age_v1, c(.25, .5, .75), na.rm = TRUE)
tab1 <- bind_rows(
  make_rows("Sex", baseline_final$sex_label, c("Male", "Female")),
  tibble(Characteristic = c("Age at baseline, years", "  Median (IQR)"),
         Value = c("", sprintf("%.1f (%.0f-%.0f)", age_q[2], age_q[1], age_q[3]))),
  make_rows("Cooking fuel", baseline_final$fuel_type, c("Wood", "Charcoal")),
  make_rows("Floor type", baseline_final$floor_cat, levels(baseline_final$floor_cat),
            c("Earthen", "Cement/unfinished hard floor", "Finished floor (tile/carpet)")),
  make_rows("Distance to main road", baseline_final$dist_road_cat,
            levels(baseline_final$dist_road_cat), c("<50 m", "50-200 m", ">200 m")),
  make_rows("Household income tertiles", baseline_final$income_tertile,
            c("Low", "Medium", "High"), show_missing = TRUE),
  make_rows("Electricity in home", baseline_final$electricity_bin,
            c("Yes", "No"), show_missing = TRUE))
print(tab1, n = Inf)
save_as_docx(
  flextable(tab1) |> 
    set_header_labels(Value = "N = 132") |> 
    autofit(),
  path = "output/Table1_Baseline.docx"
)

# ---- Table 2: GM (GSD) of 24 h PM2.5 and CO by characteristic ---------------
# Built from VISIT-LEVEL covariates in cv_joined 
t2dat <- cv_joined |>
  filter(!is.na(season), !is.na(fuel_wood)) |>
  mutate(
    sex_label = factor(sex, levels = c("Male", "Female")),
    age_num   = as.numeric(age_at_visit),
    age_cat   = case_when(age_num < 9 ~ "5-8", age_num <= 11 ~ "9-11", age_num > 11 ~ "12-14"),
    fuel_type = ifelse(fuel_wood == 1, "Wood", "Charcoal"),
    dist_cat  = case_when(dist_road == 1 ~ "<50 m", dist_road %in% 2:3 ~ "50-200 m",
                          dist_road %in% 4:5 ~ ">200 m"),
    floor_cat = case_when(floor_type == "Earth/Sand" ~ "Earthen",
                          floor_type %in% c("Cement", "Vinyl or asphalt strips", "Other") ~ "Basic hard",
                          floor_type %in% c("Ceramic tile", "Carpet") ~ "Tiled/finished"),
    hh_income = as.numeric(hh_income),
    income_tertile = case_when(is.na(hh_income) ~ NA_character_,
                               hh_income <= quantile(hh_income, .33, na.rm = TRUE) ~ "Low",
                               hh_income <= quantile(hh_income, .67, na.rm = TRUE) ~ "Medium",
                               TRUE ~ "High"),
    electricity_cat = case_when(electricity == 1 ~ "Yes", electricity == 2 ~ "No"),
    waste_cat = case_when(waste_disposal == "Open burning in the household" ~ "Household burning",
                          waste_disposal == "Open burning in a designated area in thecommunity" ~ "Community burning"),
    smoker_cat = case_when(smoker_compound == 1 ~ "Yes", smoker_compound == 0 ~ "No"))

t2_vars <- c(Sex = "sex_label", `Age (years)` = "age_cat", Season = "season",
             `Cooking fuel` = "fuel_type", `Distance to road` = "dist_cat",
             `Incense burning` = "burn_incense", `Herbs burning` = "burn_herbs",
             `Floor type` = "floor_cat", `Household income` = "income_tertile",
             Electricity = "electricity_cat", `Waste disposal` = "waste_cat",
             `Smoker in compound` = "smoker_cat")


make_table2 <- function(offset_pm = 0, offset_co = 0.01) {
  imap_dfr(t2_vars, function(v, lab) {
    t2dat |> filter(!is.na(.data[[v]])) |>
      group_by(Category = as.character(.data[[v]])) |>
      summarise(N = n(),
                PM = sprintf("%.2f (%.2f)", exp(mean(log(pm_mean + offset_pm))),
                             exp(sd(log(pm_mean + offset_pm)))),
                CO = sprintf("%.2f (%.2f)", exp(mean(log(co_mean + offset_co))),
                             exp(sd(log(co_mean + offset_co)))),
                .groups = "drop") |>
      mutate(Characteristic = lab, .before = 1)
  })
}
table2_correct   <- make_table2(0, 0.01)
table2_published <- make_table2(0.5, 0.5)
write_csv(table2_correct,   "output/Table2_GM_GSD_correct_offsets.csv")
write_csv(table2_published, "output/Table2_GM_GSD_as_published_offset0.5.csv")

# ---- Table S1: child-visits monitored vs retained ---------------------------
tabS1 <- pats_raw |> distinct(rfchildid, visitmp) |> count(visitmp, name = "monitored") |>
  left_join(count(cv_summary, visitmp, name = "retained"), by = "visitmp") |>
  mutate(pct = round(100 * retained / monitored)) |>
  bind_rows(tibble(visitmp = "Total", monitored = NA, retained = NA, pct = NA)) |>
  mutate(monitored = ifelse(visitmp == "Total", sum(monitored, na.rm = TRUE), monitored))
print(tabS1)    # expect 726 monitored / 567 retained (78%)
write_csv(tabS1, "output/TableS1_Monitored_Retained.csv")


# 8. FIGURES 1-3 AND S1-S4 (descriptive)

season_cols <- c(Dry = "#E67E22", Wet = "#2980B9")

# ---- Figure 1: study map ----------------------------------------------------
library(sf); library(rnaturalearth); library(rnaturalearthdata); library(ggspatial)
africa  <- ne_countries(continent = "Africa", scale = "medium", returnclass = "sf")
gambia  <- ne_countries(country = "Gambia",   scale = "medium", returnclass = "sf")
senegal <- ne_countries(country = "Senegal",  scale = "medium", returnclass = "sf")

gadm_file <- file.path(dir_in, "gadm41_GMB_1.json")
if (!file.exists(gadm_file))
  download.file("https://geodata.ucdavis.edu/gadm/gadm4.1/json/gadm41_GMB_1.json", gadm_file)
gambia_adm1 <- st_read(gadm_file, quiet = TRUE) |>
  mutate(region_name = recode(NAME_1, LowerRiver = "Lower River Region",
                              MaccarthyIsland = "Central River Region",
                              NorthBank = "North Bank Region", UpperRiver = "Upper River Region",
                              Western = "Western Region"))
region_cols <- c("Lower River Region" = "#FF8D85", "Central River Region" = "#F9E79F",
                 "Upper River Region" = "#D7BDE2", "North Bank Region" = "#D5F5E3",
                 "Western Region" = "#A8D5A2")
banjul <- st_as_sf(data.frame(name = "Banjul", lon = -16.577, lat = 13.453),
                   coords = c("lon", "lat"), crs = 4326)

p_africa <- ggplot() +
  geom_sf(data = africa, fill = "#D5D8DC", colour = "white", linewidth = .25) +
  geom_sf(data = gambia, fill = "#E74C3C", colour = "#C0392B", linewidth = .6) +
  coord_sf(xlim = c(-20, 52), ylim = c(-36, 38), expand = FALSE) +
  theme_void() +
  theme(panel.background = element_rect(fill = "#AED6F1", colour = NA),
        panel.border = element_rect(colour = "#888888", fill = NA, linewidth = .6))
ggsave(out("Figure1A_AfricaMap.png"), p_africa, width = 12, height = 5.5, dpi = 300, bg = "white")

p_gambia <- ggplot() +
  geom_sf(data = senegal, fill = "#E8E8E3", colour = "white", linewidth = .3) +
  geom_sf(data = gambia_adm1, aes(fill = region_name), colour = "#444444", linewidth = .4) +
  scale_fill_manual(values = region_cols, na.value = "#DDDDDD") +
  geom_sf(data = banjul, shape = 8, size = 4, stroke = .9) +
  geom_sf_label(data = banjul, aes(label = name), nudge_y = .04, size = 2.8, fontface = "bold") +
  geom_sf_text(data = gambia_adm1, aes(label = region_name), size = 4, check_overlap = TRUE) +
  annotate("point", x = -16.193, y = 13.206, shape = 21, size = 3.5,
           colour = "#C0392B", fill = "#E74C3C") +
  annotate("label", x = -16.193, y = 13.155, label = "Foni Bintang-Karanai District",
           colour = "#C0392B", fontface = "bold", size = 4) +
  annotate("text", x = -15.5, y = c(12.82, 13.92), label = "Senegal",
           colour = "#888888", fontface = "italic", size = 3.5) +
  annotation_scale(location = "bl", width_hint = .25) +
  annotation_north_arrow(location = "br", which_north = "true",
                         style = north_arrow_orienteering(text_size = 10)) +
  coord_sf(xlim = c(-17.15, -13.55), ylim = c(12.72, 14.05), expand = FALSE) +
  labs(x = "Longitude", y = "Latitude") +
  theme_classic(base_size = 11) +
  theme(panel.background = element_rect(fill = "#AED6F1", colour = NA), legend.position = "none")
ggsave(out("Figure1B_GambiaMap.png"), p_gambia, width = 12, height = 5.5, dpi = 300, bg = "white")


# Figure 2: monthly exposure

month_n <- cv_month |> count(month_name) |>
  mutate(month_label = factor(paste0(month_name, "\n(n=", n, ")"),
                              levels = paste0(month_name, "\n(n=", n, ")")))
cv_month2 <- left_join(cv_month, month_n |> select(month_name, month_label), by = "month_name")
shading <- data.frame(xmin = c(.5, 7.5), xmax = c(7.5, 12.5), season = c("Dry", "Wet"))

monthly_loess <- function(var, ylab, title, show_x) {
  d <- cv_month2 |> mutate(x = as.numeric(month_label))
  m <- d |> group_by(x) |>
    summarise(mean = mean(.data[[var]]), se = sd(.data[[var]]) / sqrt(n()), .groups = "drop")
  ggplot() +
    geom_rect(data = shading, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = season),
              alpha = .15) +
    geom_smooth(data = d, aes(x, .data[[var]]), method = "loess", formula = y ~ x,
                colour = "black", fill = "grey60", alpha = .3, linewidth = 1) +
    geom_errorbar(data = m, aes(x, ymin = mean - se, ymax = mean + se), width = .2) +
    geom_point(data = m, aes(x, mean), size = 2) +
    scale_fill_manual(values = season_cols) +
    scale_x_continuous(breaks = 1:12, labels = levels(cv_month2$month_label)) +
    labs(y = ylab, title = title, x = NULL) +
    theme_minimal() +
    theme(legend.position = "none", panel.grid = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.text.x = if (show_x) element_text(angle = 45, hjust = 1, colour = "black") else element_blank(),
          plot.title = element_text(face = "bold", size = 13), plot.title.position = "panel")
}
fig2 <- monthly_loess("pm_mean", expression(PM[2.5]~(mu*g/m^3)), expression(PM[2.5]), FALSE) +
  monthly_loess("co_mean", "CO (ppm)", "CO", FALSE) +
  monthly_loess("temp_mean", "Temperature (\u00b0C)", "Temperature", TRUE) +
  monthly_loess("rh_mean", "Relative humidity (%)", "Relative humidity", TRUE) +
  plot_layout(ncol = 2)
fig2
ggsave("output/Figure2_Monthly.png", fig2, width = 11, height = 8, dpi = 300, bg = "white")



############################################
# ---- Figure 3: diurnal pattern ----------------------------------------------
# hourly_cv = one value per child-visit per hour (defined earlier)
hourly_pooled <- pats_final |> group_by(hour_of_day) |>
  summarise(pm_adj = mean(pm_adj), co = mean(co), temp = mean(temp), rh = mean(rh),
            .groups = "drop")

hourly_loess <- function(var, ylab, title, show_x) {
  ggplot(hourly_cv, aes(hour_of_day, .data[[var]])) +
    geom_smooth(method = "loess", formula = y ~ x, colour = "black", fill = "grey70",
                linewidth = 1.1, alpha = .4) +
    geom_point(data = hourly_pooled, size = 2, colour = "black") +
    scale_x_continuous(name = if (show_x) "Hour of day" else NULL,
                       breaks = seq(0, 21, 3), labels = paste0(seq(0, 21, 3), ":00")) +
    labs(y = ylab, title = title) +
    theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
          axis.text.x = if (show_x) element_text(colour = "black", size = 9) else element_blank(),
          plot.title = element_text(face = "bold", size = 13), plot.title.position = "panel")
}

fig3 <- hourly_loess("pm_adj", expression(PM[2.5]~(mu*g/m^3)), expression(PM[2.5]), FALSE) +
  hourly_loess("co",   "CO (ppm)",              "CO",                FALSE) +
  hourly_loess("temp", "Temperature (\u00b0C)", "Temperature",       TRUE)  +
  hourly_loess("rh",   "Relative humidity (%)", "Relative humidity", TRUE)  +
  plot_layout(ncol = 2)
fig3

ggsave("output/Figure3_Diurnal.png", fig3, width = 11, height = 8, dpi = 300, bg = "white")



#################################
# ---- Figure S1: distribution of 24 h means ----------------------------------
dist_panel <- function(var, xlab, fill, line, xmax, digits, title) {
  m <- mean(cv_joined[[var]]); md <- median(cv_joined[[var]])
  ggplot(cv_joined, aes(.data[[var]])) +
    geom_histogram(aes(y = after_stat(density)), bins = 40, fill = fill, colour = "white", alpha = .85) +
    geom_density(colour = line, linewidth = 1) +
    geom_vline(xintercept = m,  colour = "#E74C3C", linetype = "dashed", linewidth = .8) +
    geom_vline(xintercept = md, colour = "#27AE60", linetype = "dotted", linewidth = .8) +
    annotate("text", x = m + xmax / 50, y = Inf, vjust = 1.5, hjust = 0, size = 3.2,
             colour = "#E74C3C", label = paste0("Mean = ", round(m, digits))) +
    annotate("text", x = md + xmax / 50, y = Inf, vjust = 3.2, hjust = 0, size = 3.2,
             colour = "#27AE60", label = paste0("Median = ", round(md, digits + 1))) +
    scale_x_continuous(limits = c(0, xmax)) +
    labs(x = xlab, y = "Density", title = title) + theme_classic(base_size = 11)
}
figS1 <- dist_panel("pm_mean", "PM2.5 24-hour mean (\u00b5g/m\u00b3)", "#2980B9", "#1A5276", 100, 1, "A") +
  dist_panel("co_mean", "CO 24-hour mean (ppm)", "#E67E22", "#784212", 30, 2, "B")

ggsave("output/FigureS1_Distributions.png", figS1, width = 11, height = 8, dpi = 300, bg = "white")

# ---- Figure S2: per-child exposure heatmap ----------------------------------
visit_labels <- c("V1\n(02/19-04/19)", "V2\n(05/19-07/19)", "V3\n(08/19-12/19)", "V4\n(12/19-03/20)",
                  "V5\n(02/21-06/21)", "V6\n(06/21-08/21)", "V7\n(09/21-10/21)", "V8\n(12/21-01/22)")
heat <- cv_joined |>
  complete(rfchildid, visitmp = paste0("V", 1:8)) |>
  mutate(visitmp = factor(visitmp, levels = paste0("V", 1:8), labels = visit_labels),
         child_order = reorder(rfchildid, pm_mean, FUN = function(x) mean(x, na.rm = TRUE)))
heat_theme <- theme_classic(base_size = 10) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        plot.title = element_text(size = 14, face = "bold"),
        axis.title = element_text(size = 13, face = "bold"),
        axis.text.x = element_text(size = 13, face = "bold"))
heat_panel <- function(var, cols, vals, lim, nm, title)
  ggplot(heat, aes(visitmp, child_order, fill = .data[[var]])) +
  geom_tile(colour = "white", linewidth = .15) +
  scale_fill_gradientn(colours = cols, values = rescale(vals), na.value = "grey85",
                       name = nm, limits = lim, oob = squish) +
  labs(x = "Study visit", y = "Individual children (n = 132)", title = title) + heat_theme
figS2 <- heat_panel("pm_mean", c("#EBF5FB", "#AED6F1", "#2980B9", "#1A5276", "#17202A"),
                    c(0, 10, 20, 40, 94), c(7, 94), "PM2.5\n(\u00b5g/m\u00b3)", "A") /
  heat_panel("co_mean", c("#FEF9E7", "#FAD7A0", "#E67E22", "#A04000", "#4A235A"),
             c(0, .5, 2, 5, 26.8), c(0, 26.8), "CO\n(ppm)", "B")
ggsave("output/FigureS2_heatmap.png", figS2, width = 12, height = 14, dpi = 300, bg = "white")


# ---- Figure S3: individual 24 h means by season -----------------------------
cv_season <- cv_joined |> filter(!is.na(season))
points_panel <- function(var, who, ymax, step, ylab, title)
  ggplot(cv_season, aes(season, .data[[var]], colour = season)) +
  geom_jitter(width = .2, alpha = .45, size = 1.5) +
  geom_boxplot(width = .35, outlier.shape = NA, fill = NA, colour = "black", linewidth = .6) +
  geom_hline(yintercept = who, linetype = "dashed", colour = "#E74C3C", linewidth = .9) +
  scale_colour_manual(values = season_cols) +
  scale_y_continuous(limits = c(0, ymax), breaks = seq(0, ymax, step)) +
  labs(x = "Season", y = ylab, title = title) +
  theme_classic(base_size = 11) +
  theme(legend.position = "none", axis.text.x = element_text(size = 11, face = "bold"))
figS3 <- points_panel("pm_mean", WHO_PM, 100, 20, "PM2.5 24-hour mean (\u00b5g/m\u00b3)", "A") +
  points_panel("co_mean", WHO_CO, 30, 5, "CO 24-hour mean (ppm)", "B")

ggsave("output/FigureS3_Distributions.png", figS3, width = 11, height = 8, dpi = 300, bg = "white")



# ---- Figure S4: WHO exceedance by season ------------------------------------
exceed_panel <- function(flag) {
  d <- cv_season |>
    mutate(g = ifelse(.data[[flag]], "Exceeds WHO guideline", "Within WHO guideline")) |>
    count(season, g) |> group_by(season) |> mutate(pct = 100 * n / sum(n))
  ggplot(d, aes(season, pct, fill = g)) +
    geom_col(width = .55, colour = "white") +
    geom_text(aes(label = paste0(round(pct, 1), "%")), position = position_stack(vjust = .5),
              size = 3.5, fontface = "bold", colour = "white") +
    scale_fill_manual(values = c("Exceeds WHO guideline" = "#E74C3C",
                                 "Within WHO guideline" = "#AED6F1"), name = "") +
    scale_y_continuous(limits = c(0, 105), breaks = seq(0, 100, 25),
                       labels = paste0(seq(0, 100, 25), "%")) +
    labs(x = "", y = "Child-visits (%)") + theme_classic(base_size = 11) +
    theme(axis.text.x = element_text(size = 10, face = "bold"))
}
figS4 <- (exceed_panel("pm_exceeds_WHO2021") + exceed_panel("co_exceeds_WHO")) +
  plot_layout(guides = "collect") + plot_annotation(tag_levels = "A") &
  theme(legend.position = "bottom")

ggsave("output/FigureS4_WHOExceedance.png", figS4, width = 11, height = 8, dpi = 300, bg = "white")


# 9. MODELS  (Tables 3, S2, S3, S4, S5; Figures 4, 5, S5)

aim3_data <- aim3_data |>
  mutate(season_num = as.numeric(season == "Wet"),
         income_num = as.numeric(income_tertile)) |>
  mutate(log_pm_s = c(scale(log_pm)), log_co_s = c(scale(log_co)),
         temp_s = c(scale(temp_mean)), rh_s = c(scale(rh_mean)),
         age_s = c(scale(age_at_visit_final))) |>
  group_by(rfchildid) |>
  mutate(log_pm_between = mean(log_pm), log_co_between = mean(log_co),
         temp_between = mean(temp_mean), rh_between = mean(rh_mean),
         log_pm_within = log_pm - log_pm_between, log_co_within = log_co - log_co_between,
         temp_within = temp_mean - temp_between, rh_within = rh_mean - rh_between) |>
  ungroup() |>
  mutate(across(c(log_pm_between, log_co_between, temp_between, rh_between,
                  log_pm_within, log_co_within, temp_within, rh_within),
                ~c(scale(.x)), .names = "{.col}_s"))

fit <- function(f) glmer(as.formula(f), data = aim3_data, family = binomial,
                         control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))

exp_s   <- c(`PM2.5` = "log_pm_s", CO = "log_co_s", Temperature = "temp_s", `Rel. humidity` = "rh_s")
b_s     <- c(`PM2.5` = "log_pm_between_s", CO = "log_co_between_s",
             Temperature = "temp_between_s", `Rel. humidity` = "rh_between_s")
w_s     <- c(`PM2.5` = "log_pm_within_s", CO = "log_co_within_s",
             Temperature = "temp_within_s", `Rel. humidity` = "rh_within_s")
raw_var <- c(`PM2.5` = "log_pm", CO = "log_co", Temperature = "temp_mean", `Rel. humidity` = "rh_mean")
conf    <- "age_s + season + income_tertile"
re      <- "(1 | rfchildid)"

# conventional models
m_conv_crude  <- map(exp_s, ~fit(paste("spn_bin ~", .x, "+", re)))
m_conv_adj1   <- map(exp_s, ~fit(paste("spn_bin ~", .x, "+", conf, "+", re)))
m_conv_mutual <- fit(paste("spn_bin ~", paste(exp_s, collapse = " + "), "+", conf, "+", re))
# within-between models
m_wb_crude    <- map2(b_s, w_s, ~fit(paste("spn_bin ~", .x, "+", .y, "+", re)))
m_wb_adj1     <- map2(b_s, w_s, ~fit(paste("spn_bin ~", .x, "+", .y, "+", conf, "+", re)))
m_wb_mutual   <- fit(paste("spn_bin ~", paste(c(rbind(b_s, w_s)), collapse = " + "), "+", conf, "+", re))

# IQR scaling: coefficients are per SD -> multiply by IQR/SD
iqr_vals <- map_dbl(raw_var, ~IQR(aim3_data[[.x]], na.rm = TRUE))
mult     <- iqr_vals / map_dbl(raw_var, ~sd(aim3_data[[.x]], na.rm = TRUE))
print(round(iqr_vals, 2))      # article: 0.66, 3.22, 1.63, 26.8

or_iqr <- function(model, trm, m) {
  r <- tidy(model, effects = "fixed", conf.int = TRUE) |> filter(term == trm)
  sprintf("%.2f|%.2f-%.2f|%.3f", exp(r$estimate * m), exp(r$conf.low * m),
          exp(r$conf.high * m), r$p.value)
}
split_cell <- function(s) {p <- str_split_fixed(s, "\\|", 3); tibble(OR = p[, 1], CI = p[, 2], p = p[, 3])}

# ---- Table S3 (conventional) ------------------------------------------------
tabS3 <- imap_dfr(exp_s, function(t, nm) bind_cols(
  tibble(Exposure = nm),
  split_cell(or_iqr(m_conv_crude[[nm]], t, mult[[nm]])) |> rename_with(~paste0("Crude_", .)),
  split_cell(or_iqr(m_conv_adj1[[nm]],  t, mult[[nm]])) |> rename_with(~paste0("Adj1_", .)),
  split_cell(or_iqr(m_conv_mutual,      t, mult[[nm]])) |> rename_with(~paste0("Adj2_", .))))

print(tabS3, width = Inf)

write_csv(tabS3, "output/TableS3_Conventional.csv")

# ---- Table 3 (within-between) -----------------------------------------------
tab3 <- imap_dfr(exp_s, function(t, nm) map_dfr(c(Between = "between", Within = "within"), function(comp) {
  term <- if (comp == "between") b_s[[nm]] else w_s[[nm]]
  bind_cols(tibble(Exposure = nm, Component = comp),
            split_cell(or_iqr(m_wb_crude[[nm]], term, mult[[nm]])) |> rename_with(~paste0("Crude_", .)),
            split_cell(or_iqr(m_wb_adj1[[nm]],  term, mult[[nm]])) |> rename_with(~paste0("Adj1_", .)),
            split_cell(or_iqr(m_wb_mutual,      term, mult[[nm]])) |> rename_with(~paste0("Adj2_", .)))
}))
print(tab3, width = Inf); 

write_csv(tab3, "output/Table3_withinbetween.csv")

# ---- Tables S4 / S5 (covariates; age is per SD, not per year) ---------------
cov_terms <- c("age_s", "seasonWet", "income_tertileMedium", "income_tertileHigh")
cov_tbl <- function(model) tidy(model, effects = "fixed", conf.int = TRUE, exponentiate = TRUE) |>
  filter(term %in% cov_terms) |>
  transmute(term, OR = round(estimate, 2), CI = sprintf("%.2f-%.2f", conf.low, conf.high),
            p = round(p.value, 3))
tabS4 <- bind_rows(
  imap_dfr(m_conv_adj1, ~cov_tbl(.x) |> mutate(model = "Conventional", exposure_model = .y)),
  imap_dfr(m_wb_adj1,   ~cov_tbl(.x) |> mutate(model = "Within-between", exposure_model = .y)))
tabS5 <- bind_rows(cov_tbl(m_conv_mutual) |> mutate(model = "Conventional"),
                   cov_tbl(m_wb_mutual)   |> mutate(model = "Within-between"))

write_csv(tabS4, "output/TableS4_Covariates_Adj1.csv")
write_csv(tabS5, "output/TableS5_Covariates_Mutual.csv")


# ---- Table S2: VIF ----------------------------------------------------------
vif_tbl <- function(m) {
  v <- as.data.frame(vif(m)) |> rownames_to_column("term")
  if ("GVIF^(1/(2*Df))" %in% names(v))     # df-adjusted for multi-df terms (income); GVIF otherwise
    v <- v |> mutate(VIF = ifelse(Df == 1, GVIF, `GVIF^(1/(2*Df))`)) |> select(term, VIF)
  else names(v)[2] <- "VIF"
  mutate(v, VIF = round(VIF, 2))
}
tabS2 <- bind_rows(vif_tbl(m_conv_mutual) |> mutate(model = "Conventional"),
                   vif_tbl(m_wb_mutual)   |> mutate(model = "Within-between"))
print(tabS2)
write_csv(tabS2, "output/TableS2_VIF.csv")

# ---- Figure S5: Pearson correlation matrix ----------------------------------
var_order <- c("PM2.5", "CO", "Temperature", "Humidity", "Season", "Age", "Income")
cor_dat <- aim3_data |>
  transmute(PM2.5 = log_pm, CO = log_co, Temperature = temp_mean, Humidity = rh_mean,
            Season = season_num, Age = age_at_visit_final, Income = income_num) |>
  na.omit()
cor_mat <- cor(cor_dat, method = "pearson")[var_order, var_order]
print(round(cor_mat, 2))      # RH-season 0.80; PM-RH -0.34; CO-RH -0.22
figS5 <- ggcorrplot(cor_mat, type = "lower", lab = TRUE, lab_size = 3.8,
                    colors = c("#2166AC", "white", "#B2182B"), outline.color = "white",
                    ggtheme = theme_classic(base_size = 12)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 11),
        axis.text.y = element_text(face = "bold", size = 11))

ggsave("output/FigureS5_Correlation.png", figS5, width = 8, height = 7, dpi = 300, bg = "white")

# ---- Figure 4: prevalence by exposure quartile ------------------------------
quartile_panel <- function(var, label, units, fill, title) {
  d <- aim3_data |> mutate(q = factor(ntile(.data[[var]], 4), labels = paste0("Q", 1:4))) |>
    group_by(q) |> summarise(prev = 100 * mean(spn_bin), med = round(median(.data[[var]]), 2))
  ggplot(d, aes(q, prev)) +
    geom_col(fill = fill, alpha = .8, width = .6) +
    geom_text(aes(label = paste0(round(prev, 1), "%")), vjust = -.5, size = 4.2, fontface = "bold") +
    geom_text(aes(label = paste0("Med\n", med, " ", units)), vjust = 1.5, size = 3.2,
              colour = "white", fontface = "bold") +
    scale_y_continuous(limits = c(0, 65), labels = function(x) paste0(x, "%")) +
    labs(x = paste(label, "quartile"), y = "Prevalence (%)", title = title) +
    theme_classic(base_size = 10) +
    theme(plot.title = element_text(size = 16, face = "bold"),
          axis.title = element_text(size = 13, face = "bold"),
          axis.text = element_text(size = 12, colour = "black"))
}
fig4 <- (quartile_panel("pm_mean", "PM2.5", "\u00b5g/m\u00b3", "#2980B9", "PM2.5") +
           quartile_panel("co_mean", "CO", "ppm", "#E67E22", "CO")) /
  (quartile_panel("temp_mean", "Temperature", "\u00b0C", "#C0392B", "Temperature") +
     quartile_panel("rh_mean", "Relative humidity", "%", "#27AE60", "Relative humidity"))
ggsave("output/Figure4_Quartiles.png", fig4, width = 11, height = 9, dpi = 300, bg = "white")

# ---- Figure 5: within- (A) and between-child (B) associations ---------------
overall_prev <- mean(aim3_data$spn_bin)
cat("Overall prevalence for Fig 5 line:", round(100 * overall_prev, 1), "%\n")
child_summary <- aim3_data |> group_by(rfchildid) |>
  summarise(mean_pm = mean(log_pm), mean_co = mean(log_co), mean_temp = mean(temp_mean),
            mean_rh = mean(rh_mean), prop_col = mean(spn_bin), n_visits = n(), .groups = "drop")

fig5_theme <- theme_classic(base_size = 11) +
  theme(axis.title = element_text(face = "bold"), axis.text = element_text(colour = "black"),
        plot.title = element_text(face = "bold", size = 12))
within_panel <- function(var, xlab, colour, title, show_y) {
  d <- aim3_data |> mutate(bin = cut(.data[[var]], breaks = 10, include.lowest = TRUE)) |>
    group_by(bin) |> summarise(prop = mean(spn_bin), n = n(), mid = mean(.data[[var]]), .groups = "drop") |>
    filter(n >= 5)
  ggplot(d, aes(mid, prop)) +
    geom_hline(yintercept = overall_prev, linetype = "dotted", colour = "grey50") +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
    geom_point(aes(size = n), colour = colour, alpha = .75) +
    geom_smooth(method = "loess", formula = y ~ x, colour = "black", fill = "#C8C8C8", alpha = .4) +
    scale_y_continuous(labels = percent_format(), limits = c(0, .7)) +
    scale_size_continuous(name = "n visits", range = c(2, 7)) +
    labs(x = xlab, y = if (show_y) "Proportion colonised (%)" else NULL, title = title) +
    fig5_theme + if (!show_y) theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
}
between_panel <- function(var, xlab, colour, title, show_y) {
  ggplot(child_summary, aes(.data[[var]], prop_col)) +
    geom_hline(yintercept = overall_prev, linetype = "dotted", colour = "grey50") +
    geom_point(aes(size = n_visits), colour = colour, alpha = .65) +
    geom_smooth(method = "loess", formula = y ~ x, colour = "black", fill = "#C8C8C8", alpha = .4) +
    scale_y_continuous(labels = percent_format(), limits = c(0, 1)) +
    scale_size_continuous(name = "Number of\nvisits", range = c(1.5, 6)) +
    labs(x = xlab, y = if (show_y) "Proportion of visits colonised (%)" else NULL, title = title) +
    fig5_theme + if (!show_y) theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
}
fig5A <- (within_panel("log_pm_within", "Deviation from child's mean log PM2.5", "#4A4A4A", "PM2.5", TRUE) +
            within_panel("log_co_within", "Deviation from child's mean log CO", "#696969", "CO", FALSE)) /
  (within_panel("temp_within", "Deviation from child's mean temperature (\u00b0C)", "#5C5C5C", "Temperature", TRUE) +
     within_panel("rh_within", "Deviation from child's mean humidity (% RH)", "#7A7A7A", "Relative Humidity", FALSE)) +
  plot_layout(guides = "collect")
fig5B <- (between_panel("mean_pm", "Child-level mean log PM2.5", "#4A4A4A", "PM2.5", TRUE) +
            between_panel("mean_co", "Child-level mean log CO", "#696969", "CO", FALSE)) /
  (between_panel("mean_temp", "Child-level mean temperature (\u00b0C)", "#5C5C5C", "Temperature", TRUE) +
     between_panel("mean_rh", "Child-level mean relative humidity (%)", "#7A7A7A", "Relative Humidity", FALSE)) +
  plot_layout(guides = "collect")
fig5 <- wrap_elements(full = fig5A) / wrap_elements(full = fig5B) +
  plot_annotation(tag_levels = "A")
fig5
ggsave("output/Figure5_WithinBetween.png", fig5, width = 13, height = 18, dpi = 300, bg = "white")

writeLines(
  capture.output(sessionInfo()),
  "sessionInfo.txt"
)


##########**********************END*******************************####
