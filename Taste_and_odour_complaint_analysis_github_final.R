# ============================================================
# Taste and odour complaint analysis 
# Master thesis: Borgar Rudshaug, NTNU 2026
# ============================================================
# Sys.setenv(FROST_CLIENT_ID = "86043375-7353-4916-8234-23b84bfb30e6")
# ---- 0. Setup ----

library(tidyverse)
library(sf)
library(readxl)
library(readr)
library(lubridate)
library(vegan)
library(ggrepel)
library(httr)
library(jsonlite)
library(viridis)

config <- list(
  project_dir       = "D:/NTNU/4. semester masteroppgave/Kode",
  crs_wgs84         = 4326,
  crs_network       = 25832,
  crs_analysis      = 25833,
  event_start       = as.Date("2021-02-01"),
  event_end         = as.Date("2021-04-30"),
  benna_locations   = c("Benna rC%vann", "Benna rC%vann 2 (ikke i drift)", "Benna beh.vann"),
  jonsvatnet_locations = c("VIVA beh.vann", "Fortuna ventilkammer"),
  frost_client_id   = Sys.getenv("FROST_CLIENT_ID", unset = ""),
  weather_station   = "SN68230"
)

if (dir.exists(config$project_dir)) setwd(config$project_dir)

input_files <- list(
  complaints = "complaints.xlsx",
  basins     = "basins.xlsx",
  rawwater   = "RawwaterALLE.xlsx",
  benna_prod = "benna_production.xlsx",
  junctions  = "Storoutput.csv",
  pipes      = "pipes.xlsx",
  lookup     = "lookup.xlsx"
)

output_dir <- "Figures for the report"
dir.create(output_dir, showWarnings = FALSE)
dir.create(file.path(output_dir, "figures"), showWarnings = FALSE, recursive = TRUE)

# ---- Plot text size ----
# Increase figure text and legends by 20% relative to ggplot2 default (11 pt).
plot_base_size <- 11 * 1.20

dir.create(file.path(output_dir, "tables"),  showWarnings = FALSE, recursive = TRUE)

active_vars <- c(
  "colour_mean", "turbidity_mean", "water_age_mean",
  "air_temp_mean", "precipitation_sum", "kimtall_mean"
)

categorical_group_vars <- c(
  "nearest_basin_id", "dominant_source", "smell_group_dominant"
)

# ---- Service reservoir anonymisation (display only) ----
# Real reservoir names are kept out of every figure and table at the request
# of Trondheim kommune. The name-to-code mapping lives in a private file that
# is NOT committed to the public repository. basin_code is used only for
# display; all analysis and joins keep the original nearest_basin_id.
if (file.exists("basin_key_private.R")) {
  source("basin_key_private.R")
} else {
  stop("basin_key_private.R not found. It holds the private name-to-code key.")
}

code_basin <- function(x) {
  out <- unname(basin_key[as.character(x)])
  missing <- unique(as.character(x)[is.na(out)])
  if (length(missing)) stop("Unmapped basin name(s): ", paste(missing, collapse = ", "))
  out
}

# ---- 1. Helper functions ----

clean_text <- function(x) trimws(as.character(x))

safe_excel_date <- function(x) suppressWarnings(as.Date(as.numeric(x), origin = "1899-12-30"))

safe_numeric <- function(x) suppressWarnings(as.numeric(x))

parse_date_flexible <- function(x) {
  if (inherits(x, "Date")) return(x)
  if (inherits(x, "POSIXct") || inherits(x, "POSIXlt")) return(as.Date(x))
  out <- suppressWarnings(as.Date(x))
  if (all(is.na(out)) && is.numeric(x)) out <- safe_excel_date(x)
  if (all(is.na(out))) {
    out <- suppressWarnings(as.Date(lubridate::parse_date_time(
      as.character(x), orders = c("d.m.Y", "d/m/Y", "Y-m-d", "d-m-Y", "m/d/Y")
    )))
  }
  out
}

add_nearest_feature <- function(from_sf, to_sf, prefix) {
  if (nrow(from_sf) == 0 || nrow(to_sf) == 0) {
    return(tibble(
      !!paste0("nearest_", prefix, "_idx") := integer(),
      !!paste0("dist_to_", prefix, "_m")   := numeric()
    ))
  }
  d <- st_distance(from_sf, to_sf)
  tibble(
    !!paste0("nearest_", prefix, "_idx") := apply(d, 1, which.min),
    !!paste0("dist_to_", prefix, "_m")   := as.numeric(apply(d, 1, min))
  )
}

standardise_smell_category <- function(x) {
  case_when(
    x == "Earthy, clay"    ~ "earthy",
    x == "Earthy, fishy"   ~ "earthy",
    x == "Earthy, musty"   ~ "earthy",
    x == "Musty, Earthy"   ~ "musty",
    x == "Swampy, Earthy"  ~ "swampy",
    x == "Gasoline"        ~ "gasoline",
    x == "Mud"             ~ "earthy",
    x == "Mineral"         ~ "metallic",
    x == "Old toothbrush"  ~ "musty",
    x %in% c("No", "Yes") ~ NA_character_,
    TRUE ~ str_to_lower(as.character(x))
  )
}

build_mix_index <- function(colour, benna_median, jonsvatnet_median) {
  case_when(
    is.na(colour) | is.na(benna_median) | is.na(jonsvatnet_median) ~ NA_real_,
    benna_median == jonsvatnet_median ~ NA_real_,
    TRUE ~ (colour - benna_median) / (jonsvatnet_median - benna_median)
  ) |> pmin(1) |> pmax(0)
}

classify_source <- function(mix_index) {
  case_when(
    is.na(mix_index) ~ NA_character_,
    mix_index < 0.33  ~ "Benna",
    mix_index > 0.66  ~ "Jonsvatnet",
    TRUE              ~ "Mixed"
  ) |> factor(levels = c("Benna", "Mixed", "Jonsvatnet"))
}

combine_herlofsenloypa <- function(df) {
  if (!"nearest_basin_id" %in% names(df)) return(df)
  if (!"nearest_basin_id_raw" %in% names(df)) {
    df <- mutate(df, nearest_basin_id_raw = as.character(nearest_basin_id))
  }
  df |> mutate(
    nearest_basin_id = case_when(
      nearest_basin_id_raw %in% herlofsen_chambers ~ herlofsen_base,
      TRUE ~ as.character(nearest_basin_id_raw)
    )
  )
}

clean_group <- function(x) {
  x <- as.character(x)
  x[x == ""] <- NA_character_
  as.factor(x)
}

subset_dist <- function(d, idx) as.dist(as.matrix(d)[idx, idx])

safe_mean   <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
safe_median <- function(x) if (all(is.na(x))) NA_real_ else median(x, na.rm = TRUE)
safe_sd     <- function(x) if (sum(!is.na(x)) < 2) NA_real_ else sd(x, na.rm = TRUE)

first_non_missing <- function(x) { x <- x[!is.na(x)]; if (length(x) == 0) NA else x[1] }

dominant_non_missing <- function(x) {
  x <- as.character(x)
  x <- x[!is.na(x) & x != ""]
  if (length(x) == 0) return(NA_character_)
  names(sort(table(x), decreasing = TRUE))[1]
}

extract_adonis_row <- function(res, variable, type = "categorical") {
  tibble(
    variable = variable, type = type,
    Df = res$Df[1], R2 = res$R2[1],
    F_stat = res$F[1], p_value = res$`Pr(>F)`[1]
  )
}

safe_betadisper <- function(distance_object, group_vector, variable_name) {
  group_vector <- droplevels(clean_group(group_vector))
  keep <- !is.na(group_vector)
  if (sum(keep) < 10 || nlevels(group_vector[keep]) < 2) return(NULL)
  if (any(table(group_vector[keep]) < 2)) return(NULL)
  d_sub <- subset_dist(distance_object, keep)
  g_sub <- droplevels(group_vector[keep])
  bd <- vegan::betadisper(d_sub, g_sub)
  list(variable = variable_name, betadisper = bd, test = vegan::permutest(bd, permutations = 999))
}

# ---- 2. Data import ----

read_complaints <- function(path, crs_wgs84, crs_analysis) {
  df <- read_excel(path) |>
    rename(
      complaint_id = Sak, date_reported = Date, coord_raw = `Google koordinat`,
      location = Where, source = Source, taste_category = Taste,
      smell_category = Smell, particles = Particles, discolouration = Discolouration
    ) |>
    select(-any_of(c("lat", "long"))) |>
    separate(coord_raw, into = c("lat", "lon"), sep = ",", convert = TRUE) |>
    mutate(
      lat = safe_numeric(trimws(lat)), lon = safe_numeric(trimws(lon)),
      date_reported = parse_date_flexible(date_reported),
      complaint_id = as.character(complaint_id),
      smell_category = factor(smell_category), taste_category = factor(taste_category)
    ) |>
    filter(!is.na(lat), !is.na(lon))
  sf <- st_as_sf(df, coords = c("lon", "lat"), crs = crs_wgs84) |> st_transform(crs_analysis)
  list(df = df, sf = sf)
}

read_basins <- function(path, crs_network, crs_analysis) {
  df <- read_excel(path) |>
    select(Basin, lat, lon) |>
    rename(basin_id = Basin, easting = lat, northing = lon) |>
    mutate(easting = safe_numeric(easting), northing = safe_numeric(northing)) |>
    filter(!is.na(easting), !is.na(northing))
  sf <- st_as_sf(df, coords = c("easting", "northing"), crs = crs_network) |> st_transform(crs_analysis)
  list(df = df, sf = sf)
}

read_rawwater <- function(path) {
  read_excel(path, skip = 1) |>
    rename(
      source = Anlegg, date = Date, ecoli = `E-coli`, colour = Fargetall,
      kimtall = `Kimtall 22`, turbidity = Turbiditet, ph = pH, toc = TOC,
      conductivity = Konduktivitet, calcium = Kalsium, magnesium = Magnesium,
      alkalinity = Alkalitet, aluminium = Aluminium, ammonium = Ammonium
    ) |>
    fill(source, .direction = "down") |>
    filter(!is.na(date)) |>
    mutate(
      date = as.Date(date), colour = safe_numeric(colour), ph = safe_numeric(ph),
      toc = safe_numeric(toc), turbidity = safe_numeric(turbidity), kimtall = safe_numeric(kimtall)
    )
}

read_benna_production <- function(path) {
  read_excel(path, col_names = FALSE) |>
    set_names(c("date_serial", "volume_m3")) |>
    filter(!is.na(date_serial), date_serial != "Dag") |>
    mutate(date = safe_excel_date(date_serial), volume_m3 = safe_numeric(volume_m3)) |>
    select(date, volume_m3) |>
    filter(!is.na(date), !is.na(volume_m3))
}

read_junctions <- function(path, crs_network, crs_analysis) {
  raw <- read_delim(path, delim = ";", skip = 1, col_names = TRUE, trim_ws = TRUE, show_col_types = FALSE) |>
    rename_with(str_trim)
  df <- raw |>
    rename(
      junction_id = mw_Junction_MUID, x = mw_Junction_GeomX, y = mw_Junction_GeomY,
      result_name = `Result Name`, age_min = `Result Min`, age_max = `Result Max`
    ) |>
    filter(str_detect(result_name, regex("Water Quality", ignore_case = TRUE))) |>
    mutate(
      x = safe_numeric(x), y = safe_numeric(y),
      age_min = safe_numeric(age_min), age_max = safe_numeric(age_max),
      age_mean = (age_min + age_max) / 2
    ) |>
    filter(!is.na(x), !is.na(y), !is.na(age_mean))
  sf <- st_as_sf(df, coords = c("x", "y"), crs = crs_network) |> st_transform(crs_analysis)
  list(raw = raw, df = df, sf = sf)
}

read_pipes <- function(path, junctions_raw, crs_network, crs_analysis) {
  all_junctions <- junctions_raw |>
    rename(junction_id = mw_Junction_MUID, x = mw_Junction_GeomX, y = mw_Junction_GeomY) |>
    mutate(x = safe_numeric(x), y = safe_numeric(y),
           junction_id = iconv(junction_id, from = "UTF-8", to = "UTF-8", sub = "")) |>
    select(junction_id, x, y) |> filter(!is.na(x), !is.na(y)) |>
    distinct(junction_id, .keep_all = TRUE)
  
  df <- read_excel(path) |>
    transmute(
      pipe_id = ID, from_node = clean_text(`From node`), to_node = clean_text(`To Node`),
      length_m = safe_numeric(Length), diameter_mm = safe_numeric(Diameter),
      material = clean_text(Material), construction_date = safe_excel_date(`Construction date`)
    ) |>
    filter(!is.na(pipe_id), !is.na(material), material != "") |>
    mutate(
      pipe_age_years = as.numeric(difftime(Sys.Date(), construction_date, units = "days")) / 365.25,
      material_group = case_when(
        material %in% c("SJK", "SJG") ~ "Cast iron",
        material %in% c("PVC", "RDEL", "PERC") ~ "Plastic (PVC)",
        material %in% c("PE", "PE50", "PE80", "PE100", "PE125") ~ "Plastic (PE)",
        material == "BET" ~ "Concrete",
        material %in% c("ST", "STK", "GS") ~ "Steel",
        TRUE ~ "Other"
      )
    ) |>
    left_join(all_junctions, by = c("from_node" = "junction_id")) |>
    filter(!is.na(x), !is.na(y))
  sf <- st_as_sf(df, coords = c("x", "y"), crs = crs_network) |> st_transform(crs_analysis)
  list(df = df, sf = sf)
}

read_lookup_tables <- function(path) {
  category_lookup <- read_excel(path) |>
    rename(type = Type, category = Category, descriptor = Descriptor, possible_cause = `Possible cause`) |>
    mutate(across(everything(), str_trim))
  
  list(category_lookup = category_lookup, basin_name_lookup = basin_name_lookup)
}

# ---- 3. Derived datasets ----

build_monthly_colour <- function(rawwater_all, benna_locations, jonsvatnet_locations) {
  rawwater_all |>
    filter(source %in% c(benna_locations, jonsvatnet_locations)) |>
    mutate(
      source_group = case_when(
        source %in% benna_locations ~ "Benna",
        source %in% jonsvatnet_locations ~ "Jonsvatnet"
      ),
      year_month = floor_date(date, "month")
    ) |>
    group_by(source_group, year_month) |>
    summarise(median_colour = median(colour, na.rm = TRUE), .groups = "drop") |>
    pivot_wider(names_from = source_group, values_from = median_colour) |>
    rename(Benna_median = Benna, Jonsvatnet_median = Jonsvatnet)
}

build_basins_mix <- function(rawwater_all, monthly_colour) {
  basins_quality <- rawwater_all |>
    filter(str_detect(source, regex("hb", ignore_case = TRUE))) |>
    mutate(year_month = floor_date(date, "month"))
  
  basins_mix <- basins_quality |>
    left_join(monthly_colour, by = "year_month") |>
    mutate(
      mix_index = build_mix_index(colour, Benna_median, Jonsvatnet_median),
      dominant_source = classify_source(mix_index)
    )
  
  basins_mix_monthly <- basins_mix |>
    group_by(source, year_month) |>
    summarise(
      mix_index = mean(mix_index, na.rm = TRUE),
      colour_mean = mean(colour, na.rm = TRUE),
      ph_mean = mean(ph, na.rm = TRUE),
      turbidity_mean = mean(turbidity, na.rm = TRUE),
      toc_mean = mean(toc, na.rm = TRUE),
      kimtall_mean = mean(kimtall, na.rm = TRUE),
      dominant_source = classify_source(mean(mix_index, na.rm = TRUE)),
      .groups = "drop"
    )
  list(raw = basins_mix, monthly = basins_mix_monthly)
}

attach_spatial_neighbours <- function(complaints_df, complaints_sf, basins_df, basins_sf,
                                      pipes_df, pipes_sf, junctions_df, junctions_sf) {
  nearest_basin    <- add_nearest_feature(complaints_sf, basins_sf, "basin")
  nearest_pipe     <- add_nearest_feature(complaints_sf, pipes_sf, "pipe")
  nearest_junction <- add_nearest_feature(complaints_sf, junctions_sf, "junction")
  
  complaints_df |>
    bind_cols(nearest_basin) |>
    mutate(nearest_basin_id = basins_df$basin_id[nearest_basin_idx]) |>
    bind_cols(nearest_pipe) |>
    mutate(
      pipe_material    = pipes_df$material_group[nearest_pipe_idx],
      pipe_diameter_mm = pipes_df$diameter_mm[nearest_pipe_idx],
      pipe_age_years   = pipes_df$pipe_age_years[nearest_pipe_idx]
    ) |>
    bind_cols(nearest_junction) |>
    mutate(
      water_age_mean = junctions_df$age_mean[nearest_junction_idx],
      water_age_max  = junctions_df$age_max[nearest_junction_idx]
    )
}

attach_pressure_zones <- function(complaints_df, junctions_raw) {
  junction_zones <- junctions_raw |>
    rename(junction_id = mw_Junction_MUID, pressure_zone = mw_Junction_Pressure_ZoneID) |>
    mutate(
      junction_id   = iconv(junction_id, to = "UTF-8", sub = ""),
      pressure_zone = iconv(pressure_zone, to = "UTF-8", sub = ""),
      pressure_zone = trimws(pressure_zone),
      pressure_zone = na_if(pressure_zone, "\xa0"),
      pressure_zone = na_if(pressure_zone, "")
    ) |>
    select(junction_id, pressure_zone) |>
    distinct(junction_id, .keep_all = TRUE)
  
  complaints_df <- complaints_df |>
    mutate(pressure_zone = junction_zones$pressure_zone[nearest_junction_idx])
  list(complaints_df = complaints_df, junction_zones = junction_zones)
}

attach_lookup_and_time_series <- function(complaints_df, monthly_colour, basin_name_lookup,
                                          basins_mix_monthly, benna_production) {
  complaints_df |>
    mutate(year_month = floor_date(date_reported, "month")) |>
    left_join(monthly_colour, by = "year_month") |>
    left_join(basin_name_lookup, by = c("nearest_basin_id" = "basin_id")) |>
    left_join(basins_mix_monthly, by = c("source_name" = "source", "year_month"),
              relationship = "many-to-one") |>
    left_join(benna_production |> rename(benna_volume_m3 = volume_m3),
              by = c("date_reported" = "date"))
}

attach_smell_groups <- function(complaints_df, category_lookup) {
  smell_lookup <- category_lookup |>
    filter(type == "Smell") |>
    transmute(descriptor = str_to_lower(str_trim(descriptor)), category) |>
    distinct(descriptor, .keep_all = TRUE)
  complaints_df |>
    mutate(smell_normalised = standardise_smell_category(smell_category)) |>
    left_join(smell_lookup, by = c("smell_normalised" = "descriptor")) |>
    rename(smell_group = category)
}

attach_taste_groups <- function(complaints_df, category_lookup) {
  taste_lookup <- category_lookup |>
    filter(type == "Taste") |>
    transmute(descriptor = str_to_lower(str_trim(descriptor)), category) |>
    distinct(descriptor, .keep_all = TRUE)
  complaints_df |>
    mutate(taste_normalised = str_to_lower(str_trim(as.character(taste_category)))) |>
    left_join(taste_lookup, by = c("taste_normalised" = "descriptor")) |>
    rename(taste_group = category)
}

flag_event_period <- function(complaints_df, event_start, event_end, event_basins) {
  complaints_df |>
    mutate(event_2021 = date_reported >= event_start &
             date_reported <= event_end &
             nearest_basin_id %in% event_basins)
}

# ---- 4. Frost API ----

get_frost_daily_weather <- function(client_id, station, start_date, end_date) {
  response <- httr::GET(
    "https://frost.met.no/observations/v0.jsonld",
    query = list(
      sources = station,
      referencetime = paste0(start_date, "/", end_date),
      elements = "mean(air_temperature P1D),sum(precipitation_amount P1D)",
      timeoffsets = "default", levels = "default", qualities = "0,1,2,3,4"
    ),
    httr::authenticate(client_id, "")
  )
  httr::stop_for_status(response)
  json <- jsonlite::fromJSON(httr::content(response, as = "text", encoding = "UTF-8"), flatten = TRUE)
  json$data |>
    tidyr::unnest(observations) |>
    transmute(
      date = as.Date(referenceTime), element_id = elementId, value = value
    ) |>
    mutate(variable = case_when(
      element_id == "mean(air_temperature P1D)"  ~ "air_temp_daily",
      element_id == "sum(precipitation_amount P1D)" ~ "precipitation_daily",
      TRUE ~ element_id
    )) |>
    select(date, variable, value) |>
    pivot_wider(names_from = variable, values_from = value) |>
    arrange(date)
}

# ============================================================
# 5. Main preprocessing pipeline
# ============================================================

complaints    <- read_complaints(input_files$complaints, config$crs_wgs84, config$crs_analysis)
basins        <- read_basins(input_files$basins, config$crs_network, config$crs_analysis)
rawwater_all  <- read_rawwater(input_files$rawwater)
benna_production <- read_benna_production(input_files$benna_prod)
junctions     <- read_junctions(input_files$junctions, config$crs_network, config$crs_analysis)
pipes         <- read_pipes(input_files$pipes, junctions$raw, config$crs_network, config$crs_analysis)
lookups       <- read_lookup_tables(input_files$lookup)

monthly_colour <- build_monthly_colour(rawwater_all, config$benna_locations, config$jonsvatnet_locations)
basins_mix_obj <- build_basins_mix(rawwater_all, monthly_colour)

complaints_df <- attach_spatial_neighbours(
  complaints$df, complaints$sf, basins$df, basins$sf,
  pipes$df, pipes$sf, junctions$df, junctions$sf
)

complaints_df <- combine_herlofsenloypa(complaints_df)
complaints_df <- complaints_df |> mutate(basin_code = code_basin(nearest_basin_id))

complaints_df <- attach_lookup_and_time_series(
  complaints_df, monthly_colour, lookups$basin_name_lookup,
  basins_mix_obj$monthly, benna_production
)

pz_obj <- attach_pressure_zones(complaints_df, junctions$raw)
complaints_df  <- pz_obj$complaints_df

complaints_df <- complaints_df |>
  attach_smell_groups(lookups$category_lookup) |>
  attach_taste_groups(lookups$category_lookup) |>
  flag_event_period(config$event_start, config$event_end, event_basins) |>
  mutate(complaint_date = date_reported, year_month = floor_date(complaint_date, "month"))

write_csv(complaints_df, file.path(output_dir, "tables", "complaints_df_processed.csv"))

# ============================================================
# 6. Descriptive summaries
# ============================================================

complaints_baseline <- complaints_df |> filter(!event_2021)

monthly_complaints <- complaints_baseline |>
  count(year_month, name = "n_complaints") |> arrange(year_month)

basin_complaints <- complaints_baseline |>
  count(nearest_basin_id, name = "n_complaints", sort = TRUE) |>
  mutate(basin_code = code_basin(nearest_basin_id),
         pct = round(n_complaints / sum(n_complaints) * 100, 1))

source_complaints <- complaints_baseline |>
  filter(!is.na(dominant_source)) |>
  count(dominant_source, name = "n_complaints") |>
  mutate(pct = round(n_complaints / sum(n_complaints) * 100, 1))

smell_group_counts <- complaints_baseline |>
  count(smell_group, name = "n_complaints", sort = TRUE) |>
  mutate(pct = round(n_complaints / sum(n_complaints) * 100, 1))

taste_group_counts <- complaints_baseline |>
  count(taste_group, name = "n_complaints", sort = TRUE) |>
  mutate(pct = round(n_complaints / sum(n_complaints) * 100, 1))

# Complaint rate per junction (proxy for consumers served)
junction_basin_idx <- add_nearest_feature(junctions$sf, basins$sf, "basin")
junctions_per_basin <- tibble(
  nearest_basin_id = basins$df$basin_id[junction_basin_idx$nearest_basin_idx]
) |>
  mutate(nearest_basin_id = case_when(
    str_detect(nearest_basin_id, herlofsen_pattern) ~ herlofsen_base,
    TRUE ~ nearest_basin_id
  )) |>
  count(nearest_basin_id, name = "n_junctions")

complaint_rate <- complaints_baseline |>
  count(nearest_basin_id, name = "n_complaints") |>
  left_join(junctions_per_basin, by = "nearest_basin_id") |>
  mutate(basin_code = code_basin(nearest_basin_id),
         complaints_per_1000 = n_complaints / n_junctions * 1000) |>
  arrange(desc(complaints_per_1000))

basin_source <- complaints_baseline |>
  filter(!is.na(dominant_source)) |>
  count(nearest_basin_id, dominant_source) |>
  slice_max(n, by = nearest_basin_id, with_ties = FALSE) |>
  select(nearest_basin_id, basin_dominant_source = dominant_source)

rate_by_source <- complaint_rate |>
  left_join(basin_source, by = "nearest_basin_id") |>
  filter(!is.na(basin_dominant_source)) |>
  summarise(
    n_complaints = sum(n_complaints),
    n_junctions = sum(n_junctions),
    complaints_per_1000 = sum(n_complaints) / sum(n_junctions) * 1000,
    .by = basin_dominant_source
  )

# Save descriptive tables
write_csv(monthly_complaints,  file.path(output_dir, "tables", "monthly_complaints_baseline.csv"))
write_csv(basin_complaints,    file.path(output_dir, "tables", "complaints_by_basin_baseline.csv"))
write_csv(source_complaints,   file.path(output_dir, "tables", "complaints_by_dominant_source.csv"))
write_csv(smell_group_counts,  file.path(output_dir, "tables", "complaints_by_smell_group.csv"))
write_csv(taste_group_counts,  file.path(output_dir, "tables", "complaints_by_taste_group.csv"))
write_csv(complaint_rate,      file.path(output_dir, "tables", "complaint_rate_by_basin.csv"))
write_csv(rate_by_source,      file.path(output_dir, "tables", "complaint_rate_by_source.csv"))

# Descriptive figures
p_monthly <- ggplot(monthly_complaints, aes(x = year_month, y = n_complaints)) +
  geom_col(fill = "#3b82f6", alpha = 0.85) +
  scale_x_date(date_breaks = "3 months", date_labels = "%b %Y",
               guide = guide_axis(angle = 45)) +
  scale_y_continuous(breaks = seq(0, 50, by = 5)) +
  labs(x = NULL, y = "Number of complaints") +
  theme_bw(base_size = plot_base_size)
ggsave(file.path(output_dir, "figures", "monthly_complaints_baseline.png"), p_monthly, width = 9, height = 5, dpi = 600)



p_basin <- ggplot(basin_complaints, aes(x = reorder(basin_code, n_complaints), y = n_complaints)) +
  geom_col(fill = "#3b82f6", alpha = 0.85) + coord_flip() +
  scale_y_continuous(breaks = scales::breaks_pretty(n = 10)) +
  labs(x = "Service reservoir", y = "Number of complaints") +
  theme_bw(base_size = plot_base_size)
ggsave(file.path(output_dir, "figures", "complaints_by_basin_baseline.png"), p_basin, width = 9, height = 6, dpi = 600)



# ---- Monthly complaints coloured by dominant source ----

complaints_by_source_monthly <- complaints_baseline |>
  filter(!is.na(dominant_source)) |>
  count(year_month, dominant_source, name = "n_complaints")

p_monthly_source <- ggplot(complaints_by_source_monthly,
                           aes(x = year_month, y = n_complaints, fill = dominant_source)) +
  geom_col(alpha = 0.85) +
  scale_fill_manual(values = c(Benna = "#ef4444", Mixed = "#f59e0b", Jonsvatnet = "#3b82f6")) +
  scale_x_date(date_breaks = "3 months", date_labels = "%b %Y",
               guide = guide_axis(angle = 45)) +
  scale_y_continuous(breaks = seq(0, 50, by = 5)) +
  labs(x = NULL, y = "Number of complaints", fill = "Dominant source") +
  theme_bw(base_size = plot_base_size)
ggsave(file.path(output_dir, "figures", "monthly_complaints_by_source.png"),
       p_monthly_source, width = 9, height = 5, dpi = 600)

# ---- Smell classification coverage ----

n_baseline <- nrow(complaints_baseline)
n_classified <- sum(!is.na(complaints_baseline$smell_group))
n_unclassified <- n_baseline - n_classified

cat("\n=== Smell classification coverage ===\n")
cat("Total baseline complaints:", n_baseline, "\n")
cat("Classified smell group:", n_classified,
    "(", round(n_classified / n_baseline * 100, 1), "%)\n")
cat("No classified smell group:", n_unclassified,
    "(", round(n_unclassified / n_baseline * 100, 1), "%)\n")

# ---- Water age: network junctions vs complaint locations ----

junction_ages <- tibble(
  water_age = junctions$df$age_mean,
  group = "Network junctions"
)

complaint_ages <- complaints_baseline |>
  filter(!is.na(water_age_mean)) |>
  transmute(water_age = water_age_mean, group = "Complaint locations")

age_combined <- bind_rows(junction_ages, complaint_ages)

network_median <- median(junction_ages$water_age, na.rm = TRUE)
complaint_median <- median(complaint_ages$water_age, na.rm = TRUE)

p_basin <- ggplot(basin_complaints, aes(x = reorder(basin_code, n_complaints), y = n_complaints)) +
  geom_col(fill = "#3b82f6", alpha = 0.85) + coord_flip() +
  scale_y_continuous(breaks = scales::breaks_pretty(n = 10)) +
  labs(x = "Service reservoir", y = "Number of complaints") +
  theme_bw(base_size = plot_base_size)
ggsave(file.path(output_dir, "figures", "complaints_by_basin_baseline.png"), p_basin, width = 9, height = 6, dpi = 600)

cat("\n=== Water age comparison ===\n")
cat("Network junction median:", round(network_median, 1), "h\n")
cat("Complaint location median:", round(complaint_median, 1), "h\n")
cat("Complaints above network median:",
    round(mean(complaint_ages$water_age > network_median, na.rm = TRUE) * 100, 1), "%\n")
cat("Complaints above network 75th percentile:",
    round(mean(complaint_ages$water_age > quantile(junction_ages$water_age, 0.75, na.rm = TRUE),
               na.rm = TRUE) * 100, 1), "%\n")
wilcox_age <- wilcox.test(complaint_ages$water_age, junction_ages$water_age,
                          alternative = "greater")
cat("Wilcoxon W:", wilcox_age$statistic,
    " p:", format.pval(wilcox_age$p.value),
    " n_complaints:", nrow(complaint_ages),
    " n_junctions:", nrow(junction_ages), "\n")

# ---- Smell group by mix index and water age ----

smell_mix_df <- complaints_df |>
  filter(!is.na(smell_group), !is.na(mix_index)) |>
  mutate(dataset = ifelse(event_2021, "2021 event", "Baseline"))

p_smell_mix <- ggplot(smell_mix_df,
                      aes(x = mix_index, y = smell_group)) +
  geom_boxplot(alpha = 0.3, outlier.shape = NA) +
  geom_jitter(aes(shape = dataset, colour = dataset),
              height = 0.15, size = 2, alpha = 0.7) +
  scale_colour_manual(values = c("Baseline" = "grey40", "2021 event" = "#ef4444")) +
  scale_shape_manual(values = c("Baseline" = 16, "2021 event" = 17)) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.1)) +
  labs(x = "Mix index (0 = Benna, 1 = Jonsvatnet)",
       y = NULL, shape = NULL, colour = NULL) +
  theme_bw(base_size = plot_base_size) +
  theme(legend.position = "bottom")
ggsave(file.path(output_dir, "figures", "smell_group_by_mix_index.png"),
       p_smell_mix, width = 8, height = 5, dpi = 600)

network_median_wa <- median(junctions$df$age_mean, na.rm = TRUE)
smell_age_df <- complaints_df |>
  filter(!is.na(smell_group), !is.na(water_age_mean)) |>
  mutate(dataset = ifelse(event_2021, "2021 event", "Baseline"))

p_smell_age <- ggplot(smell_age_df,
                      aes(x = water_age_mean, y = smell_group)) +
  geom_boxplot(alpha = 0.3, outlier.shape = NA) +
  geom_jitter(aes(shape = dataset, colour = dataset),
              height = 0.15, size = 2, alpha = 0.7) +
  scale_colour_manual(values = c("Baseline" = "grey40", "2021 event" = "#ef4444")) +
  scale_shape_manual(values = c("Baseline" = 16, "2021 event" = 17)) +
  scale_x_continuous(breaks = seq(0, 100, by = 5)) +
  labs(x = "Modelled water age (hours)",
       y = NULL, shape = NULL, colour = NULL) +
  theme_bw(base_size = plot_base_size) +
  theme(legend.position = "bottom")

ggsave(file.path(output_dir, "figures", "smell_group_by_water_age.png"),
       p_smell_age, width = 8, height = 5, dpi = 600)
# ---- Shared map extent (covers every point, 5% margin) ----
bb <- sf::st_bbox(c(sf::st_geometry(basins$sf),
                    sf::st_geometry(complaints$sf)))
xpad <- 0.05 * as.numeric(bb["xmax"] - bb["xmin"])
ypad <- 0.05 * as.numeric(bb["ymax"] - bb["ymin"])
map_xlim <- c(bb["xmin"] - xpad, bb["xmax"] + xpad)
map_ylim <- c(bb["ymin"] - ypad, bb["ymax"] + ypad)

map_theme <- theme_void(base_size = plot_base_size) +
  theme(
    plot.background = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
    plot.title = element_text(
      hjust = 0,
      margin = margin(b = 8)
    ),
    plot.margin = margin(t = 10, r = 18, b = 10, l = 18),
    legend.position = "right",
    legend.title = element_text(size = plot_base_size),
    legend.text = element_text(size = plot_base_size * 0.9)
  )

# ---- Complaint locations and service reservoirs ----
p_map <- ggplot() +
  geom_sf(
    data = complaints$sf,
    aes(colour = "Complaint location"),
    size = 1.8,
    alpha = 0.55
  ) +
  geom_sf(
    data = basins$sf,
    aes(colour = "Service reservoir"),
    size = 3
  ) +
  scale_colour_manual(
    name = NULL,
    values = c(
      "Complaint location" = "#3b82f6",
      "Service reservoir" = "#ef4444"
    )
  ) +
  coord_sf(xlim = map_xlim, ylim = map_ylim, expand = FALSE, datum = NA) +
  labs(title = "Complaint locations and service reservoirs") +
  map_theme

ggsave(
  file.path(output_dir, "figures", "complaint_locations_and_basins.png"),
  p_map,
  width = 7,
  height = 9,
  dpi = 300,
  bg = "white"
)

p_age_map <- ggplot() +
  geom_sf(
    data = complaints_map_sf,
    aes(colour = water_age_mean),
    size = 2.4,
    alpha = 0.85
  ) +
  scale_colour_viridis_c(
    option = "inferno",
    name = "Water age\n(hours)"
  ) +
  coord_sf(xlim = map_xlim, ylim = map_ylim, expand = FALSE, datum = NA) +
  labs(title = "Complaint locations coloured by water age") +
  map_theme +
  theme(
    legend.key.height = unit(1.2, "cm"),
    legend.key.width = unit(0.35, "cm")
  )

ggsave(
  file.path(output_dir, "figures", "complaint_map_water_age.png"),
  p_age_map,
  width = 7,
  height = 9,
  dpi = 300,
  bg = "white"
)
# ---- Mix index vs complaints for three key basins ----

# focus_basins is defined in basin_key_private.R

mix_monthly <- basins_mix_obj$monthly |>
  left_join(lookups$basin_name_lookup, by = c("source" = "source_name")) |>
  filter(basin_id %in% focus_basins) |>
  mutate(basin_label = code_basin(basin_id))

complaints_monthly_basin <- complaints_df |>
  filter(nearest_basin_id %in% focus_basins) |>
  count(nearest_basin_id, year_month, name = "n_complaints") |>
  mutate(basin_label = code_basin(nearest_basin_id)) |>
  select(-nearest_basin_id)

p_mix <- ggplot() +
  geom_line(data = mix_monthly,
            aes(x = year_month, y = mix_index, colour = "Mix index"),
            linewidth = 0.6) +
  geom_col(data = complaints_monthly_basin,
           aes(x = year_month,
               y = n_complaints / max(n_complaints, na.rm = TRUE),
               fill = "Complaints"),
           alpha = 0.5, width = 25) +
  facet_wrap(~basin_label, ncol = 1, scales = "fixed") +
  scale_colour_manual(name = NULL, values = c("Mix index" = "#3b82f6")) +
  scale_fill_manual(name = NULL, values = c("Complaints" = "#ef4444")) +
  scale_x_date(date_breaks = "3 months", date_labels = "%b %Y",
               guide = guide_axis(angle = 45)) +
  scale_y_continuous(
    name = "Mix index (0 = Benna, 1 = Jonsvatnet)",
    sec.axis = sec_axis(~ . * max(complaints_monthly_basin$n_complaints, na.rm = TRUE),
                        name = "Complaints",
                        breaks = seq(0, 30, by = 5))
  ) +
  labs(x = NULL) +
  theme_bw(base_size = plot_base_size) +
  theme(strip.text = element_text(face = "bold"),
        legend.position = "bottom")
ggsave(file.path(output_dir, "figures", "mix_vs_complaints_combined.png"),
       p_mix, width = 9, height = 8, dpi = 600)
# ============================================================
# 7. Fetch and join Frost weather data
# ============================================================

frost_start <- min(complaints_df$complaint_date, na.rm = TRUE)
frost_end   <- max(complaints_df$complaint_date, na.rm = TRUE)

weather_daily <- get_frost_daily_weather(
  config$frost_client_id, config$weather_station, frost_start, frost_end
)

weather_monthly <- weather_daily |>
  mutate(year_month = floor_date(date, "month")) |>
  group_by(year_month) |>
  summarise(
    air_temp_mean      = safe_mean(air_temp_daily),
    precipitation_sum  = if (all(is.na(precipitation_daily))) NA_real_ else sum(precipitation_daily, na.rm = TRUE),
    precipitation_mean = safe_mean(precipitation_daily),
    wet_days           = sum(precipitation_daily > 0, na.rm = TRUE),
    heavy_rain_days    = sum(precipitation_daily >= 10, na.rm = TRUE),
    .groups = "drop"
  )

complaints_df <- complaints_df |>
  select(-any_of(c("air_temp_mean", "precipitation_sum", "precipitation_mean", "wet_days", "heavy_rain_days"))) |>
  left_join(weather_monthly, by = "year_month")

write_csv(weather_daily,   file.path(output_dir, "tables", "frost_weather_daily.csv"))
write_csv(weather_monthly, file.path(output_dir, "tables", "frost_weather_monthly.csv"))

# ============================================================
# 8. Basin-month aggregation
# ============================================================

data_pre <- complaints_df |> filter(!event_2021)

data_ord <- data_pre |>
  group_by(nearest_basin_id, year_month) |>
  summarise(
    n_complaints       = n(),
    complaint_date_min = min(complaint_date, na.rm = TRUE),
    complaint_date_max = max(complaint_date, na.rm = TRUE),
    colour_mean        = safe_mean(colour_mean),
    turbidity_mean     = safe_mean(turbidity_mean),
    water_age_mean     = safe_mean(water_age_mean),
    water_age_median   = safe_median(water_age_mean),
    water_age_sd       = safe_sd(water_age_mean),
    kimtall_mean       = safe_mean(kimtall_mean),
    air_temp_mean      = first_non_missing(air_temp_mean),
    precipitation_sum  = first_non_missing(precipitation_sum),
    precipitation_mean = first_non_missing(precipitation_mean),
    wet_days           = first_non_missing(wet_days),
    heavy_rain_days    = first_non_missing(heavy_rain_days),
    dominant_source        = dominant_non_missing(dominant_source),
    mix_index_mean         = safe_mean(mix_index),
    smell_group_dominant   = dominant_non_missing(smell_group),
    n_smell_group_nonmissing = sum(!is.na(smell_group)),
    taste_group_dominant   = dominant_non_missing(taste_group),
    n_taste_group_nonmissing = sum(!is.na(taste_group)),
    .groups = "drop"
  ) |>
  mutate(
    smell_group_dominant = clean_group(smell_group_dominant),
    taste_group_dominant = clean_group(taste_group_dominant),
    nearest_basin_id     = clean_group(nearest_basin_id),
    dominant_source      = clean_group(dominant_source)
  )

write_csv(data_ord, file.path(output_dir, "tables", "ordination_data_basin_month_before_complete_cases.csv"))

# Complete case filtering
data_ord <- data_ord |> filter(if_all(all_of(active_vars), ~ !is.na(.x)))

for (v in categorical_group_vars) {
  if (v %in% names(data_ord)) data_ord[[v]] <- clean_group(data_ord[[v]])
}

data_ord <- data_ord |> mutate(basin_code = code_basin(nearest_basin_id))



group_sizes <- map_dfr(
  categorical_group_vars[categorical_group_vars %in% names(data_ord)],
  function(v) {
    data_ord |>
      mutate(group = as.character(.data[[v]])) |>
      count(group, name = "n") |>
      mutate(
        label = v,
        included_in_PERMANOVA = !is.na(group)
      ) |>
      select(label, group, n, included_in_PERMANOVA)
  }
)

write_csv(group_sizes, file.path(output_dir, "tables", "group_sizes.csv"))

# ============================================================
# 9. Scaling, distance and PCA
# ============================================================

set.seed(42)

mat_raw    <- data_ord |> select(all_of(active_vars))
mat_scaled <- as.data.frame(scale(mat_raw))
colnames(mat_scaled) <- active_vars
dist_euc <- dist(mat_scaled, method = "euclidean")

pca_res     <- prcomp(mat_scaled, center = FALSE, scale. = FALSE)
pca_summary <- summary(pca_res)
pca_var     <- pca_summary$importance[2, 1:3] * 100

pca_scores <- as.data.frame(pca_res$x[, 1:3])
colnames(pca_scores) <- c("PC1", "PC2","PC3")
pca_scores <- bind_cols(data_ord, pca_scores)

loadings <- as.data.frame(pca_res$rotation[, 1:3])
loadings$variable <- rownames(loadings)

arrow_scale <- min(
  diff(range(pca_scores$PC1)) / diff(range(loadings$PC1)),
  diff(range(pca_scores$PC2)) / diff(range(loadings$PC2)),
  diff(range(pca_scores$PC3)) / diff(range(loadings$PC3))
) * 0.7


var_labels <- c(
  colour_mean       = "Colour",
  turbidity_mean    = "Turbidity",
  water_age_mean    = "Water age",
  air_temp_mean     = "Air temperature",
  precipitation_sum = "Precipitation",
  kimtall_mean      = "HPC"
)

loadings_for_plot <- loadings |>
  mutate(PC1s  = PC1 * arrow_scale,
         PC2s  = PC2 * arrow_scale,
         label = unname(var_labels[variable]))



# PCA biplot coloured by nearest basin
p_bi <- ggplot() +
  geom_point(data = pca_scores, aes(x = PC1, y = PC2, colour = basin_code),
             size = 2, alpha = 0.7) +
  geom_segment(data = loadings_for_plot,
               aes(x = 0, y = 0, xend = PC1s, yend = PC2s),
               arrow = arrow(length = unit(0.2, "cm")), colour = "grey30", linewidth = 0.5) +
  geom_text_repel(data = loadings_for_plot,
                  aes(x = PC1s, y = PC2s, label = label), size = 3.2, colour = "grey20") +
  labs(title = "PCA biplot by service reservoir",
       x = paste0("PC1 (", round(pca_var[1], 1), "%)"),
       y = paste0("PC2 (", round(pca_var[2], 1), "%)"),
       colour = "Service reservoir") +
  theme_bw(base_size = plot_base_size)
ggsave(file.path(output_dir, "figures", "pca_biplot.png"), p_bi, width = 9, height = 6, dpi = 300)

# PCA coloured by dominant source (with convex hulls)
source_df <- pca_scores |> filter(!is.na(dominant_source))
hull_df <- source_df |>
  group_by(dominant_source) |> filter(n() >= 3) |>
  slice(chull(PC1, PC2)) |> ungroup()

p_source <- ggplot(source_df, aes(x = PC1, y = PC2, colour = dominant_source, fill = dominant_source)) +
  geom_polygon(data = hull_df, alpha = 0.10, linetype = "dashed", linewidth = 0.6, show.legend = FALSE) +
  geom_point(size = 2, alpha = 0.7) +
  labs(title = "PCA by dominant source",
       x = paste0("PC1 (", round(pca_var[1], 1), "%)"),
       y = paste0("PC2 (", round(pca_var[2], 1), "%)"),
       colour = "Dominant source") +
  guides(fill = "none") +
  theme_bw(base_size = plot_base_size)
ggsave(file.path(output_dir, "figures", "pca_dominant_source.png"), p_source, width = 9, height = 6, dpi = 300)

# ============================================================
# 9b. NMDS (supplementary ordination check)
# ============================================================

set.seed(42)
nmds_res <- vegan::metaMDS(dist_euc, k = 2, trymax = 100, trace = FALSE)

cat("NMDS stress (2D):", round(nmds_res$stress, 3), "\n")
writeLines(
  paste0("NMDS stress (2D, Euclidean, k=2): ", round(nmds_res$stress, 4)),
  file.path(output_dir, "tables", "nmds_stress.txt")
)

nmds_scores <- as.data.frame(vegan::scores(nmds_res, display = "sites"))
colnames(nmds_scores) <- c("NMDS1", "NMDS2")
nmds_scores <- bind_cols(data_ord, nmds_scores)
nmds_subtitle <- paste0("Stress = ", round(nmds_res$stress, 3))

# NMDS coloured by nearest basin
p_nmds_basin <- ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2, colour = basin_code)) +
  geom_point(size = 2, alpha = 0.7) +
  labs(title = "NMDS by service reservoir", subtitle = nmds_subtitle, colour = "Service reservoir") +
  theme_bw(base_size = plot_base_size)
ggsave(file.path(output_dir, "figures", "nmds_nearest_basin_id.png"),
       p_nmds_basin, width = 9, height = 6, dpi = 300)

# NMDS coloured by dominant source (with convex hulls)
nmds_source_df <- nmds_scores |> filter(!is.na(dominant_source))
nmds_hull_df <- nmds_source_df |>
  group_by(dominant_source) |> filter(n() >= 3) |>
  slice(chull(NMDS1, NMDS2)) |> ungroup()

p_nmds_source <- ggplot(nmds_source_df,
                        aes(x = NMDS1, y = NMDS2, colour = dominant_source, fill = dominant_source)) +
  geom_polygon(data = nmds_hull_df, alpha = 0.10, linetype = "dashed",
               linewidth = 0.6, show.legend = FALSE) +
  geom_point(size = 2, alpha = 0.7) +
  labs(title = "NMDS by dominant source", subtitle = nmds_subtitle, colour = "Dominant source") +
  guides(fill = "none") +
  theme_bw(base_size = plot_base_size)
ggsave(file.path(output_dir, "figures", "nmds_dominant_source.png"),
       p_nmds_source, width = 9, height = 6, dpi = 300)

# NMDS coloured by water age
p_nmds_age <- ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2, colour = water_age_mean)) +
  geom_point(size = 2, alpha = 0.85) +
  scale_colour_viridis_c(option = "inferno", name = "Water age\n(hours)") +
  labs(title = "NMDS by water age", subtitle = nmds_subtitle) +
  theme_bw(base_size = plot_base_size)
ggsave(file.path(output_dir, "figures", "nmds_water_age_mean.png"),
       p_nmds_age, width = 9, height = 6, dpi = 300)

# ============================================================
# 10. PERMANOVA and betadisper
# ============================================================

adonis_results    <- tibble()
betadisper_results <- list()

for (v in categorical_group_vars[categorical_group_vars %in% names(data_ord)]) {
  idx       <- !is.na(data_ord[[v]])
  group_vec <- droplevels(data_ord[[v]][idx])
  
  if (sum(idx) >= 10 && nlevels(group_vec) >= 2) {
    d_sub    <- subset_dist(dist_euc, idx)
    meta_sub <- tibble(group = group_vec)
    
    set.seed(42)
    res <- vegan::adonis2(d_sub ~ group, data = meta_sub, permutations = 999, by = "margin")
    adonis_results <- bind_rows(adonis_results, extract_adonis_row(res, variable = v))
    
    bd <- safe_betadisper(dist_euc, data_ord[[v]], v)
    if (!is.null(bd)) betadisper_results[[v]] <- bd
  }
}

# Combined model: basin + source
comb_idx <- complete.cases(data_ord[, c("nearest_basin_id", "dominant_source"), drop = FALSE])
if (sum(comb_idx) >= 10) {
  meta_comb <- data_ord[comb_idx, c("nearest_basin_id", "dominant_source"), drop = FALSE] |>
    mutate(across(everything(), clean_group))
  if (all(map_int(meta_comb, nlevels) >= 2)) {
    set.seed(42)
    adonis_combined <- vegan::adonis2(
      subset_dist(dist_euc, comb_idx) ~ nearest_basin_id + dominant_source,
      data = meta_comb, permutations = 999, by = "margin"
    )
  }
}

adonis_results <- adonis_results |> arrange(p_value, desc(R2))

# ============================================================
# 11. Spearman correlations
# ============================================================

cor_mat <- mat_raw |> cor(use = "pairwise.complete.obs", method = "spearman")

pairs_df <- expand.grid(var1 = active_vars, var2 = active_vars, stringsAsFactors = FALSE) |>
  filter(var1 < var2) |>
  rowwise() |>
  mutate(
    n_complete = sum(complete.cases(data_ord[[var1]], data_ord[[var2]])),
    rho = {
      cc <- complete.cases(data_ord[[var1]], data_ord[[var2]])
      if (sum(cc) >= 5) cor(data_ord[[var1]][cc], data_ord[[var2]][cc], method = "spearman") else NA_real_
    },
    p_value = {
      cc <- complete.cases(data_ord[[var1]], data_ord[[var2]])
      if (sum(cc) >= 5) cor.test(data_ord[[var1]][cc], data_ord[[var2]][cc],
                                 method = "spearman", exact = FALSE)$p.value else NA_real_
    }
  ) |>
  ungroup() |>
  arrange(p_value)

# ============================================================
# 12. Export all tables
# ============================================================

rationale <- tribble(
  ~variable,           ~description,                  ~rationale,
  "colour_mean",       "Colour",                      "Source water origin and NOM indicator.",
  "turbidity_mean",    "Turbidity",                   "Particles, treatment performance and resuspension indicator.",
  "water_age_mean",    "Mean modelled water age",     "Residence time proxy from hydraulic model.",
  "air_temp_mean",     "Mean air temperature",        "Seasonal variation from Frost API.",
  "precipitation_sum", "Monthly precipitation sum",   "Runoff and source water quality proxy from Frost API.",
  "kimtall_mean",      "Heterotrophic plate count",   "Microbial activity and biological stability indicator."
)

write_csv(rationale,       file.path(output_dir, "tables", "active_variable_rationale.csv"))
write_csv(loadings,        file.path(output_dir, "tables", "pca_loadings.csv"))
write_csv(adonis_results,  file.path(output_dir, "tables", "adonis2_results.csv"))
write_csv(pairs_df,        file.path(output_dir, "tables", "spearman_active_vars.csv"))
write_csv(data_ord,        file.path(output_dir, "tables", "ordination_data_basin_month_complete_cases.csv"))

write_csv(
  as.data.frame(cor_mat) |> rownames_to_column("variable"),
  file.path(output_dir, "tables", "spearman_correlation_matrix.csv")
)

if (exists("adonis_combined") && !is.null(adonis_combined)) {
  write_csv(
    as.data.frame(adonis_combined) |> rownames_to_column("term"),
    file.path(output_dir, "tables", "adonis2_combined_basin_source.csv")
  )
}

if (length(betadisper_results) > 0) {
  betadisper_summary <- map_dfr(names(betadisper_results), function(v) {
    tab <- betadisper_results[[v]]$test$tab
    tibble(variable = v, Df = tab$Df[1], SumSq = tab$`Sum Sq`[1],
           MeanSq = tab$`Mean Sq`[1], F_stat = tab$F[1], p_value = tab$`Pr(>F)`[1])
  })
  write_csv(betadisper_summary, file.path(output_dir, "tables", "betadisper_results.csv"))
}