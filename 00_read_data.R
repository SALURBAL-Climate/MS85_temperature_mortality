library(tidyverse)
library(here)
library(glue)

# Read metadata describing cities -----------------------------------------

# Paths to read locations
path_overall     <- "//files.drexel.edu/encrypted/SOPH/UHC/Projects/SALURBAL HEAT Study/"
path_metadata    <- glue("{path_overall}/Data/Files/Derived Data/Meta Data/")
path_climatedata <- glue("{path_overall}/Data/Files/Climate Zone/Processed Data/")

# Path to updated data mortality/temperature data
path_all_cause <- "//files.drexel.edu/colleges/SOPH/Shared/UHC/Projects/Wellcome_Trust/Manuscripts/MS85_Update/2025_08_20/data/clean/"

# Read metadata
metadata0 <- haven::read_sas(glue("{path_metadata}/city_level_temp.sas7bdat"))
hasmortality <- haven::read_sas(glue("{path_metadata}/hasmortality_indic.sas7bdat"))
clusters <- read_csv(glue("{path_overall}/Analysis/Clusters/Results/ECDFClusters.csv"))
climate_labels <- read_csv(glue("{path_climatedata}/climate_labels2.csv"))
inAICcalculations <- read_csv(glue("{path_metadata}/inAICcalculations.csv"), col_select = c(nsalid1, inAIC))

# Aggregate metadata
metadata <-
  metadata0 |> 
  left_join(hasmortality, by = "nsalid1") |> 
  left_join(inAICcalculations, by = "nsalid1") |> 
  left_join(climate_labels, by = "KP_climate_zone") |> 
  left_join(clusters, by = "nsalid1") |> 
  #left_join(citycodes, by = join_by(nsalid1 == SALID1)) |> 
  filter(hasMortality == 1,
         inAIC == 1) |> 
  mutate(across(contains("ECDF"), as_factor),
         across(contains("KP_"), as_factor),
         across(c(nsalid1, country_num, country, hasMortality), as_factor)) |> 
  select(-c(SALID1, `_TYPE_`, `_FREQ_`)) |> 
  select(nsalid1, everything())

# Save aggregated metadata
saveRDS(metadata,   here("data", "metadata.rds"))

# Read city specific data -------------------------------------------------

# Raw file locations of data
list_datafiles <- 
  list(all_cause = path_all_cause,     # All cause mortality
       respiratory = path_respiratory, # Respiratory deaths
       cardio = path_cardio) |>        # Cardiovascular deaths
  map(\(p) paste0(p, "/c", metadata$nsalid1, ".sas7bdat")) # Files for each city

# Function to read in sas7bdat file, convert to dataframe, and format
process_data <- function(path_to_data) {
  path_to_data |>
    haven::read_sas() |>
    mutate(country = toupper(country), 
           nsalid1 = as.factor(nsalid1),
           date = ymd(allDate),
           sex  = factor(male, levels = c(0, 1), labels = c("Female", "Male")),
           age = factor(thisage, levels = c("Infant", "1-4", "5-19", "20-34", "35-49",  "50-64", "65+")),
           deaths = deaths,
           tmean = ADtemp_pw,
           pop = pop_count, .keep = "none") |>
    arrange(across(c(country, nsalid1, date, sex, age))) |> # Order rows
    filter(!is.na(tmean)) # Remove days with missing temperature
}

# Create dataframe containing each city/mortality type and order columns
df <- 
  pmap(list_datafiles,
       \(all_cause, respiratory, cardio) {
         # Variables that uniquely specify a row
         each_row <- c("nsalid1", "date", "sex", "age")
         
         # Process data from a single city for all cause, respiratory, 
         # and cardiovascular deaths
         df1 <- map(all_cause, process_data)   |> list_rbind()
         df2 <- map(respiratory, process_data) |> list_rbind() |> select(all_of(each_row), respiratory = deaths)
         df3 <- map(cardio, process_data)      |> list_rbind() |> select(all_of(each_row), cardio = deaths)
         
         # Return a single dataframe with 3 death columns:
         # all cause, respiratory, and cardio
         df1 |> 
           left_join(df2, by = each_row, keep = FALSE) |> 
           left_join(df3, by = each_row, keep = FALSE)
       }) |> 
  list_rbind() |> # Combine dataframes for each city into one dataframe
  select(c(country, nsalid1, date, sex, age, # Unique determinant of each row
           deaths, respiratory, cardio,      # Deaths by type
           tmean, pop)) |> # Pop weighted mean daily temperature and log pop 
  arrange(across(c(nsalid1, date))) |> # Arrange by city and date
  mutate(country = factor(ifelse(country %in% c("PA", "GT", "SV", "CR"), "CA", country))) # Group Central American countries

# Write multi-city object to disk
write_csv(df, here("data", "mort_temp.csv"))
saveRDS(df,   here("data", "mort_temp.rds"))