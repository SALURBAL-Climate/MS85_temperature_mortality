library(tidyverse)
library(here)
library(glue)

# Helper functions --------------------------------------------------------

# Function to read in sas7bdat file, convert to dataframe, and format
process_data <- function(path_to_data) {
  path_to_data |>
    haven::read_sas() |>
    mutate(country = toupper(country), 
           nsalid1 = nsalid1,
           date = ymd(allDate),
           sex  = factor(male, levels = c(0, 1), labels = c("Female", "Male")),
           age = factor(thisage, levels = c("Infant", "1-4", "5-19", "20-34", "35-49",  "50-64", "65+")),
           deaths = deaths,
           # noninjury = noninjury,
           # cardio = cardio,
           # respinf = respinf,
           # respdis = respdis,
           tmean = L1ADtemp_pw, .keep = "none") |>
    filter(!is.na(tmean)) # Remove days with missing temperature
}

# Read city specific data -------------------------------------------------

# Path to updated data mortality/temperature data
path_all_cause <- "//files.drexel.edu/colleges/SOPH/Shared/UHC/Projects/Wellcome_Trust/Manuscripts/MS85_Update/2025_08_20/data/clean/"

# Files for each city
list_datafiles <- list.files(path_all_cause,         # path to data
                             pattern = "*.sas7bdat", # filename pattern
                             full.names = TRUE)      # include path to file

# Read files into a single data frame
df <- list_datafiles |> 
  map(process_data) |> 
  list_rbind() |> 
  mutate(across(c(country, nsalid1), factor),
         # Combine Central American countries
         country = fct_collapse(country, CA = c("PA", "GT", "SV", "CR"))) |> 
  arrange(across(c(country, nsalid1, date, sex, age))) # Order rows

# Write cleaned data frame to disk
saveRDS(df,   here("data", "mort_temp.rds"))
