library(tidyverse)
library(glue)
library(here)
library(SALURhelper)
library(dlnm)
library(gnm)
library(splines)

# Helper functions --------------------------------------------------------

# Analyze a given city with specified parameters
analyze_city <- function(df_city, death_var = "deaths") {
  # Specify knots for given city
  n_lag <- 21
  lag_knots <- logknots(n_lag, nk = 3)

  # Temperature percentiles are calculate manually
  pred_knots <- quantile(df_city$tmean, c(.1, .75, .9), na.rm = TRUE)

  # Define cross-basis
  cbt <- crossbasis(
    df_city$tmean,
    lag = n_lag,
    argvar = list(fun = "ns", knots = pred_knots),
    arglag = list(fun = "ns", knots = lag_knots),
    group = df_city$group
  )

  # Define reduced cross-basis for future meta regression
  cbt_red <- onebasis(df_city$tmean, fun = "ns", knots = pred_knots)

  # Fit the model
  model <- gnm(
    "{death_var} ~ cbt" |> glue() |> as.formula(), # Select death variable
    family = quasipoisson(),
    eliminate = strata,
    data = df_city
  )

  # Object this function returns, note that model objects can be extremely large.
  # This list contains only the relevant information calculated from the model fitting.
  list(
    cb = cbt,
    cb_red = cbt_red,
    pred_seq = seq(min(df_city$tmean), max(df_city$tmean), length.out = 100), # Number sequence for crosspred `at` parameter
    coef = coef(model),
    vcov = vcov(model),
    link = model$family$link,
    death_var = death_var
  )
}

# Get reduced prediction for a given city
city_pred <- function(city_model) {
  # Make the initial prediction to get estimate the MMT
  pred <- crossreduce(
    basis = city_model$cb,
    coef = city_model$coef,
    vcov = city_model$vcov,
    model.link = city_model$link,
    at = city_model$pred_seq
  )

  # Make the final prediction centered at the MMT
  crossreduce(
    basis = city_model$cb,
    coef = city_model$coef,
    vcov = city_model$vcov,
    model.link = city_model$link,
    at = city_model$pred_seq,
    cen = get_MMT(pred)
  )
}


# Prepare data ------------------------------------------------------------

# Read in data
df <- readRDS(here("data", "mort_temp.rds"))

# Deaths from all ages/sexes, split by city
list_city_df <-
  split(df, df$nsalid1) |> # Split data by city
  map(\(df_city) {
    df_city |>
      # Add the strata variable for conditional model fitting
      mutate(
        strata = factor(glue(
          "{year}:{month}:{dow}:{sex}:{age}",
          year = year(date),
          month = month(date),
          dow = wday(date, label = TRUE)
        )),
        group = factor(glue("{sex}:{age}"))
      )
  })

# Deaths from ages 65+/all sexes, split by city
list_city_df_65 <-
  split(df, df$nsalid1) |> # Split data by city
  map(\(df_city) {
    df_city |>
      # Restrict data to only ages 65 or older
      filter(age == "65+") |>
      mutate(
        strata = factor(glue(
          "{year}:{month}:{dow}:{sex}",
          year = year(date),
          month = month(date),
          dow = wday(date, label = TRUE)
        )),
        group = sex
      )
  })

# Get vector of city identifiers
city_names <- names(list_city_df)

# Analyze all cities ------------------------------------------------------

# Our collection of data sets and response variables
data_set <- list(all_ages = list_city_df, over_65 = list_city_df_65)
death_var <- c("deaths") # noninjury, cardio, respinf, respdis

# A data frame organizing all analyses combinations
df_analyses <- expand_grid(data_set, death_var)

# A character vector of each analyses' name
analyses_names <-
  expand_grid(age_cat = names(data_set), death_var = death_var) |>
  pmap_chr(\(age_cat, death_var) glue("{age_cat}_{death_var}"))

# Analyze all causes of death and age categories for each city
list_city_model <-
  # Map over age categories and death variable
  pmap(df_analyses, \(data_set, death_var) {
    # Map over chosen list of city data frames to model chosen death var
    data_set |>
      map(\(df_city) {
        analyze_city(df_city, death_var)
      })
  }) |>
  set_names(analyses_names)

# Predict temperature-mortality curve for each city -----------------------

# Predict the reduced ERF for each city in our age/death categories
list_city_pred <- list_city_model |>
  map(\(age_death_cat) {
    age_death_cat |>
      map(\(city_model) {
        city_pred(city_model)
      })
  })
