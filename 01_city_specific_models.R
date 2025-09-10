library(tidyverse)
library(glue)
library(here)
library(SALURhelper)
library(dlnm); library(gnm); library(splines)

# Helper functions --------------------------------------------------------

# Analyze a given city with specified parameters
analyze_city <- function(df_city, death_var = "deaths", df_meta = NULL) {
  
  # Specify knots for given city
  n_lag <- 21
  lag_knots <- logknots(n_lag, nk = 3)
  
  # If metadata temperature percentiles are provided, use those, otherwise, 
  # temperature percentiles are calculate manually
  pred_knots <- quantile(df_city$tmean, c(.1, .75, .9), na.rm = TRUE) # manual
  if(!is.null(df_meta)) pred_knots <- c(df_meta$p10, df_meta$p75, df_meta$p90) # pre-calcualted
  
  # Define cross-basis
  cbt <- crossbasis(df_city$tmean,
                    lag = n_lag,
                    argvar = list(fun = "ns", knots = pred_knots),
                    arglag = list(fun = "ns", knots = lag_knots), 
                    group = df_city$group)
  
  # Define reduced cross-basis for future meta regression
  cbt_red <- onebasis(df_city$tmean, fun = "ns", knots = pred_knots)
  
  # Fit the model
  model <- gnm("{death_var} ~ cbt" |> glue() |> as.formula(), # Select death variable
               family = quasipoisson(),
               eliminate = strata,
               data = df_city)
  
  # Object this function returns, note that model objects can be extremely large.
  # This list contains only the relevant information calculated from the model fitting.
  list(cb = cbt,
       cb_red = cbt_red,
       pred_seq = seq(min(df_city$tmean),
                      max(df_city$tmean), 
                      length.out = 100), # Number sequence for crosspred `at` parameter
       coef = coef(model), 
       vcov = vcov(model),
       link = model$family$link,
       death_var = death_var)
}

# Get reduced prediction for a given city
city_pred <- function(city_model) {
  # Make the initial prediction to get estimate the MMT
  pred <- crossreduce(basis = city_model$cb, 
                      coef  = city_model$coef,
                      vcov  = city_model$vcov,
                      model.link = city_model$link,
                      at = city_model$pred_seq)
  
  # Make the final prediction centered at the MMT
  crossreduce(basis = city_model$cb, 
              coef  = city_model$coef,
              vcov  = city_model$vcov,
              model.link = city_model$link, 
              at = city_model$pred_seq, cen = get_MMT(pred))
}


# Prepare data ------------------------------------------------------------

# Read in data
df <- readRDS(here("data", "mort_temp.rds"))
df_metadata <- readRDS(here("data", "metadata.rds"))

# Extract city names  
city_names <- unique(df$nsalid1)

# Deaths from all ages/sexes, split by city
list_city_df <- df |> 
  group_by(nsalid1) |> 
  group_split() |> # Split data by city
  map(\(df_city) {
    df_city |> 
      # Get totals for all cause, respiratory, and cardio deaths by city-day
      # summarize(across(c(deaths, respiratory, cardio), sum), 
      #           across(c(country, tmean, pop), first), .by = c(nsalid1, date)) |> 
      # Add the strata variable for conditional model fitting (year:month:dow)
      mutate(strata = factor(paste(year(date), 
                                   month(date), 
                                   wday(date, label = TRUE),
                                   sex, 
                                   age, sep = ":")),
             group = factor(paste(sex, age, sep = ":")))
  }) |> set_names(city_names) # Label this list of dataframes with the city names

# Deaths from ages 65+/all sexes, split by city
list_city_df_65 <- df |> 
  group_by(nsalid1) |> 
  group_split() |> # Split data by city
  map(\(df_city) {
    df_city |> 
      # Restrict data to only ages 65 or older
      filter(age == "65+") |> 
      # Get totals for all cause, respiratory, and cardio deaths by city-day
      # summarize(across(c(deaths, respiratory, cardio), sum), 
      #           across(c(country, tmean, pop), first), .by = c(nsalid1, date)) |> 
      # Add the strata variable for conditional model fitting (year:month:dow)
      mutate(strata = factor(paste(year(date), 
                                   month(date), 
                                   wday(date, label = TRUE),
                                   sex, sep = ":")),
             group = sex)
  }) |> set_names(city_names) # Label this list of dataframes with the city names

# Meta data, split by city
list_metadata <- df_metadata |> 
  group_by(nsalid1) |> 
  group_split() |> 
  set_names(city_names)

# Analyze all cities ------------------------------------------------------

# Our collection of data sets and response variables
data_set  <- list(list_city_df, list_city_df_65)
death_var <- c("deaths", "respiratory", "cardio")

# All combinations that will be used in our analysis
df_analyses <- expand_grid(data_set, death_var)

name_analyses <- 
  expand_grid(x = c("all_ages", "over_65"),
              y = death_var) |> pmap_chr(\(x, y) paste(x, y, sep = "_"))

# Analyze each city for all causes of death and age categories
list_city_model <-
  pmap(df_analyses, 
       \(data_set, death_var) {
         # Map over each city in chosen data set
         map2(data_set, list_metadata,
              \(df_city, df_meta) {
                analyze_city(df_city, death_var, df_meta = NULL)
           })
        }) |> set_names(name_analyses)

# Predict temperature-mortality curve for each city -----------------------

# Predict the reduced ERF for each city in our age/death categories
list_city_pred <- list_city_model |> 
  map(\(age_death_cat) {
    age_death_cat |> 
      map(\(city_model) {
        city_pred(city_model)
      })
  })
