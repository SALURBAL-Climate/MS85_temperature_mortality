# source("02_meta_analysis.R")

# Total number of deaths
list_total_deaths <-
  pmap(df_analyses,
       \(data_set, death_var) {
         data_set |> 
           map_dbl(\(df_city) {
             # Sum within city
             sum(df_city[[death_var]])
           }) |> sum() # Sum across cities
    }) |> set_names(analyses_names)

# Get attributable number for each city, temp range, and analysis
list_blup_ANs <-
  # Cycle through analyses for each cause of death/age-cat combination
  pmap(list(list_model = list_city_model,
            list_BLUP  = list_blup_meta,
            list_pred  = list_blup_pred),
       \(list_model, list_BLUP, list_pred) {
         pmap(
           # Information to get AN for each city
           list(df_city = list_city_df,
                model   = list_model,
                BLUP    = list_BLUP,
                pred    = list_pred),
           # Calculation of said AN
           \(df_city, model, BLUP, pred) {
             # Specify temperature ranges of interest
             temp_ranges <- list(
               range_total        = range(df_city$tmean),
               range_all_heat     = c(get_MMT(pred), max(df_city$tmean)),
               range_extreme_heat = quantile(df_city$tmean, c(.95, 1)),
               range_all_cold     = c(min(df_city$tmean), get_MMT(pred)),
               range_extreme_cold = quantile(df_city$tmean, c(0, .05))
             )
             
             temp_ranges |> 
               map(\(range) {
                 # Get AN (for EDF estimate)
                 AN <- attrdl(df_city$tmean,
                              basis = model$cb,
                              cases = df_city[[model$death_var]],
                              coef = BLUP$blup,
                              vcov = BLUP$vcov,
                              model.link = "log",
                              type = "an",
                              dir = "forw",
                              cen = get_MMT(pred),
                              range = range)
                 
                 # Get 1000 simulated ANs (for lower/upper bounds of EDF estimate)
                 sim_ANs <- attrdl(df_city$tmean,
                                   basis = model$cb,
                                   cases = df_city[[model$death_var]],
                                   coef = BLUP$blup,
                                   vcov = BLUP$vcov,
                                   model.link = "log",
                                   type = "an",
                                   dir = "forw",
                                   cen = get_MMT(pred),
                                   range = range,
                                   sim = TRUE,
                                   nsim = 1000)
                 
                 list(AN = AN,
                      sim_ANs = sim_ANs)
               })
           })
       }, .progress = TRUE)

# Save list of ANs estimated from BLUP
saveRDS(list_blup_ANs, here("results", "list_blup_ANs.rds"))
readRDS(here("results", "list_blup_ANs.rds"))

# Estimate overall ANs for each analysis/temp range
df_blup_AN <- list_blup_ANs |> 
  map(\(age_death_cat) {
    age_death_cat |> 
      map(\(city) {
        city |> 
          map_dbl(\(temp_range) {
            temp_range$AN
          }) |> as_tibble_row()
      }) |> list_rbind(names_to = "nsalid1")
  })

# Estimate EDF for each analysis/temp range
df_AF <- 
  map2(df_blup_AN, list_total_deaths,
       \(age_death_cat, death_total) {
         age_death_cat |> 
           summarize(across(contains("range"), \(x) 100*sum(x)/death_total)) |> 
           mutate(conf = "center", .before = range_total)
  }) |> list_rbind(names_to = "analysis")

# Estimate overall simulated ANs for each analysis/temp range
df_blup_sim_ANs <- list_blup_ANs |> 
  map(\(age_death_cat) {
    age_death_cat |> 
      map(\(city) {
        city |> 
          map(\(temp_range) {
            temp_range$sim_ANs
          }) |> as_tibble() |> mutate(iter = 1:1000, .before = range_total)
      }) |> list_rbind(names_to = "nsalid1")
  })


# Estimate EDF confidence intervals for each analysis/temp range
df_AF_conf <- 
  map2(df_blup_sim_ANs, list_total_deaths,
       \(age_death_cat, total_deaths) {
         age_death_cat |> 
           summarize(across(contains("range"), \(x) 100*sum(x)/total_deaths), .by = iter) |> 
           reframe(conf = c("lower", "upper"),
                   across(contains("range"), \(x) quantile(x, c(.025, .975))))
       }) |> list_rbind(names_to = "analysis")

# Final estimate
df_EDF <- 
  bind_rows(df_AF, df_AF_conf) |> 
  arrange(analysis)

saveRDS(df_EDF, here("results", "EDFs.rds"))
