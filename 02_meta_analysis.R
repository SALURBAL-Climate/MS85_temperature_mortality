#source("01_city_specific_models.R")
library(mixmeta)

# Set up meta-regression --------------------------------------------------

# Predictors by city with pre-calculated metadata
# df_meta <- df_metadata |>
#   mutate(nsalid1 = nsalid1,
#          p50 = p50,
#          temp_range = max - min,
#          range_high = max - p50,
#          range_low  = p50 - min,
#          country = country_num,
#          KP_level1_Interpretation = KP_level1_Interpretation, 
#          .keep = "none")

# predictors by city with manually calculated metadata
df_meta <- df |>
  summarize(p50 = median(tmean),
            temp_range = max(tmean) - min(tmean),
            range_high = max(tmean) - p50,
            range_low  = min(tmean) - p50,
            country = first(country), .by = nsalid1) |>
  left_join(df_metadata, by = join_by(nsalid1 == nsalid1), suffix = c("", "_meta"))

# List of matrices of coefficients (response of meta-regression)
list_mat_coef <- list_city_pred |> 
  map(\(age_death_cat) {
    age_death_cat |> 
      map_dfr(\(x) x$coefficients) |> 
      as.matrix()
  })

# List of list of list of variance-covariance matrices
list_list_vcov <- list_city_pred |> 
  map(\(age_death_cat) {
    age_death_cat |> 
      map(\(x) x$vcov)
  })

# Compute meta-regression for each death/age category
list_model_meta <- 
  pmap(list(mat_coef  = list_mat_coef,
            list_vcov = list_list_vcov),
       \(mat_coef, list_vcov) {
         mixmeta(mat_coef ~
                   ns(p50, df = 2) +
                   #ns(temp_range, df = 2) +
                   ns(range_high, df = 2) +
                   ns(range_low, df = 2) +
                   country,
                   #KP_level1_Interpretation, 
                 S = list_vcov,
                 data = df_meta)
       })

# Best linear unbiased predictors for each city
list_blup_meta <- map(list_model_meta, \(x) blup(x, vcov = TRUE) |> set_names(city_names))

# Predicted ERF for each city
list_blup_pred <-
  pmap(list(city_models = list_city_model,
            blup_metas  = list_blup_meta),
       \(city_models, blup_metas) {
         map2(city_models, blup_metas, 
              \(model, BLUP) {
                pred <- crosspred(basis = model$cb_red,
                                  coef  = BLUP$blup,
                                  vcov  = BLUP$vcov,
                                  model.link = model$link,
                                  at = model$pred_seq)
                
                crosspred(basis = model$cb_red, 
                          coef  = BLUP$blup,
                          vcov  = BLUP$vcov,
                          model.link = model$link, 
                          at = model$pred_seq, cen = get_MMT(pred))
              })
       })
