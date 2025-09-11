# General Libraries
library(tidyverse)
library(here)
library(patchwork)
library(lubridate)
library(glue)
library(sf)

# ggplot2 theming
theme_set(theme_bw())
theme_update(plot.title = element_text(hjust = 0.5),
             plot.subtitle = element_text(hjust = 0.5),
             strip.background = element_blank(),
             strip.text = element_text(size = 12))


# Prep for spatial plots --------------------------------------------------
# Read shape files
df_l1_sf <- read_sf("//files.drexel.edu/colleges/SOPH/Shared/UHC/Projects/Wellcome_Trust/Data Methods Core/Geodatabases/SALURBAL/Shapefiles/SALURBAL_L1.shp") |> mutate(SALID1 = factor(SALID1))
df_la_sf <- read_sf(here("data", "LA_outline", "LA_outline.shp"))

# Default ranges of longitude and latitude
range_long <- c(-120, -38)
range_lat  <- c(-55, 30)

# Default fig width/height
fig_width  <- 4
fig_height <- 9

# Figure 1 ----------------------------------------------------------------

# Read data for temperature
df <- readRDS(here("data", "mort_temp.rds")) 
# df_metadata <- readRDS(here("data", "metadata.rds")) 
# 
# df_temp <- df_metadata |> 
#   mutate(nsalid1 = factor(SALID1),
#          tmean = mean, .keep = "none")

# Calculate annual mean temperatures
df_temp <- df |> 
  # Add year variable for grouping
  mutate(year = year(date)) |> 
  # Extract each all temperature without duplicates
  summarize(tmean = first(tmean), .by = c(country, nsalid1, year, date)) |> 
  # Get yearly mean temperature for each L1
  summarize(tmean = mean(tmean), .by = c(country, nsalid1, year)) |> 
  # Get average yearly mean temperature for each L1
  summarize(tmean = mean(tmean), .by = c(country, nsalid1))

# Combine shape and temperature information
df_l1_temp <- right_join(df_l1_sf, df_temp, by = join_by(SALID1 == nsalid1))

# Create ggplot object
ggplot() +
  # Plot country outlines
  geom_sf(data = df_la_sf) +
  # Plot city temperatures
  geom_sf(aes(fill=tmean), shape = 21, data = df_l1_temp) +
  # Specify the x & y coordinates of the plot
  coord_sf(xlim = range_long, ylim = range_lat) +
  # Select the color scheme
  scale_fill_viridis_c(option = "plasma") +
  # Add labels
  labs(fill = "Temp. (°C)", x = "Longitude", y = "Latitude") +
  # Minor adjustments to the legend
  theme(legend.position = c(.1, .16),
        legend.background = element_rect(color = "grey"),
        legend.title = element_text(size = 6),
        legend.key.height = unit(.15, 'in'),
        legend.key.width = unit(.15, 'in'))

# Save the figure 1 as a png file
ggsave("results/figure_1.png", width = fig_width, height = fig_height, units = 'in')

# Figure 2 ----------------------------------------------------------------
# source("02_meta_analysis.R")

fig_2_cities <- c("8101112", "102190", "204191", "105113", "204141", "103116")
fig_2_titles <- c("Buenos Aires, Argentina",
                  "Rio de Janeiro, Brazil",
                  "Mérida, Mexico",
                  "Lima, Peru",
                  "Ciudad de México, Mexico",
                  "Los Ángeles, Chile")

fig_2_plots <- map2(fig_2_cities, fig_2_titles, \(ID, title) {
  # Extract P01, P05, MMT, P95, and P99
  p01_99 <- quantile(list_city_df[[ID]]$tmean, c(0.01, 0.99))
  p05_95 <- quantile(list_city_df[[ID]]$tmean, c(0.05, 0.95))
  MMT    <- get_MMT(list_blup_pred$all_ages_deaths[[ID]])
  x_lims <- c(min(list_city_df[[ID]]$tmean), max(list_city_df[[ID]]$tmean))
  mean_temp <- round(mean(list_city_df[[ID]]$tmean), 1)
  
  # Plot RRs
  plot_RR <- list_blup_pred$all_ages_deaths[[ID]] |> 
    get_df_RRs() |> 
    mutate(Temperature = exposure,
           type = ifelse(Temperature >= MMT, "Hot temperatures", "Cold temperatures")) |> 
    ggplot(aes(x = Temperature)) +
    geom_hline(yintercept = 1) + 
    geom_line(aes(y = RR, color = type)) +
    geom_ribbon(aes(ymin = RR_low, ymax = RR_high), alpha = .2) + 
    scale_y_continuous(limits = c(.75, 3.5),
                       breaks = c(0.8, 1, 1.25, 1.5, 2)) +
    scale_x_continuous(limits = x_lims) +
    coord_trans(y = "log") +
    scale_color_manual(values = c("Hot temperatures" = "red", "Cold temperatures" = "blue")) +
    labs(title = title, subtitle = glue("Mean temperature {mean_temp}°C"), color = "") +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank())
  
  # Plot temperature histogram
  plot_hist <- list_city_df[[ID]] |> 
    ggplot(aes(tmean)) +
    geom_histogram(bins = 30, color = "black", fill = "grey") +
    scale_y_continuous(position = "right") +
    scale_x_continuous(limits = x_lims) +
    labs(x = "Temperature (°C)", y = "")
  
  # Combine plots and add vertical lines at temp percentiles
  plot_RR / 
    plot_hist &
    geom_vline(xintercept = p01_99, linetype = "dotdash") &
    geom_vline(xintercept = p05_95, linetype = "dashed") &
    geom_vline(xintercept = MMT,    linetype = "dotted")
  })

# Create figure 2
wrap_plots(fig_2_plots, guides = "collect") & 
  theme(legend.position = "bottom", legend.box.background = element_rect(color = "black"))

# Save figure 2
ggsave("results/figure_2.png", width = 10, height = 8, units = 'in')

# Figure 3 ----------------------------------------------------------------
# source("03_RRs.R")

# Create data frame with RR slopes for all cause all ages
df_RR_slope <- list_red_95_99$all_ages_deaths |> 
  mutate(RR = exp(log_RR),
         RR_factor = case_when(
           RR <= 1                ~ "≤ 1.000",
           RR > 1 & RR < 1.05     ~ "1.001 - 1.049",
           RR >= 1.05 & RR < 1.1  ~ "1.050 - 1.099",
           RR >= 1.10 & RR < 1.15 ~ "1.100 - 1.149",
           RR >= 1.15 ~ "1.150 +",
           .default = "ISSUE"
         ) |> factor())

# Combine shape and temperature information
df_l1_slope <- right_join(df_l1_sf, df_RR_slope, by = join_by(SALID1 == nsalid1))

# Create plot
ggplot() +
  geom_sf(data = df_la_sf) +
  geom_sf(aes(fill=RR_factor), shape = 21, data = filter(df_l1_slope, exposure == "heat")) +
  coord_sf(xlim = range_long, ylim = range_lat) +
  scale_fill_viridis_d(option = "plasma", direction = -1) +
  labs(fill = "RR per 1° C above P95", x = "Longitude", y = "Latitude") +
  theme(legend.position = c(.18, .16),
        legend.background = element_rect(color = "grey"),
        legend.title = element_text(size = 6),
        legend.key.height = unit(.15, 'in'),
        legend.key.width = unit(.15, 'in'))

# Save the figure 1 as a png file
ggsave("results/figure_3.png", width = fig_width, height = fig_height, units = 'in')

# Extended Data: Figure 1 -------------------------------------------------


# Extended Data: Figure 2 -------------------------------------------------


