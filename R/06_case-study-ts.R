# Plot the VPD/AF time series for the case study fire
library(sf)
library(tidyverse)
library(vroom)
library(patchwork)

fired <- st_read("data/fired/events_w_attributes.gpkg") %>%
  st_transform(4326)

case_id <- 74238

fire <- fired %>%
  filter(id == case_id)

events <- vroom("out/goes-af-era5.csv")

d <- events %>%
  filter(id == case_id) %>%
  filter(datetime_utc < as.Date("2017-10-20"), 
         datetime_utc > as.Date("2017-10-07"))



case_study_df <- d %>%
  select(datetime_utc, n, vpd_kPa) %>%
  rename(`GOES Active Fire detections` = n, 
         `VPD (kPa)` = vpd_kPa) %>%
  pivot_longer(-datetime_utc)

case_study_df %>%
  mutate(fire_name = "Tubbs", id = case_id) %>%
  write_csv("tubbs-hourly.csv")



# Time series plot --------------------------------------------------------

thresholds <- read_csv("out/zero-goes-af-vpd-thresholds.csv")

event_names <- tibble(
  id = c(64121, 73403, 74238), 
  Event = c("Rough Fire, California", 
           "Rice Ridge Fire, Montana", 
           "Tubbs Fire, California"
           )
)

long_events <- events %>%
  group_by(id) %>%
  filter(n > 0) %>%
  summarize(n = n()) %>%
  top_n(15, n) %>%
  arrange(n)

long_pts <- events %>%
  filter(id %in% c(64121, 73403, 74238)) %>%
  left_join(select(as_tibble(fired), id, lc_name)) %>%
  left_join(thresholds) %>%
  group_by(id) %>%
  mutate(naf = sum(n)) %>%
  ungroup %>%
  left_join(event_names) %>%
  mutate(id = fct_reorder(as.factor(id), -naf))

last_date_df <- long_pts %>%
  group_by(id) %>%
  filter(n > 0) %>%
  summarize(last_date = max(datetime_utc) + 60 * 60 * 24 * 1, 
            first_date = min(datetime_utc) - 60 * 60 * 24 * 1)

split_events <- long_pts %>%
  select(datetime_utc, n, id, vpd_kPa, Event, vpd_thresh_kpa, 
         solar_radiation_w_per_m2) %>%
  left_join(last_date_df) %>%
  filter(datetime_utc <= last_date, 
         datetime_utc >= first_date) %>%
  pivot_longer(cols = c("n", "vpd_kPa")) %>%
  mutate(value = ifelse(!(name == "n" & value == 0), value, NA)) %>% 
  mutate(name = ifelse(name == "n", "AF counts", "VPD (kPa)")) %>%
  split(.$id)


plot_list <- split_events %>%
  lapply(function(x) {
    event_name <- unique(x$Event)
    
    pct_night <- x %>%
      filter(name == "AF counts", !is.na(value)) %>%
      mutate(day = solar_radiation_w_per_m2 > .001) %>%
      group_by(day) %>%
      summarize(n_af = sum(value)) %>%
      ungroup()
    
    pct_night <- pct_night$n_af[!pct_night$day] / sum(pct_night$n_af)
    title <- paste0(event_name, ": ", 100 * round(pct_night, 2), 
                    "% nighttime fire detections")
    
    ggplot(x, aes(datetime_utc, value, color = solar_radiation_w_per_m2 > .001)) + 
      geom_point(size = .6, alpha= .8) + 
      geom_line(alpha = 0.3, color = "black") + 
      theme_minimal() + 
      facet_wrap(~name, scales = "free_y", ncol = 1, strip.position = "left") + 
      xlab("") + 
      ylab("") + 
      theme(panel.grid.minor = element_blank(), 
            legend.position = "none") +
      scale_color_manual(values = c("#2166ac", "red")) + 
      ggtitle(title) + 
      geom_hline(aes(yintercept = vpd_thresh_kpa),
                 data = distinct(x, vpd_thresh_kpa) %>%
                   mutate(name = "VPD (kPa)"),
                 linetype = "dashed") + 
      theme(plot.title = element_text(size=12))
  })

p <- wrap_plots(plot_list, ncol = 1)
p
ggsave("fig/case-study-ts.pdf", plot = p, width = 8, height = 6.5)
ggsave("fig/case-study-ts.png", plot = p, width = 8, height = 6.5)

