
# Visualize GAM predictions -----------------------------------------------

library(tidyverse)
library(vroom)
library(ggrepel)
library(ggthemes)
library(patchwork)
library(sf)

events <- vroom("out/events.csv")

max_vpd <- events %>%
  group_by(lc_name) %>%
  summarize(max_vpd = max(vpd_kPa))

predictions <- list.files(pattern = "-preds.csv$") %>%
  lapply(vroom) %>%
  bind_rows %>%
  left_join(max_vpd) %>%
  mutate(in_range = vpd_kPa <= max_vpd) %>%
  filter(vpd_kPa <= quantile(events$vpd_kPa, 0.99))

predictions

fire_areas <- st_read("data/fired/events_w_attributes.gpkg") %>%
  st_transform(4326) %>%
  as.data.frame %>%
  select(-geom) %>%
  distinct(id, total_area_km2, lc_name, lc) %>%
  as_tibble

partial_df <- predictions %>%
  mutate(new_log_mu = Xb + log(mean(fire_areas$total_area_km2))) %>%
  group_by(lc_name, vpd_kPa, in_range) %>%
  summarize(lo = quantile(new_log_mu, .025), 
            hi = quantile(new_log_mu, .975),
            mu = mean(new_log_mu), 
            n = n()) %>%
  ungroup

# Write a CSV file with land-cover specific colors
library(scales)
show_col(ptol_pal()(11))
tibble(lc_name = unique(partial_df$lc_name), 
       color = ptol_pal()(11)) %>%
  write_csv("landcover-colors.csv")




# Compute vpd thresholds --------------------------------------------------

zero_df <- predictions %>%
  filter(vpd_kPa <= 3) %>%
  group_by(lc_name, vpd_kPa) %>%
  summarize(median_pr_zero = median(pr_zero)) %>%
  ungroup()

thresh_df <- zero_df %>%
  mutate(abs_diff = abs(median_pr_zero - 0.95)) %>%
  group_by(lc_name) %>%
  filter(abs_diff == min(abs_diff))

zero_df %>%
  ggplot(aes(vpd_kPa, median_pr_zero)) + 
  geom_line() + 
  facet_wrap(~lc_name) + 
  geom_point(data = thresh_df)

thresholds <- thresh_df %>%
  ungroup %>%
  rename(vpd_thresh_kpa = vpd_kPa) %>%
  select(-abs_diff, -median_pr_zero) %>%
  mutate(probability_of_zero = 0.95) 

thresholds %>%
  write_csv(paste0(Sys.Date(), "_zero-goes-af-vpd-thresholds.csv"))


partial_df <- partial_df %>%
  left_join(thresholds) %>%
  mutate(lc_name = paste0(lc_name, " (", round(vpd_thresh_kpa, 1), ')'))

partial_plot <- partial_df %>%
  ggplot(aes(vpd_kPa, exp(mu), color = lc_name)) +
  geom_histogram(bins = 200, 
                 alpha = .2,
                 data = events %>%
                   filter(vpd_kPa <= quantile(events$vpd_kPa, 0.99)),
                 inherit.aes = FALSE,
                 aes(x = vpd_kPa,
                     y = exp(ifelse(..count.. <= 1, 
                                    NA, 
                                    ..count.. * .0002)))) +
  geom_path(aes(alpha = in_range)) +
  scale_y_log10(breaks = c(0.01, 1, 100), 
                labels = c(0.01, 1, 100)) +
  geom_text_repel(aes(label = lc_name, x = vpd_kPa), 
                  data = partial_df %>%
                    group_by(lc_name) %>%
                    filter(vpd_kPa == max(vpd_kPa)) %>%
                    ungroup, 
                  nudge_x      = 0.4,
                  direction    = "y",
                  hjust        = 0,
                  segment.size = 0.1, size = 3, 
                  segment.alpha = .5,
                  segment.color = "black") + 
  scale_x_continuous(breaks = c(0:4), limits = c(0, 7)) +
  theme_minimal() + 
  theme(legend.position = "none", 
        panel.grid.minor = element_blank(), 
        panel.grid.major.y = element_blank(), 
        axis.title.x = element_text(hjust = .24)) + 
  xlab("Vapor pressure deficit (kPa)") + 
  ylab("Active fire detections per hour") + 
  scale_alpha_manual(values = c(.1, 1)) + 
  scale_color_ptol()
partial_plot

dir.create("fig", showWarnings = FALSE)
ggsave(plot = partial_plot, 
       filename = "fig/vpd-partial-effects.pdf", width = 7.5, height = 3.5)
ggsave(plot = partial_plot, 
       filename = "fig/vpd-partial-effects.png", width = 7.5, height = 3.5)





# Partial effects plots, with uncertainty
bayes_df <- predictions %>%
  mutate(new_mu = exp(Xb + log(mean(fire_areas$total_area_km2))))

bayes_plot <- bayes_df %>%
  droplevels() %>%
  left_join(thresholds) %>%
  ggplot(aes(vpd_kPa, new_mu, color = lc_name)) +
  geom_path(aes(group = j), alpha = .01) +
  scale_y_log10(breaks = c(0.01, 1, 100), 
                labels = c(0.01, 1, 100)) + 
  theme_minimal() + 
  theme(legend.position = "none", 
        panel.grid.minor = element_blank()) + 
  xlab("Vapor pressure deficit (kPa)") + 
  ylab("Active fire detections per hour") + 
  scale_alpha_manual(values = c(.1, 1)) + 
  scale_color_ptol() + 
  facet_wrap(~fct_reorder(lc_name, vpd_thresh_kpa), 
              nrow = 4, , labeller = label_wrap_gen()) + 
  geom_path(data = partial_df %>%
              mutate(lc_name = trimws(gsub("\\(.*", "", lc_name))) %>%
              left_join(thresholds), 
            aes(y = exp(mu)), size = 1) + 
  coord_cartesian(ylim = c(1e-4, 1e4)) + 
  geom_rug(data = filter(events, lc_name %in% bayes_df$lc_name) %>%
             filter(vpd_kPa <= quantile(events$vpd_kPa, 0.99)) %>%
             left_join(thresholds), 
           aes(x = vpd_kPa, color = lc_name), 
           inherit.aes = FALSE, alpha = .008)
bayes_plot
ggsave("fig/bayes_plot.png", bayes_plot, width = 5, height = 6)
ggsave("fig/bayes_plot.pdf", bayes_plot, width = 5, height = 6)
