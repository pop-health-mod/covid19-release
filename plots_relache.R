# Plot conversion for the Relâche paper
# July 2020
setwd("~/Desktop/covid19")
devtools::load_all("~/Desktop/covid19")

# Loading data from the last calibration
library(dplyr)
library(ggplot2)
load("output/output_relache_rw_2020-07-01")

# ---- Plotting the fit ----

# We must first create an appropriate dataframe for plotting
nh <- out[["nh"]] %>%
  mutate(date = start_date + (SEIR_Const$iter_time * SEIR_Const$dt) - SEIR_Const$dt) %>%
  mutate(outcome = "Hospitalisations")

# We then create a data
r0 <- tibble(median = out[["r0_med"]], lci = out[["r0_lci"]],
             uci = out[["r0_uci"]], l50 = out[["r0_l50"]],
             u50 = out[["r0_u50"]],
             date = start_date + (SEIR_Const$iter_time * SEIR_Const$dt) - SEIR_Const$dt,
             outcome = "Rt")

# We then bind the data frame into a single data instance
df <- rbind(nh, r0)

# Empirical data
df_empi <- tibble(date = start_date + SEIR_Const$y_time * SEIR_Const$dt,
                  hosp = SEIR_Data$nh) %>%
  filter(!is.na(hosp)) %>%
  mutate(outcome = "Hospitalisations")

# Data for the ylines
data_y <- tibble(xint = 1,
                 outcome = "Rt")

# Data for segments and text
int_1 <- start_date + SEIR_Const$intervention1_ind * SEIR_Const$dt - SEIR_Const$dt
data_seg1 <- tibble(x = int_1,
                   ymin = df$median[which(df$date == int_1 & df$outcome == "Rt")],
                   ymax = df$median[which(df$date == int_1 & df$outcome == "Rt")] + 1,
                   outcome = "Rt")

int_2 <- start_date + SEIR_Const$intervention2_ind * SEIR_Const$dt - SEIR_Const$dt
data_seg2 <- tibble(x = int_2,
                    ymin = df$median[which(df$date == int_2 & df$outcome == "Rt")],
                    ymax = df$median[which(df$date == int_2 & df$outcome == "Rt")] + 1,
                    outcome = "Rt")

# Data for panel title
data_panel <- tibble(x = c(min(df$date), min(df$date)),
                     y = c(170, 3.25),
                     outcome = c("Hospitalisations", "Rt"),
                     label = c("A) New hospitalisations (Daily)",
                               "B) Reproductive number - R(t)"))

# Data for custom legend
data_legend <- tibble(x = rep(min(df$date), 4),
                      y = c(0.5, 0.75, 35, 45),
                      outcome = c(rep("Rt", 2), rep("Hospitalisations", 2)),
                      labels = c("50% CrI", "95% CrI", "50% CrI", "95% CrI"),
                      coloring = c(rgb(151, 148, 97, 153, maxColorValue = 255),
                                   rgb(151, 148, 97, 102, maxColorValue = 255),
                                   rgb(205, 87, 51, 153, maxColorValue = 255),
                                   rgb(205, 87, 51, 102, maxColorValue = 255)))


# See if we can get a quick two-panel plot with the appropriate structure
plot1 <- ggplot(df %>% filter(date <= max(df_empi$date)), aes(x = date)) +
  theme_minimal() +
  geom_line(aes(y = median, col = outcome), show.legend = F, size = 1) +
  geom_ribbon(aes(ymin = lci, ymax = uci, fill = outcome),
              alpha = 0.4, show.legend = F) +
  geom_ribbon(aes(ymin = l50, ymax = u50, fill = outcome),
              alpha = 0.6, show.legend = F) +
  geom_point(data = df_empi, aes(x = date, y = hosp),
             col = "grey20", fill = "grey20", alpha = 0.5) +
  geom_hline(data = data_y, aes(yintercept = xint), color = "grey40",
             lty = 2) +
  geom_segment(data = data_seg1, aes(x = x, xend = x, y = ymin, yend = ymax),
               col = "grey40", alpha = 0.8, lty=3) +
  geom_segment(data = data_seg2, aes(x = x, xend = x, y = ymin, yend = ymax),
               col = "grey40", alpha = 0.8, lty=3) +
  geom_text(data = data_seg1, aes(x = x, y = ymax), label = "First measures",
            col = "grey20", alpha = 0.8, hjust = 0, nudge_x = 0.5, vjust = 0,
            size = 3,fontface = "bold") +
  geom_segment(data = data_seg2, aes(x = x, xend = x, y = ymin, yend = ymax),
               col = "grey40", alpha = 0.8, lty=3) +
  geom_text(data = data_seg2, aes(x = x, y = ymax), label = "Law on essential services",
            col = "grey20", alpha = 0.8, hjust = 0, nudge_x = 0.5, vjust = 0,
            size = 3, fontface = "bold") +
  geom_text(data = data_panel, aes(x = x, y = y, label = label), size = 4,
            col = "grey40", hjust = 0, nudge_x = 0.5, vjust = 0) +
  geom_point(data = data_legend, aes(x = x, y = y), col = data_legend$coloring,
             fill = data_legend$coloring, shape = 22, size = 2) +
  geom_text(data = data_legend, aes(x = x, y = y, label = labels), size = 3,
            col = "grey40", hjust = 0, nudge_x = 2, vjust = 0.5) +
  facet_wrap(~outcome, nrow = 2, scales = "free_y",
             strip.position = "left",
             labeller = as_labeller(c("Hospitalisations" = "Hospitalisations",
                                      Rt = "R(t)"))) +
  scale_color_manual(values = nationalparkcolors::park_palette("ArcticGates")[c(5,4)]) +
  scale_fill_manual(values = nationalparkcolors::park_palette("ArcticGates")[c(5,4)]) +
  scale_x_date(breaks = "2 weeks", date_labels = "%b %d") +
  expand_limits(y = 0) +
  labs(x = "Date") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, family = "Source Sans Pro"),
        strip.text = element_text(size = 10, family = "Source Sans Pro"),
        axis.title.x = element_text(size = 10, family = "Source Sans Pro"),
        strip.background = element_blank(),
        strip.placement = "outside",
        axis.title.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        text = element_text(family = "Source Sans Pro"))
plot1

ggsave(sprintf("quebec_ggfit_%s.png", Sys.Date()),
       path = "fig/",
       device = "png", height = 7, width = 6, units = "in",
       dpi = 320)

# ---- Plotting the output of the relâche paper ----
# Loading output from the last relâche run
today <- as.Date("2020-07-01") - 1
load( paste("./output/output_relache_rw_", today + 1, sep = ""))

start_date <- as.Date("2020-02-15")
start <- 0
dt <- 0.05
end <- as.numeric(today - start_date) +  7
time <- seq(start, end, by = dt)
niter <- (end - start) / dt + 1

import <- SEIR_Data$num_import
max_time <- as.numeric(today + 1 - start_date)

# import ontario
data(ont)
ont <- as.data.frame(ont)
ont$date <- as.Date(rownames(ont))
ont$nc_travel_minus <- NA
ont$nc_travel_minus[1:(nrow(ont) - 2)] <- ont$nc_travel[3:nrow(ont)]
ont$nc_travel_plus[3:nrow(ont)] <- ont$nc_travel[1:(nrow(ont) - 2)]
# imported cases after mandatory quarantine assume not infectious
ont$nc_travel[ont$date >= as.Date("2020-03-25")] <- 0
ont$nc_travel_minus[ont$date >= as.Date("2020-03-25")] <- 0
ont$nc_travel_plus[ont$date >= as.Date("2020-03-25")] <- 0
ont <- subset(ont, date >= start_date)
ont$time <- seq(1, length(ont$time), by = 1)

import_ont <- rep(0, niter)
import_ont[1] <- 1
import_ont[(ont$time[which(ont$nc_travel > 0)]) / dt] <- ont$nc_travel[which(ont$nc_travel > 0)]

import_ont_minus <- rep(0, niter)
import_ont_minus[1] <- 1
import_ont_minus[(ont$time[which(ont$nc_travel_minus > 0)]) / dt] <- ont$nc_travel_minus[which(ont$nc_travel_minus > 0)]

import_ont_plus <- rep(0, niter)
import_ont_plus[1] <- 1
import_ont_plus[(ont$time[which(ont$nc_travel_plus > 0)]) / dt] <- ont$nc_travel_plus[which(ont$nc_travel_plus > 0)]


# ---- we re-run to get the outcome ----
scenario <- run_scenario_relache(SEIR_Data, SEIR_Const, fit,
                                 import_scenario = SEIR_Data$num_import, rw = TRUE)
out <- get_out_scenario_relache(scenario, max_time)

# Dataframe to join
out_nh <- out$cum_nh %>%
  mutate(date = seq(start_date + 1, as.Date(start_date + max_time), by = "day")) %>%
  mutate(scenario = "Québec")

# ---- ontario - symptoms ----
scenario_ont <- run_scenario_relache(SEIR_Data, SEIR_Const, fit,
                                     import_scenario = import_ont, rw = TRUE)
out_ont <- get_out_scenario_relache(scenario_ont, max_time)

# Dataframe
out_ont_nh <- out_ont$cum_nh %>%
  mutate(date = seq(start_date + 1, as.Date(start_date + max_time), by = "day")) %>%
  mutate(scenario = "Ontario")

# ---- ontario - symptoms ----
scenario_ont_m <- run_scenario_relache(SEIR_Data, SEIR_Const, fit,
                                       import_scenario = import_ont_minus, rw = TRUE)
out_ont_m <- get_out_scenario_relache(scenario_ont_m, max_time)

out_ont_m_nh <- out_ont_m$cum_nh %>%
  mutate(date = seq(start_date + 1, as.Date(start_date + max_time), by = "day")) %>%
  mutate(scenario = "Ontario minus 1")

# ---- ontario - symptoms ----
scenario_ont_p <- run_scenario_relache(SEIR_Data, SEIR_Const, fit,
                                       import_scenario = import_ont_plus, rw = TRUE)
out_ont_p <- get_out_scenario_relache(scenario_ont_p, max_time)

out_ont_p_nh <- out_ont_p$cum_nh %>%
  mutate(date = seq(start_date + 1, as.Date(start_date + max_time), by = "day")) %>%
  mutate(scenario = "Ontario plus 1")

# ---- no more importation after March 8th ----
march02_ind <- as.numeric(as.Date("2020-03-02") - start_date) / dt + 1
import_0302 <- SEIR_Data$num_import
import_0302[march02_ind:niter] <- 0
scenario_0302 <- run_scenario_relache(SEIR_Data, SEIR_Const, fit,
                                      import_scenario = import_0302, rw = TRUE)
out_0302 <- get_out_scenario_relache(scenario_0302, max_time)

out_0302_nh <- out_0302$cum_nh %>%
  mutate(date = seq(start_date + 1, as.Date(start_date + max_time), by = "day")) %>%
  mutate(scenario = "No importation March 2")

# ---- no more importation after March 8th ----
march08_ind <- as.numeric(as.Date("2020-03-07") - start_date) / dt + 1
import_0308 <- SEIR_Data$num_import
import_0308[march08_ind:niter] <- 0
scenario_0308 <- run_scenario_relache(SEIR_Data, SEIR_Const, fit,
                                      import_scenario = import_0308, rw = TRUE)
out_0308 <- get_out_scenario_relache(scenario_0308, max_time)

out_0308_nh <- out_0308$cum_nh %>%
  mutate(date = seq(start_date + 1, as.Date(start_date + max_time), by = "day")) %>%
  mutate(scenario = "No importation March 8")

# ---- no more importation after March 14th ----
march15_ind <- as.numeric(as.Date("2020-03-14") - start_date) / dt + 1
import_0315 <- SEIR_Data$num_import
import_0315[march15_ind:niter] <- 0
scenario_0315 <- run_scenario_relache(SEIR_Data, SEIR_Const, fit,
                                      import_scenario = import_0315, rw = TRUE)
out_0315 <- get_out_scenario_relache(scenario_0315, max_time)

out_0315_nh <- out_0315$cum_nh %>%
  mutate(date = seq(start_date + 1, as.Date(start_date + max_time), by = "day")) %>%
  mutate(scenario = "No importation March 15")

# We build a data frame from which we can plot the data.
df_relache <- bind_rows(out_nh, out_0302_nh, out_0308_nh, out_0315_nh) %>%
  mutate(dashed = ifelse(scenario == "Québec", "Dashed", "Not Dashed"))

# We can then plot the data
ggplot(df_relache) +
  theme_minimal() +
  geom_line(aes(x = date, y = median, color = scenario, lty = dashed),
            show.legend = F) +
  geom_ribbon(aes(x = date, ymin = lci, ymax = uci, fill = scenario),
              show.legend = F, alpha = 0.7) +
  geom_point(data = df_relache %>% filter(date == max(date)),
             aes(x = date, y = median, color = scenario), show.legend = F) +
  scale_x_date(limits=c(as.Date("2020-03-01"), Sys.Date() + 30), expand=c(0,1),
               breaks = "2 weeks", date_labels = "%b %d") +
  scale_color_manual(values = nationalparkcolors::park_palette("ArcticGates")) +
  scale_fill_manual(values = nationalparkcolors::park_palette("ArcticGates")) +
  xlab("Date") + ylab("Cumulative number of hospitalisations") +
  theme(strip.text = element_text(face = "bold"),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        axis.text.y = element_text(size = 8),
        axis.title.x = element_text(size=9),
        axis.title.y = element_text(size=9))


