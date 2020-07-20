setwd("~/Google Drive/McGill/Research/Covid/covid19")
devtools::load_all("~/Google Drive/McGill/Research/Covid/covid19")


today <- as.Date("2020-07-01") - 1
load( paste("./output/output_relache_rw_", today + 1, sep = ""))

start_date <- as.Date("2020-02-27")
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
ont$nc_travel_minus <- 0
ont$nc_travel_minus[1:(nrow(ont) - 1)] <- ont$nc_travel[2:nrow(ont)]


data(qc)
qc <- as.data.frame(qc)
qc$date <- as.Date(rownames(qc))
qc$nc_travel_symp_minus <- 0
qc$nc_travel_symp_minus[1:(nrow(qc) - 1)] <- qc$nc_travel_symp[2:nrow(qc)]
qc <- subset(qc, date >= start_date)
qc$time <- seq(1, length(qc$time), by = 1)
sum(qc$nc_travel_symp[qc$date < as.Date("2020-03-25")])


qc$nc_travel_symp_arrival_March8_m <- 0
qc$nc_travel_symp_arrival_March8_m[1:(nrow(qc) - 1)] <- qc$nc_travel_symp_arrival_March8[2:nrow(qc)]
qc$nc_travel_symp_arrival_March15_m <- 0
qc$nc_travel_symp_arrival_March15_m[1:(nrow(qc) - 1)] <- qc$nc_travel_symp_arrival_March15[2:nrow(qc)]

# differences ON-QC
sum(qc$nc_travel_symp[qc$date < as.Date("2020-02-27")])
sum(ont$nc_travel[ont$date < as.Date("2020-02-27")])

plot(ont$nc_travel ~ ont$date, type = 'l', col = 'firebrick4', lwd = 2,
     xlim = as.Date(c("2020-02-01", "2020-04-15")), ylim = c(0, 200))
  lines(qc$nc_travel_symp_arrival ~ qc$date, col = "cadetblue3", lwd = 2, lty = 3)
  lines(qc$nc_travel_symp_arrival_March15 ~ qc$date, col = "#979461", lwd = 2, lty = 1)
  lines(qc$nc_travel_symp_arrival_March8 ~ qc$date, col = "antiquewhite3", lwd = 2, lty = 1)
  lines(qc$nc_travel_symp ~ qc$date, col = "cadetblue4", lwd = 2)

# imported cases after mandatory quarantine assume not infectious
ont$nc_travel_minus[ont$date >= as.Date("2020-03-25") + 4] <- 0
ont <- subset(ont, date >= start_date)
ont$time <- seq(1, length(ont$time), by = 1)
sum(ont$nc_travel[ont$date < as.Date("2020-03-25")])

# ---- we re-run to get the outcome ----
scenario <- run_scenario_relache(SEIR_Data, SEIR_Const, fit,
                         import_scenario = SEIR_Data$num_import, rw = TRUE)
out <- get_out_scenario_relache(scenario, max_time)

# ---- ontario - symptoms ----
import_ont_minus <- rep(0, niter)
import_ont_minus[1] <- 0
import_ont_minus[ont$time / dt] <- ont$nc_travel_minus
scenario_ont <- run_scenario_relache(SEIR_Data, SEIR_Const, fit,
                         import_scenario = import_ont_minus, rw = TRUE)
out_ont <- get_out_scenario_relache(scenario_ont, max_time)

# ---- same as ontario after March 8th----
import_ont_0308 <- rep(0, niter)
import_ont_0308[1] <- 0
qc$travel_ont_0308_toplot <- qc$nc_travel_symp
qc$travel_ont_0308 <- qc$nc_travel_symp_minus
qc$travel_ont_0308[qc$date > as.Date("2020-03-08") & qc$date < min(c(max(qc$date), max(ont$date)))] <- ont$nc_travel_minus[ont$date > as.Date("2020-03-08") & ont$date < min(c(max(qc$date), max(ont$date)))]
qc$travel_ont_0308_toplot[qc$date > as.Date("2020-03-08") & qc$date < min(c(max(qc$date), max(ont$date)))] <- ont$nc_travel[ont$date > as.Date("2020-03-08") & ont$date < min(c(max(qc$date), max(ont$date)))]
import_ont_0308[qc$time / dt] <- qc$travel_ont_0308
scenario_ont_0308 <- run_scenario_relache(SEIR_Data, SEIR_Const, fit,
                         import_scenario = import_ont_0308, rw = TRUE)
out_ont_0308 <- get_out_scenario_relache(scenario_ont_0308, max_time)

# ---- no more importation after March 8th ----
import_0308 <- rep(0, niter)
import_0308[1] <- 0
import_0308[qc$time / dt] <- qc$nc_travel_symp_arrival_March8_m
scenario_0308 <- run_scenario_relache(SEIR_Data, SEIR_Const, fit,
                         import_scenario = import_0308, rw = TRUE)
out_0308 <- get_out_scenario_relache(scenario_0308, max_time)

# ---- no more importation after March 14th ----
import_0315 <- rep(0, niter)
import_0315[1] <- 0
import_0315[qc$time / dt] <- qc$nc_travel_symp_arrival_March15_m
scenario_0315 <- run_scenario_relache(SEIR_Data, SEIR_Const, fit,
                         import_scenario = import_0315, rw = TRUE)
out_0315 <- get_out_scenario_relache(scenario_0315, max_time)

# ---- impact ----
get_out_impact(scenario, scenario_0315, max_time)
get_out_impact(scenario, scenario_0308, max_time)
get_out_impact(scenario, scenario_ont_0308, max_time)
get_out_impact(scenario, scenario_ont, max_time)

# ---- graphs ----
time_date <- seq(start_date + 1, as.Date(start_date + max_time), by = "day")
  x_lim <- c(start_date, today + 1)
  begin_date <- lubridate::ceiling_date(start_date, unit = "month") - 15
  first_date <- lubridate::floor_date(seq(as.Date("2020-02-15"), x_lim[2], by = "month"), unit = "month") + 14
  mid_date <- lubridate::ceiling_date(seq(as.Date("2020-02-15"), x_lim[2], by = "month"), unit = "month")
  sel_date <- sort(c(begin_date, first_date, mid_date))

  x_lim1 <- c(start_date, as.Date("2020-03-29"))
  mid_date1 <- lubridate::ceiling_date(seq(x_lim1[1], x_lim1[2], by = "week"), unit = "week")
  sel_date1 <- sort(mid_date1)

# Cumulative graphs
fam <- "Source Sans Pro"
# extrafont::font_import(pattern = fam, prompt = F)
# extrafont::loadfonts(quiet = T)

png(paste("./fig/", "CMAJ_F1", '.png', sep = ""), width = 8.5, height = 5, unit = "in", res = 320)
par(mfrow = c(1, 1), family = fam, yaxs = "i", mar = c(4, 4.5, 2, 0))

nf <- layout(matrix(c(1,2), ncol = 2, nrow = 1),
             widths = c(lcm(9), lcm(12)), heights = lcm(10),
             respect = TRUE)

# plotting scenarios
y_lim1 <- c(0, 250)
plot(qc$nc_travel_symp  ~ qc$date, type = 'l', col = "#678096", lwd = 1,
     xlim = x_lim1,
     ylim = y_lim1,
     ylab = "", xlab = "", axes = FALSE)
mtext("A) Observed and counterfactual scenarios", side = 3, line = 1, cex = 1)
polygon(x = c(qc$date, rev(qc$date)),
          y = c(rep(0, nrow(qc)), rev(qc$nc_travel_symp)),
          col = transparent_col(color = "#678096", percent = 80), border = NA)
abline(h = seq(0, y_lim1[2], by = 25), col = "lightgray", lwd = 0.5)
axis.Date(1, at = sel_date1, las = 2, format = "%b %d", cex.axis = 0.9)
axis(2, at = seq(0, y_lim1[2], by = 25), labels = seq(0, y_lim1[2], by = 25), las = 1)
mtext("Daily imported cases by symptoms onset date", side = 2, line = 3, cex = 1)

lines(I(qc$nc_travel_symp + 1) ~ qc$date, col = "#678096", lwd = 2)
lines(I(qc$nc_travel_symp_arrival_March8 - 0.9) ~ qc$date, col = "antiquewhite3", lwd = 2)
lines(ont$nc_travel ~ ont$date, col = "#A12A19", lwd = 2)
lines(I(qc$travel_ont_0308_toplot - 1.75) ~ qc$date, col = "#979461", lwd = 2, lty = 5)


legend("topleft", legend = c("Québec's imported cases (observed)", "Scenario: Ontario's imported cases",
                             "Scenario: No importations after March 8th", "Scenario: Same as Ontario after March 8th"),
       col = c(transparent_col(color = "#678096", percent = 70),
               "#A12A19", "antiquewhite4", "#979461"),
       pch = c(15, NA, NA, NA), pt.cex = c(1.5, NA, NA, NA),
       lty = c(NA, 1, 1, 5), lwd = c(NA, 2, 2, 2), bty = 'n', cex = 0.75)


# plotting results
y_lim <- c(0, ceiling(max(out$cum_nh$uci) / 1000) * 1000)
pt <- max(out$cum_nh$median)
pt_ont_0308 <- max(out_ont_0308$cum_nh$median)
pt_0308 <- max(out_0308$cum_nh$median)
pt_ont <- max(out_ont$cum_nh$median)

plot(out$cum_nh$median  ~ time_date, type = 'l', col = "#678096", lwd = 2, xlim = c(x_lim[1], x_lim[2] + 60), ylim = y_lim,
     ylab = "", xlab = "", axes = FALSE)
mtext("B) Impact of scenarios on hospitalisations", side = 3, line = 1, cex = 1)
#abline(h = seq(0, y_lim[2], by = 1000), col = "lightgray", lwd = 0.5)
segments(x0 = 0, y0 = seq(0, y_lim[2], by = 1000), x1 = x_lim[2], y1 = seq(0, y_lim[2], by = 1000), col = "lightgray", lwd = 0.5)
polygon(x = c(time_date, rev(time_date)),
          y = c(out$cum_nh$lci, rev(out$cum_nh$uci)),
          col = transparent_col(color = "#678096", percent = 50), border = NA)
polygon(x = c(time_date, rev(time_date)),
          y = c(out_ont_0308$cum_nh$lci, rev(out_ont_0308$cum_nh$uci)),
          col = transparent_col(color = "#979461", percent = 70), border = NA)
polygon(x = c(time_date, rev(time_date)),
          y = c(out_0308$cum_nh$lci, rev(out_0308$cum_nh$uci)),
          col = transparent_col(color = "antiquewhite2", percent = 10), border = NA)
polygon(x = c(time_date, rev(time_date)),
          y = c(out_ont$cum_nh$lci, rev(out_ont$cum_nh$uci)),
          col = transparent_col(color = "#A12A19", percent = 60), border = NA)

lines(out_ont_0308$cum_nh$median  ~ time_date, col = '#979461', lwd = 2, lty = 1)
lines(out_0308$cum_nh$median  ~ time_date, col = 'antiquewhite3', lwd = 2, lty = 1)
lines(out_ont$cum_nh$median  ~ time_date, col = '#A12A19', lwd = 2, lty = 1)
lines(out$cum_nh$median  ~ time_date, col = "#678096", lwd = 2)
points(x = x_lim[2], y = pt, pch = 21, cex = 1, col = "grey75", bg = "#678096")
points(x = x_lim[2], y = pt_ont_0308, pch = 21, cex = 1, col = "grey75", bg = "#979461")
points(x = x_lim[2], y = pt_0308, pch = 21, cex = 1, col = "grey75", bg = "antiquewhite2")
points(x = x_lim[2], y = pt_ont, pch = 21, cex = 1, col = "grey75", bg = "#A12A19")
text(x = x_lim[2] + 1, y = pt, labels = "Québec's imported cases \n (observed scenario)", pos = 4, cex = 0.5)
text(x = x_lim[2] + 1, y = pt_ont + 60, labels = "If same as Ontario", pos = 4, cex = 0.5)
text(x = x_lim[2] + 1, y = pt_ont_0308 - 60, labels = "If same as Ontario \n after March 8th", pos = 4, cex = 0.5)
text(x = x_lim[2] + 1, y = pt_0308, labels = "If no importations \n after March 8th", pos = 4, cex = 0.5)
  axis.Date(1, at = sel_date, las = 2, format = "%b %d", cex.axis = 0.9)
  axis(2, at = seq(0, y_lim[2], by = 1000), labels = seq(0, y_lim[2], by = 1000), las = 1)
  mtext("Cumulative number of hospitalisations", side = 2, line = 3, cex = 1)
dev.off()

