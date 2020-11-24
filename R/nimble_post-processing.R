is_wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol

get_outcomes <- function(fit, SEIR_Data, SEIR_Const) {

  trace <- as.data.frame(do.call(rbind, fit))
  if (any(is.na(trace))) { print("NA in some outcomes, using na.rm = TRUE") }

  # new cases
  nom <- "prd_nc"
    med <- apply(trace[, grepl(nom, names(trace))], MARGIN = 2, median, na.rm = TRUE)
    lci <- apply(trace[, grepl(nom, names(trace))], MARGIN = 2, function(x) quantile(x, probs = c(0.025), na.rm = TRUE))
    uci <- apply(trace[, grepl(nom, names(trace))], MARGIN = 2, function(x) quantile(x, probs = c(0.975), na.rm = TRUE))
    l50 <- apply(trace[, grepl(nom, names(trace))], MARGIN = 2, function(x) quantile(x, probs = c(0.25), na.rm = TRUE))
    u50 <- apply(trace[, grepl(nom, names(trace))], MARGIN = 2, function(x) quantile(x, probs = c(0.75), na.rm = TRUE))
    nc <- data.frame(median = med, lci = lci, uci = uci, l50 = l50, u50 = u50)

  # number daily new hospitalization
  nom <- "prd_nh"
    med <- apply(trace[, grepl(nom, names(trace))], MARGIN = 2, median, na.rm = TRUE)
    lci <- apply(trace[, grepl(nom, names(trace))], MARGIN = 2, function(x) quantile(x, probs = c(0.025), na.rm = TRUE))
    uci <- apply(trace[, grepl(nom, names(trace))], MARGIN = 2, function(x) quantile(x, probs = c(0.975), na.rm = TRUE))
    l50 <- apply(trace[, grepl(nom, names(trace))], MARGIN = 2, function(x) quantile(x, probs = c(0.25), na.rm = TRUE))
    u50 <- apply(trace[, grepl(nom, names(trace))], MARGIN = 2, function(x) quantile(x, probs = c(0.75), na.rm = TRUE))
    nh <- data.frame(median = med, lci = lci, uci = uci, l50 = l50, u50 = u50)

  # deaths
  nom <- "prd_dx"
    med <- apply(trace[, grepl(nom, names(trace))], MARGIN = 2, median, na.rm = TRUE)
    lci <- apply(trace[, grepl(nom, names(trace))], MARGIN = 2, function(x) quantile(x, probs = c(0.025), na.rm = TRUE))
    uci <- apply(trace[, grepl(nom, names(trace))], MARGIN = 2, function(x) quantile(x, probs = c(0.975), na.rm = TRUE))
    l50 <- apply(trace[, grepl(nom, names(trace))], MARGIN = 2, function(x) quantile(x, probs = c(0.25), na.rm = TRUE))
    u50 <- apply(trace[, grepl(nom, names(trace))], MARGIN = 2, function(x) quantile(x, probs = c(0.75), na.rm = TRUE))
    dx <- data.frame(median = med, lci = lci, uci = uci, l50 = l50, u50 = u50)

  # total recovered
  nom <- "prd_rx"
    med <- apply(trace[, grepl(nom, names(trace))], MARGIN = 2, median, na.rm = TRUE)
    lci <- apply(trace[, grepl(nom, names(trace))], MARGIN = 2, function(x) quantile(x, probs = c(0.025), na.rm = TRUE))
    uci <- apply(trace[, grepl(nom, names(trace))], MARGIN = 2, function(x) quantile(x, probs = c(0.975), na.rm = TRUE))
    l50 <- apply(trace[, grepl(nom, names(trace))], MARGIN = 2, function(x) quantile(x, probs = c(0.25), na.rm = TRUE))
    u50 <- apply(trace[, grepl(nom, names(trace))], MARGIN = 2, function(x) quantile(x, probs = c(0.75), na.rm = TRUE))
    rx <- data.frame(median = med, lci = lci, uci = uci, l50 = l50, u50 = u50)

  # calculating R0
  nom <- "beta_iter"
  trace_beta <- trace[, grepl(nom, names(trace))]
  prop_s <-     1 - trace[, grepl("prd_rx", names(trace))] / SEIR_Data[["init"]][1]
  # if we are fitting SEIR
  if (length(SEIR_Const$rho) == 0) {
      duration_inf <- SEIR_Const[["sigma_c"]]
  }
  # if we are fitting SEAIR
  if (length(SEIR_Const$rho) == 1) {
      duration_inf <- 1 / (1 / SEIR_Const[["eps"]] + 1 / SEIR_Const[["sigma_c"]])
  }
  trace_r0 <-  (SEIR_Const[["p_symp"]] * trace_beta / duration_inf +
              (1 - SEIR_Const[["p_symp"]]) * trace_beta * SEIR_Const[["infa"]] / duration_inf) * prop_s

    r0_med <- apply(trace_r0, MARGIN = 2, median)
    r0_lci <- apply(trace_r0, MARGIN = 2, function(x) quantile(x, probs = c(0.025)))
    r0_uci <- apply(trace_r0, MARGIN = 2, function(x) quantile(x, probs = c(0.975)))
    r0_l50 <- apply(trace_r0, MARGIN = 2, function(x) quantile(x, probs = c(0.25)))
    r0_u50 <- apply(trace_r0, MARGIN = 2, function(x) quantile(x, probs = c(0.75)))

  val <- list(nc = nc, nh = nh, dx = dx, rx = rx,
              r0_med = r0_med, r0_lci = r0_lci, r0_uci = r0_uci,
              r0_l50 = r0_l50, r0_u50 = r0_u50)
}


plot_forecast <- function(out, SEIR_Data, SEIR_Const,
                          start_date = as.Date("2020-01-22"),
                          today = today, forecast = 7,
                          file_name = "Quebec", plot_title,
                          save = FALSE, watermark = TRUE, plot_max = FALSE,
                          width = 4, height = 6.7, unit = 'in', quality = 75,
                          res = 320, report = FALSE) {

  old_locale <- Sys.getlocale("LC_TIME")
  Sys.setlocale("LC_TIME", "fr_FR")

  if (report){
    fam <- "Raleway"
    load_font(fam)
  } else {
    fam <- "Arial"
  }

  plot_forecast_wrapper(out, SEIR_Data, SEIR_Const,
                        start_date, today, forecast,
                        plot_title, watermark, plot_max, fam)

  if (save == TRUE) {
  jpeg(paste("./fig/", file_name, '.jpg', sep = ""), width = width, height = height, unit = unit, quality = quality, res = res)
  plot_forecast_wrapper(out, SEIR_Data, SEIR_Const,
                        start_date, today, forecast,
                        plot_title, watermark, plot_max, fam)
  dev.off()
  }
  Sys.setlocale("LC_TIME", old_locale)
}

plot_forecast_wrapper <- function(out, SEIR_Data, SEIR_Const,
                                  start_date, today, forecast,
                                  plot_title, watermark, plot_max, fam) {
  nc <- out[["nc"]]
  #hx <- out[["hx"]]
  nh <- out[["nh"]]
  dx <- out[["dx"]]
  r0_med <- out[["r0_med"]]
  r0_lci <- out[["r0_lci"]]
  r0_uci <- out[["r0_uci"]]
  r0_l50 <- out[["r0_l50"]]
  r0_u50 <- out[["r0_u50"]]

  time_date <- start_date + (SEIR_Const$iter_time * SEIR_Const$dt) - SEIR_Const$dt
  data_nc <- data.frame(nc = SEIR_Data$nc)
  data_nc$time_date <- start_date + SEIR_Const$y_time * SEIR_Const$dt
  data_nc <- subset(data_nc, !is.na(nc))

  # x limit axis
  x_lim <- c(start_date, today + forecast)
  begin_date <- lubridate::ceiling_date(start_date, unit = "month")
  first_date <- lubridate::floor_date(seq(as.Date("2020-03-01"), x_lim[2], by = "month"), unit = "month") + 14
  mid_date <- lubridate::ceiling_date(seq(as.Date("2020-03-01"), x_lim[2], by = "month"), unit = "month")
  sel_date <- sort(c(begin_date, first_date, mid_date))
  if (max(sel_date) > (x_lim[2] - 5)) { sel_date[length(sel_date)] <- x_lim[2] }


  par(mfrow = c(4, 1), mar = c(2, 4, 1, 1), family = fam)
  if (all(is.na(data_nc$nc))) {
    #y_lim_nc <- c(0, nc$uci[time_date %in% x_lim[2]])
    y_lim_nc <- c(0, max(nc$uci[time_date < (x_lim[2] - 1)]))
    y_lim_nh <- y_lim_nc / 5
    y_lim_dx <- y_lim_nc / 20
  } else {
    y_lim_nc <- c(0, max(data$nc* 1.5, na.rm = TRUE))
    y_lim_nh <- y_lim_nc / 5
    y_lim_dx <- y_lim_nc / 20
  }
  if (plot_max) {
    y_lim_nc <- c(0, max(nc$uci[time_date < (x_lim[2] - 1)]))
    y_lim_nh <- c(0, max(nh$uci[time_date < (x_lim[2] - 1)]))
    y_lim_dx <- c(0, max(dx$uci[time_date < (x_lim[2] - 1)]))
    }

  plot(nc$median ~ time_date, type = 'n', xlim = x_lim, ylim = y_lim_nc,
       ylab = "", xlab = "", main = plot_title, axes = FALSE)
  axis.Date(1, at = sel_date)
  axis(2, at = seq(0, y_lim_nc[2], by = 200), labels = seq(0, y_lim_nc[2], by = 200), las = 1)
  polygon(x = c(time_date, rev(time_date)),
          y = c(nc$lci, rev(nc$uci)),
          col = rgb(217, 168, 33, 150, max = 255), border = NA)
  polygon(x = c(time_date, rev(time_date)),
          y = c(nc$l50, rev(nc$u50)),
          col = rgb(217, 168, 33, 150, max = 255), border = NA)
  lines(nc$median ~ time_date, col = rgb(217, 168, 33, 255, max = 255), lwd = 2)
  points(data_nc$nc ~ data_nc$time_date, pch = 16, col = rgb(100, 100, 100, 200, max = 255))
  abline(v = max(data_nc$time_date + 1), col = "grey50", lty = 2)
  text(x = max(data_nc$time_date + 1.1), y = 0, "Scénario", pos = 4, srt = 0,
       col = rgb(100, 100, 100, 176, max = 255))
  text(x = min(time_date), y = y_lim_nc[2] * 0.95, "A) Nombre de cas confirmés (par jour)", pos = 4)
  mtext("Cas confirmés", side = 2, line = 3, cex = 0.65)
  legend("bottomleft", legend = c("IQR", "95% CrI"),
         pch = 15, pt.cex = 2, bty = 'n',
         col = c(rgb(217, 168, 33, 220, max = 255), rgb(217, 168, 33, 140, max = 255)))
  if (watermark) {
    text(x = grconvertX(0.5, from = "npc"),  # align to center of plot X axis
         y = grconvertY(0.5, from = "npc"), # align to center of plot Y axis
         labels = "Préliminaire", # our watermark
         cex = 2, font = 2, # large, bold font - hard to miss
         col = rgb(1, 0, 0, 0.4), # translucent (0.2 = 20%) red color
         srt = 35) # srt = angle of text: 45 degree angle to X axis
  }

  data_nh <- data.frame(nh = SEIR_Data$nh)
  data_nh$time_date <- start_date + SEIR_Const$y_time * SEIR_Const$dt
  data_nh <- subset(data_nh, !is.na(nh))
  plot(nh$median ~ time_date, type = 'n', xlim = x_lim, ylim = y_lim_nh,
       ylab = "", xlab = "", axes = FALSE)
  axis.Date(1, at = sel_date)
  axis(2, at = seq(0, y_lim_nh[2], by = 20), labels = seq(0, y_lim_nh[2], by = 20), las = 1)
  polygon(x = c(time_date, rev(time_date)),
          y = c(nh$lci, rev(nh$uci)),
          col = rgb(245, 160, 142, 125, max = 255), border = NA)
  polygon(x = c(time_date, rev(time_date)),
          y = c(nh$l50, rev(nh$u50)),
          col = rgb(245, 160, 142, 175, max = 255), border = NA)
  lines(nh$median ~ time_date, col = rgb(210, 140, 130, 255, max = 255), lwd = 2)
  points(data_nh$nh ~ data_nh$time_date, pch = 16,
         col = rgb(100, 100, 100, 200, max = 255))
  abline(v = max(data_nh$time_date + 1), col = "grey50", lty = 2)
  text(x = max(data_nh$time_date + 1.1), y = 0, "Scénario", pos = 4, srt = 0,
       col = rgb(100, 100, 100, 176, max = 255))
  text(x = min(time_date), y = y_lim_nh[2] * 0.95, "B) Nouvelles hospitalisations (par jour)", pos = 4)
  mtext("Hospitalisations", side = 2, line = 3, cex = 0.65)
  if (watermark) {
    text(x = grconvertX(0.5, from = "npc"),  # align to center of plot X axis
         y = grconvertY(0.5, from = "npc"), # align to center of plot Y axis
         labels = "Préliminaire", # our watermark
         cex = 2, font = 2, # large, bold font - hard to miss
         col = rgb(1, 0, 0, 0.4), # translucent (0.2 = 20%) red color
         srt = 35) # srt = angle of text: 45 degree angle to X axis
  }

  data_dx <- data.frame(dx = SEIR_Data$dx)
  data_dx$time_date <- start_date + SEIR_Const$y_time * SEIR_Const$dt
  data_dx <- subset(data_dx, !is.na(dx))
  #data_dx$cum_dx <- cumsum(data_dx$dx)
  #data_dx <- subset(data_dx, cum_dx > 0)
  plot(dx$median ~ time_date, type = 'n', xlim = x_lim, ylim = y_lim_dx,
       ylab = "", xlab = "", axes = FALSE)
  axis.Date(1, at = sel_date)
  axis(2, at = seq(0, y_lim_dx[2], by = 5), labels = seq(0, y_lim_dx[2], by = 5), las = 1)
  polygon(x = c(time_date, rev(time_date)),
          y = c(dx$lci, rev(dx$uci)),
          col = rgb(139, 58, 58, 120, max = 255), border = NA)
  polygon(x = c(time_date, rev(time_date)),
          y = c(dx$l50, rev(dx$u50)),
          col = rgb(139, 58, 58, 120, max = 255), border = NA)
  lines(dx$median ~ time_date,col = rgb(139, 58, 58, 255, max = 255), lwd = 2)
  abline(v = max(data_dx$time_date + 1), col = "grey50", lty = 2)
  text(x = max(data_dx$time_date + 1.1), y = 0, "Scénario", pos = 4, srt = 0,
       col = rgb(100, 100, 100, 176, max = 255))
  points(data_dx$dx ~ data_dx$time_date, pch = 16,
         col = rgb(100, 100, 100, 200, max = 255))
  text(x = min(time_date), y = y_lim_dx[2] * 0.95, "C) Nombre de décès (par jour)", pos = 4)
  mtext("Décès", side = 2, line = 3, cex = 0.65)
  if (watermark) {
    text(x = grconvertX(0.5, from = "npc"),  # align to center of plot X axis
         y = grconvertY(0.5, from = "npc"), # align to center of plot Y axis
         labels = "Préliminaire", # our watermark
         cex = 2, font = 2, # large, bold font - hard to miss
         col = rgb(1, 0, 0, 0.4), # translucent (0.2 = 20%) red color
         srt = 35) # srt = angle of text: 45 degree angle to X axis
  }

  y_lim_r <- c(0, max(5, max(r0_uci) * 1.2))
  plot(r0_med ~ time_date, type = 'n', xlim = x_lim, ylim = y_lim_r,
       ylab = "", xlab = "", axes = FALSE)
  axis.Date(1, at = sel_date)
  axis(2, at = c(0, 0.5, seq(1, y_lim_r[2], by = 1)), labels = c(0, 0.5, seq(1, y_lim_r[2], by = 1)), las = 1)
  abline(h = c(0, 0.5, 2, 3, 4), col = "lightgray", lwd = 0.5)
  polygon(x = c(time_date, rev(time_date)),
          y = c(r0_lci, rev(r0_uci)),
          col = rgb(115, 160, 140, 150, max = 255), border = NA)
  polygon(x = c(time_date, rev(time_date)),
          y = c(r0_l50, rev(r0_u50)),
          col = rgb(115, 160, 140, 150, max = 255), border = NA)
  lines(r0_med ~ time_date, col = rgb(115, 160, 140, 255, max = 255), lwd = 2)
  abline(h = 1, col = "grey25", lty = 3)
  abline(v = today + 1, col = "grey50", lty = 2)
  text(x = today + 1.1, y = 0, "Scénario", pos = 4, srt = 0,
       col = rgb(100, 100, 100, 176, max = 255))
  text(x = min(time_date), y = y_lim_r[2] * 0.95, "D) Taux de reproduction - R(t)", pos = 4)
  mtext("Taux", side = 2, line = 3, cex = 0.65)
  int_1 <- start_date + SEIR_Const$intervention1_ind * SEIR_Const$dt - SEIR_Const$dt
  int_2 <- start_date + SEIR_Const$intervention2_ind * SEIR_Const$dt - SEIR_Const$dt
  y0 <- r0_uci[which(time_date == int_1)] + 0.25
  y1 <- r0_uci[which(time_date == int_1)] + 1
  arrows(x0 = int_1, y0,
         x1 = int_1, y1,
         code = 1, length = 0.05, col = "grey50")
  text(x = int_1 + 0.5, y = y1, pos = 4, "Premières mesures", col = "grey50")
  y0 <- r0_uci[which(time_date == int_2)] + 0.25
  y1 <- r0_uci[which(time_date == int_2)] + 1
  arrows(x0 = int_2, y0,
         x1 = int_2, y1,
         code = 1, length = 0.05, col = "grey50")
  text(x = int_2 + 0.5, y = y1, pos = 4, "Loi services essentiels", col = "grey50")
  if (watermark) {
    text(x = grconvertX(0.5, from = "npc"),  # align to center of plot X axis
         y = grconvertY(0.5, from = "npc"), # align to center of plot Y axis
         labels = "Préliminaire", # our watermark
         cex = 2, font = 2, # large, bold font - hard to miss
         col = rgb(1, 0, 0, 0.4), # translucent (0.2 = 20%) red color
         srt = 35) # srt = angle of text: 45 degree angle to X axis
  }

}


# plot hospit only

plot_hospit <- function(out, SEIR_Data, SEIR_Const,
                          start_date = as.Date("2020-01-22"),
                          today = today, forecast = 7,
                          file_name = "Quebec",
                          plot_title = "Québec",
                          ylim_hosp = c(0, 140),
                          save = FALSE, watermark = TRUE,
                          width = 5, height = 6.5, unit = 'in', quality = 75,
                          res = 320, report = FALSE) {

  old_locale <- Sys.getlocale("LC_TIME")
  Sys.setlocale("LC_TIME", "fr_FR")

  if (report){
    fam <- "Raleway"
    load_font(fam)
  } else {
    fam <- "Arial"
  }

  plot_hospit_wrapper(out, SEIR_Data, SEIR_Const,
                        start_date, today, forecast,
                        plot_title, ylim_hosp, watermark, fam)

  if (save == TRUE) {
  #jpeg(paste("./fig/", file_name, '.jpg', sep = ""), width = width, height = height, unit = unit, quality = quality, res = res)
  png(paste("./fig/", file_name, '.png', sep = ""), width = width, height = height, unit = unit, res = res)
    plot_hospit_wrapper(out, SEIR_Data, SEIR_Const,
                        start_date, today, forecast,
                        plot_title, ylim_hosp, watermark, fam)
  dev.off()
  }

  Sys.setlocale("LC_TIME", old_locale)

}

plot_hospit_wrapper <- function(out, SEIR_Data, SEIR_Const,
                                  start_date, today, forecast,
                                  plot_title, ylim_hosp, watermark, fam) {

  nh <- out[["nh"]]
  r0_med <- out[["r0_med"]]
  r0_lci <- out[["r0_lci"]]
  r0_uci <- out[["r0_uci"]]
  r0_l50 <- out[["r0_l50"]]
  r0_u50 <- out[["r0_u50"]]

  time_date <- start_date + (SEIR_Const$iter_time * SEIR_Const$dt) - SEIR_Const$dt
  data_nc <- data.frame(nc = SEIR_Data$nc)
  data_nc$time_date <- start_date + SEIR_Const$y_time * SEIR_Const$dt
  data_nc <- subset(data_nc, !is.na(nc))


  # x limit axis
  x_lim <- c(start_date, today + forecast)
  begin_date <- lubridate::ceiling_date(start_date, unit = "month")
  first_date <- lubridate::floor_date(seq(as.Date("2020-03-01"), x_lim[2], by = "month"), unit = "month") + 14
  mid_date <- lubridate::ceiling_date(seq(as.Date("2020-03-01"), x_lim[2], by = "month"), unit = "month")
  sel_date <- sort(c(begin_date, first_date, mid_date))
  if (max(sel_date) > (x_lim[2] - 5)) { sel_date[length(sel_date)] <- x_lim[2] }

  # first graph
  par(mfrow = c(2, 1), mar = c(4, 4, 1, 1), family = fam)
  data_nh <- data.frame(nh = SEIR_Data$nh)
  data_nh$time_date <- start_date + SEIR_Const$y_time * SEIR_Const$dt
  data_nh <- subset(data_nh, !is.na(nh))
  plot(nh$median ~ time_date, type = 'n', xlim = x_lim, ylim = ylim_hosp,
       ylab = "", xlab = "", main = plot_title, axes = FALSE)
  axis.Date(1, at = sel_date, las = 2, format = "%b %d", cex.axis = 0.9)
  axis(2, at = seq(0, ylim_hosp[2], by = 20), labels = seq(0, ylim_hosp[2], by = 20), las = 1)
  abline(h = seq(0, ylim_hosp[2], by = 20), col = "lightgray", lwd = 0.5)
  polygon(x = c(time_date, rev(time_date)),
          y = c(nh$lci, rev(nh$uci)),
          col = rgb(245, 160, 142, 125, max = 255), border = NA)
  polygon(x = c(time_date, rev(time_date)),
          y = c(nh$l50, rev(nh$u50)),
          col = rgb(245, 160, 142, 175, max = 255), border = NA)
  lines(nh$median ~ time_date, col = rgb(210, 140, 130, 255, max = 255), lwd = 2)
  points(data_nh$nh ~ data_nh$time_date, pch = 16, cex = 0.6,
         col = rgb(100, 100, 100, 200, max = 255))
  abline(v = max(data_nh$time_date + 1), col = "grey50", lty = 2)
  text(x = max(data_nh$time_date + 1.1), y = 0, "Prédiction", pos = 4, srt = 90,
       col = rgb(100, 100, 100, 176, max = 255), cex = 0.7)
  text(x = min(time_date), y = ylim_hosp[2] * 0.95, "A) Nouvelles hospitalisations (par jour)", pos = 4)
  mtext("Hospitalisations", side = 2, line = 2.5, cex = 1)
  if (watermark) {
    text(x = grconvertX(0.5, from = "npc"),  # align to center of plot X axis
         y = grconvertY(0.5, from = "npc"), # align to center of plot Y axis
         labels = "Préliminaire", # our watermark
         cex = 1.5, font = 2, # large, bold font - hard to miss
         col = rgb(1, 0, 0, 0.3), # translucent (0.2 = 20%) red color
         srt = 35) # srt = angle of text: 45 degree angle to X axis
  }
  legend(x = time_date[1], y = ylim_hosp[2] * 0.9, legend = c("ICr 50%", "ICr 95%"),
         pch = 15, pt.cex = 1.5, cex = 0.8, bty = 'n',
         col = c(rgb(245, 160, 142, 200, max = 255), rgb(245, 160, 142, 125, max = 255)))

  y_lim_r <- c(0, max(5, max(r0_uci) * 1.2))
  plot(r0_med ~ time_date, type = 'n', xlim = x_lim, ylim = y_lim_r, axes = FALSE,
       xlab = "", ylab = "")
  axis.Date(1, at = sel_date, las = 2, format = "%b %d", cex.axis = 0.9)
  axis(2, at = c(0, 0.5, seq(1, y_lim_r[2], by = 1)), labels = c(0, 0.5, seq(1, y_lim_r[2], by = 1)), las = 1)
  abline(h = c(0, 0.5, 2, 3, 4), col = "lightgray", lwd = 0.5)
  polygon(x = c(time_date, rev(time_date)),
          y = c(r0_lci, rev(r0_uci)),
          col = rgb(115, 160, 140, 150, max = 255), border = NA)
  polygon(x = c(time_date, rev(time_date)),
          y = c(r0_l50, rev(r0_u50)),
          col = rgb(115, 160, 140, 150, max = 255), border = NA)
  lines(r0_med ~ time_date, col = rgb(115, 160, 140, 255, max = 255), lwd = 2)
  abline(h = 1, col = "grey25", lty = 3)
  abline(v = max(data_nh$time_date + 1), col = "grey50", lty = 2)
  text(x = max(data_nh$time_date) + 1.1, y = 0, "Prédiction", pos = 4, srt = 90,
       col = rgb(100, 100, 100, 176, max = 255), cex = 0.7)
  text(x = min(time_date), y = y_lim_r[2] * 0.95, "B) Taux de reproduction - R(t)", pos = 4)
  mtext("Taux", side = 2, line = 2.5, cex = 1)
  int_1 <- start_date + SEIR_Const$intervention1_ind * SEIR_Const$dt - SEIR_Const$dt
  int_2 <- start_date + SEIR_Const$intervention2_ind * SEIR_Const$dt - SEIR_Const$dt
  y0 <- r0_uci[which(time_date == int_1)] + 0.25
  y1 <- r0_uci[which(time_date == int_1)] + 1
  arrows(x0 = int_1, y0,
         x1 = int_1, y1,
         code = 1, length = 0.05, col = "grey50")
  text(x = int_1 + 0.5, y = y1 + 0.25, pos = 4, "Premières mesures", col = "grey50")
  if (watermark) {
    text(x = grconvertX(0.5, from = "npc"),  # align to center of plot X axis
         y = grconvertY(0.5, from = "npc"), # align to center of plot Y axis
         labels = "Préliminaire", # our watermark
         cex = 1.5, font = 2, # large, bold font - hard to miss
         col = rgb(1, 0, 0, 0.3), # translucent (0.2 = 20%) red color
         srt = 35) # srt = angle of text: 45 degree angle to X axis
  }
    legend(x = time_date[1], y = 1, legend = c("ICr 50%", "ICr 95%"),
         pch = 15, pt.cex = 1.5, cex = 0.8, bty = 'n',
         col = c(rgb(115, 160, 140, 200, max = 255), rgb(115, 160, 140, 125, max = 255)))

}

# Font loading wrapper
load_font <- function(pattern) {
  if (pattern %in% names(pdfFonts())) return
  if (!(pattern %in% extrafont::fonts())) extrafont::font_import(pattern = pattern, prompt = F)
  extrafont::loadfonts(quiet = T)
}
