
transparent_col <- function(color, percent = 50, name = NULL) {
## Mark Gardener 2015
## www.dataanalytics.org.uk
## Get RGB values for named color
rgb.val <- col2rgb(color)
## Make new color using input color as base and alpha set by transparency
t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
             max = 255,
             alpha = (100 - percent) * 255 / 100,
             names = name)
## Save the color
return(t.col)
}

run_scenario_relache <- function(SEIR_Data, SEIR_Const, fit,
                         import_scenario = SEIR_Data$num_import,
                         rw = FALSE) {

trace <- as.data.frame(do.call(rbind, fit))
prd_nc <- matrix(NA, nrow = nrow(trace), ncol = SEIR_Const$Z)
prd_nh <- matrix(NA, nrow = nrow(trace), ncol = SEIR_Const$Z)
prd_r0 <- matrix(NA, nrow = nrow(trace), ncol = SEIR_Const$Z)
prd_rx <- matrix(NA, nrow = nrow(trace), ncol = SEIR_Const$Z)
prd_h <- matrix(NA, nrow = nrow(trace), ncol = SEIR_Const$Z)
prd_dx <- matrix(NA, nrow = nrow(trace), ncol = SEIR_Const$Z)

cum_nc <- cum_nh <- cum_dx <- matrix(NA, nrow = nrow(trace), ncol = SEIR_Const$Z)

for (i in 1:nrow(trace)) {

if (rw == FALSE) {
# we fill in beta_exp
beta_exp <- rep(exp(trace$b[i]), SEIR_Const$niter)
for (j in 1:SEIR_Const$n_beta) {
  nom_reduction <- paste("reduction[", j, "]", sep = '')
  beta_exp[SEIR_Const$beta_ind[j]:(SEIR_Const$beta_ind[j + 1] - 1)] <- exp(trace$b[i]) * trace[i, names(trace) %in% nom_reduction]
}
  nom_reduction <- paste("reduction[", SEIR_Const$n_beta, "]", sep = '')
  beta_exp[SEIR_Const$beta_ind[SEIR_Const$n_beta]:SEIR_Const$niter] <- exp(trace$b[i]) *  trace[i, names(trace) %in% nom_reduction]
}

#if we are doing random walk
if (rw == TRUE) {
  beta_exp <- rep(exp(trace[i, names(trace) %in% "b[1]"]), SEIR_Const$niter)
  for (j in 2:SEIR_Const$n_beta) {
    nom_rw <- paste("b[", j, "]", sep = '')
    beta_exp[SEIR_Const$beta_ind[j - 1]:SEIR_Const$beta_ind[j]] <- exp(trace[i, names(trace) %in% nom_rw])
  }
    #nom_rw <- paste("b[", SEIR_Const$n_beta, "]", sep = '')
    #last_b_rw <- exp(trace[i, names(trace) %in% nom_rw] + rnorm(1, mean = 0, sd = trace[i, names(trace) %in% "sd_rw"]))
    beta_exp[(SEIR_Const$beta_ind[SEIR_Const$n_beta]):SEIR_Const$niter] <- trace[i, names(trace) %in% "b_last"]#last_b_rw
}


num_import_t0 <- import_scenario
if (is.null(trace$import_t0[i])) {
  num_import_t0[1] <- 1
} else {
  num_import_t0[1] <- trace$import_t0[i]
}


# SEIR model
if (length(SEIR_Const$rho) == 0) {
mod <-  seirCpp(beta = beta_exp,
              infa = SEIR_Const$infa,          # reduction of infectious if asymptomatic
              p_symp = SEIR_Const$p_symp,       # proportion of cases that are symptomatic
              eps = SEIR_Const$eps,       # rate from exposed to infectious
              p_hosp = SEIR_Const$p_hosp,      # probability of hospitalization if symptomatic
              sigma_h = SEIR_Const$sigma_h,     # rate from symptomatic to hospitalization
              sigma_c = SEIR_Const$sigma_c,     # rate from infectious to recovered
              sigma_d = SEIR_Const$sigma_d,       # rate from hospitalized to recovered (alive)
              conf_rate = SEIR_Data$conf_rate,   # confirmation rate among those with symptoms (delays of 5 days from symptom onset)
              mu = SEIR_Const$mu, # proportion dying in hospitals of covid-mortality
              num_import = num_import_t0,
              prp_detect = SEIR_Const$prp_detect,
              p_conf = SEIR_Data$p_conf,     # scaling fator for case confirmation.
              init = SEIR_Data$init,
              start = SEIR_Const$start, end = SEIR_Const$end, dt = SEIR_Const$dt)

    duration_inf <- SEIR_Const[["sigma_c"]]
}

# SEAIR model
if (length(SEIR_Const$rho) == 1) {
mod <-  seairCpp(beta = beta_exp,
              infa = SEIR_Const$infa,          # reduction of infectious if asymptomatic
              p_symp = SEIR_Const$p_symp,       # proportion of cases that are symptomatic
              rho = SEIR_Const$rho,
              eps = SEIR_Const$eps,       # rate from exposed to infectious
              p_hosp = SEIR_Const$p_hosp,      # probability of hospitalization if symptomatic
              sigma_h = SEIR_Const$sigma_h,     # rate from symptomatic to hospitalization
              sigma_c = SEIR_Const$sigma_c,     # rate from infectious to recovered
              sigma_d = SEIR_Const$sigma_d,       # rate from hospitalized to recovered (alive)
              conf_rate = SEIR_Data$conf_rate,   # confirmation rate among those with symptoms (delays of 5 days from symptom onset)
              mu = SEIR_Const$mu, # proportion dying in hospitals of covid-mortality
              num_import = num_import_t0,
              prp_detect = SEIR_Const$prp_detect,
              p_conf = SEIR_Data$p_conf,     # scaling fator for case confirmation.
              init = SEIR_Data$init,
              start = SEIR_Const$start, end = SEIR_Const$end, dt = SEIR_Const$dt)

  duration_inf <- 1 / (1 / SEIR_Const[["eps"]] + 1 / SEIR_Const[["sigma_c"]])
}

# putting everything together
  prop_s <- 1 - (mod[, "R"][SEIR_Const$iter_time] / SEIR_Data[["init"]][1])

  r_nc <- trace$theta_nc[i]
  p_nc <- r_nc / (r_nc + mod[, "NC"][SEIR_Const$iter_time] + 1e-6)
  prd_nc[i, ] <- rnbinom(length(p_nc), prob = p_nc, size = r_nc)
  cum_nc[i, ] <- cumsum(prd_nc[i, ])

  r_nh <- trace$theta_nh[i]
  p_nh <- r_nh / (r_nh + mod[, "NH"][SEIR_Const$iter_time] + 1e-6)
  prd_nh[i, ] <- rnbinom(length(p_nh), prob = p_nh, size = r_nh)
  cum_nh[i, ] <- cumsum(prd_nh[i, ])

  r_dx <- trace$theta_dx[i]
  p_dx <- r_dx / (r_dx + mod[, "DX"][SEIR_Const$iter_time] + 1e-6)
  prd_dx[i, ] <- rnbinom(length(p_dx), prob = p_dx, size = r_dx)
  cum_dx[i, ] <- cumsum(prd_dx[i, ])

  prd_r0[i, ] <- (SEIR_Const[["p_symp"]] * mod[, "beta_iter"][SEIR_Const$iter_time] / duration_inf +
                 (1 - SEIR_Const[["p_symp"]]) * mod[, "beta_iter"][SEIR_Const$iter_time] * SEIR_Const[["infa"]] / duration_inf) * prop_s
  prd_rx[i, ] <- mod[, "R"][SEIR_Const$iter_time]
  prd_h[i, ] <- mod[, "H"][SEIR_Const$iter_time]
}

out_sc_trace <- list(prd_nc = as.data.frame(prd_nc),
                     prd_nh = as.data.frame(prd_nh),
                     prd_dx = as.data.frame(prd_dx),
                     cum_nc = as.data.frame(cum_nc),
                     cum_nh = as.data.frame(cum_nh),
                     cum_dx = as.data.frame(cum_dx),
                prd_r0 = as.data.frame(prd_r0),
                prd_rx = as.data.frame(prd_rx),
                prd_h  = as.data.frame(prd_h))
return(out_sc_trace)
}

get_out_scenario_relache <- function(out_sc_trace, max_time) {

    trace_nc <- out_sc_trace[["prd_nc"]][1:max_time]
    med <- apply(trace_nc, MARGIN = 2, median)
    lci <- apply(trace_nc, MARGIN = 2, function(x) quantile(x, probs = c(0.025)))
    uci <- apply(trace_nc, MARGIN = 2, function(x) quantile(x, probs = c(0.975)))
    l50 <- apply(trace_nc, MARGIN = 2, function(x) quantile(x, probs = c(0.25)))
    u50 <- apply(trace_nc, MARGIN = 2, function(x) quantile(x, probs = c(0.75)))
    nc <- data.frame(median = med, lci = lci, uci = uci, l50 = l50, u50 = u50)

    trace_cum_nc <- out_sc_trace[["cum_nc"]][1:max_time]
    med <- apply(trace_cum_nc, MARGIN = 2, median)
    lci <- apply(trace_cum_nc, MARGIN = 2, function(x) quantile(x, probs = c(0.025)))
    uci <- apply(trace_cum_nc, MARGIN = 2, function(x) quantile(x, probs = c(0.975)))
    l50 <- apply(trace_cum_nc, MARGIN = 2, function(x) quantile(x, probs = c(0.25)))
    u50 <- apply(trace_cum_nc, MARGIN = 2, function(x) quantile(x, probs = c(0.75)))
    cum_nc <- data.frame(median = med, lci = lci, uci = uci, l50 = l50, u50 = u50)

    trace_nh <- out_sc_trace[["prd_nh"]][1:max_time]
    med <- apply(trace_nh, MARGIN = 2, median)
    lci <- apply(trace_nh, MARGIN = 2, function(x) quantile(x, probs = c(0.025)))
    uci <- apply(trace_nh, MARGIN = 2, function(x) quantile(x, probs = c(0.975)))
    l50 <- apply(trace_nh, MARGIN = 2, function(x) quantile(x, probs = c(0.25)))
    u50 <- apply(trace_nh, MARGIN = 2, function(x) quantile(x, probs = c(0.75)))
    nh <- data.frame(median = med, lci = lci, uci = uci, l50 = l50, u50 = u50)

    trace_cum_nh <- out_sc_trace[["cum_nh"]][1:max_time]
    med <- apply(trace_cum_nh, MARGIN = 2, median)
    lci <- apply(trace_cum_nh, MARGIN = 2, function(x) quantile(x, probs = c(0.025)))
    uci <- apply(trace_cum_nh, MARGIN = 2, function(x) quantile(x, probs = c(0.975)))
    l50 <- apply(trace_cum_nh, MARGIN = 2, function(x) quantile(x, probs = c(0.25)))
    u50 <- apply(trace_cum_nh, MARGIN = 2, function(x) quantile(x, probs = c(0.75)))
    cum_nh <- data.frame(median = med, lci = lci, uci = uci, l50 = l50, u50 = u50)

    trace_dx <- out_sc_trace[["prd_dx"]][1:max_time]
    med <- apply(trace_dx, MARGIN = 2, median)
    lci <- apply(trace_dx, MARGIN = 2, function(x) quantile(x, probs = c(0.025)))
    uci <- apply(trace_dx, MARGIN = 2, function(x) quantile(x, probs = c(0.975)))
    l50 <- apply(trace_dx, MARGIN = 2, function(x) quantile(x, probs = c(0.25)))
    u50 <- apply(trace_dx, MARGIN = 2, function(x) quantile(x, probs = c(0.75)))
    dx <- data.frame(median = med, lci = lci, uci = uci, l50 = l50, u50 = u50)

    trace_cum_dx <- out_sc_trace[["cum_dx"]][1:max_time]
    med <- apply(trace_cum_dx, MARGIN = 2, median)
    lci <- apply(trace_cum_dx, MARGIN = 2, function(x) quantile(x, probs = c(0.025)))
    uci <- apply(trace_cum_dx, MARGIN = 2, function(x) quantile(x, probs = c(0.975)))
    l50 <- apply(trace_cum_dx, MARGIN = 2, function(x) quantile(x, probs = c(0.25)))
    u50 <- apply(trace_cum_dx, MARGIN = 2, function(x) quantile(x, probs = c(0.75)))
    cum_dx <- data.frame(median = med, lci = lci, uci = uci, l50 = l50, u50 = u50)

    trace_r0 <- out_sc_trace[["prd_r0"]][1:max_time]
    r0_med <- apply(trace_r0, MARGIN = 2, median)
    r0_lci <- apply(trace_r0, MARGIN = 2, function(x) quantile(x, probs = c(0.025)))
    r0_uci <- apply(trace_r0, MARGIN = 2, function(x) quantile(x, probs = c(0.975)))
    r0_l50 <- apply(trace_r0, MARGIN = 2, function(x) quantile(x, probs = c(0.25)))
    r0_u50 <- apply(trace_r0, MARGIN = 2, function(x) quantile(x, probs = c(0.75)))

    trace_rx <- out_sc_trace[["prd_rx"]][1:max_time]
    med <- apply(trace_rx, MARGIN = 2, median)
    lci <- apply(trace_rx, MARGIN = 2, function(x) quantile(x, probs = c(0.025)))
    uci <- apply(trace_rx, MARGIN = 2, function(x) quantile(x, probs = c(0.975)))
    l50 <- apply(trace_rx, MARGIN = 2, function(x) quantile(x, probs = c(0.25)))
    u50 <- apply(trace_rx, MARGIN = 2, function(x) quantile(x, probs = c(0.75)))
    rx <- data.frame(median = med, lci = lci, uci = uci, l50 = l50, u50 = u50)

    trace_h <- out_sc_trace[["prd_h"]][1:max_time]
    med <- apply(trace_h, MARGIN = 2, median)
    lci <- apply(trace_h, MARGIN = 2, function(x) quantile(x, probs = c(0.025)))
    uci <- apply(trace_h, MARGIN = 2, function(x) quantile(x, probs = c(0.975)))
    l50 <- apply(trace_h, MARGIN = 2, function(x) quantile(x, probs = c(0.25)))
    u50 <- apply(trace_h, MARGIN = 2, function(x) quantile(x, probs = c(0.75)))
    h <- data.frame(median = med, lci = lci, uci = uci, l50 = l50, u50 = u50)

  out_sc <- list(nc = nc, nh = nh, dx = dx, rx = rx, h = h,
                 cum_nc = cum_nc, cum_nh = cum_nh, cum_dx = cum_dx,
              r0_med = r0_med, r0_lci = r0_lci, r0_uci = r0_uci,
              r0_l50 = r0_l50, r0_u50 = r0_u50)

  return(out_sc)
}


get_out_impact <- function(out_1, out_2, max_time, results = "paf") {

  if (!(results %in% c("paf", "ratio"))){
    break("error - paf or ratio")
  }
  if (results == "paf"){
    trace_1 <- out_1[["cum_nh"]][max_time]
    trace_2 <- out_2[["cum_nh"]][max_time]
    impact <- (trace_1 - trace_2 + 1e-6) / (trace_1 + 1e-6) * 100
    med <- apply(impact, MARGIN = 2, median)
    lci <- apply(impact, MARGIN = 2, function(x) quantile(x, probs = c(0.025)))
    uci <- apply(impact, MARGIN = 2, function(x) quantile(x, probs = c(0.975)))
    l50 <- apply(impact, MARGIN = 2, function(x) quantile(x, probs = c(0.25)))
    u50 <- apply(impact, MARGIN = 2, function(x) quantile(x, probs = c(0.75)))
    out_impact <- data.frame(median = med, lci = lci, uci = uci, l50 = l50, u50 = u50)

    return(round(out_impact))
  }

  if (results == "ratio"){
    trace_1 <- out_1[["cum_nh"]][max_time]
    trace_2 <- out_2[["cum_nh"]][max_time]
    impact <- trace_1 / trace_2
    med <- apply(impact, MARGIN = 2, median)
    lci <- apply(impact, MARGIN = 2, function(x) quantile(x, probs = c(0.025)))
    uci <- apply(impact, MARGIN = 2, function(x) quantile(x, probs = c(0.975)))
    l50 <- apply(impact, MARGIN = 2, function(x) quantile(x, probs = c(0.25)))
    u50 <- apply(impact, MARGIN = 2, function(x) quantile(x, probs = c(0.75)))
    out_impact <- data.frame(median = med, lci = lci, uci = uci, l50 = l50, u50 = u50)
    return(round(out_impact, 1))
  }
}

run_scenario_relache_measures <- function(SEIR_Data, SEIR_Const, fit,
                                          import_scenario = SEIR_Data$num_import,
                                          start_date, date_measures, scenario_date,
                                          rw = FALSE) {

  SEIR_Const_mod <- SEIR_Const
  trace <- as.data.frame(do.call(rbind, fit))
  prd_nc <- matrix(NA, nrow = nrow(trace), ncol = SEIR_Const$Z)
  prd_nh <- matrix(NA, nrow = nrow(trace), ncol = SEIR_Const$Z)
  prd_r0 <- matrix(NA, nrow = nrow(trace), ncol = SEIR_Const$Z)
  prd_rx <- matrix(NA, nrow = nrow(trace), ncol = SEIR_Const$Z)
  prd_h <- matrix(NA, nrow = nrow(trace), ncol = SEIR_Const$Z)
  prd_dx <- matrix(NA, nrow = nrow(trace), ncol = SEIR_Const$Z)

  cum_nc <- cum_nh <- cum_dx <- matrix(NA, nrow = nrow(trace), ncol = SEIR_Const$Z)

  for (i in 1:nrow(trace)) {
    if (rw == FALSE) {
      # we fill in beta_exp
      beta_exp <- rep(exp(trace$b[i]), SEIR_Const$niter)
      for (j in 1:SEIR_Const$n_beta) {
        nom_reduction <- paste("reduction[", j, "]", sep = '')
        beta_exp[SEIR_Const$beta_ind[j]:(SEIR_Const$beta_ind[j + 1] - 1)] <- exp(trace$b[i]) * trace[i, names(trace) %in% nom_reduction]
      }
      nom_reduction <- paste("reduction[", SEIR_Const$n_beta, "]", sep = '')
      beta_exp[SEIR_Const$beta_ind[SEIR_Const$n_beta]:SEIR_Const$niter] <- exp(trace$b[i]) *  trace[i, names(trace) %in% nom_reduction]
    }

    # If we are doing random walk
    if (rw == TRUE) {
      SEIR_Const_mod$beta_ind <- SEIR_Const$beta_ind - as.numeric((date_measures - scenario_date)) / SEIR_Const$dt
      if (max(SEIR_Const_mod$beta_ind) > SEIR_Const$niter) {
        SEIR_Const_mod$beta_ind[length(SEIR_Const_mod)] <- SEIR_Const$niter
      }
      beta_exp <- rep(exp(trace[i, names(trace) %in% "b[1]"]), SEIR_Const$niter)
      for (j in 2:SEIR_Const$n_beta) {
        nom_rw <- paste("b[", j, "]", sep = '')
        beta_exp[SEIR_Const_mod$beta_ind[j - 1]:SEIR_Const_mod$beta_ind[j]] <- exp(trace[i, names(trace) %in% nom_rw])
      }
    }

    num_import_t0 <- import_scenario
    if (is.null(trace$import_t0[i])) {
      num_import_t0[1] <- 1
    } else {
      num_import_t0[1] <- trace$import_t0[i]
    }


    # SEIR model
    if (length(SEIR_Const$rho) == 0) {
      mod <-  seirCpp(beta = beta_exp,
                      infa = SEIR_Const$infa,          # reduction of infectious if asymptomatic
                      p_symp = SEIR_Const$p_symp,       # proportion of cases that are symptomatic
                      eps = SEIR_Const$eps,       # rate from exposed to infectious
                      p_hosp = SEIR_Const$p_hosp,      # probability of hospitalization if symptomatic
                      sigma_h = SEIR_Const$sigma_h,     # rate from symptomatic to hospitalization
                      sigma_c = SEIR_Const$sigma_c,     # rate from infectious to recovered
                      sigma_d = SEIR_Const$sigma_d,       # rate from hospitalized to recovered (alive)
                      conf_rate = SEIR_Data$conf_rate,   # confirmation rate among those with symptoms (delays of 5 days from symptom onset)
                      mu = SEIR_Const$mu, # proportion dying in hospitals of covid-mortality
                      num_import = num_import_t0,
                      prp_detect = SEIR_Const$prp_detect,
                      p_conf = SEIR_Data$p_conf,     # scaling fator for case confirmation.
                      init = SEIR_Data$init,
                      start = SEIR_Const$start, end = SEIR_Const$end, dt = SEIR_Const$dt)

      duration_inf <- SEIR_Const[["sigma_c"]]
    }

    # SEAIR model
    if (length(SEIR_Const$rho) == 1) {
      mod <-  seairCpp(beta = beta_exp,
                       infa = SEIR_Const$infa,          # reduction of infectious if asymptomatic
                       p_symp = SEIR_Const$p_symp,       # proportion of cases that are symptomatic
                       rho = SEIR_Const$rho,
                       eps = SEIR_Const$eps,       # rate from exposed to infectious
                       p_hosp = SEIR_Const$p_hosp,      # probability of hospitalization if symptomatic
                       sigma_h = SEIR_Const$sigma_h,     # rate from symptomatic to hospitalization
                       sigma_c = SEIR_Const$sigma_c,     # rate from infectious to recovered
                       sigma_d = SEIR_Const$sigma_d,       # rate from hospitalized to recovered (alive)
                       conf_rate = SEIR_Data$conf_rate,   # confirmation rate among those with symptoms (delays of 5 days from symptom onset)
                       mu = SEIR_Const$mu, # proportion dying in hospitals of covid-mortality
                       num_import = num_import_t0,
                       prp_detect = SEIR_Const$prp_detect,
                       p_conf = SEIR_Data$p_conf,     # scaling fator for case confirmation.
                       init = SEIR_Data$init,
                       start = SEIR_Const$start, end = SEIR_Const$end, dt = SEIR_Const$dt)

      duration_inf <- 1 / (1 / SEIR_Const[["eps"]] + 1 / SEIR_Const[["sigma_c"]])
    }

    # putting everything together
    prop_s <- 1 - (mod[, "R"][SEIR_Const$iter_time] / SEIR_Data[["init"]][1])

    r_nc <- trace$theta_nc[i]
    p_nc <- r_nc / (r_nc + mod[, "NC"][SEIR_Const$iter_time] + 1e-6)
    prd_nc[i, ] <- rnbinom(length(p_nc), prob = p_nc, size = r_nc)
    cum_nc[i, ] <- cumsum(prd_nc[i, ])

    r_nh <- trace$theta_nh[i]
    p_nh <- r_nh / (r_nh + mod[, "NH"][SEIR_Const$iter_time] + 1e-6)
    prd_nh[i, ] <- rnbinom(length(p_nh), prob = p_nh, size = r_nh)
    cum_nh[i, ] <- cumsum(prd_nh[i, ])

    r_dx <- trace$theta_dx[i]
    p_dx <- r_dx / (r_dx + mod[, "DX"][SEIR_Const$iter_time] + 1e-6)
    prd_dx[i, ] <- rnbinom(length(p_dx), prob = p_dx, size = r_dx)
    cum_dx[i, ] <- cumsum(prd_dx[i, ])

    prd_r0[i, ] <- (SEIR_Const[["p_symp"]] * mod[, "beta_iter"][SEIR_Const$iter_time] / duration_inf +
                      (1 - SEIR_Const[["p_symp"]]) * mod[, "beta_iter"][SEIR_Const$iter_time] * SEIR_Const[["infa"]] / duration_inf) * prop_s
    prd_rx[i, ] <- mod[, "R"][SEIR_Const$iter_time]
    prd_h[i, ] <- mod[, "H"][SEIR_Const$iter_time]
  }

  out_sc_trace <- list(prd_nc = as.data.frame(prd_nc),
                       prd_nh = as.data.frame(prd_nh),
                       prd_dx = as.data.frame(prd_dx),
                       cum_nc = as.data.frame(cum_nc),
                       cum_nh = as.data.frame(cum_nh),
                       cum_dx = as.data.frame(cum_dx),
                       prd_r0 = as.data.frame(prd_r0),
                       prd_rx = as.data.frame(prd_rx),
                       prd_h  = as.data.frame(prd_h))
  return(out_sc_trace)
}



