combine_par_mcmc <- function(fit_par) {
  # for example, true if 4 parallel runs with 2 chains each, false if 8 runs with 1 chain each.
  has_mult_chain_per_run <- length(fit_par[[1]]) < length(fit_par[[1]][[1]])
  n_par <- length(fit_par)

  if (has_mult_chain_per_run) {
    n_chain <- length(fit_par[[1]])
    index <- 0
    out_lst <- list()
    for (i in 1:n_par) {
      for (j in 1:n_chain) {
        index <- index + 1
        trace_fit <- as.data.frame(do.call(rbind, fit_par[[i]][j]))
        out_lst[[index]] <- coda::mcmc(trace_fit)
      }
    }
    return(coda::mcmc.list(out_lst))
  }

  # Single chain per parallel run.
  return(coda::mcmc.list(lapply(lapply(fit_par, as.data.frame), coda::mcmc)))
}




parallel_nimble <- function(SEIR_Data, SEIR_Const,
                            n_iterations,
                            n_chains,
                            warm_up,
                            n_thin) {
require(nimble)
FUN_SEIRc <- nimble::nimbleRcall(function(
        beta =        double(1),
        infa =        double(0),        # reduction of infectious if asymptomatic
        p_symp =      double(0),      # proportion of cases that are symptomatic
        eps =         double(0),         # rate from exposed to infectious
        p_hosp =      double(0),      # probabiliyt of hospitalization if symptomatic
        sigma_h =     double(0),     # rate from symptomatic to hospitalization
        sigma_c =     double(0),     # rate from infectious to recovered
        sigma_d =     double(0),     # rate from hospitalized to recovered (alive)
        conf_rate =   double(1),   # confirmation rate among those with symptoms (delays of 5 days from symptom onset)
        mu =          double(0),          # covid-mortality rate among hospitalized patients
        num_import =  double(1),
        prp_detect =  double(0),
        p_conf     =  double(1),
        init       =  double(1),
        start = double(0), end = double(0), dt = double(0)) {},
        Rfun = "seirCpp", returnType = double(2))

# R0 post lock-down: 0.62 (0.37-0.89)
# https://cmmid.github.io/topics/covid19/current-patterns-transmission/comix-impact-of-physical-distance-measures-on-transmission-in-the-UK.html
# R0 post lock-down: 1.4 accross all countries (a 64% reduction from pre-lockdown level)
# https://www.imperial.ac.uk/media/imperial-college/medicine/mrc-gida/2020-03-30-COVID19-Report-13.pdf
# Wuhan lockdown and traffic ban R0 = 1.25 (required centralized quarantine for R to be 0.32)
# https://www.medrxiv.org/content/10.1101/2020.03.03.20030593v1.full.pdf
# https://cdn1.sph.harvard.edu/wp-content/uploads/sites/21/2020/03/COVID-19-03-16-2020-Lin.pdf


# Write your nimble model here - the likelihood + priors
SEIR_Code <- nimble::nimbleCode({
  # prior
  # exp(log(2.2 * sigma_c) - 1.96 * 0.3) / sigma_c; exp(log(2.2 * sigma_c) + 1.96 * 0.3) / sigma_c
  # R0 in symptomatic is: 2.2 / (-infa * p_symp + infa + p_symp)
  b ~ dnorm(log(2.44 * sigma_c), sd = 0.3)
  b_exp[1:(intervention1_ind - 1)] <- exp(b)
  for (i in 1:n_beta) {
    reduction[i] ~ dunif(0, 1)
    b_exp[beta_ind[i]:(beta_ind[i + 1] - 1)] <- exp(b) * reduction[i]
    }
  b_exp[beta_ind[n_beta]:niter] <- exp(b) * reduction[n_beta]

  import_t0 ~ dunif(0, 30)
  num_import[1] <- import_t0

  # we get the outcomes
  mod_y[1:niter, 1:ncol_val] <- FUN_SEIRc(b_exp[1:niter],
                          infa, p_symp, eps,
                          p_hosp, sigma_h, sigma_c, sigma_d,
                          conf_rate[1:niter], mu,
                          num_import[1:niter], prp_detect, p_conf[1:niter],
                          init[1:14],
                          start, end, dt)

  for (i in 1:N) {
    # negative binomial likelihood
    nc[i] ~ dnegbin(prob = p_nc[i], size = r_nc)
    nh[i] ~ dnegbin(prob = p_nh[i], size = r_nh)
    dx[i] ~ dnegbin(prob = p_dx[i], size = r_dx)
    p_nc[i] <- r_nc / (r_nc + (mod_y[y_time[i], ind_nc]) + 1e-6)
    p_nh[i] <- r_nh / (r_nh + (mod_y[y_time[i], ind_nh]) + 1e-6)
    p_dx[i] <- r_dx / (r_dx + (mod_y[y_time[i], ind_dx]) + 1e-6)
    }
    r_nc <- theta_nc
    r_nh <- theta_nh
    r_dx <- theta_dx

  # dispersion for negative binomial
  theta_nc ~ dgamma(shape = 10, rate = 10/1000) # dgamma(b, b/a) where a is the mean and var is a^(2/b)
  theta_nh ~ dgamma(shape = 10, rate = 10/1000)
  theta_dx ~ dgamma(shape = 10, rate = 10/1000)

  for (z in 1:Z) {
    prd_nc[z] ~ dnegbin(prd_p_nc[z], prd_r_nc)
    prd_p_nc[z] <- prd_r_nc / (prd_r_nc + mod_y[iter_time[z], ind_nc] + 1e-6)

    prd_nh[z] ~ dnegbin(prd_p_nh[z], prd_r_nh)
    prd_p_nh[z] <- prd_r_nh / (prd_r_nh + mod_y[iter_time[z], ind_nh] + 1e-6)

    prd_dx[z] ~ dnegbin(prd_p_dx[z], prd_r_dx)
    prd_p_dx[z] <- prd_r_dx / (prd_r_dx + mod_y[iter_time[z], ind_dx] + 1e-6)

    prd_rx[z] <- mod_y[iter_time[z], ind_rx]
    beta_iter[z] <- mod_y[iter_time[z], ind_b]
  }
    prd_r_nc <- theta_nc
    prd_r_nh <- theta_nh
    prd_r_dx <- theta_dx
})

SEIR_Inits <- list(b = log(runif(1, min = 1.5, max = 3) / (1 / SEIR_Const[["sigma_c"]]) ) )

SEIR_Model <- nimble::nimbleModel(SEIR_Code, SEIR_Const, SEIR_Data, inits = SEIR_Inits)


SEIR_Model <- nimble::nimbleModel(SEIR_Code, SEIR_Const, SEIR_Data)

# we compile the model
cSEIR <- nimble::compileNimble(SEIR_Model)

SEIRConf <- nimble::configureMCMC(cSEIR, thin = n_thin, autoblock = TRUE)
SEIRConf$addMonitors(c("b", "import_t0",
                       "reduction",
                       "theta_nc", "theta_nh", "theta_dx",
                       "prd_nc", "prd_nh", "prd_dx", "prd_rx", "beta_iter"))
#SEIRConf$removeSamplers(c("b", "reduction", "import_t0"))
#SEIRConf$addSampler(target = c("b", "reduction[1]", "import_t0"), type = 'AF_slice')
#for (i in 2:SEIR_Const$n_beta) {
#  nom_red <- paste("reduction[", i, "]", sep = "")
#  SEIRConf$addSampler(target = nom_red, type = 'slice')
#}
SEIRConf$removeSamplers(c("b", "import_t0", "reduction"))
SEIRConf$addSampler(target = c("b", "import_t0", "reduction"), type = 'AF_slice')

SEIRMCMC <- nimble::buildMCMC(SEIRConf)
CSEIRMCMC <- nimble::compileNimble(SEIRMCMC, project = SEIR_Model, resetFunctions = TRUE)

fit <- nimble::runMCMC(CSEIRMCMC, niter = n_iterations, nchains = n_chains,
                   nburnin = warm_up, samplesAsCodaMCMC = TRUE)

return(fit)

}





parallel_rw <- function(SEIR_Data, SEIR_Const,
                            n_iterations,
                            n_chains,
                            warm_up,
                            n_thin) {
require(nimble)
FUN_SEIRc <- nimble::nimbleRcall(function(
        beta =        double(1),
        infa =        double(0),        # reduction of infectious if asymptomatic
        p_symp =      double(0),      # proportion of cases that are symptomatic
        eps =         double(0),         # rate from exposed to infectious
        p_hosp =      double(0),      # probabiliyt of hospitalization if symptomatic
        sigma_h =     double(0),     # rate from symptomatic to hospitalization
        sigma_c =     double(0),     # rate from infectious to recovered
        sigma_d =     double(0),     # rate from hospitalized to recovered (alive)
        conf_rate =   double(1),   # confirmation rate among those with symptoms (delays of 5 days from symptom onset)
        mu =          double(0),          # covid-mortality rate among hospitalized patients
        num_import =  double(1),
        prp_detect =  double(0),
        p_conf     =  double(1),
        init       =  double(1),
        start = double(0), end = double(0), dt = double(0)) {},
        Rfun = "seirCpp", returnType = double(2))

# R0 post lock-down: 0.62 (0.37-0.89)
# https://cmmid.github.io/topics/covid19/current-patterns-transmission/comix-impact-of-physical-distance-measures-on-transmission-in-the-UK.html
# R0 post lock-down: 1.4 accross all countries (a 64% reduction from pre-lockdown level)
# https://www.imperial.ac.uk/media/imperial-college/medicine/mrc-gida/2020-03-30-COVID19-Report-13.pdf
# Wuhan lockdown and traffic ban R0 = 1.25 (required centralized quarantine for R to be 0.32)
# https://www.medrxiv.org/content/10.1101/2020.03.03.20030593v1.full.pdf
# https://cdn1.sph.harvard.edu/wp-content/uploads/sites/21/2020/03/COVID-19-03-16-2020-Lin.pdf


SEIR_Data$ones <- rep(1, SEIR_Const$n_beta)

# Write your nimble model here - the likelihood + priors
SEIR_Code_RW <- nimble::nimbleCode({
  # prior
  # exp(log(2.2 * sigma_c) - 1.96 * 0.3) / sigma_c; exp(log(2.2 * sigma_c) + 1.96 * 0.3) / sigma_c
  # R0 in symptomatic is: 2.2 / (-infa * p_symp + infa + p_symp)
  b[1] ~ dnorm(log(2.6 * sigma_c), sd = 0.2)
  #b[1] ~ T(dnorm(mean = log(2.44 * sigma_c), sd = 0.1), log(1 * sigma_c), log(8 * sigma_c))
  b_exp[1:(beta_ind[1] - 1)] <- exp(b[1])
  for (i in 2:n_beta) {
    #b[i] ~ T(dnorm(b[i - 1], sd_rw), log(0.1 * sigma_c), log(8 * sigma_c))
    b[i] ~ dnorm(mean = b[i - 1], sd = sd_rw)
    b_exp[beta_ind[i - 1]:(beta_ind[i] - 1)] <- exp(b[i])
  }
  sd_rw ~ T(dt(mu = 0, sigma = 25, df = 1), 0, 100) # half-cauchy prior
  #sd_rw ~ dgamma(shape = 1, rate = 1 / 1) # dgamma(shape = 5, rate = 1 / 5) # dgamma(b, b/a) where a is the mean and var is a^(2/b)
  b_exp[beta_ind[n_beta - 1]:niter] <- exp(b[n_beta])

  b_last_log ~ dnorm(mean = b[n_beta], sd = sd_rw)
  b_last <- exp(b_last_log)

  import_t0 ~ dunif(1, 30)
  num_import[1] <- import_t0
  # we get the outcomes
  mod_y[1:niter, 1:ncol_val] <- FUN_SEIRc(b_exp[1:niter],
                          infa, p_symp, eps,
                          p_hosp, sigma_h, sigma_c, sigma_d,
                          conf_rate[1:niter], mu,
                          num_import[1:niter], prp_detect, p_conf[1:niter],
                          init[1:14],
                          start, end, dt)
  for (i in 1:N) {
    # negative binomial likelihood
    nc[i] ~ dnegbin(prob = p_nc[i], size = r_nc)
    nh[i] ~ dnegbin(prob = p_nh[i], size = r_nh)
    dx[i] ~ dnegbin(prob = p_dx[i], size = r_dx)
    p_nc[i] <- r_nc / (r_nc + (mod_y[y_time[i], ind_nc]) + 1e-6)
    p_nh[i] <- r_nh / (r_nh + (mod_y[y_time[i], ind_nh]) + 1e-6)
    p_dx[i] <- r_dx / (r_dx + (mod_y[y_time[i], ind_dx]) + 1e-6)
    }
    r_nc <- theta_nc
    r_nh <- theta_nh
    r_dx <- theta_dx
  # dispersion for negative binomial
  theta_nc ~ dgamma(shape = 10, rate = 10 / 1000) # dgamma(b, b/a) where a is the mean and var is a^(2/b)
  theta_nh ~ dgamma(shape = 10, rate = 10 / 1000)
  theta_dx ~ dgamma(shape = 10, rate = 10 / 1000)
  for (z in 1:Z) {
    prd_nc[z] ~ dnegbin(prd_p_nc[z], r_nc)
    prd_p_nc[z] <- r_nc / (r_nc + mod_y[iter_time[z], ind_nc] + 1e-6)
    prd_nh[z] ~ dnegbin(prd_p_nh[z], r_nh)
    prd_p_nh[z] <- r_nh / (r_nh + mod_y[iter_time[z], ind_nh] + 1e-6)
    prd_dx[z] ~ dnegbin(prd_p_dx[z], r_dx)
    prd_p_dx[z] <- r_dx / (r_dx + mod_y[iter_time[z], ind_dx] + 1e-6)
    prd_rx[z] <- mod_y[iter_time[z], ind_rx]
    beta_iter[z] <- mod_y[iter_time[z], ind_b]
  }
})

# 0.7 / (-SEIR_Const$infa * SEIR_Const$p_symp + SEIR_Const$infa + SEIR_Const$p_symp)
b_inits <- log(c(runif(1, min = 2, max = 3),
             runif(1, min = 1, max = 2),
             runif(SEIR_Const$n_beta - 2, min = 0.5, max = 1)) / (1 / SEIR_Const[["sigma_c"]]) )
SEIR_Inits <-  list(b = b_inits,
                  #import_t0 = 15,
                   theta_nc = 1000,
                   theta_nh = 1000,
                   theta_dx = 1000,
                   sd_rw = runif(1, min = 0.1, max = 1))

SEIR_Model_RW <- nimble::nimbleModel(SEIR_Code_RW, SEIR_Const, SEIR_Data, inits = SEIR_Inits)

# we compile the model
cSEIR_RW <- nimble::compileNimble(SEIR_Model_RW)

SEIRConf_RW <- nimble::configureMCMC(cSEIR_RW, thin = n_thin, autoblock = TRUE)
SEIRConf_RW$addMonitors(c("b", "import_t0",
                       "sd_rw", "b_last",
                       "theta_nc", "theta_nh", "theta_dx",
                       "prd_nc", "prd_nh", "prd_dx", "prd_rx", "beta_iter"))
SEIRConf_RW$removeSamplers(c("b", "import_t0", "sd_rw"))
SEIRConf_RW$addSampler(target = c("b", "import_t0", "sd_rw"), type = 'AF_slice')

SEIRMCMC_RW <- nimble::buildMCMC(SEIRConf_RW)
CSEIRMCMC_RW <- nimble::compileNimble(SEIRMCMC_RW, project = SEIR_Model_RW, resetFunctions = TRUE)

fit <- nimble::runMCMC(CSEIRMCMC_RW, niter = n_iterations, nchains = n_chains,
                   nburnin = warm_up, samplesAsCodaMCMC = TRUE)

return(fit)

}

#' Adjust parallelization for high perf. computing environments
get_parallel_params <- function() {
  is_high_core_count <- parallel::detectCores() >= 32
  if (is_high_core_count) print("Detected High Perf. Environment")
  return(list(
    "n_runs" = ifelse(is_high_core_count, 8, 4),
    "n_chains" = ifelse(is_high_core_count, 1, 2)
  ))
}
