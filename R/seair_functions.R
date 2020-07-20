seair <- function(param, init, start = 0, end = 250, dt) {

  val <- seairCpp(beta = exp(param[["beta"]]),
         infa =       param[["infa"]],        # reduction of infectious if asymptomatic
         p_symp =     param[["p_symp"]],      # proportion of cases that are symptomatic
         rho =        param[["rho"]],
         eps =        param[["eps"]],         # rate from exposed to infectious
         p_hosp =     param[["p_hosp"]],      # probabiliyt of hospitalization if symptomatic
         sigma_h =    param[["sigma_h"]],     # rate from symptomatic to hospitalization
         sigma_c =    param[["sigma_c"]],     # rate from infectious to recovered
         sigma_d =    param[["sigma_d"]],     # rate from hospitalized to recovered (alive)
         conf_rate =  param[["conf_rate"]],   # confirmation rate among those with symptoms (delays of 5 days from symptom onset)
         mu =         param[["mu"]],          # covid-mortality rate among hospitalized patients
         num_import = param[["num_import"]],
         prp_detect = param[["prp_detect"]], # proportion of imported cases detected
         p_conf     = param[["p_conf"]],
         init       = init,
         start = start, end = end, dt = dt)
  val <- as.data.frame(val)
  return(val)
}


Rcpp::cppFunction('NumericMatrix seairCpp(
    NumericVector beta,
    double infa,
    double p_symp,
    double rho,
    double eps,
    double p_hosp,
    double sigma_h,
    double sigma_c,
    double sigma_d,
    NumericVector conf_rate,
    double mu,
    IntegerVector num_import,
    double prp_detect,
    NumericVector p_conf,
    NumericVector init,
    int start, int end, double dt) {

  double sigma_eps = 1 / (1 / sigma_c + 1 / eps);
  double S_init = init[0];
  double E1_init = init[1];
  double E2_init = init[2];
  double P1_init = init[3];
  double P2_init = init[4];
  double I1a_init = init[5];
  double I2a_init = init[6];
  double I1s_init = init[7];
  double I2s_init = init[8];
  double I1h_init = init[9];
  double I2h_init = init[10];
  double H_init = init[11];
  double R_init = init[12];
  double D_init = init[13];
  double O_init = init[14];
  double C_init = init[15];
  double zero = 0;

  double N_init = S_init +
                E1_init + E2_init +
                P1_init + P2_init +
                I1s_init + I2s_init +
                I1a_init + I2a_init +
                I1h_init + I2h_init +
                H_init + R_init;

  int niter = (end - start) / dt + 1;
  NumericVector time(niter);
  time[0] = start;

  // we index beta over the nubmer of iterations
  NumericVector beta_iter = beta;

  NumericVector S(niter, S_init);
  NumericVector E1(niter, E1_init);
  NumericVector E2(niter, E2_init);
  NumericVector P1(niter, P1_init);
  NumericVector P2(niter, P2_init);
  NumericVector I1s(niter, I1s_init);
  NumericVector I2s(niter, I2s_init);
  NumericVector I1a(niter, I1a_init);
  NumericVector I2a(niter, I2a_init);
  NumericVector I1h(niter, I1h_init);
  NumericVector I2h(niter, I2h_init);
  NumericVector H(niter, H_init);
  NumericVector R(niter, R_init);
  NumericVector D(niter, D_init);
  NumericVector N(niter, N_init);

  NumericVector O(niter, O_init);
  NumericVector C(niter, C_init);
  NumericVector NC(niter, zero);
  NumericVector DX(niter, zero);
  NumericVector NH(niter, zero);

  //int StoE;
  for(int i = 0; i < (niter - 1); ++i) {
    double S_i = S[i];
    double E1_i = E1[i];
    double E2_i = E2[i];
    double P1_i = P1[i];
    double P2_i = P2[i];
    double I1s_i = I1s[i];
    double I2s_i = I2s[i];
    double I1a_i = I1a[i];
    double I2a_i = I2a[i];
    double I1h_i = I1h[i];
    double I2h_i = I2h[i];
    double H_i = H[i];
    double R_i = R[i];
    double D_i = D[i];
    double N_i = N[i];

    double O_i = O[i];
    double C_i = C[i];

    // Susceptible
    S[i + 1] = S_i + dt * (- beta_iter[i] * S_i *
              (P1_i + P2_i + I1s_i + I2s_i + ((I1a_i + I2a_i) * infa) + I1h_i + I2h_i) / N_i);

    // Erlang distributed exposed period
    E1[i + 1] = E1_i + dt * (+ beta_iter[i] * S_i *
                (P1_i + P2_i + I1s_i + I2s_i + ((I1a_i + I2a_i) * infa) + I1h_i + I2h_i) / N_i - (rho * 2) * E1_i);
    E2[i + 1] = E2_i + dt * (+ (rho * 2) * E1_i - (rho * 2) * E2_i);

    // Pre-clinical (erlang distributed period)
    P1[i + 1] = P1_i + dt * (+ p_symp * (rho * 2) * E2_i - (eps * 2) * P1_i) + num_import[i];
    P2[i + 1] = P2_i + dt * (+ (eps * 2) * P1_i          - (eps * 2) * P2_i);

    // Erlang distributed infectious period
    // Infectious symptomatic (mild symptoms - no hospitalization)
    I1s[i + 1] = I1s_i + dt * (+ (1 - p_hosp) * (eps * 2) * P2_i - (sigma_c * 2) * I1s_i);
    I2s[i + 1] = I2s_i + dt * (+ (sigma_c * 2) * I1s_i           - (sigma_c * 2) * I2s_i);

    // Infectious severe symptomatic (they will be hospitalized)
    I1h[i + 1] = I1h_i + dt * (+ p_hosp * (eps * 2) * P2_i - (sigma_h * 2) * I1h_i);
    I2h[i + 1] = I2h_i + dt * (+ (sigma_h * 2) * I1h_i     - (sigma_h * 2) * I2h_i);

    // Infectious asymptomatic
    I1a[i + 1] = I1a_i + dt * (+ (1 - p_symp) * (rho * 2) * E2_i - (sigma_eps * 2) * I1a_i) + num_import[i] * (1 / p_symp);
    I2a[i + 1] = I2a_i + dt * (+ (sigma_eps * 2) * I1a_i         - (sigma_eps * 2) * I2a_i);

    // Hospitalized
    H[i + 1]  = H_i  + dt * ((sigma_h * 2) * I2h_i - sigma_d * H_i);

    // Recovered
    R[i + 1]  = R_i  + dt * ((sigma_eps * 2) * I2a_i +
                             (sigma_c * 2) * I2s_i +
                              sigma_d * H_i);

    // Dead
    D[i + 1]  = D_i  + dt * (sigma_d * mu * H_i);

    // Total population
    N[i + 1]  = S[i + 1] +
                E1[i + 1] + E2[i + 1] +
                P1[i + 1] + P2[i + 1] +
                I1s[i + 1] + I2s[i + 1] +
                I1a[i + 1] + I2a[i + 1] +
                I1h[i + 1] + I2h[i + 1] +
                H[i + 1] + R[i + 1];

    // Onset of symptoms
    O[i + 1]  = O_i + dt * ((eps * 2) * P2_i - conf_rate[i] * O_i) + num_import[i] * prp_detect;

    // Cumulative confirmed cases (after a reporting delay)
    C[i + 1]  = C_i + dt * (p_conf[i] * conf_rate[i] * O_i);

    // Daily confirmed cases
    NC[i + 1] = p_conf[i] * conf_rate[i] * O_i;

    // Daily deaths
    DX[i + 1] = sigma_d * mu * H_i;

    // New daily hospitalization
    NH[i + 1] = (sigma_h * 2) * I2h_i;

    time[i + 1] = time[i] + dt;
  }

  NumericVector E = E1 + E2;
  NumericVector P = P1 + P2;
  NumericVector Is = I1s + I2s;
  NumericVector Ia = I1a + I2a;
  NumericVector Ih = I1h + I2h;

    DataFrame seir_dat = DataFrame::create(Named("S") = S,
            Named("E")  = E, Named("P") = P,
            Named("Is") = Is, Named("Ia") = Ia, Named("Ih") = Ih,
            Named("R")  = R, Named("H") = H,
            Named("D")  = D,
            Named("O")  = O,  Named("C")  = C,
            Named("NC") = NC, Named("DX") = DX,
            Named("N")  = N,  Named("NH") = NH,
            Named("beta_iter") = beta_iter,
            Named("time") = time);

  int seir_size = seir_dat.size();
  NumericMatrix seir(seir_dat.nrows(), seir_size);
  for(int i = 0; i < seir_size; i++) {
    seir(_, i) = NumericVector(seir_dat[i]);
  }

  colnames(seir) = CharacterVector::create(
      "S", "E", "P",
      "Is", "Ia", "Ih",
      "R", "H",
      "D",
      "O", "C",
      "NC", "DX",
      "N", "NH", "beta_iter", "time" );
  return(seir);
}')


parallel_seair <- function(SEIR_Data, SEIR_Const,
                            n_iterations,
                            n_chains,
                            warm_up,
                            n_thin) {
require(nimble)
FUN_SEAIRc <- nimble::nimbleRcall(function(
        beta =        double(1),
        infa =        double(0),        # reduction of infectious if asymptomatic
        p_symp =      double(0),      # proportion of cases that are symptomatic
        rho =         double(0),
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
        Rfun = "seairCpp", returnType = double(2))

# R0 post lock-down: 0.62 (0.37-0.89)
# https://cmmid.github.io/topics/covid19/current-patterns-transmission/comix-impact-of-physical-distance-measures-on-transmission-in-the-UK.html
# R0 post lock-down: 1.4 accross all countries (a 64% reduction from pre-lockdown level)
# https://www.imperial.ac.uk/media/imperial-college/medicine/mrc-gida/2020-03-30-COVID19-Report-13.pdf
# Wuhan lockdown and traffic ban R0 = 1.25 (required centralized quarantine for R to be 0.32)
# https://www.medrxiv.org/content/10.1101/2020.03.03.20030593v1.full.pdf
# https://cdn1.sph.harvard.edu/wp-content/uploads/sites/21/2020/03/COVID-19-03-16-2020-Lin.pdf


# Write your nimble model here - the likelihood + priors
SEAIR_Code <- nimble::nimbleCode({
  # prior
  # exp(log(2.6 * sigma_c) - 1.96 * 0.2) / sigma_c; exp(log(2.6 * sigma_c) + 1.96 * 0.2) / sigma_c
  # R0 in symptomatic is: 2.6 / (-infa * p_symp + infa + p_symp)
  b ~ dnorm(log(2.6 * sigma_c), sd = 0.2)
  b_exp[1:(intervention1_ind - 1)] <- exp(b)
  for (i in 1:n_beta) {
    reduction[i] ~ dunif(0, 1)
    b_exp[beta_ind[i]:(beta_ind[i + 1] - 1)] <- exp(b) * reduction[i]
    }
  b_exp[beta_ind[n_beta]:niter] <- exp(b) * reduction[n_beta]

  import_t0 ~ dunif(0, max_t0)
  num_import[1] <- import_t0

  # we get the outcomes
  mod_y[1:niter, 1:ncol_val] <- FUN_SEAIRc(b_exp[1:niter],
                          infa, p_symp, rho, eps,
                          p_hosp, sigma_h, sigma_c, sigma_d,
                          conf_rate[1:niter], mu,
                          num_import[1:niter], prp_detect, p_conf[1:niter],
                          init[1:16],
                          start, end, dt)

  for (i in 1:N) {
    # negative binomial likelihood
    nc[i] ~ dnegbin(prob = p_nc[i], size = r_nc)
    nh[i] ~ dnegbin(prob = p_nh[i], size = r_nh)
    dx[i] ~ dnegbin(prob = p_dx[i], size = r_dx)
    p_nc[i] <- r_nc / (r_nc + (mod_y[y_time[i], ind_nc]) + 1e-6)
    p_nh[i] <- r_nh / (r_nh + (mod_y[y_time[i], ind_nh]) + 1e-6)
    p_dx[i] <- r_dx / (r_dx + (mod_y[y_time[i], ind_dx]) + 1e-6)
    r_nc <- theta_nc
    r_nh <- theta_nh
    r_dx <- theta_dx
    }

  # dispersion for negative binomial
  theta_nc ~ dgamma(shape = 10, rate = 10/1000) # dgamma(b, b/a) where a is the mean and var is a^(2/b)
  theta_nh ~ dgamma(shape = 10, rate = 10/1000)
  theta_dx ~ dgamma(shape = 10, rate = 10/1000)

  for (z in 1:Z) {
    prd_nc[z] ~ dnegbin(prd_p_nc[z], prd_r_nc)
    prd_p_nc[z] <- prd_r_nc / (prd_r_nc + mod_y[iter_time[z], ind_nc] + 1e-6)
    prd_r_nc <- theta_nc

    prd_nh[z] ~ dnegbin(prd_p_nh[z], prd_r_nh)
    prd_p_nh[z] <- prd_r_nh / (prd_r_nh + mod_y[iter_time[z], ind_nh] + 1e-6)
    prd_r_nh <- theta_nh

    prd_dx[z] ~ dnegbin(prd_p_dx[z], prd_r_dx)
    prd_p_dx[z] <- prd_r_dx / (prd_r_dx + mod_y[iter_time[z], ind_dx] + 1e-6)
    prd_r_dx <- theta_dx

    prd_rx[z] <- mod_y[iter_time[z], ind_rx]
    beta_iter[z] <- mod_y[iter_time[z], ind_b]
  }
})

SEAIR_Inits <- list(b = log(runif(1, min = 1.5, max = 2.75) / (1 / SEIR_Const[["eps"]] + 1 / SEIR_Const[["sigma_c"]]) ),
                   import_t0 = runif(1, min = 0, max = SEIR_Const$max_t0),
                   theta_nc = runif(1, min = 10, max = 1000),
                   theta_nh = runif(1, min = 10, max = 1000),
                   theta_dx = runif(1, min = 10, max = 1000))

SEAIR_Model <- nimble::nimbleModel(SEAIR_Code, SEIR_Const, SEIR_Data, inits = SEAIR_Inits)

# we compile the model
cSEAIR <- nimble::compileNimble(SEAIR_Model)

SEAIRConf <- nimble::configureMCMC(cSEAIR, thin = n_thin, autoblock = TRUE)
SEAIRConf$addMonitors(c("b", "import_t0",
                       "reduction",
                       "theta_nc", "theta_nh", "theta_dx",
                       "prd_nc", "prd_nh", "prd_dx", "prd_rx", "beta_iter"))
#SEIRConf$removeSamplers(c("b", "reduction", "import_t0"))
#SEIRConf$addSampler(target = c("b", "reduction[1]", "import_t0"), type = 'AF_slice')
#for (i in 2:SEIR_Const$n_beta) {
#  nom_red <- paste("reduction[", i, "]", sep = "")
#  SEIRConf$addSampler(target = nom_red, type = 'slice')
#}
SEAIRConf$removeSamplers(c("b", "reduction"))
SEAIRConf$addSampler(target = c("b", "reduction"), type = 'AF_slice')

SEAIRMCMC <- nimble::buildMCMC(SEAIRConf)
CSEAIRMCMC <- nimble::compileNimble(SEAIRMCMC, project = SEAIR_Model, resetFunctions = TRUE)

fit <- nimble::runMCMC(CSEAIRMCMC, niter = n_iterations, nchains = n_chains,
                   nburnin = warm_up, samplesAsCodaMCMC = TRUE)

return(fit)

}

parallel_seair_rw <- function(SEIR_Data, SEIR_Const,
                            n_iterations,
                            n_chains,
                            warm_up,
                            n_thin) {

if (is.null( SEIR_Const[["prior_sd_rw"]])) {
  SEIR_Const[["prior_sd_rw"]] <- c(0, 25)
}
require(nimble)
FUN_SEAIRc <- nimble::nimbleRcall(function(
        beta =        double(1),
        infa =        double(0),        # reduction of infectious if asymptomatic
        p_symp =      double(0),      # proportion of cases that are symptomatic
        rho =         double(0),
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
        Rfun = "seairCpp", returnType = double(2))


# Write your nimble model here - the likelihood + priors
SEAIR_Code <- nimble::nimbleCode({
  # prior
  # exp(log(2.6 * SEIR_Const$sigma_c) - 1.96 * 0.2) / SEIR_Const$sigma_c; exp(log(2.6 * SEIR_Const$sigma_c) + 1.96 * 0.2) / sigma_c
  # R0 in symptomatic is: 2.6 / (-SEIR_Const$infa * SEIR_Const$p_symp + SEIR_Const$infa + SEIR_Const$p_symp)
  b[1] ~ dnorm(log(2.6 * sigma_c), sd = 0.2)
  b_exp[1:(beta_ind[1] - 1)] <- exp(b[1])
  for (i in 2:n_beta) {
    b[i] ~ dnorm(mean = b[i - 1], sd = sd_rw)
    b_exp[beta_ind[i - 1]:(beta_ind[i] - 1)] <- exp(b[i])
  }
  sd_rw ~ T(dt(mu = prior_sd_rw[1], sigma = prior_sd_rw[2], df = 1), 0, 100) # half-cauchy prior
  #sd_rw ~ dgamma(shape = 1, rate = 1 / 1) # dgamma(shape = 5, rate = 1 / 5) # dgamma(b, b/a) where a/b is the mean and var is a/b^2
  # plot(dgamma(seq(0, 10, 0.01), shape = 0.1, scale = 0.1) ~ seq(0, 10, 0.01), type = "l")
  # hist(rgamma(10000, shape = 100, scale = 0.1), breaks = 100)
  b_exp[beta_ind[n_beta - 1]:niter] <- exp(b[n_beta])

  b_last_log ~ dnorm(mean = b[n_beta], sd = sd_rw)
  b_last <- exp(b_last_log)

  import_t0 ~ dunif(0, max_t0)
  num_import[1] <- import_t0

  # we get the outcomes
  mod_y[1:niter, 1:ncol_val] <- FUN_SEAIRc(b_exp[1:niter],
                          infa, p_symp, rho, eps,
                          p_hosp, sigma_h, sigma_c, sigma_d,
                          conf_rate[1:niter], mu,
                          num_import[1:niter], prp_detect, p_conf[1:niter],
                          init[1:16],
                          start, end, dt)

  for (i in 1:N) {
    # negative binomial likelihood
    nc[i] ~ dnegbin(prob = p_nc[i], size = r_nc)
    nh[i] ~ dnegbin(prob = p_nh[i], size = r_nh)
    dx[i] ~ dnegbin(prob = p_dx[i], size = r_dx)
    p_nc[i] <- r_nc / (r_nc + (mod_y[y_time[i], ind_nc]) + 1e-6)
    p_nh[i] <- r_nh / (r_nh + (mod_y[y_time[i], ind_nh]) + 1e-6)
    p_dx[i] <- r_dx / (r_dx + (mod_y[y_time[i], ind_dx]) + 1e-6)
    r_nc <- theta_nc
    r_nh <- theta_nh
    r_dx <- theta_dx
    }

  # dispersion for negative binomial
  theta_nc ~ dgamma(shape = 10, rate = 10/1000) # dgamma(b, b/a) where a is the mean and var is a^(2/b)
  theta_nh ~ dgamma(shape = 10, rate = 10/1000)
  theta_dx ~ dgamma(shape = 10, rate = 10/1000)

  for (z in 1:Z) {
    prd_nc[z] ~ dnegbin(prd_p_nc[z], prd_r_nc)
    prd_p_nc[z] <- prd_r_nc / (prd_r_nc + mod_y[iter_time[z], ind_nc] + 1e-6)
    prd_r_nc <- theta_nc

    prd_nh[z] ~ dnegbin(prd_p_nh[z], prd_r_nh)
    prd_p_nh[z] <- prd_r_nh / (prd_r_nh + mod_y[iter_time[z], ind_nh] + 1e-6)
    prd_r_nh <- theta_nh

    prd_dx[z] ~ dnegbin(prd_p_dx[z], prd_r_dx)
    prd_p_dx[z] <- prd_r_dx / (prd_r_dx + mod_y[iter_time[z], ind_dx] + 1e-6)
    prd_r_dx <- theta_dx

    prd_rx[z] <- mod_y[iter_time[z], ind_rx]
    beta_iter[z] <- mod_y[iter_time[z], ind_b]
  }
})

# 0.7 / (-SEIR_Const$infa * SEIR_Const$p_symp + SEIR_Const$infa + SEIR_Const$p_symp)
b_inits <- log(c(runif(1, min = 2, max = 3),
             runif(1, min = 1, max = 2),
             runif(SEIR_Const$n_beta - 2, min = 0.5, max = 1)) / (1 / SEIR_Const[["eps"]] + 1 / SEIR_Const[["sigma_c"]]) )
SEAIR_Inits <-  list(b = b_inits,
                   import_t0 = runif(1, min = 0, max = SEIR_Const$max_t0),
                   theta_nc = runif(1, min = 10, max = 1000),
                   theta_nh = runif(1, min = 10, max = 1000),
                   theta_dx = runif(1, min = 10, max = 1000),
                   sd_rw = runif(1, min = 0.2, max = 0.8))

SEAIR_Model <- nimble::nimbleModel(SEAIR_Code, SEIR_Const, SEIR_Data, inits = SEAIR_Inits)

# we compile the model
cSEAIR <- nimble::compileNimble(SEAIR_Model)

SEAIRConf <- nimble::configureMCMC(cSEAIR, thin = n_thin, autoblock = TRUE)
SEAIRConf$addMonitors(c("b", "import_t0",
                       "sd_rw", "b_last",
                       "theta_nc", "theta_nh", "theta_dx",
                       "prd_nc", "prd_nh", "prd_dx", "prd_rx", "beta_iter"))
SEAIRConf$removeSamplers(c("b", "sd_rw"))
SEAIRConf$addSampler(target = c("b", "sd_rw"), type = 'AF_slice')


SEAIRMCMC <- nimble::buildMCMC(SEAIRConf)
CSEAIRMCMC <- nimble::compileNimble(SEAIRMCMC, project = SEAIR_Model, resetFunctions = TRUE)

fit <- nimble::runMCMC(CSEAIRMCMC, niter = n_iterations, nchains = n_chains,
                   nburnin = warm_up, samplesAsCodaMCMC = TRUE)

return(fit)

}

get_outcomes_seair <- function(fit, SEIR_Data, SEIR_Const) {

  trace <- as.data.frame(do.call(rbind, fit))

  # new cases
  nom <- "prd_nc"
    med <- apply(trace[, grepl(nom, names(trace))], MARGIN = 2, median)
    lci <- apply(trace[, grepl(nom, names(trace))], MARGIN = 2, function(x) quantile(x, probs = c(0.025)))
    uci <- apply(trace[, grepl(nom, names(trace))], MARGIN = 2, function(x) quantile(x, probs = c(0.975)))
    l50 <- apply(trace[, grepl(nom, names(trace))], MARGIN = 2, function(x) quantile(x, probs = c(0.25)))
    u50 <- apply(trace[, grepl(nom, names(trace))], MARGIN = 2, function(x) quantile(x, probs = c(0.75)))
    nc <- data.frame(median = med, lci = lci, uci = uci, l50 = l50, u50 = u50)

  # number daily new hospitalization
  nom <- "prd_nh"
    med <- apply(trace[, grepl(nom, names(trace))], MARGIN = 2, median)
    lci <- apply(trace[, grepl(nom, names(trace))], MARGIN = 2, function(x) quantile(x, probs = c(0.025)))
    uci <- apply(trace[, grepl(nom, names(trace))], MARGIN = 2, function(x) quantile(x, probs = c(0.975)))
    l50 <- apply(trace[, grepl(nom, names(trace))], MARGIN = 2, function(x) quantile(x, probs = c(0.25)))
    u50 <- apply(trace[, grepl(nom, names(trace))], MARGIN = 2, function(x) quantile(x, probs = c(0.75)))
    nh <- data.frame(median = med, lci = lci, uci = uci, l50 = l50, u50 = u50)

  # deaths
  nom <- "prd_dx"
    med <- apply(trace[, grepl(nom, names(trace))], MARGIN = 2, median)
    lci <- apply(trace[, grepl(nom, names(trace))], MARGIN = 2, function(x) quantile(x, probs = c(0.025)))
    uci <- apply(trace[, grepl(nom, names(trace))], MARGIN = 2, function(x) quantile(x, probs = c(0.975)))
    l50 <- apply(trace[, grepl(nom, names(trace))], MARGIN = 2, function(x) quantile(x, probs = c(0.25)))
    u50 <- apply(trace[, grepl(nom, names(trace))], MARGIN = 2, function(x) quantile(x, probs = c(0.75)))
    dx <- data.frame(median = med, lci = lci, uci = uci, l50 = l50, u50 = u50)

  # total recovered
  nom <- "prd_rx"
    med <- apply(trace[, grepl(nom, names(trace))], MARGIN = 2, median)
    lci <- apply(trace[, grepl(nom, names(trace))], MARGIN = 2, function(x) quantile(x, probs = c(0.025)))
    uci <- apply(trace[, grepl(nom, names(trace))], MARGIN = 2, function(x) quantile(x, probs = c(0.975)))
    l50 <- apply(trace[, grepl(nom, names(trace))], MARGIN = 2, function(x) quantile(x, probs = c(0.25)))
    u50 <- apply(trace[, grepl(nom, names(trace))], MARGIN = 2, function(x) quantile(x, probs = c(0.75)))
    rx <- data.frame(median = med, lci = lci, uci = uci, l50 = l50, u50 = u50)


  # calculating R0
  duration_inf <- 1 / (1 / SEIR_Const[["sigma_c"]] + 1 / SEIR_Const[["eps"]])
  nom <- "beta_iter"
  trace_beta <- trace[, grepl(nom, names(trace))]
  trace_r0 <- SEIR_Const[["p_symp"]] * trace_beta / duration_inf +
              (1 - SEIR_Const[["p_symp"]]) * trace_beta * SEIR_Const[["infa"]] / duration_inf
    r0_med <- apply(trace_r0, MARGIN = 2, median)
    r0_lci <- apply(trace_r0, MARGIN = 2, function(x) quantile(x, probs = c(0.025)))
    r0_uci <- apply(trace_r0, MARGIN = 2, function(x) quantile(x, probs = c(0.975)))
    r0_l50 <- apply(trace_r0, MARGIN = 2, function(x) quantile(x, probs = c(0.25)))
    r0_u50 <- apply(trace_r0, MARGIN = 2, function(x) quantile(x, probs = c(0.75)))

  val <- list(nc = nc, nh = nh, dx = dx, rx = rx,
              r0_med = r0_med, r0_lci = r0_lci, r0_uci = r0_uci,
              r0_l50 = r0_l50, r0_u50 = r0_u50)
}

