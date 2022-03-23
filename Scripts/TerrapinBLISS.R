library(bayesplot) # Plotting for Bayesian Models
library(coda) # Output Analysis and Diagnostics for MCMC
library(ggridges) # Ridgeline Plots in 'ggplot2'
library(here) # A Simpler Way to Find Your Files
library(knitr) # A General-Purpose Package for Dynamic Report Generation in R
library(lubridate) # Make Dealing with Dates a Little Easier
library(posterior) # Tools for Working with Posterior Distributions
library(PresenceAbsence) # Presence-Absence Model Evaluation.
library(rjags) # Bayesian Graphical Models using MCMC
library(scales) # Scale Functions for Visualization
library(tidyverse) # Easily Install and Load the 'Tidyverse'


d <- read.csv(here::here("Terrapin_occupancy/Presence_both_years.csv"))
s1 <- str_split(d$Start_1, ":", simplify = T)
s1.d <- as.numeric(s1[,1]) + as.numeric(s1[,2])/60
s2 <- str_split(d$Start_2, ":", simplify = T)
s2.d <- as.numeric(s2[,1]) + as.numeric(s2[,2])/60
s3 <- str_split(d$Start_3, ":", simplify = T)
s3.d <- as.numeric(s3[,1]) + as.numeric(s3[,2])/60

# Prep the data
nSites <- nrow(d)
visits <- 3
y <- as.matrix(d[,3:5])
dimnames(y) <- NULL
start.time <- array(
  data = scale(c(s1.d, s2.d, s3.d)),
  dim = c(nSites, visits)
)
precip <- array(
  data = as.numeric(c(d$Precip1, d$Precip2, d$Precip3)),
  dim = c(nSites, visits)
)
ag <- as.matrix(d[,67:82]); dimnames(ag) <- NULL
ma <- as.matrix(d[,83:98]); dimnames(ma) <- NULL
gh <- as.matrix(d[,99:114]); dimnames(gh) <- NULL
lu <- as.matrix(d[,131:146]); dimnames(lu) <- NULL
ha <- as.matrix(d[,147:162]); dimnames(ha) <- NULL
cp <- as.matrix(d[,164]); dimnames(cp) <- NULL

ag <- apply(ag, 2, rescale)
ma <- apply(ma, 2, rescale)
gh <- apply(gh, 2, rescale)
lu <- apply(lu, 2, rescale)
ha <- apply(ha, 2, rescale)
cp <- apply(cp, 2, rescale)

sp.cov <- array(
  data = c(ag, ma, gh, lu, ha),
  dim = c(nSites, ncol(ag), 5),
  dimnames = NULL
)
dt.cov <- array(
  data = c(start.time, precip),
  dim = c(165, 3, 2)
)
det.hist <- as.matrix(d[,3:5]); dimnames(det.hist) <- NULL

jags.data2 <- list(
  jags_occurrence_data = list(
    nsites = 165,
    ncovariates_occurrence = 6,
    occurrence_covariates = sp.cov,
    occurrence_crabpots = cp
  ),
  jags_detection_data = list(
    ncovariates_detection = 2,
    nsurveys = 3,
    detection_covariates = dt.cov,
    yi = det.hist
  ),
  jags_inits_occurrence = list(
    list(y = as.vector(d$Presence)),
    list(y = as.vector(d$Presence)),
    list(y = as.vector(d$Presence)),
    list(y = as.vector(d$Presence))
  ),
  covariate_data = list(
    nsites = 165,
    nsurveys = 3,
    ncovariates_occurrence = 6,
    ncovariates_detection = 2,
    nscales_per_site = ncol(ag)
  )
)

occ_BLISS = function(data, jags_file, niterations = 20000, burnin = 5000, chains = 1) {
  message(paste("occ_BLISS( niterations =", niterations, ")"))
  
  start_time = lubridate::now()
  
  # initialize JAGS
  jags_model = rjags::jags.model(
    file = jags_file,
    data = c(data$jags_occurrence_data, data$jags_detection_data),
    inits = data$jags_inits_occurrence,
    n.chains = chains,
    n.adapt = max(100, ceiling(.1 * niterations)),
    quiet = FALSE
  )
  
  # burn-in
  if (burnin > 0) {
    message(paste("burn-in:", burnin, "iterations"))
    rjags::jags.samples(
      jags_model,
      variable.names = c("alpha"),
      n.iter = burnin,
      thin = 1
    )
  }
  
  # posterior simulation
  message(paste("main MCMC simulation:", niterations, " iterations"))
  samples_jags = rjags::jags.samples(
    jags_model,
    variable.names = c(
      "alpha",
      "occurrence_scale",
      "beta",
      "logit_psi",
      "logit_d_fixed",
      "y",
      "yi"
    ),
    n.iter = niterations,
    thin = 1
  )
  
  # get selected scales
  selected_scales = rep(NA, data$covariate_data$ncovariates_occurrence)
  for (i in 1:data$covariate_data$ncovariates_occurrence-1) {
    tb_mcmc_scales_i = table(samples_jags$occurrence_scale[i,,])
    
    selected_scales[i] = as.integer(names(which.max(tb_mcmc_scales_i)))
  }
  
  message("the selected scales are")
  print(selected_scales)
  
  finish_time = lubridate::now()
  
  message(paste0("Total elapsed time: ",
                 round(as.numeric(finish_time-start_time),2),
                 " minutes"))
  
  # return result
  list(
    niterations = niterations,
    burnin = burnin,
    nchains = chains,
    selected_scales = selected_scales,
    selected_model_results = samples_jags,
    start_time = start_time,
    finish_time = finish_time
  )
}

# run BLISS

# occ2 <- occ_BLISS(
#   data = jags.data2,
#   jags_file = here::here("Scripts/TerrapinBLISSJAGSModel2.txt"),
#   niterations = 3000,
#   burnin = 1000,
#   chains = 4
# )
# saveRDS(
#   object = occ2,
#   file = here::here("Output/Models/TerrapinBLISSModel.rds")
# )

occ2 <- readRDS(file = here::here("Output/Models/TerrapinBLISSModel.rds"))


# Assess model fit
bayesplot::mcmc_trace(
  coda::as.mcmc.list(occ2$selected_model_results$alpha),
) + ggtitle("Occurrence Coefficients")
bayesplot::mcmc_trace(
  coda::as.mcmc.list(occ2$selected_model_results$beta),
) + ggtitle("Detection Coefficients")
bayesplot::mcmc_trace(
  coda::as.mcmc.list(occ2$selected_model_results$occurrence_scale),
) + ggtitle("Best Scales")


# extract some results
alpha_posterior_sample = data.frame(
  chain = rep(1:occ2$nchains, each = occ2$niterations),
  iteration = rep(1:occ2$niterations, times = occ2$nchains),
  alpha = rbind(t(occ2$selected_model_results$alpha[,,1]),
                t(occ2$selected_model_results$alpha[,,2]),
                t(occ2$selected_model_results$alpha[,,3]),
                t(occ2$selected_model_results$alpha[,,4]))
)
names(alpha_posterior_sample) <- c(
  "chain", "iteration", "Ag", "Marsh", "Ghost Pots", 
  "Low Urban", "Armoring", "Active Pots"
)
scales_posterior_sample = data.frame(
  chain = rep(1:occ2$nchains, each = occ2$niterations),
  iteration = rep(1:occ2$niterations, times = occ2$nchains),
  scale = rbind(t(occ2$selected_model_results$occurrence_scale[,,1]),
                t(occ2$selected_model_results$occurrence_scale[,,2]),
                t(occ2$selected_model_results$occurrence_scale[,,3]),
                t(occ2$selected_model_results$occurrence_scale[,,4]))
)
names(scales_posterior_sample) <- c(
  "chain", "iteration", "Ag", "Marsh", "Ghost Pots", 
  "Low Urban", "Armoring"
)
site_posterior_sample = data.frame(
  chain = rep(1:occ2$nchains, each = occ2$niterations),
  iteration = rep(1:occ2$niterations, times = occ2$nchains),
  sites = rbind(t(occ2$selected_model_results$y[,,1]),
                t(occ2$selected_model_results$y[,,2]),
                t(occ2$selected_model_results$y[,,3]),
                t(occ2$selected_model_results$y[,,4]))
)
beta_posterior_sample = data.frame(
  chain = rep(1:occ2$nchains, each = occ2$niterations),
  iteration = rep(1:occ2$niterations, times = occ2$nchains),
  beta =  rbind(t(occ2$selected_model_results$beta[,,1]),
                t(occ2$selected_model_results$beta[,,2]),
                t(occ2$selected_model_results$beta[,,3]),
                t(occ2$selected_model_results$beta[,,4]))
)
names(beta_posterior_sample) <- c(
  "chain", "iteration", "Start Time", "Precipitation"
)

# Posterior summaries
posterior::summarise_draws(alpha_posterior_sample[,3:8])
posterior::summarise_draws(scales_posterior_sample[,3:7])
posterior::summarise_draws(beta_posterior_sample[,3:4])

# Plot the posteriors
## Occurrence
bayesplot::mcmc_areas(alpha_posterior_sample[,3:8], prob = 0.89) + 
  labs(x = "Posterior Coefficient Estimates")
## Detection
bayesplot::mcmc_areas(beta_posterior_sample[,3:4], prob = 0.89) + 
  labs(x = "Posterior Coefficient Estimates")
## Scale
bayesplot::mcmc_hist(scales_posterior_sample[,3:7], binwidth = 1) + 
  labs(x = "Posterior Scale Estimates")


# Posterior means
post <- data.frame(
  site = 1:165,
  Obs = d$Presence,
  Mean = apply(site_posterior_sample[,-c(1:2)], 2, mean)
)
rownames(post) <- NULL
knitr::kable(post, digits = 3)

# Optimal threshold for classification
opt.thresh <- PresenceAbsence::optimal.thresholds(
  DATA = post,
  threshold = 1001 # number of divisions between 0 and 1
)

# Confusion matrix
PresenceAbsence::cmx(
  DATA = post,
  threshold = opt.thresh$Mean[2] # Sens=Spec
)

# Add predicted classification
post2 <- data.frame(
  post,
  Class = ifelse(post$Mean < opt.thresh$Mean[2], 0, 1)
)
knitr::kable(post2, digits = 3)

# Summary Stats
sstats <- tibble(
  "Type" = c("Naive", "Bayesian", "Paper"),
  "Occupancy Estimate" = c(
    mean(post2$Obs),
    mean(post2$Class),
    0.440
  ),
  "Occupied Sites" = c(
    sum(post2$Obs),
    sum(post2$Class),
    73
  )
)
knitr::kable(sstats, digits = 3)

