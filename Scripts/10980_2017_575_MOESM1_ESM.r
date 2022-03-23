# function for Bayesian latent indicator scale selection (BLISS)
# this function requires JAGS 4.2; and the R packages rjags and lubridate
run_BLISS = function(simulation_data, jags_file, niterations = 20000, burnin = 5000) {
  message(paste("run_BLISS( niterations =", niterations, ")"))

  start_time = lubridate::now()

  # initialize JAGS
  jags_model = rjags::jags.model(
    file = jags_file,
    data = c(simulation_data$jags_abundance_data, simulation_data$jags_detection_data),
    inits = simulation_data$jags_inits_abundance,
    n.chains = 1,
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
      "abundance_scale",
      "beta",
      "sigma.observer",
      "u.observer",
      "log_lambda",
      "logit_p_fixed",
      "logit_p_random",
      "p",
      "N"
    ),
    n.iter = niterations,
    thin = 1
  )

  # get selected scales
  selected_scales = rep(NA, simulation_data$covariate_data$ncovariates_abundance)
  for (i in 1:simulation_data$covariate_data$ncovariates_abundance) {
    tb_mcmc_scales_i = table(samples_jags$abundance_scale[i,,])

    selected_scales[i] = as.integer(names(which.max(tb_mcmc_scales_i)))
  }

  message("the selected scales are")
  print(selected_scales)

  finish_time = lubridate::now()

  # return result
  list(
    niterations = niterations,
    burnin = burnin,
    selected_scales = selected_scales,
    selected_model_results = samples_jags,
    start_time = start_time,
    finish_time = finish_time
  )
}


## run study
# load data
load("simulation_data_example.RData")

# run BLISS
BLISS_result = run_BLISS(
  simulation_data = simulation_data_example,
  jags_file = here::here("jags_example.txt",
  niterations = 20000,
  burnin = 5000
)

# extract some results
alpha_posterior_sample = data.frame(
  iteration = 1:BLISS_result$niterations,
  alpha = t(BLISS_result$selected_model_results$alpha[,,1])
)
scales_posterior_sample = data.frame(
  iteration = 1:BLISS_result$niterations,
  scale = t(BLISS_result$selected_model_results$abundance_scale[,,1])
)
N_posterior_sample = data.frame(
  iteration = 1:BLISS_result$niterations,
  N = t(BLISS_result$selected_model_results$N[,,1])
)

## create some figures (requires the dplyr, tidyr, ggplot2 packages)

# load dplyr
library(dplyr)

# alpha trace plots
alpha_posterior_sample %>%
  tidyr::gather("variable", "value", -iteration) %>%
  ggplot2::ggplot(ggplot2::aes(iteration, value)) +
  ggplot2::facet_wrap(~variable) +
  ggplot2::geom_line()

# alpha posterior density plots
alpha_posterior_sample %>%
  tidyr::gather("variable", "value", -iteration) %>%
  ggplot2::ggplot(ggplot2::aes(value)) +
  ggplot2::facet_wrap(~variable) +
  ggplot2::geom_density()

# scale trace plots
scales_posterior_sample %>%
  tidyr::gather("variable", "value", -iteration) %>%
  ggplot2::ggplot(ggplot2::aes(iteration, value)) +
  ggplot2::facet_wrap(~variable) +
  ggplot2::geom_line()

# scale histograms
scales_posterior_sample %>%
  tidyr::gather("variable", "value", -iteration) %>%
  ggplot2::ggplot(ggplot2::aes(value)) +
  ggplot2::facet_wrap(~variable) +
  ggplot2::geom_histogram(binwidth = 1)

# abundance trace plots
N_posterior_sample[,c(1, 1 + sort(base::sample(ncol(N_posterior_sample) - 1, 6)))] %>%
  tidyr::gather("variable", "value", -iteration) %>%
  ggplot2::ggplot(ggplot2::aes(iteration, value)) +
  ggplot2::facet_wrap(~variable) +
  ggplot2::geom_line()

# abundance histograms
N_posterior_sample[,c(1, 1 + sort(base::sample(ncol(N_posterior_sample) - 1, 6)))] %>%
  tidyr::gather("variable", "value", -iteration) %>%
  ggplot2::ggplot(ggplot2::aes(value)) +
  ggplot2::facet_wrap(~variable) +
  ggplot2::geom_histogram(binwidth = 1)
