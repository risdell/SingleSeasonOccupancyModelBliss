# Bayesian Single-season Occupancy Model example using BLISS

I have been wanting to update my occupancy modeling from the old frequentist technique to Bayesian. I recently found `nimble`, and that is great, but there was still the nagging discomfort about model selection when working with multiple spatial scales for various coefficients. I finally came across the BLISS method from [Stuber et al. 2017][2], and it was exactly what I was looking for! It uses a latent variable that selects a scale from the candidate set, and removes the issues associated with highly-correlated land-cover variables. The paper is a really great read. To try and figure out how their code works, I went back to an old data set of mine that I published [(Isdell et al. 2014)][5] using a single-season occupancy model implemented in `unmarked` ([Chandler et al. 2021][3]) in R. 

The data consist of three replicate presence/absence surveys for diamondback terrapins at 165 sites in the southern Chesapeake Bay. Precipitation (Y/N) and start time were recorded during each survey for the detection probability, and the number of active crab pot buoys within 250 m of a site were also counted. Spatial covariates were derived from GIS layers at multiple scales (0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 3, 4, 5, 6, 7, 8, 9, 10 km) for agricultural, marsh, and low intensity development, as well as ghost crab pots and shoreline armoring. Active pots were scaled across the three surveys (mean/sd) for a single site-estimate. Below are some of the highlights and modifications that I made to the original code from [Stuber et al. 2017][2] to adapt it to a single-season occupancy model running on multiple parallel chains in JAGS. Grab the repo for the full data and code.

## Data structure
```r
str(jags.data)
List of 4
 $ jags_occurrence_data :List of 4
  ..$ nsites                : num 165
  ..$ ncovariates_occurrence: num 6
  ..$ occurrence_covariates : num [1:165, 1:16, 1:5] 0 0 0 0 0.0316 ...
  ..$ occurrence_crabpots   : num [1:165, 1] 0 0.156 0.739 0 0.437 ...
 $ jags_detection_data  :List of 4
  ..$ ncovariates_detection: num 2
  ..$ nsurveys             : num 3
  ..$ detection_covariates : num [1:165, 1:3, 1:2] 2.739 2.432 2.24 2.703 -0.407 ...
  ..$ yi                   : int [1:165, 1:3] 0 0 0 0 0 1 0 0 0 0 ...
 $ jags_inits_occurrence:List of 4
  ..$ :List of 1
  .. ..$ y: int [1:165] 0 1 0 0 0 1 0 0 0 0 ...
  ..$ :List of 1
  .. ..$ y: int [1:165] 0 1 0 0 0 1 0 0 0 0 ...
  ..$ :List of 1
  .. ..$ y: int [1:165] 0 1 0 0 0 1 0 0 0 0 ...
  ..$ :List of 1
  .. ..$ y: int [1:165] 0 1 0 0 0 1 0 0 0 0 ...
 $ covariate_data       :List of 5
  ..$ nsites                : num 165
  ..$ nsurveys              : num 3
  ..$ ncovariates_occurrence: num 6
  ..$ ncovariates_detection : num 2
  ..$ nscales_per_site      : int 16
```

Above, the `occurrence_covariates` is an array of the values of the five spatial covariates at 16 scales for each of the 165 sites. Because the active pots were derived from on-site counts rather than a GIS layer, they are included separately as `occurrence_crabpots`. Similarly, the `detection_covariates` is an array of start times (centered and scaled) and precipitation records per site per survey, and `yi` is the detection history for each site. The `jags_inits_occurrence` list is a list of initial values supplied to JAGS consisting of the recorded presence values at the site level (i.e., if a terrapin was sighted on any survey, the site recieves a 1). There is one list per parallel chain.

## JAGS model
```jags
model {
    prec_coef <- 1

    # Priors for covariates of ecological model
    for (i in 1:(ncovariates_occurrence)) {
        alpha[i] ~ dnorm(0.0, prec_coef)
    }
    alpha.int ~ dnorm(0, 1)

    occurrence_scale[1] ~ dcat(c(0.0625, 0.0625, 0.0625, 0.0625, 0.0625, 
                                 0.0625, 0.0625, 0.0625, 0.0625, 0.0625, 
                                 0.0625, 0.0625, 0.0625, 0.0625, 0.0625, 
                                 0.0625))
    occurrence_scale[2] ~ dcat(c(0.0625, 0.0625, 0.0625, 0.0625, 0.0625, 
                                 0.0625, 0.0625, 0.0625, 0.0625, 0.0625, 
                                 0.0625, 0.0625, 0.0625, 0.0625, 0.0625, 
                                 0.0625))
    occurrence_scale[3] ~ dcat(c(0.0625, 0.0625, 0.0625, 0.0625, 0.0625, 
                                 0.0625, 0.0625, 0.0625, 0.0625, 0.0625, 
                                 0.0625, 0.0625, 0.0625, 0.0625, 0.0625, 
                                 0.0625))                        
    occurrence_scale[4] ~ dcat(c(0.0625, 0.0625, 0.0625, 0.0625, 0.0625, 
                                 0.0625, 0.0625, 0.0625, 0.0625, 0.0625, 
                                 0.0625, 0.0625, 0.0625, 0.0625, 0.0625, 
                                 0.0625))
    occurrence_scale[5] ~ dcat(c(0.0625, 0.0625, 0.0625, 0.0625, 0.0625, 
                                 0.0625, 0.0625, 0.0625, 0.0625, 0.0625, 
                                 0.0625, 0.0625, 0.0625, 0.0625, 0.0625, 
                                 0.0625))                            
    # Priors for covariates of observation model
    for (i in 1:ncovariates_detection) {
        beta[i] ~ dnorm(0.0, prec_coef)
    }
    beta.int ~ dnorm(0, 1)

    # Random effects for the observation model
    # sigma.observer ~ dexp(1)
    # tau.observer <- 1 / (sigma.observer * sigma.observer)
    # for (o in 1:nobservers) {
    #     u.observer[o] ~ dnorm(0, tau.observer)
    # }

    # Ecological model
    for (s in 1:nsites) {
      # Ecological model
        logit_psi[s] <- alpha.int + 
        alpha[1] * occurrence_covariates[s, occurrence_scale[1], 1] + 
        alpha[2] * occurrence_covariates[s, occurrence_scale[2], 2] + 
        alpha[3] * occurrence_covariates[s, occurrence_scale[3], 3] + 
        alpha[4] * occurrence_covariates[s, occurrence_scale[4], 4] + 
        alpha[5] * occurrence_covariates[s, occurrence_scale[5], 5] + 
        alpha[6] * occurrence_crabpots[s, 1]
        # Set limits
        logit_psi.lim[s] <- min(999,max(-999,logit_psi[s]))
        psi[s] <- 1/(1 + exp(-logit_psi.lim[s]))

        y[s] ~ dbern(psi[s])
	# Observation model
    	for (i in 1:nsurveys) {
          logit_p[s,i] <- beta.int + beta %*% detection_covariates[s,i,]
          logit_p.lim[s,i] <- min(999,max(-999,logit_p[s,i]))
          p[s,i] <- 1/(1 + exp(-logit_p.lim[s,i]))
          yi[s,i] ~ dbern(p[s,i]*y[s])
      
      presid[s,i] <- (yi[s,i] - p[s,i])^2 # squared residual of obs dat
      y.new[s,i] ~ dbern(p[s,i]*y[s]) # sim data
      presid.new[s,i] <- (y.new[s,i] - p[s,i])^2 # sq resid of sim data
    	}
    } 
  mean.p <- exp(beta.int)/(1 + exp(beta.int))  # mean detection probability
  occ.sites <- sum(y[]) # number of occupied sites per iteration
}
```
Using the same methods as the original BLISS paper ([Stuber et al. 2017][2]), uniform priors were set for each of the candidate scales. Observer was excluded from this model because there was a single observer for the entire study, but see [here][4] for the original version of the JAGS model and how to incoporate multiple observers. The model was further modified to include the detection as a sub-model of the occurrence model. 

## BLISS function

```r
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
      "alpha", "occurrence_scale", "beta", "logit_psi", 
      "y", "yi", "p", "presid", "y.new", "presid.new", 
      "mean.p", "occ.sites"
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
```
Only minor modifications were made to the [Stuber et al. 2017][2] BLISS function. Variable names were updated to reflect the model, and the number of chains was adapted to allow for multiple chains.

## Run the model
```r
occ <- occ_BLISS(
  data = jags.data2,
  jags_file = here::here("Scripts/TerrapinBLISSJAGSModel3.txt"),
  niterations = 3000,
  burnin = 1000,
  chains = 4
)
```

## Check fit
Once you have results, you can assess the model fit.
```r
bayesplot::mcmc_trace(
  coda::as.mcmc.list(occ$selected_model_results$alpha),
) + ggtitle("Occurrence Coefficients")
bayesplot::mcmc_trace(
  coda::as.mcmc.list(occ$selected_model_results$beta),
) + ggtitle("Detection Coefficients")
bayesplot::mcmc_trace(
  coda::as.mcmc.list(occ$selected_model_results$occurrence_scale),
) + ggtitle("Best Scales")
```
![ExampleTracePlot](https://user-images.githubusercontent.com/48960489/159991269-26579e76-5cad-4f13-9ea4-e8bd06ec2d94.jpg)


## Extract the posteriors
Extract the draws for plots and summaries. I'm positive there's a more efficient way to do this, but, as a quick side project, I just adapted the methods from the original paper.
```r
# Occurrence covariates
alpha_posterior_sample = data.frame(
  chain = rep(1:occ$nchains, each = occ$niterations),
  iteration = rep(1:occ$niterations, times = occ$nchains),
  alpha = rbind(t(occ$selected_model_results$alpha[,,1]),
                t(occ$selected_model_results$alpha[,,2]),
                t(occ$selected_model_results$alpha[,,3]),
                t(occ$selected_model_results$alpha[,,4]))
)
names(alpha_posterior_sample) <- c(
  "chain", "iteration", "Ag", "Marsh", "Ghost Pots", 
  "Low Urban", "Armoring", "Active Pots"
)
# Scales
scales_posterior_sample = data.frame(
  chain = rep(1:occ$nchains, each = occ$niterations),
  iteration = rep(1:occ$niterations, times = occ$nchains),
  scale = rbind(t(occ$selected_model_results$occurrence_scale[,,1]),
                t(occ$selected_model_results$occurrence_scale[,,2]),
                t(occ$selected_model_results$occurrence_scale[,,3]),
                t(occ$selected_model_results$occurrence_scale[,,4]))
)
names(scales_posterior_sample) <- c(
  "chain", "iteration", "Ag", "Marsh", "Ghost Pots", 
  "Low Urban", "Armoring"
)
# Detection covariates
beta_posterior_sample = data.frame(
  chain = rep(1:occ$nchains, each = occ$niterations),
  iteration = rep(1:occ$niterations, times = occ$nchains),
  beta =  rbind(t(occ$selected_model_results$beta[,,1]),
                t(occ$selected_model_results$beta[,,2]),
                t(occ$selected_model_results$beta[,,3]),
                t(occ$selected_model_results$beta[,,4]))
)
names(beta_posterior_sample) <- c(
  "chain", "iteration", "Start Time", "Precipitation"
)
# Number of occupied sites
total_sites_posterior_sample = data.frame(
  chain = rep(1:occ$nchains, each = occ$niterations),
  iteration = rep(1:occ$niterations, times = occ$nchains),
  total = c(as.vector(occ$selected_model_results$occ.sites[,,1]),
            as.vector(occ$selected_model_results$occ.sites[,,2]),
            as.vector(occ$selected_model_results$occ.sites[,,3]),
            as.vector(occ$selected_model_results$occ.sites[,,4]))
)
# Detection probability
mean_p_posterior_sample = data.frame(
  chain = rep(1:occ$nchains, each = occ$niterations),
  iteration = rep(1:occ$niterations, times = occ$nchains),
  prob = c(as.vector(occ$selected_model_results$mean.p[,,1]),
            as.vector(occ$selected_model_results$mean.p[,,2]),
            as.vector(occ$selected_model_results$mean.p[,,3]),
            as.vector(occ$selected_model_results$mean.p[,,4]))
)
```

## Quick summaries
```r
knitr::kable(posterior::summarise_draws(alpha_posterior_sample[,3:8]), digits = 3)
knitr::kable(posterior::summarise_draws(scales_posterior_sample[,3:7]), digits = 3)
knitr::kable(posterior::summarise_draws(beta_posterior_sample[,3:4]), digits = 3)
```

|variable    |   mean| median|    sd|   mad|     q5|    q95|  rhat| ess_bulk| ess_tail|
|:-----------|------:|------:|-----:|-----:|------:|------:|-----:|--------:|--------:|
|Ag          | -1.400| -1.382| 0.775| 0.788| -2.708| -0.153| 1.001| 2099.921| 4258.567|
|Marsh       |  3.163|  3.151| 0.762| 0.773|  1.935|  4.435| 1.001| 2194.200| 4208.497|
|Ghost Pots  | -0.445| -0.465| 0.812| 0.811| -1.761|  0.906| 1.000| 2579.304| 4730.546|
|Low Urban   | -0.295| -0.283| 0.856| 0.848| -1.730|  1.100| 1.000| 2607.460| 5161.259|
|Armoring    | -2.073| -2.066| 0.816| 0.808| -3.426| -0.729| 1.001| 2301.367| 3698.624|
|Active Pots | -1.398| -1.400| 0.750| 0.747| -2.640| -0.158| 1.001| 3345.458| 5629.992|

## Plot the posteriors
```r
## Occurrence
bayesplot::mcmc_areas(alpha_posterior_sample[,3:8], prob = 0.89) + 
  labs(x = "Posterior Coefficient Estimates")
## Detection
bayesplot::mcmc_areas(beta_posterior_sample[,3:4], prob = 0.89) + 
  labs(x = "Posterior Coefficient Estimates")
## Scale
bayesplot::mcmc_hist(scales_posterior_sample[,3:7], binwidth = 1) + 
  labs(x = "Posterior Scale Estimates")
```
![ExamplePosteriorPlot](https://user-images.githubusercontent.com/48960489/159991906-82bc30c7-8dd1-4e76-8921-e28dc96ad7c5.jpg)

And here is the really cool part about BLISS. You can see how the posterior is distributed across the scales for each parameter. In the plots below, the x-axis is the index of each scale (1-16). The BLISS function returns which scale recieves the highest density, but it operates over all of them. 
![ExampleScalesPlot](https://user-images.githubusercontent.com/48960489/159992335-e3a0feb9-e08b-4569-8244-4e66c69ac478.jpg)


The posteriors for number of occupied sites and detection probability are always useful.
![ExampleSitesAndDetProbPosteriors](https://user-images.githubusercontent.com/48960489/159991836-248e7f56-38dd-4072-870a-5d5d40d924f6.jpg)

With your model in place, apply it to your study area.
![ExampleTerrapinMap](https://github.com/risdell/SingleSeasonOccupancyModelBliss/blob/8d8c40711f066b0620c4bcefb93b18fc96a2817a/Output/Figures/TerrapinOccupancy.png)

There's plenty more to do, but this should be enough to help you get started. Download the repo for the full analysis.


[1]: <https://books.google.com/books?id=pf-w-JAUOd0C&hl=en> 
[2]: <https://doi.org/10.1007/s10980-017-0575-y> 
[3]: <https://cran.r-project.org/web/packages/unmarked/index.html>
[4]: <https://static-content.springer.com/esm/art%3A10.1007%2Fs10980-017-0575-y/MediaObjects/10980_2017_575_MOESM2_ESM.txt>
[5]: <https://onlinelibrary.wiley.com/doi/full/10.1111/ddi.12289>