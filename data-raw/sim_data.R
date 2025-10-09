# Data Generation Process -------------------------------------------------

#Calculate scale parameters from Weibull distribution given median and shape
get_weibull_scale <- function(median, shape) {
  scale = log(2)/(median^shape);
  return(scale)
}

#Generate Weibull distributed survival times (no covariates)
generate_weibull <- function(n, shape, median, scale=NULL) {
  # Generate uniform random variables
  u <- runif(n, 0, 1)

  if (is.null(scale)) {
    #Use Weibull shape and desired median survival time to get scale
    scale = get_weibull_scale(median, shape)
  }
  # Generate survival times from Weibull distribution using inverse CDF
  t <- (-log(u)/scale)^(1/shape)

  return(t)
}

#Generate Cox-Weibull distributed survival times
generate_cox_weibull <- function(n, shape, median, X, beta) {
  # Generate uniform random variables
  u <- runif(n, 0, 1)
  #Use Weibull shape and desired median survival time to get scale
  scale = get_weibull_scale(median, shape)
  # Generate survival times from a Cox PH model with Weibull baseline hazard using inverse CDF
  t <- (-log(u)/(scale*exp(beta*X)))^(1/shape)
  return(t)
}

#Simulate from a multi-state model with three states: healthy (0), ill (1), dead (2)
#n- sample size (numeric)
#beta- effect (logHR) of T1 on lambda_12. When beta = 0, the MSM satisfies the Markov assumption. Default = 0. (numeric)
#shapeij- Weibull distribution shape for the i -> j transition (numeric > 0)
#medianij - median for the sojourn time between i -> j, used to calculate the Weibull scale for the i -> j transition (numeric > 0)
#return_latent_data- indicates if the function should return the underlying survival times in addition to the censored data that is observed. (logical)
#theta- upper bound for the uniform distribution for censoring (numeric)
simulate_illness_death <- function(n, beta = 0,
                                   shape01, median01,
                                   shape02, median02,
                                   shape12, median12,
                                   return_latent_data = FALSE,
                                   return_censored_data = TRUE,
                                   theta = NULL) {

  #Determine first event (illness or death)
  sojourn01 <- generate_weibull(n, shape01, median01)
  sojourn02 <- generate_weibull(n, shape02, median02)

  #Generate time from illness to death only for those with illness
  sojourn12 <- rep(NA, n)
  ill_ids <- which(sojourn01 < sojourn02)
  adjusted_scale_12 <- get_weibull_scale(median12, shape12) * exp(beta * sojourn01[ill_ids])
  sojourn12[ill_ids] <- generate_weibull(length(ill_ids), shape12, median = NULL, scale = adjusted_scale_12)

  temp = data.frame(sojourn01 = sojourn01, sojourn02 = sojourn02, sojourn12 = sojourn12)

  #Time of entry into each state and event indicators
  temp$id = 1:n
  temp$t1 = ifelse(temp$sojourn01 < temp$sojourn02, temp$sojourn01, temp$sojourn02)
  temp$event1 = ifelse(temp$sojourn01 < temp$sojourn02, 1, 0)
  temp$t2 = ifelse(temp$sojourn01 < temp$sojourn02, temp$sojourn01 + temp$sojourn12, temp$sojourn02)
  temp$event2 = 1
  temp[,c("id", "t1", "event1", "t2", "event2", "sojourn01", "sojourn02", "sojourn12")]

  #Add in uniform censoring times
  temp_censoring <- temp
  temp_censoring$C = runif(n, 0, theta)
  temp_censoring$event1 = ifelse(temp_censoring$C < temp_censoring$t1, 0, temp_censoring$event1)
  temp_censoring$t1 = ifelse(temp_censoring$C < temp_censoring$t1, temp_censoring$C, temp_censoring$t1)
  temp_censoring$event2 = ifelse(temp_censoring$C < temp_censoring$t2, 0, temp_censoring$event2)
  temp_censoring$t2 = ifelse(temp_censoring$C < temp_censoring$t2, temp_censoring$C, temp_censoring$t2)

  temp_censoring <- temp_censoring[, c("id", "t1", "event1", "t2", "event2", "sojourn01", "sojourn02", "sojourn12")]

  if (return_latent_data == FALSE & return_censored_data == TRUE) {
    return(d.observed = temp_censoring)
  } else if (return_censored_data == FALSE & return_latent_data == TRUE) {
    return(d.latent = temp)
  } else {
    return(list(d.observed = data.frame(temp_censoring), d.latent = data.frame(temp)))
  }

}
# Generate the dataset ------------------------------------------------
sim.data <- simulate_illness_death(n = 100, beta = log(1),
                                 shape01 = 1.3, median01 = 3,
                                 shape02 = 1.5, median02 = 8,
                                 shape12 = 2, median12 = 3,
                                 return_latent_data = FALSE,
                                 return_censored_data = TRUE,
                                 theta = 12)

usethis::use_data(sim.data, overwrite = TRUE)
