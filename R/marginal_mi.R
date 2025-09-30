
#Clock reset at entry into illness state. Impute the sojourn time from illness to death using KMI where the
#risk set is everyone who has survived in the illness state beyond the censoring time.
#should be appropriate for semi-markov models, but not for extended semi-markov models

#' Layer 2 Imputation: Marginal Approach
#'
#' @param d A data frame with one row per subject and the columns event1, t1, event2, t2 on which mici::mici.impute was previously ran
#'
#' @returns A data frame with imputed times for event2 where event2 was censored
#' @export
marginal_mi <- function(d) {

  #Add an ID column
  d$id <- seq(1, nrow(d))

  # Create data for transition 1->2 (ill to death)
  u <- d[d$event1 == 1, ]  # All who became ill
  u$sojourn12 <- u$t2 - u$t1

  # Identify individuals need imputation for t2
  dc <- d[d$event2 == 1, ]
  dd <- d[d$event2 == 0, ]
  xt <- dd$t2 - dd$t1

  # Fit Kaplan-Meier for transition from ill to death
  km_summary <- summary(survival::survfit(survival::Surv(sojourn12, event2) ~ 1, data=u,timefix = FALSE))
  surv_probs <- km_summary$surv[length(km_summary$surv)]
  surv_times <- km_summary$time
  prob_diffs <- -diff(c(1, km_summary$surv))

  # Handle tail probability (if survival doesn't reach 0)
  if (surv_probs > 0) {
    prob_diffs <- c(prob_diffs, surv_probs)
    surv_times <- c(surv_times, max(u$sojourn12) + 1)
  }

  # Impute times for censored individuals
  cts <- NULL
  for (jj in 1:length(xt)) {
    # Find times greater than censoring time
    sub <- surv_times > xt[jj]
    # Sample time from illness to death
    if (sum(sub) > 1) {
      cts[jj] <- resample(surv_times[sub], 1,replace=TRUE, prob = prob_diffs[sub])
    } else if(sum(sub) == 1){
      cts[jj] <- surv_times[sub]
    } else {
      cts[jj] <- max(u$sojourn12) + 1 #shadow event time
    }
  }

  # Update death time and event indicator
  dd$event2 <- 1
  dd$t2 <- dd$t1 + cts
  ipd <- dplyr::bind_rows(dd, dc) %>% dplyr::arrange(id) %>% dplyr::select(-id)

  return(ipd)
}
