#' Layer 2 Imputation: Cox Model (Conditional) Approach
#'
#' @param d A data frame with one row per subject and the columns event1, t1, event2, t2 on which mici::mici.impute was previously ran
#'
#' @returns A data frame with imputed times for event2 where event2 was censored
#' @export
cox_mi <- function(d) {

  #create ID column
  d$id <- seq(1:nrow(d))

  # Create data for transition 1->2 (ill to death)
  u <- d[d$event1 == 1, ]  # All who became ill
  u$sojourn12 <- u$t2 - u$t1

  # Identify individuals need imputation for t2
  dc <- d[d$event2 == 1, ]
  dd <- d[d$event2 == 0, ]
  xt <- dd$t2 - dd$t1

  #Estimate Survival function for ill to death sojourn time using coxph model with t1 as a covariate
  #each person has their own survival curve based on their time to illness
  cox_model <- survival::coxph(survival::Surv(t2-t1, event2) ~ t1, data=u)
  surv_summary <- summary(survival::survfit(cox_model, newdata = dd))
  surv_probs <- surv_summary$surv[dim(surv_summary$surv)[1], ]
  surv_times <- surv_summary$time
  prob_diffs <- apply(surv_summary$surv, 2, function(x) -diff(c(1, x)))

  # Handle tail probability (if survival doesn't reach 0)
  if (max(surv_probs) > 0) {
    prob_diffs <- rbind(prob_diffs, surv_probs)
    surv_times <- c(surv_times, max(u$sojourn12) + 1)
  }

  # Impute times for censored individuals
  cts <- NULL
  for (jj in 1:length(xt)) {
    # Find times greater than censoring time
    sub <- surv_times > xt[jj]
    # Sample time from illness to death
    if (sum(sub) > 1) {
      cts[jj] <- resample(surv_times[sub], 1,replace=TRUE, prob = prob_diffs[sub, jj])
    } else if(sum(sub) == 1){
      cts[jj] <- surv_times[sub]
    } else {
      cts[jj] <- max(u$sojourn12) + 1 #shadow event time
    }
  }
  # Update data
  dd$event2 <- 1
  dd$t2 <- dd$t1 + cts
  ipd <- dplyr::bind_rows(dd, dc) %>% dplyr::arrange(id) %>% dplyr::select(-id)

  return(ipd)
}

