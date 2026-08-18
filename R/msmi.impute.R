# Cox Imputation ----------------------------------------------------------

#' Layer 2 Imputation: Cox Model (Conditional) Approach
#'
#' @param d A data frame with one row per subject and the columns event1, t1, event2, t2 on which mici::mici.impute was previously ran
#'
#' @returns A data frame with imputed times for event2 where event2 was censored
cox_mi <- function(d, bootstrap = TRUE) {
  #create ID column
  d$id <- seq_along(1:nrow(d))

  # Create data for transition 1->2 (ill to death)
  u <- d[d$event1 == 1, ] # All who became ill
  u$sojourn23 <- u$t2 - u$t1

  # Identify individuals need imputation for t2
  dc <- d[d$event2 == 1, ]
  dd <- d[d$event2 == 0, ]
  xt <- dd$t2 - dd$t1

  # Guard: nothing to impute for this iteration
  if (nrow(dd) == 0) {
    return(d %>% dplyr::select(-id))
  }

  #Take boostrap sample
  if (bootstrap) {
    boot <- u[sample(1:nrow(u), replace = TRUE), ]
  } else {
    boot <- u
  }

  #Define shadow time
  t_shadow <- max(u$sojourn23) + 1

  #Gaurd: bootstrap sample has no observed transitions: impute a shadow time
  if (nrow(boot[boot$event2 == 1, ]) == 0) {
    dd$event2 <- 1
    dd$t2 <- dd$t1 + t_shadow

    ipd <- dplyr::bind_rows(dd, dc) %>%
      dplyr::arrange(id) %>%
      dplyr::select(-id)

    return(ipd)
  } else if (nrow(boot[boot$event2 == 1, ]) <= 5) {
    warning(
      "There are very few uncensored observations after the first round of imputation for this dataset.
                Cox imputation models may have trouble converging. Changing to the marginal imputation method instead."
    )

    km_summary <- summary(survival::survfit(
      survival::Surv(sojourn23, event2) ~ 1,
      data = boot,
      timefix = FALSE
    ))
    surv_probs <- km_summary$surv[length(km_summary$surv)]
    surv_times <- km_summary$time
    prob_diffs <- -diff(c(1, km_summary$surv))

    # Handle tail probability (if survival doesn't reach 0)
    if (surv_probs > 0) {
      prob_diffs <- c(prob_diffs, surv_probs)
      surv_times <- c(surv_times, t_shadow) #shadow time from original dataset (not boostrap sample)
    }

    # Impute times for censored individuals
    cts <- NULL
    for (jj in seq_along(xt)) {
      # Find times greater than censoring time
      sub <- surv_times > xt[jj]

      # Sample time from illness to death
      if (sum(sub) > 1) {
        cts[jj] <- resample(
          surv_times[sub],
          1,
          replace = TRUE,
          prob = prob_diffs[sub]
        )
      } else if (sum(sub) == 1) {
        cts[jj] <- surv_times[sub]
      } else {
        cts[jj] <- t_shadow
      }
    }

    # Update death time and event indicator
    dd$event2 <- 1
    dd$t2 <- dd$t1 + cts
    ipd <- dplyr::bind_rows(dd, dc) %>%
      dplyr::arrange(id) %>%
      dplyr::select(-id)

    return(ipd)
  } else {
    #Estimate Survival function for ill to death sojourn time using coxph model with t1 as a covariate
    #each person has their own survival curve based on their time to illness
    cox_model <- survival::coxph(
      survival::Surv(t2 - t1, event2) ~ t1,
      data = boot
    )
    surv_summary <- summary(survival::survfit(cox_model, newdata = dd))

    if (nrow(dd) > 1) {
      surv_probs <- surv_summary$surv[dim(surv_summary$surv)[1], ]
    } else if (nrow(dd) == 1) {
      surv_probs <- surv_summary$surv[length(surv_summary$surv)]
    }

    surv_times <- surv_summary$time
    surv_times_list <- replicate(
      length(surv_probs),
      surv_times,
      simplify = FALSE
    )

    if (nrow(dd) > 1) {
      prob_diffs <- apply(surv_summary$surv, 2, function(x) -diff(c(1, x)))
      prob_diffs_list <- split(prob_diffs, col(prob_diffs))
    } else if (nrow(dd) == 1) {
      prob_diffs <- -diff(c(1, surv_summary$surv))
      prob_diffs_list <- list(prob_diffs)
    }

    # Handle tail probability (if survival doesn't reach 0)
    tail.probs <- function(tail, probs) {
      if (tail > 0) {
        probs <- c(probs, tail)
      }
      return(probs)
    }

    tail.times <- function(tail, times) {
      if (tail > 0) {
        times <- c(times, t_shadow)
      }
      return(times)
    }

    prob_diffs <- Map(tail.probs, tail = surv_probs, probs = prob_diffs_list) #output is a list of length nrow(dd)

    surv_times <- Map(tail.times, surv_probs, surv_times_list) #output is a list of length nrow(dd)

    # Impute times for censored individuals
    cts <- NULL
    for (jj in seq_along(xt)) {
      # Find times greater than censoring time
      sub <- surv_times[[jj]] > xt[jj]

      #ensure at least one positive probability
      if (sum(prob_diffs[[jj]][sub]) == 0) {
        prob_diffs[[jj]][sub][1] <- 1
      }

      # Sample time from illness to death
      if (sum(sub) > 1) {
        cts[jj] <- resample(
          surv_times[[jj]][sub],
          1,
          replace = TRUE,
          prob = prob_diffs[[jj]][sub]
        )
      } else if (sum(sub) == 1) {
        cts[jj] <- surv_times[[jj]][sub]
      } else {
        cts[jj] <- t_shadow
      }
    }
    # Update data
    dd$event2 <- 1
    dd$t2 <- dd$t1 + cts
    ipd <- dplyr::bind_rows(dd, dc) %>%
      dplyr::arrange(id) %>%
      dplyr::select(-id)

    return(ipd)
  }
}


# Marginal Imputation -----------------------------------------------------
#Idea use clock forward timescale instead of reset?

#' Layer 2 Imputation: Marginal Approach
#'
#' @param d A data frame with one row per subject and the columns event1, t1, event2, t2 on which mici::mici.impute was previously ran

#' @returns A data frame with imputed times for event2 where event2 was censored
marginal_mi <- function(d, bootstrap = TRUE) {
  #Add an ID column
  d$id <- seq(1, nrow(d))

  # Create data for transition 1->2 (ill to death)
  u <- d[d$event1 == 1, ] # All who became ill
  u$sojourn23 <- u$t2 - u$t1

  # Identify individuals need imputation for t2
  dc <- d[d$event2 == 1, ]
  dd <- d[d$event2 == 0, ]
  xt <- dd$t2 - dd$t1

  # Guard: nothing to impute for this iteration
  if (nrow(dd) == 0) {
    return(d %>% dplyr::select(-id))
  }

  #random boostrap sample of the original dataset to use for risk set construction
  if (bootstrap) {
    boot <- u[sample(1:nrow(u), replace = TRUE), ]
  } else {
    boot <- u
  }

  # Impute shadow time if no observed transitions
  t_shadow <- max(u$sojourn23) + 1
  if (nrow(boot[boot$event2 == 1, ]) == 0) {
    dd$event2 <- 1
    dd$t2 <- dd$t1 + t_shadow

    ipd <- dplyr::bind_rows(dd, dc) %>%
      dplyr::arrange(id) %>%
      dplyr::select(-id)

    return(ipd)
  } else {
    # Fit Kaplan-Meier for transition from ill to death
    km_summary <- summary(survival::survfit(
      survival::Surv(sojourn23, event2) ~ 1,
      data = boot,
      timefix = FALSE
    ))
    surv_probs <- km_summary$surv[length(km_summary$surv)]
    surv_times <- km_summary$time
    prob_diffs <- -diff(c(1, km_summary$surv))

    # Handle tail probability (if survival doesn't reach 0)
    if (surv_probs > 0) {
      prob_diffs <- c(prob_diffs, surv_probs)
      surv_times <- c(surv_times, t_shadow) #shadow time from original dataset (not boostrap sample)
    }

    # Impute times for censored individuals
    cts <- NULL
    for (jj in seq_along(xt)) {
      # Find times greater than censoring time
      sub <- surv_times > xt[jj]

      # Sample time from illness to death
      if (sum(sub) > 1) {
        cts[jj] <- resample(
          surv_times[sub],
          1,
          replace = TRUE,
          prob = prob_diffs[sub]
        )
      } else if (sum(sub) == 1) {
        cts[jj] <- surv_times[sub]
      } else {
        cts[jj] <- t_shadow
      }
    }

    # Update death time and event indicator
    dd$event2 <- 1
    dd$t2 <- dd$t1 + cts
    ipd <- dplyr::bind_rows(dd, dc) %>%
      dplyr::arrange(id) %>%
      dplyr::select(-id)

    return(ipd)
  }
}

# Wrapper Function --------------------------------------------------------

#Provide a wrapper function that seamlessly transitions between the first and second imputation layers
#and allows the user to choose which method to use for the second layer

#' Create multiple imputed datasets for data arising from multi-state models subject to censoring
#'
#' @param dat a dataframe with one row per subject and columns corresponding to the time and event indicators for each state transition.
#'  The column names for the time and event indicators should be in the format specified by prefix.states, specifically
#'  event<i> and t<i> for i = 1, ..., n.states-1
#' @param M an integer, the number of imputations
#' @param prefix.states a character vector of length 2, specify the prefix for the event and time columns in the data in that order
#' @param method a character string, either "marginal" or "cox" indicating which method to use for the second layer of imputation
#' @param seed an integer, the seed for random number generation
#' @examples
#' msmi.impute(sim.data, M = 5, prefix.states = c("event", "t"), method = "marginal")
#' @returns A list of length M, where each element is a data frame with imputed times for censored events
#' @export
msmi.impute <- function(
  dat,
  M,
  prefix.states = c("event", "t"),
  method = "marginal",
  seed = sample(1:.Machine$integer.max, size = 1),
  bootstrap = TRUE
) {
  #Check inputs
  if (length(prefix.states) != 2) {
    stop("prefix.states must be a character vector of length 2")
  }
  if (method != "marginal" & method != "cox") {
    stop("method must be either 'marginal' or 'cox'")
  }
  if (!(M == floor(M))) {
    stop("M must be an integer")
  }

  #set seed
  set.seed(seed)

  #Create standardized dataset d from user provided dataframe dat
  d <- data.frame(row.names = 1:nrow(dat))

  for (i in 1:2) {
    if (
      !all(
        c(paste0(prefix.states[1], i), paste0(prefix.states[2], i)) %in%
          colnames(dat)
      )
    ) {
      stop(paste0(
        "Columns ",
        paste0(prefix.states[1], i),
        " and ",
        paste0(prefix.states[2], i),
        " must be present in the data"
      ))
    }
    d[paste0("event", i)] <- dat[, paste0(prefix.states[1], i)]
    d[paste0("t", i)] <- dat[, paste0(prefix.states[2], i)]
  }
  #d now has columns event1, t1, event2, t2

  #prepare data in competing events structure for time to first event
  d.comp <- d %>%
    dplyr::mutate(
      t.first = pmin(t1, t2),
      event.first = dplyr::case_when(
        t1 < t2 & event1 == 1 ~ 1,
        !(t1 < t2) & event2 == 1 ~ 2,
        TRUE ~ 0
      )
    )
  #multiple imputations for time to first event using mici.impute
  d.imp1 <- mici.impute(
    "t.first",
    "event.first",
    data = d.comp,
    M = M,
    bootstrap
  )

  #put data back into the original format for second layer of imputation
  d.imp1 <- purrr::map(d.imp1, function(x) {
    idx2 <- x$ftype == 2 & x$event.first == 0 #indices for originally censored people who were imputed to have event type 2 first
    idx1 <- x$ftype == 1 & x$event.first == 0 #indices for originally censored people who were imputed to have event type 1 first

    x$event1[idx2] <- 0 #those imputed to have event type 2 first did not have event 1
    x$event1[idx1] <- 1 #those imputed to have event type 1 first did have event 1
    x$t1[idx2 | idx1] <- x$ftime[idx2 | idx1] #everyone who was censored and has their t1 updated to the imputed time
    x$event2[idx2] <- 1 #those imputed to have event type 2 first did have event 2
    x$event2[idx1] <- 0 #those imputed to have event type 1 first did not have event 2
    x$t2[idx2 | idx1] <- x$ftime[idx2 | idx1] #everyone who was censored has their t2 updated to the imputed time

    x %>% dplyr::select(t1, event1, t2, event2)
  })

  #second layer of imputation for time to second event

  dd <- d[d$event2 == 0, ]
  if (nrow(dd) == 0) {
    warning(
      "There was no censoring in the second layer of imputation for this dataset.
            Imputation was not needed."
    )
    return(d.imp1)
  }

  u <- d[d$event1 == 1 & d$event2 == 1, ]
  if (nrow(u) == 0) {
    warning(
      "Note: There are no uncensored transitions between the first and second event in this dataset.
            Imputing a shadow time for all censored observations."
    )
  }

  if (method == "marginal") {
    d.imp2 <- purrr::map(d.imp1, function(x) {
      marginal_mi(x, bootstrap)
    })
  } else if (method == "cox") {
    d.imp2 <- purrr::map(d.imp1, function(x) {
      cox_mi(x, bootstrap)
    })
  }

  return(d.imp2)
}
