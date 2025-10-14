
# Cox Imputation ----------------------------------------------------------


#' Layer 2 Imputation: Cox Model (Conditional) Approach
#'
#' @param d A data frame with one row per subject and the columns event1, t1, event2, t2 on which mici::mici.impute was previously ran
#' @param m An integer, the nest id for the nested imputation structure
#' @param r An integer, the impute id for the nested imputation structure
#'
#' @returns A data frame with imputed times for event2 where event2 was censored
nested_cox_mi <- function(d, m, r) {

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
  surv_times_list <- replicate(length(surv_probs), surv_times, simplify = FALSE)

  prob_diffs <- apply(surv_summary$surv, 2, function(x) -diff(c(1, x)))
  prob_diffs_list <- split(prob_diffs, col(prob_diffs))

  # Handle tail probability (if survival doesn't reach 0)
  tail.probs <- function(tail, probs) {
    if (tail > 0) {
      probs <- c(probs, tail)
    }
    return(probs)
  }

  tail.times <- function(tail, times) {
    if (tail > 0) {
      times <- c(times, max(u$sojourn12) + 1)
    }
    return(times)
  }

  prob_diffs <- Map(tail.probs, tail = surv_probs, probs = prob_diffs_list) #output is a list of length nrow(dd)

  surv_times <- Map(tail.times, surv_probs, surv_times_list) #output is a list of length nrow(dd)

  # Impute times for censored individuals
  cts <- NULL
  for (jj in 1:length(xt)) {
    # Find times greater than censoring time
    sub <- surv_times[[jj]] > xt[jj]
    # Sample time from illness to death
    if (sum(sub) > 1) {
      cts[jj] <- resample(surv_times[[jj]][sub], 1,replace=TRUE, prob = prob_diffs[[jj]][sub])
    } else if(sum(sub) == 1){
      cts[jj] <- surv_times[[jj]][sub]
    } else {
      cts[jj] <- max(u$sojourn12) + 1 #shadow event time
    }
  }
  # Update data
  dd$event2 <- 1
  dd$t2 <- dd$t1 + cts
  ipd <- dplyr::bind_rows(dd, dc) %>% dplyr::arrange(id) %>% dplyr::select(-id)
  ipd$nest.id <- m
  ipd$impute.id <- r

  return(ipd)
}


# Marginal Imputation -----------------------------------------------------
#Idea use clock forward timescale instead of reset?

#' Layer 2 Imputation: Marginal Approach
#'
#' @param d A data frame with one row per subject and the columns event1, t1, event2, t2 on which mici::mici.impute was previously ran
#' @param m An integer, the nest id for the nested imputation structure
#' @param r An integer, the impute id for the nested imputation structure
#'
#' @returns A data frame with imputed times for event2 where event2 was censored
nested_marginal_mi <- function(d, m, r) {

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
  ipd$nest.id <- m
  ipd$impute.id <- r

  return(ipd)
}


# Wrapper Function --------------------------------------------------------

#Provide a wrapper function that seamlessly transitions between the first and second imputation layers
#and allows the user to choose which method to use for the second layer

#' Create multiple imputed datasets for data arising from multi-state models subject to censoring
#'
#' @param dat a dataframe with one row per subject and columns corresponding to the time and event indicators for each state transition.
#'  The column names for the time and event indicators should be in the format specified by prefix.states, specifically
#'  event<i> and t<i> for i = 1, ..., n.states-1
#' @param M an integer, the number of imputation nests (1st layer)
#' @param R an integer, the number of imputations within each nest (2n layer)
#' @param n.states an integer, the number of states in the multistate model
#' @param prefix.states a character vector of length 2, specify the prefix for the event and time columns in the d in that order
#' @param method a character string, either "marginal" or "cox", indicating which method to use for the second layer of imputation
#' @param seed an integer, the seed for random number generation
#' @examples
#' msmi.nested.impute(sim.data, M = 5, R=3, n.states = 3, prefix.states = c("event", "t"), method = "marginal")
#' @returns A nested list with MxR total entries, where each element is a data frame with imputed times for censored events
#' @export
msmi.nested.impute <- function(dat, M, R = 1, n.states = 3, prefix.states = c("event", "t"), method = "marginal", seed = sample(1:.Machine$integer.max, size=1)) {

  #Check inputs
  if (n.states != 3) {
    stop("Currently msmi only supports 3-state models")
  }
  if (length(prefix.states) != 2) {
    stop("prefix.states must be a character vector of length 2")
  }
  if (method != "marginal" & method != "cox") {
    stop("method must be either 'marginal' or 'cox'")
  }
  if (!(M == floor(M)) | !(M > 0)) {
    stop("M must be a positive integer")
  }
  if (!(M == floor(M)) | !(R > 0)) {
    stop("M must be a positive integer.")
  }


  #set seed
  set.seed(seed)

  #Create standardized dataset d from user provided dataframe dat
  d <- data.frame(row.names = 1:nrow(dat))

  for (i in 1:(n.states-1)) {
    if (!all(c(paste0(prefix.states[1], i), paste0(prefix.states[2], i)) %in% colnames(dat))) {
      stop(paste0("Columns ", paste0(prefix.states[1], i), " and ", paste0(prefix.states[2], i), " must be present in the data"))
    }
    d[[paste0("event", i)]]  <- dat[, paste0(prefix.states[1], i)]
    d[[paste0("t", i)]]  <- dat[, paste0(prefix.states[2], i)]

  }
  #d now has columns eventi, ti for i = 1, 2, ..., n.states-1


  #prepare data in competing events structure for time to first event
  d.comp <- d %>% dplyr::mutate(t.first = pmin(t1, t2),
                                event.first = dplyr::case_when(t1 < t2 & event1 == 1 ~ 1,
                                                               !(t1 < t2) & event2 == 1 ~ 2,
                                                               TRUE ~ 0))
  #create M multiple imputations for time to first event using mici.impute
  d.imp1 <- mici.impute(t.first, event.first, data=d.comp, scheme = "KMI", M = M)

  #put data back into the original format for second layer of imputation
  d.imp1 <- purrr::map(d.imp1, function(x) {

    idx2 <- x$ftype == 2 & x$event.first == 0 #indices for originally censored people who were imputed to have event type 2 first
    idx1 <- x$ftype == 1 & x$event.first == 0 #indices for originally censored people who were imputed to have event type 1 first

    x$event1[idx2] <- 0  #those imputed to have event type 2 first did not have event 1
    x$event1[idx1] <- 1  #those imputed to have event type 1 first did have event 1
    x$t1[idx2 | idx1] <- x$ftime[idx2 | idx1] #everyone who was censored and has their t1 updated to the imputed time
    x$event2[idx2] <- 1  #those imputed to have event type 2 first did have event 2
    x$event2[idx1] <- 0  #those imputed to have event type 1 first did not have event 2
    x$t2[idx2 | idx1] <- x$ftime[idx2 | idx1] #everyone who was censored has their t2 updated to the imputed time

    x %>% dplyr::select(t1, event1, t2, event2)
  })

  #for each imputation in first layer, create R nested imputations for for time to second event
  if (method == "marginal") {
    d.imp2 <- purrr::imap(d.imp1, function(d, m) {
      purrr::map(1:R, function(r) {
        nested_marginal_mi(d, m, r)
      })
    })
  }

  if (method == "cox") {
    d.imp2 <- purrr::imap(d.imp1, function(d, m) {
      purrr::map(1:R, function(r) {
        nested_marginal_cox(d, m, r)
      })
    })
  }

  return(d.imp2)
}



