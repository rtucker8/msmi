# Empirical Probability in State ------------------------------------------


#' Calculate empirical state occupation probabilities for a single imputed dataset
#'
#' @param df A data frame with imputed event times and event indicators
#' @param times Vector of times at which to calculate state occupation probabilities
#'
#' @returns A data frame with empirical state occupation probabilities at specified times
get_empirical_probs <- function(df, times) {

  # Initialize an empty list to store results for each time
  results <- lapply(times, function(t) {
    # logical conditions for each state
    is_healthy <- df$t1 > t & df$t2 > t
    is_ill <- df$t1 <= t & df$t2 > t
    is_dead <- df$t2 <= t


    # Calculate proportions in each state
    n <- nrow(df)
    p1 <- sum(is_healthy) / n
    p2<- sum(is_ill) / n
    p3 <- sum(is_dead) / n

    tibble(
      time = t,
      p1 = p1,
      p2= p2,
      p3= p3
    )
  })

  bind_rows(results)
}


# Calculate Transition Matrix and Give Confidence Intervals ---------------

#' Calculate Empirical Transition Matrix from Multiple Imputed Datasets for an Illness-Death Model
#'
#' @param imp_obj Output from msmi.impute()
#' @param times Vector of times at which to calculate state occupation probabilities
#' @param int.type Type of confidence region to compute: "multinomial_logistic" or "dirichlet" (default is "multinomial_logistic")
#' @param alpha Significance level for the confidence region (default is 0.95)
#' @param vcov Logical indicating whether to return the variance-covariance matrices for the state occupation probabilities at each time (default is FALSE)
#'
#' @returns A list with the following components:
#'  \item{mi_estimate}{A data frame with the multiple imputation point estimates
#' of state occupation probabilities at each time}
#' \item{int.type}{The type of confidence interval used}
#' \item{alpha}{The significance level for the confidence region}
#' \item{vcov}{if vcov = TRUE, A list of variance-covariance matrices for the state occupation probabilities at each time in the unconstrained space}
#' \item{cr_list}{A list of confidence regions (defined via point clouds) for the state occupation probabilities at each time
#'  in the (1) unconstrained 2D space and (2) the probability space.}
#'
#' @export
#'
#' @examples
#' imps <- msmi.impute(dat = sim.data, M = 5, n.states = 3,
#'             prefix.states = c("event", "t"), method = "marginal")
#' msmi.tprobs(imp_obj = imps, times = seq(1, 5, 1))
msmi.tprobs <- function(imp_obj = NULL,
                        times = NULL,
                        int.type = "multinomial_logistic",
                        alpha = 0.95,
                        vcov = FALSE) {

  #For each imputed dataset, calculate the state occupation probabilities at each time of interest
  empirical_probs <- purrr::map(imp_obj, function(df) {
    get_empirical_probs(df, times)
  })

  #Rubin's Rules Point Estimate: Average empirical probabilities across imputations for each simulated dataset
  mi_estimate <- dplyr::bind_rows(empirical_probs, .id = "imputation") %>%
      dplyr::group_by(time) %>%
      dplyr::summarise(dplyr::across(dplyr::starts_with("p"), \(x) mean(x, na.rm = TRUE)), .groups = "drop") %>%
    as.data.frame()

  rownames(mi_estimate) <- mi_estimate$time

  #Calculate Confidence Intervals


  #transform probability simplex to unconstrained space with dimension k-1 at each time point
  #TO NOTE: transformation is undefined on the boundary of p, i.e. p=0 or p=1, so we add epsilon at the boundary to make it numerically stable
  #TO NOTE: var(theta1) = INF at latter timepoints when p1 = 0 which causes problems; consequence of WALD interval
  if (int.type == "dirichlet") {

  } else if (int.type == "multinomial_logistic") {

    #constants for Rubin's Rules
    n <- nrow(imp_obj[[1]]) #sample size
    M <- length(imp_obj) #number of imputations
    k <- 3  #number of states, fixed at 3 for now

    #create Mxk matrix of estimates for each time point
    ps <- dplyr::bind_rows(empirical_probs) %>%
      dplyr::group_by(time) %>%
      dplyr::group_split()

    p_list <- purrr::map(ps, function(x) {
      as.matrix(x[ , names(x) != "time"])
    })

    x_list <- purrr::map(p_list, function(p) {
      t(apply(p, 1, multinomial_logit))
    })

  #calculate mean estimate on unconstrained scale at each time point
  x_est <- purrr::map(x_list, function(x) {
    colMeans(x)
  })
  names(x_est) <- times

  #put unconstrained estimates into a dataframe
  unconstrained_estimate <- data.frame(
    time = as.numeric(names(x_est)),
    do.call(rbind, x_est)
  )
  names(unconstrained_estimate)[2:3] <- c("theta1", "theta2")
  rownames(unconstrained_estimate) <- unconstrained_estimate$time

  #between variance: 1/(M-1)*sum[(x-x_est)'(x-x_est)] {Note: inner product because x_est is a row vector in R}
  dif_list <- purrr::map2(x_list, x_est, function(x, est) {
    sweep(x, 2, est, FUN = "-")
  })

  between_var_list <- purrr::map(dif_list, function(d) {
    (t(d) %*% d) / (M - 1)
  })

  compute_within_var_mat <- function(row, n) {
    mat = matrix(0, k - 1, k - 1)
    ref = row[k]

    var_mat = diag(1/row[1:(k-1)])
    cov_mat = matrix(1/ref, nrow = k-1, ncol = k-1)
    mat = (var_mat + cov_mat)/n

    return(mat)
  }

  within_var_list <- purrr::map(p_list, function(p) {
    purrr::map(seq_len(M), ~ compute_within_var_mat(p[.x, ], n)) %>%
      purrr::reduce(`+`) / M
  })

  #total variance
  total_var_list <- purrr::map2(between_var_list, within_var_list, function(b, w) {
    w + (1 + 1/M)*b
  })

  #compute joint confidence regions at each time point
  conf.width <- 1 - alpha
  crit <- stats::qchisq(conf.width, df = k - 1, lower.tail = FALSE)

  #function to get a point cloud that lies within the Wald region
  wald_region <- function(theta_hat, Sigma, n_draw = 5000) {

    draws <- MASS::mvrnorm(n_draw, mu = theta_hat, Sigma = Sigma)

    quad <- apply(draws, 1, function(z) {
      d <- z - theta_hat
      emulator::quad.form.inv(Sigma, d)
    })

    return(draws[quad <= crit, , drop = FALSE])
  }

  cr_list <- purrr::map2(x_est, total_var_list, function(theta_hat, Sigma) {
    theta_draws <- wald_region(theta_hat, Sigma)
    return(list(unconstrained = theta_draws, p.space = t(apply(theta_draws,  1, multinomial_logit_inverse))))
  })

  names(cr_list) <- times
  names(total_var_list) <- times

  } else {
    stop("int.type must be one of 'multinomial_logistic' or 'dirichlet'")
  }

  if (vcov == FALSE) {
    return(list(mi_estimate = mi_estimate, unconstrained_estimate = unconstrained_estimate, int.type = int.type, alpha = alpha, cr_list = cr_list))
  } else if (vcov == TRUE) {
    return(list(mi_estimate = mi_estimate, unconstrained_estimate = unconstrained_estimate, int.type = int.type, alpha = alpha, vcov = total_var_list, cr_list = cr_list))
  }

}
