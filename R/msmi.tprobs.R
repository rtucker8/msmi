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
    p2 <- sum(is_ill) / n
    p3 <- sum(is_dead) / n

    tibble(
      time = t,
      p1 = p1,
      p2 = p2,
      p3 = p3
    )
  })

  bind_rows(results)
}

#' Calculate Multinomial Distribution Log-Likelihood
#'
#' @param p a k-dimensional vector characterizing the cell probabilities for a multinomial distribution
#' @param counts a k-dimensional vector giving the observed counts for each category arising from a multinomial distribution
#'
#' @returns the multinomial log-likelihood, scalar
multinom_loglik <- function(p, counts) {

  if(any(p < 0))
    return(-Inf)

  if(abs(sum(p) - 1) > 1e-10)
    return(-Inf)

  # impossible point
  if(any(counts > 0 & p == 0))
    return(-Inf)

  idx <- counts > 0

  sum(counts[idx] * log(p[idx]))
}

#' Calculate Multinomial Deviance
#'
#' @param p a k-dimensional vector characterizing the cell probabilities for a multinomial distribution
#' @param counts a k-dimensional vector giving the observed counts for each category arising from a multinomial distribution
#'
#' @returns the multinomial deviance for a candidate probability vector compared to the MLE, scalar
multinom_deviance <- function(p, counts) {

  phat <- counts / sum(counts)

  2 * (multinom_loglik(phat, counts) - multinom_loglik(p, counts))
}

simplex_grid <- function(n) {
  g <- expand.grid(i = 0:n, j = 0:n)
  g <- subset(g, i + j <= n)

  transform(
    g,
    x1 = i / n,
    x2 = j / n,
    x3 = (n - i - j) / n
  )[c("x1", "x2", "x3")]
}


#' Calculate Meng-Rubin LRT Statistic
#'
#' @param L_counts a list of k-dimensional vectors given the observed counts for each state under each imputation at a desired point in time
#' @returns
meng.rubin <- function(L_counts) {

  #candidate values for confidence region
  df_p <- as.matrix(simplex_grid(1000))

  #list of deviance values for candidate points under each imputation's estimate
  s.e <- map(L_counts, function(counts) {

    d <- apply(df_p, 1, function(p) multinom_deviance(p, counts))
    dplyr::bind_cols(df_p, tibble::tibble(d = d))
  })

  #average deviance across imputations for each candidate point
  Dbar <- dplyr::bind_cols(df_p, tibble::tibble(Dbar =  rowMeans(sapply(s.e, function(df) df$d))))

}

# Calculate Transition Matrix and Give Confidence Intervals ---------------

#' Calculate Empirical Transition Matrix from Multiple Imputed Datasets for an Illness-Death Model
#'
#' @param imp_obj Output from msmi.impute()
#' @param times Vector of times at which to calculate state occupation probabilities
#' @param int.type Type of confidence region to compute: "trWald" or "bayes" (default is "trWald")
#' @param prior List of prior parameters for the Dirichlet distribution when int.type = "bayes" (default = NULL)
#' @param alpha Significance level for the confidence region (default is 0.05)
#' @param vcov Logical indicating whether to return the variance-covariance matrices for the state occupation probabilities at each time (default is FALSE)
#'
#' @returns A list with the following components:
#'  \item{mi_estimate}{A data frame with the multiple imputation point estimates
#' of state occupation probabilities at each time}
#' \item{unconstrained_estimate}{(if int.type = "trWald") A data frame with the multiple imputation point estimates}
#' \item{int.type}{The type of confidence interval used}
#' \item{alpha}{The significance level for the confidence region}
#' \item{vcov}{(if vcov = TRUE and int.type = "trWald") A list of variance-covariance matrices for the state occupation probabilities at each time in the unconstrained space}
#' \item{cr_list}{A list of confidence regions (defined via convex hulls) for the state occupation probabilities at each time
#'  in the (1) (if int.type = "trWald") unconstrained 2D space and (2) the probability space.}
#'
#' @export
#'
#' @examples
#' imps <- msmi.impute(dat = sim.data, M = 5,
#'             prefix.states = c("event", "t"), method = "marginal")
#' msmi.tprobs(imp_obj = imps, times = seq(1, 5, 1))
msmi.tprobs <- function(
  imp_obj = NULL,
  times = NULL,
  int.type = "trWald",
  prior = NULL,
  alpha = 0.05,
  vcov = FALSE
) {
  #For each imputed dataset, calculate the state occupation probabilities at each time of interest
  empirical_probs <- purrr::map(imp_obj, function(df) {
    get_empirical_probs(df, times)
  })

  #Rubin's Rules Point Estimate: Average empirical probabilities across imputations for each simulated dataset
  mi_estimate <- dplyr::bind_rows(empirical_probs, .id = "imputation") %>%
    dplyr::group_by(time) %>%
    dplyr::summarise(
      dplyr::across(dplyr::starts_with("p"), \(x) mean(x, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    as.data.frame()

  rownames(mi_estimate) <- mi_estimate$time

  #Calculate Confidence Intervals


  if (int.type == "bayes") {
    #interpretation, adding alpha_i successes in each category at timepoint of interest

    if (is.null(prior)) {
      stop("When int.type = 'bayes', prior must be specified.")
    }

    #calculate posterior parameters for Dirichlet distribution at each time point and for each imputation
    counts <- dplyr::bind_rows(empirical_probs, .id = "imp") %>%
      dplyr::mutate(dplyr::across(
        dplyr::starts_with("p"),
        ~ .x * nrow(imp_obj[[1]])
      ))
    posterior <- counts %>%
      dplyr::mutate(
        a1 = p1 + prior[1],
        a2 = p2 + prior[2],
        a3 = p3 + prior[3]
      ) %>%
      dplyr::select(imp, time, a1, a2, a3) %>%
      dplyr::group_by(time) %>%
      dplyr::group_split()

    n_draws <- 1000

    cr_list <- purrr::map(posterior, function(df) {
      # number of imputations
      M <- nrow(df)

      # store alpha vectors as a list
      alpha_list <- lapply(seq_len(M), function(i) {
        as.numeric(df[i, c("a1", "a2", "a3")])
      })

      #draw from each posterior and pool
      p_draws <- purrr::map_dfr(alpha_list, function(a) {
        draws <- gtools::rdirichlet(n_draws, a) %>% as.data.frame()
        return(draws)
      })

      #mixture density at each draw
      dens_mat <- sapply(alpha_list, function(a) {
        gtools::ddirichlet(p_draws, a)
      })
      mix_density <- rowMeans(dens_mat)
      p_mix <- cbind(p_draws, density = mix_density)

      #HPD cutoff
      cutoff <- stats::quantile(p_mix$density, probs = alpha)

      #return convex hull for HPD region for this time (could not be ideal since HPD region not necessarily contiguous)
      pts <- p_mix[p_mix$density >= cutoff, , drop = FALSE]
      pts <- pts[, colnames(pts) != "density"]

      chull <- grDevices::chull(pts)
      p_hull <- pts[chull, ]

      return(list(p.space = pts[chull, ]))
    })

    names(cr_list) <- times

    return(list(
      mi_estimate = mi_estimate,
      int.type = int.type,
      alpha = alpha,
      cr_list = cr_list
    ))

  } else if (int.type %in% c("agresticoull", "trWald")) {

    n = nrow(imp_obj[[1]]) #sample size
    M <- length(imp_obj) #number of imputations
    k <- 3 #number of states, fixed at 3 for now

    if (int.type == "agresticoull") {
      #add 4/3 pseudo-counts to each category at each time point and for each imputation
      pseudo <- dplyr::bind_rows(empirical_probs, .id = "imp") %>%
        dplyr::mutate(dplyr::across(dplyr::starts_with("p"), ~ (.x*n +(4/3))/(n+4))) %>%
        dplyr::rename(p1_adj = p1,
                      p2_adj = p2,
                      p3_adj = p3) %>%
        dplyr::group_by(time) %>%
        dplyr::group_split()

      #adjusted estimates at each time for each imputation
      p_list <- purrr::map(pseudo, function(x) {
        m <- as.matrix(x[, c("p1_adj", "p2_adj", "p3_adj")])
        return(m)
      })
    } else if (int.type == "trWald") {

      ps <- dplyr::bind_rows(empirical_probs) %>%
        dplyr::group_by(time) %>%
        dplyr::group_split()

      p_list <- purrr::map(ps, function(x) {
        m <- as.matrix(x[, names(x) != "time"])
        if (any(colMeans(m) == 0)) {
          # throw a warning because Wald method will have collapsed uncertainty at this time point
          message(
            "At least one estimated state occupation probability is 0 at time ",
            unique(x$time),
            ". Wald-type confidence regions are undefined. Try using int.type = 'bayes' instead."
          )
          return(NULL)
        }
        return(m)
      })
    }

    #transform proportions on the 2-simplex to an unconstrained space
    x_list <- purrr::map(p_list, function(p) {
      if (is.null(p)) {
        return(NULL)
      }
      t(apply(p, 1, multinomial_logit))
    })

    #calculate mean estimate on unconstrained scale at each time point
    x_est <- purrr::map(x_list, function(x) {
      if (is.null(x)) {
        return(NULL)
      }
      colMeans(x)
    })
    names(x_est) <- times

    #put unconstrained estimates into a dataframe
    unconstrained_estimate <- purrr::imap_dfr(x_est, function(x, t) {
      if (is.null(x)) {
        return(NULL)
      }
      data.frame(time = as.numeric(t), theta1 = x[1], theta2 = x[2])
    })

    rownames(unconstrained_estimate) <- unconstrained_estimate$time

    #between variance: 1/(M-1)*sum[(x-x_est)'(x-x_est)] {Note: inner product because x_est is a row vector in R}
    dif_list <- purrr::map2(x_list, x_est, function(x, est) {
      if (is.null(x) || is.null(est)) {
        return(NULL)
      }
      sweep(x, 2, est, FUN = "-")
    })

    between_var_list <- purrr::map(dif_list, function(d) {
      if (is.null(d)) {
        return(NULL)
      }
      (t(d) %*% d) / (M - 1)
    })

    compute_within_var_mat <- function(row, n) {
      mat = matrix(0, k - 1, k - 1)
      ref = row[k]

      var_mat = diag(1 / row[1:(k - 1)])
      cov_mat = matrix(1 / ref, nrow = k - 1, ncol = k - 1)
      mat = (var_mat + cov_mat) / n

      return(mat)
    }

    within_var_list <- purrr::map(p_list, function(p) {
      if (is.null(p)) {
        return(NULL)
      }
      purrr::map(seq_len(M), ~ compute_within_var_mat(p[.x, ], n)) %>%
        purrr::reduce(`+`) /
        M
    })

    #total variance
    total_var_list <- purrr::map2(
      between_var_list,
      within_var_list,
      function(b, w) {
        if (is.null(b) || is.null(w)) {
          return(NULL)
        }
        w + (1 + 1 / M) * b
      }
    )

    cr_list <- purrr::map2(x_est, total_var_list, function(theta_hat, Sigma) {
      if (any(!is.finite(Sigma))) {
        return(NULL)
      }
      if (is.null(theta_hat) || is.null(Sigma)) {
        return(NULL)
      }
      thetas <- mixtools::ellipse(
        mu = theta_hat,
        sigma = Sigma,
        alpha = alpha,
        npoints = 500,
        draw = FALSE
      )
      ps <- t(apply(thetas, 1, multinomial_logit_inverse))

      return(list(unconstrained = thetas, p.space = ps))
    })

    names(cr_list) <- times
    names(total_var_list) <- times

    #create return
    #TO DO: return adjusted estimates for Agresti-Coull approach
    out <- list(
      mi_estimate = mi_estimate,
      unconstrained_estimate = unconstrained_estimate,
      int.type = int.type,
      alpha = alpha,
      cr_list = cr_list
    )

    if (vcov) {
      out$vcov <- total_var_list
      out$w.in <- within_var_list
    }

    return(out)

  } else {
    stop("int.type must be one of 'trWald', 'agresticoull', or 'bayes'")
  }
}
