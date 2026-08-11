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


# Calculate Transition Matrix and Give Confidence Intervals ---------------

#' Calculate Empirical Transition Matrix from Multiple Imputed Datasets for an Illness-Death Model
#'
#' @param imp_obj Output from msmi.impute()
#' @param times Vector of times at which to calculate state occupation probabilities
#' @param method The interval estimation method to use: "trWald" or "goodman" (default is "trWald")
#' @param alpha Significance level for the confidence region (default is 0.05)
#'
#' @returns A list with the following components:
#'  \item{mi_estimate}{A data frame with the multiple imputation point estimates
#' of state occupation probabilities at each time}
#' \item{method}{The method used for confidence region construction}
#' \item{alpha}{The significance level for the confidence region}
#' \item{cr_list}{A list of confidence regions (defined via convex hulls) for the state occupation probabilities at each time
#'  in the probability space.}
#'
#' @examples
#' imps <- msmi.impute(dat = sim.data, M = 5,
#'             prefix.states = c("event", "t"), method = "marginal")
#' msmi.tprobs.v2(imp_obj = imps, times = seq(1, 5, 1))
#' @export

msmi.tprobs.v2 <- function(
  imp_obj = NULL,
  times = NULL,
  method = "trWald",
  alpha = 0.05
) {
  # Check the inputs
  if (!(method %in% c("trWald", "goodman"))) {
    stop(
      "method must be one of 'trWald' or 'goodman'"
    )
  }

  if (alpha <= 0 | alpha >= 1 | !is.numeric(alpha)) {
    stop("alpha must be a number between 0 and 1")
  }

  if (!is.numeric(times) | !is.vector(times)) {
    stop("times must be a numeric vector")
  }

  #Define constants
  n <- nrow(imp_obj[[1]]) #sample size
  M <- length(imp_obj) #number of imputations
  k <- 3 #number of states, fixed at 3

  ############################################################
  # Calculate point estimate (average over all imputations)  #
  ############################################################

  #For each imputed dataset, calculate the state occupation probabilities at each time of interest
  empirical_probs <- purrr::map_dfr(imp_obj, function(df) {
    get_empirical_probs(df, times)
  })

  # if (method == "agresticoull") {
  #   #for each imputation: add 4/k successes to each category and 4 trials total
  #   empirical_probs <- empirical_probs %>%
  #     dplyr::mutate(dplyr::across(
  #       dplyr::starts_with("p"),
  #       ~ (.x * n + (4 / 3)) / (n + 4)
  #     ))
  # }

  #multiple imputation point estimate
  mi_estimate <- empirical_probs %>%
    dplyr::group_by(time) %>%
    dplyr::summarise(
      dplyr::across(dplyr::starts_with("p"), \(x) mean(x, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    as.data.frame()

  rownames(mi_estimate) <- mi_estimate$time

  #############################################################
  # Construct confidence region (pooled over all imputations) #
  #############################################################

  if (method == "goodman") {
    pooled_goodman <- purrr::map(times, function(t) {
      #Matrix of estimates across imputations
      d <- empirical_probs %>%
        dplyr::filter(time == t) %>%
        dplyr::select(dplyr::starts_with("p")) %>%
        as.matrix()

      #Point estimate (average, same as mi_estiamte)
      mean_q <- colMeans(d)

      # Between and within impuation variance
      var.between <- apply(d, 2, var)
      var.within <- colMeans(d * (1 - d) / n)
      r.m <- (1 + (1 / M)) * (var.between / var.within)

      #Reference distribution
      v <- (M - 1) * (1 + (1 / r.m))^2

      #Lott and Reiter's Ad Hoc Solution for undefined r.m
      r.m[is.nan(r.m)] <- 0
      v[is.nan(v)] <- Inf

      #Critical value
      crit_value <- qt(1 - (alpha / (2 * k)), df = v)

      #Interval construction
      t_n <- (crit_value^2) / n
      tr_n <- (crit_value^2) * r.m / n
      q_sum <- 2 * mean_q + t_n + tr_n
      one_sum <- 2 * (1 + t_n + tr_n)

      centerpoint <- q_sum / one_sum
      summand <- sqrt((q_sum^2 / one_sum^2) - (2 * (mean_q^2) / one_sum))

      lower <- centerpoint - summand
      upper <- centerpoint + summand

      tibble(est = mean_q, lower = lower, upper = upper)
    })
    names(pooled_goodman) <- as.character(times)

    cr_list <- purrr::map(pooled_goodman, function(d) {
      rectangle_simplex(lower = d$lower, upper = d$upper)
    })
  } else if (method == "trWald") {
    #else if (method %in% c("agresticoull", "trWald"))
    cr_list <- purrr::map(times, function(t) {
      #Matrix of estimates across imputations
      p <- empirical_probs %>%
        dplyr::filter(time == t) %>%
        dplyr::select(dplyr::starts_with("p")) %>%
        as.matrix()

      # throw a warning if Wald method will have collapsed uncertainty at this time point
      if (any(colMeans(p) == 0)) {
        message(
          "At least one estimated state occupation probability is 0 at time ",
          as.character(t),
          ". Wald-type confidence regions are undefined.\nTry using method = 'goodman' instead."
        )
        return(NULL)
      }
      if (any(p == 0)) {
        message(
          "At least one imputation had a state occupation probability on the boundary (probability = 0) at time ",
          as.character(t),
          ". Wald-type confidence regions are undefined.\nTry using method = 'goodman' instead."
        )
        return(NULL)
      }

      #transform proportions on the 2-simplex to an unconstrained space
      x <- t(apply(p, 1, multinomial_logit))

      #mean of estimates on the unconstrained scale
      mean_x <- colMeans(x)

      #between variance: 1/(M-1)*sum[(x-x_est)'(x-x_est)] {Note: inner product because x_est is a row vector in R}
      dif <- sweep(x, 2, mean_x, FUN = "-")
      var.between <- (t(dif) %*% dif) / (M - 1)

      #vcov matrix on the unconstrained scale in terms of probabilities
      compute_within_var_mat <- function(row, n) {
        mat <- matrix(0, k - 1, k - 1)
        ref <- row[k]
        var_mat <- diag(1 / row[1:(k - 1)])
        cov_mat <- matrix(1 / ref, nrow = k - 1, ncol = k - 1)
        mat <- (var_mat + cov_mat) / n

        return(mat)
      }

      #Within variance
      var.within <- purrr::map(
        seq_len(M),
        ~ compute_within_var_mat(p[.x, ], n)
      ) %>%
        purrr::reduce(`+`) /
        M

      # #Rubin's Total variance
      # var.total <- var.within + (1 + 1 / M) * var.between

      # #Naive Pooled Wald
      # #Build eliptical region in the 2D space and transform points on the boundary of the region to the simplex
      # naiive <- t(apply(
      #   mixtools::ellipse(
      #     mu = mean_x,
      #     sigma = var.total,
      #     alpha = alpha,
      #     npoints = 500,
      #     draw = FALSE
      #   ),
      #   1,
      #   multinomial_logit_inverse
      # ))

      #D1 Method: Pooled Multivariate Wald Test

      #D1 Test Statistic
      r <- (1 + 1 / M) * sum(diag(var.between %*% solve(var.within))) / (k - 1)
      T_bar <- (1 + r) * var.within

      #Reference distribution
      t <- (k - 1) * (M - 1)
      nu <- 4 + (t - 4) * (1 + (1 - 2 * t^(-1)) * r^(-1))^2
      F_value <- qf(1 - alpha, df1 = k - 1, df2 = nu) * (k - 1)

      #Build ellipse based on the critical value
      angles <- seq(0, 2 * pi, length.out = 500)
      circle <- cbind(cos(angles), sin(angles))
      A <- chol(T_bar)

      ellipse <- sweep(
        sqrt(F_value) * circle %*% A,
        2,
        mean_x,
        "+"
      )

      D1 <- t(apply(
        ellipse,
        1,
        multinomial_logit_inverse
      ))

      D1
    })
  }
  # } else if (method == "score") {
  #   cr_list <- purrr::map(times, function(t) {
  #     #Matrix of estimates across imputations
  #     x <- empirical_probs %>%
  #       dplyr::filter(time == t) %>%
  #       dplyr::select(dplyr::starts_with("p")) %>%
  #       as.matrix()
  #     x <- x * n

  #     p_hat <- colMeans(x) / n

  #     score_test <- function(p) {
  #       #Calculate score test statistic for each imputed dataset
  #       U <- sweep(x[, 1:(k - 1), drop = FALSE], 2, p[1:(k - 1)], "/") -
  #         x[, k] / p[k]

  #       I <- (n / p[k]) * matrix(1, k - 1, k - 1)
  #       diag(I) <- diag(I) + n / p[1:(k - 1)]

  #       Iinv <- solve(I)

  #       statistics <- apply(U, 1, function(u) {
  #         drop(t(u) %*% Iinv %*% u)
  #       })

  #       q_bar <- mean(statistics)

  #       #Combine the test statistics using the D2 approach
  #       r_2 <- (1 + (1 / M)) * var(sqrt(statistics))
  #       num <- q_bar / (k - 1) - (M + 1) * r_2 / (M - 1)
  #       denom <- 1 + r_2
  #       test_stat <- num / denom

  #       #Compare pooled statistic to reference distribution
  #       nu <- ((k - 1)^(-3 / M)) * (M - 1) * (1 + 1 / r_2)^2
  #       p.value <- pf(test_stat, df1 = k - 1, df2 = nu, lower.tail = FALSE)

  #       # return(list(
  #       #   test_stat = test_stat,
  #       #   nu = nu,
  #       #   p.value = p.value
  #       # ))

  #       return(p.value)
  #     }

  #     #Adaptive Grid Search
  #     eps <- 1e-6
  #     grid <- seq(0.0001, 0.9999, by = 0.0015)

  #     ## Grid indexed by (i,j)
  #     G <- expand.grid(
  #       i = seq_along(grid),
  #       j = seq_along(grid)
  #     )

  #     G$p1 <- grid[G$i]
  #     G$p2 <- grid[G$j]
  #     G$p3 <- 1 - G$p1 - G$p2

  #     ## Keep only simplex
  #     G <- G[G$p3 > eps, ]

  #     ## Fast lookup from (i,j) -> row number
  #     key <- paste(G$i, G$j)
  #     lookup <- setNames(seq_len(nrow(G)), key)

  #     ## Status:
  #     ## NA = not evaluated
  #     ## TRUE = accepted
  #     ## FALSE = rejected
  #     status <- rep(NA, nrow(G))

  #     ## Start at grid point nearest phat
  #     start <- which.min(
  #       (G$p1 - p_hat[1])^2 +
  #         (G$p2 - p_hat[2])^2
  #     )

  #     queue <- start

  #     while (length(queue) > 0) {
  #       idx <- queue[1]
  #       queue <- queue[-1]

  #       ## already evaluated
  #       if (!is.na(status[idx])) {
  #         next
  #       }

  #       p <- c(G$p1[idx], G$p2[idx], G$p3[idx])

  #       status[idx] <- score_test(p) >= alpha

  #       ## Only expand from accepted points
  #       if (status[idx]) {
  #         i <- G$i[idx]
  #         j <- G$j[idx]

  #         nbrs <- rbind(
  #           c(i + 1, j),
  #           c(i - 1, j),
  #           c(i, j + 1),
  #           c(i, j - 1),
  #           c(i + 1, j - 1),
  #           c(i - 1, j + 1)
  #         )

  #         nbr_keys <- paste(nbrs[, 1], nbrs[, 2])

  #         nbr_idx <- lookup[nbr_keys]
  #         nbr_idx <- nbr_idx[!is.na(nbr_idx)]

  #         ## Add only unevaluated neighbors
  #         nbr_idx <- nbr_idx[is.na(status[nbr_idx])]

  #         queue <- c(queue, nbr_idx)
  #       }
  #     }

  #     accepted <- as.matrix(
  #       G[which(status), c("p1", "p2", "p3")]
  #     )

  #     hull <- chull(
  #       accepted[, 1],
  #       accepted[, 2]
  #     )

  #     boundary_hull <- accepted[hull, ]
  #   })
  # }

  names(cr_list) <- times

  # Create return object
  list(
    mi_estimate = mi_estimate,
    method = method,
    alpha = alpha,
    cr_list = cr_list
  )
}
