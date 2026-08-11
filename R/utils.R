#Clear notes in check() log
id <- t1 <- t2 <- t.first <- event.first <- event1 <- event2 <- p1 <- p2 <- p3 <- time <- theta1 <- theta2 <- a1 <- a2 <- a3 <- imp <- NULL

#wrapper for the sample function that will take a sample of size 1 if x has length 1
resample <- function(x, ...) x[sample.int(length(x), ...)]

#Multinomial logit transformation
multinomial_logit_inverse <- function(x) {
  k <- length(x) + 1
  p <- numeric(k)

  denom <- 1 + sum(exp(x))

  for (i in 1:(k - 1)) {
    p[i] <- exp(x[i]) / denom
  }
  p[k] <- 1 / denom

  return(p)
}

#Code to safely compute sum(exp(x)) without overflow
# #log_sum_exp <- function(x) {
# max_x <- max(x)
# return(max_x + log(sum(exp(x - max_x))))
# }
#
# # To get the actual sum safely:
# safe_sum_exp <- function(x) {
#   exp(log_sum_exp(x))
# }

multinomial_logit <- function(p) {
  #add eps for numerical stability and renormalize
  eps <- 1e-8
  p[p == 0] <- eps
  p <- p / sum(p)

  k <- length(p)
  x <- numeric(k - 1)

  for (i in 1:(k - 1)) {
    x[i] <- log(p[i] / p[k])
  }

  return(x)
}

#find polygon that intersects the hyperrectangle on the simplex
rectangle_simplex <- function(lower, upper, tol = 1e-10) {
  lines <- list(
    c(1, 0, lower[1]), # p1 = lower1
    c(1, 0, upper[1]), # p1 = upper1
    c(0, 1, lower[2]), # p2 = lower2
    c(0, 1, upper[2]), # p2 = upper2
    c(1, 1, 1 - upper[3]), # p1+p2 = 1-upper3
    c(1, 1, 1 - lower[3]) # p1+p2 = 1-lower3
  )

  pts <- list()

  for (i in 1:5) {
    for (j in (i + 1):6) {
      A <- rbind(lines[[i]][1:2], lines[[j]][1:2])

      if (abs(det(A)) < tol) {
        next
      }

      b <- c(lines[[i]][3], lines[[j]][3])

      x <- solve(A, b)

      p <- c(x, 1 - sum(x))

      if (
        all(p >= lower - tol) &&
          all(p <= upper + tol)
      ) {
        pts[[length(pts) + 1]] <- p
      }
    }
  }

  pts <- unique(do.call(rbind, lapply(pts, round, 10)))

  ## order vertices
  cen <- colMeans(pts[, 1:2])

  ang <- atan2(pts[, 2] - cen[2], pts[, 1] - cen[1])

  pts[order(ang), ]
}
