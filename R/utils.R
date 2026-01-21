#Clear notes in check() log
id <- t1 <- t2 <- t.first <- event.first <- event1 <- event2 <- p1 <- p2 <- p3 <- time <- theta1 <- theta2 <- NULL

#wrapper for the sample function that will take a sample of size 1 if x has length 1
resample <- function(x, ...) x[sample.int(length(x), ...)]

#Multinomial logit transformation
multinomial_logit_inverse <- function(x) {
  k <- length(x) + 1
  p <- numeric(k)

  denom <- 1 + sum(exp(x))

  for (i in 1:(k-1)) {
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
  p = p / sum(p)

  k <- length(p)
  x <- numeric(k-1)

  for (i in 1:(k-1)) {
    x[i] <- log(p[i] / p[k])
  }

  return(x)
}
