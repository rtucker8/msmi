#' World Health Organization TB data
#'
#' A subset of data from the World Health Organization Global Tuberculosis
#' Report ...
#'
#' @format ## `sim.data`
#' A data frame with 100 rows and 8 columns representing simulated data from an illness-death model (Healthy(0), Ill(1), Dead(3)) subject to censoring:
#' \describe{
#'   \item{id}{Unique person identifier}
#'   \item{t1, t2}{Follow up times for event type 1 and event type 2 respectively}
#'   \item{event1, event2}{Event indicators for event type 1 and event type 2 respectively (1 = event, 0 = censored)}
#'   \item{sojourn01}{Survival time for transition from state 0 to state 1 (healthy to ill)}
#'   \item{sojourn02}{Survival time for transition from state 0 to state 2 (healthy to death)}
#'   \item{sojourn12}{Survival time for transition from state 1 to state 2 (ill to death)}
#'   ...
#' }
"sim.data"
