#Provide a wrapper function that seamlessly transitions between the first and second imputation layers
  #and allows the user to choose which method to use for the second layer


#' Create multiple imputed datasets for data arising from multi-state models subject to censoring
#'
#' @param d a dataframe with one row per subject and columns corresponding to the time and event indicators for each state transition
#' @param M an integer, the number of imputations
#' @param n.states an integer, the number of states in the multistate model
#' @param prefix.states a character vector of length 2, specify the prefix for the event and time columns in the d
#' @param method a character string, either "marginal" or "cox", indicating which method to use for the second layer of imputation
#'
#' @returns A list of length M, where each element is a data frame with imputed times for censored events
#' @export
msmi.impute <- function(d, M, n.states = 3, prefix.states = c("event", "t"), method = "marginal") {

  #prepare data in competing events structure for time to first event
  d.comp <- d %>% dplyr::mutate(t.first = pmin(t1, t2),
                                event.first = dplyr::case_when(t1 < t2 & event1 == 1 ~ 1,
                                                               t2 < t1 & event2 == 1 ~ 2,
                                                               TRUE ~ 0))
  #impute time to first event using mici::mici.impute
  d.imp1 <- mici::mici.impute(t.first, event.first, data=d.comp, scheme = "KMI", M = M)

  #put data back into the original format for second layer of imputation


}
