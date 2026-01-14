#' Plot Confidence Regions State Occupation Probabilities
#'
#' @param p_draws A matrix of points in the confidence region (each row is a point, each column is a state)
#' @param p_hat (Optional) A data frame with one row and three columns (p1, p2, p3) representing the point estimate to overlay on the plot
#'
#' @returns A ggtern plot object showing the confidence region and optional point estimate
#' @export
#'
#' @examples
#' imps <- msmi.impute(dat = sim.data, M = 5, n.states = 3,
#'                     prefix.states = c("event", "t"), method = "marginal")
#' tprobs <- msmi.tprobs(imp_obj = imps, times = c(1,2))
#' #Plot the confidence region for state occupation probabilities at time = 1
#' plot_region( p_draws = tprobs$cr_list[[1]], p_hat   = tprobs$mi_estimate[1, c("p1", "p2", "p3")])

plot_region <- function(p_draws,
                        p_hat = NULL) {

  stopifnot(ncol(p_draws) == 3)

  df <- as.data.frame(p_draws)
  colnames(df) <- c("p1", "p2", "p3")

  p <- ggtern::ggtern(df, ggplot2::aes(x=p1, y=p2, z=p3)) +
    ggplot2::geom_point(size = 0.6, alpha = 0.3) +
    ggtern::theme_rgbw() +
    ggtern::Tlab("Healthy") +
    ggtern::Llab("Ill") +
    ggtern::Rlab("Dead") +
    ggplot2::ggtitle("95% Confidence Region for State Occupation Probabilies")

  if (!is.null(p_hat)) {

    p <- p +
      ggplot2::geom_point(
        data = p_hat,
        ggplot2::aes(x=p1, y=p2, z=p3),
        color = "red",
        size = 0.5,
        inherit.aes = FALSE
      )
  }

  return(p)
}

#TO DO: Fix the warning message
  # "In geom_point(data = p_hat, aes(p1, p2, p3), color = "red", size = 0.5,  :
  # Ignoring unknown aesthetics: z"
#TO DO: Fix the scale of the axes to go from 0 to 1
