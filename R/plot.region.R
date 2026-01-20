#' Plot Confidence Regions State Occupation Probabilities
#'
#' @param p_draws A matrix of points in the confidence region (each row is a point, each column is a state) or a list with two elements:
#'                'inside' and 'outside', each a matrix of points inside and outside the confidence region respectively
#' @param p_hat (Optional) A data frame with one row and three columns (p1, p2, p3) representing the point estimate to overlay on the plot
#' @param time (Optional) The time point corresponding to the state occupation probabilities, used for plot title
#'
#' @returns A ggtern plot object showing the confidence region and optional point estimate
#' @export
#'
#' @examples
#' imps <- msmi.impute(dat = sim.data, M = 5, n.states = 3,
#'                     prefix.states = c("event", "t"), method = "marginal")
#' tprobs <- msmi.tprobs(imp_obj = imps, times = c(1,3))
#' #Plot the confidence region for state occupation probabilities at time = 1
#' plot_region( p_draws = tprobs$cr_list[['3']][['inside']], p_hat   = tprobs$mi_estimate['3', c("p1", "p2", "p3")], time = 3)

plot_region <- function(p_draws,
                        p_hat = NULL,
                        time = NULL) {

  # allow either a matrix (inside only) or a list(inside/outside)
  if (is.list(p_draws)) {
    inside  <- p_draws$inside
    outside <- p_draws$outside
  } else {
    inside  <- p_draws
    outside <- NULL
  }

  stopifnot(ncol(inside) == 3)

  df_in <- as.data.frame(inside)
  colnames(df_in) <- c("p1", "p2", "p3")
  df_in$region <- "Inside"

  if (!is.null(outside)) {
    stopifnot(ncol(outside) == 3)
    df_out <- as.data.frame(outside)
    colnames(df_out) <- c("p1", "p2", "p3")
    df_out$region <- "Outside"

    df <- rbind(df_in, df_out)
  } else {
    df <- df_in
  }

  if (!is.null(time)) {
    title <- paste0(
      "95% Confidence Region for State Occupation Probabilities\n",
      "at time t = ", as.character(time)
    )
  } else {
    title <- "95% Confidence Region for State Occupation Probabilities"
  }

  p <- ggtern::ggtern(
    df,
    ggplot2::aes(x = p1, y = p2, z = p3, color = region)
  ) +
    ggplot2::geom_point(size = 0.6, alpha = 0.35) +
    ggtern::theme_rgbw() +
    ggtern::Tlab("Healthy") +
    ggtern::Llab("Ill") +
    ggtern::Rlab("Dead") +
    ggplot2::scale_color_manual(
      values = c("Inside" = "blue", "Outside" = "red")
    ) +
    ggplot2::ggtitle(title)

  if (!is.null(p_hat)) {
    p <- p +
      ggplot2::geom_point(
        data = p_hat,
        ggplot2::aes(x = p1, y = p2, z = p3),
        color = "black",
        size = 1.2,
        inherit.aes = FALSE
      )
  }

  return(p)
}
#TO DO: Fix the warning message
  # "In geom_point(data = p_hat, aes(p1, p2, p3), color = "red", size = 0.5,  :
  # Ignoring unknown aesthetics: z"
#TO DO: Fix the scale of the axes to go from 0 to 1
