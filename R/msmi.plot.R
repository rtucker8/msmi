#' Plot Confidence Regions State Occupation Probabilities
#'
#' @param tprob_obj The list returned by msmi.tprobs()
#' @param time  The time point of interest at which to generate the plot
#'
#' @returns A plot showing the confidence region and point estimate in the probability simplex.

#' @export
#'
#' @examples
#' imps <- msmi.impute(dat = sim.data, M = 5,
#'                     prefix.states = c("event", "t"), method = "marginal")
#' tprobs <- msmi.tprobs(imp_obj = imps, times = c(1,3))
#' #Plot the confidence region for state occupation probabilities at time = 1
#' msmi.plot(tprob_obj = tprobs, time = 1)

msmi.plot <- function(tprob_obj, time) {
  #Convex Hull for confidence region
  p_hull <- tprob_obj$cr_list[[as.character(time)]]

  if (is.null(p_hull)) {
    stop(
      "No confidence region available for this time point. Wald regions have collapsed uncertainty at the boundary."
    )
  }

  p_hull <- as.data.frame(p_hull)
  colnames(p_hull) <- c("p1", "p2", "p3")

  #Point estimate to overlay
  p_hat <- tprob_obj$mi_estimate[as.character(time), c("p1", "p2", "p3")]

  title <- paste0(
    as.character((1 - tprob_obj$alpha) * 100),
    "% Confidence Region for State Occupation Probabilities at time t = ",
    as.character(time)
  )

  #Make plot
  suppressWarnings({
    #Suppress warnings about aesthetic z that happens when ggtern and ggplot2 interact
    p <- ggtern::ggtern(p_hull, ggplot2::aes(x = p1, y = p2, z = p3)) +
      ggplot2::geom_polygon(
        data = p_hull,
        ggplot2::aes(x = p1, y = p2, z = p3),
        fill = "#FF000044",
        alpha = 0.3
      ) +
      ggtern::theme_bw() +
      ggtern::theme_showarrows() +
      ggplot2::xlab("Healthy") +
      ggplot2::ylab("Ill") +
      ggtern::zlab("Dead") +
      ggplot2::ggtitle(title)

    #Add point estimate
    p <- p +
      ggplot2::geom_point(
        data = p_hat,
        ggplot2::aes(x = p1, y = p2, z = p3),
        color = "black",
        size = 2,
        inherit.aes = FALSE
      )
  })

  #Return object
  return(p)
}
