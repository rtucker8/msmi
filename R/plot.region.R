#' Plot Confidence Regions State Occupation Probabilities
#'
#' @param tprob_obj The list returned by msmi::tprobs()
#' @param time  The time point of interest at which to generate the plot
#' @param unconstrained Logical indicating whether to provide an additional plot in the unconstrained 2D space. Default is FALSE.
#'
#' @returns If unconstrained == FALSE, A plot showing the confidence region and point estimate in the probability simplex.
#' If unconstrained == TRUE, A list of plots showing the confidence region and point estimate in the (1) unconstrained 2D space and the (2) probability simplex.
#' @export
#'
#' @examples
#' imps <- msmi.impute(dat = sim.data, M = 5, n.states = 3,
#'                     prefix.states = c("event", "t"), method = "marginal")
#' tprobs <- msmi.tprobs(imp_obj = imps, times = c(1,3))
#' #Plot the confidence region for state occupation probabilities at time = 1
#' plot_region(tprob_obj = tprobs, time = 1)

plot_region <- function(tprob_obj, time, unconstrained = FALSE) {

  if (unconstrained == TRUE) {

    #Plot in the unconstrained 2D space
    theta_hat <- tprob_obj$unconstrained_estimate[as.character(time), c("theta1", "theta2")]

    #Convex hull for confidence ellipse
    r_hull <- tprob_obj$cr_list[[as.character(time)]][["unconstrained"]]
    r_hull <- as.data.frame(r_hull)
    colnames(r_hull) <- c("theta1", "theta2")

    r <- ggplot2::ggplot() +
      ggplot2::geom_vline(xintercept = 0, linewidth = 0.5, alpha = 0.8) +
      ggplot2::geom_hline(yintercept = 0, linewidth = 0.5, alpha = 0.8) +
      ggplot2::geom_polygon(data = r_hull, ggplot2::aes(x = theta1, y = theta2), fill="#FF000044",alpha = 0.3) +
      ggplot2::theme_bw() +
      ggplot2::labs(x = "Theta 1", y = "Theta 2") +
      ggplot2::xlim(min(r_hull$theta1)-2, max(r_hull$theta1)+2) +
      ggplot2::ylim(min(r_hull$theta2)-2, max(r_hull$theta2)+2) +
      ggplot2::ggtitle(paste0(as.character((1-tprob_obj$alpha)*100), "% Confidence Region in Unconstrained Space at time t = ", as.character(time)))

    r <- r + ggplot2::geom_point(data = theta_hat,
                                 ggplot2::aes(x = theta1, y = theta2),
                                 color = "black",
                                 size = 2,
                                 inherit.aes = FALSE
    )

  }

#Plot in the probability simplex

  #Convex Hull for confidence region
  p_hull <- tprob_obj$cr_list[[as.character(time)]][["p.space"]]
  p_hull <- as.data.frame(p_hull)
  colnames(p_hull) <- c("p1", "p2", "p3")

  #Point estimate to overlay
  p_hat <- tprob_obj$mi_estimate[as.character(time), c("p1", "p2", "p3")]

  title <- paste0(as.character((1-tprob_obj$alpha)*100), "% Confidence Region for State Occupation Probabilities at time t = ", as.character(time))

 #Make plot
 suppressWarnings({ #Suppress warnings about aesthetic z that happens when ggtern and ggplot2 interact
  p <- ggtern::ggtern(p_hull, ggplot2::aes(x = p1, y = p2, z = p3)) +
    ggplot2::geom_polygon(data = p_hull, ggplot2::aes(x = p1, y = p2, z = p3), fill="#FF000044",alpha = 0.3) +
    ggtern::theme_bw() +
    ggtern::theme_showarrows() +
    ggplot2::xlab("Healthy") +
    ggplot2::ylab("Ill") +
    ggtern::zlab("Dead") +
    ggplot2::ggtitle(title)

  #Add point estimate
  p <- p +
    ggplot2::geom_point(data = p_hat,
                        ggplot2::aes(x = p1, y = p2, z = p3),
                        color = "black",
                        size = 2,
                        inherit.aes = FALSE
    )

 })

  if (unconstrained == FALSE) {
    return(p)
  } else if (unconstrained == TRUE) {
    return(list(r,p))
  }
}

