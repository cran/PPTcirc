#' Summary for the prior projected Polya tree simulations
#'
#' @description Mean, quantiles 2.5\% and 97.5\% of the mean direction and concentration.
#'
#' @param priorppt.circ  object returned by  \code{dsimpriorppt} function.
#' @param units units of the support: "radians", "degrees" or "hours".
#'
#' @examples z <- dsimpriorppt(mu = c(-1,0), nsim = 5, units = "hours")
#' priorppt.summary(z)
#'
#' @return Table of descriptive statistics for mean direction and concentration.
#' @export
#'

priorppt.summary <- function(priorppt.circ, units ="radians"){

  stats <- priorppt.circ$stats
  mean.direction <- circular::circular(stats$mean.direction, units = units, type = "angles", zero = 0, rotation = "counter", template = "none", modulo= "asis" )
  concentration <- circular::circular(stats$concentration, units = units, type = "angles", zero = 0, rotation = "counter", template="none", modulo = "asis")

  mean_ <- c(circ.mean2(mean.direction),circ.mean2(concentration))
  quantile2.5 <-  c(circular::quantile.circular(mean.direction, 0.025),
                    circular::quantile.circular(concentration, 0.025))
  quantile97.5 <-  c(circular::quantile.circular(mean.direction, 0.975),
                     circular::quantile.circular(concentration, 0.975))

  ppt.moments <- rbind(mean_, quantile2.5, quantile97.5)
  rownames(ppt.moments) <- c('mean', 'quantile2.5',"quantile97.5")
  colnames(ppt.moments) <- c("mean.direction","concentration")
  return(ppt.moments)
}
