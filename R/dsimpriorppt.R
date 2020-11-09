#' Prior projected Polya tree distribution
#'
#' @description Simulates paths of prior projected Polya tree distributions centered
#' around a projected normal distribution.
#'
#' @usage dsimpriorppt(nsim = 5, mm = 4,mu = c(0, 0),
#' sig = 1, ll = 100, aa = 1, delta = 1.1, units = "radians")
#' @param nsim integer indicating the number of simulations.
#' @param mm integer indicating the number of finite levels of the Polya tree.
#' @param mu mean vector of the projected bivariate normal distribution.
#' @param sig standard deviation of the projected bivariate normal distribution. We advise to always use sig = 1.
#' @param ll 	number of equally spaced points at which the projected distribution will be evaluated.
#' @param aa alpha. Precision parameter of the  Polya tree.
#' @param delta controls of the speed at which the variances of the branching probabilities move down in the tree,  rho(m)=m^delta.
#' @param units units of the support: "radians", "degrees" or "hours".
#'
#' @examples z <- dsimpriorppt(mu = c(5,5), nsim = 5, units = "radians")
#' priorppt.plot(z, plot.type = "line")
#' summary(z$stats)
#'
#' @seealso \code{\link[PPTcirc]{priorppt.plot}}, \code{\link[PPTcirc]{priorppt.summary}}
#'
#' @return An object with class priorppt.circ whose underlying structure is a list containing the following components:
#' \item{x}{points where the density is evaluated.}
#' \item{ppt.sims}{simulated density paths of the prior projected Polya tree.}
#' \item{stats}{descriptive statistics: mean direction and concentration of each simulated density.}
#' @export
#'
#'
#' @references  Nieto-Barajas, L.E. & Nunez-Antonio, G. (2019). Projected Polya tree. https://arxiv.org/pdf/1902.06020.pdf

dsimpriorppt <- function(nsim=5 , mm=4, mu = c(0,0), sig=1, ll = 100, aa=1,
                         delta=1.1, units= "radians"){
  #Prog_bar <- txtProgressBar(min = 0, max = nsim, style = 3)
  dt <- 2*pi/ll
  t <- seq(from = 0, to= 2*pi*(1-1/ll), by = dt)
  FT <-  matrix(data = NA, nrow = ll, ncol = nsim, byrow = FALSE, dimnames = NULL)
  vt <- mt <- vector("numeric", length = nsim)

  for(i in 1:nsim){
    #Generate simulations
    FT[,i] <- dpptcirc(mm,mu = mu, sig=sig, ll = ll, aa=aa, delta=delta)
    #Generate stats
    alpha1 <- dt*sum(cos(t[1:ll])* FT[,i])
    beta1 <- dt*sum(sin(t[1:ll])* FT[,i])
    vt[i] <- sqrt(alpha1^2+ beta1^2)
    mt[i] <- ppt.tan(alpha1, beta1)
    #setTxtProgressBar(Prog_bar,i)

  }


    if (units == "degrees") {
      mt.cc <- mt/pi *180
      vt.cc <- vt/pi*180
      t.cc <- t/pi*180
    }
    else if (units == "hours") {
      mt.cc <- mt/pi*12
      vt.cc <- vt/pi*12
      t.cc <- t/pi*12
    }
    else if (units== "radians"){
      mt.cc <- mt
      vt.cc <- vt
      t.cc <- t
    }
  FT <- as.data.frame(FT)
  colnames(FT) <- paste0(rep("sim", nsim),1:nsim )
  #close(Prog_bar)
  return(structure(list(x = t.cc ,ppt.sim = FT,  stats=data.frame(mean.direction =mt.cc, concentration=vt.cc)), class = "priorppt.circ"))
}
