#' @title Plot method for posterior projected Polya tree
#'
#' @description Plots posterior projected Polya tree estimates.
#'
#' @usage postppt.plot(postppt.circ,
#' plot.type = c("circle", "line", "summary", "a.sim", "mu.sim", "cpos"),
#' interval = TRUE, control.circular = list(),
#' shrink = 1, tol = 0.04,sep = 0.025, ylim = NULL, xlim = NULL, breaks = 12)
#'
#' @param postppt.circ object returned by the \code{dsimpostppt} function.
#' @param plot.type type of plot to be drawn:
#' "circle" for circular plot,
#' "line" for linear plot,
#' "summary" for boxplot of mean direction and concentration,
#' "cpos" for cpos scatter plot,
#' "a.sim" for summary plots of simulated alphas and
#' "mu.sim" for summary plots of simulated mu1 and mu2.
# }
#'
#' @param interval logical. If TRUE 95\% credible intervals will be shown in the circular and linear plots.
#' @param control.circular  attributes of circular object in order to draw the circle.See \code{\link[circular]{circular}}.
#' @param shrink parameter that controls the size of the plotted circle. Default is 1. Larger values shrink the circle, while smaller values enlarge the circle.
#' @param tol proportion of white space at the margins of plot.
#' @param sep constant used to specify the distance between stacked points. Default is 0.025;smaller values will create smaller spaces
#' @param ylim range to be encompassed by "y" axis.
#' @param xlim range to be encompassed by "x" axis.
#' @param breaks one of: a vector giving the breakpoints between histogram cells,
#' a function to compute the vector of breakpoints,
#' a single number giving the number of cells for the histogram,
#' a character string naming an algorithm to compute the number of cells,
#' a function to compute the number of cells.
#'
#' @importFrom methods is
#'
#' @examples \donttest{ z2 <- dsimpostppt(deer, units = "radians", it = 10, ti =1, bi=0, ha = 1)
#' postppt.plot(z2, plot.type= "line" , shrink = 1.4, tol = 1.2, ylim = c(0,0.6))
#' postppt.summary(z2)
#' postppt.plot(z2, plot.type= "cpos" )
#' postppt.plot(z2, plot.type= "circle" , shrink = 1.4, tol = 1.2)}
#'
#' @seealso \link[base]{plot}, \code{\link[circular]{plot.density.circular}}
#' @export
#'

postppt.plot <- function(postppt.circ, plot.type=c("circle", "line", "summary", "a.sim", "mu.sim", "cpos"), interval=TRUE,
                         control.circular = list(),shrink=1,
                         tol = 0.04, sep = 0.025, ylim=NULL, xlim= NULL,
                         breaks = 12){

  if (!is(postppt.circ,"postppt.circ"))
    stop("object must be class postppt.circ")

  datacircularp <- list(type="angles", units="radians", template="none", modulo="asis", zero=0, rotation="counter")
  dc <- control.circular
  if (is.null(dc$type))
    dc$type <- datacircularp$type
  if (is.null(dc$units))
    dc$units <- datacircularp$units
  if (is.null(dc$template))
    dc$template <- datacircularp$template
  if (is.null(dc$modulo))
    dc$modulo <- datacircularp$modulo
  if (is.null(dc$zero))
    dc$zero <- datacircularp$zero
  if (is.null(dc$rotation))
    dc$rotation <- datacircularp$rotation

  if (is.null(xlim))
    xlim <- range(c(postppt.circ$x, postppt.circ$data))

  if(plot.type == "line"){

    hist(postppt.circ$data,freq=FALSE,col="grey80",border=FALSE,ylim=ylim,
         xlab="theta",ylab="f(theta)",main="", xlim = xlim,breaks=breaks )
    lines(postppt.circ$x,postppt.circ$predictive,lty=1)

    if(interval==TRUE){
      lines(postppt.circ$x, postppt.circ$quantile2.5, lty =2)
      lines(postppt.circ$x, postppt.circ$quantile97.5, lty = 2)
    }

  }else if(plot.type == "circle"){
    data.c <- circular(postppt.circ$data, units = dc$units, template = dc$template, zero = dc$zero, modulo = dc$modulo, rotation = dc$rotation)
    x.c <- circular(postppt.circ$x, units = dc$units, template = dc$template, zero = dc$zero, modulo = dc$modulo, rotation = dc$rotation)
    plot(data.c, stack = TRUE ,col="grey80",
         type= "p", tol = tol,
         shrink=shrink, sep =sep)
    lines(x.c,  postppt.circ$predictive, type = "l")
    if(interval==TRUE){
      lines(x.c, postppt.circ$quantile2.5, type= "l",lty =2)
      lines(x.c, postppt.circ$quantile97.5,type="l", lty = 2)
    }
  }else if(plot.type == "summary"){
    boxplot(postppt.circ$stats, main = "Boxplot of moments of projected polya tree")

  }else if(plot.type == "a.sim"){
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar))
    par(mfrow=c(1,2))
    hist(postppt.circ$aa.sims,main = "Histogram of alpha",xlab="", col= "grey80")
    plot(postppt.circ$aa.sims, type="l", xlab= "Iteration", ylab = "Value", main = "Trace plot")
    acf(postppt.circ$aa.sims, main = "Alpha simulations")
    erg.mean  <- (cumsum(postppt.circ$aa.sims) / seq_along(postppt.circ$aa.sims))
    plot(erg.mean,type="l", xlab = "Iteration", ylab = " ", main = "Ergodic mean")
  }else if(plot.type == "mu.sim"){

    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar))
    par(mfrow=c(1,2))
    hist(postppt.circ$mu.sims[,1],main = "Histogram of mu1",xlab="", col= "grey80")
    plot(postppt.circ$mu.sims[,1], type="l", xlab= "Iteration", ylab = "mu1", main = "Trace plot")
    acf(postppt.circ$mu.sims[,1], main = "Mu1 simulations")
    erg.mean  <- (cumsum(postppt.circ$mu.sims[,1]) / seq_along(postppt.circ$mu.sims[,1]))
    plot(erg.mean,type="l", xlab = "Iteration", ylab = "mu1", main = "Ergodic mean")

    hist(postppt.circ$mu.sims[,2],main = "Histogram of mu2",xlab="", col= "grey80")
    plot(postppt.circ$mu.sims[,2], type="l", xlab= "Iteration", ylab = "mu2", main = "Trace plot")
    acf(postppt.circ$mu.sims[,2], main = "Mu2 simulations")
    erg.mean  <- (cumsum(postppt.circ$mu.sims[,2]) / seq_along(postppt.circ$mu.sims[,2]))
    plot(erg.mean,type="l", xlab = "Iteration", ylab = "mu2", main = "Ergodic mean")
    }else if(plot.type == "cpos"){

    plot(postppt.circ$cpo, type = "p", ylab = "cpo", main = "CPO graphic")

  }else print("Error: Please specify a type parameter type (line, circle, summary, a.sim, mu.sim or cpos)  of graph)")


}
