
#' @title Plot method for prior projected Polya tree
#'
#' @description Plots density paths of simulated prior projected Polya tree, mean direction and concentration.
#'
#' @importFrom circular circular
#' @importFrom methods is
#' @importFrom grDevices rainbow
#' @usage priorppt.plot(priorppt.circ, n.path="all",
#' plot.type = c("circle", "line", "summary"),control.circular = list(),
#' shrink=1, tol = 0.04,ylim)
#' @param priorppt.circ object returned by \code{dsimpriorppt} function.
#' @param n.path "all" plots all the simulated paths or numeric parameter indicates the simulation path of the priorppt.circ object that will be plot.
#' @param plot.type type of plot to be drawn:
#' "circle" for circular plot,
#' "line" for linear plot and
#' "summary" for boxplot of mean direction and concentration.
#' @param control.circular  attributes of circular object in order to draw the circle.See \code{\link[circular]{circular}}.
#' @param shrink parameter that controls the size of the plotted circle. Default is 1. Larger values shrink the circle, while smaller values enlarge the circle.
#' @param tol proportion of white space at the margins of plot.
#' @param ylim range to be encompassed by "y" axis.
#'
#' @examples z <- dsimpriorppt(mu = c(0,1), nsim = 5, units = "degrees")
#' priorppt.plot(z, plot.type = "circle",shrink =0.5, tol = 4)
#' priorppt.plot(z, plot.type = "line")
#' priorppt.plot(z, plot.type = "summary")
#'
#' @seealso \link[base]{plot}, \code{\link[circular]{plot.density.circular}}
#'
#' @return Circular plot of simulated paths when plot.type = "circle". Linear plot of simulated paths for plot.type = "line".
#' Boxplot of mean direction and concentration for plot.type = "summary"
#'
#' @export
#'
priorppt.plot <- function(priorppt.circ, n.path="all",plot.type = c("circle", "line", "summary"),
                          control.circular = list(),
                          shrink=1, tol = 0.04,ylim=NULL){


  if (!is(priorppt.circ,"priorppt.circ")){
    stop("object must be class priorppt.circ")
  }


  nsim <- ncol(priorppt.circ$ppt.sim)
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

  if(plot.type == "circle") {
    t <- priorppt.circ$x
    ll <- length(t)
    if(t[ll]<=6.3){
     units <- "radians"
    }else if(t[ll]<=24){
      units <- "hours"
    }else if(t[ll]<=360){
      units <- "degrees"
    }
    df <- data.frame(priorppt.circ$ppt.sim)
    x <- circular::circular(t, type = "angles",units = units, template = "none",  modulo = "asis", zero=0, rotation= "counter")

    if(is.character(n.path)){
      cl <- rainbow(nsim)
      plot(circular::circular(0, type = "angles",units = units, template = "none",  modulo = "asis", zero=0, rotation= "counter"), shrink=shrink, tol =tol, cex = 0.5)
      invisible(lapply(1:nsim, function(i) lines(x, df[,i],col = cl[i],type = 'l')))

      #suppressWarnings(matplot(x, df[,-c(1)], shrink=shrink,type = "l",ylab = "", tol = tol, main="Paths of projected Polya tree", xlab=""))


    }
    else if(is.double(n.path)){

      if(n.path <= nsim){
        y1 <- priorppt.circ$ppt.sim[,n.path]
        cl <- rainbow(nsim)
        plot(circular::circular(0, type = "angles",units = units, template = "none",  modulo = "asis", zero=0, rotation= "counter"), shrink=shrink, tol=tol, cex = 0.5)
        invisible(lapply(1:length(y1), function(i) lines(x, y1,col = cl[i],type = 'l')))
      }else{
        print("Error: Please specify the number of a simulated path" )
      }



    }
    else print("Error: Please specify number parameter (the number of the simulation or 'all')")
  }


  else  if(plot.type == "line"){
    df <- data.frame(priorppt.circ$x, priorppt.circ$ppt.sim)

    if(n.path == "all"){
      matplot(df[,1], df[,-c(1)], type="l", lty=1:nsim, pch=1, col=1:nsim,  xlab= "theta", ylab = "f(theta)")
      }
    else if(is.double(n.path)){

      if(n.path <= nsim){
        plot(df$priorppt.circ.x, df$sim1, type="l", xlab= "theta", ylab = "f(theta)")

      }else{
        print("Error: Please specify the number of a simulated path" )
      }


    }

    else print("Error: Please specify number parameter (one or all)")

  }
  else  if(plot.type == "summary"){

    boxplot(priorppt.circ$stats, main = "Boxplot of moments of projected Polya tree",
            ylim = ylim)

  }else print("Error: Please specify a type parameter (line, circular or summary type of graph)")





}

