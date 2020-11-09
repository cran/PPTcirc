#' Summary statistics for the post projected Polya tree
#'
#' @description Extracts mean, quantiles 2.5\% and 97.5\% of the mean direction and concentration.
#' @param postppt.circ object returned by \code{dsimpostppt} function.
#' @examples z1 <- dsimpostppt(tapir, units = "radians", it = 5, ti =1, bi=0)
#' postppt.summary(z1)
#' @return table of descriptive statistics.
#' @export
#'
postppt.summary <- function(postppt.circ){

  if (class(postppt.circ)!="postppt.circ")
    stop("object must be class postppt.circ")

  t <- postppt.circ$x
  ll <- length(t)
  if(t[ll]<=6.3){
    units <- "radians"
  }else if(t[ll]<=24){
    units <- "hours"
  }else if(t[ll]<=360){
    units <- "degrees"
  }

  stats <- postppt.circ$stats
  mean.direction <- circular(stats[,1], units = units , type = "angles",
                             modulo = "asis", template= "none", rotation = "counter", zero =0)
  concentration <- circular(stats[,2], units = units , type = "angles",
                            modulo = "asis", template= "none", rotation = "counter", zero =0)

  if(units == "radians"){
    mean_ <- c(circ.mean2(mean.direction),
               circ.mean2(concentration))
  }else if(units == "degrees"){
   mean_ <- c(circ.mean2(stats[,1]/180 * pi), circ.mean2(stats[,2]/180*pi))
  }else if(units == "hours"){
    mean_ <- c(circ.mean2(stats[,1]/12 * pi), circ.mean2(stats[,2]/12*pi))
  }

  quantile2.5 <-  c(circular::quantile.circular(mean.direction, 0.025),
                    circular::quantile.circular(concentration, 0.025))
  quantile97.5 <-  c(circular::quantile.circular(mean.direction, 0.975),
                     circular::quantile.circular(concentration, 0.975))

  ppt.moments <- rbind(mean_, quantile2.5, quantile97.5)
  rownames(ppt.moments) <- c('mean', 'quantile2.5',"quantile97.5" )
  colnames(ppt.moments) <- c("mean.direction","concentration" )

  if(is.null(postppt.circ$mu.sims) & is.null(postppt.circ$aa.sims)){
    ppt.moments <- ppt.moments
  }else if(is.null(postppt.circ$mu.sims) &   !is.null(postppt.circ$aa.sims)){
    sum.a <- c(mean(postppt.circ$aa.sims), quantile(postppt.circ$aa.sims, 0.025), quantile(postppt.circ$aa.sims, 0.975))
    ppt.moments <- cbind(ppt.moments, sum.a)
    colnames(ppt.moments) <- c("mean.direction","concentration" ,"alpha")
  }else if(!is.null(postppt.circ$mu.sims) &   is.null(postppt.circ$aa.sims)){
    sum.mu1 <- c(mean(postppt.circ$mu.sims[,1]), quantile(postppt.circ$mu.sims[,1], 0.025), quantile(postppt.circ$mu.sims[,1], 0.975))
    sum.mu2 <- c(mean(postppt.circ$mu.sims[,2]), quantile(postppt.circ$mu.sims[,2], 0.025), quantile(postppt.circ$mu.sims[,2], 0.975))
    ppt.moments <- cbind(ppt.moments, sum.mu1, sum.mu2)
    colnames(ppt.moments) <- c("mean.direction","concentration" ,"mu1", "mu2")
  }else{
    sum.a <- c(mean(postppt.circ$aa.sims), quantile(postppt.circ$aa.sims, 0.025), quantile(postppt.circ$aa.sims, 0.975))
    sum.mu1 <- c(mean(postppt.circ$mu.sims[,1]), quantile(postppt.circ$mu.sims[,1], 0.025), quantile(postppt.circ$mu.sims[,1], 0.975))
    sum.mu2 <- c(mean(postppt.circ$mu.sims[,2]), quantile(postppt.circ$mu.sims[,2], 0.025), quantile(postppt.circ$mu.sims[,2], 0.975))
    ppt.moments <- cbind(ppt.moments, sum.a, sum.mu1, sum.mu2)
    colnames(ppt.moments) <- c("mean.direction","concentration" ,"alpha", "mu1", "mu2")

  }
  return(ppt.moments)
}
