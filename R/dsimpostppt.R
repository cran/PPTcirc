#' @title Posterior projected Polya Tree distribution
#'
#'
#' @importFrom graphics boxplot hist lines matplot par
#' @importFrom stats acf dgamma dnorm qnorm quantile rgamma rnorm runif
#' @importFrom progress progress_bar
#'
#' @description Performs posterior inference for a given a circular dataset with the Projected Polya Tree via a MCMC algorithm.
#' @param datafile the data from which the estimate is to be computed. The object is circular or will be coerced to circular.
#' @param units units of the support: "radians", "degrees" or "hours".
#' @param mm number of finite levels of the Polya tree
#' @param mu mean vector of the projected bivariate normal centering distribution.
#' @param sig precision of the projected bivariate normal centering distribution.
#' @param aa alpha. Standard deviation parameter of the projected Polya tree.
#' @param delta controls of the speed at which the variances of the branching probabilities move down in the tree,  rho(m)=m^delta.
#' @param it number of iterations for MCMC.
#' @param bi number of burn in iterations for MCMC.
#' @param ti thinning parameter of the MCMC chain.
#' @param kapa tunning parameter in the MH proposal distribution for the latent resultants R.
#' @param ha logical. If TRUE alpha will be assigned Ga(c0,c1) hyper-prior distribution.
#' @param hm logical. If TRUE mu will be assigned N(mu0,taum) independent hyper-prior distributions for each coordinate.
#' @param c0,c1 shape and rate hyper-parameters of the gamma prior distribution for alpha. These will be used only when ha=1.
#' @param iota tunning parameter in the MH proposal distribution for alpha.
#' @param mu0,taum mean and precision hyper-parameters of the independent normal prior distribution for each coordinate of mu. These will be used only when hm=1.
#' @param control.circular the attribute used to coerced the resulting. object. See circular.
#'
#' @usage dsimpostppt(datafile,units = c("radians", "degrees", "hours"),
#' mm = 4, mu = c(0, 0), sig = 1, aa = 1, delta = 1.1,
#' it = 500, bi = 50, ti = 2, kapa = 0.5, ha = 0, hm = 0,
#' c0 = 1, c1 = 2, iota = 6, mu0 = 0, taum = 1, control.circular = list())
#'
#' @return An object of class postppt.circ whose underlying structure is a list containing the following components:
#' \item{x}{points where the density is evaluated.}
#' \item{predictive}{predicitive density estimated with the projected Polya tree.}
#' \item{quantile2.5 quantile97.5}{lower and upper 95\% credible interval limits.}
#' \item{stats}{descriptive statistics: mean direction and concentration of each MCMC density.}
#' \item{cpo}{conditional predictive ordinate statistic for the data.}
#' \item{LMPL}{logarithm of the pseudo marginal likelihood statistic.}
#' \item{aa.sims}{vector of simulated alphas when ha=1.}
#' \item{mu.sims}{matrix of simulated bivariate means when hm=1.}
#' \item{acceptancerate}{Acceptance rate of MH step for the latent resultants.}
#' \item{acceptancerate_aa}{Acceptance rate of MH step for alpha.}
#' \item{data}{original dataset.}
#' @export
#'
#' @seealso \code{\link[PPTcirc]{postppt.plot}}, \code{\link[PPTcirc]{postppt.summary}}
#'
#' @examples data(tapir)
#' #It is advised to increase the number of iterations for a better fitting
#' z1 <- dsimpostppt(tapir, units = "radians", it = 5, ti =1, bi=0, ha = 1, hm =1)
#' class(z1)
#' length(z1$acceptancerate)
#' z1$acceptancerate
#' \donttest{
#' postppt.summary(z1)
#' postppt.plot(z1, plot.type= "line" , ylim = c(0,0.8))
#' }
#' @references  Nieto-Barajas, L.E. & Nunez-Antonio, G. (2019). Projected Polya tree. https://arxiv.org/pdf/1902.06020.pdf
#'
dsimpostppt <- function(datafile, units = c("radians", "degrees", "hours"), mm=4,mu=c(0,0),sig=1,aa=1,delta=1.1,
                        it=500,bi=50,ti=2,kapa=0.5,
                        ha=0,hm=0, c0=1, c1=2, iota=6, mu0=0, taum=1,
                        control.circular=list()){

  #---------------------------------------------
  units <- match.arg(units)
  datafile <- as.vector(datafile)

  if(is.null(units)){
    stop("'units' must be 'radians' or 'deegrees' or 'hours")
  }
  if (!is.null(units)) {
    if (units == "degrees") {
      th <- datafile/180 * pi
    }
    else if (units == "hours") {
      th <- datafile/12 * pi
    }
    else if (units== "radians"){
      th <- datafile
    }

  }

  #---------------------------------------------
  jm <- 2^mm;km <- 2^mm;n<-length(datafile)
  rh <- NULL; q <- NULL
  aux <- c(NA,NA, NA, NA)
  xh <- matrix(rep(NA, 2*n), ncol = 2, byrow = TRUE)
  vsim <- seq(from = bi+1,to = it, by = ti)
  name <- deparse(substitute(datafile))
  jj <- NULL; kk <- NULL
  pb <- matrix(data = NA, nrow = jm, ncol = km, byrow = FALSE, dimnames = NULL)
  xh1 <- c(NULL, NULL); ft <- NULL; fh <- NULL; y <- c(NA, NA)
  a1 <- matrix(data = NA, nrow = mm, ncol = 4, byrow = FALSE, dimnames = NULL)
  aratea <- 0
  vt <- NULL
  mt <- NULL

  #Initial values
  arate <- rep(0, n)
  cpo <- rep(0,n)
  l0 <- 100
  ll <- 100
  FT <-  matrix(data = NA, nrow = ceiling((it-bi)/ti), ncol = ll, byrow = FALSE, dimnames = NULL)

  dt <- 2*pi/(ll-1)
  t<-seq(0,2*pi,length.out=ll)

  dr <- (sqrt(mu[1]^2+ mu[2]^2)+4)/ll
  r<-seq(dr, sqrt(mu[1]^2+ mu[2]^2)+4,length.out=ll)

  #Definition of alpha parameters
  a <-  matrix(data = NA, nrow = mm, ncol = 4, byrow = FALSE, dimnames = NULL)
  for(m in 1:mm){
    a[m,1:4]<- aa*(m^delta)*rep(1,4)
  }

  #Identifying grid location for th
  ih <- NULL
  for(i in 1:n){
    j<- 1
    while(t[j]<th[i]){
      j <- j+1
    }
    ih[i] <- j
  }

  #Deinition of partitions
  b <- pt.bipartition(mm, l0, mu, sig)
  b1 <- b$b1
  b2 <- b$b2

  #Sampling auxiliary lengths
  rh<-runif(n,0,3)
  xh[,1]<-rh*cos(th)
  xh[,2]<-rh*sin(th)

  #Computing statistics
  nc_df <- post.statistics(mm,xh, b1, b2,n)


  #-------------- MCMC ---------------------------
  #Progress bar
  progbar <- progress_bar$new(total = it)

  #Defining returns
  ss <- 1
  sim_aa <- vector("numeric", it)
  sim_mu <- matrix(data = NA, nrow = it, ncol = 2, byrow = FALSE, dimnames = NULL)


  for(s in 1:it){

    progbar$tick()
    Sys.sleep(1 / 100)

    #Simulating from y(m,j,k)
    y_df <- post.pt.branchprob(mm, a, nc_df)

    #Simulating from aa
    if(ha == 1){
      aa_sim <- pt.simaa(aa,a,c0,c1,iota,delta,mm, aratea, y_df)
      aa <- aa_sim$aa
      a <- aa_sim$a
      aratea <- aa_sim$aratea
    }

    #Computing PT probabilities
    pb <- pt.bipartprob(mm, y_df, pb)

    #Simulating from rh
    for(i in 1:n){
      rh1 <- rgamma(1, kapa, kapa/rh[i])
      xh1[1] <- rh1*cos(th[i])
      xh1[2] <- rh1*sin(th[i])
      q <- log(pt.density(jm, pb[1:jm,1:jm], b1[mm, 1:(jm+1)], b2[mm, 1:(jm+1)], mu, sig, xh1)) + log(rh1)
      q <- q - log(pt.density(jm, pb[1:jm,1:jm], b1[mm, 1:(jm+1)], b2[mm, 1:(jm+1)], mu, sig, xh[i, 1:2])) - log(rh[i])
      q <- q + dgamma(rh[i], kapa, kapa/rh1, log = TRUE)
      q <- q - dgamma(rh1, kapa, kapa/rh[i], log = TRUE)

      u <- runif(1)
      if(log(u)< q){
        rh[i] <- rh1
        xh[i, 1:2] <- xh1[1:2]
        arate[i] <- arate[i] + 1
      }
    }

    #Simulating from mu
    mu <- pt.simmu(n,taum,mu0,rh, th)

    #Re-defining partitions
    if(hm==1){
      b <- pt.bipartition(mm, l0, mu, sig)
      b1 <- b$b1
      b2 <- b$b2
    }

    #Re-computing statistics
    nc_df <- post.statistics(mm,xh, b1, b2,n)

    #Evaluating marginal density
    for(i in 1:ll){
      ft[i] <- 0
      for(l in 1:ll){
        y[1] <- r[l]*cos(t[i])
        y[2] <- r[l]*sin(t[i])
        ft[i] <- ft[i] + pt.density(jm, pb[1:jm, 1:jm], b1[mm, 1:(jm+1)], b2[mm, 1:(jm+1)], mu, sig,y)*r[l]*dr
      }
    }



    if(s %in% vsim){

      FT[ss,] <- ft
      ss <- ss+1

    }

    if(ha==1){
      sim_aa[s] <- aa
    }
    if(hm==1){
      sim_mu[s,] <- mu
    }


    #Mean direction and concentration
    alpha1 <- dt*sum(cos(t[1:ll])* ft[1:ll])
    beta1 <- dt*sum(sin(t[1:ll])* ft[1:ll])
    vt[s] <- sqrt(alpha1^2+ beta1^2)
    mt[s] <- ppt.tan(alpha1, beta1)

    #Computing cPO
    if(s %in% vsim){
      fh <- NULL
      for(i in 1:n){
        fh[i] <- ft[ih[i]]
        cpo[i] <- (cpo[i]*(ss-1)+1/fh[i])/ss
      }
    }



    #setTxtProgressBar(pb_, s)
  }

  #Computing gof measures
  ss <- ss-1
  lmpl <- -sum(log(cpo))


  #Writing acceptance rate
  if(ha==1){
    acceptrate1 <- arate/it
    acceptrate2 <- aratea/it
  }else{
    acceptrate1 <- arate/it
  }

  #The return of the function
  predictive <- colMeans(FT)
  quantile2.5 <-  apply(FT, 2, quantile, p =0.025)
  quantile97.5 <- apply(FT, 2, quantile, p =0.975)


  #close(pb_)

  #--------------------------------------------

  if (!is.null(units)) {
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

  }

  stats <- data.frame(mean.direction = mt.cc, concentration = vt.cc)

 #--------------------------------------------

  if(ha==1 & hm==1){
    return(structure(list(x = t.cc, predictive = predictive, quantile2.5 = quantile2.5,
                quantile97.5 = quantile97.5, stats = stats,
                cpo = cpo, LMPL = lmpl,
                aa.sims = sim_aa, mu.sims = sim_mu, acceptancerate= acceptrate1,
                acceptancerate_aa = acceptrate2, data = datafile),class = "postppt.circ"))
  }else if(ha ==1 & hm==0){
    return(structure(list(x = t.cc, predictive = predictive, quantile2.5 = quantile2.5,
                quantile97.5 = quantile97.5, stats = stats,
                cpo = cpo, LMPL = lmpl,  aa.sims = sim_aa, acceptancerate= acceptrate1,
                acceptancerate_aa = acceptrate2, data = datafile), class = "postppt.circ"))
  }else if(ha == 0 & hm == 1){
    return(structure(list(x = t.cc, predictive = predictive, quantile2.5 = quantile2.5,
                quantile97.5 = quantile97.5, stats = stats,
                cpo = cpo, LMPL = lmpl, mu.sims = sim_mu,
                acceptancerate= acceptrate1, data = datafile), class = "postppt.circ"))
  }else{
    return(structure(list(x = t.cc, predictive = predictive, quantile2.5 = quantile2.5,
                quantile97.5 = quantile97.5, stats = stats,
                cpo = cpo, LMPL = lmpl, acceptancerate= acceptrate1, data= datafile), class = "postppt.circ"))
  }
}
