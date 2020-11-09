r.dirichlet <- function (n, alpha)
{
  ll <- length(alpha)
  x <- matrix(rgamma(ll * n, alpha), ncol = ll, byrow = TRUE)
  sim <- x %*% rep(1, ll)
  rd <- x/as.vector(sim)

  return(rd)
}

d.dirichlet <- function (x, alpha)
{
  aux_dirichlet <- function(x, alpha) {
    logD <- sum(lgamma(alpha)) - lgamma(sum(alpha))
    s <- (alpha - 1) * log(x)
    s <- ifelse(alpha == 1 & x == 0, -Inf, s)
    exp(sum(s) - logD)
  }
  if (!is.matrix(x))
    if (is.data.frame(x))
      x <- as.matrix(x)
    else x <- t(x)
    if (!is.matrix(alpha))
      alpha <- matrix(alpha, ncol = length(alpha), nrow = nrow(x),
                      byrow = TRUE)
    if (any(dim(x) != dim(alpha)))
      stop("Mismatch between dimensions of 'x' and 'alpha'.")
    pd <- vector(length = nrow(x))
    for (i in 1:nrow(x)) pd[i] <- aux_dirichlet(x[i, ], alpha[i,
    ])
    pd[apply(x, 1, function(z) any(z < 0 | z > 1))] <- 0
    pd[apply(x, 1, function(z) all.equal(sum(z), 1) != TRUE)] <- 0

    return(pd)
}

pt.density <- function(jm, pb, b1,b2, mu, sig, x){
  j <- 1
  if(x[1]<b1[1]){
    j <- 1
  }else if(x[1]> b1[jm]){
    j <- jm
  }
  else{
    while(b1[j+1]<x[1] & !is.na(b1[j+1]) ){
      j <- j +1
    }
  }
  k <- 1
  if(x[2]<b2[1]){
    k <- 1
  }else if(x[2]>b2[jm]){
    k <- jm
  }else{
    while(b2[k+1]< x[2] & !is.na(b2[k+1])){
      k <- k+1
    }
  }
  pt.density <-  pb[j,k]*jm^2*dnorm(x[1], mu[1], sig)*dnorm(x[2], mu[2], sig)
  return(pt.density)
}

pt.bipartition <- function(mm, l0, mu, sig){
  b1 <- matrix(data = NA, nrow = mm, ncol = (2^mm)+1, byrow = FALSE, dimnames = NULL)
  b2 <- matrix(data = NA, nrow = mm, ncol = (2^mm)+1, byrow = FALSE, dimnames = NULL)
  for(m in 1:mm){
    for(j in 0:2^m){
      if(j == 0){
        b1[m,j+1] <- -l0*sig + mu[1]
        b2[m,j+1] <- -l0*sig + mu[2]
      }else if(j == 2^m){
        b1[m,j+1] <- l0*sig + mu[1]
        b2[m,j+1] <- l0*sig + mu[2]
      }else{
        b1[m,j+1] <- qnorm(j/(2^m), mu[1], sig)
        b2[m,j+1] <- qnorm(j/(2^m), mu[2], sig)
      }
    }
  }
  return(list(b1=b1,b2=b2))
}

pt.branchprob <- function(mm, a){
  y_df <- data.frame( y=rep(NA,(1-4^(mm+1))/(-3)-1), m_=rep(NA,(1-4^(mm+1))/(-3)-1), j_=rep(NA,(1-4^(mm+1))/(-3)-1), k_=rep(NA,(1-4^(mm+1))/(-3)-1) )
  i <- 1
  for(m in 1:mm){
    for(j in 1:2^(m-1)){
      for(k in 1:2^(m-1)){
        x <- r.dirichlet(1,a[m,1:4])
        y_df$y[i] <- x[1]
        y_df$m_[i] <- m
        y_df$j_[i] <- 2*j-1
        y_df$k_[i] <- 2*k-1
        i <- i + 1
        y_df$y[i] <- x[2]
        y_df$m_[i] <- m
        y_df$j_[i] <- 2*j-1
        y_df$k_[i] <- 2*k
        i <- i + 1
        y_df$y[i] <- x[3]
        y_df$m_[i] <- m
        y_df$j_[i] <- 2*j
        y_df$k_[i] <- 2*k-1
        i <- i + 1
        y_df$y[i] <- x[4]
        y_df$m_[i] <- m
        y_df$j_[i] <- 2*j
        y_df$k_[i] <- 2*k
        i <- i + 1
      }
    }
  }
  return(y_df)
}

pt.bipartprob <- function(mm, y_df, pb ){
  jj <- kk <- NULL; jm <- 2^mm
  pb <- matrix(data = NA, nrow = jm, ncol = jm, byrow = FALSE, dimnames = NULL)
  for(j in 1:jm){
    for(k in 1:jm){
      jj <- j
      kk <- k
      pb[j,k] <- y_df[(y_df$m_ == mm & y_df$j_ == jj & y_df$k_ == kk),1]
      for(l in 2:mm){
        jj <- ceiling((jj)/2)
        kk <- ceiling((kk)/2)
        pb[j,k] <- pb[j,k]*y_df[(y_df$m_ == mm-l+1 & y_df$j_ == jj & y_df$k == kk),1]
      }
    }
  }
  return(pb)
}

dpptcirc <-function(mm=4,mu = c(0,0), sig=1, ll = 100, aa=1, delta=1.1){

  l0 <-  100;dt <- 2*pi/ll;dr <- (sqrt(mu[1]^2+ mu[2]^2)+4)/ll
  t <- seq(from = 0, to= 2*pi*(1-1/ll), by = 2*pi/ll)
  r <- seq(from = 0, to = (sqrt(mu[1]^2+ mu[2]^2)+4)*(1-1/ll), by = dr )
  jm <- 2^mm; km <- 2^mm


  #Definition of alpha parameters
  a <-  matrix(data = NA, nrow = mm, ncol = 4, byrow = FALSE, dimnames = NULL)
  for(m in 1:mm){
    a[m,1:4]<- aa*(m^delta)*rep(1,4)
  }

  #Definition of partitions
  b <- pt.bipartition(mm, l0, mu, sig)
  b1 <- b$b1
  b2 <- b$b2

  #Simulate branch probabilities
  y_df <- pt.branchprob(mm, a)

  #Calculate partition probability
  pb <- pt.bipartprob(mm, y_df, pb)

  #Calculate density
  ft <- NULL
  y <- c(NA, NA)
  for(i in 1:ll){
    ft[i] <- 0
    for(l in 1:ll){
      y[1] <- r[l]*cos(t[i])
      y[2] <- r[l]*sin(t[i])
      ft[i] <- ft[i] + pt.density(jm, pb[1:jm, 1:jm], b1[mm, 0:jm], b2[mm, 0:jm], mu, sig,y)*r[l]*dr
    }
  }
  return(ft)
}

ppt.tan <- function(a,b){
  if(a>0 & b>=0){
    ppt.tan <- atan(b/a)
  }else if(a == 0 & b>0 ){
    ppt.tan <- pi/2
  }else if(a<0){
    ppt.tan <- atan(b/a)+ pi
  }else if(a>=0 & b<0){
    ppt.tan <- atan(b/a) + 2*pi
  }else{
    print("Error, undefinied angle")
    ppt.tan <- NA
  }
  return(ppt.tan)
}



post.statistics <- function(mm,xh, b1, b2,n){
  nc_df <- data.frame( nc=rep(NA,(1-4^(mm+1))/(-3)-1), m_=rep(NA,(1-4^(mm+1))/(-3)-1), j_=rep(NA,(1-4^(mm+1))/(-3)-1), k_=rep(NA,(1-4^(mm+1))/(-3)-1) )
  l <- 1
  for(m in 1:mm){
    for(j in c(2:(2^m+1))){
      for(k in c(2:(2^m+1))){
        nc_df$nc[l] <- 0
        for(i in 1:n){
          nc_df$nc[l] <- nc_df$nc[l] +  (b1[m,j-1]<xh[i,1] & xh[i,1]<=b1[m,j])*(b2[m,k-1]<xh[i,2] & xh[i,2]<= b2[m,k])
        }
        nc_df$m_[l] <- m
        nc_df$j_[l] <- j-1
        nc_df$k_[l] <- k-1
        l <- l+1

      }
    }
  }
  return(nc_df)
}

post.pt.branchprob <- function(mm, a, nc_df){
  y_df <- data.frame( y=rep(NA,(1-4^(mm+1))/(-3)-1), m_=rep(NA,(1-4^(mm+1))/(-3)-1), j_=rep(NA,(1-4^(mm+1))/(-3)-1), k_=rep(NA,(1-4^(mm+1))/(-3)-1) )
  ncc <- c(NA, NA, NA, NA)
  i <- 1
  for(m in 1:mm){
    for(j in 1: 2^(m-1)){
      for(k in 1:2^(m-1)){
        ncc[1] <- nc_df$nc[nc_df$m_ == m & nc_df$j_ == 2*j-1 & nc_df$k_ == 2*k-1]
        ncc[2] <- nc_df$nc[nc_df$m_ == m & nc_df$j_ == 2*j-1 & nc_df$k_ == 2*k]
        ncc[3] <- nc_df$nc[nc_df$m_ == m & nc_df$j_ == 2*j & nc_df$k_ == 2*k-1]
        ncc[4] <- nc_df$nc[nc_df$m_ == m & nc_df$j_ == 2*j & nc_df$k_ == 2*k]

        x <- r.dirichlet(1,a[m,1:4] + ncc)
        y_df$y[i] <- x[1]
        y_df$m_[i] <- m
        y_df$j_[i] <- 2*j-1
        y_df$k_[i] <- 2*k-1
        i <- i + 1
        y_df$y[i] <- x[2]
        y_df$m_[i] <- m
        y_df$j_[i] <- 2*j-1
        y_df$k_[i] <- 2*k
        i <- i + 1
        y_df$y[i] <- x[3]
        y_df$m_[i] <- m
        y_df$j_[i] <- 2*j
        y_df$k_[i] <- 2*k-1
        i <- i + 1
        y_df$y[i] <- x[4]
        y_df$m_[i] <- m
        y_df$j_[i] <- 2*j
        y_df$k_[i] <- 2*k
        i <- i + 1
      }
    }
  }

  return(y_df)
}
pt.simmu <- function(n,taum,mu0,rh, th){
  mu <- c(NA, NA)
  tau1 <- n + taum
  mu1 <- taum*mu0+sum(rh[1:n]*cos(th[1:n]))
  mu1 <- mu1/(n+taum)
  mu[1] <- rnorm(1,mu1, 1/tau1 )
  mu1=taum*mu0+sum(rh[1:n]*sin(th[1:n]))
  mu1=mu1/(n+taum)
  mu[2]=rnorm(1, mu1,1/tau1)
  return(mu)
}

pt.simaa <- function(aa,a,c0,c1,iota,delta,mm, aratea, y_df){
  aux <- rep(NA, 4)
  a1 <- matrix(data = NA, nrow = mm, ncol = 4, byrow = FALSE, dimnames = NULL)
  aa1 <- rgamma(1, shape=iota, rate =iota/aa)
  q <- dgamma(aa1, shape=c0, rate=c1, log = TRUE) - dgamma(aa, shape=c0, rate=c1, log = TRUE)
  q <- q + dgamma(aa, iota, iota/aa1, log = TRUE) - dgamma(aa1, iota, iota/aa, log = TRUE)
  for(m in 1:mm){
    a1[m, 1:4] <- aa1*(m^delta)*rep(1,4)
    for(j in 1:2^(m-1)){
      for(k in 1:2^(m-1)){
        aux[1] <- y_df$y[y_df$m_ == m & y_df$j_ == (2*j-1) & y_df$k_ == (2*k-1)]
        aux[2] <- y_df$y[y_df$m_ == m & y_df$j_ == (2*j-1) & y_df$k_ == (2*k)]
        aux[3] <- y_df$y[y_df$m_ == m & y_df$j_ == (2*j) & y_df$k_ == (2*k-1)]
        aux[4] <- y_df$y[y_df$m_ == m & y_df$j_ == (2*j) & y_df$k_ == (2*k)]
        q <-  q + log(d.dirichlet(aux,a1[m,1:4])) - log(d.dirichlet(aux,a[m,1:4]))
      }
    }
  }
  u <- runif(1,0,1)
  if(log(u)<=q){
    aa <- aa1
    a <- a1
    aratea <- aratea+1
  }
  return(list(aa=aa, a=a,aratea=aratea))
}


circ.mean2 <- function (x)
{
  sinr <- sum(sin(x))
  cosr <- sum(cos(x))
  circmean <- ppt.tan(cosr, sinr)
  circmean
}


