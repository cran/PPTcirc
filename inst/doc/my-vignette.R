## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  warnings = FALSE
)

## ----setup, warnings = FALSE--------------------------------------------------
library(PPTcirc)

## ---- echo = TRUE, fig.width=4.5, fig.height=3--------------------------------
z1 <- dsimpriorppt(nsim = 10, mm = 4,  mu = c(0,0), sig = 1, ll = 100,
                   aa = 1, delta = 1.1, units = "radians")
class(z1)
priorppt.plot(z1, plot.type = "line")
priorppt.plot(z1, plot.type = "circle", tol = 8, shrink = .13)
priorppt.summary(z1)

## ----fig.width=4, fig.height=3------------------------------------------------
mat_mu <- matrix(c(0,0,-1,0,0,1,2,2), nrow = 4, ncol = 2, byrow = TRUE)
for(i in 1:4){
  z <- dsimpriorppt(nsim = 10, mm = 4,  mu = mat_mu[i,], sig = 1, ll = 100,
                    aa = 1, delta = 1.1, units = "radians")
  priorppt.plot(z, plot.type = "line")
}

## -----------------------------------------------------------------------------
peccary <- data("peccary")


## ----fig.height=3, fig.width=4.5----------------------------------------------
set.seed(2021)
rm(list=ls())
model_1 <- dsimpostppt(peccary, units = "radians", mu = c(0,0),
                       it = 100, ti=2, bi=10)
model_2 <- dsimpostppt(peccary, units = "radians", mu = c(0,0),
                       it = 100, ti=2, bi=10, ha = 1, hm = 1)
postppt.plot(model_1, plot.type= "line" , ylim = c(0,0.8))
postppt.plot(model_2, plot.type= "line" , ylim = c(0,0.8))

## -----------------------------------------------------------------------------
model_1$LMPL
model_2$LMPL

## ---- fig.height=4, fig.width=5-----------------------------------------------

postppt.summary(model_1)
postppt.plot(model_1, plot.type = "summary")


## ----fig.height=4, fig.width=6------------------------------------------------
postppt.plot(model_2, plot.type = "a.sim" )
postppt.plot(model_2, plot.type = "mu.sim")

