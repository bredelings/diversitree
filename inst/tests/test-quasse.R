source("helper-diversitree.R")

context("QuaSSE")

test_that("quasse", {

## Build tree:
lambda <- function(x) sigmoid.x(x, 0.1, 0.2,  0, 2.5)
mu <- function(x) constant.x(x, 0.03)
char <- make.brownian.with.drift(0, 0.025)

## set.seed(1)
## phy <- tree.quasse(c(lambda, mu, char), max.taxa=15, x0=0,
##                    single.lineage=FALSE, verbose=TRUE)
load("phy.Rdata")

pars <- c(.1, .2, 0, 2.5, .03, 0, .01)

sd <- 1/200
control.C.1 <- list(dt.max=1/200)
control.C.2 <- c(control.C.1, tips.combined=TRUE)
control.M.1 <- list(method="mol")
control.R.1 <- list(dt.max=1/200, method="fftR")

if (check.fftC(FALSE)) {
lik.C.1 <- make.quasse(phy, phy$tip.state, sd, sigmoid.x, constant.x,
                       control.C.1)
##lik.C.2 <- make.quasse(phy, phy$tip.state, sd, sigmoid.x, constant.x,
##                       control.C.2)
lik.M.1 <- make.quasse(phy, phy$tip.state, sd, sigmoid.x, constant.x,
                       control.M.1)
lik.R.1 <- make.quasse(phy, phy$tip.state, sd, sigmoid.x, constant.x,
                       control.R.1)

## Conditioned references remove the historical extra dx^2 factor; dx=5*trait range/1024.
expect_that(lik.C.1(pars), equals(-52.90819039933757))
## expect_that(lik.R.1(pars), equals(-52.908190399391778)) # slow...
expect_that(lik.M.1(pars), equals(-52.876474686224441))

## Now, test root treatment:
expect_that(lik.C.1(pars, root=ROOT.FLAT), equals(-55.620545547334686))
expect_that(lik.C.1(pars, condition.surv=FALSE), equals(-55.068978927659984))
root.f <- function(x)
  dnorm(x, mean(phy$tip.state), sd(phy$tip.state))
expect_that(lik.C.1(pars, root=ROOT.GIVEN, root.f=root.f),
            equals(-53.969183708607389))

## With drift:
pars2 <- pars
pars2[6] <- 0.01
expect_that(lik.C.1(pars2), equals(-52.884261834967347))
}
})
