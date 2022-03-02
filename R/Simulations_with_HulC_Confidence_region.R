## This file contains codes for the simulation of plots for the different choices
## of the true valuation distribution function F, and their corresponding HulC
## confidence regions.


#------------------------------------------------------------------------------------
## Remove everything from the Global environment
#------------------------------------------------------------------------------------
rm(list = ls())


#------------------------------------------------------------------------------------
## Install the required packages into R if it's not already installed.
#------------------------------------------------------------------------------------
if(!require(rmutil)) install.packages("rmutil")
if(!require(ggplot2)) install.packages("ggplot2")
if(!require(GoFKernel)) install.packages("GoFKernel")
if(!require(transport)) install.packages("transport")


#------------------------------------------------------------------------------------
## Load the required packages into the R environment.
#------------------------------------------------------------------------------------
library(rmutil)     # For generating from pareto distribution
library(ggplot2)    # For the graphics/plots in R.
library(GoFKernel)  # For the inverse function of a CDF
library(transport)  # For calculating Wasserstein distance between two distribution functions.


#------------------------------------------------------------------------------------
## Source all the functions from "Functions.R" and "HulC.R"
#------------------------------------------------------------------------------------
source("Functions.R")
source("HulC.R")


#------------------------------------------------------------------------------------
## User provided values
#------------------------------------------------------------------------------------
lambda.true <- 1  # The true arrival rate of bidders in a second price auction.
auction.window <- 100  # It represents how long the auction is open for bidders.
N.auction <- 200  # Number of independent auctions of identical copies of an item.


#------------------------------------------------------------------------------------
## Different choices of true valuation distribution functions and the corresponding 
## choices of parameters and reserve prices are listed below:

## Uncomment any of the following choices to see the plots for that particular choice
## of the true valuation distribution F.
#------------------------------------------------------------------------------------

## Choice 1: Uniform distribution.
# method.in <- "unif"
# para1 <- 1
# para2 <- 20
# x.med.vec <- seq(0.1, 20, by= 0.995)
# reserve.price <- runif(N.auction, 0.1, 3)
# reserve.price.Cutoff <- 1


## Choice 2: Pareto distribution.
# method.in <- "pareto"
# para1 <- 3
# para2 <- 100
# x.med.vec <- seq(0, 20, by= 1)
# reserve.price <- runif(N.auction, 0.001, 0.1)
# reserve.price.Cutoff <- 0.05


## Choice 3: Gamma distribution.
method.in <- "gamma"
para1 <- 10
para2 <- 2
x.med.vec <- seq(0, 12, by= 12/20)
reserve.price <- runif(N.auction, 0.1, 3)
reserve.price.Cutoff <- 1


## Choice 4: Beta distribution.
# method.in <- "beta"
# para1 <- 2
# para2 <- 2
# x.med.vec <- seq(0, 1, by= 1/20)
# reserve.price <- runif(N.auction, 0.001, 0.2)
# reserve.price.Cutoff <- 0.05


## Choice 5: Piecewise uniform distribution.
# method.in <- "piecewise_unif"
# para1 <- 2
# para2 <- 4
# x.med.vec <- seq(0, 4, by = 4/20)
# reserve.price <- runif(N.auction, 0.1, 1.5)
# reserve.price.Cutoff <- 1


#------------------------------------------------------------------------------------
## Main body of codes
#------------------------------------------------------------------------------------
raw.data.list <- Data.Gen.2nd.Price.Raw(N.auction, auction.window, lambda.true, 
                                        reserve.price, method = method.in, para1, para2)
data <- Data.Gen.2nd.Price.Processed(raw.data.list)
LL <- length(data$pooled.observed.bids)
Init.Values <- initialization.2nd.price(data, reserve.price.Cutoff)
theta.init <- (1-Init.Values$F.y)/ c(1,(1-Init.Values$F.y)[-length(Init.Values$F.y)])


## MLE of F (without any modification). 
## Uncomment the below line to see how it is unstable near the boundaries, especially near 0 in our case.
# MLE.old <- MLE.2ndprice(data, lambda.in = Init.Values$lambda, theta.in = theta.init, tol = 1e-5)


## MLE with the modifications that the first few theta_i_{MLE}'s  
## defined to be the corresponding theta_initial_i's.
MLE <- MLE.2ndprice.init(data, lambda.in = Init.Values$lambda, theta.in = theta.init,
                         tol = 1e-5)


## HulC confidence regions
Hulc.Conf <- HulC1d_Ebay(data.in = raw.data.list, eval.points = x.med.vec, alpha=.1, 
                         MLE.Cutoff = reserve.price.Cutoff)


#------------------------------------------------------------------------------------
## Plot the F_{initial} estimate and F_{MLE}.
#------------------------------------------------------------------------------------

## The following command plots the F_{initial} estimate of F.
plot(Init.Values$F.x, Init.Values$F.y, type="l", col = "blue",
     main = paste("Plot of true F (black) Vs. F_init(blue) Vs. F_MLE (red), N.auction =", 
                  N.auction, ", True F =", data$true.F),
     xlab = "x", ylab = "F(x)")


## Uncomment the following command to check the unstable MLE of F (unstable near 0)
# lines(MLE.old$F.x, MLE.old$F.y, type="l", col = "green")


## The following command plots the MLE of F on top of F_{initial}.
lines(MLE$F.x, MLE$F.y, type="l", col = "red")


#------------------------------------------------------------------------------------
## Plot the true F.
## Uncomment one of the following lines according to the user's choice of the true F.
#------------------------------------------------------------------------------------

## Uncomment the following command if true F follows beta distribution.
# lines(data$pooled.data[,1], pbeta(data$pooled.data[,1], shape1 = para1, shape2 = para2), 
#       col = "black", type = "l") 


## Uncomment the following command if true F follows gamma distribution.
lines(data$pooled.data[,1], pgamma(data$pooled.data[,1], shape = para1, rate = para2), 
      col = "black", type = "l") 


## Uncomment the following command if true F follows uniform distribution.
# lines(data$pooled.data[,1], punif(data$pooled.data[,1], min = para1, max = para2), 
#       col = "black", type = "l") 


## Uncomment the following command if true F follows pareto distribution.
# lines(data$pooled.data[,1], ppareto(data$pooled.data[,1], m = para1, s = para2), 
#       col = "black", type = "l") 


## Uncomment both of the following two commands if true F follows piecewise uniform distribution.
# F0.piecewise.unif_pooled.data <- sapply(data$pooled.data[,1], 
#                                         function(y) {piecewise.unif.cdf(y, para1 = para1, para2 = para2)} )
# lines(data$pooled.data[,1], F0.piecewise.unif_pooled.data, 
#       col = "black", type = "l")


#------------------------------------------------------------------------------------
## HulC Confidence band for both F_{initial} and F_{MLE}.
#------------------------------------------------------------------------------------

points(x.med.vec, as.vector( Hulc.Conf$CI.Init[1,]), col="blue", lwd=1)
points(x.med.vec, as.vector( Hulc.Conf$CI.Init[2,]), col="blue", lwd=1)
points(x.med.vec, as.vector( Hulc.Conf$CI.MLE[1,]), col="red", lwd=1)
points(x.med.vec, as.vector( Hulc.Conf$CI.MLE[2,]), col="red", lwd=1)

plot(Hulc.Conf$Band.CI.Init$lwr, add=TRUE, col="blue")
plot(Hulc.Conf$Band.CI.Init$upr, add=TRUE, col="blue")
plot(Hulc.Conf$Band.CI.MLE$lwr, add=TRUE, col="red", lwd=1)
plot(Hulc.Conf$Band.CI.MLE$upr, add=TRUE, col="red", lwd=1)

