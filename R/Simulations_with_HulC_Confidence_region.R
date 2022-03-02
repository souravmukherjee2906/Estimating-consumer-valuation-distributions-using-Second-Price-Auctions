rm(list = ls())   # Remove everything from the Global environment
library(rmutil)   # For generating from pareto distribution
library(ggplot2)
library(GoFKernel)  # For the inverse function of a CDF
library(transport)  # For calculating Wasserstein distance between two distribution functions.

source("EbayFunctions_Ver2.R")
source("HulC.R")

lambda.true <- 1
# repetition_number <- 100        # replication numbers of each such event i.e. selling an identical stuff on ebay with multiple independent auctions.
auction.window <- 100
N.auction <- 200


## Different distributions and the corresponding choice of parameters and reserve prices are listed below:


# method.in <- "unif"
# para1 <- 1
# para2 <- 20
# x.med.vec <- seq(0.1, 20, by= 0.995)
# reserve.price <- runif(N.auction, 0.1, 3)
# reserve.price.Cutoff <- 1


# method.in <- "pareto"
# para1 <- 3
# para2 <- 100
# x.med.vec <- seq(0, 20, by= 1)
# reserve.price <- runif(N.auction, 0.001, 0.1)
# reserve.price.Cutoff <- 0.05


method.in <- "gamma"
para1 <- 10
para2 <- 2
x.med.vec <- seq(0, 12, by= 12/20)
reserve.price <- runif(N.auction, 0.1, 3)
reserve.price.Cutoff <- 1


# method.in <- "beta"
# para1 <- 2
# para2 <- 2
# x.med.vec <- seq(0, 1, by= 1/20)
# reserve.price <- runif(N.auction, 0.001, 0.2)
# reserve.price.Cutoff <- 0.05


# method.in <- "piecewise_unif"
# para1 <- 2
# para2 <- 4
# x.med.vec <- seq(0, 4, by = 4/20)
# reserve.price <- runif(N.auction, 0.1, 1.5)
# reserve.price.Cutoff <- 1


# N.auction <- N.auction.vec[jj]

raw.data.list <- Data.Gen.2nd.Price.Raw(N.auction, auction.window, lambda.true, reserve.price, method = method.in, para1, para2)
data <- Data.Gen.2nd.Price.Processed(raw.data.list)
LL <- length(data$pooled.observed.bids)
Init.Values <- initialization.2nd.price(data, reserve.price.Cutoff)
theta.init <- (1-Init.Values$F.y)/ c(1,(1-Init.Values$F.y)[-length(Init.Values$F.y)])

MLE.old <- MLE.2ndprice(data, lambda.in = Init.Values$lambda, theta.in = theta.init, tol = 1e-5)
MLE <- MLE.2ndprice.init(data, lambda.in = Init.Values$lambda, theta.in = theta.init, tol = 1e-5)





Hulc.Conf <- HulC1d_Ebay(data.in = raw.data.list, eval.points = x.med.vec, alpha=.1, MLE.Cutoff = reserve.price.Cutoff)


# Init.cov <- Init.cov+ (Hulc.Conf[[ii]]$CI.Init[1,] < pbeta(x.med.vec, shape1 = para1, shape2 = para2))*
# (Hulc.Conf[[ii]]$CI.Init[2,] > pbeta(x.med.vec, shape1 = para1, shape2 = para2))

# MLE.Cov <- MLE.Cov+  (Hulc.Conf[[ii]]$CI.MLE[1,] < pbeta(x.med.vec, shape1 = para1, shape2 = para2))*
# (Hulc.Conf[[ii]]$CI.MLE[2,] > pbeta(x.med.vec, shape1 = para1, shape2 = para2))



# Init.cov  <- Init.cov  /repetition_number
# MLE.Cov<- MLE.Cov/repetition_number
# save.image(paste("(U)N.auction=" , N.auction, "__repetition_number=", repetition_number, "__breaks=20__p",method.in,"__Hulc.RData"))
# lines(eval.points, as.vector(testing$CI.MLE[1,]), col="blue", lwd=2)
# lines(eval.points, as.vector(testing$CI.MLE[2,]), col="blue", lwd=2)



# rm(list=ls())
# library(ggplot2)
# filename<- "(U)N.auction=200__repetition_number=100__breaks=20__ppiecewise_unif__Hulc"
# load(paste0(filename, ".RData"))

#---------------------------------------------------------------------------------------------------
## Plotting the confidence interval and the MLEs
# par(mfrow = c(1, 1))


plot(Init.Values$F.x, Init.Values$F.y, type="l", col = "blue",
     main = paste("Plot of true F (black) Vs. F_init(blue) Vs. F_MLE_OLD (green), N.auction =", N.auction, ", True F =", data$true.F),
     xlab = "x", ylab = "F(x)")


plot(Init.Values$F.x, Init.Values$F.y, type="l", col = "blue",
     main = paste("Plot of true F (black) Vs. F_init(blue) Vs. F_MLE (red), N.auction =", N.auction, ", True F =", data$true.F),
     xlab = "x", ylab = "F(x)")

lines(MLE.old$F.x, MLE.old$F.y, type="l", col = "green")
lines(MLE$F.x, MLE$F.y, type="l", col = "red")

# lines(data$pooled.data[,1], pbeta(data$pooled.data[,1], shape1 = para1, shape2 = para2), col = "black", type = "l")
lines(data$pooled.data[,1], pgamma(data$pooled.data[,1], shape = para1, rate = para2), col = "black", type = "l")
# lines(data$pooled.data[,1], punif(data$pooled.data[,1], min = para1, max = para2), col = "black", type = "l")
# lines(data$pooled.data[,1], ppareto(data$pooled.data[,1], m = para1, s = para2), col = "black", type = "l")
# F0.piecewise.unif_pooled.data <- sapply(data$pooled.data[,1], function(y) {piecewise.unif.cdf(y, para1 = para1, para2 = para2)} )
# lines(data$pooled.data[,1], F0.piecewise.unif_pooled.data, col = "black", type = "l")


# lines(seq(0,4, by=.001),sapply(seq(0,4, by=.001), function(y){piecewise.unif.cdf(y, para1 = para1, para2 = para2)}), col="black")




# Confidence band based on each partition of the dataset
# for (jj in 1:Hulc.Conf$B){
#   lines(x.med.vec, Hulc.Conf$extra.info[[1]][,jj], col="blue", lty=2)
# }


# HulC Confidence band for the whole dataset
points(x.med.vec, as.vector( Hulc.Conf$CI.Init[1,]), col="blue", lwd=1)
points(x.med.vec, as.vector( Hulc.Conf$CI.Init[2,]), col="blue", lwd=1)
points(x.med.vec, as.vector( Hulc.Conf$CI.MLE[1,]), col="red", lwd=1)
points(x.med.vec, as.vector( Hulc.Conf$CI.MLE[2,]), col="red", lwd=1)

plot(Hulc.Conf$Band.CI.Init$lwr, add=TRUE, col="blue")
plot(Hulc.Conf$Band.CI.Init$upr, add=TRUE, col="blue")
plot(Hulc.Conf$Band.CI.MLE$lwr, add=TRUE, col="red", lwd=1)
plot(Hulc.Conf$Band.CI.MLE$upr, add=TRUE, col="red", lwd=1)



## Wasserstein distance between True F and F.mle
# diff.F.init.y <- c(Init.Values$F.y[1], diff(Init.Values$F.y))
# weights.init <- diff.F.init.y[which(diff.F.init.y != 0)]
# x.init <- Init.Values$F.x[which(diff.F.init.y != 0)]
# 
# diff.F.MLE.y <- c(MLE$F.y[1], diff(MLE$F.y))
# weights.MLE <- diff.F.MLE.y[which(diff.F.MLE.y != 0)]
# x.MLE <- MLE$F.x[which(diff.F.MLE.y != 0)]
# 
# dist_F.init_F.mle <- wasserstein1d(a = x.init, b = x.MLE, p = 1, wa = weights.init, wb =  weights.MLE)
# 
# print(paste("Wasserstein distance between F.init and F.mle for the whole Xbox dataset =", dist_F.init_F.mle))

## Wasserstein distance between True F and F.init

