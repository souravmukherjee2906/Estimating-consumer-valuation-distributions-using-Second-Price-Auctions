rm(list = ls())   # Remove everything from the Global environment
library(rmutil)   # For generating from pareto distribution
library(ggplot2)
library(GoFKernel)
library(transport)  # For calculating Wasserstein distance between two distribution functions.


source("EbayFunctions_Ver2.R")
source("HulC.R")

raw.unorganized.data <- read.csv("Xbox_7day_auctions.csv", header = T)
# raw.unorganized.data$bid <- raw.unorganized.data$bid + runif(length(raw.unorganized.data$bid), min = 0, max = 0.01)
raw.data.list <- Xbox.data.extract.function(raw.unorganized.data)
  
processed.data <- Data.Gen.2nd.Price.Processed(raw.data.list)

x.med.vec <- seq(0, max(processed.data$pooled.data[,1]), by = max(processed.data$pooled.data[,1])/20)
# Init.cov <- rep(0, length(x.med.vec))
# MLE.Cov <- rep(0, length(x.med.vec))


reserve.price.Cutoff <- ceiling(quantile(processed.data$reserve.price, probs = 0.35))

Init.Values <- initialization.2nd.price(data = processed.data, reserve.price.Cutoff)
theta.init <- (1-Init.Values$F.y)/ c(1,(1-Init.Values$F.y)[-length(Init.Values$F.y)])
MLE <- MLE.2ndprice.init(processed.data, lambda.in = Init.Values$lambda, theta.in = theta.init, tol = 1e-5)
Hulc.Conf <- HulC1d_Ebay(data.in = raw.data.list, eval.points = x.med.vec, alpha=.1, MLE.Cutoff = reserve.price.Cutoff)


## Wasserstein distance between F.mle and F.init
diff.F.init.y <- c(Init.Values$F.y[1], diff(Init.Values$F.y))
weights.init <- diff.F.init.y[which(diff.F.init.y != 0)]
x.init <- Init.Values$F.x[which(diff.F.init.y != 0)]

diff.F.MLE.y <- c(MLE$F.y[1], diff(MLE$F.y))
weights.MLE <- diff.F.MLE.y[which(diff.F.MLE.y != 0)]
x.MLE <- MLE$F.x[which(diff.F.MLE.y != 0)]

dist_F.init_F.mle <- wasserstein1d(a = x.init, b = x.MLE, p = 1, wa = weights.init, wb =  weights.MLE)

print(paste("Wasserstein distance between F.init and F.mle for the whole Xbox dataset =", dist_F.init_F.mle))


## Plotting the F.mle and F.init for the Xbox Auction dataset
plot(Init.Values$F.x, Init.Values$F.y, type="l", col = "blue", main = "Plot of F_init (blue) Vs F_MLE (red)",
     xlab = "x", ylab = "F(x)")
lines(MLE$F.x, MLE$F.y, col = "red")

# HulC Confidence band for the whole dataset
points(x.med.vec, as.vector( Hulc.Conf$CI.Init[1,]), col="blue", lwd=1)
points(x.med.vec, as.vector( Hulc.Conf$CI.Init[2,]), col="blue", lwd=1)
points(x.med.vec, as.vector( Hulc.Conf$CI.MLE[1,]), col="red", lwd=1)
points(x.med.vec, as.vector( Hulc.Conf$CI.MLE[2,]), col="red", lwd=1)

plot(Hulc.Conf$Band.CI.Init$lwr, add=TRUE, col="blue")
plot(Hulc.Conf$Band.CI.Init$upr, add=TRUE, col="blue")
plot(Hulc.Conf$Band.CI.MLE$lwr, add=TRUE, col="red", lwd=1)
plot(Hulc.Conf$Band.CI.MLE$upr, add=TRUE, col="red", lwd=1)

# Confidence band based on each partition of the dataset
for (jj in 1:Hulc.Conf$B){
  lines(x.med.vec, Hulc.Conf$extra.info[[1]][,jj], col="red", lty=2)
}


# Plotting the F.f and F.s at several evaluation points
evalpt <- seq(0, 400, by =.1)
lines(evalpt, Init.Values$F.f(evalpt), type="l", col = "green")

plot.vec <- rep(0,length(evalpt) )  
for (ii in 1:length(evalpt)){
	plot.vec[ii]<- Init.Values$F.s(evalpt[ii])
}
lines(evalpt, plot.vec, type="l", col = "red")  





  
## Plot of the Xbox Auction dataset

# raw.pooled.data.time.sort <- raw.pooled.data[order(raw.pooled.data$bidtime),]
# plot(raw.pooled.data.time.sort$bidtime, raw.pooled.data.time.sort$bid, type = "h", xlab = "BidTime", ylab = "Bid", main = "Xbox 7 day Auction")
  