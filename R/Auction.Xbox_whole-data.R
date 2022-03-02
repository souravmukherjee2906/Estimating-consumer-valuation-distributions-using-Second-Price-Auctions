## This file contains codes for the data analysis of Xbox 7 day auctions dataset
## found in the "dataset" folder of this repository.


#------------------------------------------------------------------------------------
## Remove everything from the Global environment
#------------------------------------------------------------------------------------
rm(list = ls())


#------------------------------------------------------------------------------------
## Install the required packages into R if it's not already installed.
#------------------------------------------------------------------------------------
if(!require(GoFKernel)) install.packages("GoFKernel")
if(!require(transport)) install.packages("transport")


#------------------------------------------------------------------------------------
## Load the required packages into the R environment.
#------------------------------------------------------------------------------------
library(GoFKernel)  # For the inverse function of a CDF.
library(transport)  # For calculating Wasserstein distance between two distribution functions.


#------------------------------------------------------------------------------------
## Source all the functions from "Functions.R" and "HulC.R"
#------------------------------------------------------------------------------------
source("EbayFunctions_Ver2.R")
source("HulC.R")


#------------------------------------------------------------------------------------
## Read and load the dataset
#------------------------------------------------------------------------------------
raw.unorganized.data <- read.csv("Xbox_7day_auctions.csv", header = T)


#------------------------------------------------------------------------------------
## Main body of codes
#------------------------------------------------------------------------------------
raw.data.list <- Xbox.data.extract.function(raw.unorganized.data)
processed.data <- Data.Gen.2nd.Price.Processed(raw.data.list)

x.med.vec <- seq(0, max(processed.data$pooled.data[,1]), 
		 by = max(processed.data$pooled.data[,1])/20)
reserve.price.Cutoff <- ceiling(quantile(processed.data$reserve.price, probs = 0.35))

Init.Values <- initialization.2nd.price(data = processed.data, reserve.price.Cutoff)
theta.init <- (1-Init.Values$F.y)/ c(1,(1-Init.Values$F.y)[-length(Init.Values$F.y)])
MLE <- MLE.2ndprice.init(processed.data, lambda.in = Init.Values$lambda, 
			 theta.in = theta.init, tol = 1e-5)
Hulc.Conf <- HulC1d_Ebay(data.in = raw.data.list, eval.points = x.med.vec, alpha=.1, 
			 MLE.Cutoff = reserve.price.Cutoff)


#------------------------------------------------------------------------------------
## Wasserstein distance between F_{MLE} and F_{initial}  (estimates of F)
#------------------------------------------------------------------------------------
diff.F.init.y <- c(Init.Values$F.y[1], diff(Init.Values$F.y))
weights.init <- diff.F.init.y[which(diff.F.init.y != 0)]
x.init <- Init.Values$F.x[which(diff.F.init.y != 0)]

diff.F.MLE.y <- c(MLE$F.y[1], diff(MLE$F.y))
weights.MLE <- diff.F.MLE.y[which(diff.F.MLE.y != 0)]
x.MLE <- MLE$F.x[which(diff.F.MLE.y != 0)]

dist_F.init_F.mle <- wasserstein1d(a = x.init, b = x.MLE, p = 1, wa = weights.init, 
				   wb =  weights.MLE)
print(paste("Wasserstein distance between F.init and F.mle for the whole Xbox dataset =", 
	    dist_F.init_F.mle))


#------------------------------------------------------------------------------------
## Plotting the F_{MLE} and F_{initial} (estimates of F) for the Xbox Auction dataset
#------------------------------------------------------------------------------------
plot(Init.Values$F.x, Init.Values$F.y, type="l", col = "blue", 
     main = "Plot of F_init (blue) Vs F_MLE (red)",
     xlab = "x", ylab = "F(x)")
lines(MLE$F.x, MLE$F.y, col = "red")


#------------------------------------------------------------------------------------
# HulC Confidence bands for both F_{MLE} and F_{initial}
#------------------------------------------------------------------------------------
points(x.med.vec, as.vector( Hulc.Conf$CI.Init[1,]), col="blue", lwd=1)
points(x.med.vec, as.vector( Hulc.Conf$CI.Init[2,]), col="blue", lwd=1)
points(x.med.vec, as.vector( Hulc.Conf$CI.MLE[1,]), col="red", lwd=1)
points(x.med.vec, as.vector( Hulc.Conf$CI.MLE[2,]), col="red", lwd=1)

plot(Hulc.Conf$Band.CI.Init$lwr, add=TRUE, col="blue")
plot(Hulc.Conf$Band.CI.Init$upr, add=TRUE, col="blue")
plot(Hulc.Conf$Band.CI.MLE$lwr, add=TRUE, col="red", lwd=1)
plot(Hulc.Conf$Band.CI.MLE$upr, add=TRUE, col="red", lwd=1)

