rm(list = ls())   # Remove everything from the Global environment
library(rmutil)   # For generating from pareto distribution
library(ggplot2)
library(GoFKernel)
library(transport)  # For calculating Wasserstein distance between two distribution functions.

source("EbayFunctions_Ver2.R")
source("HulC.R")

raw.unorganized.data <- read.csv("Xbox_7day_auctions.csv", header = T)
raw.data.list <- Xbox.data.extract.function(raw.unorganized.data)


#---------------------------------------------------------------------------------------------------
## Comparing G_lambda(F.init), G_lambda(F.mle) with ecdf(selling prices)

# Training dataset
proportion.of.train <- 1/2
train.idx <- sample.int(raw.data.list$N.auction, ceiling(proportion.of.train * raw.data.list$N.auction))

train.data.list <- raw.data.list
train.data.list$Raw.biddata.list <- train.data.list$Raw.biddata.list[sort(train.idx)]
train.data.list$N.auction <- length(train.idx)
train.data.list$reserve.price <- train.data.list$reserve.price[sort(train.idx)]
train.processed.data <- Data.Gen.2nd.Price.Processed(train.data.list)

reserve.price.Cutoff.train <- ceiling(quantile(train.processed.data$reserve.price, probs = 0.25))

Init.Values.train <- initialization.2nd.price(data = train.processed.data, reserve.price.Cutoff.train)
theta.init.train <- (1-Init.Values.train$F.y)/ c(1,(1-Init.Values.train$F.y)[-length(Init.Values.train$F.y)])
MLE.train <- MLE.2ndprice.init(train.processed.data, lambda.in = Init.Values.train$lambda, theta.in = theta.init.train, tol = 1e-5)


## Comparing the ecdf of selling prices with the function: G_lambda(F_0(.))
G.lambda.fun <- function(Fx){
  lambda.tau <- Init.Values.train$lambda * train.processed.data$auction.window
  numerator <- ((exp(lambda.tau*Fx) - 1)*((lambda.tau*(1 - Fx)) + 1)) - (lambda.tau*Fx)
  denominator <- exp(lambda.tau) - lambda.tau - 1
  return(numerator/denominator)
}

G.lambda.F.init.train.vec <- sapply(Init.Values.train$F.y, FUN = G.lambda.fun)
G.lambda.F.MLE.train.vec <- sapply(MLE.train$F.y, FUN = G.lambda.fun)

## Testing dataset
test.data.list <- raw.data.list
test.data.list$Raw.biddata.list <- test.data.list$Raw.biddata.list[-sort(train.idx)]
test.data.list$N.auction <- raw.data.list$N.auction - length(train.idx)
test.data.list$reserve.price <- test.data.list$reserve.price[-sort(train.idx)]
test.processed.data <- Data.Gen.2nd.Price.Processed(test.data.list)

# reserve.price.Cutoff.test <- ceiling(quantile(test.processed.data$reserve.price, probs = 0.25))
# Init.Values.test <- initialization.2nd.price(data = test.processed.data, reserve.price.Cutoff.test)
F.s.test.function <- ecdf(test.processed.data$selling.price)

## Rohit Code
plot(F.s.test.function)
lines(Init.Values.train$F.x, G.lambda.F.init.train.vec, col="blue")
lines(MLE.train$F.x, G.lambda.F.MLE.train.vec, col="red", type="l")
# plot(MLE.train$F.x, MLE.train$F.y, col="red")
# plot(Init.Values.train$F.x, MLE.train$F.y, col="red")


## Performing Kolmogorov-Smirnov statistic

# K.S.init <- ks.test(x = F.s.test.data.vec, y = G.lambda.F.init.train.At.test.data.vec)
# print(K.S.init)
# print(paste("K.S.distance.init = ", K.S.init$statistic))
# K.S.MLE <- ks.test(x = F.s.test.data.vec, y = G.lambda.F.MLE.train.At.test.data.vec)
# print(K.S.MLE)
# print(paste("K.S.distance.MLE = ", K.S.MLE$statistic))

## finding L2-norm
# library(wavethresh)
# l2norm.init <- wavethresh::l2norm(F.s.test.data.vec, G.lambda.F.init.train.At.test.data.vec)
# l2norm.MLE <- wavethresh::l2norm(F.s.test.data.vec, G.lambda.F.MLE.train.At.test.data.vec)
# print(paste("l2norm.init = ", l2norm.init))
# print(paste("l2norm.MLE = ", l2norm.MLE))

## finding L-infinity norm
# linfnorm.init <- wavethresh::linfnorm(F.s.test.data.vec, G.lambda.F.init.train.At.test.data.vec)
# linfnorm.MLE <- wavethresh::linfnorm(F.s.test.data.vec, G.lambda.F.MLE.train.At.test.data.vec)
# print(paste("L-infinity norm.init = ", linfnorm.init))
# print(paste("L-infinity norm.MLE = ", linfnorm.MLE))




#---------------------------------------------------------------------------------------------------
## Comparing the three distribution functions: F.init.train, F.mle.train, and F.mle.test.

nrep <- 1000
proportion.of.train.vec <- c(1/2, 2/3)

wasserstein.distance.mat <- matrix(0, nrow = length(proportion.of.train.vec), ncol = 3)
colnames(wasserstein.distance.mat) <- c("proportion.of.training.data", "distance-F.init.train-To-F.mle.test",
                                        "distance-F.mle.train-To-F.mle.test")
wasserstein.distance.mat[,1] <- proportion.of.train.vec

# par(mfrow = c(2,1))

for(pp in 1:length(proportion.of.train.vec)){
  dist.mat <- matrix(0, nrow = nrep, ncol = 2)
  
  for(ii in 1:nrep){
    train.idx <- sample.int(raw.data.list$N.auction, ceiling(proportion.of.train.vec[pp] * raw.data.list$N.auction))
    train.data.list <- raw.data.list
    train.data.list$Raw.biddata.list <- train.data.list$Raw.biddata.list[sort(train.idx)]
    train.data.list$N.auction <- length(train.idx)
    train.data.list$reserve.price <- train.data.list$reserve.price[sort(train.idx)]
    train.processed.data <- Data.Gen.2nd.Price.Processed(train.data.list)
    
    reserve.price.Cutoff.train <- ceiling(quantile(train.processed.data$reserve.price, probs = 0.25))
    Init.Values.train <- initialization.2nd.price(data = train.processed.data, reserve.price.Cutoff.train)
    theta.init.train <- (1-Init.Values.train$F.y)/ c(1,(1-Init.Values.train$F.y)[-length(Init.Values.train$F.y)])
    MLE.train <- MLE.2ndprice.init(train.processed.data, lambda.in = Init.Values.train$lambda, theta.in = theta.init.train, tol = 1e-5)
    
    test.data.list <- raw.data.list
    test.data.list$Raw.biddata.list <- test.data.list$Raw.biddata.list[-sort(train.idx)]
    test.data.list$N.auction <- raw.data.list$N.auction - length(train.idx)
    test.data.list$reserve.price <- test.data.list$reserve.price[-sort(train.idx)]
    test.processed.data <- Data.Gen.2nd.Price.Processed(test.data.list)
    
    reserve.price.Cutoff.test <- ceiling(quantile(test.processed.data$reserve.price, probs = 0.25))
    Init.Values.test <- initialization.2nd.price(data = test.processed.data, reserve.price.Cutoff.test)
    theta.init.test <- (1-Init.Values.test$F.y)/ c(1,(1-Init.Values.test$F.y)[-length(Init.Values.test$F.y)])
    MLE.test <- MLE.2ndprice.init(test.processed.data, lambda.in = Init.Values.test$lambda, theta.in = theta.init.test, tol = 1e-5)
    
    # Plotting of F.init.train vs. F.mle.train vs. F.mle.test
    # plot(Init.Values.train$F.x, Init.Values.train$F.y, type = "l", col = "blue",
    #      main = paste("F.init.train [blue] vs. F.mle.train [red] vs. F.mle.test [orange], 
    #                   for proportion of training = ", round(proportion.of.train.vec[pp], digits = 4)),
    #      xlab = "x", ylab = "F(x)")
    # lines(MLE.train$F.x, MLE.train$F.y, col = "red")
    # lines(MLE.test$F.x, MLE.test$F.y, col = "orange")
    
    
    ## Calculating the Wasserstein distance between two distribution functions
    diff.F.init.train.y <- c(Init.Values.train$F.y[1], diff(Init.Values.train$F.y))
    weights.init.train <- diff.F.init.train.y[which(diff.F.init.train.y != 0)]
    x.init.train <- Init.Values.train$F.x[which(diff.F.init.train.y != 0)]
    
    diff.F.MLE.train.y <- c(MLE.train$F.y[1], diff(MLE.train$F.y))
    weights.MLE.train <- diff.F.MLE.train.y[which(diff.F.MLE.train.y != 0)]
    x.MLE.train <- MLE.train$F.x[which(diff.F.MLE.train.y != 0)]
    
    diff.F.MLE.test.y <- c(MLE.test$F.y[1], diff(MLE.test$F.y))
    weights.MLE.test <- diff.F.MLE.test.y[which(diff.F.MLE.test.y != 0)]
    x.MLE.test <- MLE.test$F.x[which(diff.F.MLE.test.y != 0)]
    
    dist.mat[ii,1] <- wasserstein1d(a = x.init.train, b = x.MLE.train, p = 1,
                                    wa = weights.init.train,
                                    wb = weights.MLE.train )
    
    dist.mat[ii,2] <- wasserstein1d(a = x.MLE.train, b = x.MLE.test, p = 1,
                                    wa = weights.MLE.train,
                                    wb = weights.MLE.test )
  }
  
  # boxplot.matrix(dist.mat, use.cols = T, main = paste("Boxplots with", round(proportion.of.train.vec[pp], digits = 4),
  #                                                     "proportion of training data"),
  #                xlab = "Boxplots of W(F.init.train, F.mle.test) and W(F.mle.train, F.mle.test) values",
  #                ylab = "Wasserstein distances")
  
  wasserstein.distance.mat[pp,2] <- mean(dist.mat[,1])
  wasserstein.distance.mat[pp,3] <- mean(dist.mat[,2])
}

# par(mfrow = c(1,1))
print(wasserstein.distance.mat)


