## This file contains codes for the model performance evaluation of
## both F_{initial} and F_{MLE} (estimates of F, obtained from our
## proposed non-parametric methodology) on the Xbox 7 day auctions dataset 
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
## Read and load the dataset. Then, extract the necessary informations (in raw format)
## from the dataset.
#------------------------------------------------------------------------------------
raw.unorganized.data <- read.csv("Xbox_7day_auctions.csv", header = T)
raw.data.list <- Xbox.data.extract.function(raw.unorganized.data)


#------------------------------------------------------------------------------------
## Train the model with the training dataset.
## The training dataset is obtained by splitting the whole Xbox 7 day 
## auctions dataset (e.g., 1:1 and 2:1 split for the training:test split)
## and taking that corresponding proportion as our training data.
## We then compare G_lambda(F.init(.)), G_lambda(F.mle(.)) with the empirical CDF
## of final selling prices, on the training dataset.
#------------------------------------------------------------------------------------
## Training dataset
proportion.of.train <- 1/2
train.idx <- sample.int(raw.data.list$N.auction, 
                        ceiling(proportion.of.train * raw.data.list$N.auction))
train.data.list <- raw.data.list
train.data.list$Raw.biddata.list <- train.data.list$Raw.biddata.list[sort(train.idx)]
train.data.list$N.auction <- length(train.idx)
train.data.list$reserve.price <- train.data.list$reserve.price[sort(train.idx)]
train.processed.data <- Data.Gen.2nd.Price.Processed(train.data.list)

## Implementing our proposed non-parametric methodology on the training dataset
## to get the estimates F_{initial} and F_{MLE} of true F.
reserve.price.Cutoff.train <- ceiling(quantile(train.processed.data$reserve.price, 
                                               probs = 0.25) )
Init.Values.train <- initialization.2nd.price(data = train.processed.data, 
                                              reserve.price.Cutoff.train)
theta.init.train <- (1-Init.Values.train$F.y)/ c(1,(1-Init.Values.train$F.y)[-length(Init.Values.train$F.y)])
MLE.train <- MLE.2ndprice.init(train.processed.data, 
                               lambda.in = Init.Values.train$lambda, 
                               theta.in = theta.init.train, tol = 1e-5)

## G_lambda(F(.))
G.lambda.fun <- function(Fx){
  lambda.tau <- Init.Values.train$lambda * train.processed.data$auction.window
  numerator <- ((exp(lambda.tau*Fx) - 1)*((lambda.tau*(1 - Fx)) + 1)) - (lambda.tau*Fx)
  denominator <- exp(lambda.tau) - lambda.tau - 1
  return(numerator/denominator)
}

## Vector of values of G_lambda(F.init(.)) evaluated at training dataset.
G.lambda.F.init.train.vec <- sapply(Init.Values.train$F.y, FUN = G.lambda.fun)
## Vector of values of G_lambda(F.mle(.)) evaluated at training dataset.
G.lambda.F.MLE.train.vec <- sapply(MLE.train$F.y, FUN = G.lambda.fun)

## ecdf of final selling prices of all the auctions in the training dataset.
F.s.train.function <- ecdf(train.processed.data$selling.price)

## Plot of G_lambda(F.init(.)) vs. G_lambda(F.mle(.)) vs. the empirical CDF
## of final selling prices, evaluated on the training dataset.
plot(F.s.train.function)
lines(Init.Values.train$F.x, G.lambda.F.init.train.vec, col="blue")
lines(MLE.train$F.x, G.lambda.F.MLE.train.vec, col="red", type="l")


#------------------------------------------------------------------------------------
## Evaluate the model performance.
## We construct both F.init.train and F.mle.train based on the training dataset, and
## compare each of them with F.mle.test, constructed based on the testing data.

## F.init.train represents the estimate F_{initial} constructed based on the training dataset.
## F.mle.train represents the estimate F_{MLE} constructed based on the training dataset.
## F.mle.test represents the estimate F_{MLE} constructed based on the testing dataset.

## So, basically we compare the three distribution functions: F.init.train, F.mle.train,
## and F.mle.test.
## We also find the Wasserstein distance between F.mle.test and both of F.init.train,
## F.mle.test, averaged over 1000 replications of the random splits with same proportion.
#------------------------------------------------------------------------------------
nrep <- 1000  # Number of replications.
proportion.of.train.vec <- c(1/2, 2/3)  # Vector of split proportions.
wasserstein.distance.mat <- matrix(0, nrow = length(proportion.of.train.vec), ncol = 3)
colnames(wasserstein.distance.mat) <- c("proportion.of.training.data", 
                                        "distance-F.init.train-To-F.mle.test",
                                        "distance-F.mle.train-To-F.mle.test")
wasserstein.distance.mat[,1] <- proportion.of.train.vec

par(mfrow = c(2,1))

for(pp in 1:length(proportion.of.train.vec)){
  dist.mat <- matrix(0, nrow = nrep, ncol = 2)
  for(ii in 1:nrep){
    ## Training dataset
    train.idx <- sample.int(raw.data.list$N.auction, 
                            ceiling(proportion.of.train.vec[pp] * raw.data.list$N.auction))
    train.data.list <- raw.data.list
    train.data.list$Raw.biddata.list <- train.data.list$Raw.biddata.list[sort(train.idx)]
    train.data.list$N.auction <- length(train.idx)
    train.data.list$reserve.price <- train.data.list$reserve.price[sort(train.idx)]
    train.processed.data <- Data.Gen.2nd.Price.Processed(train.data.list)
    
    ## Implementing our proposed non-parametric methodology on the training dataset
    ## to get the estimates F_{initial} and F_{MLE} based on training data.
    reserve.price.Cutoff.train <- ceiling(quantile(train.processed.data$reserve.price, 
                                                   probs = 0.25))
    Init.Values.train <- initialization.2nd.price(data = train.processed.data, 
                                                  reserve.price.Cutoff.train)
    theta.init.train <- (1-Init.Values.train$F.y)/ c(1,(1-Init.Values.train$F.y)[-length(Init.Values.train$F.y)])
    MLE.train <- MLE.2ndprice.init(train.processed.data, 
                                   lambda.in = Init.Values.train$lambda, 
                                   theta.in = theta.init.train, tol = 1e-5)
    
    ## Testing dataset
    test.data.list <- raw.data.list
    test.data.list$Raw.biddata.list <- test.data.list$Raw.biddata.list[-sort(train.idx)]
    test.data.list$N.auction <- raw.data.list$N.auction - length(train.idx)
    test.data.list$reserve.price <- test.data.list$reserve.price[-sort(train.idx)]
    test.processed.data <- Data.Gen.2nd.Price.Processed(test.data.list)
    
    ## Implementing our proposed non-parametric methodology on the testing dataset
    ## to get the estimates F_{initial} and F_{MLE} based on testing data.
    reserve.price.Cutoff.test <- ceiling(quantile(test.processed.data$reserve.price, 
                                                  probs = 0.25))
    Init.Values.test <- initialization.2nd.price(data = test.processed.data, 
                                                 reserve.price.Cutoff.test)
    theta.init.test <- (1-Init.Values.test$F.y)/ c(1,(1-Init.Values.test$F.y)[-length(Init.Values.test$F.y)])
    MLE.test <- MLE.2ndprice.init(test.processed.data, 
                                  lambda.in = Init.Values.test$lambda, 
                                  theta.in = theta.init.test, tol = 1e-5)
    
    ## Plotting of F.init.train vs. F.mle.train vs. F.mle.test
    plot(Init.Values.train$F.x, Init.Values.train$F.y, type = "l", col = "blue",
         main = paste("F.init.train [blue] vs. F.mle.train [red] vs. F.mle.test [orange],
                      for proportion of training = ", 
                      round(proportion.of.train.vec[pp], digits = 4) ),
         xlab = "x", ylab = "F(x)")
    lines(MLE.train$F.x, MLE.train$F.y, col = "red")
    lines(MLE.test$F.x, MLE.test$F.y, col = "orange")
    
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
  
  wasserstein.distance.mat[pp,2] <- mean(dist.mat[,1])
  wasserstein.distance.mat[pp,3] <- mean(dist.mat[,2])
}

par(mfrow = c(1,1))

## Matrix/Table containing the Wasserstein distances between F.mle.test and both of F.init.train,
## F.mle.test, averaged over 1000 replications of the random splits with same proportion.
print(wasserstein.distance.mat)

