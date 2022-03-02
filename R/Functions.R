### source this file in other .R files in the repository. ###


#------------------------------------------------------------------------------------
##  CDF of piecewise uniform distribution.
##  Considering two uniform distributions: Unif(1, para1] and Unif(para1 + 1, para2]
##  for the corresponding piecewise uniform distribution.
#------------------------------------------------------------------------------------

piecewise.unif.cdf <- function(y, para1 = 2, para2 = 4){
  # Considering two uniform distributions: Unif(1, para1] and Unif(para1 + 1, para2] 
  if((para1<=1) || (para2 <= para1+1)) stop("This function works only when para1 > 1 
                                            and para2 > para1 + 1, since, by default,
                                            we are considering the two uniform
                                            distributions to be Unif(1, para1] and
                                            Unif(para1 + 1, para2]")
  
  if(y <= 1){
    temp <- 0
  }else if ((1 < y) && (y <= para1)){
    temp <- 0.5 * ((y-1)/(para1-1))
  }else if ((para1 < y) && (y <= (para1 + 1))){
    temp <- 0.5
  }else if (((para1 + 1) < y) && (y <= para2)){
    temp <- 0.5 + 0.5 * ((y - para1 - 1)/ (para2 - para1 - 1) )
  }else{
    temp <- 1
  }
  
  return(temp)
}


#------------------------------------------------------------------------------------
##  Inverse CDF of piecewise uniform distribution.
##  Considering two uniform distributions: Unif(1, para1] and Unif(para1 + 1, para2]
##  for the corresponding piecewise uniform distribution.
#------------------------------------------------------------------------------------

piecewise.unif.invcdf <- function(u, para1 = 2, para2 = 4){
  if((u<0) || (u >1)){
    stop("This function works only when the first argument u lies within the closed
         interval [0,1]")
  }else if(u == 0.5){
    temp <- para1
  }else{
    temp <- inverse(function(y) {piecewise.unif.cdf(y, para1, para2)},
                    lower = 1, upper = para2)(u)
  }
  return(temp)
}


#------------------------------------------------------------------------------------
##  Functions to generate the data from independent auctions  
##  where $N.auction$ is not fixed but the total time window  
##  \tau is fixed.                                            
#------------------------------------------------------------------------------------

Data.Gen.2nd.Price.Raw <- function(N.auction, auction.window, lambda.true,
                                   reserve.price,
                                   method = c("unif", "pareto", "gamma",
                                              "beta", "piecewise_unif"),
                                   para1, para2){
  arrival.count <- rep(0, N.auction)
  if(is.null(reserve.price)) reserve.price <- rep(0, N.auction)
  
  # List containing each auction data separately.
  Raw.biddata.list <- vector("list", N.auction)
  
  for (jj in 1: N.auction) {
    tts <- rexp( (2* auction.window* lambda.true), rate = lambda.true)
    # Total number of bids for the jth auction.
    arrival.count[jj] <- sum(cumsum(tts) <= auction.window)  
    totts <- cumsum(tts)[1:arrival.count[jj]]
    if(method == "unif"){
      # All the private values of the bidders for the jth process
      yy <- runif(n = arrival.count[jj], min=para1, max=para2)  
    }else if (method == "pareto"){
      # All the private values of the  bidders for the jth process
      yy <- rpareto(n = arrival.count[jj], m = para1, s = para2) 
    } else if( method == "gamma"){
      yy <- rgamma(n = arrival.count[jj], shape = para1, rate=para2) 
    }else if (method == "beta"){
      # All the private values of the bidders for the jth process
      yy <- rbeta(n = arrival.count[jj], shape1 = para1, shape2 = para2)  
    }else if (method == "piecewise_unif"){
      uu <- runif(n = arrival.count[jj])
      yy <- sapply(uu, 
                   function(u) {piecewise.unif.invcdf(u, para1 = para1, para2 = para2)})
    }
    
    Raw.biddata.matrix <- matrix(0, nrow = arrival.count[jj], ncol = 4)
    colnames(Raw.biddata.matrix) <- c("auction.index", "bid", "bidtime", "openbid")
    Raw.biddata.matrix[,1] <- rep(jj, arrival.count[jj])
    Raw.biddata.matrix[,2] <- yy
    Raw.biddata.matrix[,3] <- totts
    Raw.biddata.matrix[,4] <- rep(reserve.price[jj], arrival.count[jj])
    
    Raw.biddata.list[[jj]] <- Raw.biddata.matrix
  } # END of for(jj in 1: N.auction)
  
  Raw.data.list <- list(Raw.biddata.list = Raw.biddata.list,
                        N.auction = N.auction,
                        auction.window = auction.window,
                        reserve.price = reserve.price,
                        True.Method = method) 
  Raw.data.list$class = "SecondPriceAuction.Rawdata"
  
  return(Raw.data.list)
}



Data.Gen.2nd.Price.Processed <- function(Raw.data.list){
  if(Raw.data.list$class != "SecondPriceAuction.Rawdata") stop("This fucntion works 
                                                               only with data in the 
                                                               class 'SecondPriceAuction.Rawdata'.
                                                               See e.g., 'Data.Gen.2nd.Price.Raw'
                                                               where we generate such data.")
  
  N.auction <- Raw.data.list$N.auction
  reserve.price <- Raw.data.list$reserve.price
  auction.window <- Raw.data.list$auction.window
  
  pooled.data <- matrix(0, ncol=2, nrow=0)
  pooled.observed.bids <- NULL
  colnames(pooled.data) <- c("pooled.rec.value", "pooled.wait.time")
  # first.jumps is the vector of first change of current selling price in all the auctions.
  first.jumps <- selling.price <- rep(0, N.auction)
  
  # Vector of waiting times when the current selling price jumps for the first time
  # to a strictly larger value than the corresponding reserve price.
  T0_vec <- rep(0, N.auction)
  # Sold.Auctions.indexes[which(Sold.Auctions.indexes!=0)] is the vector of all the auction indexes where the item has been sold.
  Sold.Auctions.indexes <- rep(0, N.auction)
  obs.bids<- rep(0, N.auction)
  # List containing values and intermediate times of all the standing price
  # changes for all the auctions.
  data.list <- vector("list", N.auction)
  
  for (jj in 1: N.auction) {
    
    # Creating the list containing values and intermediate times of all the standing price
    # changes for each auction.
    if(sum(reserve.price[jj]<= Raw.data.list$Raw.biddata.list[[jj]][,2]) == 0){
      obs.bids[jj] <- 0
      xx <- NA
      selling.price[jj] <- NA
      # First jump value / first change in standing price for the jth auction.
      first.jumps[jj] <- NA
      T0_vec[jj] <- auction.window
      data.list[[jj]] <- list(res.price = reserve.price[jj], sold.ind = 0,
                              obs.bids = NULL, obs.bid.times = NULL)
      
    } else if( sum(reserve.price[jj]<= Raw.data.list$Raw.biddata.list[[jj]][,2]) == 1){
      obs.bids[jj] <- 1
      selling.price[jj] <- reserve.price[jj]
      first.jumps[jj] <- NA
      T0_vec[jj] <- auction.window
      Sold.Auctions.indexes[jj] <- jj
      data.list[[jj]] <- list(res.price = reserve.price[jj], sold.ind = 1,
                              obs.bids = NULL, obs.bid.times = NULL)
    } else if( sum(reserve.price[jj]<= Raw.data.list$Raw.biddata.list[[jj]][,2]) >= 2){
      xx <- c(reserve.price[jj], rep(0,(nrow(Raw.data.list$Raw.biddata.list[[jj]]) - 1)))
      for (ii in 2: nrow(Raw.data.list$Raw.biddata.list[[jj]])){
        if(xx[ii-1] >= Raw.data.list$Raw.biddata.list[[jj]][,2][ii]){
          xx[ii] <- xx[ii-1]
        } else { 
          obs.bids[jj] <- obs.bids[jj]+1
          temp.two.max <- sort(c(reserve.price[jj],
                                 Raw.data.list$Raw.biddata.list[[jj]][,2][1:ii]),
                               decreasing = TRUE)[-1]
          xx[ii] <- temp.two.max[1]
        }
      }
      inds <- 2:nrow(Raw.data.list$Raw.biddata.list[[jj]])
      # Indexes when selling price changes
      inds <- inds[diff(xx)>0]
      # Vector (T_0, T_0 + T_1, ...., \sum_{i=0}^{i=M-1}T_i)'
      cum.jump.times <- Raw.data.list$Raw.biddata.list[[jj]][,3][inds]
      # Vector (T_0, T_1, ...., T_M)'
      inter.arrival.time <- diff(c(0,cum.jump.times, auction.window))
      # Pooled data of reserve prices, standing prices, inter arrival times (including T_0) for all the auctions.
      pooled.data <- rbind(pooled.data, 
                           cbind(c(reserve.price[jj], xx[inds]),inter.arrival.time) )
      # Pooled data of reserve prices, standing prices, inter arrival times (including T_0) for all the auctions.
      pooled.observed.bids <- c(pooled.observed.bids, xx[inds])
      selling.price[jj] <- tail(xx[inds], n=1)
      # First jump value / first change in standing price for the jth auction.
      first.jumps[jj] <- xx[inds][1]
      T0_vec[jj] <- inter.arrival.time[1]
      Sold.Auctions.indexes[jj] <- jj
      data.list[[jj]] <- list(res.price = reserve.price[jj], sold.ind = 1,
                              obs.bids = pooled.observed.bids,
                              obs.bid.times = inter.arrival.time[-1])
    }
  } # END of for(jj in 1: N.auction)
  
  pooled.data <- pooled.data[order(pooled.data[,1]),]
  ret <- list(data.list = data.list,
              pooled.observed.bids = pooled.observed.bids, # Set of Xi's
              selling.price = selling.price,
              T0_vec = T0_vec,
              Sold.Auctions.indexes = Sold.Auctions.indexes,
              pooled.data = pooled.data,
              first.jumps = first.jumps,
              N.auction = Raw.data.list$N.auction,
              auction.window = Raw.data.list$auction.window,
              reserve.price = Raw.data.list$reserve.price,
              true.F = Raw.data.list$True.Method)
  
  ret$class = "SecondPriceAuction.Processed"
  return(ret)
}


#------------------------------------------------------------------------------------
##  Define the general version of 2nd Price likelihood function  
##  (considering all the situations)                             
#------------------------------------------------------------------------------------

General_Log_Likelihood_fn_2ndPrice <- function( data, lambda, theta){
  # Below gives the vector of auction indexes for which the item is sold.
  ss.non.zero <- data$Sold.Auctions.indexes[which(data$Sold.Auctions.indexes!=0)] 
  T0.sold.vec <- data$T0_vec[ss.non.zero]
  # Total Number of observed bids for all the auctions.
  LL <- length(data$pooled.observed.bids)  
  # Vector/Set of ranks/positions of the selling prices of all the auctions (where item is sold) in the pooled data set.
  SS.sold <- match(data$selling.price[ss.non.zero], data$pooled.data[ ,1])
  # This is vector: (uo, u1, u2,...,uL).
  RankXi <- sort(c(0, match(data$pooled.observed.bids, data$pooled.data[ ,1]))) 
  # RankXi is the vector of ranks of all the standing price changes
  # in the pooled data Z with the first element being u0 = 0.
  term_without_theta <- ((LL+length(ss.non.zero))*log(lambda)) + (LL*log(2)) + sum(log(T0.sold.vec))
  term_with_theta_1st_term <- -(lambda * sum(data$pooled.data[ ,2] * cumprod(theta))) + sum(sapply(1:(length(theta)), function(kk) {sum(SS.sold>=kk)}) * log(theta))
 
  if(LL == 0){
    term_with_theta_2nd_term <- 0
  }else{
    partial.prod.theta <- sapply(1:(length(RankXi)-1), 
                                 function(i) {prod(theta[(RankXi[i]+1) : RankXi[i+1]])})
    term_with_theta_2nd_term <- sum(log(1 - partial.prod.theta)) + sum((LL - sapply(1:(length(theta)), function(kk) {sum(RankXi < kk)})) * log(theta))
  }
  # Above if condition takes care of the case when any of the auctions do not have any observed bid i.e. LL = 0 case.
  return(term_without_theta + term_with_theta_1st_term + term_with_theta_2nd_term)
}


#------------------------------------------------------------------------------------
##  Coordinate maximization without any modifications  
#------------------------------------------------------------------------------------

cordwise.maximization <- function( data, lambda, theta){
  # Total number of standing prices from all the auctions.
  LL <- length(data$pooled.observed.bids)
  # This gives the vector of auction indexes for which the item is sold.
  ss.non.zero <- data$Sold.Auctions.indexes[which(data$Sold.Auctions.indexes!=0)]
  
  # Maximizing theta_i's one by one keeping other parameters fixed.
  LK <- length(theta)
  if(length(ss.non.zero)==0){   # if item is not sold in any of the auctions.
    theta <- rep(0, LK)
  }else{
    # This is vector: (u0, u1, u2,...,uL)
    RankXi <- sort(c(0, match(data$pooled.observed.bids, data$pooled.data[ ,1])))
    # Vector/Set of ranks/positions of the selling prices of all the auctions
    # (where item is sold) in the pooled data set.
    SS.sold <- match(data$selling.price[ss.non.zero], data$pooled.data[ ,1])   
    for (ii in 1:LK){
      if(ii==1){
        AA <- lambda * sum(data$pooled.data[ ,2] * c(1, cumprod(theta[-1])))
      }else if (ii==LK){
        AA <- lambda * prod(theta[1:(ii-1)]) * data$pooled.data[LK,2]
      }else {
        AA <- lambda * prod(theta[1:(ii-1)]) * sum(data$pooled.data[ii:LK,2] * c(1,cumprod( theta[(ii+1):LK])) )
      }
      
      if(is.nan(AA)){
        print(paste("AA", AA))
        break}
      
      if(ii<= tail(RankXi, n=1)){
        u_l1 <- RankXi[sum(RankXi < ii)]     # This is u_{l-1}
        u_l <- RankXi[sum(RankXi < ii) + 1]  # This is u_l
        BB <- sum(SS.sold>=ii) + ((1*(LL>0))*(LL - sum(RankXi < ii)) )
        if(u_l1 +1 == u_l){
          CC <- 1
        }else{
          CC <- prod(theta[-ii][(u_l1 + 1):(u_l - 1)])
        }
        
        denomenator <- 2*AA*CC
        numerator <- (AA+(BB*CC)+((1*(LL>0))*CC)) - sqrt((AA+(BB*CC)+((1*(LL>0))*CC))^2 - (4*AA*BB*CC))
        theta[ii] <- min(c(1,numerator/denomenator))
      }else{
        theta[ii] <- min(c(1,(sum(SS.sold>=ii)/AA)))
      }
    } # END of for(ii in 1:LK) loop.
  }
  return(list(lambda= lambda, theta = theta))
} # END of cordwise.maximization function.


#------------------------------------------------------------------------------------
##  Coordinate-wise maximization with first few (index.cutoff many) theta_i's  
##  defined to be theta_initial_i's in every iteration.                        
#------------------------------------------------------------------------------------

cordwise.maximization.init <- function( data, lambda, theta){
  # Total number of standing prices from all the auctions.
  LL <- length(data$pooled.observed.bids)  
  # This gives the vector of auction indexes for which the item is sold.
  ss.non.zero <- data$Sold.Auctions.indexes[which(data$Sold.Auctions.indexes!=0)]  
  
  # Maximizing theta_i's one by one keeping other parameters fixed
  LK <- length(theta)
  
  # index of the min(observed bids) value in the pooled data.
  index.cutoff <- max(which(data$pooled.data[,1] == min(data$pooled.observed.bids)) )  
  # In above, we took max() to avoid the situation where min(standing prices) occurs in 
  # more than one consecutive indexes in the pooled data.
  
  if(length(ss.non.zero)==0){   # if item is not sold in any of the auctions.
    theta <- rep(0, LK)
  }else{
    # This is vector: (u0, u1, u2,...,uL)
    RankXi <- sort(c(0, match(data$pooled.observed.bids, data$pooled.data[ ,1]))) 
    # Vector/Set of ranks/positions of the selling prices of all the auctions
    # (where item is sold) in the pooled data set.
    SS.sold <- match(data$selling.price[ss.non.zero], data$pooled.data[ ,1])   
    for (ii in 1:LK){
      if((index.cutoff < ii) && (ii < LK)){
        AA <- lambda * prod(theta[1:(ii-1)]) * sum(data$pooled.data[ii:LK,2] * c(1,cumprod( theta[(ii+1):LK])) )
      }
      if(ii==LK){
        AA <- lambda * prod(theta[1:(ii-1)]) * data$pooled.data[LK,2]
      }
      
      if(index.cutoff < ii){
      if(is.nan(AA)){
        print(paste("AA", AA))
        break}
      }
      
      if(ii <= index.cutoff){
        theta[ii] <- theta[ii]
      }else if ((index.cutoff < ii) && (ii<= tail(RankXi, n=1))){
        u_l1 <- RankXi[sum(RankXi < ii)]     # This is u_{l-1}
        u_l <- RankXi[sum(RankXi < ii) + 1]  # This is u_l
        BB <- sum(SS.sold>=ii) + ((1*(LL>0))*(LL - sum(RankXi < ii)) )
        if(u_l1 +1 == u_l){
          CC <- 1
        }else{
          CC <- prod(theta[-ii][(u_l1 + 1):(u_l - 1)])
        }
        
        denomenator <- 2*AA*CC
        numerator <- (AA+(BB*CC)+((1*(LL>0))*CC)) - sqrt((AA+(BB*CC)+((1*(LL>0))*CC))^2 - (4*AA*BB*CC))
        theta[ii] <- min(c(1,numerator/denomenator)) 
      }else{
        theta[ii] <- min(c(1,(sum(SS.sold>=ii)/AA)))
      }
    } # END of for(ii in 1:LK) loop.
  }
  return(list(lambda= lambda, theta = theta))
} # END of cordwise.maximization.init function.


#------------------------------------------------------------------------------------
## Function for finding Non-parametric MLE of the consumer valuation distribution 
## function F (without any modification to theta_MLE).                           
#------------------------------------------------------------------------------------

MLE.2ndprice <- function(data, lambda.in, theta.in, tol){
  Lik.in <- -1e100
  Lik.out <- General_Log_Likelihood_fn_2ndPrice(data, lambda.in, theta.in)
  Lik.path<- NULL

  while(Lik.out > Lik.in + tol){
    Lik.in <- Lik.out
    Lik.path <- c(Lik.path, Lik.in )
    par.out <- cordwise.maximization(data, lambda.in, theta.in)
    # lambda.in <- par.out$lambda
    theta.in <- par.out$theta
    Lik.out <- General_Log_Likelihood_fn_2ndPrice(data, lambda.in, theta.in)
    if(is.na(Lik.out) ){
      print(paste("in", Lik.in, "out", Lik.out))
      Lik.out= -1e100
    }
    # print(paste(theta.in[1:2]))
  }
 G.mle <- cumprod(theta.in)
 F.mle <- 1 - G.mle
 return(list(Lik = Lik.out, theta.mle = theta.in, lambda = lambda.in,
             F.x = data$pooled.data[,1], F.y = F.mle, Lik.path = Lik.path))
}



#------------------------------------------------------------------------------------
##  Function for finding the values (at the observed data points) of the 
##  Non-parametric MLE of the consumer valuation distribution function F,
##  with the modifications that the with first few (index.cutoff many)    
##  theta_i's defined to be theta_initial_i's.
#------------------------------------------------------------------------------------

MLE.2ndprice.init <- function(data, lambda.in, theta.in, tol){
  Lik.in <- -1e100
  Lik.out <- General_Log_Likelihood_fn_2ndPrice(data, lambda.in, theta.in)
  Lik.path <- NULL
  
  while(Lik.out > Lik.in + tol){
    Lik.in <- Lik.out
    Lik.path <- c(Lik.path, Lik.in )
    par.out <- cordwise.maximization.init(data, lambda.in, theta.in)
    theta.in <- par.out$theta
    Lik.out <- General_Log_Likelihood_fn_2ndPrice(data, lambda.in, theta.in)
    
    if(is.na(Lik.out) ){
      print(paste("in", Lik.in, "out", Lik.out))
      Lik.out <- -1e100
    }
  }
  
  G.mle <- cumprod(theta.in)
  F.mle <- 1 - G.mle
  return(list(Lik = Lik.out, theta.mle = theta.in, lambda = lambda.in,
              F.x = data$pooled.data[,1], F.y = F.mle, Lik.path = Lik.path))
}


#------------------------------------------------------------------------------------
##  Function to compute Initial values of Lambda and F for the coordinate ascent  
##  algorithm.                                                                     
#------------------------------------------------------------------------------------

initialization.2nd.price <- function(data, reserve.price.Cutoff){
  # Lambda initialization
  temp.fun <- function(lambda){
    N.ran <- rpois(10^4, lambda * data$auction.window)
    k.lam <- 2* mean(log(N.ran[N.ran>=2]))+ 2*0.5772156649-2
    mu.n <- length(data$pooled.observed.bids)/length(data$selling.price)
    return(abs(k.lam- mu.n))
  }
  lambda.initial <- optimize(temp.fun, interval = c(0, 10))$minimum
  
  # F initialization
  G.lambda.fun <- function(Fx){
    lambda.tau <- lambda.initial * data$auction.window
    numerator <- ((exp(lambda.tau*Fx) - 1)*((lambda.tau*(1 - Fx)) + 1)) - (lambda.tau*Fx)
    denominator <- exp(lambda.tau) - lambda.tau - 1
    return(numerator/denominator)
  }
  
  selling.price <- data$selling.price[intersect(which(data$reserve.price < reserve.price.Cutoff), 
                                                which(is.na(data$selling.price)==FALSE))]
  # Removing the NA terms and those auction's first jump values whose reserve prices are comparatively higher.
  # We are taking those auctions whose reserve prices are less than an user-specified reserve.price.Cutoff.
  
  Selling.price.cdf.inverse <- inverse(G.lambda.fun, lower = 0, upper = 1)
  G.n <- ecdf(selling.price)
  F.s <- function(x){
    return(Selling.price.cdf.inverse(G.n(x)))
  }
  
  first.jumps.vec <- data$first.jumps[intersect(which(data$reserve.price < reserve.price.Cutoff),
                                                which(is.na(data$first.jumps)==FALSE))]
  # Removing the NA terms and those auction's first jump values whose reserve prices are comparatively higher.
  # We are taking those auctions whose reserve prices are less than an user-specified reserve.price.Cutoff.
  
  G.n.first.jumps <- ecdf(first.jumps.vec)
  F.f <- function(x){
    return(1 - sqrt(1- G.n.first.jumps(x)))
  }

  m1 <- min(selling.price)
  m2 <- max(first.jumps.vec)
  seq.vec <- seq(0, min(m1,m2), by = 0.0001)
  Fc <- sapply(seq.vec, F.f) - F.s(m1)
  c <- seq.vec[sum(1*(Fc<=0))]
  
  F.initial <- function(x){
    first.term <- (F.f(x) * (1*(x<=c))) + (F.s(x) * (1*(x>m1)))
    second.term <- (F.f(c) + (((F.s(m1) - F.f(c))/(m1 - c)) * (x - c)) ) * ((1*(x<=m1)) - (1*(x<=c)))
    return(first.term + second.term)
  }  
  Fis <- sapply(data$pooled.data[ ,1], F.initial )
  x.temp<- data$pooled.data[ ,1]
  x.jump <- c(0, x.temp[diff(c(0,Fis))>0])
  x.jump[length(x.jump)] <- 1.1*max(x.temp)
  Fis.jump <- c(0, Fis[diff(c(0,Fis))>0])
  Fis.linear <- approx (x.jump, y = Fis.jump, x.temp, method = "linear", yleft=0, yright=.99999)
  
  return(list(lambda = lambda.initial, F.x = Fis.linear$x, F.y = Fis.linear$y,
              F.f= F.f, F.s= F.s ))
}


#------------------------------------------------------------------------------------
##  Function to extract raw data from a real Xbox Auctions data set.  
#------------------------------------------------------------------------------------

Xbox.data.extract.function <- function(raw.unorganized.data){
  colnames(raw.unorganized.data) <- c("auctionid", "bid", "bidtime", "bidder",
                                      "bidderrate", "openbid", "price")
  
  # Adding small random noise to all the bids to avoid equality among any two bid values.
  raw.unorganized.data$bid <- raw.unorganized.data$bid + runif(length(raw.unorganized.data$bid),
                                                               min = 0, max = 0.01)
  auction.window <- ceiling(max(raw.unorganized.data$bidtime))
  unique.auction.data <- unique(cbind(raw.unorganized.data$auctionid,
                                      raw.unorganized.data$openbid))
  N.auction <- nrow(unique.auction.data)
  # Vector of total number of bids in all the auctions.
  arrival.count<- rep(0, N.auction)  
  
  reserve.price <- unique.auction.data[,2]
  if(is.null(reserve.price)){ reserve.price <- rep(0, N.auction)}
  
  # List containing each auction data separately.
  Raw.biddata.list <- vector("list", N.auction) 
  
  for(jj in 1:N.auction){
    totts <- raw.unorganized.data$bidtime[which(raw.unorganized.data$auctionid == unique.auction.data[,1][jj])]
    # Total number of bids for the jth auction.
    arrival.count[jj] <- sum(raw.unorganized.data$auctionid == unique.auction.data[,1][jj])  
    # All the private values of the bidders for the jth process
    yy <-  raw.unorganized.data$bid[which(raw.unorganized.data$auctionid == unique.auction.data[,1][jj])]  
    
    Raw.biddata.matrix <- matrix(0, nrow = arrival.count[jj], ncol = 4)
    colnames(Raw.biddata.matrix) <- c("auction.index", "bid", "bidtime", "openbid")
    Raw.biddata.matrix[,1] <- rep(jj, arrival.count[jj])
    Raw.biddata.matrix[,2] <- yy
    Raw.biddata.matrix[,3] <- totts
    Raw.biddata.matrix[,4] <- rep(reserve.price[jj], arrival.count[jj])
    
    Raw.biddata.list[[jj]] <- Raw.biddata.matrix
  } # END of for(jj in 1: N.auction)
  
  Raw.data.list <- list(Raw.biddata.list = Raw.biddata.list,
                        N.auction = N.auction,
                        auction.window = auction.window,
                        reserve.price = reserve.price,
                        True.Method = "unknown")
  Raw.data.list$class = "SecondPriceAuction.Rawdata"
  
  return(Raw.data.list)
}


#------------------------------------------------------------------------------------
##  Function to compute the Non-parametric MLE of the consumer valuation  
##  distribution function F.
#------------------------------------------------------------------------------------
MLE.from.raw.data.2ndprice <- function(raw.data.list, reserve.price.Cutoff){
  if(raw.data.list$class != "SecondPriceAuction.Rawdata"){
    stop("This fucntion works only with data in the class
         'SecondPriceAuction.Rawdata'.
         See e.g., 'Data.Gen.2nd.Price.Raw' where we generate
         such data.")
  }

  data <- Data.Gen.2nd.Price.Processed(raw.data.list)
  LL <- length(data$pooled.observed.bids)
  Init.Values<- initialization.2nd.price(data, reserve.price.Cutoff)  
  theta.init <- (1-Init.Values$F.y)/ c(1,(1-Init.Values$F.y)[-length(Init.Values$F.y)])
  MLE<- MLE.2ndprice.init(data, lambda.in = Init.Values$lambda, theta.in = theta.init,
                          tol = 1e-5)
  ret <- list(Init.Values = Init.Values,
              MLE = MLE)
  ret$class = "SecondPriceAuctionInitAndMLE"
  return(ret)
}

