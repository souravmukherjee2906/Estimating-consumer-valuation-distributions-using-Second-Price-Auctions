#------------------------------------------------------------------------------------
##  Source the necessary functions from "auxiliary_functions.R" file. 
#------------------------------------------------------------------------------------
source("auxiliary_functions.R")

#------------------------------------------------------------------------------------
##  HulC1d_Ebay() uses asymptotic median bias value to construct
##  convex hull confidence interval for a univariate 
##  parameter. This is Algorithm 1 of the paper.
##  data is a data frame.
##  estimate is a function that takes a data frame as input 
##  and returns a one-dimensional estimate.
##  alpha is the level.
##  Delta is the median bias of the estimate(). 
##  randomize is a logical. If TRUE then the number of splits
##  is randomized. If FALSE, then the larger number of splits
##  is used.
#------------------------------------------------------------------------------------
HulC1d_Ebay <- function(data.in, eval.points, alpha = 0.1, MLE.Cutoff = 0){
  randomize = TRUE
  # alpha = 0.05; Delta = 0; randomize = TRUE
  nn <- length(data.in$Raw.biddata.list)
  Delta <- 0
	
  corrected.alpha <- alpha/length(eval.points) # Bonferroni Correction
  B1 <- solve_for_B(alpha = corrected.alpha, Delta = Delta, t = 0)
  B <- B1
	
  if(randomize){
    p1 <- (1/2 + Delta)^B1 + (1/2 - Delta)^B1
    B0 <- B1 - 1
    p0 <- (1/2 + Delta)^B0 + (1/2 - Delta)^B0
    U <- runif(1)
    tau <- (corrected.alpha - p1)/(p0 - p1)
    B <- B0*(U <= tau)+ B1*(U > tau)
  }
  if(B > nn){
    print(paste0("Delta = ", Delta, ", No. of splits = ", B, ", Sample size = ", nn))
    stop("Error: not enough samples for splitting!")
  }
  MLE_est <- matrix(0,nrow=length(eval.points),ncol= B)
  Init_est <- matrix(0,nrow=length(eval.points),ncol= B)
  # Right now we are using the same split but we will later use different splits.
  TMP <- split(sample(nn), sort((1:nn)%%B))
  idx <- 1
  for(idx in 1:B){
    # partioning real data
    temp.data  <- NULL
    temp.data$Raw.biddata.list<- data.in$Raw.biddata.list[TMP[[idx]]]
    temp.data$reserve.price <- data.in$reserve.price[TMP[[idx]]]
    temp.data$N.auction <- length(TMP[[idx]])
    temp.data$True.Method <- data.in$True.Method
    temp.data$reserve.price <- data.in$reserve.price
    temp.data$auction.window = data.in$auction.window
    temp.data$class =  "SecondPriceAuction.Rawdata"
    F.local<- MLE.from.raw.data.2ndprice(temp.data, MLE.Cutoff)
    
    MLE_est[,idx] <- approx(x = F.local$MLE$F.x, y = F.local$MLE$F.y, xout = eval.points,
			    method = "linear", yleft = 0, yright = 1)$y
    
    Init_est[,idx] <- approx(x = F.local$Init.Values$F.x, y = F.local$Init.Values$F.y,
			     xout = eval.points, method = "linear", yleft = 0, yright = 1)$y
  }
	
  CI.Init <- apply(Init_est, 1, range)
  CI.MLE <- apply(MLE_est, 1, range)
  rownames(CI.MLE)<- rownames(CI.Init)<- c("lwr", "upr")
  CI.Init[1,] <- pmax(CI.Init[1,] - log(2)*B/nn, 0) # expansion by log(2)/m (split size) to account for binomial dist
  CI.Init[2,] <- pmin(CI.Init[2,] + log(2)*B/nn, 1) # expansion to account for binomial dist
  CI.MLE[1, ] <- pmax(CI.MLE[1, ] - log(2)*B/nn, 0) # expansion by log(2)/m (split size) to account for binomial dist
  CI.MLE[2, ] <- pmin(CI.MLE[2, ] + log(2)*B/nn, 1) # expansion to account for binomial dist
	
  # we can use the monotonicity of F_0 to improve provide a confidence band for the whole function, 
  # suppose $\ell(t_1) \le f_0(t_1) \le u(t_1)$ and $\ell(t_2) \le f_0(t_2) \le u(t_2)$ for two points
  # $t_1 \le t_2\in[0,1]$ in the domain, then using the information $F_0(\cdot)$ is increasing and between
  # 0 and 1, we can conclude that $\overline{\ell}(t) \le f_0(t) \le \overline{u}(t)$ for all $t\in[0, 1]$, where
  # $$\begin{equation*}
  # \overline{\ell}(t) = \begin{cases}0,&amp;\mbox{for }t &lt; t_1,\\
  # \ell(t_1),&amp;\mbox{for }t_1 \le t &lt; t_2,\\
  # \ell(t_2), &amp;\mbox{for }t_2 \le t \le 1,\end{cases}
  # \quad\mbox{and}\quad 
  # \overline{u}(t) = \begin{cases}u(t_1),&amp;\mbox{for }t \le t_1,\\
  # u(t_2),&amp;\mbox{for }t_1 &lt; t \le t_2,\\
  # 1, &amp;\mbox{for }t_2 &lt; t \le 1.\end{cases}
  # \label{eq:interval-to-band-monotone}
  # \end{equation*}$$
  
  extra.info <- list(MLE_est,Init_est) # all the estimators based on the data splits.
  
  MakeBand <- function(temp){
    rval<- NULL
    rval$lwr <- approxfun(eval.points, cummax(temp[1,]), method = "constant", 
			  yleft = 0, yright = max(temp[1,]), f = 0, ties = "ordered")
    class(rval$lwr) <- c( "stepfun")
    rval$upr <- approxfun(eval.points, rev(cummin(rev(temp[2,]))), method = "constant", 
			  yleft =  min(temp[2,]), yright = 1, f = 1, ties = "ordered")
    class(rval$upr) <- c("stepfun")
    return(rval)
  }
	
  Band.CI.Init <- MakeBand(CI.Init)
  Band.CI.MLE <- MakeBand(CI.MLE)
  ret <- list(CI.Init = CI.Init, CI.MLE = CI.MLE, B = B, Band.CI.Init = Band.CI.Init, 
	      Band.CI.MLE = Band.CI.MLE, extra.info = extra.info)
  return(ret)
}

