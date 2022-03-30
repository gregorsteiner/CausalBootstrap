

######### Functions ######### 


# function provided by Aronow, Green & Lee (2014)
sharp.var <- function(yt, yc, N=length(c(yt,yc)), upper=TRUE) {
  m <- length(yt)
  n <- m + length(yc)
  FPvar <- function(x,N) (N-1)/(N*(length(x)-1)) * sum((x - mean(x))^2)
  yt <- sort(yt)
  if(upper == TRUE) yc <- sort(yc) else
    yc <- sort(yc,decreasing=TRUE)
  p_i <- unique(sort(c(seq(0,n-m,1)/(n-m),seq(0,m,1)/m))) -
    .Machine$double.eps^.5
  p_i[1] <- .Machine$double.eps^.5
  yti <- yt[ceiling(p_i*m)]
  yci <- yc[ceiling(p_i*(n-m))]
  p_i_minus <- c(NA,p_i[1: (length(p_i)-1)])
  return(((N-m)/m * FPvar(yt,N) + (N-(n-m))/(n-m) * FPvar(yc,N)
          + 2*sum(((p_i-p_i_minus)*yti*yci)[2:length(p_i)])
          - 2*mean(yt)*mean(yc))/(N-1))
}



# function for the AGL variance
Var_AGL <- function(Y, D, n0, n1, N){
  # get S1^2 and S0^2
  S1.2 <- var(Y[D == 1])
  S0.2 <- var(Y[D == 0])
  
  # S01^2
  S01.2 <- S0.2 + S1.2 - 2 * sharp.var(Y[D == 1], Y[D == 0], upper = TRUE)
  
  # in total
  Var <- S0.2 / n0 + S1.2 / n1 - S01.2 / N
  return(Var)
  
}


# bootstrap function
causal_boot <- function(data, dep.var, treatment, N, B = 100, alpha = 0.05){
  
  # define treatment and outcome
  D <- data[, treatment]
  Y <- data[, dep.var]
  
  # remove NAs
  na.ind <- !is.na(Y) & !is.na(D)
  Y <- Y[na.ind]
  D <- D[na.ind]
  
  # some parameters
  p <- mean(D) # proportion of treated observations
  n <- nrow(data) # observations in the sample
  n1 <- sum(D)
  n0 <- n - n1
  N0 <- ceiling((n0 / n) * N)
  N1 <- N - N0
  
  
  # estimators
  tau.hat <- mean(Y[D == 1], na.rm = TRUE) - mean(Y[D == 0], na.rm = TRUE)
  sigma.hat <- sqrt(Var_AGL(Y, D, n0, n1, N))
  
  # generate empirical population
  dat.boot <- do.call(rbind, Map(function(d, n.int, N.int){
    # select treatment
    dat <- data.frame(Y = Y[D == d], D = D[D == d])
    
    # order in an increasing fashion
    dat <- dat[order(dat$Y), ]
    
    # define rank
    dat$U <- 1:nrow(dat) / n.int
    
    # loop over rows and include copies of them
    for (i in 1:(n.int - 1)) {
      # compute number of inclusions
      M <- ceiling(dat[i+1, "U"] * N.int) - ceiling(dat[i, "U"] * N.int)
      
      # assign
      if(i == 1) store <- cbind(Y = rep(dat[i, "Y"], M), D = rep(dat[i, "D"], M))
      
      else{
        # try to generate replications, if an error occurs, we simply skip this one
        append <- try(cbind(Y = rep(dat[i, "Y"], M), D = rep(dat[i, "D"], M)), silent = TRUE)
        
        # if the append object is a numeric matrix, we add it to the store object
        if(is.matrix(append) & is.numeric(append)) store <- rbind(store, append)
      }
      
    }
    
    # return
    return(store)
    
  }, c(0, 1), list(n0, n1), list(N0, N1)))
  
  
  # get empirical cdfs (F_1 and F_0)
  F1 <- ecdf(dat.boot[dat.boot[, "D"] == 1, "Y"])
  F0 <- ecdf(dat.boot[dat.boot[, "D"] == 0, "Y"])
  
  
  # impute missing potential outcomes
  dat.boot <- cbind(dat.boot,
                    "Y0" = ifelse(dat.boot[, "D"] == 0, dat.boot[, "Y"], quantile(F0, F1(dat.boot[, "Y"]))), 
                    "Y1" = ifelse(dat.boot[, "D"] == 1, dat.boot[, "Y"], quantile(F1, F0(dat.boot[, "Y"]))))
  
  T.store <- numeric(B)
  
  for (i in 1:B) {
    # randomly draw from population
    dat.boot.it <- dat.boot[sample(1:nrow(dat.boot), size = n), ]
    
    
    # randomly assign treatment (TRUE is treatment)
    Treat <- sample(c(FALSE, TRUE), size = n, replace = TRUE, prob = c(1 - p, p))
    # add to df
    dat.boot.it <- cbind(dat.boot.it,
                         "Treat.Boot" = as.numeric(Treat),
                         "Y.Boot" = ifelse(Treat, dat.boot.it[, "Y1"], dat.boot.it[, "Y0"]))
    
    # sample sizes
    n1.boot <- sum(Treat)
    n0.boot <- n - n1.boot
    
    # compute bootstrap sample statistics
    tau.hat.boot <- mean(dat.boot.it[Treat, "Y.Boot"]) - mean(dat.boot.it[!Treat, "Y.Boot"])
    sigma.hat.boot <- sqrt(Var_AGL(Y = dat.boot.it[, "Y.Boot"],
                                   D = dat.boot.it[, "Treat.Boot"],
                                   n0 = n0.boot, n1 = n1.boot, N = N))
    
    
    # compute t-ratio
    T.store[i] <- (tau.hat.boot - tau.hat) / (sigma.hat.boot / sqrt(n))
    
  }
  
  # get ecdf of t-ratios
  G <- ecdf(T.store)
  
  # compute and return CI
  CI <- c("Lower" = tau.hat - quantile(G, 1 - alpha, names = FALSE) * sigma.hat / sqrt(n),
          "Upper" = tau.hat - quantile(G, alpha, names = FALSE) * sigma.hat / sqrt(n))
  return(list("Confidence Interval" = CI,
              "ATE Estimator" = tau.hat,
              "Estimated Upper Bound SD" = sigma.hat))
  
}


