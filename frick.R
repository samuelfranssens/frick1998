# The first section of this script replicates tables 1 and 2 from 
# Frick, 1998, A better stopping rule for conventional statistical tests,
# available at: http://datacolada.org/wp-content/uploads/2015/11/5364-Frick-BRM-1998-A-better-stopping-rule-for-conventional-statistical-tests.pdf
# and brought to my attention by: https://rolfzwaan.blogspot.be/2015/05/p20-what-now-adventures-of-good-ship.html.
# Note that, unlike the original article, this scripts looks at the difference between two independent means (article: difference between two dependent means).
# All sample sizes are n's per cell.

# As you can read in the article, the power achieved by COAST depends on the lower bound that is chosen.
# The second section of this script calculates how many observations one would need to achieve about .80 power with COAST.
# Result: a lower bound of about 55% of the N required by a fixed N stopping rule seems enough to achieve .80 power.

# Both simulations also save estimates of Cohen's d.

# The script starts with a function to calculate p values and other statistics: ----------------------
differences <- function(x1,x2,ci = 0.95) {
  mean1 <- mean(x1)
  mean2 <- mean(x2)
  sd1 <- sd(x1)
  sd2 <- sd(x2)
  n1 <- length(x1)
  n2 <- length(x2)
  
  df=n1+n2-2
  sd.pooled <- sqrt(( sd1*sd1 * (n1-1) + sd2*sd2 * (n2-1))/(df))
  
  difference <- mean1 - mean2
  # lwr <- difference + qt((1-ci)/2, df = df) * sd.pooled * sqrt(1/n1 + 1/n2)
  # upr <- difference + qt((1+ci)/2, df = df) * sd.pooled * sqrt(1/n1 + 1/n2)
  
  d.est <- (difference) / sd.pooled
  # d.lwr <- cohend(d.est,n1,n2,"lwr")
  # d.upr <- cohend(d.est,n1,n2,"upr")
  
  p <- 2 * pt( -abs( d.est / sqrt(1/n1+1/n2)  ), df=df )
  
  return(cbind(difference,p,d.est))
}


# SECTION 1: ----------------------------------------------------------------------------
# TABLE 1: alpha for different lower bounds of COAST ------------------------------------
# TABLE 2: power of COAST and a comparison with the fixed-sample stopping rule ----------

number.of.experiments <- 1000
lower.bounds <- c(20,40,60)
effect.sizes <- c(0,.2,.4,.6)
simulations <- number.of.experiments * length(lower.bounds) * length(effect.sizes)

table <- data.frame( matrix( data=NA, ncol=5, nrow=simulations) )
names(table) <- c("effectsize","lowerbound","N","p","d")

begin <- Sys.time() # how long does the simulation take
k <- 1 # iteration

for (lowerbound in lower.bounds){
  for (i in 1:number.of.experiments){
    for (effect.size in effect.sizes){
      
      # create two independent normal distributions with true d = 0
      condition1 <- rnorm (lowerbound,mean=0,sd=1) 
      condition2 <- rnorm (lowerbound,mean=effect.size,sd=1)
      
      # calculate p-value and decide whether to continue testing
      result <- data.frame(differences(condition1,condition2))
      while(result$p > .01 & result$p < .36){   
        condition1 <- c(condition1, rnorm(1,mean=0,sd=1)) # test after every 1 additional subject per condition
        condition2 <- c(condition2, rnorm(1,mean=effect.size,sd=1))
        result <- data.frame(differences(condition1,condition2))
      }
      
      # save results
      table$effectsize[k] <- effect.size
      table$lowerbound[k] <- lowerbound
      table$N[k] <- length(condition1)
      table$p[k] <- result$p
      table$d[k] <- result$d.est
      
      # keep an eye on progress
      percent <- simulations/100
      if(k %% percent == 0){ print(paste(k/percent,"% completed",sep=""))}
      k <- k+1
    }
  }
}

duration <- Sys.time()-begin # how long does the simulation take
# write.csv(table,"table12.csv") # save the results of the simulation to a file

# return power (or alpha for d = 0) per lower bound & effect size
results <- data.frame(matrix(data=NA,nrow=length(lower.bounds)*length(effect.sizes),ncol=3))
names(results) <- c("effectsize","lowerbound","power")

library(pwr) # to calculate power

j <- 1 # iteration
for (i in lower.bounds){
  for (k in effect.sizes){
    table2 <- subset(table,table$lowerbound==i & table$effectsize==k)
    
    results$effectsize[j] <- k
    results$lowerbound[j] <- i
    results$power[j] <- sum(table2$p<.01) / length(table2$p)
    results$average.N[j] <- ceiling(mean(table2$N))
    results$average.d[j] <- round(mean(table2$d)*-1,2)
    
    if(k>0) # don't calculate power when true d = 0
    { 
    x  <- pwr.t.test(d=k,power=results$power[j],sig.level=.05,type="two.sample",alternative="two.sided") 
    results$fixedN[j] <- ceiling(x$n) # returns the sample size needed to achieve the same power as COAST with fixed N
    }
    else {results$fixedN[j] <- NA}
    
    j <- j+1
  }
}

results$efficiency <- round((results$fixedN/results$average.N-1)*100,2) # fixed N vs COAST
results[order(results$effectsize),]
# write.csv(results[order(results$effectsize),],"table12_summary.csv")


# SECTION 2: -------------------------------------------------------------------
# Percentage of fixed N required as lower bound to achieve certain power -------
number.of.experiments <- 5000
effect.sizes <- c(.2,.25,.30,.35,.4,.45,.50) # cohen's d
lower.bounds <- c(.40,.45,.50,.55,.60,.65,.70) # percentage of fixed N sample size
desiredpower <- .80

simulations <- number.of.experiments * length(lower.bounds) * length(effect.sizes)

table <- data.frame( matrix( data=NA, ncol=6, nrow=simulations) )
names(table) <- c("effectsize","lowerbound","lowerbound.actual","N","p","d")

begin <- Sys.time() # how long does the simulation take
k <- 1 # iteration

for (effect.size in effect.sizes){
  power.n <- pwr.t.test(d=effect.size,power=desiredpower,sig.level=.05,type="two.sample",alternative="two.sided") 
  requiredsamplesize <- ceiling(power.n$n) # required sample size to achieve desired power with fixed N
  
  for (lowerbound in lower.bounds){
    for (i in 1:number.of.experiments){ # think of i as one experiment
      lowerbound2 <- ceiling(lowerbound*requiredsamplesize) # actual lower bound (= percentage of required sample size)
      # create two independent normal distributions with true d = 0
      condition1 <- rnorm (lowerbound2,mean=0,sd=1) 
      condition2 <- rnorm (lowerbound2,mean=effect.size,sd=1)
      
      # calculate p-value and decide whether to continue testing
      result <- data.frame(differences(condition1,condition2))
      while(result$p > .01 & result$p < .36){   
        condition1 <- c(condition1, rnorm(1,mean=0,sd=1)) # test after every 1 additional subject per condition
        condition2 <- c(condition2, rnorm(1,mean=effect.size,sd=1))
        result <- data.frame(differences(condition1,condition2))
      }
      
      # save results
      table$effectsize[k] <- effect.size
      table$lowerbound[k] <- lowerbound
      table$lowerbound.actual[k] <- lowerbound2
      table$N[k] <- length(condition1)
      table$p[k] <- result$p
      table$d[k] <- result$d.est
      
      # keep an eye on progress
      percent <- simulations/100
      if(k %% percent == 0){ print(paste(k/percent,"% completed",sep=""))}
      k <- k+1
    }
  }
}

duration <- Sys.time()-begin # how long does the simulation take
# write.csv(table,"lowerbound.csv") # save the results of the simulation to a file
# This takes a long time to simulate, at least on my laptop. 
# Feel free to download the file from here: https://www.dropbox.com/s/8eawogo5vjiqg64/lowerbound.csv?dl=0


# return power per lower bound & effect size
x <- data.frame(matrix(data=NA,nrow=length(lower.bounds)*length(effect.sizes),ncol=5))
names(x) <- c("effectsize","fixedN","lb.percent","lb.actual","power")

library(pwr)
j <- 1 # iteration
for (i in lower.bounds){
  for (k in effect.sizes){
    table2 <- subset(table,table$lowerbound==i & table$effectsize==k)
    
    x$effectsize[j] <- k
    x$lb.percent[j] <- i
    x$lb.actual[j] <- mean(table2$lowerbound.actual) # actual lower bound
    
    y  <- pwr.t.test(d=k,power=desiredpower,sig.level=.05,type="two.sample",alternative="two.sided") # required sample size to achieve desired power w fixed N
    x$fixedN[j] <- ceiling(y$n) # save required sample size
    
    x$power[j]     <- sum(table2$p<.01) / length(table2$p) # observed power w optional stopping
    x$average.N[j] <- ceiling(mean(table2$N)) # average N w optional stopping
    x$average.d[j] <- round(mean(table2$d) *-1,2) # average N w optional stopping
    
    # at N = lower bound, we can
    # find p<.01 and reject H0
    # find p>.36 and accept H0
    # find .01<p<.36 and continue testing
    x$nr_of_experiments[j] <- number.of.experiments
    x$lb_stop[j]     <- length(table2$p[table2$N==x$lb.actual[j]])        # how many times do we stop at lower bound
    x$lb_stop_reject[j]   <- sum(table2$p[table2$N==x$lb.actual[j]]<.01)  # how many times do we stop & reject
    x$lb_continue[j] <- x$nr_of_experiments[j] - x$lb_stop[j]             # how many times do we continue at lower bound
    x$lb_continue_reject[j] <- sum(table2$p[table2$N>x$lb.actual[j]]<.01) # how many times do we continue & reject
    
    j <- j+1
  }
}

x$efficiency  <- round((x$fixedN/x$average.N-1)*100,2) # fixed N vs COAST
results <- x[order(x$effectsize),]
results$effectsize <- factor(results$effectsize)
# write.csv(results,"lowerbound_summary.csv")