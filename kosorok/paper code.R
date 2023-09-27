#Two-Sample Quantile Tests under General Conditions (Kosorok 1999)
#The null hypothesis tested is that the quantiles, with quantile probability p, are equal but that other aspects of the distributions may differ between the two samples

#https://bios.unc.edu/~kosorok/qtest.sas

require(survival)
library(dplyr)
#library(psignifit)

SurvQtest <- function(DATA1, Time1, DATA2, Time2, p, Q = 0.25, ALPHA = 0.05, OUT = "_QTEST") {
  # INPUT:
  # DATA1 (Required) = Random sample from F. It consists of two columns. First column contains the time to event. Second column contains censoring value( 1 = event, 0 = censoring ).
  # DATA2 (Required) = Random sample from G. It consists of two columns. First column contains the time to event. Second column contains censoring value( 1 = event, 0 = censoring ).
  # TIME1 (Required) = Variable containing time to event or last follow-up in any units in DATA1.
  # TIME2 (Required) = Variable containing time to event or last follow-up in any units in DATA2.
  # P (Required) = The probability of which the quantile would be tested. The null hypothesis is that the pth quantile of the F distribution is same as the pth quantile of the G distribution. .
  # Q (Default is 0.25) = Bandwidth adjustment term for kernel density estimation.
  # ALPHA (Default is 0.05) = Type I error.
  # OUT (Default name is _QTEST) = Data set name where the output would be stored.
  # OUTPUT:
  # result = Dataframe with NUM_F, NUM_G, XI_F, XI_G, DIFF, STD_ERR, LOWER, UPPER, CHISQ, PVALUE

  # Load data1 and data2 (assuming they are data frames)
  data1 <- DATA1
  data2 <- DATA2
  
  ##########
  # In R, the equivalent functionality to the proc lifetest in SAS is provided by the survfit function from the survival package (compute KM)
  
  # Run lifetest on data1
  out1 <- survfit(Surv(x, event == 0) ~ x, data = data1)
  #out1 <- survfit(Surv(x, event == 0) ~ 1, data = data1)
  
  #print(out1)
  
  out1_df <- data.frame(x = out1$time, event = data1$event, survival = out1$surv)
  #out1_df <- data.frame(time = out1$time, x = rep(unique(data1$x), each = length(out1$time)))
  
  # Run lifetest on data2
  out2 <- survfit(Surv(x, event == 0) ~ x, data = data2)
  #out2 <- survfit(Surv(x, event == 0) ~ 1, data = data2)
  out2_df <- data.frame(x = out2$time, event = data2$event, survival = out2$surv)
  #out2_df <- data.frame(time = out2$time, event = rep(unique(data2$x), each = length(out2$time)))
  
  # Filter out1_df and out2_df
  out1_df <- out1_df %>%
    filter(x > 0 & event== 0)
  
  out2_df <- out2_df %>%
    filter(x > 0 & event == 0)
  
  # Calculate the minimum survival times for data1 and data2
  mymin3 <- min(out1$time)
  mymin4 <- min(out2$time)
  
  # Calculate m = 1 - max(mymin3, mymin4)
  m <- 1 - max(mymin3, mymin4)
  print(m)
  
  # Check if m < P and display an error message if true
  if (m < p) {
    cat('ERROR: Probability <P> is not appropriate. Try another one.\n')
    stop("ERROR - detected in the input data to the function <SurvQtest>.")
  }
  
  # Sort data sets
  out1 <- arrange(data1, Time1)
  out2 <- arrange(data2, Time2)
  
  # Read data sets into matrices
  x <- as.matrix(out1)
  y <- as.matrix(out2)
  
  n1 <- nrow(x)
  n2 <- nrow(y)
  
  xtmp1 <- matrix(0, n1, 2) 
  # 2 column matrix: 
  #the 1st indicates whether the value in the 2nd column of the corresponding row in out1 is equal to 1
  #the 2nd represents the count of rows in out1 where the value in the first column is greater than or equal to the value in the corresponding row
  for (i in 1:n1) {
    xtmp1[i, 1] <- as.numeric(x[i, 2] == 1)
    xtmp1[i, 2] <- sum(x[, 1] >= x[i, 1])
  }
  
  xtmp2 <- numeric(n1)
  # vector of length n1 (=number of rows in out1), each element corresponds to a row in the ou1 dataset
  # computes proportion of times the value in the second column of the corresponding row of out1 is not equal to 1, among the times where the value in the first column is greater than or equal to the value in the corresponding row
  for (i in 1:n1) {
    xtmp2[i] <- 1 - xtmp1[i, 1] / xtmp1[i, 2]
  }
  
  xtmp3 <- numeric(n1)
  # vector of length n1
  #vector where each element is the cumulative product of specific elements from the xtmp2 vector
  for (i in 1:n1) {
    xtmp3[i] <- 1
    for (k in 1:i) {
      xtmp3[i] <- xtmp3[i] * xtmp2[k]
    }
  }
  
  xtmp4 <- numeric(n1)
  # vector of length n1
  # each element has the difference between 1 and the corresponding element of xtmp3
  for (i in 1:n1) {
    xtmp4[i] <- 1 - xtmp3[i]
  }
  
  xtmp5 <- numeric(n1)
  # vector of length n1
  # each element corresponds to the difference between consecutive elements from xtmp4
  xtmp5[1] <- 1 - xtmp4[1]
  for (i in 2:n1) {
    xtmp5[i] <- xtmp4[i] - xtmp4[i - 1]
  }
  
  # single numeric value representing the minimum index in the xtmp4 vector where the value is greater than or equal to the specified threshold p
  index1 <- min(which(xtmp4 >= p))
  
  # value from the first column of the matrix x at the row index1
  xxi <- x[index1, 1]
  
  # single numeric value representing the element from vector xtmp3 at position index1
  xhats <- xtmp3[index1]
  
  # vector where each element is: (difference between 2nd and 1st columns of xtmp1 for each row) times (values from the 2nd column of xtmp1 for each row)
  xden <- (xtmp1[, 2] - xtmp1[, 1]) * xtmp1[, 2]
  
  xtmp6 <- numeric(n1)
  # vector of length n1
  # each element corresponds to the division og the element in the 1st column of xtmp1 by the corresponding element from xden
  for (i in 1:n1) {
    if (xden[i] != 0) {
      xtmp6[i] <- xtmp1[i, 1] / xden[i]
    }
  }
  
  # Calculate phi: single numeric value
  phi <- (xhats^2) * sum((x[, 1] <= xxi) * xtmp6)
  
  
  # now we compute the equivalent vectors and matrix considering out 2 and y (instead of out1 and x)
  
  ytmp1 <- matrix(0, n2, 2)
  # matrix of size n2 by 2
  # similar to xtmp1 but for out2 and y (instead of out1 and x)
  for (i in 1:n2) {
    ytmp1[i, 1] <- as.numeric(y[i, 2] == 1)
    ytmp1[i, 2] <- sum(y[, 1] >= y[i, 1])
  }
  
  ytmp2 <- numeric(n2)
  for (i in 1:n2) {
    ytmp2[i] <- 1 - ytmp1[i, 1] / ytmp1[i, 2]
  }
  
  ytmp3 <- numeric(n2)
  for (i in 1:n2) {
    ytmp3[i] <- 1
    for (k in 1:i) {
      ytmp3[i] <- ytmp3[i] * ytmp2[k]
    }
  }
  
  ytmp4 <- numeric(n2)
  for (i in 1:n2) {
    ytmp4[i] <- 1 - ytmp3[i]
  }
  
  ytmp5 <- numeric(n2)
  ytmp5[1] <- 1 - ytmp4[1]
  for (i in 2:n2) {
    ytmp5[i] <- ytmp4[i] - ytmp4[i - 1]
  }
  
  index2 <- min(which(ytmp4 >= p))
  
  yxi <- y[index2, 1]
  yhats <- ytmp3[index2]
  
  yden <- (ytmp1[, 2] - ytmp1[, 1]) * ytmp1[, 2]
  
  ytmp6 <- numeric(n2)
  for (i in 1:n2) {
    if (yden[i] != 0) {
      ytmp6[i] <- ytmp1[i, 1] / yden[i]
    }
  }
  
  # Calculate gamma: single numeric value
  gamma <- (yhats^2) * sum((y[, 1] <= yxi) * ytmp6)
  
  #single numeric value representing the index in the vector where the value is greater than or equal to the specified threshold p
  index3 <- min(which(xtmp4 >= q))
  
  #retrieves the value from the x vector at the index specified by index3, in the first column.
  xq_xi <- x[index3, 1]
  
  
  #####################################################################
  # Computations:
  
  # Calculate h1
  h1 <- (n1^(-1/5)) * 4 * xq_xi
  
  # Calculate xcom (vector)
  xcom <- numeric(n1)
  # computing the function K(x) based on conditions on x
  for (i in 1:n1) {
    if (-1 <= (xxi - x[i, 1]) / h1 & (xxi - x[i, 1]) / h1 <= 0) {
      xcom[i] <- (((xxi - x[i, 1]) / h1) + 1)
    } else if (0 < (xxi - x[i, 1]) / h1 & (xxi - x[i, 1]) / h1 <= 1) {
      xcom[i] <- (1 - ((xxi - x[i, 1]) / h1))
    }
  }
  
  ########################3
  # Verify the computations (I just translated from the original code)
  
  # a1: vector that returns the minimum value element-wise between (xxi+h1) and (maximum value in the 1st column of x)
  a1 <- pmin(xxi + h1, pmax(x[, 1]))
  
  # b1: scalar value
  b1 <- max(xxi, 0)
  
  # c1: scalar value that returns the minimum value between xxi and the maximum value in the 1st column of the x vector
  c1 <- pmin(xxi, pmax(x[, 1]))
  
  # d1: scalar value representing the maximum between xxi-h1 and 0
  d1 <- max(xxi - h1, 0)
  
  # wf: scalar value
  wf <- (1 / h1) * ((1 + xxi / h1) * (a1 - b1) - (a1^2 - b1^2) / (2 * h1) + (1 - xxi / h1) * (c1 - d1) + (c1^2 - d1^2) / (2 * h1))
  
  # Calculate f
  f <- (1 / wf) * (1 / h1) * sum(xcom * xtmp5)
  
  index4 <- which(ytmp4 >= q)[1]
  yq_xi <- y[index4, 1]
  h2 <- (n2^(-1/5)) * 4 * yq_xi
  ycom <- ((yxi - y[, 1]) / h2 + 1) * (-1 <= (yxi - y[, 1]) / h2 & (yxi - y[, 1]) / h2 <= 0) +
    (1 - ((yxi - y[, 1]) / h2)) * (0 < (yxi - y[, 1]) / h2 & (yxi - y[, 1]) / h2 <= 1)
  
  a2 <- min(yxi + h2, max(y[, 1]))
  b2 <- max(yxi, 0)
  c2 <- min(yxi, max(y[, 1]))
  d2 <- max(yxi - h2, 0)
  wg <- (1 / h2) * ((1 + yxi / h2) * (a2 - b2) - (a2^2 - b2^2) / (2 * h2) + (1 - yxi / h2) * (c2 - d2) + (c2^2 - d2^2) / (2 * h2))
  g <- (1 / wg) * (1 / h2) * sum(ycom * ytmp5)
  
  beta <- (xxi - yxi)
  psi <- phi / (f * f) + gamma / (g * g)
  stat <- (beta^2) / psi
  
  mystat <- data.frame(n1 = n1, n2 = n2, xxi = xxi, yxi = yxi, psi = psi, stat = stat)
  write.table(mystat, file = "mystat.txt", row.names = FALSE, col.names = TRUE)
  
  # print result:
  alpha <- ALPHA
  NUM_F <- n1
  NUM_G <- n2
  XI_F <- xxi
  XI_G <- yxi
  DIFF <- xxi - yxi
  z <- qnorm(1 - alpha / 2)  # Z-score for the specified alpha level
  STD_ERR <- sqrt(psi)
  LOWER <- DIFF - z * sqrt(psi)
  UPPER <- DIFF + z * sqrt(psi)
  CHISQ <- stat
  PVALUE <- 1 - pchisq(stat, df = 1)
  
  result <- data.frame(
    NUM_F = NUM_F,
    NUM_G = NUM_G,
    XI_F = XI_F,
    XI_G = XI_G,
    DIFF = DIFF,
    STD_ERR = STD_ERR[1], # multiple values
    LOWER = LOWER[1], # multiple values
    UPPER = UPPER[1], # multiple values
    CHISQ = CHISQ[1], # multiple values
    PVALUE = PVALUE[1] # multiple values
  )
  
  write.table(result, file = "output.txt", row.names = FALSE, col.names = TRUE)
 
  return(result)
}

##########################
# Simulations

##########################
# Parameters
n0 = 500
n1 = 500
shape0 = 2.0
scale0 = 1
shape1 = 2.0
scale1 = 1
rate = 1/15

# Generating survival times, censoring status and censoring times
set.seed(1)
truetimes0_unsort = rweibull(n0,shape=shape0,scale=scale0) # True times: censored and uncensored (=time to event)
set.seed(1)
censtimes0_unsort = rexp(n0,rate=rate) # Censoring times
set.seed(1)
status0_unsort = as.numeric(truetimes0_unsort <= censtimes0_unsort) # Censoring status: 1 if uncensored, 0 if censored

alltimes0_unsort = pmin(truetimes0_unsort, censtimes0_unsort) # Observed times: min between true times and censored times

DATA1 <- data.frame("x" = truetimes0_unsort, "event" = status0_unsort)
Time1 <- alltimes0_unsort

set.seed(2)
truetimes1_unsort = rweibull(n1,shape=shape1,scale=scale1) # True times: censored and uncensored
set.seed(2)
censtimes1_unsort = rexp(n1,rate=1/15) # Censoring times
set.seed(2)
status1_unsort = as.numeric(truetimes1_unsort <= censtimes1_unsort) # Censoring status: 1 if uncensored, 0 if censored

alltimes1_unsort = pmin(truetimes1_unsort, censtimes1_unsort) # Observed times: min between true times and censored times

DATA2 <- data.frame("x" = truetimes1_unsort, "event" = status1_unsort)
Time2 <- alltimes1_unsort

# Quantile that will be compared: p
p = 0.5

q = 0.25

ALPHA = 0.05

# Running our function for testing the equality of quantile:
res = SurvQtest(DATA1, Time1, DATA2, Time2, p=p)
View(res)

# Percentage of data that is censored in each group:
taux_cens1 = (n0 - sum(DATA1$event))/n0

taux_cens2 = (n1 - sum(DATA2$event))/n1

#########
# Plotting the survival functions and quantile

#########
library(ggplot2)
library(ggfortify)
DATA1_PLOT = data.frame(DATA1)
DATA1_PLOT$groups = rep(1, n0)

DATA2_PLOT = data.frame(DATA2)
DATA2_PLOT$groups = rep(2, n0)

DATA_PLOT = rbind(DATA1_PLOT, DATA2_PLOT)

r_fit = survfit(Surv(DATA_PLOT$x, DATA_PLOT$event)~DATA_PLOT$groups, data=DATA_PLOT)
autoplot(r_fit ) + geom_hline(yintercept = p,  linetype="dashed") + 
  labs(x="Time", y=" Survival probability", color = 'Groups', fill = "Groups") 

