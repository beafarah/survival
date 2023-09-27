require(survival)
n0 = 100
n1 = 100
n = n0+n1

shape0 = 0.5
scale0 = 5
shape1 = 0.5
scale1 = 5

truetimes0_unsort = rweibull(n0,shape=shape0,scale=scale0) # censored and uncensored
censtimes0_unsort = rexp(n0,rate=1/15) 
status0_unsort = as.numeric(truetimes0_unsort <= censtimes0_unsort) # 1 if uncensored

alltimes0_unsort = pmin(truetimes0_unsort, censtimes0_unsort)
Tsort0 = sort(alltimes0_unsort, index.return = TRUE)
alltimes0 = Tsort0$x # observed times (min between death and censored)
status0 = status0_unsort[Tsort0$ix]
group0 = rep(0, n0)

truetimes1_unsort = rweibull(n1,shape=shape1,scale=scale1) # censored and uncensored
censtimes1_unsort = rexp(n1,rate=1/15) 
status1_unsort = as.numeric(truetimes1_unsort <= censtimes1_unsort)

alltimes1_unsort = pmin(censtimes1_unsort, truetimes1_unsort)

Tsort1 = sort(alltimes1_unsort, index.return = TRUE)
alltimes1 = Tsort1$x
status1 = status1_unsort[Tsort1$ix] # corresponding status for the sorted times
group1 = rep(1, n1)

taus_sort = sort(c(alltimes0, alltimes1), index.return = TRUE)
taus = taus_sort$x

status_unsort = c(status0, status1) # all censoring status
status = status_unsort[taus_sort$ix]

groups_unsort = c(group0, group1)
groups = groups_unsort[taus_sort$ix]

data <- data.frame(taus, status, groups)

# computing number of deaths and number at risk for each group for each time tauj:
m = length(taus)

d0 = c()
d1 = c()
cens0 = c()
cens1 = c()

for(i in 1:length(taus)){
  val0 = 0
  # number of deaths at each tauj: uncensored and betweeen tauj and tau_{j-1}
  countd0 = sum(as.numeric(alltimes0[status0==1]==taus[i])) # uncensored times are less than tauj
  d0 <- c(d0,countd0)
  
  countc0 = sum(as.numeric(alltimes0[status0==0]==taus[i])) # uncensored times are less than tauj
  cens0 <- c(cens0,countc0)
  
  countd1 = sum(as.numeric(alltimes1[status1==1]==taus[i])) # uncensored times are less than tauj
  d1 <- c(d1,countd1)
  
  countc1 = sum(as.numeric(alltimes1[status1==0]==taus[i])) # uncensored times are less than tauj
  cens1 <- c(cens1,countc1) 
}

d = d0+d1

m = length(taus)

# computing number at risk

risk0 = rep(n0,m)
risk1 = rep(n1,m)
risk0[2:m] = risk0[2:m] - cumsum(d0+cens0)[1:(m-1)]
risk1[2:m] = risk1[2:m] - cumsum(d1+cens1)[1:(m-1)]

Y = risk0+risk1 #seq(n,1,by=-1)

E = sum(d*risk1/Y)

#Calculate the test statistic
#Z = sum(d1*risk0/Y - d0*risk1/Y)
#V = sum(d*risk0*risk1/Y^2)
#X = Z^2/V
#pchisq(X,1,lower.tail = F)

#https://web.stanford.edu/~lutian/coursepdf/unitweek3.pdf
O = sum(d1)

V = ((Y - risk1)/(Y-1))*risk1*(d/Y)*(1-d/Y)

Z = (O-E)/sqrt(sum(na.omit(V))) # follows (asymptotically) a N(0,1)

Z2 = Z^2 # test statistic that follows (asymptotically) a chi squared with one degree of liberty
#Get the p-value: smallest alpha for which H0 is rejected. H0: hazard functions are equal 
pchisq(Z2,1,lower.tail = F)
print(Z2)
survdiff(Surv(taus, status)~groups, data=data)

# test level
# comparing with chi squared distributions:
alpha = 0.05
critical_level_chisq = pchisq(1-alpha/2,1)*sum(na.omit(V))

# rejects when Z2*V>critical_level
Z2*sum(na.omit(V))


###########################################################################################################
# testing
delta = 2
mu = 1/2
# with no repeated times:
#alltimes: censored and uncensored!

# verify the seq!!!!
#scalier_fun0 <- stepfun(alltimes0,c(1,seq((n0-1)/n0, 1/n0, by = -1/n0))) #verify the steps
scalier_fun0 <- stepfun(alltimes0,c(1,seq((n0-1)/n0, 1/n0, length.out= n0))) #verify the steps
scalier_fun1 <- stepfun(alltimes1,c(1,seq((n1-1)/n1, 1/n1, length.out = n1))) 
scalier_fun <- stepfun(taus,c(1,seq((n-1)/n, 1/n, length.out = n))) 

# uncensored times of deaths
timed0_uncens = alltimes0[status0 == 1]
timed1_uncens = alltimes1[status1 == 1]
timed_uncens = taus[status == 1]

var_h1 = (1-mu)*delta*mean((scalier_fun0(timed1_uncens)* scalier_fun1(timed1_uncens))/(scalier_fun(timed1_uncens)^2 )) + mu * mean(scalier_fun0(timed1_uncens)/(scalier_fun(timed1_uncens)^2))

var_h0 = V # previously computed variance (under H0)

z_quant = qchisq(1-alpha/2, df = 1) # quantile of order 1-alpha/2 of chisq

c = sqrt(n/(n0*n1))

centering_term = var_h0*(delta - 1)

# power: 1-beta
beta = 1 - pchisq((var_h0*z_quant - centering_term)/(var_h1) , df = 1) - pchisq((-var_h0*z_quant - centering_term)/(var_h1) , df = 1)
mean(na.omit(beta))