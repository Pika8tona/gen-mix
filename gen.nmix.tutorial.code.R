# Tutorial on fitting the generalized N-mixture models in R

# adapted from David Dail's (March 2010) gen.nmix tutorial

# load the R file gen.nmix.R into an active R session, using the command

source("<pathway>/gen.nmix.R")

# verify the 3 files contained in this archive
  
ls()

# the 2 functions require the logit and reverse-logit transform, built by:
  
expit = function(x) { return((exp(x))/(1+exp(x))) }
logit = function(x) { return(log(x/(1-x))) }

# build covariance matrices for lambda and p

data = mallard.data
elev <- data[,"elev"]
length <- data[,"length"]
forest <- data[,"forest"]
n.it <- data[,c("count1","count2","count3")]
R <- nrow(n.it)
T <- ncol(n.it)
DATE <- data[,c("date1","date2","date3")]
IVEL <- data[,c("ivel1","ivel2","ivel3")]
DATE2 <- DATE^2


# reformat X and Z

X.const = rep(1,R)
X.lam = cbind(X.const,elev,forest)

p.date.lin = as.vector(t(DATE))
Z.const = rep(1,R*T)
Z.p = cbind(Z.const,p.date.lin)


# reconstruct the "sampling date" matrix 
  
DATE.vec = c(DATE[,1],DATE[,2],DATE[,3])
stds = sort(na.exclude(unique(DATE.vec)))
diffs = numeric(length(stds)-1)
for(i in 1:(length(stds)-1)){ diffs[i] = stds[i+1] - stds[i] }

diffs.days = round(diffs/min(diffs))
day.unique = numeric(length(diffs.days)+1)

day.unique[1]=1
for(i in 1:(length(diffs.days))){ day.unique[i+1] = day.unique[i] + diffs.days[i] }
day.unique = c(day.unique,NA)
DATE.mat = cbind(DATE.vec,seq(1:length(DATE.vec)))
DATE.sort = DATE.mat[order(DATE.mat[,1]),]
dup = duplicated(DATE.sort[,1])

j=1
for(i in 1:length(DATE.vec)){
  if(dup[i]==TRUE) DATE.sort[i,1]=DATE.sort[i-1,1]
  if(dup[i]!=TRUE) { 
    DATE.sort[i,1]= day.unique[j]
    j=j+1
  }
}

DATE.mat2 = DATE.sort[order(DATE.sort[,2]),]
DATE.2 = matrix(DATE.mat2,nrow=239,ncol=3,byrow=FALSE)


# construct a "date" matrix that provides the primary period number associated with each sample

DATE.3 = ceiling(DATE.2/30)


# specify the search method for optim()
  
Method = "BFGS"

# get info related to optim()

help(optim)

# specify the cut-off value K that will approx the infinite summations in the N-mixture models
# K = 200 was used in Dail and Madsen (2011), but K = 40 gives similar param ests in less time
  
K.lim = 40


# fit the orig N-mixture model w/negative binomial prior and no covariates (null model)
  
model.null = optim(c(0,0,0), nmix.mig, method=Method, hessian=TRUE, n=n.it,
                       X=X.const, Z=Z.const, migration="none", prior="NB", Date=DATE.3, K=K.lim)

# make sure optim() converged ("0" is good)

model.null$conv

# check the stability of the MLEs by calculating the condition number, the ratio of the largest 
# to the smallest eigenvalues of the hessian matrix (<=0 is bad--maybe overfitted)
  
ev.null = eigen(model.null$hessian)$values
cn.null = max(ev.null)/min(ev.null)

# get the MLEs

model.null$par

# back-transform the MLEs (optim() finds MLEs for log(lambda), logit(p) and log(dispersion param))
  
lambda.est = exp(model.null$par[1])
p.est = expit(model.null$par[2])
c(lambda.est, p.est)

# find the asymptotic standard errors of the param ests by solving the Hessian matrix
  
se = sqrt(diag(solve(model.null$hess)))

# get the min value of the negative log likelihood (optim() is minimizing the negative loglike)

nll = model.null$val

# calculate the AIC score
  
aic = nll + 2*length(model.null$par)

# fit the N-mixture model with covariates for lambda and p, starting at the MLEs from the null
# model and inserting "0"s in the start vec for the number of covars in the X and Z matrices

model.closed = optim(c(model.null$par[1],0,0,model.null$par[2],0,model.null$par[3]),
                          nmix.mig, method=Method, hessian=TRUE, n=n.it, X=X.lam, Z=Z.p, 
                          migration="none", prior="NB",Date=DATE.3,K=K.lim)

model.closed$conv
ev.closed = eigen(model.closed$hessian)$values
cn.closed = max(ev.closed)/min(ev.closed)


# fit the generalized model, starting at the MLEs from the closed model

model.open = optim(c(model.closed$par[1:5],-2,2,model.closed$par[6]), nmix.mig, 
                        method=Method, hessian=TRUE, n=n.it, X=X.lam, Z=Z.p, migration="constant", 
                        prior="NB", Date=DATE.3, K=K.lim)

model.open$conv
ev.open = eigen(model.open$hessian)$values
cn.open = max(ev.open)/min(ev.open)


# obtain the p-value from the closure test

t.stat = 2*(model.closed$val - model.open$val)
hess.mig = model.open$hess[6:7,6:7]
obs.inf = -hess.mig
prop = acos(obs.inf[1,2]/(sqrt(obs.inf[1,1]*obs.inf[2,2])))/(2*pi)
prop.0 = 0.5-prop
prop.1 = 0.5
prop.2 = prop
p.value = prop.0*(0) + prop.1*(1-pchisq(t.stat,1)) + prop.2*(1-pchisq(t.stat,2))

# obtain abundance ests and asymptotic SEs for every primary period and model, plus 
# 95% CIs for pop dyn params; supply number of primary periods "T"

ests.closed = ests(model.closed$par,model.closed$hess, migration="none", n=n.it,
                        X=X.lam, Z=Z.p, T=3, prior="NB")
ests.open = ests(model.open$par, model.open$hess, migration="constant", n=n.it, 
                        X=X.lam, Z=Z.p, T=3, prior="NB")

# the first T values are the estimated total abundance, and the next T values are the 
#asymptotic standard errors associated with each estimate, respectively

# the next 2 values are the back-transformed, mean estimates for lambda (across sites) 
# and detection probability (across sites and surveys)

# the last 2 values are the back-transformed estimates of gamma and omega (mig and survival)

gamma = ests.open[9]
omega = ests.open[10]

# obtain asymptotic SEs of the dynamic parameter estimates (before back-transformation)

se = sqrt(diag(solve(model.open$hess)))
se.gamma = se[6]
se.omega = se[7]

# compute asymptotic 95% confidence intervals for gamma and omega

gamma.ci = exp( c( model.open$par[6]- 1.96*se.gamma , model.open$par[6] +
1.96*se.gamma))
omega.ci = expit( c( model.open$par[7] - 1.96*se.omega, model.open$par[7] +
1.96*se.omega))
