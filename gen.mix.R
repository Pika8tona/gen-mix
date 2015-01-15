## gen.mix.R

# contains "nmix.mig", "ests", and "Mallard.data"
# nmix.mig is a function that returns the negative log likelihood for the generalized model
# ests is a function to return estimates of total abundance during each primary period
# Mallard.data is a point count data set originally published by Kery, Royle, Schmid (2005) 


nmix.mig = function(vars,n,X,Z,migration="none",prior="poisson",Date=matrix(rep(1:length(n[1,]),length(n[,1])),nrow=length(n[,1]),ncol=length(n[1,]),byrow=TRUE), K=200) {
  # function returns the negative log likelihood of the generalized N-mixture model. 
  # David Dail, revised March 2, 2010
  #
  # migration is one these: "none", "constant", "autoreg", "trend1","reshuf".
  # migration="none" assigns r=0, w=1, so G~Pois(0) [ie, P(G=0)=1 ] and S~Binom(Nit-1,1) [ie, P(S=Nit-1) = 1 ];
  # migration="constant" has G~Pois(e^r0);
  # migration="autoreg" has G~Pois(e^r0*N.it-1).
  # migration="trend1" has r = (1-w)*lambda[1] ... needs no covariates for lambda
  # migration="reshuf" has r=lambda, w=0 ... needs no covariates for lambda
  #
  # X is the covariate matrix for log(lambda.i), including a column of 1's if an intercept is desired
  # Z is the covariate matrix for logit(p.it) [detection probability], which is a (R*T) X (no. of covariates) matrix.
  # The covariates for p.it need to be recorded across rows (sites): X11, X12, ..., X1T, X21, X22..., X2T, ..., XRT; 
  # need to include a column of 1's to include an intercept. 
  #
  # prior is either "poisson" or "NB".
  # 
  # w is survival probability, r is Poisson entering migration rate (assumed independent here)
  # Note that covariates for the migration parameters (r,w) are not accommodated in this program code.
  # 
  # vars is B0, B1,... Bk for each of: log(lambda.i), logit(p.it); then log(r0),logit(w) if migration != "none"; then log(dispersion).
  #
  # n is the matrix of number observed: sites are rows, and successive sampling occasions go across the columns,
  #  with NA (missing) entered at the end of row if less samples are obtained for that site (row) than the
  # maximum obtained for any site.
  #
  # Date is a matrix of the sampling dates (positive integers); it is best to have the first sampling occasion to be 1,
  # with NA listed (only when n.it is also missing) for missing data at the end of each row.
  #
  # K is the cut-off for the summations over N.it; the likelihood can be sensitive to this choice.
  # It is good practice to choose a value that is "high enough" so that increasing K does not change the resulting likelihood.
  # (Higher K's will take more computing power & time).
  
  n.it= n
  R = length(n[,1])
  nsite=R
  #T.i is the number of sampling occasions per site
  T.i = numeric(R)
  for(i in 1:R){
    T.i[i] = sum( !is.na(n.it[i,]))
  }
  T=max(T.i)
  #Delta.it is the time difference (no. of primary periods) between sampling occasions - could be 0
  Delta.it = matrix(0, nrow=R, ncol=(T-1))
  for(i in 1:R){
    for(t in 2:T){
      Delta.it[i,(t-1)] = Date[i,t] - Date[i,(t-1)]
    }
  }
  
  J.it= matrix(1,nrow=R,ncol=T)
  for(i in 1:R){
    for(t in 1:T){
      if(!is.finite(n.it[i,t])) J.it[i,t]=0
    }
  }
  
  T = max(T.i)
  X.lam.i=as.matrix(X)  # in null case this is a 1xR vector of ones
  X.theta.it=as.matrix(Z) # in null case this is a 1xRT vector of ones
  
  # obtain the parameter values for each site & sampling occasion
  
  lam.i = exp(X.lam.i%*%vars[1:length(X.lam.i[1,])])
  p.it.vec = expit(X.theta.it%*%vars[(1+length(X.lam.i[1,])):(length(X.lam.i[1,])+length(X.theta.it[1,]))])
  p.it = matrix(p.it.vec,nrow=R,ncol=T ,byrow=TRUE)
  
  if(prior!="poisson") {
    disp=exp(vars[length(vars)])
    disp2= max(disp, 0.0000001)
  }
  p.mat = matrix(p.it,nrow=R,ncol=T,byrow=TRUE)
  K.num = 0:K
  
  if(migration == "constant"){
    r0 = vars[1+(length(X.lam.i[1,])+length(X.theta.it[1,]))]
    w = expit(vars[2+(length(X.lam.i[1,])+length(X.theta.it[1,]))])
    r1 = 0
    r.it = exp(r0)*(0:K)^(r1)
  }
  if(migration == "autoreg") {
    r0 = vars[1+(length(X.lam.i[1,])+length(X.theta.it[1,]))]
    w = expit(vars[2+(length(X.lam.i[1,])+length(X.theta.it[1,]))])
    r1 = 1 # allows G.it ~ Pois(e^r0 * N.it-1)
    r.it = exp(r0)*(0:K)^(r1)
  }
  if(migration == "trend1") {
    r0 = 0
    w = expit(vars[1+(length(X.lam.i[1,])+length(X.theta.it[1,]))])
    r1 = 0 # allows G.it ~ Pois(e^r0 * N.it-1)
    r.it = (1-w)*lam.i[1] # requires no covariates for lambda
  }
  if(migration == "reshuf") {
    r0 = 0
    w = 0
    r1 = 0 # allows G.it ~ Pois(e^r0 * N.it-1)
    r.it = lam.i[1] # requires no covariates for lambda
  }
  
  if(migration == "none"){
    r0 = 0
    w = 1
    r.it = rep(0,(K+1))
  }
  
  # some preliminaries we'll need for each site, may as well build them now:
  
  Pois1 = rep(0:-K,(K+1)^2) + rep(0:K, each=(K+1)^2)
  r.it.rep = rep(rep(r.it,each=(K+1)),(K+1))
  Pois.mat = matrix(dpois(Pois1,r.it.rep),nrow=(K+1)^2,ncol=(K+1),byrow=TRUE)
  ####
  
  Bin1 = rep(0:(K),(K+1)^2)
  Bin2 = rep(rep(0:(K),each=(K+1)),(K+1))
  Bin.mat = matrix(dbinom(Bin1,Bin2,w),nrow=(K+1)^2, ncol=(K+1),byrow=TRUE)
  
  Pt.Agoal2 = matrix((Pois.mat*Bin.mat)%*%rep(1,(K+1)),nrow=(K+1),ncol=(K+1),byrow=FALSE)
  Pt.A2.array = array(0,dim=c(1+max(na.omit(Delta.it)),(K+1),(K+1)))
  Pt.A2 = diag((K+1))
  Pt.A2.array[1,,]=Pt.A2
  # using the Chapman-Kolmogorov equations:
  max.delta = max(max(na.omit(Delta.it)),max(na.omit(Date[,1]-1)))
  for(i in 2:(1+max.delta)){
    Pt.A2 = Pt.A2%*%Pt.Agoal2
    Pt.A2.array[i,,]=Pt.A2
  }
  
  probs.i = numeric(R)
  
  #Construct the likelihood, for each site (follows the supplemental appendix of Dail & Madsen.
  for (i in 1:R){
    
    #Make A
    
    T=T.i[i]
    p.t = p.it[i,]
    lam = lam.i[i]  
    delta.t = Delta.it[i,]
    n.t = n.it[i,]  
    J.t = J.it[i,]
    if(T>0){
      #make g_{1}(N_{iT})
      
      if(J.t[T]==0)break #this shouldn't happen...should only have observations in the data matrix!
      Bin.first = rep(c(n.t[(T-(J.t[T]-1)):T]), each=(K+1))
      Bin.index = rep(K.num,J.t[T])
      B=matrix(dbinom(Bin.first,Bin.index,p.t[T],log=TRUE),nrow=(K+1),ncol=(J.t[T]),byrow=FALSE)
      #B[,is.na(n[i,])]<-0
      g1.T = exp(B%*%rep(1,length(J.t[T])))
      g1.T.mat = matrix(rep(g1.T,(K+1)),nrow=K+1,ncol=K+1,byrow=TRUE)
      # make g_{3}
      # if T=1 then what? skip this part, set A.goal = g1.T.
      if(T>1){
        if(delta.t[T-1] == 0) { 
          A.goal = matrix(g1.T,nrow=1,ncol=(K+1))
        }
        if(delta.t[T-1] != 0) {
          g3.T = Pt.A2.array[delta.t[T-1]+1,,]
          A.goal = matrix((g1.T.mat * g3.T)%*%rep(1,(K+1)),nrow=1,ncol=(K+1),byrow=TRUE)
        }
      }
      if(T==1){
        A.goal = rep(1,(K+1))
      }
      
      #######
      # part B
      
      B.goal=A.goal
      
      hold.j = 1
      if(T>2) {
        for(reps in 1:(T-2)) {
          counter = T-1-reps
          Bin.first = rep(c(n.t[(T-hold.j-(J.t[counter+1]-1)):(T-hold.j)]), each=(K+1))
          Bin.index = rep(K.num,J.t[counter+1])
          B=matrix(dbinom(Bin.first,Bin.index,p.t[counter+1],log=FALSE),nrow=(K+1),ncol=(J.t[T]),byrow=FALSE)
          g1.T = apply(B,1,prod)
          g1.T.mat = matrix(rep(g1.T,(K+1)),nrow=K+1,ncol=K+1,byrow=TRUE)
          hold.j = hold.j+1
          # make g_{3}
          if(delta.t[counter]==0){
            B.goal = B.goal*g1.T
          }      
          if(delta.t[counter]!=0){
            g3.T = Pt.A2.array[delta.t[counter+1]+1,,]
            #get g.star(N_{iT-1})
            
            pt.B3 = matrix(rep(B.goal,(K+1)),nrow=(K+1),ncol=(K+1),byrow=TRUE)
            B.goal = matrix((g1.T.mat * g3.T * pt.B3)%*%rep(1,(K+1)),nrow=1,ncol=(K+1),byrow=TRUE)
          }
        }
      }
      
      
      ####
      # part C
      
      Bin.first = rep(c(n.t[(1:(T-hold.j))]), each=(K+1))
      Bin.index = rep(K.num,J.t[1])
      B=matrix(dbinom(Bin.first,Bin.index,p.t[1],log=FALSE),nrow=(K+1),ncol=(J.t[T]),byrow=FALSE)
      g1.T = apply(B,1,prod)
      
      # now we have g_{2}
      if(prior!="NB"){
        pt.C2.i = matrix(dpois((0:K),rep(lam.i[i],K+1)),nrow=1,ncol=(K+1),byrow=TRUE)
      }
      if (prior=="NB"){
        pt.C2.i = matrix(dnbinom((0:K),mu=rep(lam.i[i],K+1),size=disp2),nrow=1,ncol=(K+1),byrow=TRUE)
      }
      
      pt.C2 = pt.C2.i/sum(pt.C2.i)
      pt.C3 = B.goal
      g1.T = as.vector(g1.T)
      pt.C2 = as.vector(pt.C2)
      pt.C3 = as.vector(pt.C3)
      # if the first sampling period for this site was the first overall period:
      if(Date[i,1] == 1){    
        probs.i[i] = sum((g1.T)*pt.C2*(pt.C3))
      } 
      
      # if the first sampling period for this site was NOT the first overall period:
      if(Date[i,1]!=1){
        g1.T.mat = matrix(rep(g1.T,(K+1)),nrow=K+1,ncol=K+1,byrow=TRUE)
        pt.B3 = matrix(rep(B.goal,(K+1)),nrow=(K+1),ncol=(K+1),byrow=TRUE)
        g4 = matrix((g1.T.mat* Pt.A2.array[Date[i,1],,]*pt.B3)%*%rep(1,(K+1)),nrow=1,ncol=(K+1),byrow=TRUE)
        g4 = as.vector(g4)
        probs.i[i] = sum(g4*pt.C2)
      }
      
    } # end of if(T>0)
    if((sum(J.t)==0)|(T==0)) probs.i[i]=1  
    
  } # end of site loop (i)
  
  q = -1*sum(log(probs.i))
  if(!is.finite(q)) return(1000000000)
  return(q)
}


ests = function(maxlikes,input.hess,migration="none", n, X, Z, T=length(n.it[1,]), prior="poisson"){
  #David Dail, revised March 2, 2010
  #
  #function to calculate the estimated values of N.t's (total abundances) for every sampling period
  #Will also provide the asymptotic SE's, using the Delta method
  #
  #maxlikes is the estimated MLE values of the parameters
  #input.hess is the hessian matrix of the likelihood evaluated at the MLEs (provided by optim())
  #migration="none" means r=0, w=1, closed popn assumption
  #migration="constant" means constant r (ie, log(r) = r0)
  #migration="autoreg" means r depends on N.it-1 (ie, log(r) = r0 + log(N.it-1)
  #if migration type is reshuffle or trend1, use migration="none" here.
  #
  #T is the number of primary sampling occasions
  #X and Z are the covariate matrices for lambda.i and p.it, respectively (constructed the same way
  # as in nmix.mig()
  #
  #prior="poisson" or "NB" (for negative binomial)
  
  X=as.matrix(X)
  Z=as.matrix(Z)
  betas = length(X[1,])
  phis = length(Z[1,])
  n.it=n
  R=length(n.it[,1])
  lam.i = exp(X%*%maxlikes[1:length(X[1,])])
  p.it.vec = expit(Z%*%maxlikes[(1+length(X[1,])):(length(X[1,])+length(Z[1,]))])
  p.it = matrix(p.it.vec,nrow=R,ncol=length(n.it[1,]),byrow=TRUE)
  ave.pit = mean(na.omit(p.it))
  
  if(prior=="NB") disp=exp(maxlikes[length(maxlikes)])
  
  if(migration=="none"){ #if there is not migration... a little simpler.
    inv.hess=abs(solve(input.hess)[1:betas,1:betas])
    lambda = lam.i
    par.ests = c(mean(lambda),ave.pit,0,1)
    N.hat = sum(lambda)
    partial.f.1 = N.hat
    partial.f = partial.f.1
    if(betas>1) {
      partial.f2.i = numeric(R)
      for (i in 1:R) {
        partial.f2.i[i] = X[i,2] * exp(X[i,]%*%maxlikes[1:betas])
      }
      partial.f2 = sum(partial.f2.i)
      partial.f = c(partial.f, partial.f2)
    }
    if(betas>2) {
      partial.f3.i = numeric(R)
      for (i in 1:R) {
        partial.f3.i[i] = X[i,3] * exp(X[i,]%*%maxlikes[1:betas])
        
      }
      partial.f3 = sum(partial.f3.i)
      partial.f = c(partial.f, partial.f3)
    }
    k.b=4
    while(k.b<=betas) {
      partial.f3.i = numeric(R)
      for (i in 1:R) {
        partial.f3.i[i] = X[i,k.b] * exp(X[i,]%*%maxlikes[1:betas])
        
      }
      partial.f3 = sum(partial.f3.i)
      partial.f = c(partial.f, partial.f3)
      k.b=k.b+1
    }
    
    
    
    var.Nhat = t(partial.f)%*%inv.hess%*%(partial.f)
    
    Nt = N.hat
    var.Nt = var.Nhat
    se.Nt = sqrt(abs(var.Nhat))
    N.ave = N.hat
    var.Nave = var.Nhat/T
    sd.Nave=sqrt(abs(var.Nave))
    
    trend=0
    var.trend = NA
    sd.trend=NA
    N.t = rep(Nt, T)
    seN.t = rep(se.Nt,T)
  }
  
  if(migration!="none") {  #if there is migration...
    Beta = maxlikes[1:betas]
    inv.hess=solve(input.hess)[c(1:betas,(betas+phis+1),(betas+phis+2)),c(1:betas,(betas+phis+1),(betas+phis+2))]
    
    w=expit(maxlikes[betas+phis+2])
    r=exp(maxlikes[betas+phis+1])
    lambda=lam.i
    par.ests = c(mean(lambda),ave.pit,r,w)
    ##hess and f.t need to be in this order: B0, B1, B2, ..., Bp, r, w
    
    N.t = numeric(T)
    seN.t = numeric(T)
    for(t in 1:T){   
      
      if( migration=="constant"){
        Nt.i = lambda*(w^(t-1)) + rep( r*(w^(t-1)-1)/(w-1), R)
        Nt = sum(Nt.i)
        St.i = lambda*((1-w^t)/(1-w)) - rep((r/(1-w))*(1-w^t)/(1-w),R) + rep(r*t/(1-w),R) 
        St = sum(St.i)
        N.ave = 1/t*St
        trend = rep(w,R) + rep(r,R)/lambda
        f.Nt.i = matrix(0,nrow=R,ncol=(betas+2))
        f.St.i = matrix(0,nrow=R,ncol=(betas+2))
        f.trend.i = matrix(0,nrow=R,ncol=(betas+2))
        
        for(i in 1:R){
          f.Nt.i[i,1]=lambda[i]*w^(t-1)
          if(betas>1) f.Nt.i[i,2] = X[i,2]* lambda[i]*w^(t-1)
          if(betas>2) f.Nt.i[i,3] = X[i,3]* lambda[i]*w^(t-1)
          k.b=4
          while(k.b<=betas) {
            f.Nt.i[i,k.b] = X[i,k.b]*lambda[i]*w^(t-1)
            k.b = k.b+1
          }
          
          f.Nt.i[i,betas+1] = r*(w^(t-1)-1)/(w-1)
          f.Nt.i[i,betas+2] = (lambda[i]*(w^(t-1)*(t-1))*(w-w^2)/w + 
                                 r*w^(t-1)*(t-1)*(w-w^2)/(w*(w-1)) - r*(w^(t-1)-1)*(w-w^2)/(w-1)^2 )
          
          f.St.i[i,1]=lambda[i]*(1-w^t)/(1-w)
          if(betas>1) f.St.i[i,2] = X[i,2]*lambda[i]*(1-w^t)/(1-w)
          if(betas>2) f.St.i[i,3] = X[i,3]*lambda[i]*(1-w^t)/(1-w)
          k.b=4
          while(k.b<=betas) {
            f.St.i[i,k.b] = X[i,k.b]*lambda[i]*(1-w^t)/(1-w)
            k.b = k.b+1
          }
          
          f.St.i[i,betas+1] = -1*r*(1-w^t)/(1-w)^2 + r*t/(1-w)
          f.St.i[i,betas+2] = (-1*lambda[i]*w^t*t*(w-w^2)/(w*(1-w)) 
                               - lambda[i]*(1-w^t)*(-w+w^2)/(1-w)^2 + lambda[i]*w^t*t*(w-w^2)/(w*(1-w)^2) 
                               + 2*r*(1-w^t)*(-w+w^2)/(1-w)^3 - r*t*(-w+w^2)/(1-w)^2 )
          
          f.trend.i[i,1] = -r/(lambda[i]) 
          if(betas>1) f.trend.i[i,2] = -r*X[i,2]/(lambda[i])
          if(betas>2) f.trend.i[i,3] = -r*X[i,3]/(lambda[i])
          k.b=4
          while(k.b<=betas) {
            f.trend.i[i,k.b] = -r*X[i,k.b]/(lambda[i])
            k.b = k.b+1
          }
          
          
          f.trend.i[i,betas+1] = r/(lambda[i])
          f.trend.i[i,betas+2] = w-w^2
        } # end of site= 1,...,R loop
        
        f.Nt = colSums(f.Nt.i)
        f.St = colSums(f.St.i)
        #f.trend = colSums(f.trend.i)
        
        var.Nt = t(f.Nt) %*% inv.hess %*% (f.Nt)
        var.St = t(f.St) %*% inv.hess %*% (f.St)
        var.Nave = (1/t)^2*var.St
        se.Nt = sqrt(abs(var.Nt))
        se.Nave = sqrt(abs(var.Nave))
        trend.i = rep(w,R) + rep(r,R)/lambda
        var.trend.i = numeric(R)
        
        for(i in 1:R){
          var.trend.i[i] = t(f.trend.i[i,])%*%inv.hess%*%(f.trend.i[i,])
        }
        trend=sum(((trend.i-rep(1,R))^2)/var.trend.i)
        
        trend.i.se=sqrt(abs(var.trend.i))
        
      } #end of migration="constant"
      
      
      if( migration=="autoreg") {
        Nt.i = lambda*(w+r)^(t-1)
        Nt = sum(Nt.i)
        St.i = lambda*(((w+r)^t-1)/(w+r-1))
        St = sum(St.i)
        trend = w+r
        
        N.ave = 1/t*St
        f.Nt.i = matrix(0,nrow=R,ncol=(betas+2))
        f.St.i = matrix(0,nrow=R,ncol=(betas+2))
        for(i in 1:R){
          f.Nt.i[i,1]=lambda[i]*(w+r)^(t-1)
          if(betas>1) f.Nt.i[i,2] = X[i,2]* lambda[i]*(w+r)^(t-1)
          if(betas>2) f.Nt.i[i,3] = X[i,3]* lambda[i]*(w+r)^(t-1)
          k.b=4
          while(k.b<=betas) {
            f.Nt.i[i,k.b] = X[i,k.b]* lambda[i]*(w+r)^(t-1)
            k.b = k.b+1
          }
          
          f.Nt.i[i,betas+1] = lambda[i]*(w+r)^(t-1)*(t-1)*r/(w+r)
          f.Nt.i[i,betas+2] = lambda[i]*(w+r)^(t-1)*(t-1)*(w-w^2)/(w+r)
          
          f.St.i[i,1]=lambda[i]*(((w+r)^t )-1)/(w+r-1)
          if(betas>1) f.St.i[i,2] = X[i,2]*lambda[i]*(((w+r)^t )-1)/(w+r-1)
          if(betas>2) f.St.i[i,3] = X[i,3]*lambda[i]*(((w+r)^t )-1)/(w+r-1)
          k.b=4
          while(k.b<=betas) {
            f.St.i[i,k.b] = X[i,k.b]* lambda[i]*(((w+r)^t )-1)/(w+r-1)
            k.b = k.b+1
          }
          
          f.St.i[i,betas+1] =( lambda[i]*(w+r)^t*t*r/(w+r)/
                                 ( w+r-1) - (lambda[i]*(w+r)^t - lambda[i])*r/(w+r-1)^2)
          
          f.St.i[i,betas+2] = (lambda[i]*(w+r)^t*t*(w-w^2)/(w+r)/(w+r-1) - 
                                 (lambda[i]*(w+r)^t - lambda[i])*(w-w^2)/ (w+r-1)^2 )
        } #end of site=1,...,R loop
        f.Nt = colSums(f.Nt.i)
        f.St = colSums(f.St.i)
        
        
        var.Nt = t(f.Nt) %*% inv.hess %*% (f.Nt)
        se.Nt = sqrt(abs(var.Nt))
        var.St = t(f.St) %*% inv.hess %*% (f.St)
        var.Nave = (1/t)^2*var.St
        se.Nave = sqrt(abs(var.Nave))
        f.trend = c(rep(0,betas), r, (w-w^2))
        
        var.trend=t(f.trend)%*%inv.hess%*%(f.trend)
        se.trend = sqrt(abs(var.trend))
      } #end of migration="autoreg"
      
      N.t[t] = Nt
      seN.t[t] = se.Nt
    } # end of loop through t= 1,...,T
    
  } # end of if migration!="none" 
  
  #could slightly modify to return N.ave, se.Nave, trend, se.trend also.
  
  ans.ret=c(N.t,seN.t,par.ests)
  return(ans.ret)
  
}

## Below is the Mallard data set that was originally published in the online Supplement
# of Kery, Royle, Schmid, Modeling avian abundance from replicated counts using 
# binomial mixture models, Ecological Applications, 2005.

"mallard.data" <-
  structure(c(0, 0, 3, 0, 3, 0, 0, 0, 0, 0, 0, NA, 0, 0, 0, 0, 
              0, 1, 0, 1, 0, 0, 0, 0, 0, 10, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 
              0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
              0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, NA, 4, 0, 0, 0, 0, 0, 0, 0, 
              0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 
              0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, NA, 
              0, 0, 0, 0, 1, 4, 0, 0, 0, 0, 0, 0, 1, 0, 3, 0, 0, 0, 2, 1, 0, 
              0, 0, 1, 0, 0, 0, NA, 0, 0, 2, 0, 0, 0, 0, 4, 0, 0, 0, 0, 0, 
              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 1, 0, 0, 
              0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 
              1, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 
              0, 0, 0, 0, 0, 0, 0, NA, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 
              0, 12, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
              0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 
              0, 0, 0, NA, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 
              NA, 0, 0, 0, NA, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 
              0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, NA, 0, 0, 0, 0, 1, 4, 0, 
              0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 2, 1, 0, 0, 0, 0, 0, 0, 0, NA, 
              0, 0, 1, 0, 0, 0, 0, 3, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
              0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 
              0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 2, 0, 0, 0, 0, 0, 0, 0, 
              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 3, 0, 0, 0, 0, 0, 0, NA, 
              0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 7, 0, 0, 0, 0, 1, 0, 0, 
              0, 0, 0, NA, 0, 0, NA, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, NA, 
              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, NA, 1, 0, 0, 0, NA, 4, NA, NA, 
              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, NA, 0, 0, 0, NA, 0, 0, 
              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, NA, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
              1, 0, 0, 0, 0, NA, 0, 0, 0, 0, 0, 0, 0, 0, NA, 0, 0, 0, 1, 0, 
              1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, NA, NA, 0, 0, NA, 0, 0, 
              0, 3, NA, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, NA, 0, 0, 0, 0, 0, 
              0, 0, 0, 0, 3, 1, 0, 0, 0, 2, NA, 0, 0, 0, 0, 0, 0, 0, 0, 1, 
              0, 0, 0, 1, 0, 0, NA, 1, 0, 2, 0, NA, NA, NA, NA, 0, NA, NA, 
              NA, 0, 1, 0, 0, 0, NA, 0, 0, NA, NA, NA, NA, NA, NA, NA, NA, 
              NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, -0.506, -0.934, 
              -1.136, -0.819, 0.638, -1.329, -1.448, -0.321, -0.231, -1.097, 
              -0.224, NA, 0.417, -1.354, -1.117, -0.278, -0.483, -1.07, -1.115, 
              -0.862, -0.242, -0.462, -0.358, 0.417, -0.458, 1.866, -0.264, 
              0.427, -0.42, 0.952, 0.904, -1.697, -1.067, -0.837, 0.568, 0.223, 
              0.977, 0.952, 0.395, 0.338, -0.727, -0.819, 0.24, -0.139, -1.508, 
              -0.101, -0.539, 0.897, -0.374, -0.308, -1.199, -0.984, -0.042, 
              0.073, -1.097, 0.334, 0.662, -0.206, -1.568, 0.568, -1.097, -0.394, 
              -1.212, 0.919, -0.019, -0.123, -0.735, -0.672, -0.476, 0.429, 
              0.666, -0.417, -0.803, -0.172, -0.785, -0.792, -0.795, 0.647, 
              -0.308, 0.298, 1.382, -0.394, -0.669, 0.476, -0.941, -0.601, 
              0.735, 1.076, -0.212, -1.513, -0.958, 0.516, -0.692, -0.172, 
              0.09, -0.731, 0.481, -1.428, -0.004, -0.443, -0.394, -0.366, 
              -0.312, 0.152, 0.417, 0.558, -0.194, -0.431, -0.348, 0.984, -1.014, 
              -1.148, -0.605, 0.88, -0.554, -0.355, -0.524, 0.355, -0.729, 
              -1.393, -0.01, 0.006, -0.997, -1.108, -0.513, 1.973, 3.066, -0.413, 
              -0.834, -1.377, -0.647, -0.984, -0.681, -1.148, -1.014, 0.104, 
              1.168, -1.122, 0.557, -1.097, -0.652, 0.568, -0.469, -1.088, 
              -0.086, NA, -0.735, -0.042, -1.249, -0.373, 1.013, -0.081, 0.342, 
              -0.451, 1.225, -0.862, -1.196, -0.482, 0.806, 0.806, 0.735, -0.06, 
              0.467, -1.485, 1.502, -0.264, 0.615, 3.066, -0.113, 2.119, 0.099, 
              0.152, -0.35, 0.318, -0.51, 0.436, -0.086, -0.362, 0.502, 0.768, 
              2.759, -0.451, 1.102, 0.066, 1.289, 0.028, 0.994, 1.734, 0.443, 
              -0.139, 0.735, -0.448, 0.4, 0.659, 1.147, -0.194, 0.21, -0.745, 
              -0.272, -0.344, -0.897, -0.819, 1.422, -0.643, 2.355, 1.379, 
              -0.072, 1.193, 1.609, -0.291, 1.982, 0.323, -1.356, 0.107, -0.036, 
              1.039, 0.493, 1.28, 0.201, 2.094, -0.203, 1.71, -0.506, 0.152, 
              1.329, 0.251, -0.913, 0.092, 0.735, -1.228, 5.355, 0.066, 5.494, 
              -0.378, -0.618, 0.264, -0.126, 3.713, 1.056, -0.506, -0.991, 
              -1.339, -0.927, 0.88, -1.042, -1.562, -0.557, -0.231, -1.021, 
              0.058, NA, 0.284, -1.014, -0.224, -0.182, -0.884, -1.258, -0.735, 
              -1.115, -0.242, -0.615, -0.212, 0.682, -1.542, 1.883, 0.568, 
              0.702, -0.264, 1.409, 0.913, -1.753, -1.543, -0.472, 0.152, -0.63, 
              0.042, 1.123, -0.01, 0.586, -1.097, -0.738, 0.24, -0.139, -1.644, 
              -0.228, -0.688, 0.491, -0.469, -0.308, -0.904, -0.589, -0.237, 
              0.546, -1.216, 0.425, 0.735, 0.101, -1.52, 0.568, -1.051, -0.394, 
              -1.274, 0.459, -0.362, 0.207, -0.228, -0.672, -0.61, 0.984, 0.666, 
              -0.985, -0.903, -0.172, -0.889, -0.669, -1.014, 0.592, 0.305, 
              0.298, 0.605, -0.212, -0.546, 0.476, -0.941, NA, 0.56, 0.863, 
              0.334, NA, -1.143, 0.759, -0.615, 0.368, 0.214, -0.466, 0.199, 
              -1.132, 0.256, -0.264, 0.061, -1.143, -0.047, 0.152, 0.735, 0.965, 
              -0.293, -0.166, 0.035, 0.984, -1.055, -1.058, -0.264, 1.245, 
              -0.025, -0.482, -0.472, 0.694, -0.729, -1.498, -0.01, 0.443, 
              -0.505, -0.872, -0.666, 1.791, 3.066, -0.283, -0.655, -1.377, 
              0.058, -0.539, -1.027, -1.148, -0.664, 0.104, 1.507, -0.334, 
              0.314, -0.889, -0.652, 0.984, -0.182, -0.716, -0.026, NA, -0.482, 
              -0.042, -1.417, -0.314, 1.145, -0.373, 0.215, -0.551, 1.225, 
              -0.862, -0.935, -0.101, 0.806, 0.806, 0.851, -0.06, 0.073, -1.245, 
              1.502, -0.264, 0.615, 3.267, -0.378, 1.609, 0.099, 0.298, -0.652, 
              0.318, 0.086, 0.01, -0.264, -0.362, 0.385, 0.6, 1.379, -0.2, 
              0.912, -0.019, 1.289, -0.902, 1.123, 1.234, 0.268, 0.152, 0.735, 
              0.323, 1.268, 0.532, 0.507, -0.095, 0.385, -0.488, -0.855, -0.282, 
              -1.208, -1.305, 1.497, -0.908, 2.355, 0.919, -0.072, 1.262, 1.344, 
              0.025, 2.185, 0.323, -1.048, -0.117, 0.058, 1.039, 0.427, 0.908, 
              0.201, 2.094, -0.203, 1.71, -0.882, 0.006, 1.161, 0.547, -0.577, 
              0.092, -0.348, -0.615, 5.98, 0.323, 5.494, 0.615, -0.453, 0.264, 
              0.152, 3.713, 2.262, -0.506, -1.162, -1.61, -1.197, 1.042, -0.899, 
              -1.676, -0.636, -0.001, -0.832, -0.224, NA, 0.549, -1.159, -0.788, 
              0.009, -0.819, -1.211, -0.482, -0.988, -0.399, -0.462, 0.735, 
              0.639, -1.434, 0.923, -0.264, 0.482, -0.108, 1.237, 1.054, -1.585, 
              -1.437, -0.639, 0.152, -0.203, NA, 1.18, 0.152, NA, -1.004, -0.981, 
              0.682, 0.298, -1.474, -0.101, -0.737, 0.491, -0.66, -0.308, -1.157, 
              -0.836, NA, 0.231, -1.037, 0.789, 0.516, -0.001, -1.615, 0.568, 
              -1.189, -0.212, -1.491, NA, -0.819, 0.702, -0.228, -0.672, -0.52, 
              0.984, NA, NA, -0.903, 0.152, -1.097, -0.546, -0.722, 0.207, 
              -0.104, 0.006, 0.346, -0.394, -0.423, 0.638, -0.941, NA, 0.035, 
              1.076, -0.03, NA, -0.264, 0.82, -0.615, -0.064, 0.214, -0.554, 
              0.622, -1.28, 0.672, 0.271, -0.394, -1.273, NA, -0.2, -0.007, 
              0.558, -0.144, -0.219, -0.098, 0.984, -1.014, -1.148, -0.491, 
              1.245, 0.24, -0.608, -0.837, 1.033, -0.322, -1.217, 0.152, 0.37, 
              -0.423, -0.99, -0.615, 1.609, NA, 0.108, -0.852, -1.424, 0.246, 
              -0.539, -0.819, -1.148, -0.664, 0.582, 1.101, -0.161, 0.88, -0.847, 
              -0.652, 2.649, -0.182, -0.902, 0.509, NA, NA, 0.152, -1.361, 
              NA, 1.278, -0.442, 0.659, -0.551, NA, -1.115, -1.153, -0.101, 
              0.806, 0.806, 1.201, 0.099, -0.084, -1.325, 1.076, -0.383, NA, 
              2.965, -0.554, 1.609, 0.258, 0.103, -0.602, 0.735, -0.179, 0.436, 
              -0.443, -0.362, 0.385, 1.105, 1.839, -0.099, NA, -0.448, 1.502, 
              -0.282, 1.123, 1.234, 0.56, -0.431, 0.516, -0.448, 1.144, 0.279, 
              1.289, 0.053, 0.735, -0.328, NA, -0.282, -0.936, -1.062, 1.557, 
              NA, NA, NA, NA, 1.331, NA, NA, NA, 0.323, -1.528, 0.017, 0.105, 
              1.039, NA, 0.747, 0.498, NA, NA, NA, NA, NA, NA, NA, NA, NA, 
              NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, -1.761, -2.904, -1.69, 
              -2.19, -1.833, -2.619, -2.69, -2.119, -2.047, -2.333, -1.69, 
              -1.047, -1.833, -0.904, -2.547, -1.333, -1.833, -2.333, -2.404, 
              -1.333, 0.024, -1.69, -1.976, -0.476, -2.404, -2.619, -2.761, 
              -1.19, -0.619, -1.047, -2.19, -2.904, -1.761, -1.333, -2.119, 
              -0.761, 1.096, 0.453, -2.119, 1.596, -2.19, -0.476, -1.833, -2.333, 
              -1.833, 0.096, -1.833, -1.904, -1.833, -1.976, -2.119, -2.19, 
              -0.547, -2.619, -1.69, -2.333, -2.404, -0.476, -2.69, -2.904, 
              -2.404, -1.976, -2.404, 2.667, -2.261, -2.333, -2.404, -2.404, 
              -2.547, -1.833, 1.596, 2.453, -1.904, -1.976, -2.047, -2.261, 
              0.739, -1.833, -0.261, -2.19, -1.833, -2.904, -2.69, -2.404, 
              -1.833, 0.453, -1.904, 0.524, -2.404, -2.904, -2.761, -2.333, 
              -2.19, -2.404, -0.19, -2.261, -2.261, -1.904, -0.261, -0.976, 
              -2.19, -2.404, 1.31, -0.761, -0.261, -0.404, -2.333, -2.547, 
              -1.976, -0.833, -1.619, -1.19, -2.333, 1.239, -1.976, -2.333, 
              -1.976, -2.333, -2.333, -2.333, -0.333, -2.333, -1.833, -2.333, 
              -2.547, 2.167, 1.881, -2.547, -2.69, -1.976, -2.69, 1.096, -2.761, 
              -2.761, -2.404, -1.976, -2.261, -2.547, 0.596, -2.69, -0.904, 
              -1.119, -2.404, -2.333, -2.833, -2.333, 0.739, -1.047, -1.976, 
              -0.761, 1.167, -0.19, -2.69, -1.904, 0.739, -2.547, -2.333, -1.833, 
              -1.761, -1.761, -2.404, -2.404, -0.547, -2.761, 0.953, -2.333, 
              -0.404, 0.096, -0.476, -0.547, -1.69, -2.333, -1.619, -2.119, 
              -2.19, -2.261, -2.69, -1.619, -1.904, -2.904, 0.524, -2.333, 
              1.596, -0.619, -1.047, -0.976, -1.833, -0.476, -2.904, -2.904, 
              -1.833, -1.833, -0.047, -0.69, -1.619, -2.19, -2.904, -2.119, 
              1.739, -1.833, -1.833, -2.19, -1.119, 0.667, 0.667, 0.667, 0.596, 
              -1.976, 1.596, 1.667, 2.381, 0.667, -0.19, 0.453, 0.381, 0.096, 
              1.453, -0.261, 0.81, 2.096, 0.453, 1.596, 0.167, 2.381, 2.31, 
              2.096, 0.596, 1.81, 1.524, 2.381, 1.381, 2.239, 2.096, 1.167, 
              0.667, 1.167, 0.667, 1.739, 1.81, 0.31, -1.047, -0.476, -0.69, 
              0.167, 0.167, -1.19, -0.476, -0.547, -1.119, 0.453, 0.739, -0.333, 
              1.096, -0.261, 0.31, -1.047, -0.976, -0.833, -0.261, 1.524, -0.261, 
              -0.476, 0.596, -0.976, -0.547, -0.619, 0.453, 0.524, 0.453, -0.619, 
              -0.976, -0.547, -0.261, -0.119, -0.19, 2.881, 1.381, -0.261, 
              2.096, -1.261, 1.31, -0.761, -0.619, -0.833, 0.667, -1.047, -0.261, 
              -0.119, -0.404, -0.619, -0.547, 0.81, -0.833, -1.19, -0.261, 
              -1.333, 1.096, -0.761, -1.047, -1.261, -0.976, 0.096, 3.239, 
              -0.619, -0.261, -0.976, -0.619, -1.19, -0.261, 2.667, 3.024, 
              -0.976, -0.476, -0.547, -0.261, 1.739, -0.904, 1.596, 0.453, 
              -0.761, -1.19, -1.261, -0.761, -0.833, 1.596, -0.833, 1.31, -0.976, 
              -0.904, -1.404, -0.547, -0.547, -1.19, 1.167, -1.69, -0.833, 
              -0.904, 1.096, -0.476, -1.19, -0.904, 3.239, 1.381, 0.381, 0.453, 
              -0.833, -1.619, -0.547, 0.524, -1.119, -0.976, -0.261, 1.596, 
              -0.19, -0.547, -0.547, -0.904, -0.976, -1.047, 2.667, -0.833, 
              -0.261, -0.404, -0.547, 3.096, 3.024, -1.19, -0.547, 0.31, -0.833, 
              1.596, -1.333, -1.047, -0.904, -0.619, -1.19, -0.547, 1.524, 
              -0.19, -0.476, -0.119, -0.19, -0.619, -0.976, -0.476, 2.167, 
              -0.476, -0.547, 0.096, 2.167, 1.31, -0.547, -0.619, 2.453, -0.261, 
              -1.047, -0.476, -1.119, -1.119, -0.833, -0.833, 1.239, -1.761, 
              2.381, -0.404, 1.596, 1.596, 0.596, -0.261, -0.619, 0.167, -0.261, 
              -0.547, -0.619, -0.476, -1.69, -1.047, 0.096, -0.261, 2.167, 
              -1.047, 3.024, 0.667, -0.261, 1.381, -0.619, 0.524, -1.047, -1.619, 
              -0.404, -0.261, 1.667, 0.881, -1.047, -1.19, -0.547, 0.024, 2.31, 
              0.096, -0.19, -0.904, -0.119, 1.667, 1.881, 1.596, 1.096, -0.904, 
              2.596, 2.239, 2.881, 1.453, 0.596, 1.31, 0.524, 1.024, 2.453, 
              0.881, 1.453, 3.239, 2.881, 2.167, 2.096, 3.453, 2.524, 3.096, 
              1.453, 2.167, 1.739, 3.239, 1.81, 2.667, 2.667, 2.167, 1.596, 
              3.167, 2.667, 2.81, 3.31, 1.381, 0.596, 1.453, 1.239, 1.381, 
              1.381, 1.596, 1.453, 1.167, -0.261, 1.453, 2.596, 1.953, 1.239, 
              1.31, 1.524, 0.524, 2.31, 0.667, 1.524, 3.239, 1.596, 1.239, 
              1.667, 0.524, 0.524, 0.881, 2.596, 2.31, 2.096, 0.31, 1.453, 
              -0.404, 1.31, 0.596, 0.596, NA, 2.167, 1.596, NA, 1.096, 1.953, 
              0.596, 1.453, 0.596, 1.596, 0.596, 2.596, 2.167, 1.596, 0.524, 
              1.524, NA, 0.596, -0.547, 0.524, 0.381, 1.667, 0.381, 0.096, 
              0.453, -0.476, 1.667, NA, 0.167, 0.667, 1.024, -0.19, 0.596, 
              0.667, NA, NA, 0.81, 1.31, 0.596, 1.524, 2.881, 0.167, 2.596, 
              1.453, -0.119, -0.476, 0.31, 0.524, 0.453, NA, 0.596, 1.667, 
              0.596, 1.096, 0.024, 1.381, 1.096, -0.261, 1.667, 0.381, 0.167, 
              1.096, 1.739, 0.739, 0.524, -0.119, NA, 3.167, 1.596, 1.524, 
              0.667, 0.524, 0.596, 1.667, 1.596, 0.096, 1.096, 1.953, 0.524, 
              1.453, 0.596, 0.667, 0.667, 0.596, 2.739, 0.596, 0.667, 0.667, 
              0.167, 3.524, NA, 0.881, -0.261, 0.81, 0.881, 2.096, 1.381, 0.524, 
              1.096, 1.31, 0.524, 0.739, 1.667, 0.81, 0.524, 1.524, 2.024, 
              0.881, 0.739, 0.667, NA, 1.096, 0.596, 1.31, 3.167, 2.167, 1.453, 
              0.596, NA, 2.167, 0.667, 0.667, -0.261, -0.261, 1.596, 1.167, 
              1.739, -0.119, 3.81, 0.596, NA, 3.381, 1.739, 0.596, 0.381, 1.31, 
              0.596, 0.453, -0.261, 1.453, -0.619, -0.833, 1.881, 0.453, 2.596, 
              -0.547, NA, 2.096, 0.453, 2.596, 0.096, 3.381, 0.096, -0.404, 
              0.667, 1.953, 2.953, 1.739, -0.476, 0.81, 0.596, 1.739, NA, 1.596, 
              0.667, 0.524, 1.31, NA, NA, NA, NA, 0.381, NA, NA, NA, 2.453, 
              1.31, 1.453, 1.381, 2.024, NA, 1.596, 2.953, NA, NA, NA, NA, 
              NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 
              -1.173, -1.127, -0.198, -0.105, -1.034, -0.848, -0.91, -1.003, 
              -0.058, -0.631, 0.066, -0.709, -1.111, -0.26, -0.647, 0.189, 
              0.143, -1.235, -0.972, -1.219, 0.313, -1.003, 0.019, 0.484, -0.399, 
              -0.863, -0.987, -0.554, 0.066, 0.004, -0.693, -0.925, -1.065, 
              0.143, -0.383, -1.034, 1.01, 1.18, -0.136, 1.289, -0.6, 0.855, 
              -0.089, -0.863, -0.91, 1.025, -0.383, -0.987, -1.173, 0.344, 
              -0.445, 0.437, 1.041, 0.375, -0.894, -0.972, -0.662, 0.948, -0.77, 
              -0.291, -1.018, -0.647, -1.25, 2.372, -0.817, -0.801, -0.693, 
              -1.111, -1.266, -1.436, 1.18, 1.923, 0.22, -0.6, -1.096, -1.188, 
              1.18, 0.607, 1.087, 0.174, 0.251, -0.461, -0.724, -0.848, -1.127, 
              1.83, 0.731, 1.35, -0.6, -0.693, -1.065, -0.523, -1.127, -0.26, 
              0.422, -0.507, 0.097, -0.941, 0.499, 0.375, -0.616, -1.096, 1.397, 
              1.304, 0.886, 0.7, -0.832, -1.003, -0.801, 0.282, 0.174, -0.786, 
              -1.08, 1.273, -0.585, -0.26, -0.77, -1.065, -1.204, -1.328, 0.344, 
              -0.91, -1.157, -1.188, -1.034, 1.165, 0.592, -0.739, -1.142, 
              -1.08, -1.065, 1.041, -1.065, -1.08, -1.142, -0.554, -0.786, 
              -0.879, 0.809, -0.43, 0.112, -0.414, -0.987, -1.127, -1.096, 
              -1.096, 1.087, 0.344, -1.096, -0.538, 1.196, 0.747, 0.05, -1.235, 
              1.025, -0.352, -0.461, -1.08, -0.972, -0.972, -1.297, -0.972, 
              0.22, -1.421, 1.289, -0.074, 1.134, 0.7, 0.561, 0.251, 0.004, 
              -0.538, -0.987, 0.36, 0.22, -0.445, -0.925, -1.049, 0.066, -0.724, 
              0.994, -1.018, 1.103, 0.731, 0.329, 0.329, 0.174, 0.979, -0.074, 
              -1.08, -0.213, -1.188, 1.32, 1.087, -0.352, -1.127, -0.352, -0.167, 
              1.242, 0.267, -0.987, -1.173, 0.375, 1.645, 1.242, 1.784, 1.35, 
              0.004, 1.289, 1.614, 1.753, 1.567, 0.917, 1.01, 1.072, 1.025, 
              1.923, 0.994, 1.041, 1.66, 1.66, 1.97, 1.815, 2.279, 2.434, 1.97, 
              1.196, 1.505, 1.041, 1.35, 1.041, 1.35, 1.505, 1.35, 1.196, 1.196, 
              1.66, 1.66, 1.815, 0.801, 0.115, -0.479, 0.315, -1.102, 0.741, 
              0.115, -1.007, -0.913, 1.556, -1.626, 1.647, -0.399, 0.685, 0.801, 
              0.741, -0.322, 0.801, 1.175, -0.245, -1.007, 0.503, -0.734, 0.381, 
              -0.479, -1.301, -0.563, 0.252, 0.444, 0.115, 0.801, 0.185, 0.381, 
              0.444, -0.563, -0.647, 0.252, 0.115, -1.102, -0.168, 0.857, -1.102, 
              -1.406, -0.734, 1.944, 1.175, 0.626, -0.479, 0.741, 0.503, 1.175, 
              0.626, -0.322, -1.007, -0.025, -1.514, -0.734, 0.503, 0.741, 
              -0.563, 0.857, -0.095, 2.255, -3.336, 0.115, 0.252, 1.175, -0.245, 
              0.965, -0.563, 0.115, -0.647, 0.566, 0.315, 1.86, 1.273, -0.734, 
              0.252, 0.503, -0.734, -0.322, -0.095, 1.273, 0.315, -0.095, 2.063, 
              0.045, -0.647, -0.095, 0.444, 0.857, -0.095, 1.51, 0.315, -0.168, 
              -1.406, 0.801, 0.626, 0.444, -0.025, -1.514, -0.322, -0.399, 
              0.566, 0.381, -0.479, 0.626, 0.381, -1.199, -1.199, 1.224, 0.965, 
              1.556, -0.095, -1.406, -0.245, 0.444, -0.479, -0.479, 1.818, 
              0.315, -0.734, 1.273, 1.416, 0.503, -0.095, -1.738, 1.07, 0.965, 
              0.741, 0.801, 0.626, -0.563, 0.965, 0.045, 0.741, -0.479, 0.315, 
              -1.102, 1.224, 0.566, -0.563, 0.741, -0.168, -0.025, 1.322, 1.175, 
              -0.322, 0.185, 0.045, -0.399, 0.045, -0.245, 0.566, -0.913, 1.175, 
              1.07, -0.245, -0.025, -0.025, 0.045, 0.381, -1.007, 1.371, -0.647, 
              -0.025, -0.399, -1.857, -1.406, -0.734, 0.381, 0.685, 0.566, 
              -1.199, -0.399, -0.647, -0.025, 0.115, 0.045, 0.185, -0.913, 
              0.566, -0.245, -1.301, -0.647, -0.168, -0.322, -1.199, 0.045, 
              0.685, -0.734, -1.301, -0.168, -0.245, -0.647, 0.626, 0.045, 
              2.14, 0.381, -0.168, 1.465, 1.322, -0.822, 1.017, -0.647, -0.913, 
              -0.822, -0.563, -0.399, -0.245, -0.479, 1.122, 1.902, 0.965, 
              0.801, -0.245, 0.252, -0.168, 0.626, -1.102, -0.647, -0.479, 
              2.217, -3.158, 0.185, 0.626, 0.185, -0.025, -1.199, 0.503, -4.406, 
              -1.301, -4.945, -0.399, 0.252, 0.185, -0.563, -2.109, -1.857, 
              -1.156, -0.501, -0.101, 0.008, -1.193, 0.917, -1.083, -0.792, 
              0.553, 0.808, 1.79, -0.72, -0.647, 0.844, -0.32, -0.101, 1.135, 
              0.372, -0.138, -0.974, 1.572, 0.081, -0.611, 0.699, -1.047, -1.156, 
              -0.283, 0.772, 1.135, 0.299, 0.19, -1.083, -1.156, 0.626, 0.117, 
              0.517, -0.72, 0.99, 1.717, -1.265, -1.083, 0.844, 2.299, -0.065, 
              0.663, -0.647, -0.683, 0.917, -1.193, 0.626, 0.844, 0.553, -0.829, 
              -0.21, -0.683, 1.754, 1.608, 0.59, -0.974, 0.372, -1.083, 1.063, 
              -0.756, -1.265, -1.265, -0.174, -0.647, 1.099, -0.174, -1.265, 
              -1.083, -1.265, 1.135, -0.974, 1.463, 1.208, 1.354, 2.19, 0.808, 
              1.063, 0.335, 0.335, -0.029, -0.501, -0.429, -1.229, 0.299, -0.174, 
              -0.32, -1.011, 1.899, 0.517, 0.153, 2.117, 0.481, -0.865, 1.39, 
              -0.829, 0.262, -0.32, -0.792, -1.193, -0.72, 0.699, 1.426, 0.372, 
              0.954, -0.065, 0.699, 1.245, 0.444, -1.083, 0.044, 0.917, 1.354, 
              -0.21, -1.011, -1.011, -0.792, -1.193, -0.683, 0.99, -0.938, 
              0.081, 1.317, 0.226, -0.574, -0.392, -0.902, -1.265, 1.281, 0.372, 
              -0.611, -0.538, -1.156, -0.538, 0.044, -0.574, 0.008, 1.899, 
              0.59, -0.501, -0.865, -0.21, 0.735, 0.699, -0.101, -0.101, -1.156, 
              0.844, 0.044, 1.863, 1.463, -0.429, -1.265, 1.317, 1.245, -1.047, 
              1.172, 1.172, 0.262, 1.354, 1.608, 0.444, -0.32, 0.772, -0.611, 
              -0.065, -0.501, 0.153, 1.245, 0.262, 2.299, 0.626, 2.008, 0.553, 
              0.699, -0.938, 0.59, 0.008, 0.808, -0.611, -0.829, -0.029, 1.208, 
              1.026, 0.481, -0.065, 2.008, -0.065, 0.081, -0.974, 0.044, 1.135, 
              0.772, -1.265, 0.226, 1.026, -1.265, 0.008, -0.756, -1.229, -0.283, 
              -1.265, -1.047, -1.265, -1.229, 0.844, -0.938, -0.756, -1.193, 
              -0.174, 2.117, 1.717, 1.681, 1.717, -1.265, 1.245, 0.663, -1.265, 
              -1.265, -1.265, -1.265, -1.265, -1.265, -1.265, -1.265, -1.265, 
              -1.229, -1.265, -1.265, -0.72, -1.265, -1.265, -1.193, -1.265, 
              -1.265, -1.265, -1.265), .Dim = as.integer(c(239, 12)), .Dimnames = list(
                NULL, c("count1", "count2", "count3", "ivel1", "ivel2", "ivel3", 
                        "date1", "date2", "date3", "elev", "length", "forest")))


