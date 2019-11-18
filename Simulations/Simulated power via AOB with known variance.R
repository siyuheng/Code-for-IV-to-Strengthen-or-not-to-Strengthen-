library(mvtnorm)
library(nbpMatching)

#############
##Functions##
#############


# Paul's functions
# rank based Mahalanobis distance between each pair
smahal=function(X){
  X<-as.matrix(X)
  n<-dim(X)[1]
  k<-dim(X)[2]
  for (j in 1:k) X[,j]<-rank(X[,j])
  cv<-cov(X)
  vuntied<-var(1:n)
  rat<-sqrt(vuntied/diag(cv))
  cv<-diag(rat)%*%cv%*%diag(rat)
  out<-matrix(NA,n,n)
  library(MASS)
  icov<-ginv(cv)
  for (i in 1:n) out[i,]<-mahalanobis(X,X[i,],icov,inverted=T)
  out
}


#Calipers
calipers <- function(distmat,variable,tolerance=0.2){
  sd.distmat <- median(sd(distmat))
  sd.var <- sd(variable)
  penalty <- matrix(variable,nrow=length(variable),ncol=length(variable))
  penalty <- penalty - t(penalty)
  distmat[abs(penalty)> tolerance*sd.var]  <- distmat[abs(penalty)> tolerance*sd.var]  + abs((penalty/sd.var)[abs(penalty)> tolerance*sd.var]*.5*sd.distmat)
  distmat
}

#Wrap function for matching function.  Creates sinks.
matches <- function(select1row,X,sinks,dif.travel.time,mindist){
  cols.w.var <-  which(apply(X[select1row,],2,var)>0)        #Only give the columns with variation, otherwise "smahal" freaks out
  distmat=smahal(X[select1row,cols.w.var]);                  #Create distance matrix
  distmat=round(distmat*1000);              #Make into integer distance and add constant so sinks are smallest.
  
  #Add calipers
  #distmat <- calipers(distmat,X[select1row,"conditionI"],tolerance=0.01)
  #distmat <- calipers(distmat,X[select1row,"SES.missing"],tolerance=0.1)
  #distmat <- calipers(distmat,X[select1row,"income"],tolerance=0.2)
  #distmat <- calipers(distmat,X[select1row,"black"],tolerance=0.1)
  
  #Make sure we don't match babies who have similar encouragement.
  distmat.adjust <- min(sd(distmat))
  dtt <- dif.travel.time[select1row]              #My fingers are getting tired.
  store <- matrix(dtt,nrow=length(dtt),ncol=length(dtt))
  store <- abs(store - t(store))
  distmat[store<mindist]  <-   distmat[store<mindist]  + ifelse((max(store[store<mindist])-store[store<mindist])>2,2,(max(store[store<mindist])-store[store<mindist]))*distmat.adjust
  
  #Create artifical sinks
  size <- dim(distmat)[2]
  num.sinks <- size*sinks                                  #The "sinks" variable tells you what fraction of the babies should be matched to sinks - see line 14 above
  num.sinks <- 2*ceiling((size+num.sinks)/2) - size        #Make sure we have an even number of (babies + sinks) to match
  total <- size + num.sinks
  distmat <- cbind(rbind(distmat,matrix(0,nrow=num.sinks,ncol=size)),matrix(0,nrow=total,ncol=num.sinks));     #Add all of the sinks
  if(num.sinks>0){distmat[(size+1):total,(size+1):total] <- max(distmat)*3}                               #If there are sinks, make sure sinks don't match to sinks
  
  #Do the matching
  distmat=round(distmat);              #Make into integer distance 
  distmat=distancematrix(distmat)
  matching=nonbimatch(distmat)$matches;
  dif.travel.time.orig=dtt[matching$Group1.Row]
  dif.travel.time.matched=dtt[matching$Group2.Row]
  encouraged=dif.travel.time.orig<dif.travel.time.matched;                                                   #There are going to be NAs because of the sinks
  encouraged[(size+1):total]=FALSE;
  encouraged[is.na(dif.travel.time.matched)]=FALSE;                                                          #The NAs tell you that there's been a baby matched to a sink.
  
  #Some diff.travel.time will be the same within a pair.  Chose the first half to be encouraged and second to be unencouraged.  If odd number just drop last pair.
  same.diff <- sum(dif.travel.time.orig==dif.travel.time.matched,na.rm=TRUE)
  dif.travel.time.matched[is.na(dif.travel.time.matched)] <- dif.travel.time.orig[is.na(dif.travel.time.matched)]+1
  if(same.diff>0){
    encouraged[which(dif.travel.time.orig==dif.travel.time.matched)] <- sample(c(rep(TRUE,floor(same.diff/2)),rep(FALSE,ceiling(same.diff/2))),size=same.diff)
  }
  
  encouraged.numbers=matching$Group1.Row[which(encouraged)];
  encouraged.match=matching$Group2.Row[which(encouraged)];
  return(cbind(encouraged.numbers,encouraged.match))
}

expit <- function(x) return(exp(x)/(1+exp(x)))


# two endpoints of SI
low_high_endpoint <- function(dt, delta, lambda_0, lambda_1){
  
  I=nrow(dt)
  
  estimate_bias=rep(0,5)
  
  compliance = mean(dt$D.e) - mean(dt$D.u)
  
  for (j in 1:5){
    E_treat_U = mean(rbinom(I, 1, prob = expit(lambda_0 + lambda_1*dt$Z.e)))
    E_control_U = mean(rbinom(I, 1, prob = expit(lambda_0 + lambda_1*dt$Z.u)))
    #cat(E_treat_U, E_control_U, mean(expit(lambda_0 + lambda_1*dt$Z.e)), 
    #    mean(expit(lambda_0 + lambda_1*dt$Z.u)), '\n')
    estimate_bias[j] = delta*abs((E_treat_U - E_control_U)/compliance)
  }

  wald = (mean(dt$Y.e) - mean(dt$Y.u))/compliance
  bias = mean(estimate_bias)
  var_within=2/(I*(compliance)^2)
  var_across=(1+1/5)*var(estimate_bias)
  sd=sqrt(var_within+var_across)
  v=(5-1)*(1+var_within/var_across)^2
  return(c(wald - abs(bias)-qt(0.975, df = v)*sd, wald + abs(bias) + qt(0.975, df = v)*sd, bias, var_within, var_across, sd))
}


# Solve for lambda_0 for a given tau and lambda_1
solve_lambda_0 <- function(tau, lambda_1){
  func <- function(x, lambda_0, lambda_1) expit(lambda_0 + lambda_1*x)*(1/sqrt(2*pi))*exp(-x^2)
  func_2 <- function(lambda_0, lambda_1, tau) integrate(func, lower = 0, upper = 10, lambda_0 = lambda_0, lambda_1 = lambda_1)$value - 
    integrate(func, lower = -10, upper = 0, lambda_0 = lambda_0, lambda_1 = lambda_1)$value - 2*tau
  
  tryCatch({uniroot(func_2, c(0, 10), extendInt = c('yes'), lambda_1 = lambda_1, tau = tau)$root},
           error=function(error_message) {
             message(" Error in integrate(func, lower = 0, upper = 10, lambda_0 = lambda_0, lambda_1 = lambda_1) : 
                     non-finite function value ")
             return(NA)
           })
}

power_AOB_varytau<-function(sinks,mindist,delta,n,S,beta,lambda_1,tau_max,lambda_0_D,lambda_1_D,Mean,SD){
  count_1=0
  count_2=0
  ASI_low=rep(0, S)
  ASI_high=rep(0, S)
  bias_AOB=rep(0, S)
  diff_D=rep(0, S)
  diff_X=rep(0, S)
  withinvar=rep(0, S)
  acrossvar=rep(0, S)
  totalsd=rep(0, S)
  for(i in 1:S){
    Z=rnorm(n, mean = Mean, sd = SD)
    X_1=rnorm(n)
    X_2=rnorm(n)
    X_3=rnorm(n)
    sigma_epsilon<-matrix(c(1,0.5,0.5,1), nrow = 2, ncol = 2)
    epsilon=rmvnorm(n, mean = c(0,0), sigma = sigma_epsilon)
    
    prob_D=expit(lambda_0_D+lambda_1_D*Z)
    D=rbinom(n, size=1, prob = prob_D)
    Y=beta*D+0.2*X_1+0.5*log(abs(X_2))+0.3*sin(X_3)+epsilon[,1]
    
    data_unmatched<-data.frame(Z, X_1, X_2, X_3, D, Y)
    nomatches <- ceiling(dim(data_unmatched)[1]/2)       # Total number of matches that will be made. If you don't adjust this number you'll have more rows in your output than are necessary.
    
    X.reduced<-data_unmatched[, 2:4]
    dif.travel.time=data_unmatched[,1]
    select1row<-c(1:nrow(data_unmatched))
    matched <- matches(select1row,X.reduced,sinks,dif.travel.time,mindist)
    encouraged.numbers <- matched[,2]
    encouraged.match   <- matched[,1]
    data_encouraged<-data_unmatched[encouraged.numbers,]
    data_control<-data_unmatched[encouraged.match,]
    colnames(data_encouraged)<-c("Z.e", "X_1.e", "X_2.e", "X_3.e", "D.e", "Y.e")
    colnames(data_control)<-c("Z.u", "X_1.u", "X_2.u", "X_3.u", "D.u", "Y.u")
    data_matched<-cbind(data_encouraged, data_control)
    
    tau=tau_max
    lambda_0 = solve_lambda_0(tau, lambda_1)
    mat=0
    if (is.na(lambda_0)){
      mat = -1
    } 
    else {
      result= low_high_endpoint(data_matched, delta, lambda_0, lambda_1)
      ASI_low[i] = result[1]
      ASI_high[i] = result[2]
      if (ASI_low[i] < 0 && ASI_high[i]>0) mat = 1 
      bias_AOB[i]=result[3]
      withinvar[i]=result[4]
      acrossvar[i]=result[5]
      totalsd[i]=result[6]
      diff_D[i]=mean(data_matched$D.e, na.rm = TRUE)-mean(data_matched$D.u, na.rm = TRUE)
      diff_X[i]=mean(0.2*data_matched$X_1.e+0.5*log(abs(data_matched$X_2.e))+0.3*sin(data_matched$X_3.e), na.rm = TRUE)-mean(0.2*data_matched$X_1.u+0.5*log(abs(data_matched$X_2.u))+0.3*sin(data_matched$X_3.u), na.rm = TRUE)
    }
    if (mat==0){
      count_1=count_1+1
    }
    if (mat!=-1){
      count_2=count_2+1
    }
    if(i%%100==0){
      print(i)
    }
  }
  vec_return<-rep(0, 9)
  vec_return[1]=count_1/count_2
  vec_return[2]=sum(bias_AOB)/count_2
  vec_return[3]=sum(diff_D)/count_2
  vec_return[4]=sum(diff_X)/count_2
  vec_return[5]=sum(withinvar)/count_2
  vec_return[6]=sum(acrossvar)/count_2
  vec_return[7]=sum(totalsd)/count_2
  vec_return[8]=sum(ASI_low)/count_2
  vec_return[9]=sum(ASI_high)/count_2
  return(vec_return)
}



n=1000
S=2000
beta_0=0.8
delta_0=0.5
beta=4
delta=10
Mean=0
SD=1
tau_max=0.01


power_IV<-matrix(0, nrow = 100, ncol = 9)

power_IV[1,]<-power_AOB_varytau(sinks=0,mindist=0,delta_0,n,S,beta_0,lambda_1=1,tau_max,lambda_0_D=0,lambda_1_D=1,Mean,SD)
power_IV[2,]<-power_AOB_varytau(sinks=0,mindist=0,delta_0,n,S,beta_0,lambda_1=1.5,tau_max,lambda_0_D=0,lambda_1_D=1,Mean,SD)
power_IV[3,]<-power_AOB_varytau(sinks=0,mindist=0,delta_0,n,S,beta_0,lambda_1=2,tau_max,lambda_0_D=0,lambda_1_D=1,Mean,SD)
power_IV[4,]<-power_AOB_varytau(sinks=0,mindist=0,delta_0,n,S,beta_0,lambda_1=2.5,tau_max,lambda_0_D=0,lambda_1_D=1,Mean,SD)
power_IV[5,]<-power_AOB_varytau(sinks=0,mindist=0,delta_0,n,S,beta_0,lambda_1=3,tau_max,lambda_0_D=0,lambda_1_D=1,Mean,SD)
power_IV[6,]<-power_AOB_varytau(sinks=0,mindist=0,delta_0,n,S,beta_0,lambda_1=3.5,tau_max,lambda_0_D=0,lambda_1_D=1,Mean,SD)
power_IV[7,]<-power_AOB_varytau(sinks=0,mindist=0,delta_0,n,S,beta_0,lambda_1=1,tau_max,lambda_0_D=0,lambda_1_D=1.2,Mean,SD)
power_IV[8,]<-power_AOB_varytau(sinks=0,mindist=0,delta_0,n,S,beta_0,lambda_1=1.5,tau_max,lambda_0_D=0,lambda_1_D=1.2,Mean,SD)
power_IV[9,]<-power_AOB_varytau(sinks=0,mindist=0,delta_0,n,S,beta_0,lambda_1=2,tau_max,lambda_0_D=0,lambda_1_D=1.2,Mean,SD)
power_IV[10,]<-power_AOB_varytau(sinks=0,mindist=0,delta_0,n,S,beta_0,lambda_1=2.5,tau_max,lambda_0_D=0,lambda_1_D=1.2,Mean,SD)
power_IV[11,]<-power_AOB_varytau(sinks=0,mindist=0,delta_0,n,S,beta_0,lambda_1=3,tau_max,lambda_0_D=0,lambda_1_D=1.2,Mean,SD)
power_IV[12,]<-power_AOB_varytau(sinks=0,mindist=0,delta_0,n,S,beta_0,lambda_1=3.5,tau_max,lambda_0_D=0,lambda_1_D=1.2,Mean,SD)
power_IV[13,]<-power_AOB_varytau(sinks=0,mindist=0,delta,n,S,beta,lambda_1=6,tau_max,lambda_0_D=0,lambda_1_D=4,Mean,SD)
power_IV[14,]<-power_AOB_varytau(sinks=0,mindist=0,delta,n,S,beta,lambda_1=9,tau_max,lambda_0_D=0,lambda_1_D=4,Mean,SD)
power_IV[15,]<-power_AOB_varytau(sinks=0,mindist=0,delta,n,S,beta,lambda_1=12,tau_max,lambda_0_D=0,lambda_1_D=4,Mean,SD)
power_IV[16,]<-power_AOB_varytau(sinks=0,mindist=0,delta,n,S,beta,lambda_1=15,tau_max,lambda_0_D=0,lambda_1_D=4,Mean,SD)
power_IV[17,]<-power_AOB_varytau(sinks=0,mindist=0,delta,n,S,beta,lambda_1=18,tau_max,lambda_0_D=0,lambda_1_D=4,Mean,SD)
power_IV[18,]<-power_AOB_varytau(sinks=0,mindist=0,delta,n,S,beta,lambda_1=21,tau_max,lambda_0_D=0,lambda_1_D=4,Mean,SD)
power_IV[19,]<-power_AOB_varytau(sinks=0,mindist=0,delta,n,S,beta,lambda_1=6,tau_max,lambda_0_D=0,lambda_1_D=5,Mean,SD)
power_IV[20,]<-power_AOB_varytau(sinks=0,mindist=0,delta,n,S,beta,lambda_1=9,tau_max,lambda_0_D=0,lambda_1_D=5,Mean,SD)
power_IV[21,]<-power_AOB_varytau(sinks=0,mindist=0,delta,n,S,beta,lambda_1=12,tau_max,lambda_0_D=0,lambda_1_D=5,Mean,SD)
power_IV[22,]<-power_AOB_varytau(sinks=0,mindist=0,delta,n,S,beta,lambda_1=15,tau_max,lambda_0_D=0,lambda_1_D=5,Mean,SD)
power_IV[23,]<-power_AOB_varytau(sinks=0,mindist=0,delta,n,S,beta,lambda_1=18,tau_max,lambda_0_D=0,lambda_1_D=5,Mean,SD)
power_IV[24,]<-power_AOB_varytau(sinks=0,mindist=0,delta,n,S,beta,lambda_1=21,tau_max,lambda_0_D=0,lambda_1_D=5,Mean,SD)




power_IV_st<-matrix(0, nrow = 100, ncol = 9)

power_IV_st[1,]<-power_AOB_varytau(sinks=0.5,mindist=1.3,delta_0,n,S,beta_0,lambda_1=1,tau_max,lambda_0_D=0,lambda_1_D=1,Mean,SD)
power_IV_st[2,]<-power_AOB_varytau(sinks=0.5,mindist=1.3,delta_0,n,S,beta_0,lambda_1=1.5,tau_max,lambda_0_D=0,lambda_1_D=1,Mean,SD)
power_IV_st[3,]<-power_AOB_varytau(sinks=0.5,mindist=1.3,delta_0,n,S,beta_0,lambda_1=2,tau_max,lambda_0_D=0,lambda_1_D=1,Mean,SD)
power_IV_st[4,]<-power_AOB_varytau(sinks=0.5,mindist=1.3,delta_0,n,S,beta_0,lambda_1=2.5,tau_max,lambda_0_D=0,lambda_1_D=1,Mean,SD)
power_IV_st[5,]<-power_AOB_varytau(sinks=0.5,mindist=1.3,delta_0,n,S,beta_0,lambda_1=3,tau_max,lambda_0_D=0,lambda_1_D=1,Mean,SD)
power_IV_st[6,]<-power_AOB_varytau(sinks=0.5,mindist=1.3,delta_0,n,S,beta_0,lambda_1=3.5,tau_max,lambda_0_D=0,lambda_1_D=1,Mean,SD)
power_IV_st[7,]<-power_AOB_varytau(sinks=0.5,mindist=1.3,delta_0,n,S,beta_0,lambda_1=1,tau_max,lambda_0_D=0,lambda_1_D=1.2,Mean,SD)
power_IV_st[8,]<-power_AOB_varytau(sinks=0.5,mindist=1.3,delta_0,n,S,beta_0,lambda_1=1.5,tau_max,lambda_0_D=0,lambda_1_D=1.2,Mean,SD)
power_IV_st[9,]<-power_AOB_varytau(sinks=0.5,mindist=1.3,delta_0,n,S,beta_0,lambda_1=2,tau_max,lambda_0_D=0,lambda_1_D=1.2,Mean,SD)
power_IV_st[10,]<-power_AOB_varytau(sinks=0.5,mindist=1.3,delta_0,n,S,beta_0,lambda_1=2.5,tau_max,lambda_0_D=0,lambda_1_D=1.2,Mean,SD)
power_IV_st[11,]<-power_AOB_varytau(sinks=0.5,mindist=1.3,delta_0,n,S,beta_0,lambda_1=3,tau_max,lambda_0_D=0,lambda_1_D=1.2,Mean,SD)
power_IV_st[12,]<-power_AOB_varytau(sinks=0.5,mindist=1.3,delta_0,n,S,beta_0,lambda_1=3.5,tau_max,lambda_0_D=0,lambda_1_D=1.2,Mean,SD)
power_IV_st[13,]<-power_AOB_varytau(sinks=0.5,mindist=1.3,delta,n,S,beta,lambda_1=6,tau_max,lambda_0_D=0,lambda_1_D=4,Mean,SD)
power_IV_st[14,]<-power_AOB_varytau(sinks=0.5,mindist=1.3,delta,n,S,beta,lambda_1=9,tau_max,lambda_0_D=0,lambda_1_D=4,Mean,SD)
power_IV_st[15,]<-power_AOB_varytau(sinks=0.5,mindist=1.3,delta,n,S,beta,lambda_1=12,tau_max,lambda_0_D=0,lambda_1_D=4,Mean,SD)
power_IV_st[16,]<-power_AOB_varytau(sinks=0.5,mindist=1.3,delta,n,S,beta,lambda_1=15,tau_max,lambda_0_D=0,lambda_1_D=4,Mean,SD)
power_IV_st[17,]<-power_AOB_varytau(sinks=0.5,mindist=1.3,delta,n,S,beta,lambda_1=18,tau_max,lambda_0_D=0,lambda_1_D=4,Mean,SD)
power_IV_st[18,]<-power_AOB_varytau(sinks=0.5,mindist=1.3,delta,n,S,beta,lambda_1=21,tau_max,lambda_0_D=0,lambda_1_D=4,Mean,SD)
power_IV_st[19,]<-power_AOB_varytau(sinks=0.5,mindist=1.3,delta,n,S,beta,lambda_1=6,tau_max,lambda_0_D=0,lambda_1_D=5,Mean,SD)
power_IV_st[20,]<-power_AOB_varytau(sinks=0.5,mindist=1.3,delta,n,S,beta,lambda_1=9,tau_max,lambda_0_D=0,lambda_1_D=5,Mean,SD)
power_IV_st[21,]<-power_AOB_varytau(sinks=0.5,mindist=1.3,delta,n,S,beta,lambda_1=12,tau_max,lambda_0_D=0,lambda_1_D=5,Mean,SD)
power_IV_st[22,]<-power_AOB_varytau(sinks=0.5,mindist=1.3,delta,n,S,beta,lambda_1=15,tau_max,lambda_0_D=0,lambda_1_D=5,Mean,SD)
power_IV_st[23,]<-power_AOB_varytau(sinks=0.5,mindist=1.3,delta,n,S,beta,lambda_1=18,tau_max,lambda_0_D=0,lambda_1_D=5,Mean,SD)
power_IV_st[24,]<-power_AOB_varytau(sinks=0.5,mindist=1.3,delta,n,S,beta,lambda_1=21,tau_max,lambda_0_D=0,lambda_1_D=5,Mean,SD)




