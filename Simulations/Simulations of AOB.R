library(mvtnorm)
library(nbpMatching)

#############
##Functions##
#############
#Andreas Buja's Principal Component Analysis code
# site <- "http://www-stat.wharton.upenn.edu/~buja/STAT-541/"
# source(paste(site, "collinearity-pca.R", sep=""))


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
  cols.w.var <-  which(apply(X[select1row, ],2,var)>0)        #Only give the columns with variation, otherwise "smahal" freaks out
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


bias_AOB<-function(sinks, mindist, n, S, beta, p, beta_1, beta_2, C, Mean, Sd){
  bias<-rep(0, 4)
  bias_D<-rep(0, S)
  bias_U<-rep(0, S)
  bias_X<-rep(0, S)
  bias_total<-rep(0, S)
  for(i in 1:S){
    Z=rnorm(n, mean = Mean, sd = Sd)
    X_1=rnorm(n)
    X_2=rnorm(n)
    e_3=rnorm(n)
    sigma_epsilon<-matrix(c(1,0.5,0.5,1), nrow = 2, ncol = 2)
    epsilon=rmvnorm(n, mean = c(0,0), sigma = sigma_epsilon)
    
    U=Z^p+e_3
    D=as.numeric(beta_1*Z^3+beta_2*Z+0.2*X_1+0.4*X_2+epsilon[,2]>C)
    R=beta*D+sin(X_1)+X_2^3+U+epsilon[,1]
    
    data_unmatched<-data.frame(Z, X_1, X_2, D, U, R)
    nomatches <- ceiling(dim(data_unmatched)[1]/2)       # Total number of matches that will be made. If you don't adjust this number you'll have more rows in your output than are necessary.
    
    X.reduced<-data_unmatched[, c(2,3)]
    dif.travel.time=data_unmatched$Z
    select1row<-c(1:nrow(data_unmatched))
    matched <- matches(select1row,X.reduced,sinks,dif.travel.time,mindist)
    encouraged.numbers <- matched[,1]
    encouraged.match   <- matched[,2]
    data_encouraged<-data_unmatched[encouraged.numbers,]
    data_control<-data_unmatched[encouraged.match,]
    
    bias_D[i]=mean(data_encouraged$D)-mean(data_control$D)
    bias_U[i]=mean(data_encouraged$U)-mean(data_control$U)
    bias_total[i]=(mean(data_encouraged$R)-mean(data_control$R))/bias_D[i]-beta
    if (i%%200==0){
      print(i)
    }
  }
  bias[1]=mean(bias_D)
  bias[2]=mean(bias_U)
  bias[3]=mean(bias_U/bias_D)
  bias[4]=mean(bias_total)
  return(bias)
}

bias_IV<-matrix(0, nrow = 100, ncol = 4)

n=400
S=20000

 
bias_IV[1,]=bias_AOB(sinks=0, mindist=0, n, S, beta = 5, p = 1, beta_1 = 0, beta_2 = 1, C = 0, Mean = 0, Sd = 1)
bias_IV[2,]=bias_AOB(sinks=0.5, mindist=8, n, S, beta = 5, p = 1, beta_1 = 0, beta_2 = 1, C = 0, Mean = 0, Sd = 1)



bias_AOB_2<-function(sinks, mindist, n, S, beta, beta_1, beta_2, C, Mean, Sd){
  bias<-rep(0, 4)
  bias_D<-rep(0, S)
  bias_U<-rep(0, S)
  bias_X<-rep(0, S)
  bias_total<-rep(0, S)
  for(i in 1:S){
    Z=rnorm(n, mean = Mean, sd = Sd)
    X_1=rnorm(n)
    X_2=rnorm(n)
    e_3=rnorm(n)
    sigma_epsilon<-matrix(c(1,0.5,0.5,1), nrow = 2, ncol = 2)
    epsilon=rmvnorm(n, mean = c(0,0), sigma = sigma_epsilon)
    
    U=1/(Z-Mean)+e_3
    D=as.numeric(beta_1*Z^3+beta_2*Z+0.2*X_1+0.4*X_2+epsilon[,2]>C)
    R=beta*D+sin(X_1)+X_2^3+U+epsilon[,1]
    
    data_unmatched<-data.frame(Z, X_1, X_2, D, U, R)
    nomatches <- ceiling(dim(data_unmatched)[1]/2)       # Total number of matches that will be made. If you don't adjust this number you'll have more rows in your output than are necessary.
    
    X.reduced<-data_unmatched[, c(2,3)]
    dif.travel.time=data_unmatched$Z
    select1row<-c(1:nrow(data_unmatched))
    matched <- matches(select1row,X.reduced,sinks,dif.travel.time,mindist)
    encouraged.numbers <- matched[,1]
    encouraged.match   <- matched[,2]
    data_encouraged<-data_unmatched[encouraged.numbers,]
    data_control<-data_unmatched[encouraged.match,]
    
    bias_D[i]=mean(data_encouraged$D)-mean(data_control$D)
    bias_U[i]=mean(data_encouraged$U)-mean(data_control$U)
    bias_total[i]=(mean(data_encouraged$R)-mean(data_control$R))/bias_D[i]-beta
    if (i%%200==0){
      print(i)
    }
  }
  bias[1]=mean(bias_D)
  bias[2]=mean(bias_U)
  bias[3]=mean(bias_U/bias_D)
  bias[4]=mean(bias_total)
  return(bias)
}


bias_IV[3,]=bias_AOB_2(sinks=0, mindist=0, n, S, beta = 5, beta_1 = 1, beta_2 = 1, C = 4, Mean = 1, Sd = 5)
bias_IV[4,]=bias_AOB_2(sinks=0.5, mindist=8, n, S, beta = 5, beta_1 = 1, beta_2 = 1, C = 4, Mean = 1, Sd = 5)


bias_AOB_3<-function(sinks, mindist, n, S, beta, beta_1, beta_2, C, Mean, Sd){
  bias<-rep(0, 4)
  bias_D<-rep(0, S)
  bias_U<-rep(0, S)
  bias_X<-rep(0, S)
  bias_total<-rep(0, S)
  for(i in 1:S){
    Z=rnorm(n, mean = Mean, sd = Sd)
    X_1=rnorm(n)
    X_2=rnorm(n)
    e_3=rnorm(n)
    sigma_epsilon<-matrix(c(1,0.5,0.5,1), nrow = 2, ncol = 2)
    epsilon=rmvnorm(n, mean = c(0,0), sigma = sigma_epsilon)
    
    U=exp(Z)+e_3
    D=as.numeric(beta_1*Z^3+beta_2*Z+0.2*X_1+0.4*X_2+epsilon[,2]>C)
    R=beta*D+sin(X_1)+X_2^3+U+epsilon[,1]
    
    data_unmatched<-data.frame(Z, X_1, X_2, D, U, R)
    nomatches <- ceiling(dim(data_unmatched)[1]/2)       # Total number of matches that will be made. If you don't adjust this number you'll have more rows in your output than are necessary.
    
    X.reduced<-data_unmatched[, c(2,3)]
    dif.travel.time=data_unmatched$Z
    select1row<-c(1:nrow(data_unmatched))
    matched <- matches(select1row,X.reduced,sinks,dif.travel.time,mindist)
    encouraged.numbers <- matched[,1]
    encouraged.match   <- matched[,2]
    data_encouraged<-data_unmatched[encouraged.numbers,]
    data_control<-data_unmatched[encouraged.match,]
    
    bias_D[i]=mean(data_encouraged$D)-mean(data_control$D)
    bias_U[i]=mean(data_encouraged$U)-mean(data_control$U)
    bias_total[i]=(mean(data_encouraged$R)-mean(data_control$R))/bias_D[i]-beta
    if (i%%200==0){
      print(i)
    }
  }
  bias[1]=mean(bias_D)
  bias[2]=mean(bias_U)
  bias[3]=mean(bias_U/bias_D)
  bias[4]=mean(bias_total)
  return(bias)
}


bias_IV[5,]=bias_AOB_3(sinks=0, mindist=0, n, S, beta = 5, beta_1 = 1, beta_2 = 1, C = 4, Mean = 1, Sd = 5)
bias_IV[6,]=bias_AOB_3(sinks=0.5, mindist=8, n, S, beta = 5, beta_1 = 1, beta_2 = 1, C = 4, Mean = 1, Sd = 5)




