# Programs for strengthening IV project, applying it to data
# Code for matching for Donut Hole IV project

#################################Original IV, sink=0###################

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

#Data storage function
store.matches <- function(matchdataframe,encouraged.numbers, encouraged.match, select1row, dif.travel.time, X, high_level_NICU, medu, death_fatel, losI, matchcount){
  noe=length(encouraged.numbers);
  n.X=ncol(X)
  matchdataframe[(matchcount+1):(matchcount+noe),1]=(dif.travel.time[select1row])[encouraged.numbers];
  matchdataframe[(matchcount+1):(matchcount+noe),2:(n.X+1)]=(X[select1row,])[encouraged.numbers,];
  matchdataframe[(matchcount+1):(matchcount+noe),n.X+2]=(high_level_NICU[select1row])[encouraged.numbers];
  matchdataframe[(matchcount+1):(matchcount+noe),n.X+3]=(losI[select1row])[encouraged.numbers];
  matchdataframe[(matchcount+1):(matchcount+noe),n.X+4]=(medu[select1row])[encouraged.numbers];  
  matchdataframe[(matchcount+1):(matchcount+noe),n.X+5]=(death_fatel[select1row])[encouraged.numbers];
  
  matchdataframe[(matchcount+1):(matchcount+noe),n.X+6]=(dif.travel.time[select1row])[encouraged.match];
  matchdataframe[(matchcount+1):(matchcount+noe),(n.X+7):(2*n.X+6)]=(X[select1row,])[encouraged.match,];
  matchdataframe[(matchcount+1):(matchcount+noe),2*n.X+7]=(high_level_NICU[select1row])[encouraged.match];
  matchdataframe[(matchcount+1):(matchcount+noe),2*n.X+8]=(losI[select1row])[encouraged.match];
  matchdataframe[(matchcount+1):(matchcount+noe),2*n.X+9]=(medu[select1row])[encouraged.match];
  matchdataframe[(matchcount+1):(matchcount+noe),2*n.X+10]=(death_fatel[select1row])[encouraged.match];
  return(matchdataframe)
}


data.names     <- paste("C:/Users/siyuheng/Desktop/Real Data Analysis for Strengthening IV/Data/Preemie Data/preemies",sep="",1995:2005,".csv")
matches.names  <- paste("C:/Users/siyuheng/Desktop/Real Data Analysis for Strengthening IV/Data/Original IV no medu/matched_originalIV_no_medu",sep="",1995:2005,".csv")


for(k in 1:length(data.names)){
  ########
  ##data##
  ########
  preemiedata=read.csv(data.names[k])
  dim(preemiedata)
  
  
  #####################
  ##Adjustable Inputs##
  #####################
  sinks <- 0.00                             # Fraction of data to be matched to sinks
  nomatches <- ceiling(dim(preemiedata)[1]/2)       # Total number of matches that will be made. If you don't adjust this number you'll have more rows in your output than are necessary.
  mindist <- 0                            # Minimum distance babies are allowed to be apart and not get a penalty
  
  #############################
  ##Prepare data and matrices##
  #############################
  
  # Create missing data indicators
  # Impute gestational age and birthweight using a regression of one
  # variable on the other
  bthwght=preemiedata$bthwght;
  bthwghtreg=lm(bthwght~preemiedata$gestage_weeks);
  bthwght[is.na(bthwght)]=coef(bthwghtreg)[1]+coef(bthwghtreg)[2]*preemiedata$gestage_weeks[is.na(bthwght)];
  gestage=preemiedata$gestage_weeks;
  gestagereg=lm(gestage~preemiedata$bthwght);
  gestage[is.na(gestage)]=coef(gestagereg)[1]+coef(gestagereg)[2]*preemiedata$bthwght[is.na(gestage)]; #There's no record with both missing gestational age and birth weight
  gestdiabetes=preemiedata$Gestational_DiabetesM  #gestational diabetes is complete
  singlebirth=as.numeric(preemiedata$multiples==1); #Single birth is complete
  parity=preemiedata$parity;
  parity[is.na(parity)]=mean(parity,na.rm=TRUE);
  
  #Mother's covariates
  ageM=preemiedata$ageM; #Mother's age is complete
  medu=preemiedata$medu;
  medu[is.na(medu)]=mean(medu, na.rm = TRUE)
  race=preemiedata$raceM;
  white <- as.numeric(race==1)
  white[is.na(white)]=0
  race.mis <- as.numeric(is.na(race))
  
  #Mother's neighborhood
  below_poverty=preemiedata$below_poverty
  below_poverty[is.na(below_poverty)]=mean(below_poverty, na.rm = TRUE)
  
  X.reduced=cbind(bthwght,gestage,gestdiabetes,singlebirth,parity,ageM,white,race.mis,below_poverty);
  
  
  
  #create high level NICU indicator
  high_level_NICU <- preemiedata$vol_level_big_2500
  # Set up strata on race of mother, age of mother and gestational age and match within strata
  dif.travel.time=preemiedata$diff_travel_2500;
  losI = preemiedata$losI
  death_fatel = preemiedata$death_fetal
  
  
  ########################
  ##Initialize Variables##
  ########################
  
  X.names <- colnames(X.reduced)
  
  
  varnames=c("dif.travel.time.e",paste(X.names,sep="",".e"),"high_level_NICU.e", "losI.e", "medu.e", "death_fatel.e",
             "dif.travel.time.u",paste(X.names,sep="",".u"),"high_level_NICU.u", "losI.u", "medu.u", "death_fatel.u")
  
  matchdataframe=data.frame(matrix(rep(NA,nomatches*length(varnames)),ncol=length(varnames)));
  names(matchdataframe)=varnames;
  matchcount <- 0
  
  ##############################
  ##Run the matching by strata##    New strata - adaptive
  ##############################
  
  # gest. age <=28
  select1=as.logical(gestage<=28*(is.na(dif.travel.time)==FALSE));
  select1row.age.strata=which(select1)
  length(select1row.age.strata)
  num.strata <- max(round(length(select1row.age.strata)/1000),2)
  num.strata
  probs <- (0:num.strata)/num.strata
  quans <- quantile(bthwght[select1row.age.strata],probs=probs)
  for(i in 1:(num.strata-1)){
    select1=as.logical((gestage<=28)*(is.na(dif.travel.time)==FALSE)*(bthwght>=quans[i])*(bthwght<quans[i+1]));
    select1row=which(select1)
    #create matches
    matched <- matches(select1row,X.reduced,sinks,dif.travel.time,mindist)
    encouraged.numbers <- matched[,1]
    encouraged.match   <- matched[,2]
    # Save matches
    matchdataframe <- store.matches(matchdataframe,encouraged.numbers, encouraged.match, select1row, dif.travel.time, X.reduced, high_level_NICU, medu, death_fatel, losI, matchcount)
    matchcount <- length(encouraged.numbers) + matchcount
    print(matchcount)
  }
  #do for the last birth weight stratum
  i <- num.strata
  select1=as.logical((gestage<=28)*(is.na(dif.travel.time)==FALSE)*(bthwght>=quans[i])*(bthwght<quans[i+1]));
  select1row=which(select1)
  #create matches
  matched <- matches(select1row,X.reduced,sinks,dif.travel.time,mindist)
  encouraged.numbers <- matched[,1]
  encouraged.match   <- matched[,2]
  # Save matches
  matchdataframe <- store.matches(matchdataframe,encouraged.numbers, encouraged.match, select1row, dif.travel.time, X.reduced, high_level_NICU, medu, death_fatel, losI, matchcount)
  matchcount <- length(encouraged.numbers) + matchcount
  matchcount
  
  
  # Gestational Age for 28<gestage<=32
  select1=as.logical((gestage>28)*(gestage<=32)*(is.na(dif.travel.time)==FALSE));
  select1row.age.strata=which(select1)
  length(select1row.age.strata)
  num.strata <- max(round(length(select1row.age.strata)/1000),2)
  num.strata
  probs <- (0:num.strata)/num.strata
  quans <- quantile(bthwght[select1row.age.strata],probs=probs)
  for(i in 1:(num.strata-1)){
    select1=as.logical((gestage>28)*(gestage<=32)*(is.na(dif.travel.time)==FALSE)*(bthwght>=quans[i])*(bthwght<quans[i+1]));
    select1row=which(select1)
    #create matches
    matched <- matches(select1row,X.reduced,sinks,dif.travel.time,mindist)
    encouraged.numbers <- matched[,1]
    encouraged.match   <- matched[,2]
    # Save matches
    matchdataframe <- store.matches(matchdataframe,encouraged.numbers, encouraged.match, select1row, dif.travel.time, X.reduced, high_level_NICU, medu, death_fatel, losI, matchcount)
    matchcount <- length(encouraged.numbers) + matchcount
    print(matchcount)
  }
  #do for the last birth weight stratum
  i <- num.strata
  select1=as.logical((gestage>28)*(gestage<=32)*(is.na(dif.travel.time)==FALSE)*(bthwght>=quans[i])*(bthwght<=quans[i+1]));
  select1row=which(select1)
  #create matches
  matched <- matches(select1row,X.reduced,sinks,dif.travel.time,mindist)
  encouraged.numbers <- matched[,1]
  encouraged.match   <- matched[,2]
  # Save matches
  matchdataframe <- store.matches(matchdataframe,encouraged.numbers, encouraged.match, select1row, dif.travel.time, X.reduced, high_level_NICU, medu, death_fatel, losI, matchcount)
  matchcount <- length(encouraged.numbers) + matchcount
  matchcount
  
  
  
  
  
  
  
  
  
  # Gestational Age for 32<gestage<37
  select1=as.logical((gestage>32)*(gestage<37)*(is.na(dif.travel.time)==FALSE));
  select1row.age.strata=which(select1)
  length(select1row.age.strata)
  num.strata <- max(round(length(select1row.age.strata)/1000),2)
  num.strata
  probs <- (0:num.strata)/num.strata
  quans <- quantile(bthwght[select1row.age.strata],probs=probs)
  for(i in 1:(num.strata-1)){
    select1=as.logical((gestage>32)*(gestage<37)*(is.na(dif.travel.time)==FALSE)*(bthwght>=quans[i])*(bthwght<quans[i+1]));
    select1row=which(select1)
    #create matches
    matched <- matches(select1row,X.reduced,sinks,dif.travel.time,mindist)
    encouraged.numbers <- matched[,1]
    encouraged.match   <- matched[,2]
    # Save matches
    matchdataframe <- store.matches(matchdataframe,encouraged.numbers, encouraged.match, select1row, dif.travel.time, X.reduced, high_level_NICU, medu, death_fatel, losI, matchcount)
    matchcount <- length(encouraged.numbers) + matchcount
    print(matchcount)
  }
  i <- num.strata
  select1=as.logical((gestage>32)*(gestage<37)*(is.na(dif.travel.time)==FALSE)*(bthwght>=quans[i])*(bthwght<=quans[i+1]));
  select1row=which(select1)
  #create matches
  matched <- matches(select1row,X.reduced,sinks,dif.travel.time,mindist)
  encouraged.numbers <- matched[,1]
  encouraged.match   <- matched[,2]
  # Save matches
  matchdataframe <- store.matches(matchdataframe,encouraged.numbers, encouraged.match, select1row, dif.travel.time, X.reduced, high_level_NICU, medu, death_fatel, losI, matchcount)
  matchcount <- length(encouraged.numbers) + matchcount
  matchcount
  
  
  
  # Gestational Age for gestage==37
  select1=as.logical((gestage==37)*(is.na(dif.travel.time)==FALSE));
  select1row.age.strata=which(select1)
  length(select1row.age.strata)
  num.strata <- max(round(length(select1row.age.strata)/1000),2)
  num.strata
  probs <- (0:num.strata)/num.strata
  quans <- quantile(bthwght[select1row.age.strata],probs=probs)
  for(i in 1:(num.strata-1)){
    select1=as.logical((gestage==37)*(is.na(dif.travel.time)==FALSE)*(bthwght>=quans[i])*(bthwght<quans[i+1]));
    select1row=which(select1)
    #create matches
    matched <- matches(select1row,X.reduced,sinks,dif.travel.time,mindist)
    encouraged.numbers <- matched[,1]
    encouraged.match   <- matched[,2]
    # Save matches
    matchdataframe <- store.matches(matchdataframe,encouraged.numbers, encouraged.match, select1row, dif.travel.time, X.reduced, high_level_NICU, medu, death_fatel, losI, matchcount)
    matchcount <- length(encouraged.numbers) + matchcount
    print(matchcount)
  }
  #do for the last birth weight stratum
  i <- num.strata
  select1=as.logical((gestage==37)*(is.na(dif.travel.time)==FALSE)*(bthwght>=quans[i])*(bthwght<=quans[i+1]));
  select1row=which(select1)
  #create matches
  matched <- matches(select1row,X.reduced,sinks,dif.travel.time,mindist)
  encouraged.numbers <- matched[,1]
  encouraged.match   <- matched[,2]
  # Save matches
  matchdataframe <- store.matches(matchdataframe,encouraged.numbers, encouraged.match, select1row, dif.travel.time, X.reduced, high_level_NICU, medu, death_fatel, losI, matchcount)
  matchcount <- length(encouraged.numbers) + matchcount
  matchcount
  
  
  # Save data
  write.csv(matchdataframe[1:nomatches,],file=matches.names[k]);
  print(k)
}







#####################################Strengthened IV sink=0.5###############################

data.names     <- paste("C:/Users/siyuheng/Desktop/Real Data Analysis for Strengthening IV/Data/Preemie Data/preemies",sep="",1995:2005,".csv")
matches.names  <- paste("C:/Users/siyuheng/Desktop/Real Data Analysis for Strengthening IV/Data/Stronger IV no medu/matched_strongerIV_no_medu",sep="",1995:2005,".csv")

for(k in 1:length(data.names)){
  ########
  ##data##
  ########
  preemiedata=read.csv(data.names[k])
  dim(preemiedata)
  
  
  #####################
  ##Adjustable Inputs##
  #####################
  sinks <- 0.5                             # Fraction of data to be matched to sinks
  nomatches <- ceiling(dim(preemiedata)[1]/2)       # Total number of matches that will be made. If you don't adjust this number you'll have more rows in your output than are necessary.
  mindist <- 25                            # Minimum distance babies are allowed to be apart and not get a penalty
  
  #############################
  ##Prepare data and matrices##
  #############################
  
  # Create missing data indicators
  # Impute gestational age and birthweight using a regression of one
  # variable on the other
  bthwght=preemiedata$bthwght;
  bthwghtreg=lm(bthwght~preemiedata$gestage_weeks);
  bthwght[is.na(bthwght)]=coef(bthwghtreg)[1]+coef(bthwghtreg)[2]*preemiedata$gestage_weeks[is.na(bthwght)];
  gestage=preemiedata$gestage_weeks;
  gestagereg=lm(gestage~preemiedata$bthwght);
  gestage[is.na(gestage)]=coef(gestagereg)[1]+coef(gestagereg)[2]*preemiedata$bthwght[is.na(gestage)]; #There's no record with both missing gestational age and birth weight
  gestdiabetes=preemiedata$Gestational_DiabetesM  #gestational diabetes is complete
  singlebirth=as.numeric(preemiedata$multiples==1); #Single birth is complete
  parity=preemiedata$parity;
  parity[is.na(parity)]=mean(parity,na.rm=TRUE);
  
  #Mother's covariates
  ageM=preemiedata$ageM; #Mother's age is complete
  medu=preemiedata$medu;
  medu[is.na(medu)]=mean(medu, na.rm = TRUE)
  race=preemiedata$raceM;
  white <- as.numeric(race==1)
  white[is.na(white)]=0
  race.mis <- as.numeric(is.na(race))
  
  #Mother's neighborhood
  below_poverty=preemiedata$below_poverty
  below_poverty[is.na(below_poverty)]=mean(below_poverty, na.rm = TRUE)
  
  X.reduced=cbind(bthwght,gestage,gestdiabetes,singlebirth,parity,ageM,white,race.mis,below_poverty);
  
  
  
  #create high level NICU indicator
  high_level_NICU <- preemiedata$vol_level_big_2500
  # Set up strata on race of mother, age of mother and gestational age and match within strata
  dif.travel.time=preemiedata$diff_travel_2500;
  losI = preemiedata$losI
  death_fatel = preemiedata$death_fetal
  
  
  ########################
  ##Initialize Variables##
  ########################
  
  X.names <- colnames(X.reduced)
  
  
  varnames=c("dif.travel.time.e",paste(X.names,sep="",".e"),"high_level_NICU.e", "losI.e", "medu.e", "death_fatel.e",
             "dif.travel.time.u",paste(X.names,sep="",".u"),"high_level_NICU.u", "losI.u", "medu.u", "death_fatel.u")
  
  matchdataframe=data.frame(matrix(rep(NA,nomatches*length(varnames)),ncol=length(varnames)));
  names(matchdataframe)=varnames;
  matchcount <- 0
  
  ##############################
  ##Run the matching by strata##    New strata - adaptive
  ##############################
  
  # gest. age <=28
  select1=as.logical(gestage<=28*(is.na(dif.travel.time)==FALSE));
  select1row.age.strata=which(select1)
  length(select1row.age.strata)
  num.strata <- max(round(length(select1row.age.strata)/1000),2)
  num.strata
  probs <- (0:num.strata)/num.strata
  quans <- quantile(bthwght[select1row.age.strata],probs=probs)
  for(i in 1:(num.strata-1)){
    select1=as.logical((gestage<=28)*(is.na(dif.travel.time)==FALSE)*(bthwght>=quans[i])*(bthwght<quans[i+1]));
    select1row=which(select1)
    #create matches
    matched <- matches(select1row,X.reduced,sinks,dif.travel.time,mindist)
    encouraged.numbers <- matched[,1]
    encouraged.match   <- matched[,2]
    # Save matches
    matchdataframe <- store.matches(matchdataframe,encouraged.numbers, encouraged.match, select1row, dif.travel.time, X.reduced, high_level_NICU, medu, death_fatel, losI, matchcount)
    matchcount <- length(encouraged.numbers) + matchcount
    print(matchcount)
  }
  #do for the last birth weight stratum
  i <- num.strata
  select1=as.logical((gestage<=28)*(is.na(dif.travel.time)==FALSE)*(bthwght>=quans[i])*(bthwght<quans[i+1]));
  select1row=which(select1)
  #create matches
  matched <- matches(select1row,X.reduced,sinks,dif.travel.time,mindist)
  encouraged.numbers <- matched[,1]
  encouraged.match   <- matched[,2]
  # Save matches
  matchdataframe <- store.matches(matchdataframe,encouraged.numbers, encouraged.match, select1row, dif.travel.time, X.reduced, high_level_NICU, medu, death_fatel, losI, matchcount)
  matchcount <- length(encouraged.numbers) + matchcount
  matchcount
  
  
  # Gestational Age for 28<gestage<=32
  select1=as.logical((gestage>28)*(gestage<=32)*(is.na(dif.travel.time)==FALSE));
  select1row.age.strata=which(select1)
  length(select1row.age.strata)
  num.strata <- max(round(length(select1row.age.strata)/1000),2)
  num.strata
  probs <- (0:num.strata)/num.strata
  quans <- quantile(bthwght[select1row.age.strata],probs=probs)
  for(i in 1:(num.strata-1)){
    select1=as.logical((gestage>28)*(gestage<=32)*(is.na(dif.travel.time)==FALSE)*(bthwght>=quans[i])*(bthwght<quans[i+1]));
    select1row=which(select1)
    #create matches
    matched <- matches(select1row,X.reduced,sinks,dif.travel.time,mindist)
    encouraged.numbers <- matched[,1]
    encouraged.match   <- matched[,2]
    # Save matches
    matchdataframe <- store.matches(matchdataframe,encouraged.numbers, encouraged.match, select1row, dif.travel.time, X.reduced, high_level_NICU, medu, death_fatel, losI, matchcount)
    matchcount <- length(encouraged.numbers) + matchcount
    print(matchcount)
  }
  #do for the last birth weight stratum
  i <- num.strata
  select1=as.logical((gestage>28)*(gestage<=32)*(is.na(dif.travel.time)==FALSE)*(bthwght>=quans[i])*(bthwght<=quans[i+1]));
  select1row=which(select1)
  #create matches
  matched <- matches(select1row,X.reduced,sinks,dif.travel.time,mindist)
  encouraged.numbers <- matched[,1]
  encouraged.match   <- matched[,2]
  # Save matches
  matchdataframe <- store.matches(matchdataframe,encouraged.numbers, encouraged.match, select1row, dif.travel.time, X.reduced, high_level_NICU, medu, death_fatel, losI, matchcount)
  matchcount <- length(encouraged.numbers) + matchcount
  matchcount
  
  
  
  # Gestational Age for 32<gestage<37
  select1=as.logical((gestage>32)*(gestage<37)*(is.na(dif.travel.time)==FALSE));
  select1row.age.strata=which(select1)
  length(select1row.age.strata)
  num.strata <- max(round(length(select1row.age.strata)/1000),2)
  num.strata
  probs <- (0:num.strata)/num.strata
  quans <- quantile(bthwght[select1row.age.strata],probs=probs)
  for(i in 1:(num.strata-1)){
    select1=as.logical((gestage>32)*(gestage<37)*(is.na(dif.travel.time)==FALSE)*(bthwght>=quans[i])*(bthwght<quans[i+1]));
    select1row=which(select1)
    #create matches
    matched <- matches(select1row,X.reduced,sinks,dif.travel.time,mindist)
    encouraged.numbers <- matched[,1]
    encouraged.match   <- matched[,2]
    # Save matches
    matchdataframe <- store.matches(matchdataframe,encouraged.numbers, encouraged.match, select1row, dif.travel.time, X.reduced, high_level_NICU, medu, death_fatel, losI, matchcount)
    matchcount <- length(encouraged.numbers) + matchcount
    print(matchcount)
  }
  i <- num.strata
  select1=as.logical((gestage>32)*(gestage<37)*(is.na(dif.travel.time)==FALSE)*(bthwght>=quans[i])*(bthwght<=quans[i+1]));
  select1row=which(select1)
  #create matches
  matched <- matches(select1row,X.reduced,sinks,dif.travel.time,mindist)
  encouraged.numbers <- matched[,1]
  encouraged.match   <- matched[,2]
  # Save matches
  matchdataframe <- store.matches(matchdataframe,encouraged.numbers, encouraged.match, select1row, dif.travel.time, X.reduced, high_level_NICU, medu, death_fatel, losI, matchcount)
  matchcount <- length(encouraged.numbers) + matchcount
  matchcount
  
  
  
  # Gestational Age for gestage==37
  select1=as.logical((gestage==37)*(is.na(dif.travel.time)==FALSE));
  select1row.age.strata=which(select1)
  length(select1row.age.strata)
  num.strata <- max(round(length(select1row.age.strata)/1000),2)
  num.strata
  probs <- (0:num.strata)/num.strata
  quans <- quantile(bthwght[select1row.age.strata],probs=probs)
  for(i in 1:(num.strata-1)){
    select1=as.logical((gestage==37)*(is.na(dif.travel.time)==FALSE)*(bthwght>=quans[i])*(bthwght<quans[i+1]));
    select1row=which(select1)
    #create matches
    matched <- matches(select1row,X.reduced,sinks,dif.travel.time,mindist)
    encouraged.numbers <- matched[,1]
    encouraged.match   <- matched[,2]
    # Save matches
    matchdataframe <- store.matches(matchdataframe,encouraged.numbers, encouraged.match, select1row, dif.travel.time, X.reduced, high_level_NICU, medu, death_fatel, losI, matchcount)
    matchcount <- length(encouraged.numbers) + matchcount
    print(matchcount)
  }
  #do for the last birth weight stratum
  i <- num.strata
  select1=as.logical((gestage==37)*(is.na(dif.travel.time)==FALSE)*(bthwght>=quans[i])*(bthwght<=quans[i+1]));
  select1row=which(select1)
  #create matches
  matched <- matches(select1row,X.reduced,sinks,dif.travel.time,mindist)
  encouraged.numbers <- matched[,1]
  encouraged.match   <- matched[,2]
  # Save matches
  matchdataframe <- store.matches(matchdataframe,encouraged.numbers, encouraged.match, select1row, dif.travel.time, X.reduced, high_level_NICU, medu, death_fatel, losI, matchcount)
  matchcount <- length(encouraged.numbers) + matchcount
  matchcount
  
  
  # Save data
  write.csv(matchdataframe[1:nomatches,],file=matches.names[k]);
  print(k)
}



###########################Combine all the matched data#####################
matches.names_originalIV  <- paste("C:/Users/siyuheng/Desktop/Real Data Analysis for Strengthening IV/Data/Original IV no medu/matched_originalIV_no_medu",sep="",1995:2005,".csv")
matches.names_strongerIV  <- paste("C:/Users/siyuheng/Desktop/Real Data Analysis for Strengthening IV/Data/Stronger IV no medu/matched_strongerIV_no_medu",sep="",1995:2005,".csv")

matched_allyears_originalIV<-NULL
matched_allyears_strongerIV<-NULL

for (i in 1:length(matches.names_originalIV)){
  tmp=read.csv(matches.names_originalIV[i])
  matched_allyears_originalIV=rbind(matched_allyears_originalIV, tmp)
  print(i)
}

write.csv(matched_allyears_originalIV, "C:/Users/siyuheng/Desktop/Real Data Analysis for Strengthening IV/Data/Original IV no medu/allyears_matched_originalIV_no_medu.csv")

for (i in 1:length(matches.names_strongerIV)){
  tmp=read.csv(matches.names_strongerIV[i])
  matched_allyears_strongerIV=rbind(matched_allyears_strongerIV, tmp)
  print(i)
}

write.csv(matched_allyears_strongerIV, "C:/Users/siyuheng/Desktop/Real Data Analysis for Strengthening IV/Data/Stronger IV no medu/allyears_matched_strongerIV_no_medu.csv")

