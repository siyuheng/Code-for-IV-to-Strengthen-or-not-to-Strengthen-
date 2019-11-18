# Programs for strengthening IV project, applying it to data
# Code for matching for Donut Hole IV project

#####################################Strengthened IV via debiased matching###############################

library(gurobi)


##############Make the balance table for 2005 years###################
originalIV<-read.csv("C:/Users/siyuheng/Desktop/Real Data Analysis for Strengthening IV/Data/Debiased IV/incomplete_matched_originalIV2005.csv")
near_mean=rep(0,7)
far_mean=rep(0,7)
for (i in 1:length(near_mean)){
  near_mean[i]=mean(originalIV[,1+i], na.rm = TRUE) #Check the number of the column every time!!!
  far_mean[i]=mean(originalIV[,12+i], na.rm = TRUE)
}

cov_dif = abs(near_mean[2:7]-far_mean[2:7]) ###Check the order of the covariates!
IV_diff = 2*abs(near_mean[1]-far_mean[1])
near_far=0 ####Indicator for near or far receives the encouragement


#Data storage function
store.matches <- function(matchdataframe,encouraged.numbers, encouraged.match, select1row, dif.travel.time, X, high_level_NICU, death, death_fatel, losI, matchcount){
  noe=length(encouraged.numbers);
  n.X=ncol(X)
  matchdataframe[(matchcount+1):(matchcount+noe),1]=(dif.travel.time[select1row])[encouraged.numbers];
  matchdataframe[(matchcount+1):(matchcount+noe),2:(n.X+1)]=(X[select1row,])[encouraged.numbers,];
  matchdataframe[(matchcount+1):(matchcount+noe),n.X+2]=(high_level_NICU[select1row])[encouraged.numbers];
  matchdataframe[(matchcount+1):(matchcount+noe),n.X+3]=(losI[select1row])[encouraged.numbers];
  matchdataframe[(matchcount+1):(matchcount+noe),n.X+4]=(death[select1row])[encouraged.numbers];  
  matchdataframe[(matchcount+1):(matchcount+noe),n.X+5]=(death_fatel[select1row])[encouraged.numbers];
  
  matchdataframe[(matchcount+1):(matchcount+noe),n.X+6]=(dif.travel.time[select1row])[encouraged.match];
  matchdataframe[(matchcount+1):(matchcount+noe),(n.X+7):(2*n.X+6)]=(X[select1row,])[encouraged.match,];
  matchdataframe[(matchcount+1):(matchcount+noe),2*n.X+7]=(high_level_NICU[select1row])[encouraged.match];
  matchdataframe[(matchcount+1):(matchcount+noe),2*n.X+8]=(losI[select1row])[encouraged.match];
  matchdataframe[(matchcount+1):(matchcount+noe),2*n.X+9]=(death[select1row])[encouraged.match];
  matchdataframe[(matchcount+1):(matchcount+noe),2*n.X+10]=(death_fatel[select1row])[encouraged.match];
  return(matchdataframe)
}


model <- list()
IV_debias<-function(X, IV, delta, phi){
  L=nrow(X)
  B=matrix(0, nrow = (L+1 + 2*length(delta)), ncol = (L*(L-1))/2)
  for (l in 1:L){
    count_1=0
    if (l>1){
      for (m in 1:(l-1)){
        B[l, l-m+count_1]=1
        count_1=count_1+L-m
      }
    }
    if (l<L){
      for (m in (l+1):L){
        B[l, m-l+count_1]=1
      }
    }
  }
  count_2=0
  for (l in 1:(L-1)){
    for (m in (l+1):L){
      B[L+1, count_2+m-l]=IV[l]-IV[m]-phi
    }
    count_2=count_2+L-l
  }
  for (i in 1:length(delta)){
    count_3=0
    for (l in 1:(L-1)){
      for (m in (l+1):L){
        B[L+2*i, count_3+m-l]=X[l, i]-X[m, i]-delta[i]
        B[L+2*i+1, count_3+m-l]=X[m, i]-X[l, i]-delta[i]
      }
      count_3=count_3+L-l
    }
  }
  model$A          <- B
  model$obj        <- rep(1, (L*(L-1))/2)
  model$modelsense <- 'max'
  model$rhs        <- c( rep(1, L), rep(0, 2*length(delta)+1))
  model$sense      <- c(rep('<=', L), c('>='), rep('<=', 2*length(delta)))
  model$vtype      <- 'B'
  params <- list(OutputFlag=0)
  result <- gurobi(model, params)
  return(result$x)
}

translate_index<-function(a, L, IV, direction){
  index=which(a!=0)
  index_l=rep(0, length(index))
  index_m=rep(0, length(index))
  for (i in 1:length(index)){
    count_l=L-1
    index_l[i]=1
    if (index[i]<L){
      index_m[i]=index[i]+1
    }
    else {
      while(index[i]>count_l){
        index_l[i]=index_l[i]+1
        count_l=count_l+L-index_l[i]
      }
      count_l=count_l-L+index_l[i]
      index_m[i]=index_l[i]+index[i]-count_l
    }
  }
  treated_index=rep(0, length(index))
  control_index=rep(0, length(index))
  if (direction==1){
    for (i in 1:length(index)){
      if (IV[index_l[i]]>IV[index_m[i]]){
        treated_index[i]=index_l[i]
        control_index[i]=index_m[i]
      }
      else {
        treated_index[i]=index_m[i]
        control_index[i]=index_l[i]
      }
    }
  }
  else {
    for (i in 1:length(index)){
      if (IV[index_l[i]]<IV[index_m[i]]){
        treated_index[i]=index_l[i]
        control_index[i]=index_m[i]
      }
      else {
        treated_index[i]=index_m[i]
        control_index[i]=index_l[i]
      }
    }
  }
  A=rbind(treated_index, control_index)
  return(A)
}



data.names     <- paste("C:/Users/siyuheng/Desktop/Real Data Analysis for Strengthening IV/Data/Preemie Data/preemies2005.csv")
matches.names  <- paste("C:/Users/siyuheng/Desktop/Real Data Analysis for Strengthening IV/Data/Debiased IV/incomplete_matched_debiasedIV2005.csv")

#for(k in 1:length(data.names)){
########
##data##
########

preemiedata=read.csv(data.names)
dim(preemiedata)


#####################
##Adjustable Inputs##
#####################
nomatches <- ceiling(dim(preemiedata)[1]/2)       # Total number of matches that will be made. If you don't adjust this number you'll have more rows in your output than are necessary.

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

X.reduced=cbind(bthwght,gestage,gestdiabetes,singlebirth,parity,ageM);

#create high level NICU indicator
high_level_NICU <- preemiedata$vol_level_big_2500
# Set up strata on race of mother, age of mother and gestational age and match within strata
dif.travel.time=preemiedata$diff_travel_2500;
losI = preemiedata$losI
death = preemiedata$death
death_fatel = preemiedata$death_fetal


########################
##Initialize Variables##
########################

X.names <- colnames(X.reduced)

varnames=c("dif.travel.time.e",paste(X.names,sep="",".e"),"high_level_NICU.e", "losI.e", "death.e", "death_fatel.e",
           "dif.travel.time.u",paste(X.names,sep="",".u"),"high_level_NICU.u", "losI.u", "death.u", "death_fatel.u")

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
  X.reduced.select=X.reduced[select1row, ]
  dif.travel.time.select=dif.travel.time[select1row]
  pair_indicator<-IV_debias(X.reduced[select1row, ], dif.travel.time.select, cov_dif, IV_diff)
  matched<-translate_index(pair_indicator, length(select1row), dif.travel.time.select, near_far)
  encouraged.numbers <- matched[1,]
  encouraged.match   <- matched[2,]
  # Save matches
  matchdataframe <- store.matches(matchdataframe,encouraged.numbers, encouraged.match, select1row, dif.travel.time, X.reduced, high_level_NICU, death, death_fatel, losI, matchcount)
  matchcount <- length(encouraged.numbers) + matchcount
  print(matchcount)
}
#do for the last birth weight stratum
i <- num.strata
select1=as.logical((gestage<=28)*(is.na(dif.travel.time)==FALSE)*(bthwght>=quans[i])*(bthwght<quans[i+1]));
select1row=which(select1)
#create matches
X.reduced.select=X.reduced[select1row, ]
dif.travel.time.select=dif.travel.time[select1row]
pair_indicator<-IV_debias(X.reduced[select1row, ], dif.travel.time.select, cov_dif, IV_diff)
matched<-translate_index(pair_indicator, length(select1row), dif.travel.time.select, near_far)
encouraged.numbers <- matched[1,]
encouraged.match   <- matched[2,]
# Save matches
matchdataframe <- store.matches(matchdataframe,encouraged.numbers, encouraged.match, select1row, dif.travel.time, X.reduced, high_level_NICU, death, death_fatel, losI, matchcount)
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
  X.reduced.select=X.reduced[select1row, ]
  dif.travel.time.select=dif.travel.time[select1row]
  pair_indicator<-IV_debias(X.reduced[select1row, ], dif.travel.time.select, cov_dif, IV_diff)
  matched<-translate_index(pair_indicator, length(select1row), dif.travel.time.select, near_far)
  encouraged.numbers <- matched[1,]
  encouraged.match   <- matched[2,]
  # Save matches
  matchdataframe <- store.matches(matchdataframe,encouraged.numbers, encouraged.match, select1row, dif.travel.time, X.reduced, high_level_NICU, death, death_fatel, losI, matchcount)
  matchcount <- length(encouraged.numbers) + matchcount
  print(matchcount)
}
#do for the last birth weight stratum
i <- num.strata
select1=as.logical((gestage>28)*(gestage<=32)*(is.na(dif.travel.time)==FALSE)*(bthwght>=quans[i])*(bthwght<quans[i+1]));
select1row=which(select1)
#create matches
X.reduced.select=X.reduced[select1row, ]
dif.travel.time.select=dif.travel.time[select1row]
pair_indicator<-IV_debias(X.reduced[select1row, ], dif.travel.time.select, cov_dif, IV_diff)
matched<-translate_index(pair_indicator, length(select1row), dif.travel.time.select, near_far)
encouraged.numbers <- matched[1,]
encouraged.match   <- matched[2,]
# Save matches
matchdataframe <- store.matches(matchdataframe,encouraged.numbers, encouraged.match, select1row, dif.travel.time, X.reduced, high_level_NICU, death, death_fatel, losI, matchcount)
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
  X.reduced.select=X.reduced[select1row, ]
  dif.travel.time.select=dif.travel.time[select1row]
  pair_indicator<-IV_debias(X.reduced[select1row, ], dif.travel.time.select, cov_dif, IV_diff)
  matched<-translate_index(pair_indicator, length(select1row), dif.travel.time.select, near_far)
  encouraged.numbers <- matched[1,]
  encouraged.match   <- matched[2,]
  # Save matches
  matchdataframe <- store.matches(matchdataframe,encouraged.numbers, encouraged.match, select1row, dif.travel.time, X.reduced, high_level_NICU, death, death_fatel, losI, matchcount)
  matchcount <- length(encouraged.numbers) + matchcount
  print(matchcount)
}
i <- num.strata
select1=as.logical((gestage>32)*(gestage<37)*(is.na(dif.travel.time)==FALSE)*(bthwght>=quans[i])*(bthwght<=quans[i+1]));
select1row=which(select1)
#create matches
X.reduced.select=X.reduced[select1row, ]
dif.travel.time.select=dif.travel.time[select1row]
pair_indicator<-IV_debias(X.reduced[select1row, ], dif.travel.time.select, cov_dif, IV_diff)
matched<-translate_index(pair_indicator, length(select1row), dif.travel.time.select, near_far)
encouraged.numbers <- matched[1,]
encouraged.match   <- matched[2,]
# Save matches
matchdataframe <- store.matches(matchdataframe,encouraged.numbers, encouraged.match, select1row, dif.travel.time, X.reduced, high_level_NICU, death, death_fatel, losI, matchcount)
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
  X.reduced.select=X.reduced[select1row, ]
  dif.travel.time.select=dif.travel.time[select1row]
  pair_indicator<-IV_debias(X.reduced[select1row, ], dif.travel.time.select, cov_dif, IV_diff)
  matched<-translate_index(pair_indicator, length(select1row), dif.travel.time.select, near_far)
  encouraged.numbers <- matched[1,]
  encouraged.match   <- matched[2,]
  # Save matches
  matchdataframe <- store.matches(matchdataframe,encouraged.numbers, encouraged.match, select1row, dif.travel.time, X.reduced, high_level_NICU, death, death_fatel, losI, matchcount)
  matchcount <- length(encouraged.numbers) + matchcount
  print(matchcount)
}
#do for the last birth weight stratum
i <- num.strata
select1=as.logical((gestage==37)*(is.na(dif.travel.time)==FALSE)*(bthwght>=quans[i])*(bthwght<quans[i+1]));
select1row=which(select1)
#create matches
X.reduced.select=X.reduced[select1row, ]
dif.travel.time.select=dif.travel.time[select1row]
pair_indicator<-IV_debias(X.reduced[select1row, ], dif.travel.time.select, cov_dif, IV_diff)
matched<-translate_index(pair_indicator, length(select1row), dif.travel.time.select, near_far)
encouraged.numbers <- matched[1,]
encouraged.match   <- matched[2,]
# Save matches
matchdataframe <- store.matches(matchdataframe,encouraged.numbers, encouraged.match, select1row, dif.travel.time, X.reduced, high_level_NICU, death, death_fatel, losI, matchcount)
matchcount <- length(encouraged.numbers) + matchcount
matchcount


# Save data
write.csv(matchdataframe[1:nomatches,],file=matches.names)


##############Make the balance table for 2005 years###################
debiasedIV<-read.csv("C:/Users/siyuheng/Desktop/Real Data Analysis for Strengthening IV/Data/Debiased IV/incomplete_matched_debiasedIV2005.csv")
near_mean_debiased=rep(0,7)
far_mean_debiased=rep(0,7)
for (i in 1:length(near_mean_debiased)){
  near_mean_debiased[i]=mean(debiasedIV[,1+i], na.rm = TRUE) #Check the number of the column every time!!!
  far_mean_debiased[i]=mean(debiasedIV[,12+i], na.rm = TRUE)
}

cov_dif_debiased = abs(near_mean_debiased[2:7]-far_mean_debiased[2:7]) ###Check the order of the covariates!
IV_diff_debiased = 2*abs(near_mean_debiased[1]-far_mean_debiased[1])


cov_dif
cov_dif_debiased
IV_diff
IV_diff_debiased