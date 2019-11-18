unmatched_data<-read.csv("C:/Users/siyuheng/Desktop/Real Data Analysis for Strengthening IV/Data/Preemie Data/allyears_preemie_data_imputed.csv")
matched_data_0<-read.csv("C:/Users/siyuheng/Desktop/Real Data Analysis for Strengthening IV/Data/Debiased IV/incomplete_matched_originalIV2005.csv")
matched_data_1<-read.csv("C:/Users/siyuheng/Desktop/Real Data Analysis for Strengthening IV/Data/Debiased IV/incomplete_matched_strengthenedIV2005.csv")
matched_data_2<-read.csv("C:/Users/siyuheng/Desktop/Real Data Analysis for Strengthening IV/Data/Debiased IV/incomplete_matched_debiasedIV2005.csv")

balance_table<-matrix(0, nrow = 6, ncol = 6)
colnames(balance_table)<-c("vanilla 1", "vanilla 2", "strengthened 1", "strengthened 2", "debiased 1", "debiased 2")
rownames(balance_table)<-colnames(unmatched_data)[7:12]

std_0<-rep(0,7)
std_1<-rep(0,7)
std_2<-rep(0,7)

for (i in 1:7){
  std_0[i]=sqrt((var(matched_data_0[,i+1], na.rm = TRUE)+var(matched_data_0[,i+12], na.rm = TRUE))/2)
}

for (i in 1:7){
  std_1[i]=sqrt((var(matched_data_1[,i+1], na.rm = TRUE)+var(matched_data_1[,i+12], na.rm = TRUE))/2)
}

for (i in 1:7){
  std_2[i]=sqrt((var(matched_data_2[,i+1], na.rm = TRUE)+var(matched_data_2[,i+12], na.rm = TRUE))/2)
}

com_rate<-rep(0,3)
com_rate[1]<-abs(mean(matched_data_0$high_level_NICU.e, na.rm = TRUE)-mean(matched_data_0$high_level_NICU.u, na.rm = TRUE))
com_rate[2]<-abs(mean(matched_data_1$high_level_NICU.e, na.rm = TRUE)-mean(matched_data_1$high_level_NICU.u, na.rm = TRUE))
com_rate[3]<-abs(mean(matched_data_2$high_level_NICU.e, na.rm = TRUE)-mean(matched_data_2$high_level_NICU.u, na.rm = TRUE))
for (i in 1:6){
  balance_table[i,1]=abs(mean(matched_data_0[,i+2], na.rm = TRUE)-mean(matched_data_0[,i+13], na.rm = TRUE))/std_0[i+1]
  balance_table[i,2]=abs(mean(matched_data_0[,i+2], na.rm = TRUE)-mean(matched_data_0[,i+13], na.rm = TRUE))/(std_0[i+1]*com_rate[1])
  balance_table[i,3]=abs(mean(matched_data_1[,i+2], na.rm = TRUE)-mean(matched_data_1[,i+13], na.rm = TRUE))/std_1[i+1]
  balance_table[i,4]=abs(mean(matched_data_1[,i+2], na.rm = TRUE)-mean(matched_data_1[,i+13], na.rm = TRUE))/(std_1[i+1]*com_rate[2])
  balance_table[i,5]=abs(mean(matched_data_2[,i+2], na.rm = TRUE)-mean(matched_data_2[,i+13], na.rm = TRUE))/std_2[i+1]
  balance_table[i,6]=abs(mean(matched_data_2[,i+2], na.rm = TRUE)-mean(matched_data_2[,i+13], na.rm = TRUE))/(std_2[i+1]*com_rate[3])
}
num_pairs<-rep(0,3)
num_pairs[1]=sum(!is.na(matched_data_0$bthwght.e))
num_pairs[2]=sum(!is.na(matched_data_1$bthwght.e))
num_pairs[3]=sum(!is.na(matched_data_2$bthwght.e))
balance_travel_time<-rep(0, 3)
balance_travel_time[1]=abs(mean(matched_data_0$dif.travel.time.e, na.rm = TRUE)-mean(matched_data_0$dif.travel.time.u, na.rm = TRUE))/std_0[1]
balance_travel_time[2]=abs(mean(matched_data_1$dif.travel.time.e, na.rm = TRUE)-mean(matched_data_1$dif.travel.time.u, na.rm = TRUE))/std_1[1]
balance_travel_time[3]=abs(mean(matched_data_2$dif.travel.time.e, na.rm = TRUE)-mean(matched_data_2$dif.travel.time.u, na.rm = TRUE))/std_2[1]
