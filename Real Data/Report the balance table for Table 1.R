unmatched_data<-read.csv("C:/Users/siyuheng/Desktop/Real Data Analysis for Strengthening IV/Data/Preemie Data/allyears_preemie_data_imputed.csv")
matched_data<-read.csv("C:/Users/siyuheng/Desktop/Real Data Analysis for Strengthening IV/Data/Original IV/allyears_matched_originalIV.csv")
balance_table<-matrix(0, nrow = 11, ncol = 4)
colnames(balance_table)<-c("near mean", "far mean", "std", "std diff")
rownames(balance_table)<-colnames(unmatched_data)[c(2, 7:16)]
balance_table[1,3]=sqrt((var(matched_data$dif.travel.time.e, na.rm = TRUE)+var(matched_data$dif.travel.time.u, na.rm = TRUE))/2)
for (i in 2:11){
  balance_table[i,3]=sqrt((var(matched_data[,i+2], na.rm = TRUE)+var(matched_data[,i+17], na.rm = TRUE))/2)
}
for (i in 1:11){
  balance_table[i,1]=mean(matched_data[,i+2], na.rm = TRUE)
  balance_table[i,2]=mean(matched_data[,i+17], na.rm = TRUE)
}
balance_table[,4]=abs(balance_table[,1]-balance_table[,2])/balance_table[,3]
sample_size=sum(!is.na(matched_data$bthwght.e))

strength<-abs(mean(matched_data$high_level_NICU.e, na.rm = TRUE)-mean(matched_data$high_level_NICU.u, na.rm = TRUE))

matched_data_st<-read.csv("C:/Users/siyuheng/Desktop/Real Data Analysis for Strengthening IV/Data/Stronger IV/allyears_matched_strongerIV.csv")
balance_table_st<-matrix(0, nrow = 11, ncol = 4)
colnames(balance_table)<-c("near mean", "far mean", "std", "std diff")
rownames(balance_table)<-colnames(unmatched_data)[c(2, 7:16)]
balance_table_st[1,3]=sqrt((var(matched_data_st$dif.travel.time.e, na.rm = TRUE)+var(matched_data_st$dif.travel.time.u, na.rm = TRUE))/2)
for (i in 2:11){
  balance_table_st[i,3]=sqrt((var(matched_data_st[,i+2], na.rm = TRUE)+var(matched_data_st[,i+17], na.rm = TRUE))/2)
}
for (i in 1:11){
  balance_table_st[i,1]=mean(matched_data_st[,i+2], na.rm = TRUE)
  balance_table_st[i,2]=mean(matched_data_st[,i+17], na.rm = TRUE)
}
balance_table_st[,4]=abs(balance_table_st[,1]-balance_table_st[,2])/balance_table_st[,3]
sample_size_st=sum(!is.na(matched_data_st$bthwght.e))
strength_st<-abs(mean(matched_data_st$high_level_NICU.e, na.rm = TRUE)-mean(matched_data_st$high_level_NICU.u, na.rm = TRUE))

ratio=(sample_size/sample_size_st)*(strength/strength_st)^2

