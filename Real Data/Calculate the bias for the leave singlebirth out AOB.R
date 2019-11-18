
##########AOB for original IV#####
originalIV<-read.csv("C:/Users/siyuheng/Desktop/Real Data Analysis for Strengthening IV/Data/Original IV no single birth/allyears_matched_originalIV_no_singlebirth.csv")
diff_U_original=abs(mean(originalIV$singlebirth.e, na.rm = TRUE)-mean(originalIV$singlebirth.u, na.rm = TRUE))
diff_D_original=abs(mean(originalIV$high_level_NICU.e, na.rm = TRUE)-mean(originalIV$high_level_NICU.u, na.rm = TRUE))
AOB_original=diff_U_original/diff_D_original

#########AOB for strengthened IV#########
strongerIV<-read.csv("C:/Users/siyuheng/Desktop/Real Data Analysis for Strengthening IV/Data/Stronger IV no single birth/allyears_matched_strongerIV_no_singlebirth.csv")
diff_U_stronger=abs(mean(strongerIV$singlebirth.e, na.rm = TRUE)-mean(strongerIV$singlebirth.u, na.rm = TRUE))
diff_D_stronger=abs(mean(strongerIV$high_level_NICU.e, na.rm = TRUE)-mean(strongerIV$high_level_NICU.u, na.rm = TRUE))
AOB_stronger=diff_U_stronger/diff_D_stronger

ratio=AOB_stronger/AOB_original
