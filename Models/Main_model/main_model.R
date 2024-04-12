library(dplyr)
library(tidyr)
# library(tidyverse)
library(data.table)
library(Hmisc)
library(haven)
library(survival)		#Survival analysis
# library(survminer) 		#Survival analysis
library(rms)
library(survMisc)
library(Greg)


data_add <- fread("/Users/yuchen/Desktop/Cancer/Data/data_add.csv", header=TRUE, sep=",")

data_add$Alcohol_status = as.factor(data_add$Alcohol_status)
data_add$Smoking_status = as.factor(data_add$Smoking_status)
data_add$IPAQ_activity = as.factor(data_add$IPAQ_activity)
data_add$ageg = as.factor(data_add$ageg)
data_add$sex = as.factor(data_add$sex)
data_add$Fish_oil_baseline = as.factor(data_add$Fish_oil_baseline)
data_add$omega3_pctg = as.factor(data_add$omega3_pctg)
data_add$omega6_pctg = as.factor(data_add$omega6_pctg)
data_add$omega_ratiog = as.factor(data_add$omega_ratiog)

data_add$omega3_pctstd = scale(data_add$omega3_pct,center=TRUE,scale=TRUE)
data_add$omega6_pctstd = scale(data_add$omega6_pct,center=TRUE,scale=TRUE)
data_add$omegar_std = scale(data_add$omega_ratio,center=TRUE,scale=TRUE)

data_add$skin_color = as.factor(data_add$skin_color)
data_add$bright_sun = as.factor(data_add$bright_sun)
data_add$uv_protect = as.factor(data_add$uv_protect)
data_add$sunburn = as.factor(data_add$sunburn)
data_add$sunlamp = as.factor(data_add$sunlamp)

data_add$HRT = as.factor(data_add$HRT)
data_add$OC = as.factor(data_add$OC)
data_add$meno = as.factor(data_add$meno)
data_add$hyster = as.factor(data_add$hyster)

data1 = data_add %>% filter(cancer_sub == "01_head" | is.na(cancer_sub) == TRUE)
data2 = data_add %>% filter(cancer_sub == "02_esopha" | is.na(cancer_sub) == TRUE)
data3 = data_add %>% filter(cancer_sub == "03_stomach" | is.na(cancer_sub) == TRUE)
data4 = data_add %>% filter(cancer_sub == "04_colon" | is.na(cancer_sub) == TRUE)
data5 = data_add %>% filter(cancer_sub == "05_rectum" | is.na(cancer_sub) == TRUE)
data6 = data_add %>% filter(cancer_sub == "06_hepatob" | is.na(cancer_sub) == TRUE)
data7 = data_add %>% filter(cancer_sub == "07_pancrea" | is.na(cancer_sub) == TRUE)
data8 = data_add %>% filter(cancer_sub == "08_lung" | is.na(cancer_sub) == TRUE)
data9 = data_add %>% filter(cancer_sub == "09_mal_melan" | is.na(cancer_sub) == TRUE)
data10 = data_add %>% filter(cancer_sub == "10_connect" | is.na(cancer_sub) == TRUE)
data11 = data_add %>% filter(cancer_sub == "11_breast" | is.na(cancer_sub) == TRUE)
data12 = data_add %>% filter(cancer_sub == "12_uterus" | is.na(cancer_sub) == TRUE)
data13 = data_add %>% filter(cancer_sub == "13_ovary" | is.na(cancer_sub) == TRUE)
data14 = data_add %>% filter(cancer_sub == "14_prostate" | is.na(cancer_sub) == TRUE)
data15 = data_add %>% filter(cancer_sub == "15_kidney" | is.na(cancer_sub) == TRUE)
data16 = data_add %>% filter(cancer_sub == "16_bladder" | is.na(cancer_sub) == TRUE)
data17 = data_add %>% filter(cancer_sub == "17_brain" | is.na(cancer_sub) == TRUE)
data18 = data_add %>% filter(cancer_sub == "18_thyroid" | is.na(cancer_sub) == TRUE)
data19 = data_add %>% filter(cancer_sub == "19_lymphoma" | is.na(cancer_sub) == TRUE)

data = list(data_add, data1, data2, data3, data4, data5, data6, data7, data8, data9, data10, data11,
            data12, data13, data14, data15, data16, data17, data18, data19)

# omega-3 model 1
#######################
# m1_5_o3
## omega3_pctg + strata(ageg) + strata(sex)
#######################

m1_5_o3 = function(dataset) {
  fit = summary(coxph(Surv(time_to_cancer_year, status)~ omegar_std + strata(ageg) + strata(sex), data = dataset))
  res = round(fit$coefficients[1,5],3)
  return(res)
}

n = 20
datalist = list()

for (i in 1:n) {
  dat <- m1_5_o3(dataset = data[[i]])
  datalist[[i]] <- dat # add it to your list
}
m1_5_o3_r <- bind_cols(datalist)
m1_5_o3_r <- as.data.frame(t(as.matrix(m1_5_o3_r)))
m1_5_o3_r$dataid = c(1:20)
m1_5_o3_r$model = 1


# omega-3 model 2
#######################
# m2_5_o3
## omega3_pctg + strata(ageg) + strata(sex) + TDI + eth_group_new + BMI + Smoking_status + Alcohol_status + IPAQ_activity
#######################

m2_5_o3 = function(dataset) {
  fit = summary(coxph(Surv(time_to_cancer_year, status)~ omegar_std + strata(ageg) + strata(sex) + TDI + eth_group_new + BMI + Smoking_status + Alcohol_status + IPAQ_activity, data = dataset))
  res = round(fit$coefficients[1,5],3)
  return(res)
}

n = 20
datalist = list()

for (i in 1:n) {
  dat <- m2_5_o3(dataset = data[[i]])
  datalist[[i]] <- dat # add it to your list
}
m2_5_o3_r <- bind_cols(datalist)
m2_5_o3_r <- as.data.frame(t(as.matrix(m2_5_o3_r)))
m2_5_o3_r$dataid = c(1:20)
m2_5_o3_r$model = 2


# merge the results
or_5 = rbind(m1_5_o3_r, m2_5_o3_r)
or_5 = or_5 %>% arrange(dataid)

write.csv(or_5, file="/Users/yuchen/Desktop/Cancer/Data/or_con_p.csv",row.names = F)



### FDR adjusted p-values

# omega-3
m1_5_o3_ad = m1_5_o3_r[-1,]
m1_5_o3_ad$p_ad = p.adjust(m1_5_o3_ad$V1, method = "fdr", n = 19)
m1_5_o3_ad$p_ad = round(m1_5_o3_ad$p_ad,3)

m2_5_o3_ad = m2_5_o3_r[-1,]
m2_5_o3_ad$p_ad = p.adjust(m2_5_o3_ad$V1, method = "fdr", n = 19)
m2_5_o3_ad$p_ad = round(m2_5_o3_ad$p_ad,3)

# merge the results
or_5_ad = rbind(m1_5_o3_ad, m2_5_o3_ad)
or_5_ad = or_5_ad %>% arrange(dataid)

write.csv(or_5_ad, file="/Users/yuchen/Desktop/Cancer/Data/or_conp_ad.csv",row.names = F)



### Additionally adjusted model

# 2. Esophagus
## + gas_dises_yn + wh_ratio

### omega-3
fit2_o3_cat = summary(coxph(Surv(time_to_cancer_year, status)~ omegar_std + strata(ageg) + strata(sex) + TDI + eth_group_new + BMI + Smoking_status + 
                              Alcohol_status + IPAQ_activity + gas_dises_yn + wh_ratio, data = data2))
cat2_o3_r = round(fit2_o3_cat$coefficients[1,5],3)


##----------------------------------------------------------------------------------------------------------##

# 4.	Colon
## + diabet + aspir_yn + meat_yn + wh_ratio + fh_bowel

### omega-3
fit4_o3_cat = summary(coxph(Surv(time_to_cancer_year, status)~ omegar_std + strata(ageg) + strata(sex) + TDI + eth_group_new + BMI + Smoking_status + 
                              Alcohol_status + IPAQ_activity + diabet + aspir_yn + meat_yn + wh_ratio + fh_bowel, data = data4))
cat4_o3_r = round(fit4_o3_cat$coefficients[1,5],3)




##----------------------------------------------------------------------------------------------------------##

# 5.	Rectum
## + diabet + aspir_yn + meat_yn + wh_ratio + fh_bowel

### omega-3
fit5_o3_cat = summary(coxph(Surv(time_to_cancer_year, status)~ omegar_std + strata(ageg) + strata(sex) + TDI + eth_group_new + BMI + Smoking_status + 
                              Alcohol_status + IPAQ_activity + diabet + aspir_yn + meat_yn + wh_ratio + fh_bowel, data = data5))
cat5_o3_r = round(fit5_o3_cat$coefficients[1,5],3)



##----------------------------------------------------------------------------------------------------------##

# 7.	Pancreas
## + diabet

### omega-3
fit7_o3_cat = summary(coxph(Surv(time_to_cancer_year, status)~ omegar_std + strata(ageg) + strata(sex) + TDI + eth_group_new + BMI + Smoking_status + 
                              Alcohol_status + IPAQ_activity + diabet, data = data7))
cat7_o3_r = round(fit7_o3_cat$coefficients[1,5],3)




##----------------------------------------------------------------------------------------------------------##

# 8.	Lung
## + fh_lung

### omega-3
fit8_o3_cat = summary(coxph(Surv(time_to_cancer_year, status)~ omegar_std + strata(ageg) + strata(sex) + TDI + eth_group_new + BMI + Smoking_status + 
                              Alcohol_status + IPAQ_activity + fh_lung, data = data8))
cat8_o3_r = round(fit8_o3_cat$coefficients[1,5],3)



##----------------------------------------------------------------------------------------------------------##

# 9.	Malignant melanoma
## + skin_color + bright_sun + uv_protect + sunburn + sunlamp

### omega-3
fit9_o3_cat = summary(coxph(Surv(time_to_cancer_year, status)~ omegar_std + strata(ageg) + strata(sex) + TDI + eth_group_new + BMI + Smoking_status + 
                              Alcohol_status + IPAQ_activity + skin_color + bright_sun + uv_protect + sunburn + sunlamp, data = data9))
cat9_o3_r = round(fit9_o3_cat$coefficients[1,5],3)




##----------------------------------------------------------------------------------------------------------##

# 11.	Breast
## + HRT + OC + num_birth + age_mena + meno + hyster + fh_breast


### omega-3
fit11_o3_cat = summary(coxph(Surv(time_to_cancer_year, status)~ omegar_std + strata(ageg) + strata(sex) + TDI + eth_group_new + BMI + Smoking_status + 
                               Alcohol_status + IPAQ_activity + HRT + OC + num_birth + age_mena + meno + hyster + fh_breast, data = data11))
cat11_o3_r = round(fit11_o3_cat$coefficients[1,5],3)





##----------------------------------------------------------------------------------------------------------##

# 12.	Uterus
## + HRT + OC + num_birth + age_mena + meno + hyster

### omega-3
fit12_o3_cat = summary(coxph(Surv(time_to_cancer_year, status)~ omegar_std + strata(ageg) + strata(sex) + TDI + eth_group_new + BMI + Smoking_status + 
                               Alcohol_status + IPAQ_activity + HRT + OC + num_birth + age_mena + meno + hyster, data = data12))
cat12_o3_r = round(fit12_o3_cat$coefficients[1,5],3)





##----------------------------------------------------------------------------------------------------------##

# 13.	Ovary
## + HRT + OC + num_birth + age_mena + meno + hyster

### omega-3
fit13_o3_cat = summary(coxph(Surv(time_to_cancer_year, status)~ omegar_std + strata(ageg) + strata(sex) + TDI + eth_group_new + BMI + Smoking_status + 
                               Alcohol_status + IPAQ_activity + HRT + OC + num_birth + age_mena + meno + hyster, data = data13))
cat13_o3_r = round(fit13_o3_cat$coefficients[1,5],3)




##----------------------------------------------------------------------------------------------------------##

# 14.	Prostate
## + fh_prostat

### omega-3
fit14_o3_cat = summary(coxph(Surv(time_to_cancer_year, status)~ omegar_std + strata(ageg) + strata(sex) + TDI + eth_group_new + BMI + Smoking_status + 
                               Alcohol_status + IPAQ_activity + fh_prostat, data = data14))
cat14_o3_r = round(fit14_o3_cat$coefficients[1,5],3)



##----------------------------------------------------------------------------------------------------------##


## Results summary


o3_5_sum = rbind(cat2_o3_r, cat4_o3_r, cat5_o3_r, cat7_o3_r, cat8_o3_r,
                 cat9_o3_r, cat11_o3_r, cat12_o3_r, cat13_o3_r, cat14_o3_r)

o6_5_sum = rbind(cat2_o6_r, cat4_o6_r, cat5_o6_r, cat7_o6_r, cat8_o6_r,
                 cat9_o6_r, cat11_o6_r, cat12_o6_r, cat13_o6_r, cat14_o6_r)

write.csv(o3_5_sum, file="/Users/yuchen/Desktop/Cancer/Data/o3_conp_add.csv",row.names = F)
write.csv(o6_5_sum, file="/Users/yuchen/Desktop/Cancer/Data/o6_conp_add.csv",row.names = F)



