# test for the overall effects

library(dplyr)
library(tidyr)
#library(tidyverse)
library(data.table)
library(Hmisc)
library(haven)
library(survival)		#Survival analysis
#library(survminer) 		#Survival analysis
library(rms)
library(survMisc)
library(Greg)
library(ggplot2)

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


# omega-3

## model 1

#######################
## omega3_pctg + strata(ageg) + strata(sex)
#######################

m1_anova_o3 = function(dataset) {
  fit1 = coxph(Surv(time_to_cancer_year, status)~ omega3_pctg + strata(ageg) + strata(sex), data = dataset)
  fit2 = coxph(Surv(time_to_cancer_year, status)~ strata(ageg) + strata(sex), data = dataset)
  res = round(anova(fit1, fit2)[2,4],3)
  return(res)
}

n = 20
datalist = list()

for (i in 1:n) {
  dat <- m1_anova_o3(dataset = data[[i]])
  datalist[[i]] <- dat # add it to your list
}
m1_anova_o3_r <- bind_cols(datalist)
m1_anova_o3_r <- as.data.frame(t(as.matrix(m1_anova_o3_r)))
m1_anova_o3_r$dataid = c(1:20)
m1_anova_o3_r$model = 1

## model 2
#######################
# m2_cat_o3
## omega3_pctg + strata(ageg) + strata(sex) + TDI + eth_group_new + BMI + Smoking_status + Alcohol_status + IPAQ_activity
#######################

m2_anova_o3 = function(dataset) {
  fit1 = coxph(Surv(time_to_cancer_year, status)~ omega3_pctg + strata(ageg) + strata(sex) + TDI + eth_group_new + BMI + Smoking_status + Alcohol_status + IPAQ_activity, data = dataset)
  fit2 = coxph(Surv(time_to_cancer_year, status)~ strata(ageg) + strata(sex) + TDI + eth_group_new + BMI + Smoking_status + Alcohol_status + IPAQ_activity, data = dataset)
  res = round(anova(fit1, fit2)[2,4],3)
  return(res)
}

n = 20
datalist = list()

for (i in 1:n) {
  dat <- m2_anova_o3(dataset = data[[i]])
  datalist[[i]] <- dat # add it to your list
}
m2_anova_o3_r <- bind_cols(datalist)
m2_anova_o3_r <- as.data.frame(t(as.matrix(m2_anova_o3_r)))
m2_anova_o3_r$dataid = c(1:20)
m2_anova_o3_r$model = 2

# merge the results
anova_o3 = rbind(m1_anova_o3_r, m2_anova_o3_r)
anova_o3 = anova_o3 %>% arrange(dataid)

write.csv(anova_o3, file="/Users/yuchen/Desktop/Cancer/Data/anova_o3.csv",row.names = F)

# omega-6

## model 1

#######################
## omega6_pctg + strata(ageg) + strata(sex)
#######################

m1_anova_o6 = function(dataset) {
  fit1 = coxph(Surv(time_to_cancer_year, status)~ omega6_pctg + strata(ageg) + strata(sex), data = dataset)
  fit2 = coxph(Surv(time_to_cancer_year, status)~ strata(ageg) + strata(sex), data = dataset)
  res = round(anova(fit1, fit2)[2,4],3)
  return(res)
}

n = 20
datalist = list()

for (i in 1:n) {
  dat <- m1_anova_o6(dataset = data[[i]])
  datalist[[i]] <- dat # add it to your list
}
m1_anova_o6_r <- bind_cols(datalist)
m1_anova_o6_r <- as.data.frame(t(as.matrix(m1_anova_o6_r)))
m1_anova_o6_r$dataid = c(1:20)
m1_anova_o6_r$model = 1

## model 2
#######################
# m2_cat_o6
## omega6_pctg + strata(ageg) + strata(sex) + TDI + eth_group_new + BMI + Smoking_status + Alcohol_status + IPAQ_activity
#######################

m2_anova_o6 = function(dataset) {
  fit1 = coxph(Surv(time_to_cancer_year, status)~ omega6_pctg + strata(ageg) + strata(sex) + TDI + eth_group_new + BMI + Smoking_status + Alcohol_status + IPAQ_activity, data = dataset)
  fit2 = coxph(Surv(time_to_cancer_year, status)~ strata(ageg) + strata(sex) + TDI + eth_group_new + BMI + Smoking_status + Alcohol_status + IPAQ_activity, data = dataset)
  res = round(anova(fit1, fit2)[2,4],3)
  return(res)
}

n = 20
datalist = list()

for (i in 1:n) {
  dat <- m2_anova_o6(dataset = data[[i]])
  datalist[[i]] <- dat # add it to your list
}
m2_anova_o6_r <- bind_cols(datalist)
m2_anova_o6_r <- as.data.frame(t(as.matrix(m2_anova_o6_r)))
m2_anova_o6_r$dataid = c(1:20)
m2_anova_o6_r$model = 2

# merge the results
anova_o6 = rbind(m1_anova_o6_r, m2_anova_o6_r)
anova_o6 = anova_o6 %>% arrange(dataid)

write.csv(anova_o6, file="/Users/yuchen/Desktop/Cancer/Data/anova_o6.csv",row.names = F)


### FDR adjusted p-values

# omega-3
m1_anova_o3_ad = m1_anova_o3_r[-1,]
m1_anova_o3_ad$p_ad = p.adjust(m1_anova_o3_ad$V1, method = "fdr", n = 19)
m1_anova_o3_ad$p_ad = round(m1_anova_o3_ad$p_ad,3)

m2_anova_o3_ad = m2_anova_o3_r[-1,]
m2_anova_o3_ad$p_ad = p.adjust(m2_anova_o3_ad$V1, method = "fdr", n = 19)
m2_anova_o3_ad$p_ad = round(m2_anova_o3_ad$p_ad,3)

# merge the results
anova_o3_ad = rbind(m1_anova_o3_ad, m2_anova_o3_ad)
anova_o3_ad = anova_o3_ad %>% arrange(dataid)

write.csv(anova_o3_ad, file="/Users/yuchen/Desktop/Cancer/Data/anova_o3_ad.csv",row.names = F)

# omega-6

m1_anova_o6_ad = m1_anova_o6_r[-1,]
m1_anova_o6_ad$p_ad = p.adjust(m1_anova_o6_ad$V1, method = "fdr", n = 19)
m1_anova_o6_ad$p_ad = round(m1_anova_o6_ad$p_ad,3)

m2_anova_o6_ad = m2_anova_o6_r[-1,]
m2_anova_o6_ad$p_ad = p.adjust(m2_anova_o6_ad$V1, method = "fdr", n = 19)
m2_anova_o6_ad$p_ad = round(m2_anova_o6_ad$p_ad,3)

# merge the results
anova_o6_ad = rbind(m1_anova_o6_ad, m2_anova_o6_ad)
anova_o6_ad = anova_o6_ad %>% arrange(dataid)

write.csv(anova_o6_ad, file="/Users/yuchen/Desktop/Cancer/Data/anova_o6_ad.csv",row.names = F)

### Additionally adjusted model

# 2. Esophagus
## + gas_dises_yn + wh_ratio

### omega-3 0.187
fit2_o3 = coxph(Surv(time_to_cancer_year, status)~ omega3_pctg + strata(ageg) + strata(sex) + TDI + eth_group_new + BMI + Smoking_status + 
                          Alcohol_status + IPAQ_activity + gas_dises_yn + wh_ratio, data = data2)
fit2_o3r = coxph(Surv(time_to_cancer_year, status)~ strata(ageg) + strata(sex) + TDI + eth_group_new + BMI + Smoking_status + 
                  Alcohol_status + IPAQ_activity + gas_dises_yn + wh_ratio, data = data2)
anova(fit2_o3,fit2_o3r)

### omega-6 0.006
fit2_o6 = coxph(Surv(time_to_cancer_year, status)~ omega6_pctg + strata(ageg) + strata(sex) + TDI + eth_group_new + BMI + Smoking_status + 
                  Alcohol_status + IPAQ_activity + gas_dises_yn + wh_ratio, data = data2)
fit2_o6r = coxph(Surv(time_to_cancer_year, status)~ strata(ageg) + strata(sex) + TDI + eth_group_new + BMI + Smoking_status + 
                   Alcohol_status + IPAQ_activity + gas_dises_yn + wh_ratio, data = data2)
anova(fit2_o6,fit2_o6r)

##----------------------------------------------------------------------------------------------------------##

# 4.	Colon
## + diabet + aspir_yn + meat_yn + wh_ratio + fh_bowel

####### continuous
### omega-3 0.009
fit4_o3 = coxph(Surv(time_to_cancer_year, status)~ omega3_pctg + strata(ageg) + strata(sex) + TDI + eth_group_new + BMI + Smoking_status + 
                          Alcohol_status + IPAQ_activity + diabet + aspir_yn + meat_yn + wh_ratio + fh_bowel, data = data4)
fit4_o3r = coxph(Surv(time_to_cancer_year, status)~ strata(ageg) + strata(sex) + TDI + eth_group_new + BMI + Smoking_status + 
                  Alcohol_status + IPAQ_activity + diabet + aspir_yn + meat_yn + wh_ratio + fh_bowel, data = data4)
anova(fit4_o3,fit4_o3r)

### omega-6 0.105
fit4_o6 = coxph(Surv(time_to_cancer_year, status)~ omega6_pctg + strata(ageg) + strata(sex) + TDI + eth_group_new + BMI + Smoking_status + 
                  Alcohol_status + IPAQ_activity + diabet + aspir_yn + meat_yn + wh_ratio + fh_bowel, data = data4)
fit4_o6r = coxph(Surv(time_to_cancer_year, status)~ strata(ageg) + strata(sex) + TDI + eth_group_new + BMI + Smoking_status + 
                   Alcohol_status + IPAQ_activity + diabet + aspir_yn + meat_yn + wh_ratio + fh_bowel, data = data4)
anova(fit4_o6,fit4_o6r)


##----------------------------------------------------------------------------------------------------------##

# 5.	Rectum
## + diabet + aspir_yn + meat_yn + wh_ratio + fh_bowel

####### continuous
### omega-3 0.475
fit5_o3 = coxph(Surv(time_to_cancer_year, status)~ omega3_pctg + strata(ageg) + strata(sex) + TDI + eth_group_new + BMI + Smoking_status + 
                  Alcohol_status + IPAQ_activity + diabet + aspir_yn + meat_yn + wh_ratio + fh_bowel, data = data5)
fit5_o3r = coxph(Surv(time_to_cancer_year, status)~ strata(ageg) + strata(sex) + TDI + eth_group_new + BMI + Smoking_status + 
                   Alcohol_status + IPAQ_activity + diabet + aspir_yn + meat_yn + wh_ratio + fh_bowel, data = data5)
anova(fit5_o3,fit5_o3r)

### omega-6 0.125
fit5_o6 = coxph(Surv(time_to_cancer_year, status)~ omega6_pctg + strata(ageg) + strata(sex) + TDI + eth_group_new + BMI + Smoking_status + 
                  Alcohol_status + IPAQ_activity + diabet + aspir_yn + meat_yn + wh_ratio + fh_bowel, data = data5)
fit5_o6r = coxph(Surv(time_to_cancer_year, status)~ strata(ageg) + strata(sex) + TDI + eth_group_new + BMI + Smoking_status + 
                   Alcohol_status + IPAQ_activity + diabet + aspir_yn + meat_yn + wh_ratio + fh_bowel, data = data5)
anova(fit5_o6,fit5_o6r)

##----------------------------------------------------------------------------------------------------------##

# 7.	Pancreas
## + diabet

### omega-3 0.312
fit7_o3 = coxph(Surv(time_to_cancer_year, status)~ omega3_pctg + strata(ageg) + strata(sex) + TDI + eth_group_new + BMI + Smoking_status + 
                          Alcohol_status + IPAQ_activity + diabet, data = data7)
fit7_o3r = coxph(Surv(time_to_cancer_year, status)~ strata(ageg) + strata(sex) + TDI + eth_group_new + BMI + Smoking_status + 
                  Alcohol_status + IPAQ_activity + diabet, data = data7)
anova(fit7_o3,fit7_o3r)

### omega-6 0.012
fit7_o6 = coxph(Surv(time_to_cancer_year, status)~ omega6_pctg + strata(ageg) + strata(sex) + TDI + eth_group_new + BMI + Smoking_status + 
                  Alcohol_status + IPAQ_activity + diabet, data = data7)
fit7_o6r = coxph(Surv(time_to_cancer_year, status)~ strata(ageg) + strata(sex) + TDI + eth_group_new + BMI + Smoking_status + 
                   Alcohol_status + IPAQ_activity + diabet, data = data7)
anova(fit7_o6,fit7_o6r)


##----------------------------------------------------------------------------------------------------------##

# 8.	Lung
## + fh_lung

### omega-3 0.001
fit8_o3 = coxph(Surv(time_to_cancer_year, status)~ omega3_pctg + strata(ageg) + strata(sex) + TDI + eth_group_new + BMI + Smoking_status + 
                          Alcohol_status + IPAQ_activity + fh_lung, data = data8)
fit8_o3r = coxph(Surv(time_to_cancer_year, status)~ strata(ageg) + strata(sex) + TDI + eth_group_new + BMI + Smoking_status + 
                  Alcohol_status + IPAQ_activity + fh_lung, data = data8)
anova(fit8_o3,fit8_o3r)

### omega-6  0.000
fit8_o6 = coxph(Surv(time_to_cancer_year, status)~ omega6_pctg + strata(ageg) + strata(sex) + TDI + eth_group_new + BMI + Smoking_status + 
                  Alcohol_status + IPAQ_activity + fh_lung, data = data8)
fit8_o6r = coxph(Surv(time_to_cancer_year, status)~ strata(ageg) + strata(sex) + TDI + eth_group_new + BMI + Smoking_status + 
                   Alcohol_status + IPAQ_activity + fh_lung, data = data8)
anova(fit8_o6,fit8_o6r)


##----------------------------------------------------------------------------------------------------------##

# 9.	Malignant melanoma
## + skin_color + bright_sun + uv_protect + sunburn + sunlamp

### omega-3 0.686
fit9_o3 = coxph(Surv(time_to_cancer_year, status)~ omega3_pctg + strata(ageg) + strata(sex) + TDI + eth_group_new + BMI + Smoking_status + 
                          Alcohol_status + IPAQ_activity + skin_color + bright_sun + uv_protect + sunburn + sunlamp, data = data9)
fit9_o3r = coxph(Surv(time_to_cancer_year, status)~ strata(ageg) + strata(sex) + TDI + eth_group_new + BMI + Smoking_status + 
                  Alcohol_status + IPAQ_activity + skin_color + bright_sun + uv_protect + sunburn + sunlamp, data = data9)
anova(fit9_o3,fit9_o3r)

### omega-6 0.014
fit9_o6 = coxph(Surv(time_to_cancer_year, status)~ omega6_pctg + strata(ageg) + strata(sex) + TDI + eth_group_new + BMI + Smoking_status + 
                  Alcohol_status + IPAQ_activity + skin_color + bright_sun + uv_protect + sunburn + sunlamp, data = data9)
fit9_o6r = coxph(Surv(time_to_cancer_year, status)~ strata(ageg) + strata(sex) + TDI + eth_group_new + BMI + Smoking_status + 
                   Alcohol_status + IPAQ_activity + skin_color + bright_sun + uv_protect + sunburn + sunlamp, data = data9)
anova(fit9_o6,fit9_o6r)

##----------------------------------------------------------------------------------------------------------##

# 11.	Breast
## + HRT + OC + num_birth + age_mena + meno + hyster + fh_breast

### omega-3 0.108
fit11_o3 = coxph(Surv(time_to_cancer_year, status)~ omega3_pctg + strata(ageg) + strata(sex) + TDI + eth_group_new + BMI + Smoking_status + 
                           Alcohol_status + IPAQ_activity + HRT + OC + num_birth + age_mena + meno + hyster + fh_breast, data = data11)
fit11_o3r = coxph(Surv(time_to_cancer_year, status)~ strata(ageg) + strata(sex) + TDI + eth_group_new + BMI + Smoking_status + 
                   Alcohol_status + IPAQ_activity + HRT + OC + num_birth + age_mena + meno + hyster + fh_breast, data = data11)
anova(fit11_o3,fit11_o3r)

### omega-6 0.991
fit11_o6 = coxph(Surv(time_to_cancer_year, status)~ omega6_pctg + strata(ageg) + strata(sex) + TDI + eth_group_new + BMI + Smoking_status + 
                   Alcohol_status + IPAQ_activity + HRT + OC + num_birth + age_mena + meno + hyster + fh_breast, data = data11)
fit11_o6r = coxph(Surv(time_to_cancer_year, status)~ strata(ageg) + strata(sex) + TDI + eth_group_new + BMI + Smoking_status + 
                    Alcohol_status + IPAQ_activity + HRT + OC + num_birth + age_mena + meno + hyster + fh_breast, data = data11)
anova(fit11_o6,fit11_o6r)

##----------------------------------------------------------------------------------------------------------##

# 12.	Uterus
## + HRT + OC + num_birth + age_mena + meno + hyster

### omega-3 0.817
fit12_o3 = coxph(Surv(time_to_cancer_year, status)~ omega3_pctg + strata(ageg) + strata(sex) + TDI + eth_group_new + BMI + Smoking_status + 
                           Alcohol_status + IPAQ_activity + HRT + OC + num_birth + age_mena + meno + hyster, data = data12)
fit12_o3r = coxph(Surv(time_to_cancer_year, status)~ strata(ageg) + strata(sex) + TDI + eth_group_new + BMI + Smoking_status + 
                   Alcohol_status + IPAQ_activity + HRT + OC + num_birth + age_mena + meno + hyster, data = data12)
anova(fit12_o3,fit12_o3r)

### omega-6 0.022
fit12_o6 = coxph(Surv(time_to_cancer_year, status)~ omega6_pctg + strata(ageg) + strata(sex) + TDI + eth_group_new + BMI + Smoking_status + 
                   Alcohol_status + IPAQ_activity + HRT + OC + num_birth + age_mena + meno + hyster, data = data12)
fit12_o6r = coxph(Surv(time_to_cancer_year, status)~ strata(ageg) + strata(sex) + TDI + eth_group_new + BMI + Smoking_status + 
                    Alcohol_status + IPAQ_activity + HRT + OC + num_birth + age_mena + meno + hyster, data = data12)
anova(fit12_o6,fit12_o6r)

##----------------------------------------------------------------------------------------------------------##

# 13.	Ovary
## + HRT + OC + num_birth + age_mena + meno + hyster

### omega-3 0.105
fit13_o3 = coxph(Surv(time_to_cancer_year, status)~ omega3_pctg + strata(ageg) + strata(sex) + TDI + eth_group_new + BMI + Smoking_status + 
                   Alcohol_status + IPAQ_activity + HRT + OC + num_birth + age_mena + meno + hyster, data = data13)
fit13_o3r = coxph(Surv(time_to_cancer_year, status)~ strata(ageg) + strata(sex) + TDI + eth_group_new + BMI + Smoking_status + 
                    Alcohol_status + IPAQ_activity + HRT + OC + num_birth + age_mena + meno + hyster, data = data13)
anova(fit13_o3,fit13_o3r)

### omega-6 0.151
fit13_o6 = coxph(Surv(time_to_cancer_year, status)~ omega6_pctg + strata(ageg) + strata(sex) + TDI + eth_group_new + BMI + Smoking_status + 
                   Alcohol_status + IPAQ_activity + HRT + OC + num_birth + age_mena + meno + hyster, data = data13)
fit13_o6r = coxph(Surv(time_to_cancer_year, status)~ strata(ageg) + strata(sex) + TDI + eth_group_new + BMI + Smoking_status + 
                    Alcohol_status + IPAQ_activity + HRT + OC + num_birth + age_mena + meno + hyster, data = data13)
anova(fit13_o6,fit13_o6r)

##----------------------------------------------------------------------------------------------------------##

# 14.	Prostate
## + fh_prostat

### omega-3 0.346
fit14_o3 = coxph(Surv(time_to_cancer_year, status)~ omega3_pctg + strata(ageg) + strata(sex) + TDI + eth_group_new + BMI + Smoking_status + 
                           Alcohol_status + IPAQ_activity + fh_prostat, data = data14)
fit14_o3r = coxph(Surv(time_to_cancer_year, status)~ strata(ageg) + strata(sex) + TDI + eth_group_new + BMI + Smoking_status + 
                   Alcohol_status + IPAQ_activity + fh_prostat, data = data14)
anova(fit14_o3,fit14_o3r)

### omega-6 0.005
fit14_o6 = coxph(Surv(time_to_cancer_year, status)~ omega6_pctg + strata(ageg) + strata(sex) + TDI + eth_group_new + BMI + Smoking_status + 
                   Alcohol_status + IPAQ_activity + fh_prostat, data = data14)
fit14_o6r = coxph(Surv(time_to_cancer_year, status)~ strata(ageg) + strata(sex) + TDI + eth_group_new + BMI + Smoking_status + 
                    Alcohol_status + IPAQ_activity + fh_prostat, data = data14)
anova(fit14_o6,fit14_o6r)




