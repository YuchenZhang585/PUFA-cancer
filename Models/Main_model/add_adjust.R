#################
## Additional adjusted
#################

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

data_add = data_add %>% mutate(omega3_pctg = ntile(omega3_pct, 5))
data_add = data_add %>% mutate(omega6_pctg = ntile(omega6_pct, 5))
data_add = data_add %>% mutate(omega_ratiog = ntile(omega_ratio, 5))

data_add$Alcohol_status = as.factor(data_add$Alcohol_status)
data_add$Smoking_status = as.factor(data_add$Smoking_status)
data_add$IPAQ_activity = as.factor(data_add$IPAQ_activity)
data_add$ageg = as.factor(data_add$ageg)
data_add$sex = as.factor(data_add$sex)
data_add$Fish_oil_baseline = as.factor(data_add$Fish_oil_baseline)
data_add$omega3_pctg = as.factor(data_add$omega3_pctg)
data_add$omega6_pctg = as.factor(data_add$omega6_pctg)
data_add$omega_ratiog = as.factor(data_add$omega_ratiog)

data_add = data_add %>% group_by(omega3_pctg) %>% mutate(omega3_m = median(omega3_pct))
data_add = data_add %>% group_by(omega6_pctg) %>% mutate(omega6_m = median(omega6_pct))
data_add = data_add %>% group_by(omega_ratiog) %>% mutate(omegar_m = median(omega_ratio))

data_add$omega3_pctstd = scale(data_add$omega3_pct,center=TRUE,scale=TRUE)
data_add$omega6_pctstd = scale(data_add$omega6_pct,center=TRUE,scale=TRUE)
data_add$omegar_std = scale(data_add$omega_ratio,center=TRUE,scale=TRUE)

data_add$wh_ratio = data_add$waist/data_add$hip
data_add = data_add %>% mutate(gas_dises_yn = ifelse(is.na(gas_dises)==FALSE, 1, 0))
data_add = data_add %>% mutate(aspir_yn = ifelse(aspirin==1, 1, 0))
data_add = data_add %>% mutate(meat_yn = ifelse(meat>0, 1, 0))
data_add = data_add %>% mutate(fh_bowel = ifelse(f_history==4 | m_history==4, 1, 0))
data_add = data_add %>% mutate(fh_lung = ifelse(f_history==3 | m_history==3, 1, 0))
data_add = data_add %>% mutate(fh_breast = ifelse(f_history==5 | m_history==5, 1, 0))
data_add = data_add %>% mutate(fh_prostat = ifelse(f_history==13 | m_history==13, 1, 0))

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


##----------------------------------------------------------------------------------------------------------##


# 2. Esophagus
## + gas_dises_yn + wh_ratio

####### continuous
### omega-3
fit2_o3 = summary(coxph(Surv(time_to_cancer_year, status)~ omega3_pctstd + strata(ageg) + strata(sex) + TDI + eth_group_new + BMI + Smoking_status + 
                       Alcohol_status + IPAQ_activity + gas_dises_yn + wh_ratio, data = data2))
contin2_o3 = as.data.frame(fit2_o3$conf.int)[1,-2]

### omega-6
fit2_o6 = summary(coxph(Surv(time_to_cancer_year, status)~ omega6_pctstd + strata(ageg) + strata(sex) + TDI + eth_group_new + BMI + Smoking_status + 
                          Alcohol_status + IPAQ_activity + gas_dises_yn + wh_ratio, data = data2))
contin2_o6 = as.data.frame(fit2_o6$conf.int)[1,-2]

### omega-ratio
fit2_or = summary(coxph(Surv(time_to_cancer_year, status)~ omegar_std + strata(ageg) + strata(sex) + TDI + eth_group_new + BMI + Smoking_status + 
                          Alcohol_status + IPAQ_activity + gas_dises_yn + wh_ratio, data = data2))
contin2_or = as.data.frame(fit2_or$conf.int)[1,-2]

####### quintiles and test for trend

### omega-3
fit2_o3_cat = summary(coxph(Surv(time_to_cancer_year, status)~ omega3_pctg + strata(ageg) + strata(sex) + TDI + eth_group_new + BMI + Smoking_status + 
                              Alcohol_status + IPAQ_activity + gas_dises_yn + wh_ratio, data = data2))
cat2_o3 = as.matrix(fit2_o3_cat$conf.int)[1:4,-2]
cat2_o3_r = c(t(cat2_o3))

fit2_o3_tre = summary(coxph(Surv(time_to_cancer_year, status)~ omega3_m + strata(ageg) + strata(sex) + TDI + eth_group_new + BMI + Smoking_status + 
                              Alcohol_status + IPAQ_activity + gas_dises_yn + wh_ratio, data = data2))
p_val = fit2_o3_tre$coefficients[1,5]

cat2_o3_r = c(cat2_o3_r, p_val)


### omega-6
fit2_o6_cat = summary(coxph(Surv(time_to_cancer_year, status)~ omega6_pctg + strata(ageg) + strata(sex) + TDI + eth_group_new + BMI + Smoking_status + 
                              Alcohol_status + IPAQ_activity + gas_dises_yn + wh_ratio, data = data2))
cat2_o6 = as.matrix(fit2_o6_cat$conf.int)[1:4,-2]
cat2_o6_r = c(t(cat2_o6))

fit2_o6_tre = summary(coxph(Surv(time_to_cancer_year, status)~ omega6_m + strata(ageg) + strata(sex) + TDI + eth_group_new + BMI + Smoking_status + 
                              Alcohol_status + IPAQ_activity + gas_dises_yn + wh_ratio, data = data2))
p_val = fit2_o6_tre$coefficients[1,5]

cat2_o6_r = c(cat2_o6_r, p_val)

### omega-ratio
fit2_or_cat = summary(coxph(Surv(time_to_cancer_year, status)~ omega_ratiog + strata(ageg) + strata(sex) + TDI + eth_group_new + BMI + Smoking_status + 
                              Alcohol_status + IPAQ_activity + gas_dises_yn + wh_ratio, data = data2))
cat2_or = as.matrix(fit2_or_cat$conf.int)[1:4,-2]
cat2_or_r = c(t(cat2_or))
fit2_or_tre = summary(coxph(Surv(time_to_cancer_year, status)~ omegar_m + strata(ageg) + strata(sex) + TDI + eth_group_new + BMI + Smoking_status + 
                              Alcohol_status + IPAQ_activity + gas_dises_yn + wh_ratio, data = data2))
p_val = fit2_or_tre$coefficients[1,5]

cat2_or_r = c(cat2_or_r, p_val)

##----------------------------------------------------------------------------------------------------------##

# 4.	Colon
## + diabet + aspir_yn + meat_yn + wh_ratio + fh_bowel

####### continuous
### omega-3
fit4_o3 = summary(coxph(Surv(time_to_cancer_year, status)~ omega3_pctstd + strata(ageg) + strata(sex) + TDI + eth_group_new + BMI + Smoking_status + 
                          Alcohol_status + IPAQ_activity + diabet + aspir_yn + meat_yn + wh_ratio + fh_bowel, data = data4))
contin4_o3 = as.data.frame(fit4_o3$conf.int)[1,-2]

### omega-6
fit4_o6 = summary(coxph(Surv(time_to_cancer_year, status)~ omega6_pctstd + strata(ageg) + strata(sex) + TDI + eth_group_new + BMI + Smoking_status + 
                          Alcohol_status + IPAQ_activity + diabet + aspir_yn + meat_yn + wh_ratio + fh_bowel, data = data4))
contin4_o6 = as.data.frame(fit4_o6$conf.int)[1,-2]

### omega-ratio
fit4_or = summary(coxph(Surv(time_to_cancer_year, status)~ omegar_std + strata(ageg) + strata(sex) + TDI + eth_group_new + BMI + Smoking_status + 
                          Alcohol_status + IPAQ_activity + diabet + aspir_yn + meat_yn + wh_ratio + fh_bowel, data = data4))
contin4_or = as.data.frame(fit4_or$conf.int)[1,-2]

####### quintiles and test for trend

### omega-3
fit4_o3_cat = summary(coxph(Surv(time_to_cancer_year, status)~ omega3_pctg + strata(ageg) + strata(sex) + TDI + eth_group_new + BMI + Smoking_status + 
                              Alcohol_status + IPAQ_activity + diabet + aspir_yn + meat_yn + wh_ratio + fh_bowel, data = data4))
cat4_o3 = as.matrix(fit4_o3_cat$conf.int)[1:4,-2]
cat4_o3_r = c(t(cat4_o3))

fit4_o3_tre = summary(coxph(Surv(time_to_cancer_year, status)~ omega3_m + strata(ageg) + strata(sex) + TDI + eth_group_new + BMI + Smoking_status + 
                              Alcohol_status + IPAQ_activity + diabet + aspir_yn + meat_yn + wh_ratio + fh_bowel, data = data4))
p_val = fit4_o3_tre$coefficients[1,5]

cat4_o3_r = c(cat4_o3_r, p_val)


### omega-6
fit4_o6_cat = summary(coxph(Surv(time_to_cancer_year, status)~ omega6_pctg + strata(ageg) + strata(sex) + TDI + eth_group_new + BMI + Smoking_status + 
                              Alcohol_status + IPAQ_activity + diabet + aspir_yn + meat_yn + wh_ratio + fh_bowel, data = data4))
cat4_o6 = as.matrix(fit4_o6_cat$conf.int)[1:4,-2]
cat4_o6_r = c(t(cat4_o6))

fit4_o6_tre = summary(coxph(Surv(time_to_cancer_year, status)~ omega6_m + strata(ageg) + strata(sex) + TDI + eth_group_new + BMI + Smoking_status + 
                              Alcohol_status + IPAQ_activity + diabet + aspir_yn + meat_yn + wh_ratio + fh_bowel, data = data4))
p_val = fit4_o6_tre$coefficients[1,5]

cat4_o6_r = c(cat4_o6_r, p_val)

### omega-ratio
fit4_or_cat = summary(coxph(Surv(time_to_cancer_year, status)~ omega_ratiog + strata(ageg) + strata(sex) + TDI + eth_group_new + BMI + Smoking_status + 
                              Alcohol_status + IPAQ_activity + diabet + aspir_yn + meat_yn + wh_ratio + fh_bowel, data = data4))
cat4_or = as.matrix(fit4_or_cat$conf.int)[1:4,-2]
cat4_or_r = c(t(cat4_or))

fit4_or_tre = summary(coxph(Surv(time_to_cancer_year, status)~ omegar_m + strata(ageg) + strata(sex) + TDI + eth_group_new + BMI + Smoking_status + 
                              Alcohol_status + IPAQ_activity + diabet + aspir_yn + meat_yn + wh_ratio + fh_bowel, data = data4))
p_val = fit4_or_tre$coefficients[1,5]

cat4_or_r = c(cat4_or_r, p_val)


##----------------------------------------------------------------------------------------------------------##

# 5.	Rectum
## + diabet + aspir_yn + meat_yn + wh_ratio + fh_bowel

####### continuous
### omega-3
fit5_o3 = summary(coxph(Surv(time_to_cancer_year, status)~ omega3_pctstd + strata(ageg) + strata(sex) + TDI + eth_group_new + BMI + Smoking_status + 
                          Alcohol_status + IPAQ_activity + diabet + aspir_yn + meat_yn + wh_ratio + fh_bowel, data = data5))
contin5_o3 = as.data.frame(fit5_o3$conf.int)[1,-2]

### omega-6
fit5_o6 = summary(coxph(Surv(time_to_cancer_year, status)~ omega6_pctstd + strata(ageg) + strata(sex) + TDI + eth_group_new + BMI + Smoking_status + 
                          Alcohol_status + IPAQ_activity + diabet + aspir_yn + meat_yn + wh_ratio + fh_bowel, data = data5))
contin5_o6 = as.data.frame(fit5_o6$conf.int)[1,-2]

### omega-ratio
fit5_or = summary(coxph(Surv(time_to_cancer_year, status)~ omegar_std + strata(ageg) + strata(sex) + TDI + eth_group_new + BMI + Smoking_status + 
                          Alcohol_status + IPAQ_activity + diabet + aspir_yn + meat_yn + wh_ratio + fh_bowel, data = data5))
contin5_or = as.data.frame(fit5_or$conf.int)[1,-2]

####### quintiles and test for trend

### omega-3
fit5_o3_cat = summary(coxph(Surv(time_to_cancer_year, status)~ omega3_pctg + strata(ageg) + strata(sex) + TDI + eth_group_new + BMI + Smoking_status + 
                              Alcohol_status + IPAQ_activity + diabet + aspir_yn + meat_yn + wh_ratio + fh_bowel, data = data5))
cat5_o3 = as.matrix(fit5_o3_cat$conf.int)[1:4,-2]
cat5_o3_r = c(t(cat5_o3))

fit5_o3_tre = summary(coxph(Surv(time_to_cancer_year, status)~ omega3_m + strata(ageg) + strata(sex) + TDI + eth_group_new + BMI + Smoking_status + 
                              Alcohol_status + IPAQ_activity + diabet + aspir_yn + meat_yn + wh_ratio + fh_bowel, data = data5))
p_val = fit5_o3_tre$coefficients[1,5]

cat5_o3_r = c(cat5_o3_r, p_val)


### omega-6
fit5_o6_cat = summary(coxph(Surv(time_to_cancer_year, status)~ omega6_pctg + strata(ageg) + strata(sex) + TDI + eth_group_new + BMI + Smoking_status + 
                              Alcohol_status + IPAQ_activity + diabet + aspir_yn + meat_yn + wh_ratio + fh_bowel, data = data5))
cat5_o6 = as.matrix(fit5_o6_cat$conf.int)[1:4,-2]
cat5_o6_r = c(t(cat5_o6))

fit5_o6_tre = summary(coxph(Surv(time_to_cancer_year, status)~ omega6_m + strata(ageg) + strata(sex) + TDI + eth_group_new + BMI + Smoking_status + 
                              Alcohol_status + IPAQ_activity + diabet + aspir_yn + meat_yn + wh_ratio + fh_bowel, data = data5))
p_val = fit5_o6_tre$coefficients[1,5]

cat5_o6_r = c(cat5_o6_r, p_val)

### omega-ratio
fit5_or_cat = summary(coxph(Surv(time_to_cancer_year, status)~ omega_ratiog + strata(ageg) + strata(sex) + TDI + eth_group_new + BMI + Smoking_status + 
                              Alcohol_status + IPAQ_activity + diabet + aspir_yn + meat_yn + wh_ratio + fh_bowel, data = data5))
cat5_or = as.matrix(fit5_or_cat$conf.int)[1:4,-2]
cat5_or_r = c(t(cat5_or))

fit5_or_tre = summary(coxph(Surv(time_to_cancer_year, status)~ omegar_m + strata(ageg) + strata(sex) + TDI + eth_group_new + BMI + Smoking_status + 
                              Alcohol_status + IPAQ_activity + diabet + aspir_yn + meat_yn + wh_ratio + fh_bowel, data = data5))
p_val = fit5_or_tre$coefficients[1,5]

cat5_or_r = c(cat5_or_r, p_val)

##----------------------------------------------------------------------------------------------------------##

# 7.	Pancreas
## + diabet

####### continuous
### omega-3
fit7_o3 = summary(coxph(Surv(time_to_cancer_year, status)~ omega3_pctstd + strata(ageg) + strata(sex) + TDI + eth_group_new + BMI + Smoking_status + 
                          Alcohol_status + IPAQ_activity + diabet, data = data7))
contin7_o3 = as.data.frame(fit7_o3$conf.int)[1,-2]

### omega-6
fit7_o6 = summary(coxph(Surv(time_to_cancer_year, status)~ omega6_pctstd + strata(ageg) + strata(sex) + TDI + eth_group_new + BMI + Smoking_status + 
                          Alcohol_status + IPAQ_activity + diabet, data = data7))
contin7_o6 = as.data.frame(fit7_o6$conf.int)[1,-2]

### omega-ratio
fit7_or = summary(coxph(Surv(time_to_cancer_year, status)~ omegar_std + strata(ageg) + strata(sex) + TDI + eth_group_new + BMI + Smoking_status + 
                          Alcohol_status + IPAQ_activity + diabet, data = data7))
contin7_or = as.data.frame(fit7_or$conf.int)[1,-2]

####### quintiles and test for trend

### omega-3
fit7_o3_cat = summary(coxph(Surv(time_to_cancer_year, status)~ omega3_pctg + strata(ageg) + strata(sex) + TDI + eth_group_new + BMI + Smoking_status + 
                              Alcohol_status + IPAQ_activity + diabet, data = data7))
cat7_o3 = as.matrix(fit7_o3_cat$conf.int)[1:4,-2]
cat7_o3_r = c(t(cat7_o3))

fit7_o3_tre = summary(coxph(Surv(time_to_cancer_year, status)~ omega3_m + strata(ageg) + strata(sex) + TDI + eth_group_new + BMI + Smoking_status + 
                              Alcohol_status + IPAQ_activity + diabet, data = data7))
p_val = fit7_o3_tre$coefficients[1,5]

cat7_o3_r = c(cat7_o3_r, p_val)


### omega-6
fit7_o6_cat = summary(coxph(Surv(time_to_cancer_year, status)~ omega6_pctg + strata(ageg) + strata(sex) + TDI + eth_group_new + BMI + Smoking_status + 
                              Alcohol_status + IPAQ_activity + diabet, data = data7))
cat7_o6 = as.matrix(fit7_o6_cat$conf.int)[1:4,-2]
cat7_o6_r = c(t(cat7_o6))

fit7_o6_tre = summary(coxph(Surv(time_to_cancer_year, status)~ omega6_m + strata(ageg) + strata(sex) + TDI + eth_group_new + BMI + Smoking_status + 
                              Alcohol_status + IPAQ_activity + diabet, data = data7))
p_val = fit7_o6_tre$coefficients[1,5]

cat7_o6_r = c(cat7_o6_r, p_val)

### omega-ratio
fit7_or_cat = summary(coxph(Surv(time_to_cancer_year, status)~ omega_ratiog + strata(ageg) + strata(sex) + TDI + eth_group_new + BMI + Smoking_status + 
                              Alcohol_status + IPAQ_activity + diabet, data = data7))
cat7_or = as.matrix(fit7_or_cat$conf.int)[1:4,-2]
cat7_or_r = c(t(cat7_or))

fit7_or_tre = summary(coxph(Surv(time_to_cancer_year, status)~ omegar_m + strata(ageg) + strata(sex) + TDI + eth_group_new + BMI + Smoking_status + 
                              Alcohol_status + IPAQ_activity + diabet, data = data7))
p_val = fit7_or_tre$coefficients[1,5]

cat7_or_r = c(cat7_or_r, p_val)


##----------------------------------------------------------------------------------------------------------##

# 8.	Lung
## + fh_lung

####### continuous
### omega-3
fit8_o3 = summary(coxph(Surv(time_to_cancer_year, status)~ omega3_pctstd + strata(ageg) + strata(sex) + TDI + eth_group_new + BMI + Smoking_status + 
                          Alcohol_status + IPAQ_activity + fh_lung, data = data8))
contin8_o3 = as.data.frame(fit8_o3$conf.int)[1,-2]

### omega-6
fit8_o6 = summary(coxph(Surv(time_to_cancer_year, status)~ omega6_pctstd + strata(ageg) + strata(sex) + TDI + eth_group_new + BMI + Smoking_status + 
                          Alcohol_status + IPAQ_activity + fh_lung, data = data8))
contin8_o6 = as.data.frame(fit8_o6$conf.int)[1,-2]

### omega-ratio
fit8_or = summary(coxph(Surv(time_to_cancer_year, status)~ omegar_std + strata(ageg) + strata(sex) + TDI + eth_group_new + BMI + Smoking_status + 
                          Alcohol_status + IPAQ_activity + fh_lung, data = data8))
contin8_or = as.data.frame(fit8_or$conf.int)[1,-2]

####### quintiles and test for trend

### omega-3
fit8_o3_cat = summary(coxph(Surv(time_to_cancer_year, status)~ omega3_pctg + strata(ageg) + strata(sex) + TDI + eth_group_new + BMI + Smoking_status + 
                              Alcohol_status + IPAQ_activity + fh_lung, data = data8))
cat8_o3 = as.matrix(fit8_o3_cat$conf.int)[1:4,-2]
cat8_o3_r = c(t(cat8_o3))

fit8_o3_tre = summary(coxph(Surv(time_to_cancer_year, status)~ omega3_m + strata(ageg) + strata(sex) + TDI + eth_group_new + BMI + Smoking_status + 
                              Alcohol_status + IPAQ_activity + fh_lung, data = data8))
p_val = fit8_o3_tre$coefficients[1,5]

cat8_o3_r = c(cat8_o3_r, p_val)


### omega-6
fit8_o6_cat = summary(coxph(Surv(time_to_cancer_year, status)~ omega6_pctg + strata(ageg) + strata(sex) + TDI + eth_group_new + BMI + Smoking_status + 
                              Alcohol_status + IPAQ_activity + fh_lung, data = data8))
cat8_o6 = as.matrix(fit8_o6_cat$conf.int)[1:4,-2]
cat8_o6_r = c(t(cat8_o6))

fit8_o6_tre = summary(coxph(Surv(time_to_cancer_year, status)~ omega6_m + strata(ageg) + strata(sex) + TDI + eth_group_new + BMI + Smoking_status + 
                              Alcohol_status + IPAQ_activity + fh_lung, data = data8))
p_val = fit8_o6_tre$coefficients[1,5]

cat8_o6_r = c(cat8_o6_r, p_val)

### omega-ratio
fit8_or_cat = summary(coxph(Surv(time_to_cancer_year, status)~ omega_ratiog + strata(ageg) + strata(sex) + TDI + eth_group_new + BMI + Smoking_status + 
                              Alcohol_status + IPAQ_activity + fh_lung, data = data8))
cat8_or = as.matrix(fit8_or_cat$conf.int)[1:4,-2]
cat8_or_r = c(t(cat8_or))

fit8_or_tre = summary(coxph(Surv(time_to_cancer_year, status)~ omegar_m + strata(ageg) + strata(sex) + TDI + eth_group_new + BMI + Smoking_status + 
                              Alcohol_status + IPAQ_activity + fh_lung, data = data8))
p_val = fit8_or_tre$coefficients[1,5]

cat8_or_r = c(cat8_or_r, p_val)


##----------------------------------------------------------------------------------------------------------##

# 9.	Malignant melanoma
## + skin_color + bright_sun + uv_protect + sunburn + sunlamp

####### continuous
### omega-3
fit9_o3 = summary(coxph(Surv(time_to_cancer_year, status)~ omega3_pctstd + strata(ageg) + strata(sex) + TDI + eth_group_new + BMI + Smoking_status + 
                          Alcohol_status + IPAQ_activity + skin_color + bright_sun + uv_protect + sunburn + sunlamp, data = data9))
contin9_o3 = as.data.frame(fit9_o3$conf.int)[1,-2]

### omega-6
fit9_o6 = summary(coxph(Surv(time_to_cancer_year, status)~ omega6_pctstd + strata(ageg) + strata(sex) + TDI + eth_group_new + BMI + Smoking_status + 
                          Alcohol_status + IPAQ_activity + skin_color + bright_sun + uv_protect + sunburn + sunlamp, data = data9))
contin9_o6 = as.data.frame(fit9_o6$conf.int)[1,-2]

### omega-ratio
fit9_or = summary(coxph(Surv(time_to_cancer_year, status)~ omegar_std + strata(ageg) + strata(sex) + TDI + eth_group_new + BMI + Smoking_status + 
                          Alcohol_status + IPAQ_activity + skin_color + bright_sun + uv_protect + sunburn + sunlamp, data = data9))
contin9_or = as.data.frame(fit9_or$conf.int)[1,-2]

####### quintiles and test for trend

### omega-3
fit9_o3_cat = summary(coxph(Surv(time_to_cancer_year, status)~ omega3_pctg + strata(ageg) + strata(sex) + TDI + eth_group_new + BMI + Smoking_status + 
                              Alcohol_status + IPAQ_activity + skin_color + bright_sun + uv_protect + sunburn + sunlamp, data = data9))
cat9_o3 = as.matrix(fit9_o3_cat$conf.int)[1:4,-2]
cat9_o3_r = c(t(cat9_o3))

fit9_o3_tre = summary(coxph(Surv(time_to_cancer_year, status)~ omega3_m + strata(ageg) + strata(sex) + TDI + eth_group_new + BMI + Smoking_status + 
                              Alcohol_status + IPAQ_activity + skin_color + bright_sun + uv_protect + sunburn + sunlamp, data = data9))
p_val = fit9_o3_tre$coefficients[1,5]

cat9_o3_r = c(cat9_o3_r, p_val)


### omega-6
fit9_o6_cat = summary(coxph(Surv(time_to_cancer_year, status)~ omega6_pctg + strata(ageg) + strata(sex) + TDI + eth_group_new + BMI + Smoking_status + 
                              Alcohol_status + IPAQ_activity + skin_color + bright_sun + uv_protect + sunburn + sunlamp, data = data9))
cat9_o6 = as.matrix(fit9_o6_cat$conf.int)[1:4,-2]
cat9_o6_r = c(t(cat9_o6))

fit9_o6_tre = summary(coxph(Surv(time_to_cancer_year, status)~ omega6_m + strata(ageg) + strata(sex) + TDI + eth_group_new + BMI + Smoking_status + 
                              Alcohol_status + IPAQ_activity + skin_color + bright_sun + uv_protect + sunburn + sunlamp, data = data9))
p_val = fit9_o6_tre$coefficients[1,5]

cat9_o6_r = c(cat9_o6_r, p_val)

### omega-ratio
fit9_or_cat = summary(coxph(Surv(time_to_cancer_year, status)~ omega_ratiog + strata(ageg) + strata(sex) + TDI + eth_group_new + BMI + Smoking_status + 
                              Alcohol_status + IPAQ_activity + skin_color + bright_sun + uv_protect + sunburn + sunlamp, data = data9))
cat9_or = as.matrix(fit9_or_cat$conf.int)[1:4,-2]
cat9_or_r = c(t(cat9_or))

fit9_or_tre = summary(coxph(Surv(time_to_cancer_year, status)~ omegar_m + strata(ageg) + strata(sex) + TDI + eth_group_new + BMI + Smoking_status + 
                              Alcohol_status + IPAQ_activity + skin_color + bright_sun + uv_protect + sunburn + sunlamp, data = data9))
p_val = fit9_or_tre$coefficients[1,5]

cat9_or_r = c(cat9_or_r, p_val)


##----------------------------------------------------------------------------------------------------------##

# 11.	Breast
## + HRT + OC + num_birth + age_mena + meno + hyster + fh_breast

####### continuous
### omega-3
fit11_o3 = summary(coxph(Surv(time_to_cancer_year, status)~ omega3_pctstd + strata(ageg) + strata(sex) + TDI + eth_group_new + BMI + Smoking_status + 
                          Alcohol_status + IPAQ_activity + HRT + OC + num_birth + age_mena + meno + hyster + fh_breast, data = data11))
contin11_o3 = as.data.frame(fit11_o3$conf.int)[1,-2]

### omega-6
fit11_o6 = summary(coxph(Surv(time_to_cancer_year, status)~ omega6_pctstd + strata(ageg) + strata(sex) + TDI + eth_group_new + BMI + Smoking_status + 
                          Alcohol_status + IPAQ_activity + HRT + OC + num_birth + age_mena + meno + hyster + fh_breast, data = data11))
contin11_o6 = as.data.frame(fit11_o6$conf.int)[1,-2]

### omega-ratio
fit11_or = summary(coxph(Surv(time_to_cancer_year, status)~ omegar_std + strata(ageg) + strata(sex) + TDI + eth_group_new + BMI + Smoking_status + 
                          Alcohol_status + IPAQ_activity + HRT + OC + num_birth + age_mena + meno + hyster + fh_breast, data = data11))
contin11_or = as.data.frame(fit11_or$conf.int)[1,-2]

####### quintiles and test for trend

### omega-3
fit11_o3_cat = summary(coxph(Surv(time_to_cancer_year, status)~ omega3_pctg + strata(ageg) + strata(sex) + TDI + eth_group_new + BMI + Smoking_status + 
                              Alcohol_status + IPAQ_activity + HRT + OC + num_birth + age_mena + meno + hyster + fh_breast, data = data11))
cat11_o3 = as.matrix(fit11_o3_cat$conf.int)[1:4,-2]
cat11_o3_r = c(t(cat11_o3))

fit11_o3_tre = summary(coxph(Surv(time_to_cancer_year, status)~ omega3_m + strata(ageg) + strata(sex) + TDI + eth_group_new + BMI + Smoking_status + 
                              Alcohol_status + IPAQ_activity + HRT + OC + num_birth + age_mena + meno + hyster + fh_breast, data = data11))
p_val = fit11_o3_tre$coefficients[1,5]

cat11_o3_r = c(cat11_o3_r, p_val)


### omega-6
fit11_o6_cat = summary(coxph(Surv(time_to_cancer_year, status)~ omega6_pctg + strata(ageg) + strata(sex) + TDI + eth_group_new + BMI + Smoking_status + 
                              Alcohol_status + IPAQ_activity + HRT + OC + num_birth + age_mena + meno + hyster + fh_breast, data = data11))
cat11_o6 = as.matrix(fit11_o6_cat$conf.int)[1:4,-2]
cat11_o6_r = c(t(cat11_o6))

fit11_o6_tre = summary(coxph(Surv(time_to_cancer_year, status)~ omega6_m + strata(ageg) + strata(sex) + TDI + eth_group_new + BMI + Smoking_status + 
                              Alcohol_status + IPAQ_activity + HRT + OC + num_birth + age_mena + meno + hyster + fh_breast, data = data11))
p_val = fit11_o6_tre$coefficients[1,5]

cat11_o6_r = c(cat11_o6_r, p_val)

### omega-ratio
fit11_or_cat = summary(coxph(Surv(time_to_cancer_year, status)~ omega_ratiog + strata(ageg) + strata(sex) + TDI + eth_group_new + BMI + Smoking_status + 
                              Alcohol_status + IPAQ_activity + HRT + OC + num_birth + age_mena + meno + hyster + fh_breast, data = data11))
cat11_or = as.matrix(fit11_or_cat$conf.int)[1:4,-2]
cat11_or_r = c(t(cat11_or))

fit11_or_tre = summary(coxph(Surv(time_to_cancer_year, status)~ omegar_m + strata(ageg) + strata(sex) + TDI + eth_group_new + BMI + Smoking_status + 
                              Alcohol_status + IPAQ_activity + HRT + OC + num_birth + age_mena + meno + hyster + fh_breast, data = data11))
p_val = fit11_or_tre$coefficients[1,5]

cat11_or_r = c(cat11_or_r, p_val)


##----------------------------------------------------------------------------------------------------------##

# 12.	Uterus
## + HRT + OC + num_birth + age_mena + meno + hyster

####### continuous
### omega-3
fit12_o3 = summary(coxph(Surv(time_to_cancer_year, status)~ omega3_pctstd + strata(ageg) + strata(sex) + TDI + eth_group_new + BMI + Smoking_status + 
                           Alcohol_status + IPAQ_activity + HRT + OC + num_birth + age_mena + meno + hyster, data = data12))
contin12_o3 = as.data.frame(fit12_o3$conf.int)[1,-2]

### omega-6
fit12_o6 = summary(coxph(Surv(time_to_cancer_year, status)~ omega6_pctstd + strata(ageg) + strata(sex) + TDI + eth_group_new + BMI + Smoking_status + 
                           Alcohol_status + IPAQ_activity + HRT + OC + num_birth + age_mena + meno + hyster , data = data12))
contin12_o6 = as.data.frame(fit12_o6$conf.int)[1,-2]

### omega-ratio
fit12_or = summary(coxph(Surv(time_to_cancer_year, status)~ omegar_std + strata(ageg) + strata(sex) + TDI + eth_group_new + BMI + Smoking_status + 
                           Alcohol_status + IPAQ_activity + HRT + OC + num_birth + age_mena + meno + hyster, data = data12))
contin12_or = as.data.frame(fit12_or$conf.int)[1,-2]

####### quintiles and test for trend

### omega-3
fit12_o3_cat = summary(coxph(Surv(time_to_cancer_year, status)~ omega3_pctg + strata(ageg) + strata(sex) + TDI + eth_group_new + BMI + Smoking_status + 
                               Alcohol_status + IPAQ_activity + HRT + OC + num_birth + age_mena + meno + hyster, data = data12))
cat12_o3 = as.matrix(fit12_o3_cat$conf.int)[1:4,-2]
cat12_o3_r = c(t(cat12_o3))

fit12_o3_tre = summary(coxph(Surv(time_to_cancer_year, status)~ omega3_m + strata(ageg) + strata(sex) + TDI + eth_group_new + BMI + Smoking_status + 
                               Alcohol_status + IPAQ_activity + HRT + OC + num_birth + age_mena + meno + hyster, data = data12))
p_val = fit12_o3_tre$coefficients[1,5]

cat12_o3_r = c(cat12_o3_r, p_val)


### omega-6
fit12_o6_cat = summary(coxph(Surv(time_to_cancer_year, status)~ omega6_pctg + strata(ageg) + strata(sex) + TDI + eth_group_new + BMI + Smoking_status + 
                               Alcohol_status + IPAQ_activity + HRT + OC + num_birth + age_mena + meno + hyster, data = data12))
cat12_o6 = as.matrix(fit12_o6_cat$conf.int)[1:4,-2]
cat12_o6_r = c(t(cat12_o6))

fit12_o6_tre = summary(coxph(Surv(time_to_cancer_year, status)~ omega6_m + strata(ageg) + strata(sex) + TDI + eth_group_new + BMI + Smoking_status + 
                               Alcohol_status + IPAQ_activity + HRT + OC + num_birth + age_mena + meno + hyster, data = data12))
p_val = fit12_o4_tre$coefficients[1,5]

cat12_o6_r = c(cat12_o6_r, p_val)

### omega-ratio
fit12_or_cat = summary(coxph(Surv(time_to_cancer_year, status)~ omega_ratiog + strata(ageg) + strata(sex) + TDI + eth_group_new + BMI + Smoking_status + 
                               Alcohol_status + IPAQ_activity + HRT + OC + num_birth + age_mena + meno + hyster, data = data12))
cat12_or = as.matrix(fit12_or_cat$conf.int)[1:4,-2]
cat12_or_r = c(t(cat12_or))

fit12_or_tre = summary(coxph(Surv(time_to_cancer_year, status)~ omegar_m + strata(ageg) + strata(sex) + TDI + eth_group_new + BMI + Smoking_status + 
                               Alcohol_status + IPAQ_activity + HRT + OC + num_birth + age_mena + meno + hyster, data = data12))
p_val = fit12_or_tre$coefficients[1,5]

cat12_or_r = c(cat12_or_r, p_val)



##----------------------------------------------------------------------------------------------------------##

# 13.	Ovary
## + HRT + OC + num_birth + age_mena + meno + hyster

####### continuous
### omega-3
fit13_o3 = summary(coxph(Surv(time_to_cancer_year, status)~ omega3_pctstd + strata(ageg) + strata(sex) + TDI + eth_group_new + BMI + Smoking_status + 
                           Alcohol_status + IPAQ_activity + HRT + OC + num_birth + age_mena + meno + hyster, data = data13))
contin13_o3 = as.data.frame(fit13_o3$conf.int)[1,-2]

### omega-6
fit13_o6 = summary(coxph(Surv(time_to_cancer_year, status)~ omega6_pctstd + strata(ageg) + strata(sex) + TDI + eth_group_new + BMI + Smoking_status + 
                           Alcohol_status + IPAQ_activity + HRT + OC + num_birth + age_mena + meno + hyster , data = data13))
contin13_o6 = as.data.frame(fit13_o6$conf.int)[1,-2]

### omega-ratio
fit13_or = summary(coxph(Surv(time_to_cancer_year, status)~ omegar_std + strata(ageg) + strata(sex) + TDI + eth_group_new + BMI + Smoking_status + 
                           Alcohol_status + IPAQ_activity + HRT + OC + num_birth + age_mena + meno + hyster, data = data13))
contin13_or = as.data.frame(fit13_or$conf.int)[1,-2]

####### quintiles and test for trend

### omega-3
fit13_o3_cat = summary(coxph(Surv(time_to_cancer_year, status)~ omega3_pctg + strata(ageg) + strata(sex) + TDI + eth_group_new + BMI + Smoking_status + 
                               Alcohol_status + IPAQ_activity + HRT + OC + num_birth + age_mena + meno + hyster, data = data13))
cat13_o3 = as.matrix(fit13_o3_cat$conf.int)[1:4,-2]
cat13_o3_r = c(t(cat13_o3))

fit13_o3_tre = summary(coxph(Surv(time_to_cancer_year, status)~ omega3_m + strata(ageg) + strata(sex) + TDI + eth_group_new + BMI + Smoking_status + 
                               Alcohol_status + IPAQ_activity + HRT + OC + num_birth + age_mena + meno + hyster, data = data13))
p_val = fit13_o3_tre$coefficients[1,5]

cat13_o3_r = c(cat13_o3_r, p_val)


### omega-6
fit13_o6_cat = summary(coxph(Surv(time_to_cancer_year, status)~ omega6_pctg + strata(ageg) + strata(sex) + TDI + eth_group_new + BMI + Smoking_status + 
                               Alcohol_status + IPAQ_activity + HRT + OC + num_birth + age_mena + meno + hyster, data = data13))
cat13_o6 = as.matrix(fit13_o6_cat$conf.int)[1:4,-2]
cat13_o6_r = c(t(cat13_o6))

fit13_o6_tre = summary(coxph(Surv(time_to_cancer_year, status)~ omega6_m + strata(ageg) + strata(sex) + TDI + eth_group_new + BMI + Smoking_status + 
                               Alcohol_status + IPAQ_activity + HRT + OC + num_birth + age_mena + meno + hyster, data = data13))
p_val = fit13_o6_tre$coefficients[1,5]

cat13_o6_r = c(cat13_o6_r, p_val)

### omega-ratio
fit13_or_cat = summary(coxph(Surv(time_to_cancer_year, status)~ omega_ratiog + strata(ageg) + strata(sex) + TDI + eth_group_new + BMI + Smoking_status + 
                               Alcohol_status + IPAQ_activity + HRT + OC + num_birth + age_mena + meno + hyster, data = data13))
cat13_or = as.matrix(fit13_or_cat$conf.int)[1:4,-2]
cat13_or_r = c(t(cat13_or))

fit13_or_tre = summary(coxph(Surv(time_to_cancer_year, status)~ omegar_m + strata(ageg) + strata(sex) + TDI + eth_group_new + BMI + Smoking_status + 
                               Alcohol_status + IPAQ_activity + HRT + OC + num_birth + age_mena + meno + hyster, data = data13))
p_val = fit13_or_tre$coefficients[1,5]

cat13_or_r = c(cat13_or_r, p_val)



##----------------------------------------------------------------------------------------------------------##

# 14.	Prostate
## + fh_prostat

####### continuous
### omega-3
fit14_o3 = summary(coxph(Surv(time_to_cancer_year, status)~ omega3_pctstd + strata(ageg) + strata(sex) + TDI + eth_group_new + BMI + Smoking_status + 
                           Alcohol_status + IPAQ_activity + fh_prostat, data = data14))
contin14_o3 = as.data.frame(fit14_o3$conf.int)[1,-2]

### omega-6
fit14_o6 = summary(coxph(Surv(time_to_cancer_year, status)~ omega6_pctstd + strata(ageg) + strata(sex) + TDI + eth_group_new + BMI + Smoking_status + 
                           Alcohol_status + IPAQ_activity + fh_prostat, data = data14))
contin14_o6 = as.data.frame(fit14_o6$conf.int)[1,-2]

### omega-ratio
fit14_or = summary(coxph(Surv(time_to_cancer_year, status)~ omegar_std + strata(ageg) + strata(sex) + TDI + eth_group_new + BMI + Smoking_status + 
                           Alcohol_status + IPAQ_activity + fh_prostat, data = data14))
contin14_or = as.data.frame(fit14_or$conf.int)[1,-2]

####### quintiles and test for trend

### omega-3
fit14_o3_cat = summary(coxph(Surv(time_to_cancer_year, status)~ omega3_pctg + strata(ageg) + strata(sex) + TDI + eth_group_new + BMI + Smoking_status + 
                               Alcohol_status + IPAQ_activity + fh_prostat, data = data14))
cat14_o3 = as.matrix(fit14_o3_cat$conf.int)[1:4,-2]
cat14_o3_r = c(t(cat14_o3))

fit14_o3_tre = summary(coxph(Surv(time_to_cancer_year, status)~ omega3_m + strata(ageg) + strata(sex) + TDI + eth_group_new + BMI + Smoking_status + 
                               Alcohol_status + IPAQ_activity + fh_prostat, data = data14))
p_val = fit14_o3_tre$coefficients[1,5]

cat14_o3_r = c(cat14_o3_r, p_val)


### omega-6
fit14_o6_cat = summary(coxph(Surv(time_to_cancer_year, status)~ omega6_pctg + strata(ageg) + strata(sex) + TDI + eth_group_new + BMI + Smoking_status + 
                               Alcohol_status + IPAQ_activity + fh_prostat, data = data14))
cat14_o6 = as.matrix(fit14_o6_cat$conf.int)[1:4,-2]
cat14_o6_r = c(t(cat14_o6))

fit14_o6_tre = summary(coxph(Surv(time_to_cancer_year, status)~ omega6_m + strata(ageg) + strata(sex) + TDI + eth_group_new + BMI + Smoking_status + 
                               Alcohol_status + IPAQ_activity + fh_prostat, data = data14))
p_val = fit14_o6_tre$coefficients[1,5]

cat14_o6_r = c(cat14_o6_r, p_val)

### omega-ratio
fit14_or_cat = summary(coxph(Surv(time_to_cancer_year, status)~ omega_ratiog + strata(ageg) + strata(sex) + TDI + eth_group_new + BMI + Smoking_status + 
                               Alcohol_status + IPAQ_activity + fh_prostat, data = data14))
cat14_or = as.matrix(fit14_or_cat$conf.int)[1:4,-2]
cat14_or_r = c(t(cat14_or))

fit14_or_tre = summary(coxph(Surv(time_to_cancer_year, status)~ omegar_m + strata(ageg) + strata(sex) + TDI + eth_group_new + BMI + Smoking_status + 
                               Alcohol_status + IPAQ_activity + fh_prostat, data = data14))
p_val = fit14_or_tre$coefficients[1,5]

cat14_or_r = c(cat14_or_r, p_val)



##----------------------------------------------------------------------------------------------------------##


## Results summary

con_o3_sum = rbind(contin2_o3, contin4_o3, contin5_o3, contin7_o3, contin8_o3,
                   contin9_o3, contin11_o3, contin12_o3, contin13_o3, contin14_o3)

con_o6_sum = rbind(contin2_o6, contin4_o6, contin5_o6, contin7_o6, contin8_o6,
                   contin9_o6, contin11_o6, contin12_o6, contin13_o6, contin14_o6)

con_or_sum = rbind(contin2_or, contin4_or, contin5_or, contin7_or, contin8_or,
                   contin9_or, contin11_or, contin12_or, contin13_or, contin14_or)
  
cat_o3_sum = rbind(cat2_o3_r, cat4_o3_r, cat5_o3_r, cat7_o3_r, cat8_o3_r,
                   cat9_o3_r, cat11_o3_r, cat12_o3_r, cat13_o3_r, cat14_o3_r)

cat_o6_sum = rbind(cat2_o6_r, cat4_o6_r, cat5_o6_r, cat7_o6_r, cat8_o6_r,
                   cat9_o6_r, cat11_o6_r, cat12_o6_r, cat13_o6_r, cat14_o6_r)

cat_or_sum = rbind(cat2_or_r, cat4_or_r, cat5_or_r, cat7_or_r, cat8_or_r,
                   cat9_or_r, cat11_or_r, cat12_or_r, cat13_or_r, cat14_or_r)

write.csv(con_o3_sum, file="/Users/yuchen/Desktop/Cancer/Data/con_o3_sum.csv",row.names = F)
write.csv(con_o6_sum, file="/Users/yuchen/Desktop/Cancer/Data/con_o6_sum.csv",row.names = F)
write.csv(con_or_sum, file="/Users/yuchen/Desktop/Cancer/Data/con_or_sum.csv",row.names = F)
write.csv(cat_o3_sum, file="/Users/yuchen/Desktop/Cancer/Data/cat_o3_sum.csv",row.names = F)
write.csv(cat_o6_sum, file="/Users/yuchen/Desktop/Cancer/Data/cat_o6_sum.csv",row.names = F)
write.csv(cat_or_sum, file="/Users/yuchen/Desktop/Cancer/Data/cat_or_sum.csv",row.names = F)




