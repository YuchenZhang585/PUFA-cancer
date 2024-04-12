library(smoothHR)
library(survival)
library(rms)
library(Hmisc) 
library(dplyr)
library(data.table)

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

new_mydata= data_add %>% select(ID, omega_ratio2, omega3_pct, omega6_pct, ageg, TDI, eth_group_new, sex, BMI, Smoking_status, Alcohol_status, IPAQ_activity, 
                                time_to_cancer_year, status, cancer_sub) 
new_mydata= new_mydata %>% filter(omega_ratio2 <= 40)


data_prostate = new_mydata %>% filter(cancer_sub == "14_prostate" | is.na(cancer_sub) == TRUE) %>% select(-cancer_sub)
data_breast = new_mydata %>% filter(cancer_sub == "11_breast" | is.na(cancer_sub) == TRUE)  %>% select(-cancer_sub)
data_head = new_mydata %>% filter(cancer_sub == "01_head" | is.na(cancer_sub) == TRUE)  %>% select(-cancer_sub)

new_mydata= new_mydata %>% select(-cancer_sub)

# Determine the optimal number of nodes knots,

for (i in 3:7) {
    fit <- cph(Surv(time_to_cancer_year, status) ~ rcs(omega6_pct,i) + strat(ageg) + strat(sex) + TDI + 
                   eth_group_new + BMI + Smoking_status + Alcohol_status + IPAQ_activity, data = new_mydata, x=TRUE)
    tmp <- extractAIC(fit)
    if(i == 3) {AIC = tmp[2]; nk = 3}
    if(tmp[2] < AIC) {AIC = tmp[2]; nk = i} 
}

# The above code is based on AIC Criteria filter minimum AIC Corresponding knots. The results show that the optimal knots by 3.

################ Overall cancer incidence - ratio  #####################

# Start fitting the model

fit = cph(Surv(time_to_cancer_year, status) ~ rcs(omega_ratio2,3) + strat(ageg) + strat(sex) + TDI + 
              eth_group_new + BMI + Smoking_status + Alcohol_status + IPAQ_activity, data = new_mydata, x=TRUE)
anova(fit)

# omega_ratio2     1.77      2   0.4125
# Nonlinear       1.76      1   0.1848

## p-non-linear
p <-round(anova(fit)[,3],3)

# define the reference number
refvalue <- 4

# Package the data , And specify the reference value .
ddist <<- datadist(new_mydata)
ddist$limits$omega_ratio2[2]<-refvalue
options(datadist="ddist")

pred_HR<-Predict(fit,omega_ratio2,ref.zero=TRUE,fun=exp)


tiff(file.path("/Users/yuchen/Desktop/all_ratio.tiff"), units="in", width=4.5, height=4.5, res=300)

#  Sets the background color of the density curve 
grey <- "#D3D3D3"
    #  Draw the left and right double axes baseplot
par(mar = c(5, 4, 4, 4) + 0.3)
par(xpd=NA)

ylim.bot <- min(pred_HR$lower)
ylim.top <- max(pred_HR$upper)

#  Draw a density chart first to avoid obscuring the line chart below 
dens <- density(new_mydata$omega_ratio2) #  Calculate density 
plot(dens$x,dens$y, col=ggplot2::alpha(grey,0.5), type="l", xlab = "", ylab = "",xaxt="n",yaxt="n")
polygon(dens$x,dens$y,col = ggplot2::alpha(grey,0.5),border = ggplot2::alpha(grey,0.5)) #  Color transparent anti covering line 
axis(side=4, at = pretty(range(dens$y))[-length(pretty(range(dens$y)))])
mtext("Fraction of population (Density)", side=4, line=2)


par(new=TRUE) #  New canvas 
plot(pred_HR$omega_ratio2,pred_HR$yhat,
     xlab = "Ratio of omega-6/omega-3 PUFAs",ylab = "Hazard Ratio (ref = 4)",
     type = "l", ylim = c(0.5,3), xlim = c(0,40),
     col="red",lwd=2) 
lines(pred_HR$omega_ratio2,pred_HR$lower,lty=2,lwd=1.5)
lines(pred_HR$omega_ratio2,pred_HR$upper,lty=2,lwd=1.5)
lines(x=range(pred_HR$omega_ratio2),y=c(1,1),lty=3,col="grey40",lwd=1.3) 
points(as.numeric(refvalue),1,pch=16,cex=1)
title(main = "Overall Cancer Incidence")


legend("topright",lty = c(1,2),col = c("red","black"),
       c("Estimation","95% CI"),
       bty="n",cex=0.8)
dev.off()

################ Prostate cancer incidence - ratio  #####################

# Start fitting the model

fit = cph(Surv(time_to_cancer_year, status) ~ rcs(omega_ratio2,3) + strat(ageg) + TDI + 
              eth_group_new + BMI + Smoking_status + Alcohol_status + IPAQ_activity, data = data_prostate, x=TRUE)
anova(fit)

## p-non-linear
p <-round(anova(fit)[,3],3)

# omega_ratio2   11.61       2   0.0030
# Nonlinear      7.43       1   0.0064

# define the reference number
refvalue <- 4

# Package the data , And specify the reference value .
ddist <<- datadist(data_prostate)
ddist$limits$omega_ratio2[2]<-refvalue
options(datadist="ddist")

pred_HR<-Predict(fit,omega_ratio2,ref.zero=TRUE,fun=exp)


tiff(file.path("/Users/yuchen/Desktop/prostate_ratio.tiff"), units="in", width=4.5, height=4.5, res=300)

#  Sets the background color of the density curve 
grey <- "#D3D3D3"
    #  Draw the left and right double axes baseplot
par(mar = c(5, 4, 4, 4) + 0.3)
par(xpd=NA)

ylim.bot <- min(pred_HR$lower)
ylim.top <- max(pred_HR$upper)

#  Draw a density chart first to avoid obscuring the line chart below 
dens <- density(new_mydata$omega_ratio2) #  Calculate density 
plot(dens$x,dens$y, col=ggplot2::alpha(grey,0.5), type="l", xlab = "", ylab = "",xaxt="n",yaxt="n")
polygon(dens$x,dens$y,col = ggplot2::alpha(grey,0.5),border = ggplot2::alpha(grey,0.5)) #  Color transparent anti covering line 
axis(side=4, at = pretty(range(dens$y))[-length(pretty(range(dens$y)))])
mtext("Fraction of population (Density)", side=4, line=2)


par(new=TRUE) #  New canvas 
plot(pred_HR$omega_ratio2,pred_HR$yhat,
     xlab = "Ratio of omega-6/omega-3 PUFAs",ylab = "Hazard Ratio (ref = 4)",
     type = "l", ylim = c(0.5,3), xlim = c(0,40),
     col="red",lwd=2) 
lines(pred_HR$omega_ratio2,pred_HR$lower,lty=2,lwd=1.5)
lines(pred_HR$omega_ratio2,pred_HR$upper,lty=2,lwd=1.5)
lines(x=range(pred_HR$omega_ratio2),y=c(1,1),lty=3,col="grey40",lwd=1.3) 
points(as.numeric(refvalue),1,pch=16,cex=1)
title(main = "Prostate Cancer Incidence")

legend("topright",lty = c(1,2),col = c("red","black"),
       c("Estimation","95% CI"),
       bty="n",cex=0.8)
dev.off()

################ Breast cancer incidence - ratio  #####################

# Start fitting the model

fit = cph(Surv(time_to_cancer_year, status) ~ rcs(omega_ratio2,3) + strat(ageg)  + TDI + 
              eth_group_new + BMI + Smoking_status + Alcohol_status + IPAQ_activity, data = data_breast, x=TRUE)
anova(fit)

## p-non-linear
p <-round(anova(fit)[,3],3)

# omega_ratio2    5.81       2   0.0548
# Nonlinear      0.90       1   0.3438

# define the reference number
refvalue <- 4

# Package the data , And specify the reference value .
ddist<-datadist(data_breast)
ddist$limits$omega_ratio2[2]<-refvalue
options(datadist="ddist")

pred_HR<-Predict(fit,omega_ratio2,ref.zero=TRUE,fun=exp)


tiff(file.path("/Users/yuchen/Desktop/breast_ratio.tiff"), units="in", width=4.5, height=4.5, res=300)

#  Sets the background color of the density curve 
grey <- "#D3D3D3"
    #  Draw the left and right double axes baseplot
par(mar = c(5, 4, 4, 4) + 0.3)
par(xpd=NA)

ylim.bot <- min(pred_HR$lower)
ylim.top <- max(pred_HR$upper)

#  Draw a density chart first to avoid obscuring the line chart below 
dens <- density(new_mydata$omega_ratio2) #  Calculate density 
plot(dens$x,dens$y, col=ggplot2::alpha(grey,0.5), type="l", xlab = "", ylab = "",xaxt="n",yaxt="n")
polygon(dens$x,dens$y,col = ggplot2::alpha(grey,0.5),border = ggplot2::alpha(grey,0.5)) #  Color transparent anti covering line 
axis(side=4, at = pretty(range(dens$y))[-length(pretty(range(dens$y)))])
mtext("Fraction of population (Density)", side=4, line=2)


par(new=TRUE) #  New canvas 
plot(pred_HR$omega_ratio2,pred_HR$yhat,
     xlab = "Ratio of omega-6/omega-3 PUFAs",ylab = "Hazard Ratio (ref = 4)",
     type = "l", ylim = c(0.5,3), xlim = c(0,40),
     col="red",lwd=2) 
lines(pred_HR$omega_ratio2,pred_HR$lower,lty=2,lwd=1.5)
lines(pred_HR$omega_ratio2,pred_HR$upper,lty=2,lwd=1.5)
lines(x=range(pred_HR$omega_ratio2),y=c(1,1),lty=3,col="grey40",lwd=1.3) 
points(as.numeric(refvalue),1,pch=16,cex=1)
title(main = "Breast Cancer Incidence")

legend("topright",lty = c(1,2),col = c("red","black"),
       c("Estimation","95% CI"),
       bty="n",cex=0.8)
dev.off()


################ Head cancer incidence - ratio  #####################

# Start fitting the model

fit = cph(Surv(time_to_cancer_year, status) ~ rcs(omega_ratio2,3) + strat(ageg) + strat(sex) + TDI + 
              eth_group_new + BMI + Smoking_status + Alcohol_status + IPAQ_activity, data = data_head, x=TRUE)
anova(fit)

## p-non-linear
p <-round(anova(fit)[,3],3)

# omega_ratio2    0.21       2   0.8986
# Nonlinear      0.11       1   0.7443

# define the reference number
refvalue <- 4

# Package the data , And specify the reference value .
ddist <<- datadist(data_head)
ddist$limits$omega_ratio2[2]<-refvalue
options(datadist="ddist")

pred_HR<-Predict(fit,omega_ratio2,ref.zero=TRUE,fun=exp)


tiff(file.path("/Users/yuchen/Desktop/head_ratio.tiff"), units="in", width=4.5, height=4.5, res=300)

#  Sets the background color of the density curve 
grey <- "#D3D3D3"
    #  Draw the left and right double axes baseplot
par(mar = c(5, 4, 4, 4) + 0.3)
par(xpd=NA)

ylim.bot <- min(pred_HR$lower)
ylim.top <- max(pred_HR$upper)

#  Draw a density chart first to avoid obscuring the line chart below 
dens <- density(new_mydata$omega_ratio2) #  Calculate density 
plot(dens$x,dens$y, col=ggplot2::alpha(grey,0.5), type="l", xlab = "", ylab = "",xaxt="n",yaxt="n")
polygon(dens$x,dens$y,col = ggplot2::alpha(grey,0.5),border = ggplot2::alpha(grey,0.5)) #  Color transparent anti covering line 
axis(side=4, at = pretty(range(dens$y))[-length(pretty(range(dens$y)))])
mtext("Fraction of population (Density)", side=4, line=2)


par(new=TRUE) #  New canvas 
plot(pred_HR$omega_ratio2,pred_HR$yhat,
     xlab = "Ratio of omega-6/omega-3 PUFAs",ylab = "Hazard Ratio (ref = 4)",
     type = "l", ylim = c(0.5,3), xlim = c(0,40),
     col="red",lwd=2) 
lines(pred_HR$omega_ratio2,pred_HR$lower,lty=2,lwd=1.5)
lines(pred_HR$omega_ratio2,pred_HR$upper,lty=2,lwd=1.5)
lines(x=range(pred_HR$omega_ratio2),y=c(1,1),lty=3,col="grey40",lwd=1.3) 
points(as.numeric(refvalue),1,pch=16,cex=1)
title(main = "Head & Neck Cancer Incidence")

legend("topright",lty = c(1,2),col = c("red","black"),
       c("Estimation","95% CI"),
       bty="n",cex=0.8)
dev.off()



############################################################################
#           Omega - 3
############################################################################

################ Overall cancer incidence  #####################

# Start fitting the model

fit = cph(Surv(time_to_cancer_year, status) ~ rcs(omega3_pct,3) + strat(ageg) + strat(sex) + TDI + 
              eth_group_new + BMI + Smoking_status + Alcohol_status + IPAQ_activity, data = new_mydata, x=TRUE)
anova(fit)

# omega3_pct       1.56      2   0.4575
# Nonlinear       0.24      1   0.6238

## p-non-linear
p <-round(anova(fit)[,3],3)

# define the reference number
refvalue <- 1

# Package the data , And specify the reference value .
ddist <<- datadist(new_mydata)
ddist$limits$omega3_pct[2]<-refvalue
options(datadist="ddist")

pred_HR<-Predict(fit,omega3_pct,ref.zero=TRUE,fun=exp)


tiff(file.path("/Users/yuchen/Desktop/all_3.tiff"), units="in", width=4.5, height=4.5, res=300)

#  Sets the background color of the density curve 
grey <- "#D3D3D3"
    #  Draw the left and right double axes baseplot
par(mar = c(5, 4, 4, 4) + 0.3)
par(xpd=NA)

ylim.bot <- min(pred_HR$lower)
ylim.top <- max(pred_HR$upper)

#  Draw a density chart first to avoid obscuring the line chart below 
dens <- density(new_mydata$omega3_pct) #  Calculate density 
plot(dens$x,dens$y, col=ggplot2::alpha(grey,0.5), type="l", xlab = "", ylab = "",xaxt="n",yaxt="n")
polygon(dens$x,dens$y,col = ggplot2::alpha(grey,0.5),border = ggplot2::alpha(grey,0.5)) #  Color transparent anti covering line 
axis(side=4, at = pretty(range(dens$y))[-length(pretty(range(dens$y)))])
mtext("Fraction of population (Density)", side=4, line=2)


par(new=TRUE) #  New canvas 
plot(pred_HR$omega3_pct,pred_HR$yhat,
     xlab = "Omega-3 percentage",ylab = "Hazard Ratio (ref = 1%)",
     type = "l", ylim = c(0.5,3),
     col="red",lwd=2) 
lines(pred_HR$omega3_pct,pred_HR$lower,lty=2,lwd=1.5)
lines(pred_HR$omega3_pct,pred_HR$upper,lty=2,lwd=1.5)
lines(x=range(pred_HR$omega3_pct),y=c(1,1),lty=3,col="grey40",lwd=1.3) 
points(as.numeric(refvalue),1,pch=16,cex=1)
title(main = "Overall Cancer Incidence")


legend("topright",lty = c(1,2),col = c("red","black"),
       c("Estimation","95% CI"),
       bty="n",cex=0.8)
dev.off()

################ Prostate cancer incidence  #####################

# Start fitting the model

fit = cph(Surv(time_to_cancer_year, status) ~ rcs(omega3_pct,3) + strat(ageg) + strat(sex) + TDI + 
              eth_group_new + BMI + Smoking_status + Alcohol_status + IPAQ_activity, data = data_prostate, x=TRUE)
anova(fit)

# omega3_pct      3.26       2   0.1961
# Nonlinear      0.38       1   0.5399

## p-non-linear
p <-round(anova(fit)[,3],3)

# define the reference number
refvalue <- 1

# Package the data , And specify the reference value .
ddist <<- datadist(data_prostate)
ddist$limits$omega3_pct[2]<-refvalue
options(datadist="ddist")

pred_HR<-Predict(fit,omega3_pct,ref.zero=TRUE,fun=exp)


tiff(file.path("/Users/yuchen/Desktop/prostate_3.tiff"), units="in", width=4.5, height=4.5, res=300)

#  Sets the background color of the density curve 
grey <- "#D3D3D3"
    #  Draw the left and right double axes baseplot
par(mar = c(5, 4, 4, 4) + 0.3)
par(xpd=NA)

ylim.bot <- min(pred_HR$lower)
ylim.top <- max(pred_HR$upper)

#  Draw a density chart first to avoid obscuring the line chart below 
dens <- density(data_prostate$omega3_pct) #  Calculate density 
plot(dens$x,dens$y, col=ggplot2::alpha(grey,0.5), type="l", xlab = "", ylab = "",xaxt="n",yaxt="n")
polygon(dens$x,dens$y,col = ggplot2::alpha(grey,0.5),border = ggplot2::alpha(grey,0.5)) #  Color transparent anti covering line 
axis(side=4, at = pretty(range(dens$y))[-length(pretty(range(dens$y)))])
mtext("Fraction of population (Density)", side=4, line=2)


par(new=TRUE) #  New canvas 
plot(pred_HR$omega3_pct,pred_HR$yhat,
     xlab = "Omega-3 percentage",ylab = "Hazard Ratio (ref = 1%)",
     type = "l", ylim = c(0.5,3),
     col="red",lwd=2) 
lines(pred_HR$omega3_pct,pred_HR$lower,lty=2,lwd=1.5)
lines(pred_HR$omega3_pct,pred_HR$upper,lty=2,lwd=1.5)
lines(x=range(pred_HR$omega3_pct),y=c(1,1),lty=3,col="grey40",lwd=1.3) 
points(as.numeric(refvalue),1,pch=16,cex=1)
title(main = "Prostate Cancer Incidence")


legend("topright",lty = c(1,2),col = c("red","black"),
       c("Estimation","95% CI"),
       bty="n",cex=0.8)
dev.off()

################ Breast cancer incidence  #####################

# Start fitting the model

fit = cph(Surv(time_to_cancer_year, status) ~ rcs(omega3_pct,3) + strat(ageg) + strat(sex) + TDI + 
              eth_group_new + BMI + Smoking_status + Alcohol_status + IPAQ_activity, data = data_breast, x=TRUE)
anova(fit)

# omega3_pct      0.71       2   0.7017
# Nonlinear      0.02       1   0.8914

## p-non-linear
p <-round(anova(fit)[,3],3)

# define the reference number
refvalue <- 1

# Package the data , And specify the reference value .
ddist <<- datadist(data_breast)
ddist$limits$omega3_pct[2]<-refvalue
options(datadist="ddist")

pred_HR<-Predict(fit,omega3_pct,ref.zero=TRUE,fun=exp)


tiff(file.path("/Users/yuchen/Desktop/breast_3.tiff"), units="in", width=4.5, height=4.5, res=300)

#  Sets the background color of the density curve 
grey <- "#D3D3D3"
    #  Draw the left and right double axes baseplot
par(mar = c(5, 4, 4, 4) + 0.3)
par(xpd=NA)

ylim.bot <- min(pred_HR$lower)
ylim.top <- max(pred_HR$upper)

#  Draw a density chart first to avoid obscuring the line chart below 
dens <- density(data_breast$omega3_pct) #  Calculate density 
plot(dens$x,dens$y, col=ggplot2::alpha(grey,0.5), type="l", xlab = "", ylab = "",xaxt="n",yaxt="n")
polygon(dens$x,dens$y,col = ggplot2::alpha(grey,0.5),border = ggplot2::alpha(grey,0.5)) #  Color transparent anti covering line 
axis(side=4, at = pretty(range(dens$y))[-length(pretty(range(dens$y)))])
mtext("Fraction of population (Density)", side=4, line=2)


par(new=TRUE) #  New canvas 
plot(pred_HR$omega3_pct,pred_HR$yhat,
     xlab = "Omega-3 percentage",ylab = "Hazard Ratio (ref = 1%)",
     type = "l", ylim = c(0.5,3),
     col="red",lwd=2) 
lines(pred_HR$omega3_pct,pred_HR$lower,lty=2,lwd=1.5)
lines(pred_HR$omega3_pct,pred_HR$upper,lty=2,lwd=1.5)
lines(x=range(pred_HR$omega3_pct),y=c(1,1),lty=3,col="grey40",lwd=1.3) 
points(as.numeric(refvalue),1,pch=16,cex=1)
title(main = "Breast Cancer Incidence")


legend("topright",lty = c(1,2),col = c("red","black"),
       c("Estimation","95% CI"),
       bty="n",cex=0.8)
dev.off()

################ Head cancer incidence  #####################

# Start fitting the model

fit = cph(Surv(time_to_cancer_year, status) ~ rcs(omega3_pct,3) + strat(ageg) + strat(sex) + TDI + 
              eth_group_new + BMI + Smoking_status + Alcohol_status + IPAQ_activity, data = data_head, x=TRUE)
anova(fit)

# omega3_pct      1.32       2   0.5167
# Nonlinear      1.06       1   0.3026

## p-non-linear
p <-round(anova(fit)[,3],3)

# define the reference number
refvalue <- 1

# Package the data , And specify the reference value .
ddist <<- datadist(data_head)
ddist$limits$omega3_pct[2]<-refvalue
options(datadist="ddist")

pred_HR<-Predict(fit,omega3_pct,ref.zero=TRUE,fun=exp)


tiff(file.path("/Users/yuchen/Desktop/head_3.tiff"), units="in", width=4.5, height=4.5, res=300)

#  Sets the background color of the density curve 
grey <- "#D3D3D3"
    #  Draw the left and right double axes baseplot
par(mar = c(5, 4, 4, 4) + 0.3)
par(xpd=NA)

ylim.bot <- min(pred_HR$lower)
ylim.top <- max(pred_HR$upper)

#  Draw a density chart first to avoid obscuring the line chart below 
dens <- density(data_head$omega3_pct) #  Calculate density 
plot(dens$x,dens$y, col=ggplot2::alpha(grey,0.5), type="l", xlab = "", ylab = "",xaxt="n",yaxt="n")
polygon(dens$x,dens$y,col = ggplot2::alpha(grey,0.5),border = ggplot2::alpha(grey,0.5)) #  Color transparent anti covering line 
axis(side=4, at = pretty(range(dens$y))[-length(pretty(range(dens$y)))])
mtext("Fraction of population (Density)", side=4, line=2)


par(new=TRUE) #  New canvas 
plot(pred_HR$omega3_pct,pred_HR$yhat,
     xlab = "Omega-3 percentage",ylab = "Hazard Ratio (ref = 1%)",
     type = "l", ylim = c(0.5,3),
     col="red",lwd=2) 
lines(pred_HR$omega3_pct,pred_HR$lower,lty=2,lwd=1.5)
lines(pred_HR$omega3_pct,pred_HR$upper,lty=2,lwd=1.5)
lines(x=range(pred_HR$omega3_pct),y=c(1,1),lty=3,col="grey40",lwd=1.3) 
points(as.numeric(refvalue),1,pch=16,cex=1)
title(main = "Head & Neck Cancer Incidence")


legend("topright",lty = c(1,2),col = c("red","black"),
       c("Estimation","95% CI"),
       bty="n",cex=0.8)
dev.off()




############################################################################
#           Omega - 6
############################################################################

################ Overall cancer incidence  #####################

# Start fitting the model

fit = cph(Surv(time_to_cancer_year, status) ~ rcs(omega6_pct,3) + strat(ageg) + strat(sex) + TDI + 
              eth_group_new + BMI + Smoking_status + Alcohol_status + IPAQ_activity, data = new_mydata, x=TRUE)
anova(fit)

# omega6_pct       5.27      2   0.0716
# Nonlinear       0.00      1   0.9722

## p-non-linear
p <-round(anova(fit)[,3],3)

# define the reference number
refvalue <- 20

# Package the data , And specify the reference value .
ddist <<- datadist(new_mydata)
ddist$limits$omega6_pct[2]<-refvalue
options(datadist="ddist")

pred_HR<-Predict(fit,omega6_pct,ref.zero=TRUE,fun=exp)


tiff(file.path("/Users/yuchen/Desktop/all_6.tiff"), units="in", width=4.5, height=4.5, res=300)

#  Sets the background color of the density curve 
grey <- "#D3D3D3"
    #  Draw the left and right double axes baseplot
par(mar = c(5, 4, 4, 4) + 0.3)
par(xpd=NA)

ylim.bot <- min(pred_HR$lower)
ylim.top <- max(pred_HR$upper)

#  Draw a density chart first to avoid obscuring the line chart below 
dens <- density(new_mydata$omega6_pct) #  Calculate density 
plot(dens$x,dens$y, col=ggplot2::alpha(grey,0.5), type="l", xlab = "", ylab = "",xaxt="n",yaxt="n")
polygon(dens$x,dens$y,col = ggplot2::alpha(grey,0.5),border = ggplot2::alpha(grey,0.5)) #  Color transparent anti covering line 
axis(side=4, at = pretty(range(dens$y))[-length(pretty(range(dens$y)))])
mtext("Fraction of population (Density)", side=4, line=2)


par(new=TRUE) #  New canvas 
plot(pred_HR$omega6_pct,pred_HR$yhat,
     xlab = "Omega-6 percentage",ylab = "Hazard Ratio (ref = 1%)",
     type = "l", ylim = c(0.5,3),
     col="red",lwd=2) 
lines(pred_HR$omega6_pct,pred_HR$lower,lty=2,lwd=1.5)
lines(pred_HR$omega6_pct,pred_HR$upper,lty=2,lwd=1.5)
lines(x=range(pred_HR$omega6_pct),y=c(1,1),lty=3,col="grey40",lwd=1.3) 
points(as.numeric(refvalue),1,pch=16,cex=1)
title(main = "Overall Cancer Incidence")


legend("topright",lty = c(1,2),col = c("red","black"),
       c("Estimation","95% CI"),
       bty="n",cex=0.8)
dev.off()

################ Prostate cancer incidence  #####################

# Start fitting the model

fit = cph(Surv(time_to_cancer_year, status) ~ rcs(omega6_pct,3) + strat(ageg) + strat(sex) + TDI + 
              eth_group_new + BMI + Smoking_status + Alcohol_status + IPAQ_activity, data = data_prostate, x=TRUE)
anova(fit)

# omega6_pct      3.69       2   0.1584
# Nonlinear      3.41       1   0.0646

## p-non-linear
p <-round(anova(fit)[,3],3)

# define the reference number
refvalue <- 20

# Package the data , And specify the reference value .
ddist <<- datadist(data_prostate)
ddist$limits$omega6_pct[2]<-refvalue
options(datadist="ddist")

pred_HR<-Predict(fit,omega6_pct,ref.zero=TRUE,fun=exp)


tiff(file.path("/Users/yuchen/Desktop/prostate_6.tiff"), units="in", width=4.5, height=4.5, res=300)

#  Sets the background color of the density curve 
grey <- "#D3D3D3"
    #  Draw the left and right double axes baseplot
par(mar = c(5, 4, 4, 4) + 0.3)
par(xpd=NA)

ylim.bot <- min(pred_HR$lower)
ylim.top <- max(pred_HR$upper)

#  Draw a density chart first to avoid obscuring the line chart below 
dens <- density(data_prostate$omega6_pct) #  Calculate density 
plot(dens$x,dens$y, col=ggplot2::alpha(grey,0.5), type="l", xlab = "", ylab = "",xaxt="n",yaxt="n")
polygon(dens$x,dens$y,col = ggplot2::alpha(grey,0.5),border = ggplot2::alpha(grey,0.5)) #  Color transparent anti covering line 
axis(side=4, at = pretty(range(dens$y))[-length(pretty(range(dens$y)))])
mtext("Fraction of population (Density)", side=4, line=2)


par(new=TRUE) #  New canvas 
plot(pred_HR$omega6_pct,pred_HR$yhat,
     xlab = "Omega-6 percentage",ylab = "Hazard Ratio (ref = 1%)",
     type = "l", ylim = c(0.5,3),
     col="red",lwd=2) 
lines(pred_HR$omega6_pct,pred_HR$lower,lty=2,lwd=1.5)
lines(pred_HR$omega6_pct,pred_HR$upper,lty=2,lwd=1.5)
lines(x=range(pred_HR$omega6_pct),y=c(1,1),lty=3,col="grey40",lwd=1.3) 
points(as.numeric(refvalue),1,pch=16,cex=1)
title(main = "Prostate Cancer Incidence")


legend("topright",lty = c(1,2),col = c("red","black"),
       c("Estimation","95% CI"),
       bty="n",cex=0.8)
dev.off()

################ Breast cancer incidence  #####################

# Start fitting the model

fit = cph(Surv(time_to_cancer_year, status) ~ rcs(omega6_pct,3) + strat(ageg) + strat(sex) + TDI + 
              eth_group_new + BMI + Smoking_status + Alcohol_status + IPAQ_activity, data = data_breast, x=TRUE)
anova(fit)

# omega6_pct      2.86       2   0.2387
# Nonlinear      0.09       1   0.7596

## p-non-linear
p <-round(anova(fit)[,3],3)

# define the reference number
refvalue <- 20

# Package the data , And specify the reference value .
ddist <<- datadist(data_breast)
ddist$limits$omega6_pct[2]<-refvalue
options(datadist="ddist")

pred_HR<-Predict(fit,omega6_pct,ref.zero=TRUE,fun=exp)


tiff(file.path("/Users/yuchen/Desktop/breast_6.tiff"), units="in", width=4.5, height=4.5, res=300)

#  Sets the background color of the density curve 
grey <- "#D3D3D3"
    #  Draw the left and right double axes baseplot
par(mar = c(5, 4, 4, 4) + 0.3)
par(xpd=NA)

ylim.bot <- min(pred_HR$lower)
ylim.top <- max(pred_HR$upper)

#  Draw a density chart first to avoid obscuring the line chart below 
dens <- density(data_breast$omega6_pct) #  Calculate density 
plot(dens$x,dens$y, col=ggplot2::alpha(grey,0.5), type="l", xlab = "", ylab = "",xaxt="n",yaxt="n")
polygon(dens$x,dens$y,col = ggplot2::alpha(grey,0.5),border = ggplot2::alpha(grey,0.5)) #  Color transparent anti covering line 
axis(side=4, at = pretty(range(dens$y))[-length(pretty(range(dens$y)))])
mtext("Fraction of population (Density)", side=4, line=2)


par(new=TRUE) #  New canvas 
plot(pred_HR$omega6_pct,pred_HR$yhat,
     xlab = "Omega-6 percentage",ylab = "Hazard Ratio (ref = 1%)",
     type = "l", ylim = c(0.5,3),
     col="red",lwd=2) 
lines(pred_HR$omega6_pct,pred_HR$lower,lty=2,lwd=1.5)
lines(pred_HR$omega6_pct,pred_HR$upper,lty=2,lwd=1.5)
lines(x=range(pred_HR$omega6_pct),y=c(1,1),lty=3,col="grey40",lwd=1.3) 
points(as.numeric(refvalue),1,pch=16,cex=1)
title(main = "Breast Cancer Incidence")


legend("topright",lty = c(1,2),col = c("red","black"),
       c("Estimation","95% CI"),
       bty="n",cex=0.8)
dev.off()

################ Head cancer incidence  #####################

# Start fitting the model

fit = cph(Surv(time_to_cancer_year, status) ~ rcs(omega6_pct,3) + strat(ageg) + strat(sex) + TDI + 
              eth_group_new + BMI + Smoking_status + Alcohol_status + IPAQ_activity, data = data_head, x=TRUE)
anova(fit)

# omega6_pct      4.64       2   0.0982
# Nonlinear      1.12       1   0.2895

## p-non-linear
p <-round(anova(fit)[,3],3)

# define the reference number
refvalue <- 20

# Package the data , And specify the reference value .
ddist <<- datadist(data_head)
ddist$limits$omega6_pct[2]<-refvalue
options(datadist="ddist")

pred_HR<-Predict(fit,omega6_pct,ref.zero=TRUE,fun=exp)


tiff(file.path("/Users/yuchen/Desktop/head_6.tiff"), units="in", width=4.5, height=4.5, res=300)

#  Sets the background color of the density curve 
grey <- "#D3D3D3"
    #  Draw the left and right double axes baseplot
par(mar = c(5, 4, 4, 4) + 0.3)
par(xpd=NA)

ylim.bot <- min(pred_HR$lower)
ylim.top <- max(pred_HR$upper)

#  Draw a density chart first to avoid obscuring the line chart below 
dens <- density(data_head$omega6_pct) #  Calculate density 
plot(dens$x,dens$y, col=ggplot2::alpha(grey,0.5), type="l", xlab = "", ylab = "",xaxt="n",yaxt="n")
polygon(dens$x,dens$y,col = ggplot2::alpha(grey,0.5),border = ggplot2::alpha(grey,0.5)) #  Color transparent anti covering line 
axis(side=4, at = pretty(range(dens$y))[-length(pretty(range(dens$y)))])
mtext("Fraction of population (Density)", side=4, line=2)


par(new=TRUE) #  New canvas 
plot(pred_HR$omega6_pct,pred_HR$yhat,
     xlab = "Omega-6 percentage",ylab = "Hazard Ratio (ref = 1%)",
     type = "l", ylim = c(0.5,3),
     col="red",lwd=2) 
lines(pred_HR$omega6_pct,pred_HR$lower,lty=2,lwd=1.5)
lines(pred_HR$omega6_pct,pred_HR$upper,lty=2,lwd=1.5)
lines(x=range(pred_HR$omega6_pct),y=c(1,1),lty=3,col="grey40",lwd=1.3) 
points(as.numeric(refvalue),1,pch=16,cex=1)
title(main = "Head & Neck Cancer Incidence")


legend("topright",lty = c(1,2),col = c("red","black"),
       c("Estimation","95% CI"),
       bty="n",cex=0.8)
dev.off()


