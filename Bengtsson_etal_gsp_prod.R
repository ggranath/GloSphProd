##########################################################################################
# Script to analyse data and produce the results presetned in                            #
# " Environmental drivers of Sphagnum growth in peatlands across the Holarctic region"   #
#                                                                                        #
# Contact: Gustaf.Granath@gmail.com                                                      #
##########################################################################################

# Load extra functions (R2, plot etc) ####
source("gsp_extra_func.R")

# Load data
dat <- read.csv("Bengtsson_etal_gsp_prod.csv")

# Remove two site in China (inconsistent sampling) and one site that do not have N concentration
rem <- which(dat$Site == "Hani"| dat$Site == "Mangui"| dat$Site == "Verh-Tarka")
dat <- dat[-rem, ]


#Make df in long format
# and fix variables
library(tidyr)
data_long <- gather(dat, year, LI, c(LI13, LI14), factor_key=TRUE)
levels(data_long$year)<-c("yr2013", "yr2014")

data_long$grow_days <- c(dat$days2013, dat$days2014)
data_long$LI_d <- c(dat$LI_per_day2013, dat$LI_per_day2014)
data_long$prod <- c(dat$Prod13, dat$Prod14)
data_long$prod_d <- c(dat$prod_per_day2013, dat$prod_per_day2014)
data_long$temp <- c(dat$temp13, dat$temp14)
data_long$evap <- c(dat$evap13, dat$evap14)
data_long$evap_d <- c(dat$ev_d13, dat$ev_d14)
data_long$prec <- c(dat$precip13, dat$precip14)
data_long$prec_d <- c(dat$prec_d13, dat$prec_d14)
data_long$norain <- c(dat$norain13.mean, dat$norain14.mean)
data_long$par <- c(dat$par13, dat$par14)
data_long$par_d <- c(dat$par_d13, dat$par_d14)
data_long$cover <- c(dat$cover13, dat$cover14)
data_long$hwt <- c(dat$HWT13, dat$HWT14)
data_long$BD <- c(dat$BD13, dat$BD14)
data_long$numden <- c(dat$num.den13, dat$num.den14)
data_long$ev_pre_d <- data_long$prec_d - data_long$evap_d # precipitation minus evaporation per day
data_long$ev_pre <- data_long$prec - data_long$evap # precipitation minus evaporation 

# Add site mean nitrogen concentration and NP ratio
data_long$id <- 1:nrow(data_long)
data_long <- data_long %>% 
  group_by(Site, Species) %>% 
  mutate(Nmean = mean(N_per, na.rm=TRUE),
         NPmean = mean(NP, na.rm=TRUE)) %>% arrange(id) %>% as.data.frame()

# Loop to: scale for each response variable after removing NA rows,
# and save each data frame in a list named 'resp. for later use
responses <- c("LI_d", "LI", "prod_d", "prod", "BD", "numden")
resp <-list()
for (i in 1:length(responses)) {
temp <- data_long[complete.cases(data_long[ , c(responses[i], "prec_d", "evap_d", "par_d",
                                                       "prec", "temp", "evap", "par", "hwt", "cover", "norain", 
                                                       "ndep_Lam13", "Nmean", "NPmean"),]),]

temp.sc <- scale(temp[ , c("prec_d", "evap_d", "par_d", "ev_pre", "ev_pre_d",
                                  "prec", "temp", "evap", "par", "hwt", "cover", "norain", 
                                  "ndep_Lam13", "Nmean", "NPmean"),])
colnames(temp.sc) <- c("pr_d", "ev_d", "pa_d","ev_pre","ev_pre_d",
                      "pr", "tem", "ev", "pa", "wt", "cov", "nor", "ndeL", "Nm", "NPm")
resp[[i]] <- as.data.frame(cbind(temp, temp.sc))
}
names(resp) <- responses # name list

# Run models ####

# Mixed models. Sites as random, 
# Using lmer and Powertransform to transform data so residuals 
# are normaly distributed, predictors are scaled
library(lme4)
library(car) #needed for powerTransform()
library(nlme)

#___LI per day####
# first get lambda for Box-cox transformation, response is too small, so *100
# make all values positive for easier transformation
min.val <- min(resp[["LI_d"]]$LI_d, na.rm=T)*-1+0.0001 # get value to add for positive response
mod1 <- lmer(((LI_d+min.val)*100) ~ Species+year+ Species*pr_d  + Species*tem + Species*ev_d + 
               Species*pa_d + Species*wt + Species*cov + 
             Species*nor + Species*ndeL + Species*Nm+
             (year-1|Site),
            data = resp[["LI_d"]])
plot(mod1) #increasing variance
summary(l1<-powerTransform(mod1, family = "bcPower"))#bcnPower neccessary with neg. values. Lambda=0.30

#transform response
resp[["LI_d"]]$LId_trans<-((resp[["LI_d"]]$LI_d+min.val)*100)^(l1$lambda)
mod1.2<-lmer(LId_trans ~ Species+year+Species*pr_d  + Species*tem + Species*ev_d + Species*pa_d + 
               Species*wt + Species*cov + 
               Species*nor + Species*ndeL + Species*Nm+
               (year-1|Site),
         data = resp[["LI_d"]])

        #data_long.na2<- mod1.2@frame # save data w/o NAs
#Anova to check if interactions are interesting library(car)
Anova(mod1.2, test.statistic = "F")
summary(mod1.2)
plot(mod1.2)

#to get coefficients per species
mod1.2.1<-lmer(LId_trans ~ Species+year+Species/pr_d  + Species/tem + Species/ev_d + Species/pa_d + 
                 Species/wt + Species/cov + 
               Species/nor + Species/ndeL + Species/Nm+
               (year-1|Site),
             data = resp[["LI_d"]])
summary(mod1.2.1)

#LI_d null model
mod1n<-lmer((LId_trans)~1 + (1|Site),
            data = resp[["LI_d"]])
#plot(mod1n)
summary(mod1n)

#______Excluding outliers from model with 2-way interactions ####
outliersLId<-resp[["LI_d"]][which(residuals(mod1.2)<(-1.5)),]#which are the outliers?
data_long1<-resp[["LI_d"]][which(residuals(mod1.2)>(-1.5)),]
modxx1<-lmer((LId_trans) ~ Species/year+ Species/pr_d  + Species/tem + Species/ev_d + Species/pa_d + 
               Species/wt + Species/cov + 
              Species/nor + Species/nde + Species/Nm + (1|Site), 
             data = data_long1, na.action = na.omit)

modx1<-lmer((LId_trans) ~ Species*year+ Species*pr_d  + Species*tem + Species*ev_d + Species*pa_d + 
              Species*wt + Species*cov + 
              Species*nor + Species*nde + Species*Nm + (1|Site),
            data = data_long1, na.action = na.omit)
plot(modx1)

#Anova to check if interactions are interesting library(car)
summary(modxx1)
Anova(modx1, test.statistic = "F")
summary(modx1)

#Model with only main effects
mod1.5<-lmer(LId_trans ~ Species + year + pr_d + tem +ev_d + pa_d + wt + cov + nor+ ndeL + Nm + 
              (year-1|Site), 
               data = resp[["LI_d"]])
plot(mod1.5)
summary(mod1.5)
Anova(mod1.5, test.statistic = "F")

# Check residuals
res.test <- mod4.2@frame
res.test$resid <- residuals(mod4.2,type="pearson")
# numeric predictors
res.test %>%
  select_if(is.numeric) %>%
  gather(preds, value, 2:(ncol(.)-1)) %>% 
  ggplot(aes(x=value, y=resid)) +
  geom_point() +
  facet_wrap(~ preds, scales = "free")
# categorical predictors also
res.test %>%
  select(-Site) %>%
  gather(preds, value, 2:(ncol(.)-1)) %>% 
  ggplot(aes(x=value, y=resid)) +
  geom_point() +
  facet_wrap(~ preds, scales = "free")

# check influential points
ggplot(data.frame(lev=hatvalues(mod4.2),pearson=residuals(mod4.2,type="pearson")),
       aes(x=lev,y=pearson)) +
  geom_point() +
  theme_bw()

#identify points with high leverage, set a cutoff
res.test[which(hatvalues(mod2.2) >= .25),]


#______Excluding outliers from model with no interactions####
outliersLId2 <- resp[["LI_d"]][which(residuals(mod1.5)<(-1.5)),]

data_long2<-resp[["LI_d"]][which(residuals(mod1.5)>(-1.5)),]

modx2<-lmer((LId_trans) ~ Species + year + pr_d  + tem + ev_d + pa_d + wt + cov + nor+ nde + Nm+
              (1|Site), 
            data = data_long2, na.action = na.omit)
plot(modx2)
summary(modx2)
Anova(modx2, test.statistic = "F")##end outlier test


#______Prec-evap mixed model####
# Test if a index works better than evap and precip as separate variables
min.val <- min(resp[["LI_d"]]$LI_d, na.rm=T)*-1+0.0001 # get value to add for positive response
mod1_evpre <- lmer(((LI_d+min.val)*100) ~ Species+year+ Species*ev_pre_d  + Species*tem + 
               Species*pa_d + Species*wt + Species*cov + 
               Species*nor + Species*ndeL + Species*Nm+
               (year-1|Site),
             data = resp[["LI_d"]])
plot(mod1_evpre) #increasing variance
summary(l1<-powerTransform(mod1_evpre, family = "bcPower"))#bcnPower neccessary with neg. values. Lambda=0.32

#transform response
resp[["LI_d"]]$LId_trans<-((resp[["LI_d"]]$LI_d+min.val)*100)^(l1$lambda)
mod2_evpre<-lmer(LId_trans ~ Species+year+Species*ev_pre_d  + Species*tem + Species*pa_d + 
               Species*wt + Species*cov + 
               Species*nor + Species*ndeL + Species*Nm+
               (year-1|Site),
             data = resp[["LI_d"]])
#Anova to check if interactions are interesting library(car)
Anova(mod2_evpre, test.statistic = "F")

#to get coefficients per species
mod2.1_evpre<-lmer(LId_trans ~ Species+year+Species/ev_pre_d  + Species/tem +  Species/pa_d + 
                 Species/wt + Species/cov + 
                 Species/nor + Species/ndeL + Species/Nm+
                 (year-1|Site),
               data = resp[["LI_d"]])
summary(mod2.1_evpre)

#Model with only main effects
mod3_evpre<-lmer(LId_trans ~ Species + year + ev_pre_d + tem + pa_d + wt + cov + nor+ ndeL + Nm + 
               (year-1|Site), 
             data = resp[["LI_d"]])
plot(mod3_evpre)
summary(mod3_evpre)
Anova(mod3_evpre, test.statistic = "F")


#___NPP per day####
# first get lambda for Box-cox transformation, response is too small, so *100
# make all values positive for easier transformation

min.val <- min(resp[["prod_d"]]$prod_d)*-1+0.001 # get value to add for positive response
mod2<-lmer(( (prod_d+min.val)*100) ~ Species+year+ Species*pr_d  + Species*tem + Species*ev_d + 
             Species*pa_d + Species*wt + Species*cov + 
             Species*nor + Species*ndeL + Species*Nm+
             (year-1|Site), data = resp[["prod_d"]])

plot(mod2)
summary(l2<-powerTransform(mod2, family = "bcPower"))#bcnPower neccessary with neg. values. Lambda=0.32

#transform response
resp[["prod_d"]]$prod.d_trans <- ((resp[["prod_d"]]$prod_d+min.val)*100)^(l2$lambda)
mod2.2<-lmer(prod.d_trans ~ Species+year+Species*pr_d  + Species*tem + Species*ev_d + 
               Species*pa_d + Species*wt + Species*cov + 
               Species*nor + Species*ndeL + Species*Nm+
               (year-1|Site), 
             data= resp[["prod_d"]])

#data_long.prod.d.na<- mod2.2@frame # save data w/o NAs
plot(mod2.2)
#to get coefficients
mod2.2.1<-lmer(prod.d_trans ~ Species+year+Species/pr_d  + Species/tem + Species/ev_d + Species/pa_d + 
                 Species/wt + Species/cov + 
                 Species/nor + Species/ndeL + Species/Nm+
                 (year-1|Site), 
               data = resp[["prod_d"]])
summary(mod2.2.1)
Anova(mod2.2, test.statistic = "F")
summary(mod2.2)

#prod_d null model
mod2n<-lmer((prod.d_trans)~1 + (1|Site),
            data = resp[["prod_d"]])
#plot(mod2n)
#summary(mod2n)

#Run without interactions to get coefficients
mod2.5<-lmer(prod.d_trans ~ Species + year + pr_d  + tem + ev_d + pa_d + wt + cov + nor+ ndeL + Nm+
               (year-1|Site), 
             data = resp[["prod_d"]])
plot(mod2.5)
summary(mod2.5)
Anova(mod2.5, test.statistic = "F")

#_______Prec-evap Mixed model####
min.val <- min(resp[["prod_d"]]$prod_d)*-1+0.001 # get value to add for positive response
modE_P<-lmer(( (prod_d+min.val)*100) ~ Species+year+ Species*ev_pre_d  + Species*tem + 
             Species*pa_d + Species*wt + Species*cov + 
             Species*nor + Species*ndeL + Species*Nm+
             (year-1|Site), data = resp[["prod_d"]])

plot(modE_P)
summary(l2<-powerTransform(modE_P, family = "bcPower"))#bcnPower neccessary with neg. values. Lambda=0.32

#transform response
resp[["prod_d"]]$prod.d_trans <- ((resp[["prod_d"]]$prod_d+min.val)*100)^(l2$lambda)

modE_P2<-lmer(prod.d_trans ~ Species+year+Species*ev_pre_d  + Species*tem + 
               Species*pa_d + Species*wt + Species*cov + 
               Species*nor + Species*ndeL + Species*Nm+
               (year-1|Site), 
             data= resp[["prod_d"]])

Anova(modE_P2, test.statistic = "F")


modE_P4<-lmer(prod.d_trans ~ Species+year + Species/ev_pre_d  + Species/tem + Species/pa_d + 
                Species/wt + Species/cov + Species/nor + Species/ndeL + Species/Nm+
                (year-1|Site), data = resp[["prod_d"]])
summary(modE_P4)

#Run without interactions to get coefficients
modE_P5<-lmer(prod.d_trans ~ Species + year + ev_pre_d  + tem + pa_d + wt + cov + nor+ ndeL + Nm+ 
               (year-1|Site),
              data = resp[["prod_d"]])
plot(modE_P5)
summary(modE_P5)
Anova(modE_P5, test.statistic = "F")

#___LI per season ####
min.val <- min(resp[["LI"]]$LI)*-1+0.001 # get value to add for positive response
mod3<-lmer((LI+min.val) ~ Species+year + Species*pr + Species*tem + Species*ev + Species*pa + 
            Species*wt + Species*cov + Species*nor + Species*ndeL + Species*Nm +
           (year-1|Site), data = resp[["LI"]])
plot(mod3)

summary(l3<-powerTransform(mod3, family = "bcPower"))#bcnPower neccessary with neg. values lambda=0.38

resp[["LI"]]$LI_trans <- (resp[["LI"]]$LI+min.val)^(l3$lambda)

mod3.2<-lmer(LI_trans ~ Species+year + Species*pr + Species*tem + Species*ev + Species*pa + 
               Species*wt + Species*cov + Species*nor + Species*ndeL + Species*Nm +
               (year-1|Site), data = resp[["LI"]])

#LI only PAR variable
mod3.2par<-lmer(LI_trans ~ Species+year + Species*pa +
               (year-1|Site), data = resp[["LI"]])
summary(mod3.2par)

#LI season null model
mod3n<-lmer((LI_trans)~1 + (year-1|Site),
            data = resp[["LI"]])
summary(mod3.2)
summary(mod3n)

#Calculate within and between site variance 
mod3n<-lmer((LI_trans)~1 + (1|Site),
            data = resp[["LI"]],subset = Species=="S.magellanicum")

# code to test if random slope is better
mod3.2<-lmer(LI_trans ~ Species+year + Species*pr + Species*tem + Species*ev + Species*pa + 
               Species*wt + Species*cov + Species*nor + Species*ndeL + Species*Nm +
               (year-1|Site), data = resp[["LI"]])
chistat <- anova(mod3.2b, mod3.2)[2,"Chisq"] 
0.5 * pchisq(chistat, 1, lower.tail=FALSE) + 0.5 * pchisq(chistat, 2, lower.tail=FALSE)

plot(mod3.2)

# get coefs
mod3.2.1<-lmer(LI_trans ~ Species+year + Species/pr + Species/tem + Species/ev + Species/pa + 
                 Species/wt + Species/cov + Species/nor + Species/ndeL + Species/Nm +
               (year-1|Site), data = resp[["LI"]])

summary(mod3.2.1)
Anova(mod3.2, test.statistic = "F")
summary(mod3.2)


#no interactions
mod3.4<-lmer(LI_trans ~ Species + year + pr + tem + ev + pa + wt + cov + nor + ndeL + Nm 
             + (year-1|Site), 
             data = resp[["LI"]])

plot(mod3.4)
summary(mod3.4)
Anova(mod3.4, test.statistic = "F")

#______Excluding outliers from model with no interactions####
outliersLIseas<-resp[["LI"]][which(residuals(mod3.2)<(-1.5)),]

data_long2<-resp[["LI"]][which(residuals(mod3.2)>(-1.5)),]

modx3.1<-lmer(LI_trans ~ Species/year + Species/pr + Species/tem + Species/ev + Species/pa + 
                Species/wt + Species/cov + Species/nor + Species/nde + Species/Nm +
                (1|Site), data = data_long2, na.action = na.omit)
plot(modx3.1)
summary(modx3.1)
Anova(modx3.1, test.statistic = "F")

#no interactions
modx3<-lmer(LI_trans ~ Species + year + pr + tem + ev + pa + wt + cov + nor + nde + Nm +
              (1|Site), data = data_long2, na.action = na.omit)
plot(modx3)
summary(modx3)
Anova(modx3, test.statistic = "F")


#______Prec-evap mixed model####
min.val <- min(resp[["LI"]]$LI)*-1+0.001 # get value to add for positive response
LI_P_E1<-lmer((LI+min.val) ~ Species+year + Species*ev_pre + Species*tem + Species*pa + 
             Species*wt + Species*cov + Species*nor + Species*ndeL + Species*Nm +
             (year-1|Site), data = resp[["LI"]])
plot(LI_P_E1)

summary(l3<-powerTransform(LI_P_E1, family = "bcPower"))#bcnPower neccessary with neg. values lambda=0.34

resp[["LI"]]$LI_trans <- (resp[["LI"]]$LI+min.val)^(l3$lambda)

LI_P_E<-lmer(LI_trans ~ Species+year + Species*ev_pre  + Species*tem + Species*pa_d + 
               Species*wt + Species*cov + Species*nor + Species*ndeL + Species*Nm+
               (year-1|Site), data = resp[["LI"]])

plot(LI_P_E)

LI_P_E2<-lmer(LI_trans ~ Species+year + Species/ev_pre  + Species/tem + Species/pa_d + 
                Species/wt + Species/cov + Species/nor + Species/ndeL + Species/Nm+
                (year-1|Site), data = resp[["LI"]])
summary(LI_P_E2)
Anova(LI_P_E, test.statistic = "F")

#Run without interactions to get coefficients
LI_P_E3<-lmer(LI_trans ~ Species + year + ev_pre  + tem + pa_d + wt + cov + nor+ ndeL + Nm+ 
               (year-1|Site),
             data = resp[["LI"]])
plot(LI_P_E3)
summary(LI_P_E3)
Anova(LI_P_E3, test.statistic = "F")

#___NPP per season ####
min.val <- min(resp[["prod"]]$prod)*-1+0.001 # get value to add for positive response
mod4<-lmer((prod+min.val) ~ Species+year + Species*pr + Species*tem + Species*ev + Species*pa + 
             Species*wt + Species*cov + Species*nor + Species*ndeL + Species*Nm +
             (year-1|Site), data = resp[["prod"]])
plot(mod4)
summary(l4<-powerTransform(mod4, family = "bcPower"))#bcnPower neccessary with neg. values 
resp[["prod"]]$prod_trans<-(resp[["prod"]]$prod+min.val)^(l4$lambda)

mod4.2<-lmer(prod_trans ~ Species+year + Species*pr + Species*tem + Species*ev + Species*pa + 
               Species*wt + Species*cov + Species*nor + Species*ndeL + Species*Nm +
               (year-1|Site), data = resp[["prod"]][-(which(resp[["prod"]]$orgOrder ==177)),])
plot(mod4.2)

# coefs
mod4.2.1<-lmer(prod_trans ~ Species+year + Species/pr + Species/tem + Species/ev + Species/pa + 
                 Species/wt + Species/cov + Species/nor + Species/ndeL + Species/Nm +
                 (year-1|Site), data = resp[["prod"]])

summary(mod4.2.1)
Anova(mod4.2, test.statistic = "F")
summary(mod4.2)

# no interactions
mod4.4<-lmer(prod_trans ~ Species + year + pr + tem + ev + pa + 
               wt + cov + nor + ndeL + Nm + (year-1|Site), 
             data = resp[["prod"]][-(which(resp[["prod"]]$orgOrder ==177)),])

plot(mod4.4)
summary(mod4.4)
Anova(mod4.4, test.statistic = "F")

#prod season null model
mod4n<-lmer((prod_trans)~1 + (year-1|Site),
            data = resp[["prod"]][-(which(resp[["prod"]]$orgOrder ==177)),])
#plot(mod4n)
summary(mod4n)

#Calculate site variance
#Overall
mod4n<-lmer((prod_trans)~1 + (1|Site),
            data = resp[["prod"]][-(which(resp[["prod"]]$orgOrder ==177)),])
summary(mod4n)
mod4n<-lmer((prod_trans)~1 + (1|Site),
            data = resp[["prod"]][-(which(resp[["prod"]]$orgOrder ==177)),],subset = Species=="S.fuscum")

#______Excluding outliers from model with no interactions####
outliersprodseas<-resp[["prod"]][which(residuals(mod4.2)<(-1.5)),]

data_long4<-resp[["prod"]][which(residuals(mod4.2)>(-1.5)),]

modx4.1<-lmer(prod_trans ~ Species/year + Species/pr + Species/tem + Species/ev + 
                Species/pa + Species/wt + Species/cov + Species/nor + Species/nde + Species/Nm +
                (year-1|Site), 
              data = data_long4)
plot(modx4.1)
summary(modx4.1)
Anova(modx4.1, test.statistic = "F")

#no interactions
modx4<-lmer(prod_trans ~ Species + year + pr + tem + ev + pa + wt + cov + nor + nde + Nm +
              (year-1|Site), 
            data = data_long4)
plot(modx4)
summary(modx4)
Anova(modx4, test.statistic = "F")#


#______Prec-Evap mixed model####
min.val <- min(resp[["prod"]]$prod)*-1+0.001 # get value to add for positive response
NPP_E_P<-lmer((prod+min.val) ~ Species+year + Species*ev_pre + Species*tem + Species*ev + Species*pa + 
             Species*wt + Species*cov + Species*nor + Species*ndeL + Species*Nm +
             (year-1|Site), data = resp[["prod"]])
plot(NPP_E_P)
summary(l4<-powerTransform(NPP_E_P, family = "bcPower"))#bcnPower neccessary with neg. values 
resp[["prod"]]$prod_trans<-(resp[["prod"]]$prod+min.val)^(l4$lambda)

NPP_E_P2<-lmer(prod_trans ~ Species+year + Species*ev_pre + Species*tem + Species*pa + 
               Species*wt + Species*cov + Species*nor + Species*ndeL + Species*Nm +
               (year-1|Site), data = resp[["prod"]])
plot(NPP_E_P2)

NPP_E_P3<-lmer(prod_trans ~ Species+year + Species/ev_pre + Species/tem + Species/pa + 
                 Species/wt + Species/cov + Species/nor + Species/ndeL + Species/Nm +
                 (year-1|Site), data = resp[["prod"]])

summary(NPP_E_P3)
Anova(NPP_E_P2, test.statistic = "F")

#no interactions
NPP_E_P4<-lmer(prod_trans ~ Species + year + ev_pre + tem + pa + 
               wt + cov + nor + ndeL + Nm + (year-1|Site), 
             data = resp[["prod"]])

plot(NPP_E_P4)
summary(NPP_E_P4)
Anova(NPP_E_P4, test.statistic = "F")


#___Bulk density####
mod5<-lmer(BD ~ Species*year + Species*pr + Species*tem + Species*ev + Species*pa + 
             Species*wt + Species*cov + Species*nor + Species*nde + Species*Nm +
             (1|Site), data = resp[["BD"]])
plot(mod5)
summary(l5<-powerTransform(mod5, family = "bcPower"))#bcnPower neccessary with neg. values lambda=0.27
resp[["BD"]]$BD_trans<-resp[["BD"]]$BD^(l5$lambda)

mod5.2<-lmer(BD_trans ~ Species*year + Species*pr + Species*tem + Species*ev + Species*pa + 
               Species*wt + Species*cov + Species*nor + Species*nde + Species*Nm +
               (1|Site), data = resp[["BD"]])
plot(mod5.2)

mod5.2.1<-lmer(BD_trans ~ Species/year + Species/pr + Species/tem + Species/ev + Species/pa + 
                 Species/wt + Species/cov + Species/nor + Species/nde + Species/Nm +
                 (1|Site), data = resp[["BD"]])

summary(mod5.2.1)
Anova(mod5.2, test.statistic = "F")
summary(mod5.2)


#no interactions
mod5.4<-lmer(BD_trans ~ Species + year + pr + tem + ev + pa + 
               wt + cov + nor + nde + Nm + (1|Site), 
             data = resp[["BD"]])

plot(mod5.4)
summary(mod5.4)
Anova(mod5.4, test.statistic = "F")

#BD season null model
mod5n<-lmer((BD_trans)~1 + (1|Site),
            data = resp[["BD"]])
plot(mod5n)
summary(mod5n)


#______Excluding outliers from model with no interactions####
outliersBDseas<-resp[["BD"]][which(residuals(mod5.2)<(-1.0)),]

data_long5<-resp[["BD"]][which(residuals(mod5.2)>(-1.0)),]

modx5.1<-lmer(BD_trans ~ Species/year + Species/pr + Species/tem + Species/ev + Species/pa + 
                Species/wt + Species/cov + Species/nor + Species/nde + Species/Nm +
                (1|Site), 
              data = data_long5)

mod5.1<-lmer(BD_trans ~ Species*year + Species*pr + Species*tem + Species*ev + Species*pa + 
               Species*wt + Species*cov + Species*nor + Species*nde + Species*Nm +
                (1|Site), 
              data = data_long5)

plot(mod5.1)
summary(modx5.1)
Anova(mod5.1, test.statistic = "F")

#no inrteractions
modx5<-lmer(BD_trans ~ Species + year + pr + tem + ev + pa + wt + cov + nor + nde + Nm +
              (1|Site), 
            data = data_long5, na.action = na.omit)
plot(modx5)
summary(modx5)
Anova(modx5, test.statistic = "F")##slut p?? outlier test


#MIXED MODEL with prec.-evap.
mod5.5<-lmer(BD_trans ~ Species*year + Species*ev_pre  + Species*tem + Species*pa_d + 
               Species*wt + Species*cov + Species*nor + Species*nde + Species*Nm+
               (1|Site),
             data = resp[["BD"]])
plot(mod5.5)

mod5.5x<-lmer(BD_trans ~ Species/year + Species/ev_pre  + Species/tem + Species/pa_d + 
                Species/wt + Species/cov + Species/nor + Species/nde + Species/Nm+
                (1|Site),
              data = resp[["BD"]])

summary(mod5.5)
Anova(mod5.5, test.statistic = "F")
summary(mod5.5x)

#Run without interactions to get coefficients
mod5.6<-lmer(BD_trans ~ Species + year + ev_pre  + tem + pa_d + wt + cov + nor+ nde + Nm+ 
               (1|Site),
             data = resp[["BD"]])
plot(mod5.6)
summary(mod5.6)
Anova(mod5.6, test.statistic = "F")
#data_long.seas2.na<- mod3.6@frame # save data w/o NAs

#___Numerical density####
mod6<-lmer(numden ~ Species*year + Species*pr + Species*tem + Species*ev + Species*pa + 
             Species*wt + Species*cov + Species*nor + Species*nde + Species*Nm +
             (1|Site), data = resp[["numden"]])
plot(mod6)
summary(l6<-powerTransform(mod6, family = "bcPower"))#bcnPower neccessary with neg. values lambda=0.1366
resp[["numden"]]$den_trans<-(resp[["numden"]]$numden*100)^(l6$lambda)

mod6.2<-lmer(den_trans ~ Species*year + Species*pr + Species*tem + Species*ev + Species*pa + 
               Species*wt + Species*cov + Species*nor + Species*nde + Species*Nm +
               (1|Site), data = resp[["numden"]])
plot(mod6.2)

mod6.2.1<-lmer(den_trans ~ Species/year + Species/pr + Species/tem + Species/ev + Species/pa + 
                 Species/wt + Species/cov + Species/nor + Species/nde + Species/Nm +
                 (1|Site), resp[["numden"]])

summary(mod6.2.1)
Anova(mod6.2, test.statistic = "F")
summary(mod6.2)


#no interactions
mod6.4<-lmer(den_trans ~ Species + year + pr + tem + ev + pa + 
               wt + cov + nor + nde + Nm + (1|Site), 
             data = resp[["numden"]])

plot(mod6.4)
summary(mod6.4)
Anova(mod6.4, test.statistic = "F")

#den season null model
mod6n<-lmer((den_trans)~1 + (1|Site),
            data = resp[["numden"]])
plot(mod6n)
summary(mod6n)


#______Excluding outliers from model with no interactions####
outliersdenseas<-resp[["numden"]][which(residuals(mod6.2)<(-0.55)),]

data_long6<-data_long.na6[which(residuals(mod6.2)>(-0.55)),]

modx6.1<-lmer(den_trans ~ Species/year + Species/pr + Species/tem + Species/ev + Species/pa + 
                Species/wt + Species/cov + Species/nor + Species/nde + Species/Nm +
                (1|Site), 
              data = data_long6)

mod6.1<-lmer(den_trans ~ Species*year + Species*pr + Species*tem + Species*ev + Species*pa + 
               Species*wt + Species*cov + Species*nor + Species*nde + Species*Nm +
               (1|Site), 
             data = data_long6)

plot(mod6.1)
summary(modx6.1)
Anova(mod6.1, test.statistic = "F")

#no inrteractions
modx6<-lmer(den_trans ~ Species + year + pr + tem + ev + pa + wt + cov + nor + nde + Nm +
              (1|Site), 
            data = data_long6)
plot(modx6)
summary(modx6)
Anova(modx6, test.statistic = "F")##end outlier test


#MIXED MODEL with prec.-evap.
mod6.5<-lmer(den_trans ~ Species*year + Species*ev_pre  + Species*tem + Species*pa_d + 
               Species*wt + Species*cov + Species*nor + Species*nde + Species*Nm+
               (1|Site),
             data = resp[["numden"]])
plot(mod6.5)

mod6.5x<-lmer(den_trans ~ Species/year + Species/ev_pre  + Species/tem + Species/pa_d + 
                Species/wt + Species/cov + Species/nor + Species/nde + Species/Nm+
                (1|Site),
              data = resp[["numden"]])

summary(mod6.5)
Anova(mod6.5, test.statistic = "F")
summary(mod6.5x)

#Run without interactions to get coefficients
mod6.6<-lmer(den_trans ~ Species + year + ev_pre  + tem + pa_d + wt + cov + nor+ nde + Nm+ 
               (1|Site),
             data = resp[["numden"]])
plot(mod6.6)
summary(mod6.6)
Anova(mod6.6, test.statistic = "F")

# Variances, R2 tables ####

# get R2 values for (1) between, and within site, (2) marginal (total variation explained by fixed effects)
# and conditional (total variation fixed + random effects)

#___LI season####
m1 <- lmer(LI_trans ~ Species+year + Species*pr + Species*tem + Species*ev + Species*pa + 
             Species*wt + Species*cov + Species*nor + Species*ndeL + Species*Nm +
             (year|Site), data = resp[["LI"]])
m2 <- lmer(LI_trans ~ 1 +
             (year|Site), data = resp[["LI"]])
r2.fnc(m1, m2)

#__Main effects only
m1<-lmer(LI_trans ~ Species+year+pr  + tem +ev +pa + 
           wt + cov + nor + ndeL + Nm+
           (year|Site), data = resp[["LI"]])
m2 <- lmer(LI_trans ~ 1 +
             (year|Site), data = resp[["LI"]])
r2.fnc(m1, m2)

# only PAR model
m1 <-lmer(LI_trans ~ Species+year + Species*pa +
                  (year|Site), data = resp[["LI"]])
m2 <-lmer(LI_trans ~ Species+year + 
            (year|Site), data = resp[["LI"]])
r2.fnc(m1, m2)

#___LI day####
m1<-lmer(LId_trans ~ Species+year+Species*pr_d  + Species*tem + Species*ev_d + Species*pa_d + 
               Species*wt + Species*cov + 
               Species*nor + Species*ndeL + Species*Nm+
               (year|Site), data = resp[["LI_d"]])
m2 <- lmer(LId_trans ~ 1 +
             (year|Site), data = resp[["LI_d"]])
r2.fnc(m1, m2)

#__Main effects only
m1<-lmer(LId_trans ~ Species+year+pr_d  + tem +ev_d +pa_d + 
          wt + cov + nor + ndeL + Nm+
           (year|Site), data = resp[["LI_d"]])
m2 <- lmer(LId_trans ~ 1 +
             (year|Site), data = resp[["LI_d"]])
r2.fnc(m1, m2)

#___NPP season####
m1 <-lmer(prod_trans ~ Species+year + Species*pr + Species*tem + Species*ev + Species*pa + 
               Species*wt + Species*cov + Species*nor + Species*ndeL + Species*Nm +
               (year|Site), data = resp[["prod"]][-(which(resp[["prod"]]$orgOrder ==177)),])
m2 <-lmer(prod_trans ~ 1 +
              (year|Site), data = resp[["prod"]][-(which(resp[["prod"]]$orgOrder ==177)),])
r2.fnc(NPP_E_P4, m2)

#__Main effects only
m1<-lmer(prod_trans ~ Species+year + pr + tem + ev + pa + 
           wt + cov + nor + ndeL + Nm +
           (year|Site), data = resp[["prod"]][-(which(resp[["prod"]]$orgOrder ==177)),])
m2 <- lmer(prod_trans ~ 1 +
             (year|Site), data = resp[["prod"]][-(which(resp[["prod"]]$orgOrder ==177)),])
r2.fnc(m1, m2)

#___NPP day####
m1 <-lmer(prod.d_trans ~ Species+year+Species*pr_d  + Species*tem + Species*ev_d + 
            Species*pa_d + Species*wt + Species*cov + 
            Species*nor + Species*ndeL + Species*Nm+
            (year|Site), 
          data= resp[["prod_d"]][-(which(resp[["prod_d"]]$orgOrder ==177)),])
m2 <-lmer(prod.d_trans ~ 1 +
            (year|Site), data = resp[["prod_d"]][-(which(resp[["prod_d"]]$orgOrder ==177)),])
r2.fnc(m1, m2)

#__Main effects only
m1 <-lmer(prod.d_trans ~ Species+year+pr_d  + tem + ev_d + 
            pa_d + wt + cov + 
            nor + ndeL + Nm+
            (year|Site), 
          data= resp[["prod_d"]][-(which(resp[["prod_d"]]$orgOrder ==177)),])
m2 <-lmer(prod.d_trans ~ 1 +
            (year|Site), data = resp[["prod_d"]][-(which(resp[["prod_d"]]$orgOrder ==177)),])
r2.fnc(m1, m2)

#__R2 CIs####
# run resampling R2 function to get CIs of R2 
# set sims=1000, this takes a while though!

ee <- r2.resamp(mod4.2, # m1 is the model of interest
                sims = 1000,  # number of simulations
                only.whole=FALSE) # TRUE=only whole model CIs, if FALSE also each variable(SLOW!!!)

# get Table with CI and median R2 values
ee.vec <- sapply(ee, function (x) apply(x,1:2, median))
#ee.vec[is.na(ee.vec)] <- 0
ee.ci <- apply(ee.vec,1, function (x) paste(round(quantile(x, c(0.05,0.50,0.95), na.rm =T), digits=3), 
                                            collapse = ",", sep = "") )
out = data.frame(row.names = as.character(rownames(ee[[1]])), 
                 b.site = ee.ci[1:nrow(ee[[1]])], 
                 w.site = ee.ci[(nrow(ee[[1]])+1):(2*nrow(ee[[1]]))], 
                 marginal = ee.ci[(2*nrow(ee[[1]])+1):(3*nrow(ee[[1]]))],
                 conditional = ee.ci[(3*nrow(ee[[1]])+1):(4*nrow(ee[[1]]))])
out



# Plot response Fig 5-6####
library(interactions)
library(jtools)

#__NPP####
toplot <- plot.model.fnc(mod=mod4.2,           # model name          
           back = TRUE,      # if you want to backtransform the response (FALSE if not) 
           x.scale.orig = TRUE,   # should x-axis (predictor) be on the original scale?
           orig.data = resp[["prod"]],  # orignal data frame used in fitting the model resp[["prod"]][-(which(resp[["prod"]]$orgOrder ==177)),] or resp[["LI"]]
           lam = l4$lambda,       # give transformation lambda value (lam=) for backtransformation, 
           # does not account for the *100 for daily values
           shift = min.val,  #min.val value given to get positive responses. If not set as 0.
           mult = 1,       # if response was multiplied for larger values. If not set as 1.
           ylab = expression(NPP ~ (g ~m^{-2} ~yr^{-1})), # set y-axis label 
           vars.plot = c("Species", "year",  # you always need these # IF all variables in model, set vars=NULL 
                         "pr", "tem", "ev", "cov", "nor")) # NPP: "pr", "tem", "ev", "cov", "nor"

png("npp.png", width=18, height=18, units="cm", res=300)
grid.arrange(grobs=toplot, ncol=2) # plot the plots
dev.off()

#__LI####
toplot <- plot.model.fnc(mod=mod3.2,           # model name          
                         back = TRUE,      # if you want to backtransform the response (FALSE if not) 
                         x.scale.orig = TRUE,   # should x-axis (predictor) be on the original scale?
                         orig.data = resp[["LI"]],  # orignal data frame used in fitting the model resp[["prod"]][-(which(resp[["prod"]]$orgOrder ==177)),] or resp[["LI"]]
                         lam = l3$lambda,       # give transformation lambda value (lam=) for backtransformation, 
                         # does not account for the *100 for daily values
                         shift = min.val,  #min.val value given to get positive responses. If not set as 0.
                         mult = 1,       # if response was multiplied for larger values. If not set as 1.
                         ylab = expression(LI ~ (mm ~yr^{-1})), # set y-axis label LI,  LI ~ (mm ~yr^{-1})
                         vars.plot = c("Species", "year", "cov",  # you always need these # IF all variables in model, set vars=NULL 
                                       "pr", "tem", "pa", "cov", "nor", "Nm")) # LI: "pr", "tem", "ev", "cov", "nor"

png("li.png", width=18, height=18, units="cm", res=300)
grid.arrange(grobs=toplot, ncol=2) # plot the plots
dev.off()

# Plot Biome and sites####
# Need to add worldclim variables to data file
#The biome plot using worldclim data for our sites
#names(dat3)
#tem per month 167:178
#prec per month 179:190
#calculate rowMeans, so av. per year
temp_year<-rowMeans(subset(dat, select=167:178),na.rm = TRUE)
prec_year<-(rowSums(subset(dat, select=179:190),na.rm = TRUE))/10
temp_year_site <- aggregate(temp_year ~ site.verified, dat, na.rm=TRUE,na.action="na.pass", FUN=mean)
prec_year_site<- aggregate(prec_year ~ site.verified, dat, na.rm=TRUE,na.action="na.pass", FUN=mean)
year_site<-merge(prec_year_site, temp_year_site)
#plot(temp_year_site$temp_year,prec_year_site$prec_year/10, ylim=c(16,440), xlim=c(-17, 33))


# install.packages("devtools")
# you need the 'plotbiomes' package from github
# devtools::install_github("valentinitnelav/plotbiomes")
library(ggplot2)
library(plotbiomes)
whittaker_base_plot() +
  geom_point(data = year_site,
             aes(x = temp_year,
                 y = prec_year),
             size   = 3,
             shape  = 21,
             stroke = 1,
             alpha  = 0.5) +
  theme_bw()

#Site level data and tables ####
# aggregate data frame by returning means for numeric variables
# means based on site and species

# only include bulk density data from patches where there are NPP data as well.
dat$BD13 <- ifelse(is.na(dat$Prod13), NA,  dat$BD13)
dat$BD14 <- ifelse(is.na(dat$Prod14), NA,  dat$BD14)


Avs <- aggregate(cbind(Beging_season_vasc_2013, End_season_vasc_2013, HWT_begin_season_2013, HWT_end_season_2013, Beging_season_vasc_2014,End_season_vasc_2014,HWT_begin_season_2014,HWT_end_season_2014,
                       Patch_coord_lat, Patch_coord_lon, days2013, days2014, LI13, LI14, Prod13, Prod14, LI_per_day2013, LI_per_day2014, prod_per_day2013, prod_per_day2014, BD13, BD14,
                       cover13, cover14, HWT13, HWT14, N_per, C_per, CN, P_per, ndep_Lam13,
                       temp13, temp14, evap13, evap14, ev_d13, ev_d14, precip13, precip14, prec_d13, prec_d14, norain13.mean, norain14.mean, norain13.max, norain14.max, par13, par14, par_d13, par_d14,
                       Patch_coord_lat, googlemap_elevation_masl) ~ Site + Species, dat,na.rm=TRUE, na.action="na.pass", FUN=mean)


#Average per site
Avs_site <- aggregate(cbind(Beging_season_vasc_2013, End_season_vasc_2013, HWT_begin_season_2013, HWT_end_season_2013, 
                            Beging_season_vasc_2014,End_season_vasc_2014,
                            HWT_begin_season_2014,HWT_end_season_2014,
                            cover13, cover14, HWT13, HWT14,Patch_coord_lat, Patch_coord_lon, 
                            ndep_Lam13, temp13, temp14, 
                            evap13, evap14, precip13, precip14, norain13.mean, norain14.mean, par13, par14,
                            googlemap_elevation_masl, LI13, LI14, Prod13, Prod14) ~ 
                      Site, dat, na.rm=TRUE,na.action="na.pass", FUN=mean)


#Average per site, no matter which species (for table S1,S2)
Avs_site2 <- aggregate(cbind(ndep_Lam13, googlemap_elevation_masl) ~ site.verified, dat, na.rm=TRUE,na.action="na.pass", FUN=mean)
#library(openxlsx)
#write.xlsx(Avs_site2,"Avs_site.xlsx") 

# Save site averages based on data long
av.long<-aggregate(cbind(LI_d, prod_d, LI, prod, prec_d, temp, evap_d, par_d,
                         hwt, cover, norain, ndep_Lamarque2013, Nmean) ~ 
                     Site, data_long, na.rm=TRUE,na.action="na.pass", FUN=mean)

# Correlations####
#_ndep and n tissue and ndep and vasc. plant cover####
cor.test(av.long$ndep_Lamarque2013,av.long$Nmean)#0.35, df 97
cor.test(av.long$ndep_Lamarque2013,av.long$cover)#0.25, df 95

Avs_site$mean_cover <- apply(Avs_site[, c("cover13","cover14")],1, mean, na.rm=T)#average cover between the years
cor.test(Avs_site$ndep_Lam13, Avs_site$mean_cover)#0.23, p=0.045, df=74
cor.test(Avs_site$ndep_Lam13, Avs_site$Nmean)

#_LI and NPP correlated?####
with(data_long[data_long$Species=="S.fuscum",], cor.test(prod, LI))#0.59, df 574
with(data_long[data_long$Species=="S.magellanicum",], cor.test(prod, LI))#0.57, df 598

#_Beginning and end season####
cor.test(data_long$HWT_begin_season_2013,data_long$HWT_end_season_2013)#cor=0.75, df=1106 
cor.test(data_long$HWT_begin_season_2014,data_long$HWT_end_season_2014)#=0.70, df 1132
cor.test(data_long$Beging_season_vasc_2013,data_long$End_season_vasc_2013)#0.92, df 1036
cor.test(data_long$Beging_season_vasc_2014,data_long$End_season_vasc_2014)#0.86, df 896

#On site level.
cor.test(Avs_site$HWT_begin_season_2013, Avs_site$HWT_end_season_2013)#r=0.79, df 83
cor.test(Avs_site$HWT_begin_season_2014, Avs_site$HWT_end_season_2014)#r=0.79, df 85
cor.test(Avs_site$Beging_season_vasc_2013, Avs_site$End_season_vasc_2013)#r=0.92, df 76
cor.test(Avs_site$Beging_season_vasc_2014, Avs_site$End_season_vasc_2014)#r=0.88, df 64

#_Within species correlations####
Av.fu<-Avs[which(Avs$Species== "S.fuscum"),]#S. fuscum df
levels(as.factor(Av.fu$Site))#85 Sites
cor.test(Av.fu$HWT13, Av.fu$HWT14)#r=0,74, df 67
cor.test(Av.fu$cover13, Av.fu$cover14)#r=0,89, df 61
cor.test(Av.fu$LI13, Av.fu$LI14)#r=0,68, df=68, p<0.0001
cor.test(Av.fu$Prod13, Av.fu$Prod14)#r=0,58, df=65, p<0.0001
cor.test(Av.fu$LI13, Av.fu$Prod13)#r=0,51, df=74, p<0.0001
cor.test(Av.fu$LI14, Av.fu$Prod14)#r=0,55, df=73, p<0.0001


Av.mg<-Avs[which(Avs$Species=='S.magellanicum'),]#S. magellanicum df
levels(sitmg<-as.factor(Av.mg$Site))#91
cor.test(Av.mg$HWT13, Av.mg$HWT14)#r=0,78, df 78
cor.test(Av.mg$cover13, Av.mg$cover14)#r=0,85, df 68
cor.test(Av.mg$LI13, Av.mg$LI14)#r=0,69, df=77, p<0.0001
cor.test(Av.mg$Prod13, Av.mg$Prod14)#r=0,48, df=73, p<0.0001
cor.test(Av.mg$LI13, Av.mg$Prod13)#r=0,56, df=79, p<0.0001
cor.test(Av.mg$LI14, Av.mg$Prod14)#r=0,57, df=82, p<0.0001


#_corr plots####
library(corrplot)
M <- cor(data.frame(Avs[c("days2013", "days2014", "LI13", "LI14", "Prod13", "Prod14", "LI_per_day2013", "LI_per_day2014",   
                                "prod_per_day2013", "prod_per_day2014", "BD13", "BD14", "cover13",  "cover14", "HWT13", "HWT14",  "N_per",  "C_per", "CN",
                                "P_per",  "ndep_Lam13", "temp13", "temp14", "evap13", "evap14", "ev_d13", "ev_d14", "precip13", "precip14", "prec_d13", "prec_d14", "norain13.mean",   
                                 "norain14.mean", "norain13.max", "norain14.max", "par13", "par14", "par_d13", "par_d14", "Patch_coord_lat")] ), use = "complete.obs")
par(mfrow=c(1,1))
corrplot(M, method = "number", number.cex = 0.7)

#Split per year and day/season to see better
M_per_d13 <- cor(data.frame(Avs[c("LI_per_day2013", "prod_per_day2013", "BD13", "cover13",  "HWT13", "N_per", "C_per", "CN",
                          "P_per",  "ndep_Lam13", "temp13", "ev_d13", "prec_d13", "norain13.mean",   
                           "norain13.max", "par_d13", "Patch_coord_lat")] ), use = "complete.obs")
M_per_d14 <- cor(data.frame(Avs[c( "LI_per_day2014", "prod_per_day2014", "BD14",   "cover14",  "HWT14",  "N_per",  "C_per", "CN",
                                  "P_per",  "ndep_Lam13",  "temp14",  "ev_d14",  "prec_d14", "norain14.mean",  "norain14.max", 
                                  "par_d14", "Patch_coord_lat")] ), use = "complete.obs")
M_season13 <- cor(data.frame(Avs[c("LI13", "Prod13",  "BD13", "cover13", "HWT13", "N_per",  "C_per", "CN",
                          "P_per",  "ndep_Lam13", "temp13", "evap13", "precip13", "norain13.mean",   
                           "norain13.max",  "par13", "Patch_coord_lat")] ), use = "complete.obs")
M_season14 <- cor(data.frame(Avs[c("days2014", "LI14", "Prod14", "BD14",  "cover14", "HWT14",  "N_per",  "C_per", "CN",
                                   "P_per",  "ndep_Lam13",  "temp14", "evap14", "precip14", "norain14.mean", "norain14.max", "par14", "Patch_coord_lat")] ), use = "complete.obs")

par(mfrow=c(1,1))
corrplot(M, method = "number", number.cex = 0.7)
corrplot(M_per_d13, method = "number", number.cex = 0.7)
corrplot(M_per_d14, method = "number", number.cex = 0.7)
corrplot(M_season13, method = "number", number.cex = 0.7)
corrplot(M_season14, method = "number", number.cex = 0.7)

#Pairs for each year
#2013
pairs(data.frame(Avs2[c("LI_per_day2013", "prod_per_day2013", "BD13",
                       "BD14", "temp13", "evap13",  "precip13",
                       "norain13.mean", "norain13.max", "par13", "ndep_Lam13", "Beging_season_vasc_2013","End_season_vasc_2013",
                       "HWT_begin_season_2013",  "HWT_end_season_2013", "CN", "P_per", "Patch_coord_lat")] ), main= "2013 excluding China & Japan")
#2014
pairs(data.frame(Avs[c("LI_per_day2014", "prod_per_day2014", 
                       "BD14",  "temp14",  "evap14",  "precip14",
                        "norain14.mean", "norain14.max", "par14", "ndep_Lam13", 
                       "HWT_begin_season_2014", "HWT_end_season_2014", "CN", "P_per", "Patch_coord_lat")] ),main= "2014 excluding China & Japan")

###Bar plots Figure 4####
# using datafram Avs

#function for se:
std<-function(x) sd(x, na.rm=T)/sqrt(sum(!is.na(x)))

LI13average<-tapply(Avs$LI13,list(Avs$Species),mean,na.rm=T)
LI13max<-tapply(Avs$LI13,list(Avs$Species),max, na.rm=T)
LIsampN<- cbind(LI13_N = colSums(table(Avs$LI13, Avs$Species)), LI14_N = colSums(table(Avs$LI14, Avs$Species)))

LI13min<-tapply(Avs$LI13,list(Avs$Species),min, na.rm=T)
LI13SE<-tapply(Avs$LI13,list(Avs$Species),std)
LI14average<-tapply(Avs$LI14,list(Avs$Species),mean,na.rm=T)
LI14max<-tapply(Avs$LI14,list(Avs$Species),max, na.rm=T)
LI14min<-tapply(Avs$LI14,list(Avs$Species),min, na.rm=T)
LI14SE<-tapply(Avs$LI14,list(Avs$Species),std)
LItable<-cbind(LI13average, LI13SE, LI13min, LI13max, LI14average, LI14SE, LI14min, LI14max, LIsampN)

BD13<-Avs$BD13
BD14<-Avs$BD14
BDsampN<- cbind(BD13_N = colSums(table(Avs$BD13, Avs$Species)), BD14_N = colSums(table(Avs$BD14, Avs$Species)))

BD13average<-tapply(BD13,list(Avs$Species),mean,na.rm=T)
BD13max<-tapply(BD13,list(Avs$Species),max, na.rm=T)
BD13min<-tapply(BD13,list(Avs$Species),min, na.rm=T)
BD13SE<-tapply(BD13,list(Avs$Species),std)
BD14average<-tapply(BD14,list(Avs$Species),mean,na.rm=T)
BD14max<-tapply(BD14,list(Avs$Species),max, na.rm=T)
BD14min<-tapply(BD14,list(Avs$Species),min, na.rm=T)
BD14SE<-tapply(BD14,list(Avs$Species),std)
BDtable<-cbind(BD13average, BD13SE, BD13min, BD13max, BD14average, BD14SE, BD14min, BD14max, BDsampN)
#cbind(Prod14, Avs$Species)

Prod13<-Avs$Prod13
Prod14<-Avs$Prod14
ProdSampN<- cbind(Prod13_N =colSums(table(Avs$Prod13, Avs$Species)), Prod14_N = colSums(table(Avs$Prod14, Avs$Species)))

Prod13average<-tapply(Prod13,list(Avs$Species),mean,na.rm=T)
Prod13max<-tapply(Prod13,list(Avs$Species),max, na.rm=T)
Prod13min<-tapply(Prod13,list(Avs$Species),min, na.rm=T)
Prod13SE<-tapply(Prod13,list(Avs$Species),std)
Prod14average<-tapply(Prod14,list(Avs$Species),mean,na.rm=T)
Prod14max<-tapply(Prod14,list(Avs$Species),max, na.rm=T)
Prod14min<-tapply(Prod14,list(Avs$Species),min, na.rm=T)
Prod14SE<-tapply(Prod14,list(Avs$Species),std)
Prodtable<-cbind(Prod13average, Prod13SE, Prod13min, Prod13max, Prod14average, Prod14SE, Prod14min, Prod14max, ProdSampN)

#__Bar plot LI season####
par(mfrow=c(3,1), mar = c(3,5,2,2) + 0.8, cex=0.9)#bottom, left, top and right margins
meansLI<-c(LItable[1,1],LItable[1,5],LItable[2,1],LItable[2,5])
ses<-c(LItable[1,2],LItable[1,6],LItable[2,2],LItable[2,6])
stapel<-barplot(c(meansLI),beside=T, ylim = c(0,25), col=c("lightcyan3","lightcoral","lightcyan3","lightcoral"),
                ylab=expression("LI"  ~ (mm ~yr^{-1})), las=1.1,cex.axis = 1.5, cex.lab=1.6, tck=-0.02, space=c(1,0,0.6,0))
arrows(stapel, meansLI, y1=c(meansLI-ses,meansLI+ses), angle=90, length=0)

#X axis labels
lab <- c(expression(paste(italic("S. fuscum"))),expression(paste(italic("S. magellanicum"))))
ticks<-c(2,4.8)
text(x = ticks, y=-1,  labels = lab, srt = 0, xpd = NA, pos=1, cex=1.5) #adj=c(0.9,-1),offset=0.04,
text(1.55, 2, "n = 80", cex=1.4)
text(2.55, 2, "n = 74", cex=1.4)
text(4.1, 2, "n = 86", cex=1.4)
text(5.1, 2, "n = 83", cex=1.4)
legend(1,30,horiz=T,xpd=T,c("2013", "2014"),title="Year", bty="n",fill=c("lightcyan3","lightcoral"),cex=1.5, x.intersp=0.3, text.width = 0.75 )
text(1, 30,xpd=T, "A", cex=1.4, adj=7)

#__Bar plot BD ####
meansBD<-c(BDtable[1,1],BDtable[1,5],BDtable[2,1],BDtable[2,5])
stapel<-barplot(c(meansBD),beside=T, ylim = c(0,20), col=c("lightcyan3","lightcoral","lightcyan3","lightcoral"),
                ylab=expression("Bulk density"  ~ (kg ~m^{-3})), cex.axis = 1.5, cex.lab=1.5, las=1, tck=-0.02, space=c(1,0,0.6,0))#

ses<-c(BDtable[1,2],BDtable[1,6],BDtable[2,2],BDtable[2,6])
arrows(stapel, meansBD, y1=c(meansBD-ses,meansBD+ses), angle=90, length=0)

#X axis labels
ticks<-c(2,4.6)
lab <- c(expression(paste(italic("S. fuscum"))),expression(paste(italic("S. magellanicum"))))
text(x = ticks, y=-1,  labels = lab, srt = 0, xpd = NA, pos=1, cex=1.5) #adj=c(0.9,-1),offset=0.04,
text(1.55, 2, "n = 77", cex=1.4)
text(2.55, 2, "n = 74", cex=1.4)
text(4.1, 2, "n = 82", cex=1.4)
text(5.1, 2, "n = 83", cex=1.4)
text(1, 23,xpd=T, "B", cex=1.4, adj=7)

#__Bar plot prod####
#par(mfrow=c(1,1), mar = c(5,4,2,3) + 0.8)#bottom, left, top and right margins
meansProd<-c(Prodtable[1,1],Prodtable[1,5],Prodtable[2,1],Prodtable[2,5])
stapel<-barplot(c(meansProd),beside=T, ylim = c(0,250), col=c("lightcyan3","lightcoral","lightcyan3","lightcoral"),
                ylab=expression("NPP"  ~ (g ~m^{-2} ~yr^{-1})), las=1,cex.axis = 1.4, cex.lab=1.5, tck=-0.02, space=c(1,0,0.6,0))
ses<-c(Prodtable[1,2],Prodtable[1,6],Prodtable[2,2],Prodtable[2,6])
arrows(stapel, meansProd, y1=c(meansProd-ses,meansProd+ses), angle=90, length=0)

#X axis labels
ticks<-c(2,4.6)
lab <- c(expression(paste(italic("S. fuscum"))),expression(paste(italic("S. magellanicum"))))
text(x = ticks, y=-1,  labels = lab, srt = 0, xpd = NA, pos=1, cex=1.5)
text(1.5, 20, "n = 77", cex=1.4)
text(2.5, 20, "n = 74", cex=1.4)
text(4.1, 20, "n = 82", cex=1.4)
text(5.1, 20, "n = 83", cex=1.4)
text(1.1, 260,xpd=T, "C", cex=1.4, adj=7)



#Tables with descriptive stats for Supplementary material, from averaged per site df####
std<-function(x) sd(x, na.rm=T)/sqrt(sum(!is.na(x)))
names(Avs)
LI13average<-tapply(Avs$LI13,list(Avs$Species),mean,na.rm=T)
LI13SE<-tapply(Avs$LI13,list(Avs$Species),std)
LI13min<-tapply(Avs$LI13,list(Avs$Species),min, na.rm=T)
LI13max<-tapply(Avs$LI13,list(Avs$Species),max, na.rm=T)
NPP13average<-tapply(Avs$Prod13,list(Avs$Species),mean,na.rm=T)
NPP13SE<-tapply(Avs$Prod13,list(Avs$Species),std)
NPP13min<-tapply(Avs$Prod13,list(Avs$Species),min, na.rm=T)
NPP13max<-tapply(Avs$Prod13,list(Avs$Species),max, na.rm=T)
Prec13average<-tapply(Avs$precip13,list(Avs$Species),mean,na.rm=T)
Prec13SE<-tapply(Avs$precip13,list(Avs$Species),std)
Prec13min<-tapply(Avs$precip13,list(Avs$Species),min, na.rm=T)
Prec13max<-tapply(Avs$precip13,list(Avs$Species),max, na.rm=T)
Temp13average<-tapply(Avs$temp13,list(Avs$Species),mean,na.rm=T)
Temp13SE<-tapply(Avs$temp13,list(Avs$Species),std)
Temp13min<-tapply(Avs$temp13,list(Avs$Species),min, na.rm=T)
Temp13max<-tapply(Avs$temp13,list(Avs$Species),max, na.rm=T)
Evap13average<-tapply(Avs$evap13,list(Avs$Species),mean,na.rm=T)
Evap13SE<-tapply(Avs$evap13,list(Avs$Species),std)
Evap13min<-tapply(Avs$evap13,list(Avs$Species),min, na.rm=T)
Evap13max<-tapply(Avs$evap13,list(Avs$Species),max, na.rm=T)
PAR13average<-tapply(Avs$par13,list(Avs$Species),mean,na.rm=T)
PAR13SE<-tapply(Avs$par13,list(Avs$Species),std)
PAR13min<-tapply(Avs$par13,list(Avs$Species),min, na.rm=T)
PAR13max<-tapply(Avs$par13,list(Avs$Species),max, na.rm=T)
HWT13average<-tapply(Avs$HWT13,list(Avs$Species),mean,na.rm=T)
HWT13SE<-tapply(Avs$HWT13,list(Avs$Species),std)
HWT13min<-tapply(Avs$HWT13,list(Avs$Species),min, na.rm=T)
HWT13max<-tapply(Avs$HWT13,list(Avs$Species),max, na.rm=T)
Cover13average<-tapply(Avs$cover13,list(Avs$Species),mean,na.rm=T)
Cover13SE<-tapply(Avs$cover13,list(Avs$Species),std)
Cover13min<-tapply(Avs$cover13,list(Avs$Species),min, na.rm=T)
Cover13max<-tapply(Avs$cover13,list(Avs$Species),max, na.rm=T)
Nor13average<-tapply(Avs$norain13.mean,list(Avs$Species),mean,na.rm=T)
Nor13SE<-tapply(Avs$norain13.mean,list(Avs$Species),std)
Nor13min<-tapply(Avs$norain13.mean,list(Avs$Species),min, na.rm=T)
Nor13max<-tapply(Avs$norain13.mean,list(Avs$Species),max, na.rm=T)

LI14average<-tapply(Avs$LI14,list(Avs$Species),mean,na.rm=T)
LI14SE<-tapply(Avs$LI14,list(Avs$Species),std)
LI14min<-tapply(Avs$LI14,list(Avs$Species),min, na.rm=T)
LI14max<-tapply(Avs$LI14,list(Avs$Species),max, na.rm=T)
NPP14average<-tapply(Avs$Prod14,list(Avs$Species),mean,na.rm=T)
NPP14SE<-tapply(Avs$Prod14,list(Avs$Species),std)
NPP14min<-tapply(Avs$Prod14,list(Avs$Species),min, na.rm=T)
NPP14max<-tapply(Avs$Prod14,list(Avs$Species),max, na.rm=T)
Prec14average<-tapply(Avs$precip14,list(Avs$Species),mean,na.rm=T)
Prec14SE<-tapply(Avs$precip14,list(Avs$Species),std)
Prec14min<-tapply(Avs$precip14,list(Avs$Species),min, na.rm=T)
Prec14max<-tapply(Avs$precip14,list(Avs$Species),max, na.rm=T)
Temp14average<-tapply(Avs$temp14,list(Avs$Species),mean,na.rm=T)
Temp14SE<-tapply(Avs$temp14,list(Avs$Species),std)
Temp14min<-tapply(Avs$temp14,list(Avs$Species),min, na.rm=T)
Temp14max<-tapply(Avs$temp14,list(Avs$Species),max, na.rm=T)
Evap14average<-tapply(Avs$evap14,list(Avs$Species),mean,na.rm=T)
Evap14SE<-tapply(Avs$evap14,list(Avs$Species),std)
Evap14min<-tapply(Avs$evap14,list(Avs$Species),min, na.rm=T)
Evap14max<-tapply(Avs$evap14,list(Avs$Species),max, na.rm=T)
PAR14average<-tapply(Avs$par14,list(Avs$Species),mean,na.rm=T)
PAR14SE<-tapply(Avs$par14,list(Avs$Species),std)
PAR14min<-tapply(Avs$par14,list(Avs$Species),min, na.rm=T)
PAR14max<-tapply(Avs$par14,list(Avs$Species),max, na.rm=T)
HWT14average<-tapply(Avs$HWT14,list(Avs$Species),mean,na.rm=T)
HWT14SE<-tapply(Avs$HWT14,list(Avs$Species),std)
HWT14min<-tapply(Avs$HWT14,list(Avs$Species),min, na.rm=T)
HWT14max<-tapply(Avs$HWT14,list(Avs$Species),max, na.rm=T)
Cover14average<-tapply(Avs$cover14,list(Avs$Species),mean,na.rm=T)
Cover14SE<-tapply(Avs$cover14,list(Avs$Species),std)
Cover14min<-tapply(Avs$cover14,list(Avs$Species),min, na.rm=T)
Cover14max<-tapply(Avs$cover14,list(Avs$Species),max, na.rm=T)
Nor14average<-tapply(Avs$norain14.mean,list(Avs$Species),mean,na.rm=T)
Nor14SE<-tapply(Avs$norain14.mean,list(Avs$Species),std)
Nor14min<-tapply(Avs$norain14.mean,list(Avs$Species),min, na.rm=T)
Nor14max<-tapply(Avs$norain14.mean,list(Avs$Species),max, na.rm=T)

Ndepaverage<-tapply(Avs$ndep_Lam13,list(Avs$Species),mean,na.rm=T)
NdepSE<-tapply(Avs$ndep_Lam13,list(Avs$Species),std)
Ndepmin<-tapply(Avs$ndep_Lam13,list(Avs$Species),min, na.rm=T)
Ndepmax<-tapply(Avs$ndep_Lam13,list(Avs$Species),max, na.rm=T)
Ntissaverage<-tapply(Avs$N_per,list(Avs$Species),mean,na.rm=T)
NtissSE<-tapply(Avs$N_per,list(Avs$Species),std)
Ntissmin<-tapply(Avs$N_per,list(Avs$Species),min, na.rm=T)
Ntissmax<-tapply(Avs$N_per,list(Avs$Species),max, na.rm=T)
Eleaverage<-tapply(Avs$googlemap_elevation_masl,list(Avs$Species),mean,na.rm=T)
EleSE<-tapply(Avs$googlemap_elevation_masl,list(Avs$Species),std)
Elemin<-tapply(Avs$googlemap_elevation_masl,list(Avs$Species),min, na.rm=T)
Elemax<-tapply(Avs$googlemap_elevation_masl,list(Avs$Species),max, na.rm=T)

tLI<-cbind(LI13average, LI13SE, LI13min, LI13max, LI14average, LI14SE, LI14min, LI14max)
tNPP<-cbind(NPP13average, NPP13SE, NPP13min, NPP13max, NPP14average, NPP14SE, NPP14min, NPP14max)
tPrec<-cbind(Prec13average, Prec13SE, Prec13min, Prec13max, Prec14average, Prec14SE, Prec14min, Prec14max)
tTemp<-cbind(Temp13average, Temp13SE, Temp13min, Temp13max, Temp14average, Temp14SE, Temp14min, Temp14max)
tEvap<-cbind(Evap13average, Evap13SE, Evap13min, Evap13max, Evap14average, Evap14SE, Evap14min, Evap14max)
tPAR<-cbind(PAR13average, PAR13SE, PAR13min, PAR13max, PAR14average, PAR14SE, PAR14min, PAR14max)
tHWT<-cbind(HWT13average, HWT13SE, HWT13min, HWT13max, HWT14average, HWT14SE, HWT14min, HWT14max)
tCover<-cbind(Cover13average, Cover13SE, Cover13min, Cover13max, Cover14average, Cover14SE, Cover14min, Cover14max)
tNOR<-cbind(Nor13average, Nor13SE, Nor13min, Nor13max, Nor14average, Nor14SE, Nor14min, Nor14max)
tNdep<-cbind(Ndepaverage, NdepSE, Ndepmin, Ndepmax)
tNtiss<-cbind(Ntissaverage, NtissSE, Ntissmin, Ntissmax)
tElevation<-cbind(Eleaverage, EleSE, Elemin, Elemax)

desc.table<-cbind(tLI, tNPP, tPrec, tTemp, tEvap, tPAR, tHWT, tCover, tNOR, tNdep, tNtiss, tElevation)
#library(openxlsx)
#write.xlsx(desc.table,"Descriptive.xlsx") 
#desc.table<-cbind(tNdep, tNtiss, tElevation)

#For how many sites do we have data?####
#Cover
length(which(Avs_site$cover13>0))# & Avs_site$cover13>0)) 90 2013,  2014, 83 
length(which(Avs_site$cover13>0 | Avs_site$cover14 >0)) #total 76 with data for both years, either 97

#HWT
length(which(Avs_site$HWT13>0))# 94 2013, 90 2014
length(which(Avs_site$HWT13>0 & Avs_site$HWT14>0)) #total 85 with data for both years, 97 either

#Data for both spring and fall?
#Cover 
length(which(Avs_site$Beging_season_vasc_2013>0 & Avs_site$End_season_vasc_2013>0)) #total 78 with data for both spring and fall in 2013
91-78#=13 sites with data that only has data for one period 
length(which(Avs_site$End_season_vasc_2013>0))#87 begining, 82 end season

length(which(Avs_site$Beging_season_vasc_2014>0 & Avs_site$End_season_vasc_2014>0)) #total 66 with data for both spring and fall in 2013
83-66#= 17 sites with data for only one period 
length(which(Avs_site$Beging_season_vasc_2014>0))#76 begining, 73 end season

#HWT
length(which(Avs_site$HWT_begin_season_2013>0 & Avs_site$HWT_end_season_2013>0)) #total 85 with data for both spring and fall in 2013
95-85#=10 sites that only has data for one period 
length(which(Avs_site$HWT_begin_season_2013>0))#91 begining, 89 end season

length(which(Avs_site$HWT_begin_season_2014>0 & Avs_site$HWT_end_season_2014>0)) #total 87 with data for both spring and fall in 2013
90-87#= 3 sites with data for only one period 
length(which(Avs_site$HWT_end_season_2014>0))#89 begining, 88 end season

#Overall Avs_site : averages between species per site calculated 
mean(Avs_site$LI13, na.rm=T)#16.95
mean(Avs_site$LI14, na.rm=T)#18.12
mean(Avs_site$Prod13, na.rm=T)#194.70
mean(Avs_site$Prod14, na.rm=T)#182.10
std<-function(x) sd(x, na.rm=T)/sqrt(sum(!is.na(x)))
std(Avs_site$LI13)#1.09
std(Avs_site$LI14)#1.19
std(Avs_site$Prod13)#12.13
std(Avs_site$Prod14)#11.12

