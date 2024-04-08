## Simulation Parameters using 6 Covariates
### Scenario with 50\% Vaccine Coverage 
## March 6 2024

# Set library path
.libPaths(c("/home/users/landrew2/R_lib"))
library(readr)
library(dplyr)
#file_path <- "/home/users/landrew2/Dissertation/Project1/"



#### Clean Data ######
## Remove those with missing race/ethnicity

# 
# 
# simtest0<- covid.rct.df %>% filter(!is.na(RACEc)&!is.na(ETHNICc))%>%
#   select(USUBJID,SUBJID,AGE,SEX,REGION,RISKGR1,POC,RACE,RACEc,
#          ETHNICc, AGEGR1,TRT01PN,CNSR,CASE, CALTIME, VACCINE) %>%
#   mutate(
#     # Combine Hawaiian and other race categories
#     RACE6c = factor(ifelse(RACEc=="NATIVE HAWAIIAN OR OTHER PACIFIC ISLANDER","OTHER",RACEc)),
#     RACE6c=relevel(RACE6c, ref="WHITE"),
#     ETHNICc = relevel(factor(ETHNICc), ref="NOT HISPANIC OR LATINO"),
#     RISKGR1 = relevel(factor(RISKGR1), ref="Not At Risk"),
#     CALTIMEs2 =  I(CALTIME>=90)*(CALTIME-90), #ifelse(CALTIME>=90, CALTIME-90,0)
#     CALTIMEs3 = I(CALTIME>=135)*(CALTIME-135),# ifelse(CALTIME>=135, CALTIME-135,0)) 
#     CALTIMEs4 = I(CALTIME>=180)*(CALTIME-180),
#     FEMALE = ifelse(SEX=="F",1,0)) 
# 
# simtest <- simtest0 %>% select(SUBJID, AGE, FEMALE, RISKGR1, 
#                                CALTIME, CALTIMEs2,CALTIMEs3,CALTIMEs4)

#write.csv(simtest, "~/Desktop/UW Stuff/IS TND/Code/simtest.csv")
#simtest<-read_csv("~/Desktop/UW Stuff/IS TND/Code/simtest.csv")


simtest0<-read_csv(paste0(file_path,"simtest.csv"))
#simtest0 <- simtest
simtest <- simtest0 %>% mutate(
  RISKGR10 = RISKGR1,
  # relevel Comorbidities variable
  RISKGR1 = relevel(factor(RISKGR1), ref="Not At Risk")
)
table(simtest$RISKGR10)
table(simtest$RISKGR1)

dim(simtest)


####### Vaccination A Parameters #########

#a_formula <- "~FEMALE +RISKGR1+ CALTIME"
a_formula <- "~ FEMALE + RISKGR1+ FEMALE:RISKGR1+ CALTIME+CALTIMEs2+CALTIMEs3"

### Sex
# female prob of vacc 52%, male 48%
#fem_a_odds <- 52/48
#male_a_odds <- 48/52
sex_a_or <- 3 #2 #3 #male_a_odds/fem_a_odds
sex_a_beta <- log(sex_a_or)
names(sex_a_beta) <- "FEMALE"

### Comorbidities
#comorb_a_odds <- 60/40
#nocomorb_a_odds <- 45/55
comorb_a_or <- 4#2 #3 #4 #comorb_a_odds/nocomorb_a_odds
comorb_a_beta <- log(comorb_a_or)
names(comorb_a_beta) <- "RISKGR1At Risk"

### Sex and Comorbidity Interaction
sexxcomorb_a_beta0 <- 0
names(sexxcomorb_a_beta0) <-"FEMALE:RISKGR1At Risk"
sexxcomorb_a_ror <-  0.25 #3 # 4#3 
sexxcomorb_a_beta <-  log(sexxcomorb_a_ror) # log of ratio of two ORs 
names(sexxcomorb_a_beta) <-"FEMALE:RISKGR1At Risk"

### Calendar Time
# at time 0 (September 1), 40% vacc prob, so 0.67 vacc odds
sept_a_prob <- 0.4
sept_a_odds <- sept_a_prob/(1-sept_a_prob)
# at time 90, beginning of December
dec_a_prob <- 0.5
dec_a_odds <- dec_a_prob/(1-dec_a_prob)
septdec_a_or <- dec_a_odds/sept_a_odds
septdec_a_beta<-log(septdec_a_or)/90
# at time 135 (mid january)
jan_a_prob <- 0.6
jan_a_odds <- jan_a_prob/(1-jan_a_prob)
decjan_a_or <- jan_a_odds/dec_a_odds
decjan_a_beta0 <- 0
decjan_a_beta<- log(decjan_a_or)/45-septdec_a_beta
# at time 210 (mid March), 60% vacc prob
mar_a_prob <- 0.30  
mar_a_odds <- mar_a_prob/(1-mar_a_prob)
# log (vacc OR comparing march 210 to sept 0)=210beta
janmar_a_or <- mar_a_odds/jan_a_odds
janmar_a_beta0 <- 0
janmar_a_beta<- log(janmar_a_or)/75-decjan_a_beta-septdec_a_beta

time_a_betas0<- c(septdec_a_beta, decjan_a_beta0, janmar_a_beta0)
names(time_a_betas0) <- c("CALTIME","CALTIMEs2","CALTIMEs3")
time_a_betas<- c(septdec_a_beta, decjan_a_beta,janmar_a_beta)
names(time_a_betas) <- c("CALTIME","CALTIMEs2","CALTIMEs3")


### Intercept (Odds of Vacc for Reference Group)
# vacc odds for noncomorbid male in September
# choose intercept to give 50% vaccination
ref_a_odds <- .20/(1-.20)    
names(ref_a_odds) <-"(Intercept)" 

#a_betas <- c(log(c(ref_a_odds,sex_a_or,comorb_a_or)),septmar_a_beta)

# Scenario 1: 0 interactions and 0 splines parameters 
(a_betas1 <- c(log(ref_a_odds), time_a_betas0,
              comorb_a_beta, sex_a_beta, sexxcomorb_a_beta0))
# Scenario 2: 0 splines parameters 
(a_betas2 <- c(log(ref_a_odds), time_a_betas0,
              comorb_a_beta, sex_a_beta,sexxcomorb_a_beta))
# Scenario 3: 
(a_betas3 <- c(log(ref_a_odds), time_a_betas,
              comorb_a_beta, sex_a_beta,sexxcomorb_a_beta))



# # Probability of Vaccination by time
# # for reference group (white, mw, nonhispanic, nocomorbid, young person)
# data.frame(x=0:210)%>%mutate(
#   y=plogis(log(ref_a_odds) +septdec_a_beta*x+
#            decjan_a_beta*I(x>=90)*(x-90)+
#              janmar_a_beta*I(x>=135)*(x-135))) %>%
#   ggplot(aes(x,y))+geom_line()+theme_bw()+
#   labs(x="Calendar Time (Sept to April)", y="Prob of Vaccination")
# data.frame(FEM=c(0,0,1,1),
#            COM = c(0,1,0,1))%>%mutate(
#              y=plogis(log(ref_a_odds) +sex_a_beta*FEM+
#                         comorb_a_beta*COM+
#                         sexxcomorb_a_beta*FEM*COM)) %>%
#   ggplot(aes(FEM,y, col=factor(COM)))+geom_line()+theme_bw()+
#   labs(x="Sex", y="Prob of Vaccination")
# 
# data.frame(FEM=c(0,0,1,1),
#            COM = c(0,1,0,1))%>%mutate(
#   y=plogis(log(ref_a_odds) +sex_a_beta*FEM+
#            comorb_a_beta*COM+
#              log(0.25)*FEM*COM)) %>%
#   ggplot(aes(COM,y, col=factor(FEM)))+geom_line()+theme_bw()+
#   labs(x="Comorbidities", y="Prob of Vaccination")
# 
# data.frame(FEM=c(0,0,1,1),
#            COM = c(0,1,0,1))%>%mutate(
#              y=plogis(log(ref_y_odds) +sex_y_beta*FEM+
#                         comorb_y_beta*COM+
#                         log(0.2)*FEM*COM)) %>%
#   ggplot(aes(COM,y, col=factor(FEM)))+geom_line()+theme_bw()+
#   labs(x="Comorbidities", y="Prob of Infection")

######## SARs-CoV-2 Infection Y Parameters ########

y_formula <- "~ A+ FEMALE + RISKGR1+ FEMALE:RISKGR1+ CALTIME+CALTIMEs2+CALTIMEs3"


### Sex (female vs male (reference))
sex_y_or <- 3 # 2 # 3 #0.9 # 0.94  # for noncomorbid
names(sex_y_or) <- "FEMALE"
sex_y_beta<- log(sex_y_or)

### Comorbidities
comorb_y_or <- 4 #1.8 # 1.06 # for male
names(comorb_y_or) <-"RISKGR1At Risk"
comorb_y_beta <- log(comorb_y_or)

### Sex and Comorbidity Interaction
sexxcomorb_y_beta0 <- 0
names(sexxcomorb_y_beta0) <-"FEMALE:RISKGR1At Risk"
sexxcomorb_y_ror <- 0.2 #3 #4 #3
sexxcomorb_y_beta <-  log(sexxcomorb_y_ror)  # log of ratio of two ORs  
names(sexxcomorb_y_beta) <-"FEMALE:RISKGR1At Risk"
# OR of infection comparing female vs male comorbid person
# exp(sex_y_beta+sexxcomorb_y_beta)
# # OR of infection comparing female vs male noncomorbid person
# exp(sex_y_beta)
# # OR of infection comparing comorbid and noncomorbid female
# exp(comorb_y_beta+sexxcomorb_y_beta)
# # OR of infection comparing comorbid and noncomorbid male
# exp(comorb_y_beta)
# # Prob of Infection for noncomorb male
# plogis(log(ref_y_odds) +comorb_y_beta*0+  sex_y_beta*0 +comorb_y_beta*sex_y_beta*0*1)
# # Prob of Infection for comorbid male
# plogis(log(ref_y_odds) +comorb_y_beta*1+  sex_y_beta*0 +comorb_y_beta*sex_y_beta*0*1)
# # Prob of Infection for noncomorb female
# plogis(log(ref_y_odds) +comorb_y_beta*0+  sex_y_beta*1 +comorb_y_beta*sex_y_beta*1*0)
# # Prob of Infection for comorbid female
# plogis(log(ref_y_odds) +comorb_y_beta*1+  sex_y_beta*1 +sexxcomorb_y_beta*1*1)


### Calendar Time
# at time 0 (Sept 1), prob of SARS-CoV-2 infection is 295000/330mil/0.6
sept_y_prob <- 0.001  #295000/330000000/0.6
sept_y_odds <- sept_y_prob/(1-sept_y_prob)
dec_y_prob <- 0.0066  #295000/330000000/0.6
dec_y_odds <- dec_y_prob/(1-dec_y_prob)
septdec_y_or <- dec_y_odds/sept_y_odds
septdec_y_beta<-log(septdec_y_or)/90
# at time 135 (mid january)
jan_y_prob <- 0.00859 #1.7 mil/330 mill/0.6
jan_y_odds <- jan_y_prob/(1-jan_y_prob)
decjan_y_or <- jan_y_odds/dec_y_odds
decjan_y_beta0 <- 0
decjan_y_beta<- log(decjan_y_or)/45-septdec_y_beta
# at time 210 (end of March, basically April)
mar_y_prob <- sept_y_prob
mar_y_odds <- sept_y_odds
janmar_y_or <- mar_y_odds/jan_y_odds
janmar_y_beta0 <- 0
janmar_y_beta <- log(janmar_y_or)/75-decjan_y_beta-septdec_y_beta

time_y_betas0<- c(septdec_y_beta, decjan_y_beta0, janmar_y_beta0)
names(time_y_betas0) <- c("CALTIME","CALTIMEs2","CALTIMEs3")
time_y_betas<- c(septdec_y_beta, decjan_y_beta, janmar_y_beta)
names(time_y_betas) <- c("CALTIME","CALTIMEs2","CALTIMEs3")

# # Prob of Infection by Time for Reference Group  
# data.frame(x=0:210)%>%mutate(
#   y=plogis(log(ref_y_odds) +septdec_y_beta*x+
#                decjan_y_beta*I(x>=90)*(x-90)+
#                  janmar_y_beta*I(x>=135)*(x-135))) %>%
#   ggplot(aes(x,y))+geom_line()+theme_bw()+
#   labs(x="Calendar Time (Sept to April)", y="Prob of SARS Infection")
# data.frame(x=0:210)%>%mutate(
#   y=plogis( -0.69805573  +0.01966002*x+
#               -0.14851899 *I(x>=90)*(x-90)+
#               0.08598985 *I(x>=135)*(x-135))) %>%
#   ggplot(aes(x,y))+geom_line()+theme_bw()+
#   labs(x="Calendar Time (Sept to April)", y="Prob of SARS Infection")

# 
# data.frame(x=0:210)%>%mutate(
#   y=plogis(log(ref_y_odds) +septdec_y_beta*x+ 
#                  0.0000000000000001*x^3)) %>%
#   ggplot(aes(x,y))+geom_line()+theme_bw()+
#   labs(x="Calendar Time (Sept to April)", y="Prob of SARS Infection")


### Vaccination
vacc_y_or <- 1 # should have 1 - (0%, 30%, 50%, 80%)
vacc_y_beta <- log(vacc_y_or)
names(vacc_y_beta) <-"A" 

### Odds of SARS-CoV-2 Infection for Reference Group
#  noncomorbid unvacc, male in MW in September
ref_y_prob <- 0.04 #0.02 #0.005
ref_y_odds <- ref_y_prob/(1-ref_y_prob)  
names(ref_y_odds) <-"(Intercept)" 

# Scenario 1: 0 interactions and 0 splines parameters 
(y_betas1 <- c(log(ref_y_odds), time_y_betas0, #vacc_y_beta,
              comorb_y_beta, sex_y_beta,sexxcomorb_y_beta0))
# Scenario 2: 0 splines parameters 
(y_betas2 <- c(log(ref_y_odds), time_y_betas0, #vacc_y_beta,
              comorb_y_beta, sex_y_beta,sexxcomorb_y_beta))
# Scenario 3: 
(y_betas3 <- c(log(ref_y_odds), time_y_betas, #vacc_y_beta,
             comorb_y_beta, sex_y_beta,sexxcomorb_y_beta))


######## Other Infections W Parameters ########

w_formula <- "~ A+ FEMALE + RISKGR1+ CALTIME+CALTIMEs2+CALTIMEs4"

### Sex (female vs male (reference))
sex_w_or <- 2 #3 #2#3 #1/0.9 # 0.94
names(sex_w_or) <- "FEMALE"
sex_w_beta<- log(sex_w_or)

### Comorbidities (same as COVID)
comorb_w_or <- 2#3# 2 #3 #1.8 # 1.06
names(comorb_w_or) <-"RISKGR1At Risk"
comorb_w_beta <- log(comorb_w_or)

### Calendar Time (second cutpoint different than SARS)
# at time 0 (Sept 1), prob of non-SARS infection
sept_w_prob <- 0.07/0.3 
sept_w_odds <- sept_w_prob/(1-sept_w_prob)
dec_w_prob <- 0.14/0.3 
dec_w_odds <- dec_w_prob/(1-dec_w_prob)
septdec_w_or <- dec_w_odds/sept_w_odds
septdec_w_beta<- log(septdec_w_or)/90
# at time 180 (~ March )
mar_w_prob <- 0.17/0.3
mar_w_odds <- mar_w_prob/(1-mar_w_prob)
decmar_w_or <- mar_w_odds/dec_w_odds
decmar_w_beta<- log(decmar_w_or)/90-septdec_w_beta
# June (just to get slope accurate)
jun_w_prob <- 0.08/0.3
jun_w_odds <- jun_w_prob/(1-jun_w_prob)
marjun_w_or <- jun_w_odds/mar_w_odds
marjun_w_beta <- log(marjun_w_or)/120-decmar_w_beta-septdec_w_beta

time_w_betas<- c(septdec_w_beta, decmar_w_beta, marjun_w_beta)
names(time_w_betas) <- c("CALTIME","CALTIMEs2","CALTIMEs4")

# # Prob of Non-SARS Infection by Time for Reference Group
# data.frame(x=0:210)%>%mutate(
#   y=plogis(log(ref_w_odds) +septdec_w_beta*x+
#                  decmar_w_beta*I(x>=90)*(x-90)+
#                  marjun_w_beta*I(x>=180)*(x-180))) %>%
#   ggplot(aes(x,y))+geom_line()+theme_bw()+
#   labs(x="Calendar Time (Sept to April)", y="Prob of non-SARS Infection")


### Vaccination 
vacc_w_or <- 1 # vaccine has no effect on other infections
vacc_w_beta <- log(vacc_w_or)
names(vacc_w_beta) <- "A" 


### Odds of Non-SARS Infection for Reference Group
# nonhisp white, noncomorbid 0 year old, unvacc male in MW in September
#ref_w_prob <- 0.1 #0.15
ref_w_odds <- 0.1 #ref_w_prob/(1-ref_w_prob)  
names(ref_w_odds) <-"(Intercept)" 
w_betas <- c(log(ref_w_odds), time_w_betas, #vacc_w_beta,
              comorb_w_beta,sex_w_beta)
#### Clinical Definition and Symptoms C #####

c_formula <- "~Y+W+ A+ Y*FEMALE +Y*RISKGR1  + W*FEMALE+ W*RISKGR1"

### intercept (odds of having symptoms when young, male, noncommorbid, no infections)
#ref_c_prob <- 0.1
ref_c_odds <- 0.1 #ref_c_prob/(1-ref_c_prob)
ref_c_beta <- log(ref_c_odds)
names(ref_c_beta) <- "(Intercept)"

### Vaccination 
vacc_c_or <- 1 # assume vaccine has no effect on symptoms
vacc_c_beta <- log(vacc_c_or)
names(vacc_c_beta) <- "A" 


### No infection 

## Sex
noinf_sex_c_beta <- 0
names(noinf_sex_c_beta) <- "FEMALE"

## Comorbidities
noinf_comorb_c_prob <- 0.13
noinf_comorb_c_odds <- noinf_comorb_c_prob/(1-noinf_comorb_c_prob)
noinf_nocomorb_c_prob <- 0.08
noinf_nocomorb_c_odds <- noinf_nocomorb_c_prob/(1-noinf_nocomorb_c_prob)
noinf_comorb_c_or <- 2 #noinf_comorb_c_odds/noinf_nocomorb_c_odds
names(noinf_comorb_c_or) <- "RISKGR1At Risk"
noinf_comorb_c_beta<- log(noinf_comorb_c_or)

### Other Infections  
## Main effect term (noncomorbid, young, male)
# 30% of common cold infections are symptomatic,
other_c_prob <- 0.3
other_c_odds <- other_c_prob/(1-other_c_prob)
noinf_c_prob <- 0.1
noinf_c_odds <- noinf_c_prob/(1-noinf_c_prob)
other_c_or <- 4 #other_c_odds/noinf_c_odds
other_c_beta <- log(other_c_or)
names(other_c_beta) <- "W"

## Sex (from challenge study)
other_fem_c_prob <- 0.96
other_fem_c_odds <- other_fem_c_prob/(1-other_fem_c_prob)
other_male_c_prob <- 0.80
other_male_c_odds <- other_male_c_prob/(1-other_male_c_prob)
other_sex_c_or <-other_fem_c_odds/other_male_c_odds
other_sex_c_beta<- log(other_sex_c_or)-noinf_sex_c_beta
names(other_sex_c_beta) <- "W:FEMALE"

## Comorbidities
other_comorb_c_or <- 1.06 # used same as COVID
names(other_comorb_c_or) <- "W:RISKGR1At Risk"
other_comorb_c_beta<- log(other_comorb_c_or)-noinf_comorb_c_beta


### SARS-CoV-2 Infection 
## Main effect term (noncomorbid,women, young)
# 60% of SARS-CoV-2 infections are symptomatic
sars_c_prob <- 0.6
sars_c_odds <- sars_c_prob/(1-sars_c_prob)
noinf_c_prob <- 0.1
noinf_c_odds <- noinf_c_prob/(1-noinf_c_prob)
sars_c_or <- sars_c_odds/noinf_c_odds
sars_c_beta <- log(sars_c_or)
names(sars_c_beta) <- "Y"

## Sex 
sars_sex_c_or <- 1/0.93
names(sars_sex_c_or) <- "Y:FEMALE"
sars_sex_c_beta<- log(sars_sex_c_or)-noinf_sex_c_beta

## Comorbidities
sars_comorb_c_or <- 1.06
names(sars_comorb_c_or) <- "Y:RISKGR1At Risk"
sars_comorb_c_beta<- log(sars_comorb_c_or)-noinf_comorb_c_beta

c_betas <- c(ref_c_beta, #vacc_c_beta, 
             noinf_sex_c_beta, noinf_comorb_c_beta, 
             other_c_beta,  other_sex_c_beta,other_comorb_c_beta, 
             sars_c_beta, sars_sex_c_beta, sars_comorb_c_beta)

##### Missing At Random Vaccination Delta Dependent on Comorbidities ######
comorb_mar_d_formula <- "~ RISKGR1"

comorb_mar_d_or <-  0.02 #0.1 #0.5
names(comorb_mar_d_or) <- "RISKGR1At Risk"
comorb_mar_d_beta<- log(comorb_mar_d_or)

ref_comorb_mar20_d_odds <- 0.38 #0.36 #0.3
names(ref_comorb_mar20_d_odds) <- "(Intercept)"
ref_comorb_mar20_d_beta<- log(ref_comorb_mar20_d_odds)

ref_comorb_mar50_d_odds <-2.1 #1.75  #1.2 #0.3
names(ref_comorb_mar50_d_odds) <- "(Intercept)"
ref_comorb_mar50_d_beta<- log(ref_comorb_mar50_d_odds)

comorb_mar20_d_betas <- c(ref_comorb_mar20_d_beta, comorb_mar_d_beta)
comorb_mar50_d_betas <- c(ref_comorb_mar50_d_beta, comorb_mar_d_beta)

# # # Prob Not at risk
# plogis(ref_comorb_mar20_d_beta)
# #
# # # Prob at risk
# plogis(ref_comorb_mar20_d_beta+ comorb_mar_d_beta)
# #
# # # marginal Prob missingness is
# plogis(ref_comorb_mar20_d_beta)*(1-0.2765)+ plogis(ref_comorb_mar20_d_beta+ comorb_mar_d_beta)*0.2765
# 
# # # Prob Not at risk
# plogis(ref_comorb_mar50_d_beta)
# #
# # # Prob at risk
# plogis(ref_comorb_mar50_d_beta+ comorb_mar_d_beta)
# #
# # # marginal missingness is
# plogis(ref_comorb_mar50_d_beta)*(1-0.2765)+ plogis(ref_comorb_mar50_d_beta+ comorb_mar_d_beta)*0.2765

#hist(raw.df$SampMargCom)
#median(raw.df$SampMargCom)

##### Missing At Random Vac/Antibodies Delta Dependent on SARS-CoV-2 Status ######
# y_mar_d_formula <- "~ Y"
# 
# y_mar_d_or <-  0.1 #0.5
# names(y_mar_d_or) <- "Y"
# y_mar_d_beta<- log(y_mar_d_or)
# 
# ref_y_mar20_d_odds <- 0.36 #0.3
# names(ref_y_mar20_d_odds) <- "(Intercept)"
# ref_y_mar20_d_beta<- log(ref_y_mar20_d_odds)
# 
# ref_y_mar50_d_odds <-1.75  #1.2 #0.3
# names(ref_y_mar50_d_odds) <- "(Intercept)"
# ref_y_mar50_d_beta<- log(ref_y_mar50_d_odds)
# 
# y_mar20_d_betas <- c(ref_y_mar20_d_beta, y_mar_d_beta)
# y_mar50_d_betas <- c(ref_y_mar50_d_beta, y_mar_d_beta)

# # # Prob Not at risk
# plogis(ref_y_mar20_d_beta) 
# # 
# # # Prob at risk
# plogis(ref_y_mar20_d_beta+ y_mar_d_beta)
# # 
# # # marginal Prob missingness is
# plogis(ref_y_mar20_d_beta)*(1-0.2765)+ plogis(ref_y_mar20_d_beta+ y_mar_d_beta)*0.2765
# 
# # # Prob Not at risk
# plogis(ref_y_mar50_d_beta) 
# # 
# # # Prob at risk
# plogis(ref_y_mar50_d_beta+ y_mar_d_beta)
# # 
# # # marginal missingness is
# plogis(ref_y_mar50_d_beta)*(1-0.646)+ plogis(ref_y_mar50_d_beta+ y_mar_d_beta)*0.646
# 
# 
# hist(raw.df$SampMargY, breaks=500)
# median(raw.df$SampMargY)
