# This is the plotting & analysis script for:
#age-1 Pacific cod comparing body condition and otoliths between Eras
#all cook bay and Era A is 2013; Era B is 2023

# Libraries
library(patchwork)
library(tidyverse)
library(lubridate)
library(ggplot2)
library(dplyr)
library(plyr)
library(mgcv)
library(rstan)
library(brms)
library(bayesplot)

# This is the data import script
# Read in the data and rename any columns that need renaming
cod1 <- read_csv("C:/Users/alask/Documents/Git/seine-data/data/cod_esass_data.csv")
head(cod1)
str(cod1)
tail(cod1)

##can see that last row are all NA. Must remove these last 3 rows
cod1<-filter(cod1, !is.na(TL))
tail(cod1)
#looks good. all data rows and col present

#### PLOTTING ####
##if need to rename columns, you would type 
#codlen <- rename(codlen,"TL" = "Length (mm)") 

#len <- len %>%
# mutate(year_fac = as.factor(year),
#         TL = as.numeric(TL),
#         month_fac = as.factor(month))

ggplot(cod1, aes(TL, fill=(Era_year))) +
  geom_histogram(alpha=0.5) +
  xlab("age-1 Pacific cod TL")+
  theme_bw()+
  theme(legend.position = "bottom")

ggplot(cod1, aes(TL, fill=(Era_year))) +
  geom_density(alpha=0.5) +
  xlab("age-1 Pacific cod TL")+
  theme_bw()+
  theme(legend.position = "bottom")

a <- lm(formula = TL ~ Era_year, data = cod1)
summary(a)


cod1 <- cod1 %>%
  mutate(Era = as.factor(Era_year))

#now begin to plot HSI and Kwet by Era
b <- lm(formula = Kwet ~ Era_year, data = cod1)
summary(b)
c <- gam(Kwet ~ s(TL, k=4) + Era, data = cod1)
summary(c)

plot_Kwet <- ggplot(data = cod1, aes(x = factor(Era_year), y = Kwet, fill = Era)) +
  geom_boxplot()+
  theme_minimal()+
  labs(x = "Era", y = "Fulton's condition factor, Kwet")+
  theme(legend.position ="bottom") 

plot(plot_Kwet)


d <- lm(formula = Kwet ~ TL + Era, data = cod1)

summary(d)

MuMIn::AICc(b,c,d)
#c is lowest, so use the GAM
plot(c)
summary(c) #no difference between ERA

d <- ggplot(data = cod1, 
            aes(x = TL,
                y = Kwet, color = Era)) +
  geom_point(size = 2.5, alpha = 0.8) + 
  theme_bw()+
  labs(y = "Fulton's K wet", x = "Total length (mm)")
d

#talk with Mike and now test Kdry, HSIwet, and HSIdry
plot_Kdry <- ggplot(data = cod1, aes(x = factor(Era_year), y = Kdry, fill = Era)) +
  geom_boxplot()+
  theme_minimal()+
  labs(x = "Era", y = "Fulton's condition factor, Kdry")+
  theme(legend.position ="bottom") 

plot(plot_Kdry)

b <- lm(formula = Kdry ~ Era_year, data = cod1)
summary(b)
c <- gam(Kdry ~ s(TL, k=4) + Era, data = cod1)
summary(c)

d <- lm(formula = Kdry ~ TL + Era, data = cod1)
MuMIn::AICc(b,c,d)
#c and d are same. best relationship is linear!
summary (d)
#no difference between era R2 = .36


e <- ggplot(data = cod1, 
            aes(x = TL,
                y = Kdry, color = Era)) +
  geom_point(size = 2.5, alpha = 0.8) + 
  theme_bw()+
  labs(y = "Fulton's K dry", x = "Total length (mm)")
e

plot_HSIwet <- ggplot(data = cod1, aes(x = factor(Era_year), y = HSI_wet, fill = Era)) +
  geom_boxplot()+
  theme_minimal()+
  labs(x = "Era", y = "Hepatosomatic Index HSI_wet")+
  theme(legend.position ="bottom") 

plot(plot_HSIwet)

f <- ggplot(data = cod1, 
            aes(x = TL,
                y = HSI_wet, color = Era)) +
  geom_point(size = 2.5, alpha = 0.8) + 
  theme_bw()+
  labs(y = "Hepatosomatic Indext HSI_wet", x = "Total length (mm)")
f


b <- lm(formula = HSI_wet ~ Era_year, data = cod1)
summary(b)
c <- gam(HSI_wet ~ s(TL, k=4) + Era, data = cod1)
summary(c)

d <- lm(formula = HSI_wet ~ TL + Era, data = cod1)
MuMIn::AICc(b,c,d)
#again, c and d are same. so use linear model. yes difference w TL and ERA
summary(d)
plot(d)


plot_HSIdry <- ggplot(data = cod1, aes(x = factor(Era_year), y = HSI_dry, fill = Era)) +
  geom_boxplot()+
  theme_minimal()+
  labs(x = "Era", y = "Hepatosomatic Index HSI_dry")+
  theme(legend.position ="bottom") 

plot(plot_HSIdry)

g <- ggplot(data = cod1, 
            aes(x = TL,
                y = HSI_dry, color = Era)) +
  geom_point(size = 2.5, alpha = 0.8) + 
  theme_bw()+
  labs(y = "Hepatosomatic Indext HSI_dry", x = "Total length (mm)")
g


b <- lm(formula = HSI_dry ~ Era_year, data = cod1)
summary(b)
c <- gam(HSI_dry ~ s(TL, k=4) + Era, data = cod1)
summary(c)

d <- lm(formula = HSI_dry ~ TL + Era, data = cod1)
e <- lm(HSI_dry ~ TL*Era, data = cod1)
MuMIn::AICc(b,c,d,e)
#again, c and d are same. so use linear model. yes difference w TL and ERA
summary(d)
plot(d)


#################now begin to plot otolith data
#first see if difference between weights at day 2 and day 7
oto<-lm(formula = Loto_Wgt2 ~ Loto_Wgt7, data = cod1)
summary(oto)
Roto<-lm(formula = Roto_Wgt2 ~ Roto_Wgt7, data = cod1)
summary(Roto)

p <- ggplot(data = cod1, #here chose one dataframe that has all datapoints considered 
            aes(x = Roto_Wgt7,
                y = Roto_Wgt2)) +
  geom_point(size = 2.5, alpha = 0.8) + 
  theme_bw()+
  labs(y = "otolith weight 7 days", x = "otolith weight 2 days")+
  geom_smooth(method = "lm", data=cod1, aes(x=Loto_Wgt7, y = Loto_Wgt2), 
              color= "red") + 
  geom_smooth(method = "lm", data=cod1, aes(x=Roto_Wgt7,y=Roto_Wgt2), color = "grey", alpha = 0.3)
p
#very clear that no weight difference between day 2 or day 7


########otoliths
#Test if relationship between fish length and otolith length
TL_Loto <- lm(Loto_Len ~ TL, data = cod1, na.omit = T)
summary(TL_Loto)
plot(TL_Loto)
TL_Roto<-lm(formula = Roto_Len ~ TL, data = cod1)
summary(TL_Roto)

#try to filter data and then use in LL plot
codA <- filter(cod1, Era == "A", preserve = TRUE)
head(codA)
codB <- filter(cod1, Era == "B", preserve = TRUE)
head(codB)

LL <- ggplot(data = cod1,
            aes(x = TL,
                y = Loto_Len)) +
  geom_point(size = 2.5, alpha = 0.8) + 
  theme_bw()+
  labs(y = "Left otolith length", x = "Total Length (mm)") +
  geom_smooth(method = "lm", data=codB, aes(x=TL,y=Loto_Len, color = "blue"), color = "blue", alpha = 0.3, fill = "blue") +
  geom_point(data = codB, color = "blue", size = 2.5, alpha = 0.8)  +
  geom_smooth(method = "lm", data=codA, aes(x=TL,y=Loto_Len, color = "red"), color = "red", alpha = 0.3, fill = "red") +
  geom_point(data = codA, color = "red", size = 2.5, alpha = 0.8)  
LL


RL <- ggplot(data = cod1,
             aes(x = TL,
                 y = Roto_Len)) +
  geom_point(size = 2.5, alpha = 0.8) + 
  theme_bw()+
  labs(y = "Right otolith length", x = "Total Length (mm)") +
  geom_smooth(method = "lm", data=codB, aes(x=TL,y=Roto_Len, color = "blue"), color = "blue", alpha = 0.3, fill = "blue") +
  geom_point(data = codB, color = "blue", size = 2.5, alpha = 0.8)  +
  geom_smooth(method = "lm", data=codA, aes(x=TL,y=Roto_Len, color = "red"), color = "red", alpha = 0.3, fill = "red") +
  geom_point(data = codA, color = "red", size = 2.5, alpha = 0.8)  

RL



#test if otolith height (related to oto increment diameter) varies by ERA
plot_LotoHgt <- ggplot(data = cod1, aes(x = factor(Era_year), y = Loto_hgt, fill = Era)) +
  geom_boxplot()+
  theme_minimal()+
  labs(x = "Era", y = "Height Left otolith")+
  theme(legend.position ="bottom") 

plot(plot_LotoHgt)

plot_RotoHgt <- ggplot(data = cod1, aes(x = factor(Era_year), y = Roto_hgt, fill = Era)) +
  geom_boxplot()+
  theme_minimal()+
  labs(x = "Era", y = "Height Right otolith")+
  theme(legend.position ="bottom") 

plot(plot_RotoHgt)

LotoH <- lm(formula = Loto_hgt ~ TL + Era, data = cod1)
summary(LotoH)
RotoH <- lm(formula = Roto_hgt ~ TL + Era, data = cod1)
summary(RotoH)
RotoH <- lm(formula = Roto_hgt ~ Era * TL, data = cod1)
summary(RotoH)

#Test if relationship between L and R otolith weight. Ideally choose one to use
##then test if otolith weight varies by Era, as oto weight represents growth

plot_Lotowgt <- ggplot(data = cod1, aes(x = factor(Era_year), y = Loto_Wgt2, fill = Era)) +
  geom_boxplot()+
  theme_minimal()+
  labs(x = "Era", y = "Weight Left otolith")+
  theme(legend.position ="bottom") 

plot(plot_Lotowgt)

plot_Rotowgt <- ggplot(data = cod1, aes(x = factor(Era_year), y = Roto_Wgt2, fill = Era)) +
  geom_boxplot()+
  theme_minimal()+
  labs(x = "Era", y = "Weight Right otolith")+
  theme(legend.position ="bottom") 

plot(plot_Rotowgt)
#but is there a relationship between fish length and otolith weight? 

oto_eraL <- lm(Loto_Wgt2 ~ TL + Era, data = cod1, na.omit = T)
summary(oto_eraL)
oto_eraR<-lm(Roto_Wgt2 ~ TL + Era, data = cod1, na.omit = T)
summary(oto_eraR)
#otolith weights are related to fish length but not to ERA

LW <- ggplot(data = cod1,
             aes(x = TL,
                 y = Loto_Wgt2)) +
  geom_point(size = 2.5, alpha = 0.8) + 
  theme_bw()+
  labs(y = "Left otolith weight", x = "Total Length (mm)") +
  geom_smooth(method = "lm", data=codB, aes(x=TL,y=Loto_Wgt2, color = "blue"), color = "blue", alpha = 0.3, fill = "blue") +
  geom_point(data = codB, color = "blue", size = 2.5, alpha = 0.8)  +
  geom_smooth(method = "lm", data=codA, aes(x=TL,y=Loto_Wgt2, color = "red"), color = "red", alpha = 0.3, fill = "red") +
  geom_point(data = codA, color = "red", size = 2.5, alpha = 0.8)  
LW


RW <- ggplot(data = cod1,
             aes(x = TL,
                 y = Roto_Wgt2)) +
  geom_point(size = 2.5, alpha = 0.8) + 
  theme_bw()+
  labs(y = "Right otolith weight", x = "Total Length (mm)") +
  geom_smooth(method = "lm", data=codB, aes(x=TL,y=Roto_Wgt2, color = "blue"), color = "blue", alpha = 0.3, fill = "blue") +
  geom_point(data = codB, color = "blue", size = 2.5, alpha = 0.8)  +
  geom_smooth(method = "lm", data=codA, aes(x=TL,y=Roto_Wgt2, color = "red"), color = "red", alpha = 0.3, fill = "red") +
  geom_point(data = codA, color = "red", size = 2.5, alpha = 0.8)  

RW


#ggsave(filename = "C:/Users/alask/Documents/Git/KDSP_gadids/output/cod_kwet_by_age_year.png", 
#       width = 8,
#       height = 6, units = "in")
head(cod1)

codLW <- ggplot(data = cod1,
             aes(x = TL,
                 y = fresh_body_wgt_wet)) +
  geom_point(size = 2.5, alpha = 0.8) + 
  theme_bw()+
  labs(y = "Body Weight (wet)", x = "Total Length (mm)") +
  geom_smooth(method = "lm", data=codB, aes(x=TL,y=fresh_body_wgt_wet, color = "blue"), color = "blue", alpha = 0.3, fill = "blue") +
  geom_point(data = codB, color = "blue", size = 2.5, alpha = 0.8)  +
  geom_smooth(method = "lm", data=codA, aes(x=TL,y=fresh_body_wgt_wet, color = "red"), color = "red", alpha = 0.3, fill = "red") +
  geom_point(data = codA, color = "red", size = 2.5, alpha = 0.8)  

codLW
