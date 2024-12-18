# This is the plotting & analysis script for:
#age-1 Pacific cod comparing diet between Eras
#all cook bay and Era A is 2013; Era B is 2023
#### NMS Ordinations by IRI of prey item in stomach
##  GC, weight percent, values were then 4th root
###one empty stomach from era A removed. currently this is NOT each stomach adds to 100%.

#### Load Libraries ####
library('vegan')
library("MASS")
library("labdsv")
library("indicspecies")
library("tidyverse")
library("ggrepel")
library("goeveg")
library("corrplot")
library("rcartocolor")
library("dplyr")
library("tidyr")
library("ggeasy")




# This is the data import script
# Read in the data and rename any columns that need renaming
prey <- read_csv("C:/Users/alask/Documents/Git/seine-data/data/cod_age1_diet.csv")
head(prey)
str(prey)
tail(prey)

##can see that last row are all NA. Must remove these last 3 rows
prey<-filter(prey, !is.na(TL))
tail(prey)
#looks good. all data rows and col present

### Create Environmental and Species Matrices ####
prey_1 <- prey[,c(3:14)]
head(prey_1)

preyEnvData <- prey[,c(1:2, 15:19)]
head(preyEnvData)

### Run a Stress Plot using the package "goeveg" where k = # dimensions
dimcheckMDS(prey_1,distance = "bray",k = 6,trymax = 300,autotransform = FALSE)
#converged quickly. good solution is 2 stress dimensions = 0.16, 3D = 0.089
#output has graph with line at 0.20 being good stress, so use 2 D as 0.16 < 0.2


#### Run the NMS Ordination dimensions k = # dimensions####
prey_wgtMDS<-metaMDS(prey_1, distance="bray", k=2, trymax=300, autotransform=FALSE)
## Set up your two NMS Principal Axes 
NMDS1 <- prey_wgtMDS$points[,1]
NMDS2 <- prey_wgtMDS$points[,2]


#Add them to the ordination data plot
codprey<-cbind(prey, NMDS1, NMDS2)
head(codprey)
str(codprey)
codprey <- codprey %>%
  mutate(ERA = as.factor(Era))

##save file w NMDS 3 dimensions. 
write.csv(codprey, file = "output/nmds_cod1_diet.csv")

#### Plot the NMS in Base R ###
plot(prey_wgtMDS, type = 't', display = c('species'))

### Plot the NMS in ggplot ##
head(preyEnvData)

#isolate factors
preyEnvDataFact <- preyEnvData[,c(2,4,5)]
head(preyEnvDataFact)
en <- envfit(prey_wgtMDS, preyEnvDataFact, permutations = 999, na.rm = T)
head(en)

#Only continuous environmental variables
head(preyEnvData)
preyEnvDataCont <- preyEnvData[,c(6,7)]
head(preyEnvDataCont)
#names(preyEnvDataCont)[1]<-paste("Day of Year")
#names(preyEnvDataCont)[2]<-paste("Temp")
#names(preyEnvDataCont)[3]<-paste("Salinity")
#names(preyEnvDataCont)[4]<-paste("TL (mm)")
#head(preyEnvDataCont)

en1 <- envfit(prey_wgtMDS, preyEnvDataCont, permutations = 999, na.rm = T)
head(en1)
en1 <- as.data.frame(en1)

#### Get the vectors the correct length
en_coord_cont = as.data.frame(scores(en1, "vectors")) * ordiArrowMul(en1)


#### Plot the NMS with environmental variable overlay in Base R ###

plot(prey_wgtMDS, type = 't', display = c('species'))
plot(en1)

#### Plot the NMS with environmental variable overlay in ggplot #####

safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                             "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")

Cod_Env <- ggplot(data=codprey, aes(NMDS1, NMDS2))+
  geom_point(data=codprey, aes(NMDS1, NMDS2, color=ERA), show.legend=F, position=position_jitter(.1))+
  stat_ellipse(aes(fill=ERA, color = ERA), alpha=.25,type='t',size =1, geom="polygon")+ 
  theme_classic()+
  scale_linetype_manual(values = "solid") +
  scale_fill_manual(values=c("#332288", "#888888", "#CC6677")) +
  scale_color_manual(values=c("#332288", "#888888", "#661100", "#000000")) +
  #geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), data = en1, size = 1, alpha = 0.5, colour = "black") +
  #geom_text(data = en1, aes(x=NMDS1, y = NMDS2), colour = "black", fontface = "bold", label = row.names(en_coord_cont), position = position_jitter(.01))+
  theme(axis.text=element_text(size=15), axis.title=element_text(size=14,face="bold")) 

plot(Cod_Env)

### Set up Prey joint-plot to overlay prey names on Season elipses
en_prey <- envfit(prey_wgtMDS, prey_1, permutations = 999, na.rm = T)
head(prey_1)
head(en_prey)

# Make a subset of Prey Variables with R^2 > 0.1 that are important to final NMDS by season
#select the species that have R^2 with significance
en_prey1 <- envfit(NMDS1, prey_1, permutations = 999, na.rm = T)
head(prey_1)
head(en_prey1)

en_prey2 <- envfit(NMDS2, prey_1, permutations = 999, na.rm = T)
head(en_prey2)

#first make the plot of env var and species
head(prey_1)
preywgt_cut <- prey_1[,c(1, 2, 3, 6, 7, 8, 10, 12)]
head(preywgt_cut)
#preywgt_cut <- rename(preywgt_cut,"Poly" = "poly")
#preywgt_cut <- rename(preywgt_cut,"Gamm" = "Gammarid")
#preywgt_cut <- rename(preywgt_cut,"Caprellid" = "Caprellidae4")
#preywgt_cut <- rename(preywgt_cut,"Shrimp" = "shrimp")

head(preywgt_cut)

en_prey_cut <- envfit(prey_wgtMDS, preywgt_cut, permutations = 999, na.rm = T)
head(en_prey_cut)
en_prey_cut <- as.data.frame(en_prey_cut)
prey_coord = as.data.frame(scores(en_prey, "vectors")) * ordiArrowMul(en_prey)
#stuck because above line will not work

#try below
prey_coord = as.data.frame(scores(en_prey_cut)$vectors * ordiArrowMul(en_prey_cut))
#not getting above line to work and need a dataframe
#going to try fortify method for envfit object

fortify(en_prey_cut)
prey_coord = as.data.frame(scores(en_prey, "vectors")) * ordiArrowMul(en_prey)
#above wont' work either


plot(prey_wgtMDS, type = 't', display = c('species'))
plot(en_prey)
#those plots overlay. not sure why below will not overlay 

#### Plot the NMS in ggplot with prey vector overlay for Month ####
Wgt_prey <- ggplot(data=codprey, aes(NMDS1, NMDS2))+
  geom_point(data=codprey, aes(NMDS1, NMDS2, color=ERA), show.legend=F, position=position_jitter(.1))+
  stat_ellipse(aes(fill=ERA, color = ERA), alpha=.2,type='t',size =1, geom="polygon")+ 
  theme_classic()+
  scale_linetype_manual(values = "solid") +
  #scale_fill_manual(values=c("#332288", "#888888", "#CC6677")) +
  #scale_color_manual(values=c("#332288", "#888888", "#661100", "#000000")) + 
  #geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), data = en_coord_cont, size = 1, alpha = 0.5, colour = "black", arrow = arrow()) +
  #geom_text(data = en_prey, aes(x=NMDS1, y = NMDS2), colour = "black", fontface = "bold", label = row.names(en_prey), position=position_jitter(0.25))+
  theme(axis.text=element_text(size=15), axis.title=element_text(size=14,face="bold")) 

plot(Wgt_prey)

plot(Wgt_prey)
plot(en_prey, add = TRUE)



##using: jkzorz.github/io/2020/04/04/NMDS-extras.html
en_coord_cont = as.data.frame(scores(en1, "vectors")) * ordiArrowMul(en1)
en_coord_cont = as.data.frame(scores(en1, "factors")) * ordiArrowMul(en1)

en = envfit(prey_wgtMDS, preyEnvDataFact, permutations = 999, na.rm = T)
head(en)
plot(en_prey_cut) #this is what I want.

gg = ggplot(data = codprey, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(data = codprey, aes(colour = ERA), size = 3, alpha = 0.5) + 
  stat_ellipse(aes(fill=ERA, color = ERA), alpha=.2,type='t',size =1, geom="polygon")+ 
  scale_colour_manual(values = c("orange", "steelblue")) + 
  theme(axis.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "grey30"), 
        axis.ticks = element_blank(), axis.text = element_blank(), legend.key = element_blank(), 
        legend.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        legend.text = element_text(size = 9, colour = "grey30")) +
  labs(colour = "ERA")

gg


plot(prey_wgtMDS, type = 't', display = c('species')) 
par(new=TRUE)
plot(en_prey) 
par(new=TRUE)
plot(gg, add = TRUE)
plot(en_prey, add = TRUE)



plot(gg) + print(en_prey_cut, add=TRUE) 

#this is what I want (above) but it doesn't work. Maybe print one on top?
#en_prey_cut is not a dataframe and cannot be coerced into one



