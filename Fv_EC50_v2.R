# Copyright (c) 2015 Jie Wang
# calculation of the EC50 for fluopyram
# install.packages("drc")
library("drc")
options(digits = 4,scipen = 999)
# set working directory
setwd("/Users/wjidea/GoogleDrive/Graduate_Study/0_Research_projects/5_Fungicide_project/2014_Bayer_fluopyram/2015_02_Data_analysis_full/input_data")

# list all the files in the preset format to load
files <- dir("./", pattern = "^fungal_set[0-9]*$")
data_list <- lapply(files, read.table, header = TRUE,sep="\t")
#data_list[[6]][data_list[[6]]$Strain != "34938",]
data_list <- na.omit(data_list)
str(data_list)
# Scripts used to test which set of data is usable for the drm calculation
#   fluo_test <- data_list[[2]]
#   fluo_test$diameter1 <- with(fluo_test, 2*sqrt(Area1/3.1415926))
#   drm(diameter1 ~ Conc, 
#       format(Strain, trim=TRUE),
#       fct=LL.4(names=c("slope","Lower Limit","Upper Limit","EC50")),
#       data=fluo_test)

# A new set of data list to store the set of data that can be used to calculate in
# "drm" with the diamter 1 data
# the other list will have a convergence fail issue
data_list_diam1 <- list(data_list[[1]],data_list[[2]],data_list[[4]],data_list[[5]],
                        data_list[[8]],data_list[[9]],data_list[[13]],data_list[[14]],
                        data_list[[15]])

############# Calculate all the EC50 in for-loop ##############
# Run for loop to loop through all the data files in the list
# CAUTION! some of the relative growth area data may not work
for (fluo in data_list ){
  # Remove the space the variable names
  fluo$Strain <- gsub("\\s*$", "", fluo$Strain, perl=T)
  # back calculate the radius
  fluo$diameter2 <- with(fluo, 2*sqrt(Area2/3.1415926))
  fluo$diameter1 <- with(fluo, 2*sqrt(Area1/3.1415926))
  # calculate mycelia growth rate per day in centermeters
  fluo$diameterR <- (fluo$diameter2 - fluo$diameter1) / 14
  fluo.m1 <- drm(diameterR ~ Conc, 
                 format(Strain, trim=TRUE),
                 fct=LL.4(names=c("slope","Lower Limit","Upper Limit","EC50")),
                 data=fluo)
  ED_table <- ED(fluo.m1,c(50), interval="delta")
  write.table(format(ED_table,digits=3,trim=TRUE),"EC50_diameter_diamR.txt",
               append = TRUE,sep = "\t",col.names = FALSE, eol = "\n",
               quote = FALSE,row.names=TRUE)
}

###########################
### Figure 1 Histogram ####
###########################
library(ggplot2)
EC50_sum <- na.omit(read.table(file="/Users/wjidea/GoogleDrive/Graduate_Study/0_Research_projects/5_Fungicide_project/2014_Bayer_fluopyram/2015_02_Data_analysis_full/input_data/EC50_diameter_diamR.txt",header=TRUE,sep="\t"))
#EC50_p05 <- EC50_sum[EC50_sum$p.value< 0.05,]
#with(EC50_sum,hist(EC50,))
EC50_sum$Year <- as.factor(EC50_sum$Year)
lm_EC <- with(EC50_sum,lm(EC50~Year))
anova(lm_EC)

ggplot(EC50_sum, aes(x=Estimate)) + 
  geom_histogram(aes(y=(..count..)),binwidth=.5, 
                 colour="dark green", fill="white",size=0.9) +
  geom_vline(aes(xintercept=mean(Estimate, na.rm=T)),colour="red", 
             linetype="dashed", size=1) + theme_grey() + geom_rug() + 
  scale_x_continuous(breaks=0:10) + 
  labs(y="Frequency", x=expression('Fluopyram EC'[50]*' (mg/L)')) +
  theme(axis.text = element_text(size=16), axis.title = element_text(size=22),
        axis.text.x=element_text(colour="black"), 
        axis.text.y=element_text(colour="black")) 

#####################################################
#### Figure 2 Growth Curve fit and special cases ####
#####################################################

par(mfrow=c(1,2))
# Specific examples to show good and bad fitness of the model
# possiblely due to the inconsistent growth of the colony
fluo1 <- data_list[[3]]
fluo1$Strain <- gsub("\\s*$", "", fluo1$Strain, perl=T)
fluo1$diameter2 <- with(fluo1, 2*sqrt(Area2/3.1415926))
fluo1$diameter1 <- with(fluo1, 2*sqrt(Area1/3.1415926))
fluo1$diameterR <- (fluo1$diameter2 - fluo1$diameter1)/14
drc_func1 <- drm(diameterR~Conc,
                 format(Strain,trim=TRUE),
                 fct=LL.4(names=c("slope","Lower Limit","Upper Limit","EC50")),
                 data=fluo1)

plot(drc_func1,type = "all",
     broken = TRUE,
     col = TRUE,
     # Spaces in the variable names are essential to id those variables
     level = c("MIBer_A5","STJ_3a  "),  
     cex.legend = 0.9,
     xlab = expression('Log10 Transformed Fluopyram EC'[50]*' (mg/L)'),
     ylab = "Mycelial Growth Rate (cm/d)",
     cex.axis = 1,
     main = "Hormetic effect",
     legendPos = c(95,0.14),
     ylim = c(0,0.14))
# MIBer_A5 annotation
abline(v=3.789, col="black", lty=1)
text(x = 6.7, y = .12,labels = expression('EC'[50]*' = 3.682'),cex=.8)
# STJ_3a annotation
abline(v=3.539,col="red",lty=2)
text(x = 6.7, y = .08,labels = expression('EC'[50]*' = 3.539'),cex=.8,col="red")

# Normal examples to show good fitness of the LL.4 model
# possiblely due to the inconsistent growth of the colony
fluo2 <- data_list[[4]]
fluo2$Strain <- gsub("\\s*$","", fluo2$Strain, perl=T)
fluo2$diameter2 <- with(fluo2,2*sqrt(Area2/3.1415926))
fluo2$diameter1 <- with(fluo2,2*sqrt(Area1/3.1415926))
fluo2$diameterR <- (fluo2$diameter2 - fluo2$diameter1) / 14
drc_func2 <- drm(diameterR~Conc,format(Strain,trim=TRUE),data=fluo2,
                fct=LL.4(names=c("slope","Lower Limit","Upper Limit","EC50")))
ED(drc_func2,c(50),interval="delta")
plot(drc_func2,type="all",
     broken=TRUE,
     col=TRUE,
     level=c("KSSH_A2 ","MIVB_A5 "),
     cex.legend =0.9, 
     xlab=expression('Log10 Transformed Fluopyram EC'[50]*' (mg/L)'),
     ylab="Mycelial Growth Rate (cm/d)",
     main="Log-logistic Model",
     legendPos=c(95,0.14),
     ylim=c(0,0.14))
abline(v=2.172,col="black",lty=1)
text(x = 4, y = .12,labels = expression('EC'[50]*' = 2.172'),cex=.8)
abline(v=3.793,col="red",lty=2)
text(x = 7, y = .08,labels = expression('EC'[50]*' = 3.793'),cex=.8,col="red")

###################################
### Figure 3 Methods Validation ###
###################################
# Read in EC50 summary data to 
diam1 <- read.table("/Users/wjidea/Google\ Drive/Graduate_Study/Research_projects/5_Fungicide_project/2014_Bayer_fluopyram/2015_02_Data_analysis_full/input_data/EC50_diameter_diam1.txt",sep="\t")
diam1$V1 <- gsub("\\s*:50$","", diam1$V1, perl=T)
diam1$V6 <- "3_days"

diam2 <- read.table("/Users/wjidea/Google\ Drive/Graduate_Study/Research_projects/5_Fungicide_project/2014_Bayer_fluopyram/2015_02_Data_analysis_full/input_data/EC50_diameter_diam2.txt",sep="\t")
diam2$V1 <- gsub("\\s*:50$","", diam2$V1, perl=T)
diam2$V6 <- "10_days"

diamR <- read.table("/Users/wjidea/Google\ Drive/Graduate_Study/Research_projects/5_Fungicide_project/2014_Bayer_fluopyram/2015_02_Data_analysis_full/input_data/EC50_diameter_diamR.txt",sep="\t")
diamR$V1 <- gsub("\\s*:50$","", diamR$V1, perl=T)
diamR$V6 <- "Relative"

diam <- rbind(diam1,diam2,diamR)
colnames(diam) <- c("Strain", "EC50", "StdErr", "Lower", "Upper", "Methods")
dim(diam)
library(beeswarm)
# change the order of the methods by making the 3 days measurement in the beiginning
diam$Methods <- factor(as.factor(diam$Methods),
                       levels=levels(as.factor(diam$Methods))[c(2,1,3)] )

boxplot(EC50 ~ Methods, data = diam,outline=FALSE, 
        names = c("3 Days", "10 Days", "Relative"),
        ylab = "EC50 (mg/L)",
        xlab = "Measurement methods",
        cex.axis = 1)
beeswarm(EC50 ~ Methods, data = diam, method="hex",
         col = 2:4, pch = 16, add = TRUE)
diam[diam$Methods != "3_days",]
summary(aov(EC50 ~ Methods, data= diam[diam$Methods != "Relative",]))


###########################
### Growth Rate Figure ####
###########################
# extract growth rate from control groups 
growth_rate <- data.frame(Strain=character(0),Conc=numeric(0),Rep=integer(0),
                          Area1=numeric(0), Area2=numeric(0), RelArea=numeric(0),
                          Set=character(0), diameter2=numeric(0), 
                          diameter1=numeric(0), diameterR=numeric(0))

for (fluo in data_list){
  # Remove the space the variable names
  fluo$Strain <- gsub("\\s*$","", fluo$Strain, perl=T)
  # back calculate the radius
  fluo$diameter2 <- with(fluo,2*sqrt(Area2/3.1415926))
  fluo$diameter1 <- with(fluo,2*sqrt(Area1/3.1415926))
  # calculate mycelia growth rate per day in centermeters
  fluo$diameterR <- (fluo$diameter2 - fluo$diameter1)/14
  # select the control groups
  growth_rate <- rbind(growth_rate,fluo[fluo$Conc == 0,])
}

growth_rate[order(growth_rate$diameterR, decreasing=TRUE)]
hist(growth_rate$diameterR, 
     breaks = 15, 
     xlim=c(0,0.2), 
     main="Histogram of mycelial growth rate",
     xlab="Mycelial Growth Rate (cm/d)")
rug(growth_rate$diameterR)
summary(growth_rate$diameterR)

#########################
##### J's Sand Box ######
#########################

# Codes to calculate EC50 for nonFv species EC50
nonfv <- read.table("nonfv_isolates.txt",header=TRUE, sep = "\t")
nonfv$Strain <- gsub("\\s*$", "", nonfv$Strain, perl=T)
# back calculate the radius
nonfv$diameter2 <- with(nonfv, 2*sqrt(Area2/3.1415926))
nonfv$diameter1 <- with(nonfv, 2*sqrt(Area1/3.1415926))
nonfv$diameterR <- (nonfv$diameter2 - nonfv$diameter1) / 14
fluo.nonfv <- drm(diameter2 ~ Conc, 
               format(Strain, trim=TRUE),
               fct=LL.4(names=c("slope","Lower Limit","Upper Limit","EC50")),
               data=nonfv[nonfv$Strain!=22387,])
plot(fluo.nonfv,cex.legend = 0.7)
ED(fluo.nonfv,c(50),interval="delta")
summary(fluo.nonfv)
nonfv_mod <- nonfv[nonfv$Strain!=22387, ]
nonfv_mod <- nonfv_mod[nonfv_mod$Strain != 34938,]
# nonfv_mod <- nonfv_mod[nonfv_mod$Strain != 54362,]
# nonfv_mod <- nonfv_mod[nonfv_mod$Strain != 22395,]

fluo.nonfv_mod <- drm(diameterR ~ Conc, 
                  format(Strain, trim=TRUE),
                  fct=LL.4(names=c("slope","Lower Limit","Upper Limit","EC50")),
                  data=nonfv_mod[nonfv_mod$Strain != 34938,])
plot(fluo.nonfv_mod,cex.legend = 0.7,col=TRUE,level = c(22276,31096,22395))

# examples for hormetic response to fungicide treatment at low concentrations
hormetic <- drm(diameterR~Conc,format(Strain,trim=TRUE),
                fct=CRS.4c(),data=fluo1)
plot(hormetic,type="all",broken=TRUE,col=TRUE,level=c("MIBer_A5","STJ_3a  "),
     cex.legend = .8,xlab="Log10-Transformed Concentrations", 
     ylab="Mycelial Growth Rate (cm/d)",cex.axis=1)
ED(hormetic,c(50),interval="delta")

fluo2 <- data_list[[4]]
fluo2$Strain <- gsub("\\s*$","", fluo2$Strain, perl=T)
fluo2$diameter2 <- with(fluo2,2*sqrt(Area2/3.1415926))
fluo2$diameter1 <- with(fluo2,2*sqrt(Area1/3.1415926))
fluo2$diameterR <- (fluo2$diameter2 - fluo2$diameter1) / 14
drc_func2 <- drm(diameterR~Conc,format(Strain,trim=TRUE),data=fluo2,
                 fct=LL.4(names=c("slope","Lower Limit","Upper Limit","EC50")))
ED(drc_func2,c(50),interval="delta")
plot(drc_func2,type="all",
     broken=TRUE,
     col=TRUE,
     level=c("KSSH_A2 ","MIVB_A5 "),
     cex.legend =0.9, 
     xlab=expression('Log10 Transformed Fluopyram EC'[50]*' (mg/L)'),
     ylab="Mycelial Growth Rate (cm/d)",
     main="Log-logistic Model",
     legendPos=c(95,0.14),
     ylim=c(0,0.14))
abline(v=2.172,col="black",lty=1)
text(x = 4, y = .12,labels = expression('EC'[50]*' = 2.172'),cex=.8)
abline(v=3.793,col="red",lty=2)
text(x = 7, y = .08,labels = expression('EC'[50]*' = 3.793'),cex=.8,col="red")

# test code for package broom
library(broom)
lmfit <- lm(EC50 ~ Methods, data= diam[diam$Methods != "Relative",])
summary(lmfit)
tidy(lmfit)
head(augment(lmfit))
glance(lmfit)
td <- tidy(lmfit, conf.int = TRUE)
ggplot(td, aes(estimate, term, color = term)) +
  geom_point() +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high)) +
  geom_vline()


# rjags
set.seed(1337)
y <- rnorm(n = 20, mean = 10, sd = 5)
mean(y)

library("rjags")
install.packages("rjags")
# The model specification
model_string <- "model{
for(i in 1:length(y)) {
y[i] ~ dnorm(mu, tau)
}
mu ~ dnorm(0, 0.0001)
sigma ~ dlnorm(0, 0.0625)
tau <- 1 / pow(sigma, 2)
}"

# Running the model
model <- jags.model(textConnection(model_string), data = list(y = y), n.chains = 3, n.adapt= 10000)
update(model, 10000); # Burnin for 10000 samples
mcmc_samples <- coda.samples(model, variable.names=c("mu", "sigma"), n.iter=20000)
plot(mcmc_samples)



# BEST model

library(BEST)
y1 <- c(5.77, 5.33, 4.59, 4.33, 3.66, 4.48)
y2 <- c(3.88, 3.55, 3.29, 2.59, 2.33, 3.59)
priors <- list(muM = 6, muSD = 2)
BESTout <- BESTmcmc(y1, y2, priors=priors, parallel=FALSE)
plot(BESTout)
plot(BESTout, compVal=1, ROPE=c(-0.1,0.1))
plot(BESTout, which="sd")
