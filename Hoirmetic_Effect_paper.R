# load packages
library("drc")
library("ggplot2")
library("dplyr")
library("stringr")

# generate Rdata file 
hormeticData <- read.csv("./Data/hormeticData.csv")
#unique(hormeticData$setStrain)

# function to calculate EC50
calEC50 <- function(hData, fct){
  # subset data by strain set
  subData.m1 <- drm(growthRate ~ Conc, 
                    format(setStrain, trim=TRUE),
                    fct=fct,
                    data=hData)
  ED_table <- ED(subData.m1,c(50), interval="delta", bound = FALSE)
  return(ED_table)
}

# calculate EC50 using BC.4
EC50RawHorm.BC4 <- data.frame(straiName = character(0), EC50 = numeric(0), 
                          StdErr = numeric(0), Lower = numeric(0), 
                          Upper = numeric(0))
for (strainSet in levels(factor(hormeticData$setStrain)) ){
  subData.4 <- filter(hormeticData, setStrain == strainSet)
  tableEC50 <- calEC50(subData.4, BC.4())
  EC50RawHorm.BC4 <- rbind(EC50RawHorm.BC4, tableEC50)
}
EC50RawHorm.BC4$method <- "BC4"

# calculate EC50 using LL.4
EC50RawHorm.LL4 <- data.frame(straiName = character(0), EC50 = numeric(0), 
                              StdErr = numeric(0), Lower = numeric(0), 
                              Upper = numeric(0))
for (strainSet in levels(factor(hormeticData$setStrain)) ){
  subData.4 <- filter(hormeticData, setStrain == strainSet)
  tableEC50 <- calEC50(subData.4, LL.4())
  EC50RawHorm.LL4 <- rbind(EC50RawHorm.LL4, tableEC50)
}
EC50RawHorm.LL4$method <- "LL4"

# calculate EC50 using W2.4
EC50RawHorm.W24 <- data.frame(straiName = character(0), EC50 = numeric(0), 
                              StdErr = numeric(0), Lower = numeric(0), 
                              Upper = numeric(0))
for (strainSet in levels(factor(hormeticData$setStrain)) ){
  subData.4 <- filter(hormeticData, setStrain == strainSet)
  tableEC50 <- calEC50(subData.4, W2.4())
  EC50RawHorm.W24 <- rbind(EC50RawHorm.W24, tableEC50)
}
EC50RawHorm.W24$method <- "W2.4"

# calculate EC50 using W2.4
EC50RawHorm.CRS4 <- data.frame(straiName = character(0), EC50 = numeric(0), 
                              StdErr = numeric(0), Lower = numeric(0), 
                              Upper = numeric(0))
for (strainSet in levels(factor(hormeticData$setStrain)) ){
  subData.4 <- filter(hormeticData, setStrain == strainSet)
  tableEC50 <- calEC50(subData.4, CRS.4c())
  EC50RawHorm.CRS4 <- rbind(EC50RawHorm.CRS4, tableEC50)
}
EC50RawHorm.CRS4$method <- "CRS4"

# merge all data into full table
# summarise all data together
allHormData <- rbind(EC50RawHorm.BC4, EC50RawHorm.LL4, EC50RawHorm.W24,
                     EC50RawHorm.CRS4)
allHormData <- allHormData %>% rename(EC50 = Estimate)

# check how many isolates with EC50 over 30 ppm
table(allHormData[allHormData$EC50>30,]$method)

# plot all the histograms all together
# for easy of data visulization, the EC50 values above 30 ppm
# were remove for the models BC4 and CRS4
# generate histogram for each model
ggplot(allHormData, aes(x=EC50)) + 
  geom_histogram(aes(y=(..count..)), color="dark blue", fill="white",
                 binwidth = 0.5, size=0.9) +
  # geom_vline(aes(xintercept=mean(EC50, na.rm=T)),colour="red", 
  #            linetype="dashed", size=1) +
  scale_x_continuous(breaks=seq(0,30,5), limits = c(0,30)) + 
  theme_classic() +
  labs(y="Frequency", x=expression('Fluopyram EC'[50]*' (mg/L)')) +
  facet_wrap(~method)


# make a model selection table with the statistics for each model for this model
hData.m1 <- drm(growthRate ~ Conc, 
                format(setStrain, trim=TRUE),
                fct=LL.4(),
                data=hormeticData)

hormModelSelection <- mselect(hData.m1, list(W2.4(), CRS.4c(), BC.4()))
tableModelSel <- hormModelSelection[,c(1,2,4)]




