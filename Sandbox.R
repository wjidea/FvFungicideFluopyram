as.numeric(levels(coniDataFormat$Conc))[coniDataFormat$Conc]

TestData <- subset(fullPlateData, subset = fullPlateData$setStrain == "9@VB-1")

# back calculate the radius
TestData$diameter2 <- with(TestData, 2*sqrt(Area2/3.1415926))
TestData$diameter1 <- with(TestData, 2*sqrt(Area1/3.1415926))
# calculate mycelia growth rate per day in centermeters
TestData$diameterR <- (TestData$diameter2 - TestData$diameter1) / 14
TestData.m1 <- drm(diameterR ~ Conc, 
                  format(setStrain, trim=TRUE),
                  fct=LL.4(names=c("slope","Lower Limit","Upper Limit","EC50")),
                  data=TestData)
ED_table <- ED(TestData.m1,c(50), interval="delta")
testSummary <- summary(TestData.m1)
ED_table_p <- cbind(ED_table,pvalue = testSummary$coefficients[4,4])


# import demographic data for MI and IL experiment
setwd(dataPath)
MIStrainInfo <- read.csv("StrainInfoMI.csv", header = TRUE )
MIStrainInfoSub <- subset(MIStrainInfo, 
                          select = c("isolate","state","county","year"))
ILStrainInfo <- read.csv("StrainInfoIL.csv", header = TRUE )
ILStrainInfoSub <- subset(ILStrainInfo,
                          select = c("IsolateID", "State", "County", "year"))
colnames(ILStrainInfoSub) <- c("isolate","state","county","year")
strainInfoAll <- rbind(ILStrainInfoSub[1:53,], MIStrainInfoSub)

EC50SumFilterPlate$strainName <- do.call(str_replace_all, 
                                 list(EC50SumFilterPlate$Strain, "_", "-"))
EC50SumFilterConidia$strainName <- do.call(str_replace_all, 
                                 list(EC50SumFilterConidia$Strain, "_", "-"))

EC50SumPlateGeoYear <- merge(EC50SumFilterPlate, strainInfoAll, 
                             by.x = "strainName", by.y = "isolate")
EC50SumConidiaGeoYear <- merge(EC50SumFilterConidia, ILStrainInfoSub, 
                               by.x = "strainName", by.y = "isolate")
dim(EC50SumPlateGeoYear)
dim(EC50SumFilterPlate)

EC50SumPlateGeoYear[!(EC50SumPlateGeoYear$strainName %in% EC50SumFilterPlate$strainName),]
EC50SumFilterPlate[!(EC50SumFilterPlate$strainName %in% EC50SumPlateGeoYear$strainName),]

EC50SumConidiaGeoYear[!(EC50SumConidiaGeoYear$strainName %in% EC50SumFilterConidia$strainName),]
EC50SumFilterConidia[!(EC50SumFilterConidia$strainName %in% EC50SumConidiaGeoYear$strainName),]


# Function to calculate accurate EC50 values
accuEC <- function(data, method = c(),...){
  

}



# correct strain name
EC50SumConidiaGeoYear$strainName[EC50SumConidiaGeoYear$strainName == "Mont-1"] <- "Mont1"
# merge two EC50 df
plateEC50Simple <- select(EC50SumPlateGeoYear, strainName, EC50, StdErr, 
                          Lower, Upper) %>% rename(plateEC50 = EC50)
conidiaEC50Simple <- select(EC50SumConidiaGeoYear, strainName, EC50, StdErr, 
                          Lower, Upper) %>% rename(conidiaEC50 = EC50)
dupNamePlate <- unique(plateEC50Simple$strainName[duplicated(plateEC50Simple$strainName)])
dupNameConidia <- unique(conidiaEC50Simple$strainName[duplicated(conidiaEC50Simple$strainName)])

plateEC50SimpleUniq <- plateEC50Simple[!plateEC50Simple$strainName %in% dupNamePlate,]
conidiaEC50SimpleUniq <- conidiaEC50Simple[!conidiaEC50Simple$strainName %in% dupNameConidia,]

plateEC50SimpleDupl  <- 
  plateEC50Simple[plateEC50Simple$strainName %in% dupNamePlate,] %>% 
  group_by(strainName) %>% 
  summarise(plateECMean = mean(plateEC50), 
            StdErr = sd(plateEC50) / sqrt(length(plateEC50)),
            Lower = plateECMean - StdErr,
            Upper = plateECMean + StdErr) %>%
  rename(plateEC50 = plateECMean)
  
conidiaEC50SimpleDupl <- 
  conidiaEC50Simple[conidiaEC50Simple$strainName %in% dupNameConidia,] %>% 
  group_by(strainName) %>% 
  summarise(conidiaECMean = mean(conidiaEC50), 
            StdErr = sd(conidiaEC50) / sqrt(length(conidiaEC50)),
            Lower = conidiaECMean - StdErr,
            Upper = conidiaECMean + StdErr) %>%
  rename(conidiaEC50 = conidiaECMean)

plateEC50 <- rbind(plateEC50SimpleUniq, plateEC50SimpleDupl) %>%
  rename(pEC50 = plateEC50, pSE = StdErr, pL = Lower, pU = Upper)
conidiaEC50 <- rbind(conidiaEC50SimpleUniq, conidiaEC50SimpleDupl) %>%
  rename(cEC50 = conidiaEC50, cSE = StdErr, cL = Lower, cU = Upper)
rownames(plateEC50) <- 1:length(rownames(plateEC50))
rownames(conidiaEC50) <- 1:length(rownames(conidiaEC50))

compStrainsNames <- merge(plateEC50, conidiaEC50, by = "strainName")$strainName
compPdata <- plateEC50[plateEC50$strainName %in% compStrainsNames, ]
compPdata$assay <- "plate"
colnames(compPdata) <- c("strainName", "EC50", "SE", "Lower", "Upper", "assay")
compCdata <- conidiaEC50[conidiaEC50$strainName %in% compStrainsNames,]
compCdata$assay <- "conidia"
colnames(compCdata) <- c("strainName", "EC50", "SE", "Lower", "Upper", "assay")
compData <- rbind(compPdata, compCdata)

dodge <- position_dodge(width=0.9)
ggplot(compData, aes(x = strainName, y = EC50, fill=assay)) +
  geom_bar(aes(fill = assay), stat="identity", 
           position = position_dodge(width = 0.9)) +
  geom_errorbar(aes(ymax = EC50 + SE, ymin=EC50 - SE), 
                position = position_dodge(width = 0.9), width=0.25) +
  scale_fill_grey(start = 0.4, end = 0.7) +
  theme_classic() + labs(x = "Strain names", y = "EC50 estimates") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))





























