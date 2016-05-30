# unused code

# find duplicated run strain names
dupNamePlate <- unique(
  plateEC50Simple$Strain[duplicated(plateEC50Simple$Strain)]
)
dupNameConidia <- unique(
  conidiaEC50Simple$Strain[duplicated(conidiaEC50Simple$Strain)]
)
# find unique run strain names
plateEC50SimpleUniq <- 
  plateEC50Simple[!plateEC50Simple$Strain %in% dupNamePlate,]
conidiaEC50SimpleUniq <- 
  conidiaEC50Simple[!conidiaEC50Simple$Strain %in% dupNameConidia,]

# calculate mean, stderr, Lower and Upper for duplicated strain
plateEC50SimpleDupl  <- 
  plateEC50Simple[plateEC50Simple$Strain %in% dupNamePlate,] %>% 
  group_by(Strain) %>% 
  summarise(plateECMean = mean(plateEC50), 
            StdErr = sd(plateEC50) / sqrt(length(plateEC50)),
            Lower = plateECMean - StdErr,
            Upper = plateECMean + StdErr) %>%
  rename(plateEC50 = plateECMean)

conidiaEC50SimpleDupl <- 
  conidiaEC50Simple[conidiaEC50Simple$Strain %in% dupNameConidia,] %>% 
  group_by(Strain) %>% 
  summarise(conidiaECMean = mean(conidiaEC50), 
            StdErr = sd(conidiaEC50) / sqrt(length(conidiaEC50)),
            Lower = conidiaECMean - StdErr,
            Upper = conidiaECMean + StdErr) %>%
  rename(conidiaEC50 = conidiaECMean)





###########
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


# bad fit strain with BC model
selStrain1 <- c("MIBer_A5", "STJ3a")
selStrain2 <- c("MIBer_A5")
selStrain1Data <- fullPlateData[fullPlateData$Strain %in% selStrain1,]
selStrain1Model <- drm(growthRate~Conc,
                       format(Strain,trim=TRUE),
                       fct=LL.4(names=c("slope","Lower","Upper","EC50")),
                       data=selStrain1Data)
plot(selStrain1Model,type = "all",
     broken = TRUE,
     col = TRUE,
     cex.legend = 0.9,
     xlab = expression('Log10 Transformed Fluopyram Concentration (mg/L)'),
     ylab = "Mycelial Growth Rate (cm/d)",
     cex.axis = 1,
     main = "Hormetic effect")

selStrain1.m0 <- drm(growthRate~Conc,
                     format(Strain,trim=TRUE), 
                     fct=LL.3(),
                     data=selStrain1Data)
selStrain1.m1 <- drm(growthRate~Conc,
                       format(Strain,trim=TRUE),
                       fct=LL.4(names=c("slope","Lower","Upper","EC50")),
                       data=selStrain1Data)
selStrain1.m2 <- drm(growthRate~Conc,
                     format(Strain,trim=TRUE),
                     fct=LL.5(),
                     data=selStrain1Data)
selStrain1.m3 <- drm(growthRate~Conc,
                     format(Strain,trim=TRUE), fct=CRS.4c(),
                     data=selStrain1Data)
selStrain1.m4 <- drm(growthRate~Conc,
                     format(Strain,trim=TRUE), fct=BC.4(),
                     data=selStrain1Data)
selStrain1.m5 <- drm(growthRate~Conc,
                     format(Strain,trim=TRUE), fct=W1.4(),
                     data=selStrain1Data)
selStrain1.m6 <- drm(growthRate~Conc,
                     format(Strain,trim=TRUE), fct=W2.4(),
                     data=selStrain1Data)
anova(selStrain1.m4, selStrain1.m1)
modelFit(selStrain1.m4)

plot(selStrain1.m4)
plot(selStrain1.m1,add = TRUE, col="Red")

plot(selStrain1.m3,type = "all",
     broken = TRUE,
     col = TRUE,
     cex.legend = 0.9,
     xlab = expression('Log10 Transformed Fluopyram Concentration (mg/L)'),
     ylab = "Mycelial Growth Rate (cm/d)",
     cex.axis = 1)

modelRSS <- c(m0=modelFit(selStrain1.m0)[2,2],m1=modelFit(selStrain1.m1)[2,2],
m2=modelFit(selStrain1.m2)[2,2], m3=modelFit(selStrain1.m3)[2,2],
m4=modelFit(selStrain1.m4)[2,2], m5=modelFit(selStrain1.m5)[2,2],
m6=modelFit(selStrain1.m6)[2,2])
sort(modelRSS)

ED(selStrain1.m1,50)
ED(selStrain1.m4,50)

mselect(selStrain1.mTest, list(LL.3(), LL.5(), W1.4(), W2.4(), CRS.4c(), BC.4()))

selStrain1.mTest <- drm(growthRate~Conc,format(Strain,trim=TRUE),
                     data=selStrain1Data, fct = LL.4())

selStrain2 <- c("12@KSSH_A2","12@MIVB_A5")
selStrain2Data <- fullPlateData[fullPlateData$setStrain %in% selStrain2,]
selStrain2Model <- drm(growthRate~Conc,
                       format(Strain,trim=TRUE),
                       fct=LL.4(names=c("slope","Lower","Upper","EC50")),
                       data=selStrain2Data)
modelFit(selStrain2Model)
coef(selStrain2Model) 


# hormetic effect 
subData.test <- filter(fullPlateDataExclude50, setStrain == "2@MIVB_A5")
subData.mTest <- drm(growthRate ~ Conc,
                     format(setStrain, trim=TRUE),
                     fct=LL.4(),
                     data=subData.test)
plot(subData.mTest )

plot(subData.mTest,type = "all",
     broken = TRUE,
     col = TRUE,
     cex.legend = 0.9,
     #level=c("KSSH_A2","MIVB_A5"),
     xlab = expression('Log10 Transformed Fluopyram Concentration (mg/L)'),
     ylab = "Mycelial Growth Rate (cm/d)",
     cex.axis = 1,
     main = paste(subData.mTest$data[,4][1]),
     legendPos = c(95,0.14),
     ylim = c(0,0.14))

# mselect(subData.m3, list(LL.3(),BC.4()))

subData.test <- filter(coniNoNaDataExclude50, Isolate == "13Fv165")
subData.mTest <- drm(germinate ~ Conc,
                     format(setStrain, trim=TRUE),
                     fct=LL.4(),
                     data=subData.test)
plot(subData.mTest, type = "all")

mselect(subData.mTest, list(LL.3(), LL.4(), W1.4(), W2.4(), CRS.4c(), BC.4(), BC.5()))

# check for relative growth inhibition at highest concentration

######## 2016-05-28 ###################
# revised Table S1
# mycelial assay
fullPlateData$Strain <- do.call(str_replace_all, 
                                list(fullPlateData$Strain, "_", "-"))
subPlate50Ppm <- fullPlateData[fullPlateData$Conc == 50,]

plateRelInhPerc <- subPlate50Ppm %>% group_by(Strain) %>%
  summarise(N = length(relGrowthRate),
            relInhibition = mean(relGrowthRate),
            relSE = sd(relGrowthRate) / sqrt(N))

newECTablePlate1 <- merge(EC50PlateTable, plateRelInhPerc, by="Strain")

isolateNotInFinalTable <- 
  unique(fullPlateData$Strain)[!unique(fullPlateData$Strain) 
                               %in% EC50PlateTable$Strain]

isolateNotinFinalDF <- fullPlateData[fullPlateData$Strain %in% 
                                       isolateNotInFinalTable,]
plateRelInhPerc_1 <-
  isolateNotinFinalDF[isolateNotinFinalDF$Conc == 50,] %>%
  group_by(Strain) %>%
  summarise(N = length(relGrowthRate),
            relInhibition = mean(relGrowthRate),
            relSE = sd(relGrowthRate) / sqrt(N))

# Measurment error at INMO-C4, with relative growth inhibition at 200% at 50 ppm
plateRelInhPerc_1 <- plateRelInhPerc_1[plateRelInhPerc_1$Strain != "INMO-C4",]

plateRelInhPerc_1_geo <- merge(y=plateRelInhPerc_1, x=strainInfoAll, 
                               by.y = "Strain", by.x = "isolate")

naDF <- data.frame(EC50=12*NA, StdErr=12*NA,Lower=12*NA, Upper=12*NA)
plateRelInhPerc_2 <- cbind(plateRelInhPerc_1_geo[,1:4],
                           naDF,plateRelInhPerc_1_geo[,5:7])

names(plateRelInhPerc_2) <- names(newECTablePlate1)
newECTablePlate <- rbind(newECTablePlate1, plateRelInhPerc_2)



############### Conidia Table S2 ######################

subConidia20Ppm <- coniNoNaData[coniNoNaData$Conc == 20,]

coniNoNaData$Isolate <- do.call(str_replace_all, 
                                list(coniNoNaData$Isolate, "_", "-"))

conidiaIsolateNo20ppm <- coniNoNaData[!(coniNoNaData$Isolate %in% 
                                          conidiaRelInhPerc$Isolate),]

conidiaRelInhPerc <- subConidia20Ppm %>% group_by(Isolate) %>%
  summarise(N = length(relGerminate),
            relInhibition = mean(relGerminate),
            relSE = sd(relGerminate) / sqrt(N))
# conidiaRelInhPerc[conidiaRelInhPerc$Isolate %in% "13Fv190",]

newECTableConidia1 <- merge(EC50ConidialTable, conidiaRelInhPerc, 
                            by.x="Strain", by.y="Isolate")

conidiaRelInhPerc_1 <-
  conidiaIsolateNo20ppm %>% group_by(Isolate) %>%
  summarise(N = length(relGerminate),
            relInhibition = mean(relGerminate),
            relSE = sd(relGerminate) / sqrt(N))

conidiaRelInhPerc_1_geo <- merge(y=conidiaRelInhPerc_1, x=strainInfoAll, 
                                by.x = "isolate", by.y = "Isolate")

naDFby4 <- data.frame(EC50=4*NA, StdErr=4*NA,Lower=4*NA, Upper=4*NA)
conidiaRelInhPerc_2 <- cbind(conidiaRelInhPerc_1_geo[,1:4],
                           naDFby4,conidiaRelInhPerc_1_geo[,5:7])
names(conidiaRelInhPerc_2) <- names(newECTableConidia1)
newECTableConidia <- rbind(newECTableConidia1, conidiaRelInhPerc_2)

