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
strainInfoAll <- rbind(ILStrainInfoSub, MIStrainInfoSub)

EC50SumFilterPlate$strainName <- do.call(str_replace_all, 
                                 list(EC50SumFilterPlate$Strain, "_", "-"))
EC50SumFilterConidia$strainName <- do.call(str_replace_all, 
                                 list(EC50SumFilterConidia$Strain, "_", "-"))

EC50SumPlateGeoYear <- merge(EC50SumFilterPlate, strainInfoAll, 
                             by.x = "strainName", by.y = "isolate")
EC50SumConidiaGeoYear <- merge(EC50SumFilterConidia, ILStrainInfoSub, 
                               by.x = "strainName", by.y = "isolate")


accuEC <- function(data, method = c(),...){
  

}
