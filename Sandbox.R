






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


