# function to calculate the sub data from a for-loop 
drm.1 <- function(subdata, conc, measure, straiName,drmMethod = c(),...) {
  # step to evaluate the lowest boundary reach the limit
  maxConc <- max(as.numeric.factor(subdata$measure))
  minConc <- min(as.numeric.factor(subdata$measure))
  halfGrowth <-  mean(subdata$measure[conc == minConc])/2 # half of control Growth
  minGrowth <- min(subdata$measure[conc == maxConc])
  if (minGrowth > halfGrowth){
    return(data.frame(row.names = sprintf("@NoEC50 %s", straiName), 
               EC50 = maxConc, StdErr = NA, Lower = NA,
               Upper = NA))
  }

# step to determine if there is an hormetic effect

# step to calculate EC50 with the proper model
  
  model.1 <- drm(growthRate ~ Conc, 
                 format(straiName, trim=TRUE),
                 fct=LL.4(names=c("slope",
                                  "Lower Limit",
                                  "Upper Limit","EC50")),
                 data=subData)
# return EC50 tables
  
  
} 

for (strainSet in levels(fullPlateData$setStrain) ){
  # subset data by strain set
  subData <- subset(fullPlateData, 
                    subset = fullPlateData$setStrain == strainSet)
  maxConc <- max(as.numeric.factor(subdata$growthRate))
  minConc <- min(as.numeric.factor(subdata$growthRate))
  halfGrowth <-  mean(subdata$growthRate[conc == minConc])/2 # half of control Growth
  minGrowth <- min(subdata$growthRate[conc == maxConc])
  if (minGrowth > halfGrowth){
    noEC50Df <- data.frame(row.names = sprintf("@NoEC50 %s", straiName), 
                           EC50 = maxConc, StdErr = NA, Lower = NA,
                           Upper = NA)
    
  }
  else {
    ED_table <- ED(subData.m1,c(50), interval="delta")
    ModelSummary <- summary(subData.m1)
    # include pvalue along with EC50 calculation
    ED_table_p <- cbind(ED_table,pvalue = ModelSummary$coefficients[4,4])
    EC50Raw <- rbind(EC50Raw, ED_table_p)
    if (write2File){ # if I want to write tab delimited file
      write.table(format(ED_table_p,digits=3,trim=TRUE), eol = "\n",
                  "./data/EC50_diameter_diamR.txt", sep = "\t",
                  col.names = FALSE, quote = FALSE, row.names=TRUE)}
  }
}




testDf <- EC50SumFilterPlate[,1:4]

testDf.1 <-
  rbind(data.frame(row.names = "TestStrain", EC50 = 30, StdErr = NA,
                   Lower = NA, Upper = NA),testDf, deparse.level = 0)

EC50SumFilterPlate[1,]

dd <- 10

rbind(1:4, c = 2, "a++" = 10, dd, deparse.level = 0)

data.frame(row.names = sprintf("@NoEC50 %s", straiName), EC50 = sprintf(">%d", 30), StdErr = NA,
           Lower = NA, Upper = NA)

sprintf("@NoEC50 %s", "MIIC-14-01")
