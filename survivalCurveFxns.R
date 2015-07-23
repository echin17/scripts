library("coxphw")
library("OIsurv")

#####################################################################################################################
#use the function "createSurvivalCurves" to generate a survival curve, significance statistics,
# and schoenfeld residuals graphs (to test proportional-hazards assumption)
# createSurvivalCurves input:
# listOfNames = list of names for the columns
# fileLoc = location of file containing data of variables you are interested in exploring
#         File Format: Data from each variable must be in columns
#               C1: outcome     C2: days alive    C3-Cn: covariate data
# NOTE: continuous data is broken down into groups. You must alter the function itself if you want to use it for
#       different continuous variables
#
# Created by Elizabeth Chin. 07/23/215. Please email lizchin317ATgmailDOTcom if you have any questions
#####################################################################################################################

getSurvival <- function(dat, name, gp)
{
  outcome <- dat[,1]
  days <- dat[,2]
  var <- dat[,3]
  csurv <- Surv(days,outcome)~var
  if(is.na(gp[1])){
    gp = unique(dat[,3])
  }
  
  jpeg(paste0(name,"_KM.jpg"))
  plot(survfit(csurv), main=paste0("K-M estimate with 95% CI for ",name),
       xlab="time",ylab="survival function", xlim=c(0,12), col=1:length(gp),lwd=2)
  legend("bottomleft",legend=gp, col=1:length(gp), lty=1,lwd=2)
  dev.off()
  
  time.dep <- coxph(csurv, na.action=na.exclude)
  time.dep.zph <- cox.zph(time.dep, transform = 'km') 
  
  sink(paste0(name,"_log.txt"),append=FALSE, split=FALSE)
  print("Log-rank:")
  print(survdiff(csurv, rho=0))
  print("")
  print("Wilcoxon test:")
  print(survdiff(csurv, rho=1))
  print("")
  print("Cox Proportional Hazards Test:")
  print(coxphw(csurv, data=dat, template="PH"))
  print("")
  print("Average Regression Hazard Adjustment (Xu et al., 2000):")
  print(coxphw(csurv, data=dat, template="ARE"))
  print("")
  print("Weighted Cox Regression (Schemper et al., 2009):")
  print(coxphw(csurv, data=dat, template="AHR"))
  print("")
  print("Schoenfeld Residuals")
  print(time.dep.zph)
  sink()
  
  jpeg(paste0(name,"_SD.jpg"))
  plot(time.dep.zph[1], main=paste0(name," Schoenfeld Residuals"))
  abline(h=0, lty=3)
  dev.off()  
}

getGroup <- function(var, intervals)
{
  options(warn=-1)
  var <- as.numeric(as.character(var))
  options(warn=0)
  groups <- rep(NA,length(var))
  intervals <- as.numeric(intervals)
  numGroups <- length(intervals) - 1
  groups[var < intervals[1]] <- 1
  if(numGroups > 0)
  {
    for(i in (1:numGroups))
    {
      groups[var >= intervals[i] & var < intervals[i+1]] <- as.integer(i + 1)
    }    
  }
  groups[var >= tail(intervals, n=1)] <- numGroups + 2
  return(groups)
}

intervalToString <- function(intervals, units)
{
  numGroups <- length(intervals) + 1
  groups <- rep("",numGroups)
  groups[1] <- paste("Under", intervals[1],units, sep= " ")
  if(numGroups > 2)
  {
    for(i in (1:(numGroups-2)))
    {
      groups[i+1] <- paste0(intervals[i], "-", intervals[i+1], " ",units) 
    }    
  }
  groups[numGroups] <- paste(tail(intervals, n=1),units, "and over", sep=" ")
  return(groups)
}

createSurvivalCurves <- function(listOfNames, fileLoc)
{
  topVar <- read.csv(fileLoc, stringsAsFactors=FALSE)
  colnames(topVar) <- c("OUTCOME","NODA",listOfNames)
  for (i in (1:length(listOfNames)))
  {
    print(listOfNames[i])
    var <- topVar[,i+2]
    if(grepl("Max_Systolic_BP", listOfNames[i]))
    {
      iv <- c(120,140,160)
      gp <- intervalToString(iv,"mmHg")
      var <- getGroup(var, iv)
    } else if(grepl("Temperature", listOfNames[i])){
      iv <- c(36.5,37.2,38.9,39.4)
      gp <- intervalToString(iv, "C")
      var <- getGroup(var, iv)
    } else if(grepl("Distance", listOfNames[i])){
      iv <- c(75,150,225,300)
      gp <- intervalToString(iv, "km")
      var <- getGroup(var, iv)
    } else if(grepl("Respiratory_Rate",listOfNames[i])){
      iv <- c(12,18,25,30)
      gp <- intervalToString(iv, "bpm")
      var <- getGroup(var, iv)
    } else if(grepl("Age",listOfNames[i])){
      iv <- c(25,40,55)
      gp <- intervalToString(iv, "years")
      var <- getGroup(var, iv)
    } else if(grepl("Blood_Urea_Nitrogen", listOfNames[i])){
      iv <- c(7,20)
      gp <- intervalToString(iv, "mg/dL")
      var <- getGroup(var, iv)
    } else if(grepl("Creatinine", listOfNames[i])){
      iv <- c(0.5,1.7,7.7)
      gp <- intervalToString(iv, "mg/dL")
      var <- getGroup(var, iv)
    } else if(grepl("White_Blood_Cell", listOfNames[i])){
      iv <- c(3.5,10.5)
      gp <- intervalToString(iv, "wbc/L")
      var <- getGroup(var, iv)
    } else if(grepl("Granulocytes", listOfNames[i])){
      iv <- c(0.5,1.7,7.7)
      gp <- intervalToString(iv, "10^3/mm^3")
      var <- getGroup(var, iv)
    } else if(grepl("Glutamic_Oxaloacetic_Transaminase", listOfNames[i])){
      iv <- c(5.0,45,0,250.0)
      gp <- intervalToString(iv, "U/L")
      var <- getGroup(var, iv)
    } else {
      gp <- NA
      var <- topVar[,2+i]
    }
    name <- listOfNames[i]
    outcome <- topVar$OUTCOME
    days <- topVar$NODA
    lengths <- max(c(length(outcome), length(days), length(var)))
    length(outcome) <- lengths
    length(days) <- lengths
    length(var) <- lengths
    options(warn=-1)
    dat <- cbind(as.numeric(as.character(outcome)), as.numeric(as.character(days)), as.numeric(as.character(var)))
    options(warn=0)
    dat <- as.data.frame(na.omit(dat))
    colnames(dat) <- c("outcome","days","var")
    getSurvival(dat, name, gp)
  }
}
