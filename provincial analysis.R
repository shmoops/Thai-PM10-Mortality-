###############################################################################
# Code originally for the analysis in:
#
#   "Reducing and meta-analyzing estimates from distributed lag non-linear models"
#   Gasparrini and Armstrong 
#   BMC Medical Research Methodology - 2013
#   http://www.ag-myresearch.com/2013_gasparrini_bmcmrm.html
#
# *modified to suite the purposes of our research (particulate matter)
###############################################################################

###############################################################################
# CREATE 2 OBJECTS:
#
# 1) A VECTOR WITH NAMES OF REGIONS OF THAILAND
#
# 2) A LIST WITH THE DATA FOR EACH REGION, INCLUDING:
#   - DATE, YEAR, MONTH, DAY, TIME, DAY OF THE YEAR, DAY OF THE WEEK
#   - REGION NUMBERS AND NAMES
#   - DAILY MEAN OF PM10, TEMPERATURE AND RELATIVE HUMIDITY
#   - MORTALITY (ALL-CAUSE, CARDIOVASCULAR, AND RESPIRATORY)
#
###############################################################################

# LOAD PACKAGES (ASSUMED ALREADY INSTALLED)
library(dlnm) ; library(mvmeta) ; library(splines)

# LOAD THE DATASET
dim(mvmeta_MONSTER_north)
head(mvmeta_MONSTER_north)

# REGIONS
regions <- as.character(unique(mvmeta_MONSTER_north$regnames))

# CREATE A LIST WITH THE REGIONAL SERIES
data <- lapply(regions,function(x) mvmeta_MONSTER_north[mvmeta_MONSTER_north$regnames==x,])
names(data) <- regions
m <- length(regions)

# PM10 RANGES
ranges <- t(sapply(data, function(x) range(x$PM10, na.rm=TRUE)))


####################################################################
# FUNCTION TO COMPUTE THE Q-AIC IN QUASI-POISSON MODELS
fqaic <- function(model) {
  loglik <- sum(dpois(model$y,model$fitted.values,log=TRUE))
  phi <- summary(model)$dispersion
  qaic <- -2*loglik + 2*summary(model)$df[3]*phi
  return(qaic)
}
####################################################################

# Provincial Statistics:
 
for(i in seq(data)) {
  
  # LOAD
sub <- data[[i]]


print(paste("Region: ", sub$regnames[1]))

print(paste("Start Date: ", min(sub$date)))
  
print(paste("NAD: ", mean(sub$nad, na.rm = TRUE), "sd: ", sd(sub$nad, na.rm = TRUE)))

print(paste("CD: ", mean(sub$cd, na.rm = TRUE), "sd: ", sd(sub$cd, na.rm = TRUE)))

print(paste("RD: ", mean(sub$rd, na.rm = TRUE), "sd: ", sd(sub$rd, na.rm = TRUE)))

print(paste("PM10: ", mean(sub$PM10, na.rm = TRUE), "sd: ", sd(sub$PM10, na.rm = TRUE)))

print(paste("Temp: ", mean(sub$temp, na.rm = TRUE), "sd: ", sd(sub$temp, na.rm = TRUE)))

print(paste("Rh: ", mean(sub$rh, na.rm = TRUE), "sd: ", sd(sub$rh, na.rm = TRUE)))

cat("\n")

}

#Average [PM10] across all regions:
region_averages= vector(mode = "numeric", length= length(data))

for(i in seq(data)) {
  
  # LOAD
  sub <- data[[i]]
  region_averages[i]= mean(sub$PM10, na.rm = TRUE) 
  
}

mean(region_averages)


####################################################################
# For loop to produce the Figure 1 for every region
####################################################################

# ARGUMENTS AND LISTS FOR CROSS-BASIS DEFINITION
# bounty <- colMeans(ranges)
bound <- c(10, 160)
varknots <- equalknots(bound, fun = "bs", degree = 2,df=6)
argvar <- list(fun="bs",degree=2,knots=varknots,bound=bound)

lagknots <- logknots(10,df=5, intercept = TRUE)
arglag <- list(fun="ns", knots=lagknots)

lag <- c(0,10)
lag1 <- c(0,1)
lag2 <- c(0,2)
lag3 <- c(0,3) 

# GENERATING PROVINCE SPECIFIC SURFACE AND PEDICTOR SPECIFIC SLICES

for(i in seq(data)) {
  
  # LOAD
  sub <- data[[i]]

  
  # BASIS FOR PM10:
  # - QUADRATIC SPLINE FOR PREDICTOR, WITH SPECIFIC KNOT SELECTION
  # - NATURAL CUBIC SPLINE FOR LAG, WITH DF AT EQUALLY-SPACED LOG-VALUES
  # WARNING FOR PREDICTION BEYOND BOUNDARIES SUPPRESSED
  suppressWarnings(
  cb <- crossbasis(sub$PM10,lag=10,argvar=argvar,arglag=arglag)
)

   # RUN THE MODEL
  # THE BELOW IF STATEMENT IS FOR THE ANALYSIS OF NORTHERN THAILAND BECAUSE CERTAIN PROVINCES IN THE NORTH
  # DO NOT CONTAIN DATA FOR CERTAIN CO-POLLUTANTS OR YEARS
  if(length(data)== 8){
  
  if(i == 2|| i==3){

    model <- glm(nad ~ cb + ns(rh, 4) + ns(temp, 4)+ ns(O3, 4) + dow + ns(time,df=10*8), family=quasipoisson(),sub)
  
    }else if (i == 5|| i==6){
    
      model <- glm(nad ~ cb + ns(NOX,4) + ns(rh, 4) + ns(temp, 4)+ ns(O3, 4) +ns(SO2,4)+ dow + ns(time,df=10*7), family=quasipoisson(),sub)

     }else if(i == 7|| i==8){
 
        model <- glm(nad ~ cb + ns(NOX,4) + ns(rh, 4) + ns(temp, 4)+ ns(O3, 4) +ns(SO2,4)+ dow + ns(time,df=10*6), family=quasipoisson(),sub)
    
         } else{
       
       model <- glm(nad ~ cb + ns(NOX,4) + ns(rh, 4) + ns(temp, 4)+ ns(O3, 4) +ns(SO2,4)+ dow + ns(time,df=10*10), family=quasipoisson(),sub)
       
     }
  
  }else{
    
model <- glm(nad ~ cb + ns(NOX,4) + ns(rh, 4) + ns(temp, 4)+ ns(O3, 4) +ns(SO2,4)+ dow + ns(time,df=10*10), family=quasipoisson(),sub)
    
  }
  
  # PREDICTION USING:
  #   crosspred FOR BI-DIMENSIONAL RELATIONSHIP
  #   crossreduce FOR UNI-DIMENSIONAL SUMMARIES
  # (NB: CENTERING AT SPECIFIC PM10 VALUE)
  # (NB: ALL THE ESTIMATES ARE ALREADY REPORTED BY crosspred ALONE)
  
  cp <- crosspred(cb,model,from=bound[1],to=bound[2],by=2, bylag=1,cen=10)
  
  crall <- crossreduce(cb,model, type = "overall", lag= lag1, from=bound[1],to=bound[2],  by=0.2,cen=10)

  crlag <- crossreduce(cb,model,type="lag",value=1, from=bound[1],to=bound[2],
                       bylag=0.2,cen=10)
 
  crvar <- crossreduce(cb,model,type="var",value=100, from=bound[1],to=bound[2],
                       bylag=0.2,cen=10)
  
  # PLOT
  par(mar=c(1.5,1,0,0)+0.1,cex.axis=0.9,cex.lab=1)
  layout(matrix(rep(1:4,each=2),2,4,byrow=TRUE))
  
  # 3D PLOT WITH DIFFERENT NON-DEFAULT PERSPECTIVE AND GREY SURFACE LINES
  d3 <- plot(cp,xlab=("PM10 (??g/m³)"),zlab=paste("NAD RR   ", sub$regnames[1], sep=""),phi=20,theta=205,ltheta=170,
             shade=0.4, cex.lab=1.5)
  
  # LINES IN THE SURFACE CORRESPONDING TO THE EFFECTS IN THE PLOTS BELOW
  lines(trans3d(x=100,y=0:10,z=cp$matRRfit[as.character(96),],
                pmat=d3),lwd=2,col=2)
  lines(trans3d(x=cp$predvar,y=1,z=cp$matRRfit[,"lag1"],
                pmat=d3),lwd=2,col=2)
  
  par(mar=c(5,4,1,1)+0.1,mgp=c(2.5,1,0))
  
  # PLOTS FOR PREDICTOR-SPECIFIC, LAG-SPECIFIC AND OVERALL CUMULATIVE SUMMARIES
  
  plot(crvar,xlab="Lag",ylab="RR", ci.level = 0.9, col=2,lwd=2)
  mtext(text=paste("Predictor-specific association at PM10 concentration (",100,
                   "??g/m³)",sep=""))
  plot(crlag,xlab="PM10 (??g/m³)",ylab="RR",ci.level = 0.9, col=2,ylim=c(.8,1.3),lwd=2)
  mtext(text="Lag-specific association at lag 1")
  plot(crall,xlab="PM10 (??g/m³)",ylab="RR",ci.level = 0.9, ylim=c(.8,1.5),col=2,lwd=2)
  mtext(text="Overall cumulative association for 1 lag")
  
  
}

##############################################################################################
# Regional RR CALCULATIONS AT DIFFERENT DISCRETE LAG DAYS AND OVERALL CUMULATIVE ASSOCIATIONS
##############################################################################################

for(i in seq(data)) {
  
  # LOAD
  sub <- data[[i]]
  
  suppressWarnings(
    cb <- crossbasis(sub$PM10,lag=lag,argvar=argvar,arglag=arglag)
    
  )
  
  # Creating Cross Bases for Moving Averages
  suppressWarnings(
    cb1 <- crossbasis(sub$PM10ma1,lag=lag,argvar=argvar,arglag=arglag)
    
  )
  suppressWarnings(
    cb2 <- crossbasis(sub$PM10ma2,lag=lag,argvar=argvar,arglag=arglag)
    
  )
  suppressWarnings(
    cb3 <- crossbasis(sub$PM10ma3,lag=lag,argvar=argvar,arglag=arglag)
    
  )
  # RUN THE MODELS
  
  if(length(data)== 8){
    
  if(i == 2|| i==3){
    
    model <- glm(cd ~ cb + ns(rh, 4) + ns(temp, 4)+ ns(O3, 4) + dow + ns(time,df=10*8), family=quasipoisson(),sub)
    model1 <- glm(cd ~ cb1 + ns(rh, 4) + ns(temp, 4)+ ns(O3, 4) + dow + ns(time,df=10*8), family=quasipoisson(),sub)
    model2 <- glm(cd ~ cb2 + ns(rh, 4) + ns(temp, 4)+ ns(O3, 4) + dow + ns(time,df=10*8), family=quasipoisson(),sub)
    model3 <- glm(cd ~ cb3 + ns(rh, 4) + ns(temp, 4)+ ns(O3, 4) + dow + ns(time,df=10*8), family=quasipoisson(),sub)
    
  }else if (i == 5|| i==6){
    
    model <- glm(cd ~ cb + ns(NOX,4) + ns(rh, 4) + ns(temp, 4)+ ns(O3, 4) +ns(SO2,4)+ dow + ns(time,df=10*7), family=quasipoisson(),sub)
    model1 <- glm(cd ~ cb1 + ns(NOX,4) + ns(rh, 4) + ns(temp, 4)+ ns(O3, 4) +ns(SO2,4)+ dow + ns(time,df=10*7), family=quasipoisson(),sub)
    model2 <- glm(cd ~ cb2 + ns(NOX,4) + ns(rh, 4) + ns(temp, 4)+ ns(O3, 4) +ns(SO2,4)+ dow + ns(time,df=10*7), family=quasipoisson(),sub)
    model3 <- glm(cd ~ cb3 + ns(NOX,4) + ns(rh, 4) + ns(temp, 4)+ ns(O3, 4) +ns(SO2,4)+ dow + ns(time,df=10*7), family=quasipoisson(),sub)
    
  }else if(i == 7|| i==8){
    
    model <- glm(cd ~ cb + ns(NOX,4) + ns(rh, 4) + ns(temp, 4)+ ns(O3, 4) +ns(SO2,4)+ dow + ns(time,df=10*6), family=quasipoisson(),sub)
    model1 <- glm(cd ~ cb1 + ns(NOX,4) + ns(rh, 4) + ns(temp, 4)+ ns(O3, 4) +ns(SO2,4)+ dow + ns(time,df=10*6), family=quasipoisson(),sub)
    model2 <- glm(cd ~ cb2 + ns(NOX,4) + ns(rh, 4) + ns(temp, 4)+ ns(O3, 4) +ns(SO2,4)+ dow + ns(time,df=10*6), family=quasipoisson(),sub)
    model3 <- glm(cd ~ cb3 + ns(NOX,4) + ns(rh, 4) + ns(temp, 4)+ ns(O3, 4) +ns(SO2,4)+ dow + ns(time,df=10*6), family=quasipoisson(),sub)
    
  } else{
    
    model <- glm(cd ~ cb + ns(NOX,4) + ns(rh, 4) + ns(temp, 4)+ ns(O3, 4) +ns(SO2,4)+ dow + ns(time,df=10*10), family=quasipoisson(),sub)
    model1 <- glm(cd ~ cb1 + ns(NOX,4) + ns(rh, 4) + ns(temp, 4)+ ns(O3, 4) +ns(SO2,4)+ dow + ns(time,df=10*10), family=quasipoisson(),sub)
    model2 <- glm(cd ~ cb2+ ns(NOX,4) + ns(rh, 4) + ns(temp, 4)+ ns(O3, 4) +ns(SO2,4)+ dow + ns(time,df=10*10), family=quasipoisson(),sub)
    model3 <- glm(cd ~ cb3 + ns(NOX,4) + ns(rh, 4) + ns(temp, 4)+ ns(O3, 4) +ns(SO2,4)+ dow + ns(time,df=10*10), family=quasipoisson(),sub)
    
  }
  }else{
    
    model <- glm(cd ~ cb + ns(NOX,4) + ns(rh, 4) + ns(temp, 4)+ ns(O3, 4) +ns(SO2,4)+ dow + ns(time,df=10*10), family=quasipoisson(),sub)
    model1 <- glm(cd ~ cb1 + ns(NOX,4) + ns(rh, 4) + ns(temp, 4)+ ns(O3, 4) +ns(SO2,4)+ dow + ns(time,df=10*10), family=quasipoisson(),sub)
    model2 <- glm(cd ~ cb2 + ns(NOX,4) + ns(rh, 4) + ns(temp, 4)+ ns(O3, 4) +ns(SO2,4)+ dow + ns(time,df=10*10), family=quasipoisson(),sub)
    model3 <- glm(cd ~ cb3 + ns(NOX,4) + ns(rh, 4) + ns(temp, 4)+ ns(O3, 4) +ns(SO2,4)+ dow + ns(time,df=10*10), family=quasipoisson(),sub)
    
  }
  

# Cumulative Overall:
    
crall <- crossreduce(cb,model, type = "overall", lag =lag, from= bound[1], to= bound[2], cen=10)
    
crall1 <- crossreduce(cb,model,type = "overall", lag = lag1,from= bound[1], to= bound[2],cen=10)

crall2 <- crossreduce(cb,model,type = "overall", lag = lag2,from= bound[1], to= bound[2],cen=10)

crall3 <- crossreduce(cb,model,type = "overall", lag = lag3,from= bound[1], to= bound[2],cen=10)


# Predictor Specific:

crlag0 <- crossreduce(cb,model,type="lag",value=0,from=bound[1],to=bound[2],
                     bylag=0.2,cen=10)

crlag1 <- crossreduce(cb,model,type="lag",value=1,from=bound[1],to=bound[2],
                     bylag=0.2,cen=10)

crlag2 <- crossreduce(cb,model,type="lag",value=2,from=bound[1],to=bound[2],
                     bylag=0.2,cen=10)

crlag3 <- crossreduce(cb,model,type="lag",value=3,from=bound[1],to=bound[2],
                     bylag=0.2,cen=10)

# Predictor Specific, MOVING AVERAGES:

crmov1 <- crossreduce(cb1,model1,type="lag",value=0,from=bound[1],to=bound[2],
                      bylag=0.2,cen=10)

crmov2 <- crossreduce(cb2,model2,type="lag",value=0,from=bound[1],to=bound[2],
                      bylag=0.2,cen=10)

crmov3 <- crossreduce(cb3,model3,type="lag",value=0,from=bound[1],to=bound[2],
                      bylag=0.2,cen=10)


# Overall Cumulative Plots

par(mfrow=c(2,2))

plot(crall,ylab="RR",ylim=c(.8,1.5),xlab="PM10 (??g/m³)")
mtext(text=paste("Overall Cumulative association at for lag period 0-10 ", sub$regnames[1]),cex=0.7)

plot(crall1,ylab="RR",ylim=c(.8,1.5),xlab="PM10 (??g/m³)")
mtext(text="Overall Cumulative association at for lag period 0-1",cex=0.7)

plot(crall2,ylab="RR",ylim=c(.8,1.5),xlab="PM10 (??g/m³)")
mtext(text="Overall Cumulative association at for lag period 0-2",cex=0.7)

plot(crall3,ylab="RR",ylim=c(.8,1.5),xlab="PM10 (??g/m³)")
mtext(text="Overall Cumulative association at for lag period 0-3",cex=0.7)


# Predictor Specific Plots
par(mfrow=c(2,2))

plot(crlag0,xlab="PM10 (??g/m³)",ylab="RR", col=2,ylim=c(.8,1.3),lwd=2, cex.lab= 1.5)
mtext(text=paste("Lag-specific association at lag 0"))

plot(crlag1,xlab="PM10 (??g/m³)",ylab="RR", col=2,ylim=c(.8,1.3),lwd=2, cex.lab= 1.5)
mtext(text="Lag-specific association at lag 1")

plot(crlag2,xlab="PM10 (??g/m³)",ylab="RR", col=2,ylim=c(.8,1.3),lwd=2, cex.lab= 1.5)
mtext(text="Lag-specific association at lag 2")

plot(crlag3,xlab="PM10 (??g/m³)",ylab="RR", col=2,ylim=c(.8,1.3),lwd=2, cex.lab= 1.5)
mtext(text="Lag-specific association at lag 3")


# Predictor Specific Plots, MOVING AVERAGES

par(mfrow=c(2,2))

plot(crmov1,xlab="PM10 (??g/m³)",ylab="RR", col=2,ylim=c(.8,1.3),lwd=2, cex.lab= 1.5)
mtext(text=paste("Moving average association for 0-1"))

plot(crmov2,xlab="PM10 (??g/m³)",ylab="RR", col=2,ylim=c(.8,1.3),lwd=2, cex.lab= 1.5)
mtext(text="Moving average association for 0-2")

plot(crmov3,xlab="PM10 (??g/m³)",ylab="RR", col=2,ylim=c(.8,1.3),lwd=2, cex.lab= 1.5)
mtext(text="Moving average association for 0-3")


# DEFINING THE ARRAYS THAT WILL BE USED IN THE RR CALCULATION

# Cumulative:

array1.1 = vector(mode = "numeric", length= 16)
array2.1 = vector(mode = "numeric", length= 16)
array3.1 = vector(mode = "numeric", length= 16)

array1.2 = vector(mode = "numeric", length= 16)
array2.2 = vector(mode = "numeric", length= 16)
array3.2 = vector(mode = "numeric", length= 16)

array1.3 = vector(mode = "numeric", length= 16)
array2.3 = vector(mode = "numeric", length= 16)
array3.3 = vector(mode = "numeric", length= 16)

# Lag specific:

array1.4 = vector(mode = "numeric", length= 16)
array2.4 = vector(mode = "numeric", length= 16)
array3.4 = vector(mode = "numeric", length= 16)

array1.5 = vector(mode = "numeric", length= 16)
array2.5 = vector(mode = "numeric", length= 16)
array3.5 = vector(mode = "numeric", length= 16)

array1.6 = vector(mode = "numeric", length= 16)
array2.6 = vector(mode = "numeric", length= 16)
array3.6 = vector(mode = "numeric", length= 16)

array1.7 = vector(mode = "numeric", length= 16)
array2.7 = vector(mode = "numeric", length= 16)
array3.7 = vector(mode = "numeric", length= 16)


# Moving averages:
array1.8 = vector(mode = "numeric", length= 16)
array2.8 = vector(mode = "numeric", length= 16)
array3.8 = vector(mode = "numeric", length= 16)

array1.9 = vector(mode = "numeric", length= 16)
array2.9 = vector(mode = "numeric", length= 16)
array3.9 = vector(mode = "numeric", length= 16)

array1.11 = vector(mode = "numeric", length= 16)
array2.11 = vector(mode = "numeric", length= 16)
array3.11 = vector(mode = "numeric", length= 16)



# GOING THROUGH THE RANGE OF 10-160 [PM10] TO FIND THE RR AT 10 UNIT INCREMENTS

for(i in 1:16){
  
  j = (i*10)
  j.1= as.character(j)
  
  # Cumulative:
  
  array1.1[i]= round(with(crall1,cbind(RRlow)[j.1,]),6)
  array2.1[i]= round(with(crall1,cbind(RRfit)[j.1,]),6)
  array3.1[i]= round(with(crall1,cbind(RRhigh)[j.1,]),6)
  
  array1.2[i]= round(with(crall2,cbind(RRlow)[j.1,]),6)
  array2.2[i]= round(with(crall2,cbind(RRfit)[j.1,]),6)
  array3.2[i]= round(with(crall2,cbind(RRhigh)[j.1,]),6)
  
  array1.3[i]= round(with(crall3,cbind(RRlow)[j.1,]),6)
  array2.3[i]= round(with(crall3,cbind(RRfit)[j.1,]),6)
  array3.3[i]= round(with(crall3,cbind(RRhigh)[j.1,]),6)
  
  # Lag specific:
  
  array1.4[i]= round(with(crlag0,cbind(RRlow)[j.1,]),6)
  array2.4[i]= round(with(crlag0,cbind(RRfit)[j.1,]),6)
  array3.4[i]= round(with(crlag0,cbind(RRhigh)[j.1,]),6)
  

  array1.5[i]= round(with(crlag1,cbind(RRlow)[j.1,]),6)
  array2.5[i]= round(with(crlag1,cbind(RRfit)[j.1,]),6)
  array3.5[i]= round(with(crlag1,cbind(RRhigh)[j.1,]),6)
  
  
  array1.6[i]= round(with(crlag2,cbind(RRlow)[j.1,]),6)
  array2.6[i]= round(with(crlag2,cbind(RRfit)[j.1,]),6)
  array3.6[i]= round(with(crlag2,cbind(RRhigh)[j.1,]),6)
  
  array1.7[i]= round(with(crlag3,cbind(RRlow)[j.1,]),6)
  array2.7[i]= round(with(crlag3,cbind(RRfit)[j.1,]),6)
  array3.7[i]= round(with(crlag3,cbind(RRhigh)[j.1,]),6)
  
  # Moving averages:
  
  array1.8[i]= round(with(crmov1,cbind(RRlow)[j.1,]),6)
  array2.8[i]= round(with(crmov1,cbind(RRfit)[j.1,]),6)
  array3.8[i]= round(with(crmov1,cbind(RRhigh)[j.1,]),6)
  
  array1.9[i]= round(with(crmov2,cbind(RRlow)[j.1,]),6)
  array2.9[i]= round(with(crmov2,cbind(RRfit)[j.1,]),6)
  array3.9[i]= round(with(crmov2,cbind(RRhigh)[j.1,]),6)
  
  array1.11[i]= round(with(crmov3,cbind(RRlow)[j.1,]),6)
  array2.11[i]= round(with(crmov3,cbind(RRfit)[j.1,]),6)
  array3.11[i]= round(with(crmov3,cbind(RRhigh)[j.1,]),6)
  
  }

# Cumulative:

array_difference1.1 = vector(mode = "numeric", length= 15)
array_difference2.1 = vector(mode = "numeric", length= 15)
array_difference3.1 = vector(mode = "numeric", length= 15)


array_difference1.2 = vector(mode = "numeric", length= 15)
array_difference2.2 = vector(mode = "numeric", length= 15)
array_difference3.2 = vector(mode = "numeric", length= 15)


array_difference1.3 = vector(mode = "numeric", length= 15)
array_difference2.3 = vector(mode = "numeric", length= 15)
array_difference3.3 = vector(mode = "numeric", length= 15)

# Lag specific:

array_difference1.4 = vector(mode = "numeric", length= 15)
array_difference2.4 = vector(mode = "numeric", length= 15)
array_difference3.4 = vector(mode = "numeric", length= 15)


array_difference1.5 = vector(mode = "numeric", length= 15)
array_difference2.5 = vector(mode = "numeric", length= 15)
array_difference3.5 = vector(mode = "numeric", length= 15)


array_difference1.6 = vector(mode = "numeric", length= 15)
array_difference2.6 = vector(mode = "numeric", length= 15)
array_difference3.6 = vector(mode = "numeric", length= 15)


array_difference1.7 = vector(mode = "numeric", length= 15)
array_difference2.7 = vector(mode = "numeric", length= 15)
array_difference3.7 = vector(mode = "numeric", length= 15)

# Moving averages:

array_difference1.8 = vector(mode = "numeric", length= 15)
array_difference2.8 = vector(mode = "numeric", length= 15)
array_difference3.8 = vector(mode = "numeric", length= 15)


array_difference1.9 = vector(mode = "numeric", length= 15)
array_difference2.9 = vector(mode = "numeric", length= 15)
array_difference3.9 = vector(mode = "numeric", length= 15)

array_difference1.11 = vector(mode = "numeric", length= 15)
array_difference2.11 = vector(mode = "numeric", length= 15)
array_difference3.11 = vector(mode = "numeric", length= 15)


# STORING THE CHANGE OF RR PER 10 UNIT INCREASE OF PM10

for(i in 1:15){
 
  # Cumulative:
  
  array_difference1.1[i] = (array1.1[i+1]-array1.1[i])  
  array_difference2.1[i] = (array2.1[i+1]-array2.1[i])  
  array_difference3.1[i] = (array3.1[i+1]-array3.1[i])  
 
  array_difference1.2[i] = (array1.2[i+1]-array1.2[i])  
  array_difference2.2[i] = (array2.2[i+1]-array2.2[i])  
  array_difference3.2[i] = (array3.2[i+1]-array3.2[i])  
  
  array_difference1.3[i] = (array1.3[i+1]-array1.3[i])  
  array_difference2.3[i] = (array2.3[i+1]-array2.3[i])  
  array_difference3.3[i] = (array3.3[i+1]-array3.3[i])  
    
  # Lag specific:
  
  array_difference1.4[i] = (array1.4[i+1]-array1.4[i])  
  array_difference2.4[i] = (array2.4[i+1]-array2.4[i])  
  array_difference3.4[i] = (array3.4[i+1]-array3.4[i])  

  array_difference1.5[i] = (array1.5[i+1]-array1.5[i])  
  array_difference2.5[i] = (array2.5[i+1]-array2.5[i])  
  array_difference3.5[i] = (array3.5[i+1]-array3.5[i])  
  
  array_difference1.6[i] = (array1.6[i+1]-array1.6[i])  
  array_difference2.6[i] = (array2.6[i+1]-array2.6[i])  
  array_difference3.6[i] = (array3.6[i+1]-array3.6[i])  
  
  array_difference1.7[i] = (array1.7[i+1]-array1.7[i])  
  array_difference2.7[i] = (array2.7[i+1]-array2.7[i])  
  array_difference3.7[i] = (array3.7[i+1]-array3.7[i])  
  
  
  # Moving averages:
  array_difference1.8[i] = (array1.8[i+1]-array1.8[i])  
  array_difference2.8[i] = (array2.8[i+1]-array2.8[i])  
  array_difference3.8[i] = (array3.8[i+1]-array3.8[i])  
  
  array_difference1.9[i] = (array1.9[i+1]-array1.9[i])  
  array_difference2.9[i] = (array2.9[i+1]-array2.9[i])  
  array_difference3.9[i] = (array3.9[i+1]-array3.9[i])  
  
  array_difference1.11[i] = (array1.11[i+1]-array1.11[i])  
  array_difference2.11[i] = (array2.11[i+1]-array2.11[i])  
  array_difference3.11[i] = (array3.11[i+1]-array3.11[i])  
  
  }


# USING THE EQUATION RR = e^(BETA*DELTA PM10), WE CONVERT THE MEAN CHANGE
# IN RR PER 10 UNIT INCREASE OF PM10 TO BETA (THE HEALTH EFFECT ESTIMATE)
# THIS IS DONE FOR THE UPPER AND LOWER CONFIDENCE INTERVALS AS WELL

print(sub$regnames[1])

cat("\n")

print("Cumulative 0-1:")
print(paste("RRlow:", mean(array_difference1.1)+1 , " Beta:", log(mean(array_difference1.1)+1)/10))
print(paste("RRfit:", mean(array_difference2.1)+1, " Beta:",log(mean(array_difference2.1)+1)/10))
print(paste("RRhigh:", mean(array_difference3.1)+1, " Beta:", log(mean(array_difference3.1)+1)/10))
cat("\n")

print("Cumulative 0-2:")
print(paste("RRlow:", mean(array_difference1.2)+1, " Beta:", log(mean(array_difference1.2)+1)/10))
print(paste("RRfit:", mean(array_difference2.2)+1, " Beta:", log(mean(array_difference2.2)+1)/10))
print(paste("RRhigh:", mean(array_difference3.2)+1, " Beta:", log(mean(array_difference3.2)+1)/10))
cat("\n")

print("Cumulative 0-3:")
print(paste("RRlow:", mean(array_difference1.3)+1, " Beta:", log(mean(array_difference1.3)+1)/10))
print(paste("RRfit:", mean(array_difference2.3)+1, " Beta:", log(mean(array_difference2.3)+1)/10))
print(paste("RRhigh:", mean(array_difference3.3)+1, " Beta:", log(mean(array_difference3.3)+1)/10))

cat("\n")
cat("\n")

print("Lag 0:")
print(paste("RRlow:", mean(array_difference1.4)+1 , " Beta:", log(mean(array_difference1.4)+1)/10))
print(paste("RRfit:", mean(array_difference2.4)+1, " Beta:",log(mean(array_difference2.4)+1)/10))
print(paste("RRhigh:", mean(array_difference3.4)+1, " Beta:", log(mean(array_difference3.4)+1)/10))

cat("\n")

print("Lag 1:")
print(paste("RRlow:", mean(array_difference1.5)+1 , " Beta:", log(mean(array_difference1.5)+1)/10))
print(paste("RRfit:", mean(array_difference2.5)+1, " Beta:",log(mean(array_difference2.5)+1)/10))
print(paste("RRhigh:", mean(array_difference3.5)+1, " Beta:", log(mean(array_difference3.5)+1)/10))
cat("\n")

print("Lag 2:")
print(paste("RRlow:", mean(array_difference1.6)+1 , " Beta:", log(mean(array_difference1.6)+1)/10))
print(paste("RRfit:", mean(array_difference2.6)+1, " Beta:",log(mean(array_difference2.6)+1)/10))
print(paste("RRhigh:", mean(array_difference3.6)+1, " Beta:", log(mean(array_difference3.6)+1)/10))
cat("\n")

print("Lag 3:")
print(paste("RRlow:", mean(array_difference1.7)+1 , " Beta:", log(mean(array_difference1.7)+1)/10))
print(paste("RRfit:", mean(array_difference2.7)+1, " Beta:",log(mean(array_difference2.7)+1)/10))
print(paste("RRhigh:", mean(array_difference3.7)+1, " Beta:", log(mean(array_difference3.7)+1)/10))
cat("\n")
cat("\n")

print("Moving 0-1:")
print(paste("RRlow:", mean(array_difference1.8)+1 , " Beta:", log(mean(array_difference1.8)+1)/10))
print(paste("RRfit:", mean(array_difference2.8)+1, " Beta:",log(mean(array_difference2.8)+1)/10))
print(paste("RRhigh:", mean(array_difference3.8)+1, " Beta:", log(mean(array_difference3.8)+1)/10))
cat("\n")

print("Moving 0-2:")
print(paste("RRlow:", mean(array_difference1.9)+1 , " Beta:", log(mean(array_difference1.9)+1)/10))
print(paste("RRfit:", mean(array_difference2.9)+1, " Beta:",log(mean(array_difference2.9)+1)/10))
print(paste("RRhigh:", mean(array_difference3.9)+1, " Beta:", log(mean(array_difference3.9)+1)/10))
cat("\n")

print("Moving 0-3:")
print(paste("RRlow:", mean(array_difference1.11)+1 , " Beta:", log(mean(array_difference1.11)+1)/10))
print(paste("RRfit:", mean(array_difference2.11)+1, " Beta:",log(mean(array_difference2.11)+1)/10))
print(paste("RRhigh:", mean(array_difference3.11)+1, " Beta:", log(mean(array_difference3.11)+1)/10))
cat("\n")

}
