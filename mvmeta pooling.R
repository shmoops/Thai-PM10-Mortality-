###############################################################################
#   Code originally for the analysis in:
#
#   "Reducing and meta-analyzing estimates from distributed lag non-linear models"
#   Gasparrini and Armstrong 
#   BMC Medical Research Methodology - 2013
#   http://www.ag-myresearch.com/2013_gasparrini_bmcmrm.html
#
# *modified to suite the purposes of our research (particulate matter)
###############################################################################

# LOAD PACKAGES (ASSUMED ALREADY INSTALLED)
library(dlnm) ; library(mvmeta) ; library(splines)

# REGIONS
regions <- as.character(unique(mvmeta_MONSTER_central$regnames))

# CREATE A LIST WITH THE REGIONAL SERIES
data <- lapply(regions,function(x) mvmeta_MONSTER_central[mvmeta_MONSTER_central$regnames==x,])
names(data) <- regions
m <- length(regions)

# ARGUMENTS AND LISTS FOR CROSS-BASIS DEFINITION
bound <- c(10, 160)
varknots <- equalknots(bound, fun = "bs", degree = 2,df=6)
argvar <- list(fun="bs",degree=2,knots=varknots,bound=bound)

lagknots <- logknots(10,df=5, intercept = TRUE)
arglag <- list(fun="ns", knots=lagknots)

# DIFFERENT MAX LAGS
lag <- c(0,10)
lag1 <- c(0,1)
lag2 <- c(0,2)
lag3 <- c(0,3) 

####################################################################
# BUILT OBJECTS WHERE RESULTS WILL BE STORED
#   y- IS THE MATRIX FOR THE OUTCOME PARAMETERS
#   S- IS THE LISTS OF (CO)VARIANCE MATRICES

# OVERALL CUMULATIVE AND PREDICTOR-SPECIFIC SUMMARIES
yall <- matrix(NA,length(data),6,dimnames=list(regions,paste("b",seq(6),sep="")))
ymov1 <- ymov2 <- ymov3 <- ylag0 <- ylag1 <-ylag2 <-ylag3 <-yall1 <-yall2 <-yall3 <-yall


# (CO)VARIANCE MATRICES
Sall <- vector("list",length(data))
names(Sall) <- regions
Smov1 <- Smov2 <- Smov3 <- Slag0 <- Slag1<- Slag2<- Slag3<- Sall1<- Sall2<- Sall3<- Sall

####################################################################
# RUN THE MODEL FOR EACH CITY

# WARNING FOR PREDICTION BEYOND BOUNDARIES SUPPRESSED
options(warn=-1)

# LOOP FOR CITIES
for(i in seq(data)) {
  
  # PRINT
  cat(i,"")
  
  # LOAD
  sub <- data[[i]]
  
  # DEFINE THE CROSS-BASES
  cb <- crossbasis(sub$PM10,lag=lag,argvar=argvar,arglag=arglag)
  cb1 <- crossbasis(sub$PM10ma1,lag=lag,argvar=argvar,arglag=arglag)
  cb2 <- crossbasis(sub$PM10ma2,lag=lag,argvar=argvar,arglag=arglag)
  cb3 <- crossbasis(sub$PM10ma3,lag=lag,argvar=argvar,arglag=arglag)
  
  # RUN THE FIRST-STAGE MODELS
  if(length(data)== 8){
    
    if(i == 2|| i==3){
      
      mfirst <- glm(nad ~ cb + ns(rh, 4) + ns(temp, 4)+ ns(O3, 4) + dow + ns(time,df=10*8), family=quasipoisson(),sub)
      mfirst1 <- glm(nad ~ cb1 + ns(rh, 4) + ns(temp, 4)+ ns(O3, 4) + dow + ns(time,df=10*8), family=quasipoisson(),sub)
      mfirst2 <- glm(nad ~ cb2 + ns(rh, 4) + ns(temp, 4)+ ns(O3, 4) + dow + ns(time,df=10*8), family=quasipoisson(),sub)
      mfirst3 <- glm(nad ~ cb3 + ns(rh, 4) + ns(temp, 4)+ ns(O3, 4) + dow + ns(time,df=10*8), family=quasipoisson(),sub)
      
    }else if (i == 5|| i==6){
      
      mfirst <- glm(nad ~ cb + ns(NOX,4) + ns(rh, 4) + ns(temp, 4)+ ns(O3, 4) +ns(SO2,4)+ dow + ns(time,df=10*7), family=quasipoisson(),sub)
      mfirst1 <- glm(nad ~ cb1 + ns(NOX,4) + ns(rh, 4) + ns(temp, 4)+ ns(O3, 4) +ns(SO2,4)+ dow + ns(time,df=10*7), family=quasipoisson(),sub)
      mfirst2 <- glm(nad ~ cb2 + ns(NOX,4) + ns(rh, 4) + ns(temp, 4)+ ns(O3, 4) +ns(SO2,4)+ dow + ns(time,df=10*7), family=quasipoisson(),sub)
      mfirst3 <- glm(nad ~ cb3 + ns(NOX,4) + ns(rh, 4) + ns(temp, 4)+ ns(O3, 4) +ns(SO2,4)+ dow + ns(time,df=10*7), family=quasipoisson(),sub)
      
    }else if(i == 7|| i==8){
      
      mfirst <- glm(nad ~ cb + ns(NOX,4) + ns(rh, 4) + ns(temp, 4)+ ns(O3, 4) +ns(SO2,4)+ dow + ns(time,df=10*6), family=quasipoisson(),sub)
      mfirst1 <- glm(nad ~ cb1 + ns(NOX,4) + ns(rh, 4) + ns(temp, 4)+ ns(O3, 4) +ns(SO2,4)+ dow + ns(time,df=10*6), family=quasipoisson(),sub)
      mfirst2 <- glm(nad ~ cb2 + ns(NOX,4) + ns(rh, 4) + ns(temp, 4)+ ns(O3, 4) +ns(SO2,4)+ dow + ns(time,df=10*6), family=quasipoisson(),sub)
      mfirst3 <- glm(nad ~ cb3 + ns(NOX,4) + ns(rh, 4) + ns(temp, 4)+ ns(O3, 4) +ns(SO2,4)+ dow + ns(time,df=10*6), family=quasipoisson(),sub)
      
    } else{
      
      mfirst <- glm(nad ~ cb + ns(NOX,4) + ns(rh, 4) + ns(temp, 4)+ ns(O3, 4) +ns(SO2,4)+ dow + ns(time,df=10*10), family=quasipoisson(),sub)
      mfirst1 <- glm(nad ~ cb1 + ns(NOX,4) + ns(rh, 4) + ns(temp, 4)+ ns(O3, 4) +ns(SO2,4)+ dow + ns(time,df=10*10), family=quasipoisson(),sub)
      mfirst2 <- glm(nad ~ cb2 + ns(NOX,4) + ns(rh, 4) + ns(temp, 4)+ ns(O3, 4) +ns(SO2,4)+ dow + ns(time,df=10*10), family=quasipoisson(),sub)
      mfirst3 <- glm(nad ~ cb3 + ns(NOX,4) + ns(rh, 4) + ns(temp, 4)+ ns(O3, 4) +ns(SO2,4)+ dow + ns(time,df=10*10), family=quasipoisson(),sub)
      
    }
  }else{
    
    mfirst <- glm(nad ~ cb + ns(NOX,4) + ns(rh, 4) + ns(temp, 4)+ ns(O3, 4) +ns(SO2,4)+ dow + ns(time,df=10*10), family=quasipoisson(),sub)
    mfirst1 <- glm(nad ~ cb1 + ns(NOX,4) + ns(rh, 4) + ns(temp, 4)+ ns(O3, 4) +ns(SO2,4)+ dow + ns(time,df=10*10), family=quasipoisson(),sub)
    mfirst2 <- glm(nad ~ cb2 + ns(NOX,4) + ns(rh, 4) + ns(temp, 4)+ ns(O3, 4) +ns(SO2,4)+ dow + ns(time,df=10*10), family=quasipoisson(),sub)
    mfirst3 <- glm(nad ~ cb3 + ns(NOX,4) + ns(rh, 4) + ns(temp, 4)+ ns(O3, 4) +ns(SO2,4)+ dow + ns(time,df=10*10), family=quasipoisson(),sub)
    
  }
  
  
  
  ####################################################################
  # REDUCTION TO SUMMARY ASSOCIATIONS
  
  # TO OVERALL CUMULATIVE SUMMARY
  # NB: CENTERING NOT REALLY NEEDED HERE, AS COEF-VCOV (EXTRACTED BELOW) IN THE
  #   VAR SPACE DO NOT DEPEND ON CENTERING VALUE
  crall <- crossreduce(cb, mfirst, type = "overall", lag = lag , from= bound[1], to= bound[2],cen=10)
  crall1 <- crossreduce(cb, mfirst, type = "overall", lag = lag1 , from= bound[1], to= bound[2],cen=10)
  crall2 <- crossreduce(cb, mfirst, type = "overall", lag = lag2 , from= bound[1], to= bound[2],cen=10)
  crall3 <- crossreduce(cb, mfirst, type = "overall", lag = lag3 , from= bound[1], to= bound[2],cen=10)
  
  # To lag-specific associations
  crlag0 <- crossreduce(cb, mfirst, type = "lag", val= 0, cen = 10)
  crlag1 <- crossreduce(cb, mfirst, type = "lag", val= 1, cen = 10)
  crlag2 <- crossreduce(cb, mfirst, type = "lag", val= 2, cen = 10)
  crlag3 <- crossreduce(cb, mfirst, type = "lag", val= 3, cen = 10)
  
  # TO MOVEING AVERAGE ASSOCIATIONS
 
  crmov1 <- crossreduce(cb1,mfirst1,type="lag",value=0,from=bound[1],to=bound[2],
                        bylag=0.2,cen=10)
  
  crmov2 <- crossreduce(cb2,mfirst2,type="lag",value=0,from=bound[1],to=bound[2],
                        bylag=0.2,cen=10)
  
  crmov3 <- crossreduce(cb3,mfirst3,type="lag",value=0,from=bound[1],to=bound[2],
                        bylag=0.2,cen=10)
  
  
  ####################################################################
  # STORE THE RESULTS
  
  # OVERALL CUMULATIVE SUMMARY FOR THE MAIN MODEL
  yall[i,] <- coef(crall)
  yall1[i,] <- coef(crall1)
  yall2[i,] <- coef(crall2)
  yall3[i,] <- coef(crall3)  
  
  Sall[[i]] <- vcov(crall)
  Sall1[[i]] <- vcov(crall1)
  Sall2[[i]] <- vcov(crall2)
  Sall3[[i]] <- vcov(crall3)
 
  # PREDICTOR SPECIFIC SUMMARY FOR THE MAIN MODEL
  
  ylag0[i,] <- coef(crlag0)
  ylag1[i,] <- coef(crlag1)
  ylag2[i,] <- coef(crlag2)
  ylag3[i,] <- coef(crlag3)

  Slag0[[i]] <- vcov(crlag0)
  Slag1[[i]] <- vcov(crlag1)
  Slag2[[i]] <- vcov(crlag2)
  Slag3[[i]] <- vcov(crlag3)
  
  # MOVING AVERAGE SUMMARY
  ymov1[i,] <- coef(crmov1)
  ymov2[i,] <- coef(crmov2)
  ymov3[i,] <- coef(crmov3)
  
  Smov1[[i]] <- vcov(crmov1)
  Smov2[[i]] <- vcov(crmov2)
  Smov3[[i]] <- vcov(crmov3)
  
}

####################################################################

# RESET WARNING
options(warn=0)

####################################################################
# SECOND STAGE
# - RUN THE MULTIVARIATE META-ANALYTICAL MODELS WITH mvmeta
# - CREATE BASIS VARIABLES USING onebasis, TO BE USED FOR PREDICTION
# - OBTAIN PREDICTIONS THROUGH crosspred (dlnm)
####################################################################

####################################################################
# PERFORM MULTIVARIATE META-ANALYSIS

# LOAD THE PACKAGES (mvmeta PACKAGE IS ASSUMED TO BE INSTALLED)
library(mvmeta)

# SELECT THE ESTIMATION METHOD
method <- "reml"
# IN THE CURRENT VERSION, SET control=list(showiter=T) TO 
#   INSPECT THE OPTIMIZATION SEARCH

# OVERALL CUMULATIVE SUMMARY FOR THE MAIN MODEL
mvall <- mvmeta(yall~1,Sall,method=method)
mvall1 <- mvmeta(yall1~1,Sall1,method=method)
mvall2 <- mvmeta(yall2~1,Sall2,method=method)
mvall3 <- mvmeta(yall3~1,Sall3,method=method)

summary(mvall)

# PREDICTOR-SPECIFIC SUMMARY FOR lags 0-3 (MAIN MODEL)
mvall_lag0 <- mvmeta(ylag0~1,Slag0, method = method)
mvall_lag1 <- mvmeta(ylag1~1,Slag1, method = method)
mvall_lag2 <- mvmeta(ylag2~1,Slag2, method = method)
mvall_lag3 <- mvmeta(ylag3~1,Slag3, method = method)


# MOVING AVERAGES SUMMARY 
mvall_mov1 <- mvmeta(ymov1~1,Smov1, method = method)
mvall_mov2<- mvmeta(ymov2~1,Smov2, method = method)
mvall_mov3 <- mvmeta(ymov3~1,Smov3, method = method)

# Q TEST AND I-SQUARE
(qall <- qtest(mvall_mov3))
round(((qall$Q-qall$df)/qall$Q)[1]*100,1)

####################################################################
# CREATE BASES FOR PREDICTION

# BASES OF PM10 AND LAG USED TO PREDICT, EQUAL TO THAT USED FOR ESTIMATION
# COMPUTED USING THE ATTRIBUTES OF THE CROSS-BASIS USED IN ESTIMATION
xvar <- seq(bound[1],bound[2],by=0.1)
bvar <- do.call("onebasis",c(list(x=xvar),attr(cb,"argvar")))
xlag <- 0:100/10
blag <- do.call("onebasis",c(list(x=xlag),attr(cb,"arglag")))

####################################################################
# provincial FIRST-STAGE SUMMARIES

# CUMULATIVE

regall <- lapply(seq(nrow(yall)),function(i) crosspred(bvar,coef=yall[i,],
                                                       vcov=Sall[[i]],model.link="log",cen=10))
regall1 <- lapply(seq(nrow(yall1)),function(i) crosspred(bvar,coef=yall1[i,],
                                                         vcov=Sall1[[i]],model.link="log",cen=10))
regall2 <- lapply(seq(nrow(yall2)),function(i) crosspred(bvar,coef=yall2[i,],
                                                         vcov=Sall2[[i]],model.link="log",cen=10))
regall3 <- lapply(seq(nrow(yall3)),function(i) crosspred(bvar,coef=yall3[i,],
                                                         vcov=Sall3[[i]],model.link="log",cen=10))
# DISCRETE LAGS

regall_lag0 <- lapply(seq(nrow(ylag0)),function(i) crosspred(bvar,coef=ylag0[i,],
                                                             vcov=Slag0[[i]],model.link="log",cen=10))
regall_lag1 <- lapply(seq(nrow(ylag1)),function(i) crosspred(bvar,coef=ylag1[i,],
                                                             vcov=Slag1[[i]],model.link="log",cen=10))
regall_lag2 <- lapply(seq(nrow(ylag2)),function(i) crosspred(bvar,coef=ylag2[i,],
                                                             vcov=Slag2[[i]],model.link="log",cen=10))
regall_lag3 <- lapply(seq(nrow(ylag3)),function(i) crosspred(bvar,coef=ylag3[i,],
                                                             vcov=Slag3[[i]],model.link="log",cen=10))

# MOVING AVERAGE

regall_mov1 <- lapply(seq(nrow(ymov1)),function(i) crosspred(bvar,coef=ymov1[i,],
                                                             vcov=Smov1[[i]],model.link="log",cen=10))
regall_mov2 <- lapply(seq(nrow(ymov2)),function(i) crosspred(bvar,coef=ymov2[i,],
                                                             vcov=Smov2[[i]],model.link="log",cen=10))
regall_mov3 <- lapply(seq(nrow(ymov3)),function(i) crosspred(bvar,coef=ymov3[i,],
                                                             vcov=Smov3[[i]],model.link="log",cen=10))

####################################################################
# PREDICTION FOR A GRID OF PM10 AND LAG VALUES

# OVERALL CUMULATIVE SUMMARY ASSOCIATION FOR MAIN MODEL
cpall <- crosspred(bvar,coef=coef(mvall),vcov=vcov(mvall),
                   model.link="log",by=0.1,from=bound[1],to=bound[2],cen=10)
cpall1 <- crosspred(bvar,coef=coef(mvall1),vcov=vcov(mvall1),
                    model.link="log",by=0.1,from=bound[1],to=bound[2],cen=10)
cpall2 <- crosspred(bvar,coef=coef(mvall2),vcov=vcov(mvall2),
                    model.link="log",by=0.1,from=bound[1],to=bound[2],cen=10)
cpall3 <- crosspred(bvar,coef=coef(mvall3),vcov=vcov(mvall3),
                    model.link="log",by=0.1,from=bound[1],to=bound[2],cen=10)


# PREDICTOR-SPECIFIC SUMMARIES FOR 0-3 lag (MAIN MODEL)

cpall_lag0 <- crosspred(bvar,coef=coef(mvall_lag0),vcov=vcov(mvall_lag0),
                        model.link="log",by=0.1,from=bound[1],to=bound[2],cen=10)

cpall_lag1 <- crosspred(bvar,coef=coef(mvall_lag1),vcov=vcov(mvall_lag1),
                        model.link="log",by=0.1,from=bound[1],to=bound[2],cen=10)

cpall_lag2 <- crosspred(bvar,coef=coef(mvall_lag2),vcov=vcov(mvall_lag2),
                        model.link="log",by=0.1,from=bound[1],to=bound[2],cen=10)

cpall_lag3 <- crosspred(bvar,coef=coef(mvall_lag3),vcov=vcov(mvall_lag3),
                        model.link="log",by=0.1,from=bound[1],to=bound[2],cen=10)

# MOVING AVERAGES

cpall_mov1 <- crosspred(bvar,coef=coef(mvall_mov1),vcov=vcov(mvall_mov1),
                        model.link="log",by=0.1,from=bound[1],to=bound[2],cen=10)

cpall_mov2 <- crosspred(bvar,coef=coef(mvall_mov2),vcov=vcov(mvall_mov2),
                        model.link="log",by=0.1,from=bound[1],to=bound[2],cen=10)

cpall_mov3 <- crosspred(bvar,coef=coef(mvall_mov3),vcov=vcov(mvall_mov3),
                        model.link="log",by=0.1,from=bound[1],to=bound[2],cen=10)


####################################################################
# PLOTS:
for (i in 1){
# OVERALL CUMULATIVE SUMMARY ASSOCIATION

# PLOTS:
par(mfrow=c(1,1))
plot(cpall,type="n",ylab="RR",ylim=c(.8,1.5),xlab="PM10 (??g/m³)", cex.lab=1.5)
for(i in seq(regall)) lines(regall[[i]],ptype="overall",col=grey(0.5),lty=2)
abline(h=1)
lines(cpall,col=2,lwd=2)
mtext("Overall Cumulative Pooled Estimates MAX LAG = 10",cex=1.5)
legend ("top",c("Pooled (with 95%CI)","Provincial"),
        lty=c(1,2),lwd=1.5,col=c(2,grey(0.7)),bty="n",inset=0.1)

plot(cpall1,type="n",ylab="RR",ylim=c(.8,1.5),xlab="PM10 (??g/m³)", cex.lab=1.5)
for(i in seq(regall1)) lines(regall1[[i]],ptype="overall",col=grey(0.5),lty=2)
abline(h=1)
lines(cpall1,col=2,lwd=2)
mtext("Overall Cumulative Pooled Estimates MAX LAG = 1",cex=1.5)
legend ("top",c("Pooled (with 95%CI)","Provincial"),
        lty=c(1,2),lwd=1.5,col=c(2,grey(0.7)),bty="n",inset=0.1)

plot(cpall2,type="n",ylab="RR",ylim=c(.8,1.5),xlab="PM10 (??g/m³)", cex.lab=1.5)
for(i in seq(regall2)) lines(regall2[[i]],ptype="overall",col=grey(0.5),lty=2)
abline(h=1)
lines(cpall2,col=2,lwd=2)
mtext("Overall Cumulative Pooled Estimates MAX LAG = 2",cex=1.5)
legend ("top",c("Pooled (with 95%CI)","Provincial"),
        lty=c(1,2),lwd=1.5,col=c(2,grey(0.7)),bty="n",inset=0.1)

plot(cpall3,type="n",ylab="RR",ylim=c(.8,1.5),xlab="PM10 (??g/m³)", cex.lab=1.5)
for(i in seq(regall3)) lines(regall3[[i]],ptype="overall",col=grey(0.5),lty=2)
abline(h=1)
lines(cpall3,col=2,lwd=2)
mtext("Overall Cumulative Pooled Estimates MAX LAG = 3",cex=1.5)
legend ("top",c("Pooled (with 95%CI)","Provincial"),
        lty=c(1,2),lwd=1.5,col=c(2,grey(0.7)),bty="n",inset=0.1)


# POINTS OF MINIMUM MORTALITY
cpall$predvar[which.min(cpall$allRRfit)]
round(sum(mvmeta_MONSTER_central$PM10<10, na.rm = TRUE)/nrow(mvmeta_MONSTER_central)*100,1)


round(with(cpall,cbind(allRRfit,allRRlow,allRRhigh)["45",]),3)
round(with(cpall,cbind(allRRfit,allRRlow,allRRhigh)["55",]),3)



####################################################################

# PREDICTOR-SPECIFIC SUMMARIES

# PLOTS:
par(mfrow= c(1,1))
plot(cpall_lag0,type="n",ylab="Non-accidental RR",ylim=c(.95,1.12),xlab="PM10(??g/m³)", cex.lab=1.5, cex.axis= 1.5)
for(i in seq(regall_lag0)) lines(regall_lag0[[i]],ptype="overall",col=grey(0.5),lty=2)
abline(h=1)
lines(cpall_lag0,col=1,lwd=2)
legend ("top",c("Pooled (with 95%CI)","Provincial Relationships"),
        lty=c(1,2),lwd=1.5,col=c(1,grey(0.7)),bty="n",inset=0.1)
# mtext(text=paste("Predictor-specific summary for Lag = (",0,
#                  "days)",sep=""),cex=1)

plot(cpall_lag1,type="n",ylab="RR",ylim=c(.95,1.12),xlab="PM10(??g/m³)",cex.lab=1.5, cex.axis= 1.5, font.axis = 5)
for(i in seq(regall_lag1)) lines(regall_lag1[[i]],ptype="overall",col=grey(0.5),lty=2)
abline(h=1)
lines(cpall_lag1,col=1,lwd=2)
legend ("topleft",c("Pooled (with 95%CI)","provincial"),
        lty=c(1,2),lwd=2,col=c(1,grey(0.7)),bty="n",inset=0.1, cex = 1.5)
# mtext(text=paste("Predictor-specific summary for Lag = (",1,
#                  "day)",sep=""),cex=1)

plot(cpall_lag2,type="n",ylab="RR",ylim=c(.95,1.12),xlab="PM10(??g/m³)",cex.lab=1.5)
for(i in seq(regall_lag2)) lines(regall_lag2[[i]],ptype="overall",col=grey(0.5),lty=2)
abline(h=1)
lines(cpall_lag2,col=2,lwd=2)
legend ("top",c("Pooled (with 95%CI)","provincial"),
        lty=c(1,2),lwd=1.5,col=c(2,grey(0.7)),bty="n",inset=0.1)
mtext(text=paste("Predictor-specific summary for Lag = (",2,
                 "days)",sep=""),cex=1)

plot(cpall_lag3,type="n",ylab="RR",ylim=c(.95,1.12),xlab="PM10(??g/m³)",cex.lab=1.5)
for(i in seq(regall_lag3)) lines(regall_lag3[[i]],ptype="overall",col=grey(0.5),lty=2)
abline(h=1)
lines(cpall_lag3,col=2,lwd=2)
legend ("top",c("Pooled (with 95%CI)","provincial"),
        lty=c(1,2),lwd=1.5,col=c(2,grey(0.7)),bty="n",inset=0.1)
mtext(text=paste("Predictor-specific summary for Lag = (",3,
                 "days)",sep=""),cex=1)


##############################################################

# MOVING AVERAGES:

plot(cpall_mov1,type="n",ylab="RR",ylim=c(.88,1.3),xlab="PM10(??g/m³)",cex.lab=1.5)
for(i in seq(regall_mov1)) lines(regall_mov1[[i]],ptype="overall",col=grey(0.5),lty=2)
abline(h=1)
lines(cpall_mov1,col=2,lwd=2)
legend ("top",c("Pooled (with 95%CI)","provincial"),
        lty=c(1,2),lwd=1.5,col=c(2,grey(0.7)),bty="n",inset=0.1)
mtext(text=("Moving average 0-1"),cex=1)

plot(cpall_mov2,type="n",ylab="RR",ylim=c(.88,1.3),xlab="PM10(??g/m³)",cex.lab=1.5)
for(i in seq(regall_mov2)) lines(regall_mov2[[i]],ptype="overall",col=grey(0.5),lty=2)
abline(h=1)
lines(cpall_mov2,col=2,lwd=2)
legend ("top",c("Pooled (with 95%CI)","provincial"),
        lty=c(1,2),lwd=1.5,col=c(2,grey(0.7)),bty="n",inset=0.1)
mtext(text=("Moving average 0-2"),cex=1)

plot(cpall_mov3,type="n",ylab="RR",ylim=c(.85,1.35),xlab="PM10(??g/m³)",cex.lab=1.5)
for(i in seq(regall_mov3)) lines(regall_mov3[[i]],ptype="overall",col=grey(0.5),lty=2)
abline(h=1)
lines(cpall_mov3,col=2,lwd=2)
legend ("top",c("Pooled (with 95%CI)","provincial"),
        lty=c(1,2),lwd=1.5,col=c(2,grey(0.7)),bty="n",inset=0.1)
mtext(text=("Moving average 0-3"),cex=1)


}
#############################################################
# Beta Calculations:
#############################################################

for(i in 1){

# DEFINING THE ARRAYS THAT WILL BE USED IN THE RR CALCULATION

# Cumulative:

array1.1 = vector(mode = "numeric", length= 7)
array2.1 = vector(mode = "numeric", length= 7)
array3.1 = vector(mode = "numeric", length= 7)

array1.2 = vector(mode = "numeric", length= 7)
array2.2 = vector(mode = "numeric", length= 7)
array3.2 = vector(mode = "numeric", length= 7)

array1.3 = vector(mode = "numeric", length= 7)
array2.3 = vector(mode = "numeric", length= 7)
array3.3 = vector(mode = "numeric", length= 7)

# Lag specific:

array1.4 = vector(mode = "numeric", length= 7)
array2.4 = vector(mode = "numeric", length= 7)
array3.4 = vector(mode = "numeric", length= 7)

array1.5 = vector(mode = "numeric", length= 7)
array2.5 = vector(mode = "numeric", length= 7)
array3.5 = vector(mode = "numeric", length= 7)

array1.6 = vector(mode = "numeric", length= 7)
array2.6 = vector(mode = "numeric", length= 7)
array3.6 = vector(mode = "numeric", length= 7)

array1.7 = vector(mode = "numeric", length= 7)
array2.7 = vector(mode = "numeric", length= 7)
array3.7 = vector(mode = "numeric", length= 7)


# Moving averages:
array1.8 = vector(mode = "numeric", length= 7)
array2.8 = vector(mode = "numeric", length= 7)
array3.8 = vector(mode = "numeric", length= 7)

array1.9 = vector(mode = "numeric", length= 7)
array2.9 = vector(mode = "numeric", length= 7)
array3.9 = vector(mode = "numeric", length= 7)

array1.11 = vector(mode = "numeric", length= 7)
array2.11 = vector(mode = "numeric", length= 7)
array3.11 = vector(mode = "numeric", length= 7)



# GOING THROUGH THE RANGE OF 10-160 [PM10] TO FIND THE RR AT 10 UNIT INCREMENTS

for(i in 1:7){
  
  j = ((i*10)+10)
  j.1= as.character(j)
  
  # Cumulative:
  
  array1.1[i]= round(with(cpall1,cbind(allRRlow)[j.1,]),6)
  array2.1[i]= round(with(cpall1,cbind(allRRfit)[j.1,]),6)
  array3.1[i]= round(with(cpall1,cbind(allRRhigh)[j.1,]),6)
  
  array1.2[i]= round(with(cpall2,cbind(allRRlow)[j.1,]),6)
  array2.2[i]= round(with(cpall2,cbind(allRRfit)[j.1,]),6)
  array3.2[i]= round(with(cpall2,cbind(allRRhigh)[j.1,]),6)
  
  array1.3[i]= round(with(cpall3,cbind(allRRlow)[j.1,]),6)
  array2.3[i]= round(with(cpall3,cbind(allRRfit)[j.1,]),6)
  array3.3[i]= round(with(cpall3,cbind(allRRhigh)[j.1,]),6)
  
  # Lag specific:
  
  array1.4[i]= round(with(cpall_lag0,cbind(allRRlow)[j.1,]),6)
  array2.4[i]= round(with(cpall_lag0,cbind(allRRfit)[j.1,]),6)
  array3.4[i]= round(with(cpall_lag0,cbind(allRRhigh)[j.1,]),6)
  
  
  array1.5[i]= round(with(cpall_lag1,cbind(allRRlow)[j.1,]),6)
  array2.5[i]= round(with(cpall_lag1,cbind(allRRfit)[j.1,]),6)
  array3.5[i]= round(with(cpall_lag1,cbind(allRRhigh)[j.1,]),6)
  
  
  array1.6[i]= round(with(cpall_lag2,cbind(allRRlow)[j.1,]),6)
  array2.6[i]= round(with(cpall_lag2,cbind(allRRfit)[j.1,]),6)
  array3.6[i]= round(with(cpall_lag2,cbind(allRRhigh)[j.1,]),6)
  
  array1.7[i]= round(with(cpall_lag3,cbind(allRRlow)[j.1,]),6)
  array2.7[i]= round(with(cpall_lag3,cbind(allRRfit)[j.1,]),6)
  array3.7[i]= round(with(cpall_lag3,cbind(allRRhigh)[j.1,]),6)
  
  # Moving averages:
  
  array1.8[i]= round(with(cpall_mov1,cbind(allRRlow)[j.1,]),6)
  array2.8[i]= round(with(cpall_mov1,cbind(allRRfit)[j.1,]),6)
  array3.8[i]= round(with(cpall_mov1,cbind(allRRhigh)[j.1,]),6)
  
  array1.9[i]= round(with(cpall_mov2,cbind(allRRlow)[j.1,]),6)
  array2.9[i]= round(with(cpall_mov2,cbind(allRRfit)[j.1,]),6)
  array3.9[i]= round(with(cpall_mov2,cbind(allRRhigh)[j.1,]),6)
  
  array1.11[i]= round(with(cpall_mov3,cbind(allRRlow)[j.1,]),6)
  array2.11[i]= round(with(cpall_mov3,cbind(allRRfit)[j.1,]),6)
  array3.11[i]= round(with(cpall_mov3,cbind(allRRhigh)[j.1,]),6)
  
}

# Cumulative:

array_difference1.1 = vector(mode = "numeric", length= 6)
array_difference2.1 = vector(mode = "numeric", length= 6)
array_difference3.1 = vector(mode = "numeric", length= 6)


array_difference1.2 = vector(mode = "numeric", length= 6)
array_difference2.2 = vector(mode = "numeric", length= 6)
array_difference3.2 = vector(mode = "numeric", length= 6)


array_difference1.3 = vector(mode = "numeric", length= 6)
array_difference2.3 = vector(mode = "numeric", length= 6)
array_difference3.3 = vector(mode = "numeric", length= 6)

# Lag specific:

array_difference1.4 = vector(mode = "numeric", length= 6)
array_difference2.4 = vector(mode = "numeric", length= 6)
array_difference3.4 = vector(mode = "numeric", length= 6)


array_difference1.5 = vector(mode = "numeric", length= 6)
array_difference2.5 = vector(mode = "numeric", length= 6)
array_difference3.5 = vector(mode = "numeric", length= 6)


array_difference1.6 = vector(mode = "numeric", length= 6)
array_difference2.6 = vector(mode = "numeric", length= 6)
array_difference3.6 = vector(mode = "numeric", length= 6)


array_difference1.7 = vector(mode = "numeric", length= 6)
array_difference2.7 = vector(mode = "numeric", length= 6)
array_difference3.7 = vector(mode = "numeric", length= 6)

# Moving averages:

array_difference1.8 = vector(mode = "numeric", length= 6)
array_difference2.8 = vector(mode = "numeric", length= 6)
array_difference3.8 = vector(mode = "numeric", length= 6)


array_difference1.9 = vector(mode = "numeric", length= 6)
array_difference2.9 = vector(mode = "numeric", length= 6)
array_difference3.9 = vector(mode = "numeric", length= 6)

array_difference1.11 = vector(mode = "numeric", length= 6)
array_difference2.11 = vector(mode = "numeric", length= 6)
array_difference3.11 = vector(mode = "numeric", length= 6)


# STORING THE CHANGE OF RR PER 10 UNIT INCREASE OF PM10

for(i in 1:6){
  
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

print("Pooled Statistics")
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

####################################################################

