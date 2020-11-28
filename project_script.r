#STAT523_Project
#Md Muhtasim Billah


###getting started
#set working directory
setwd("/Users/muhtasim/Desktop/STAT523/Project")
#read in the dataset
mydata=read.csv("data.csv",header=T,colClasses=c("factor","factor","numeric"))
#manage plot window
par(mfrow=c(1,2))
#means by factor plot
plot.design(mydata, main = "Mean by Factor")
#interaction plot
interaction.plot(mydata$RecepLen,mydata$ReacCutoff,mydata$MaxRbonds,
                 xlab="Factor A : RecepLen", ylab="Average Max R_Bonds",
                 trace.label="Factor B : \n ReacCutoff", main= "Interaction Plot")
#two way additive ANOVA model
model=aov(MaxRbonds ~ RecepLen + ReacCutoff,data=mydata)
# ANOVA table for summary
summary(model)
#estimated parameter
coef(model) 
#manage plot window
par(mfrow=c(1,2))
#Tukey's simultaneous difference
TukeyHSD(model,conf.level=.95)
#Tukey's underscore method
source(url("http://math.wsu.edu/math/faculty/jpascual/stat423/R/two-way-T-method.R"))
twoway.t.method(model)
#Diagnostic plots for ANOVA
par(mfrow=c(1,2))
#residuals vs fitted plot
plot(model$fitted, model$res, main="Residuals vs Fitted Plot", 
     xlab = "Fitted Values", ylab = "Residuals")
#QQ plot
qqnorm(model$res)
qqline(model$res)

#Alternative options
#source(url("http://math.wsu.edu/math/faculty/jpascual/stat423/R/chap13fcns.R")) 
#residplots.lm(model,std.resid=F)

#Breusch-Pagan test
library(lmtest)
bptest(model)


##For simple linear regression
#data setup
regdata=read.csv("regdata.csv",header=T)
cutoff=regdata[,1] 
max_Rbonds=regdata[,2]
#scatter plot
plot(cutoff,max_Rbonds,xlab="Reaction Cutoff",ylab="Max R_Bonds")
#SLR
model1=lm(max_Rbonds~cutoff,data=regdata)
abline(model1)
#model parameters
s=summary(model1)
#coefficients(model1)
cf=confint(model1,level=0.95)
#prediction intervals
y=max_Rbonds
x=cutoff
predict(lm(y ~ x))
#creating new values for prediction
new <- data.frame(x = seq(0.0, 6.0, 0.25))
predict(lm(y ~ x), new, se.fit = TRUE)
pred.w.plim <- predict(lm(y ~ x), new, interval = "prediction")
pred.w.clim <- predict(lm(y ~ x), new, interval = "confidence")
#prediction plots
matplot(new$x, cbind(pred.w.clim, pred.w.plim[,-1]),  type = "l", 
        ylab = "Max R_Bonds", xlab = "Reaction Cutoff")
#Diagnostics plots for SLR
#Normal probability plot of Residuals
par(mfrow=c(1,2))
#plot of residual vs x and fitted values
plot(model1$fitted, model1$res, main = "Residuals vs Fitted Plot", 
     xlab = "Fitted Values", ylab = "Residuals")
qqnorm(model1$res)
qqline(model1$res)
#Diagnostic tests
#test for normality
library(stats)
shapiro.test(model1$res)
#breusch pagan test
#might need to install lmtest
library(lmtest)
bptest(model1)
#durbin-watson test
dwtest(model1)