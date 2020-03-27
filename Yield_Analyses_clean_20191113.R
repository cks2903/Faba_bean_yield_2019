
###################################################
#      Yield analysis NorFab 2019                 #
###################################################

setwd("/Volumes/ST_MBG-PMg/Cathrine/NorFab/Article/2. version Dec2019/")

library(lme4)
library(nlme)
library(ggplot2)
library(gridExtra)
library(RLRsim)
library(agricolae)
library(afex)
library(cowplot)
library(dmm)
library(car)
library(lmerTest)
library(nlme)
library(cowplot)
library(multcompView)
library(dplyr)



scaleFUN <- function(x) sprintf("%.2f", x)


# Load data
d=read.table("YieldDataJan2020.txt",sep="\t",header=T)
head(d)
d$Yield..g..Pr.kvm.=d$Yield..kg..Pr.kvm.*1000
head(d)
# Variance by cultivar name
completeFun <- function(data, desiredCols) {
  completeVec <- complete.cases(data[, desiredCols])
  return(data[completeVec, ])
}
newd=completeFun(d,"Yield..g..Pr.kvm.")

variances=aggregate(newd$Yield..g..Pr.kvm.,list(newd$Cultivar),var)
mean=aggregate(newd$Yield..g..Pr.kvm.,list(newd$Cultivar),mean)
sd_yieldmeans=aggregate(newd$Yield..g..Pr.kvm.,list(newd$Cultivar),sd)
sd_yieldmeans
repcultivar=count(newd, vars = Cultivar)
se_yieldMeans=sd_yieldmeans[,2]/sqrt(repcultivar$n) #Se of cultivar means

total=cbind(variances,mean)
colnames(total)=c("Cultivar","YieldVar","Cultivar2","MeanYield")
total$MeanYield

# Produce plot yield vs. variance
plot=ggplot(total, aes(x=MeanYield, y=YieldVar))+geom_point(col="#00ba38",cex=3.5,alpha=0.8) +
  geom_text(label=total$Cultivar,size=4.5, vjust=0.5) +
  #scale_y_continuous(labels=scaleFUN) +
  #scale_x_continuous(labels=scaleFUN) +
  ylab("Variance of log(yield)") +
  xlab("Mean yield g/sqm") +
  ylim(c(min(total$YieldVar),max(total$YieldVar) )) +
  xlim(c(min(total$MeanYield),max(total$MeanYield) )) +
  geom_vline(xintercept=mean(total$MeanYield),colour="black",linetype=2) +
  geom_hline(yintercept=mean(total$YieldVar),colour="black",linetype=2)

plot+theme_bw()

# Remove plants with no observations
rowstoremove=which(is.na(d$Yield..g..Pr.kvm.))
filtered=d[-rowstoremove,]


#### Analyses

# Make a Year-location interaction column
filtered$YearLoc=interaction(filtered$Year,filtered$Location)

#Calculating mean environmental yield
filtered$Environmental_yield_mean=filtered$Yield..g..Pr.kvm.

NS2016=subset(filtered,filtered$YearLoc=="1.N")
mean1=mean(NS2016$Yield..g..Pr.kvm.)
NS2017=subset(filtered,filtered$YearLoc=="2.N")
mean2=mean(na.omit(NS2017$Yield..g..Pr.kvm.))
NS2018=subset(filtered,filtered$YearLoc=="3.N")
mean3=mean(na.omit(NS2018$Yield..g..Pr.kvm.))
S2017=subset(filtered,filtered$YearLoc=="2.S")
mean4=mean(S2017$Yield..g..Pr.kvm.)
S2018=subset(filtered,filtered$YearLoc=="3.S")
mean5=mean(S2018$Yield..g..Pr.kvm.)
SF2017=subset(filtered,filtered$YearLoc=="2.SF")
mean6=mean(SF2017$Yield..g..Pr.kvm.)
SF2018=subset(filtered,filtered$YearLoc=="3.SF")
mean7=mean(SF2018$Yield..g..Pr.kvm.)
V2016=subset(filtered,filtered$YearLoc=="1.V")
mean8=mean(V2016$Yield..g..Pr.kvm.)
V2018=subset(filtered,filtered$YearLoc=="3.V")
mean9=mean(V2018$Yield..g..Pr.kvm.)

means_pr_trial=c(mean5,mean7,mean9,mean3,mean1,mean6,mean8,mean4,mean2)
means_pr_trial=as.data.frame(means_pr_trial)
rownames(means_pr_trial)=c("Sejet_J-2018","Sejet_F-2018","Finland-2018","Nordic_Seed-2018","Nordic_Seed-2016","Sejet_F-2017","Finland-2016","Sejet_J-2017","Nordic_Seed-2017")

filtered$Environmental_yield_mean[which(filtered$YearLoc=="3.S")]=means_pr_trial[1,1]
filtered$Environmental_yield_mean[which(filtered$YearLoc=="3.SF")]=means_pr_trial[2,1]
filtered$Environmental_yield_mean[which(filtered$YearLoc=="3.V")]=means_pr_trial[3,1]
filtered$Environmental_yield_mean[which(filtered$YearLoc=="3.N")]=means_pr_trial[4,1]
filtered$Environmental_yield_mean[which(filtered$YearLoc=="1.N")]=means_pr_trial[5,1]
filtered$Environmental_yield_mean[which(filtered$YearLoc=="2.SF")]=means_pr_trial[6,1]
filtered$Environmental_yield_mean[which(filtered$YearLoc=="1.V")]=means_pr_trial[7,1]
filtered$Environmental_yield_mean[which(filtered$YearLoc=="2.S")]=means_pr_trial[8,1]
filtered$Environmental_yield_mean[which(filtered$YearLoc=="2.N")]=means_pr_trial[9,1]

length(unique(filtered$Environmental_yield_mean))
filtered$Scaled_Environmental_yield_mean=scale(filtered$Environmental_yield_mean) #scale to get variance=1, mean=0


#Fitting LMM on seed yield

M6<- lmer(Yield..g..Pr.kvm. ~ (1|Location)+(1|Year)+(1|Cultivar:Location)+(1|Cultivar:Year)+(1|YearLoc)+ (Scaled_Environmental_yield_mean|Cultivar), data = filtered,REML=T) 
summary(M6)
Residuals <- residuals(M6)
shapiro.test(Residuals) 
qqnorm(Residuals); qqline(Residuals) ##A few genotypes perform poorly compared to environmental mean. But data is not wrong
plot(M6, abline=c(0, 0))
lattice::qqmath(M6,id=0.1,idLabels=~.obs) 

######## Testing significance of all random effects
# do a likelihood ratio test
M6_<- lmer(Yield..g..Pr.kvm. ~ (1|Location)+(1|Year)+(1|Cultivar:Location)+(1|Cultivar:Year)+(1|YearLoc)+ (Scaled_Environmental_yield_mean|Cultivar), data = filtered,REML=F) 
summary(M6_)

M6_<- lmer(Yield..g..Pr.kvm. ~ (1|Location)+(1|Year)+(1|YearLoc)+ (1|Cultivar), data = filtered,REML=F) 
M6_0<- lmer(Yield..g..Pr.kvm. ~ (1|Location)+(1|Year)+(1|YearLoc), data = filtered,REML=F) 

summary(M6_)
summary(M6_0) 
anova(M6_,M6_0) #Cultivar has a significant effect ***

M6__<- lmer(Yield..g..Pr.kvm. ~ (1|Cultivar) + (1|Location)+(1|Year)+(1|Cultivar:Year), data = filtered,REML=F) 
M6__0<- lmer(Yield..g..Pr.kvm. ~ (1|Cultivar) + (1|Year)+(1|Cultivar:Year), data = filtered,REML=F) 
anova(M6__,M6__0) #Location has a significant effect ***

M6_<- lmer(Yield..g..Pr.kvm. ~ (1|Cultivar) + (1|Location) + (1|Year) + (1|Cultivar:Location), data = filtered,REML=F) 
summary(M6_)
M6_0<- lmer(Yield..g..Pr.kvm. ~ (1|Cultivar) + (1|Location) + (1|Cultivar:Location), data = filtered,REML=F) 
summary(M6_0)
anova(M6_,M6_0) #Year has a significant effect ***

M6_<- lmer(Yield..g..Pr.kvm. ~ (1|Location)+(1|Year)+(1|Cultivar:Location)+(1|Cultivar:Year)+(1|YearLoc)+ (Scaled_Environmental_yield_mean|Cultivar), data = filtered,REML=F) 
summary(M6)
M6_0<- lmer(Yield..g..Pr.kvm. ~ (1|Location)+(1|Year)+(1|Cultivar:Year)+(1|YearLoc)+ (Scaled_Environmental_yield_mean|Cultivar), data = filtered,REML=F) 
summary(M6_0)
anova(M6_,M6_0) #Cultivar:location has significant effect ***

M6_<- lmer(Yield..g..Pr.kvm. ~ (1|Location)+(1|Year)+(1|Cultivar:Location)+(1|Cultivar:Year)+(1|YearLoc)+ (Scaled_Environmental_yield_mean|Cultivar), data = filtered,REML=F) 
summary(M6)
M6_0<- lmer(Yield..g..Pr.kvm. ~ (1|Location)+(1|Year)+(1|Cultivar:Location)+(1|YearLoc)+ (Scaled_Environmental_yield_mean|Cultivar), data = filtered,REML=F) 
summary(M6_0)
anova(M6_,M6_0) #Cultivar:Year has significant effect *

M6_<- lmer(Yield..g..Pr.kvm. ~ (1|Location)+(1|Year)+(1|Cultivar:Location)+(1|Cultivar:Year)+(1|YearLoc)+ (Scaled_Environmental_yield_mean|Cultivar), data = filtered,REML=F) 
summary(M6)
M6_0<- lmer(Yield..g..Pr.kvm. ~ (1|Location)+(1|Year)+(1|Cultivar:Location)+ (Scaled_Environmental_yield_mean|Cultivar), data = filtered,REML=F) 
summary(M6_0)
anova(M6_,M6_0) #Location:Year has significant effect **

M6_<- lmer(Yield..g..Pr.kvm. ~ (1|Location)+(1|Year)+(1|Cultivar:Location)+(1|Cultivar:Year)+(1|YearLoc)+ (Scaled_Environmental_yield_mean|Cultivar), data = filtered,REML=F) 
summary(M6)
M6_0<- lmer(Yield..g..Pr.kvm. ~ (1|Location)+(1|Year)+(1|Cultivar:Location)+(1|Cultivar:Year)+(1|YearLoc)+ (1|Cultivar), data = filtered,REML=F) 
summary(M6_0)
anova(M6_,M6_0) #Cultivar by FW regression has significant effect  ***

# Calculating BLUPs for all random effect
Cultivar_effects=ranef(M6)$Cultivar[1] 

Highyield=subset(Cultivar_effects,Cultivar_effects>0) #Subset cultivars with positive yield effect

Variances=tapply(filtered$Log_Yield, filtered$Cultivar, var)
Variances=data.matrix(Variances)
Variances=Variances[-6,1]
Variances=data.matrix(Variances)
Variances=Variances[-13,1]
Variances=data.matrix(Variances)


############ Calculate CV%, coefficient of variation 

CV_mean=tapply(filtered$Yield..g..Pr.kvm., filtered$Cultivar, mean)

CV_stdev=tapply(filtered$Yield..g..Pr.kvm., filtered$Cultivar, sd)

CV_var=tapply(filtered$Yield..g..Pr.kvm., filtered$Cultivar, var)

CVpercent=CV_stdev/CV_mean
CVpercent




### Finlay-Wilkinson regression, produce Figure 2

# Calculate mean of every cultivar for each trial
NS2016_CultivarMeans=aggregate(NS2016$Yield..g..Pr.kvm.,list(NS2016$Cultivar),mean)
NS2016_CultivarMeans=as.data.frame(NS2016_CultivarMeans)
NS2017_CultivarMeans=aggregate(NS2017$Yield..g..Pr.kvm.,list(NS2017$Cultivar),mean)
NS2017_CultivarMeans=as.data.frame(NS2017_CultivarMeans)
NS2018_CultivarMeans=aggregate(NS2018$Yield..g..Pr.kvm.,list(NS2018$Cultivar),mean)
NS2018_CultivarMeans=as.data.frame(NS2018_CultivarMeans)
S2017_CultivarMeans=aggregate(S2017$Yield..g..Pr.kvm.,list(S2017$Cultivar),mean)
S2017_CultivarMeans=as.data.frame(S2017_CultivarMeans)
S2018_CultivarMeans=aggregate(S2018$Yield..g..Pr.kvm.,list(S2018$Cultivar),mean)
S2018_CultivarMeans=as.data.frame(S2018_CultivarMeans)
SF2017_CultivarMeans=aggregate(SF2017$Yield..g..Pr.kvm.,list(SF2017$Cultivar),mean)
SF2017_CultivarMeans=as.data.frame(SF2017_CultivarMeans)
SF2018_CultivarMeans=aggregate(SF2018$Yield..g..Pr.kvm.,list(SF2018$Cultivar),mean)
SF2018_CultivarMeans=as.data.frame(SF2018_CultivarMeans)
V2016_CultivarMeans=aggregate(V2016$Yield..g..Pr.kvm.,list(V2016$Cultivar),mean)
V2016_CultivarMeans=as.data.frame(V2016_CultivarMeans)
V2018_CultivarMeans=aggregate(V2018$Yield..g..Pr.kvm.,list(V2018$Cultivar),mean)
V2018_CultivarMeans=as.data.frame(V2018_CultivarMeans)


All_years=Reduce(function(x, y) merge(x, y, by="Group.1", all=TRUE), list(S2018_CultivarMeans, SF2018_CultivarMeans, V2018_CultivarMeans,NS2018_CultivarMeans,NS2016_CultivarMeans,SF2017_CultivarMeans,V2016_CultivarMeans,S2017_CultivarMeans,NS2017_CultivarMeans))
colnames(All_years)=c("Cultivar","S2018","SF2018","V2018","NS2018","NS2016","SF2017","V2016","S2017","NS2017")
transpose=t(All_years)

means_pr_trial

merged=cbind(means_pr_trial,transpose[2:nrow(transpose),])
colnames(merged)=c("means_pr_trial",transpose[1,])
head(merged)

pdf(file="YieldStabilityAllPlots_1.pdf",width=65,height=75,useDingbats = F)

lm_fit=lm(unfactor(`247-13`)~means_pr_trial,data=merged)
summary(lm_fit)
coef(lm_fit)
#Wald Test to test if slope is significantly different from zero
linearHypothesis(lm_fit,hypothesis.matrix=c(0,1), rhs=1.0) 
p1=ggplot(merged,aes(means_pr_trial,unfactor(`247-13`)))+geom_point(col="#00ba38",cex=10,alpha=0.8) + ylim(c(100, 1000)) + xlim(c(100,1000)) + ggtitle("247-13") +geom_abline(intercept = 0, slope = 1,col="black",cex=2) + geom_abline(col="#00ba38",intercept=coef(lm_fit)[1],slope=coef(lm_fit)[2],cex=4) +
  theme(axis.text.x=element_text(size=55),
        axis.text.y=element_text(size=55),
        plot.title = element_text(size=75),) +
  geom_text(size=20,aes(x=750,y=150,hjust=0,vjust=0,label="0.96"))+ theme(
    panel.background = element_rect(fill="white",
                                    colour = "white",
                                    size = 0.5, linetype = "solid"),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                    colour = "gray93"), 
    panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                    colour = "gray93")
    
  ) 

                             
lm_fit=lm(unfactor(`749-13`)~means_pr_trial,data=merged)                                  
summary(lm_fit)
linearHypothesis(lm_fit,hypothesis.matrix=c(0,1), rhs=1.0) 
p2=ggplot(merged,aes(means_pr_trial,unfactor(`749-13`)))+geom_point(col="#00ba38",cex=10,alpha=0.8) + ylim(c(100, 1000)) + xlim(c(100,1000)) + ggtitle("749-13") +geom_abline(intercept = 0, slope = 1,col="black",cex=2) + geom_abline(col="#00ba38",intercept=coef(lm_fit)[1],slope=coef(lm_fit)[2],cex=4) +
  theme(axis.text.x=element_text(size=55),
        axis.text.y=element_text(size=55),
        plot.title = element_text(size=75),) +
  geom_text(size=20,aes(x=750,y=150,hjust=0,vjust=0,label="0.91"))+ theme(
    panel.background = element_rect(fill="white",
                                    colour = "white",
                                    size = 0.5, linetype = "solid"),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                    colour = "gray93"), 
    panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                    colour = "gray93")
    
  )
                                     
lm_fit=lm(unfactor(`Alexia`)~means_pr_trial,data=merged)                                  
summary(lm_fit)
linearHypothesis(lm_fit,hypothesis.matrix=c(0,1), rhs=1.0) 
p3=ggplot(merged,aes(means_pr_trial,unfactor(`Alexia`)))+geom_point(col="#00ba38",cex=10,alpha=0.8) + ylim(c(100, 1000)) + xlim(c(100,1000)) + ggtitle("Alexia") +geom_abline(intercept = 0, slope = 1,col="black",cex=2) + geom_abline(col="#00ba38",intercept=coef(lm_fit)[1],slope=coef(lm_fit)[2],cex=4) +
  theme(axis.text.x=element_text(size=55),
        axis.text.y=element_text(size=55),
        plot.title = element_text(size=75),) +
  geom_text(size=20,aes(x=750,y=150,hjust=0,vjust=0,label="0.97"))+ theme(
    panel.background = element_rect(fill="white",
                                    colour = "white",
                                    size = 0.5, linetype = "solid"),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                    colour = "gray93"), 
    panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                    colour = "gray93")
    
  )                                    
                                                                                                                               
lm_fit=lm(unfactor(`Banquise`)~means_pr_trial,data=merged)                                  
summary(lm_fit)
linearHypothesis(lm_fit,hypothesis.matrix=c(0,1), rhs=1.0) 
p4=ggplot(merged,aes(means_pr_trial,unfactor(`Banquise`)))+geom_point(col="#00ba38",cex=10,alpha=0.8) + ylim(c(100, 1000)) + xlim(c(100,1000)) + ggtitle("Banquise") +geom_abline(intercept = 0, slope = 1,col="black",cex=2) + geom_abline(col="#00ba38",intercept=coef(lm_fit)[1],slope=coef(lm_fit)[2],cex=4) +
  theme(axis.text.x=element_text(size=55),
        axis.text.y=element_text(size=55),
        plot.title = element_text(size=75),) +
  geom_text(size=20,aes(x=750,y=150,hjust=0,vjust=0,label="1.04"))+ theme(
    panel.background = element_rect(fill="white",
                                    colour = "white",
                                    size = 0.5, linetype = "solid"),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                    colour = "gray93"), 
    panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                    colour = "gray93")
    
  )                                                                                                                                 

lm_fit=lm(unfactor(`Boxer`)~means_pr_trial,data=merged)                                  
summary(lm_fit)
linearHypothesis(lm_fit,hypothesis.matrix=c(0,1), rhs=1.0) 
p5=ggplot(merged,aes(means_pr_trial,unfactor(`Boxer`)))+geom_point(col="#00ba38",cex=10,alpha=0.8) + ylim(c(100, 1000)) + xlim(c(100,1000)) + ggtitle("Boxer") +geom_abline(intercept = 0, slope = 1,col="black",cex=2) + geom_abline(col="#00ba38",intercept=coef(lm_fit)[1],slope=coef(lm_fit)[2],cex=4) +
  theme(axis.text.x=element_text(size=55),
        axis.text.y=element_text(size=55),
        plot.title = element_text(size=75),) +
  geom_text(size=20,aes(x=750,y=150,hjust=0,vjust=0,label="1.12"))+ theme(
    panel.background = element_rect(fill="white",
                                    colour = "white",
                                    size = 0.5, linetype = "solid"),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                    colour = "gray93"), 
    panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                    colour = "gray93")
    
  )                                                                                                                                
                                                                                                                            
lm_fit=lm(unfactor(`Fanfare`)~means_pr_trial,data=merged)                                  
summary(lm_fit)
linearHypothesis(lm_fit,hypothesis.matrix=c(0,1), rhs=1.0) 
p6=ggplot(merged,aes(means_pr_trial,unfactor(`Fanfare`)))+geom_point(col="#00ba38",cex=10,alpha=0.8) + ylim(c(100, 1000)) + xlim(c(100,1000)) + ggtitle("Fanfare") +geom_abline(intercept = 0, slope = 1,col="black",cex=2) + geom_abline(col="#00ba38",intercept=coef(lm_fit)[1],slope=coef(lm_fit)[2],cex=4) +
  theme(axis.text.x=element_text(size=55),
        axis.text.y=element_text(size=55),
        plot.title = element_text(size=75),) +
  geom_text(size=20,aes(x=750,y=150,hjust=0,vjust=0,label="1.14"))+ theme(
    panel.background = element_rect(fill="white",
                                    colour = "white",
                                    size = 0.5, linetype = "solid"),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                    colour = "gray93"), 
    panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                    colour = "gray93")
    
  )  


lm_fit=lm(unfactor(`Fuego`)~means_pr_trial,data=merged)                                  
summary(lm_fit)
linearHypothesis(lm_fit,hypothesis.matrix=c(0,1), rhs=1.0) 
p7=ggplot(merged,aes(means_pr_trial,unfactor(`Fuego`)))+geom_point(col="#00ba38",cex=10,alpha=0.8) + ylim(c(100, 1000)) + xlim(c(100,1000)) + ggtitle("Fuego") +geom_abline(intercept = 0, slope = 1,col="black",cex=2) + geom_abline(col="#00ba38",intercept=coef(lm_fit)[1],slope=coef(lm_fit)[2],cex=4) +
  theme(axis.text.x=element_text(size=55),
        axis.text.y=element_text(size=55),
        plot.title = element_text(size=75),) +
  geom_text(size=20,aes(x=750,y=150,hjust=0,vjust=0,label="1.06"))+ theme(
    panel.background = element_rect(fill="white",
                                    colour = "white",
                                    size = 0.5, linetype = "solid"),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                    colour = "gray93"), 
    panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                    colour = "gray93")
    
  )  

lm_fit=lm(unfactor(`Gloria`)~means_pr_trial,data=merged)                                  
summary(lm_fit)
linearHypothesis(lm_fit,hypothesis.matrix=c(0,1), rhs=1.0) 
p8=ggplot(merged,aes(means_pr_trial,unfactor(`Gloria`)))+geom_point(col="#00ba38",cex=10,alpha=0.8) + ylim(c(100, 1000)) + xlim(c(100,1000)) + ggtitle("Gloria") +geom_abline(intercept = 0, slope = 1,col="black",cex=2) + geom_abline(col="#00ba38",intercept=coef(lm_fit)[1],slope=coef(lm_fit)[2],cex=4) +
  theme(axis.text.x=element_text(size=55),
        axis.text.y=element_text(size=55),
        plot.title = element_text(size=75),) +
  geom_text(size=20,aes(x=750,y=150,hjust=0,vjust=0,label="1.02"))+ theme(
    panel.background = element_rect(fill="white",
                                    colour = "white",
                                    size = 0.5, linetype = "solid"),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                    colour = "gray93"), 
    panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                    colour = "gray93")
    
  )  

lm_fit=lm(unfactor(`Gracia`)~means_pr_trial,data=merged)                                  
summary(lm_fit)
linearHypothesis(lm_fit,hypothesis.matrix=c(0,1), rhs=1.0)  
p9=ggplot(merged,aes(means_pr_trial,unfactor(`Gracia`)))+geom_point(col="#00ba38",cex=10,alpha=0.8) + ylim(c(100, 1000)) + xlim(c(100,1000)) + ggtitle("Gracia") +geom_abline(intercept = 0, slope = 1,col="black",cex=2) + geom_abline(col="#00ba38",intercept=coef(lm_fit)[1],slope=coef(lm_fit)[2],cex=4) +
  theme(axis.text.x=element_text(size=55),
        axis.text.y=element_text(size=55),
        plot.title = element_text(size=75),) +
  geom_text(size=20,aes(x=750,y=150,hjust=0,vjust=0,label="1.12 *"))+ theme(
    panel.background = element_rect(fill="white",
                                    colour = "white",
                                    size = 0.5, linetype = "solid"),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                    colour = "gray93"), 
    panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                    colour = "gray93")
    
  )  

lm_fit=lm(unfactor(`Kontu`)~means_pr_trial,data=merged)                                  
summary(lm_fit)
linearHypothesis(lm_fit,hypothesis.matrix=c(0,1), rhs=1.0)  
p10=ggplot(merged,aes(means_pr_trial,unfactor(`Kontu`)))+geom_point(col="#00ba38",cex=10,alpha=0.8) + ylim(c(100, 1000)) + xlim(c(100,1000)) + ggtitle("Kontu") +geom_abline(intercept = 0, slope = 1,col="black",cex=2) + geom_abline(col="#00ba38",intercept=coef(lm_fit)[1],slope=coef(lm_fit)[2],cex=4) +
  theme(axis.text.x=element_text(size=55),
        axis.text.y=element_text(size=55),
        plot.title = element_text(size=75),) +
  geom_text(size=20,aes(x=750,y=150,hjust=0,vjust=0,label="0.72 *"))+ theme(
    panel.background = element_rect(fill="white",
                                    colour = "white",
                                    size = 0.5, linetype = "solid"),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                    colour = "gray93"), 
    panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                    colour = "gray93")
    
  )  

lm_fit=lm(unfactor(`Lynx`)~means_pr_trial,data=merged)                                  
summary(lm_fit)
linearHypothesis(lm_fit,hypothesis.matrix=c(0,1), rhs=1.0)  
p11=ggplot(merged,aes(means_pr_trial,unfactor(`Lynx`)))+geom_point(col="#00ba38",cex=10,alpha=0.8) + ylim(c(100, 1000)) + xlim(c(100,1000)) + ggtitle("Lynx") +geom_abline(intercept = 0, slope = 1,col="black",cex=2) + geom_abline(col="#00ba38",intercept=coef(lm_fit)[1],slope=coef(lm_fit)[2],cex=4) +
  theme(axis.text.x=element_text(size=55),
        axis.text.y=element_text(size=55),
        plot.title = element_text(size=75),) +
  geom_text(size=20,aes(x=750,y=150,hjust=0,vjust=0,label="1.13"))+ theme(
    panel.background = element_rect(fill="white",
                                    colour = "white",
                                    size = 0.5, linetype = "solid"),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                    colour = "gray93"), 
    panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                    colour = "gray93")
    
  )

lm_fit=lm(unfactor(`Mistral`)~means_pr_trial,data=merged)                                  
summary(lm_fit)
linearHypothesis(lm_fit,hypothesis.matrix=c(0,1), rhs=1.0)  
p12=ggplot(merged,aes(means_pr_trial,unfactor(`Mistral`)))+geom_point(col="#00ba38",cex=10,alpha=0.8) + ylim(c(100, 1000)) + xlim(c(100,1000)) + ggtitle("Mistral") +geom_abline(intercept = 0, slope = 1,col="black",cex=2) + geom_abline(col="#00ba38",intercept=coef(lm_fit)[1],slope=coef(lm_fit)[2],cex=4) +
  theme(axis.text.x=element_text(size=55),
        axis.text.y=element_text(size=55),
        plot.title = element_text(size=75),) +
  geom_text(size=20,aes(x=750,y=150,hjust=0,vjust=0,label="0.92"))+ theme(
    panel.background = element_rect(fill="white",
                                    colour = "white",
                                    size = 0.5, linetype = "solid"),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                    colour = "gray93"), 
    panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                    colour = "gray93")
    
  )


lm_fit=lm(unfactor(`Pyramid`)~means_pr_trial,data=merged)                                  
summary(lm_fit)
linearHypothesis(lm_fit,hypothesis.matrix=c(0,1), rhs=1.0)  
p13=ggplot(merged,aes(means_pr_trial,unfactor(`Pyramid`)))+geom_point(col="#00ba38",cex=10,alpha=0.8) + ylim(c(100, 1000)) + xlim(c(100,1000)) + ggtitle("Pyramid") +geom_abline(intercept = 0, slope = 1,col="black",cex=2) + geom_abline(col="#00ba38",intercept=coef(lm_fit)[1],slope=coef(lm_fit)[2],cex=4) +
  theme(axis.text.x=element_text(size=55),
        axis.text.y=element_text(size=55),
        plot.title = element_text(size=75),) +
  geom_text(size=20,aes(x=750,y=150,hjust=0,vjust=0,label="1.11 *"))+ theme(
    panel.background = element_rect(fill="white",
                                    colour = "white",
                                    size = 0.5, linetype = "solid"),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                    colour = "gray93"), 
    panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                    colour = "gray93")
    
  )  

lm_fit=lm(unfactor(`Snowdrop`)~means_pr_trial,data=merged)                                  
summary(lm_fit)
linearHypothesis(lm_fit,hypothesis.matrix=c(0,1), rhs=1.0)  
p14=ggplot(merged,aes(means_pr_trial,unfactor(`Snowdrop`)))+geom_point(col="#00ba38",cex=10,alpha=0.8) + ylim(c(100, 1000)) + xlim(c(100,1000)) + ggtitle("Snowdrop") +geom_abline(intercept = 0, slope = 1,col="black",cex=2) + geom_abline(col="#00ba38",intercept=coef(lm_fit)[1],slope=coef(lm_fit)[2],cex=4) +
  theme(axis.text.x=element_text(size=55),
        axis.text.y=element_text(size=55),
        plot.title = element_text(size=75),) +
  geom_text(size=20,aes(x=750,y=150,hjust=0,vjust=0,label="0.62 **"))+ theme(
    panel.background = element_rect(fill="white",
                                    colour = "white",
                                    size = 0.5, linetype = "solid"),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                    colour = "gray93"), 
    panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                    colour = "gray93")
    
  )  

lm_fit=lm(unfactor(`SSNS-1`)~means_pr_trial,data=merged)                                  
summary(lm_fit)
linearHypothesis(lm_fit,hypothesis.matrix=c(0,1), rhs=1.0)  
p15=ggplot(merged,aes(means_pr_trial,unfactor(`SSNS-1`)))+geom_point(col="#00ba38",cex=10,alpha=0.8) + ylim(c(100, 1000)) + xlim(c(100,1000)) + ggtitle("SSNS-1") +geom_abline(intercept = 0, slope = 1,col="black",cex=2) + geom_abline(col="#00ba38",intercept=coef(lm_fit)[1],slope=coef(lm_fit)[2],cex=4) +
  theme(axis.text.x=element_text(size=55),
        axis.text.y=element_text(size=55),
        plot.title = element_text(size=75),) +
  geom_text(size=20,aes(x=750,y=150,hjust=0,vjust=0,label="0.90 *"))+ theme(
    panel.background = element_rect(fill="white",
                                    colour = "white",
                                    size = 0.5, linetype = "solid"),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                    colour = "gray93"), 
    panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                    colour = "gray93")
    
  )  

lm_fit=lm(unfactor(`Taifun`)~means_pr_trial,data=merged)                                  
summary(lm_fit)
linearHypothesis(lm_fit,hypothesis.matrix=c(0,1), rhs=1.0)  
p16=ggplot(merged,aes(means_pr_trial,unfactor(`Taifun`)))+geom_point(col="#00ba38",cex=10,alpha=0.8) + ylim(c(100, 1000)) + xlim(c(100,1000)) + ggtitle("Taifun") +geom_abline(intercept = 0, slope = 1,col="black",cex=2) + geom_abline(col="#00ba38",intercept=coef(lm_fit)[1],slope=coef(lm_fit)[2],cex=4) +
  theme(axis.text.x=element_text(size=55),
        axis.text.y=element_text(size=55),
        plot.title = element_text(size=75),) +
  geom_text(size=20,aes(x=750,y=150,hjust=0,vjust=0,label="1.01"))+ theme(
    panel.background = element_rect(fill="white",
                                    colour = "white",
                                    size = 0.5, linetype = "solid"),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                    colour = "gray93"), 
    panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                    colour = "gray93")
                                  
  )
 

lm_fit=lm(unfactor(`Vertigo`)~means_pr_trial,data=merged)                                  
summary(lm_fit)
linearHypothesis(lm_fit,hypothesis.matrix=c(0,1), rhs=1.0)  
p17=ggplot(merged,aes(means_pr_trial,unfactor(`Vertigo`)))+geom_point(col="#00ba38",cex=10,alpha=0.8) + ylim(c(100, 1000)) + xlim(c(100,1000)) + ggtitle("Vertigo") +geom_abline(intercept = 0, slope = 1,col="black",cex=2) + geom_abline(col="#00ba38",intercept=coef(lm_fit)[1],slope=coef(lm_fit)[2],cex=4) +
  theme(axis.text.x=element_text(size=55),
        axis.text.y=element_text(size=55),
        plot.title = element_text(size=75),) +
  geom_text(size=20,aes(x=750,y=150,hjust=0,vjust=0,label="1.19 *"))+ theme(
    panel.background = element_rect(fill = "white",
                                    colour = "white",
                                    size = 0.5, linetype = "solid"),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                    colour = "gray93"), 
    panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                    colour = "gray93")
  )

plot_grid(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17, nrow=5,align="v")

dev.off()




####################################################################
#                                                                  #
#         Correlation of yield and protein content of seeds        #
#                                                                  #
####################################################################
options(scipen = 999)


proteindata=which(filtered$Protein!="NA")
Filtered_protein=filtered[proteindata,]
nrow(Filtered_protein)


# Get environmental mean for protein

#Yield*protein means
mean1_protein=mean(NS2016$Protein)
mean2_protein=mean(NS2017$Protein)
mean3_protein=mean(NS2018$Protein)
mean4_protein=mean(S2017$Protein)
mean5_protein=mean(V2016$Protein)

means_pr_trial_protein=c(mean1_protein,mean2_protein,mean3_protein,mean4_protein,mean5_protein)
means_pr_trial_protein=as.data.frame(means_pr_trial_protein)
rownames(means_pr_trial_protein)=c("Nordic_Seed-2016","Nordic_Seed-2017","Nordic_Seed-2018","Sejet_J-2017","FN_2016")

Filtered_protein$Environmental_protein_mean=Filtered_protein$Protein

Filtered_protein$Environmental_protein_mean[which(Filtered_protein$YearLoc=="1.N")]=means_pr_trial_protein[1,1]
Filtered_protein$Environmental_protein_mean[which(Filtered_protein$YearLoc=="2.N")]=means_pr_trial_protein[2,1]
Filtered_protein$Environmental_protein_mean[which(Filtered_protein$YearLoc=="3.N")]=means_pr_trial_protein[3,1]
Filtered_protein$Environmental_protein_mean[which(Filtered_protein$YearLoc=="2.S")]=means_pr_trial_protein[4,1]
Filtered_protein$Environmental_protein_mean[which(Filtered_protein$YearLoc=="1.V")]=means_pr_trial_protein[5,1]

length(unique(Filtered_protein$Environmental_protein_mean))
Filtered_protein$Scaled_Environmental_protein_mean=scale(Filtered_protein$Environmental_protein_mean) #scale to get variance=1, mean=0





# Calculate phenotypic correlation between yield and protein content
correlation(Filtered_protein$Protein,Filtered_protein$Yield..g..Pr.kvm.) #-0.08, not significant


# Fit model to protein data
CM<- lmer(Protein ~ (1|Location)+(1|Year)+(1|Cultivar:Location)+(1|Cultivar:Year)+ (Scaled_Environmental_protein_mean|Cultivar), data = Filtered_protein,REML=T) 
summary(CM)

lol<- lmer(Protein ~ (1|Location)+(1|Year)+(1|YearLoc), data = Filtered_protein,REML=T) 
summary(lol)

Residuals <- residuals(CM)
shapiro.test(Residuals) #fits
qqnorm(Residuals); qqline(Residuals)
plot(CM, abline=c(0, 0))


#Genetic correlations between protein and yield
Cultivar_effects_protein=ranef(CM)$Cultivar[1] 

correlation(Cultivar_effects_protein,Cultivar_effects) #BLUPs of yield on all data, -0.61 **

plot(Cultivar_effects_protein[,1]~Cultivar_effects_yield[,1])
new=cbind(rownames(Cultivar_effects_protein),as.numeric(as.character(Cultivar_effects_protein[,1])),as.numeric(as.character(Cultivar_effects[,1])))
colnames(new)=c("cultivar","BLUPS_Protein","BLUPS_Yield")
new=as.data.frame(new)


# Test significance and size of different parameters in mixed model of protein

######## Testing significance of all random effects
# do a likelihood ratio test
CM_ML<- lmer(Protein ~ (1|Location)+(1|Year)+(1|Cultivar:Location)+(1|Cultivar:Year)+ (Scaled_Environmental_protein_mean|Cultivar), data = Filtered_protein,REML=F) 
summary(CM_ML)

#Test significans of interaction effects
CM_0<- lmer(Protein ~ (1|Location)+(1|Year)+(1|Cultivar:Year)+ (Scaled_Environmental_protein_mean|Cultivar), data = Filtered_protein,REML=F) 
anova(CM_ML,CM_0) #Cultivar:Location has a significant effect ***

CM_0<- lmer(Protein ~ (1|Location)+(1|Year)+(1|Cultivar:Location)+ (Scaled_Environmental_protein_mean|Cultivar), data = Filtered_protein,REML=F) 
anova(CM_ML,CM_0) #Cultivar:Year has a significant effect ***

CM_0<- lmer(Protein ~ (1|Location)+(1|Year)+(1|Cultivar:Location)+(1|Cultivar:Year)+ (1|Cultivar), data = Filtered_protein,REML=F) 
anova(CM_ML,CM_0) #Cultivar by FW regression do not have a significant effect

CM_<- lmer(Protein ~ (1|Cultivar) + (1|Location)+(1|Year)+(1|Cultivar:Year), data = Filtered_protein,REML=F) 
CM_0<- lmer(Protein ~ (1|Cultivar) + (1|Year)+(1|Cultivar:Year), data = Filtered_protein,REML=F) 
anova(CM_,CM_0)#Location has a significant effect, ***

CM_<- lmer(Protein ~ (1|Cultivar) + (1|Location)+(1|Year)+(1|Cultivar:Year), data = Filtered_protein,REML=F) 
CM_0<- lmer(Protein ~ (1|Cultivar) + (1|Location)+(1|Cultivar:Year), data = Filtered_protein,REML=F) 
anova(CM_,CM_0)#Year has a significant effect, **

CM_<- lmer(Protein ~ (1|Location)+(1|Year) + (1|Cultivar), data = Filtered_protein,REML=F) 
CM_0<- lmer(Protein ~ (1|Location)+(1|Year), data = Filtered_protein,REML=F) 
anova(CM_,CM_0)#Cultivar has a significant effect, ***

####Fit model on protein yield
Filtered_protein$Proteinyield=(Filtered_protein$Protein/100)*Filtered_protein$Yield..g..Pr.kvm.
Environmental_proteinyield_mean=aggregate(Filtered_protein$Proteinyield,list(Filtered_protein$YearLoc),mean)
Filtered_protein$Environmental_proteinyield_mean=Filtered_protein$Proteinyield
Filtered_protein$Environmental_proteinyield_mean[which(Filtered_protein$YearLoc=="1.N")]=Environmental_proteinyield_mean[1,2]
Filtered_protein$Environmental_proteinyield_mean[which(Filtered_protein$YearLoc=="2.N")]=Environmental_proteinyield_mean[2,2]
Filtered_protein$Environmental_proteinyield_mean[which(Filtered_protein$YearLoc=="3.N")]=Environmental_proteinyield_mean[3,2]
Filtered_protein$Environmental_proteinyield_mean[which(Filtered_protein$YearLoc=="2.S")]=Environmental_proteinyield_mean[4,2]
Filtered_protein$Environmental_proteinyield_mean[which(Filtered_protein$YearLoc=="1.V")]=Environmental_proteinyield_mean[5,2]

Filtered_protein$Scaled_Environmental_proteinyield_mean=scale(Filtered_protein$Environmental_proteinyield_mean) #scale to get variance=1, mean=0

Model_proteinyield<- lmer(Proteinyield ~ (1|Location)+(1|Year)+(1|Cultivar:Location)+(1|Cultivar:Year)+(1|YearLoc)+ (Scaled_Environmental_protein_mean|Cultivar), data = Filtered_protein,REML=T) 
summary(Model_proteinyield)

#test if Finlay-Wilkinson regressions are significant 
Model_proteinyield_<- lmer(Proteinyield ~ (1|Location)+(1|Year)+(1|Cultivar:Location)+(1|Cultivar:Year)+(1|YearLoc)+ (Scaled_Environmental_protein_mean|Cultivar), data = Filtered_protein,REML=F) 
Model_proteinyield_0<- lmer(Proteinyield ~ (1|Location)+(1|Year)+(1|Cultivar:Location)+(1|Cultivar:Year)+(1|YearLoc)+ (1|Cultivar), data = Filtered_protein,REML=F) 
summary(Model_proteinyield_0)
anova(Model_proteinyield_,Model_proteinyield_0) #The GxE FW interaction of proteinyield is not significant, really.

#Produce Figure 3
library("varhandle")
ggplot(new, aes(x=unfactor(BLUPS_Yield), y=unfactor(BLUPS_Protein)))+geom_point(col="#00ba38",cex=3.5) +
  geom_text(label=new$cultivar,size=5.5,vjust = 0, nudge_y = 0.1) +
  stat_smooth(method="lm",col="#00ba38") +
  labs(x="Cultivar BLUPs Yield", y = "Cultivar BLUPs Protein") + theme_bw()



# Test significance of protein-yield cor. when leaving some out
correlation(as.numeric(as.character(new$BLUPS_Protein)),as.numeric(as.character(new$BLUPS_Yield))) 
correlation(Cultivar_effects_protein,Cultivar_effects) #BLUPs of yield on all data, -0.61 **

new_withoutkontuandsd=new[-c(10,14),]
ggplot(new_withoutkontuandsd, aes(x=unfactor(BLUPS_Yield), y=unfactor(BLUPS_Protein)))+geom_point(col="#00ba38",cex=3.5) +
  geom_text(label=new_withoutkontuandsd$cultivar,size=5.5,vjust = 0, nudge_y = 0.1) +
  stat_smooth(method="lm",col="#00ba38") + 
  labs(x="Cultivar BLUPs Yield", y = "Cultivar BLUPs Protein") + theme_bw()

correlation(as.numeric(as.character(new_withoutkontuandsd$BLUPS_Protein)),as.numeric(as.character(new_withoutkontuandsd$BLUPS_Yield))) #-0.63 *

new_withoutkontuandsdandssns1=new[-c(10,14,15),]
ggplot(new_withoutkontuandsdandssns1, aes(x=unfactor(BLUPS_Yield), y=unfactor(BLUPS_Protein)))+geom_point(col="#00ba38",cex=3.5) +
  geom_text(label=new_withoutkontuandsdandssns1$cultivar,size=5.5,vjust = 0, nudge_y = 0.1) +
  stat_smooth(method="lm",col="#00ba38") +
  labs(x="Cultivar BLUPs Yield", y = "Cultivar BLUPs Protein") + theme_bw()

correlation(as.numeric(as.character(new_withoutkontuandsdandssns1$BLUPS_Protein)),as.numeric(as.character(new_withoutkontuandsdandssns1$BLUPS_Yield))) #-0.54 *

new_withouthighprotein=new[-c(8,12,15),]
ggplot(new_withouthighprotein, aes(x=unfactor(BLUPS_Yield), y=unfactor(BLUPS_Protein)))+geom_point(col="#00ba38",cex=3.5) +
  geom_text(label=new_withouthighprotein$cultivar,size=5.5,vjust = 0, nudge_y = 0.1) +
  stat_smooth(method="lm",col="#00ba38") + 
  labs(x="Cultivar BLUPs Yield", y = "Cultivar BLUPs Protein") + theme_bw()

correlation(as.numeric(as.character(new_withouthighprotein$BLUPS_Protein)),as.numeric(as.character(new_withouthighprotein$BLUPS_Yield))) #-0.73 **


new_rm_5lowestyield=new[-c(10,14,15,8,12),]
ggplot(new_rm_5lowestyield, aes(x=unfactor(BLUPS_Yield), y=unfactor(BLUPS_Protein)))+geom_point(col="#00ba38",cex=3.5) +
  geom_text(label=new_rm_5lowestyield$cultivar,size=5.5,vjust = 0, nudge_y = 0.1) +
  stat_smooth(method="lm",col="#00ba38") +
  labs(x="Cultivar BLUPs Yield", y = "Cultivar BLUPs Protein") + theme_bw()

correlation(as.numeric(as.character(new_rm_5lowestyield$BLUPS_Protein)),as.numeric(as.character(new_rm_5lowestyield$BLUPS_Yield))) #0.14 (not signif)

lowyield=new[-c(1,3,2,4,5,6,7,9,11,13,16,17),]
ggplot(lowyield, aes(x=unfactor(BLUPS_Yield), y=unfactor(BLUPS_Protein)))+geom_point(col="#00ba38",cex=3.5) +
  geom_text(label=lowyield$cultivar,size=5.5,vjust = 0, nudge_y = 0.1) +
  stat_smooth(method="lm",col="#00ba38") +
  labs(x="Cultivar BLUPs Yield", y = "Cultivar BLUPs Protein") + theme_bw()

correlation(as.numeric(as.character(lowyield$BLUPS_Protein)),as.numeric(as.character(lowyield$BLUPS_Yield))) #0.97 (**)

## Make Figure 1, panel 1
CVpercent
total$CV=c(40.59527,37.12608,40.05453,41.90156,42.78845,41.12105,37.52563,45.72125,42.89822,44.01177,40.12525,39.70240,39.88650,36.86544,40.03375,39.39686,44.05776)
plot=ggplot(total, aes(x=MeanYield, y=CV))+geom_point(col="#00ba38",cex=3.5,alpha=0.8) +
  geom_text(label=total$Cultivar,size=4.5, vjust=0.5) +
  #scale_y_continuous(labels=scaleFUN) +
  #scale_x_continuous(labels=scaleFUN) +
  ylab("Coefficient of variance (%)") + ylim(36.5,47) +
  xlab("Mean yield g/sqm") +
  geom_vline(xintercept=mean(filtered$Yield..g..Pr.kvm.),colour="black",linetype=2) +
  geom_hline(yintercept=mean(total$CV),colour="black",linetype=2)

plot+theme_bw()

# Panel 1 Figure 1
e <- ggplot(newd, aes(x = reorder(Cultivar, Yield..g..Pr.kvm., fun = mean), y = Yield..g..Pr.kvm.)) + 
  geom_boxplot(fill = "#00ba38",alpha=.8) +geom_jitter(position=position_jitter(0.2),cex=1)  +
  labs(x="Cultivar", y="Seed yield (g/m2)")+theme_bw(base_size = 14) 

e


e1<-ggplot(data=total, aes(x=reorder(Cultivar, MeanYield), y=MeanYield)) +
  geom_bar(stat="identity", fill="#00ba38",alpha=.8)+ ylim()
  theme_minimal()
e1

### Calculate correlations between stability parameters
total$regression=c(0.96,0.91,0.97,1.04,1.12,1.14,1.06,1.02,1.12,0.72,1.13,0.92,1.11,0.62,0.90,1.01,1.19)

correlation(total$CV,total$MeanYield) # -0.06 not significant
correlation(total$regression,total$MeanYield) #0.91 ***

correlation(total$CV,total$regression) #0.33 (not significant)




###Calculate correlations between environmental mean yield and weather conditions. 
#For climatic table in article

climate=read.table("/Volumes/ST_MBG-PMg/Cathrine/NorFab/Article/3. Version March2020/Climatedata.txt",sep="\t",header=T)

correlation(climate$Days_of_growth,climate$Mean_yield) # 0.85, **
correlation(climate$Avg_Temp,climate$Mean_yield) # -0.69, *
correlation(climate$Avg_Prec,climate$Mean_yield) # 0.81, **
correlation(climate$Days_with_prec,climate$Mean_yield) # 0.81, **
correlation(climate$Avg_hum,climate$Mean_yield) # 0.64, not significant
correlation(climate$Avg_sun,climate$Mean_yield) # -0.75, *
correlation(climate$Dry_days,climate$Mean_yield) # -0.76, *
correlation(climate$Mean_yield,climate$Mean_yield) # 1, ***

correlation(climate$Days_of_growth,climate$Mean_protein) # 0.35, not significant
correlation(climate$Avg_Temp,climate$Mean_protein) # -0.25, not significant
correlation(climate$Avg_Prec,climate$Mean_protein) # 0.32, not significant
correlation(climate$Days_with_prec,climate$Mean_protein) # 0.05, not significant
correlation(climate$Avg_hum,climate$Mean_protein) # -0.31, not significant
correlation(climate$Avg_sun,climate$Mean_protein) # -0.08, not significant
correlation(climate$Dry_days,climate$Mean_protein) # 0, not significant
correlation(climate$Mean_yield,climate$Mean_protein) # 0.22, not significant



###Calculate correlations between traits of cultivars
#For cultivar description table in article
CultivarYieldMeans=aggregate(filtered$Yield..g..Pr.kvm.,list(filtered$Cultivar),mean)
CultivarProteinMeans=aggregate(Filtered_protein$Protein,list(Filtered_protein$Cultivar),mean)

sd_proteinmeans=aggregate(Filtered_protein$Protein,list(Filtered_protein$Cultivar),sd)
sd_proteinmeans
repcultivar_prot=count(Filtered_protein, vars = Cultivar)
se_ProteinMeans=sd_proteinmeans[,2]/sqrt(repcultivar_prot$n) #Se of cultivar means
round(se_ProteinMeans,1)

Filtered_protein$Proteinyield=Filtered_protein$Yield..g..Pr.kvm.*Filtered_protein$Protein/100
Proteinyieldmean=aggregate(Filtered_protein$Proteinyield,list(Filtered_protein$Cultivar),mean)


correlation(CultivarYieldMeans[,2],CultivarProteinMeans[,2]) #-0.68 **
correlation(CultivarYieldMeans[,2],Proteinyieldmean[,2]) #0.93 ***



# Panel 2 Supporting Figure 1

CV_mean=tapply(Filtered_protein$Protein, Filtered_protein$Cultivar, mean)

CV_stdev=tapply(Filtered_protein$Protein, Filtered_protein$Cultivar, sd)

CV_var=tapply(Filtered_protein$Protein, Filtered_protein$Cultivar, var)

CVpercent=CV_stdev/CV_mean
CVpercent
CVs=c(4.611475,3.108916,3.692058,4.954016,3.648544,3.421201,4.335544,3.965863,4.730999,4.875816,3.645966,6.554386,5.334901,4.228979,7.383690,3.549022,5.177707)
total1=cbind(CultivarProteinMeans,CVs)
colnames(total1)[2]="MeanProtein"
colnames(total1)[1]="Cultivar"

plot=ggplot(total1, aes(x=MeanProtein, y=CVs))+geom_point(col="#00ba38",cex=3.5,alpha=0.8) +
  geom_text(label=total1$Cultivar,size=4.5, vjust=0.5) +
  #scale_y_continuous(labels=scaleFUN) +
  #scale_x_continuous(labels=scaleFUN) +
  ylab("Coefficient of variation (%)") + ylim(2,8) +
  xlab("Protein content (%)") +
  geom_vline(xintercept=mean(Filtered_protein$Protein),colour="black",linetype=2) +
  geom_hline(yintercept=mean(total1$CVs),colour="black",linetype=2)

plot+theme_bw()

# Save as default size suggested

# Panel 1 Supporting Figure 1
e <- ggplot(Filtered_protein, aes(x = reorder(Cultivar, Protein, fun = mean), y = Protein)) + 
  geom_boxplot(fill = "#00ba38",alpha=.8) +geom_jitter(position=position_jitter(0.2),cex=1)  +
  labs(x="Cultivar", y="Protein content (%)")+theme_bw(base_size = 14) 

e #Save as default size suggested


# compare cultivar means for protein 
setwd("/Volumes/ST_MBG-PMg/Cathrine/NorFab/Article/1. version/Publication ready/Final_20191112")
data=read.table("YieldData_nodash.txt",sep="\t",header=T)
rowstoremove=which(is.na(data$Protein))
data=data[-rowstoremove,]

fit <- aov(Protein~ factor(Cultivar)  , data=data)

results <- TukeyHSD(fit, ordered=TRUE)
multcompLetters4(fit, results)
plot(TukeyHSD(fit))








Description=read.table("/Volumes/ST_MBG-PMg/Cathrine/NorFab/Article/3. Version March2020/Genotype_information.txt",sep="\t",header=T)
Genotypeinformation=cbind(CultivarYieldMeans,Description[,2:6])
colnames(Genotypeinformation)[1:2]=c("Cultivar","Yield_mean")
head(Genotypeinformation)
Genotypeinformation=cbind(Genotypeinformation,CultivarProteinMeans[,2])
colnames(Genotypeinformation)[8]=c("Protein_mean")


correlation(Genotypeinformation$Yield_mean,Genotypeinformation$TGW) #0.81, ***
correlation(Genotypeinformation$Yield_mean,Genotypeinformation$Days_to_Maturation) #0.68, **
correlation(Genotypeinformation$Yield_mean,Genotypeinformation$Days_to_flowering) #0.09961, not significant
correlation(Genotypeinformation$Yield_mean,Genotypeinformation$Days_of_flowering) #-0.22, not significant
correlation(Genotypeinformation$Yield_mean,Genotypeinformation$zt) #-0.37, not significant


correlation(Genotypeinformation$Protein_mean,Genotypeinformation$TGW) #-0.60, *
correlation(Genotypeinformation$Protein_mean,Genotypeinformation$Days_to_Maturation) #-0.23, not significant
correlation(Genotypeinformation$Protein_mean,Genotypeinformation$Days_to_flowering) #0.075, not significant
correlation(Genotypeinformation$Protein_mean,Genotypeinformation$Days_of_flowering) #-0.08, not significant
correlation(Genotypeinformation$Protein_mean,Genotypeinformation$zt) #0.59, *

correlation(Genotypeinformation$Yield_mean,Genotypeinformation$Protein_mean) #-0.68, **



#Anova test to compare seed yield means 
setwd("/Volumes/ST_MBG-PMg/Cathrine/NorFab/Article/1. version/Publication ready/Final_20191112")
data=read.table("YieldData_nodash.txt",sep="\t",header=T)
rowstoremove=which(is.na(data$Yield..kg..Pr.kvm.))
data=data[-rowstoremove,]

data$Yield..g..Pr.kvm.=data$Yield..kg..Pr.kvm.*1000

fit <- aov(Yield..g..Pr.kvm.~ factor(Cultivar)  , data=data)

results <- TukeyHSD(fit, ordered=TRUE)
multcompLetters4(fit, results)
plot(TukeyHSD(fit))

#n*=29

#Mean comparison climate data
climatedata=read.table("/Volumes/ST_MBG-PMg/Cathrine/NorFab/Article/3. Version March2020/Climate_data.txt",head=T)
head(climatedata)

temp <- aov(Temp~ factor(TRID)  , data=climatedata)
resultstemp <- TukeyHSD(temp, ordered=TRUE)
multcompLetters4(temp, resultstemp)
plot(TukeyHSD(temp))

temp_year <- aov(Temp~ factor(year)  , data=climatedata)
resultstemp_temp <- TukeyHSD(temp_year, ordered=TRUE)
multcompLetters4(temp_year, resultstemp_temp)
plot(TukeyHSD(temp_year))

temp_loc <- aov(Temp~ factor(loc)  , data=climatedata)
temp_loc_res <- TukeyHSD(temp_loc, ordered=TRUE)
multcompLetters4(temp_loc, temp_loc_res)
plot(TukeyHSD(temp_loc))

prec <- aov(Prec~ factor(TRID)  , data=climatedata)
resultsPrec <- TukeyHSD(prec, ordered=TRUE)
multcompLetters4(prec, resultsPrec)
plot(TukeyHSD(prec))

prec_year <- aov(Prec~ factor(year)  , data=climatedata)
resultsPrec_year <- TukeyHSD(prec_year, ordered=TRUE)
multcompLetters4(prec_year, resultsPrec_year)
plot(TukeyHSD(prec_year))

prec_loc <- aov(Prec~ loc  , data=climatedata)
resultsPrec_loc <- TukeyHSD(prec_loc, ordered=TRUE)
multcompLetters4(prec_loc, resultsPrec_loc)
plot(TukeyHSD(prec_loc))


Humidity <- aov(Humidity~ factor(TRID)  , data=climatedata)
resultsHum <- TukeyHSD(Humidity, ordered=TRUE)
multcompLetters4(Humidity, resultsHum)
plot(TukeyHSD(Humidity))

Humidity_year <- aov(Humidity~ factor(year)  , data=climatedata)
Res_Humidity_year <- TukeyHSD(Humidity_year, ordered=TRUE)
multcompLetters4(Humidity_year, Res_Humidity_year)
plot(TukeyHSD(Humidity_year))

Humidity_loc <- aov(Humidity~ factor(loc)  , data=climatedata)
Res_Humidity_loc <- TukeyHSD(Humidity_loc, ordered=TRUE)
multcompLetters4(Humidity_loc, Res_Humidity_loc)
plot(TukeyHSD(Humidity_loc))

sunshineDur <- aov(sunshineDur~ factor(TRID)  , data=climatedata)
resultsSun <- TukeyHSD(sunshineDur, ordered=TRUE)
multcompLetters4(sunshineDur, resultsSun)
plot(TukeyHSD(sunshineDur))

sunshineDur_year <- aov(sunshineDur~ factor(year)  , data=climatedata)
sunshineDur_year_results <- TukeyHSD(sunshineDur_year, ordered=TRUE)
multcompLetters4(sunshineDur_year, sunshineDur_year_results)
plot(TukeyHSD(sunshineDur_year))

sunshineDur_loc <- aov(sunshineDur~ factor(loc)  , data=climatedata)
sunshineDur_loc_results <- TukeyHSD(sunshineDur_loc, ordered=TRUE)
multcompLetters4(sunshineDur_loc, sunshineDur_loc_results)
plot(TukeyHSD(sunshineDur_loc))
