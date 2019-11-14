
###################################################
#     Initial yield analysis NorFab 2019          #
###################################################

setwd("/Volumes/ST_MBG-PMg/Cathrine/NorFab/Phenotypic_analysis/Yield_Stability_Analyses_20190610")

library(lme4)
library(nlme)
library(ggplot2)
library(gridExtra)
library("RLRsim")
library("agricolae")
library(afex)
library("cowplot")
library("dmm")

scaleFUN <- function(x) sprintf("%.2f", x)


# Load data
d=read.table("YieldData_short.txt",sep="\t",header=T)
head(d)
d$Yield..g..Pr.kvm.=d$Yield..kg..Pr.kvm.*1000
d$Log_Yield=log(d$Yield..g..Pr.kvm.)
head(d)
# Variance by cultivar name
# We do not differentiate between locations and years
completeFun <- function(data, desiredCols) {
  completeVec <- complete.cases(data[, desiredCols])
  return(data[completeVec, ])
}
newd=completeFun(d,"Yield..g..Pr.kvm.")

variances=aggregate(newd$Log_Yield,list(newd$Cultivar),var)
mean=aggregate(newd$Yield..g..Pr.kvm.,list(newd$Cultivar),mean)
total=cbind(variances,mean)
colnames(total)=c("Cultivar","YieldVar","Cultivar2","MeanYield")


# Produce plot for Figure 1A
plot=ggplot(total, aes(x=MeanYield, y=YieldVar))+geom_point(col="#00ba38",cex=3.5,alpha=0.8) +
  geom_text(label=total$Cultivar,size=4.5, vjust=0.5) +
  #scale_y_continuous(labels=scaleFUN) +
  #scale_x_continuous(labels=scaleFUN) +
  ylab("Variance of log(yield)") +
  xlab("Mean yield g/sqm") +
  ylim(c(0.18, 0.275)) +
  xlim(c(275, 500)) +
  geom_vline(xintercept=mean(total$Yield..g..Pr.kvm.),colour="black",linetype=2) +
  geom_hline(yintercept=mean(total$YieldVar),colour="black",linetype=2)

plot+theme_bw()

# Remove plants with no observations

rowstoremove=which(is.na(d$Yield..g..Pr.kvm.))
filtered=d[-rowstoremove,]


#### Analyses


#check if a year has very low yield or very high yield.
#If this is the case it will be skew the overall mean and possibly reduce the crop genotype yield rankings

#there are different means for each year and for each locations
  
ggplot(filtered, aes(x=Location, y=Yield..g..Pr.kvm.,group=Location)) + 
  geom_boxplot(fill="steelblue", alpha=0.5)+ geom_jitter(shape=16, position=position_jitter(0.2))+theme_classic() +
  labs(title="Yield in different locations",x="Location", y = "Yield g/sqrm") +
  theme(axis.text.x = element_text(angle = 90)) +
  geom_hline(yintercept=mean(na.omit(filtered$Yield..g..Pr.kvm.)), linetype="dashed", 
             color = "black", linetype=2)

ggplot(filtered, aes(x=Year, y=Yield..g..Pr.kvm.,group=Year)) + 
  geom_boxplot(fill="steelblue", alpha=0.5) + geom_jitter(shape=16, position=position_jitter(0.2))+theme_classic() +
  labs(title="Yield in different years",x="Year", y = "Yield g/sqrm") +
  theme(axis.text.x = element_text(angle = 90)) +
  geom_hline(yintercept=mean(na.omit(filtered$Yield..g..Pr.kvm.)), linetype="dashed", 
             color = "black", linetype=2)

ggplot(filtered, aes(x=Year, y=Yield..g..Pr.kvm.,group=Year)) + 
  geom_boxplot(fill="steelblue", alpha=0.5) + geom_jitter(shape=16, position=position_jitter(0.2),aes(colour=Location))+theme_classic() +
  labs(title="Yield in different years",x="Year", y = "Yield g/sqrm") +
  theme(axis.text.x = element_text(angle = 90)) +
  geom_hline(yintercept=mean(na.omit(filtered$Yield..g..Pr.kvm.)), linetype="dashed", 
             color = "black", linetype=2)



# Make a Year-location interaction column
filtered$YearLoc=interaction(filtered$Year,filtered$Location)


##########################################
# Calculating average neighbour values   #
##########################################

for (i in 1:nrow(filtered)){
  nabo=which(filtered$Plot==filtered$Ne1[i])
  if (length(nabo)==1){
    filtered$Ne1value[i]=filtered$Yield..g..Pr.kvm.[nabo]
    
  }
  else{
    filtered$Ne1value[i]=NA
  }
}


for (i in 1:nrow(filtered)){
  nabo=which(filtered$Plot==filtered$Ne2[i])
  if (length(nabo)==1){
    filtered$Ne2value[i]=filtered$Yield..g..Pr.kvm.[nabo]
    
  }
  else{
    filtered$Ne2value[i]=NA
  }
}


for (i in 1:nrow(filtered)){
  nabo=which(filtered$Plot==filtered$Ne3[i])
  if (length(nabo)==1){
    filtered$Ne3value[i]=filtered$Yield..g..Pr.kvm.[nabo]
  }
  else{
    filtered$Ne3value[i]=NA
  }
}

for (i in 1:nrow(filtered)){
  nabo=which(filtered$Plot==filtered$Ne4[i])
  if (length(nabo)==1){
    filtered$Ne4value[i]=filtered$Yield..g..Pr.kvm.[nabo]
  }
  else{
    filtered$Ne4value[i]=NA
  }
}

for (i in 1:nrow(filtered)){
  nabo=which(filtered$Plot==filtered$Ne5[i])
  if (length(nabo)==1){
    filtered$Ne5value[i]=filtered$Yield..g..Pr.kvm.[nabo]
  }
  else{
    filtered$Ne5value[i]=NA
  }
}


for (i in 1:nrow(filtered)){
  nabo=which(filtered$Plot==filtered$Ne6[i])
  if (length(nabo)==1){
    filtered$Ne6value[i]=filtered$Yield..g..Pr.kvm.[nabo]
  }
  else{
    filtered$Ne6value[i]=NA
  }
}

for (i in 1:nrow(filtered)){
  nabo=which(filtered$Plot==filtered$Ne7[i])
  if (length(nabo)==1){
    filtered$Ne7value[i]=filtered$Yield..g..Pr.kvm.[nabo]
  }
  else{
    filtered$Ne7value[i]=NA
  }
}

for (i in 1:nrow(filtered)){
  nabo=which(filtered$Plot==filtered$Ne8[i])
  if (length(nabo)==1){
    filtered$Ne8value[i]=filtered$Yield..g..Pr.kvm.[nabo]
  }
  else{
    filtered$Ne8value[i]=NA
  }
}


filtered$Average_neighbour_score=rowMeans(filtered[,21:28],na.rm =TRUE)



#Fitting LMM on seed yield
M6<- lmer(Yield..g..Pr.kvm. ~ (1|Cultivar)+(1|Location)+(1|Year)+Average_neighbour_score+(1|Cultivar:Location)+(1|Cultivar:Year), data = filtered,REML=T) 
summary(M6)

Residuals <- residuals(M6)
shapiro.test(Residuals) #Fits
qqnorm(Residuals); qqline(Residuals)
plot(M6, abline=c(0, 0))

lattice::qqmath(M6,id=0.1,idLabels=~.obs)

anova(M6) #This F score corresponds to neighbour having a significant effect
VarCorr(M6)




# Calculating BLUPs for all random effect
library("nlme")
Cultivar_effects=ranef(M6)$Cultivar 

Highyield=subset(Cultivar_effects,Cultivar_effects>0) #Subset cultivars with positive yield effect

Variances=tapply(filtered$Log_Yield, filtered$Cultivar, var)
Variances=data.matrix(Variances)
Variances=Variances[-6,1]
Variances=data.matrix(Variances)
Variances=Variances[-13,1]
Variances=data.matrix(Variances)


############ Calculate CV%, coefficient of variation 

CV_mean=tapply(filtered$Yield..g..Pr.kvm., filtered$Cultivar, mean)
CV_mean=CV_mean[-6]
CV_mean=CV_mean[-10]
CV_mean=CV_mean[-12]

CV_stdev=tapply(filtered$Yield..g..Pr.kvm., filtered$Cultivar, sd)
CV_stdev=CV_stdev[-6]
CV_stdev=CV_stdev[-10]
CV_stdev=CV_stdev[-12]
CV_var=tapply(filtered$Yield..g..Pr.kvm., filtered$Cultivar, var)

CV_var=CV_var[-6]
CV_var=CV_var[-10]
CV_var=CV_var[-12]

CVpercent=CV_stdev/CV_mean
CVpercent




### Finlay-Wilkinson regression, produce Figure 2

#Yield means
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

means_pr_trial=c(mean5,mean7,mean9,mean3,mean1,mean8,mean6,mean4,mean2)
means_pr_trial=as.data.frame(means_pr_trial)
rownames(means_pr_trial)=c("Sejet_J-2018","Sejet_F-2018","Finland-2018","Nordic_Seed-2018","Nordic_Seed-2016","Finland-2016","Sejet_F-2017","Sejet_J-2017","Nordic_Seed-2017")

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


All_years=Reduce(function(x, y) merge(x, y, by="Group.1", all=TRUE), list(S2018_CultivarMeans, SF2018_CultivarMeans, V2018_CultivarMeans,NS2018_CultivarMeans,NS2016_CultivarMeans,V2016_CultivarMeans,SF2017_CultivarMeans,S2017_CultivarMeans,NS2017_CultivarMeans))
colnames(All_years)=c("Cultivar","S2018","SF2018","V2018","NS2018","NS2016","V2016","SF2017","S2017","NS2017")
transpose=t(All_years)

means_pr_trial

merged=cbind(means_pr_trial,transpose[2:nrow(transpose),])
colnames(merged)=c("means_pr_trial",transpose[1,])
head(merged)


library("cowplot")
pdf(file="YieldStabilityAllPlots_1.pdf",width=65,height=75,useDingbats = F)

lm_fit=lm(unfactor(`247-13`)~means_pr_trial,data=merged)
summary(lm_fit)
p1=ggplot(merged,aes(means_pr_trial,unfactor(`247-13`)))+geom_point(col="#00ba38",cex=10,alpha=0.8) + ylim(c(100, 1000)) + xlim(c(100,1000)) + ggtitle("Cultivar 247-13") +geom_abline(intercept = 0, slope = 1,col="black",cex=2) + geom_abline(col="#00ba38",intercept=4.75409,slope=0.96696,cex=4) +
  theme(axis.text.x=element_text(size=55),
        axis.text.y=element_text(size=55),
        plot.title = element_text(size=75),) +
  geom_text(size=20,aes(x=850,y=150,hjust=0,vjust=0,label="0.97"))+ theme(
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
p2=ggplot(merged,aes(means_pr_trial,unfactor(`749-13`)))+geom_point(col="#00ba38",cex=10,alpha=0.8) + ylim(c(100, 1000)) + xlim(c(100,1000)) + ggtitle("Cultivar 749-13") +geom_abline(intercept = 0, slope = 1,col="black",cex=2) + geom_abline(col="#00ba38",intercept=37.44494,slope=0.91350,cex=4) +
  theme(axis.text.x=element_text(size=55),
        axis.text.y=element_text(size=55),
        plot.title = element_text(size=75),) +
  geom_text(size=20,aes(x=850,y=150,hjust=0,vjust=0,label="0.91"))+ theme(
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
p3=ggplot(merged,aes(means_pr_trial,unfactor(`Alexia`)))+geom_point(col="#00ba38",cex=10,alpha=0.8) + ylim(c(100, 1000)) + xlim(c(100,1000)) + ggtitle("Cultivar Alexia") +geom_abline(intercept = 0, slope = 1,col="black",cex=2) + geom_abline(col="#00ba38",intercept=10.86308,slope=0.97304,cex=4) +
  theme(axis.text.x=element_text(size=55),
        axis.text.y=element_text(size=55),
        plot.title = element_text(size=75),) +
  geom_text(size=20,aes(x=850,y=150,hjust=0,vjust=0,label="0.97"))+ theme(
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
p4=ggplot(merged,aes(means_pr_trial,unfactor(`Banquise`)))+geom_point(col="#00ba38",cex=10,alpha=0.8) + ylim(c(100, 1000)) + xlim(c(100,1000)) + ggtitle("Cultivar Banquise") +geom_abline(intercept = 0, slope = 1,col="black",cex=2) + geom_abline(col="#00ba38",intercept=3.9112,slope=1.0394,cex=4) +
  theme(axis.text.x=element_text(size=55),
        axis.text.y=element_text(size=55),
        plot.title = element_text(size=75),) +
  geom_text(size=20,aes(x=850,y=150,hjust=0,vjust=0,label="1.04"))+ theme(
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
p5=ggplot(merged,aes(means_pr_trial,unfactor(`Boxer`)))+geom_point(col="#00ba38",cex=10,alpha=0.8) + ylim(c(100, 1000)) + xlim(c(100,1000)) + ggtitle("Cultivar Boxer") +geom_abline(intercept = 0, slope = 1,col="black",cex=2) + geom_abline(col="#00ba38",intercept=-24.78261,slope=1.11656,cex=4) +
  theme(axis.text.x=element_text(size=55),
        axis.text.y=element_text(size=55),
        plot.title = element_text(size=75),) +
  geom_text(size=20,aes(x=850,y=150,hjust=0,vjust=0,label="1.12"))+ theme(
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
p6=ggplot(merged,aes(means_pr_trial,unfactor(`Fanfare`)))+geom_point(col="#00ba38",cex=10,alpha=0.8) + ylim(c(100, 1000)) + xlim(c(100,1000)) + ggtitle("Cultivar Fanfare") +geom_abline(intercept = 0, slope = 1,col="black",cex=2) + geom_abline(col="#00ba38",intercept=-9.77691,slope=1.15259,cex=4) +
  theme(axis.text.x=element_text(size=55),
        axis.text.y=element_text(size=55),
        plot.title = element_text(size=75),) +
  geom_text(size=20,aes(x=850,y=150,hjust=0,vjust=0,label="1.15"))+ theme(
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
p7=ggplot(merged,aes(means_pr_trial,unfactor(`Fuego`)))+geom_point(col="#00ba38",cex=10,alpha=0.8) + ylim(c(100, 1000)) + xlim(c(100,1000)) + ggtitle("Cultivar Fuego") +geom_abline(intercept = 0, slope = 1,col="black",cex=2) + geom_abline(col="#00ba38",intercept=15.77001,slope=1.06541,cex=4) +
  theme(axis.text.x=element_text(size=55),
        axis.text.y=element_text(size=55),
        plot.title = element_text(size=75),) +
  geom_text(size=20,aes(x=850,y=150,hjust=0,vjust=0,label="1.07"))+ theme(
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
p8=ggplot(merged,aes(means_pr_trial,unfactor(`Gloria`)))+geom_point(col="#00ba38",cex=10,alpha=0.8) + ylim(c(100, 1000)) + xlim(c(100,1000)) + ggtitle("Cultivar Gloria") +geom_abline(intercept = 0, slope = 1,col="black",cex=2) + geom_abline(col="#00ba38",intercept=-47.47847,slope=1.03937,cex=4) +
  theme(axis.text.x=element_text(size=55),
        axis.text.y=element_text(size=55),
        plot.title = element_text(size=75),) +
  geom_text(size=20,aes(x=850,y=150,hjust=0,vjust=0,label="1.04"))+ theme(
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
p9=ggplot(merged,aes(means_pr_trial,unfactor(`Gracia`)))+geom_point(col="#00ba38",cex=10,alpha=0.8) + ylim(c(100, 1000)) + xlim(c(100,1000)) + ggtitle("Cultivar Gracia") +geom_abline(intercept = 0, slope = 1,col="black",cex=2) + geom_abline(col="#00ba38",intercept=-29.57231,slope=1.11590,cex=4) +
  theme(axis.text.x=element_text(size=55),
        axis.text.y=element_text(size=55),
        plot.title = element_text(size=75),) +
  geom_text(size=20,aes(x=850,y=150,hjust=0,vjust=0,label="1.12"))+ theme(
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
p10=ggplot(merged,aes(means_pr_trial,unfactor(`Kontu`)))+geom_point(col="#00ba38",cex=10,alpha=0.8) + ylim(c(100, 1000)) + xlim(c(100,1000)) + ggtitle("Cultivar Kontu") +geom_abline(intercept = 0, slope = 1,col="black",cex=2) + geom_abline(col="#00ba38",intercept=0.8570,slope=0.6969,cex=4) +
  theme(axis.text.x=element_text(size=55),
        axis.text.y=element_text(size=55),
        plot.title = element_text(size=75),) +
  geom_text(size=20,aes(x=850,y=150,hjust=0,vjust=0,label="0.70"))+ theme(
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
p11=ggplot(merged,aes(means_pr_trial,unfactor(`Lynx`)))+geom_point(col="#00ba38",cex=10,alpha=0.8) + ylim(c(100, 1000)) + xlim(c(100,1000)) + ggtitle("Cultivar Lynx") +geom_abline(intercept = 0, slope = 1,col="black",cex=2) + geom_abline(col="#00ba38",intercept=6.55917,slope=1.12644,cex=4) +
  theme(axis.text.x=element_text(size=55),
        axis.text.y=element_text(size=55),
        plot.title = element_text(size=75),) +
  geom_text(size=20,aes(x=850,y=150,hjust=0,vjust=0,label="1.13"))+ theme(
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
p12=ggplot(merged,aes(means_pr_trial,unfactor(`Mistral`)))+geom_point(col="#00ba38",cex=10,alpha=0.8) + ylim(c(100, 1000)) + xlim(c(100,1000)) + ggtitle("Cultivar Mistral") +geom_abline(intercept = 0, slope = 1,col="black",cex=2) + geom_abline(col="#00ba38",intercept=11.13386,slope=0.92487,cex=4) +
  theme(axis.text.x=element_text(size=55),
        axis.text.y=element_text(size=55),
        plot.title = element_text(size=75),) +
  geom_text(size=20,aes(x=850,y=150,hjust=0,vjust=0,label="0.92"))+ theme(
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
p13=ggplot(merged,aes(means_pr_trial,unfactor(`Pyramid`)))+geom_point(col="#00ba38",cex=10,alpha=0.8) + ylim(c(100, 1000)) + xlim(c(100,1000)) + ggtitle("Cultivar Pyramid") +geom_abline(intercept = 0, slope = 1,col="black",cex=2) + geom_abline(col="#00ba38",intercept=4.27386,slope=1.10906,cex=4) +
  theme(axis.text.x=element_text(size=55),
        axis.text.y=element_text(size=55),
        plot.title = element_text(size=75),) +
  geom_text(size=20,aes(x=850,y=150,hjust=0,vjust=0,label="1.11"))+ theme(
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
p14=ggplot(merged,aes(means_pr_trial,unfactor(`Snowdrop`)))+geom_point(col="#00ba38",cex=10,alpha=0.8) + ylim(c(100, 1000)) + xlim(c(100,1000)) + ggtitle("Cultivar Snowdrop") +geom_abline(intercept = 0, slope = 1,col="black",cex=2) + geom_abline(col="#00ba38",intercept=66.55529,slope=0.59872,cex=4) +
  theme(axis.text.x=element_text(size=55),
        axis.text.y=element_text(size=55),
        plot.title = element_text(size=75),) +
  geom_text(size=20,aes(x=850,y=150,hjust=0,vjust=0,label="0.60"))+ theme(
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
p15=ggplot(merged,aes(means_pr_trial,unfactor(`SSNS-1`)))+geom_point(col="#00ba38",cex=10,alpha=0.8) + ylim(c(100, 1000)) + xlim(c(100,1000)) + ggtitle("Cultivar SSNS-1") +geom_abline(intercept = 0, slope = 1,col="black",cex=2) + geom_abline(col="#00ba38",intercept=3.3778,slope=0.8953,cex=4) +
  theme(axis.text.x=element_text(size=55),
        axis.text.y=element_text(size=55),
        plot.title = element_text(size=75),) +
  geom_text(size=20,aes(x=850,y=150,hjust=0,vjust=0,label="0.90"))+ theme(
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
p16=ggplot(merged,aes(means_pr_trial,unfactor(`Taifun`)))+geom_point(col="#00ba38",cex=10,alpha=0.8) + ylim(c(100, 1000)) + xlim(c(100,1000)) + ggtitle("Cultivar Taifun") +geom_abline(intercept = 0, slope = 1,col="black",cex=2) + geom_abline(col="#00ba38",intercept=13.86413,slope=1.01178,cex=4) +
  theme(axis.text.x=element_text(size=55),
        axis.text.y=element_text(size=55),
        plot.title = element_text(size=75),) +
  geom_text(size=20,aes(x=850,y=150,hjust=0,vjust=0,label="1.01"))+ theme(
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
p17=ggplot(merged,aes(means_pr_trial,unfactor(`Vertigo`)))+geom_point(col="#00ba38",cex=10,alpha=0.8) + ylim(c(100, 1000)) + xlim(c(100,1000)) + ggtitle("Cultivar Vertigo") +geom_abline(intercept = 0, slope = 1,col="black",cex=2) + geom_abline(col="#00ba38",intercept=-34.30419,slope=1.18961,cex=4) +
  theme(axis.text.x=element_text(size=55),
        axis.text.y=element_text(size=55),
        plot.title = element_text(size=75),) +
  geom_text(size=20,aes(x=850,y=150,hjust=0,vjust=0,label="1.19"))+ theme(
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




### Yield protein variance, producing Figure 4

#Yield*protein means
NS2016$yield_x_protein=NS2016$Yield..g..Pr.kvm.*NS2016$Protein/100
mean1=mean(NS2016$yield_x_protein)
NS2017$yield_x_protein=NS2017$Yield..g..Pr.kvm.*NS2017$Protein/100
mean2=mean(NS2017$yield_x_protein)
NS2018$yield_x_protein=NS2018$Yield..g..Pr.kvm.*NS2018$Protein/100
mean3=mean(NS2018$yield_x_protein)
S2017$yield_x_protein=S2017$Yield..g..Pr.kvm.*S2017$Protein/100
mean4=mean(S2017$yield_x_protein)
V2016$yield_x_protein=V2016$Yield..g..Pr.kvm.*V2016$Protein/100
mean5=mean(V2016$yield_x_protein)



means_pr_trial=c(mean3,mean1,mean5,mean4,mean2)
means_pr_trial=as.data.frame(means_pr_trial)
rownames(means_pr_trial)=c("Nordic_Seed-2018","Nordic_Seed-2016","FN_2016","Sejet_J-2017","Nordic_Seed-2017")

# Calculate mean of every cultivar for each trial
NS2016_CultivarMeans=aggregate(NS2016$yield_x_protein,list(NS2016$Cultivar),mean)
NS2016_CultivarMeans=as.data.frame(NS2016_CultivarMeans)
NS2017_CultivarMeans=aggregate(NS2017$yield_x_protein,list(NS2017$Cultivar),mean)
NS2017_CultivarMeans=as.data.frame(NS2017_CultivarMeans)
NS2018_CultivarMeans=aggregate(NS2018$yield_x_protein,list(NS2018$Cultivar),mean)
NS2018_CultivarMeans=as.data.frame(NS2018_CultivarMeans)
S2017_CultivarMeans=aggregate(S2017$yield_x_protein,list(S2017$Cultivar),mean)
S2017_CultivarMeans=as.data.frame(S2017_CultivarMeans)
V2016_CultivarMeans=aggregate(V2016$yield_x_protein,list(V2016$Cultivar),mean)
V2016_CultivarMeans=as.data.frame(V2016_CultivarMeans)


All_years=Reduce(function(x, y) merge(x, y, by="Group.1", all=TRUE), list(NS2018_CultivarMeans,NS2016_CultivarMeans,V2016_CultivarMeans,S2017_CultivarMeans,NS2017_CultivarMeans))
colnames(All_years)=c("Cultivar","NS2018","V2016","NS2016","S2017","NS2017")
transpose=t(All_years)

means_pr_trial


merged=cbind(means_pr_trial,transpose[2:nrow(transpose),])
colnames(merged)=c("means_pr_trial",transpose[1,])
head(merged)

library("cowplot")
pdf(file="ProteinYieldStabilityAllPlots_1_20191114.pdf",width=65,height=75,useDingbats = F)


lm_fit=lm(unfactor(`247-13`)~means_pr_trial,data=merged)
summary(lm_fit)
p1=ggplot(merged,aes(means_pr_trial,unfactor(`247-13`)))+geom_point(col="#00ba38",cex=10,alpha=0.8) + ylim(c(0, 300)) + xlim(c(0,300)) + ggtitle("247-13") +geom_abline(intercept = 0, slope = 1,col="black",cex=2) + geom_abline(col="#00ba38",intercept=-11.07123,slope=1.03384,cex=4) +
  theme(axis.text.x=element_text(size=55),
        axis.text.y=element_text(size=55),
        plot.title = element_text(size=75),) +
  geom_text(size=20,aes(x=250,y=25,hjust=0,vjust=0,label="1.03"))+ theme(
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
p2=ggplot(merged,aes(means_pr_trial,unfactor(`749-13`)))+geom_point(col="#00ba38",cex=10,alpha=0.8) + ylim(c(0, 300)) + xlim(c(0,300)) + ggtitle("749-13") +geom_abline(intercept = 0, slope = 1,col="black",cex=2) + geom_abline(col="#00ba38",intercept=16.90689,slope=0.84496,cex=4) +
  theme(axis.text.x=element_text(size=55),
        axis.text.y=element_text(size=55),
        plot.title = element_text(size=75),) +
  geom_text(size=20,aes(x=250,y=25,hjust=0,vjust=0,label="0.84"))+ theme(
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
p3=ggplot(merged,aes(means_pr_trial,unfactor(`Alexia`)))+geom_point(col="#00ba38",cex=10,alpha=0.8) + ylim(c(0, 300)) + xlim(c(0,300)) + ggtitle("Alexia") +geom_abline(intercept = 0, slope = 1,col="black",cex=2) + geom_abline(col="#00ba38",intercept=-3.8225 ,slope=1.0253,cex=4) +
  theme(axis.text.x=element_text(size=55),
        axis.text.y=element_text(size=55),
        plot.title = element_text(size=75),) +
  geom_text(size=20,aes(x=250,y=25,hjust=0,vjust=0,label="1.03"))+ theme(
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
p4=ggplot(merged,aes(means_pr_trial,unfactor(`Banquise`)))+geom_point(col="#00ba38",cex=10,alpha=0.8) + ylim(c(0, 300)) + xlim(c(0,300)) + ggtitle("Banquise") +geom_abline(intercept = 0, slope = 1,col="black",cex=2) + geom_abline(col="#00ba38",intercept=-15.5672,slope=1.0626,cex=4) +
  theme(axis.text.x=element_text(size=55),
        axis.text.y=element_text(size=55),
        plot.title = element_text(size=75),) +
  geom_text(size=20,aes(x=250,y=25,hjust=0,vjust=0,label="1.06"))+ theme(
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
p5=ggplot(merged,aes(means_pr_trial,unfactor(`Boxer`)))+geom_point(col="#00ba38",cex=10,alpha=0.8) + ylim(c(0, 300)) + xlim(c(0,300)) + ggtitle("Boxer") +geom_abline(intercept = 0, slope = 1,col="black",cex=2) + geom_abline(col="#00ba38",intercept=-25.2708,slope=1.1825,cex=4) +
  theme(axis.text.x=element_text(size=55),
        axis.text.y=element_text(size=55),
        plot.title = element_text(size=75),) +
  geom_text(size=20,aes(x=250,y=25,hjust=0,vjust=0,label="1.18"))+ theme(
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
p6=ggplot(merged,aes(means_pr_trial,unfactor(`Fanfare`)))+geom_point(col="#00ba38",cex=10,alpha=0.8) + ylim(c(0, 300)) + xlim(c(0,300)) + ggtitle("Fanfare") +geom_abline(intercept = 0, slope = 1,col="black",cex=2) + geom_abline(col="#00ba38",intercept=-2.9345  ,slope=1.1417  ,cex=4) +
  theme(axis.text.x=element_text(size=55),
        axis.text.y=element_text(size=55),
        plot.title = element_text(size=75),) +
  geom_text(size=20,aes(x=250,y=25,hjust=0,vjust=0,label="1.14"))+ theme(
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
p7=ggplot(merged,aes(means_pr_trial,unfactor(`Fuego`)))+geom_point(col="#00ba38",cex=10,alpha=0.8) + ylim(c(0, 300)) + xlim(c(0,300)) + ggtitle("Fuego") +geom_abline(intercept = 0, slope = 1,col="black",cex=2) + geom_abline(col="#00ba38",intercept=18.6347,slope=0.9597,cex=4) +
  theme(axis.text.x=element_text(size=55),
        axis.text.y=element_text(size=55),
        plot.title = element_text(size=75),) +
  geom_text(size=20,aes(x=250,y=25,hjust=0,vjust=0,label="0.96"))+ theme(
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
p8=ggplot(merged,aes(means_pr_trial,unfactor(`Gloria`)))+geom_point(col="#00ba38",cex=10,alpha=0.8) + ylim(c(0, 300)) + xlim(c(0,300)) + ggtitle("Gloria") +geom_abline(intercept = 0, slope = 1,col="black",cex=2) + geom_abline(col="#00ba38",intercept=-24.0104,slope=1.1923,cex=4) +
  theme(axis.text.x=element_text(size=55),
        axis.text.y=element_text(size=55),
        plot.title = element_text(size=75),) +
  geom_text(size=20,aes(x=250,y=25,hjust=0,vjust=0,label="1.19"))+ theme(
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
p9=ggplot(merged,aes(means_pr_trial,unfactor(`Gracia`)))+geom_point(col="#00ba38",cex=10,alpha=0.8) + ylim(c(0, 300)) + xlim(c(0,300)) + ggtitle("Gracia") +geom_abline(intercept = 0, slope = 1,col="black",cex=2) + geom_abline(col="#00ba38",intercept=-7.24819,slope=1.11538,cex=4) +
  theme(axis.text.x=element_text(size=55),
        axis.text.y=element_text(size=55),
        plot.title = element_text(size=75),) +
  geom_text(size=20,aes(x=250,y=25,hjust=0,vjust=0,label="1.12"))+ theme(
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
p10=ggplot(merged,aes(means_pr_trial,unfactor(`Kontu`)))+geom_point(col="#00ba38",cex=10,alpha=0.8) + ylim(c(0, 300)) + xlim(c(0,300)) + ggtitle("Kontu") +geom_abline(intercept = 0, slope = 1,col="black",cex=2) + geom_abline(col="#00ba38",intercept=-16.4106,slope=0.8401,cex=4) +
  theme(axis.text.x=element_text(size=55),
        axis.text.y=element_text(size=55),
        plot.title = element_text(size=75),) +
  geom_text(size=20,aes(x=250,y=25,hjust=0,vjust=0,label="0.84"))+ theme(
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
p11=ggplot(merged,aes(means_pr_trial,unfactor(`Lynx`)))+geom_point(col="#00ba38",cex=10,alpha=0.8) + ylim(c(0, 300)) + xlim(c(0,300)) + ggtitle("Lynx") +geom_abline(intercept = 0, slope = 1,col="black",cex=2) + geom_abline(col="#00ba38",intercept=11.3656 ,slope=1.0494,cex=4) +
  theme(axis.text.x=element_text(size=55),
        axis.text.y=element_text(size=55),
        plot.title = element_text(size=75),) +
  geom_text(size=20,aes(x=250,y=25,hjust=0,vjust=0,label="1.05"))+ theme(
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
p12=ggplot(merged,aes(means_pr_trial,unfactor(`Mistral`)))+geom_point(col="#00ba38",cex=10,alpha=0.8) + ylim(c(0, 300)) + xlim(c(0,300)) + ggtitle("Mistral") +geom_abline(intercept = 0, slope = 1,col="black",cex=2) + geom_abline(col="#00ba38",intercept=-0.1765,slope=1.0044,cex=4) +
  theme(axis.text.x=element_text(size=55),
        axis.text.y=element_text(size=55),
        plot.title = element_text(size=75),) +
  geom_text(size=20,aes(x=250,y=25,hjust=0,vjust=0,label="1.00"))+ theme(
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
p13=ggplot(merged,aes(means_pr_trial,unfactor(`Pyramid`)))+geom_point(col="#00ba38",cex=10,alpha=0.8) + ylim(c(0, 300)) + xlim(c(0,300)) + ggtitle("Pyramid") +geom_abline(intercept = 0, slope = 1,col="black",cex=2) + geom_abline(col="#00ba38",intercept=20.9683,slope=0.9448,cex=4) +
  theme(axis.text.x=element_text(size=55),
        axis.text.y=element_text(size=55),
        plot.title = element_text(size=75),) +
  geom_text(size=20,aes(x=250,y=25,hjust=0,vjust=0,label="0.94"))+ theme(
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
p14=ggplot(merged,aes(means_pr_trial,unfactor(`Snowdrop`)))+geom_point(col="#00ba38",cex=10,alpha=0.8) + ylim(c(0, 300)) + xlim(c(0,300)) + ggtitle("Snowdrop") +geom_abline(intercept = 0, slope = 1,col="black",cex=2) + geom_abline(col="#00ba38",intercept=69.7920,slope=0.3315,cex=4) +
  theme(axis.text.x=element_text(size=55),
        axis.text.y=element_text(size=55),
        plot.title = element_text(size=75),) +
  geom_text(size=20,aes(x=250,y=25,hjust=0,vjust=0,label="0.33"))+ theme(
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
p15=ggplot(merged,aes(means_pr_trial,unfactor(`SSNS-1`)))+geom_point(col="#00ba38",cex=10,alpha=0.8) + ylim(c(0, 300)) + xlim(c(0,300)) + ggtitle("SSNS-1") +geom_abline(intercept = 0, slope = 1,col="black",cex=2) + geom_abline(col="#00ba38",intercept=-1.9893,slope=0.9406,cex=4) +
  theme(axis.text.x=element_text(size=55),
        axis.text.y=element_text(size=55),
        plot.title = element_text(size=75),) +
  geom_text(size=20,aes(x=250,y=25,hjust=0,vjust=0,label="0.94"))+ theme(
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
p16=ggplot(merged,aes(means_pr_trial,unfactor(`Taifun`)))+geom_point(col="#00ba38",cex=10,alpha=0.8) + ylim(c(0, 300)) + xlim(c(0,300)) + ggtitle("Taifun") +geom_abline(intercept = 0, slope = 1,col="black",cex=2) + geom_abline(col="#00ba38",intercept=12.34668,slope=0.93737,cex=4) +
  theme(axis.text.x=element_text(size=55),
        axis.text.y=element_text(size=55),
        plot.title = element_text(size=75),) +
  geom_text(size=20,aes(x=250,y=25,hjust=0,vjust=0,label="0.94"))+ theme(
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
p17=ggplot(merged,aes(means_pr_trial,unfactor(`Vertigo`)))+geom_point(col="#00ba38",cex=10,alpha=0.8) + ylim(c(0, 300)) + xlim(c(0,300)) + ggtitle("Vertigo") +geom_abline(intercept = 0, slope = 1,col="black",cex=2) + geom_abline(col="#00ba38",intercept=-24.9007,slope=1.2725,cex=4) +
  theme(axis.text.x=element_text(size=55),
        axis.text.y=element_text(size=55),
        plot.title = element_text(size=75),) +
  geom_text(size=20,aes(x=250,y=25,hjust=0,vjust=0,label="1.27"))+ theme(
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

correlation(Filtered_protein$Protein,Filtered_protein$Yield..g..Pr.kvm.)


CM<- lmer(Protein ~ (1|Cultivar)+(1|Location)+(1|Year)+Average_neighbour_score+(1|Cultivar:Location)+(1|Cultivar:Year), data = Filtered_protein,REML=T) 
summary(CM)

Residuals <- residuals(CM)
shapiro.test(Residuals) 
qqnorm(Residuals); qqline(Residuals)
plot(CM, abline=c(0, 0))


VarCorr(CM)


CM1<- lmer(Yield..g..Pr.kvm. ~ (1|Cultivar)+(1|Location)+(1|Year)+Average_neighbour_score+(1|Cultivar:Location)+(1|Cultivar:Year), data = Filtered_protein,REML=T) 
summary(CM1)

Residuals <- residuals(CM1)
shapiro.test(Residuals) 
qqnorm(Residuals); qqline(Residuals)
plot(CM1, abline=c(0, 0))

#Check how much the genetic correlation is by taking pearsson correlation of BLUPs in the different models.
library("nlme")
Cultivar_effects_protein=ranef(CM)$Cultivar 
Cultivar_effects_yield=ranef(CM1)$Cultivar 

correlation(Cultivar_effects_protein,Cultivar_effects_yield) #Only BLUP of yield on the data we have protein data on 

correlation(Cultivar_effects_protein,Cultivar_effects) #BLUPs of yield on all data

plot(Cultivar_effects_protein[,1]~Cultivar_effects_yield[,1])
new=cbind(rownames(Cultivar_effects_protein),as.numeric(as.character(Cultivar_effects_protein[,1])),as.numeric(as.character(Cultivar_effects_yield[,1])))
colnames(new)=c("cultivar","BLUPS_Protein","BLUPS_Yield")
new=as.data.frame(new)




#Produce Figure 3
library("varhandle")
ggplot(new, aes(x=unfactor(BLUPS_Yield), y=unfactor(BLUPS_Protein)))+geom_point(col="#00ba38",cex=3.5) +
  geom_text(label=new$cultivar,size=5.5,vjust = 0, nudge_y = 0.1) +
  stat_smooth(method="lm",col="#00ba38") +
  labs(x="Cultivar BLUPs Yield", y = "Cultivar BLUPs Protein") + theme_bw()
  




## Make Figure 1B
CVpercent
total$CV=c(41.9,37.1,41.0,43.2,44.0,42.7,38.7,47.5,43.9,42.8,41.2,40.9,41.1,36.1,41.0,40.4,45.1)
plot=ggplot(total, aes(x=MeanYield, y=CV))+geom_point(col="#00ba38",cex=3.5,alpha=0.8) +
  geom_text(label=total$Cultivar,size=4.5, vjust=0.5) +
  #scale_y_continuous(labels=scaleFUN) +
  #scale_x_continuous(labels=scaleFUN) +
  ylab("Coefficient of variance (%)") +
  xlab("Mean yield g/sqm") +
  geom_vline(xintercept=mean(filtered$Yield..g..Pr.kvm.),colour="black",linetype=2) +
  geom_hline(yintercept=mean(total$CV),colour="black",linetype=2)

plot+theme_bw()


### Calculate correlations between stability parameters
total$regression=c(0.97,0.91,0.97,1.04,1.12,1.15,1.07,1.04,1.12,0.70,1.13,0.92,1.11,0.60,0.90,1.01,1.19)

correlation(total$CV,total$MeanYield)
correlation(total$YieldVar,total$MeanYield)
correlation(total$regression,total$MeanYield)

correlation(total$CV,total$YieldVar)
correlation(total$CV,total$regression)
correlation(total$YieldVar,total$regression)


