###############################################################
#### Meta-Analysis by  Lisa.Pegram                  
###############################################################

### LIBRARIES ####

# 

# LOAD THE PACKAGES INTO THE CURRENT R SESSION
library(mvtnorm)        
# TO SIMULATE MULTIVARIATE DATA
library(compute.es)     
# TO COMPUTE EFFECT SIZES
library(MAd)           
# META-ANALYSIS PACKAGE
library(metafor)        
# META-ANALYSIS PACKAGE
library(ggplot2)        
# GRAPHICS PACKAGE
library(gridExtra)      
# ARRANGE GRID FOR GRAPHICS
library(Gmisc)          
# TABLE GENERATION
library(scales)         
# ADDITIONAL GRAPHICS


id<- c(1,2.1,2.1,4.1,5.1,10) #6 checked
t<-c(1.12,-1.01,-1.00,-1.45,-.36,NA)#6 checked
df<-c(19,26,26,36,31,19)#6 checked
n<-c(21) #6 checked
p<-c(0.28,.31,.32,0.156,.72,0.005) #6 checked
s_mean<-c(7.533,8.31,8.31,8.778,7.1324,12.72)#6 checked
h_mean<-c(8.833,8.43,8.67,8.075,7.4063,16.45) #6 checked
s_sd<-c(3.72,1.01,1.01,1.3956,1.94478,3.31)#6 checked
h_sd<-c(1.60,.85,.78,1.5751,2.36621,1.69)#6 checked
s_n<-c(15,16,16,20,17,11)#6 checked
h_n<-c(6,14,12,18,16,10)#6 checked



############################################################
###########################ENCODING #########################
############################################################



outE <- data.frame(id,t,df,n,p,s_mean,h_mean,s_sd,h_sd,s_n,h_n)

# OUTCOME 1 EFFECT SIZES #null replace
res1 <- mes(m.1=outE[1,6], m.2=outE[1,7], sd.1=outE[1,8], sd.2=outE[1,9],n.1= outE[1,10], n.2= outE[1,11], level=95, dig=2, id=outE[1,1], data=NULL)[,c(1,2,3,12,13,18)]
names(res1) <- c("id", "s_n", "h_n", "g", "var.g", "pval.g")
res1$outcome <- "g1" 

res1final=c(res1[,"id"],res1[,"g"],res1[,"var.g"])
names(res1final) <- c("id", "g", "var.g")
# OUTCOME 2.1 EFFECT SIZES #Aggregated

res2A <- mes(m.1=outE[2,6], m.2=outE[2,7], sd.1=outE[2,8], sd.2=outE[2,9],n.1= outE[2,10], n.2= outE[2,11], level=95, dig=2, id=outE[2,1], data=NULL)[,c(2,3,12,13,18)]
names(res2A) <- c("s_n", "h_n", "g", "var.g", "pval.g")
res2A$outcome <- "g2A" 

res2A.2 <- mes(m.1=outE[3,6], m.2=outE[3,7], sd.1=outE[3,8], sd.2=outE[3,9],n.1= outE[3,10], n.2= outE[3,11], level=95, dig=2, id=outE[3,1], data=NULL)[,c(2,3,12,13,18)]
names(res2A.2) <- c("s_n", "h_n", "g", "var.g", "pval.g")
res2A.2$outcome <- "g2A2" 

dat.sim.2ONLY <-rbind(res2A, res2A.2)

RES2agg<- agg(
  id = 2, es = g, var = var.g,
  n.1 = s_n, n.2 = h_n, cor = 0.5,
  method = "BHHR", data = dat.sim.2ONLY
)  #s_n, h_n pval.g outcome 

names(RES2agg) <- c("id", "g", "var.g")
RES2agg$id <- "2" 


# OUTCOME 4.1 EFFECT SIZES
res4A <- mes(m.1=outE[4,6], m.2=outE[4,7], sd.1=outE[4,8], sd.2=outE[4,9],n.1= outE[4,10], n.2= outE[4,11], level=95, dig=2, id=outE[4,1], data=NULL)[,c(1,2,3,12,13,18)]
names(res4A) <- c("id", "s_n", "h_n", "g", "var.g", "pval.g")
res4A$outcome <- "g4A"

res4Afinal=c(res4A[,"id"],res4A[,"g"],res4A[,"var.g"])
names(res4Afinal) <- c("id", "g", "var.g")

# OUTCOME 5 EFFECT SIZES

res5A <- mes(m.1=outE[5,6], m.2=outE[5,7], sd.1=outE[5,8], sd.2=outE[5,9],n.1= outE[5,10], n.2= outE[5,11], level=95, dig=2, id=outE[5,1], data=NULL)[,c(1,2,3,12,13,18)]
names(res5A) <- c("id", "s_n", "h_n", "g", "var.g", "pval.g")
res5A$outcome <- "g5A"

res5Afinal=c(res5A[,"id"],res5A[,"g"],res5A[,"var.g"])
names(res5Afinal) <- c("id", "g", "var.g")


# OUTCOME 10 EFFECT SIZES wE NEED TO AGGREGATE USED ANCOVA FOR VALUE
res10 <- mes(m.1=outE[6,6], m.2=outE[6,7], sd.1=outE[6,8], sd.2=outE[6,9],n.1= outE[6,10], n.2= outE[6,11], level=95, dig=2, id=outE[6,1], data=NULL)[,c(1,2,3,12,13,18)]
names(res10) <- c("id", "s_n", "h_n", "g", "var.g", "pval.g")
res10$outcome <- "g10"

res10final=c(res10[,"id"],res10[,"g"],res10[,"var.g"])
names(res10final) <- c("id", "g", "var.g")

dat.sim.final <-rbind(res1final, RES2agg, res4Afinal, res5Afinal, res10final)

sapply(dat.sim.final,class)


#OMNIBUS MODEL FOR RANDOM EFFECTS LIKELIHOOD MODEL
m0 <- mareg(g~1, var=var.g, method="REML", data=dat.sim.final)
summary(m0)
#A q-value of 5% means that 5% of significant results will result in false positives.
#small efect Small effect (cannot be discerned by the naked eye) <= 0.2 and our value is .184


#' MODEL COMPARISONS
# likelihood ratio test 
# alternative to Wald test, this conducts a full vs reduced model
# comparison.  Need to add the other 
m0a <- mareg(g~1, var=var.g, data=dat.sim.es, method = "ML") # EMPTY OMNIBUS MODEL
#m3a <- mareg(es~dose+stress, var=var, data=dat.sim.final, method = "ML")  # FULL MODEL
summary(m0a)
# NOTE: IN THIS SIMULATED DATA EXTRACTING A RANDOM SELECTION 
# OF K=8 (FROM A K=1000 POPULATION), EACH MODERATOR ACCOUNTS FOR
# ALL OF THE BETWEEN-STUDY VARIABILITY (REDUCING TAU^2 TO 0).
# THIS IS UNLIKELY TO HAPPEN WHEN EXAMINING REAL DATA WITH LARGER K.
#anova(m0a, m3a) 

####' DIAGNOSTICS
funnel(m0) # EXAMINE PUBLICATION BIAS VISUALLY 
fsn(yi=m0$yi, vi=m0$vi)  # NO PUB BIAS no publication effect found that out by symmetry estimator of unpublished studies in meta-analysis
#The Rosenthal method (sometimes called a ‘file drawer analysis’) calculates the number of studies
#averaging null results that would have to be added to the given set of observed outcomes to reduce
#the combined significance level (p-value) to a target alpha level (e.g., .05). The calculation is based
#on Stouffer’s method to combine p-values and is described in Rosenthal (1979).
#how many missing studies we would need to retrieve and incorporate in the, doesn't have a significant effect so we don't need it
#analysis before the p-value became nonsignificant.

# FIGURES FOR SIMULATED DATA, Plotted well 
fig1a <- ggplot(res2B, aes(x=g)) + 
  geom_histogram(binwidth=.1, colour="black", fill="blue") +
  geom_vline(aes(xintercept=median(g, na.rm=T)),   
             color="red", linetype="dashed", size=1)+
  ggtitle("Outcome One in Study Population")+
  xlab("g")+
  theme_bw()  

fig1b <- ggplot(res3an2, aes(x=g)) + 
  geom_histogram(binwidth=.1, colour="black", fill="yellow") +
  geom_vline(aes(xintercept=median(g, na.rm=T)),   
             color="red", linetype="dashed", size=1)+
  ggtitle("Outcome Two in Study Population")+
  xlab("g")+
  theme_bw()

#grid.arrange(fig1a, fig1b, ncol=2, top="Effect Sizes Generated from Simulation")
#library(ggplot2)
#p1 = qplot(1:10,rnorm(10))
#p2 = qplot(1:10,rnorm(10))

#library(gridExtra)
#grid.arrange(p1, p2, ncol=2,top="Main Title")
#grid.arrange( grobs = list(p1,p2,...),...

# VISUALIZE THE HETEROGENEITY (fig2)
confint(m0,digits=2)

# FOR THE TUTORIAL GRAPHIC, CHANGE ID VARIABLES FOR BETTER DISPLAY 
dat.sim.final$id <- paste("Study", 1:length(dat.sim.final$id))

forest(m0)
text(-2.4 , 9.7, "Study",   pos=4,  cex=.7)
text(3.81,  9.7, "Observed g [95% CI]", pos=2, cex=.7)

# DOSE MOD (CONTINUOUS)

# USING FUNCTION IN MAd PACKAGE (HERE IT DOESNT MATCH UP WELL SO WE
# WILL USE THE CODE BELOW plotcon CODE)
plotcon(g=es, var=es, mod=dose, data=dat.sim.final, modname = "Treatment dose", 
        title = "Scatterplot examining dose moderator")

# LOOKS BETTER:
fig3a <- ggplot(data=dat.sim.final,aes(x=dose, y=es, size=1/var), na.rm=TRUE) +
  geom_point(alpha=.5) +   
  geom_abline(intercept=m1$b[1], slope=m1$b[2]) + 
  #ylim(0,1.25)+
  expand_limits(x = 5, y=0) +
  xlab("Treatment dose") + ylab("Effect Size") +  
  ggtitle("(a) Scatterplot examining dose moderator") +
  theme_bw() +
  theme(legend.position = "none")+
  scale_x_continuous(breaks= pretty_breaks())

# STRESS MOD (CATEGORICAL)
dat.sim.final$stress <- relevel(dat.sim.final$stress, ref="low")

# USING FUNCTION IN MAd PACKAGE
plotcat(g=es, var=es, mod=stress, data=dat.sim.final, modname = "Stress",
        title = "Boxplot examining stress moderator") 

# FOR ADDED FLEXIBILITY USE:
fig3b <- ggplot(dat.sim.final,  aes(factor(stress), es, weight = 1/var), na.rm=TRUE) + 
  geom_boxplot(aes(weight = 1/var), outlier.size=3,na.rm=TRUE) + 
  geom_jitter(aes(shape=factor(stress), color=factor(stress), size=1/var), alpha=.5) + 
  xlab("Stress") + 
  ylab("Effect Size")  +
  ggtitle("(b) Boxplot examining stress moderator") +
  theme_bw() +
  theme(legend.position = "none")


fig3c <- ggplot(dat.sim.final,  aes(factor(stress), dose), na.rm=TRUE) + 
  geom_boxplot(outlier.size=3,na.rm=TRUE) + 
  geom_jitter(aes(shape=factor(stress), color=factor(stress), size=1/var), alpha=.5) + 
  coord_flip()+
  xlab("Stress") + 
  ylab("Treatmen dose")  +
  ggtitle("(c) Confounding between moderators") +
  theme_bw() +
  theme(legend.position = "none")

fig3d <- ggplot(data=dat.sim.final,aes(x=dose, y=es, size=1/var), na.rm=TRUE) +
  geom_point(alpha=.5, aes(shape=stress, color=factor(stress))) +   
  geom_abline(intercept=m1$b[1], slope=m1$b[2]) + 
  #ylim(0,1.25)+
  expand_limits(x = 5, y=0) +
  xlab("Treatment dose") + ylab("Effect Size") +  
  ggtitle("(d) Effect sizes by dose and stress moderators") +
  theme_bw() +
  theme(legend.position = "none")+
  scale_x_continuous(breaks= pretty_breaks())

grid.arrange(fig3a, fig3b, fig3c, fig3d, ncol=2, main="Visualizing continuous and categorical moderator variables")
