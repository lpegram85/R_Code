###############################################################
#### Meta-Analysis by  Lisa.Pegram                  
###############################################################

### LIBRARIES ####

# 
# INSTALL PACKAGES (IF NOT ALREADY IN THE LIBRARY)

install.packages("mvtnorm") 
install.packages("compute.es") 
install.packages("MAd") 
install.packages("metafor") 
install.packages("ggplot2") 
install.packages("gridExtra") 
install.packages("Gmisc") 

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

id<- c(2,3.2,4,5,6,7,9) #7 checked
t<-c(1.49,0.43,-0.6,NA,-1.94,1.13,NA)#7 checked
df<-c(26,24,36,31,27,38,22)#7 checked
n<-c(28,26,38,33,29,40,24) #7 checked
p<-c(0.148,0.67,0.55, 0.16,0.06,0.27,0.23) #7 checked
s_mean<-c(8.31,8.13,7,7.65,15.54,15.23,16.74)#7 checked
h_mean<-c(8.67,7.89,6.6,8.38,17.11,14,15.97) #7 checked
s_sd<-c(1.01,1.52,1.97,0.36,2.13,3.2,1.4)#7 checked
h_sd<-c(.78,1.19,2.14,0.35,2.2,3.65,1.62)#7 checked
s_n<-c(16,12,18,17,15,23,13)#7 checked
h_n<-c(12,14,20,16,14,23,11)#7 checked

#learning<-c(“Encoding”)
#university<-c(“CalStateLA”)


############################################################
###########################DECODING#########################
############################################################
out <- data.frame(id,t,df,n,p,s_mean,h_mean,s_sd,h_sd,s_n,h_n)

# OUTCOME 2 EFFECT SIZES #null replace
res2B <- mes(m.1=out[1,6], m.2=out[1,7], sd.1=out[1,8], sd.2=out[1,9],n.1= out[1,10], n.2= out[1,11], level=95, dig=2, id=out[1,1], data=NULL)[,c(1,2,3,12,13,18)]
names(res2B) <- c("id", "s_n", "h_n", "g", "var.g", "pval.g")
res2B$outcome <- "g2B" 

#OUTCOME 3.2 EFFECT SIZES

res3an2 <- mes(m.1=out[2,6], m.2=out[2,7], sd.1=out[2,8], sd.2=out[2,9],n.1= out[2,10], n.2= out[2,11], level=95, dig=2, id=out[2,1], data=NULL)[,c(1,2,3,12,13,18)]
names(res3an2) <- c("id", "s_n", "h_n", "g", "var.g", "pval.g")
res3an2$outcome <- "g3.2" 

# OUTCOME 4 EFFECT SIZES
res4B <- mes(m.1=out[3,6], m.2=out[3,7], sd.1=out[3,8], sd.2=out[3,9],n.1= out[3,10], n.2= out[3,11], level=95, dig=2, id=out[3,1], data=NULL)[,c(1,2,3,12,13,18)]
names(res4B) <- c("id", "s_n", "h_n", "g", "var.g", "pval.g")
res4B$outcome <- "g4B"

# OUTCOME 5 EFFECT SIZES

res5B <- mes(m.1=out[4,6], m.2=out[4,7], sd.1=out[4,8], sd.2=out[4,9],n.1= out[4,10], n.2= out[4,11], level=95, dig=2, id=out[4,1], data=NULL)[,c(1,2,3,12,13,18)]
names(res5B) <- c("id", "s_n", "h_n", "g", "var.g", "pval.g")
res5B$outcome <- "g5"

# OUTCOME 6 EFFECT SIZES v1

#res6 <- mes(m.1=out[5,6], m.2=out[5,7], sd.1=out[5,8], sd.2=out[5,9],n.1= out[5,10], n.2= out[5,1], level=95, dig=2, id=out[5,1], data=NULL)[,c(1,2,3,12,13,18)]
#names(res6) <- c("id", "s_n", "h_n", "g", "var.g", "pval.g")
#res6$outcome <- "g6"

# OUTCOME 6 EFFECT SIZES wE NEED TO AGGREGATE USED  FOR VALUE
#d to g function NEED THE n.1 and n.2 assumed 15 and 14
res6<- des(d=.73, n.1=15, n.2=14, level = 95, dig = 2, verbose = TRUE, id=NULL, data=NULL)[,c(1,2,3,12,13,18)]
names(res6) <- c("id", "s_n", "h_n", "g", "var.g", "pval.g")
res6$outcome <- "g6"

# OUTCOME 7 EFFECT SIZES wE NEED TO AGGREGATE USED ANCOVA FOR VALUE
#a.fes(f, n.1, n.2, R, q, level=95, cer = 0.2, dig = 2, verbose = TRUE, id=NULL, data=NULL) NEED THE n.1 and n.2 assumed 23 and 23
res7<- des(d=.36, n.1=23, n.2=23, level = 95,  dig = 2, verbose = TRUE, id=NULL, data=NULL)[,c(1,2,3,12,13,18)]
names(res7) <- c("id", "s_n", "h_n", "g", "var.g", "pval.g")
res7$outcome <- "g7"

# OUTCOME 8 EFFECT SIZES wE NEED TO AGGREGATE USED ANCOVA FOR VALUE
#a.fes(f, n.1, n.2, R, q, level=95, cer = 0.2, dig = 2, verbose = TRUE, id=NULL, data=NULL) NEED THE n.1 and n.2 assumed 23 and 23
res8<- des(d=.44, n.1=23, n.2=23, level = 95,  dig = 2, verbose = TRUE, id=NULL, data=NULL)[,c(1,2,3,12,13,18)]
names(res8) <- c("id", "s_n", "h_n", "g", "var.g", "pval.g")
res8$outcome <- "g8"

# OUTCOME 5 EFFECT SIZES

res9 <- mes(m.1=out[7,6], m.2=out[7,7], sd.1=out[7,8], sd.2=out[7,9],n.1= out[7,10], n.2= out[7,11], level=95, dig=2, id=out[7,1], data=NULL)[,c(1,2,3,12,13,18)]
names(res9) <- c("id", "s_n", "h_n", "g", "var.g", "pval.g")
res9$outcome <- "g9"

dat.sim.es <-rbind(res2B, res3an2, res4B, res5B, res6, res7, res8,res9) # BIND EACH DATA SET TOGETHER 
sapply(dat.sim.es,class)

#OMNIBUS MODEL FOR RANDOM EFFECTS LIKELIHOOD MODEL
m0 <- mareg(g~1, var=var.g, data=dat.sim.es)
summary(m0)
#A q-value of 5% means that 5% of significant results will result in false positives.
#small efect Small effect (cannot be discerned by the naked eye) <= 0.2 and our value is .061


#' MODEL COMPARISONS
# likelihood ratio test 
# alternative to Wald test, this conducts a full vs reduced model
# comparison.  Need to add the other 
m0a <- mareg(g~1, var=var.g, data=dat.sim.es, method = "ML") # EMPTY OMNIBUS MODEL
#m3a <- mareg(es~dose+stress, var=var, data=dat.sim.final, method = "ML")  # FULL MODEL
summary(m0a)

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

# FOR THE TUTORIAL GRAPHIC, CHANGE ID VARIABLES FOR BETTER DISPLAY 
dat.sim.final$id <- paste("Study", 1:length(dat.sim.final$id))

forest(m0)
text(-2.4 , 9.7, "Study",   pos=4,  cex=.7)
text(3.81,  9.7, "Observed g [95% CI]", pos=2, cex=.7)
confint(m0,digits=2)

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

