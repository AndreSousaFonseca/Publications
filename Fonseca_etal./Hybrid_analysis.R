# Load libraries:
#Load libraries
library(MASS)       # Perform general functions
library(dplyr)      # Better handle data
library(pROC)       # Perofrm Roc curve 
#library(nortest)    # Perform normality test
library(mixsmsn)    # Mixture model
library(AID)        # Perform box-cox transformation
library(doParallel) # Parallel processing
library(sn)         # Regression skew
library(tidyr)      # Handle data for ggplot2
library(ggplot2)    # Plot the data
library(lmtest)     # Perform the likelihood test
library(caret)      # Use for CV prerformance
library(SuperLearner) # Apply Superlearner to make predictions
library(OptimalCutpoints) # Obtain the optimal cutpoints
#library(naniar)     # View missing data
library(ggrepel)
set.seed(42)       # Have reproducible results

# Specify the number of cores
numCores <- detectCores()
registerDoParallel(numCores)  # use multicore, set to the number of our cores

#Load and pre-process the data
KEN <- read.csv("KEN.csv", header=T, check.names=FALSE)

# Attribute subject id to rownames
rownames(KEN) <- KEN$SubjectID

# select only patients with parasite at screening
my.data <- subset(KEN, KEN$ParasiteAtScreen==1)

# Ignore variables that were not used in the original paper
ignore <- c("SubjectID","NEpisodes","TimeToFirstEpisode", "ParasiteAtScreen", "Age", "SchizontReactivity", "Status")

# Specify the predictive variables
X <- my.data[, !colnames(my.data) %in% ignore] 

# Specify the outcome variable
Y <- data.frame( "Status" = my.data[, colnames(my.data) %in% "Status"]); rownames(Y) <- my.data$SubjectID

#  Specify parameters for the box-cox function
n <- ncol(X)
outcome <-  as.factor(Y$Status)
out <- ifelse(outcome == "protected",1,0)

# Run the algorithm across all aantibodies
#bc.vals <- foreach (i=X[,1:n], .combine=rbind) %dopar% {
#  bc <- AID::boxcoxlm(x= matrix(out, ncol = 1, nrow = length(out)), y= i, lambda = seq(-10,10,0.01))
#  vec <-  data.frame("lambda" = bc$lambda.hat, "p.value" = bc$p.value)
#}
#Save data
#saveRDS(bc.vals,"suplementary1")
bc.vals <- readRDS("suplementary1")

# Attribute rownames
rownames(bc.vals) <- colnames(X)

# Add a column called "Sig" wich will allow us to collor the dots according to their significance level
bc.vals$sigificance <- as.factor(ifelse(bc.vals$lambda > -4 & bc.vals$p.value > 0.05, "p.value > 0.05", "p.value < 0.05"))
bc.vals$names <-  rownames(bc.vals)

ggplot(bc.vals,  aes(lambda,-log10(p.value), color = sigificance))+ 
  geom_point(size = 3) +
  geom_vline(xintercept = 4, linetype="dashed", color = "grey") + theme_classic() +
  geom_vline(xintercept = -4, linetype="dashed", color = "grey") + theme_classic() +
  scale_y_continuous(labels= c(1, 0.01, 0.0001, 0.000001)) +
  scale_x_continuous(position = 'top') +
  scale_color_manual(values=c("grey", "#3498db")) +
  geom_hline(yintercept= -log10(0.05), linetype="dashed", color = "grey") + theme_classic() +
  ylab("p.value") + 
  geom_text_repel(data = subset(bc.vals, p.value > 0.05),
    aes(label = names), size = 3,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines"))

# Plot
#ggplot(bc.vals, aes(-log10(p.value), lambda, color = sigificance))+ geom_point(size = 2) +
#  geom_hline(yintercept= 4, linetype="dashed", color = "grey") + theme_classic() +
#  geom_hline(yintercept= -4, linetype="dashed", color = "grey") + theme_classic() +
#  scale_color_manual(values=c("grey", "#F8766D"))  + scale_x_continuous(labels= c(0,0.01, 0.0001, 0.000001),position = 'top') +
#  geom_vline(xintercept= -log10(0.05), linetype="dashed", color = "grey") + theme_classic() +
#  xlab("p.value") 

## Select the normally distributed antibodies
bc.normal <- subset(bc.vals, bc.vals$p.value > 0.05)
print(bc.normal)
nrow(bc.normal)
min(bc.normal$lambda)
max(bc.normal$lambda)

## Obtain the transformed data for the 6 antibodies  
# Subset only the 6 antibodies from the original dataset
tf.data <- X[, colnames(X) %in% rownames(bc.normal)]

# Obtain the transformed datasets for the 6 antibodies
tf.datasets <- data.frame("Status" = Y$Status)

for (i in 1:length(tf.data)) {
  bc.val <- (tf.data[,i]^bc.normal$lambda[i] - 1)/ bc.normal$lambda[i]
  tf.datasets <- cbind(tf.datasets, bc.val)
}

# Attribute the names
colnames(tf.datasets)[2:ncol(tf.datasets)] <-  colnames(tf.data)


                        # Perform the T-test
# Remove the Status variable
t.data <- tf.datasets[, !colnames(tf.datasets) %in% "Status"]

# Writte the T-test function
t_test <- function(x) {
  ab.t <- t.test(x ~ Y$Status, alternative = "two.sided", var.equal = F )$p.value
  return(ab.t)
}

# Perform the t-test
t.vals <- data.frame(sapply(t.data, t_test))
colnames(t.vals) <- "p.value"


                                  # Analyze the non-normally distributed antibodies throught Mixture Models
# Select the non normally distributed antibodies
non.data <- X[, ! colnames(X) %in% colnames(t.data)]
mm.ab <- log10(non.data)
print(dim(mm.ab))

## Prepare for the Mixture Model
#my.p.mix.skew.t<-function(x,pii,mu,sigma2,shape,nu){
#  
#  param<-cbind(pii,mu,sigma2,shape,nu)
#  prob<-sum(apply(param,1,function(y,param)param[1]*sn::pst(x,param[2],sqrt(param[3]),param[4],param[5]),y=x))
#  return(prob)
#  
#}
#
## Prepare the pearson function
#my.gof.smsn<-function(data,pii,mu,sigma2,shape=0,nu=Inf,lag=0.1,test='pearson'){
#  num.comp<-length(pii)
#  num.param<-3*num.comp-1
#  if(sum(shape!=0)!=0){
#    num.param<-num.param+num.comp
#  }
#  
#  if(nu!=Inf){
#    num.param<-num.param+1
#  }
#  
#  q1<-unique(quantile(data,seq(lag,1-lag,lag)))
#  quantiles1<-c(-Inf,q1,Inf)
#  obs<-table(cut(data,breaks=quantiles1))
#  num.bin<-length(obs)
#  p<-sapply(q1,my.p.mix.skew.t,pii=pii,mu=mu,sigma2=sigma2,shape=shape,nu=nu)
#  prob<-c(p[1],sapply(2:length(q1),function(x,p)p[x]-p[x-1],p=p))
#  prob<-c(prob,1-sum(prob))
#  n<-sum(obs)
#  exp<-n*prob
#  print(rbind(obs,exp))
#  
#  if(test=='pearson'){
#    cat("Pearson's goodness-of-fit test\n")
#    stat<-sum(((obs-exp)^2)/exp)		
#  } else {
#    cat("LRT goodness-of-fit test\n")
#    stat<-(-2)*sum(obs*(log(prob)-log(obs/n)))
#  }
#  df<-num.bin-num.param-1
#  p.value<-1-pchisq(stat,df)
#  cat("Statistic=",stat," (",df," d.f.)\n",sep='')
#  cat("P-value=",round(p.value,3),"\n",sep='')
#  return(p.value)
#}
#
#
## Run the mixture model
#l=ncol(mm.ab)
#output <- c()
#
#set.seed(42)
#data.out <- foreach(antibody = mm.ab[1:l], .combine=rbind) %dopar%{
#  output <- c()
#  #---- NORMAL DISTRIBUTION ----#
#  
#  mm.n <- function(ab) {return(tryCatch(mixsmsn::smsn.mix(ab, nu =3, g = 1, family = "Normal", iter.max = 1000), error=function(e) NULL))}
#  m.n <- replicate(10, mm.n(antibody), simplify = FALSE)
#  n.l <- cbind(lapply(m.n, function(y)y$aic))
#  p.n <- data.frame("AIC" =n.l)
#  p.n$iter <- 1:nrow(p.n) 
#  p.n$AIC <- as.numeric(as.character(p.n$AIC))
#  p.n <-  p.n[which.min(p.n$AIC),]
#  
#  
#  if (nrow(p.n) != 0) {
#    
#    best.fit.norm1 <- m.n[[p.n$iter]]
#    
#    p1.pearson<-my.gof.smsn(antibody,
#                            pii=best.fit.norm1$pii,mu=best.fit.norm1$mu,
#                            sigma2=best.fit.norm1$sigma2,shape=0,lag=0.1,test='pearson')
#    
#    output<- rbind(output, c('Normal',1,round(best.fit.norm1$bic,2),
#                             round(best.fit.norm1$aic,2),round(p1.pearson,3)))
#  }else{
#    output<- rbind(output, c('Normal', 1 ,999,999, 0))
#  }
#  n.l <- c()
#  
#  # Mixture of 2 Normal Distribution #
#  
#  mm.n2 <- function(ab) {return(tryCatch(mixsmsn::smsn.mix(ab, nu =3, g = 2, family = "Normal", iter.max = 1000), error=function(e) NULL))}
#  m.n2 <- replicate(10, mm.n2(antibody), simplify = FALSE)
#  n.l2 <- cbind(lapply(m.n2, function(y)y$aic))
#  p.n2 <- data.frame("AIC" =n.l2)
#  p.n2$iter <- 1:nrow(p.n2) 
#  p.n2$AIC <- as.numeric(as.character(p.n2$AIC))
#  p.n2 <-  p.n2[which.min(p.n2$AIC),]
#  
#  if (nrow(p.n2) != 0) {
#    
#    best.fit.norm2 <- m.n2[[p.n2$iter]]
#    
#    p1.pearson<-my.gof.smsn(antibody,
#                            pii=best.fit.norm2$pii,mu=best.fit.norm2$mu,
#                            sigma2=best.fit.norm2$sigma2,shape=c(0,0),lag=0.1,test='pearson')
#    
#    
#    output <- rbind(output, c("mixNormal",2,round(best.fit.norm2$bic,2), 
#                              round(best.fit.norm2$aic,2),round(p1.pearson,3)))
#  }else{
#    output<- rbind(output,c('mixNormal', 2 ,999,999, 0))
#  }
#  n.l2 <-c ()
#  
#  #--- normal skew distribution ---#
#  
#  mm.sn <- function(ab) {return(tryCatch(mixsmsn::smsn.mix(ab, nu =3, g = 1, family = "Skew.normal", iter.max = 1000), error=function(e) NULL))}
#  m.sn <- replicate(10, mm.sn(antibody), simplify = FALSE)
#  sn.l <- cbind(lapply(m.sn, function(y)y$aic))
#  p.sn <- data.frame("AIC" =sn.l)
#  p.sn$iter <- 1:nrow(p.sn) 
#  p.sn$AIC <- as.numeric(as.character(p.sn$AIC))
#  p.sn <-  p.sn[which.min(p.sn$AIC),]
#  
#  if (nrow(p.sn) != 0) {
#    
#    best.fit.sn1 <- m.sn[[p.sn$iter]]
#    
#    p1.pearson<-my.gof.smsn(antibody,
#                            pii=best.fit.sn1$pii,mu=best.fit.sn1$mu,
#                            sigma2=best.fit.sn1$sigma2,shape=best.fit.sn1$shape,lag=0.1,test='pearson')
#    
#    
#    output<-rbind(output,c('SkewNormal',1,round(best.fit.sn1$bic,2),
#                           round(best.fit.sn1$aic,2),round(p1.pearson,3)))
#    
#    
#  }else{
#    output<- rbnd(output,c('SkewNormal', 1 ,999,999,0))
#  }
#  sn.l <-c ()
#  
#  # Mixture of 2 Normal Skew Distribution #
#  
#  mm.sn2 <- function(ab) {return(tryCatch(mixsmsn::smsn.mix(ab, nu =3, g = 2, family = "Skew.normal", iter.max = 1000), error=function(e) NULL))}
#  m.sn2 <- replicate(10, mm.sn2(antibody), simplify = FALSE)
#  sn2.l <- cbind(lapply(m.sn2, function(y)y$aic))
#  p.sn2 <- data.frame("AIC" =sn2.l)
#  p.sn2$iter <- 1:nrow(p.sn2) 
#  p.sn2$AIC <- as.numeric(as.character(p.sn2$AIC))
#  p.sn2 <-  p.sn2[which.min(p.sn2$AIC),]
#  
#  if (nrow(p.sn2) != 0) {
#    
#    best.fit.sn2 <- m.sn2[[p.sn2$iter]]
#    
#    p1.pearson<-my.gof.smsn(antibody,
#                            pii=best.fit.sn2$pii,mu=best.fit.sn2$mu,
#                            sigma2=best.fit.sn2$sigma2,shape=best.fit.sn2$shape,lag=0.1,test='pearson')
#    
#    output<-rbind(output,c('mixSkewNormal',2,round(best.fit.sn2$bic,2),
#                           round(best.fit.sn2$aic,2),round(p1.pearson,3)))
#    
#    
#  }else{
#    output<- rbind(output,c('mixSkewNormal', 2 ,999,999,0))
#  }
#  sn2.l <-c ()
#  
#  
#  #--- T DISTRIBUTION ---#
#  
#  mm.t <- function(ab) {return(tryCatch(mixsmsn::smsn.mix(ab, nu =3, g = 1, family = "t", iter.max = 1000), error=function(e) NULL))}
#  m.t <- replicate(10, mm.t(antibody), simplify = FALSE)
#  t.l <- cbind(lapply(m.t, function(y)y$aic))
#  p.t <- data.frame("AIC" =t.l)
#  p.t$iter <- 1:nrow(p.t) 
#  p.t$AIC <- as.numeric(as.character(p.t$AIC))
#  p.t <-  p.t[which.min(p.t$AIC),]
#  
#  
#  if (nrow(p.t) != 0) {
#    
#    best.fit.t1 <- m.t[[ p.t$iter]]
#    
#    p1.pearson<-my.gof.smsn(antibody,
#                            pii=best.fit.t1$pii,mu=best.fit.t1$mu,
#                            sigma2=best.fit.t1$sigma2,shape=0,nu=best.fit.t1$nu,lag=0.1,test='pearson')
#    
#    
#    output<-rbind(output,c('T',1,round(best.fit.t1$bic,2),
#                           round(best.fit.t1$aic,2),round(p1.pearson,3)))
#    
#  }else{
#    output<- rbind(output, c('T', 1 ,999,999,0))
#  }
#  t.l <-c ()
#  # Mixture of 2 T Distribution #
#  
#  mm.t2 <- function(ab) {return(tryCatch(mixsmsn::smsn.mix(ab, nu =3, g = 2, family = "t", iter.max = 1000), error=function(e) NULL))}
#  m.t2 <- replicate(10, mm.t2(antibody), simplify = FALSE)
#  t2.l <- cbind(lapply(m.t2, function(y)y$aic))
#  p.t2 <- data.frame("AIC" =t2.l)
#  p.t2$iter <- 1:nrow(p.t2) 
#  p.t2$AIC <- as.numeric(as.character(p.t2$AIC))
#  p.t2 <-  p.t2[which.min(p.t2$AIC),]
#  
#  
#  if (nrow(p.t2) != 0) {
#    
#    best.fit.t2 <- m.t2[[p.t2$iter]]
#    
#    p1.pearson<-my.gof.smsn(antibody,
#                            pii=best.fit.t2$pii,mu=best.fit.t2$mu,
#                            sigma2=best.fit.t2$sigma2,shape=c(0,0),nu=best.fit.t2$nu,lag=0.1,test='pearson')
#    
#    
#    output<-rbind(output,c('mixT',2,round(best.fit.t2$bic,2),
#                           round(best.fit.t2$aic,2),round(p1.pearson,3)))
#    
#    
#  }else{
#    output<- rbind(output,c('mixT', 2 ,999,999,0))
#  }
#  
#  t2.l <-c ()
#  
#  #--- SKEW T DISTRIBUTION ---#
#  
#  mm.st <- function(ab) {return(tryCatch(mixsmsn::smsn.mix(ab, nu =3, g = 1, family = "Skew.t", iter.max = 1000),
#                                         error=function(e) NULL))}
#  m.st <- replicate(10, mm.st(antibody), simplify = FALSE)
#  st.l <- cbind(lapply(m.st, function(y)y$aic))
#  p.st <- data.frame("AIC" =st.l)
#  p.st$iter <- 1:nrow(p.st) 
#  p.st$AIC <- as.numeric(as.character(p.st$AIC))
#  p.st <-  p.st[which.min(p.st$AIC),]
#  
#  
#  if (nrow(p.st) != 0) {
#    
#    best.fit.st1 <- m.st[[ p.st$iter]]
#    
#    p1.pearson <-my.gof.smsn(antibody,
#                             pii=best.fit.st1$pii,mu=best.fit.st1$mu,
#                             sigma2=best.fit.st1$sigma2,shape=best.fit.st1$shape,
#                             nu=best.fit.st1$nu,lag=0.1,test='pearson')
#    
#    
#    output<-rbind(output,c('SkewT',1,round(best.fit.st1$bic,2),
#                           round(best.fit.st1$aic,2),round(p1.pearson,3)))
#    
#  }else{
#    output<- rbind(output, c('SkewT', 1 ,999,999,0))
#  }
#  st.l <-c ()
#  # Mixture of 2 Skew T Distribution #
#  
#  mm.st2 <- function(ab) {return(tryCatch(mixsmsn::smsn.mix(ab, nu =3, g = 2, family = "Skew.t", iter.max = 1000),
#                                          error=function(e) NULL))}
#  m.st2 <- replicate(10, mm.st2(antibody), simplify = FALSE)
#  st2.l <- cbind(lapply(m.st2, function(y)y$aic))
#  p.st2 <- data.frame("AIC" =st2.l)
#  p.st2$iter <- 1:nrow(p.st2) 
#  p.st2$AIC <- as.numeric(as.character(p.st2$AIC))
#  p.st2 <-  p.st2[which.min(p.st2$AIC),]
#  
#  
#  if (nrow(p.st2) != 0) {
#    
#    best.fit.st2 <- m.st2[[ p.st2$iter]]
#    
#    p1.pearson<-my.gof.smsn(antibody,
#                            pii=best.fit.st2$pii,mu=best.fit.st2$mu,
#                            sigma2=best.fit.st2$sigma2,shape=best.fit.st2$shape,nu=best.fit.st2$nu,lag=0.1,test='pearson')
#    
#    
#    output<-rbind(output,c('mixSkewT',2,round(best.fit.st2$bic,2),
#                           round(best.fit.st2$aic,2),round(p1.pearson,3)))
#    
#    
#  }else{
#    output<- rbind(output, c('mixSkewT', 2 ,999,999,0))
#  }
#  st2.l <-c ()
#  
#  output <-as.data.frame(output)
#  return(output)
#}
#
#
###Atribute names to data.out
#colnames(data.out)<-c('Model','Num.Comp','BIC','AIC','GOF.Pearson.P.10%')
#saveRDS(data.out, "paremetric.RDS")
data.out <-  readRDS("paremetric.RDS")

#Obtain the antibody names and multiply by 8 (the number of models performed)
ab <- colnames(mm.ab)
ab <- replicate(8, ab)
l <- nrow(ab)

ab.names <- c()
for (i in 1:l){
  ab.names <- c(ab.names,ab[i,])
}

data.out <- data.frame("Antibody"=ab.names, data.out)

# Transform character data to numeric data
for (i in 3:ncol(data.out)) {
  data.out[,i] <-  as.numeric(data.out[,i])
}

# Select the antibodies with the lowest aic while significant for the Pearson an LRT test
l.abs <- unique(data.out$Antibody)

best.auc <- function(x) {
  v.ab <-subset(data.out, data.out$Antibody == x)
  f10 <- subset(v.ab , v.ab$GOF.Pearson.P.10. >= 0.01) 
  
  if(nrow(f10) < 1){
    f.ab <-v.ab[1,]
    f.ab$AIC <- 999
  }
  else{
    f.ab <- subset(f10, f10$AIC == min(f10$AIC))
  }
  
  return(f.ab)
}

all.values <- data.frame()

for( i in l.abs){
  all.values <-rbind( all.values, best.auc(i))
}

# Obtain the non valid antibodies
f.non.all.values <- subset(all.values, all.values$AIC == 999 )
print(paste("A total of", dim(f.non.all.values)[1], "antibodies were NOT adjusted  by the mixture models"))

# Obtain the valid antibodies
f.all.values <- subset(all.values, all.values$AIC != 999 )
print(paste("A total of", dim(f.all.values)[1], "antibodies were well adjusted by the mixture models"))

# See the table of components for the valid ones
table(f.all.values$Num.Comp)
table(f.all.values$Model)


# Plot the data
ggplot(f.all.values, aes(AIC, GOF.Pearson.P.10.,  label=Antibody)) +
  geom_point(size= 3,alpha=0.7) +
  geom_text_repel(max.overlaps = 100) +
  theme_classic()



                              #~ Perform Mann-Whitney ~#
mw.data <- mm.ab[, colnames(mm.ab) %in% f.non.all.values$Antibody]

# Perform Mann-Whitney for the 2 antibodies that were not fited by the mixture models
mw.test <- c()
for(ab in mw.data){
  mw <- wilcox.test(ab ~ Y$Status)$p.value
  mw.test <- c(mw.test, mw)
}

mw.test <-  data.frame(mw.test)
mw.test$antibody <-  colnames(mw.data)
colnames(mw.test) <- c("p.value", "antibody")
print(paste("All", dim(subset(mw.test, mw.test$p.value<0.5))[1], " antibodies were statistically significant for the Mann-Whitney test"))

                            #~ Analyze the single population mixture models ~#
mm1 <- subset(f.all.values, f.all.values$Num.Comp == 1)
mm1.ab <- mm.ab[,colnames(mm.ab) %in% mm1$Antibody]

lik <- function(x){
  sn1 <- selm(x ~ 1, data =mm1.ab)
  sn2 <- selm(x ~ Y$Status , data = mm1.ab)
  a <- lmtest::lrtest(sn1,sn2)
  a$`Pr(>Chisq)`[2]
}

lk.vals <- lapply(mm1.ab, lik) %>% bind_rows()
lk.vals <- t(lk.vals)
lk.vals <- data.frame(lk.vals)
colnames(lk.vals)  <- "p.value"
print(paste(dim(subset(lk.vals,lk.vals$p.value< 0.05))[1], "out of the", table(f.all.values$Num.Comp)[1], "antibodies best fitted by a single component were stastistically significant"))


        
                      #~ Work on the 2 sub-populations Mixture-models ~#

#Here we are using the Chi-square test to achieve the optimal cutpoint
#Select the antibodies  with two components
mm.k2 <- subset( f.all.values, f.all.values$Num.Comp == 2)
X <- mm.ab[, colnames(mm.ab) %in% mm.k2$Antibody]

f.pval <- data.frame()                        # Data.frame to store p.values from single antibody
all.pval <- data.frame()                      # Data.frame to store best p.value fo each antibody
sens_spe.vals <-  data.frame()                # Dataf.rame to store sensitivity and specificity values
new.data <-  data.frame("Status" = Y$Status)  # Data.frame to store the dichotomized values

# Function with the algorithm itself
for (i in 1:ncol(X)) {
  ab <- X[,i]   # Define column (antibody)
  
  for (j in 1:length(ab)) {   # For each value of the antibody:
    trs <-  ab[j]             # Define the antibody value as a threshold
    ss <- ifelse(ab > trs, "Seropositive", "Seronegative") # Patients with Values above the treshold are sero +  and below sero -
    tab <- table(ss, Y$Status)  # Obtain the number of individuals for each class
    
    if( nrow(tab) < 2 | ncol(tab) <2 ){ # Put every important value into a dataframe
      p.val <- data.frame("p.value"= NA, "threshold" = trs, "Seronegative" = table(ss)[1], "Seropositive" = table(ss)[2], "Seropositivity" = (table(ss)[2])/121 )
    }else{
      p.val <- data.frame("p.value"= chisq.test(tab)$p.value, "threshold" =trs, "Seronegative"= table(ss)[1], "Seropositive" = table(ss)[2],"Seropositivity" = (table(ss)[2])/121)
    }
    f.pval <- rbind(f.pval, p.val)  # Bind the above information for every possible threshold
    
  }
  b.pval <- f.pval[which.min(f.pval$p.value),]    # Select the threshold for which the p-value is lowest ( higher discriminatory power)
  all.pval <- rbind(all.pval, b.pval)             # Bind this information together
  new.val <-  ifelse(ab > b.pval$threshold, 1, 2) # Dichotomize the data
  new.data <- cbind(new.data, new.val)            # Build the data.frame with the dichotomized values
  tab <- table(new.val, Y$Status)                 # Create a 2x2 table to obtain sensitivity and specificity
  sn.sp <- data.frame("Sensitivity" = tab[1]/sum(tab[,1]),
                      "Specificity" = tab[4]/sum(tab[,2])) # Get sensitivity and specificity
  sens_spe.vals <- rbind(sens_spe.vals, sn.sp)  # Bind everything on a data.frame
  
  f.pval <-  data.frame() # Clear the data.frame
  b.pval <- data.frame()  # Clear the data.frame
}

# Attirbute the respective col/rownames
rownames(all.pval) <- colnames(X)
colnames(new.data)[2:ncol(new.data)] <- colnames(X)
rownames(sens_spe.vals) <-  colnames(X)


# Obtain all the p.values across all test in a single table
p.values_data <- rbind(data.frame("Antibody"= rownames(t.vals),
                                  "p.value" = t.vals$p.value, 
                                  "Test"= "T-test"),
                       data.frame("Antibody" = mw.test$antibody,
                                  "p.value" = mw.test$p.value,
                                  "Test" = "Man-Whitney"),
                       data.frame("Antibody"= rownames(lk.vals),
                                  "p.value" = lk.vals$p.value, 
                                  "Test" = "Likelihood"),
                       data.frame("Antibody"= rownames(all.pval),
                                  "p.value" = all.pval$p.value, 
                                  "Test" = "Chi-squared"))

dim(p.values_data)



# Plot the p.values graph
p.values_data$Significance <- ifelse(p.values_data$p.value <0.05, "p-value < 0.05", "p-value > 0.05")

# Plot using ggplot2
ggplot(p.values_data , aes(x=Antibody, y= log10(p.value), fill = Significance ,color = Significance, shape = Test)) +
  geom_point(size =4) + geom_hline(yintercept = log10(0.05), lty=2, col = "black") +
  coord_cartesian( ylim = c(-5, 0)) +
  theme_bw() + 
  scale_fill_manual(values=c("#3498db","#ff5555")) +
  scale_color_manual(values=c("#3498db","#ff5555")) +
  scale_y_continuous( breaks = c( 0, -1,-2,-3,-4,-5), labels= c(1, 0.1, 0.01, 0.001 , 0.0001, 0.00001)) +
  ylab("p.value") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
                          axis.text=element_text(size=12),
                          panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank())

# Adjust the p.value
p.values_ad<- p.adjust(p.values_data$p.value, method = "BY", n=length(p.values_data$p.value))

p.values_ad <-  data.frame("Antibody"= p.values_data$Antibody, 
                           "p.value" = p.values_ad,
                           "Test" = p.values_data$Test,
                           "Significance" = ifelse(p.values_ad <0.05,"p-value < 0.05", "p-value > 0.05"))

ggplot(p.values_ad , aes(x=Antibody, y= log10(p.value), fill = Significance ,color = Significance, shape = Test)) +
  geom_point(size =4) + geom_hline(yintercept = log10(0.05), lty=2, col = "black") +
  coord_cartesian( ylim = c(-5, 0)) +
  theme_bw() + 
  scale_fill_manual(values=c("#3498db","#ff5555")) +
  scale_color_manual(values=c("#3498db","#ff5555")) +
  scale_y_continuous( breaks = c( 0, -1,-2,-3,-4,-5), labels= c(1, 0.1, 0.01, 0.001 , 0.0001, 0.00001)) +
  ylab("p.value") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
                          axis.text=element_text(size=12),
                          panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank())



# Subset the statistically significant antibodies
p.value.sig <- subset(p.values_ad, p.values_ad$p.value < 0.05 )
table(p.value.sig$Test)
print(paste("After the p.value correction with the BY test we have a total of", dim(p.value.sig)[1], "significant antibodies"))

# Obtain the final with all the significant antibodies
mw.abs <-  mw.data[, colnames(mw.data) %in% p.value.sig$Antibody]
l.abs <- mm1.ab[, colnames(mm1.ab) %in% p.value.sig$Antibody]
p.abs <- new.data[, colnames(new.data) %in% p.value.sig$Antibody]

for (i in 1:ncol(p.abs)) {
  p.abs[,i] <-  ifelse(p.abs[,i] == 1,1,0)
  
}

# Bind everything
data <- cbind("Status"= Y$Status, mw.abs, l.abs, p.abs)

# Use the SuperLearner
x_train <- data[, !colnames(data) %in% "Status" ]
y_train <- ifelse(data$Status == "protected",1,0)

# Implment the SuperLearner
set.seed(42)
#sl = CV.SuperLearner(Y = y_train, X = x_train, family = binomial(),
#                     V = nrow(x_train),
#                     SL.library = c("SL.glm", "SL.randomForest","SL.lda", "SL.qda", "SL.xgboost"))
# Review results
#saveRDS(sl, "parametric2.RDS")
sl <- readRDS("parametric2.RDS")


summary(sl)

#Obtain the weights
review_weights <- function(cv_sl) {
  meta_weights = coef(cv_sl)
  means = colMeans(meta_weights)
  sds = apply(meta_weights, MARGIN = 2,  FUN = sd)
  mins = apply(meta_weights, MARGIN = 2, FUN = min)
  maxs = apply(meta_weights, MARGIN = 2, FUN = max)
  sl_stats = cbind("mean(weight)" = means, "sd" = sds, "min" = mins, "max" = maxs)
  sl_stats[order(sl_stats[, 1], decreasing = TRUE), ]
}

# WEIGHTS
print(review_weights(sl), digits = 3)


# Obtain the roc curve 
roc <- roc(y_train , sl$SL.predict)
print(roc)

# Obtain the cut-offs
X <- data.frame("X"= sl$SL.predict, "Y" = y_train )

# Obtain the optimal cutpoints (point closest to the top left and sensitivity = specificity)
opt.roc01 <- optimal.cutpoints("X", "Y", tag.healthy = 0, method = "ROC01" , data = X)
opt.spse <- optimal.cutpoints("X", "Y", tag.healthy = 0, method = "SpEqualSe", data = X)

print("ROC01")
print(summary(opt.roc01))

print("SpEqualSe")
print(summary(opt.spse ))

#Obtain the summary
s.roc01 <- summary(opt.roc01) 
s.sesp <- summary(opt.spse)   

# Obtain the thresold point
t.roc01 <- s.roc01$p.table$Global$ROC01[[1]][1]
t.sesp <- s.sesp$p.table$Global$SpEqualSe[[1]][1]


#PLot in the roc curve
plot(roc, legacy.axes = T )

#PLot the points
points(x= s.roc01$p.table$Global$ROC01[[1]][3], y=s.roc01$p.table$Global$ROC01[[1]][2], col ="#F8766D", pch = 19, cex=1.5)
text(x=s.roc01$p.table$Global$ROC01[[1]][3] + 0.1, y =s.roc01$p.table$Global$ROC01[[1]][2] +0.05, paste("ROC01:",round(s.roc01$p.table$Global$ROC01[[1]][1],3)), col= "#F8766D")

points(x= s.sesp$p.table$Global$SpEqualSe[[1]][3], y=s.sesp$p.table$Global$SpEqualSe[[1]][2], col ="#619CFF", pch = 19, cex=1.5)
text(x=s.sesp$p.table$Global$SpEqualSe[[1]][3] - 0.20, y =s.sesp$p.table$Global$SpEqualSe[[1]][2] -0.075, paste("SpEqualSe:",round(s.sesp$p.table$Global$SpEqualSe[[1]][1],3)), col= "#619CFF")


# Obtain classification labels
sl_loo_roc01_labels <- ifelse(sl$SL.predict >= t.roc01 , 1, 0)
sl_loo_sesp_labels <- ifelse(sl$SL.predict >= t.sesp , 1, 0)

#Plot confusion matrix for "ROC01"
cf <- confusionMatrix( as.factor(sl_loo_roc01_labels), as.factor(y_train), positive = "1")$table
tb <- matrix(data =c(cf[4],cf[3],cf[2],cf[1]), 2)
fourfoldplot(tb, color = c("#FFCCCC", "#ff5555"),
             conf.level = 0,std ="ind.max")

#Plot confusion matrix for "SpEqualSe"
cf <- confusionMatrix(as.factor(sl_loo_sesp_labels), as.factor(y_train),positive = "1")$table
tb <- matrix(data =c(cf[4],cf[3],cf[2],cf[1]), 2)
fourfoldplot(tb, color = c("#CCFFFF", "#18bc9c"),
             conf.level = 0,std ="ind.max")


# Obtain model metrics according to hte 'ROC01' threshold
library.pred <-  sl$library.predict[,-4]

measure_individual_performance <- function(sl){
  output_list <- list()
  for(i in 1:ncol(library.pred)){
    prediction_vector <- library.pred[,i] # Obtain the predictive values fore each model
    prediction_labels <- ifelse(prediction_vector > s.roc01$p.table$Global$ROC01[[1]][1], 1, 0)
    conf <- caret::confusionMatrix(data = factor(prediction_labels, levels = c(1, 0)),
                                   reference = factor(y_train, levels = c(1, 0)),
                                   mode = 'everything',
                                   positive = "1")
    sl_roc <- roc(y_train~prediction_vector)
    metrics <- c(
      AUC = sl_roc$auc[1],
      conf$overall[1],
      conf$byClass["Sensitivity"],
      conf$byClass["Specificity"],
      conf$byClass["Precision"],
      conf$byClass["Recall"],
      conf$byClass["F1"]
    )
    output_list[colnames(sl$library.predict)[i]] <- list(metrics)
  }
  output_list
}
print("For the 'ROC01' metric the results are the following:")
print(data.frame(measure_individual_performance(sl)))
```

# Obtain model metrics according to hte 'SpEqualSe' threshold
```{r}
measure_individual_performance <- function(sl){
  output_list <- list()
  for(i in 1:ncol(library.pred)){
    prediction_vector <- library.pred[, i] # Obtain the predictive values fore each model
    prediction_labels <- ifelse(prediction_vector > s.sesp$p.table$Global$SpEqualSe[[1]][1], 1, 0)
    conf <- caret::confusionMatrix(data = factor(prediction_labels, levels = c(1, 0)),
                                   reference = factor(y_train, levels = c(1, 0)),
                                   mode = 'everything',
                                   positive = "1")
    sl_roc <- roc(y_train~prediction_vector)
    metrics <- c(
      AUC = sl_roc$auc[1],
      conf$overall[1],
      conf$byClass["Sensitivity"],
      conf$byClass["Specificity"],
      conf$byClass["Precision"],
      conf$byClass["Recall"],
      conf$byClass["F1"]
    )
    output_list[colnames(sl$library.predict)[i]] <- list(metrics)
  }
  output_list
}
print("For the 'SpEqualSe' metric the results are the following:")
print(data.frame(measure_individual_performance(sl)))
```

