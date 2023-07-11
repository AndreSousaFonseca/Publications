# Load libraries:
library(MASS)             # Perform basic data analysis
library(dplyr)            # Better handle data
library(pROC)             # Perform an plot the ROC curve and obtain the AIC   
library(tidyr)            # Better handle data
library(ggplot2)          # Plot the data
library(caret)            # Perform Confusion Matrices
library(doParallel)       # Use parallel processing for faster run times
library(OptimalCutpoints) # Obtain the optimal outpoints
library(SuperLearner)     # Perform predictive analysis
library(ggrepel)          # Avoid overlapping text in graph
set.seed(42)              # Specify a seed

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

# Perform chi-squared test
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

# Gather the useafull information to be plotted
p.vals <-  data.frame("Antibody" = rownames(all.pval),
                      "p.value" = all.pval$p.value,
                      "Significance" = ifelse(all.pval$p.value <0.05, "p-value < 0.05","p-value > 0.05"))
# Plot using ggplot2
ggplot(p.vals, aes(x=Antibody, y= log10(p.value), fill = Significance ,color = Significance)) +
  geom_point(shape = 21, size =4) + geom_hline(yintercept = log10(0.05), lty=2, col = "black") +
  coord_cartesian( ylim = c(-5, 0)) +
  theme_bw() + 
  scale_fill_manual(values=c("#3498db","#ff5555")) +
  scale_color_manual(values=c("#3498db","#ff5555")) +
  scale_y_continuous( breaks = c( 0, -1,-2,-3,-4,-5), labels= c(1, 0.1, 0.01, 0.001 , 0.0001, 0.00001)) +
  ylab("p.value") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
                          axis.text=element_text(size=12),
                          panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank())

# Plot the sensitivity against specificity scatter plot for each antibody
# Add a column called "Antibody" wih will have the Antibody names
sens_spe.vals$Antibody <-  rownames(sens_spe.vals)

# Plot the sensitivity against specificity scatter plot for each antibody
ggplot(sens_spe.vals, aes(1-Specificity, Sensitivity,  label=Antibody)) + geom_point(size= 3,alpha=0.7) + geom_text_repel(max.overlaps = 100) + theme_classic() + ylim(0,1) + xlim(0,1) + geom_abline(intercept = 0, slope = 1, col= "darkgrey", lty = 2) 

# Gather the data to plot a ggplot2 barplot
ss.plot <-  gather(sens_spe.vals, "metric", "value", -Antibody)

# Plot
ggplot(ss.plot, aes(x= Antibody, y = value, fill=metric)) + 
  geom_bar(position="dodge", stat="identity", alpha =0.5)+  theme_test() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ylim(0,1) 

# Obtain the corrected p-values
p.values_ad <- p.adjust(all.pval$p.value, method = "BY", n=length(all.pval$p.value))

# Print the  number of significant antibodies after correction
print(paste("There are", sum(p.values_ad <0.05), "significant antibodies after Benjamini-Yekutieli correction"))

# Gather the information in a data.frame with the antibody names
adjusted.p <-  data.frame("Antibody" = rownames(all.pval),
                          "p.value" = p.values_ad,
                          "Significance" = ifelse(p.values_ad <0.05, "p-value < 0.05","p-value > 0.05"))

ggplot(adjusted.p, aes(x=Antibody, y= log10(p.value), fill = Significance ,color = Significance)) +
  geom_point(shape = 21, size =4) + geom_hline(yintercept = log10(0.05), lty=2, col = "black") +
  coord_cartesian( ylim = c(-5, 0)) +
  theme_bw() + 
  scale_fill_manual(values=c("#3498db","#ff5555")) +
  scale_color_manual(values=c("#3498db","#ff5555")) +
  scale_y_continuous( breaks = c( 0, -1,-2,-3,-4,-5), labels= c(1, 0.1, 0.01, 0.001 , 0.0001, 0.00001)) +
  ylab("p.value") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
                          axis.text=element_text(size=12),
                          panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank())

# Select only the statistically significant antibodies after correction
sig.abs <- subset(adjusted.p , adjusted.p$p.value <0.05)
dim(sig.abs)

# subset in the dichotomized (new.data) dataframe this significant antibodies
data <- new.data[, colnames(new.data) %in% sig.abs$Antibody]

# Transform the data from "1/2" to "1/0" to perform the SuperLEarner (Personal Preference)
for(i in 1:ncol(data)){
  data[,i] <-  ifelse(data[,i] == 1,1,0)
}

# Add the variable Status
data <- data.frame("Status" = Y$Status, data)

x_train <- data[, !colnames(data) %in% "Status" ]
y_train <- ifelse(data$Status == "protected",1,0)

# Implment the SuperLearner
set.seed(42)
#sl = CV.SuperLearner(Y = y_train, X = x_train, family = binomial(),
#                     V = nrow(x_train),
#                     SL.library = c("SL.glm", "SL.randomForest","SL.lda", "SL.qda", "SL.xgboost"))
# Caching result
#saveRDS(sl, "opt.RDS")
sl<- readRDS("opt.RDS")

# Review results
summary(sl)
sl

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
ci(roc)



# Obtain the cut-offs
X <- data.frame("X"= sl$SL.predict, "Y" = y_train )

# Obtain the optimal cutpoints (point closest to the top left and sensitivity = specificity)
opt.roc01 <- optimal.cutpoints("X", "Y", tag.healthy = 0, method = "ROC01" , data = X)
opt.spse <- optimal.cutpoints("X", "Y", tag.healthy = 0, method = "SpEqualSe", data = X)

print("ROC01")
print(summary(opt.roc01))

print("SpEqualSe")
print(summary(opt.spse))

#Obtain the summary
s.roc01 <- summary(opt.roc01) 
s.sesp <- summary(opt.spse)   

# Obtain the thresold point
t.roc01 <- s.roc01$p.table$Global$ROC01[[1]][1]
t.sesp <- s.sesp$p.table$Global$SpEqualSe[[1]][1]


#PLot in the roc curve
plot(roc, legacy.axes = T )

#PLot the points
sl_loo_roc01_labels <- ifelse(sl$SL.predict >= t.roc01 , 1, 0)

#Plot confusion matrix for "ROC01"
cf <- confusionMatrix( as.factor(sl_loo_roc01_labels), as.factor(y_train), positive = "1")$table
tb <- matrix(data =c(cf[4],cf[3],cf[2],cf[1]), 2)
fourfoldplot(tb, color = c("#FFCCCC", "#ff5555"),
             conf.level = 0,std ="ind.max")


# Obtain the predictions for each individaul model
sl$library.predict



roc_glm <- roc(y_train , sl$library.predict[,1])
plot(roc_glm, col= "#f8a800", lwd=2,  legacy.axes = T)

roc_rf <- roc(y_train , sl$library.predict[,2])
plot(roc_rf, col= "#18bc9c", add=T, lwd=2)

roc_xgb <- roc(y_train , sl$library.predict[,5])
plot(roc_xgb, col= "#3498db", add=T, lwd=2)


plot(roc,  add=T, lwd=1, lty =1)




