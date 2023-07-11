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

# Calculate the wilcoxon p.values
output <-  c()
for(i in 1:ncol(X)){
  output <- c(output,wilcox.test(X[,i]~ Y$Status)$p.value)
}

# Create a data frame with the wilcoxon p.values
p.vals <- data.frame("Antibody" = colnames(X),
                     "p.value" = output, 
                     "Significance" = ifelse(output < 0.05, "p-value < 0.05", "p-value > 0.05"))

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

# Adjust the p.values using Benjamini-Yekutieli
adjusted.p <-p.adjust(output, method = "BY", n=length(output))
sum(adjusted.p <0.05)
p.adjust <- data.frame("Antibody" = colnames(X),
                       "p.value" = adjusted.p, 
                       "Significance" = ifelse(adjusted.p < 0.05, "p-value < 0.05", "p-value > 0.05") )

# Plot using ggplot2
ggplot(p.adjust , aes(x=Antibody, y= log10(p.value), fill = Significance ,color = Significance)) +
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


# Perform Spearman correlation on all antibodies
m <-  matrix(data= NA, nrow = 36, ncol = 36)
l <-  ncol(X)
for(i in 1:l){
  for(j in 1:l){
    ab1 <-  X[,i]
    ab2 <-  X[,j]
    m[i,j] <- cor(ab1,ab2,use="everything", method = "spearman")
  }
}

# Convert to data.frame
m_df <- as.data.frame(m)

# Rename columns and rows
colnames(m_df) <-  colnames(X) ; rownames(m_df) <-  colnames(X)

# Remove lower and diagonial data.frames
m_df[lower.tri(m_df,diag = T)] <- NA

# Bind everything in a single column
spdata<- data.frame()
for(i in 1:ncol(m_df)){
  values <- c(m_df[,i])
  spdata <- rbind(spdata, data.frame(values))
}

# remove na data
spdata <- na.omit(spdata)

# Plot the boxplot
boxplot(spdata)

# Subset the significant antibodies
sig.abs <- subset(p.adjust, p.adjust$p.value <0.05)
sig.abs <- X[, colnames(X) %in% sig.abs$Antibody]

# Create the predictor and outcome variables
x_train <- sig.abs
y_train <- ifelse(Y$Status == "protected",1,0)

# Implment the SuperLearner
set.seed(42)
#sl = CV.SuperLearner(Y = y_train, X = x_train, family = binomial(),
#                    V = nrow(x_train),
#                    SL.library = c("SL.glm", "SL.randomForest","SL.lda", "SL.qda", "SL.xgboost"))
# Caching result
#saveRDS(sl,"pragmatic.RDS")
sl <-  readRDS("pragmatic.RDS")

# Review results.
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
points(x= s.roc01$p.table$Global$ROC01[[1]][3], y=s.roc01$p.table$Global$ROC01[[1]][2], col ="#F8766D", pch = 19, cex=1.5)
text(x=s.roc01$p.table$Global$ROC01[[1]][3] + 0.1, y =s.roc01$p.table$Global$ROC01[[1]][2] +0.05, paste("ROC01:",round(s.roc01$p.table$Global$ROC01[[1]][1],3)), col= "#F8766D")

points(x= s.sesp$p.table$Global$SpEqualSe[[1]][3], y=s.sesp$p.table$Global$SpEqualSe[[1]][2], col ="#619CFF", pch = 19, cex=1.5)
text(x=s.sesp$p.table$Global$SpEqualSe[[1]][3] - 0.20, y =s.sesp$p.table$Global$SpEqualSe[[1]][2] -0.01, paste("SpEqualSe:",round(s.sesp$p.table$Global$SpEqualSe[[1]][1],3)), col= "#619CFF")

sl_loo_roc01_labels <- ifelse(sl$SL.predict >= t.roc01 , 1, 0)
sl_loo_sesp_labels <- ifelse(sl$SL.predict >= t.sesp , 1, 0)

#Plot confusion matrix for "ROC01"
cf <- confusionMatrix( as.factor(sl_loo_roc01_labels), as.factor(y_train), positive = "1")$table
tb <- matrix(data =c(cf[4],cf[3],cf[2],cf[1]), 2)
fourfoldplot(tb, color = c("#FFCCCC", "#ff5555"),
             conf.level = 0,std ="ind.max")

#Plot confusion matrix for "SpEqualSe"
cf2 <-confusionMatrix(as.factor(sl_loo_sesp_labels), as.factor(y_train),positive = "1")$table
tb <- matrix(data =c(cf2[4],cf2[3],cf2[2],cf2[1]), 2)
fourfoldplot(tb, color = c("#CCFFFF", "#18bc9c"),
             conf.level = 0,std ="ind.max")

measure_individual_performance <- function(sl){
  output_list <- list()
  for(i in 1:ncol(sl$library.predict)){
    prediction_vector <- sl$library.predict[, i] # Obtain the predictive values fore each model
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

measure_individual_performance <- function(sl){
  output_list <- list()
  for(i in 1:ncol(sl$library.predict)){
    prediction_vector <- sl$library.predict[, i] # Obtain the predictive values fore each model
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


wei <-data.frame(measure_individual_performance(sl))[1,]
colnames(wei) <-  c("LRM","RF", "LDA", "QDA", "XGB")
wei.gg <- gather(wei,"Model", "Values", LRM:XGB)

ggplot(wei.gg, aes(Model,Values, fill= Model)) +
  geom_bar(stat="identity") + coord_cartesian(ylim = c(0.7, 0.73)) +
  scale_fill_manual(values = c("#f8a800","#f87f2a","#ff5555","#18bc9c","#3498db"))+
  theme_light() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
  
  

               