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
library(ranger)
library(ggwaffle)
library(hrbrthemes)
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

# Specify datasete for Random Forest
rf_data <-  data.frame(X,Y)
rf_data$Status <-  as.factor(rf_data$Status)

# Random Forest
#control <- trainControl(method='repeatedcv', 
#                        number=10, 
#                        repeats=100,
#                        search = "random",
#                        allowParallel = TRUE,
#                        savePredictions = T,
#                        classProbs = T,
#                        summaryFunction=twoClassSummary
#                        )
#
### Random search ##
#rf_random <- train(Status ~ .,
#                   data = rf_data,
#                   method = 'rf',
#                   trControl = control)

# Save the data and load it
#saveRDS(rf_random,"rf_preliminary.RDS")
rf_final <-  readRDS("rf_preliminary.RDS")
rf_final$results

# plot
plot(rf_final)
rf_final$bestTune

max(rf_final$results$ROC)

selectedIndices <- rf_final$pred$mtry == 1
plot.roc(rf_final$pred$obs[selectedIndices],
         rf_final$pred$protected[selectedIndices])

roc.curve <- roc(rf_final$pred$obs[selectedIndices],rf_final$pred$protected[selectedIndices] )
roc.curve
plot(roc.curve, col = "#18bc9c", lwd = 4, legacy.axes = T)

plot(varImp(rf_final))

vi <- varImp(rf_final)
vi <-vi$importance
vi


varImp

importance(rf_final)

# Order the vi data
ordered.vi <-  data.frame(vi, "Antibody"= rownames(vi))
ordered.vi <- ordered.vi [order(ordered.vi$Overall, decreasing = T),]
ordered.vi$pos <- 1:nrow(ordered.vi)

waffle_data <-  data.frame()

# Specify number of iterations
j <-  nrow(ordered.vi)
for(i in 1:j){
  ab <- rownames(ordered.vi)[i] # Colect the names
  vec <- data.frame("antibody" = ab, "position" = ordered.vi$pos[i] , "values" = seq(1,100,2))# Create dataframe
  vec$cond <-  as.factor(ifelse(vec$values <= ordered.vi$Overall[i] , "bellow", "above"))
  
  waffle_data <-  rbind(waffle_data, vec)
  vec <-  data.frame()
}

ggplot(waffle_data, aes(values, reorder(antibody, -position), colour = cond)) + 
  geom_waffle(tile_shape = 'circle', size = 3) +
  scale_colour_manual(values = c("#18bc9c","grey95")) +
  ylab("Antibody") +
  xlab("Importance") +
  theme_minimal() +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_x_continuous(breaks = c(0,20,40,60,80,100))+
  geom_vline(xintercept = 20, linetype = "dashed", color= "grey20", size = 1)
  
#ggplot(waffle_data, aes(values, reorder(antibody, -position), colour = cond)) + 
#  geom_waffle(tile_shape = 'circle', size = 3) +
#  scale_colour_manual(values = c("grey","grey95")) +
#  ylab("Antibody") +
#  theme_minimal() +
#  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(waffle_data, aes(values,  reorder(antibody, -position), fill = cond)) + 
  geom_waffle(size = 1) +
  scale_colour_waffle() +
  scale_fill_manual(values = c("orange","grey95")) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())


