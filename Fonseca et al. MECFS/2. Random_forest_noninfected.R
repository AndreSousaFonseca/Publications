# Random Forest (HCs vs. non-EBV infect individuals)

# Load libraries
library(caret)        # Cross validation analysis
library(caTools)      # Split the data
library(pROC)         # Roc curve   
library(tidyr)        # Arrange the data
library(dplyr)        # For data cleaning
library(gplots)       # Heatmap
library(MASS)         # Use several fucntions to handle the data 
library(AID)          # Variety of functions that facilitate the analyzes
library(ggplot2)      # Plot graphs
library(readxl)       # Read the excell files
library(doParallel)       # Use parallel processing for faster run times
library(OptimalCutpoints) # Obtain the optimal outpoints
library(SuperLearner)     # Perform predictive analysis
library(ggrepel)          # Avoid overlapping text in graph
library(ranger)
library(writexl)
set.seed(42)        # Set a seed for reproducibility

# Specificy the number of cores for parallel processing
numCores <- detectCores()     # Detect the number of cores
registerDoParallel(numCores)  # Use multicore, set to the number to our cores

# Load and process the data
dataset <- read_excel("Data_CFS.xlsx")
dataset <- t(dataset)               # Transpose the dataset (antibodies in columns and patients in rows)
dataset <- as.data.frame(dataset)   # Set as data.frame

# Set the colnames as the first row
colnames(dataset) <- dataset[1,]

# Treat the colnames  (omit the " 00_ " in the begining of the EBV virus)
colnames(dataset) <- substring(colnames(dataset), 4)

# Check how many unique antibodt name swe have
print(paste("There are a total of",length(unique(colnames(dataset))), "unique antibody names"))

# Obtain the indexes for the duplicated antigens
duplicated_index <- duplicated(colnames(dataset))

# Split the data into duplicated and non-duplicated
duplicated_data <- dataset[,duplicated_index]
unique_data <- dataset[,!duplicated_index]

# Obtain the unique names of the repeated antibodies
duplicated_abs <- unique(colnames(duplicated_data))

# Let's now get the duplicated antigens and add a number in front to 
# distinguish them. Bare in mind that the duplicated antigens are
# already the repetitions of the value in unique_data. To prove this:
which(colnames(unique_data) == duplicated_abs[1]) 
# Indeed the names in duplicated_ags are in unique_data

# So let's add the number in front of the antibodies
different_data <- data.frame("Index" = c(1:nrow(unique_data)))

for (ab in duplicated_abs) { # For each antibody
  dataset <- data.frame(duplicated_data[,colnames(duplicated_data) == ab]) #obtain all duplicates
  n <- ncol(dataset) # Check the number of replicates
  
  for( k in 1:n){   # For each replicate add a "_" and a number
    colnames(dataset)[k] <-  paste0(ab[k],"_", k)
  }
  different_data <-  cbind(different_data, dataset) #Bind everything to the different_data
}

# Drop the index column
different_data <-  different_data[, !colnames(different_data) %in% "Index"]

# Bind the different_data and the unique_data together
dataset <-cbind(unique_data, different_data)
dim(dataset)

# To confirm that everything was performed well
length(unique(colnames(dataset)))

# Store the sequence information
seq_data <- dataset

# Remove the Name, Sequence, AG876, B95-8, GD1, Cao, Raji, P3HR-1 rows
dataset <- dataset[! rownames(dataset) %in% c("Name", "Sequence", "AG876", "B95-8", "GD1", "Cao", "Raji", "P3HR-1") ,]

# Transform the the values from cathegorical to numerical
n <-  ncol(dataset)
for (i in 1:n) {
  dataset[,i] <-  as.numeric(dataset[,i])
}

# Add the meta data
# Load the meta data
meta <- read_excel("epi_data_age_gender.xlsx")
metacs <- read.csv("Data_Epi_Data_CFS_Infection.csv")

# Merge both meta datasets
metas <- merge(x = meta , metacs, by= c("ID", "age", "sex"), all = TRUE)

# Chech if the ID matches the rownames of our dataset
which(rownames(dataset) !=  metas$ID) 

# Order the metas dataset to match the rows of dataset
metas <- metas[ match(rownames(dataset), metas$ID),]
which(rownames(dataset) !=  metas$ID) 
# No everthing is ok, we can bidn everything

# Add everything into a single data.frame
dataset <-data.frame("Age" = metas$age,
                     "sex" = metas$sex,
                     "Infection" =metas$infection.triggered.onset,
                     dataset
)


# Add a varaible identifying healthy and CFS patients
dataset$Status <- as.factor(ifelse(startsWith(rownames(dataset), "HEALTHY") == TRUE, "healthy", "cfs"))


# Subset the 50 healthy and 38 non-infected
data.h <- subset(dataset, dataset$Status == "healthy")
data.c <- subset(dataset, dataset$Status == "cfs")
data.c <- subset(data.c, data.c$Infection != "ja" | is.na(data.c$Infection))

data <- rbind(data.h, data.c)
dim(data)

# Drop unnecessary variables and perform data splitting
dataset <-  data[, !colnames(data) %in% c("Age", "sex", "Infection")]

# Perform a 90% training and 10% testing validation 
spliting <- sample.split(dataset$Status, SplitRatio = 0.9)
train_set <-  subset(dataset, spliting == TRUE) # Obtain the train set
test_set <-  subset(dataset, spliting == FALSE) # obtain the test set

# Obtain the train and test sets antibody values only
training <-  train_set[, !colnames(train_set) %in% "Status"]
testing <- test_set[, !colnames(test_set) %in% "Status"]


# Perform the Random Forest
# Specify a dataset that will contain all the relevant data
#info_data <-  data.frame()
#feature_imp <-  data.frame()

# Specify hyperparameters
#num_trees <- c(100, 500, 750, 1000, 2000) # Define the number of trees
#num_mtry <- sample(1:100, 50)             # Define the number of variables to possibly split at in each node
#min_node <- seq(1,10)                     # Define the minimal node size
#iteration <- 1                            # Set a iteration value

# Specify paremeters
#l <- length(num_trees)
#k <- ncol(training)

# define variables used in RandomForest
#info_tree <-  c()
#info_mtry <-  c()
#info_node <-  c()
#info_acc <-  c()
#info_imp <-  c()
#info_ab <-  c()

# Run the random forest loop
#rf_info_data <- foreach(ntree = num_trees[1:l], .combine = rbind, .packages="ranger" ) %dopar% {
#  for(mtry  in  num_mtry){
#    for(node  in  min_node){
#    # Set a seed for reproducibility
#      set.seed(42)
#      #Perform ranger
#      rangerforest <- ranger(Status ~ . ,
#                             data = train_set,
#                             num.trees = ntree,
#                             mtry = mtry,
#                             min.node.size = node,
#                             importance = "impurity",
#                             num.threads = numCores,
#                             classification =  TRUE)
#      
#      # Bind the information obtained to the info data
#      info_tree <-  c(info_tree, rep(rangerforest$num.trees, k))    # Obtain the num_tree value
#      info_mtry <-  c(info_mtry, rep(rangerforest$mtry, k))         # Obtain the mtry value
#      info_node <-  c(info_node, rep(rangerforest$min.node.size,k)) # Obtain the node values
#      
#      accuracy <- (rangerforest$confusion.matrix[1] +
#                     rangerforest$confusion.matrix[4]) / rangerforest$num.samples
#      info_acc <- c(info_acc, rep(accuracy,k))                      # Obtain the accuracy values
#      
#      ab.name <- data.frame(importance(rangerforest))               
#      info_ab <-  c(info_ab, rownames(ab.name))                     # Obtain the antibody names
#      info_imp <- c(info_imp, importance(rangerforest))             # Obtain hte importance values
#    }
#  }
#  info <- data.frame(
#    "Antibody" = info_ab,
#    "Importance"= info_imp,
#    "Num.tree" = info_tree,
#    "Mtry" = info_mtry,
#    "Node.size" = info_node,
#    "Accuracy" =info_acc)
#}
#
## Obtain a itterarion identifier.
#iteration <- c()
#l <- nrow(subset(rf_info_data, rf_info_data$Antibody == rf_info_data$Antibody[1]))
#k <- length(unique(rf_info_data$Antibody))
# Each RandomForest run gives me an output of k = 3054.
# Therefore I will have to replicate a single interarion 3054 times.
# Then we performed 2500 runs. So we need to have a total of 2500 different 
# iterations.
#for(i in 1:l){
#  iteration <-  c(iteration, rep(i, k))
#}
# Add the iterarion ot the results ouputed by the Random Forest
#rf_info_data$Iteration <- iteration

#saveRDS(rf_info_data, "rf_info_data_no.RDS") # Save the results form Random Froest
rf_info_data <- readRDS("rf_info_data_no.RDS") # Load the results from Random Forest
# Obtain the mean importancce for each antibody
#names <-  unique(rf_info_data$Antibody)
#l <-  length(names)
#Obtain the mean importance value for each antibody
#mean_imp <- foreach(name = names[1:l], .combine= "rbind") %dopar% {
#  v <- subset(rf_info_data, rf_info_data$Antibody == name)
#  s <- v$Importance
#  mean <-  data.frame("Antibody" = name, "Importance_Mean" = mean(s))
#}

# Save the data
saveRDS(mean_imp , "rf.mean_no.RDS")
mean_imp <- readRDS("rf.mean_no.RDS")
mean_imp$Protein <- stringr::str_extract(mean_imp$Antibody, "[^_]*_[^_]*")
length(unique(mean_imp$Protein)) # We confirm that we have 14 proteins

mean_imp$Antibody <- substr(mean_imp$Antibody, 5, 20)
mean_imp$Protein <- substr(mean_imp$Protein, 5, 20)


                                            # PLOTS #

#Plot te data for the importance. Any tendency of a protein being more important?
col <- c(viridis::viridis(n=14))
#col <- dichromat::colorschemes$BluetoOrangeRed.14

sup2 <-ggplot(mean_imp, aes(x = Antibody , y=Importance_Mean, color = as.factor(Protein))) +
  geom_point(size = 2, alpha = 0.7) +
  ylim(c(0.,0.4)) +
  theme_bw() + ylab("Mean Importance") +
  guides(color = guide_legend(title = "Protein")) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_color_manual(values =col)+
  ggtitle("HCs against ME/CFS with non-infectious or unknown disease trigger") +
  theme(plot.title = element_text(hjust = 0.5))

saveRDS(sup2, "sup2(B).RDS")

sup1 <-ggplot(mean_imp, aes(x = Protein , y=log10(Importance_Mean), fill = as.factor(Protein))) +
  geom_boxplot()+theme_bw() +
  ylim(c(-3,0)) +
  theme_bw() + ylab("Mena importance") +
  guides(fill = guide_legend(title = "Protein")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) + 
scale_fill_manual(values =col)+ theme(legend.position = "none") +
  ggtitle("HCs against ME/CFS with non-infectious or unknown disease trigger") +
  theme(plot.title = element_text(hjust = 0.5))

sup1
saveRDS(sup1, "sup1(B).RDS")


#ggplot(mean_imp, aes(x = Protein , y=Importance_Mean, colour = as.factor(Protein))) +
#  geom_boxplot()+theme_bw() +
#  ylim(c(0.,0.4)) +
#  theme_bw() + ylab("Importance") +
#  guides(colour = guide_legend(title = "Protein")) +
#  theme(panel.grid.major = element_blank(),
#        panel.grid.minor = element_blank()) #+ 
#  scale_colour_manual(values =col)



# Mean protein importance plot
col <- c(viridis::viridis(n=14))
proteins <-  unique(mean_imp$Protein)
k <-  length(proteins)
mean_prot <- foreach(prot = proteins[1:k], .combine="rbind") %dopar% {
  prt <-  subset(mean_imp, mean_imp$Protein == prot)
  mean <-  data.frame("Protein" = prot, "Importance_Mean" = mean(prt$Importance))
}
ggplot(mean_prot, aes(x = Protein , y=Importance_Mean, fill = as.factor(Protein))) +
  geom_bar(stat="identity") +theme_bw() +
  ylab("Mean Importance")+
  guides(fill = guide_legend(title = "Protein")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) + 
  scale_fill_manual(values =col)


# Is data normally distributed?
shapiro.test(mean_imp$Importance_Mean)

# Perform the Kruskal wallis
kw_res <- kruskal.test(Importance_Mean ~ Protein, mean_imp)
kw_res

dt <- dunnTest(Importance_Mean ~ Protein,
               data=mean_imp,
               method="bonferroni")

dt
write.csv(dt$res,"table_portein-wise_no.csv")

# As the p-value is greater than the significance level 0.05,
# we can conclude that there are no significant differences 
# between the groups in the model summary.



# Density plot
density.value <- density(mean_imp$Importance_Mean)
density.data <- data.frame("x" = density.value$x, "y"= density.value$y)

dens.values <- data.frame()
for (values in mean_imp$Importance_Mean) {
  min <- which(abs(density.data$x - values)==min(abs(density.data$x-values)))
  min.val <- data.frame(density.data[min,])
  dens.values <-  rbind(dens.values, min.val)
}
dens.values$names <-  mean_imp$Antibody

# Draw the 5 topmost important antibody barplot mean values
best.results <- mean_imp[order(mean_imp$Importance_Mean, decreasing = T),][1:10,]
best.results <- best.results[order(best.results$Importance_Mean, decreasing = F),]
best.results$pos <- c(1:nrow(best.results)) 

best.results$col <- rev(c(viridis::plasma(n=10)))

for(i in 1:length(best.results$Antibody)){
  if(endsWith(best.results$Antibody[i], "_1") == T){
    best.results$Antibody[i] <- gsub("_1", "*", best.results$Antibody[i])
  }
}


plb <- ggplot(best.results, aes(x = reorder(Antibody, pos), y = Importance_Mean, fill = as.factor(pos))) +  
  geom_bar(stat = "identity" ) + coord_flip() + guides(fill="none") + ylim(c(0, 0.4)) + 
  theme_classic()  + xlab("Antibody") + ylab("Mean Importance")+
  scale_fill_manual(values = best.results$col)
plb 


#ab_box <- rf_info_data[rf_info_data$Antibody %in% best.results$Antibody,]
#plg <- ggplot(ab_box, aes(x=reorder(Antibody, Importance), Importance, fill = reorder(Antibody, Importance))) +
#  geom_boxplot() + 
#  coord_flip() + xlab("Antibody") +
#  scale_fill_manual(values = best.results$co) + theme_classic()  + ylim(c(0, 1)) +
#  theme(legend.position = "none")  
#
#plg


# Obtain the most important
best.values <- mean_imp[order(mean_imp$Importance_Mean, decreasing = T),][1:10,"Antibody"]
best.results <- dens.values[dens.values$names %in% best.values,]
best.results <-  best.results[order(best.results$x, decreasing = T),]

# Obtain the density plot
plt <- ggplot(mean_imp, aes(x =Importance_Mean))+
  geom_histogram(aes(y =..density..), colour="lightgrey", fill="grey80",
                 alpha=0.5, position="identity", bins = 100)  +
  xlim(c(0,0.4)) + ylim(c(0,100)) + ylab("Density") +
  geom_density( color="grey25", alpha=0.5, lty = 1, lwd =1) +
  geom_jitter(data = best.results, aes(x,y), size = 3, color = c(viridis::plasma(n=10))) +
  theme_classic() + xlab("Mean Importance")

plt

library(ggpubr)
#ggarrange(plt, plb, labels = c("A", "B"),
#          ncol = 2, nrow = 1)

plt1 <- plt + annotation_custom(ggplotGrob(plb),
                        xmin = 0.2, ymin = 10,
                        xmax = 0.4, ymax = 100)


#plt + annotation_custom(ggplotGrob(plg),
#                        xmin = 0.1, ymin = 10,
#                        xmax = 0.2, ymax = 90)


# Order the dataset according to the importance and mean_importance
mean_imp <- readRDS("rf.mean_yes.RDS")
mean_imp$Protein <- stringr::str_extract(mean_imp$Antibody, "[^_]*_[^_]*")
top.imp <- mean_imp[order(mean_imp$Importance_Mean, decreasing = T),]

# We can see that the order of the antibody importance changes when we use a weighted
# mean to obtain the overall importance
# Obtain the dataets with the top important antiboides by their mean importance
i <- 2
train.subsets <-  list()
ncols <- 1 

while(ncols < 100){
  i <- i +1
  top <- top.imp[1:i,"Antibody"] # Select the top most antibodies~
  sub <-  training[, colnames(training) %in% top]
  cor <-  cor(sub,  use = "everything",method =  "spearman")
  drop <-  findCorrelation(cor, cutoff = .8)
  drop <-  colnames(cor)[drop]
  dataset <-  sub[,! colnames(sub) %in% drop]
  train.subsets <- c(train.subsets, list(dataset))
  ncols <-  ncol(dataset)
}


#Obtain the test_subsets sets
test_subsets <-  list()
for(subset in train.subsets){
  test <- test_set[, colnames(test_set) %in% colnames(subset)]
  test_subsets <- c(test_subsets, list(test))
}

# Apply the SuperLEarner for predictive analysis:
y_train <- ifelse(train_set$Status == "cfs", 1,0)
y_test <-  ifelse(test_set$Status == "cfs",1,0)

#slearn.subsets <-  list()
#for (data in train.subsets) {
#  sr = SuperLearner(Y = y_train, X = data, family = binomial(), 
#                    SL.library = c("SL.glmnet", "SL.randomForest","SL.lda", "SL.qda", "SL.xgboost"))
#  slearn.subsets <- c(slearn.subsets, list(sr))
#}
#saveRDS(slearn.subsets, "slearn.subsets.no.RDS")
slearn.subsets <- readRDS("slearn.subsets.no.RDS")


### Obtain the Train Set Accuracies ###

accuracies_training <-  data.frame()
for (model in slearn.subsets) {
  #Obtain the roc curve AUC
  roc <- roc(as.factor(y_train) , model$SL.predict[,1])
  
  # Obtain the predicted values 
  pred <-  ifelse(model$SL.predict[,1] > 0.5, 1 ,0)
  cf50 <- confusionMatrix( as.factor(pred), as.factor(y_train), positive = "1")
  
  
  #Obtain the pred values 
  X <- data.frame("X"= model$SL.predict[,1], "Y" = y_train )
  
  # Using the roc01
  opt.roc01 <- optimal.cutpoints("X", "Y", tag.healthy = 0, method = "ROC01", data = X)
  s.roc01 <- summary(opt.roc01) 
  roc01 <- s.roc01$p.table$Global$ROC01[[1]][1]
  pred <-  ifelse(model$SL.predict[,1] > roc01, 1 ,0)
  cf.roc <- confusionMatrix(as.factor(pred), as.factor(y_train), positive = "1")
  
  # Using the sesp
  opt.spse <- optimal.cutpoints("X", "Y", tag.healthy = 0, method = "SpEqualSe", data = X)
  s.sesp <- summary(opt.spse)
  sesp <-  s.sesp$p.table$Global$SpEqualSe[[1]][1]
  pred <-  ifelse(model$SL.predict[,1] > sesp, 1 ,0)
  cf.sesp <- confusionMatrix(as.factor(pred), as.factor(y_train), positive = "1")
  
  accs <- c(roc$auc,cf50$overall[1], cf.roc$overall[1], cf.sesp$overall[1], length(model$varNames))
  accuracies_training <-  rbind(accuracies_training, accs)
  
}

# Attribute column names
colnames(accuracies_training) <-  c("AUC","Fifty","ROC01","SpEqualSE", "Number_antibody")



# Obtain the predictive values for each model
models_training <-  data.frame()

for (model in slearn.subsets) {
  # Obtain the values for the models
  models_pred <- data.frame(model$library.predict)
  
  # Obtain each individual preedicitons 
  glm <- models_pred$SL.glmnet_All
  rf <- models_pred$SL.randomForest_All
  lda <- models_pred$SL.lda_All
  xgb <- models_pred$SL.xgboost_All
  
  
  # Obtain the predicted vlues
  glm.pred <-  ifelse(glm > 0.5, 1 ,0)
  rf.pred <-  ifelse(rf > 0.5, 1 ,0)
  lda.pred <-  ifelse(lda > 0.5, 1 ,0)
  xgb.pred <- ifelse(xgb > 0.5, 1 ,0)
  
  
  cf.glm <- confusionMatrix( as.factor(glm.pred), as.factor(y_train), positive = "1")
  cf.rf <- confusionMatrix( as.factor(rf.pred), as.factor(y_train), positive = "1")
  cf.lda <- confusionMatrix( as.factor(lda.pred), as.factor(y_train), positive = "1")
  cf.xbg <- confusionMatrix( as.factor(xgb.pred), as.factor(y_train), positive = "1")
  
  
  models_accs <- c(cf.glm$overall[1], cf.rf$overall[1], cf.lda$overall[1], cf.xbg$overall[1] ,length(model$varNames))
  models_training <-  rbind(models_training, models_accs)
}

colnames(models_training) <-  c("Glmnet","Random-Forest","LDA","XGB", "Number_antibody")

# Gather the data
models_training <- gather(models_training, "Method","Values", Glmnet : XGB)

# Plot the accuracies (50%)
sup3 <- ggplot(models_training, aes(x =Number_antibody, y =Values , color = Method)) + geom_line(size=1) + theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  xlab("Number of antibodies") + ylab("Accuracy") +
  geom_hline(yintercept=0.50,
             linetype="dashed", 
             color = "blue") +  geom_hline(yintercept=0.85,
                                           linetype="dashed", 
                                           color = "black") +
  scale_color_manual(values = c("#FD636B", "#ECB22E", "#2EB67D", "#36C5F0")) + labs(color='Model') +
  scale_y_continuous(breaks=c(0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1), limits =c(0.4,1)) + 
  ggtitle("HCs against ME/CFS with non-infectious or unknown disease trigger") +
  theme(plot.title = element_text(hjust = 0.5))
#1200 X  700

saveRDS(sup3, "sup3(B).RDS")


### Obtain the Test Set Accuracies ####
r <- 1
accuracies_testing <-  data.frame()

for (model in slearn.subsets) {
  
  predicted <- predict(model, test_subsets[[r]], onlySL = TRUE)
  
  # Obtain the roc curve
  roc <- roc(as.factor(y_test) , predicted$pred[,1])
  
  # Obtain the predicted values
  pred <-  ifelse(predicted$pred > 0.5, 1 ,0)
  cf50 <- confusionMatrix(as.factor(pred), as.factor(y_test), positive = "1")
  
  # Obtain the pred values
  X <- data.frame("X"= model$SL.predict[,1], "Y" = y_train )
  
  # Using the roc01
  opt.roc01 <- optimal.cutpoints("X", "Y", tag.healthy = 0, method = "ROC01", data = X)
  s.roc01 <- summary(opt.roc01) 
  roc01 <- s.roc01$p.table$Global$ROC01[[1]][1]
  pred <-  ifelse(predicted$pred > roc01, 1 ,0)
  cf.roc <- confusionMatrix(as.factor(pred), as.factor(y_test), positive = "1")
  
  # Using the sesp
  opt.spse <- optimal.cutpoints("X", "Y", tag.healthy = 0, method = "SpEqualSe", data = X)
  s.sesp <- summary(opt.spse)
  sesp <-  s.sesp$p.table$Global$SpEqualSe[[1]][1]
  pred <-  ifelse(predicted$pred > sesp, 1 ,0)
  cf.sesp <- confusionMatrix(as.factor(pred), as.factor(y_test), positive = "1")
  
  r <-  r + 1
  
  accs <- c( roc$auc, cf50$overall[1], cf.roc$overall[1], cf.sesp$overall[1], length(model$varNames))
  accuracies_testing <-  rbind(accuracies_testing, accs)
  
}

# Attribute column names
colnames(accuracies_testing) <-  c("AUC","Fifty","ROC01","SpEqualSE", "Number_antibody")


# Bind both the train an test dataset to analyze if they complete the 85% criteria
acc.grahps <-  rbind( data.frame(accuracies_testing, "Subset" = "Test" ),data.frame(accuracies_training, "Subset" = "Train"))

plt2 <- ggplot(acc.grahps, aes(x= Number_antibody, y= ROC01, color = Subset ))+ geom_line(size=1) +
  geom_hline(yintercept=0.85,
             linetype="dashed", 
             color = "black")+
  geom_hline(yintercept=0.50,
             linetype="dashed", 
             color = "blue")+
  theme_bw() +
  ylab("Accuracy") +
  xlab("Number of antibodies") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  geom_point(aes(x=57, y=0.89), colour="grey9") +
  annotate("text", x=57 , y= 0.91, label= "Optimal classifier",  colour="grey9") +
  scale_y_continuous(breaks=c(0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1), limits = c(0.4,1)) +
  scale_color_manual(values = c("#F79044FF", "#9512A1FF")) 

library(cowplot)

#plot_grid(plt1, plt2, labels=c("A", "B"), c(0.6,0.3), ncol = 1, nrow = 2)
pltf <- plot_grid(plt1, plt2, ncol = 2, nrow = 1)
pltf


title <- ggdraw() + 
  draw_label("HCs against ME/CFS with non-infectious or unknown disease trigger",
    x=0.30,
    hjust = 0, size =15
  ) +
  theme(
    
    plot.margin = margin(0, 0, 0, 0)
  )


fig1 <-  plot_grid(title, pltf,
                   ncol = 1,
                   rel_heights = c(0.1, 1))


fig1
saveRDS(fig1, "fig1(B).RDS")

accuracies_testing$id <-  seq(1:nrow(accuracies_testing))
write.csv(slearn.subsets[[71]]$varNames, "no_antibodies.csv")
read.csv("no_antibodies.csv")


