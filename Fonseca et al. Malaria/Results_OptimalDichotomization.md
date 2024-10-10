---
title: "ME/CFS Normality"
author: "Andre"
date: '2022-04-11'
output:
  html_document:
    theme: united
    toc_depht: 4
    toc_float:
      collapsed: false 
      smoth_scroll: false
    toc: true
    code_download: true
    code_folding: hide
---
45+8

```{r}
knitr::opts_chunk$set(echo = TRUE, warning = FALSE)
```

# Load libraries
```{r}
library(doParallel)   # Parallel processing
library(data.table)   # Work in data_table form
library(naniar)       # View missing data
library(ggplot2)      # Plot data
library(caret)        # Cross validationa anlysis
library(ResourceSelection) # Perform the Hosmel- Lemeshow
library(caTools)      # Split the data
library(pROC)         # Roc curve   
library(glmnet)       # For ridge regression
library(tidyr)        # Arrange the data
library(dplyr)        # For data cleaning
library(gplots)       # Heatmap
library(psych)        # Phi coefficient
library(MASS)
library(pROC)
library(mixtools)   # Gaussian mixture model
library(nortest)
library(moments)    #perform agostino and jarque-bera test
library(teigen)     #perform mm for t-distributions
library(AID)        # Variety of functions that facilitate the analyzes
library(ggplot2)
library(multimode)
library(readxl)
```

## Specify the number of cores
```{r cores, message=FALSE} 
numCores <- detectCores()
print(numCores)
registerDoParallel(numCores)  # use multicore, set to the number of our cores
```

## Set a seed
```{r}
set.seed(42)
```

# Load the data
```{r}
# Load the dataset
dataset <- read_excel("Data_CFS.xlsx")
```

```{r}
# Transpose the dataset
CFS <- t(dataset)
CFS <- as.data.frame(CFS)
head(CFS)

# Set the colnames as the first row
colnames(CFS) <- CFS[1,]

#R emove the Name, Sequence, AG876, B95-8, GD1, Cao, Raji, P3HR-1 rows
CFS <- CFS[! rownames(CFS) %in% c("Name", "Sequence", "AG876", "B95-8", "GD1", "Cao", "Raji", "P3HR-1") ,]

# Treat the colnames  ( omit the " 00_ " in the begining of the EBV virus)
colnames(CFS) <- substring(colnames(CFS), 4)

length(unique(colnames(CFS)))

# Obtain the indexes for the duplicated antigens
duplicated_index <- duplicated(colnames(CFS))

# Subset the data whih has and which has not duplicated values
duplicated_data <- CFS[,duplicated_index]
unique_data <- CFS[,!duplicated_index]

duplicated_ags <- unique(colnames(duplicated_data))

# Let's now get the duplicated antigens and add a number in front to 
# distinguish them. Bare in mind that the duplicated ags are already
# the repeatitions of the value in unique_data

different_data <- data.frame("Index" = c(1:nrow(CFS)))

for (ag in duplicated_ags) {
  # Obtain the dataset of the duplicated ags
  dataset <- data.frame(duplicated_data[,colnames(duplicated_data) == ag])
  n <- ncol(dataset)
    
  for( k in 1:n){
    # Obtain the names with the indexes on the end
    colnames(dataset)[k] <-  paste0(ag[k],"_", k)
  }
    # Bind everything
  different_data <-  cbind(different_data, dataset)
}

# Remove the index column
different_data <-  different_data[, !colnames(different_data) %in% "Index"]

# Bind the different_data and the unique_data together
dataset <-cbind(unique_data, different_data)
dim(dataset)

# To confirm that everything was eprforemd well
length(unique(colnames(dataset)))


# The IDS are only for the CFS individuals
# Therefore we need to devide in healthy and CFS
data.healthy <- CFS[ grepl( "HEAL" , rownames(dataset)), ]
data.cfs <- CFS[ grepl( "CFS" , rownames(dataset)), ]

# Load the dataset with the ids
IDS <- read.csv("Data_Epi_Data_CFS_Infection.csv", row.names = 1)

#Order both datasets
which(rownames(data.cfs) != rownames(IDS))

#Add na column for the heatly dataset
data.cfs <- cbind(IDS,data.cfs)

#Obtain the names of the absent variables in the healthy individuals
names(data.cfs)[1:3]
data.cfs <- data.frame("Status"= rep("cfs",nrow(data.cfs)), data.cfs)

# Create a matrix that wll be filled with NAs just to create the dull dataset
na.mat <- matrix(NA, 50,3)
na.mat <- data.frame(na.mat)
colnames(na.mat) <- colnames(data.cfs)[2:4]

data.healthy <- cbind(na.mat, data.healthy)
data.healthy <- data.frame("Status" = rep("healthy", nrow(data.healthy)), data.healthy)

#Bind both datasets
which(colnames(data.cfs)!= colnames(data.healthy))
CFS <- rbind(data.healthy, data.cfs) 
head(CFS)

# Devide the data into informative data and ...
info.data <- CFS[,colnames(CFS) %in% c("Status","infection.triggered.onset","age", "sex")]

#CFS values
CFS <-  CFS[,!colnames(CFS) %in% c("Status","infection.triggered.onset","age", "sex")]

# Replicate CFS dataset and transform data into numeric 
CFS.num <- CFS
n <-  ncol(CFS.num)
for (i in 1:n) {
 CFS.num[,i] <-  as.numeric(CFS.num[,i])
}

#Obtain the all data
CFS <-  data.table(info.data, CFS.num)
```


# Subset only the 50 healthy and 54 infection prior to disease individuals
```{r}
# Obtain the heathy individuals
m.healthy <- CFS[Status == "healthy",colnames(CFS.num), with = FALSE]

#Obtain the infected disease individuals
m.cfs <- CFS[infection.triggered.onset == "ja", colnames(CFS.num), with = FALSE]

#Add Status variable
m.healthy$Status <- rep("healthy", nrow(m.healthy))
m.cfs$Status <- rep("cfs", nrow(m.cfs)) 

# Bind both datasets ( healthy and cfs with infection)
#which(colnames(m.healthy)!= colnames(m.cfs))
data <- rbind(m.healthy, m.cfs)
data <-  data.frame(data)

# Dichotomize the response variable
data$Status <-  as.factor(data$Status)
```


# Data splitting
```{r}
dataset <-  data

spliting <- sample.split(dataset$Status, SplitRatio = 0.9)
train_set <-  subset(dataset, spliting == TRUE)
test_set <-  subset(dataset, spliting == FALSE)
training <-  train_set[, !colnames(train_set) %in% "Status"]

#write.csv(train_set, "train_set_MALI.csv")
#write.csv(test_set, "test_set_MALI.csv")
```


# Obtain the best thresholds
```{r}
target <-  train_set$Status

#l <-  ncol(training)
#best.pval <- foreach(antibody = training[1:l], name = #colnames(training)[1:l],.combine=rbind) %dopar%{
#  df <- data.frame(antibody, target)
#  colnames(df) <-  c("Ab", "Status")
#  k <- length(antibody)                                 
#  f.pval <- data.frame()
#  df <- df[order(df$Ab, decreasing = F), ]

  # Perform the loop
#  for (i in 1:k) {
#    tr <- df$Ab[i]
#    s <- ifelse( df$Ab >= tr , "healthy", "cfs")
#    sp <- table(s)
    
#    tab <- table(s, df$Status)
    
#    if( nrow(tab) < 2 | ncol(tab) <2 ){
#      p.val <- data.frame("p.value"= NA,"threshold" = tr, "Seronegative" = sp[1], #"Seropositive" = sp[2],
#                          "Seropositivity" = sp[2]/k, "antibody" = name )   
#    }else{
#      p.val <- data.frame("p.value"= chisq.test(tab)$p.value, "threshold" =tr, #"Seronegative" =sp[1],
#                          "Seropositive" = sp[2],"Seropositivity" = sp[2]/k, "antibody" #= name)
#    }
#    f.pval <- rbind(f.pval, p.val)
#   
#  }
#  b.pval <- f.pval[which.min(f.pval$p.value),]
#  return(b.pval)
#}

# Export as a CSV
#write.csv(best.pval, "results_pragmatic_best_pvals.csv")
```


# Load the parametric data p-value
```{r}
pval.data <- read.csv( "results_pragmatic_best_pvals.csv", row.names = 1)
```

# 1-Obtain statsitically significant antiboldies
```{r}
#Statistically signifcant antibodies without correction
treshold_list <-  list()

treshold <-  c(0.05,0.01, 0.001)

for(trs in treshold){
  sig.pval <-  subset(pval.data, pval.data$p.value <trs)
  sig.data <- as.data.frame(training[, colnames(training) %in% sig.pval$antibody])
  
  if(nrow(sig.data)==1){
    colnames(sig.data) <-  sig.pval$antibody
  }
  
  treshold_list <- c(treshold_list, list(sig.data))

  #Statistically signifcant antibodies with Bonferroni correction
  bf.pval <- pval.data
  bf.pval$p.value <- p.adjust(bf.pval$p.value, method = "bonferroni", n = length(bf.pval$p.value))
  sig.bf.pval <- subset(bf.pval, bf.pval$p.value <trs)
  sig.bf.data <- as.data.frame(training[, colnames(training) %in% sig.bf.pval$antibody])
    
  if(nrow(sig.bf.pval)==1){
    colnames(sig.bf.data ) <-  sig.pval$antibody
  }
  treshold_list <- c(treshold_list, list(sig.bf.data ))
  
  #Statistically signifcant antibodies with Benjamini-Yekutieli correction
  by.pval <- pval.data
  by.pval$p.value <- p.adjust(by.pval$p.value, method = "BY", n = length(by.pval$p.value))
  sig.by.pval <- subset(by.pval, by.pval$p.value < trs)
  sig.by.data <- as.data.frame(training[, colnames(training) %in% sig.by.pval$antibody])
    
  if(nrow(sig.by.pval)==1){
    colnames(sig.by.data) <-  sig.by.pval$antibody
  }
  
  treshold_list<- c(treshold_list, list(sig.by.data))
}

names(treshold_list) <-  c("0.05", "0.05 BF", "0.05 BY",
                           "0.01", "0.01 BF", "0.01 BY",
                           "0.001", "0.001 BF", "0.001 BY")


treshold_lists <- list()

for(t in treshold_list){
  if (dim(t)[2] > 1) {
    treshold_lists  <-  c(treshold_lists , list(t))
  }
  
}
```


# 2) Obtain the binary values
```{r}
binary_dataset <- list()

for(data in treshold_lists){
  df.abs <- data.frame("Status" = train_set$Status) # Create empty dataset
  info <-  pval.data[pval.data$antibody %in% colnames(data),] # Obtain the treshold for the abs
  which(info$antibody != colnames(data)) # See if the order differs
  
  for(j in 1:ncol(data)){
      tr  <- info$threshold[j]
      dt <-  data[,j] # Obtain data
      ab.dt <-  as.factor(ifelse(dt >= tr, 1, 0)) # Convert to factor
      df.abs <- cbind(df.abs, ab.dt)
      ab.dt <- c()
  }
  
    df.abs <- data.frame(df.abs[, !colnames(df.abs) %in% "Status"])
    colnames(df.abs) <- colnames(data)
    
    binary_dataset <- c(binary_dataset, list(df.abs))
  
}
```

We have obtained the 9 datasets:
- 0.05 non, bf, by
- 0.01 non, bf
- 0.001 non


# Obtain the test dataset
```{r}
binary_test <- list()

for(data in treshold_lists){
  df.abs <- data.frame("Status" = test_set$Status) # Create empty dataset
  info <-  pval.data[pval.data$antibody %in% colnames(data),] # Obtain the treshold for the abs
  
  # Obtain the test dataset
  test.abs <- test_set[colnames(test_set) %in% colnames(data)]
  which(info$antibody != colnames(test.abs)) # See if the order differs
  
  for(j in 1:ncol(test.abs)){
      tr  <- info$threshold[j]
      dt <-  test.abs[,j] # Obtain data
      ab.dt <-  as.factor(ifelse(dt >= tr, 1, 0)) # Convert to factor
      df.abs <- cbind(df.abs, ab.dt)
      ab.dt <- c()
  }
  
    df.abs <- data.frame(df.abs[, !colnames(df.abs) %in% "Status"])
    colnames(df.abs) <- colnames(test.abs)
    
    binary_test <- c(binary_test, list(df.abs))
  
}
```



# 3 - Remove highly correlated variables
```{r}
#Obtain the phi coefficients for each antibody pair
#phi.dataset <-  list()
#
#for (data in binary_dataset){
#  cor.matrix <- matrix(NA, ncol = ncol(data), nrow = ncol(data)) # Specify empty #matrix
#  
#  for(i in 1:ncol(data)){
#    for(j in 1:ncol(data)){
#      c <- data[,i]
#      r <- data[,j]
#      fd <- table(c,r)
#      p <- phi(fd)
#      cor.matrix[i,j] <-  p
#    }
#  }
#  cor.matrix <- as.data.frame(cor.matrix)
#  colnames(cor.matrix) <- colnames(data);   rownames(cor.matrix) <- #colnames(data)
#  phi.dataset <- c(phi.dataset, list(cor.matrix))
#}
#
## Save the data.frames outside R
#i <- 1
#for (data in phi.dataset){
#  write.csv(data, paste0("thephi_data_",i,".csv"))
#  i <- i+1
#}

# Load the data
phi.data <-  list()
for (i in 1:length(binary_dataset)) {
  loaded <- read.csv(paste0("C:/Users/andre/Desktop/Second paper/datasets/thephi_data_",i,".csv"), row.names = 1)
  phi.data <- c(phi.data, list(loaded))
}
```


```{r}
# Drop phi values above 0.7
binary_non_dataset <-  list()

for(i in 1:length(phi.data)){
  iteration <- phi.data[[i]]

  if (ncol(iteration) > 1) {
    iter <-  as.matrix(iteration)
    drop = findCorrelation(iter, cutoff = .7)
    drop = colnames(iteration)[drop]
    dt <- binary_dataset[[i]]
    
    iteration.cor <-  dt[,! colnames(dt) %in% drop]
    binary_non_dataset <- c(binary_non_dataset, list(iteration.cor))
  }
  
  else{
    dt <- binary_dataset[[i]]
    binary_non_dataset <- c(binary_non_dataset, list(dt))
  }
}

```


# Obtain the test dataset
```{r}
binary_non_test <-  list()

for(l in 1:length(binary_non_dataset)){
  t.names <- colnames(binary_non_dataset[[l]])
  t.data <- binary_test[[l]]
  binary_non_test <- c(binary_non_test, list(t.data[, colnames(t.data) %in% t.names]))

}
```


# Prepare the dataset for Predictive analysis
```{r}
# Bind cor and non cor analysi
all_data <- c(binary_dataset, binary_non_dataset)
all_test <- c(binary_test, binary_non_test)


# Add statust to all panels
all_datasets <-  list()
for (l in 1:length(all_data)) {
  all_datasets[[l]] <- data.frame(all_data[[l]], "Status"=train_set$Status)
  
}


data_names <- c("0.05", "0.05 BF", "0.05 BY",
                "0.01", "0.01 BF", 
                "0.001", 
                "0.05 non.cor", "0.05 BF non.cor", "0.05 BY non.cor",
                "0.01 non.cor", "0.01 BF non.cor",
                "0.001 non.cor")
```

# Predictive Analysis #
# Elastic net regression
```{r binary Elastic-net regression}
# Set the lambdas range
lambdas_to_try <- seq(0.001,1, length.out = 100)

objControl <- trainControl(method='repeatedcv', 
                           number=10,
                           repeats = 10,
                           returnResamp='all',
                           classProbs = TRUE,
                           savePredictions = "all",
                           verboseIter = FALSE,
                           allowParallel=TRUE,
                           summaryFunction = twoClassSummary,
                           preProc = c("center","scale"))

elastic <- function(x) {
  set.seed(42)
  train(Status~ .,
  data = x,
  method = "glmnet",
  trControl=objControl,
  family="binomial",
  metric = 'ROC',
  na.action = na.omit,
  tuneGrid = expand.grid(alpha = seq(0,1,0.1), lambda = lambdas_to_try))
}
  

el <- list()
el <- lapply(all_datasets, elastic)

for( i in 1:length(el)){
  vimp <- varImp(el[[i]])
  vimp <- vimp$importance
  len <-  dim(subset(vimp, vimp$Overall != 0))[1]
  print( paste( data_names[i],"-", len ,"antibodies" , ":",  max(el[[i]]$results$ROC)))
}

for( i in 1:length(el)){
  print(data_names[i])
  predictions <- predict(el[[i]], all_test[[i]])
  print(confusionMatrix(as.factor(test_set$Status), predictions))
}

```


# Random Forest
```{r random Forest}
objControl <- trainControl(method='repeatedcv', 
                           number=10,
                           repeats = 10,
                           returnResamp='all',
                           classProbs = TRUE,
                           search = "random",
                           savePredictions = "all",
                           verboseIter = FALSE ,
                           allowParallel=TRUE,
                           summaryFunction = twoClassSummary)

## Random search ##
my.rf <-  function(x){
  set.seed(42)
  train(Status ~ .,
  data = x ,
  method = 'rf',
  metric = 'ROC',
  trControl = objControl)
}


rf <- list()
rf <- lapply(all_datasets, my.rf)

for( i in 1:length(rf)){
  vimp <- varImp(rf[[i]])
  vimp <- vimp$importance
  len <-  dim(subset(vimp, vimp$Overall != 0))[1]
  print( paste(data_names[i],"-", len ,"antibodies" , ":",  max(rf[[i]]$results$ROC)))
}
  

for( i in 1:length(rf)){
  print(data_names[i])
  predictions <- predict(rf[[i]], all_test[[i]])
  print(confusionMatrix(as.factor(test_set$Status), predictions))
}

```


# XGBoost
```{r binary xgb}
xgbGrid <- expand.grid(nrounds = c(100,200),
                       max_depth = c(10, 15, 20, 25),
                       colsample_bytree = seq(0.5, 0.9, length.out = 5),
                       eta = 0.1,
                       gamma=0,
                       min_child_weight = 1,
                       subsample = 1)


my.xgb <-  function(x){
  set.seed(42)
  train(Status ~ .,
  data = x ,
  method = 'xgbTree',
  metric = 'ROC',
  trControl = objControl,
  tuneGrid = xgbGrid)
}


xgb <- list()
xgb <- lapply(all_datasets, my.xgb)


for( i in 1:length(xgb)){
  vimp <- varImp(xgb[[i]])
  vimp <- vimp$importance
  len <-  dim(subset(vimp, vimp$Overall != 0))[1]
  print( paste(data_names[i],"-", len ,"antibodies" , ":",  max(xgb[[i]]$results$ROC)))
}

for( i in 1:length(xgb)){
  print(data_names[i])
  predictions <- predict(xgb[[i]], all_test[[i]])
  print(confusionMatrix(as.factor(test_set$Status), predictions))
}
```


