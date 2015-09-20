# This is a RF Forest Exploring the Predictors to the Variable Degree
# author: Georgi D. Gospodinov
# date: "September 20, 2015"
# 
# Data Sources:
#
# https://data.hdx.rwlabs.org/dataset/scnepal-agency-data
# master_hlcit.csv
# 
# Relevant materials and data can be found at:
# 
# https://www.dropbox.com/sh/tb9854hzcof7x23/AACEDTGk8EmYQ6r4ukSFLBspa?dl = 0
#
# in the folder /Displacement Network Model/
#
# 
#
# LOAD PACKAGES
library(plyr)
library(dplyr)
library(igraph)
library(RColorBrewer)
library(randomForest)
library(readr)
library(xgboost)
library(quantreg)
library(caret)
library(Hmisc)
library(ROCR)
library(tidyr)
library(scales)
library(rpart)
library(pROC)
library(bootstrap)
library(MASS)
#
#
#


# SET FILE SOURCE PATH
DIR <- "/Users/ggospodinov/Desktop/UN_OCHA_project/data/"


#
#
# AND THIS PART: SET UP RF PARAMETERS
ntree <- 501
mtry <- 2
target <- "degree"


#
#
# DEFINE FUNCTIONS
#
#
#


# WRITE OBJECT FUNCTION
writeObj <- function(obj, file_name) {
  save(obj, file=file_name)
  return(file_name)
}


# READ OBJECT FUNCTION
readObj <- function(file_name) {
  obj_name <- load(file_name)
  obj <- get(obj_name)
  return(obj)
}


# PERFORMANCE STATS FUNCTION
aprf <- function(tp, fp, fn, tn) {
  acc <- sum(tp,tn)/sum(tp,tn,fp,fn)
  prec <- tp/sum(tp,fp)
  recall <- tp/sum(tp,fn)
  fscore <- 2*prec*recall/sum(prec,recall)
  res <- data.frame(accuracy=acc, precision=prec, recall=recall, fscore=fscore)
  return(res)
}


# AGGREGATE STATS FUNCTION
aggStats <- function(df,
                     actual_col_name,
                     pred_col_name,
                     target_class) {
  
  group_columns <- c(actual_col_name, pred_col_name)
  dots <- lapply(group_columns, as.symbol)
  grouped <- df
  #t <- grouped %>% dplyr::group_by_(.dots=dots) %>% dplyr::summarize(count=n())
  grouped_df <- dplyr::group_by_(.data = grouped, .dots=dots)
  t <- dplyr::summarize(.data = grouped_df, count=n())
  
  tp <- t$count[t[[actual_col_name]]==target_class & t[[pred_col_name]]==target_class]
  fp <- t$count[t[[actual_col_name]]!=target_class & t[[pred_col_name]]==target_class]
  fn <- t$count[t[[actual_col_name]]==target_class & t[[pred_col_name]]!=target_class]
  tn <- t$count[t[[actual_col_name]]!=target_class & t[[pred_col_name]]!=target_class]
  
  return(aprf(tp, fp, fn, tn))
}


# RF TRAIN FUNCTION
rf.train <- function(df, target, ntree, mtry, subset)
{
  df2 <- df[,which(names(df) %in% c(target, subset))]
  fo <- paste0(target,"~.")
  rf.fit <- randomForest(as.formula(fo), 
                         data = df2,
                         ntree = ntree,
                         mtry = mtry,
                         na.action = na.omit)
  return(rf.fit)
}


# RF PERFORMANCE EVALUATION AND STATS FUNCTION
rf.performance <- function(df, rf.fit, target, subset)
{
  df2 <- df[,which(names(df) %in% subset)]
  yhat.train <- predict(rf.fit, df2, type="response")
  df.final1 <- data.frame(df[,which(names(df) %in% target)], yhat.train)
  df.final1 <- df.final1[complete.cases(df.final1),]
  names(df.final1) <- c("actual","predicted")
  actual_col_name <- "actual"
  pred_col_name <- "predicted"
  target_class <- 1
  return(aggStats(df.final1,actual_col_name, pred_col_name, target_class))
}


# FUNCTION TO SET NAs TO MEAN
nas_to_mean <- function(df){
  for (i in 1:(dim(df)[2])){
    df[is.na(df[,i]),i] <- mean(df[!is.na(df[,i]),i])
  }
  return(df)
}


# FUNCTION FOR VARIABLE SELECTION USING RANDOM FOREST
weightsRF <- function(df,target,ntree,mtry) {
  fo <- as.formula(paste0(target,"~ ."))
  rfmodel <- randomForest(fo, df,maxnodes=5,ntree=ntree,na.action=na.omit,mtry=mtry,nodesize=5, importance = TRUE)
  df.rfImportance <- cbind.data.frame(variable=names(rfmodel$importance[,1]),importance=rfmodel$importance[,1])
  df.rfImportance <- df.rfImportance[order(-df.rfImportance[,2]),]
  weights_rf <- cbind.data.frame(df.rfImportance[,1],df.rfImportance[,2])
  colnames(weights_rf) <- c("features","importance")
  weights_rf$features <- as.character(weights_rf$features)
  return(weights_rf)
}


#
#
#
#
#
#
#
#
# READ IN THE DATA AND EXCLUDE INDIVIDUAL_ID VAR
#
#
#
#
#
#
#
#

#
#
#
#
#
#
# READ IN THE DATA
#
#
#
#
#
#


# READ IN AID AND SEVERITY COMPLETE DATA TABLE
aid_sev <- readObj(file_name = paste0(DIR,"aid_and_severity.df"))

# READ MODELING TABLE
aid_sev_modeling <- readObj(file_name=paste0(DIR,"aid_sev_modeling.df"))


# SET NAs TO MEAN
aid_sev_modeling <- nas_to_mean(aid_sev_modeling)


# CREATE INITIAL FEATURE RANKING
feature_weights <- weightsRF(df = aid_sev_modeling,
                             target = target,
                             ntree = ntree,
                             mtry = mtry)
initial_features <- feature_weights$features

# print(feature_weights)
#        features importance
# 1      severity  130.47852
# 2 vulnerability   83.61868
# 3        hazard   79.03830
# 4       poverty   77.74036
# 5       housing   63.98707
# 6      exposure   10.19198


# FEATURES
var_list <- c("hazard","exposure","housing","poverty","vulnerability","severity") 


# SPLIT DATA IN 70-30 TRAIN/TEST
train.index <- createDataPartition(aid_sev_modeling$degree, p = .7, list = FALSE)
train <- aid_sev_modeling[ train.index,]
test  <- aid_sev_modeling[-train.index,]


# RUN THE LINEAR REGRESSION PREDICTOR WITH THE SELECTED FEATURES
fit <- lm(degree ~ ., data=train)

# GENERATE PLOTS
layout(matrix(c(1,2,3,4),2,2)) # optional 4 graphs/page 
plot(fit)

# Assessing R2 shrinkage using 10-Fold Cross-Validation 
theta.fit <- function(x,degree){lsfit(x,degree)}
theta.predict <- function(fit,x){cbind(1,x)%*%fit$coef} 

# matrix of predictors
X <- as.matrix(train[,var_list])
# vector of predicted values
y <- as.matrix(train[c("degree")]) 

results <- crossval(X,y,theta.fit,theta.predict,ngroup=10)
cor(y, fit$fitted.values)**2 # raw R2 
cor(y,results$cv.fit)**2 # cross-validated R2

# STEPWISE REGRESSION
# Stepwise Regression
library(MASS)
fit <- lm(degree~.,data=train)
step <- stepAIC(fit, direction="both")
step$anova # display results




# All Subsets Regression
library(leaps)
attach(train)
leaps<-regsubsets(degree~.,data=train,nbest=10)
# view results 
summary(leaps)
# plot a table of models showing variables in each model.
# models are ordered by the selection statistic.
plot(leaps,scale="r2")
# plot statistic by subset size 
library(car)
subsets(leaps, statistic="rsq")


library(relaimpo)
calc.relimp(fit,type=c("lmg","last","first","pratt"),
            rela=TRUE)

# Bootstrap Measures of Relative Importance (1000 samples) 
boot <- boot.relimp(fit, b = 1000, type = c("lmg", 
                                            "last", "first", "pratt"), rank = TRUE, 
                    diff = TRUE, rela = TRUE)
booteval.relimp(boot) # print result
plot(booteval.relimp(boot,sort=TRUE)) # plot result



