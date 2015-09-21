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
library(leaps)
library(car)
library(relaimpo)
library(foreign)
library(reshape2)
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

aid <- scale(aid_sev_modeling)
# CREATE INITIAL FEATURE RANKING
feature_weights <- weightsRF(df = aid_sev_modeling,
                             target = target,
                             ntree = ntree,
                             mtry = mtry)
initial_features <- feature_weights$features

# print(feature_weights)
#        features importance
# 1      severity 193.406731
# 2        hazard 111.277955
# 3       poverty  77.736398
# 4 vulnerability  67.782081
# 5       housing  45.624705
# 6      exposure  -5.480277


# CORRELATION ANALYSIS
dat <- melt(round(cor(aid_sev_modeling),2))
dat$Var1 <- factor(dat$Var1,levels=colnames(aid_sev_modeling))
dat$Var2 <- factor(dat$Var2,levels=rev(colnames(aid_sev_modeling)))
ggplot(data = dat,
       aes(x = Var1, y = Var2)) + 
  geom_tile(aes(fill = value),colour = "white") + 
  geom_text(aes(label = sprintf("%1.2f",value)),vjust = 1) + 
  scale_fill_gradient(low = "white",high = "darkorange") + 
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  ggtitle("Correaltion Analysis of Severity Variables") + 
  theme(plot.title = element_text(face="bold"))

# HERE ARE THE VALUES
# round(cor(aid_sev_modeling),2)
# hazard exposure housing poverty vulnerability severity degree
# hazard          1.00     0.07    0.17   -0.24         -0.12     0.72  -0.06
# exposure        0.07     1.00    0.12   -0.12         -0.04     0.52  -0.05
# housing         0.17     0.12    1.00   -0.31          0.38     0.37  -0.05
# poverty        -0.24    -0.12   -0.31    1.00          0.76    -0.23   0.09
# vulnerability  -0.12    -0.04    0.38    0.76          1.00     0.03   0.06
# severity        0.72     0.52    0.37   -0.23          0.03     1.00  -0.07
# degree         -0.06    -0.05   -0.05    0.09          0.06    -0.07   1.00



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

# SUMMARY OF THE REGRESSION STATS
print(summary(fit))

# Call:
#   lm(formula = degree ~ ., data = train)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -154.61  -84.95  -21.24   75.12  299.13 
# 
# Coefficients: (1 not defined because of singularities)
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept)    139.432     32.243   4.324 1.82e-05 ***
#   hazard          -9.503      7.142  -1.330    0.184    
# exposure       -14.314     15.234  -0.940    0.348    
# housing         -1.395      4.333  -0.322    0.748    
# poverty          3.035      2.755   1.102    0.271    
# vulnerability       NA         NA      NA       NA    
# severity        11.409     24.678   0.462    0.644    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 102.3 on 542 degrees of freedom
# Multiple R-squared:  0.0134,	Adjusted R-squared:  0.004297 
# F-statistic: 1.472 on 5 and 542 DF,  p-value: 0.1971

# LOOK AT THE RESIDUALS, AUTOCORRELATION ANALYSIS
acf(fit$residuals)

# STANDARD RESIDUALS
plot(rstandard(fit))






#
#
#
#
#
#
#
# FURTHER INVESTIGATIONS
#
#
#
#
#
#


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
fit <- lm(degree~.,data=train)
step <- stepAIC(fit, direction="both")
step$anova # display results



# ROBUST REGRESSION ANALYSIS
summary(ols <- lm(degree ~ ., data = train))
opar <- par(mfrow = c(2,2), oma = c(0, 0, 1.1, 0))
plot(ols, las = 1)
par(opar)
d1 <- cooks.distance(ols)
r <- stdres(ols)
a <- cbind(train, d1, r)
a[d1 > 4/dim(train)[1], ]
rabs <- abs(r)
a <- cbind(train, d1, r, rabs)
asorted <- a[order(-rabs), ]
summary(rr.hubber <- rlm(degree ~ hazard+exposure+housing+poverty+severity, data = train))
hweights <- data.frame(severity = train$severity,degree = train$degree,resid = rr.huber$resid, weight = rr.huber$w)
hweights2 <- hweights[order(rr.huber$w), ]
hweights2[1:15, ]
rr.bisquare <- rlm(degree ~ hazard+exposure+housing+poverty+severity, data=train, psi = psi.bisquare)
summary(rr.bisquare)


#       severity degree    resid    weight
# 681  1.1760796    447 322.5770 0.4987031
# 579  0.7316551    443 309.9230 0.5190625
# 846  0.5702592    426 289.7085 0.5552798
# 562  0.2008899    410 267.1058 0.6022627
# 1830 1.1063028    366 240.2861 0.6694975
# 823  1.1834204    360 236.1161 0.6813207
# 2104 0.5773126    371 235.0135 0.6845122
# 627  0.5353111    362 225.4604 0.7135142
# 1815 0.8919543    352 223.0202 0.7213219
# 1331 0.7076701    353 220.8907 0.7282704
# 1744 0.7791218    351 218.9036 0.7348917
# 1141 0.7853564    350 218.4414 0.7364438
# 1410 0.7289234    338 205.0274 0.7846284
# 729  1.0206061    332 204.5952 0.7862914
# 1417 0.7365237    330 197.0624 0.8163440





