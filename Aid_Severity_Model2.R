# This is Joint Table to Create the Need Attribute and Aid Attribute
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
library(igraph)
library(RColorBrewer)
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
#
# LOAD DATA
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

# LOAD LAT/LON COORDINATES (OF CENTROIDS) AND HLCIT CODES
hlcit <- read.csv(paste0(DIR,"master_hlcit.csv"))
colnames(hlcit) <- c("lon","lat","vdc_name","vname","hlcit_code")
hlcit$hlcit_code <- as.factor(hlcit$hlcit_code)
hlcit$vname <- as.character(hlcit$vname)
hlcit$vdc_name <- as.character(hlcit$vdc_name)
hlcit <- rm_space(hlcit,"hlcit_code")
hlcit$hlcit_code <- as.numeric(levels(hlcit$hlcit_code))[hlcit$hlcit_code]


# COLUMN NAMES FOR agency_relief.csv ARE:
aid_data <- read.csv(paste0(DIR,"agency_relief.csv"), sep=",")


# READ IN AID AND SEVERITY COMPLETE DATA TABLE
aid_sev <- readObj(file_name = paste0(DIR,"aid_and_severity.df"))

# READ MODELING TABLE
aid_sev_modeling <- readObj(file_name=paste0(DIR,"aid_sev_modeling.df"))


# SET NAs TO MEAN
aid_sev_modeling <- nas_to_mean(aid_sev_modeling)


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
# AGENCY-VDC AID NETWORK
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


# CHANGE FOMRAT TO CHARACTER FOR VDC AND AGENCY NAMES
aid_data$vdc <- trim(as.character(aid_data$vdc))
aid_data$impl_ag <- trim(as.character(aid_data$impl_ag))

# FILTER OUT THE EMPTY ENTRIES
aid_data <- aid_data[nchar(aid_data$vdc)>0 & nchar(aid_data$impl_ag)>0,]
aid_data <- rm_space(aid_data,"hlcit")
aid_data$hlcit <- as.numeric(levels(aid_data$hlcit))[aid_data$hlcit]
for (k in 1:dim(aid_data)[1]){
  aid_data$vdc[k] <- hlcit$vdc_name[which(hlcit$hlcit_code %in% aid_data$hlcit[k])[1]]
}

# SINCE aid_data$hlcit AND aid_data$vdc HAVE 239 ROWS OF SIMULTANEOUS NAs,
# WE DROP THESE ROWS
aid_data <- aid_data[!is.na(aid_data$hlcit),]


# SELECT UNIQUE AGENCIES AND TARGET VDC
ag <- unique(aid_data$impl_ag)
hlc <- unique(aid_data$hlcit)
all <- union(ag,hlc)

# EXTRAPOLATE LHCIT NUMBERS FROM LAT AND LON FOR VDCS FROM HLCIT_MASTER
hl <- vector()
xc <- vector()
yc <- vector()
for (k in 1:length(hlc)){
  hl[k] <- hlc[k]
  xc[k] <- hlcit$lat[which(hlcit$hlcit_code==hlc[k])[1]]
  yc[k] <- hlcit$lon[which(hlcit$hlcit_code==hlc[k])[1]]
}
koords<-cbind(xc,yc)

xa <- 20+5*runif(length(ag))
ya <- 80+12*runif(length(ag))
koords1 <-cbind(xa,ya)
koords2<- rbind(koords1,koords)

# DEFINE THE AGENCY-VDC AID NETWORK ADJACENCY MATRIX
aid_m <- matrix(0,nrow=length(all),ncol=length(all))
for (i in 1:length(ag)){
  for (j in 1:length(hlc)){
    aid_m[[i,length(ag)+j]] <- 
      dim(aid_data[aid_data$impl_ag==ag[i] & aid_data$hlcit==hlc[j],c(3,5)])[1]
  }
}

# BUILD THE AGENCY-VDC AID NETWORK
av <- graph.adjacency(aid_m,
                      mode = "directed",
                      weighted = TRUE)

# COLOR VERTICES REPRESENTING AGENCIES (GREEN) AND VDCs (BLUE) WHERE AID WAS SENT
V(av)$color <- rep("green",length(all))
for (k in 1:length(all)){
  if(is.element(all[k],hlc)){
    V(av)$color[k] <- "SkyBlue2"
  }  
}

for (k in 1:dim(aid_m)[1]){
  if(k-1<length(ag)){
    V(av)$size[k] <- 3
    V(av)$name[k] <- ag[k]
  } else {
    V(av)$size[k] <- 2
    V(av)$name[k] <- NA}
}

# PLOT THE AGENCY-VDC AID NETWORK
plot(av,
     layout = layout.fruchterman.reingold(av, 
                                          niter = 200),
     vertex.color = V(av)$color,
     vertex.size = V(av)$size,
     vertex.label = NA, 
     vertex.label.color = "black", 
     vertex.label.font = 1, 
     vertex.label.cex = 0.2, 
     edge.width = 0.2*sqrt(E(av)$weight),
     edge.arrow.size = 0.2,
     edge.curved = TRUE,
     edge.color = gray.colors(1),
     main = "Weighted Agency-VDC Aid Abstract Network (VDC Level)")
legend("topleft",
       c("Implementing Aid Agency","Aid Target VDC"),
       fill = c("green","SkyBlue2"),
       bty = "n")






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
# PROJECT EACH GRAPH COMPONENT WITH APPROPRIATE CONNECTIONS
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


# NOTE FOR THE AGENCY NETWORK PROJECTION, WE ARE COUNTING THE NUMBER OF
# INSTANCES WHEN TWO AGENCIES SUPPLIED AID TO THE SAME VDC
# WE CAN ALSO AUGMENT THIS MEASURE BY ACCOUNTING FOR THE NUMBER OF DIFFERENT 
# COMMON INSTANCES OF AID WITHIN A VDC, A MORE GRANULAR APPROACH
ag_m <- matrix(0,nrow=length(ag),ncol=length(ag))
for (i in 1:length(ag)){
  for (j in 1:length(ag)){
    common <- aid_m[c(i,j),(length(ag)+1):dim(aid_m)[1]]
    common[common>0] <- 1
    ag_m[[i,j]] <- sum(t(common[1,])*t(common[2,]))
  }
}

# REMOVE SELF LOOPS
for (k in 1:dim(ag_m)[1]){ag_m[[k,k]] <- 0}

# DEFINE AGENCY GRAPH
agg <- as.undirected(graph.adjacency(ag_m,weighted=TRUE))

# SET THE GRAPH COLOR, LABELS, AND NODE SIZE
V(agg)$color <- rep("green",length(ag))
V(agg)$name <- ag
V(agg)$size <- sqrt(degree(agg)-min(degree(agg)))

# PLOT AGENCY GRAPH AND FILTER
plot(agg,
     layout = layout.fruchterman.reingold(agg, 
                                          niter = 200),
     vertex.color = V(agg)$color,
     vertex.size = V(agg)$size,
     vertex.label = V(agg)$name, 
     vertex.label.color = "darkgreen", 
     vertex.label.font = 1, 
     vertex.label.cex = 0.75, 
     edge.width = 0.25*E(agg)$weight,
     edge.curved = TRUE,
     edge.color = gray.colors(1),
     main = "Aid Agency Association Network (Node Size = Sqrt(Degree)
     (Node Degree = number of agencies with shared VDC aid targets)")


# CHANGE THE NODE SIZE TO REFLECT NUBER OF SHARED VDCs USING SQRT
V(agg)$size <- sqrt(graph.strength(agg)-min(graph.strength(agg)))

# PLOT AGENCY GRAPH AND FILTER
plot(agg,
     layout = layout.fruchterman.reingold(agg, 
                                          niter = 200),
     vertex.color = V(agg)$color,
     vertex.size = V(agg)$size,
     vertex.label = V(agg)$name, 
     vertex.label.color = "darkgreen", 
     vertex.label.font = 1, 
     vertex.label.cex = 0.75, 
     edge.width = 0.05*E(agg)$weight,
     edge.curved = TRUE,
     edge.color = gray.colors(1),
     main = "Aid Agency Association Network (Node Size = Sqrt(Weighted Degree)
     (Node Weighted Degree = number of agencies with shared VDC aid targets, 
     weighted by the number of VDCs per agency)")


# CHANGE THE NODE SIZE TO REFLECT NUMBER OF SHARED VDCs USING LOG
V(agg)$size <- log(graph.strength(agg))

# PLOT AGENCY GRAPH AND FILTER
plot(agg,
     layout = layout.fruchterman.reingold(agg, 
                                          niter = 200),
     vertex.color = V(agg)$color,
     vertex.size = V(agg)$size,
     vertex.label = V(agg)$name, 
     vertex.label.color = "darkgreen", 
     vertex.label.font = 1, 
     vertex.label.cex = 0.75, 
     edge.width = 0.05*E(agg)$weight,
     edge.curved = TRUE,
     edge.color = gray.colors(1),
     main = "Aid Agency Association Network (Node Size = Log(Weighted Degree)
     (Node Weighted Degree = number of agencies with shared VDC aid targets, 
     weighted by the number of VDCs per agency)")


# HEAT MAP ACCORDING TO DEGREE
agg <- as.undirected(graph.adjacency(ag_m,weighted=TRUE))
V(agg)$color <- rep("green",length(ag))
V(agg)$name <- ag
V(agg)$size <- sqrt(degree(agg)-min(degree(agg)))

for (k in 1:length(degree(agg))){
  V(agg)$color[k] <- rev(heat.colors(1+max(degree(agg))))[degree(agg)[k]+1]
}

# PLOT AGENCY GRAPH AND FILTER
plot(agg,
     layout = layout.fruchterman.reingold(agg,
                                          niter = 20000),
     vertex.color = V(agg)$color,
     vertex.size = 5,
     vertex.label = V(agg)$name, 
     vertex.label.color = "darkgreen",
     vertex.label.font = 0.03, 
     vertex.label.cex = 0.5, 
     edge.width = 0.05*E(agg)$weight,
     edge.curved = TRUE,
     edge.color = gray.colors(1),
     main = "Aid Agency Association Network (Node Size = Log(Weighted Degree)
     (Node Weighted Degree = number of agencies with shared VDC aid targets, 
     weighted by the number of VDCs per agency)")



#
#
#
#
#
#
#
# ANALYSIS OF THE AID AGENCY NETWORK
#
#
#
#
#
#
#
#


# RANGE OF NUMBER OF DISTINCT AID INSTANCES FOR EACH AGENCY
summary(as.data.frame(table(aid_data$impl_ag))[,2])

# NOTE: THIS IS NOT THE SAME AS
# summary(graph.strength(av))
# SINCE BOTH AGENCIES AND VDCs ARE INCLUDED IN THIS

# RANGE OF NUMBER OF DISTINCT VDCs OF AID FOR EACH AGENCY
unique_aid <- unique(cbind.data.frame(aid_data$impl_ag,aid_data$hlcit))
colnames(unique_aid) <- c("impl_ag","hlcit")
summary(as.data.frame(table(unique_aid$impl_ag))[,2])

# NOTE: THIS IS NOT THE SAME AS
# summary(degree(av))
# SINCE BOTH AGENCIES AND VDCs ARE INCLUDED IN THIS

# PLOT RELIEF AGENCY WEIGHTED DEGREE DISTRIBUTION (DISTINCT TYPES OF AID)
plot(sort(as.data.frame(table(aid_data$impl_ag))[,2]),
     col = adjustcolor(rgb(1,0,1,1)),
     pch = 19,
     xlab = "Agency index",
     ylab = "Numer of Distinct Aid Activities",
     main = "Sorted Agencies by Number of Distinct Aid Activities")
par(new = T)
lines(x = c(0,length(ag)),y = rep(mean(as.data.frame(table(aid_data$impl_ag))[,2]),2), col ="black", lwd=4)
text(x = 25,y = 75,paste("MEAN =",mean(as.data.frame(table(aid_data$impl_ag))[,2])),col="black",cex=2.5)


histP1(as.data.frame(table(aid_data$impl_ag))[,2],
       breaks=100,
       col = adjustcolor(rgb(1,0,1,1)),
       xlab="Agency Network Aid Action Numbers", 
       main="Agency Network Number of Aid Actions Distribution
       (VDC Overlap Counts Dsitribution)")


# PLOT RELIEF AGENCY DEGREE DISTRIBUTION (DISTINCT VDCs)
plot(sort(as.data.frame(table(unique_aid$impl_ag))[,2]),
     col = adjustcolor(rgb(1,0,1,1)),
     pch = 19,
     xlab = "Agency index",
     ylab = "Numer of Distinct Aid Activities",
     main = "Sorted Agencies by Number of Distinct VDC")

hist(as.data.frame(table(unique_aid$impl_ag))[,2], breaks=100,
     col=adjustcolor(rgb(1,0,1,1)),
     xlab="Agency Network Aid Action Numbers", 
     main="Agency Network Number of Targeted VDCs Distribution
     (VDC Overlap Counts Distribution)")



# ANALYSIS OF THE AGENCY NETWORK ITSELF: OVERLAP OF AGENCY EFFORTS
agg <- as.undirected(graph.adjacency(ag_m,weighted=TRUE))

# SET THE GRAPH COLOR, LABELS, AND NODE SIZE
V(agg)$color <- rep("green",length(ag))
V(agg)$name <- ag
V(agg)$size <- sqrt(degree(agg)-min(degree(agg)))

# THIS IS THE NUMBER OF AGENCIES WITH COMMON VDC TARGETS AS A GIVEN AGENCY
summary(degree(agg))

# THIS IS THE WEIGHTED NUMBER OF THE ABOVE AGENCIES, SO THE NUMBER OF SHARED VDC
# TARGETS IS ACCOUNTED FOR BETWEEN EACH PAIR OF AGENCIES
summary(graph.strength(agg))


# PLOT THE NUMBER OF DISTINCT AGENCIES THAT SHARE TARGETS WITH A GIVEN AGENCY
plot(sort(degree(agg)),
     col = adjustcolor(rgb(1,0,1,1)),
     pch = 19,
     xlab = "Agency index",
     ylab = "Numer of Agencies",
     main = "Number of Agencies with Shared Target with an Agency (Sorted)")


histP1(degree(agg),
       breaks = 50,
       col = "green",
       xlab = "Agency Network Degree Values", 
       main = "Agency Network Degree Distribution
       (Distribution of the Number of Agencies with Common Targets as a Given Agency)")


# PLOT THE NUMBER OF DISTINCT AGENCIES THAT SHARE TARGETS WITH A GIVEN AGENCY
# WEIGHTED BY THE NUMBER OF SHARED VDC DISTRICT BETWEEN EACH PAIR OF AGENCIES
plot(sort(graph.strength(agg)),
     col = adjustcolor(rgb(1,0,1,1)),
     pch = 19,
     xlab = "Agency index",
     ylab = "Numer of Agencies",
     main = "Number of Agencies with Shared Target with an Agency (Sorted)
     (Weighted By The Number of Shared VDCs)")


histP1(graph.strength(agg), 
       breaks = 50,
       col = adjustcolor(rgb(1,0,1,1)),
       xlab = "Weighted Degree Values", 
       main = "Agency Network Weighted Degree Distribution
       (Weighted VDC Overlap Counts Dsitribution)")


# GRAPH DENSITY IS THE RATIO OF THE NUMBER OF EDGES AND THE NUMBER OF POSSIBLE EDGES
# TYPICALLY ON THE ORDER OF 1-10%
100*graph.density(agg)

# CLUSTERS ARE CONNECTED COMPONENTS, WE HAVE 4 in the UNFILTERED AGENCY-VDC NETWORK 
clusters(agg)$no

# SORTED CLUSTERS BY SIZE, NOTE THAT FILTRATIONS RESULT IN INCREASED NUMBER OF CLUSTERS AND A DROP IN CLUSTER SIZE
sort(clusters(agg)$csize,decreasing=TRUE)

# GLOBAL CLUSTERING COEFFICIENT (TRANSITIVITY) IS THE RATIO FO TRIANGLES AND CONNECTED TRIPLES
transitivity(agg)
cut75 <- quantile(as.vector(ag_m[ag_m>0]),0.75)
agg_f <- filter(cutoff = cut75,
                edge_matrix = ag_m,
                vertex_colors = V(agg)$color,
                vertex_names = ag,
                vertex_size = V(agg)$size)
agg_f <- as.undirected(agg_f)
transitivity(agg_f)

# RELATIVE MAXIMAL CLUSTER SIZE (AS % OF NUMBER OF NODES) 
max(clusters(agg)$csize)/vcount(agg)

# RELATIVE NUMBER OF ISOLATED NODES (AS % OF NUMBER OF NODES)  
sum(degree(agg)==0)/vcount(agg)

# PATH DISTRIBUTION: This shows the different lengths of shortest paths (geodesics) in our network. 
agg <- as.undirected(graph.adjacency(ag_m,weighted=TRUE))
V(agg)$color <- rep("green", length(ag))
agg <- giant_comp(graph = agg,
                  vertex_color = V(agg)$color,
                  vertex_names = ag)
sh <- shortest.paths(agg)
is.na(sh) <- sapply(sh,is.infinite)
sh[1:5,1:5]
paths <- na.omit(as.vector(sh))
length(paths)
summary(paths)
plot(sort(paths),
     xlab = "Path Index", 
     ylab = "Path Length", 
     main = "Paths (sorted by length)", 
     pch = 20,
     col = adjustcolor(rgb(1,0,1,1)))


histP1(paths,
       breaks = 15,
       col = adjustcolor(rgb(1,0,1,1)),
       xlab = "Path Length Values",
       main = "Path Length Distribution for g")



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
#
#
#
#
# VDC AID TARGET NETWORK
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
#
#
#
#
#


# ANALYSIS OF THE VDC AID TARGET NETWORK
u_hlcit <- as.character(unique(aid_data$hlcit))

# BUILD THE SHARED AGENCY ASSOCIATION NETWORK FOR THE VDCs
aid_m_hlc <- matrix(0, nrow = length(u_hlcit), ncol = length(u_hlcit))
for (i in 1:length(u_hlcit)){
  for (j in 1:length(u_hlcit)){
    common <- aid_m[1:length(ag),c(i,j)]
    common[common>0] <- 1
    aid_m_hlc[i,j] <- sum((common[,1])*(common[,2]))
  }
}


# REMOVE SELF LOOPS
for (k in 1:dim(aid_m_hlc)[1]){aid_m_hlc[[k,k]] <- 0}

# DEFINE AGENCY GRAPH
vgg <- as.undirected(graph.adjacency(aid_m_hlc,
                                     weighted = TRUE))

# SET THE GRAPH PROPERTIES
V(vgg)$color <- rep("SkyBlue2",length(u_hlcit))
V(vgg)$name <- u_hlcit
V(vgg)$size <- log(exp(1)+degree(vgg)/max(degree(vgg)))


# SAVE HLCIT DEGREE FOR MODELING PURPOSES
hlcit_degree <- cbind.data.frame(u_hlcit,as.numeric(degree(vgg)))
colnames(hlcit_degree) <- c("hlcit","hlcit_degree")
write.csv(hlcit_degree,file = paste0(DIR,"hlcit_degree.csv"))
writeObj(hlcit_degree,file = paste0(DIR,"hlcit_degree.df"))


# PLOT AGENCY GRAPH AND FILTER
plot(vgg,
     layout = layout.fruchterman.reingold(vgg,
                                          niter = 200),
     vertex.color = "SkyBlue2",
     vertex.size = 2,
     vertex.label = NA, 
     vertex.label.color = "black",
     vertex.label.font = 1, 
     vertex.label.cex = 0.5, 
     edge.width = 0.05*E(vgg)$weight,
     edge.curved = TRUE,
     edge.color = gray.colors(1),
     main = "VDC Aid Association Network
     (node degree = number of VDCs with same aid agency)")


# REMOVE ISOLATED
vgg <- drop_isolated(graph = vgg,
                     vertex_colors = V(vgg)$color,
                     vertex_names = V(vgg)$name)

# GET THE GIANT CONENCTED COMPONENT (TWO CLSUTERS ONLY)
vgg <- giant_comp(graph = vgg, 
                  vertex_colors = V(vgg)$color,
                  vertex_names = V(vgg)$name,
                  vertex_size = V(vgg)$size)

plot(vgg,
     layout = layout.fruchterman.reingold(vgg,
                                          niter = 200),
     vertex.color = "SkyBlue2",
     vertex.size = 2,
     vertex.label = NA, 
     vertex.label.color = "black",
     vertex.label.font = 1, 
     vertex.label.cex = 0.5, 
     edge.width = 0.03*E(vgg)$weight,
     edge.curved = TRUE,
     edge.color = gray.colors(1),
     main = "VDC Aid Association Network Giant Component
     (node degree = number of VDCs with same aid agency)")


# FILTER BY EDGE WEIGHT LEVELS, NOTE THAT LOWER CUTOFF VALUES
vgg <- as.undirected(graph.adjacency(aid_m_hlc,weighted=TRUE))
V(vgg)$color <- rep("SkyBlue2",length(u_hlcit))
V(vgg)$name <- u_hlcit
V(vgg)$size <- log(3+degree(vgg)/max(degree(vgg)))
cut15 <- quantile(as.vector(aid_vdc[aid_vdc>0]),0.15)
vgg_f <- filter(cutoff = cut15,
                edge_matrix = aid_m_hlc,
                vertex_colors = V(vgg)$color,
                vertex_names = V(vgg)$name,
                vertex_size = V(vgg)$size)
vgg_f <- as.undirected(vgg_f)
vgg_f <- giant_comp(graph = vgg_f,
                    vertex_colors = V(vgg_f)$color,
                    vertex_names = V(vgg_f)$name,
                    vertex_size = V(vgg_f)$size)
plot(as.undirected(vgg_f),
     layout=layout.fruchterman.reingold(vgg_f,
                                        niter = 200),
     vertex.color = V(vgg_f)$color,
     vertex.size = 2,
     vertex.label = NA, 
     vertex.label.color = "black", 
     vertex.label.font = 1, 
     vertex.label.cex = 1, 
     edge.width = 0.05*(E(vgg_f)$weight),
     edge.curved = TRUE,
     edge.color = gray.colors(1),
     main = "VDC Aid Association Network With 15% Threshold Filtration
     (node degree = number of VDCs with same aid agency)")

# THE NEXT CUT IS TOO BIG OF A JUMP, WHICH WILL BE EVIDENT IN THE DEGREE DISTRIBUTION
vgg <- as.undirected(graph.adjacency(aid_m_hlc,weighted=TRUE))
V(vgg)$color <- rep("SkyBlue2",length(u_hlcit))
V(vgg)$name <- u_hlcit
V(vgg)$size <- log(3+degree(vgg)/max(degree(vgg)))
cut85 <- quantile(as.vector(aid_vdc[aid_vdc>0]),0.85)
vgg_f <- filter(cutoff = cut85,
                edge_matrix = aid_m_hlc,
                vertex_colors = V(vgg)$color,
                vertex_names = V(vgg)$name,
                vertex_size = V(vgg)$size)
vgg_f <- as.undirected(vgg_f)
vgg_f <- giant_comp(graph = vgg_f,
                    vertex_colors = V(vgg_f)$color,
                    vertex_names = V(vgg_f)$name,
                    vertex_size = V(vgg_f)$size)
plot(as.undirected(vgg_f),
     layout=layout.fruchterman.reingold(vgg_f,
                                        niter = 200),
     vertex.color = V(vgg_f)$color,
     vertex.size = 2,
     vertex.label = NA, 
     vertex.label.color = "black", 
     vertex.label.font = 1, 
     vertex.label.cex = 1, 
     edge.width = 0.05*(E(vgg_f)$weight),
     edge.curved = TRUE,
     edge.color=gray.colors(1),
     main = "VDC Aid Association Network With 85% Threshold Filtration
     (node degree = number of VDCs with same aid agency)")


# ANALYSIS OF THE TARGET VDC NETWORK ITSELF:
vgg <- as.undirected(graph.adjacency(aid_m_hlc,
                                     weighted = TRUE))
V(vgg)$color <- rep("SkyBlue2",length(u_hlcit))
V(vgg)$name <- u_hlcit
V(vgg)$size <- log(3+degree(vgg)/max(degree(vgg)))

# THIS IS THE NUMBER OF VDCS A GIVEN VDC HAS AN AGENCY IN COMMON
summary(degree(vgg))

# THIS IS THE WEIGHTED NUMBER OF THE ABOVE INSTANCES, SO THE NUMBER OF SHARED
# AGENCIES IS ACCOUNTED FOR BETWEEN EACH PAIR OF VDCs
summary(graph.strength(vgg))

# PLOT THE NUMBER OF DISTINCT AGENCIES THAT SHARE TARGETS WITH A GIVEN AGENCY
plot(sort(degree(vgg)),
     col = adjustcolor(rgb(1,0,1,1)),
     pch = 19,
     xlab = "VDC index",
     ylab = "Numer of VDCs",
     main = "Number of VDCs with Shared Agency with a Given VDC (Sorted)")

histP2(degree(vgg),
       breaks = 130,
       col = adjustcolor(rgb(1,0,1,1)),
       xlab = "VDC Aid Association Network Degree Values", 
       main = "VDC Aid Association Network Degree Distribution")


# PLOT THE NUMBER OF DISTINCT AGENCIES THAT SHARE TARGETS WITH A GIVEN AGENCY
# WEIGHTED BY THE NUMBER OF SHARED VDC DISTRICT BETWEEN EACH PAIR OF AGENCIES
plot(sort(graph.strength(vgg)),
     col = adjustcolor(rgb(1,0,1,1)),
     pch = 19,
     xlab = "Agency index",
     ylab = "Numer of Agencies",
     main = "Number of Agencies with Shared Target with an Agency (Sorted)
     (Weighted By The Number of Shared VDCs)")


histP2(graph.strength(vgg), 
       breaks = 150,
       col = "SkyBlue2",
       xlab = "Weighted Degree Values", 
       main = "Agency Network Weighted Degree Distribution
       (Weighted VDC Overlap Counts Distribution)")


# GRAPH DENSITY IS THE RATIO OF THE NUMBER OF EDGES AND THE NUMBER OF POSSIBLE EDGES
# TYPICALLY ON THE ORDER OF 1-10%
100*graph.density(vgg)

# CLUSTERS ARE CONNECTED COMPONENTS
clusters(vgg)$no

# SORTED CLUSTERS BY SIZE, NOTE THAT FILTRATIONS RESULT IN INCREASED NUMBER OF CLUSTERS AND A DROP IN CLUSTER SIZE
sort(clusters(vgg)$csize,decreasing=TRUE)



# GLOBAL CLUSTERING COEFFICIENT (TRANSITIVITY) IS THE RATIO FO TRIANGLES AND CONNECTED TRIPLES
vgg <- as.undirected(graph.adjacency(aid_m_hlc,weighted=TRUE))
V(vgg)$color <- rep("SkyBlue2",length(u_hlcit))
V(vgg)$name <- u_hlcit
V(vgg)$size <- log(3+degree(vgg)/max(degree(vgg)))
cut15 <- quantile(as.vector(aid_vdc[aid_vdc>0]),0.15)
vgg_f <- filter(cutoff = cut15,
                edge_matrix = aid_m_hlc,
                vertex_colors = V(vgg)$color,
                vertex_names = V(vgg)$name,
                vertex_size = V(vgg)$size)
vgg_f <- as.undirected(vgg_f)
transitivity(vgg_f)

# RELATIVE MAXIMAL CLUSTER SIZE (AS % OF NUMBER OF NODES) 
max(clusters(vgg)$csize)/vcount(vgg)

# RELATIVE NUMBER OF ISOLATED NODES (AS % OF NUMBER OF NODES)  
sum(degree(vgg)==0)/vcount(vgg)

# PATH DISTRIBUTION: This shows the different lengths of shortest paths (geodesics) in our network. 
vgg <- as.undirected(graph.adjacency(aid_m_hlc,weighted=TRUE))
V(vgg)$color <- rep("SkyBlue2",length(u_hlcit))
V(vgg)$name <- u_hlcit
V(vgg)$size <- log(3+degree(vgg)/max(degree(vgg)))
cut15 <- quantile(as.vector(aid_vdc[aid_vdc>0]),0.15)
vgg_f <- filter(cutoff = cut15,
                edge_matrix = aid_m_hlc,
                vertex_colors = V(vgg)$color,
                vertex_names = V(vgg)$name,
                vertex_size = V(vgg)$size)
sh<-shortest.paths(vgg)
is.na(sh)<-sapply(sh,is.infinite)
sh[1:5,1:5]
paths<-na.omit(as.vector(sh))
length(paths)
summary(paths)
plot(sort(paths),
     xlab="Path Index", 
     ylab="Path Length", 
     main="Paths (sorted by length)", 
     pch=20,
     col=adjustcolor(rgb(1,0,1/2,1)))
hist(paths,
     breaks=15,
     col=adjustcolor(rgb(1,0,1,1)),
     xlab="Path Length Values",
     main="Path Length Distribution for g")



#
#
#
#
# HEAT MAPS ACCORDING TO DEGREE DISTRIBUTION
# AND WEIGHTED DEGREE DISTRIBUTION
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
# HEAT MAP FOR VDCs ACCORDING TO DEGREE
#
#
#
#
#

xa1 <- 20+5*runif(length(ag))
ya1 <- 80+12*runif(length(ag))
koords11 <-cbind(xa1,ya1)
koords21<- rbind(koords11,koords)

# BUILD THE AGENCY-VDC AID NETWORK
av <- graph.adjacency(aid_m,
                      mode = "directed",
                      weighted = TRUE)
vgg <- as.undirected(graph.adjacency(aid_vdc,weighted=TRUE))


# COLOR VERTICES REPRESENTING AGENCIES (GREEN) AND VDCs (BLUE) WHERE AID WAS SENT
V(av)$color <- rep("green",length(all))
for (k in 1:length(all)){
  if(is.element(all[k],vd)){
    V(av)$color[k] <- rev(heat.colors(1+as.integer(max(degree(vgg)))))[as.integer(degree(vgg)[which(vd %in% all[k])])+1]
  }  
}

# PLOT THE AGENCY-VDC AID NETWORK
for (k in 1:dim(aid_m)[1]){
  if(k-1<length(ag)){
    V(av)$size[k] <- 3
    V(av)$name[k] <- ag[k]
  } else {
    V(av)$size[k] <- 1
    V(av)$name[k] <- NA}
}

plot(av,
     layout = koords21,
     vertex.color = V(av)$color,
     vertex.size = V(av)$size,
     vertex.label = V(av)$name, 
     vertex.label.color = "darkgreen", 
     vertex.label.font = 1, 
     vertex.label.cex = 0.75, 
     edge.width = 0.05*sqrt(E(av)$weight),
     edge.arrow.size = 0.2,
     edge.curved = TRUE,
     edge.color = gray.colors(1),
     main = "Nepal Agency-VDC Aid Relief Geo-Network")
legend("topright",
       c("Implementing Aid Agency","VDC With Geo-Coordinates"),
       fill = c("green","darkorange"),
       bty = "n")



xa1 <- 20+5*runif(length(ag))
ya1 <- 82+10*runif(length(ag))
koords11 <-cbind(xa1,ya1)
koords21<- rbind(koords11,koords)
plot(av,
     layout = koords21,
     vertex.color = V(av)$color,
     vertex.size = V(av)$size,
     vertex.label = V(av)$name, 
     vertex.label.color = "darkgreen", 
     vertex.label.font = 1, 
     vertex.label.cex = 0.75, 
     edge.width = 0.05*sqrt(E(av)$weight),
     edge.arrow.size = 0.2,
     edge.curved = TRUE,
     edge.color = gray.colors(1),
     main = "Nepal Agency-VDC Aid Relief Geo-Network")
legend("topleft",
       c("Implementing Aid Agencies","VDCs With Geo-Coordinates"),
       fill = c("green","darkorange"),
       bty = "n")






V(av)$color <- rep("green",length(all))
for (k in 1:length(all)){
  if(is.element(all[k],vd)){
    V(av)$color[k] <- "SkyBlue2"
  }  
}


plot(av,
     layout = koords21,
     vertex.color = V(av)$color,
     vertex.size = V(av)$size,
     vertex.label = V(av)$name, 
     vertex.label.color = "darkgreen", 
     vertex.label.font = 1, 
     vertex.label.cex = 0.75, 
     edge.width = 0.05*sqrt(E(av)$weight),
     edge.arrow.size = 0.2,
     edge.curved = TRUE,
     edge.color = gray.colors(1),
     main = "Nepal Agency-VDC Aid Relief Geo-Network")
legend("topleft",
       c("Implementing Aid Agencies","VDCs With Geo-Coordinates"),
       fill = c("green","SkyBlue2"),
       bty = "n")









































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





