# This is the Modeling TAblef or Need and Aid Analysis
# author: Georgi D. Gospodinov
# date: "September 27, 2015"
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
library(mgcv)
library(stats)
library(cluster)
library(fpc)
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
library(ggplot2)
library(plotrix)
#
#
#


# SET FILE SOURCE PATH
DIR <- "/Users/ggospodinov/Desktop/UN_OCHA_project/data/"


#
#
# DEFINE FUNCTIONS
#
#
#


# DEFINE HISTOGRAM FORMATTING FUNCTIONS
histP <- function(x,breaks, ...) {
  H <- hist(x, plot = FALSE, breaks=breaks)
  H$density <- with(H, 100 * density* diff(breaks)[1])
  labs <- paste(round(H$density), "%", sep="")
  plot(H, freq = FALSE, labels = labs, ylim=c(0, 1.08*max(H$density)),...)
}

histP1 <- function(x,breaks, ...) {
  H <- hist(x, plot = FALSE, breaks=breaks)
  H$density <- with(H, 100 * density* diff(breaks)[1])
  labs <- ifelse(round(H$density)>0,paste(round(H$density), "%", sep=""),NA)
  plot(H, freq = FALSE, labels = labs, ylim=c(0, 1.08*max(H$density)),...)
}

histP2 <- function(x,breaks, ...) {
  H <- hist(x, plot = FALSE, breaks=breaks)
  H$density <- with(H, 100 * density* diff(breaks)[1])
  labs <- ifelse(round(H$density)>5,paste(round(H$density), "%", sep=""),NA)
  plot(H, freq = FALSE, labels = labs, ylim=c(0, 1.08*max(H$density)),...)
}


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


# FUNCITON TO REMOVE ALL PUNCTUATION CHARACTERS FROM LEVEL NAMES OF VARIABLES
rm_punct <- function(df){
  for (i in 1:dim(df)[2]){
    if (class(df[,i])=="factor" || class(df[,i])=="character" ){
      level_names <- unique(levels(df[,i]))
      df[,i] <- mapvalues(df[,i], from=level_names,to=gsub("[[:punct:]]","",level_names))
    }
  }
  return(df)
}


# FUNCITON THAT DROPS LOOPS
drop_loops <- function(graph, vertex_colors, vertex_names){
  g <- simplify(graph,remove.loops=TRUE)
  g_f <- delete.vertices(g,V(g)[degree(g)==0])
  v_g_f <- setdiff(V(g),V(g)[degree(g)==0])
  # filter names and color
  V(g_f)$name <- vertex_names[v_g_f]
  V(g_f)$color <- vertex_colors[v_g_f]
  return(g_f)
}


# FUNCITON TO REMOVE ALL SPACES FROM LEVEL NAMES OF VARIABLES
rm_space <- function(df){
  for (i in 1:dim(df)[2]){
    if (class(df[,i])=="factor" || class(df[,i])=="character" ){
      level_names <- unique(levels(df[,i]))
      df[,i] <- mapvalues(df[,i], from=level_names,to=gsub("[[:space:]]","",level_names))
    }
  }
  return(df)
}

# FUNCTION TO SET NAs TO MEAN
nas_to_mean <- function(df){
  for (i in 1:(dim(df)[2])){
    df[is.na(df[,i]),i] <- mean(df[!is.na(df[,i]),i])
  }
  return(df)
}

# FUNCTION THAT TRIMS LEADING WHITESPACE
trim.leading <- function (x)  sub("^\\s+", "", x)


# FUNCTION THAT TRIMS TRAILING WHITESPACE
trim.trailing <- function (x) sub("\\s+$", "", x)


# FUNCTION THAT TRIMS LEADING OR TRAILING WHITESPACE
trim <- function (x) gsub("^\\s+|\\s+$", "", x)


# FUNCITON TO REMOVE ALL SPACES FROM LEVEL NAMES OF A VARIABLE
rm_space <- function(df,col_name){
  level_names <- unique(levels(df[,which(names(df) %in% col_name)]))
  df[,which(names(df) %in% col_name)] <- mapvalues(df[,which(names(df) %in% col_name)], 
                                                   from = level_names,
                                                   to = gsub("[[:space:]]","",level_names))
  return(df)
}


# FUNCTION TO GET GEODESIC DISTANCE IN KM
gcd.slc <- function(long1, lat1, long2, lat2) {
  R <- 6371 # Earth mean radius [km]
  d <- acos(sin(lat1)*sin(lat2) + cos(lat1)*cos(lat2) * cos(long2-long1)) * R
  return(d) # Distance in km
}


# BUILD HISTOGRAM PLOTS BY 2-CLUSTERS BY VARIABLE
histogram2 <- function(df1,df2,var_name,breaks){
  # extract data
  data1 <- df1[[var_name]]
  data2 <- df2[[var_name]]
  data1 <- (data1-min(data1))/max(data1-min(data1))
  data2 <- (data2-min(data2))/max(data2-min(data2))
  cnt1 <- hist(data1, breaks = breaks, plot = FALSE)$counts
  brk1 <- hist(data1, breaks = breaks, plot = FALSE)$breaks
  cnt2 <- hist(data2, breaks = breaks, plot = FALSE)$counts
  brk2 <- hist(data2, breaks = breaks, plot = FALSE)$breaks
  pct1 <- 100*cnt1/sum(cnt1)
  pct2 <- 100*cnt2/sum(cnt2)
  bardata <- cbind(pct1,pct2)
  barplot(t(bardata), 
          beside = T,
          xlab = paste("VDC Index For",var_name),
          ylab = "Relative Percentages",
          main = paste("Distribution Comparison For",var_name),
          col = c("black","red"))
  axis(1, 
       at = 3*(1:length(brk1)-1)-1,
       labels = NA,
       cex.axis = 1,
       las = 1)
  legend("topright",
         legend = c("Cluster 1","Cluster 2"),
         fill = c("black","red"),
         bty = "n",
         cex = 1.25)
}


# BUILD HISTOGRAM PLOTS BY 3-CLUSTERS BY VARIABLE
histogram3 <- function(df1,df2,df3,var_name,breaks){
  # extract data
  data1 <- df1[[var_name]]
  data2 <- df2[[var_name]]
  data3 <- df3[[var_name]]
  data1 <- (data1-min(data1))/max(data1-min(data1))
  data2 <- (data2-min(data2))/max(data2-min(data2))
  data3 <- (data3-min(data3))/max(data3-min(data3))
  cnt1 <- hist(data1, breaks = breaks, plot = FALSE)$counts
  brk1 <- hist(data1, breaks = breaks, plot = FALSE)$breaks
  cnt2 <- hist(data2, breaks = breaks, plot = FALSE)$counts
  brk2 <- hist(data2, breaks = breaks, plot = FALSE)$breaks
  cnt3 <- hist(data3, breaks = breaks, plot = FALSE)$counts
  brk3 <- hist(data3, breaks = breaks, plot = FALSE)$breaks
  pct1 <- 100*cnt1/sum(cnt1)
  pct2 <- 100*cnt2/sum(cnt2)
  pct3 <- 100*cnt3/sum(cnt3)
  bardata <- cbind(pct1,pct2,pct3)
  barplot(t(bardata), 
          beside = T,
          xlab = paste("VDC Index For",var_name),
          ylab = "Relative Percentages",
          main = paste("Distribution Comparison For",var_name),
          col = c("black","red","green"))
  axis(1, 
       at = 4*(1:(length(brk1)-1))-2,
       labels = NA,
       cex.axis = 1,
       las = 1)
  legend("topright",
         legend = c("VDC Cluster 1","VDC Cluster 2","VDC Cluster 3"),
         fill = c("black","red","green"),
         bty = "n",
         cex = 1.25)
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


# QUAKE STATS
quake_stats <- read.csv(paste0(DIR,"quake_stats.csv"))

# NEPAL AFFECTED POPULATIN (DISTRICT LEVEL)
affected_pop <- read.csv(paste0(DIR,"affected_pop.csv"))

# NEPAL HAZARD SCORE
n_hazard <- read.csv(paste0(DIR,"hazard_score.csv"))

# NEPAL HEALTH FACILITIES
nepal_h <- read.csv(paste0(DIR,"nepal_health.csv"))

# LOAD POPULATION CENSUS TABLE
popt <- read.csv(paste0(DIR,"npl-popt.csv"))

# READ IN PCODE TO HLCIT TABLE
p_to_h <- read.csv(paste0(DIR,"p_codes_to_hlcit.csv"))

# RAW SEVERITY DATA
sev <- read.csv(paste0(DIR,"severity.csv"))

# LOAD LAT/LON COORDINATES (OF CENTROIDS) AND HLCIT CODES
centroids <- read.csv(paste0(DIR,"centroids.csv"))
hlcit <- read.csv(paste0(DIR,"master_hlcit.csv"))
colnames(hlcit) <- c("lon","lat","vdc_name","vname","hlcit_code")
hlcit$hlcit_code <- as.factor(hlcit$hlcit_code)
hlcit$vname <- as.character(hlcit$vname)
hlcit$vdc_name <- as.character(hlcit$vdc_name)
hlcit <- rm_space(hlcit,"hlcit_code")
hlcit$hlcit_code <- as.numeric(levels(hlcit$hlcit_code))[hlcit$hlcit_code]

# READ SEVERITY TABLE
severity_data <- readObj(file_name = paste0(DIR,"severity_mapvalues.df"))

# READ DISASTER AID RELIEF
aid_data <- read.csv(paste0(DIR,"agency_relief.csv"), sep=",")

# READ IN AID AND SEVERITY COMPLETE DATA TABLE
aid_sev <- readObj(file_name = paste0(DIR,"aid_and_severity.df"))

# READ MODELING TABLE
aid_sev_modeling <- readObj(file_name=paste0(DIR,"aid_sev_modeling.df"))


# SET NAs TO MEAN
aid_sev_modeling <- nas_to_mean(aid_sev_modeling)

# CONVERT EGREES TO RADIANS
deg2rad <- function(deg) return(deg*pi/180)


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
# Severity = (Hazard x Exposure × Vulnerability)^ 1/3
#
#
#
# NEED FOR AID  IS PROPORTIONAL TO
# severity
# 1/dist_to_epicenter
# 1/density of health care facilities
# 
#
#
# WITH FIXED APRIL 25 LT AND LON
#
#
#
#
#
#
#
#
#
# BUILD THE NEED ATTRIBUTE MODELING TABLE
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

# TRANSFORM HOSPITAL TYPE VARIABLE
nepal_h$hf_type <- mapvalues(x = nepal_h$hf_type,
                             from = c("Sub Center","DPHO","Sub Health Post","Health Post","Hospital","Supply Center","Primary Health Center",
                                      "District Cold Room","Laxmipur","DIstrict Cold Room","Ayurvedic Aushadhalaya","DAHC","","Distict Cold Room",
                                      "Zonal Hospital","Private Hospital","Health Care Center","Health Center","Refugee Camp","Primary Health Post",
                                      "Central Hospital","District Center","D(P)HO","District Ayurvedic HC","RMS"),
                             to = c("Four","Five","One","Two","Four","Two","Three","Two","Two","Two","Three","Four",NA,"Two","Five","Four","Three","Three","Two","Two","Five","Five","Five","Four","Four"))
nepal_h$hf_type <- mapvalues(x = nepal_h$hf_type,
                             from = c("One","Two","Three","Four","Five"),
                             to = 1:5)
nepal_h$hf_type <- as.numeric(levels(nepal_h$hf_type))[nepal_h$hf_type]
nepal_h$hf_type[is.na(nepal_h$hf_type)] <- 2


# TRANSFORM DISASTER AID DATA
aid_data$vdc <- trim(as.character(aid_data$vdc))
aid_data <- rm_space(aid_data,"vdc")
aid_data$vdc <- substr(as.character(aid_data$vdc),1,12)
aid_data$impl_ag <- trim(as.character(aid_data$impl_ag))
aid_data <- rm_space(aid_data,"impl_ag")
aid_data$impl_ag <- substr(as.character(aid_data$impl_ag),1,12)
aid_data <- aid_data[nchar(aid_data$vdc)>0 & nchar(aid_data$impl_ag)>0,]
aid_data <- rm_space(aid_data,"hlcit")
aid_data$hlcit <- as.numeric(levels(aid_data$hlcit))[aid_data$hlcit]
for (k in 1:dim(aid_data)[1]){
  aid_data$vdc[k] <- hlcit$vdc_name[which(hlcit$hlcit_code %in% aid_data$hlcit[k])[1]]
}
aid_data <- aid_data[!is.na(aid_data$hlcit),]
aid_data[is.na(aid_data$no_hh),]$no_hh <- mean(aid_data[!is.na(aid_data$no_hh),]$no_hh)
aid_data[is.na(aid_data$ave_cost),]$ave_cost <- mean(aid_data[!is.na(aid_data$ave_cost),]$ave_cost)
aid_data$no_actions <- log(1+aid_data$no_actions)


# TRANSFORM ALL HLCIT CODES
p_to_h <- rm_space(p_to_h,"hlcit")
p_to_h$hlcit <- as.numeric(levels(p_to_h$hlcit))[p_to_h$hlcit]


# ADD HLCIT CODES TO SEVERITY TABLE
sev_hlcit <- sev
for (k in 1:dim(sev_hlcit)[1]){
  if (sev_hlcit$p_code[k] %in% p_to_h$p_code){
    sev_hlcit$hlcit[k] <- mean(p_to_h[p_to_h$p_code %in% sev_hlcit$p_code[k],]$hlcit)
  }
}


# ADD LAT AND LON OF CENTROIDS TO SEVERITY TABLE
need_attribute_table <- sev_hlcit
for (k in 1:dim(need_attribute_table)[1]){
  if (need_attribute_table$hlcit[k] %in% hlcit$hlcit_code){
    need_attribute_table$lon[k] <- mean(hlcit[hlcit$hlcit_code %in% need_attribute_table$hlcit[k],]$lon)
    need_attribute_table$lat[k] <- mean(hlcit[hlcit$hlcit_code %in% need_attribute_table$hlcit[k],]$lat)
  } else {
    need_attribute_table$lon[k] <- NA
    need_attribute_table$lat[k] <- NA
  }
}


# ADDING POPULATION DENSITY WITH NEED_ATTRIBUTE_TABLE
popt1 <- popt[,c("P_CODE","Popden2011")]
colnames(popt1) <- c("p_code","pop_density")
for (k in 1:dim(need_attribute_table)[1]){
  if (need_attribute_table$p_code[k] %in% popt1$p_code){
    need_attribute_table$pop_density[k] <- mean(popt1[popt1$p_code %in% need_attribute_table$p_code[k],]$pop_density)
  }
}


# DROP TWO DUPLICATING VDCs WITH DIFFERENT NAMES
need_attribute_table <- need_attribute_table[-c(3931,3937),]


# ADD HELATHCARE FACILITIES
for (k in 1:dim(need_attribute_table)[1]){
  if (need_attribute_table$p_code[k] %in% nepal_h$vdc_code1){
    h_type <- nepal_h[nepal_h$vdc_code1 %in% need_attribute_table$p_code[k],]$hf_type
    need_attribute_table$hc_cnt[k] <- length(h_type)
    need_attribute_table$hc_wt_cnt[k] <- sum(h_type)
  } else {
    need_attribute_table$hc_cnt[k] <- NA
    need_attribute_table$hc_wt_cnt[k] <- NA
  }
}

# RESOLVE THE NAs
need_attribute_table$hc_cnt[is.na(need_attribute_table$hc_cnt)] <- median(need_attribute_table$hc_cnt[!is.na(need_attribute_table$hc_cnt)])
need_attribute_table$hc_wt_cnt[is.na(need_attribute_table$hc_wt_cnt)] <- median(need_attribute_table$hc_wt_cnt[!is.na(need_attribute_table$hc_wt_cnt)])



# ADD NEPAL HAZARD SCORE
for (k in 1:dim(need_attribute_table)[1]){
  if (need_attribute_table$p_code[k] %in% n_hazard$p_code){
    need_attribute_table$hazard_score[k] <- mean(n_hazard[n_hazard$p_code %in% need_attribute_table$p_code[k],]$index_cnt_vdc)
  }
}

# INSERT DISTANCE TO EPICENTER
q_lat <- 28.1473
q_lon <- 84.7079
for (k in 1:dim(need_attribute_table)[1]){
  need_attribute_table$dist_epicenter[k] <- sqrt((need_attribute_table$lon[k]-q_lon)^2+(need_attribute_table$lat[k]-q_lat)^2)
}

# INTERMEDIATE SAVE
write.csv(need_attribute_table,file=paste0(DIR,"need_attribute_table.csv"))
writeObj(need_attribute_table,file=paste0(DIR,"need_attribute_table.df"))




#
#
# READ THE MODELIN TABLE
#
#

need_attribute_table <- readObj(file_name = paste0(DIR,"need_attribute_table.df"))


# DEFINE NEED


need_attribute_table$severity <- scale(need_attribute_table$severity)-
  min(scale(need_attribute_table$severity))
need_attribute_table$hc_wt_cnt <- scale(need_attribute_table$hc_wt_cnt)-
  min(scale(need_attribute_table$hc_wt_cnt))
need_attribute_table$dist_epicenter <- scale(need_attribute_table$dist_epicenter)-
  min(scale(need_attribute_table$dist_epicenter))


need_attribute_table$need <- 3*(need_attribute_table$severity)/
  (((need_attribute_table$hc_wt_cnt/6+0.1)^(1/3))*(need_attribute_table$dist_epicenter+0.1)^(1/3))



summary(need_attribute_table$need)
# OBSERVE THREE DIFFERENT GROUPS
hist(need_attribute_table$need[need_attribute_table$need < 67], 
       breaks = 100,
       col = rev(heat.colors(1+1.5*max(as.integer(need_attribute_table$need[need_attribute_table$need < 60])))),
       xlab = "VDC NEED FOR AID RANK",
       ylab = " FREQUENCY OF NEED RANK OCCURENCE",
       main = " DISTRIBUTION OF NEED FOR ALL VDCs")






# A VIEW OF THE DISTRIBUTION WITH THE GROUPING BY NEED
hist(need_attribute_table$need[need_attribute_table$need < 43], 
       breaks = 100,
       col = rev(heat.colors(1+1.85*max(as.integer(need_attribute_table$need[need_attribute_table$need < 43])))),
       xlab = "VDC NEED FOR AID RANK",
       ylab = " FREQUENCY OF NEED RANK OCCURENCE",
       main = " DISTRIBUTION OF NEED FOR ALL VDCs")
rect(xleft = 0,
     xright = 4.75,
     ybottom = 0,
     ytop = 810,
     border = "blue",
     density = 7,
     col = "blue",
     lwd = 0.8)
text(3,600,
      "LOW
NEED
63.9%",
      col = "blue",
      cex = 1.5,
     font = 2)
rect(xleft = 4.75,
     xright = 20,
     ybottom = 0,
     ytop = 300,
     border = "blue",
     density = 7,
     col = "blue",
     lwd = 0.8)
text(12,200,
     "MED
NEED
30.7%",
     col = "blue",
     cex = 1.5,
     font = 2)
rect(xleft = 20,
     xright = 40,
     ybottom = 0,
     ytop = 200,
     border = "blue",
     density = 7,
     col = "blue",
     lwd = 0.8)
text(30,75,
     "HIGH NEED 5.4%",
     col = "blue",
     cex = 1.5,
     font = 2)




# SEPARATE THE DIFFERNET CLUSTERS:
low_need <- need_attribute_table[need_attribute_table$need <= 4.75,]
med_need <- need_attribute_table[need_attribute_table$need > 4.75 & need_attribute_table$need <= 20,]
high_need<- need_attribute_table[need_attribute_table$need > 20,]

# COMPUTE THE PERCENTAGES
length(unique(low_need$hlcit))/length(unique(need_attribute_table$hlcit))
length(unique(med_need$hlcit))/length(unique(need_attribute_table$hlcit))
length(unique(high_need$hlcit))/length(unique(need_attribute_table$hlcit))








#
#
#
#
#
#
#
# VISUALIZE THE THREE CLUSTER GROUPS
#
#
#
#
#
#
#
#
#

# SELECT UNIQUE AGENCIES AND TARGET VDC
ag <- unique(aid_data$impl_ag)
hlc <- unique(aid_data$hlcit)
all <- union(ag,hlc)
koords<-cbind(hlcit$lon,hlcit$lat)
u_hl <- unique(hlcit$hlcit_code)

# DEFINE THE AGENCY-VDC AID NETWORK ADJACENCY MATRIX
aid_m <- matrix(0,
                nrow = length(u_hl),
                ncol = length(u_hl))

# BUILD THE AGENCY-VDC AID NETWORK
av <- graph.adjacency(aid_m)


# COLOR VERTICES REPRESENTING AGENCIES (GREEN) AND VDCs (BLUE) WHERE AID WAS SENT
V(av)$color <- rep(NA,length(u_hl))
for (k in 1:length(unique(hlcit$hlcit_code))){
  if(is.element(unique(hlcit$hlcit_code)[k],high_need$hlcit) & 
     is.element(unique(hlcit$hlcit_code)[k],aid_data$hlcit)){
    V(av)$color[k] <- "red"
  }  
  if(is.element(unique(hlcit$hlcit_code)[k],med_need$hlcit)& 
     is.element(unique(hlcit$hlcit_code)[k],aid_data$hlcit)){
    V(av)$color[k] <- "orange"
  } 
  if(is.element(unique(hlcit$hlcit_code)[k],low_need$hlcit) & 
     is.element(unique(hlcit$hlcit_code)[k],aid_data$hlcit) &
     (need_attribute_table$dist_epicenter[k] < quantile(need_attribute_table$dist_epicenter,0.85))){
    V(av)$color[k] <- "yellow"
  } 
  if (need_attribute_table$lat[k]>q_lat+0.7*q_lat){
    V(av)$color[k] <- NA
  }
}


plot(av,
     layout = koords,
     vertex.color = V(av)$color,
     vertex.size = 2,
     vertex.label = NA, 
     vertex.label.color = "darkgreen", 
     vertex.label.font = 1, 
     vertex.label.cex = 0.75
)



# COLOR VERTICES REPRESENTING AGENCIES (GREEN) AND VDCs (BLUE) WHERE AID WAS SENT
V(av)$color <- rep(NA,length(u_hl))
for (k in 1:length(unique(hlcit$hlcit_code))){
  if(is.element(unique(hlcit$hlcit_code)[k],high_need$hlcit)){
    V(av)$color[k] <- "red"
  }  
  if(is.element(unique(hlcit$hlcit_code)[k],med_need$hlcit)){
    V(av)$color[k] <- "orange"
  } 
  if(is.element(unique(hlcit$hlcit_code)[k],low_need$hlcit)){
    V(av)$color[k] <- "yellow"
  }
  if (need_attribute_table[need_attribute_table$hlcit %in% hlcit$hlcit_code[k],]$dist_epicenter>
      quantile(need_attribute_table$dist_epicenter,0.70)){
    V(av)$color[k] <- "white"
  }
}


plot(av,
     layout = koords,
     vertex.color = V(av)$color,
     vertex.size = 2,
     vertex.label = NA, 
     vertex.label.color = "darkgreen", 
     vertex.label.font = 1, 
     vertex.label.cex = 0.75
)



# CONTINUOUS HEAT MAP






V(av)$color <- rep(NA,length(u_hl))
for (k in 1:length(unique(hlcit$hlcit_code))){
  V(av)$color[k] <- rev(heat.colors(1+max(as.integer(1.2*need_attribute_table$need))))[as.integer(need_attribute_table$need)[k]+1]
  
}


plot(av,
     layout = koords,
     vertex.color = V(av)$color,
     vertex.size = 2,
     vertex.label = NA, 
     vertex.label.color = "darkgreen", 
     vertex.label.font = 1, 
     vertex.label.cex = 0.75
)




#
#
#
#
#
#
#
#
#
# LAST ONE
#
#
#
#
#
#
#
#
#





# EXTRAPOLATE LHCIT NUMBERS FROM LAT AND LON FOR VDCS FROM HLCIT_MASTER
ag <- unique(aid_data$impl_ag)
hlc <- unique(aid_data$hlcit)
all <- union(ag,hlc)

# EXTRAPOLATE LHCIT NUMBERS FROM LAT AND LON FOR VDCS FROM HLCIT_MASTER
hl <- vector()
xc <- vector()
yc <- vector()
for (k in 1:length(hlc)){
  hl[k] <- hlc[k]
  xc[k] <- hlcit$lon[which(hlcit$hlcit_code==hlc[k])[1]]
  yc[k] <- hlcit$lat[which(hlcit$hlcit_code==hlc[k])[1]]
}
koords3 <- cbind(xc,yc)

xa <- 71+9*runif(length(ag))
ya <- 25+8*runif(length(ag))
koords1 <- cbind(xa,ya)
koords2<- rbind(koords1,koords3)


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


for (k in 1:dim(aid_m)[1]){
  if(k-1<length(ag)){
    V(av)$size[k] <- 3
    V(av)$name[k] <- ag[k]
  } else {
    V(av)$size[k] <- 1
    V(av)$name[k] <- NA}
}


# PLOT THE AGENCY-VDC AID NETWORK
V(av)$color <- rep("green",length(all))
for (k in 1:length(unique(hlcit$hlcit_code))){
  if(is.element(unique(hlcit$hlcit_code)[k],high_need$hlcit) & 
     is.element(unique(hlcit$hlcit_code)[k],aid_data$hlcit)){
    V(av)$color[k] <- "red"
  }  
  if(is.element(unique(hlcit$hlcit_code)[k],med_need$hlcit)& 
     is.element(unique(hlcit$hlcit_code)[k],aid_data$hlcit)){
    V(av)$color[k] <- "orange"
  } 
  if(is.element(unique(hlcit$hlcit_code)[k],low_need$hlcit) & 
     is.element(unique(hlcit$hlcit_code)[k],aid_data$hlcit) &
     (need_attribute_table$dist_epicenter[k] < quantile(need_attribute_table$dist_epicenter,0.85))){
    V(av)$color[k] <- "yellow"
  } 
  if (need_attribute_table$lat[k]>q_lat+0.7*q_lat){
    V(av)$color[k] <- NA
  }
}



plot(av,
     layout = koords2,
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
       c("Implementing Aid Agency","VDC With Geo-Coordinates"),
       fill = c("green","darkorange"),
       bty = "n")




















#
#
#
#
#
#
#
#
# AID RECEIVED
#
#
#
#
#
#
#
#
#




# IDENTIFY VDC FROM EACH CLUSTER THAT RECEIVED AID:
aid_hlcit <- unique(aid_data$hlcit)
aid_hlcit_low <- aid_hlcit[aid_hlcit %in% unique(low_need$hlcit)]
aid_hlcit_med <- aid_hlcit[aid_hlcit %in% unique(med_need$hlcit)]
aid_hlcit_high <- aid_hlcit[aid_hlcit %in% unique(high_need$hlcit)]


# CHECK RELATIVE SIZES
length(aid_hlcit_3a)
length(aid_hlcit_3b)
length(aid_hlcit_3c)
dim(clust_3a)[1]
dim(clust_3b)[1]
dim(clust_3c)[1]


# CREATE THE THREE AID TABLES
aid_data_3a <- aid_data[aid_data$hlcit %in% aid_hlcit_3a,]
aid_data_3b <- aid_data[aid_data$hlcit %in% aid_hlcit_3b,]
aid_data_3c <- aid_data[aid_data$hlcit %in% aid_hlcit_3c,]


# COUNT WEIGHTED SUM OF INSTANCES OF AID
length(aid_hlcit_low)/dim(clust_3ah)[1]
length(aid_hlcit_med)/dim(clust_3bh)[1]
length(aid_hlcit_high)/dim(clust_3ch)[1]


length(unique(low_need$hlcit))/length(unique(need_attribute_table$hlcit))
length(unique(med_need$hlcit))/length(unique(need_attribute_table$hlcit))
length(unique(high_need$hlcit))/length(unique(need_attribute_table$hlcit))

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



















# BUILD COORDS FOR ALL VDCS
u_hl <- unique(hlcit$hlcit_code)
kds <- cbind(hlcit$lon,hlcit$lat)

# DEFINE THE AGENCY-VDC AID NETWORK ADJACENCY MATRIX
aid_add <- matrix(0,
                nrow = length(u_hl),
                ncol = length(u_hl))
for (i in 1:length(u_hl)){
    aid_add[[i,i]] <- 1
}


# BUILD THE AGENCY-VDC AID NETWORK
av_add <- graph.adjacency(aid_add)

# DROP SELF-LOOPS
av_add <- drop_loops(graph = av_add,
                 vertex_colors = NA,
                 vertex_names = NA)


plot(av_add,
     layout = kds,
     vertex.color = NA,
     vertex.size = 1,
     vertex.label = NA)

uav <- graph.disjoint.union(av_add,av)

plot(uav,
     layout=rbind(kds, koords2),
     vertex.label = NA,
     vertex.size = 1,
     edge.arrow.size = 0.1,
     edge.width = 0.1)


















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
# BUILD THE AID ATTRIBUTE MODELING TABLE
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


















# MERGE WITH AGENCY-VDC AID TABLE
for (k in 1:dim(sev)[1]){
  if (sev$vdc[k] %in% hlcit$vname){
    sev$hlcit[k] <- hlcit[which(hlcit$vname==sev$vdc[k])[1],]$hlcit
  } else {
    sev$hlcit[k] <- NA
  }
}

# EXPORT TRANSFORMED SEVERITY TABLE
write.csv(sev,file=paste0(DIR,"severity_mapvalues.csv"))
writeObj(sev,file=paste0(DIR,"severity_mapvalues.df"))


# RESOLVE REMAINING VDC NAMES AND HLCIT CODES
resolve_vdc <- sev[is.na(sev$hlcit),]$vdc
resolve_which <- which(is.na(sev$hlcit))
for (k in 1: length(resolve_vdc)){
  if (sev$vdc[resolve_which[k]] %in% hlcit$vdc_name){
    sev$hlcit[resolve_which[k]] <- hlcit[which(hlcit$vdc_name==sev$vdc[resolve_which[k]])[1],]$hlcit
  } else {
    sev$hlcit[resolve_which[k]] <- NA
  }
}

# NOW WE MERGE WITH THE SEVERITY INDEX DATA
# GET UNIQUE HLCIT FROM AID DATA
aid_sev <- aid_data
for (k in 1:dim(aid_data)[1]){
  aid_sev$hazard[k] <- mean(sev[sev$hlcit %in% aid_data$hlcit[k],]$hazard)
  aid_sev$exposure[k] <- mean(sev[sev$hlcit %in% aid_data$hlcit[k],]$exposure)
  aid_sev$housing[k] <- mean(sev[sev$hlcit %in% aid_data$hlcit[k],]$housing)
  aid_sev$poverty[k] <- mean(sev[sev$hlcit %in% aid_data$hlcit[k],]$poverty)
  aid_sev$vulnerability[k] <- mean(sev[sev$hlcit %in% aid_data$hlcit[k],]$vulnerability)
  aid_sev$severity[k] <- mean(sev[sev$hlcit %in% aid_data$hlcit[k],]$severity)
}

# ADD HLCIT AND VDC DEGREE VARIABLE
vdc_degree <- readObj(paste0(DIR,"vdc_degree.df"))
hlcit_degree <- readObj(paste0(DIR,"hlcit_degree.df"))
for (k in 1:dim(aid_sev)[1]){
  if (aid_sev$hlcit[k] %in% hlcit_degree$hlcit){
    aid_sev$degree[k] <- hlcit_degree[hlcit_degree$hlcit %in% aid_sev$hlcit[k],]$hlcit_degree
  } else {
    aid_sev$degree[k] <- NA
  }
}

# EXPORT THE DATA
write.csv(aid_sev,file=paste0(DIR,"aid_and_severity.csv"))
writeObj(aid_sev,file=paste0(DIR,"aid_and_severity.df"))


# EXPORT JUST SEVERITY DATA ON THE AID RECEIVING VDCs
aid_hlcit_list <- unique(aid_sev$hlcit)
var_list <- c("hazard","exposure","housing","poverty","vulnerability","severity","degree")   
aid_sev_modeling <- unique(aid_sev[aid_sev$hlcit %in% aid_hlcit_list,var_list])

# SAVE MODELING TABLE
write.csv(aid_sev_modeling,file=paste0(DIR,"aid_sev_modeling.csv"))
writeObj(aid_sev_modeling,file=paste0(DIR,"aid_sev_modeling.df"))









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
# BUILD THE AGENCY-VDC AID NETWORK
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


















