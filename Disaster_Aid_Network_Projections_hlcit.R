# Nepal Disaster Aid Distribution Network Analysis and Projections BY HLCIT CODES
# author: Georgi D. Gospodinov
# date: "September 19, 2015"
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
# FUNCTION TO DISPLAY RELATIVE PERCENTAGES FOR HSITOGRAM COLUMNS
histP <- function(x,breaks, ...) {
  H <- hist(x, plot = FALSE, breaks=breaks)
  H$density <- with(H, 100 * density* diff(breaks)[1])
  labs <- paste(round(H$density), "%", sep="")
  plot(H, freq = FALSE, labels = labs, ylim=c(0, 1.08*max(H$density)),...)
}


# FUNCTION TO DISPLAY RELATIVE PERCENTAGE VALUES GREATER THAN 0
histP1 <- function(x,breaks, ...) {
  H <- hist(x, plot = FALSE, breaks=breaks)
  H$density <- with(H, 100 * density* diff(breaks)[1])
  labs <- ifelse(round(H$density)>0,paste(round(H$density), "%", sep=""),NA)
  plot(H, freq = FALSE, labels = labs, ylim=c(0, 1.08*max(H$density)),...)
}


# FUNCITON TO DISPLAY RELATIVE PERCENTAGE VALUES GREATER THAN 5, CAN BE HARDCODED
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

# FUNCTION THAT DROPS ISOLATED VERTICES
drop_isolated <- function(graph, vertex_colors, vertex_names) {
  
  # get the definitions
  g <- graph
  
  # filter to degree > 0 eliminate isolated vertices
  g_f <- delete.vertices(g,V(g)[degree(g)==0])
  v_g_f <- setdiff(V(g),V(g)[degree(g)==0])
  
  # filter names and color
  V(g_f)$name <- vertex_names[v_g_f]
  V(g_f)$color <- vertex_colors[v_g_f]
  return(g_f)
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


# DEFINE EDGE-FILTRATION FUNCTION FOR THE NETWORKS
filter <- function(cutoff,edge_matrix,vertex_colors,vertex_names, vertex_size) {
  
  # get the definitions
  cut <- cutoff
  adj <- edge_matrix
  adj[adj<cut] <- 0
  adj_0 <- adj
  
  # define the filtered graph
  g <- graph.adjacency(adj_0,mode="directed",weighted=TRUE)
  V(g)$color <- vertex_colors
  
  # filter to degree > 0 eliminate isolated vertices
  g_f <- delete.vertices(g,V(g)[degree(g)==0])
  v_g_f <- setdiff(V(g),V(g)[degree(g)==0])
  V(g_f)$name <- vertex_names[v_g_f]
  V(g_f)$color <- vertex_colors[v_g_f]
  V(g_f)$size <- vertex_size[v_g_f]
  
  return(g_f)
}


# DEFINE DEGREE FILTRATION FUNCTION FOR THE NETWORKS
filter_deg <- function(cutoff,edge_matrix,vertex_colors,vertex_names) {
  
  # get the definitions
  cut <- cutoff
  adj <- edge_matrix
  
  # set the cut-off
  adj_0 <- adj
  adj_0[adj_0>0] <- 1
  for (i in 1:dim(adj)[1]){
    if (sum(adj_0[,i])<cutoff){
      adj <- adj[,-i]
      adj <- adj[-i,]
    }
  }
  
  # define the filtered graph
  g <- graph.adjacency(adj,mode="directed",weighted=TRUE)
  V(g)$color <- vertex_color
  
  # filter to degree > 0 eliminate isolated vertices
  g_f <- delete.vertices(g,V(g)[degree(g)==0])
  v_g_f <- setdiff(V(g),V(g)[degree(g)==0])
  V(g_f)$name <- vertex_names[v_g_f]
  
  # color the filtered graph with sources and endpoints for the directed edges
  V(g_f)$color <- V(g)$colors[v_g_f]
  
  
  return(g_f)
}

# DEFINE WEIGHTED DEGREE FILTRATION USING graph.strength

# DEFINE A FILTRATION OF GRAPH TO DISPLAY LARGEST CLUSTER (GIANT COMPONENT)
giant_comp <- function(graph, vertex_colors, vertex_names, vertex_size){
  
  # get the definitions
  g <- graph
  
  # identify the largest cluster
  clusters <- as.data.frame(table(clusters(g)$membership))
  ind <- as.numeric(clusters[which(clusters$Freq==max(clusters$Freq)),]$Var1)
  vertices <- which(clusters(g)$membership==ind)
  vertices_complement <- which(clusters(g)$membership!=ind)
  
  # filter to only include the giant component
  g_f <- delete.vertices(g,V(g)[which(V(g) %in% vertices_complement)])
  v_g_f <- setdiff(V(g),V(g)[which(V(g) %in% vertices_complement)])
  V(g_f)$name <- vertex_names[v_g_f]
  
  # color the filtered graph with sources and endpoints for the directed edges
  V(g_f)$color <- vertex_colors[v_g_f]
  V(g_f)$size <- vertex_size[v_g_f]
  
  return(g_f)
}


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
for (k in 1:length(vd)){
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
legend("topleft",
       c("Implementing Aid Agency","VDC With Geo-Coordinates"),
       fill = c("green","SkyBlue2"),
       bty = "n")

#
#
#
#
#
#
#
# ANALYSIS OF THE AGENCY-VDC AID NETWORK
#
#
#
#
#
#
#

# DEFINE THE AGENCY-VDC AID NETWORK
av <- graph.adjacency(aid_m,
                      mode = "directed",
                      weighted = TRUE)

# COLOR VERTICES REPRESENTING AGENCIES (GREEN) AND VDCs (BLUE) WHERE AID WAS SENT
V(av)$color <- rep("green",length(all))
for (k in 1:length(all)){
  if(is.element(all[k],hlc)){
    V(av)$color[k]<-"SkyBlue2"
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


# THIS IS THE NUMBER OF AGENCIES PER VDC
summary(degree(av,mode = "in"))

# THIS IS THE NUMBER OF VDCs PER AGENCY
summary(degree(av,mode = "out"))

# THIS IS THE WEIGHTED NUMBER OF AGENCIES,
# WEIGHED BY NUMBER OF AID ACTIONS PER VDC
summary(graph.strength(av,mode = "in"))


# THIS IS THE WEIGHTED NUMBER OF VDCs PER AGENCY
# WEIGHED BY THE NUMBER OF AID ACTIONS
summary(graph.strength(av,mode = "out"))

# PLOT THE AGENCY AND VDC DEGREE DISTRIBUTIONS 
plot(sort(degree(av,mode = "in")),
     col = "green",
     pch = 19,
     xaxt = "n",
     yaxt = "n",
     ann = FALSE)
par(new=T)
plot(sort(degree(av,mode = "out")),
     col = "SkyBlue2",
     pch = 19,
     xlab = "Agency/VDC Index",
     ylab = "Agency/VDC Degree",
     main = "Agency/VDC Degree (Sorted)")
legend("topleft",
       c("Number of VDCs per Agency", "Number of Agencies per VDC"),
       fill = c("green","SkyBlue2"),
       cex = 1.8,
       bty = "n")
text(300,200, paste("  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
                    0.000   0.000   0.000   4.782   0.000 337.000 "),
     col = "green",
     cex = 1.6,
     font = 2)
text(300,150, paste("   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
                    0.000   1.000   3.000   4.782   6.000  64.000"),
     col = "SkyBlue2",
     cex = 1.6,
     font = 2)


# THE DEGREE DISTRIBUTION (AGENCY-VDC CONNECTIONS)
deg1 <- hist(degree(av,mode = "in"), breaks = 20)$counts
deg1n <- 100*deg1/sum(deg1)
deg2 <- hist(degree(av,mode = "out"), breaks = 20)$counts
deg2n <- 100*deg2/sum(deg2)
maxn <- max(length(deg1n), length(deg2n))
d1 <- append(deg1n,rep(0,maxn-length(deg1n)))
d2 <- append(deg2n,rep(0,maxn-length(deg2n)))
ddata <- cbind(d1,d2)
barplot(t(ddata), 
        beside = T, 
        xlab = "Number of Agency/VDC Connections for a VDC/Agency", 
        ylab = "Relative Agency/VDC Connection Percentages", 
        main = "Distribution of Relative Agency/VDC Connections",
        col = c("green","SkyBlue2"))
axis(1, 
     at = 4*(1:maxn)-2,
     labels = as.character(1:maxn),
     cex.axis = 1,
     las = 1)
legend("topright",
       legend = c("Number of VDCs per Agency", "Number of Agencies per VDC"),
       col = c("green","SkyBlue2"), 
       bty = "n",pch = 15, 
       cex = 1.5)


# THE WEIGHTED DEGREE DISTRIBUTION (AGENCY-VDC CONNECTIONS)
plot(sort(graph.strength(av,mode = "in")),
     col = "green",
     pch = 19,
     xaxt = "n",
     yaxt = "n",
     ann = FALSE)
par(new=T)
plot(sort(graph.strength(av,mode = "out")),
     col = "SkyBlue2",
     pch = 19,
     xlab = "Weighted Agency/VDC Index",
     ylab = "Weighted Agency/VDC Degree",
     main = "Weighted Number of Agency/VDC Connections (Sorted)")
legend("topleft",
       c("Weighted Number of VDCs per Agency", "Weighted Number of Agencies per VDC"),
       fill = c("green","SkyBlue2"),
       bty = "n")


# THE WEIGHTED DEGREE DISTRIBUTION (VDC-VDC CONNECTIONS)
deg1 <- hist(graph.strength(av,mode = "in"), breaks = 20)$counts
deg1n <- 100*deg1/sum(deg1)
deg2 <- hist(graph.strength(av,mode = "out"), breaks = 20)$counts
deg2n <- 100*deg2/sum(deg2)
maxn <- max(length(deg1n), length(deg2n))
d1 <- append(deg1n,rep(0,maxn-length(deg1n)))
d2 <- append(deg2n,rep(0,maxn-length(deg2n)))
wddata <- cbind(d1,d2)
barplot(t(wddata), 
        beside = T, 
        xlab = "Weighted Number of Agency/VDC Connections for a VDC/Agency",
        ylab = "Weighted Relative Agency/VDC Connection Percentages",
        main = "Distribution of Weighted Relative Agency/VDC Connections",
        col = c("green","SkyBlue2"))
axis(1, 
     at = 4*(1:maxn)-2,
     labels = as.character(1:maxn),
     cex.axis = 1,
     las = 1)
legend("topright",
       legend = c("Weighted Number of VDCs per Agency", "Weighted Number of Agencies per VDC"),
       col = c("green","SkyBlue2"), 
       bty = "n",pch = 15, 
       cex = 1.5)

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
     edge.width = 0.25*E(agg)$weight,
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
     edge.width = 0.15*E(agg)$weight,
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
     col = adjustcolor(rgb(1,0,1/2,1)))
hist(paths,
     breaks = 15,
     col = adjustcolor(rgb(1,0,1/2,1)),
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

# CLUSTERS ARE CONNECTED COMPONENTS, WE HAVE 4 in the UNFILTERED AGENCY-VDC NETWORK 
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


