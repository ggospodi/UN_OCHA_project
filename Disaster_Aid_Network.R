# Nepal Disaster Aid Distribution Network Analysis
# author: Georgi D. Gospodinov
# date: "Augist 11, 2015"
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
# Agency-VDC Aid NETWORK AT VDC LEVEL
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
vd <- unique(aid_data$vdc)
all <- union(ag,vd)

# EXTRAPOLATE LHCIT NUMBERS FROM LAT AND LON FOR VDCS FROM HLCIT_MASTER
hl <- vector()
xc <- vector()
yc <- vector()
for (k in 1:length(vd)){
    hl[k] <- hlcit$hlcit_code[which(hlcit$vdc_name==vd[k])[1]]
    xc[k] <- hlcit$lat[which(hlcit$vdc_name==vd[k])[1]]
    yc[k] <- hlcit$lon[which(hlcit$vdc_name==vd[k])[1]]
    }
koords<-cbind(xc,yc)

xa <- 20+5*runif(length(ag))
ya <- 80+12*runif(length(ag))
koords1 <-cbind(xa,ya)
koords2<- rbind(koords1,koords)

# DEFINE THE AGENCY-VDC AID NETWORK ADJACENCY MATRIX
aid_m <- matrix(0,nrow=length(all),ncol=length(all))
for (i in 1:length(ag)){
  for (j in 1:length(vd)){
    aid_m[[i,length(ag)+j]] <- 
      dim(aid_data[aid_data$impl_ag==ag[i] & aid_data$vdc==vd[j],c(3,5)])[1]
  }
}

# BUILD THE AGENCY-VDC AID NETWORK
av <- graph.adjacency(aid_m,
                      mode = "directed",
                      weighted = TRUE)

# COLOR VERTICES REPRESENTING AGENCIES (GREEN) AND VDCs (BLUE) WHERE AID WAS SENT
V(av)$color <- rep("green",length(all))
for (k in 1:length(all)){
  if(is.element(all[k],vd)){
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
                                          niter = 200,
                                          area = 2000*vcount(av)),
     vertex.color = V(av)$color,
     vertex.size = V(av)$size,
     vertex.label = NA, 
     vertex.label.color = "black", 
     vertex.label.font = 2, 
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
     layout = koords2,
     vertex.color = V(av)$color,
     vertex.size = V(av)$size,
     vertex.label = V(av)$name, 
     vertex.label.color = "darkgreen", 
     vertex.label.font = 2, 
     vertex.label.cex = 0.75, 
     edge.width = 0.05*sqrt(E(av)$weight),
     edge.arrow.size = 0.2,
     edge.curved = TRUE,
     edge.color = gray.colors(1),
     main = "Nepal Agency-VDC Aid Relief Geo-Network")
legend("topright",
       c("Implementing Aid Agency","VDC With Geo-Coordinates"),
       fill = c("green","SkyBlue2"),
       bty = "n")


# EDGE-FILTRATION BY EDGE WEIGHT OF THE AGENCY-VDC AID NETWORK: CUT-OFF = 25% PERCENTILE
for (k in 1:dim(aid_m)[1]){
  if(k-1<length(ag)){
    V(av)$size[k] <- 3
    V(av)$name[k] <- ag[k]
  } else {
    V(av)$size[k] <- 1
    V(av)$name[k] <- NA}
}
cut25 <- quantile(as.vector(aid_m[aid_m>0]),0.25)
av_f <- filter(cutoff = cut25,
             edge_matrix = aid_m,
             vertex_colors = V(av)$color,
             vertex_names = V(av)$name,
             vertex_size = V(av)$size)

# DISPLAY THE EDGE-FILTERED GRAPH
plot(av_f,
     layout = layout.fruchterman.reingold(av_f,
                                          niter = 200,
                                          area = 2000*vcount(av_f)),
     vertex.color = V(av_f)$color,
     vertex.size = V(av_f)$size,
     vertex.label = NA, 
     vertex.label.color = "black", 
     vertex.label.font = 2, 
     vertex.label.cex = 0.7, 
     edge.width = 0.7*sqrt(E(av_f)$weight),
     edge.arrow.size = 0.2,
     edge.curved = TRUE,
     edge.color = gray.colors(1),
     main = "Filtered Nepal Agency-VDC Aid Relief Abstract Network")
legend("topright",
       c("Implementing Aid Agency","Aid Target VDC"),
       fill = c("green","SkyBlue2"),
       bty = "n")

# PLOT THE AGENCY-VDC AID GEO-NETWORK
koords2_f <- koords2[which(V(av)$name %in% V(av_f)$name),]
plot(av_f,
     layout = koords2_f,
     vertex.color = V(av_f)$color,
     vertex.size = V(av_f)$size,
     vertex.label = V(av_f)$name, 
     vertex.label.color = "darkgreen", 
     vertex.label.font = 2, 
     vertex.label.cex = 0.75, 
     edge.width = 0.05*sqrt(E(av_f)$weight),
     edge.arrow.size = 0.2,
     edge.curved = TRUE,
     edge.color = gray.colors(1),
     main = "Filtered Nepal Agency-VDC Aid Relief Geo-Network")
legend("topright",
       c("Implementing Aid Agency","VDC With Geo-Coordinates"),
       fill = c("green","SkyBlue2"),
       bty = "n")


# EDGE-FILTRATION BY EDGE WEIGHT OF THE AGENCY-VDC AID NETWORK: CUT-OFF = 50% PERCENTILE
av <- graph.adjacency(aid_m,
                      mode = "directed",
                      weighted = TRUE)

# COLOR VERTICES REPRESENTING AGENCIES (GREEN) AND VDCs (BLUE) WHERE AID WAS SENT
V(av)$color <- rep("green",length(all))
for (k in 1:length(all)){
  if(is.element(all[k],vd)){
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
cut50 <- quantile(as.vector(aid_m[aid_m>0]),0.5)
av_f <- filter(cutoff = cut50,
             edge_matrix = aid_m,
             vertex_colors = V(av)$color,
             vertex_names = V(av)$name,
             vertex_size = V(av)$size)

# DISPLAY THE EDGE-FILTERED GRAPH
plot(av_f,
     layout = layout.fruchterman.reingold(av_f,
                                          niter = 200,
                                          area = 2000*vcount(av_f)),
     vertex.color = V(av_f)$color,
     vertex.size = V(av_f)$size,
     vertex.label = V(av_f)$name, 
     vertex.label.color = "black", 
     vertex.label.font = 2, 
     vertex.label.cex = 0.7, 
     edge.width = 0.3*sqrt(E(av_f)$weight),
     edge.arrow.size = 0.2,
     edge.curved = TRUE,
     edge.color = gray.colors(1),
     main = "Filtered Nepal Agency-VDC Aid Relief Abstract Network")
legend("topright",
       c("Implementing Aid Agency","Aid Target VDC"),
       fill = c("green","SkyBlue2"),
       bty = "n")

# PLOT THE AGENCY-VDC AID GEO-NETWORK
av <- graph.adjacency(aid_m,
                      mode = "directed",
                      weighted = TRUE)
V(av)$color <- rep("green",length(all))
for (k in 1:length(all)){
  if(is.element(all[k],vd)){
    V(av)$color[k] <- "SkyBlue2"
  }  
}
for (k in 1:dim(aid_m)[1]){
  if(k-1<length(ag)){
    V(av)$size[k] <- 3
    V(av)$name[k] <- all[k]
  } else {
    V(av)$size[k] <- 1
    V(av)$name[k] <- all[k]}
}
cut50 <- quantile(as.vector(aid_m[aid_m>0]),0.5)
av_f <- filter(cutoff = cut50,
               edge_matrix = aid_m,
               vertex_colors = V(av)$color,
               vertex_names = V(av)$name,
               vertex_size = V(av)$size)
koords2_f <- koords2[which(V(av)$name %in% V(av_f)$name),]
V(av_f)$name[which(all %in% V(av_f)$name)>length(ag)] <- NA
plot(av_f,
     layout = koords2_f,
     vertex.color = V(av_f)$color,
     vertex.size = V(av_f)$size,
     vertex.label = V(av_f)$name, 
     vertex.label.color = "darkgreen", 
     vertex.label.font = 2, 
     vertex.label.cex = 0.75, 
     edge.width = 0.05*sqrt(E(av_f)$weight),
     edge.arrow.size = 0.2,
     edge.curved = TRUE,
     edge.color = gray.colors(1),
     main = "Filtered Nepal Agency-VDC Aid Relief Geo-Network")
legend("topright",
       c("Implementing Aid Agency","VDC With Geo-Coordinates"),
       fill = c("green","SkyBlue2"),
       bty = "n")


# EDGE-FILTRATION BY EDGE WEIGHT OF THE AGENCY-VDC AID NETWORK: CUT-OFF = 75% PERCENTILE
for (k in 1:dim(aid_m)[1]){
  if(k-1<length(ag)){
    V(av)$size[k] <- 3
    V(av)$name[k] <- ag[k]
  } else {
    V(av)$size[k] <-2
    V(av)$name[k] <- NA}
}
cut75 <- quantile(as.vector(aid_m[aid_m>0]),0.75)
av_f <- filter(cutoff = cut75,
             edge_matrix = aid_m,
             vertex_colors = V(av)$color,
             vertex_names = all,
             vertex_size = V(av)$size)

# DISPLAY THE EDGE-FILTERED GRAPH
plot(av_f,
     layout = layout.fruchterman.reingold(av_f,
                                          niter = 200,
                                          area = 2000*vcount(av_f)),
     vertex.color = V(av_f)$color,
     vertex.size = V(av_f)$size,
     vertex.label = NA, 
     vertex.label.color = "black", 
     vertex.label.font = 2, 
     vertex.label.cex = 0.7, 
     edge.width = 0.7*sqrt(E(av_f)$weight),
     edge.arrow.size = 0.3,
     edge.curved = TRUE,
     edge.color = gray.colors(1),
     main = "Filtered Nepal Agency-VDC Aid Relief Abstract Network")
legend("topright",
       c("Implementing Aid Agency","Aid Target VDC"),
       fill = c("green","SkyBlue2"),
       bty = "n")

# PLOT THE AGENCY-VDC AID GEO-NETWORK
av <- graph.adjacency(aid_m,
                      mode = "directed",
                      weighted = TRUE)
V(av)$color <- rep("green",length(all))
for (k in 1:length(all)){
  if(is.element(all[k],vd)){
    V(av)$color[k] <- "SkyBlue2"
  }  
}
for (k in 1:dim(aid_m)[1]){
  if(k-1<length(ag)){
    V(av)$size[k] <- 3
    V(av)$name[k] <- all[k]
  } else {
    V(av)$size[k] <- 1
    V(av)$name[k] <- all[k]}
}
cut75 <- quantile(as.vector(aid_m[aid_m>0]),0.75)
av_f <- filter(cutoff = cut75,
               edge_matrix = aid_m,
               vertex_colors = V(av)$color,
               vertex_names = V(av)$name,
               vertex_size = V(av)$size)
koords2_f <- koords2[which(V(av)$name %in% V(av_f)$name),]
V(av_f)$name[which(all %in% V(av_f)$name)>length(ag)] <- NA
plot(av_f,
     layout = koords2_f,
     vertex.color = V(av_f)$color,
     vertex.size = V(av_f)$size,
     vertex.label = V(av_f)$name, 
     vertex.label.color = "darkgreen", 
     vertex.label.font = 2, 
     vertex.label.cex = 0.75, 
     edge.width = 0.05*sqrt(E(av_f)$weight),
     edge.arrow.size = 0.2,
     edge.curved = TRUE,
     edge.color = gray.colors(1),
     main = "Filtered Nepal Agency-VDC Aid Relief Geo-Network")
legend("topright",
       c("Implementing Aid Agency","VDC With Geo-Coordinates"),
       fill = c("green","SkyBlue2"),
       bty = "n")


# EDGE-FILTRATION BY EDGE WEIGHT OF THE AGENCY-VDC AID NETWORK: CUT-OFF = 85% PERCENTILE
for (k in 1:dim(aid_m)[1]){
  if(k-1<length(ag)){
    V(av)$size[k] <- 4
    V(av)$name[k] <- ag[k]
  } else {
    V(av)$size[k] <-2
    V(av)$name[k] <- NA}
}
cut85 <- quantile(as.vector(aid_m[aid_m>0]),0.85)
av_f <- filter(cutoff = cut85,
               edge_matrix = aid_m,
               vertex_colors = V(av)$color,
               vertex_names = V(av)$name,
               vertex_size = V(av)$size)

# DISPLAY THE EDGE-FILTERED GRAPH
plot(av_f,
     layout = layout.fruchterman.reingold(av_f,
                                          niter = 200,
                                          area = 2000*vcount(av_f)),
     vertex.color = V(av_f)$color,
     vertex.size = V(av_f)$size,
     vertex.label = V(av_f)$name, 
     vertex.label.color = "black", 
     vertex.label.font = 2, 
     vertex.label.cex = 0.7, 
     edge.width = 0.7*sqrt(E(av_f)$weight),
     edge.arrow.size = 0.3,
     edge.curved = TRUE,
     edge.color = gray.colors(1),
     main = "Filtered Nepal Agency-VDC Aid Relief Abstract Network")
legend("topright",
       c("Implementing Aid Agency","Aid Target VDC"),
       fill = c("green","SkyBlue2"),
       bty = "n")

# PLOT THE AGENCY-VDC AID GEO-NETWORK
av <- graph.adjacency(aid_m,
                      mode = "directed",
                      weighted = TRUE)
V(av)$color <- rep("green",length(all))
for (k in 1:length(all)){
  if(is.element(all[k],vd)){
    V(av)$color[k] <- "SkyBlue2"
  }  
}
for (k in 1:dim(aid_m)[1]){
  if(k-1<length(ag)){
    V(av)$size[k] <- 4
    V(av)$name[k] <- all[k]
  } else {
    V(av)$size[k] <- 1
    V(av)$name[k] <- all[k]}
}
cut85 <- quantile(as.vector(aid_m[aid_m>0]),0.85)
av_f <- filter(cutoff = cut85,
               edge_matrix = aid_m,
               vertex_colors = V(av)$color,
               vertex_names = V(av)$name,
               vertex_size = V(av)$size)
koords2_f <- koords2[which(V(av)$name %in% V(av_f)$name),]
V(av_f)$name[which(all %in% V(av_f)$name)>length(ag)] <- NA
plot(av_f,
     layout = koords2_f,
     vertex.color = V(av_f)$color,
     vertex.size = V(av_f)$size,
     vertex.label = V(av_f)$name, 
     vertex.label.color = "darkgreen", 
     vertex.label.font = 2, 
     vertex.label.cex = 0.75, 
     edge.width = 0.1*sqrt(E(av_f)$weight),
     edge.arrow.size = 0.2,
     edge.curved = TRUE,
     edge.color = gray.colors(1),
     main = "Filtered Nepal Agency-VDC Aid Relief Geo-Network")
legend("topright",
       c("Implementing Aid Agency","VDC With Geo-Coordinates"),
       fill = c("green","SkyBlue2"),
       bty = "n")


# DISPLAY THE LARGEST CLUSTER (GIANT COMPONENT) FOR THE
# FILTERED AGENCY-VDC AID NETWORK WITH 
# FILTER THRESHOLD = 75% PERCENTILE
av <- graph.adjacency(aid_m,
                      mode = "directed",
                      weighted = TRUE)
V(av)$color <- rep("green",length(all))
for (k in 1:length(all)){
  if(is.element(all[k],vd)){
    V(av)$color[k] <- "SkyBlue2"
  }  
}
for (k in 1:dim(aid_m)[1]){
  if(k-1<length(ag)){
    V(av)$size[k] <- 3
    V(av)$name[k] <- all[k]
  } else {
    V(av)$size[k] <- 2
    V(av)$name[k] <- all[k]}
}
cut75 <- quantile(as.vector(aid_m[aid_m>0]),0.75)
av_f <- filter(cutoff = cut75,
               edge_matrix = aid_m,
               vertex_colors = V(av)$color,
               vertex_names = V(av)$name,
               vertex_size = V(av)$size)
koords2_f <- koords2[which(V(av)$name %in% V(av_f)$name),]

av_f_c <- giant_comp(graph = av_f,
                     vertex_colors = V(av_f)$color,
                     vertex_names = V(av_f)$name,
                     vertex_size = V(av_f)$size)
koords2_f_c <- koords2[which(V(av)$name %in% V(av_f_c)$name),]
V(av_f_c)$name[which(all %in% V(av_f_c)$name)>length(ag)] <- NA
# DISPLAY THE EDGE-FILTERED GRAPH AGAIN
plot(av_f_c,
     layout = layout.fruchterman.reingold(av_f_c,
                                          niter = 200,
                                          area = 2000*vcount(av_f_c)),
     vertex.color = V(av_f_c)$color,
     vertex.size = V(av_f_c)$size,
     vertex.label = V(av_f_c)$name, 
     vertex.label.color = "darkgreen", 
     vertex.label.font = 1, 
     vertex.label.cex = 1, 
     edge.width = 0.3*sqrt(E(av_f_c)$weight),
     edge.arrow.size = 0.3,
     edge.curved = TRUE,
     edge.color = gray.colors(1),
     main = "Giant Component for Filtered Nepal Agency-VDC Aid Relief Network")
legend("topright",
       c("Implementing Aid Agency","VDC With Geo-Coordinates"),
       fill = c("green","SkyBlue2"),
       bty = "n")

# DISPLAY THE GEO-NETWORK GIANT COMPONENT
plot(av_f_c,
     layout = koords2_f_c,
     vertex.color = V(av_f_c)$color,
     vertex.size = V(av_f_c)$size,
     vertex.label = V(av_f_c)$name, 
     vertex.label.color = "darkgreen", 
     vertex.label.font = 1, 
     vertex.label.cex = 1, 
     edge.width = 0.1*sqrt(E(av_f_c)$weight),
     edge.arrow.size = 0.3,
     edge.curved = TRUE,
     edge.color = gray.colors(1),
     main = "Giant Component for Filtered Nepal Agency-VDC Aid Relief Network")
legend("topright",
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
# ANALYSIS OF THE AGENCY-VDC AID NETWORK ITSELF
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
  if(is.element(all[k],vd)){
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
# AUTHORITY SCORE: THIS IS A MEASURE FOR DIRECTED NETWORKS, THE NUMBER OF NODES THAT ARE
# HUBS AND POINT TO A GIVEN NODE. IT IS DEFIEND AS THE PRINCIPAL EIGENVECTOR FOR t(A)*A, 
# WHERE A SANDS FOR THE ADJACENCY MATRIX OF THE NETWORK. FOR UNDIRECTED NETWORKS, 
# AUTHORITY SCORE IS EQUAL TO THE HUB SCORE
#
# NOTE: THIS IS LIKELY NOT THE RIGHT TOOL TO APPLY HERE aT THE VDC LEVEL
# BUT IF WE OBTAIN MORE GRANUALR DATA, WE WILL BE ABLE TO USE IT
#
#
#
av <- graph.adjacency(aid_m,
                      mode = "directed",
                      weighted = TRUE)
V(av)$color <- rep("green",length(all))
for (k in 1:length(all)){
  if(is.element(all[k],vd)){
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
au <- authority.score(graph = av,
                      weights = E(av)$weight)$vector
au <- 100*sqrt(au)
plot(sort(au, decreasing = TRUE),
     col = adjustcolor(rgb(1,0,1,1)), 
     xlab = "Node Index", 
     ylab = "Authority Score", 
     main = "Sorted Agency-VDC Aid Network Authority Score Values", 
     pch = 19)
histP1(au,
       breaks = 50,
       col = adjustcolor(rgb(1,0,1,1)),
       xlab = "Authority Score Values",
       main = "Agency-VDC Aid Network Authority Score Distribution")


# PLOT HEAT MAP ON VERTICES ACCORDING TO AUTHORITY SCORE
for (k in 1:length(au)){
  if (k>length(ag)){
    V(av)$color[k] <- rev(heat.colors(1+max(as.integer(au))))[1+as.integer(au[k])]
  } else {
    V(av)$color[k] <- "green"
  }
}

plot(av,
     layout = layout.fruchterman.reingold(av, 
                                          niter = 200, 
                                          area = 2000*vcount(av)),
     vertex.color = V(av)$color,
     vertex.size = V(av)$size,
     vertex.label = NA, 
     vertex.label.color = "black",
     vertex.label.font = 1, 
     vertex.label.cex = 0.85, 
     edge.width = 0.2*sqrt(E(av)$weight),
     edge.arrow.size = 0.2,
     edge.curved = TRUE,
     edge.color = gray.colors(1),
     main = "Weighted Agency-VDC Aid Network Authority Score")
legend("topleft",
       c("Highest Authority Score","Lowest Authority Score"),
       fill = c("red","White"),
       bty = "n")

# PLOT THE AUTHORITY SCORE HEAT MAP OF THE AGENCY-VDC AID NETWORK
for (k in 1:dim(aid_m)[1]){
  if(k-1<length(ag)){
    V(av)$size[k] <- 3
    V(av)$name[k] <- all[k]
  } else {
    V(av)$size[k] <- 2
    V(av)$name[k] <- all[k]}
}
av_coords <- koords2[which(all %in% V(av)$name),]

plot(av,
     layout = av_coords,
     vertex.color = V(av)$color,
     vertex.size = V(av)$size,
     vertex.label = NA, 
     vertex.label.color = "black",
     vertex.label.font = 1, 
     vertex.label.cex = 0.75, 
     edge.width = 0.1*sqrt(E(av)$weight),
     edge.arrow.size = 0.2,
     edge.curved = TRUE,
     edge.color = gray.colors(1),
     main = "Weighted Agency-VDC Aid Geo-Network Authority Scores")
legend("topright",
       c("Highest Authority Score","Lowest Authority Score"),
       fill = c("red","White"),
       bty = "n")


#
#
#
# HUB SCORE: MEASURE THE NUMBER OF AUTHORITY NODES THAT A GIVEN
# HUB NODE POINTS TO
#
# NOTE: THIS IS LIKELY NOT THE RIGHT TOOL TO APPLY HERE aT THE VDC LEVEL
# BUT IF WE OBTAIN MORE GRANUALR DATA, WE WILL BE ABLE TO USE IT
#
#
#

av <- graph.adjacency(aid_m,
                      mode = "directed",
                      weighted = TRUE)
V(av)$color <- rep("green",length(all))
for (k in 1:length(all)){
  if(is.element(all[k],vd)){
    V(av)$color[k] <- "SkyBlue2"
  }  
}
for (k in 1:dim(aid_m)[1]){
  if(k-1<length(ag)){
    V(av)$size[k] <- 4
    V(av)$name[k] <- ag[k]
  } else {
    V(av)$size[k] <- 2
    V(av)$name[k] <- NA}
}
hb <- hub.score(graph = av,
                      weights = E(av)$weight)$vector
hb <- 100*(hb)^(1/3)
plot(sort(hb),
     col = adjustcolor(rgb(1,0,1,1)), 
     xlab = "Node Index", 
     ylab = "Hub Score", 
     main = "Sorted Agency-VDC Aid Network Hub Score Values", 
     pch = 19)

histP1(hb,
       breaks = 50,
       col = adjustcolor(rgb(1,0,1,1)),
       xlab = "Hub Score Values",
       main = "Agency-VDC Aid Network Hub Score Distribution")


# PLOT HEAT MAP ON VERTICES ACCORDING TO HUB SCORE
for (k in 1:length(hb)){
  if (k<length(ag)+1){
    V(av)$color[k] <- rev(heat.colors(1+max(as.integer(hb))))[1+as.integer(hb[k])]
  } else {
    V(av)$color[k] <- "SkyBlue2"
  }
}

plot(av,
     layout = layout.fruchterman.reingold(av, 
                                          niter = 200, 
                                          area = 2000*vcount(av)),
     vertex.color = V(av)$color,
     vertex.size = V(av)$size,
     vertex.label = V(av)$name, 
     vertex.label.color = "black",
     vertex.label.font = 1, 
     vertex.label.cex = 0.5, 
     edge.width = 0.1*sqrt(E(av)$weight),
     edge.arrow.size = 0.2,
     edge.curved = TRUE,
     edge.color = gray.colors(1),
     main = "Weighted Agency-VDC Aid Network Hub Score")
legend("topleft",
       c("Highest Hub Score","Lowest Hub Score"),
       fill = c("red","White"),
       bty = "n")

# PLOT THE HUB SCORE HEAT MAP OF THE AGENCY-VDC AID NETWORK
for (k in 1:dim(aid_m)[1]){
  if(k-1<length(ag)){
    V(av)$size[k] <- 5
    V(av)$name[k] <- all[k]
  } else {
    V(av)$size[k] <- 2
    V(av)$name[k] <- all[k]}
}
hb_coords <- koords2[which(all %in% V(av)$name),]
for (k in 1:dim(aid_m)[1]){
  if(k-1<length(ag)){
    V(av)$size[k] <- 6
    V(av)$name[k] <- ag[k]
  } else {
    V(av)$size[k] <- 1
    V(av)$name[k] <- NA}
}
plot(av,
     layout = hb_coords,
     vertex.color = V(av)$color,
     vertex.size = V(av)$size,
     vertex.label = V(av)$name, 
     vertex.label.color = "darkgreen",
     vertex.label.font = 1, 
     vertex.label.cex = 0.5, 
     edge.width = 0.05*sqrt(E(av)$weight),
     edge.arrow.size = 0.2,
     edge.curved = TRUE,
     edge.color = gray.colors(1),
     main = "Weighted Agency-VDC Aid Geo-Network Hub Scores")
legend("topright",
       c("Highest Hub Score","Lowest Hub Score"),
       fill = c("red","White"),
       bty = "n")



#
#
# COMMUNITY STRUCTURES FOR THE AGENCY-VDC AID NETWORK
#
#

#
#
#
# SPINGLASS COMMUNITY DETECTION
#
#
# Spinglass community detection aims to find 
# communities using a spin-glass model 
# and simulated annealing. 
# Edge directions are ignored.
#
#
#


# DEFINE THE WEIGHTED DISPLACEMENT GRAPH
av <- graph.adjacency(aid_m,mode="directed",weighted=TRUE)

# COLOR VERTICES REPRESENTING AGENCIES (GREEN) AND VDCs (BLUE) WHERE AID WAS SENT
V(av)$color<-rep("green",length(all))
for (k in 1:length(all)){
  if(is.element(all[k],vd)){
    V(av)$color[k]<-"SkyBlue2"
  }  
}

# SET UP THE VERTEX SIZES
for (k in 1:dim(aid_m)[1]){
  if(k-1<length(ag)){
    V(av)$size[k] <-3
    V(av)$name[k] <- ag[k]
  } else {
    V(av)$size[k] <-2
    V(av)$name[k] <- NA}
}


# DEfine the SPINGLASS COMMUNITY STRUCTURE
sp <- spinglass.community(graph = av,
                          weights = E(av)$weights,
                          spins = 20)

# PLOT THE COMMUNITIES OF THE Agency-VDC Aid NETWORK
plot(sp,
     av,
     layout = layout.fruchterman.reingold(av, 
                                          niter = 200, 
                                          area = 2000*vcount(av)),
     vertex.color = sp,
     vertex.size = 2,
     vertex.label = NA, 
     vertex.label.color = "black",
     vertex.label.font = 1, 
     vertex.label.cex = 0.85, 
     edge.width = 0.3*sqrt(E(av)$weight),
     edge.arrow.size = 0.3,
     edge.curved = TRUE,
     edge.color = gray.colors(1),
     main = "Weighted Agency-VDC Aid Network Spinglass Communities")

# PLOT THE SPINGLASS COMMUNITIES FOR THE GEO-NETWORK
plot(sp,
     av,
     layout = koords2,
     vertex.color = sp,
     vertex.size = 2,
     vertex.label = NA, 
     vertex.label.color = "black",
     vertex.label.font = 1, 
     vertex.label.cex = 0.85, 
     edge.width = 0.3*sqrt(E(av)$weight),
     edge.arrow.size = 0.3,
     edge.curved = TRUE,
     edge.color = gray.colors(1),
     main = "Weighted Agency-VDC Aid Network Spinglass Communities")


# SPINGLASS COMMUNITIES AND FILTER AT CUTOFF = 90%
cut90 <- quantile(as.vector(aid_m[aid_m>0]),0.90)
av_f <- filter(cutoff = cut90,
               edge_matrix = aid_m,
               vertex_colors = V(av)$color,
               vertex_names = all,
               vertex_size = V(av)$size)
av_f_c <- giant_comp(graph = av_f,
                     vertex_colors = V(av_f)$color,
                     vertex_names = V(av_f)$name,
                     vertex_size = V(av)$size)

# DEFINE THE SPINGLASS COMMUNITY STRUCTURE
sp_f_c <- spinglass.community(graph = av_f_c,
                            weights = E(av_f_c)$weights,
                            spins = 20)

# PLOT THE COMMUNITIES OF THE Agency-VDC Aid NETWORK
plot(sp_f_c,
     av_f_c,
     layout = layout.fruchterman.reingold(av_f_c, 
                                          niter = 200, 
                                          area = 2000*vcount(av_f_c)),
     vertex.color = sp_f_c,
     vertex.size = 2,
     vertex.label = NA, 
     vertex.label.color = "black",
     vertex.label.font = 1, 
     vertex.label.cex = 0.85, 
     edge.width = 0.3*sqrt(E(av_f_c)$weight),
     edge.arrow.size = 0.3,
     edge.curved = TRUE,
     edge.color = gray.colors(1),
     main = "Weighted Agency-VDC Aid Network Spinglass Communities")


# TO PLOT THE FILTERED GEO-NETWORK, RE-DEFINE THE WEIGHTED DISPLACEMENT GRAPH
av <- graph.adjacency(aid_m,mode="directed",weighted=TRUE)

# COLOR VERTICES REPRESENTING AGENCIES (GREEN) AND VDCs (BLUE) WHERE AID WAS SENT
V(av)$color<-rep("green",length(all))
for (k in 1:length(all)){
  if(is.element(all[k],vd)){
    V(av)$color[k]<-"SkyBlue2"
  }  
}

# SET UP THE VERTEX NAMES, SIZES, COORDINATES
for (k in 1:dim(aid_m)[1]){
  if(k-1<length(ag)){
    V(av)$size[k] <-3
    V(av)$name[k] <- all[k]
  } else {
    V(av)$size[k] <-2
    V(av)$name[k] <- all[k]}
}

koords2_f_c <- koords2[which(V(av)$name %in% V(av_f_c)$name),]

# PLOT THE SPINGLASS COMMUNITIES FOR THE GEO-NETWORK
plot(sp_f_c,
     av_f_c,
     layout = koords2_f_c,
     vertex.color = sp_f_c,
     vertex.size = 2,
     vertex.label = NA, 
     vertex.label.color = "black",
     vertex.label.font = 1, 
     vertex.label.cex = 0.85, 
     edge.width = 0.3*sqrt(E(av_f_c)$weight),
     edge.arrow.size = 0.3,
     edge.curved = TRUE,
     edge.color = gray.colors(1),
     main = "Weighted Agency-VDC Aid Network Spinglass Communities")

# SOME BASIC SPINGLASS COMMUNITY STATS
plot(sort(sp$membership), 
     col = adjustcolor(rgb(1,0,1,1)), 
     xlab = "Node Index in the Network", 
     ylab = "Spinglass Community Values", 
     main = "Sorted Spinglass Community Values for Nepal Agency-VDC Aid Network", 
     pch = 19)

histP1(sp$membership,
       breaks = 60,
       col = adjustcolor(rgb(1,0,1,1)),
       xlab = "Spinglass Community Values",
       main = "Spinglass Community Distribution for Nepal Displacement Network")



#
#
#
# WALKTRAP COMMUNITY DETECTION
#
#
# This is an approach based on random walks. The general idea is that if you perform random walks on the graph, 
# then the walks are more likely to stay within the same community because there are only a few edges that lead 
# outside a given community. Walktrap runs short random walks of 3-4-5 steps (depending on one of its parameters) 
# and uses the results of these random walks to merge separate communities in a bottom-up manner like fastgreedy.community. 
# Again, you can use the modularity score to select where to cut the dendrogram. It is a bit slower than the fast greedy 
# approach but also a bit more accurate (according to the original publication).
#
#
#


# DEFINE THE WEIGHTED DISPLACEMENT GRAPH
av <- graph.adjacency(aid_m,mode="directed",weighted=TRUE)

# COLOR VERTICES REPRESENTING AGENCIES (GREEN) AND VDCs (BLUE) WHERE AID WAS SENT
V(av)$color<-rep("green",length(all))
for (k in 1:length(all)){
  if(is.element(all[k],vd)){
    V(av)$color[k]<-"SkyBlue2"
  }  
}

# SET UP THE VERTEX SIZES
for (k in 1:dim(aid_m)[1]){
  if(k-1<length(ag)){
    V(av)$size[k] <-3
    V(av)$name[k] <- ag[k]
  } else {
    V(av)$size[k] <-2
    V(av)$name[k] <- NA}
}


# DEFINE THE WALKTRAP COMMUNITY STRUCTURE
wk <- walktrap.community(graph = av,
                           membership = TRUE,
                           weights = E(av)$weights)

# PLOT THE COMMUNITIES OF THE Agency-VDC Aid NETWORK
plot(wk,
     av,
     layout = layout.fruchterman.reingold(av, 
                                          niter = 200, 
                                          area = 2000*vcount(av)),
     vertex.color = wk,
     vertex.size = V(av)$size,
     vertex.label = NA, 
     vertex.label.color = "black",
     vertex.label.font = 1, 
     vertex.label.cex = 0.85, 
     edge.width = 0.3*sqrt(E(av)$weight),
     edge.arrow.size = 0.3,
     edge.curved = TRUE,
     edge.color = gray.colors(1),
     main = "Weighted Agency-VDC Aid Network Fastgreedy Communities")

# PLOT THE WALKTRAP COMMUNITIES FOR THE GEO-NETWORK
plot(wk,
     av,
     layout = koords2,
     vertex.color = wk,
     vertex.size = V(av)$size,
     vertex.label = NA, 
     vertex.label.color = "black",
     vertex.label.font = 1, 
     vertex.label.cex = 0.85, 
     edge.width = 0.3*sqrt(E(av)$weight),
     edge.arrow.size = 0.3,
     edge.curved = TRUE,
     edge.color = gray.colors(1),
     main = "Weighted Agency-VDC Aid Network Walktrap Communities")


# WALKTRAP COMMUNITIES AND FILTER AT CUTOFF = 70%
cut70 <- quantile(as.vector(aid_m[aid_m>0]),0.70)
av_f <- filter(cutoff = cut70,
               edge_matrix = aid_m,
               vertex_colors = V(av)$color,
               vertex_names = all)
av_f_c <- giant_comp(graph = av_f,
                     vertex_colors = V(av_f)$color,
                     vertex_names = V(av_f)$name)

# DEFINE THE WALKTRAP COMMUNITY STRUCTURE
wk_f_c <- walktrap.community(graph = av_f_c,
                              weights = E(av_f_c)$weights)

# PLOT THE COMMUNITIES OF THE Agency-VDC Aid NETWORK
plot(wk_f_c,
     av_f_c,
     layout = layout.fruchterman.reingold(av_f_c, 
                                          niter = 200, 
                                          area = 2000*vcount(av_f_c)),
     vertex.color = wk_f_c,
     vertex.size = 2,
     vertex.label = NA, 
     vertex.label.color = "black",
     vertex.label.font = 1, 
     vertex.label.cex = 0.85, 
     edge.width = 0.3*sqrt(E(av_f_c)$weight),
     edge.arrow.size = 0.3,
     edge.curved = TRUE,
     edge.color = gray.colors(1),
     main = "Weighted Agency-VDC Aid Network Walktrap Communities")


# TO PLOT THE FILTERED GEO-NETWORK, RE-DEFINE THE WEIGHTED DISPLACEMENT GRAPH
av <- graph.adjacency(aid_m,mode="directed",weighted=TRUE)

# COLOR VERTICES REPRESENTING AGENCIES (GREEN) AND VDCs (BLUE) WHERE AID WAS SENT
V(av)$color<-rep("green",length(all))
for (k in 1:length(all)){
  if(is.element(all[k],vd)){
    V(av)$color[k]<-"SkyBlue2"
  }  
}

# SET UP THE VERTEX NAMES, SIZES, COORDINATES
for (k in 1:dim(aid_m)[1]){
  if(k-1<length(ag)){
    V(av)$size[k] <-3
    V(av)$name[k] <- all[k]
  } else {
    V(av)$size[k] <-2
    V(av)$name[k] <- all[k]}
}

koords2_f_c <- koords2[which(V(av)$name %in% V(av_f_c)$name),]

# PLOT THE WALKTRAP COMMUNITIES FOR THE GEO-NETWORK
plot(wk_f_c,
     av_f_c,
     layout = koords2_f_c,
     vertex.color = wk_f_c,
     vertex.size = 2,
     vertex.label = NA, 
     vertex.label.color = "black",
     vertex.label.font = 1, 
     vertex.label.cex = 0.85, 
     edge.width = 0.3*sqrt(E(av_f_c)$weight),
     edge.arrow.size = 0.3,
     edge.curved = TRUE,
     edge.color = gray.colors(1),
     main = "Weighted Agency-VDC Aid Network Walktrap Communities")


# WALKTRAP COMMUNITIES AND FILTER AT CUTOFF = 80%
cut80 <- quantile(as.vector(aid_m[aid_m>0]),0.80)
av_f <- filter(cutoff = cut80,
               edge_matrix = aid_m,
               vertex_colors = V(av)$color,
               vertex_names = all)
av_f_c <- giant_comp(graph = av_f,
                     vertex_colors = V(av_f)$color,
                     vertex_names = V(av_f)$name)

# DEFINE THE WALKTRAP COMMUNITY STRUCTURE
wk_f_c <- walktrap.community(graph = av_f_c,
                             weights = E(av_f_c)$weights)

# PLOT THE COMMUNITIES OF THE Agency-VDC Aid NETWORK
plot(wk_f_c,
     av_f_c,
     layout = layout.fruchterman.reingold(av_f_c, 
                                          niter = 200, 
                                          area = 2000*vcount(av_f_c)),
     vertex.color = wk_f_c,
     vertex.size = 2,
     vertex.label = NA, 
     vertex.label.color = "black",
     vertex.label.font = 1, 
     vertex.label.cex = 0.85, 
     edge.width = 0.3*sqrt(E(av_f_c)$weight),
     edge.arrow.size = 0.3,
     edge.curved = TRUE,
     edge.color = gray.colors(1),
     main = "Weighted Agency-VDC Aid Network Walktrap Communities")


# TO PLOT THE FILTERED GEO-NETWORK, RE-DEFINE THE WEIGHTED DISPLACEMENT GRAPH
av <- graph.adjacency(aid_m,mode="directed",weighted=TRUE)

# COLOR VERTICES REPRESENTING AGENCIES (GREEN) AND VDCs (BLUE) WHERE AID WAS SENT
V(av)$color<-rep("green",length(all))
for (k in 1:length(all)){
  if(is.element(all[k],vd)){
    V(av)$color[k]<-"SkyBlue2"
  }  
}

# SET UP THE VERTEX NAMES, SIZES, COORDINATES
for (k in 1:dim(aid_m)[1]){
  if(k-1<length(ag)){
    V(av)$size[k] <-3
    V(av)$name[k] <- all[k]
  } else {
    V(av)$size[k] <-2
    V(av)$name[k] <- all[k]}
}

koords2_f_c <- koords2[which(V(av)$name %in% V(av_f_c)$name),]

# PLOT THE WALKTRAP COMMUNITIES FOR THE GEO-NETWORK
plot(wk_f_c,
     av_f_c,
     layout = koords2_f_c,
     vertex.color = wk_f_c,
     vertex.size = 2,
     vertex.label = NA, 
     vertex.label.color = "black",
     vertex.label.font = 1, 
     vertex.label.cex = 0.85, 
     edge.width = 0.3*sqrt(E(av_f_c)$weight),
     edge.arrow.size = 0.3,
     edge.curved = TRUE,
     edge.color = gray.colors(1),
     main = "Weighted Agency-VDC Aid Network Walktrap Communities")


# WALKTRAP COMMUNITIES AND FILTER AT CUTOFF = 90%
cut90 <- quantile(as.vector(aid_m[aid_m>0]),0.90)
av_f <- filter(cutoff = cut90,
               edge_matrix = aid_m,
               vertex_colors = V(av)$color,
               vertex_names = all)
av_f_c <- giant_comp(graph = av_f,
                     vertex_colors = V(av_f)$color,
                     vertex_names = V(av_f)$name)

# DEFINE THE WALKTRAP COMMUNITY STRUCTURE
wk_f_c <- walktrap.community(graph = av_f_c,
                             weights = E(av_f_c)$weights)

# PLOT THE COMMUNITIES OF THE Agency-VDC Aid NETWORK
plot(wk_f_c,
     av_f_c,
     layout = layout.fruchterman.reingold(av_f_c, 
                                          niter = 200, 
                                          area = 2000*vcount(av_f_c)),
     vertex.color = wk_f_c,
     vertex.size = 2,
     vertex.label = NA, 
     vertex.label.color = "black",
     vertex.label.font = 1, 
     vertex.label.cex = 0.85, 
     edge.width = 0.3*sqrt(E(av_f_c)$weight),
     edge.arrow.size = 0.3,
     edge.curved = TRUE,
     edge.color = gray.colors(1),
     main = "Weighted Agency-VDC Aid Network Walktrap Communities")


# TO PLOT THE FILTERED GEO-NETWORK, RE-DEFINE THE WEIGHTED DISPLACEMENT GRAPH
av <- graph.adjacency(aid_m,mode="directed",weighted=TRUE)

# COLOR VERTICES REPRESENTING AGENCIES (GREEN) AND VDCs (BLUE) WHERE AID WAS SENT
V(av)$color<-rep("green",length(all))
for (k in 1:length(all)){
  if(is.element(all[k],vd)){
    V(av)$color[k]<-"SkyBlue2"
  }  
}

# SET UP THE VERTEX NAMES, SIZES, COORDINATES
for (k in 1:dim(aid_m)[1]){
  if(k-1<length(ag)){
    V(av)$size[k] <-3
    V(av)$name[k] <- all[k]
  } else {
    V(av)$size[k] <-2
    V(av)$name[k] <- all[k]}
}

koords2_f_c <- koords2[which(V(av)$name %in% V(av_f_c)$name),]

# PLOT THE WALKTRAP COMMUNITIES FOR THE GEO-NETWORK
plot(wk_f_c,
     av_f_c,
     layout = koords2_f_c,
     vertex.color = wk_f_c,
     vertex.size = 2,
     vertex.label = NA, 
     vertex.label.color = "black",
     vertex.label.font = 1, 
     vertex.label.cex = 0.85, 
     edge.width = 0.3*sqrt(E(av_f_c)$weight),
     edge.arrow.size = 0.3,
     edge.curved = TRUE,
     edge.color = gray.colors(1),
     main = "Weighted Agency-VDC Aid Network Walktrap Communities")


# WALKTRAP COMMUNITIES AND FILTER AT CUTOFF = 93%
cut93 <- quantile(as.vector(aid_m[aid_m>0]),0.93)
av_f <- filter(cutoff = cut93,
               edge_matrix = aid_m,
               vertex_colors = V(av)$color,
               vertex_names = all)
av_f_c <- giant_comp(graph = av_f,
                     vertex_colors = V(av_f)$color,
                     vertex_names = V(av_f)$name)

# DEFINE THE WALKTRAP COMMUNITY STRUCTURE
wk_f_c <- walktrap.community(graph = av_f_c,
                             weights = E(av_f_c)$weights)

# PLOT THE COMMUNITIES OF THE Agency-VDC Aid NETWORK
plot(wk_f_c,
     av_f_c,
     layout = layout.fruchterman.reingold(av_f_c, 
                                          niter = 200, 
                                          area = 2000*vcount(av_f_c)),
     vertex.color = wk_f_c,
     vertex.size = 2,
     vertex.label = NA, 
     vertex.label.color = "black",
     vertex.label.font = 1, 
     vertex.label.cex = 0.85, 
     edge.width = 0.3*sqrt(E(av_f_c)$weight),
     edge.arrow.size = 0.3,
     edge.curved = TRUE,
     edge.color = gray.colors(1),
     main = "Weighted Agency-VDC Aid Network Walktrap Communities")


# TO PLOT THE FILTERED GEO-NETWORK, RE-DEFINE THE WEIGHTED DISPLACEMENT GRAPH
av <- graph.adjacency(aid_m,mode="directed",weighted=TRUE)

# COLOR VERTICES REPRESENTING AGENCIES (GREEN) AND VDCs (BLUE) WHERE AID WAS SENT
V(av)$color<-rep("green",length(all))
for (k in 1:length(all)){
  if(is.element(all[k],vd)){
    V(av)$color[k]<-"SkyBlue2"
  }  
}

# SET UP THE VERTEX NAMES, SIZES, COORDINATES
for (k in 1:dim(aid_m)[1]){
  if(k-1<length(ag)){
    V(av)$size[k] <-3
    V(av)$name[k] <- all[k]
  } else {
    V(av)$size[k] <-2
    V(av)$name[k] <- all[k]}
}

koords2_f_c <- koords2[which(V(av)$name %in% V(av_f_c)$name),]

# PLOT THE WALKTRAP COMMUNITIES FOR THE GEO-NETWORK
plot(wk_f_c,
     av_f_c,
     layout = koords2_f_c,
     vertex.color = wk_f_c,
     vertex.size = 2,
     vertex.label = NA, 
     vertex.label.color = "black",
     vertex.label.font = 1, 
     vertex.label.cex = 0.85, 
     edge.width = 0.3*sqrt(E(av_f_c)$weight),
     edge.arrow.size = 0.3,
     edge.curved = TRUE,
     edge.color = gray.colors(1),
     main = "Weighted Agency-VDC Aid Network Walktrap Communities")


# SOME BASIC WALKTRAP COMMUNITY STATS
plot(sort(wk$membership), 
     col = adjustcolor(rgb(1,0,1,1)), 
     xlab = "Node Index in the Network", 
     ylab = "Walktrap Community Values", 
     main = "Sorted Walktrap Community Values for Nepal Agency-VDC Aid Network", 
     pch = 19)

histP1(wk$membership,
       breaks = 60,
       col = adjustcolor(rgb(1,0,1,1)),
       xlab = "Walktrap Community Values",
       main = "Walktrap Community Size Distribution for Nepal Agency-VDC Aid Network")
