# Nepal Disaster Aid Distribution Network Analysis and Projections
# author: Georgi D. Gospodinov
# date: "Augist 15, 2015"
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
# AGENCY-VDC AID NETWORK
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
    V(av)$size[k] <- 2
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
    V(av)$size[k] <- 2
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
    V(av)$size[k] <- 3
    V(av)$name[k] <- ag[k]
  } else {
    V(av)$size[k] <- 2
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
    V(av)$size[k] <- 3
    V(av)$name[k] <- all[k]
  } else {
    V(av)$size[k] <- 2
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
# This is an approach based on random walks. The general idea is that if you perform random walks 
# on the graph,vthen the walks are more likely to stay within the same community because there are 
# only a few edges that lead outside a given community. Walktrap runs short random walks of 3-4-5 
# steps (depending on one of its parameters) and uses the results of these random walks to merge 
# separate communities in a bottom-up manner like fastgreedy.community. Again, you can use the 
# modularity score to select where to cut the dendrogram. It is a bit slower than the fast greedy 
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
    V(av)$size[k] <- 3
    V(av)$name[k] <- ag[k]
  } else {
    V(av)$size[k] <- 2
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
               vertex_names = all,
               vertex_size = V(av)$size)
av_f_c <- giant_comp(graph = av_f,
                     vertex_colors = V(av_f)$color,
                     vertex_names = V(av_f)$name,
                     vertex_size = V(av_f)$size)

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
    V(av)$size[k] <- 3
    V(av)$name[k] <- all[k]
  } else {
    V(av)$size[k] <- 2
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
               vertex_names = all,
               vertex_size = V(av)$size)
av_f_c <- giant_comp(graph = av_f,
                     vertex_colors = V(av_f)$color,
                     vertex_names = V(av_f)$name,
                     vertex_size = V(av_f)$size)

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
    V(av)$size[k] <- 3
    V(av)$name[k] <- all[k]
  } else {
    V(av)$size[k] <- 2
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
               vertex_names = all,
               vertex_size = V(av)$size)
av_f_c <- giant_comp(graph = av_f,
                     vertex_colors = V(av_f)$color,
                     vertex_names = V(av_f)$name,
                     vertex_size = V(av_f)$size)

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
    V(av)$size[k] <- 3
    V(av)$name[k] <- all[k]
  } else {
    V(av)$size[k] <- 2
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
               vertex_names = all,
               vertex_size = V(av)$size)
av_f_c <- giant_comp(graph = av_f,
                     vertex_colors = V(av_f)$color,
                     vertex_names = V(av_f)$name,
                     vertex_size = V(av_f)$size)

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
    V(av)$size[k] <- 3
    V(av)$name[k] <- all[k]
  } else {
    V(av)$size[k] <- 2
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


















# PROJECT EACH GRAPH COMPONENT WITH APPROPRIATE CONNECTIONS



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
                                        niter = 200,
                                        area = 2000*vcount(agg)),
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
                                          niter = 200,
                                          area = 2000*vcount(agg)),
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

# CHANGE THE NODE SIZE TO REFLECT NUBER OF SHARED VDCs USING LOG
V(agg)$size <- log(graph.strength(agg))

# PLOT AGENCY GRAPH AND FILTER
plot(agg,
     layout = layout.fruchterman.reingold(agg, 
                                          niter = 200,
                                          area = 2000*vcount(agg)),
     vertex.color = V(agg)$color,
     vertex.size = V(agg)$size,
     vertex.label = V(agg)$name, 
     vertex.label.color = "darkgreen", 
     vertex.label.font = 1, 
     vertex.label.cex = 0.75, 
     edge.width = 0.25*E(agg)$weight,
     edge.curved = TRUE,
     edge.color = gray.colors(1),
     main = "Aid Agency Association Network (Node Size = Log(Weighted Degree)
(Node Weighted Degree = number of agencies with shared VDC aid targets, 
     weighted by the number of VDCs per agency)")




# BEFORE WE FILTER, MULTILEVEL COMMUNITY DETECTION
mc <- multilevel.community(agg)
plot(mc,
     agg, 
     vertex.size = V(agg)$size,
     edge.width = 0.15*E(agg)$weight,
     vertex.label.cex = 0.8,
     vertex.label = V(agg)$name,
     main = "ML Community Detection for Aid Agency Association Network")

# MULTILEVEL COMMUNITIES BREAKDOWN
as.data.frame(sort(membership(mc)))


# MULTILEVEL COMMUNITY DETECTION WITH 75% FILTRATION
cut75 <- quantile(as.vector(ag_m[ag_m>0]),0.75)
agg_f<-filter(cutoff = cut75,
              edge_matrix = ag_m,
              vertex_colors = V(agg)$color,
              vertex_names = V(agg)$name,
              vertex_size = V(agg)$size)

plot(as.undirected(agg_f),
     layout = layout.fruchterman.reingold(agg_f,
                                          niter = 200,
                                          area = 2000*vcount(agg_f)),
     vertex.color = V(agg_f)$color,
     vertex.size = V(agg_f)$size,
     vertex.label = V(agg_f)$name, 
     vertex.label.color = "darkgreen", 
     vertex.label.font = 1, 
     vertex.label.cex = 1, 
     edge.width = 0.25*(E(agg_f)$weight),
     edge.curved = TRUE,
     edge.color = gray.colors(1),
     main = "75% Level Filtration of the Aid Agency Association Network (Node Size = Log(Weighted Degree))")

mc_f <- multilevel.community(as.undirected(agg_f))
plot(mc_f,
     as.undirected(agg_f), 
     vertex.size = V(agg_f)$size,
     edge.width = 0.15*E(agg_f)$weight,
     vertex.label.cex = 1,
     vertex.label = V(agg_f)$name,
     edge.curved = FALSE,
     main = "ML Communities for the Aid Agency Association Network (75% Filtration)")


# MULTILEVEL COMMUNITY DETECTION WITH 85% FILTRATION
cut85 <- quantile(as.vector(ag_m[ag_m>0]),0.85)
agg_f <- filter(cutoff = cut85,
                edge_matrix = ag_m,
                vertex_colors = V(agg)$color,
                vertex_names = V(agg)$name,
                vertex_size = V(agg)$size)

plot(as.undirected(agg_f),
     layout = layout.fruchterman.reingold(agg_f,
                                          niter = 200,
                                          area = 2000*vcount(agg_f)),
     vertex.color = V(agg_f)$color,
     vertex.size = V(agg_f)$size,
     vertex.label = V(agg_f)$name, 
     vertex.label.color = "darkgreen", 
     vertex.label.font = 1, 
     vertex.label.cex = 1, 
     edge.width = 0.25*(E(agg_f)$weight),
     edge.curved = TRUE,
     edge.color = gray.colors(1),
     main = "ML Communities for the Aid Agency Association Network (85% Filtration)")


mc_f <- multilevel.community(as.undirected(agg_f))
plot(mc_f,
     as.undirected(agg_f), 
     vertex.size = V(agg_f)$size,
     edge.width = 0.15*E(agg_f)$weight,
     vertex.label.cex = 1,
     vertex.label = V(agg_f)$name,
     edge.curved = FALSE,
     main = "ML Communities for the Aid Agency Association Network (85% Filtration)")


# MULTILEVEL COMMUNITY DETECTION WITH 90% FILTRATION
cut90 <- quantile(as.vector(ag_m[ag_m>0]),0.90)
agg_f <- filter(cutoff = cut90,
                edge_matrix = ag_m,
                vertex_colors = V(agg)$color,
                vertex_names = V(agg)$name,
                vertex_size = V(agg)$size)

plot(as.undirected(agg_f),
     layout = layout.fruchterman.reingold(agg_f,
                                          niter = 200,
                                          area = 2000*vcount(agg_f)),
     vertex.color = V(agg_f)$color,
     vertex.size = V(agg_f)$size,
     vertex.label = V(agg_f)$name, 
     vertex.label.color = "darkgreen", 
     vertex.label.font = 1, 
     vertex.label.cex = 1, 
     edge.width = 0.15*(E(agg_f)$weight),
     edge.curved = TRUE,
     edge.color = gray.colors(1),
     main = "ML Communities for the Aid Agency Association Network (90% Filtration)")


mc_f <- multilevel.community(as.undirected(agg_f))
plot(mc_f,
     as.undirected(agg_f), 
     vertex.size = V(agg_f)$size,
     edge.width = 0.15*E(agg_f)$weight,
     vertex.label.cex = 1,
     vertex.label = V(agg_f)$name,
     edge.curved = FALSE,
     main = "ML Communities for the Aid Agency Association Network (90% Filtration)")


# MULTILEVEL COMMUNITY DETECTION WITH 95% FILTRATION
cut95 <- quantile(as.vector(ag_m[ag_m>0]),0.95)
agg_f <- filter(cutoff = cut95,
                edge_matrix = ag_m,
                vertex_colors = V(agg)$color,
                vertex_names = V(agg)$name,
                vertex_size = V(agg)$size)

plot(as.undirected(agg_f),
     layout = layout.fruchterman.reingold(agg_f,
                                          niter = 200,
                                          area = 2000*vcount(agg_f)),
     vertex.color = V(agg_f)$color,
     vertex.size = V(agg_f)$size,
     vertex.label = V(agg_f)$name, 
     vertex.label.color = "darkgreen", 
     vertex.label.font = 1, 
     vertex.label.cex = 1, 
     edge.width = 0.15*(E(agg_f)$weight),
     edge.curved = TRUE,
     edge.color = gray.colors(1),
     main = "ML Communities for the Aid Agency Association Network (95% Filtration)")


mc_f <- multilevel.community(as.undirected(agg_f))
plot(mc_f,
     as.undirected(agg_f), 
     vertex.size = V(agg_f)$size,
     edge.width = 0.15*E(agg_f)$weight,
     vertex.label.cex = 1,
     vertex.label = V(agg_f)$name,
     edge.curved = FALSE,
     main = "ML Communities for the Aid Agency Association Network (95% Filtration)")

# MULTILEVEL COMMUNITY DETECTION WITH 97% FILTRATION
cut97 <- quantile(as.vector(ag_m[ag_m>0]),0.97)
agg_f <- filter(cutoff = cut97,
                edge_matrix = ag_m,
                vertex_colors = V(agg)$color,
                vertex_names = V(agg)$name,
                vertex_size = V(agg)$size)

plot(as.undirected(agg_f),
     layout = layout.fruchterman.reingold(agg_f,
                                          niter = 200,
                                          area = 2000*vcount(agg_f)),
     vertex.color = V(agg_f)$color,
     vertex.size = V(agg_f)$size,
     vertex.label = V(agg_f)$name, 
     vertex.label.color = "darkgreen", 
     vertex.label.font = 1, 
     vertex.label.cex = 1, 
     edge.width = 0.15*(E(agg_f)$weight),
     edge.curved = TRUE,
     edge.color = gray.colors(1),
     main = "ML Communities for the Aid Agency Association Network (97% Filtration)")


mc_f <- multilevel.community(as.undirected(agg_f))
plot(mc_f,
     as.undirected(agg_f), 
     vertex.size = V(agg_f)$size,
     edge.width = 0.15*E(agg_f)$weight,
     vertex.label.cex = 1,
     vertex.label = V(agg_f)$name,
     edge.curved = FALSE,
     main = "ML Communities for the Aid Agency Association Network (97% Filtration)")

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
unique_aid <- unique(cbind.data.frame(aid_data$impl_ag,aid_data$vdc))
colnames(unique_aid) <- c("impl_ag","vdc")
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
                vertex_names = ag)
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
sh<-shortest.paths(agg)
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
     col=adjustcolor(rgb(1,0,1/2,1)),
     xlab="Path Length Values",
     main="Path Length Distribution for g")



# BETWEENNESS CENTRALITY: THE NUMBER OF GEODESICS GOING THROUGH A NODE
bc <- betweenness(agg,v=V(agg), directed=FALSE)
plot(sort(bc, decreasing=TRUE),
     col=adjustcolor(rgb(1/2,0,0,1/2)), 
     xlab="Node Index", 
     ylab="Betweenness Centrality", 
     main="Sorted Relief Agency Betweenness Centrality Values", 
     pch=19)
histP2(bc,
       breaks=100,
       col=adjustcolor(rgb(1/2,0,0,1/2)),
       xlab="Betweenness Centrality Values",
       main="Relief Agency Betweenness Centrality Distribution")

# PLOT HEAT MAP ON VERTICES ACCORDING TO BETWEENNESS CENTRALITY
bc_int <- as.integer(round(bc,0))
for (k in 1:length(bc_int)){
  V(agg)$color[k] <- rev(heat.colors(1+as.integer(max(bc_int))))[as.integer(bc_int[k])+1]
}

# PLOT AGENCY GRAPH AND FILTER
plot(agg,
     layout = layout.fruchterman.reingold(agg, niter=200, area=2000*vcount(agg)),
     vertex.color = V(agg)$color,
     vertex.size = 7,
     vertex.label = V(agg)$name, 
     vertex.label.color = "black",
     vertex.label.font = 1.5, 
     vertex.label.cex = 1, 
     edge.width = 0.5*E(agg)$weight,
     edge.curved = TRUE,
     edge.color = gray.colors(1))

# FIND THE TOP 10% BETWEENNES NODES
top_bc <- bc[which(bc > quantile(bc,0.9))]
top_bc

# FIND THE TOP 5% BETWENNES NODES
top_bc <- bc[which(bc > quantile(bc,0.95))]
top_bc

# REMOVE ISOLATED
agg <- drop_isolated(graph = agg,
                     vertex_colors = V(agg)$color,
                     vertex_names = V(agg)$name)

# GET THE GIANT CONNECTED COMPONENT (TWO CLSUTERS ONLY)
agg <- giant_comp(graph = agg,
                  vertex_colors = V(agg)$color,
                  vertex_names = V(agg)$name)











# SET THE GRAPH COLOR ACCORDING TO BC
agg <- as.undirected(graph.adjacency(ag_m,weighted=TRUE))
# SET THE GRAPH COLOR
V(agg)$color <- rep("green",length(ag))
V(agg)$name <- ag
# REMOVE ISOLATED
agg <- drop_isolated(graph = agg,
                     vertex_colors = V(agg)$color,
                     vertex_names = V(agg)$name)

# GET THE GIANT CONNECTED COMPONENT (TWO CLSUTERS ONLY)
agg <- giant_comp(graph = agg,
                  vertex_colors = V(agg)$color,
                  vertex_names = V(agg)$name)



bc<-betweenness(agg,v=V(agg), directed=FALSE)
bc_int <- as.integer(round(bc,0))
for (k in 1:length(bc_int)){
  V(agg)$color[k] <- rev(heat.colors(1+as.integer(max(bc_int))))[as.integer(bc_int[k])+1]
}
plot(agg,
     layout = layout.fruchterman.reingold(agg, niter=2000, area=20000*vcount(agg)),
     vertex.color = V(agg)$color,
     vertex.size = 5,
     vertex.label = NA, 
     vertex.label.color = "black",
     vertex.label.font = 1, 
     vertex.label.cex = 0.5, 
     edge.width = 0.1*E(agg)$weight,
     edge.curved = TRUE,
     edge.color = gray.colors(1),
     main="Aid Agency Network Betweenness Centrality Heat Map")


# FILTER AND REPEAT:
agg <- as.undirected(graph.adjacency(ag_m,weighted=TRUE))
V(agg)$color <- rep("green",length(ag))
V(agg)$name <- ag
cut75 <- quantile(as.vector(ag_m[ag_m>0]),0.75)
agg_f <- filter(cutoff = cut75,
                edge_matrix = ag_m,
                vertex_colors = V(agg)$color,
                vertex_names = V(agg)$name)
agg_f <- as.undirected(agg_f)
agg_f <- giant_comp(graph = agg_f,
                  vertex_color = V(agg_f)$color,
                  vertex_names = V(agg_f)$name)
bc<-betweenness(agg_f,v=V(agg_f), directed=FALSE)
bc_int <- as.integer(round(bc,0))
for (k in 1:length(bc_int)){
  V(agg_f)$color[k] <- rev(heat.colors(1+as.integer(max(bc_int))))[as.integer(bc_int[k])+1]
}
plot(agg_f,
     layout = layout.fruchterman.reingold(agg_f, niter=200, area=2000*vcount(agg_f)),
     vertex.color = V(agg_f)$color,
     vertex.size = 5,
     vertex.label = NA, 
     vertex.label.color = "black",
     vertex.label.font = 1, 
     vertex.label.cex = 0.5, 
     edge.width = 0.1*E(agg_f)$weight,
     edge.curved = TRUE,
     edge.color = gray.colors(1))

# FILTER AND REPEAT:
agg <- as.undirected(graph.adjacency(ag_m,weighted=TRUE))
V(agg)$color <- rep("green",length(ag))
V(agg)$name <- ag
cut90 <- quantile(as.vector(ag_m[ag_m>0]),0.90)
agg_f <- filter(cutoff = cut90,
                edge_matrix = ag_m,
                vertex_colors = V(agg)$color,
                vertex_names = V(agg)$name)
agg_f <- as.undirected(agg_f)
agg_f <- giant_comp(graph = agg_f,
                    vertex_color = V(agg_f)$color,
                    vertex_names = V(agg_f)$name)
bc<-betweenness(agg_f,v=V(agg_f), directed=FALSE)
bc_int <- as.integer(round(bc,0))
for (k in 1:length(bc_int)){
  V(agg_f)$color[k] <- rev(heat.colors(1+as.integer(max(bc_int))))[as.integer(bc_int[k])+1]
}
plot(agg_f,
     layout = layout.fruchterman.reingold(agg_f, niter=200, area=2000*vcount(agg_f)),
     vertex.color = V(agg_f)$color,
     vertex.size = 5,
     vertex.label = NA, 
     vertex.label.color = "black",
     vertex.label.font = 1, 
     vertex.label.cex = 0.5, 
     edge.width = 0.1*E(agg_f)$weight,
     edge.curved = TRUE,
     edge.color = gray.colors(1))








# EDGE-BETWEENNESS CENTRALITY: THE NUMBER OF GEODESICS GOING THROUGH AN EDGE
agg <- as.undirected(graph.adjacency(ag_m,weighted=TRUE))
V(agg)$color <- rep("green",length(ag))
V(agg)$name <- ag
ec <- edge.betweenness(graph = agg,
                       e = E(agg), 
                       directed = FALSE,
                       weights = E(agg)$weight)
plot(sort(ec, decreasing=TRUE),
     col=adjustcolor(rgb(1/2,0,0,1/2)), 
     xlab="Node Index", 
     ylab="Edge-Betweenness Centrality", 
     main="Sorted Relief Agency Edge-Betweenness Centrality Values", 
     pch=19)
histP2(ec,
       breaks=100,
       col=adjustcolor(rgb(1/2,0,0,1/2)),
       xlab="Edge-Betweenness Centrality Values",
       main="Relief Agency Edge-Betweenness Centrality Distribution")

# PLOT HEAT MAP ON VERTICES ACCORDING TO EDGE-BEWEENNESS CENTRALITY
ec_int <- as.integer(round(ec,0))
for (k in 1:length(ec_int)){
  E(agg)$color[k] <- rev(heat.colors(1+as.integer(max(ec_int))))[as.integer(ec_int[k])+1]
}

# PLOT AGENCY GRAPH AND FILTER
plot(agg,
     layout = layout.fruchterman.reingold(agg, niter=200, area=2000*vcount(agg)),
     vertex.color = "green",
     vertex.size = 7,
     vertex.label = V(agg)$name, 
     vertex.label.color = "black",
     vertex.label.font = 1.5, 
     vertex.label.cex = 1, 
     edge.width = 0.5*E(agg)$weight,
     edge.curved = TRUE,
     edge.color = E(agg)$color)

# FIND THE TOP 10% EDGE-BETWEENNES EDGES
top_ec <- ec[which(ec > quantile(ec,0.9))]
E(agg)[which(ec %in% top_ec)]

# FIND THE TOP 5% EDGE-BETWEENNES EDGES
top_ec <- ec[which(ec > quantile(ec,0.99))]
E(agg)[which(ec %in% top_ec)]

# REMOVE ISOLATED
agg <- drop_isolated(graph = agg,
                     vertex_colors = V(agg)$color,
                     vertex_names = V(agg)$name)

# GET THE GIANT CONNECTED COMPONENT (TWO CLSUTERS ONLY)
agg <- giant_comp(graph = agg,
                  vertex_colors = V(agg)$color,
                  vertex_names = V(agg)$name)

# SET THE GRAPH COLOR ACCORDING TO EC
ec <- edge.betweenness(graph = agg,
                       e = E(agg), 
                       directed = FALSE,
                       weights = E(agg)$weight)
ec_int <- as.integer(round(ec,0))
for (k in 1:length(ec_int)){
  E(agg)$color[k] <- rev(heat.colors(1+as.integer(max(ec_int))))[as.integer(ec_int[k])+1]
  }
plot(agg,
     layout = layout.fruchterman.reingold(agg, niter=200, area=2000*vcount(agg)),
     vertex.color = "green",
     vertex.size = 7,
     vertex.label = V(agg)$name, 
     vertex.label.color = "black",
     vertex.label.font = 1.5, 
     vertex.label.cex = 1, 
     edge.width = 0.5*E(agg)$weight,
     edge.curved = TRUE,
     edge.color = E(agg)$color)

# FILTER AND REPEAT:
agg <- as.undirected(graph.adjacency(ag_m,weighted=TRUE))
V(agg)$color <- rep("green",length(ag))
V(agg)$name <- ag
cut85 <- quantile(as.vector(ag_m[ag_m>0]),0.85)
agg_f <- filter(cutoff = cut85,
                edge_matrix = ag_m,
                vertex_colors = V(agg)$color,
                vertex_names = V(agg)$name)
agg_f <- as.undirected(agg_f)
agg_f <- giant_comp(graph = agg_f,
                    vertex_color = V(agg_f)$color,
                    vertex_names = V(agg_f)$name)
ec <- edge.betweenness(graph = agg_f,
                       e = E(agg_f), 
                       directed = FALSE,
                       weights = E(agg_f)$weight)
ec_int <- as.integer(round(ec,0))
for (k in 1:length(ec_int)){
  E(agg_f)$color[k] <- rev(heat.colors(1+as.integer(max(ec_int))))[as.integer(ec_int[k])+1]
  }
plot(agg_f,
     layout = layout.fruchterman.reingold(agg_f, niter=200, area=2000*vcount(agg_f)),
     vertex.color = "green",
     vertex.size = 7,
     vertex.label = V(agg_f)$name, 
     vertex.label.color = "black",
     vertex.label.font = 1.5, 
     vertex.label.cex = 1, 
     edge.width = 0.5*E(agg_f)$weight,
     edge.curved = TRUE,
     edge.color = E(agg_f)$color)

# FILTER AND REPEAT:
agg <- as.undirected(graph.adjacency(ag_m,weighted=TRUE))
V(agg)$color <- rep("green",length(ag))
V(agg)$name <- ag
cut90 <- quantile(as.vector(ag_m[ag_m>0]),0.90)
agg_f <- filter(cutoff = cut90,
                edge_matrix = ag_m,
                vertex_colors = V(agg)$color,
                vertex_names = V(agg)$name)
agg_f <- as.undirected(agg_f)
agg_f <- giant_comp(graph = agg_f,
                    vertex_color = V(agg_f)$color,
                    vertex_names = V(agg_f)$name)
ec <- edge.betweenness(graph = agg_f,
                       e = E(agg_f), 
                       directed = FALSE,
                       weights = E(agg_f)$weight)
ec_int <- as.integer(round(ec,0))
for (k in 1:length(ec_int)){
  E(agg_f)$color[k] <- rev(heat.colors(1+as.integer(max(ec_int))))[as.integer(ec_int[k])+1]
}
plot(agg_f,
     layout = layout.fruchterman.reingold(agg_f, niter=200, area=2000*vcount(agg_f)),
     vertex.color = "green",
     vertex.size = 7,
     vertex.label = V(agg_f)$name, 
     vertex.label.color = "black",
     vertex.label.font = 1.5, 
     vertex.label.cex = 1, 
     edge.width = 0.5*E(agg_f)$weight,
     edge.curved = TRUE,
     edge.color = E(agg_f)$color)

# FILTER AND REPEAT:
agg <- as.undirected(graph.adjacency(ag_m,weighted=TRUE))
V(agg)$color <- rep("green",length(ag))
V(agg)$name <- ag
cut95 <- quantile(as.vector(ag_m[ag_m>0]),0.95)
agg_f <- filter(cutoff = cut95,
                edge_matrix = ag_m,
                vertex_colors = V(agg)$color,
                vertex_names = V(agg)$name)
agg_f <- as.undirected(agg_f)
agg_f <- giant_comp(graph = agg_f,
                    vertex_color = V(agg_f)$color,
                    vertex_names = V(agg_f)$name)
ec <- edge.betweenness(graph = agg_f,
                       e = E(agg_f), 
                       directed = FALSE,
                       weights = E(agg_f)$weight)
ec_int <- as.integer(round(ec,0))
for (k in 1:length(ec_int)){
  E(agg_f)$color[k] <- rev(heat.colors(1+as.integer(max(ec_int))))[as.integer(ec_int[k])+1]
}
plot(agg_f,
     layout = layout.fruchterman.reingold(agg_f, niter=200, area=2000*vcount(agg_f)),
     vertex.color = "green",
     vertex.size = 7,
     vertex.label = V(agg_f)$name, 
     vertex.label.color = "black",
     vertex.label.font = 1.5, 
     vertex.label.cex = 1, 
     edge.width = 0.5*E(agg_f)$weight,
     edge.curved = TRUE,
     edge.color = E(agg_f)$color)



# EDGE-BETWEENNESS COMMUNITY
# The edge betweenness score of an edge measures the number of shortest paths through it, see edge.betweenness for details. 
# The idea of the edge betweenness based community structure detection is that it is likely that edges connecting separate modules 
# have high edge betweenness as all the shortest paths from one module to another must traverse through them. 
# So if we gradually remove the edge with the highest edge betweenness score we will get a hierarchical map, a rooted tree, 
# called a dendrogram of the graph. The leafs of the tree are the individual vertices and the root of the tree represents the whole graph.
#
# edge.betweenness.community performs this algorithm by calculating the edge betweenness of the graph, 
# removing the edge with the highest edge betweenness score, then recalculating edge betweenness of the edges 
# and again removing the one with the highest score, etc.

# FILTER AND CLUSTER WITH CUTOFF = 0.75
agg <- as.undirected(graph.adjacency(ag_m,weighted=TRUE))
V(agg)$color <- rep("green",length(ag))
V(agg)$name <- ag
cut75 <- quantile(as.vector(ag_m[ag_m>0]),0.75)
agg_f <- filter(cutoff = cut75,
                edge_matrix = ag_m,
                vertex_colors = V(agg)$color,
                vertex_names = V(agg)$name)
agg_f <- as.undirected(agg_f)
agg_f <- giant_comp(graph = agg_f,
                    vertex_color = V(agg_f)$color,
                    vertex_names = V(agg_f)$name)
ebc <- edge.betweenness.community(graph = agg_f)
plot(ebc,
     agg_f, 
     vertex.size=5,
     edge.width=0.1*E(agg_f)$weight,
     main="Example: ML Communities",
     vertex.label.cex=0.8,
     vertex.label=V(agg_f)$name)

# FILTER AND CLUSTER WITH CUTOFF = 0.85
agg <- as.undirected(graph.adjacency(ag_m,weighted=TRUE))
V(agg)$color <- rep("green",length(ag))
V(agg)$name <- ag
cut85 <- quantile(as.vector(ag_m[ag_m>0]),0.85)
agg_f <- filter(cutoff = cut85,
                edge_matrix = ag_m,
                vertex_colors = V(agg)$color,
                vertex_names = V(agg)$name)
agg_f <- as.undirected(agg_f)
agg_f <- giant_comp(graph = agg_f,
                    vertex_color = V(agg_f)$color,
                    vertex_names = V(agg_f)$name)
ebc <- edge.betweenness.community(graph = agg_f)
plot(ebc,
     agg_f, 
     vertex.size=5,
     edge.width=0.02*E(agg_f)$weight,
     main="Example: ML Communities",
     vertex.label.cex=1.5,
     vertex.label=V(agg_f)$name)

# FILTER AND CLUSTER WITH CUTOFF = 0.90
agg <- as.undirected(graph.adjacency(ag_m,weighted=TRUE))
V(agg)$color <- rep("green",length(ag))
V(agg)$name <- ag
cut90 <- quantile(as.vector(ag_m[ag_m>0]),0.90)
agg_f <- filter(cutoff = cut90,
                edge_matrix = ag_m,
                vertex_colors = V(agg)$color,
                vertex_names = V(agg)$name)
agg_f <- as.undirected(agg_f)
agg_f <- giant_comp(graph = agg_f,
                    vertex_color = V(agg_f)$color,
                    vertex_names = V(agg_f)$name)
ebc <- edge.betweenness.community(graph = agg_f)
plot(ebc,
     agg_f, 
     vertex.size=3,
     edge.width=0.02*E(agg_f)$weight,
     main="Aid Agency Betweenness Community Structure",
     vertex.label.cex=1.2,
     vertex.label=V(agg_f)$name)

# FILTER AND CLUSTER WITH CUTOFF = 0.95
agg <- as.undirected(graph.adjacency(ag_m,weighted=TRUE))
V(agg)$color <- rep("green",length(ag))
V(agg)$name <- ag
cut95 <- quantile(as.vector(ag_m[ag_m>0]),0.95)
agg_f <- filter(cutoff = cut95,
                edge_matrix = ag_m,
                vertex_colors = V(agg)$color,
                vertex_names = V(agg)$name)
agg_f <- as.undirected(agg_f)
agg_f <- giant_comp(graph = agg_f,
                    vertex_color = V(agg_f)$color,
                    vertex_names = V(agg_f)$name)
ebc <- edge.betweenness.community(graph = agg_f)
plot(ebc,
     agg_f, 
     vertex.size=5,
     edge.width=0.5*E(agg_f)$weight,
     main="Example: ML Communities",
     vertex.label.cex=0.8,
     vertex.label=V(agg_f)$name)



# BETWEENNESS ESTIMATE: This measure calculates betweenness by considering only paths of a certain length 
# that is smaller than or equal to the cutoff value. Similarly for edge betweenness estimates. 
# In our analysis, we will use cutoff lengths 10-30 (refer to the path length distribution analysis above). 
# Here, we use a cutoff=3 just to illustrate the application.

be<-betweenness.estimate(vgg,v=V(vgg), directed=FALSE,3)
plot(sort(be, decreasing=TRUE), xlab="Node Index", ylab="Betweenness Centrality Estimate for g", main="Betweenness Centrality Estimate for g",col=adjustcolor(rgb(0,1/2,0,1)))
hist(be,breaks=200,col=adjustcolor(rgb(0,1/2,0,1)),xlab="Betweenness Centrality Estimate Values (c=3) for g",main="Positive Betweenness Centrality Estimate for g")





# BETWEENNESS ESTIMATE CENTRALITY: THE NUMBER OF GEODESICS GOING THROUGH A NODE
agg <- as.undirected(graph.adjacency(ag_m,weighted=TRUE))
V(agg)$color <- rep("green",length(ag))
V(agg)$name <- ag
bce <- betweenness.estimate(graph = agg,
                            vids = V(agg),
                            weights = E(agg)$weight,
                            directed=FALSE,
                            4)
plot(sort(bce, decreasing=TRUE),
     col=adjustcolor(rgb(1/2,0,0,1/2)), 
     xlab="Node Index", 
     ylab="Betweenness Centrality Estimate", 
     main="Sorted Relief Agency Betweenness Centrality Values", 
     pch=19)
histP2(bce,
       breaks=100,
       col=adjustcolor(rgb(1/2,0,0,1/2)),
       xlab="Betweenness Estimate Centrality Values",
       main="Relief Agency Betweenness Estimate Centrality Distribution")

# PLOT HEAT MAP ON VERTICES ACCORDING TO BETWEENNESS CENTRALITY
bce_int <- as.integer(round(bce,0))
for (k in 1:length(bce_int)){
  V(agg)$color[k] <- rev(heat.colors(1+as.integer(max(bce_int))))[as.integer(bce_int[k])+1]
}

# PLOT AGENCY GRAPH AND FILTER
plot(agg,
     layout = layout.fruchterman.reingold(agg, niter=200, area=2000*vcount(agg)),
     vertex.color = V(agg)$color,
     vertex.size = 7,
     vertex.label = V(agg)$name, 
     vertex.label.color = "black",
     vertex.label.font = 1.5, 
     vertex.label.cex = 1, 
     edge.width = 0.5*E(agg)$weight,
     edge.curved = TRUE,
     edge.color = gray.colors(1))

# FIND THE TOP 10% BETWEENNES NODES
top_bce <- bce[which(bce > quantile(bce,0.9))]
top_bce

# FIND THE TOP 5% BETWENNES NODES
top_bce <- bce[which(bce > quantile(bce,0.95))]
top_bce

# REMOVE ISOLATED
agg <- drop_isolated(graph = agg,
                     vertex_colors = V(agg)$color,
                     vertex_names = V(agg)$name)

# GET THE GIANT CONNECTED COMPONENT (TWO CLSUTERS ONLY)
agg <- giant_comp(graph = agg,
                  vertex_colors = V(agg)$color,
                  vertex_names = V(agg)$name)

# SET THE GRAPH COLOR ACCORDING TO BC
bce <- betweenness.estimate(graph = agg,
                            vids = V(agg),
                            weights = E(agg)$weight,
                            directed=FALSE,
                            cutoff = 4)
bce_int <- as.integer(round(bce,0))
for (k in 1:length(bce_int)){
  V(agg)$color[k] <- rev(heat.colors(1+as.integer(max(bce_int))))[as.integer(bce_int[k])+1]
}
plot(agg,
     layout = layout.fruchterman.reingold(agg, niter=200, area=2000*vcount(agg)),
     vertex.color = V(agg)$color,
     vertex.size = 5,
     vertex.label = NA, 
     vertex.label.color = "black",
     vertex.label.font = 1, 
     vertex.label.cex = 0.5, 
     edge.width = 0.1*E(agg)$weight,
     edge.curved = TRUE,
     edge.color = gray.colors(1))


# EDGE BETWEENNESS ESTIMATE: This betweenness estimate is defined in a similar way as betweenness estimate, 
# we show it here with cutoff=4.

agg <- as.undirected(graph.adjacency(ag_m,weighted=TRUE))
V(agg)$color <- rep("green",length(ag))
V(agg)$name <- ag
ebe <- edge.betweenness.estimate(graph = agg,
                                weights =E(agg)$weight,
                                directed=FALSE,
                                cutoff = 4)
plot(sort(ebe, decreasing=TRUE),
     col=adjustcolor(rgb(1/2,0,0,1/2)), 
     xlab="Node Index", 
     ylab="Edge-Betweenness Centrality", 
     main="Sorted Relief Agency Edge-Betweenness Centrality Values", 
     pch=19)
histP2(ebe,
       breaks=100,
       col=adjustcolor(rgb(1/2,0,0,1/2)),
       xlab="Edge-Betweenness Centrality Values",
       main="Relief Agency Edge-Betweenness Centrality Distribution")

# PLOT HEAT MAP ON VERTICES ACCORDING TO EDGE-BEWEENNESS CENTRALITY
ebe_int <- as.integer(round(ebe,0))
for (k in 1:length(ebe_int)){
  E(agg)$color[k] <- rev(heat.colors(1+as.integer(max(ebe_int))))[as.integer(ebe_int[k])+1]
}

# PLOT AGENCY GRAPH AND FILTER
plot(agg,
     layout = layout.fruchterman.reingold(agg, niter=200, area=2000*vcount(agg)),
     vertex.color = "green",
     vertex.size = 7,
     vertex.label = V(agg)$name, 
     vertex.label.color = "black",
     vertex.label.font = 1.5, 
     vertex.label.cex = 1, 
     edge.width = 0.5*E(agg)$weight,
     edge.curved = TRUE,
     edge.color = E(agg)$color)

# FIND THE TOP 10% EDGE-BETWEENNES EDGES
top_ebe <- ebe[which(ebe > quantile(ebe,0.9))]
E(agg)[which(ebe %in% top_ebe)]

# FIND THE TOP 5% EDGE-BETWEENNES EDGES
top_ebe <- ebe[which(ebe > quantile(ebe,0.99))]
E(agg)[which(ebe %in% top_ebe)]

# REMOVE ISOLATED
agg <- drop_isolated(graph = agg,
                     vertex_colors = V(agg)$color,
                     vertex_names = V(agg)$name)

# GET THE GIANT CONNECTED COMPONENT (TWO CLSUTERS ONLY)
agg <- giant_comp(graph = agg,
                  vertex_colors = V(agg)$color,
                  vertex_names = V(agg)$name)

# SET THE GRAPH COLOR ACCORDING TO EC
ebe <- edge.betweenness.estimate(graph = agg,
                                 weights =E(agg)$weight,
                                 directed=FALSE,
                                 cutoff = 4)
ebe_int <- as.integer(round(ebe,0))
for (k in 1:length(ec_int)){
  E(agg)$color[k] <- rev(heat.colors(1+as.integer(max(ebe_int))))[as.integer(ebe_int[k])+1]
}
plot(agg,
     layout = layout.fruchterman.reingold(agg, niter=200, area=2000*vcount(agg)),
     vertex.color = "green",
     vertex.size = 7,
     vertex.label = V(agg)$name, 
     vertex.label.color = "black",
     vertex.label.font = 1.5, 
     vertex.label.cex = 1, 
     edge.width = 0.5*E(agg)$weight,
     edge.curved = TRUE,
     edge.color = E(agg)$color)


# WE CAN FURTHER FILTER AND APPLY THIS DEPENDING ON MODELING GOALS




# CLOSENESS CENTRALITY: This measure takes into account the distribution of distances to other nodes from a given node. 
# It is defined as the reciprocal of the farness of a node, where farness is defined as the sum of its distances to all 
# other nodes. Closeness can be regarded as a measure of how long it will take to spread information 
# (or an efect of an event) from a node to all other nodes. To demonstrate this concept, we compute the closeness 
# centrality for the unweighted network.
agg <- as.undirected(graph.adjacency(ag_m,weighted=TRUE))
V(agg)$color <- rep("green",length(ag))
V(agg)$name <- ag
cl<-clusters(agg)
agg1<-induced.subgraph(agg, which(cl$membership == which.max(cl$csize)))
cc<-closeness(graph = agg1,vids = V(agg1),weights = E(agg1)$weight)
plot(sort(cc/max(cc), decreasing=TRUE), col=adjustcolor(rgb(1/2,0,1,1)), xlab="Node Id in the Giant Conencted Component (gg1)", ylab="Normalized Closeness Centrality", main="Closeness Centrality for the Giant Component (gg1)")
hist(cc/max(cc),breaks=200,col=adjustcolor(rgb(1/2,0,1,1)),xlab="Normalized Closeness Centrality Values for gg1",main="Normalized Closeness Centrality Distribution for gg1")

cce<-closeness.estimate(graph = agg1,vids = V(agg1),weights = E(agg1)$weight,cutoff = 7)
plot(sort(cce/max(cce), decreasing=TRUE), col=adjustcolor(rgb(1/2,0,1,1)), xlab="Node Id in the Giant Conencted Component (gg1)", ylab="Normalized Closeness Centrality", main="Closeness Centrality for the Giant Component (gg1)")
hist(cce/max(cce),breaks=200,col=adjustcolor(rgb(1/2,0,1,1)),xlab="Normalized Closeness Centrality Values for gg1",main="Normalized Closeness Centrality Distribution for gg1")





# EIGENVECTOR CENTRALITY: This is a measure of the influence of a node in the network. 
# It assigns relative scores to all nodes in the network based on the concept that connections to high-scoring nodes 
# contribute more to the score of the given node than equal conenctions to low-scoring nodes. 
# A variant of egenvector centrality is Google's PageRank algorithm.
agg <- as.undirected(graph.adjacency(ag_m,weighted=TRUE))
V(agg)$color <- rep("green",length(ag))
V(agg)$name <- ag
clu<-clusters(agg)
agg1<-induced.subgraph(agg, which(clu$membership == which.max(clu$csize)))
ec<-evcent(agg1)$vector
plot(sort(ec, decreasing=TRUE), col=adjustcolor(rgb(0,0,1,1/2)), xlab="Node Id in the Network (1:200)", ylab="Closeness Centrality Values", main="Essential (first 200 nodes) Closeness Centrality for g", pch=20)
hist(ec,breaks=100,col=adjustcolor(rgb(0,0,1,1/2)),xlab="Closeness Centrality Values",main="Essential Closeness Centrality Distribution")


# AUTHORITY SCORE: This is a measure for DIRECTED NETWORKS, and it measures the number of nodes that are hubs and point 
# to a given node. It is defined as the principle eigenvector values for t(A)*A, where A stands for the adjacency 
# matrix of the network. For undirected networks like ours, the adjacency matrix is symmetric, so the authority score 
# is equivalent to the hub score. In subsequent analyses, we will be looking at directed extensions of this network 
# model, so we are including these two scores in the analysis for completeness.
agg <- as.undirected(graph.adjacency(ag_m,weighted=TRUE))
V(agg)$color <- rep("green",length(ag))
V(agg)$name <- ag
au<-authority.score(agg)$vector
plot(sort(au, decreasing=TRUE), col=adjustcolor(rgb(0,0,1,1/2)), xlab="Node Id in the Network (1:200)", ylab="Authority Score Values", main="Essential (first 200 nodes) Authority Scores for g", pch=20)
hist(au,breaks=100,col=adjustcolor(rgb(0,0,1,1/2)),xlab="Authority Score Values",main="Essential Authority Score Distribution")



# HUB SCORE: This is a measure FOR DIRECTED NETWORKS and it measures the number of authority nodes that a given hub node points to.
agg <- as.undirected(graph.adjacency(ag_m,weighted=TRUE))
V(agg)$color <- rep("green",length(ag))
V(agg)$name <- ag
hb<-hub.score(agg)$vector
plot(sort(hb, decreasing=TRUE), col=adjustcolor(rgb(0,0,1,1/2)), xlab="Node Id in the Network (1:200)", ylab="Hub Score Values", main="Essential (first 200 nodes) Hub Scores for g", pch=20)
hist(hb,breaks=100,col=adjustcolor(rgb(0,0,1,1/2)),xlab="Hub Score Values",main="Essential Hub Score Distribution")


# CLUSTERING COEFFICIENTS: This is a measure of the clustering of the network,defined by the ratio of the number of closed triplets 
# and the number of connected triplets of vertices. We computed it in the previous report, but here we include the local and weighted version 
# of the clustering coefficients. Clustering is particularly relevant to social netowrks where nodes tend to create tightly knit groups charaterized 
# by a high density of ties, this likelihood is greater than the average probability of an edge between two randomly selected nodes.
agg <- as.undirected(graph.adjacency(ag_m,weighted=TRUE))
V(agg)$color <- rep("green",length(ag))
V(agg)$name <- ag
tr<-transitivity(agg, type="local")
plot(sort(tr), col=adjustcolor(rgb(0,0,1,1/2)), xlab="Node Id in the Network (1:3200)", ylab="Clustering Coefficient Values", main="Essential Clustering Coefficients for g", pch=20)
hist(tr,breaks=100,col=adjustcolor(rgb(0,0,1,1/2)),xlab="Clustering Coefficient Values",main="Clustering Coefficient Distribution for g")


# The weighted analogue of the clustering coefficient
trw<-transitivity(agg, type="weighted")
plot(sort(trw), col=adjustcolor(rgb(0,0,1,1/2)), xlab="Node Id in the Network (1:200)", ylab="Clustering Coefficient Values", main="Essential Clustering Coefficients for g", pch=20)
hist(trw,breaks=100,col=adjustcolor(rgb(0,0,1,1/2)),xlab="Clustering Coefficient Values",main="Clustering Coefficient Distribution for g")


# COMMUNITY STRUCTURES: This is a way of performing funcitonal clustering in complex networks. We have already looked at the connected components, 
# this is an elementary community detection based on connectivity.
strongclusters<-clusters(agg)$membership
plot(agg,vertex.color=strongclusters, layout=layout.fruchterman.reingold,vertex.size=4, edge.color="black", edge.width=E(agg)$weight,vertex.label=NA,main="Clustering for Store Network g200")

# ADD SOME FILTERING AND TRY AGAIN

mc<-multilevel.community(agg)
plot(sort(mc$membership), col=adjustcolor(rgb(0,0,1,1/2)), xlab="Node Id in the Network", ylab="Multilevel Community Values", main="Multilevel Community Values for g", pch=20)
hist(mc$membership,breaks=100,col=adjustcolor(rgb(0,0,1,1/2)),xlab="Multilevel Community Values",main="Multilevel Community Distribution")


# Next, we show the walktrap community algorithm.
wc<-walktrap.community(agg)
plot(sort(wc$membership), col=adjustcolor(rgb(0,0,1,1/2)), xlab="Node Id in the Network", ylab="Walktrap Community Values", main="Walktrap Community Values for g", pch=20)
hist(wc$membership,breaks=100,col=adjustcolor(rgb(0,0,1,1/2)),xlab="Walktrap Community Values",main="Walktrap Community Distribution")


plot(wc,agg,vertex.size=4, vertex.label=NA,edge.width=E(agg)$weight,main="Walktrap Community Detection for g200")
plot(agg, vertex.color=membership(wc), vertex.size=6, edge.color="black", edge.width=E(agg)$weight,vertex.label=NA,main="Walktrap Community Detection for g200")
















# REMINDER: RUN IT FOR THE TOTAL AGENCY-VDC NETWORK
au<-authority.score(av)$vector
plot(sort(au, decreasing=TRUE), col=adjustcolor(rgb(0,0,1,1/2)), xlab="Node Id in the Network (1:200)", ylab="Authority Score Values", main="Essential (first 200 nodes) Authority Scores for g", pch=20)
hist(au,breaks=100,col=adjustcolor(rgb(0,0,1,1/2)),xlab="Authority Score Values",main="Essential Authority Score Distribution")

# HUB SCORE: This is a measure FOR DIRECTED NETWORKS and it measures the number of authority nodes that a given hub node points to.
agg <- as.undirected(graph.adjacency(ag_m,weighted=TRUE))
V(agg)$color <- rep("green",length(ag))
V(agg)$name <- ag
hb<-hub.score(agg)$vector
plot(sort(hb, decreasing=TRUE), col=adjustcolor(rgb(0,0,1,1/2)), xlab="Node Id in the Network (1:200)", ylab="Hub Score Values", main="Essential (first 200 nodes) Hub Scores for g", pch=20)
hist(hb,breaks=100,col=adjustcolor(rgb(0,0,1,1/2)),xlab="Hub Score Values",main="Essential Hub Score Distribution")

# RUN IT FOR THE TOTAL AGENCY-VDC NETWORK
hb<-hub.score(av)$vector
plot(sort(hb, decreasing=TRUE), col=adjustcolor(rgb(0,0,1,1/2)), xlab="Node Id in the Network (1:200)", ylab="Hub Score Values", main="Essential (first 200 nodes) Hub Scores for g", pch=20)
hist(hb,breaks=100,col=adjustcolor(rgb(0,0,1,1/2)),xlab="Hub Score Values",main="Essential Hub Score Distribution")
































# VDC AID TARGET NETWORK:


# ANALYSIS OF THE VDC AID TARGET NETWORK
u_vdc <- as.character(unique(aid_data$vdc))

# BUILD THE SHARED AGENCY ASSOCIATION NETWORK FOR THE VDCs
aid_vdc <- matrix(0, nrow = length(u_vdc), ncol = length(u_vdc))
for (i in 1:length(u_vdc)){
  for (j in 1:length(u_vdc)){
    common <- aid_m[1:length(ag),c(i,j)]
    common[common>0] <-1
    aid_vdc[i,j]<-sum((common[,1])*(common[,2]))
  }
}


# REMOVE SELF LOOPS
for (k in 1:dim(aid_vdc)[1]){aid_vdc[[k,k]] <- 0}

# DEFINE AGENCY GRAPH
vgg <- as.undirected(graph.adjacency(aid_vdc,weighted=TRUE))

# SET THE GRAPH COLOR
V(vgg)$color <- rep("SkyBlue2",length(u_vdc))
V(vgg)$name <- u_vdc

# PLOT AGENCY GRAPH AND FILTER
plot(vgg,
     layout = layout.fruchterman.reingold(vgg, niter=200, area=2000*vcount(vgg)),
     vertex.color = "SkyBlue2",
     vertex.size = 2,
     vertex.label = NA, 
     vertex.label.color = "black",
     vertex.label.font = 1, 
     vertex.label.cex = 0.5, 
     edge.width = 0.5*E(vgg)$weight,
     edge.curved = TRUE,
     edge.color = gray.colors(1))

# REMOVE ISOLATED
vgg <- drop_isolated(graph = vgg,
                     vertex_colors = V(vgg)$color,
                     vertex_names = V(vgg)$name)

# GET THE GIANT CONENCTED COMPONENT (TWO CLSUTERS ONLY)
vgg <- giant_comp(graph = vgg,
                  vertex_colors = V(vgg)$color,
                  vertex_names = V(vgg)$name)

plot(vgg,
     layout = layout.fruchterman.reingold(vgg, niter=200, area=2000*vcount(vgg)),
     vertex.color = "SkyBlue2",
     vertex.size = 2,
     vertex.label = NA, 
     vertex.label.color = "black",
     vertex.label.font = 1, 
     vertex.label.cex = 0.5, 
     edge.width = 0.1*E(vgg)$weight,
     edge.curved = TRUE,
     edge.color = gray.colors(1),
     main="VDC Aid Association Network")


# FILTER BY EDGE WEIGHT LEVELS, NOTE THAT LOWER CUTOFF VALUES
# DO NOT RESULT IN SIGNIFICANT FILTRATION

# THIS INITIAL FILTRATION IS REDUNDANT SINCE MINIMAL DEGREE IS 1
vgg <- as.undirected(graph.adjacency(aid_vdc,weighted=TRUE))
V(vgg)$color <- rep("SkyBlue2",length(u_vdc))
V(vgg)$name <- u_vdc
transitivity(vgg)
cut85 <- quantile(as.vector(aid_vdc[aid_vdc>0]),0.85)
vgg_f <- filter(cutoff = cut85,
                edge_matrix = aid_vdc,
                vertex_colors = V(vgg)$color,
                vertex_names = V(vgg)$name)
vgg_f <- as.undirected(vgg_f)
vgg_f <- giant_comp(graph = vgg_f,
                    vertex_colors = V(vgg_f)$color,
                    vertex_names = V(vgg_f)$name)
plot(as.undirected(vgg_f),
     layout=layout.fruchterman.reingold(vgg_f, niter=200, area=2000*vcount(vgg_f)),
     vertex.color = V(vgg_f)$color,
     vertex.size=2,
     vertex.label=NA, 
     vertex.label.color="black", 
     vertex.label.font=1, 
     vertex.label.cex=1, 
     edge.width=0.25*(E(vgg_f)$weight),
     edge.curved=TRUE,
     edge.color=gray.colors(1))

# THE NEXT CUT IS TOO BIG OF A JUMP, WHICH WILL BE EVIDENT IN THE DEGREE DISTRIBUTION
vgg <- as.undirected(graph.adjacency(aid_vdc,weighted=TRUE))
vgg <- as.undirected(graph.adjacency(aid_vdc,weighted=TRUE))
V(vgg)$color <- rep("SkyBlue2",length(u_vdc))
V(vgg)$name <- u_vdc
transitivity(vgg)
cut95 <- quantile(as.vector(aid_vdc[aid_vdc>0]),0.95)
vgg_f <- filter(cutoff = cut95,
                edge_matrix = aid_vdc,
                vertex_colors = V(vgg)$color,
                vertex_names = V(vgg)$name)
vgg_f <- as.undirected(vgg_f)
vgg_f <- giant_comp(graph = vgg_f,
                    vertex_colors = V(vgg_f)$color,
                    vertex_names = V(vgg_f)$name)
plot(as.undirected(vgg_f),
     layout=layout.fruchterman.reingold(vgg_f, niter=200, area=2000*vcount(vgg_f)),
     vertex.color = V(vgg_f)$color,
     vertex.size=2,
     vertex.label=NA, 
     vertex.label.color="black", 
     vertex.label.font=1, 
     vertex.label.cex=1, 
     edge.width=0.25*(E(vgg_f)$weight),
     edge.curved=TRUE,
     edge.color=gray.colors(1))

# FURTHER FILTRATION
vgg <- as.undirected(graph.adjacency(aid_vdc,weighted=TRUE))
V(vgg)$color <- rep("SkyBlue2",length(u_vdc))
V(vgg)$name <- u_vdc
transitivity(vgg)
cut99.5 <- quantile(as.vector(aid_vdc[aid_vdc>0]),0.995)
vgg_f <- filter(cutoff = cut99.5,
                edge_matrix = aid_vdc,
                vertex_colors = V(vgg)$color,
                vertex_names = V(vgg)$name)
vgg_f <- as.undirected(vgg_f)
vgg_f <- giant_comp(graph = vgg_f,
                    vertex_colors = V(vgg_f)$color,
                    vertex_names = V(vgg_f)$name)
plot(as.undirected(vgg_f),
     layout=layout.fruchterman.reingold(vgg_f, niter=200, area=2000*vcount(vgg_f)),
     vertex.color = V(vgg_f)$color,
     vertex.size=4,
     vertex.label=NA, 
     vertex.label.color="black", 
     vertex.label.font=1, 
     vertex.label.cex=1, 
     edge.width=sqrt(E(vgg_f)$weight),
     edge.curved=TRUE,
     edge.color=gray.colors(1))


# ANALYSIS OF THE TARGET VDC NETWORK ITSELF:

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
vgg <- as.undirected(graph.adjacency(aid_vdc,weighted=TRUE))
V(vgg)$color <- rep("SkyBlue2",length(u_vdc))
V(vgg)$name <- u_vdc
transitivity(vgg)
cut75 <- quantile(as.vector(aid_vdc[aid_vdc>0]),0.75)
vgg_f <- filter(cutoff = cut75,
                edge_matrix = aid_vdc,
                vertex_colors = V(vgg)$color,
                vertex_names = V(vgg)$name)
vgg_f <- as.undirected(vgg_f)
transitivity(vgg_f)

# RELATIVE MAXIMAL CLUSTER SIZE (AS % OF NUMBER OF NODES) 
max(clusters(vgg)$csize)/vcount(vgg)

# RELATIVE NUMBER OF ISOLATED NODES (AS % OF NUMBER OF NODES)  
sum(degree(vgg)==0)/vcount(vgg)

# PATH DISTRIBUTION: This shows the different lengths of shortest paths (geodesics) in our network. 
vgg <- as.undirected(graph.adjacency(aid_vdc,weighted=TRUE))
V(vgg)$color <- rep("SkyBlue2",length(u_vdc))
V(vgg)$name <- u_vdc
vgg <- giant_comp(graph = vgg,
                  vertex_color = V(vgg)$color,
                  vertex_names = V(vgg)$name)
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
     col=adjustcolor(rgb(1,0,1/2,1)),
     xlab="Path Length Values",
     main="Path Length Distribution for g")



# BETWEENNESS CENTRALITY: THE NUMBER OF GEODESICS GOING THROUGH A NODE
bc <- betweenness(vgg,v=V(vgg), directed=FALSE)
plot(sort(bc, decreasing=TRUE),
     col=adjustcolor(rgb(1/2,0,0,1/2)), 
     xlab="Node Index", 
     ylab="Betweenness Centrality", 
     main="Sorted Relief Agency Betweenness Centrality Values", 
     pch=19)
histP2(bc,
       breaks=100,
       col=adjustcolor(rgb(1/2,0,0,1/2)),
       xlab="Betweenness Centrality Values",
       main="Relief Agency Betweenness Centrality Distribution")

# PLOT HEAT MAP ON VERTICES ACCORDING TO BETWEENNESS CENTRALITY
bc_int <- as.integer(round(bc,0))
for (k in 1:length(bc_int)){
  V(vgg)$color[k] <- rev(heat.colors(1+as.integer(max(bc_int))))[as.integer(bc_int[k])+1]
}

# PLOT AGENCY GRAPH AND FILTER
plot(vgg,
     layout = layout.fruchterman.reingold(vgg, niter=200, area=2000*vcount(vgg)),
     vertex.color = V(vgg)$color,
     vertex.size = 4,
     vertex.label = NA, 
     vertex.label.color = "black",
     vertex.label.font = 1.5, 
     vertex.label.cex = 1, 
     edge.width = 0.5*E(vgg)$weight,
     edge.curved = TRUE,
     edge.color = gray.colors(1))

# FIND THE TOP 10% BETWEENNES NODES
top_bc <- bc[which(bc > quantile(bc,0.9))]
top_bc

# FIND THE TOP 5% BETWENNES NODES
top_bc <- bc[which(bc > quantile(bc,0.95))]
top_bc

# REMOVE ISOLATED
vgg <- drop_isolated(graph = vgg,
                     vertex_colors = V(vgg)$color,
                     vertex_names = V(vgg)$name)

# GET THE GIANT CONNECTED COMPONENT (TWO CLSUTERS ONLY)
vgg <- giant_comp(graph = vgg,
                  vertex_colors = V(vgg)$color,
                  vertex_names = V(vgg)$name)

# SET THE GRAPH COLOR ACCORDING TO BC

bc<-betweenness(vgg,v=V(vgg), directed=FALSE)
bc_int <- as.integer(round(bc,0))
for (k in 1:length(bc_int)){
  V(vgg)$color[k] <- rev(heat.colors(1+as.integer(max(bc_int))))[as.integer(bc_int[k])+1]
}
plot(vgg,
     layout = layout.fruchterman.reingold(vgg, niter=2000, area=20000*vcount(vgg)),
     vertex.color = V(vgg)$color,
     vertex.size = 3,
     vertex.label = NA, 
     vertex.label.color = "black",
     vertex.label.font = 1, 
     vertex.label.cex = 0.5, 
     edge.width = 0.1*E(vgg)$weight,
     edge.curved = TRUE,
     edge.color = gray.colors(1),
     main="VDC Aid Network Betweenness Centrality Heat Map")






















# FILTER AND REPEAT:
vgg <- as.undirected(graph.adjacency(aid_vdc,weighted=TRUE))
V(vgg)$color <- rep("green",length(u_vdc))
V(vgg)$name <- u_vdc
cut85 <- quantile(as.vector(aid_vdc[aid_vdc>0]),0.85)
vgg_f <- filter(cutoff = cut85,
                edge_matrix = aid_vdc,
                vertex_colors = V(vgg)$color,
                vertex_names = V(vgg)$name)
vgg_f <- as.undirected(vgg_f)
vgg_f <- giant_comp(graph = vgg_f,
                    vertex_color = V(vgg_f)$color,
                    vertex_names = V(vgg_f)$name)
bc<-betweenness(vgg_f,v=V(vgg_f), directed=FALSE)
bc_int <- as.integer(round(bc,0))
for (k in 1:length(bc_int)){
  V(vgg_f)$color[k] <- rev(heat.colors(1+as.integer(max(bc_int))))[as.integer(bc_int[k])+1]
}
plot(vgg_f,
     layout = layout.fruchterman.reingold(vgg_f, niter=200, area=2000*vcount(vgg_f)),
     vertex.color = V(vgg_f)$color,
     vertex.size = 5,
     vertex.label = NA, 
     vertex.label.color = "black",
     vertex.label.font = 1, 
     vertex.label.cex = 0.5, 
     edge.width = 0.1*E(vgg_f)$weight,
     edge.curved = TRUE,
     edge.color = gray.colors(1))

# FILTER AND REPEAT:
vgg <- as.undirected(graph.adjacency(aid_vdc,weighted=TRUE))
V(vgg)$color <- rep("green",length(u_vdc))
V(vgg)$name <- u_vdc
cut95 <- quantile(as.vector(aid_vdc[aid_vdc>0]),0.95)
vgg_f <- filter(cutoff = cut95,
                edge_matrix = aid_vdc,
                vertex_colors = V(vgg)$color,
                vertex_names = V(vgg)$name)
vgg_f <- as.undirected(vgg_f)
vgg_f <- giant_comp(graph = vgg_f,
                    vertex_color = V(vgg_f)$color,
                    vertex_names = V(vgg_f)$name)
bc<-betweenness(graph = vgg_f,
                v=V(vgg_f), 
                directed=FALSE,
                weights = E(vgg_f)$weight)
bc_int <- as.integer(round(bc,0))
for (k in 1:length(bc_int)){
  V(vgg_f)$color[k] <- rev(heat.colors(1+as.integer(max(bc_int))))[as.integer(bc_int[k])+1]
}
plot(vgg_f,
     layout = layout.fruchterman.reingold(vgg_f, niter=200, area=2000*vcount(vgg_f)),
     vertex.color = V(vgg_f)$color,
     vertex.size = 5,
     vertex.label = NA, 
     vertex.label.color = "black",
     vertex.label.font = 1, 
     vertex.label.cex = 0.5, 
     edge.width = 0.1*E(vgg_f)$weight,
     edge.curved = TRUE,
     edge.color = gray.colors(1))

# FILTER AND REPEAT:
vgg <- as.undirected(graph.adjacency(aid_vdc,weighted=TRUE))
V(vgg)$color <- rep("green",length(u_vdc))
V(vgg)$name <- u_vdc
cut97 <- quantile(as.vector(aid_vdc[aid_vdc>0]),0.97)
vgg_f <- filter(cutoff = cut97,
                edge_matrix = aid_vdc,
                vertex_colors = V(vgg)$color,
                vertex_names = V(vgg)$name)
vgg_f <- as.undirected(vgg_f)
vgg_f <- giant_comp(graph = vgg_f,
                    vertex_color = V(vgg_f)$color,
                    vertex_names = V(vgg_f)$name)
bc<-betweenness(graph = vgg_f,
                v=V(vgg_f), 
                directed=FALSE,
                weights = E(vgg_f)$weight)
bc_int <- as.integer(round(bc,0))
for (k in 1:length(bc_int)){
  V(vgg_f)$color[k] <- rev(heat.colors(1+as.integer(max(bc_int))))[as.integer(bc_int[k])+1]
}
plot(vgg_f,
     layout = layout.fruchterman.reingold(vgg_f, niter=200, area=2000*vcount(vgg_f)),
     vertex.color = V(vgg_f)$color,
     vertex.size = 5,
     vertex.label = NA, 
     vertex.label.color = "black",
     vertex.label.font = 1, 
     vertex.label.cex = 0.5, 
     edge.width = 0.1*E(vgg_f)$weight,
     edge.curved = TRUE,
     edge.color = gray.colors(1))



# FILTER AND REPEAT:
vgg <- as.undirected(graph.adjacency(aid_vdc,weighted=TRUE))
V(vgg)$color <- rep("green",length(u_vdc))
V(vgg)$name <- u_vdc
cut99.8 <- quantile(as.vector(aid_vdc[aid_vdc>0]),0.998)
vgg_f <- filter(cutoff = cut99.8,
                edge_matrix = aid_vdc,
                vertex_colors = V(vgg)$color,
                vertex_names = V(vgg)$name)
vgg_f <- as.undirected(vgg_f)
vgg_f <- giant_comp(graph = vgg_f,
                    vertex_color = V(vgg_f)$color,
                    vertex_names = V(vgg_f)$name)
bc<-betweenness(graph = vgg_f,
                v=V(vgg_f), 
                directed=FALSE,
                weights = E(vgg_f)$weight)
bc_int <- as.integer(round(bc,0))
for (k in 1:length(bc_int)){
  V(vgg_f)$color[k] <- rev(heat.colors(1+as.integer(max(bc_int))))[as.integer(bc_int[k])+1]
}
plot(vgg_f,
     layout = layout.fruchterman.reingold(vgg_f, niter=200, area=2000*vcount(vgg_f)),
     vertex.color = V(vgg_f)$color,
     vertex.size = 12,
     vertex.label = V(vgg_f)$name, 
     vertex.label.color = "black",
     vertex.label.font = 1, 
     vertex.label.cex = 1.3, 
     edge.width = 0.7*E(vgg_f)$weight,
     edge.curved = TRUE,
     edge.color = gray.colors(1))





# EDGE-BETWEENNESS CENTRALITY: THE NUMBER OF GEODESICS GOING THROUGH AN EDGE
vgg <- as.undirected(graph.adjacency(aid_vdc,weighted=TRUE))
V(vgg)$color <- rep("green",length(u_vdc))
V(vgg)$name <- u_vdc
ec <- edge.betweenness(graph = vgg,
                       e = E(vgg), 
                       directed = FALSE,
                       weights = E(vgg)$weight)
plot(sort(ec, decreasing=TRUE),
     col=adjustcolor(rgb(1/2,0,0,1/2)), 
     xlab="Node Index", 
     ylab="Edge-Betweenness Centrality", 
     main="Sorted Relief Agency Edge-Betweenness Centrality Values", 
     pch=19)
histP2(ec,
       breaks=100,
       col=adjustcolor(rgb(1/2,0,0,1/2)),
       xlab="Edge-Betweenness Centrality Values",
       main="Relief Agency Edge-Betweenness Centrality Distribution")

# PLOT HEAT MAP ON VERTICES ACCORDING TO EDGE-BEWEENNESS CENTRALITY
ec_int <- as.integer(round(ec,0))
for (k in 1:length(ec_int)){
  E(vgg)$color[k] <- rev(heat.colors(1+as.integer(max(ec_int))))[as.integer(ec_int[k])+1]
}

# PLOT AGENCY GRAPH AND FILTER
plot(vgg,
     layout = layout.fruchterman.reingold(vgg, niter=200, area=2000*vcount(vgg)),
     vertex.color = "green",
     vertex.size = 2,
     vertex.label = NA, 
     vertex.label.color = "black",
     vertex.label.font = 1.5, 
     vertex.label.cex = 1, 
     edge.width = 0.5*E(vgg)$weight,
     edge.curved = TRUE,
     edge.color = E(vgg)$color)

# FIND THE TOP 10% EDGE-BETWEENNES EDGES
top_ec <- ec[which(ec > quantile(ec,0.999))]
E(vgg)[which(ec %in% top_ec)]

# FIND THE TOP 5% EDGE-BETWEENNES EDGES
top_ec <- ec[which(ec > quantile(ec,0.9995))]
E(vgg)[which(ec %in% top_ec)]

# REMOVE ISOLATED
vgg <- drop_isolated(graph = vgg,
                     vertex_colors = V(vgg)$color,
                     vertex_names = V(vgg)$name)

# GET THE GIANT CONNECTED COMPONENT (TWO CLSUTERS ONLY)
vgg <- giant_comp(graph = vgg,
                  vertex_colors = V(vgg)$color,
                  vertex_names = V(vgg)$name)

# SET THE GRAPH COLOR ACCORDING TO EC
ec <- edge.betweenness(graph = vgg,
                       e = E(vgg), 
                       directed = FALSE,
                       weights = E(vgg)$weight)
ec_int <- as.integer(round(ec,0))
for (k in 1:length(ec_int)){
  E(vgg)$color[k] <- rev(heat.colors(1+as.integer(max(ec_int))))[as.integer(ec_int[k])+1]
}
plot(vgg,
     layout = layout.fruchterman.reingold(vgg, niter=200, area=2000*vcount(vgg)),
     vertex.color = "green",
     vertex.size = 2,
     vertex.label = NA, 
     vertex.label.color = "black",
     vertex.label.font = 1.5, 
     vertex.label.cex = 1, 
     edge.width = 0.5*E(vgg)$weight,
     edge.curved = TRUE,
     edge.color = E(vgg)$color)

# FILTER AND REPEAT:
vgg <- as.undirected(graph.adjacency(aid_vdc,weighted=TRUE))
V(vgg)$color <- rep("green",length(u_vdc))
V(vgg)$name <- u_vdc
cut85 <- quantile(as.vector(aid_vdc[aid_vdc>0]),0.85)
vgg_f <- filter(cutoff = cut85,
                edge_matrix = aid_vdc,
                vertex_colors = V(vgg)$color,
                vertex_names = V(vgg)$name)
vgg_f <- as.undirected(vgg_f)
vgg_f <- giant_comp(graph = vgg_f,
                    vertex_color = V(vgg_f)$color,
                    vertex_names = V(vgg_f)$name)
ec <- edge.betweenness(graph = vgg_f,
                       e = E(vgg_f), 
                       directed = FALSE,
                       weights = E(vgg_f)$weight)
ec_int <- as.integer(round(ec,0))
for (k in 1:length(ec_int)){
  E(vgg_f)$color[k] <- rev(heat.colors(1+as.integer(max(ec_int))))[as.integer(ec_int[k])+1]
}
plot(vgg_f,
     layout = layout.fruchterman.reingold(vgg_f, niter=200, area=2000*vcount(vgg_f)),
     vertex.color = "green",
     vertex.size = 2,
     vertex.label = NA, 
     vertex.label.color = "black",
     vertex.label.font = 1.5, 
     vertex.label.cex = 1, 
     edge.width = 0.5*E(vgg_f)$weight,
     edge.curved = TRUE,
     edge.color = E(vgg_f)$color)

# FILTER AND REPEAT:
vgg <- as.undirected(graph.adjacency(aid_vdc,weighted=TRUE))
V(vgg)$color <- rep("green",length(u_vdc))
V(vgg)$name <- u_vdc
cut95 <- quantile(as.vector(aid_vdc[aid_vdc>0]),0.95)
vgg_f <- filter(cutoff = cut95,
                edge_matrix = aid_vdc,
                vertex_colors = V(vgg)$color,
                vertex_names = V(vgg)$name)
vgg_f <- as.undirected(vgg_f)
vgg_f <- giant_comp(graph = vgg_f,
                    vertex_color = V(vgg_f)$color,
                    vertex_names = V(vgg_f)$name)
ec <- edge.betweenness(graph = vgg_f,
                       e = E(vgg_f), 
                       directed = FALSE,
                       weights = E(vgg_f)$weight)
ec_int <- as.integer(round(ec,0))
for (k in 1:length(ec_int)){
  E(vgg_f)$color[k] <- rev(heat.colors(1+as.integer(max(ec_int))))[as.integer(ec_int[k])+1]
}
plot(vgg_f,
     layout = layout.fruchterman.reingold(vgg_f, niter=200, area=2000*vcount(vgg_f)),
     vertex.color = "green",
     vertex.size = 2,
     vertex.label = NA, 
     vertex.label.color = "black",
     vertex.label.font = 1.5, 
     vertex.label.cex = 1, 
     edge.width = 0.5*E(vgg_f)$weight,
     edge.curved = TRUE,
     edge.color = E(vgg_f)$color)

# FILTER AND REPEAT:
vgg <- as.undirected(graph.adjacency(aid_vdc,weighted=TRUE))
V(vgg)$color <- rep("green",length(u_vdc))
V(vgg)$name <- u_vdc
cut99 <- quantile(as.vector(aid_vdc[aid_vdc>0]),0.99)
vgg_f <- filter(cutoff = cut99,
                edge_matrix = aid_vdc,
                vertex_colors = V(vgg)$color,
                vertex_names = V(vgg)$name)
vgg_f <- as.undirected(vgg_f)
vgg_f <- giant_comp(graph = vgg_f,
                    vertex_color = V(vgg_f)$color,
                    vertex_names = V(vgg_f)$name)
ec <- edge.betweenness(graph = vgg_f,
                       e = E(vgg_f), 
                       directed = FALSE,
                       weights = E(vgg_f)$weight)
ec_int <- as.integer(round(ec,0))
for (k in 1:length(ec_int)){
  E(vgg_f)$color[k] <- rev(heat.colors(1+as.integer(max(ec_int))))[as.integer(ec_int[k])+1]
}
plot(vgg_f,
     layout = layout.fruchterman.reingold(vgg_f, niter=200, area=2000*vcount(vgg_f)),
     vertex.color = "green",
     vertex.size = 2,
     vertex.label = NA, 
     vertex.label.color = "black",
     vertex.label.font = 1.5, 
     vertex.label.cex = 1, 
     edge.width = 0.5*E(vgg_f)$weight,
     edge.curved = TRUE,
     edge.color = E(vgg_f)$color)


# FILTER AND REPEAT:
vgg <- as.undirected(graph.adjacency(aid_vdc,weighted=TRUE))
V(vgg)$color <- rep("green",length(u_vdc))
V(vgg)$name <- u_vdc
cut99.9 <- quantile(as.vector(aid_vdc[aid_vdc>0]),0.999)
vgg_f <- filter(cutoff = cut99.9,
                edge_matrix = aid_vdc,
                vertex_colors = V(vgg)$color,
                vertex_names = V(vgg)$name)
vgg_f <- as.undirected(vgg_f)
vgg_f <- giant_comp(graph = vgg_f,
                    vertex_color = V(vgg_f)$color,
                    vertex_names = V(vgg_f)$name)
ec <- edge.betweenness(graph = vgg_f,
                       e = E(vgg_f), 
                       directed = FALSE,
                       weights = E(vgg_f)$weight)
ec_int <- as.integer(round(ec,0))
for (k in 1:length(ec_int)){
  E(vgg_f)$color[k] <- rev(heat.colors(1+as.integer(max(ec_int))))[as.integer(ec_int[k])+1]
}
plot(vgg_f,
     layout = layout.fruchterman.reingold(vgg_f, niter=200, area=2000*vcount(vgg_f)),
     vertex.color = "green",
     vertex.size = 10,
     vertex.label = V(vgg_f)$name, 
     vertex.label.color = "black",
     vertex.label.font = 1.5, 
     vertex.label.cex = 1.3, 
     edge.width = 0.9*E(vgg_f)$weight,
     edge.curved = TRUE,
     edge.color = E(vgg_f)$color)


# EDGE-BETWEENNESS COMMUNITY
# The edge betweenness score of an edge measures the number of shortest paths through it, see edge.betweenness for details. 
# The idea of the edge betweenness based community structure detection is that it is likely that edges connecting separate modules 
# have high edge betweenness as all the shortest paths from one module to another must traverse through them. 
# So if we gradually remove the edge with the highest edge betweenness score we will get a hierarchical map, a rooted tree, 
# called a dendrogram of the graph. The leafs of the tree are the individual vertices and the root of the tree represents the whole graph.
#
# edge.betweenness.community performs this algorithm by calculating the edge betweenness of the graph, 
# removing the edge with the highest edge betweenness score, then recalculating edge betweenness of the edges 
# and again removing the one with the highest score, etc.

# FILTER AND CLUSTER WITH CUTOFF = 0.75
vgg <- as.undirected(graph.adjacency(aid_vdc,weighted=TRUE))
V(vgg)$color <- rep("green",length(u_vdc))
V(vgg)$name <- u_vdc
cut75 <- quantile(as.vector(aid_vdc[aid_vdc>0]),0.75)
vgg_f <- filter(cutoff = cut75,
                edge_matrix = aid_vdc,
                vertex_colors = V(vgg)$color,
                vertex_names = V(vgg)$name)
vgg_f <- as.undirected(vgg_f)
vgg_f <- giant_comp(graph = vgg_f,
                    vertex_color = V(vgg_f)$color,
                    vertex_names = V(vgg_f)$name)
ebc <- edge.betweenness.community(graph = vgg_f)
plot(ebc,
     vgg_f, 
     vertex.size=2,
     edge.width=0.5*E(vgg_f)$weight,
     main="Example: ML Communities",
     vertex.label.cex=0.8,
     vertex.label=NA)

# FILTER AND CLUSTER WITH CUTOFF = 0.85
vgg <- as.undirected(graph.adjacency(aid_vdc,weighted=TRUE))
V(vgg)$color <- rep("green",length(u_vdc))
V(vgg)$name <- u_vdc
cut85 <- quantile(as.vector(aid_vdc[aid_vdc>0]),0.85)
vgg_f <- filter(cutoff = cut85,
                edge_matrix = aid_vdc,
                vertex_colors = V(vgg)$color,
                vertex_names = V(vgg)$name)
vgg_f <- as.undirected(vgg_f)
vgg_f <- giant_comp(graph = vgg_f,
                    vertex_color = V(vgg_f)$color,
                    vertex_names = V(vgg_f)$name)
ebc <- edge.betweenness.community(graph = vgg_f)
plot(ebc,
     vgg_f, 
     vertex.size=2,
     edge.width=0.5*E(vgg_f)$weight,
     main="Example: ML Communities",
     vertex.label.cex=0.8,
     vertex.label=NA)

# FILTER AND CLUSTER WITH CUTOFF = 0.95
vgg <- as.undirected(graph.adjacency(aid_vdc,weighted=TRUE))
V(vgg)$color <- rep("green",length(u_vdc))
V(vgg)$name <- u_vdc
cut95 <- quantile(as.vector(aid_vdc[aid_vdc>0]),0.95)
vgg_f <- filter(cutoff = cut95,
                edge_matrix = aid_vdc,
                vertex_colors = V(vgg)$color,
                vertex_names = V(vgg)$name)
vgg_f <- as.undirected(vgg_f)
vgg_f <- giant_comp(graph = vgg_f,
                    vertex_color = V(vgg_f)$color,
                    vertex_names = V(vgg_f)$name)
ebc <- edge.betweenness.community(graph = vgg_f)
plot(ebc,
     vgg_f, 
     vertex.size=3,
     edge.width=0.5*E(vgg_f)$weight,
     main="Example: ML Communities",
     vertex.label.cex=0.8,
     vertex.label=NA)

# FILTER AND CLUSTER WITH CUTOFF = 0.99
vgg <- as.undirected(graph.adjacency(aid_vdc,weighted=TRUE))
V(vgg)$color <- rep("green",length(u_vdc))
V(vgg)$name <- u_vdc
cut99 <- quantile(as.vector(aid_vdc[aid_vdc>0]),0.99)
vgg_f <- filter(cutoff = cut99,
                edge_matrix = aid_vdc,
                vertex_colors = V(vgg)$color,
                vertex_names = V(vgg)$name)
vgg_f <- as.undirected(vgg_f)
vgg_f <- giant_comp(graph = vgg_f,
                    vertex_color = V(vgg_f)$color,
                    vertex_names = V(vgg_f)$name)
ebc <- edge.betweenness.community(graph = vgg_f)
plot(ebc,
     vgg_f, 
     vertex.size=4,
     edge.width=0.05*E(vgg_f)$weight,
     main="Filtered VDC Aid Network Betweenness Community Structure",
     vertex.label.cex=0.8,
     vertex.label=NA)


# BETWEENNESS ESTIMATE: This measure calculates betweenness by considering only paths of a certain length 
# that is smaller than or equal to the cutoff value. Similarly for edge betweenness estimates. 
# In our analysis, we will use cutoff lengths 10-30 (refer to the path length distribution analysis above). 
# Here, we use a cutoff=3 just to illustrate the application.


# BETWEENNESS ESTIMATE CENTRALITY: THE NUMBER OF GEODESICS GOING THROUGH A NODE
vgg <- as.undirected(graph.adjacency(aid_vdc,weighted=TRUE))
V(vgg)$color <- rep("green",length(u_vdc))
V(vgg)$name <- u_vdc
bce <- betweenness.estimate(graph = vgg,
                            vids = V(vgg),
                            weights = E(vgg)$weight,
                            directed=FALSE,
                            4)
plot(sort(bce, decreasing=TRUE),
     col=adjustcolor(rgb(1/2,0,0,1/2)), 
     xlab="Node Index", 
     ylab="Betweenness Centrality Estimate", 
     main="Sorted Relief Agency Betweenness Centrality Values", 
     pch=19)
histP2(bce,
       breaks=100,
       col=adjustcolor(rgb(1/2,0,0,1/2)),
       xlab="Betweenness Estimate Centrality Values",
       main="Relief Agency Betweenness Estimate Centrality Distribution")

# PLOT HEAT MAP ON VERTICES ACCORDING TO BETWEENNESS CENTRALITY
bce_int <- as.integer(round(bce,0))
for (k in 1:length(bce_int)){
  V(vgg)$color[k] <- rev(heat.colors(1+as.integer(max(bce_int))))[as.integer(bce_int[k])+1]
}

# PLOT AGENCY GRAPH AND FILTER
plot(vgg,
     layout = layout.fruchterman.reingold(vgg, niter=200, area=2000*vcount(vgg)),
     vertex.color = V(vgg)$color,
     vertex.size = 2,
     vertex.label = NA, 
     vertex.label.color = "black",
     vertex.label.font = 1.5, 
     vertex.label.cex = 1, 
     edge.width = 0.5*E(vgg)$weight,
     edge.curved = TRUE,
     edge.color = gray.colors(1))

# FIND THE TOP 10% BETWEENNES NODES
top_bce <- bce[which(bce > quantile(bce,0.9))]
top_bce

# FIND THE TOP 5% BETWENNES NODES
top_bce <- bce[which(bce > quantile(bce,0.95))]
top_bce

# REMOVE ISOLATED
vgg <- drop_isolated(graph = vgg,
                     vertex_colors = V(vgg)$color,
                     vertex_names = V(vgg)$name)

# GET THE GIANT CONNECTED COMPONENT (TWO CLSUTERS ONLY)
vgg <- giant_comp(graph = vgg,
                  vertex_colors = V(vgg)$color,
                  vertex_names = V(vgg)$name)

# SET THE GRAPH COLOR ACCORDING TO BC
bce <- betweenness.estimate(graph = vgg,
                            vids = V(vgg),
                            weights = E(vgg)$weight,
                            directed=FALSE,
                            cutoff = 4)
bce_int <- as.integer(round(bce,0))
for (k in 1:length(bce_int)){
  V(vgg)$color[k] <- rev(heat.colors(1+as.integer(max(bce_int))))[as.integer(bce_int[k])+1]
}
plot(vgg,
     layout = layout.fruchterman.reingold(vgg, niter=200, area=2000*vcount(vgg)),
     vertex.color = V(vgg)$color,
     vertex.size = 5,
     vertex.label = NA, 
     vertex.label.color = "black",
     vertex.label.font = 1, 
     vertex.label.cex = 0.5, 
     edge.width = 0.1*E(vgg)$weight,
     edge.curved = TRUE,
     edge.color = gray.colors(1))


# EDGE BETWEENNESS ESTIMATE: This betweenness estimate is defined in a similar way as betweenness estimate, 
# we show it here with cutoff=4.
vgg <- as.undirected(graph.adjacency(aid_vdc,weighted=TRUE))
V(vgg)$color <- rep("green",length(u_vdc))
V(vgg)$name <- u_vdc
ebe <- edge.betweenness.estimate(graph = vgg,
                                 weights =E(vgg)$weight,
                                 directed=FALSE,
                                 cutoff = 4)
plot(sort(ebe, decreasing=TRUE),
     col=adjustcolor(rgb(1/2,0,0,1/2)), 
     xlab="Node Index", 
     ylab="Edge-Betweenness Centrality", 
     main="Sorted Relief Agency Edge-Betweenness Centrality Values", 
     pch=19)
histP2(ebe,
       breaks=100,
       col=adjustcolor(rgb(1/2,0,0,1/2)),
       xlab="Edge-Betweenness Centrality Values",
       main="Relief Agency Edge-Betweenness Centrality Distribution")

# PLOT HEAT MAP ON VERTICES ACCORDING TO EDGE-BEWEENNESS CENTRALITY
ebe_int <- as.integer(round(ebe,0))
for (k in 1:length(ebe_int)){
  E(vgg)$color[k] <- rev(heat.colors(1+as.integer(max(ebe_int))))[as.integer(ebe_int[k])+1]
}

# PLOT AGENCY GRAPH AND FILTER
plot(vgg,
     layout = layout.fruchterman.reingold(vgg, niter=200, area=2000*vcount(vgg)),
     vertex.color = "green",
     vertex.size = 3,
     vertex.label = NA, 
     vertex.label.color = "black",
     vertex.label.font = 1.5, 
     vertex.label.cex = 1, 
     edge.width = 0.5*E(vgg)$weight,
     edge.curved = TRUE,
     edge.color = E(vgg)$color)

# FIND THE TOP 10% EDGE-BETWEENNES EDGES
top_ebe <- ebe[which(ebe > quantile(ebe,0.9))]
E(vgg)[which(ebe %in% top_ebe)]

# FIND THE TOP 5% EDGE-BETWEENNES EDGES
top_ebe <- ebe[which(ebe > quantile(ebe,0.999))]
E(vgg)[which(ebe %in% top_ebe)]

# REMOVE ISOLATED
vgg <- drop_isolated(graph = vgg,
                     vertex_colors = V(vgg)$color,
                     vertex_names = V(vgg)$name)

# GET THE GIANT CONNECTED COMPONENT (TWO CLSUTERS ONLY)
vgg <- giant_comp(graph = vgg,
                  vertex_colors = V(vgg)$color,
                  vertex_names = V(vgg)$name)

# SET THE GRAPH COLOR ACCORDING TO EC
ebe <- edge.betweenness.estimate(graph = vgg,
                                 weights =E(vgg)$weight,
                                 directed=FALSE,
                                 cutoff = 4)
ebe_int <- as.integer(round(ebe,0))
for (k in 1:length(ebe_int)){
  E(vgg)$color[k] <- rev(heat.colors(1+as.integer(max(ebe_int))))[as.integer(ebe_int[k])+1]
}
plot(vgg,
     layout = layout.fruchterman.reingold(vgg, niter=200, area=2000*vcount(vgg)),
     vertex.color = "green",
     vertex.size = 2,
     vertex.label = NA, 
     vertex.label.color = "black",
     vertex.label.font = 1.5, 
     vertex.label.cex = 1, 
     edge.width = 3*E(vgg)$weight,
     edge.curved = TRUE,
     edge.color = E(vgg)$color)





# WE CAN FURTHER FILTER AND APPLY THIS DEPENDING ON MODELING GOALS




# CLOSENESS CENTRALITY: This measure takes into account the distribution of distances to other nodes from a given node. 
# It is defined as the reciprocal of the farness of a node, where farness is defined as the sum of its distances to all 
# other nodes. Closeness can be regarded as a measure of how long it will take to spread information 
# (or an efect of an event) from a node to all other nodes. To demonstrate this concept, we compute the closeness 
# centrality for the unweighted network.
vgg <- as.undirected(graph.adjacency(aid_vdc,weighted=TRUE))
V(vgg)$color <- rep("green",length(u_vdc))
V(vgg)$name <- u_vdc
cl<-clusters(vgg)
vgg1<-induced.subgraph(vgg, which(cl$membership == which.max(cl$csize)))
cc<-closeness(graph = vgg1,vids = V(vgg1),weights = E(vgg1)$weight)
plot(sort(cc/max(cc), decreasing=TRUE), col=adjustcolor(rgb(1/2,0,1,1)), xlab="Node Id in the Giant Conencted Component (gg1)", ylab="Normalized Closeness Centrality", main="Closeness Centrality for the Giant Component (gg1)")
hist(cc/max(cc),breaks=200,col=adjustcolor(rgb(1/2,0,1,1)),xlab="Normalized Closeness Centrality Values for gg1",main="Normalized Closeness Centrality Distribution for gg1")

# DISPLAY THE GRAPH WITH THE APPORPRIATE COLORS FOR THE VERTICES
vgg <- as.undirected(graph.adjacency(aid_vdc,weighted=TRUE))
V(vgg)$color <- rep("green",length(u_vdc))
V(vgg)$name <- u_vdc
cc <- closeness(graph = vgg,
              vids = V(vgg), 
              weights = E(vgg)$weight,
              normalized = TRUE)
cc_int <- as.integer(round(10000*cc,0))
cc_int <- cc_int-min(cc_int)

# FIND THE TOP 10% BETWEENNES NODES
top_cc <- cc[which(cc > quantile(cc,0.9))]
top_cc

# FIND THE TOP 10% BETWEENNES NODES
top_cc <- cc[which(cc > quantile(cc,0.95))]
top_cc

# FIND THE TOP 10% BETWEENNES NODES
top_cc <- cc[which(cc > quantile(cc,0.99))]
top_cc

# REMOVE ISOLATED
vgg <- drop_isolated(graph = vgg,
                     vertex_colors = V(vgg)$color,
                     vertex_names = V(vgg)$name)

# GET THE GIANT CONNECTED COMPONENT (TWO CLSUTERS ONLY)
vgg <- giant_comp(graph = vgg,
                  vertex_colors = V(vgg)$color,
                  vertex_names = V(vgg)$name)

# NOTE: MUST RECALCULATE CC AFER GRAPH TRANSFORMATIONS
cc <- closeness(graph = vgg,
                vids = V(vgg), 
                weights = E(vgg)$weight,
                normalized = TRUE)
cc_int <- as.integer(round(10000*cc,0))
cc_int <- cc_int-min(cc_int)

for (k in 1:length(cc_int)){
  V(vgg)$color[k] <- rev(heat.colors(1+as.integer(max(cc_int))))[as.integer(cc_int[k])+1]
}
plot(vgg,
     layout = layout.fruchterman.reingold(vgg, niter=200, area=2000*vcount(vgg)),
     vertex.color = V(vgg)$color,
     vertex.size = 2,
     vertex.label = NA, 
     vertex.label.color = "black",
     vertex.label.font = 1, 
     vertex.label.cex = 0.5, 
     edge.width = 0.1*E(vgg)$weight,
     edge.curved = TRUE,
     edge.color = gray.colors(1))




# FILTER AND REPEAT:
vgg <- as.undirected(graph.adjacency(aid_vdc,weighted=TRUE))
V(vgg)$color <- rep("green",length(u_vdc))
V(vgg)$name <- u_vdc
cut95 <- quantile(as.vector(aid_vdc[aid_vdc>0]),0.95)
vgg_f <- filter(cutoff = cut95,
                edge_matrix = aid_vdc,
                vertex_colors = V(vgg)$color,
                vertex_names = V(vgg)$name)
vgg_f <- as.undirected(vgg_f)
vgg_f <- giant_comp(graph = vgg_f,
                    vertex_color = V(vgg_f)$color,
                    vertex_names = V(vgg_f)$name)
cc <- closeness(graph = vgg_f,
              vids = V(vgg_f), 
              weights = E(vgg_f)$weight,
              normalized = TRUE)
cc_int <- as.integer(round(10000*cc,0))
cc_int <- cc_int-min(cc_int)
for (k in 1:length(cc_int)){
  V(vgg_f)$color[k] <- rev(heat.colors(1+as.integer(max(bc_int))))[as.integer(bc_int[k])+1]
}
plot(vgg_f,
     layout = layout.fruchterman.reingold(vgg_f, niter=200, area=2000*vcount(vgg_f)),
     vertex.color = V(vgg_f)$color,
     vertex.size = 4,
     vertex.label = NA, 
     vertex.label.color = "black",
     vertex.label.font = 1, 
     vertex.label.cex = 0.5, 
     edge.width = 0.1*E(vgg_f)$weight,
     edge.curved = TRUE,
     edge.color = gray.colors(1))


# FILTER AND REPEAT:
vgg <- as.undirected(graph.adjacency(aid_vdc,weighted=TRUE))
V(vgg)$color <- rep("green",length(u_vdc))
V(vgg)$name <- u_vdc
cut99 <- quantile(as.vector(aid_vdc[aid_vdc>0]),0.99)
vgg_f <- filter(cutoff = cut99,
                edge_matrix = aid_vdc,
                vertex_colors = V(vgg)$color,
                vertex_names = V(vgg)$name)
vgg_f <- as.undirected(vgg_f)
vgg_f <- giant_comp(graph = vgg_f,
                    vertex_color = V(vgg_f)$color,
                    vertex_names = V(vgg_f)$name)
cc <- closeness(graph = vgg_f,
              vids = V(vgg_f), 
              weights = E(vgg_f)$weight,
              normalized = TRUE)
cc_int <- as.integer(round(10000*cc,0))
cc_int <- cc_int-min(cc_int)
for (k in 1:length(cc_int)){
  V(vgg_f)$color[k] <- rev(heat.colors(1+as.integer(max(bc_int))))[as.integer(bc_int[k])+1]
}
plot(vgg_f,
     layout = layout.fruchterman.reingold(vgg_f, niter=200, area=2000*vcount(vgg_f)),
     vertex.color = V(vgg_f)$color,
     vertex.size = 4,
     vertex.label = NA, 
     vertex.label.color = "black",
     vertex.label.font = 1, 
     vertex.label.cex = 0.5, 
     edge.width = 0.1*E(vgg_f)$weight,
     edge.curved = TRUE,
     edge.color = gray.colors(1))




cce<-closeness.estimate(graph = vgg1,vids = V(vgg1),weights = E(vgg1)$weight,cutoff = 7)
plot(sort(cce/max(cce), decreasing=TRUE), col=adjustcolor(rgb(1/2,0,1,1)), xlab="Node Id in the Giant Conencted Component (gg1)", ylab="Normalized Closeness Centrality", main="Closeness Centrality for the Giant Component (gg1)")
hist(cce/max(cce),breaks=200,col=adjustcolor(rgb(1/2,0,1,1)),xlab="Normalized Closeness Centrality Values for gg1",main="Normalized Closeness Centrality Distribution for gg1")


# WE CAN FILTER FURTHER AT A LATER STAGE



# EIGENVECTOR CENTRALITY: This is a measure of the influence of a node in the network. 
# It assigns relative scores to all nodes in the network based on the concept that connections to high-scoring nodes 
# contribute more to the score of the given node than equal conenctions to low-scoring nodes. 
# A variant of egenvector centrality is Google's PageRank algorithm.
vgg <- as.undirected(graph.adjacency(aid_vdc,weighted=TRUE))
V(vgg)$color <- rep("green",length(u_vdc))
V(vgg)$name <- u_vdc
clu<-clusters(vgg)
vgg1<-induced.subgraph(vgg, which(clu$membership == which.max(clu$csize)))
ec <- evcent(graph = vgg1)$vector
plot(sort(ec, decreasing=TRUE), col=adjustcolor(rgb(0,0,1,1/2)), xlab="Node Id in the Network (1:200)", ylab="Closeness Centrality Values", main="Essential (first 200 nodes) Closeness Centrality for g", pch=20)
hist(ec,breaks=100,col=adjustcolor(rgb(0,0,1,1/2)),xlab="Closeness Centrality Values",main="Essential Closeness Centrality Distribution")


# DISPLAY THE GRAPH WITH THE APPORPRIATE COLORS FOR THE VERTICES
vgg <- as.undirected(graph.adjacency(aid_vdc,weighted=TRUE))
V(vgg)$color <- rep("green",length(u_vdc))
V(vgg)$name <- u_vdc
ec <- evcent(graph = vgg)$vector
ec_int <- as.integer(round(1000*ec,0))
ec_int <- ec_int-min(ec_int)

# FIND THE TOP 10% BETWEENNES NODES
top_ec <- ec[which(ec > quantile(ec,0.9))]
top_ec

# FIND THE TOP 10% BETWEENNES NODES
top_ec <- ec[which(ec > quantile(ec,0.95))]
top_ec

# FIND THE TOP 10% BETWEENNES NODES
top_ec <- ec[which(ec > quantile(ec,0.99))]
top_ec

# REMOVE ISOLATED
vgg <- drop_isolated(graph = vgg,
                     vertex_colors = V(vgg)$color,
                     vertex_names = V(vgg)$name)

# GET THE GIANT CONNECTED COMPONENT (TWO CLSUTERS ONLY)
vgg <- giant_comp(graph = vgg,
                  vertex_colors = V(vgg)$color,
                  vertex_names = V(vgg)$name)

# NOTE: MUST RECALCULATE EC AFTER GRAPH TRANSFORMATIONS
ec <- evcent(graph = vgg)$vector
ec_int <- as.integer(round(1000*ec,0))
ec_int <- ec_int-min(ec_int)

for (k in 1:length(ec_int)){
  V(vgg)$color[k] <- rev(heat.colors(1+as.integer(max(ec_int))))[as.integer(ec_int[k])+1]
}
plot(vgg,
     layout = layout.fruchterman.reingold(vgg, niter=200, area=2000*vcount(vgg)),
     vertex.color = V(vgg)$color,
     vertex.size = 2,
     vertex.label = NA, 
     vertex.label.color = "black",
     vertex.label.font = 1, 
     vertex.label.cex = 0.5, 
     edge.width = 0.1*E(vgg)$weight,
     edge.curved = TRUE,
     edge.color = gray.colors(1))

# FILTER AND REPEAT:
vgg <- as.undirected(graph.adjacency(aid_vdc,weighted=TRUE))
V(vgg)$color <- rep("green",length(u_vdc))
V(vgg)$name <- u_vdc
cut95 <- quantile(as.vector(aid_vdc[aid_vdc>0]),0.95)
vgg_f <- filter(cutoff = cut95,
                edge_matrix = aid_vdc,
                vertex_colors = V(vgg)$color,
                vertex_names = V(vgg)$name)
vgg_f <- as.undirected(vgg_f)
vgg_f <- giant_comp(graph = vgg_f,
                    vertex_color = V(vgg_f)$color,
                    vertex_names = V(vgg_f)$name)
ec <- evcent(graph = vgg_f)$vector
ec_int <- as.integer(round(1000*ec,0))
ec_int <- ec_int-min(ec_int)

for (k in 1:length(ec_int)){
  V(vgg_f)$color[k] <- rev(heat.colors(1+as.integer(max(bc_int))))[as.integer(bc_int[k])+1]
}
plot(vgg_f,
     layout = layout.fruchterman.reingold(vgg_f, niter=200, area=2000*vcount(vgg_f)),
     vertex.color = V(vgg_f)$color,
     vertex.size = 4,
     vertex.label = NA, 
     vertex.label.color = "black",
     vertex.label.font = 1, 
     vertex.label.cex = 0.5, 
     edge.width = 0.1*E(vgg_f)$weight,
     edge.curved = TRUE,
     edge.color = gray.colors(1))


# FILTER AND REPEAT:
vgg <- as.undirected(graph.adjacency(aid_vdc,weighted=TRUE))
V(vgg)$color <- rep("green",length(u_vdc))
V(vgg)$name <- u_vdc
cut99 <- quantile(as.vector(aid_vdc[aid_vdc>0]),0.99)
vgg_f <- filter(cutoff = cut99,
                edge_matrix = aid_vdc,
                vertex_colors = V(vgg)$color,
                vertex_names = V(vgg)$name)
vgg_f <- as.undirected(vgg_f)
vgg_f <- giant_comp(graph = vgg_f,
                    vertex_color = V(vgg_f)$color,
                    vertex_names = V(vgg_f)$name)
ec <- evcent(graph = vgg)$vector
ec_int <- as.integer(round(1000*ec,0))
ec_int <- ec_int-min(ec_int)
for (k in 1:length(ec_int)){
  V(vgg_f)$color[k] <- rev(heat.colors(1+as.integer(max(bc_int))))[as.integer(bc_int[k])+1]
}
plot(vgg_f,
     layout = layout.fruchterman.reingold(vgg_f, niter=200, area=2000*vcount(vgg_f)),
     vertex.color = V(vgg_f)$color,
     vertex.size = 4,
     vertex.label = NA, 
     vertex.label.color = "black",
     vertex.label.font = 1, 
     vertex.label.cex = 0.5, 
     edge.width = 0.1*E(vgg_f)$weight,
     edge.curved = TRUE,
     edge.color = gray.colors(1))




# AUTHORITY SCORE: This is a measure for DIRECTED NETWORKS, and it measures the number of nodes that are hubs and point 
# to a given node. It is defined as the principle eigenvector values for t(A)*A, where A stands for the adjacency 
# matrix of the network. For undirected networks like ours, the adjacency matrix is symmetric, so the authority score 
# is equivalent to the hub score. In subsequent analyses, we will be looking at directed extensions of this network 
# model, so we are including these two scores in the analysis for completeness.
vgg <- as.undirected(graph.adjacency(aid_vdc,weighted=TRUE))
V(vgg)$color <- rep("green",length(u_vdc))
V(vgg)$name <- u_vdc
au <- authority.score(graph = vgg,
                      weights = E(vgg)$weight)$vector
plot(sort(au, decreasing=TRUE), col=adjustcolor(rgb(0,0,1,1/2)), xlab="Node Id in the Network (1:200)", ylab="Authority Score Values", main="Essential (first 200 nodes) Authority Scores for g", pch=20)
hist(au,breaks=100,col=adjustcolor(rgb(0,0,1,1/2)),xlab="Authority Score Values",main="Essential Authority Score Distribution")


# DISPLAY THE GRAPH WITH THE APPORPRIATE COLORS FOR THE VERTICES
vgg <- as.undirected(graph.adjacency(aid_vdc,weighted=TRUE))
V(vgg)$color <- rep("green",length(u_vdc))
V(vgg)$name <- u_vdc
au <- authority.score(graph = vgg,
                      weights = E(vgg)$weight)$vector
au_int <- as.integer(round(10000*au,0))
au_int <- au_int-min(au_int)

# FIND THE TOP 10% BETWEENNES NODES
top_au <- au[which(au > quantile(au,0.9))]
top_au

# FIND THE TOP 10% BETWEENNES NODES
top_au <- au[which(au > quantile(au,0.95))]
top_au

# FIND THE TOP 10% BETWEENNES NODES
top_au <- au[which(au > quantile(au,0.99))]
top_au

# REMOVE ISOLATED
vgg <- drop_isolated(graph = vgg,
                     vertex_colors = V(vgg)$color,
                     vertex_names = V(vgg)$name)

# GET THE GIANT CONNECTED COMPONENT (TWO CLSUTERS ONLY)
vgg <- giant_comp(graph = vgg,
                  vertex_colors = V(vgg)$color,
                  vertex_names = V(vgg)$name)

# NOTE: MUST RECALCULATE au AFER GRAPH TRANSFORMATIONS
au <- authority.score(graph = vgg,
                      weights = E(vgg)$weight)$vector
au_int <- as.integer(round(10000*au,0))
au_int <- au_int-min(au_int)

for (k in 1:length(au_int)){
  V(vgg)$color[k] <- rev(heat.colors(1+as.integer(max(au_int))))[as.integer(au_int[k])+1]
}
plot(vgg,
     layout = layout.fruchterman.reingold(vgg, niter=200, area=2000*vcount(vgg)),
     vertex.color = V(vgg)$color,
     vertex.size = 2,
     vertex.label = NA, 
     vertex.label.color = "black",
     vertex.label.font = 1, 
     vertex.label.cex = 0.5, 
     edge.width = 0.1*E(vgg)$weight,
     edge.curved = TRUE,
     edge.color = gray.colors(1))

# FILTER AND REPEAT:
vgg <- as.undirected(graph.adjacency(aid_vdc,weighted=TRUE))
V(vgg)$color <- rep("green",length(u_vdc))
V(vgg)$name <- u_vdc
cut95 <- quantile(as.vector(aid_vdc[aid_vdc>0]),0.95)
vgg_f <- filter(cutoff = cut95,
                edge_matrix = aid_vdc,
                vertex_colors = V(vgg)$color,
                vertex_names = V(vgg)$name)
vgg_f <- as.undirected(vgg_f)
vgg_f <- giant_comp(graph = vgg_f,
                    vertex_color = V(vgg_f)$color,
                    vertex_names = V(vgg_f)$name)
au <- authority.score(graph = vgg_f,
                      weights = E(vgg_f)$weight)$vector
au_int <- as.integer(round(10000*au,0))
au_int <- au_int-min(au_int)
for (k in 1:length(au_int)){
  V(vgg_f)$color[k] <- rev(heat.colors(1+as.integer(max(bc_int))))[as.integer(bc_int[k])+1]
}
plot(vgg_f,
     layout = layout.fruchterman.reingold(vgg_f, niter=200, area=2000*vcount(vgg_f)),
     vertex.color = V(vgg_f)$color,
     vertex.size = 4,
     vertex.label = NA, 
     vertex.label.color = "black",
     vertex.label.font = 1, 
     vertex.label.cex = 0.5, 
     edge.width = 0.1*E(vgg_f)$weight,
     edge.curved = TRUE,
     edge.color = gray.colors(1))

# FILTER AND REPEAT:
vgg <- as.undirected(graph.adjacency(aid_vdc,weighted=TRUE))
V(vgg)$color <- rep("green",length(u_vdc))
V(vgg)$name <- u_vdc
cut99 <- quantile(as.vector(aid_vdc[aid_vdc>0]),0.99)
vgg_f <- filter(cutoff = cut99,
                edge_matrix = aid_vdc,
                vertex_colors = V(vgg)$color,
                vertex_names = V(vgg)$name)
vgg_f <- as.undirected(vgg_f)
vgg_f <- giant_comp(graph = vgg_f,
                    vertex_color = V(vgg_f)$color,
                    vertex_names = V(vgg_f)$name)
au <- authority.score(graph = vgg,
                      weights = E(vgg)$weight)$vector
au_int <- as.integer(round(10000*au,0))
au_int <- au_int-min(au_int)
for (k in 1:length(au_int)){
  V(vgg_f)$color[k] <- rev(heat.colors(1+as.integer(max(bc_int))))[as.integer(bc_int[k])+1]
}
plot(vgg_f,
     layout = layout.fruchterman.reingold(vgg_f, niter=200, area=2000*vcount(vgg_f)),
     vertex.color = V(vgg_f)$color,
     vertex.size = 4,
     vertex.label = NA, 
     vertex.label.color = "black",
     vertex.label.font = 1, 
     vertex.label.cex = 0.5, 
     edge.width = 0.1*E(vgg_f)$weight,
     edge.curved = TRUE,
     edge.color = gray.colors(1))



# APPLY TO AV NETWORK


# HUB SCORE: This is a measure FOR DIRECTED NETWORKS and it measures the number of authority nodes that a given hub node points to.
vgg <- as.undirected(graph.adjacency(aid_vdc,weighted=TRUE))
V(vgg)$color <- rep("green",length(ag))
V(vgg)$name <- u_vdc
hb<-hub.score(vgg)$vector
plot(sort(hb, decreasing=TRUE), col=adjustcolor(rgb(0,0,1,1/2)), xlab="Node Id in the Network (1:200)", ylab="Hub Score Values", main="Essential (first 200 nodes) Hub Scores for g", pch=20)
hist(hb,breaks=100,col=adjustcolor(rgb(0,0,1,1/2)),xlab="Hub Score Values",main="Essential Hub Score Distribution")

# APPLY TO AV NETWORK

# CLUSTERING COEFFICIENTS: This is a measure of the clustering of the network,defined by the ratio of the number of closed triplets 
# and the number of connected triplets of vertices. Here we include the local and weighted version 
# of the clustering coefficients. Clustering is particularly relevant to social netowrks where nodes tend to create tightly knit groups charaterized 
# by a high density of ties, this likelihood is greater than the average probability of an edge between two randomly selected nodes.
vgg <- as.undirected(graph.adjacency(aid_vdc,weighted=TRUE))
V(vgg)$color <- rep("green",length(u_vdc))
V(vgg)$name <- u_vdc
tr <- transitivity(graph = vgg,
                   vids = V(vgg), 
                   type="local")
plot(sort(tr), col=adjustcolor(rgb(0,0,1,1/2)), xlab="Node Id in the Network (1:3200)", ylab="Clustering Coefficient Values", main="Essential Clustering Coefficients for g", pch=20)
hist(tr,breaks=100,col=adjustcolor(rgb(0,0,1,1/2)),xlab="Clustering Coefficient Values",main="Clustering Coefficient Distribution for g")

# DISPLAY THE GRAPH WITH THE APPORPRIATE COLORS FOR THE VERTICES
vgg <- as.undirected(graph.adjacency(aid_vdc,weighted=TRUE))
V(vgg)$color <- rep("green",length(u_vdc))
V(vgg)$name <- u_vdc
# REMOVE ISOLATED
vgg <- drop_isolated(graph = vgg,
                     vertex_colors = V(vgg)$color,
                     vertex_names = V(vgg)$name)

# GET THE GIANT CONNECTED COMPONENT (TWO CLSUTERS ONLY)
vgg <- giant_comp(graph = vgg,
                  vertex_colors = V(vgg)$color,
                  vertex_names = V(vgg)$name)
# CALCULATE LOCAL TRANSITIVITY
tr <- transitivity(graph = vgg,
                   vids = V(vgg),
                   weights =  E(vgg)$weight,
                   type="local")
tr_int <- 100*tr
tr_int <- tr_int-min(tr_int)

for (k in 1:length(tr_int)){
  V(vgg)$color[k] <- rev(heat.colors(1+as.integer(max(tr_int))))[as.integer(tr_int[k])+1]
}
plot(vgg,
     layout = layout.fruchterman.reingold(vgg, niter=200, area=2000*vcount(vgg)),
     vertex.color = V(vgg)$color,
     vertex.size = 2,
     vertex.label = NA, 
     vertex.label.color = "black",
     vertex.label.font = 1, 
     vertex.label.cex = 0.5, 
     edge.width = 0.1*E(vgg)$weight,
     edge.curved = TRUE,
     edge.color = gray.colors(1))


# FILTER AND REPEAT:
vgg <- as.undirected(graph.adjacency(aid_vdc,weighted=TRUE))
V(vgg)$color <- rep("green",length(u_vdc))
V(vgg)$name <- u_vdc
cut95 <- quantile(as.vector(aid_vdc[aid_vdc>0]),0.95)
vgg_f <- filter(cutoff = cut95,
                edge_matrix = aid_vdc,
                vertex_colors = V(vgg)$color,
                vertex_names = V(vgg)$name)
vgg_f <- as.undirected(vgg_f)
vgg_f <- giant_comp(graph = vgg_f,
                    vertex_color = V(vgg_f)$color,
                    vertex_names = V(vgg_f)$name)
tr <- transitivity(graph = vgg_f,
                   vids = V(vgg_f),
                   weights =  E(vgg_f)$weight,
                   type="local")
tr_int <- 100*tr
tr_int <- tr_int-min(tr_int)
for (k in 1:length(tr_int)){
  V(vgg_f)$color[k] <- rev(heat.colors(1+as.integer(max(bc_int))))[as.integer(bc_int[k])+1]
}
plot(vgg_f,
     layout = layout.fruchterman.reingold(vgg_f, niter=200, area=2000*vcount(vgg_f)),
     vertex.color = V(vgg_f)$color,
     vertex.size = 4,
     vertex.label = NA, 
     vertex.label.color = "black",
     vertex.label.font = 1, 
     vertex.label.cex = 0.5, 
     edge.width = 0.1*E(vgg_f)$weight,
     edge.curved = TRUE,
     edge.color = gray.colors(1))


# FILTER AND REPEAT:
vgg <- as.undirected(graph.adjacency(aid_vdc,weighted=TRUE))
V(vgg)$color <- rep("green",length(u_vdc))
V(vgg)$name <- u_vdc
cut99 <- quantile(as.vector(aid_vdc[aid_vdc>0]),0.99)
vgg_f <- filter(cutoff = cut99,
                edge_matrix = aid_vdc,
                vertex_colors = V(vgg)$color,
                vertex_names = V(vgg)$name)
vgg_f <- as.undirected(vgg_f)
vgg_f <- giant_comp(graph = vgg_f,
                    vertex_color = V(vgg_f)$color,
                    vertex_names = V(vgg_f)$name)
tr <- transitivity(graph = vgg_f,
                   vids = V(vgg_f),
                   weights =  E(vgg_f)$weight,
                   type="local")
tr_int <- 100*tr
tr_int <- tr_int-min(tr_int)
for (k in 1:length(tr_int)){
  V(vgg_f)$color[k] <- rev(heat.colors(1+as.integer(max(bc_int))))[as.integer(bc_int[k])+1]
}
plot(vgg_f,
     layout = layout.fruchterman.reingold(vgg_f, niter=200, area=2000*vcount(vgg_f)),
     vertex.color = V(vgg_f)$color,
     vertex.size = 4,
     vertex.label = NA, 
     vertex.label.color = "black",
     vertex.label.font = 1, 
     vertex.label.cex = 0.5, 
     edge.width = 0.1*E(vgg_f)$weight,
     edge.curved = TRUE,
     edge.color = gray.colors(1))


# The weighted analogue of the clustering coefficient
trw<-transitivity(vgg, type="weighted")
plot(sort(trw), col=adjustcolor(rgb(0,0,1,1/2)), xlab="Node Id in the Network (1:200)", ylab="Clustering Coefficient Values", main="Essential Clustering Coefficients for g", pch=20)
hist(trw,breaks=100,col=adjustcolor(rgb(0,0,1,1/2)),xlab="Clustering Coefficient Values",main="Clustering Coefficient Distribution for g")

# DISPLAY THE GRAPH WITH THE APPORPRIATE COLORS FOR THE VERTICES
vgg <- as.undirected(graph.adjacency(aid_vdc,weighted=TRUE))
V(vgg)$color <- rep("green",length(u_vdc))
V(vgg)$name <- u_vdc
# REMOVE ISOLATED
vgg <- drop_isolated(graph = vgg,
                     vertex_colors = V(vgg)$color,
                     vertex_names = V(vgg)$name)
# GET THE GIANT CONNECTED COMPONENT (TWO CLSUTERS ONLY)
vgg <- giant_comp(graph = vgg,
                  vertex_colors = V(vgg)$color,
                  vertex_names = V(vgg)$name)
trw <- transitivity(graph = vgg,
                    type = "weighted",
                    vids = V(vgg))
trw_int <- as.integer(round(100*trw,0))
trw_int <- trw_int-min(trw_int)


trw_int <- trw
for (k in 1:length(trw_int)){
  V(vgg)$color[k] <- rev(heat.colors(1+as.integer(max(trw_int))))[as.integer(trw_int[k])+1]
}
plot(vgg,
     layout = layout.fruchterman.reingold(vgg, niter=200, area=2000*vcount(vgg)),
     vertex.color = V(vgg)$color,
     vertex.size = 2,
     vertex.label = NA, 
     vertex.label.color = "black",
     vertex.label.font = 1, 
     vertex.label.cex = 0.5, 
     edge.width = 0.1*E(vgg)$weight,
     edge.curved = TRUE,
     edge.color = gray.colors(1))

# FILTER AND REPEAT:
vgg <- as.undirected(graph.adjacency(aid_vdc,weighted=TRUE))
V(vgg)$color <- rep("green",length(u_vdc))
V(vgg)$name <- u_vdc
cut95 <- quantile(as.vector(aid_vdc[aid_vdc>0]),0.95)
vgg_f <- filter(cutoff = cut95,
                edge_matrix = aid_vdc,
                vertex_colors = V(vgg)$color,
                vertex_names = V(vgg)$name)
vgg_f <- as.undirected(vgg_f)
vgg_f <- giant_comp(graph = vgg_f,
                    vertex_color = V(vgg_f)$color,
                    vertex_names = V(vgg_f)$name)
trw <- transitivity(graph = vgg_f,
                    type = "weighted",
                    vids = V(vgg_f))
trw_int <- as.integer(round(100*trw,0))
trw_int <- trw_int-min(trw_int)
for (k in 1:length(trw_int)){
  V(vgg_f)$color[k] <- rev(heat.colors(1+as.integer(max(bc_int))))[as.integer(bc_int[k])+1]
}
plot(vgg_f,
     layout = layout.fruchterman.reingold(vgg_f, niter=200, area=2000*vcount(vgg_f)),
     vertex.color = V(vgg_f)$color,
     vertex.size = 4,
     vertex.label = NA, 
     vertex.label.color = "black",
     vertex.label.font = 1, 
     vertex.label.cex = 0.5, 
     edge.width = 0.1*E(vgg_f)$weight,
     edge.curved = TRUE,
     edge.color = gray.colors(1))


# FILTER AND REPEAT:
vgg <- as.undirected(graph.adjacency(aid_vdc,weighted=TRUE))
V(vgg)$color <- rep("green",length(u_vdc))
V(vgg)$name <- u_vdc
cut99 <- quantile(as.vector(aid_vdc[aid_vdc>0]),0.99)
vgg_f <- filter(cutoff = cut99,
                edge_matrix = aid_vdc,
                vertex_colors = V(vgg)$color,
                vertex_names = V(vgg)$name)
vgg_f <- as.undirected(vgg_f)
vgg_f <- giant_comp(graph = vgg_f,
                    vertex_color = V(vgg_f)$color,
                    vertex_names = V(vgg_f)$name)
trw <- transitivity(graph = vgg_f,
                    type = "weighted",
                    vids = V(vgg_f))
trw_int <- as.integer(round(100*trw,0))
trw_int <- trw_int-min(trw_int)
for (k in 1:length(trw_int)){
  V(vgg_f)$color[k] <- rev(heat.colors(1+as.integer(max(bc_int))))[as.integer(bc_int[k])+1]
}
plot(vgg_f,
     layout = layout.fruchterman.reingold(vgg_f, niter=200, area=2000*vcount(vgg_f)),
     vertex.color = V(vgg_f)$color,
     vertex.size = 4,
     vertex.label = NA, 
     vertex.label.color = "black",
     vertex.label.font = 1, 
     vertex.label.cex = 0.5, 
     edge.width = 0.1*E(vgg_f)$weight,
     edge.curved = TRUE,
     edge.color = gray.colors(1))




# COMMUNITY STRUCTURES: This is a way of performing funcitonal clustering in complex networks. 
# We have already looked at the connected components, 
# this is an elementary community detection based on connectivity.
vgg <- as.undirected(graph.adjacency(aid_vdc,weighted=TRUE))
V(vgg)$color <- rep("green",length(u_vdc))
V(vgg)$name <- u_vdc
cut99 <- quantile(as.vector(aid_vdc[aid_vdc>0]),0.998)
vgg_f <- filter(cutoff = cut99,
                edge_matrix = aid_vdc,
                vertex_colors = V(vgg)$color,
                vertex_names = V(vgg)$name)
vgg_f <- as.undirected(vgg_f)
strongclusters <- clusters(vgg_f)$membership
plot(vgg_f,
     vertex.color = strongclusters, 
     layout = layout.fruchterman.reingold,
     vertex.size = 4, 
     edge.color = "black", 
     edge.width = E(vgg_f)$weight,
     vertex.label = NA,
     main="Clustering for Store Network g200")




# ADD SOME FILTERING AND TRY AGAIN

mc<-multilevel.community(vgg)
plot(sort(mc$membership), col=adjustcolor(rgb(0,0,1,1/2)), xlab="Node Id in the Network", ylab="Multilevel Community Values", main="Multilevel Community Values for g", pch=20)
hist(mc$membership,breaks=100,col=adjustcolor(rgb(0,0,1,1/2)),xlab="Multilevel Community Values",main="Multilevel Community Distribution")


# MULTILEVEL COMMUNITY DETECTION:

mc<-multilevel.community(vgg)
plot(mc,vgg, vertex.size=2,edge.width=0.1*E(vgg)$weight,
     main="Example: ML Communities",
     vertex.label=NA)


# FILTER BY EDGE WEIGHT LEVELS

# FIRST FILTRATION CUTOFF LEVEL DOES NOT REALLY INDUCE ANY FILTRATION
vgg <- as.undirected(graph.adjacency(aid_vdc,weighted=TRUE))
V(vgg)$color <- rep("SkyBlue2",length(u_vdc))
V(vgg)$name <- u_vdc
cut90 <- quantile(as.vector(aid_vdc[aid_vdc>0]),0.90)
vgg_f <- filter(cutoff = cut90,
                edge_matrix = aid_vdc,
                vertex_colors = V(vgg)$color,
                vertex_names = V(vgg)$name)
vgg_f <- as.undirected(vgg_f)
vgg_f <- giant_comp(graph = vgg_f,
                    vertex_colors = V(vgg_f)$color,
                    vertex_names = V(vgg_f)$name)

mc_f <- multilevel.community(as.undirected(vgg_f))
plot(mc_f,
     as.undirected(vgg_f), 
     vertex.size=2,
     edge.width=sqrt(E(vgg_f)$weight),
     main="Example: ML Communities",
     vertex.label=NA)


# FURTHER FILTRATION, HUGE JUMP
vgg <- as.undirected(graph.adjacency(aid_vdc,weighted=TRUE))
V(vgg)$color <- rep("SkyBlue2",length(u_vdc))
V(vgg)$name <- u_vdc
cut97 <- quantile(as.vector(aid_vdc[aid_vdc>0]),0.97)
vgg_f <- filter(cutoff = cut97,
                edge_matrix = aid_vdc,
                vertex_colors = V(vgg)$color,
                vertex_names = V(vgg)$name)
vgg_f <- as.undirected(vgg_f)
vgg_f <- giant_comp(graph = vgg_f,
                    vertex_colors = V(vgg_f)$color,
                    vertex_names = V(vgg_f)$name)

mc_f <- multilevel.community(as.undirected(vgg_f))
plot(mc_f,
     as.undirected(vgg_f), 
     vertex.size=3,
     edge.width=sqrt(E(vgg_f)$weight),
     main="Example: ML Communities",
     vertex.label=NA)


# FURTHER FILTRATION
vgg <- as.undirected(graph.adjacency(aid_vdc,weighted=TRUE))
V(vgg)$color <- rep("SkyBlue2",length(u_vdc))
V(vgg)$name <- u_vdc
cut99.8 <- quantile(as.vector(aid_vdc[aid_vdc>0]),0.998)
vgg_f <- filter(cutoff = cut99.8,
                edge_matrix = aid_vdc,
                vertex_colors = V(vgg)$color,
                vertex_names = V(vgg)$name)
vgg_f <- as.undirected(vgg_f)
vgg_f <- giant_comp(graph = vgg_f,
                    vertex_colors = V(vgg_f)$color,
                    vertex_names = V(vgg_f)$name)

mc_f <- multilevel.community(as.undirected(vgg_f))
plot(mc_f,
     as.undirected(vgg_f), 
     vertex.size = 7,
     edge.width=sqrt(E(vgg_f)$weight),
     main="Example: ML Communities",
     vertex.label = V(vgg_f)$name)






# Next, we show the walktrap community algorithm.
vgg <- as.undirected(graph.adjacency(aid_vdc,weighted=TRUE))
V(vgg)$color <- rep("SkyBlue2",length(u_vdc))
V(vgg)$name <- u_vdc
wc <- walktrap.community(graph = vgg,
                       weights = E(vgg)$weight)
plot(sort(wc$membership), col=adjustcolor(rgb(0,0,1,1/2)), xlab="Node Id in the Network", ylab="Walktrap Community Values", main="Walktrap Community Values for g", pch=20)
hist(wc$membership,breaks=100,col=adjustcolor(rgb(0,0,1,1/2)),xlab="Walktrap Community Values",main="Walktrap Community Distribution")

# MULTILEVEL COMMUNITY DETECTION:
vgg <- as.undirected(graph.adjacency(aid_vdc,weighted=TRUE))
V(vgg)$color <- rep("SkyBlue2",length(u_vdc))
V(vgg)$name <- u_vdc
wc <- walktrap.community(graph = vgg,
                         weights = E(vgg)$weight)
plot(wc,vgg, vertex.size=2,edge.width=0.1*E(vgg)$weight,
     main="Example: ML Communities",
     vertex.label=NA)


# FILTER BY EDGE WEIGHT LEVELS

# FIRST FILTRATION CUTOFF LEVEL DOES NOT REALLY INDUCE ANY FILTRATION
vgg <- as.undirected(graph.adjacency(aid_vdc,weighted=TRUE))
V(vgg)$color <- rep("SkyBlue2",length(u_vdc))
V(vgg)$name <- u_vdc
cut90 <- quantile(as.vector(aid_vdc[aid_vdc>0]),0.90)
vgg_f <- filter(cutoff = cut90,
                edge_matrix = aid_vdc,
                vertex_colors = V(vgg)$color,
                vertex_names = V(vgg)$name)
vgg_f <- as.undirected(vgg_f)
vgg_f <- giant_comp(graph = vgg_f,
                    vertex_colors = V(vgg_f)$color,
                    vertex_names = V(vgg_f)$name)
wc <- walktrap.community(graph = vgg_f,
                         weights = E(vgg_f)$weight)
plot(wc_f,
     as.undirected(vgg_f), 
     vertex.size=2,
     edge.width=sqrt(E(vgg_f)$weight),
     main="Example: ML Communities",
     vertex.label=NA)


# FURTHER FILTRATION, HUGE JUMP
vgg <- as.undirected(graph.adjacency(aid_vdc,weighted=TRUE))
V(vgg)$color <- rep("SkyBlue2",length(u_vdc))
V(vgg)$name <- u_vdc
cut97 <- quantile(as.vector(aid_vdc[aid_vdc>0]),0.97)
vgg_f <- filter(cutoff = cut97,
                edge_matrix = aid_vdc,
                vertex_colors = V(vgg)$color,
                vertex_names = V(vgg)$name)
vgg_f <- as.undirected(vgg_f)
vgg_f <- giant_comp(graph = vgg_f,
                    vertex_colors = V(vgg_f)$color,
                    vertex_names = V(vgg_f)$name)
wc <- walktrap.community(graph = vgg_f,
                         weights = E(vgg_f)$weight)
plot(wc_f,
     as.undirected(vgg_f), 
     vertex.size=3,
     edge.width=sqrt(E(vgg_f)$weight),
     main="Example: ML Communities",
     vertex.label=NA)


# FURTHER FILTRATION
vgg <- as.undirected(graph.adjacency(aid_vdc,weighted=TRUE))
V(vgg)$color <- rep("SkyBlue2",length(u_vdc))
V(vgg)$name <- u_vdc
cut99.8 <- quantile(as.vector(aid_vdc[aid_vdc>0]),0.998)
vgg_f <- filter(cutoff = cut99.8,
                edge_matrix = aid_vdc,
                vertex_colors = V(vgg)$color,
                vertex_names = V(vgg)$name)
vgg_f <- as.undirected(vgg_f)
vgg_f <- giant_comp(graph = vgg_f,
                    vertex_colors = V(vgg_f)$color,
                    vertex_names = V(vgg_f)$name)
wc <- walktrap.community(graph = vgg_f,
                         weights = E(vgg_f)$weight)
plot(wc_f,
     as.undirected(vgg_f), 
     vertex.size = 7,
     edge.width=sqrt(E(vgg_f)$weight),
     main="Example: ML Communities",
     vertex.label = V(vgg_f)$name)














# FOR THE AGENCIES: SIMILARITY MEASURE BASED ON AID TYPE AND QUANTITY


# COMPARE VDC PROJECTION DISPLACEMENT TRACKING WITH SEVERITY INDEX
# COMPARE AID DISTRIBUTION WITH SEVERITY INDEX
# COMPARE DISPLACEMENT TRACKING WITH 





CLIQUE ANALYSIS:
```{r, echo=FALSE,results='markup',warning=FALSE,message=FALSE,fig.width=12,fig.height=6, dpi=200,out.width='1200px',out.height='600px'}
# Largest clique
largest.cliques(g1)[1]
# Maximal clique
mc<-maximal.cliques(g1)
mc[length(mc)] 
maximal.cliques.count(g1)
clique.number(g1)
````
SMALL-WORLD ANALYSIS:
```{r, echo=FALSE,results='markup',warning=FALSE,message=FALSE}
cat("We compute the giant component.")
cl<-clusters(g1)
gg1<-induced.subgraph(g1, which(cl$membership == which.max(cl$csize)))
cat("We measure average path length for the giant coponent.")
disthist<-path.length.hist(gg1, directed=TRUE)$res
diameter<-length(disthist)
avdistg<-weighted.mean(1:diameter, disthist)
cat("We compute the tail component.")
tg1<-induced.subgraph(g1, which(cl$membership != which.max(cl$csize)))
cat("We measure avearge path length for the tail component.")
disthist<-path.length.hist(tg1, directed=TRUE)$res
diameter<-length(disthist)
avdist<-weighted.mean(1:diameter, disthist)
cat("We measure the clustering coefficient.")
clustering_coef<-transitivity(gg1, type="global")
cat("We create erdos&renyi random network.")	
er<-erdos.renyi.game(vcount(gg1),ecount(gg1), type="gnm",directed=TRUE)
# We measure avearge path length for the erdos-renyi network.")
disthist_er<-path.length.hist(er, directed=TRUE)$res
diameter_er<-length(disthist_er)
avdist_er<-weighted.mean(1:diameter_er, disthist_er)     
cat("We measure the clustering coefficient for the erdos-renyi network for comparison.")
clustering_coef_er<-transitivity(er, type="global")
sindex<-(clustering_coef/clustering_coef_er) / (avdistg/avdist_er)
cat("We compute the small-world index of the store network. We will take a closer look at this and other properties in subsequent sections.")
sindex
```
library(plyr)
#We next compare the dept-level multi-basket (more than one dept) size distribution for renewal versus nonrenewal members after the first year. Since the distributions are exponential and hard to compare, we use the slope of the linear fit to the semi-log distribution to better see the difference. 
dy2<-as.data.frame(unique(cbind(dy1$visit_nbr,dy1$dept_nbr)))
dn2<-as.data.frame(unique(cbind(dn1$visit_nbr,dn1$dept_nbr)))
colnames(dy2)<-c("visit_nbr","dept_nbr")
colnames(dn2)<-c("visit_nbr","dept_nbr")
bsizey<-length(unique(dy2$visit_nbr))
bsizen<-length(unique(dn2$visit_nbr))
yvisits<-unique(dy2$visit_nbr)
nvisits<-unique(dn2$visit_nbr)
dy2$visit_num<-mapvalues(dy2$visit_nbr, from=yvisits,to=1:length(yvisits))
dn2$visit_num<-mapvalues(dn2$visit_nbr, from=nvisits,to=1:length(nvisits))
baskety<-vector()
basketn<-vector()
ty<-table(dy2$visit_num)
tn<-table(dn2$visit_num)
for (k in 1:bsizey){
  baskety[k]<-ty[k][[1]] 
}
for (k in 1:bsizen){
  basketn[k]<-tn[k][[1]] 
}
tby<-vector()
tbn<-vector()
for (k in 1:length(table(baskety))){
  tby[k]<-table(baskety)[k][[1]] 
}
for (k in 1:length(table(basketn))){
  tbn[k]<-table(basketn)[k][[1]] 
}
dy<-cbind.data.frame(as.vector(1:length(tby)),log(tby))
colnames(dy)<-c("x","y")
dn<-cbind.data.frame(as.vector(1:length(tbn)),log(tbn))
colnames(dn)<-c("x","y")
fity<-lm(y~x, data=dy)
fitn<-lm(y~x, data=dn)
par(mfrow=c(1,2))
hist(baskety[baskety>1], breaks=90,col=adjustcolor(rgb(1,0,1,1)),xlab="Sams Renewal Dept Index", main="Reneweal Multi-Basket Sizes at Dept Level",labels=ifelse(as.numeric(as.vector(hist(baskety[baskety>1],breaks=90,plot=FALSE)$counts))<1,"",hist(baskety[baskety>1],breaks=90,plot=FALSE)$counts))
hist(basketn[basketn>1], breaks=80,col=adjustcolor(rgb(0,0,1,1/2)),xlab="Sams Non-Renewal Dept Index", main="Non-Reneweal Multi-Basket Sizes at Dept Level",labels=ifelse(as.numeric(as.vector(hist(basketn[basketn>1],breaks=80,plot=FALSE)$counts))<1,"",hist(basketn[basketn>1],breaks=80,plot=FALSE)$counts))
par(mfrow=c(1,2))
plot(dy$x,dy$y,xlab="Sams Renewal Dept Index", ylab="Log[Basket Size Frequency]", pch=16, main="Linear Fit for Semi-Log Renewal Multi-Baskets")
lines(abline(coef=coef(fity), col="red"),xlim=range(0:length(dy)),ylim=range(0:max(log(dy))))
legend(15,7.5,legend=rbind(c("Intercept:","Slope:"),round(coef(fity),5)))
plot(dn$x,dn$y,xlab="Sams Non-Renewal Dept Index", ylab="Log[Basket Size Frequency]", pch=16, main="Linear Fit for Semi-Log Non-Renewal Multi-Baskets")
lines(abline(coef=coef(fitn), col="red"),xlim=range(0:length(dn)),ylim=range(0:max(log(dn))))




# TO USE HERE FOR EACH PROJECTION:






#
#
#
#
#
#
#
# CENTRALITY ANALYSIS OF THE NETWORK
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
# EIGENVECTOR CENTRALITY: THIS IS A MEASURE OF THE INFLUENCE OF A NODE IN A NETWORK.
# IT ASSIGNS RELATIVE SCORES TO ALL NODES BASED ON THE CONCEPT THAT CONNECTIONS
# TO HIGH-SCORING NODES CONTRIBUTE MORE TO THE SCORE OF A GIVEN NODE THAN CONNECTIONS
# TO LOW-SCORING NODES. A VARIANT OF EIGNEVECTOR CENTRALITY IS GOOGLE'S PAGERANK.
#
#
# NOTE: THIS IS LIKELY NOT THE RIGHT TOOL TO APPLY HERE aT THE VDC LEVEL
# BUT IF WE OBTAIN MORE GRANUALR DATA, WE WILL BE ABLE TO USE IT
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
    V(av)$name[k] <- all[k]
  } else {
    V(av)$size[k] <- 2
    V(av)$name[k] <- all[k]}
}
av <- drop_isolated(graph = av,
                    vertex_colors = V(av)$color,
                    vertex_names = V(av)$name)
av <- drop_loops(graph = av,
                 vertex_colors = V(av)$color,
                 vertex_names = V(av)$name)
ec <- evcent(av,
             directed = TRUE,
             weights = E(av)$weight)$vector
ec <- 10*ec
plot(sort(ec, decreasing = TRUE),
     col = adjustcolor(rgb(1,0,1,1)), 
     xlab = "Node Index", 
     ylab = "Eigenvector Centrality", 
     main = "Sorted VDC Network Eigenvector Centrality Values", 
     pch = 19)
histP2(ec,
       breaks = 50,
       col = adjustcolor(rgb(1,0,1,1)),
       xlab = "Eigenvector Centrality Values",
       main = "VDC Network Eigenvector Centrality Distribution")


# PLOT HEAT MAP ON VERTICES ACCORDING TO EIGENVECTOR CENTRALITY
for (k in 1:length(ec)){
  V(av)$color[k] <- rev(heat.colors(1+as.integer(max(ec))))[as.integer(ec[k])+1]
}
av_coords <- koords[which(vd %in% V(av)$name),]
plot(av,
     layout = layout.fruchterman.reingold(av, 
                                          niter = 200, 
                                          area = 2000*vcount(av)),
     vertex.color = V(av)$color,
     vertex.size = 9,
     vertex.label = V(av)$name, 
     vertex.label.color = "black",
     vertex.label.font = 1, 
     vertex.label.cex = 0.85, 
     edge.width = 0.2*sqrt(E(av)$weight),
     edge.arrow.size = 0.5,
     edge.curved = TRUE,
     edge.color = gray.colors(1),
     main = "Weighted Agency-VDC Aid Network Flow (VDC Level)")
legend("topleft",
       c("Highest Eigenvector Centrality","Lowest Eigenvector Centrality"),
       fill = c("red","White"),
       bty = "n")

# PLOT THE HEAT MAP OF THE AGENCY-VDC AID NETWORK
plot(av,
     layout = av_coords,
     vertex.color = V(av)$color,
     vertex.size = 4,
     vertex.label = NA, 
     vertex.label.color = "black",
     vertex.label.font = 1, 
     vertex.label.cex = 0.75, 
     edge.width = 0.2*sqrt(E(av)$weight),
     edge.arrow.size = 0.5,
     edge.curved = TRUE,
     edge.color = gray.colors(1),
     main = "Weighted Agency-VDC Aid Geo-Network Flow (VDC Level)")
legend("topright",
       c("Highest Eigenvector Centrality","Lowest Eigenvector Centrality"),
       fill = c("red","White"),
       bty = "n")
