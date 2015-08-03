# Nepal Disaster Relief Distribution and Displacement Tracking Network Analysis
# author: Georgi D. Gospodinov
# date: "July 30, 2015"
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
filter <- function(cutoff,edge_matrix,vertex_colors,vertex_names) {
  
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
giant_comp <- function(graph, vertex_colors, vertex_names){
  
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
  
  return(g_f)
}


# LOAD AGENCY RELIEF DATA:
aid_data <- read.csv("agency_relief.csv",sep=",")


# LOAD LAT/LON COORDINATES (OF CENTROIDS) AND HLCIT CODES
hlcit <- read.csv("master_hlcit.csv")
colnames(hlcit) <- c("lon","lat","vdc_name","vname","hlcit_code")
hlcit$hlcit_code <- as.factor(hlcit$hlcit_code)
hlcit$vname <- as.character(hlcit$vname)
hlcit$vdc_name <- as.character(hlcit$vdc_name)
hlcit <- rm_space(hlcit,"hlcit_code")
hlcit$hlcit_code <- as.numeric(levels(hlcit$hlcit_code))[hlcit$hlcit_code]



# AGENCY RELIEF NETWORK AT VDC LEVEL BELOW:


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

# SELECT UNIQUE AGENCIES AND TARGET VDC
ag <- unique(aid_data$impl_ag)
vd <- unique(aid_data$vdc)
all <- union(ag,vd)


vdc <- unique(aid_data$vdc)

# EXTRAPOLATE LHCIT NUMBERS FROM LAT AND LON FOR VDCS FROM HLCIT_MASTER
hl <- vector()
xc <- vector()
yc <- vector()
for (k in 1:length(vdc)){
    hl[k] <- hlcit$hlcit_code[which(hlcit$vdc_name==vdc[k])[1]]
    xc[k] <- hlcit$lat[which(hlcit$vdc_name==vdc[k])[1]]
    yc[k] <- hlcit$lon[which(hlcit$vdc_name==vdc[k])[1]]
    }
koords<-cbind(xc,yc)



# DEFINE THE AGENCY-VDC RELIEF AID NETWORK ADJACENCY MATRIX
aid_m <- matrix(0,nrow=length(all),ncol=length(all))
for (i in 1:length(ag)){
  for (j in 1:length(vd)){
    aid_m[[i,length(ag)+j]] <- 
      dim(aid_data[aid_data$impl_ag==ag[i] & aid_data$vdc==vd[j],c(3,6)])[1]
  }
}

# BUILD THE AGENCY-VDC RELIEF EFFORT AID NETWORK
av <- graph.adjacency(aid_m,mode="directed",weighted=TRUE)

# COLOR VERTICES REPRESENTING AGENCIES (GREEN) AND VDCs (BLUE) WHERE AID WAS SENT
V(av)$color<-rep("green",length(all))
for (k in 1:length(all)){
  if(is.element(all[k],vd)){
    V(av)$color[k]<-"SkyBlue2"
  }  
}

# PLOT THE AGENCY-VDC AID NETWORK
plot(av,
     layout=layout.fruchterman.reingold(av, niter=200, area=2000*vcount(av)),
     vertex.color=V(av)$color,
     vertex.size=2,
     vertex.label=NA, 
     vertex.label.color="black", 
     vertex.label.font=2, 
     vertex.label.cex=0.7, 
     edge.width=0.5*sqrt(E(av)$weight),
     edge.arrow.size=0.3,
     edge.curved=TRUE,
     edge.color=gray.colors(1))

