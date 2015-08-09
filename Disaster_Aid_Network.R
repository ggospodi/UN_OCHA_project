# Nepal Disaster Relief Distribution and Displacement Tracking Network Analysis
# author: Georgi D. Gospodinov
# date: "Augist 5, 2015"
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
# AGENCY RELIEF NETWORK AT VDC LEVEL BELOW:
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

# DEFINE THE AGENCY-VDC RELIEF AID NETWORK ADJACENCY MATRIX
aid_m <- matrix(0,nrow=length(all),ncol=length(all))
for (i in 1:length(ag)){
  for (j in 1:length(vd)){
    aid_m[[i,length(ag)+j]] <- 
      dim(aid_data[aid_data$impl_ag==ag[i] & aid_data$vdc==vd[j],c(3,5)])[1]
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


V(av)$color<-rep("green",length(all))
for (k in 1:length(all)){
  if(is.element(all[k],vd)){
    V(av)$color[k]<-"SkyBlue2"
  }  
}

for (k in 1:dim(aid_m)[1]){
  if(k-1<length(ag)){
    V(av)$size[k] <-3
    V(av)$name[k] <- ag[k]
  } else {
    V(av)$size[k] <-2
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
     edge.color = gray.colors(1))


# PLOT THE AGENCY-VDC AID NETWORK
for (k in 1:dim(aid_m)[1]){
  if(k-1<length(ag)){
    V(av)$size[k] <-3
    V(av)$name[k] <- ag[k]
  } else {
    V(av)$size[k] <-1
    V(av)$name[k] <- NA}
}

plot(av,
     layout=koords2,
     vertex.color=V(av)$color,
     vertex.size=V(av)$size,
     vertex.label=V(av)$name, 
     vertex.label.color="darkgreen", 
     vertex.label.font=2, 
     vertex.label.cex=0.75, 
     edge.width=0.05*sqrt(E(av)$weight),
     edge.arrow.size=0.2,
     edge.curved=TRUE,
     edge.color=gray.colors(1),
     main="Nepal Agency Aid Relief Geo-Network")
legend("topright",c("Implementing Aid Agency ","VDC with Geo-Coords"),fill=c("green","SkyBlue2"),bty="n")
































# EDGE-FILTRATION BY EDGE WEIGHT OF THE AGENCY-VDC AID NETWORK: CUT-OFF = 25% percentile
cut25 <- quantile(as.vector(aid_m[aid_m>0]),0.25)
av_f <- filter(cutoff = cut25,
             edge_matrix = aid_m,
             vertex_colors = V(av)$color,
             vertex_names = all)

# DISPLAY THE EDGE-FILTERED GRAPH
plot(av_f,
     layout=layout.fruchterman.reingold(av_f, niter=200, area=2000*vcount(av_f)),
     vertex.color=V(av_f)$color,
     vertex.size=4,
     vertex.label=NA, 
     vertex.label.color="black", 
     vertex.label.font=2, 
     vertex.label.cex=0.7, 
     edge.width=0.7*sqrt(E(av_f)$weight),
     edge.arrow.size=0.5,
     edge.curved=TRUE,
     edge.color=gray.colors(1))

# EDGE-FILTRATION BY EDGE WEIGHT OF THE AGENCY-VDC AID NETWORK: CUT-OFF = 50% percentile
cut50 <- quantile(as.vector(aid_m[aid_m>0]),0.5)
av_f <- filter(cutoff = cut50,
             edge_matrix = aid_m,
             vertex_colors = V(av)$color,
             vertex_names = all)

# DISPLAY THE EDGE-FILTERED GRAPH
plot(av_f,
     layout=layout.fruchterman.reingold(av_f, niter=200, area=2000*vcount(av_f)),
     vertex.color=V(av_f)$color,
     vertex.size=3,
     vertex.label=NA, 
     vertex.label.color="black", 
     vertex.label.font=2, 
     vertex.label.cex=0.7, 
     edge.width=0.3*(E(av_f)$weight),
     edge.arrow.size=0.5,
     edge.curved=TRUE,
     edge.color=gray.colors(1))

# EDGE-FILTRATION BY EDGE WEIGHT OF THE AGENCY-VDC AID NETWORK: CUT-OFF = 75% percentile
cut75 <- quantile(as.vector(aid_m[aid_m>0]),0.75)
av_f <- filter(cutoff = cut75,
             edge_matrix = aid_m,
             vertex_colors = V(av)$color,
             vertex_names = all)

# DISPLAY THE EDGE-FILTERED GRAPH
plot(av_f,
     layout=layout.fruchterman.reingold(av_f, niter=200, area=2000*vcount(av_f)),
     vertex.color=V(av_f)$color,
     vertex.size=3,
     vertex.label=NA, 
     vertex.label.color="black", 
     vertex.label.font=2, 
     vertex.label.cex=0.7, 
     edge.width=0.3*(E(av_f)$weight),
     edge.arrow.size=0.6,
     edge.curved=TRUE,
     edge.color=gray.colors(1))

# DISPLAY THE LARGEST CLUSTER (GIANT COMPONENT):
av_f_c <- giant_comp(graph = av_f,
                     vertex_colors = V(av_f)$color,
                     vertex_names = V(av_f)$name)

# PLOT THE WEIGHTED DISPLACEMENT GRAPH
plot(av_f_c,
     layout=layout.fruchterman.reingold(av_f_c, niter=200, area=2000*vcount(av_f_c)),
     vertex.color=V(av_f_c)$color,
     vertex.size=4,
     vertex.label=NA, 
     vertex.label.color="black", 
     vertex.label.font=2, 
     vertex.label.cex=1, 
     edge.width=sqrt(E(av_f_c)$weight),
     edge.arrow.size=0.6,
     edge.curved=TRUE,
     edge.color=gray.colors(1))

# FILTRATION PLUS GIANT CONNECTED COMPONENT, CUTOFF = 90% quantile
cut90 <- quantile(as.vector(aid_m[aid_m>0]),0.9)
av_f <- filter(cutoff = cut90,
               edge_matrix = aid_m,
               vertex_colors = V(av)$color,
               vertex_names = all)
av_f_c <- giant_comp(graph = av_f,
                     vertex_colors = V(av_f)$color,
                     vertex_names = V(av_f)$name)
plot(av_f_c,
     layout=layout.fruchterman.reingold(av_f_c, niter=200, area=2000*vcount(av_f_c)),
     vertex.color=V(av_f_c)$color,
     vertex.size=6,
     vertex.label=NA, 
     vertex.label.color="black", 
     vertex.label.font=2, 
     vertex.label.cex=1, 
     edge.width=sqrt(E(av_f_c)$weight),
     edge.arrow.size=0.6,
     edge.curved=TRUE,
     edge.color=gray.colors(1))

# FILTRATION PLUS GIANT CONNECTED COMPONENT, CUTOFF = 95% quantile
cut95 <- quantile(as.vector(aid_m[aid_m>0]),0.95)
av_f<-filter(cutoff = cut95,
             edge_matrix = aid_m,
             vertex_colors = V(av)$color,
             vertex_names = all)
av_f_c<-giant_comp(graph = av_f,
                   vertex_colors = V(av_f)$color,
                   vertex_names = V(av_f)$name)
plot(av_f_c,
     layout=layout.fruchterman.reingold(av_f_c, niter=200, area=2000*vcount(av_f_c)),
     vertex.color=V(av_f_c)$color,
     vertex.size=12,
     vertex.label=V(av_f_c)$name, 
     vertex.label.color="black", 
     vertex.label.font=1, 
     vertex.label.cex=1.2, 
     edge.width=0.5*(E(av_f_c)$weight),
     edge.arrow.size=0.8,
     edge.curved=TRUE,
     edge.color=gray.colors(1))
