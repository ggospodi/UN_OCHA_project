# Nepal Disaster Relief Distribution and Displacement Tracking Network Analysis
# author: Georgi D. Gospodinov
# date: "August 7, 2015"
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
     layout=layout.fruchterman.reingold(av, niter=200, area=2000*vcount(av)),
     vertex.color=V(av)$color,
     vertex.size=V(av)$size,
     vertex.label=NA, 
     vertex.label.color="black", 
     vertex.label.font=2, 
     vertex.label.cex=0.2, 
     edge.width=0.2*sqrt(E(av)$weight),
     edge.arrow.size=0.2,
     edge.curved=TRUE,
     edge.color=gray.colors(1))




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














































# EDGE-FILTRATION BY EDGE WEIGHT OF THE WEIGHTED DISPLACEMENT GRAPH: CUT-OFF = 25% quantile
cut50 <- quantile(as.vector(dtm[dtm>0]),0.50)
gd <- graph.adjacency(dtm,mode="directed",weighted=TRUE)
V(gd)$name <- vdc
V(gd)$color <- rep("SkyBlue2",length(vdc))
for (k in 1:length(vdc)){
  o_count <- length(which(dt_data$idp_origin_vdc==vdc[k]))+length(which(dt_data$idp2_origin_vdc==vdc[k]))
  d_count <- length(which(dt_data$vdc==vdc[k]))
  if(o_count>d_count){
    V(gd)$color[k] <- "green"
  } 
}
gd_f <- filter(cutoff = cut50,
               edge_matrix = dtm,
               vertex_colors = V(gd)$color,
               vertex_names = V(gd)$name)
gd_f <- drop_loops(graph = gd_f,
                   vertex_colors = V(gd_f)$color,
                   vertex_names = V(gd_f)$name)

# DISPLAY THE EDGE-FILTERED GRAPH
gd_f_coords <- koords[which(vdc %in% V(gd_f)$name),]
plot(gd_f,
     layout=gd_f_coords,
     vertex.color=V(gd_f)$color,
     vertex.size=5, 
     vertex.label=V(gd_f)$name,
     vertex.label.color="black", 
     vertex.label.font=1, 
     vertex.label.cex=0.7, 
     edge.width=0.2*sqrt(E(gd_f)$weight),
     edge.arrow.size=0.7,
     edge.curved=TRUE,
     edge.color=gray.colors(1),
     main="Filtered Nepal Displacement Network Flow (VDC Level)")
legend("topright",c("Origins of Displacement",
                      "Destinations of Displacement"),fill=c("green","SkyBlue2"),bty="n")

plot(gd_f,
     layout=layout.fruchterman.reingold(gd_f, niter=200, area=2000*vcount(gd_f)),
     vertex.color=V(gd_f)$color,
     vertex.size=10,
     vertex.label=V(gd_f)$name, 
     vertex.label.color="black", 
     vertex.label.font=1, 
     vertex.label.cex=1.2, 
     edge.width=0.4*sqrt(E(gd_f)$weight),
     edge.arrow.size=01,
     edge.curved=TRUE,
     edge.color=gray.colors(1),
     main="Filtered Nepal Displacement Network Flow (VDC Level)")

# EDGE-FILTRATION BY EDGE WEIGHT OF THE WEIGHTED DISPLACEMENT GRAPH: CUT-OFF = 50% quantile
cut50 <- quantile(as.vector(dtm[dtm>0]),0.5)
gd <- graph.adjacency(dtm,mode="directed",weighted=TRUE)
V(gd)$name <- vdc
V(gd)$color <- rep("SkyBlue2",length(vdc))
for (k in 1:length(vdc)){
  o_count <- length(which(dt_data$idp_origin_vdc==vdc[k]))+length(which(dt_data$idp2_origin_vdc==vdc[k]))
  d_count <- length(which(dt_data$vdc==vdc[k]))
  if(o_count>d_count){
    V(gd)$color[k] <- "green"
  } 
}
gd_f <- filter(cutoff = cut50,
               edge_matrix = dtm,
               vertex_colors = V(gd)$color,
               vertex_names = V(gd)$name)
gd_f <- drop_loops(graph = gd_f,
                   vertex_colors = V(gd_f)$color,
                   vertex_names = V(gd_f)$name)

# DISPLAY THE EDGE-FILTERED GRAPH
plot(gd_f,
     layout=layout.fruchterman.reingold(gd_f, niter=200, area=2000*vcount(gd_f)),
     vertex.color=V(gd_f)$color,
     vertex.size=10,
     vertex.label=V(gd_f)$name, 
     vertex.label.color="black", 
     vertex.label.font=1, 
     vertex.label.cex=1.4, 
     edge.width=0.4*sqrt(E(gd_f)$weight),
     edge.arrow.size=01,
     edge.curved=TRUE,
     edge.color=gray.colors(1))

# EDGE-FILTRATION BY EDGE WEIGHT OF THE WEIGHTED DISPLACEMENT GRAPH: CUT-OFF = 75% quantile
cut75 <- quantile(as.vector(dtm[dtm>0]),0.75)
gd <- graph.adjacency(dtm,mode="directed",weighted=TRUE)
V(gd)$name <- vdc
V(gd)$color <- rep("SkyBlue2",length(vdc))
for (k in 1:length(vdc)){
  o_count <- length(which(dt_data$idp_origin_vdc==vdc[k]))+length(which(dt_data$idp2_origin_vdc==vdc[k]))
  d_count <- length(which(dt_data$vdc==vdc[k]))
  if(o_count>d_count){
    V(gd)$color[k] <- "green"
  } 
}
gd_f <- filter(cutoff = cut75,
               edge_matrix = dtm,
               vertex_colors = V(gd)$color,
               vertex_names = V(gd)$name)
gd_f <- drop_loops(graph = gd_f,
                   vertex_colors = V(gd_f)$color,
                   vertex_names = V(gd_f)$name)

# DISPLAY THE EDGE-FILTERED GRAPH
plot(gd_f,
     layout=layout.fruchterman.reingold(gd_f, niter=200, area=2000*vcount(gd_f)),
     vertex.color=V(gd_f)$color,
     vertex.size=14,
     vertex.label=V(gd_f)$name, 
     vertex.label.color="black", 
     vertex.label.font=1, 
     vertex.label.cex=1.4, 
     edge.width=0.5*sqrt(E(gd_f)$weight),
     edge.arrow.size=1.4,
     edge.curved=TRUE,edge.color=gray.colors(1))






# ANALYSIS OF THE vdcENCY NETWORK ITSELF: OVERLAP OF vdcENCY EFFORTS


# DEFINE THE WEIGHTED DISPLACEMENT GRAPH
gd <- graph.adjacency(dtm,mode="directed",weighted=TRUE)

# SET VERTEX COLORS
V(gd)$color <- rep("SkyBlue2",length(vdc))
for (k in 1:length(vdc)){
  o_count <- length(which(dt_data$idp_origin_vdc==vdc[k]))+
    length(which(dt_data$idp2_origin_vdc==vdc[k]))
  d_count <- length(which(dt_data$vdc==vdc[k]))
  if(o_count>d_count){
    V(gd)$color[k] <- "green"
  } 
}

# SET THE VERTEX LABELS
V(gd)$name <- vdc

# PLOT THE WEIGHTED DISPLACEMENT GRAPH
gd <- drop_isolated(graph = gd,
                    vertex_colors = V(gd)$color,
                    vertex_names = V(gd)$name)


# DROP LOOPS ONLY VERTICES AS WELL 
gd <- drop_loops(graph = gd,
                 vertex_colors = V(gd)$color,
                 vertex_names = V(gd)$name)

# SELECT THE COORDINATES
gd_coords <- koords[which(vdc %in% V(gd)$name),]
plot(gd,
     layout=gd_coords,
     vertex.color=V(gd)$color,
     vertex.size=4, 
     vertex.label=V(gd)$name,
     vertex.label.color="black", 
     vertex.label.font=1, 
     vertex.label.cex=0.5, 
     edge.width=0.15*sqrt(E(gd)$weight),
     edge.arrow.size=0.5,
     edge.curved=TRUE,
     edge.color=gray.colors(1),
     main="Weighted Nepal Displacement Network Flow (VDC Level)")
legend("top",c("Displacement Origin","Displacement Destination"),fill=c("green","SkyBlue2"),bty="n")


# THIS IS THE NUMBER OF vdcENCIES WITH COMMON VDC TARGETS AS A GIVEN vdcENCY
summary(degree(gd))

# THIS IS THE WEIGHTED NUMBER OF THE ABOVE vdcENCIES, SO THE NUMBER OF SHARED VDC
# TARGETS IS ACCOUNTED FOR BETWEEN EACH PAIR OF vdcENCIES
summary(graph.strength(gd))

# PLOT THE NUMBER OF DISTINCT vdcENCIES THAT SHARE TARGETS WITH A GIVEN vdcENCY
plot(sort(degree(gd)),
     col = adjustcolor(rgb(1,0,1,1)),
     pch = 19,
     xlab = "vdcency index",
     ylab = "Numer of vdcencies",
     main = "Number of vdcencies with Shared Target with an vdcency (Sorted)")

histP1(degree(gd),
       breaks = 50,
       col = adjustcolor(rgb(1,0,1,1)),
       xlab = "vdcency Network Degree Values", 
       main = "vdcency Network Degree Distribution
       (Distribution of the Number of vdcencies with Common Targets as a Given vdcency)")

# PLOT THE NUMBER OF DISTINCT vdcENCIES THAT SHARE TARGETS WITH A GIVEN vdcENCY
# WEIGHTED BY THE NUMBER OF SHARED VDC DISTRICT BETWEEN EACH PAIR OF vdcENCIES
plot(sort(graph.strength(gd)),
     col = adjustcolor(rgb(1,0,1,1)),
     pch = 19,
     xlab = "vdcency index",
     ylab = "Numer of vdcencies",
     main = "Number of vdcencies with Shared Target with an vdcency (Sorted)
     (Weighted By The Number of Shared VDCs)")

histP1(graph.strength(gd), 
       breaks = 50,
       col = adjustcolor(rgb(1,0,1,1)),
       xlab = "Weighted Degree Values", 
       main = "vdcency Network Weighted Degree Distribution
       (Weighted VDC Overlap Counts Dsitribution)")


# GRAPH DENSITY IS THE RATIO OF THE NUMBER OF EDGES AND THE NUMBER OF POSSIBLE EDGES
# TYPICALLY ON THE ORDER OF 1-10%
100*graph.density(gd)

# CLUSTERS ARE CONNECTED COMPONENTS, WE HAVE 4 in the UNFILTERED vdcENCY-VDC NETWORK 
clusters(gd)$no

# SORTED CLUSTERS BY SIZE, NOTE THAT FILTRATIONS RESULT IN INCREASED NUMBER OF CLUSTERS AND A DROP IN CLUSTER SIZE
sort(clusters(gd)$csize,decreasing=TRUE)

# GLOBAL CLUSTERING COEFFICIENT (TRANSITIVITY) IS THE RATIO FO TRIANGLES AND CONNECTED TRIPLES
transitivity(gd)
cut75 <- quantile(as.vector(vdc_m[vdc_m>0]),0.75)
gd_f <- filter(cutoff = cut75,
               edge_matrix = vdc_m,
               vertex_colors = V(gd)$color,
               vertex_names = vdc)
gd_f <- as.undirected(gd_f)
transitivity(gd_f)

# RELATIVE MAXIMAL CLUSTER SIZE (AS % OF NUMBER OF NODES) 
max(clusters(gd)$csize)/vcount(gd)

# RELATIVE NUMBER OF ISOLATED NODES (AS % OF NUMBER OF NODES)  
sum(degree(gd)==0)/vcount(gd)

# PATH DISTRIBUTION: This shows the different lengths of shortest paths (geodesics) in our network. 
sh<-shortest.paths(gd)
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
bc <- betweenness(gd,v=V(gd), directed=FALSE)
plot(sort(bc, decreasing=TRUE),
     col=adjustcolor(rgb(1/2,0,0,1/2)), 
     xlab="Node Index", 
     ylab="Betweenness Centrality", 
     main="Sorted Relief vdcency Betweenness Centrality Values", 
     pch=19)
histP2(bc,
       breaks=100,
       col=adjustcolor(rgb(1/2,0,0,1/2)),
       xlab="Betweenness Centrality Values",
       main="Relief vdcency Betweenness Centrality Distribution")

# PLOT HEAT MAP ON VERTICES ACCORDING TO BETWEENNESS CENTRALITY
bc_int <- as.integer(round(bc,0))/5
for (k in 1:length(bc_int)){
  V(gd)$color[k] <- rev(heat.colors(1+as.integer(max(bc_int))))[as.integer(bc_int[k])+1]
}

# PLOT vdcENCY GRAPH AND FILTER
plot(gd,
     layout = layout.fruchterman.reingold(gd, niter=2000, area=20000*vcount(gd)),
     vertex.color = V(gd)$color,
     vertex.size = 9,
     vertex.label = V(gd)$name, 
     vertex.label.color = "black",
     vertex.label.font = 1, 
     vertex.label.cex = 0.75, 
     edge.width=0.2*sqrt(E(gd)$weight),
     edge.arrow.size=0.5,
     edge.curved = TRUE,
     edge.color = gray.colors(1))
plot(gd,
     layout = gd_coords,
     vertex.color = V(gd)$color,
     vertex.size = 4,
     vertex.label = NA, 
     vertex.label.color = "black",
     vertex.label.font = 1, 
     vertex.label.cex = 0.75, 
     edge.width = 0.2*sqrt(E(gd)$weight),
     edge.arrow.size=0.5,
     edge.curved = TRUE,
     edge.color = gray.colors(1),
     main="Weighted Nepal Displacement Network Flow (VDC Level)")
legend("top",c("Highest BC","Lowest BC"),fill=c("red","White"),bty="n")



# FIND THE TOP 10% BETWEENNES NODES
top_bc <- bc[which(bc > quantile(bc,0.9))]
top_bc

# FIND THE TOP 5% BETWENNES NODES
top_bc <- bc[which(bc > quantile(bc,0.95))]
top_bc

# REMOVE ISOLATED
gd <- drop_isolated(graph = gd,
                    vertex_colors = V(gd)$color,
                    vertex_names = V(gd)$name)

# GET THE GIANT CONNECTED COMPONENT (TWO CLSUTERS ONLY)
gd <- giant_comp(graph = gd,
                 vertex_colors = V(gd)$color,
                 vertex_names = V(gd)$name)

# SET THE GRAPH COLOR ACCORDING TO BC

bc<-betweenness(gd,v=V(gd), directed=FALSE)
bc_int <- as.integer(round(bc,0))
for (k in 1:length(bc_int)){
  V(gd)$color[k] <- rev(heat.colors(1+as.integer(max(bc_int))))[as.integer(bc_int[k])+1]
}
plot(gd,
     layout = layout.fruchterman.reingold(gd, niter=200, area=2000*vcount(gd)),
     vertex.color = V(gd)$color,
     vertex.size = 5,
     vertex.label = NA, 
     vertex.label.color = "black",
     vertex.label.font = 1, 
     vertex.label.cex = 0.5, 
     edge.width = 0.1*E(gd)$weight,
     edge.curved = TRUE,
     edge.color = gray.colors(1))






# FILTER AND REPEAT:
# DEFINE THE WEIGHTED DISPLACEMENT GRAPH
gd <- graph.adjacency(dtm,mode="directed",weighted=TRUE)

# SET VERTEX COLORS
V(gd)$color <- rep("SkyBlue2",length(vdc))
for (k in 1:length(vdc)){
  o_count <- length(which(dt_data$idp_origin_vdc==vdc[k]))+
    length(which(dt_data$idp2_origin_vdc==vdc[k]))
  d_count <- length(which(dt_data$vdc==vdc[k]))
  if(o_count>d_count){
    V(gd)$color[k] <- "green"
  } 
}

# SET THE VERTEX LABELS
V(gd)$name <- vdc


cut25 <- quantile(as.vector(dtm[dtm>0]),0.25)
gd_f <- filter(cutoff = cut25,
               edge_matrix = dtm,
               vertex_colors = V(gd)$color,
               vertex_names = V(gd)$name)
# PLOT THE WEIGHTED DISPLACEMENT GRAPH
gd_f <- drop_isolated(graph = gd_f,
                    vertex_colors = V(gd_f)$color,
                    vertex_names = V(gd_f)$name)


# DROP LOOPS ONLY VERTICES AS WELL 
gd_f <- drop_loops(graph = gd_f,
                 vertex_colors = V(gd_f)$color,
                 vertex_names = V(gd_f)$name)
bc<-betweenness(gd_f,v=V(gd_f), directed=FALSE)
bc_int <- as.integer(round(bc,0))
for (k in 1:length(bc_int)){
  V(gd_f)$color[k] <- rev(heat.colors(1+as.integer(max(bc_int))))[as.integer(bc_int[k])+1]
}
plot(gd_f,
     layout = layout.fruchterman.reingold(gd_f, niter=200, area=2000*vcount(gd_f)),
     vertex.color = V(gd_f)$color,
     vertex.size = 5,
     vertex.label = NA, 
     vertex.label.color = "black",
     vertex.label.font = 1, 
     vertex.label.cex = 0.5, 
     edge.width = 0.1*sqrt(E(gd_f)$weight),
     edge.curved = TRUE,
     edge.color = gray.colors(1))














# EDGE-BETWEENNESS CENTRALITY: THE NUMBER OF GEODESICS GOING THROUGH AN EDGE

# DEFINE THE WEIGHTED DISPLACEMENT GRAPH
gd <- graph.adjacency(dtm,mode="directed",weighted=TRUE)

# SET VERTEX COLORS
V(gd)$color <- rep("SkyBlue2",length(vdc))
for (k in 1:length(vdc)){
  o_count <- length(which(dt_data$idp_origin_vdc==vdc[k]))+
    length(which(dt_data$idp2_origin_vdc==vdc[k]))
  d_count <- length(which(dt_data$vdc==vdc[k]))
  if(o_count>d_count){
    V(gd)$color[k] <- "green"
  } 
}

# SET THE VERTEX LABELS
V(gd)$name <- vdc

# PLOT THE WEIGHTED DISPLACEMENT GRAPH
gd <- drop_isolated(graph = gd,
                    vertex_colors = V(gd)$color,
                    vertex_names = V(gd)$name)


# DROP LOOPS ONLY VERTICES AS WELL 
gd <- drop_loops(graph = gd,
                 vertex_colors = V(gd)$color,
                 vertex_names = V(gd)$name)


ec <- edge.betweenness(graph = gd,
                       e = E(gd), 
                       directed = FALSE,
                       weights = E(gd)$weight)
plot(sort(ec, decreasing=TRUE),
     col=adjustcolor(rgb(1/2,0,0,1/2)), 
     xlab="Node Index", 
     ylab="Edge-Betweenness Centrality", 
     main="Sorted Relief vdcency Edge-Betweenness Centrality Values", 
     pch=19)
histP2(ec,
       breaks=100,
       col=adjustcolor(rgb(1/2,0,0,1/2)),
       xlab="Edge-Betweenness Centrality Values",
       main="Relief vdcency Edge-Betweenness Centrality Distribution")

# PLOT HEAT MAP ON VERTICES ACCORDING TO EDGE-BEWEENNESS CENTRALITY
ec_int <- 0.25*as.integer(round(ec,0))
for (k in 1:length(ec_int)){
  E(gd)$color[k] <- rev(heat.colors(1+as.integer(max(ec_int))))[as.integer(ec_int[k])+1]
}

# PLOT vdcENCY GRAPH AND FILTER
plot(gd,
     layout = layout.fruchterman.reingold(gd, niter=2000, area=20000*vcount(gd)),
     vertex.color = V(gd)$color,
     vertex.size=9, 
     vertex.label=V(gd)$name,
     vertex.label.color="black", 
     vertex.label.font=1, 
     vertex.label.cex=0.85, 
     edge.width=0.35*sqrt(E(gd)$weight),
     edge.arrow.size=0.5,
     edge.curved=TRUE,
     edge.color = E(gd)$color)


gd_coords <- koords[which(vdc %in% V(gd)$name),]
plot(gd,
     layout=gd_coords,
     vertex.color=V(gd)$color,
     vertex.size=4, 
     vertex.label=NA,
     vertex.label.color="black", 
     vertex.label.font=1, 
     vertex.label.cex=0.5, 
     edge.width=0.25*sqrt(E(gd)$weight),
     edge.arrow.size=0.5,
     edge.curved=TRUE,
     edge.color=E(gd)$color,
     main="Weighted Nepal Displacement Network Flow (VDC Level)")
legend("top",c("Highest EBC","Lowest EBC"),fill=c("red","white"),bty="n")





# FIND THE TOP 10% EDGE-BETWEENNES EDGES
top_ec <- ec[which(ec > quantile(ec,0.9))]
E(gd)[which(ec %in% top_ec)]

# FIND THE TOP 5% EDGE-BETWEENNES EDGES
top_ec <- ec[which(ec > quantile(ec,0.99))]
E(gd)[which(ec %in% top_ec)]


# EDGE-BETWEENNESS COMMUNITY
# The edge betweenness score of an edge measures the number of shortest paths through it, see edge.betweenness for details. 
# The idea of the edge betweenness based community structure detection is that it is likely that edges connecting separate modules 
# have high edge betweenness as all the shortest paths from one module to another must traverse through them. 
# So if we gradually remove the edge with the highest edge betweenness score we will get a hierarchical map, a rooted tree, 
# called a dendrogram of the graph. The leafs of the tree are the individual vertices and the root of the tree represents the whole graph.
#
# edge.betweenness.community performs this algorithm by calculating the edge betweenness of the graph, 
# removing the edge with the highest edge betweenness score, then recalculating edge betweenness of the edges 
# and vdcain removing the one with the highest score, etc.

# FILTER AND CLUSTER WITH CUTOFF = 0.75
gd <- graph.adjacency(dtm,mode="directed",weighted=TRUE)

# SET VERTEX COLORS
V(gd)$color <- rep("SkyBlue2",length(vdc))
for (k in 1:length(vdc)){
  o_count <- length(which(dt_data$idp_origin_vdc==vdc[k]))+
    length(which(dt_data$idp2_origin_vdc==vdc[k]))
  d_count <- length(which(dt_data$vdc==vdc[k]))
  if(o_count>d_count){
    V(gd)$color[k] <- "green"
  } 
}

# SET THE VERTEX LABELS
V(gd)$name <- vdc

# PLOT THE WEIGHTED DISPLACEMENT GRAPH
gd <- drop_isolated(graph = gd,
                    vertex_colors = V(gd)$color,
                    vertex_names = V(gd)$name)


# DROP LOOPS ONLY VERTICES AS WELL 
gd <- drop_loops(graph = gd,
                 vertex_colors = V(gd)$color,
                 vertex_names = V(gd)$name)
ebc <- edge.betweenness.community(graph = gd)
plot(ebc,
     gd, 
     layout=gd_coords,
     vertex.size=4, 
     vertex.label=NA,
     vertex.label.color="black", 
     vertex.label.font=1, 
     vertex.label.cex=0.85, 
     edge.width=0.1*sqrt(E(gd)$weight),
     edge.arrow.size=0.5,
     edge.curved=TRUE,
main="Weighted Nepal Displacement Network Communities (EB)")







gd_coords <- coords[which(vdc %in% V(gd)$name),]
plot(gd,
     layout=gd_coords,
     vertex.color=V(gd)$color,
     vertex.size=4, 
     vertex.label=NA,
     vertex.label.color="black", 
     vertex.label.font=1, 
     vertex.label.cex=0.5, 
     edge.width=0.25*sqrt(E(gd)$weight),
     edge.arrow.size=0.5,
     edge.curved=TRUE,
     edge.color=E(gd)$color,
     main="Weighted Nepal Displacement Network Flow (VDC Level)")
legend("top",c("Highest EBC","Lowest EBC"),fill=c("red","white"),bty="n")







# CLOSENESS CENTRALITY: This measure takes into account the distribution of distances to other nodes from a given node. 
# It is defined as the reciprocal of the farness of a node, where farness is defined as the sum of its distances to all 
# other nodes. Closeness can be regarded as a measure of how long it will take to spread information 
# (or an efect of an event) from a node to all other nodes. To demonstrate this concept, we compute the closeness 
# centrality for the unweighted network.
gd <- as.undirected(graph.adjacency(vdc_m,weighted=TRUE))
V(gd)$color <- rep("green",length(vdc))
V(gd)$name <- vdc
cl<-clusters(gd)
gd1<-induced.subgraph(gd, which(cl$membership == which.max(cl$csize)))
cc<-closeness(graph = gd1,vids = V(gd1),weights = E(gd1)$weight)
plot(sort(cc/max(cc), decreasing=TRUE), col=adjustcolor(rgb(1/2,0,1,1)), xlab="Node Id in the Giant Conencted Component (gg1)", ylab="Normalized Closeness Centrality", main="Closeness Centrality for the Giant Component (gg1)")
hist(cc/max(cc),breaks=200,col=adjustcolor(rgb(1/2,0,1,1)),xlab="Normalized Closeness Centrality Values for gg1",main="Normalized Closeness Centrality Distribution for gg1")

cce<-closeness.estimate(graph = gd1,vids = V(gd1),weights = E(gd1)$weight,cutoff = 7)
plot(sort(cce/max(cce), decreasing=TRUE), col=adjustcolor(rgb(1/2,0,1,1)), xlab="Node Id in the Giant Conencted Component (gg1)", ylab="Normalized Closeness Centrality", main="Closeness Centrality for the Giant Component (gg1)")
hist(cce/max(cce),breaks=200,col=adjustcolor(rgb(1/2,0,1,1)),xlab="Normalized Closeness Centrality Values for gg1",main="Normalized Closeness Centrality Distribution for gg1")





# EIGENVECTOR CENTRALITY: This is a measure of the influence of a node in the network. 
# It assigns relative scores to all nodes in the network based on the concept that connections to high-scoring nodes 
# contribute more to the score of the given node than equal conenctions to low-scoring nodes. 
# A variant of egenvector centrality is Google's PvdceRank algorithm.
gd <- graph.adjacency(dtm,mode="directed",weighted=TRUE)

# SET VERTEX COLORS
V(gd)$color <- rep("SkyBlue2",length(vdc))
for (k in 1:length(vdc)){
  o_count <- length(which(dt_data$idp_origin_vdc==vdc[k]))+
    length(which(dt_data$idp2_origin_vdc==vdc[k]))
  d_count <- length(which(dt_data$vdc==vdc[k]))
  if(o_count>d_count){
    V(gd)$color[k] <- "green"
  } 
}

# SET THE VERTEX LABELS
V(gd)$name <- vdc

# PLOT THE WEIGHTED DISPLACEMENT GRAPH
gd <- drop_isolated(graph = gd,
                    vertex_colors = V(gd)$color,
                    vertex_names = V(gd)$name)


# DROP LOOPS ONLY VERTICES AS WELL 
gd <- drop_loops(graph = gd,
                 vertex_colors = V(gd)$color,
                 vertex_names = V(gd)$name)
clu<-clusters(gd)
gd1<-induced.subgraph(gd, which(clu$membership == which.max(clu$csize)))
ec<-evcent(gd1)$vector
plot(sort(ec, decreasing=TRUE), col=adjustcolor(rgb(0,0,1,1/2)), xlab="Node Id in the Network (1:200)", ylab="Closeness Centrality Values", main="Essential (first 200 nodes) Closeness Centrality for g", pch=20)
hist(ec,breaks=100,col=adjustcolor(rgb(0,0,1,1/2)),xlab="Closeness Centrality Values",main="Essential Closeness Centrality Distribution")


# AUTHORITY SCORE: This is a measure for DIRECTED NETWORKS, and it measures the number of nodes that are hubs and point 
# to a given node. It is defined as the principle eigenvector values for t(A)*A, where A stands for the adjacency 
# matrix of the network. For undirected networks like ours, the adjacency matrix is symmetric, so the authority score 
# is equivalent to the hub score. In subsequent analyses, we will be looking at directed extensions of this network 
# model, so we are including these two scores in the analysis for completeness.
gd <- graph.adjacency(dtm,mode="directed",weighted=TRUE)

# SET VERTEX COLORS
V(gd)$color <- rep("SkyBlue2",length(vdc))
for (k in 1:length(vdc)){
  o_count <- length(which(dt_data$idp_origin_vdc==vdc[k]))+
    length(which(dt_data$idp2_origin_vdc==vdc[k]))
  d_count <- length(which(dt_data$vdc==vdc[k]))
  if(o_count>d_count){
    V(gd)$color[k] <- "green"
  } 
}

# SET THE VERTEX LABELS
V(gd)$name <- vdc

# PLOT THE WEIGHTED DISPLACEMENT GRAPH
gd <- drop_isolated(graph = gd,
                    vertex_colors = V(gd)$color,
                    vertex_names = V(gd)$name)


# DROP LOOPS ONLY VERTICES AS WELL 
gd <- drop_loops(graph = gd,
                 vertex_colors = V(gd)$color,
                 vertex_names = V(gd)$name)
au<-authority.score(gd)$vector
plot(sort(au, decreasing=TRUE), col=adjustcolor(rgb(0,0,1,1/2)), xlab="Node Id in the Network (1:200)", ylab="Authority Score Values", main="Essential (first 200 nodes) Authority Scores for g", pch=20)
hist(au,breaks=100,col=adjustcolor(rgb(0,0,1,1/2)),xlab="Authority Score Values",main="Essential Authority Score Distribution")



# HUB SCORE: This is a measure FOR DIRECTED NETWORKS and it measures the number of authority nodes that a given hub node points to.
gd <- graph.adjacency(dtm,mode="directed",weighted=TRUE)

# SET VERTEX COLORS
V(gd)$color <- rep("SkyBlue2",length(vdc))
for (k in 1:length(vdc)){
  o_count <- length(which(dt_data$idp_origin_vdc==vdc[k]))+
    length(which(dt_data$idp2_origin_vdc==vdc[k]))
  d_count <- length(which(dt_data$vdc==vdc[k]))
  if(o_count>d_count){
    V(gd)$color[k] <- "green"
  } 
}

# SET THE VERTEX LABELS
V(gd)$name <- vdc

# PLOT THE WEIGHTED DISPLACEMENT GRAPH
gd <- drop_isolated(graph = gd,
                    vertex_colors = V(gd)$color,
                    vertex_names = V(gd)$name)


# DROP LOOPS ONLY VERTICES AS WELL 
gd <- drop_loops(graph = gd,
                 vertex_colors = V(gd)$color,
                 vertex_names = V(gd)$name)
hb<-hub.score(gd)$vector
plot(sort(hb, decreasing=TRUE), col=adjustcolor(rgb(0,0,1,1/2)), xlab="Node Id in the Network (1:200)", ylab="Hub Score Values", main="Essential (first 200 nodes) Hub Scores for g", pch=20)
hist(hb,breaks=100,col=adjustcolor(rgb(0,0,1,1/2)),xlab="Hub Score Values",main="Essential Hub Score Distribution")


# CLUSTERING COEFFICIENTS: This is a measure of the clustering of the network,defined by the ratio of the number of closed triplets 
# and the number of connected triplets of vertices. We computed it in the previous report, but here we include the local and weighted version 
# of the clustering coefficients. Clustering is particularly relevant to social netowrks where nodes tend to create tightly knit groups charaterized 
# by a high density of ties, this likelihood is greater than the avervdce probability of an edge between two randomly selected nodes.
gd <- as.undirected(graph.adjacency(vdc_m,weighted=TRUE))
V(gd)$color <- rep("green",length(vdc))
V(gd)$name <- vdc
tr<-transitivity(gd, type="local")
plot(sort(tr), col=adjustcolor(rgb(0,0,1,1/2)), xlab="Node Id in the Network (1:3200)", ylab="Clustering Coefficient Values", main="Essential Clustering Coefficients for g", pch=20)
hist(tr,breaks=100,col=adjustcolor(rgb(0,0,1,1/2)),xlab="Clustering Coefficient Values",main="Clustering Coefficient Distribution for g")


# The weighted analogue of the clustering coefficient
trw<-transitivity(gd, type="weighted")
plot(sort(trw), col=adjustcolor(rgb(0,0,1,1/2)), xlab="Node Id in the Network (1:200)", ylab="Clustering Coefficient Values", main="Essential Clustering Coefficients for g", pch=20)
hist(trw,breaks=100,col=adjustcolor(rgb(0,0,1,1/2)),xlab="Clustering Coefficient Values",main="Clustering Coefficient Distribution for g")


# COMMUNITY STRUCTURES: This is a way of performing funcitonal clustering in complex networks. We have already looked at the connected components, 
# this is an elementary community detection based on connectivity.
strongclusters<-clusters(gd)$membership
plot(gd,vertex.color=strongclusters, layout=layout.fruchterman.reingold,vertex.size=4, edge.color="black", edge.width=E(gd)$weight,vertex.label=NA,main="Clustering for Store Network g200")

# ADD SOME FILTERING AND TRY vdcAIN

mc<-multilevel.community(gd)
plot(sort(mc$membership), col=adjustcolor(rgb(0,0,1,1/2)), xlab="Node Id in the Network", ylab="Multilevel Community Values", main="Multilevel Community Values for g", pch=20)
hist(mc$membership,breaks=100,col=adjustcolor(rgb(0,0,1,1/2)),xlab="Multilevel Community Values",main="Multilevel Community Distribution")


# Next, we show the walktrap community algorithm.
wc<-walktrap.community(gd)
plot(sort(wc$membership), col=adjustcolor(rgb(0,0,1,1/2)), xlab="Node Id in the Network", ylab="Walktrap Community Values", main="Walktrap Community Values for g", pch=20)
hist(wc$membership,breaks=100,col=adjustcolor(rgb(0,0,1,1/2)),xlab="Walktrap Community Values",main="Walktrap Community Distribution")


plot(wc,gd,vertex.size=4, vertex.label=NA,edge.width=E(gd)$weight,main="Walktrap Community Detection for g200")
plot(gd, vertex.color=membership(wc), vertex.size=6, edge.color="black", edge.width=E(gd)$weight,vertex.label=NA,main="Walktrap Community Detection for g200")



# quickyl explore severity

sev <- read.csv(paste0(DIR,"severity.csv"))
sev$vdc <- as.character(sev$vdc)
sev$vdc <- mapvalues(sev$vdc,
                     from = c("Agara","Baruneshwor","Betini","BhaktapurN.P.","BhimesworN.P.",
                              "ChandeniMandan","Chhatara","GunsiBhadaure","HetaudaN.P.","JaisithokMandan",
                              "JhangajholiRalmata","Jhyaku","JyamdiMandan","KakurThakur","KathmanduN.P.",
                              "LalitpurN.P.","Sangu","KirtipurN.P.","TokhaChandeswori","Thulogoun",
                              "Talkududechour","Sankhu","Puranagau","PokhariNarayansthan",
                              "Pukhulachhi","NaikapPuranoBhanjya","Mankha","Mahankal","MadhyapurThimiN.P.",
                              "Lamidada","Laharepouwa","Daxinkali","Bajrayogini","Budanilkantha","Fulpingkatti",
                              "Orang"),
                     to = c("Agra","Barudeshwor","Beteni","Bhaktapur Municipality","Bhimeswor Municipality",
                            "Chandeni Mandan","Chautara","Gunsi","Hetauda Municipality","Jaisithok Mandan",
                            "Jhangajholi Ratmata","Jhyanku","Jyamdi Mandan","Kakur Thakur","Kathmandu Metropolitan",
                            "Lalitpur Sub Metropolitan","Sangkhu","Kirtipur Municipality","Tokhachandeshwari",
                            "Thulo Gaun","Talkudunde Chaur","Sangkhu Suntol","Puranagaun","Pokhari Narayansthan",
                            "Pukulachhi","Naikap Naya","Mangkha","Mahangkal","Madhyapur Thimi Municipality",
                            "Lamidanda","Laharepauwa","Dakshinkali","Sangkhu Bajrayogini","Budhanilkantha",
                            "Phulpingkatti","Worang"))

sev$vdc[218] <- "Betini"
sev$vdc[624] <-"Lamidada"


bc
































































# AGENCY RELIEF NETWORK AT VDC LEVEL BELOW:




# CHANGE FOMRAT TO CHARACTER FOR VDC AND AGENCY NAMES
aid_data$vdc <- trim(as.character(aid_data$vdc))
aid_data$impl_agency <- trim(as.character(aid_data$impl_agency))

# FILTER OUT THE EMPTY ENTRIES
aid_data <- aid_data[nchar(aid_data$vdc)>0 & nchar(aid_data$impl_agency)>0,]
aid_data <- rm_space(aid_data,"vdc_code")
aid_data$vdc_code <- as.numeric(levels(aid_data$vdc_code))[aid_data$vdc_code]
for (k in 1:dim(aid_data)[1]){
  aid_data$vdc[k] <- hlcit$vdc_name[which(hlcit$hlcit_code %in% aid_data$vdc_code[k])[1]]
}

# SELECT UNIQUE AGENCIES AND TARGET VDC
ag <- unique(aid_data$impl_agency)
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
      dim(aid_data[aid_data$impl_agency==ag[i] & aid_data$vdc==vd[j],c(3,6)])[1]
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
