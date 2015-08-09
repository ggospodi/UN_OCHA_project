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












# EDIT









#
#
#
#
#
#
#
# ANALYSIS OF THE AGENCY-VDC NETWORK ITSELF
#
#
#
#
#
#
#


# DEFINE THE WEIGHTED DISPLACEMENT GRAPH



# THIS IS THE NUMBER OF EDGES FROM EACH NODE
# IGNORING DIRECTION (TOTAL DEGREE)
summary(degree(av))


# THIS IS THE IN DEGREE SUMMARY
summary(degree(av,mode = "in"))


# THIS IS THE OUT DEGREE SUMMARY
summary(degree(av,mode = "out"))


# THIS IS THE WEIGHTED NUMBER OF THE ABOVE vdcENCIES, SO THE NUMBER OF SHARED VDC
# TARGETS IS ACCOUNTED FOR BETWEEN EACH PAIR OF vdcENCIES
summary(graph.strength(av))


# AGAIN, THE INWARD WEIGHTED DEGREE
summary(graph.strength(av,mode = "in"))


# AGAIN, THE OUTWARD WEIGHTED DEGREE
summary(graph.strength(av,mode = "out"))


# PLOT THE NUMBER OF DISTINCT VDC-VDC CONNECTIONS
plot(sort(degree(av)),
     col = adjustcolor(rgb(1,0,1,1)),
     pch = 19,
     xaxt = "n",
     yaxt = "n",
     ann = FALSE)
par(new=T)
plot(sort(degree(av,mode = "in")),
     col = "green",
     pch = 19,
     xaxt = "n",
     yaxt = "n",
     ann = FALSE)
par(new=T)
plot(sort(degree(av,mode = "out")),
     col = "blue",
     pch = 19,
     xlab = "VDC Index",
     ylab = "Number of VDC Transitions",
     main = "Number of VDC Transitions To and From a Given VDC (Sorted)")
legend("topleft",
       c("Total VDC-VDC Transitions","VDC In-Transitions", "VDC Out-Transitions"),
       fill = c(adjustcolor(rgb(1,0,1,1)),"green","blue"),
       bty = "n")


# THE DEGREE DISTRIBUTION (VDC-VDC CONNECTIONS)
deg1 <- hist(degree(av), breaks = 20)$counts
deg1n <- 100*deg1/sum(deg1)
deg2 <- hist(degree(av,mode = "in"), breaks = 20)$counts
deg2n <- 100*deg2/sum(deg2)
deg3 <- hist(degree(av,mode = "out"), breaks = 20)$counts
deg3n <- 100*deg3/sum(deg3)
maxn <- max(length(deg1n), length(deg2n),length(deg3n))
d1 <- append(deg1n,rep(0,maxn-length(deg1n)))
d2 <- append(deg2n,rep(0,maxn-length(deg2n)))
d3 <- append(deg3n,rep(0,maxn-length(deg3n)))
ddata <- cbind(d1,d2,d3)
barplot(t(ddata), 
        beside = T, 
        xlab = "Number of VDC-VDC Connection", 
        ylab = "Relative VDC-VDC Conneciton %", 
        main = "Distribution of VDC Transitions To and From a Given VDC",
        col = c(adjustcolor(rgb(1,0,1,1)),"green","blue"))
axis(1, 
     at = 4*(1:maxn)-2,
     labels = as.character(1:maxn),
     cex.axis = 1,
     las = 1)
legend("topright",
       legend = c("Total VDC-VDC Transitions","VDC In-Transitions", "VDC Out-Transitions"),
       col = c(adjustcolor(rgb(1,0,1,1)),"green","blue"), 
       bty = "n",pch = 15, 
       cex = 1.5)


# PLOT THE WEIGHTED NUMBER OF DISTINCT VDC-VDC CONNECTIONS
plot(sort(graph.strength(av)),
     col = adjustcolor(rgb(1,0,1,1)),
     pch = 19,
     xaxt = "n",
     yaxt = "n",
     ann = FALSE)
par(new=T)
plot(sort(graph.strength(av,mode = "in")),
     col = "green",
     pch = 19,
     xaxt = "n",
     yaxt = "n",
     ann = FALSE)
par(new=T)
plot(sort(graph.strength(av,mode = "out")),
     col = "blue",
     pch = 19,
     xlab = "VDC Index",
     ylab = "Weighted Number of VDC Transitions",
     main = "Weighted Number of VDC Transitions To and From a Given VDC (Sorted)")
legend("topleft",
       c("Total Weighted VDC-VDC Transitions","Weighted VDC In-Transitions", "Weighted VDC Out-Transitions"),
       fill = c(adjustcolor(rgb(1,0,1,1)),"green","blue"),
       bty = "n")


# THE WEIGHTED DEGREE DISTRIBUTION (VDC-VDC CONNECTIONS)
deg1 <- hist(graph.strength(av), breaks = 20)$counts
deg1n <- 100*deg1/sum(deg1)
deg2 <- hist(graph.strength(av,mode = "in"), breaks = 20)$counts
deg2n <- 100*deg2/sum(deg2)
deg3 <- hist(graph.strength(av,mode = "out"), breaks = 20)$counts
deg3n <- 100*deg3/sum(deg3)
maxn <- max(length(deg1n), length(deg2n),length(deg3n))
d1 <- append(deg1n,rep(0,maxn-length(deg1n)))
d2 <- append(deg2n,rep(0,maxn-length(deg2n)))
d3 <- append(deg3n,rep(0,maxn-length(deg3n)))
wddata <- cbind(d1,d2,d3)
barplot(t(wddata), 
        beside = T, 
        xlab = "Weighted Number of VDC-VDC Connection", 
        ylab = "Weighted Relative VDC-VDC Conneciton %", 
        main = "Distribution of Weighted VDC Transitions To and From a Given VDC",
        col = c(adjustcolor(rgb(1,0,1,1)),"green","blue"))
axis(1, 
     at = 4*(1:maxn)-2,
     labels = as.character(1:maxn),
     cex.axis = 1,
     las = 1)
legend("topright",
       legend = c("Total Weighted VDC-VDC Transitions","Weighted VDC In-Transitions", "Weighted VDC Out-Transitions"),
       col = c(adjustcolor(rgb(1,0,1,1)),"green","blue"), 
       bty = "n",pch = 15, 
       cex = 1.5)


# GRAPH DENSITY IS THE RATIO OF THE NUMBER OF EDGES AND THE NUMBER OF POSSIBLE EDGES
# TYPICALLY ON THE ORDER OF 1-10%
100*graph.density(av)


# CLUSTERS ARE CONNECTED COMPONENTS
clusters(av)$no


# SORTED CLUSTERS BY SIZE
sort(clusters(av)$csize,decreasing = TRUE)


# GLOBAL CLUSTERING COEFFICIENT (TRANSITIVITY) IS THE RATIO OF TRIANGLES AND CONNECTED TRIPLES
transitivity(av)
cut75 <- quantile(as.vector(dtm[dtm>0]),0.75)
av_f <- filter(cutoff = cut75,
               edge_matrix = dtm,
               vertex_colors = V(av)$color,
               vertex_names = V(av)$name)
av_f <- as.undirected(av_f)
transitivity(av_f)


# RELATIVE MAXIMAL CLUSTER SIZE (AS % OF NUMBER OF NODES) 
max(clusters(av)$csize)/vcount(av)


# RELATIVE NUMBER OF ISOLATED NODES (AS % OF NUMBER OF NODES)  
sum(degree(av) == 0)/vcount(av)


# PATH DISTRIBUTION: This shows the different lengths of shortest paths (geodesics) in our network. 
sh <- shortest.paths(graph = av, 
                     mode = "all",
                     weights = E(av)$weight)
is.na(sh) <- sapply(sh,is.infinite)
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
       breaks = 50,
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

av <- graph.adjacency(dtm,
                      mode = "directed",
                      weighted = TRUE)
V(av)$color <- rep("SkyBlue2",length(vdc))
for (k in 1:length(vdc)){
  o_count <- length(which(dt_data$idp_origin_vdc == vdc[k]))+
    length(which(dt_data$idp2_origin_vdc == vdc[k]))
  d_count <- length(which(dt_data$vdc == vdc[k]))
  if(o_count>d_count){
    V(av)$color[k] <- "green"
  } 
}
V(av)$name <- vdc
av <- drop_isolated(graph = av,
                    vertex_colors = V(av)$color,
                    vertex_names = V(av)$name)
av <- drop_loops(graph = av,
                 vertex_colors = V(av)$color,
                 vertex_names = V(av)$name)
ec <- evcent(av)$vector
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
av_coords <- koords[which(vdc %in% V(av)$name),]
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
     main = "Weighted Agency Aid Network Flow (VDC Level)")
legend("topleft",
       c("Highest Eigenvector Centrality","Lowest Eigenvector Centrality"),
       fill = c("red","White"),
       bty = "n")

# PLOT THE HEAT MAP OF THE AGENCY AID NETWORK
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
     main = "Weighted Agency Aid Geo-Network Flow (VDC Level)")
legend("topright",
       c("Highest Eigenvector Centrality","Lowest Eigenvector Centrality"),
       fill = c("red","White"),
       bty = "n")

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

av <- graph.adjacency(dtm,
                      mode = "directed",
                      weighted = TRUE)
V(av)$color <- rep("SkyBlue2",length(vdc))
for (k in 1:length(vdc)){
  o_count <- length(which(dt_data$idp_origin_vdc == vdc[k]))+
    length(which(dt_data$idp2_origin_vdc == vdc[k]))
  d_count <- length(which(dt_data$vdc == vdc[k]))
  if(o_count>d_count){
    V(av)$color[k] <- "green"
  } 
}
V(av)$name <- vdc
av <- drop_isolated(graph = av,
                    vertex_colors = V(av)$color,
                    vertex_names = V(av)$name)
av <- drop_loops(graph = av,
                 vertex_colors = V(av)$color,
                 vertex_names = V(av)$name)
au <- authority.score(av)$vector
au <- 100*au
plot(sort(au, decreasing = TRUE),
     col = adjustcolor(rgb(1,0,1,1)), 
     xlab = "Node Index", 
     ylab = "Eigenvector Centrality", 
     main = "Sorted VDC Network Authority Score Values", 
     pch = 19)
histP2(au,
       breaks = 50,
       col = adjustcolor(rgb(1,0,1,1)),
       xlab = "Authority Score Values",
       main = "VDC Network Authority Score Distribution")


# PLOT HEAT MAP ON VERTICES ACCORDING TO AUTHORITY SCORE
for (k in 1:length(au)){
  V(av)$color[k] <- rev(heat.colors(1+as.integer(max(au))))[as.integer(au[k])+1]
}
av_coords <- koords[which(vdc %in% V(av)$name),]
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
     main = "Weighted Agency Aid Network Flow (VDC Level)")
legend("topleft",
       c("Highest Authority Score","Lowest Authority Score"),
       fill = c("red","White"),
       bty = "n")

# PLOT THE HEAT MAP OF THE AGENCY AID NETWORK
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
     main = "Weighted Agency Aid Geo-Network Flow (VDC Level)")
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

av <- graph.adjacency(dtm,
                      mode = "directed",
                      weighted = TRUE)
V(av)$color <- rep("SkyBlue2",length(vdc))
for (k in 1:length(vdc)){
  o_count <- length(which(dt_data$idp_origin_vdc == vdc[k]))+
    length(which(dt_data$idp2_origin_vdc == vdc[k]))
  d_count <- length(which(dt_data$vdc == vdc[k]))
  if(o_count>d_count){
    V(av)$color[k] <- "green"
  } 
}
V(av)$name <- vdc
av <- drop_isolated(graph = av,
                    vertex_colors = V(av)$color,
                    vertex_names = V(av)$name)
av <- drop_loops(graph = av,
                 vertex_colors = V(av)$color,
                 vertex_names = V(av)$name)
hb <- hub.score(av)$vector
hb <- 100*hb
plot(sort(hb, decreasing = TRUE),
     col = adjustcolor(rgb(1,0,1,1)), 
     xlab = "Node Index", 
     ylab = "Eigenvector Centrality", 
     main = "Sorted VDC Network Hub Score Values", 
     pch = 19)
histP2(hb,
       breaks = 50,
       col = adjustcolor(rgb(1,0,1,1)),
       xlab = "Hub Score Values",
       main = "VDC Network Hub Score Distribution")


# PLOT HEAT MAP ON VERTICES ACCORDING TO AUTHORITY SCORE
for (k in 1:length(hb)){
  V(av)$color[k] <- rev(heat.colors(1+as.integer(max(hb))))[as.integer(hb[k])+1]
}
av_coords <- koords[which(vdc %in% V(av)$name),]
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
     main = "Weighted Agency Aid Network Flow (VDC Level)")
legend("topleft",
       c("Highest Hub Score","Lowest Hub Score"),
       fill = c("red","White"),
       bty = "n")

# PLOT THE HUB SCORES OF THE AGENCY AID NETWORK
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
     main = "Weighted Agency Aid Geo-Network Flow (VDC Level)")
legend("topright",
       c("Highest Hub Score","Lowest Hub Score"),
       fill = c("red","White"),
       bty = "n")

# PLOT HEAT MAP ON VERTICES ACCORDING TO AUTHORITY SCORE
# for (k in 1:length(hb)){
#   V(av)$color[k] <- rev(heat.colors(1+as.integer(max(hb))))[as.integer(hb[k])+1]
# }
# av_coords <- koords[which(vdc %in% V(av)$name),]
# plot(av,
#      layout = layout.fruchterman.reingold(av, 
#                                           niter = 200, 
#                                           area = 2000*vcount(av)),
#      vertex.color = V(av)$color,
#      vertex.size = 9,
#      vertex.label = V(av)$name, 
#      vertex.label.color = "black",
#      vertex.label.font = 1, 
#      vertex.label.cex = 0.85, 
#      edge.width = 0.2*sqrt(E(av)$weight),
#      edge.arrow.size = 0.5,
#      edge.curved = TRUE,
#      edge.color = gray.colors(1),
#      main = "Weighted Agency Aid Network Flow (VDC Level)")
# legend("topleft",
#        c("Highest Hub Score","Lowest Hub Score"),
#        fill = c("red","White"),
#        bty = "n")
# plot(av,
#      layout = av_coords,
#      vertex.color = V(av)$color,
#      vertex.size = 4,
#      vertex.label = NA, 
#      vertex.label.color = "black",
#      vertex.label.font = 1, 
#      vertex.label.cex = 0.75, 
#      edge.width = 0.2*sqrt(E(av)$weight),
#      edge.arrow.size = 0.5,
#      edge.curved = TRUE,
#      edge.color = gray.colors(1),
#      main = "Weighted Agency Aid Geo-Network Flow (VDC Level)")
# legend("topright",
#        c("Highest Hub Score","Lowest Hub Score"),
#        fill = c("red","White"),
#        bty = "n")
#
#
#
#
# COMMUNITY STRUCTURES FOR THE DISPLACEMENT NETWORK
#
#
#
#
#
# MULTILEVEL COMMUNITY DETECTION
#
# NOTE: THIS APPLEIS TO UNDIRECTED GRAPHS ONLY
#
# Multilevel community detection is based on the following approach. 
# Assume that we start with a weighted network of N nodes. 
# First, we assign a different community to each node of the network. 
# So, in this initial partition there are as many communities as there are nodes. 
# Then, for each node i we consider the neighbours j of i and we evaluate 
# the gain of modularity that would take place by removing i from its community 
# and by placing it in the community of j. The node i is then placed in the community 
# for which this gain is maximum (in case of a tie we use a breaking rule), 
# but only if this gain is positive. If no positive gain is possible, 
# i stays in its original community. This process is applied repeatedly 
# and sequentially for all nodes until no further improvement can be achieved. 
# Note that a node may be, and often is, considered several times. Also, 
# the output of the algorithm depends on the order in which the nodes 
# are considered, although it can be shown that the order has little effect.
#
#
av <- graph.adjacency(dtm,
                      mode = "directed",
                      weighted = TRUE)
V(av)$color <- rep("SkyBlue2",length(vdc))
for (k in 1:length(vdc)){
  o_count <- length(which(dt_data$idp_origin_vdc == vdc[k]))+
    length(which(dt_data$idp2_origin_vdc == vdc[k]))
  d_count <- length(which(dt_data$vdc == vdc[k]))
  if(o_count>d_count){
    V(av)$color[k] <- "green"
  } 
}
V(av)$name <- vdc
av <- drop_isolated(graph = av,
                    vertex_colors = V(av)$color,
                    vertex_names = V(av)$name)
av <- drop_loops(graph = av,
                 vertex_colors = V(av)$color,
                 vertex_names = V(av)$name)
av_coords <- koords[which(vdc %in% V(av)$name),]
av <- as.undirected(av)
mc <- multilevel.community(graph = av,
                           weights = E(av)$weights)
plot(mc,
     av,
     layout = layout.fruchterman.reingold(av, 
                                          niter = 200, 
                                          area = 2000*vcount(av)),
     vertex.color = mc,
     vertex.size = 9,
     vertex.label = V(av)$name, 
     vertex.label.color = "black",
     vertex.label.font = 1, 
     vertex.label.cex = 0.85, 
     edge.width = 0.3*sqrt(E(av)$weight),
     edge.arrow.size = 0.5,
     edge.curved = TRUE,
     edge.color = gray.colors(1),
     main = "Weighted Agency Aid Network Multilevel Communities")

# PLOT THE COMMUNITIES OF THE AGENCY AID NETWORK
plot(mc,
     av,
     layout = av_coords,
     vertex.color = strongclusters,
     vertex.size = 4,
     vertex.label = NA, 
     vertex.label.color = "black",
     vertex.label.font = 1, 
     vertex.label.cex = 0.75, 
     edge.width = 0.2*sqrt(E(av)$weight),
     edge.arrow.size = 0.5,
     edge.curved = TRUE,
     edge.color = gray.colors(1),
     main = "Weighted Agency Aid Geo-Network Multilevel Communities")


# SOME BASIC MULTILEVEL COMMUNITY STATS
plot(sort(mc$membership), 
     col = adjustcolor(rgb(1,0,1,1)), 
     xlab = "Node Index (VDC) in the Network", 
     ylab = "Multilevel Community Values", 
     main = "Sorted Multilevel Community Values for Nepal Displacement Network", 
     pch = 19)
histP1(mc$membership,
       breaks = 60,
       col = adjustcolor(rgb(1,0,1,1)),
       xlab = "Multilevel Community Values",
       main = "Multilevel Community Distribution for Nepal Displacement Network")

#
#
#
# WALKTRAP COMMUNITY DETECTION
#
#
# Walktrap community detection aims to find 
# densely connected subgraphs using random walks with the 
# premise that short random walks should be contained within the same cluster.
#
#
#

av <- graph.adjacency(dtm,
                      mode = "directed",
                      weighted = TRUE)
V(av)$color <- rep("SkyBlue2",length(vdc))
for (k in 1:length(vdc)){
  o_count <- length(which(dt_data$idp_origin_vdc == vdc[k]))+
    length(which(dt_data$idp2_origin_vdc == vdc[k]))
  d_count <- length(which(dt_data$vdc == vdc[k]))
  if(o_count>d_count){
    V(av)$color[k] <- "green"
  } 
}
V(av)$name <- vdc
av <- drop_isolated(graph = av,
                    vertex_colors = V(av)$color,
                    vertex_names = V(av)$name)
av <- drop_loops(graph = av,
                 vertex_colors = V(av)$color,
                 vertex_names = V(av)$name)
av_coords <- koords[which(vdc %in% V(av)$name),]
av <- as.undirected(av)
wc <- walktrap.community(graph = av,
                         weights = E(av)$weights)
plot(wc,
     av,
     layout = layout.fruchterman.reingold(av, 
                                          niter = 200, 
                                          area = 2000*vcount(av)),
     vertex.color = mc,
     vertex.size = 9,
     vertex.label = V(av)$name, 
     vertex.label.color = "black",
     vertex.label.font = 1, 
     vertex.label.cex = 0.85, 
     edge.width = 0.3*sqrt(E(av)$weight),
     edge.arrow.size = 0.5,
     edge.curved = TRUE,
     edge.color = gray.colors(1),
     main = "Weighted Agency Aid Network Walktrap Communities")

# PLOT THE COMMUNITIES OF THE AGENCY AID NETWORK
plot(wc,
     av,
     layout = av_coords,
     vertex.color = strongclusters,
     vertex.size = 4,
     vertex.label = NA, 
     vertex.label.color = "black",
     vertex.label.font = 1, 
     vertex.label.cex = 0.75, 
     edge.width = 0.2*sqrt(E(av)$weight),
     edge.arrow.size = 0.5,
     edge.curved = TRUE,
     edge.color = gray.colors(1),
     main = "Weighted Agency Aid Geo-Network Walktrap Communities")


# SOME BASIC WALKTRAP STATS
plot(sort(wc$membership), 
     col = adjustcolor(rgb(1,0,1,1)), 
     xlab = "Node Index (VDC) in the Network", 
     ylab = "Walktrap Community Values", 
     main = "Sorted Walktrap Community Values for Nepal Displacement Network", 
     pch = 19)
histP1(wc$membership,
       breaks = 60,
       col = adjustcolor(rgb(1,0,1,1)),
       xlab = "Walktrap Community Values",
       main = "Walktrap Community Distribution for Nepal Displacement Network")
