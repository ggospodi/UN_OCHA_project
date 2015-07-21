# Nepal Disaster Relief Distribution and Displacement Tracking Network Analysis
# author: Georgi D. Gospodinov
# date: "July 11, 2015"
# 
# Data Sources:
#
# Tables: CCCM_Nepal_DTM_R2.csv
# agency_relief.csv
# centroids.csv
# 
# This report contains the initial displacement tracking network model construction and some analytics.
# 
# 1. NEPAL DISPLACEMENT TRACKING NETWORK CONSTRUCTION AND ANALYSIS
# 2. Nepal Disaster Relief Distribution Network Construction and Analysis
# 3. Nepal Disaster Agency Network Construction and Analysis
# 4. Severity Index Correlation With Disaster Agency NEtwork. 
# 
# 
# Definition of the Nepal Displacement Tracking Network: 
#   
# NOTE: This report is intended to only demonstrate the construction of the networks and some of the analytical tools. 
# In subsequent reports, we will develop the analytics further and address the actionable advances that this apporach offers.
#
#
#
# LOAD PACKAGES
library(igraph)
library(RColorBrewer)

# SET FILE SOURCE PATH
DIR <- "/Users/ggospodinov/Desktop/UN_OCHA_project/data/"

# DEFINE FUNCTIONS


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




# SECTION 1: NEPAL DISPLACEMENT TRACKING NETWORK CONSTRUCTION AND ANALYSIS
#
#
# COLUMN NAMES FOR CCCM_Nepal_DTM_R2.csv ARE:
#
# 1.1.c.1 Site ID (SSID)
# 1.1.d.1 Site Name
# 1.1.a.2 Survey Round
# 1.1.a.1 Survey date (DD.MM.YYYY)
# 1.1.e.2 Zone
# 1.1.e.3 District
# 1.1.e.4 VDC
# 1.1.f.2 Site GPS Latitude
# 1.1.f.1 Site GPS Longitude
# 1.4.a.2 Site open?  
# 1.4.c.1 Closing Date (DD.MM.YYYY)  
# 1.3.a.1 Site Classification  
# 1.3.b.1 Site Type  
# If other, please specify  
# 1.4.a.1 Site start/open date (DD.MM.YYYY)  
# 1.4.b.1 Site Expected Closing Date  
# 1.1.i.1 Accessibility to site  
# 1.3.d.1 Ownership of land of site  
# 1.3.c.1 What is the most common type of shelter  
# If other or No Answer, please specify  
# 1.2.b.1 Is there a Site Management Committee (SMC) at the site?  
# 1.2.b.10 % of women participating in the Site Management Committee (SMC)  
# 1.2.b.2 Is the site management committee (SMC) made up from the community at the site?  
# 1.2.b.7 SMC member name or focal point at site  
# 1.2.b.8 SMC member phone number  
# 1.2.c.1 Is there a Site Management Agency (SMA) at the site?  
# 1.2.c.2 Type of SMA  
# If other, Specify.  
# 1.2.c.3 Name of SMA  
# 1.2.c.4 SMA member name or focal point at site  
# 1.2.c.5 SMA member phone number  
# 1.2.a.1 Is there any registration activity/IDP List Maintained  
# 1.2.n.1 Is WASH support being provided at the site ?  
# If yes who provides the service?  
# 1.2.o.1 Is HEALTH support being provided at the site  
# If yes who provides the service?  
# 1.2.p.1 Is SHELTER/NFI support provided at the site ?  
# If yes who provides the service?  
# 1.2.q.1 Is FOOD support provided at the site ?  
# If yes who provides the service?  
# 1.2.r.1 Is PROTECTION support provided at the site ?  
# If yes who provides the service?  
# 1.2.s.1 Is EDUCATION support provided at the site ?  
# If yes who provides the service?  
# 1.2.t.1 Is LIVELIHOOD support being provided at the site?  
# If yes who provides the service?  
# 1.5.b.2 Place of Origin of the largest IDP group in camp (Zone)  
# 1.5.b.2 Place of Origin of the largest IDP group in camp (District)  
# 1.5.b.2 Place of Origin of the largest IDP group in camp (VDC)  
# 1.5.b.2 Place of Origin of the largest IDP group in camp (Ward)  
# 1.5.c.2 Place of Origin of the second largest IDP group (Zone)  
# 1.5.c.2 Place of Origin of the second largest IDP group District)  
# 1.5.c.2 Place of Origin of the second largest IDP group (VDC)  
# 1.5.c.2 Place of Origin of the second largest IDP group (Ward)  
# 2.1.a.1 Total Number of IDP Families/HHs  
# 2.1.c.1 Number of Males by age: <1 year  
# 2.1.c.1 Number of Males by age: 1-5 year  
# 2.1.c.1 Number of Males by age: 6-17 year  
# 2.1.c.1 Number of Males by age: 18-59 year  
# 2.1.c.1 Number of Males by age: 60+ year  
# 2.1.b.2 Total number of IDP Male Individuals  
# 2.1.d.1 Number of Females by age: <1 year  
# 2.1.d.1 Number of Females by age: 1-5 year  
# 2.1.d.1 Number of Females by age: 6-17 year  
# 2.1.d.1 Number of Females by age: 18-59 year  
# 2.1.d.1 Number of Females by age: 60+ year  
# 2.1.b.2 Total number of IDP Female Individuals  
# 2.1.b.2 Total number of IDP  
# 2.2.c.1 Number of pregnant women  
# 2.2.d.1 Number of breastfeeding mothers  
# 2.2.g.1 Number of persons with Disabilities  
# 2.2.f.1 Number of Persons with Chronic Diseases or Serious Medical Conditions  
# 2.2.n.1 Number of single female-headed households  
# 2.2.p.1 Number of child-headed households  
# 2.2.x.3 Number of elderly-headed households  
# 2.2.u.1 Number of members of marginalized caste/ethnicity  
# 2.3.c.1 Date of arrival of first IDP group  
# 2.3.c.2 Date of arrival of last IDP group  
# 2.3.e.1 Area of intended return for largest IDP group  
# 2.3.b.4 Have IDPs previously been displaced  
# 2.3.e.7 What is preventing the largest IDP group from returning home?  
# If other, please specify  
# 2.3.g.1 Estimated % of IDPs sleeping in the site  
# 2.3.f.1 Is there relocation plan for the IDPs?  
# 3.1.a.1 Percentage of HH living outside (no shelter)  
# 3.1.b.1 Percentage of HH living in tents  
# 3.1.c.1 Percentage of HH living in makeshift shelters  
# 3.1.d.1 Percentage of HH living indoors (solid walls)  
# 3.2.a.1 Percentage of HH with access to electricity  
# 3.2.b.1 Percentage of HH with access to safe cooking facilities  
# 3.2.c.1 Percentage of HH with private living area  
# 3.7.c.2 Percentage of HH with mosquito nets  
# 3.7.j.1 Most needed type of NFI  
# If other or No Answer, please specify  
# 3.7.j.2 Second most needed type of NFI  
# If other or No Answer, please specify  
# 3.7.j.3 Third most needed type of NFI  
# If other or No Answer, please specify  
# 3.8.b.1 Is there a need for shelter repair materials  
# 3.8.b.1 Is there a need for tools for shelter construction?  
# 4.1.a.1 Location of site's main water source (walking, one-way)  
# 4.1.b.1 Main non-drinking water source provided or available  
# If other, please specify  
# 4.1.f.1 What is the main drinking water source provided or available  
# If other, please specify  
# 4.2.a.1 Average amount of water available per day and per person  
# 4.3.e.1 What is the main problem with the water?  
# If other, please specify  
# 4.4.j.1 Condition of most of the latrines  
# 4.4.a.1 Number of functioning toilets available on-site  
# 4.8.b.1 Do toilets and bathrooms have locks from the inside  
# 4.1.g.1 Is the drinking water potable ?  
# 4.4.b.1 Availability of separate male and female toilets  
# 4.5.b.1 Availability of separate male and female bathing areas  
# 4.6.a.1 Garbage Disposal Type  
# If other, please specify  
# 4.4.n.1 Evidence of open defecation?  
# 4.7.g.1 Availability of hand-washing station filled in with water and soap close to the toilets?  
# 4.7.g.2 Evidence of hand-washing practices?  
# 5.1.a.1 Is there access to food (distribution, vouchers, trade, fishing…)  
# 5.1.e.1 Is there access to a market near from the site?  
# 5.1.d.1 Most common source for obtaining food  
# If other, please specify  
# 6.1.a.1 Screening for malnutrition conducted in the area? (i.e. weight, height, mid-upper arm circumference screening)  
# 6.1.b.1 Availability of supplementary feeding for pregnant and lactating mothers  
# 6.1.c.1 Availability of supplementary feeding for children  
# 7.1.b.1 What is the most prevalent health problem at the site  
# If other, please specify  
# 7.1.b.2 What is the second most prevalent health problem at the site  
# If other, please specify  
# 7.1.b.3 What is the third most prevalent health problem at the site  
# If other, please specify  
# 7.2.a.1 Access to health facilities  
# 7.2.a.4 Do most women utilize health facilities?  
# 7.2.b.1 Location of health facilities/services  
# 7.2.c.2 Who provides health facilities/services  
# If other, please specify  
# 7.2.a.3 Regular access to medicine  
# 8.1.b.1 Access to formal/informal education services for children from displaced HHs  
# 8.2.b.1 Distance to nearest education facility  
# 8.3.a.2 Of children in site attending school, what % are girls?  
# 8.3.a.2 Of children in site attending school, what % are boys?  
# 8.2.a.1 Location of formal/informal education facilities/services for children from displaced HHs  
# 9.1.a.1 Occupation/trade of majority of displaced households (coping mechanism)  
# If other, please specify  
# 9.2.h.1 Percentage of HH in the site with a source of income  
# 9.2.i.1 Access to income generating activities  
# 9.2.k.1 Do the majority of HHs in the site receive remittances?  
# 9.3.a.2 Is there livestock on site  
# 9.1.c.1 Do IDPs have access to land for cultivation?  
# 11.3.b.1 Do any recruiters come to the site for day labour, domestic work or working abroad?  
# 11.3.a.2 if yes, to which city/country  
# 10.1.a.2 Is there security on-site/settlement areas  
# 10.1.e.1 Are security incidents reported in the site  
# 10.1.b.1 Who provides main security in the site  
# If other, please specify  
# 10.2.g.1 Most common type of security incidents reported/known occurring in the settlement area  
# If others, please specify  
# 10.2.g.5 Are security incidents commonly reported in the community before the earthquake?  
# 10.2.i.4 Most reported problem in receiving support  
# If other, please specify  
# 10.1.c.1 Number of children friendly spaces  
# 10.1.d.1 Number of women friendly spaces  
# 10.1.f.1 Do the majority of people have identification card or other documentation?  
# 10.2.j.2 Reporting/referral mechanism for GBV survivors  
# 10.2.t.2 Are women segregated during menstruation?  
# If yes, where do they go?  
# 10.2.t.1 Are there issues of discrimination towards minority groups?  
# 10.3.a.1 Do men feel safe in site?  
# 10.3.a.2 Do children feel safe in site  
# 10.3.a.3 Do women feel safe at the site  
# 10.1.s.1 Is there adequate lighting in the majority of communal point (WASH facilities, public spaces…)  
# 11.1.a.1 Where do residents mostly get their information from?  
# If other, please specify  
# 11.2.c.1 Are complaints being reported?  
# 11.2.c.2 If yes, to whom?  
# 11.1.c.1 What is the main topic on which the community is requesting information on?  
# If other, please specify  
# 11.1.h.1 Is everyone aware that donations do not need to be exchanged for anything?  
# Site classification


# LOAD VDC CENTROIDS FOR VISUALIZATION PURPOSES
centroids <- read.csv(paste0(DIR,"centroids.csv"))


# Attempts to call the file directly from online HDX server:
# library(XLConnect)
# data1<-readWorksheetFromFile("http://data.hdx.rwlabs.org/dataset/io/CCCM Nepal Displacement Tracking Matrix.xlsx",sheet=1)
# library(xlsx)
# data1<-read.xlsx("https://www.dropbox.com/s/6powpj6wsp9r9aw/De-identified%20SPUS.xlsx", sheetIndex=1)


# READ IN THE DISPLACEMENT FILE DATA
dt_data <- read.csv(paste0(DIR,"CCCM_Nepal_DTM_R2.csv"), sep=",")


# SET VAIRABLE NAMES THAT ARE SHORTER AND DO NOT CONTAIN SPACES AND OTHER SYMBOLS
nam <- c("ssid","site_name","survey_round","survey_date","zone","district","vdc","lat","lon","is_open","closing_date",
       "site_class","site_type","if_other","site_start","expect_close","access","owner","shelter_type",
       "specify","is_smc","women_smc","smc_locals","smc_members","member_phone","is_sma","sma_type",
       "if_other2","sma_name","sma_member","sma_phone","is_idp","is_wash","wash_provider","is_health",
       "health_provider","is_nfi","nfi_provider","is_food","food_provider","is_protect","protector",
       "is_educate","educator","is_livelihood","livelihood_provider","idp_origin_zone","idp_origin_district",
       "idp_origin_vdc","idp_origin_ward","idp2_origin_zone","idp2_origin_district","idp2_origin_vdc",
       "idp2_origin_ward","idp_hh","males_age_1","males_age_1-5","males_age_6-17","males_age_18-59",
       "males_age_over60","total_idp_males","females_age_1","females_age_1-5","females_age_6-17",
       "females_age_18-59","females_age_over60","total_idp_females","total_idp","preg","breastfeed",
       "disabled","chronical","female_hh","child_hh","elder_hh","marginal","idp_arrival","idp_last_arrival",
       "return_area","prev_idp","why_idp","if_other3","idp_pct_sleep","is_relocation","hh_pct_out","hh_pct_tents",
       "hh_pct_mkshift","hh_pct_indoors","hh_pct_electricity","hh_pct_cooking","hh_pct_private","hh_pct_nets",
       "nfi_most_need","if_other4","nfi2_most_need","if_other5","nfi3_most_need","if_other6","is_repairs",
       "is_tools","water_loc","non_drinking_water","if_other7","drinking_water","if_other8","ave_water_per_day",
       "water_prob","if_other9","latrines","no_wc","wc_locks","water_potable","m_f_wc","m_f_bath","garbage",
       "if_other10","open_wc","hand_wash","is_hand_wash","food_access","market_access","food_src","if_other11",
       "screening","food_for_preg","food_child","illness","if_other12","illness2","if_other13","illness3",
       "if_other14","access_hosp","women_use_hosp","hosp_loc","who_doc","if_other15","access_meds","access_school",
       "dist_school","girls_pct_school","boys_pct_school","school_loc","coping_mech","if_other16","hh_pct_income",
       "paid_work","is_remittance","is_livestock","is_land","recruit","to_which","is_security","reported",
       "main_security","if_other17","common_crimes","if_other18","before_reported","biggest_prob","if_other19",
       "child_friendly","women_friendly","id","gbv_survivors","w_segregated","where","discrimination","men_safe",
       "children_safe","women_safe","lights","info","if_other20","complaints","to_whom","info_need","if_other21",
       "donations","site_class2")
colnames(dt_data) <- nam


# DROP THE FIRST THREE ROWS DUE TO EXTRA NOTATION
dt_data <- dt_data[4:dim(dt_data)[1],]


# FORMAT THE POPULATION COUNT (IN HH) OF DISPLACEMENT
dt_data$idp_hh <- as.numeric(levels(dt_data$idp_hh))[dt_data$idp_hh]


# COLUMNS FOR THE INITIAL MODEL:
# 1.1.c.1 Site ID (SSID)
# 1.1.d.1 Site Name
# 1.1.a.2 Survey Round
# 1.1.a.1 Survey date (DD.MM.YYYY)
# 1.1.e.2 Zone
# 1.1.e.3 District
# 1.1.e.4 VDC
# 1.1.f.2 Site GPS Latitude
# 1.1.f.1 Site GPS Longitude
# 1.5.b.2 Place of Origin of the largest IDP group in camp (Zone)  
# 1.5.b.2 Place of Origin of the largest IDP group in camp (District)  
# 1.5.b.2 Place of Origin of the largest IDP group in camp (VDC)  
# 1.5.b.2 Place of Origin of the largest IDP group in camp (Ward)  
# 1.5.c.2 Place of Origin of the second largest IDP group (Zone)  
# 1.5.c.2 Place of Origin of the second largest IDP group District)  
# 1.5.c.2 Place of Origin of the second largest IDP group (VDC)  
# 1.5.c.2 Place of Origin of the second largest IDP group (Ward)  
# 2.1.a.1 Total Number of IDP Families/HHs 


# SELECT INITIAL COLUMNS, SAVE THE FILE
names <- c("ssid","site_name","survey_round","survey_date","zone","district",
         "vdc","lat","lon","idp_origin_zone","idp_origin_district","idp_origin_vdc",
         "idp_origin_ward","idp2_origin_zone","idp2_origin_district","idp2_origin_vdc",
         "idp2_origin_ward","idp_hh")
dt_data <- dt_data[,names]
write.csv(dt_data,file=paste0(DIR,"Nepal_DTM_Network_model1.csv"))


# COLUMN NAMES FOR agency_relief.csv ARE:
aid_data <- read.csv(paste0(DIR,"agency_relief.csv"), sep=",")
colnames(aid_data) <- c("district_code","vdc_code","impl_agency","src_agency","district","vdc")


# READ IN DISPLACEMENT DATA
dt_data <- read.csv(paste0(DIR,"Nepal_DTM_Network_model1.csv"))
dt_data <- dt_data[,2:dim(dt_data)[2]]


# FORMAT THE VDC RELATED COLUMNS AS CHARACTER
dt_data$vdc <- as.character(dt_data$vdc)
dt_data$idp_origin_vdc <- as.character(dt_data$idp_origin_vdc)
dt_data$idp2_origin_vdc <- as.character(dt_data$idp2_origin_vdc)


# FILTER TO USE DESTINATION VDC LEVELS AS WELL AS ORIGIN (1 OR 2) VDC ENTRIES THAT ARE NON-EMPTY
dt_data <- dt_data[nchar(dt_data$vdc)>0 & (nchar(dt_data$idp_origin_vdc) + nchar(dt_data$idp2_origin_vdc))>0,]


# CREATE A LIST OF UNIQUE VDC DESTINATION NAMES
vdcs <- unique(dt_data$vdc)


# CREATE TWO LISTS OF VDC ORIGIN NAMES, FOR THE LARGEST (1) AND SECOND LARGEST (2) POPULATION BY SHELTER VDC
vdc1_o <- unique(dt_data$idp_origin_vdc[nchar(dt_data$idp_origin_vdc)>0])
vdc2_o <- unique(dt_data$idp2_origin_vdc[nchar(dt_data$idp2_origin_vdc)>0])


# UNION OF ALL LEVEL 1 AND 2 ORIGIN VDC NAMES
vdc_o <- union(vdc1_o,vdc2_o)


# OVERALL LIST OF VDC CODES, BOTH ORIGIN AND DESTINATION
vdc <- unique(trim(union(vdc1_o,vdcs)))


# GET COORDINATES OF THE VDC CENTROIDS FROM CHRIS AFTER CONVERTING AND TRIMMING WHITESPACES FROM NAMES
# centroids$name <- trim(as.character(centroids$name))
# xc <- vector()
# yc <- vector()
# for (k in 1:length(vdc)){
#     xc[k] <- centroids$lat[which(centroids$name==vdc[k])[1]]
#     yc[k] <- centroids$lon[which(centroids$name==vdc[k])[1]]
#     }
# coords<-cbind(xc,yc)


# VDC-LEVEL ADJACENCY MATRIX FOR THE DISPLACEMENT GRAPH, AT THE LEVEL OF DISPLACEMENT TRACK
vdc_m <- matrix(0,nrow=length(vdc),ncol=length(vdc))


# CREATE A MATRIX THAT TRACKS THE NUMBERS OF DISPALCED POPULATIONS AS WELL (THIS COULD BE USED TO DERIVE EDGE WEIGHTS)
dtm <- matrix(nrow=length(vdc),ncol=length(vdc))


# COMPUTE THE ENTRIES OF BOTH THE TRACKING DISPLACEMENT MATRIX AND THE WEIGHTED DISPLACEMENT TRACKING MATRIX
for (i in 1:length(vdc)){
  for (j in 1:length(vdc)){
    vdc_m[[i,j]]<-length(dt_data[dt_data$idp_origin_vdc==vdc[i] & dt_data$vdc==vdc[j],1])+
      length(dt_data[dt_data$idp2_origin_vdc==vdc[i] & dt_data$vdc==vdc[j],1])
    dtm[[i,j]]<-sum(dt_data[dt_data$idp_origin_vdc==vdc[i] & dt_data$vdc==vdc[j],18])+
      sum(dt_data[dt_data$idp2_origin_vdc==vdc[i] & dt_data$vdc==vdc[j],18])
    }
  }


# BUILD THE DIRECTED WEIGHTED VDC NETWORK
gv <- graph.adjacency(vdc_m,mode="directed",weighted=TRUE)


# COLOR VDC NAMES OF ORIGIN (GREEN) AND VDC NAMES OF DESTINATION (BLUE)
V(gv)$color <- rep("SkyBlue2",length(vdc))
for (k in 1:length(vdc)){
  o_count <- length(which(dt_data$idp_origin_vdc==vdc[k]))+length(which(dt_data$idp2_origin_vdc==vdc[k]))
  d_count <- length(which(dt_data$vdc==vdc[k]))
  if(o_count>d_count){
    V(gv)$color[k] <- "green"
  } 
}


# SET THE VERTEX LABELS
V(gv)$name <- vdc


# PLOT THE VDC DISPLACEMENT TRACKING GRAPH
# plot(gv,
#      layout=layout.fruchterman.reingold(gv, niter=20, area=2000*vcount(gv)),
#      vertex.color=V(gv)$color,vertex.size=9, vertex.label=V(gv)$name,
#      vertex.label.color="black", vertex.label.font=2, vertex.label.cex=0.7, 
#      edge.width=0.3*(E(gv)$weight),edge.arrow.size=0.7,edge.curved=FALSE,edge.color=gray.colors(1))

# DROP ISOLATED VERTICES (NO DISPLACEMENT)
gv <- drop_isolated(graph = gv,
                    vertex_colors = V(gv)$color,
                    vertex_names = vdc)
plot(gv,
     layout=layout.fruchterman.reingold(gv, niter=20, area=2000*vcount(gv)),
     vertex.color=V(gv)$color,
     vertex.size=9, 
     vertex.label=V(gv)$name,
     vertex.label.color="black", 
     vertex.label.font=2, 
     vertex.label.cex=0.7, 
     edge.width=0.3*(E(gv)$weight),
     edge.arrow.size=0.5,
     edge.curved=FALSE,
     edge.color=gray.colors(1))

# DROP LOOPS ONLY VERTICES AS WELL ()
gv <- drop_loops(graph = gv,
                 vertex_colors = V(gv)$color,
                 vertex_names = vdc)

# RESULTING CLEANED UP GRAPH SHOWING NONTRIVIAL MIGRATION
plot(gv,
     layout=layout.fruchterman.reingold(gv, niter=20, area=2000*vcount(gv)),
     vertex.color=V(gv)$color,
     vertex.size=9, 
     vertex.label=V(gv)$name,
     vertex.label.color="black", 
     vertex.label.font=2, 
     vertex.label.cex=0.7, 
     edge.width=2*(E(gv)$weight),
     edge.arrow.size=0.6,
     edge.curved=TRUE,
     edge.color=gray.colors(1))

# DEFINE THE WEIGHTED DISPLACEMENT GRAPH
gd <- graph.adjacency(dtm,mode="directed",weighted=TRUE)

# SET VERTEX COLORS
V(gd)$color <- rep("SkyBlue2",length(vdc))
for (k in 1:length(vdc)){
  o_count <- length(which(dt_data$idp_origin_vdc==vdc[k]))+length(which(dt_data$idp2_origin_vdc==vdc[k]))
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
plot(gd,
     layout=layout.fruchterman.reingold(gd, niter=20, area=2000*vcount(gd)),
     vertex.color=V(gd)$color,
     vertex.size=9, 
     vertex.label=V(gd)$name,
     vertex.label.color="black", 
     vertex.label.font=2, 
     vertex.label.cex=0.7, 
     edge.width=0.2*sqrt(E(gd)$weight),
     edge.arrow.size=0.7,
     edge.curved=TRUE,
     edge.color=gray.colors(1))

# DROP LOOPS ONLY VERTICES AS WELL 
gd <- drop_loops(graph = gd,
                 vertex_colors = V(gd)$color,
                 vertex_names = V(gd)$name)

# RESULTING CLEANED UP GRAPH SHOWING NONTRIVIAL MIGRATION
plot(gd, 
     layout=layout.fruchterman.reingold(gd, niter=20, area=2000*vcount(gd)),
     vertex.color=V(gd)$color,
     vertex.size=9, 
     vertex.label=V(gd)$name,
     vertex.label.color="black", 
     vertex.label.font=2, 
     vertex.label.cex=0.7, 
     edge.width=0.3*sqrt(E(gd)$weight),
     edge.arrow.size=0.7,
     edge.curved=TRUE,
     edge.color=gray.colors(1))


# DISPLAY THE LARGEST CLUSTER (GIANT COMPONENT):
gd_c <- giant_comp(graph = gd,
                   vertex_colors = V(gd)$color,
                   vertex_names = V(gd)$name)

# PLOT THE WEIGHTED DISPLACEMENT GRAPH
plot(gd_c,
     layout=layout.fruchterman.reingold(gd_c, niter=200, area=2000*vcount(gd_c)),
     vertex.color=V(gd_c)$color,
     vertex.size=12,
     vertex.label=V(gd_c)$name, 
     vertex.label.color="black",
     vertex.label.font=2, 
     vertex.label.cex=1, 
     edge.width=0.5*sqrt(E(gd_c)$weight),
     edge.arrow.size=1.0,edge.curved=TRUE,edge.color=gray.colors(1))

# EDGE-FILTRATION BY EDGE WEIGHT OF THE WEIGHTED DISPLACEMENT GRAPH: CUT-OFF = 25% quantile
cut25 <- quantile(as.vector(dtm[dtm>0]),0.25)
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
gd_f <- filter(cutoff = cut25,
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
     vertex.label.cex=1, 
     edge.width=0.3*sqrt(E(gd_f)$weight),
     edge.arrow.size=0.8,
     edge.curved=TRUE,
     edge.color=gray.colors(1))


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




# AGENCY RELIEF NETWORK AT VDC LEVEL BELOW:




# CHANGE FOMRAT TO CHARACTER FOR VDC AND AGENCY NAMES
aid_data$vdc <- trim(as.character(aid_data$vdc))
aid_data$impl_agency <- trim(as.character(aid_data$impl_agency))

# FILTER OUT THE EMPTY ENTRIES
aid_data <- aid_data[nchar(aid_data$vdc)>0 & nchar(aid_data$impl_agency)>0,]

# SELECT UNIQUE AGENCIES AND TARGET VDC
ag <- unique(aid_data$impl_agency)
vd <- unique(aid_data$vdc)
all <- union(ag,vd)

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




# PROJECT EACH GRAPH COMPONENT WITH APPROPRIATE CONNECTIONS



# NOTE FOR THE AGENCY NETWORK PROJECTION, WE ARE COUNTING THE NUMBER OF
# INSTANCES WHEN TWO AGENCIES SUPPLIED AID TO THE SAME VDC
# WE CAN ALSO AUGMENT THIS MEASURE BY ACCOUNTING FOR THE NUMBER OF DIFFERENT 
# COMMON INSTANCES OF AID WITHIN A VDC, A MORE GRANULAR APPROACH
ag_m <- matrix(0,nrow=length(ag),ncol=length(ag))
for (i in 1:length(ag)){
  for (j in 1:length(ag)){
    common <- aid_m[c(i,j),(length(ag)+1):dim(aid_m)[1]]
    common[common>0] <-1
    ag_m[[i,j]]<-sum(t(common[1,])*t(common[2,]))
  }
}

# REMOVE SELF LOOPS
for (k in 1:dim(ag_m)[1]){ag_m[[k,k]] <- 0}

# DEFINE AGENCY GRAPH
agg <- as.undirected(graph.adjacency(ag_m,weighted=TRUE))

# SET THE GRAPH COLOR
V(agg)$color <- rep("green",length(ag))
V(agg)$name <- ag

# PLOT AGENCY GRAPH AND FILTER
plot(agg,
     layout=layout.fruchterman.reingold(agg, niter=200, area=2000*vcount(agg)),
     vertex.color="green",
     vertex.size=6,
     vertex.label=ag, 
     vertex.label.color="black", 
     vertex.label.font=1, 
     vertex.label.cex=0.5, 
     edge.width=0.5*E(agg)$weight,
     edge.curved=TRUE,
     edge.color=gray.colors(1))

# BEFORE WE FILTER, MULTILEVEL COMMUNITY DETECTION:
mc<-multilevel.community(agg)
plot(mc,
     agg, 
     vertex.size=5,
     edge.width=0.15*E(agg)$weight,
     main="Example: ML Communities",
     vertex.label.cex=0.8,
     vertex.label=V(agg)$name)

# FILTRATION, CUTOFF = 75% quantile
cut75 <- quantile(as.vector(ag_m[ag_m>0]),0.75)
agg_f<-filter(cutoff = cut75,
              edge_matrix = ag_m,
              vertex_colors = V(agg)$color,
              vertex_names = ag)

plot(as.undirected(agg_f),
     layout=layout.fruchterman.reingold(agg_f, niter=200, area=2000*vcount(agg_f)),
     vertex.color="green",
     vertex.size=10,
     vertex.label=V(agg_f)$name, 
     vertex.label.color="black", 
     vertex.label.font=1, 
     vertex.label.cex=1, 
     edge.width=(E(agg_f)$weight),
     edge.curved=TRUE,
     edge.color=gray.colors(1))

# MULTILEVEL COMMUNITY DETECTION WITH FILTRATIONS:
mc_f <- multilevel.community(as.undirected(agg_f))
plot(mc_f,
     as.undirected(agg_f), 
     vertex.size=10,
     edge.width=0.5*E(agg_f)$weight,
     main="Example: ML Communities",
     vertex.label.cex=1,
     vertex.label=V(agg_f)$name)

cut85 <- quantile(as.vector(ag_m[ag_m>0]),0.85)
agg_f<-filter(cutoff = cut85,
              edge_matrix = ag_m,
              vertex_colors = V(agg)$color,
              vertex_names = ag)

plot(as.undirected(agg_f),
     layout=layout.fruchterman.reingold(agg_f, niter=200, area=2000*vcount(agg_f)),
     vertex.color="green",
     vertex.size=10,
     vertex.label=V(agg_f)$name, 
     vertex.label.color="black", 
     vertex.label.font=1, 
     vertex.label.cex=1, 
     edge.width=(E(agg_f)$weight),
     edge.curved=TRUE,
     edge.color=gray.colors(1))

c_f <- multilevel.community(as.undirected(agg_f))
plot(mc_f,
     as.undirected(agg_f), 
     vertex.size=10,
     edge.width=0.5*E(agg_f)$weight,
     main="Example: ML Communities",
     vertex.label.cex=1,
     vertex.label=V(agg_f)$name)

cut90 <- quantile(as.vector(ag_m[ag_m>0]),0.90)
agg_f<-filter(cut90,
              edge_matrix = ag_m,
              vertex_colors = V(agg)$color,
              vertex_names = ag)

plot(as.undirected(agg_f),
     layout=layout.fruchterman.reingold(agg_f, niter=200, area=2000*vcount(agg_f)),
     vertex.color="green",
     vertex.size=10,
     vertex.label=V(agg_f)$name, 
     vertex.label.color="black", 
     vertex.label.font=1, 
     vertex.label.cex=1, 
     edge.width=0.5*(E(agg_f)$weight),
     edge.curved=TRUE,
     edge.color=gray.colors(1))

mc_f <- multilevel.community(as.undirected(agg_f))
plot(mc_f,
     as.undirected(agg_f), 
     vertex.size=10,
     edge.width=0.5*E(agg_f)$weight,
     main="Example: ML Communities",
     vertex.label.cex=1,
     vertex.label=V(agg_f)$name)

cut95 <- quantile(as.vector(ag_m[ag_m>0]),0.95)
agg_f<-filter(cut95,
              edge_matrix = ag_m,
              vertex_colors = V(agg)$color,
              vertex_names = ag)

plot(as.undirected(agg_f),
     layout=layout.fruchterman.reingold(agg_f, niter=200, area=2000*vcount(agg_f)),
     vertex.color="green",
     vertex.size=10,
     vertex.label=V(agg_f)$name, 
     vertex.label.color="black", 
     vertex.label.font=1, 
     vertex.label.cex=1, 
     edge.width=0.5*(E(agg_f)$weight),
     edge.curved=TRUE,
     edge.color=gray.colors(1))

mc_f <- multilevel.community(as.undirected(agg_f))
plot(mc_f,
     as.undirected(agg_f), 
     vertex.size=10,
     edge.width=0.5*E(agg_f)$weight,
     main="Example: ML Communities",
     vertex.label.cex=1,
     vertex.label=V(agg_f)$name)

cut97 <- quantile(as.vector(ag_m[ag_m>0]),0.97)
agg_f<-filter(cut97,
              edge_matrix = ag_m,
              vertex_colors = V(agg)$color,
              vertex_names = ag)

plot(as.undirected(agg_f),
     layout=layout.fruchterman.reingold(agg_f, niter=200, area=2000*vcount(agg_f)),
     vertex.color="green",
     vertex.size=10,
     vertex.label=V(agg_f)$name, 
     vertex.label.color="black", 
     vertex.label.font=1, 
     vertex.label.cex=1, 
     edge.width=0.5*(E(agg_f)$weight),
     edge.curved=TRUE,
     edge.color=gray.colors(1))

mc_f <- multilevel.community(as.undirected(agg_f))
plot(mc_f,
     as.undirected(agg_f), 
     vertex.size=10,
     edge.width=0.5*E(agg_f)$weight,
     main="Example: ML Communities",
     vertex.label.cex=1,
     vertex.label=V(agg_f)$name)

# CHECKING THE GIANT CONNECTED COMPONENT WE SEE IT IS CONNECTED VERY WELL
# agg_c <- as.undirected(giant_comp(agg,V(agg)$name))
# 
# plot(agg_c,
#      layout=layout.fruchterman.reingold(agg_c, niter=200, area=2000*vcount(agg_c)),
#      vertex.color="green",vertex.size=10,vertex.label=V(agg_c)$name, 
#      vertex.label.color="black", vertex.label.font=1, vertex.label.cex=1, 
#      edge.width=0.5*(E(agg_c)$weight),edge.curved=TRUE,edge.color=gray.colors(1))





# ANALYSIS OF AGENCY NETWORK: 




# RANGE OF NUMBER OF DISTINCT AID INSTANCES FOR EACH AGENCY
summary(as.data.frame(table(aid_data$impl_agency))[,2])

# NOTE: THIS IS NOT THE SAME AS
# summary(graph.strength(av))
# SINCE BOTH AGENCIES AND VDCs ARE INCLUDED IN THIS

# RANGE OF NUMBER OF DISTINCT VDCs OF AID FOR EACH AGENCY
unique_aid <- unique(cbind.data.frame(aid_data$impl_agency,aid_data$vdc))
colnames(unique_aid) <- c("impl_agency","vdc")
summary(as.data.frame(table(unique_aid$impl_agency))[,2])

# NOTE: THIS IS NOT THE SAME AS
# summary(degree(av))
# SINCE BOTH AGENCIES AND VDCs ARE INCLUDED IN THIS

# PLOT RELIEF AGENCY WEIGHTED DEGREE DISTRIBUTION (DISTINCT TYPES OF AID)
plot(sort(as.data.frame(table(aid_data$impl_agency))[,2]),
     col = adjustcolor(rgb(1,0,1,1)),
     pch = 19,
     xlab = "Agency index",
     ylab = "Numer of Distinct Aid Activities",
     main = "Sorted Agencies by Number of Distinct Aid Activities")
par(new = T)
lines(x = c(0,length(ag)),y = rep(mean(as.data.frame(table(aid_data$impl_agency))[,2]),2), col ="black", lwd=4)
text(x = 25,y = 75,paste("MEAN =",mean(as.data.frame(table(aid_data$impl_agency))[,2])),col="black",cex=2.5)

histP1(as.data.frame(table(aid_data$impl_agency))[,2],
       breaks=100,
       col = adjustcolor(rgb(1,0,1,1)),
       xlab="Agency Network Aid Action Numbers", 
       main="Agency Network Number of Aid Actions Distribution
  (VDC Overlap Counts Dsitribution)")

# PLOT RELIEF AGENCY DEGREE DISTRIBUTION (DISTINCT VDCs)
plot(sort(as.data.frame(table(unique_aid$impl_agency))[,2]),
     col = adjustcolor(rgb(1,0,1,1)),
     pch = 19,
     xlab = "Agency index",
     ylab = "Numer of Distinct Aid Activities",
     main = "Sorted Agencies by Number of Distinct VDC")

hist(as.data.frame(table(unique_aid$impl_agency))[,2], breaks=100,
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
       col = adjustcolor(rgb(1,0,1,1)),
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

bc<-betweenness(agg,v=V(agg), directed=FALSE)
bc_int <- as.integer(round(bc,0))
for (k in 1:length(bc_int)){
  V(agg)$color[k] <- rev(heat.colors(1+as.integer(max(bc_int))))[as.integer(bc_int[k])+1]
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
     edge.width=0.5*E(agg_f)$weight,
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
     edge.width=0.5*E(agg_f)$weight,
     main="Example: ML Communities",
     vertex.label.cex=0.8,
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
     vertex.size=5,
     edge.width=0.5*E(agg_f)$weight,
     main="Example: ML Communities",
     vertex.label.cex=0.8,
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
u_vdc <- as.character(unique(unique_aid$vdc))

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
     edge.color = gray.colors(1))

# BEFORE WE FILTER, COMMUNITY DETECTION:

mc<-multilevel.community(vgg)
plot(mc,vgg, vertex.size=2,edge.width=0.1*E(vgg)$weight,
     main="Example: ML Communities",
     vertex.label.cex=0.8,
     vertex.label=NA)


# FILTER BY EDGE WEIGHT LEVELS

# THIS INITIAL FILTRATION IS REDUNDNAT SINCE MINIMAL DEGREE IS 1
vgg <- as.undirected(graph.adjacency(aid_vdc,weighted=TRUE))
V(vgg)$color <- rep("SkyBlue2",length(u_vdc))
V(vgg)$name <- u_vdc
cut <- 1
vgg_f <- filter(cutoff = cut,
                edge_matrix = aid_vdc,
                vertex_colors = V(vgg)$color,
                vertex_names = u_vdc)
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

mc_f <- multilevel.community(as.undirected(vgg_f))
plot(mc_f,as.undirected(vgg_f), vertex.size=2,edge.width=0.5*E(vgg_f)$weight,
     main="Example: ML Communities",
     vertex.label.cex=1,
     vertex.label=NA)


# THE NEXT CUT IS TOO BIG OF A JUMP, WHICH WILL BE EVIDENT IN THE DEGREE DISTRIBUTION
vgg <- as.undirected(graph.adjacency(aid_vdc,weighted=TRUE))
V(vgg)$color <- rep("SkyBlue2",length(u_vdc))
V(vgg)$name <- u_vdc
cut <- 2
vgg_f <- filter(cutoff = cut,
                edge_matrix = aid_vdc,
                vertex_colors = V(vgg)$color,
                vertex_names = u_vdc)
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

mc_f <- multilevel.community(as.undirected(vgg_f))
plot(mc_f,as.undirected(vgg_f), vertex.size=3,edge.width=0.5*E(vgg_f)$weight,
     main="Example: ML Communities",
     vertex.label.cex=1,
     vertex.label=NA)


vgg <- as.undirected(graph.adjacency(aid_vdc,weighted=TRUE))
V(vgg)$color <- rep("SkyBlue2",length(u_vdc))
V(vgg)$name <- u_vdc
cut <- 3
vgg_f <- filter(cutoff = cut,
                edge_matrix = aid_vdc,
                vertex_colors = V(vgg)$color,
                vertex_names = V(vgg)$name)
vgg_f <- giant_comp(graph = vgg_f,
                    vertex_colors = V(vgg_f)$color,
                    vertex_names = V(vgg_f)$name)
plot(as.undirected(vgg_f),
     layout=layout.fruchterman.reingold(vgg_f, niter=200, area=2000*vcount(vgg_f)),
     vertex.color = V(vgg_f)$color,
     vertex.size=5,
     vertex.label=V(vgg_f)$name, 
     vertex.label.color="black", 
     vertex.label.font=1, 
     vertex.label.cex=1, 
     edge.width=(E(vgg_f)$weight),
     edge.curved=TRUE,
     edge.color=gray.colors(1))

mc_f <- multilevel.community(as.undirected(vgg_f))
plot(mc_f,
     as.undirected(vgg_f), 
     vertex.size=3,
     edge.width=0.5*E(vgg_f)$weight,
     main="Example: ML Communities",
     vertex.label.cex=1,
     vertex.label=V(vgg_f)$name)


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
       xlab = "Target VDC Network Degree Values", 
       main = "Target VDC Network Degree Distribution")


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
       col = adjustcolor(rgb(1,0,1,1)),
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







# CHECK FROM HERE ON








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
     vertex.size = 7,
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
     edge.width=0.25*E(vgg_f)$weight,
     main="Example: ML Communities",
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
ec<-evcent(vgg1)$vector
plot(sort(ec, decreasing=TRUE), col=adjustcolor(rgb(0,0,1,1/2)), xlab="Node Id in the Network (1:200)", ylab="Closeness Centrality Values", main="Essential (first 200 nodes) Closeness Centrality for g", pch=20)
hist(ec,breaks=100,col=adjustcolor(rgb(0,0,1,1/2)),xlab="Closeness Centrality Values",main="Essential Closeness Centrality Distribution")


# DISPLAY THE GRAPH WITH THE APPORPRIATE COLORS FOR THE VERTICES
vgg <- as.undirected(graph.adjacency(aid_vdc,weighted=TRUE))
V(vgg)$color <- rep("green",length(u_vdc))
V(vgg)$name <- u_vdc
ec <- closeness(graph = vgg,
                vids = V(vgg), 
                weights = E(vgg)$weight,
                normalized = TRUE)
ec_int <- as.integer(round(10000*ec,0))
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
ec <- closeness(graph = vgg,
                vids = V(vgg), 
                weights = E(vgg)$weight,
                normalized = TRUE)
ec_int <- as.integer(round(10000*ec,0))
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
ec <- closeness(graph = vgg_f,
                vids = V(vgg_f), 
                weights = E(vgg_f)$weight,
                normalized = TRUE)
ec_int <- as.integer(round(10000*ec,0))
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
ec <- closeness(graph = vgg_f,
                vids = V(vgg_f), 
                weights = E(vgg_f)$weight,
                normalized = TRUE)
ec_int <- as.integer(round(10000*ec,0))
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




# COMMUNITY STRUCTURES: This is a way of performing funcitonal clustering in complex networks. We have already looked at the connected components, 
# this is an elementary community detection based on connectivity.
strongclusters<-clusters(vgg)$membership
plot(vgg,vertex.color=strongclusters, layout=layout.fruchterman.reingold,vertex.size=4, edge.color="black", edge.width=E(vgg)$weight,vertex.label=NA,main="Clustering for Store Network g200")

# ADD SOME FILTERING AND TRY AGAIN

mc<-multilevel.community(vgg)
plot(sort(mc$membership), col=adjustcolor(rgb(0,0,1,1/2)), xlab="Node Id in the Network", ylab="Multilevel Community Values", main="Multilevel Community Values for g", pch=20)
hist(mc$membership,breaks=100,col=adjustcolor(rgb(0,0,1,1/2)),xlab="Multilevel Community Values",main="Multilevel Community Distribution")


# Next, we show the walktrap community algorithm.
wc<-walktrap.community(vgg)
plot(sort(wc$membership), col=adjustcolor(rgb(0,0,1,1/2)), xlab="Node Id in the Network", ylab="Walktrap Community Values", main="Walktrap Community Values for g", pch=20)
hist(wc$membership,breaks=100,col=adjustcolor(rgb(0,0,1,1/2)),xlab="Walktrap Community Values",main="Walktrap Community Distribution")


plot(wc,vgg,vertex.size=4, vertex.label=NA,edge.width=E(vgg)$weight,main="Walktrap Community Detection for g200")
plot(vgg, vertex.color=membership(wc), vertex.size=6, edge.color="black", edge.width=E(vgg)$weight,vertex.label=NA,main="Walktrap Community Detection for g200")




















































# FOR THE AGENCIES: SIMILARITY MEASURE BASED ON AID TYPE AND QUANTITY

# FOR THE VDC NETWORK: DIRECTED DISPLACEMENT TRACKING GRAPH
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






SECTION 2: Nepal Displacement Tracking Network Analytics Examples.
library("igraph")
library("stats")
# We define all the graphs that we will be using in this research.
g100<-graph.adjacency(edgem[1:100,1:100],mode="undirected",weighted=TRUE)
g200<-graph.adjacency(edgem[1:200,1:200],mode="undirected",weighted=TRUE)
g201<-graph.adjacency(edges[1:200,1:200],mode="undirected")
g400<-graph.adjacency(edgem[1:400,1:400],mode="undirected",weighted=TRUE)
g401<-graph.adjacency(edges[1:400,1:400],mode="undirected")
g1000<-graph.adjacency(edgem[1:1000,1:1000],mode="undirected",weighted=TRUE)
g1001<-graph.adjacency(edges[1:1000,1:1000],mode="undirected")
g2000<-graph.adjacency(edgem[1:2000,1:2000],mode="undirected",weighted=TRUE)
g2001<-graph.adjacency(edges[1:2000,1:2000],mode="undirected")
g<-graph.adjacency(edgem,mode="undirected",weighted=TRUE)
g1<-graph.adjacency(edges,mode="undirected")
#We next plot the examples of store networks for stores 1-100 and 1-200. Note that since many stores are not presented, there are many edges missing.")
plot(g100, layout=layout.fruchterman.reingold,vertex.size=8, edge.color="black", edge.width=E(g100)$weight,vertex.label.cex=0.8,main="Store Network for g100",vertex.label=ids[1:100])
plot(g200, layout=layout.fruchterman.reingold,vertex.size=5, edge.color="black", edge.width=E(g200)$weight,vertex.label.cex=0.3,main="Store Network for g200",vertex.label=ids[1:200])


SECTION 3: Nepal Displacement Tracking Network Analytics for the Complete Store Network

#Now we will show a few basic properties of the Wal Mart Store Trade Area network. We will use various sub-networks (g200, g400, g1000, g2000) to visualize some of the results.")
plot(g, layout=layout.fruchterman.reingold,vertex.label=NA,vertex.size=2,main="Wal Mart Store Network",vertex.label=ids)

NUMBER OF NODES: This is the number of stores (nodes) in this network.

length(V(g))
NUMBER OF EDGES: This is the number of store-to-store connections in this network. Each connection between stores indicates shared trade areas between the stores.
length(E(g))

GRAPH DENSITY: This is a measure of the density of edges on the graph (Number of edges/Number of possible edges) in this network.
graph.density(g)

DEGREE OF A NODE: This shows the number of connections to the node. That is, the number of different stores that a store shares trade areas with. In further analysis, we will take into account the different stores' percent sales that come from a shared trade area. 
summary(degree(g))
NUMBER OF CONNECTED COMPONENTS: The number of connected components or "islands" or "clusters"
clusters(g)$no
sort(clusters(g)$csize,decreasing=TRUE)[1:80]
GIANT CONNECTED COMPONENT: The largest cluster or island in the network.
#Notice the giant connected component (cluster) of about 3500 stores. This is subject of investigation: when did it arise, what is the store demographic distribution within it, how are stores clustered within it (clique anlaysis).") 
sort(clusters(g)$csize, decreasing=TRUE)[1]
#The giant connected component as a percentage:")
(sort(clusters(g)$csize, decreasing=TRUE)[1])/vcount(g)
NUMBER OF ISOLATED NODES:
sum(degree(g)==0)
#Number of isolated nodes as a percentage:")
(sum(degree(g)==0))/vcount(g)
GLOBAL CLUSTERING COEFFICIENT: This is a measure of the clustering of the network (Number of of triangles of closed triplets/ Number of total triples (both open and closed)), it is also called the transitivity measure. For weighted networks like ours, there are more than four refinements of this measure, and we will choose each depending on the question we are looking to answer in subsequent reports. For now, we compute the generic one (ignoring weights):
transitivity(g)
#For comparison, transitivity of a random graph of the same size:")
ge<-erdos.renyi.game(vcount(g),ecount(g),type="gnm")
transitivity(ge)
EXAMPLES OF CLUSTERING: g200, g400, g2000
sort(clusters(g201)$csize,decreasing=TRUE)
strongclusters<-clusters(g200)$membership
plot(g200,vertex.color=strongclusters, layout=layout.fruchterman.reingold,vertex.size=4, edge.color="black", edge.width=E(g200)$weight,vertex.label=NA,main="Clustering for Store Network g200")
sort(clusters(g400)$csize,decreasing=TRUE)
strongclusters<-clusters(g400)$membership
plot(g400,vertex.color=strongclusters, layout=layout.fruchterman.reingold,vertex.size=2, edge.color="black", edge.width=E(g400)$weight,vertex.label=NA,main="Clustering for Store Network g400")
sort(clusters(g2001)$csize, decreasing=TRUE)[1:90]
strongclusters<-clusters(g2001)$membership
plot(g2001,vertex.color=strongclusters, layout=layout.fruchterman.reingold,vertex.size=2, vertex.label=NA, edge.color="black",main="Clustering for Store Network g2000")
EDGE-CONNECTIVITY: This is a measure of whether the network is connected or not. Here, the value is 0 since the network is disconnected.
edge.connectivity(g)
DIAMETER OF THE NETWORK: The longest path (among all connected components). We will use the unweighted graph first so that we can count the longest unweighted path in an intuitive way. For the weighted graph, we will consider a path as the weighted sum of paths. For the definition of a weigthed path, we take the sum of the weights of each of the edges along the path. We illustrate the diameters of g400 and g2000 as examples.
#We compute the diameter of the complete network first.")
diameter(g1)
#Now we illustrate the examples for g400 and g200.")
d<-get.diameter(g401)
E(g401)$color<-"black"
E(g401)$width<-1
E(g401,path=d)$color<-"red"
E(g401,path=d)$width<-2
V(g401)$label.color<-"blue"
V(g401)$color<-"SkyBlue2"
V(g401)[d]$label.color<-"black"
V(g401)[d]$color<-"red"
plot(g401,layout=layout.fruchterman.reingold,vertex.label.dist=0,vertex.size=5,vertex.label.cex=0.4, main="Diameter for the Largest Connected Component of g400 (unweighted)",vertex.label=ids[1:400])
d<-get.diameter(g2001)
E(g2001)$color<-"black"
E(g2001,path=d)$color<-"red"
E(g2001,path=d)$width<-2
V(g2001)$color<-"SkyBlue2"
V(g2001)[d]$color<-"red"
plot(g2001,layout=layout.fruchterman.reingold,vertex.label=NA,vertex.size=2, main="Diameter for the Largest Connected Component of g2000")
DEGREE DISTRIBUTION: This studies the distribution of node connections across the network.The functionality in igraph has a degree distribution function which gives cumulative normalized degree distribution, but we will use the absolute degree distribution frequence definition (count how many times a certain node degree occurs.) 
degree.distribution(g)
# We illustrate the degree distribution of the subgraphs of g200 and g400 for comparison.
par(mfrow=c(1,2))
hist(degree(g201), breaks=60, col=adjustcolor(rgb(1,0,1,1)), xlab="Node Degree", main="Node Degree Distribution Histogram for g200")
# Just for comparison, we will compute the degree distribution of g400.
hist(degree(g401), breaks=100, col=adjustcolor(rgb(1,0,1,1)), xlab="Node Degree", main="Node Degree Distribution Histogram for g400")
#Note that the plots are rescaled, so in g400, we have more than twice as many nodes with degree 1 or 2, for example. The main point to note is the similar shape of the degree distribution. This is a property of networks (real-life networks) which is known as scale-free, that is, we notice similar network structure at different scales.")
#Typically, the degree distributions are compared in a semi-log graph, because resutls from complex network science show certian networks have exponential degree distribution (and ours is one of them).")
d21<-hist(degree(g201))$counts
d22<-d21[which(!d21==0)]
d200<-cbind.data.frame(as.vector(1:length(d22)),log(d22))
colnames(d200)<-c("x","y")
fit2<-lm(y~x, data=d200)
d41<-hist(degree(g401))$counts
d42<-d41[which(!d41==0)]
d400<-cbind.data.frame(as.vector(1:length(d42)),log(d42))
colnames(d400)<-c("x","y")
fit4<-lm(y~x, data=d400)
#Notice that the slopes (the exponent of the degree distribution) are very close in value (compare with the slope of the degree distribution of the whole network below).")
par(mfrow=c(1,2))
plot(d200$x,d200$y,xlab="Degree", ylab="Log[Degree Frequency]", pch=16, main="Semi-Log Degree Distribution with a Linear Fit for g200")
lines(abline(coef=coef(fit2), col="red"),xlim=range(0:length(d22)),ylim=range(0:max(log(d22))))
legend(5.5,4.25,legend=rbind(c("Intercept:","Slope:"),coef(fit2)))
plot(d400$x,d400$y,xlab="Degree", ylab="Log[Degree Frequency]", pch=16, main="Semi-Log Degree Distribution with a Linear Fit for g400")
lines(abline(coef=coef(fit4), col="red"),xlim=range(0:length(d42)),ylim=range(0:max(log(d42))))
legend(6,4.5,legend=rbind(c("Intercept:","Slope:"),coef(fit4)))
#For the degree distribution analysis of the complete network, we use the unweighted version, for simplicity. When the weights are taken into consideration, degree is a weighted sum of all connections, each one weighted by its, well, weight. In some cases, we normalize all weights, and that makes some of the figures difficult to interpret. That is why, for simplicity, we will just cound the connections, and use that as the definition of a degree.")
d1<-hist(degree(g1))$counts
d2<-d1[which(!d1==0)]
dg<-cbind.data.frame(1:length(d2),log(d2))
colnames(dg)<-c("x","y")
fit<-lm(y~x, data=dg)
par(mfrow=c(1,2))
hist(degree(g), breaks=400, col=adjustcolor(rgb(1,0,1,1)), xlab="Node Degree", main="Node Degree Distribution for the Unweighted Store Network")
plot(dg$x,dg$y,xlab="Degree", ylab="Log[Degree Frequency]", pch=16, main="Degree Distribution for g")
lines(abline(coef=coef(fit), col="red"),xlim=range(0:length(d2)),ylim=range(0:max(log(d2))))
legend(8,7.5,legend=rbind(c("Intercept:","Slope:"),coef(fit)))
MINIMUM WEIGHT SPANNING TREE: This is the minimum number of edges that we can keep so that connected components are preserved, vertices are preserved, and the weight is minimized. This can be useful in order to show us the most essential shared trade areas and their adjacent stores. We observe that most of the edges are needed in order to preserve the connectivity of the graph. The real power of this tool comes in when we have higher connectivity and essential features (the "spine of the graph") need to be detected. We will use g400 and g1000 to illustrate the concept.
E(g401)$id<-seq_len(ecount(g401))
mst<-minimum.spanning.tree(g401)
E(g401)$color<-"black"
E(g401)$width<-1
E(g401)$color[E(mst)$id]<-"red"
E(g401)$width[E(mst)$id]<-1
V(g401)$color<-"SkyBlue2"
plot(g401,layout=layout.fruchterman.reingold,vertex.label=NA,vertex.size=2,main="Minimum Spanning Tree for the Largest Connected Component of g400")
#We compute the minimum weight spanning tree for g1000.")
E(g1001)$id<-seq_len(ecount(g1001))
mst<-minimum.spanning.tree(g1001)
E(g1001)$color<-"black"
E(g1001)$color[E(mst)$id]<-"red"
V(g1001)$color<-"SkyBlue2"
plot(g1001,layout=layout.fruchterman.reingold,vertex.label=NA,vertex.size=2,main="Minimum Spanning Tree for the Largest Connected Component of the Network g1000")
