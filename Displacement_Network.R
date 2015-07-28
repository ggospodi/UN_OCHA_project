# Nepal Displacement Tracking Network Analysis
# author: Georgi D. Gospodinov
# date: "July 27, 2015"
# 
# Data Sources:
#
# Tables:
# 
# https://data.hdx.rwlabs.org/dataset/nepal-earthquake-severity-index
# https://data.hdx.rwlabs.org/dataset/io (sheet 3)
# https://data.hdx.rwlabs.org/dataset/population-movements-after-the-nepal-earthquake-v-3-up-to-11th-june-2015
#
# 
# Relevant materials and data cna be found at:
# 
# https://www.dropbox.com/sh/tb9854hzcof7x23/AACEDTGk8EmYQ6r4ukSFLBspa?dl = 0
#
# in the folder /Displacement Network
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
  df[,which(names(df) %in% col_name)] <- mapvalues(df[,which(names(df) %in% col_name)], from=level_names,to=gsub("[[:space:]]","",level_names))
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


# INITIAL DATA INPUT AND TRANSFORMATIONS
#
#
# COLUMN NAMES FOR CCCM_Nepal_DTM_R2.csv
#
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


# LOAD LAT/LON COORDINATES (OF CENTROIDS) AND HLCIT CODES
hlcit <- read.csv(paste0(DIR,"master_hlcit.csv"))
colnames(hlcit) <- c("lon","lat","vdc_name","vname","hlcit_code")
hlcit$hlcit_code <- as.factor(hlcit$hlcit_code)
hlcit$vname <- as.character(hlcit$vname)
hlcit$vdc_name <- as.character(hlcit$vdc_name)
hlcit <- rm_space(hlcit,"hlcit_code")
hlcit$hlcit_code <- as.numeric(levels(hlcit$hlcit_code))[hlcit$hlcit_code]


# READ IN THE DISPLACEMENT FILE DATA
dt_data <- read.csv(paste0(DIR,"CCCM_Nepal_DTM_R2.csv"), sep = ",")


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
write.csv(dt_data,file = paste0(DIR,"Nepal_DTM_Network_model1.csv"))


# READ IN DISPLACEMENT DATA
dt_data <- read.csv(paste0(DIR,"Nepal_DTM_Network_model1.csv"))
dt_data <- dt_data[,2:dim(dt_data)[2]]


# FORMAT THE VDC RELATED COLUMNS AS CHARACTER
dt_data$vdc <- trim(as.character(dt_data$vdc))
dt_data$idp_origin_vdc <- trim(as.character(dt_data$idp_origin_vdc))
dt_data$idp2_origin_vdc <- trim(as.character(dt_data$idp2_origin_vdc))


# FILTER TO USE DESTINATION VDC LEVELS AS WELL AS ORIGIN (1 OR 2) VDC ENTRIES THAT ARE NON-EMPTY
dt_data <- dt_data[nchar(dt_data$vdc)>0 & (nchar(dt_data$idp_origin_vdc) + nchar(dt_data$idp2_origin_vdc))>0,]


# RESOLVE VDC NAMES ACCORDING TO THE HLCIT MASTER LIST:
dt_data$vdc <- mapvalues(dt_data$vdc,
                         from = c("Barahbise","Barpak","Chandeni","Charikot","Gokarneshwar","Kathmandu",
                                  "Kathmandu Municipality","Lalitpur","Mankha","Manmaijn","Naikap Naya",
                                  "Puranagaun","Sanga","Sangkhu Suntol","Tokhasaraswati"),
                         to = c("Barhabise","Warpak","Chandeni Mandan","Narikot","Gokarneswor",
                                "Kathmandu Metropolitan","Kathmandu Metropolitan","Lalitpur Sub Metropolitan",
                                "Mangkha","Manmaiju","NaikapNayaBhanjyang","PuranogaunDapcha",
                                "Sangla","Sangkhu","TokhaSarswoti"))
dt_data$idp_origin_vdc <- mapvalues(dt_data$idp_origin_vdc,
                                    from = c("barabise","Bhimeshwor Municipality","Chautara Municipality","Maanka"),
                                    to = c("Barhabise","Bhimeswor Municipality","Chautara","Mahangkal"))
dt_data$idp2_origin_vdc <- mapvalues(dt_data$idp2_origin_vdc,
                                     from = c("barabise","Bhimeshwor Municipality","Jhaku"),
                                     to = c("Barhabise","Bhimeswor Municipality","Jhyanku"))


# CREATE A LIST OF UNIQUE VDC DESTINATION NAMES
vdc <- unique(c(dt_data$vdc, dt_data$idp_origin_vdc,dt_data$idp2_origin_vdc))
vdc <- vdc[-which(vdc == "")]
vdc1 <- intersect(vdc,hlcit$vdc_name)
vdc2 <- setdiff(vdc,hlcit$vdc_name)


# RESOLVE THE ONLY VDC NAME LEFT THAT IS NOT ON RECORD
index <- which(dt_data$vdc == setdiff(vdc2,hlcit$vname))
closest <- vector()
for (k in 1:dim(hlcit)[1]){
  closest[k] <- (dt_data$lat[index]-hlcit$lat[k])^2+(dt_data$lon[index]-hlcit$lon[k])^2
}
dt_data$vdc[index] <- hlcit$vname[which(closest == min(closest))]
vdc <- unique(c(dt_data$vdc, dt_data$idp_origin_vdc,dt_data$idp2_origin_vdc))
vdc <- vdc[-which(vdc == "")]
vdc1 <- intersect(vdc,hlcit$vdc_name)
vdc2 <- setdiff(vdc,hlcit$vdc_name)


# EXTRAPOLATE LHCIT NUMBERS FROM LAT AND LON FOR VDCS FROM HLCIT_MASTER
hl <- vector()
xc <- vector()
yc <- vector()
for (k in 1:length(vdc)){
  if (vdc[k] %in% vdc1){
    hl[k] <- hlcit$hlcit_code[which(hlcit$vdc_name == vdc[k])[1]]
    xc[k] <- hlcit$lat[which(hlcit$vdc_name == vdc[k])[1]]
    yc[k] <- hlcit$lon[which(hlcit$vdc_name == vdc[k])[1]]
  }
  if (vdc[k] %in% vdc2){
    hl[k] <- hlcit$hlcit_code[which(hlcit$vname == vdc[k])[1]]
    xc[k] <- hlcit$lat[which(hlcit$vname == vdc[k])[1]]
    yc[k] <- hlcit$lon[which(hlcit$vname == vdc[k])[1]]
  }
}
koords <- cbind(xc,yc)


# FURTHER NARROW DOWN COLUMNS
dt_data <- dt_data[,c("vdc","idp_origin_vdc","idp2_origin_vdc","idp_hh")]


# VDC-LEVEL ADJACENCY MATRIX FOR THE DISPLACEMENT GRAPH, AT THE LEVEL OF DISPLACEMENT TRACK
vdc_m <- matrix(0,nrow = length(vdc),ncol = length(vdc))


# CREATE A MATRIX THAT TRACKS THE NUMBERS OF DISPLACED POPULATIONS AS WELL 
# (THIS COULD BE USED TO DERIVE EDGE WEIGHTS)
dtm <- matrix(nrow = length(vdc),ncol = length(vdc))


# COMPUTE THE ENTRIES OF BOTH THE TRACKING DISPLACEMENT MATRIX 
# AND THE WEIGHTED DISPLACEMENT TRACKING MATRIX
for (i in 1:length(vdc)){
  for (j in 1:length(vdc)){
    vdc_m[[i,j]] <- length(dt_data[dt_data$idp_origin_vdc == vdc[i] & dt_data$vdc == vdc[j],1])+
      length(dt_data[dt_data$idp2_origin_vdc == vdc[i] & dt_data$vdc == vdc[j],1])
    dtm[[i,j]] <- (2/3)*sum(dt_data[dt_data$idp_origin_vdc == vdc[i] & dt_data$vdc == vdc[j],]$idp_hh)+
      (1/3)*sum(dt_data[dt_data$idp2_origin_vdc == vdc[i] & dt_data$vdc == vdc[j],]$idp_hh)
  }
}


# BUILD THE DIRECTED WEIGHTED VDC NETWORK
gv <- graph.adjacency(vdc_m,
                      mode = "directed",
                      weighted = TRUE)


# COLOR VDC NAMES OF ORIGIN (GREEN) AND VDC NAMES OF DESTINATION (BLUE)
V(gv)$color <- rep("SkyBlue2",length(vdc))
for (k in 1:length(vdc)){
  o_count <- length(which(dt_data$idp_origin_vdc == vdc[k]))+
    length(which(dt_data$idp2_origin_vdc == vdc[k]))
  d_count <- length(which(dt_data$vdc == vdc[k]))
  if(o_count>d_count){
    V(gv)$color[k] <- "green"
  } 
}


# SET THE VERTEX LABELS
V(gv)$name <- vdc


# DROP ISOLATED VERTICES (NO DISPLACEMENT)
gv <- drop_isolated(graph = gv,
                    vertex_colors = V(gv)$color,
                    vertex_names= vdc)
plot(gv,
     layout = layout.fruchterman.reingold(gv,
                                          niter = 200,
                                          area = 2000*vcount(gv)),
     vertex.color = V(gv)$color,
     vertex.size = 9, 
     vertex.label = V(gv)$name,
     vertex.label.color = "black", 
     vertex.label.font = 1, 
     vertex.label.cex = 0.9, 
     edge.width = 0.2*(E(gv)$weight),
     edge.arrow.size = 0.5,
     edge.curved = TRUE,
     edge.color = gray.colors(1),
     main = "Abstract Nepal Displacement Network Flow (VDC Level, with Self-Loops)")
legend("topright",
       c("Displacement Origin","Displacement Destination"),
       fill = c("green","SkyBlue2"),
       bty = "n")


# DROP LOOPS ONLY VERTICES AS WELL ()
gv <- drop_loops(graph = gv,
                 vertex_colors = V(gv)$color,
                 vertex_names = vdc)

# RESULTING CLEANED UP GRAPH SHOWING NONTRIVIAL MIGRATION
plot(gv,
     layout = layout.fruchterman.reingold(gv,
                                          niter = 200,
                                          area = 2000*vcount(gv)),
     vertex.color = V(gv)$color,
     vertex.size = 9, 
     vertex.label = V(gv)$name,
     vertex.label.color = "black", 
     vertex.label.font = 1, 
     vertex.label.cex = 0.9, 
     edge.width = 2*(E(gv)$weight),
     edge.arrow.size = 0.5,
     edge.curved = TRUE,
     edge.color = gray.colors(1),
     main = "Abstract Nepal Displacement Network Flow (VDC Level, No Self-Loops)")
legend("topright",
       c("Origins of Displacement","Destinations of Displacement"),
       fill = c("green","SkyBlue2"),
       bty = "n")


# PLOT THE GRAPH WITH SELECTED COORDINATES
gv_coords <- koords[which(vdc %in% V(gv)$name),]
plot(gv,
     layout = gv_coords,
     vertex.color = V(gv)$color,
     vertex.size = 4, 
     vertex.label = V(gv)$name,
     vertex.label.color = "black", 
     vertex.label.font = 1, 
     vertex.label.cex = 0.7, 
     edge.width = (E(gv)$weight),
     edge.arrow.size = 0.4,
     edge.curved = TRUE,
     edge.color = gray.colors(1),
     main = "Nepal Displacement Geo-Network Flow (VDC Level)")
legend("topright",
       c("Origins of Displacement","Destinations of Displacement"),
       fill = c("green","SkyBlue2"),
       bty = "n")


# DEFINE THE WEIGHTED DISPLACEMENT GRAPH
gd <- graph.adjacency(dtm,
                      mode = "directed",
                      weighted = TRUE)

# SET VERTEX COLORS
V(gd)$color <- rep("SkyBlue2",length(vdc))
for (k in 1:length(vdc)){
  o_count <- length(which(dt_data$idp_origin_vdc == vdc[k]))+
    length(which(dt_data$idp2_origin_vdc == vdc[k]))
  d_count <- length(which(dt_data$vdc == vdc[k]))
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


# RESULTING CLEANED UP GRAPH SHOWING NONTRIVIAL MIGRATION
plot(gd, 
     layout = layout.fruchterman.reingold(gd,
                                          niter = 200,
                                          area = 2000*vcount(gd)),
     vertex.color = V(gd)$color,
     vertex.size = 9, 
     vertex.label = V(gd)$name,
     vertex.label.color = "black", 
     vertex.label.font = 2, 
     vertex.label.cex = 0.7, 
     edge.width = 0.02*(E(gd)$weight),
     edge.arrow.size = 0.7,
     edge.curved = TRUE,
     edge.color = gray.colors(1),
     main = "Weighted Abstract Nepal Displacement Network Flow (VDC Level)")
legend("topright",
       c("Displacement Origin","Displacement Destination"),
       fill = c("green","SkyBlue2"),
       bty = "n")


# SELECT THE COORDINATES
gd_coords <- koords[which(vdc %in% V(gd)$name),]
plot(gd,
     layout = gd_coords,
     vertex.color = V(gd)$color,
     vertex.size = 4, 
     vertex.label = V(gd)$name,
     vertex.label.color = "black", 
     vertex.label.font = 1, 
     vertex.label.cex = 0.6, 
     edge.width = 0.2*sqrt(E(gd)$weight),
     edge.arrow.size = 0.2,
     edge.curved = TRUE,
     edge.color = gray.colors(1),
     main = "Weighted Nepal Displacement Geo-Network Flow (VDC Level)")
legend("topright",
       c("Displacement Origin","Displacement Destination"),
       fill = c("green","SkyBlue2"),
       bty = "n")


# DISPLAY THE LARGEST CLUSTER (GIANT COMPONENT):
gd_c <- giant_comp(graph = gd,
                   vertex_colors = V(gd)$color,
                   vertex_names = V(gd)$name)


# PLOT THE WEIGHTED DISPLACEMENT GRAPH
gd_c_coords <- koords[which(V(gd)$name %in% V(gd_c)$name),]
plot(gd_c,
     layout = gd_c_coords,
     vertex.color = V(gd_c)$color,
     vertex.size = 5, 
     vertex.label = V(gd_c)$name,
     vertex.label.color = "black", 
     vertex.label.font = 1, 
     vertex.label.cex = 0.7, 
     edge.width = 0.2*sqrt(E(gd_c)$weight),
     edge.arrow.size = 0.5,
     edge.curved = TRUE,
     edge.color = gray.colors(1),
     main = "Weighted Nepal Displacement Network Flow (VDC Level)")
legend("top",
       c("Displacement Origin","Displacement Destination"),
       fill = c("green","SkyBlue2"),
       bty = "n")
plot(gd_c,
     layout = layout.fruchterman.reingold(gd_c,
                                          niter = 200,
                                          area = 2000*vcount(gd_c)),
     vertex.color = V(gd_c)$color,
     vertex.size = 7, 
     vertex.label = V(gd_c)$name,
     vertex.label.color = "black", 
     vertex.label.font = 1, 
     vertex.label.cex = 1, 
     edge.width = 0.2*sqrt(E(gd_c)$weight),
     edge.arrow.size = 0.7,
     edge.curved = TRUE,
     edge.color = gray.colors(1),
     main = "Weighted Nepal Displacement Network Flow (VDC Level)")
legend("top",
       c("Displacement Origin","Displacement Destination"),
       fill = c("green","SkyBlue2"),
       bty = "n")


# EDGE-FILTRATION BY EDGE WEIGHT OF THE WEIGHTED DISPLACEMENT GRAPH: CUT-OFF = 50% quantile
cut50 <- quantile(as.vector(dtm[dtm>0]),0.50)
gd <- graph.adjacency(dtm,
                      mode = "directed",
                      weighted = TRUE)
V(gd)$name <- vdc
V(gd)$color <- rep("SkyBlue2",length(vdc))
for (k in 1:length(vdc)){
  o_count <- length(which(dt_data$idp_origin_vdc == vdc[k]))
    + length(which(dt_data$idp2_origin_vdc == vdc[k]))
  d_count <- length(which(dt_data$vdc == vdc[k]))
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
gd_f_coords <- koords[which(vdc %in% V(gd_f)$name),]
plot(gd_f,
     layout = gd_f_coords,
     vertex.color = V(gd_f)$color,
     vertex.size = 5, 
     vertex.label = V(gd_f)$name,
     vertex.label.color = "black", 
     vertex.label.font = 1, 
     vertex.label.cex = 0.7, 
     edge.width = 0.2*sqrt(E(gd_f)$weight),
     edge.arrow.size = 0.5,
     edge.curved = TRUE,
     edge.color = gray.colors(1),
     main = "Filtered Abstract Nepal Displacement Network (Cut-Off = 50%))")
legend("topright",
       c("Origins of Displacement","Destinations of Displacement"),
       fill = c("green","SkyBlue2"),
       bty = "n")
plot(gd_f,
     layout = layout.fruchterman.reingold(gd_f,
                                          niter = 200,
                                          area = 2000*vcount(gd_f)),
     vertex.color = V(gd_f)$color,
     vertex.size = 10,
     vertex.label = V(gd_f)$name, 
     vertex.label.color = "black", 
     vertex.label.font = 1, 
     vertex.label.cex = 1.2, 
     edge.width = 0.4*sqrt(E(gd_f)$weight),
     edge.arrow.size = 01,
     edge.curved = TRUE,
     edge.color = gray.colors(1),
     main = "Filtered Nepal Displacement Geo-Network (Cut-Off = 50%))")
legend("topright",
       c("Origins of Displacement","Destinations of Displacement"),
       fill = c("green","SkyBlue2"),
       bty = "n")


# EDGE-FILTRATION BY EDGE WEIGHT OF THE WEIGHTED DISPLACEMENT GRAPH: CUT-OFF = 75% quantile
cut75 <- quantile(as.vector(dtm[dtm>0]),0.75)
gd <- graph.adjacency(dtm,
                      mode = "directed",
                      weighted = TRUE)
V(gd)$name <- vdc
V(gd)$color <- rep("SkyBlue2",length(vdc))
for (k in 1:length(vdc)){
  o_count <- length(which(dt_data$idp_origin_vdc == vdc[k]))
  + length(which(dt_data$idp2_origin_vdc == vdc[k]))
  d_count <- length(which(dt_data$vdc == vdc[k]))
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
gd_f_coords <- koords[which(vdc %in% V(gd_f)$name),]
plot(gd_f,
     layout = layout.fruchterman.reingold(gd_f,
                                          niter = 200,
                                          area = 2000*vcount(gd_f)),
     vertex.color = V(gd_f)$color,
     vertex.size = 14,
     vertex.label = V(gd_f)$name, 
     vertex.label.color = "black", 
     vertex.label.font = 1, 
     vertex.label.cex = 1.4, 
     edge.width = 0.05*(E(gd_f)$weight),
     edge.arrow.size = 1.4,
     edge.curved = TRUE,edge.color = gray.colors(1),
     main = "Filtered Abstract Nepal Displacement Network (Cut-Off = 75%))")
legend("topright",
       c("Origins of Displacement","Destinations of Displacement"),
       fill = c("green","SkyBlue2"),
       bty = "n")
plot(gd_f,
     layout = layout.fruchterman.reingold(gd_f,
                                          niter = 200,
                                          area = 2000*vcount(gd_f)),
     vertex.color = V(gd_f)$color,
     vertex.size = 14,
     vertex.label = V(gd_f)$name, 
     vertex.label.color = "black", 
     vertex.label.font = 1, 
     vertex.label.cex = 1.4, 
     edge.width = 0.05*(E(gd_f)$weight),
     edge.arrow.size = 1.4,
     edge.curved = TRUE,edge.color = gray.colors(1),
     main = "Filtered Nepal Displacement Geo-Network (Cut-Off = 75%))")
legend("topright",
       c("Origins of Displacement","Destinations of Displacement"),
       fill = c("green","SkyBlue2"),
       bty = "n")

#
#
# ANALYSIS OF THE VDC NETWORK ITSELF
#
#

# DEFINE THE WEIGHTED DISPLACEMENT GRAPH
gd <- graph.adjacency(dtm,
                      mode = "directed",
                      weighted = TRUE)


# SET VERTEX COLORS
V(gd)$color <- rep("SkyBlue2",length(vdc))
for (k in 1:length(vdc)){
  o_count <- length(which(dt_data$idp_origin_vdc == vdc[k]))+
    length(which(dt_data$idp2_origin_vdc == vdc[k]))
  d_count <- length(which(dt_data$vdc == vdc[k]))
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
     layout = gd_coords,
     vertex.color = V(gd)$color,
     vertex.size = 4, 
     vertex.label = V(gd)$name,
     vertex.label.color = "black", 
     vertex.label.font = 1, 
     vertex.label.cex = 0.5, 
     edge.width = 0.15*sqrt(E(gd)$weight),
     edge.arrow.size = 0.5,
     edge.curved = TRUE,
     edge.color = gray.colors(1),
     main = "Weighted Nepal Displacement Network Flow (VDC Level)")
legend("topright",
       c("Displacement Origin","Displacement Destination"),
       fill = c("green","SkyBlue2"),
       bty = "n")

# REGENERATE THE GRAPH AGAIN
gd <- graph.adjacency(dtm,
                      mode = "directed",
                      weighted = TRUE)


# SET VERTEX COLORS
V(gd)$color <- rep("SkyBlue2",length(vdc))
for (k in 1:length(vdc)){
  o_count <- length(which(dt_data$idp_origin_vdc == vdc[k]))+
    length(which(dt_data$idp2_origin_vdc == vdc[k]))
  d_count <- length(which(dt_data$vdc == vdc[k]))
  if(o_count>d_count){
    V(gd)$color[k] <- "green"
  } 
}


# SET THE VERTEX LABELS
V(gd)$name <- vdc

# THIS IS THE NUMBER OF EDGES FROM EACH NODE
# IGNORING DIRECTION (TOTAL DEGREE)
summary(degree(gd))


# THIS IS THE IN DEGREE SUMMARY
summary(degree(gd,mode = "in"))


# THIS IS THE OUT DEGREE SUMMARY
summary(degree(gd,mode = "out"))


# THIS IS THE WEIGHTED NUMBER OF THE ABOVE vdcENCIES, SO THE NUMBER OF SHARED VDC
# TARGETS IS ACCOUNTED FOR BETWEEN EACH PAIR OF vdcENCIES
summary(graph.strength(gd))


# AGAIN, THE INWARD WEIGHTED DEGREE
summary(graph.strength(gd,mode = "in"))


# AGAIN, THE OUTWARD WEIGHTED DEGREE
summary(graph.strength(gd,mode = "out"))


# PLOT THE NUMBER OF DISTINCT VDC-VDC CONNECTIONS
plot(sort(degree(gd)),
     col = adjustcolor(rgb(1,0,1,1)),
     pch = 19,
     xaxt = "n",
     yaxt = "n",
     ann = FALSE)
par(new=T)
plot(sort(degree(gd,mode = "in")),
     col = "green",
     pch = 19,
     xaxt = "n",
     yaxt = "n",
     ann = FALSE)
par(new=T)
plot(sort(degree(gd,mode = "out")),
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
deg1 <- hist(degree(gd), breaks = 20)$counts
deg1n <- 100*deg1/sum(deg1)
deg2 <- hist(degree(gd,mode = "in"), breaks = 20)$counts
deg2n <- 100*deg2/sum(deg2)
deg3 <- hist(degree(gd,mode = "out"), breaks = 20)$counts
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
plot(sort(graph.strength(gd)),
     col = adjustcolor(rgb(1,0,1,1)),
     pch = 19,
     xaxt = "n",
     yaxt = "n",
     ann = FALSE)
par(new=T)
plot(sort(graph.strength(gd,mode = "in")),
     col = "green",
     pch = 19,
     xaxt = "n",
     yaxt = "n",
     ann = FALSE)
par(new=T)
plot(sort(graph.strength(gd,mode = "out")),
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
deg1 <- hist(graph.strength(gd), breaks = 20)$counts
deg1n <- 100*deg1/sum(deg1)
deg2 <- hist(graph.strength(gd,mode = "in"), breaks = 20)$counts
deg2n <- 100*deg2/sum(deg2)
deg3 <- hist(graph.strength(gd,mode = "out"), breaks = 20)$counts
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
100*graph.density(gd)


# CLUSTERS ARE CONNECTED COMPONENTS
clusters(gd)$no


# SORTED CLUSTERS BY SIZE
sort(clusters(gd)$csize,decreasing = TRUE)


# GLOBAL CLUSTERING COEFFICIENT (TRANSITIVITY) IS THE RATIO OF TRIANGLES AND CONNECTED TRIPLES
transitivity(gd)
cut75 <- quantile(as.vector(dtm[dtm>0]),0.75)
gd_f <- filter(cutoff = cut75,
               edge_matrix = dtm,
               vertex_colors = V(gd)$color,
               vertex_names = V(gd)$name)
gd_f <- as.undirected(gd_f)
transitivity(gd_f)


# RELATIVE MAXIMAL CLUSTER SIZE (AS % OF NUMBER OF NODES) 
max(clusters(gd)$csize)/vcount(gd)


# RELATIVE NUMBER OF ISOLATED NODES (AS % OF NUMBER OF NODES)  
sum(degree(gd) == 0)/vcount(gd)


# PATH DISTRIBUTION: This shows the different lengths of shortest paths (geodesics) in our network. 
sh <- shortest.paths(graph = gd, 
                     mode = "all",
                     weights = E(gd)$weight)
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
# CENTRALITY ANALYSIS OF THE NETWORK
#
#
#
#
#
# BETWEENNESS CENTRALITY: THE NUMBER OF GEODESICS GOING THROUGH A NODE
#
#
gd <- drop_isolated(graph = gd,
                    vertex_colors = V(gd)$color,
                    vertex_names = V(gd)$name)
gd <- drop_loops(graph = gd,
                 vertex_colors = V(gd)$color,
                 vertex_names = V(gd)$name)
bc <- betweenness(graph = gd,
                  v = V(gd),
                  directed = TRUE,
                  weights = E(gd)$weight) 
plot(sort(bc, decreasing = TRUE),
     col = adjustcolor(rgb(1,0,1,1)), 
     xlab = "Node Index", 
     ylab = "Betweenness Centrality", 
     main = "Sorted VDC Network Betweenness Centrality Values", 
     pch = 19)
histP2(bc,
       breaks = 50,
       col = adjustcolor(rgb(1,0,1,1)),
       xlab = "Betweenness Centrality Values",
       main = "VDC Network Betweenness Centrality Distribution")


# PLOT HEAT MAP ON VERTICES ACCORDING TO BETWEENNESS CENTRALITY
bc_int <- as.integer(round(bc,0))/5
for (k in 1:length(bc_int)){
  V(gd)$color[k] <- rev(heat.colors(1+as.integer(max(bc_int))))[as.integer(bc_int[k])+1]
}
gd_coords <- koords[which(vdc %in% V(gd)$name),]
plot(gd,
     layout = layout.fruchterman.reingold(gd, niter = 200, area = 2000*vcount(gd)),
     vertex.color = V(gd)$color,
     vertex.size = 9,
     vertex.label = V(gd)$name, 
     vertex.label.color = "black",
     vertex.label.font = 1, 
     vertex.label.cex = 0.85, 
     edge.width = 0.2*sqrt(E(gd)$weight),
     edge.arrow.size = 0.5,
     edge.curved = TRUE,
     edge.color = gray.colors(1),
     main = "Weighted Nepal Displacement Network Flow (VDC Level)")
legend("topleft",
       c("Highest Betweenness Centrality","Lowest Betweenness Centrality"),
       fill = c("red","White"),
       bty = "n")
plot(gd,
     layout = gd_coords,
     vertex.color = V(gd)$color,
     vertex.size = 4,
     vertex.label = NA, 
     vertex.label.color = "black",
     vertex.label.font = 1, 
     vertex.label.cex = 0.75, 
     edge.width = 0.2*sqrt(E(gd)$weight),
     edge.arrow.size = 0.5,
     edge.curved = TRUE,
     edge.color = gray.colors(1),
     main = "Weighted Nepal Displacement Network Flow (VDC Level)")
legend("top",
       c("Highest Betweenness Centrality","Lowest Betweenness Centrality"),
       fill = c("red","White"),
       bty = "n")


# FIND THE TOP 10% BETWEENNES NODES
top_bc <- bc[which(bc > quantile(bc,0.9))]
top_bc


# FIND THE TOP 5% BETWENNES NODES
top_bc <- bc[which(bc > quantile(bc,0.95))]
top_bc


# FILTRATION AND BETWEENNESS CENTRALITY
gd <- graph.adjacency(dtm,
                      mode = "directed",
                      weighted = TRUE)
V(gd)$color <- rep("SkyBlue2",length(vdc))
for (k in 1:length(vdc)){
  o_count <- length(which(dt_data$idp_origin_vdc == vdc[k]))+
    length(which(dt_data$idp2_origin_vdc == vdc[k]))
  d_count <- length(which(dt_data$vdc == vdc[k]))
  if(o_count>d_count){
    V(gd)$color[k] <- "green"
  } 
}
V(gd)$name <- vdc


# BETWEENNESS CENTRALITY AND FILTER AT CUTOFF = 25%
cut25 <- quantile(as.vector(dtm[dtm>0]),0.25)
gd_f <- filter(cutoff = cut25,
               edge_matrix = dtm,
               vertex_colors = V(gd)$color,
               vertex_names = V(gd)$name)
gd_f <- drop_isolated(graph = gd_f,
                    vertex_colors = V(gd_f)$color,
                    vertex_names = V(gd_f)$name)
gd_f <- drop_loops(graph = gd_f,
                 vertex_colors = V(gd_f)$color,
                 vertex_names = V(gd_f)$name)
bc <- betweenness(graph = gd_f,
                  v = V(gd_f),
                  directed = TRUE,
                  weights = E(gd_f)$weight) 
bc_int <- as.integer(round(bc,0))
for (k in 1:length(bc_int)){
  V(gd_f)$color[k] <- rev(heat.colors(1+as.integer(max(bc_int))))[as.integer(bc_int[k])+1]
}
gd_f_coords <- koords[which(V(gd)$name %in% V(gd_f)$name),]
plot(gd_f,
     layout = layout.fruchterman.reingold(gd_f, niter = 200, area = 2000*vcount(gd_f)),
     vertex.color = V(gd_f)$color,
     vertex.size = 8,
     vertex.label = V(gd_f)$name, 
     vertex.label.color = "black",
     vertex.label.font = 1, 
     vertex.label.cex = 0.9, 
     edge.width = 0.4*sqrt(E(gd_f)$weight),
     edge.curved = TRUE,
     edge.arrow.size = 0.75,
     edge.color = gray.colors(1),
     main = "Filtered (Threshold = 25%) Weighted Abstract Nepal Displacement Network")
legend("topleft",
       c("Highest Betweenness Centrality","Lowest Betweenness Centrality"),
       fill = c("red","White"),
       bty = "n")
plot(gd_f,
     layout = gd_f_coords,
     vertex.color = V(gd_f)$color,
     vertex.size = 4,
     vertex.label = NA, 
     vertex.label.color = "black",
     vertex.label.font = 1, 
     vertex.label.cex = 0.75, 
     edge.width = 0.2*sqrt(E(gd_f)$weight),
     edge.arrow.size = 0.5,
     edge.curved = TRUE,
     edge.color = gray.colors(1),
     main = "Filtered (Threshold = 25%) Weighted Nepal Displacement Geo-Network")
legend("top",
       c("Highest Betweenness Centrality","Lowest Betweenness Centrality"),
       fill = c("red","White"),
       bty = "n")


# BETWEENNESS CENTRALITY AND FILTER AT CUTOFF = 50%
cut50 <- quantile(as.vector(dtm[dtm>0]),0.50)
gd_f <- filter(cutoff = cut50,
               edge_matrix = dtm,
               vertex_colors = V(gd)$color,
               vertex_names = V(gd)$name)
gd_f <- drop_isolated(graph = gd_f,
                      vertex_colors = V(gd_f)$color,
                      vertex_names = V(gd_f)$name)
gd_f <- drop_loops(graph = gd_f,
                   vertex_colors = V(gd_f)$color,
                   vertex_names = V(gd_f)$name)
bc <- betweenness(graph = gd_f,
                  v = V(gd_f),
                  directed = TRUE,
                  weights = E(gd_f)$weight) 
bc_int <- as.integer(round(bc,0))
for (k in 1:length(bc_int)){
  V(gd_f)$color[k] <- rev(heat.colors(1+as.integer(max(bc_int))))[as.integer(bc_int[k])+1]
}
gd_f_coords <- koords[which(V(gd)$name %in% V(gd_f)$name),]
plot(gd_f,
     layout = layout.fruchterman.reingold(gd_f, niter = 200, area = 2000*vcount(gd_f)),
     vertex.color = V(gd_f)$color,
     vertex.size = 8,
     vertex.label = V(gd_f)$name, 
     vertex.label.color = "black",
     vertex.label.font = 1, 
     vertex.label.cex = 0.9, 
     edge.width = 0.4*sqrt(E(gd_f)$weight),
     edge.curved = TRUE,
     edge.arrow.size = 0.75,
     edge.color = gray.colors(1),
     main = "Filtered (Threshold = 50%) Weighted Abstract Nepal Displacement Network")
legend("topleft",
       c("Highest Betweenness Centrality","Lowest Betweenness Centrality"),
       fill = c("red","White"),
       bty = "n")
plot(gd_f,
     layout = gd_f_coords,
     vertex.color = V(gd_f)$color,
     vertex.size = 6,
     vertex.label = V(gd_f)$name, 
     vertex.label.color = "black",
     vertex.label.font = 1, 
     vertex.label.cex = 0.9, 
     edge.width = 0.2*sqrt(E(gd_f)$weight),
     edge.arrow.size = 0.5,
     edge.curved = TRUE,
     edge.color = gray.colors(1),
     main = "Filtered (Threshold = 50%) Weighted Nepal Displacement Geo-Network")
legend("top",
       c("Highest Betweenness Centrality","Lowest Betweenness Centrality"),
       fill = c("red","White"),
       bty = "n")


# BETWEENNESS CENTRALITY AND FILTER AT CUTOFF = 75%
cut75 <- quantile(as.vector(dtm[dtm>0]),0.75)
gd_f <- filter(cutoff = cut75,
               edge_matrix = dtm,
               vertex_colors = V(gd)$color,
               vertex_names = V(gd)$name)
gd_f <- drop_isolated(graph = gd_f,
                      vertex_colors = V(gd_f)$color,
                      vertex_names = V(gd_f)$name)
gd_f <- drop_loops(graph = gd_f,
                   vertex_colors = V(gd_f)$color,
                   vertex_names = V(gd_f)$name)
bc <- betweenness(graph = gd_f,
                  v = V(gd_f),
                  directed = TRUE,
                  weights = E(gd_f)$weight) 
bc_int <- as.integer(round(bc,0))
for (k in 1:length(bc_int)){
  V(gd_f)$color[k] <- rev(heat.colors(1+as.integer(max(bc_int))))[as.integer(bc_int[k])+1]
}
gd_f_coords <- koords[which(V(gd)$name %in% V(gd_f)$name),]
plot(gd_f,
     layout = layout.fruchterman.reingold(gd_f, niter = 200, area = 2000*vcount(gd_f)),
     vertex.color = V(gd_f)$color,
     vertex.size = 8,
     vertex.label = V(gd_f)$name, 
     vertex.label.color = "black",
     vertex.label.font = 1, 
     vertex.label.cex = 0.9, 
     edge.width = 0.4*sqrt(E(gd_f)$weight),
     edge.curved = TRUE,
     edge.arrow.size = 0.75,
     edge.color = gray.colors(1),
     main = "Filtered (Threshold = 75%) Weighted Abstract Nepal Displacement Network")
legend("topleft",
       c("Highest Betweenness Centrality","Lowest Betweenness Centrality"),
       fill = c("red","White"),
       bty = "n")
plot(gd_f,
     layout = gd_f_coords,
     vertex.color = V(gd_f)$color,
     vertex.size = 8,
     vertex.label = V(gd_f)$name, 
     vertex.label.color = "black",
     vertex.label.font = 1, 
     vertex.label.cex = 1, 
     edge.width = 0.3*sqrt(E(gd_f)$weight),
     edge.arrow.size = 0.75,
     edge.curved = TRUE,
     edge.color = gray.colors(1),
     main = "Filtered (Threshold = 75%) Weighted Nepal Displacement Geo-Network")
legend("top",
       c("Highest Betweenness Centrality","Lowest Betweenness Centrality"),
       fill = c("red","White"),
       bty = "n")

#
#
# EDGE-BETWEENNESS CENTRALITY: THE NUMBER OF GEODESICS GOING THROUGH AN EDGE
#
#
gd <- graph.adjacency(dtm,
                      mode = "directed",
                      weighted = TRUE)
V(gd)$color <- rep("SkyBlue2",length(vdc))
for (k in 1:length(vdc)){
  o_count <- length(which(dt_data$idp_origin_vdc == vdc[k]))+
    length(which(dt_data$idp2_origin_vdc == vdc[k]))
  d_count <- length(which(dt_data$vdc == vdc[k]))
  if(o_count>d_count){
    V(gd)$color[k] <- "green"
  } 
}
V(gd)$name <- vdc
gd <- drop_isolated(graph = gd,
                    vertex_colors = V(gd)$color,
                    vertex_names = V(gd)$name)
gd <- drop_loops(graph = gd,
                 vertex_colors = V(gd)$color,
                 vertex_names = V(gd)$name)
ec <- edge.betweenness(graph = gd,
                       e = E(gd),
                       directed = TRUE,
                       weights = E(gd)$weight)
plot(sort(ec, decreasing = TRUE),
     col = adjustcolor(rgb(1,0,1,1)), 
     xlab = "Node Index", 
     ylab = "Edge Betweenness Centrality", 
     main = "Sorted VDC Network Edge Betweenness Centrality Values", 
     pch = 19)
histP2(ec,
       breaks = 50,
       col = adjustcolor(rgb(1,0,1,1)),
       xlab = "Edge Betweenness Centrality Values",
       main = "VDC Network Edge Betweenness Centrality Distribution")


# PLOT HEAT MAP ON VERTICES ACCORDING TO EDGE BETWEENNESS CENTRALITY
for (k in 1:length(ec)){
  E(gd)$color[k] <- rev(heat.colors(1+as.integer(max(ec))))[as.integer(ec[k])+1]
}
gd_coords <- koords[which(vdc %in% V(gd)$name),]
plot(gd,
     layout = layout.fruchterman.reingold(gd, niter = 200, area = 2000*vcount(gd)),
     vertex.color = V(gd)$color,
     vertex.size = 9,
     vertex.label = V(gd)$name, 
     vertex.label.color = "black",
     vertex.label.font = 1, 
     vertex.label.cex = 0.85, 
     edge.width = 0.2*sqrt(E(gd)$weight),
     edge.arrow.size = 0.5,
     edge.curved = TRUE,
     edge.color = E(gd)$color,
     main = "Weighted Nepal Displacement Network Flow (VDC Level)")
legend("topleft",
       c("Highest Edge Betweenness Centrality","Lowest Edge Betweenness Centrality"),
       fill = c("red","White"),
       bty = "n")
plot(gd,
     layout = gd_coords,
     vertex.color = V(gd)$color,
     vertex.size = 4,
     vertex.label = NA, 
     vertex.label.color = "black",
     vertex.label.font = 1, 
     vertex.label.cex = 0.75, 
     edge.width = 0.2*sqrt(E(gd)$weight),
     edge.arrow.size = 0.5,
     edge.curved = TRUE,
     edge.color = E(gd)$color,
     main = "Weighted Nepal Displacement Geo-Network Flow (VDC Level)")
legend("topright",
       c("Highest Edge Betweenness Centrality","Lowest Edge Betweenness Centrality"),
       fill = c("red","White"),
       bty = "n")


# FILTRATION AND EDGE BETWEENNESS CENTRALITY
gd <- graph.adjacency(dtm,
                      mode = "directed",
                      weighted = TRUE)
V(gd)$color <- rep("SkyBlue2",length(vdc))
for (k in 1:length(vdc)){
  o_count <- length(which(dt_data$idp_origin_vdc == vdc[k]))+
    length(which(dt_data$idp2_origin_vdc == vdc[k]))
  d_count <- length(which(dt_data$vdc == vdc[k]))
  if(o_count>d_count){
    V(gd)$color[k] <- "green"
  } 
}
V(gd)$name <- vdc


# EDGE BETWEENNESS CENTRALITY AND FILTER AT CUTOFF = 25%
cut25 <- quantile(as.vector(dtm[dtm>0]),0.25)
gd_f <- filter(cutoff = cut25,
               edge_matrix = dtm,
               vertex_colors = V(gd)$color,
               vertex_names = V(gd)$name)
gd_f <- drop_isolated(graph = gd_f,
                      vertex_colors = V(gd_f)$color,
                      vertex_names = V(gd_f)$name)
gd_f <- drop_loops(graph = gd_f,
                   vertex_colors = V(gd_f)$color,
                   vertex_names = V(gd_f)$name)
ec <- edge.betweenness(graph = gd_f,
                       e = E(gd_f),
                       directed = TRUE,
                       weights = E(gd_f)$weight)
for (k in 1:length(ec)){
  E(gd_f)$color[k] <- rev(heat.colors(1+as.integer(max(ec))))[as.integer(ec[k])+1]
}
gd_f_coords <- koords[which(V(gd)$name %in% V(gd_f)$name),]
plot(gd_f,
     layout = layout.fruchterman.reingold(gd_f, niter = 200, area = 2000*vcount(gd_f)),
     vertex.color = V(gd_f)$color,
     vertex.size = 8,
     vertex.label = V(gd_f)$name, 
     vertex.label.color = "black",
     vertex.label.font = 1, 
     vertex.label.cex = 0.9, 
     edge.width = 0.2*sqrt(E(gd_f)$weight),
     edge.curved = TRUE,
     edge.arrow.size = 0.75,
     edge.color = E(gd_f)$color,
     main = "Filtered (Threshold = 25%) Weighted Abstract Nepal Displacement Network")
legend("topleft",
       c("Highest Betweenness Centrality","Lowest Betweenness Centrality"),
       fill = c("red","White"),
       bty = "n")
plot(gd_f,
     layout = gd_f_coords,
     vertex.color = V(gd_f)$color,
     vertex.size = 4,
     vertex.label = V(gd_f)$name, 
     vertex.label.color = "black",
     vertex.label.font = 1, 
     vertex.label.cex = 0.75, 
     edge.width = 0.2*sqrt(E(gd_f)$weight),
     edge.arrow.size = 0.5,
     edge.curved = TRUE,
     edge.color = E(gd_f)$color,
     main = "Filtered (Threshold = 25%) Weighted Nepal Displacement Geo-Network")
legend("topright",
       c("Highest Betweenness Centrality","Lowest Betweenness Centrality"),
       fill = c("red","White"),
       bty = "n")


# EDGE BETWEENNESS CENTRALITY AND FILTER AT CUTOFF = 50%
cut50 <- quantile(as.vector(dtm[dtm>0]),0.50)
gd_f <- filter(cutoff = cut50,
               edge_matrix = dtm,
               vertex_colors = V(gd)$color,
               vertex_names = V(gd)$name)
gd_f <- drop_isolated(graph = gd_f,
                      vertex_colors = V(gd_f)$color,
                      vertex_names = V(gd_f)$name)
gd_f <- drop_loops(graph = gd_f,
                   vertex_colors = V(gd_f)$color,
                   vertex_names = V(gd_f)$name)
ec <- edge.betweenness(graph = gd_f,
                       e = E(gd_f),
                       directed = TRUE,
                       weights = E(gd_f)$weight)
for (k in 1:length(ec)){
  E(gd_f)$color[k] <- rev(heat.colors(1+as.integer(max(ec))))[as.integer(ec[k])+1]
}
gd_f_coords <- koords[which(V(gd)$name %in% V(gd_f)$name),]
plot(gd_f,
     layout = layout.fruchterman.reingold(gd_f, niter = 200, area = 2000*vcount(gd_f)),
     vertex.color = V(gd_f)$color,
     vertex.size = 8,
     vertex.label = V(gd_f)$name, 
     vertex.label.color = "black",
     vertex.label.font = 1, 
     vertex.label.cex = 0.9, 
     edge.width = 0.2*sqrt(E(gd_f)$weight),
     edge.curved = TRUE,
     edge.arrow.size = 0.75,
     edge.color = E(gd_f)$color,
     main = "Filtered (Threshold = 50%) Weighted Abstract Nepal Displacement Network")
legend("topleft",
       c("Highest Betweenness Centrality","Lowest Betweenness Centrality"),
       fill = c("red","White"),
       bty = "n")
plot(gd_f,
     layout = gd_f_coords,
     vertex.color = V(gd_f)$color,
     vertex.size = 5,
     vertex.label = V(gd_f)$name, 
     vertex.label.color = "black",
     vertex.label.font = 1, 
     vertex.label.cex = 0.85, 
     edge.width = 0.2*sqrt(E(gd_f)$weight),
     edge.arrow.size = 0.5,
     edge.curved = TRUE,
     edge.color = E(gd_f)$color,
     main = "Filtered (Threshold = 50%) Weighted Nepal Displacement Geo-Network")
legend("topright",
       c("Highest Betweenness Centrality","Lowest Betweenness Centrality"),
       fill = c("red","White"),
       bty = "n")


# EDGE BETWEENNESS CENTRALITY AND FILTER AT CUTOFF = 75%
cut75 <- quantile(as.vector(dtm[dtm>0]),0.75)
gd_f <- filter(cutoff = cut75,
               edge_matrix = dtm,
               vertex_colors = V(gd)$color,
               vertex_names = V(gd)$name)
gd_f <- drop_isolated(graph = gd_f,
                      vertex_colors = V(gd_f)$color,
                      vertex_names = V(gd_f)$name)
gd_f <- drop_loops(graph = gd_f,
                   vertex_colors = V(gd_f)$color,
                   vertex_names = V(gd_f)$name)
ec <- edge.betweenness(graph = gd_f,
                       e = E(gd_f),
                       directed = TRUE,
                       weights = E(gd_f)$weight)
for (k in 1:length(ec)){
  E(gd_f)$color[k] <- rev(heat.colors(1+as.integer(max(ec))))[as.integer(ec[k])+1]
}
gd_f_coords <- koords[which(V(gd)$name %in% V(gd_f)$name),]
plot(gd_f,
     layout = layout.fruchterman.reingold(gd_f, niter = 200, area = 2000*vcount(gd_f)),
     vertex.color = V(gd_f)$color,
     vertex.size = 8,
     vertex.label = V(gd_f)$name, 
     vertex.label.color = "black",
     vertex.label.font = 1, 
     vertex.label.cex = 0.9, 
     edge.width = 0.2*sqrt(E(gd_f)$weight),
     edge.curved = TRUE,
     edge.arrow.size = 0.75,
     edge.color = E(gd_f)$color,
     main = "Filtered (Threshold = 75%) Weighted Abstract Nepal Displacement Network")
legend("topleft",
       c("Highest Edge Betweenness Centrality","Lowest Edge Betweenness Centrality"),
       fill = c("red","White"),
       bty = "n")
plot(gd_f,
     layout = gd_f_coords,
     vertex.color = V(gd_f)$color,
     vertex.size = 6,
     vertex.label = V(gd_f)$name, 
     vertex.label.color = "black",
     vertex.label.font = 1, 
     vertex.label.cex = 0.95, 
     edge.width = 0.25*sqrt(E(gd_f)$weight),
     edge.arrow.size = 0.5,
     edge.curved = TRUE,
     edge.color = E(gd_f)$color,
     main = "Filtered (Threshold = 75%) Weighted Nepal Displacement Geo-Network")
legend("topright",
       c("Highest Edge Betweenness Centrality","Lowest Edge Betweenness Centrality"),
       fill = c("red","White"),
       bty = "n")

#
#
# CLOSENES CENTRALITY: THIS MEASURE TAKES INTO ACCOUNT THE DISTRIBUTION OF DISTANCES
# TO OTHER NODES FROM A GIVEN NODE. IT IS DEFINED AS THE RCIPROCAL OF THE FARNESS OF A NODE.
# FARNESS OF A NODE IS DEFINED AS THE SUM OF ITS DISTANCES TO ALL OTHER NODES. CLOSENESS CAN BE
# REGARDED AS A MEASURE OF HOW LONG IT WILL TAKE TO SPREAD INFORMATION FROM A NODE TO ALL OTHER
# NODES.
#
#
gd <- graph.adjacency(dtm,
                      mode = "directed",
                      weighted = TRUE)
V(gd)$color <- rep("SkyBlue2",length(vdc))
for (k in 1:length(vdc)){
  o_count <- length(which(dt_data$idp_origin_vdc == vdc[k]))+
    length(which(dt_data$idp2_origin_vdc == vdc[k]))
  d_count <- length(which(dt_data$vdc == vdc[k]))
  if(o_count>d_count){
    V(gd)$color[k] <- "green"
  } 
}
V(gd)$name <- vdc
gd <- drop_isolated(graph = gd,
                    vertex_colors = V(gd)$color,
                    vertex_names = V(gd)$name)
gd <- drop_loops(graph = gd,
                 vertex_colors = V(gd)$color,
                 vertex_names = V(gd)$name)
# cl <- clusters(gd)
# gd <- induced.subgraph(gd, which(cl$membership == which.max(cl$csize)))
cc <- closeness(graph = gd,
                vids = V(gd),
                weights = E(gd)$weight)
cc <- 100*(cc/max(cc))^9
plot(sort(cc, decreasing = TRUE),
     col = adjustcolor(rgb(1,0,1,1)), 
     xlab = "Node Index", 
     ylab = "Closeness Centrality", 
     main = "Sorted VDC Network Closeness Centrality Values", 
     pch = 19)
histP2(cc,
       breaks = 50,
       col = adjustcolor(rgb(1,0,1,1)),
       xlab = "Closeness Centrality Values",
       main = "VDC Network Closeness Centrality Distribution")


# PLOT HEAT MAP ON VERTICES ACCORDING TO CLOSENESS CENTRALITY
for (k in 1:length(cc)){
  V(gd)$color[k] <- rev(heat.colors(1+as.integer(max(cc))))[as.integer(cc[k])+1]
}
gd_coords <- koords[which(vdc %in% V(gd)$name),]
plot(gd,
     layout = layout.fruchterman.reingold(gd, 
                                          niter = 200, 
                                          area = 2000*vcount(gd)),
     vertex.color = V(gd)$color,
     vertex.size = 9,
     vertex.label = V(gd)$name, 
     vertex.label.color = "black",
     vertex.label.font = 1, 
     vertex.label.cex = 0.85, 
     edge.width = 0.2*sqrt(E(gd)$weight),
     edge.arrow.size = 0.5,
     edge.curved = TRUE,
     edge.color = gray.colors(1),
     main = "Weighted Nepal Displacement Network Flow (VDC Level)")
legend("topleft",
       c("Highest Closeness Centrality","Lowest Closeness Centrality"),
       fill = c("red","White"),
       bty = "n")
plot(gd,
     layout = gd_coords,
     vertex.color = V(gd)$color,
     vertex.size = 4,
     vertex.label = NA, 
     vertex.label.color = "black",
     vertex.label.font = 1, 
     vertex.label.cex = 0.75, 
     edge.width = 0.2*sqrt(E(gd)$weight),
     edge.arrow.size = 0.5,
     edge.curved = TRUE,
     edge.color = gray.colors(1),
     main = "Weighted Nepal Displacement Geo-Network Flow (VDC Level)")
legend("topright",
       c("Highest Closeness Centrality","Lowest Closeness Centrality"),
       fill = c("red","White"),
       bty = "n")


# FILTRATION AND CLOSENESS CENTRALITY
gd <- graph.adjacency(dtm,
                      mode = "directed",
                      weighted = TRUE)
V(gd)$color <- rep("SkyBlue2",length(vdc))
for (k in 1:length(vdc)){
  o_count <- length(which(dt_data$idp_origin_vdc == vdc[k]))+
    length(which(dt_data$idp2_origin_vdc == vdc[k]))
  d_count <- length(which(dt_data$vdc == vdc[k]))
  if(o_count>d_count){
    V(gd)$color[k] <- "green"
  } 
}
V(gd)$name <- vdc


# CLOSENESS CENTRALITY AND FILTER AT CUTOFF = 25%
cut25 <- quantile(as.vector(dtm[dtm>0]),0.25)
gd_f <- filter(cutoff = cut25,
               edge_matrix = dtm,
               vertex_colors = V(gd)$color,
               vertex_names = V(gd)$name)
gd_f <- drop_isolated(graph = gd_f,
                      vertex_colors = V(gd_f)$color,
                      vertex_names = V(gd_f)$name)
gd_f <- drop_loops(graph = gd_f,
                   vertex_colors = V(gd_f)$color,
                   vertex_names = V(gd_f)$name)
cc <- closeness(graph = gd_f,
                vids = V(gd_f),
                weights = E(gd_f)$weight)
cc <- 100*(cc/max(cc))^9
for (k in 1:length(cc)){
  V(gd_f)$color[k] <- rev(heat.colors(1+as.integer(max(cc))))[as.integer(cc[k])+1]
}
gd_f_coords <- koords[which(V(gd)$name %in% V(gd_f)$name),]
plot(gd_f,
     layout = layout.fruchterman.reingold(gd_f,
                                          niter = 200,
                                          area = 2000*vcount(gd_f)),
     vertex.color = V(gd_f)$color,
     vertex.size = 8,
     vertex.label = V(gd_f)$name, 
     vertex.label.color = "black",
     vertex.label.font = 1, 
     vertex.label.cex = 0.9, 
     edge.width = 0.2*sqrt(E(gd_f)$weight),
     edge.curved = TRUE,
     edge.arrow.size = 0.65,
     edge.color = grey.colors(1),
     main = "Filtered (Threshold = 25%) Weighted Abstract Nepal Displacement Network")
legend("topleft",
       c("Highest Closeness Centrality","Lowest Closeness Centrality"),
       fill = c("red","White"),
       bty = "n")
plot(gd_f,
     layout = gd_f_coords,
     vertex.color = V(gd_f)$color,
     vertex.size = 6,
     vertex.label = V(gd_f)$name, 
     vertex.label.color = "black",
     vertex.label.font = 1, 
     vertex.label.cex = 0.75, 
     edge.width = 0.2*sqrt(E(gd_f)$weight),
     edge.arrow.size = 0.5,
     edge.curved = TRUE,
     edge.color = grey.colors(1),
     main = "Filtered (Threshold = 25%) Weighted Nepal Displacement Geo-Network")
legend("topright",
       c("Highest Closeness Centrality","Lowest Closeness Centrality"),
       fill = c("red","White"),
       bty = "n")


# CLOSENESS CENTRALITY AND FILTER AT CUTOFF = 50%
cut50 <- quantile(as.vector(dtm[dtm>0]),0.50)
gd_f <- filter(cutoff = cut50,
               edge_matrix = dtm,
               vertex_colors = V(gd)$color,
               vertex_names = V(gd)$name)
gd_f <- drop_isolated(graph = gd_f,
                      vertex_colors = V(gd_f)$color,
                      vertex_names = V(gd_f)$name)
gd_f <- drop_loops(graph = gd_f,
                   vertex_colors = V(gd_f)$color,
                   vertex_names = V(gd_f)$name)
cc <- closeness(graph = gd_f,
                vids = V(gd_f),
                weights = E(gd_f)$weight)
cc <- 100*(cc/max(cc))^9
for (k in 1:length(cc)){
  V(gd_f)$color[k] <- rev(heat.colors(1+as.integer(max(cc))))[as.integer(cc[k])+1]
}
gd_f_coords <- koords[which(V(gd)$name %in% V(gd_f)$name),]
plot(gd_f,
     layout = layout.fruchterman.reingold(gd_f,
                                          niter = 200,
                                          area = 2000*vcount(gd_f)),
     vertex.color = V(gd_f)$color,
     vertex.size = 8,
     vertex.label = V(gd_f)$name, 
     vertex.label.color = "black",
     vertex.label.font = 1, 
     vertex.label.cex = 0.9, 
     edge.width = 0.2*sqrt(E(gd_f)$weight),
     edge.curved = TRUE,
     edge.arrow.size = 0.65,
     edge.color = grey.colors(1),
     main = "Filtered (Threshold = 50%) Weighted Abstract Nepal Displacement Network")
legend("topleft",
       c("Highest Closeness Centrality","Lowest Closeness Centrality"),
       fill = c("red","White"),
       bty = "n")
plot(gd_f,
     layout = gd_f_coords,
     vertex.color = V(gd_f)$color,
     vertex.size = 6,
     vertex.label = V(gd_f)$name, 
     vertex.label.color = "black",
     vertex.label.font = 1, 
     vertex.label.cex = 0.75, 
     edge.width = 0.2*sqrt(E(gd_f)$weight),
     edge.arrow.size = 0.5,
     edge.curved = TRUE,
     edge.color = grey.colors(1),
     main = "Filtered (Threshold = 50%) Weighted Nepal Displacement Geo-Network")
legend("topright",
       c("Highest Closeness Centrality","Lowest Closeness Centrality"),
       fill = c("red","White"),
       bty = "n")


# EDGE BETWEENNESS CENTRALITY AND FILTER AT CUTOFF = 75%
cut75 <- quantile(as.vector(dtm[dtm>0]),0.75)
gd_f <- filter(cutoff = cut75,
               edge_matrix = dtm,
               vertex_colors = V(gd)$color,
               vertex_names = V(gd)$name)
gd_f <- drop_isolated(graph = gd_f,
                      vertex_colors = V(gd_f)$color,
                      vertex_names = V(gd_f)$name)
gd_f <- drop_loops(graph = gd_f,
                   vertex_colors = V(gd_f)$color,
                   vertex_names = V(gd_f)$name)
cc <- closeness(graph = gd_f,
                vids = V(gd_f),
                weights = E(gd_f)$weight)
cc <- 100*(cc/max(cc))^9
for (k in 1:length(cc)){
  V(gd_f)$color[k] <- rev(heat.colors(1+as.integer(max(cc))))[as.integer(cc[k])+1]
}
gd_f_coords <- koords[which(V(gd)$name %in% V(gd_f)$name),]
plot(gd_f,
     layout = layout.fruchterman.reingold(gd_f,
                                          niter = 200,
                                          area = 2000*vcount(gd_f)),
     vertex.color = V(gd_f)$color,
     vertex.size = 12,
     vertex.label = V(gd_f)$name, 
     vertex.label.color = "black",
     vertex.label.font = 1, 
     vertex.label.cex = 1.3, 
     edge.width = 0.3*sqrt(E(gd_f)$weight),
     edge.curved = TRUE,
     edge.arrow.size = 1,
     edge.color = grey.colors(1),
     main = "Filtered (Threshold = 75%) Weighted Abstract Nepal Displacement Network")
legend("topleft",
       c("Highest Closeness Centrality","Lowest Closeness Centrality"),
       fill = c("red","White"),
       bty = "n")
plot(gd_f,
     layout = gd_f_coords,
     vertex.color = V(gd_f)$color,
     vertex.size = 8,
     vertex.label = V(gd_f)$name, 
     vertex.label.color = "black",
     vertex.label.font = 1, 
     vertex.label.cex = 1, 
     edge.width = 0.3*sqrt(E(gd_f)$weight),
     edge.arrow.size = 0.75,
     edge.curved = TRUE,
     edge.color = grey.colors(1),
     main = "Filtered (Threshold = 75%) Weighted Nepal Displacement Geo-Network")
legend("topright",
       c("Highest Closeness Centrality","Lowest Closeness Centrality"),
       fill = c("red","White"),
       bty = "n")
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
gd <- graph.adjacency(dtm,
                      mode = "directed",
                      weighted = TRUE)
V(gd)$color <- rep("SkyBlue2",length(vdc))
for (k in 1:length(vdc)){
  o_count <- length(which(dt_data$idp_origin_vdc == vdc[k]))+
    length(which(dt_data$idp2_origin_vdc == vdc[k]))
  d_count <- length(which(dt_data$vdc == vdc[k]))
  if(o_count>d_count){
    V(gd)$color[k] <- "green"
  } 
}
V(gd)$name <- vdc
gd <- drop_isolated(graph = gd,
                    vertex_colors = V(gd)$color,
                    vertex_names = V(gd)$name)
gd <- drop_loops(graph = gd,
                 vertex_colors = V(gd)$color,
                 vertex_names = V(gd)$name)
ec <- evcent(gd)$vector
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
  V(gd)$color[k] <- rev(heat.colors(1+as.integer(max(ec))))[as.integer(ec[k])+1]
}
gd_coords <- koords[which(vdc %in% V(gd)$name),]
plot(gd,
     layout = layout.fruchterman.reingold(gd, 
                                          niter = 200, 
                                          area = 2000*vcount(gd)),
     vertex.color = V(gd)$color,
     vertex.size = 9,
     vertex.label = V(gd)$name, 
     vertex.label.color = "black",
     vertex.label.font = 1, 
     vertex.label.cex = 0.85, 
     edge.width = 0.2*sqrt(E(gd)$weight),
     edge.arrow.size = 0.5,
     edge.curved = TRUE,
     edge.color = gray.colors(1),
     main = "Weighted Nepal Displacement Network Flow (VDC Level)")
legend("topleft",
       c("Highest Eigenvector Centrality","Lowest Eigenvector Centrality"),
       fill = c("red","White"),
       bty = "n")
plot(gd,
     layout = gd_coords,
     vertex.color = V(gd)$color,
     vertex.size = 4,
     vertex.label = NA, 
     vertex.label.color = "black",
     vertex.label.font = 1, 
     vertex.label.cex = 0.75, 
     edge.width = 0.2*sqrt(E(gd)$weight),
     edge.arrow.size = 0.5,
     edge.curved = TRUE,
     edge.color = gray.colors(1),
     main = "Weighted Nepal Displacement Geo-Network Flow (VDC Level)")
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

gd <- graph.adjacency(dtm,
                      mode = "directed",
                      weighted = TRUE)
V(gd)$color <- rep("SkyBlue2",length(vdc))
for (k in 1:length(vdc)){
  o_count <- length(which(dt_data$idp_origin_vdc == vdc[k]))+
    length(which(dt_data$idp2_origin_vdc == vdc[k]))
  d_count <- length(which(dt_data$vdc == vdc[k]))
  if(o_count>d_count){
    V(gd)$color[k] <- "green"
  } 
}
V(gd)$name <- vdc
gd <- drop_isolated(graph = gd,
                    vertex_colors = V(gd)$color,
                    vertex_names = V(gd)$name)
gd <- drop_loops(graph = gd,
                 vertex_colors = V(gd)$color,
                 vertex_names = V(gd)$name)
au <- authority.score(gd)$vector
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
  V(gd)$color[k] <- rev(heat.colors(1+as.integer(max(au))))[as.integer(au[k])+1]
}
gd_coords <- koords[which(vdc %in% V(gd)$name),]
plot(gd,
     layout = layout.fruchterman.reingold(gd, 
                                          niter = 200, 
                                          area = 2000*vcount(gd)),
     vertex.color = V(gd)$color,
     vertex.size = 9,
     vertex.label = V(gd)$name, 
     vertex.label.color = "black",
     vertex.label.font = 1, 
     vertex.label.cex = 0.85, 
     edge.width = 0.2*sqrt(E(gd)$weight),
     edge.arrow.size = 0.5,
     edge.curved = TRUE,
     edge.color = gray.colors(1),
     main = "Weighted Nepal Displacement Network Flow (VDC Level)")
legend("topleft",
       c("Highest Authority Score","Lowest Authority Score"),
       fill = c("red","White"),
       bty = "n")
plot(gd,
     layout = gd_coords,
     vertex.color = V(gd)$color,
     vertex.size = 4,
     vertex.label = NA, 
     vertex.label.color = "black",
     vertex.label.font = 1, 
     vertex.label.cex = 0.75, 
     edge.width = 0.2*sqrt(E(gd)$weight),
     edge.arrow.size = 0.5,
     edge.curved = TRUE,
     edge.color = gray.colors(1),
     main = "Weighted Nepal Displacement Geo-Network Flow (VDC Level)")
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
gd <- graph.adjacency(dtm,
                      mode = "directed",
                      weighted = TRUE)
V(gd)$color <- rep("SkyBlue2",length(vdc))
for (k in 1:length(vdc)){
  o_count <- length(which(dt_data$idp_origin_vdc == vdc[k]))+
    length(which(dt_data$idp2_origin_vdc == vdc[k]))
  d_count <- length(which(dt_data$vdc == vdc[k]))
  if(o_count>d_count){
    V(gd)$color[k] <- "green"
  } 
}
V(gd)$name <- vdc
gd <- drop_isolated(graph = gd,
                    vertex_colors = V(gd)$color,
                    vertex_names = V(gd)$name)
gd <- drop_loops(graph = gd,
                 vertex_colors = V(gd)$color,
                 vertex_names = V(gd)$name)
hb <- hub.score(gd)$vector
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
  V(gd)$color[k] <- rev(heat.colors(1+as.integer(max(hb))))[as.integer(hb[k])+1]
}
gd_coords <- koords[which(vdc %in% V(gd)$name),]
plot(gd,
     layout = layout.fruchterman.reingold(gd, 
                                          niter = 200, 
                                          area = 2000*vcount(gd)),
     vertex.color = V(gd)$color,
     vertex.size = 9,
     vertex.label = V(gd)$name, 
     vertex.label.color = "black",
     vertex.label.font = 1, 
     vertex.label.cex = 0.85, 
     edge.width = 0.2*sqrt(E(gd)$weight),
     edge.arrow.size = 0.5,
     edge.curved = TRUE,
     edge.color = gray.colors(1),
     main = "Weighted Nepal Displacement Network Flow (VDC Level)")
legend("topleft",
       c("Highest Hub Score","Lowest Hub Score"),
       fill = c("red","White"),
       bty = "n")
plot(gd,
     layout = gd_coords,
     vertex.color = V(gd)$color,
     vertex.size = 4,
     vertex.label = NA, 
     vertex.label.color = "black",
     vertex.label.font = 1, 
     vertex.label.cex = 0.75, 
     edge.width = 0.2*sqrt(E(gd)$weight),
     edge.arrow.size = 0.5,
     edge.curved = TRUE,
     edge.color = gray.colors(1),
     main = "Weighted Nepal Displacement Geo-Network Flow (VDC Level)")
legend("topright",
       c("Highest Hub Score","Lowest Hub Score"),
       fill = c("red","White"),
       bty = "n")

#
#
#
# CLUSTERING COEFFICIENT: THIS MEASURES THE PROBABILITY THAT TWO 
# VERTICES ARE CONNECTED
#
# NOTE: SINC ENO CLSOED TRIPLETS ARE PRESENT, WE GET 0 OR NA VALUES ONLY
# THIS IS LIKELY NOT THE RIGHT TOOL TO APPLY HERE aT THE VDC LEVEL
# BUT IF WE OBTAIN MORE GRANUALR DATA, WE WILL BE ABLE TO USE IT
#
gd <- graph.adjacency(dtm,
                      mode = "directed",
                      weighted = TRUE)
V(gd)$color <- rep("SkyBlue2",length(vdc))
for (k in 1:length(vdc)){
  o_count <- length(which(dt_data$idp_origin_vdc == vdc[k]))+
    length(which(dt_data$idp2_origin_vdc == vdc[k]))
  d_count <- length(which(dt_data$vdc == vdc[k]))
  if(o_count>d_count){
    V(gd)$color[k] <- "green"
  } 
}
V(gd)$name <- vdc
gd <- drop_isolated(graph = gd,
                    vertex_colors = V(gd)$color,
                    vertex_names = V(gd)$name)
gd <- drop_loops(graph = gd,
                 vertex_colors = V(gd)$color,
                 vertex_names = V(gd)$name)
tr <- transitivity(graph = gd,type = "local")
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
# for (k in 1:length(hb)){
#   V(gd)$color[k] <- rev(heat.colors(1+as.integer(max(hb))))[as.integer(hb[k])+1]
# }
# gd_coords <- koords[which(vdc %in% V(gd)$name),]
# plot(gd,
#      layout = layout.fruchterman.reingold(gd, 
#                                           niter = 200, 
#                                           area = 2000*vcount(gd)),
#      vertex.color = V(gd)$color,
#      vertex.size = 9,
#      vertex.label = V(gd)$name, 
#      vertex.label.color = "black",
#      vertex.label.font = 1, 
#      vertex.label.cex = 0.85, 
#      edge.width = 0.2*sqrt(E(gd)$weight),
#      edge.arrow.size = 0.5,
#      edge.curved = TRUE,
#      edge.color = gray.colors(1),
#      main = "Weighted Nepal Displacement Network Flow (VDC Level)")
# legend("topleft",
#        c("Highest Hub Score","Lowest Hub Score"),
#        fill = c("red","White"),
#        bty = "n")
# plot(gd,
#      layout = gd_coords,
#      vertex.color = V(gd)$color,
#      vertex.size = 4,
#      vertex.label = NA, 
#      vertex.label.color = "black",
#      vertex.label.font = 1, 
#      vertex.label.cex = 0.75, 
#      edge.width = 0.2*sqrt(E(gd)$weight),
#      edge.arrow.size = 0.5,
#      edge.curved = TRUE,
#      edge.color = gray.colors(1),
#      main = "Weighted Nepal Displacement Geo-Network Flow (VDC Level)")
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
# EDGE-BETWEENNESS COMMUNITY
# The edge betweenness score of an edge measures the number of shortest paths through it, see edge.betweenness for details. 
# The idea of the edge betweenness based community structure detection is that it is likely that edges connecting separate modules 
# have high edge betweenness as all the shortest paths from one module to another must traverse through them. 
# So if we gradually remove the edge with the highest edge betweenness score we will get a hierarchical map, a rooted tree, 
# called a dendrogram of the graph. The leafs of the tree are the individual vertices and the root of the tree represents the whole graph.
# edge.betweenness.community performs this algorithm by calculating the edge betweenness of the graph, 
# removing the edge with the highest edge betweenness score, then recalculating edge betweenness of the edges 
# and again removing the one with the highest score, etc.
#
#
#
gd <- graph.adjacency(dtm,
                      mode = "directed",
                      weighted = TRUE)
V(gd)$color <- rep("SkyBlue2",length(vdc))
for (k in 1:length(vdc)){
  o_count <- length(which(dt_data$idp_origin_vdc == vdc[k]))+
    length(which(dt_data$idp2_origin_vdc == vdc[k]))
  d_count <- length(which(dt_data$vdc == vdc[k]))
  if(o_count>d_count){
    V(gd)$color[k] <- "green"
  } 
}
V(gd)$name <- vdc
gd <- drop_isolated(graph = gd,
                    vertex_colors = V(gd)$color,
                    vertex_names = V(gd)$name)
gd <- drop_loops(graph = gd,
                 vertex_colors = V(gd)$color,
                 vertex_names = V(gd)$name)
ebc <- edge.betweenness.community(graph = gd,
                                  weights = E(gd)$weight,
                                  directed = TRUE)
gd_coords <- coords[which(vdc %in% V(gd)$name),]
plot(ebc,
     gd, 
     layout = gd_coords,
     vertex.size = 4, 
     vertex.label = NA,
     vertex.label.color = "black", 
     vertex.label.font = 1, 
     vertex.label.cex = 0.85, 
     edge.width = 0.1*sqrt(E(gd)$weight),
     edge.arrow.size = 0.5,
     edge.curved = TRUE,
main = "Weighted Nepal Displacement Geo-Network Communities (EB)")
plot(ebc,
     gd,
     layout = layout.fruchterman.reingold(gd, 
                                          niter = 200, 
                                          area = 2000*vcount(gd)),
     vertex.color = V(gd)$color,
     vertex.size = 9,
     vertex.label = V(gd)$name, 
     vertex.label.color = "black",
     vertex.label.font = 1, 
     vertex.label.cex = 0.85, 
     edge.width = 0.2*sqrt(E(gd)$weight),
     edge.arrow.size = 0.5,
     edge.curved = TRUE,
     edge.color = gray.colors(1),
     main = "Weighted Nepal Displacement Network Communities (EB)")

# SOME BASIC EDGE-BETWEENNESS COMMUNITY STATS
plot(sort(ebc$membership), 
     col = adjustcolor(rgb(1,0,1,1)), 
     xlab = "Node Index (VDC) in the Network", 
     ylab = "Edge-Betweenness Community Values", 
     main = "Sorted Edge-Betweenness Community Values for Nepal Displacement Network", 
     pch = 19)
histP1(ebc$membership,
       breaks = 60,
       col = adjustcolor(rgb(1,0,1,1)),
       xlab = "Edge-Betweenness Community Values",
       main = "Edge-Betweenness Community Distribution for Nepal Displacement Network")
#
#
#
# SIMPLEST CLUSTERING/COMMUNITY DETECTION
#
# This clustering is based on the connectivity of the graph.
# In a way, the msot basic one, but, as it can be observed, 
# not necessarily the most informative or functional as
# there are occasional conenctions that ar enonessential
# but would have a significant imapct for this method.
#
#
#
gd <- graph.adjacency(dtm,
                      mode = "directed",
                      weighted = TRUE)
V(gd)$color <- rep("SkyBlue2",length(vdc))
for (k in 1:length(vdc)){
  o_count <- length(which(dt_data$idp_origin_vdc == vdc[k]))+
    length(which(dt_data$idp2_origin_vdc == vdc[k]))
  d_count <- length(which(dt_data$vdc == vdc[k]))
  if(o_count>d_count){
    V(gd)$color[k] <- "green"
  } 
}
V(gd)$name <- vdc
gd <- drop_isolated(graph = gd,
                    vertex_colors = V(gd)$color,
                    vertex_names = V(gd)$name)
gd <- drop_loops(graph = gd,
                 vertex_colors = V(gd)$color,
                 vertex_names = V(gd)$name)
gd_coords <- koords[which(vdc %in% V(gd)$name),]
strongclusters <- clusters(gd)$membership
plot(gd,
     layout = layout.fruchterman.reingold(gd, 
                                          niter = 200, 
                                          area = 2000*vcount(gd)),
     vertex.color = strongclusters,
     vertex.size = 9,
     vertex.label = V(gd)$name, 
     vertex.label.color = "black",
     vertex.label.font = 1, 
     vertex.label.cex = 0.85, 
     edge.width = 0.2*sqrt(E(gd)$weight),
     edge.arrow.size = 0.5,
     edge.curved = TRUE,
     edge.color = gray.colors(1),
     main = "Weighted Nepal Displacement Network Clusters")
plot(gd,
     layout = gd_coords,
     vertex.color = strongclusters,
     vertex.size = 4,
     vertex.label = NA, 
     vertex.label.color = "black",
     vertex.label.font = 1, 
     vertex.label.cex = 0.75, 
     edge.width = 0.2*sqrt(E(gd)$weight),
     edge.arrow.size = 0.5,
     edge.curved = TRUE,
     edge.color = gray.colors(1),
     main = "Weighted Nepal Displacement Geo-Network Flow Clusters")


# SOME BASIC CLUSTERING COMMUNITY STATS
plot(sort(clusters(gd)$membership), 
     col = adjustcolor(rgb(1,0,1,1)), 
     xlab = "Node Index (VDC) in the Network", 
     ylab = "Clustering Community Values", 
     main = "Sorted Clustering Community Values for Nepal Displacement Network", 
     pch = 19)
histP1(clusters(gd)$membership,
       breaks = 60,
       col = adjustcolor(rgb(1,0,1,1)),
       xlab = "Clustering Community Values",
       main = "Clustering Community Distribution for Nepal Displacement Network")
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
gd <- graph.adjacency(dtm,
                      mode = "directed",
                      weighted = TRUE)
V(gd)$color <- rep("SkyBlue2",length(vdc))
for (k in 1:length(vdc)){
  o_count <- length(which(dt_data$idp_origin_vdc == vdc[k]))+
    length(which(dt_data$idp2_origin_vdc == vdc[k]))
  d_count <- length(which(dt_data$vdc == vdc[k]))
  if(o_count>d_count){
    V(gd)$color[k] <- "green"
  } 
}
V(gd)$name <- vdc
gd <- drop_isolated(graph = gd,
                    vertex_colors = V(gd)$color,
                    vertex_names = V(gd)$name)
gd <- drop_loops(graph = gd,
                 vertex_colors = V(gd)$color,
                 vertex_names = V(gd)$name)
gd_coords <- koords[which(vdc %in% V(gd)$name),]
gd <- as.undirected(gd)
mc <- multilevel.community(graph = gd,
                           weights = E(gd)$weights)
plot(mc,
     gd,
     layout = layout.fruchterman.reingold(gd, 
                                          niter = 200, 
                                          area = 2000*vcount(gd)),
     vertex.color = mc,
     vertex.size = 9,
     vertex.label = V(gd)$name, 
     vertex.label.color = "black",
     vertex.label.font = 1, 
     vertex.label.cex = 0.85, 
     edge.width = 0.3*sqrt(E(gd)$weight),
     edge.arrow.size = 0.5,
     edge.curved = TRUE,
     edge.color = gray.colors(1),
     main = "Weighted Nepal Displacement Network Multilevel Communities")
plot(mc,
     gd,
     layout = gd_coords,
     vertex.color = strongclusters,
     vertex.size = 4,
     vertex.label = NA, 
     vertex.label.color = "black",
     vertex.label.font = 1, 
     vertex.label.cex = 0.75, 
     edge.width = 0.2*sqrt(E(gd)$weight),
     edge.arrow.size = 0.5,
     edge.curved = TRUE,
     edge.color = gray.colors(1),
     main = "Weighted Nepal Displacement Geo-Network Multilevel Communities")


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
gd <- graph.adjacency(dtm,
                      mode = "directed",
                      weighted = TRUE)
V(gd)$color <- rep("SkyBlue2",length(vdc))
for (k in 1:length(vdc)){
  o_count <- length(which(dt_data$idp_origin_vdc == vdc[k]))+
    length(which(dt_data$idp2_origin_vdc == vdc[k]))
  d_count <- length(which(dt_data$vdc == vdc[k]))
  if(o_count>d_count){
    V(gd)$color[k] <- "green"
  } 
}
V(gd)$name <- vdc
gd <- drop_isolated(graph = gd,
                    vertex_colors = V(gd)$color,
                    vertex_names = V(gd)$name)
gd <- drop_loops(graph = gd,
                 vertex_colors = V(gd)$color,
                 vertex_names = V(gd)$name)
gd_coords <- koords[which(vdc %in% V(gd)$name),]
gd <- as.undirected(gd)
wc <- walktrap.community(graph = gd,
                           weights = E(gd)$weights)
plot(wc,
     gd,
     layout = layout.fruchterman.reingold(gd, 
                                          niter = 200, 
                                          area = 2000*vcount(gd)),
     vertex.color = mc,
     vertex.size = 9,
     vertex.label = V(gd)$name, 
     vertex.label.color = "black",
     vertex.label.font = 1, 
     vertex.label.cex = 0.85, 
     edge.width = 0.3*sqrt(E(gd)$weight),
     edge.arrow.size = 0.5,
     edge.curved = TRUE,
     edge.color = gray.colors(1),
     main = "Weighted Nepal Displacement Network Walktrap Communities")
plot(wc,
     gd,
     layout = gd_coords,
     vertex.color = strongclusters,
     vertex.size = 4,
     vertex.label = NA, 
     vertex.label.color = "black",
     vertex.label.font = 1, 
     vertex.label.cex = 0.75, 
     edge.width = 0.2*sqrt(E(gd)$weight),
     edge.arrow.size = 0.5,
     edge.curved = TRUE,
     edge.color = gray.colors(1),
     main = "Weighted Nepal Displacement Geo-Network Walktrap Communities")


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
