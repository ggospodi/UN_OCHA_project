# Building a TAble that combines the Aid and the Severity Index tables
# author: Georgi D. Gospodinov
# date: "September 18, 2015"
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


#
#
#
#
# SEVERITY ANALYSIS ON THE VDC PROJECTION OF THE AGENCY-VDC AID NETWORK
# 
#
#
#
#


# INCORPORATE SEVERITY SCORES

sev <- read.csv(paste0(DIR,"severity.csv"))
sev$vdc <- as.character(sev$vdc)
sev$vdc <- mapvalues(sev$vdc,
                     from = c("Akhikarichaur","AmargadhiN.P.","Andheri (Narayansthan)","Arbabijaya",
                              "Arjun Chaupari","Arthar Dadakharka","Awal Parajul","Baad Bhanjyang",
                              "Bada Khola","Badabhairab","Badagama","Badaka Diyale","BageswariTitrona",
                              "BaglungN.P.","Bahadurgunj","Baikunthe","Baklauri","BanepaN. P.",
                              "Banskhaur","Bekhsimle","BeriyaBirta","BeriyaBirta(Wa. Pu.)",
                              "Bhadauretamagi","BhadgaunSinawari","BhadrapurN.P.","BharatpurN. P.",
                              "Bhinggithe","Bicharichautara","BidurN.P.","Bijayshwari(Chaurjahari)",
                              "Bhusedanda",
                              "Binun",
                              "Chande",
                              "Dahafatgaun",
                              "Dhaibung",
                              "Khadang",
                              "Bindyabasini",
                              "BiratnagarN.P.",
                              "BirendranagarN.P.",
                              "BirgunjN. P.",
                              "Bodebarsain","BodhaPokharathok",
                              "Bungadovan","ButawalN.P.","Byas N. P.",
                              "Chaphamandaun","Chhatedhunga","Chimadi","Chok chisapani",
                              "DamakN.P.",
                              "Danda Bazar","DandaParajul","DasharathchandaN.P.","Dewaldibyapur",
                              "Dhamia","DhangadhiN.P.","DhankutaN.P.","DhapukSimalBhanjyan",
                              "Dhar(Institutional)",
                              "Dharampaniya","DharanN. P.","Dhikurpokhari","Dholagohe","DhulikhelN. P.",
                              "DhullaJaidi",
                              "DipayalSilgadhiN.P.","Drabesh","Duhilamati","Dulegauda","Durling",
                              "Falametar","GairiBisouna Deupur","Gamhariya Parwaha","GaurN.P.",
                              "Ghanghalgau",
                              "Ghasikuwa","Gobargadha","GulariyaN.P.","Guthi Prasauni","Hadiryapaltuwa",
                              "Harpatgunj","Hemja","IlamN.P.","InaruwaN .P.",
                              "InarwaFulbariya",
                              "ItahariN. P.",
                              "Jaiamangalapur","Jaira","JaleshworN.P.","JanakpurN.P.","Jhouwaguthi",
                              "KalaiyaN.P.",
                              "Kalati Bhumidanda","KamalamiN.P.","Kanhu","Kanhushivapur","Kaphalkot",
                              "Kapilbastu N. P.","KauwabanKataiya","KavreNitya Chandeswor",
                              "KhandbariN.P.",
                              "Khaptad National Par","Kheen","Khilungdeurali","Khokhsarparbaha",
                              "KholegaunKhanigau","Khukhani","Kinhu","Ko. Madhepura","Koshi Tappu Wildlife",
                              "Kotbhairab","Kushma","Kusmi","LahanN.P.","Lalapatti","Lalutantikot","Laphagaun",
                              "Lekhnath N. P.","Lumbini Development","Machapuchre","MagaPauwa","Maharajgunj",
                              "MahendranagarN.P.","Mainamaini","MalangawaN.P.","Malhaniya Gamharia",
                              "MalhaniyaKhori","Malikabota","Malkot(Institutional)","Malm","Mashtabandali",
                              "Mashtamandaun",
                              "Matsyapokhari","MechinagarN.P.","Mumrakot","Nachnechaur",
                              "NarayanN.P.","NasikasthanSangha","NawamiDanda","NepalgunjN.P.","Nyawal",
                              "Paiyunnyap","Palungmainadi","PanautiN. P.","Panchawati","Paschimkusaha",
                              "Patlepani","Patto","Phoi Mahadev","PiparpatiJabdi","Pokhara N. P.","Prajvunpata",
                              "PrithbinarayanN.P.","Purbakusaha","Purtimkanda","PutalibazarN.P.",
                              "PyuthanKhalanga","RajbirajN. P.","RakhuPiple","RamgramN. P.","RamgunjBelgachhi",
                              "Ramjadeurali","Ranipokhari (Resing)","RatnanagarN.P.","Rayalebhir",
                              "Royal Bardiya Nation","Royal Chitawan Natio","Royal Shuklaphanta",
                              "Rumale(Khaumale)","Sabhapokhari","Sahebgunj","Salmechakala","Saraswar",
                              "SarasyunKhark","ShankarSaraiya","ShankhupatiChour","Shibaratha","Shyamgha",
                              "SiddharthNagarN.P.","Sinwang","SirahaN.P.","SisaKhani","SiswaBelhi","Siuna",
                              "Sreesiya(Nau.Ta.Ja)",
                              "Sukitaya",
                              "Sunaulabazar","Sundhara (Ghiring)","Sunhoo",
                              "Tanahusur","TansenN.P.","Thulolumpek","TikapurN.P.","TribhuwanNagarN. P.",
                              "TriyugaN. P.","TulsipurN. P.","WalingN.P.",
                              "Agara","Aabukhaireni","Baruneshwor","Betini","BhaktapurN.P.","BhimesworN.P.",
                              "ChandeniMandan","Chhatara","GunsiBhadaure","HetaudaN.P.","JaisithokMandan",
                              "JhangajholiRalmata","Jhyaku","JyamdiMandan","KakurThakur","KathmanduN.P.",
                              "LalitpurN.P.","Sangu","KirtipurN.P.","TokhaChandeswori","Thulogoun",
                              "Talkududechour","Sankhu","Puranagau","PokhariNarayansthan",
                              "Pukhulachhi","NaikapPuranoBhanjya","Mankha","Mahankal","MadhyapurThimiN.P.",
                              "Lamidada","Laharepouwa","Daxinkali","Bajrayogini","Budanilkantha","Fulpingkatti",
                              "Orang"),
                     to = c("Adhikarichaur","Amargadhi Municipality","Narayansthan","Arba Bijaypur",
                            "Arjunchaupari","Arthar Dandakharka","Awalparajul","Badbhanjyang",
                            "Badakhola","Badabhairav","Badgama","Baadkadiyale","Bageshwari Titarauna",
                            "Baglung Municipality","Bahadurganj","Baikhunthe","Bakloura","Banepa Municipality",
                            "Banskhor","Bekhsimle Ghartigaon","Bairiyabirta (Nau.Ta.Ja.)","Bairiyanbirta (Wa.Pu.)",
                            "Bhadauretamago","Bhadgaun Sinuwari","Bhadrapur Municipality","Bharatpur Municipality",
                            "Bhingithe","Bichari Chautara","Bidur Municipality",
                            "Chourjahari",
                            "Bhasedawa",
                            "Binauna",
                            "Change",
                            "Daha Phalgaun",
                            "Darbung",
                            "Kharang",
                            "Bindhyabasini",
                            "Biratnagar Sub Metropolitan",
                            "Birendranagar Municipality",
                            "Birgunj Sub Metropolitan",
                            "Barashine (Bode)","Bodha Pokhara Thok",
                            "Bungadobhan","Butawal Municipality","Byas Municipality",
                            "Chaphamadaun","Chatedhunga","Chimdi","Chokchisapani",
                            "Damak Municipality",
                            "Dadhabazaar","Dandaparajul","Dasharath Chanda Municipality","Dewaldebhyapur",
                            "Dhamja","Dhangadhi Municipality","Dhankuta Municipality","Dhapuk Simalbhanjyang",
                            "Dhari",
                            "Dharmapaniya","Dharan Municipality","Dhikur Pokhari","Dholagoha","Dhulikhel Municipality",
                            "Jaidi","Dipayal Silgadhi Municipality","Darbesa","Dudilabhati","Dulegaunda","Durlung",
                            "Phalametar","Gairi Bisauna Deupur","Gamariya Parawaha","Gaur Municipality",
                            "Khagalgaun",
                            "Ghansikuwa","Gobargada","Gulariya Municipality","Guthiparsauni","Hardiyapaltuwa",
                            "Harpatganj","Hyangja","Ilam Municipality","Inarwa Municipality",
                            "Inarwaphulwariya",
                            "Itahari Municipality",
                            "Jayamangalapur","Jair","Jaleshwor Municipality","Janakpur Municipality","Jhauwaguthi",
                            "Kalaiya Municipality",
                            "Kapali Bhumaedanda","Kamalami Municipality","Kanku","Kahunshivapur","Kaphal Kot",
                            "Kapilbastu Municipality","Kauwaban Kataiya","Kabhrenitya Chandeshwari",
                            "Khandbari Municipality",
                            "Khaptad National Park","Khin","Khilu Deurali","Khoksarparwaha",
                            "Khanigaun","Kakani","Kihun","Komadhepura","Koshi Tappu Wildlife Reserve",
                            "Kotbhairav","Kusma","Kusma","Lahan Municipality","Lalapatthi","Lalu","Lafagaun",
                            "Lekhnath Municipality","Lumbini Development Area","Machhapuchchhre","Magapauwa","Maharajganj",
                            "Mahendranagar Municipality","Mainamaine","Malangawa Municipality","Malahaniya Gamariya",
                            "Malhaniya Khori","Malikabota (Hatsinja)","Malkot","Malma","Masyawandali","Mashtamadaun",
                            "Matshyapokhari","Mechinagar Municipality","Mumra","Kristi Nachnechaur",
                            "Narayan Municipality","Nasikasthan Sanga","Nawamidanda","Nepalgunj Municipality",
                            "Nwali",
                            "Paiyunthanthap","Palung Mainadi","Panauti Municipality","Panchabatti","Kusahapaschim",
                            "Patelepani","Pato","Phoimahadev","Khutawa Jabdi","Pokhara Sub Metropolitan","Paiyunpata",
                            "Prithbinarayan Municipality","Purba Kusahha","Purtimdanda","Putalibazar Municipality",
                            "Khalanga","Rajbiraj Municipality","Piple","Ramgram Municipality","Ramgunjbelgachhiya",
                            "Ramja Deurali","Resing Ranipokhari","Ratnanagar Municipality","Rayal",
                            "Royal Bardiya National Park","Royal chitwan National Park","Royal Shukla Phanta National Park",
                            "Khamale","Shavapokhari","Shahebgunj",
                            "Salme","Sarashwor","Sarsyunkharka","Shankarsaraiya","Sangkhupatichaur","Shiba","Syamgha",
                            "Siddhartha Nagar Municipality","Sunwal","Siraha Municipality","Sisakhani","Shishwa","Syuna",
                            "Sirsiya",
                            "Sukatinya",
                            "Sunaulabajar","Sundhara (Thiring )","Sangkhu","Tanahunsur","Tansen Municipality",
                            "Thulo Lumpek","Tikapur Municipality","Tribhuwan Nagar Municipality","Triyuga Municipality",
                            "Tulsipur Municipality","Waling Municipality",
                            "Agra","Anbukhaireni","Barudeshwor","Beteni","Bhaktapur Municipality","Bhimeswor Municipality",
                            "Chandeni Mandan","Chautara","Gunsi","Hetauda Municipality","Jaisithok Mandan",
                            "Jhangajholi Ratmata","Jhyanku","Jyamdi Mandan","Kakur Thakur","Kathmandu Metropolitan",
                            "Lalitpur Sub Metropolitan","Sangkhu","Kirtipur Municipality","Tokhachandeshwari",
                            "Thulo Gaun","Talkudunde Chaur","Sangkhu Suntol","Puranagaun","Pokhari Narayansthan",
                            "Pukulachhi","Naikap Naya","Mangkha","Mahangkal","Madhyapur Thimi Municipality",
                            "Lamidanda","Laharepauwa","Dakshinkali","Sangkhu Bajrayogini","Budhanilkantha",
                            "Phulpingkatti","Worang"))

sev$vdc[218] <- "Betini"
sev$vdc[624] <-"Lamidada"

# MERGE WITH AGENCY-VDC AID TABLE
for (k in 1:dim(sev)[1]){
  if (sev$vdc[k] %in% hlcit$vname){
    sev$hlcit[k] <- hlcit[which(hlcit$vname==sev$vdc[k])[1],]$hlcit
  } else {
    sev$hlcit[k] <- NA
  }
}

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

# EXPORT THE DATA
write.csv(aid_sev,file=paste0(DIR,"aid_and_severity.csv"))

# EXPORT JUST SEVERITY DATA ON THE AID RECEIVING VDCs
sev_list <- c("hazard","exposure","housing","poverty","vulnerability","severity")
sev_aid <- sev[sev$hlcit %in% unique(aid_data$hlcit),sev_list]
