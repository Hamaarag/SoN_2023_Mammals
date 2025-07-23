library(dplyr)
library(readr)
library(stringr)
library(lubridate)

# handle T2 first
# read CSV exported from eMammal and generate an RDS data file
SeqData <- read_csv('Data/Image_sequence_metadata_20210523125503.csv')
DepData <- read_csv('Data/deployment_metadata_20210523125503.csv')

# fix time zone and sort by time of photo----
SeqData$`Begin Time` <-  lubridate::ymd_hms(SeqData$`Begin Time`,tz = "Asia/Jerusalem")
SeqData$`End Time` <-  lubridate::ymd_hms(SeqData$`End Time`,tz = "Asia/Jerusalem")
SeqData <- SeqData %>% arrange(`Deployment ID`,`Begin Time`)

# fix species names----
idx <- which(!is.na(str_match(SeqData$`Species Name`,"Hystrix cristata")))
SeqData$`Species Name`[idx] <- "Hystrix indica"
SeqData$`Common Name`[idx] <- "Indian crested porcupine"
idx <- which(!is.na(str_match(SeqData$`Species Name`,"Dipodomys")))
SeqData$`Species Name`[idx] <- "Jaculus species"
SeqData$`Common Name`[idx] <- "Jaculus species"

# remove wrong observations (these were also fixed in eMammal after the Image sequence metadata used here was created)
# single and wrong observation of Oryx leucoryx
idx <- which(!is.na(str_match(SeqData$`Species Name`,"Oryx")))

# remove empty/technician observations with wrong dates----
# single observation in Inland sands with date 2013
idx <- c(idx,which(year(SeqData$`Begin Time`)==2013))
# single observation from Nov 16 2020 in the Sfar
idx <- c(idx,which(year(SeqData$`Begin Time`)==2020))
# wrong dates (2016) from Forest
idx <- c(idx,which(year(SeqData$`Begin Time`)==2016 & grepl("Aderet",SeqData$`Deployment Name`)))
idx <- c(idx,which(year(SeqData$`Begin Time`)==2016 & grepl("Meron.+7",SeqData$`Deployment Name`)))
SeqData <- SeqData[-idx,]

# fix date where necessary
# fix year, day and month correct
idx <- which(year(SeqData$`Begin Time`)==2016 & grepl("Bat Shlomo",SeqData$`Deployment Name`))
idx <- c(idx,which(year(SeqData$`Begin Time`)==2016 & grepl("Meron.+8",SeqData$`Deployment Name`)))
idx <- c(idx,which(year(SeqData$`Begin Time`)==2016 & grepl("Ofer.+9",SeqData$`Deployment Name`)))
SeqData$`Begin Time`[idx] <- str_replace(SeqData$`Begin Time`[idx],"2016","2017")
SeqData$`End Time`[idx] <- str_replace(SeqData$`End Time`[idx],"2016","2017")

# fix times of 8 cameras with wrong dates in Har Hanegev----
# the explicit times expressed here were extracted from the Fulcrum form
idx <- which(year(SeqData$`Begin Time`)==2017 & grepl("Bislach Far 9",SeqData$`Deployment Name`))
dT <- as_datetime("2019-01-14T15:02:00",tz="Asia/Jerusalem")-SeqData$`Begin Time`[idx[length(idx)]] # last image is camera trapper
SeqData$`Begin Time`[idx] <- SeqData$`Begin Time`[idx]+dT
SeqData$`End Time`[idx] <- SeqData$`End Time`[idx]+dT

idx <- which(year(SeqData$`Begin Time`)==2017 & grepl("Bislach Near 4",SeqData$`Deployment Name`))
dT <- as_datetime("2019-01-04T10:37:00",tz="Asia/Jerusalem")-SeqData$`Begin Time`[idx[1]] # only one image
SeqData$`Begin Time`[idx] <- SeqData$`Begin Time`[idx]+dT
SeqData$`End Time`[idx] <- SeqData$`End Time`[idx]+dT

idx <- which(year(SeqData$`Begin Time`)==2017 & grepl("Bislach Near 9",SeqData$`Deployment Name`))
dT <- as_datetime("2019-01-04T10:27:08",tz="Asia/Jerusalem")-SeqData$`Begin Time`[idx[1]] # first & last images are camera trapper
SeqData$`Begin Time`[idx] <- SeqData$`Begin Time`[idx]+dT
SeqData$`End Time`[idx] <- SeqData$`End Time`[idx]+dT

idx <- which(year(SeqData$`Begin Time`)==2017 & grepl("Sde Boker Far 8",SeqData$`Deployment Name`))
dT <- as_datetime("2019-01-09T10:30:00",tz="Asia/Jerusalem")-SeqData$`Begin Time`[idx[length(idx)]] # last image is camera trapper
SeqData$`Begin Time`[idx] <- SeqData$`Begin Time`[idx]+dT
SeqData$`End Time`[idx] <- SeqData$`End Time`[idx]+dT

idx <- which(year(SeqData$`Begin Time`)==2017 & grepl("Sde Boker Near 8",SeqData$`Deployment Name`))
dT <- as_datetime("2018-12-30T14:34:00",tz="Asia/Jerusalem")-SeqData$`Begin Time`[idx[1]] # first & last images are camera trapper
SeqData$`Begin Time`[idx] <- SeqData$`Begin Time`[idx]+dT
SeqData$`End Time`[idx] <- SeqData$`End Time`[idx]+dT

idx <- which(year(SeqData$`Begin Time`)==2017 & grepl("Yeruham Far 3",SeqData$`Deployment Name`))
dT <- as_datetime("2018-12-11T12:05:00",tz="Asia/Jerusalem")-SeqData$`Begin Time`[idx[1]] # first & last images are camera trapper
SeqData$`Begin Time`[idx] <- SeqData$`Begin Time`[idx]+dT
SeqData$`End Time`[idx] <- SeqData$`End Time`[idx]+dT

idx <- which(year(SeqData$`Begin Time`)==2017 & grepl("Yeruham Far 8",SeqData$`Deployment Name`))
dT <- as_datetime("2018-12-11T11:33:00",tz="Asia/Jerusalem")-SeqData$`Begin Time`[idx[1]] # first & last images are camera trapper
SeqData$`Begin Time`[idx] <- SeqData$`Begin Time`[idx]+dT
SeqData$`End Time`[idx] <- SeqData$`End Time`[idx]+dT

idx <- which(year(SeqData$`Begin Time`)==2017 & grepl("Yeruham Near 5",SeqData$`Deployment Name`))
dT <- as_datetime("2018-12-11T13:10:00",tz="Asia/Jerusalem")-SeqData$`Begin Time`[idx[1]] # first & last images are camera trapper
SeqData$`Begin Time`[idx] <- SeqData$`Begin Time`[idx]+dT
SeqData$`End Time`[idx] <- SeqData$`End Time`[idx]+dT

# fix lat lon of Zuriel Far 7 (fixed in eMammal)----
idx <- which(grepl("Zuriel.+7$", SeqData$`Deployment Name`))
SeqData$`Actual Lon`[idx] <- 35.31376

# extract point name----
parsed_point_name <- str_match(SeqData$`Deployment Name`,
                               "^T2_(\\w+\\s??[:alpha:]*)\\s?\\d? (Far|Near|KKL Plantings)\\s?((Semi-)?Shifting)? (\\d_?\\d?)$")

# rearrange variables----
S <- SeqData %>% mutate(Treatment = if_else(is.na(Treatment),Subproject,Treatment)) %>% 
  mutate(unit = if_else(grepl("Dry Shrub",Treatment),"Mediterranean-desert transition zone",
                        if_else(grepl("Grassland",Treatment),"Herbaceous and dwarf-shrub vegetation",
                                if_else(grepl("Maquis",Treatment),"Mediterranean maquis",
                                        if_else(grepl("Conifer",Treatment),"Planted conifer forests",
                                                if_else(grepl("Highlands",Treatment),"Negev highlands",
                                                        if_else(grepl("Arid",Treatment),"Arid desert",
                                                                "Inland sands"))))))) %>% 
  mutate (subunit = if_else(grepl("Judea",Treatment),"Judea",
                            if_else(grepl("Galilee",Treatment),"Galilee",
                                    if_else(grepl("Carmel",Treatment),"Carmel",
                                            if_else(grepl("Grassland",Treatment),"Grassland",
                                                    if_else(str_detect(Treatment,"Long.+Maquis"),"Galilee",NA_character_)))))) %>% 
  mutate (site = parsed_point_name[,2],
          settlements = if_else(unit == "Planted conifer forests","Far",
                                if_else(unit %in% c("Arid desert","Inland sands"),NA_character_,
                                        parsed_point_name[,3])),
          agriculture = if_else(unit %in% c("Arid desert","Inland sands"),parsed_point_name[,3],NA_character_),
          dunes = parsed_point_name[,4],
          point_number = parsed_point_name[,6]) %>%
  mutate (point_name = paste(site,if_else(is.na(settlements),agriculture,settlements),if_else(is.na(dunes),"",dunes),point_number,sep=" ")) %>%
  mutate (subunit=if_else(unit=="Planted conifer forests",
                          if_else(site=="Amatzia","Judea",
                                  if_else(site %in% c("Ramot Naftali","Manara"),"Galilee",
                                          if_else(site %in% c("Bat Shlomo","Elyakim"),"Carmel",subunit))),
                          subunit)) %>% 
  
  rename (SciName = `Species Name`,
          CommonName = `Common Name`,
          Lat = `Actual Lat`,
          Lon = `Actual Lon`,
          eMammalDepID = `Deployment ID`,
          eMammalSeqID = `Sequence ID`,
          BeginTime = `Begin Time`,
          EndTime = `End Time`) %>% 
  mutate (is_mammal = if_else(str_detect(CommonName,
                                         "(Hare|Gazelle|Hyaena|Canid|Wolf|Felid|Fox|Rodent| Cat|Ass|Jackal|Hedgehog|Jaculus|Dog|Gerbil|Ibex|Porcupine|Rat|Boar|Cow|Badger|Horse|Mouse|Sheep|Goat|Mule|Jird|Oryx|Dipodomys|Mongoose|Apodemus|Marten|Deer|Donkey|Camel|Caracal)"),
                              TRUE,FALSE)) %>%
  mutate (is_mammal_species = if_else(is_mammal & !str_detect(CommonName,"(Unknown|species)"),TRUE,FALSE)) %>% 
  mutate (is_mammal_species_med_large = if_else(is_mammal_species & !str_detect(CommonName,"(Hedgehog|Gerbil|Rat|Mouse)"),
                                                TRUE,FALSE)) %>% 
  mutate (is_mammal_species_domestic = if_else(is_mammal_species & str_detect(CommonName,"([Dd]omestic|[Mm]ule|[Cc]amel)"),
                                                TRUE,FALSE)) %>%
  select(unit,subunit,site,settlements,agriculture,dunes,point_number,point_name,SciName,CommonName,is_mammal,is_mammal_species,
         is_mammal_species_med_large,is_mammal_species_domestic,Age,Sex,Count,BeginTime,EndTime,Lat,Lon,eMammalDepID,eMammalSeqID) %>% 
  arrange(unit,subunit,site,settlements,agriculture,dunes,point_number,BeginTime)

# save----
saveRDS(S,'Data/T2_mammal_observations.RDS')
write_csv(S,'Data/T2_mammal_observations.csv')
