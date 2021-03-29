library(dplyr)
library(tidyr)
library(reshape2)
library(ggplot2)
library(smooth)

# Download data
dl <- "bash get_data.sh"
system(dl)

COVID19 <- readRDS("COVID-19.rds")

cities <- read.csv("candidate_cities.txt",sep = "\t", header = FALSE)
names(cities) <- c("Name","ID")

# Extract all cities
covcit <- COVID19[COVID19$ID %in% cities$ID,]
# Add city column
covcitid <- covcit %>% left_join(cities, by = "ID")
# Sum duplicate rows based on "City,Date,Type"
covcitid <- aggregate(Cases_New ~Date+Type+Name, data=covcitid, sum, na.rm=TRUE)
covcitid$Name <- gsub("Birmingham", "Birmingham Alabama", covcitid$Name)

# Add UK data
man <- read.csv("manchester_region_deaths.csv",sep = ",", header = TRUE)
bir <- read.csv("birmingham_region_deaths.csv",sep = ",", header = TRUE)
lon <- read.csv("london_region_deaths.csv",sep = ",", header = TRUE)
uk <- rbind(man, bir,lon)
uk <- uk[,c(1,4,7,8)]
names(uk) <- c("date", "Name","Confirmed","Deaths")
# Convert to long
uk <- melt(uk, id = c("date","Name"))
names(uk) <- c("Date","Name","Type","Cases_New")
# Concatenate
covcitid <- rbind(covcitid, uk)
# Sort
covcitid <-covcitid[with(covcitid, order(Name,Type, Date)),]
# Replace NA with zero
covcitid$Cases_New[is.na(covcitid$Cases_New)] <- 0
# Write out raw master table
write.table(covcitid, file = "covid_cities.csv", sep = ",", quote = FALSE,row.names = FALSE)

# Drop outlier cities and days
cleancov <- covcitid
cleancov <- cleancov[cleancov$Type %in% c("Deaths"),]
cleancov <- cleancov[!(cleancov$Name %in% c("Tianjin","Beijing","Xinjiang", "Shanghai","Kyoto")),]
cleancov <- cleancov[!(cleancov$Name %in% c("Tokyo") & cleancov$Cases_New == 299),]
cleancov <- cleancov[!(cleancov$Name %in% c("Mexico City") & cleancov$Cases_New == 1618),]
cleancov <- cleancov[!(cleancov$Name %in% c("Mexico City") & cleancov$Date == "2020-11-21" & cleancov$Cases_New == 1618),]
cleancov <- cleancov[!(cleancov$Name %in% c("Delhi") & cleancov$Cases_New == 905),]
cleancov <- cleancov[!(cleancov$Name %in% c("Stockholm") & cleancov$Cases_New == 2136),]
cleancov <- cleancov[!(cleancov$Name %in% c("Moscow") & cleancov$Cases_New == 2553),]
cleancov <- cleancov[!(cleancov$Name %in% c("Madrid") & cleancov$Cases_New == 8779),]
cleancov <- cleancov[!(cleancov$Name %in% c("Amsterdam") & cleancov$Cases_New == 817),]
cleancov <- cleancov[!(cleancov$Name %in% c("Osaka") & cleancov$Cases_New == 82),]
cleancov <- cleancov[!(cleancov$Name %in% c("Lima") & cleancov$Cases_New == 1473),]

# Correct extreme values manually
# Take the three neighboring values before and after and use mean to replace extreme value
cleancov$Cases_New[which(cleancov$Name == "Mexico City" & cleancov$Date == "2020-10-05")] <- 48
cleancov$Cases_New[which(cleancov$Name == "Mexico City" & cleancov$Date == "2020-11-21")] <- 75
cleancov$Cases_New[which(cleancov$Name == "Mexico City" & cleancov$Date == "2020-11-22")] <- 75
cleancov$Cases_New[which(cleancov$Name == "Mexico City" & cleancov$Date == "2020-11-27")] <- 86
cleancov$Cases_New[which(cleancov$Name == "Mexico City" & cleancov$Date == "2020-11-29")] <- 87
cleancov$Cases_New[which(cleancov$Name == "Lima" & cleancov$Date == "2020-07-23")] <- 122
cleancov$Cases_New[which(cleancov$Name == "Lima" & cleancov$Date == "2020-08-14")] <- 65
cleancov$Cases_New[which(cleancov$Name == "Santiago" & cleancov$Date == "2020-06-04")] <- 56
cleancov$Cases_New[which(cleancov$Name == "Santiago" & cleancov$Date == "2020-06-08")] <- 94
cleancov$Cases_New[which(cleancov$Name == "Santiago" & cleancov$Date == "2020-07-17")] <- 48


# Set negative death counts to 0
cleancov$Cases_New[cleancov$Cases_New<0] <- 0
names(cleancov) <- c("Date","Type", "Name", "Cases_New_raw")
ccov <- cleancov %>% group_by(Name,Type)  %>%
  mutate(Cases_New = cma(Cases_New_raw,7)$fitted)
ccout <- ccov[,c(1,2,3,5)]
write.table(ccout, file = "covid_cities_corrected.csv", sep = ",", quote = FALSE,row.names = FALSE)

# Make plots

#ggplot(cleancov, aes(x=Date, y=Cases_New)) +
#  geom_line() +
#  scale_x_date(date_labels = "%b")+
#  facet_wrap(~ Name,scales = "free") +
#  theme(axis.text=element_text(size=6),
#        axis.title=element_text(size=6,))

#ggplot(ccov) +
#  geom_line(data=ccov,mapping=aes(x=Date, y=Cases_New, size = 1),col="red", size=1, alpha = 0.9) +
#  geom_line(data=ccov,mapping=aes(x=Date, y=Cases_New_raw, size = 1),col="blue", size=1, alpha = 0.25) +            
#  scale_x_date(date_labels = "%b")+
#  facet_wrap(~ Name,scales = "free") +
#  theme(axis.text=element_text(size=6),
#        axis.title=element_text(size=6,))
#
## Show only non-US cities
#clean <- theme(panel.grid.major = element_blank(), strip.background = element_blank(),
#               panel.grid.minor = element_blank(),
#               strip.text = element_text(face = "bold",size=20),
#               legend.title = element_blank(),
#               panel.background = element_blank(), axis.line = element_line(colour = "black"),
#               text = element_text(size=20))
#nonus <- ccov[ccov$Name %in% c("Delhi","Hiroshima","Lima","Melbourne","Madrid","Rome","Sao Paulo","Rio de Janeiro","Santiago","Mexico City","Moscow","Tokyo"),]
#ggplot(nonus) +
#  geom_line(data=nonus,mapping=aes(x=Date, y=Cases_New, size = 1),col="red", size=1, alpha = 0.9) +
#  geom_line(data=nonus,mapping=aes(x=Date, y=Cases_New_raw, size = 1),col="blue", size=1, alpha = 0.25) +            
#  scale_x_date(date_labels = "%b")+
#  facet_wrap(~ Name,scales = "free") +
#  clean+
#  ylab("Deaths")
