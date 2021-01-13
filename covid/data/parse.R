library(dplyr)
library(tidyr)
library(reshape2)
library(ggplot2)
library(smooth)

remove_outliers <- function(x, na.rm = TRUE, ...) {
  qnt <- quantile(x, probs=c(.02, .98), na.rm = na.rm, ...)
  H <- 5 * IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - H)] <- NA
  y[x > (qnt[2] + H)] <- NA
  y
}

# Download data
dl <- "bash get_data.sh" 
system(dl)

COVID19 <- readRDS("COVID-19.rds")
cities <- read.csv("candidate_cities.txt",sep = "\t", header = FALSE)
names(cities) <- c("Name","ID")
# Input table of col1=cityname and col2=ID
# Each City may have multiple IDs in which case data for these needs to be merged
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

# Correct raw death data
cleancov <- covcitid
# Exclude cities with low COVID incidence
cleancov <- cleancov[!(cleancov$Name %in% c("Tianjin","Hong Kong","Hiroshima","Beijing","Xinjiang", "Shanghai","Mexico City", "Kyoto")),]

# Remove negative values
cleancov$Cases_New[cleancov$Cases_New<0] <- 0
options(scipen=999)
# Get mean cases
cc_sum <- covcitid %>%
  group_by(Name,Type) %>%
  summarise(
    n = n(),
    meanCases = mean(Cases_New, na.rm = TRUE)
  )

# Find extreme outlier days
dat_noOutliers<- cleancov %>%
  group_by(Name,Type) %>%
  mutate(outlier = remove_outliers(Cases_New)) %>%
 mutate(outlier = ifelse(Cases_New <10, Cases_New, outlier))
# Exclude outliers that are very low (these are likely correct values)

# Remove outliers
# merge with mean values and replace extreme outliers with mean and apply 5-day centered moving average with 33% increase
corrected <- dat_noOutliers %>% left_join(cc_sum, by = c("Name","Type")) %>% 
  mutate(Cases_New = ifelse(is.na(outlier), meanCases, Cases_New)) %>% 
  group_by(Name,Type)  %>% 
  mutate(Cases_New = cma(Cases_New,7)$fitted) %>%
  mutate(Cases_New = ifelse(Type == "Deaths", Cases_New*1.3, Cases_New))

corrected <- corrected[ ,!(names(corrected) %in% c("outlier","n","meanCases"))]
corrected$Cases_New <- round(corrected$Cases_New,2)
          
write.table(corrected, file = "covid_cities_corrected.csv", sep = ",", quote = FALSE, row.names = FALSE)

# Generate Population Size table
#system("wget https://github.com/CSSEGISandData/COVID-19_Unified-Dataset/blob/master/COVID-19_LUT.csv")
#geodat <- read.csv("COVID-19_LUT.csv",sep = ",", header = TRUE)
#citsiz <- geodat %>% left_join(cities, by = "ID") %>% filter(!is.na(Name)) 
#citsiz <- aggregate(Population ~Name, data=citsiz, sum, na.rm=TRUE)
#ukcit<-data.frame(c("Manchester","London","Birmingham_UK"),c(553230,8170000,1149000))
#names(ukcit)<-c("Name","Population")
#citsiz <- rbind(citsiz, ukcit)
#write.table(citsiz, file = "covid_cities_populations.csv", sep = ",", quote = FALSE,row.names = FALSE)


#### PLOT ####


# Get unsmoothed for plotting comparison
uncor <- dat_noOutliers %>% left_join(cc_sum, by = c("Name","Type")) %>% 
  mutate(Cases_New = ifelse(is.na(outlier), meanCases, Cases_New))
uncor <- uncor[ ,!(names(uncor) %in% c("outlier","n","meanCases"))]
# Rename
names(uncor) <- c("Date","Type","Name","Cases_New_raw")
# Join
compar <- corrected %>% left_join(uncor, by = c("Date","Name","Type"))
compar <- melt(compar, id=c("Date","Name","Type"))

dc <- subset(compar , Type == "Deaths")
names(dc) <- c("Date","Name","Type","Deaths","Count")
ggplot(dc, aes(x=Date, y=Count,color = Deaths)) + 
  geom_line(size = 1) + 
  facet_wrap(~ Name,scales = "free") +
  theme(axis.text=element_text(size=6),
        axis.title=element_text(size=6,))

library(scales)
ggplot(dc) + 
  geom_line(data=dc[dc$Deaths=="Cases_New",],mapping=aes(x=Date, y=Count, size = 1),col="red", size=1, alpha = 0.9) +
  geom_line(data=dc[dc$Deaths=="Cases_New_raw",],mapping=aes(x=Date, y=Count, size = 1),col="blue", size=1, alpha = 0.25) +            
  geom_vline(xintercept=as.numeric(as.Date("2020-06-01")), linetype=4) +
  scale_x_date(date_labels = "%b")+
  facet_wrap(~ Name,scales = "free") +
  theme(axis.text=element_text(size=6),
        axis.title=element_text(size=6,)) +
  ggtitle("Covid-19 Deaths per City")
ggsave("covid_city_deaths_cor_raw.pdf")

# Compare AZ (no DST) with CA
dst <- dc[dc$Name %in% c("Tucson","Phoenix","San Diego","Los Angeles"),]
ggplot(dst) + 
  geom_line(data=dst[dst$Deaths=="Cases_New",],mapping=aes(x=Date, y=Count, size = 1),col="red", size=1, alpha = 0.9) +
  geom_line(data=dst[dst$Deaths=="Cases_New_raw",],mapping=aes(x=Date, y=Count, size = 1),col="blue", size=1, alpha = 0.25) +            
  geom_vline(xintercept=as.numeric(as.Date("2020-03-14")), linetype=4) +
  geom_vline(xintercept=as.numeric(as.Date("2020-11-7")), linetype=4) +
  scale_x_date(date_labels = "%b")+
  facet_wrap(. ~ Name,scales = "free",ncol = 1) +
  theme(axis.text=element_text(size=6),
        axis.title=element_text(size=6,)) +
  ggtitle("Covid-19 Deaths per City")
                      
corrected $Date <- as.Date(corrected$Date)
d <- subset(corrected , Type == "Deaths")
ggplot(d, aes(x=Date, y=Cases_New,color = Type)) + 
  geom_line(size = 1) + 
  facet_wrap(~ Name,scales = "free") +
  theme(axis.text=element_text(size=6),
        axis.title=element_text(size=6,))

h <- subset(corrected , Type == "Hospitalized")
ggplot(h, aes(x=Date, y=Cases_New,color = Type)) + 
  geom_line(size = 1) + 
  facet_wrap(~ Name,scales = "free") +
  theme(axis.text=element_text(size=6),
        axis.title=element_text(size=6,))

c <- subset(corrected , Type == "Confirmed")
ggplot(c, aes(x=Date, y=Cases_New,color = Type)) + 
  geom_line(size = 1) + 
  facet_wrap(~ Name,scales = "free") +
  theme(axis.text=element_text(size=6),
        axis.title=element_text(size=6,))

cleancov$Date <- as.Date(cleancov$Date)
d <- subset(cleancov, Type == "Deaths")
c <- subset(cleancov, Type == "Confirmed")
ggplot(d, aes(x=Date, y=Cases_New,color = Type)) + 
  geom_line(size = 1) + 
  facet_wrap(~ Name,scales = "free") +
  theme(axis.text=element_text(size=6),
        axis.title=element_text(size=6,))
ggplot(c, aes(x=Date, y=Cases_New,color = Type)) + 
  geom_line(size = 1, color = "steelblue") + 
  facet_wrap(~ Name,scales = "free") +
  theme(axis.text=element_text(size=6),
        axis.title=element_text(size=6,))







