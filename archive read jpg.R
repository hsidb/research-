#X55stations_temp_precip <- read_excel("Downloads/climatedata/MeteoFrance_55stations/55stations_temp&precip.xlsx")
X55stations_temp_precip <- read_excel("climatedata/MeteoFrance_55stations/55stations_temp&precip.xlsx",
col_types = c("numeric", "text", "text",
"text", "text", "text", "numeric",
"numeric", "numeric", "numeric",
"numeric"))

X55stations_temp_precip <- read_excel("climatedata/MeteoFrance_55stations/55stations_temp&precip.xlsx")
X55stations_temp_precip$LAT1 = gsub('°', ' ', X55stations_temp_precip$LAT)
X55stations_temp_precip$LAT1 = gsub('\'', '.', X55stations_temp_precip$LAT1)
X55stations_temp_precip$LAT1 = gsub('\"N', ' ', X55stations_temp_precip$LAT1)
X55stations_temp_precip$LAT1 = measurements::conv_unit(X55stations_temp_precip$LAT1, from = 'deg_dec_min', to = 'dec_deg')

X55stations_temp_precip$LON1 = gsub('°', ' ', X55stations_temp_precip$LON)
X55stations_temp_precip$LON1 = gsub('\'', '.', X55stations_temp_precip$LON1)
X55stations_temp_precip$LON1 = gsub('\"E', ' ', X55stations_temp_precip$LON1)
X55stations_temp_precip$LON1 = gsub('\"W', '-', X55stations_temp_precip$LON1)
#X55stations_temp_precip$LON1 = gsub('-*', '^-', X55stations_temp_precip$LON1)
merged_weather1$altitude = gsub('m', ' ', merged_weather1$ALT)
merged_weather1$altitude <- as.numeric(merged_weather1$altitude)

X55stations_temp_precip$LON1 = measurements::conv_unit(X55stations_temp_precip$LON1, from = 'deg_dec_min', to = 'dec_deg')
X55stations_temp_precip$LON1 <- as.numeric(X55stations_temp_precip$LON1) 
X55stations_temp_precip$LAT1 <- as.numeric(X55stations_temp_precip$LAT1)

weather_station_list <- unique(X55stations_temp_precip[c("NOM", "LAT1", "LON1")])
weather_station_list$LON1 = measurements::conv_unit(weather_station_list$LON, from = 'deg_dec_min', to = 'dec_deg')

#View(X55stations_temp_precip)

write.table(weather_stations_list, "weather_station_list.xlsx", sep = " ", row.names = FALSE, col.names = TRUE)

#X55stations_temp_precip$date <- as.numeric(X55stations_temp_precip$DATE)
#strptime('2/5/2010', "%m/%d/%Y")
#X55stations_temp_precip$date <- strptime('X55stations_temp_precip$date', "%d%m%Y")
#sapply(X55stations_temp_precip, function(x) length(unique(x)))
#create date column yourself as 01011944 etc with range 

#summary(X55stations_temp_precip)
#sapply(X55stations_temp_precip, function(x) length(unique(x)))
#sapply(X55stations_temp_precip, function(x) sum(is.na(x)))
#summary(X55stations_temp_precip$DATE)
install.packages("pdftools")
library(pdftools)

setwd("~/Documents/Pauline")
getwd()
amiens_glisy_jan44 <- pdf_text("Amiens-Glisy Jan 44.pdf")
text <- amiens_glisy_jan44
text2 <- strsplit(text, "\n")
head(text2[[1]])

setwd("~/Documents/Pauline/text recognized")
getwd()
amiens <- pdf_text("Amiens-Glisy.pdf")
text <- amiens
text2 <- strsplit(text, "\n")
head(text2[[1]])

# img <- readJPEG("MeteoFrance_GermanArchives/amienss/Amiens/amiens2/FRAN_0122_480658_L.jpg", native = TRUE)
# if(exists("rasterImage")){
#   plot(1:2, type='n')
#   rasterImage(img,1,1,2,2)
# }

setwd("~/Documents/Pauline/MeteoFrance_GermanArchives/avords")
getwd()
avords <- pdf_text("avords_43.pdf")
text <- avords
text2 <- strsplit(text, "\n")
(text2[[1]])
 # library(tm)
 # read <- readPDF(control = list(text = "-layout"))
 
#avords <- do.call(rbind.data.frame, text2)
#document <- Corpus(URISource("avords_43.pdf"), readerControl = list(reader = read))
#doc <- content(document[[1]])
#head(doc)
#write.csv(doc, file = "text.csv")

 
#merged_weather3 <-  merged_weather3[with(merged_weather3, order(Date)),]

X55stations_temp_precip <- X55stations_temp_precip[with(X55stations_temp_precip, order(NOM, DATE)),] 

x$long = gsub('°', ' ', x$long)

#convert from decimal minutes to decimal degrees
x$long = measurements::conv_unit(x$long, from = 'deg_dec_min', to = 'dec_deg')
# library(readxl)
# > Amiens_Glisy_Jan_44 <- read_excel("Amiens-Glisy Jan 44.xlsx", 
#                                     +     col_types = c("text"))
# > View(Amiens_Glisy_Jan_44)

summary <- read_excel("summary.xlsx", sheet = "Sheet2")
summary$LAT1 = gsub('°', ' ', summary$latitude)
summary$LAT1 = gsub('N', ' ', summary$LAT1)
#summary$LAT1 = gsub('\'', ' ', summary$LAT1)


summary$LAT1 = measurements::conv_unit(summary$lat, from = 'deg_dec_min', to = 'dec_deg')
summary$LON1 = measurements::conv_unit(summary$long, from = 'deg_dec_min', to = 'dec_deg')

summary$LON1 = gsub('°', ' ', summary$LON)
summary$LON1 = gsub('\'', '.', summary$LON1)
summary$LON1 = gsub('\"E', ' ', summary$LON1)
summary$LON1 = gsub('\"W', '-', summary$LON1)
#X55stations_temp_precip$LON1 = gsub('-*', '^-', X55stations_temp_precip$LON1)


merged_weather3$altitude <- as.numeric(gsub('m', ' ', merged_weather3$ALT))

X55stations_temp_precip$LON1 = measurements::conv_unit(X55stations_temp_precip$LON1, from = 'deg_dec_min', to = 'dec_deg')
X55stations_temp_precip$LON1 <- as.numeric(X55stations_temp_precip$LON1) 
X55stations_temp_precip$LAT1 <- as.numeric(X55stations_temp_precip$LAT1)

weather_station_list <- unique(summary[c("NOM", "LAT1", "LON1")])
weather_station_list$LON1 = measurements::conv_unit(weather_station_list$LON, from = 'deg_dec_min', to = 'dec_deg')

geocoded_attentats_done$Date <- as.Date(with(geocoded_attentats_done,paste(Year,Month,Day, sep = "-")),"%Y-%m-%d")

X55stations_temp_precip$date <- as.character(X55stations_temp_precip$DATE)
X55stations_temp_precip$date <- as.Date(X55stations_temp_precip$date, format = "%d%m%Y")
X55stations_temp_precip$date1 <- as.character.Date(X55stations_temp_precip$DATE)
X55stations_temp_precip$date1 <- as.Date(X55stations_temp_precip$date1, format = "%d%m%Y")

#subset weather data to match attacks
X55_subset <- subset(X55stations_temp_precip, date >= "1943-01-02" & date <= "1944-07-20")
X55_subset11 <- subset(X55stations_temp_precip, a >= "1943-01-02" & date <= "1944-07-20")

sapply(X55_subset,function(x) sum(is.na(x)))
sapply(X55_subset, function(x) length(unique(x)))

# g <- with(X55_subset, aggregate(RR, list(date = date), sum, na.rm = TRUE))
# X55_subset1 = merge(X55_subset, g, by.x= "date", by.y= "date")
# names(X55_subset1)[names(X55_subset1)=="x"] <- "rain_per_day"
# 
# h <- with(X55_subset, aggregate(RR, list(date = date), mean, na.rm = TRUE))
# X55_subset1 = merge(X55_subset1, h, by.x= "date", by.y= "date")
# names(X55_subset1)[names(X55_subset1)=="x"] <- "mean_rain_per_day"
# 
# X55_subset1$week <- format(X55_subset1$date,format = "%Y %W")
# h <- with(X55_subset1, aggregate(RR, list(week = week, NOM = NOM), sum, na.rm = TRUE))
# X55_subset2 = merge(X55_subset1, h)
# names(X55_subset2)[names(X55_subset2)=="x"] <- "sum_rain_week"
# 
# 
# h <- with(merged_weather3, aggregate(TMEAN, list(Date = Date), mean, na.rm = TRUE))
# merged_weather3 = merge(merged_weather3, h)
# names(merged_weather3)[names(merged_weather3)=="x"] <- "mean_temp"



plot(X55_subset$DATE, X55_subset$RR)

data <- merge(activity_per_day, X55_subset1[c(1,16,17)], by.x= "Date", by.y = "date", all.x = TRUE)

plot.ts(data[c(2:4)])

plot(data$ave_rain_per_station
     ,data$activities)


#calculate distance between two points
library(geosphere)
geocoded_attentats_unique_locations <- unique(geocoded_attentats_done[,c(5,7,8)])
mat <- distm(geocoded_attentats_unique_locations[,c('lng','lat')], weather_station_list[,c('LON1','LAT1')], fun=distVincentyEllipsoid)
summary$lat <- as.numeric(summary$lat)
summary$long <- as.numeric(summary$long)


mat <- distm(geocoded_attentats_unique_locations[,c('lng', 'lat')], summary[,c('long','lat')])


# assign the name to the point in list1 based on shortest distance in the matrix
geocoded_attentats_unique_locations$NOM <- weather_station_list$NOM[max.col(-mat)]
geocoded_attentats_unique_locations$distance <- max.col(-mat)
# 
geocoded_attentats_unique_locations$archive_station <- summary$X__1[max.col(-mat)]
geocoded_attentats_unique_locations$dist_archive <- max.col(-mat)
geocoded_attentats_unique_locations <- merge(geocoded_attentats_unique_locations,summary, by.x = "archive_station", by.y = "X__1", all.x = TRUE)
               
#geocoded_attentats_done <- geocoded_attentats_done[with(geocoded_attentats_done, order("NOM", )),]
attentats_weatherstations <- geocoded_attentats_done[c("Location", "NOM")]
attentats_weatherstations <- unique(attentats_weatherstations)


sapply(attentats_weatherstations,function(x) sum(is.na(x)))
sapply(attentats_weatherstations, function(x) length(unique(x)))

merged_attentats <- merge(geocoded_attentats_done,geocoded_attentats_unique_locations, by.x = "Location", by.y = "Location", all.x =  TRUE)
merged_attentats <- 
merged_weather <- merge(geocoded_attentats_done,X55_subset, by.x = c("NOM", "Date"), by.y = c("NOM", "date"))

merged_weather1 <- merge(geocoded_attentats_done,X55_subset2, by.x = c("NOM", "Date"), by.y = c("NOM", "date"))

merged_weather3 <- merge(merged_attentats,X55_subset2, by.x = c("NOM", "Date"), by.y = c("NOM", "date"))

write.csv(merged_weather, file = "merged_weather.csv")
geocoded_attentats_done$weeksdays <- weekdays(geocoded_attentats_done$Date)

r$sum_theft <- with(geocoded_attentats_done,aggregate(theft, list(weeksdays = weeksdays), sum, na.rm = TRUE))
r$sum_sabotage <- with(geocoded_attentats_done,aggregate(sabotage, list(weeksdays = weeksdays), sum, na.rm = TRUE))
r$sum_attentat <- with(geocoded_attentats_done,aggregate(attentat, list(weeksdays = weeksdays), sum, na.rm = TRUE))


names(r)[2] <- "sum_of_activities"
names(data)[3]<-"new_name"
write.csv(r, file = "activities_per_weekday.csv")

 h <- with(merged_weather3, aggregate(RR, list(Date = Date), mean, na.rm = TRUE))
 merged_weather3 = merge(merged_weather3, h)
 names(merged_weather3)[names(merged_weather3)=="x"] <- "mean_rain_per_day"
# 
 merged_weather3$TMEAN <- (merged_weather3$TX + merged_weather3$TN) /2
# 
 h <- with(merged_weather4, aggregate(TMEAN, list(Date = Date), mean, na.rm = TRUE))
 merged_weather4 = merge(merged_weather4, h)
 names(merged_weather4)[names(merged_weather4)=="x"] <- "mean_temp"

merged_weather3$norain_indicator <- ifelse(merged_weather3$RR == 0.0, 1, 0)
h <- with(merged_weather3, aggregate(activities, list(NOM = NOM, Date = Date), sum, na.rm = TRUE))
merged_weather3 = merge(merged_weather3, h)
names(merged_weather3)[names(merged_weather3)=="x"] <- "sum_activities_day_weatherstation"

h <- with(merged_weather1, aggregate(activities, list(NOM = NOM, week = week), sum, na.rm = TRUE))

plot.ts(merged_weather1$sum_activities_perweek_weatherstation)
merged_weather3$rain_indicator <- ifelse(merged_weather3$sum_rain_week <= 10, 1, 0)

hist(merged_weather1$sum_rain_week)

training_data = merged_weather1[which(merged_weather1$Date <= "1944-03-01"),]
test_data = merged_weather1[which(merged_weather1$Date > "1944-01-01"),]

training_data = merged_weather3[which(merged_weather3$Date <= "1944-03-01"),]
test_data = merged_weather3[which(merged_weather3$Date > "1944-01-01"),]


m <- lm(sum_activities_perweek_weatherstation ~ RR + TN + TX + RR:TN + sum_rain_week  + rain_indicator + altitude, merged_weather1)
m <- lm(activities ~ RR + TN + TX + RR:TN + RR:TX + sum_rain_week +  altitude + altitude:TN + altitude:RR + altitude:TX , training_data)

merged_weather3$weeksdays <- as.factor(merged_weather3$weeksdays)
merged_weather3$Department <- as.factor(merged_weather3$Department)

m1 <- lm(activities ~ RR + TN + TX + RR:TN + RR:TX  + factor(weeksdays) + week, merged_weather3) # R = 0.65
m2 <- lm(theft_per_day ~ RR + TN + TX + RR:TN + RR:TX +  factor(weeksdays) + week , merged_weather3) # R = 0.70
#m2a <- lm(activities ~ RR + TN + TX + RR:TN + RR:TX  + factor(weeksdays) + week, merged_weather4)
m3 <- lm(fires_per_day ~ RR + TN + TX + RR:TN + RR:TX +  factor(weeksdays) + week , merged_weather3) # R= 0.57
m4 <- lm(sabotage_per_day ~ RR + TN + TX + RR:TN + RR:TX +  factor(weeksdays) + week  , merged_weather3) # R = 0.45
m5 <- lm(attentat_per_day ~ RR + TN + TX + RR:TN + RR:TX +  factor(weeksdays) + week  , merged_weather3)  # R = 0.36

m6  <- lm(theft_per_day ~ RR + TN + TX + RR:TN + RR:TX + sum_rain_week  + factor(weeksdays) + week   , merged_weather3)

m7 <- lm(activities ~ RR + TN + TX + RR:TN + RR:TX + altitude, merged_weather3) # R = 0.1

m8 <- lm(sum_activities_day_weatherstation ~ RR + TN + TX + RR:TN + RR:TX + altitude  , merged_weather3)

m9 <- lm(theft_per_day ~ RR + TN + TX + RR:TN + RR:TX + altitude+ sum_rain_week + rain_indicator + rain + factor(weeksdays) + week , merged_weather3)


summary(m12)
summary(m3)
summary(m4)
summary(m8)

# sum rain per week per station
plot(m2)
# sum rain per department per day

p <- predict(m, test_data)

#normalize <- function(x) { return((x - min(x))/ (max(x) - min(x)))}

#merged_weather1_norm <- as.data.frame(lapply(merged_weather3, normalize))
plot(merged_weather3$Date, merged_weather3$attentat_per_day)
plot.ts(merged_weather3[,c("theft_per_day", "fires_per_day", "sabotage_per_day", "attentat_per_day")])

library(foreign)
library(grid)
library(arm)
library(texreg)

screenreg(list(m1,m2,m3,m4,m5), single.row = FALSE)

library(caret)
library(corrplot)
ranking <- varImp(m1, scale = FALSE)


geocoded_attentats_unique_locations1 <- merge(geocoded_attentats_unique_locations, points_with_cantons_Hagen[c(4,41:46)], by.x = "Location", by.y = "Location", all.x =  TRUE)
geocoded_attentats_unique_locations1 <- unique(geocoded_attentats_unique_locations1)
#geocoded_attentats_unique_locations1 <- geocoded_attentats_unique_locations1[!is.na(geocoded_attentats_unique_locations1$archive_station),]


sapply(geocoded_attentats_unique_locations1,function(x) sum(is.na(x)))
sapply(geocoded_attentats_unique_locations1,function(x) length(unique(x)))
sapply(merged_weather3,function(x) sum(is.na(x)))
sapply(merged_weather3,function(x) length(unique(x)))


merged_weather4 <- merge(merged_weather3,geocoded_attentats_unique_locations1[c(1,2,7,18,20)],by.x = "Location", by.y = "Location", all.x = TRUE)
#merged_weather6 <- merge(merged_weather4, geocoded_attentats_unique_locations1[c(1,2,7)],by.x = "Location", by.y = "Location")

g <- with(merged_weather4, aggregate(RR, list(Date = Date), sum, na.rm = TRUE))
merged_weather4 = merge(merged_weather4, g, by.x= "Date", by.y= "Date")
names(merged_weather4)[names(merged_weather4)=="x"] <- "sum_rain_per_day"

h <- with(merged_weather4, aggregate(RR, list(Date = Date), mean, na.rm = TRUE))
merged_weather4 = merge(merged_weather4, h, by.x= "Date", by.y= "Date")
names(merged_weather4)[names(merged_weather4)=="x"] <- "mean_rain_per_day"

h <- with(merged_weather4, aggregate(RR, list(week = week), mean, na.rm = TRUE))
merged_weather4 = merge(merged_weather4, h, by.x= "week", by.y= "week")
names(merged_weather4)[names(merged_weather4)=="x"] <- "mean_rain_per_week"

h <- with(merged_weather4, aggregate(RR, list(week = week), sum, na.rm = TRUE))
merged_weather4 = merge(merged_weather4, h, by.x= "week", by.y= "week")
names(merged_weather4)[names(merged_weather4)=="x"] <- "sum_rain_per_week"

# variation in daily activities can't be explained by weather parameters 

m9 <- lm(theft_per_day ~    sum_rain_per_day    + norain_indicator +  mean_temp     , merged_weather4) # R = 0.24
m91 <- lm(activities ~    sum_rain_per_day    + norain_indicator +  mean_temp     , merged_weather4) # R = 0.23
m92 <- lm(fires_per_day ~    mean_rain_per_day    + norain_indicator +  mean_temp     , merged_weather4) # R = 0.02
m93 <- lm(attentat_per_day ~    sum_rain_per_day    + norain_indicator +  mean_temp     , merged_weather4) # R = 0.05
m94 <- lm(sabotage_per_day ~    sum_rain_per_day    + norain_indicator +  mean_temp     , merged_weather4) # R = 0.02

#m10 <- lm(theft_per_day ~ RR + TN + TX + RR:TN + RR:TX +  rain_indicator , merged_weather3)
#m11 <- lm(theft_per_day ~ mean_rain_per_week + sum_rain_per_day + mean_rain_per_day + rain_indicator  + mean_temp, merged_weather4)
m18 <- lm(theft_per_day ~   mean_rain_per_day +  mean_temp  + factor(weeksdays) + week, merged_weather4) # R = 0.71
m16 <- lm(activities ~   sum_rain_per_day + mean_rain_per_day + factor(weeksdays) + week , merged_weather4)  # R = 0.71
m14 <- lm(fires_per_day ~   mean_rain_per_day  + norain_indicator + mean_temp + factor(weeksdays) + week, merged_weather4) # R 0.57
m17 <- lm(attentat_per_day ~    mean_temp  + sum_rain_per_week  + factor(weeksdays) + week, merged_weather4) #R = 0.37
m15 <- lm(sabotage_per_day ~   sum_rain_per_week  +  mean_temp + factor(weeksdays) + week, merged_weather4) # R = 0.44


# variation in theft_per_day is explained by time parameters
m18a <- lm(theft_per_day ~   factor(weeksdays) + week, merged_weather4) # R = 0.69



m13 <- lm(theft_per_day ~  sum_rain_per_day +  rain_indicator + mean_temp , merged_weather4) #R = 0.24
ranking <- varImp(m14, scale = FALSE)
merged_weather5 <- merged_weather4[c(2,71:76)]
merged_weather5 <- unique(merged_weather5)
activity_per_day2 <- merge(activity_per_day, merged_weather5, by.x = "Date", by.y = "Date") 
plot.ts(activity_per_day2[c(-1)])

# local activities can't be explained by variation in weather and time

m19 <- lm(theft ~  RR + TN + TX + RR:TN + RR:TX + factor(weeksdays) + week, merged_weather3) # R = 0.12
#m20 <- lm(theft ~   RR + TN + TX + RR:TN + RR:TX  + factor(weeksdays) + week , merged_weather3)
m21 <- lm(sabotage ~  RR + TN + TX + RR:TN + RR:TX + factor(weeksdays) + week, merged_weather3) # R = 0.08
m22 <- lm(fires ~  RR + TN + TX + RR:TN + RR:TX + factor(weeksdays) + week, merged_weather3) # R = 0.09
m23 <- lm(attentat ~   RR + TN + TX + RR:TN + RR:TX +  factor(weeksdays) + week, merged_weather3) # R = 0.02


summary(m19)

merged_weather6 <- merge(merged_weather4, geocoded_attentats_unique_locations1[c(1,2,7)],by.x = "Location", by.y = "Location")

merged_weather7 <- subset(merged_weather6, archive_station == "Amiens")
merged_weather8 <- subset(merged_weather7, Department ==  "Somme")

anova(m14, test = 'Chisq')
library(car)

vif(m23)

princomp(activity_per_day2[-1])

h <- with(X55_subset2, aggregate(RR, list(date = date), sum, na.rm = TRUE))
X55_subset2 = merge(X55_subset2, h, by.x= "date", by.y= "date")
names(X55_subset2)[names(X55_subset2)=="x"] <- "sum_rain_per_day"

h <- with(X55_subset2, aggregate(RR, list(date = date), mean, na.rm = TRUE))
X55_subset2 = merge(X55_subset2, h, by.x= "date", by.y= "date")
names(X55_subset2)[names(X55_subset2)=="x"] <- "mean_rain_per_day"

X55_subset2$TMEAN <- (X55_subset2$TX - X55_subset2$TN)/2


h <- with(X55_subset2, aggregate(TMEAN, list(date = date), mean, na.rm = TRUE))
X55_subset2 = merge(X55_subset2, h, by.x= "date", by.y= "date")
names(X55_subset2)[names(X55_subset2)=="x"] <- "mean_temp_per_day"

weather_summary <- X55_subset2[c(1,16:20)]
weather_summary <- unique(weather_summary)

merged_weather9 <- merge(merged_weather4, weather_summary, by.x = "Date", by.y = "date", all.x = TRUE)

screenreg(list(m9, m91, m92, m93, m94), single.row = TRUE)

