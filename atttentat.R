geocoded_attentats_done <- read_excel("geocoded_attentats_done.xlsx",
col_types = c("text", "text", "text",
"text", "text", "numeric", "numeric",
"text", "text", "text", "text", "text",
"numeric", "text", "numeric", "numeric",
"numeric", "numeric", "numeric",
"numeric", "text", "text", "text",
"text", "text", "text", "text", "text",
"text", "numeric", "numeric", "numeric",
"numeric", "numeric"))

#create new date column
a <- seq(as.Date("1939/01/01"), by = "day", length.out = 2191)
b <- rep(a,)

geocoded_attentats_done$Date <- as.Date(with(geocoded_attentats_done,paste(Year,Month,Day, sep = "-")),"%Y-%m-%d")

g <- with(geocoded_attentats_done, aggregate(attentat, list(Date = Date), sum, na.rm = TRUE))
geocoded_attentats_done = merge(geocoded_attentats_done, g, by.x= "Date", by.y= "Date")
names(geocoded_attentats_done)[names(geocoded_attentats_done)=="x"] <- "attentat_per_day"

g <- with(geocoded_attentats_done, aggregate(theft, list(Date = Date), sum, na.rm = TRUE))
geocoded_attentats_done = merge(geocoded_attentats_done, g, by.x= "Date", by.y= "Date")
names(geocoded_attentats_done)[names(geocoded_attentats_done)=="x"] <- "theft_per_day"

g <- with(geocoded_attentats_done, aggregate(fires, list(Date = Date), sum, na.rm = TRUE))
geocoded_attentats_done = merge(geocoded_attentats_done, g, by.x= "Date", by.y= "Date")
names(geocoded_attentats_done)[names(geocoded_attentats_done)=="x"] <- "fires_per_day"

g <- with(geocoded_attentats_done, aggregate(sabotage, list(Date = Date), sum, na.rm = TRUE))
geocoded_attentats_done = merge(geocoded_attentats_done, g, by.x= "Date", by.y= "Date")
names(geocoded_attentats_done)[names(geocoded_attentats_done)=="x"] <- "sabotage_per_day"

#sum all daily activities
geocoded_attentats_done$activities <- rowSums(geocoded_attentats_done[,c(39,40,41,42)])

activity_per_day <- geocoded_attentats_done[,c("Date","activities")]
activity_per_day <- unique(activity_per_day)
plot(activity_per_day$Date, activity_per_day$activities)

sapply(geocoded_attentats_unique_locations,function(x) sum(is.na(x)))
sapply(geocoded_attentats_unique_locations, function(x) length(unique(x)))

