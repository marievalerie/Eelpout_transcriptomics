library(ncdf4)
library(lubridate)


#################### BALTIC SEA ##########################

setwd('Trier/UPB - poolseq time line/copernicus_env_data/DarsserOrt')

nc_file <- nc_open('2021-model.nc')
print(nc_file)

#file strucutre
names(nc_file)

#name of variable
names(nc_file$var)

#names of dimensions
names(nc_file$dim)

#the list values coule be also accessed by their position
depth <- nc_file$dim[[1]]$vals ##is the same as:
depth <- nc_file$dim[["depth"]]$vals


latitude <- nc_file$dim[["latitude"]]$vals
longitude <- nc_file$dim[["longitude"]]$vals
time <- nc_file$dim[["time"]]$vals


latitude
longitude
depth


##time is given as timestamp -> units as secs. since defined time point
t_units <- ncatt_get(nc_file, "time", "units")
t_units


names(nc_file$var)
##get potential temperature at ocean bottom
T_array <- ncvar_get(nc_file,nc_file$var[[1]])
T_array



#variable's attributes
ncatt_get(nc_file, "bottomT", "long_name")   #long name
ncatt_get(nc_file, "bottomT", "units")       #measure unit
fillvalue <- ncatt_get(nc_file, "bottomT", "_FillValue")  #(optional)  

##convert the time 
t_ustr <- strsplit(t_units$value, " ")
t_dstr <- strsplit(unlist(t_ustr)[3], "-")
date <- ymd(t_dstr) + dseconds(time)
date


##create csv
coords <- as.matrix(expand.grid(longitude, latitude, date))
temperature <- ncvar_get(nc_file, "bottomT", collapse_degen=FALSE)

nc_df <- data.frame(cbind(coords, temperature))
names(nc_df) <- c("lon", "lat",  "time", "potential sea floor temp ")
nc_df

#csv_fname <- "2021-sea_floor_temperatures_extracted_from_model.csv"
#write.table(nc_df, csv_fname, row.names=FALSE, sep=";")
nc_close(nc_file)


##initialize dataframe all -> only one time
#df_all <- nc_df
df_all<- rbind(df_all, nc_df)


###plot temp
library(stringr)
library(ggplot2)

df <- df_all[, c(3,4)]


df$year <- str_split_fixed(df$time, "-", 3)[,1]
df$year <- as.factor(df$year)
df$year

df$date <- str_split_fixed(df$time, "-", 2)[,2]
df$date <- as.factor(df$date)
df$date

df$temp <- round(as.numeric(df$`potential sea floor temp `), 2)


col <- c("#FFDB6D", "#C4961A", "#F4EDCA", 
         "#D16103", "#C3D7A4", "#52854C", "#4E84C4", "#293352")

ggplot(df, aes(date, temp, color = year, group = year)) +  
  geom_line()+
  scale_colour_manual(values=col)


ggplot(df, aes(year, temp, color = year)) +
  geom_boxplot()

write.table(df, 'temp_data_baltic_sea.csv', row.names=FALSE, sep=";")




#################### NORTH SEA ##########################

setwd('../NorthSea')

nc_file <- nc_open('2021-model.nc')
print(nc_file)



latitude <- nc_file$dim[["latitude"]]$vals
longitude <- nc_file$dim[["longitude"]]$vals
time <- nc_file$dim[["time"]]$vals


latitude
longitude


##time is given as timestamp -> units as secs. since defined time point
t_units <- ncatt_get(nc_file, "time", "units")
t_units

#variable's attributes
ncatt_get(nc_file, "bottomT", "long_name")   #long name
ncatt_get(nc_file, "bottomT", "units")       #measure unit
fillvalue <- ncatt_get(nc_file, "bottomT", "_FillValue")  #(optional)  

##convert the time 
t_ustr <- strsplit(t_units$value, " ")
t_dstr <- strsplit(unlist(t_ustr)[3], "-")
date <- ymd(t_dstr) + dseconds(time)
date


##create csv
coords <- as.matrix(expand.grid(longitude, latitude, date))
temperature <- ncvar_get(nc_file, "bottomT", collapse_degen=FALSE)

nc_df <- data.frame(cbind(coords, temperature))
names(nc_df) <- c("lon", "lat",  "time", "potential sea floor temp ")
nc_df

csv_fname <- "2021-sea_floor_temperatures_extracted_from_model.csv"
write.table(nc_df, csv_fname, row.names=FALSE, sep=";")
nc_close(nc_file)


##initialize dataframe all -> only one time
#df_all <- nc_df
df_all<- rbind(df_all, nc_df)


df <- df_all[, c(3,4)]


df$year <- str_split_fixed(df$time, "-", 3)[,1]
df$year <- as.factor(df$year)
df$year

df$date <- str_split_fixed(df$time, "-", 2)[,2]
df$date <- as.factor(df$date)
df$date

df$temp <- round(as.numeric(df$`potential sea floor temp `), 2)


col <- c("#FFDB6D", "#C4961A", "#F4EDCA", 
         "#D16103", "#C3D7A4", "#52854C", "#4E84C4", "#293352")

ggplot(df, aes(date, temp, color = year, group = year)) +  
  geom_line()+
  scale_colour_manual(values=col)


ggplot(df, aes(year, temp, color = year)) +
  geom_boxplot()

write.table(df, 'temp_data_north_sea.csv', row.names=FALSE, sep=";")
