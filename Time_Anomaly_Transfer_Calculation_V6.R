###########################################################
#Time for Anomaly Transfer
#Created By: Debasish Mishra
#Date: June 25th, 2023
###########################################################

#Objective: Looking at the relationship between time for anomaly transfer and Wasserstein distance.

#Let's call in the required libraries first:
library(pacman)
pacman::p_load(tidyverse, terra, sf, sp, tmap, spData, spDataLarge, lubridate, quantreg, raster)

#Load in the required datasets:

#1. Wasserstein Distance 3-month running Spatraster
WD_monthly_rast<-rast("G:/Shared drives/VZRG_TAMU/Debasish/Derived Rasters/Wasserstein Distance between Soil Moisture and Evapotranspiration Anomalies/Global_WD_RunningMean_Jan2010_Dec2019.nc")
names(WD_monthly_rast)<-seq(as.Date("2010/1/1"), by = "month", length.out = 120)

#2. Chemical Potential Gradient (J/kg)
#water_pot_rast<-rast("G:/Shared drives/VZRG_TAMU/Debasish/Derived Rasters/Chemical Water Potential/Water_Chemical_Potential_Gradient_JperKg.nc")
water_pot_rast<-rast("G:/Shared drives/VZRG_TAMU/Debasish/Derived Rasters/Chemical Water Potential/Chemical_water_potentail_gradient_Monthly_Mean_Jan2010_Dec2019_JperKg.nc")
names(water_pot_rast)<-seq(as.Date("2010/1/1"), by = "month", length.out = 120)

#3. Water-Energy Classification Spatraster
WEC_rast<- rast("G:/Shared drives/VZRG_TAMU/Debasish/Other Modelled Datasets/WEC Climate Classification/WEC15_raster_2021.nc")
#Reclassify into 5 broad divisions
WEC_class<-c(1, 1,  
             2, 1,  
             3, 1,   
             4, 2,  
             5, 2,  
             6, 2,   
             7, 3,   
             8, 3,   
             9, 3,   
             10, 4,  
             11, 4,  
             12, 4,  
             13, 5,  
             14, 5,  
             15, 5)  
rclWEC <- matrix(WEC_class, ncol=2, byrow=TRUE)
WEC_rast2 <- classify(WEC_rast, rclWEC, include.lowest=TRUE)

#Resample Water Potential and WEC raster according to WD raster
water_pot_rast<-resample(water_pot_rast, WD_monthly_rast)
WEC_rast2<-resample(WEC_rast2, WD_monthly_rast)


########################################################################
#Function to extract and compute time of anomaly transfer

#Create a spatial point for college station
college_station<-st_sfc(st_point(x = c(-96.33, 30.62), dim = "XY"), crs = "EPSG:4326")

cs_mu_ts<-as.data.frame(t(terra::extract(water_pot_rast, vect(college_station))[,-1]))
cs_WD_ts<-as.data.frame(t(terra::extract(WD_monthly_rast, vect(college_station))[,-1]))

cs_mu_ts$Date<-as.Date(row.names(cs_mu_ts))
cs_WD_ts$Date<-as.Date(row.names(cs_WD_ts))

names(cs_mu_ts)<-c("mu_grad", "Date")
names(cs_WD_ts)<-c("WD", "Date")

cs_ts<-inner_join(cs_WD_ts, cs_mu_ts)
# plot(as.Date(cs_ts$Date), abs(cs_ts$mu_grad), type = "l", lty = 1)
# plot(as.Date(cs_ts$Date), cs_ts$WD^2, type = "l", lty = 1)
# plot(abs(cs_ts$mu_grad), (cs_ts$WD^2), type = "l", lty = 1)
# plot(lm(cs_ts$WD^2 ~ abs(cs_ts$mu_grad)))



#Texas_polygon<-us_states %>% dplyr::filter(NAME == "Texas")

#cs_ts$Date[month(cs_ts$Date) %in% 4:5]

# plot(abs(cs_ts_sub$mu_grad), (cs_ts_sub$WD^2))


Inst_Slope<-c() #Instantaneous Slope vector
months_vector<-list(c("03","04","05"), c("06", "07", "08"), c("09", "10", "11"), c("12", "01", "02"))
length(months_vector)


slope_fun<-function(x){
  
  get_dates_with_months <<- function(date_vector, months_vector) {
    new_dates <- list()
    for(i in 1:length(date_vector)) {
      month <- strsplit(as.character(date_vector[i]), '-')[[1]][2]
      if(month %in% months_vector) { new_dates[i] <- list(date_vector[[i]]) }
    }
    res <- Filter(Negate(is.null), new_dates)
    return(res)
  }
  
  
  print("Stage 1")
  #ts_val<-as.numeric(t(terra::extract(texas_rast, vect(college_station))[,-1]))
  ts_val<- as.numeric(x)
  mu_ts <- ts_val[1:(length(ts_val)/2)]
  WD_ts <- ts_val[(1+(length(ts_val)/2)):(length(ts_val))]
  
  ts_df<- data.frame(mu_ts, WD_ts)
  ts_df$Date <- as.Date(seq(as.Date("2010/1/1"), by = "month", length.out = 120))
  ts_df<- ts_df[complete.cases(ts_df), ]
  names(ts_df) <- c("mu", "WD", "Date")
  
  tryCatch( if(nrow(ts_df) > 30) {
    print("Stage 2")
    Inst_Slope<-c()
    #Inst_Slope2<-c()
    months_vector<-list(c("03","04","05"), c("06", "07", "08"), c("09", "10", "11"), c("12", "01", "02"))
    
    for (i in 1:length(months_vector)) {
      print("Stage 3")
      #Loop over Months vector
      dateList <- get_dates_with_months(ts_df$Date, months_vector[[i]])
      #dateList <-get_dates_with_months(ts_df$Date, months_vector[[2]])
      dateVector <- Reduce(c, dateList)
  
      ts_df_sub <- ts_df[ts_df$Date %in% dateVector,]
  
        if(nrow(ts_df_sub) > 10) {
          print("Stage 4")
          
          #Quantile slope
          rqfit <- rq((WD) ~ abs(mu), data = ts_df_sub)  # Quantile Regression
          Inst_Slope[i] <- as.numeric(as.data.frame(rqfit$coefficients)[2,])  #Slope
          print(Inst_Slope[i])
          
        } else {
          
          print("Stage 5")
          
          Inst_Slope[i] <- NA
          #Inst_Slope2[i] <- NA
        }
    } 
    print("Stage 6")
    print(round(Inst_Slope*(10e10)))
    return(as.numeric((Inst_Slope)))
  } else {
    
    print("Stage 7")
    return(rep(NA, 4))
    
  }, error = function(e){return(rep(NA, 4))}) 
}

#Let's try the function on texas 
South_US<-us_states %>% dplyr::filter(NAME %in% c("Texas", "Oklahoma", "Louisiana", "Arkansas")) %>% st_union()
South_US<-st_transform(South_US, crs = 4326)

mu_WD_rast<-c(water_pot_rast, WD_monthly_rast)

South_US_rast<-crop(mu_WD_rast, vect(South_US))
South_US_rast<-mask(South_US_rast, vect(South_US))
#plot(South_US_rast[[20]])

slope_Out1 = app(mu_WD_rast, fun = slope_fun, cores = 1)
plot(slope_Out1)
writeCDF(slope_Out1, "Seasonal_Slope_with_WD_and_Chem_Water_Pot_Grad.nc")
############################


slope_fun2<-function(x){
  
  get_dates_with_months <<- function(date_vector, months_vector) {
    new_dates <- list()
    for(i in 1:length(date_vector)) {
      month <- strsplit(as.character(date_vector[i]), '-')[[1]][2]
      if(month %in% months_vector) { new_dates[i] <- list(date_vector[[i]]) }
    }
    res <- Filter(Negate(is.null), new_dates)
    return(res)
  }
  
  
  print("Stage 1")
  #ts_val<-as.numeric(t(terra::extract(South_US_rast, vect(college_station))[,-1]))
  ts_val<- as.numeric(x)
  mu_ts <- ts_val[1:(length(ts_val)/2)]
  WD_ts <- ts_val[(1+(length(ts_val)/2)):(length(ts_val))]
  
  ts_df<- data.frame(mu_ts, WD_ts)
  ts_df$Date <- as.Date(seq(as.Date("2010/1/1"), by = "month", length.out = 120))
  ts_df<- ts_df[complete.cases(ts_df), ]
  names(ts_df) <- c("mu", "WD", "Date")
  
  tryCatch( if(nrow(ts_df) > 30) {
    print("Stage 2")
    Inst_Slope<-c()
    #Inst_Slope2<-c()
    months_vector<-list(c("03","04","05"), c("06", "07", "08"), c("09", "10", "11"), c("12", "01", "02"))
    
    for (i in 1:length(months_vector)) {
      print("Stage 3")
      #Loop over Months vector
      dateList <- get_dates_with_months(ts_df$Date, months_vector[[i]])
      #dateList <-get_dates_with_months(ts_df$Date, months_vector[[2]])
      dateVector <- Reduce(c, dateList)
      
      ts_df_sub <- ts_df[ts_df$Date %in% dateVector,]
      
      if(nrow(ts_df_sub) > 10) {
        print("Stage 4")
        
        #Piecewise slope
        Piecewise_Slope<-c()
        chunk<-3
        r  <- rep(1:ceiling(nrow(ts_df_sub)/chunk), each=chunk)[1:nrow(ts_df_sub)]
        d <- split(ts_df_sub, r)
          for (s in 1:length(d)) {
            
            print("Stage 4.1")
            #s=2
            df<-as.data.frame(d[s])
            names(df)<-c("mu", "WD", "Date")
            Lmfit <- lm((WD) ~ abs(mu), data = df)
            Seasonal_slope <- as.numeric(as.data.frame(Lmfit$coefficients)[2,])
            Piecewise_Slope<- c(Piecewise_Slope, Seasonal_slope)
            print(Piecewise_Slope)
          }
        
        Inst_Slope2dens <- density(Piecewise_Slope)
        Inst_Slope[i] <- as.numeric(Inst_Slope2dens$x[which.max(Inst_Slope2dens$y)])
        
        print(Inst_Slope[i])
        
      } else {
        
        print("Stage 5")
        
        Inst_Slope[i] <- NA
        #Inst_Slope2[i] <- NA
      }
    } 
    print("Stage 6")
    print(round(Inst_Slope*(10e10)))
    return(as.numeric((Inst_Slope)))
  } else {
    
    print("Stage 7")
    return(rep(NA, 4))
    
  }, error = function(e){return(rep(NA, 4))}) 
}


slope_Out2 = app(mu_WD_rast, fun = slope_fun2, cores = 1)
hist(slope_Out2)
plot(slope_Out2)
writeCDF(slope_Out2, "Seasonal_Slope2_with_WD_and_Chem_Water_Pot_Grad.nc")
slope_Out1

############################
#Seasonal Analysis

mu_seasonal<-rast("G:/Shared drives/VZRG_TAMU/Debasish/Derived Rasters/Chemical Water Potential/Chemical_water_potentail_gradient_Seasonal_Mean_2010_2019_JperKg.nc")
names(mu_seasonal)<-c("mu_MAM", "mu_JJA", "mu_SON", "mu_DJF")

WD_rast<-rast("G:/Shared drives/VZRG_TAMU/Debasish/Derived Rasters/Wasserstein Distance between Soil Moisture and Evapotranspiration Anomalies/Global_WD_seasonal_Mean_2010_2019.nc")
names(WD_rast)<-c("WD_MAM", "WD_JJA", "WD_SON", "WD_DJF")

#slope_Out2<-rast("G:/Shared drives/VZRG_TAMU/Debasish/Derived Rasters/Chemical Water Potential/Seasonal_Slope2_for_WD_and_Chem_Water_Pot_Grad.nc")
slope_Out2<-rast("G:/Shared drives/VZRG_TAMU/Debasish/Derived Rasters/Chemical Water Potential/Seasonal_Slope2_with_WD_and_Chem_Water_Pot_Grad.nc")

names(slope_Out1)<-c("qrslope_MAM", "qrslope_JJA", "qrslope_SON", "qrslope_DJF")
names(slope_Out2)<-c("plslope_MAM", "plslope_JJA", "plslope_SON", "plslope_DJF")

mu_seasonal<-resample(mu_seasonal, WD_rast)
slope_Out1<-resample(slope_Out1, WD_rast)
slope_Out2<-resample(slope_Out2, WD_rast)

slope_Out1[slope_Out1 == 0]<-NA
slope_Out2[slope_Out2 == 0]<-NA

#plot(slope_Out2)
#Focal aggregation (7*7 window)
slope_Out2_focal = focal(slope_Out2, w = matrix(1, nrow = 7, ncol = 7), fun = median, na.rm=TRUE)

slope_Out2_focal_mean<-mean(slope_Out2_focal, na.rm = T)
slope_Out2_focal_sd<-stdev(slope_Out2_focal, na.rm = T)

Africa_sp<-world %>% dplyr::filter(continent == "Africa") %>% st_union()

slope_africa<-crop(slope_Out2_focal_mean, Africa_sp)
slope_africa<-mask(slope_africa, vect(Africa_sp))

slope_sd_africa<-crop(slope_Out2_focal_sd, Africa_sp)
slope_sd_africa<-mask(slope_sd_africa, vect(Africa_sp))

#tmap
library(tmaptools)
#tmaptools::palette_explorer()
data(World, land)
Africa_tm<-World %>% dplyr::filter(continent == "Africa")
US_tm<-World %>% dplyr::filter(name == "United States")

summary(slope_Out2_focal_sd*(10000))

slope_Out2_avg<-mean(slope_Out2*(10000), na.rm=T)

tmplot<-tm_shape(slope_Out2_focal*(10000), ylim = c(-60,88))+ #raster.warp = FALSE, projection="+proj=robin"
  #tm_graticules(alpha = 0.2, labels.size = 1) +
  tm_raster(title = expression(paste(gamma, " [10"^4*" kg/J]")),  #expression(paste(Delta, "SM"[AIF], "[10"^-5*" kg/J]"))
            legend.hist = F, style = "fixed", #"fixed"
            breaks = c(-Inf, -1.0, 0, 1.0, Inf),  #c(-Inf,-1,0,1,Inf)
            labels = c("< -1.00", "-1.00 to 0.00", "0.00 to 1.00", "> 1.00"),
            legend.show = T, 
            legend.hist.title = "Pixel Count", 
            palette = "PuOr",  #, palette = "-viridis", "-Spectral", "-RdBu", "BrBG"
            colorNA = NULL)+   
  tm_shape(World) +   #World, Africa_tm, US_tm, projection="+proj=robin"
  tm_borders(col = "black")+
  tm_layout(
    # inner.margins=c(.04,.03, .02, .01), outer.margins =c(0.25,0.25,0.25,0.25),
  #           earth.boundary = F, bg.color = "transparent",
  #           earth.datum = c(4326), space.color = "transparent",
  #           earth.boundary.color = "black",earth.boundary.lwd = 1.0,
            frame = F, 
            legend.show = T, legend.outside = T,
            legend.position = c("LEFT", "BOTTOM"),  #c("LEFT", "BOTTOM")
            legend.stack = "vertical",
            legend.hist.width = 0.25, legend.just = c("right","center"),
            legend.title.fontface = "bold",
            legend.title.size = 1.0, legend.text.size = 1.0,
            title.fontface = "bold",
            panel.labels = c("MAM", "JJA", "SON", "DJF"),
            panel.label.fontface = c("italic", "bold"),
            panel.label.size = 1.0,
            frame.lwd = 0,
            legend.format = options(digits=2))+
  tm_facets(ncol = 2, free.scales.raster = F)


tmap_save(tm = tmplot, 
          filename = "Median_Slope_between_WD_Pot_Grad.png",
          asp = 2.5, device = png, bg = "transparent",
          dpi = 600)

################################################################
#Time analysis
time_dimless1<-exp(-(WD_rast)/(abs(mu_seasonal)*slope_Out1))
time_dimless2<-exp(-(WD_rast)/(abs(mu_seasonal)*slope_Out2))
# time_dimless3<-exp(-(WD_rast^2)/(abs(mu_seasonal)*slope_Out1))
# time_dimless2<-exp(-(WD_rast^2)/(abs(mu_seasonal)*slope_Out2))
plot(time_dimless2)
summary(time_dimless1); summary(time_dimless2)

#Eliminate inf values
time_dimless1[time_dimless1 == Inf]<-NA
time_dimless2[time_dimless2 == Inf]<-NA


names(time_dimless1)<-c("tau_MAM", "tau_JJA", "tau_SON", "tau_DJF")
names(time_dimless2)<-c("tau_MAM", "tau_JJA", "tau_SON", "tau_DJF")

boxplot(time_dimless2, ymax = 1000)

hist(clamp(time_dimless2, 0, 100))
time_dimless2[time_dimless2 > 10000]<-NA

time_dimless2_focal = focal(time_dimless2, w = matrix(1, nrow = 7, ncol = 7), fun = median, na.rm=TRUE)
names(time_dimless2_focal)<-c("tau_MAM", "tau_JJA", "tau_SON", "tau_DJF")

summary(slope_Out2_focal); summary(time_dimless2_focal)
hist(clamp(time_dimless2_focal, 0, 100))
time_dimless2_focal[time_dimless2_focal > 1000]<-NA

writeCDF(time_dimless1, "Relative_Time_Dimensionaless_Using_Quantile_Regression_with_WD.nc")
writeCDF(time_dimless2, "Relative_Time_Dimensionaless_Using_Piecewise_Regression_with_WD.nc")

time_dimless2_focal_mean<-mean(time_dimless2_focal, na.rm=T)
time_dimless2_focal_sd<-stdev(time_dimless2_focal, na.rm = T)

time_dimless_Africa<-crop(time_dimless2_focal_mean, Africa_sp)
time_dimless_Africa<-mask(time_dimless_Africa, vect(Africa_sp))

time_dimless_sd_Africa<-crop(time_dimless2_focal_sd, Africa_sp)
time_dimless_sd_Africa<-mask(time_dimless_sd_Africa, vect(Africa_sp))

data("World")
library(tmaptools)
tmplot2<-tm_shape(time_dimless2_focal, ylim = c(-60,88) )+  #, ylim = c(-60,88), raster.warp = FALSE, projection="+proj=robin"
  #tm_grid(alpha = 0.2, labels.size = 1) +
  tm_raster(title = expression(paste(tau/tau0,  "[-]")),  #expression(paste(Delta, "SM"[AIF], "[10"^-5*" kg/J]"))
          legend.hist = F, 
          style = "quantile",   #"fixed"
          #breaks = c(0,0.1,1,10,Inf),  #c(-Inf,-1,0,1,Inf)
          #labels = c("0.00 to 0.10", "0.10 to 1.00", "1.00 to 10.00", "10.00 +"),
          legend.show = T, 
          legend.hist.title = "Pixel Count", 
          palette = get_brewer_pal("-PuOr", n = 5),  #, palette = "-viridis", "-Spectral", "-RdBu", "BrBG", "PuOr"
          colorNA = NULL)+   
  tm_shape(World) +   #World,Africa_tm , projection="+proj=robin"
  tm_borders(col = "black")+
  tm_layout(#scale = 1, 
    # inner.margins=c(.04,.03, .02, .01), outer.margins =c(0.25,0.25,0.25,0.25),
    # earth.boundary = F, bg.color = "transparent",
    # earth.datum = c(4326), space.color = "transparent",
    # earth.boundary.color = "black",earth.boundary.lwd = 1.0,
            frame = F, 
            legend.show = T,
            #legend.outside = T,
            #legend.position = c("LEFT", "BOTTOM"),  #c("LEFT", "BOTTOM")
            legend.stack = "vertical",
            #legend.hist.width = 0.25,
            legend.just = c("right","center"),
            legend.title.fontface = "bold",
            legend.title.size = 1.5,
            legend.text.size = 1.0,
            title.fontface = "bold",
            panel.labels = c("MAM", "JJA", "SON", "DJF"),
            panel.label.fontface = c("italic", "bold"),
            panel.label.size = 1.0,
            frame.lwd = 0,
            legend.format = options(digits=2) )+
  tm_facets(ncol = 2, free.scales.raster = F)

tmap_save(tm = tmplot2, 
          filename = "Seasonal_Relative_time.png",
          asp = 2.5, #device = png, bg = "transparent",
          dpi = 600)
