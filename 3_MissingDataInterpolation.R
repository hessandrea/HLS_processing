#this script reads the HLS index time series
#it fits an rlm tothe time series
#and then interpolates missing data based on the rlm 

rm(list = ls(all = TRUE))

#load libraries
#remove.packages("dplyr")
#install.packages('dplyr')
library(tidyverse)
library(raster)
library(lubridate)
library(MASS)
library(mgcv)

out_csv <- "./data_crunch/SatData/HLS/indices_robustified_bufferm_12knots.csv"
out_plot <-  "./docs/graphs/HLS/imputed_indices_10knots_buffer_byyear.png"
# load data ---------------------------------------------------------------
dat <- read_csv("./data_Crunch/SatData/HLS/HLS_time_series_buffer_cloudmasked.csv") %>%
  filter(value != "NA") %>%
  mutate(date = as.Date(obs_date), 
         trap = as.factor(trap))
head(dat)

dat_dup_means <- dat %>%
  group_by(trap, date, sensor, index) %>%
  summarize(value = mean(value, na.rm = TRUE))


dat_wide <- dat_dup_means %>%
  group_by(trap, date, sensor)%>%
  pivot_wider(names_from = index, values_from = value)
head(dat_wide)
names(dat_wide) <- c( "trap", "date","sensor","ndvi_buffer","ndwi_ga_buffer", "ndwi_mc_buffer", "ndwi_xu_buffer", "ndvi" ,       "ndwi_ga"   ,  "ndwi_mc"    , "ndwi_xu"    )



#plot data
ggplot(dat_wide) +
  geom_point(aes(x=date, y=ndwi_xu_buffer), color = "blue", size=0.5) +
  geom_point(aes(x=date, y=ndvi_buffer), color = "darkgreen", size=0.5) +
  facet_wrap(~trap)+
  xlab("")+
  ylab("NDVI_bufferzone and NDMI_xu_bufferzone")+
  ggtitle("Original values")+
  scale_x_date(date_labels = "%m-%Y")+
  theme(axis.text.x=element_text(angle=60, hjust=1))+
  ylim(-1, 1)


# create data frame of all possible trap and day combinations
str(unique(dat_wide$date))   #I have data on 131 days
all_dates <- seq(as.Date("2019/1/1"), as.Date("2020/12/31"), by = "day")
str(all_dates)   #I need data for all 1431 days

all_traps <- unique(dat_wide$trap)
str(all_traps)

trap_day_combinations <- data.frame(expand.grid(as.character(all_traps),all_dates))
names(trap_day_combinations) <- c("trap", "date")
trap_day_combinations <-trap_day_combinations %>%
  mutate(year = lubridate::year(date))
str(trap_day_combinations)


#join dataframe with all possible trap-day combinations to recorded data
#this will create na values on all days of the year where there is no sat data recorded
spectral_all_days <- right_join(dat_wide, trap_day_combinations, by = (c("trap", "date"))) 
summary(spectral_all_days)


#### running rlm by trap

plot_dat_template <- spectral_all_days %>%
  dplyr::select(-c("trap"))

smp_siz = floor(0.5*nrow(spectral_all_days))  # creates a value for dividing the data into train and test. In this case the value is defined as 75% of the number of rows in the dataset
smp_siz  # shows the value of the sample size
set.seed(123)   # set seed to ensure you always have same random numbers generated
train_ind = sample(seq_len(nrow(spectral_all_days)),size = smp_siz)  # Randomly identifies therows equal to sample size ( defined in previous instruction) from  all the rows of Smarket dataset and stores the row number in train_ind
train <- spectral_all_days[train_ind,] #creates the training dataset with row numbers stored in train_ind
test <- spectral_all_days[-train_ind,]  



bind_dat <- data.frame(matrix(ncol = 6, nrow = 0))
colnames(bind_dat) <- c("trap", "date", 
                        "imputed_ndvi", "imputed_ndvi_buffer", 
                        "imputed_ndwi_ga", "imputed_ndwi_ga_buffer")



# running the interpolation model -----------------------------------------

unique_traps <- unique(spectral_all_days$trap)
unique_years <- unique(year(spectral_all_days$date))

for (y in 1:length(unique_years)) {
  my_year = unique_years[y]
  spectral_year_sub <- spectral_all_days %>%
    filter(year(date)== my_year)
  
  for (i in 1:length(unique_traps)) {
    this_trap <- unique_traps[i]
    spectral_sub <- spectral_year_sub %>%
      filter(trap == this_trap) 
    
    
    spectral_sub$numdate <- as.numeric(spectral_sub$date)
    doy <- as.numeric(format(spectral_sub$date, "%j"))
    my_knots <- seq(min(doy),max(doy),length=10)
    my_spline <- cSplineDes(doy, my_knots)
    
    plot(doy,my_spline[,1],type="l")
    
    df <- data.frame(my_spline, spectral_sub$numdate) 
    new_df <- df %>%
      dplyr::rename(my_spline1 = X1, my_spline2 = X2,  my_spline3 = X3, my_spline4 = X4, my_spline5 = X5, my_spline6 = X6,
                    my_spline7 = X7, my_spline8 = X8,  my_spline9 = X9,
                    #number of my_splines is knot - 1
                    numdate = spectral_sub.numdate)
    
    ### running the rlm()
    ndvi_rlm <- rlm(ndvi ~ 0 + numdate + my_spline, data = spectral_sub)
    ndvi_buffer_rlm <- rlm(ndvi_buffer ~ 0 + numdate + my_spline, data = spectral_sub)
    ndwi_ga_rlm <- rlm(ndwi_ga ~ 0 + numdate + my_spline, data = spectral_sub)
    ndwi_ga_buffer_rlm <- rlm(ndwi_ga_buffer ~ 0 + numdate + my_spline, data = spectral_sub)

    
    
    #rlm(ndvi ~ numdate + my_spline[,1:2], data=spectral)
    
    ndvi_robustexpectation <- predict(ndvi_rlm, newdata = new_df)
    ndvi_buffer_robustexpectation <- predict(ndvi_buffer_rlm, newdata = new_df)
    ndwi_ga_robustexpectation <- predict(ndwi_ga_rlm, newdata = new_df)
    ndwi_ga_buffer_robustexpectation <- predict(ndwi_ga_buffer_rlm, newdata = new_df)
    
    
    #creatng robustscale for each index
    ndvi_robustscale <- ndvi_rlm$s
    ndvi_observed <- spectral_sub$ndvi
    ndvi_buffer_robustscale <- ndvi_buffer_rlm$s
    ndvi_buffer_observed <- spectral_sub$ndvi_buffer
    
    ndwi_ga_robustscale <- ndwi_ga_rlm$s
    ndwi_ga_observed <- spectral_sub$ndwi_ga
    ndwi_ga_buffer_robustscale <- ndwi_ga_buffer_rlm$s
    ndwi_ga_buffer_observed <- spectral_sub$ndwi_ga_buffer

    
  
    # filling in actual observations into robustified data frame
    ndvi_robustifiedobservation <- ndvi_observed
    ndvi_buffer_robustifiedobservation <- ndvi_buffer_observed
    ndwi_ga_robustifiedobservation <- ndwi_ga_observed
    ndwi_ga_buffer_robustifiedobservation <- ndwi_ga_buffer_observed

    
    #calculating robust score for each index
    ndvi_robustzscore <- (ndvi_observed - ndvi_robustexpectation) / ndvi_robustscale
    ndvi_buffer_robustzscore <- (ndvi_buffer_observed - ndvi_buffer_robustexpectation) / ndvi_buffer_robustscale
    ndwi_ga_robustzscore <- (ndwi_ga_observed - ndwi_ga_robustexpectation) / ndwi_ga_robustscale
    ndwi_ga_buffer_robustzscore <- (ndwi_ga_buffer_observed - ndwi_ga_buffer_robustexpectation) / ndwi_ga_buffer_robustscale

    
    
    #robust_df <- data.frame(ndvi_observed, ndvi_robustzscore,  ndvi_robustexpectation, 
    #                        ndmi_observed, ndmi_robustzscore,  ndmi_robustexpectation)
    
    #filtering out index observations with a robustscore above or below 0
    ndvi_robustifiedobservation[ndvi_robustzscore >= 3] <- NA
    ndvi_robustifiedobservation[ndvi_robustzscore <= -3] <- NA
    ndvi_robustifiedobservation[is.na(ndvi_robustifiedobservation)] <- NA
    ndvi_buffer_robustifiedobservation[ndvi_buffer_robustzscore >= 3] <- NA
    ndvi_buffer_robustifiedobservation[ndvi_buffer_robustzscore <= -3] <- NA
    ndvi_buffer_robustifiedobservation[is.na(ndvi_buffer_robustifiedobservation)] <- NA
    
    ndwi_ga_robustifiedobservation[ndwi_ga_robustzscore >= 3] <- NA
    ndwi_ga_robustifiedobservation[ndwi_ga_robustzscore <= -3] <- NA
    ndwi_ga_robustifiedobservation[is.na(ndwi_ga_robustifiedobservation)] <- NA
    ndwi_ga_buffer_robustifiedobservation[ndwi_ga_buffer_robustzscore >= 3] <- NA
    ndwi_ga_buffer_robustifiedobservation[ndwi_ga_buffer_robustzscore <= -3] <- NA
    ndwi_ga_buffer_robustifiedobservation[is.na(ndwi_ga_buffer_robustifiedobservation)] <- NA

    
    
    robust_df <- data.frame(ndvi_robustifiedobservation, ndvi_buffer_robustifiedobservation, 
                            ndwi_ga_robustifiedobservation, ndwi_ga_buffer_robustifiedobservation)
    
    #robust_df <- robust_df %>%
    #  mutate(imputed_ndvi = ifelse(ndvi_robustzscore <= 3,  ndvi_observed, NA), #
    #         imputed_ndvi = ifelse(ndvi_robustzscore >= -3, ndvi_observed, NA), #n
    #         imputed_ndmi = ifelse(ndmi_robustzscore <= 3, ndmi_observed, NA), #
    #         imputed_ndmi = ifelse(ndmi_robustzscore >= -3, ndmi_observed, NA)) #
    
    
    robust_df <- within(robust_df, ndvi_robustifiedobservation[is.na(ndvi_robustifiedobservation)] <- ndvi_robustexpectation[is.na(ndvi_robustifiedobservation)])
    robust_df <- within(robust_df, ndvi_buffer_robustifiedobservation[is.na(ndvi_buffer_robustifiedobservation)] <- ndvi_buffer_robustexpectation[is.na(ndvi_buffer_robustifiedobservation)])
    
    robust_df <- within(robust_df, ndwi_ga_robustifiedobservation[is.na(ndwi_ga_robustifiedobservation)] <- ndwi_ga_robustexpectation[is.na(ndwi_ga_robustifiedobservation)])
    robust_df <- within(robust_df, ndwi_ga_buffer_robustifiedobservation[is.na(ndwi_ga_buffer_robustifiedobservation)] <- ndwi_ga_buffer_robustexpectation[is.na(ndwi_ga_buffer_robustifiedobservation)])
    

    summary(robust_df)
    str(robust_df)
    
    robust_df$trap <- this_trap
    robust_df$date <- spectral_sub$date
    robust_df <- robust_df %>%
      dplyr::select(c("trap", "date",
                      "ndvi_robustifiedobservation", "ndvi_buffer_robustifiedobservation", 
                      "ndwi_ga_robustifiedobservation", "ndwi_ga_buffer_robustifiedobservation"))
    
    bind_dat <- rbind(bind_dat, robust_df)
  }
}

unique(bind_dat$trap)
str(bind_dat$date)
str(bind_dat)

plot_dat <- merge(plot_dat_template, bind_dat, by = c("trap", "date")) %>%
  filter(between(date, as.Date("2019-03-01"), as.Date("2019-11-30")) |
           between(date, as.Date("2020-03-01"), as.Date("2020-11-30")))
summary(plot_dat)

par(mfrow = c(1,1))

ggplot(plot_dat) +
  geom_point(aes(x=date, y=ndvi_buffer_robustifiedobservation, color = "darkgreen"), size=0.5) +
  geom_point(aes(x=date, y=ndwi_ga_buffer_robustifiedobservation, color = "darkgoldenrod"), size=0.5) +
  facet_wrap(~trap)+
  xlab("")+
  ylab("indices")+
  ggtitle("HLS data, imputed to daily, 12knots")+
  scale_x_date(date_labels = "%m-%Y")+
  theme(axis.text.x=element_text(angle=60, hjust=1))+
  scale_colour_manual(name = 'indices', 
                      values =c("blue"="blue", 'darkgoldenrod'='darkgoldenrod','darkgreen'='darkgreen'), 
                      labels = c('NDWI','NDMI', 'NDVI'))+
  guides(color = guide_legend(override.aes = list(size = 2)))

ggsave(filename = out_plot, plot=last_plot(), 
       width = 30, height = 25, unit = "cm")


out_dat <- plot_dat %>%
  mutate(doy = as.numeric(format(plot_dat$date, "%j")))

write.csv(out_dat, out_csv)
