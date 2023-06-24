
library(tidyverse)

#------ Set Plotting theme----------------
gg_theme <- function () { 
  theme_bw() %+replace% 
    theme(
      axis.text = element_text(colour = "black", size = 11),
      axis.title = element_text(size=13),
      axis.ticks = element_line(colour = "black"),
      panel.grid = element_blank(),
      strip.background = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA),
      axis.line = element_line(colour = "black")
      
    )
}


library(rsoi)

##################################################
# Download Antarctic Oscillation (= SAM) data
##################################################
# Description
# Projection of the monthly 700 hPa anomaly height field south of 20?S on the first EOF obtained
# from the monthly 700 hPa height anomaly.
# References
# https://www.cpc.ncep.noaa.gov/products/precip/CWlink/daily_ao_index/aao/aao.shtml

aao = download_aao(use_cache = FALSE)

#################################################
# Download Southern Oscillation Index data
#################################################

# Description
# The Southern Oscillation Index is defined as the standardized difference between barometric readings at Darwin, Australia and Tahiti. The Oceanic Nino Index is average sea surface temperature
# in the Nino 3.4 region (120W to 170W) averaged over three months. Phases are categorized by
# Oceanic Nino Index:
# Month: Month of record
# Year: Year of record?
# SOI: Southern Oscillation Index
# SOI_3MON_AVG: 3 Month Average Southern Oscillation Index
#References
# https://www.ncdc.noaa.gov/teleconnections/enso/indicators/soi/

soi = download_soi(use_cache = FALSE)

#-------------------------------------------------------------------------------
# Data formatting 
#-------------------------------------------------------------------------------

# SOI
# Select data from March to September (pre-breeding foraging trip) only, for all years
min(soi$Year)
soi.ses = subset(soi, Year < 2020)
soi.ses = subset(soi.ses, Month != "Jan")
soi.ses = subset(soi.ses, Month != "Feb")
soi.ses = subset(soi.ses, Month != "Oct")
soi.ses = subset(soi.ses, Month != "Nov")
soi.ses = subset(soi.ses, Month != "Dec")
soi.ses 

soi.ses.year =  soi.ses %>%
  group_by(Year) %>%
  dplyr::summarise(Mean = mean(SOI, na.rm=TRUE))

## Save Plot 
 pdf("./plots/FigureS5.pdf",
    useDingbats = FALSE, width = 10, height = 7)
plot(soi.ses.year$Year, soi.ses.year$Mean, type = "b", xlab = "Year", ylab = "SOI Index", pch = 16) 
abline(h = 0)
dev.off()


# SAM
# Select data from March to September (pre-breeding foraging trip) only, for study years
min(aao$Year)
sam.ses = subset(aao , Year < 2020)
sam.ses = subset(sam.ses, Month != "Jan")
sam.ses = subset(sam.ses, Month != "Feb")
sam.ses = subset(sam.ses, Month != "Oct")
sam.ses = subset(sam.ses, Month != "Nov")
sam.ses = subset(sam.ses, Month != "Dec")
sam.ses 
sam.ses.year =  sam.ses %>%
  group_by(Year) %>%
  dplyr::summarise(Mean = mean(AAO, na.rm=TRUE))

## Save Plot 
pdf("./plots/FigureS6.pdf",
     useDingbats = FALSE, width = 10, height = 7)
plot(sam.ses.year$Year, sam.ses.year$Mean, type = "b", xlab = "Year", ylab = "SAM Index", pch = 16)
abline(h = 0)
dev.off()

# add female population size on 15 October
popN <- read.csv("./data/ses_15_oct_counts.csv")  # import total island count data
head(popN)

## Save Plot 
pdf("./plots/FigureS4.pdf",
     useDingbats = FALSE, width = 10, height = 7)
plot(popN$year, popN$N, type = "b", xlab = "Year", ylab = "Peak breeding season count", pch = 16)
abline(h = 0)
dev.off()

