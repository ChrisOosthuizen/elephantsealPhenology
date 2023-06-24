
# Mixed model analysis of female elephant seal breeding (arrival) date 
# from mark recapture observations at Marion Island.
# Proceedings of the Royal Society B
# DOI: 10.1098/rspb.2023-1170
# Reproductive phenology is a repeatable, heritable trait linked to the 
# timing of other life history events in a migratory marine predator

# Chris Oosthuizen
# June 2023

# ------------------------------------------------------------------------------------------
# Load packages
#------------------------------------------------------------------------------------------

library(tidyverse)
library(lubridate)
library(lme4)
library(broom)
library(broom.mixed)
library(chron)
library(sjPlot)
library(ggeffects)  

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

# ------------------------------------------------------------------------------------------
# Import cleaned census data
# ------------------------------------------------------------------------------------------

# import the first breeding observation (social class 9/3) for a seal
tags = readRDS("./data/tagdat.rds")
tags

# Get year of first breeding
tags =  tags %>%
  group_by(ID) %>%
  mutate(first.breeding.year = year[1])%>% 
  ungroup()  

tags$AFB = tags$first.breeding.year - tags$birthyear

# Is the animal breeding for the first time (recruitment) or not?
# order bY date, ID
tags = tags[order(tags$Date, tags$ID),]

# This counts 1 value for every new individual; the value does not increase for old individuals
tags = tags %>% 
      mutate(new.indiv = cumsum(!duplicated(ID)))

# This tells you if this is the first time you see an animal in a breeding season, or a repeat 
tags <- tags %>%
  mutate(first.record = if_else(duplicated(new.indiv) == FALSE, 1, 0))

# Assign variable state to FTB or B
tags <- tags %>%
  group_by(ID, year) %>%
  mutate(state = case_when(any(first.record == 1) ~ "FTB",
                                any(first.record == 0) ~ "B")) %>%
  ungroup() 

# order by date
tags = tags[order(tags$ID, tags$Date),]

#Add the first breeding date of a seal
tags =  tags %>%
  group_by(ID) %>%
  mutate(first.breeding.date = julian[1]) %>%   
  ungroup()

# sanity check: look at 'age'
table(tags$state, tags$Age)  # field age (contain some errors)
table(tags$state, tags$age)  # calculated age (good)

## clean up slightly 
tags = tags %>% dplyr::select(-c(Age, new.indiv, first.record))
tags

# Remove unrealistic outliers: 
# first.breeding.date of RR(1)478 was 319 (far outlier to the right). Exclude this animal
tags = subset(tags, ID != "RR(1)478")
# also exclude records before September (unrealistic)
tags = subset(tags, julian > 245)
# also exclude a record from 14 December (unrealistic)
tags = subset(tags, julian < 348)

# Look at sample size per age
tapply(tags$julian, tags$age, length)  

# Make and age class with maximum age == 21
tags$ageclass = pmin(tags$age, 21)
tapply(tags$julian, tags$ageclass, length)  
tags$ageclass.factor = as.factor(tags$ageclass)

# Look at sample size per season
tapply(tags$julian, tags$year, length)  
# Sample size in 1980s is low. Assess over the period 1989-2019
tags = subset(tags, year >1988)
tapply(tags$julian, tags$year, length)  

# add SAM and SOI environmental covariates
covar <- readRDS("./data/SOI_SAM_covar.RDS")   # import SAM and SOI values per year (averages) obtained from rsoi package
tags = left_join(tags, covar, by = "year")  # merge: left_join keeps all rows of tags, and matches 'year' from tags and covar, adding 2 columns

# add female population size on 15 October
popN <- read.csv("./data/ses_15_oct_counts.csv")  # import total island count data
popN$Nstd = (popN$N-(mean(popN$N)))/sd(popN$N)  # standardize to mean 0, sd 1

## clean up slightly 
popN = popN %>% dplyr::select(c(year,Nstd))
popN = subset(popN , year >1988)

tags = left_join(tags, popN, by = "year")  # merge: left_join keeps all rows of tags, and matches 'year' from tags and popN

# Creat age * state interaction variable (don't do in in the model, as there are not data for each level (no old FTB))
tags = tags %>%
  mutate(age.state = ageclass) %>%
  mutate(age.state = case_when(
    .$ageclass == "4" & .$state == "FTB" ~ 104,
    .$ageclass == "5" & .$state == "FTB" ~ 105,
    .$ageclass == "6" & .$state == "FTB" ~ 106,
    .$ageclass == "7" & .$state == "FTB" ~ 107,
    TRUE ~ as.numeric(as.character(.$ageclass))))

# Update the numerical codes from above code to more meaningful text codes 
tags$age.state = factor(tags$age.state, 
                             levels = c(3:21, 104:107),  # these are the numerical codes above
                             labels = c(3:21, "4FTB","5FTB","6FTB","7FTB"))
# all age 3 animals are FTB, so we have an age.state = 3, and do not need an additional age.state = 3FTB

# Check data
table(tags$state, tags$ageclass)  # calculated age (good)
table(tags$state, tags$age.state)  # age * state interaction 

# exclude 1998 due to low observation effort
tags = subset(tags, year != 1998)

#-----------------------------------------------------------------------------------------
# Plot 'discovery curve' across all years
# Colours of the discovery curve represent whether it is a new or repeat observation in 
# every breeding season 
#-----------------------------------------------------------------------------------------

# First group data by year
tags <- group_by(tags, year) %>%
  mutate(new.indiv.year = cumsum(!duplicated(ID))) %>%
  mutate(first.record.year = if_else(duplicated(new.indiv.year) == FALSE, 1, 0))

# order by date
tags = tags[order(tags$Date, tags$ID),]

# How many observations per individual per breeding season? (tags was GROUPED above, by tags,year)
tags %>% count(ID)

# Add discovery.curve column: This counts 1 value for every new individual; the value does not increase for old individuals
tags = tags %>% 
  mutate(discovery.curve = cumsum(!duplicated(ID)))

# Add first.record.year column: this tells you if this is the first time you see an animal in a breeding season, or a repeat 
tags <- tags %>%
  mutate(first.record.year = if_else(duplicated(discovery.curve) == FALSE, 1, 0))

# Plot discovery curve with lines per year
dc = ggplot(tags, aes(x=julian, y= discovery.curve))+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 8),limits = c(0, 300))+
  scale_x_continuous(breaks = seq(250,330, by = 10),limits = c(245, 330),
                     labels =c('07 Sep','17 Sep','27 Sep','07 Oct','17 Oct','27 Oct','06 Nov','16 Nov','26 Nov'))+
  gg_theme()+
  xlab("Date") +
  ylab("Total number of individuals identified")+
  theme(axis.text.x = element_text(angle = -90, vjust = 0.5, hjust = 0))+
  geom_line(aes(group = as.factor(year), col = as.factor(year)))
#  geom_line(aes(group = as.factor(year)), color = 'black')
  
dc

# You have to ungroup 'tags', otherwise the following code does not work.
tags = tags %>%
  ungroup() 


#-----------------------------------------------------------------------------------------
#Lmer Model selection for breeding date
#-----------------------------------------------------------------------------------------

# select only the first record per year per breeding season
tags1 = subset(tags, first.record.year == 1)

# How many observations per individual per year (max to min sorted)?
tags1 %>% count(ID, year, sort = TRUE)

# Add column that gives the number of years that an individual is in the data
tags1 = tags1 %>% 
  group_by(ID) %>%
  mutate(n.years = n_distinct(year))

# Check that column n is correct (use summarise - i.e. one row per individual)
n.breeding = tags %>% 
  group_by(ID) %>%
  summarise(n = n_distinct(year))%>%
  ungroup() 

arrange(n.breeding, desc(n)) # Yes, all the old individuals are breeding ~ 15 times
# summarise - this is not the data used in the paper, so will differ at this point
mean(n.breeding$n)
sum(n.breeding$n)
n_distinct(n.breeding$ID)

# average breeding attempts if breeding more than once:
n.breeding2 = subset(n.breeding, n.breeding$n > 1)
mean(n.breeding2$n)

# Now one can subset the data to include only females with multiple breeding occasions
# Don't do this subsetting - e.g., Martin et al 2011, Clermont et al 2018) 
# tags1 = subset(tags1, n.years > 1)

#------------------------------------------------
#### Covariates - within-subject centering ####
#------------------------------------------------
# 'within-subject centering' of covariates (Van de Pol & Wright 2009 Animal Behaviour 77:753-758)
# do this for tags1, otherwise the means may differ (depends on the input data;
# tags and tags1 will have different means as the observations differ.)

#--------
# SOI
#--------
# The between-individual effect is fitted as the average x value for each individual
tags1 =  tags1 %>%
  group_by(ID) %>%
  mutate(aveSOI = mean(SOI))
# The within individual effect is fitted as the deviation from the individual mean 
# for a given observation  (x - mean(x)) 
tags1$deltaSOI = tags1$SOI - tags1$aveSOI

#--------
# SAM
#--------
# The between-individual effect is fitted as the average x value for each individual
tags1 =  tags1 %>%
  group_by(ID) %>%
  mutate(aveSAM = mean(SAM))
# The within individual effect is fitted as the deviation from the individual mean 
# for a given observation  (x - mean(x)) 
tags1$deltaSAM = tags1$SAM - tags1$aveSAM

#-------------------------
# Population density
#-------------------------
# The between-individual effect is fitted as the average x value for each individual
tags1 =  tags1 %>%
  group_by(ID) %>%
  mutate(aveNstd = mean(Nstd))
# The within individual effect is fitted as the deviation from the individual mean 
# for a given observation  (x - mean(x)) 
tags1$deltaNstd = tags1$Nstd - tags1$aveNstd

# Wean date has only 1 value per individual - don't fit it as a within-between covariate.

#-------------------------------------------------------------------------------
# Moult
#-------------------------------------------------------------------------------
# import moult observations

moultdat = readRDS("./data/moultdat.rds")

# add sealyear
moultdat = moultdat %>% mutate(sealyear = ifelse(julian > 200, year, year-1))

moultdat$sealyear4B = moultdat$sealyear+1

moultdat$julian365 = moultdat$julian + 365 
moultdat$moultjulian = moultdat$julian  

moultdat = moultdat %>%
  mutate(
    moultjulian = replace(moultjulian,       ## mutate (and replace) the column julian.end
                          moultjulian  < 200,    # where this is the case (> 200 days difference)
                          julian365[julian  < 200]))   # if the above is TRUE, replace julian.start into the new julian.start column

# Make it days since 1 November 
moultdat$moultjulian = moultdat$moultjulian - 305    # 305 is 1 November in julian days 

plot(moultdat$moultjulian)  # looks OK (some records before 1 Nov, which is now day 0)

# Select the last moult entry for every seal, for every year
moultdat = moultdat %>%
  group_by(ID, sealyear4B) %>%
  slice(n()) %>%
  ungroup() 

plot(moultdat$moultjulian)  # looks OK (some records before 1 Nov, which is now day 0 are unlikely 
                            # but I will keep them here)

#----------------------------------------------------------------------------
# Now merge with tags1 data
#----------------------------------------------------------------------------
names(moultdat)

## clean up slightly 
moultdat = moultdat %>% dplyr::select(c(Date, ID, julian, sealyear4B, moultjulian))
names(moultdat) = c("moultdate", "ID", 'moultjulian', 'sealyear4B', 'moultjulian1NOV')

# Now merge
tags1moult <- merge(x = moultdat, y = tags1,  by.x=c("ID", "sealyear4B"), by.y=c("ID", "year"))
dim(tags1moult)

# remove a moult outlier: we know this specific female FB (2008) moult at any time of the year:
tmp = subset(tags1moult, tags1moult$birthyear == 2008)
plot(tmp$moultjulian)  # looks OK (some records before 1 Nov, which is now day 0)
# view(tmp)  # FB253

tags1moult = subset(tags1moult, ID != "FB(1)253")

#--------------------------------------------------------
# Summary statistics
#--------------------------------------------------------
dim(tags1moult)
n_distinct(tags1moult$ID) # how many mothers
tags1moult$tally = paste0(tags1moult$ID, tags1moult$year.random)
n_distinct(tags1moult$tally) # how many breeding attempts
mean(tags1moult$n.years)
sd(tags1moult$n.years)
range(tags1moult$n.years)

quantile(tags1moult$julian)
quantile(tags1moult$julian, c(0.025, 0.975)) # 95% 
quantile(tags1moult$julian, c(0.1, 0.9))   # 80% 

subset(tags1moult, tags1moult$julian == 266)
subset(tags1moult, tags1moult$julian == 292)
subset(tags1moult, tags1moult$julian == 277)  # median

#===============================================================================================
# Model selection with centering
#===============================================================================================

# wean date was not centered, as it does not vary with time (ave wean = wean ij, thus delta wean = 0)

# Z-standardize some variables first (avoid convergence errors with random slope models)
# moult date (julian)
tags1moult$moultstd = datawizard::standardize(tags1moult$moultjulian1NOV)
# wean date (julian)
tags1moult$weanstd = datawizard::standardize(tags1moult$wean.date)
# time trend (years)
tags1moult$lineartrend = datawizard::standardize(as.numeric(tags1moult$year.random))

#-----------------------------------------
# Now draw Median arrival date histogram
#-----------------------------------------
# Select any year
tags1moulth = subset(tags1moult, tags1moult$year == 2009)

# order by date
tags1moulth = tags1moulth[order(tags1moulth$Date, tags1moulth$ID),]

# How many observations per individual in this breeding season?
tags1moulth %>% count(ID)

# Add discovery.curve column: This counts 1 value for every new individual; the value does not increase for old individuals
tags1moulth = tags1moulth %>% 
   mutate(discovery.curve = cumsum(!duplicated(ID)))

# Add first.record.year column: this tells you if this is the first time you see an animal in a breeding season, or a repeat 
tags1moulth <- tags1moulth %>%
  mutate(first.record.year = if_else(duplicated(discovery.curve) == FALSE, 1, 0))

# find INDEX (position) of minimum entry in column discovery.curve:
idx <- match(median(tags1moulth$discovery.curve), tags1moulth$discovery.curve) # gives you the first entry of the median
idx

# if does not fall on an even number (median(dat1$discovery.curve)), this works:
# idx <- match(median(tags1moulth$discovery.curve)+0.5, tags1moulth$discovery.curve) # gives you the first entry of the median
# idx

idx2 <- which(idx %in% tags1moulth$discovery.curve) # gives you all indices of the median
#idx2 <- which(median(tags1moulth$discovery.curve) %in% tags1moulth$discovery.curve) # gives you all indices of the median
idx2  # should be 1

# What is the Date where [idx] is (ie. the date of the median in discover.curve column:
tags1moulth$Date[idx]

# median dates
# Subtracted a day fom leap years so that all julian days are the same for the breeding season period 

median.arrival.date = data.frame(median = c(
"1989-10-05 UTC",
"1990-10-01 UTC",
"1991-10-04 UTC",
"1992-10-02 UTC",#-1
"1993-10-02 UTC",
"1994-10-02 UTC",
"1995-10-06 UTC",
"1996-10-07 UTC",#-1
"1997-10-03 UTC",
"1999-10-06 UTC",
"2000-10-04 UTC",#-1 
"2001-10-07 UTC",
"2002-10-04 UTC",
"2003-10-04 UTC",
"2004-10-01 UTC",#-1
"2005-10-06 UTC",
"2006-10-04 UTC",
"2007-10-03 UTC",
"2008-10-03 UTC",#-1
"2009-10-06 UTC",
"2010-10-05 UTC",
"2011-10-05 UTC",
"2012-10-02 UTC",#-1 
"2013-10-04 UTC",
"2014-10-05 UTC",
"2015-10-02 UTC",
"2016-10-04 UTC", #-1 
"2017-10-04 UTC",
"2018-10-04 UTC",
"2019-10-05 UTC"))

# format date as POSIXct
median.arrival.date$dt = strptime(as.character(median.arrival.date$median),"%Y-%m-%d") 

# add julian date with lubridate
median.arrival.date$julian <- yday(median.arrival.date$dt) 

# Now plot a histogram. Important, set the same xlim values here as for the plot you want to add
# this histogram to.
xhist.arrival <- ggplot(median.arrival.date, aes(julian)) +     # data = summary of medians, not raw data
  geom_histogram(binwidth = 1)+
  xlim(c(245, 330))+
  gg_theme()

library(cowplot)
xhist.arrival
p2 <- insert_xaxis_grob(dc, xhist.arrival, grid::unit(1, "in"), position = "top")
ggdraw(p2)

## Save Plot 
 pdf("./plots/FigureS3.pdf",
     useDingbats = FALSE, width = 6, height = 5)
 ggdraw(p2)
 dev.off()

#-----------------------------------------
##### Model covariates ####
#-----------------------------------------

names(tags1moult)

# First we have a look at the correlations between covariates:
cov = tags1moult %>% dplyr::select(SOI, SAM, Nstd, lineartrend, moultstd, wean.date, age)
names(cov) = c('SOI', 'SAM', 'Pop.N', 'Linear.trend', 'Moult.date', 'Wean.date', 'Age')

GGally::ggcorr(cov, nbreaks = 10, label = T, low = "red3", high = "green3", method = c("pairwise", "pearson"),
       label_round = 2, name = "Correlation Scale", label_alpha = F, hjust = 0.75, label_size = 6,
       legend.size = 16, size = 6)


## Save Plot 
pdf("./plots/FigureS7.pdf",
     useDingbats = FALSE, width = 9, height = 9)
GGally::ggcorr(cov, nbreaks = 10, label = T, low = "red3", high = "green3", method = c("pairwise", "pearson"),
               label_round = 2, name = "Correlation Scale", label_alpha = F, hjust = 0.75, label_size = 6,
               legend.size = 16, size = 6)
dev.off()


#-----------------------------------------
##### Model selection ####
#-----------------------------------------

#-----------------------------------------
# Model set 1
# All fixed covariates
# Random intercepts for ID and year.
# Select random slopes on ID
#-----------------------------------------

full.interc <- lmer(julian ~  age.state + moultstd + weanstd +
                     aveSAM + deltaSAM +
                     aveSOI + deltaSOI + 
                     aveNstd + deltaNstd +
                     lineartrend + 
                     (1|ID) + (1|year.random), tags1moult, REML = T)

full.slope.moult <- lmer(julian ~ age.state + moultstd + weanstd +
                     aveSAM + deltaSAM +
                     aveSOI + deltaSOI + 
                     aveNstd + deltaNstd +
                     lineartrend + 
                     (moultstd|ID) + (1|year.random), tags1moult, REML = T)

# Pers. comm. Niels Dingemanse:
# This fits SAM while controlling for differences in mean SAM:
full.slope.SAM_mean <- lmer(julian ~ age.state + moultstd + weanstd +
                     aveSAM + SAM +
                     aveSOI + deltaSOI + 
                     aveNstd + deltaNstd +
                     lineartrend + 
                     (SAM|ID) + (1|year.random), tags1moult, REML = T)
# singular

# Pers. comm. Niels Dingemanse:
# This fits the random slopes of SAM *within* individuals. A problem is that if the effect
# of SAM is non-linear, there will be the appearance of random slopes where none exists. 
# See the Suppl Mat of Dingemanse & Dochtermann 2013, JAE).

full.slope.SAM_within <- lmer(julian ~ age.state + moultstd + weanstd +
                     aveSAM + deltaSAM +
                     aveSOI + deltaSOI + 
                     aveNstd + deltaNstd +
                     lineartrend + 
                     (deltaSAM|ID) + (1|year.random), tags1moult, REML = T)
# singular

# Pers. comm. Niels Dingemanse:
# This fits SAM while controlling for differences in mean SAM:
full.slope.SOI_mean <- lmer(julian ~ age.state + moultstd + weanstd +
                     aveSAM + deltaSAM +
                     aveSOI + SOI + 
                     aveNstd + deltaNstd +
                     lineartrend + 
                     (SOI|ID) + (1|year.random), tags1moult, REML = T)
# singular

# Pers. comm. Niels Dingemanse:
# This fits the random slopes of SOI *within* individuals. 
full.slope.SOI_within <- lmer(julian ~ age.state + moultstd + weanstd +
                     aveSAM + deltaSAM +
                     aveSOI + deltaSOI + 
                     aveNstd + deltaNstd +
                     lineartrend + 
                     (deltaSOI|ID) + (1|year.random), tags1moult, REML = T)

full.slope.Nstd_mean <- lmer(julian ~ age.state + moultstd + weanstd +
                     aveSAM + deltaSAM +
                     aveSOI + deltaSOI + 
                     aveNstd + Nstd +
                     lineartrend + 
                     (Nstd|ID) + (1|year.random), tags1moult, REML = T)
# singular

full.slope.Nstd_within <- lmer(julian ~ age.state + moultstd + weanstd +
                     aveSAM + deltaSAM +
                     aveSOI + deltaSOI + 
                     aveNstd + deltaNstd +
                     lineartrend + 
                     (deltaNstd|ID) + (1|year.random), tags1moult, REML = T,
                      control = lmerControl(
                           optimizer ='optimx', optCtrl=list(method='L-BFGS-B')))



AIC(full.interc , full.slope.moult, 
    full.slope.SAM_mean, full.slope.SAM_within, 
    full.slope.SOI_mean, full.slope.SOI_within,
    full.slope.Nstd_mean, full.slope.Nstd_within)

# https://m-clark.github.io/posts/2020-03-16-convergence/
# library(mixedup)
# summarize_model(full.slope.moult, ci = FALSE)

#----------------------------------
# Test convergence with different optimizers
#----------------------------------
library(optimx)

full.slope.SAM_within_all = allFit(full.slope.SAM_within)
full.slope.SOI_within_all = allFit(full.slope.SOI_within)
full.slope.Nstd_within_all = allFit(full.slope.Nstd_within)

full.slope.SAM_summary = summary(full.slope.SAM_within_all)
full.slope.SOI_summary = summary(full.slope.SOI_within_all)
full.slope.Nstd_summary = summary(full.slope.Nstd_within_all)

full.slope.SAM_summary
full.slope.SOI_summary
full.slope.Nstd_summary

# Should do likelihood ratio tests: 
# Models are being fitted with ML
# does inclusion of random slope for moult improves model?
anova(full.interc, full.slope.moult)    # Yes
AIC(full.interc) - AIC(full.slope.moult)

# does inclusion of random slope for SOI improves model?
anova(full.interc, full.slope.SOI_within)      # No
anova(full.interc, full.slope.SOI_mean)      # No
AIC(full.slope.SOI_within) - AIC(full.slope.moult) 

# does inclusion of random slope for SAM improves model?
anova(full.interc, full.slope.SAM_within)   # No
anova(full.interc, full.slope.SAM_mean)     # No
AIC(full.slope.SAM_mean) - AIC(full.slope.moult) 

# does inclusion of random slope for N improves model?
anova(full.interc, full.slope.Nstd_mean)    # Almost relative to null, but model is singular.
anova(full.interc, full.slope.Nstd_within)  # No
AIC(full.slope.Nstd_within) - AIC(full.slope.moult) 

# Fit of moult slope model vs Nstd_mean model: 
anova(full.slope.moult, full.slope.Nstd_mean)  # No 

#--------------------------------------------------------------------- 
# Continue with model with random slope for moult (full.slope.moult)
# Now do model selection on fixed effects.
#---------------------------------------------------------------------

full.slope.moult <-  lmer(julian ~ age.state + moultstd + weanstd +
                     aveSAM + deltaSAM +
                     aveSOI + deltaSOI + 
                     aveNstd + deltaNstd +
                     lineartrend + 
                     (moultstd|ID) + (1|year.random), tags1moult, REML = F)

nowean <-  lmer(julian ~ age.state + moultstd + 
                     aveSAM + deltaSAM +
                     aveSOI + deltaSOI + 
                     aveNstd + deltaNstd +
                     lineartrend + 
                     (moultstd|ID) + (1|year.random), tags1moult, REML = F)


noagestate <-  lmer(julian ~ moultstd + weanstd +
                     aveSAM + deltaSAM +
                     aveSOI + deltaSOI + 
                     aveNstd + deltaNstd +
                     lineartrend + 
                     (moultstd|ID) + (1|year.random), tags1moult, REML = F)

nomoult <-  lmer(julian ~ age.state + weanstd +
                     aveSAM + deltaSAM +
                     aveSOI + deltaSOI + 
                     aveNstd + deltaNstd +
                     lineartrend + 
                     (moultstd|ID) + (1|year.random), tags1moult, REML = F)

null <- lmer(julian ~ 1 + (moultstd|ID)+ (1|year.random), tags1moult, REML = F)


state <-  lmer(julian ~ state + moultstd + weanstd +
                     aveSAM + deltaSAM +
                     aveSOI + deltaSOI + 
                     aveNstd + deltaNstd +
                     lineartrend + 
                     (moultstd|ID) + (1|year.random), tags1moult, REML = F)


ageF <-  lmer(julian ~ ageclass.factor + moultstd + weanstd +
                     aveSAM + deltaSAM +
                     aveSOI + deltaSOI + 
                     aveNstd + deltaNstd +
                     lineartrend + 
                     (moultstd|ID) + (1|year.random), tags1moult, REML = F)

AIC(full.slope.moult, nowean, noagestate, nomoult, null, state, ageF)

# Removing wean date, moult and age or breeding state resulted in large increase in AIC,
# indicating the importance of these variables.

noSAM <-  lmer(julian ~ age.state + moultstd + weanstd +
                   #  aveSAM + deltaSAM +
                     aveSOI + deltaSOI + 
                     aveNstd + deltaNstd +
                     lineartrend + 
                     (moultstd|ID) + (1|year.random), tags1moult, REML = F)

noSOI <-  lmer(julian ~ age.state + moultstd + weanstd +
                     aveSAM + deltaSAM +
                   #  aveSOI + deltaSOI + 
                     aveNstd + deltaNstd +
                     lineartrend + 
                     (moultstd|ID) + (1|year.random), tags1moult, REML = F)

noNstd <-  lmer(julian ~ age.state + moultstd + weanstd +
                     aveSAM + deltaSAM +
                     aveSOI + deltaSOI + 
                   # aveNstd + deltaNstd +
                     lineartrend + 
                     (moultstd|ID) + (1|year.random), tags1moult, REML = F)

notrend <-  lmer(julian ~ age.state + moultstd + weanstd +
                     aveSAM + deltaSAM +
                     aveSOI + deltaSOI + 
                    aveNstd + deltaNstd +
                   #  lineartrend + 
                     (moultstd|ID) + (1|year.random), tags1moult, REML = F)

noE <-  lmer(julian ~ age.state + moultstd + weanstd +
                   # aveSAM + deltaSAM +
                   #  aveSOI + deltaSOI + 
                   #  aveNstd + deltaNstd +
                   #   lineartrend + 
                   (moultstd|ID) + (1|year.random), tags1moult, REML = F)

AIC(full.slope.moult, noSAM, noSOI, noNstd, notrend, noE)

# Model parsimony improved when environmental factors and density were removed.

AIC(full.slope.moult, nowean, noagestate, nomoult, null, state, ageF)

models <- c("full.slope.moult", "null", "nomoult", "noagestate", "nowean", 
            "state", "ageF","noSAM", "noSOI", 
            "noNstd", "notrend", 'noE')

aics <- AIC(full.slope.moult, null, nomoult, noagestate, nowean, 
            state, ageF,noSAM, noSOI, 
            noNstd, notrend, noE)

delta.aics <- aics$AIC - min(aics$AIC) # To work out the change in delta AIC #see definitions in book (from Anderson (2008))
exp.delta <- exp(-0.5*delta.aics)
wi <- exp.delta/sum(exp.delta)      # these are the Akaike weights for each model #The probability that model is the actual (fitted) k-l best model in the set (Anderson 2008) (See Burnham et al (2011) paper) 
delta.null <- aics$AIC - aics$AIC[1] # AIC difference to null model

(modtable <- data.frame(models, numpar=aics$df,  aics, delta.aics, wi, delta.null))

# The AIC score of models that did not contain moult date, wean date or age increased 
# by at least 12 units relative to the full model, indicating strong model support for these variables.
# Dropping state from the age.state variable increased by 2 AIC units.
# In contrast, model parsimony tended to improve when environmental factors 
# and density were removed from the model.

# Pull out the deviance 
dev.full.slope.moult =  glance(full.slope.moult) 
dev.null =  glance(null)
dev.nomoult =  glance(nomoult) 
dev.age.state =  glance(noagestate)
dev.nowean =  glance(nowean) 
dev.state =  glance(state) 
dev.ageF =  glance(ageF)
dev.noSAM =  glance(noSAM)
dev.noSOI =  glance(noSOI) 
dev.noNstd =  glance(noNstd)
dev.notrend =  glance(notrend)
dev.noE =  glance(noE)

# --
dev.full.slope.moult = dev.full.slope.moult$deviance
dev.null =  dev.null$deviance
dev.nomoult =  dev.nomoult$deviance
dev.age.state =  dev.age.state$deviance
dev.nowean  = dev.nowean$deviance
dev.state =  dev.state$deviance
dev.ageF =  dev.ageF$deviance
dev.noSAM = dev.noSAM$deviance
dev.noSOI =  dev.noSOI$deviance
dev.noNstd =  dev.noNstd$deviance
dev.notrend =  dev.notrend$deviance
dev.noE =  dev.noE$deviance

dev = c(dev.full.slope.moult,   
dev.null, dev.nomoult,dev.age.state,dev.nowean,dev.state, dev.ageF, dev.noSAM,
dev.noSOI,dev.noNstd,dev.notrend, dev.noE ) 

(modtable <- data.frame(models, numpar=aics$df,  delta.aics, wi, delta.null, dev))
modtable = modtable[with(modtable, order(delta.aics)), ]

write.csv(modtable, "./output/model.table_with_centering_May_2023.csv")


# Refit with REML 
null <- lmer(julian ~ 1 + (moultstd|ID)+ (1|year.random), tags1moult, REML = T)

best <- lmer(julian ~ age.state + moultstd + weanstd + (moultstd|ID) + (1|year.random), tags1moult, REML = T)
best

# Confidence intervals:           
confint.merMod(best, oldNames = FALSE)

# --------------------------------------------------------------------------
# fit a linear slope on age to show the slope / trend of decrease with age
# --------------------------------------------------------------------------

tags1moult$ageslope =
        pmax(tags1moult$ageclass, 8)
   
table(tags1moult$ageslope)
table(tags1moult$ageclass)

adults = subset(tags1moult, tags1moult$state == "B")

best3 <- lmer(julian ~ ageslope + moultstd + weanstd + (moultstd|ID) + 
                     (1|year.random), adults, REML = T)
summary(best3)
# Confidence intervals:           
# confint.merMod(best3, oldNames = FALSE)

#----------------------------------------------------------
# Variance partitioning
# https://cran.r-project.org/web/packages/sjPlot/vignettes/tab_mixed.html
# Summary of Mixed Models as HTML Table  Daniel Ludecke  2020-05-23
#----------------------------------------------------------
# tab_model(null)
# tab_model(best)

tot = 44.559 + 42.318 + 3.08 + 2.582
(42.318 / tot)*100
(44.559 / tot)*100
(3.080 / tot)*100

# now pull out the variances for the null model
var.null = tidy(null, effects="ran_pars")
(var.nullResidual = (var.null[5,4])^2)
(var.nullID = (var.null[1,4])^2) 
(var.nullYear = (var.null[4,4])^2) 
(var.nullmoultID = (var.null[3,4])^2) 
#Total phenotypic variance:
(var.null.total = var.nullID + var.nullmoultID + var.nullYear + var.nullResidual)

# proportion of variance accounted for by random effects:
var.nullYear/var.null.total  # proportion of var accounted for by year random effect
var.nullID/var.null.total   # proportion of var accounted for by ID random effect

# And calculate the PCV as follows:
# PCV between null model and age.state
# ((var.nullResidual - var.moult.baselineResidual)/var.nullResidual*100)
# ((var.nullYear - var.moult.baselineYear)/var.nullYear*100)
# ((var.nullID -var.moult.baselineID)/var.nullID*100)


# ====================================================================================
# Figure: Age class plus state
# ====================================================================================
# https://ggplot2-book.org/scale-colour.html#colour-blindness
colorBlindness::displayAllColors(viridis::viridis(6))
library(scales)
show_col(viridis_pal()(6))
viridis_purple = '#440154FF'
viridis_green =  '#7AD151FF'

pr <- ggeffect(best, terms = c("age.state"), type = "fe")

#----------------------------------------------
# the ggpredict data is in a dataframe pr and can be used to plot with ggplot
head(pr)
str(pr)

# first create columns with the same names as in the Results data frame
pr$ageclass.pred = as.factor(pr$x)
pr$julian = pr$predicted

# now plot the raw data. 
# Note that the y-axis is very narrow around 15 Oct 
# (so that the model predictions are shown in a sensible way)
# and thus the plot does not include all raw data.

z = as.data.frame(pr)
pr.plot1 = filter(pr,
                    ageclass.pred != "3" &
                    ageclass.pred != "4FTB" &
                    ageclass.pred != "5FTB" &
                    ageclass.pred != "6FTB" &
                    ageclass.pred != "7FTB")

pr.plot1

pr.plot2 = filter(z,
                    ageclass.pred == "3" |
                    ageclass.pred == "4FTB" |
                    ageclass.pred == "5FTB" |
                    ageclass.pred == "6FTB" |
                    ageclass.pred == "7FTB")

pr.plot2

pr.plot2$ageclass.pred = c(1, 1.9, 2.9, 3.9, 4.9)


g = ggplot(data = tags1moult, aes(x = factor(ageclass), y= julian))+
  stat_sum(aes(colour = state, size=..n..), alpha = 0.2, shape=16)+
  scale_colour_manual(values = c(viridis_green, viridis_purple))+
  scale_size_area(breaks=c(1,15,30),"N")+
  #  scale_y_continuous(breaks = scales::pretty_breaks(n = 8),limits = c(275, 290))+
  #  scale_x_discrete(breaks = unique(tags1moult$ageclass)) 
  scale_y_continuous(breaks = c(seq(265,285, by = 2)), 
                labels =c("22 Sep", "24 Sep","26 Sep", "28 Sep","30 Sep",
                          "2 Oct","4 Oct", "6 Oct","8 Oct","10 Oct","12 Oct"),
                     limits = c(265, 285))+
  scale_x_discrete(breaks = c(seq(3,21, by = 2))) 

# now add model estimates and 95% confidence intervals
g = g + geom_pointrange(data = pr.plot1, size = 0.6, shape = 21 , stroke = 0.5,  
                        fill = viridis_green, colour = "black", 
                        aes(x = factor(ageclass.pred), y = julian,
                                                   ymin = conf.low, ymax = conf.high)) +
  gg_theme()+
  theme(legend.position = "top")+
  xlab("Age") +
  ylab("Date")


g = g + geom_pointrange(data = pr.plot2, size = 0.6, stroke = 0.5,shape = 22, 
                        fill = viridis_purple, colour = "black", 
                        aes(x = ageclass.pred-0.1, y = julian,
                                                         ymin = conf.low, ymax = conf.high)) 
  
g

fig_age.state = g

# ## Save Plot 
# pdf("./plots/2022_Age class and state model predictions.pdf",
#     useDingbats = FALSE, width = 3.9, height = 3.9)
# print(fig_age.state)
# dev.off()

#-----------------------------------------------------------------------

# Add the linear regression from best4 model to the plot

# ====================================================================================
# Figure: Age class plus state
# ====================================================================================
A <- ggeffect(best3, terms = c("ageslope"), type = "fe")
head(A)

# first create columns with the same names as in the Results data frame
A$ageclass.pred = as.factor(A$x)
A$julian = A$predicted

A = as.data.frame(A)
A

g2 = g + 
  geom_line(data = A, size = 0.4, colour = viridis_green, 
                   aes(x = ageclass.pred, y = julian, group = 1))
  
g2 = g2 + 
  geom_line(data = A, size = 0.4, colour = viridis_green, linetype = "twodash", 
                   aes(x = ageclass.pred, y = conf.low, group = 1))
g2 = g2 + 
  geom_line(data = A, size = 0.4, colour = viridis_green, linetype = "twodash", 
                   aes(x = ageclass.pred, y = conf.high, group = 1))

g2 

## Save Plot 
pdf("./plots/Figure4.pdf",
     useDingbats = FALSE, width = 4, height = 4)
print(g2)
dev.off()

# ====================================================================================
# Figure: Wean trend with geom_hex
# ====================================================================================

# Make values for the x-axis on the z-standardized scale that corresponds to the real scale
#plot(tags1moult$wean.date, tags1moult$weanstd)
xaxs = lm(weanstd ~ wean.date, data = tags1moult)
summary(tags1moult$wean.date)
x.scale = predict(xaxs, data.frame(wean.date = seq(280, 340, by = 10)))
x.scale


# Set up colour palettes:

viridis <- viridis::viridis(11)
viridis_hcl <- colorspace::sequential_hcl(11,
                                          h = c(300, 75), c = c(35, 95), l = c(15, 90), power = c(0.8, 1.2))

plasma <- viridis::plasma(11)
plasma_hcl <- colorspace::sequential_hcl(11,
                                         h = c(-100, 100), c = c(60, 100), l = c(15, 95), power = c(2, 0.9))

pal <- colorspace::sequential_hcl(5, "Heat", rev = T)

# ----------------------------------------
# Use best model output:
# ----------------------------------------
pr <- ggeffect(best, terms = c("weanstd"), type = "fe")

# first create columns with the same names as in the Results data frame
pr$weanstd = pr$x
pr$julian = pr$predicted

# now plot the raw data. 
# Note that the y-axis is very narrow around 15 Oct 
# (so that the model predictions are shown in a sensible way)
# and thus the plot does not include all raw data.

g = ggplot(data = tags1moult, aes(x = weanstd, y= julian))+
  geom_hex(bins = 40) +
 #  scale_fill_gradientn(colours = terrain.colors(7, rev = T))+             # https://colorspace.r-forge.r-project.org/articles/hcl_palettes.html
  #scale_fill_gradientn(colours = colorspace::heat_hcl(7, rev = T))+
  # scale_fill_gradientn(colours = colorspace::diverge_hcl(7))+
  # scale_fill_gradientn(colours = plasma_hcl)+
  # scale_fill_gradientn(colours = viridis_hcl)+
  # scale_fill_gradientn(colours = pal)+
  #scale_fill_viridis(option="cividis", direction = -1, begin = 0, end = 0.98, alpha = 1)+
  scale_fill_viridis_c(direction = -1,begin = 0, end = 1, alpha = 1, option = 'rocket', name = "(N)")+
  scale_size_continuous(breaks=c(1,10,50),range=c(1,6),"N")+
  scale_y_continuous(breaks = c(seq(260,310, by = 10)), 
                     labels =c("17 Sep","27 Sep","07 Oct","17 Oct","27 Oct","06 Nov"),
                     limits = c(260, 310))+
#  scale_x_continuous(breaks = c(seq(280, 340, by = 10)), # real scale
   scale_x_continuous(breaks = x.scale,  # z-stand. scale 
                     labels =c("07 Oct",	"17 Oct",	"27 Oct",	"06 Nov","16 Nov",	"26 Nov",	"06 Dec"),
                     # limits = c(280, 340)) # real scale
                     limits = c(-3.2, 4.2))  # z-stand. scale

# now add model estimates and 95% confidence intervals
#g = g + geom_ribbon(data = pr, aes(ymin=conf.low,ymax=conf.high),
#                   alpha = 0.6, colour = NA, fill = "green")+

g = g + geom_line(data = pr,size = 0.8, linetype = "twodash", aes(x = weanstd, y = conf.low)) + 
  geom_line(data = pr,size = 0.8, linetype = "twodash", aes(x = weanstd, y = conf.high)) + 
  gg_theme()+
  theme(axis.text.x = element_text(angle = -90, vjust = 0.5, hjust = 0))+
  #theme(legend.position = "top")+
  theme(legend.title=element_blank(),
        legend.position = c(-0.02, 1.02),
        legend.direction="horizontal",
        legend.justification = c("left", "top"))+
  xlab("Weaning date (as pup)") +
  ylab("Breeding arrival date") 

g = g + geom_line(data = pr, size = 1.2, aes(x = weanstd, y = julian))
g

fig_wean = g

# Save png to disk
# ggsave("./plots/2022_geomhex_wean model predictions.png", device = 'png', width = 4, height = 4, dpi = 600, units = 'in', limitsize = T)

# ## Save Plot 
# pdf("./plots/2022_geomhex_wean model predictions.pdf",
#     useDingbats = FALSE, width = 4, height = 4.5)
# print(fig_wean)
# dev.off()

# ====================================================================================
# Figure: Moult trend with geom_hex
# ====================================================================================

moult.xaxs = lm(moultstd ~ moultjulian1NOV, data = tags1moult)
summary(tags1moult$moultjulian1NOV)
moult.x.scale = predict(moult.xaxs, data.frame(moultjulian1NOV = 
                                                 c(44, 61, 75, 92, 106)))
moult.x.scale

# ----------------------------------------
# Use best model output:
# ----------------------------------------
pr <- ggeffect(best, terms = c("moultstd"), type = "fe")

# first create columns with the same names as in the Results data frame
pr$moultstd = pr$x
pr$julian = pr$predicted

# now plot the raw data. 
# Note that the y-axis is very narrow around 15 Oct 
# (so that the model predictions are shown in a sensible way)
# and thus the plot does not include all raw data.

g = ggplot(data = tags1moult, aes(x = moultstd, y= julian))+
  geom_hex(bins = 40) +
#  scale_fill_gradientn(colours = terrain.colors(7, rev = T))+             # https://colorspace.r-forge.r-project.org/articles/hcl_palettes.html
  #scale_fill_gradientn(colours = colorspace::heat_hcl(7, rev = T))+
  # scale_fill_gradientn(colours = colorspace::diverge_hcl(7))+
  # scale_fill_gradientn(colours = plasma_hcl)+
  # scale_fill_gradientn(colours = viridis_hcl)+
  # scale_fill_gradientn(colours = pal)+
 # scale_fill_viridis(option="cividis", direction = -1, begin = 0, end = 0.98, alpha = 1)+
  scale_fill_viridis_c(direction = -1,begin = 0, end = 1, alpha = 1, option = 'rocket')+
  
  scale_size_continuous(breaks=c(1,10,50),range=c(1,6),"N")+
  scale_y_continuous(breaks = c(seq(260,310, by = 10)), 
                     labels =c("17 Sep","27 Sep","07 Oct","17 Oct","27 Oct","06 Nov"),
                     limits = c(260, 310))+
#  scale_x_continuous(breaks = c(seq(280, 340, by = 10)), # real scale
   scale_x_continuous(breaks = moult.x.scale,  # z-stand. scale 
                     labels =c("15 Dec",
                               "01 Jan",	
                               "15 Jan",	
                               "01 Feb",
                               "15 Feb"),
                     # limits = c(280, 340)) # real scale
                     limits = c(-2.6, 2))  # z-stand. scale

g = g + geom_line(data = pr, size = 0.8, linetype = "twodash", aes(x = moultstd, y = conf.low)) + 
  geom_line(data = pr,size = 0.8, linetype = "twodash", aes(x = moultstd, y = conf.high)) + 
  
  gg_theme()+
  theme(axis.text.x = element_text(angle = -90, vjust = 0.5, hjust = 0))+
#  theme(legend.position = "top")+
  theme(legend.title=element_blank(),
        legend.position = c(-0.02, 1.02),
        legend.direction="horizontal",
        legend.justification = c("left", "top"))+
  xlab("Moult departure date") +
  ylab("Breeding arrival date")

g = g + geom_line(data = pr, size = 1.2, aes(x = moultstd, y = julian))
g

fig_moult = g

# Save png to disk
# ggsave("./plots/2022_geomhex_moult model predictions.png", device = 'png', width = 4, height = 4, dpi = 600, units = 'in', limitsize = T)

## Save Plot 
# pdf("./plots/2022_geomhex_moult model predictions23.pdf",
#     useDingbats = FALSE, width = 4, height = 4.5)
# print(fig_moult)
# dev.off()


# Plot all figures together

cowplot::plot_grid(fig_wean, 
                  fig_moult,   
                  labels = "AUTO",
                  nrow = 1,
                  align = "v")

## Save Plot 
pdf("./plots/Figure2.pdf",
     useDingbats = FALSE, width = 8, height = 4.5)
 cowplot::plot_grid(fig_wean, 
                    fig_moult,   
                    labels = "AUTO",
                    nrow = 1,
                    align = "v")
dev.off()

#=============================================================================
#----------------------------------------------------------------------------------------------------
# Correlating FTB arrival date with B arrival dates
# Not part of the model selection, as you are otherise using some response data as a predictor too
#=============================================================================
#----------------------------------------------------------------------------------------------------
# First select only experienced breeders:
EBtags = subset(tags1moult, tags1moult$state == "B") 

# fit lmer: test whether experienced breeder (this data set) arrival is predicted by 1st breeding date 
FTB.trend <- lmer(julian ~ first.breeding.date + (1|ID)+ (1|year.random), EBtags, REML = T)
 
# Confidence intervals:           
# tpr <- profile(FTB.trend)
# (confint(tpr) -> CIpr)

# Perhaps fit the full model (other important covariates). It has very little impact on the slope
FTB.trend2 <- lmer(julian ~ first.breeding.date + 
                            weanstd + 
                            age.state + 
                            moultstd + 
                            (moultstd|ID) + (1|year.random), 
                       EBtags, REML = T)


# ====================================================================================
# Figure: First breeding trend with geom_hex
# ====================================================================================
pr <- ggeffect(FTB.trend2, terms = c("first.breeding.date"), type = "fe")

# first create columns with the same names as in the Results data frame
pr$first.breeding.date = pr$x
pr$julian = pr$predicted

# now plot the raw data. 
# Note that the y-axis is very narrow around 15 Oct 
# (so that the model predictions are shown in a sensible way)
# and thus the plot does not include all raw data.

g = ggplot(data = EBtags, aes(x = first.breeding.date, y= julian))+
  #geom_bin2d(bins = 50) +
  geom_hex(bins = 40) +
#  scale_fill_gradientn(colours = terrain.colors(7, rev = T))+  
#  scale_fill_gradientn(colours = terrain.colors(7))+             # https://colorspace.r-forge.r-project.org/articles/hcl_palettes.html
#  scale_fill_gradientn(colours = colorspace::heat_hcl(10))+
#  scale_fill_gradientn(colours = colorspace::diverge_hcl(7))+
  scale_fill_viridis_c(direction = -1,begin = 0, end = 1, alpha = 1, option = 'rocket', breaks = c(5, 10, 15, 20,25))+    # could insert into bracket: breaks = c(0, 5, 10, 15, 20)
  scale_size_continuous(breaks=c(1,10,50),range=c(1,6),"N")+
  scale_y_continuous(breaks = c(seq(250,310, by = 10)), 
                     labels =c("07 Sep", "17 Sep","27 Sep","07 Oct","17 Oct","27 Oct","06 Nov"),
                     limits = c(250, 310))+
  scale_x_continuous(breaks = c(seq(250, 310, by = 10)), 
                     labels =c("07 Sep", "17 Sep","27 Sep","07 Oct","17 Oct","27 Oct","06 Nov"),
                     limits = c(250, 310)) 

# now add model estimates and 95% confidence intervals
#g = g + geom_ribbon(data = pr, aes(ymin=conf.low,ymax=conf.high),
#                    alpha = 0.6, colour = NA, fill = "darkorchid4")+  
g = g + geom_line(data = pr,size = 0.8, linetype = "twodash", aes(x = first.breeding.date, y = conf.low)) + 
  geom_line(data = pr,size = 0.8, linetype = "twodash", aes(x = first.breeding.date, y = conf.high)) + 
  gg_theme()+
  theme(axis.text.x = element_text(angle = -90, vjust = 0.5, hjust = 0))+
  theme(legend.title=element_blank(),
        legend.position = c(-0.01, 1),
        legend.direction="horizontal",
        legend.justification = c("left", "top"))+
#  guides(fill = guide_colourbar(barwidth = 0.5, barheight = 10))+
  guides(color=guide_legend(override.aes=list(fill=NA)))+
    xlab("First breeding arrival date") +
  ylab("Subsequent arrival date(s)")

g = g + geom_line(data = pr, size = 1.2, aes(x = first.breeding.date, y = julian))
g

fig_first = g

# This comment is only applicable to the case where you do this analysis on all data (and not 
# only on Exp breeders, as here)
# The blue FTB lines are on the 1:1 line (because it plots the observations against itself)
# The series of 'higher' blue lines are because FTB are seen multiple times in the breeding season
# (They appear in 7 day bands).
# These bands will disappear if you only use the first breeding encounter.
# This comparision should perhaps happen without FTB (the result is the same)

# Save png to disk
# ggsave("./plots/2022_geom_hex_First breeding predictions.png", device = 'png', width = 4, height = 4, dpi = 600, units = 'in', limitsize = T)
# 
## Save Plot 
pdf("./plots/Figure5.pdf",
     useDingbats = FALSE, width = 4, height = 4)
print(fig_first)
dev.off()


# ====================================================================================
# Figure: random effects
# ====================================================================================

# null model
t <- tidy(null, effects="ran_vals")
t <- tidy(best, effects="ran_vals")

# Plot random effects by year
tt = subset(t, group == "year.random")
ggplot(tt,aes(level,estimate))+
  geom_pointrange(aes(ymin=estimate-1.96*std.error,
                      ymax=estimate+1.96*std.error))+
  gg_theme()


# Plot random effects by year
tt = subset(t, group == "year.random")

# plot in order of conditional mode:
tt = tt[order(tt$estimate),]
tt$plotorder = seq(1:length(tt$estimate))

# plot in order of years:
tt = tt[order(tt$level),]
tt$plotorder = seq(1:length(tt$level))

ry =  ggplot(tt, aes(plotorder, estimate))+
    scale_x_continuous(breaks = c(seq(1,30, by = 1)),
    labels = (tt$level)) +
  #scale_y_continuous(limits = c(-25, 27.5), breaks = c(seq(-20,20, by = 10)))+  # same scale as ID
  scale_y_continuous(limits = c(-5, 7.5), breaks = c(seq(-10,10, by = 5)))+
      geom_pointrange(aes(ymin=estimate-1.96*std.error,
    ymax=estimate+1.96*std.error),
    size = 0.8, # line width 
    fatten = 1.3, # point size
    color= "black", fill = "darkorchid4", shape=21)+
  geom_point(aes(plotorder, estimate), colour = "darkorchid4", size = 0.8)+
  
  gg_theme()+
  ylab("Conditional mode") +
  xlab("Year")+
  coord_flip()+ geom_hline(yintercept = 0, linetype = 'dotted')

ry

## Save Plot 
# pdf("./plots/2022 year random effect plot.pdf",
#     useDingbats = FALSE, width = 4, height = 4.5)
# print(ry)
# dev.off()

# ------------------------------
# plot moult random slopes
# ------------------------------

# If you add other covariates, then you add population slopes.
# If you add the year random effect, then your blue lines get added annual variation.
# Thus, fit only moult date, as a fixed and random effect on individual to illustrate the shape

# to make the plots more readable, plot cohorts per BIRTH cohort!

# Make x axis variable (because it was standardized)
xaxs = lm(moultstd ~ moultjulian1NOV, data = tags1moult)
x.scale = predict(xaxs, data.frame(moultjulian1NOV = c(1,15,31,45,62,76,93,107,121,135)))
x.scale
# 1	01-Nov
# 15	15-Nov
# 31	01-Dec
# 45	15-Dec
# 62	01-Jan
# 76	15-Jan
# 93	01-Feb
# 107	15-Feb
# 121	01-Mar
# 135	15-Mar

# Code from Ben Bolker (in Leandri MSc folder - search for Bolker)

p1 <- lmer(julian ~ moultstd + (moultstd|ID) , 
         tags1moult, 
         REML = T, na.action = "na.fail")

tags1moult$pop.pred <- predict(p1, re.form=NA)  ## population level
tags1moult$ind.pred <- predict(p1) ## individual level

gg <- ggplot(tags1moult, aes(moultstd, julian)) +
         geom_line(aes(y = pop.pred,  group=ID), colour="black", size = 1.3, alpha = 1) + 
         geom_line(aes(y = ind.pred,  group=ID), colour="darkorchid4", size = 0.05, alpha = 0.3) +
  
  gg_theme() + 
  xlab("Moult departure date") +
  ylab("Breeding arrival date")+
  theme(axis.text.x = element_text(angle = -90, vjust = 0.5, hjust = 0)) +
  scale_x_continuous(breaks = x.scale,      # z-standardized scale 
                     labels =c("01 Nov","15 Nov","01 Dec",	"15 Dec",
                               "01 Jan","15 Jan","01 Feb",	"15 Feb",	
                               "01 Mar", "15 Mar"),
                     limits = c(-4, 3.8)) +
scale_y_continuous(breaks = c(seq(260,300, by = 10)), 
                     labels =c("17 Sep","27 Sep","07 Oct","17 Oct","27 Oct"),
                     limits = c(260, 300)) + 
        geom_label(aes(x = -3.02, y = 298, label = "1983 - 2016"), fontface = "plain", size = 4.2,
                  label.size = 0)
gg

# ## Save Plot 
# pdf("./plots/2022 random effect slopes_1986_2016_purple.pdf",
#     useDingbats = FALSE, width = 4, height = 4.5)
# print (gg)
# dev.off()


## Save Plot 
pdf("./plots/Figure3.pdf",
     useDingbats = FALSE, width = 9, height = 5.5)
cowplot::plot_grid( gg,ry, labels = "AUTO",
                    nrow = 1,
                    align = "v")
dev.off()

# Now subset predicted values per cohort
g83 = tags1moult %>% dplyr::filter(between(birthyear, 1983, 1989))
g90 = tags1moult %>% dplyr::filter(between(birthyear, 1990, 1994))
g95 = tags1moult %>% dplyr::filter(between(birthyear, 1995, 1999))
g00 = tags1moult %>% dplyr::filter(between(birthyear, 2000, 2004))
g05 = tags1moult %>% dplyr::filter(between(birthyear, 2005, 2009))
g10 = tags1moult %>% dplyr::filter(between(birthyear, 2010, 2016))
  

gg83 <- ggplot(g83, aes(moultstd, julian)) +
         geom_line(aes(y = pop.pred,  group=ID), colour="black", size = 1.3, alpha = 1) + 
         geom_line(aes(y = ind.pred,  group=ID), colour="darkorchid4", size = 0.05, alpha = 0.3) +
  gg_theme() + 
  xlab("Moult departure date") +
  ylab("Breeding arrival date")+
  theme(axis.text.x = element_text(angle = -90, vjust = 0.5, hjust = 0)) +
  scale_x_continuous(breaks = x.scale,      # z-standardized scale 
                     labels =c("01 Nov","15 Nov","01 Dec",	"15 Dec",
                               "01 Jan","15 Jan","01 Feb",	"15 Feb",	
                               "01 Mar", "15 Mar"),
                     limits = c(-4, 3.8)) +
scale_y_continuous(breaks = c(seq(260,300, by = 10)), 
                     labels =c("17 Sep","27 Sep","07 Oct","17 Oct","27 Oct"),
                     limits = c(260, 300)) + 
        geom_label(aes(x = -3.02, y = 298, label = "1983 - 1989"), fontface = "plain", size = 4.2,
                  label.size = 0)
#gg83


gg90 <- ggplot(g90, aes(moultstd, julian)) +
         geom_line(aes(y = pop.pred,  group=ID), colour="black", size = 1.3, alpha = 1) + 
         geom_line(aes(y = ind.pred,  group=ID), colour="darkorchid4", size = 0.05, alpha = 0.3) +
  gg_theme() + 
  xlab("Moult departure date") +
  ylab("Breeding arrival date")+
  theme(axis.text.x = element_text(angle = -90, vjust = 0.5, hjust = 0)) +
  scale_x_continuous(breaks = x.scale,      # z-standardized scale 
                     labels =c("01 Nov","15 Nov","01 Dec",	"15 Dec",
                               "01 Jan","15 Jan","01 Feb",	"15 Feb",	
                               "01 Mar", "15 Mar"),
                     limits = c(-4, 3.8)) +
scale_y_continuous(breaks = c(seq(260,300, by = 10)), 
                     labels =c("17 Sep","27 Sep","07 Oct","17 Oct","27 Oct"),
                     limits = c(260, 300)) + 
        geom_label(aes(x = -3.02, y = 298, label = "1990 - 1994"), fontface = "plain", size = 4.2,
                  label.size = 0)
#gg90


gg95 <- ggplot(g95, aes(moultstd, julian)) +
         geom_line(aes(y = pop.pred,  group=ID), colour="black", size = 1.3, alpha = 1) + 
         geom_line(aes(y = ind.pred,  group=ID), colour="darkorchid4", size = 0.05, alpha = 0.3) +
  gg_theme() + 
  xlab("Moult departure date") +
  ylab("Breeding arrival date")+
  theme(axis.text.x = element_text(angle = -90, vjust = 0.5, hjust = 0)) +
  scale_x_continuous(breaks = x.scale,      # z-standardized scale 
                     labels =c("01 Nov","15 Nov","01 Dec",	"15 Dec",
                               "01 Jan","15 Jan","01 Feb",	"15 Feb",	
                               "01 Mar", "15 Mar"),
                     limits = c(-4, 3.8)) +
scale_y_continuous(breaks = c(seq(260,300, by = 10)), 
                     labels =c("17 Sep","27 Sep","07 Oct","17 Oct","27 Oct"),
                     limits = c(260, 300)) + 
        geom_label(aes(x = -3.02, y = 298, label = "1995 - 1999"), fontface = "plain", size = 4.2,
                  label.size = 0)
#gg95



gg00 <- ggplot(g00, aes(moultstd, julian)) +
         geom_line(aes(y = pop.pred,  group=ID), colour="black", size = 1.3, alpha = 1) + 
         geom_line(aes(y = ind.pred,  group=ID), colour="darkorchid4", size = 0.05, alpha = 0.3) +
  gg_theme() + 
  xlab("Moult departure date") +
  ylab("Breeding arrival date")+
  theme(axis.text.x = element_text(angle = -90, vjust = 0.5, hjust = 0)) +
  scale_x_continuous(breaks = x.scale,      # z-standardized scale 
                     labels =c("01 Nov","15 Nov","01 Dec",	"15 Dec",
                               "01 Jan","15 Jan","01 Feb",	"15 Feb",	
                               "01 Mar", "15 Mar"),
                     limits = c(-4, 3.8)) +
scale_y_continuous(breaks = c(seq(260,300, by = 10)), 
                     labels =c("17 Sep","27 Sep","07 Oct","17 Oct","27 Oct"),
                     limits = c(260, 300)) + 
        geom_label(aes(x = -3.02, y = 298, label = "2000 - 2004"), fontface = "plain", size = 4.2,
                  label.size = 0)
#gg00


gg05 <- ggplot(g05, aes(moultstd, julian)) +
         geom_line(aes(y = pop.pred,  group=ID), colour="black", size = 1.3, alpha = 1) + 
         geom_line(aes(y = ind.pred,  group=ID), colour="darkorchid4", size = 0.05, alpha = 0.3) +
  gg_theme() + 
  xlab("Moult departure date") +
  ylab("Breeding arrival date")+
  theme(axis.text.x = element_text(angle = -90, vjust = 0.5, hjust = 0)) +
  scale_x_continuous(breaks = x.scale,      # z-standardized scale 
                     labels =c("01 Nov","15 Nov","01 Dec",	"15 Dec",
                               "01 Jan","15 Jan","01 Feb",	"15 Feb",	
                               "01 Mar", "15 Mar"),
                     limits = c(-4, 3.8)) +
scale_y_continuous(breaks = c(seq(260,300, by = 10)), 
                     labels =c("17 Sep","27 Sep","07 Oct","17 Oct","27 Oct"),
                     limits = c(260, 300)) + 
        geom_label(aes(x = -3.02, y = 298, label = "2005 - 2009"), fontface = "plain", size = 4.2,
                  label.size = 0)
#gg05



gg10 <- ggplot(g10, aes(moultstd, julian)) +
         geom_line(aes(y = pop.pred,  group=ID), colour="black", size = 1.3, alpha = 1) + 
         geom_line(aes(y = ind.pred,  group=ID), colour="darkorchid4", size = 0.05, alpha = 0.3) +
  gg_theme() + 
  xlab("Moult departure date") +
  ylab("Breeding arrival date")+
  theme(axis.text.x = element_text(angle = -90, vjust = 0.5, hjust = 0)) +
  scale_x_continuous(breaks = x.scale,      # z-standardized scale 
                     labels =c("01 Nov","15 Nov","01 Dec",	"15 Dec",
                               "01 Jan","15 Jan","01 Feb",	"15 Feb",	
                               "01 Mar", "15 Mar"),
                     limits = c(-4, 3.8)) +
scale_y_continuous(breaks = c(seq(260,300, by = 10)), 
                     labels =c("17 Sep","27 Sep","07 Oct","17 Oct","27 Oct"),
                     limits = c(260, 300)) + 
        geom_label(aes(x = -3.02, y = 298, label = "2010 - 2016"), fontface = "plain", size = 4.2,
                  label.size = 0)
#gg10



cowplot::plot_grid(gg83,
                   gg90,
                   gg95,
                   gg00,
                   gg05,
                   gg10,
                   labels = "AUTO",
                   nrow = 3,
                   align = "v")


## Save Plot
pdf("./plots/FigureS9.pdf",
    useDingbats = FALSE, width = 8, height = 10)

cowplot::plot_grid(gg83,
                   gg90,
                   gg95,
                   gg00,
                   gg05,
                   gg10,
                   labels = "AUTO",
                   nrow = 3,
                   align = "v")
dev.off()


#-----------------------------------------------------------------------------------------
# Calculate Repeatability (R)
#-----------------------------------------------------------------------------------------

# ----------------------------------------------
# REPEATABILITY: nakagawa and Schielzeth 2010
# https://cran.r-project.org/web/packages/rptR/vignettes/rptR.html
# ----------------------------------------------

library(rptR)
# If you want to run the random slope model, then you have to have 
# moultstd as a fixed effect also, otherwise the rpt function gives an error
# about subscript out of bounds
null.rpt = rpt(julian ~ moultstd + (moultstd|ID)+ (1|year.random), grname = c("ID", "year.random"), 
           data = tags1moult, datatype = "Gaussian",
           nboot = 100, npermut = 0, adjusted = FALSE)

null.rpt

best.rpt = rpt(julian ~ age.state + moultstd + weanstd +
              (moultstd|ID)+ (1|year.random), grname = c("ID", "year.random"), 
           data = tags1moult, datatype = "Gaussian",
           nboot = 100, npermut = 0, adjusted = FALSE)
best.rpt 

# The results show substantial repeatable variation among ID
# with low variation among years.

print(null.rpt)
plot(null.rpt, grname = "ID", type = "boot", cex.main = 0.8)
plot(null.rpt, grname = "year.random", type = "boot", cex.main = 0.8)

summary(null.rpt$mod)

#Adjusted repeatabilities for Gaussian data
adj.rpt = rpt(julian ~ age.state + moultstd + weanstd + (moultstd|ID) + (1|year.random),
              grname = c("ID", "year.random"), 
              data = tags1moult, datatype = "Gaussian",
              nboot = 1000, npermut = 0, parallel=TRUE)

print(adj.rpt)
plot(adj.rpt, grname = "year.random", scale= "link", cex.main=0.8, main= "Year variance", las=1)
plot(adj.rpt, grname = "ID", scale= "link", cex.main=0.8, main= "ID variance", las=1)


#---------------------------------------------------------------------------------------------
# Heritability in arrival date
#---------------------------------------------------------------------------------------------

# now import mother-pup association data from supersmall records
ss <- read.csv("./data/mom_pup.csv")  
head(ss)

# Create ID columns:
# format the tag number as 000 
ss$mom_tag = paste(formatC(ss$mom_tag, width=3, flag="0"))
ss$pup_tag = paste(formatC(ss$pup_tag, width=3, flag="0"))

# concatenate with Cohort
ss$momID = paste(ss$mom_cohort, ss$mom_tag, sep = "")
ss$pupID = paste(ss$pup_cohort, ss$pup_tag, sep = "")

## clean up slightly 
ss = ss %>% dplyr::select(c(momID, pupID))
head(ss)  # these are the mom-pup pairs from supersmall data

dim(ss)
n_distinct(ss$momID) # how many mothers
n_distinct(ss$pupID) # all pup IDs are unique (good)

# -----------------------------------------------------------------
# Merge mom-pup SS data with tag data (to get dates of weaning and first breeding)
# This is a bit of a hack...
# -----------------------------------------------------------------

# assign new names to merge by ID:
names(ss) = c("ID", "pup")

# get wean date for moms
m <- merge(y = ss, x = tags, by = "ID", all = F)
dim(m)

## clean up slightly 
mm = m %>% dplyr::select(c(ID,  wean.date, first.breeding.date, pup, year, birthyear, AFB))
mm = mm[!duplicated(mm[c(4)]),]
dim(mm)

head(mm)
names(mm) = c("momID", "mom.wean", "mom.first", "ID", "reproyear", "mom.birthyear", "AFB.m")
head(mm)

# Now get wean date for pups
p <- merge(y = mm, x = tags, by = "ID", all = F)
dim(p)

head(p)
## clean up slightly 
pp =  p %>% dplyr::select(c(ID,  wean.date,  first.breeding.date , momID , mom.wean, mom.first, mom.birthyear, birthyear, AFB.m, AFB ))

head(pp)
names(pp) = c("pupID", "pupwean", "pup.first", "momID", "momwean", "mom.first","mom.birthyear","pup.birthyear", "AFB.m","AFB.pup")
head(pp)
dim(pp)
pp = pp[!duplicated(pp[c(1)]),]
dim(pp)

# ===========================================================================================
# Falconer Page 168
# Calculating heritability regression on the mean parent and the mean offspring values
# rather than each of the individual values
# ===========================================================================================

# Are there females with more than 1 pup in the heritability data?
n_distinct(pp$momID)  

# yes, the pp df is 216 rows long (there is 216 mom-pup records)
# but the unique df is only 187 records, meaning that some females have more than
# one pup in this data. Calculate the mean parent and mean offspring arrival dates 
# as done by Falconer

216 - 187

N = pp %>% 
  group_by(momID) %>% 
  tally()

table(N$n)

Falconer = pp %>% 
  group_by(momID) %>% 
  mutate(pups.ave = mean(pup.first)) %>%
  mutate(moms.ave = mean(mom.first)) %>%
  distinct(momID, .keep_all = TRUE) %>%
  ungroup()

Falconer   #same nr of rows as unique moms

mean(Falconer$pup.first)
mean(Falconer$mom.first)

#linear model
lm.fit.Falconer = lm(pups.ave ~ moms.ave, data =  Falconer)
summary(lm.fit.Falconer)
anova(lm.fit.Falconer)

plot(Falconer$pups.ave, Falconer$moms.ave, xlim = c(250,330), ylim = c(250,330))
abline(coef = c(0,1))
abline(lm.fit.Falconer, lwd=2, col="red")

#add sample size to Flaconer data (records where there are more than 1 mom-pup pair with
# the same dates)

samplesize = Falconer %>% 
  group_by(pups.ave, moms.ave) %>%
  mutate(n= n()) %>% 
  distinct(pupID, .keep_all=TRUE)

samplesize   # same nr of rows as Falconer
table(samplesize$n)  # 20 records are the same date (?)

#h.plot = ggplot(data = Falconer, aes(x = pups.ave, y = moms.ave)) + 
h.plot = ggplot(data = samplesize, aes(x = pups.ave, y = moms.ave, color = as.factor(n))) + 
  geom_point(#position=position_jitter(h=0.1, w=0.1),
    #color = "tan2", alpha = 0.5, shape = 15, size = 1.5) +
    alpha = 0.9, shape = 15, size = 1.5) +
  scale_colour_manual(values = c("#4daf4a", "#984ea3"))+
  scale_y_continuous(breaks = c(seq(250,310, by = 10)), 
                     labels =c("07 Sep","17 Sep","27 Sep","07 Oct","17 Oct","27 Oct","06 Nov"),
                     limits = c(245, 320))+
  scale_x_continuous(breaks = c(seq(250,318, by = 10)), 
                     labels =c("07 Sep","17 Sep","27 Sep","07 Oct","17 Oct","27 Oct","06 Nov"),
                     limits = c(250, 308)) + 
  gg_theme() + 
  theme(axis.text.x = element_text(angle = -90, vjust = 0.5, hjust = 0))+
  theme(legend.position = "top")+
  xlab("Daughter arrival date") +
  ylab("Mother arrival date") +
  geom_ribbon(stat='smooth', method = "lm",  fill = NA,  outline.type = 'both',
              aes(), color = 'black', size = 0.8, linetype = "twodash" )+
  geom_line(stat='smooth', method = "lm", alpha = 1 ,color = "black", size = 1.2)

h.plot

n_distinct(Falconer$pups.ave)  # there is overlap of points
n_distinct(Falconer$moms.ave)  # there is overlap of points


# t = subset(samplesize, samplesize$n == 2)
# 20 records are duplicate, thus it is 10 pairs. Thus, on the graph, there is only 10
# purple blocks. Correct.

sort(samplesize$moms.ave)

# This plot excludes 3 records with mom arrival dates after 1 November
# and it excludes 1 record with a mom arrival date of 6 Sept. 

h.plot = ggplot(data = samplesize, aes(x = moms.ave,  y = pups.ave, color = as.factor(n))) + 
  geom_point(#position=position_jitter(h=0.1, w=0.1),
    #color = "tan2", alpha = 0.5, shape = 15, size = 1.5) +
    alpha = 0.95, shape = 15, size = 1.5) +
  scale_colour_manual(values = c('orangered', '#440154FF'))+
  scale_y_continuous(breaks = c(seq(250,310, by = 10)), 
                     labels =c("07 Sep","17 Sep","27 Sep","07 Oct","17 Oct","27 Oct","06 Nov"),
                     limits = c(252, 305))+
  scale_x_continuous(breaks = c(seq(250,318, by = 10)), 
                     labels =c("07 Sep","17 Sep","27 Sep","07 Oct","17 Oct","27 Oct","06 Nov"),
                     limits = c(255, 305)) + 
  gg_theme() + 
  theme(axis.text.x = element_text(angle = -90, vjust = 0.5, hjust = 0))+
  # theme(legend.position=c(0.1, 0.87), legend.title = element_blank())+
  theme(legend.position = "none") +
  # legend position inside the plot is given as between (0 and 1), not on the real axis scale
  xlab("Mother arrival date") +
  ylab("Daughter arrival date") +
  geom_ribbon(stat='smooth', method = "lm",  fill = NA,  outline.type = 'both',
              aes(), color = 'black', size = 0.8, linetype = "twodash" )+
  geom_line(stat='smooth', method = "lm", alpha = 1 ,color = "black", size = 1.2)

h.plot

## Save Plot 
pdf("./plots/Figure6.pdf",
    useDingbats = FALSE, width = 5.5, height = 5.5)
print(h.plot)
dev.off()

plot(samplesize$moms.ave, samplesize$pups.ave)
abline(lm.fit.Falconer, lwd=2, col="red")

#==============================================================================================

#----------------------------------------------
# Sanity check
#----------------------------------------------
# write data to go double check the correctness
# write.csv(pp, "mom_pup_weaning dates from R.csv")
# subset(tags1, ID == "GR(1)444")
# wean dates for moms and pups in pp object matches those in tags1 object. 
# Sanity check passed
#----------------------------------------------

###
plot(pp$pup.first, pp$mom.first, xlim = c(250,330), ylim = c(250,330))
abline(coef = c(0,1))

# correlation?
r = cor(pp$pup.first, pp$mom.first)
r

#linear model R2?
lm.fit = lm(pp$pup.first ~  pp$mom.first)
summary(lm.fit)

R2 = summary(lm.fit)$r.squared     # names(summary(lm.fit)) gives you the variables you can extract

# now extract the model p value:
lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

lmp(lm.fit)


# Print results
cat("Correlation coefficient = ", r)
## Correlation coefficient =  0.8068949
cat("Square of corr. coef. = ", r*r)
## Square of corr. coef. =  0.6510794
cat("Coefficient of determination= ", R2)
## Coefficient of determination=  0.6510794

#====================================================================================================================
# I could not get a mixed model to work with pup.first.
# Problem with singularity that seems to be associated with this variable (because pupwean, for example, works as a lmer)
pp$pup.birthyearf = as.factor(pp$pup.birthyear)
pp$momIDf = as.factor(pp$momID)
mm = lmer(pup.first ~ mom.first  + (1|momIDf)+  (1|pup.birthyearf), pp)  # singular
mm = lmer(pupwean ~ mom.first  + (1|momIDf)+  (1|pup.birthyearf), pp)  # works
#====================================================================================================================


# --------------------------------------------------
# select only those mom-pup pairs with the same AFB
# --------------------------------------------------

ppp = subset(pp, pp$AFB.m == pp$AFB.pup)
dim(ppp)
as.data.frame(ppp)

# correlation?
cor(ppp$pup.first, ppp$mom.first)

#linear model R2?
summary(lm(ppp$pup.first ~  ppp$mom.first))

AFB3 = subset(ppp, ppp$AFB.m == 3)
AFB4 = subset(ppp, ppp$AFB.m == 4)

plot(ppp$pup.first, ppp$mom.first, xlim = c(250,330), ylim = c(250,330))
points(AFB3$pup.first, AFB3$mom.first, pch = 16)
points(AFB4$pup.first, AFB4$mom.first, pch = 16, col = "blue")
abline(coef = c(0,1))

ppp = ppp[order(ppp$mom.first),]
tail(ppp)

#---------------------------------------------------
# Can I do a bootstrap?
#---------------------------------------------------

# What is the correlation of a random sample (same sample size as above ppp)

tagsAFB3 = subset(tags, tags$AFB == 3)
tagsAFB4 = subset(tags, tags$AFB == 4)

rand.pup3 = sample_n(tagsAFB3, 54)
rand.pup4 = sample_n(tagsAFB4, 60)
rand.pup = rbind(rand.pup3, rand.pup4)

rand.mom3 = sample_n(tagsAFB3, 54)
rand.mom4 = sample_n(tagsAFB4, 60)
rand.mom = rbind(rand.mom3, rand.mom4)

# correlation?
cor(rand.pup$first.breeding.date, rand.mom$first.breeding.date)

plot(rand.pup$first.breeding.date, rand.mom$first.breeding.date, xlim = c(250,330), ylim = c(250,330))
abline(coef = c(0,1))

#--------------------------------------------------------------------------------------
# BOOTSTRAP A RANDOM SAMPLE 
# this code draws 1000 data sets, and then conducts a lm on each of the random data sets
# the random data sets assume every individual can have another individual as mother.
# to see if the code works, change:
# Analysis function    to 
# for(i in 1:1){  
# then run code on this single random data set.
# Compare output of summary(lm(rand.pup$first.breeding.date ~  rand.mom$first.breeding.date))   and
# mean(rand.out) # mean R2
# mean(rand.out2) # mean values
# That should give the same answers
#--------------------------------------------------------------------------------------

output.R2 = list()
output.coef = list()
output.pvalue = list()

repeats = 1000

# Analysis function
for(i in 1:repeats){                          
  
  tagsAFB3 = subset(tags, tags$AFB == 3)
  tagsAFB4 = subset(tags, tags$AFB == 4)
  
  rand.pup3 = sample_n(tagsAFB3, 54)
  rand.pup4 = sample_n(tagsAFB4, 60)
  rand.pup = rbind(rand.pup3, rand.pup4)
  rand.pup$pup.first = rand.pup$first.breeding.date
  
  rand.mom3 = sample_n(tagsAFB3, 54)
  rand.mom4 = sample_n(tagsAFB4, 60)
  rand.mom = rbind(rand.mom3, rand.mom4)
  rand.mom$mom.first = rand.mom$first.breeding.date
  
  # correlation?
  # temp = cor(rand.pup$first.breeding.date, rand.mom$first.breeding.date)
  
  #linear model R2?
  temp =  lm(rand.pup$first.breeding.date ~  rand.mom$first.breeding.date)
  
  #   plot(rand.pup$first.breeding.date, rand.mom$first.breeding.date,  xlim = c(250,330), ylim = c(250,330))
  #  abline(temp)
  
  #  output[[i]] = temp
  
  output.R2[[i]] =  summary(temp)$r.squared
  output.coef[[i]] =  summary(temp)$coefficients[2,1]
  output.pvalue[[i]] = lmp(temp)
  rm(temp)
}

#output  # is a list: class(output)
# to work with output data, unlist:

output.R2 = unlist(output.R2, use.names=FALSE)
output.coef = unlist(output.coef, use.names=FALSE)
output.pvalue = unlist(output.pvalue, use.names=FALSE)

mean(output.R2) # mean R2
mean(output.coef) # mean values
mean(output.pvalue) # mean values

# to get a confidence interval of 95%, we would select the value at the 2.5% percentile 
# as the lower bound and the 97.5% percentile as the upper bound.
# if we calculated 1,000 statistics from 1,000 bootstrap samples, then the 
# lower bound would be the 25th value and the upper bound would be the 975th value, 
# assuming the list of statistics was ordered.

output.coef <- sort(output.coef, decreasing = FALSE)

lcl <- output.coef[ceiling(repeats/40)]            # the [ceiling(repeats/40)] selects the 
ucl <- output.coef[repeats-ceiling(repeats/40)]    # 2.5 and 97.5 percentiles

lcl #lower confidence limits
ucl #upper confidence limits

hist(output.R2)
hist(output.coef)
hist(output.pvalue)

lower95 =quantile(output.coef, c(0.025))
upper95 = quantile(output.coef, c(0.975))

quantile(output.coef, c(0.025, 0.975))
quantile(output.pvalue, c(0.025, 0.975))

output.coef = as.data.frame(output.coef)
# Histogram with density instead of count on y-axis
bootplot = ggplot(output.coef, aes(x=output.coef)) + 
  geom_histogram(aes(y=..density..),      
                 bins=50,
                 colour="black", fill="grey")+
  # geom_density(alpha=.2, fill="blue") +
  gg_theme() +
  geom_vline(aes(xintercept=lm.fit$coefficients[2]), 
             color="blue", linetype="dashed", size=1)+
  xlab("Slope coefficient") +
  ylab("Density") +
  geom_vline(aes(xintercept=lower95), 
             color="grey", linetype="dotted", size=1)+
  
  geom_vline(aes(xintercept=upper95), 
             color="grey", linetype="dotted", size=1)

bootplot

## Save Plot 
pdf("./plots/FigureS10.pdf",
    useDingbats = FALSE, width = 6, height = 5)
ggdraw(bootplot)
dev.off()

# ====================================================================
# - HOW TO DO BOOTSTRAP FROM MOM-PUP PAIRS *SUBSAMPLE MOM-PUP DATA* 
# ====================================================================

# Bootstrap the observed data (subsample mom-pup pairs)

library(boot)

# function to obtain R-Squared from the data
rsq <- function(formula, data, indices) {
  d <- data[indices,] # allows boot to select sample
  fit <- lm(formula, data=d)
  return(summary(fit)$r.square)
}

# bootstrapping with 1000 replications

results <- boot(data= pp, statistic = rsq,
                R=1000, formula = pup.first ~ mom.first)

# view results
plot(results)

# get 95% confidence interval
boot.ci(results, type="bca")


# function to obtain model coefficient from the data
rcoef <- function(formula, data, indices) {
  d <- data[indices,] # allows boot to select sample
  fit <- lm(formula, data=d)
  return(summary(fit)$coefficients[2,1])
}

# bootstrapping with 1000 replications
results <- boot(data= pp, statistic = rcoef,
                R = 1000, formula = pup.first ~ mom.first)

# view results
results
plot(results)

# get 95% confidence interval
boot.ci(results, type="bca")

bootresults = results$t
bootresults = as.data.frame(bootresults)

gg = cbind(output.coef,bootresults)
names(gg) = c("booot.coef"," booot.results")

gg = gg%>%
  gather(booot,output.coef)

str(gg)

ggsim = subset(gg, gg$booot == "booot.results")

# Overlaid histograms are created by setting the argument position="identity".
# We have also set the alpha parameter as alpha=.5 for transparency.
ggplot(gg, aes(x= output.coef, fill= booot)) + 
  geom_histogram(aes(y=..density..),      
                 bins=50,
                 colour="black",  alpha=0.4, position= "identity")+
  # geom_density(alpha=.2, fill="blue") +
  geom_vline(aes(xintercept=mean(ggsim$output.coef)), color="blue", linetype="dashed", size=1)+
  gg_theme()


#--------------------------------------
# Model checking
#--------------------------------------

library(performance)
# checking model assumptions
#dev.new()

## Save Plot 
pdf("./plots/FigureS8.pdf",
   useDingbats = FALSE, width = 10, height = 6)

check_model(best, panel = T, dot_size = 1, alpha = 0.15, line_size = 0.6,
            check = c(#"qq",
                      "normality", 
                      "linearity", 
                      "homogeneity"
                     #"pp_check"
                     ))
            
dev.off()

# pdf("./plots/RE model diagnostics.pdf",
#     useDingbats = FALSE, width = 12, height = 6)
check_model(best, panel = T, dot_size = 1, alpha = 0.15, line_size = 0.6,
            check = c("reqq"
            ))

#dev.off()

check_model(best, panel = T)
check_collinearity(best)
check_normality(best)
check_heteroscedasticity(best)

plot(best)
# residuals-vs-fitted plot with some outliers (Ben Bolker)
plot(best, id = 0.01)


#-----------------------------------
# Figure 1
#-----------------------------------

head(tags1moult)

f = dplyr::select(tags1moult, julian, moultjulian1NOV, moultjulian, state)
head(f)

b = dplyr::select(tags1moult, julian)
b$class = 'breeding'
head(b)

m = dplyr::select(tags1moult, moultjulian)
names(m) = 'julian'
m$class = 'moulting'
head(m)

d = rbind(b,m)
head(d)
tail(d)

d = d %>%
  group_by(julian) %>%
  mutate(freq = n()) %>%
  ungroup()

head(d)

# missing days
fullyear = 1:365
fullyear = data.frame(fullyear[!fullyear %in% d$julian])
names(fullyear) = 'julian'
head(fullyear)
fullyear$freq = 0
fullyear$class = NA
d = rbind(d,fullyear)
head(d)

# Remove any duplicates
d <- d %>%
  filter(duplicated(julian) == FALSE)

# Make the plot
p <- ggplot(d, aes(x= julian, y=freq, fill= class)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
  
  # This add the bars with a blue color
  geom_bar(stat="identity") +
  
  # Limits of the plot = very important. 
  # The negative value controls the size of the inner circle, 
  # the positive one is useful to add size over each bar
  ylim(-600,900) +
  
  # Custom the theme: no axis title and no cartesian grid
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-2,4), "cm")     # This remove unnecessary margin around plot
  ) +
  
  # This makes the coordinate polar instead of cartesian.
  coord_polar(start = 180)
p

p = p + guides(fill="none")

# # Save png to disk
# ggsave("./plots/Haulouts.png", device = 'png', bg = "white",
#        width = 50, height = 50, dpi = 800, units = 'in', limitsize = T)
# 
# ## Save Plot 
pdf("./plots/Figure1.pdf",
     useDingbats = FALSE, width = 5, height = 5)
 print(p)
 dev.off()

