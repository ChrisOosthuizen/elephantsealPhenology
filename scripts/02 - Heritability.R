# run 01 - Repeatability of breeding phenology
# until start of mixed model section

# now import mother-pup association data from supersmall records
ss <- read.csv("mom_pup.csv")  
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
  ggplot_theme_rrr() + 
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
             alpha = 0.85, shape = 15, size = 1.1) +
    scale_colour_manual(values = c("#4daf4a", "#984ea3"))+
  scale_y_continuous(breaks = c(seq(250,310, by = 10)), 
                     labels =c("07 Sep","17 Sep","27 Sep","07 Oct","17 Oct","27 Oct","06 Nov"),
                     limits = c(252, 305))+
  scale_x_continuous(breaks = c(seq(250,318, by = 10)), 
                     labels =c("07 Sep","17 Sep","27 Sep","07 Oct","17 Oct","27 Oct","06 Nov"),
                     limits = c(255, 305)) + 
  ggplot_theme_rrr() + 
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
pdf("./plots/2022_heritability_plot5.pdf",
    useDingbats = FALSE, width = 3.9, height = 3.9)
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

quantile(output.coef, c(0.025, 0.975))
quantile(output.pvalue, c(0.025, 0.975))

output.coef = as.data.frame(output.coef)
# Histogram with density instead of count on y-axis
ggplot(output.coef, aes(x=output.coef)) + 
  geom_histogram(aes(y=..density..),      
                 bins=50,
                 colour="black", fill="white")+
 # geom_density(alpha=.2, fill="blue") +
   ggplot_theme_rrr() +
  geom_vline(aes(xintercept=lm.fit$coefficients[2]), 
              color="blue", linetype="dashed", size=1)
  


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

# Overlaid histograms are created by setting the argument position="identity".
# We have also set the alpha parameter as alpha=.5 for transparency.
ggplot(gg, aes(x= output.coef, fill= booot)) + 
  geom_histogram(aes(y=..density..),      
                 bins=50,
                 colour="black",  alpha=0.4, position= "identity")+
  # geom_density(alpha=.2, fill="blue") +
  geom_vline(aes(xintercept=mean(output.coef)), color="blue", linetype="dashed", size=1)+
  ggplot_theme_rrr()







# =====================================================
# OLD - SCRAPS
# =====================================================

ma = left_join(ss, tags, by = c("momID" = "ID"))
str(ma)

# sanity check
length(ss$momID)
length(ma$momID)

## clean up slightly 
ma = ma %>% dplyr::select(c(momID,  pupID, wean.date, first.breeding.date, julian, year, age))
dim(ma)
ma
subset(ss, momID == "LB(2)001")
subset(ma, momID == "LB(2)001")
tail(subset(ma, momID == "LB(2)001"))

kind = left_join(ss, tags, by = c("pupID" = "ID"))
kind = kind %>% dplyr::select(c(momID,  pupID, wean.date, first.breeding.date, julian, year, age))
kind

left_join(ma, kind, by = c("momID" = "momID"))

str(ma)

