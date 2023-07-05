library(suncalc)
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

# Macquarie Island is not given in the data frame site name, so I can add the label 
# to the other side of the point in ggplot

day_lat = read.csv('./data/daylength.csv')

day_lat$peak = as.Date(day_lat$peak, format = "%Y/%m/%d")

# Plot latitude against peak female haulout

day_lat = day_lat[1:7 ,]

p.day_lat = ggplot(data = day_lat,
              aes(x = peak, y = lat, label=site)) + 
  gg_theme() + 
  ylab("Latitude") +
  xlab("Peak female breeding haulout") + 
  geom_smooth(method='lm', formula= y~x, color = "black", fill = "grey51") +
  geom_point(colour = "red", size = 3) + 
  geom_text(hjust=-0.1, vjust=-0.5)+
  scale_x_date(limits = as.Date(c('2021-10-05','2021-11-01')))+
  annotate("text", x=as.Date('2021-10-16'), y=-54, label= "Macquarie Is.") + 
  annotate("text", x=as.Date('2021-10-15'), y=-56, label= "North of the Polar Front", size = 3) +
  annotate("text", x=as.Date('2021-10-28'), y=-55.2, label= "South of the Polar Front", size = 3)+ 
  annotate("text", x=as.Date('2021-10-15'), y=-57, label= "Early breeding peak for latitude", size = 3) +
  annotate("text", x=as.Date('2021-10-28'), y=-56.2, label= "Late breeding peak for latitude", size = 3) +
  annotate("text", x=as.Date('2021-10-28'), y=-58.2, label= "South of the SBACC", size = 3) +
  annotate("text", x=as.Date('2021-10-28'), y=-59.2, label= "Early breeding peak for latitude", size = 3)

#annotate("text", x=as.Date('2021-10-28'), y=-58.2, label= "Southern Boundary of the Antarctic Circumpolar Current", size = 3) +
  
 # geom_hline(yintercept=-55, linetype="dashed", color = "grey51")
  p.day_lat

## Save Plot 
pdf("./plots/FigureS1.pdf",
    useDingbats = FALSE, width = 6, height = 5)
cowplot::ggdraw(p.day_lat)
dev.off()

