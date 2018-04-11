# Draw a lin graph from strainLog.txt
# 

library(zoo)
library(ggplot2)
library(reshape2)

df = read.table('./strainLog.txt', header = F, stringsAsFactors = F)
df$V1 = NULL
titles = c('Date', 'Total', gsub(":[0-9]+$", "", df[1,])[-c(1,2)])
df <- as.data.frame(lapply(df, FUN = function(x) as.numeric(gsub("[0-9]+_[0-9]+:", "", x))))
colnames(df) = titles

percent = (df/df$Total)*100
percent$Date = df$Date
percent$Total = NULL

smooth_factor = 50

smooth_percent = as.data.frame(rollmean(percent,smooth_factor))

# Set average percent threshold for organisms to display
smooth_percent_filtered = smooth_percent[,!colMeans(smooth_percent) < 0.5]

# Somehow last line doesn't work for me in Ubuntu
#
# smooth_percent = percent[smooth_factor:nrow(percent),]
# for (i in rownames(smooth_percent)) {
#   i = as.numeric(i)
#   for (j in colnames(smooth_percent)[2:ncol(smooth_percent)]){
#     smooth_percent[i,j] = mean(percent[(i-smooth_factor):i,j])
#   }
# }
# smooth_percent$Date = percent$Date
# smooth_percent = rbind(percent[1:(smooth_factor-1),],smooth_percent)

smooth_percent_long <- melt(smooth_percent_filtered, id="Date")

ggplot(data=smooth_percent_long,
       aes(x=Date, y=value, colour=variable)) +
  ylab('Percent') +
  geom_line()

ggsave("dynamics.png", plot = last_plot())
