library(ggplot2)
library(cowplot)


# sliding window function
slideFunct <- function(data, window, step){
  total <- length(data)
  spots <- seq(from=1, to=(total-window), by=step)
  result <- vector(length = length(spots))
  for(i in 1:length(spots)){
    result[i] <- mean(data[spots[i]:(spots[i]+window)])
  }
  return(result)
}


# 

# read data using bjorns script
inp.tbl <- read.table(text = system("/data/bnf/scripts/plot_load.pl",intern=TRUE) )

# format data
colnames(inp.tbl)<-c("day","time","load")
inp.tbl$daytime <- paste(inp.tbl$day,inp.tbl$time)
#convert date into posix type
inp.tbl$daytime<-as.POSIXct(inp.tbl$daytime, format = "%Y-%m-%d %H:%M")


# add sliding average
inp.tbl$load[inp.tbl$load > 40] <- 40
inp.tbl$loadslide <- c(inp.tbl$load[1:6],slideFunct(inp.tbl$load,6))

inp.tblLastweek <- tail(inp.tbl,7*288)

inp.tblLastweek$day<-factor(inp.tblLastweek$day)



# plot
png("lennart_load.png",w=630,h=315)
ggplot(tail(inp.tbl,7*288),aes(x=daytime,y=load,col=day))+xlab("")+ theme(legend.position="none")+
  scale_colour_brewer(palette="Set2")+ylab("")+
  geom_line(aes(x=daytime,y=loadslide,col=day))+scale_x_datetime(date_breaks = "days",date_labels = "%d/%m")

dev.off()


system("scp lennart_load.png pi@10.0.224.47:/home/pi/seqmon")
