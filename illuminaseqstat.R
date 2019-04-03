
library(dplyr)
library(cowplot)
library(RColorBrewer)

setwd("/data/bnf/proj/seqmon/")

manad <- factor(c("jan","feb","mar","apr","maj","jun","jul","aug","sep","okt","nov","dec"))
manad <- factor(manad, levels=c("jan","feb","mar","apr","maj","jun","jul","aug","sep","okt","nov","dec"))


colorPal <- c("#7BC0F7", "#3B8AD9", "#F18226", "#FFDB69", "#61737B", "#A6B3B3", "#E24B26", brewer.pal(name = "Set1",9))

plotmachine <- function(seqstat, FlowCell_filter, expected_clusters,show_y_axis, show_x_axis, show_trendline ){
  
  # function to produce a line graph for a single machine
  # - seqstat = input data frame
  # - FlowCell_filter = which flowcell to plot
  # - expected_clusters = where to to put the gray line
  # - show_y_axis = TRUE/FALSE show the axis label
  # - show_x_axis = TRUE/FALSE show the axis label
  # - show_trendline = TRUE/FALSE show trendline
  
  filter.data <- dplyr::filter(seqstat,FlowCell==FlowCell_filter)
  
  # Use millions
  expected_clusters<-expected_clusters/1e6
  filter.data$PF_readcount<-filter.data$PF_readcount/1e6
  
  breaks_to_plot <- filter.data$Date
  breaks_to_plot <- as.character(seq(from=1,to=length(filter.data$Date)))
  
  if(length(breaks_to_plot)>50){
   breaks_to_plot[c(FALSE,TRUE,TRUE,TRUE)]<-""
  }
  
 
  filter.data<-  filter.data %>% group_by(machine) %>% mutate(id = row_number())
  
  g1<-ggplot(filter.data,aes(x=id,y=PF_readcount,col=machine,group=machine))+geom_point()+ylim(low=0,high=max(filter.data$PF_readcount)*1.1)+#geom_smooth(method = "lm",se = F)+
    geom_hline(yintercept=expected_clusters,col="grey",show.legend=T)+xlab("")+
    ylab("Clusters passing filter [M]")+ggtitle(FlowCell_filter) + scale_color_brewer(palette="Dark2")
    
  g1
  
  # if multiple lines in graph, show legend 
  if(length(unique(filter.data$machine)) > 1){
    g1<-g1+theme(legend.title=element_blank(),legend.position = "top")
  }else{
    g1 <- g1 + theme(legend.position="none")
  }
  
  if(!show_y_axis){
    g1<-g1+ylab("")
  }
  
  if(show_x_axis){
    g1<-g1+xlab("Run number")
  }
  
  
  if(show_trendline){
    g1<-g1+geom_smooth()
  }
  
  
  
  return(g1)
  
  
}


####################################################################
#
# Read data
#
####################################################################


#
# SEQUENCINGDATA
#


options(stringsAsFactors = F)
illumina.seqstat <- read.csv("/data/bnf/proj/seqmon/illumina_run_stats.tsv",sep="\t",header=T)
illumina.seqstat$Date <- substr(illumina.seqstat$run,1,6)

iontorrent.seqstat <- read.csv("/data/bnf/proj/seqmon/s5_seqmon.tsv",sep="\t",header=F)
iontorrent.seqstat$Date <- iontorrent.seqstat$V1


iontorrent.seqstat <- data.frame(machine="S5",run=iontorrent.seqstat$V1,raw_readcount=iontorrent.seqstat$V2, PF_readcount=iontorrent.seqstat$V2,
                                 yield=1,   Date=iontorrent.seqstat$Date)

seqstat<-rbind(illumina.seqstat,iontorrent.seqstat)

illumina.clusterdens <- read.csv("illumina_clusterdensity_stats.tsv",sep="\t")

#
# Cluster density



#
# SAMPLEDATA
#

sampleData <- read.csv(paste0("/data/bnf/proj/clarityReport/data/production/sample_data.fixed.tsv"),sep="\t",row.names=NULL)
sampleData$Date<-as.Date(sampleData$Date)
sampleData$Week <- format(as.POSIXct(sampleData$Date), format = "%W-%Y")
sampleData$Week2 <- format(as.POSIXct(sampleData$Date), format = "%Y\n%W")
sampleData$Month <- format(as.POSIXct(sampleData$Date), format = "%m-%Y")
sampleData$Month2 <- format(as.POSIXct(sampleData$Date), format = "%Y\n%m")


sampleData<-sampleData[sampleData$Tissue!="Benmärg",]
sampleData<-sampleData[!duplicated(sampleData),]
sampleData<-sampleData[sampleData$Classification == "Rutinprov",]
sampleData<-sampleData[sampleData$Department != "TEST",]
sampleData<-sampleData[sampleData$Project    != "The_TEST_Project",]

sampleData$Workflow<-plyr::revalue(sampleData$Analysis,  c(#myeloisk panel
                                                           "Myeloisk Panel - Parad"="Myeloisk Panel", 
                                                           "Myeloisk Panel - Oparad"="Myeloisk Panel",
                                                           "Myeloisk Panel - Oparad - KLL"="Myeloisk Panel",
                                                           "Myeloisk Panel - Oparad - MPN"="Myeloisk Panel",
                                                           "Myeloisk Panel - Oparad - AML"="Myeloisk Panel",		
                                                           #tumörexom
                                                           "SureSelectXTHS - Paired Tumor Exome"="Övriga",
                                                           "SureSelect XTHS Clinical Exome - Tumor Paired Exome"="Övriga",
                                                           "SureSelectXTHS - Unpaired Tumor Exome" = "Övriga",
                                                           "SureSelectXTHS - Single Exome"="Exom",
                                                           "SureSelectXTHS - Trio Exome"="Exom",
                                                           "TruSeq Stranded mRNA"="RNA-seq",
                                                           "TruSeq Stranded mRNA - Bladder"="RNA-seq"
                                                           ,"TruSeq Stranded mRNA - Fusion"="RNA-seq",
                                                           "AmpliSeq CancerHotspot"="AmpliSeq Cancer", 
                                                           "AmpliSeq ColonLung" = "AmpliSeq Cancer", 
                                                           "Oncomine Focus Assay" = "Oncomine Cancer",
                                                           "Clarigo NIPT Analys" = "NIPT",
                                                           "AmpliSeq CF"="Övriga"),
                                   warn_missing = F)



sampleData[grep("Myeloisk",sampleData$Workflow),"Workflow"]<-"Myeloisk Panel"
sampleData[grep("SWEA",sampleData$Workflow),"Workflow"]<-"BRCA"
sampleData[grep("BRCA",sampleData$Workflow),"Workflow"]<-"BRCA"
sampleData[grep("Liquid Biopsy",sampleData$Workflow),"Workflow"]<-"Oncomine Liquid Biopsy"




#
# TATDATA
#

tatData <- read.csv(paste0("/data/bnf/proj/clarityReport/data/production/tat_data.fixed.tsv"),sep="\t",row.names=NULL)
tatData <-filter(tatData,Sample.classification=="Rutinprov - Godkänt")
tatData<-merge(sampleData,tatData,by="Sample",all.x=F)
tatData$TAT <- as.integer(as.Date(tatData$Date_process)-as.Date(tatData$Date))




####################################################################
#
# Prepare data
#
####################################################################


seqstat$FlowCell<-as.character(seqstat$machine)

seqstat[(seqstat$machine=="NextSeq1" | seqstat$machine=="NextSeq2")  & seqstat$PF_readcount > 3e08,"FlowCell"]<-"NextSeq High"
seqstat[(seqstat$machine=="NextSeq1" | seqstat$machine=="NextSeq2")  & seqstat$PF_readcount < 3e08,"FlowCell"]<-"NextSeq Mid"
seqstat[seqstat$machine=="S5"  & seqstat$PF_readcount > 10e6,"FlowCell"]<-"S5 - 540"



####################################################################
#
# Plot 1: line graphs per machine
#
####################################################################


plot_grid(plotmachine(seqstat,"MiniSeq",2.5e07,TRUE,FALSE,FALSE),plotmachine(seqstat,"MiSeq",2.5e07,FALSE,FALSE,FALSE), 
          plotmachine(seqstat,"NextSeq Mid",1.3e08,TRUE,TRUE,TRUE),plotmachine(seqstat,"NextSeq High",4e08,FALSE,TRUE,TRUE) , nrow = 2)

ggsave("illumina_sequencing_output.png",width = 13.4,height = 7.24)



####################################################################
#
# Plot 2: TAT bar graph
#
####################################################################

#head(tatData)


tatData$Analysis<-plyr::revalue(tatData$Analysis, c(#myeloisk panel
  "Myeloisk Panel - Parad"="Myeloisk Panel", 
  "Myeloisk Panel - Oparad"="Myeloisk Panel",
  "Myeloisk Panel - Oparad - KLL"="Myeloisk Panel",
  "Myeloisk Panel - Oparad - MPN"="Myeloisk Panel",
  "Myeloisk Panel - Oparad - AML"="Myeloisk Panel",		
  #tumörexom
  "SureSelectXTHS - Paired Tumor Exome"="Övriga",
  "SureSelect XTHS Clinical Exome - Tumor Paired Exome"="Övriga",
  "SureSelectXTHS - Unpaired Tumor Exome" = "Övriga",
  "SureSelectXTHS - Single Exome"="Exom",
  "SureSelectXTHS - Trio Exome"="Exom",
  "TruSeq Stranded mRNA"="RNA-seq",
  "TruSeq Stranded mRNA - Bladder"="RNA-seq"
  ,"TruSeq Stranded mRNA - Fusion"="RNA-seq",
  "AmpliSeq CancerHotspot"="AmpliSeq Cancer", 
  "AmpliSeq ColonLung" = "AmpliSeq Cancer", 
  "Oncomine Focus Assay" = "Oncomine Cancer",
  "Clarigo NIPT Analys" = "NIPT",
  "AmpliSeq CF"="Övriga"))

a1<-aggregate(TAT~Week2+Analysis,tatData,mean)

a1 <- filter(a1,TAT <50)
a1 <- filter(a1,Analysis !=  "Övriga")
a1$Week2 <- substr(a1$Week2,3,8)
a1[a1$Workflow=="Myeloisk Panel" & a1$TAT>15 ,"TAT"]<-NA


labels <- a1$Week2
labels[-seq(0,length(labels),by=5)]<-""


tatplot<-ggplot(a1,aes(x=Week2,y=TAT,fill=Analysis))+geom_bar(stat = "identity")+scale_fill_manual(values = colorPal[c(1,3,5,7)])+
  ylab("Svarstid [dagar]")+xlab("Vecka")+  theme(legend.title=element_blank(),legend.position = "top")+
  scale_x_discrete(labels=labels)
tatplot


ggsave("tat.png",w=8,h=9)




####################################################################
#
# Plot 3A: Workflows per month
#
####################################################################

#workflow_month<-reshape::melt(table(sampleData$Month,sampleData$Workflow))
#workflow_month$manad <- manad[as.integer(substr(workflow_month$Var.1,1,2))]
#workflow_month$År <- substr(workflow_month$Var.1,4,8)
#g1<-ggplot(workflow_month,aes(x=manad,y=value,fill=År,group=ar))+geom_bar(stat="identity")+facet_wrap(~Var.2,scales="free")+
#  xlab("")+ylab("Inkomna prover")

sampleData.cmd<-sampleData[sampleData$Department == "Klinisk Genetik" | 
                             sampleData$Department == "Klinisk Patologi" |
                             sampleData$Department == "Sahlgrenska Klinisk Genetik" |
                             sampleData$Department == "CMD - Diana Karpman" |
                             sampleData$Department == "Klinisk immunologi och transfusionsmedicin" ,
                         ]

workflow_month<-reshape::melt(table(sampleData.cmd$Month2,sampleData.cmd$Workflow))
workflow_month$Var.1<-substr(workflow_month$Var.1,3,12)

g1<-ggplot(workflow_month,aes(x=Var.1,y=value,fill=Var.2))+ geom_hline(yintercept=c(50,100,150,200,250,300,350,400,450,500,550), linetype="dotted")+geom_bar(stat="identity")+
  scale_fill_manual(values = colorPal)+theme(legend.title=element_blank(),legend.position = "top")+xlab("")+ylab("Inkomna prover")

g1

ggsave("workflows_per_month.png")



####################################################################
#
# Plot 3B: Workflows per month
#
####################################################################


options(stringsAsFactors = F)
workflow_month$Year<-paste0("20",substr(workflow_month$Var.1,1,2))
workflow_month$Month<-substr(workflow_month$Var.1,4,5)
workflow_month<-workflow_month[workflow_month$value>0,]

workflow_month$Var.2<-as.character(workflow_month$Var.2)
workflow_month[workflow_month$Var.2=="Oncomine Cancer","Var.2"]<-"AmpliSeq Cancer"

workflow_month<-workflow_month[workflow_month$Var.2 == "AmpliSeq Cancer" | workflow_month$Var.2 == "Myeloisk Panel" | workflow_month$Var.2 == "NIPT" | workflow_month$Var.2 == "Exom",]


labels <- workflow_month$Var.1
labels[seq(0,length(labels),by=2)]<-""

ggplot(workflow_month,aes(x=Var.1,y=value,fill=Year,group=Year))+facet_wrap(~Var.2,scales="free_y")+
  geom_bar(stat="identity")+  scale_fill_manual(values = colorPal)+theme(legend.title=element_blank(),legend.position = "top")+xlab("")+ylab("Inkomna prover")+scale_color_manual(values = colorPal[c(1,3,5,7)])+
  scale_x_discrete(labels=labels)

ggsave("samples_per_month_per_analysis.png",w=16,h=8,dpi = 300)



####################################################################
#
# Plot 4: combine plot 2 & 3
#
####################################################################


plot_grid(g1,tatplot,rel_widths = c(1.3,1))
ggsave("workflows_per_month_TAT.png",scale=1.4, w=9, h=5)



####################################################################
#
# Plot 5: cluster density
#
####################################################################


seqstat.clusterdens<-merge(seqstat,illumina.clusterdens,by.x="Date",by.y="Datum")

seqstat.clusterdens<-filter(seqstat.clusterdens,FlowCell == "NextSeq High" | FlowCell == "NextSeq Mid")

g1<-ggplot(seqstat.clusterdens,aes(x=Kluisterdensittet,y=PF_readcount/1e6,col=FlowCell))+geom_point()+scale_color_manual(values = brewer.pal("Blues",n=9)[c(5,7,9)])+
  theme(legend.title=element_blank(),legend.position = "none")+xlab("Klusterdensitet [K/mm2]")+ylab("Godkända kluster [M]")+facet_wrap(~FlowCell,scales="free")+geom_smooth(n=10)
g1

ggsave("klusterdensitet.png",scale=1.4, w=9, h=5)




####################################################################
#
# Transfer to seqmon
#
####################################################################

system("scp illumina_sequencing_output.png pi@10.0.224.47:/home/pi/seqmon")
system("scp workflows_per_month.png pi@10.0.224.47:/home/pi/seqmon")
system("scp samples_per_month_per_analysis.png pi@10.0.224.47:/home/pi/seqmon")
system("scp tat.png pi@10.0.224.47:/home/pi/seqmon")
system("scp workflows_per_month_TAT.png pi@10.0.224.47:/home/pi/seqmon")
system("scp klusterdensitet.png pi@10.0.224.47:/home/pi/seqmon")
