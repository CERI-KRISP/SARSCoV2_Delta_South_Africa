
library(ggplot2)
library(ape)
library(repr)
library("readxl")
library('gridExtra')
#library(tidyverse)
library(dplyr)
library(hrbrthemes)
library(ggpubr)
library(cowplot)
library(ggthemes)
library(viridis)
library(ggrepel)
library("ggsci")
library(ggalt)
library("Hmisc")
library("scales")

require(tidyverse)
library("readxl")
library("lubridate")
require("ggalt")
library("wesanderson")



library(ggplot2)
library(ape)
library(repr)
library("readxl")
library('gridExtra')
#library(tidyverse)
library(dplyr)
library(hrbrthemes)
library(ggpubr)
library(cowplot)
library(ggthemes)
library(viridis)
library(ggrepel)
library("ggsci")
library(ggalt)
library("Hmisc")
library("scales")




# load pacakges
library(ggtree)
library(tidyverse)
library(tidytree)
library(ape)
library(treeio)

#PLOTS MAY BE FURTHER EDITTED IN ILLUSTRATOR/POWERPOINT TO ADD LABELS OR CHANGE LAYOUTS

######################################
########### FIG 1 ####################
######################################

###### PANEL 1A


data3<-read_excel('SA_provincial_daily_cases.xlsx')

data3$days<-as.Date(cut(data3$date,
                        breaks = "day",
                        start.on.monday = FALSE))

data3$date<-as.Date(cut(data3$date,
                        breaks = "week",
                        start.on.monday = FALSE))


data3 <- data3 %>% 
  dplyr::mutate(total_daily_7day = zoo::rollmean(SA, k = 7, fill = NA)) %>% 
  dplyr::ungroup()


R_estimates<-read_excel('Re.xlsx')
R_estimates$date<-as.Date(cut(R_estimates$date,
                              breaks = "day",
                              start.on.monday = FALSE))

R_estimates_SA=subset(subset(subset(R_estimates, region=='ZAF'),data_type=="Confirmed cases"),estimate_type=="Cori_slidingWindow")

excess_deaths<-read_excel('Excess_deaths.xlsx')
excess_deaths$date<-as.Date(cut(excess_deaths$date,
                                breaks = "week",
                                start.on.monday = FALSE))

p1A<-ggplot()+
  theme_excel_new()+
  geom_hline(yintercept=10000, color='cyan4', linetype=2, size=0.5,alpha=0.7) +
  geom_line(data = R_estimates_SA, aes(x = date, y = median_R_mean*10000, color = "Re"), size=0.5) +
  geom_ribbon(data = R_estimates_SA,aes(x=date, ymin=median_R_lowHPD*10000, ymax=median_R_highHPD*10000), fill='cyan4', alpha=0.2) +
  #geom_line(data=excess_deaths,aes(x = date,y=RSA,color='Deaths'), size=1)+ ### Uncomment this to add death data
  geom_line(data=data3,aes(x = days,y=total_daily_7day,color='Cases'), size=1)+
  theme(axis.title.x = element_blank())+
  theme(axis.title.y = element_text(color="black", size=15, face="bold"))+
  theme(axis.text.y = element_text(color="black", size=12))+
  theme(axis.text.x = element_text(color="black", size=10))+
  scale_x_date(date_labels = "%b\n%Y",date_breaks = "2 month",limits=as.Date(c("2020/08/01","2021/08/15")))+
  #scale_colour_manual(values=c('slategray4', 'indianred3','cyan4'), name="")+ ### Use this if adding death data
  scale_colour_manual(values=c('slategray4','cyan4'), name="")+
  xlab(' ')+
  theme(legend.position = 'bottom')+
  theme(legend.title = element_text(size=12))+
  theme(legend.text = element_text(size=12))+
  
  theme(plot.title = element_text(hjust = 0.5))+
  scale_y_continuous(
    
    # Features of the first axis
    name = "Daily Cases",
    
    # Add a second axis and specify its features
    sec.axis = sec_axis(~./10000, name="Re"), limits=c(0,2.2*10000)
    
  )+
  ggtitle('Epidemiological Modeling in South Africa')

p1A

###### PANEL 1B


data2<-read_excel('SA_all_data_2Sept_gooddates.xlsx')
data2$Nextstrain_variants<-factor(data2$Nextstrain_variants,levels = c("Other Lineages","20H (Beta, V2)","20I (Alpha, V1)","21A (Delta)","21B (Kappa)","C.1.2"))

data2$division<-factor(data2$division,levels = c("Mpumalanga","North West","Northern Cape","Free State","Eastern Cape","Limpopo","Western Cape","Gauteng","KwaZulu-Natal"))

data2$days<-as.Date(cut(data2$date,breaks = "day",start.on.monday = FALSE))
data2$date<-as.Date(cut(data2$date,breaks = "week",start.on.monday = FALSE))
data2$date2<-as.Date(cut(data2$date,breaks = "2 week",start.on.monday = FALSE))
data2<- data2 %>% filter(division!="South Africa")

myColor3 <- randomcoloR::distinctColorPalette(k = 13)


### Extended Data Figure showing distribution of Delta lineages
P_pangolin_delta<-data2 %>%
  mutate(date=as.POSIXct(date)) %>% #convert date to date
   ggplot(data=subset(subset(data2, Nextstrain_variants=='21A (Delta)'), pango_lineage!='None'),
         mapping = aes(x = date,fill=pango_lineage))+
  geom_bar(position='fill',width=5, size=0.2, color='black')+
  theme_classic()+
  theme(axis.text.x = element_text(color="black", size=16))+
  xlab("Sampling Date")+ 
  scale_x_date(date_labels = "%b\n%Y",date_breaks = "1 month",limits=as.Date(c("2021/03/01","2021/08/15")))+
  theme(axis.title.x = element_text(color="black", size=12, face="bold"))+
  theme(axis.text.x = element_text(color="black", size=8))+
  theme(axis.title.y = element_text(color="black", size=12, face="bold"))+
  theme(axis.text.y = element_text(color="black", size=12))+
   scale_fill_manual(values = myColor3, name='Delta Sublineages')+
   theme(legend.text = element_text(size=11))+
  theme(legend.title = element_text(size=12))+
  xlab('Date')+
  ylab('Proportion of Genomes')+
  ggtitle('Delta Lineages in South Africa')

P_pangolin_delta



P1b<-ggplot(data=subset(data2, !is.na(Nextstrain_variants)), aes(x = days,fill=Nextstrain_variants))+
  geom_area(stat='bin', position='fill', binwidth=3)+
  theme_classic()+
  theme(axis.text.x = element_text(color="black", size=16))+
  xlab("Sampling Date")+ 
  scale_x_date(date_labels = "%b\n%Y",date_breaks = "2 month",limits=as.Date(c("2020/08/01","2021/08/15")))+
  theme(axis.title.x = element_blank())+
  theme(axis.text.x = element_text(color="black", size=10))+
  theme(axis.title.y = element_text(color="black", size=12, face="bold"))+
  theme(axis.text.y = element_text(color="black", size=10))+
   scale_fill_manual(values=c('grey90','aquamarine3','dodgerblue3','indianred3','gold3'), name='Variants')+
  theme(legend.text = element_text(size=10))+
  theme(legend.title = element_text(size=10))+
  theme(legend.position = "bottom")+
  theme(plot.margin = unit(c(2,2,0,0), "lines"))+
  ylab('Proportion of Genomes')+
  scale_y_continuous(labels = scales::percent)+
  ggtitle('SARS-CoV-2 Genomes in South Africa')+
  theme(plot.title = element_text(hjust = 0.5))

P1b


### Extended Data Figure showing the absolute number of genomes samples in South Africa, coloured by variants
P_variants_absolute<-ggplot(data=subset(data2, !is.na(Nextstrain_variants)), aes(x = days,fill=Nextstrain_variants))+
  geom_area(stat='bin', binwidth=6)+
  theme_classic()+
  xlab("Sampling Date")+ 
  scale_x_date(date_labels = "%b\n%Y",date_breaks = "2 month",limits=as.Date(c("2020/08/01","2021/08/15")))+
  theme(axis.title.x = element_blank())+
  theme(axis.text.x = element_text(color="black", size=10))+
  theme(axis.title.y = element_text(color="black", size=12, face="bold"))+
  theme(axis.text.y = element_text(color="black", size=10))+
   scale_fill_manual(values=c('grey90','aquamarine3','dodgerblue3','indianred3','gold3'), name='Variants')+
  theme(legend.text = element_text(size=10))+
  theme(legend.title = element_text(size=10))+
  theme(legend.position = "bottom")+
  theme(plot.margin = unit(c(2,2,0,0), "lines"))+
  ylab('Genome Count')+
  ggtitle('SARS-CoV-2 Genomes in South Africa')+
  theme(plot.title = element_text(hjust = 0.5))

P_variants_absolute

##### Panel 1C - Additional Scripts Provided (maps.R, cases_maps.R)

##### Panel 1D - Additional Script Provided (variants_multinomial.R) for which the data can be downloaded from GISAID according the the instructions in variants_multinomial_data_download_GISAID.png


######################################
########### FIG 2 ####################
######################################



# load pacakges
library(ggtree)
library(tidyverse)
library(tidytree)
library(ape)
library(treeio)

# read in the tree and metadata
tree<-read.tree('timetree.nwk')
#The following metadata is annoted with cluster info obtained with Phylotype (pbm folder)
metadata_df <- read_excel('SADelta_5k_Delta_ref_5k_metadata_no_duplicates_n10582_clusters.xlsx')
# transform dates
metadata_df$date<-as.Date(cut(metadata_df$date,
                              breaks = "week",
                              start.on.monday = FALSE))


p2A<-ggtree(tree, mrsd="2021-08-20", as.Date=TRUE,color='grey80',show.legend = FALSE,size=0.2) %<+% metadata_df + theme_tree2() +
  scale_colour_manual(values=c("grey80",'red2')) +
  scale_fill_manual(values=c("red4",'red2'))+
 scale_x_date(date_labels = "%b-%Y",date_breaks = "2 month") +
  theme(axis.text=element_text(size=10)) +
  theme(axis.text.x = element_text(size=6,hjust = 1,vjust=0.5, angle=90))+
   geom_tippoint(size = .2,color='grey60',show.legend = FALSE) +
   geom_tippoint(aes(
    subset=(grepl('SouthAfrica',label,fixed=TRUE)==TRUE)), fill='red2',size=2, stroke=0.8, shape=21, alpha=1, color='black') +
  geom_tippoint(aes(subset=(grepl('SouthAfrica',label,fixed=TRUE)==TRUE)),size = 1,color='red2',show.legend = FALSE, alpha=0.4) +
   geom_tippoint(aes(subset=(grepl('SouthAfrica',label,fixed=TRUE)==TRUE) & cluster=='11637'), fill='lightpink2',stroke=0.8, shape=21,size = 2,alpha=1,show.legend = FALSE,color='black') +
  geom_tippoint(aes(subset=(grepl('SouthAfrica',label,fixed=TRUE)==TRUE) & cluster=='11637'),size = 1,color='lightpink2',show.legend = FALSE, alpha=0.4) +
  geom_tippoint(aes(subset=(!grepl('SouthAfrica',label,fixed=TRUE)==TRUE)),size = 0.6,color='grey70',show.legend = FALSE, alpha=0.4) +
  theme(
    legend.position      = "top",
    legend.direction = 'horizontal',
    legend.justification = c(0,1),
    legend.key.size      = unit(.5,'lines'),
    legend.title         = element_text(size = 6),  
    legend.text          = element_text(size = 6),  
    axis.text            = element_text(size = 6.4),
    legend.spacing.y     = unit(-.2,"cm"),
    legend.background    = element_rect(fill = alpha("white",0)))+
  expand_limits(y = 11000)
p2A

### Extended Data Figure showing the ML tree, with monophyletic clusters highlighted in pink including the biggest one (Cluster A) in darker pink

p<-ggtree(tree, mrsd="2021-08-20", as.Date=TRUE,color='grey80',show.legend = FALSE,size=0.2) %<+% metadata_df + theme_tree2() +
  scale_colour_manual(values=c("grey80",'red2')) +
  scale_fill_manual(values=c("red4",'red2'))+
  scale_x_date(date_labels = "%b-%Y",date_breaks = "2 month") +
  theme(axis.text=element_text(size=10)) +
  theme(axis.text.x = element_text(size=6,hjust = 1,vjust=0.5, angle=90))+
  geom_tippoint(size = .2,color='grey60',show.legend = FALSE) +
   geom_tippoint(aes(
     subset=(grepl('SouthAfrica',label,fixed=TRUE)==TRUE)), fill='red2',size=2, stroke=0.8, shape=21, alpha=1, color='black') +
  geom_tippoint(aes(subset=(grepl('SouthAfrica',label,fixed=TRUE)==TRUE)),size = 1,color='red2',show.legend = FALSE) +
  
  geom_tippoint(aes(subset=(grepl('SouthAfrica',label,fixed=TRUE)==TRUE) & !is.na(cluster)), fill='lightpink2',stroke=0.8, shape=21,size = 2,alpha=1,show.legend = FALSE,color='black') +
   geom_tippoint(aes(subset=(grepl('SouthAfrica',label,fixed=TRUE)==TRUE) & !is.na(cluster)),size = 1,color='lightpink2',show.legend = FALSE) +
  geom_tippoint(aes(subset=(grepl('SouthAfrica',label,fixed=TRUE)==TRUE) & cluster=='11637'), fill='lightpink4',stroke=0.8, shape=21,size = 2,alpha=1,show.legend = FALSE,color='black') +
  
  geom_tippoint(aes(subset=(grepl('SouthAfrica',label,fixed=TRUE)==TRUE) & cluster=='11637'),size = 1,color='lightpink4',show.legend = FALSE) +
  theme(
    legend.position      = "top",
    legend.direction = 'horizontal',
    legend.justification = c(0,1),
    legend.key.size      = unit(.5,'lines'),
    legend.title         = element_text(size = 6),  
    legend.text          = element_text(size = 6),  
    axis.text            = element_text(size = 6.4),
    legend.spacing.y     = unit(-.2,"cm"),
    legend.background    = element_rect(fill = alpha("white",0)))+
  expand_limits(y = 11000)
p
#ggsave('Delta_tree_v2_9.pdf', width = 2.5, height = 15, units = "cm",limitsize = FALSE)
#ggsave('Delta_bigtree_v2_15.pdf', width = 15, height = 350, units = "cm",limitsize = FALSE)




##### Panel 2B
library(lubridate)
importexport<-read_excel('Delta_ImportExport_v2.xlsx') ### Data generated from Mugration analysis. Script and input files (timetree in nexus format) available in folder: AncestralStateReconstruction.

importexport$EventTime<-as.numeric(importexport$EventTime)


importexport$date <- as.Date(format(date_decimal(importexport$EventTime), "%Y-%m-%d"))

importexport$date


importexport$days<-as.Date(cut(importexport$date,
                              breaks = "day",
                              start.on.monday = FALSE))

importexport$date<-as.Date(cut(importexport$date,
                              breaks = "2 week",
                              start.on.monday = FALSE))
importexport$date2<-as.Date(cut(importexport$date,
                               breaks = "1 month",
                               start.on.monday = FALSE))


p2B<-importexport %>%
  mutate(date=as.POSIXct(date)) %>%
  ggplot()+
  geom_area(stat='bin',data=subset(data2,Nextstrain_variants=="21A (Delta)"), mapping = aes(x = days, fill='All SA Delta Genomes'),color='black', size=0.4,binwidth=15,position = position_stack(reverse=TRUE))+
  geom_area(stat='bin',binwidth=10,data=subset(importexport, Destination=='South Africa'), mapping = aes(x = days,fill='Inferred Introduction Events'),color='black',position = position_stack(reverse=TRUE))+
  theme_classic()+
  theme(axis.text.x = element_text(color="black", size=16))+
  theme(axis.title.x = element_blank())+
  theme(axis.text.x = element_text(color="black", size=10))+
  theme(axis.title.y = element_text(color="black", size=12, face="bold"))+
  theme(axis.text.y = element_text(color="black", size=10))+
  scale_fill_manual(values=c('coral3','grey80'))+
  scale_x_date(date_labels = "%b\n%Y",date_breaks = "1 month", limits=as.Date(c('2021-03-01','2021-08-01')))+
   theme(legend.text = element_text(size=10))+
  theme(legend.title = element_text(size=14))+
  theme(plot.margin = unit(c(2,2,0,0), "lines"))+
  theme(legend.position = "bottom")+
  labs(fill=" ")+
  xlab('Date')+
  ylab('Counts (log)')+
  scale_y_log10()
pane2B


library(sf)
library(raster)
library(spData)
library(tmap)
library(leaflet)
library(cartogram)
library(ggnewscale)

### Extended Data Figure showing the inferred origins of the Delta variant introduced into South Africa


df_imports <- subset(importexport, Destination=='South Africa') %>% count(Origin)




world = world %>% 
  left_join(df_imports, by = c("name_long" = "Origin"))



panelA<-ggplot(subset(world,!is.na(pop))) +
  geom_sf(aes(geometry = geom),color='grey95', fill = 'grey95')+
  geom_sf(aes(geometry = geom, fill=n),color='black', size=0.1)+
  
  #geom_sf(data=subset(africa, name=='Mauritania'),aes(geometry = geom),color='grey40', fill = 'white')+
  #geom_polygon(data=subset(africa, name=='Western Sahara'), aes(geometry = geom), colour="red")+
  
  #scale_x_continuous(trans = log2_trans())+
  scale_fill_distiller(palette = "OrRd", direction = 1,na.value = "white",name='Inferred\nIntroduction\nEvents', trans='log',breaks = c(32,16,8,4,2), labels = c(32,16,8,4,2)) + 
  # scale_fill_gradient(name = "count", trans = "log",
  #                     breaks = my_breaks, labels = my_breaks)
  #geom_sf_label(aes(label = Count),alpha=0.5, size=3,label.padding = unit(0.1, "lines"))+
  #geom_sf_text(data=subset(subset(africa,name!='Democratic Republic of the Congo'), area_km2>50000),aes(label=ifelse(Count>10,name,''), size=area_km2),alpha=0.85,vjust=2, color='black', show.legend = FALSE)+
  #geom_sf_text(data=subset(subset(africa,name=='Democratic Republic of the Congo'), area_km2>50000),aes(label=ifelse(Count>10,'DRC',''), size=area_km2),alpha=0.85,vjust=2, color='black', show.legend = FALSE)+
  
  theme(legend.position = 'bottom')+
  # scale_size_continuous(range = c(2, 4))+
  theme_void()
panelA 

##### Panel 2C - Additional Script Provided (Fig2C_Seraphim.R in BEAST folder)

library("readxl")
library(lubridate)
library("scales")
library(ggplot2)


ggplotRegression <- function(fit){
  
  ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
    geom_point() +
    stat_smooth(method = "lm", col = "red") +
    labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                       #"Intercept =",signif(fit$coef[[1]],5 ),
                       "R = ",signif(summary(fit)$coef, 5),
                       " Slope =",signif(fit$coef[[2]], 5),
                       " P =",signif(summary(fit)$coef[2,4], 5)))
}

###Extended Data Figure showing the tempest plot for the Delta ML tree

tempest_data<-read_excel('SA_Delta_tree_tempest.xlsx')
tempest_data$date2<-date_decimal(tempest_data$date)

tempest_data$date2<-as.Date(cut(tempest_data$date2,
                                breaks = "day",
                                start.on.monday = FALSE))

ggplotRegression(lm(distance ~ date2, data = tempest_data))


p_tempest<-ggplot(tempest_data, aes(date2,distance))+
  geom_point(color='black',size=1.5,shape=21, stroke=0.3, aes(fill=Location), alpha=0.8)+
  geom_smooth(method=lm,se=T, color='black')+
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.text=element_text(size=10))+
  theme(axis.text.x =element_text(size=10))+
  theme(axis.title=element_text(size=10))+
  scale_x_date(date_labels = "%b\n%Y",breaks='2 month')+
  scale_fill_manual(values=c('grey75','red2'), name='Location')+
  ylab("Root-to-tip Distance")+
  xlab("Date")+
  scale_y_continuous(labels = scientific)+
  annotate(geom="text", x=as.Date("2020/12/01"),y=0.0017,label="r = 0.3",
           color="black",size=4)+
  annotate(geom="text", x=as.Date("2020/12/01"),y=0.00155,label="r2 = 0.09",
           color="black",size=4)+
  ggtitle('Delta tree')

p_tempest



###Extended Data Figure showing the tempest plot for Cluster A

tempest_data<-read_excel('clusterA_n286_tempest.xlsx')
tempest_data$date2<-date_decimal(tempest_data$date)

tempest_data$date2<-as.Date(cut(tempest_data$date2,
                                breaks = "day",
                                start.on.monday = FALSE))


ggplotRegression(lm(distance ~ date2, data = tempest_data))


p_tempest<-ggplot(tempest_data, aes(date2,distance))+
  geom_point(color='black',size=3,shape=21, stroke=0.3, fill='lightpink2')+
  geom_smooth(method=lm,se=T, color='black')+
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.text=element_text(size=10))+
  theme(axis.text.x =element_text(size=10))+
  theme(axis.title=element_text(size=10))+
  scale_x_date(date_labels = "%b\n%Y",breaks='1 month')+
  scale_fill_manual(values=c('grey75','red2'), name='Location')+
  ylab("Root-to-tip Distance")+
  xlab("Date")+
  scale_y_continuous(labels = scientific)+
  annotate(geom="text", x=as.Date("2021/06/01"),y=0.0004,label="r = 0.3",
           color="black",size=4)+
  annotate(geom="text", x=as.Date("2021/06/01"),y=0.00038,label="r2 = 0.09",
           color="black",size=4)+
  ggtitle('Cluster A - n=286')

p_tempest

