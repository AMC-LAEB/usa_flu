library(scales)
library(stringr)
library(sleekts)
library(ape)

setwd("/Users/simondejong/US_phylo")
source("plotFunctions.R")

ili = get_incidence_data()
dat = getClusterData_subsamp("time_0.083_pct_0.1_2", 0.05)

mtdt = dat[[1]]
cluster_df = dat[[2]]
cluster_df_country = dat[[3]]

cluster_df_country$Time_Since_Onset = as.numeric(cluster_df_country$RelativeOnset_Date)-cluster_df_country$EpidemicOnset
cluster_df_country$Time_Since_Sample = cluster_df_country$FirstSamp-cluster_df_country$EpidemicOnset
cluster_df_country$Time_Since_TMRCA = cluster_df_country$TMRCA-cluster_df_country$EpidemicOnset

thm = theme(axis.text = element_text(color="black",size=7),
            panel.spacing.x = unit(1, "mm"),
            axis.line=element_line(color='black',size = 0.2),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            legend.position='none',
            axis.title= element_text(size=7),
            strip.background=element_rect(colour="black",fill="white"),
            strip.text.x = element_text(size = 7),
            legend.text=element_text(size=7),
            legend.title =element_text(size=7),
            legend.key=element_rect(fill="white"),
            axis.ticks = element_line(size=0.1),
            legend.key.width = unit(.5, "line"),
            legend.spacing.y = unit(.001, 'cm'),
            legend.box.spacing = unit(0, "pt"),
            legend.margin=margin(1,1,1,1))

plotDistribution <- function(){
  prop_df = data.frame(subtype=NA,season=NA,seqprop=NA,clprop=NA,ncluster=NA)
  for (subtype in unique(cluster_df_country$Subtype)){
    for (season in unique(cluster_df_country[cluster_df_country$Subtype==subtype,]$Season)){
      if (is.na(season)){next}
      cldf = cluster_df_country[cluster_df_country$Subtype==subtype & cluster_df_country$Season==season,]
      if (nrow(cldf[cldf$Cluster!="UNCLUSTERED",])==0){next}
      prop_df = rbind(prop_df,data.frame(subtype=subtype,season=season,seqprop=0,clprop=0,ncluster=0))
      cldf = cldf[!(is.na(cldf[,1])),]
      cldf = cldf[order(cldf$Seqcount,decreasing=T),]
      cldf = cldf[!(cldf$Cluster=="UNCLUSTERED"),]
      seqprop = cumsum(cldf$Seqcount)/sum(cldf$Seqcount)
      prop_df = rbind(prop_df,data.frame(subtype=subtype,season=season,seqprop=seqprop,clprop=(1:length(seqprop))/length(seqprop),ncluster=(1:length(seqprop))))
    }
  }
  
  prop_df = prop_df[!(is.na(prop_df$subtype)),]
  prop_df$subtype = factor(prop_df$subtype,levels=c("H3N2","H1N1","Yam","Vic"))
  
  which_to_keep = c()
  for (i in 1:nrow(prop_df)){
    subt= prop_df[i,]$subtype
    ssn = prop_df[i,]$season
    prop =epidemic_compositions[epidemic_compositions$Season==ssn & 
                                  epidemic_compositions$Subtype==c("A/H3N2","A/H1N1pdm09","B/Yam","B/Vic")[match(subt,c("H3N2","H1N1","Yam","Vic"))],]$Prop
    if (prop > 0.2){which_to_keep = c(which_to_keep,i)}
  }
  
  prop_df = prop_df[which_to_keep,]
  distribution_plot = ggplot(prop_df[!(is.na(prop_df$subtype)) & prop_df$ncluster>0,],aes(group=season))+ thm +  
    geom_line(aes(x=ncluster,y=seqprop,col=subtype,group=interaction(subtype,season)),lwd=0.4,lty=2)  + 
    xlab("Number of lineages") + ylab("Cumulative proportion\nof sequences") + 
    theme(legend.title=element_blank(),
          legend.key.size = unit(.1, "cm"),
          legend.position = c(0.85,0.25)) + 
    scale_color_brewer(palette='Set1') + 
    scale_x_log10() + annotation_logticks(side='b',outside=T,short=unit(0.7,'mm'),mid=unit(0.7,'mm'),long=unit(0,'mm'),size=0.1) + 
    coord_cartesian(clip='off') 
  
  return(distribution_plot)
}


shannonEntropy <- function(counts.table, log.base = getOption("GeneFamilies.entropy.log.base",
                                                              base::exp(1))) {
  if (length(counts.table) <= 1)
    return(0)
  c.t.s <- sum(counts.table)
  -sum(sapply(counts.table, function(x) ifelse(x>0,x/c.t.s * log(x/c.t.s, base = log.base),0)))#/log(length(counts.table),
}


plotNumberOfStates <- function(){
  
  prop_df_50 = prop_df[prop_df$seqprop<0.5,]
  print(quantile(table(paste0(prop_df_50$season,prop_df_50$subtype))))
  prop_df_50 = prop_df_50[order(prop_df_50$ncluster,decreasing=T),]
  prop_df_50 = prop_df_50[!(duplicated(interaction(prop_df_50$subtype,prop_df_50$season))),]
  prop_df_50$ncluster = prop_df_50$ncluster+1
  prop_df_50$Onset = NA
  for (i in 1:nrow(prop_df_50)){
    prop_df_50[i,]$Onset = getOnsetCountry(prop_df_50[i,]$subtype,prop_df_50[i,]$season)
  }
  
  isolates = read.csv("data/National_PH.csv")
  isolates$Date = NA
  for (year in unique(isolates$YEAR)){
    isolates[isolates$YEAR==year,]$Date = year + isolates[isolates$YEAR==year,]$WEEK/max(isolates[isolates$YEAR==year,]$WEEK)
  }
  isolates$Season = NA
  for (year in unique(isolates$YEAR)){
    isolates[isolates$Date>year+0.5 & isolates$Date<year+1.5,]$Season = year
  }
  isolates = isolates[,c("Season","A..2009.H1N1.","A..H3.", "A..Subtyping.not.Performed.",   "B", "BVic", "BYam")]
  isolates$H1 = isolates$A..2009.H1N1. + (isolates$A..2009.H1N1./(isolates$A..H3.+isolates$A..2009.H1N1.)) * isolates$A..Subtyping.not.Performed.
  isolates$H3 = isolates$A..H3. + (isolates$A..H3./(isolates$A..H3.+isolates$A..2009.H1N1.)) * isolates$A..Subtyping.not.Performed.
  isolates$Vic = isolates$BVic + (isolates$BVic/(isolates$BVic+isolates$BYam)) * isolates$B
  isolates$Yam = isolates$BYam + (isolates$BYam/(isolates$BYam+isolates$BVic)) * isolates$B
  isolates$Prop_H3 = isolates$H3 / (rowSums(isolates[,c("H1","H3","Vic","Yam")]))
  isolates$Prop_H1 = isolates$H1 / (rowSums(isolates[,c("H1","H3","Vic","Yam")]))
  isolates$Prop_Vic = isolates$Vic / (rowSums(isolates[,c("H1","H3","Vic","Yam")]))
  isolates$Prop_Yam = isolates$Yam / (rowSums(isolates[,c("H1","H3","Vic","Yam")]))
  prop_df_50$intens = NA
  prop_df_50$entropy = NA
  for (i in 1:nrow(prop_df_50)){
    if (prop_df_50[i,]$season==2014){next}
    if (prop_df_50[i,]$subtype == "H3N2"){tests <- isolates[isolates$Season==prop_df_50[i,]$season,]$Prop_H3}#A..H3.}
    if (prop_df_50[i,]$subtype == "H1N1"){tests <- isolates[isolates$Season==prop_df_50[i,]$season,]$Prop_H1}#A..2009.H1N1.}
    if (prop_df_50[i,]$subtype == "Vic"){tests <- isolates[isolates$Season==prop_df_50[i,]$season,]$Prop_Vic}#$BVic}
    if (prop_df_50[i,]$subtype == "Yam"){tests <- isolates[isolates$Season==prop_df_50[i,]$season,]$Prop_Yam}#BYam}
    # tests[is.na(tests)] = 0
    tests = tests/sum(tests,na.rm=T)
    intens = 1/(-sum(tests*log(tests),na.rm=T))
    prop_df_50[i,]$intens = intens
    mt = mtdt[mtdt$Season==prop_df_50[i,]$season & mtdt$Subtype==prop_df_50[i,]$subtype,]$Cluster
    prop_df_50[i,]$entropy = shannonEntropy(table(mt))/log(sum(table(mt)))
  }
  prop_df_50$intens = (prop_df_50$intens - min(prop_df_50$intens,na.rm=T)) / (max(prop_df_50$intens,na.rm=T) - min(prop_df_50$intens,na.rm=T))
  cor.test(prop_df_50$intens,prop_df_50$entropy)
  cor.test(prop_df_50$Onset,prop_df_50$entropy)
  
  cluster_df$Circulating = cluster_df$Percentage>0.05
  number_of_states = aggregate(Circulating~Cluster+Season+Subtype,FUN=sum,data=cluster_df[cluster_df$Cluster!="UNCLUSTERED",])
  number_of_states = data.frame(N = rep(1:50,4),Subtype = c(rep("H3N2",50),rep("H1N1",50),rep("Yam",50),rep("Vic",50)),ncluster = 
                                  c(sapply(c("A/H3N2","A/H1N1pdm09","B/Yam","B/Vic"),function(y)sapply(1:50,function(x)nrow(number_of_states[number_of_states$Subtype==y & number_of_states$Circulating>=x,])))))
  number_of_states$Subtype = factor(number_of_states$Subtype,levels=c("H3N2","H1N1","Yam","Vic"))
  
  number_of_states_hist = ggplot(number_of_states) +thm + 
    theme(strip.background=element_rect(colour="white",fill="white"),
          strip.text.x = element_text(size = 5),
          legend.position = c(0.8,0.8),
          legend.key.size = unit(.1, "cm"),
          legend.title=element_blank()) + 
    geom_line(aes(x=N,color=Subtype,fill=Subtype,y=ncluster),lty=2,lwd=0.4) +  
    theme(legend.position="NA") + xlab("Number of states") + ylab("Lineages in at least\nnumber of states ") +  
    scale_fill_brewer(palette='Set1') + 
    scale_y_log10() + coord_cartesian(clip='off') +
    annotation_logticks(side='l',outside=T,short=unit(0.7,'mm'),mid=unit(0.7,'mm'),long=unit(0,'mm'),size=0.1) + 
    scale_color_brewer(palette='Set1')
  
  length(unique(cluster_df[cluster_df$Circulating==1 & cluster_df$Cluster != "UNCLUSTERED",]$Cluster))
  sum(number_of_states[number_of_states$N==10,]$ncluster)
  sum(number_of_states[number_of_states$N==25,]$ncluster)
  return(number_of_states_hist)
}


plotOnset <- function(){
  
  onset_plot = ggplot(cluster_df_country[cluster_df_country$Cluster!="UNCLUSTERED" & 
                                           !(is.infinite(cluster_df_country$Time_Since_Onset)) &
                                           cluster_df_country$Substantial_Subtype==1,],
                      aes(x=Time_Since_Onset,as.numeric(Percentage))) + thm + 
    geom_point(aes(bg=Subtype),position=position_jitter(width=0.005),pch=21,col='black',stroke=0.1,cex=1) + 
    xlab("Relative timing of\nestablishment (y)") + scale_fill_brewer(palette='Set1') + scale_color_brewer(palette='Set1') + 
    ylab("Lineage relative\nsize (%)") + scale_y_log10(breaks=c(0.001,0.01,0.1,1,10,100), labels=trans_format("log10", math_format(10^.x))) +
    stat_smooth(geom = "ribbon", aes(group=Subtype,fill=Subtype),method='lm',lty=1,lwd=0.15,span=5,color='black',alpha=1,level=0.5)+
    annotation_logticks(side='l',outside=T,short=unit(0.7,'mm'),mid=unit(0.7,'mm'),long=unit(0,'mm'),size=0.1) + coord_cartesian(clip='off') +
    theme(legend.position = c(0.85,0.85),
          legend.key.size = unit(.1, "cm"),
          legend.title = element_blank())
  
  for (subtype in unique(cluster_df_country$Subtype)){
    m_onset = cor.test(cluster_df_country[cluster_df_country$Cluster!="UNCLUSTERED" & 
                                            cluster_df_country$Subtype==subtype & 
                                            !(is.infinite(cluster_df_country$Time_Since_Onset)) &
                                            cluster_df_country$Substantial_Subtype==1,]$Time_Since_Onset,
                       log(as.numeric(cluster_df_country[cluster_df_country$Cluster!="UNCLUSTERED" & 
                                                           !(is.infinite(cluster_df_country$Time_Since_Onset)) &
                                                           cluster_df_country$Subtype==subtype & 
                                                           cluster_df_country$Substantial_Subtype==1,]$Percentage),10),method='spearman')
    print(m_onset$estimate)
    print(m_onset$p.value)
  }
  return(onset_plot)
}

plotSampling <- function(){
  sampling_plot = ggplot(cluster_df_country[cluster_df_country$Cluster!="UNCLUSTERED" & !(is.infinite(cluster_df_country$Time_Since_Sample)) & cluster_df_country$Substantial_Subtype==1,],aes(x=Time_Since_Sample,as.numeric(Percentage))) + 
    thm + geom_point(aes(bg=Subtype),pch=21,col='black',stroke=0.1,cex=1) +
    xlab("Relative timing of\nfirst sampling (y)") + scale_fill_brewer(palette='Set1',labels=c("A/H3N2","A/H1N1pdm09","B/Yamagata","B/Victoria")) +
    scale_color_brewer(palette='Set1',labels=c("A/H3N2","A/H1N1pdm09","B/Yamagata","B/Victoria")) + 
    ylab("Lineage relative\nsize (%)") + scale_y_log10(breaks=c(0.001,0.01,0.1,1,10,100), labels=trans_format("log10", math_format(10^.x))) +
    stat_smooth(geom = "ribbon", aes(group=Subtype,fill=Subtype),method='lm',lty=1,lwd=0.15,span=5,color='black',alpha=1,level=0.5)+
    annotation_logticks(side='l',outside=T,short=unit(0.7,'mm'),mid=unit(0.7,'mm'),long=unit(0,'mm'),size=0.1) + coord_cartesian(clip='off') +
    theme(legend.position = c(0.85,0.85),
          legend.key.size = unit(.1, "cm"),
          legend.title = element_blank()) +
    scale_x_continuous(breaks=c(-0.5,0,0.5))
  
  m_sample = cor.test(cluster_df_country[cluster_df_country$Cluster!="UNCLUSTERED" & 
                                           !(is.infinite(cluster_df_country$Time_Since_Sample)) &
                                           cluster_df_country$Substantial_Subtype==1,]$Time_Since_Sample,
                      log(as.numeric(cluster_df_country[cluster_df_country$Cluster!="UNCLUSTERED" & 
                                                          !(is.infinite(cluster_df_country$Time_Since_Sample)) &
                                                          cluster_df_country$Substantial_Subtype==1,]$Percentage),10),method='spearman')
  return(sampling_plot)
}


getMixedModel <- function(){
  m_onset_h3 = lmer(log(Percentage,10) ~ Time_Since_Onset +  (1 | Season) ,data=cluster_df_country[cluster_df_country$Cluster!="UNCLUSTERED" & 
                                                                                                     !(is.infinite(cluster_df_country$Time_Since_Onset)) &
                                                                                                     cluster_df_country$Substantial_Subtype==1&
                                                                                                     cluster_df_country$Subtype=="H3N2",])
  m_onset_h1 = lmer(log(Percentage,10) ~ Time_Since_Onset +  (1 | Season) ,data=cluster_df_country[cluster_df_country$Cluster!="UNCLUSTERED" & 
                                                                                                     !(is.infinite(cluster_df_country$Time_Since_Onset)) &
                                                                                                     cluster_df_country$Substantial_Subtype==1 &
                                                                                                     cluster_df_country$Subtype=="H1N1",])
  m_onset_yam = lmer(log(Percentage,10) ~ Time_Since_Onset +  (1 | Season) ,data=cluster_df_country[cluster_df_country$Cluster!="UNCLUSTERED" & 
                                                                                                      !(is.infinite(cluster_df_country$Time_Since_Onset)) &
                                                                                                      cluster_df_country$Substantial_Subtype==1 &
                                                                                                      cluster_df_country$Subtype=="Yam",])
  m_onset_vic = lm(log(Percentage,10) ~ Time_Since_Onset ,data=cluster_df_country[cluster_df_country$Cluster!="UNCLUSTERED" & 
                                                                                    !(is.infinite(cluster_df_country$Time_Since_Onset)) &
                                                                                    cluster_df_country$Substantial_Subtype==1 &
                                                                                    cluster_df_country$Subtype=="Vic",])
  
  pct_changes = (10^(c(fixef(m_onset_h3)[2],fixef(m_onset_h1)[2],fixef(m_onset_yam)[2],coefficients(m_onset_vic)[2])/52) - 1)*100
}


plotTimingByThreshold <- function(){
  
  inc_thresholds = c(0.001,0.005,0.01,0.02,0.05,0.1)
  alt_thresh_cdc = lapply(inc_thresholds,function(x)getClusterData_subsamp(clusterpath, x)[[3]])
  all_country_dfs = list()
  for (i in 1:length(alt_thresh_cdc)){
    df = alt_thresh_cdc[[i]]
    df$thresh = inc_thresholds[i]
    df$Time_Since_Onset = as.numeric(df$RelativeOnset_Date)-df$EpidemicOnset
    df$Time_Since_Sample = df$FirstSamp-df$EpidemicOnset
    df$Time_Since_TMRCA = df$TMRCA-df$EpidemicOnset
    all_country_dfs[[i]] = df
  }
  all_country_data = do.call(rbind,all_country_dfs)
  
  table((all_country_data %>% filter(Cluster!="UNCLUSTERED" &!(is.infinite(Time_Since_Onset)) & Substantial_Subtype==1))$thresh)
  
  all_country_data$thresh = factor(all_country_data$thresh,labels=c("0.1% (N = 942)",
                                                                    "0.5% (N = 884)",
                                                                    "1% (N=840)",
                                                                    "2% (N=785)",
                                                                    "5% (N=569)",
                                                                    "10% (N=323)"))
  
  timing_plot_per_threshold = ggplot(all_country_data %>% filter(Cluster!="UNCLUSTERED" &!(is.infinite(Time_Since_Onset)) & Substantial_Subtype==1),
                                     aes(x=Time_Since_Onset,as.numeric(Percentage))) + thm + theme(panel.grid.major = element_line(color = "grey90",size = 0.1,linetype = 1)) +  
    geom_point(aes(bg=Subtype),position=position_jitter(width=0.005),pch=21,col='black',stroke=0.1,cex=1) + 
    xlab("Relative timing of\nestablishment (y)") + scale_fill_brewer(palette='Set1') + scale_color_brewer(palette='Set1') + 
    ylab("Lineage relative\nsize (%)") + scale_y_log10(breaks=c(0.001,0.01,0.1,1,10,100), labels=trans_format("log10", math_format(10^.x))) +
    stat_smooth(aes(group=Subtype),col='white',method='lm',se=F,lty=1,lwd=0.4) + stat_smooth(aes(group=Subtype,col=Subtype),method='lm',se=F,lty=1,lwd=0.3) +
    annotation_logticks(side='l',outside=T,short=unit(0.7,'mm'),mid=unit(0.7,'mm'),long=unit(0,'mm'),size=0.1) + coord_cartesian(clip='off') +
    theme(legend.position = c(0.4,0.14)) + theme(legend.key.size = unit(.1, "cm"),legend.title = element_blank()) + facet_wrap(thresh~.)
  ggarrange(legplot,timing_plot_per_threshold + theme(axis.text = element_text(color="black",size=8),
                                                      strip.text.x = element_text(size = 8),
                                                      legend.text=element_text(size=8),
                                                      strip.background = element_blank(),
                                                      legend.position='none',
                                                      axis.title = element_text(size=8)),heights=c(0.1,0.9),ncol=1)
  
  
  ggsave("SuppFig_Timing.pdf",width=2080,units="px",height=1500,dpi=320)
}


table(table(cluster_df[cluster_df$Circulating==T & cluster_df$Cluster!="UNCLUSTERED",]$Cluster)>10)


mrg = theme(plot.margin = margin(c(6,0,4,0),'cm'),
            plot.title=element_text(size=7,margin=margin(-3,0,0,0)),
            axis.ticks.x = element_line(size=0.2),
            axis.line.x = element_line(size=0.1))
tree_1 = plotTree(clusterpath,"H1N1",2019,mtdt,c(),F) + mrg 
tree_1 = addIsolationsToTree(tree_1,2019,"H1N1")
tree_2 = plotTree(clusterpath,"H3N2",2016,mtdt,c(),F) + mrg
tree_2 = addIsolationsToTree(tree_2,2016,"H3N2")
tree_3 = plotTree(clusterpath,"H3N2",2017,mtdt,c(),F) + mrg
tree_3 = addIsolationsToTree(tree_3,2017,"H3N2")
tree_4 = plotTree(clusterpath,"H1N1",2015,mtdt,c(),F) + mrg
tree_4 = addIsolationsToTree(tree_4,2015,"H1N1")
tree_5 = plotTree(clusterpath,"Vic",2016,mtdt,c(),F) + mrg
tree_5 = addIsolationsToTree(tree_5,2016,"Vic")
tree_6 = plotTree(clusterpath,"Yam",2017,mtdt,c(),F) + mrg
tree_6 = addIsolationsToTree(tree_6,2017,"Yam")

leg = get_legend(sampling_plot+theme(legend.position = 'top',legend.spacing.x = unit(.15,"cm")))
legplot =as_ggplot(leg)

distribution_plot = plotDistribution()
number_of_states_plot = plotNumberOfStates()
sampling_plot = plotSampling()
onset_plot = plotOnset()

ggarrange(ggarrange(tree_1,tree_2,ncol=1),ggarrange(tree_3,tree_4,ncol=1),ggarrange(tree_5,tree_6,ncol=1),
          ggarrange(ggarrange(distribution_plot+theme(legend.position='none'),number_of_states_plot+theme(legend.position='none'),nrow=1),
                    ggarrange(sampling_plot+theme(legend.position='none'),onset_plot+theme(legend.position='none'),nrow=1),legplot,nrow=3,heights=c(0.5,0.5,0.05)),ncol=4,widths=c(0.23,0.23,0.23,0.36+0.36))

ggsave("Figure_1_resubmit.pdf",width=2080,units="px",height=1000,dpi=320)

source("SuppFig_PhyParams.R")
plotTimingByThreshold()
