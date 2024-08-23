source("plotFunctions.R")

### Get proportion of sequences that is accounted for by a proportion of clusters

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
  prop =epidemic_compositions[epidemic_compositions$Season==ssn & epidemic_compositions$Subtype==subt,]$Prop
  if (prop > 0.1){which_to_keep = c(which_to_keep,i)}
}

prop_df = prop_df[which_to_keep,]
distribution_plot = ggplot(prop_df[!(is.na(prop_df$subtype)) & prop_df$ncluster>0,],aes(group=season))+thm+ theme(panel.grid.major = element_line(color = "grey90",size = 0.1,linetype = 1)) +  
  geom_line(aes(x=ncluster,y=seqprop,col=subtype,group=interaction(subtype,season)),lwd=0.2,lty=2)  + 
  xlab("Number of lineages") + ylab("Cumulative proportion\nof sequences") + theme(legend.position='top') + 
  scale_color_brewer(palette='Set1') + theme(legend.title=element_blank()) + theme(legend.position = c(0.85,0.25)) + theme(legend.key.size = unit(.1, "cm")) +
  scale_x_log10() + annotation_logticks(side='b',outside=T,short=unit(0.7,'mm'),mid=unit(0.7,'mm'),long=unit(0,'mm'),size=0.1) + coord_cartesian(clip='off') 


prop_df_50 = prop_df[prop_df$seqprop<0.5,]
print(quantile(table(paste0(prop_df_50$season,prop_df_50$subtype))))


# Get the distibution of the number of states that see substantial circulation of a transmission lineage

cluster_df$Circulating = cluster_df$Percentage>0.05
number_of_states = aggregate(Circulating~Cluster+Season+Subtype,FUN=sum,data=cluster_df[cluster_df$Cluster!="UNCLUSTERED",])
number_of_states = data.frame(N = rep(1:50,4),Subtype = c(rep("H3N2",50),rep("H1N1",50),rep("Yam",50),rep("Vic",50)),ncluster = 
                                c(sapply(c("A/H3N2","A/H1N1pdm09","B/Yam","B/Vic"),function(y)sapply(1:50,function(x)nrow(number_of_states[number_of_states$Subtype==y & number_of_states$Circulating>=x,])))))
number_of_states$Subtype = factor(number_of_states$Subtype,levels=c("H3N2","H1N1","Yam","Vic"))

number_of_states_hist = ggplot(number_of_states) +thm + 
  theme(panel.grid.major = element_line(color = "grey90",size = 0.1,linetype = 1)) + 
  geom_line(aes(x=N,color=Subtype,fill=Subtype,y=ncluster),lty=2,lwd=0.2) +  
  theme(legend.position="NA") + xlab("Number of states") + ylab("Lineages in at least\nnumber of states ") +  
  scale_fill_brewer(palette='Set1') + 
  theme(strip.background=element_rect(colour="white",fill="white"),
        strip.text.x = element_text(size = 5)) + #scale_y_continuous(breaks=c(0,10,100,1000),trans=scales::pseudo_log_trans(base = 10)) +
  coord_cartesian(clip='off') + scale_y_log10() + annotation_logticks(side='l',outside=T,short=unit(0.7,'mm'),mid=unit(0.7,'mm'),long=unit(0,'mm'),size=0.1) + coord_cartesian(clip='off') + 
  scale_color_brewer(palette='Set1')+ theme(legend.title=element_blank()) + theme(legend.position = c(0.8,0.8)) + theme(legend.key.size = unit(.1, "cm"))

length(unique(cluster_df[cluster_df$Circulating==1 & cluster_df$Cluster != "UNCLUSTERED",]$Cluster))
sum(number_of_states[number_of_states$N==10,]$ncluster)
sum(number_of_states[number_of_states$N==25,]$ncluster)


# Get correlation between emergence and establishment timing and lineage size

cluster_df_country$Time_Since_Onset = as.numeric(cluster_df_country$RelativeOnset_Date)-cluster_df_country$EpidemicOnset
cluster_df_country$Time_Since_Sample = cluster_df_country$FirstSamp-cluster_df_country$EpidemicOnset
cluster_df_country$Time_Since_TMRCA = cluster_df_country$TMRCA-cluster_df_country$EpidemicOnset

onset_plot = ggplot(cluster_df_country[cluster_df_country$Cluster!="UNCLUSTERED" & 
                                         !(is.infinite(cluster_df_country$Time_Since_Onset)) &
                                         cluster_df_country$Substantial_Subtype==1,],
                    aes(x=Time_Since_Onset,as.numeric(Percentage))) + thm + theme(panel.grid.major = element_line(color = "grey90",size = 0.1,linetype = 1)) +  
  geom_point(aes(bg=Subtype),position=position_jitter(width=0.005),pch=21,col='black',stroke=0.1,cex=0.5) + 
  xlab("Relative timing of\nestablishment (y)") + scale_fill_brewer(palette='Set1') + scale_color_brewer(palette='Set1') + 
  ylab("Lineage relative\nsize (%)") + scale_y_log10(breaks=c(0.001,0.01,0.1,1,10,100), labels=trans_format("log10", math_format(10^.x))) +
  stat_smooth(aes(group=Subtype),col='white',method='lm',se=F,lty=1,lwd=0.4) + stat_smooth(aes(group=Subtype,col=Subtype),method='lm',se=F,lty=1,lwd=0.3) +
  annotation_logticks(side='l',outside=T,short=unit(0.7,'mm'),mid=unit(0.7,'mm'),long=unit(0,'mm'),size=0.1) + coord_cartesian(clip='off') +
  theme(legend.position = c(0.85,0.85)) + theme(legend.key.size = unit(.1, "cm"),legend.title = element_blank())

sampling_plot = ggplot(cluster_df_country[cluster_df_country$Cluster!="UNCLUSTERED" & !(is.infinite(cluster_df_country$Time_Since_Sample)) & cluster_df_country$Substantial_Subtype==1,],aes(x=Time_Since_Sample,as.numeric(Percentage))) + 
  thm + theme(panel.grid.major = element_line(color = "grey90",size = 0.1,linetype = 1)) +  geom_point(aes(bg=Subtype),pch=21,col='black',stroke=0.1,cex=0.5) +
  xlab("Relative timing of\nfirst sampling (y)") + scale_fill_brewer(palette='Set1') + scale_color_brewer(palette='Set1') + 
  ylab("Lineage relative\nsize (%)") + scale_y_log10(breaks=c(0.001,0.01,0.1,1,10,100), labels=trans_format("log10", math_format(10^.x))) +
  stat_smooth(aes(group=Subtype),col='white',method='lm',se=F,lty=1,lwd=0.4) + stat_smooth(aes(group=Subtype,col=Subtype),method='lm',se=F,lty=1,lwd=0.3) +
  annotation_logticks(side='l',outside=T,short=unit(0.7,'mm'),mid=unit(0.7,'mm'),long=unit(0,'mm'),size=0.1) + coord_cartesian(clip='off') +
  theme(legend.position = c(0.85,0.85)) + theme(legend.key.size = unit(.1, "cm"),legend.title = element_blank()) +
  scale_x_continuous(breaks=c(-0.5,0,0.5))

 
m_onset = cor.test(cluster_df_country[cluster_df_country$Cluster!="UNCLUSTERED" & 
                                        !(is.infinite(cluster_df_country$Time_Since_Onset)) &
                                        cluster_df_country$Substantial_Subtype==1,]$Time_Since_Onset,
                   log(as.numeric(cluster_df_country[cluster_df_country$Cluster!="UNCLUSTERED" & 
                                                       !(is.infinite(cluster_df_country$Time_Since_Onset)) &
                                                       cluster_df_country$Substantial_Subtype==1,]$Percentage),10),method='spearman')
m_onset$estimate
m_onset$p.value

m_firstsample = cor.test(cluster_df_country[cluster_df_country$Cluster!="UNCLUSTERED" & 
                                              !(is.infinite(cluster_df_country$Time_Since_Sample)) &
                                              cluster_df_country$Substantial_Subtype==1,]$Time_Since_Sample,
                         log(as.numeric(cluster_df_country[cluster_df_country$Cluster!="UNCLUSTERED" & 
                                                             !(is.infinite(cluster_df_country$Time_Since_Sample)) &
                                                             cluster_df_country$Substantial_Subtype==1,]$Percentage),10))
m_firstsample$estimate
m_firstsample$p.value


table(table(cluster_df[cluster_df$Circulating==T & cluster_df$Cluster!="UNCLUSTERED",]$Cluster)>10)



# Get trees to plot
mrg = theme(plot.margin = margin(c(4,0,4,0),'cm'),plot.title=element_text(size=5,margin=margin(-3,0,0,0)))#,axis.text.x = element_blank(),axis.line.x = element_blank(),axis.ticks.x = element_blank())
tree_1 = plotTree(clusterpath,"H1N1",2019,mtdt,c(),F) + mrg #+ theme(axis.text.x = element_blank(),axis.line.x = element_blank(),axis.ticks.x = element_blank())
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


ggarrange(ggarrange(tree_1,tree_2,ncol=1),ggarrange(tree_3,tree_4,ncol=1),ggarrange(tree_5,tree_6,ncol=1),ggarrange(distribution_plot,sampling_plot,ncol=1,heights=c(0.5,0.5)),ggarrange(number_of_states_hist,onset_plot,ncol=1),nrow=1,widths=c(0.23,0.23,0.23,0.38,0.38))
ggsave("Figure_1.pdf",width=125,units="mm",height=55)

tmrca_plot
ggsave("SuppFig_TMRCA.pdf",width=50,units="mm",height=60)

