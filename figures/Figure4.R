
plotTrueAndSimulated_2019_vic <- function(fctr){
  
  leg = as_ggplot(get_legend(plotCluster_true(3,0,0,c()) + theme(legend.position='top',legend.text=element_text(size=5),legend.title=element_text(size=5),legend.key.height=unit(4,'pt'),legend.key.width=unit(15,'pt')) +
                               guides(fill = guide_colourbar(title.position="top", title.hjust = 0.5))))
  
  season = 2019
  subtype = "Vic"
  cdc_2019 = cluster_df_country[cluster_df_country$Season==2019 & cluster_df_country$Subtype=="Vic",]
  cdc_2019 = cdc_2019[order(cdc_2019$Seqcount,decreasing=T),]
  clusters = cdc_2019$Cluster[1:7]
  clusters = clusters[!(clusters=="UNCLUSTERED")]
  onsetstates = c("California","Florida","Texas","Louisiana","Nevada","Washington")
  #indexonset = unlist(lapply(clusters,function(x)min(cluster_df[cluster_df$Cluster==x,]$Onset,na.rm=T)))
  indexstates = match(onsetstates,states)-1
  
  mrg = theme(plot.margin = margin(c(-1,-3,-3,-4),'cm'))
  
  indexonset = c(21, 22, 23, 21, 22, 25)
  
  sim = simulate_fun_separate(nday,indexstates,indexonset*7,r0,d_commuting,d_airtravel,pops_raw,1)

  pcts = getPcts(sim)
  
  sim_plots = list()
  for (i in 1:length(clusters)){
    sim_plots[[i]] = plotCluster_sim(pcts[[1]][,i],pcts[[2]][,i],indexonset[i]-min(indexonset),fctr,clusters)+mrg
  }
  
  sims = ggarrange(plotlist=sim_plots,nrow=1)
  
  true_plots = list()
  for (i in 1:length(clusters)){
    true_plots[[i]] = plotCluster_true(clusters[i],indexonset[i]-min(indexonset),fctr,clusters)+mrg
  }
  
  true_ann = ggplot()+theme_void() #+ annotate("text",x=0,y=0,label='Data',size=5/.pt) + xlim(c(-0.5,0.5)) + ylim(c(-0.5,0.5))
  sim_ann =ggplot()+theme_void() #+ annotate("text",x=0,y=0,label='Simulated',size=5/.pt) + xlim(c(-0.5,0.5)) + ylim(c(-0.5,0.5))
  
  titleplots = list()
  for (i in 1:length(clusters)){
    titleplots[[i]] = ggplot()+theme_void() + annotate("text",x=0,y=0,label=i,size=5/.pt) + xlim(c(-0.5,0.5)) + ylim(c(-0.5,0.5))
  }
  titlerow = ggarrange(ggplot()+theme_void(),ggarrange(plotlist=titleplots,nrow=1),widths=c(0.07,0.95))
  
  trues = ggarrange(plotlist=true_plots,nrow=1)
  treeplot = plotTree_vert(clusterpath,"Vic",2019,mtdt,clusters,c(0.5,2),F)[[1]]
  ggarrange(treeplot,titlerow,ggarrange(true_ann/sim_ann,ggarrange(trues,sims,ncol=1),nrow=1,widths=c(0.07,0.95)),ggarrange(ggplot()+theme_void(),leg,ggplot()+theme_void(),nrow=1,widths=c(0.6,0.3,0.4)),nrow=4,heights=c(0.45,0.05,0.7,0.15))
  
  ggsave("Figure_4.pdf",width=180,height=75,units='mm')  #sims = (plotCluster_sim(pcts[[1]][,2],pcts[[2]][,2],0,20)|plotCluster_sim(pcts[[1]][,1],pcts[[2]][,1],2,20)) + plot_annotation(title='True') &theme(text=element_text(size=6,hjust=0.5))
  return(plt)
}



plotTrueAndSimulated_2019_vic_comm <- function(fctr){
  
  leg = as_ggplot(get_legend(plotCluster_true(3,0,0,c()) + theme(legend.position='top',legend.text=element_text(size=5),legend.title=element_text(size=5),legend.key.height=unit(4,'pt'),legend.key.width=unit(15,'pt')) +
                               guides(fill = guide_colourbar(title.position="top", title.hjust = 0.5))))
  
  season = 2019
  subtype = "Vic"
  cdc_2019 = cluster_df_country[cluster_df_country$Season==2019 & cluster_df_country$Subtype=="Vic",]
  cdc_2019 = cdc_2019[order(cdc_2019$Seqcount,decreasing=T),]
  clusters = cdc_2019$Cluster[1:7]
  clusters = clusters[!(clusters=="UNCLUSTERED")]
  onsetstates = c("California","Florida","Texas","Louisiana","Nevada","Washington")
  indexonset = unlist(lapply(clusters,function(x)min(cluster_df[cluster_df$Cluster==x,]$Onset,na.rm=T)))
  indexstates = match(onsetstates,states)-1
  
  mrg = theme(plot.margin = margin(c(-1,-3,-3,-4),'cm'))
  
  
  indexonset = c(21, 22, 23, 21, 22, 25)
  
  sim = simulate_fun_separate(nday,indexstates,indexonset*7,r0,mob[[4]],d_airtravel*0,pops_raw,1)
  
  pcts = getPcts(sim)
  
  sim_plots = list()
  for (i in 1:length(clusters)){
    sim_plots[[i]] = plotCluster_sim(pcts[[1]][,i],pcts[[2]][,i],indexonset[i]-min(indexonset),fctr,clusters)+mrg
  }
  
  sims = ggarrange(plotlist=sim_plots,nrow=1)
  
  true_plots = list()
  for (i in 1:length(clusters)){
    true_plots[[i]] = plotCluster_true(clusters[i],indexonset[i]-min(indexonset),fctr,clusters)+mrg
  }
  
  true_ann = ggplot()+theme_void() #+ annotate("text",x=0,y=0,label='Data',size=5/.pt) + xlim(c(-0.5,0.5)) + ylim(c(-0.5,0.5))
  sim_ann =ggplot()+theme_void() #+ annotate("text",x=0,y=0,label='Simulated',size=5/.pt) + xlim(c(-0.5,0.5)) + ylim(c(-0.5,0.5))
  
  titleplots = list()
  for (i in 1:length(clusters)){
    titleplots[[i]] = ggplot()+theme_void() + annotate("text",x=0,y=0,label=i,size=5/.pt) + xlim(c(-0.5,0.5)) + ylim(c(-0.5,0.5))
  }
  titlerow = ggarrange(ggplot()+theme_void(),ggarrange(plotlist=titleplots,nrow=1),widths=c(0.07,0.95))
  
  trues = ggarrange(plotlist=true_plots,nrow=1)
  treeplot = plotTree_vert(clusterpath,"Vic",2019,mtdt,clusters,c(0.5,2),F)[[1]]
  ggarrange(true_ann/sim_ann,ggarrange(trues,sims,ncol=1),nrow=1,widths=c(0.07,0.95))
  
  ggsave("Figure_4_comm.pdf",width=180,height=50,units='mm')  #sims = (plotCluster_sim(pcts[[1]][,2],pcts[[2]][,2],0,20)|plotCluster_sim(pcts[[1]][,1],pcts[[2]][,1],2,20)) + plot_annotation(title='True') &theme(text=element_text(size=6,hjust=0.5))
  return(plt)
}


plotTrueAndSimulated_2019_vic_air <- function(fctr){
  
  leg = as_ggplot(get_legend(plotCluster_true(3,0,0,c()) + theme(legend.position='top',legend.text=element_text(size=5),legend.title=element_text(size=5),legend.key.height=unit(4,'pt'),legend.key.width=unit(15,'pt')) +
                               guides(fill = guide_colourbar(title.position="top", title.hjust = 0.5))))
  
  season = 2019
  subtype = "Vic"
  cdc_2019 = cluster_df_country[cluster_df_country$Season==2019 & cluster_df_country$Subtype=="Vic",]
  cdc_2019 = cdc_2019[order(cdc_2019$Seqcount,decreasing=T),]
  clusters = cdc_2019$Cluster[1:7]
  clusters = clusters[!(clusters=="UNCLUSTERED")]
  onsetstates = c("California","Florida","Texas","Louisiana","Nevada","Washington")
  indexonset = unlist(lapply(clusters,function(x)min(cluster_df[cluster_df$Cluster==x,]$Onset,na.rm=T)))
  indexstates = match(onsetstates,states)-1
  
  mrg = theme(plot.margin = margin(c(-1,-3,-3,-4),'cm'))
  
  
  indexonset = c(21, 22, 23, 21, 22, 25)
  
  sim = simulate_fun_separate(nday,indexstates,indexonset*7,r0,mob[[3]],d_airtravel*0,pops_raw,1)
  
  pcts = getPcts(sim)
  
  sim_plots = list()
  for (i in 1:length(clusters)){
    sim_plots[[i]] = plotCluster_sim(pcts[[1]][,i],pcts[[2]][,i],indexonset[i]-min(indexonset),fctr,clusters)+mrg
  }
  
  sims = ggarrange(plotlist=sim_plots,nrow=1)
  
  true_plots = list()
  for (i in 1:length(clusters)){
    true_plots[[i]] = plotCluster_true(clusters[i],indexonset[i]-min(indexonset),fctr,clusters)+mrg
  }
  
  true_ann = ggplot()+theme_void() #+ annotate("text",x=0,y=0,label='Data',size=5/.pt) + xlim(c(-0.5,0.5)) + ylim(c(-0.5,0.5))
  sim_ann =ggplot()+theme_void() #+ annotate("text",x=0,y=0,label='Simulated',size=5/.pt) + xlim(c(-0.5,0.5)) + ylim(c(-0.5,0.5))
  
  titleplots = list()
  for (i in 1:length(clusters)){
    titleplots[[i]] = ggplot()+theme_void() + annotate("text",x=0,y=0,label=i,size=5/.pt) + xlim(c(-0.5,0.5)) + ylim(c(-0.5,0.5))
  }
  titlerow = ggarrange(ggplot()+theme_void(),ggarrange(plotlist=titleplots,nrow=1),widths=c(0.07,0.95))
  
  trues = ggarrange(plotlist=true_plots,nrow=1)
  treeplot = plotTree_vert(clusterpath,"Vic",2019,mtdt,clusters,c(0.5,2),F)[[1]]
  ggarrange(true_ann/sim_ann,ggarrange(trues,sims,ncol=1),nrow=1,widths=c(0.07,0.95))
  
  ggsave("Figure_4_air.pdf",width=180,height=50,units='mm')  #sims = (plotCluster_sim(pcts[[1]][,2],pcts[[2]][,2],0,20)|plotCluster_sim(pcts[[1]][,1],pcts[[2]][,1],2,20)) + plot_annotation(title='True') &theme(text=element_text(size=6,hjust=0.5))
  return(plt)
}

plotTrueAndSimulated_2019_vic(15)
plotTrueAndSimulated_2019_vic_comm(15)
plotTrueAndSimulated_2019_vic_air(15)


params = read.csv("/Users/simondejong/US_phylo/2019_Vic.log",sep='\t')
date_decimal(mean(params$cluster_965.age.root[100:1000]))
date_decimal(mean(params$cluster_967.age.root[100:1000]))

date_decimal(quantile(params$cluster_965.age.root[100:1000],c(0.025,0.5,0.975)))
date_decimal(quantile(params$cluster_967.age.root[100:1000],c(0.025,0.5,0.975)))
date_decimal(min(mtdt[mtdt$Cluster==967,]$Date))


