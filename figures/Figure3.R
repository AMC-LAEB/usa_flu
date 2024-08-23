library(ggpubr)
library(usmap)
source("plotFunctions.R")
source("simulateMobility.R")

plotTrueAndSimulated_2018 <- function(fctr){
  
  clusters = c(670,671)
  states = c(state.name[-c(2,11)],"District of Columbia")
  
  indexstates = match(c("Georgia","Nebraska"),states)-1
  
  indexonset = cluster_df_country$RelativeOnset_Country[match(clusters,cluster_df_country$Cluster)]
  indexonset = as.numeric(indexonset) - min(as.numeric(indexonset))
    
  sim_onlycommuting = simulate_fun_separate(nday,indexstates,indexonset*7,r0,mob[[4]],d_airtravel*0,pops_raw,1)
  sim_onlyairtravel = simulate_fun_separate(nday,indexstates,indexonset*7,r0,mob[[3]],d_airtravel*0,pops_raw,1)
  
  pcts = getPcts(sim)
  pcts_commuting = getPcts(sim_onlycommuting)
  pcts_airtravel = getPcts(sim_onlyairtravel)
  
  mrg = theme(plot.margin = margin(c(4,-3,-3,-4),'cm'))
  
  sims = ggarrange(plotCluster_sim(pcts[[1]][,1],pcts[[2]][,1],indexonset[1],fctr,clusters)+mrg,
                   plotCluster_sim(pcts[[1]][,2],pcts[[2]][,2],indexonset[2],fctr,clusters)+mrg,ncol=2)
  
  trues = ggarrange(plotCluster_true(clusters[1],indexonset[1],fctr,clusters)+mrg,
                    plotCluster_true(clusters[2],indexonset[2],fctr,clusters)+mrg,ncol=2)
  
  sims_commuting = ggarrange(plotCluster_sim(pcts_commuting[[1]][,1],pcts_commuting[[2]][,1],indexonset[1],fctr,clusters)+mrg,
                             plotCluster_sim(pcts_commuting[[1]][,2],pcts_commuting[[2]][,2],indexonset[2],fctr,clusters)+mrg,ncol=2)
  
  sims_airtravel = ggarrange(plotCluster_sim(pcts_airtravel[[1]][,1],pcts_airtravel[[2]][,1],indexonset[1],fctr,clusters)+mrg,
                             plotCluster_sim(pcts_airtravel[[1]][,2],pcts_airtravel[[2]][,2],indexonset[2],fctr,clusters)+mrg,ncol=2)
  
  plt = ggarrange(trues,sims_commuting,sims_airtravel,ncol=1)
  return(plt)
}



plotTrueAndSimulated_2017 <- function(fctr){
  
  clusters = c(96,104)
  states = c(state.name[-c(2,11)],"District of Columbia")
  indexstates = match(c("Mississippi","Oregon"),states)-1
  indexonset = cluster_df_country$RelativeOnset_Country[match(clusters,cluster_df_country$Cluster)]
  indexonset = as.numeric(indexonset) - min(as.numeric(indexonset))

  sim = simulate_fun_separate(nday,indexstates,indexonset*7,r0,d_commuting,d_airtravel,pops_raw,1)
  sim_onlycommuting = simulate_fun_separate(nday,indexstates,indexonset*7,r0,mob[[4]],d_airtravel*0,pops_raw,1)
  sim_onlyairtravel = simulate_fun_separate(nday,indexstates,indexonset*7,r0,mob[[3]],d_airtravel*0,pops_raw,1)
  
  pcts = getPcts(sim)
  pcts_commuting = getPcts(sim_onlycommuting)
  pcts_airtravel = getPcts(sim_onlyairtravel)
  
  mrg = theme(plot.margin = margin(c(4,-3,-3,-4),'cm'))
  
  sims = ggarrange(plotCluster_sim(pcts[[1]][,1],pcts[[2]][,1],indexonset[1],fctr,clusters)+mrg,
                   plotCluster_sim(pcts[[1]][,2],pcts[[2]][,2],indexonset[2],fctr,clusters)+mrg,ncol=2)
  
  trues = ggarrange(plotCluster_true(clusters[1],indexonset[1],fctr,clusters)+mrg,
                    plotCluster_true(clusters[2],indexonset[2],fctr,clusters)+mrg,ncol=2)
  
  sims_commuting = ggarrange(plotCluster_sim(pcts_commuting[[1]][,1],pcts_commuting[[2]][,1],indexonset[1],fctr,clusters)+mrg,
                             plotCluster_sim(pcts_commuting[[1]][,2],pcts_commuting[[2]][,2],indexonset[2],fctr,clusters)+mrg,ncol=2)
  
  sims_airtravel = ggarrange(plotCluster_sim(pcts_airtravel[[1]][,1],pcts_airtravel[[2]][,1],indexonset[1],fctr,clusters)+mrg,
                             plotCluster_sim(pcts_airtravel[[1]][,2],pcts_airtravel[[2]][,2],indexonset[2],fctr,clusters)+mrg,ncol=2)
  
  plt = ggarrange(trues,sims_commuting,sims_airtravel,ncol=1)
  return(plt)
}


plotFigure3 <- function(fctr){
  leg = as_ggplot(get_legend(plotCluster_true(3,0,0,c()) + theme(legend.position='top',
                                                                 legend.text=element_text(size=5),
                                                                 legend.title=element_text(size=5),
                                                                 legend.key.height=unit(4,'pt'),
                                                                 legend.key.width=unit(15,'pt')) +
                               guides(fill = guide_colourbar(title.position="top", title.hjust = 0.5))))
  
  treeplot1 = plotTree_vert(clusterpath,"H3N2",2018,mtdt,c(670,671),c(0.5,2.3))[[1]]
  treeplot2 = plotTree_vert(clusterpath,"H1N1",2017,mtdt,c(96,104),c(0.5,2.3))[[1]]
  
  plt1 = plotTrueAndSimulated_2018(fctr)
  plt2 = plotTrueAndSimulated_2017(fctr) 
  
  ggarrange(
    ggarrange(
      ggarrange(treeplot1,plt1,heights=c(0.25,0.8),ncol=1),ggplot() + theme_void(),ggarrange(treeplot2,plt2,heights=c(0.25,0.8),ncol=1),ncol=3,widths=c(0.4,0.05,0.4)),
    ggarrange(ggplot()+theme_void(),leg,ggplot()+theme_void(),nrow=1),ncol=1,heights=c(0.9,0.15))
  ggsave("Figure_3.pdf",width=100,height=75,units='mm')  #sims = (plotCluster_sim(pcts[[1]][,2],pcts[[2]][,2],0,20)|plotCluster_sim(pcts[[1]][,1],pcts[[2]][,1],2,20)) + plot_annotation(title='True') &theme(text=element_text(size=6,hjust=0.5))
}

plotFigure3(12)

plotTimingEffect_1 <- function(fctr){
  
  clusters = c(670,671)
  states = c(state.name[-c(2,11)],"District of Columbia")
  indexstates = match(c("Georgia","Nebraska"),states)-1
  
  indexonset = cluster_df_country$RelativeOnset_Country[match(clusters,cluster_df_country$Cluster)]
  indexonset = as.numeric(indexonset) - min(as.numeric(indexonset))
  sim = simulate_fun_separate(nday,indexstates,indexonset*7,r0,mob[[4]],d_airtravel*0,pops_raw,1)
  
  indexonset_later = indexonset
  indexonset_later[2] = indexonset_later[2] + 4
  indexonset_later = indexonset_later - min(indexonset_later)
  sim_later = simulate_fun_separate(nday,indexstates,indexonset_later*7,r0,mob[[4]],d_airtravel*0,pops_raw,1)
  
  indexonset_earlier = indexonset
  indexonset_earlier[1] = indexonset_earlier[1] + 4
  indexonset_earlier = indexonset_earlier - min(indexonset_earlier)
  sim_earlier = simulate_fun_separate(nday,indexstates,indexonset_earlier*7,r0,mob[[4]],d_airtravel*0,pops_raw,1)
  
  pcts = getPcts(sim)
  pcts_earlier = getPcts(sim_earlier)
  pcts_later = getPcts(sim_later)
  
  normal_sum = (colSums(pcts[[1]],na.rm=T))
  earlier_sum =(colSums(pcts_earlier[[1]],na.rm=T))
  later_sum = (colSums(pcts_later[[1]],na.rm=T))
  
  print(table(pcts[[1]][,2]>0.1))
  print(table(pcts_earlier[[1]][,2]>0.1))
  print(table(pcts_later[[1]][,2]>0.1))
  
  mrg = theme(plot.margin = margin(c(4,-3,-3,-4),'cm'))
  
  sims_earlier = ggarrange(plotCluster_sim(pcts_earlier[[1]][,1],pcts_earlier[[2]][,1],indexonset_earlier[1],fctr,clusters)+mrg,
                           plotCluster_sim(pcts_earlier[[1]][,2],pcts_earlier[[2]][,2],indexonset_earlier[2],fctr,clusters)+mrg,ncol=2)
  
  sims_later = ggarrange(plotCluster_sim(pcts_later[[1]][,1],pcts_later[[2]][,1],indexonset_later[1],fctr,clusters)+mrg,
                         plotCluster_sim(pcts_later[[1]][,2],pcts_later[[2]][,2],indexonset_later[2],fctr,clusters)+mrg,ncol=2)
  
  ggarrange(sims_earlier,sims_later,ncol=1)
  
  ggsave("SuppFigTiming1.pdf",width=90,units="mm",height=60)
  
}


plotTimingEffect_2 <- function(fctr){
  
  clusters = c(96,104)
  states = c(state.name[-c(2,11)],"District of Columbia")
  indexstates = match(c("Mississippi","Oregon"),states)-1
  
  indexonset = cluster_df_country$RelativeOnset_Country[match(clusters,cluster_df_country$Cluster)]
  indexonset = as.numeric(indexonset) - min(as.numeric(indexonset))
  sim = simulate_fun_separate(nday,indexstates,indexonset*7,r0,mob[[4]],d_airtravel*0,pops_raw,1)
  
  indexonset_later = indexonset
  indexonset_later[2] = indexonset_later[2] + 4
  indexonset_later = indexonset_later - min(indexonset_later)
  sim_later = simulate_fun_separate(nday,indexstates,indexonset_later*7,r0,mob[[4]],d_airtravel*0,pops_raw,1)
  
  indexonset_earlier = indexonset
  indexonset_earlier[1] = indexonset_earlier[1] + 4
  indexonset_earlier = indexonset_earlier - min(indexonset_earlier)
  sim_earlier = simulate_fun_separate(nday,indexstates,indexonset_earlier*7,r0,mob[[4]],d_airtravel*0,pops_raw,1)
  
  pcts = getPcts(sim)
  pcts_earlier = getPcts(sim_earlier)
  pcts_later = getPcts(sim_later)
  
  normal_sum = (colSums(pcts[[1]],na.rm=T))
  earlier_sum =(colSums(pcts_earlier[[1]],na.rm=T))
  later_sum = (colSums(pcts_later[[1]],na.rm=T))

  
  mrg = theme(plot.margin = margin(c(4,-3,-3,-4),'cm'))
  
  sims_earlier = ggarrange(plotCluster_sim(pcts_earlier[[1]][,1],pcts_earlier[[2]][,1],indexonset_earlier[1],fctr,clusters)+mrg,
                           plotCluster_sim(pcts_earlier[[1]][,2],pcts_earlier[[2]][,2],indexonset_earlier[2],fctr,clusters)+mrg,ncol=2)
  
  sims_later = ggarrange(plotCluster_sim(pcts_later[[1]][,1],pcts_later[[2]][,1],indexonset_later[1],fctr,clusters)+mrg,
                         plotCluster_sim(pcts_later[[1]][,2],pcts_later[[2]][,2],indexonset_later[2],fctr,clusters)+mrg,ncol=2)
  
  ggarrange(sims_earlier,sims_later,ncol=1)
  
  ggsave("SuppFigTiming2.pdf",width=90,units="mm",height=60)
  
}

plotTimingEffect_1(12)
plotTimingEffect_2(12)


