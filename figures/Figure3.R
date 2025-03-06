library(ggpubr)
library(usmap)
source("plotFunctions.R")
source("simulateMobility.R")

getOptimalOnsets <- function(clusters,indexonset,seqcounts,mob1,mob2){
  
  sim = simulate_fun_separate(nday,indexstates,indexonset*7,r0,mob1,mob2,pops_raw,1)
  pcts = getPcts(sim,seqcounts)
  
  diff = sapply(1:length(clusters),function(x)sum((pcts[[1]][,x]-(seqcounts[,x]/rowSums(seqcounts))),na.rm=T))
  poss_range = lapply(1:length(clusters),function(x)c(0,1,2)*sign(diff[x]))
  diffs = expand.grid(poss_range)
  
  cur_ll = -100000
  cur_best = NA
  
  for (i in 1:nrow(diffs)){
    
    indons = (indexonset+unlist(diffs[i,]))
    indons = as.numeric(indons) - min(as.numeric(indons))
    sim = simulate_fun_separate(nday,indexstates,indons*7,r0,mob1,mob2,pops_raw,1)
    pcts = getPcts(sim,seqcounts)
    
    if (pcts[[3]] > cur_ll){
      cur_ll = pcts[[3]]
      cur_best = i
    }
    
  }
  
  return(diffs[cur_best,])
}


plotTrueAndSimulated <- function(fctr, clusters, indexstates){
  
 
  states = c(state.name[-c(2,11)],"District of Columbia")
  
  indexonset = as.numeric(cluster_df_country$RelativeOnset_Country[match(clusters,cluster_df_country$Cluster)])
  indexonset = as.numeric(indexonset) - min(as.numeric(indexonset))
  
  seqcounts = t(sapply(states,function(y)sapply(clusters,function(x)
    ifelse(length(cluster_df[cluster_df$Cluster==x & cluster_df$State==y,]$Seqcount)==1,
           cluster_df[cluster_df$Cluster==x & cluster_df$State==y,]$Seqcount,
           0))))
  
  onsets = t(sapply(states,function(y)sapply(clusters,function(x)cluster_df[cluster_df$Cluster==x & cluster_df$State==y,]$Onset)))
  
  time_adj_commuting = getOptimalOnsets(clusters,indexonset,seqcounts,mob[[4]],mob[[4]]*0)
  time_adj_air = getOptimalOnsets(clusters,indexonset,seqcounts,mob[[3]],mob[[4]]*0)
  
  onset_comm = unlist(indexonset + time_adj_commuting)
  onset_air = unlist(indexonset + time_adj_air)
  
  onset_comm = onset_comm-min(onset_comm)
  onset_air = onset_air-min(onset_air)
  
  sim_onlycommuting = simulate_fun_separate(nday,indexstates,onset_comm*7,r0,mob[[4]],d_airtravel*0,pops_raw,1)
  sim_onlyairtravel = simulate_fun_separate(nday,indexstates,onset_air*7,r0,mob[[3]],d_airtravel*0,pops_raw,1)
  sim_onlycommuting_later = simulate_fun_separate(nday,indexstates,(onset_comm+c(4,0))*7,r0,mob[[4]],d_airtravel*0,pops_raw,1)
  sim_onlycommuting_earlier = simulate_fun_separate(nday,indexstates,(onset_comm+c(0,4))*7,r0,mob[[4]],d_airtravel*0,pops_raw,1)
    
  pcts_commuting = getPcts(sim_onlycommuting,seqcounts)
  pcts_airtravel = getPcts(sim_onlyairtravel,seqcounts)
  pcts_earlier = getPcts(sim_onlycommuting_earlier,seqcounts)
  pcts_later = getPcts(sim_onlycommuting_later,seqcounts)
  
  print(mean((as.numeric(pcts_commuting[[2]]+(lm(as.numeric(onsets)~1+offset(as.numeric(pcts_commuting[[2]]))))$coefficients[1])-(as.numeric(onsets)))^2,na.rm=T))
  print(mean((as.numeric(pcts_airtravel[[2]]+(lm(as.numeric(onsets)~1+offset(as.numeric(pcts_airtravel[[2]]))))$coefficients[1])-(as.numeric(onsets)))^2,na.rm=T))
  
  trues = ggarrange(plotCluster_true(clusters[1],indexonset[1],fctr,clusters)+mrg,
                    plotCluster_true(clusters[2],indexonset[2],fctr,clusters)+mrg,ncol=2)

  sims_commuting = ggarrange(plotCluster_sim(pcts_commuting[[1]][,1],pcts_commuting[[2]][,1],onset_comm[1],fctr,clusters)+mrg,
                             plotCluster_sim(pcts_commuting[[1]][,2],pcts_commuting[[2]][,2],onset_comm[2],fctr,clusters)+mrg,ncol=2)
  
  sims_airtravel = ggarrange(plotCluster_sim(pcts_airtravel[[1]][,1],pcts_airtravel[[2]][,1],onset_air[1],fctr,clusters)+mrg,
                             plotCluster_sim(pcts_airtravel[[1]][,2],pcts_airtravel[[2]][,2],onset_air[2],fctr,clusters)+mrg,ncol=2)
  
  sims_earlier = ggarrange(plotCluster_sim(pcts_earlier[[1]][,1],pcts_earlier[[2]][,1],onset_air[1],fctr,clusters)+mrg,
                             plotCluster_sim(pcts_earlier[[1]][,2],pcts_earlier[[2]][,2],onset_air[2],fctr,clusters)+mrg,ncol=2)
  
  sims_later = ggarrange(plotCluster_sim(pcts_later[[1]][,1],pcts_later[[2]][,1],onset_air[1],fctr,clusters)+mrg,
                             plotCluster_sim(pcts_later[[1]][,2],pcts_later[[2]][,2],onset_air[2],fctr,clusters)+mrg,ncol=2)
  
  plt = ggarrange(trues,sims_commuting,sims_airtravel,ncol=1)
  plt_timing = ggarrange(sims_earlier,sims_later,ncol=1)
  return(list(plt,plt_timing))
}



plotFigure3 <- function(fctr){
  leg = as_ggplot(get_legend(plotCluster_true(3,0,0,c()) + theme(legend.position='top',
                                                                 legend.text=element_text(size=7),
                                                                 legend.title=element_text(size=7),
                                                                 legend.key.height=unit(2.5,'pt'),
                                                                 legend.key.width=unit(15,'pt')) +
                               guides(fill = guide_colourbar(title.position="top", title.hjust = 0.5))))
  
  treeplot1 = plotTree_vert(clusterpath,"H3N2",2018,mtdt,c(670,671),c(0.5,2.4))[[1]]
  treeplot2 = plotTree_vert(clusterpath,"H1N1",2017,mtdt,c(96,104),c(0.5,2.4))[[1]]
  
  treeplot1 = treeplot1 + scale_y_continuous(expand=c(.03,0))
  treeplot2 = treeplot2 + scale_y_continuous(expand=c(.03,0))
  
  plt_2018 = plotTrueAndSimulated(fctr,c(670,671), match(c("Georgia","Nebraska"),states)-1)
  plt_2017 = plotTrueAndSimulated(fctr,c(96,104), match(c("Mississippi","Oregon"),states)-1)

  plt1 = plt_2018[[1]]
  plt2 = plt_2017[[1]]
  plt_timing1 = plt_2018[[2]]
  plt_timing2 = plt_2017[[2]]
  plt_timing1/plt_timing2
  ggsave("SuppFig_Timing_Both.pdf",width=90,units="mm",height=120)
  
  ggarrange(
    ggarrange(
      ggarrange(treeplot1,plt1,heights=c(0.3,0.8),ncol=1),ggplot() + theme_void(),ggarrange(treeplot2,plt2,heights=c(0.3,0.8),ncol=1),ncol=3,widths=c(0.4,0.05,0.4)),
    ggarrange(ggplot()+theme_void(),leg,ggplot()+theme_void(),nrow=1,widths=c(0.5,0.4,0.4)),ncol=1,heights=c(0.9,0.12))
  ggsave("Figure_3_resubmit.pdf",width=1700,height=1700*100/140,units='px',dpi=320)  
  
}

plotFigure3(13)




