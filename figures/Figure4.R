leg = as_ggplot(get_legend(plotCluster_true(3,0,0,c()) + theme(legend.position='top',legend.text=element_text(size=7),legend.title=element_text(size=7),legend.key.height=unit(2.5,'pt'),legend.key.width=unit(15,'pt')) +
                             guides(fill = guide_colourbar(title.position="top", title.hjust = 0.5))))

season = 2019
subtype = "Vic"
cdc_2019 = cluster_df_country[cluster_df_country$Season==2019 & cluster_df_country$Subtype=="Vic",]
cdc_2019 = cdc_2019[order(cdc_2019$Seqcount,decreasing=T),]
clusters = cdc_2019$Cluster[1:7]
clusters = clusters[!(clusters=="UNCLUSTERED")]

onsetstates = c("California","Florida","Texas","Louisiana","Nevada","Washington")
indexstates = match(onsetstates,states)-1

seqcounts = t(sapply(states,function(y)sapply(clusters,function(x)
  ifelse(length(cluster_df[cluster_df$Cluster==x & cluster_df$State==y,]$Seqcount)==1,
         cluster_df[cluster_df$Cluster==x & cluster_df$State==y,]$Seqcount,
         0))))

onsets = t(sapply(states,function(y)sapply(clusters,function(x)cluster_df[cluster_df$Cluster==x & cluster_df$State==y,]$Onset)))

mrg = theme(plot.margin = margin(c(-1,-3,-3,-4),'cm'))

indexonset = unlist(lapply(clusters,function(x)min(cluster_df[cluster_df$Cluster==x,]$Onset,na.rm=T)))
indexonset = as.numeric(indexonset) - min(as.numeric(indexonset))

treeplot = plotTree_vert(clusterpath,"Vic",2019,mtdt,clusters,c(0.5,2.1),F)[[1]]
treeplot = treeplot + scale_y_continuous(expand=c(.01,0)) + theme()


simulate_2019_vic <- function(fctr){
  
  output_list = list()
  
  plot_list = list()
  
  mob_combs = list(list(d_commuting,d_airtravel),list(mob[[4]],d_airtravel*0),list(mob[[3]],d_airtravel*0))
  
  for (mob_list_idx in 1:length(mob_combs)){
    
    mob_list = mob_combs[[mob_list_idx]]
    
    time_adj = getOptimalOnsets(clusters,indexonset,seqcounts,mob_list[[1]],mob_list[[2]])
    
    indexonset_adj = unlist(indexonset+time_adj)
    
    if (mob_list_idx == 1){time_adj_combined = indexonset_adj}
    
    sim = simulate_fun_separate(nday,indexstates,(indexonset_adj)*7,r0,mob_list[[1]],mob_list[[2]],pops_raw,1)
    
    pcts = getPcts(sim,seqcounts)
    
    mse = (mean((as.numeric(pcts[[2]]+(lm(as.numeric(onsets)~1+offset(as.numeric(pcts[[2]]))))$coefficients[1])-(as.numeric(onsets)))^2,na.rm=T))
    
    pcts[[4]] = mse
    
    output_list[[mob_list_idx]] = pcts
    
    sim_plots = list()
    for (i in 1:length(clusters)){
      sim_plots[[i]] = plotCluster_sim(pcts[[1]][,i],pcts[[2]][,i],indexonset_adj[i]-min(indexonset_adj),fctr,clusters)+mrg
    }
    
    plots =  ggarrange(plotlist=sim_plots,nrow=1)
    
    plot_list[[mob_list_idx]] = ggarrange(ggplot()+theme_void(),plots,nrow=1,widths=c(0.07,0.95))
    
  }
  
  true_plots = list()
  for (i in 1:length(clusters)){
    true_plots[[i]] = plotCluster_true(clusters[i],indexonset[i]-min(indexonset),fctr,clusters)+mrg
  }
  plots =  ggarrange(plotlist=true_plots,nrow=1)
  plot_list[[4]] = ggarrange(ggplot()+theme_void(),plots,nrow=1,widths=c(0.07,0.95))
  
  titleplots = list()
  for (i in 1:length(clusters)){
    titleplots[[i]] = ggplot()+theme_void() + annotate("text",x=0,y=0,label=i,size=7/.pt) + xlim(c(-0.5,0.5)) + ylim(c(-0.5,0.5))
  }
  titlerow = ggarrange(ggplot()+theme_void(),ggarrange(plotlist=titleplots,nrow=1),widths=c(0.07,0.95))
  
  ggarrange(treeplot,titlerow,plot_list[[4]],plot_list[[1]],ggarrange(ggplot()+theme_void(),leg,ggplot()+theme_void(),nrow=1,widths=c(0.5,0.3,0.4)),nrow=5,heights=c(0.45,0.05,0.35,0.35,0.20))
  ggsave("Figure_4_resubmit.pdf",width=2080,height=900,units='px',dpi=320)  
  
  output_list[[4]] = fitGrav(time_adj_combined, T)
  
  sim_plots = list()
  for (i in 1:length(clusters)){
    sim_plots[[i]] = plotCluster_sim(output_list[[4]][[1]][,i],output_list[[4]][[2]][,i],time_adj_combined[i]-min(time_adj_combined),fctr,clusters)+mrg
  }
  
  plots =  ggarrange(plotlist=sim_plots,nrow=1)
  plot_list[[5]] = ggarrange(ggplot()+theme_void(),plots,nrow=1,widths=c(0.07,0.95))
  
  
  output_list[[5]] = fitGrav(time_adj_combined, F)
  
  sim_plots = list()
  for (i in 1:length(clusters)){
    sim_plots[[i]] = plotCluster_sim(output_list[[5]][[1]][,i],output_list[[5]][[2]][,i],time_adj_combined[i]-min(time_adj_combined),fctr,clusters)+mrg
  }
  
  plots =  ggarrange(plotlist=sim_plots,nrow=1)
  plot_list[[6]] = ggarrange(ggplot()+theme_void(),plots,nrow=1,widths=c(0.07,0.95))
  
  plot_list2 = list(plot_list[[4]],plot_list[[1]],plot_list[[2]],plot_list[[3]],plot_list[[6]],plot_list[[5]])
  
  ggarrange(ggarrange(plotlist=plot_list2,nrow=6),ggarrange(ggplot()+theme_void(),leg,ggplot()+theme_void(),nrow=1,widths=c(0.5,0.3,0.4)),ncol=1,heights=c(0.95,0.08))
  ggsave("SuppFig_2019.pdf",width=2080,height=1500,units='px',dpi=320)  
  
  return(plt)
}


fitGrav <- function(indexonset_adj, single_efflux){
  
  distparams = c()
  popparams = c()
  
  lls = c()
  
  distance = read.csv("/Users/simondejong/usa/centroid_dists.csv",h=F)
  
  cur_max_ll = -10000
  
  for (distparam in seq(-0.5,3,0.1)){
    
    print(distparam)
    
    for (popparam in seq(-0.5,2,0.1)){
      
      distparams = c(distparams,distparam)
      popparams = c(popparams,popparam)
      
      dmat = (distance^-distparam)*(distance>0)
      colnames(dmat) = state.abb
      mmat = do.call(cbind,lapply(1:50,function(x)dmat[,x]/sum(dmat[,x],na.rm=T)))
      mmat = mmat[-c(2,11),-c(2,11)]
      mmat = cbind(mmat,0)
      mmat = rbind(mmat,0)
      mmat=do.call(cbind,lapply(1:49,function(x)(pops_raw^popparam*mmat[,x])))
      colnames(mmat) = colnames(mob[[4]])
      rownames(mmat) = colnames(mob[[4]])
      mmat[is.nan(mmat)] = 0
      
      for (i in 1:49){mmat[,i] = mmat[,i]/sum(mmat[,i])}
      if (single_efflux){
        for (i in 1:49){mmat[,i] = mmat[,i]*median(colSums(mob[[4]]+mob[[3]]))}
      } else {
        for (i in 1:49){mmat[,i] = mmat[,i]*(sum(mob[[4]][,i]) + sum(mob[[3]][,i]))}
      }
      
      mmat[is.nan(mmat)] = 0
      
      sim_grav = simulate_fun_separate(nday,indexstates,indexonset_adj*7,r0,mmat,d_airtravel*0,pops_raw,1)
      pcts = getPcts(sim_grav,seqcounts)
      
      mse = (mean((as.numeric(pcts[[2]]+(lm(as.numeric(onsets)~1+offset(as.numeric(pcts[[2]]))))$coefficients[1])-(as.numeric(onsets)))^2,na.rm=T))
      pcts[[4]] = mse
      
      ll =  pcts[[3]]
    
      if (ll > cur_max_ll){max_pcts = pcts;cur_max_ll=ll}
      
      lls = c(lls,ll)
    }
  }
  
  lls_save = lls
  lls = lls-max(lls)
  
  llr = -2 * lls
  chisq = 3.84
  ci_idxes = which(llr<chisq)
  
  if (!(single_efflux)){
    
    print(quantile(popparams[ci_idxes]))
    print(quantile(distparams[ci_idxes]))
    
    print(popparams[which.max(lls)])
    print(distparams[which.max(lls)])
    
    grav_df = data.frame(pop=popparams,dist=distparams,ll=(1/abs(log(exp(lls)/sum(exp(lls))))),inci=((1:length(lls))%in%ci_idxes))
    
    ggplot(grav_df,aes(x=popparams,y=distparams,fill=inci,size=ll)) + geom_point(pch=21) + scale_fill_brewer(guide='none',palette='Set1') + theme_bw() +
      theme(legend.position='none') + xlab(expression(tau)) + ylab(expression(rho)) + theme(axis.text = element_text(size=10))
    
    ggsave("SuppFig_grav.pdf",width=120,units="mm",height=120)
  }
  
  return(max_pcts)
}

simulate_2019_vic(15)

params = read.csv("/Users/simondejong/US_phylo/2019_Vic.log",sep='\t')
date_decimal(mean(params$cluster_965.age.root[100:1000]))
date_decimal(mean(params$cluster_967.age.root[100:1000]))

date_decimal(quantile(params$cluster_965.age.root[100:1000],c(0.025,0.5,0.975)))
date_decimal(quantile(params$cluster_967.age.root[100:1000],c(0.025,0.5,0.975)))
date_decimal(min(mtdt[mtdt$Cluster==967,]$Date))


