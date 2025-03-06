library(scales)

jumpFrequency_analysis <- function(all_d){ 
  
  mob = getMobilityMatrix()
  
  states = c(state.name[-c(2,11)],"District Of Columbia")
  states_abb = c(state.abb[-c(2,11)],"DC")
  
  relative_jump_frequencies = lapply(states,function(y) ((sapply(states,function(x)sum((all_d[all_d[,2]==y,3]==x),na.rm=T))+sapply(states,function(x)sum((all_d[all_d[,3]==y,2]==x),na.rm=T))))/(nrow(all_d[all_d[,2]==y,])+nrow(all_d[all_d[,3]==y,])))
  relative_jump_frequencies = do.call(rbind,relative_jump_frequencies)
  rownames(relative_jump_frequencies) = colnames(relative_jump_frequencies)
  
  relative_jump_frequencies_long = relative_jump_frequencies
  relative_jump_frequencies_long = melt(relative_jump_frequencies_long)
  colnames(relative_jump_frequencies_long) = c("state1","state2",'prop')
  relative_jump_frequencies_long$dir = paste0(states_abb[match(relative_jump_frequencies_long$state2,states)]," \u2192 ",states_abb[match(relative_jump_frequencies_long$state1,states)])
  
  comm = mob[[4]]
  comm = apply(comm,2,function(x)x/sum(x,na.rm=T))
  
  relative_jump_frequencies_long$comm_rel = melt(t(comm))[,3]
  
  air = mob[[3]]
  air = apply(air,2,function(x)x/sum(x,na.rm=T))
  
  relative_jump_frequencies_long$airtravel = melt(air)[,3]
  
  dist = mob[[6]]
  relative_jump_frequencies_long$dist = unlist(as.vector(dist))
    
  adjacency_mat = mob[[5]]
  relative_jump_frequencies_long$adj = NA
  for (stateidx in 1:length(states)){
    state = states[stateidx]
    relative_jump_frequencies_long[relative_jump_frequencies_long$state2==state,]$adj = as.vector(unlist(adjacency_mat[stateidx,]))
  }
  
  pl_1 = ggplot(relative_jump_frequencies_long[order(relative_jump_frequencies_long$prop,decreasing=T),][1:25,],
                aes(x=reorder(dir,prop),y=prop))+ geom_point(pch=22,color='black',aes(fill=factor(adj)),cex=2,stroke=0.1) + theme_bw()+ 
    theme(axis.text.x = element_text(color='black',size=6,angle = 90, vjust = 0.5, hjust=1), panel.grid.major = element_line(size=0.2),
          panel.border=element_rect(linewidth=0.2),
          axis.ticks = element_line(size=0.2),axis.text.y=element_text(color='black')) + 
    theme(axis.title.x=element_blank(),axis.text.y=element_text(size=7),axis.title = element_text(size=7),axis.line=element_blank()) +
    ylab("Relative jump\ncontribution") + scale_fill_manual(values=c("#69b3a2", "#404080")) + theme(legend.position='none') +
    scale_y_continuous(expand=expansion(mult=0.1))
  
  
  jumpmobcorplot_comm = ggplot(relative_jump_frequencies_long[relative_jump_frequencies_long$prop>0,]) + 
    geom_point(aes(fill=factor(adj),y=as.numeric(prop),x=(as.numeric(comm_rel))),pch=21,cex=1,stroke=0.05) + 
    scale_y_log10(breaks=c(0.01,0.1),labels = label_log(digits = 2)) + scale_x_log10(breaks=c(10^seq(-5,0)),labels=label_log(digits = 2)) +
    annotation_logticks(side='bl',outside=T,short=unit(0.7,'mm'),mid=unit(0.7,'mm'),long=unit(0,'mm'),size=0.2) +
    coord_cartesian(clip='off') +
    thm + ylab("Relative jump contribution") + xlab("Relative commuting contribution") + scale_fill_manual(values=c("#69b3a2", "#404080"))
  
  print(cor.test(log(relative_jump_frequencies_long[relative_jump_frequencies_long$prop>0,]$prop),log(relative_jump_frequencies_long[relative_jump_frequencies_long$prop>0,]$comm_rel),method='spearman'))
  print(cor.test(log(relative_jump_frequencies_long[relative_jump_frequencies_long$prop>0,]$prop),log(relative_jump_frequencies_long[relative_jump_frequencies_long$prop>0,]$airtravel),method='spearman'))
    
  jumpmobcorplot_air = ggplot(relative_jump_frequencies_long[relative_jump_frequencies_long$prop>0,]) + 
    geom_point(aes(fill=factor(adj),y=as.numeric(prop),x=(as.numeric(airtravel))),pch=21,cex=1,stroke=0.05) + 
    scale_y_log10(breaks=c(0.01,0.1),labels=label_log(digits = 2)) + scale_x_log10(breaks=c(1e-7, 1e-5,1e-3,1e-1,1e1),labels=label_log(digits = 2)) +
    annotation_logticks(side='bl',outside=T,short=unit(0.7,'mm'),mid=unit(0.7,'mm'),long=unit(0,'mm'),size=0.2) +
    coord_cartesian(clip='off') +
    thm + ylab("Relative jump contribution") + xlab("Relative air travel contribution") + #+ stat_smooth(linewidth=0.3,lty=2,se=F,col='darkgrey',method='lm',aes(x=as.numeric(jumps),y=(as.numeric(mobility))))
    scale_fill_manual(values=c("#69b3a2", "#404080"))
  
  
  jumpmobcorplot_dist = ggplot(relative_jump_frequencies_long[relative_jump_frequencies_long$prop>0,]) + 
    geom_point(aes(fill=factor(adj),y=as.numeric(prop),x=(as.numeric(dist))),pch=21,cex=.7,stroke=0.1) + 
    scale_y_log10(breaks=c(0.01,0.1),labels=label_log(digits = 2)) + scale_x_log10(breaks=c(1,10,100,1000,10000),labels=label_log(digits = 2)) +
    annotation_logticks(side='bl',outside=T,short=unit(0.7,'mm'),mid=unit(0.7,'mm'),long=unit(0,'mm'),size=0.2) +
    coord_cartesian(clip='off') +
    thm + ylab("Relative jump contribution") + xlab("Relative air travel contribution") + #+ stat_smooth(linewidth=0.3,lty=2,se=F,col='darkgrey',method='lm',aes(x=as.numeric(jumps),y=(as.numeric(mobility))))
    scale_fill_manual(values=c("#69b3a2", "#404080"))
  
  relative_connectedness = relative_jump_frequencies
  for (i in 1:nrow(relative_connectedness)){relative_connectedness[,i] = relative_connectedness[,i]/mean(relative_connectedness[,i],na.rm=T)}
  relative_connectedness = (relative_connectedness + t(relative_connectedness))/2
  relative_connectedness[lower.tri(relative_connectedness)] = NA
  relative_connectedness_long = melt(relative_connectedness)
  colnames(relative_connectedness_long) = c("state1","state2",'prop')
  relative_connectedness_long$dir = paste0(states_abb[match(relative_connectedness_long$state2,states)]," ↔ ",states_abb[match(relative_connectedness_long$state1,states)])
  comm = mob[[4]]
  comm = apply(comm,2,function(x)x/sum(x,na.rm=T))
  relative_connectedness_long$comm_rel = melt(comm)[,3]
  
  relative_connectedness_long$adj = NA
  for (stateidx in 1:length(states)){
    state = states[stateidx]
    relative_connectedness_long[relative_connectedness_long$state2==state,]$adj = as.vector(unlist(adjacency_mat[stateidx,]))
  }
  
  relative_connectedness_long$dist = unlist(as.vector(dist))
  
  m_dist = (cor.test(log(relative_connectedness_long[relative_connectedness_long$prop>0,]$prop),log(relative_connectedness_long[relative_connectedness_long$prop>0,]$dist),method='pearson'))
  pl_2 = ggplot(relative_connectedness_long[order(relative_connectedness_long$prop,decreasing=T),][1:25,],aes(x=reorder(dir,prop),y=prop))+ 
    geom_point(pch=22,color='black',aes(fill=factor(adj)),cex=2,stroke=0.1) + theme_bw()+ 
    theme(axis.text.x = element_text(color='black',size=6,angle = 90, vjust = 0.5, hjust=1), panel.grid.major = element_line(size=0.2),
          panel.border=element_rect(linewidth=0.2),
          axis.ticks = element_line(size=0.2),axis.text.y=element_text(color='black')) + 
    theme(axis.title.x=element_blank(),axis.text.y=element_text(size=7),axis.title = element_text(size=7),axis.line=element_blank()) +
    ylab("Normalized jump\nfrequency") + scale_fill_manual(values=c("#69b3a2", "#404080")) + theme(legend.position='none') +
    scale_y_continuous(expand=expansion(mult=0.1))
  
  diffhist = ggplot(relative_connectedness_long[relative_connectedness_long$prop>0,],aes(fill=factor(adj),x=prop)) + geom_density( color='black', alpha=0.6, position = 'identity',linewidth=0) +     
    scale_fill_manual(values=c("#69b3a2", "#404080"),labels=c('Non-adjoining','Adjoining')) + thm + theme(legend.position='top') + xlab("Normalized jump frequency") + ylab("Density") + 
    theme(legend.direction='vertical') + theme(legend.title = element_blank()) + theme(legend.key.size = unit(.3,"cm"))
  
  return(list(pl_1,pl_2,jumpmobcorplot_comm,jumpmobcorplot_air,diffhist,jumpmobcorplot_dist,m_dist))
}

plts = jumpFrequency_analysis(all_sourcesink)

leg = as_ggplot(get_legend(plts[[5]] + theme(legend.position='top',
                                             legend.text=element_text(size=7),
                                             legend.direction='horizontal',
                                             legend.key.height=unit(7,'pt'),
                                             legend.key.width=unit(7,'pt'),
                                             legend.title=element_blank())))

ggarrange(ggarrange((plts[[3]]/plts[[4]]),(plts[[1]]/plts[[2]]),plts[[5]]+theme(legend.position='none'),widths=c(0.35,0.55,0.3),ncol=3),leg,nrow=2,heights=c(0.9,0.1)) 
ggsave("Figure_5_resubmit.pdf",width=2080,height=1000,units='px',dpi=320)

by_subtype = lapply(unique(cluster_df_country$Subtype),function(x)jumpFrequency_analysis(all_sourcesink[all_sourcesink$cluster%in%cluster_df_country[cluster_df_country$Subtype==x,]$Cluster,]))
(by_subtype[[1]][[5]]+ggtitle("H3N2"))+(by_subtype[[2]][[5]]+ggtitle("H1N1pdm09"))+(by_subtype[[3]][[5]]+ggtitle("B/Vic"))+(by_subtype[[4]][[5]]+ggtitle("B/Yam"))+plot_layout(guides='collect') & 
  theme(legend.position = 'right',plot.title = element_text(size=6))
ggsave("SuppFig_adjoining.pdf",width=1,units="mm",height=120)

single_season_cor_list = list()
idx = 1
combs = unique(interaction(cluster_df_country$Season,cluster_df_country$Subtype,sep='-'))
for (i in 1:length(combs)){
  comb = combs[i]
  lst = list()
  lst$subtype = substr(comb,6,nchar(as.character(comb)))
  lst$season = substr(comb,1,4)
  if (lst$season %in%c(2020,2021,2023)){next}
  df = all_sourcesink[all_sourcesink$cluster%in%cluster_df_country[cluster_df_country$Subtype==lst$subtype & cluster_df_country$Season==lst$season,]$Cluster,]
  if (nrow(df)<10000){next}
  lst$df = jumpFrequency_analysis(df[sample(nrow(df),10000,replace=T),])
  single_season_cor_list[[idx]] = lst
  idx = idx + 1
}

per_season_df = as.data.frame(do.call(rbind,lapply(single_season_cor_list,function(x)c(x$df[[7]]$conf.int[1],x$df[[7]]$estimate,x$df[[7]]$conf.int[2]))))
per_season_df = cbind(per_season_df,unlist(lapply(single_season_cor_list,function(x)x$subtype)))
per_season_df = cbind(per_season_df,unlist(lapply(single_season_cor_list,function(x)x$season)))
per_season_df = as.data.frame(per_season_df)
colnames(per_season_df) = c("lo","estimate","hi","subtype","season")
per_season_df$season = factor(per_season_df$season,levels=unique(per_season_df$season),labels=paste0(unique(per_season_df$season),"/",as.numeric(unique(per_season_df$season))+1))

ggplot(per_season_df) + geom_point(aes(x=interaction(season,subtype,sep="\n"),y=estimate)) + geom_errorbar(aes(x=interaction(season,subtype,sep="\n"),ymin=lo,ymax=hi),linewidth=0.3) +
  theme_bw() + thm + xlab("Season-Subtype") + ylab("Correlation")
ggsave("SuppFig_distancebyseason.pdf",width=160,units="mm",height=60)


