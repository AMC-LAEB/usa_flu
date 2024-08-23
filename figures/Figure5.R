source("simulateMobility.R")
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
  relative_jump_frequencies_long$dir = paste0(states_abb[match(relative_jump_frequencies_long$state2,states)],"\u2192",states_abb[match(relative_jump_frequencies_long$state1,states)])
  
  comm = mob[[4]]
  comm = apply(comm,2,function(x)x/sum(x,na.rm=T))
  
  relative_jump_frequencies_long$comm_rel = melt(t(comm))[,3]
  
  air = mob[[3]]
  air = apply(air,2,function(x)x/sum(x,na.rm=T))
  
  relative_jump_frequencies_long$airtravel = melt(air)[,3]

  adjacency_mat = mob[[5]]
  relative_jump_frequencies_long$adj = NA
  for (stateidx in 1:length(states)){
    state = states[stateidx]
    relative_jump_frequencies_long[relative_jump_frequencies_long$state2==state,]$adj = as.vector(unlist(adjacency_mat[stateidx,]))
  }
  
  pl_1 = ggplot(relative_jump_frequencies_long[order(relative_jump_frequencies_long$prop,decreasing=T),][1:25,],aes(x=reorder(dir,prop),y=prop))+ geom_point(pch=21,aes(col=factor(adj)),cex=0.8) + theme_bw()+ 
    theme(axis.text.x = element_text(color='black',size=5,angle = 90, vjust = 0.5, hjust=1)) + theme(axis.title.x=element_blank(),axis.text.y=element_text(size=5),axis.title = element_text(size=5)) +
    ylab("Relative jump contribution") + scale_color_manual(values=c("#69b3a2", "#404080")) + theme(legend.position='none')
  
  
  jumpmobcorplot_comm = ggplot(relative_jump_frequencies_long[relative_jump_frequencies_long$prop>0,]) + 
    geom_point(aes(fill=factor(adj),y=as.numeric(prop),x=(as.numeric(comm_rel))),pch=21,cex=.7,stroke=0.1) + 
    scale_y_log10(breaks=c(0.01,0.1),labels = label_log(digits = 2)) + scale_x_log10(breaks=c(10^seq(-5,0)),labels=label_log(digits = 2)) +
    annotation_logticks(side='bl',outside=T,short=unit(0.7,'mm'),mid=unit(0.7,'mm'),long=unit(0,'mm'),size=0.1) +
    coord_cartesian(clip='off') +
    thm + ylab("Relative jump contribution") + xlab("Relative commuting contribution") + scale_fill_manual(values=c("#69b3a2", "#404080"))
  
  print(cor.test(log(relative_jump_frequencies_long[relative_jump_frequencies_long$prop>0,]$prop),log(relative_jump_frequencies_long[relative_jump_frequencies_long$prop>0,]$comm_rel),method='spearman'))
  print(cor.test(log(relative_jump_frequencies_long[relative_jump_frequencies_long$prop>0,]$prop),log(relative_jump_frequencies_long[relative_jump_frequencies_long$prop>0,]$airtravel),method='spearman'))
  
  jumpmobcorplot_air = ggplot(relative_jump_frequencies_long[relative_jump_frequencies_long$prop>0,]) + 
    geom_point(aes(fill=factor(adj),y=as.numeric(prop),x=(as.numeric(airtravel))),pch=21,cex=.7,stroke=0.1) + 
    scale_y_log10(breaks=c(0.01,0.1),labels=label_log(digits = 2)) + scale_x_log10(breaks=c(1e-7, 1e-5,1e-3,1e-1,1e1),labels=label_log(digits = 2)) +
    annotation_logticks(side='bl',outside=T,short=unit(0.7,'mm'),mid=unit(0.7,'mm'),long=unit(0,'mm'),size=0.1) +
    coord_cartesian(clip='off') +
    thm + ylab("Relative jump contribution") + xlab("Relative air travel contribution") + #+ stat_smooth(linewidth=0.3,lty=2,se=F,col='darkgrey',method='lm',aes(x=as.numeric(jumps),y=(as.numeric(mobility))))
    scale_fill_manual(values=c("#69b3a2", "#404080"))
  

  relative_connectedness = relative_jump_frequencies
  for (i in 1:nrow(relative_connectedness)){relative_connectedness[,i] = relative_connectedness[,i]/mean(relative_connectedness[,i],na.rm=T)}
  relative_connectedness = (relative_connectedness + t(relative_connectedness))/2
  relative_connectedness[lower.tri(relative_connectedness)] = NA
  relative_connectedness_long = melt(relative_connectedness)
  colnames(relative_connectedness_long) = c("state1","state2",'prop')
  relative_connectedness_long$dir = paste0(states_abb[match(relative_connectedness_long$state2,states)],"↔",states_abb[match(relative_connectedness_long$state1,states)])
  comm = mob[[4]]
  comm = apply(comm,2,function(x)x/sum(x,na.rm=T))
  relative_connectedness_long$comm_rel = melt(comm)[,3]
  
  relative_connectedness_long$adj = NA
  for (stateidx in 1:length(states)){
    state = states[stateidx]
    relative_connectedness_long[relative_connectedness_long$state2==state,]$adj = as.vector(unlist(adjacency_mat[stateidx,]))
  }
  
  pl_2 = ggplot(relative_connectedness_long[order(relative_connectedness_long$prop,decreasing=T),][1:25,],aes(x=reorder(dir,prop),y=prop))+ geom_point(pch=21,aes(col=factor(adj)),cex=0.8) + theme_bw()+ 
    theme(axis.text.x = element_text(color='black',size=5,angle = 90, vjust = 0.5, hjust=1)) + theme(axis.title.x=element_blank(),axis.text.y=element_text(size=5),axis.title = element_text(size=5)) +
    ylab("Normalized jump frequency") + scale_color_manual(values=c("#69b3a2", "#404080")) + theme(legend.position='none')
  
  
  diffhist = ggplot(relative_connectedness_long[relative_connectedness_long$prop>0,],aes(fill=factor(adj),x=prop)) + geom_density( color="#e9ecef", alpha=0.6, position = 'identity') +     
    scale_fill_manual(values=c("#69b3a2", "#404080"),labels=c('Non-adjoining','Adjoining')) + thm + theme(legend.position='top') + xlab("Normalized jump frequency") + ylab("Density") + 
    theme(legend.direction='vertical') + theme(legend.title = element_blank()) + theme(legend.key.size = unit(.3,"cm"))
  
  return(list(pl_1,pl_2,jumpmobcorplot_comm,jumpmobcorplot_air,diffhist))
}

plts = jumpFrequency_analysis(all_sourcesink)
ggarrange((plts[[3]]/plts[[4]]),(plts[[1]]/plts[[2]]),plts[[5]],widths=c(0.35,0.5,0.3),ncol=3) 
ggsave("Figure_5.pdf",width=120,units="mm",height=60)



