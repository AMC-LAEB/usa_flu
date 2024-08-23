library(ape)
library(ggplot2)
library(ggtree)
library(lubridate)
library(RColorBrewer)

plotTree <- function(clusterpath, subtype, year, mtdt,cltoplot=NA,vert=F){
  
  axissize=0.2
  
  trees = list()
  idx = 1
  
  large_clusters = unique(mtdt[mtdt$Subtype==subtype & mtdt$Season==year,]$File)
  large_clusters = large_clusters[!(is.na(large_clusters))]
  
  large_clusters = gsub(".txt",".tre",large_clusters)
  large_clusters = gsub("cluster","tree",large_clusters)
  
  
  for (i in large_clusters){
    
    t = read.nexus(paste0("clustering_output","/",clusterpath,"/",i))
    t$tip.label[grepl("cluster",t$tip.label)] = sub('_[^_]*$', '', t$tip.label[grepl("cluster",t$tip.label)])
    t = keep.tip(t,mtdt[mtdt$Strain%in%t$tip.label & mtdt$Date>year+0.5,]$Strain)
    maxdate = max(mtdt[mtdt$Strain %in% t$tip.label,]$Date,na.rm=T)
    height = max(node.depth.edgelength(t),na.rm=T)
    rootlen = maxdate-height-year+0.5
    
    if (rootlen < 0){
      
      tipskept = t$tip.label
      nodeorder = order(node.depth.edgelength(t))
      norder = nodeorder[nodeorder>length(t$tip.label)]
      norder = norder[abs(maxdate-(node.depth.edgelength(t)[norder]-height))-year+0.5 > 0]
      branchidx = 1
      
      while (length(norder)>0){
        
        subtree = extract.clade(t,norder[1])
        rootlen = maxdate-max(node.depth.edgelength(subtree),na.rm=T)-year+0.5# - 2
        subtree$root.edge = rootlen
        trees[[idx]] = subtree
        idx = idx + 1
        norder = norder[2:length(norder)]
        norder = norder[!(norder %in% subtree$edge[,2])]
      }
    }
    
    else {
      t$root.edge = rootlen
      trees[[idx]] = t
      idx = idx + 1
    }
  }
  
  w = rtree(1)
  w$edge.length=0
  t = bind.tree(w,trees[[1]])
  if (length(trees) > 1){
    for (tr in 2:length(trees)){
      t = bind.tree(t,trees[[tr]])
    }
  }

  tokeep = sample(1:length(t$tip.label),500,replace=T)#seq(1,length(t$tip.label),as.integer(length(t$tip.label)/500))
  t = drop.tip(t,setdiff(1:length(t$tip.label),tokeep))
  t = drop.tip(t,mtdt[mtdt$Strain %in% t$tip.label & mtdt$Date > year + 1.4,]$Strain)
  
  tr = t
  y =ggtree(t)
  y$data$isroot = 'black'
  y$data$isroot[y$data$parent == min(y$data$parent)] = 'black'#white'
  y$data$isroot[y$data$parent == min(y$data$parent)] = 'black'#white'
  y$data$size = 1
  y$data$size[y$data$label == "t1"] = NA
  
  y$data$cluster = 0
  for (i in 1:length(y$data$cluster)){
    if (length(mtdt[mtdt$Strain == y$data$label[i],]$Cluster) == 0){
      next
    }
    y$data$cluster[i] = mtdt[mtdt$Strain == y$data$label[i],]$Cluster[1]
  }
  y$data$clustersave = y$data$cluster
  y$data$cluster = match(y$data$cluster,unique(y$data$cluster))
 
  fig2 = y + aes(color=I(isroot),size=I(0.03))+ geom_tippoint(data=y$data[as.numeric(y$data$clustersave)>0 & (!(is.na(y$data$clustersave))) & y$data$clustersave!="UNCLUSTERED",],aes(fill=I(brewer.pal(100,"Paired")[as.numeric(clustersave)%%12+1]),size=I(size)*.7),pch=21,stroke=0.05) + theme_tree2() + 
    theme(axis.line=element_line(color='black',size = axissize),
          axis.text.x = element_text(size=7,color='black'),
          #plot.title = element_text(size = 7,hjust=0.5)) + ggtitle(paste0(c("A/H3N2","A/H1N1pdm09","B/Victoria","B/Yamagata")[match(subtype,c("H3N2","H1N1","Vic","Yam"))],'\n',year,"/",year+1))
          plot.title = element_text(size = 7,hjust=0.5)) + ggtitle(paste0(c("A/H3","A/H1pdm","B/Vic","B/Yam")[match(subtype,c("H3N2","H1N1","Vic","Yam"))],", ","'",substr(year,3,4),"/'",substr(year+1,3,4)))
  
  
  idx = 1
  for (clustertoplot in cltoplot){
    mrca = getMRCA(t,mtdt[mtdt$Cluster==clustertoplot & mtdt$Strain%in%t$tip.label,]$Strain)
    cl = y$data$cluster[y$data$label%in%mtdt[mtdt$Cluster==clustertoplot & mtdt$Strain%in%t$tip.label,]$Strain][1]
    fig2 = fig2 + geom_cladelab(node=mrca, label=idx, align=TRUE,  offset = .05, textcolor='black', barcolor=cl+10,fontsize=5/.pt)
    idx = idx + 1
  }
  
  fig2_plot = fig2 
  xmax = max(fig2_plot$data$x)
  xmin = xmax-1.5
  maxdate = max(mtdt[mtdt$Strain%in%fig2_plot$data$label,]$Date)
  max_x = max(fig2_plot$data$x[!(is.na(fig2_plot$data$label))])
 
  fig2_plot = fig2_plot + xlim(c(xmin,xmax))
  fig2_plot = fig2_plot + scale_x_continuous(labels=c(paste0("Jul '",substr(year,3,4)),paste0("Jul '",substr(year+1,3,4))),breaks=max_x+year+1+c(-0.5,0.5)-maxdate,limits=c(.8,2)) + theme(axis.text.x = element_text(hjust=c(0.5,1)))
  fig2_plot = fig2_plot + theme(axis.text.x = element_text(size=5),plot.title=element_text(size=5))+ theme(plot.title=element_text(vjust=.5,size=5))
  
  return(fig2_plot+ theme(axis.text=element_text(size=5)) + ylim(c(1,max(fig2_plot$data$y)+1)))
}

plotTree_vert <- function(clusterpath, subtype, year, mtdt, cltoplot=NA, lims, vert=F){
  
  trees = list()
  idx = 1
  
  large_clusters = unique(mtdt[mtdt$Subtype==subtype & mtdt$Season==year,]$File)
  large_clusters = large_clusters[!(is.na(large_clusters))]
  
  large_clusters = gsub(".txt",".tre",large_clusters)
  large_clusters = gsub("cluster","tree",large_clusters)
  
  axissize=0.2
  for (i in large_clusters){
    
    t = read.nexus(paste0("clustering_output","/",clusterpath,"/",i))

    t$tip.label[grepl("cluster",t$tip.label)] = sub('_[^_]*$', '', t$tip.label[grepl("cluster",t$tip.label)])
    
    t = keep.tip(t,mtdt[mtdt$Strain %in% t$tip.label & mtdt$Date > year+0.5,]$Strain)
    
    maxdate = max(mtdt[mtdt$Strain %in% t$tip.label,]$Date,na.rm=T)
    height = max(node.depth.edgelength(t),na.rm=T)
    rootlen = maxdate-height-year+0.5
    
    if (rootlen < 0){
      
      tipskept = t$tip.label
      nodeorder = order(node.depth.edgelength(t))
      norder = nodeorder[nodeorder>length(t$tip.label)]
      norder = norder[abs(maxdate-(node.depth.edgelength(t)[norder]-height))-year+0.5 > 0]
      branchidx = 1
      
      while (length(norder)>0){
        
        subtree = extract.clade(t,norder[1])
        rootlen = maxdate-max(node.depth.edgelength(subtree),na.rm=T)-year+0.5# - 2
        subtree$root.edge = rootlen
        trees[[idx]] = subtree
        idx = idx + 1
        norder = norder[2:length(norder)]
        norder = norder[!(norder %in% subtree$edge[,2])]
      }
    }
    
    else {
      t$root.edge = rootlen
      trees[[idx]] = t
      idx = idx + 1
    }
  }
  
  w = rtree(1)
  w$edge.length=0
  t = bind.tree(w,trees[[1]])
  if (length(trees) > 1){
    for (tr in 2:length(trees)){
      t = bind.tree(t,trees[[tr]])
    }
  }
  
  tokeep = sample(1:length(t$tip.label),500,replace=T)#seq(1,length(t$tip.label),as.integer(length(t$tip.label)/500))
  t = drop.tip(t,setdiff(1:length(t$tip.label),tokeep))
  t = drop.tip(t,mtdt[mtdt$Strain %in% t$tip.label & decimal_date(ymd(mtdt$Date)) > year + 1.4,]$Strain)
  
  tr = t
  y =ggtree(t) #+ geom_cladelab(node=mrca, label="1", align=TRUE,  offset = .1, textcolor=519+10, barcolor=519+10,fontsize=5/.pt)#,layout="circular")
  y$data$isroot = 'black'
  y$data$isroot[y$data$parent == min(y$data$parent)] = 'black'#white'
  y$data$isroot[y$data$parent == min(y$data$parent)] = 'black'#white'
  y$data$size = 1
  y$data$size[y$data$label == "t1"] = NA
  
  y$data$cluster = 0
  for (i in 1:length(y$data$cluster)){
    if (length(mtdt[mtdt$Strain == y$data$label[i],]$Cluster) == 0){
      next
    }
    y$data$cluster[i] = mtdt[mtdt$Strain == y$data$label[i],]$Cluster[1]
  }
  y$data$clustersave = y$data$cluster
  y$data$cluster = match(y$data$cluster,unique(y$data$cluster))
  
  fig2 = y + coord_flip() + scale_x_reverse() + aes(color=I(isroot),size=I(0.1))+ geom_tippoint(aes(fill=I(cluster+10),size=I(size)*1.2),pch=21,stroke=0.1) + theme_tree2() + 
    theme(axis.line=element_line(color='black',size = axissize),
          axis.text.y = element_text(size=7,color='black'),
          plot.title = element_text(size = 7,hjust=0.5)) + ggtitle(paste0(c("A/H3N2","A/H1N1pdm09","B/Victoria","B/Yamagata")[match(subtype,c("H3N2","H1N1","Vic","Yam"))],'\n',year,"/",year+1))
  #plot.title = element_text(size = 7,hjust=0.5)) + ggtitle(paste0(c("A/H3","A/H1pdm","B/Vic","B/Yam")[match(subtype,c("H3N2","H1N1","Vic","Yam"))],", ","'",substr(year,3,4),"/'",substr(year+1,3,4)))
  
  idx = 1
  colors = c()
  for (clustertoplot in cltoplot){
    mrca = getMRCA(t,mtdt[mtdt$Cluster==clustertoplot & mtdt$Strain%in%t$tip.label,]$Strain)
    cl = y$data$cluster[y$data$label%in%mtdt[mtdt$Cluster==clustertoplot & mtdt$Strain%in%t$tip.label,]$Strain][1]
    fig2 = fig2 + geom_cladelab(node=mrca, label=idx, align=TRUE,  offset = .05,offset.text=0.1, textcolor='black', barcolor=cl+10,fontsize=5/.pt)
    colors = c(colors,cl+10)
    idx = idx + 1
  }
  
  #ll = pl + xlim(0,2.1) + theme(plot.margin=unit(c(-3,-0.5,-3,-0.5),"cm")) 
  fig2_plot = fig2 #+ xlim(0,2.1) #+ theme(plot.margin=unit(c(-3,-0.5,-3,-0.5),"cm")) 
  xmax = max(fig2_plot$data$x)#floor(maxdate) + 0.3 - maxdate + max(fig2_plot$data$x)
  xmin = xmax-1.5
  maxdate = max(mtdt[mtdt$Strain%in%fig2_plot$data$label,]$Date)
  max_x = max(fig2_plot$data$x[!(is.na(fig2_plot$data$label))])
  
  fig2_plot = fig2_plot  + scale_x_reverse(limits=(rev(lims)),breaks=(max_x+year+1+c(-1,-0.5,0,0.5)-maxdate),labels=(year+1+c(-1,-0.5,0,0.5))) +
    theme(axis.line.y=element_line(size=0.2),axis.text.y=element_text(size=5),axis.ticks.y=element_line(size=0.2),axis.line.x = element_blank(),axis.ticks.x=element_blank(),axis.text.x = element_blank())
  
  return(list(fig2_plot+ theme(axis.text=element_text(size=5),plot.title=element_blank()),colors))
}


plotCluster_true <- function(idx,onsetdiff,fctr,clusters){
  
  df = cluster_df[cluster_df$Cluster==idx,]
  absent_states = state.name[!(state.name %in%df$State)]  
  for (i in absent_states){
    row = cluster_df[which(cluster_df$State==i)[1],]
    row$Percentage = 0
    row$Onset = NA
    df = rbind(df,row)
  }
  
  df = df[!(df$State%in%c("Hawaii","Alaska")),]
  df$Onset = df$Onset - min(df$Onset,na.rm=T) + onsetdiff
  df = df[!(is.na(df$Lat)),]
  df[,c("Lon","Lat")] = usmap_transform(df[,c("Lon","Lat")])[,c(3,4)]
  df$Percentage = as.numeric(df$Percentage)
  df$Percentage[df$Percentage<0.05] = 0
  
  
  df$Totpct = NA
  for (state in unique(df$State)){
    
    df[df$State==state,]$Totpct = sum(cluster_df[cluster_df$State==state & cluster_df$Cluster%in%clusters,]$Percentage,na.rm=T)
  }
  
  lims = c(0,15)
  c = plot_usmap(size=0.03,exclude=c('HI','AK'))  + geom_point(data = df[df$Totpct>0,], aes(x = Lon, y = Lat, size = sqrt(Totpct/20/pi)*fctr*2),bg='lightgrey', alpha = 0.2, pch=21,stroke=0.1)+
    geom_point(data = df[df$Percentage>0,], aes(x = Lon, y = Lat, size = sqrt(Percentage/20/pi)*fctr*2, bg = as.numeric(Onset)), alpha = 1, pch=21,stroke=0.1)  +  
    scale_fill_distiller(direction=-1,type='seq',palette='YlOrRd',name="Onset week",limits=lims,na.value='grey75') + scale_size_identity() 
  
  c = c + theme(legend.position = 'none')#+ theme(plot.margin=margin) 
  return(c)
}


plotCluster_sim = function(prop,timing,onsetdiff,fctr,clusters){
  
  df = data.frame(State=states,Percentage=prop,Onset=timing)
  df$Onset = df$Onset - min(df$Onset,na.rm=T) + onsetdiff + 1
  absent_states = states[!(states %in%df$State)]
  
  df = df[!(df$State%in%c("Hawaii","Alaska","New Hampshire")),]
  for (i in absent_states){
    df = rbind(df,c(0,i,0,0,0,NA,0,0,0))
  }
  df = df[!(df$State%in%c("Hawaii","Alaska","New Hampshire","District of Columbia")),]
  df = df[!(df$State%in%c("Hawaii","Alaska","New Hampshire","District Of Columbia")),]
  
  df$Percentage = as.numeric(df$Percentage)
  df$Percentage[df$Percentage<0.05] = 0
  df$Onset = as.numeric(df$Onset)
  lims = c(0,15)
  
  df$Totpct = NA
  for (state in unique(df$State)){
    
    df[df$State==state,]$Totpct = sum(cluster_df[cluster_df$State==state & cluster_df$Cluster%in%clusters,]$Percentage,na.rm=T)
  }
  
  df$Percentage = df$Percentage * df$Totpct
  
  
  df$Lon = NA
  df$Lat = NA
  for (stateidx in 1:50){
    state = state.name[stateidx]
    if (nrow(df[df$State == state,]) >= 1){
      df[df$State == state,]$Lon = state.center$x[stateidx]
      df[df$State == state,]$Lat = state.center$y[stateidx]
    }
  }
  
  df[,c("Lon","Lat")] = usmap_transform(df[,c("Lon","Lat")])[,c(3,4)]
  
  c = plot_usmap(size=0.03,exclude=c('HI','AK'))  + geom_point(data = df[df$Totpct>0,], aes(x = Lon, y = Lat, size = sqrt(Totpct/20/pi)*fctr*2),bg='lightgrey', alpha = 0.2, pch=21,stroke=0.1) +
    geom_point(data = df[df$Percentage>0,], aes(x = Lon, y = Lat, size = sqrt(Percentage/20/pi)*fctr*2, bg = as.numeric(Onset)), alpha = 1, pch=21,stroke=0.1)  +  
    scale_fill_distiller(direction=-1,type='seq',palette='YlOrRd',name="Onset week",limits=lims,na.value='grey75') + scale_size_identity() 
  
  c = c+ theme(legend.position = 'none')# + theme(plot.margin=margin) 
  return(c)
}



addIsolationsToTree <- function(treeplot, season, subtype){
  column = c("A..2009.H1N1.",'A..H3.','BYam','BVic')[match(subtype,c("H1N1","H3N2","Yam","Vic"))]
  isolates = read.csv("data/US_PH_labs.csv")
  isolates$Date = NA
  for (year in unique(isolates$YEAR)){
    isolates[isolates$YEAR==year,]$Date = year + isolates[isolates$YEAR==year,]$WEEK/max(isolates[isolates$YEAR==year,]$WEEK)
  }
  isolates$Season = NA
  for (i in seq(2015,2019)){
    isolates[isolates$Date>i+0.5 & isolates$Date<i+1.5,]$Season = i
  }
  islts = isolates[isolates$Season==season,]
  islts = islts[complete.cases(islts),]
  islts[,column] = cumsum(islts[,column])
  islts[,column] = (islts[,column] / max(islts[,column])) * max(treeplot$data$y)
  
  islts$D = islts$Date - islts$Season + 0.5
  islts$Val = islts[,column]
  treeplot =  treeplot +geom_line(data=islts,aes(x=D,y=Val),lty=1,inherit.aes=F,col='black',lwd=0.1) +  geom_ribbon(data=islts,aes(x=D,ymin=min(treeplot$data$y),ymax=Val),inherit.aes=F,fill='lightgrey',alpha=0.3)
  return(treeplot)
}
