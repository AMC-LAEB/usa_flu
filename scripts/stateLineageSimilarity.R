library(dplyr)
library(ggrepel)
library(reshape2)
library(vegan)
library(MASS)
library(ggpubr)
library(usmap)


getSingleDist <- function(season=NA, subtype=NA){
  
  mt = mtdt
  tokeep = unlist(lapply(states,function(x)length(which(tabulate(findInterval(mtdt[mtdt$State==x,]$Season,c(2014:2019,2022)))<10))<1))
  
  if (!(is.na(season))){
    mt = mtdt[mtdt$Season==season & mtdt$Subtype==subtype,]# & (cluster_df$Subtype%in%c("H1N1","H3N2")),]
    tokeep = unlist(lapply(states,function(x)nrow(mtdt[mtdt$State==x & mtdt$Season==season & mtdt$Subtype==subtype,])>=10))
    
  }
  
  file = file[file$State%in%states,]
  file=file[!(file$State%in%c("Alaska","Hawaii")),]
  
  dfs = vector(mode='list')
  idx = 1
  for (i in states[tokeep]){
    sampl = c()
    for (ssn in c(2014:2019,2022)){
      mt_i = mtdt[mtdt$State==i & mtdt$Season==ssn & !(mtdt$Cluster=="UNCLUSTERED"),]
      n = nrow(mt_i)
      if (n > 20){
        samp = sample(mt_i$Cluster,20,replace=F)
      } else {
        samp = mt_i$Cluster
      }
      sampl = c(sampl,samp)
    }
    tb = tabulate(as.numeric(sampl),nbins=max(as.numeric(mt$Cluster),na.rm=T))
    dfs[[idx]] = tb
    idx = idx + 1
  }
  
  if (length(dfs)<40){return(NA)}
  
  df = do.call(rbind,dfs)
  
  dists = vegdist((df),method='bray') 
  
  m = as.matrix(dists)
  
  colnames(m) = states[tokeep]
  rownames(m) = colnames(m)
  
  for (i in 1:ncol(m)){
    m[i,i] = 0
  }
  return(m)
}

getSimMat <- function(){
  
  s = replicate(50,getSingleDist())
  m = apply(s,c(1,2),mean)
  
  cmd = isoMDS(m)$points
  
  colnames(cmd) <- c("Dim.1", "Dim.2")
  cmd = as.data.frame(cmd)
  cmd$division = c(state.region[-c(2,11)],"Northeast")[match(rownames(cmd),states)]
  cmd$hhs = hhsregions[match(rownames(cmd),state.name)]
  cmd$state = rownames(cmd)
  return(list(m,cmd))
}

plotMDS <- function(cmd){
  cmd$state = rownames(cmd)
  cmd2 = cmd %>% group_by(division) %>% slice(chull(Dim.1,Dim.2))
  
  mds_plot = ggplot(cmd, aes(x = -1*Dim.1, y = -1*Dim.2,fill=factor(division),col=factor(division)),
                    label = factor(rownames(cmd)),
                    size = 4,
                    repel = TRUE) + scale_fill_brewer(palette='Set2') + scale_fill_brewer(palette='Set2',labels=c('Northeast','South','Midwest','West'))+
    thm + geom_text(aes(color=division,label=state.abb[match(state,state.name)]),size=4/.pt) + thm +
    theme(legend.position='top',legend.direction='horizontal') + geom_polygon(data = cmd2, alpha = 0.2,linewidth=0,aes(fill = factor(division),colour = factor(division))) + 
    scale_color_brewer(palette='Set2',guide='none')
  
  return(mds_plot)
}


getClosestAdjacent <- function(m){
  adjacency_mat = read.csv("data/adjacency_matrix.csv")
  for (i in 1:nrow(m)){m[i,i]=NA}
  closest_states = colnames(m)[apply(m,2,which.min)]
  adjacent_states = apply(adjacency_mat,1,function(x)state.abb[which(x==1)])
  is_adjacent = sapply(1:ncol(m),function(x)closest_states[x] %in% adjacent_states[[match(colnames(m)[x],state.abb)]])
  return(table(is_adjacent))
}

getHighestSimilarities <- function(m,N){
  sapply(1:N,function(x)print(which(m==sort(m,decreasing=F)[x*2],arr.ind=TRUE)))
}

get_distance_metrics = function(){
  distance = read.csv("/Users/simondejong/usa/centroid_dists.csv",h=F)
  return(list(distance))
}

plotDistMet <- function(m){
  
  distance_metrics = get_distance_metrics()
  dist_met_dfs = list()
  
  for (dist_met in 1:1){
    dm = distance_metrics[[dist_met]]
    idxes = match(colnames(m),state.name)
    dm = dm[idxes,idxes]
    
    rank_dist = (as.vector(unlist(apply(dm,2,rank))))
    rank_sim = (as.vector(apply(m,2,rank)))
    zer = which(unlist(dm)!=0)
    rank_dist = rank_dist[zer]
    rank_sim = rank_sim[zer]
    
    if (dist_met>1){
      rank_dist = 51 - rank_dist
    }
    
    nstate = dim(m)[1]
    plot_rank_df = data.frame(Metric=c("Centroid distance","Commuting","Air travel")[dist_met],
                              Hi=sapply(seq(1,nstate,1),function(x)quantile(rank_dist[rank_sim>=x & rank_sim<=x],c(0.75))),
                              Med=sapply(seq(1,nstate,1),function(x)quantile(rank_dist[rank_sim>=x & rank_sim<=x],c(0.50))),
                              Lo=sapply(seq(1,nstate,1),function(x)quantile(rank_dist[rank_sim>=x & rank_sim<=x],c(0.25))),
                              rank=1:nstate)
    
    
    dist_met_dfs[[dist_met]] = plot_rank_df
  }
  
  plot_rank_df = do.call(rbind,dist_met_dfs)
  
  thm = theme(axis.text = element_text(color="black",size=5),
              panel.spacing.x = unit(1, "mm"),
              axis.line=element_line(color='black',size = 0.2),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(),
              panel.background = element_blank(),
              legend.position='none',
              axis.title= element_text(size=5),
              strip.background=element_rect(colour="white",fill="white"),
              strip.text.x = element_text(size = 5),
              legend.text = element_text(size = 5),
              legend.title=element_blank(),
              legend.spacing.x = unit(0.05,"cm"),
              legend.key=element_rect(fill="white"))
  
  
  distmet_plot = ggplot(plot_rank_df[plot_rank_df$Metric=="Centroid distance",]) + geom_errorbar(width=1,aes(x=rank,ymin=Lo,ymax=Hi),position=position_dodge(width=1.5),linewidth=0.2) +
    thm + xlab("Similarity rank") + ylab("Distance rank") + thm + 
    stat_smooth(aes(x=rank,y=Med),col='black',se=T,linewidth=0.5)
  return(distmet_plot)
}

getSeasonPValues <- function(){
  
  distance = read.csv("/Users/simondejong/usa/centroid_dists.csv",h=F)
  
  mantel_test = list()
  idx = 1
  for (i in unique(cluster_df_country$Season)){
    season = i
    if (season==2021){next}
    for (j in unique(cluster_df_country[cluster_df_country$Season==i,]$Subtype)){
      subtype = j
      mat = getSingleDist(season,subtype)
      if ("District Of Columbia" %in% colnames(mat)){
        mat = mat[-c(match("District Of Columbia",colnames(mat))),-c(match("District Of Columbia",colnames(mat)))]
      }
      if (is.na(mat)){next}
      mantel_t = mantel(mat,log(distance[match(colnames(mat),state.name),match(colnames(mat),state.name)]),method='spearman')
      mantel_t$subtype = subtype
      mantel_t$season = season
      mantel_test[[idx]] = mantel_t
      idx = idx + 1
    }
  }
  
  p_values = sapply(mantel_test,function(x)(x$signif))
  return(p_values)
}

