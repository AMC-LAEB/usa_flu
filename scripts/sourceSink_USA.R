allJumps <- function(folder, hhs=F){
  
  burnin = 100000
  
  clusters = list.files(folder)
  
  outlist = vector(mode='list',length=length(clusters))
  
  idx = 1
  
  for (cluster in clusters){
    
    if (hhs){
      if (!(grepl("hhs",cluster))){next}
    } else {
      if (grepl("hhs",cluster)){next}
    }
    if (!(grepl("log",cluster))){next}
    
    if (hhs){
      cluster_name = gsub('.hhs.geo.jumpHistory.log',"",cluster)
    } else {
      cluster_name = gsub('.geo.jumpHistory.log',"",cluster)
    }
    
    path = paste0("/Users/simondejong/US_phylo/",folder,"/",cluster)
    
    outpath = gsub(".geo.jumpHistory.log",".jumpTimes.txt",path)
    
    checkpath = file.exists(outpath)
    
    system(paste0("../collect_times ", burnin," < ",path," > ",outpath))
    
    f = read.csv(outpath,sep='\t')
    
    f$cluster = gsub("cluster_","",cluster_name)
    f$root = NA
    rts = getRootByState(f)
    rts = rts[!(is.na(rts$Cluster)),]
    for (i in unique(rts$Sample)){
      f[f$state==i,]$root = rts[rts$Sample==i,]$Root[1]
    }
    outlist[[idx]] = f
    idx = idx + 1
    
  }
  
  return(outlist)
}

getRootByState <- function(f){
  root_df = data.frame(Cluster=NA,Sample=NA,Root=NA)
  for (cluster in unique(f$cluster)){
    for (sample in unique(f[f$cluster==cluster,]$state)){
      df = f[f$cluster==cluster & f$state==sample,]
      first10 = table(df[order(df$time,decreasing=T)[1:10],]$from)
      rootloc = names(first10)[which.max(first10)]
      root_df = rbind(root_df,c(cluster,sample,rootloc))
    }
  }
  return(root_df)
}

getPropFromHHS <- function(state, output){
  
  clusters = unique(cluster_df[cluster_df$State==state,]$Cluster)
  clusters = clusters[!(clusters=="UNCLUSTERED")]
  all_tb = rep(0,10)
  for (cluster in clusters){
    season = cluster_df[cluster_df$Cluster==cluster & cluster_df$State==state,]$Season
    if (nrow(mtdt[mtdt$State==state & mtdt$Season==season,])<10){next}
    prop = nrow(cluster_df[cluster_df$Cluster==cluster & cluster_df$State==state,])/nrow(cluster_df[cluster_df$Season==season & cluster_df$State==state,])
    tb = tabulate(as.numeric(gsub("HHS","",output[output$cluster==cluster,]$root)),nbins=10)
    if (sum(tb)>0){
      tb = tb/sum(tb)
    }
    tb = tb * prop
    all_tb = all_tb + tb
  }
  names(all_tb) = 1:10
  return(all_tb/sum(all_tb))
}


table(cluster_df_country[cluster_df_country$Cluster!="UNCLUSTERED",]$Percentage_Season>=0.5)

all_sourcesink = do.call(rbind,allJumps("jumpHistories",F))
all_sourcesink_hhs = do.call(rbind,allJumps("jumpHistories",T))
all_sourcesink_hhs_subsamp = do.call(rbind,allJumps("jumpHistories_subsamp",T))

state_props = sapply(states,function(x)getPropFromHHS(x,all_sourcesink_hhs))
state_props_subsamp = sapply(states,function(x)getPropFromHHS(x,all_sourcesink_hhs_subsamp))

hhsregions <- c(4,10,9,6,9,8,1,3,4,4,9,10,5,5,7,7,4,6,1,3,1,5,5,4,7,8,7,9,1,2,6,2,4,8,5,6,10,3,1,4,8,4,6,8,1,3,10,3,5,8)
hhsreg = c(hhsregions[-c(2,11)],3)

mlt = melt(state_props_subsamp)
colnames(mlt) = c("source_hhs","state","prop")
mlt$to_hhs = hhsreg[match(mlt$state,states)]
mlt$st = c(state.abb[-c(2,11)],"DC")[match(mlt$state,states)]
ggplot(mlt,aes(x=st,y=prop)) + geom_point(aes(fill=factor(to_hhs)),pch=21,stroke=0.1,cex=2) + 
  scale_fill_brewer(palette='Set3',name="HHS region") + facet_wrap(.~factor(source_hhs),scales='free',ncol=2,nrow=5) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + theme(axis.text=element_text(size=5)) + xlab("State") + 
  ylab("Proportion") + theme(axis.title=element_text(size=6))
ggsave("Figure_SS_subsamp.pdf",width=180,units="mm",height=200)

mlt = melt(state_props)
colnames(mlt) = c("source_hhs","state","prop")
mlt$to_hhs = hhsreg[match(mlt$state,states)]
mlt$st = c(state.abb[-c(2,11)],"DC")[match(mlt$state,states)]
ggplot(mlt,aes(x=st,y=prop)) + geom_point(aes(fill=factor(to_hhs)),pch=21,stroke=0.1,cex=2) + 
  scale_fill_brewer(palette='Set3',name="HHS region") + facet_wrap(.~factor(source_hhs),scales='free',ncol=2,nrow=5) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + theme(axis.text=element_text(size=5)) + xlab("State") + 
  ylab("Proportion") + theme(axis.title=element_text(size=6))
ggsave("Figure_SS_nonsubsamp.pdf",width=180,units="mm",height=200)


adjacency_mat = mob[[5]]

m1 = as.matrix(dist(t(state_props_subsamp[,1:48])))
m2 = distance[match(colnames(state_props_subsamp[,1:48]),state.name),match(colnames(state_props_subsamp[,1:48]),state.name)]
for (i in 1:48){
  for (j in 1:48){
    if (hhsreg[i]==hhsreg[j]){
      m1[i,j] = NA
      m2[i,j] = NA
    }
  }
}

mantel(m1,m2,na.rm=T)


getLocalPct <- function(state,df){
  hhs_local = hhsreg[match(state,states)]
  return(sum(df[hhs_local,match(state,colnames(df))],na.rm=T))
}

local_prop = sapply(states,function(x)getLocalPct(x,state_props))
quantile(local_prop)

local_prop_subsamp = sapply(states,function(x)getLocalPct(x,state_props_subsamp))
quantile(local_prop_subsamp)

cor.test(state_props,state_props_subsamp,method='spearman')
quantile(c(local_prop,local_prop_subsamp))

quantile(sapply(1:length(local_prop),function(x)mean(c(local_prop[x],local_prop_subsamp[x]))))
