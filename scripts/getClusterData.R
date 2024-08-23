library(stringr)
library(sleekts)
library(ape)

setwd("/Users/simondejong/US_phylo")


clusterpath = "time_0.083_pct_0.1_2"


get_incidence_data <- function(plot=F){
  
  # Read ILI data
  ili = read.csv("~/united_states/ILINet.csv")
  ili = ili[ili$YEAR >= 2014,]
  ili[ili=="X"] = NA
  ili = ili[,c("REGION","YEAR","WEEK","X.UNWEIGHTED.ILI")]
  
  ili$YEARWEEK = as.numeric(ili$YEAR) + as.numeric(ili$WEEK) / 52
  ili[ili$YEAR%in%c(2014,2020),]$YEARWEEK = as.numeric(ili[ili$YEAR%in%c(2014,2020),]$YEAR) + as.numeric(ili[ili$YEAR%in%c(2014,2020),]$WEEK) / 53
  
  
  colnames(ili) = c("Region","Year","Week","Incidence","YearWeek")
  ili = ili[,c("Region","Incidence","YearWeek","Year","Week")]
  iliYearWeek = as.numeric(ili$YearWeek)
  ili = ili[ili$Region %in% state.name,]
  ili$Incidence = as.numeric(ili$Incidence) * 0.01
  ili[is.na(ili)] = 0
  
  #Read confirmed cases data
  
  #First, for year > 2015
  epidata = read.csv("~/united_states/WHO_NREVSS_Clinical_Labs.csv")
  epidata[epidata=="X"] = NA
  epidata = epidata[,c("REGION","YEAR","WEEK","PERCENT.A","PERCENT.B")]
  epidata$YEARWEEK = as.numeric(epidata$YEAR) + as.numeric(epidata$WEEK)/52
  epidata[epidata$YEAR%in%c(2014,2020),]$YEARWEEK = as.numeric(epidata[epidata$YEAR%in%c(2014,2020),]$YEAR) + as.numeric(epidata[epidata$YEAR%in%c(2014,2020),]$WEEK)/53
  
  colnames(epidata) = c("Region","Year","Week","A","B","YearWeek")
  epidata[is.na(epidata)] = 0
  epidata[,c(2,3)] = as.numeric(as.matrix(epidata[,c(2,3)]))
  epidata$A = as.numeric(epidata$A)
  epidata$B = as.numeric(epidata$B)
  epidata$A = epidata$A * 0.01
  epidata$B = epidata$B * 0.01
  
  #For year < 2015
  epidatapre2015 = read.csv("~/united_states/WHO_NREVSS_Combined_prior_to_2015_16.csv")
  epidatapre2015[epidatapre2015=="X"] = NA
  epidatapre2015[,4:ncol(epidatapre2015)] = as.numeric(as.matrix(epidatapre2015[,4:ncol(epidatapre2015)]))
  
  epidatapre2015$YearWeek = as.numeric(epidatapre2015$YEAR) + as.numeric(epidatapre2015$WEEK)/52
  epidatapre2015[epidatapre2015$YEAR%in%c(2014,2020),]$YearWeek = as.numeric(epidatapre2015[epidatapre2015$YEAR%in%c(2014,2020),]$YEAR) + as.numeric(epidatapre2015[epidatapre2015$YEAR%in%c(2014,2020),]$WEEK)/53
  
  epidatapre2015$A = rowSums(epidatapre2015[,7:11],na.rm=T)/epidatapre2015$TOTAL.SPECIMENS 
  epidatapre2015$B = epidatapre2015$B/epidatapre2015$TOTAL.SPECIMENS 
  #epidatapre2015 = epidatapre2015[epidatapre2015$REGION != "District of Columbia",]
  epidatapre2015[is.na(epidatapre2015)] = 0
  epidatapre2015 = epidatapre2015[,c("REGION","YEAR","WEEK","A","B","YearWeek")]
  colnames(epidatapre2015) = c("Region","Year","Week","A","B","YearWeek")
  
  #Combine epidata datasets
  epidata = rbind(epidatapre2015,epidata)
  epidata$Season = 0
  for (i in 2014:2022){
    #epidata[epidata$YearWeek > i + 0.5 & epidata$YearWeek < i + 1.5,]$Season = i
    epidata[epidata$Week >= 26 & epidata$Year==i,]$Season = i
    epidata[epidata$Week < 26 & epidata$Year==i+1,]$Season = i
  }
  epidata = epidata[epidata$Season %in% 2014:2022,]
  epidata[epidata$YearWeek<2015.769 & epidata$YearWeek>2015.5,]$A = 0
  epidata[epidata$YearWeek<2015.769 & epidata$YearWeek>2015.5,]$B = 0
  epidata = epidata[epidata$Region %in% ili$Region,]
  
  #ili = ili[ili$YearWeek %in% epidata$YearWeek,]
  ili = ili[paste0(ili$Year,ili$Week) %in% paste0(epidata$Year,epidata$Week),]
  ili$Season = epidata$Season
  ili = ili[ili$Season!="2020",]
  
  ili$ILIPLUS_A = NA
  ili$ILIPLUS_B = NA
  for (season in unique(ili$Season)){
    states = unique(ili[ili$Season==season,]$Region)
    for (state in states){
      ili[ili$Season==season & ili$Region==state,]$ILIPLUS_A = sleek(as.numeric(ili[ili$Season==season & ili$Region==state,]$Incidence)*epidata[epidata$Season==season & epidata$Region==state,]$A)
      ili[ili$Season==season & ili$Region==state,]$ILIPLUS_B = sleek(as.numeric(ili[ili$Season==season & ili$Region==state,]$Incidence)*epidata[epidata$Season==season & epidata$Region==state,]$B)
      
    }
  }
  
  ili$EpiWeek = NA
  for (region in unique(ili$Region)){
    for (season in c(2015:2019,2021,2022)){
      firstidx = which(ili[ili$Region==region & ili$Season==season,]$Week>=26)[1]
      ili[ili$Region==region & ili$Season==season,][firstidx:(firstidx+51),]$EpiWeek = 1:52
    }
    season = 2014
    firstidx = which(ili[ili$Region==region & ili$Season==season,]$Week>=40)[1]
    ili[ili$Region==region & ili$Season==season,][firstidx:(firstidx+37),]$EpiWeek = 1:38
  }
  
  return(ili)
}

# Get incidebce data
ili = get_incidence_data()


getClusterData_subsamp <- function(clusterpath){
  
  # # Read sequence metadata files
  # # 
  # mtdt = setNames(cbind("H3N2",read.csv("metadata/H3N2_mt.csv")),c("Subtype","Strain","Date"))
  # mtdt = rbind(mtdt,setNames(cbind("H1N1",read.csv("metadata/H1N1_mt.csv")),c("Subtype","Strain","Date")))
  # mtdt = rbind(mtdt,setNames(cbind("Yam",read.csv("metadata/Yam_mt.csv")),c("Subtype","Strain","Date")))
  # mtdt = rbind(mtdt,setNames(cbind("Vic",read.csv("metadata/Vic_mt.csv")),c("Subtype","Strain","Date")))
  # 
  # mtdt$Season = NA
  # for (i in 2013:2023){
  #   mtdt[mtdt$Date > i + 0.5 & mtdt$Date < i + 1.5,]$Season = i
  # }
  # 
  # mtdt$Week = ceiling((mtdt$Date-floor(mtdt$Date))*52)
  # mtdt$Week[mtdt$Week==0] = 1
  # 
  # mtdt$State = gsub("_"," ",str_match(mtdt$Strain, "/\\s*(.*?)\\s*/")[,2])
  # mtdt[mtdt$State=="PENNSYLVANIA",]$State = "Pennsylvania"
  # 
  # mtdt$EpiWeek = ili$EpiWeek[match(paste0(mtdt$Season,mtdt$Week),paste0(ili$Season,ili$Week))]
  # 
  # mtdt_all = mtdt
  # 
  # nrow(mtdt[!(mtdt$State%in%c("Alaska","Hawaii")),])
  
  mtdt = read.csv("subsampled_mtdt.tsv")
  
  # Get cluster for each sequence
  
  mtdt$Cluster = "UNCLUSTERED"
  mtdt$File = NA
  mtdt$TMRCA = NA
  
  clusters = list.files(paste0("clustering_output/",clusterpath))
  trees_paths = clusters[grepl(".tre",clusters)]
  clusters = clusters[grepl(".txt",clusters)]
  
  idx = 1
  for (clidx in 1:length(clusters)){
    season = str_match(clusters[clidx], '([^_]+)(?:_[^_]+){1}$')[2]
    if (season==2021){season=2022}
    file = read.table(paste0("clustering_output/",clusterpath,"/",clusters[clidx]),h=T)
    tree = read.nexus(paste0("clustering_output/",clusterpath,"/",trees_paths[clidx]))
    tree$tip.label[grepl("cluster",tree$tip.label)] = sub('_[^_]*$', '', tree$tip.label[grepl("cluster",tree$tip.label)])
    clusters_cl = unique(file$CLUSTER)
    for (clidx_2 in 1:length(clusters_cl)){
      if (nrow(mtdt[mtdt$Season==season&mtdt$Strain%in%file[file$CLUSTER==clusters_cl[clidx_2],]$TAXA,])==0){next}
      
      ##if (season<2020){
      #  mtdt[mtdt$Season==season&mtdt$Strain%in%file[file$CLUSTER==clusters_cl[clidx_2],]$TAXA,]$Cluster = idx
      #  mtdt[mtdt$Season==season&mtdt$Strain%in%file[file$CLUSTER==clusters_cl[clidx_2],]$TAXA,]$File = clusters[clidx]
      #} else {
        mtdt[mtdt$Season>=season&mtdt$Strain%in%file[file$CLUSTER==clusters_cl[clidx_2],]$TAXA,]$Cluster = idx
        mtdt[mtdt$Season>=season&mtdt$Strain%in%file[file$CLUSTER==clusters_cl[clidx_2],]$TAXA,]$File = clusters[clidx]
     # }
        
      mrcA = getMRCA(tree,file[file$CLUSTER==clusters_cl[clidx_2],]$TAXA)
      maxdate = max(mtdt[mtdt$Strain%in%file[file$CLUSTER==clusters_cl[clidx_2],]$TAXA,]$Date)
      depths = node.depth.edgelength(tree)
      tmrca = maxdate - (max(depths) - depths[mrcA])
      mtdt[mtdt$Season==season & mtdt$Strain%in%file[file$CLUSTER==clusters_cl[clidx_2],]$TAXA,]$TMRCA = tmrca
      
      idx = idx + 1
    }
  }
  
  # Add season and week

  
  mtdt$Week = ceiling((mtdt$Date-floor(mtdt$Date))*52)
  mtdt$Week[mtdt$Week==0] = 1
  
  mtdt$State = gsub("_"," ",str_match(mtdt$Strain, "/\\s*(.*?)\\s*/")[,2])
  mtdt[mtdt$State=="PENNSYLVANIA",]$State = "Pennsylvania"
  
  mtdt$EpiWeek = ili$EpiWeek[match(paste0(mtdt$Season,mtdt$Week),paste0(ili$Season,ili$Week))]
  
  
  statepops = read.csv("/Users/simondejong/usa/statepops.csv")
  statepops = statepops[statepops$year=="2012",]
  statepops = statepops[statepops$ages=="total",]
  statepops = statepops[statepops$state.region %in% state.abb,]
  pops = statepops$population
  
  
  # Get weekly incidence by cluster
  
  subtype_dfs = list()
  subtype_dfs_trajec = list()
  
  for (type in 1:2){
    
    if (type==1){
      subtypes = c("H3N2","H1N1")
    } else {
      subtypes = c("Yam","Vic")
    }
    by_week_state_year = vector(mode="list")
    by_week_state_year_trajec = vector(mode="list")
    
    
    for (state_idx in 1:length(unique(mtdt$State))){
      
      by_state = vector(mode="list")
      by_state_trajec = vector(mode="list")
      
      state = unique(mtdt$State)[state_idx]
      popsize = pops[match(state,state.name)]
      
      for (season_idx in 1:8){
        
        season = c(2014:2019,2021,2022)[season_idx]
        
        mt = mtdt[mtdt$Season == season & mtdt$Subtype %in% subtypes & mtdt$State==state,]
        
        clusters = unique(mt$Cluster)
        
        seqs_per_week = 
          do.call(rbind,
                  lapply(clusters,
                         function(x)tabulate(mt[mt$Cluster==x,]$EpiWeek,nbins=52)))
        
        if (is.null(seqs_per_week)){next}
        
        if (season==2014){
          seqs_per_week_smoothed = do.call(rbind,lapply(1:nrow(seqs_per_week),function(y)colSums(do.call(rbind,lapply(1:38,function(x)dnorm(1:38,x,6)*seqs_per_week[y,][x])),na.rm=T)))
        } else {
          seqs_per_week_smoothed = do.call(rbind,lapply(1:nrow(seqs_per_week),function(y)colSums(do.call(rbind,lapply(1:52,function(x)dnorm(1:52,x,6)*seqs_per_week[y,][x])),na.rm=T)))
        }
        
        for (wk in 1:ncol(seqs_per_week_smoothed)){
          seqs_per_week_smoothed[,wk] = seqs_per_week_smoothed[,wk]/sum(seqs_per_week_smoothed[,wk],na.rm=T)
        }
        if (type == 1){
          incidence = ili[ili$Season==season & ili$Region==state & !(is.na(ili$EpiWeek)),]$ILIPLUS_A
        } else {
          incidence = ili[ili$Season==season & ili$Region==state & !(is.na(ili$EpiWeek)),]$ILIPLUS_B
        }
        
        cluster_incs_fit = t(t(seqs_per_week_smoothed)*incidence)
        
        df_trajec = data.frame(cluster_incs_fit)
        df_trajec = cbind(season,df_trajec)
        df_trajec = cbind(clusters,df_trajec)
        df_trajec = cbind(state,df_trajec)
        df_trajec = cbind(type,df_trajec)
        colnames(df_trajec)[1:4] = c("Type","State","Cluster","Season")
        
        seqcounts = as.vector(table(mt$Cluster)[match(unique(mt$Cluster),names(table(mt$Cluster)))])
        
        total_inc = sum(cluster_incs_fit)
        
        onsets = apply(cluster_incs_fit,1,function(x)which(cumsum(x)>0.05*total_inc)[1])
        
        percentages = seqcounts/sum(seqcounts)#apply(cluster_incs_fit,1,sum)/total_inc
        
        percentages_absolute = apply(cluster_incs_fit,1,sum)/total_inc
        
        subtype = mt$Subtype[match(clusters,mt$Cluster)]
        
        df = data.frame(Subtype = subtype, State=state,Season=season,Cluster=clusters,Percentage=percentages,Onset=onsets,Seqcount=seqcounts)
        
        for (subtype in subtypes){
          
          df[df$Subtype==subtype,]$Percentage = df[df$Subtype==subtype,]$Percentage / sum(df[df$Subtype==subtype,]$Percentage)
          
        }
        
        by_state[[season_idx]] = df
      }
      
      by_week_state_year[[state_idx]] = by_state
      
    }
    
    full_df = do.call(rbind,lapply(1:length(by_week_state_year),function(x)do.call(rbind,by_week_state_year[[x]])))
    subtype_dfs[[type]] = full_df
  }
  
  # Combine all outputs
  cluster_df = do.call(rbind,subtype_dfs)
  
  cluster_df$Lon = NA
  cluster_df$Lat = NA
  
  for (stateidx in 1:50){
    state = state.name[stateidx]
    if (nrow(cluster_df[cluster_df$State == state,]) > 1){
      cluster_df[cluster_df$State == state,]$Lon = state.center$x[stateidx]
      cluster_df[cluster_df$State == state,]$Lat = state.center$y[stateidx]
    }
  }
  
  # Get state-specific relative onsets
  
  cluster_df$RelativeOnset_State = NA
  for (subtype in unique(cluster_df$Subtype)){
    
    for (season in unique(cluster_df$Season)){
      
      for (state in unique(cluster_df$State)){
        cluster_df[cluster_df$Subtype==subtype & cluster_df$Season==season & cluster_df$State==state,]$RelativeOnset_State = 
          cluster_df[cluster_df$Subtype==subtype & cluster_df$Season==season & cluster_df$State==state,]$Onset - 
          min(cluster_df[cluster_df$Subtype==subtype & cluster_df$Season==season & cluster_df$State==state,]$Onset,na.rm=T)
        
      }
    }
  }
  
  # Get country_specific relative onsets
  
  cluster_df_country = data.frame(Subtype=NA,Season=NA,Cluster=NA,Percentage=NA,RelativeOnset=NA,RelativeOnset_Country=NA,Seqcount=NA)
  
  for (subtype in unique(cluster_df$Subtype)){
    
    for (season in unique(cluster_df$Season)){
      
      df = cluster_df[cluster_df$Subtype==subtype & cluster_df$Season==season,]
      
      for (cl in unique(df$Cluster)){
        
        ons = min(df[df$Cluster==cl & !(df$State %in% c("Alaska","Hawaii")),]$Onset,na.rm=T)
        onset = min(df[df$Cluster==cl,]$Onset,na.rm=T) - min(df$Onset,na.rm=T)
        sum_pct = sum(df[df$Cluster==cl,]$Percentage,na.rm=T)
        sum_seqcount = sum(df[df$Cluster==cl,]$Seqcount,na.rm=T)
        cluster_df_country = rbind(cluster_df_country,c(subtype,season,cl,sum_pct,ons,onset,sum_seqcount))
        
      }
    }
  }
  
  cluster_df_country = cluster_df_country[!(is.na(cluster_df_country[,1])),]
  cluster_df_country$Clustered = 1
  cluster_df_country[cluster_df_country$Cluster=="UNCLUSTERED",]$Clustered = 0
  cluster_df_country$Seqcount = as.numeric(cluster_df_country$Seqcount)
  cluster_df_country$Percentage = as.numeric(cluster_df_country$Percentage)*2
  
  cluster_df$Subtype = factor(cluster_df$Subtype,levels=c('H3N2','H1N1','Yam','Vic'),labels=c("A/H3N2","A/H1N1pdm09","B/Yam","B/Vic"))
  
  cluster_df_country$Percentage_Season = NA
  cluster_df_country$Substantial_Subtype = 0
  for (season in unique(cluster_df_country$Season)){
    for (subtype in unique(cluster_df_country$Subtype)){
      if (nrow(cluster_df_country[cluster_df_country$Season==season & cluster_df_country$Subtype==subtype,])==0){next}
      prop = epidemic_compositions[epidemic_compositions$Season==season & epidemic_compositions$Subtype==subtype,]$Prop
      cluster_df_country[cluster_df_country$Season==season & cluster_df_country$Subtype==subtype,]$Percentage_Season = 
        cluster_df_country[cluster_df_country$Season==season & cluster_df_country$Subtype==subtype,]$Percentage * prop
      if (prop > 0.1){
        cluster_df_country[cluster_df_country$Season==season & cluster_df_country$Subtype==subtype,]$Substantial_Subtype = 1
      }
    }
  }
  
  cluster_df_country$Subtype=factor(cluster_df_country$Subtype,levels=c("H3N2","H1N1","Yam","Vic"))
  
  cluster_df_country = cluster_df_country[!(is.na(cluster_df_country[,1])),]
  cluster_df_country$TMRCA = NA
  for (i in 1:nrow(cluster_df_country)){
    cl = cluster_df_country[i,]$Cluster
    if (cl == "UNCLUSTERED"){next}
    tmrca = mtdt[mtdt$Cluster==cl,]$TMRCA[1]
    cluster_df_country[i,]$TMRCA = tmrca
  }
  
  
  getOnsetCountry <- function(subtype,season){
    column = c("A..2009.H1N1.",'A..H3.','BYam','BVic')[match(subtype,c("H1N1","H3N2","Yam","Vic"))]
    if (season > 2014){
      ssns = c(2015:2019,2022)
      isolates = read.csv("data/US_PH_labs.csv")
    } else {
      if (subtype %in% c("Yam","Vic")){return(NA)}
      ssns = c(2014)
      isolates = read.csv("data/National_pre2015.csv")
    }
    isolates$Date = NA
    for (year in unique(isolates$YEAR)){
      isolates[isolates$YEAR==year,]$Date = year + isolates[isolates$YEAR==year,]$WEEK/max(isolates[isolates$YEAR==year,]$WEEK)
    }
    isolates$Season = NA
    for (i in c(ssns)){
      isolates[isolates$Date>i+0.5 & isolates$Date<i+1.5,]$Season = i
    }
    islts = isolates[isolates$Season==season,]
    islts = islts[complete.cases(islts),]
    
    onset_idx = which(cumsum(islts[,column])>0.05*sum(islts[,column]))[1]
    wk = islts[onset_idx,]$Date
    return(wk)
  } 
  
  cluster_df_country$FirstSamp = sapply(cluster_df_country$Cluster,function(x)min(mtdt[mtdt$Cluster==x,]$Date))
  cluster_df_country$EpidemicOnset = NA
  for (ssn in unique(cluster_df_country$Season)){
    for (subtype in unique(cluster_df_country$Subtype)){
      wk = getOnsetCountry(subtype,ssn)
      if (is.na(wk)){next}
      if (nrow(cluster_df_country[cluster_df_country$Season==ssn & cluster_df_country$Subtype==subtype,])>0){
        cluster_df_country[cluster_df_country$Season==ssn & cluster_df_country$Subtype==subtype,]$EpidemicOnset = wk
      }
    }
  }
  
  cluster_df_country$RelativeOnset_Date = NA
  cluster_df_country[cluster_df_country$Season!=2014,]$RelativeOnset_Date = as.numeric(cluster_df_country[cluster_df_country$Season!=2014,]$Season)+0.5 + 
    as.numeric(cluster_df_country[cluster_df_country$Season!=2014,]$RelativeOnset)/52
  cluster_df_country[cluster_df_country$Season==2014,]$RelativeOnset_Date = as.numeric(cluster_df_country[cluster_df_country$Season==2014,]$Season)+(40/52) +
    as.numeric(cluster_df_country[cluster_df_country$Season==2014,]$RelativeOnset)/52
  
  
  return(list(mtdt,cluster_df,cluster_df_country))
  
}



dat = getClusterData_subsamp(clusterpath)
mtdt = dat[[1]]
cluster_df = dat[[2]]
cluster_df_country = dat[[3]]



