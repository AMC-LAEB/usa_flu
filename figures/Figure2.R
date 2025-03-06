source("stateLineageSimilarity.R")
library(ggdendro)
library(RColorBrewer)

thm = theme(axis.text = element_text(color="black",size=6),
            panel.spacing.x = unit(1, "mm"),
            axis.line=element_line(color='black',size = 0.1),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            legend.position='none',
            axis.title= element_text(size=7),
            strip.background=element_rect(colour="black",fill="white"),
            strip.text.x = element_text(size = 7),
            legend.text=element_text(size=7),
            legend.title =element_text(size=7),
            legend.key=element_rect(fill="white"),
            legend.key.width = unit(.5, "line"),
            legend.spacing.y = unit(.001, 'cm'),
            legend.box.spacing = unit(0, "pt"),
            legend.margin=margin(1,1,1,1))

sim_mat = getSimMat()

distmet = plotDistMet(sim_mat[[1]]) + thm
mds = plotMDS(sim_mat[[2]]) + xlab("MDS Dim. 1") + ylab("MDS Dim. 2")+ 
  theme(legend.key.height=unit(5,'pt'),
        legend.key.width=unit(5,'pt'),
        legend.spacing.x = unit(.1, 'cm'),
        legend.text = element_text(size=6))

m = sim_mat[[1]]
colnames(m) = state.abb[match(colnames(m),state.name)]
rownames(m) = colnames(m)

for (i in 1:ncol(m)){m[i,i]=NA}

getHighestSimilarities(m,10)
getClosestAdjacent(m)
getSeasonPValues()

getDendrogram <- function(){
  dend = ggdendrogram(hclust(as.dist(m)),size=0.1) + theme(axis.text.y = element_text(size=7,hjust=0.5,color='black'),axis.line.y = element_line(size=0.1),axis.ticks.y = element_line(size=0.2)) + 
    coord_cartesian(ylim=c(0.65,0.91)) + scale_y_continuous(breaks=seq(0.65,0.85,0.1)) 
  dend$layers[[2]]$aes_params$size <- .2
  f <- function(pal) brewer.pal(brewer.pal.info[pal, "maxcolors"], pal)
  cols <- f("Set2")[c(2,4,1,3)]
  
  dend = dend +
    theme(axis.text.x = element_text(size=6,color=cols[match(state.region[match(ggplot_build(dend)$layout$panel_params[[1]]$x$get_labels(),state.abb)],unique(state.region))]),
          axis.text.y=element_text(color='black'))
  
  return(dend)
}

dend = getDendrogram()

getMaps <- function(){
  
  mrg = theme(plot.margin = margin(c(-1,-3,-3,-4),'cm'))
  maplist = list()
  cluster_to_plot = c(747,225,17,414,226,422,984,498,306,1180,309,274)
  cdc = cluster_df_country
  for (cl in 1:length(cluster_to_plot)){
    ssn = cdc[cdc$Cluster==cluster_to_plot[cl],]$Season
    sbt = cdc[cdc$Cluster==cluster_to_plot[cl],]$Subtype
    others = cdc[cdc$Season==ssn & cdc$Subtype==sbt & cdc$Cluster!="UNCLUSTERED",]$Cluster
    maplist[[cl]] = plotCluster_true(cluster_to_plot[cl],0,13,c(),c(0,12)) +mrg
  }
  return(maplist)
}

maplist = getMaps()

leg = as_ggplot(get_legend(plotCluster_true(cluster_to_plot[cl],0,20,others) + 
                             theme(legend.direction='vertical',legend.position='top',
                                   legend.text=element_text(size=6),legend.title=element_text(size=6),
                                   legend.key.height=unit(12,'pt'),legend.key.width=unit(5,'pt')) +
                             guides(fill = guide_colourbar(title='Week',title.position="top", title.hjust = 0.5))))

leg2 = as_ggplot(get_legend(mds))
ggarrange(
  ggarrange(leg2,dend+theme(plot.margin=margin(0,-4,-6,1)),ggarrange(mds+theme(legend.position = 'none'),distmet,nrow=1,widths=c(0.6,0.4)),ncol=1,heights=c(0.05,0.3,0.5)),
  ggarrange(
    ggarrange(plotlist=maplist,
              ncol=3,nrow=4),
    ggarrange(ggplot()+theme_void(),leg,ggplot()+theme_void(),ncol=1,heights=c(0.4,0.3,0.3)),
    nrow=1,widths=c(0.9,.13)),
  nrow=1,widths=c(0.26,0.28))

#ggsave("Figure_2_resubmit.pdf",width=180,units="mm",height=75)
ggsave("Figure_2_resubmit.pdf",width=2080,units="px",height=850,dpi=320)

plotMDSSingleSeason()
plotSimilarityCor()