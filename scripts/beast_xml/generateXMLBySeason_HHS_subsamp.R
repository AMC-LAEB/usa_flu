# all_seqs_to_remove = c()
# for (season in c(2014:2019,2021,2022)){
#   for (subtype in unique(mtdt_all$Subtype)){
#     if (nrow(mtdt_all[mtdt_all$Season==season & mtdt_all$Subtype==subtype,])<1){next}
#     maxsampcount = getMaxSampCounts(season,subtype,mtdt_all,0)
#     for (state in 1:length(states)){
#       if (!(is.na(maxsampcount[state]))){
#         if (length(mtdt[mtdt$Season==season & mtdt$Subtype == subtype & mtdt$State==states[state],]$Strain)>maxsampcount[state]){
#           seqs_to_keep = sample(mtdt[mtdt$Season==season & mtdt$Subtype == subtype & mtdt$State==states[state],]$Strain,maxsampcount[state],replace=F)
#           seqs_to_remove = mtdt[mtdt$Season==season & mtdt$Subtype == subtype & mtdt$State==states[state] &!(mtdt$Strain%in%seqs_to_keep),]$Strain
#           all_seqs_to_remove = c(all_seqs_to_remove,seqs_to_remove)
#         }
#       }
#     }
#   }
# }
# write.csv(mtdt[!(mtdt$Strain%in%all_seqs_to_remove),],"/Users/simondejong/US_phylo/subsampled_mtdt_phylo.tsv",quote=F,row.names=F,sep='\t')
# 

# states = c(state.name[-c(2,11)],"District Of Columbia")
# hhsregions <- c(4,10,9,6,9,8,1,3,4,4,9,10,5,5,7,7,4,6,1,3,1,5,5,4,7,8,7,9,1,2,6,2,4,8,5,6,10,3,1,4,8,4,6,8,1,3,10,3,5,8)
# hhsreg = c(hhsregions[-c(2,11)],3)
# 
# all_seqs_to_remove = c()
# for (season in c(2014:2019,2021,2022)){
#   for (subtype in unique(mtdt_all$Subtype)){
#     if (nrow(mtdt[mtdt$Season==season & mtdt$Subtype==subtype,])<1){next}
#     
#     tab = tabulate(hhsreg[match(mtdt[mtdt$Season==season & mtdt$Subtype==subtype ,]$State,states)],nbin=10)
#     maxsampcount = as.integer(quantile(tab,c(0.25)))
#     for (hhsregion in 1:10){
#       
#       if (tab[hhsregion]>maxsampcount){
#         strains = mtdt[!(mtdt$State %in% c("Hawaii","Alaska")) & mtdt$Season==season & mtdt$Subtype == subtype & hhsreg[match(mtdt$State,states)] == hhsregion,]$Strain
#         seqs_to_keep = sample(strains[!(is.na(strains))],maxsampcount,replace=F)
#         seqs_to_remove = mtdt[mtdt$Season==season & mtdt$Subtype == subtype & hhsreg[match(mtdt$State,states)] == hhsregion &!(mtdt$Strain%in%seqs_to_keep),]$Strain
#         all_seqs_to_remove = c(all_seqs_to_remove,seqs_to_remove)
#         
#       }
#     }
#   }
# }
# write.csv(mtdt[!(mtdt$Strain%in%all_seqs_to_remove),],"/Users/simondejong/US_phylo/subsampled_mtdt_hhs_phylo.tsv",quote=F,row.names=F,sep='\t')



generateSingleStrain <- function(taxon,date,state){return(paste0(
  "<taxon id='",taxon,"'>","<date value='",date,
  "' direction='forwards' units='years' uncertainty='0.0'/>",
  "<attr name='state'>",state,"</attr>",
  "</taxon>"))}


generateTaxa <- function(name,taxa,dates,states){
  taxa_xml = do.call(paste0,lapply(1:length(taxa),function(x)generateSingleStrain(taxa[x],dates[x],states[x])))
  taxa_xml = paste0("\n<taxa id = '",name,"'>",taxa_xml)
  taxa_xml = paste0(taxa_xml,"</taxa>\n")
  
  return(taxa_xml)
}

generateStartingTree <- function(name,tree){
  xml_tree = paste0("<newick id ='",name,".startingTree' usingDates='true'>")
  xml_tree = paste0(xml_tree,write.tree(tree,file=""),"</newick>")
  return(xml_tree)
}


generateDataTree <- function(name,tree){
  xml_tree = paste0("<newick id ='",name,".dataTree' usingDates='false' usingHeights='true'>")
  xml_tree = paste0(xml_tree,write.tree(tree,file=""),"</newick>")
  return(xml_tree)
}


generateModelBlockFirst <- function(name,mtdt,maxdate,nsite,mod){
  cluster = gsub("cluster_","",name)
  
  firstseqdate = min(mtdt[mtdt$Cluster==cluster,]$Date)
  cutoff = maxdate - (firstseqdate-1/12)
  numgridpoints = ceiling((maxdate-(firstseqdate-1/12))*53)
  skygriddim = numgridpoints+1
  
  
  str = ('<constrainedTreeModel id = "tree1.treeModel">
    <tree idref="tree1.startingTree"/>
      <constraintsTree>
      <tree idref="tree1.dataTree"/>
        </constraintsTree>
        </constrainedTreeModel>
        <!-- Statistic for root height of the tree       -->
        <treeHeightStatistic id="tree1.treeModel.rootHeight">
          <treeModel idref="tree1.treeModel"/>
            </treeHeightStatistic>
            <!-- Statistic for sum of the branch lengths of the tree (tree length)       -->
            <treeLengthStatistic id="tree1.treeLength">
              <treeModel idref="tree1.treeModel"/>
                </treeLengthStatistic>
                <!-- Statistic for time of most recent common ancestor of tree               -->
                <tmrcaStatistic id="tree1.age(root)" absolute="true">
                  <treeModel idref="tree1.treeModel"/>
                    </tmrcaStatistic>
         
         
         	<thorneyTreeLikelihood id="tree1.treeLikelihood">
		<constrainedTreeModel idref="tree1.treeModel"/>
				<strictClockBranchLengthLikelihood id="cluster_727.branchLengthLikelihood" scale="14402.0">
			<parameter idref="clock.rate"/>
		</strictClockBranchLengthLikelihood>

		
		<constrainedBranchLengthProvider scale="13427.0">
			<constrainedTreeModel idref="tree1.treeModel"/>
			<dataTree>
				<tree idref="tree1.dataTree"/>
			</dataTree>
		</constrainedBranchLengthProvider>
	</thorneyTreeLikelihood>
	
	
'
  )
  
  if (mod=="skygrid"){
    str = paste0(str,'
                          <gmrfSkyGridLikelihood id="tree1.skygrid">
  <populationSizes>
  <!-- skygrid.logPopSize is in log units unlike other popSize                 -->
  <parameter id="tree1.skygrid.logPopSize" dimension="skygriddim" value="1.0"/>
  </populationSizes>
  <precisionParameter>
  <parameter id="skygrid.precision" value="0.1" lower="0.0"/>
  </precisionParameter>
  <numGridPoints>
  <parameter id="tree1.skygrid.numGridPoints" value="ngridpoint"/>
  </numGridPoints>
  <cutOff>
  <parameter id="tree1.skygrid.cutOff" value="cutoff"/>
  </cutOff>

  <intervals>
  <bigFastTreeIntervals>
  <treeModel idref="tree1.treeModel"/>
  </bigFastTreeIntervals>
  </intervals>

  </gmrfSkyGridLikelihood>
                 ')
  } else {
    
    str = paste0(str,'
                 <exponentialGrowth id="tree1.exponential" units="years">
		<populationSize>
			<parameter id="tree1.exponential.popSize" value="1.0" lower="0.0"/>
		</populationSize>
		<growthRate>
			<parameter id="tree1.exponential.growthRate" value="0.0"/>
		</growthRate>
	</exponentialGrowth>
	
	
	 <coalescentLikelihood id="tree1.coalescent">
           	<model>
		<exponentialGrowth idref="tree1.exponential"/>
	          </model>
           
           <intervals>
           <bigFastTreeIntervals>
           <treeModel idref="tree1.treeModel"/>
           </bigFastTreeIntervals>
           </intervals>
           
           </coalescentLikelihood>
')
  }
  
  
  str = gsub("tree1",name,str)
  str = gsub("cutoff",cutoff,str)
  str = gsub("ngridpoint",numgridpoints,str)
  str = gsub("skygriddim",skygriddim,str)
  str = gsub("nsite",nsite,str)
  
  return(str)
}



generateModelBlockSecond <- function(name,mtdt,maxdate,nsite,mod){
  
  cluster = gsub("cluster_","",name)
  
  firstseqdate = min(mtdt[mtdt$Cluster==cluster,]$Date)
  cutoff = maxdate - (firstseqdate-1/12)
  numgridpoints = ceiling((maxdate-(firstseqdate-1/12))*53)
  skygriddim = numgridpoints+1
  
  str = ('
  
  	<constrainedTreeModel id = "tree2.treeModel">
		<tree idref="tree2.startingTree"/>
		<constraintsTree>
			<tree idref="tree2.dataTree"/>
		</constraintsTree>
	</constrainedTreeModel>
		<!-- Statistic for root height of the tree       -->
	<treeHeightStatistic id="tree2.treeModel.rootHeight">
	<treeModel idref="tree2.treeModel"/>
	</treeHeightStatistic>
			<!-- Statistic for sum of the branch lengths of the tree (tree length)       -->
	<treeLengthStatistic id="tree2.treeLength">
	<treeModel idref="tree2.treeModel"/>
	</treeLengthStatistic>
			<!-- Statistic for time of most recent common ancestor of tree               -->
	<tmrcaStatistic id="tree2.age(root)" absolute="true">
	<treeModel idref="tree2.treeModel"/>
	</tmrcaStatistic>
	
         
         	<thorneyTreeLikelihood id="tree2.treeLikelihood">
		<constrainedTreeModel idref="tree2.treeModel"/>
		
		<strictClockBranchLengthLikelihood id="tree2.branchLengthLikelihood" scale="nsite.0">
			<parameter idref="clock.rate"/>
		</strictClockBranchLengthLikelihood>
		
		<constrainedBranchLengthProvider scale="13427.0">
			<constrainedTreeModel idref="tree2.treeModel"/>
			<dataTree>
				<tree idref="tree2.dataTree"/>
			</dataTree>
		</constrainedBranchLengthProvider>
	</thorneyTreeLikelihood>

	
	
'
  )
  
  if (mod=="skygrid"){
    str = paste0(str,'
                          <gmrfSkyGridLikelihood id="tree2.skygrid">
  <populationSizes>
  <!-- skygrid.logPopSize is in log units unlike other popSize                 -->
  <parameter id="tree2.skygrid.logPopSize" dimension="skygriddim" value="1.0"/>
  </populationSizes>
  <precisionParameter>
  <parameter idref="skygrid.precision"/>
  </precisionParameter>
  <numGridPoints>
  <parameter id="tree2.skygrid.numGridPoints" value="ngridpoint"/>
  </numGridPoints>
  <cutOff>
  <parameter id="tree2.skygrid.cutOff" value="cutoff"/>
  </cutOff>

  <intervals>
  <bigFastTreeIntervals>
  <treeModel idref="tree2.treeModel"/>
  </bigFastTreeIntervals>
  </intervals>

  </gmrfSkyGridLikelihood>
                 ')
  } else {
    
    str = paste0(str,'
                 <exponentialGrowth id="tree2.exponential" units="years">
		<populationSize>
			<parameter id="tree2.exponential.popSize" value="1.0" lower="0.0"/>
		</populationSize>
		<growthRate>
			<parameter id="tree2.exponential.growthRate" value="0.0"/>
		</growthRate>
	</exponentialGrowth>
	
	
	 <coalescentLikelihood id="tree2.coalescent">
           	<model>
		<exponentialGrowth idref="tree2.exponential"/>
	          </model>
           
           <intervals>
           <bigFastTreeIntervals>
           <treeModel idref="tree2.treeModel"/>
           </bigFastTreeIntervals>
           </intervals>
           
           </coalescentLikelihood>
')
  }
  str = gsub("tree2",name,str)
  str = gsub("cutoff",cutoff,str)
  str = gsub("ngridpoint",numgridpoints,str)
  str = gsub("skygriddim",skygriddim,str)
  str = gsub("nsite",nsite,str)
  
  return(str)
  
}


generateOperatorSingle <- function(name, bifurc,mod){
  print(bifurc)
  str = '



		<nodeHeightOperator type="uniform" weight="574">
			<treeModel idref="tree1.treeModel"/>
		</nodeHeightOperator>
		
		<nodeHeightOperator type="scaleRoot" weight="114.8" scaleFactor="0.75" >
			<treeModel idref="tree1.treeModel"/>
		</nodeHeightOperator>
		


  '
  
  if (!(bifurc)){str = paste0(str,'
                              <uniformSubtreePruneRegraft weight="574">
			<constrainedTreeModel idref="tree1.treeModel"/>
		</uniformSubtreePruneRegraft>
				<narrowExchange weight="114.8">
			<constrainedTreeModel idref="tree1.treeModel"/>
		</narrowExchange>
		
		<wideExchange weight="114.8">
			<constrainedTreeModel idref="tree1.treeModel"/>
		</wideExchange>
		
		<wilsonBalding weight="114.8">
			<constrainedTreeModel idref="tree1.treeModel"/>
		</wilsonBalding>
')}
  
  if (mod=="skygrid"){
    str = paste0(str,'
      		<gmrfSkygridBlockUpdateOperator scaleFactor="1.1" weight="500">
			<gmrfSkygridLikelihood idref="tree1.skygrid"/>
		</gmrfSkygridBlockUpdateOperator>
		
				<scaleOperator scaleFactor="0.75" weight="100">
			<parameter idref="skygrid.precision"/>
		</scaleOperator>
		<randomWalkOperator windowSize="1.0" weight="100">
			<parameter idref="tree1.skygrid.logPopSize"/>
		</randomWalkOperator>           
                 ')
  } else {
    str = paste0(str,'
                 				<scaleOperator scaleFactor="0.75" weight="50">
			<parameter idref="tree1.exponential.popSize"/>
		</scaleOperator>
		<randomWalkOperator windowSize="1.0" weight="50">
			<parameter idref="tree1.exponential.growthRate"/>
		</randomWalkOperator>
                 ')
  }
  
  str = gsub("tree1",name,str)
  return(str)
}

generateAllOperators <- function(names,bifurcs,lengths){
  str = paste0("<operators id = 'operators' optimizationSchedule = 'default'>")
  str = paste0(str,'
                		<scaleOperator scaleFactor="0.75" weight="100">
			<parameter idref="clock.rate"/>
		</scaleOperator>
               ')
  for (i in 1:length(names)){
    str = paste0(str,generateOperatorSingle(names[i],bifurcs[[i]],c("skygrid","coal")[(lengths[[i]]<50)+1]))
  }
  str = paste0(str,"</operators>")
  return(str)
}

generateMCMC <- function(names,season,subtype,lengths){
  print(lengths)
  str = paste0('<mcmc id="mcmc" chainLength="500000000" autoOptimize="true">
                 <joint id="joint">
                 <prior id="prior">
               <gammaPrior id="clock.prior" shape="0.001" scale="1000" offset="0.0">
					<parameter idref="clock.rate"/>
				</gammaPrior>
               
               ')
  
  
  if (any(unlist(lengths)>=50)){
    str = paste0(str,'
    <gammaPrior id="skygrid.precision.prior" shape="0.001" scale="1000" offset="0.0">
      <parameter idref="skygrid.precision"/>
        </gammaPrior>
        '
    )
  }
  
  
  for (i in 1:length(names)){
    mod = c("skygrid","coal")[(lengths[[i]]<50)+1]
    if (mod=="skygrid"){
      str = paste0(str,gsub("tree1",names[[i]],'
				<gmrfSkyGridLikelihood idref="tree1.skygrid"/>	'))
    } else {
      
      str = paste0(str,gsub("tree1",names[[i]],'
    <oneOnXPrior>
					<parameter idref="tree1.exponential.popSize"/>
				</oneOnXPrior>
				<laplacePrior mean="0.0" scale="1.0">
					<parameter idref="tree1.exponential.growthRate"/>
				</laplacePrior>
				<coalescentLikelihood idref="tree1.coalescent"/>	'))
    }
  }
  
  
  
  str = paste0(str,'</prior>
			<likelihood id="likelihood">')
  for (i in names){
    str = paste0(str,gsub("tree1",i,'<thorneyTreeLikelihood idref="tree1.treeLikelihood"/>'))
  }
  str = paste0(str,'
               </likelihood>
		</joint>
		<operators idref="operators"/>

		<!-- write log to screen                                                     -->
		<log id="screenLog" logEvery="1000000">
			<column label="Joint" dp="4" width="12">
				<joint idref="joint"/>
			</column>
			<column label="Prior" dp="4" width="12">
				<prior idref="prior"/>
			</column>
			<column label="Likelihood" dp="4" width="12">
				<likelihood idref="likelihood"/>
			</column>
			
						<column label="clock.rate" sf="6" width="12">
				<parameter idref="clock.rate"/>
			</column>
		</log>
		
		<!-- write log to file                                                       -->
		<log id="fileLog" logEvery="500000" fileName="',season,"_",subtype,'.log" overwrite="false">
			<joint idref="joint"/>
			<prior idref="prior"/>
			<likelihood idref="likelihood"/>
			<parameter idref="clock.rate"/>
			
               ')
  for (i in 1:length(names)){
    if (lengths[[i]]>=50){
      str = paste0(str,gsub("tree1",names[[i]],'
                          <parameter idref="tree1.treeModel.rootHeight"/>
			<tmrcaStatistic idref="tree1.age(root)"/>
			<treeLengthStatistic idref="tree1.treeLength"/>

			<parameter idref="skygrid.precision"/>
            <parameter idref="tree1.skygrid.logPopSize"/>
            <parameter idref="tree1.skygrid.cutOff"/>
            <gmrfSkyGridLikelihood idref="tree1.skygrid"/>
		
			<treeDataLikelihood idref="tree1.treeLikelihood"/>
                          
                          '))
    } else {
      str = paste0(str,gsub("tree1",names[[i]],'
                          <parameter idref="tree1.treeModel.rootHeight"/>
			<tmrcaStatistic idref="tree1.age(root)"/>
			<treeLengthStatistic idref="tree1.treeLength"/>

	  	<coalescentLikelihood idref="tree1.coalescent"/>
		
			<treeDataLikelihood idref="tree1.treeLikelihood"/>
                          
                          '))
      
    }
  }
  
  str = paste0(str,'	</log>')
  
  for (i in names){
    str = paste0(str,gsub("tree1",i,'
                          
                          <!-- write tree log to file                                                  -->
		<logTree id="tree1.treeFileLog" logEvery="5000000" nexusFormat="true" fileName="tree1.trees" sortTranslationTable="true">
			<treeModel idref="tree1.treeModel"/>

			<joint idref="joint"/>
		</logTree>
                          
                          '))
  }
  
  str = paste0(str,'
               </mcmc>
	
	<report>
		<property name="timer">
			<mcmc idref="mcmc"/>
		</property>
	</report>
	</beast>
               ')
  
  return(str)
  
}


generateHeader <- function(){
  return(paste0('<?xml version="1.0" standalone="yes"?>

<!-- Generated by BEAUTi v1.10.4                                             -->
<!--       by Alexei J. Drummond, Andrew Rambaut and Marc A. Suchard         -->
<!--       Department of Computer Science, University of Auckland and        -->
<!--       Institute of Evolutionary Biology, University of Edinburgh        -->
<!--       David Geffen School of Medicine, University of California, Los Angeles-->
<!--       http://beast.community/                                           -->



<beast version="1.10.4">
'))
}


generateFullXML_subsamp <- function(clusterdata,season,subtype){
  
  names = lapply(clusterdata,function(x)x$cluster)#[14:21]
  taxa = lapply(clusterdata,function(x)x$tips)#[14:21]
  dates = lapply(clusterdata,function(x)x$dates)#[14:21]
  states = lapply(clusterdata,function(x)x$states)#[14:21]
  trees = lapply(clusterdata,function(x)x$tree)#[14:21]
  
  lengths = lapply(taxa,length)
  
  bifurcs = names
  for (i in 1:length(trees)){
    Nnode = trees[[i]]$Nnode
    Ntip = length(trees[[i]]$tip.label)
    if (Ntip==Nnode+1 | Ntip==2){bifurcs[[i]]=T} else{bifurcs[[i]]=F}
    
  }
  
  maxdate = max(unlist(dates))
  
  nsite = c(13427,13432,14398,14402)[match(subtype,c("H3N2","H1N1","Yam","Vic"))]
  
  str = generateHeader()
  for (i in 1:length(names)){
    str = paste0(str,generateTaxa(names[[i]],taxa[[i]],dates[[i]],states[[i]]))
  }
  for (i in 1:length(names)){
    str = paste0(str,
                 generateDataTree(names[[i]],trees[[i]]),
                 generateStartingTree(names[[i]],trees[[i]]))
  }
  
  str = paste0(str,'	<parameter id="clock.rate" value="0.002357" lower="0.0"/>')
  
  firstlarge = which(unlist(lengths)>=50)[1]
  if (!(is.na(firstlarge))){
    str = paste0(str,generateModelBlockFirst(names[[firstlarge]],mtdt,maxdate,nsite,c("skygrid","coal")[(lengths[[firstlarge]]<50)+1]))
  }
  
  rest = (1:length(trees))
  if (!(is.na(firstlarge))){rest=rest[-firstlarge]}
  
  if (length(rest)>1){
    for (i in 1:length(rest)){
      str = paste0(str,generateModelBlockSecond(names[[rest[i]]],mtdt,maxdate,nsite,c("skygrid","coal")[(lengths[[rest[i]]]<50)+1]))
    }
  }
  
  if (length(rest)==1){
    str = paste0(str,generateModelBlockFirst(names[[rest[i]]],mtdt,maxdate,nsite,c("skygrid","coal")[(lengths[[rest[i]]]<50)+1]))
    
  }
  
  str = paste0(str,generateAllOperators(names,bifurcs,lengths))
  str = paste0(str,generateMCMC(names,season,subtype,lengths))
  
  if (!(dir.exists(paste0("beast_xmls_subsamp/",season,"_",subtype)))){
    dir.create(paste0("beast_xmls_subsamp/",season,"_",subtype))
  }
  
  cat(str,file=paste0("beast_xmls_subsamp/",season,"_",subtype,'/',season,"_",subtype,'.xml'))
  
}



getTreeData <- function(cluster,subtype=NA,season=NA){
  
  taxa_to_keep = c()
  mtdt_subsampled = read.csv("/Users/simondejong/US_phylo/subsampled_mtdt_hhs_phylo.tsv",h=T,sep=',')
  mtdt2 = mtdt[mtdt$Strain%in%mtdt_subsampled$Strain,]
  tree = mtdt[mtdt$Cluster==cluster,]$File[1]
  tree = regmatches(tree,gregexpr(paste0("(?<=)",subtype,".*"),tree,perl=TRUE))[[1]]
  tree = gsub(".txt",".treefile",tree)
  treepath = tree
  tree = read.tree(paste0("divergence_trees/",tree))
  tree = keep.tip(tree,tree$tip.label[which(tree$tip.label %in% mtdt2[mtdt2$Cluster==cluster & !(mtdt2$State %in% c("Alaska","Hawaii")),]$Strain)])
  if (is.null(tree)){return()}
  tree = di2multi(tree,tol=0.00005)
  tips = tree$tip.label
  dates = sapply(tree$tip.label,function(x)mtdt2[mtdt2$Strain==x,]$Date)
  states = sapply(tree$tip.label,function(x)mtdt2[mtdt2$Strain==x,]$State)
  
  return(list(cluster=paste0("cluster_",cluster),tips=tips,dates=dates,states=states,tree=tree))
}


source("generateXML_HHS_subsamp.R")
#source("generateXML_State.R")

seasons = c(2014:2019,2022)

for (i in 1:length(seasons)){
  season = seasons[i]
  for (subtype in c("H3N2","H1N1","Vic","Yam")){
    
    clusters = cluster_df_country[cluster_df_country$Season == season & 
                                    cluster_df_country$Subtype==subtype & 
                                    cluster_df_country$Percentage_Season>0.5,]$Cluster
    
    clusters = clusters[clusters!="UNCLUSTERED"]
    
    if (length(clusters) == 0){next}
    
    l = lapply(clusters,getTreeData,subtype=subtype,season=season)
    lengths = unlist(lapply(l,function(x)length(x$tips)))
    l = l[which(lengths>2)]
    
    generateFullXML_subsamp(l,season,subtype)
    generateBEASTGeoXML_HHS_subsamp (l,season,subtype)
   # generateBEASTGeoXML(l,season,subtype)
    
  }
}
