###############################################
######MEN construction and the test of H1######
###############################################
#Ctrl+Shift+O:catalogue of this script
#R version 4.2.2

#Load required packages
library(dplyr)
library(stringr)
library(Matrix)
library(RMThreshold)
library(igraph)
library(reshape2)
library(tidyverse)
library(pulsar)
library(psych)
library(vegan)

#1.Generating the correlation matrix-----------------
#(1)Calculating Spearman correlations----------
#Using the LE of PG as examples
f_gin_le_corr <- corr.test(f_gin_le_0.01%>%t,method='spearman',ci=F,adjust='fdr')
#f_gin_le_0.01:ASV table of the LE of PG, only containing ASVs with relative abundance > 0.01%
#In this table,row represents ASVs,col represents samples

#(2)Filtering correlations based on FDR----------
f_gin_le_cor_p <- f_gin_le_corr$p
f_gin_le_cor_p[lower.tri(f_gin_le_cor_p)]=0#The upper tri of the p value matrix was adjusted
f_gin_le_cor_p <- f_gin_le_cor_p+t(f_gin_le_cor_p)

f_gin_le_matrix <- f_gin_le_corr$r
f_gin_le_matrix[f_gin_le_cor_p>0.05]=0

#(3)Selecting thresholds using RMT method----------
f_gin_le_rmt <- rm.get.threshold(f_gin_le_matrix,nr.thresholds = 100,
                               wait.seconds = 1,unfold.method = 'spline',
                               discard.zeros = T,interval = c(0.4,1.0))#0.6242

f_gin_le_matrix_filter <- rm.denoise.mat(f_gin_le_matrix,0.6242)%>%rm.discard.zeros()

#(4)Constructing MENs
f_gin_le_ig <- graph.adjacency(f_gin_le_matrix_filter,
                               mode = 'upper',
                               diag = F,
                               weighted = T)
#2.MEN analysis-----------------
#(1)Calculating topological properties for true MENs and randomly generated networks-----------
my_netproper_new <- function(ig_name){
  result <- data.frame(matrix(NA,nrow=length(ig_name),ncol = 22))
  rownames(result)=ig_name
  #calculating properties for true MENs
  for(i in 1:length(ig_name)){
    sample_ig <- get(ig_name[i])
    E(sample_ig)$weight <- NA
    result[i,'edge'] = length(E(sample_ig))#Number of edges
    result[i,'node'] = length(V(sample_ig))#Number of nodes
    result[i,'density'] = edge_density(sample_ig,loops = FALSE)#graph density
    result[i,'mean_degree'] <- mean(igraph::degree(sample_ig))#mean degree
    result[i,'path_length'] <- average.path.length(sample_ig)#avgL
    result[i,'transitivity'] <- transitivity(sample_ig) #clustering coefficient
    result[i,'deg_centralization'] <- centralization.degree(sample_ig)$centralization #centralization.degree度中心性
    greedy <- cluster_fast_greedy(sample_ig,weights = NULL)#detecting modules
    result[i,'greedy_modularity'] <- modularity(greedy)#calculating modularity
    result[i,'global_efficiency'] <- global_efficiency(sample_ig,directed = F)#global efficiency
    result[i,'natural_connectivity'] <- natural.connectivity(as_adjacency_matrix(sample_ig))#natural connectivity
    
    #Calculating properties for random networks
    #avgL,clustering coefficient,deg_centralization,greedy_modularity,global_efficiency,natural_connectivity
    random_path_length=c()
    random_transitivity=c()
    random_deg_centralization=c()
    random_modularity=c()
    random_efficiency=c()
    random_natural = c()
    for(j in 1:100){
      random_ig <- sample_gnm(length(V(sample_ig)),length(E(sample_ig)))
      random_path_length = c(random_path_length,average.path.length(random_ig))
      random_transitivity=c(random_transitivity,transitivity(random_ig))
      random_deg_centralization=c(random_deg_centralization,centralization.degree(random_ig)$centralization)
      random_modularity = c(random_modularity,modularity(cluster_fast_greedy(random_ig,weights = NULL)))
      random_efficiency = c(random_efficiency,global_efficiency(random_ig,directed = F))
      random_natural = c(random_natural,natural.connectivity(as_adjacency_matrix(random_ig)))
    }
    result[i,'mean_random_path_length'] = mean(random_path_length)
    result[i,'sd_random_path_length'] = sd(random_path_length)
    result[i,'mean_random_transitivity'] = mean(random_transitivity)
    result[i,'sd_random_transitivity'] = sd(random_transitivity)
    result[i,'mean_random_deg_centralization'] = mean(random_deg_centralization)
    result[i,'sd_random_deg_centralization'] = sd(random_deg_centralization)
    result[i,'mean_random_modularity'] = mean(random_modularity)
    result[i,'sd_random_modularity'] = sd(random_modularity)
    result[i,'mean_random_efficiency'] = mean(random_efficiency)
    result[i,'sd_random_efficiency'] = sd(random_efficiency)
    result[i,'mean_random_natural'] = mean(random_natural)
    result[i,'sd_random_natural'] = sd(random_natural)
    
  }
  return(result)
}

f_gin_le_net_proper <- my_netproper_new("f_gin_le_ig")

#(2)Calculating stability for networks with 20% of nodes being removed----------
my_stability <- function(ig_name,thresh=0.20){
  result <- data.frame(matrix(NA,nrow=length(ig_name),ncol=4))
  colnames(result) = c('ori_rm_nc','ori_rm_nc_sd','ran_rm_nc','ran_rm_nc_sd')
  rownames(result) = ig_name
  for(i in 1:length(ig_name)){
    print(i)
    ig = get(ig_name[i])
    ig_matrix <- as_adjacency_matrix(ig)%>%as.matrix
    
    node_num = length(V(ig))
    ori_result <- c()
    for(j in 1:100){
      #Extracting the ids of nodes which will be removed
      remove_node = sample(1:node_num,node_num*thresh)
      remove_ig_matrix = ig_matrix[-remove_node,-remove_node]
      remove_ig_matrix = rm.discard.zeros(remove_ig_matrix)
      ori_result <- c(ori_result,natural.connectivity(remove_ig_matrix))
    }
    result[i,1] = mean(ori_result)
    result[i,2] = sd(ori_result)
    
    #Calculation for random graphs
    random_result <- c()
    for(t in 1:100){
      random_ig = sample_gnm(length(V(ig)),length(E(ig)))
      random_matrix = as_adjacency_matrix(random_ig)%>%as.matrix
      for(x in 1:100){
        random_remove = sample(1:node_num,node_num*thresh)
        random_remove_matrix = random_matrix[-random_remove,-random_remove]
        random_remove_matrix = rm.discard.zeros(random_remove_matrix)
        random_result <- c(random_result,natural.connectivity(random_remove_matrix))
      }
    }
    result[i,3] = mean(random_result)
    result[i,4] = sd(random_result)
  }
  return(result)
}

f_gin_le_stability <- my_stability("f_gin_le_ig")

#(3)Identifying the roles of nodes in MENs------------
f_gin_le_ig_edge <- data.frame(as_edgelist(f_gin_le_ig),E(f_gin_le_ig)$weight)

my_node_attribute<-function(data,edge_table){
  #data:ig object
  #edge_table:edge table
  E(data)$weight <- NA 
  #(a)topological properties of nodes
  bet<-as.matrix(betweenness(data))
  degree<-as.matrix(degree(data))
  close<-as.matrix(closeness(data))
  eigen <- eigen_centrality(data)$vector %>% as.matrix
  module.greedy<-cluster_fast_greedy(data,weights=NULL)
  modularity.greedy<-modularity(module.greedy)
  module.greedy.table<-as.matrix(membership(module.greedy))
  attri<-cbind(degree,bet,close,eigen,module.greedy.table)
  attri<-as.data.frame(attri)
  names(attri)<-c("degree","bet","close",'eigen',"greedy.module")

  #(b)Zi and Pi of nodes
  my.zi<-function(data1,data2){
    #data1:node table
    #data2:edge table
    test<-data1
    test.edge<-data2
    test$Source<-row.names(test)
    test.edge<-merge(test.edge,test[,c(5,6)],by="Source")
    test$Target<-row.names(test)
    test.edge<-merge(test.edge,test[,c(5,7)],by="Target")
    test<-test[test$degree!=0,]
    max<-max(test$greedy.module)
    test.edge<-test.edge[,c(1,2,5,6)]
    #Extracting node information from each module
    null<-list()
    for(x in 1:max){
      mo_x<-test.edge[test.edge[,3]==x & test.edge[,4]==x,]
      null<-c(null,list(mo_x))
    }
    
    haha<-list()
    haha1<-data.frame()
    for(x in 1:max){
      t1<-data.frame(null[[x]]$Target,1)
      names(t1)<-c("node","number")
      t2<-data.frame(null[[x]]$Source,1)
      names(t2)<-c("node","number")
      mo_x<-rbind(t1,t2)
      mo_x<-aggregate(mo_x[,2],list(factor(mo_x[,1])),sum)
      mo_x[,2]<-(mo_x[,2]-mean(mo_x[,2]))/sd(mo_x[,2])
      haha<-c(haha,list(mo_x))
      haha1<-rbind(haha1,haha[[x]])
    }
    names(haha1)<-c("OTU","zi")
    return(haha1)
  }
  my.pi<-function(data1,data2){
    require(tidyverse)
    require(reshape2)
    test.node<-data1
    test.edge<-data2[,1:2]
    test.node$Source<-row.names(test.node)
    test.edge<-merge(test.edge,test.node[,c(5,6)],by="Source")
    names(test.edge)<-c("Source","Target","Smodule")
    test.node$Target<-row.names(test.node)
    test.edge<-merge(test.edge,test.node[,c(5,7)],by="Target")
    names(test.edge)<-c("Target","Source","Smodule","Tmodule")
    t1<-test.edge[,c(1,3)]
    t2<-test.edge[,c(2,4)]
    names(t1)<-c("ASV","Tmodule")
    names(t2)<-c("ASV","Tmodule")
    t<-rbind(t1,t2)
    t$calculate<-paste(t[,1],t[,2],sep="#module")
    t.new<-data.frame(t[,3],1)
    t.new<-aggregate(t.new[,2],list(factor(t.new[,1])),sum)
    ttt<-separate(t.new,1,c("ASV","Tmodule"),"#")
    ttt<-cast(ttt,ASV~Tmodule)
    ttt[is.na(ttt)]<-0
    ncol<-ncol(ttt)
    nrow<-nrow(ttt)
    for(x in 1:nrow){
      sum<-sum(ttt[x,2:ncol])
      n<-0
      for(y in 2:ncol){
        n1<-(ttt[x,y]/sum)^2
        n<-n+n1
      }
      ttt[x,ncol+1]<-1-n
    }
    return(ttt)
  }
  #merging data
  zi = my.zi(attri,edge_table)
  pi = my.pi(attri,edge_table)%>% select(.,ASV,pi=starts_with('V'))
  
  zp <- merge(zi,pi,by='ASV')
  name <- c("Row.names","degree","bet","close","eigen","greedy.module","zi","pi")
  
  attri <- merge(attri,zp,by.x='row.names',by.y='ASV');attri=attri[name]
  
  return(attri)
}

f_gin_le_ig_node <- my_node_attribute(f_gin_le_ig,f_gin_le_ig_edge)

#(3)The role of different factors in shaping MENs------------
my_lted <- function(edge_table,asv_table,bs_factor,site,net_thresh){
  #edge_table:edge table
  #asv_table:relative abundance table
  #bs_factor:edaphic factors
  #site:longitude and latitude
  #net_thresh为构建网络所使用的阈值
  #edaphic factors::pH,OC,TN,AK,AP,NIN,AMN,Mg,Ca,MC,S
  space_dist <- site[,1:2] %>% distm(.,fun = distVincentyEllipsoid)%>%as.dist
  edge_table$pH_direction=NA;edge_table$pH_cov=NA
  edge_table$OC_direction=NA;edge_table$OC_cov=NA
  edge_table$TN_direction=NA;edge_table$TN_cov=NA
  edge_table$AK_direction=NA;edge_table$AK_cov=NA
  edge_table$AP_direction=NA;edge_table$AP_cov=NA
  edge_table$NIN_direction=NA;edge_table$NIN_cov=NA
  edge_table$AMN_direction=NA;edge_table$AMN_cov=NA
  edge_table$Mg_direction=NA;edge_table$Mg_cov=NA
  edge_table$Ca_direction=NA;edge_table$Ca_cov=NA
  edge_table$MC_direction=NA;edge_table$MC_cov=NA
  edge_table$S_direction=NA;edge_table$S_cov=NA
  
  edge_table$source_space=NA;edge_table$target_space=NA;edge_table$space_cov=NA
  edge_table$env_cov=NA
  
  for(i in 1:nrow(edge_table)){
    print(paste(i,'in',nrow(edge_table),sep=' '))
    asv_source <- asv_table[edge_table[i,'Source'],1:27]%>%as.numeric
    asv_target <- asv_table[edge_table[i,'Target'],1:27]%>%as.numeric
    #pH
    asv_source_pH <- cor.test(asv_source,bs_factor$pH,method='spearman')
    asv_target_pH <- cor.test(asv_target,bs_factor$pH,method='spearman')
    if(asv_source_pH$p.value>0.05){asv_source_pH_rho=0}
    if(asv_source_pH$p.value<=0.05){asv_source_pH_rho=asv_source_pH$estimate}
    if(asv_target_pH$p.value>0.05){asv_target_pH_rho=0}
    if(asv_target_pH$p.value<=0.05){asv_target_pH_rho=asv_target_pH$estimate}
    if(abs(asv_source_pH_rho + asv_target_pH_rho) < max(abs(asv_source_pH_rho),abs(asv_target_pH_rho))){edge_table[i,'pH_direction']='oppo'}
    if(abs(asv_source_pH_rho + asv_target_pH_rho) > max(abs(asv_source_pH_rho),abs(asv_target_pH_rho))){edge_table[i,'pH_direction']='same'}
    if(abs(asv_source_pH_rho + asv_target_pH_rho) == max(abs(asv_source_pH_rho),abs(asv_target_pH_rho))){edge_table[i,'pH_direction']='all_zero'}
    if(edge_table[i,'pn']=='positive'& edge_table[i,'pH_direction']=='same' & abs(asv_source_pH_rho * asv_target_pH_rho)>=net_thresh){edge_table[i,'pH_cov']='True'}
    if(edge_table[i,'pn']=='negative'& edge_table[i,'pH_direction']=='oppo' & abs(asv_source_pH_rho * asv_target_pH_rho)>=net_thresh){edge_table[i,'pH_cov']='True'}
    #OC
    asv_source_OC <- cor.test(asv_source,bs_factor$OC,method='spearman')
    asv_target_OC <- cor.test(asv_target,bs_factor$OC,method='spearman')
    if(asv_source_OC$p.value>0.05){asv_source_OC_rho=0}
    if(asv_source_OC$p.value<=0.05){asv_source_OC_rho=asv_source_OC$estimate}
    if(asv_target_OC$p.value>0.05){asv_target_OC_rho=0}
    if(asv_target_OC$p.value<=0.05){asv_target_OC_rho=asv_target_OC$estimate}
    if(abs(asv_source_OC_rho + asv_target_OC_rho) < max(abs(asv_source_OC_rho),abs(asv_target_OC_rho))){edge_table[i,'OC_direction']='oppo'}
    if(abs(asv_source_OC_rho + asv_target_OC_rho) > max(abs(asv_source_OC_rho),abs(asv_target_OC_rho))){edge_table[i,'OC_direction']='same'}
    if(abs(asv_source_OC_rho + asv_target_OC_rho) == max(abs(asv_source_OC_rho),abs(asv_target_OC_rho))){edge_table[i,'OC_direction']='all_zero'}
    if(edge_table[i,'pn']=='positive'& edge_table[i,'OC_direction']=='same' & abs(asv_source_OC_rho * asv_target_OC_rho)>=net_thresh){edge_table[i,'OC_cov']='True'}
    if(edge_table[i,'pn']=='negative'& edge_table[i,'OC_direction']=='oppo' & abs(asv_source_OC_rho * asv_target_OC_rho)>=net_thresh){edge_table[i,'OC_cov']='True'}
    #TN
    asv_source_TN <- cor.test(asv_source,bs_factor$TN,method='spearman')
    asv_target_TN <- cor.test(asv_target,bs_factor$TN,method='spearman')
    if(asv_source_TN$p.value>0.05){asv_source_TN_rho=0}
    if(asv_source_TN$p.value<=0.05){asv_source_TN_rho=asv_source_TN$estimate}
    if(asv_target_TN$p.value>0.05){asv_target_TN_rho=0}
    if(asv_target_TN$p.value<=0.05){asv_target_TN_rho=asv_target_TN$estimate}
    if(abs(asv_source_TN_rho + asv_target_TN_rho) < max(abs(asv_source_TN_rho),abs(asv_target_TN_rho))){edge_table[i,'TN_direction']='oppo'}
    if(abs(asv_source_TN_rho + asv_target_TN_rho) > max(abs(asv_source_TN_rho),abs(asv_target_TN_rho))){edge_table[i,'TN_direction']='same'}
    if(abs(asv_source_TN_rho + asv_target_TN_rho) == max(abs(asv_source_TN_rho),abs(asv_target_TN_rho))){edge_table[i,'TN_direction']='all_zero'}
    if(edge_table[i,'pn']=='positive'& edge_table[i,'TN_direction']=='same' & abs(asv_source_TN_rho * asv_target_TN_rho)>=net_thresh){edge_table[i,'TN_cov']='True'}
    if(edge_table[i,'pn']=='negative'& edge_table[i,'TN_direction']=='oppo' & abs(asv_source_TN_rho * asv_target_TN_rho)>=net_thresh){edge_table[i,'TN_cov']='True'}
    #AK
    asv_source_AK <- cor.test(asv_source,bs_factor$AK,method='spearman')
    asv_target_AK <- cor.test(asv_target,bs_factor$AK,method='spearman')
    if(asv_source_AK$p.value>0.05){asv_source_AK_rho=0}
    if(asv_source_AK$p.value<=0.05){asv_source_AK_rho=asv_source_AK$estimate}
    if(asv_target_AK$p.value>0.05){asv_target_AK_rho=0}
    if(asv_target_AK$p.value<=0.05){asv_target_AK_rho=asv_target_AK$estimate}
    if(abs(asv_source_AK_rho + asv_target_AK_rho) < max(abs(asv_source_AK_rho),abs(asv_target_AK_rho))){edge_table[i,'AK_direction']='oppo'}
    if(abs(asv_source_AK_rho + asv_target_AK_rho) > max(abs(asv_source_AK_rho),abs(asv_target_AK_rho))){edge_table[i,'AK_direction']='same'}
    if(abs(asv_source_AK_rho + asv_target_AK_rho) == max(abs(asv_source_AK_rho),abs(asv_target_AK_rho))){edge_table[i,'AK_direction']='all_zero'}
    if(edge_table[i,'pn']=='positive'& edge_table[i,'AK_direction']=='same' & abs(asv_source_AK_rho * asv_target_AK_rho)>=net_thresh){edge_table[i,'AK_cov']='True'}
    if(edge_table[i,'pn']=='negative'& edge_table[i,'AK_direction']=='oppo' & abs(asv_source_AK_rho * asv_target_AK_rho)>=net_thresh){edge_table[i,'AK_cov']='True'}
    #AP
    asv_source_AP <- cor.test(asv_source,bs_factor$AP,method='spearman')
    asv_target_AP <- cor.test(asv_target,bs_factor$AP,method='spearman')
    if(asv_source_AP$p.value>0.05){asv_source_AP_rho=0}
    if(asv_source_AP$p.value<=0.05){asv_source_AP_rho=asv_source_AP$estimate}
    if(asv_target_AP$p.value>0.05){asv_target_AP_rho=0}
    if(asv_target_AP$p.value<=0.05){asv_target_AP_rho=asv_target_AP$estimate}
    if(abs(asv_source_AP_rho + asv_target_AP_rho) < max(abs(asv_source_AP_rho),abs(asv_target_AP_rho))){edge_table[i,'AP_direction']='oppo'}
    if(abs(asv_source_AP_rho + asv_target_AP_rho) > max(abs(asv_source_AP_rho),abs(asv_target_AP_rho))){edge_table[i,'AP_direction']='same'}
    if(abs(asv_source_AP_rho + asv_target_AP_rho) == max(abs(asv_source_AP_rho),abs(asv_target_AP_rho))){edge_table[i,'AP_direction']='all_zero'}
    if(edge_table[i,'pn']=='positive'& edge_table[i,'AP_direction']=='same' & abs(asv_source_AP_rho * asv_target_AP_rho)>=net_thresh){edge_table[i,'AP_cov']='True'}
    if(edge_table[i,'pn']=='negative'& edge_table[i,'AP_direction']=='oppo' & abs(asv_source_AP_rho * asv_target_AP_rho)>=net_thresh){edge_table[i,'AP_cov']='True'}
    #NIN
    asv_source_NIN <- cor.test(asv_source,bs_factor$NIN,method='spearman')
    asv_target_NIN <- cor.test(asv_target,bs_factor$NIN,method='spearman')
    if(asv_source_NIN$p.value>0.05){asv_source_NIN_rho=0}
    if(asv_source_NIN$p.value<=0.05){asv_source_NIN_rho=asv_source_NIN$estimate}
    if(asv_target_NIN$p.value>0.05){asv_target_NIN_rho=0}
    if(asv_target_NIN$p.value<=0.05){asv_target_NIN_rho=asv_target_NIN$estimate}
    if(abs(asv_source_NIN_rho + asv_target_NIN_rho) < max(abs(asv_source_NIN_rho),abs(asv_target_NIN_rho))){edge_table[i,'NIN_direction']='oppo'}
    if(abs(asv_source_NIN_rho + asv_target_NIN_rho) > max(abs(asv_source_NIN_rho),abs(asv_target_NIN_rho))){edge_table[i,'NIN_direction']='same'}
    if(abs(asv_source_NIN_rho + asv_target_NIN_rho) == max(abs(asv_source_NIN_rho),abs(asv_target_NIN_rho))){edge_table[i,'NIN_direction']='all_zero'}
    if(edge_table[i,'pn']=='positive'& edge_table[i,'NIN_direction']=='same' & abs(asv_source_NIN_rho * asv_target_NIN_rho)>=net_thresh){edge_table[i,'NIN_cov']='True'}
    if(edge_table[i,'pn']=='negative'& edge_table[i,'NIN_direction']=='oppo' & abs(asv_source_NIN_rho * asv_target_NIN_rho)>=net_thresh){edge_table[i,'NIN_cov']='True'}
    #AMN
    asv_source_AMN <- cor.test(asv_source,bs_factor$AMN,method='spearman')
    asv_target_AMN <- cor.test(asv_target,bs_factor$AMN,method='spearman')
    if(asv_source_AMN$p.value>0.05){asv_source_AMN_rho=0}
    if(asv_source_AMN$p.value<=0.05){asv_source_AMN_rho=asv_source_AMN$estimate}
    if(asv_target_AMN$p.value>0.05){asv_target_AMN_rho=0}
    if(asv_target_AMN$p.value<=0.05){asv_target_AMN_rho=asv_target_AMN$estimate}
    if(abs(asv_source_AMN_rho + asv_target_AMN_rho) < max(abs(asv_source_AMN_rho),abs(asv_target_AMN_rho))){edge_table[i,'AMN_direction']='oppo'}
    if(abs(asv_source_AMN_rho + asv_target_AMN_rho) > max(abs(asv_source_AMN_rho),abs(asv_target_AMN_rho))){edge_table[i,'AMN_direction']='same'}
    if(abs(asv_source_AMN_rho + asv_target_AMN_rho) == max(abs(asv_source_AMN_rho),abs(asv_target_AMN_rho))){edge_table[i,'AMN_direction']='all_zero'}
    if(edge_table[i,'pn']=='positive'& edge_table[i,'AMN_direction']=='same' & abs(asv_source_AMN_rho * asv_target_AMN_rho)>=net_thresh){edge_table[i,'AMN_cov']='True'}
    if(edge_table[i,'pn']=='negative'& edge_table[i,'AMN_direction']=='oppo' & abs(asv_source_AMN_rho * asv_target_AMN_rho)>=net_thresh){edge_table[i,'AMN_cov']='True'}
    #Mg
    asv_source_Mg <- cor.test(asv_source,bs_factor$Mg,method='spearman')
    asv_target_Mg <- cor.test(asv_target,bs_factor$Mg,method='spearman')
    if(asv_source_Mg$p.value>0.05){asv_source_Mg_rho=0}
    if(asv_source_Mg$p.value<=0.05){asv_source_Mg_rho=asv_source_Mg$estimate}
    if(asv_target_Mg$p.value>0.05){asv_target_Mg_rho=0}
    if(asv_target_Mg$p.value<=0.05){asv_target_Mg_rho=asv_target_Mg$estimate}
    if(abs(asv_source_Mg_rho + asv_target_Mg_rho) < max(abs(asv_source_Mg_rho),abs(asv_target_Mg_rho))){edge_table[i,'Mg_direction']='oppo'}
    if(abs(asv_source_Mg_rho + asv_target_Mg_rho) > max(abs(asv_source_Mg_rho),abs(asv_target_Mg_rho))){edge_table[i,'Mg_direction']='same'}
    if(abs(asv_source_Mg_rho + asv_target_Mg_rho) == max(abs(asv_source_Mg_rho),abs(asv_target_Mg_rho))){edge_table[i,'Mg_direction']='all_zero'}
    if(edge_table[i,'pn']=='positive'& edge_table[i,'Mg_direction']=='same' & abs(asv_source_Mg_rho * asv_target_Mg_rho)>=net_thresh){edge_table[i,'Mg_cov']='True'}
    if(edge_table[i,'pn']=='negative'& edge_table[i,'Mg_direction']=='oppo' & abs(asv_source_Mg_rho * asv_target_Mg_rho)>=net_thresh){edge_table[i,'Mg_cov']='True'}
    #Ca
    asv_source_Ca <- cor.test(asv_source,bs_factor$Ca,method='spearman')
    asv_target_Ca <- cor.test(asv_target,bs_factor$Ca,method='spearman')
    if(asv_source_Ca$p.value>0.05){asv_source_Ca_rho=0}
    if(asv_source_Ca$p.value<=0.05){asv_source_Ca_rho=asv_source_Ca$estimate}
    if(asv_target_Ca$p.value>0.05){asv_target_Ca_rho=0}
    if(asv_target_Ca$p.value<=0.05){asv_target_Ca_rho=asv_target_Ca$estimate}
    if(abs(asv_source_Ca_rho + asv_target_Ca_rho) < max(abs(asv_source_Ca_rho),abs(asv_target_Ca_rho))){edge_table[i,'Ca_direction']='oppo'}
    if(abs(asv_source_Ca_rho + asv_target_Ca_rho) > max(abs(asv_source_Ca_rho),abs(asv_target_Ca_rho))){edge_table[i,'Ca_direction']='same'}
    if(abs(asv_source_Ca_rho + asv_target_Ca_rho) == max(abs(asv_source_Ca_rho),abs(asv_target_Ca_rho))){edge_table[i,'Ca_direction']='all_zero'}
    if(edge_table[i,'pn']=='positive'& edge_table[i,'Ca_direction']=='same' & abs(asv_source_Ca_rho * asv_target_Ca_rho)>=net_thresh){edge_table[i,'Ca_cov']='True'}
    if(edge_table[i,'pn']=='negative'& edge_table[i,'Ca_direction']=='oppo' & abs(asv_source_Ca_rho * asv_target_Ca_rho)>=net_thresh){edge_table[i,'Ca_cov']='True'}
    #MC
    asv_source_MC <- cor.test(asv_source,bs_factor$MC,method='spearman')
    asv_target_MC <- cor.test(asv_target,bs_factor$MC,method='spearman')
    if(asv_source_MC$p.value>0.05){asv_source_MC_rho=0}
    if(asv_source_MC$p.value<=0.05){asv_source_MC_rho=asv_source_MC$estimate}
    if(asv_target_MC$p.value>0.05){asv_target_MC_rho=0}
    if(asv_target_MC$p.value<=0.05){asv_target_MC_rho=asv_target_MC$estimate}
    if(abs(asv_source_MC_rho + asv_target_MC_rho) < max(abs(asv_source_MC_rho),abs(asv_target_MC_rho))){edge_table[i,'MC_direction']='oppo'}
    if(abs(asv_source_MC_rho + asv_target_MC_rho) > max(abs(asv_source_MC_rho),abs(asv_target_MC_rho))){edge_table[i,'MC_direction']='same'}
    if(abs(asv_source_MC_rho + asv_target_MC_rho) == max(abs(asv_source_MC_rho),abs(asv_target_MC_rho))){edge_table[i,'MC_direction']='all_zero'}
    if(edge_table[i,'pn']=='positive'& edge_table[i,'MC_direction']=='same' & abs(asv_source_MC_rho * asv_target_MC_rho)>=net_thresh){edge_table[i,'MC_cov']='True'}
    if(edge_table[i,'pn']=='negative'& edge_table[i,'MC_direction']=='oppo' & abs(asv_source_MC_rho * asv_target_MC_rho)>=net_thresh){edge_table[i,'MC_cov']='True'}
    #S
    asv_source_S <- cor.test(asv_source,bs_factor$S,method='spearman')
    asv_target_S <- cor.test(asv_target,bs_factor$S,method='spearman')
    if(asv_source_S$p.value>0.05){asv_source_S_rho=0}
    if(asv_source_S$p.value<=0.05){asv_source_S_rho=asv_source_S$estimate}
    if(asv_target_S$p.value>0.05){asv_target_S_rho=0}
    if(asv_target_S$p.value<=0.05){asv_target_S_rho=asv_target_S$estimate}
    if(abs(asv_source_S_rho + asv_target_S_rho) < max(abs(asv_source_S_rho),abs(asv_target_S_rho))){edge_table[i,'S_direction']='oppo'}
    if(abs(asv_source_S_rho + asv_target_S_rho) > max(abs(asv_source_S_rho),abs(asv_target_S_rho))){edge_table[i,'S_direction']='same'}
    if(abs(asv_source_S_rho + asv_target_S_rho) == max(abs(asv_source_S_rho),abs(asv_target_S_rho))){edge_table[i,'S_direction']='all_zero'}
    if(edge_table[i,'pn']=='positive'& edge_table[i,'S_direction']=='same' & abs(asv_source_S_rho * asv_target_S_rho)>=net_thresh){edge_table[i,'S_cov']='True'}
    if(edge_table[i,'pn']=='negative'& edge_table[i,'S_direction']=='oppo' & abs(asv_source_S_rho * asv_target_S_rho)>=net_thresh){edge_table[i,'S_cov']='True'}
    
    #Spatial distance
    asv_source_dist = data.frame(asv=asv_source,other=0)%>%vegdist()
    asv_source_dist[is.na(asv_source_dist)]=0
    asv_target_dist = data.frame(asv=asv_target,other=0)%>%vegdist()
    asv_target_dist[is.na(asv_target_dist)]=0
    source_space = mantel(asv_source_dist,space_dist,method = 'spearman')
    target_space = mantel(asv_target_dist,space_dist,method = 'spearman')
    if(source_space$signif>0.05){edge_table[i,'source_space']=0}
    if(source_space$signif<=0.05){edge_table[i,'source_space']=source_space$statistic}
    if(target_space$signif>0.05){edge_table[i,'target_space']=0}
    if(target_space$signif<=0.05){edge_table[i,'target_space']=target_space$statistic}
  }
  edge_table[is.na(edge_table)]='FALSE'
  for(i in 1:nrow(edge_table)){
    edge_table[i,'env_cov']=str_count(edge_table[i,],'True')%>%sum
  }
  for(i in 1:nrow(edge_table)){
    if(edge_table[i,'source_space']>=0.6 & edge_table[i,'target_space']>=0.6){edge_table[i,'space_cov']='Yes'}
  }
  edge_table$space_cov[is.na(edge_table$space_cov)]='No'
  return(edge_table)
}

f_gin_le_ig_edge_lted <- my_lted(f_gin_le_ig_edge,f_gin_le_0.01,gin_bs_factor,gin_site,0.6242)