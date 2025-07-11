
CalNetAttrib <- function(LinkNodeLs){
  edge <- LinkNodeLs[[1]]
  node <- LinkNodeLs[[2]]
  if (any(duplicated(node$id))){
    dupid <- node[which(duplicated(node$id)),]$id
    dupid <- unique(dupid)
    nodeundup <- node[!node$id %in% dupid,]
    for (id in dupid) {
      dupnode <- node[node$id %in% id,]
      dupnode <- dupnode[which(dupnode$ReInt==max(dupnode$ReInt)),][1,]
      nodeundup <- rbind(nodeundup,dupnode)
    }
    node <- nodeundup
  }
  edge <- edge[!(is.na(edge$Formula1)|is.na(edge$Formula2)),]
  network <- graph_from_data_frame(edge, directed = F,vertices = node)
  #network <- graph_from_data_frame(edge, directed = F)
  V(network)$taxa <- V(network)$name
  data("dataset")
  t1 <- trans_network$new(dataset = dataset, cor_method = "spearman", filter_thres = 0.001)
  t1$cal_network(COR_p_thres = 0.01, COR_cut = 0.7)
  networknew <- simplify(network)#删除multiedges和loop
  t1$res_network<-networknew
  t1$cal_module(method = "cluster_fast_greedy")
  t1$cal_network_attr()
  t1$get_node_table(node_roles = TRUE)
  p <- t1$plot_taxa_roles(use_type = 1)
  #Path <- sprintf("./output/%snodeRoles.PDF",SampleName)
  #ggsave(Path,dpi = 300,width = 11.69, height = 8.27, units = "in")
  network_attrib <- t1$res_network_attr
  node_roles <- t1$res_node_table
  node_roles <-node_roles[,1:10]
  colnames(node_roles)[10]<-"Node_Role"
  NetAttrib <- list(network_attrib,node_roles)
  names(NetAttrib) <- c("network_attrib","node_roles")
  #write.csv(NetAttrib$node_roles,sprintf("./output/%snodeRoles.csv",SampleName),row.names = FALSE)
  Res<- list(NetAttrib,p)
  return(Res)
}
# Keystone molecule 分布
#NetAtrribList <- list(AADNetAttrib$node_roles)

#NodeList <- list(AADnode)
#names(NodeList)<- c("AADnode")
#AAD <- read_xlsx("D:/Network/kegg/DOM-Microb/BeijingDigester/AAD.xlsx",col_names = T)
#FtMsFile <- list(AAD)
#names(FtMsFile)<-c("AAD")
KeystoneMolecule <- function(NetAtrribList,NodeList,FtMsFile){
  Res_keynode <- tibble(name = NA,degree=NA,betweenness=NA,module=NA,z=NA,p=NA,taxa_roles=NA,Label=NA,Formula=NA,ReInt=NA,group=NA,`O/C`=NA,`H/C`=NA,`N/C`=NA,`S/C`=NA,NOSC=NA,`DBE/C`=NA,Group=NA,groupKVD=NA,NewChemgroup=NA,sample=NA)
  for (i in 1:length(NetAtrribList)) {
    keynode <- NetAtrribList[[i]][!NetAtrribList[[i]]$taxa_roles %in% 'Peripheral nodes' & !NetAtrribList[[i]]$taxa_roles %in% NA, ]
    keynode <- keynode[c('name','degree','betweenness','module','z','p','taxa_roles')]
    Node <- NodeList[[i]][c('id','Label','group')]
    if (is.character(keynode$name[1])){
      Node$name <- as.character(Node$id)
    }else{
      Node$name <- Node$id
    }
    keynode1 <- left_join(keynode,Node,by='name')
    molecule <- FtMsFile[[i]]
    if (is.character(keynode1$name[1])){
      molecule$name <- as.character(round(molecule$`Neutral m/z`,digits = 6))
    }else{
      molecule$name <- round(molecule$`Neutral m/z`,digits = 6)
    }
    keynode2 <- left_join(keynode1,molecule,by='name')
    keynode2$sample <-names(FtMsFile)[i]
    keynode2 <- keynode2[c('name','degree','betweenness','module','z','p','taxa_roles','Label','Formula','ReInt','group','O/C',"H/C","N/C","S/C","NOSC","DBE/C","Group","groupKVD","NewChemgroup","sample")]
    Res_keynode <- rbind(Res_keynode,keynode2)
  }
  Res_keynode <- Res_keynode[2:nrow(Res_keynode),]
  #p <- ggscatter(Res_keynode,'O/C','N/C',color = 'sample',shape = 'sample', label = 'Formula',repel = T)+color_palette(c("#672494", "#9F2786", "#B94D63", "#BC792B", "#AAA217", "#83C86D", "#5CE6BE", "#B0F4FA"))+ theme_prism(base_size = 13,border = T)+scale_shape_manual(values = c(0,1,2,3,4,5,6,20))
  #ggpar(p,legend = 'top')
  #ggsave("./Keystone.pdf",width=8.5,height=7,dpi=300)
  #write.csv(Res_keynode,"./output/KeyNode.csv",row.names = F)
  return(Res_keynode)
}

PlotNodeAttribPropor <- function(NetworkFiles){
  count_elementGroup <- NetworkFiles$node_attrib %>% count(ElementGroup)
  count_elementGroup$Relat <- round(count_elementGroup$n*100/sum(count_elementGroup$n),digits = 2)
  count_elementGroup <- count_elementGroup[order(count_elementGroup$Relat,decreasing = T),]
  p1 <- ggpar(ggbarplot(count_elementGroup,'ElementGroup','Relat', fill = 'ElementGroup',size =8/70,alpha=0.6, width = 0.4,label = TRUE,lab.size = 1,lab.pos = 'in')+
          fill_palette('jco')+xlab(NULL)+ylab("Proportion (%)")+
          theme(axis.text.x = element_text(vjust=1,size = 6))+
          theme_prism(base_size = 7,border = T,base_fontface = "plain",base_line_size = 8/50,base_rect_size = 8/70)+
          guides(fill=guide_legend(keywidth=0.4,keyheight=0.3,size=8/70)),legend = 'top',x.text.angle = 20)
  #ggsave(sprintf("./output/%sNetNodeElemRatio.pdf",SampleName),width=8,height=7,dpi=300,units = 'cm')
  
  count_KVDgroup <- NetworkFiles$node_attrib %>% count(groupKVD)
  count_KVDgroup$Relat <- round(count_KVDgroup$n*100/sum(count_KVDgroup$n),digits = 2)
  count_KVDgroup <- count_KVDgroup[order(count_KVDgroup$Relat,decreasing = T),]
  p2 <- ggpar(ggbarplot(count_KVDgroup,'groupKVD','Relat', fill = 'groupKVD',size =8/70,alpha=0.6, width = 0.4,label = TRUE,lab.size = 1,lab.pos = 'in')+
          fill_palette('jco')+xlab(NULL)+ylab("Proportion (%)")+
          theme(axis.text.x = element_text(vjust=1,size = 5.5))+
          theme_prism(base_size = 7,border = T,base_fontface = "plain",base_line_size = 8/50,base_rect_size = 8/70)+
          guides(fill=guide_legend(keywidth=0.4,keyheight=0.3,size=8/70)),legend = 'top',x.text.angle = 20)
  #ggsave(sprintf("./output/%sNetNodeKVDRatio.pdf",SampleName),width=8,height=8,dpi=300,units = 'cm')
  
  count_KVD2group <- NetworkFiles$node_attrib %>% count(NewChemgroup)
  count_KVD2group$Relat <- round(count_KVD2group$n*100/sum(count_KVD2group$n),digits = 2)
  count_KVD2group <- count_KVD2group[order(count_KVD2group$Relat,decreasing = T),]
  p3 <- ggpar(ggbarplot(count_KVD2group,'NewChemgroup','Relat', fill = 'NewChemgroup',size =8/70,alpha=0.6, width = 0.4,label = TRUE,lab.size = 1,lab.pos = 'in')+
          fill_palette('jco')+xlab(NULL)+ylab("Proportion (%)")+
          theme(axis.text.x = element_text(vjust=1,size = 5.5))+
          theme_prism(base_size = 7,border = T,base_fontface = "plain",base_line_size = 8/50,base_rect_size = 8/70)+
          guides(fill=guide_legend(keywidth=0.4,keyheight=0.3,size=8/70)),legend = 'top',x.text.angle = 20)
  #ggsave(sprintf("./output/%sNetNodeKVD2Ratio.pdf",SampleName),width=8,height=8,dpi=300,units = 'cm')
  
  
  count_elementReact <- NetworkFiles$link_mass %>% count(ElementReact)
  count_elementReact$Relat <- round(count_elementReact$n*100/sum(count_elementReact$n),digits = 2)
  count_elementReact <- count_elementReact[order(count_elementReact$Relat,decreasing = T),]
  count_elementReactNew <- tibble(ElementReact=NA,n=NA,Relat=NA)
  for (i in 1:nrow(count_elementReact)) {
    if (!any(count_elementReactNew$ElementReact %in% paste(str_split(count_elementReact$ElementReact[i],pattern = '-')[[1]][2],str_split(count_elementReact$ElementReact[i],pattern = '-')[[1]][1],sep = '-' ))){
      if (!str_split(count_elementReact$ElementReact[i],pattern = '-')[[1]][2]==str_split(count_elementReact$ElementReact[i],pattern = '-')[[1]][1]){
        tmp <- paste(str_split(count_elementReact$ElementReact[i],pattern = '-')[[1]][2],str_split(count_elementReact$ElementReact[i],pattern = '-')[[1]][1],sep = '-' )
        if (any(count_elementReact$ElementReact %in% tmp)){
          data_tmp <- count_elementReact[count_elementReact$ElementReact %in% tmp,]
          n_tmp <- sum(count_elementReact$n[i],data_tmp$n)
          Relat_tmp <- sum(count_elementReact$Relat[i],data_tmp$Relat)
          count_elementReactNewtmp <- tibble(ElementReact=count_elementReact$ElementReact[i],n=n_tmp,Relat=Relat_tmp)
        }
      }else{
        count_elementReactNewtmp <- tibble(ElementReact=count_elementReact$ElementReact[i],n=count_elementReact$n[i],Relat=count_elementReact$Relat[i])
      }
      count_elementReactNew <- rbind(count_elementReactNew,count_elementReactNewtmp)
    }
  }
  count_elementReactNew <- count_elementReactNew[2:nrow(count_elementReactNew),]
  count_elementReactplot <- unique(count_elementReactNew[1:10,])
  
  p4 <- ggpar(ggbarplot(count_elementReactplot,'ElementReact','Relat', fill = 'ElementReact',size =8/70,alpha=0.6, width = 0.4,label = TRUE,lab.size = 1,lab.pos = 'in')+
          fill_palette('jco')+xlab(NULL)+ylab("Proportion (%)")+
          theme(axis.text.x = element_text(vjust=1,size = 5.5))+
          theme_prism(base_size = 7,border = T,base_fontface = "plain",base_line_size = 8/50,base_rect_size = 8/70)+
          guides(fill=guide_legend(keywidth=0.4,keyheight=0.3,size=8/70)),legend = 'top',x.text.angle = 15)
  #ggsave(sprintf("./output/%sNetLinkElemRatio.pdf",SampleName),width=8,height=8,dpi=300,units = 'cm')
  
  count_KVDReact <- NetworkFiles$link_mass %>% count(KVDReact)
  count_KVDReact$Relat <- round(count_KVDReact$n*100/sum(count_KVDReact$n),digits = 2)
  count_KVDReact <- count_KVDReact[order(count_KVDReact$Relat,decreasing = T),]
  count_KVDReactNew <- tibble(KVDReact=NA,n=NA,Relat=NA)
  for (i in 1:nrow(count_KVDReact)) {
    if (!any(count_KVDReactNew$KVDReact %in% paste(str_split(count_KVDReact$KVDReact[i],pattern = '-')[[1]][2],str_split(count_KVDReact$KVDReact[i],pattern = '-')[[1]][1],sep = '-' ))){
      if (!str_split(count_KVDReact$KVDReact[i],pattern = '-')[[1]][2]==str_split(count_KVDReact$KVDReact[i],pattern = '-')[[1]][1]){
        tmp <- paste(str_split(count_KVDReact$KVDReact[i],pattern = '-')[[1]][2],str_split(count_KVDReact$KVDReact[i],pattern = '-')[[1]][1],sep = '-' )
        if (any(count_KVDReact$KVDReact %in% tmp)){
          data_tmp <- count_KVDReact[count_KVDReact$KVDReact %in% tmp,]
          n_tmp <- sum(count_KVDReact$n[i],data_tmp$n)
          Relat_tmp <- sum(count_KVDReact$Relat[i],data_tmp$Relat)
          count_KVDReactNewtmp <- tibble(KVDReact=count_KVDReact$KVDReact[i],n=n_tmp,Relat=Relat_tmp)
        }
      }else{
        count_KVDReactNewtmp <- tibble(KVDReact=count_KVDReact$KVDReact[i],n=count_KVDReact$n[i],Relat=count_KVDReact$Relat[i])
      }
      count_KVDReactNew <- rbind(count_KVDReactNew,count_KVDReactNewtmp)
    }
  }
  count_KVDReactNew <- count_KVDReactNew[2:nrow(count_KVDReactNew),]
  count_KVDReactplot <- unique(count_KVDReactNew[1:8,])
  count_KVDReactplot <- count_KVDReactplot[order(count_KVDReactplot$Relat,decreasing = T),]
  p5 <- ggpar(ggbarplot(count_KVDReactplot,'KVDReact','Relat', fill = 'KVDReact',size =8/70,alpha=0.6, width = 0.4,label = TRUE,lab.size = 1,lab.pos = 'in')+
          fill_palette('jco')+xlab(NULL)+ylab("Proportion (%)")+
          theme(axis.text.x = element_text(vjust=1,size = 5.5))+
          theme_prism(base_size = 7,border = T,base_fontface = "plain",base_line_size = 8/50,base_rect_size = 8/70)+
          guides(fill=guide_legend(keywidth=0.4,keyheight=0.3,size=8/70)),legend = 'top',x.text.angle = 15)
  #ggsave(sprintf("./output/%sNetLinkKVDRatio.pdf",SampleName),width=8,height=8,dpi=300,units = 'cm')
  
  count_KVD2React <- NetworkFiles$link_mass %>% count(KVD2React)
  count_KVD2React$Relat <- round(count_KVD2React$n*100/sum(count_KVD2React$n),digits = 2)
  count_KVD2React <- count_KVD2React[order(count_KVD2React$Relat,decreasing = T),]
  count_KVD2ReactNew <- tibble(KVD2React=NA,n=NA,Relat=NA)
  for (i in 1:nrow(count_KVD2React)) {
    if (!any(count_KVD2ReactNew$KVD2React %in% paste(str_split(count_KVD2React$KVD2React[i],pattern = '-')[[1]][2],str_split(count_KVD2React$KVD2React[i],pattern = '-')[[1]][1],sep = '-' ))){
      if (!str_split(count_KVD2React$KVD2React[i],pattern = '-')[[1]][2]==str_split(count_KVD2React$KVD2React[i],pattern = '-')[[1]][1]){
        tmp <- paste(str_split(count_KVD2React$KVD2React[i],pattern = '-')[[1]][2],str_split(count_KVD2React$KVD2React[i],pattern = '-')[[1]][1],sep = '-' )
        if (any(count_KVD2React$KVD2React %in% tmp)){
          data_tmp <- count_KVD2React[count_KVD2React$KVD2React %in% tmp,]
          n_tmp <- sum(count_KVD2React$n[i],data_tmp$n)
          Relat_tmp <- sum(count_KVD2React$Relat[i],data_tmp$Relat)
          count_KVD2ReactNewtmp <- tibble(KVD2React=count_KVD2React$KVD2React[i],n=n_tmp,Relat=Relat_tmp)
        }
      }else{
        count_KVD2ReactNewtmp <- tibble(KVD2React=count_KVD2React$KVD2React[i],n=count_KVD2React$n[i],Relat=count_KVD2React$Relat[i])
      }
      count_KVD2ReactNew <- rbind(count_KVD2ReactNew,count_KVD2ReactNewtmp)
    }
  }
  count_KVD2ReactNew <- count_KVD2ReactNew[2:nrow(count_KVD2ReactNew),]
  count_KVD2Reactplot <- unique(count_KVD2ReactNew[1:8,])
  count_KVD2Reactplot <- count_KVD2Reactplot[order(count_KVD2Reactplot$Relat,decreasing = T),]
  p6 <- ggpar(ggbarplot(count_KVD2Reactplot,'KVD2React','Relat', fill = 'KVD2React',size =8/70,alpha=0.6, width = 0.4,label = TRUE,lab.size = 1,lab.pos = 'in')+
          fill_palette('jco')+xlab(NULL)+ylab("Proportion (%)")+
          theme(axis.text.x = element_text(vjust=1,size = 5.5))+
          theme_prism(base_size = 7,border = T,base_fontface = "plain",base_line_size = 8/50,base_rect_size = 8/70)+
          guides(fill=guide_legend(keywidth=0.4,keyheight=0.3,size=8/70)),legend = 'top',x.text.angle = 15)
  # ggsave(sprintf("./output/%sNetLinkKVD2Ratio.pdf",SampleName),width=8,height=8.5,dpi=300,units = 'cm')
  Res <- list(p1,p2,p3,p4,p5,p6)
  names(Res) <- c("NetNodeElemRatio","NetNodeChemicClassRatio","NetNodeChemicClass2Ratio","NetLinkElemRatio","NetLinkChemicClassRatio","NetLinkChemicClass2Ratio")
  return(Res)
}

PMDHeatmap <- function(NetLinkList,PMD_info,cluster_col=TRUE,TopNum=60){
  #dir.create("./Transformation")
  #colnames(PMD_info) <- c('PMD','group')
  if(any(str_detect(colnames(PMD_info),'Group'))){
    colnames(PMD_info)[which(colnames(PMD_info)=="Group")]<-"group"
  }
  if(any(PMD_info$group %in% 'C1')){
    PMD_info[PMD_info$group %in% 'C1',]$PMD <- 12.0
  }
  if(any(PMD_info$group %in% 'C')){
    PMD_info[PMD_info$group %in% 'C',]$PMD <- 12.0
  }
  PMD_info$PMD1 <- trunc(PMD_info$PMD*10000)/10000
  #PMD_info[PMD_info$group %in% 'C1',]$PMD1 <- PMD_info[PMD_info$group %in% 'C1',]$PMD
  PMD_res <- tibble(group=NA,n=NA,Relat=NA,sample=NA,PMD=NA)
  for (i in 1:length(NetLinkList)) {
    NetLinkList[[i]]$PMD1 <- trunc(NetLinkList[[i]]$PMD*10000)/10000
    data_tmp <- left_join(NetLinkList[[i]],PMD_info,by='PMD1')
    count_pmd <- data_tmp%>%count(group)
    count_pmd$Relat <- count_pmd$n/sum(count_pmd$n)
    count_pmd$sample <- names(NetLinkList)[i]
    count_pmd <- left_join(count_pmd,PMD_info[,c('group','PMD')],by='group')
    PMD_res <- rbind(PMD_res,count_pmd)
  }
  PMD_res <- PMD_res[2:nrow(PMD_res),]
  #write.csv(PMD_res,"./Transformation/PMDcount.csv",row.names = FALSE)
  PMD_resspread <- spread(PMD_res[,c('PMD','sample','n')],sample,n,fill = 0)
  #write.csv(PMD_resspread,"./Transformation/PMDcountSpread.csv",row.names = FALSE)
  PMD_res1 <- PMD_res[c('group','Relat','sample')]
  PMD_res2 <- spread(PMD_res1,sample,Relat)
  PMD_res2[is.na(PMD_res2)]<-0.0
  PMD_res3 <-as.data.frame(PMD_res2)
  rownames(PMD_res3)<-as.character(PMD_res3$group)
  PMD_res3 <- PMD_res3[,2:ncol(PMD_res3)]
  PMD_res3$sum <- rowSums(PMD_res3)
  PMD_res3 <- PMD_res3[order(PMD_res3$sum,decreasing = T),]
  if(nrow(PMD_res3)>TopNum){
    TopPMD_res <- PMD_res3[1:TopNum,1:(ncol(PMD_res3)-1)]*100
  }else{
    TopPMD_res <- PMD_res3[,1:(ncol(PMD_res3)-1)]*100
  }
  
  #p_tol <- pheatmap(TopPMD_res, cluster_cols = cluster_col,treeheight_col = 4,treeheight_row=15 ,fontsize=6, lwd=8/60, cellheight = 4.6, cellwidth = 4.6, fontsize_row = 5, border_color = "gray")
  p_tol <- pheatmap(TopPMD_res, cluster_cols = cluster_col, cellheight = 10, border_color = "gray",silent = TRUE,treeheight_col=30,treeheight_row=25)
  #ggsave("./Transformation/transformationHeatmap.pdf",plot = p_tol,width=7,height=10,dpi=300,units = 'cm')
  TopPMD_resBar <- round(PMD_res3[1:9,1:(ncol(PMD_res3)-1)]*100,digits = 2)
  TopPMD_resBar <- rownames_to_column(TopPMD_resBar,"Group")
  TopPMD_resBar <- gather(TopPMD_resBar,"Sample","RelativeAbundance",-1)
  p_tol_bar <- ggpar(ggbarplot(TopPMD_resBar,"Sample","RelativeAbundance",width = 0.4,alpha=0.6,size = 8/80,fill = "Group",order = names(NetLinkList),facet.by = "Group")+
          fill_palette(palette = "jco")+ylab("Proportion of frequency (%)")+xlab(NULL)+
          theme_prism(base_size = 8,base_fontface = 'plain',base_rect_size = 8/60,base_line_size = 8/60,border = TRUE)+
          theme(legend.position = 'top',legend.box.spacing = unit(0.01,'mm'),legend.key.height = unit(0.3,'cm'),strip.background = element_rect(fill = "grey"),strip.placement = "outside"),yticks.by = 0.5,xtickslab.rt = 20)
  #ggsave("./Transformation/TOPPMDbarplot.pdf",width=15,height=12,dpi=300,units = "cm")
  
  PMD_S <- PMD_res3[str_detect(rownames(PMD_res3),'S'),1:(ncol(PMD_res3)-1)]*100
  if(nrow(PMD_S)>=2){
    p_S <- pheatmap(PMD_S, cluster_cols = cluster_col, cellheight = 10, border_color = "gray",silent = TRUE,treeheight_col=30,treeheight_row=25)
    #ggsave("./Transformation/S_transformationHeatmap.pdf",plot = p_S, width=7,height=10,dpi=300,units = 'cm')
    PMD_SBar <-PMD_res3[str_detect(rownames(PMD_res3),'S'),1:(ncol(PMD_res3)-1)]*100
    PMD_SBar <- round(PMD_SBar[1:9,],2)
    PMD_SBar <- rownames_to_column(PMD_SBar,"Group")
    PMD_SBar <- gather(PMD_SBar,"Sample","RelativeAbundance",-1)
    p_S_bar <-ggpar(ggbarplot(PMD_SBar,"Sample","RelativeAbundance",width = 0.4,alpha=0.6,size = 8/80,fill = "Group",order = names(NetLinkList),facet.by = "Group")+
                      fill_palette(palette = "jco")+ylab("Proportion of frequency (%)")+xlab(NULL)+
                      theme_prism(base_size = 8,base_fontface = 'plain',base_rect_size = 8/60,base_line_size = 8/60,border = TRUE)+
                      theme(legend.position = 'top',legend.box.spacing = unit(0.01,'mm'),legend.key.height = unit(0.3,'cm'),strip.background = element_rect(fill = "grey"),strip.placement = "outside"),yticks.by = 0.2,xtickslab.rt = 20)
    #ggsave("./Transformation/TOPPMDSbarplot.pdf",width=15,height=12,dpi=300,units = "cm")
  }else{
    p_S <- NULL
    p_S_bar <- NULL
  }
  PMD_N <- PMD_res3[str_detect(rownames(PMD_res3),'N'),1:(ncol(PMD_res3)-1)]*100
  if(nrow(PMD_N)>=2){
    p_N <- pheatmap(PMD_N, cluster_cols = cluster_col, cellheight = 10, border_color = "gray",silent = TRUE,treeheight_col=30,treeheight_row=25)
    #ggsave("./Transformation/N_transformationHeatmap.pdf",plot = p_N,width=7,height=10,dpi=300,units = 'cm')
    
    PMD_NBar <-PMD_res3[str_detect(rownames(PMD_res3),'N'),1:(ncol(PMD_res3)-1)]*100
    PMD_NBar <- round(PMD_NBar[1:9,],2)
    PMD_NBar <- rownames_to_column(PMD_NBar,"Group")
    PMD_NBar <- gather(PMD_NBar,"Sample","RelativeAbundance",-1)
    p_N_bar <-ggpar(ggbarplot(PMD_NBar,"Sample","RelativeAbundance",width = 0.4,alpha=0.6,size = 8/80,fill = "Group",order = names(NetLinkList),facet.by = "Group")+
                      fill_palette(palette = "jco")+ylab("Proportion of frequency (%)")+xlab(NULL)+
                      theme_prism(base_size = 8,base_fontface = 'plain',base_rect_size = 8/60,base_line_size = 8/60,border = TRUE)+
                      theme(legend.position = 'top',legend.box.spacing = unit(0.01,'mm'),legend.key.height = unit(0.3,'cm'),strip.background = element_rect(fill = "grey"),strip.placement = "outside"),yticks.by = 0.2,xtickslab.rt = 20)
    #ggsave("./Transformation/TOPPMDNbarplot.pdf",width=15,height=12,dpi=300,units = "cm")
  }else{
    p_N <- NULL
    p_N_bar <- NULL
  }
  PMD_P <- log2(PMD_res3[str_detect(rownames(PMD_res3),'P'),1:(ncol(PMD_res3)-1)]+1)
  if(nrow(PMD_P)>=2){
    p_P <- pheatmap(PMD_P, cluster_cols = cluster_col, cellheight = 10, border_color = "gray",silent = TRUE,treeheight_col=30,treeheight_row=25)
    #ggsave("./P_transformationHeatmap.pdf",plot = p_P,width=8.27,height=11.69,dpi=300,units = 'in')
    
    PMD_PBar <-PMD_res3[str_detect(rownames(PMD_res3),'P'),1:(ncol(PMD_res3)-1)]*100
    PMD_PBar <- round(PMD_PBar[1:9,],2)
    PMD_PBar <- rownames_to_column(PMD_PBar,"Group")
    PMD_PBar <- gather(PMD_PBar,"Sample","RelativeAbundance",-1)
    p_P_bar <- ggpar(ggbarplot(PMD_PBar,"Sample","RelativeAbundance",width = 0.4,alpha=0.6,size = 8/80,fill = "Group",order = names(NetLinkList),facet.by = "Group")+
                       fill_palette(palette = "jco")+ylab("Proportion of frequency (%)")+xlab(NULL)+
                       theme_prism(base_size = 8,base_fontface = 'plain',base_rect_size = 8/60,base_line_size = 8/60,border = TRUE)+
                       theme(legend.position = 'top',legend.box.spacing = unit(0.01,'mm'),legend.key.height = unit(0.3,'cm'),strip.background = element_rect(fill = "grey"),strip.placement = "outside"),yticks.by = 0.2,xtickslab.rt = 20)
    #ggsave("./Transformation/TOPPMDPbarplot.pdf",width=14,height=11,dpi=300,units = "cm")
  }else{
    p_P <- NULL
    p_P_bar <- NULL
  }
  
  #write.csv(t(PMD_res3),"./transformation.csv",row.names = T)
  res <-as.data.frame(t(PMD_res3))
  res$sample <- rownames(res)
  res_p <- list(p_tol,p_N,p_S,p_P)
  names(res_p) <- c("p_tol","p_N","p_S","p_P")
  #write_xlsx(res,"./transformation.xlsx",col_names = T)
  return(res_p)
}

