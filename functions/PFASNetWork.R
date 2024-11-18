
PFASDataProcess <- function(data,checkdup=FALSE){
  if(str_detect(str_flatten(colnames(data)),"Measured Mass")){
    colnames(data)[which(colnames(data)=="Measured Mass")]<-'Measured m/z'
  }else if(str_detect(str_flatten(colnames(data)),"mz")){
    colnames(data)[which(colnames(data)=="mz")]<-'Measured m/z'
  }else if(str_detect(str_flatten(colnames(data)),"m/z")){
    colnames(data)[which(colnames(data)=="m/z")]<-'Measured m/z'
  }
  data$`Measured m/z` <- trunc(data$`Measured m/z`*100000)/100000
  
  if(str_detect(str_flatten(colnames(data)),"mf")){
    colnames(data)[which(colnames(data)=="mf")]<-'Formula'
  }
  
  if(str_detect(str_flatten(colnames(data)),"S/N")){
    data <- filter(data,`S/N`>= 4)
  }
  if(any(colnames(data) %in% "I")&!any(colnames(data) %in% "Intensity")){
    colnames(data)[which(colnames(data)=="I")]<-'Intensity'
  }
  if(str_detect(str_flatten(colnames(data)),"Intensity")&(!any((colnames(data) %in% "RA")))){
    data$RA <- data$Intensity*100/(max(data$Intensity))
  }
  # 1.0072765为H-e的质量
  if(!any(colnames(data) %in% "KMD_CF2")){
    data$`KMD_CF2` <- (data$`Measured m/z`+1.0072765)*(50/49.9968)-trunc((data$`Measured m/z`+1.0072765)*(50/49.9968))
    #data$`Ave_KMD_CF2` <- trunc(data$`KMD_CF2`*1000)/1000
    data$`Kendrick_Mass_CF2` <- (data$`Measured m/z`+1.0072765)*(50/49.9968)
  }
  
  if(str_detect(str_flatten(colnames(data)),"Formula")){
    data <- data[c("Measured m/z","S/N","Intensity","RA","Formula","KMD_CF2","Kendrick_Mass_CF2")]
  }else{
    data <- data[c("Measured m/z","S/N","Intensity","RA","KMD_CF2","Kendrick_Mass_CF2")]
  }
  
  if(checkdup==TRUE){
    if(any(duplicated(data$`Measured m/z`))){
      dupid <- which(duplicated(data$`Measured m/z`))
      undupdata <- data[!data$`Measured m/z` %in% data$`Measured m/z`[dupid],]
      for (id in dupid) {
        dupdata <- data[data$`Measured m/z` %in% data$`Measured m/z`[id],]
        mergedata <- dupdata[1,]
        undupdata <- rbind(undupdata,mergedata)
      }
      data <-unique(undupdata)
      
    }
  }
  return(data)
}

f_calc_masses_dif <- function(peaklist){
  diff_matrix <- outer(peaklist$`Measured m/z`,peaklist$`Measured m/z`,'-')
  diff_matrix <- abs(diff_matrix)
  diff_matrix[!upper.tri(diff_matrix,diag = TRUE)] <- 0
  diff_DF <-  as.data.frame(diff_matrix)
  colnames(diff_DF)<- peaklist$`Measured m/z`
  diff_DF$Mass1 <- peaklist$`Measured m/z`
  pmd_result <- gather(diff_DF,Mass2,PMD,-ncol(diff_DF))
  pmd_result <- setDT(pmd_result)
  pmd_result <- pmd_result[PMD>0]
  pmd_result$Mass2 <- as.numeric(pmd_result$Mass2)
  return(pmd_result)
}

FindPMDNew2 <- function(peaklist,PMDDF,masserror=1,ndistance="norestrict"){
  pmdlist <- f_calc_masses_dif(peaklist)
  # sort pmdlist by mass
  setkey(pmdlist,PMD)
  # assign pmdlist masses
  s <- pmdlist$PMD
  n <- length(s)
  
  # assign pmd masses
  PMDDF <- setDT(PMDDF)
  setkey(PMDDF)
  d <- PMDDF$PMD
 
  pmd_find <- tibble(Mass1=NA,Mass2=NA,PMD_cal=NA,PMD=NA)
  for (i in 1:length(d)) {
    if(ndistance=="norestrict"){
      tmp <- pmdlist%>% mutate(DiffPMD=abs(1-PMD/(d[i]*round(PMD/d[i])))*1000000)%>% mutate(DiffPMD1=abs(PMD-d[i]))
    }else{
      pmdlist <- filter(pmdlist,round(PMD/d[i])<= as.numeric(ndistance))
      tmp <- pmdlist%>% mutate(DiffPMD=abs(1-PMD/(d[i]*round(PMD/d[i])))*1000000)%>% mutate(DiffPMD1=abs(PMD-d[i]))
    }
    
    if(masserror<0.005){
      tmp <- filter(tmp,DiffPMD1<= masserror)
    }else{
      tmp <- filter(tmp,DiffPMD<= masserror)
    }
    if(!nrow(tmp)==0){
      PMD_findtmp <- tmp[,c("Mass1","Mass2","PMD")]
      PMD_findtmp$PMDRef <- d[i]
      colnames(PMD_findtmp) <- colnames(pmd_find)
      pmd_find <- rbind(pmd_find,PMD_findtmp)
    }
  }
  pmd_find <- unique(pmd_find[2:nrow(pmd_find),])
  return(pmd_find)
}  

FindpfasPMD_betw_samp <- function(DF1,DF2,PMDDF,masserror=1,cal_nontrans='yes'){
  pmdlist <- f_calc_masses_dif_betw_sam(DF1,DF2)
  # sort pmdlist by mass
  setkey(pmdlist,PMD)
  # assign pmdlist masses
  s <- pmdlist$PMD
  n <- length(s)
  
  # assign pmd masses
  PMDDF <- setDT(PMDDF)
  setkey(PMDDF)
  d <- PMDDF$PMD
  # ranges
  r.d <- range(d)
  r.s <- range(s)
  #if (r.s[1] < r.d[1]){
    #pmdlist <- pmdlist[PMD>=(r.d[1]-0.01),]
    # assign peaklist masses
    #s <- pmdlist$PMD
    #n <- length(s)
    # ranges
    #r.s <- range(s)
    
  #}
  
  if (r.s[2] > r.d[2]){
    pmdlist <- pmdlist[PMD<=(r.d[2]+0.01),]
    # assign peaklist masses
    s <- pmdlist$PMD
    n <- length(s)
    # ranges
    r.s <- range(s)
    
  }
  
  pmd_find <- tibble(Mass1=NA,Mass2=NA,PMD_cal=NA,PMD=NA)
  for (i in 1:length(d)) {
    tmp <- pmdlist%>% mutate(DiffPMD=abs((d[i]-PMD)/d[i])*1000000)
    tmp <- filter(tmp,DiffPMD<= masserror)
    if(!nrow(tmp)==0){
      PMD_findtmp <- tmp[,c("Mass1","Mass2","PMD")]
      PMD_findtmp$PMDRef <- d[i]
      colnames(PMD_findtmp) <- colnames(pmd_find)
      pmd_find <- rbind(pmd_find,PMD_findtmp)
    }
  }
  pmd_find <- unique(pmd_find[2:nrow(pmd_find),])
  cal_nontrans=cal_nontrans
  if(cal_nontrans=='yes'){
    tmp2 <- filter(pmdlist,PMD <= 0.003)
    print(nrow(tmp2))
    if(!nrow(tmp2)==0){
      PMD_findtmp2 <- tmp2[,c("Mass1","Mass2","PMD")]
      PMD_findtmp2$PMDRef <- '<0.003'
      colnames(PMD_findtmp2) <- colnames(pmd_find)
      pmd_find <- rbind(pmd_find,PMD_findtmp2)
    }
  }
  
  return(pmd_find)
}



PFASNetwork <- function(FtMsDf,PMD,type="Res",mserror=0.5,plotgraph=FALSE,plotFormula=FALSE,ndistance="norestrict"){
  
  Links <- FindPMDNew2(FtMsDf[c('Measured m/z')],PMD=PMD,masserror=mserror,ndistance=ndistance)
  Nodes <- unique(as.data.frame(c(Links$Mass1,Links$Mass2)))
  colnames(Nodes)<-"Measured m/z"
  NodesNew <- unique(left_join(Nodes,FtMsDf,by='Measured m/z'))
  #NodesNew <- filter(NodesNew,!is.na(NodesNew$`S/N`))
  if(any(duplicated(NodesNew$`Measured m/z`))){
    dupid <- which(duplicated(NodesNew$`Measured m/z`))
    undupdata <- NodesNew[!NodesNew$`Measured m/z` %in% NodesNew$`Measured m/z`[dupid],]
    for (id in dupid) {
      dupdata <- NodesNew[NodesNew$`Measured m/z` %in% NodesNew$`Measured m/z`[id],]
      mergedata <- dupdata[1,]
      mergedata$Formula <- str_flatten(dupdata$Formula,collapse =';')
      undupdata <- rbind(undupdata,mergedata)
    }
    NodesNew <-unique(undupdata)
  }
  colnames(Links)[1:2]<-c('source','target')
  colnames(NodesNew)[1]<-c('id')
  
  if(plotFormula==TRUE){
    Network <- graph_from_data_frame(Links,vertices = NodesNew,directed=F)
  }else{
    Network <- graph_from_data_frame(Links,directed=F,vertices = NodesNew)
  }
  
  if(plotgraph==TRUE){
    if(plotFormula==TRUE){
      pdf(sprintf("./Netoutput/%s_%sFormula.pdf",type,names(FtMsDf)[i]))
      plot(Network, layout=layout_with_fr, vertex.size=4,vertex.label=vertex_attr(Network)$Formula,
           vertex.label.dist=0.5, vertex.label.color='blue',vertex.color="grey",vertex.label.cex=0.4,
           edge.label=edge_attr(Network)$PMD_cal,edge.label.cex=0.4,edge.label.color="red" )
      dev.off()
    }else{
      pdf(sprintf("./Netoutput/%s_%s.pdf",type,names(FtMsDf)[i]))
      plot(Network, layout=layout_with_fr, vertex.size=4,
           vertex.label.dist=0.5,vertex.label.color='blue', vertex.color="grey",vertex.label.cex=0.6,
           edge.label=edge_attr(Network)$PMD_cal,edge.label.cex=0.4,edge.label.color="red" )
      dev.off()
    }
  }
  Network<-simplify(Network)
  clust <- cluster_fast_greedy(Network)
  module <- membership(clust)
  Dfmodule <- as.data.frame(names(module))
  colnames(Dfmodule)[1]<-'id'
  Dfmodule$module <- module
  node_deg <- degree(Network)
  node_deg <- as.data.frame(node_deg)
  colnames(node_deg)[1] <- "degree"
  Dfmodule$degree <- node_deg$degree
  com_size <- sizes(clust)
  Dfmodule$module_size <- as.data.frame(com_size[Dfmodule$module])$Freq
  NodesNew$id <- as.character(NodesNew$id)
  NodesNew <- left_join(NodesNew,Dfmodule,by='id')
  NodesNew <- setDT(NodesNew)
  setkey(NodesNew,degree)
  NodesNew <- setDF(NodesNew)
  NodesNew <- setorder(NodesNew,-degree)
  NodesNew$id <- as.double(NodesNew$id)
  New_network <- graph_from_data_frame(Links,vertices = NodesNew)
  NodesNew <- as.data.frame(NodesNew)
  Netlisttmp <- list(Links,NodesNew,New_network)
  names(Netlisttmp) <- c('Links','Nodes','Network')
  return(Netlisttmp)
}

PFASNetworkCrossSamp <- function(DF1,DF2,PMD,type="Res",mserror=0.5,cal_nontrans='yes',plotgraph=FALSE,plotFormula=FALSE){
  
  Links <- FindpfasPMD_betw_samp(DF1,DF2,PMD=PMD,masserror=mserror,cal_nontrans=cal_nontrans)
  Nodes <- unique(as.data.frame(c(Links$Mass1,Links$Mass2)))
  colnames(Nodes)<-"Measured m/z"
  if(any(colnames(DF1)%in%'id')){
    DF1$`Measured m/z` <- DF1$id
    DF2$`Measured m/z` <- DF2$id
  }
  DFCombine <- rbind(DF1,DF2)
  NodesNew <- unique(left_join(Nodes,DFCombine,by='Measured m/z'))
  if(any(duplicated(NodesNew$`Measured m/z`))){
    dupid <- which(duplicated(NodesNew$`Measured m/z`))
    undupdata <- NodesNew[!NodesNew$`Measured m/z` %in% NodesNew$`Measured m/z`[dupid],]
    for (id in dupid) {
      dupdata <- NodesNew[NodesNew$`Measured m/z` %in% NodesNew$`Measured m/z`[id],]
      mergedata <- dupdata[1,]
      mergedata$Formula <- str_flatten(dupdata$Formula,collapse =';')
      undupdata <- rbind(undupdata,mergedata)
    }
    NodesNew <-unique(undupdata)
  }
  colnames(Links)[1:2]<-c('source','target')
  colnames(NodesNew)[1]<-c('id')
  
  if(plotFormula==TRUE){
    Network <- graph_from_data_frame(Links,vertices = NodesNew,directed=F)
  }else{
    Network <- graph_from_data_frame(Links,directed=F)
  }
  
  if(plotgraph==TRUE){
    if(plotFormula==TRUE){
      pdf(sprintf("./Netoutput/%s_%sFormula.pdf",type,names(FtMsDf)[i]))
      plot(Network, layout=layout_with_fr, vertex.size=4,vertex.label=vertex_attr(Network)$Formula,
           vertex.label.dist=0.5, vertex.label.color='blue',vertex.color="grey",vertex.label.cex=0.4,
           edge.label=edge_attr(Network)$PMD_cal,edge.label.cex=0.4,edge.label.color="red" )
      dev.off()
    }else{
      pdf(sprintf("./Netoutput/%s_%s.pdf",type,names(FtMsDf)[i]))
      plot(Network, layout=layout_with_fr, vertex.size=4,
           vertex.label.dist=0.5,vertex.label.color='blue', vertex.color="grey",vertex.label.cex=0.6,
           edge.label=edge_attr(Network)$PMD_cal,edge.label.cex=0.4,edge.label.color="red" )
      dev.off()
    }
  }
  Network<-simplify(Network)
  clust <- cluster_fast_greedy(Network)
  module <- membership(clust)
  Dfmodule <- as.data.frame(names(module))
  colnames(Dfmodule)[1]<-'id'
  Dfmodule$module <- module
  node_deg <- degree(Network)
  node_deg <- as.data.frame(node_deg)
  colnames(node_deg)[1] <- "degree"
  Dfmodule$degree <- node_deg$degree
  com_size <- sizes(clust)
  Dfmodule$module_size <- as.data.frame(com_size[Dfmodule$module])$Freq
  NodesNew$id <- as.character(NodesNew$id)
  NodesNew <- left_join(NodesNew,Dfmodule,by='id')
  NodesNew <- setDT(NodesNew)
  setkey(NodesNew,degree)
  NodesNew <- setDF(NodesNew)
  NodesNew <- setorder(NodesNew,-degree)
  NodesNew$id <- as.numeric(NodesNew$id)
  New_network <- graph_from_data_frame(Links,vertices = NodesNew)
  NodesNew <- as.data.frame(NodesNew)
  Netlisttmp <- list(Links,NodesNew,New_network)
  names(Netlisttmp) <- c('Links','Nodes','Network')
  return(Netlisttmp)
}

#在筛查基础上的样品间分子差网络关系
PFASNetworkCrossSamp2 <- function(DF1,DF2,PMD,type="Res",mserror=0.5,cal_nontrans='yes'){
  
  Links <- FindpfasPMD_betw_samp(DF1,DF2,PMD=PMD,masserror=mserror,cal_nontrans=cal_nontrans)
  Nodes <- unique(as.data.frame(c(Links$Mass1,Links$Mass2)))
  colnames(Nodes)<-"Measured m/z"
  if(any(colnames(DF1)%in%'id')){
    DF1$`Measured m/z` <- DF1$id
    DF2$`Measured m/z` <- DF2$id
  }
  DFCombine <- rbind(DF1,DF2)
  NodesNew <- unique(left_join(Nodes,DFCombine,by='Measured m/z'))
  if(any(duplicated(NodesNew$`Measured m/z`))){
    dupid <- which(duplicated(NodesNew$`Measured m/z`))
    undupdata <- NodesNew[!NodesNew$`Measured m/z` %in% NodesNew$`Measured m/z`[dupid],]
    for (id in dupid) {
      dupdata <- NodesNew[NodesNew$`Measured m/z` %in% NodesNew$`Measured m/z`[id],]
      mergedata <- dupdata[1,]
      mergedata$Formula <- str_flatten(unique(dupdata$Formula),collapse =';')
      undupdata <- rbind(undupdata,mergedata)
    }
    NodesNew <-unique(undupdata)
  }
  colnames(Links)[1:2]<-c('source','target')
  colnames(NodesNew)[1]<-c('id')
  NodesNew <-NodesNew[,c('id', 'S/N','Intensity','RA','module',	'degree','module_size','Formula')]
  Netlisttmp <- list(Links,NodesNew)
  names(Netlisttmp) <- c('Links','Nodes')
  return(Netlisttmp)
}




