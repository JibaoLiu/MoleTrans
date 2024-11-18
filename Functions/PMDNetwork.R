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

f_calc_masses_dif_betw_sam <- function(DF1,DF2){
  if(any(colnames(DF1)%in%'id')){
    DF1$`Measured m/z` <- DF1$id
    DF2$`Measured m/z` <- DF2$id
  }
  
  diff_matrix <- outer(unique(DF1$`Measured m/z`),unique(DF2$`Measured m/z`),'-')
  diff_matrix <- abs(diff_matrix)
  #diff_matrix[!upper.tri(diff_matrix,diag = TRUE)] <- 0
  diff_DF <-  as.data.frame(diff_matrix)
  colnames(diff_DF)<- unique(DF2$`Measured m/z`)
  diff_DF$Mass1 <- unique(DF1$`Measured m/z`)
  diff_DF <- unique(diff_DF)
  pmd_result <- gather(diff_DF,Mass2,PMD,-ncol(diff_DF))
  pmd_result <- setDT(pmd_result)
  pmd_result <- pmd_result[PMD>1e-3]
  pmd_result$Mass2 <- as.numeric(pmd_result$Mass2)
  return(pmd_result)
}

FindPMD <- function(peaklist,PMDDF,mserror=2){
  pmdlist <- f_calc_masses_dif(peaklist)
  # sort pmdlist by mass
  setkey(pmdlist,PMD)
  # assign pmdlist masses
  s <- pmdlist$PMD
  n <- length(s)
  
  # assign pmd masses
  PMDDF <- setDT(PMDDF)
  setkey(PMDDF)
  d <- as.numeric(PMDDF$PMD)
  # ranges
  r.d <- range(d)
  r.s <- range(s)
  if (r.s[1] < r.d[1]){
    pmdlist <- pmdlist[PMD>=(r.d[1]-0.01),]
    # assign peaklist masses
    s <- pmdlist$PMD
    n <- length(s)
    # ranges
    r.s <- range(s)
    
  }
  
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
    tmp <- pmdlist%>% mutate(DiffPMD=abs(PMD-d[i])*1000000/PMD)%>% mutate(DiffPMD1=abs(PMD-d[i]))
    if(mserror<0.005){
      tmp <- filter(tmp,DiffPMD1<= mserror)
    }else{
      tmp <- filter(tmp,DiffPMD<= mserror)
    }
    
    #tmp <- filter(tmp,DiffPMD<=mserror)
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

FindPMD_betw_samp <- function(DF1,DF2,PMDDF,multiDistance=FALSE){
  pmdlist <- f_calc_masses_dif_betw_sam(DF1,DF2)
  # sort pmdlist by mass
  setkey(pmdlist,PMD)
  # assign pmdlist masses
  s <- pmdlist$PMD
  n <- length(s)
  
  # assign pmd masses
  PMDDF <- setDT(PMDDF)
  setkey(PMDDF)
  d <- as.numeric(PMDDF$PMD)
  if(multiDistance==FALSE){
    # ranges
    r.d <- range(d)
    r.s <- range(s)
    if (r.s[1] < r.d[1]){
      pmdlist <- pmdlist[PMD>=(r.d[1]-0.01),]
      # assign peaklist masses
      s <- pmdlist$PMD
      n <- length(s)
      # ranges
      r.s <- range(s)
      
    }
    
    if (r.s[2] > r.d[2]){
      pmdlist <- pmdlist[PMD<=(r.d[2]+0.01),]
      # assign peaklist masses
      s <- pmdlist$PMD
      n <- length(s)
      # ranges
      r.s <- range(s)
      
    }
  }
  pmd_find <- tibble(Mass1=NA,Mass2=NA,PMD_cal=NA,PMD=NA)
  for (i in 1:length(d)) {
    if(multiDistance==FALSE){
      tmp <- pmdlist%>% mutate(DiffPMD=abs(PMD-d[i])*1000000/PMD) 
      tmp <- filter(tmp,DiffPMD<=0.01)
    }else{
      tmp <- pmdlist%>% mutate(DiffPMD=abs(1-PMD/(d[i]*round(PMD/d[i])))*1000000)
      tmp <- filter(tmp,DiffPMD<=0.01)
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

ToLinkNode <- function(pmd_result,featureDF){
  
  node_1 <- pmd_result[,c("Mass1","Formula1")]
  colnames(node_1) <- c("Mass","Formula")
  node_2 <- pmd_result[,c("Mass2","Formula2")] 
  colnames(node_2) <- c("Mass","Formula")
  node <- unique(rbind(node_1,node_2))
  colnames(node) <- c("id","Label")
  if(any(colnames(featureDF)%in%'Compund_Name')&any(colnames(featureDF)%in%'SMILES')){
    Feature <- featureDF[,c("Measured m/z","Formula","RA","O/C","N/C","AImod","NOSC","DBE/C","Group","ReInt","groupKVD","NewChemgroup","DeltaG","MLB","Compund_Name","SMILES")]
    #用trunc直接截取，避免四舍五入遇到5的问题
    Feature$`Measured m/z`<- trunc(Feature$`Measured m/z`*10000)/10000
    colnames(Feature) <- c("id_round","Formula","RA","O/C","N/C","AImod","NOSC","DBE/C","ElementGroup","ReInt","groupKVD","NewChemgroup","DeltaG","MLB","Compund_Name","SMILES")
    
  }else{
    Feature <- featureDF[,c("Measured m/z","Formula","RA","O/C","N/C","AImod","NOSC","DBE/C","Group","ReInt","groupKVD","NewChemgroup","DeltaG","MLB")]
    #用trunc直接截取，避免四舍五入遇到5的问题
    Feature$`Measured m/z`<- trunc(Feature$`Measured m/z`*10000)/10000
    colnames(Feature) <- c("id_round","Formula","RA","O/C","N/C","AImod","NOSC","DBE/C","ElementGroup","ReInt","groupKVD","NewChemgroup","DeltaG","MLB")
  }
 
  node_tmp <- node
  node_tmp$id_round <- trunc(node_tmp$id*10000)/10000
  node_attrib <- left_join(node_tmp,Feature,by="id_round")
  # id按照mass的小数点3位匹配不够，会有重复值
  #依据元素组成和化合物类别来区分反应类别，为Link添加反应标签
  node_attrib <- unique(node_attrib)
  ElementReact <- tibble(ElementReact=NA)
  KVDReact <- tibble(KVDReact=NA)
  KVD2React <- tibble(KVD2React=NA)
  for (i in 1:nrow(pmd_result)) {
    ElementGP1=node_attrib[node_attrib$id %in% pmd_result$Mass1[i],]$ElementGroup[1]
    ElementGP2=node_attrib[node_attrib$id %in% pmd_result$Mass2[i],]$ElementGroup[1]
    ElementReactLab <- paste(ElementGP1,ElementGP2,sep = "-")
    ElementReact_tmp <- tibble(ElementReact=ElementReactLab)
    KVDGP1=node_attrib[node_attrib$id %in% pmd_result$Mass1[i],]$groupKVD[1]
    KVDGP2=node_attrib[node_attrib$id %in% pmd_result$Mass2[i],]$groupKVD[1]
    KVDReactLab <- paste(KVDGP1,KVDGP2,sep = "-")
    KVDReact_tmp <- tibble(KVDReact=KVDReactLab)
    KVD2GP1=node_attrib[node_attrib$id %in% pmd_result$Mass1[i],]$NewChemgroup[1]
    KVD2GP2=node_attrib[node_attrib$id %in% pmd_result$Mass2[i],]$NewChemgroup[1]
    KVD2ReactLab <- paste(KVD2GP1,KVD2GP2,sep = "-")
    KVD2React_tmp <- tibble(KVD2React=KVD2ReactLab)
    
    ElementReact <- rbind(ElementReact,ElementReact_tmp)
    KVDReact <-rbind(KVDReact,KVDReact_tmp)
    KVD2React <-rbind(KVD2React,KVD2React_tmp)
  }
  ElementReact <- ElementReact[2:nrow(ElementReact),]
  KVDReact <- KVDReact[2:nrow(KVDReact),]
  KVD2React <- KVD2React[2:nrow(KVD2React),]
  pmd_result <- cbind(pmd_result,ElementReact)
  pmd_result <- cbind(pmd_result,KVDReact)
  pmd_result <- cbind(pmd_result,KVD2React)
  
  link_mass <- pmd_result[,c("Mass1","Mass2","Formula1","Formula2","PMD","PMD_cal","ElementReact","KVDReact","KVD2React")]
  colnames(link_mass) <- c("source","target","Formula1","Formula2","PMD","PMD_cal","ElementReact","KVDReact","KVD2React")
  link_formula <- pmd_result[,c("Formula1","Formula2","PMD","PMD_cal","ElementReact","KVDReact","KVD2React")]
  colnames(link_formula) <- c("source","target","PMD","PMD_cal","ElementReact","KVDReact","KVD2React")
  node_attrib <- filter(node_attrib,!is.na(node_attrib$Label))
  LinKNodeLs <- list(link_mass,link_formula,node_attrib)
  names(LinKNodeLs)<-c("link_mass","link_formula","node_attrib")
  return(LinKNodeLs)
}


ToLinkNode_CrossNet <- function(pmd_result,featureDF){
  
  node_1 <- pmd_result[,c("Mass1","Formula1")]
  colnames(node_1) <- c("Mass","Formula")
  node_2 <- pmd_result[,c("Mass2","Formula2")] 
  colnames(node_2) <- c("Mass","Formula")
  node <- unique(rbind(node_1,node_2))
  colnames(node) <- c("id","Label")
  
  Feature <- featureDF[,c("Measured m/z","Formula","O/C","N/C","AImod","NOSC","DBE/C","Group","groupKVD","NewChemgroup","DeltaG","MLB")]
  #用trunc直接截取，避免四舍五入遇到5的问题
  Feature$`Measured m/z`<- trunc(Feature$`Measured m/z`*1000000)/1000000
  colnames(Feature) <- c("id_round","Formula","O/C","N/C","AImod","NOSC","DBE/C","ElementGroup","groupKVD","NewChemgroup","DeltaG","MLB")
  node_tmp <- node
  node_tmp$id_round <- trunc(node_tmp$id*1000000)/1000000
  node_attrib <- left_join(node_tmp,Feature,by="id_round")
  # id按照mass的小数点3位匹配不够，会有重复值
  #依据元素组成和化合物类别来区分反应类别，为Link添加反应标签
  node_attrib <- unique(node_attrib)
  ElementReact <- tibble(ElementReact=NA)
  KVDReact <- tibble(KVDReact=NA)
  KVD2React <- tibble(KVD2React=NA)
  for (i in 1:nrow(pmd_result)) {
    ElementGP1=node_attrib[node_attrib$id %in% pmd_result$Mass1[i],]$ElementGroup[1]
    ElementGP2=node_attrib[node_attrib$id %in% pmd_result$Mass2[i],]$ElementGroup[1]
    ElementReactLab <- paste(ElementGP1,ElementGP2,sep = "-")
    ElementReact_tmp <- tibble(ElementReact=ElementReactLab)
    KVDGP1=node_attrib[node_attrib$id %in% pmd_result$Mass1[i],]$groupKVD[1]
    KVDGP2=node_attrib[node_attrib$id %in% pmd_result$Mass2[i],]$groupKVD[1]
    KVDReactLab <- paste(KVDGP1,KVDGP2,sep = "-")
    KVDReact_tmp <- tibble(KVDReact=KVDReactLab)
    KVD2GP1=node_attrib[node_attrib$id %in% pmd_result$Mass1[i],]$NewChemgroup[1]
    KVD2GP2=node_attrib[node_attrib$id %in% pmd_result$Mass2[i],]$NewChemgroup[1]
    KVD2ReactLab <- paste(KVD2GP1,KVD2GP2,sep = "-")
    KVD2React_tmp <- tibble(KVD2React=KVD2ReactLab)
    
    ElementReact <- rbind(ElementReact,ElementReact_tmp)
    KVDReact <-rbind(KVDReact,KVDReact_tmp)
    KVD2React <-rbind(KVD2React,KVD2React_tmp)
  }
  ElementReact <- ElementReact[2:nrow(ElementReact),]
  KVDReact <- KVDReact[2:nrow(KVDReact),]
  KVD2React <- KVD2React[2:nrow(KVD2React),]
  pmd_result <- cbind(pmd_result,ElementReact)
  pmd_result <- cbind(pmd_result,KVDReact)
  pmd_result <- cbind(pmd_result,KVD2React)
  pmd_result$n <- round(pmd_result$PMD_cal/pmd_result$PMD)
  link_mass <- pmd_result[,c("Mass1","Mass2","Formula1","Formula2","PMD","PMD_cal","n","ElementReact","KVDReact","KVD2React")]
  colnames(link_mass) <- c("source","target","Formula1","Formula2","PMD","PMD_cal","n","ElementReact","KVDReact","KVD2React")
  link_formula <- pmd_result[,c("Formula1","Formula2","PMD","PMD_cal","n","ElementReact","KVDReact","KVD2React")]
  colnames(link_formula) <- c("source","target","PMD","PMD_cal","n","ElementReact","KVDReact","KVD2React")
  node_attrib <- unique(node_attrib)
  #node_attrib <- filter(node_attrib,!duplicated(node_attrib$id))
  if(any(duplicated(node_attrib$id_round))){
    node_attrib <- node_attrib[-which(duplicated(node_attrib$id_round)),]
  }
  node_attrib <- filter(node_attrib,!is.na(node_attrib$Label))
  LinKNodeLs <- list(link_mass,link_formula,node_attrib)
  names(LinKNodeLs)<-c("link_mass","link_formula","node_attrib")
  return(LinKNodeLs)
}


MassToNetworkFile <- function(FtmsDF,PMD,mserror=2){

  if(str_detect(str_flatten(colnames(FtmsDF)),"Measured Mass")){
    colnames(FtmsDF)[1]<-'Measured m/z'
  }
  if(!str_detect(str_flatten(colnames(FtmsDF)),"Neutral m/z")){
    FtmsDF$`Neutral m/z` <- FtmsDF$`Measured m/z`+1.0072765
  }
  if(str_detect(str_flatten(colnames(FtmsDF)),"%")){
    FtmsDF$RA <- FtmsDF$`RA(%)`
  }
  if(!str_detect(str_flatten(colnames(FtmsDF)),"15N")){
    FtmsDF$`15N` <- 0
  }
  if(!str_detect(str_flatten(colnames(FtmsDF)),"13C")){
    FtmsDF$`13C` <- FtmsDF$`C13`
  }
  if(!str_detect(str_flatten(colnames(FtmsDF)),"2H")){
    FtmsDF$`2H` <- 0
  }
  if(!str_detect(str_flatten(colnames(FtmsDF)),"18O")){
    FtmsDF$`18O` <- 0
  }
  if(!str_detect(str_flatten(colnames(FtmsDF)),"34S")){
    FtmsDF$`34S` <- 0
  }
  if(!str_detect(str_flatten(colnames(FtmsDF)),"37Cl")){
    FtmsDF$`37Cl` <- 0
  }
  if(!str_detect(str_flatten(colnames(FtmsDF)),"81Br")){
    FtmsDF$`81Br` <- 0
  }
  if (!str_detect(str_flatten(colnames(FtmsDF)),"Formula")){
    if(str_detect(str_flatten(colnames(FtmsDF)),"Cl")){
      FtmsDF$Formula <- paste('C',FtmsDF$C+FtmsDF$`13C`,'H',FtmsDF$H+FtmsDF$`2H`,'O',FtmsDF$O+FtmsDF$`18O`,'N',FtmsDF$N+FtmsDF$`15N`,'P',FtmsDF$P,'S',FtmsDF$S+FtmsDF$`34S`,'Cl',FtmsDF$Cl+FtmsDF$`37Cl`,'Br',FtmsDF$Br+FtmsDF$`81Br`,'I',FtmsDF$I,sep = '')%>%str_remove_all('N0')%>% str_remove_all('P0') %>% str_remove_all('S0')%>% str_remove_all('Cl0')%>% str_remove_all('Br0')%>% str_remove_all('I0')
    }else{
      FtmsDF$Formula <- paste('C',FtmsDF$C+FtmsDF$`13C`,'H',FtmsDF$H+FtmsDF$`2H`,'O',FtmsDF$O+FtmsDF$`18O`,'N',FtmsDF$N+FtmsDF$`15N`,'P',FtmsDF$P,'S',FtmsDF$S+FtmsDF$`34S`,sep = '')%>%str_remove_all('N0')%>% str_remove_all('P0') %>% str_remove_all('S0')
    }
  }
  if (!str_detect(str_flatten(colnames(FtmsDF)),"ReInt")){
    FtmsDF$ReInt <- FtmsDF$Intensity*100/sum(FtmsDF$Intensity)
  }
  if (!str_detect(str_flatten(colnames(FtmsDF)),"groupKVD")){
    FtmsDF <- KVD(FtmsDF)
  }

  if (!str_detect(str_flatten(colnames(FtmsDF)),"Group")){
    FtmsDF <- ElementGroup(FtmsDF)
  }
  FtmsDF$`Measured m/z`<- trunc(FtmsDF$`Measured m/z`*1000000)/1000000
 
  PMD_result <- FindPMD(FtmsDF,PMD,mserror=mserror)
  Kegg_formula <- FtmsDF[,c("Measured m/z","Formula")]#Input ms，formula
  Kegg_formula$Mass <- Kegg_formula$`Measured m/z`
  Kegg_formula_tmp <- Kegg_formula[,c("Mass","Formula")]
  
  colnames(Kegg_formula_tmp)<- c(colnames(PMD_result)[1],"Formula1")
  PMD_result <- left_join(PMD_result,Kegg_formula_tmp,by="Mass1")
  colnames(Kegg_formula_tmp)<- c(colnames(PMD_result)[2],"Formula2")
  PMD_result <- left_join(PMD_result,Kegg_formula_tmp,by="Mass2")
  PMD_result <- unique(PMD_result)
  
  LinkNodeLs <- ToLinkNode(PMD_result,FtmsDF)
  return(LinkNodeLs)
}


MassToNetworkFile_betw_sam <- function(FtmsDF1,FtmsDF2,PMD,multiDistance=FALSE){
  
  if(str_detect(str_flatten(colnames(FtmsDF1)),"Measured Mass")){
    colnames(FtmsDF1)[1]<-'Measured m/z'
    colnames(FtmsDF2)[1]<-'Measured m/z'
  }
  if(!str_detect(str_flatten(colnames(FtmsDF1)),"Neutral m/z")){
    FtmsDF1$`Neutral m/z` <- FtmsDF1$`Measured m/z`+1.0072765
    FtmsDF2$`Neutral m/z` <- FtmsDF2$`Measured m/z`+1.0072765
  }
  if(str_detect(str_flatten(colnames(FtmsDF1)),"%")){
    FtmsDF1$RA <- FtmsDF1$`RA(%)`
    FtmsDF2$RA <- FtmsDF2$`RA(%)`
  }
  if(!str_detect(str_flatten(colnames(FtmsDF1)),"15N")){
    FtmsDF1$`15N` <- 0
    FtmsDF2$`15N` <- 0
  }
  if(!str_detect(str_flatten(colnames(FtmsDF1)),"13C")){
    FtmsDF1$`13C` <- FtmsDF1$`C13`
    FtmsDF2$`13C` <- FtmsDF2$`C13`
  }
  if(!str_detect(str_flatten(colnames(FtmsDF1)),"2H")){
    FtmsDF1$`2H` <- 0
    FtmsDF2$`2H` <- 0
  }
  if(!str_detect(str_flatten(colnames(FtmsDF1)),"18O")){
    FtmsDF1$`18O` <- 0
    FtmsDF2$`18O` <- 0
  }
  if(!str_detect(str_flatten(colnames(FtmsDF1)),"34S")){
    FtmsDF1$`34S` <- 0
    FtmsDF2$`34S` <- 0
  }
  if(!str_detect(str_flatten(colnames(FtmsDF1)),"37Cl")){
    FtmsDF1$`37Cl` <- 0
    FtmsDF2$`37Cl` <- 0
  }
  if(!str_detect(str_flatten(colnames(FtmsDF1)),"81Br")){
    FtmsDF1$`81Br` <- 0
    FtmsDF2$`81Br` <- 0
  }
  if (!str_detect(str_flatten(colnames(FtmsDF1)),"Formula")|!str_detect(str_flatten(colnames(FtmsDF1)),"accuratemass")){
    if(any(colnames(FtmsDF1)%in%"Br")){
      if(any(colnames(FtmsDF1)%in%"D")){
        FtmsDF1 <- FtmsDF1%>%filter(`13C`==0)%>%filter(`2H`==0)%>%filter(`18O`==0)%>%filter(`15N`==0)%>%filter(`34S`==0)%>%filter(`37Cl`==0)
        FtmsDF2 <- FtmsDF2%>%filter(`13C`==0)%>%filter(`2H`==0)%>%filter(`18O`==0)%>%filter(`15N`==0)%>%filter(`34S`==0)%>%filter(`37Cl`==0)
        FtmsDF1 <- setDT(FtmsDF1)
        FtmsDF2 <- setDT(FtmsDF2)
        FtmsDF1 <- FtmsDF1[,accuratemass:=C*12.0 + 
                             H * 1.00782503223 + 
                             D * 2.014102 + 
                             N * 14.00307400443 + 
                             O * 15.99491461957 + 
                             P * 30.97376199842 +
                             S * 31.9720711744 +
                             Cl * 34.968852682 +
                             Br * 78.9183376 +
                             `81Br` * 80.916291 +
                             I * 126.9044719]
        FtmsDF2 <- FtmsDF2[,accuratemass:=C*12.0 + 
                             H * 1.00782503223 + 
                             D * 2.014102 + 
                             N * 14.00307400443 + 
                             O * 15.99491461957 + 
                             P * 30.97376199842 +
                             S * 31.9720711744 +
                             Cl * 34.968852682 +
                             Br * 78.9183376 +
                             `81Br` * 80.916291 +
                             I * 126.9044719]
        
      }else{
        FtmsDF1$Formula <- paste('C',FtmsDF1$C+FtmsDF1$`13C`,'H',FtmsDF1$H+FtmsDF1$`2H`,'O',FtmsDF1$O+FtmsDF1$`18O`,'N',FtmsDF1$N+FtmsDF1$`15N`,'P',FtmsDF1$P,'S',FtmsDF1$S+FtmsDF1$`34S`,'Cl',FtmsDF1$Cl+FtmsDF1$`37Cl`,'Br',FtmsDF1$Br+FtmsDF1$`81Br`,'I',FtmsDF1$I,sep = '')%>%str_remove_all('N0')%>% str_remove_all('P0') %>% str_remove_all('S0')%>% str_remove_all('Cl0')%>% str_remove_all('Br0')%>% str_remove_all('I0')
        FtmsDF2$Formula <- paste('C',FtmsDF2$C+FtmsDF2$`13C`,'H',FtmsDF2$H+FtmsDF2$`2H`,'O',FtmsDF2$O+FtmsDF2$`18O`,'N',FtmsDF2$N+FtmsDF2$`15N`,'P',FtmsDF2$P,'S',FtmsDF2$S+FtmsDF2$`34S`,'Cl',FtmsDF2$Cl+FtmsDF2$`37Cl`,'Br',FtmsDF2$Br+FtmsDF2$`81Br`,'I',FtmsDF2$I,sep = '')%>%str_remove_all('N0')%>% str_remove_all('P0') %>% str_remove_all('S0')%>% str_remove_all('Cl0')%>% str_remove_all('Br0')%>% str_remove_all('I0')
        FtmsDF1 <- FtmsDF1%>%filter(`13C`==0)%>%filter(`2H`==0)%>%filter(`18O`==0)%>%filter(`15N`==0)%>%filter(`34S`==0)%>%filter(`37Cl`==0)%>%filter(`81Br`==0)
        FtmsDF2 <- FtmsDF2%>%filter(`13C`==0)%>%filter(`2H`==0)%>%filter(`18O`==0)%>%filter(`15N`==0)%>%filter(`34S`==0)%>%filter(`37Cl`==0)%>%filter(`81Br`==0)
        
        FtmsDF1 <- setDT(FtmsDF1)
        FtmsDF2 <- setDT(FtmsDF2)
        FtmsDF1 <- FtmsDF1[,accuratemass:=C*12.0 + 
                             H * 1.00782503223 + 
                             N * 14.00307400443 + 
                             O * 15.99491461957 + 
                             P * 30.97376199842 +
                             S * 31.9720711744 +
                             Cl * 34.968852682 +
                             Br * 78.9183376 +
                             `81Br` * 80.916291 +
                             I * 126.9044719]
        FtmsDF2 <- FtmsDF2[,accuratemass:=C*12.0 + 
                             H * 1.00782503223 + 
                             N * 14.00307400443 + 
                             O * 15.99491461957 + 
                             P * 30.97376199842 +
                             S * 31.9720711744 +
                             Cl * 34.968852682 +
                             Br * 78.9183376 +
                             `81Br` * 80.916291 +
                             I * 126.9044719]
      }
    }else{
      FtmsDF1$Formula <- paste('C',FtmsDF1$C+FtmsDF1$`13C`,'H',FtmsDF1$H+FtmsDF1$`2H`,'O',FtmsDF1$O+FtmsDF1$`18O`,'N',FtmsDF1$N+FtmsDF1$`15N`,'P',FtmsDF1$P,'S',FtmsDF1$S+FtmsDF1$`34S`,sep = '')%>%str_remove_all('N0')%>% str_remove_all('P0') %>% str_remove_all('S0')
      FtmsDF2$Formula <- paste('C',FtmsDF2$C+FtmsDF2$`13C`,'H',FtmsDF2$H+FtmsDF2$`2H`,'O',FtmsDF2$O+FtmsDF2$`18O`,'N',FtmsDF2$N+FtmsDF2$`15N`,'P',FtmsDF2$P,'S',FtmsDF2$S+FtmsDF2$`34S`,sep = '')%>%str_remove_all('N0')%>% str_remove_all('P0') %>% str_remove_all('S0')
      FtmsDF1 <- FtmsDF1%>%filter(`13C`==0)%>%filter(`2H`==0)%>%filter(`18O`==0)%>%filter(`15N`==0)%>%filter(`34S`==0)
      FtmsDF2 <- FtmsDF2%>%filter(`13C`==0)%>%filter(`2H`==0)%>%filter(`18O`==0)%>%filter(`15N`==0)%>%filter(`34S`==0)
      FtmsDF1 <- setDT(FtmsDF1)
      FtmsDF2 <- setDT(FtmsDF2)
      FtmsDF1 <- FtmsDF1[,accuratemass:=C*12.0 + 
                           H * 1.00782503223 + 
                           N * 14.00307400443 + 
                           O * 15.99491461957 + 
                           P * 30.97376199842 +
                           S * 31.9720711744]
      FtmsDF2 <- FtmsDF2[,accuratemass:=C*12.0 + 
                           H * 1.00782503223 + 
                           N * 14.00307400443 + 
                           O * 15.99491461957 + 
                           P * 30.97376199842 +
                           S * 31.9720711744]
      
    }
  }
  if (!str_detect(str_flatten(colnames(FtmsDF1)),"ReInt")){
    FtmsDF1$ReInt <- FtmsDF1$Intensity*100/sum(FtmsDF1$Intensity)
    FtmsDF2$ReInt <- FtmsDF2$Intensity*100/sum(FtmsDF2$Intensity)
  }
  if (!str_detect(str_flatten(colnames(FtmsDF1)),"groupKVD")){
    FtmsDF1 <- KVD(FtmsDF1)
    FtmsDF2 <- KVD(FtmsDF2)
  }
  
  if (!str_detect(str_flatten(colnames(FtmsDF1)),"Group")){
    FtmsDF1 <- ElementGroup(FtmsDF1)
    FtmsDF2 <- ElementGroup(FtmsDF2)
  }
  FtmsDF1$`Measured m/z`<- trunc(FtmsDF1$accuratemass*1000000)/1000000
  FtmsDF2$`Measured m/z`<- trunc(FtmsDF2$accuratemass*1000000)/1000000
  FtmsDF1 <- unique(FtmsDF1[,c("Measured m/z","Formula","O/C","N/C","AImod","NOSC","DBE/C","Group","groupKVD","NewChemgroup","DeltaG","MLB")])
  FtmsDF2 <- unique(FtmsDF2[,c("Measured m/z","Formula","O/C","N/C","AImod","NOSC","DBE/C","Group","groupKVD","NewChemgroup","DeltaG","MLB")])
  FtmsDF1 <- filter(FtmsDF1,!duplicated(FtmsDF1$`Measured m/z`))
  FtmsDF2 <- filter(FtmsDF2,!duplicated(FtmsDF2$`Measured m/z`))
  PMD_result <- FindPMD_betw_samp(FtmsDF1,FtmsDF2,PMD,multiDistance=multiDistance)
  formula1 <- unique(FtmsDF1[,c("Measured m/z","Formula")])#Input ms，formula
  formula2 <- unique(FtmsDF2[,c("Measured m/z","Formula")])
  colnames(formula1) <- c("Mass1","Formula1")
  colnames(formula2) <- c("Mass2","Formula2")
  PMD_result <- left_join(PMD_result,formula1,by="Mass1")
  PMD_result <- left_join(PMD_result,formula2,by="Mass2")
  PMD_result$label1 <- abs(PMD_result$Mass1-PMD_result$Mass2)+log(PMD_result$Mass1)*log(PMD_result$Mass2)
  PMD_result <- PMD_result[!duplicated(PMD_result[,ncol(PMD_result)]),]
  PMD_result <- unique(PMD_result)[,-ncol(PMD_result)]
  FtmsDF <- unique(rbind(FtmsDF1,FtmsDF2))
  LinkNodeLs <- ToLinkNode_CrossNet(PMD_result,FtmsDF)
  return(LinkNodeLs)
}

PMDFormulaConsensusFilter <- function(netlink,netnode,PMD){
  Formula1 <- FormuTrans(select(netlink,c('Formula1')))
  Formula2 <- FormuTrans(select(netlink,c('Formula2')))
  Formula1$rowname <- rownames(Formula1)
  Formula2$rowname <- rownames(Formula2)
  new_df <- full_join(Formula1,Formula2,by="rowname")
  new_df$C.p <- new_df$C.y-new_df$C.x
  new_df$H.p <- new_df$H.y-new_df$H.x
  new_df$O.p <- new_df$O.y-new_df$O.x
  new_df$N.p <- new_df$N.y-new_df$N.x
  new_df$P.p <- new_df$P.y-new_df$P.x
  new_df$S.p <- new_df$S.y-new_df$S.x
  colnames(new_df)<-str_remove_all(colnames(new_df),'.p')
  new_df$TransGroup <- paste('C',new_df$C,'H',new_df$H,'O',new_df$O,'N',new_df$N,'P',new_df$P,'S',new_df$S,sep = '') %>% str_remove_all('C0')%>% str_remove_all('H0') %>% str_remove_all('O0') %>% str_remove_all('N0')%>% str_remove_all('P0') %>% str_remove_all('S0')
  new_df$PMDmass <-abs(new_df$C*12+new_df$H*1.007825+new_df$O*15.994915+new_df$N*14.003074+new_df$S*31.972071+new_df$P*30.973762)
  PMDinfo <-  PMD
  colnames(PMDinfo)[which(colnames(PMDinfo) == "Group")] <- "Formula"
  PMDinfo1 <-  FormuTrans(select(PMDinfo,c('Formula')))
  PMDinfo1$PMDmass_ref <- abs(PMDinfo1$C*12+PMDinfo1$H*1.007825+PMDinfo1$O*15.994915+PMDinfo1$N*14.003074+PMDinfo1$S*31.972071+PMDinfo1$P*30.973762)
  netlink1 <- cbind(netlink,select(new_df,c('TransGroup','PMDmass')))
  netlink1 <- filter(netlink1, netlink1$PMDmass %in% PMDinfo1$PMDmass_ref)
  netnode1 <- filter(netnode,netnode$id %in% c(netlink1$source,netlink1$target))
  LinKNodeLs <- list(netlink1,netlink,netnode1)
  names(LinKNodeLs)<-c("link_mass","link_origin","node_attrib")
  return(LinKNodeLs)
}



FindSubNetFromNode <- function(NetEdge,NetNode,SelectNode,Length){
  Selectid<-NetNode[NetNode$Label %in% SelectNode,]$id
  if(any(colnames(NetEdge) %in% "TransGroup")){
    SelectEdge <- tibble(source=NA,target=NA,Formula1=NA,Formula2=NA,PMD=NA,PMD_cal=NA,ElementReact=NA,KVDReact=NA,KVD2React=NA,TransGroup=NA,PMDmass=NA)
  }else if(any(colnames(NetEdge) %in% "n")){
    SelectEdge <- tibble(source=NA,target=NA,Formula1=NA,Formula2=NA,PMD=NA,PMD_cal=NA,n=NA,ElementReact=NA,KVDReact=NA,KVD2React=NA)
  }else{
    SelectEdge <- tibble(source=NA,target=NA,Formula1=NA,Formula2=NA,PMD=NA,PMD_cal=NA,ElementReact=NA,KVDReact=NA,KVD2React=NA)
  }
    
  SelectNode <- NetNode[NetNode$Label %in% SelectNode,]
  for (i in 1:Length) {
    SelectEdge_tmp <- NetEdge[NetEdge$source %in% Selectid|NetEdge$target %in% Selectid,]
    Selectid <- unique(c(SelectEdge_tmp$source,SelectEdge_tmp$target))
    SelectNode_tmp <- NetNode[NetNode$id %in% Selectid,]
    SelectEdge <- rbind(SelectEdge,SelectEdge_tmp)
    SelectNode <- rbind(SelectNode,SelectNode_tmp)
  }
  SelectEdge <- unique(SelectEdge[2:nrow(SelectEdge),])
  SelectNode <- unique(SelectNode)
  SelectNet <- list(SelectEdge,SelectNode)
  names(SelectNet)<-c("Edge","Node")
  return(SelectNet)
}

CrossNetPMDFreq <- function(NetLink,PMD_info){
  if(any(PMD_info$Group %in% 'C1')){
    PMD_info[PMD_info$Group %in% 'C1',]$PMD <- 12.0
  }
  if(any(PMD_info$Group %in% 'C')){
    PMD_info[PMD_info$Group %in% 'C',]$PMD <- 12.0
  }
  PMD_info$PMD1 <- trunc(PMD_info$PMD*10000)/10000
  NetLink$PMD1 <- trunc(NetLink$PMD*10000)/10000
  data_tmp <- left_join(NetLink,PMD_info,by='PMD1')
  count_pmd <- data_tmp%>%count(Group)
  count_pmd$`Frequency(%)` <- round(count_pmd$n*100/sum(count_pmd$n),digits = 2)
  return(count_pmd)
}

PMDRadarPlot <- function(dataTB,NumTop=10){
  if(nrow(dataTB)>NumTop){
    dataTB <- dataTB[order(dataTB$`Frequency(%)`,decreasing = TRUE),]
    Others <- tibble(Group=NA,n=NA,`Frequency(%)`=NA)
    Others[1,1] <- "Others"
    Others[1,2] <- sum(dataTB[(NumTop+1):nrow(dataTB),]$n)
    Others[1,3] <- sum(dataTB[(NumTop+1):nrow(dataTB),]$`Frequency(%)`)
    dataTB <- rbind(dataTB[1:NumTop,],Others)
  }
  
  rownames(dataTB) <- dataTB$Group
  #dataTB$Min <- 0.0
  #dataTB$Max <- 100.0
  dataTBtmp <- select(dataTB,'Frequency(%)')
  dataTBtmp <- as.data.frame(sapply(dataTBtmp,as.numeric))
  dataTBtmp$`Frequency` <-dataTBtmp$`Frequency(%)`/100
  rownames(dataTBtmp) <-  dataTB$Group
  dataTBtmp <- select(dataTBtmp,'Frequency')
  dataTBtmp <- as.data.frame(t(dataTBtmp))
  dataTBtmp <- rownames_to_column(dataTBtmp)
  p <- ggradar(dataTBtmp,values.radar = c("0%","50%","100%"),base.size = 12,grid.label.size = 5,
          group.colours="#00AFBB",group.point.size = 0.8, group.line.width = 0.5,axis.label.size = 4,
          fill=TRUE,fill.alpha = 0.5)+theme(plot.margin = ggplot2::margin(2,2,2,2,'mm'))
  return(p)
}
DiffNetwork <- function(FtmsDF1,FtmsDF2, PMD,masserror=0.00001){
  Net2 <- MassToNetworkFile_betw_sam(FtmsDF1,FtmsDF2,PMD)
  Net1 <- MassToNetworkFile_betw_sam(FtmsDF1,FtmsDF1,PMD)
  Net1Link <- Net1$link_mass
  Net2Link <- Net2$link_mass
  #Net1Link <- Net1Link[,3:9]
  #Net2Link <- Net2Link[,3:9]
  colnames(Net1Link)[1:4]<-c('mz1','mz2','source','target')
  colnames(Net2Link)[1:4]<-c('mz1','mz2','source','target')
  Net1Link <- Net1Link[,c("source","target", "mz1", "mz2","PMD","PMD_cal","ElementReact","KVDReact","KVD2React")]
  Net2Link <- Net2Link[,c("source","target", "mz1", "mz2","PMD","PMD_cal","ElementReact","KVDReact","KVD2React")]
  Graph1 <- graph_from_data_frame(Net1Link,directed = FALSE)
  Graph2 <- graph_from_data_frame(Net2Link,directed = FALSE)
  #Graph1 <- graph_from_data_frame(Net1$link_mass,directed = FALSE)
  #Graph2 <- graph_from_data_frame(Net2$link_mass,directed = FALSE)
  GraphDiff <- Graph2 %m% Graph1
  info_tmp <- edge.attributes(GraphDiff)
  link_info <- as.data.frame(info_tmp[[1]])
  colnames(link_info) <- colnames(Net2$link_mass)[3]
  for (i in 2:length(info_tmp)) {
    link_info <- cbind(link_info,as.data.frame(info_tmp[[i]]))
    colnames(link_info)[i] <- colnames(Net2$link_mass)[2+i]
  }
  link_mass <- as_edgelist(GraphDiff) %>% as.data.frame()%>% cbind(link_info)
  colnames(link_mass) <- colnames(Net2Link)
  colnames(link_mass)[1:4] <- c('Formula1','Formula2','source','target')
  link_mass <- link_mass[c("source","target", 'Formula1','Formula2',"PMD","PMD_cal","ElementReact","KVDReact","KVD2React")]
  node_id <- unique(c(link_mass$source,link_mass$target))
  link_formula <- link_mass[,3:ncol(link_mass)]
  colnames(link_formula)[1:2]<-c('source','target')
  node_attrib <- filter(Net2$node_attrib,id %in% node_id)
  res <- list(link_mass,link_formula,node_attrib)
  names(res)<-c('link_mass','link_formula','node_attrib')
  return(res)
}
DiffNetwork2 <- function(FtmsDF1,FtmsDF2, PMD,masserror=0.00001){
  FtmsDF1_2 <- rbind(FtmsDF1,FtmsDF2) 
  Net2 <- MassToNetworkFile_betw_sam(FtmsDF1_2,FtmsDF2,PMD)
  Net1 <- MassToNetworkFile_betw_sam(FtmsDF1,FtmsDF1,PMD)
  
  Net1Link <- Net1$link_mass
  Net2Link <- Net2$link_mass
  #Net1Link <- Net1Link[,3:9]
  #Net2Link <- Net2Link[,3:9]
  colnames(Net1Link)[1:4]<-c('mz1','mz2','source','target')
  colnames(Net2Link)[1:4]<-c('mz1','mz2','source','target')
  Net1Link <- Net1Link[,c("source","target", "mz1", "mz2","PMD","PMD_cal","ElementReact","KVDReact","KVD2React")]
  Net2Link <- Net2Link[,c("source","target", "mz1", "mz2","PMD","PMD_cal","ElementReact","KVDReact","KVD2React")]
  Graph1 <- graph_from_data_frame(Net1Link,directed = FALSE)
  Graph2 <- graph_from_data_frame(Net2Link,directed = FALSE)
  #Graph1 <- graph_from_data_frame(Net1$link_mass,directed = FALSE)
  #Graph2 <- graph_from_data_frame(Net2$link_mass,directed = FALSE)
  GraphDiff <- Graph2 %m% Graph1
  info_tmp <- edge.attributes(GraphDiff)
  link_info <- as.data.frame(info_tmp[[1]])
  colnames(link_info) <- colnames(Net2$link_mass)[3]
  for (i in 2:length(info_tmp)) {
    link_info <- cbind(link_info,as.data.frame(info_tmp[[i]]))
    colnames(link_info)[i] <- colnames(Net2$link_mass)[2+i]
  }
  link_mass <- as_edgelist(GraphDiff) %>% as.data.frame()%>% cbind(link_info)
  colnames(link_mass) <- colnames(Net2Link)
  colnames(link_mass)[1:4] <- c('Formula1','Formula2','source','target')
  link_mass <- link_mass[c("source","target", 'Formula1','Formula2',"PMD","PMD_cal","ElementReact","KVDReact","KVD2React")]
  node_id <- unique(c(link_mass$source,link_mass$target))
  link_formula <- link_mass[,3:ncol(link_mass)]
  colnames(link_formula)[1:2]<-c('source','target')
  node_attrib <- filter(Net2$node_attrib,id %in% node_id)
  res <- list(link_mass,link_formula,node_attrib)
  names(res)<-c('link_mass','link_formula','node_attrib')
  return(res)
}

plot_molecule <- function(molecule, name = NULL, smile=NULL,sma = NULL, ...){
  #' molecule an object as returned by rcdk::load.molecules or rcdk::parse.smiles()
  #' name a character for the name of the molecule, 
  #' sma a character witht the smarts string as passed onto get.depictor()
  #' ... other arguments for get.depictor()
  
  # Image aesthetics 
  dep <- get.depictor(
    width = 1000, height = 1000,
    zoom = 7, sma = sma, ...
  )
  molecule_sdf <- view.image.2d(molecule[[1]], depictor = dep)
  
  
  ## Remove extra margins around the molecule
  par(mar=c(0,0,0,0))
  plot(NA, 
       xlim=c(1, 10), ylim=c(1, 10), 
       # Remove the black bounding boxes around the molecule
       axes = F)
  rasterImage(molecule_sdf, 1,1, 10,10)
  # Annotate the molecule
  text(x = 5.5, y = 1.1,  smile)
  text(x= 5.5, y = 10, name)
}


