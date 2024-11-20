dataPrecess <- function(data,Deu_DT=FALSE){
  if(str_detect(str_flatten(colnames(data)),"Measured Mass")){
    colnames(data)[1]<-'Measured m/z'
  }
  if(str_detect(str_flatten(colnames(data)),"%")){
    data$RA <- data$`RA(%)`
  }
  if(!str_detect(str_flatten(colnames(data)),"Neutral m/z")){
    data$`Neutral m/z` <- data$`Measured m/z`+1.0072765
  }
  if(!str_detect(str_flatten(colnames(data)),"15N")){
    data$`15N` <- 0
  }
  if(!str_detect(str_flatten(colnames(data)),"13C")){
    data$`13C` <- data$`C13`
  }

  if(!any(colnames(data) %in% "2H")){
    data$`2H` <- 0
  }
  if(!any(colnames(data) %in% "P")){
    data$P <- 0
  }
  if(!any(colnames(data) %in% "S/C")){
    data <- mutate(data,`S/C`=S/C)
  }
  if(!str_detect(str_flatten(colnames(data)),"18O")){
    data$`18O` <- 0
  }
  if(!str_detect(str_flatten(colnames(data)),"33S")){
    data$`33S` <- 0
  }
  if(!str_detect(str_flatten(colnames(data)),"34S")){
    data$`34S` <- 0
  }
  if(!str_detect(str_flatten(colnames(data)),"37Cl")){
    data$`37Cl` <- 0
  }
  if(!str_detect(str_flatten(colnames(data)),"81Br")){
    data$`81Br` <- 0
  }
  if(Deu_DT=="deu_dt"){
    data$D <- data$`2H`
    data$`2H` <- 0
    data$`H/C` <- (data$H+data$D)/(data$C+data$`13C`)
    data$`O/C` <- (data$O+data$`18O`)/(data$C+data$`13C`)
  }
  
  if (str_detect(str_flatten(colnames(data)),"Formula")){
    if(any(data$Formula %in% "Addcut Formula")){
      if(str_detect(str_flatten(colnames(data)),"Cl")){
        data$Formula <- paste('C',data$C+data$`13C`,'H',data$H+data$`2H`,'O',data$O+data$`18O`,'N',data$N+data$`15N`,'P',data$P,'S',data$S+data$`33S`+data$`34S`,'Cl',data$Cl+data$`37Cl`,'Br',data$Br+data$`81Br`,'I',data$I,sep = '')%>%str_remove_all('N0')%>% str_remove_all('P0') %>% str_remove_all('S0')%>% str_remove_all('Cl0')%>% str_remove_all('Br0')%>% str_remove_all('I0')
      }else{
        data$Formula <- paste('C',data$C+data$`13C`,'H',data$H+data$`2H`,'O',data$O+data$`18O`,'N',data$N+data$`15N`,'P',data$P,'S',data$S+data$`33S`+data$`34S`,sep = '')%>%str_remove_all('N0')%>% str_remove_all('P0') %>% str_remove_all('S0')
      }
    }
  }
  if (!str_detect(str_flatten(colnames(data)),"Formula")){
    if(Deu_DT=="deu_dt"){
      data$Formula <- paste('C',data$C+data$`13C`,'H',data$H+data$`2H`,'D',data$D,'O',data$O+data$`18O`,'N',data$N+data$`15N`,'P',data$P,'S',data$S+data$`33S`+data$`34S`,'Cl',data$Cl+data$`37Cl`,'Br',data$Br+data$`81Br`,'I',data$I,sep = '')%>%str_remove_all('D0')%>%str_remove_all('N0')%>% str_remove_all('P0') %>% str_remove_all('S0')%>% str_remove_all('Cl0')%>% str_remove_all('Br0')%>% str_remove_all('I0')
    }else{
      if(str_detect(str_flatten(colnames(data)),"Cl")){
        if(!any(colnames(data) %in% "Br")){
          data$`Br` <- 0
        }
        if(!any(colnames(data) %in% "I")){
          data$`I` <- 0
        }
        
        data$Formula <- paste('C',data$C+data$`13C`,'H',data$H+data$`2H`,'O',data$O+data$`18O`,'N',data$N+data$`15N`,'P',data$P,'S',data$S+data$`33S`+data$`34S`,'Cl',data$Cl+data$`37Cl`,'Br',data$Br+data$`81Br`,'I',data$I,sep = '')%>%str_remove_all('N0')%>% str_remove_all('P0') %>% str_remove_all('S0')%>% str_remove_all('Cl0')%>% str_remove_all('Br0')%>% str_remove_all('I0')
        #data <- data %>% filter(`13C`==0) %>% filter(`2H`==0) %>% filter(`18O`==0)%>% filter(`15N`==0)%>% filter(`33S`==0)%>% filter(`34S`==0)%>% filter(`37Cl`==0)%>% filter(`81Br`==0)
      }else{
        data$Formula <- paste('C',data$C+data$`13C`,'H',data$H+data$`2H`,'O',data$O+data$`18O`,'N',data$N+data$`15N`,'P',data$P,'S',data$S+data$`33S`+data$`34S`,sep = '')%>%str_remove_all('N0')%>% str_remove_all('P0') %>% str_remove_all('S0')
        #data <- data %>% filter(`13C`==0) %>% filter(`2H`==0) %>% filter(`18O`==0)%>% filter(`15N`==0)%>% filter(`33S`==0)%>% filter(`34S`==0)
      }
    }
  }
  if(str_detect(str_flatten(colnames(data)),"Cl")){
    
    if(any(colnames(data) %in% "33S")){
      data <- data %>% filter(`13C`==0) %>% filter(`2H`==0) %>% filter(`18O`==0)%>% filter(`15N`==0)%>% filter(`33S`==0)%>% filter(`34S`==0)%>% filter(`37Cl`==0)%>% filter(`81Br`==0)
    }else{
      data <- data %>% filter(`13C`==0) %>% filter(`2H`==0) %>% filter(`18O`==0)%>% filter(`15N`==0)%>% filter(`34S`==0)%>% filter(`37Cl`==0)%>% filter(`81Br`==0)
    }
  }else{
    if(any(colnames(data) %in% "33S")){
      data <- data %>% filter(`13C`==0) %>% filter(`2H`==0) %>% filter(`18O`==0)%>% filter(`15N`==0)%>% filter(`33S`==0)%>% filter(`34S`==0)
    }else{
      data <- data %>% filter(`13C`==0) %>% filter(`2H`==0) %>% filter(`18O`==0)%>% filter(`15N`==0)%>% filter(`34S`==0)
    }
  }
  
  
  
  if (!str_detect(str_flatten(colnames(data)),"ReInt")){
    data$ReInt <- data$Intensity*100/sum(data$Intensity)
  }
  if (!str_detect(str_flatten(colnames(data)),"DeltaG")){
    data$DeltaG <- 60.3-(data$NOSC)*28.5
  }
  if (!any(colnames(data) %in% "DBE")|!any(colnames(data) %in% "DBE/C")){
    data$DBE <- (2*data$C+data$N+data$P-data$H+2)/2
    data$`DBE-O` <- data$DBE-data$O
    data$`DBE/C` <- data$DBE/data$C
  }
  
  if (!any(colnames(data) %in% "(DBE-O)/C")){
    data$`(DBE-O)/C`<-(data$DBE-data$O)/data$C
  }
  if (!any(colnames(data) %in% "(DBE-0.5O)/C")){
    data$`(DBE-0.5O)/C`<-(data$DBE-0.5*data$O)/data$C
  }
  
  if (!str_detect(str_flatten(colnames(data)),"groupKVD")){
    data <- KVD(data)
  }
  if (!str_detect(str_flatten(colnames(data)),"NewChemgroup")){
    data <- KVDNewGroup(data)
  }
  if (!str_detect(str_flatten(colnames(data)),"MLB")){
    data <- MLB(data)
  }
  if(str_detect(str_flatten(colnames(data)),"Group")&any(data$Group %in% 'Addcut Formula')){
    data<-data%>%select(-Group)
    data<- ElementGroup(data)
  }
  
  if (!str_detect(str_flatten(colnames(data)),"Group")){
    data <- ElementGroup(data)
  }
  if (str_detect(str_flatten(colnames(data)),"Formula")){
    if(!(any(str_detect(data$Formula,"Cl"))|any(str_detect(data$Formula,"Br"))|any(str_detect(data$Formula,"I")))){
      data <- as.data.frame(data)
      data1 <- data['Formula']
      colnames(data1)[1]<-'MolForm'
      df <- main_run(data1['MolForm'])
      data$lambda <-round(df$lambda,digits = 4)
    }
  }
  if (!str_detect(str_flatten(colnames(data)),"Xc")){
    data <- Xc(data,m=0.5,n=0.5)
  }
  data$`Measured m/z`<- trunc(data$`Measured m/z`*1000000)/1000000
  if(any(duplicated(data$`Measured m/z`))){
    data <- filter(data,!duplicated(data$`Measured m/z`))
  }
  data<-setDT(data)
  setkey(data,`Measured m/z`)
  data<-setDF(data)
  return(data)
}

dataPrecess_Formularity <- function(data,Deu_DT=FALSE){
  if(str_detect(str_flatten(colnames(data)),"Mass")){
    colnames(data)[1]<-'Measured m/z'
  }
  data <- filter(data,!C==0)
  if(str_detect(str_flatten(colnames(data)),"%")){
    data$RA <- data$`RA(%)`
  }
  if(!str_detect(str_flatten(colnames(data)),"Neutral m/z")){
    data$`Neutral m/z` <- data$`Measured m/z`+1.0072765
  }
  if(!str_detect(str_flatten(colnames(data)),"15N")){
    data$`15N` <- 0
  }
  if(!str_detect(str_flatten(colnames(data)),"13C")){
    data$`13C` <- data$`C13`
  }
  
  if(!any(colnames(data) %in% "2H")){
    data$`2H` <- 0
  }
  if(!any(colnames(data) %in% "P")){
    data$P <- 0
  }
  if(!any(colnames(data) %in% "S/C")){
    data <- mutate(data,`S/C`=S/C)
  }
  if(!str_detect(str_flatten(colnames(data)),"18O")){
    data$`18O` <- 0
  }
  if(!str_detect(str_flatten(colnames(data)),"33S")){
    data$`33S` <- 0
  }
  if(!str_detect(str_flatten(colnames(data)),"34S")){
    data$`34S` <- 0
  }
  if(!str_detect(str_flatten(colnames(data)),"37Cl")){
    data$`37Cl` <- 0
  }
  if(!str_detect(str_flatten(colnames(data)),"81Br")){
    data$`81Br` <- 0
  }
  
  if(Deu_DT=="deu_dt"){
    data$D <- data$`2H`
    data$`2H` <- 0
    data$`H/C` <- (data$H+data$D)/(data$C+data$`13C`)
    data$`O/C` <- (data$O+data$`18O`)/(data$C+data$`13C`)
  }
  
  if (str_detect(str_flatten(colnames(data)),"Formula")){
    if(any(data$Formula %in% "Addcut Formula")){
      if(str_detect(str_flatten(colnames(data)),"Cl")){
        data$Formula <- paste('C',data$C+data$`13C`,'H',data$H+data$`2H`,'O',data$O+data$`18O`,'N',data$N+data$`15N`,'P',data$P,'S',data$S+data$`33S`+data$`34S`,'Cl',data$Cl+data$`37Cl`,'Br',data$Br+data$`81Br`,'I',data$I,sep = '')%>%str_remove_all('N0')%>% str_remove_all('P0') %>% str_remove_all('S0')%>% str_remove_all('Cl0')%>% str_remove_all('Br0')%>% str_remove_all('I0')
      }else{
        data$Formula <- paste('C',data$C+data$`13C`,'H',data$H+data$`2H`,'O',data$O+data$`18O`,'N',data$N+data$`15N`,'P',data$P,'S',data$S+data$`33S`+data$`34S`,sep = '')%>%str_remove_all('N0')%>% str_remove_all('P0') %>% str_remove_all('S0')
      }
    }
  }
  if (!str_detect(str_flatten(colnames(data)),"Formula")){
    if(Deu_DT=="deu_dt"){
      data$Formula <- paste('C',data$C+data$`13C`,'H',data$H+data$`2H`,'D',data$D,'O',data$O+data$`18O`,'N',data$N+data$`15N`,'P',data$P,'S',data$S+data$`33S`+data$`34S`,'Cl',data$Cl+data$`37Cl`,'Br',data$Br+data$`81Br`,'I',data$I,sep = '')%>%str_remove_all('D0')%>%str_remove_all('N0')%>% str_remove_all('P0') %>% str_remove_all('S0')%>% str_remove_all('Cl0')%>% str_remove_all('Br0')%>% str_remove_all('I0')
    }else{
      if(any(colnames(data) %in% "Cl")){
        data$Formula <- paste('C',data$C+data$`13C`,'H',data$H+data$`2H`,'O',data$O+data$`18O`,'N',data$N+data$`15N`,'P',data$P,'S',data$S+data$`33S`+data$`34S`,'Cl',data$Cl+data$`37Cl`,'Br',data$Br+data$`81Br`,'I',data$I,sep = '')%>%str_remove_all('N0')%>% str_remove_all('P0') %>% str_remove_all('S0')%>% str_remove_all('Cl0')%>% str_remove_all('Br0')%>% str_remove_all('I0')
        #data <- data %>% filter(`13C`==0) %>% filter(`2H`==0) %>% filter(`18O`==0)%>% filter(`15N`==0)%>% filter(`33S`==0)%>% filter(`34S`==0)%>% filter(`37Cl`==0)%>% filter(`81Br`==0)
      }else{
        data$Formula <- paste('C',data$C+data$`13C`,'H',data$H+data$`2H`,'O',data$O+data$`18O`,'N',data$N+data$`15N`,'P',data$P,'S',data$S+data$`33S`+data$`34S`,sep = '')%>%str_remove_all('N0')%>% str_remove_all('P0') %>% str_remove_all('S0')
        #data <- data %>% filter(`13C`==0) %>% filter(`2H`==0) %>% filter(`18O`==0)%>% filter(`15N`==0)%>% filter(`33S`==0)%>% filter(`34S`==0)
      }
    }
  }
  if(any(colnames(data) %in% "Cl")){
    if(any(colnames(data) %in% "33S")){
      data <- data %>% filter(`13C`==0) %>% filter(`2H`==0) %>% filter(`18O`==0)%>% filter(`15N`==0)%>% filter(`33S`==0)%>% filter(`34S`==0)%>% filter(`37Cl`==0)
    }else{
      data <- data %>% filter(`13C`==0) %>% filter(`2H`==0) %>% filter(`18O`==0)%>% filter(`15N`==0)%>% filter(`34S`==0)%>% filter(`37Cl`==0)
    }
  }else{
    if(any(colnames(data) %in% "33S")){
      data <- data %>% filter(`13C`==0) %>% filter(`2H`==0) %>% filter(`18O`==0)%>% filter(`15N`==0)%>% filter(`33S`==0)%>% filter(`34S`==0)
    }else{
      data <- data %>% filter(`13C`==0) %>% filter(`2H`==0) %>% filter(`18O`==0)%>% filter(`15N`==0)%>% filter(`34S`==0)
    }
  }
  
  
  
  if (!str_detect(str_flatten(colnames(data)),"ReInt")){
    data$ReInt <- data$RA*100/sum(data$RA)
  }
  if (!str_detect(str_flatten(colnames(data)),"Intensity")){
    data$Intensity <- data$ReInt
  }
  if (!str_detect(str_flatten(colnames(data)),"DBE")){
    data$DBE <- (2*data$C+data$N+data$P-data$H+2)/2
    data$`DBE-O` <- data$DBE-data$O
    data$`DBE/C` <- data$DBE/data$C
  }
  if (!str_detect(str_flatten(colnames(data)),"N/C")){
    data$`N/C` <- data$N/data$C
    data$`P/C` <- data$P/data$C
  }
  if (!str_detect(str_flatten(colnames(data)),"S/N")){
    data$`S/N` <- 6
  }
  
  if (!str_detect(str_flatten(colnames(data)),"NOSC")){
    if(!any(colnames(data) %in% "Br")){
      data$NOSC <- 4-(4*data$C+data$H-2*data$O-3*data$N-2*data$S+5*data$P)/data$C
    }else if(any(colnames(data) %in% "Br")){
      data$NOSC <- 4-(4*data$C+data$H-2*data$O-3*data$N-2*data$S+5*data$P-data$Cl-data$Br-data$`81Br`-data$I)/data$C
    }
  }
  
  
  if (!str_detect(str_flatten(colnames(data)),"AImod")){
    if(!any(colnames(data) %in% "Br")){
      data$AImod <- (1+data$C-0.5*data$O-data$S-0.5*(data$N+data$P+data$H))/(data$C-0.5*data$O-data$N-data$S-data$P)
      data$DBEAI <- (1+data$C-0.5*data$O-data$S-0.5*(data$N+data$P+data$H))
      data$CAI <-(data$C-0.5*data$O-data$N-data$S-data$P)
      datatmp1 <- filter(data,DBEAI<=0)
      datatmp2 <- filter(data,CAI<=0)
      datatmp3 <- filter(data,DBEAI>0&CAI>0)
      if(nrow(datatmp1)>0){
        datatmp1$AImod <- 0
      }
      if(nrow(datatmp2)>0){
        datatmp2$AImod <- 0
      }
      data <- rbind(datatmp1,datatmp2)%>%rbind(datatmp3)
      
    }else if(any(colnames(data) %in% "Br")){
      data$AImod <- (1+data$C-0.5*data$O-data$S-0.5*(data$N+data$P+data$H+data$Cl+data$Br+data$`81Br`))/(data$C-0.5*data$O-data$N-data$S-data$P)
      data$DBEAI <-(1+data$C-0.5*data$O-data$S-0.5*(data$N+data$P+data$H+data$Cl+data$Br+data$`81Br`))
      data$CAI <-(data$C-0.5*data$O-data$N-data$S-data$P)
      datatmp1 <- filter(data,DBEAI<=0)
      datatmp2 <- filter(data,CAI<=0)
      datatmp3 <- filter(data,DBEAI>0&CAI>0)
      if(nrow(datatmp1)>0){
        datatmp1$AImod <- 0
      }
      if(nrow(datatmp2)>0){
        datatmp2$AImod <- 0
      }
      data <- rbind(datatmp1,datatmp2)%>%rbind(datatmp3)
    }
  }
  
  if (!str_detect(str_flatten(colnames(data)),"DeltaG")){
    data$DeltaG <- 60.3-(data$NOSC)*28.5
  }
  if (!any(colnames(data) %in% "(DBE-O)/C")){
    data$`(DBE-O)/C`<-(data$DBE-data$O)/data$C
  }
  
  if (!any(colnames(data) %in% "(DBE-0.5O)/C")){
    data$`(DBE-0.5O)/C`<-(data$DBE-0.5*data$O)/data$C
  }
  
  if (!str_detect(str_flatten(colnames(data)),"groupKVD")){
    data <- KVD(data)
  }
  if (!str_detect(str_flatten(colnames(data)),"NewChemgroup")){
    data <- KVDNewGroup(data)
  }
  if (!str_detect(str_flatten(colnames(data)),"MLB")){
    data <- MLB(data)
  }
  if(str_detect(str_flatten(colnames(data)),"Group")&any(data$Group %in% 'Addcut Formula')){
    data<-data%>%select(-Group)
    data<- ElementGroup(data)
  }
  
  if (!str_detect(str_flatten(colnames(data)),"Group")){
    data <- ElementGroup(data)
  }
  if (str_detect(str_flatten(colnames(data)),"Formula")){
    if(!(any(str_detect(data$Formula,"Cl"))|any(str_detect(data$Formula,"Br"))|any(str_detect(data$Formula,"I")))){
      data <- as.data.frame(data)
      data1 <- data['Formula']
      colnames(data1)[1]<-'MolForm'
      df <- main_run(data1['MolForm'])
      data$lambda <-round(df$lambda,digits = 4)
    }
  }
  if (!str_detect(str_flatten(colnames(data)),"Xc")){
    data <- Xc(data,m=0.5,n=0.5)
  }
  data$`Measured m/z`<- trunc(data$`Measured m/z`*1000000)/1000000
  if(any(duplicated(data$`Measured m/z`))){
    data <- filter(data,!duplicated(data$`Measured m/z`))
  }
  data<-setDT(data)
  setkey(data,`Measured m/z`)
  data<-setDF(data)
  data <- data[,c("Measured m/z","C","H","O","N","C13","S","P","Error_ppm","RA","O/C","H/C","S/N","Neutral m/z","15N","13C","2H","S/C","18O","33S","34S","37Cl","81Br","Formula","ReInt","Intensity","DBE","DBE-O","DBE/C","N/C","P/C","NOSC","DeltaG","(DBE-0.5O)/C","AImod","groupKVD","NewChemgroup","MLB","Group","lambda","Xc")]
  return(data)
}


dataPrecess_Other <- function(data,Deu_DT=FALSE){
  if(any(colnames(data) %in% "mz")|any(colnames(data) %in% "m/z")){
    if(any(colnames(data) %in% "mz")){
      colnames(data)[which(colnames(data)=="mz")]<-'Measured m/z'
    }else if(any(colnames(data) %in% "m/z")){
      colnames(data)[which(colnames(data)=="m/z")]<-'Measured m/z'
    }
    data$`Measured m/z` <- as.numeric(data$`Measured m/z`)
  }
  if(str_detect(str_flatten(colnames(data)),"%")){
    data$RA <- data$`RA(%)`
  }
  if(!str_detect(str_flatten(colnames(data)),"Neutral m/z")){
    data$`Neutral m/z` <- data$`Measured m/z`+1.0072765
  }
  if(!any(colnames(data)%in%"C")){
    datatemp <- FormuTrans(data['Formula'])
    data <- cbind(data,datatemp[,-1])
  }
  data <- data[!duplicated(data$Formula),]
  
  if(!str_detect(str_flatten(colnames(data)),"15N")){
    data$`15N` <- 0
  }
  if(!str_detect(str_flatten(colnames(data)),"13C")){
    if(any(colnames(data)%in%"C13")){
      data$`13C` <- data$`C13`
    }else{
      data$`13C` <- 0
    }
    
  }

  
  
  if(!any(colnames(data) %in% "2H")){
    data$`2H` <- 0
  }
  if(!any(colnames(data) %in% "P")){
    data$P <- 0
  }
  if (!str_detect(str_flatten(colnames(data)),"H/C")){
    data$`H/C` <- data$H/data$C
    data$`O/C` <- data$O/data$C
  }
  if(!any(colnames(data) %in% "S/C")){
    data <- mutate(data,`S/C`=S/C)
  }
  if(!str_detect(str_flatten(colnames(data)),"18O")){
    data$`18O` <- 0
  }
  if(!str_detect(str_flatten(colnames(data)),"33S")){
    data$`33S` <- 0
  }
  if(!str_detect(str_flatten(colnames(data)),"34S")){
    data$`34S` <- 0
  }
  if(!str_detect(str_flatten(colnames(data)),"37Cl")){
    data$`37Cl` <- 0
  }
  if(!str_detect(str_flatten(colnames(data)),"81Br")){
    data$`81Br` <- 0
  }
  if(Deu_DT=="deu_dt"){
    data$D <- data$`2H`
    data$`2H` <- 0
    data$`H/C` <- (data$H+data$D)/(data$C+data$`13C`)
    data$`O/C` <- (data$O+data$`18O`)/(data$C+data$`13C`)
  }
  
  if (str_detect(str_flatten(colnames(data)),"Formula")){
    if(any(data$Formula %in% "Addcut Formula")){
      if(str_detect(str_flatten(colnames(data)),"Cl")){
        data$Formula <- paste('C',data$C+data$`13C`,'H',data$H+data$`2H`,'O',data$O+data$`18O`,'N',data$N+data$`15N`,'P',data$P,'S',data$S+data$`33S`+data$`34S`,'Cl',data$Cl+data$`37Cl`,'Br',data$Br+data$`81Br`,'I',data$I,sep = '')%>%str_remove_all('N0')%>% str_remove_all('P0') %>% str_remove_all('S0')%>% str_remove_all('Cl0')%>% str_remove_all('Br0')%>% str_remove_all('I0')
      }else{
        data$Formula <- paste('C',data$C+data$`13C`,'H',data$H+data$`2H`,'O',data$O+data$`18O`,'N',data$N+data$`15N`,'P',data$P,'S',data$S+data$`33S`+data$`34S`,sep = '')%>%str_remove_all('N0')%>% str_remove_all('P0') %>% str_remove_all('S0')
      }
    }
  }
  if (!str_detect(str_flatten(colnames(data)),"Formula")){
    if(Deu_DT=="deu_dt"){
      data$Formula <- paste('C',data$C+data$`13C`,'H',data$H+data$`2H`,'D',data$D,'O',data$O+data$`18O`,'N',data$N+data$`15N`,'P',data$P,'S',data$S+data$`33S`+data$`34S`,'Cl',data$Cl+data$`37Cl`,'Br',data$Br+data$`81Br`,'I',data$I,sep = '')%>%str_remove_all('D0')%>%str_remove_all('N0')%>% str_remove_all('P0') %>% str_remove_all('S0')%>% str_remove_all('Cl0')%>% str_remove_all('Br0')%>% str_remove_all('I0')
    }else{
      if(any(colnames(data) %in% "Cl")){
        data$Formula <- paste('C',data$C+data$`13C`,'H',data$H+data$`2H`,'O',data$O+data$`18O`,'N',data$N+data$`15N`,'P',data$P,'S',data$S+data$`33S`+data$`34S`,'Cl',data$Cl+data$`37Cl`,'Br',data$Br+data$`81Br`,'I',data$I,sep = '')%>%str_remove_all('N0')%>% str_remove_all('P0') %>% str_remove_all('S0')%>% str_remove_all('Cl0')%>% str_remove_all('Br0')%>% str_remove_all('I0')
        #data <- data %>% filter(`13C`==0) %>% filter(`2H`==0) %>% filter(`18O`==0)%>% filter(`15N`==0)%>% filter(`33S`==0)%>% filter(`34S`==0)%>% filter(`37Cl`==0)%>% filter(`81Br`==0)
      }else{
        data$Formula <- paste('C',data$C+data$`13C`,'H',data$H+data$`2H`,'O',data$O+data$`18O`,'N',data$N+data$`15N`,'P',data$P,'S',data$S+data$`33S`+data$`34S`,sep = '')%>%str_remove_all('N0')%>% str_remove_all('P0') %>% str_remove_all('S0')
        #data <- data %>% filter(`13C`==0) %>% filter(`2H`==0) %>% filter(`18O`==0)%>% filter(`15N`==0)%>% filter(`33S`==0)%>% filter(`34S`==0)
      }
    }
  }
  if(any(colnames(data) %in% "Cl")){
    if(any(colnames(data) %in% "33S")){
      data <- data %>% filter(`13C`==0) %>% filter(`2H`==0) %>% filter(`18O`==0)%>% filter(`15N`==0)%>% filter(`33S`==0)%>% filter(`34S`==0)%>% filter(`37Cl`==0)
    }else{
      data <- data %>% filter(`13C`==0) %>% filter(`2H`==0) %>% filter(`18O`==0)%>% filter(`15N`==0)%>% filter(`34S`==0)%>% filter(`37Cl`==0)
    }
  }else{
    if(any(colnames(data) %in% "33S")){
      data <- data %>% filter(`13C`==0) %>% filter(`2H`==0) %>% filter(`18O`==0)%>% filter(`15N`==0)%>% filter(`33S`==0)%>% filter(`34S`==0)
    }else{
      data <- data %>% filter(`13C`==0) %>% filter(`2H`==0) %>% filter(`18O`==0)%>% filter(`15N`==0)%>% filter(`34S`==0)
    }
  }
  
  if (!str_detect(str_flatten(colnames(data)),"ReInt")){
    data$ReInt <- data$Intensity*100/sum(data$Intensity)
  }
  
  if(!any(colnames(data) %in% "RA")){
    data$RA <- data$Intensity*100/max(data$Intensity)
  }

  if (!str_detect(str_flatten(colnames(data)),"DBE")){
    if(any(colnames(data) %in% "Cl")){
      data$DBE <- (2*data$C+data$N+data$P-data$H-data$Cl-data$Br-data$I+2)/2
    }else{
      data$DBE <- (2*data$C+data$N+data$P-data$H+2)/2
    }
    data$`DBE-O` <- data$DBE-data$O
    data$`DBE/C` <- data$DBE/data$C
  }
  if (!str_detect(str_flatten(colnames(data)),"N/C")){
    data$`N/C` <- data$N/data$C
    data$`P/C` <- data$P/data$C
  }
  if (!str_detect(str_flatten(colnames(data)),"S/N")){
    data$`S/N` <- 6
  }
  
  if (!str_detect(str_flatten(colnames(data)),"NOSC")){
    if(!any(colnames(data) %in% "Br")){
      data$NOSC <- 4-(4*data$C+data$H-2*data$O-3*data$N-2*data$S+5*data$P)/data$C
    }else if(any(colnames(data) %in% "Br")){
      data$NOSC <- 4-(4*data$C+data$H-2*data$O-3*data$N-2*data$S+5*data$P-data$Cl-data$Br-data$`81Br`-data$I)/data$C
    }
  }
  
  if (!str_detect(str_flatten(colnames(data)),"AImod")){
    if(!any(colnames(data) %in% "Br")){
      data$AImod <- (1+data$C-0.5*data$O-data$S-0.5*(data$N+data$P+data$H))/(data$C-0.5*data$O-data$N-data$S-data$P)
      data$DBEAI <- (1+data$C-0.5*data$O-data$S-0.5*(data$N+data$P+data$H))
      data$CAI <-(data$C-0.5*data$O-data$N-data$S-data$P)
      datatmp1 <- filter(data,DBEAI<=0)
      datatmp2 <- filter(data,CAI<=0)
      datatmp3 <- filter(data,DBEAI>0&CAI>0)
      if(nrow(datatmp1)>0){
        datatmp1$AImod <- 0
      }
      if(nrow(datatmp2)>0){
        datatmp2$AImod <- 0
      }
      data <- rbind(datatmp1,datatmp2)%>%rbind(datatmp3)
      
    }else if(any(colnames(data) %in% "Br")){
      data$AImod <- (1+data$C-0.5*data$O-data$S-0.5*(data$N+data$P+data$H+data$Cl+data$Br+data$`81Br`+data$I))/(data$C-0.5*data$O-data$N-data$S-data$P)
      data$DBEAI <-(1+data$C-0.5*data$O-data$S-0.5*(data$N+data$P+data$H+data$Cl+data$Br+data$`81Br`+data$I))
      data$CAI <-(data$C-0.5*data$O-data$N-data$S-data$P)
      datatmp1 <- filter(data,DBEAI<=0)
      datatmp2 <- filter(data,CAI<=0)
      datatmp3 <- filter(data,DBEAI>0&CAI>0)
      if(nrow(datatmp1)>0){
        datatmp1$AImod <- 0
      }
      if(nrow(datatmp2)>0){
        datatmp2$AImod <- 0
      }
      data <- rbind(datatmp1,datatmp2)%>%rbind(datatmp3)
    }
  }
  
  if (!str_detect(str_flatten(colnames(data)),"DeltaG")){
    data$DeltaG <- 60.3-(data$NOSC)*28.5
  }
  if (!any(colnames(data) %in% "(DBE-O)/C")){
    data$`(DBE-O)/C`<-(data$DBE-data$O)/data$C
  }
  
  if (!any(colnames(data) %in% "(DBE-0.5O)/C")){
    data$`(DBE-0.5O)/C`<-(data$DBE-0.5*data$O)/data$C
  }
  
  if (!str_detect(str_flatten(colnames(data)),"groupKVD")){
    data <- KVD(data)
  }
  if (!str_detect(str_flatten(colnames(data)),"NewChemgroup")){
    data <- KVDNewGroup(data)
  }
  if (!str_detect(str_flatten(colnames(data)),"MLB")){
    data <- MLB(data)
  }
  if(str_detect(str_flatten(colnames(data)),"Group")&any(data$Group %in% 'Addcut Formula')){
    data<-data%>%select(-Group)
    data<- ElementGroup(data)
  }
  
  if (!str_detect(str_flatten(colnames(data)),"Group")){
    data <- ElementGroup(data)
  }
  if (str_detect(str_flatten(colnames(data)),"Formula")){
    if(!(any(str_detect(data$Formula,"Cl"))|any(str_detect(data$Formula,"Br"))|any(str_detect(data$Formula,"I")))){
      data <- as.data.frame(data)
      data1 <- data['Formula']
      colnames(data1)[1]<-'MolForm'
      df <- main_run(data1['MolForm'])
      data$lambda <-round(df$lambda,digits = 4)
    }
  }
  if (!str_detect(str_flatten(colnames(data)),"Xc")){
    data <- Xc(data,m=0.5,n=0.5)
  }
  data$`Measured m/z`<- trunc(data$`Measured m/z`*1000000)/1000000
  if(any(duplicated(data$`Measured m/z`))){
    data <- filter(data,!duplicated(data$`Measured m/z`))
  }
  data<-setDT(data)
  setkey(data,`Measured m/z`)
  data<-setDF(data)
  if(any(colnames(data) %in% "lambda")){
    data <- data[,c("Measured m/z","C","H","O","N","S","P","Cl","Br","81Br","I","RA","O/C","H/C","S/N","Neutral m/z","15N","13C","2H","S/C","18O","33S","34S","37Cl","Formula","ReInt","Intensity","DBE","DBE-O","DBE/C","N/C","P/C","NOSC","DeltaG","(DBE-0.5O)/C","AImod","groupKVD","NewChemgroup","MLB","Group","lambda","Xc")]
  }else{
    data <- data[,c("Measured m/z","C","H","O","N","S","P","Cl","Br","81Br","I","RA","O/C","H/C","S/N","Neutral m/z","15N","13C","2H","S/C","18O","33S","34S","37Cl","Formula","ReInt","Intensity","DBE","DBE-O","DBE/C","N/C","P/C","NOSC","DeltaG","(DBE-0.5O)/C","AImod","groupKVD","NewChemgroup","MLB","Group","Xc")]
  }
  
  return(data)
}



ElementGroup <- function(FtmsDf) {
  Group <- tibble(Group=NA)
  for (i in 1:nrow(FtmsDf)) {
    Grouplabel <- str_flatten(str_sub(FtmsDf$Formula[i],str_locate(FtmsDf$Formula[i],c('C','H','D','O','N','S','P','Cl','Br','I')))[!is.na(str_sub(FtmsDf$Formula[i],str_locate(FtmsDf$Formula[i],c('C','H','D','O','N','S','P','Cl','Br','I'))))])
    Group_tmp <- tibble(Group=Grouplabel)
    Group <- rbind(Group,Group_tmp)
  }
  Group <- Group[2:nrow(Group),]
  FtmsDf <- cbind(FtmsDf,Group)
  return(FtmsDf)
}

MLB <- function(userdata){
  labileMLB <- filter(userdata,`H/C`>=1.5)
  nonlabileMLB <- filter(userdata,`H/C`<1.5)
  labileMLB$MLB <- 'labileMLB'
  nonlabileMLB$MLB <- 'nonlabileMLB'
  MLBdata <- rbind(labileMLB,nonlabileMLB)
  return(MLBdata)
}

NOSC <- function(userdata){
  FormuTransDf <- FormuTrans(userdata['Formula'])
  FormuTransDf$NOSC <- -((4*FormuTransDf$C+FormuTransDf$H-3*FormuTransDf$N-2*FormuTransDf$O+5*FormuTransDf$P-2*FormuTransDf$S)/FormuTransDf$C)+4
  data <- cbind(userdata[,-which(colnames(userdata)=='Formula')],FormuTransDf)
  #data$DeltaG <- 60.3-(data$NOSC)*28.5
}

Xc <- function(userdata,m=1,n=1){
  #FormuTransDf <- FormuTrans(userdata['Formula'])
  userdata <- as.data.frame(userdata)
  datatemp <- userdata[,c('Formula','C')]
  datatemp$XcInter = userdata$DBE-(m*(userdata$O+userdata$`18O`)+n*(userdata$S+userdata$`34S`))
  datatemp <- rownames_to_column(datatemp,'rowname')
  datatemp1 = filter(datatemp,XcInter <= 0)
  if(!nrow(datatemp1)==0){
    datatemp1$Xc = 0
  }
  datatemp2 = filter(datatemp,XcInter > 0)
  if(!nrow(datatemp2)==0){
    datatemp2$Xc = round((3*datatemp2$XcInter-2)/datatemp2$XcInter,digits = 4)
  }
  datatempNew = rbind(datatemp1,datatemp2)
  userdata <- rownames_to_column(userdata,'rowname')
  Xcdata <- left_join(userdata,datatempNew[c('rowname','Xc')],by='rowname')
  Xcdata <-Xcdata[,-which(colnames(Xcdata)=='rowname')]
  return(Xcdata)
}

KVD <- function(userdata){
  lipid <- filter(userdata,`H/C`>1.5&`H/C`<=2&`O/C`>0&`O/C`<=0.3)
  protein<-filter(userdata,`H/C`>1.5&`H/C`<=2.2&`O/C`>0.3&`O/C`<=0.67)
  carbohydrates<-filter(userdata,`H/C`>1.5&`H/C`<=2.4&`O/C`>0.67&`O/C`<=1.2)
  unsaturated <-filter(userdata,`H/C`>0.7&`H/C`<=1.5&`O/C`>0&`O/C`<=0.1)
  lignins <-filter(userdata,`H/C`>0.7&`H/C`<=1.5&`O/C`>0.1&`O/C`<=0.67)
  aromatic <-filter(userdata,`H/C`>0.2&`H/C`<=0.7&`O/C`>0&`O/C`<=0.67)
  tannin<-filter(userdata,`H/C`>0.6&`H/C`<=1.5&`O/C`>0.67&`O/C`<=1)
  Res <- list(lipid,protein,carbohydrates,unsaturated,lignins,aromatic,tannin)
  names(Res)<-c('Lipid','Aliphatic/protein','Carbohydrates','Unsaturatedhydrocarbons','Lignin/CRAM','AromaticStructures','Tannin')
  nrowRes <- unlist(lapply(Res, nrow))
  if(any(nrowRes==0)==TRUE){
    ResNew <- Res[-which(nrowRes==0)]
  }else{
    ResNew <- Res
  }
  for (i in 1:length(ResNew)) {
    ResNew[[i]]$groupKVD <- names(ResNew)[i]
  }
  KVDdata <- rbindlist(ResNew)
  Others <- userdata[!userdata$Intensity%in%KVDdata$Intensity,]
  if(!nrow(Others)==0){
    Others$groupKVD <- "Others"
    KVDdata <- rbind(KVDdata,Others)
  }
  
  return(KVDdata)
}

KVDNewGroup <- function(userdata){
  SaturatedFattyAcidLike <- filter(userdata,`H/C`>=2)
  PeptideLike <- filter(userdata,`H/C`>=1.5 & `H/C` < 2 & (N>0|`15N`>0))
  UnsaturatedAliphaticLike <- filter(userdata,`H/C`>=1.5 & `H/C` < 2 & (N==0&`15N`==0))
  HighlyUnsaturatedLike <- filter(userdata,`H/C`< 1.5&AImod <= 0.5)
  PolyphenolLike <- filter(userdata, AImod > 0.5 & AImod <= 0.66)
  CondensedAromaticLike <- filter(userdata,AImod > 0.66)
  SaturatedFattyAcidLike$NewChemgroup <- "SaturatedFattyAcidLike"
  PeptideLike$NewChemgroup <- "PeptideLike"
  UnsaturatedAliphaticLike$NewChemgroup <- "UnsaturatedAliphaticLike"
  HighlyUnsaturatedLike$NewChemgroup <- "HighlyUnsaturatedLike"
  PolyphenolLike$NewChemgroup <- "PolyphenolLike"
  CondensedAromaticLike$NewChemgroup <- "CondensedAromaticLike"
  NewChemgroup <- rbind(SaturatedFattyAcidLike,PeptideLike) %>% rbind(UnsaturatedAliphaticLike) %>% rbind(HighlyUnsaturatedLike) %>% rbind(PolyphenolLike) %>% rbind(CondensedAromaticLike)
  if (any(colnames(userdata) %in% "Intensity")){
    Others <- userdata[!userdata$Intensity%in%NewChemgroup$Intensity,] 
    if (!nrow(Others)==0){
      Others$NewChemgroup <- "Others"
      KVDdata <- rbind(NewChemgroup,Others) 
    }else{
      KVDdata <- NewChemgroup
    }
  }else{
    KVDdata <-NewChemgroup
  }
  return(KVDdata)
}

FormuTrans <- function(formula_ls){
  message("Transforing Formula to Element:\n" )
  if (ncol(formula_ls)>1){
    colnames(formula_ls) <- c("Reaction","Formula")
    Formula_df <- tibble(Formula = NA, C = NA, H=NA, N=NA, O=NA, P=NA, S=NA, Cl=NA, Br = NA, I=NA, F=NA)
    for (i in 1:nrow(formula_ls)) {
      formula <- formula_ls$Formula[i]
      formula_split <- str_split(formula, pattern = "")
      split_num <- as.numeric(formula_split[[1]])
      
      if(any(str_detect(formula_split[[1]],"-"))){
        element_id <- which(is.na(split_num)&!str_detect(formula_split[[1]],"-"))
        #num_id <-which(!is.na(split_num))
        #neg_id <- which(str_detect(formula_split[[1]],"-"))
        element_id <- c(element_id,length(formula_split[[1]])+1)
        formula_split_up <- c()
        for (i in 1:(length(element_id)-1)) {
          tmp <- c()
          for (j in 1:(element_id[i+1]-element_id[i]-1)) {
            tmp<- paste(tmp,formula_split[[1]][element_id[i]+j],sep = "")
        }
          formula_split_up <-c(formula_split_up,formula_split[[1]][element_id[i]],tmp)
        }
        formula_split <-list(formula_split_up)
      }
      
      split_num <- as.numeric(formula_split[[1]])
      split_num1<-split_num
      split_num[is.na(split_num)] <- 0
      if (!(str_detect(formula,"n")|str_detect(formula,"o")|str_detect(formula,"a")|str_detect(formula,"\\(")|str_detect(formula,"\\.")|is.na(formula))){
        if (str_detect(formula,"C")){
          C_index <- which(formula_split[[1]]%in%"C")[1]
          if (is.na(split_num1[C_index+1])){
            C_num <- 1
          }
          else if (!is.na(split_num1[C_index+2])){
            C_num <- split_num[C_index+1]*10+split_num[C_index+2]
          } else {C_num <- split_num[C_index+1]}
        } else {C_num <- 0}
        if (str_detect(formula,"H")){
          H_index <- which(formula_split[[1]]%in%"H")
          if (is.na(split_num1[H_index+1])){
            H_num <- 1
          }
          else if (!is.na(split_num1[H_index+2])){
            H_num <- split_num[H_index+1]*10+split_num[H_index+2]
          } else {H_num <- split_num[H_index+1]}
        } else {H_num <- 0}
        if (str_detect(formula,"O")){
          O_index <- which(formula_split[[1]]%in%"O")
          if (is.na(split_num1[O_index+1])){
            O_num <- 1
          }
          else if (!is.na(split_num1[O_index+2])){
            O_num <- split_num[O_index+1]*10+split_num[O_index+2]
          } else {O_num <- split_num[O_index+1]}
        } else {O_num <- 0}
        if (str_detect(formula,"N")){
          N_index <- which(formula_split[[1]]%in%"N")
          if (is.na(split_num1[N_index+1])){
            N_num <- 1
          }
          else if (!is.na(split_num1[N_index+2])){
            N_num <- split_num[N_index+1]*10+split_num[N_index+2]
          } else {N_num <- split_num[N_index+1]}
        } else {N_num <- 0}
        if (str_detect(formula,"P")){
          P_index <- which(formula_split[[1]]%in%"P")
          if (is.na(split_num1[P_index+1])){
            P_num <- 1
          }
          else if (!is.na(split_num1[P_index+2])){
            P_num <- split_num[P_index+1]*10+split_num[P_index+2]
          } else {P_num <- split_num[P_index+1]}
        } else {P_num <- 0}
        if (str_detect(formula,"S")){
          S_index <- which(formula_split[[1]]%in%"S")
          if (is.na(split_num1[S_index+1])){
            S_num <- 1
          }
          else if (!is.na(split_num1[S_index+2])){
            S_num <- split_num[S_index+1]*10+split_num[S_index+2]
          } else {S_num <- split_num[S_index+1]}
        } else {S_num <- 0}
        if (str_detect(formula,"l")){
          Cl_index <- which(formula_split[[1]]%in%"l")
          if (is.na(split_num1[Cl_index+1])){
            Cl_num <- 1
          }
          else if (!is.na(split_num1[Cl_index+2])){
            Cl_num <- split_num[Cl_index+1]*10+split_num[Cl_index+2]
          } else {Cl_num <- split_num[Cl_index+1]}
        } else {Cl_num <- 0}
        if (str_detect(formula,"r")){
          Br_index <- which(formula_split[[1]]%in%"r")
          if (is.na(split_num1[Br_index+1])){
            Br_num <- 1
          }
          else if (!is.na(split_num1[Br_index+2])){
            Br_num <- split_num[Br_index+1]*10+split_num[Br_index+2]
          } else {Br_num <- split_num[Br_index+1]}
        } else {Br_num <- 0}
        if (str_detect(formula,"I")){
          I_index <- which(formula_split[[1]]%in%"I")
          if (is.na(split_num1[I_index+1])){
            I_num <- 1
          }
          else if (!is.na(split_num1[I_index+2])){
            I_num <- split_num[I_index+1]*10+split_num[I_index+2]
          } else {I_num <- split_num[I_index+1]}
        } else {I_num <- 0}
        if (str_detect(formula,"F")){
          F_index <- which(formula_split[[1]]%in%"F")
          if (is.na(split_num1[F_index+1])){
            F_num <- 1
          }
          else if (!is.na(split_num1[F_index+2])){
            F_num <- split_num[F_index+1]*10+split_num[F_index+2]
          } else {F_num <- split_num[F_index+1]}
        } else {F_num <- 0}
        Formula_temp <- tibble(Reaction = formula_ls$Reaction[i], Formula = formula_ls$Formula[i], C = C_num, H=H_num, N=N_num, O=O_num, P=P_num, S=S_num,Cl=Cl_num, Br=Br_num, I=I_num, F=F_num)
      } else {
        Formula_temp <- tibble(Reaction = formula_ls$Reaction[i], Formula = formula_ls$Formula[i], C = NA, H=NA, N=NA, O=NA, P=NA, S=NA, Cl=NA, Br=NA, I=NA, F=NA)
      }
      Formula_df <- rbind(Formula_df,Formula_temp)
    }
    Formula_df <- Formula_df[2:nrow(Formula_df),]
  } else {
    colnames(formula_ls) <- "Formula"
    Formula_df <- tibble(Formula = NA, C = NA, H=NA, N=NA, O=NA, P=NA, S=NA, Cl=NA, Br = NA, I=NA, F=NA)
    for (i in 1:nrow(formula_ls)) {
      formula <- formula_ls$Formula[i]
      formula_split <- str_split(formula, pattern = "")
      split_num <- as.numeric(formula_split[[1]])
      
      if(any(str_detect(formula_split[[1]],"-"))){
        element_id <- which(is.na(split_num)&!str_detect(formula_split[[1]],"-"))
        #num_id <-which(!is.na(split_num))
        #neg_id <- which(str_detect(formula_split[[1]],"-"))
        element_id <- c(element_id,length(formula_split[[1]])+1)
        formula_split_up <- c()
        for (i in 1:(length(element_id)-1)) {
          tmp <- c()
          for (j in 1:(element_id[i+1]-element_id[i]-1)) {
            tmp<- paste(tmp,formula_split[[1]][element_id[i]+j],sep = "")
          }
          formula_split_up <-c(formula_split_up,formula_split[[1]][element_id[i]],tmp)
        }
        formula_split <-list(formula_split_up)
      }
      
      split_num <- as.numeric(formula_split[[1]])
      split_num1<-split_num
      split_num[is.na(split_num)] <- 0
      if (!(str_detect(formula,"n")|str_detect(formula,"o")|str_detect(formula,"a")|str_detect(formula,"\\(")|str_detect(formula,"\\.")|is.na(formula))){
        if (str_detect(formula,"C")){
          C_index <- which(formula_split[[1]]%in%"C")[1]
          if (is.na(split_num1[C_index+1])){
            C_num <- 1
          }
          else if (!is.na(split_num1[C_index+2])){
            C_num <- split_num[C_index+1]*10+split_num[C_index+2]
          } else {C_num <- split_num[C_index+1]}
        } else {C_num <- 0}
        if (str_detect(formula,"H")){
          H_index <- which(formula_split[[1]]%in%"H")
          if (is.na(split_num1[H_index+1])){
            H_num <- 1
          }
          else if (!is.na(split_num1[H_index+2])){
            H_num <- split_num[H_index+1]*10+split_num[H_index+2]
          } else {H_num <- split_num[H_index+1]}
        } else {H_num <- 0}
        if (str_detect(formula,"O")){
          O_index <- which(formula_split[[1]]%in%"O")
          if (is.na(split_num1[O_index+1])){
            O_num <- 1
          }
          else if (!is.na(split_num1[O_index+2])){
            O_num <- split_num[O_index+1]*10+split_num[O_index+2]
          } else {O_num <- split_num[O_index+1]}
        } else {O_num <- 0}
        if (str_detect(formula,"N")){
          N_index <- which(formula_split[[1]]%in%"N")
          if (is.na(split_num1[N_index+1])){
            N_num <- 1
          }
          else if (!is.na(split_num1[N_index+2])){
            N_num <- split_num[N_index+1]*10+split_num[N_index+2]
          } else {N_num <- split_num[N_index+1]}
        } else {N_num <- 0}
        if (str_detect(formula,"P")){
          P_index <- which(formula_split[[1]]%in%"P")
          if (is.na(split_num1[P_index+1])){
            P_num <- 1
          }
          else if (!is.na(split_num1[P_index+2])){
            P_num <- split_num[P_index+1]*10+split_num[P_index+2]
          } else {P_num <- split_num[P_index+1]}
        } else {P_num <- 0}
        if (str_detect(formula,"S")){
          S_index <- which(formula_split[[1]]%in%"S")
          if (is.na(split_num1[S_index+1])){
            S_num <- 1
          }
          else if (!is.na(split_num1[S_index+2])){
            S_num <- split_num[S_index+1]*10+split_num[S_index+2]
          } else {S_num <- split_num[S_index+1]}
        } else {S_num <- 0}
        if (str_detect(formula,"l")){
          Cl_index <- which(formula_split[[1]]%in%"l")
          if (is.na(split_num1[Cl_index+1])){
            Cl_num <- 1
          }
          else if (!is.na(split_num1[Cl_index+2])){
            Cl_num <- split_num[Cl_index+1]*10+split_num[Cl_index+2]
          } else {Cl_num <- split_num[Cl_index+1]}
        } else {Cl_num <- 0}
        if (str_detect(formula,"r")){
          Br_index <- which(formula_split[[1]]%in%"r")
          if (is.na(split_num1[Br_index+1])){
            Br_num <- 1
          }
          else if (!is.na(split_num1[Br_index+2])){
            Br_num <- split_num[Br_index+1]*10+split_num[Br_index+2]
          } else {Br_num <- split_num[Br_index+1]}
        } else {Br_num <- 0}
        if (str_detect(formula,"I")){
          I_index <- which(formula_split[[1]]%in%"I")
          if (is.na(split_num1[I_index+1])){
            I_num <- 1
          }
          else if (!is.na(split_num1[I_index+2])){
            I_num <- split_num[I_index+1]*10+split_num[I_index+2]
          } else {I_num <- split_num[I_index+1]}
        } else {I_num <- 0}
        if (str_detect(formula,"F")){
          F_index <- which(formula_split[[1]]%in%"F")
          if (is.na(split_num1[F_index+1])){
            F_num <- 1
          }
          else if (!is.na(split_num1[F_index+2])){
            F_num <- split_num[F_index+1]*10+split_num[F_index+2]
          } else {F_num <- split_num[F_index+1]}
        } else {F_num <- 0}
        Formula_temp <- tibble(Formula = formula, C = C_num, H=H_num, N=N_num, O=O_num, P=P_num, S=S_num,Cl=Cl_num, Br=Br_num, I=I_num, F=F_num)
      } else {
        Formula_temp <- tibble(Reaction = formula_ls$Reaction[i], Formula = formula_ls$Formula[i], C = NA, H=NA, N=NA, O=NA, P=NA, S=NA, Cl=NA, Br=NA, I=NA, F=NA)
      }
      Formula_df <- rbind(Formula_df,Formula_temp)
    }
    Formula_df <- Formula_df[2:nrow(Formula_df),]
  }
  return(Formula_df)
}

VKDGrah1 <- function(data,pointsize=TRUE){
    if (pointsize){
      P1Element1 <- ggscatter(data,'O/C','H/C',fill = 'Group',shape = 21,color = 'black',size = 'ReInt',alpha = 0.6,stroke=0.2)+
        fill_palette('jco')+scale_size(range = c(0.8,3))+
        theme_prism(base_size = 7,border = TRUE,base_fontface = "plain",base_line_size = 8/50,base_rect_size = 8/50)+facet_wrap(vars(Group),ncol = 4)+
        theme(strip.background = element_rect(fill = "grey"),strip.placement = "outside",legend.key.height=unit(0.8,"line"),legend.margin = ggplot2::margin(0,0,0,0),legend.box.margin=ggplot2::margin(0,-10,-5,-10),plot.margin = ggplot2::margin(2,2,2,2,'mm'))+
        guides(fill=guide_legend(override.aes = list(size=2)))
      
    }else{
      P1Element1 <- ggscatter(data,'O/C','H/C',fill = 'Group',shape = 21,color = 'black',size = 0.8,alpha = 0.6,stroke=0.2)+
        fill_palette('jco')+
        theme_prism(base_size = 7,border = TRUE,base_fontface = "plain",base_line_size = 8/50,base_rect_size = 8/50)+facet_wrap(vars(Group),ncol = 4)+
        theme(strip.background = element_rect(fill = "grey"),strip.placement = "outside",legend.key.height=unit(0.8,"line"),legend.margin = ggplot2::margin(0,0,0,0),legend.box.margin=ggplot2::margin(0,-10,-5,-10),plot.margin = ggplot2::margin(2,2,2,2,'mm'))+
        guides(fill=guide_legend(override.aes = list(size=2)))
    }
    P1Element1 <- ggpar(P1Element1,xlim = c(0,1.2),ylim = c(0,2.5),legend = 'top')
    P1Element1Up <- P1Element1+geom_hline(aes(yintercept=1.5),linetype="longdash",color="darkgrey",linewidth=0.1)+
      geom_linerange(aes(xmin=0,xmax=0.3,y=2.0),linetype="longdash",color="darkgrey",linewidth=0.1)+
      geom_linerange(aes(x=0,ymin=1.5,ymax=2.0),linetype="longdash",color="darkgrey",linewidth=0.1)+
      geom_linerange(aes(x=0.3,ymin=1.5,ymax=2.0),linetype="longdash",color="darkgrey",linewidth=0.1)+
      geom_linerange(aes(x=0.3,ymin=1.5,ymax=2.2),linetype="longdash",color="darkgrey",linewidth=0.1)+
      geom_linerange(aes(x=0.67,ymin=1.5,ymax=2.2),linetype="longdash",color="darkgrey",linewidth=0.1)+
      geom_linerange(aes(xmin=0.3,xmax=0.67,y=2.2),linetype="longdash",color="darkgrey",linewidth=0.1)+
      geom_linerange(aes(x=0.67,ymin=1.5,ymax=2.4),linetype="longdash",color="darkgrey",linewidth=0.1)+
      geom_linerange(aes(x=1.2,ymin=1.5,ymax=2.4),linetype="longdash",color="darkgrey",linewidth=0.1)+
      geom_linerange(aes(xmin=0.67,xmax=1.2,y=2.4),linetype="longdash",color="darkgrey",linewidth=0.1)+
      geom_linerange(aes(x=0,ymin=0.7,ymax=1.5),linetype="longdash",color="darkgrey",linewidth=0.1)+
      geom_linerange(aes(x=0.1,ymin=0.7,ymax=1.5),linetype="longdash",color="darkgrey",linewidth=0.1)+
      geom_linerange(aes(xmin=0,xmax=0.1,y=0.7),linetype="longdash",color="darkgrey",linewidth=0.1)+
      geom_linerange(aes(x=0.1,ymin=0.7,ymax=1.5),linetype="longdash",color="darkgrey",linewidth=0.1)+
      geom_linerange(aes(xmin=0.1,xmax=0.67,y=0.7),linetype="longdash",color="darkgrey",linewidth=0.1)+
      geom_linerange(aes(x=0,ymin=0.2,ymax=0.7),linetype="longdash",color="darkgrey",linewidth=0.1)+
      geom_linerange(aes(x=0.67,ymin=0.2,ymax=0.7),linetype="longdash",color="darkgrey",linewidth=0.1)+
      geom_linerange(aes(xmin=0,xmax=0.67,y=0.2),linetype="longdash",color="darkgrey",linewidth=0.1)+
      geom_linerange(aes(x=0.67,ymin=0.6,ymax=1.5),linetype="longdash",color="darkgrey",linewidth=0.1)+
      geom_linerange(aes(x=1,ymin=0.6,ymax=1.5),linetype="longdash",color="darkgrey",linewidth=0.1)+
      geom_linerange(aes(xmin=0.67,xmax=1,y=0.6),linetype="longdash",color="darkgrey",linewidth=0.1)
   
  return(P1Element1Up)
}

VKDGrah2 <- function(data,pointsize=TRUE){
  if (pointsize){
    P1Element1 <- ggscatter(data,'O/C','H/C',fill = 'Group',shape = 21,color = 'black',size = 'ReInt',alpha = 0.6,stroke=0.2)+
      fill_palette('jco')+scale_size(range = c(0.8,3))+
      theme_prism(base_size = 7,border = TRUE,base_fontface = "plain",base_line_size = 8/50,base_rect_size = 8/50)+facet_wrap(vars(Group),ncol = 4)+
      theme(strip.background = element_rect(fill = "grey"),strip.placement = "outside",legend.key.height=unit(0.8,"line"),legend.margin = ggplot2::margin(0,0,0,0),legend.box.margin=ggplot2::margin(0,-10,-5,-10),plot.margin = ggplot2::margin(2,2,2,2,'mm'))+
      guides(fill=guide_legend(override.aes = list(size=2)))
    
  }else{
    P1Element1 <- ggscatter(data,'O/C','H/C',fill = 'Group',shape = 21,color = 'black',size = 0.8,alpha = 0.6,stroke=0.2)+
      fill_palette('jco')+
      theme_prism(base_size = 7,border = TRUE,base_fontface = "plain",base_line_size = 8/50,base_rect_size = 8/50)+facet_wrap(vars(Group),ncol = 4)+
      theme(strip.background = element_rect(fill = "grey"),strip.placement = "outside",legend.key.height=unit(0.8,"line"),legend.margin = ggplot2::margin(0,0,0,0),legend.box.margin=ggplot2::margin(0,-10,-5,-10),plot.margin = ggplot2::margin(2,2,2,2,'mm'))+
      guides(fill=guide_legend(override.aes = list(size=2)))
  }
  P1Element1 <- ggpar(P1Element1,xlim = c(0,1.2),ylim = c(0,2.5),legend = 'top')
  P1Element1Up <- P1Element1+geom_hline(aes(yintercept=2),linetype="longdash",color="darkgrey",linewidth=0.1)+
    geom_hline(aes(yintercept=1.5),linetype="longdash",color="darkgrey",linewidth=0.1)+
    geom_abline(intercept = 1.05,slope = -0.38,linetype="longdash",color="darkgrey",linewidth=0.1)+
    geom_abline(intercept = 0.8,slope = -0.38,linetype="longdash",color="darkgrey",linewidth=0.1)
  
  return(P1Element1Up)
}

VKDGrah3 <- function(data,pointsize=TRUE){
  if (pointsize){
    P1Element1_size_legend <- ggscatter(data,'O/C','H/C',shape = 21,color = 'black',fill = 'groupKVD',size = 'ReInt',alpha = 0.6,stroke=0.2)+
      scale_fill_manual(values = paletteer_dynamic("cartography::multi.pal", 20))+scale_size(range = c(0.8,3))+xlim(0,1.2)+ylim(0,2.5)+
      theme_prism(base_size = 7,border = TRUE,base_fontface = "plain",base_line_size = 8/50,base_rect_size = 8/50)+
      theme(legend.position = c(0.7, 0.95),legend.direction = 'horizontal')+guides(fill='none') 
    P1Element1_size_legend <- P1Element1_size_legend+geom_hline(aes(yintercept=1.5),linetype="longdash",color="darkgrey",linewidth=0.1)+
      geom_linerange(aes(xmin=0,xmax=0.3,y=2.0),linetype="longdash",color="darkgrey",linewidth=0.1)+
      geom_linerange(aes(x=0,ymin=1.5,ymax=2.0),linetype="longdash",color="darkgrey",linewidth=0.1)+
      geom_linerange(aes(x=0.3,ymin=1.5,ymax=2.0),linetype="longdash",color="darkgrey",linewidth=0.1)+
      geom_linerange(aes(x=0.3,ymin=1.5,ymax=2.2),linetype="longdash",color="darkgrey",linewidth=0.1)+
      geom_linerange(aes(x=0.67,ymin=1.5,ymax=2.2),linetype="longdash",color="darkgrey",linewidth=0.1)+
      geom_linerange(aes(xmin=0.3,xmax=0.67,y=2.2),linetype="longdash",color="darkgrey",linewidth=0.1)+
      geom_linerange(aes(x=0.67,ymin=1.5,ymax=2.4),linetype="longdash",color="darkgrey",linewidth=0.1)+
      geom_linerange(aes(x=1.2,ymin=1.5,ymax=2.4),linetype="longdash",color="darkgrey",linewidth=0.1)+
      geom_linerange(aes(xmin=0.67,xmax=1.2,y=2.4),linetype="longdash",color="darkgrey",linewidth=0.1)+
      geom_linerange(aes(x=0,ymin=0.7,ymax=1.5),linetype="longdash",color="darkgrey",linewidth=0.1)+
      geom_linerange(aes(x=0.1,ymin=0.7,ymax=1.5),linetype="longdash",color="darkgrey",linewidth=0.1)+
      geom_linerange(aes(xmin=0,xmax=0.1,y=0.7),linetype="longdash",color="darkgrey",linewidth=0.1)+
      geom_linerange(aes(x=0.1,ymin=0.7,ymax=1.5),linetype="longdash",color="darkgrey",linewidth=0.1)+
      geom_linerange(aes(xmin=0.1,xmax=0.67,y=0.7),linetype="longdash",color="darkgrey",linewidth=0.1)+
      geom_linerange(aes(x=0,ymin=0.2,ymax=0.7),linetype="longdash",color="darkgrey",linewidth=0.1)+
      geom_linerange(aes(x=0.67,ymin=0.2,ymax=0.7),linetype="longdash",color="darkgrey",linewidth=0.1)+
      geom_linerange(aes(xmin=0,xmax=0.67,y=0.2),linetype="longdash",color="darkgrey",linewidth=0.1)+
      geom_linerange(aes(x=0.67,ymin=0.6,ymax=1.5),linetype="longdash",color="darkgrey",linewidth=0.1)+
      geom_linerange(aes(x=1,ymin=0.6,ymax=1.5),linetype="longdash",color="darkgrey",linewidth=0.1)+
      geom_linerange(aes(xmin=0.67,xmax=1,y=0.6),linetype="longdash",color="darkgrey",linewidth=0.1)
    
    P1Element1_fill_legend <- ggscatter(data,'O/C','H/C',shape = 21,color = 'black',fill = 'groupKVD',alpha = 0.6,stroke=0.2)+
      scale_fill_manual(values = paletteer_dynamic("cartography::multi.pal", 20))+scale_size(range = c(0.8,3))+
      theme_prism(base_size = 7,border = TRUE,base_fontface = "plain",base_line_size = 8/50,base_rect_size = 8/50)+
      theme(legend.position = 'top',legend.key.height=unit(0.2,"line"),legend.margin = ggplot2::margin(2,0,0,0),legend.box.margin=ggplot2::margin(2,-10,0,-10),plot.margin = unit(c(0.2,1,0.6,0.6),'cm'))+
      guides(fill=guide_legend(override.aes = list(size=2)))
    
    guide_fill <- get_legend(P1Element1_fill_legend)
    
    P1Element1Up <- plot_grid(guide_fill,P1Element1_size_legend, 
              nrow = 2, rel_heights = c(1, 10))
  }else{
    P1Element1 <- ggscatter(data,'O/C','H/C',fill = 'groupKVD',shape = 21,color = 'black',size = 0.8,alpha = 0.6,stroke=0.2)+
      scale_fill_manual(values = paletteer_dynamic("cartography::multi.pal", 20))+xlim(0,1.2)+ylim(0,2.5)+
      theme_prism(base_size = 7,border = TRUE,base_fontface = "plain",base_line_size = 8/50,base_rect_size = 8/50)+
      theme(legend.position = 'top',legend.key.height=unit(0.2,"line"),legend.margin = ggplot2::margin(2,0,0,0),legend.box.margin=ggplot2::margin(2,-10,0,-10),plot.margin = unit(c(0.2,1,0.6,0.6),'cm'))+
      guides(fill=guide_legend(override.aes = list(size=2)))
    P1Element1Up <- P1Element1+geom_hline(aes(yintercept=1.5),linetype="longdash",color="darkgrey",linewidth=0.1)+
      geom_linerange(aes(xmin=0,xmax=0.3,y=2.0),linetype="longdash",color="darkgrey",linewidth=0.1)+
      geom_linerange(aes(x=0,ymin=1.5,ymax=2.0),linetype="longdash",color="darkgrey",linewidth=0.1)+
      geom_linerange(aes(x=0.3,ymin=1.5,ymax=2.0),linetype="longdash",color="darkgrey",linewidth=0.1)+
      geom_linerange(aes(x=0.3,ymin=1.5,ymax=2.2),linetype="longdash",color="darkgrey",linewidth=0.1)+
      geom_linerange(aes(x=0.67,ymin=1.5,ymax=2.2),linetype="longdash",color="darkgrey",linewidth=0.1)+
      geom_linerange(aes(xmin=0.3,xmax=0.67,y=2.2),linetype="longdash",color="darkgrey",linewidth=0.1)+
      geom_linerange(aes(x=0.67,ymin=1.5,ymax=2.4),linetype="longdash",color="darkgrey",linewidth=0.1)+
      geom_linerange(aes(x=1.2,ymin=1.5,ymax=2.4),linetype="longdash",color="darkgrey",linewidth=0.1)+
      geom_linerange(aes(xmin=0.67,xmax=1.2,y=2.4),linetype="longdash",color="darkgrey",linewidth=0.1)+
      geom_linerange(aes(x=0,ymin=0.7,ymax=1.5),linetype="longdash",color="darkgrey",linewidth=0.1)+
      geom_linerange(aes(x=0.1,ymin=0.7,ymax=1.5),linetype="longdash",color="darkgrey",linewidth=0.1)+
      geom_linerange(aes(xmin=0,xmax=0.1,y=0.7),linetype="longdash",color="darkgrey",linewidth=0.1)+
      geom_linerange(aes(x=0.1,ymin=0.7,ymax=1.5),linetype="longdash",color="darkgrey",linewidth=0.1)+
      geom_linerange(aes(xmin=0.1,xmax=0.67,y=0.7),linetype="longdash",color="darkgrey",linewidth=0.1)+
      geom_linerange(aes(x=0,ymin=0.2,ymax=0.7),linetype="longdash",color="darkgrey",linewidth=0.1)+
      geom_linerange(aes(x=0.67,ymin=0.2,ymax=0.7),linetype="longdash",color="darkgrey",linewidth=0.1)+
      geom_linerange(aes(xmin=0,xmax=0.67,y=0.2),linetype="longdash",color="darkgrey",linewidth=0.1)+
      geom_linerange(aes(x=0.67,ymin=0.6,ymax=1.5),linetype="longdash",color="darkgrey",linewidth=0.1)+
      geom_linerange(aes(x=1,ymin=0.6,ymax=1.5),linetype="longdash",color="darkgrey",linewidth=0.1)+
      geom_linerange(aes(xmin=0.67,xmax=1,y=0.6),linetype="longdash",color="darkgrey",linewidth=0.1)
  }
  return(P1Element1Up)
}

VKDGrah4 <- function(data,pointsize=TRUE){
  if (pointsize){
    P1Element1_size_legend <- ggscatter(data,'O/C','H/C',shape = 21,color = 'black',fill = 'NewChemgroup',size = 'ReInt',alpha = 0.6,stroke=0.2)+
      scale_fill_manual(values = paletteer_dynamic("cartography::multi.pal", 20))+scale_size(range = c(0.8,3))+xlim(0,1.2)+ylim(0,2.5)+
      theme_prism(base_size = 7,border = TRUE,base_fontface = "plain",base_line_size = 8/50,base_rect_size = 8/50)+theme(legend.position = c(0.7, 0.95),legend.direction = 'horizontal')+guides(fill='none')
    P1Element1_size_legend <- P1Element1_size_legend+geom_hline(aes(yintercept=2),linetype="longdash",color="darkgrey",linewidth=0.1)+
      geom_hline(aes(yintercept=1.5),linetype="longdash",color="darkgrey",linewidth=0.1)+
      geom_abline(intercept = 1.05,slope = -0.38,linetype="longdash",color="darkgrey",linewidth=0.1)+
      geom_abline(intercept = 0.8,slope = -0.38,linetype="longdash",color="darkgrey",linewidth=0.1)
    P1Element1_fill_legend <- ggscatter(data,'O/C','H/C',shape = 21,color = 'black',fill = 'NewChemgroup',alpha = 0.6,stroke=0.2)+
      scale_fill_manual(values = paletteer_dynamic("cartography::multi.pal", 20))+scale_size(range = c(0.8,3))+
      theme_prism(base_size = 7,border = TRUE,base_fontface = "plain",base_line_size = 8/50,base_rect_size = 8/50)+
      theme(legend.position = 'top',legend.key.height=unit(0.2,"line"),legend.margin = ggplot2::margin(2,0,0,0),legend.box.margin=ggplot2::margin(2,-10,0,-10),plot.margin = unit(c(0.2,1,0.6,0.6),'cm'))+
      guides(fill=guide_legend(override.aes = list(size=2)))
    guide_fill <- get_legend(P1Element1_fill_legend)
    P1Element1Up <- plot_grid(guide_fill,P1Element1_size_legend, 
                              nrow = 2, rel_heights = c(1, 10))
  }else{
    P1Element1 <- ggscatter(data,'O/C','H/C',fill = 'NewChemgroup',shape = 21,color = 'black',size = 0.8,alpha = 0.6,stroke=0.2)+
      scale_fill_manual(values = paletteer_dynamic("cartography::multi.pal", 20))+
      theme_prism(base_size = 7,border = TRUE,base_fontface = "plain",base_line_size = 8/50,base_rect_size = 8/50)+
      theme(legend.position = 'top',legend.key.height=unit(0.2,"line"),legend.margin = ggplot2::margin(2,0,0,0),legend.box.margin=ggplot2::margin(2,-10,0,-10),plot.margin = unit(c(0.2,1,0.6,0.6),'cm'))+
      guides(fill=guide_legend(override.aes = list(size=2)))
    P1Element1Up <- P1Element1+geom_hline(aes(yintercept=2),linetype="longdash",color="darkgrey",linewidth=0.1)+
      geom_hline(aes(yintercept=1.5),linetype="longdash",color="darkgrey",linewidth=0.1)+
      geom_abline(intercept = 1.05,slope = -0.38,linetype="longdash",color="darkgrey",linewidth=0.1)+
      geom_abline(intercept = 0.8,slope = -0.38,linetype="longdash",color="darkgrey",linewidth=0.1)
  }
  return(P1Element1Up)
}

HighlightVKPlot <- function(FTMSDTLS,biochemcal_type="NewChemgroup",pointsize=FALSE,highlight_filter=NULL,deu_dt=FALSE){
  if(is.null(highlight_filter)){
    highlight_filter=names(FTMSDTLS)[1]
  }
  
  for (i in 1:length(FTMSDTLS)) {
    FTMSDTLS[[i]]$Sample <- names(FTMSDTLS)[i]
    if(deu_dt=="deu_dt"){
      FTMSDTLS[[i]]$H <- FTMSDTLS[[i]]$H+FTMSDTLS[[i]]$D
    }
    
  }
  datatmp0 <- bind_rows(FTMSDTLS)
  #datatmp1 <- datatmp0
  #datatmp1$Sample <- "All"
  #data <- rbind(datatmp1,datatmp0)
  if(biochemcal_type=="NewChemgroup"){
    if (pointsize){
      P1Element1_size_legend <- ggscatter(datatmp0,'O/C','H/C',shape = 21,color = 'black',fill = 'Sample',size = 'ReInt',alpha = 0.8,stroke=0.2)+gghighlight(Sample %in% highlight_filter)+
        scale_fill_manual(values = paletteer_dynamic("cartography::multi.pal", 20))+scale_size(range = c(0.8,3))+xlim(0,1.2)+ylim(0,2.5)+
        theme_prism(base_size = 7,border = TRUE,base_fontface = "plain",base_line_size = 8/50,base_rect_size = 8/50)+theme(legend.position = c(0.7, 0.95),legend.direction = 'horizontal')+guides(fill='none')
      P1Element1_size_legend <- P1Element1_size_legend+geom_hline(aes(yintercept=2),linetype="longdash",color="darkgrey",linewidth=0.1)+
        geom_hline(aes(yintercept=1.5),linetype="longdash",color="darkgrey",linewidth=0.1)+
        geom_abline(intercept = 1.05,slope = -0.38,linetype="longdash",color="darkgrey",linewidth=0.1)+
        geom_abline(intercept = 0.8,slope = -0.38,linetype="longdash",color="darkgrey",linewidth=0.1)
      P1Element1_fill_legend <- ggscatter(datatmp0,'O/C','H/C',shape = 21,color = 'black',fill = 'Sample',alpha = 0.6,stroke=0.2)+
        scale_fill_manual(values = paletteer_dynamic("cartography::multi.pal", 20))+scale_size(range = c(0.8,3))+
        theme_prism(base_size = 7,border = TRUE,base_fontface = "plain",base_line_size = 8/50,base_rect_size = 8/50)+
        theme(legend.position = 'top',legend.key.height=unit(0.2,"line"),legend.margin = ggplot2::margin(2,0,0,0),legend.box.margin=ggplot2::margin(2,-10,0,-10),plot.margin = unit(c(0.2,1,0.6,0.6),'cm'))+
        guides(fill=guide_legend(override.aes = list(size=2)))
      guide_fill <- get_legend(P1Element1_fill_legend)
      P1Element1Up <- plot_grid(guide_fill,P1Element1_size_legend, 
                                nrow = 2, rel_heights = c(1, 10))
    }else{
      P1Element1 <- ggscatter(datatmp0,'O/C','H/C',fill = 'Sample', shape = 21,color = 'black',size = 1,alpha = 0.8,stroke=0.2)+gghighlight(Sample %in% highlight_filter)+
        scale_fill_manual(values = paletteer_dynamic("cartography::multi.pal", 20))+
        theme_prism(base_size = 7,border = TRUE,base_fontface = "plain",base_line_size = 8/50,base_rect_size = 8/50)+
        theme(legend.position = 'top',legend.key.height=unit(0.2,"line"),legend.margin = ggplot2::margin(2,0,0,0),legend.box.margin=ggplot2::margin(2,-10,0,-10),plot.margin = unit(c(0.2,1,0.6,0.6),'cm'))+
        guides(fill=guide_legend(override.aes = list(size=2)))
      P1Element1Up <- P1Element1+geom_hline(aes(yintercept=2),linetype="longdash",color="darkgrey",linewidth=0.1)+
        geom_hline(aes(yintercept=1.5),linetype="longdash",color="darkgrey",linewidth=0.1)+
        geom_abline(intercept = 1.05,slope = -0.38,linetype="longdash",color="darkgrey",linewidth=0.1)+
        geom_abline(intercept = 0.8,slope = -0.38,linetype="longdash",color="darkgrey",linewidth=0.1)
    }
  }else{
    if (pointsize){
      P1Element1_size_legend <- ggscatter(datatmp0,'O/C','H/C',shape = 21,color = 'black',fill = 'Sample',size = 'ReInt',alpha = 0.8,stroke=0.2)+gghighlight(Sample %in% highlight_filter)+
        scale_fill_manual(values = paletteer_dynamic("cartography::multi.pal", 20))+scale_size(range = c(0.8,3))+xlim(0,1.2)+ylim(0,2.5)+
        theme_prism(base_size = 7,border = TRUE,base_fontface = "plain",base_line_size = 8/50,base_rect_size = 8/50)+theme(legend.position = c(0.7, 0.95),legend.direction = 'horizontal')+guides(fill='none')
      P1Element1_size_legend <- P1Element1_size_legend+geom_hline(aes(yintercept=1.5),linetype="longdash",color="darkgrey",linewidth=0.1)+
        geom_linerange(aes(xmin=0,xmax=0.3,y=2.0),linetype="longdash",color="darkgrey",linewidth=0.1)+
        geom_linerange(aes(x=0,ymin=1.5,ymax=2.0),linetype="longdash",color="darkgrey",linewidth=0.1)+
        geom_linerange(aes(x=0.3,ymin=1.5,ymax=2.0),linetype="longdash",color="darkgrey",linewidth=0.1)+
        geom_linerange(aes(x=0.3,ymin=1.5,ymax=2.2),linetype="longdash",color="darkgrey",linewidth=0.1)+
        geom_linerange(aes(x=0.67,ymin=1.5,ymax=2.2),linetype="longdash",color="darkgrey",linewidth=0.1)+
        geom_linerange(aes(xmin=0.3,xmax=0.67,y=2.2),linetype="longdash",color="darkgrey",linewidth=0.1)+
        geom_linerange(aes(x=0.67,ymin=1.5,ymax=2.4),linetype="longdash",color="darkgrey",linewidth=0.1)+
        geom_linerange(aes(x=1.2,ymin=1.5,ymax=2.4),linetype="longdash",color="darkgrey",linewidth=0.1)+
        geom_linerange(aes(xmin=0.67,xmax=1.2,y=2.4),linetype="longdash",color="darkgrey",linewidth=0.1)+
        geom_linerange(aes(x=0,ymin=0.7,ymax=1.5),linetype="longdash",color="darkgrey",linewidth=0.1)+
        geom_linerange(aes(x=0.1,ymin=0.7,ymax=1.5),linetype="longdash",color="darkgrey",linewidth=0.1)+
        geom_linerange(aes(xmin=0,xmax=0.1,y=0.7),linetype="longdash",color="darkgrey",linewidth=0.1)+
        geom_linerange(aes(x=0.1,ymin=0.7,ymax=1.5),linetype="longdash",color="darkgrey",linewidth=0.1)+
        geom_linerange(aes(xmin=0.1,xmax=0.67,y=0.7),linetype="longdash",color="darkgrey",linewidth=0.1)+
        geom_linerange(aes(x=0,ymin=0.2,ymax=0.7),linetype="longdash",color="darkgrey",linewidth=0.1)+
        geom_linerange(aes(x=0.67,ymin=0.2,ymax=0.7),linetype="longdash",color="darkgrey",linewidth=0.1)+
        geom_linerange(aes(xmin=0,xmax=0.67,y=0.2),linetype="longdash",color="darkgrey",linewidth=0.1)+
        geom_linerange(aes(x=0.67,ymin=0.6,ymax=1.5),linetype="longdash",color="darkgrey",linewidth=0.1)+
        geom_linerange(aes(x=1,ymin=0.6,ymax=1.5),linetype="longdash",color="darkgrey",linewidth=0.1)+
        geom_linerange(aes(xmin=0.67,xmax=1,y=0.6),linetype="longdash",color="darkgrey",linewidth=0.1)
      P1Element1_fill_legend <- ggscatter(datatmp0,'O/C','H/C',shape = 21,color = 'black',fill = 'Sample',alpha = 0.6,stroke=0.2)+
        scale_fill_manual(values = paletteer_dynamic("cartography::multi.pal", 20))+scale_size(range = c(0.8,3))+
        theme_prism(base_size = 7,border = TRUE,base_fontface = "plain",base_line_size = 8/50,base_rect_size = 8/50)+
        theme(legend.position = 'top',legend.key.height=unit(0.2,"line"),legend.margin = ggplot2::margin(2,0,0,0),legend.box.margin=ggplot2::margin(2,-10,0,-10),plot.margin = unit(c(0.2,1,0.6,0.6),'cm'))+
        guides(fill=guide_legend(override.aes = list(size=2)))
      guide_fill <- get_legend(P1Element1_fill_legend)
      P1Element1Up <- plot_grid(guide_fill,P1Element1_size_legend, 
                                nrow = 2, rel_heights = c(1, 10))
    }else{
      P1Element1 <- ggscatter(datatmp0,'O/C','H/C',fill = 'Sample', shape = 21,color = 'black',size = 1,alpha = 0.8,stroke=0.2)+gghighlight(Sample %in% highlight_filter)+
        scale_fill_manual(values = paletteer_dynamic("cartography::multi.pal", 20))+
        theme_prism(base_size = 7,border = TRUE,base_fontface = "plain",base_line_size = 8/50,base_rect_size = 8/50)+
        theme(legend.position = 'top',legend.key.height=unit(0.2,"line"),legend.margin = ggplot2::margin(2,0,0,0),legend.box.margin=ggplot2::margin(2,-10,0,-10),plot.margin = unit(c(0.2,1,0.6,0.6),'cm'))+
        guides(fill=guide_legend(override.aes = list(size=2)))
      P1Element1Up <- P1Element1+geom_hline(aes(yintercept=1.5),linetype="longdash",color="darkgrey",linewidth=0.1)+
        geom_linerange(aes(xmin=0,xmax=0.3,y=2.0),linetype="longdash",color="darkgrey",linewidth=0.1)+
        geom_linerange(aes(x=0,ymin=1.5,ymax=2.0),linetype="longdash",color="darkgrey",linewidth=0.1)+
        geom_linerange(aes(x=0.3,ymin=1.5,ymax=2.0),linetype="longdash",color="darkgrey",linewidth=0.1)+
        geom_linerange(aes(x=0.3,ymin=1.5,ymax=2.2),linetype="longdash",color="darkgrey",linewidth=0.1)+
        geom_linerange(aes(x=0.67,ymin=1.5,ymax=2.2),linetype="longdash",color="darkgrey",linewidth=0.1)+
        geom_linerange(aes(xmin=0.3,xmax=0.67,y=2.2),linetype="longdash",color="darkgrey",linewidth=0.1)+
        geom_linerange(aes(x=0.67,ymin=1.5,ymax=2.4),linetype="longdash",color="darkgrey",linewidth=0.1)+
        geom_linerange(aes(x=1.2,ymin=1.5,ymax=2.4),linetype="longdash",color="darkgrey",linewidth=0.1)+
        geom_linerange(aes(xmin=0.67,xmax=1.2,y=2.4),linetype="longdash",color="darkgrey",linewidth=0.1)+
        geom_linerange(aes(x=0,ymin=0.7,ymax=1.5),linetype="longdash",color="darkgrey",linewidth=0.1)+
        geom_linerange(aes(x=0.1,ymin=0.7,ymax=1.5),linetype="longdash",color="darkgrey",linewidth=0.1)+
        geom_linerange(aes(xmin=0,xmax=0.1,y=0.7),linetype="longdash",color="darkgrey",linewidth=0.1)+
        geom_linerange(aes(x=0.1,ymin=0.7,ymax=1.5),linetype="longdash",color="darkgrey",linewidth=0.1)+
        geom_linerange(aes(xmin=0.1,xmax=0.67,y=0.7),linetype="longdash",color="darkgrey",linewidth=0.1)+
        geom_linerange(aes(x=0,ymin=0.2,ymax=0.7),linetype="longdash",color="darkgrey",linewidth=0.1)+
        geom_linerange(aes(x=0.67,ymin=0.2,ymax=0.7),linetype="longdash",color="darkgrey",linewidth=0.1)+
        geom_linerange(aes(xmin=0,xmax=0.67,y=0.2),linetype="longdash",color="darkgrey",linewidth=0.1)+
        geom_linerange(aes(x=0.67,ymin=0.6,ymax=1.5),linetype="longdash",color="darkgrey",linewidth=0.1)+
        geom_linerange(aes(x=1,ymin=0.6,ymax=1.5),linetype="longdash",color="darkgrey",linewidth=0.1)+
        geom_linerange(aes(xmin=0.67,xmax=1,y=0.6),linetype="longdash",color="darkgrey",linewidth=0.1)
    }
  }
  return(P1Element1Up)
  
}

RelativeAbund <- function(FtMSDataList,SampleOrder=NA){
  if (any(is.na(SampleOrder))){
    SampleOrder <- names(FtMSDataList)
  }
  res_countElem <- tibble(Group=NA, n=NA, Relat=NA, Sample=NA)
  res_countKVD <- tibble(groupKVD=NA, n=NA, Relat=NA,Sample=NA)
  res_countKVD2 <- tibble(NewChemgroup=NA, n=NA, Relat=NA,Sample=NA)
  res_sumIntElem <- tibble(Group=NA, SumInt=NA, Relat=NA,Sample=NA)
  res_sumIntKVD <- tibble(groupKVD=NA, SumInt=NA, Relat=NA,Sample=NA)
  res_sumIntKVD2 <- tibble(NewChemgroup=NA, SumInt=NA, Relat=NA,Sample=NA)
  
  for (i in 1:length(FtMSDataList)) {
    data_countElem <- data.frame(Group=c("CHO","CHON","CHONP","CHONS","CHONSP","CHOP","CHOS","CHOSP"))
    data_countKVD <- data.frame(groupKVD=c("Aliphatic/protein","AromaticStructures","Carbohydrates","Lignin/CRAM","Lipid","Others","Tannin","Unsaturatedhydrocarbons"))
    data_countKVD2 <- data.frame(NewChemgroup=c("SaturatedFattyAcidLike","PeptideLike","UnsaturatedAliphaticLike","HighlyUnsaturatedLike","PolyphenolLike","CondensedAromaticLike"))
    data_sumIntElem <- data.frame(Group=c("CHO","CHON","CHONP","CHONS","CHONSP","CHOP","CHOS","CHOSP"))
    data_sumIntKVD <- data.frame(groupKVD=c("Aliphatic/protein","AromaticStructures","Carbohydrates","Lignin/CRAM","Lipid","Others","Tannin","Unsaturatedhydrocarbons"))
    data_sumIntKVD2 <- data.frame(NewChemgroup=c("SaturatedFattyAcidLike","PeptideLike","UnsaturatedAliphaticLike","HighlyUnsaturatedLike","PolyphenolLike","CondensedAromaticLike"))
    count_element <- FtMSDataList[[i]]%>%count(Group)
    count_element$Relat <- round(count_element$n*100/sum(count_element$n),2)
    data_countElem <-full_join(data_countElem, count_element,by="Group")
    data_countElem$Sample <- names(FtMSDataList)[i]
    count_KVD <- FtMSDataList[[i]]%>%count(groupKVD)
    count_KVD$Relat <- round(count_KVD$n*100/sum(count_KVD$n),2)
    data_countKVD <-full_join(data_countKVD, count_KVD,by="groupKVD")
    data_countKVD$Sample <- names(FtMSDataList)[i]
    count_KVD2 <- FtMSDataList[[i]]%>%count(NewChemgroup)
    count_KVD2$Relat <- round(count_KVD2$n*100/sum(count_KVD2$n),2)
    data_countKVD2 <-full_join(data_countKVD2, count_KVD2,by="NewChemgroup")
    data_countKVD2$Sample <- names(FtMSDataList)[i]
    
    sum_element <- aggregate(FtMSDataList[[i]]$Intensity,by=list(Group=FtMSDataList[[i]]$Group),sum)
    sum_element$Relat <- round(sum_element$x*100/sum(sum_element$x),2)
    data_sumIntElem <-full_join(data_sumIntElem, sum_element,by="Group")
    colnames(data_sumIntElem)[2]<-'SumInt'
    data_sumIntElem$Sample <- names(FtMSDataList)[i]
    sum_KVD <- aggregate(FtMSDataList[[i]]$Intensity,by=list(groupKVD=FtMSDataList[[i]]$groupKVD),sum)
    sum_KVD$Relat <- round(sum_KVD$x*100/sum(sum_KVD$x),2)
    data_sumIntKVD <-full_join(data_sumIntKVD, sum_KVD,by="groupKVD")
    colnames(data_sumIntKVD)[2]<-'SumInt'
    data_sumIntKVD$Sample <- names(FtMSDataList)[i]
    sum_KVD2 <- aggregate(FtMSDataList[[i]]$Intensity,by=list(NewChemgroup=FtMSDataList[[i]]$NewChemgroup),sum)
    sum_KVD2$Relat <- round(sum_KVD2$x*100/sum(sum_KVD2$x),2)
    data_sumIntKVD2 <-full_join(data_sumIntKVD2, sum_KVD2,by="NewChemgroup")
    colnames(data_sumIntKVD2)[2]<-'SumInt'
    data_sumIntKVD2$Sample <- names(FtMSDataList)[i]
    res_countElem <- rbind(res_countElem, data_countElem)
    res_countKVD <- rbind(res_countKVD, data_countKVD)
    res_countKVD2 <- rbind(res_countKVD2, data_countKVD2)
    res_sumIntElem <- rbind(res_sumIntElem, data_sumIntElem)
    res_sumIntKVD <- rbind(res_sumIntKVD, data_sumIntKVD)
    res_sumIntKVD2 <- rbind(res_sumIntKVD2, data_sumIntKVD2)
  }
  res_countElem <-res_countElem[2:nrow(res_countElem),]
  p1<-ggpar(ggbarplot(res_countElem,'Sample','Relat', fill = 'Group',size =8/70,alpha=0.6, width = 0.4, order = SampleOrder,label=TRUE, lab.size = 1,lab.pos = 'in')+
              xlab(NULL)+ylab("Count proportion (%)")+
              theme(axis.text.x = element_text(vjust=1,size = 5.5))+
              #geom_label_repel(aes(label = Relat),segment.size=0.1, label.r=0.02,fontface = "plain",label.size = 8/80,size=0.6, label.padding = 8/90,box.padding = 8/90)+
              theme_prism(base_size = 9,border = T,base_fontface = "plain",base_line_size = 8/50,base_rect_size = 8/70)+
              guides(fill=guide_legend(keywidth=0.4,keyheight=0.3,size=8/70)),legend = 'top',x.text.angle = 20)
  #ggsave("./output/res_countElem.pdf",width=8.5,height=7,dpi=300,units = "cm")
  p2<-ggpar(ggbarplot(res_countElem,'Sample','n',fill = 'Group',size =8/70,alpha=0.6, width = 0.4, order = SampleOrder,label = TRUE,lab.size = 1,lab.pos = 'in')+
              xlab(NULL)+ylab("Count")+
              theme(axis.text.x = element_text(vjust=1,size = 5.5))+
              #geom_label_repel(aes(label = n),segment.size=0.1, label.r=0.02,fontface = "plain",label.size = 8/80,size=0.6, label.padding = 8/90,box.padding = 8/90)+
              theme_prism(base_size = 9,border = T,base_fontface = "plain",base_line_size = 8/50,base_rect_size = 8/70)+
              guides(fill=guide_legend(keywidth=0.4,keyheight=0.3,size=8/70,direction = 'horizontal')),legend = 'top',x.text.angle = 20)
  #ggsave("./output/res_countNumElem.pdf",width=8.5,height=7,dpi=300,units = "cm")
  p1_t <- p1+ theme(axis.text.x=element_blank())
  p1_2 <- guide_area()/p1_t/p2+plot_layout(heights=c(0.1,2,2), guides = 'collect')
  #ggsave("./output/res_countNumProElem.pdf",width=9,height=8,dpi=300,units = "cm")
  res_countKVD <-res_countKVD[2:nrow(res_countKVD),]
  p3 <- ggpar(ggbarplot(res_countKVD,'Sample','Relat', fill = 'groupKVD',size =8/70,alpha=0.6, width = 0.4, order = SampleOrder,label = TRUE,lab.size = 1,lab.pos = 'in')+
                fill_palette('jco')+xlab(NULL)+ylab("Count proportion (%)")+
                theme(axis.text.x = element_text(vjust=1,size = 5.5))+
                theme_prism(base_size = 9,border = T,base_fontface = "plain",base_line_size = 8/50,base_rect_size = 8/70)+
                guides(fill=guide_legend(keywidth=0.4,keyheight=0.3,size=8/70)),legend = 'top',x.text.angle = 20)
  #ggsave("./output/res_countKVD.pdf",width=8.5,height=7,dpi=300,units = 'cm')
  p4 <- ggpar(ggbarplot(res_countKVD,'Sample','n', fill = 'groupKVD',size =8/70,alpha=0.6, width = 0.4, order = SampleOrder,label = TRUE,lab.size = 1,lab.pos = 'in')+
                fill_palette('jco')+xlab(NULL)+ylab("Count")+
                theme(axis.text.x = element_text(vjust=1,size = 5.5))+
                theme_prism(base_size = 9,border = T,base_fontface = "plain",base_line_size = 8/50,base_rect_size = 8/70)+
                guides(fill=guide_legend(keywidth=0.4,keyheight=0.3,size=8/70)),legend = 'top',x.text.angle = 20)
  #ggsave("./output/res_countNumKVD.pdf",width=8.5,height=7,dpi=300,units = 'cm')
  p3_t <- p3+ theme(axis.text.x=element_blank())
  p3_4 <- guide_area()/p3_t/p4+plot_layout(heights=c(0.1,2,2), guides = 'collect')
  #ggsave("./output/res_countNumProKVD.pdf",width=9,height=8,dpi=300,units = "cm")
  
  res_countKVD2 <-res_countKVD2[2:nrow(res_countKVD2),]
  p5 <- ggpar(ggbarplot(res_countKVD2,'Sample','Relat', fill = 'NewChemgroup',size =8/70,alpha=0.6, width = 0.4, order = SampleOrder,label = TRUE,lab.size = 1,lab.pos = 'in')+
                fill_palette('jco')+xlab(NULL)+ylab("Count proportion (%)")+
                theme(axis.text.x = element_text(vjust=1,size = 5.5))+
                theme_prism(base_size = 9,border = T,base_fontface = "plain",base_line_size = 8/50,base_rect_size = 8/70)+
                guides(fill=guide_legend(keywidth=0.4,keyheight=0.3,size=8/70)),legend = 'top',x.text.angle = 20)
  #ggsave("./output/res_countKVD2.pdf",width=8.5,height=7,dpi=300,units = 'cm')
  p6 <- ggpar(ggbarplot(res_countKVD2,'Sample','n', fill = 'NewChemgroup',size =8/70,alpha=0.6, width = 0.4, order = SampleOrder,label = TRUE,lab.size = 1,lab.pos = 'in')+
                fill_palette('jco')+xlab(NULL)+ylab("Count")+
                theme(axis.text.x = element_text(vjust=1,size = 5.5))+
                theme_prism(base_size = 9,border = T,base_fontface = "plain",base_line_size = 8/50,base_rect_size = 8/70)+
                guides(fill=guide_legend(keywidth=0.4,keyheight=0.3,size=8/70)),legend = 'top',x.text.angle = 20)
  #ggsave("./output/res_countNumKVD2.pdf",width=8.5,height=7,dpi=300,units = 'cm')
  p5_t <- p5+ theme(axis.text.x=element_blank())
  p5_6 <- guide_area()/p5_t/p6+plot_layout(heights=c(0.1,2,2), guides = 'collect')
  #ggsave("./output/res_countNumProKVD2.pdf",width=9,height=8,dpi=300,units = "cm")
  
  res_sumIntElem <- res_sumIntElem[2:nrow(res_sumIntElem),]
  p7 <- ggpar(ggbarplot(res_sumIntElem,'Sample','Relat', fill = 'Group', size =8/70,alpha=0.6, width = 0.4, order = SampleOrder,label = TRUE,lab.size = 1,lab.pos = 'in')+
                xlab(NULL)+ylab("Relative abundance (%)")+
                theme(axis.text.x = element_text(vjust=1,size = 5.5))+
                #geom_label_repel(aes(label = Relat),segment.size=0.1, label.r=0.02,fontface = "plain",label.size = 8/80,size=0.6, label.padding = 8/90,box.padding = 8/90)+
                theme_prism(base_size = 9,border = T,base_fontface = "plain",base_line_size = 8/50,base_rect_size = 8/70)+
                guides(fill=guide_legend(keywidth=0.4,keyheight=0.3,size=8/70,direction = 'horizontal')),legend = 'top',x.text.angle = 20)
  #ggsave("./output/res_sumIntElem.pdf",width=8.5,height=7,dpi=300,units = 'cm')
  res_sumIntKVD <- res_sumIntKVD[2:nrow(res_sumIntKVD),]
  p8 <- ggpar(ggbarplot(res_sumIntKVD,'Sample','Relat', fill = 'groupKVD', size =8/70,alpha=0.6, width = 0.4, order = SampleOrder,label = TRUE,lab.size = 1,lab.pos = 'in')+
                fill_palette('npg')+xlab(NULL)+ylab("Relative abundance (%)")+
                theme(axis.text.x = element_text(vjust=1,size = 5.5))+
                #geom_label_repel(aes(label = Relat),segment.size=0.1, label.r=0.02,fontface = "plain",label.size = 8/80,size=0.6, label.padding = 8/90,box.padding = 8/90)+
                theme_prism(base_size = 9,border = T,base_fontface = "plain",base_line_size = 8/50,base_rect_size = 8/70)+
                guides(fill=guide_legend(keywidth=0.4,keyheight=0.3,size=8/70,direction = 'horizontal')),legend = 'top',x.text.angle = 20)
  #ggsave("./output/res_sumIntKVD.pdf",width=8.5,height=7,dpi=300,units = 'cm')
  res_sumIntKVD2 <- res_sumIntKVD2[2:nrow(res_sumIntKVD2),]
  p9 <- ggpar(ggbarplot(res_sumIntKVD2,'Sample','Relat', fill = 'NewChemgroup', size =8/70,alpha=0.6, width = 0.4, order = SampleOrder,label = TRUE,lab.size = 1,lab.pos = 'in')+
                fill_palette('npg')+xlab(NULL)+ylab("Relative abundance (%)")+
                theme(axis.text.x = element_text(vjust=1,size = 5.5))+
                #geom_label_repel(aes(label = Relat),segment.size=0.1, label.r=0.02,fontface = "plain",label.size = 8/80,size=0.6, label.padding = 8/90,box.padding = 8/90)+
                theme_prism(base_size = 9,border = T,base_fontface = "plain",base_line_size = 8/50,base_rect_size = 8/70)+
                guides(fill=guide_legend(keywidth=0.4,keyheight=0.3,size=8/70,direction = 'horizontal')),legend = 'top',x.text.angle = 20)
  #ggsave("./output/res_sumIntKVD2.pdf",width=8.5,height=7,dpi=300,units = 'cm')
  #write.csv(res_countElem,"./Output/res_countElem.csv",row.names = F)
  #write.csv(res_countKVD,"./Output/res_countKVD.csv",row.names = F)
  #write.csv(res_countKVD2,"./Output/res_countKVD2.csv",row.names = F)
  #write.csv(res_sumIntElem,"./Output/res_sumIntElem.csv",row.names = F)
  #write.csv(res_sumIntKVD,"./Output/res_sumIntKVD.csv",row.names = F)
  #write.csv(res_sumIntKVD2,"./Output/res_sumIntKVD2.csv",row.names = F)
  res_RelaAbund <- list(p1,p2,p1_2,p3,p4,p3_4,p5,p6,p5_6,p7,p8,p9)
  names(res_RelaAbund) <- c("p1","p2","p1_2","p3","p4","p3_4","p5","p6","p5_6","p7","p8","p9")
  return(res_RelaAbund)
}



XcGraph <- function(FtMSDataList){
  
  Fig_ls <- list()
  Fig2_ls <- list()
  FigGroup_ls <- list()
  FigCount_ls <- list()
  XcCount <- tibble(total = NA, XcLower2.5=NA, Xc0=NA,Xc2.5=NA, Xc2.7=NA, Xc2.8=NA)
  for (i in 1:length(FtMSDataList)) {
    nXc0 <- nrow(FtMSDataList[[i]])
    tempCountXc1 <- filter(FtMSDataList[[i]],Xc < 2.5)
    nXc1 <- nrow(tempCountXc1)
    tempCountXc2 <- filter(FtMSDataList[[i]],Xc == 0)
    nXc2 <- nrow(tempCountXc2)
    tempCountXc3 <- filter(FtMSDataList[[i]],Xc == 2.5)
    nXc3 <- nrow(tempCountXc3)
    tempCountXc4 <- filter(FtMSDataList[[i]],Xc == 2.7143)
    nXc4 <- nrow(tempCountXc4)
    tempCountXc5 <- filter(FtMSDataList[[i]],Xc == 2.8)
    nXc5 <- nrow(tempCountXc5)
    XcCount <- tibble(total = nXc0, XcLower2.5=nXc1,Xc0=nXc2,Xc2.5=nXc3, Xc2.7=nXc4, Xc2.8=nXc5)
    #write.csv(XcCount,sprintf("./output/Xc/%s_XcCount.csv",names(FtMSDataList)[i]),row.names = FALSE)
    
    #FormuTransDf <- FormuTrans(FtMSDataList[[i]]["Formula"])
    #FtMSDataList[[i]]$CarbonNumber <- FtMSDataList[[i]]$C+FtMSDataList[[i]]$`18C`
    FtMSDataList[[i]] <- FtMSDataList[[i]]%>% mutate(CarbonNumber=C+`13C`)
    FtMSDataList[[i]]$XctoC <- paste(FtMSDataList[[i]]$Xc,FtMSDataList[[i]]$CarbonNumber,sep = "-")
    FtMSDataList[[i]] <- FtMSDataList[[i]]%>% add_count(XctoC)
    FtMSDataList[[i]]$Xc <- round(FtMSDataList[[i]]$Xc,digits = 4)
    FtMSDataList[[i]]$Xc <- as.character(FtMSDataList[[i]]$Xc)
    Fig_ls[[i]]<-ggpar(ggscatter(FtMSDataList[[i]],'CarbonNumber','Xc',size = 'ReInt',fill = 'groupKVD',shape=21, color = 'black',alpha=0.8,stroke=0.2)+
                         scale_size(range = c(1,3.5))+fill_palette('jco') +
                         theme_prism(base_size = 5.5,border = TRUE,base_fontface = "plain",base_rect_size = 8/50,base_line_size = 8/50)+
                         theme(legend.position = 'top',legend.box.spacing = unit(0.001,'cm'),legend.box.margin = ggplot2::margin(0.01, .01, .01, .01, "cm"))+
                         guides(fill=guide_legend(override.aes = list(size=2))),xticks.by = 2)
    Fig2_ls[[i]]<-ggpar(ggscatter(FtMSDataList[[i]],'CarbonNumber','Xc',size = 'ReInt',fill = 'NewChemgroup',shape=21, color = 'black',alpha=0.8,stroke=0.2)+
                          scale_size(range = c(1,3.5))+fill_palette('jco') +
                          theme_prism(base_size = 5.5,border = TRUE,base_fontface = "plain",base_rect_size = 8/50,base_line_size = 8/50)+
                          theme(legend.position = 'top',legend.box.spacing = unit(0.001,'cm'),legend.box.margin = ggplot2::margin(0.01, .01, .01, .01, "cm"))+
                          guides(fill=guide_legend(override.aes = list(size=2))),xticks.by = 2)
    FigGroup_ls[[i]]<-ggpar(ggscatter(FtMSDataList[[i]],'CarbonNumber','Xc',size = 'ReInt',fill = 'Group',shape=21, color = 'black',alpha=0.8,stroke=0.2)+
                              scale_size(range = c(1,3.5))+
                              fill_palette('jco') +theme_prism(base_size = 5.5,border = TRUE,base_fontface = "plain",base_rect_size = 8/50,base_line_size = 8/50)+
                              theme(legend.position = 'top',legend.box.spacing = unit(0.001,'cm'),legend.box.margin = ggplot2::margin(0.01, .01, .01, .01, "cm"))+
                              guides(fill=guide_legend(override.aes = list(size=2))),xticks.by = 2)
    FigCount_ls[[i]]<-ggpar(ggscatter(FtMSDataList[[i]],'CarbonNumber','Xc',size = 1.5,shape=15, color = 'n',alpha=1)+
                              geom_text(aes(label=n),color="white",size=1) +
                              scale_color_gradientn(colors=rev(paletteer_c("ggthemes::Classic Orange-Blue", 15)))+
                              theme_prism(base_size = 5.3,border = TRUE,base_fontface = 'plain',base_rect_size = 8/50,base_line_size = 8/50),xticks.by = 2,xlim = c(5,50))+
      theme(legend.key.width = unit(0.2,'cm'),legend.spacing = unit(0.01,'cm'),legend.box.spacing =unit(0.01,'cm') )
    names(Fig_ls)[i]<- names(FtMSDataList)[i]
    names(Fig2_ls)[i]<- names(FtMSDataList)[i]
    names(FigGroup_ls)[i]<- names(FtMSDataList)[i]
    names(FigCount_ls)[i]<- names(FtMSDataList)[i]
    #ggsave(sprintf("./output/Xc/%s_Xc.PDF",names(Fig_ls[i])),Fig_ls[[i]],dpi=300,width = 8.5,height = 7,units = "cm")
    #ggsave(sprintf("./output/Xc/%s_Xc2.PDF",names(Fig2_ls[i])),Fig2_ls[[i]],dpi=300,width = 8.5,height = 7,units = "cm")
    #ggsave(sprintf("./output/Xc/%s_GroupXc.PDF",names(FigGroup_ls[i])),FigGroup_ls[[i]],dpi=300,width = 8.5,height = 7,units = "cm")
    #ggsave(sprintf("./output/Xc/%s_CountXc.PDF",names(FigCount_ls[i])),FigCount_ls[[i]],dpi=300,width = 8.5,height = 6.6,units = "cm")
  }
  return(FigCount_ls[[1]])
}


FtMsUpsetData <- function(FtMSDataList){
  FormularList <- list()
  for (i in 1:length(FtMSDataList)) {
    FormularList$data <- FtMSDataList[[i]]$Formula
    if(str_detect(names(FtMSDataList)[i],".xlsx")){
      names(FormularList)[i]<- str_remove(names(FtMSDataList)[i],".xlsx")
    }else{
      names(FormularList)[i] <-names(FtMSDataList)[i]
    }
    #names(FormularList)[i]<- names(FtMSDataList)[i]
  }
  return(FormularList)
}

PCAData <- function(FtMsDf,Select_TopRANum=1000){
  RADataLS <- list()
  for (i in 1:length(FtMsDf)) {
    if(any(colnames(FtMsDf[[i]]) %in% 'RA(%)')){
      FtMsDf[[i]]$RA <- FtMsDf[[i]]$`RA(%)`
    }
    
    FtMsDf[[i]] <- FtMsDf[[i]][order(FtMsDf[[i]]$RA,decreasing = T),]
    FtMsDf[[i]] <-FtMsDf[[i]][1:floor(0.7*nrow(FtMsDf[[i]])),]
    dupid <- which(duplicated(FtMsDf[[i]]$Formula))
    nodupdata <- FtMsDf[[i]][!FtMsDf[[i]]$Formula %in% FtMsDf[[i]]$Formula[dupid],]
    RAdata <- nodupdata[c('Formula','RA')]
    for (id in dupid) {
      RAdatatemp <- FtMsDf[[i]][FtMsDf[[i]]$Formula%in%FtMsDf[[i]]$Formula[id],][c('Formula','RA')]
      #RAdatatemp$RA <- sum(RAdatatemp$RA)
      RAdatatemp$RA <- max(RAdatatemp$RA)
      RAdatatemp <- unique(RAdatatemp)
      RAdata <- rbind(RAdata,RAdatatemp)
    }
    RADataLS[[i]]<-RAdata
    if(str_detect(names(FtMsDf)[i],".xlsx")){
      names(RADataLS)[i]<- str_remove(names(FtMsDf)[i],".xlsx")
    }else{
      names(RADataLS)[i]<-names(FtMsDf)[i]
    }
  }
  if(is.na(Select_TopRANum)){
    PCAdata <- full_join(RADataLS[[1]][c('Formula','RA')],RADataLS[[2]][c('Formula','RA')],by='Formula')
    PCAdata <- unique(PCAdata)
    if (length(RADataLS)>2){
      for (i in 3:length(RADataLS)) {
        PCAdata <- full_join(PCAdata,RADataLS[[i]][c('Formula','RA')],by='Formula')
        PCAdata <- unique(PCAdata)
      }
    }
  }else{
    RADataLS[[1]] <- RADataLS[[1]][order( RADataLS[[1]]$RA,decreasing = T),][1:Select_TopRANum,]
    RADataLS[[2]] <- RADataLS[[2]][order( RADataLS[[2]]$RA,decreasing = T),][1:Select_TopRANum,]
    PCAdata <- full_join(RADataLS[[1]][c('Formula','RA')],RADataLS[[2]][c('Formula','RA')],by='Formula')
    PCAdata <- unique(PCAdata)
    if (length(RADataLS)>2){
      for (i in 3:length(RADataLS)) {
        RADataLS[[i]] <- RADataLS[[i]][order( RADataLS[[i]]$RA,decreasing = T),][1:Select_TopRANum,]
        PCAdata <- full_join(PCAdata,RADataLS[[i]][c('Formula','RA')],by='Formula')
        PCAdata <- unique(PCAdata)
      }
    }
  }
  colnames(PCAdata)<-c('Formula',names(RADataLS)) 
  PCAdata <- unique(PCAdata)
  rownames(PCAdata)<-PCAdata$Formula
  PCAdata<-PCAdata[,2:ncol(PCAdata)]
  PCAdata[is.na(PCAdata)] <- 0
  PCAdata <- t(PCAdata)
  res.pca <- PCA(PCAdata, graph = FALSE)
  return(res.pca)
}



HeatmapData <- function(FtMsDf,group=NA,topnum=60,Select_TopRANum=1000){
  #把因为有同位素的相同分子式的RA进行合并
  RADataLS <- list()
  for (i in 1:length(FtMsDf)) {
    if(any(colnames(FtMsDf[[i]]) %in% 'RA(%)')){
      FtMsDf[[i]]$RA <- FtMsDf[[i]]$`RA(%)`
    }
    FtMsDf[[i]] <- FtMsDf[[i]][order(FtMsDf[[i]]$RA,decreasing = T),]
    FtMsDf[[i]] <-FtMsDf[[i]][1:floor(0.7*nrow(FtMsDf[[i]])),]
    dupid <- which(duplicated(FtMsDf[[i]]$Formula))
    nodupdata <- FtMsDf[[i]][!FtMsDf[[i]]$Formula %in% FtMsDf[[i]]$Formula[dupid],]
    RAdata <- nodupdata[c('Formula','RA')]
    for (id in dupid) {
      RAdatatemp <- FtMsDf[[i]][FtMsDf[[i]]$Formula%in%FtMsDf[[i]]$Formula[id],][c('Formula','RA')]
      RAdatatemp$RA <- max(RAdatatemp$RA)
      RAdatatemp <- unique(RAdatatemp)
      RAdata <- rbind(RAdata,RAdatatemp)
    }
    RADataLS[[i]]<-RAdata
    if(str_detect(names(FtMsDf)[i],".xlsx")){
      names(RADataLS)[i]<- str_remove(names(FtMsDf)[i],".xlsx")
    }else{
      names(RADataLS)[i]<-names(FtMsDf)[i]
    }
  }
  
  if(is.na(Select_TopRANum)){
    PCAdata <- full_join(RADataLS[[1]][c('Formula','RA')],RADataLS[[2]][c('Formula','RA')],by='Formula')
    if (length(RADataLS)>2){
      for (i in 3:length(RADataLS)) {
        PCAdata <- full_join(PCAdata,RADataLS[[i]][c('Formula','RA')],by='Formula')
        PCAdata <- unique(PCAdata)
      }
    }
  }else{
    RADataLS[[1]] <- RADataLS[[1]][order( RADataLS[[1]]$RA,decreasing = T),][1:Select_TopRANum,]
    RADataLS[[2]] <- RADataLS[[2]][order( RADataLS[[2]]$RA,decreasing = T),][1:Select_TopRANum,]
    PCAdata <- full_join(RADataLS[[1]][c('Formula','RA')],RADataLS[[2]][c('Formula','RA')],by='Formula')
    if (length(RADataLS)>2){
      for (i in 3:length(RADataLS)) {
        RADataLS[[i]] <- RADataLS[[i]][order( RADataLS[[i]]$RA,decreasing = T),][1:Select_TopRANum,]
        PCAdata <- full_join(PCAdata,RADataLS[[i]][c('Formula','RA')],by='Formula')
        PCAdata <- unique(PCAdata)
      }
    }
  }
  
  
  colnames(PCAdata)<-c('Formula',names(RADataLS)) 
  #PCAdata <- as.data.frame(PCAdata)
  PCAdata <- unique(PCAdata)
  
  rownames<-as.character(PCAdata$Formula)
  PCAdata<-PCAdata[,2:ncol(PCAdata)]
  
  rownames(PCAdata) <- as.character(rownames[1:length(rownames)])
  PCAdata[is.na(PCAdata)] <- 0
  PCAdata$sum <- rowSums(PCAdata)
  PCAdata <- PCAdata[order(PCAdata$sum,decreasing = T),]
  #PCAdata <- as.data.frame(t(PCAdata))
  TopPCAdata <- log2(PCAdata[1:topnum,1:(ncol(PCAdata)-1)]+1)
  TopPCAdata <- as.data.frame(TopPCAdata)
  return(TopPCAdata)
}

