WeightedIndex <- function(FtMSDataList){
  data_DBEw <- data.frame()
  `data_DBE-Ow` <- data.frame()
  `data_(DBE-0.5O)/Cw` <- data.frame()
  `data_(DBE-O)/Cw` <- data.frame()
  data_MWw <- data.frame()
  `data_O/Cw` <- data.frame()
  `data_H/Cw` <- data.frame()
  `data_N/Cw` <- data.frame()
  `data_S/Cw` <- data.frame()
  `data_P/Cw` <- data.frame()
  `data_DBE/Cw` <- data.frame()
  data_AImodw <- data.frame()
  data_NOSCw <- data.frame()
  data_DeltaGw <- data.frame()
  data_MLB <- data.frame()
  data_lambdaw <- data.frame()
  data_diversity_shannon <- data.frame()
  data_diversity_simpson <- data.frame()
  data_diversity_func_NOSC <- data.frame()
  data_diversity_func_Unsaturation <- data.frame()
  data_diversity_func_ThermoDyn <- data.frame()
  for (i in 1:length(FtMSDataList)) {
    
    if (any(colnames(FtMSDataList[[i]]) %in% "D")){
      if (!any(str_detect(colnames(FtMSDataList[[i]]),'33S'))){
        FtMSDataList[[i]]$MolecularWeight <- FtMSDataList[[i]]$C*12+FtMSDataList[[i]]$`13C`*13.003355+FtMSDataList[[i]]$H*1.007825+FtMSDataList[[i]]$`2H`*2.014102+ FtMSDataList[[i]]$`D`*2.014102+FtMSDataList[[i]]$O*15.994915+FtMSDataList[[i]]$`18O`*17.999160+FtMSDataList[[i]]$N*14.003074+FtMSDataList[[i]]$`15N`*15.000109+FtMSDataList[[i]]$S*31.972071+FtMSDataList[[i]]$`34S`*33.967867+FtMSDataList[[i]]$P*30.973762
      }else{
        FtMSDataList[[i]]$MolecularWeight <- FtMSDataList[[i]]$C*12+FtMSDataList[[i]]$`13C`*13.003355+FtMSDataList[[i]]$H*1.007825+FtMSDataList[[i]]$`2H`*2.014102+ FtMSDataList[[i]]$`D`*2.014102+FtMSDataList[[i]]$O*15.994915+FtMSDataList[[i]]$`18O`*17.999160+FtMSDataList[[i]]$N*14.003074+FtMSDataList[[i]]$`15N`*15.000109+FtMSDataList[[i]]$S*31.972071+FtMSDataList[[i]]$`33S`*32.971458+FtMSDataList[[i]]$`34S`*33.967867+FtMSDataList[[i]]$P*30.973762
      }
      
    }else{
      if (!any(str_detect(colnames(FtMSDataList[[i]]),'33S'))){
        FtMSDataList[[i]]$MolecularWeight <- FtMSDataList[[i]]$C*12+FtMSDataList[[i]]$`13C`*13.003355+FtMSDataList[[i]]$H*1.007825+FtMSDataList[[i]]$`2H`*2.014102+FtMSDataList[[i]]$O*15.994915+FtMSDataList[[i]]$`18O`*17.999160+FtMSDataList[[i]]$N*14.003074+FtMSDataList[[i]]$`15N`*15.000109+FtMSDataList[[i]]$S*31.972071+FtMSDataList[[i]]$`34S`*33.967867+FtMSDataList[[i]]$P*30.973762
      }else{
        FtMSDataList[[i]]$MolecularWeight <- FtMSDataList[[i]]$C*12+FtMSDataList[[i]]$`13C`*13.003355+FtMSDataList[[i]]$H*1.007825+FtMSDataList[[i]]$`2H`*2.014102+FtMSDataList[[i]]$O*15.994915+FtMSDataList[[i]]$`18O`*17.999160+FtMSDataList[[i]]$N*14.003074+FtMSDataList[[i]]$`15N`*15.000109+FtMSDataList[[i]]$S*31.972071+FtMSDataList[[i]]$`33S`*32.971458+FtMSDataList[[i]]$`34S`*33.967867+FtMSDataList[[i]]$P*30.973762
      }
      
    }
    DBEw <- as.data.frame(sum(FtMSDataList[[i]]$DBE*FtMSDataList[[i]]$Intensity)/sum(FtMSDataList[[i]]$Intensity))
    colnames(DBEw)[1]<-"DBEw"
    data_DBEw <- rbind(data_DBEw,DBEw)
    `DBE-Ow` <- as.data.frame(sum(FtMSDataList[[i]]$`DBE-O`*FtMSDataList[[i]]$Intensity)/sum(FtMSDataList[[i]]$Intensity))
    colnames(`DBE-Ow`)[1]<-"DBE-Ow"
    `data_DBE-Ow` <- rbind(`data_DBE-Ow`,`DBE-Ow`)
    MWw <- as.data.frame(sum(FtMSDataList[[i]]$MolecularWeight*FtMSDataList[[i]]$Intensity)/sum(FtMSDataList[[i]]$Intensity))
    colnames(MWw)[1]<-"MWw"
    data_MWw <- rbind(data_MWw,MWw)
    `O/Cw` <- as.data.frame(sum(FtMSDataList[[i]]$`O/C`*FtMSDataList[[i]]$Intensity)/sum(FtMSDataList[[i]]$Intensity))
    colnames(`O/Cw`)[1]<-"O/Cw"
    `data_O/Cw` <- rbind(`data_O/Cw`,`O/Cw`)
    `H/Cw` <- as.data.frame(sum(FtMSDataList[[i]]$`H/C`*FtMSDataList[[i]]$Intensity)/sum(FtMSDataList[[i]]$Intensity))
    colnames(`H/Cw`)[1]<-"H/Cw"
    `data_H/Cw` <- rbind(`data_H/Cw`,`H/Cw`)
    `N/Cw` <- as.data.frame(sum(FtMSDataList[[i]]$`N/C`*FtMSDataList[[i]]$Intensity)/sum(FtMSDataList[[i]]$Intensity))
    colnames(`N/Cw`)[1]<-"N/Cw"
    `data_N/Cw` <- rbind(`data_N/Cw`,`N/Cw`)
    `S/Cw` <- as.data.frame(sum(FtMSDataList[[i]]$`S/C`*FtMSDataList[[i]]$Intensity)/sum(FtMSDataList[[i]]$Intensity))
    colnames(`S/Cw`)[1]<-"S/Cw"
    `data_S/Cw` <- rbind(`data_S/Cw`,`S/Cw`)
    `P/Cw` <- as.data.frame(sum(FtMSDataList[[i]]$`P/C`*FtMSDataList[[i]]$Intensity)/sum(FtMSDataList[[i]]$Intensity))
    colnames(`P/Cw`)[1]<-"P/Cw"
    `data_P/Cw` <- rbind(`data_P/Cw`,`P/Cw`)
    `DBE/Cw` <- as.data.frame(sum(FtMSDataList[[i]]$`DBE/C`*FtMSDataList[[i]]$Intensity)/sum(FtMSDataList[[i]]$Intensity))
    colnames(`DBE/Cw`)[1]<-"DBE/Cw"
    `data_DBE/Cw` <- rbind(`data_DBE/Cw`,`DBE/Cw`)
    `(DBE-O)/Cw` <- as.data.frame(sum(FtMSDataList[[i]]$`(DBE-O)/C`*FtMSDataList[[i]]$Intensity)/sum(FtMSDataList[[i]]$Intensity))
    colnames(`(DBE-O)/Cw`)[1]<-"(DBE-O)/Cw"
    `data_(DBE-O)/Cw` <- rbind(`data_(DBE-O)/Cw`,`(DBE-O)/Cw`)
    `(DBE-0.5O)/Cw` <- as.data.frame(sum(FtMSDataList[[i]]$`(DBE-0.5O)/C`*FtMSDataList[[i]]$Intensity)/sum(FtMSDataList[[i]]$Intensity))
    colnames(`(DBE-0.5O)/Cw`)[1]<-"(DBE-0.5O)/Cw"
    `data_(DBE-0.5O)/Cw` <- rbind(`data_(DBE-0.5O)/Cw`,`(DBE-0.5O)/Cw`)
    AImodw <- as.data.frame(sum(FtMSDataList[[i]]$AImod*FtMSDataList[[i]]$Intensity)/sum(FtMSDataList[[i]]$Intensity))
    colnames(AImodw)[1]<-"AImodw"
    data_AImodw <- rbind(data_AImodw,AImodw)
    NOSCw <- as.data.frame(sum(FtMSDataList[[i]]$NOSC*FtMSDataList[[i]]$Intensity)/sum(FtMSDataList[[i]]$Intensity))
    colnames(NOSCw)[1]<-"NOSCw"
    data_NOSCw <- rbind(data_NOSCw,NOSCw)
    DeltaGw <- as.data.frame(sum(FtMSDataList[[i]]$DeltaG*FtMSDataList[[i]]$Intensity)/sum(FtMSDataList[[i]]$Intensity))
    colnames(DeltaGw)[1]<-"DeltaGw"
    data_DeltaGw <- rbind(data_DeltaGw,DeltaGw)
    if(any(colnames(FtMSDataList[[i]]) %in% "lambda")){
      datatmp <- FtMSDataList[[i]]%>%filter(!is.na(lambda))
      lambdaw <- as.data.frame(sum(datatmp$lambda*datatmp$Intensity)/sum(datatmp$Intensity))
      colnames(lambdaw)[1]<-"lambdaw"
      data_lambdaw <- rbind(data_lambdaw,lambdaw)
    }
    countMLB <- FtMSDataList[[i]] %>% count(MLB)
    countMLB$Prop <- countMLB$n*100/sum(countMLB$n)
    MLBtemp <- as.data.frame(countMLB$Prop[1])
    colnames(MLBtemp)[1]<-'MLBindex'
    data_MLB <- rbind(data_MLB,MLBtemp)
    diversity_shannon <- as.data.frame(vegan::diversity(FtMSDataList[[i]]$RA,index = 'shannon'))
    colnames(diversity_shannon)[1] <- "diversity_shannon"
    data_diversity_shannon <- rbind(data_diversity_shannon,diversity_shannon)
    
    diversity_simpson <- as.data.frame(vegan::diversity(FtMSDataList[[i]]$RA,index = 'simpson'))
    colnames(diversity_simpson)[1] <- "diversity_simpson"
    data_diversity_simpson <- rbind(data_diversity_simpson,diversity_simpson)
    mat <- FtMSDataList[[i]] %>% select('Measured m/z',RA) %>% column_to_rownames(var = "Measured m/z") %>% t()
    rx_trait <- FtMSDataList[[i]] %>% select('Measured m/z',NOSC) %>% column_to_rownames(var = 'Measured m/z')
    diversity_func_NOSC <- as.data.frame(as.numeric(rao.diversity(mat,rx_trait)$FunRao))
    colnames(diversity_func_NOSC)[1] <- "diversity_Func_NOSC"
    data_diversity_func_NOSC <- rbind(data_diversity_func_NOSC,diversity_func_NOSC)
    
    ia_trait <- FtMSDataList[[i]] %>% select('Measured m/z',DBE,AImod) %>% column_to_rownames(var = 'Measured m/z')
    diversity_func_Unsaturation <- as.data.frame(as.numeric(rao.diversity(mat,ia_trait)$FunRao))
    colnames(diversity_func_Unsaturation)[1] <- "diversity_Func_Unsaturation"
    data_diversity_func_Unsaturation <- rbind(data_diversity_func_Unsaturation,diversity_func_Unsaturation)
    if(any(colnames(FtMSDataList[[i]]) %in% "lambda")){
      lambda_trait <- FtMSDataList[[i]] %>% select('Measured m/z',lambda) %>% column_to_rownames(var = 'Measured m/z')
      diversity_func_ThermoDyn <- as.data.frame(as.numeric(rao.diversity(mat,lambda_trait)$FunRao))
      colnames(diversity_func_ThermoDyn)[1] <- "diversity_Func_ThermoDyn"
      data_diversity_func_ThermoDyn <- rbind(data_diversity_func_ThermoDyn,diversity_func_ThermoDyn)
    }

  }
  if(any(colnames(FtMSDataList[[i]]) %in% "lambda")&(nrow(data_lambdaw)==length(FtMSDataList))){
    res <- cbind(data_DBEw,`data_DBE-Ow`,data_MWw,`data_O/Cw`,`data_H/Cw`,`data_N/Cw`,`data_S/Cw`,`data_P/Cw`,`data_DBE/Cw`,`data_(DBE-O)/Cw`,`data_(DBE-0.5O)/Cw`,data_AImodw,data_NOSCw,data_DeltaGw,data_lambdaw,data_MLB,data_diversity_shannon,diversity_simpson,data_diversity_func_NOSC,data_diversity_func_Unsaturation,data_diversity_func_ThermoDyn)
  }else{
    res <- cbind(data_DBEw,`data_DBE-Ow`,data_MWw,`data_O/Cw`,`data_H/Cw`,`data_N/Cw`,`data_S/Cw`,`data_P/Cw`,`data_DBE/Cw`,`data_(DBE-O)/Cw`,`data_(DBE-0.5O)/Cw`,data_AImodw,data_NOSCw,data_DeltaGw,data_MLB,data_diversity_shannon,diversity_simpson,data_diversity_func_NOSC,data_diversity_func_Unsaturation)
  }
  
  res <- round(res,digits = 3)
  rownames(res) <- names(FtMSDataList)
  return(res)
}
