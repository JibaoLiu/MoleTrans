
# 新的方法，直接从kegg 的 reaction pair 获取反应对信息，用于计算PMD
GetkeggReact_PMD <- function(PathID){
  message("Getting KEGG Reaction online:\n" )
  path_info <- keggGet(PathID)
  compound_info <- path_info[[1]]$COMPOUND
  cpd_id <- names(compound_info)
  i1=0
  Reaction_df <- tibble(Reactant=NA, Product=NA, Reaction=NA, Name=NA, Enzyme=NA, Orthology=NA, Pathway=NA)
  for (id in cpd_id){
    
    id_info <- keggGet(id)
    Reaction_ls <- id_info[[1]]$REACTION
    for (Rid in Reaction_ls){
      Rid_split <- str_split(Rid,pattern = " ")
      for (Rid_sub in Rid_split[[1]]){
        
        Sys.sleep(0.5)
        
        R_info <- keggGet(Rid_sub)
        if (str_detect(str_flatten(names(R_info[[1]])),"PATHWAY")|str_detect(str_flatten(names(R_info[[1]])),"ORTHOLOGY")){
          R_pathid <- names(R_info[[1]]$PATHWAY)
          R_Orthid <- names(R_info[[1]]$ORTHOLOGY)
          if (!is.null(R_Orthid)){
            Orth_info <- c()
            for (R_Orthidsup in R_Orthid) {
              Orth_infotmp <- keggGet(R_Orthidsup) #似乎存在问题
              Orth_infotmpname <- str_flatten(names(unlist(lapply(Orth_infotmp,'[[',"PATHWAY"))))
              Orth_info <- paste(Orth_info,Orth_infotmpname)
              i1=i1+2
              if(i1%%10==0){
                Sys.sleep(1.5)
              }
            }
            
          }else{
            Orth_info <- 'None'
          }
          
          if (str_detect(str_flatten(R_pathid),str_remove(PathID,"ko"))|str_detect(Orth_info,str_remove(PathID,"ko"))){
            if (str_detect(str_flatten(names(R_info[[1]])),"RCLASS")){
              RClass <- R_info[[1]]$RCLASS
              for (i in 1:length(RClass)) {
                RClassSplit <- unlist(str_split(RClass[i],pattern = " "))
                RPair <- RClassSplit[which(str_detect(RClassSplit,"_"))]
                RPairSplit <- unlist(str_split(RPair,pattern = '_'))
                if (!str_detect(str_flatten(names(R_info[[1]])),"NAME")){
                  if (str_detect(str_flatten(names(R_info[[1]])),"ORTHOLOGY")){
                    Reaction_df_tmp <- tibble(Reactant=RPairSplit[1],Product=RPairSplit[2],Reaction=R_info[[1]]$ENTRY,Name=NA, Enzyme=str_flatten(R_info[[1]]$ORTHOLOGY,collapse = ";"),Orthology=str_flatten(names(R_info[[1]]$ORTHOLOGY),collapse = ";"), Pathway=PathID)
                  }else{
                    if (str_detect(str_flatten(names(R_info[[1]])),"ENZYME")){
                      Reaction_df_tmp <- tibble(Reactant=RPairSplit[1],Product=RPairSplit[2],Reaction=R_info[[1]]$ENTRY,Name=NA, Enzyme=paste("EC",R_info[[1]]$ENZYME,sep = ":"),Orthology=NA, Pathway=PathID)
                    }else{
                      Reaction_df_tmp <- tibble(Reactant=RPairSplit[1],Product=RPairSplit[2],Reaction=R_info[[1]]$ENTRY,Name=NA, Enzyme=NA,Orthology=NA, Pathway=PathID)
                    }
                  }
                  
                } else {
                  if (str_detect(str_flatten(names(R_info[[1]])),"ORTHOLOGY")){
                    Reaction_df_tmp <- tibble(Reactant=RPairSplit[1],Product=RPairSplit[2],Reaction=R_info[[1]]$ENTRY,Name=R_info[[1]]$NAME, Enzyme=str_flatten(R_info[[1]]$ORTHOLOGY,collapse = ";"),Orthology=str_flatten(names(R_info[[1]]$ORTHOLOGY),collapse = ";"), Pathway=PathID)
                  }else{
                    if (str_detect(str_flatten(names(R_info[[1]])),"ENZYME")){
                      Reaction_df_tmp <- tibble(Reactant=RPairSplit[1],Product=RPairSplit[2],Reaction=R_info[[1]]$ENTRY,Name=R_info[[1]]$NAME, Enzyme=paste("EC",R_info[[1]]$ENZYME,sep = ":"), Orthology=NA, Pathway=PathID)
                    }else{
                      Reaction_df_tmp <- tibble(Reactant=RPairSplit[1],Product=RPairSplit[2],Reaction=R_info[[1]]$ENTRY,Name=R_info[[1]]$NAME, Enzyme=NA,Orthology=NA, Pathway=PathID)
                    }
                  }
                }
                Reaction_df <- rbind(Reaction_df,Reaction_df_tmp)
                Sys.sleep(0.2)
              }
            }
          }
        }
      }
    }
  }
  Reaction_df <-Reaction_df[2:nrow(Reaction_df),]
  Reaction_df <- unique(Reaction_df)
  return(Reaction_df)
}

GetkeggFormula_PMD <- function(Cid_df) {
  message("Getting KEGG Formula online:\n" )
  colnames(Cid_df)<-"Cid"
  i1=0
  Formula_df <- tibble(Cid=NA,Formula=NA)
  for (i in 1:nrow(Cid_df)) {
    cid <- Cid_df$Cid[i]
    c_info <- keggGet(cid)
    if (!str_detect(str_flatten(names(c_info[[1]])),"FORMULA")){
      Formula_df_tmp <- tibble(Cid=cid,Formula=NA)
    } else {
      Formula <- c_info[[1]]$FORMULA
      Formula_df_tmp <- tibble(Cid=cid,Formula=Formula)
    }
    Formula_df <- rbind(Formula_df,Formula_df_tmp)
    i1=i1+2
    if(i1%%10==0){
      Sys.sleep(1.5)
    }
    Sys.sleep(0.5)
  }
  Formula_df <- Formula_df[2:nrow(Formula_df),]
  return(Formula_df)
}

FormuTrans_keggpmd <- function(formula_ls){
  message("Transforing Formula to Element:\n" )
  if (ncol(formula_ls)>1){
    colnames(formula_ls) <- c("Reaction","Formula")
    Formula_df <- tibble(Reaction = NA, Formula = NA, C = NA, H=NA, N=NA, O=NA, P=NA, S=NA)
    for (i in 1:nrow(formula_ls)) {
      formula <- formula_ls$Formula[i]
      formula_split <- str_split(formula, pattern = "")
      split_num <- as.numeric(formula_split[[1]])
      split_num1<-split_num
      split_num[is.na(split_num)] <- 0
      if (!(str_detect(formula,"l")|str_detect(formula,"n")|str_detect(formula,"o")|str_detect(formula,"a")|str_detect(formula,"\\(")|str_detect(formula,"\\.")|is.na(formula))){
        if (str_detect(formula,"C")){
          C_index <- which(formula_split[[1]]%in%"C")
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
        Formula_temp <- tibble(Reaction = formula_ls$Reaction[i], Formula = formula_ls$Formula[i], C = C_num, H=H_num, N=N_num, O=O_num, P=P_num, S=S_num)
      } else {
        Formula_temp <- tibble(Reaction = formula_ls$Reaction[i], Formula = formula_ls$Formula[i], C = NA, H=NA, N=NA, O=NA, P=NA, S=NA)
      }
      Formula_df <- rbind(Formula_df,Formula_temp)
    }
    Formula_df <- Formula_df[2:nrow(Formula_df),]
  } else {
    colnames(formula_ls) <- "Formula"
    Formula_df <- tibble(Formula = NA, C = NA, H=NA, N=NA, O=NA, P=NA, S=NA)
    for (i in 1:nrow(formula_ls)) {
      formula <- formula_ls$Formula[i]
      formula_split <- str_split(formula, pattern = "")
      split_num <- as.numeric(formula_split[[1]])
      split_num1<-split_num
      split_num[is.na(split_num)] <- 0
      if (!(str_detect(formula,"l")|str_detect(formula,"n")|str_detect(formula,"o")|str_detect(formula,"a")|str_detect(formula,"\\(")|str_detect(formula,"\\.")|is.na(formula))){
        if (str_detect(formula,"C")){
          C_index <- which(formula_split[[1]]%in%"C")
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
        Formula_temp <- tibble(Formula = formula, C = C_num, H=H_num, N=N_num, O=O_num, P=P_num, S=S_num)
      } else {
        Formula_temp <- tibble(Reaction = formula_ls$Reaction[i], Formula = formula_ls$Formula[i], C = NA, H=NA, N=NA, O=NA, P=NA, S=NA)
      }
      Formula_df <- rbind(Formula_df,Formula_temp)
    }
    Formula_df <- Formula_df[2:nrow(Formula_df),]
  }
  return(Formula_df)
}

CalCAR_kegg <- function (PMD_df) {
  message("Caculating CAR:\n" )
  PMD_df$Nc <- apply(PMD_df[,c("C.x","C.y")],1,min)+apply(PMD_df[,c("N.x","N.y")],1,min)+apply(PMD_df[,c("O.x","O.y")],1,min)+apply(PMD_df[,c("P.x","P.y")],1,min)+apply(PMD_df[,c("S.x","S.y")],1,min)
  PMD_df$Nr <- PMD_df$C.x+PMD_df$N.x+PMD_df$O.x+PMD_df$P.x+PMD_df$S.x
  PMD_df$Np <- PMD_df$C.y+PMD_df$N.y+PMD_df$O.y+PMD_df$P.y+PMD_df$S.y
  PMD_df$CAR <- ((PMD_df$Nc/PMD_df$Nr+PMD_df$Nc/PMD_df$Np)/2)*(1-abs(PMD_df$Nc/PMD_df$Nr-PMD_df$Nc/PMD_df$Np))
  return(PMD_df)
}

CalPMD_kegg <- function(PMD_df) {
  message("Caculating PMD:\n" )
  new_df <- PMD_df
  new_df$Mass.x <- new_df$C.x*12+new_df$H.x*1.007825+new_df$O.x*15.994915+new_df$N.x*14.003074+new_df$S.x*31.972071+new_df$P.x*30.973762
  new_df$Mass.y <- new_df$C.y*12+new_df$H.y*1.007825+new_df$O.y*15.994915+new_df$N.y*14.003074+new_df$S.y*31.972071+new_df$P.y*30.973762
  new_df$PMD <- round(abs(new_df$Mass.x-new_df$Mass.y),digits = 6)
  new_df$C.p <- new_df$C.y-new_df$C.x
  new_df$H.p <- new_df$H.y-new_df$H.x
  new_df$O.p <- new_df$O.y-new_df$O.x
  new_df$N.p <- new_df$N.y-new_df$N.x
  new_df$P.p <- new_df$P.y-new_df$P.x
  new_df$S.p <- new_df$S.y-new_df$S.x
  #按照CAR筛选PMD
  PMD_df <- unique(filter(new_df, new_df$CAR>0.45&new_df$CAR<1&!(new_df$C.p==0&abs(new_df$H.p)==1&new_df$O.p==0&new_df$N.p==0&new_df$P.p==0&new_df$S.p==0)))
  #PMD只删除H2，H+的情况
  # PMD_df <- filter(new_df, new_df$CAR>0.2 &new_df$CAR<1&!(new_df$C.p==0&abs(new_df$H.p)==1&new_df$O.p==0&new_df$N.p==0&new_df$P.p==0&new_df$S.p==0))
  # 可根据CAR取值的大小，得到是否完全还原所有代谢网络
  PMD_unique <- unique(round(PMD_df["PMD"],digits = 6))
  PMD_info <- unique(PMD_df[c("PMD","C.p","H.p","O.p","N.p","P.p","S.p","Reaction.x")])
  colnames(PMD_info)<-str_remove_all(colnames(PMD_info),'.p')
  colnames(PMD_info)[ncol(PMD_info)]<-"Reaction"
  PMD_info$Group <- paste('C',PMD_info$C,'H',PMD_info$H,'O',PMD_info$O,'N',PMD_info$N,'P',PMD_info$P,'S',PMD_info$S,sep = '') %>% str_remove_all('C0')%>% str_remove_all('H0') %>% str_remove_all('O0') %>% str_remove_all('N0')%>% str_remove_all('P0') %>% str_remove_all('S0')
  PMD_info <-mutate(PMD_info,sum_index=C+H+O+N+P+S)
  PMD_info <-setDT(PMD_info)
  PMD_info <-setorder(PMD_info,-sum_index)
  undupid <- which(!duplicated(PMD_info$PMD))
  PMD_infotmp <- PMD_info[1,]
  for (i in undupid) {
    PMD_infotmptmp <- PMD_info[PMD_info$PMD %in% PMD_info$PMD[i],]
    PMD_infotmptmptmp <- PMD_infotmptmp[1,]
    PMD_infotmptmptmp$Reaction <- str_flatten(PMD_infotmptmp$Reaction,collapse = ";")
    PMD_infotmp <- rbind(PMD_infotmp,PMD_infotmptmptmp)
  }
  
  PMD_info <-PMD_infotmp[2:nrow(PMD_infotmp),-10]
  PMD_info <-setDF(PMD_info)
  PMD_list <- list(PMD_unique,PMD_info,PMD_df)
  names(PMD_list) <- c("PMD_unique","PMD_info","PMD_Orig")
  return(PMD_list)
}

KEGG_PMD <- function(PathID){
  React_df <- GetkeggReact_PMD(PathID)
  Formula1 <- GetkeggFormula_PMD(React_df["Reactant"])
  Formula2 <- GetkeggFormula_PMD(React_df["Product"])
  React_df <- cbind(React_df,Formula1["Formula"])
  colnames(React_df)[ncol(React_df)]<-"Formula1"
  React_df <- cbind(React_df,Formula2["Formula"])
  colnames(React_df)[ncol(React_df)]<-"Formula2"
  Reactant_df <- FormuTrans_keggpmd(React_df[,c("Reaction","Formula1")])
  Product_df <- FormuTrans_keggpmd(React_df[,c("Reaction","Formula2")])
  Reactant_df$rowname <- rownames(Reactant_df)
  Product_df$rowname <- rownames(Product_df)
  new_df <- full_join(Reactant_df,Product_df,by="rowname")
  new_df <- CalCAR_kegg(new_df)
  PMD_df <- CalPMD_kegg(new_df)
  PMD_df$PMD_info$PathID <- PathID
  Reaction_info <- React_df
  Reaction_info$label <- paste(Reaction_info$Formula1,Reaction_info$Formula2,sep = '-')
  PMD_info <- PMD_df$PMD_Orig[c("Reaction.x","Formula.x","Formula.y","PMD","C.p","H.p","O.p","N.p","P.p","S.p")]
  PMD_info$label <- paste(PMD_info$Formula.x,PMD_info$Formula.y,sep = '-')
  PMD_KO <- left_join(PMD_info,Reaction_info,by="label")
  PMD_KO <- PMD_KO[c("Reaction","Name","PMD","Formula1","Formula2","Enzyme","Orthology","Pathway","C.p","H.p","O.p","N.p","P.p","S.p")]
  PMD_KO <- unique(PMD_KO)
  PMD_ls <- list(React_df,new_df,PMD_df,PMD_KO)
  names(PMD_ls)<-c("Reaction_df","new_df","PMD_df","PMD_KO")
  return(PMD_ls)
}

KEGG_PMD_LS <- function(PathIDInput){
  PathIDVect <- str_split(PathIDInput,pattern = ';')[[1]]
  PMDresLs <- list()
  PathTotal <- tibble(PMD=NA,C=NA,H=NA,O=NA,N=NA,P=NA,S=NA,Reaction=NA,Group=NA,PathID=NA)
  for (PathID in PathIDVect) {
    PMDrestmp <- KEGG_PMD(PathID)
    PMDresLs[[PathID]]<-PMDrestmp$PMD_df$PMD_info
    PathTotal <- rbind(PathTotal,PMDrestmp$PMD_df$PMD_info)
  }
  PathTotal <-unique(PathTotal[2:nrow(PathTotal),])
  PMDresLs[["PathTotal"]]<-PathTotal
  return(PMDresLs)
}

