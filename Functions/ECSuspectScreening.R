
FormulaSplit <- function(formula_ls){
  message("Transforing Formula to Element:\n" )
  if (ncol(formula_ls)>1){
    colnames(formula_ls) <- c("Reaction","Formula")
    Formula_df <- tibble(Formula = NA, C = NA, H=NA, N=NA, O=NA, P=NA, S=NA, Cl=NA, Br = NA, I=NA, F=NA)
    for (i in 1:nrow(formula_ls)) {
      formula <- formula_ls$Formula[i]
      formula_split <- str_split(formula, pattern = "")
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

CleanDatabaseData <- function(Dataset){
  if(str_detect(str_flatten(colnames(Dataset)),"Formula")){
    colnames(Dataset)[which(str_detect(colnames(Dataset),'Formula'))]<-'Formula'
  }
  if(str_detect(str_flatten(colnames(Dataset)),"Monoiso_Mass")){
    DataFilter <- filter(Dataset,!str_detect(Formula,'Si')) %>%
      filter(str_detect(Formula,'O')) %>%
      filter(!str_detect(Formula,'R')) %>%
      filter(!str_detect(Formula,'D')) %>%
      filter(!str_detect(Formula,'\\[')) %>%
      filter(!str_detect(Formula,'X')) %>%
      filter(!str_detect(Formula,'Se')) %>%
      filter(!str_detect(Formula,'Sn')) %>%
      filter(!str_detect(Formula,'n')) %>%
      filter(!str_detect(Formula,'T')) %>%
      filter(!str_detect(Formula,"(CF2)"))%>%
      filter(!str_detect(Formula,'Sb')) %>%
      filter(!str_detect(Formula,'Pb')) %>%
      filter(!str_detect(Formula,'Pt')) %>%
      filter(!str_detect(Formula,'A')) %>%
      filter(!str_detect(Formula,'Cs')) %>%
      filter(!str_detect(Formula,'Cu')) %>%
      filter(!str_detect(Formula,'d')) %>%
      filter(!str_detect(Formula,'e')) %>%
      filter(!str_detect(Formula,'C39H68N7O17N3S')) %>%
      filter(!str_detect(Formula,'g')) %>%
      filter(!str_detect(Formula,'a')) %>%
      filter(!str_detect(Formula,'K')) %>%
      filter(!str_detect(Formula,'Y')) %>%
      filter(!str_detect(Formula,'Bi')) %>%
      filter(!str_detect(Formula,'Ni')) %>%
      filter(!str_detect(Formula,'Li')) %>%
      filter(!str_detect(Formula,'L')) %>%
      filter(!str_detect(Formula,'E')) %>%
      filter(!str_detect(Formula,'Sr')) %>%
      filter(!str_detect(Formula,'Sm')) %>%
      filter(!str_detect(Formula,'SX')) %>%
      filter(!str_detect(Formula,'Sg')) %>%
      filter(!str_detect(Formula,'Sc')) %>%
      filter(!str_detect(Formula,'Os')) %>%
      filter(!str_detect(Formula,'Hf')) %>%
      filter(str_detect(Formula,'C')) %>%
      filter(!(str_detect(Formula,'B')&!str_detect(Formula,'r'))) %>%
      filter(!is.na(Monoiso_Mass)) %>%
      filter(!is.na(Formula))
  }else{
    DataFilter <- filter(Dataset,!str_detect(Formula,'Si')) %>%
      filter(str_detect(Formula,'O')) %>%
      filter(!str_detect(Formula,'R')) %>%
      filter(!str_detect(Formula,'D')) %>%
      filter(!str_detect(Formula,'X')) %>%
      filter(!str_detect(Formula,'\\[')) %>%
      filter(!str_detect(Formula,'Se')) %>%
      filter(!str_detect(Formula,'Sn')) %>%
      filter(!str_detect(Formula,'n')) %>%
      filter(!str_detect(Formula,'T')) %>%
      filter(!str_detect(Formula,"(CF2)"))%>%
      filter(!str_detect(Formula,'Sb')) %>%
      filter(!str_detect(Formula,'Pb')) %>%
      filter(!str_detect(Formula,'Pt')) %>%
      filter(!str_detect(Formula,'A')) %>%
      filter(!str_detect(Formula,'Cs')) %>%
      filter(!str_detect(Formula,'Cu')) %>%
      filter(!str_detect(Formula,'d')) %>%
      filter(!str_detect(Formula,'e')) %>%
      filter(!str_detect(Formula,'C39H68N7O17N3S')) %>%
      filter(!str_detect(Formula,'g')) %>%
      filter(!str_detect(Formula,'a')) %>%
      filter(!str_detect(Formula,'K')) %>%
      filter(!str_detect(Formula,'Y')) %>%
      filter(!str_detect(Formula,'Bi')) %>%
      filter(!str_detect(Formula,'Ni')) %>%
      filter(!str_detect(Formula,'Li')) %>%
      filter(!str_detect(Formula,'L')) %>%
      filter(!str_detect(Formula,'E')) %>%
      filter(!str_detect(Formula,'Sr')) %>%
      filter(!str_detect(Formula,'Sg')) %>%
      filter(!str_detect(Formula,'Sm')) %>%
      filter(!str_detect(Formula,'SX')) %>%
      filter(!str_detect(Formula,'Os')) %>%
      filter(!str_detect(Formula,'Hf')) %>%
      filter(!str_detect(Formula,'Sc')) %>%
      filter(str_detect(Formula,'C')) %>%
      filter(!(str_detect(Formula,'B')&!str_detect(Formula,'r'))) %>%
      filter(!is.na(Formula))
  }
  
  
  FormTran <- FormulaSplit(DataFilter['Formula'])
  DataForm <- left_join(DataFilter,FormTran,by='Formula')
  DataForm <- unique(DataForm)
  DataTable <- setDT(DataForm)
  DataTable <- DataTable[,mass:=C*12.0 + 
                           H * 1.00782503223 + 
                           N * 14.00307400443 + 
                           O * 15.99491461957 + 
                           P * 30.97376199842 +
                           S * 31.9720711744 +
                           Cl * 34.968852682 +
                           F * 18.99840316273 +
                           Br * 78.9183376 +
                           I * 126.9044719]
  DataTable <- DataTable[!is.na(C),]
  setkey(DataTable,'mass')
  return(DataTable)
}
  
  

databasetoTD <- function(databaseDF){
  ultramassmf <- setDT(databaseDF)
  setkey(ultramassmf, "mass")
  libversion <- 1
  ultramassmf[, vkey:=(libversion*10^12+1):(libversion*10^12+length(ultramassmf$mass))]
  return(ultramassmf)
}

prepare_peaklist <- function(FtMsDataLs, status.file_id = "available",status.peak_id = "not available",pol='neg',ion="H",ppm_dev=0.5){
  peaklist <- FtMsDataLs[[1]]
  peaklist <- peaklist[c("Measured m/z","Intensity","RA","S/N")]
  colnames(peaklist) <- c('mz','i_magnitude','RA',"s_n")
  
  peaklist <- setDT(peaklist)
  peaklist <- setkey(peaklist,'mz')
  if(status.file_id == "not available"){
    peaklist[, file_id:=rep("Sample_01",dim(peaklist)[1])]
  } else {
    peaklist[, file_id:=rep(names(FtMsDataLs)[1],dim(peaklist)[1])]
  }
  
  setkeyv(peaklist,c("file_id","mz"))
  
  peaklist <- data.table(peaklist)
  if(status.peak_id == "not available"){
    peaklist[, peak_id:=(1:(dim(peaklist)[1]))]
  }
  
  if(length(unique(peaklist$peak_id))!=dim(peaklist)[1]){
    peaklist[,peak_id:=(1:(dim(peaklist)[1]))]
  }
  
  peaklist <- f_um_calc_neutral_masses(peaklist,pol, ion=ion, ppm_dev)
  return(peaklist)
}




f_um_calc_neutral_masses <- function(peaklist, pol, ion="H", ppm_dev){
  
  if (pol=="neg")
  {
    peaklist[, m:=mz+1.0072763]
    peaklist[, m_min:=mz + 1.0072763-((mz + 1.0072763)*ppm_dev/1000000)]
    peaklist[, m_max := mz + 1.0072763 + ((mz + 1.0072763)*ppm_dev/1000000)]
  }
  
  if (pol=="pos") 
  {
    if (ion=="H"){
      peaklist[, m:=mz-1.0072763]
      peaklist[, m_min:= mz-1.0072763-((mz-1.0072763)*ppm_dev/1000000)]
      peaklist[, m_max:=mz-1.0072763 + ((mz-1.0072763)*ppm_dev/1000000)]
    }else if(ion=="Na"){
      peaklist[, m:=mz-22.9892215]
      peaklist[, m_min:= mz-22.9892215-((mz-22.9892215)*ppm_dev/1000000)]
      peaklist[, m_max:=mz-22.9892215 + ((mz-22.9892215)*ppm_dev/1000000)]
    }
  }
  
  if (pol=="neu") 
  {
    peaklist[, m:=mz]
    peaklist[, m_min:= mz-(mz*ppm_dev/1000000)]
    peaklist[, m_max:=mz+(mz*ppm_dev/1000000)]
  }
  return(peaklist)
}

# load public database

loadDatabase <- function(DatabaseDir){
  ultramassmf <- setDT(read.csv(DatabaseDir,header = TRUE))
  setkey(ultramassmf, "mass")
  libversion <- 1
  ultramassmf[, vkey:=(libversion*10^12+1):(libversion*10^12+length(ultramassmf$mass))]
  return(ultramassmf)
}


f_formula_query <- function(peaklist,
                            ultramassmf,
                            lib_version = 10 ^ 12,
                            status.s_n ="available",
                            status.peak_id = "available",
                            status.file_id= "available",
                            s_nTh = 6){
  print("***************************************************")
  print("Assigning molecular formulas to masses...")
  # make sure that peak magnitude is numeric
  peaklist$i_magnitude<-as.numeric(peaklist$i_magnitude)  
  
  # sort peaklist by mass
  setkey(peaklist,m)
  # assign peaklist masses
  s <- peaklist$m
  n <- length(s)
  dev <- (peaklist$m_max-peaklist$m_min)/2
  # assign library masses
  d <- ultramassmf$mass
  # ranges
  r.d <- range(d)
  r.s <- range(s)
  if (r.s[2] > r.d[2]){
    peaklist <- peaklist[m<=r.d[2],]
    # assign peaklist masses
    s <- peaklist$m
    n <- length(s)
    # ranges
    r.s <- range(s)
    
  }
  
  # check if sorted
  d.d <- diff(d)
  all(d.d >= 0)
  s.d <- diff(s)
  all(s.d >= 0)
  
  # insert s into d by index in d
  # check if range of samples (s) is within range of db (d)
  stopifnot (r.s[1] > r.d[1] && r.s[2] < r.d[2])
  # create vectors for indexing lower (pl) and upper limit (pt)
  pl = rep(NA,n)
  pt = rep(NA,n)
  
  # define start/end index
  #Kl = which.min(abs(d - dev[1]))
  Kl = which.min(abs(d-(s[1]-dev[1])))
  if (d[Kl]>s[1]-dev[1]){
    Kl=Kl-1
  }
  # Kt = which.min(abs(d + dev[1]))
  Kt = which.min(abs(d+(s[1]-dev[1])))
  if (d[Kt]>s[1] + dev[1]){
    Kt=Kt-1
  }
  
  # find elements in database whith lowest deviation
  
  for (I in 1:n){
    S=s[I]
    while(S - dev[I] > d[Kl]) {Kl = Kl + 1}
    while(S + dev[I] > d[Kt]) {Kt = Kt + 1}
    if (S + dev[I] > d[Kl] & S - dev[I] < d[Kt - 1]) {
      pl[I] = Kl;pt[I] = Kt - 1
    }
  }
  
  # create datatable of peaks and corresponding indexes
  #if(status.s_n!="not available"){matches <- data.table(peaklist[,.(peak_id,file_id,i_magnitude,mz,s_n)], s, pl, pt)}
  #if(status.s_n=="not available"){matches <- data.table(peaklist[,.(peak_id,file_id,i_magnitude,mz)], s, pl, pt)}
  # create datatable of peaks and corresponding indexes
  if(status.s_n!="not available"){matches <- data.table(peaklist[,.(peak_id,file_id,i_magnitude,mz,s_n,RA)], s, pl, pt)}
  if(status.s_n=="not available"){matches <- data.table(peaklist[,.(peak_id,file_id,i_magnitude,mz,RA)], s, pl, pt)}
  
  # remove non matched peaks
  matches <- na.omit(matches)
  if (nrow(matches)==0){
    return(NA)
  } else {
    # create list for indeces
    L <- as.list(1:length(matches$s))
    # fill list with indeces of respective ranges
    P.range <- lapply(L,function(x) matches$pl[x]:matches$pt[x])
    # create list of matched keys of db
    P.Keys <- sapply(P.range,function(x) as.list(ultramassmf$vkey[x]));
    
    if(status.s_n!="not available"){
      # unlist matches into new datatable
      matchesunl <- data.table(peak_id = unlist(lapply(L, function (x) rep(matches$peak_id[x], sapply(P.Keys[x],length)))),
                               file_id = unlist(lapply(L, function (x) rep(matches$file_id[x], sapply(P.Keys[x],length)))),
                               i_magnitude = unlist(lapply(L, function (x) rep(matches$i_magnitude[x], sapply(P.Keys[x],length)))),
                               mz = unlist(lapply(L, function (x) rep(matches$mz[x], sapply(P.Keys[x],length)))),
                               s_n = unlist(lapply(L, function (x) rep(matches$s_n[x], sapply(P.Keys[x],length)))),
                               RA = unlist(lapply(L, function (x) rep(matches$RA[x], sapply(P.Keys[x],length)))),
                               m = unlist(lapply(L, function (x) rep(matches$s[x], sapply(P.Keys[x],length))))
      )
    }
    
    if(status.s_n=="not available"){
      # unlist matches into new datatable
      matchesunl <- data.table(peak_id = unlist(lapply(L, function (x) rep(matches$peak_id[x], sapply(P.Keys[x],length)))),
                               file_id = unlist(lapply(L, function (x) rep(matches$file_id[x], sapply(P.Keys[x],length)))),
                               i_magnitude = unlist(lapply(L, function (x) rep(matches$i_magnitude[x], sapply(P.Keys[x],length)))),
                               mz = unlist(lapply(L, function (x) rep(matches$mz[x], sapply(P.Keys[x],length)))),
                               RA = unlist(lapply(L, function (x) rep(matches$RA[x], sapply(P.Keys[x],length)))),
                               m = unlist(lapply(L, function (x) rep(matches$s[x], sapply(P.Keys[x],length))))
      )
      
    }
    # Bind matches with corresponding rows in ultramass indexed by vkey
    df <- data.table(matchesunl, ultramassmf[(unlist(P.Keys)-lib_version)])
    df <- unique(df)
    df <- df[s_n>=s_nTh,]
    return(df)
  }
}

MainStepScreen <- function(FtMsDataLs,ultramassmf,pol='neg',ion='H',ppm_dev=0.5,s_nTh=6) {
  ResList <- list()
  for (i in 1:length(FtMsDataLs)){
    #prepare_peaklist data.table
    datatable <- prepare_peaklist(FtMsDataLs[i],pol=pol,ion = ion,ppm_dev=ppm_dev)
    #formula_query
    df <- f_formula_query(datatable,
                          ultramassmf,
                          lib_version = 10 ^ 12,
                          status.s_n ="available",
                          status.peak_id = "available",
                          status.file_id= "available",
                          s_nTh = s_nTh)
    if (all(is.na(df))){
      next
    }else{
      colnames(df)[3]<-'Intensity'
      df2<-unique(df[,!c('Name','Smiles','vkey','Source')])
      df2$Formula <- paste('C',df2$C,'H',df2$H,'Cl',df2$Cl,'Br',df2$Br,'I',df2$I,'F',df2$F,'N',df2$N,'O',df2$O,'S',df2$S,'P',df2$P,sep = '')%>%str_remove_all('N0')%>% str_remove_all('P0') %>% str_remove_all('S0')%>% str_remove_all('Cl0')%>% str_remove_all('Br0')%>% str_remove_all('I0')%>% str_remove_all('F0')%>% str_remove_all('O0')
      df2 <- unique(df2)
      #RawDf <-setDT(FtMsDataLs[[i]])[!`m/z` %in% df2$mz,]
      ResList[[i]]<- df2
      names(ResList)[i]<-names(FtMsDataLs)[i]
    }
  }
  return(ResList)
}




