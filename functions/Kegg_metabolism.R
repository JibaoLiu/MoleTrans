
get.kegg.byId<-
function(keggId) {
  kegg = data.frame()
  i = 1
  while(i <= length(keggId)) {
    
    cat('processing', keggId[i], '\n')
    query <- keggGet(keggId[i:(i+9)])
    
    for(l in 1:length(query)) {
      
      keggRow = query[[l]]
      
      for(j in names(keggRow)) {
        if(j == 'DBLINKS') {
          for(k in 1:length(keggRow$DBLINKS)) {
            db = unlist(strsplit(keggRow$DBLINKS[k], ': '))[1]
            id = unlist(strsplit(keggRow$DBLINKS[k], ': '))[2]
            keggRow[[db]] = id
          }
        } else if (j == 'PATHWAY') {
          for(k in 1:length(keggRow$PATHWAY)) {
            keggRow$PATHWAY[k] = paste(names(keggRow$PATHWAY[k]), keggRow$PATHWAY[k], sep=': ')
          }
          keggRow$PATHWAY = paste(keggRow$PATHWAY, collapse='///')
        } else if (j == 'REFERENCE') {
          keggRow$REFERENCE = paste(keggRow$REFERENCE[[1]]$REFERENCE, collapse='///')
        } else {
          if(length(keggRow[[j]]) > 1) {
            keggRow[[j]] = paste(keggRow[[j]], collapse='///')
          }
        }
      }
      keggRow[['DBLINKS']] = NULL
      keggRow = as.data.frame(keggRow, stringsAsFactors=FALSE)
      kegg = rbind.fill(kegg, keggRow)
      kegg[is.na(kegg)] = ''
    }
    i = i + 10
    Sys.sleep(2)
  }
  return(kegg)
}


get.kegg.all <-
  function() {
    cmp <- keggList("compound")
    reactionEntry = keggList("reaction")
    
    cmpId = names(cmp)
    cmpId = sub('cpd:', '', cmpId)
    
    reactionEntry = names(reactionEntry)
    reactionEntry = sub('rn:', '', reactionEntry)
    
    keggReaction = get.kegg.byId(reactionEntry)
    keggReaction[is.na(keggReaction)] = ""
    
    keggCompound = get.kegg.byId(cmpId)
    keggCompound[is.na(keggCompound)] = ""
    
    # reference
    referIndex = grep('.+', keggReaction$REFERENCE)
    referId = keggReaction[grep('.+', keggReaction$REFERENCE), 'ENTRY']
    referIdUnique = unique(keggReaction[grep('.+', keggReaction$REFERENCE), 'ENTRY'])
    
    redundantIndex = c()
    for(i in referIdUnique) {
      index = grep(i, referId)
      index = referIndex[index[-1]]
      redundantIndex = c(redundantIndex, index)
    }
    
    if(length(redundantIndex) > 0) {
      keggReaction_unique = keggReaction[-redundantIndex,]
    } else {
      keggReaction_unique = keggReaction
    }
    
    result = list()
    result[['reaction']] = keggReaction_unique
    result[['compound']] = keggCompound
    cat('# of reactions:', nrow(keggReaction_unique), '\n')
    cat('# of compounds:', nrow(keggCompound), '\n')
    return(result)
  }

prepare_peaklist_kegg <- function(FtMsDataLs, status.file_id = "available",status.peak_id = "not available",pol='neg',ion="H",ppm_dev=0.5){
  peaklist <- FtMsDataLs[[1]]
  if(any(colnames(peaklist) %in% "S/N")){
    peaklist <- peaklist[,c("Measured m/z","Intensity","RA","S/N","Formula")]
    colnames(peaklist) <- c('mz','i_magnitude','RA',"s_n","Formula")
  }else{
    peaklist <- peaklist[,c("Measured m/z","Intensity","RA","Formula")]
    colnames(peaklist) <- c('mz','i_magnitude','RA',"Formula")
    
  }
  
  
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

Prepare_KEGGDB <- function(database){
  PreKeggDatabase <- setDT(database)
  if(any(colnames(PreKeggDatabase) %in% "EXACT_MASS")){
    PreKeggDatabase <- filter(PreKeggDatabase, !EXACT_MASS=="")
    colnames(PreKeggDatabase)[which(colnames(PreKeggDatabase)=="EXACT_MASS")] <- 'mass'
    
  }
  PreKeggDatabase$mass <- as.numeric(PreKeggDatabase$mass)
  setkey(PreKeggDatabase, "mass")
  libversion <- 1
  PreKeggDatabase[, vkey:=(libversion*10^12+1):(libversion*10^12+length(PreKeggDatabase$mass))]
  return(PreKeggDatabase)
}

f_formula_query_kegg <- function(peaklist,
                            PreKeggDatabase,
                            lib_version = 10 ^ 12,
                            status.s_n ="available",
                            status.peak_id = "available",
                            status.file_id= "available"){
  print("***************************************************")
  #print("Assigning molecular formulas to masses...")
  # make sure that peak magnitude is numeric
  peaklist$i_magnitude<-as.numeric(peaklist$i_magnitude)  
  
  # sort peaklist by mass
  setkey(peaklist,m)
  # assign peaklist masses
  s <- peaklist$m
  n <- length(s)
  dev <- (peaklist$m_max-peaklist$m_min)/2
  # assign library masses
  d <- PreKeggDatabase$mass
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
  if(status.s_n!="not available"){matches <- data.table(peaklist[,.(peak_id,file_id,i_magnitude,mz,Formula,s_n,RA)], s, pl, pt)}
  if(status.s_n=="not available"){matches <- data.table(peaklist[,.(peak_id,file_id,i_magnitude,mz,Formula,RA)], s, pl, pt)}
  
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
    P.Keys <- sapply(P.range,function(x) as.list(PreKeggDatabase$vkey[x]));
    
    if(status.s_n!="not available"){
      # unlist matches into new datatable
      matchesunl <- data.table(peak_id = unlist(lapply(L, function (x) rep(matches$peak_id[x], sapply(P.Keys[x],length)))),
                               file_id = unlist(lapply(L, function (x) rep(matches$file_id[x], sapply(P.Keys[x],length)))),
                               i_magnitude = unlist(lapply(L, function (x) rep(matches$i_magnitude[x], sapply(P.Keys[x],length)))),
                               mz = unlist(lapply(L, function (x) rep(matches$mz[x], sapply(P.Keys[x],length)))),
                               Formula = unlist(lapply(L, function (x) rep(matches$Formula[x], sapply(P.Keys[x],length)))),
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
                               Formula = unlist(lapply(L, function (x) rep(matches$Formula[x], sapply(P.Keys[x],length)))),
                               RA = unlist(lapply(L, function (x) rep(matches$RA[x], sapply(P.Keys[x],length)))),
                               m = unlist(lapply(L, function (x) rep(matches$s[x], sapply(P.Keys[x],length))))
      )
      
    }
    # Bind matches with corresponding rows in ultramass indexed by vkey
    df <- data.table(matchesunl,PreKeggDatabase[(unlist(P.Keys)-lib_version)])
    df <- unique(df)
    #df <- df[s_n>=s_nTh,]
    df <- df[,3:16]
    df$`MassError(ppm)` = round(abs(df$mass-df$m)*1000000/df$mass,digits = 4)
    colnames(df)[1] <- "Intensity"
    return(df)
  }
}

MainStepScreen_KEGG <- function(FtMsDataLs,KeggDatabase,pol='neg',ion='H',ppm_dev=0.5) {
  ResList <- list()
  PreKeggDatabase <- Prepare_KEGGDB(KeggDatabase)
  for (i in 1:length(FtMsDataLs)){
    #prepare_peaklist data.table
    datatable <- prepare_peaklist_kegg(FtMsDataLs[i],pol=pol,ion = ion,ppm_dev=ppm_dev)
    #formula_query
    df <- f_formula_query_kegg(datatable,
                          PreKeggDatabase,
                          lib_version = 10 ^ 12,
                          status.s_n ="available",
                          status.peak_id = "available",
                          status.file_id= "available")
    df<- as.data.frame(df)
    if (all(is.na(df))){
      next
    }else{
      #colnames(df)[3]<-'Intensity'
      #df2<-unique(df[,!c('Name','Smiles','vkey','Source')])
      #df2$Formula <- paste('C',df2$C,'H',df2$H,'Cl',df2$Cl,'Br',df2$Br,'I',df2$I,'F',df2$F,'N',df2$N,'O',df2$O,'S',df2$S,'P',df2$P,sep = '')%>%str_remove_all('N0')%>% str_remove_all('P0') %>% str_remove_all('S0')%>% str_remove_all('Cl0')%>% str_remove_all('Br0')%>% str_remove_all('I0')%>% str_remove_all('F0')%>% str_remove_all('O0')
      #df2 <- unique(df2)
      #RawDf <-setDT(FtMsDataLs[[i]])[!`m/z` %in% df2$mz,]
      
      ResList[[i]]<- df
      names(ResList)[i]<-names(FtMsDataLs)[i]
    }
  }
  return(ResList)
}




FindKeggMetabo_Form <- function(FtMsDf){
  message("Find Kegg Formula\n" )
  if(!str_detect(str_flatten(colnames(FtMsDf)),"Neutral m/z")){
    FtMsDf$`Neutral m/z` <- FtMsDf$`Measured m/z`+1.0072765
  }
  mslist <- FtMsDf['Neutral m/z']
  formulalist <- FtMsDf['Formula']
  colnames(mslist)<-"Mass"
  colnames(formulalist)<-"Formula"
  kegg_formula <- tibble(Cpd_id=NA, Mass=NA,KEGGformula=NA, OrigFormula=NA, Label=NA,`MassError(ppm)`=NA)
  for (i in 1:nrow(mslist)){
    message(i)
    mass <- trunc(mslist$Mass[i]*1000)/1000
    cpd_id <- names(keggFind('compound', mass, 'exact_mass'))
    if (!is.null(cpd_id)){
      cpd_info <- keggGet(cpd_id[1])
      if (str_detect(cpd_info[[1]]$FORMULA,"l")|str_detect(cpd_info[[1]]$FORMULA,"o")|str_detect(cpd_info[[1]]$FORMULA,"n")|str_detect(cpd_info[[1]]$FORMULA,"\\.")|str_detect(cpd_info[[1]]$FORMULA,"a")|str_detect(cpd_info[[1]]$FORMULA,"\\(")){
        temp <- tibble(Mass=round(mslist$Mass[i],digits=6), KEGGformula = formulalist$Formula[i], OrigFormula=formulalist$Formula[i], Label='Found_but_unconsistant')
      } else {
        temp <- tibble(Mass=round(mslist$Mass[i],digits=6), KEGGformula = cpd_info[[1]]$FORMULA, OrigFormula=formulalist$Formula[i], Label='Found')
      }
      cpd_id <- str_remove_all(cpd_id,"cpd:")
      res_tmp <- as.data.frame(cpd_id)
      colnames(res_tmp) <- "Cpd_id"
      res_tmp$Mass=temp$Mass
      res_tmp$KEGGformula=temp$KEGGformula
      res_tmp$OrigFormula=temp$OrigFormula
      res_tmp$Label = temp$Label
      formula_element=FormuTrans(temp[,'KEGGformula'])
      formula_element = formula_element%>%mutate(exactmass=C*12.0 + 
                                                   H * 1.00782503223 + 
                                                   N * 14.00307400443 + 
                                                   O * 15.99491461957 + 
                                                   P * 30.97376199842 +
                                                   S * 31.9720711744)
      res_tmp$`MassError(ppm)` = round(abs(res_tmp$Mass-formula_element$exactmass)*1000000/formula_element$exactmass,digits = 4)
      
      kegg_formula <- rbind(kegg_formula,res_tmp)
    }
    
    Sys.sleep(0.5)
  } 
  kegg_formula <- kegg_formula[2:nrow(kegg_formula),]
  return(kegg_formula)
}