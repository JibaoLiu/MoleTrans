
############################################
# Machine Learning Using tidymodels    #####
############################################
#两两样品比较，设定标签
CompareLable <- function(FtMSData1,FtMSData2){
  CompleteRemove <- FtMSData1[[1]][!FtMSData1[[1]]$Formula %in% FtMSData2[[1]]$Formula,]
  CompleteRemove$Label <- 'Disappeared'
  Produced <- FtMSData2[[1]][!FtMSData2[[1]]$Formula %in% FtMSData1[[1]]$Formula,]
  Produced$Label <- 'Newly_appeared'
  
  #不考虑C13
  Common <- FtMSData1[[1]][!duplicated(FtMSData1[[1]]$Formula),][FtMSData1[[1]][!duplicated(FtMSData1[[1]]$Formula),]$Formula %in% FtMSData2[[1]][!duplicated(FtMSData2[[1]]$Formula),]$Formula, ]
  Common$Int1 <- FtMSData1[[1]][!duplicated(FtMSData1[[1]]$Formula),][FtMSData1[[1]][!duplicated(FtMSData1[[1]]$Formula),]$Formula %in% Common$Formula,]$ReInt
  Common$Int2 <- FtMSData2[[1]][!duplicated(FtMSData2[[1]]$Formula),][FtMSData2[[1]][!duplicated(FtMSData2[[1]]$Formula),]$Formula %in% Common$Formula,]$ReInt
  Common$fc <- Common$Int2/Common$Int1
  PartiallyRemove <- filter(Common,fc < 0.5)[,1:(ncol(Common)-3)]#降低的是前驱物的相对峰强
  PartiallyRemove$Label <- "Partially Remove"
  Common2 <- FtMSData2[[1]][!duplicated(FtMSData2[[1]]$Formula),][FtMSData2[[1]][!duplicated(FtMSData2[[1]]$Formula),]$Formula %in% FtMSData1[[1]][!duplicated(FtMSData1[[1]]$Formula),]$Formula, ]
  Common2$Int1 <- FtMSData1[[1]][!duplicated(FtMSData1[[1]]$Formula),][FtMSData1[[1]][!duplicated(FtMSData1[[1]]$Formula),]$Formula %in% Common2$Formula,]$ReInt
  Common2$Int2 <- FtMSData2[[1]][!duplicated(FtMSData2[[1]]$Formula),][FtMSData2[[1]][!duplicated(FtMSData2[[1]]$Formula),]$Formula %in% Common2$Formula,]$ReInt
  Common2$fc <- Common2$Int2/Common2$Int1
  Increased <- filter(Common2,fc > 2)[,1:(ncol(Common2)-3)]#Reint是增加物的相对峰强
  Increased$Label <- "Increased"
  
  Retained <- filter(Common,fc>=0.5 & fc <= 2)[,1:(ncol(Common)-3)]
  Retained$Label <- "Retained"
  bDOM <- rbind(CompleteRemove,PartiallyRemove)
  bDOM$Label <- 'bDOM'
  rDOM <- rbind(Produced,Increased)%>% rbind(Retained)
  rDOM$Label <- 'rDOM'
  ResList <- list(CompleteRemove,PartiallyRemove,Retained,Increased,Produced,bDOM,rDOM)
  names(ResList)<-c('Completely Remove','Partially Remove','Retained','Increased','Produced','bDOM','rDOM')
  return(ResList)
}

DataPreforML <- function(FtMsDf1,FtMsDf2,balanceData=FALSE){
  CompareRes <- CompareLable(FtMsDf1,FtMsDf2)
  #合并为三个标签
  ResistantData <- rbind(CompareRes$Retained,CompareRes$Increased) %>% rbind(CompareRes$`Partially Remove`)
  ResistantData$Label <- 'Resistant'
  ResistantData$Label_Code <- 0
  CompareRes$`Completely Remove`$Label_Code <- 1
  CompareRes$Produced$Label_Code <- 2
  DataSet <- rbind(CompareRes$`Completely Remove`,CompareRes$Produced)%>%rbind(ResistantData)
  DataSet <- as.data.frame(DataSet)
  #DataSetNew <- DataSet[c('Measured m/z','O/C','H/C','N/C','S/C','P/C','AImod','NOSC','DBE-O','DBE/C','X/C','lambda','Label','Label_Code')]
  #colnames(DataSetNew)<-c('Measuredmz','OC','HC','NC','SC','PC','AImod','NOSC','DBEO','DBEC','X/C','lambda','Label','LabelCode')
  if(str_detect(str_flatten(colnames(DataSet)),"lambda")){
    DataSetNew <- DataSet[c('Measured m/z','O/C','H/C','N/C','S/C','P/C','AImod','NOSC','DBE-O','DBE/C','(DBE-0.5O)/C','lambda','Label','Label_Code')]
    colnames(DataSetNew)<-c('Measured_mz','O_C','H_C','N_C','S_C','P_C','AImod','NOSC','DBE_O','DBE_C','(DBE_0.5O)_C','lambda','Label','LabelCode')
  }else if(any(colnames(DataSet) %in% "P/C")){
    DataSetNew <- DataSet[c('Measured m/z','O/C','H/C','N/C','S/C','P/C','AImod','NOSC','DBE-O','DBE/C','(DBE-0.5O)/C','Label','Label_Code')]
    colnames(DataSetNew)<-c('Measured_mz','O_C','H_C','N_C','S_C','P_C','AImod','NOSC','DBE_O','DBE_C','(DBE_0.5O)_C','Label','LabelCode')
  }else{
    DataSetNew <- DataSet[c('Measured m/z','O/C','H/C','N/C','S/C','AImod','NOSC','DBE-O','DBE/C','(DBE-0.5O)/C','Label','Label_Code')]
    colnames(DataSetNew)<-c('Measured_mz','O_C','H_C','N_C','S_C','AImod','NOSC','DBE_O','DBE_C','(DBE_0.5O)_C','Label','LabelCode')
  }
  
  DataSetNew <-setDT(DataSetNew)
  setkey(DataSetNew,'Measured_mz')
  DataSetNew <-setDF(DataSetNew)
  if(balanceData==TRUE){
    DataSetNew$LabelCode <- as.factor(DataSetNew$LabelCode)
    trainup<-upSample(x=DataSetNew[,-ncol(DataSetNew)],
                      y=DataSetNew$LabelCode)
    trainup$Class <- as.numeric(trainup$Class)-1
    colnames(trainup)[ncol(trainup)] <- "LabelCode"
    x_train <- as.matrix(subset(trainup,select = -c(Label,LabelCode)))
    y_train <- as.matrix(trainup['LabelCode'])
    res <- list(x_train,y_train,trainup)
  }else{
    x_train <- as.matrix(subset(DataSetNew,select = -c(Label,LabelCode)))
    y_train <- as.matrix(DataSetNew['LabelCode'])
    res <- list(x_train,y_train,DataSetNew)
  }
  names(res)<-c("x_train","y_train","FormulaBehavior")
  return(res)
  #y_train <- to_categorical(y_train)
}

Venn <- function(FtMsDf1,FtMsDf2){
  CompareRes <- CompareLable(FtMsDf1,FtMsDf2)
  #合并为三个标签
  ResistantData <- rbind(CompareRes$Retained,CompareRes$Increased) %>% rbind(CompareRes$`Partially Remove`)
  ResistantData$Label <- 'Resistant'
  ResistantData$Label_Code <- 0
  CompareRes$`Completely Remove`$Label_Code <- 1
  CompareRes$Produced$Label_Code <- 2
  DataSet <- rbind(CompareRes$`Completely Remove`,CompareRes$Produced)%>%rbind(ResistantData)
  `Newly_appeared` <- filter(DataSet,Label=='Newly_appeared')$"Formula"
  `Resistant` <- filter(DataSet,Label=='Resistant')$"Formula"
  `Disappeared` <- filter(DataSet,Label=='Disappeared')$"Formula"
  a <- list(`Newly_appeared`=Newly_appeared,`Resistant`=Resistant,`Disappeared`=Disappeared)
  fig <- ggvenn(a,set_name_size = 3,text_size=3, stroke_size = 0.5)
  return(fig)
}

Venn_Merge <- function(FtMsDfLs){
  a <- list()
  for (i in 1:length(FtMsDfLs)) {
    a[[i]] <- FtMsDfLs[[i]]$Formula
    names(a)[i] <- names(FtMsDfLs)[i]
  }
  fig <- ggvenn(a,set_name_size =3,text_size = 3,stroke_size = 0.5)
  return(fig)
}

Sample_Merge <- function(FtMsDfLs){
  data_Merge <- FtMsDfLs[[1]]
  inner_formula <- data_Merge$Formula
  if(length(FtMsDfLs)>1){
    for (i in 2:length(FtMsDfLs)) {
      inner_formula <- intersect(inner_formula,FtMsDfLs[[i]]$Formula)
      data_Merge <- data_Merge[data_Merge$Formula %in% inner_formula,]
    }
  }else{
    data_Merge <- FtMsDfLs[[1]]
  }
  return(data_Merge)
}


ml_dataset <- function(x_train,y_train, vfold=10){
  data_set <- cbind(as.data.frame(x_train),y_train)
  data_set <- as_tibble(data_set)
  data_set$LabelCode <- as.factor(data_set$LabelCode)
  trees_split <- initial_split(data_set, strata = LabelCode)
  train_set <- training(trees_split)
  test_set <- testing(trees_split)
  set.seed(123)
  folds <- vfold_cv(train_set, strata = LabelCode, v = vfold)
  recipe <- recipe(LabelCode~.,data = train_set) 
  res <- list(recipe, train_set, test_set,folds,trees_split)
  names(res) <- c("recipe", "train_set", "test_set","folds","tree_split")
  return(res)
}





conf_matrix <- function(df.true, df.pred, title = "", true.lab ="True Class", pred.lab ="Predicted Class",
                        high.col = 'red', low.col = 'white') {
  #convert input vector to factors, and ensure they have the same levels
  df.true <- as.factor(df.true)
  df.pred <- factor(df.pred, levels = levels(df.true))
  
  #generate confusion matrix, and confusion matrix as a pecentage of each true class (to be used for color) 
  df.cm <- table(True = df.true, Pred = df.pred)
  df.cm.col <- df.cm / rowSums(df.cm)
  
  #convert confusion matrices to tables, and binding them together
  df.table <- reshape2::melt(df.cm)
  df.table.col <- reshape2::melt(df.cm.col)
  df.table <- left_join(df.table, df.table.col, by =c("True", "Pred"))
  
  #calculate accuracy and class accuracy
  acc.vector <- c(diag(df.cm)) / c(rowSums(df.cm))
  class.acc <- data.frame(Pred = "Class Acc.", True = names(acc.vector), value = acc.vector)
  acc <- sum(diag(df.cm)) / sum(df.cm)
  
  #plot
  ggplot() +
    geom_tile(aes(x=Pred, y=True, fill=value.y),
              data=df.table, size=0.2, color=grey(0.5)) +
    geom_tile(aes(x=Pred, y=True),
              data=df.table[df.table$True==df.table$Pred, ], size=1, color="black", fill = 'transparent') +
    #scale_x_discrete(position = "top",  limits = c(levels(df.table$Pred), "Class Acc.")) +
    #scale_y_discrete(limits = rev(unique(levels(df.table$Pred)))) +
    labs(x=pred.lab, y=true.lab, fill=NULL,
         title= paste0(title, "\nAccuracy ", round(100*acc, 1), "%")) +
    geom_text(aes(x=Pred, y=True, label=value.x),
              data=df.table, size=4, colour="black") +
    #geom_text(data = class.acc, aes(Pred, True, label = paste0(round(100*value), "%"))) +
    scale_fill_gradient(low=low.col, high=high.col, labels = scales::percent,
                        limits = c(0,1), breaks = c(0,0.5,1)) +
    guides(size=F) +
    theme_bw() +
    theme(panel.border = element_blank(), legend.position = "bottom",
          axis.text = element_text(color='black'), axis.ticks = element_blank(),
          panel.grid = element_blank()) +
    coord_fixed()
  
}

shapvalue_cal <- function(train_set, test_set,final_wf,
                          classification_type="multiclass"){
  X <- test_set[sample(nrow(test_set), min(nrow(test_set)/2,200)), -ncol(test_set)]
  bg_X <- train_set[sample(nrow(train_set), min(nrow(train_set)/2,100)), ]
  if(classification_type=="multiclass"){
    # multi classes
    pred_fun <- function(mod, X) predict(mod, X, type="prob")
  }else{
    pred_fun <- function(mod, X) {
      pre <- predict(mod, X)
      pre$.pred_class <- as.numeric(pre$.pred_class)-1
      pre
    }
  }
  fit <- final_wf %>%
    fit(data = train_set) 
  shap_rf <- kernelshap(fit, X, bg_X =bg_X,pred_fun = pred_fun)
  sv_rf <- shapviz(shap_rf)
  return(sv_rf)
}

shapley_plot <- function(sv_rf,type="type1",x_var, color_var){
  if(type=="type1"){
    p <-sv_importance(sv_rf, show_numbers = TRUE)
  }else if(type=="type2"){
    p <-sv_importance(sv_rf, show_numbers = TRUE, bar_type="stack")
  }else if(type=="type3"){
    p <-sv_importance(sv_rf, kind = "beeswarm",bee_width=0.4,max_display=18L,show_numbers = FALSE)
  }else if(type=="type4"){
    p <- sv_dependence(sv_rf, x_var, color_var = color_var)
  }
  return(p)
}

shapplot_cal_single <- function(train_set, DataSet,int,final_wf,
                          classification_type="multiclass"){
  #X <- test_set[sample(nrow(test_set), min(nrow(test_set)/2,200)), -ncol(test_set)]
  if(str_detect(str_flatten(colnames(DataSet)),"lambda")){
    test_set <- DataSet[c('Measured m/z','O/C','H/C','N/C','S/C','P/C','AImod','NOSC','DBE-O','DBE/C','(DBE-0.5O)/C','lambda')]
    colnames(test_set)<-c('Measured_mz','O_C','H_C','N_C','S_C','P_C','AImod','NOSC','DBE_O','DBE_C','(DBE_0.5O)_C','lambda')
  }else{
    test_set <- DataSet[c('Measured m/z','O/C','H/C','N/C','S/C','P/C','AImod','NOSC','DBE-O','DBE/C','(DBE-0.5O)/C')]
    colnames(test_set)<-c('Measured_mz','O_C','H_C','N_C','S_C','P_C','AImod','NOSC','DBE_O','DBE_C','(DBE_0.5O)_C')
  }
  test_set <- test_set[,colnames(train_set[,-ncol(train_set)])]
  single_row <- as.data.frame(test_set)[int,]
  bg_X <- train_set[sample(nrow(train_set), min(nrow(train_set)/2,100)), ]
  if(classification_type=="multiclass"){
    # multi classes
    pred_fun <- function(mod, X) predict(mod, X, type="prob")
  }else{
    pred_fun <- function(mod, X) {
      pre <- predict(mod, X)
      pre$.pred_class <- as.numeric(pre$.pred_class)-1
      pre
    }
  }
  res_plot <- final_wf %>%
    fit(data = train_set) |>
    kernelshap(single_row, bg_X = bg_X, pred_fun = pred_fun) |> 
    shapviz() |>
    sv_waterfall()
  return(res_plot)
}



