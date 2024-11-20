create_piechart<- function(data){
  rate1=sapply(unique(data$Group),function(x){sum(data$Group==x)})
  rate2= as.data.frame(rate1)
  sample_data <- data.frame(
    type = rownames(rate2),
    n = rate2[,1]
  )
  pg <- (sample_data %>% arrange(desc(n)) %>%
      mutate(ypos = cumsum(n) - n / 2) %>%
      mutate(per = n / sum(n)) %>%
      mutate(label = paste0(type, "\n", scales::percent(per, 0.1))) %>%
        ggpie("per",size=8/60,fill = 'type',label='label', lab.pos = 'out',lab.font= c(2, "plain", "black"), color = 'grey',alpha=0.6)
        #+ggrepel::geom_text_repel()
        +theme( plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm"))
  )
  pg <- ggpar(pg, legend = 'none')
  return(pg)
  
}



create_piechart2<- function(data,biochemicaltype="NewChemgroup"){
  if(biochemicaltype=="NewChemgroup"){
    rate1=sapply(unique(data$NewChemgroup),function(x){sum(data$NewChemgroup==x)})
  }else{
    rate1=sapply(unique(data$groupKVD),function(x){sum(data$groupKVD==x)})
  }
  rate2= as.data.frame(rate1)
  sample_data <- data.frame(
    type = rownames(rate2),
    n = rate2[,1]
  )
  pg <- (sample_data %>% arrange(desc(n)) %>%
           mutate(ypos = cumsum(n) - n / 2) %>%
           mutate(per = n / sum(n)) %>%
           mutate(label = paste0(type, "\n", scales::percent(per, 0.1))) %>%
           ggpie("per",size=8/60,label = 'label',fill = 'type',lab.pos = 'out',lab.font = c(2,"plain", 'black'),color = 'grey',alpha=0.6)
           #+fill_palette('npg')
           +theme( plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm"))
  )
  pg <- ggpar(pg, legend = 'none')
  return(pg) 
}

