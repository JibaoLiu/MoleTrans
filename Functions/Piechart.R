create_piechart<- function(data){
  rate1=sapply(unique(data$Group),function(x){sum(data$Group==x)})
  rate2= as.data.frame(rate1)
  sample_data <- data.frame(
    type = rownames(rate2),
    n = rate2[,1]
  )
  
  pie_data <- sample_data %>% arrange(desc(n)) %>%
    mutate(per = n *100 / sum(n))
    #mutate(ypos = cumsum(per) - per / 2) %>%
    #mutate(label = paste0(type, "\n", scales::percent(per, 0.1))) 
  
  positionDF <- pie_data |>
    arrange(desc(type)) |>
    mutate(
      csum = cumsum(per),
      pos = (csum + dplyr::lag(csum, default = 0)) / 2
    )

  
  baseplot <- pie_data |>
    ggplot(aes(x = "", y = per, fill = type)) +
    geom_col(color = "grey", width = .5) +
    geom_text_repel(
      data = positionDF,
      aes(
        x = 1.2,
        y = pos,
        label = paste(
          type,
          scales::number(per, accuracy = .1, suffix = " %")
        )
      ),
      size = 3,
      color = "black",
      fontface = "italic",
      nudge_x = 0.3,
      show.legend = FALSE
    )
  
  themedChart <- baseplot +
    #scale_fill_manual(values = color_palette) +
    #labs(caption = "Collisions by light type - 2020 to 2024") +
    theme_void() +
    theme(
      plot.caption = element_text(
        hjust = 0.5, vjust = 10,
        size = 10, face = "bold"
      ),
      plot.margin = unit(c(1, 1, 1, 1), "cm"),
      panel.background = element_rect(fill = "white"),
      legend.position = "none"
    )
  
  p <- themedChart +
        coord_polar("y", start = pi / 2)
  
  return(p)
  
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
  pie_data <- sample_data %>% arrange(desc(n)) %>%
    mutate(per = n *100 / sum(n))
  #mutate(ypos = cumsum(per) - per / 2) %>%
  #mutate(label = paste0(type, "\n", scales::percent(per, 0.1))) 
  
  positionDF <- pie_data |>
    arrange(desc(type)) |>
    mutate(
      csum = cumsum(per),
      pos = (csum + dplyr::lag(csum, default = 0)) / 2
    )
  
  
  baseplot <- pie_data |>
    ggplot(aes(x = "", y = per, fill = type)) +
    geom_col(color = "grey", width = .5) +
    geom_text_repel(
      data = positionDF,
      aes(
        x = 1.2,
        y = pos,
        label = paste(
          type,
          scales::number(per, accuracy = .1, suffix = " %")
        )
      ),
      size = 3,
      color = "black",
      fontface = "italic",
      nudge_x = 0.3,
      show.legend = FALSE
    )
  
  themedChart <- baseplot +
    #scale_fill_manual(values = color_palette) +
    #labs(caption = "Collisions by light type - 2020 to 2024") +
    theme_void() +
    theme(
      plot.caption = element_text(
        hjust = 0.5, vjust = 10,
        size = 10, face = "bold"
      ),
      plot.margin = unit(c(1, 1, 1, 1), "cm"),
      panel.background = element_rect(fill = "white"),
      legend.position = "none"
    )
  
  pg <- themedChart +
    coord_polar("y", start = pi / 2)

  return(pg) 
}

