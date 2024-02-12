#SASP Volcano

volcano.SASP <- function(x){
  require(ggplot2)
  require(stringr)
  df <- all_comparisons[all_comparisons$Comparison..group1.group2. == x,]
  df$Color <- ifelse(df$AVG.Log2.Ratio >= 0.58 & df$Qvalue < myQval, "Red", "Grey")
  df$Color <- ifelse(df$AVG.Log2.Ratio <= -0.58 & df$Qvalue < myQval, "Blue", df$Color)
  df <- df %>% filter(Genes %in% SASP$Genes)
  down <- paste(count(df[df$Color == "Blue",]), "Down", sep = "_")
  up <- paste(count(df[df$Color == "Red",]), "Up", sep = "_")
  temp <- ggplot(data = df, aes(x = AVG.Log2.Ratio, y = minuslogqval, col = Color, text = Genes)) + 
    geom_point() +
    # geom_text_repel(aes(x = AVG.Log2.Ratio, y = minuslogqval)) +  
    geom_vline(xintercept = c(-0.58, 0.58), col = "black", linetype = "dashed") +
    geom_hline(yintercept = -log10(myQval), col = "black", linetype = "dashed") +
    scale_color_manual(values = mycolors) +
    ylab("-Log10(q-value)") +
    scale_x_continuous(name = "Log2(fold change)", limits = c(-7,7), labels = c(-6,-3,0,3,6), breaks = c(-6,-3,0,3,6)) +
    ylim(c(0,250))+
    #facet_wrap(~factor(Comparison..group1.group2.), nrow = 2) +
    theme_classic() +
    theme(axis.title =element_text(size = 20, color = "black"),
          axis.text = element_text(size = 18, color = "black"),
          legend.position = "none")
  ggsave(filename = paste("SASP",
                          str_split(x, " / ")[[1]][1], 
                          "vs", 
                          str_split(x, " / ")[[1]][2],
                          "volcano",
                          myname, down, up,
                          sep = "_"), 
         device = tiff, path = "output", dpi = "print", width = 5, height = 5, units = 'in')  
}
sapply(unique(all_comparisons$Comparison..group1.group2.), volcano.ECM)


#Extracellular Volcano

volcano.ECM <- function(x){
  require(ggplot2)
  require(stringr)
  require(dplyr)
  df <- all_comparisons[all_comparisons$Comparison..group1.group2. == x,]
  df$Color <- ifelse(df$AVG.Log2.Ratio >= 0.58 & df$Qvalue < myQval, "Red", "Grey")
  df$Color <- ifelse(df$AVG.Log2.Ratio <= -0.58 & df$Qvalue < myQval, "Blue", df$Color)
  df <- df %>% filter(grepl('extracellular', GO.Cellular.Component, ignore.case = T))
  down <- paste(count(df[df$Color == "Blue",]), "Down", sep = "_")
  up <- paste(count(df[df$Color == "Red",]), "Up", sep = "_")
  temp <- ggplot(data = df, aes(x = AVG.Log2.Ratio, y = minuslogqval, col = Color, text = Genes)) + 
    geom_point() +
    # geom_text_repel(aes(x = AVG.Log2.Ratio, y = minuslogqval)) +  
    geom_vline(xintercept = c(-0.58, 0.58), col = "black", linetype = "dashed") +
    geom_hline(yintercept = -log10(myQval), col = "black", linetype = "dashed") +
    scale_color_manual(values = mycolors) +
    ylab("-Log10(q-value)") +
    scale_x_continuous(name = "Log2(fold change)", limits = c(-7,7), labels = c(-6,-3,0,3,6), breaks = c(-6,-3,0,3,6)) +
    ylim(c(0,250))+
    #facet_wrap(~factor(Comparison..group1.group2.), nrow = 2) +
    theme_classic() +
    theme(axis.title =element_text(size = 20, color = "black"),
          axis.text = element_text(size = 18, color = "black"),
          legend.position = "none")
  ggsave(filename = paste("Extracellular",
                          str_split(x, " / ")[[1]][1], 
                          "vs", 
                          str_split(x, " / ")[[1]][2],
                          "volcano",
                          myname, down, up,
                          sep = "_"), 
         device = tiff, path = "output", dpi = "print", width = 5, height = 5, units = 'in')  
}
