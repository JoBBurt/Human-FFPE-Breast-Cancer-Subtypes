library(matrixStats)
# Barplots for Specific Proteins
progression <- c("DF", "Meta", "TN", "HER2", "TP", "ERPR")

df <- as.data.frame(my_pro)

mf <- sapply(split.default(df, names(df)), rowMeans)
mf <- mf[,progression]
mf <- as.data.frame(mf)
sf <- data.frame(Gene = rownames(my_pro),
                 Cont = rowSds(my_pro[,c(rep("Cont", 14))]),
                 MPBC = rowSds(my_pro[,c(rep("MPBC", 14))]),
                 TNBC = rowSds(my_pro[,c(rep("TNBC", 14))]),
                 HER2 = rowSds(my_pro[,c(rep("HER2", 14))]),
                 LumB = rowSds(my_pro[,c(rep("LumB", 14))]),
                 LumA = rowSds(my_pro[,c(rep("LumA", 14))])
)

mf  <-mf %>%  mutate(Gene = rownames(mf)) %>% 
  pivot_longer(!Gene, names_to = "Group", values_to = "Mean")

sf <- sf %>% pivot_longer(!Gene, names_to = "Group", values_to = "SD")

plotFrame <- full_join(mf, sf)
plotFrame$Group <- factor(plotFrame$Group, levels = progression)



mybar <- function(myGene){
    
  plotFrame %>% filter(Gene == myGene) %>%
    ggplot(aes(x = Group, y = Mean/1000, fill = Group))+
    geom_bar(stat = "identity", fill = mycol) +
    geom_errorbar(aes(ymax = Mean/1000 + SD/1000, ymin = Mean/1000 - SD/1000), width = 0.5) +
    ggtitle(myGene) +
    xlab("Subtype") +
    scale_x_discrete(labels = c("Cont", "MPBC", "TNBC", "HER2", "LumB", "LumA"))+
    ylab("Average Abundance (x10^3)") +
    theme(legend.position = "none",axis.text = element_text(size = 16, color = "black"), 
          axis.title = element_text(size = 18, color = "black"),
          legend.text = element_text(size = 14),
          axis.text.x = element_text(angle = 90),
          title = element_text(size = 20, color = "black")) +
    theme_bw()
  
  ggsave(filename = paste(myGene, "_Barplot.tiff"), 
         device = tiff, path = "barplot_of_interest", dpi = "print", width = 5, height = 5, units = 'in')  
  
}

mycol <- c("red", "darkgreen", "blue", "deepskyblue", "orange", "purple")
mybar("SERPINH1")
mybar("PRELP")
mybar("POSTN")
mybar("HSPG2")
mybar("TIMP1")
mybar("GPC1")
mybar("TGFB1")
mybar("SELENBP1")


mybar("EPX")
