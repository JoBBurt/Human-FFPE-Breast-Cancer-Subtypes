library(ggplot2)

#Import data
data <- read.csv("output/Output_Data/ECM_abundance.csv", header = TRUE, stringsAsFactors = FALSE)

data$Condition <- factor(data$Condition, levels = c("DF", "Meta", "TN", "HER2", "TP", "ERPR"))
data$Category <- factor(data$Category, 
                        levels = c("Other","Collagens", "ECM Glycoproteins", "Proteoglycans",
                                   "ECM-affiliated Proteins", "ECM Regulators", "Secreted Factors"
                                   ))

tiff("output/Histo_RelativeAbundance_MatrisomeDB_JB8-14_2023_0510.tiff", res = 300, height = 5, width = 7, units = "in")
ggplot(data = data, aes(x = Condition, y = Values, fill = Category)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c(alpha("black", 0.5),
                              alpha('#FD6467',1), alpha('#FD6467',0.5), alpha('#FD6467',0.2),
                               alpha('#5B1A18',1), alpha('#5B1A18',0.5), alpha('#5B1A18',0.2))) +
  scale_y_continuous(name = "Relative abundance", labels = scales::percent) +
  scale_x_discrete(labels = c("DF", "MBC", "TNBC", "HER2", "LumB", "LumA")) +
  theme_classic() +
  theme(axis.title.y = element_text(size = 20, color = "black"),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 18, color = "black"),
        axis.text.x = element_text(angle = 90, size = 18, color = "black"),
        legend.text = element_text(size = 16, color = "black"),
        legend.title = element_blank(),
        legend.position = "right")
dev.off()
