library(ggplot2)
library(dplyr)
library(scales)
library(hrbrthemes)
library(tidyr)

dq <- all_comparisons[all_comparisons$Absolute.AVG.Log2.Ratio >= 0.58,]
dq<-all_comparisons
#1.0 Common Genes in Other CIACs
common <- c("DSP","NID1","LAMB2","LAMA5","LAMC1","LAMA4","KRT8","POSTN","DES",
            "MYH11","FLNA","COL6A1","COL6A2","VWA1","COL6A3","FLNC","APCS",
            "PRELP","FGG","HSP90AA1","COL14A1","MPO","PFN1")
cf <- dq %>% filter(Genes %in% common)
cf$Condition.Numerator <- factor(cf$Condition.Numerator, levels = c("Meta", "TN", "HER2", "ERPR", "TP"))

tiff("output/Heatmap_common_2022_1127.tiff", res = 300, height = 5, width = 7, units = "in")
ggplot(cf, aes(x=Condition.Numerator, y=Genes, fill=AVG.Log2.Ratio)) +
  coord_equal(ratio = 1/2) + # aspect ratio y/x
  geom_tile(color = "white") +
  scale_y_discrete(limits = rev(levels(cf$Genes))) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", na.value = "white", 
                       "Log2(FC)",
                       limits=c(-6,6), oob=squish) + # squish and limits are used in controlling color saturation # library(scales)
  geom_point(aes(x=Condition.Numerator, y=Genes, size = -log10(Qvalue))) +
  scale_size_binned(range = c(0.5,3.5))+
  xlab("Subtype") +
  ylab("Proteins") +
  theme_ipsum() +
  theme(axis.text = element_text(size = 16, color = "black"), 
        axis.title = element_text(size = 18, color = "black"),
        legend.text = element_text(size = 14),
        axis.text.x = element_text(angle = 90))
dev.off()

#2.0 Basement Membrane
bm <- read.csv("Basement_Membrane_Data/Basement_Membrane.csv",
               stringsAsFactors = F)
bm <- bm %>% filter(Category == "Basement membrane")

bf <- dq %>% filter(Genes %in% bm$Gene.symbol)
bf$Condition.Numerator <- factor(bf$Condition.Numerator, levels = c("Meta", "TN", "HER2", "ERPR", "TP"))

tiff("output/Heatmap_BasementMembrane_2022_1128.tiff", res = 300, height = 8, width = 7, units = "in")
ggplot(bf, aes(x=Condition.Numerator, y=Genes, fill=AVG.Log2.Ratio)) +
  coord_equal(ratio = 1/2) + # aspect ratio y/x
  geom_tile(color = "white") +
  scale_y_discrete(limits = rev(levels(bf$Genes))) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", na.value = "white", 
                       "Log2(FC)",
                       limits=c(-6,6), oob=squish) + # squish and limits are used in controlling color saturation # library(scales)
  geom_point(aes(x=Condition.Numerator, y=Genes, size = -log10(Qvalue))) +
  scale_size_binned(range = c(0.25,2))+
  xlab("Subtype") +
  ylab("Proteins") +
  theme_ipsum() +
  theme(axis.text = element_text(size = 16, color = "black"), 
        axis.title = element_text(size = 18, color = "black"),
        legend.text = element_text(size = 14),
        axis.text.x = element_text(angle = 90),
        axis.text.y = element_text(size = 10))
dev.off()

#3.0 Small Leucine-Rich Proteins
slrp <- read.csv("data/small_leucine_rich_proteins.csv", stringsAsFactors = F)

sf <- dq %>% filter(Genes %in% slrp$Genes)
sf$Condition.Numerator <- factor(sf$Condition.Numerator, levels = c("Meta", "TN", "HER2", "ERPR", "TP"))

tiff("output/Heatmap_SLRPs_2022_1128.tiff", res = 300, height = 4, width = 7, units = "in")
ggplot(sf, aes(x=Condition.Numerator, y=Genes, fill=AVG.Log2.Ratio)) +
  coord_equal(ratio = 1/2) + # aspect ratio y/x
  geom_tile(color = "white") +
  scale_y_discrete(limits = rev(levels(sf$Genes))) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", na.value = "white", 
                       "Log2(FC)",
                       limits=c(-6,6), oob=squish) + # squish and limits are used in controlling color saturation # library(scales)
  geom_point(aes(x=Condition.Numerator, y=Genes, size = -log10(Qvalue))) +
  scale_size_binned(range = c(0.5,3))+
  xlab("Subtype") +
  ylab("Proteins") +
  theme_ipsum() +
  theme(axis.text = element_text(size = 16, color = "black"), 
        axis.title = element_text(size = 18, color = "black"),
        legend.text = element_text(size = 14),
        axis.text.x = element_text(angle = 90),
        axis.text.y = element_text(size = 10))
dev.off()

#4.0 Keratins

kf <- dplyr::filter(dq, grepl('KRT', Genes))
kf$Condition.Numerator <- factor(kf$Condition.Numerator, levels = c("Meta", "TN", "HER2", 
                                                                    "ERPR", "TP"))
kf$Genes <- factor(kf$Genes, levels = c("KRT1","KRT2","KRT3","KRT4","KRT5","KRT6A",
                                        "KRT6B","KRT6C","KRT7","KRT8","KRT9","KRT10",
                                        "KRT13","KRT14","KRT15","KRT16","KRT17","KRT18",
                                        "KRT19","KRT20","KRT31","KRT73","KRT75","KRT77",
                                        "KRT78","KRT80","KRT81","KRT85"))
tiff("output/Heatmap_Keratinss_2022_1128.tiff", res = 300, height = 7, width = 7, units = "in")
ggplot(kf, aes(x=Condition.Numerator, y=Genes, fill=AVG.Log2.Ratio)) +
  coord_equal(ratio = 1/2) + # aspect ratio y/x
  geom_tile(color = "white") +
  scale_y_discrete(limits = rev(levels(kf$Genes))) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", na.value = "white", 
                       "Log2(FC)",
                       limits=c(-6,6), oob=squish) + # squish and limits are used in controlling color saturation # library(scales)
  geom_point(aes(x=Condition.Numerator, y=Genes, size = -log10(Qvalue))) +
  scale_size_binned(range = c(0.25,3))+
  xlab("Subtype") +
  ylab("Proteins") +
  theme_ipsum() +
  theme(axis.text = element_text(size = 16, color = "black"), 
        axis.title = element_text(size = 18, color = "black"),
        legend.text = element_text(size = 14),
        axis.text.x = element_text(angle = 90),
        axis.text.y = element_text(size = 10))
dev.off()

#5.0 Conserved ECM Signatures
sig <- c("COL4A1","POSTN","NID1","NID2","AGRN","HSPG2","PRELP","DCN","OGN",
         "SERPINH1","MMP8","MMP9","MMP14","MUC1","MUC2","MUC5AC","MUC5B","MUC16")
cf <- dq %>% filter(Genes %in% sig)
cf$Condition.Numerator <- factor(cf$Condition.Numerator, levels = c("Meta", "TN", "HER2", "ERPR", "TP"))

tiff("output/Heatmap_signatures_2022_1128.tiff", res = 300, height = 5, width = 7, units = "in")
ggplot(cf, aes(x=Condition.Numerator, y=Genes, fill=AVG.Log2.Ratio)) +
  coord_equal(ratio = 1/2) + # aspect ratio y/x
  geom_tile(color = "white") +
  scale_y_discrete(limits = rev(levels(cf$Genes))) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", na.value = "white", 
                       "Log2(FC)",
                       limits=c(-6,6), oob=squish) + # squish and limits are used in controlling color saturation # library(scales)
  geom_point(aes(x=Condition.Numerator, y=Genes, size = -log10(Qvalue))) +
  scale_size_binned(range = c(0.25,3))+
  xlab("Subtype") +
  ylab("Proteins") +
  theme_ipsum() +
  theme(axis.text = element_text(size = 16, color = "black"), 
        axis.title = element_text(size = 18, color = "black"),
        legend.text = element_text(size = 14),
        axis.text.x = element_text(angle = 90))
dev.off()

#6.0 Desmosome

desmosome <- c("CDSN","DSC", "DSC1", "DSG", "DSG1", "DSG2", "DSG3", "DSG4", "DSP",
               "JUP", "DP3", "CTNNG", "PLEC", "PKP", "PKP1", "PKP2", "PKP3", "PPL")
cf <- dq %>% filter(Genes %in% desmosome)
cf$Condition.Numerator <- factor(cf$Condition.Numerator, levels = c("Meta", "TN", "HER2", "ERPR", "TP"))

tiff("output/Heatmap_Desmosome_2022_1128.tiff", res = 300, height = 5, width = 7, units = "in")
ggplot(cf, aes(x=Condition.Numerator, y=Genes, fill=AVG.Log2.Ratio)) +
  coord_equal(ratio = 1/2) + # aspect ratio y/x
  geom_tile(color = "white") +
  scale_y_discrete(limits = rev(levels(cf$Genes))) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", na.value = "white", 
                       "Log2(FC)",
                       limits=c(-6,6), oob=squish) + # squish and limits are used in controlling color saturation # library(scales)
  geom_point(aes(x=Condition.Numerator, y=Genes, size = -log10(Qvalue))) +
  scale_size_binned(range = c(0.25,3))+
    xlab("Subtype") +
  ylab("Proteins") +
  theme_ipsum() +
  theme(axis.text = element_text(size = 16, color = "black"), 
        axis.title = element_text(size = 18, color = "black"),
        legend.text = element_text(size = 14),
        axis.text.x = element_text(angle = 90))
dev.off()

#6.0 SASP

cf <- dq %>% filter(Genes %in% SASP$Genes)
cf$Condition.Numerator <- factor(cf$Condition.Numerator, levels = c("Meta", "TN", "HER2", "ERPR", "TP"))

tiff("output/Heatmap_SASP_2022_1128.tiff", res = 300, height = 20, width = 7, units = "in")
ggplot(cf, aes(x=Condition.Numerator, y=Genes, fill=AVG.Log2.Ratio)) +
  coord_equal(ratio = 1/2) + # aspect ratio y/x
  geom_tile(color = "white") +
  scale_y_discrete(limits = rev(levels(cf$Genes))) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", na.value = "white", 
                       "Log2(FC)",
                       limits=c(-6,6), oob=squish) + # squish and limits are used in controlling color saturation # library(scales)
  geom_point(aes(x=Condition.Numerator, y=Genes, size = -log10(Qvalue))) +
  scale_size_binned(range = c(0.25,3))+
  xlab("Subtype") +
  ylab("Proteins") +
  theme_ipsum() +
  theme(axis.text = element_text(size = 16, color = "black"), 
        axis.title = element_text(size = 18, color = "black"),
        legend.text = element_text(size = 14),
        axis.text.x = element_text(angle = 90))
dev.off()

