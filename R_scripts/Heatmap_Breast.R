library(ggplot2)
library(dplyr)
library(scales)
library(hrbrthemes)
library(tidyr)

dq <- all_comparisons[all_comparisons$Absolute.AVG.Log2.Ratio >= 0.58 & all_comparisons$Qvalue <= myQval,]
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
bf$Condition.Numerator <- factor(bf$Condition.Numerator, levels = c("Meta", "TN", "HER2", "TP", "ERPR"))

tiff("output/Heatmap_BasementMembrane_2023_0510.tiff", res = 300, height = 12, width = 7, units = "in")
ggplot(bf, aes(x=Condition.Numerator, y=Genes, fill=AVG.Log2.Ratio)) +
  coord_equal(ratio = 1/2) + # aspect ratio y/x
  geom_tile(color = "white") +
  scale_y_discrete(limits = rev(levels(bf$Genes))) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", na.value = "white", 
                       "Log2(FC)",
                       limits=c(-6,6), oob=squish) + # squish and limits are used in controlling color saturation # library(scales)
 # geom_point(aes(x=Condition.Numerator, y=Genes, size = -log10(Qvalue))) +
 # scale_size_binned(range = c(0.25,2))+
  scale_x_discrete(labels = c("MBC", "TNBC", "HER2", "LumB", "LumA")) +
  xlab("Subtype") +
  ylab("Proteins") +
  theme_ipsum() +
  theme(axis.text = element_text(size = 18, color = "black"), 
        axis.title = element_text(size = 18, color = "black"),
        legend.text = element_text(size = 18),
        axis.text.x = element_text(angle = 90),
        axis.text.y = element_text(size = 12))
dev.off()

#3.0 Small Leucine-Rich Proteins
slrp <- read.csv("data/small_leucine_rich_proteins.csv", stringsAsFactors = F)

sf <- dq %>% filter(Genes %in% slrp$Genes)
sf$Condition.Numerator <- factor(sf$Condition.Numerator, levels = c("Meta", "TN", "HER2", "TP", "ERPR"))

tiff("output/Heatmap_SLRPs_2023_0510.tiff", res = 300, height = 5, width = 5, units = "in")
ggplot(sf, aes(x=Condition.Numerator, y=Genes, fill=AVG.Log2.Ratio)) +
  coord_equal(ratio = 1/2) + # aspect ratio y/x
  geom_tile(color = "white") +
  scale_y_discrete(limits = rev(levels(sf$Genes))) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", na.value = "white", 
                       "Log2(FC)",
                       limits=c(-6,6), oob=squish) + # squish and limits are used in controlling color saturation # library(scales)
  #geom_point(aes(x=Condition.Numerator, y=Genes, size = -log10(Qvalue))) +
  #scale_size_binned(range = c(0.5,3))+
  xlab("Subtype") +
  scale_x_discrete(labels = c("MBC", "TNBC", "HER2", "LumB", "LumA")) +
  ylab("Proteins") +
  theme_ipsum() +
  theme(axis.text = element_text(size = 18, color = "black"), 
        axis.title = element_text(size = 18, color = "black"),
        legend.text = element_text(size = 18),
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
cf$Condition.Numerator <- factor(cf$Condition.Numerator, levels = c("Meta", "TN", "HER2", "TP", "ERPR"))

tiff("output/Heatmap_signatures_2023_0512.tiff", res = 300, height = 7, width = 7, units = "in")
ggplot(cf, aes(x=Condition.Numerator, y=Genes, fill=AVG.Log2.Ratio)) +
  coord_equal(ratio = 1/2) + # aspect ratio y/x
  geom_tile(color = "white") +
  scale_y_discrete(limits = rev(levels(cf$Genes))) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", na.value = "white", 
                       "Log2(FC)",
                       limits=c(-6,6), oob=squish) + # squish and limits are used in controlling color saturation # library(scales)
  #geom_point(aes(x=Condition.Numerator, y=Genes, size = -log10(Qvalue))) +
  #scale_size_binned(range = c(0.25,3))+
  xlab("Subtype") +
  scale_x_discrete(labels = c("MBC", "TNBC", "HER2", "LumB", "LumA")) +
  ylab("Proteins") +
  theme_ipsum() +
  theme(axis.text = element_text(size = 18, color = "black"), 
        axis.title = element_text(size = 18, color = "black"),
        legend.text = element_text(size = 18),
        axis.text.x = element_text(angle = 90))
dev.off()

#6.0 Desmosome

desmosome <- c("CDSN","DSC", "DSC1", "DSG", "DSG1", "DSG2", "DSG3", "DSG4", "DSP",
               "JUP", "DP3", "CTNNG", "PLEC", "PKP", "PKP1", "PKP2", "PKP3", "PPL")
cf <- dq %>% filter(Genes %in% desmosome)
cf$Condition.Numerator <- factor(cf$Condition.Numerator, levels = c("Meta", "TN", "HER2", "TP", "ERPR"))

tiff("output/Heatmap_Desmosome_2023_0510.tiff", res = 300, height = 5, width =5, units = "in")
ggplot(cf, aes(x=Condition.Numerator, y=Genes, fill=AVG.Log2.Ratio)) +
  coord_equal(ratio = 1/2) + # aspect ratio y/x
  geom_tile(color = "white") +
  scale_y_discrete(limits = rev(levels(cf$Genes))) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", na.value = "white", 
                       "Log2(FC)",
                       limits=c(-6,6), oob=squish) + # squish and limits are used in controlling color saturation # library(scales)
  #geom_point(aes(x=Condition.Numerator, y=Genes, size = -log10(Qvalue))) +
  #scale_size_binned(range = c(0.25,3))+
  xlab("Subtype") +
  scale_x_discrete(labels = c("MBC", "TNBC", "HER2", "LumB", "LumA")) +
  ylab("Proteins") +
  theme_ipsum() +
  theme(axis.text = element_text(size = 18, color = "black"), 
        axis.title = element_text(size = 18, color = "black"),
        legend.text = element_text(size = 18),
        axis.text.x = element_text(angle = 90))
dev.off()

#6.0 SASP
SASP <- c("TNC", "VCAN", "POSTN", "FN1", "FKBP10", "MSN", "HSP90B1", "HSP90AB1", "HSP90AA1", "SERPINH1", "GDI2",
          "VIM", "NID2", "HSPG2", "IGFBP2", "IGFBP4", "KRT2", "KRT10")
cf <- dq %>% filter(Genes %in% SASP)
cf$Condition.Numerator <- factor(cf$Condition.Numerator, levels = c("Meta", "TN", "HER2", "TP", "ERPR"))
cf$Genes <- factor(cf$Genes, levels = SASP)

tiff("output/Heatmap_SASP_2022_0510.tiff", res = 300, height = 7, width = 7, units = "in")
ggplot(cf, aes(x=Condition.Numerator, y=Genes, fill=AVG.Log2.Ratio)) +
  coord_equal(ratio = 1/2) + # aspect ratio y/x
  geom_tile(color = "white") +
  scale_y_discrete(limits = rev(levels(cf$Genes))) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", na.value = "white", 
                       "Log2(FC)",
                       limits=c(-6,6), oob=squish) + # squish and limits are used in controlling color saturation # library(scales)
 # geom_point(aes(x=Condition.Numerator, y=Genes, size = -log10(Qvalue))) +
#  scale_size_binned(range = c(0.25,3))+
  xlab("Subtype") +
  scale_x_discrete(labels = c("MBC", "TNBC", "HER2", "LumB", "LumA")) +
  ylab("Proteins") +
  theme_ipsum() +
  theme(axis.text = element_text(size = 18, color = "black"), 
        axis.title = element_text(size = 18, color = "black"),
        legend.text = element_text(size = 18),
        axis.text.x = element_text(angle = 90))
dev.off()


#6.1 SASP Poster

cf <- dq %>% filter(Genes %in% SASP$Genes)
cf$Condition.Numerator <- factor(cf$Condition.Numerator, levels = c("Meta", "TN", "HER2", "ERPR", "TP"))
myGenes <- c("TNC", "VCAN", "POSTN", "FN1", "FKBP10", "MSN", "HSP90B1", "HSP90AB1", "HSP90AA1", "SERPINH1", 
             "GDI2", "VIM", "NID2", "HSPG2", "IGFBP2", "IGFBP4", "KRT2", "KRT10")
cf<- cf %>% filter(Genes %in% myGenes)
cf$Genes<- factor(cf$Genes, levels = myGenes)

tiff("output/Heatmap_SASP_2022_1128_selected.tiff", res = 300, height = 6, width = 7, units = "in")
ggplot(cf, aes(x=Condition.Numerator, y=Genes, fill=AVG.Log2.Ratio)) +
  coord_equal(ratio = 1/2) + # aspect ratio y/x
  geom_tile(color = "white") +
  scale_y_discrete(limits = rev(levels(cf$Genes))) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", na.value = "white", 
                       "Log2(FC)",
                       limits=c(-6,6), oob=squish) + # squish and limits are used in controlling color saturation # library(scales)
  geom_point(aes(x=Condition.Numerator, y=Genes, size = -log10(Qvalue))) +
  scale_size_binned(range = c(0.1,3))+
  xlab("Subtype") +
  ylab("Proteins") +
  theme_ipsum() +
  theme(axis.text = element_text(size = 16, color = "black"), 
        axis.title = element_text(size = 18, color = "black"),
        legend.text = element_text(size = 14),
        axis.text.x = element_text(angle = 90))
dev.off()

# Overlap of significantly altered proteins
dq <- all_comparisons[all_comparisons$Absolute.AVG.Log2.Ratio >= 0.58 & 
                        all_comparisons$Qvalue <= 0.001,]
overlap <- unique(dq[with(dq, ave(seq_along(Comparison..group1.group2.), Genes, FUN = length) == 
                            length(unique(Comparison..group1.group2.))),])
neg_overlap <- overlap[overlap$AVG.Log2.Ratio < 0,]
tiff("output/Heatmaps/Heatmap_Overlap_2023_0510.tiff", res = 300,units = "in")
ggplot(neg_overlap, aes(x=Condition.Numerator, y=Genes, fill=AVG.Log2.Ratio)) +
  coord_equal(ratio = 1/2) + # aspect ratio y/x
  geom_tile(color = "white") +
  scale_y_discrete(limits = rev(levels(neg_overlap$Genes))) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", na.value = "white", 
                       "Log2(FC)",
                       limits=c(-6,6), oob=squish) + # squish and limits are used in controlling color saturation # library(scales)
  #geom_point(aes(x=Condition.Numerator, y=Genes, size = -log10(Qvalue))) +
  #scale_size_binned(range = c(0.1,3))+
  xlab("Subtype") +
  ylab("Proteins") +
  theme_ipsum() +
  theme(axis.text = element_text(size = 16, color = "black"), 
        axis.title = element_text(size = 18, color = "black"),
        legend.text = element_text(size = 14),
        axis.text.x = element_text(angle = 90))
dev.off()

#Ovelap Heatmap

library(openxlsx)
overlap <- read.xlsx("data/Overlap/Overlap_same_direction.xlsx")

oMat <- as.matrix(overlap[,2:ncol(overlap)])
rownames(oMat) <- overlap$Protein
colnames(oMat) <- c("MPBC", "TNBC", "HER2", "LumB", "LumA")
library(gplots)

# Compute the pairwise distances between rows
dist_matrix <- dist(overlap, method = "euclidean")

# Cluster the rows based on their distances
row_clusters <- hclust(dist_matrix, method = "ward.D")

# Plot the heatmap
heatmap.2(oMat,
          dendrogram = "row",
          Rowv = as.dendrogram(row_clusters),
          col = bluered(75),
          margins = c(5,10),
          key = TRUE,
          keysize = 1,
          trace = "none")

# New Overlap Heatmap

overlap <- read.csv(file = "data/Overlap/overlap_all.csv", stringsAsFactors = F)
cf <- dq %>% filter(Genes %in% overlap$Genes)
cf <- cf %>% select(c("Genes", "AVG.Log2.Ratio", "Condition.Numerator"))
cf$Condition.Numerator <- factor(cf$Condition.Numerator, levels = c("Meta", "TN", "HER2", "ERPR", "TP"))
cf$AVG.Log2.Ratio <- as.numeric(cf$AVG.Log2.Ratio)


df <- cf %>% 
  group_by(Genes, Condition.Numerator) %>% 
  summarise(AVG.Log2.Ratio = mean(AVG.Log2.Ratio)) %>% 
  pivot_wider(names_from = Condition.Numerator, values_from = AVG.Log2.Ratio)

oMat <- as.matrix(df[,2:6])
rownames(oMat) <- df$Genes
colnames(oMat) <- c("MPBC", "TNBC", "HER2", "LumB", "LumA")
dist_matrix <- dist(oMat, method = "euclidean")

# Cluster the rows based on their distances
row_clusters <- hclust(dist_matrix, method = "ward.D")
pdf("heatmap.pdf")
heatmap.2(oMat,
          dendrogram = "row",
          Rowv = as.dendrogram(row_clusters),
          col = bluered(75),
          margins = c(5,10),
          key = TRUE,
          keysize = 1,
          trace = "none")
dev.off()

## Look at sign. altered proteins with specific pattern of fold changes
keep_rows <- oMat[, "MPBC"] >= oMat[, "TNBC"] &
             oMat[, "TNBC"] >= oMat[, "HER2"] &
             oMat[, "HER2"] >= oMat[, "LumB"] &
             oMat[, "LumB"] >= oMat[, "LumA"]

fMat <- oMat[keep_rows,]

keep_rows <- oMat[, "MPBC"] <= oMat[, "TNBC"] &
             oMat[, "TNBC"] <= oMat[, "HER2"] &
             oMat[, "HER2"] <= oMat[, "LumB"] &
             oMat[, "LumB"] <= oMat[, "LumA"]

fMat <- oMat[keep_rows,]
fMat <- rbind(fMat, dMat)
dist_matrix <- dist(fMat, method = "euclidean")
row_clusters <- hclust(dist_matrix, method = "ward.D")
custom_order <- c("MPBC", "TNBC", "HER2", "LumB", "LumA")
pdf("heatmap_fade.pdf")
heatmap.2(fMat,
          dendrogram = "row",
          Rowv = as.dendrogram(row_clusters),
          col = bluered(75),
          margins = c(5,10),
          key = TRUE,
          keysize = 1,
          trace = "none")
dev.off()


#6.0 Desmosome

common <- c("MYO6","PPK2", "ERC1", "CHGA", "KRT17", "CES1", "CPB1", "CES1P1", "AGR3",
               "HSD11B2")
cf <- dq %>% filter(Genes %in% common)
cf$Condition.Numerator <- factor(cf$Condition.Numerator, levels = c("Meta", "TN", "HER2", "ERPR", "TP"))

tiff("output/Heatmap_Conserved_2023_0324.tiff", res = 300, height = 4, width = 7, units = "in")
ggplot(cf, aes(x=Condition.Numerator, y=Genes, fill=AVG.Log2.Ratio)) +
  coord_equal(ratio = 1/2) + # aspect ratio y/x
  geom_tile(color = "white") +
  scale_y_discrete(limits = rev(levels(cf$Genes))) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", na.value = "white", 
                       "Log2(FC)",
                       limits=c(-4,4), oob=squish) + # squish and limits are used in controlling color saturation # library(scales)
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

#modules of interest

royalblue <- c("C1RL", "CRTAP", "CTHRC1", "IGFBP3", "PLOD2", "RPL22C1")
green <- c("ACTBL2", "CPVL", "CTSB", "PLOD3", "PTX3", "QKI", "SH3BGRL3")
white <- c("DEFA1", "EPX", "MMP9", "MPO", "PRTN3")
darkorange <- c("MIEN1", "PLIN2", "TAOK2", "GRN", "GRB7")
lightcyan <- c("BNIP1", "ERBB2", "H3F3A", "RPL7", "RPL8", "RPL15", "RPL18", "RPL18A",
               "RPL21", "RPL30", "RPL32", "RPL34", "RPL35A")
saddlebrown <- c("RAB13", "SYPL1", "TMEM63A")
paleturquoise <- c("APOO", "ANXA9", "ABAT", "CRIP1", "DCTN3", "DCTPP1", "HSPB1", 
                   "MAPT", "SLC25A24", "SSH3", "YES1")
greyup <- c("ABHD14B", "CA2", "CLIC6", "CPB1", "DCCK1", "ECM1")
greydown <- c("ACYP2", "APP", "COQ7", "DCAF13", "GGH", "LSG1", "MFGE8", "TRIP11", "VWA1")
cyanup <- c("AGR2", "CYB5R1", "HSD17B8", "ISOC1", "L2HGDH", "PLEKHF1", "PVALB", "QDPR",
            "RHOT2", "RMND1", "SNAP29", "STARD10")
cyandown <- c("4-Sep", "ACOT8", "ANK2", "ARFIP2", "ARML10", "ATP55", "CBR4", "CDC73",
              "COX5A", "COX5B", "DCD", "EARS2", "GCTP", "HADHB", "LRRFID1", "MIPEP",
              "MRP52", "MT-CO2", "MTDH", "NADK2", "NDUFA3", "NDUFA5", "NDUFA9",
              "NDUFA10", "NDUFA13", "NDUFB4", "NDUFB5", "NDUFB10", "NDUFS1", "NDUFS2",
              "NDUFS5", "NDUFS7", "NDUFS8", "NDUFU1", "NEDD4", "NUMA1", "TMEM11", "TUFM",
              "VDAC2")

all <- c(royalblue, green, white, darkorange, lightcyan, saddlebrown, paleturquoise,
         greyup, greydown, cyanup, cyandown)

cf <- dq %>% filter(Genes %in% royalblue)
cf$Condition.Numerator <- factor(cf$Condition.Numerator, levels = c("Meta", "TN", "HER2", "TP", "ERPR"))

tiff("Heatmap_royalblue_2023_0522.tiff", res = 300, height = 4, width = 4, units = "in")
cf %>%
  ggplot(., aes(x=Condition.Numerator, y=Genes, fill=AVG.Log2.Ratio)) +
  coord_equal(ratio = 1/2) + # aspect ratio y/x
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", na.value = "white", 
                       "Log2(FC)",
                       limits=c(-4,4), oob=squish) + # squish and limits are used in controlling color saturation # library(scales)
 # geom_point(aes(x=Condition.Numerator, y=Genes, size = -log10(Qvalue))) +
  #scale_size_binned(range = c(0.25,3))+
  xlab("Subtype") +
  scale_x_discrete(labels = c("MBC", "TNBC", "HER2", "LumB", "LumA")) +
  ylab("Proteins") +
  theme_ipsum() +
  theme(axis.text = element_text(size = 16, color = "black"), 
        axis.title = element_text(size = 18, color = "black"),
        legend.text = element_text(size = 14),
        axis.text.x = element_text(angle = 90))
dev.off()

cf <- dq %>% filter(Genes %in% green)
cf$Condition.Numerator <- factor(cf$Condition.Numerator, levels = c("Meta", "TN", "HER2", "TP", "ERPR"))

tiff("Heatmap_green_2023_0522.tiff", res = 300, height = 4, width = 4, units = "in")
cf %>%
  ggplot(., aes(x=Condition.Numerator, y=Genes, fill=AVG.Log2.Ratio)) +
  coord_equal(ratio = 1/2) + # aspect ratio y/x
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", na.value = "white", 
                       "Log2(FC)",
                       limits=c(-4,4), oob=squish) + # squish and limits are used in controlling color saturation # library(scales)
  # geom_point(aes(x=Condition.Numerator, y=Genes, size = -log10(Qvalue))) +
  #scale_size_binned(range = c(0.25,3))+
  xlab("Subtype") +
  scale_x_discrete(labels = c("MBC", "TNBC", "HER2", "LumB", "LumA")) +
  ylab("Proteins") +
  theme_ipsum() +
  theme(axis.text = element_text(size = 16, color = "black"), 
        axis.title = element_text(size = 18, color = "black"),
        legend.text = element_text(size = 14),
        axis.text.x = element_text(angle = 90))
dev.off()

cf <- dq %>% filter(Genes %in% white)
cf$Condition.Numerator <- factor(cf$Condition.Numerator, levels = c("Meta", "TN", "HER2", "TP", "ERPR"))

tiff("Heatmap_white_2023_0522.tiff", res = 300, height = 4, width = 4, units = "in")
cf %>%
  ggplot(., aes(x=Condition.Numerator, y=Genes, fill=AVG.Log2.Ratio)) +
  coord_equal(ratio = 1/2) + # aspect ratio y/x
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", na.value = "white", 
                       "Log2(FC)",
                       limits=c(-4,4), oob=squish) + # squish and limits are used in controlling color saturation # library(scales)
  # geom_point(aes(x=Condition.Numerator, y=Genes, size = -log10(Qvalue))) +
  #scale_size_binned(range = c(0.25,3))+
  xlab("Subtype") +
  scale_x_discrete(labels = c("MBC", "TNBC", "HER2", "LumB", "LumA")) +
  ylab("Proteins") +
  theme_ipsum() +
  theme(axis.text = element_text(size = 16, color = "black"), 
        axis.title = element_text(size = 18, color = "black"),
        legend.text = element_text(size = 14),
        axis.text.x = element_text(angle = 90))
dev.off()

cf <- dq %>% filter(Genes %in% all)
cf$Condition.Numerator <- factor(cf$Condition.Numerator, levels = c("Meta", "TN", "HER2", "TP", "ERPR"))
cf <- cf %>% mutate(Module = ifelse(Genes %in% royalblue, "royalblue", 
                              ifelse(Genes %in% green, "green",
                               ifelse(Genes %in% white, "white", 
                                ifelse(Genes %in% darkorange, "darkorange",
                                 ifelse(Genes %in% lightcyan, "lightcyan",
                                  ifelse(Genes %in% saddlebrown, "saddlebrown", 
                                    ifelse(Genes %in% paleturquoise, "paleturquoise",
                                     ifelse(Genes %in% greyup, "greyup",
                                      ifelse(Genes %in% greydown, "greydown",
                                       ifelse(Genes %in% cyanup, "cyanup", 
                                        ifelse(Genes %in% cyandown, "cyandown", NA))))))))))))


tiff("Heatmap_all_2023_0523.tiff", res = 300, height = 11, width = 4, units = "in")
cf %>%
  ggplot(., aes(x=Condition.Numerator, y=Genes, fill=AVG.Log2.Ratio, group = Module)) +
  coord_equal(ratio = 1/2) + # aspect ratio y/x
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", na.value = "white", 
                       "Log2(FC)",
                       limits=c(-4,4), oob=squish) + # squish and limits are used in controlling color saturation # library(scales)
  # geom_point(aes(x=Condition.Numerator, y=Genes, size = -log10(Qvalue))) +
  #scale_size_binned(range = c(0.25,3))+
  xlab("Subtype") +
  scale_x_discrete(labels = c("MBC", "TNBC", "HER2", "LumB", "LumA")) +
  ylab("Proteins") +
  theme_ipsum() +
  theme(axis.text = element_text(size = 16, color = "black"), 
        axis.title = element_text(size = 18, color = "black"),
        legend.text = element_text(size = 14),
        axis.text.x = element_text(angle = 90))
dev.off()

tiff("Heatmap_HER2_2023_0614.tiff", res = 300, height = 4, width = 8, units = "in")
cf %>% filter(Module %in% c("darkorange", "lightcyan")) %>%
  mutate(GroupedGenes = gsub("\\.Module", "", interaction(Genes, Module))) %>%
  ggplot(., aes(x = Condition.Numerator, y = GroupedGenes, fill = AVG.Log2.Ratio)) +
  coord_equal(ratio = 1/2) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", na.value = "white", "Log2(FC)",
                       limits = c(-4, 4), oob = squish) +
  xlab("Subtype") +
  scale_x_discrete(labels = c("MPBC", "TNBC", "HER2", "LumB", "LumA")) +
  ylab("Genes") +
  theme_ipsum() +
  theme(axis.text = element_text(size = 16, color = "black"),
        axis.title = element_text(size = 18, color = "black"),
        legend.text = element_text(size = 14),
        axis.text.x = element_text(angle = 90))

dev.off()

tiff("Heatmap_LumB_2023_0614.tiff", res = 300, height = 4, width = 8, units = "in")
cf %>% filter(Module %in% c("paleturquoise")) %>%
  mutate(GroupedGenes = gsub("\\.Module", "", interaction(Genes, Module))) %>%
  ggplot(., aes(x = Condition.Numerator, y = GroupedGenes, fill = AVG.Log2.Ratio)) +
  coord_equal(ratio = 1/2) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", na.value = "white", "Log2(FC)",
                       limits = c(-4, 4), oob = squish) +
  xlab("Subtype") +
  scale_x_discrete(labels = c("MPBC", "TNBC", "HER2", "LumB", "LumA")) +
  ylab("Genes") +
  theme_ipsum() +
  theme(axis.text = element_text(size = 16, color = "black"),
        axis.title = element_text(size = 18, color = "black"),
        legend.text = element_text(size = 14),
        axis.text.x = element_text(angle = 90))

dev.off()

tiff("Heatmap_LumA_down_2023_0614.tiff", res = 300, height = 4, width = 8, units = "in")
cf %>% filter(Module %in% c("greydown")) %>%
  mutate(GroupedGenes = gsub("\\.Module", "", interaction(Genes, Module))) %>%
  ggplot(., aes(x = Condition.Numerator, y = GroupedGenes, fill = AVG.Log2.Ratio)) +
  coord_equal(ratio = 1/2) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", na.value = "white", "Log2(FC)",
                       limits = c(-4, 4), oob = squish) +
  xlab("Subtype") +
  scale_x_discrete(labels = c("MPBC", "TNBC", "HER2", "LumB", "LumA")) +
  ylab("Genes") +
  theme_ipsum() +
  theme(axis.text = element_text(size = 16, color = "black"),
        axis.title = element_text(size = 18, color = "black"),
        legend.text = element_text(size = 14),
        axis.text.x = element_text(angle = 90))

dev.off()

tiff("Heatmap_Luma_up_2023_0614.tiff", res = 300, height = 4, width = 8, units = "in")
cf %>% filter(Module %in% c("greyup")) %>%
  mutate(GroupedGenes = gsub("\\.Module", "", interaction(Genes, Module))) %>%
  ggplot(., aes(x = Condition.Numerator, y = GroupedGenes, fill = AVG.Log2.Ratio)) +
  coord_equal(ratio = 1/2) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", na.value = "white", "Log2(FC)",
                       limits = c(-4, 4), oob = squish) +
  xlab("Subtype") +
  scale_x_discrete(labels = c("MPBC", "TNBC", "HER2", "LumB", "LumA")) +
  ylab("Genes") +
  theme_ipsum() +
  theme(axis.text = element_text(size = 16, color = "black"),
        axis.title = element_text(size = 18, color = "black"),
        legend.text = element_text(size = 14),
        axis.text.x = element_text(angle = 90))

dev.off()


# Stress Extrinsic Expression Profile (SEEP) Fordyce et al 2012
SEEP <- c("HPGD", "PTGS1", "PTGR1", "PTGES3", "PTGR2", "FAM213B", "PTGES2", 
          "PTGFRN", "TNXB", "TNC", "FN1", "COL1A1", "HYOU1", "LDHA", "LDHB",
          "LDHD", "IL6", "IL8", "VNN2", "VCAM1", "TACSTD2", "STAT3", "STAT1",
          "STAT6", "STAT5A", "STAT5B", "STAT2")
cf <- filter(dq, Genes %in% SEEP)
cf$Condition.Numerator <- factor(cf$Condition.Numerator, levels = c("Meta", "TN", "HER2", "TP", "ERPR"))


tiff("Heatmap_SEEP_sign_2023_0616.tiff", res = 300, height = 4, width = 8, units = "in")
cf %>% 
  ggplot(., aes(x = Condition.Numerator, y = Genes, fill = AVG.Log2.Ratio)) +
  coord_equal(ratio = 1/2) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", na.value = "white", "Log2(FC)",
                       limits = c(-6, 6), oob = squish) +
  xlab("Subtype") +
  scale_x_discrete(labels = c("MPBC", "TNBC", "HER2", "LumB", "LumA")) +
  ylab("Genes") +
  theme_ipsum() +
  theme(axis.text = element_text(size = 16, color = "black"),
        axis.title = element_text(size = 18, color = "black"),
        legend.text = element_text(size = 14),
        axis.text.x = element_text(angle = 90))

dev.off() 

#ER Stress-Response

df <- u2osCan %>% filter(assignment == "ER")
ER <- c("UGGT1", "SERPINH1", "PPIB", "PDIA4", "PDIA3", "P4HB", "HYOU1",
        "HSPA5", "HSP90B1", "ERP44", "ERP29", "CDK5RAP3")
cf <- filter(dq, Genes %in% ER)


cf$Condition.Numerator <- factor(cf$Condition.Numerator, levels = c("Meta", "TN", "HER2", "TP", "ERPR"))

cf <- cf %>% filter(AVG.Log2.Ratio >0)

tiff("Heatmap_ER_stress_sign_2023_0616.tiff", res = 300, height = 4, width = 8, units = "in")
cf %>% 
  ggplot(., aes(x = Condition.Numerator, y = Genes, fill = AVG.Log2.Ratio)) +
  coord_equal(ratio = 1/2) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", na.value = "white", "Log2(FC)",
                       limits = c(-6, 6), oob = squish) +
  xlab("Subtype") +
  scale_x_discrete(labels = c("MPBC", "TNBC", "HER2", "LumB", "LumA")) +
  ylab("Genes") +
  theme_ipsum() +
  theme(axis.text = element_text(size = 16, color = "black"),
        axis.title = element_text(size = 18, color = "black"),
        legend.text = element_text(size = 14),
        axis.text.x = element_text(angle = 90))

dev.off() 





