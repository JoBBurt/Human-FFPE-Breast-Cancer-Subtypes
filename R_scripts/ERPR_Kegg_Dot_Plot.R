
library(clusterProfiler)
library(AnnotationHub)
library(org.Hs.eg.db)
library(ggplot2)
library(dplyr)
library(VennDiagram)


# data import
df <- read.csv("data/DAVID/All_Raw/ERPR_0p001_KEGG_All.csv", 
               stringsAsFactors = F)

###############    All Kegg Pathway Dot Plot Dotplot with ggplot##############
# filter data for plotting
data <- df %>%
grid.newpage()   

# top10 <- top_n(data, n = 10) # select the top10 pathways to display if applicable


tiff(filename="output/Pathways/ERPR_KEGG_All_0p001.tiff", width=1200, height=600)

ggplot(df, aes(x=Log2E , y=Log10FDR, 
             size=Count, fill=Bonferroni, label=Term)) + 
  geom_point(alpha = 0.8, shape = 21, color = "black") + 
  theme_classic(base_size = 28) +
  scale_fill_gradient(low = "red2",  
                      high = "blue", space = "Lab",
                      name = 'P value',
                      limit = c(min(df$Bonferroni), 
                                max(df$Bonferroni))) +
   xlab(bquote('Log'[2]~'(Enrichment)')) +
   ylab(bquote('-Log'[10]~'(FDR)')) + 
   scale_x_continuous(limits=c(-3, 3),breaks=seq(-3,3,.5)) +
   scale_y_continuous(limits=c(0, 10),breaks=seq(0,10,2)) +
   scale_size_continuous(range = c(0,10)) 
   #geom_text(size = 5, angle = 45)

dev.off()


# data import
df <- read.csv("data/DAVID/SASP_Raw/ERPR_0p001_KEGG_SASP.csv", 
               stringsAsFactors = F)

###############    Kegg SASP Pathway Dot Plot Dotplot with ggplot##############
# filter data for plotting
data <- df %>%
  grid.newpage()   

# top10 <- top_n(data, n = 10) # select the top10 pathways to display if applicable


tiff(filename="output/Pathways/ERPR_KEGG_SASP_0p001.tiff", width=1200, height=600)

ggplot(df, aes(x=Log2E , y=Log10FDR, 
               size=Count, fill=Bonferroni, label=Term)) + 
  geom_point(alpha = 0.8, shape = 21, color = "black") + 
  theme_classic(base_size = 28) +
  scale_fill_gradient(low = "red2",  
                      high = "blue", space = "Lab",
                      name = 'P value',
                      limit = c(min(df$Bonferroni), 
                                max(df$Bonferroni))) +
  xlab(bquote('Log'[2]~'(Enrichment)')) +
  ylab(bquote('-Log'[10]~'(FDR)')) + 
  scale_x_continuous(limits=c(-6, 6),breaks=seq(-6,6,1)) +
  scale_y_continuous(limits=c(0, 10),breaks=seq(0,10,2)) +
  scale_size_continuous(range = c(0,10)) 
  #geom_text(size = 5, angle = 45)

dev.off()


# data import
df <- read.csv("data/DAVID/ECM_raw/ERPR_0p001_KEGG_ECM.csv", 
               stringsAsFactors = F)

###############    Kegg SASP Pathway Dot Plot Dotplot with ggplot##############
# filter data for plotting
data <- df %>%
  grid.newpage()   

# top10 <- top_n(data, n = 10) # select the top10 pathways to display if applicable


tiff(filename="output/Pathways/ERPR_KEGG_ECM_0p001_text.tiff", width=1200, height=600)

ggplot(df, aes(x=Log2E , y=Log10FDR, 
               size=Count, fill=Bonferroni, label=Term)) + 
  geom_point(alpha = 0.8, shape = 21, color = "black") + 
  theme_classic(base_size = 28) +
  scale_fill_gradient(low = "red2",  
                      high = "blue", space = "Lab",
                      name = 'P value',
                      limit = c(min(df$Bonferroni), 
                                max(df$Bonferroni))) +
  xlab(bquote('Log'[2]~'(Enrichment)')) +
  ylab(bquote('-Log'[10]~'(FDR)')) + 
  scale_x_continuous(limits=c(-4, 4),breaks=seq(-4,4,1)) +
  scale_y_continuous(limits=c(0, 20),breaks=seq(0,20,5)) +
  scale_size_continuous(range = c(0,20)) +
  geom_text(size = 5, angle = 45)

dev.off()


