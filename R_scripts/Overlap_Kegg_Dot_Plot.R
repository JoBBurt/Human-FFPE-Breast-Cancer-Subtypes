
library(clusterProfiler)
library(AnnotationHub)
library(org.Hs.eg.db)
library(ggplot2)
library(dplyr)
library(VennDiagram)


# data import
df <- read.csv("data/DAVID/OverLap/Overlap_KEGG.csv", 
               stringsAsFactors = F)

###############    All Kegg Pathway Dot Plot Dotplot with ggplot##############
# filter data for plotting
data <- df %>%
grid.newpage()   

# top10 <- top_n(data, n = 10) # select the top10 pathways to display if applicable


tiff(filename="output/Pathways/Overlap_text.tiff", width=1200, height=600)

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
   scale_y_continuous(limits=c(0, 6),breaks=seq(0,6,1)) +
   scale_size_continuous(range = c(0,10)) +
  geom_text(size = 5, angle = 45)

dev.off()

## Unique Proteins for each condition

## Meta Unique
# data import
df <- read.csv("data/DAVID/Unique/Meta_Unique.csv", 
               stringsAsFactors = F)

###############    All Kegg Pathway Dot Plot Dotplot with ggplot##############
# filter data for plotting
data <- df %>%
   grid.newpage()   

# top10 <- top_n(data, n = 10) # select the top10 pathways to display if applicable


tiff(filename="output/Pathways/Meta_Unique_text.tiff", width=1200, height=600)

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
   scale_x_continuous(limits=c(-4, 4),breaks=seq(-4,4,.5)) +
   scale_y_continuous(limits=c(0, 12),breaks=seq(0,12,2)) +
   scale_size_continuous(range = c(0,10)) +
   geom_text(size = 5, angle = 45)

dev.off()

## TNBC Unique
# data import
df <- read.csv("data/DAVID/Unique/TN_Unique.csv", 
               stringsAsFactors = F)

###############    All Kegg Pathway Dot Plot Dotplot with ggplot##############
# filter data for plotting
data <- df %>%
   grid.newpage()   

# top10 <- top_n(data, n = 10) # select the top10 pathways to display if applicable


tiff(filename="output/Pathways/TN_Unique_text.tiff", width=1200, height=600)

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
   scale_x_continuous(limits=c(-3, 3.5),breaks=seq(-3,3.5,.5)) +
   scale_y_continuous(limits=c(0, 1),breaks=seq(0,1,0.25)) +
   scale_size_continuous(range = c(0,10)) +
   geom_text(size = 5, angle = 45)

dev.off()


## HER Unique
# data import
df <- read.csv("data/DAVID/Unique/HER_Unique.csv", 
               stringsAsFactors = F)

###############    All Kegg Pathway Dot Plot Dotplot with ggplot##############
# filter data for plotting
data <- df %>%
   grid.newpage()   

# top10 <- top_n(data, n = 10) # select the top10 pathways to display if applicable


tiff(filename="output/Pathways/HER_Unique.tiff", width=1200, height=600)

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
   scale_x_continuous(limits=c(-4.5, 4.5),breaks=seq(-4.5,4.5,.5)) +
   scale_y_continuous(limits=c(0, 17),breaks=seq(0,17,2)) +
   scale_size_continuous(range = c(0,10)) 
   #geom_text(size = 5, angle = 45)

dev.off()


## LumA Unique
# data import
df <- read.csv("data/DAVID/Unique/LumA_Unique.csv", 
               stringsAsFactors = F)

###############    All Kegg Pathway Dot Plot Dotplot with ggplot##############
# filter data for plotting
data <- df %>%
   grid.newpage()   

# top10 <- top_n(data, n = 10) # select the top10 pathways to display if applicable


tiff(filename="output/Pathways/LumA_Unique.tiff", width=1200, height=600)

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
   scale_x_continuous(limits=c(-3.5, 3.5),breaks=seq(-3.5,3.5,.5)) +
   scale_y_continuous(limits=c(0, 1),breaks=seq(0,1,.25)) +
   scale_size_continuous(range = c(0,10)) 
  # geom_text(size = 5, angle = 45)

dev.off()


## LumB Unique
# data import
df <- read.csv("data/DAVID/Unique/LumB_Unique.csv", 
               stringsAsFactors = F)

###############    All Kegg Pathway Dot Plot Dotplot with ggplot##############
# filter data for plotting
data <- df %>%
   grid.newpage()   

# top10 <- top_n(data, n = 10) # select the top10 pathways to display if applicable


tiff(filename="output/Pathways/LumB_Unique_text.tiff", width=1200, height=600)

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
   scale_x_continuous(limits=c(-4.5, 4.5),breaks=seq(-4.5,4.5,.5)) +
   scale_y_continuous(limits=c(0, 1),breaks=seq(0,1,.25)) +
   scale_size_continuous(range = c(0,10)) +
 geom_text(size = 5, angle = 45)

dev.off()