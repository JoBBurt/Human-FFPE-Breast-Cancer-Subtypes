
library(clusterProfiler)
library(AnnotationHub)
library(org.Hs.eg.db)
library(ggplot2)
library(dplyr)
library(VennDiagram)


# data import
df <- read.csv("data/DAVID/Meta_0p001_KEGG.csv", 
               stringsAsFactors = F)

###############   Meta All Kegg Pathway Dot Plot Dotplot with ggplot##############
# filter data for plotting
data <- df %>%
grid.newpage()   

# top10 <- top_n(data, n = 10) # select the top10 pathways to display if applicable


tiff(filename="output/Meta_KEGG_SASP_0p001_text.tiff", width=1200, height=600)

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
   scale_y_continuous(limits=c(0, 15),breaks=seq(0,15,5)) +
   scale_size_continuous(range = c(0,15)) +
  geom_text(size = 5, angle = 45)

dev.off()


###
S1 <- ggplot(data, aes(x=Log2E , y=Log10FDR), 
             size=Count, color=Bonferroni) + geom_point(alpha = 0.8) + 
  theme_classic()
S1 

S1 = S1+scale_color_gradient(low = "red2",  
                             high = "mediumblue", space = "Lab",
                             limit = c(min(data$q.value), 
                                       max(data$q.value)))
S1 + scale_y_discrete(name="") +
  scale_size(range = c(2, 6)) +
  theme(axis.title.x = element_text(size=16),
        axis.text.x  = element_text(size=14, angle = 90,
                                    vjust = 0.5, hjust=1),
        axis.text.y = element_text(size=10))

dev.off()

# data import
df <- read.csv("data/DAVID/Meta_0p001_KEGG_SASP.csv", 
               stringsAsFactors = F)

###############   Meta Kegg SASP Pathway Dot Plot Dotplot with ggplot##############
# filter data for plotting
data <- df %>%
  grid.newpage()   

# top10 <- top_n(data, n = 10) # select the top10 pathways to display if applicable


tiff(filename="output/Meta_KEGG_SASP_0p001.tiff", width=1200, height=600)

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
  scale_size_continuous(range = c(0,15)) 
  #geom_text(size = 5, angle = 45)

dev.off()


# data import
df <- read.csv("data/DAVID/Meta_0p001_KEGG_ECM.csv", 
               stringsAsFactors = F)

###############   Meta Kegg SASP Pathway Dot Plot Dotplot with ggplot##############
# filter data for plotting
data <- df %>%
  grid.newpage()   

# top10 <- top_n(data, n = 10) # select the top10 pathways to display if applicable


tiff(filename="output/Meta_KEGG_ECM_0p001_labeled.tiff", width=1200, height=600)

ggplot(df, aes(x=Log2E , y=Log10FDR, 
               size=Count, fill=Bonferroni, label=Label)) + 
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
  scale_size_continuous(range = c(0,15))+
  geom_text(size = 4, position = "dodge")

dev.off()


