#1.0 LIBRARIES ----

library(openxlsx)
library(tidyverse)
library(scales)


#2.0 User Inputs & Modifications ----
## Set the working directory
setwd("~/Desktop/CRUK_Storming_Cancer/Breast_FFPE/JB8-14_2rep")

## modify with your batch name
batch <- "JB8-JB14"

## What organism is it? (Options currently are human or mouse)
organism <- "human"
## modify with your search information
search_info <- "  from 7 biological replicates with 2 technical replicates per biological replicate searched with the panHuman Library"
comparison_info <- "Breast Carcinoma Subtypes vs. Disease Free after Search with PanHuman Library"

## How many technical and biological replicates are in the batch
tech_rep <- 2
bio_rep <- 7
## I order my data by condition, this could be different for your dataset
## This is the treatment factor used in the PCA and heatmap

treat <- c(rep("DF", tech_rep * bio_rep), rep("Meta",tech_rep * bio_rep), 
           rep("ERPR", tech_rep * bio_rep), rep("HER2", tech_rep * bio_rep), 
           rep("TP", tech_rep * bio_rep), rep("TN", tech_rep * bio_rep))

## There needs to be a color for each treatment factor
heatcolor <- c("red", "orange", "yellow", "green", "blue", "purple")

## mypattern is used to split the string of the file name
## "JB\\d_\\d+" splits the data file into "JB8_01" etc = Batch_SampleNumber
mypattern <- "JB\\d_\\d+"

## Which q-value filter would you like to include in the spreadsheet?
add_q0.05 <- FALSE
add_q0.01 <- FALSE
add_q0.001 <- TRUE

## Would you like a Violin Plot
### Needs colors to be specified
Violin <- TRUE
vf_man <- alpha(heatcolor,0.7)
vs_man <- heatcolor
# 3.0  User Input GET DATA ----
protein <- read.csv("data/22_1018_JB8-JB14_2rep_panHuman_v01_Report_Birgit_Protein Quant_Pivot (Pivot).csv",
                    stringsAsFactors = F)
all_comparisons <- read.csv("data/22_1018_JB8-JB14_2rep_panHuman_v01_candidates.csv",
                            stringsAsFactors = F)

# 4.0 Modify Data ----
all_comparisons$minuslogqval <- -1*log10(all_comparisons$Qvalue)
theme_set(theme_bw(base_size = 18))

## This function prepares a matrix of the protein data
prepmat <- function(df, mypattern){
  require("stringr")
  temp <- data.matrix(dplyr::select(df, contains(c("PG.Quantity"))))
  string <- colnames(temp)
  names <- str_extract(string = string, pattern = mypattern)
  colnames(temp) <- names
  rownames(temp) <- df$PG.Genes
  return(temp)
}

proIDs <- paste(comma(nrow(protein)), 
                "Protein Groups Identified with ≥ 2 Unique Peptides",
                sep = " ")
## my_pro is used for the PCA and Heatmap
my_pro <- prepmat(protein, mypattern)
colnames(my_pro) <- treat
# 5.0 Volcano Plot Function ----
mycolors <- c("Blue", "Red", "Gray")
names(mycolors) <- c("Blue", "Red", "Gray")

volcano <- function(x){
  require(ggplot2)
  require(stringr)
  df <- all_comparisons[all_comparisons$Comparison..group1.group2. == x,]
  df$Color <- ifelse(df$AVG.Log2.Ratio >= 0.58 & df$Qvalue < myQval, "Red", "Grey")
  df$Color <- ifelse(df$AVG.Log2.Ratio <= -0.58 & df$Qvalue < myQval, "Blue", df$Color)
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
    #facet_wrap(~factor(Comparison..group1.group2.), nrow = 2) +
    theme_classic() +
    theme(axis.title =element_text(size = 20, color = "black"),
          axis.text = element_text(size = 18, color = "black"),
          legend.position = "none")
  ggsave(filename = paste(str_split(x, " / ")[[1]][1], 
                          "vs", 
                          str_split(x, " / ")[[1]][2],
                          "volcano",
                          myname, down, up,
                          sep = "_"), 
         device = tiff, path = "output", dpi = "print", width = 5, height = 5, units = 'in')  
}

# 6.0 PCA, Heatmap, Correlation and Violin Plots ----

pca_heat_corr <- function(df, pca_name, heat_name, corr_name, mycol){
  require(gplots)
  require(ggplot2)
  require(ggrepel)
  require(corrplot)
  require(dplyr)
  meds<-apply(df, 2, median, na.rm=TRUE)
  nMat<-sweep(df, 2, meds/mean(meds), FUN="/")
  
  pcMat<-nMat
  pcMat<-pcMat[complete.cases(pcMat),]
  pcMat[pcMat == 0] <-1
  pcRes<-prcomp(t(log2(pcMat)), center = TRUE, scale. = TRUE)
  pcSum <- summary(pcRes)
  PC1label <- paste0("PC1, ",
                     round(100 * pcSum$importance["Proportion of Variance", "PC1"],1),
                     "% of variance")
  PC2label <- paste0("PC2, ",
                     round(100 * pcSum$importance["Proportion of Variance", "PC2"],1), 
                     "% of variance")
  
  treat <- colnames(df) # needs to be changed as data from replicates is updated
  
  treatment <- factor(treat)
  myRamp<-colorRampPalette(colors=c("#0571b0", "#92c5de", "#f7f7f7", "#f4a582", "#ca0020"))
  
  pcPlotFrame<-data.frame(treatment = treatment,
                          sample = colnames(nMat),
                          pcRes$x[,1:5])
  pcPlotFrame %>% 
    ggplot(aes(PC1, PC2,  color = treatment, shape = treatment, label = treatment))+ #label = sample
    geom_point(size=1.8) +
    scale_x_continuous(name=PC1label) +
    scale_y_continuous(name=PC2label) +
    geom_text_repel(size = 5) +
    scale_color_manual(values = mycol) +
    theme(legend.position = 'none') +
    #labs(title = "A") #+ 
    stat_ellipse(aes(color = paste0(treatment)))
  #ggsave("output/fig2a.emf", width=180, height=480, units="mm")
  ggsave(filename = pca_name, device = tiff, path = "output")
  
  myheatcol <- rep(heatcolor, each = tech_rep * bio_rep)
  tiff(file = heat_name)
  heatmap.2(t(scale(t(log10(pcMat)))), col = myRamp, trace = 'none', labRow = FALSE,
            ColSideColors = myheatcol)
  dev.off()
  
  col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
  M <- cor(df)
  
  tiff(file = paste(corr_name,".tiff",sep=""))
  corrplot(M, method="color", col=col(2000),
           tl.col = "black")
  dev.off()
  p<-corrplot(M, method="color", col=col(2000),
             tl.col = "black")
  write.csv(p$corr, paste(corr_name,".csv",sep=""))
}

pca_heat_corr(df = my_pro, pca_name = paste(batch, "pca", sep = "_"), 
         heat_name = paste("output/", batch, "_heatmap.tiff", sep = ""),
         corr_name = paste("output/", batch, "_correlation", sep = ""),
         mycol = heatcolor)

Violin_plot <- function(violin_name){
  #Import des donnees
  Import <- read.csv(paste("output/", batch, "_correlation_edited", ".csv", sep = ""), sep = ",", dec = ".", header = TRUE, na = "NA", stringsAsFactors = FALSE)
  
  tiff(violin_name, res = 300, height = 4, width = 5, units = "in")
  ggplot(Import, aes(x=factor(Group,levels = unique(treat)),
                     y=Values, fill=Group, color = Group)) +
    geom_violin(trim = T, scale = "width") +
    scale_y_continuous(name = "Pearson correlation coefficients", limits = c(0,1), 
                       labels = c(0,0.2,0.4,0.6,0.8,1), breaks = c(0,0.2,0.4,0.6,0.8,1)) +
    stat_summary(fun = mean, geom = "point", shape = 23, fill = "black", size = 3) +
    geom_jitter(shape = 16, position = position_jitter(0.1), color = alpha("black", 0.2)) +
    scale_fill_manual(values = vf_man) +
    scale_color_manual(values = vs_man) +
    theme_bw() +
    theme(axis.title.y =element_text(size = 18, color = "black"),
          axis.text.y = element_text(size = 16, color = "black"),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          legend.position = "bottom")
  dev.off()
}

if(Violin == TRUE){
  Violin_plot(violin_name = paste("output/", batch, "_violin.tiff", sep = ""))
}

# 7.0 LOPIT Plots ----

Mouse_LOPIT_ref <- function(){
  require(pRoloc)
  require(pRolocdata)
  require(tidyverse)
  require(reshape2)
  require(ggplot2)
  require(gridExtra)
  require(Rtsne)
  require(gplots)
  require(RColorBrewer)
  require(plotly)
  require(ggpubr)
  require(dplyr)
  theme_set(theme_bw(base_size = 12))
  update_geom_defaults("point", list(size = 0.7))
  ## Reference Plot
  ## https://pubmed.ncbi.nlm.nih.gov/26754106/
  data("hyperLOPIT2015")
  set.seed(11)
  p <- plot2D(hyperLOPIT2015, method = "t-SNE")
  u2os <- data.frame(Accession = row.names(p), p, assignment = as.character(fData(hyperLOPIT2015)$final.assignment), stringsAsFactors = FALSE)
  u2os$Accession <- sapply(strsplit(u2os$Accession, "-"),'[', 1)
  u2osSummary <- u2os %>% group_by(assignment) %>% count()
  textFrame <- data.frame(x = c(-5, 5, -15, -33, 34, 7, 8, -5, 44, 25, 0, 17, -20, -25), 
                          y = c(18, 5, -6, -2, -20, -39.5, -27, -41, 4, 30, 45, -3, -31, 23), 
                          text = c("40S Ribosome", "60S Ribosome", "Actin Cytoskeleton", "Cytosol", "Endoplasmic Reticulum/\nGolgi Apparatus", "Endosome", "Extracellular Matrix", "Lysosome", "Mitochondria", "Nucleus Chromatin", "Nucleus Non-Chromatin", "Peroxisome", "Plasma Membrane", "Proteasome"))
  
  mycolors <- c("#E31A1C", "#D95F02", "#70b38d", "#A6CEE3", "#B15928", "#B2DF8A","#3328b1", "#FB9A99", "#1B9E77", "#FDBF6F", "#FF7F00", "#6A3D9A", "#CAB2D6", "#dbdb4b", "#3328b1")
  
  hyperLOPIT <- u2os %>%
    mutate(annotated = assignment != "unknown") %>%
    ggplot(aes(x = Dimension.1, y = Dimension.2)) +
    geom_point(data = function(x){x[!(x$annotated), ]}, color = grey(0.9)) +
    geom_point(data = function(x){x[(x$annotated), ]}, aes(color = assignment)) +
    geom_text(data = textFrame, aes(x = x, y = y, label = text), size = 3.5) +
    scale_color_manual(values = mycolors) +
    labs(color = "Localization", x = "t-SNE Dim. 1", y = "t-SNE Dim. 2") +
    theme(#axis.text.x = element_blank(),
      #axis.text.y = element_blank(),
      #axis.ticks = element_blank(), 
      legend.position = 'none')
  
  ggsave(hyperLOPIT, file = "output/LOPIT_mouse_reference.tiff", dpi = 300)
  return(u2os)
  }  

Mouse_LOPIT <- function(x){
  require(pRoloc)
  require(pRolocdata)
  require(tidyverse)
  require(reshape2)
  require(ggplot2)
  require(gridExtra)
  require(Rtsne)
  require(gplots)
  require(RColorBrewer)
  require(plotly)
  require(ggpubr)
  require(dplyr)
  theme_set(theme_bw(base_size = 12))
  update_geom_defaults("point", list(size = 0.7)) 
 ## Fold Change Plots
    df <- all_comparisons[all_comparisons$Comparison..group1.group2. == x,]
    df$minuslogqval <- -1*log10(df$Qvalue)
    df$Color <- ifelse(df$AVG.Log2.Ratio >= 0.58 & df$Qvalue < myQval, "Red", "Grey")
    df$Color <- ifelse(df$AVG.Log2.Ratio <= -0.58 & df$Qvalue < myQval, "Blue", df$Color)
    df$Label <- NA
    u2os$inCan <- u2os$Accession %in% df$UniProtIds #1721 proteins
    u2osCan <- df %>%
      full_join(u2os, by = c("UniProtIds" = "Accession"))
    
    textFrame <- data.frame(x = c(-5, 5, -15, -33, 28, 7, 7, -5, 40, 25, 0, 17, -20, -28), 
                            y = c(18, 5, -6, -2, -20, -39.5, -26, -41, 4, 25, 45, -3, -31, 23), 
                            text = c("40S R", "60S R", "AC", "Cyt", "ER/GA", "End", "EM", "Lys", "Mito", "Nuc-Chr", "Nuc Non-Chr", "Per", "PM", "Pro"))
    #Filter out non-statistically significant data points
    myCan <- u2osCan[!is.na(u2osCan$Absolute.AVG.Log2.Ratio),]
    myCan <- myCan[myCan$Color != "Grey",]
    
    canLOPIT <- u2osCan %>%
      mutate(identified = !is.na(u2osCan$inCan)) %>%
      ggplot(aes(x = Dimension.1, y = Dimension.2)) +
      geom_point(alpha = 0.1) +
      geom_point(data = myCan, aes(x = Dimension.1, y = Dimension.2, size = Absolute.AVG.Log2.Ratio),
                 color = myCan$Color, alpha = 0.3) +
      #scale_size_manual(values = c(0.18, 0.375, .75, 1.5)) +
      labs(color = "Localization", title = x, x = "t-SNE Dim. 1", y = "t-SNE Dim. 2") +
      theme(#axis.text.x = element_blank(),
        #axis.text.y = element_blank(),
        #axis.ticks = element_blank(),
        legend.position = c(0.88, 0.83)) +
      guides(size=guide_legend(title="|log2(FC)| Ratio")) +
      geom_text(data = textFrame, aes(x = x, y = y, label = text), size = 4.5)
    ggsave(canLOPIT, path = "output", dpi = 300, 
           filename = paste(str_split(x, " / ")[[1]][1], 
                            "vs", 
                            str_split(x, " / ")[[1]][2],
                            "LOPIT_Mouse", myname,
                            sep = "_"),)
}
  
Human_LOPIT_ref <- function(){
  require(pRoloc)
  require(pRolocdata)
  require(tidyverse)
  require(reshape2)
  require(ggplot2)
  require(gridExtra)
  require(Rtsne)
  require(gplots)
  require(RColorBrewer)
  require(plotly)
  require(ggpubr)
  require(dplyr)
  theme_set(theme_bw(base_size = 12))
  update_geom_defaults("point", list(size = 0.7))
  ## Reference Plot
  # Load data from Thul et. al. 2017
  data("hyperLOPITU2OS2017")
  
  # use T-SNE algorithm to generate plot dimensions and assign points to dimensions
  set.seed(11)
  p <- plot2D(hyperLOPITU2OS2017, method = "t-SNE")
  
  # Formatting the data 
  u2os <- data.frame(Accession = row.names(p), p, assignment = as.character(fData(hyperLOPITU2OS2017)$assignment), stringsAsFactors = FALSE)
  u2os$Accession <- sapply(strsplit(u2os$Accession, "-"),'[', 1)
  
  ## Add assignment name to plot
  textFrame <- data.frame(x = c(30, 8, -25, -40, 38, -27, -31, 52, -25, 28, 14, 25), 
                          y = c(22, 30, 35, 14, 3, -30, -12, 0, -2, 8, -28, -23), 
                          text = c("Cytosol", "ER", "Golgi", "Lysosome", "Mitochondria", 
                                   "Nucleus", "Nucleus-Chromatin", "Peroxisome", "Plasma Membrane", 
                                   "Proteasome", "40S Ribosome", "60S Ribosome"))
  
  # plotting the hyperlopit data
  hyperLOPIT <- u2os %>%
    mutate(annotated = assignment != "unknown") %>%
    ggplot(aes(x = Dimension.1, y = Dimension.2)) +
    geom_point(data = function(x){x[!(x$annotated), ]}, color = grey(0.9)) +
    geom_point(data = function(x){x[(x$annotated), ]}, aes(color = assignment)) +
    geom_text(data = textFrame, aes(x = x, y = y, label = text), size = 3.5) +
    scale_color_manual(values = c(brewer.pal(12, "Paired"))) +
    labs(color = "Localization", x = "t-SNE Dim. 1", y = "t-SNE Dim. 2") +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(), 
          legend.position = 'none')
  
  # save the hyperLOPIT plot
  ggsave(hyperLOPIT, file = "output/LOPIT_Human_reference.tiff", dpi = 300)
  
  return(u2os)
}  

Human_LOPIT <- function(x){
  require(pRoloc)
  require(pRolocdata)
  require(tidyverse)
  require(reshape2)
  require(ggplot2)
  require(gridExtra)
  require(Rtsne)
  require(gplots)
  require(RColorBrewer)
  require(plotly)
  require(ggpubr)
  require(dplyr)
  theme_set(theme_bw(base_size = 12))
  update_geom_defaults("point", list(size = 0.7))
  
  ## Fold Change Plots
  df <- all_comparisons[all_comparisons$Comparison..group1.group2. == x,]
  df$minuslogqval <- -1*log10(df$Qvalue)
  df$Color <- ifelse(df$AVG.Log2.Ratio >= 0.58 & df$Qvalue < myQval, "Red", "Grey")
  df$Color <- ifelse(df$AVG.Log2.Ratio <= -0.58 & df$Qvalue < myQval, "Blue", df$Color)
  df$Label <- NA
  u2os$inCan <- u2os$Accession %in% df$UniProtIds #1721 proteins
  u2osCan <- df %>%
    full_join(u2os, by = c("UniProtIds" = "Accession"))
  
  #Filter out non-statistically significant data points
  myCan <- u2osCan[!is.na(u2osCan$Absolute.AVG.Log2.Ratio),]
  myCan <- myCan[myCan$Color != "Grey",]

  # Generate the LOPIT Plot
  canLOPIT <- u2osCan %>%
    mutate(identified = !is.na(u2osCan$inCan)) %>%
    ggplot(aes(x = Dimension.1, y = Dimension.2)) +
    geom_point(alpha = 0.1) +
    geom_point(data = myCan,
               aes(x = Dimension.1, y = Dimension.2, size = Absolute.AVG.Log2.Ratio),
               color = myCan$Color, alpha = 0.3) +
    #scale_size_manual(values = c(0.18, 0.375, .75, 1.5)) +
    labs(color = "Localization", title = paste(str_split(x, " / ")[[1]][1], 
                                               "vs", 
                                               str_split(x, " / ")[[1]][2],
                                               sep = " "),
         x = "t-SNE Dim. 1", y = "t-SNE Dim. 2") +
    theme(#axis.text.x = element_blank(),
      #axis.text.y = element_blank(),
      #axis.ticks = element_blank(),
      legend.position = c(0.88, 0.83)) +
    guides(size=guide_legend(title="|log2(FC)| Ratio")) +
    geom_text(data = textFrame, aes(x = x, y = y, label = text), size = 4.5)
  
  # Save the LOPIT Plot
  ggsave(canLOPIT, device = "tiff", path = "output", dpi = 300, 
         file = paste(str_split(x, " / ")[[1]][1], 
                          "vs", 
                          str_split(x, " / ")[[1]][2],
                          "LOPIT_Human", myname,
                          sep = "_"))
}

LOPIT_Excel <- function(x){
  require(dplyr)
  require(openxlsx)
  ## Fold Change Plots
  df <- all_comparisons[all_comparisons$Comparison..group1.group2. == x,]
  df$minuslogqval <- -1*log10(df$Qvalue)
  df$Color <- ifelse(df$AVG.Log2.Ratio >= 0.58 & df$Qvalue < myQval, "Red", "Grey")
  df$Color <- ifelse(df$AVG.Log2.Ratio <= -0.58 & df$Qvalue < myQval, "Blue", df$Color)
  df$Label <- NA
  u2os$inCan <- u2os$Accession %in% df$UniProtIds #1721 proteins
  u2osCan <- df %>%
    full_join(u2os, by = c("UniProtIds" = "Accession"))
  
  # Excel Workbook for each comparison
  canSummary <- u2osCan %>%
    mutate(identified = !is.na(u2osCan$inCan)) %>%
    filter(inCan == TRUE) %>%
    group_by(assignment, Color) %>% count() %>%
    pivot_wider(names_from = c(Color), values_from = n)
  
  #rename the columns
  names(canSummary) <- c("Assignment", "Down-Regulated","Not Significant", "Up-Regulated")
  canSummary <- canSummary %>%
    mutate(Total = `Down-Regulated`+ `Up-Regulated` + `Not Significant`)
  
  
  mysheet <- paste(str_split(x, " / ")[[1]][1], 
                   "vs.", 
                   str_split(x, " / ")[[1]][2],
                   "- q ≤", myQval)
  addWorksheet(wb, sheetName = mysheet)
  writeData(wb, sheet = mysheet, 
            x = c(paste(batch,
                        "-", str_split(mysheet, " -")[[1]][1], 
                        search_info,
                        sep = " ")), 
            startRow = 1)
  writeData(wb, sheet = mysheet, 
            x = proIDs, 
            startRow = 2)
  writeData(wb, sheet = mysheet, 
            x = paste(comma(nrow(df[df$Color != "Grey",])), 
                      "Significantly Altered Protein Groups with |log2(FC)| ≥ 0.58 & q-value ≤",
                      myQval),
            startRow = 3)
  writeData(wb, sheet = mysheet, x = canSummary, 
            startRow = 5)
  writeData(wb, sheet = mysheet, x = colSums(canSummary[,2:ncol(canSummary)]),
            startRow = 20, startCol = "B")
}

if(add_q0.05 == TRUE){
  myQval <- 0.05
  myname <- "0p05"
  sapply(unique(all_comparisons$Comparison..group1.group2.), volcano)
  if(organism == "mouse"){
    Mouse_LOPIT_ref()
    sapply(unique(all_comparisons$Comparison..group1.group2.), Mouse_LOPIT)
    wb <- createWorkbook()
    sapply(unique(all_comparisons$Comparison..group1.group2.), LOPIT_Excel)
    saveWorkbook(wb, 
                 paste(batch,"_Mouse_LOPIT",myname,".xlsx", sep = ""),
                 overwrite = TRUE)
  }
  if(organism == "human"){
    Human_LOPIT_ref()
    sapply(unique(all_comparisons$Comparison..group1.group2.), Human_LOPIT)
    wb <- createWorkbook()
    sapply(unique(all_comparisons$Comparison..group1.group2.), LOPIT_Excel)
    saveWorkbook(wb, 
                 paste(batch,"_Human_LOPIT",myname,".xlsx", sep = ""),
                 overwrite = TRUE)
  }
}
if(add_q0.01 == TRUE){
  myQval <- 0.01
  myname <- "0p01"
  sapply(unique(all_comparisons$Comparison..group1.group2.), volcano)
  if(organism == "mouse"){
    Mouse_LOPIT_ref()
    sapply(unique(all_comparisons$Comparison..group1.group2.), Mouse_LOPIT)
    wb <- createWorkbook()
    sapply(unique(all_comparisons$Comparison..group1.group2.), LOPIT_Excel)
    saveWorkbook(wb, 
                 paste(batch,"_Mouse_LOPIT",myname,".xlsx", sep = ""),
                 overwrite = TRUE)
  }
  if(organism == "human"){
    Human_LOPIT_ref()
    sapply(unique(all_comparisons$Comparison..group1.group2.), Human_LOPIT)
    wb <- createWorkbook()
    sapply(unique(all_comparisons$Comparison..group1.group2.), LOPIT_Excel)
    saveWorkbook(wb, 
                 paste(batch,"_Human_LOPIT",myname,".xlsx", sep = ""),
                 overwrite = TRUE)
  }
}
if(add_q0.001 == TRUE){
  myQval <- 0.001
  myname <- "0p001"
  sapply(unique(all_comparisons$Comparison..group1.group2.), volcano)
  if(organism == "mouse"){
    Mouse_LOPIT_ref()
    sapply(unique(all_comparisons$Comparison..group1.group2.), Mouse_LOPIT)
    wb <- createWorkbook()
    sapply(unique(all_comparisons$Comparison..group1.group2.), LOPIT_Excel)
    saveWorkbook(wb, 
                 paste(batch,"_Mouse_LOPIT",myname,".xlsx", sep = ""),
                 overwrite = TRUE)
  }
  if(organism == "human"){
    Human_LOPIT_ref()
    sapply(unique(all_comparisons$Comparison..group1.group2.), Human_LOPIT)
    wb <- createWorkbook()
    sapply(unique(all_comparisons$Comparison..group1.group2.), LOPIT_Excel)
    saveWorkbook(wb, 
                 paste(batch,"_Human_LOPIT",myname,".xlsx", sep = ""),
                 overwrite = TRUE)
  }
}
