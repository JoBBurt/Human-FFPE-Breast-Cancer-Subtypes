#1.0 LIBRARIES ----
  
library(openxlsx)
library(tidyverse)
library(scales)
library(stringr)

#2.0 User Inputs & Modifications ----
## Set the working directory
setwd("~/Desktop/CRUK_Storming_Cancer/Breast_FFPE/JB8-14_2rep")
## Set your output file
myoutput <- "output/Output_Data/Breast_FFPE_SASP_2rep_panHuman_23_0105_V01.xlsx"
## modify with your batch name
batch <- "JB8-JB14"
## modify with your search information
search_info <- " from 7 biological replicates with 2 technical replicates per biological replicate searched with the panHuman Library"
comparison_info <- "Breast Carcinoma Subtypes vs. Disease Free after Search with PanHuman Library"
## Which q-value filter would you like to include in the spreadsheet?
add_q0.05 <- FALSE
add_q0.01 <- FALSE
add_q0.001 <- TRUE
# 3.0  User Input GET DATA ----
protein <- read.csv("data/22_1018_JB8-JB14_2rep_panHuman_v01_Report_Birgit_Protein Quant_Pivot (Pivot).csv",
                    stringsAsFactors = F)
all_comparisons <- read.csv("output/Output_Data/SASP_q0p001.csv",
                            stringsAsFactors = F)
# 4.0 Format Data ----
protein <- as_tibble(protein)
colorder <- c("PG.Genes", "PG.ProteinDescriptions","PG.Qvalue",
              "PG.ProteinNames",
              "PG.UniProtIds", "PG.BiologicalProcess",
              "PG.CellularComponent", "PG.MolecularFunction", names(protein[,9:ncol(protein)]))
protein <- protein[,colorder]

all_comparisons <- as_tibble(all_comparisons)
colorder <- c("Comparison..group1.group2.", "Genes", "ProteinDescriptions", 
              "ProteinNames", "UniProtIds", "ProteinGroups", "Group", 
              "AVG.Log2.Ratio", "Qvalue", "Absolute.AVG.Log2.Ratio", 
              "Pvalue", "X..of.Ratios",
              names(all_comparisons[,13:ncol(all_comparisons)]))
all_comparisons <- all_comparisons[,colorder]

proIDs <- paste(comma(nrow(protein)), 
                "Protein Groups Identified with ≥ 2 Unique Peptides",
                sep = " ")
# 5.0 CREATE WORKBOOK ----

# * Initialize a workbook ----
wb <- createWorkbook()

# * Add a Worksheet ----
addWorksheet(wb, sheetName = "Identified Protein Groups")
addWorksheet(wb, sheetName = "All Comparisons")

# * Add Protein Identifications Data ----
writeData(wb, sheet = "Identified Protein Groups", 
          x = c(paste(batch,
                      "- Protein Group Identifications",
                      search_info,
                      sep = " ")), 
          startRow = 1)
writeData(wb, sheet = "Identified Protein Groups", 
          x = proIDs, 
          startRow = 2)
writeData(wb, sheet = "Identified Protein Groups", x = protein, 
          startRow = 4)

# * Add All Comparisons Data ---- 
writeData(wb, sheet = "All Comparisons", 
          x = c(paste(batch,
                      "- All Comparisons of",
                      comparison_info,
                      sep = " ")), 
          startRow = 1)
writeData(wb, sheet = "All Comparisons", 
          x = proIDs, 
          startRow = 2)
writeData(wb, sheet = "All Comparisons", 
          x = "Altered Protein Groups with No Filter",
          startRow = 3)
writeData(wb, sheet = "All Comparisons", x = all_comparisons, 
          startRow = 5)

# 5.0 Write Functions for q-value ----
write_all <- function(myQval){
  dq <- all_comparisons[all_comparisons$Absolute.AVG.Log2.Ratio >= 0.58 & 
                          all_comparisons$Qvalue <= myQval,]
  mySheet <- paste("Signif Altered - q ≤", myQval)
  addWorksheet(wb, sheetName = mySheet)
  writeData(wb, sheet = mySheet, 
            x = c(paste(batch,
                        "- All Comparisons of",
                        comparison_info,
                        sep = " ")), 
            startRow = 1)
  writeData(wb, sheet = mySheet, 
            x = proIDs, 
            startRow = 2)
  writeData(wb, sheet = mySheet, 
            x = paste("Significantly Altered Protein Groups with |log2(FC)| ≥ 0.58 & q-value ≤", myQval),
            startRow = 3)
  writeData(wb, sheet = mySheet, x = dq, 
            startRow = 5)
}

write_comparison <- function(x){
  df <- all_comparisons[all_comparisons$Absolute.AVG.Log2.Ratio >= 0.58 & 
                          all_comparisons$Qvalue <= myQval &
                          all_comparisons$Comparison..group1.group2. == x,]
  
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
            x = paste(comma(nrow(df)), 
                      "Significantly Altered Protein Groups with |log2(FC)| ≥ 0.58 & q-value ≤", myQval),
            startRow = 3)
  writeData(wb, sheet = mysheet, x = df, 
            startRow = 5)
}

# * Add q ≤ 0.05 data ----
if(add_q0.05 == TRUE){
  myQval <- 0.05
  write_all(myQval)
  sapply(unique(all_comparisons$Comparison..group1.group2.), write_comparison)
}

# * Add q ≤ 0.01 data ----
if(add_q0.01 == TRUE){
  myQval <- 0.01
  write_all(myQval)
  sapply(unique(all_comparisons$Comparison..group1.group2.), write_comparison)
}

# * Add q ≤ 0.01 data ----
if(add_q0.001 == TRUE){
  myQval <- 0.001
  write_all(myQval)
  sapply(unique(all_comparisons$Comparison..group1.group2.), write_comparison)
}

# * Save Workbook ----
saveWorkbook(wb, myoutput, overwrite = TRUE)

# * Open the Workbook ----
openXL(myoutput)
