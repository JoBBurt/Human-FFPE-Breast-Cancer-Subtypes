# Load required packages
library(reticulate)

# Set Python path to include the cpdbWsdlClient folder
# Replace /path/to/cpdbWsdlClient with the actual folder path
Sys.setenv(PYTHONPATH = "/Users/jobburt/Documents/R_Toolkits/cpdbWsdlClient_py3/")

# Define a function to import all functions from all Python files in the folder
import_all_functions <- function(folder_path) {
  py_run_string(paste0("import os, sys, importlib",
                       "; sys.path.append('", folder_path, "')",
                       "; modules = [f[:-3] for f in os.listdir('", folder_path, "')",
                       " if f.endswith('.py')]",
                       "; globals().update({m: importlib.import_module(m) for m in modules})"))
}

# Import all functions from the cpdbWsdlClient folder
import_all_functions("/Users/jobburt/Documents/R_Toolkits/cpdbWsdlClient_py3/")




# Set your working directory to where your TSV files are stored
input_dir <- "~/Desktop/CRUK_Storming_Cancer/Breast_FFPE/JB8-14_2rep/output/ConsensusPathDB/results/"
output_dir <- "~/Desktop/CRUK_Storming_Cancer/Breast_FFPE/JB8-14_2rep/output/ConsensusPathDB/test_01/"

# create the output directory if it doesn't exist
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# create the dotplot output directory
dotplot_output_dir <- file.path(output_dir, "dotplots")
if (!dir.exists(dotplot_output_dir)) dir.create(dotplot_output_dir)

library(tidyverse)
library(ggplot2)

process_and_dotplot <- function(input_dir, output_dir, q.val = 0.01, t.level = 5, num = 10) {
  # Get a list of all TSV files in the input directory
  files <- list.files(input_dir, pattern = "\\.tsv$", full.names = TRUE)
  
  # Loop through each file and modify the column name and add the new columns
  for (file in files) {
    # Read in the TSV file
    data <- read.delim(file)
    
    # Rename the "term_name" column to "pathway"
    colnames(data)[colnames(data) == "term_name"] <- "pathway"
    
    # Calculate the new columns
    count <- nchar(data$members_input_overlap_geneids) - nchar(gsub(";", "", data$members_input_overlap_geneids)) + 1
    GeneRatio.non.corrected <- count / data$size
    GeneRatio <- count / data$effective_size
    
    # Add the new columns to the data frame
    data$count <- count
    data$GeneRatio.non.corrected <- GeneRatio.non.corrected
    data$GeneRatio <- GeneRatio
    
    # Write the modified data frame to a new TSV file
    output_file <- file.path(output_dir, paste0("modified_", basename(file)))
    write.csv(data, output_file, row.names = FALSE)
    
    # Generate dotplots for the modified file
    generate_dotplot <- function(data, file_name, path, output_path) {
      # filter data for the current term category
      data <- data %>% filter(term_category == path, term_level >= t.level)
      
      # check if there are any rows that match the filter condition
      if (sum(data$q.value < q.val) == 0) {
        print(paste("No terms with q-value <", q.val, "in", file_name, "for", path, "category"))
        return()
      }
      
      # select the top 10 terms by GeneRatio
      top10 <- data %>%
        filter(q.value < q.val) %>%
        arrange(desc(GeneRatio)) %>%
        top_n(num)
      
      g <- ggplot(top10, aes(x = GeneRatio, y = factor(pathway, levels = rev(pathway)), 
                             size = count, color = q.value)) +
        geom_point(alpha = 0.8) +
        scale_color_gradient(low = "red2", high = "mediumblue", space = "Lab", 
                             limit = c(min(data$q.value), max(data$q.value))) +
        scale_y_discrete(name = "") +
        scale_size(range = c(2, 6)) +
        theme_classic() +
        theme(
          axis.title.x = element_text(size = 16),
          axis.text.x = element_text(size = 14, angle = 90, vjust = 0.5, hjust = 1),
          axis.text.y = element_text(size = 10)
        )
      
      # save the plot to a file
      ggsave(file.path(output_path, paste0(file_name, "_", path, ".png")), 
             g, width = 8, height = 5, dpi = 300)
    }
    
    
    # create the output directory for the dotplots
    dotplot_output_dir <- file.path(output_dir, "dotplots", basename(file))
    if (!dir.exists(dotplot_output_dir)) dir.create(dotplot_output_dir, recursive = TRUE)
    
    # generate dotplots for each term category
    for (path in c("b", "m", "c")) {
      generate_dotplot(data, tools::file_path_sans_ext(basename(file)), path, dotplot_output_dir)
    }
  }
}


process_and_dotplot(input_dir, output_dir, t.level = 5, q.val = 0.05)
