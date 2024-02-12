# Set your working directory to where your TSV files are stored
setwd("~/Desktop/CRUK_Storming_Cancer/Breast_FFPE/JB8-14_2rep/output/ConsensusPathDB/results/")

# Get a list of all TSV files in the directory
files <- list.files(pattern = "\\.tsv$")

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
  write.csv(data, paste0("modified_", file))
}

# generate dotplots

dotplot <- function(dir_path, output_dir, q.val = 0.01, t.level = 5, num = 10) {
  # define a helper function to generate the dotplot for a given term category
  generate_dotplot <- function(data, path, output_path) {
    top10 <- data %>%
      filter(q.value < q.val) %>%
      filter(term_category == path) %>%
      filter(term_level >= t.level) %>%
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
    ggsave(file.path(output_path, paste0(path, ".png")), g, width = 8, height = 5, dpi = 300)
  }
  
  # get a list of all ".tsv" files in the directory
  file_paths <- list.files(dir_path, pattern = "\\.tsv$", full.names = TRUE)
  
  # process each file
  for (file_path in file_paths) {
    # read in the file
    data <- read.delim(file_path, stringsAsFactors = FALSE)
    
    # create the output directory if it doesn't exist
    output_path <- file.path(output_dir, basename(file_path))
    if (!dir.exists(output_path)) dir.create(output_path)
    
    # generate and save the dotplots for each term category
    for (path in c("b", "m", "c")) {
      generate_dotplot(data, path, output_path)
    }
  }
}
