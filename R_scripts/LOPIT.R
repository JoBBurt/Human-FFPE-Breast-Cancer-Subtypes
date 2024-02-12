df <- overlap %>%
  mutate(AVG.Log2.Ratio = rowMeans(overlap[,2:6])) %>%
  mutate(Absolute.AVG.Log2.Ratio = abs(AVG.Log2.Ratio))
df$Color <- ifelse(df$AVG.Log2.Ratio >= 0.58, "Red", "Grey")
df$Color <- ifelse(df$AVG.Log2.Ratio <= -0.58, "Blue", df$Color)
df$Label <- NA

meta <- bind_cols(all_comparisons$ProteinGroups, 
  all_comparisons$Genes, 
  all_comparisons$ProteinDescriptions, 
  all_comparisons$ProteinNames)
names(meta) <- c("UniProtIDs", "Genes", "ProteinDescriptions", "ProteinNames")
meta <- distinct(meta)
x<- inner_join(df, meta, by = c("Protein" = "Genes"))

u2os$inCan <- u2os$Accession %in% x$UniProtIDs #0 proteins
u2osCan <- x %>%
  full_join(u2os, by = c("UniProtIDs" = "Accession"))
