suppressMessages({
  library(optparse)
  library(tidyr)
  library(dplyr)
  library(GenomicFeatures)
  library(GenomicRanges)
})

# ------------------------------
# Define command line options
# ------------------------------
option_list <- list(
  make_option(c("-c", "--circRNA_file"), type="character",
              help="Path to input circRNA file (txt)"),
  make_option(c("-p", "--mod_pvalue"), type="character",
              help="Path to modification p-value file (tsv)"),
  make_option(c("-l", "--mod_likelihood"), type="character",
              help="Path to modification likelihood file (tsv)"),
  make_option(c("-o", "--output_file"), type="character",
              help="Path to output CSV file for filtered methylation levels"),
  make_option(c("--min_depth"), type="integer", default=30,
              help="Minimum read depth to filter circRNAs [default: %default]"),
  make_option(c("--min_count"), type="integer", default=5,
              help="Minimum total count to keep a site [default: %default]")
)

opt <- parse_args(OptionParser(option_list=option_list))

# ------------------------------
# Load input data
# ------------------------------
circRNA <- read.table(opt$circRNA_file, header = FALSE)
circRNA <- unique(circRNA)

# Separate the V2 column into 'trans' and 'pos'
circRNA <- separate(circRNA, V2, c("trans", "pos"), sep = "(?<=\\d)(?=chr)", remove = FALSE)

# Filter circRNA by read depth
depth <- as.data.frame(table(circRNA$pos))
depth_filtered <- depth[depth$Freq >= opt$min_depth, ]
circRNA_filtered <- subset(circRNA, pos %in% depth_filtered$Var1)

# Load modification data
mod_p <- read.table(opt$mod_pvalue, header = FALSE)
mod_l <- read.table(opt$mod_likelihood, header = FALSE)

# Merge data
merged_df <- cbind(mod_p, mod_l)
merged_df <- merged_df[, -8:-13]
merged_df$V9 <- ifelse(merged_df$V7 == merged_df$V7.1, merged_df$V7, 0)

# Select modification sites
mod <- merged_df[merged_df$V9 == 1, ]

# Merge with filtered circRNA
merged <- merge(mod, circRNA_filtered, by.x = "V4", by.y = "V1")

# Filter records where modification sites are within circRNA regions
filtered_df <- merged[merged$V1 == merged$V3.y, ]
filtered_df <- filtered_df[filtered_df$V3.x == filtered_df$V6.y, ]
filtered_mod <- filtered_df[filtered_df$V2.x > filtered_df$V4.y & filtered_df$V2.x < filtered_df$V5.y, ]

# Count unique modification sites
unique_mod <- filtered_mod %>% count(V1, V2.x, V3.x)

# Filter reads
filtered_reads <- merged_df %>%
  semi_join(unique_mod, by = c("V1" = "V1", "V2" = "V2.x"))

# Calculate methylation levels
methylation_levels_mod <- filtered_reads %>%
  group_by(V1, V2, V3, V9) %>%
  summarise(Count = n(), .groups = "drop") %>%
  pivot_wider(names_from = V9, values_from = Count, values_fill = 0,
              names_prefix = "count_") %>%
  mutate(Total_Count = count_0 + count_1,
         freq_0 = count_0 / Total_Count,
         freq_1 = count_1 / Total_Count) %>%
  filter(Total_Count >= opt$min_count,
         count_1 != 0)

# Save output
write.table(methylation_levels_mod, file = opt$output_file,
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = ",")
