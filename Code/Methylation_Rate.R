# ---- load packages ----
suppressPackageStartupMessages({
  library(optparse)
  library(tidyr)
  library(dplyr)
})

# ---- options ----
option_list <- list(
  make_option(c("-p", "--mod_p"), type = "character", help = "Path to mod_p file (tsv)"),
  make_option(c("-l", "--mod_l"), type = "character", help = "Path to mod_l file (tsv)"),
  make_option(c("--pos"),   type = "character", help = "Path to pos file (csv)"),
  make_option(c("-o", "--output"), type = "character", help = "Path to output file"),
  make_option(c("-m", "--min_count"), type = "integer", default = 5,
              help = "Minimum total count")
)

opt <- parse_args(OptionParser(option_list = option_list))

# ---- read data ----
mod_p <- read.table(opt$mod_p, header = FALSE)
mod_l <- read.table(opt$mod_l, header = FALSE)
pos   <- read.table(opt$pos, header = FALSE, sep = ",")

colnames(pos) <- c("id","chr","start","end","strand")

# ---- merge ----
merged <- cbind(mod_p, mod_l)
merged <- merged[, -8:-13]
merged$V9 <- ifelse(merged$V7 == merged$V7.1, merged$V7, 0)
colnames(merged)[3] <- "strand"
merged <- subset(merged, strand == "+")

merged <- merged %>%
  separate(V1, into = c("chr", "range"), sep = ":", remove = FALSE) %>%
  separate(range, into = c("start", "end"), sep = "-", convert = TRUE) %>%
  mutate(width = end - start + 1)

merged_df <- merge(
  merged, pos[, c("chr", "start", "end", "strand")],
  by = c("chr", "start", "end"),
  all.x = TRUE
)

merged_df <- merged_df %>%
  mutate(
    mapped_pos = case_when(
      V2 <= width & strand.y == "+" ~ V2 + start - 1,
      V2 <= width & strand.y == "-" ~ end - V2 + 1,
      V2 >  width & strand.y == "+" ~ V2 - width + start - 1,
      V2 >  width & strand.y == "-" ~ end - (V2 - width) + 1
    )
  )

# ---- filter ----
mod_df <- merged_df[merged_df$V9 == 1, ]

mod_unique <- mod_df %>% count(chr, mapped_pos, strand.y)

filtered_reads <- merged_df %>%
  semi_join(mod_unique, by = c("chr" = "chr", "mapped_pos" = "mapped_pos"))

mod_methylation_levels <- filtered_reads %>%
  group_by(chr, mapped_pos, strand.y, V9) %>%
  summarise(Count = n(), .groups = "drop") %>%
  pivot_wider(names_from = V9, values_from = Count, values_fill = 0,
              names_prefix = "count_") %>%
  mutate(
    Total_Count = count_0 + count_1,
    freq_0 = count_0 / Total_Count,
    freq_1 = count_1 / Total_Count
  ) %>%
  filter(
    Total_Count >= opt$min_count,
    count_1 != 0
  )

# ---- rename columns ----
colnames(mod_methylation_levels)[2] <- "pos"
colnames(mod_methylation_levels)[3] <- "strand"

# ---- write output ----
write.table(
  mod_methylation_levels,
  file = opt$output,
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE,
  sep = ","
)
