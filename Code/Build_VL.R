suppressMessages({
  library(GenomicRanges)
  library(GenomicFeatures)
  library(rtracklayer)
  library(data.table)
  library(tidyr)
})

# 获取命令行参数
args <- commandArgs(trailingOnly = TRUE)
if(length(args) < 2){
  stop("Usage: Rscript script.R <TxDb_package_name> <output_file>")
}

txdb_pkg <- args[1]
output_file <- args[2]

# 加载 TxDb 包
library(txdb_pkg, character.only = TRUE)
txdb <- get(txdb_pkg)

# 提取信息
exons <- exonsBy(txdb, by = "tx")
introns <- intronsByTranscript(txdb)
transcripts <- transcripts(txdb)
meta <- mcols(transcripts)
tx_names <- meta$tx_name

# 合并 exon 和 intron
combined <- mapply(FUN = c, exons, introns, SIMPLIFY = FALSE)
names(combined) <- tx_names
combined <- endoapply(combined, sort)

# 定义函数
rebuild <- function(gr, index, output_file){
  df <- as.data.frame(gr)
  df <- df[c(1:3,5)]
  colnames(df) <- c("seq","start","end","strand")
  df_new <- df[rep(row.names(df), times = nrow(df)), ]
  df_new$start1 <- rep(df$start, times = nrow(df))
  df_new$end1 <- rep(df$end, each = nrow(df))
  df_new <- df_new[as.numeric(df_new$start1) < as.numeric(df_new$end1), ]
  df_new$transcript <- rep(index, nrow(df_new))
  df_new$name <- paste0(df_new$transcript, df_new$seq, ":", df_new$start1, "-", df_new$end1)
  df_new <- df_new[, c("name","seq","start1","end1","strand","transcript")]
  
  write.table(df_new, output_file, row.names = FALSE, col.names = FALSE, sep = ",", quote = FALSE, append = TRUE)
}

# 执行
sapply(seq_along(combined), function(index) rebuild(combined[[index]], names(combined)[index], output_file))
