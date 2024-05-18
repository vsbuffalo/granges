library(plyranges)
library(rtracklayer)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
  stop("Usage: Rscript script.R <bed_file_a> <bed_file_b> <output_file>")
}

bed_file_a <- args[1]
bed_file_b <- args[2]
output_file <- args[3]

# Load data into memory
granges_a <- read_bed(bed_file_a)
granges_b <- read_bed(bed_file_b)

# Perform the map operation
mapped_granges <- granges_a %>%
  join_overlap_inner(granges_b) %>%
  group_by(seqnames, start, end, strand) %>%
  summarise(
    min = min(score, na.rm = TRUE),
    max = max(score, na.rm = TRUE),
    mean = mean(score, na.rm = TRUE),
    sum = sum(score, na.rm = TRUE),
    median = median(score, na.rm = TRUE)
  )

mapped_df <- as.data.frame(mapped_granges)
write.table(mapped_df, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE)
