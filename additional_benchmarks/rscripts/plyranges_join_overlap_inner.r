library(plyranges)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: Rscript script.R <bed_file_a> <bed_file_b> <output_file>")
}

bed_file_a <- args[1]
bed_file_b <- args[2]
output_file <- args[3]

# load data into mem
granges_a <- read_bed(bed_file_a)
granges_b <- read_bed(bed_file_b)

intersection_result <- join_overlap_inner(granges_a, granges_b)

write_bed(intersection_result, output_file)
