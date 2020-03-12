
args <- commandArgs(trailingOnly = TRUE)
rds_fname <- args[1]

df <- readRDS(rds_fname)
sel_df <- df$reads[(df$reads$seqnames == 'chr1') & (df$reads$pos > 10e6) & (df$reads$pos < 11e6)]
write.table(sel_df, file="", sep = "\t", col.names = TRUE, row.names = FALSE)