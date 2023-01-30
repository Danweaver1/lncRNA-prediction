gene_TPMs <- read.table(file = "gene_TPM_all_samples.tsv",sep = '\t', header = TRUE,row.names=1)
TPM_cutoff <- 0.5
to_keep_full <- vector()
#keep features which have TPM of 0.5 or more in all 3 reps of at least one sample
for (i in seq(1,ncol(gene_TPMs),by=3)){
  to_keep <- which(rowSums(gene_TPMs[,i:(i+2)]>=TPM_cutoff) >= 2)
  to_keep_full <- c(to_keep_full,to_keep)
}
final_to_keep <- unique(to_keep_full)
#extract TPM data for features to keep
filtered_data <- gene_TPMs[final_to_keep, ]
#write to file
write.csv(filtered_data,"filtered_TPM.tsv" , row.names=T)

print(paste(nrow(gene_TPMs) ,"features input", sep=" "))
print(paste(length(final_to_keep),"features with TPM of",TPM_cutoff,"or more in 2 reps of a sample", sep=" "))
print(paste(nrow(filtered_data),"features output", sep=" "))
