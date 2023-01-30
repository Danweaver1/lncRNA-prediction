transcript_TPMs <- read.table(file = "salmon.isoform.TPM.not_cross_norm.matrix_subset",sep = '\t', header = TRUE,row.names=1)
TPM_cutoff <- 1
to_keep_full <- vector()
#keep features which have TPM of 0.5 or more in all 3 reps of at least one sample
for (i in seq(1,ncol(transcript_TPMs),by=3)){
  to_keep <- which(rowSums(transcript_TPMs[,i:(i+2)]>=TPM_cutoff) >= 2)
  to_keep_full <- c(to_keep_full,to_keep)
}
final_to_keep <- unique(to_keep_full)
#extract TPM data for features to keep
filtered_data <- transcript_TPMs[final_to_keep, ]
#write to file
write.csv(filtered_data,"filtered_transcripts_TPM.tsv" , row.names=T)

filtered_IDs <- rownames(filtered_data)
#write IDs to file
write.table(filtered_IDs,"filtered_IDs.csv",sep=",",  col.names=FALSE,row.names = FALSE,quote = FALSE)

print(paste(nrow(transcript_TPMs) ,"features input", sep=" "))
print(paste(length(final_to_keep),"features with TPM of",TPM_cutoff,"or more in 2 reps of a sample", sep=" "))
print(paste(nrow(filtered_data),"features output", sep=" "))
