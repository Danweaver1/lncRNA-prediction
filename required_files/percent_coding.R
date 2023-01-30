#calculates percent coding transcripts per loci 
  
noncoding_counts <- read.table(file = "noncoding_transcript_counts.tsv",sep = '\t', header = FALSE,row.names=2)
colnames(noncoding_counts) <- "noncoding"

coding_counts <- read.table(file = "coding_transcript_counts.tsv",sep = '\t', header = FALSE,row.names=2)
colnames(coding_counts) <- "coding"

#join counts
merged_cnts <- merge(noncoding_counts,coding_counts,by="row.names", all=TRUE)
#replace NA with zero
merged_cnts[is.na(merged_cnts)] <- 0
#calculate total transcript counts per locus
merged_cnts[4] <- merged_cnts[2] + merged_cnts[3]
colnames(merged_cnts)[4] <- "total"

#calculate percent coding per locus
merged_cnts[5] <- (merged_cnts$coding / merged_cnts$total)*100
colnames(merged_cnts)[5] <- "percent_coding"

write.csv(merged_cnts, "percent_coding.csv")

##filter by percent coding
cutoff=50
suppressMessages(library(dplyr))
#find loci where all transcripts are coding
coding_loci <- merged_cnts %>% filter(percent_coding == 100)

write.table(coding_loci$Row.names, "coding_loci.csv",quote=FALSE,row.names = FALSE,col.names = FALSE)

#find putative coding loci - those which have percent coding above cutoff
putative_coding_loci <- merged_cnts  %>% filter(percent_coding < 100) %>% filter(percent_coding > cutoff)

write.table(putative_coding_loci$Row.names, "putative_coding_loci.csv",quote=FALSE,row.names = FALSE,col.names = FALSE)

#find putative noncoding loci - those which have percent coding equal to or below cutoff 
noncoding_loci <- merged_cnts %>% filter(percent_coding <= cutoff)

write.table(noncoding_loci$Row.names, "noncoding_loci.csv",quote=FALSE,row.names = FALSE,col.names = FALSE)

#print results
print(paste(nrow(merged_cnts) ,"features input", sep=" "))
print(paste(nrow(coding_loci) ,"loci where all transcripts are deemed coding", sep=" "))
print(paste(nrow(putative_coding_loci)," loci where over ",cutoff,"% of transcripts are deemed coding", sep=""))
print(paste(nrow(noncoding_loci),"noncoding loci output", sep=" "))
