#code adapted from M. Love to convert counts to TPM: https://support.bioconductor.org/p/91218/
fasta.length <- read.table(file = "lengths",sep = '\t', header = FALSE,row.names=1)
#convert to kb
fasta.length[2] <- fasta.length/1000
fasta.length_v2 <- fasta.length[order(row.names(fasta.length)), ]
colnames(fasta.length_v2) <- c("bp","kb")

count <- read.table(file = "full_cnts.csv",sep = ',', header = TRUE,row.names=1)
count_v2  <- count[ order(row.names(count)), ]

x  <- count[1]/fasta.length_v2$kb
for (i in seq(2,length(count))){
  x[i]  <- count[i]/fasta.length_v2$kb
}

tpm.mat <- t( t(x) * 1e6 / colSums(x) )

write.csv(tpm.mat, "converted_TPM.csv")
