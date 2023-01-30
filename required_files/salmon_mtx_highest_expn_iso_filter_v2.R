library(plyr)
library(tidyr)

#read in salmon counts matrix
counts = read.csv("salmon.isoform.counts.matrix_subset", header = TRUE, sep="\t") #read in file with sample names, reps and conditions
row.names(counts) <- counts[,1]
counts <- counts[,-1]

#add row sum as new column
counts[,(ncol(counts)+1)] <- rowSums(counts)
colnames(counts)[ncol(counts)] <- "rowsum"

#order by rownames
counts <- counts[ order(row.names(counts)), ]

#add sample names as a column
counts[,(ncol(counts)+1)] <- row.names(counts) 
colnames(counts)[ncol(counts)] <- "ID"
counts_sum <- counts[,(ncol(counts)-1):ncol(counts)]

counts_sum2 <- separate(data = counts_sum, col = ID, into = c("ID","Locus", "Isoform"), sep = "\\.") 

#use ddply to apply function to locus subsets 
#function takes rows, orders by decreasing rowsum and takes 1st row (ie. isoform with highest rowsum!)
filtered_counts <- ddply(counts_sum2,"Locus",function(counts_sum2){
  counts_sum2[order(counts_sum2$rowsum,decreasing=TRUE),][1,]
})

#merge id locus and isoform back together
filtered_IDs <- unite(filtered_counts[2:4],"ID", sep = ".", remove = TRUE, na.rm = FALSE)

#write IDs to file
write.table(filtered_IDs,"count_filtered_IDs.csv",sep=",",  col.names=FALSE,row.names = FALSE,quote = FALSE)


print(paste(nrow(counts) ,"features input", sep=" "))
print(paste(nrow(unique(counts_sum2["Locus"])),"loci input", sep=" "))
print(paste(nrow(filtered_IDs),"features output", sep=" "))
  