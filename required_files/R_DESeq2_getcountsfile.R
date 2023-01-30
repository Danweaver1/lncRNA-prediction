suppressMessages(library(DESeq2))

sample_details = read.csv(list.files(pattern="v8_5_samples.csv"), header = TRUE, sep=",") #read in file with sample names, reps and conditions
  
condition = sample_details$condition #create variable for conditions

  
sampleTable = data.frame(sampleName = sample_details$sampleName, fileName = sample_details$sampleFile, condition = condition) #create sample info table with filenames, samples names and conditions
write.csv(sampleTable, "sampleTable.csv") #create a file of sample Table
ddsHTSeq = DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, design = ~condition)
  
counts = counts(ddsHTSeq) #create variable of counts
write.csv(counts, "full_cnts.csv") #create a file of counts
  