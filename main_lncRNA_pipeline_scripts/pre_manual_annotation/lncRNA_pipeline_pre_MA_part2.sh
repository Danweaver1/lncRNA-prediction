########## converts htseq counts to TPM values, then performs TPM filter #################################
#This part of the pipeline follows functional annotation #2 of the lncRNAs

##script which can have TPM cutoff value changed:
#htseq_TPM_transcript_filter.R

##scripts which need no changes:
#R_DESeq2_getcountsfile.R
#ML_counts_toTPM_v3.R

###data requirements:
##need samples.csv file for DE analysis in format:
#"sampleName","sampleFile","condition"
#"5Fluorocytosine_x05-1","5Fluorocytosine_x05-1_v7_7.counts","5Fluorocytosine_x05"
#"5Fluorocytosine_x05-2","5Fluorocytosine_x05-2_v7_7.counts","5Fluorocytosine_x05"
#"5Fluorocytosine_x05-3","5Fluorocytosine_x05-3_v7_7.counts","5Fluorocytosine_x05"
##htseq count files
##starting gtf and fasta files for corresponding annotation

##set main analysis directory, containing scripts to be fetched
analysis_dir="<ADD PATH HERE!>"

##set starting dir, should contain fasta & gtf file of corresponding annotation
starting_dir="<AND HERE!>"
fasta="formatted_ST_v8.5_merged.fa"
gtf="ST_v8.5_merged.gtf"

##prefix for final output files
final_output_prefix="ST_v8.5b"

##set directory for htseq counts files
htseq_dir=""
cd $starting_dir
mkdir TPM_filter_2
cd TPM_filter_2

##get lengths of fasta sequences
#code from: https://www.danielecook.com/generate-fasta-sequence-lengths/
awk '$0 ~ ">" {if (NR > 1) {print c;} c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }' ../"$fasta" > lengths

##copy htseq count files to current dir 
cp $htseq_dir/*.counts .

##run DEseq2 to get counts matrix
cp $analysis_dir/R_DESeq2_getcountsfile.R .
cp $analysis_dir/v8_5_samples.csv .
Rscript R_DESeq2_getcountsfile.R

##convert counts matrix to TPM 
cp $analysis_dir/ML_counts_toTPM_v3.R .
Rscript ML_counts_toTPM_v3.R
#Output: converted_TPM.csv

####FILTER STEP ##############
#### Uses TPM matrix as input and filters transcripts based on TPM cutoff in 2 reps of at least one sample

#get IDs from latest gtf
grep -o MSTRG\\.[0-9]*\\.[0-9]* ../$gtf | sort | uniq | sort | uniq | sort > "$gtf"_IDs

##create reduced matrix of counts for only those transcripts in current annotation version
head -1 converted_TPM.csv > converted_TPM_subset.csv

while read p
do
    grep -Fw $p converted_TPM.csv >> converted_TPM_subset.csv
done <"$gtf"_IDs

#### filter by TPM cutoff in 2 reps of at least one sample ####
##get R filter script & run
cp $analysis_dir/htseq_TPM_transcript_filter.R .
Rscript htseq_TPM_transcript_filter.R

##extract final annotation
#collapsed isoform ID list file
transcript_IDs="filtered_IDs.csv"

#check number of input IDs
echo "collapsed transcripts to fetch:" $( wc -l $transcript_IDs )
#get fasta
while read p; do grep -Fw "$p" -A 1 ../$fasta; done <$transcript_IDs > "$final_output_prefix".fa
#check number extracted
echo "number fasta sequences extracted:" $( grep '>' "$final_output_prefix".fa | wc -l)

#extract stringtie annotation
while read p
do
grep -F "$p\"" ../$gtf
done <$transcript_IDs > "$final_output_prefix".gtf

#check gtf records extracted
echo "number gtf records extracted:" $(grep -o MSTRG\\.[0-9]*\\.[0-9]* "$final_output_prefix".gtf | sort | uniq | sort | uniq | sort | wc -l)
