
####begin this analysis in dir containing stringtie expression output and expression metric matrices ##############
#### in the dir above is the starting stringtie gtf: "stringtie_merged_joined.gtf"
analysis_dir="<ADD PATH HERE!>"
cd /home/mfbx9dw5/Dropbox/Multidrug_annotation/assembly_3_stringtie/stringtie_pipeline_v4/ST_expression_200527
#set current dir ( ie. dir below $analysis_dir )

#set genome annotation path and filename
annotation="~/A1163_master_v1.combined.gtf"
#set gtf version for pre functional annotation
gtf_version=""
#set gtf version for post functional annotation
gtf_v_after_functional_annot=""
#set genome fasta path
genome_fasta_dir="<AND HERE!>"
#set genome fasta filename
genome_fasta=""
#there needs to be an index file (.fai) for the above genome fasta in the same dir

#results for functional annotation filter:
blast_dir="<AND HERE!>"
rfam_dir="<AND HERE!>"
CPC2_dir="<AND HERE!>"
#results for Salmon quantification
salmon_dir="<AND HERE!>"
#additional scripts needed should be present in $analysis_dir. They are then copied to working dir at the time of use:
#following need variables changed before running:
###run_str_loc_filter_v2.sh
###func_annot_filter1_BLAST_v2.sh
###func_annot_filter2_Rfam_CPC2_v4.sh
###extract_collapsed_isoforms_v2.sh

###following scripts do not require any changes:
#salmon_mtx_highest_expn_iso_filter_v2.R
#abundance_estimates_to_matrix.pl
#expression_matrix_200527.sh
#stringtie_expression_matrix_DWmod.pl
#percent_coding.R


#### TPM filter #1 (gene level) #######################################################################################################
cp ../expression_matrix_200527.sh .
cp ../stringtie_expression_matrix_DWmod.pl .
./expression_matrix_200527.sh 

Rscript ../TPM_cutoff_filter.R
####Rscript code: ####
##gene_TPMs <- read.table(file = "gene_TPM_all_samples.tsv",sep = '\t', header = TRUE,row.names=1)
#TPM_cutoff <- 0.5
#to_keep_full <- vector()
##keep features which have TPM of 0.5 or more in all 3 reps of at least one sample
#for (i in seq(1,ncol(gene_TPMs),by=3)){
#  to_keep <- which(rowSums(gene_TPMs[,i:(i+2)]>=TPM_cutoff) >= 2)
#  to_keep_full <- c(to_keep_full,to_keep)
#}
#final_to_keep <- unique(to_keep_full)
##extract TPM data for features to keep
#filtered_data <- gene_TPMs[final_to_keep, ]
##write to file
#write.csv(filtered_data,"filtered_TPM.tsv" , row.names=T)

#print(paste(nrow(gene_TPMs) ,"features input", sep=" "))
#print(paste(length(final_to_keep),"features with TPM of",TPM_cutoff,"or more in 2 reps of a sample", sep=" "))
#print(paste(nrow(filtered_data),"features output", sep=" "))
#### Rscript code end ####


## Extract filtered MSTRG IDs
#take ID list from filtered TPM matrix
cut -d"," -f1 filtered_TPM.tsv | sed 's/"//g' | sort > "$gtf_version".1_IDs

#extract relevant lines for all isoforms of loci from stringtie gtf
while read p; do  grep -F "\"$p\"" ../stringtie_merged_joined.gtf; done<"$gtf_version".1_IDs > "$gtf_version".1.gtf &&
#check there are 6870 MSTRG loci
grep -o MSTRG\\.[0-9]* "$gtf_version".1.gtf | sort | uniq | sort | wc -l 

#how many isoforms?
grep -o MSTRG\\.[0-9]*\\.[0-9]* "$gtf_version".1.gtf | sort | uniq | sort | uniq | sort | wc -l 

#how many isoforms did we start with? 
grep -o MSTRG\\.[0-9]*\\.[0-9]* ../stringtie_merged_joined.gtf | sort | uniq | sort | uniq | sort | wc -l 

echo "Final gtf from this step = "$gtf_version".1.gtf"

#####remove annotated features #######################################################
mkdir remove_annotated_features

#fetch the 'master' A1163 annotation 
cp $annotation remove_annotated_features/

cd remove_annotated_features/

##run gffcompare
gffcompare -o "$gtf_version".1_gff_cmp -r A1163_master_v1.combined.gtf ../"$gtf_version".1.gtf


#check class codes (overall)
grep -o "class_code \"[a-z]\"" "$gtf_version".1_gff_cmp.annotated.gtf | sort | uniq -c | sort -n


#how many xui MSTRG isoforms?
cut -f9 "$gtf_version".1_gff_cmp.annotated.gtf | grep -v 'transcript_id \"transcript' | grep 'class_code "[xui]"' | cut -d";" -f1 | cut -d" " -f2 | sed 's/"//g' | sort  | wc -l


#how many loci?
cut -f9 "$gtf_version".1_gff_cmp.annotated.gtf | grep -v 'transcript_id \"transcript' | grep 'class_code "[xui]"' | cut -d";" -f1 | cut -d" " -f2 | sed 's/"//g' | cut -d"." -f1-2 | sort | uniq | sort | wc -l 

#extract xui isoform IDs
cut -f9 "$gtf_version".1_gff_cmp.annotated.gtf | grep -v 'transcript_id \"transcript' | grep 'class_code "[xui]"' | cut -d";" -f1 | cut -d" " -f2 | sed 's/"//g' | sort > xui_MSTRG_IDs


##For actual numbers of x,u or i MSTRGs being taken:

#how many antisense MSTRGs?
cut -f9 "$gtf_version".1_gff_cmp.annotated.gtf | grep -v 'transcript_id \"transcript' | grep 'class_code "[x]"' | cut -d";" -f1 | cut -d" " -f2 | sed 's/"//g' | sort | wc -l 
#5980
#how many intergenic MSTRGs?
cut -f9 "$gtf_version".1_gff_cmp.annotated.gtf | grep -v 'transcript_id \"transcript' | grep 'class_code "[u]"' | cut -d";" -f1 | cut -d" " -f2 | sed 's/"//g' | sort | wc -l 
#5369
#how many intronic MSTRGs?
cut -f9 "$gtf_version".1_gff_cmp.annotated.gtf | grep -v 'transcript_id \"transcript' | grep 'class_code "[i]"' | cut -d";" -f1 | cut -d" " -f2 | sed 's/"//g' | sort | wc -l 
#40


##Extract xui lncRNA at the isoform level: 

while read p ; do grep -F "$p\"" "$gtf_version".1_gff_cmp.annotated.gtf |grep -v 'transcript_id \"transcript'; done<xui_MSTRG_IDs >"$gtf_version".2_xui.gtf

#check 11389 were extracted
grep -o MSTRG\\.[0-9]*\\.[0-9]* "$gtf_version".2_xui.gtf | sort | uniq | sort | uniq | sort | wc -l

#final output: version .2


####Stranded 20nt location filter ####################################################################
mkdir stranded_location_filter
cd stranded_location_filter

cp "$genome_fasta_dir"/"$genome_fasta" .
cp "$genome_fasta_dir"/"$genome_fasta".fai .
cp ../../../run_str_loc_filter_v2.sh .
chmod +x run_str_loc_filter_v2.sh
./run_str_loc_filter_v2.sh

fasta_formatter -i "$gtf_version".3.fa -o formatted_"$gtf_version".3.fa 
mv formatted_"$gtf_version".3.fa "$gtf_version".3.fa 

#final output: version .3

#### Functional annotation ############################################################################

mkdir functional_annotation_filter
cd functional_annotation_filter

##copy in the results files
cp -r $blast_dir/*.out .
cp -r $rfam_dir/*.tblout .
cp -r $CPC2_dir/*cpc2.txt .

##get first filter script - BLASTx hits - and run
cp $analysis_dir/func_annot_filter1_BLAST_v2.sh .

chmod +x func_annot_filter1_BLAST_v2.sh
./func_annot_filter1_BLAST_v2.sh

##get second filter script - rfam hits & CPC2 coding potential - and run
cp $analysis_dir/func_annot_filter2_Rfam_CPC2_v4.sh .
chmod +x func_annot_filter2_Rfam_CPC2_v4.sh
./func_annot_filter2_Rfam_CPC2_v4.sh


##create summary of functional annotation data
mkdir func_annot1_summary
cd func_annot1_summary

while read p
do
  echo "|"
  echo $p &&
  if grep -Fqw $p ../blastx_fungalnr_strandedhits; then
      grep -Fw $p ../blastx_fungalnr_strandedhits | awk -v OFS=',' '{print $2,$14,$17,$20}'
  else
      echo "-,-,-,-"
  fi &&

  if grep -Fqw $p ../rfam_results; then
    grep -Fw $p ../rfam_results | awk -v OFS=',' '{print $2,$3,$18}'
  else
    echo "-,-,-"
  fi &&

  if grep -Fqw $p ../CPC2_results; then
    grep -Fw $p ../CPC2_results | awk -v OFS=',' '{print $2,$6,$7}'
  else
    echo "-,-,-"
  fi
done<../"$gtf_version".3 >combi_results


(echo ",ID,Blastx hit accession,alignment length,%identity,evalue,Rfam target name,accession,E-value,CPC peptide length,coding probability,label"; cat combi_results) | tr "|\n" "\n," | sed 's/,//' > combi_results_final

#final output: version .4

#### Isoform collapsing ####################################################

### get Salmon quantification results matrix
cd ../"$gtf_v_after_functional_annot".4_output_files
mkdir salmon
cd salmon
cp $salmon_dir/salmon.isoform.counts.matrix .

##create reduced matrix of counts for only those transcripts in current annotation version
head -1 salmon.isoform.counts.matrix > salmon.isoform.counts.matrix_subset

while read p
do
    grep -Fw $p salmon.isoform.counts.matrix >> salmon.isoform.counts.matrix_subset
done <../"$gtf_v_after_functional_annot".4_noncoding_IDs

##get R filter script & run
cp $analysis_dir/salmon_mtx_highest_expn_iso_filter_v2.R .
Rscript salmon_mtx_highest_expn_iso_filter_v2.R


##get script to extract gtf and fa for collapsed isoforms

cp $analysis_dir/extract_collapsed_isoforms_v2.sh .
chmod +x extract_collapsed_isoforms_v2.sh
./extract_collapsed_isoforms_v2.sh

#final output: version .5

