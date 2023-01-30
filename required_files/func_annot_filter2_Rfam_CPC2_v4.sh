#This script needs:
#a starting gtf (in previous dir) for extracting filtered annotation at the end
#a starting fasta (in previous dir)  for extracting filtered annotation at the end
#prefix for final output files
#the output file from func_annot_filter1_BLAST.sh script, named as '"$transcript_IDs"_blast_filtered' (or could swap in for another transcript ID list by entering name in blast_filtered_IDs variable)
#rfam & CPC2 results files in current dir

##this script will take a set of transcript IDs,
# find if they had hits in rfam
# finalise the list of loci where none of their transcripts had hits
# extract the transcripts from these loci
# find the CPC2 results from these transcripts
# select only those transcripts deemed 'noncoding' from loci with 50% or less coding transcripts, and extract gtf lines and fasta sequences for these
##
analysis_dir="<ADD PATH HERE!>"
gtf="ST_v5.3.gtf"
final_output_prefix="ST_v8.4"
transcript_IDs=$( echo $gtf | sed 's/.gtf//')
fasta="ST_v5.3.fa"

blast_filtered_IDs=$( echo "$transcript_IDs"_blast_filtered )

#####Rfam filter ###############################

echo "####### Rfam filter #########" > Rfam_CPC2_filter_results_overview

cat *.tblout > rfam_results

#seperate rRNA, tRNA and other hits
grep rRNA rfam_results | awk {'print$4'} | sort | uniq > unique_"$transcript_IDs"_rRNA_IDs
grep tRNA rfam_results  | awk {'print$4'} | sort | uniq > unique_"$transcript_IDs"_tRNA_IDs
grep MSTRG rfam_results | grep -v rRNA | grep -v tRNA | awk {'print$4'} | sort | uniq > unique_"$transcript_IDs"_nontr_RNAhits_IDs
grep MSTRG rfam_results  | grep -v rRNA | grep -v tRNA > nontr_RNAhits_lines


#get IDs and loci for all with rfam hits
grep MSTRG rfam_results | awk {'print$4'} | sort | uniq | sort > rfam_hits_"$transcript_IDs"_IDs
echo "Transcripts with Rfam hits:" $( grep MSTRG rfam_results | awk {'print$4'} | sort | uniq | wc -l )
echo "Transcripts with Rfam hits:" $( grep MSTRG rfam_results | awk {'print$4'} | sort | uniq | wc -l ) >>Rfam_CPC2_filter_results_overview
cat rfam_hits_"$transcript_IDs"_IDs | cut -d"." -f1-2 | sort | uniq > rfam_hits_"$transcript_IDs"_loci
echo "Loci with Rfam hits:" $( cat rfam_hits_"$transcript_IDs"_IDs | cut -d"." -f1-2 | sort | uniq | wc -l )
echo "Loci with Rfam hits:" $( cat rfam_hits_"$transcript_IDs"_IDs | cut -d"." -f1-2 | sort | uniq | wc -l ) >>Rfam_CPC2_filter_results_overview

sort $blast_filtered_IDs > sorted_$blast_filtered_IDs
#get loci of isoforms without hits:
comm -23 sorted_$blast_filtered_IDs rfam_hits_"$transcript_IDs"_IDs | cut -d"." -f1-2 | sort | uniq > rfam_no_hits_"$transcript_IDs"_loci


#find loci which have no rfam hits (ie.no isoforms with hits)
#ie. checking loci which are unique to '' list
comm -13 rfam_hits_"$transcript_IDs"_loci rfam_no_hits_"$transcript_IDs"_loci > rfam_no_hits_"$transcript_IDs"_FINAL_loci
echo "Loci with no Rfam hits:" $( comm -13 rfam_hits_"$transcript_IDs"_loci rfam_no_hits_"$transcript_IDs"_loci | wc -l) >>Rfam_CPC2_filter_results_overview
echo "Loci with no Rfam hits:" $( comm -13 rfam_hits_"$transcript_IDs"_loci rfam_no_hits_"$transcript_IDs"_loci | wc -l)


#use unique loci IDs to retrieve isoform IDs for those which belong to loci with no hits:
while read p
do
  grep "$p\." $transcript_IDs
done <rfam_no_hits_"$transcript_IDs"_FINAL_loci > "$transcript_IDs"_blast_rfam_filtered &&

echo "Transcripts with no Rfam hits:" $( wc -l "$transcript_IDs"_blast_rfam_filtered ) >>Rfam_CPC2_filter_results_overview
echo "Transcripts with no Rfam hits:" $( wc -l "$transcript_IDs"_blast_rfam_filtered )

#### CPC2 filter ###############################
echo "####### CPC2 filter #########" >> Rfam_CPC2_filter_results_overview


cat *result_cpc2.txt > CPC2_results
#extract results for blast and rfam filtered isoforms only
while read p 
do
grep -wF "$p" CPC2_results 
done <"$transcript_IDs"_blast_rfam_filtered > blast_rfam_filtered_CPC2_results


awk '$7=="noncoding" {print $0}' blast_rfam_filtered_CPC2_results > noncoding
awk '$7=="coding" {print $0}' blast_rfam_filtered_CPC2_results > coding

#calculate counts of coding and noncoding transcripts per locus 
cat coding | cut -f1 | cut -d"." -f1-2 | sort | uniq -c | awk '{$1=$1};1' | sed 's/ /\t/' > coding_transcript_counts.tsv
cat noncoding | cut -f1 | cut -d"." -f1-2 | sort | uniq -c | awk '{$1=$1};1' | sed 's/ /\t/' > noncoding_transcript_counts.tsv

#run R code to calculate percent coding per locus and select loci which are equal to or below cutoff (default cutoff = 50)
Rscript $analysis_dir/percent_coding.R 
#selected loci output in "noncoding_loci.csv"

##collect transcript IDs for noncoding loci 
#get IDs of all noncoding transcripts
cut -f1 noncoding > all_noncoding_transcripts
#extract final transcripts using selected loci
while read p 
do
    grep "$p\."  all_noncoding_transcripts 
done<noncoding_loci.csv > "$final_output_prefix"_noncoding_IDs

#add final transcript and loci number to overview
echo "noncoding transcripts:" $( wc -l "$final_output_prefix"_noncoding_IDs) >> Rfam_CPC2_filter_results_overview
echo "noncoding transcripts:" $( wc -l "$final_output_prefix"_noncoding_IDs)
echo "noncoding loci:" $( wc -l noncoding_loci.csv ) >> Rfam_CPC2_filter_results_overview

#get fasta
while read p; do grep -Fw "$p" -A 1 ../$fasta; done <"$final_output_prefix"_noncoding_IDs > "$final_output_prefix".fa
#check number extracted
echo "number fasta sequences extracted:" $( grep '>' "$final_output_prefix".fa | wc -l)
echo "number fasta sequences extracted:" $( grep '>' "$final_output_prefix".fa | wc -l) >> Rfam_CPC2_filter_results_overview


#extract stringtie annotation
while read p
do
grep -F "$p\"" ../$gtf
done <"$final_output_prefix"_noncoding_IDs > "$final_output_prefix".gtf


#check number of IDs in final gtf
echo "transcripts in final gtf:" $(grep transcript "$final_output_prefix".gtf | grep -v exon | wc -l)
echo "transcripts in final gtf:" $(grep transcript "$final_output_prefix".gtf | grep -v exon | wc -l) >> Rfam_CPC2_filter_results_overview
#
mkdir "$final_output_prefix"_output_files
cp "$final_output_prefix"_noncoding_IDs "$final_output_prefix"_output_files
cp "$final_output_prefix"_loci "$final_output_prefix"_output_files
cp "$final_output_prefix".gtf "$final_output_prefix"_output_files
cp "$final_output_prefix".fa "$final_output_prefix"_output_files
