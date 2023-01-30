#This script needs:
#a starting gtf for extracting filtered annotation at the end (in previous dir)
#a starting fasta for extracting filtered annotation at the end  (in previous dir)
#prefix for final output files
#the output file from func_annot_filter1_BLAST.sh script, named as '"$transcript_IDs"_blast_filtered' (or could swap in for another transcript ID list by entering name in blast_filtered_IDs variable)
#rfam & CPC2 results files in current dir (*.tblout and *result_cpc2.txt - NB if want to include several, they will automatically be combined)

##this script will take a set of transcript IDs,
# find those without rfam hits
# find the CPC2 results from these transcripts
# select only those transcripts deemed 'noncoding' and extract gtf lines and fasta sequences for these
##


gtf="ST_v11.1.gtf"
final_output_prefix="ST_v11.2"
transcript_IDs=$( echo $gtf | sed 's/.gtf//')
fasta="ST_v11.1.fa"
blast_filtered_IDs=$( echo "$transcript_IDs"_blast_filtered )

#####Rfam filter ###############################

echo "####### Rfam filter #########" > Rfam_CPC2_filter_results_overview


#extract rfam for starting set of transcripts
cat *.tblout > rfam.tblout


#seperate rRNA, tRNA and other hits
grep rRNA rfam.tblout | awk {'print$4'} | sort | uniq > unique_"$transcript_IDs"_rRNA_IDs
grep tRNA rfam.tblout  | awk {'print$4'} | sort | uniq > unique_"$transcript_IDs"_tRNA_IDs
grep MSTRG rfam.tblout | grep -v rRNA | grep -v tRNA | awk {'print$4'} | sort | uniq > unique_"$transcript_IDs"_nontr_RNAhits_IDs
grep MSTRG rfam.tblout  | grep -v rRNA | grep -v tRNA > nontr_RNAhits_lines


#get IDs and loci for all with rfam hits
grep MSTRG rfam.tblout | awk {'print$4'} | sort | uniq | sort > rfam_hits_"$transcript_IDs"_IDs
echo "Transcripts with Rfam hits:" $( grep MSTRG rfam.tblout | awk {'print$4'} | sort | uniq | wc -l )
echo "Transcripts with Rfam hits:" $( grep MSTRG rfam.tblout | awk {'print$4'} | sort | uniq | wc -l ) >>Rfam_CPC2_filter_results_overview

sort $blast_filtered_IDs > sorted_$blast_filtered_IDs
#get loci of isoforms without hits:
comm -23 sorted_$blast_filtered_IDs rfam_hits_"$transcript_IDs"_IDs | sort | uniq > blast_rfam_no_hits_"$transcript_IDs"

echo "Transcripts with no Rfam hits:" $( wc -l blast_rfam_no_hits_"$transcript_IDs" ) >>Rfam_CPC2_filter_results_overview
echo "Transcripts with no Rfam hits:" $( wc -l blast_rfam_no_hits_"$transcript_IDs" )

#### CPC2 filter ###############################
echo "####### CPC2 filter #########" >> Rfam_CPC2_filter_results_overview


cat *result_cpc2.txt > CPC2_results_file
#extract results for blast and rfam filtered isoforms only
while read p 
do
grep -wF "$p" CPC2_results_file
done <blast_rfam_no_hits_"$transcript_IDs" > blast_rfam_filtered_CPC2_results


awk '$7=="noncoding" {print $0}' blast_rfam_filtered_CPC2_results > noncoding
awk '$7=="coding" {print $0}' blast_rfam_filtered_CPC2_results > coding

#get IDs of noncoding
cut -f1 noncoding > "$final_output_prefix"_noncoding_IDs
echo "noncoding transcripts:" $( wc -l "$final_output_prefix"_noncoding_IDs) >> Rfam_CPC2_filter_results_overview
echo "noncoding transcripts:" $( wc -l "$final_output_prefix"_noncoding_IDs)
#get loci
cut -d"." -f1-2 "$final_output_prefix"_noncoding_IDs | sort | uniq >"$final_output_prefix"_loci
echo "noncoding loci:" $( wc -l "$final_output_prefix"_loci )
echo "noncoding loci:" $( wc -l "$final_output_prefix"_loci ) >> Rfam_CPC2_filter_results_overview

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
