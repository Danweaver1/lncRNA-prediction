#need starting transcript ID list, equivalent gtf (or an earlier one, its for extracting filtered annotation at the end) and prefix for final output files
#blast results need to be in current directory (*.out file)
#modifications still functional annotation 2 include adding a search for IDs containing lettered isoforms as well as standard numbered MSTRG isoforms
gtf="ST_v11.1.gtf"
final_output_prefix="ST_v11.2"
transcript_IDs=$( echo $gtf | sed 's/.gtf//')

#search for numbered isoforms specifically - one number is required in isoform section of the ID
cat ../$gtf | grep -o MSTRG\\.[0-9]*\\.[0-9][0-9]*| sort | uniq | sort > num_"$transcript_IDs"
#search for lettered isoforms (manually reannotated)
cat ../$gtf | grep -o MSTRG\\.[0-9]*\\.[A-Z]| sort | uniq | sort > let_"$transcript_IDs"
#merge isoform IDs and remove intermediate ID files
cat let_"$transcript_IDs" num_"$transcript_IDs" | sort > "$transcript_IDs"
rm let_"$transcript_IDs"
rm num_"$transcript_IDs"

##BLAST filter #################################
echo "####### BLAST filter #########" > BLAST_filter_results_overview

#extract all hits (only taking top hit record)
grep ' [1-9][0-9]* hits' -A 1 *.out > blastx_fungalnr_strandedhits
#extract subset IDs from hits
while read p
do
  grep -woF $p blastx_fungalnr_strandedhits
done < $transcript_IDs | sort | uniq | sort > "$transcript_IDs"_blastx_fungalnr_hit_IDs

#how many from subset with hits?
echo "Transcripts with blast hits:" $( cat "$transcript_IDs"_blastx_fungalnr_hit_IDs | wc -l )
echo "Transcripts with blast hits:" $( cat "$transcript_IDs"_blastx_fungalnr_hit_IDs | wc -l ) >>BLAST_filter_results_overview

#extract all those without hits
grep " 0 hits" *.out -B 2 | grep MSTRG | cut -d" " -f3 > blastx_fungalnr_no_hits

#extract subset IDs for those without hits
while read p
do
  grep -woF $p blastx_fungalnr_no_hits
done < $transcript_IDs > "$transcript_IDs"_blast_filtered


echo "BLAST filtered transcript IDs extracted:" $( wc -l "$transcript_IDs"_blast_filtered )
echo "BLAST filtered transcript IDs extracted:" $( wc -l "$transcript_IDs"_blast_filtered )  >>BLAST_filter_results_overview
