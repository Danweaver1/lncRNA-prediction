
#need starting transcript ID list, equivalent gtf (or an earlier one, its for extracting filtered annotation at the end) and prefix for final output files
#blast results need to be in current directory
gtf="ST_v5.3.gtf"
final_output_prefix="ST_v5.4"
transcript_IDs=$( echo $gtf | sed 's/.gtf//')

cat ../$gtf | grep -o MSTRG\\.[0-9]*\\.[0-9]*| sort | uniq | sort > $transcript_IDs
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

#extract all those without hits
grep " 0 hits" *.out -B 2 | grep MSTRG | cut -d" " -f3 > blastx_fungalnr_no_hits

#extract subset IDs for those without hits
while read p
do
  grep -woF $p blastx_fungalnr_no_hits
done < $transcript_IDs > "$transcript_IDs"_blastx_fungalnr_no_hits_IDs

#how many from subset without hits?
echo "Transcripts without blast hits:" $( cat "$transcript_IDs"_blastx_fungalnr_no_hits_IDs | wc -l )

#send overview to file
echo "Transcripts with blast hits:" $( cat "$transcript_IDs"_blastx_fungalnr_hit_IDs | wc -l )>>BLAST_filter_results_overview
echo "Transcripts without blast hits:" $( cat "$transcript_IDs"_blastx_fungalnr_no_hits_IDs | wc -l ) >>BLAST_filter_results_overview


####cross reference which LOCI are in both no hits and hits sets:
#loci with hits:
cat "$transcript_IDs"_blastx_fungalnr_hit_IDs | cut -d"." -f1-2 | sort | uniq | sort > "$transcript_IDs"_blastx_fungalnr_hit_loci

echo "Loci with blast hits:" $( cat "$transcript_IDs"_blastx_fungalnr_hit_IDs | cut -d"." -f1-2 | sort | uniq | sort | wc -l)
echo "Loci with blast hits:" $( cat "$transcript_IDs"_blastx_fungalnr_hit_IDs | cut -d"." -f1-2 | sort | uniq | sort | wc -l) >>BLAST_filter_results_overview


#loci of isoforms without hits:
cat "$transcript_IDs"_blastx_fungalnr_no_hits_IDs   | cut -d"." -f1-2 | sort | uniq | sort > "$transcript_IDs"_blastx_fungalnr_no_hits_loci
echo "Loci of transcripts with no blast hits:" $( cat "$transcript_IDs"_blastx_fungalnr_no_hits_IDs   | cut -d"." -f1-2 | sort | uniq | sort | wc -l)
echo "Loci of transcripts with no blast hits:" $( cat "$transcript_IDs"_blastx_fungalnr_no_hits_IDs   | cut -d"." -f1-2 | sort | uniq | sort | wc -l) >>BLAST_filter_results_overview



#loci where none of the isoforms have hits:
#ie. checking loci which are unique to '_fungalnr_nohits_loci' list
echo "Loci where all transcripts had no blast hits:" $( comm -13 "$transcript_IDs"_blastx_fungalnr_hit_loci "$transcript_IDs"_blastx_fungalnr_no_hits_loci | wc -l)
echo "Loci where all transcripts had no blast hits:" $( comm -13 "$transcript_IDs"_blastx_fungalnr_hit_loci "$transcript_IDs"_blastx_fungalnr_no_hits_loci | wc -l) >>BLAST_filter_results_overview
comm -13 "$transcript_IDs"_blastx_fungalnr_hit_loci "$transcript_IDs"_blastx_fungalnr_no_hits_loci > "$transcript_IDs"_blastx_fungalnr_FINAL_no_hits_loci


#use unique loci IDs to retrieve isoform IDs for those which belong to loci with no hits:
while read p
do
  grep "$p\."  $transcript_IDs
done <"$transcript_IDs"_blastx_fungalnr_FINAL_no_hits_loci > "$transcript_IDs"_blast_filtered &&

echo "BLAST filtered transcript IDs extracted:" $( wc -l "$transcript_IDs"_blast_filtered )
echo "BLAST filtered transcript IDs extracted:" $( wc -l "$transcript_IDs"_blast_filtered )  >>BLAST_filter_results_overview
