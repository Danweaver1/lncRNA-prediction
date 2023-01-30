#files to extract from
gtf="../ST_v8.4.gtf"
fasta="../ST_v8.4.fa"

#collapsed isoform ID list file
transcript_IDs="count_filtered_IDs.csv"
#prefix for final output files
final_output_prefix="ST_v8.5"


#check number of input IDs
echo "collapsed transcripts to fetch:" $( wc -l $transcript_IDs )
#get fasta
while read p; do grep -Fw "$p" -A 1 $fasta; done <$transcript_IDs > "$final_output_prefix".fa
#check number extracted
echo "number fasta sequences extracted:" $( grep '>' "$final_output_prefix".fa | wc -l)




#extract stringtie annotation
while read p
do
grep -F "$p\"" ../$gtf
done <$transcript_IDs > "$final_output_prefix".gtf

#check gtf records extracted
echo "number gtf records extracted:" $(grep -o MSTRG\\.[0-9]*\\.[0-9]* "$final_output_prefix".gtf | sort | uniq | sort | uniq | sort | wc -l)
