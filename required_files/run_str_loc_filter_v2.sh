# dir =
#necessary files: an annotation to filter against, Aspergillus_fumigatusa1163.ASM15014v1.dna.nonchromosomal.fa
# a master annotation with lines only (eg. xui_lines)
#and starting gtf(gtf in previous dir)
starting_gtf="ST_v5.2_xui.gtf"
filtered_output_prefix="ST_v5.3"

starting_prefix=$( echo $starting_gtf | sed 's/.gtf//')


annotation="~/A1163_master_v1.combined.gtf"

#if ID list exists, extract the gtf lines for the subset
grep -v exon ../$starting_gtf > lines_$starting_gtf



#check for annotations which start at 20 (they will change to 0 and be removed by awk)
if [ $(awk -F '\t' '$4==20 {print $0}' lines_$starting_gtf | wc -l ) -ne 0 ]
then
  echo "need to manually change feature with coordinate value equal to cutoff"
  awk -F '\t' '$4==20 {print $0}' lines_$starting_gtf | wc -l
  exit
fi


#extend coords by 20 each side
awk -F '\t' '$4-=20' lines_$starting_gtf | sed -e 's/ /\t/g' | awk '$5+=20' | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9" "$10" "$11" "$12" "$13" "$14" "$15" "$16" "$17" "$18" "$19" "$20" "$21" "$22}' > "$starting_gtf"_ex20.gtf

#remove any neg values
sed 's/-[0-9][0-9]/0/' "$starting_gtf"_ex20.gtf | sed 's/-[0-9]/0/' | sed 's/\t0\t/\t1\t/' > "$starting_prefix"_ex20_mod.gtf


#intersect 20 nt (stranded) against annotation
bedtools intersect -s -v -a "$starting_prefix"_ex20_mod.gtf  -b $annotation > "$filtered_output_prefix".ex_out
echo "Transcripts input:"$(cat "$starting_prefix"_ex20_mod.gtf | wc -l )
echo "How many loci input?" $(cat "$starting_prefix"_ex20_mod.gtf | cut -f9 | cut -d" " -f4 | sed 's/"//g' | sed 's/;//' | sort | uniq | wc -l )
echo "Transcripts remaining after filter:" $(bedtools intersect -s -v -a "$starting_prefix"_ex20_mod.gtf  -b $annotation | wc -l )

echo "Transcripts input:"$(cat "$starting_prefix"_ex20_mod.gtf | wc -l ) > results_overview
echo "How many loci input?" $(cat "$starting_prefix"_ex20_mod.gtf | cut -f9 | cut -d" " -f4 | sed 's/"//g' | sed 's/;//' | sort | uniq | wc -l ) >> results_overview
echo "Transcripts remaining after filter:" $(bedtools intersect -s -v -a "$starting_prefix"_ex20_mod.gtf  -b $annotation | wc -l ) >> results_overview
#####

echo "how many IDs output?" $( cat "$filtered_output_prefix".ex_out | cut -f9 | cut -d" " -f2 | sed 's/"//g' | sed 's/;//' | sort | uniq | wc -l )
echo "how many IDs output?" $( cat "$filtered_output_prefix".ex_out | cut -f9 | cut -d" " -f2 | sed 's/"//g' | sed 's/;//' | sort | uniq | wc -l ) >> results_overview

IDs_output=$(cat "$filtered_output_prefix".ex_out | cut -f9 | cut -d" " -f2 | sed 's/"//g' | sed 's/;//' | sort | uniq | wc -l )
#

echo "How many loci output?" $(cat "$filtered_output_prefix".ex_out | cut -f9 | cut -d" " -f4 | sed 's/"//g' | sed 's/;//' | sort | uniq | wc -l )

echo "How many loci output?" $(cat "$filtered_output_prefix".ex_out | cut -f9 | cut -d" " -f4 | sed 's/"//g' | sed 's/;//' | sort | uniq | wc -l ) >> results_overview

#get IDs
cat "$filtered_output_prefix".ex_out | cut -f9 | cut -d" " -f2 | sed 's/"//g' | sed 's/;//' | sort >"$filtered_output_prefix"_IDs

#GET ORIGINAL ANNOTATION (without nt extension)
while read p ; do  grep -F "$p\"" ../$starting_gtf ; done <"$filtered_output_prefix"_IDs > "$filtered_output_prefix".gtf
#check extracted all

if [ $(cat "$filtered_output_prefix".gtf  | cut -f9 | cut -d" " -f2 | sed 's/"//g' | sed 's/;//' | sort | uniq | wc -l) = $IDs_output ]
then
  echo "All extracted!"
  #extract fasta
  gffread -w "$filtered_output_prefix".fa -g Aspergillus_fumigatusa1163.ASM15014v1.dna.nonchromosomal.fa "$filtered_output_prefix".gtf  &&
  #send loci IDs to file
  grep '>' "$filtered_output_prefix".fa | cut -d" " -f1 | sed 's/>//' | sed 's/.rna1//' | cut -d"." -f1-2 | sort | uniq > unique_lociIDs

  #how many isoforms per loci?
  grep '>' "$filtered_output_prefix".fa | cut -d" " -f1 | sed 's/>//' | sed 's/.rna1//' | cut -d"." -f1-2 | sort | uniq -c | sort -n > isoforms_perloci
  #count loci with 1,2,3 or over 3 isoforms
  echo "Loci with a single isoform:" $(awk '$1==1 {print $0}' isoforms_perloci | wc -l )
  echo "Loci with two isoforms:" $(awk '$1==2 {print $0}' isoforms_perloci | wc -l )
  echo "Loci with three isoforms:" $(awk '$1==3 {print $0}' isoforms_perloci | wc -l )
  echo "Loci with four or more isoforms:" $(awk '$1>3 {print $0}' isoforms_perloci | wc -l)
  echo "Loci with a single isoform:" $(awk '$1==1 {print $0}' isoforms_perloci | wc -l ) >> results_overview
  echo "Loci with two isoforms:" $(awk '$1==2 {print $0}' isoforms_perloci | wc -l ) >> results_overview
  echo "Loci with three isoforms:" $(awk '$1==3 {print $0}' isoforms_perloci | wc -l ) >> results_overview
  echo "Loci with four or more isoforms:" $(awk '$1>3 {print $0}' isoforms_perloci | wc -l) >> results_overview
else
  echo "Extraction error - incorrect number of IDs"
fi
#
