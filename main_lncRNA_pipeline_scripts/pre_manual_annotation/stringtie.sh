##de novo stringtie
for file in */*.bam
do
PREFIX=$( echo $file | cut -d"/" -f2 | sed 's/.bam//' )
Sample=$( echo $PREFIX | cut -d"-" -f1 )
stringtie $file -p 24 -l $PREFIX -o $Sample/$PREFIX/transcripts.gtf --rf 
done

##merge all sample assemblies
stringtie --merge --rf -p 8 -o stringtie_merged.gtf assembly_GTF_list.txt

##genome guided stringtie
for file in */*.bam
do
PREFIX=$( echo $file | cut -d"/" -f2 | sed 's/.bam//' )
Sample=$( echo $PREFIX | cut -d"-" -f1 )
stringtie $file -p 12 -G GG_stringtie/Aspergillus_fumigatusa1163.ASM15014v1.44.gff3 -l $PREFIX -o GG_stringtie/$Sample/$PREFIX/transcripts.gtf --rf 
done

#merge all genome guided sample assemblies
stringtie --merge --rf -p 8 -o stringtie_GG_merged.gtf -G Aspergillus_fumigatusa1163.ASM15014v1.44.gff3 assembly_GG_GTF_list.txt

#merge de novo & genome guided assemblies together
stringtie --merge --rf -p 8 -o stringtie_merged_joined.gtf -G Aspergillus_fumigatusa1163.ASM15014v1.44.gff3 list.txt
#list.txt contains:
#stringtie_GG_merged.gtf
#stringtie_merged.gtf


#compare de novo & genome guided merged stringtie transcripts with A1163 annotation
gffcompare -G -r Aspergillus_fumigatusa1163.ASM15014v1.44.gff3 stringtie_merged_joined.gtf 
