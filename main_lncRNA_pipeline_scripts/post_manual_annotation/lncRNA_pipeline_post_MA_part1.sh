#### Functional annotation #3 ############################################################################
#This part of the pipeline follows full lncRNA manual reannotation effort 

#### main pipeline directory, containing relevant scripts to be copied ####
analysis_dir="<ADD PATH HERE!>"

####begin this analysis in dir below, containing genome fasta and index (.fai) and latest manually annotated gtf and fasta files ##############
cd <AND HERE!> ####

#there needs to be an index file (.fai) for the above genome fasta in the same dir
#results for functional annotation filter:
blast_dir="<AND HERE!>"
rfam_dir="<AND HERE!>"
CPC2_dir="<AND HERE!>"

#additional scripts needed should be present in the initial working dir. They are then copied to working dir at the time of use:
#following need variables changed before running:
###func_annot2_filter1_BLAST_v2.sh
###func_annot2_filter2_Rfam_CPC2_v2.sh


mkdir functional_annot3
cd functional_annot3

##copy in the results files
cp -r $blast_dir/*.out .
cp -r $rfam_dir/*.tblout .
cp -r $CPC2_dir/*cpc2.txt .

##get first filter script - BLASTx hits - and run
#filters by transcript hits only (no filtering at loci level)
cp $analysis_dir/func_annot3_filter1_BLAST_v1.sh .


chmod +x func_annot3_filter1_BLAST_v1.sh
./func_annot3_filter1_BLAST_v1.sh


##get second filter script - rfam hits and cpc2 coding potential - and run
#filters by transcript hits only (no filtering at loci level)
cp $analysis_dir/func_annot3_filter2_Rfam_CPC2_v1.sh .
chmod +x func_annot3_filter2_Rfam_CPC2_v1.sh
./func_annot3_filter2_Rfam_CPC2_v1.sh

