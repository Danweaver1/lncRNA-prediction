
#perl script to make matrix takes directories containing stringtie output

#get directory paths
dirs=$(ls -d */*/* | tr '\n' ',')
#run perl script on directories
./stringtie_expression_matrix_DWmod.pl --expression_metric=TPM --result_dirs="$dirs" --transcript_matrix_file=NOTtranscript_TPM_all_samples.tsv --gene_matrix_file=gene_TPM_all_samples.tsv
./stringtie_expression_matrix_DWmod.pl --expression_metric=FPKM --result_dirs="$dirs" --transcript_matrix_file=NOTtranscript_FPKM_all_samples.tsv --gene_matrix_file=gene_FPKM_all_samples.tsv
./stringtie_expression_matrix_DWmod.pl --expression_metric=Coverage --result_dirs="$dirs" --transcript_matrix_file=NOTtranscript_Coverage_all_samples.tsv --gene_matrix_file=gene_Coverage_all_samples.tsv

head *all_samples.tsv
