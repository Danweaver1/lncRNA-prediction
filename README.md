# lncRNA-prediction
Pipeline to identify long noncoding RNA (lncRNA) from transcriptomic data in Aspergillus fumigatus. The pipeline is split into two sections: pre-manual annotation and post-manual annotation. 

Broadly, the pre-manual annotation pipeline performs the following steps:
- Uses StringTie to assemble transcripts
- Filters transcripts based on known annotated features, and removes any features extremely close to genes (this step is necessary due to the poor annotation of UTRs in A. fumigatus - many transcripts close to annotated genes are actually unannotated UTRs)
- Performs functional annotation using BLAST, Rfam (INFERNAL) and Coding Potential Calculator 2 (CPC2) to identify candidate lncRNA.
- Collapses loci into a single representative isoform (Most highly expressed isoform is chosen)
- Filters transcripts based on expression levels 

After these steps, predicted lncRNA were manually assessed and re-annotated where necessary. Following this manual annotation step, the following steps of the pre-manual annotation pipeline were then repeated:
- Functional annotation
- Expression filter

The pipeline is run using the following scripts:
1. **stringtie.sh** - Uses StringTie to assemble transcripts from RNAseq data.
2. **lncRNA_pipeline_pre_MA_part1.sh** - Removes known annotated features, transcripts extremely close to genes,performs functional annotation and collapses loci.
3. **lncRNA_pipeline_pre_MA_part2.sh** - Filters transcripts based on expression levels.
4. **lncRNA_pipeline_post_MA_part1.sh** - Performs functional annotation
5. **lncRNA_pipeline_post_MA_part2.sh** - Filters transcripts based on expression levels.

All additional scripts which are required by the above shell scripts are contained within 'required_files' folder.

Additional data needed for the pipeline are: 
- Transcriptomic data (fastq files, sorted BAM files and both HTSeq and Salmon-generated count data)
- BLAST search results (*.out file)
- Rfam search results (*.tblout file)
- CPC2 results (*result_cpc2.txt file)

Examples of the code used to generate these data are summarised in "additional_data_code.txt"
