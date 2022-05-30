# Transcriptomics-Final-Project: Analyzing genes associated with Alzheimer’s Disease Progression

![Screen Shot 2022-05-30 at 10 56 24 AM](https://user-images.githubusercontent.com/90015489/171018133-e7965791-c9f8-4155-bea7-6d774e2c466c.png)

## Abstract
This report investigates the kinds of genes that are differentially expressed in J20 and rTg4510 mice that are used to model the pathology of Alzheimer’s Disease (AD) based on the published work of Hill et al. in Cell Reports titled “Transcriptional Signatures of Tau and Amyloid Neuropathology”. We investigated what kinds of genes were differentially expressed when comparing the two strains. We also explored how the gene expression changed as mice age by comparing young and old mice within each strain group. Lastly, we conducted GO term enrichment and Signaling Pathway Impact Analysis (SPIA) in R. We found more genes were differentially expressed in rTg4510 mice model than J20 mice models (2223 vs . rTg4510 mice expresses very different genes as they age, whereas old and young J20 mice are more similar to each other. And lastly, two signaling pathways that could be of significance in studying tau pathology are focal adhesion and MAPK signaling pathways because it consisted of almost all the DEGs identified in rTg4510 mice samples. 

## Workflow Summary

![Screen Shot 2022-05-30 at 11 00 18 AM](https://user-images.githubusercontent.com/90015489/171018834-46cefa20-2f6a-4517-b19e-2394a69058da.png)


## Data Acquisition and Pre-processing Summary

These are the steps that we took in order to pre-process all of the files. The total number of paired end samples that were processed was 72. In order to see the detailed/full scripts that were used throughout the analysis please refer to the Appendix. We decided to include brief explanations of the data download and processing since they represented a significant portion of the project. 

### Data download

The data was downloaded from the Sequence Read Archive (SRA) using the sra-tools/2.10.9 module on the NYU Greene Cluster and the fastq-dump command. The samples for each mouse line, 42 for J20 and 30 for rtg4510, were downloaded in parallel using two corresponding slurm scripts. The accession numbers for the relevant samples for each mouse line were stored in two text files that were passed to the fastq-dump command to download the samples of interest. The command line that was used to download each of the samples was: 

fastq-dump -I --split-files --gzip $(</scratch/tt61/transcriptomics/SRR_Acc_List.txt)
fastq-dump -I --split-files --gzip $(</scratch/dq2033/transcriptomics/SRR_Acc_List.txt)

### Data processing

After the reads were downloaded, they had to be inspected for quality scores and adapters needed to be removed. First, fastp/intel/0.20.1 was used in order to trim adapter sequences. For each mouse model, a table was constructed. The table consisted of the SRA run ID, and the names for the paired-end reads. The fastp command was made so that it iterated through every line and column of the table to trim each sample and generate its corresponding output files. 

The first 2 lines of the table:
'''bash
SRR8512365 SRR8512365_1.fastq.gz SRR8512365_2.fastq.gz
SRR8512366 SRR8512366_1.fastq.gz SRR8512366_2.fastq.gz
'''
The fastp command that was used was:

fastp -i $fqdir/$fq1 \
-I $fqdir/$fq2 \
-o $fq1_fastp \
-O $fq2_fastp \
--length_required 76 \
--detect_adapter_for_pe \
--n_base_limit 50 \
--html $sample.fastp.html \
--json $sample.fastp.json

fastqc/0.11.9 was then used to generate quality control reports for each sample:

fastqc $fq1_fastp $fq2_fastp
MultiQC was used in order to group the fastqc reports for each mouse model separately and visualize the metrics for each sample all in one report. 
Alignment of Samples to reference Genome

In order to align the samples to the reference genome, STAR (Spliced Transcripts Alignment to a Reference) was used. STAR is meant for short-read whole-transcriptome sequencing data. The first step in aligning reads using STAR is to create a genome index. In order to create the genome index using star, the reference genome sequence is needed in FASTA format along with the GTF annotation file. The paper where the data came from used the GRCm38.p4 (M11 release). Both the FASTA and GTF files were downloaded from GENCODE using wget:

wget 'https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M6/gencode.vM6.chr_patch_hapl_scaff.annotation.gtf.gz'
wget 'https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M6/GRCm38.p4.genome.fa.gz'

Once the FASTA reference and GTF annotation were downloaded, the files have to be unzipped before they can be used by STAR:

gunzip gencode.vM6.chr_patch_hapl_scaff.annotation.gtf.gz
gunzip GRCm38.p4.genome.fa.gz

The star command that was used to generate the genome index was:

STAR --runThreadN 20 --runMode genomeGenerate --genomeDir /scratch/dq2033/transcriptomics/J20_samples/genome --genomeFastaFiles /scratch/dq2033/transcriptomics/J20_samples/genome/GRCm38.p4.genome.fa --sjdbGTFfile /scratch/dq2033/transcriptomics/J20_samples/genome/gencode.vM6.chr_patch_hapl_scaff.annotation.gtf

Next, the goal was to align the reads to the genome and count them simultaneously. Because the next step in the analysis is the DESeq2 workflow, it was necessary to obtain raw gene counts as the input for DESeq2. STAR can output raw gene counts with the option:

--quantmode GeneCounts

Because STAR outputs the counts in a ReadsPerGene.out.tab that has several columns according to the strandedness of the data, it was necessary to check whether the libraries were stranded or not. In order to check, an awk command line was used:

$ grep -v "N_" SRR8512365ReadsPerGene.out.tab | awk '{unst+=$2;forw+=$3;rev+=$4}END{print unst,forw,rev}'
8779853 8870916 314680

The output of this command indicated that 8870916 Reads 1 (forward) were mapped to known genes and  314680 Reads 2 (reverse) were mapped to known genes. The fact that these numbers were so different indicates that the RNA library preparation was stranded as stated in  Castanho et. al: ‘ Stranded-specific mRNA sequencing libraries were prepared using the TruSeq Stranded mRNA Sample Prep Kit (Illumina) using the Bravo Automated Liquid Handling Platform (Agilent)’.
The counts for each read were then concatenated into a single dataframe to use as input for differential gene expression analysis. 

