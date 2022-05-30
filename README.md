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

```bash
fastq-dump -I --split-files --gzip $(</scratch/tt61/transcriptomics/SRR_Acc_List.txt)
fastq-dump -I --split-files --gzip $(</scratch/dq2033/transcriptomics/SRR_Acc_List.txt)
```

### Data processing

After the reads were downloaded, they had to be inspected for quality scores and adapters needed to be removed. First, fastp/intel/0.20.1 was used in order to trim adapter sequences. For each mouse model, a table was constructed. The table consisted of the SRA run ID, and the names for the paired-end reads. The fastp command was made so that it iterated through every line and column of the table to trim each sample and generate its corresponding output files. 

The first 2 lines of the table:

```bash
SRR8512365 SRR8512365_1.fastq.gz SRR8512365_2.fastq.gz
SRR8512366 SRR8512366_1.fastq.gz SRR8512366_2.fastq.gz
```
The fastp command that was used was:

```bash
fastp -i $fqdir/$fq1 \
-I $fqdir/$fq2 \
-o $fq1_fastp \
-O $fq2_fastp \
--length_required 76 \
--detect_adapter_for_pe \
--n_base_limit 50 \
--html $sample.fastp.html \
--json $sample.fastp.json
```
fastqc/0.11.9 was then used to generate quality control reports for each sample:

```bash
fastqc $fq1_fastp $fq2_fastp
```

MultiQC was used in order to group the fastqc reports for each mouse model separately and visualize the metrics for each sample all in one report. 
Alignment of Samples to reference Genome

In order to align the samples to the reference genome, STAR (Spliced Transcripts Alignment to a Reference) was used. STAR is meant for short-read whole-transcriptome sequencing data. The first step in aligning reads using STAR is to create a genome index. In order to create the genome index using star, the reference genome sequence is needed in FASTA format along with the GTF annotation file. The paper where the data came from used the GRCm38.p4 (M11 release). Both the FASTA and GTF files were downloaded from GENCODE using wget:

```bash
wget 'https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M6/gencode.vM6.chr_patch_hapl_scaff.annotation.gtf.gz'
wget 'https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M6/GRCm38.p4.genome.fa.gz'
```

Once the FASTA reference and GTF annotation were downloaded, the files have to be unzipped before they can be used by STAR:

```bash
gunzip gencode.vM6.chr_patch_hapl_scaff.annotation.gtf.gz
gunzip GRCm38.p4.genome.fa.gz
```

The star command that was used to generate the genome index was:

```bash
STAR --runThreadN 20 --runMode genomeGenerate --genomeDir /scratch/dq2033/transcriptomics/J20_samples/genome --genomeFastaFiles /scratch/dq2033/transcriptomics/J20_samples/genome/GRCm38.p4.genome.fa --sjdbGTFfile /scratch/dq2033/transcriptomics/J20_samples/genome/gencode.vM6.chr_patch_hapl_scaff.annotation.gtf
```

Next, the goal was to align the reads to the genome and count them simultaneously. Because the next step in the analysis is the DESeq2 workflow, it was necessary to obtain raw gene counts as the input for DESeq2. STAR can output raw gene counts with the option:

```bash
--quantmode GeneCounts
```

Because STAR outputs the counts in a ReadsPerGene.out.tab that has several columns according to the strandedness of the data, it was necessary to check whether the libraries were stranded or not. In order to check, an awk command line was used:

```bash
$ grep -v "N_" SRR8512365ReadsPerGene.out.tab | awk '{unst+=$2;forw+=$3;rev+=$4}END{print unst,forw,rev}'
8779853 8870916 314680
```

The output of this command indicated that 8870916 Reads 1 (forward) were mapped to known genes and  314680 Reads 2 (reverse) were mapped to known genes. The fact that these numbers were so different indicates that the RNA library preparation was stranded as stated in  Castanho et. al: ‘ Stranded-specific mRNA sequencing libraries were prepared using the TruSeq Stranded mRNA Sample Prep Kit (Illumina) using the Bravo Automated Liquid Handling Platform (Agilent)’.
The counts for each read were then concatenated into a single dataframe to use as input for differential gene expression analysis. 

## Differential Gene Expression Analysis with DESeq2

Normalized DESeq2 read counts were plotted for each of the 5 comparisons we outlined in the methods. The plots show hierarchical clustering of mean levels of expression across all samples. The two mouse models, J20 and rTg4510 (figure 2a) show clear distinction between each model suggesting the relative expression of genes are different with minimal overlap. Comparisons of wild-type vs mutants (figure 2b-c) and young mice vs adult mice (figure 2d-e) showed samples clustered by genotype and age. J20 mice did not have many differential expressions of genes, of note were Ccdc80 (upregulated in mutant mice, FDR < 0.005, fold change >2 ), Abca8a (upregulated in mutant mice, FDR < 0.005, fold change > 2 ). Ccdc82 is involved in cell adhesion and matrix assembly and Abca8a encodes a member of the ATP-binding cassette (ABC, A-subclass) transporter family involved in homeostatic regulation of brain lipid homeostasis (1, 12). Htr1a, which encodes a major serotonin receptor and Hspa5, involved in folding and assembly of proteins,  were both downregulated in mutant mice (FDR < 0.005, fold change < 0.5) (1). These results are consistent with that found in the reference paper, except for gene Abca8a, which was reported as downregulated in mutants with an FDR of 0.02. rTg4510 had  a lot more variation in expression levels, of note were genes such as Car4, responsible for encoding carbonic anhydrase 4 (upregulated in mutant mice, FDR<0.005, fold change> 2), Blnk (upregulated in mutant mice, FDR < 0.005, fold change ~ 3 ) involved in B cell development both consistent with reference paper (1). 
J20 age comparison did not show much difference, also consistent with that reported in background paper. rTg4510 model, showed a higher range of expression levels, of note were Gfap, Cd68, Itgax genes that are key microglial markers and Trem2, a key drive of AD pathology (1).

![Screen Shot 2022-05-30 at 11 09 12 AM](https://user-images.githubusercontent.com/90015489/171020336-1d2e21eb-ed5c-48e2-bcaf-0adcf5ecefc7.png)

![Screen Shot 2022-05-30 at 11 11 06 AM](https://user-images.githubusercontent.com/90015489/171020635-1bea5e43-af48-4119-8ec9-99faede1737d.png)

Figure 2. Heatmaps of all comparisons. (a) J20 vs rTG4510 (top left), (b) rTg4510: wild-type vs mutants (top right), (c) J20: wild-type vs mutants (second from top right), (d) J20: 6m vs 12m (bottom left) and (e) rTg4510: 2m vs 8m (second from bottom left). 

## GO-Term Enrichment

### J20 wildtype vs mutants
GO-term enrichment analysis of the most differentially expressed genes in J20 wild-type vs mutants yielded very few terms within the threshold of FDR < 0.05 (Figure 3). The terms include: cytoplasmic side of plasma membrane, extrinsic side of plasma membrane and embryonic limb morphogenesis, each containing about 5-8 DE genes under each term. The DE genes are better visualized in figure 4, a category netplot that shows the relationships between the genes associated with the most significant GO terms and the fold changes of the significant genes associated with these terms. Gnai3 gene, which shows a fold reduction (8) is involved in making the inhibitory component of alpha subunit for a G protein that in turn inhibits the activity of adenylyl cyclase, AC (8). AC dysregulation is a precursor of amyloid plaque deposition (2). A fold reduction in this subunit would suggest an increase in AC activity and consequently, amyloid plaque formation (12). Another gene that is downregulated is Cav2. Cav2 encodes a major membrane protein that is implicated in essential cellular functions such as signal transduction and lipid metabolism (8). Amyloid plaque formation reports altered cholesterol metabolism and hypercholesterolemia in AD progression, thus Cav2 activity can act as a biomarker in AD pathology (11). 
### rTg4510 wildtype vs mutants
The top GO terms for rTg4510 wild-type vs mutants were: ameboidal-type cell migration, establishment of protein localization to organelle and protein import, all within FDR < 0.05. There were many genes implicated, some with prominent fold changes include: Tbx2, Sp1, Akt1 and KIT. Most of these genes are either implicated in encoding transcription factors and kinases (8). Tau pathology involves hyperphosphorylation of the tau protein thus it strengthens our hypothesis to observe that these genes are dysregulated. Figures 5 and 6 shows the dotplot and category net plots of this analysis. 

### J20: 6m vs 12m
The only comparison that we conducted that did not survive p value correction was the comparison of young and old J20 mice and WT and mutant J20 mice. This is not that surprising since we had seen that J20 mice were very tightly clustered together on the PCA graph (Fig.1). Although there were fewer differentially expressed genes when comparing old and young J20 mice, a more clean picture was obtained when comparing J20 mice in terms of mutant and wild-type genotype.
### rTg4510: 2m vs 8m
Age comparison between 2m and 8m old rTg4510 mice yielded similar terms to that observed in WT vs mutant mice however the fold changes were larger on both ends (figure 7 and 8). Some key genes of notable fold changes include: the ACVRL1 gene, which is involved in encoding activin receptor-like kinase  provides instructions for making a protein called activin receptor-like kinase, implicated in growth factor signal transduction (8). CDKN1A, a protein involved in cell cycle progression and protein import also showed a fold reduction, suggesting disruption in basal cell cycle processes that is caused by tau hyperphosphorylation (11).

### J20: 6m vs 12m
The only comparison that we conducted that did not survive p value correction was the comparison of young and old J20 mice and WT and mutant J20 mice. This is not that surprising since we had seen that J20 mice were very tightly clustered together on the PCA graph (Fig.1). Although there were fewer differentially expressed genes when comparing old and young J20 mice, a more clean picture was obtained when comparing J20 mice in terms of mutant and wild-type genotype.
### rTg4510: 2m vs 8m
Age comparison between 2m and 8m old rTg4510 mice yielded similar terms to that observed in WT vs mutant mice however the fold changes were larger on both ends (figure 7 and 8). Some key genes of notable fold changes include: the ACVRL1 gene, which is involved in encoding activin receptor-like kinase  provides instructions for making a protein called activin receptor-like kinase, implicated in growth factor signal transduction (8). CDKN1A, a protein involved in cell cycle progression and protein import also showed a fold reduction, suggesting disruption in basal cell cycle processes that is caused by tau hyperphosphorylation (11).

![Screen Shot 2022-05-30 at 11 12 32 AM](https://user-images.githubusercontent.com/90015489/171020853-ace3a242-0c1a-4c35-aa35-25b5f273cc3f.png)
![Screen Shot 2022-05-30 at 11 13 47 AM](https://user-images.githubusercontent.com/90015489/171021081-f89a6a7f-9e9d-4df3-8040-c8efcfdd3a08.png)
![Screen Shot 2022-05-30 at 11 14 11 AM](https://user-images.githubusercontent.com/90015489/171021127-adb10d84-4674-49b7-817c-c439c595b7ef.png)
![Screen Shot 2022-05-30 at 11 15 27 AM](https://user-images.githubusercontent.com/90015489/171021327-8e285e93-406d-4f69-afa8-aa03db0bb45d.png)



