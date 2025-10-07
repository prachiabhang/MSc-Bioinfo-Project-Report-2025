# MSc-Bioinfo-Project-Report-2025
MSc Bioinfo Project Report 2025


# Genetic Differentiation and Phylogenetic Analysis of the Human OPN1SW Gene Across Arctic and Equatorial Populations

 **Under the guidance of Dr. Davide Pisani**
 A big thanks to **Dr. Christopher Kay**. The following scripts were adapted and developed during the course of this MSc project, with valuable inputs and support.


This repository contains the scripts, analysis workflows, and figures 
associated with my MSc Bioinformatics dissertation at the University of Bristol (2025).


## Project Overview
The **OPN1SW** gene encodes the short-wavelength sensitive opsin (blue cone photopigment), 
critical for human colour vision. This study investigates whether OPN1SW shows signatures of **latitude-linked adaptation** 
between Arctic and equatorial populations using genomic datasets and population genetic analyses.

- **Hypothesis**: OPN1SW haplotypes show subtle but detectable divergence between populations living under contrasting light environments.
- **Datasets**:  
  - 1000 Genomes Project (VCF-based haplotypes)  
  - HGDP/SGDP Arctic cohort (FASTQ >> consensus haplotypes)  
- **Methods**: haplotype reconstruction, QC, multiple sequence alignment (MAFFT), 
  phylogenetic inference (IQ-TREE), population differentiation (Hudson’s FST).     

## Workflow Summary
1. Extract OPN1SW ±5 kb locus from 1000 Genomes (VCFs).  
2. Generate consensus haplotypes with bcftools consensus.  
3. For HGDP/SGDP Arctic samples: FASTQ >> trimming (fastp) >> alignment (bwa-mem) >> variant calling (FreeBayes) >> consensus FASTA.  
4. QC of haplotypes (Ns, GC%, duplicates).  
5. Redundancy reduction with CD-HIT.  
6. Alignment (MAFFT), phylogeny (IQ-TREE).  
7. Population differentiation with Hudson’s FST (scikit-allel): https://www.kaggle.com/code/dprdx/fst-vfinal

## List of scripts used 

## 1. Consensus_haplotype_sequence_reconstruction_pipeline.sh
A master script that takes 1000 Genomes VCFs and turns them into haplotype FASTA files.
(In short: VCF >> FASTA haplotypes.)

## 2. sbatch_Arc.sh
A simple Slurm job submission script to run the Arctic sample pipeline on the HPC cluster.
(Think of it as the “launcher” for Arctic processing.)

## 3. FASTQ_to_VCF_FAST_Arc.sh
The full pipeline for Arctic samples. It starts from raw FASTQs, trims and aligns them, calls variants, and outputs consensus FASTA haplotypes.
(FASTQ >> BAM >> VCF >> FASTA.)

## 4. run_fasta_qc_allinone.sh
Runs quality control checks on the reconstructed haplotypes. It checks sequence length, GC content, missing bases (Ns), and flags duplicates.
(Basically: QC for your haplotype FASTAs.)

## Important Figures..
-------------------------------------------------------------------------------------------------------

Table 1a. Genomic coordinates and annotation of the human OPN1SW gene (GRCh37/hg19, Ensembl release 75).
 <img width="950" height="125" alt="image" src="https://github.com/user-attachments/assets/3e3a4e86-b297-442b-bb17-d3b6c3ab7985" />

-------------------------------------------------------------------------------------------------------


Table 1b. Exon-level annotation of OPN1SW (GRCh37/hg19, Ensembl release 75).
<img width="822" height="175" alt="image" src="https://github.com/user-attachments/assets/2aa615b0-1caf-4e4a-b6db-0aaa5915d629" />

-------------------------------------------------------------------------------------------------------

Figure 1. Genomic structure of the human OPN1SW locus (GRCh37, chromosome 7). Schematic of OPN1SW (chr7:128,407,545–128,420,844, GRCh37) with ±5 kb flanks. Exons (E1-E5, dark blue), introns (I1-I4, light blue), and 5′/3′ flanking regions (grey) are shown. Local coordinates (bottom) span the 13.3 kb extracted region, genomic coordinates (top) indicate positions on chromosome 7.
<img width="965" height="200" alt="image" src="https://github.com/user-attachments/assets/3695f92c-7fb2-4532-9ec3-53639a64e8b8" />

-------------------------------------------------------------------------------------------------------
Table 2. Final dataset composition after quality control and filtering
<img width="853" height="463" alt="image" src="https://github.com/user-attachments/assets/7e1283ec-0cb5-4cbb-819e-80685fa35314" />

-------------------------------------------------------------------------------------------------------

Figure 2. Workflow for haplotype consensus sequence reconstruction at the OPN1SW locus. Pipeline illustrating the steps from raw input data to haplotype FASTA output. Chromosome-wide 1000 Genomes VCFs were subset for the OPN1SW region (chr7:128,407,545–128,420,844, GRCh37, ±5 kb) using tabix while the corresponding reference sequence was retrieved from Ensembl and indexed. Sample identifiers were extracted with bcftools query -l, and phased haplotype sequences were reconstructed per individual using bcftools consensus. The resulting haplotype FASTAs (~13 kb each) were used for further downstream analysis.
<img width="1004" height="247" alt="image" src="https://github.com/user-attachments/assets/a3975fe0-ea45-4cc3-9530-0ac0c279fb05" />

-------------------------------------------------------------------------------------------------------
Figure 3. Workflow for consensus sequence reconstruction from Arctic HGDP/SGDP samples. Paired-end FASTQ reads were retrieved from the IGSR portal and processed through a custom pipeline. Reads were adapter-trimmed and quality-filtered with fastp, then aligned to the GRCh37 reference (regional slice ±5 kb around OPN1SW) using bwa-mem. Variants within the target interval were called with FreeBayes, filtered with bcftools, and compressed/indexed with bgzip/tabix. Consensus sequences were reconstructed with bcftools consensus, encoding heterozygous sites with IUPAC ambiguity codes and applying locus-specific ploidy rules. Outputs included per-sample consensus FASTAs (~13 kb), filtered VCFs, and indexed BAMs, which were pooled with 1000 Genomes haplotypes for downstream alignment and phylogenetic analysis.
<img width="994" height="244" alt="image" src="https://github.com/user-attachments/assets/291376d1-a4a3-4320-8fe8-cf1e6618a975" />

-------------------------------------------------------------------------------------------------------

<img width="940" height="489" alt="image" src="https://github.com/user-attachments/assets/1468611c-ab4b-4ab6-b638-0e39b2a070a2" />


-------------------------------------------------------------------------------------------------------
<img width="940" height="530" alt="image" src="https://github.com/user-attachments/assets/d6bfa343-a888-443a-9074-1017a12c5bf2" />

-------------------------------------------------------------------------------------------------------
<img width="940" height="583" alt="image" src="https://github.com/user-attachments/assets/ed7c2f2e-e847-43b2-95ec-332935a04683" />

-------------------------------------------------------------------------------------------------------
Figure 6a. Per-site FST Manhattan plot with exon/intron overlay.
<img width="952" height="430" alt="image" src="https://github.com/user-attachments/assets/99347d69-4190-4927-b0e1-f2d211d0ed4f" />

-------------------------------------------------------------------------------------------------------
Figure 6b. Windowed Skyline FST with feature annotation.
<img width="940" height="470" alt="image" src="https://github.com/user-attachments/assets/86c40df5-befb-4eac-9b35-a9e2bcaf9640" />


-------------------------------------------------------------------------------------------------------
Figure 6c. Scatterplot of alternate-allele frequencies in Arctic vs. equatorial populations.
<img width="843" height="632" alt="image" src="https://github.com/user-attachments/assets/b7f86ead-5f0d-4a98-b011-4b3c6cda634b" />

-------------------------------------------------------------------------------------------------------



## Citation
If you use this repository or workflows, please cite as:

```bibtex
@misc{PrachiAbhang_OPN1SW_2025,
  author       = {Prachi R Abhang},
  title        = {Genetic Differentiation and Phylogenetic Analysis of the Human OPN1SW Gene Across Arctic and Equatorial Populations)},
  year         = {2025},
  howpublished = {\url{https://github.com/your-username/opsin-evolution}},
  note         = {MSc Bioinformatics Dissertation, University of Bristol}
}
