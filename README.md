# MSc-Bioinfo-Project-Report-2025
MSc Bioinfo Project Report 2025


# Genetic Differentiation and Phylogenetic Analysis of the Human OPN1SW Gene Across Arctic and Equatorial Populations

This repository contains the scripts, analysis workflows, and figures 
associated with my MSc Bioinformatics dissertation at the University of Bristol (2025).

## Project Overview
The **OPN1SW** gene encodes the short-wavelength–sensitive opsin (blue cone photopigment), 
critical for human colour vision. This study investigates whether OPN1SW shows signatures of **latitude-linked adaptation** 
between Arctic and equatorial populations using genomic datasets and population genetic analyses.

- **Hypothesis**: OPN1SW haplotypes show subtle but detectable divergence between populations living under contrasting light environments.
- **Datasets**:  
  - 1000 Genomes Project (VCF-based haplotypes)  
  - HGDP/SGDP Arctic cohort (FASTQ -- consensus haplotypes)  
- **Methods**: haplotype reconstruction, QC, multiple sequence alignment (MAFFT), 
  phylogenetic inference (IQ-TREE), population differentiation (Hudson’s FST).     

## Workflow Summary
1. Extract OPN1SW ±5 kb locus from 1000 Genomes (VCFs).  
2. Generate consensus haplotypes with `bcftools consensus`.  
3. For HGDP/SGDP Arctic samples: FASTQ >> trimming (fastp) >> alignment (bwa-mem) >> variant calling (FreeBayes) >> consensus FASTA.  
4. QC of haplotypes (Ns, GC%, duplicates).  
5. Redundancy reduction with CD-HIT.  
6. Alignment (MAFFT), phylogeny (IQ-TREE).  
7. Population differentiation with Hudson’s FST (scikit-allel): https://www.kaggle.com/code/dprdx1/fst-vfina

## List of scripts used 
A big thanks to **Dr. Christopher Kay**. The following scripts were adapted and developed during the course of this MSc project, with valuable input and support.

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

Table 1a. Genomic coordinates and annotation of the human OPN1SW gene (GRCh37/hg19, Ensembl release 75).
 <img width="950" height="125" alt="image" src="https://github.com/user-attachments/assets/3e3a4e86-b297-442b-bb17-d3b6c3ab7985" />


Table 1b. Exon-level annotation of OPN1SW (GRCh37/hg19, Ensembl release 75).
<img width="822" height="175" alt="image" src="https://github.com/user-attachments/assets/2aa615b0-1caf-4e4a-b6db-0aaa5915d629" />

Figure 1. Genomic structure of the human OPN1SW locus (GRCh37, chromosome 7). Schematic of OPN1SW (chr7:128,407,545–128,420,844, GRCh37) with ±5 kb flanks. Exons (E1-E5, dark blue), introns (I1-I4, light blue), and 5′/3′ flanking regions (grey) are shown. Local coordinates (bottom) span the 13.3 kb extracted region, genomic coordinates (top) indicate positions on chromosome 7.
<img width="965" height="200" alt="image" src="https://github.com/user-attachments/assets/3695f92c-7fb2-4532-9ec3-53639a64e8b8" />

<img width="853" height="463" alt="image" src="https://github.com/user-attachments/assets/6b86f37f-f55d-47fb-b024-065b5e1685d6" />

Table 2. Final dataset composition after quality control and filtering
<img width="853" height="463" alt="image" src="https://github.com/user-attachments/assets/7e1283ec-0cb5-4cbb-819e-80685fa35314" />

## Citation
If you use this repository or workflows, please cite as:

```bibtex
@misc{YourName_OPN1SW_2025,
  author       = {Prachi R Abhang},
  title        = {Genetic Differentiation and Phylogenetic Analysis of the Human OPN1SW Gene Across Arctic and Equatorial Populations)},
  year         = {2025},
  howpublished = {\url{https://github.com/your-username/opsin-evolution}},
  note         = {MSc Bioinformatics Dissertation, University of Bristol}
}
