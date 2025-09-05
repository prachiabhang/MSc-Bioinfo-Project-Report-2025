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

## 📂 Repository Structure
- `scripts/` – Bash and Python pipelines (haplotype reconstruction, QC, FST analysis).  
- `results/` – Figures and tables from the dissertation.  
- `data/` – Processed haplotypes and QC summaries (raw data from 1KG/SGDP not redistributed).  
- `docs/` – Supplementary materials and dissertation report.  

## ⚙️ Workflow Summary
![Pipeline](results/figs/Figure2_pipeline.png)

1. Extract OPN1SW ±5 kb locus from 1000 Genomes (VCFs).  
2. Generate consensus haplotypes with `bcftools consensus`.  
3. For HGDP/SGDP Arctic samples: FASTQ -> trimming (`fastp`) -> alignment (`bwa-mem`) -> variant calling (`FreeBayes`) -> consensus FASTA.  
4. QC of haplotypes (Ns, GC%, duplicates).  
5. Redundancy reduction with CD-HIT.  
6. Alignment (MAFFT), phylogeny (IQ-TREE).  
7. Population differentiation with Hudson’s FST (scikit-allel).  

## 📊 Key Figures
- **Phylogeny of Arctic vs. Equatorial haplotypes**  


- **Population differentiation (Hudson’s FST)**  
 

## 📖 Citation
If you use this repository or workflows, please cite as:

```bibtex
@misc{YourName_OPN1SW_2025,
  author       = {Your Name},
  title        = {Genetic Differentiation and Phylogenetic Analysis of the Human OPN1SW Gene Across Arctic and Equatorial Populations)},
  year         = {2025},
  howpublished = {\url{https://github.com/your-username/opsin-evolution}},
  note         = {MSc Bioinformatics Dissertation, University of Bristol}
}
