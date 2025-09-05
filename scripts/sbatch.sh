#!/bin/bash

#SBATCH --job-name=FASTQ
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=24
#SBATCH --time=96:0:00
#SBATCH --mem=32G
#SBATCH --account=BISC034804


./FASTQ_to_VCF_FAST.sh "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR132/003/ERR1329993/ERR1329993_1.fastq.gz" "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR132/003/ERR1329993/ERR1329993_2.fastq.gz" "HGDP00967" "autosomal" "7" "128767485" "128780794" "24"

./FASTQ_to_VCF_FAST.sh "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR757/ERR757810/ERR757810_1.fastq.gz" "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR757/ERR757810/ERR757810_2.fastq.gz" "HGDP00948" "autosomal" "7" "128767485" "128780794" "24"

./FASTQ_to_VCF_FAST.sh "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR133/003/ERR1330003/ERR1330003_1.fastq.gz" "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR133/003/ERR1330003/ERR1330003_2.fastq.gz" "HGDP00950" "autosomal" "7" "128767485" "128780794" "24"

./FASTQ_to_VCF_FAST.sh "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR134/002/ERR1343902/ERR1343902_1.fastq.gz" "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR134/002/ERR1343902/ERR1343902_2.fastq.gz" "HGDP00955" "autosomal" "7" "128767485" "128780794" "24"

./FASTQ_to_VCF_FAST.sh "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR134/001/ERR1344251/ERR1344251_1.fastq.gz" "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR134/001/ERR1344251/ERR1344251_2.fastq.gz" "HGDP00962" "autosomal" "7" "128767485" "128780794" "24"

./FASTQ_to_VCF_FAST.sh "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR134/009/ERR1344479/ERR1344479_1.fastq.gz" "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR134/009/ERR1344479/ERR1344479_1.fastq.gz" "HGDP00963" "autosomal" "7" "128767485" "128780794" "24"

./FASTQ_to_VCF_FAST.sh "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR134/008/ERR1344108/ERR1344108_1.fastq.gz" "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR134/008/ERR1344108/ERR1344108_2.fastq.gz" "HGDP00968" "autosomal" "7" "128767485" "128780794" "24"

./FASTQ_to_VCF_FAST.sh "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR134/000/ERR1344190/ERR1344190_1.fastq.gz" "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR134/000/ERR1344190/ERR1344190_2.fastq.gz" "HGDP00949" "autosomal" "7" "128767485" "128780794" "24"

./FASTQ_to_VCF_FAST.sh "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR134/007/ERR1344767/ERR1344767_1.fastq.gz" "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR134/007/ERR1344767/ERR1344767_2.fastq.gz" "HGDP00964" "autosomal" "7" "128767485" "128780794" "24"

./FASTQ_to_VCF_FAST.sh "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR134/004/ERR1343844/ERR1343844_1.fastq.gz" "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR134/004/ERR1343844/ERR1343844_2.fastq.gz" "HGDP00969" "autosomal" "7" "128767485" "128780794" "24"

./FASTQ_to_VCF_FAST.sh "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR757/ERR757808/ERR757808_1.fastq.gz" "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR757/ERR757808/ERR757808_2.fastq.gz" "HGDP00945" "autosomal" "7" "128767485" "128780794" "24"

./FASTQ_to_VCF_FAST.sh "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR134/006/ERR1344646/ERR1344646_1.fastq.gz" "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR134/006/ERR1344646/ERR1344646_2.fastq.gz" "HGDP00952" "autosomal" "7" "128767485" "128780794" "24"

./FASTQ_to_VCF_FAST.sh "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR134/008/ERR1344538/ERR1344538_1.fastq.gz" "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR134/008/ERR1344538/ERR1344538_2.fastq.gz" "HGDP00957" "autosomal" "7" "128767485" "128780794" "24"

./FASTQ_to_VCF_FAST.sh "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR133/004/ERR1330074/ERR1330074_1.fastq.gz" "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR133/004/ERR1330074/ERR1330074_2.fastq.gz" "HGDP00946" "autosomal" "7" "128767485" "128780794" "24"

./FASTQ_to_VCF_FAST.sh "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR134/006/ERR1344226/ERR1344226_1.fastq.gz" "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR134/006/ERR1344226/ERR1344226_2.fastq.gz" "HGDP00953" "autosomal" "7" "128767485" "128780794" "24"

./FASTQ_to_VCF_FAST.sh "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR134/002/ERR1344322/ERR1344322_1.fastq.gz" "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR134/002/ERR1344322/ERR1344322_2.fastq.gz" "HGDP00958" "autosomal" "7" "128767485" "128780794" "24"

./FASTQ_to_VCF_FAST.sh "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR134/009/ERR1343819/ERR1343819_1.fastq.gz" "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR134/009/ERR1343819/ERR1343819_2.fastq.gz" "HGDP00960" "autosomal" "7" "128767485" "128780794" "24"

./FASTQ_to_VCF_FAST.sh "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR134/007/ERR1344527/ERR1344527_1.fastq.gz" "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR134/007/ERR1344527/ERR1344527_2.fastq.gz" "HGDP00965" "autosomal" "7" "128767485" "128780794" "24"

./FASTQ_to_VCF_FAST.sh "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR134/005/ERR1344335/ERR1344335_1.fastq.gz" "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR134/005/ERR1344335/ERR1344335_2.fastq.gz" "HGDP00966" "autosomal" "7" "128767485" "128780794" "24"

./FASTQ_to_VCF_FAST.sh "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR757/ERR757809/ERR757809_1.fastq.gz" "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR757/ERR757809/ERR757809_2.fastq.gz" "HGDP00947" "autosomal" "7" "128767485" "128780794" "24"

./FASTQ_to_VCF_FAST.sh "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR134/000/ERR1344490/ERR1344490_1.fastq.gz" "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR134/000/ERR1344490/ERR1344490_2.fastq.gz" "HGDP00954" "autosomal" "7" "128767485" "128780794" "24"

./FASTQ_to_VCF_FAST.sh "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR134/007/ERR1344167/ERR1344167_1.fastq.gz" "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR134/007/ERR1344167/ERR1344167_2.fastq.gz" "HGDP00959" "autosomal" "7" "128767485" "128780794" "24"

./FASTQ_to_VCF_FAST.sh "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR134/001/ERR1343831/ERR1343831_1.fastq.gz" "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR134/001/ERR1343831/ERR1343831_2.fastq.gz" "HGDP00961" "autosomal" "7" "128767485" "128780794" "24"
