#!/bin/bash

URL1=$1
URL2=$2
SAMPLE_NAME=$3
PLOIDY_TYPE=$4          # autosomal or X
chromosome=$5           # "7"
start_position=$6       # "128767485"
end_position=$7         # "128780794"
threads=$8

# Determine ploidy
case "$PLOIDY_TYPE" in
  autosomal) PLOIDY_NUM=2 ;;
  xlinked) PLOIDY_NUM=1 ;;  # assuming male
  *) echo "Unsupported PLOIDY: $PLOIDY_TYPE. Use 'autosomal' or 'xlinked'."; exit 1 ;;
esac

# Check environment
available_kb=$(df --output=avail . | tail -1)
available_gb=$((available_kb / 1024 / 1024))

if (( available_gb < 10 )); then
  echo "Error: Less than 10GB of disk space available. Exiting."
  exit 1
elif (( available_gb < 50 )); then
  echo "Warning: Less than 50GB of disk space available."
fi

echo "# Checking script environment"

status=0
required_programs=("wget" "curl" "gunzip" "samtools" "fastp" "bwa" "freebayes" "bcftools")
for program in "${required_programs[@]}"; do
  if ! command -v "$program" &> /dev/null; then
    echo "FAIL: '$program' is not installed or not in the system's PATH"
    status=1
  else
    echo "PASS: $program is available"
  fi
done

if ping -c 1 -W 2 8.8.8.8 > /dev/null 2>&1; then
  echo "Internet connection is available."
else
  echo "No internet connection."
  status=1
fi

if [[ $status == "1" ]]; then
  exit 1
fi

# Reference handling
REFERENCE_FULL="chr${chromosome}.fa"
REFERENCE_FILE="chr${chromosome}_${start_position}_${end_position}.fa"

if [ ! -s "$REFERENCE_FULL" ]; then
  ARCHIVE_FILE="chr${chromosome}.fa.gz"
  request="ftp.ensembl.org/pub/grch37/current/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.chromosome.$chromosome.fa.gz"
  curl -s "$request" > "$ARCHIVE_FILE"

  # Replace full header with just '>chr7'
  gunzip -c "$ARCHIVE_FILE" | \
    awk -v chr="chr${chromosome}" 'BEGIN {OFS = "\n"} /^>/ {print ">" chr} !/^>/ {print}' > "$REFERENCE_FULL"

  samtools faidx "$REFERENCE_FULL"
fi

# Extract region of interest
samtools faidx "$REFERENCE_FULL" "chr${chromosome}:${start_position}-${end_position}" > "$REFERENCE_FILE"

# Ensure header is plain chr7 (already done)
sed -i "1s/.*/>chr${chromosome}/" "$REFERENCE_FILE"

samtools faidx "$REFERENCE_FILE"
bwa index "$REFERENCE_FILE"

# Output directory
OUTDIR="./$SAMPLE_NAME"
mkdir -p "$OUTDIR"

# Download reads
wget -c "$URL1" -O "$OUTDIR/${SAMPLE_NAME}_R1.fastq.gz"
wget -c "$URL2" -O "$OUTDIR/${SAMPLE_NAME}_R2.fastq.gz"

# Run fastp
fastp -i "$OUTDIR/${SAMPLE_NAME}_R1.fastq.gz" -I "$OUTDIR/${SAMPLE_NAME}_R2.fastq.gz" \
      -o "$OUTDIR/trimmed_R1.fastq.gz" -O "$OUTDIR/trimmed_R2.fastq.gz" \
      --thread $threads --detect_adapter_for_pe

# Align with BWA
bwa mem -t $threads -R "@RG\tID:$SAMPLE_NAME\tSM:$SAMPLE_NAME\tPL:ILLUMINA" "$REFERENCE_FILE" \
    "$OUTDIR/trimmed_R1.fastq.gz" "$OUTDIR/trimmed_R2.fastq.gz" | \
    samtools view -b -o "$OUTDIR/aligned.bam" -

# Sort and index
samtools sort -o "$OUTDIR/aligned.sorted.bam" "$OUTDIR/aligned.bam"
samtools index "$OUTDIR/aligned.sorted.bam"

# Call variants
freebayes -f "$REFERENCE_FILE" -p "$PLOIDY_NUM" \
  --min-base-quality 20 \
  --min-mapping-quality 30 \
  --use-best-n-alleles 2 \
  "$OUTDIR/aligned.sorted.bam" > "${SAMPLE_NAME}.raw.vcf"

# Filter variants
bcftools filter -i 'QUAL > 100 && INFO/DP > 30' -Ov -o "${SAMPLE_NAME}.filtered.vcf" "${SAMPLE_NAME}.raw.vcf"
bgzip "${SAMPLE_NAME}.filtered.vcf"
bcftools index "${SAMPLE_NAME}.filtered.vcf.gz"

# Fix CHROM field to match FASTA header ('chr7')
echo -e "chr${chromosome}:${start_position}-${end_position}\tchr${chromosome}" > rename_chrs.txt
bcftools annotate --rename-chrs rename_chrs.txt "${SAMPLE_NAME}.filtered.vcf.gz" -Oz -o "${SAMPLE_NAME}.renamed.vcf.gz"
bcftools index "${SAMPLE_NAME}.renamed.vcf.gz"

# Call consensus
cat "$REFERENCE_FILE" | bcftools consensus \
  --haplotype R \
  -s "$SAMPLE_NAME" \
  "${SAMPLE_NAME}.renamed.vcf.gz" > "${SAMPLE_NAME}.consensus.fa"

# Optional: clean up
rm -rf "$OUTDIR"
rm rename_chrs.txt

echo "✅ Finished. Output FASTA: ${SAMPLE_NAME}.consensus.fa"

