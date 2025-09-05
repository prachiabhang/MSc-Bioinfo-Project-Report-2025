#! /bin/bash

# if there are previous results, remove them
if [ -e ./sequences ]; then
  rm -rf ./sequences
fi

mkdir ./sequences

#Set some sort of stem name to make modifying this easier at later date:
chromosome="7"
start_position="128767485"
end_position="128780794"
orientation="-1"
ploidy="A"
stem_filename="$chromosome-$start_position-$end_position"

# get VCF portion
tabix -h https://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20130502/ALL.chr"$chromosome".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz "$chromosome:$start_position-$end_position" > ./sequences/$stem_filename.vcf
# compress
bgzip ./sequences/$stem_filename.vcf
VCF_FILE=./sequences/$stem_filename.vcf.gz

#index
bcftools index $VCF_FILE

# Download reference sequence from GRCh37 Ensembl endpoint
ARCHIVE_FILE="./sequences/chr${chromosome}.fa.gz"

request="ftp.ensembl.org/pub/grch37/current/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.chromosome.$chromosome.fa.gz"
curl -s "$request" > "$ARCHIVE_FILE"

REFERENCE_FILE="./sequences/chr${chromosome}.fa"
gunzip -c "$ARCHIVE_FILE" > "$REFERENCE_FILE"

REFERENCE_FILE="./sequences/chr${chromosome}.fa"
samtools faidx "$REFERENCE_FILE"

# generate samples
samples_file=./sequences/samples.txt
bcftools query -l $VCF_FILE >> $samples_file

# see if there is actually work to do
echo "🔍 Checking for actual variant genotypes in region $chromosome:$start_position-$end_position..."

variant_count=$(bcftools view -r "$chromosome:$start_position-$end_position" "$VCF_FILE" | grep -v "^#" | wc -l)
echo " - Variants in region: $variant_count"

if [ "$variant_count" -eq 0 ]; then
  echo "⚠️  No variant records found in region — bcftools consensus will produce reference-only output."
  exit
fi

# Extract region from reference once
REFERENCE_REGION="./sequences/${stem_filename}_reference_region.fa"
samtools faidx "$REFERENCE_FILE" "$chromosome:$start_position-$end_position" > "$REFERENCE_REGION"

# Now use this for consensus
for sample in $(cat "$samples_file"); do
  echo " - Sample: $sample"

  if [ "$ploidy" == "A" ]; then
		for hap in 1 2; do
		  echo " - - Haplotype $hap"
		  bcftools consensus --haplotype $hap -f "$REFERENCE_REGION" -s "$sample" \
		    -o "./sequences/${stem_filename}_${sample}_hap${hap}.fa" "$VCF_FILE"
		    
		  # Modify the FASTA header to match the filename
		  header_name="${stem_filename}_${sample}_hap${hap}"
		  sed -i "1s/^>.*/>${header_name}/" "./sequences/${stem_filename}_${sample}_hap${hap}.fa"
		done
	elif [ "$ploidy" == "X" ]; then
		# Detect ploidy by checking GT fields for this sample in the region
		gt_format=$(bcftools view -r "$chromosome:$start_position-$end_position" -s "$sample" "$VCF_FILE" \
		  | grep -v '^#' | awk -F'\t' '{print $10}' | cut -d ':' -f1 | head -n1)

		if [[ "$gt_format" == *"/"* ]]; then
		  # Diploid genotype — generate both haplotypes
		  echo " - - Diploid X locus for sample $sample (e.g., 0/1)"
		  for hap in 1 2; do
		    echo "   - Haplotype $hap"
		    output="./sequences/${stem_filename}_${sample}_hap${hap}.fa"
		    bcftools consensus --haplotype $hap -f "$REFERENCE_REGION" -s "$sample" -o "$output" "$VCF_FILE"
		    header_name="${stem_filename}_${sample}_hap${hap}"
		    sed -i "1s/^>.*/>${header_name}/" "$output"
		  done
		else
		  # Haploid genotype — only one haplotype
		  echo " - - Haploid X locus for sample $sample (e.g., 0)"
		  output="./sequences/${stem_filename}_${sample}_hap1.fa"
		  bcftools consensus --haplotype 1 -f "$REFERENCE_REGION" -s "$sample" -o "$output" "$VCF_FILE"
		  header_name="${stem_filename}_${sample}_hap1"
		  sed -i "1s/^>.*/>${header_name}/" "$output"
		fi
  fi
done

# tidy up
rm "$REFERENCE_FILE"*
rm "$VCF_FILE"*

