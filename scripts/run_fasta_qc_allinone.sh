#!/usr/bin/env bash
set -euo pipefail

# ============================================
# All-in-one FASTA concatenation + QC + Report
# ============================================
# Usage:
#   bash run_fasta_qc_allinone.sh "<path_to_fasta_dir>" [aligned_fasta]
#
# Env option:
#   MAX_N_PERC=20   # sequences with Ns% > threshold are flagged (default 20)
#
# Outputs: combined_sequences.fasta, qc_out/*, and a printed summary.

RAW_DIR="${1:-}"
ALN_FASTA="${2:-}"                     # optional aligned FASTA
COMBINED="combined_sequences.fasta"
MAX_N_PERC="${MAX_N_PERC:-20}"

if [[ -z "$RAW_DIR" || ! -d "$RAW_DIR" ]]; then
  echo "ERROR: Directory not found. Provide a valid path." >&2
  echo 'Tip: quote paths that contain spaces or "&".' >&2
  exit 1
fi

echo ">>> Scanning FASTA files in:"
echo "    $RAW_DIR"
TMP_LIST="$(mktemp)"
# Recursively find fasta-like files (case-insensitive), including gzipped
find "$RAW_DIR" -type f \
  \( -iname "*.fa" -o -iname "*.fasta" -o -iname "*.fas" -o -iname "*.fna" \
     -o -iname "*.fa.gz" -o -iname "*.fasta.gz" -o -iname "*.fas.gz" -o -iname "*.fna.gz" \) \
  -print0 | sort -z > "$TMP_LIST"

COUNT=$(tr '\0' '\n' < "$TMP_LIST" | wc -l)
echo ">>> Files detected: $COUNT"
if [[ "$COUNT" -eq 0 ]]; then
  echo "ERROR: No FASTA files found (.fa/.fasta/.fas/.fna or gzipped variants)." >&2
  rm -f "$TMP_LIST"
  exit 1
fi

echo ">>> Concatenating into $COMBINED ..."
: > "$COMBINED"
while IFS= read -r -d '' f; do
  if [[ "$f" =~ \.gz$ ]]; then
    gzip -cd -- "$f" >> "$COMBINED"
  else
    cat -- "$f" >> "$COMBINED"
  fi
done < "$TMP_LIST"
rm -f "$TMP_LIST"

SEQ_N=$(grep -c '^>' "$COMBINED" || true)
echo ">>> Sequences in combined FASTA: $SEQ_N"
echo

# -------------------------
# QC ENGINE (inline)
# -------------------------
mkdir -p qc_out

# 1) Normalize FASTA to one-seq-per-line
awk '
  /^>/ { if (NR>1) print seq; print; seq=""; next }
  { gsub(/\r/,""); seq = seq $0 }
  END { if (seq) print seq }
' "$COMBINED" > qc_out/oneLine.fasta

# 2) Per-sequence stats: LEN, Ns, Ns%, GC% (excluding N/-), non-ACGTN
awk '
  NR%2==1 { hdr=$0; sub(/^>/,"",hdr); next }
  NR%2==0 {
    seq=$0
    len=length(seq)
    nN=gsub(/N/,"&",seq)+gsub(/n/,"&",seq)

    s=seq; gsub(/-/,"",s)          # remove gaps
    t=s;   gsub(/[Nn]/,"",t)       # remove Ns for GC calc
    len_noN=length(t)
    gc=gsub(/[Gg]/,"&",t)+gsub(/[Cc]/,"&",t)

    tmp=seq; gsub(/[-ACGTNacgtn]/,"",tmp)  # unexpected symbols
    bad=length(tmp)

    nperc=(len? nN/len*100:0)
    gcperc=(len_noN? gc/len_noN*100:0)

    printf("%s\tLEN=%d\tNs=%d\tNs_perc=%.4f\tGC_perc=%.4f\tNonACGTN=%d\n",
           hdr, len, nN, nperc, gcperc, bad)
  }
' qc_out/oneLine.fasta > qc_out/fasta_qc.tsv

# 3) Flag sequences with Ns% over threshold
awk -v thr="$MAX_N_PERC" -F'\t' '
BEGIN{OFS="\t"; print "ID","Ns_perc","FLAG"}
{
  id=$1
  for(i=1;i<=NF;i++){
    if($i ~ /^Ns_perc=/){ split($i,a,"="); n=a[2]+0 }
  }
  flag=(n>thr? "EXCLUDE_N" : "OK")
  print id, n, flag
}' qc_out/fasta_qc.tsv > qc_out/flag_N.tsv

# 4) Detect exact duplicate sequences (identical strings)
awk '
  NR%2==1{h=$0; sub(/^>/,"",h); next}
  NR%2==0{
    if(!seen[$0]++){ names[$0]=names[$0] (names[$0]? "," : "") h }
  }
  END{ for(k in seen) if(seen[k]>1) printf("%d\t%s\n", seen[k], names[k]) }
' qc_out/oneLine.fasta > qc_out/duplicates.txt

# 5) Summary files
{
  echo "# Summary"
  echo "Total sequences:  $(grep -c '^>' qc_out/oneLine.fasta)"
  echo -n "Mean length:      "
  awk -F'\t' '{for(i=1;i<=NF;i++) if($i ~ /^LEN=/){split($i,a,"="); s+=a[2]; n++}}
              END{if(n) printf("%.1f\n",s/n); else print 0}' qc_out/fasta_qc.tsv
  echo -n "Median Ns%:       "
  awk -F'\t' '{for(i=1;i<=NF;i++) if($i ~ /^Ns_perc=/){split($i,a,"="); print a[2]}}' qc_out/fasta_qc.tsv \
    | sort -g | awk '{a[NR]=$1} END{if(NR%2){print a[(NR+1)/2]} else {printf("%.4f\n",(a[NR/2]+a[NR/2+1])/2)}}'
  echo -n "Median GC%:       "
  awk -F'\t' '{for(i=1;i<=NF;i++) if($i ~ /^GC_perc=/){split($i,b,"="); print b[2]}}' qc_out/fasta_qc.tsv \
    | sort -g | awk '{a[NR]=$1} END{if(NR%2){print a[(NR+1)/2]} else {printf("%.4f\n",(a[NR/2]+a[NR/2+1])/2)}}'
  echo "Duplicates (rows): $(wc -l < qc_out/duplicates.txt)"
} > qc_out/summary.txt

# -------------------------
# REPORT: print the exact items requested
# -------------------------
echo "------------------ QC REPORT ------------------"

# Median GC%
MED_GC=$(awk -F'\t' '{for(i=1;i<=NF;i++) if($i~/^GC_perc=/){split($i,a,"="); print a[2]}}' qc_out/fasta_qc.tsv \
  | sort -g | awk '{a[NR]=$1} END{if(NR%2){printf("%.4f",a[(NR+1)/2])} else {printf("%.4f",(a[NR/2]+a[NR/2+1])/2)}}')
echo "Median GC%:"
echo "$MED_GC"

# Median Ns%
MED_NS=$(awk -F'\t' '{for(i=1;i<=NF;i++) if($i~/^Ns_perc=/){split($i,a,"="); print a[2]}}' qc_out/fasta_qc.tsv \
  | sort -g | awk '{a[NR]=$1} END{if(NR%2){printf("%.4f",a[(NR+1)/2])} else {printf("%.4f",(a[NR/2]+a[NR/2+1])/2)}}')
echo
echo "Median Ns%:"
echo "$MED_NS"

# Unique lengths
UNIQUE_LENS=$(awk -F'\t' '{for(i=1;i<=NF;i++) if($i~/^LEN=/){split($i,a,"="); print a[2]}}' qc_out/fasta_qc.tsv | sort -u | tr '\n' ' ')
echo
echo "Unique sequence lengths (bp):"
echo "$UNIQUE_LENS"

# How many sequences were flagged for high Ns?
NFLAG=$(awk 'NR>1 && $3=="EXCLUDE_N"{c++} END{print c+0}' qc_out/flag_N.tsv)
echo
echo "Sequences flagged for high Ns (> ${MAX_N_PERC}%):"
echo "$NFLAG"

# List the first few flagged IDs
echo
echo "First few flagged IDs (if any):"
awk 'NR>1 && $3=="EXCLUDE_N"{print $1}' qc_out/flag_N.tsv | head || true

# How many duplicate sets?
DUPSETS=$(wc -l < qc_out/duplicates.txt)
echo
echo "Duplicate sets (rows in duplicates.txt):"
echo "$DUPSETS"

# Peek at per-sequence stats
echo
echo "Head of fasta_qc.tsv:"
head qc_out/fasta_qc.tsv

# One-line Results sentence
# Compute a representative length (~median), and detect a short outlier, if any
MED_LEN=$(awk -F'\t' '{for(i=1;i<=NF;i++) if($i~/^LEN=/){split($i,a,"="); print a[2]}}' qc_out/fasta_qc.tsv \
  | sort -g | awk '{a[NR]=$1} END{if(NR%2){printf("%.0f",a[(NR+1)/2])} else {printf("%.0f",(a[NR/2]+a[NR/2+1])/2)}}')
MIN_LEN=$(awk -F'\t' '{for(i=1;i<=NF;i++) if($i~/^LEN=/){split($i,a,"="); print a[2]}}' qc_out/fasta_qc.tsv | sort -g | head -1)

# Round GC and Ns to one decimal for the sentence
MED_GC_1D=$(awk -v x="$MED_GC" 'BEGIN{printf("%.1f", x+0)}')
MED_NS_1D=$(awk -v x="$MED_NS" 'BEGIN{printf("%.1f", x+0)}')

# Decide if we mention a truncated sequence (min < 75% of median)
SHORT_NOTE=""
awk -v min="$MIN_LEN" -v med="$MED_LEN" 'BEGIN{ if(min < 0.75*med) print "yes"; else print "no"; }' | grep -q yes && \
  SHORT_NOTE=" , except for one truncated sequence (~${MIN_LEN} bp)"

# Compose sentence
echo
echo "One-line Results (paste-ready):"
echo "Across ${SEQ_N} reconstructed haplotypes, the median GC content was ${MED_GC_1D}%, with ${MED_NS_1D}% unresolved bases (Ns), and sequence lengths were consistent at ~${MED_LEN} bp${SHORT_NOTE}."

echo "------------------------------------------------"
echo
echo ">>> QC complete."
echo "    Combined FASTA:        $COMBINED"
echo "    QC outputs directory:  qc_out/"
echo "    Key files: summary.txt, fasta_qc.tsv, flag_N.tsv, duplicates.txt"
[[ -n "$ALN_FASTA" ]] && echo "    Alignment QC: aln_gap_per_seq.tsv, aln_variable_sites.txt"
