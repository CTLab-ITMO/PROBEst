#!/bin/bash
# Script: probe_check.sh
# Description: This script filters and validates probes based on their mapping characteristics.
#              It checks for mismatches, identity thresholds, and multimapping behavior, and
#              outputs a filtered list of probes. The script is designed to handle large datasets
#              and supports customizable thresholds for filtering.
#
# Usage: ./probe_check.sh -t <true_file> -o <output> -p <path_to_filter_py> -d <contig_table> [OPTIONS] <filter_files>
#
# Required Arguments:
#   -t <true_file>:       File containing true mappings (e.g., BLAST or alignment results).
#   -o <output>:          Output file for filtered probes.
#   -p <path_to_filter_py>: Path to the Python script for final filtering.
#   -d <contig_table>:    File mapping contigs to genome names.
#   <filter_files>:       Files containing false mappings (e.g., BLAST or alignment results for negative controls).
#
# Optional Arguments:
#   -r <tmp_dir>:         Temporary directory for intermediate files (default: ./.tmp).
#   -m <max_out>:         Maximum number of probes to output (default: no limit).
#   -e <max_mismatch>:    Maximum allowed mismatches (default: 5).
#   -i <min_ident>:       Minimum identity percentage (default: 90).
#   -a <multimap_max>:    Maximum allowed multimapping occurrences (default: 0).
#   -b <negative_max>:    Maximum allowed occurrences in negative controls (default: 0).
#   -c <min_seq>:         Minimum sequence length (default: 50).
#   -d <max_seq>:         Maximum sequence length (default: 10000).
#
# Example:
#   ./probe_check.sh -t true_mappings.txt -o filtered_probes.txt -p filter.py -d contig_table.txt -e 3 -i 95 false_mappings_1.txt false_mappings_2.txt

# Parse command-line options
while getopts t:o:p:d:r:m:e:i:a:b:c:*: flag
do
    case "${flag}" in
        t) true_file=${OPTARG};;       # File with true mappings
        o) output=${OPTARG};;          # Output file for filtered probes
        p) path_to_filter_py=${OPTARG};; # Path to Python filtering script
        d) contig_table=${OPTARG};;    # Contig-to-genome mapping table
        r) tmp_dir=${OPTARG};;         # Temporary directory
        m) max_out=${OPTARG};;         # Maximum number of probes to output
        e) max_mismatch=${OPTARG};;    # Maximum allowed mismatches
        i) min_ident=${OPTARG};;       # Minimum identity percentage
        a) multimap_max=${OPTARG};;    # Maximum allowed multimapping occurrences
        b) negative_max=${OPTARG};;    # Maximum allowed occurrences in negative controls
        c) min_seq=${OPTARG};;         # Minimum sequence length
        d) max_seq=${OPTARG};;         # Maximum sequence length
    esac
done
shift $(( OPTIND - 1 ))

# Set default values for optional arguments
max_mismatch="${max_mismatch:-5}"
multimap_max="${multimap_max:-0}"
negative_max="${negative_max:-0}"
min_ident="${min_ident:-90}"
min_seq="${min_seq:-50}"
max_seq="${max_seq:-10000}"

# Functions
filter() {
    # Filter mappings based on mismatches and identity
    awk '$7 <= "'$max_mismatch'" && $6 >= "'$min_ident'"' $1
}

sort_unique() {
    # Sort and filter unique mappings
    negmax=$2
    filter $1 |\
    awk '{print $1}' |\
    sort |\
    uniq -c |\
    awk -v negmax="$negmax" '{if ($1 >= negmax) print $2}'
}

count() {
    # Count occurrences of each probe
    awk '{count[$1]++} END {for (word in count) print count[word], word}' $1 |\
    awk '{print $2 "\t" $1}'
}

# Create temporary directory
if [ -z $tmp_dir ]; then
    tmp=./.tmp
else
    tmp=$tmp_dir
fi
mkdir -p $tmp || { echo "Failed to create tmp directory"; exit 1; }

echo "Probe check initiated"

# Count false mappings
for filter_files in "$@"; do
    bn=$(basename $filter_files)
    echo " - prepare false hits: ${bn%.*}"
    sort_unique $filter_files $negative_max >> $tmp/filter_files.tmp
done

echo " - count occurrences"
# Filter true mappings
filter $true_file > $tmp/filtered.tmp

# Identify unique probes
sort $tmp/filtered.tmp > $tmp/sorted_1.tmp
count $tmp/sorted_1.tmp | sort > $tmp/count.tmp
echo " --" $(wc -l < $tmp/count.tmp) "probes successfully mapped"

# Multimapping check
echo " - multimapping check"
if [ -s "$tmp/sorted_1.tmp" ] && [ -s "$tmp/count.tmp" ]; then
    join $tmp/sorted_1.tmp $tmp/count.tmp | sort -k8 -nr > $tmp/amount_sorted.tmp
else
    echo -e "Problem with sort or join.\nOne of the files is empty or does not exist."
    exit 1
fi

sort $tmp/filtered.tmp -k2 > $tmp/sorted_2.tmp
sort $contig_table -k2 > $tmp/sorted_contig.tmp
join -1 2 -2 2 $tmp/sorted_contig.tmp $tmp/sorted_2.tmp |\
awk '{print $2 "---" $3}' > $tmp/sorted_with_genomes.tmp
count $tmp/sorted_with_genomes.tmp > $tmp/sorted_counts.tmp
sed 's/---/\t/g' $tmp/sorted_counts.tmp |\
awk '{if ($3 > "'$multimap_max'") print $2}' >> $tmp/filter_files.tmp

# Concatenate and deduplicate
sort $tmp/filter_files.tmp |\
uniq > $tmp/remove.tmp

echo " - filtration"

# Execute Python script for final filtering
python $path_to_filter_py $tmp/amount_sorted.tmp $tmp/remove.tmp $output $max_out

# Clean up temporary files
if [ -z $tmp_dir ]; then
    echo "Clear ----"
    rm -r $tmp
fi