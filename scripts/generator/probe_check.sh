# carefully: can not accept wildcard.

# contig table: contig_name \t genome_name
# prepare:
# cat group_* | sed 's/.*>//g' | sed 's/ /\t/' | sed 's/ /_/' | sed 's/ .*//g' > contigs
# cat group_* | sed 's/:>/\t/g' | sed 's/ .*//g' | sed 's/.*\///g' > contigs

# read args
#shopt -s extglob

while getopts t:o:p:d:r:m:e:i:a:b:c:*: flag
do
    case "${flag}" in
        t) true_file=${OPTARG};;
        o) output=${OPTARG};;
        p) path_to_filter_py=${OPTARG};;
        d) contig_table=${OPTARG};;
        r) tmp_dir=${OPTARG};;
        m) max_out=${OPTARG};;
        e) max_mismatch=${OPTARG};;
        i) min_ident=${OPTARG};;
        a) multimap_max=${OPTARG};;
        b) negative_max=${OPTARG};;
        c) min_seq=${OPTARG};;
        d) max_seq=${OPTARG};;
    esac
done
shift $(( OPTIND - 1 ))

# defaults
max_mismatch="${max_mismatch:-5}"
multimap_max="${multimap_max:-0}"
negative_max="${negative_max:-0}"
min_ident="${min_ident:-90}"
min_seq="${min_seq:-50}"
max_seq="${max_seq:-10000}"

# functions
filter(){
    awk '$7 <= "'$max_mismatch'" && $6 >= "'$min_ident'"' $1
}

sort_unique() {
    negmax=$2
    filter $1 |\
    awk '{print $1}' |\
    sort |\
    uniq -c |\
    awk -v negmax="$negmax" '{if ($1 >= negmax) print $2}'
    #-c | awk '{if ($1>2) print $2}'
}

count() {
    awk '{count[$1]++} END {for (word in count) print count[word], word}' $1 |\
    awk '{print $2 "\t" $1}'
}

# make tmp dir
if [ -z $tmp_dir ]; then
    tmp=./.tmp
else
    tmp=$tmp_dir
fi
mkdir -p $tmp || { echo "Failed to create tmp directory"; exit 1; }

echo "Probe check initiated"

# count false mappers
for filter_files in "$@"; do
    bn=$(basename $filter_files)
    echo " - prepare false hits: ${bn%.*}"
    sort_unique $filter_files $negative_max >> $tmp/filter_files.tmp
done

echo " - count occurences"
#filter true mappers
filter $true_file > $tmp/filtered.tmp

#unique hitters
sort $tmp/filtered.tmp > $tmp/sorted_1.tmp
count $tmp/sorted_1.tmp | sort > $tmp/count.tmp
echo " --" $(wc -l < $tmp/count.tmp) probes sucsessfully mapped

#multimapping check
echo " - multimapping check"

# test if both files created successfully
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

# concat
sort $tmp/filter_files.tmp |\
uniq > $tmp/remove.tmp

echo " - filtration"

# python exec
python $path_to_filter_py $tmp/amount_sorted.tmp $tmp/remove.tmp $output $max_out

# cleanup
if [ -z $tmp_dir ]; then
    echo "Clear ----"
    rm -r $tmp
fi