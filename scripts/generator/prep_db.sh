# opts
while getopts n:c:t: flag
do
    case "${flag}" in
        n) database_name=${OPTARG};;
        c) contig_name=${OPTARG};;
        t) tmp=${OPTARG};;
    esac
done
shift $(( OPTIND - 1 ))

# mktmp
if [ -z $tmp ]; then
    tmp=./.tmp
fi
mkdir $tmp -p

#mkdb
dbdir=$(echo $database_name | sed "s/\/.*/\//g")
mkdir $dbdir -p

#merge fasta in tmp
for fname in "$@"; do
    clear_name=$(echo $fname | sed "s/.*\///g" | sed "s/\..*//g")
    #unzip gzipped
    gzip -dq $fname
    fname2=$(echo $fname | sed "s/\.gz$//g")
    
    echo "Add $clear_name"

    sed 's/ .*//g' $fname2 >> $tmp/$clear_name.fa
    #sed -e 's/[\.,]/_/g' 
    
    grep "^>" $tmp/$clear_name.fa |\
    sed "s/^>/${clear_name}\t/" >> $contig_name
done
cat $tmp/*.fa > $tmp/merged.fa

# create base
makeblastdb -in $tmp/merged.fa -out $database_name -dbtype nucl

rm -r $tmp