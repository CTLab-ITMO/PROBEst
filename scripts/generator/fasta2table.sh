while getopts i:o: flag
do
    case "${flag}" in
        i) input=${OPTARG};;
        o) output=${OPTARG};;
    esac
done

awk 'BEGIN{RS=">"}{print ">"$1"\t"$2;}' $input |\
    sort |\
    sed 's/>\t.*//g' |\
    sed 's/\t/\n/g' |\
    awk '!/^[[:space:]]*$/' > $output