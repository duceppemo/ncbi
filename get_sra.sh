#!/bin/bash

###############################################################
#
#   This script will parallel download the SRA Illumina
#   paired-end samples matching the query string.
#
#   Define query string and output folder before running using
#   -q and -o options, respectively. Options are mandatory
#
###############################################################


function print_help()
{
    echo "\
Usage: bash get_sra.sh [-q \"Some organism\" | -l SRA.list] [-o output_folder] 

Mandatory flags:

    -q    Query text (e.g. \"Listeria monocytogenes\").

    -l    List of SRA numbers. One per line.

    -o    Output directory"

    exit 1
}


#http://wiki.bash-hackers.org/howto/getopts_tutorial
query=""
out=""
list=""

options=':q:o:l:h'

while getopts "$options" opt; do
    case "$opt" in
        q)
            query="$OPTARG"
            ;;
        o)
            export out="$OPTARG"
            ;;
        l)
            export list="$OPTARG"
            ;;
        h)
            print_help
            ;;
        \?)
            echo "Invalid option: -"$OPTARG"" >&2
            print_help
            ;;
        :)
            echo "Option -"$OPTARG" requires an argument." >&2
            print_help
            ;;
    esac
done

shift $(($OPTIND - 1))


#Checks

# Check if only one of "-q" or "-l" has been used
if [[ -n "$query" ]] && [[ -n "$list" ]]; then
    echo "You can't use both \"-q\" and \"-l\" in the same command"
    print_help
fi

# Exit if argument is missing
if [[ -n "$list" ]] && [[ ! -s "$list" ]]; then
    echo "The list file provided is empty"
    print_help
fi

# Check if the list file seem valid
if [[ -n "$list" ]]; then
    is_bad=$(head -n 1 "$list" | grep -vE "[a-zA-z0-9]")
    if [[ -n "$is_bad" ]]; then
        echo "The list file provided has an invalid format"
        print_help
    fi
fi

# Check if output folder provided
if [[ -z "$out" ]]; then
    echo "No output directory provided"
    print_help
fi


# query="Campylobacter fetus subsp. fetus"
# export out="/media/3tb_hdd/data/campylobacter/raw_reads/sra/cff"


# Create output directory if does not exist
[ -d "$out" ] || mkdir -p "$out"

if [[ "$query" ]]; then

    # Get accession numbers and store them in array
    esearch -db sra -query "$query" \
        | efetch -format runinfo \
        | awk 'NR>1' \
        | grep "GENOMIC" \
        | grep "PAIRED" \
        | grep "ILLUMINA" \
        | cut -d "," -f 1 \
        > "${out}"/SraAccList.txt
else
    cp "$list" "${out}"/SraAccList.txt
fi

# Statisticts
nResults=$(cat "${out}"/SraAccList.txt | wc -l)
if [[ "$query" ]]; then
    echo -e ""$nResults" SRA paired-end Illumina DNA samples were found for \""$query"\". Downloading..."
else
    echo -e "Downloading the provided "$nResults" SRA accessions..."
fi

# Get data from ncbi:
# Split R1 and R2 in two files
# Compress the data
# save to output directory
# Limit number of parallel download du to limited bandwith (10MB/s) ??? -P 24
cat "${out}"/SraAccList.txt \
    | parallel  -P 24 \
                --bar \
                --env out \
                --delay 0.3 \
                'fastq-dump --gzip --split-files -O $out {}'

# Renane files for pipeline compliance
rename 's/_1/_R1/' "${out}"/*.fastq.gz
rename 's/_2/_R2/' "${out}"/*.fastq.gz

# Replace 3rd line of every sequence in fastq file by "+"
# Reduces file size
find "$out" -maxdepth 1 -type f -name "*.fastq.gz" |
    parallel --bar 'zcat {} | sed "3~4 s/.*/+/" | pigz > {}.tmp'
find "$out" -maxdepth 1 -type f -name "*.fastq.gz.tmp" |
    parallel --bar 'mv {} {.}'

