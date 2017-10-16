#!/bin/bash


# Script description:

# This script will download ".fna" ".ffn" and ".gff", if availables, from NCBI RefSeq for a specific database.
# Assembly level and name filtering are also available as options.


#################
#               #
#    Options    #
#               #
#################


function displayHelp()
{
    echo -e "\

    Usage: refseqDownloader.sh -t db_type -o output_folder [-n organism(s)] [-l assemmly_level(s)]

    Mandatory flags:

        -t          RefSeq database type (bacteria, archaea, viral, fungi, protozoa, vertebrate, plant, human).

        -o          Output folder. Where downloaded files will be located.

    Optional flags:

        -h          Print this help message.

        -n          Organism name(s) to keep.
                    Comma separated with no space if multiple (e.g. "Listeria monocytogenes,Campylobacter,Salmonella").
                    Not using this option will download all species of the selected database type.

        -l          Assembly level (Chromosome, Complete Genome, Contig, Scaffold).
                    Comma separated with no space if multiple (e.g. "Complete Genome,Scaffold").
                    Not using this option will download all assembly levels of the selected database type.

        -r          Rename files according to description.

        -u          Uncompress all downloaded files.

        -j          Just fna.
    "
}


#Colored error message
BLUE='\033[1;34m'
NC='\033[0m' # No Color

db_type=''
export output=''
name=''
level=''
uncompress=0
rename=0
just_fna=0
help=''


# If the very first character of the option-string is a : (colon),
# which would normally be nonsense because there's no option letter preceding it,
# getopts switches to "silent error reporting mode"
# When you want getopts to expect an argument for an option, just place a : (colon) after the proper option flag
options=':t:o:n:l:ujrh'

while getopts "$options" opt; do
    case "$opt" in
        t)
            db_type="$OPTARG"
            ;;
        o)
            export output="$OPTARG"
            ;;
        n)
            name="$OPTARG"
            ;;
        l)
            level="$OPTARG"
            ;;
        h)
            displayHelp
            exit 1
            ;;
        r)
            rename=1
            ;;
        u)
            uncompress=1
            ;;
        j)
            just_fna=1
            ;;
        \?)
            printf ""${BLUE}"Invalid option: -"$OPTARG"\n\n"${NC}"" >&2
            # echo "Invalid option: -"$OPTARG"" >&2
            displayHelp
            exit 1
            ;;
        :)
            printf ""${BLUE}"Option -"$OPTARG" requires an argument.\n\n"${NC}"" >&2
            # echo "Option -"$OPTARG" requires an argument." >&2
            displayHelp
            exit 1
            ;;
    esac
done

shift $(($OPTIND - 1))

# Exit if option flags b or v are missing
if [[ -z "$db_type" ]] || [[ -z "$output" ]]; then
    echo "Both \"-t\" and \"-o\" options are mandatory and require arguments"
    displayHelp
    exit 1
fi


###################
#                 #
#   Downloading   #
#                 #
###################

summary_url=''

case "$db_type" in
    bacteria)   summary_url="ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt";;
    archaea)    summary_url="ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/archaea/assembly_summary.txt";;
    viral)      summary_url="ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/viral/assembly_summary.txt";;
    fungi)      summary_url="ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/fungi/assembly_summary.txt";;
    protozoa)   summary_url="ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/protozoa/assembly_summary.txt";;
    vertebrate) summary_url="ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/assembly_summary.txt";;
    plant)      summary_url="ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/plant/assembly_summary.txt";;
    human)      summary_url="ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/assembly_summary.txt";;
esac

# Check everything went fine
if [ -z "$summary_url" ]; then  # if not
    echo "Wrong database type!"
    displayHelp
    exit 1
else
    #create output directory
    [ -d "$output" ] || mkdir -p "$output"

    #create a time stamp to track the time of download
    downloadDate=$(date +%Y-%m-%d)
    summary_file=""${output}"/"${downloadDate}"_assembly-summary_"${db_type}".txt"

    if [ ! -s "$summary_file" ]; then
        #Download the refseq assembly summary file
        echo -e "Downloading \"assembly_summary.txt\" for "$db_type" from NCBI..."
        curl -# "$summary_url" > "$summary_file"
    fi
fi


#################
#               #
#   Filtering   #
#               #
#################


levels=''

#make sure "Complete Genome" is treated a a single array element
levels=($(echo "$level" | tr " " "_" | tr "," " "))

#check if correct assembly level term used
for l in "${levels[@]}"; do
    assembly_level=''
    l=$(tr "_" " " <<< "$l")

    case "$l" in  # indirect variable expansion
        Chromosome)         assembly_level="$level";;
        'Complete Genome')  assembly_level="$level";;
        Contig)             assembly_level="$level";;
        Scaffold)           assembly_level="$level";;
    esac

    if [ -z "$assembly_level" ]; then
        echo "Invalid assembly level term \""$l"\""
        displayHelp
        exit 1
    fi
done

declare -a level_filter=()

#Convert to awk search format
if [ "$level" ]; then
    if [ "${#levels[@]}" -eq 1 ]; then
        f=("\$12==\"$(echo "$levels" | tr "_" " ")\"")
    elif [ "${#levels[@]}" -gt 1 ]; then
        for l in "${levels[@]}"; do
            level_filter+=("\$12==\"$(echo "$l" | tr "_" " ")\"")
        done
        f=$(echo "${level_filter[@]}" | sed 's/" $/" || $/g')
    fi
fi

#if entered names to keep
if [ -n "$name" ]; then
    name_filter=$(echo "\""$name"\"" | tr "," "|")
    # echo "$name_filter"
fi

# File with download paths
if [ -n "$level" ] && [ -n "$name" ]; then
    echo "Downloading all "$level" sequences for "$name"..."
    cmd="cat "$summary_file" | grep -E "$name_filter" | awk -F \"\\t\" '"$f" && \$11==\"latest\" {print \$20}' > "${output}"/completeGenomePaths.txt"
elif [ -n "$level" ] && [ -z "$name" ]; then
    echo "Downloading all "$level" sequences for "$db_type"..."
    cmd="cat "$summary_file" | awk -F \"\\t\" '"$f" && \$11==\"latest\" {print \$20}' > "${output}"/completeGenomePaths.txt"
elif [ -z "$level" ] && [ -n "$name" ]; then
    echo "Downloading all available sequences for "$name"..."
    cmd="cat "$summary_file" | grep -E "$name_filter" | awk -F \"\\t\" '\$11==\"latest\" {print \$20}' > "${output}"/completeGenomePaths.txt"
else
    echo "Downloading all available sequences for "$db_type"..."
    cmd="cat "$summary_file" | awk -F \"\\t\" '\$11==\"latest\" {print \$20}' > "${output}"/completeGenomePaths.txt"
fi

# Add a check
if [ "$cmd" ]; then
    eval "$cmd"
else
    echo "Something went very wrong"
    exit 1
fi

# Create ouput subfolders
[ -d "${output}"/fna ] || mkdir -p "${output}"/fna
if [ $just_fna -eq 0 ]; then
    [ -d "${output}"/ffn ] || mkdir -p "${output}"/ffn
    [ -d "${output}"/gff ] || mkdir -p "${output}"/gff
fi

# Fucntion to parallel download
function download()
{
    fileName=$(basename "$1")

    #Determine which files to download
    fna="/"${fileName}"_genomic.fna.gz"
    #Create the full file path
    fnaDownload=$(echo "$1" | awk -v var="$fna" '{print $0var}')
    #check if already downloaded. if not, download.
    [ -s "${output}"/"${fileName}".fna.gz ] || curl -s "$fnaDownload" > "${output}"/fna/"${fileName}".fna.gz

    if [ $just_fna -eq 0 ]; then
        ffn="/"${fileName}"_cds_from_genomic.fna.gz"
        gff="/"${fileName}"_genomic.gff.gz"

        gffDownload=$(echo "$1" | awk -v var="$gff" '{print $0var}')
        ffnDownload=$(echo "$1" | awk -v var="$ffn" '{print $0var}')

        [ -s "${output}"/"${fileName}".ffn.gz ] || curl -s "$ffnDownload" > "${output}"/ffn/"${fileName}".ffn.gz
        [ -s "${output}"/"${fileName}".gff.gz ] || curl -s "$gffDownload" > "${output}"/gff/"${fileName}".gff.gz
    fi

    #rename the file with the organism name from the fasta header
    if [ "$rename" -eq 1 ]; then
        org=$(zcat "${output}"/fna/"${fileName}".fna.gz \
            | head -n 1 \
            | cut -d " " -f 2- \
            | cut -d "," -f 1 \
            | sed   -e "s%'\|(\|)\|\.%%g" \
                    -e 's/ genome//' \
                    -e 's/ chromosome//' \
                    -e 's/ complete genome//' \
                    -e 's/ chromosome sequence//' \
            | tr " " "_" \
            | tr "/" "_")

        mv "${output}"/fna/"${fileName}".fna.gz "${output}"/fna/"${org}".fna.gz
        if [ $just_fna -eq 0 ]; then
            mv "${output}"/ffn/"${fileName}".ffn.gz "${output}"/ffn/"${org}".ffn.gz
            mv "${output}"/gff/"${fileName}".gff.gz "${output}"/gff/"${org}".gff.gz
        fi
    fi
}

#make function available to parallel
export -f download  # -f is to export functions

# Download multiple refseq genomes in parallel
cat "${output}"/completeGenomePaths.txt \
    | parallel --bar \
        --delay 0.3 \
        --env download \
        --env output \
        'download {}'

# TODO -> create folder for each Genus or Species?

#uncompress?
if [ "$uncompress" -eq 1 ]; then
    find "$output" -type f -name "*.gz" \
        | parallel --bar \
            'pigz -d {}'
fi
