#!/bin/bash

version="0.1"

# Script description:

# This script will download assemblies in ".fna" ".ffn", ".faa" and ".gff" format, if available, from NCBI.
# Assembly level and name filtering are also available as options.


#################
#               #
#    Options    #
#               #
#################


function displayHelp()
{
    echo -e "\

    Usage: refseqDownloader.sh -o output_folder [-q \"organism\"] [-l \"assemmly_level\"] -r -u -n 4

    Mandatory flags:

        -o          Output folder. Where downloaded files will be located.

    Optional flags:

        -h          Print this help message.

        -q          Query. Organism name(s) to keep.
                    Comma separated with no space if multiple (e.g. "Listeria monocytogenes,Campylobacter,Salmonella").
                    Not using this option will download all species of the selected database type.

        -l          Assembly level (Chromosome, Complete Genome, Contig, Scaffold).
                    Comma separated with no space if multiple (e.g. "Complete Genome,Scaffold").
                    Not using this option will download all assembly levels of the selected database type.

        -f          File type(s) to download(fna, ffn, gff, faa).
                    Comma separated with no space if multiple (e.g. "fna,ffn,gff,faa,gbk").
                    Not using this option will download all three file types.

        -t          Download only type (representative) strains.

        -r          Rename files according to description in fna file.
                    "fna" file type must be selected to enable renaming.

        -u          Uncompress all downloaded files.

        -n          Number of CPUs to use for parallel download.
                    Default is max numbers of CPUs avaiable.
    "
}


#Colored error message
BLUE='\033[1;34m'
NC='\033[0m' # No Color

# Default values
export output=''
name=''
level=''
export seq_type=""
type=0
uncompress=0
rename=0
export cpu=$(nproc)
help=''

# If the very first character of the option-string is a : (colon),
# which would normally be nonsense because there's no option letter preceding it,
# getopts switches to "silent error reporting mode"
# When you want getopts to expect an argument for an option, just place a : (colon) after the proper option flag
options=':o:n:f:l:utrq:h'

while getopts "$options" opt; do
    case "$opt" in
        o)  export output="$OPTARG";;
        q)  name="$OPTARG";;
        l)  level="$OPTARG";;
        h)  displayHelp
            exit 1
            ;;
        r)  rename=1;;
        t)  type=1;;
        f)  export seq_type="$OPTARG";;
        u)  uncompress=1;;
        n)  export cpu="$OPTARG";;
        \?) printf ""${BLUE}"Invalid option: -"$OPTARG"\n\n"${NC}"" >&2
            # echo "Invalid option: -"$OPTARG"" >&2
            displayHelp
            exit 1;;
        :)  printf ""${BLUE}"Option -"$OPTARG" requires an argument.\n\n"${NC}"" >&2
            # echo "Option -"$OPTARG" requires an argument." >&2
            displayHelp
            exit 1;;
    esac
done

shift $(($OPTIND - 1))

# check mandatary flags and arguments
if [[ -z "$output" ]]; then
    echo "\"-o\" option is mandatory and require an argument"
    displayHelp
    exit 1
fi

# Test file type to get
if [[ -z $(echo "$seq_type" | grep -E "fna|faa|gff|faa|gbk") ]]; then
    echo "Please options are mandatory and require arguments"
    displayHelp
    exit 1
fi

# Check if number of core is exceding number of cores available
if [ "$cpu" ]; then
    if [ "$cpu" -gt "$(nproc)" ]; then
        echo "Number of cores entered ("$cpu") excedes total number of cores available ("$(nproc)")"
        print_help
    fi
fi


###################
#                 #
#   Downloading   #
#                 #
###################

summary_url="https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_genbank.txt"


# Create output directory if not exist
[ -d "$output" ] || mkdir -p "$output"

# Create a time stamp to track the time of download
downloadDate=$(date +%Y-%m-%d)
[ -e "${output}"/"${downloadDate}".txt ] && rm "${output}"/"${downloadDate}".txt
touch "${output}"/"${downloadDate}".txt

#Download the refseq assembly summary file
cd "$output"
wget "$summary_url" -q --show-progress #-P "$output"
summary_file=""${output}"/assembly_summary_genbank.txt"


#################
#               #
#   Filtering   #
#               #
#################


# Make sure "Complete Genome" is treated as a single array element
levels=($(echo "$level" | tr " " "_" | tr "," " "))

# Check if correct assembly level term used
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

# Convert to awk search format
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

# If entered names to keep
export filtered_summary=""${output}"/assembly_summary_filtered.txt"
if [[ -n "$name" ]]; then
    # Create a subset of assembly_summary using the name(s) supplied
    name_filter=$(echo "'"$name"'" | tr "," "|")
    echo "Filtering metadata using: "$name_filter""
    cmd1="cat "$summary_file" | grep -E "$name_filter" > "$filtered_summary""
    eval "$cmd1"
else
    mv "$summary_file" "$filtered_summary"
fi

# Create another subset if representive only was selected
if [[ "$type" -eq 1 ]]; then
    echo "Only keeping representative genomes"
    cmd2="cat $filtered_summary | grep 'representative genome' > ${filtered_summary}.tmp"
    eval "$cmd2"
    mv "${filtered_summary}".tmp "$filtered_summary"
fi

# File with download paths
if [ -n "$level" ] && [ -n "$name" ]; then
    echo "Downloading all "$level" sequences for "$name"..."
    cmd="cat "$filtered_summary" | awk -F \"\\t\" '"$f" && \$11==\"latest\" {print \$20}' > "${output}"/completeGenomePaths.txt"
elif [ -n "$level" ] && [ -z "$name" ]; then
    echo "Downloading all "$level" sequences for "$db_type"..."
    cmd="cat "$filtered_summary" | awk -F \"\\t\" '"$f" && \$11==\"latest\" {print \$20}' > "${output}"/completeGenomePaths.txt"
elif [ -z "$level" ] && [ -n "$name" ]; then
    echo "Downloading all available sequences for "$name"..."
    cmd="cat "$filtered_summary" | awk -F \"\\t\" '\$11==\"latest\" {print \$20}' > "${output}"/completeGenomePaths.txt"
else
    echo "Downloading all available sequences for "$db_type"..."
    cmd="cat "$filtered_summary" | awk -F \"\\t\" '\$11==\"latest\" {print \$20}' > "${output}"/completeGenomePaths.txt"
fi

# Add a check
if [ "$cmd" ]; then
    eval "$cmd"
else
    echo "Something went very wrong"
    exit 1
fi

#file types
seq_types=($(echo "$seq_type" | tr "," " "))  # Convert string to array

# check if correct sequence type terms are being used
for t in "${seq_types[@]}"; do
    ty=''

    case "$t" in  # indirect variable expansion
        fna)    ty="$t"
                [ -d "${output}"/fna ] || mkdir -p "${output}"/fna
                ;;
        ffn)    ty="$t"
                [ -d "${output}"/ffn ] || mkdir -p "${output}"/ffn
                ;;
        gff)    ty="$t"
                [ -d "${output}"/gff ] || mkdir -p "${output}"/gff
                ;;
        faa)    ty="$t"
                [ -d "${output}"/faa ] || mkdir -p "${output}"/faa
                ;;
        gbk)    ty="$t"
                [ -d "${output}"/gbk ] || mkdir -p "${output}"/gbk
                ;;
    esac

    if [ -z "$ty" ]; then
        echo "Invalid sequence type \""$t"\""
        displayHelp
        exit 1
    fi
done

# Fucntion to parallel download
function download()
{
    s_types=($(echo "$seq_type" | tr "," " "))
    fileName=$(basename "$1")

    # Determine which files to download
    # wget -P "${output}"/fna -q --show-progress "$fnaDownload"
    for t in "${s_types[@]}"; do
        case "$t" in  # Indirect variable expansion
            fna)    fna="/"${fileName}"_genomic.fna.gz"
                    fnaDownload=$(echo "$1" | awk -v var="$fna" '{print $0var}')  # Create the full file path
                    [ -s "${output}"/"${fileName}".fna.gz ] || curl -s "$fnaDownload" \
                        > "${output}"/fna/"${fileName}".fna.gz  # Check if already downloaded. if not, download.
                    ;;
            ffn)    ffn="/"${fileName}"_cds_from_genomic.fna.gz"
                    ffnDownload=$(echo "$1" | awk -v var="$ffn" '{print $0var}')
                    [ -s "${output}"/"${fileName}".ffn.gz ] || curl -s "$ffnDownload" \
                        > "${output}"/ffn/"${fileName}".ffn.gz
                    ;;
            gbk)    gbk="/"${fileName}"_genomic.gbff.gz"
                    gbkDownload=$(echo "$1" | awk -v var="$gbk" '{print $0var}')
                    [ -s "${output}"/"${fileName}".gbk.gz ] || curl -s "$gbkDownload" \
                        > "${output}"/gbk/"${fileName}".gbk.gz
                    ;;
            gff)    gff="/"${fileName}"_genomic.gff.gz"
                    gffDownload=$(echo "$1" | awk -v var="$gff" '{print $0var}')
                    [ -s "${output}"/"${fileName}".gff.gz ] || curl -s "$gffDownload" \
                        > "${output}"/gff/"${fileName}".gff.gz
                    ;;
            faa)    faa="/"${fileName}"_protein.faa.gz"
                    faaDownload=$(echo "$1" | awk -v var="$faa" '{print $0var}')
                    [ -s "${output}"/"${fileName}".faa.gz ] || curl -s "$faaDownload" \
                        > "${output}"/faa/"${fileName}".faa.gz
                    ;;
        esac
    done
}

# Make function available to parallel
export -f download  # -f is to export functions

# Download multiple assemblies in parallel
# for i in $(cat "${output}"/completeGenomePaths.txt); do
#     download "$i"
# done
# exit

cat "${output}"/completeGenomePaths.txt \
    | parallel  --bar \
                --delay 0.3 \
                --env download \
                --env output \
                --env seq_type \
                --jobs "$cpu" \
                'download {}'

# Remove empty files
find "$output" -size 0 -type f ! -name "*.txt" -exec rm {} \;

# Rename the file with the organism name from the fasta header
if [ "$rename" -eq 1 ] && [ $(ls "${output}"/fna | wc -l) -gt 1 ]; then
    echo "Renaming files..."

    function rename()
    {
        gfc=$(basename "$1" ".fna")
        name=$(echo "$gfc" | cut -d "_" -f 1,2 --output-delimiter="_")
        path=$(dirname "$1")
        
        organism_name=$(cat "$filtered_summary" | grep "$name" | awk -F $'\t' 'BEGIN {OFS = FS} {print $8, $9}')
        id=$(echo "$organism_name" | cut -f 1 | sed 's/str.*//' | tr " " "_" | tr -d ".")
        strain=$(echo "$organism_name" | cut -f 2 | cut -d "=" -f 2 | tr " " "_")
        new_name=$(echo ""${name}"_"${id}"_strain_"${strain}"" | sed -e 's/__/_/g' -e 's%/%-%g')

        mv "$1" "${path}"/"${new_name}".fna.gz
        [ -s "${output}"/ffn/"${name}".ffn.gz ] && mv "${output}"/ffn/"${name}".ffn.gz "${output}"/ffn/"${new_name}".ffn.gz
        [ -s "${output}"/gff/"${name}".gff.gz ] && mv "${output}"/gff/"${name}".gff.gz "${output}"/gff/"${new_name}".gff.gz
        [ -s "${output}"/faa/"${name}".faa.gz ] && mv "${output}"/faa/"${name}".faa.gz "${output}"/faa/"${new_name}".faa.gz
        [ -s "${output}"/gbk/"${name}".gbk.gz ] && mv "${output}"/gbk/"${name}".gbk.gz "${output}"/gbk/"${new_name}".gbk.gz
    }

    export -f rename

    find "$output" -type f -name "*.fna.gz" \
        | parallel  --bar \
                    --env rename \
                    --env filtered_summary \
                    --env output \
                    --jobs "$cpu" \
                    'rename {}'

fi

# TODO -> create folder for each Genus or Species?

# Uncompress data
if [ "$uncompress" -eq 1 ]; then
    echo "Uncompressing files..."
    find "$output" -type f -name "*.gz" \
        | parallel  --bar \
                    --jobs $((cpu/3)) \
                    'pigz -p 3 -d {}'
fi
