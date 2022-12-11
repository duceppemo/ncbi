#!/bin/bash

version="0.2.2"

# Script description:

# This script will download ".fna" ".ffn", ".faa" and ".gff", if availables, from NCBI RefSeq for a specific database.
# Assembly level and name filtering are also available as options.


#################
#               #
#    Options    #
#               #
#################


function displayHelp()
{
    echo -e "\

    Usage: refseqDownloader.sh -d db_type -o output_folder [-q \"organism\"] [-l \"assemmly_level\"] -r -u -n 4

    Mandatory flags:

        -d          RefSeq database type (bacteria, archaea, viral, fungi, protozoa, vertebrate, plant, human).

        -o          Output folder. Where downloaded files will be located.

    Optional flags:

        -h          Print this help message.

        -q          Query. Organism name(s) to keep.
                    Comma separated with no space if multiple (e.g. "Listeria monocytogenes,Campylobacter,Salmonella").
                    Not using this option will download all species of the selected database type.

        -l          Assembly level (Chromosome, Complete Genome, Contig, Scaffold).
                    Comma separated with no space if multiple (e.g. "Complete Genome,Scaffold").
                    Not using this option will download all assembly levels of the selected database type.

        -r          Rename files according to description in fna file.
                    "fna" file type must be selected to enable renaming.

        -t          Download only type (representative) strains.

        -u          Uncompress all downloaded files.

        -f          File type(s) to download(fna, ffn, gff, faa).
                    Comma separated with no space if multiple (e.g. "fna,ffn,gff,faa,gbk").
                    Not using this option will download all three file types.

        -n          Number of CPUs to use for parallel download.
                    Default is max numbers of CPUs avaiable.
    "
}


#Colored error message
BLUE='\033[1;34m'
NC='\033[0m' # No Color

# Default values
db_type=''
export output=''
name=''
level=''
type=0
uncompress=0
rename=0
export seq_type=''  # 'fna,ffn,gff,faa,gbk'
export cpu=$(nproc)
help=''

# If the very first character of the option-string is a : (colon),
# which would normally be nonsense because there's no option letter preceding it,
# getopts switches to "silent error reporting mode"
# When you want getopts to expect an argument for an option, just place a : (colon) after the proper option flag
options=':d:o:n:l:utf:rq:h'

while getopts "$options" opt; do
    case "$opt" in
        d)  db_type="$OPTARG";;
        o)  export output="$OPTARG";;
        q)  name="$OPTARG";;
        l)  level="$OPTARG";;
        h)  displayHelp
            exit 1
            ;;
        r)  rename=1;;
        t)  type=1;;
        u)  uncompress=1;;
        f)  seq_type="$OPTARG";;
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
if [[ -z "$db_type" ]] || [[ -z "$output" ]]; then
    echo "Both \"-t\" and \"-o\" options are mandatory and require arguments"
    displayHelp
    exit 1
fi

#test file type to get
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
    [ -e "${output}"/"${db_type}"_"${downloadDate}".txt ] && rm "${output}"/"${db_type}"_"${downloadDate}".txt
    touch "${output}"/"${db_type}"_"${downloadDate}".txt

    #Download the refseq assembly summary file
    cd "$output"
    wget --timestamp "$summary_url" -q --show-progress #-P "$output"
    summary_file=""${output}"/assembly_summary.txt"
fi


#################
#               #
#   Filtering   #
#               #
#################


#make sure "Complete Genome" is treated as a single array element
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
if [[ -n "$name" ]]; then
    name_filter=$(echo "'"$name"'" | tr "," "|")
fi

echo "$name_filter"

# create a subset of assembly_summary using the name(s) supplied
export filtered_summary=""${output}"/assembly_summary_filtered.txt"
cmd1="cat "$summary_file" | grep -E "$name_filter" > "$filtered_summary""
eval "$cmd1"

# create another subset if representive only was selected
cmd2="cat "$filtered_summary" | grep -E 'representative genome' > "${filtered_summary}".tmp"
eval "$cmd2"
mv "${filtered_summary}".tmp "$filtered_summary"

# File with download paths
if [ -n "$level" ] && [ -n "$name" ]; then
    echo "Downloading all "$level" sequences for "$name"..."
    cmd="cat "$filtered_summary" | awk -F \"\\t\" '"$f" && \$11==\"latest\" {print \$20}' > "${output}"/completeGenomePaths.txt"
elif [ -n "$level" ] && [ -z "$name" ]; then
    echo "Downloading all "$level" sequences for "$db_type"..."
    cmd="cat "$summary_file" | awk -F \"\\t\" '"$f" && \$11==\"latest\" {print \$20}' > "${output}"/completeGenomePaths.txt"
elif [ -z "$level" ] && [ -n "$name" ]; then
    echo "Downloading all available sequences for "$name"..."
    cmd="cat "$filtered_summary" | awk -F \"\\t\" '\$11==\"latest\" {print \$20}' > "${output}"/completeGenomePaths.txt"
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

#file types
seq_types=($(echo "$seq_type" | tr "," " "))  # convert string to array

#check if correct sequence type terms are being used
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

function check_target()
{
    if [[ $(wget -S --spider $1  2>&1 | grep 'HTTP/1.1 200 OK') ]]; then 
        echo "true"
    fi
}

# Fucntion to parallel download
function download()
{
    s_types=($(echo "$seq_type" | tr "," " "))
    fileName=$(basename "$1")

    #Determine which files to download
    for t in "${s_types[@]}"; do
        case "$t" in  # indirect variable expansion
            fna)    fna="/"${fileName}"_genomic.fna.gz"
                    fnaDownload=$(echo "$1" | awk -v var="$fna" '{print $0var}')  #Create the full file path
                    [ -s "${output}"/"${fileName}".fna.gz ] || curl -s "$fnaDownload" \
                        > "${output}"/fna/"${fileName}".fna.gz #check if already downloaded. if not, download.
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

#make function available to parallel
export -f download  # -f is to export functions

# Download multiple refseq genomes in parallel
cat "${output}"/completeGenomePaths.txt \
    | parallel  --bar \
                --delay 0.3 \
                --env download \
                --env output \
                --env seq_type \
                --jobs "$cpu" \
                'download {}'

#remove empty files
find "$output" -size 0 -type f ! -name "*.txt" -exec rm {} \;

#rename the file with the organism name from the fasta header
if [ "$rename" -eq 1 ] && [ $(ls "${output}"/fna | wc -l) -gt 1 ]; then
    echo "Renaming files..."

    function rename()
    {
            gfc=$(basename "${1%.fna}")
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

#uncompress?
if [ "$uncompress" -eq 1 ]; then
    echo "Uncompressing files..."
    find "$output" -type f -name "*.gz" \
        | parallel  --bar \
                    --jobs $((cpu/3)) \
                    'pigz -p 3 -d {}'
fi
