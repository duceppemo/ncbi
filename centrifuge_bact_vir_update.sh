#!/bin bash

today=$(date +%Y-%m-%d)

dwl="/media/6tb_raid10/db/centrifuge/download"
dwl_bact=""${dwl}"/bacteria"
dwl_vir=""${dwl}"/viruses"
dwl_human=""${dwl}"/human"

[ -d "$dwl" ] || mkdir -p "$dwl"
[ -d "$dwl_bact" ] || mkdir -p "$dwl_bact"
[ -d "$dwl_vir" ] || mkdir -p "$dwl_vir"
[ -d "$dwl_human" ] || mkdir -p "$dwl_human"

# Dowload all sequences from RefSeq
bash ~/scripts/refseqDownloader.sh \
    -t "bacteria" \
    -l "Complete Genome" \
    -o "$dwl_bact" \
    -u

bash ~/scripts/refseqDownloader.sh \
    -t "viral" \
    -l "Complete Genome" \
    -o "$dwl_vir" \
    -u

bash ~/scripts/refseqDownloader.sh \
    -t "human" \
    -o "$dwl_human" \
    -u

#Concatenate all downloaded files by taxon
declare -a mylist=("${dwl_bact}" "${dwl_vir}" "${dwl_human}")
for i in "${mylist[@]}"; do
    taxon=$(basename "$i")
    for j in $(find "${i}"/fna -type f -name "*.fna"); do
        cat "$j" | sed 's/\..//1' >> "${i}"/"${taxon}"_all.fasta  #remove version part of the accession number
        rm "$j"  # delete file downloaded to avoid double footprint
    done
done

#Accession to taxID
#Genbank files to download
echo "ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz
ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/nucl_wgs.accession2taxid.gz" \
> "${dwl}"/URLs.list

#Parallel download all the genbank files
cat "${dwl}"/URLs.list | xargs -n 1 -P 10 wget -q -P "$dwl"

#Paralle decompress all the downloaded files
find "$dwl" -maxlevel 1 -type f -name "*.gz" | parallel 'pigz -d {}'

#make files compatible with centrifuge-build
cat "${dwl}"/nucl_gb.accession2taxid \
    | sed '1d' \
    | cut -f 1,3 \
    > "${dwl}"/nucl_gb.accession2taxid.map

cat "${dwl}"/nucl_wgs.accession2taxid \
    | sed '1d' \
    | cut -f 1,3 \
    > "${dwl}"/nucl_wgs.accession2taxid.map

rm "${dwl}"/*.accession2taxid

cat "${dwl}"/*.map > "${dwl}"/acc2taxid.txt
rm "${dwl}"/*.map


#update taxonomy
[ -d "${dwl}"/taxonomy ] || mkdir -p "${dwl}"/taxonomy
wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz -P "${dwl}"/taxonomy
tar zxvf "${dwl}"/taxonomy/taxdump.tar.gz

#Build centrifuge compressed database
time centrifuge-build \
    -p $(nproc) \
    --conversion-table "${dwl}"/acc2taxid.txt \
    --taxonomy-tree "${dwl}"/taxonomy/nodes.dmp \
    --name-table "${dwl}"/taxonomy/names.dmp \
    "${dwl_bact}"/bacteria_all.fasta,"${dwl_vir}"/viruses_all.fasta,"${dwl_human}"/human_all.fasta \
    "${dwl}"/"${today}"_bact_vir_h

# nt            1371m21.118s -> 22.85h
# vir           3m20s
# bact_vir_h    4.27h

#move databases
mv "${dwl}"/"${today}"_bact_vir_h* "${dwl}"/..

#delete temp files
rm -rf "$dwl"

#Gotta update krona tools too!
bash ~/prog/KronaTools-2.7/updateTaxonomy.sh
bash ~/prog/KronaTools-2.7/updateAccessions.sh

