#!/bin bash

today=$(date +%Y-%m-%d)

dwl="/media/6tb_raid10/db/centrifuge/download"
dwl_plant=""${dwl}"/plant"

[ -d "$dwl" ] || mkdir -p "$dwl"
[ -d "$dwl_plant" ] || mkdir -p "$dwl_plant"


#############
#           #
#   Fasta   #
#           #
#############


#make accession list
perl ~/scripts/ncbiAccessionsFromTaxidFM.pl \
    -db "nucleotide" \
    --query 'txid33090[Organism:exp] AND "internal transcribed spacer"' \
    --output "${dwl}"/plant_ITS_accession.list

#Download sequences from accession list
perl ~/scripts/http_postFM.pl \
    -i "${dwl}"/plant_ITS_accession.list \
    -o "${dwl}"/plant_ITS_accession.fasta \
    -t nucleotide \
    -n $(nproc)

#remove version part of the accession number
cat "${dwl}"/plant_ITS_accession.fasta \
    | sed 's/\..//1' \
    > "${dwl}"/plant_ITS.fata  

#cleanup
rm "${dwl}"/plant_ITS_accession*


#############
#           #
#   TaxID   #
#           #
#############


#Accession to taxID
#Genbank files to download
echo "ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz
ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/nucl_wgs.accession2taxid.gz" \
> "${dwl}"/URLs.list

#Parallel download all the genbank files
cat "${dwl}"/URLs.list | xargs -n 1 -P 10 wget -q -P "$dwl"

#Paralle decompress all the downloaded files
find "$dwl" -maxdepth 1 -type f -name "*.gz" | parallel 'pigz -d {}'

#make files compatible with centrifuge-build
find "$dwl" -maxdepth 1 -type f -name "*accession2taxid" \
    | parallel "cat {} | sed '1d' | cut -f 1,3 > {}.map"

rm "${dwl}"/*.accession2taxid

cat "${dwl}"/*.map > "${dwl}"/acc2taxid.txt
rm "${dwl}"/*.map


################
#              #
#   Taxonomy   #
#              #
################


#update taxonomy
[ -d "${dwl}"/taxonomy ] || mkdir -p "${dwl}"/taxonomy
wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz -P "${dwl}"/taxonomy
tar zxvf "${dwl}"/taxonomy/taxdump.tar.gz -C "${dwl}"/taxonomy


##################
#                #
#   Centrifuge   #
#                #
##################


#Build centrifuge compressed database
time centrifuge-build \
    -p $(nproc) \
    --conversion-table "${dwl}"/acc2taxid.txt \
    --taxonomy-tree "${dwl}"/taxonomy/nodes.dmp \
    --name-table "${dwl}"/taxonomy/names.dmp \
    "${dwl}"/plant_ITS.fasta \
    "${dwl}"/"${today}"_plant_ITS

# nt            1371m21.118s -> 22.85h
# vir           3m20s
# plant_ITS     3m4s
# bact_vir_h    4.27h

#move databases
mv "${dwl}"/"${today}"_plant_ITS* "${dwl}"/..

#delete temp files
rm -rf "$dwl"


###################
#                 #
#   Krona Tools   #
#                 #
###################


#Gotta update krona tools too!
bash ~/prog/KronaTools-2.7/updateTaxonomy.sh
bash ~/prog/KronaTools-2.7/updateAccessions.sh

