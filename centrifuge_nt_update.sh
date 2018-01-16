wget ftp://ftp.ncbi.nih.gov/blast/db/FASTA/nt.gz
pigz -d nt.gz
mv -v nt nt.fa
wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz
pigz -d nucl_gb.accession2taxid.gz
cat nucl_gb.accession2taxid | sed '1d' | cut -f 2,3 > nucl_gb.accession2taxid.map


time centrifuge-build \
    -p $(nproc) \
    --bmax 1342177280 \
    --conversion-table nucl_gb.accession2taxid.map \
    --taxonomy-tree ../taxonomy/nodes.dmp \
    --name-table ../taxonomy/names.dmp \
    nt.fa \
    ../nt

# real    1371m21.118s -> 22.85h


#Gotta update krona tools too!
bash ~/prog/KronaTools-2.7/updateTaxonomy.sh
bash ~/prog/KronaTools-2.7/updateAccessions.sh

