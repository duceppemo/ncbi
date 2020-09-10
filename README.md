# ncbi
Compilation of scrips related to stuff to do on NCBI

#TODO
Create proper documentation for the scripts (including dependencies etc.)


## http_postFM.pl
Dependencies (perl modules):
- Sys::CPU
- Parallel::ForkManager

usage: perl http_postFM.pl -h

## trace_downloader.py
A script to download assemblies from Shotgut sequencing projects on GenBank with accession numbers like `AABUZE000000000.2`.

usage: python3 trace_downloader.py <accession_list.txt> <download_folder/>

where `accession_list.txt` is a text file with one accession number per line and `download_folder/` is the directory where you want your fasta files to be downlaoded.
