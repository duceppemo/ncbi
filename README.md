# ncbi
Compilation of scrips related to stuff to do on NCBI

#TODO
Create proper documentation for the scripts (including dependencies etc.)


## http_postFM.pl
Dependencies (perl modules):
- Sys::CPU
- Parallel::ForkManager
```
conda create -n ncbi -c bioconda perl-app-cpanminus perl-sys-cpu perl-parallel-forkmanager
conda activate ncbi
```
usage: perl http_postFM.pl -h

## ncbi_wgs_downloader.py
A script to download assemblies from shotgut sequencing projects on GenBank with accession numbers like `AABUZE000000000.2`.

usage: python3 ncbi_wgs_downloader.py <accession_list.txt> <download_folder/> <you_email@gmail.com>

where `accession_list.txt` is a text file with one accession number per line, `download_folder/` is the directory where you want your fasta files to be downlaoded and `your_email@gmail.com` is the email address you want to use for the handshake with ncbi server.
