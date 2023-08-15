# ncbi
Compilation of scrips related to stuff to do on NCBI

## Installation
```shell
# Create virtual environment
# `conda` can be used instead of `mamba`
mamba create -n ncbi -y biopython perl-lwp-simple perl-sys-cpu perl-parallel-forkmanager
```

## ncbi_fasta_downloader.py

```text
usage: python ncbi_fasta_downloader.py [-h] -q "txid4762[Organism:exp] AND "internal
                                       transcribed spacer"[Title] NOT
                                       uncultured[Title]" -o /output_folder/ [-s]
                                       [-t 4] [-e your.email@example.org] [-a] [-v]

Download DNA sequences from NCBI from a query string or an accession list.

optional arguments:
  -h, --help            show this help message and exit
  -q "txid4762[Organism:exp] AND ("internal transcribed spacer"[Title]) NOT uncultured[Title]", --query "txid4762[Organism:exp] AND ("internal transcribed spacer"[Title]) NOT uncultured[Title]"
                        NCBI query string OR a text file with one accession
                        number per line. Mandatory.
  -o /output_folder/, --output /output_folder/
                        Output folder. Mandatory.
  -s, --single          Output retrieved sequences in individual fasta files
                        instead of a single file.
  -t 4, --threads 4     Number of CPU. Default is 4. Optional.
  -e your.email@example.org, --email your.email@example.org
                        Your email address. Optional.
  -a , --api-key        Your NCBI API key, if you have one. Allows up to 10
                        requests per second instead of 3. Optional.
  -v, --version         show program's version number and exit
```

## ncbi_wgs_downloader.py
A script to download assemblies from shotgun sequencing projects on GenBank with accession numbers like `AABUZE000000000.2`.
```shell
# Usage
python ncbi_wgs_downloader.py \
  accession_list.txt \
  /download_folder \
  you_email@gmail.com
```

## http_postFM.pl
Legacy.
```shell
# Usage
perl http_postFM.pl -h
```