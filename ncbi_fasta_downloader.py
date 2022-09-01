#!/usr/local/env python3

from sys import argv
from Bio import Entrez
from Bio import SeqIO
import re
import http
import urllib
from time import sleep

input_list = argv[1]
output_folder = argv[2]
email = argv[3]

error_acc = output_folder + '/acc_errors.txt'
error_acc_handle = open(error_acc, 'w')

# Parse accession list
accession_list = list()
with open(input_list, 'r') as f:
    for line in f:
        line = line.rstrip()
        if not line:
            continue
        accession_list.append(line)


def fetch_wgs_contigs(accession):
    """
    1- Use NCBI Entrez to get the contigs accession numbers of a sequencing project.
    2- Download all contigs into a single file for each project.
    :param accession: string; genbank accession number from a shotgun sequencing project (e.g. NZ_AABUZE000000000.2)
    :return:
    """
    print('Processing {}... '.format(accession), end='', flush=True)
    output_fasta = output_folder + '/' + accession + '.fasta'

    # Initialize some default parameters
    Entrez.email = email  # Provide your email address
    db = 'nucleotide'  # Set search to dbVar database

    # Get WGS project info in genbank format
    handle = ''
    try:
        handle = Entrez.efetch(db=db, id=accession, rettype='fasta', retmode='text')
    except urllib.error.HTTPError:
        print('ERROR')
        error_acc_handle.write(accession + '\n')
        return

    if handle:
        try:
            with open(output_fasta, 'w') as f:
                f.write(handle.read())
                # record = SeqIO.write(handle, f, "fasta")
                handle.close()
        # except urllib.error.HTTPError:
        except http.client.IncompleteRead:
            print('ERROR')
            error_acc_handle.write(accession + '\n')
    else:
        print('ERROR')
        error_acc_handle.write(accession + '\n')
        return  # skip accession

# Loop the accession list
for i in accession_list:
    fetch_wgs_contigs(i)

error_acc_handle.close()
