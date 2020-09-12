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
        handle = Entrez.efetch(db=db, id=accession, rettype='gb', retmode='text')
    except urllib.error.HTTPError:
        print('ERROR')
        error_acc_handle.write(accession + '\n')
        return

    if handle:
        try:
            record = SeqIO.read(handle, "genbank")
            handle.close()
        # except urllib.error.HTTPError:
        except http.client.IncompleteRead:
            print('ERROR')
            error_acc_handle.write(accession + '\n')
    else:
        print('ERROR')
        error_acc_handle.write(accession + '\n')
        return  # skip accession

    # Find out contigs to get
    contigs = record.annotations['wgs']
    # version = accession.split('.')[1]
    # version = str(version).zfill(2)
    # Extract alphabetic characters at the begining of the accession
    prefix = "".join(re.split("[^a-zA-Z]*", accession))
    version = "".join(re.split("[^0-9]*", contigs[0]))[:2]
    num = "".join(re.split("[^0-9]*", contigs[0]))[2:]
    num_len = len(num)
    start = int(contigs[0][-num_len:])
    stop = int(contigs[1][-num_len:])
    # https://www.ncbi.nlm.nih.gov/nuccore/AABKDH010000003.1?report=fasta

    # List of all the files to get at once
    acc_list = list()
    for c in range(start, stop + 1, 1):
        acc = prefix + version + str(c).zfill(num_len)
        acc_list.append(acc)

    # Write sequences to file
    with open(output_fasta, 'a') as outfasta:
        try:
            handle = Entrez.efetch(db=db, id=acc_list, rettype='fasta', retmode='text')
        except urllib.error.HTTPError:
            print('ERROR')
            error_acc_handle.write(accession + '\n')
            return

        if handle:
            try:
                rec_dict = SeqIO.to_dict(SeqIO.parse(handle, 'fasta'))
                SeqIO.write(rec_dict.values(), outfasta, 'fasta')
                sleep(1)
                print('DONE')
            except http.client.IncompleteRead:
                print('ERROR')
                error_acc_handle.write(accession + '\n')
        else:
            error_acc_handle.write(accession + '\n')
            print('ERROR')


# Loop the accession list
for i in accession_list:
    fetch_wgs_contigs(i)

error_acc_handle.close()
