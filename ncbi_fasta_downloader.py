#!/usr/local/env python3
import http.client
from argparse import ArgumentParser
from multiprocessing import cpu_count
from time import time
import shutil
from urllib.request import urlopen
from pathlib import Path
import os
from Bio import Entrez
from concurrent import futures
from time import sleep
import numpy as np
from itertools import groupby


__author__ = 'duceppemo'
__version__ = '0.1'


# TODO: add download progressbar
#  (https://stackoverflow.com/questions/41106599/python-3-5-urllib-request-urlopen-progress-bar-available)


class Methods(object):

    @staticmethod
    def elapsed_time(seconds):
        """
        Transform a time value into a string

        :param seconds: Elapsed time in seconds
        :return: Formatted time string
        """

        minutes, seconds = divmod(seconds, 60)
        hours, minutes = divmod(minutes, 60)
        days, hours = divmod(hours, 24)
        periods = [('d', days), ('h', hours), ('m', minutes), ('s', seconds)]
        time_string = ''.join('{}{}'.format(int(np.round(value)), name) for name, value in periods if value)
        return time_string

    @staticmethod
    def make_folder(folder):
        """
        Create output folder.

        :param folder: string. Output folder path.
        """

        # Will create parent directories if they don't exist
        # Will not return error if it already exists
        Path(folder).mkdir(parents=True, exist_ok=True)

    @staticmethod
    def download(url, file_path):
        """
        Download the file from `url` and save it locally under `file_name`
        https://stackoverflow.com/questions/7243750/download-file-from-web-in-python-3
        https://unite.ut.ee/repository.php

        :param url: string. URL of file (remotely)
        :param file_path: string. Path of downloaded file (locally)
        """

        with open(file_path, 'wb') as out_file:
            with urlopen(url) as response:
                shutil.copyfileobj(response, out_file)

    @staticmethod
    def download_parallel(url_list, path_list, cpu):
        """
        Run the download function in parallel.

        :param url_list: list. A list of URLs
        :param path_list: list. A list of paths
        :param cpu: int or string. Number of parallel downloads
        """

        with futures.ThreadPoolExecutor(max_workers=cpu) as executor:
            args = ((u, p) for u in url_list for p in path_list)
            for results in executor.map(lambda p: Methods.download(*p), args):  # (*p) unpacks arguments
                pass

    @staticmethod
    def download_seq_from_query(query, seq_file, email, api, single):
        """
        # http://biopython.org/DIST/docs/tutorial/Tutorial.html#sec167

        :param query: string. NCBI query
        :param seq_file: string. Path to output file
        :param email: string. Email address
        :param api: string. NCBI API key
        """

        Entrez.email = email
        Entrez.api_key = api
        Entrez.sleep_between_tries = 15
        Entrez.max_tries = 3

        if os.path.exists(query) and os.path.isfile(query):  # Query is a file with one acc per line
            acc_list = list()
            with open(query, 'r') as f:
                for line in f:
                    line = line.rstrip()
                    if not line:
                        continue
                    acc_list.append(line)
            count = len(acc_list)

            print('Downloading {} records...'.format(count))

            batch_size = 300
            with open(seq_file, "w") as out_handle:
                for start in range(0, count, batch_size):
                    end = min(count, start + batch_size)
                    print("\tDownloading record %i to %i (%i)" % (start + 1, end, count))
                    try:  # Trying to deal with "http.client.IncompleteRead
                        fetch_handle = Entrez.efetch(
                            db="nucleotide",
                            rettype="fasta",
                            retmode="text",
                            id=','.join(acc_list[start:end]))
                    except (http.client.IncompleteRead, ValueError) as e:
                        print('Network error ({}). Last attempt.'.format(e))
                        fetch_handle = Entrez.efetch(
                            db="nucleotide",
                            rettype="fasta",
                            retmode="text",
                            id=','.join(acc_list[start:end]))

                    data = fetch_handle.read()
                    fetch_handle.close()
                    out_handle.write(data)
                    sleep(0.5)

        else:  # query is a string
            # Search
            print('Searching NCBI for \'{}\'.'.format(query))
            # ESearch URL can only retrieve up to 100,000 UIDs
            search_handle = Entrez.esearch(db='nucleotide', term=query, idtype="acc", usehistory='y')
            search_results = Entrez.read(search_handle)
            search_handle.close()

            count = int(search_results["Count"])
            webenv = search_results["WebEnv"]
            query_key = search_results["QueryKey"]

            # Fetch
            print('Downloading {} records...'.format(count))
            batch_size = 300
            with open(seq_file, "w") as out_handle:
                for start in range(0, count, batch_size):
                    end = min(count, start + batch_size)
                    print("\tDownloading record %i to %i (%i)" % (start + 1, end, count))
                    try:  # Trying to deal with "http.client.IncompleteRead
                        fetch_handle = Entrez.efetch(
                            db="nucleotide",
                            rettype="fasta",
                            retmode="text",
                            retstart=start,
                            retmax=batch_size,
                            webenv=webenv,
                            query_key=query_key,
                            idtype="acc")
                    except (http.client.IncompleteRead, ValueError) as e:
                        print('Network error ({}). Last attempt.'.format(e))
                        fetch_handle = Entrez.efetch(
                            db="nucleotide",
                            rettype="fasta",
                            retmode="text",
                            retstart=start,
                            retmax=batch_size,
                            webenv=webenv,
                            query_key=query_key,
                            idtype="acc")

                    data = fetch_handle.read()
                    fetch_handle.close()
                    out_handle.write(data)
                    sleep(0.5)

    @staticmethod
    def fasta_iterator(input_fasta):
        """
        Creates an iterator for multi fasta files.

        :param input_fasta: string. Path to fasta file.
        :return: iterator.
        """

        with open(input_fasta, 'rb') as in_fh:
            fasta_iter = (x[1] for x in groupby(in_fh, lambda line: str(line, 'utf-8')[0] == ">"))
            for entry in fasta_iter:
                header = str(entry.__next__(), 'utf-8')
                acc = header.split()[0][1:]
                seq = "".join(str(s, 'utf-8').strip() for s in fasta_iter.__next__())
                yield acc, header, seq

    @staticmethod
    def multi_to_single_fasta(input_fasta, output_folder):
        """
        https://www.biostars.org/p/710/
        Split a multi fasta file into multiple single fasta files, one file per entry.

        :param input_fasta: string. Input multi fasta file
        :param output_folder: string. Output folder
        """

        print('Splitting downloaded sequences into individual fasta files...')
        fasta_iter = Methods.fasta_iterator(input_fasta)
        for i in fasta_iter:
            acc, header, seq = i
            output_file = output_folder + '/' + acc + '.fasta'
            with open(output_file, 'w') as f:
                f.write('{}{}\n'.format(header, seq))

        # Remove multi fasta file
        os.remove(input_fasta)


class NcbiDownloader(object):
    def __init__(self, args):
        self.args = args

        # Args
        self.query = args.query
        self.output_folder = args.output
        self.cpu = args.threads
        self.email = args.email
        self.api = args.api_key
        self.single = args.single

        # Run
        self.run()

    def run(self):
        """
        Run all the steps
        """

        # Time tracking
        t_zero = time()

        # Check a few things before starting processing data
        self.checks()

        # Download nucleotide sequences of the query result
        seq_file = self.output_folder + '/seq.fasta'
        Methods.download_seq_from_query(self.query, seq_file, self.email, self.api, self.single)

        # Split into individual files
        if self.single:
            Methods.multi_to_single_fasta(seq_file, self.output_folder)

        end_time = time()
        interval = end_time - t_zero
        print("Done (total time: %s)." % Methods.elapsed_time(interval))

    def checks(self):
        """
        Perform a bunch of check to make sure we have everything we need before we start.
        :return:
        """

        # Check query string
        if not self.query:
            raise Exception('Your query is empty. Please supply a text file with accession numbers (one per line) or '
                            'a query string.')

        # Check output folder
        if not os.path.exists(self.output_folder):
            Methods.make_folder(self.output_folder)  # Create output folder if it does not exist

        # Check cpu and parallel processes
        cpu = cpu_count()
        if 1 > self.cpu > cpu:  # smaller than 1 or greater than available cpu
            self.cpu = cpu


if __name__ == '__main__':

    parser = ArgumentParser(description='Download DNA sequences from NCBI from a query string or an accession list.')
    parser.add_argument('-q', '--query', metavar='\"txid4762[Organism:exp] AND (\"internal transcribed spacer\"[Title]) NOT uncultured[Title]\"',
                        required=True,
                        type=str,
                        help='NCBI query string OR a text file with one accession number per line. Mandatory.')
    parser.add_argument('-o', '--output', metavar='/output_folder/',
                        required=True,
                        type=str,
                        help='Output folder. Mandatory.')
    parser.add_argument('-s', '--single',
                        required=False,
                        action='store_true',
                        help='Output retrieved sequences in individual fasta files instead of a single file.')
    parser.add_argument('-t', '--threads', metavar='4',
                        required=False, default=4,
                        type=int,
                        help='Number of CPU. Default is 4. Optional.')
    parser.add_argument('-e', '--email', metavar='your.email@example.org',
                        required=False, default='\'your.email@example.org\'',
                        type=str,
                        help='Your email address. Optional.')
    parser.add_argument('-a', '--api-key', metavar='',
                        required=False,
                        type=str,
                        help='Your NCBI API key, if you have one. Allows up to 10 requests per second instead of 3. '
                             'Optional.')
    parser.add_argument('-v', '--version', action='version',
                        version=f'{os.path.basename(__file__)}: v{__version__}')

    # Get the arguments into an object
    arguments = parser.parse_args()

    NcbiDownloader(arguments)
