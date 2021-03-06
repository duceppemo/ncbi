#!/usr/local/env python3

__author__ = 'duceppemo'
__version__ = '0.2'


import os
import sys
from nested_dict import nested_dict
from time import time


class FormatTaxo(object):

    def __init__(self, args):
        import multiprocessing

        self.args = args
        self.taxo = args.taxo
        self.accession = args.accession
        self.out = args.output

        # Number of cpu
        self.cpus = int(multiprocessing.cpu_count())

        # Main data object
        self.acc_list = list()
        self.acc_dict = nested_dict()  # test
        self.acc2taxid_dict = nested_dict()
        self.merged_dict = nested_dict()
        self.nodes_dict = nested_dict()
        self.names_dict = nested_dict()
        self.taxo_dict = nested_dict()

        # run the script
        self.run()

    def run(self):
        self.check_dependencies()
        self.parse_accession_list(self.acc_list, self.accession)
        self.parse_merged(self.merged_dict, self.taxo + '/' + 'merged.dmp')
        self.parse_nodes(self.nodes_dict, self.taxo + '/' + 'nodes.dmp')
        self.parse_names(self.names_dict, self.taxo + '/' + 'names.dmp')
        self.parse_acc2taxid(self.taxo + '/' + 'nucl_gb.accession2taxid', self.acc2taxid_dict)
        self.check_merged(self.acc2taxid_dict, self.merged_dict)
        self.make_dict(self.acc2taxid_dict, self.taxo_dict)
        self.make_lineage()
        self.make_ouput(self.taxo_dict, self.out + '/' + 'output.tsv')

    def check_dependencies(self):
        """
        Check if taxonomy files are present. If not download them
        :return:
        """

        import pathlib
        import subprocess

        # Default output path is input assembly path
        if not self.out:
            self.out = os.path.dirname(self.accession) + "/"

        # Name of required files
        taxdump_nodes = 'nodes.dmp'
        taxdump_names = 'names.dmp'
        taxdump_merged = 'merged.dmp'
        acc2taxo = 'nucl_gb.accession2taxid'
        needed = [taxdump_nodes, taxdump_names, acc2taxo]

        taxdump_url = 'https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz'
        acc2taxo_url = 'https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz'
        taxo_file_urls = [taxdump_url, acc2taxo_url]

        # Check if taxonomy folder was supplied
        if not self.taxo:
            self.taxo = self.out + 'taxonomy'

        # Check if supplied taxo folder is a folder
        if not os.path.isdir(self.taxo):
            pathlib.Path(self.taxo).mkdir(parents=True, exist_ok=True)  # If not create it

        # Check proper files are in self.taxo folder
        # Taxdump files
        if not os.path.isfile(self.taxo + '/' + taxdump_nodes)\
                or not os.path.isfile(self.taxo + '/' + taxdump_names)\
                or not os.path.isfile(self.taxo + '/' + taxdump_merged):
            if not os.path.isfile(self.taxo + '/' + 'taxdump.tar.gz'):
                self.download(taxdump_url, self.taxo + '/' + 'taxdump.tar.gz')
            else:
                # TODO: check if the file is latest and complete if already present
                pass
            os.chdir(self.taxo)
            subprocess.run(['tar', '-zxvf', self.taxo + '/' + 'taxdump.tar.gz',
                            'nodes.dmp', 'names.dmp', 'merged.dmp'])
            # os.remove(self.taxo + '/' + 'taxdump.tar.gz')

        # Accession to taxid file(s)
        if not os.path.isfile(self.taxo + '/' + taxdump_nodes):
            if not os.path.isfile(self.taxo + '/' + 'nucl_gb.accession2taxid'):
                self.download(acc2taxo_url, self.taxo + '/' + 'nucl_gb.accession2taxid.gz')
            else:
                # TODO: check if the file is latest and complete if already present
                pass
            subprocess.run(['pigz', '-d', self.taxo + '/' + 'nucl_gb.accession2taxid.gz'])

    def hbytes(self, num):
        """
        Convert bytes to KB, MB, GB or TB
        :param num: a file size in bytes
        :return: A string representing the size of the file with proper units
        """
        for x in ['bytes', 'KB', 'MB', 'GB']:
            if num < 1024.0:
                return "%3.1f%s" % (num, x)
            num /= 1024.0
        return "%3.1f%s" % (num, 'TB')

    def download(self, url, filename):
        """
        https://sumit-ghosh.com/articles/python-download-progress-bar/
        :param url: Link to download file
        :param filename: file to save on disk
        :return:
        """
        import requests

        with open(filename, 'wb') as f:
            response = requests.get(url, stream=True)
            total = response.headers.get('content-length')

            if total is None:
                f.write(response.content)
            else:
                downloaded = 0
                total = int(total)
                for data in response.iter_content(chunk_size=max(int(total / 1000), 1024 * 1024)):
                    downloaded += len(data)
                    f.write(data)
                    done = int(50 * downloaded / total)
                    sys.stdout.write('\r[{}{}] ({}/{})'.format(
                        '█' * done, '.' * (50 - done),
                        self.hbytes(downloaded), self.hbytes(total)))
                    sys.stdout.flush()
        sys.stdout.write('\n')

    def parse_accession_list(self, l, f):
        """
        Parse accession files into a list
        :param l: An empty list
        :param f: The file with the accession numbers
        :return: A populated list
        """

        # counter = 0
        # total_lines = 0

        # with open(f, 'rb') as file:
        #     while True:
        #         buffer = file.read(8192*1024)
        #         if not buffer:
        #             break
        #         total_lines += buffer.count(b'\n')

        print("Parsing {}...".format(os.path.basename(f)), end="", flush=True)
        start_time = time()
        with open(f, 'rb', 2**16) as file:
            for line in file:
                # counter += 1
                line = line.rstrip()
                if not line:
                    continue
                self.acc_dict[line] = ''
                # sys.stdout.write('\rParsing "{}"... {}%'.format(os.path.basename(f),
                #                                                 round(counter / total_lines * 100, 1)))
                # sys.stdout.flush()
            end_time = time()
            interval = end_time - start_time
            # sys.stdout.write(' ({})\n'.format(self.elapsed_time(interval)))
            print(' ({})'.format(self.elapsed_time(interval)))

    def parse_merged(self, d, f):
        """
        Check for merged taxid
        :param d: acc2taxid_dict
        :param f: merger.dmp
        :return: An updated acc2taxid_dict
        """

        # counter = 0
        # total_lines = 0
        #
        # with open(f, 'rb') as file:
        #     while True:
        #         buffer = file.read(8192 * 1024)
        #         if not buffer:
        #             break
        #         total_lines += buffer.count(b'\n')

        print("Parsing {}...".format(os.path.basename(f)), end="", flush=True)
        start_time = time()
        with open(f, 'rb', 8192 * 1024) as file:
            # counter += 1
            for line in file:
                line = line.rstrip()
                if not line:
                    continue
                old_taxid, new_taxid = line.split(b'\t|\t')
                new_taxid = new_taxid.replace(b'\t|', b'')
                d[old_taxid] = new_taxid

            # sys.stdout.write('\rParsing "{}"... {}%'.format(os.path.basename(f),
            #                                                 round(counter / total_lines * 100, 1)))
            # sys.stdout.flush()
        end_time = time()
        interval = end_time - start_time
        # sys.stdout.write(' ({})\n'.format(self.elapsed_time(interval)))
        print(' ({})'.format(self.elapsed_time(interval)))

    def parse_nodes(self, d, f):
        """
        Parse the nodes.dmp file to get the full lineage for each taxid
        :param d:
        :param f:
        :return:
        """

        # counter = 0
        # total_lines = 0

        # with open(f, 'rb') as file:
        #     while True:
        #         buffer = file.read(8192*1024)
        #         if not buffer:
        #             break
        #         total_lines += buffer.count(b'\n')

        print("Parsing {}...".format(os.path.basename(f)), end="", flush=True)
        start_time = time()
        with open(f, 'rb', 8192*1024) as file:
            for line in file:
                # counter += 1
                line = line.rstrip()
                if not line:
                    continue
                taxid, parent, rank = line.split(b'\t|\t')[:3]
                d[taxid]['parent'] = parent
                d[taxid]['rank'] = rank

                # sys.stdout.write('\rParsing "{}"... {}%'.format(os.path.basename(f),
                #                                                 round(counter / total_lines * 100, 1)))
                # sys.stdout.flush()
        end_time = time()
        interval = end_time - start_time
        # sys.stdout.write(' ({})\n'.format(self.elapsed_time(interval)))
        print(' ({})'.format(self.elapsed_time(interval)))

    def parse_names(self, d, f):
        """
        Parse names.dmp into dictionary
        :param d:
        :param f:
        :return:
        """

        # counter = 0
        # total_lines = 0
        #
        # with open(f, 'rb') as file:
        #     while True:
        #         buffer = file.read(8192*1024)
        #         if not buffer:
        #             break
        #         total_lines += buffer.count(b'\n')

        print("Parsing {}...".format(os.path.basename(f)), end="", flush=True)
        start_time = time()
        with open(f, 'rb', 8192*1024) as file:
            for line in file:
                # counter += 1
                line = line.rstrip()
                if not line:
                    continue
                if b'scientific name' in line:
                    taxid, name = line.split(b'\t|\t')[:2]
                    d[taxid] = name

                # sys.stdout.write('\rParsing "{}"... {}%'.format(os.path.basename(f),
                #                                                 round(counter / total_lines * 100, 1)))
                # sys.stdout.flush()
        end_time = time()
        interval = end_time - start_time
        # sys.stdout.write(' ({})\n'.format(self.elapsed_time(interval)))
        print(' ({})'.format(self.elapsed_time(interval)))

    def parse_acc2taxid(self, f, d):
        """
        Parse the 'nucl_gb.accession2taxid' file to create a dictionary where each accession gets assigned a taxid.
        :param f: The 'nucl_gb.accession2taxid' file
        :param d: An empty dictionary dictionary
        :return: A populated dictionary

        Reading bid data
        https://www.blopig.com/blog/2016/08/processing-large-files-using-python/
        """

        # counter = 0
        # total_lines = 0
        #
        # with open(f, 'rb') as file:
        #     while True:
        #         buffer = file.read(8192*1024)
        #         if not buffer:
        #             break
        #         total_lines += buffer.count(b'\n')

        print("Parsing {}...".format(os.path.basename(f)), end="", flush=True)
        start_time = time()
        with open(f, 'rb', 8192*1024) as file:
            # accession	accession.version	taxid	gi
            # A00002	A00002.1	9913	2
            # counter += 1
            next(file)  # skip first line
            for line in file:
                # counter += 1
                line = line.rstrip()
                if not line:
                    continue
                acc, taxid = line.split(b"\t")[1:3]
                d[acc] = taxid

                # sys.stdout.write('\rParsing "{}"... {}%'.format(os.path.basename(f),
                #                                                 round(counter/total_lines*100, 1)))
                # sys.stdout.flush()
        end_time = time()
        interval = end_time - start_time
        # sys.stdout.write(' ({})\n'.format(self.elapsed_time(interval)))
        print(' ({})'.format(self.elapsed_time(interval)))

    def check_merged(self, a_dict, m_dict):
        print("Checking for merged taxids...", end="", flush=True)
        start_time = time()
        for acc, taxid in a_dict.items():
            # counter += 1
            if taxid in m_dict:
                new_taxid = m_dict[taxid]
                a_dict[acc] = new_taxid

            # sys.stdout.write('\rCheck for merged taxids... {}%'.format(round(counter / dict_size * 100, 1)))
            # sys.stdout.flush()
        end_time = time()
        interval = end_time - start_time
        # sys.stdout.write(' ({})\n'.format(self.elapsed_time(interval)))
        print(' ({})'.format(self.elapsed_time(interval)))

    def make_dict(self, a_dict, t_dict):
        """

        :param l:
        :param acc_dict:
        :param taxo_dict:
        :return: Dictionary with tuple for each accession number with (taxid, rank, name)
        """

        print("Fetching taxid for {}".format(os.path.basename(self.accession)), end="", flush=True)
        start_time = time()
        for acc in self.acc_dict:
            if acc in a_dict:
                taxid = a_dict[acc]
                rank = self.nodes_dict[taxid]['rank']
                name = self.names_dict[taxid]
                t_dict[acc] = [(taxid, rank, name)]
            else:
                t_dict[acc] = [('NA', 'NA', 'NA')]
        end_time = time()
        interval = end_time - start_time
        print(' ({})'.format(self.elapsed_time(interval)))

    def make_lineage(self):
        """
        Format taxonomy for QIIME 1
        :param d: The final dictioany
        :param f: the output file, tsv format
        """

        print('\nMaking lineages...', end="", flush=True)
        start_time = time()

        ranks_of_interest = [b'kingdom', b'phylum', b'class', b'order', b'family', b'genus', b'species']

        for acc, list_of_tuple in self.taxo_dict.items():
            # print("acc {}: ".format(acc.decode('ascii')))
            if list_of_tuple[0][0] != 'NA':
                taxid = list_of_tuple[0][0]  # fist element of first tuple
                # print("{}({}) ".format(self.names_dict[taxid].decode('ascii'),
                #                        taxid.decode('ascii')), end="", flush=True)
                if taxid in self.nodes_dict.keys():
                    parent_taxid = self.nodes_dict[taxid]['parent']
                    if parent_taxid in self.nodes_dict.keys():
                        while parent_taxid != b'131567':  # cellular organism
                            rank = self.nodes_dict[parent_taxid]['rank']
                            name = self.names_dict[parent_taxid]
                            if rank in ranks_of_interest:
                                # print("{}({}) ".format(name.decode('ascii'),
                                #       parent_taxid.decode('ascii')), end="", flush=True)
                                self.taxo_dict[acc].insert(0, (parent_taxid, rank, name))
                            # Check new parent
                            parent_taxid = self.nodes_dict[parent_taxid]['parent']
                        # print('\n')
                    else:
                        print('Parent taxid {} not found in "nodes.dmp"'.format(parent_taxid))
                else:
                    print('Taxid {} not found in "nodes.dmp"'.format(taxid))

        end_time = time()
        interval = end_time - start_time
        print(' ({})'.format(self.elapsed_time(interval)))

    def make_ouput(self, d, f):

        print('\nMaking output file...', end="", flush=True)
        start_time = time()

        rank_shorts = {'kingdom': 'k__',
                       'phylum': 'p__',
                       'class': 'c__',
                       'order': 'o__',
                       'family': 'f__',
                       'genus': 'g__',
                       'species': 's__'}

        with open(f, 'w') as file:
            # Write header
            file.write("Accession\tTaxonomy\n")
            for acc, list_of_tuple in d.items():
                if acc == b'MG243622.1':  # because it's a subspecies!
                    test = 1
                file.write(acc.decode('ascii') + "\t")
                if list_of_tuple[0][0] != 'NA':
                    # Build taxonomy string
                    taxo_list = list()
                    for t in list_of_tuple:
                        try:  # did that because sometime the root taxid is a subspecies instead of a species
                            rank = rank_shorts[t[1].decode('ascii')]
                        except KeyError:
                            continue
                        name = t[2].decode('ascii').replace(' ', '_')
                        taxo_list.append(rank + name)
                    file.write("{}\n".format('; '.join(taxo_list)))
                else:
                    file.write("NA\n")

        end_time = time()
        interval = end_time - start_time
        print(' ({})'.format(self.elapsed_time(interval)))

    def elapsed_time(self, seconds):
        """
        Transform a time value into a string
        :param seconds: Elapsed time in seconds
        :return: Formated time string
        """

        import numpy as np

        minutes, seconds = divmod(seconds, 60)
        hours, minutes = divmod(minutes, 60)
        days, hours = divmod(hours, 24)
        periods = [('d', days), ('h', hours), ('m', minutes), ('s', seconds)]
        time_string = ''.join('{}{}'.format(int(np.round(value)), name) for name, value in periods if value)
        return time_string


if __name__ == '__main__':

    from argparse import ArgumentParser

    parser = ArgumentParser(description='Reorder assembly to have the dnaA gene first')
    parser.add_argument('-t', '--taxo', metavar='/taxonomy/',
                        required=False,
                        help='A folder where NCBI\'s taxonomy files are located'
                             '(taxdump.tar.gz and nucl_gb.accession2taxid.gz)')
    parser.add_argument('-a', '--accession', metavar='accession.list',
                        required=True,
                        help='A text file with one accession number per line')
    parser.add_argument('-o', '--output', metavar='taxonomy.tsv',
                        required=False,
                        help='Output taxonomy file. Will default to input list location if omitted.')

    # Get the arguments into an object
    arguments = parser.parse_args()

    FormatTaxo(arguments)
    
