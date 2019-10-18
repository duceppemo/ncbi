#!/usr/local/env python3

__author__ = 'duceppemo'
__version__ = '0.1'

from argparse import ArgumentParser
from Bio import SearchIO


class FilterBLastXML(object):
    def __init__(self, args):
        # Command line arguments
        self.input_file = args.input
        self.output_file = args.output
        self.sim = args.similarity

        # Data
        self.search_dict = dict()
        self.filtered_list = list()

        # Run
        self.run()

    def run(self):
        self.parse_xml()
        self.filter_xml()
        self.write_filtered_xml()
        print('Done')

    def parse_xml(self):
        print('Parsing blast xml file...')
        qresults = SearchIO.parse(self.input_file, 'blast-xml')
        print('Converting to dictionary...')
        self.search_dict = SearchIO.to_dict(qresults)

    def similarity_filter(self, hit):
        aln = hit.hsps[0].aln_annotation['similarity']
        total_aln_len = len(aln)
        total_matches = aln.count('|')
        sim_perc = round(total_matches / total_aln_len * 100, 1)
        if sim_perc >= self.sim:
            return hit

    def filter_xml(self):
        print('Filtering results to 99% similarity...')
        new_list = list()
        for ID, qresult in self.search_dict.items():
            filtered = qresult.hit_filter(self.similarity_filter)
            if filtered.hits:
                self.filtered_list.append(filtered)

    def write_filtered_xml(self):
        print('Writing output file...')
        SearchIO.write(self.filtered_list, self.output_file, format='blast-xml')


if __name__ == '__main__':
    parser = ArgumentParser(description='Filter a blast xml file to keep the hits with a minimum % similarity'
                                        'Much quicker to do than rerun the whole blast')
    parser.add_argument('-i', '--input', metavar='sample.blast90.xml',
                        required=True,
                        help='Blast xml output file ("-outfmt 5")')
    parser.add_argument('-o', '--output', metavar='sample.blast99.xml',
                        required=True,
                        help='Filtered blast xml output file')
    parser.add_argument('-s', '--similarity', metavar='99.0', type=float,
                        required=True,
                        help='% similarity')

    # Get the arguments into an object
    arguments = parser.parse_args()

    FilterBLastXML(arguments)
