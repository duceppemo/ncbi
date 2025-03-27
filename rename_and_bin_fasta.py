import os
from argparse import ArgumentParser
import pathlib
import re


__author__ = 'duceppemo'
__version__ = '0.1'


class Binner(object):
    def __init__(self, args):
        # I/O
        self.input_folder = os.path.abspath(args.input_folder)
        self.output_folder = os.path.abspath(args.output_folder)
        self.file_extensions = args.extensions

        # Data
        self.fasta_list = list()
        self.bin_dict = dict()

        # Run
        self.run()

    def run(self):
        # List fasta
        for ext in self.file_extensions:
            self.fasta_list += Binner.list_fasta(self.input_folder, ext)

        # Create output folder
        Binner.make_folder(self.output_folder)

        # Get bins into dictionary
        for fasta in self.fasta_list:
            Binner.rename_fasta(fasta, self.output_folder)

    @staticmethod
    def list_fasta(input_folder, ext):
        fasta_list = list()
        for root, directories, filenames in os.walk(input_folder):
            for filename in filenames:
                if filename.endswith(tuple(ext)):  # accept a tuple or string
                    file_path = os.path.join(root, filename)
                    fasta_list.append(file_path)
        return fasta_list

    @staticmethod
    def make_folder(folder):
        # Will create parent directories if they don't exist and will not return error if already exists
        pathlib.Path(folder).mkdir(parents=True, exist_ok=True)

    @staticmethod
    def rename_fasta(fasta, output_folder):
        if not os.path.exists(fasta):
            return
        
        with open(fasta, 'r') as f:
            header = f.readline().rstrip()  # Only read the first line

        # Try to clean the header a bit
        header = header.replace('sp.', 'sp').replace('ssp.', 'ssp').replace('subsp.', 'subsp').replace('bv.', 'bv').replace(',', '').replace('/', '-')
        header = re.sub(r'\W+>.', ' ', header)  # remove special characters
        header = re.sub(r'[:\[\]]+','', header)  # remove more special characters not included in "\W"

        # Split header into a list of words
        header_list = header.split()

        acc = ''.join(header_list[0][1:]) + '_'  # Drop the heading ">"
        genus = header_list[1] + '_'
        species = header_list[2]
        subsp = ''
        variant = ''
        biovar= ''
        strain = ''
        isolate = ''
        header_end = ' '.join(header_list[3:])

        # Check for specific taxonomic information from fasta header

        # Prefixes
        if 'TPA_asm' in header_list[1:]:  # Skip accession
            # Find its position
            i = header_list.index('TPA_asm')
            genus = 'TPA_asm_' + header_list[i + 1] + '_'
            species = header_list[i + 2]

        if 'MAG' in header_list[1:]:
            i = header_list.index('MAG')
            genus = 'MAG_' + header_list[i + 1] + '_'
            species = header_list[i + 2]

        if all(x in header_list[1:] for x in ['MAG', 'TPA_asm']):
            i = header_list.index('TPA_asm')
            genus = 'MAG_TPA_asm_' + header_list[i + 1] + '_'
            species = header_list[i + 2]

        if 'UNVERIFIED_CONTAM' in header_list[1:]:
            i = header_list.index('UNVERIFIED_CONTAM')
            genus = 'UNVERIFIED_CONTAM_' + header_list[i + 1] + '_'
            species = header_list[i + 2]

        if 'Candidatus'in header_list[1:]:
            i = header_list.index('Candidatus')
            genus = 'Candidatus_' + header_list[i + 1] + '_'
            species = header_list[i + 2]

        # Suffixes
        if any(x in header_list[1:] for x in ['subsp', 'ssp']):
            try:
                i = header_list.index('subsp')
            except ValueError:
                i = header_list.index('ssp')
            subsp = '_subsp_' + header_list[i + 1]

        if 'variant' in header_list[1:]:
            # Find its position
            i = header_list.index('variant')
            variant = '_variant_' + header_list[i + 1]

        # Suffixes
        if any(x in header_list[1:] for x in ['biovar', 'bv']):
            try:
                i = header_list.index('biovar')
            except ValueError:
                i = header_list.index('bv')
            subsp = '_biovar_' + header_list[i + 1]

        if 'strain' in header_list[1:]:
            i = header_list.index('strain')
            strain = '_strain_' + header_list[i + 1]

        if 'isolate' in header_list[1:]:
            i = header_list.index('isolate')
            isolate = '_isolate_' + header_list[i + 1]

        # Generate new name
        new_name = '{}{}{}{}{}{}{}{}'.format(acc, genus, species, subsp, variant, biovar, strain, isolate)

        # Generate bin name
        bin_name = '{}{}{}{}{}'.format(genus, species, subsp, variant, biovar)
        my_bin = output_folder + '/' + bin_name

        # Create folder
        Binner.make_folder(my_bin)

        # Rename and move file
        try:
            os.rename(fasta, my_bin + '/' + new_name + '.fasta')
        except FileNotFoundError:
            print('{} not found.'.format(fasta))
            # continue


if __name__ == "__main__":
    parser = ArgumentParser(description='Rename and bin fasta files by species based on inforamtion provided in the 1st fasta header')
    parser.add_argument('-i', '--input-folder', metavar='/path/to/input_folder',
                        required=True, type=str,
                        help='Folder that contains the fasta files. Mandatory.')
    parser.add_argument('-o', '--output-folder', metavar='/path/to/bin_folder',
                        required=True, type=str,
                        help='Folder that will hold the subfolders (bins). Mandatory')
    parser.add_argument('-e', '--extensions', metavar='fasta',
                        required=False,
                        default=['fa', 'fasta', 'fna'],
                        help='Limit binning to specific file extensions. Useful if your files are gzipped. \
                        You can use multiple comma-separated extensions (no space). \
                        Defaul is "fa,fasta,fna" . Optional.')

    # Get the arguments into an object
    arguments = parser.parse_args()

    Binner(arguments)
