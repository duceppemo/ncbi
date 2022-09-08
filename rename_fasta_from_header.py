import glob
import os
import re


input_folder = '/media/3tb/db/MbovisID/susan_wgs_2'

# Fasta list
fasta_list = glob.glob(input_folder + '/*.fasta')

for fasta in fasta_list:
    file_name = '.'.join(os.path.basename(fasta).split('.')[:-1])
    if file_name.startswith('GCA_'):
        file_name = '_'.join(file_name.split('_')[:2])

    with open(fasta, 'r') as f:
        header = f.readline().rstrip()

    header = header.replace('sp.', '').replace(',', '').replace('/', '-')
    header = re.sub(r'\W+>.', ' ', header)  # remove special characters

    header_list = header.split()
    acc = ''.join(header_list[0][1:])
    genus = header_list[1]
    species = header_list[2]
    subsp = ''
    variant = ''
    strain = ''
    header_end = ' '.join(header_list[3:])

    if 'variant' in header:
        variant = header_list[4]
        if 'strain' in header:
            strain = header_list[6]
            new_name = '{}_{}_variant_{}_strain_{}_({})'.format(genus, species, variant, strain, file_name)
        elif 'isolate' in header:
            strain = header_list[6]
            new_name = '{}_{}_variant_{}_isolate_{}_({})'.format(genus, species, variant, strain, file_name)
        else:
            try:
                strain = header_end.split()[2]
                # strain = re.split(r'strain|isolate|genome|contig|chromosome|NODE', header_end)[1].strip().split()[0]
            except IndexError:
                print('ERROR: ' + header)

            new_name = '{}_{}_variant_{}_strain_{}_({})'.format(genus, species, variant, strain, file_name)
    elif any(x in header for x in ['subsp', 'sub']):
        subsp = header_list[4]
        if 'strain' in header:
            strain = header_list[6]
            new_name = '{}_{}_subsp_{}_strain_{}_({})'.format(genus, species, subsp, strain, file_name)
        elif 'isolate' in header:
            strain = header_list[6]
            new_name = '{}_{}_subsp_{}_isolate_{}_({})'.format(genus, species, subsp, strain, file_name)
        else:
            try:
                strain = header_end.split()[2]
            except IndexError:
                print('ERROR: ' + header)

            new_name = '{}_{}_subsp_{}_strain_{}_({})'.format(genus, species, subsp, strain, file_name)
    else:
        if 'strain' in header:
            strain = header_list[4]
            new_name = '{}_{}_strain_{}_({})'.format(genus, species, strain, file_name)
        elif 'isolate' in header:
            strain = header_list[4]
            new_name = '{}_{}_isolate_{}_{}'.format(genus, species, strain, file_name)
        else:
            strain = header_list[3]
            # strain = re.split(r'strain|isolate|genome|contig|chromosome|NODE', header_end)[1].strip().split()[0]
            new_name = '{}_{}_strain_{}_({})'.format(genus, species, strain, file_name)

    # print(new_name)
    os.rename(fasta, os.path.dirname(fasta) + '/' + new_name + '.fasta')
