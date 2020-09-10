#!/usr/local/env python3

from sys import argv
import requests

input_file = argv[1]
output_folder = argv[2]

# NZ_AABUZE000000000.2
# https://sra-download.ncbi.nlm.nih.gov/traces/wgs01/wgs_aux/AA/BU/ZE/AABUZE02/AABUZE02.1.fsa_nt.gz

base = 'https://sra-download.ncbi.nlm.nih.gov/traces/wgs01/wgs_aux/'

with open(input_file, 'r') as f:
    for line in f:
        line = line.rstrip()
        if not line:
            continue
        if '_' in line:
            line = line.split('_')[1]

        f1 = line[:2]
        f2 = line[2:4]
        f3 = line[4:6]
        f4 = line[:7] + line[-1]
        url = "{}{}/{}/{}/{}/{}.1.fsa_nt.gz".format(base, f1, f2, f3, f4, f4)
        print('Downloading {}...'.format(line))
        r = requests.get(url, allow_redirects=True)

        # Download
        open(output_folder + '/' + f4 + '1.fsa_nt.gz', 'wb').write(r.content)
