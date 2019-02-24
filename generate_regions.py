#!/usr/bin/env python3

import yaml
import subprocess

def main():
    # TODO handle IO errors
    config = yaml.load(open('config.yaml', 'r'))
    faidx_filename = config['faidx_ref']
    chunk_size = int(config['chunk_size'])

    with open(faidx_filename, 'r') as faidx_file:
        for line in faidx_file:
            seq_name, seq_length = map(int, line.strip().split()[:2])
            for i, start in enumerate(range(1, seq_length, chunk_size)):
                with open('variant-calls/regions/{}.reg'.format(i+1), 'w') as out:
                    print('{}:{}-{}'.format(seq_name, start, min(start +
                        chunk_size - 1, seq_length)), file=out)

if __name__ == '__main__':
    main()
