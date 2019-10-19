#!/usr/bin/env python

import argparse
import datetime

def main(args):
    """
    Take a fasta file with a deflines like:
    >CY083910_A/Aalborg/INS132/2009_2009/11/23
    ...and parses them to a csv file with three columns (label, clade name, date)
    CY083910,A/Aalborg/INS132/2009,2009-11-23
    """
    with open(args.fasta_input, 'r') as f:
        for line in f:
            if line[0] == '>':
                defline = line[1:].strip()
                if (len(defline.split('_')) == 3):
                    defline_split = defline.split('_')
                # some clade names include '_', so they get split.
                # will reconstruct them
                elif (len(defline.split('_')) > 3):
                    defline_split = defline.split('_')
                    label = defline_split[0]
                    clade_name = "_".join(defline_split[1:-1])
                    date = defline_split[-1]
                    defline_reconstructed = [label, clade_name, date]
                    defline_split = defline_reconstructed
                label = defline_split[0]
                clade_name = defline_split[1]
                date = datetime.datetime.strptime(defline_split[2], '%Y/%m/%d').date()
                print(",".join([label, clade_name, str(date)]))
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("fasta_input", help="Input file (fasta)")
    args = parser.parse_args()
    main(args)
