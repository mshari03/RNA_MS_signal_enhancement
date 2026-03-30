#!/usr/bin/python3

"""
Last update: February 2025
Author: Max Sharin
Contact info: max.sharin@tufts.edu
GitHub project repository: https://github.com/ldascenzo/pytheas

***DESCRIPTION***
Optional step of the Pytheas in silico digest library generation. All digestion fragments are labeled with a 3' terminal oligo before 
mass information is calculated.

***OUTPUT***
-> output.4.MS2 -> the new output file used with appended labeled sequences. (unlabeled sequences are maintained)
"""

import os
import sys
import pandas as pd
import argparse
from Bio import SeqIO
from shutil import move


parser = argparse.ArgumentParser(description='List of available options')
parser.add_argument('--label_sequences', nargs='*', required=True,
                    help='Input RNA label sequence file(s). Please use fasta format and input the file names without '
                         '"=" after the option. First string of the header for each sequence will be used as id')
args = parser.parse_args()


def read_csv(input_csv='./nts_light.csv'):
    """
    Create a dictionary containing all ID : ID_ext couples, one letter and extended IDs for all nucleotides
    """
    if not os.path.isfile(input_csv):
        print("ERROR! File {} with info on nucleotides from script 2_modify.py is missing."
            "Execution terminated without generating decoys".format(input_csv))
        sys.exit(1)

    else:
        # Read the csv file with the nucleotides dictionary
        df = pd.read_csv(input_csv, usecols=['ID', 'ID_ext'])

        # Drop rows with NaN values
        df = df[pd.notnull(df['ID'])]

        return dict(zip(df.ID, df.ID_ext))
    
def get_label(args):
    """
    Open the label fasta file and return a list of seq objects that contains label/s supplied
    """
    # Separate the input files if multiple fasta files are selected, based on the running OS
    labels = []
    for fasta_file in args.label_sequences:
        # First character in file name cannot be whitespace
        if fasta_file[0].isdigit() or fasta_file[0].isalpha():
            with open(os.getcwd() + "/" + fasta_file.rstrip(), 'r') as handle:
                # Extract and process the input fasta sequences
                for seq in SeqIO.parse(handle, "fasta"):
                    with open(fasta_file.rstrip(), 'r') as handle:
                        # Extract and process the input fasta sequences
                        for seq in SeqIO.parse(handle, "fasta"):
                            labels.append(seq)
    return labels


def label_lines(input_file, labels):
    #TODO: it might be helpful to include information about which label is added to each sample
    # this would be beneficial if differential labeling is used for quantitation
    """
    Generate a list with all the labeled lines to be appended on the original file
    Each labeled line is obtained taking line per line the sequences in the input file
    and adding the 3' terminal label. 
    """
    outlines = []

    # Cycle among the lines of the input file
    for label in labels:
        with open(input_file, 'r') as infile:
            for line in infile:
                # Lines of the header are excluded
                if line[0] != '#' and line.split()[0] != 'sequence' and 'decoy' not in line:
                    seq, chem3, chem5 = line.split()[0], line.split()[3], line.split()[4]
                    # # Only sequences at least 3 nts long are considered for label addition
                    # if len(seq) > 2:
                    entry = list(seq)
                    for nt in label:
                        entry = entry + [nt]
                    new_seq = ''.join(entry)
                    new_molecule_ID = str(line.split()[7]) + "-tagged," + str(line.split()[8]) + "," + str(line.split()[9])
                    line_start = str("{} {} {}".format(new_seq, ' '.join(line.split()[1:6]), new_molecule_ID))
                    line_end = str("{}".format(' '.join(line.split()[8:10])))
                    outlines.append(
                                        "{} {} {}\n".format(line_start, "tagged", line_end))
                    #             break
        # add label oligo to list
        label_line_1  = "{} {} {} {} {} {} {} {} {} {}\n".format(label.seq, "-", "label", "OH", "OH", "1", "0", "label", "0", "0")
        outlines.append(label_line_1)
        label_line_2  = "{} {} {} {} {} {} {} {} {} {}\n".format(label.seq, "-", "label", "OH", "P", "1", "0", "label", "0", "0")
        outlines.append(label_line_2)
        label_line_3  = "{} {} {} {} {} {} {} {} {} {}\n".format(label.seq, "-", "label", "OH", "App", "1", "0", "label", "0", "0")
        outlines.append(label_line_3)
        
    return outlines


def output_lines(input_file, labels_lines):
    """
    Merge the lines of the input and the decoy, ready to be written in the output file
    """
    # Check if the input file output.3.MS2 is present in the directory
    if not os.path.isfile("./" + input_file):
        print("ERROR! Input file output.3.MS2 is missing. Execution terminated without adding labels")
        sys.exit(1)

    with open("./output.3.MS2", 'r') as infile:
        input_lines = infile.readlines()

    return input_lines + labels_lines


if __name__ == "__main__":
    open('./output.4.MS2', 'w').writelines(output_lines('output.3.MS2', label_lines('output.3.MS2',
                                                                                  get_label(args))))

    print("Done! Output file(s) -> output.4.MS2")
