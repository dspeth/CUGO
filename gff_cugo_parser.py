#!/usr/bin/env python

# import required libraries
import argparse
import sys
from pathlib import Path
import pandas as pd

# parse args
parser = argparse.ArgumentParser()

parser.add_argument("-g", "--gff_file", help="anvi'o generated GFF file to parse")
parser.add_argument("-V", "--version", action="store_true", help="show script version and exit")
parser.add_argument("-o", "--out_file", help="Define output file")

args = parser.parse_args()

# check existence of input files
if args.version:
    print("parser to extract CUGO information from anvi'o generated GFF files, version 0.1")
    print("CUGO is shorthand for colocated unidirectional gene organizaton, because the term gene cluster has become ambiguous")
    sys.exit(0)

if not args.gff_file:
    print("input gff file is required, specify with '-g' or '--gff_file'")
    sys.exit(0)
gff = Path(args.gff_file)

if not args.out_file:
    print("output file is required, specify with '-o' or '--out_file'")
    sys.exit(0)
outfile = Path(args.out_file)

if not gff.is_file():
    print("provided input GFF file does not exist")
    sys.exit(0)


# code to loop over lines of an anvio exported GFF file and write a tab delimited file with CUGO information
# CUGO is shorthand for "colocated unidirectional gene organizaton"

# define variables needed in loop
reformat_data = []
prev_line = None
prev_direction = None
prev_parent = None
prev_feat_type = None
cugo_count = 0
cugo_size = {}
cugo_size_count = 0

# Loop over all lines in file
with open(gff, "r") as GFF:
    for line in GFF:
        clean_line = line.strip().split("\t")

        # if a line is does not consist of 9 tab delimited fields, ignore it
        if len(clean_line) != 9:
            continue

        # if the feature type is not CDS (but e.g. tRNA or rRNA), ignore the line
        feat_type = clean_line[2]
        if feat_type != "CDS":
            if prev_feat_type != "CDS":
                continue
        prev_feat_type = feat_type

        # the loop needs to be aware of the lines before and after the line processed
        # my solution was to store the split line in a variable and process the stored variable at the next iteration of the loop
        if prev_line == None:
            prev_line = clean_line
        else:
            next_parent = clean_line[0]
            next_direction = clean_line[6]

            # retrieve the useful fields of the GFF as variables
            annotation = prev_line[8].split(";")
            seqID = annotation[0].split("=")[1].replace("___", "_")
            if len(annotation) == 1:
                COG_ID = "NA"
            else:
                COG_ID = annotation[1].split("=")[1]
            parent = prev_line[0]
            direction = prev_line[6]
            gene_start = prev_line[3]
            gene_end = prev_line[4]
            nuc_length = abs(int(gene_end) - int(gene_start)) + 1
            aa_length = int(nuc_length / 3)

            # genes process genes in the middle of a CUGO, where the value for both edges is "NA"
            if (direction == prev_direction == next_direction) and (parent == prev_parent == next_parent):
                cugo = prev_cugo
                cugo_size_count += 1
                cugo_start = "NA"
                cugo_end = "NA"

            # process genes at the start of a contig/scaffold
            elif parent != prev_parent:
                cugo = cugo_count
                cugo_size_count = 1
                if direction == "+":
                    cugo_start = "sequence_edge"
                    if parent != next_parent:
                        cugo_end = "sequence_edge"
                        cugo_size[cugo] = cugo_size_count
                        cugo_size_count = 0
                        cugo_count += 1
                    elif direction != next_direction:
                        cugo_end = "strand_change"
                        cugo_size[cugo] = cugo_size_count
                        cugo_size_count = 0
                        cugo_count += 1
                    else:
                        cugo_end = "NA"
                else:
                    cugo_end = "sequence_edge"
                    if parent != next_parent:
                        cugo_start = "sequence_edge"
                        cugo_size[cugo] = cugo_size_count
                        cugo_size_count = 0
                        cugo_count += 1
                    elif direction != next_direction:
                        cugo_start = "strand_change"
                        cugo_size[cugo] = cugo_size_count
                        cugo_size_count = 0
                        cugo_count += 1
                    else:
                        cugo_start = "NA"
                prev_parent = parent
                prev_direction = direction
                prev_cugo = cugo

            # process genes at the end of a scaffold
            elif parent != next_parent:
                cugo = cugo_count
                cugo_count += 1

                cugo_size_count += 1
                cugo_size[cugo] = cugo_size_count
                cugo_size_count = 0

                if direction == "+":
                    cugo_end = "sequence_edge"
                    if direction != prev_direction:
                        cugo_start = "strand_change"
                    else:
                        cugo_start = "NA"
                else:
                    cugo_start = "sequence_edge"
                    if direction != prev_direction:
                        cugo_end = "strand_change"
                    else:
                        cugo_end = "NA"

                prev_parent = parent
                prev_direction = direction
                prev_cugo = cugo

            # process genes at the beginning of a CUGO because of strand change
            elif direction != prev_direction:
                cugo = cugo_count
                cugo_size_count = 1
                if direction == "+":
                    cugo_start = "strand_change"
                    if direction != next_direction:
                        cugo_end = "strand_change"
                        cugo_size[cugo] = cugo_size_count
                        cugo_size_count = 0
                        cugo_count += 1
                    else:
                        cugo_end = "NA"
                else:
                    cugo_end = "strand_change"
                    if direction != next_direction:
                        cugo_start = "strand_change"
                        cugo_size[cugo] = cugo_size_count
                        cugo_size_count = 0
                        cugo_count += 1
                    else:
                        cugo_start = "NA"

                prev_parent = parent
                prev_direction = direction
                prev_cugo = cugo

            # process genes at the end of a CUGO because of strand change
            elif direction != next_direction:
                cugo = cugo_count
                cugo_count += 1

                cugo_size_count += 1
                cugo_size[cugo] = cugo_size_count
                cugo_size_count = 0

                if direction == "+":
                    cugo_end = "strand_change"
                    cugo_start = "NA"
                else:
                    cugo_start = "strand_change"
                    cugo_end = "NA"
                prev_parent = parent
                prev_direction = direction
                prev_cugo = cugo

            # save the relevant parameters in a list of lists, to convert to a dataframe at the end of the loop
            reformat_line = [seqID, parent, gene_start, gene_end, nuc_length, aa_length, direction, COG_ID, cugo, cugo_start, cugo_end]
            reformat_data.append(reformat_line)

            # store the current line for processing next
            prev_line = clean_line

    # get the last line, if needed
    clean_line = line.strip().split("\t")
    if len(clean_line) == 9:
        feat_type = clean_line[2]
        if feat_type == "CDS":
            annotation = clean_line[8].split(";")
            seqID = annotation[0].split("=")[1].replace("___", "_")
            if len(annotation) == 1:
                COG_ID = "NA"
            else:
                COG_ID = annotation[1].split("=")[1]
            parent = clean_line[0]
            direction = clean_line[6]
            gene_start = clean_line[3]
            gene_end = clean_line[4]
            nuc_length = abs(int(gene_end) - int(gene_start)) + 1
            aa_length = int(nuc_length / 3)

            if (direction == prev_direction) and (parent == prev_parent):
                cugo = prev_cugo
                cugo_size_count += 1
                cugo_size[cugo] = cugo_size_count
                if direction == "+":
                    cugo_start = "NA"
                    cugo_end = "sequence_edge"
                else:
                    cugo_start = "sequence_edge"
                    cugo_end = "NA"

            elif parent != prev_parent:
                cugo = cugo_count
                cugo_size_count = 1
                cugo_size[cugo] = cugo_size_count
                cugo_start = "sequence_edge"
                cugo_end = "sequence_edge"

            elif direction != prev_direction:
                cugo = cugo_count
                cugo_size_count = 1
                cugo_size[cugo] = cugo_size_count
                if direction == "+":
                    cugo_start = "strand_change"
                    cugo_end = "sequence_edge"
                else:
                    cugo_start = "sequence_edge"
                    cugo_end = "strand_change"

            reformat_line = [seqID, parent, gene_start, gene_end, nuc_length, aa_length, direction, COG_ID, cugo, cugo_start, cugo_end]
            reformat_data.append(reformat_line)

    # Finally, convert the list of lists to a dataframe, add the CUGO size dictionary, and print to file
    gff_cugo = pd.DataFrame(reformat_data, columns=["seqID", "parent_ID", "gene_start", "gene_end", "nuc_length", "aa_length", "strand", "COG_ID", "CUGO_number", "CUGO_start", "CUGO_end"])
    gff_cugo["CUGO_size"] = gff_cugo["CUGO_number"].map(cugo_size)
    gff_cugo.to_csv(outfile, sep="\t", index=False)
