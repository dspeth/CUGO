#!/usr/bin/env python

# import required libraries
import argparse
import sys
from pathlib import Path
import pandas as pd

# parse args
parser = argparse.ArgumentParser()

parser.add_argument("-i", "--prot_ids", help="list of protein accession numbers to retrieve genomic context for")
parser.add_argument("-d", "--cugo_dir", help="dir with CUGO files of the genomes containing proteins of interest")
parser.add_argument("-t", "--tmhmm_dir", help="OPTIONAL: dir with (modified) TMHMM files of the genomes containing proteins of interest")
parser.add_argument("-r", "--cugo_range", type=int, default=0, help="range of genomic context, in CUGO steps from target gene. Default 0, just CUGO of target")
parser.add_argument("-V", "--version", action="store_true", help="show script version and exit")
parser.add_argument("-o", "--out_file", help="output file, defined by user")

args = parser.parse_args()

# version supersedes all
if args.version:
    print("parser, version 0.1")
    sys.exit(0)

if not args.prot_ids or not args.cugo_dir or not args.out_file:
    print("specifying an ID list, a CUGO dir, and an Output file is required")
    print("use '-h' for full options")
    sys.exit(0)

# parse args
prot_ids = Path(args.prot_ids)
cugo_dir = Path(args.cugo_dir)
cugo_range = args.cugo_range
out_file = args.out_file
if args.tmhmm_dir:
    tmhmm_dir = Path(args.tmhmm_dir)


# check existence of files/dir
if not prot_ids.is_file():
    print("provided input file does not exist")
    sys.exit(0)

if not cugo_dir.is_dir():
    print("provided cugo dir path is not a directory")
    sys.exit(0)
if tmhmm_dir and not tmhmm_dir.is_dir():
    print("provided tmhmm dir path is not a directory")
    sys.exit(0)

count = 0
with open(prot_ids, "r") as id_list:
    for line in id_list:
        count += 1
        target_ID = line.strip()
        genome_ID = target_ID.rsplit("_", 1)[0]

        genome_cugo_df = pd.read_csv(str(cugo_dir)+"/"+genome_ID+"_cugo.tab",  sep="\t", na_filter=False)
        if tmhmm_dir:
            genome_tmhmm_df = pd.read_csv(str(tmhmm_dir)+"/"+genome_ID+"_tmhmm_clean",  sep="\t", na_filter=False)
            parse_df = pd.merge(genome_cugo_df, genome_tmhmm_df, on="prot_ID", how="left")
        else:
            parse_df = genome_cugo_df

        target_select = parse_df[parse_df["prot_ID"] == target_ID]
        target_cugo = target_select["CUGO_number"].item()
        target_parent = target_select["parent_ID"].item()
        target_strand = target_select["strand"].item()

        parent_df = parse_df[parse_df["parent_ID"] == target_parent]

        cugo_context = parent_df[(parent_df["CUGO_number"] >= (target_cugo - cugo_range)) & (parent_df["CUGO_number"] <= (target_cugo + cugo_range)) ]
        if target_strand == "-":
            cugo_context = cugo_context.iloc[::-1]
        cugo_context = cugo_context.reset_index(drop=True)
        target_index = cugo_context.index[cugo_context.prot_ID == target_ID].item()
        cugo_context.index = cugo_context.index - target_index
        cugo_context = cugo_context.transpose()

        if count != 1:
            temp = cugo_context_all
            cugo_context_all = pd.concat([temp, cugo_context], join="outer", axis=0)
        else:
            cugo_context_all = cugo_context

cugo_context_all = cugo_context_all.reset_index().rename(columns={"index" : "feat_type"})
cugo_context_all.to_csv(args.out_file, sep="\t", index=False)
