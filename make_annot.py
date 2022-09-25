#!/usr/bin/env python
from __future__ import print_function
import pandas as pd
import numpy as np
import argparse
from pybedtools import BedTool
import gzip
import os
import sys

def gene_set_to_bed(args):
    print('making gene set bed file')
    GeneSet = pd.read_csv(args.gene_set_file, header = None, names = ['GENE'])
    all_genes = pd.read_csv(args.gene_coord_file, delim_whitespace = True)
    df = pd.merge(GeneSet, all_genes, on = 'GENE', how = 'inner')
    df['START'] = np.maximum(1, df['START'] - args.windowsize)
    df['END'] = df['END'] + args.windowsize
    iter_df = [['chr'+(str(x1).lstrip('chr')), x2 - 1, x3] for (x1,x2,x3) in np.array(df[['CHR', 'START', 'END']])]
    return BedTool(iter_df).sort().merge()


def make_annot_single(df_bim, bimbed, bed_for_annot):
    if bed_for_annot.field_count() == 3:
        ## binary annotation
        annotbed = bimbed.intersect(bed_for_annot)
        score = [int(1) for x in annotbed]
    else:
        ## continuous annotation
        annotbed = bimbed.intersect(bed_for_annot, wb=True)
        score_col = annotbed.field_count()
        score = [float(x.fields[score_col-1]) for x in annotbed]

    bp = [x.start + 1 for x in annotbed]
    chrom = [x.chrom.replace("chr","") for x in annotbed]
    df_int = pd.DataFrame({'CHR': chrom,'BP': bp, 'ANNOT':score})

    df_annot = pd.merge(df_bim, df_int.drop_duplicates(), how='left', on=['CHR','BP'])
    df_annot.fillna(0, inplace=True)

    return df_annot[['ANNOT']].astype(float)


def make_annot_files(args, bed_for_annot, list_input=False, header=None):
    df_bim = pd.read_csv(args.bimfile,
                delim_whitespace=True, usecols = [0,1,2,3], names = ['CHR','SNP','CM','BP'])
    iter_bim = [['chr'+str(x1), x2 - 1, x2] for (x1, x2) in np.array(df_bim[['CHR', 'BP']])]
    bimbed = BedTool(iter_bim)
    if not list_input:
        ## if input is a single BedTool object
        df_annot = make_annot_single(df_bim, bimbed, bed_for_annot)

    else:
        ## if input is a list of BedTool objects
        df_annot=pd.DataFrame()
        for i in range(len(bed_for_annot)):
            sys.stderr.write('making annot file with ' + str(i+1) + '/' + str(len(bed_for_annot)) + 'bed files\n')
            #print('making annot file with ' + str(i+1) + '/' + str(len(bed_for_annot)) + 'bed files')
            df_annot_1 = make_annot_single(df_bim, bimbed, bed_for_annot[i])
            if header:
                df_annot_1 = df_annot_1.rename(columns={'ANNOT': header[i]})
            else:
                df_annot_1 = df_annot_1.rename(columns={'ANNOT': 'ANNOT_'+ str(i)})
            df_annot = pd.concat([df_annot, df_annot_1], axis=1)
    
    if args.baseline:
        N = len(df_bim)
        base = pd.DataFrame({'Base':np.ones(N)}).astype(int)
        df_annot = pd.concat([base, df_annot], axis=1)

    print('writing annot file')
    if args.annot_file.endswith('.gz'):
        with gzip.open(args.annot_file, 'wb') as f:
            df_annot.to_csv(f, sep = "\t", index = False)
    else:
        df_annot.to_csv(args.annot_file, sep="\t", index=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--gene-set-file', type=str, help='a file of gene names, one line per gene.')
    parser.add_argument('--gene-coord-file', type=str, default='ENSG_coord.txt', help='a file with columns GENE, CHR, START, and END, where START and END are base pair coordinates of TSS and TES. This file can contain more genes than are in the gene set. We provide ENSG_coord.txt as a default.')
    parser.add_argument('--windowsize', type=int, help='how many base pairs to add around the transcribed region to make the annotation?')
    parser.add_argument('--bed-file', type=str, help='the UCSC bed file with the regions that make up your annotation')
    parser.add_argument('--nomerge', action='store_true', default=False, help='don\'t merge the bed file; make an annot file with values proportional to the number of intervals in the bedfile overlapping the SNP.')
    parser.add_argument('--bimfile', type=str, help='plink bim file for the dataset you will use to compute LD scores.')
    parser.add_argument('--annot-file', type=str, help='the name of the annot file to output.')
    parser.add_argument('--baseline', action='store_true',  help='Output annotation with baseline (all snps). Can be disabled by [--no-baseline]')
    parser.add_argument('--no-baseline', dest='baseline', action='store_false')
    parser.set_defaults(baseline=True)

    args = parser.parse_args()

    if args.gene_set_file is not None:
        bed_for_annot = gene_set_to_bed(args)
        make_annot_files(args, bed_for_annot)
    else:
        bed_for_annot_list=[]
        annot_header=[]
        for file in args.bed_file.split(","):
            annot = os.path.basename(file).split('.')[0]
            annot_header.append(annot)
            #print('reading bed file for ' + annot )
            bed_for_annot = BedTool(file).sort()
            if not args.nomerge:
                bed_for_annot = bed_for_annot.merge()
            bed_for_annot_list.append(bed_for_annot)

        make_annot_files(args, bed_for_annot_list, True, annot_header)

