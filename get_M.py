#!/usr/bin/env python
# -*- coding: utf-8 -*-

import ldscore.ldscore as ld
import ldscore.parse as ps
import ldscore.sumstats as sumstats
import numpy as np
import pandas as pd
import argparse
import time
from functools import reduce


def sec_to_str(t):
    '''Convert seconds to days:hours:minutes:seconds'''
    [d, h, m, s, n] = reduce(lambda ll, b : divmod(ll[0], b) + ll[1:], [(t, 1), 60, 60, 24])
    f = ''
    if d > 0:
        f += '{D}d:'.format(D=d)
    if h > 0:
        f += '{H}h:'.format(H=h)
    if m > 0:
        f += '{M}m:'.format(M=m)

    f += '{S}s'.format(S=s)
    return f


class Logger(object):
    '''
    Lightweight logging.
    '''
    def __init__(self, fh):
        self.log_fh = open(fh, 'w')

    def log(self, msg):
        '''
        Print to log file and stdout with a single command.
        '''
        print(msg, file=self.log_fh)
        print(msg)


def __filter__(fname, noun, verb, merge_obj):
    merged_list = None
    if fname:
        f = lambda x,n: x.format(noun=noun, verb=verb, fname=fname, num=n)
        x = ps.FilterFile(fname)
        c = 'Read list of {num} {noun} to {verb} from {fname}'
        print(f(c, len(x.IDList)))
        merged_list = merge_obj.loj(x.IDList)
        len_merged_list = len(merged_list)
        if len_merged_list > 0:
            c = 'After merging, {num} {noun} remain'
            print(f(c, len_merged_list))
        else:
            error_msg = 'No {noun} retained for analysis'
            raise ValueError(f(error_msg, 0))

        return merged_list


def get_m(args, log):
    '''
    Calculates M and M_5_50 without calculating LD Scores.
    '''
    if args.bfile:
        snp_file, snp_obj = args.bfile+'.bim', ps.PlinkBIMFile
        ind_file, ind_obj = args.bfile+'.fam', ps.PlinkFAMFile
        array_file, array_obj = args.bfile+'.bed', ld.PlinkBEDFile

    # read bim/snp
    array_snps = snp_obj(snp_file)
    m = len(array_snps.IDList)
    log.log('Read list of {m} SNPs from {f}'.format(m=m, f=snp_file))
    if args.annot is not None:  # read --annot
        try:
            if args.thin_annot: # annot file has only annotations
                annot = ps.ThinAnnotFile(args.annot)
                n_annot, ma = len(annot.df.columns), len(annot.df)
                log.log("Read {A} annotations for {M} SNPs from {f}".format(f=args.annot,
                    A=n_annot, M=ma))
                annot_matrix = annot.df.values
                annot_colnames = annot.df.columns
                keep_snps = None
            else:
                annot = ps.AnnotFile(args.annot)
                n_annot, ma = len(annot.df.columns) - 4, len(annot.df)
                log.log("Read {A} annotations for {M} SNPs from {f}".format(f=args.annot,
                    A=n_annot, M=ma))
                annot_matrix = np.array(annot.df.iloc[:,4:])
                annot_colnames = annot.df.columns[4:]
                keep_snps = None
                if np.any(annot.df.SNP.values != array_snps.df.SNP.values):
                    raise ValueError('The .annot file must contain the same SNPs in the same'+\n                        ' order as the .bim file.')
        except Exception:
            log.log('Error parsing .annot file')
            raise

    elif args.extract is not None:  # --extract
        keep_snps = __filter__(args.extract, 'SNPs', 'include', array_snps)
        annot_matrix, annot_colnames, n_annot = None, None, 1

    else:
        annot_matrix, annot_colnames, keep_snps = None, None, None
        n_annot = 1

    # read fam
    array_indivs = ind_obj(ind_file)
    n = len(array_indivs.IDList)
    log.log('Read list of {n} individuals from {f}'.format(n=n, f=ind_file))
    # read keep_indivs
    if args.keep:
        keep_indivs = __filter__(args.keep, 'individuals', 'include', array_indivs)
    else:
        keep_indivs = None

    # read genotype array
    log.log('Reading genotypes from {fname}'.format(fname=array_file))
    geno_array = array_obj(array_file, n, array_snps, keep_snps=keep_snps,
        keep_indivs=keep_indivs, mafMin=args.maf)

    # filter annot_matrix down to only SNPs passing MAF cutoffs
    if annot_matrix is not None:
        annot_keep = geno_array.kept_snps
        annot_matrix = annot_matrix[annot_keep,:]

    log.log("Calculating M and M_5_50.")
    file_suffix = "l2"
    if annot_matrix is not None:
        M = np.atleast_1d(np.squeeze(np.asarray(np.sum(annot_matrix, axis=0))))
        ii = geno_array.maf > 0.05
        M_5_50 = np.atleast_1d(np.squeeze(np.asarray(np.sum(annot_matrix[ii,:], axis=0))))
    else:
        M = [geno_array.m]
        M_5_50 = [np.sum(geno_array.maf > 0.05)]

    # print .M
    fout_M = open(args.out + '.'+ file_suffix +'.M','w')
    print('\t'.join(map(str,M)), file=fout_M)
    fout_M.close()
    log.log('Wrote M to {f}'.format(f=args.out + '.'+ file_suffix +'.M'))

    # print .M_5_50
    fout_M_5_50 = open(args.out + '.'+ file_suffix +'.M_5_50','w')
    print('\t'.join(map(str,M_5_50)), file=fout_M_5_50)
    fout_M_5_50.close()
    log.log('Wrote M_5_50 to {f}'.format(f=args.out + '.'+ file_suffix +'.M_5_50'))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--out', required=True, type=str, help='Output filename prefix.')
    parser.add_argument('--bfile', required=True, type=str, help='Prefix for Plink .bed/.bim/.fam file')
    parser.add_argument('--annot', default=None, type=str, help='Filename prefix for annotation file.')
    parser.add_argument('--thin-annot', action='store_true', default=False, help='If your annot files have only annotations, with no SNP, CM, CHR, BP columns.')
    parser.add_argument('--maf', default=None, type=float, help='Minor allele frequency lower bound. Default is MAF > 0.')
    parser.add_argument('--keep', default=None, type=str, help='File with individuals to include.')
    parser.add_argument('--extract', default=None, type=str, help='File with SNPs to include.')

    args = parser.parse_args()
    log = Logger(args.out+'.log')
    log.log('Beginning analysis at {T}'.format(T=time.ctime()))
    start_time = time.time()

    get_m(args, log)

    log.log('Analysis finished at {T}'.format(T=time.ctime()))
    time_elapsed = round(time.time()-start_time, 2)
    log.log('Total time elapsed: {T}'.format(T=sec_to_str(time_elapsed)))
