#!/usr/bin/env python
# coding: utf-8
"""
@author      CS Lim
@create date 2020-09-15 17:40:16
@modify date 2020-09-15 20:33:55
@desc        RIBOSS pipline
"""

import os
import re
import sys
import time
import argparse
import textwrap
import csv
import random
from itertools import product
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import scipy.stats as stats
import pysam


def check_arg(args=None):
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
    description='Automating footprint size selection using Fisher-exact test',
    epilog=textwrap.dedent('''\
      examples:
        python select_footprint_size.py \\
            -b Aligned.out.bam \\
            -p Arabidopsis_thaliana.TAIR10.cdna.all_filtered.txt
        
        python select_footprint_size.py \\
            -b Aligned.out.bam \\
            -p GRCh38_latest_rna_filtered.txt \\
            -e 1

      if you find this useful, please consider citing:
        Lim, C.S. & Brown, C.M. (2020) Synergistic effects of upstream open reading frames and 
            leader exon-exon junctions on protein expression. bioRxiv.
        Lim, C.S., Wardell, S.J.T., Kleffmann, T. & Brown, C.M. (2018) The exon-intron gene 
            structure upstream of the initiation codon predicts translation efficiency.
            Nucleic Acids Res. 46: 4575-4591.
        '''))

    parser.add_argument('-b', '--bam',
                        metavar='STR',
                        help='input transcript-level BAM alignment of Ribo-seq',
                        required='True') 
    parser.add_argument('-p', '--position',
                        metavar='STR',
                        help='CDS range. _filtered.txt from prep_input_files.py.',
                        required='True')
    parser.add_argument('-o', '--output_path',
                        default='./',
                        metavar='STR',
                        help='output directory. Default = ./')              
    parser.add_argument('-d', '--downsampling',
                        type=check_downsampling,
                        default=0.2,
                        metavar='FLOAT',
                        help='fraction for downsampling BAM file. Default = 0.2')
    parser.add_argument('-l', '--length',
                        type=check_size,
                        default=25,
                        metavar='INT',
                        help='minimum footprint size cutoff (inclusive). Default = 25')
    parser.add_argument('-e', '--stringency',
                        type=check_stringency,
                        default=3,
                        metavar='INT',
                        help='footprint selection stringency (0-4). Default = 3')    
    results = parser.parse_args(args)
    return (results.bam, results.position, results.output_path, results.downsampling, results.length, results.stringency)


def check_downsampling(fraction):        
    if 0.1 <= float(fraction) <= 1:
        return float(fraction)
    else:
        raise argparse.ArgumentTypeError('Fraction out of range.')

        
def check_size(size):        
    if 20 <= int(size) <= 30:
        return int(size)
    else:
        raise argparse.ArgumentTypeError('Footprint size cutoff out of range. Recommend 20 to 30')
        
        
def check_stringency(level):        
    if 0 <= int(level) <= 4:
        return int(level)
    else:
        raise argparse.ArgumentTypeError('Footprint selection stringency out of range. Please use an integer from 0 to 4.')
        
        
def io(bamfile, outdir):
    check_dir = os.path.isdir(outdir)
    if not check_dir:
        os.makedirs(outdir)

    try:
        base = os.path.basename(bamfile) # .split(os.extsep)
        if re.search('\.bam$', bamfile):
            footprint = outdir + '/' + os.path.splitext(base)[0] + '_footprint.txt'
            offset = outdir + '/' + os.path.splitext(base)[0] + '_offset.txt'
            bar = outdir + '/' + os.path.splitext(base)[0] + '_offset_footprint.pdf'
    except Exception:
        raise argparse.ArgumentTypeError('Please provide a bam file as input!')
    return footprint, offset, bar

    
def bam_sampling(bamfile, fraction):
    start_time = time.perf_counter()
    print('Downsampling', bamfile, 'to a fraction of', fraction, '...')
    infile = pysam.AlignmentFile(bamfile, 'rb')
    sample = []
    random.seed(123)
    for read in infile:
        if read.mapping_quality==255 and random.random() < float(fraction):
            sample.append([str(read.reference_name),str(read).split('\t')])
    infile.close()
    print('Finished downsampling in', 
          round((time.perf_counter()-start_time)/60), 'min',
          (time.perf_counter()-start_time)%60, 's\n')
    return sample


def footprint_summary(sample, tsv, size, footprint):
    """
    Input from bam_sampling
    """
    start_time = time.perf_counter()
    print('Counting footprints by reading frames and fragment sizes >=', size, 'nt ...')
    df = pd.DataFrame(sample)
    df.columns = ['tid','read']
    pos = pd.DataFrame(df['read'].apply(lambda x: [x[3],x[5]]).tolist(), 
                       columns = ['pos','cigar'])
    df = pd.concat([df['tid'],pos], axis=1)
    df['footprint_len'] = df.cigar.str.split('M', expand=True)[[0]]
    df = df[~df.footprint_len.str.contains('S')]
    df['pos'] = df.pos.astype(int)
    df['footprint_len'] = df.footprint_len.astype(int)

    cds = pd.read_csv(tsv, sep='\t', header=None)
    cds.columns = ['tid','CDS_start','CDS_end']
    df = pd.merge(cds,df,on='tid')
    df['frame'] = (df.pos - df.CDS_start)%3
    df = df[(df.frame>=0) & (df.footprint_len>=int(size))]
    df = df.groupby(['footprint_len','frame']).count().reset_index()[['footprint_len','frame','tid']]
    df.columns = ['footprint_len','frame','counts']
    df.to_csv(footprint, sep='\t', index=None, quoting = csv.QUOTE_NONE, escapechar = ' ')
    print('Finished footprint counting in', 
          round((time.perf_counter()-start_time)/60), 'min',
          (time.perf_counter()-start_time)%60, 's\n')
    return df


def fisher_test(df, s, offset):
    """
    Input from footprint_summary
    """
    print('Comparing the triplet periodicity of different footprint sizes ...')
    dt = df.pivot(index='footprint_len', columns='frame', values='counts').reset_index()
    dt.columns = ['footprint_len','Frame1_ribosome','Frame2_ribosome','Frame3_ribosome']
    dt['Frame2_3_ribosome'] = dt.Frame2_ribosome + dt.Frame3_ribosome

    f = dt[['footprint_len','footprint_len']]
    f.columns = ['footprint_len','footprint_len_']
    uniques = [f[i].unique().tolist() for i in f.columns ]
    f = pd.DataFrame(product(*uniques), columns = f.columns)
    f = pd.merge(f, dt, on='footprint_len')
    f.columns = ['footprint_len_','footprint_len','Frame1_ribosome','Frame2_ribosome','Frame3_ribosome','Frame2_3_ribosome']
    f = pd.merge(f, dt, on='footprint_len')
    f = f[f.footprint_len_ != f.footprint_len]
    f['tab_x'] = f[['Frame1_ribosome_x','Frame2_3_ribosome_x']].values.tolist()
    f['tab_y'] = f[['Frame1_ribosome_y','Frame2_3_ribosome_y']].values.tolist()
    f['tab'] = f[['tab_x','tab_y']].values.tolist()
    f = f.reset_index(drop=True)
    f['Fisher_exact'] = f.tab.apply(lambda x: list(stats.fisher_exact(x, alternative='greater')))
    f['Adjusted_pvalue'] = f.Fisher_exact.apply(lambda x : 1 if x[1]*len(f) > 1 else (x[1]*len(f) if x[1]*len(f) < 0.05 else x[1]))
    f = f[f.Adjusted_pvalue<0.05].groupby('footprint_len_').count().reset_index()[['footprint_len_','tab']]

    f = f[f.tab>len(f)/(5-int(s))][['footprint_len_']]
    f.columns = ['footprint_len']
    f['P-site'] = 12

    f.to_csv(offset, sep='\t', index=None, header=None,
             quoting = csv.QUOTE_NONE, escapechar = ' ')
    return f


def plot_footprints(f, df, bar):
    """
    Input from footprint_summary and fisher_test
    """
    df = pd.merge(f, df, on='footprint_len')
    df.columns = ['Footprint length','P-site','Reading frame', 'Counts']

    sns.set(style="whitegrid")#, font_scale=2)
    sns.catplot('Reading frame', 'Counts', col='Footprint length', 
                    kind="bar", data=df, aspect=0.5)
    plt.savefig(bar)


def main():
    start_time = time.perf_counter()
    footprint,offset,bar = io(b, o)
    s = bam_sampling(b, d)
    df = footprint_summary(s, p, l, footprint)
    f = fisher_test(df, e, offset)
    plot_footprints(f, df, bar)
    print('Finished all analysis in', 
          round((time.perf_counter()-start_time)/60), 'min',
          (time.perf_counter()-start_time)%60, 's \
          \n\nOutput files:\n', footprint, '\n', offset, '\n', bar)

if __name__ == "__main__":
    b, p, o, d, l, e = check_arg(sys.argv[1:])
    main()