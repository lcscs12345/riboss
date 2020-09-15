#!/usr/bin/env python
# coding: utf-8
"""
@author      CS Lim
@create date 2020-09-15 17:40:16
@modify date 2020-09-15 21:46:10
@desc        RIBOSS pipline
"""

import os
import re
import sys
import time
import argparse
import textwrap
import csv
import pandas as pd
import scipy.stats as stats


def check_arg(args=None):
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
    description=textwrap.dedent('''\
      RIBOSS:
        Comparing the strength of open reading frames within individual transcripts,
        i.e. the degree of triplet periodicity using Fisher-exact test
        '''),
    epilog=textwrap.dedent('''\
      example:
        python riboss.py -i ribomap.base -r gencode.v35.pc_transcripts_filtered.pkl.gz
        
      if you find this useful, please consider citing:
        Lim, C.S. & Brown, C.M. (2020) Synergistic effects of upstream open reading frames and 
            leader exon-exon junctions on protein expression. bioRxiv.
        Lim, C.S., Wardell, S.J.T., Kleffmann, T. & Brown, C.M. (2018) The exon-intron gene 
            structure upstream of the initiation codon predicts translation efficiency.
            Nucleic Acids Res. 46: 4575-4591.
        '''))
    parser.add_argument('-i', '--input',
                        metavar='STR',
                        help='.base from Ribomap',
                        required='True') 
    parser.add_argument('-r', '--ref',
                        metavar='STR',
                            help='.pkl from prep_input_files.py',
                        required='True')
    parser.add_argument('-o', '--output_path',
                        default='./',
                        metavar='STR',
                        help='output directory. Default = ./')   

    results = parser.parse_args(args)
    return (results.input, results.ref, results.output_path)


def io(profile, outdir):
    check_dir = os.path.isdir(outdir)
    if not check_dir:
        os.makedirs(outdir)

    try:
        base = os.path.basename(profile)
        if re.search('\.base$', profile):
            phase = outdir + '/' + os.path.splitext(base)[0] + '_periodicity.txt'
            stat = outdir + '/' + os.path.splitext(base)[0] + '_fisher_test.txt'
    except Exception:
        raise argparse.ArgumentTypeError('Please provide a Ribomap .base file as input!')
    return phase, stat
    
    
def all_orf(seq):
    triplet = ""
    for i in range(0, len(seq)-2, 3):
        triplet = seq[i:(i+3)]
        if triplet in ["TAG", "TAA", "TGA"]:
            return seq[:i+3]
    return seq

def top_3_frames(seq):
    positions = []
    orf = ""
    for i in range(0, len(seq)):
        triplet = seq[i:i+3]
        if(triplet == "ATG"):
            orf = all_orf(seq[i:])
            positions.append([i, i + len(orf)])
    return positions

def orf_finder(data):
    """
    Input pickle files from prep_input_files.py
    """
    print('Finding ORFs in transcript isoforms ...')
    start_time = time.perf_counter()
    df = pd.read_pickle(data)
    morf = df[['tid','CDS_range']].copy()
    morf['ORF_type'] = 'mORF'
    morf.columns = ['tid','ORF_range','ORF_type']
    
    df['ORF_range'] = df.seq.apply(lambda x: top_3_frames(x))
    df = df.explode('ORF_range')
    df = df.dropna(axis=0).reset_index(drop=True)
    df = pd.concat([df, pd.DataFrame(df.ORF_range.values.tolist(), 
                                   columns = ['ORF_start','ORF_end'])], 
                  axis=1)

    uorf = df[(df.ORF_start<df.CDS_start) & 
              (df.ORF_end<=df.CDS_start)][['tid','ORF_range']].copy()
    uorf['ORF_type'] = 'uORF'
    dorf = df[df.ORF_start>=df.CDS_end][['tid','ORF_range']].copy()
    dorf['ORF_type'] = 'dORF'
    df = pd.concat([uorf,morf,dorf])
    print('Finished ORF finding in', 
          round((time.perf_counter()-start_time)/60), 'min',
          (time.perf_counter()-start_time)%60, 's\n')
    return df


def parse_ribomap(base):
    print('Parsing ribomap output ...')
    start_time = time.perf_counter()
    dt = pd.read_csv(base, sep=':\s', engine='python', header=None)
    dt = dt[(dt[0]=='tid') | (dt[0]=='ribo profile')]
    dt.columns = ['id','val']
    dt = pd.DataFrame({'tid':dt['val'].iloc[::2].values, 'rprofile':dt['val'].iloc[1::2].values})
    dt['tid'] = dt.tid.str.split('|').str[0]
    dt['rprofile'] = dt.rprofile.str.split().apply(lambda x: [float(i) for i in x])
    print('Finished parsing ribomap output in', 
          round((time.perf_counter()-start_time)/60), 'min',
          (time.perf_counter()-start_time)%60, 's\n')
    return dt


def footprint_counts(df, d):
    """
    Input from orf_finder and parse_ribomap
    """
    print('Counting footprints by reading frames ...')
    dt = pd.merge(df, d, on='tid')
    dt['tid'] = dt.tid.str.split('|').str[0] #.str.split('.').str[0]
    start_time = time.perf_counter()
    dt['frame1'] = dt.ORF_range.apply(lambda x: [i for i in range(x[0], x[1], 3)])
    dt['frame2'] = dt.ORF_range.apply(lambda x: [i for i in range(x[0]+1, x[1], 3)])
    dt['frame3'] = dt.ORF_range.apply(lambda x: [i for i in range(x[0]+2, x[1], 3)])
    dt['periodicity'] = dt[['rprofile','frame1','frame2','frame3']].values.tolist()
    dt['periodicity'] = dt.periodicity.apply(lambda x: [sum([x[0][i] for i in tuple(x[1])]),
                                                            sum([x[0][i] for i in tuple(x[2])]),
                                                            sum([x[0][i] for i in tuple(x[3])])])
    dt = pd.concat([dt, pd.DataFrame(dt.periodicity.values.tolist(), 
                                     columns = ['Frame1_ribosome','Frame2_ribosome','Frame3_ribosome'])], 
                   axis=1).drop(['rprofile','frame1','frame2','frame3','periodicity'], axis=1)
    dt = dt[dt.Frame1_ribosome + dt.Frame2_ribosome + dt.Frame3_ribosome>10]
    print('Finished footprint counting in', 
          round((time.perf_counter()-start_time)/60), 'min',
          (time.perf_counter()-start_time)%60, 's\n')
    return dt


def fisher_test(dt, phase, stat):
    """
    Input from footprint_counts
    """
    print('Comparing the triplet periodicity of mORFs with uORFs and dORFs ...')
    dt['Frame2_3_ribosome'] = dt.Frame2_ribosome + dt.Frame3_ribosome
    m = dt[dt.ORF_type=='mORF']
    o = dt[(dt.ORF_type=='uORF') | (dt.ORF_type=='dORF')]
    oo = pd.merge(o,o,on='tid')
    o1 = oo[(oo.ORF_range_x != oo.ORF_range_y) & 
            (oo.Frame1_ribosome_x/oo.Frame2_3_ribosome_x > 
             oo.Frame1_ribosome_y/oo.Frame2_3_ribosome_y)].iloc[:,0:7]
    o2 = oo[(oo.ORF_range_x != oo.ORF_range_y) & 
            (oo.Frame1_ribosome_x/oo.Frame2_3_ribosome_x == 
             oo.Frame1_ribosome_y/oo.Frame2_3_ribosome_y)].iloc[:,0:7]

    o1['oid'] = o1.tid + o1.ORF_range_x.astype('str')
    o = pd.concat([o1,o2]).drop_duplicates(subset=['oid'])
    o.drop('oid', inplace=True, axis=1)
    o.columns = ['tid','ORF_range','ORF_type','Frame1_ribosome','Frame2_ribosome','Frame3_ribosome','Frame2_3_ribosome']
    oo = pd.concat([m,o]).drop('Frame2_3_ribosome',axis=1)

    f = pd.merge(o, m, on='tid')
    f['tab_x'] = f[['Frame1_ribosome_x','Frame2_3_ribosome_x']].values.tolist()
    f['tab_y'] = f[['Frame1_ribosome_y','Frame2_3_ribosome_y']].values.tolist()
    f['tab'] = f[['tab_x','tab_y']].values.tolist()
    f['Fisher_exact'] = f.tab.apply(lambda x: list(stats.fisher_exact(x, alternative='greater')))
    f['Adjusted_pvalue'] = f.Fisher_exact.apply(lambda x : 1 if x[1]*len(f) > 1 else (x[1]*len(f) if x[1]*len(f) < 0.05 else x[1]))
    f = f[['tid','ORF_range_x','ORF_type_x','tab','Fisher_exact','Adjusted_pvalue']]
    f.columns = ['tid','ORF_range','ORF_type','Contingency_table','Fisher_exact','Adjusted_pvalue']

    oo.to_csv(phase, sep='\t', index=None, quoting = csv.QUOTE_NONE, escapechar = ' ')
    f.to_csv(stat, sep='\t', index=None, quoting = csv.QUOTE_NONE, escapechar = ' ')


def main():
    start_time = time.perf_counter()
    phase,stat = io(i, o)
    df = parse_ribomap(i)
    orf = orf_finder(r)
    df = footprint_counts(orf, df)
    fisher_test(df, phase, stat)
    print('Finished all analysis in', 
          round((time.perf_counter()-start_time)/60), 'min',
          (time.perf_counter()-start_time)%60, 's \
          \n\nOutput files:\n', phase, '\n', stat)

if __name__ == "__main__":
    i, r, o = check_arg(sys.argv[1:])
    main()

