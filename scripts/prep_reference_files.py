#!/usr/bin/env python
# coding: utf-8
"""
@author      CS Lim
@create date 2020-09-15 17:40:16
@modify date 2020-09-15 20:19:18
@desc        RIBOSS pipline
"""

import os
import sys
import re
import time
import argparse
import textwrap
import pandas as pd
import csv
import gzip
from mimetypes import guess_type
from functools import partial
from Bio import SeqIO


def check_arg(args=None):
    """
    Supported databases and example input files
    database, sequence, cds
    Ensembl, ftp://ftp.ensemblgenomes.org/pub/plants/release-48/fasta/arabidopsis_thaliana/cdna/Arabidopsis_thaliana.TAIR10.cdna.all.fa.gz, ftp://ftp.ensemblgenomes.org/pub/plants/release-48/fasta/arabidopsis_thaliana/cds/Arabidopsis_thaliana.TAIR10.cds.all.fa.gz
    FlyBase, ftp://ftp.flybase.net/releases/FB2020_04/dmel_r6.35/fasta/dmel-all-transcript-r6.35.fasta.gz, ftp://ftp.flybase.net/releases/FB2020_04/dmel_r6.35/fasta/dmel-all-CDS-r6.35.fasta.gz
    GENCODE, ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_35/gencode.v35.pc_transcripts.fa.gz, not_required
    RefSeq, ftp://ftp.ncbi.nlm.textwrapnih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_rna.gbff.gz, not_required
    """

    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, 
    description=textwrap.dedent('''\
      Preparing reference files of mRNA isoforms for RIBOSS pipeline.
        Only CDS divisible by three will be retained. Frameshifting genes will not be analysed
        '''), 
    epilog=textwrap.dedent('''\
      examples:
        python prep_input_files.py \\
            -d ensembl \\
            -s Arabidopsis_thaliana.TAIR10.cdna.all.fa.gz \\
            -c Arabidopsis_thaliana.TAIR10.cds.all.fa.gz
        
        python prep_input_files.py \\
            -d refseq \\
            -s GRCh38_latest_rna.gbff.gz \\
            -o ref

      if you find this useful, please consider citing:
        Lim, C.S. & Brown, C.M. (2020) Synergistic effects of upstream open reading frames and 
            leader exon-exon junctions on protein expression. bioRxiv.
        Lim, C.S., Wardell, S.J.T., Kleffmann, T. & Brown, C.M. (2018) The exon-intron gene 
            structure upstream of the initiation codon predicts translation efficiency.
            Nucleic Acids Res. 46: 4575-4591.
        '''))

    parser.add_argument('-d', '--database',
                        type=check_database,
                        metavar='STR',
                        help='specify either {ensembl, flybase, gencode, refseq}',
                        required='True')
    parser.add_argument('-s', '--sequence',
                        metavar='STR',
                        help='input reference transcript sequences in \
                        fasta or RefSeq gbff format',
                        required='True')
    parser.add_argument('-c', '--cds',
                        metavar='STR',
                        help='input reference CDS sequences in \
                        fasta format. Required for ensembl and flybase')
    parser.add_argument('-o', '--output_path',
                        default='./',
                        metavar='STR',
                        help='output directory. Automatically create a new directory if it does not exist. Default = ./')  
    
    results = parser.parse_args(args)
    return (results.database, results.sequence, results.cds, results.output_path)

    
def check_database(db_name):
    if db_name in ['refseq', 'ensembl', 'gencode', 'flybase']:
        return db_name
    else:
        raise argparse.ArgumentTypeError('Database not supported!')


def io(seq, outdir):
    check_dir = os.path.isdir(outdir)
    if not check_dir:
        os.makedirs(outdir)
    try:
        base = os.path.basename(seq) # .splitext(bamfile)
        if re.search('\.fasta\.gz$|\.fna\.gz$|\.fa\.gz$|\.gbff\.gz$', seq):
            base = os.path.splitext(base)[0]
            fname = outdir + '/' + os.path.splitext(base)[0] + '_filtered.fasta'
            pkl = outdir + '/' + os.path.splitext(base)[0] + '_filtered.pkl.gz'
            txt = outdir + '/' + os.path.splitext(base)[0] + '_filtered.txt'
        elif re.search('\.fasta$|\.fna$|\.fa$|\.gbff$', seq):
            fname = outdir + '/' + base + '_filtered.fasta'
            pkl = outdir + '/' + base + '_filtered.pkl.gz'
            txt = outdir + '/' + base + '_filtered.txt'
    except Exception:
        raise argparse.ArgumentTypeError('Please provide a transcript fasta file as input!')
    return pkl, txt, fname


def output_files(data, pickle, tsv, outseq):
    """
    Writing output as fasta, pickle, and txt (CDS range)
    """
    data.to_pickle(pickle)
    data[['tid','CDS_start','CDS_end']].to_csv(tsv, sep='\t', index=None, header=None,
                                              quoting = csv.QUOTE_NONE, escapechar = ' ')
    ('>' + data.tid + '\n' + data.seq).to_csv(outseq, sep='\t', index=None, header=None,
                                              quoting = csv.QUOTE_NONE, escapechar = ' ')

    
def fasta_to_dataframe(seq):
    fasta_df = pd.read_csv(seq, sep='>', lineterminator='>', header=None)
    fasta_df[['tid','seq']]=fasta_df[0].str.split('\n', 1, expand=True)
    fasta_df['tid'] = fasta_df['tid']
    fasta_df['seq'] = fasta_df['seq'].replace('\n','', regex=True)
    fasta_df.drop(0, axis=1, inplace=True)
    fasta_df = fasta_df[fasta_df.seq != '']
    final_df = fasta_df.dropna()
    return final_df


def ensembl_flybase_parser(tx_seq, cds_seq, database):
    tx = fasta_to_dataframe(tx_seq)
    tx['tid'] = tx.tid.str.split().apply(lambda x: x[0])
    cds = fasta_to_dataframe(cds_seq)
    
    if database == 'ensembl':
        cds['tid'] = cds.tid.str.split().apply(lambda x: x[0])
    elif database == 'flybase':
        cds['tid'] = cds.tid.str.split().\
        apply(lambda x: [i for i in x if re.search('FBtr',i)]).\
        apply(lambda x: re.sub(';','',x[0].split(',')[1]))
    else:
        print('\nFailed! Please check your input files and arguments.\n')

    df = pd.merge(cds, tx,on='tid')
    df['seq'] = df[['seq_x','seq_y']].values.tolist()
    df = df[df.seq.apply(lambda x: x[0] in x[1])]
    df['CDS_range'] = df.seq.apply(lambda x: [x[1].index(x[0]), x[1].index(x[0]) + len(x[0])])
    df.drop(['seq_x','seq'], inplace=True, axis=1)
    df.columns = ['tid','seq','CDS_range']
    df.reset_index(drop=True, inplace=True)
    df = pd.concat([df, 
                    pd.DataFrame(df.CDS_range.tolist(), columns = ['CDS_start','CDS_end'])], 
                   axis=1)
    df = df[(df.CDS_end-df.CDS_start)%3==0].reset_index(drop=True)
    return df
    

def gencode_parser(seq):
    df = fasta_to_dataframe(seq)
    df['CDS'] = df.tid.str.split('|').apply(lambda x: [i for i in x if i.startswith('CDS:')]).apply(pd.Series)
    df['tid'] = df.tid.str.split('|').apply(lambda x: x[0])
    df[['CDS_start','CDS_end']] = df.CDS.str.replace('CDS:','').str.split('-', expand=True)
    df.dropna(inplace=True)
    df['CDS_start'] = df.CDS_start.astype('int')
    df['CDS_start'] = df.CDS_start - 1
    df['CDS_end'] = df.CDS_end.astype('int')
    df = df[(df.CDS_end-df.CDS_start)%3==0]
    df['CDS_range'] = df[['CDS_start','CDS_end']].values.tolist()
    df = df.drop('CDS', axis=1).reset_index(drop=True)
    return df


def genbank_parser(f):
    encoding = guess_type(f)[1]  # uses file extension
    _open = partial(gzip.open, mode='rt') if encoding == 'gzip' else open

    seq = []
    with _open(f) as input_handle:
        for record in SeqIO.parse(input_handle, 'genbank'):
            if bool(re.search('NM_', record.id)) is True:
                for feat in record.features:
                    if feat.type == 'CDS':
                        try:
                            seq.append([record.id, 
                                        feat.location._start.position, 
                                        feat.location._end.position, 
                                        str(record.seq)])
                        except AttributeError:
                            print('CompoundLocation found in', record.id, '!')

    df = pd.DataFrame(seq)
    df.columns = ['tid','CDS_start','CDS_end','seq']
    df = df[(df.CDS_end-df.CDS_start)%3==0].reset_index(drop=True)
    df['CDS_range'] = df[['CDS_start','CDS_end']].values.tolist()
    return df


def main():
    start_time = time.perf_counter()
    print('Parsing', s, '...')

    if io(s, o):
        if d in ['gencode','refseq']:
            if c is not None:
                print('Warning: CDS file not needed!')
            pkl,txt,fname = io(s, o)
            if d=='gencode':
                df = gencode_parser(s)
            else:
                df = genbank_parser(s)
            output_files(df, pkl, txt, fname)                
            print('Finished in', 
                  round((time.perf_counter()-start_time)/60), 'min',
                  (time.perf_counter()-start_time)%60, 's \
                  \n\nOutput files:\n', fname, '\n', txt, '\n', pkl)

        elif c is not None:
            pkl,txt,fname = io(s, o)
            df = ensembl_flybase_parser(s, c, d)
            output_files(df, pkl, txt, fname)
            print('Finished in', 
                  round((time.perf_counter()-start_time)/60), 'min',
                  (time.perf_counter()-start_time)%60, 's. \
                  \n\nOutput files:\n', fname, '\n', txt, '\n', pkl)
        else:
            print('Failed! Please check your input files and arguments.')
    else:
        print('Failed! Please check your input files and arguments.')

        
if __name__ == "__main__":
    d, s, c, o = check_arg(sys.argv[1:])
    if d in ['ensembl', 'flybase']:
        if c is None:
            raise ValueError('Please input the CDS file!')
    main()

