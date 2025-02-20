#!/usr/bin/env python
# coding: utf-8
 
"""
@author      CS Lim
@create date 2024-09-13 15:26:12
@modify date 2025-02-19 14:22:57
@desc        RIBOSS module for finding ORFs
"""




import os, re, time, argparse, csv, gzip, logging.config, subprocess
from io import StringIO
import numpy as np
import pandas as pd 
import seaborn as sns
import matplotlib.pyplot as plt
from mimetypes import guess_type
from functools import partial
from Bio import SeqIO
from .wrapper import filename
import pyranges as pr
from tqdm import tqdm
from .wrapper import fasta_to_dataframe


DEFAULT_LOGGING = {
    'version': 1,
    'disable_existing_loggers': False,
    'loggers': {
        '': {
            'level': 'INFO',
        },
        'another.module': {
            'level': 'ERROR',
        },
    }
}

logging.config.dictConfig(DEFAULT_LOGGING)



CODON_TO_AA={'TTT':'F','TCT':'S','TAT':'Y','TGT':'C','TTC':'F','TCC':'S',\
             'TAC':'Y','TGC':'C','TTA':'L','TCA':'S','TAA':'*',\
             'TGA':'*','TTG':'L','TCG':'S','TAG':'*','TGG':'W',\
             'CTT':'L','CCT':'P','CAT':'H','CGT':'R','CTC':'L','CCC':'P',\
             'CAC':'H','CGC':'R','CTA':'L','CCA':'P','CAA':'Q','CGA':'R',\
             'CTG':'L','CCG':'P','CAG':'Q','CGG':'R','ATT':'I','ACT':'T',\
             'AAT':'N','AGT':'S','ATC':'I','ACC':'T','AAC':'N','AGC':'S',\
             'ATA':'I','ACA':'T','AAA':'K','AGA':'R','ATG':'M','ACG':'T',\
             'AAG':'K','AGG':'R','GTT':'V','GCT':'A','GAT':'D','GGT':'G',\
             'GTC':'V','GCC':'A','GAC':'D','GGC':'G','GTA':'V','GCA':'A',\
             'GAA':'E','GGA':'G','GTG':'V','GCG':'A','GAG':'E','GGG':'G'}



def translate(seq):
    seq = seq.upper()
    if seq[-3:] in ['TAA','TAG','TGA']:
        seq = seq[:-3]
        
    length = (len(seq)- len(seq)%3)
    split_func = lambda seq, n: [seq[i:i+n] for\
                                    i in range(0, length, n)]
    codons = split_func(seq, 3)
    aa = ''
    for c in codons:
        aa+=CODON_TO_AA[c]
    return aa

    


def all_orf(seq, stop_codon=False):
    seq = seq.upper()
    triplet = ""
    for i in range(0, len(seq)-2, 3):
        triplet = seq[i:(i+3)]
        if triplet in ['TAA', 'TAG', 'TGA']:
            return seq[:i+3]
    if stop_codon==False:
        return seq



def top_3_frames(seq, start_codon, stop_codon=False):
    seq = seq.upper()
    positions = []
    orf = ""
    for i in range(0, len(seq)):
        triplet = seq[i:i+3]
        if type(start_codon)==str:
            if (triplet==start_codon.upper()):
                orf = all_orf(seq[i:], stop_codon=stop_codon)
                if orf is not None:
                    l = len(orf)
                    l -= l % +3
                    positions.append([start_codon, i, i + l])
        else:
            for codon in start_codon:
                if (triplet==codon.upper()):
                    orf = all_orf(seq[i:], stop_codon=stop_codon)
                    if orf is not None:
                        l = len(orf)
                        l -= l % +3
                        positions.append([codon, i, i + l])
    return positions



def orf_finder(annotation, tx, ncrna=False, outdir=None, start_codon=["ATG", "CTG", "GTG", "TTG"], stop_codon=False):
    """
    Input:
        * annotation: gene annotation file in GTF, GFF3 or BED format 
        * tx: transcript fasta file (required)
        * ncrna: remove noncoding RNAs from analysis
        * outdir: output directory (default: None)
        * start_codon: any triplets. If a list is given, it should be sorted from high to low abundance in the species of interest (default: ATG, GTG, TTG, CTG)
    Output:
        * cds_range: as input for analyse_footprints
        * orf: as input for riboss
        * df: as input for footprint_counts
    """
    
    start_time = time.perf_counter()

    fname = filename(annotation, None, outdir)
    
    path, ext = os.path.splitext(annotation)
    genepred = path + '.gp'
    bed = path + '.bed'
    
    if 'gtf' in annotation:
        subprocess.run(['gtfToGenePred', annotation, genepred], check=True)
    elif ('gff' in annotation) | ('gff3' in annotation):
        subprocess.run(['gff3ToGenePred', annotation, genepred], check=True)
    elif 'bed' in annotation:
        subprocess.run(['bedToGenePred', annotation, genepred], check=True)
    elif 'gp' in annotation:
        pass

    subprocess.run(['genePredToBed', genepred, bed], check=True)
    gp = pd.read_csv(genepred, sep='\t', header=None)
    
    if ncrna==False:
        gp = gp[gp[5]!=gp[6]]
        
    gp[8] = gp[8].str.split(',').str[:-1].apply(lambda x: [int(i) for i in x])
    gp[9] = gp[9].str.split(',').str[:-1].apply(lambda x: [int(i) for i in x])
    gp['exons'] = gp[[8,9]].values.tolist()
    gp['exons'] = gp.exons.apply(lambda x: list(zip(x[0],x[1])))
    
    # # get transcript sequence
    # exons = gp.explode('exons')
    # exons['Start'] = exons.exons.apply(lambda x: x[0])
    # exons['End'] = exons.exons.apply(lambda x: x[1])
    # exons.rename(columns = {1:'Chromosome',2:'Strand',0:'Name'}, inplace=True)
    # seq = pr.get_transcript_sequence(pr.PyRanges(exons), group_by='Name', path=fasta)    
    
    # get individual positions as lists
    gp['exons'] = gp.exons.apply(lambda x: [list(range(i[0],i[1])) for i in x])
    gp['exons'] = gp.exons.apply(lambda x: [j for i in x for j in i])
    
    df = gp[[1,3,4,0,2,5,6,'exons']]
    df.columns = ['Chromosome','Start','End','Name','Strand','start','end','exons']

    seq = fasta_to_dataframe(tx)
    seq.columns = ['Name','Sequence']

    df = pd.merge(df, seq)
    if df.shape[0]==0:
        logging.error('Error while merging annotation and fasta sequence header! Make sure that the name column for annotation ID.')

    pbar = tqdm.pandas(desc="finding all ORFs       ", unit_scale=True, ncols=100)
    df['ORF_range'] = df.Sequence.progress_apply(lambda x: top_3_frames(x, start_codon, stop_codon=stop_codon))
    df = df.explode('ORF_range')
    df = df.dropna(axis=0).reset_index(drop=True)
    df['start_codon'] = df['ORF_range'].apply(lambda x: x[0])
    df['ORF_range'] = df['ORF_range'].apply(lambda x: [x[1],x[2]])
    df['ORF_start'] = df.ORF_range.apply(lambda x: x[0])
    df['ORF_end'] = df.ORF_range.apply(lambda x: x[1])
    df['ORF_length'] = df.ORF_end - df.ORF_start
    
    pbar = tqdm.pandas(desc="converting coordinates ", unit_scale=True, ncols=100)
    df_plus = df[df.Strand=='+'].copy()
    df_plus['genomic_coordinates'] = df_plus[['exons','ORF_range']].values.tolist()
    df_plus['genomic_range'] = df_plus.genomic_coordinates.progress_apply(lambda x: x[0][x[1][0]:x[1][1]])
    df_plus['genomic_range'] = df_plus.genomic_range.apply(lambda x: [x[0],x[-1]+1])
    
    df_minus = df[df.Strand=='-'].copy()
    df_minus['genomic_coordinates'] = df_minus.exons.progress_apply(lambda x: list(reversed(x)))
    df_minus['genomic_coordinates'] = df_minus[['genomic_coordinates','ORF_range']].values.tolist()
    df_minus['genomic_range'] = df_minus.genomic_coordinates.apply(lambda x: x[0][x[1][0]:x[1][1]])
    df_minus['genomic_range'] = df_minus.genomic_range.apply(lambda x: [x[-1],x[0]+1])
    
    df = pd.concat([df_plus,df_minus]).drop(['exons','genomic_coordinates'], axis=1)
    df['genomic_start'] = df.genomic_range.apply(lambda x: x[0])
    df['genomic_end'] = df.genomic_range.apply(lambda x: x[1])
    
    cds = df[(df.genomic_start==df.start) & (df.genomic_end==df.end)].copy()
    # Export CDS_range
    cds_range = cds
    cds_range['fasta_header'] = '>' + cds_range['Name']
    cds_range[['fasta_header','Sequence']].to_csv(fname + '.transcripts.fa', index=None, header=None, sep='\n')
    cds_range = cds_range[['Name','ORF_start','ORF_end']]
    cds_range.columns = ['tid','CDS_start','CDS_end']
    cds_range.to_csv(fname + '.cds_range.txt', sep='\t', index=None, header=None)
    cds = cds.drop(['Start','End','start','end','genomic_range','Sequence'], axis=1).copy()
    cds.rename(columns={'genomic_start':'Start','genomic_end':'End'}, inplace=True)
    cds['ORF_type'] = 'mORF'
    
    # find uORFs
    uorf = df[((df.genomic_start<df.start) & (df.Strand=='+')) | ((df.genomic_end>df.end) & (df.Strand=='-'))]
    uorf = uorf.drop(['Start','End','start','end','genomic_range','Sequence'], axis=1).copy()
    uorf.rename(columns={'genomic_start':'Start','genomic_end':'End'}, inplace=True)
    uorf['ORF_type'] = 'uORF'
    
    # find dORFs
    dorf = df[((df.genomic_start>=df.end) & (df.Strand=='+')) | ((df.genomic_end<=df.start) & (df.Strand=='-'))]
    dorf = dorf.drop(['Start','End','start','end','genomic_range','Sequence'], axis=1).copy()
    dorf.rename(columns={'genomic_start':'Start','genomic_end':'End'}, inplace=True)
    dorf['ORF_type'] = 'dORF'
    orf = pd.concat([uorf,dorf])
    
    # find ORFs in-framed with mORFs
    cdsorf = pr.PyRanges(cds[['Chromosome','Start','End','Strand']]).join(pr.PyRanges(orf), strandedness='same').df
    cdsorf['frame_plus'] = (cdsorf.Start-cdsorf.Start_b)%3
    cdsorf['frame_minus'] = (cdsorf.End-cdsorf.End_b)%3
    inframe = cdsorf[(cdsorf.frame_plus==0) | (cdsorf.frame_minus==0)]
    
    # find overlapping ORFs
    oorf = pd.concat([cdsorf, inframe, cds])
    oorf = oorf.drop_duplicates(['Name','start_codon','ORF_start','ORF_end'], keep=False)
    oorf = oorf.drop(['Start', 'End'], axis=1).rename(columns={'Start_b':'Start', 'End_b':'End'})
    oorf['ORF_type'] = 'oORF'
    
    df = pd.concat([uorf, oorf, dorf]).drop_duplicates(['Name','start_codon','ORF_start','ORF_end'])
    df = pd.concat([df, inframe]).drop_duplicates(['Name','start_codon','ORF_start','ORF_end'], keep=False)
    
    if type(start_codon)==list:
        ds = []
        for i in start_codon:
            ds.append(df[df.start_codon==i].sort_values('ORF_length', ascending=False))
        df = pd.concat(ds)
    else:
        df = df.sort_values('ORF_length').drop_duplicates(['Name','ORF_end'])
        
    df.drop_duplicates(['Name','ORF_end'], inplace=True)
    df = pd.concat([cds, df]).drop_duplicates(['Name','start_codon','ORF_start','ORF_end'])
    df = df.drop(['fasta_header','Strand_b','frame_plus','frame_minus','Start_b','End_b'], axis=1).rename(columns={'Name':'tid'})
    df['ORF_length'] = df.ORF_range.apply(lambda x: x[1]-x[0])

    # # remove in-frame ORFs again
    # morfs = pr.PyRanges(df[df.ORF_type=='mORF'])
    # norfs = pr.PyRanges(df[df.ORF_type!='mORF'])
    
    # cdsorf = norfs.join(morfs).df
    # cdsorf['frame_plus'] = (cdsorf.Start-cdsorf.Start_b)%3
    # cdsorf['frame_minus'] = (cdsorf.End-cdsorf.End_b)%3
    # inframe = cdsorf[(cdsorf.frame_plus==0) | (cdsorf.frame_minus==0)].copy()
    
    # df = pd.concat([df,inframe]).drop_duplicates(['Chromosome', 'tid', 'Strand', 'start_codon', 'ORF_start',
    #        'ORF_end', 'ORF_length', 'Start', 'End', 'ORF_type'], keep=False)
    # df = df[['Chromosome', 'tid', 'Strand', 'ORF_range', 'start_codon', 'ORF_start',
    #        'ORF_end', 'ORF_length', 'Start', 'End', 'ORF_type']].reset_index(drop=True)

    # Get actual oORFs. Remove those located within introns
    oorfs = df[df.ORF_type=='oORF'].copy()
    oorfs['Start'] = oorfs.Start.astype(int)
    oorfs['End'] = oorfs.End.astype(int)
    oorfs['oid'] = oorfs.tid + '__' + oorfs.ORF_range.str[0].astype(int).astype(str) + '-' + oorfs.ORF_range.str[1].astype(int).astype(str)
    oorfs[['Chromosome','Start','End','oid','ORF_length','Strand']].to_csv('oorfs.bed', sep='\t', header=None, index=None)
    # BED12 without UTRs
    morfs = pd.read_csv(bed, sep='\t', header=None)
    morfs[[0,6,7,3,4,5,6,7,8,9,10,11]].to_csv('morfs.bed', sep='\t', header=None, index=None)
    # actual oORFs
    results = subprocess.run(['bedtools','intersect','-a','oorfs.bed','-b','morfs.bed','-s','-split','-wo'],
                   stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)

    oorf_ids = pd.read_csv(StringIO(results.stdout), sep='\t', header=None)
    if oorf_ids.shape[0]>0:
        oorf_ids = oorf_ids[[3]]
        
    oorf_ids['tid'] = oorf_ids[3].str.split('__').str[0]
    oorf_ids['ORF_start'] = oorf_ids[3].str.split('__').str[1].str.split('-').str[0].astype(int)
    oorf_ids['ORF_end'] = oorf_ids[3].str.split('__').str[1].str.split('-').str[1].astype(int)
    oorf_ids = oorf_ids[['tid','ORF_start','ORF_end']]
    oorfs = pd.merge(oorfs,oorf_ids).drop('oid', axis=1)
    df = pd.concat([df[df.ORF_type!='oORF'],oorfs]).reset_index(drop=True)
    
    df.to_pickle(fname + '.orf_finder.pkl.gz')
    os.remove('morfs.bed')
    os.remove('oorfs.bed')
    
    logging.info('found ' + str(df.shape[0]) + ' ORFs in ' + 
          str(round((time.perf_counter()-start_time)/60)) + ' min ' +
          str(round((time.perf_counter()-start_time)%60)) + ' s')
    logging.info('saved sequences as ' + fname + '.transcripts.fa')
    logging.info('saved sequences as ' + fname + '.gp')
    logging.info('saved sequences as ' + fname + '.orf_finder.pkl.gz')
    logging.info('saved CDS range as ' + fname + '.cds_range.txt')

    return cds_range, df



def operon_distribution(op, displot_prefix, log=False):
    """
    Plot the distribution of operons.

    Input:
        * op: dataframe for operon descriptive statistics
        * displot_prefix: filename prefix for displot and jointgrid
    
    Output:
        * displot and jointgrid as PDFs
    """
    
    sns.displot(op['count'], height=3)
    plt.title('Operons predicted from transcriptome assembly')
    plt.xlabel('Number of ORFs per mRNA')
    plt.savefig(displot_prefix + '.operon_dist.pdf', bbox_inches='tight')

    plt.figure(figsize=(3,3))
    sns.boxplot(data=op, x='count', y='length')
    if log==True:
        plt.yscale('log')
    plt.xlabel('Number of ORFs per mRNA')
    plt.ylabel('mRNA length')
    plt.savefig(displot_prefix + '.operon_boxplot.pdf', bbox_inches='tight')
    
    logging.info('plotted the distribution of operons as ' + displot_prefix + '.operon_dist.pdf and ' + displot_prefix + '.operon_scatter.pdf')
    


def operon_finder(tx_assembly, bed, outdir=None, delim=None, 
                  start_codon=["ATG", "GTG", "TTG", "CTG"], stop_codon=False,
                  ncrna=False,log=False):
    """
    Predict operons from transcriptome.
    
    Input:
        * tx_assembly: transcript fasta file extracted using bedtools getfasta. Headers with genomic coordinates (required)
        * bed: converted from gff3, e.g. from NCBI Genome Assembly (required)
        * outdir: output directory (default: None)
        * delim: use :: for tx_assembly extracted using bedtools getfasta -name flag, as this appends name (column #4) to the genomic coordinates (default: None)
        * start_codon: any triplets. If a list is given, it should be sorted from high to low abundance in the species of interest (default: ATG, GTG, TTG, CTG)
        * ncrna: remove noncoding RNAs from analysis
    Output:
        * cds_range: as input for analyse_footprints
        * orf: as input for riboss
        * df: as input for footprint_counts
    """
 
    if delim!=None:
        pos = 1
    else:
        pos = 0
        
    df = fasta_to_dataframe(tx_assembly)
    df['Chromosome'] = df.tid.str.split(':').apply(lambda x: x[pos+1])
    df['Start'] = df.tid.str.split(':').apply(lambda x: x[pos+2]).str.split('-').apply(lambda x: x[0]).astype(int)
    df['End'] = df.tid.str.split(':').apply(lambda x: x[pos+2]).str.split('-').apply(lambda x: x[1]).str.split('(').apply(lambda x: x[0]).astype(int)
    df['Strand'] = df.tid.str.split(':').apply(lambda x: x[pos+2]).str.split('(').apply(lambda x: x[1]).str.replace(')','')
    df_ = df
    df_['length'] = df_.seq.apply(len)
    df_ = df_[(df_.Strand=='+') | (df_.Strand=='-')].drop('seq',axis=1)
    
    pbar = tqdm.pandas(desc="finding all ORFs       ", unit_scale=True, ncols=100)
    df['ORF_range'] = df.seq.progress_apply(lambda x: top_3_frames(x ,start_codon, stop_codon=stop_codon))
    df = df.explode('ORF_range')
    df = df.dropna(axis=0).reset_index(drop=True)
    df['start_codon'] = df['ORF_range'].apply(lambda x: x[0])
    df['ORF_range'] = df['ORF_range'].apply(lambda x: [x[1],x[2]])
    df['ORF_start'] = df.ORF_range.apply(lambda x: x[0])
    df['ORF_end'] = df.ORF_range.apply(lambda x: x[1])
    df['ORF_length'] = df.ORF_end - df.ORF_start
    
    orf = df[['Chromosome','Start','End','Strand','tid','start_codon','ORF_start','ORF_end','ORF_length']].drop_duplicates()
    orf_plus = orf[orf.Strand=='+'].copy()
    orf_minus = orf[orf.Strand=='-'].copy()
    orf_plus['Start_b'] = orf_plus.Start+orf_plus.ORF_start
    orf_plus['End_b'] = orf_plus.Start+orf_plus.ORF_end
    orf_minus['Start_b'] = orf_minus.End-orf_minus.ORF_end
    orf_minus['End_b'] = orf_minus.End-orf_minus.ORF_start
    orf = pd.concat([orf_plus, orf_minus])
    orf = orf.drop(['Start','End'], axis=1).rename(columns={'Start_b':'Start','End_b':'End'})
    orf = orf[['Chromosome','Start','End','Strand','tid','start_codon','ORF_start','ORF_end','ORF_length']]

    # Get mORFs encoded by assembled transcripts
    cds = pd.read_csv(bed, sep='\t', header=None)
    cds.columns = ['Chromosome','Start','End','Name','Score','Strand','thickStart','thickEnd','itemRgb','blockCount','blockSizes','blockStarts']
    cds['Strand'] = cds.Strand.astype(str)
    cds_ = cds   

    # Remove ORFs on ncRNAs
    if ncrna==False:
        ncorf = pr.PyRanges(orf).join(pr.PyRanges(cds_), strandedness=False).df
        ncorf = ncorf[ncorf.thickStart==ncorf.thickEnd].copy()
        ncorf = ncorf[['Chromosome','Start','End','Strand','tid','start_codon','ORF_start','ORF_end','ORF_length']]
        orf = pd.concat([orf,ncorf]).drop_duplicates(keep=False)
        orf['Strand'] = orf.Strand.astype(str)
    
    cdsorf = pr.PyRanges(cds).join(pr.PyRanges(orf), strandedness='same').df
    cds = cdsorf[(cdsorf.Start==cdsorf.Start_b) & (cdsorf.End==cdsorf.End_b) & (cdsorf.Strand==cdsorf.Strand_b)]
    cds = cds.drop_duplicates(['Chromosome','Start','End','Strand','Name']) 
    cds['Strand'] = cds.Strand.astype(str)
    cds['ORF_type'] = 'mORF'

    # find ORFs in-framed with mORFs
    cdsorf['frame_plus'] = (cdsorf.Start-cdsorf.Start_b)%3
    cdsorf['frame_minus'] = (cdsorf.End-cdsorf.End_b)%3
    inframe = cdsorf[(cdsorf.frame_plus==0) | (cdsorf.frame_minus==0)]

    # find overlapping ORFs
    oorf = pd.concat([cdsorf, inframe, cds])
    oorf = oorf[['tid','start_codon','ORF_start','ORF_end','ORF_length']].drop_duplicates(keep=False)
    oorf['ORF_type'] = 'oORF'

    # find sORFs
    sorf = pd.concat([orf, cds, inframe, oorf])
    sorf = sorf[['Chromosome','Start','End','Strand','tid','start_codon','ORF_start','ORF_end','ORF_length']]
    sorf.drop_duplicates(['tid','start_codon','ORF_start','ORF_end','ORF_length'], keep=False, inplace=True)
    sorf['Strand'] = sorf.Strand.astype(str)
    sorf['ORF_type'] = 'sORF'

    # find ORFs overlaped with partial mORFs (due to partial transcripts)
    oporf = pr.PyRanges(sorf).join(pr.PyRanges(cds_[['Chromosome','Start','End','Strand']]), strandedness='same').df
    oporf = oporf[(oporf.Strand==oporf.Strand_b)].copy()
    oporf['ORF_type'] = 'opORF'
    oporf = oporf[['tid','start_codon','ORF_start','ORF_end','ORF_length','ORF_type']]

    cds = cds[['tid','start_codon','ORF_start','ORF_end','ORF_length','ORF_type']]
    
    d = pd.concat([oorf, oporf, sorf])
    
    if type(start_codon)==list:
        ds = []
        for i in start_codon:
            ds.append(d[d.start_codon==i].sort_values('ORF_length', ascending=False))
        d = pd.concat(ds)
    else:
        d = d.sort_values('ORF_length', ascending=False)
        
    d.drop_duplicates(['tid','ORF_end'], inplace=True)
    d = pd.concat([cds, d]).drop_duplicates(['tid','start_codon','ORF_start','ORF_end','ORF_length'])
    df = pd.merge(d.drop(['Chromosome','Start','End','Strand'], axis=1), orf)

    # Export CDS_range
    fname = filename(tx_assembly, None, outdir)
    
    cds_range = df[df.ORF_type=='mORF'][['tid','ORF_start','ORF_end']].drop_duplicates()
    cds_range.columns = ['tid','CDS_start','CDS_end']
    cds_range.to_csv(fname + '.cds_range.txt', header=None, index=None, sep='\t')

    #  Plot the distribution of operons
    cdstx_ = pr.PyRanges(cds_).join(pr.PyRanges(df_), strandedness='same')
    op = cdstx_.df.value_counts(['tid','length']).reset_index()
    
    operon_distribution(op, fname, log)

    df.to_pickle(fname + '.operon_finder.pkl.gz')
    
    logging.info('saved operons and ORFs as ' + fname + '.operon_finder.pkl.gz')
    logging.info('saved CDS range as ' + fname + '.cds_range.txt')
    
    return cdstx_, cds_range, df