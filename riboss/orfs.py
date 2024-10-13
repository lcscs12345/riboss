#!/usr/bin/env python
# coding: utf-8
 
"""
@author      CS Lim
@create date 2024-09-13 15:26:12
<<<<<<< HEAD
@modify date 2024-10-13 19:40:02
=======
@modify date 2024-09-14 10:40:54
>>>>>>> refs/remotes/origin/master
@desc        RIBOSS module for finding ORFs
"""




import os, re, time, argparse, csv, gzip, logging
import pandas as pd 
import seaborn as sns
import matplotlib.pyplot as plt
from mimetypes import guess_type
from functools import partial
from Bio import SeqIO
from cgat import GTF
from .converter import TranscriptCoordInterconverter
from .wrapper import filename
import pyranges as pr
from tqdm import tqdm


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
             'TAC':'Y','TGC':'C','TTA':'L','TCA':'S','TAA':'stop',\
             'TGA':'stop','TTG':'L','TCG':'S','TAG':'stop','TGG':'W',\
             'CTT':'L','CCT':'P','CAT':'H','CGT':'R','CTC':'L','CCC':'P',\
             'CAC':'H','CGC':'R','CTA':'L','CCA':'P','CAA':'Q','CGA':'R',\
             'CTG':'L','CCG':'P','CAG':'Q','CGG':'R','ATT':'I','ACT':'T',\
             'AAT':'N','AGT':'S','ATC':'I','ACC':'T','AAC':'N','AGC':'S',\
             'ATA':'I','ACA':'T','AAA':'K','AGA':'R','ATG':'M','ACG':'T',\
             'AAG':'K','AGG':'R','GTT':'V','GCT':'A','GAT':'D','GGT':'G',\
             'GTC':'V','GCC':'A','GAC':'D','GGC':'G','GTA':'V','GCA':'A',\
             'GAA':'E','GGA':'G','GTG':'V','GCG':'A','GAG':'E','GGG':'G'}



def translate(seq):
    seq = seq[:-3]
    length = (len(seq)- len(seq)%3)
    split_func = lambda seq, n: [seq[i:i+n] for\
                                    i in range(0, length, n)]
    codons = split_func(seq, 3)
    aa = ''
    for c in codons:
        aa+=CODON_TO_AA[c]
    return aa



def check_database(db_name):
    if db_name in ['refseq', 'ensembl', 'gencode', 'flybase']:
        return db_name
    else:
        raise argparse.ArgumentTypeError('Database not supported!')


# def io(seq, outdir):
#     check_dir = os.path.isdir(outdir)
#     if not check_dir:
#         os.makedirs(outdir)
#     try:
#         base = os.path.basename(seq)
#         if re.search('.fasta.gz$|.fna.gz$|.fa.gz$|.gbff.gz$', seq):
#             base = os.path.splitext(base)[0]
#             fname = outdir + '/' + os.path.splitext(base)[0] + '_filtered.fasta'
#             pkl = outdir + '/' + os.path.splitext(base)[0] + '_filtered.pkl.gz'
#             txt = outdir + '/' + os.path.splitext(base)[0] + '_filtered.txt'
#         elif re.search('.fasta$|.fna$|.fa$|.gbff$', seq):
#             fname = outdir + '/' + base + '_filtered.fasta'
#             pkl = outdir + '/' + base + '_filtered.pkl.gz'
#             txt = outdir + '/' + base + '_filtered.txt'
#     except Exception:
#         raise argparse.ArgumentTypeError('Please provide a transcript fasta file as input!')
#     return pkl, txt, fname


    
def fasta_to_dataframe(seq):
    fasta_df = pd.read_csv(seq, sep='>', lineterminator='>', header=None)
    df = fasta_df[0].str.split('\n', n=1, expand=True)
    df[1] = df[1].replace('\n','', regex=True)
    df = df[df[1] != '']
    df = df.dropna()
    df.columns = ['tid','seq']
    return df



def ensembl_flybase_parser(tx_seq, cds_seq, db_name):
    tx = fasta_to_dataframe(tx_seq)
    tx['tid'] = tx.tid.str.split().apply(lambda x: x[0])
    cds = fasta_to_dataframe(cds_seq)
    
    if db_name == 'ensembl':
        cds['tid'] = cds.tid.str.split().apply(lambda x: x[0])
    elif db_name == 'flybase':
        cds['tid'] = cds.tid.str.split().\
        apply(lambda x: [i for i in x if re.search('FBtr',i)]).\
        apply(lambda x: re.sub(';','',x[0].split(',')[1]))
    else:
        logging.error('Failed! Please check your input files and arguments.', exc_info=True)

    df = pd.merge(cds, tx,on='tid')
    df['seq'] = df[['seq_x','seq_y']].values.tolist()
    df = df[df.seq.apply(lambda x: x[0] in x[1])]
    df['CDS_range'] = df.seq.apply(lambda x: [x[1].index(x[0]), x[1].index(x[0]) + len(x[0])])
    df.drop(['seq_x','seq'], inplace=True, axis=1)
    df.columns = ['tid','seq','CDS_range']
    df['CDS_start'] = df.CDS_range.apply(lambda x: x[0])
    df['CDS_end'] = df.CDS_range.apply(lambda x: x[1])
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
                            logging.error('CompoundLocation found in', record.id, '!', exc_info=True)

    df = pd.DataFrame(seq)
    df.columns = ['tid','CDS_start','CDS_end','seq']
    df = df[(df.CDS_end-df.CDS_start)%3==0].reset_index(drop=True)
    df['CDS_range'] = df[['CDS_start','CDS_end']].values.tolist()
    return df



def all_orf(seq):
    triplet = ""
    for i in range(0, len(seq)-2, 3):
        triplet = seq[i:(i+3)]
        if triplet in ["TAG", "TAA", "TGA"]:
            return seq[:i+3]
    return seq



def top_3_frames(seq, start_codon):
    positions = []
    orf = ""
    for i in range(0, len(seq)):
        triplet = seq[i:i+3]
        if type(start_codon)==str:
            if (triplet==start_codon):
                orf = all_orf(seq[i:])
                positions.append([start_codon, i, i + len(orf)])
        else:
            for codon in start_codon:
                if (triplet==codon):
                    orf = all_orf(seq[i:])
                    positions.append([codon, i, i + len(orf)])
    return positions



def orf_finder(start_codon, db_name, gtf, seq, cds=None, outdir=None):
    
    """
    Finding ORFs in all transcript isoforms.
    """
    
    start_time = time.perf_counter()
    
    if db_name in ['gencode','refseq']:
        if cds is not None:
            logging.warning('Reminder: CDS file not needed.')
        # pkl,txt,fname = io(seq, output)
        if db_name=='gencode':
            df = gencode_parser(seq)
        else:
            df = genbank_parser(seq)

    elif db_name in ['ensembl','flybase']:
        # pkl,txt,fname = io(seq, output)
        df = ensembl_flybase_parser(seq, cds, db_name)
    else:
        logging.error('Check the database name and files!', exc_info=True)
        
    fname = filename(seq, None, outdir)

    cds_range = df[['tid','CDS_start','CDS_end']]
    cds_range.to_csv(fname + '.cds_range.txt', sep='\t', index=None, header=None)

    df['fasta_header'] = '>' + df.tid
    df[['fasta_header','seq']].to_csv(fname + '.transcripts.fa', index=None, header=None, sep='\n')
    
    df['start_codon'] = df[['seq','CDS_start']].values.tolist()
    df['start_codon'] = df['start_codon'].apply(lambda x: x[0][x[1]:x[1]+3])
    morf = df[['tid','start_codon','CDS_range']].copy()
    morf['ORF_type'] = 'mORF'
    morf.columns = ['tid','start_codon','ORF_range','ORF_type']
    
    pbar = tqdm.pandas(desc="finding all ORFs       ", unit_scale = True)
    df['ORF_range'] = df.seq.progress_apply(lambda x: top_3_frames(x, start_codon))
    df = df.explode('ORF_range')
    df = df.dropna(axis=0).reset_index(drop=True)
    df['start_codon'] = df['ORF_range'].apply(lambda x: x[0])
    df['ORF_range'] = df['ORF_range'].apply(lambda x: [x[1],x[2]])
    df['ORF_start'] = df.ORF_range.apply(lambda x: x[0])
    df['ORF_end'] = df.ORF_range.apply(lambda x: x[1])
    
    uorf = df[(df.ORF_start<df.CDS_start) & 
              (df.ORF_end<=df.CDS_start)][['tid','start_codon','ORF_range']].copy()
    uorf['ORF_type'] = 'uORF'
    dorf = df[df.ORF_start>=df.CDS_end][['tid','start_codon','ORF_range']].copy()
    dorf['ORF_type'] = 'dORF'
    df = pd.concat([uorf,morf,dorf])
    df = orf_resolver(df, gtf)
    # df.to_pickle(pkl)

    logging.info('found ' + str(df.shape[0]) + ' ORFs in ' + 
          str(round((time.perf_counter()-start_time)/60)) + ' min ' +
          str(round((time.perf_counter()-start_time)%60)) + ' s')
    logging.info('saved sequences as ' + fname + '.transcripts.fa')
    logging.info('saved CDS range as ' + fname + '.cds_range.txt')
    
    return cds_range, fname + '.cds_range.txt', df



def orf_resolver(df, gtf):    
    # find overlapping ORFs
    u = df.groupby('tid')['ORF_range'].apply(list).reset_index()
    u['overlap'] = u['ORF_range'].apply(lambda x: [i[1] for i in x])
    u['overlap'] = u['overlap'].apply(lambda x: dict((i,x.count(i)) for i in set(x)))
    u.drop('ORF_range', axis=1, inplace=True)
    df = pd.merge(df,u)
    df['ORF_start'] = df['ORF_range'].apply(lambda x: x[0])
    df['ORF_end'] = df['ORF_range'].apply(lambda x: x[1])
    df['overlap'] = df[['ORF_end','overlap']].values.tolist()
    df['overlap'] = df['overlap'].apply(lambda x: x[1][x[0]])
    df['overlap'] = df['overlap'].apply(lambda x: x if x>1 else 0)
    
    # map transcript to genomic coordinates
    g = pd.read_csv(gtf,sep='\t',comment='#',header=None,low_memory=False)
    g = g[g[2].str.contains('transcript|RNA|pseudogene')]
    g['tid'] = g[8].str.split('; ')
    g = g.explode('tid')
    g = g[g.tid.str.contains('transcript_id')].reset_index(drop=True)
    g['tid'] = g.tid.str.split().apply(lambda x: x[1].replace('"',''))
    g['transcript_iterator'] = [transcript for transcript in GTF.transcript_iterator(GTF.iterator(gzip.open(gtf, mode="rt") if guess_type(gtf)[1] == 'gzip' else open(gtf)))]
    
    pbar = tqdm.pandas(desc="mapping start positions", unit_scale = True)
    t2g_start = pd.merge(df.groupby('tid').start.apply(list).reset_index(),g)
    t2g_start['genomic_start'] = t2g_start[['transcript_iterator','ORF_start']].values.tolist()
    t2g_start['genomic_start'] = t2g_start['genomic_start'].progress_apply(lambda x: [TranscriptCoordInterconverter(x[0]).transcript2genome(i)[0] for i in x[1]])
    pbar = tqdm.pandas(desc="mapping end positions  ", unit_scale = True)
    t2g_end = pd.merge(df.groupby('tid').end.apply(list).reset_index(),g)
    t2g_end['genomic_end'] = t2g_end[['transcript_iterator','ORF_end']].values.tolist()
    t2g_end['genomic_end'] = t2g_end['genomic_end'].progress_apply(lambda x: [TranscriptCoordInterconverter(x[0]).transcript2genome(i-3)[0] for i in x[1]])
    
    t2g = pd.merge(t2g_start[['tid','ORF_start','genomic_start']], t2g_end[['tid','ORF_end','genomic_end']])
    t2g = t2g.explode(['ORF_start','ORF_end','genomic_start','genomic_end']).drop_duplicates()
    df = pd.merge(df[['tid','start_codon','ORF_start','ORF_end','ORF_range','ORF_type','overlap']], t2g)
    
    # remove duplicated ORFs
    pbar = tqdm.pandas(desc="finding duplicated ORFs", unit_scale = True)
    orfs = df.groupby('genomic_start')['ORF_type'].progress_apply(list).reset_index()
    orfs.columns = ['genomic_start','ORF_list_start']
    df = pd.merge(df,orfs)
    df['ORF_list_unique_number'] = df.ORF_list_start.apply(lambda x: len(set(x)))
    df = pd.concat([df[df.ORF_list_unique_number==1], df[(df.ORF_list_unique_number!=1) & (df.ORF_type=='mORF')]])
    orfs = df.groupby('genomic_end')['ORF_type'].progress_apply(list).reset_index()
    orfs.columns = ['genomic_end','ORF_list_end']
    df = pd.merge(df,orfs)
    df['ORF_list_unique_number'] = df.ORF_list_end.apply(lambda x: len(set(x)))
    df = pd.concat([df[df.ORF_list_unique_number==1], df[(df.ORF_list_unique_number!=1) & (df.ORF_type=='mORF')]])

    # pbar = tqdm.pandas(desc="creating ORF sets      ", unit_scale = True)
    df['ORF_set'] = df[['ORF_list_start','ORF_list_end']].values.tolist()
    df['ORF_set'] = df['ORF_set'].apply(lambda x: set([j for i in x for j in i]))
    df.drop(['ORF_list_start','ORF_list_end','ORF_list_unique_number'], axis=1, inplace=True)

    return df



<<<<<<< HEAD
def operon_distribution(op, displot_prefix):
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
    
    sns.set_style("ticks")
    g = sns.JointGrid(data=op, x='count', y='length', height=3)
    g.plot_joint(sns.scatterplot)
    g.plot_marginals(sns.boxplot, fliersize=2)
    g.ax_joint.set_xscale('log')
    g.ax_joint.set_yscale('log')
    g.set_axis_labels(xlabel='Number of ORFs per mRNA', ylabel='mRNA length')
    plt.savefig(displot_prefix + '.operon_scatter.pdf', bbox_inches='tight')
    
    logging.info('plotted the distribution of operons as ' + displot_prefix + '.operon_dist.pdf and' + displot_prefix + '.operon_scatter.pdf')
    


=======
>>>>>>> refs/remotes/origin/master
def operon_finder(tx_assembly, bed, outdir=None, delim=None, start_codon=["ATG", "CTG", "GTG", "TTG"]):

    """
    Predict operons from transcriptome.
    
    Input:
        * tx_assembly: transcript fasta file extracted using bedtools getfasta. Headers with genomic coordinates (required)
        * bed: converted from gff3, e.g. from NCBI Genome Assembly (required)
        * outdir: output directory (default: None)
        * delim: use :: for tx_assembly extracted using bedtools getfasta -name flag, as this appends name (column #4) to the genomic coordinates (default: None)
        * start_codon: any triplets (default: ATG, ACG, CTG, GTG and TTG)
    
    Output:
        * cds_range: as input for analyse_footprints
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
    
    pbar = tqdm.pandas(desc="finding all ORFs       ", unit_scale = True)
    df['ORF_range'] = df.seq.progress_apply(lambda x: top_3_frames(x ,start_codon))
    df = df.explode('ORF_range')
    df = df.dropna(axis=0).reset_index(drop=True)
    df['start_codon'] = df['ORF_range'].apply(lambda x: x[0])
    df['ORF_range'] = df['ORF_range'].apply(lambda x: [x[1],x[2]])
    df['ORF_start'] = df.ORF_range.apply(lambda x: x[0])
    df['ORF_end'] = df.ORF_range.apply(lambda x: x[1])
    df['ORF_length'] = df.ORF_end - df.ORF_start

    tx = pr.PyRanges(df[['Chromosome','Start','End','Strand','tid','start_codon','ORF_start','ORF_end','ORF_length']].drop_duplicates())

    # Get mORFs encoded by assembled transcripts
    cds = pd.read_csv(bed, sep='\t', header=None)
    cds = cds[[0,1,2,5,3]].drop_duplicates()
    cds.columns = ['Chromosome','Start','End','Strand','Name']
    cds_ = cds
    cds = pr.PyRanges(cds)
    cdstx = cds.join(tx).df.drop_duplicates()
    cds_plus = cdstx[(cdstx.Start_b+cdstx.ORF_start==cdstx.Start) & (cdstx.Start_b+cdstx.ORF_end==cdstx.End)]
    cds_minus = cdstx[(cdstx.End_b-cdstx.ORF_start==cdstx.End) & (cdstx.End_b-cdstx.ORF_end==cdstx.Start)]
    cds = pd.concat([cds_plus,cds_minus])
    cds['ORF_type'] = 'mORF'    
    cds = cds[['tid','start_codon','ORF_start','ORF_end','ORF_type']].drop_duplicates()
    
    # Resolve ORFs
    orf = cdstx[['tid','start_codon','ORF_start','ORF_end']].drop_duplicates().copy()
    
    # find sORFs
    df = pd.merge(orf, cds, on='tid')
    sorf = df[(df.ORF_start_x<df.ORF_start_y) | (df.ORF_start_x>=df.ORF_end_y)].copy() # & (df.ORF_end_x<=df.ORF_start_y)
    sorf.drop_duplicates(['tid','start_codon_x','ORF_start_x','ORF_end_x'], inplace=True)
    sorf = sorf[['tid','start_codon_x','ORF_start_x','ORF_end_x']]
    sorf.columns = ['tid','start_codon','ORF_start','ORF_end']
    sorf['ORF_type'] = 'sORF'
    
    # find overlapping ORFs (oORFs)
    orfs = pd.concat([cds[['tid','start_codon','ORF_start','ORF_end','ORF_type']], sorf])
    d = df[['tid','start_codon_x','ORF_start_x','ORF_end_x']]
    d.columns = ['tid','start_codon','ORF_start','ORF_end']
    oorf = pd.concat([orfs, d]).drop_duplicates(['tid','start_codon','ORF_start','ORF_end'], keep=False)
    oorf = oorf[['tid','start_codon','ORF_start','ORF_end','ORF_type']]
    oorf['ORF_type'] = 'oORF'

<<<<<<< HEAD
    df = pd.concat([cds[['tid','start_codon','ORF_start','ORF_end','ORF_type']], oorf, sorf])
=======
    df = pd.concat([orfs,oorf])
>>>>>>> refs/remotes/origin/master
    df.drop_duplicates(['tid','start_codon','ORF_start','ORF_end'], inplace=True)
       
    # remove ORFs that are in-frame with mORFs
    cds.columns = ['Chromosome','start_codon','Start','End','ORF_type']
    df.columns = ['Chromosome','start_codon','Start','End','ORF_type']
    cds = pr.PyRanges(cds)    
    df = pr.PyRanges(df)
    
    dfcds = df.join(cds).df
    dfcds['frame'] = (dfcds.Start-dfcds.Start_b)%3
    inframeorf = pd.concat([dfcds[(dfcds.frame==0) & (dfcds.ORF_type!='mORF')]])
    inframeorf = inframeorf[['Chromosome','start_codon','Start','End']]
    
    df = pd.concat([df.df,inframeorf]).drop_duplicates(['Chromosome','start_codon','Start','End'], keep=False)
    df.columns = ['tid','start_codon','ORF_start','ORF_end','ORF_type']
    df['ORF_range'] = df[['ORF_start','ORF_end']].values.tolist()

    # Export CDS_range
    fname = filename(tx_assembly, None, outdir)
    
    cds_range = df[df.ORF_type=='mORF'][['tid','ORF_start','ORF_end']].drop_duplicates()
    cds_range.columns = ['tid','CDS_start','CDS_end']
    cds_range.to_csv(fname + '.cds_range.txt', header=None, index=None, sep='\t')

    #  Plot the distribution of operons
    df_ = pr.PyRanges(df_)
    cds_ = pr.PyRanges(cds_)
    cdstx_ = cds_.join(df_)
<<<<<<< HEAD
    cdstx_.length = cdstx_.df['seq'].apply(len)
    op = cdstx_.df.value_counts(['tid','length']).reset_index()
    
    operon_distribution(op, fname)

    df.to_pickle(fname + '.operon_finder.pkl.gz')

    logging.info('saved operons and ORFs as ' + fname + '.operon_finder.pkl.gz')
    logging.info('saved CDS range as ' + fname + '.cds_range.txt')
    
    return cds_range, df
=======
 
    sns.displot(cdstx_.df.value_counts('tid').reset_index()['count'], height=3)
    plt.title('Operons predicted from transcriptome assembly')
    plt.xlabel('Number of ORFs per mRNA')
    plt.savefig(fname + '.operons.pdf', bbox_inches='tight') 

    logging.info('saved the distribution of operons as ' + fname + '.operons.pdf')
    logging.info('saved CDS range as ' + fname + '.cds_range.txt')
    
    return cds_range, fname + '.cds_range.txt', df
>>>>>>> refs/remotes/origin/master
