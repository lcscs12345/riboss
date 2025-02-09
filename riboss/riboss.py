#!/usr/bin/env python
# coding: utf-8

"""
@author      CS Lim
@create date 2020-09-15 17:40:16
@modify date 2025-02-09 20:11:29
@desc        Main RIBOSS module
"""




import os, re, sys, time, argparse, textwrap, csv, logging.config, subprocess
from urllib.request import HTTPError
import numpy as np
import pandas as pd
import seaborn as sns
import seaborn.objects as so
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.ticker import MaxNLocator
import scipy.stats as stats
from scipy.stats.contingency import odds_ratio
from statsmodels.stats.multitest import multipletests
from types import SimpleNamespace
from collections import namedtuple
from tqdm import tqdm
import pyranges as pr
from Bio.Blast import NCBIWWW
from Bio import Entrez as ez
from pyfaidx import Faidx
from io import StringIO
from .chisquare import chi_square_posthoc
from .orfs import orf_finder, translate, CODON_TO_AA
from .wrapper import filename, merge_scores
from .read import read_blast_xml
import cloudpickle
import itertools


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




def base_to_bedgraph(superkingdom, base, bedgraph_prefix, profile=None, genepred=None, ncrna=False, delim=None, outdir=None):
    """
    Convert a dataframe from parse_ribomap to BedGraph format.

    Input:
        * superkingdom: Archaea, Bacteria or Eukaryota (required)
        * base: dataframe from parse_ribomap (required)
        * bedgraph_prefix: output filename prefix for BedGraph (required)
        * delim: use :: for tx_assembly extracted using bedtools getfasta -name flag, as this appends name (column #4) to the genomic coordinates (default: None)
        * outdir: output directory (default=None)
    Output:
        * BedGraph file
    """
    
    fname = filename(bedgraph_prefix, 'riboprof', outdir)
    df = base.copy()
    
    if superkingdom in ['Archaea','Bacteria']:
        if delim!=None:
            pos = 1
        else:
            pos = 0
        
        try:
            df['Chromosome'] = df.tid.str.split(':').str[pos+1]
            df['Start'] = df.tid.str.split(':').apply(lambda x: x[pos+2]).str.split('-').apply(lambda x: x[0]).astype(int)
            df['End'] = df.tid.str.split(':').apply(lambda x: x[pos+2]).str.split('-').apply(lambda x: x[1]).str.split('(').apply(lambda x: x[0]).astype(int)
            df['Strand'] = df.tid.str.split(':').apply(lambda x: x[pos+2]).str.split('(').apply(lambda x: x[1]).str.replace(')','')
            
            df['range'] = df[['Start','End']].values.tolist()
            df['range'] = df['range'].apply(lambda x: list(range(x[0],x[1])))
            df_plus = df[df.Strand=='+'].copy()
            df_minus = df[df.Strand=='-'].copy()
            df_minus['range'] = df_minus.range.apply(lambda x: list(reversed(x)))
            df = pd.concat([df_plus, df_minus])
            
            df = df.drop(['tid','Start','End'],axis=1).explode(['range',profile])
            df['End'] = df['range'] +1
            df.rename(columns={'range':'Start'}, inplace=True)
            df[profile] = df[profile].astype(int)
        except:
            logging.error('Check if delim is correct!')
        
    elif superkingdom=='Eukaryota':
        gp = pd.read_csv(genepred, sep='\t', header=None)
        
        if ncrna==False:
            gp = gp[gp[5]!=gp[6]]
            
        gp[8] = gp[8].str.split(',').str[:-1].apply(lambda x: [int(i) for i in x])
        gp[9] = gp[9].str.split(',').str[:-1].apply(lambda x: [int(i) for i in x])
        gp['exons'] = gp[[8,9]].values.tolist()
        gp['exons'] = gp.exons.apply(lambda x: list(zip(x[0],x[1])))
        gp['exons'] = gp.exons.apply(lambda x: [list(range(i[0],i[1])) for i in x])
        gp['exons'] = gp.exons.apply(lambda x: [j for i in x for j in i])
        gc = gp[[0,1,2,'exons']].rename(columns={0:'tid',1:'Chromosome',2:'Strand'})

        df = pd.merge(df,gc).copy()
        df['range'] = df.apply(lambda x: list(zip(x[profile],x['exons'])), axis=1)
        df = df.explode('range')[['tid','Chromosome','Strand','range']]
        df['Start'] = df['range'].apply(lambda x: x[1]).astype(int)
        df['End'] = df['range'].apply(lambda x: x[1]).astype(int) +1
        df[profile] = df['range'].apply(lambda x: x[0]).astype(int)
        df.drop('range', axis=1, inplace=True)               
    
    else:
        logging.error('Please check your spelling! Only Archaea, Bacteria, or Eukaryota is acceptable superkingdom.')
        
    # Split by strands
    df = df[df[profile]>0].copy()
    plus = df[df.Strand=='+'][['Chromosome','Start','End',profile]]    
    minus = df[df.Strand=='-'][['Chromosome','Start','End',profile]]
    df_merged = pd.merge(plus,minus, on=['Chromosome','Start','End'])
    
    if df_merged.shape[0]!=0:
        plus.to_csv(fname + '.plus.bg', index=None, header=None, sep='\t')
        bg = merge_scores(fname + '.plus.bg')
        f = open(fname + '.plus.bg', 'w')
        f.write('track type=bedGraph name="Plus strand" description="Ribosome profile" visibility=full color=0,0,0 priority=20\n')
        bg.to_csv(f, index=None, header=None, sep='\t', mode='a')
        f.close()
        
        minus.to_csv(fname + '.minus.bg', index=None, header=None, sep='\t')
        bg_ = merge_scores(fname + '.minus.bg')
        f = open(fname + '.minus.bg', 'w')
        f.write('track type=bedGraph name="Minus strand" description="Ribosome profile" visibility=full color=0,0,0 priority=20\n')
        bg_.to_csv(f, index=None, header=None, sep='\t', mode='a')
        f.close()

        bg = pd.concat([bg,bg_])

        if bg[bg[0].str.contains('ERROR')].shape[0]!=0:
            logging.error('Failed generating BedGraph!')        
        else:
            logging.info('saved ribosome profiles as ' + fname + '.plus.bg and ' + fname + '.minus.bg')
            
    else:
        df[['Chromosome','Start','End',profile]].to_csv(fname + '.bg', index=None, header=None, sep='\t')
        bg = merge_scores(fname + '.bg')
        f = open(fname + '.bg', 'w')
        f.write('track type=bedGraph name="riboprof" description="Ribosome profile" visibility=full color=0,0,0 priority=20\n')
        bg.to_csv(f, index=None, header=None, sep='\t', mode='a')
        f.close()
        
        if bg[bg[0].str.contains('ERROR')].shape[0]!=0:
            logging.error('Failed generating BedGraph!')
        else:
            logging.info('saved ribosome profiles as ' + fname + '.bg')
    
    bg.columns = ['Chromosome','Start','End',profile]
    
    if superkingdom in ['Archaea','Bacteria']:
        return bg
    elif superkingdom=='Eukaryota':
        return gc, bg
    else:
        logging.error('Please check your spelling! Only Archaea, Bacteria, or Eukaryota is acceptable superkingdom.')


def parse_ribomap(superkingdom, base, genepred=None, ncrna=False, delim=None, outdir=None):
    """
    Parse a ribomap/riboprof base file into a dataframe.
    
    Input:
        * superkingdom: Archaea, Bacteria or Eukaryota (required)
        * base: output base file from riboprof (required)
        * genepred: required for Eukaryota
        * ncrna: remove noncoding RNAs from analysis
        * delim: use :: for tx_assembly extracted using bedtools getfasta -name flag, as this appends name (column #4) to the genomic coordinates (default: None)
        * outdir: output directory (default=None)
    
    Output:
        * a dataframe
    """
    
    dt = pd.read_csv(base, sep=r':\s', engine='python', header=None)

    profiles = []
    for p in ['ribo profile','mRNA profile']:
        profile = list(p.split()[0])[0] + p.split()[1]
        
        prof = dt[(dt[0]=='tid') | (dt[0]==p)].copy()
        prof.columns = ['id','val']
        prof = pd.DataFrame({'tid':prof['val'].iloc[::2].values, profile:prof['val'].iloc[1::2].values})
    
        pbar = tqdm.pandas(desc="parsing ribomap output ", unit_scale=True, ncols=100)
        prof[profile] = prof[profile].str.split().progress_apply(lambda x: [float(i) for i in x])
        profiles.append(prof)

    bgs = []
    for p in profiles:
        if superkingdom in ['Archaea','Bacteria']:
            bg = base_to_bedgraph(superkingdom, p, base, profile=p.columns[1], ncrna=ncrna, delim=delim, outdir=outdir)
            bgs.append(bg)
            return bgs, profiles
            
        elif superkingdom=='Eukaryota':
            gc, bg = base_to_bedgraph(superkingdom, p, base, profile=p.columns[1], genepred=genepred, ncrna=ncrna, delim=None, outdir=outdir)
            bgs.append(bg)
            return gc, bgs, profiles
        else:
            logging.error('Please check your spelling! Only Archaea, Bacteria, or Eukaryota is acceptable superkingdom.')


def footprint_counts(df, d, gencode=None):
    """
    Count footprints by reading frames.
    
    Input:
        * df: dataframe from orf_finder/operon_finder 
        * d: dataframe from parse_ribomap
    Output:
        * dt: dataframe with lists of footprint counts by frames.
    """
    
    if 'ORF_range' not in df.columns:
        df['ORF_range'] = df[['ORF_start','ORF_end']].values.tolist()
            
    dt = pd.merge(df, d, on='tid')
    
    if gencode:
        dt['tid'] = dt.tid.str.split('|').str[0]
    
    # create frame 0
    dt['frame0'] = dt.ORF_range.apply(lambda x: [i for i in range(x[0], x[1], 3)])
    # create frame 1
    dt['frame1'] = dt.ORF_range.apply(lambda x: [i for i in range(x[0]+1, x[1], 3)])
    # create frame 2
    dt['frame2'] = dt.ORF_range.apply(lambda x: [i for i in range(x[0]+2, x[1], 3)])

    if 'rprofile' in dt.columns:
        profile = 'rprofile'
    else:
        profile = 'mprofile'

    pbar = tqdm.pandas(desc="counting footprints    ", unit_scale=True, ncols=100)
    dt['periodicity'] = dt[[profile,'frame0','frame1','frame2']].values.tolist()
    dt['periodicity'] = dt.periodicity.progress_apply(lambda x: [sum([x[0][i] for i in tuple(x[1])]),
                                                            sum([x[0][i] for i in tuple(x[2])]),
                                                            sum([x[0][i] for i in tuple(x[3])])])
    
    # get rprofile at start codons
    dt['start_pos'] = dt.ORF_range.apply(lambda x: [x[0], x[0]+1, x[0]+2])
    dt['start_pos'] = dt[[profile,'start_pos']].values.tolist()
    dt['start_rprofile'] = dt.start_pos.apply(lambda x: [x[0][i] for i in x[1]])

    # filter rprofile by frames. At least 10 footprints mapped to all 3 frames, and the number of footprints mapped to frame0 must be > frame1.
    dt['frame0'] = dt.periodicity.apply(lambda x: x[0])
    dt['frame1'] = dt.periodicity.apply(lambda x: x[1])
    dt['frame2'] = dt.periodicity.apply(lambda x: x[2])
    dt = dt[(dt.frame0 + dt.frame1 + dt.frame2>10) & (dt.frame0>dt.frame1) & (dt.frame0>dt.frame2)].reset_index(drop=True)
    
    # prioritise start_rprofile>0
    d1 = dt[dt.start_rprofile.apply(lambda x: x[0])>0].sort_values(['start_codon','ORF_start','tid'])
    d0 = dt[dt.start_rprofile.apply(lambda x: x[0])==0].sort_values(['start_codon','ORF_start','tid'])
    dt = pd.concat([d1, d0])
    dt['tab'] = dt[['frame0','frame1','frame2']].values.tolist()

    dt.drop_duplicates(['tid','ORF_type','ORF_end'], inplace=True)
    dt.drop([profile,'periodicity','start_pos','frame0','frame1','frame2'], axis=1, inplace=True)  
    
    return dt



def contingency_tables(df, tab_x='tab_x', tab_y='tab_y'):
    """
    Create 2x3 and 2x2 contingency tables.

    Input: 
        * df: dataframe with two lists (tab_x and tab_y)
        
    Output:
        * df: dataframe with contingency tables
    """
    # create 2x3 contingency tables
    df['tab'] = df[[tab_x,tab_y]].values.tolist()
    df['tab'] = df['tab'].apply(lambda x: np.nan_to_num(np.array(x)).astype(int))
    # create 2x2 contingency tables, add pseudocount of 1 to zeroes to avoid odds ratio infinity
    df['tab_x_'] = df[tab_x].apply(lambda x: [x[0],x[1]+x[2]] if (x[1]+x[2]>0) else [x[0],1])
    df['tab_y_'] = df[tab_y].apply(lambda x: [x[0],x[1]+x[2]] if (x[1]+x[2]>0) else [x[0],1])
    df['tab_'] = df[['tab_x_','tab_y_']].values.tolist()
    df['tab_'] = df['tab_'].apply(lambda x: np.nan_to_num(np.array(x)).astype(int))
    
    df['tab'] = df[['tab','tab_']].values.tolist()
    df.drop([tab_x,tab_y,'tab_x_','tab_y_'], axis=1, inplace=True)
    
    return df
    


def statistical_test(df, num_simulations=1000, padj_method='fdr_bh'):

    logging.disable(logging.INFO)
    # calculate odds ratios
    df['odds_ratio'] = df['tab'].apply(lambda x: odds_ratio(x[1]))
    
    # perform (chi-square test, num_simulations bootstrap, post-hoc test), or (boschloo-exact test) if counts were zeros
    pbar = tqdm.pandas(desc="comparing periodicity  ", ncols=100)#, mininterval=10, maxinterval=200, miniters=int(df.shape[0]/10))
    
    f = df[(df.tab.apply(lambda x: x[1][0][0]>x[1][0][1])) | (df.tab.apply(lambda x: x[1][1][0]>x[1][1][1]))].copy()
    f['result'] = f.tab.progress_apply(lambda x: chi_square_posthoc(x[0],num_simulations) if np.sum(x[0], axis=0).all()==True else stats.boschloo_exact(x[1]))
    f['statistical test'] = f['result'].apply(lambda x: 'ChiSquare' if type(x)==SimpleNamespace else 'BoschlooExact')
    
    # p-value correction
    pval = f['result'].apply(lambda x: x.pvalue).tolist()
    f['adj_pval'] = multipletests(pval, method=padj_method)[1]
    
    chi = f[f['statistical test']=='ChiSquare'].reset_index(drop=True)
    chi['statistic'] = chi['result'].apply(lambda x: x.statistic)
    chi['posthoc_statistic'] = chi['result'].apply(lambda x: x.posthoc_statistic)
    bootstrap_pval = chi['result'].apply(lambda x: x.bootstrap_pvalue).tolist()
    chi['adj_bootstrap_pval'] = multipletests(bootstrap_pval, method=padj_method)[1]
    posthoc_pval = np.concatenate(chi['result'].apply(lambda x: x.posthoc_pvalue).tolist()).ravel()
    chi['adj_posthoc_pval'] = multipletests(posthoc_pval, method=padj_method)[1].reshape(chi.shape[0], 3).tolist()
    
    ChiSquareResult = namedtuple('ChiSquareResult', ['statistic', 'adjusted_pvalue', 'adjusted_bootstrap_pvalue', 'posthoc_statistic','adjusted_posthoc_pvalue'])
    chi['result'] = chi[['statistic','adj_pval','adj_bootstrap_pval','posthoc_statistic','adj_posthoc_pval']].values.tolist()
    chi['result'] = chi['result'].apply(lambda x: ChiSquareResult(x[0],x[1],x[2],x[3],x[4]))
    
    be = f[f['statistical test']=='BoschlooExact'].reset_index(drop=True)
    be['statistic'] = be['result'].apply(lambda x: x.statistic)
    be['result'] = be[['statistic','adj_pval']].values.tolist()
    
    BoschlooExactResult = namedtuple('BoschlooExactResult', ['statistic', 'adjusted_pvalue'])
    be['result'] = be['result'].apply(lambda x: BoschlooExactResult(x[0],x[1]))
    
    if be.shape[0]>0:
        f = pd.concat([chi,be])
    else:
        f = chi
    
    xf = df[(df.tab.apply(lambda x: x[1][0][0]<x[1][0][1])) & (df.tab.apply(lambda x: x[1][1][0]<x[1][1][1]))].copy()

    logging.disable(logging.NOTSET)
    
    return f, xf



def boss(df, tx_assembly, boss_prefix, padj_method='fdr_bh', tie=False, num_simulations=1000, outdir=None):  
    """
    Compare uORFs, oORFs and dORFs to mORFs, and assign which ORFs are the bosses.
    
    Input:
        * df: dataframe from footprint_counts (required)
        * tx_assembly: transcript fasta file extracted using bedtools getfasta. Headers with genomic coordinates (required)
        * tie: if adjusted p-values between ORFs is not significant (default=False)
        * num_simulations: number of simulations (default=1000)
        * outdir: output directory (default=None)
    
    Output:
        * sig: dataframe for significant RIBOSS results
        * boss: dataframe for all RIBOSS statistics
        * RIBOSS statistics and significant results as Pickle/BED 
    """

    m = df[df.ORF_type=='mORF'].reset_index(drop=True)
    o = df[df.ORF_type!='mORF'].reset_index(drop=True)
    f = pd.merge(o, m, on='tid', how='outer')

    # ORFs win by default wihout opponents
    fo = f[f.ORF_type_y.isna()].copy()
    fo['boss'] = 'default'
    fo['odds_ratio'] = np.inf
    fo = fo[['tid','boss','start_codon_x','ORF_range_x','ORF_type_x','tab_x','odds_ratio','start_rprofile_x']]
    
    fm = f[f.ORF_type_x.isna()].copy()
    fm['boss'] = 'default'
    fm['odds_ratio'] = 0
    fm = fm[['tid','boss','start_codon_y','ORF_range_y','ORF_type_y','tab_y','odds_ratio','start_rprofile_y']]

    # ORFs with opponents
    f.dropna(inplace=True)
    f = contingency_tables(f)
    f, xf = statistical_test(f, num_simulations, padj_method=padj_method)

    # Are mORFs the boss?
    f.reset_index(drop=True, inplace=True)

    boss_x = f[(f.odds_ratio.apply(lambda x: x.statistic>1)) & (f['result'].apply(lambda x: x.adjusted_pvalue<0.05)) & (f['result'].apply(lambda x: (x.adjusted_bootstrap_pvalue<0.05) if 'adjusted_bootstrap_pvalue' in x._fields else x.adjusted_pvalue<0.05))].copy()
    boss_x['boss'] = boss_x['ORF_type_x']
    boss_y = f[(f.odds_ratio.apply(lambda x: x.statistic<1)) & (f['result'].apply(lambda x: x.adjusted_pvalue<0.05)) & (f['result'].apply(lambda x: (x.adjusted_bootstrap_pvalue<0.05) if 'adjusted_bootstrap_pvalue' in x._fields else x.adjusted_pvalue<0.05))].copy()
    boss_y['boss'] = boss_y['ORF_type_y']
    boss = pd.concat([boss_x,boss_y,f]).reset_index().drop_duplicates('index')
    boss['boss'] = boss['boss'].fillna('tie')
    boss = pd.concat([boss,xf]).fillna('lacks periodicity')
    boss['odds_ratio'] = boss.odds_ratio.apply(lambda x: x.statistic)    
    boss = pd.concat([boss,fo,fm]).reset_index(drop=True)
    boss = boss[['tid', 'boss', 'start_codon_x', 'start_codon_y', 'ORF_range_x', 'ORF_type_x', 'start_rprofile_x', 'ORF_range_y', 'ORF_type_y','start_rprofile_y', 'tab', 'odds_ratio', 'statistical test', 'result']]
    # boss.dropna(subset=['odds_ratio'],inplace=True)
    # boss.reset_index(drop=True, inplace=True)


    fname = filename(boss_prefix, 'riboss', outdir)
    boss.to_csv(fname + '.csv', index=None)
    
    with open(fname + '.boss.pkl', 'wb') as fl:
        cloudpickle.dump(boss, fl)
        logging.info('saved RIBOSS stats as ' + fname + '.boss.pkl and ' + fname + '.boss.csv')
        
    # extract significant results        
    chi = boss[(boss['statistical test']=='ChiSquare') & (boss.boss!='tie')].reset_index(drop=True)
    be = boss[(boss['statistical test']=='BoschlooExact') & (boss.boss!='tie')].reset_index(drop=True)
    
    if tie==True:
        chi = chi[(chi.result.apply(lambda x: x.adjusted_posthoc_pvalue[0])<0.05) & (chi.boss!='mORF')].reset_index(drop=True)
        be = be[be.result.apply(lambda x: x.adjusted_pvalue<0.05)].reset_index(drop=True)
    else:
        chi = chi[(chi.result.apply(lambda x: x.adjusted_posthoc_pvalue[0])<0.05) & (chi.boss!='tie') & (chi.boss!='mORF')].reset_index(drop=True)
        be = be[be.result.apply(lambda x: x.adjusted_pvalue<0.05) & (be.boss!='tie')].reset_index(drop=True)
        
    sig = pd.concat([chi,be])
    if sig.shape[0]>0:    
        sig['start'] = sig.ORF_range_x.apply(lambda x: x[0])
        sig['end'] = sig.ORF_range_x.apply(lambda x: x[1])
        sig.drop_duplicates(['tid','start','end'], inplace=True)
    
        sig.rename(columns={'tid':'Chromosome','start':'Start','end':'End'}, inplace=True)
        sig = pr.PyRanges(sig)
        sig.dna = pr.get_sequence(sig, tx_assembly)
        sig.aa = sig.dna.apply(lambda x: translate(x))
        sig.oid = sig.Chromosome.astype(str) + '__' + sig.Start.astype(int).astype(str) + '-' + sig.End.astype(int).astype(str)
        sig.fa = ('>' + sig.oid + '\n' + sig.aa).tolist()
        sig = sig.df
    
        sig.rename(columns={'Chromosome':'tid','Start':'start','End':'end'}, inplace=True)
        sig.reset_index(drop=True, inplace=True)

        with open(fname + '.sig.pkl', 'wb') as fl:
            cloudpickle.dump(sig, fl)
            logging.info('saved significant RIBOSS results (n=' + str(sig.shape[0]) + ') as ' + fname + '.sig.pkl')        

    return sig, boss



def blastp(sig, blastp_prefix, email=None, outdir=None):
    """
    BLASTP for ORFs with significantly greater triplet periodicity than mORFs.

    Input:
        * sig: significant RIBOSS hits (required)
        * blastp_prefix: output filename (required)
        * outdir: output directory (default: None)

    Output:
        * df: dataframe for BLASTP results
        * BLASTP results as xml
    """

    start_time = time.perf_counter()
    
    if email!=None:
        NCBIWWW.email=email

    logging.info('perform BLASTP for RIBOSS hits (n=' + str(sig.shape[0]) + ')')
    result_handle = NCBIWWW.qblast('blastp', 'nr', '\n'.join(sig.fa.tolist()))
    
    fname = filename(blastp_prefix, 'riboss', outdir)
    
    with open(fname + '.sig.blastp.xml', 'w') as save_to:
        save_to.write(result_handle.read())
        result_handle.close()

    logging.info('finished BLASTP in ' + 
          str(round((time.perf_counter()-start_time)/60)) + ' min ' +
          str(round((time.perf_counter()-start_time)%60)) + ' s')
    logging.info('saved BLASTP results for RIBOSS hits as ' + fname + '.sig.blastp.xml')

    df = read_blast_xml(fname + '.sig.blastp.xml')   
    
    return df



def efetch(acc, tries=5, sleep=1, email=None, api_key=None):
    """
    Fetch Identical Protein Groups (IPG) from NCBI.
    
    Input: 
        * acc: a single or a list or 1D array of NCBI accession numbers (required)
        * tries: number of attempts (default=5) 
        * sleep: pause time seconds (default=1)
        * email: email address (default=None)
        * api_key: NCBI API key https://support.nlm.nih.gov/knowledgebase/article/KA-05317/en-us (default=None)

    Output:
        * ipg: dataframe for IPG
    """
    
    start_time = time.perf_counter()
    
    if email!=None:
        ez.email = email
    if api_key!=None:
        ez.api_key = api_key
    
    logging.info('efetch IPG for ' + str(len(acc)) + ' accession numbers')
    
    if type(acc)==str:
        
        tries = tries
        for n in range(tries):
            try:
                # logging.info('efetching IPG for ' + acc)
                e = ez.efetch(db="protein", id=acc, rettype="ipg", retmode="txt")
                ipg = pd.read_csv(StringIO(e.read().decode('utf-8')), sep='\t', header=None)
            except HTTPError:
                if n < tries - 1: # i is zero indexed
                    continue
                else:
                    raise
            break
    
    elif (type(list(acc))==list) | (type(acc)==np.ndarray):
        entries = []
        for i in acc:
            
            tries = tries
            for n in range(tries):
                try:
                    # logging.info('efetching IPG for ' + i)
                    e = ez.efetch(db="protein", id=i, rettype="ipg", retmode="txt")
                    e = pd.read_csv(StringIO(e.read().decode('utf-8')), sep='\t', header=None)
                    entries.append(e)
                    time.sleep(sleep)
                except HTTPError:
                    if n < tries - 1: # i is zero indexed
                        continue
                    else:
                        raise
                break
        
        ipg = pd.concat(entries)
    
    logging.info('finished efetch in ' + 
          str(round((time.perf_counter()-start_time)/60)) + ' min ' +
          str(round((time.perf_counter()-start_time)%60)) + ' s')
    
    return ipg



def predicted_orf_profile(rp, df, utr, title, barplot_prefix):
    """
    Plot metagene ribosome profiles for unannotated ORFs.
    
    Input:
        * rp: ribosome profiles for ORFs extracted from riboprof (required)
        * df: merged dataframe from blastp and riboprof (required)
        * title: (required)
        * barplot_prefix: filename prefix for barplots (required)

    Output:
        * metagene plots as PDFs
    """

    infix = title.lower().split()[0] + '_' + '_'.join(title.lower().split()[2:])
    
    peptide_len = np.median(df.aa.apply(len))
    if round(peptide_len)==peptide_len:
        peptide_len = round(peptide_len)
    elif round(peptide_len)>peptide_len:
        peptide_len

    ymax = rp[(rp['Position from predicted start codon']>=-12) & (rp['Position from predicted start codon']<=30)]['Ribosome profile from start codon'].tolist()
    ymax = np.max(ymax) + np.max(ymax)*0.05
    
    plt.set_loglevel('WARNING')
    plt.figure(figsize=(4,3))
    ax = sns.barplot(data=rp, x='Position from predicted start codon', y='Ribosome profile from start codon', hue='Frames', width=1)
    ax.xaxis.set_major_locator(ticker.MultipleLocator(3))
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    plt.title(title + ' (n=' + str(df.shape[0]) + ', median length=' + str(peptide_len) + ' aa)')
    plt.ylabel('Metagene ribosome profile')
    plt.xlim(-12+utr, 30+utr)
    plt.ylim(0,ymax)
    plt.savefig(barplot_prefix + '.' + infix + '.start_codon.pdf', bbox_inches='tight')
            
    logging.info('saved metagene plots as ' + barplot_prefix + '.' + infix + '.start_codon.pdf')
    


def operons_to_biggenepred(orf, df, bed, fai, big_fname, delim=None):
    """
    Create UCSC bigGenePred for RIBOSS top hits for Acheae and Bacteria.

    Input:
        * orf: dataframe from operon_finder (required)
        * df: dataframe for RIBOSS hits (required)
        * bed: genome BED file (required)
        * fai: genome fasta index to generate chrom.sizes (required)
        * big_fname: output filename prefix for bigGenePred
        * delim: use :: for tx_assembly extracted using bedtools getfasta -name flag, as this appends name (column #4) to the genomic coordinates (default: None)
        
    Output:
        * bigGenePred for RIBOSS hits 
    """
    
    orf = orf.rename(columns={'start_codon':'start_codon_x','ORF_start':'start','ORF_end':'end'}).drop('ORF_length', axis=1)
    toporf = pd.merge(orf, df)
    toporf['name2'] = toporf.Chromosome.astype(str) + ':' + toporf['Start'].astype(int).astype(str) + '-' + toporf['End'].astype(int).astype(str) + '(' + toporf.Strand.astype(int).astype(str) + ')' 
    toporf['reserved'] = '255,128,0'
    toporf['blockCount'] = 1
    toporf['blockSizes'] = toporf.ORF_range_x.apply(lambda x: x[1]-x[0])
    toporf['chromStarts'] = 0
    toporf['cdsStartStat'] = 'none'
    toporf['cdsEndStat'] = 'none'
    toporf['exonFrames'] = 0
    toporf['type'] = 'none'
    toporf['geneType'] = 'none'
    
    # BLASTP hits
    ob = toporf[~toporf.title.isna()].copy()
    ob['Name'] = ob.title.str.split().str[0].str.split('|').str[1] + '__' + ob.ORF_type_x + '-BLASTP'
    
    # No hits
    on = toporf[toporf.title.isna()].copy()
    on['bits'] = 0
    
    cds = pd.read_csv(bed, sep='\t', header=None)
    # cds_ = cds
    cds = cds[[0,1,2,5,3]].drop_duplicates()
    cds.columns = ['Chromosome','Start','End','Strand','Name']
    
    # oORFs
    ono = pr.PyRanges(on).join(pr.PyRanges(cds), strandedness='same').df
    ono['Strand'] = ono.Strand.astype(str)
    ono['Name'] = ono.name2 + '__' + ono.ORF_type_x + '-Ref'
    
    # sORFs
    ons = pd.concat([ono,on]).drop_duplicates(['Chromosome','Start','End','Strand'], keep=False)
    ons['Name'] = ons.name2 + '__' + ons.ORF_type_x

    # Create a dataframe for pre-bigGenePred
    bb = pd.concat([ob,ono,ons])[['Chromosome','Start','End','Name','bits','Strand','Start','End','reserved','blockCount','blockSizes',
     'chromStarts','name2','cdsStartStat','cdsEndStat','exonFrames','type','Name','name2','geneType','ORF_type_x']]
    bb['bits'] = bb.bits.astype(int)
    bb.drop_duplicates(inplace=True)
    bb.columns = ['Chromosome','Start','End','Name','bits','Strand','thickStart','thickEnd','reserved','blockCount','blockSizes',
                  'chromStarts','name2','cdsStartStat','cdsEndStat','exonFrames','type','geneName','geneName2','geneType','ORF_type_x']
    
    # Prepare required files for making bigGenePred
    faidx = pd.read_csv(fai, sep='\t', header=None)
    chromsizes = os.path.splitext(fai)[0] + '.chrom.sizes'
    faidx[[0,1]].to_csv(chromsizes, header=None, index=None, sep='\t')

    bas = os.path.split(big_fname)[0] + '/bigGenePred.as'
    subprocess.run(['wget','https://genome.ucsc.edu/goldenpath/help/examples/bigGenePred.as',
                    '-O', bas], check=True)
    
    # Create bigGenePred by ORF_type
    for orftype in bb.ORF_type_x.unique():
        bed = big_fname + '.' + orftype + '.bed'
        bb[bb.ORF_type_x==orftype].drop('ORF_type_x',axis=1).to_csv(bed, index=None, header=None, sep='\t')
        
        pre = big_fname + '.' + orftype + '.sorted.bed'
        subprocess.run(['bedSort', bed, pre], check=True)
            
        bgp = big_fname + '.' + orftype + '.bb'
        subprocess.run(['bedToBigBed','-as=' + bas,
                        '-type=bed12+8', pre, chromsizes, bgp], check=True)
    
        os.remove(bed)        
        os.remove(pre)
        
        if os.path.exists(bgp):
            logging.info('saved bigGenePred for RIBOSS hits as ' + bgp)
        
    os.remove(chromsizes)
    os.remove(bas)




def group_ranges(L):
    """
    https://jonathanchang.org/blog/group-consecutive-number-ranges/
    Collapses a list of integers into a list of the start and end of
    consecutive runs of numbers. Returns a generator of generators.
    >>> [list(x) for x in group_ranges([1, 2, 3, 5, 6, 8])]
    [[1, 3], [5, 6], [8]]
    """
    for w, z in itertools.groupby(L, lambda x, y=itertools.count(): next(y)-x):
        grouped = list(z)
        yield (x for x in [grouped[0], grouped[-1]+1][:len(grouped)])



def orfs_to_biggenepred(orf_ranges, df, fai, big_fname, orf_range_col=None, orf_type_col=None):
    """
    Create UCSC bigGenePred for RIBOSS top hits for Eukaryota.

    Input:        
        * orf_ranges: gc dataframe from parse_ribomap (required)
        * df: dataframe for RIBOSS hits (required)
        * fai: genome fasta index to generate chrom.sizes (required)
        * big_fname: output filename prefix for bigGenePred

    Output:
        * bigGenePred for RIBOSS hits
    """

    df = df[~df[orf_type_col].isna()].copy()
    df_ranges = pd.merge(orf_ranges,df)
    dfp = df_ranges[df_ranges.Strand=='+'].copy()
    dfp['Start'] = dfp.apply(lambda x: x['exons'][x[orf_range_col][0]], axis=1)
    dfp['End'] = dfp.apply(lambda x: x['exons'][x[orf_range_col][1]-1] +1, axis=1)
    dfm = df_ranges[df_ranges.Strand=='-'].copy()
    dfm['End'] = dfm.apply(lambda x: sorted(x['exons'], reverse=True)[x[orf_range_col][0]] +1, axis=1)
    dfm['Start'] = dfm.apply(lambda x: sorted(x['exons'], reverse=True)[x[orf_range_col][1]-1], axis=1)
    df_ranges = pd.concat([dfp, dfm])

    b = df_ranges[['Chromosome','tid','Strand',orf_type_col,orf_range_col,'exons','Start','End']].copy()
    b = b.explode('exons')
    b['oid'] = b.tid + '__' + b[orf_range_col].str[0].astype(int).astype(str) + '-' + b[orf_range_col].str[1].astype(int).astype(str)
    b = b[(b.exons>=b.Start) & (b.exons<b.End)].drop(orf_range_col, axis=1).copy()
    
    bgroupby = b.groupby(['Chromosome','tid','Strand',orf_type_col,'Start','End','oid'])['exons'].apply(list).reset_index()
    bgroupby['exons'] = bgroupby.exons.apply(lambda x: [list(i) for i in group_ranges(set(x))])
    bgroupby = bgroupby.explode('exons')
    
    b1 = bgroupby[bgroupby.exons.apply(len)==1].copy()
    b1['Starts'] = b1.exons.str[0]
    b1['Ends'] = b1.exons.str[0] +1
    b1.drop('exons', axis=1, inplace=True)
    
    b2 = bgroupby[bgroupby.exons.apply(len)==2].copy()
    b2['Starts'] = b2.exons.apply(lambda x: x[0])
    b2['Ends'] = b2.exons.apply(lambda x: x[1])
    b2.drop('exons', axis=1, inplace=True)
    
    b = pd.concat([b1,b2])
    starts = b.groupby(['Chromosome','tid','Strand',orf_type_col,'Start','End','oid'])['Starts'].apply(list).reset_index()
    ends = b.groupby(['Chromosome','tid','Strand',orf_type_col,'Start','End','oid'])['Ends'].apply(list).reset_index()
    
    boss_gp = pd.merge(starts,ends)
    boss_gp['exonCount'] = boss_gp['Starts'].apply(len)
    boss_gp['Starts'] = boss_gp['Starts'].apply(lambda x: ','.join([str(i) for i in sorted(x)])) + ','
    boss_gp['Ends'] = boss_gp['Ends'].apply(lambda x: ','.join([str(i) for i in sorted(x)])) + ','
    boss_gp = boss_gp[['oid','Chromosome','Strand','Start','End','Start','End','exonCount','Starts','Ends',orf_type_col]]
    
    # Prepare required files for making bigGenePred
    faidx = pd.read_csv(fai, sep='\t', header=None)
    chromsizes = os.path.splitext(fai)[0] + '.chrom.sizes'
    faidx[[0,1]].to_csv(chromsizes, header=None, index=None, sep='\t')
    
    bpath = os.path.split(big_fname)[0]
    if (len(bpath)==0) & (os.path.exists(bpath)==False):
        bas =  'bigGenePred.as'
    elif os.path.exists(bpath)==True:
        bas = bpath + '/bigGenePred.as'
    else:
        try:
            os.makedirs(bpath)
            bas = bpath + '/bigGenePred.as'
        except Exception as e:
            logging.exception('No permission to create a new directory! Please check the output filename.')
        
    subprocess.run(['wget','https://genome.ucsc.edu/goldenpath/help/examples/bigGenePred.as',
                    '-O', bas], check=True)
    
    # Create bigGenePred by ORF_type
    for orftype in boss_gp[orf_type_col].unique():
        gp = big_fname + '.' + orftype + '.gp'
        bed = big_fname + '.' + orftype + '.bed'
        bgp = big_fname + '.' + orftype + '.bgp'
        boss_gp[boss_gp[orf_type_col]==orftype].drop(orf_type_col,axis=1).to_csv(gp, header=None, index=None, sep='\t')
        subprocess.run(['genePredToBed', gp, bed])
        subprocess.run(['genePredToBigGenePred', gp, bgp + '.pre'])
        subprocess.call('sort -k1,1 -k2,2n '.split() + [bgp + '.pre'], stdout=open(bgp, 'w'))
            
        bb = big_fname + '.' + orftype + '.bb'
        subprocess.run(['bedToBigBed','-as=' + bas,
                        '-type=bed12+8', '-tab', bgp, chromsizes, bb], check=True)
        
        if os.path.exists(bb):
            logging.info('saved bigGenePred for RIBOSS hits as ' + bb + ', ' + bed + ', and ' + gp)
            os.remove(bgp + '.pre')
            os.remove(bgp)
            
    os.remove(chromsizes)
    os.remove(bas)


    
def profile_anomaly(bedgraph, bb, bed, fasta, scatterplot_prefix=None):
    """
    Find anomalous ribosome profiles by the encoded amino acids for oORFs.
    
    Input:
        * bedgraph: dataframe for BedGraph from base_to_bedgraph (requied)
        * bb: dataframe for pre-bigGenePred from operons_to_biggenepred (required)
        * bed: genome BED file (required)
        * fasta: genome fasta file (required)
        * scatterplot_prefix: output filename prefix for scatterplot (default=None)

    Output:
        * scatterplot and dataframe for anomalous ribosome profiles
    """
    
    bb['range'] = bb[['Start','End']].values.tolist()
    bb['Start'] = bb['range'].apply(lambda x: [x[0] + 3*n for n in range(int((x[1]-x[0])/3))])
    bb = bb.explode('Start')
    bb['End'] = bb['Start']+3
    bb = bb[['Chromosome','Start','End','Name','Strand']]
    
    cds = pd.read_csv(bed, sep='\t', header=None)
    cds = cds[cds[6]!=cds[7]].copy()
    cds['Start'] = cds[[1,2]].values.tolist()
    cds['Start'] = cds['Start'].apply(lambda x: [x[0] + 3*n for n in range(int((x[1]-x[0])/3))])
    cds = cds.explode('Start')
    cds['End'] = cds['Start']+3
    cds = cds[[0,'Start','End',3,5]]
    cds.columns = ['Chromosome','Start','End','Name','Strand']
    
    bg_codon = pr.PyRanges(cds).join(pr.PyRanges(bedgraph))
    bg_codon.dna = pr.get_sequence(bg_codon, fasta)
    bg_codon.aa = bg_codon.dna.apply(lambda x: CODON_TO_AA[x])
    
    tbg_codon = pr.PyRanges(bb).join(pr.PyRanges(bg_codon.df.drop(['Start_b', 'End_b'],axis=1)))
    tbg_codon.ORF_type = tbg_codon.Name.str.split('__').str[1]
    tbg_codon = tbg_codon.df[~tbg_codon.df.aa.str.contains('stop')].copy()
    codon = pd.concat([tbg_codon,bg_codon.df])
    codon['ORF_type'] = codon.ORF_type.fillna('mORF')
    codon = codon.drop_duplicates(['Chromosome','Start','End','Strand']).copy()
    ocodon = codon[codon.ORF_type.str.contains('oORF')][['ORF_type','aa','rprofile']].groupby('aa')['rprofile'].agg(['count','min', 'max','median','mean','sem']).reset_index()
    mcodon = codon[codon.ORF_type=='mORF'][['ORF_type','aa','rprofile']].groupby('aa')['rprofile'].agg(['count','min', 'max','median','mean','sem']).reset_index()
    cc = pd.merge(ocodon,mcodon,on='aa')
    
    fig,ax = plt.subplots(figsize=(4,4))
    p = so.Plot(cc,x='mean_y', y='mean_x',text='aa')\
                .add(so.Dot(marker='o'))\
                .add(so.Text(halign='left'))\
                .label(title='oORFs vs mORFs',
                       x='Ribosome profiles by amino acids encoded',
                    y='Ribosome profiles deviating from triplet periodicity')
    plt.errorbar(x='mean_y', y='mean_x', xerr="sem_y", fmt=' ',
                yerr="sem_x", elinewidth=0.5,data=cc, label='aa', capsize=2, capthick=0.5, alpha=0.5)
    
    p.on(ax).save(scatterplot_prefix + '.riboprof_aa.pdf', bbox_inches='tight').show()
    
    if os.path.exists(scatterplot_prefix + '.riboprof_aa.pdf'):
        logging.info('saved scatterplot for anomaly ribosome profiles as ' + scatterplot_prefix + '.riboprof_aa.pdf')
        
    return cc



def riboss(superkingdom, df, riboprof_base, profile, fasta, tx_assembly, 
           genepred=None, ncrna=None, bed=None, 
           orf_range_col=None,
           utr=30, padj_method='fdr_bh', tie=False, num_simulations=1000, 
           run_blastp=False, run_efetch=False, 
           tries=5, sleep=1, 
           email=None, api_key=None, 
           delim=None, outdir=None):
    
    """
    This wrapper is the main RIBOSS function, which 
    - compares the translatability of ORFs within individual transcripts, 
    - runs BLASTP for statistically significant results, 
    - create annotation tracks and metagene plots for novel ORFs, and 
    - compare the ribosome profiles between oORFs and mORFs.
     
    Input:
        * superkingdom: Archaea, Bacteria or Eukaryota (required)
        * df: dataframe from orf_finder/operon_finder (required)
        * riboprof_base: dataframe from parse_ribomap (required)
        * tx_assembly: transcript fasta file extracted using bedtools getfasta. Headers with genomic coordinates (required)
        * fasta: genome fasta file (required)
        * bed: genome BED file (required)
        * utr: padding for metagene plot (default=30)
        * tie: if adjusted p-values between ORFs is not significant (default=False)
        * num_simulations: number of simulations (default=1000)
        * run_blastp: run BLASTP for significant RIBOSS results (default=False)
        * run_efetch: run efetch for top BLASTP hits including RefSeq (default=False)
        * tries: number of attempts (default=5)
        * sleep: pause time seconds (default=1)
        * email: email address (default=None)
        * api_key: NCBI API key https://support.nlm.nih.gov/knowledgebase/article/KA-05317/en-us (default=None)
        * delim: use :: for tx_assembly extracted using bedtools getfasta -name flag, as this appends name (column #4) to the genomic coordinates (default: None)
        * outdir: output directory (default=None)
        
    Output:
        * hits: dataframes for significant RIBOSS results with BLASTP hits
        * BLASTP hits as CSV, JSON and the top hits as pickle, and metagene plots for unannotated ORFs with BLAST hits and no hits as PDFs.
    """
    
    fname = filename(riboprof_base, 'riboss', outdir)

    basename, ext = os.path.splitext(fasta)
    if (ext=='.gz') & (os.path.isfile(fasta)):
        subprocess.run(['gunzip', fasta], check=True)
        fasta = basename
    elif os.path.isfile(basename):
        fasta = basename
    else:
        pass
        
    fai = fasta + '.fai'
    fa = Faidx(fasta)

    if superkingdom in ['Archaea','Bacteria']:
        _, bases = parse_ribomap(superkingdom, riboprof_base, ncrna=ncrna, delim=delim, outdir=outdir)
    elif superkingdom=='Eukaryota':
        gc, _, bases = parse_ribomap(superkingdom, riboprof_base, genepred=genepred, ncrna=ncrna, delim=None, outdir=outdir)
    else:
        logging.error('Please check your spelling! Only Archaea, Bacteria, or Eukaryota is acceptable superkingdrom')
        
    if profile=='Ribosome profiling':
        prof = 0
    elif profile=='RNA-seq':
        prof = 1
    else:
        logging.error('Profile must be either "Ribosome profiling" or "RNA-seq"!')
        
    dt = footprint_counts(df, bases[prof])
    sig, boss_df = boss(dt, tx_assembly, riboprof_base, padj_method=padj_method, tie=tie, num_simulations=num_simulations, outdir=outdir)

    # get bigBed, BED, and GenePred
    if superkingdom in ['Archaea','Bacteria']:
        _ = operons_to_biggenepred(df, sig, bed, fai, fname + '.sig', delim)
        if tie==True:
            _ = operons_to_biggenepred(df, boss_df[boss_df.boss!='mORF'], bed, fai, fname + '.boss_tie', delim)
        else:
            _ = operons_to_biggenepred(df, boss_df[(boss_df.boss!='mORF') & (boss_df.boss!='tie') & (boss_df.boss!='lacks periodicity')], bed, fai, fname + '.boss', delim)
    elif superkingdom=='Eukaryota':
        _ = orfs_to_biggenepred(gc, sig, fai, fname + '.sig', orf_range_col='ORF_range_x', orf_type_col='ORF_type_x')
        if tie==True:
            _ = orfs_to_biggenepred(gc, boss_df[boss_df.boss!='mORF'], fai, fname + '.boss_tie', orf_range_col='ORF_range_x', orf_type_col='ORF_type_x')
        else:
            _ = orfs_to_biggenepred(gc, boss_df[(boss_df.boss!='mORF') & (boss_df.boss!='tie') & (boss_df.boss!='lacks periodicity')], fai, fname + '.boss', orf_range_col='ORF_range_x', orf_type_col='ORF_type_x')
            
    if run_blastp==True:
        if superkingdom in ['Archaea','Bacteria']:
            refseq = 'NP_|WP_'
        elif superkingdom=='Eukaryota':
            refseq = 'NP_|XP_'

        blast = blastp(sig, riboprof_base, email, outdir)
        blast.to_pickle(fname + '.sig.blastp.pkl.gz')
        
        hits = pd.merge(sig, blast.rename(columns={'query':'oid'}), on='oid', how='outer').reset_index(drop=True)
        
        # export cleaner results as pickle
        chits = hits.drop(['boss','odds_ratio','statistical test','result','hit_id','hit_def','e_value','xml'], axis=1)
        rhits = chits.dropna()[chits.dropna().title.str.contains(refseq)]
        tophits = pd.concat([rhits, chits.sort_values('bits', ascending=False)]).drop_duplicates('oid')
        tophits.to_pickle(fname + '.tophits.pkl.gz')
        
        # if superkingdom in ['Archaea','Bacteria']:
        #     bb = operons_to_biggenepred(df, tophits, bed, fai, fname + '.tophits', delim)
        #     # _ = profile_anomaly(bedgraph, bb, bed, fasta, fname)
        # elif superkingdom=='Eukaryota':
        #     bb = orfs_to_biggenepred(gc, tophits, fai, fname + '.tophits', orf_range_col='ORF_range_x', orf_type_col='ORF_type_x')

        # plot top hits
        tb = pd.merge(tophits, bases[prof][['tid','rprofile']])
        tb['start'] = tb['start'] - utr
        tb['end'] = tb['end'] + utr
        tb['ORF_rprofile_x'] = tb[['rprofile','start','end']].values.tolist()
        tb['ORF_rprofile_x'] = tb.ORF_rprofile_x.apply(lambda x: x[0][x[1]:x[2]])
        tb['zeros'] = (np.max(tb.end - tb.start) - tb.ORF_rprofile_x.apply(len))
        tb['ORF_rprofile_x'] = tb[['ORF_rprofile_x','zeros']].values.tolist()
        tb['start_rprofile_x'] = tb['ORF_rprofile_x'].apply(lambda x: x[0] + (x[1]*[0]))
        tb['stop_rprofile_x'] = tb['ORF_rprofile_x'].apply(lambda x: (x[1]*[0]) + x[0])
        blastp_hits = tb[~pd.isnull(tb).any(axis=1)].copy()
        no_hits = tb[pd.isnull(tb).any(axis=1)].copy()

        for ot in tb.ORF_type_x.unique():
            # BLAST hits
            if blastp_hits[blastp_hits.ORF_type_x==ot].shape[0]>0:
                start_rprofile = np.sum(np.array(blastp_hits[blastp_hits.ORF_type_x==ot]['start_rprofile_x'].tolist()), axis=0)
                stop_rprofile = np.sum(np.array(blastp_hits[blastp_hits.ORF_type_x==ot]['stop_rprofile_x'].tolist()), axis=0)
                frames = [0,1,2] * int(start_rprofile.shape[0]/3)
                if start_rprofile.shape[0]-len(frames)==-1:
                    frames = frames[:-1]
                elif start_rprofile.shape[0]-len(frames)==1:
                    frames = frames + [1]
                    
                blastp_rp = pd.DataFrame({'Ribosome profile from start codon':start_rprofile,'Ribosome profile to stop codon':stop_rprofile,'Frames':frames})
                blastp_rp['Position from predicted start codon'] = blastp_rp.index -utr
                blastp_rp['Position to stop codon'] = blastp_rp.index-blastp_rp.index.stop -utr
                predicted_orf_profile(blastp_rp, blastp_hits[blastp_hits.ORF_type_x==ot], utr, str(ot) + 's with BLASTP hits', fname)
            else:
                pass
                
            # no BLAST hits
            if no_hits[no_hits.ORF_type_x==ot].shape[0]>0:
                start_rprofile = np.sum(np.array(no_hits[no_hits.ORF_type_x==ot]['start_rprofile_x'].tolist()), axis=0)
                stop_rprofile = np.sum(np.array(no_hits[no_hits.ORF_type_x==ot]['stop_rprofile_x'].tolist()), axis=0)
                frames = [0,1,2] * int(start_rprofile.shape[0]/3)
                if start_rprofile.shape[0]-len(frames)==-1:
                    frames = frames[:-1]
                elif start_rprofile.shape[0]-len(frames)==1:
                    frames = frames + [1]
                    
                no_rp = pd.DataFrame({'Ribosome profile from start codon':start_rprofile,'Ribosome profile to stop codon':stop_rprofile,'Frames':frames})
                no_rp['Position from predicted start codon'] = no_rp.index -utr
                no_rp['Position to stop codon'] = no_rp.index-no_rp.index.stop -utr
                predicted_orf_profile(no_rp, no_hits[no_hits.ORF_type_x==ot], utr, str(ot) + 's with no BLASTP hits', fname)
            else:
                pass
    
        hits.drop(['start', 'end'], axis=1, inplace=True)
        hits.to_csv(fname + '.sig.blastp.csv', index=None)
        hits.to_json(fname + '.sig.blastp.json', index=None)
        
        logging.info('saved BLASTP results for RIBOSS hits as ' + fname + '.tophits.pkl.gz, ' + fname + '.sig.blastp.csv, ' + fname + '.sig.blastp.json, and ' + fname + '.sig.blastp.pkl.gz')

        
        if (run_blastp==True) & (run_efetch==True):
            w = blast.dropna()[blast.dropna().accession.str.contains(refseq)].accession.unique()
            t = tophits.accession.dropna().unique()
            acc = np.unique(np.concatenate([w, t]))
            ipg = efetch(acc, tries, sleep, email, api_key)
            ipg.to_pickle(fname + '.sig.ipg.pkl.gz')
            
            return ipg, tophits, blast, sig, boss_df
            
        elif (run_blastp==True) & (run_efetch==False): 
            logging.info('saved BLASTP results for RIBOSS hits as ' + fname + '.tophits.pkl.gz, ' + fname + '.sig.blastp.csv, ' + fname + '.sig.blastp.json, and ' + fname + '.sig.blastp.pkl.gz')
            return tophits, blast, sig, boss_df
            
        elif (run_blastp==False) & (run_efetch==True):
             logging.error('Please enable BLASTP!')
            
    else:
        return sig, boss_df
