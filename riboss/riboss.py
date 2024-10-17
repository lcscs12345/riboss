#!/usr/bin/env python
# coding: utf-8

"""
@author      CS Lim
@create date 2020-09-15 17:40:16
@modify date 2024-10-17 18:27:07
@desc        Main RIBOSS module
"""




import os, re, sys, time, argparse, textwrap, csv, logging, subprocess
from urllib.request import HTTPError
import numpy as np
import pandas as pd
import seaborn as sns
import seaborn.objects as so
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import scipy.stats as stats
from scipy.stats.contingency import odds_ratio
from types import SimpleNamespace
from collections import namedtuple
from tqdm import tqdm
import pyranges as pr
from riboss.orfs import translate
from Bio.Blast import NCBIWWW
from Bio import Entrez as ez
from io import StringIO
from .chisquare import chi_square_posthoc
from .orfs import orf_finder, CODON_TO_AA
from .wrapper import filename, merge_scores
from .read import read_blast_xml
 


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




def base_to_bedgraph(df, bedgraph_prefix, delim=None, outdir=None):
    """
    Convert a dataframe from parse_ribomap to BedGraph format.

    Input:
        * df: dataframe from parse_ribomap (required)
        * bedgraph_prefix: output filename prefix for BedGraph (required)
        * delim: use :: for tx_assembly extracted using bedtools getfasta -name flag, as this appends name (column #4) to the genomic coordinates (default: None)
        * outdir: output directory (default=None)
    Output:
        * BedGraph file
    """
    
    if delim!=None:
        pos = 1
    else:
        pos = 0

    fname = filename(bedgraph_prefix, 'riboprof', outdir)
        
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
    df = df.drop(['tid','Start','End'],axis=1).explode(['range','rprofile'])
    df['End'] = df['range'] +1
    df.rename(columns={'range':'Start'}, inplace=True)
    df['rprofile'] = df.rprofile.astype(int)
    df = df[df.rprofile!=0].copy()

    # Split by strands
    plus = df[df.Strand=='+'][['Chromosome','Start','End','rprofile']]    
    minus = df[df.Strand=='-'][['Chromosome','Start','End','rprofile']]
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
        df[['Chromosome','Start','End','rprofile']].to_csv(fname + '.bg', index=None, header=None, sep='\t')
        bg = merge_scores(fname + '.bg')
        f = open(fname + '.bg', 'w')
        f.write('track type=bedGraph name="riboprof" description="Ribosome profile" visibility=full color=0,0,0 priority=20\n')
        bg.to_csv(f, index=None, header=None, sep='\t', mode='a')
        f.close()
        
        if bg[bg[0].str.contains('ERROR')].shape[0]!=0:
            logging.error('Failed generating BedGraph!')
        else:
            logging.info('saved ribosome profiles as ' + fname + '.bg')
    
    bg.columns = ['Chromosome','Start','End','rprofile']
    
    return bg



def parse_ribomap(superkingdom, base, delim=None, outdir=None):
    """
    Parse a ribomap/riboprof base file into a dataframe
    
    Input:
        * superkingdom: Archaea, Bacteria or Eukaryota (required)
        * base: output base file from riboprof (required)
        * delim: use :: for tx_assembly extracted using bedtools getfasta -name flag, as this appends name (column #4) to the genomic coordinates (default: None)
        * outdir: output directory (default=None)
    
    Output:
        * a dataframe
    """
    
    dt = pd.read_csv(base, sep=r':\s', engine='python', header=None)
    dt = dt[(dt[0]=='tid') | (dt[0]=='ribo profile')]
    dt.columns = ['id','val']
    dt = pd.DataFrame({'tid':dt['val'].iloc[::2].values, 'rprofile':dt['val'].iloc[1::2].values})

    pbar = tqdm.pandas(desc="parsing ribomap output ", unit_scale = True)
    dt['rprofile'] = dt.rprofile.str.split().progress_apply(lambda x: [float(i) for i in x])
    
    if superkingdom in ['Archaea','Bacteria']:
        bg = base_to_bedgraph(dt, base, delim=delim, outdir=outdir)
    # TO DO
    # elif superkingdom=='Eukaryota':
    
    return bg, dt



def footprint_counts(superkingdom, df, d, gencode=None):
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

    pbar = tqdm.pandas(desc="counting footprints    ", unit_scale = True)
    dt['periodicity'] = dt[['rprofile','frame0','frame1','frame2']].values.tolist()
    dt['periodicity'] = dt.periodicity.progress_apply(lambda x: [sum([x[0][i] for i in tuple(x[1])]),
                                                            sum([x[0][i] for i in tuple(x[2])]),
                                                            sum([x[0][i] for i in tuple(x[3])])])
    
    # get rprofile at start codons
    dt['start_pos'] = dt.ORF_range.apply(lambda x: [x[0], x[0]+1, x[0]+2])
    dt['start_pos'] = dt[['rprofile','start_pos']].values.tolist()
    dt['start_rprofile'] = dt.start_pos.apply(lambda x: [x[0][i] for i in x[1]])

    # filter rprofile by frames
    dt['frame0'] = dt.periodicity.apply(lambda x: x[0])
    dt['frame1'] = dt.periodicity.apply(lambda x: x[1])
    dt['frame2'] = dt.periodicity.apply(lambda x: x[2])    
    dt = dt[(dt.frame0 + dt.frame1 + dt.frame2>10) & (dt.frame0>dt.frame1) & (dt.frame0>dt.frame2)].reset_index(drop=True)
    # prioritise start_rprofile>0
    d1 = dt[dt.start_rprofile.apply(lambda x: x[0])>0].sort_values(['start_codon','ORF_start','tid'])
    d0 = dt[dt.start_rprofile.apply(lambda x: x[0])==0].sort_values(['start_codon','ORF_start','tid'])
    dt = pd.concat([d1, d0])
    dt['tab'] = dt[['frame0','frame1','frame2']].values.tolist()

    if superkingdom in ['Archaea','Bacteria']:
        dt.drop_duplicates(['tid','ORF_type','ORF_end'], inplace=True)
        dt.drop(['rprofile','periodicity','start_pos','frame0','frame1','frame2'], axis=1, inplace=True)
    elif superkingdom=='Eukaryota':
        dt.drop_duplicates(['tid','ORF_type','overlap','ORF_end'], inplace=True)
        dt.drop(['rprofile','genomic_start','periodicity','start_pos','frame0','frame1','frame2'], axis=1, inplace=True)
    else:
        sys.exit('Superkingdom is required! Choose either Archaea, Bacteria or Eukaryota.')    
    
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
    


def statistical_test(f, num_simulations=1000):
    
    # calculate odds ratios
    f['odds_ratio'] = f['tab'].apply(lambda x: odds_ratio(x[1]))
    
    # perform (chi-square test, num_simulations bootstrap, post-hoc test) or (boschloo-exact test)
    pbar = tqdm.pandas(desc="comparing periodicity  ", unit_scale = True)
    
    f['result'] = f.tab.progress_apply(lambda x: chi_square_posthoc(x[0],num_simulations) if np.sum(x[0], axis=0).all()==True else stats.boschloo_exact(x[1])) # # ((x[0][0][1]+x[0][0][2]>0) & x[0][1][1]+x[0][1][2]>0)
    f['statistical test'] = f['result'].apply(lambda x: 'ChiSquare' if type(x)==SimpleNamespace else 'BoschlooExact')
    
    # bonferroni correction
    f['adj_pval'] = f['result'].apply(lambda x : 1 if x.pvalue*len(f) > 1 else (x.pvalue*len(f) if x.pvalue*len(f) < 0.05 else x.pvalue))
    
    chi = f[f['statistical test']=='ChiSquare'].reset_index(drop=True)
    chi['statistic'] = chi['result'].apply(lambda x: x.statistic)
    chi['posthoc_statistic'] = chi['result'].apply(lambda x: x.posthoc_statistic)
    chi['adj_bootstrap_pval'] = chi['result'].apply(lambda x : 1 if x.bootstrap_pvalue*len(chi) > 1 else (x.bootstrap_pvalue*len(chi) if x.bootstrap_pvalue*len(chi) < 0.05 else x.bootstrap_pvalue))
    chi['adj_posthoc_pval'] = chi['result'].apply(lambda x : np.array([1 if i*(len(chi)+12) > 1 else (i*(len(chi)+12) if i*(len(chi)+12) < 0.05 else i) for i in x.posthoc_pvalue]))
    
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

    return f



def boss(superkingdom, df, tx_assembly, boss_prefix, tie=False, num_simulations=1000, outdir=None):  
    """
    Compare uORFs, oORFs and dORFs to mORFs, and assign which ORFs are the bosses.
    
    Input:
        * superkingdom: Archaea, Bacteria or Eukaryota (required)
        * df: dataframe from footprint_counts (required)
        * tx_assembly: transcript fasta file extracted using bedtools getfasta. Headers with genomic coordinates (required)
        * tie: if adjusted p-values between ORFs is not significant (default=False)
        * num_simulations: number of simulations (default=1000)
        * outdir: output directory (default=None)
    
    Output:
        * sig: dataframe for significant RIBOSS results
        * boss: dataframe for all RIBOSS statistics
        * RIBOSS statistics and significant results as CSV/JSON/BED 
    """

    m = df[df.ORF_type=='mORF'].reset_index(drop=True)
    o = df[df.ORF_type!='mORF'].reset_index(drop=True)
    f = pd.merge(o, m, on='tid', how='outer')
    noo = f[f.start_codon_x.isna()].reset_index(drop=True)
    noo['boss'] = noo['ORF_type_y']
    nom = f[f.start_codon_y.isna()].reset_index(drop=True)
    nom['boss'] = nom['ORF_type_x']
    f = pd.merge(o, m, on='tid')
    f = contingency_tables(f)
    f = statistical_test(f, num_simulations)    

    # Are mORFs the boss?
    f.reset_index(drop=True, inplace=True)
    f = f.reset_index()
    boss_x = f[(f.odds_ratio.apply(lambda x: x.statistic>1)) & (f['result'].apply(lambda x: x.adjusted_pvalue<0.05)) & (f['result'].apply(lambda x: (x.adjusted_bootstrap_pvalue<0.05) if 'adjusted_bootstrap_pvalue' in x._fields else x.adjusted_pvalue<0.05))].reset_index(drop=True)
    boss_x['boss'] = boss_x['ORF_type_x']
    boss_y = f[(f.odds_ratio.apply(lambda x: x.statistic<1)) & (f['result'].apply(lambda x: x.adjusted_pvalue<0.05)) & (f['result'].apply(lambda x: (x.adjusted_bootstrap_pvalue<0.05) if 'adjusted_bootstrap_pvalue' in x._fields else x.adjusted_pvalue<0.05))].reset_index(drop=True)
    boss_y['boss'] = boss_y['ORF_type_y']
    boss = pd.concat([boss_x,boss_y,f]).drop_duplicates('index')
    boss['boss'] = boss['boss'].fillna('tie')
    boss = pd.concat([boss,noo,nom])

    if superkingdom in ['Archaea','Bacteria']:
        boss = boss[['tid', 'boss', 'start_codon_x', 'start_codon_y', 'ORF_range_x', 'ORF_type_x', 'ORF_range_y', 'ORF_type_y', 'tab', 'odds_ratio', 'statistical test', 'result']]
    elif superkingdom=='Eukaryota':
        boss = boss[['tid', 'boss', 'start_codon_x', 'start_codon_y', 'ORF_range_x', 'ORF_type_x', 'ORF_range_y', 'ORF_type_y', 'overlap_x', 'overlap_y', 'tab', 'odds_ratio', 'statistical test', 'result']]
    else:
        sys.exit('Superkingdom is required! Choose either Archaea, Bacteria or Eukaryota.')
        
    fname = filename(boss_prefix, 'riboss', outdir)
    
    boss.dropna(inplace=True)
    boss.reset_index(drop=True, inplace=True)
    boss['odds_ratio'] = boss.odds_ratio.apply(lambda x: x.statistic)    
    boss.to_csv(fname + '.csv', index=None)
    boss.to_json(fname + '.json', index=None)
    logging.info('saved RIBOSS stats as ' + fname + '.csv and ' + fname + '.json')
    
    # extract significant results        
    chi = boss[(boss['statistical test']=='ChiSquare')].reset_index(drop=True)
    be = boss[(boss['statistical test']=='BoschlooExact')].reset_index(drop=True)
    
    if tie==True:
        chi = chi[(chi.result.apply(lambda x: x.adjusted_posthoc_pvalue[0])<0.05) & (chi.boss!='mORF')].reset_index(drop=True)
        be = be[be.result.apply(lambda x: x.adjusted_pvalue<0.05)].reset_index(drop=True)
    else:
        chi = chi[(chi.result.apply(lambda x: x.adjusted_posthoc_pvalue[0])<0.05) & (chi.boss!='tie') & (chi.boss!='mORF')].reset_index(drop=True)
        be = be[be.result.apply(lambda x: x.adjusted_pvalue<0.05) & (be.boss!='tie')].reset_index(drop=True)
        
    sig = pd.concat([chi,be])
    sig['start'] = sig.ORF_range_x.apply(lambda x: x[0])
    sig['end'] = sig.ORF_range_x.apply(lambda x: x[1])
    sig.drop_duplicates(['tid','start','end'], inplace=True)

    sig.rename(columns={'tid':'Chromosome','start':'Start','end':'End'}, inplace=True)
    sig = pr.PyRanges(sig)
    sig.dna = pr.get_sequence(sig, tx_assembly)
    sig.aa = sig.dna.apply(lambda x: translate(x))
    sig.oid = sig.Chromosome.astype(str) + '__' + sig.Start.astype(str) + '-' + sig.End.astype(str)
    sig.fa = ('>' + sig.oid + '\n' + sig.aa).tolist()
    sig = sig.df

    sig.rename(columns={'Chromosome':'tid','Start':'start','End':'end'}, inplace=True)
    sig.reset_index(drop=True, inplace=True)
    sig[['tid','start','end']].to_csv(fname + '.sig.bed', sep='\t', header=None, index=None)
    sig.to_csv(fname + '.sig.csv', index=None)
    sig.to_json(fname + '.sig.json', index=None)
    logging.info('saved significant RIBOSS results (n=' + str(sig.shape[0]) + ') as ' + fname + '.sig.csv, ' + fname + '.sig.json, and ' + fname + '.sig.bed')

    return sig, boss



def blastp(sig, blastp_prefix, email=None, outdir=None):
    """
    BLASTP for ORFs with significantly greater triplet periodicity than mORFs

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
    Plot metagene ribosome profiles for unannotated ORFs
    
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
        
    plt.set_loglevel('WARNING')
    plt.figure(figsize=(4,3))
    ax = sns.barplot(data=rp, x='Position from predicted start codon', y='Ribosome profile from start codon', hue='Frames', width=1)
    ax.xaxis.set_major_locator(ticker.MultipleLocator(3))
    plt.title(title + ' (n=' + str(df.shape[0]) + ', median peptide length=' + str(peptide_len) + ')')
    plt.ylabel('Metagene ribosome profile')
    plt.xlim(-12+utr, 30+utr)
    plt.savefig(barplot_prefix + '.' + infix + '.start_codon.pdf', bbox_inches='tight')
            
    logging.info('saved metagene plots as ' + barplot_prefix + '.' + infix + '.start_codon.pdf')
    


def tophits_to_biggenepred(orf, tophits, bed, fai, big_fname, delim=None):
    """
    Create UCSC bigGenePred for RIBOSS top hits

    Input:
        * orf: dataframe from operon_finder (required)
        * tophits: dataframe for RIBOSS top hits (required)
        * bed: genome BED file (required)
        * fai: genome fasta index to generate chrom.sizes (required)
        * big_fname: output filename prefix for bigGenePred
        * delim: use :: for tx_assembly extracted using bedtools getfasta -name flag, as this appends name (column #4) to the genomic coordinates (default: None)
        
    Output:
        * bigGenePred for tophits and dataframes for both tophits and all ORFs 
    """
    
    orf = orf.rename(columns={'start_codon':'start_codon_x','ORF_start':'start','ORF_end':'end'}).drop('ORF_length', axis=1)
    toporf = pd.merge(orf, tophits)
    toporf['name2'] = toporf.Chromosome.astype(str) + ':' + toporf['Start'].astype(str) + '-' + toporf['End'].astype(str) + '(' + toporf.Strand.astype(str) + ')' 
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
        bed = big_fname + '.' + orftype + '.tophits.bed'
        bb[bb.ORF_type_x==orftype].drop('ORF_type_x',axis=1).to_csv(bed, index=None, header=None, sep='\t')
        
        pre = big_fname + '.' + orftype + '.tophits.sorted.bed'
        subprocess.run(['bedSort', bed, pre], check=True)
            
        bgp = big_fname + '.' + orftype + '.tophits.bb'
        subprocess.run(['bedToBigBed','-as=' + bas,
                        '-type=bed12+8', pre, chromsizes, bgp], check=True)
    
        os.remove(bed)        
        os.remove(pre)
        
        if os.path.exists(bgp):
            logging.info('saved bigGenePred for RIBOSS top hits as ' + bgp)
        
    os.remove(chromsizes)
    os.remove(bas)
    
    return bb



def profile_anomaly(bedgraph, bb, bed, fasta, scatterplot_prefix=None):
    """
    Find anomalous ribosome profiles by the encoded amino acids
    
    Input:
        * bedgraph: dataframe for BedGraph from base_to_bedgraph (requied)
        * bb: dataframe for pre-bigGenePred from tophits_to_biggenepred (required)
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
    p = so.Plot(cc,x='mean_x', y='mean_y',text='aa')\
                .add(so.Dot(marker='o'))\
                .add(so.Text(halign='left'))\
                .label(title='Ribosome profiles by the encoded amino acids',
                       x='Codons deviating from triplet periodicity', y='Other codons')
    plt.errorbar(x='mean_x', y='mean_y', xerr="sem_x", fmt=' ',
                  yerr="sem_y", elinewidth=0.5,data=cc, label='aa', capsize=2, capthick=0.5)
    
    p.on(ax).show()
    plt.savefig(scatterplot_prefix + '.riboprof_aa.pdf', bbox_inches='tight')
    
    if os.path.exists(scatterplot_prefix + '.riboprof_aa.pdf'):
        logging.info('saved scatterplot for anomaly ribosome profiles as ' + scatterplot_prefix + '.riboprof_aa.pdf')
        
    return cc



def riboss(superkingdom, df, orf, riboprof_base, tx_assembly, fasta, fai, bed, 
           utr=30, tie=False, num_simulations=1000, 
           run_blastp=False, run_efetch=False, 
           tries=5, sleep=1, 
           email=None, api_key=None, 
           delim=None, outdir=None):
    
    """
    A wrapper for the functions above and construct metagene plots for unannotated ORFs.
     
    Input:
        * superkingdom: Archaea, Bacteria or Eukaryota (required)
        * df: dataframe from orf_finder/operon_finder (required)
        * riboprof_base: dataframe from parse_ribomap (required)
        * tx_assembly: transcript fasta file extracted using bedtools getfasta. Headers with genomic coordinates (required)
        * fasta: genome fasta file (required)
        * bed: genome BED file (required)
        * fai: genome fasta index (required)
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
    
    bedgraph, base = parse_ribomap(superkingdom, riboprof_base, delim=delim, outdir=outdir)

    dt = footprint_counts(superkingdom, df, base)
    sig, boss_df = boss(superkingdom, dt, tx_assembly, riboprof_base, tie, num_simulations, outdir)

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
        bb = tophits_to_biggenepred(orf, tophits, bed, fai, fname, delim)
        _ = profile_anomaly(bedgraph, bb, bed, fasta, fname)
        
        # plot top hits
        tb = pd.merge(tophits, base[['tid','rprofile']])
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
            predicted_orf_profile(blastp_rp, blastp_hits[blastp_hits.ORF_type_x==ot], utr, str(ot) + ' with BLASTP hits', fname)
            
            # no BLAST hits
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
            predicted_orf_profile(no_rp, no_hits[no_hits.ORF_type_x==ot], utr, str(ot) + ' with no BLASTP hits', fname)
    
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