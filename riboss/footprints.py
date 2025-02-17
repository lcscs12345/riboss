#!/usr/bin/env python
# coding: utf-8

"""
@author      CS Lim
@create date 2024-09-14 10:42:41
@modify date 2025-02-17 19:16:01
@desc        RIBOSS module for analysing aligned ribosome footprints
"""




import os, time, argparse, textwrap, csv, random, pysam, sys, logging.config, re
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from tqdm import tqdm
from .riboss import contingency_tables, statistical_test
from .wrapper import filename


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


# def check_downsampling(fraction):        
#     if 0.1 <= float(fraction) <= 1:
#         return float(fraction)
#     else:
#         raise argparse.ArgumentTypeError('Fraction out of range.')

        
# def check_size(size):        
#     if 25 <= int(size) <= 35:
#         return int(size)
#     else:
#         raise argparse.ArgumentTypeError('Footprint size cutoff out of range. Recommend 25 to 35')


# def io(infile, outdir):
#     check_dir = os.path.isdir(outdir)
#     if not check_dir:
#         os.makedirs(outdir)

#     try:
#         base = os.path.basename(infile) # .split(os.extsep)
#         if re.search('.bam$', infile):
#             offset = outdir + os.path.splitext(base)[0] + '_offset.txt'
#             heatmap = outdir + os.path.splitext(base)[0] + '_heatmap.pdf'
#             barplot = outdir + os.path.splitext(base)[0]
#             return offset, heatmap, barplot

#         elif re.search('.fasta$|.fas$|.fa$', infile):
#             barplot = outdir + os.path.splitext(base)[0]
#             return barplot
            
#     except Exception:
#         raise argparse.ArgumentTypeError('Please provide a bam or fasta file as input!')
    



def bam_sampling(bamfile, fraction=0.1):
    """
    Downsampling a BAM file. As input footprint_summary.

    Input:
        * bamfile: BAM file path for ribosome footprint BAM (required)
        * fraction: fraction of reads to sample (default: 0.1)
    Output:
        * sample: downsampled BAM file
    """
    
    infile = pysam.AlignmentFile(bamfile, 'rb')
    sample = []
    random.seed(123)

    with tqdm(desc='downsampling BAM       ') as pbar:
        for read in infile:
            rand = random.random()
            pbar.update()
            if read.mapping_quality==255 and rand < float(fraction):
                sample.append([str(read.reference_name),str(read).split('\t')])               
    infile.close()
    
    return sample



def footprint_summary(cds_range, sample, quality='best', offset_method='5p', adj=12, min_size=None, max_size=None, offset_prefix=None):
    """
    Compare the periodicity of ribosome footprint by sizes using odds ratio and chi-square post-hoc test.
    All vs the most abundant footprint size.
    
    Input (all required):
        * offset_method: 5p or 3p (required)
        * adj: offset value for the footprint size of 28 nt (required)
        * cds_range: dataframe or tsv file from orf_finder or operon_finder (orfs.py) (required)
        * sample: dataframe from bam_sampling (footprints.py) (required)
        * quality: quality of selected footprints. Option: best, good, or fair (required)
        * min_size: minimum footprint size for analysis (required)
        * max_size: maximum footprint size for analysis (required)
        * offset_prefix: filename prefix for footprint length(s) with offset
    Output:
        * stats: dataframe with p-values for odds ratio and chi-square post-hoc test
        * frame_stats: dataframe for heatmaps
        * selected_footprints: dataframe for footprint_periodicity
    """
    
    if isinstance(cds_range, str):
        cds = pd.read_csv(cds_range, sep='\t', header=None)
    elif isinstance(cds_range, pd.DataFrame):
        cds = cds_range
    
    cds.columns = ['tid','CDS_start','CDS_end']
        
    df = pd.DataFrame(sample)
    df.columns = ['tid','read']
    pos = pd.DataFrame(df['read'].apply(lambda x: [x[3],x[5]]).tolist(), 
                       columns = ['pos','cigar'])
    df = pd.concat([df['tid'],pos], axis=1)
    df['footprint_len'] = df.cigar.str.split('M', expand=True)[[0]]
    df = df[~df.footprint_len.str.contains('S')]
    df['pos'] = df.pos.astype(int) -1
    df['footprint_len'] = df.footprint_len.astype(int)
    
    if min_size==None:
        min_size = np.min(df.footprint_len.tolist())
        logging.info('use ' + str(min_size) + ' as minimum footprint length.')
        
    if max_size==None:
        max_size = np.max(df.footprint_len.tolist())
        if max_size<40:
            logging.info('use ' + str(max_size) + ' as the maximum footprint length.')
        elif max_size>=50:
            logging.warning('Use ' + str(max_size) + ' as the maximum footprint length!')
            
        df = df[(df.footprint_len>=int(min_size)) & (df.footprint_len<=int(max_size))]
        
    else:
        df = df[(df.footprint_len>=int(25)) & (df.footprint_len<=int(35))]


    dt = pd.merge(cds,df,on='tid')
    adj = int(adj)
    
    if offset_method=='5p':     
        dt['adj_start'] = dt.pos - dt.CDS_start + adj
        dt['adj_end'] = dt.pos - dt.CDS_end + adj
        dt['frame'] = dt['adj_start']%3        
    elif offset_method=='3p':
        dt['adj_start'] = dt.pos - dt.CDS_start + dt.footprint_len - adj
        dt['adj_end'] = dt.pos - dt.CDS_end + dt.footprint_len - adj
        dt['frame'] = dt['adj_start']%3
    else:
        sys.exit('Offset method is required! Choose either 5p, 3p or centre.')
            
    frame_stats = dt[(dt.frame>=0) & (dt['adj_start']<=90) & (dt['adj_start']>=-30)].copy()
    frame_stats = frame_stats.groupby(['footprint_len','frame']).count().reset_index()[['footprint_len','frame','tid']]
    frame_stats.columns = ['footprint_len','frame','counts']
    
    total_fp = np.sum(frame_stats['counts'].tolist())
    if total_fp<100:
        logging.warning('Total fragment counts are ' + str(total_fp) + '! Metagene plots might not be generated.')
    
    pv = frame_stats.pivot(columns='frame', index='footprint_len', values='counts').sort_values('footprint_len', ascending=False)
    pv['all counts'] = pv[[0,1,2]].values.tolist()
    pv['% counts by footprints'] = pv['all counts'].apply(lambda x: [i/np.sum(x) for i in x])
    pv['0'] = pv['% counts by footprints'].apply(lambda x: x[0])
    pv['1'] = pv['% counts by footprints'].apply(lambda x: x[1])
    pv['2'] = pv['% counts by footprints'].apply(lambda x: x[2])
    
    pv['% counts'] = pv['all counts'].apply(lambda x: [i/np.sum(pv[[0,1,2]].sum()) for i in x])
    pv[0] = pv['% counts'].apply(lambda x: x[0])
    pv[1] = pv['% counts'].apply(lambda x: x[1])
    pv[2] = pv['% counts'].apply(lambda x: x[2])
    pv.dropna(inplace=True)
    pv['frame_max'] = pv['all counts'].apply(lambda x: np.argmax(x))
    
    fplen = pd.DataFrame({'footprint_len':range(28-6,28+73), 'offset':range(adj-6,adj+73)})
    fplen['frame_max'] = fplen.offset%3
    pv = pd.merge(pv, fplen, on='footprint_len')
    
    f0x = pv.value_counts('frame_max_x').reset_index().sort_values('frame_max_x').iloc[0]['count']
    f0y = pv.value_counts('frame_max_y').reset_index().sort_values('frame_max_y').iloc[0]['count']
    if f0x>f0y:
        pv.offset = adj
    elif f0x==f0y:
        pv = pv[pv.frame_max_x==pv.frame_max_y].copy()
    else:
        pv = pv[pv.frame_max_x==pv.frame_max_y].copy()
        logging.warning('The default adj=' + str(adj) + ' for footprint offset is not suitable! Try a different offset value, for example adj=13 or adj=11.')
        
    try:
        best = frame_stats.sort_values('counts', ascending=False).head(1).footprint_len.values[0]
        if best>35:
            logging.warning('The most abundant footprint size is ' + str(best) + ' and longer than expected!')
        
        footprint_stats = frame_stats.groupby('footprint_len')['counts'].apply(list).reset_index()
        top_footprint = footprint_stats[footprint_stats.footprint_len==best]
        top_footprint = top_footprint.loc[top_footprint.index.repeat(footprint_stats.shape[0])].reset_index(drop=True)
        footprint_stats = pd.concat([top_footprint,footprint_stats], axis=1)
        footprint_stats.columns = ['top_footprint','tab_x','footprint_len','tab_y']
        footprint_stats = footprint_stats[(footprint_stats.tab_x.apply(len)==3) & (footprint_stats.tab_y.apply(len)==3)].copy()
        footprint_stats = contingency_tables(footprint_stats)

        f,xf = statistical_test(footprint_stats, 1000)
        
        if quality=='best':
            best = pd.DataFrame({'footprint_len':[best]})
            fplen = pd.merge(best, pv[['footprint_len','offset']])
        elif quality=='good':
            if f.shape[0]>0:
                good = f[(f.adj_posthoc_pval.apply(lambda x: x[0]>0.01)) | (f.odds_ratio.apply(lambda x: x.statistic)<2)].copy()[['footprint_len']]
            else:
                xf['odds_ratio'] = xf.odds_ratio.apply(lambda x: x.statistic)
                good = xf.sort_values('odds_ratio', ascending=False)[['footprint_len']]    
            fplen = pd.merge(good, pv[['footprint_len','offset']])
        else:
            sys.exit('Please provide footprint quality!')
        
        fplen.to_csv(offset_prefix + '.offset.txt', sep='\t', index=None, header=None)
        logging.info('saved selected footprint sizes with an offset as ' + offset_prefix + '.offset.txt')
        
        selected_footprints = pd.merge(dt,fplen)
    
        if offset_method=='5p':
            selected_footprints['adj_start'] = selected_footprints.pos - selected_footprints.CDS_start + selected_footprints.offset
            selected_footprints['adj_end'] = selected_footprints.pos - selected_footprints.CDS_end + selected_footprints.offset
        
        elif offset_method=='3p':
            selected_footprints['adj_start'] = selected_footprints.pos - selected_footprints.CDS_start + selected_footprints.footprint_len - selected_footprints.offset
            selected_footprints['adj_end'] = selected_footprints.pos - selected_footprints.CDS_end + selected_footprints.footprint_len - selected_footprints.offset
    
        pv.drop(['frame_max_x','frame_max_y'], axis=1, inplace=True)
        pv.set_index('footprint_len', inplace=True)
        
        return f, pv, selected_footprints
    
    except Exception as e:
        logging.exception('Reads are longer than expected! Are you sure that this is a ribosome profiling library? Here are the statistics by reading frame.')
        print(frame_stats)
        pass


def heatmap_periodicity(pv, heatmap_prefix):
    """
    Plot heatmaps by frames. Using frame_stats dataframe from footprint_summary.
    
    Input:
        * frame_stats: dataframe from footprint_summary (required)
        * heatmap_prefix: filename prefix for heatmaps (required)
    Output:
        * heatmaps as PDFs
    """
    
    if re.search('Aligned.out',heatmap_prefix):
        bam = heatmap_prefix.replace('Aligned.out','')
        bam = os.path.basename(bam)
    else:
        bam = heatmap_prefix
        bam = os.path.basename(bam)
        
    sns.set(style='ticks',font_scale=1)
    fig, axes = plt.subplots(1, 2, figsize=(4, 3.5), constrained_layout=True)
    
    plt.subplot(1, 2, 1)
    if np.max(pv[0])>0.75:
        g1 = sns.heatmap(pv[[0,1,2]], cmap="viridis", vmin=0, vmax=1, cbar_kws = {'ticks': [0,0.5,1],'location':'top'})
        g1.set(xlabel=None, ylabel=None)
    elif (np.max(pv[0])>0.5) & (np.max(pv[0])<=0.75):
        g1 = sns.heatmap(pv[[0,1,2]], cmap="viridis", vmin=0, vmax=0.75, cbar_kws =  {'ticks': [0,0.25,0.5,0.75],'location':'top'})
        g1.set(xlabel=None, ylabel=None)
    elif (np.max(pv[0])>0.25) & (np.max(pv[0])<=0.5):
        g1 = sns.heatmap(pv[[0,1,2]], cmap="viridis", vmin=0, vmax=0.5, cbar_kws =  {'ticks': [0,0.25,0.5],'location':'top'})
        g1.set(xlabel=None, ylabel=None)
    elif np.max(pv[0])<=0.25:
        g1 = sns.heatmap(pv[[0,1,2]], cmap="viridis", vmin=0, vmax=0.25, cbar_kws =  {'ticks': [0,0.25],'location':'top'})
        g1.set(xlabel=None, ylabel=None)
        
    plt.subplot(1, 2, 2)
    g2 = sns.heatmap(pv[['0','1','2']], cmap="coolwarm", vmin=0, vmax=1, cbar_kws = {'ticks': [0,0.5,1],'location':'top'})
    g2.set(xlabel=None, ylabel=None)
    
    fig.suptitle(bam, y=1.1)
    fig.supxlabel('Footprint positions on reading frames', fontsize='medium')
    fig.supylabel('Footprint size', fontsize='medium')

    plt.savefig(heatmap_prefix + '.frames.pdf', bbox_inches='tight')
    logging.info('converted mapped frames into heatmaps as ' + heatmap_prefix + '.frames.pdf')



def footprint_periodicity(selected_footprints, downsampling, barplot_prefix, ylim=None):
    """
    Plot metagene barplots by footprint sizes.
    
    Input:
        * selected_footprints: dataframe from footprint_summary (required)
        * downsampling: fraction used in bam_sampling (required)
        * barplot_prefix: filename prefix for barplots (required)
        * ylim: option to use the same ylim for all barplots (default: None)
    Output:
        * metagene plots as PDFs
    """

    if ylim!=None:
        # Roundup to next multiple of 100
        # https://stackoverflow.com/a/14092788
        ymax = selected_footprints[(selected_footprints.adj_start<=21) & (selected_footprints.adj_start>=-21)].value_counts(['footprint_len','adj_start']).iloc[0]
        ymax -= ymax % -100
        ylim=(0,ymax)

    selected_footprints = selected_footprints.sort_values('footprint_len').copy()
    selected_footprints['footprint_len'] = selected_footprints['footprint_len'].astype(str) + '-nt with offset=' + selected_footprints['offset'].astype(str)
    
    fl_start = set(selected_footprints[(selected_footprints.adj_start<=21) & (selected_footprints.adj_start>=-21)].footprint_len.unique())
    fl_end = set(selected_footprints[(selected_footprints.adj_end<=21) & (selected_footprints.adj_end>=-21)].footprint_len.unique())
    fl = pd.DataFrame({'footprint_len':list(fl_start.intersection(fl_end))})
    
    df_start = pd.merge(selected_footprints[(selected_footprints.adj_start<=21) & (selected_footprints.adj_start>=-21)],fl)
    df_end = pd.merge(selected_footprints[(selected_footprints.adj_end<=21) & (selected_footprints.adj_end>=-21)],fl)

    if re.search('Aligned.out',barplot_prefix):
        bam = barplot_prefix.replace('Aligned.out','')
        bam = os.path.basename(bam)
    else:
        bam = barplot_prefix
        bam = os.path.basename(bam)

    if (df_start.shape[0]>0) & (df_end.shape[0]>0):
        g = sns.FacetGrid(df_start, col='footprint_len', height=3, sharey=False, ylim=ylim, xlim=(-21,21)) 
        g.map_dataframe(sns.histplot, x='adj_start', hue='frame', discrete=True)
        g.set(xlabel=None)
        g.set_titles(col_template='{col_name}-nt')
        g.fig.suptitle(str(downsampling*100) + '% of footprints sampled from ' + bam, y=1.05)
        g.fig.supxlabel('Footprint positions from start codons',fontsize='medium')
        plt.savefig(barplot_prefix + '.start_codon.pdf', bbox_inches='tight')
        
        g = sns.FacetGrid(df_end, col='footprint_len', height=3, sharey=False, ylim=ylim, xlim=(-21,21)) 
        g.map_dataframe(sns.histplot, x='adj_end', hue='frame', discrete=True)
        g.set(xlabel=None)
        g.set_titles(col_template='{col_name}-nt')
        g.fig.suptitle(str(downsampling*100) + '% footprints sampled from ' + bam, y=1.05)
        g.fig.supxlabel('Footprint positions from stop codons',fontsize='medium')
        plt.savefig(barplot_prefix + '.stop_codon.pdf', bbox_inches='tight')
        
        logging.info('saved metagene plots as ' + barplot_prefix + '.start_codon.pdf and ' + barplot_prefix + '.stop_codon.pdf')
        
    else:
        logging.warning(bam + ' have insufficient mapped reads for plotting! This is likely due to longer than expected read length.')
        
    


def analyse_footprints(offset_method, adj, bam, downsampling, cds_range, quality, outdir=None, min_size=25, max_size=35, ylim=None):
    """
    A wrapper for the above functions.

    Input:
        * offset_method: 5p or 3p (required)        
        * bam: BAM file path for ribosome footprint BAM (required)
        * downsampling: fraction used in bam_sampling (required)
        * cds_range: dataframe or tsv file from orf_finder or operon_finder (orfs.py) (required)
        * quality: quality of selected footprints. Option: best, good, or fair (required)
        * outdir: (default=None)
        * min_size: minimum footprint size for analysis (default=25)
        * max_size: maximum footprint size for analysis (default=35)
        * ylim: option to use the same ylim for all barplots (default=None)
        
    Output:
        * stat: RIBOSS statistics
        * heatmaps and metagene plots as PDFs
    """

    fname = filename(bam, None, outdir)

    sample = bam_sampling(bam, downsampling)

    if quality=='best':
        try:
            stat,fp,df = footprint_summary(cds_range, sample, 'best', offset_method, adj=adj, min_size=min_size, max_size=max_size, offset_prefix=fname)
            heatmap_periodicity(fp, fname)
            footprint_periodicity(df, downsampling, fname, ylim=ylim) # metagene plots
        except Exception as e:
            logging.exception('Reads are longer than expected! Are you sure that this is a ribosome profiling library?')
            pass
    elif quality=='good':
        try:
            stat,fp,df = footprint_summary(cds_range, sample, 'good', offset_method, adj=adj, min_size=min_size, max_size=max_size, offset_prefix=fname)
            heatmap_periodicity(fp, fname)
            footprint_periodicity(df, downsampling, fname, ylim=ylim) # metagene plots
        except Exception as e:
            logging.exception('Reads are longer than expected! Are you sure that this is a ribosome profiling library?')
            pass
        

    
    return stat