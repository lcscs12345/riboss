#!/usr/bin/env python
# coding: utf-8

"""
@author      CS Lim
@create date 2020-10-10 16:49:00
@modify date 2024-10-17 17:53:55
@desc        RIBOSS module for binary wrappers
"""



import subprocess, os, sys, logging
from io import StringIO
import pandas as pd



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




def filename(infname, outfname=None, outdir=None):
    if outdir:
        os.makedirs(outdir, exist_ok=True)
        fname = os.path.basename(infname)
        fname = outdir + '/' + os.path.splitext(fname)[0]
        if outfname:
            fname = fname.replace('//','/') + '.' + outfname
        else:
            fname = fname.replace('//','/')
    else:
        if outfname:
            fname = './' + outfname
        else:
            fname = os.path.splitext(infname)[0]

    return fname


    
def transcriptome_assembly(superkingdom, genome, long_reads, short_reads=None, strandness=None, num_threads=4, trim=True,
                           min_length=100, coverage=5, single_exon_coverage=5, fraction=0.1, outdir=None):
    """
    Assembly transcriptome using long reads or a mix of short and long reads.
    
    Input:
        * superkingdom: Archaea, Bacteria or Eukaryota (required)
        * genome: genomic fasta file (required)
        * long_reads: BAM file for PacBio or Oxford Nanopore long reads (required)
        * short_reads: BAM file for Illumina long reads (default: None)
        * strandness: rf assumes a stranded library fr-firststrand, fr assumes a stranded library fr-secondstrand (default: None)
        * outdir: output directory (default: None)

    Output:
        * StringTie GTF: .gtf
        * transcript BED: .bed
        * transcript fasta: .transcripts.fa
    """
    
    if superkingdom in ['Archaea','Bacteria']:
        junction_coverage=1000
        
    elif superkingdom=='Eukaryota':
        junction_coverage=1
        
    else:
        sys.exit('Superkingdom is required! Choose either Archaea, Bacteria or Eukaryota.')
        
    fname = filename(long_reads, None, outdir)

    cmd = ['stringtie',
           '-p', str(num_threads),
           '-m', str(min_length), 
           '-c', str(coverage),
           '-s', str(single_exon_coverage),
           '-f', str(fraction), 
           '-j', str(junction_coverage), 
           '-o', fname + '.gtf']
    
    if short_reads!=None:
        cmd = cmd + ['--mix', short_reads,long_reads, '--' + strandness]
        # subprocess.run(, check=True)
    else:
        cmd = cmd + ['-L', long_reads]
        
    if trim==False:
        cmd = cmd + ['-t']
        
    subprocess.run(cmd, check=True)
    subprocess.run(['gtfToGenePred', fname + '.gtf', fname + '.gp'], check=True)
    subprocess.run(['genePredToBed', fname + '.gp', fname + '.bed'], check=True)

    track = 'track name="Transcriptome" description="Transcriptome assembly" visibility=2 colorByStrand="255,0,0 0,0,255"'
    subprocess.run(['sed','-i','1i '+ track, fname + '.bed'], check=True)
    
    subprocess.run(['bedtools', 'getfasta',
                    '-bed', fname + '.bed',
                    '-fi', genome,
                    '-s', '-name', 
                    '-fo', fname + '.transcripts.fa'])

    if os.path.exists(fname + '.gtf') & os.path.exists(fname + '.bed') & os.path.exists(fname + '.transcripts.fa'):
        logging.info('saved StringTie GTF as ' + fname + '.gtf')
        logging.info('converted transcriptome assembly to ' + fname + '.bed')
        logging.info('extracted sequences to ' + fname + '.transcripts.fa')
    else:
        logging.error('No results generated!')
        
    os.remove(fname + '.gp')

    return fname + '.transcripts.fa', fname + '.gtf'



def build_star_index(
        fasta_path,
        index,
        star='STAR',
        mode='genomeGenerate',
        num_threads=4,
        nbases=7,
        nbits=11
        ):
    """
    Wrapper to generate STAR index

    Input:
          * read1: full path to the read 1 FASTQ file (required)
          * prefix: the prefix to use for output files (required)
          * index: full path to the directory containing STAR genome index files (required)
          * star: full path to STAR (default: STAR)
          * num_threads: number of threads to use (default: 4)
          * mode: (default: genomeGenerate)
          * fasta_path: file path for fasta
          * nbases: (default: 7)
          * nbits: (default: 11)        

    Output:
          * cmd: command to execute for alignment
          * output: name of output file alignment is written to
    """
    
    program = [star]
    options = [
        '--runMode', mode,
        '--runThreadN', str(num_threads),
        '--genomeDir', index,
        '--genomeFastaFiles', fasta_path,
        '--genomeSAindexNbases', str(nbases), 
        '--genomeChrBinNbits', str(nbits)
        ]
    
    cmd = program + options
    subprocess.run(cmd, check=True)
      
    if os.path.isdir(index):
        logging.info('saved index to ' + index)
    else:
        logging.error('No index generated!')



def align_short_reads(
        read1,
        prefix,
        index,
        # clip_3p_adapter_seq=None,
        star='STAR',
        mode='alignReads',
        num_threads=4,
        read_files_command='zcat',
        seed_search_lmax=10,
        filter_multimap_score_range=0,
        filter_multimap_nmax=255,
        filter_mismatch_nmax=1,
        filter_intron_motifs='RemoveNoncanonical',
        sam_type='BAM Unsorted',
        sam_mode='NoQS',
        sam_attributes='NH NM',        
        ):
    """
    Wrapper to run STAR aligner.

    Input:
          * read1: full path to the read 1 FASTQ file (required)
          * prefix: the prefix to use for output files (required)
          * index: full path to the directory containing STAR genome index files (required)
          * star: full path to STAR (default: STAR)
          * mode: (default: alignReads)
          * num_threads: number of threads to use (default: 4)
          * read_files_command: uncompress program to run on FASTQ files (default: zcat)
          * seed_search_lmax:  (default: 10)
          * filter_multimap_score_range: (default: 0)
          * filter_multimap_nmax: (default: 255)
          * filter_mismatch_nmax: (default: 1)
          * filter_intron_motifs: filter alignments based on intron motifs (default: RemoveNoncanonical)
          * sam_type: output SAM file type (default: BAM SortedByCoordinate)
          * sam_mode: (default: NoQS)
          * sam_attributes: (default: NH NM)           

    Output:
          * output: BAM filename with prefix + Aligned.out.bam
    """
    
    program = [star]
    options = [
        '--runMode', mode,
        '--readFilesIn', read1,# read2,
        '--outFileNamePrefix', prefix,
        '--genomeDir', index,
        '--runThreadN', str(num_threads),
        '--readFilesCommand', read_files_command,
        # '--clip3pAdapterSeq', clip_3p_adapter_seq,
        '--seedSearchLmax', str(seed_search_lmax),
        '--outFilterMultimapScoreRange', str(filter_multimap_score_range),
        '--outFilterMultimapNmax', str(filter_multimap_nmax),
        '--outFilterMismatchNmax', str(filter_mismatch_nmax),
        '--outFilterIntronMotifs', filter_intron_motifs,
        '--outSAMtype', sam_type.split()[0], sam_type.split()[1],
        '--outSAMmode', sam_mode,
        '--outSAMattributes', sam_attributes.split()[0], sam_attributes.split()[1]
        ]
    
    cmd = program + options
    subprocess.run(cmd, check=True)

    bam = read1.split(os.extsep)[0] + 'Aligned.out.bam'
    
    if os.path.exists(bam):
        logging.info('saved alignment to' + bam)
    else:
        logging.error('No alignment file generated!')
        
    return bam



def quantify_transcripts(read1, fasta_path, index=None):
    """
    Build pufferfish index and count reads by mRNA isoforms using Salmon.

    Input:
        * read1: single-end read
        * fasta_path: path or fasta file
        * index: path for index

    Output:
        * [prefix]_puff and [prefix]_salmont_quant
    """
    
    build_index = ['salmon','index',
                   '-t', fasta_path, 
                   '-i', index + '_puff']

    quant_dir = read1.split(os.extsep)[0] + '_salmon_quant/'
    quant = ['salmon','quant',
             '-l', 'A', 
             '-i', index + '_puff', 
             '-r', read1,
             '-o', quant_dir, '-q']
    
    if index!=None:
        if os.path.isdir(index + '_puff/') is False:
            subprocess.run(build_index, check=True)
            subprocess.run(quant, check=True)
    else:
        subprocess.run(quant, check=True)

    if os.path.isdir(index + '_puff/'):
        logging.info('saved index to ' + index + '_puff/')
    else:
        logging.error('No index generated!')

    if os.path.isdir(quant_dir):
        logging.info('saved read counts to ' + quant_dir)
    else:
        logging.error('No quant.sf generated!')



def riboprofiler(
    offset,
    ribobam,
    mrnabam,
    fasta,
    cds_range_file,
    sf,
    out,
    tabd_cutoff=0,    
    min_fplen=25,
    max_fplen=35):
    
    """
    A wrapper for Riboprof.
    See https://github.com/Kingsford-Group/ribomap

    Input:
        * offset: text file of footprint size(s) with offset value.
        * ribobam: BAM file for ribosome profiling (required)
        * mrnabam: BAM file for paired RNA-seq (required)
        * fasta: transcript fasta file (required, e.g. tx_assembly extracted using bedtools getfasta)
        * cds_range_file: tsv file from orf_finder or operon_finder (orfs.py) (required)
        * sf: quan.sf from Salmon (required)
        * out: prefix for output files (required)
        * tabd_cutoff: cutoff for transcript abundance (default=0)
        * min_fplen: minimum footprint size for analysis (default=23)
        * max_fplen: minimum footprint size for analysis (default=35)
    
    Output:
        * file path for .base. Output include .base., .codon, _abundant.list, and _scarce.list
    """
        
    cmd = [
        'riboprof',
        '--ribobam', ribobam,
        '--mrnabam', mrnabam,
        '--fasta', fasta,
        '--cds_range', cds_range_file,
        '--sf', sf,
        '--tabd_cutoff', str(tabd_cutoff),
        '--offset', offset,
        '--out', out,
        '--useSecondary',
        '--min_fplen', str(min_fplen),
        '--max_fplen', str(min_fplen),
        ]
    
    if os.path.isdir(out)==True:
        out = out + 'riboprof'
    elif (os.path.isdir(out)!=True) & (os.path.isdir(os.path.split(out)[0])!=True):
        os.makedirs(os.path.split(out)[0])
    
    subprocess.run(cmd, check=True)

    if os.path.exists(out + '.base'):
        logging.info('saved main output as ' + out + '.base')
    else:
        logging.error('No output generated!')

    return out + '.base'



def merge_scores(bg):
    """
    Merge continuous BedGraph lines with the same scores
    https://www.biostars.org/p/267712/#267731

    Input:
        * bg: BedGraph in single-base resolution (required)

    Output:
        * dataframe for BedGraph
    """
    
    cmd = 'bedSort ' + bg  + ' stdout | groupBy -g 1,4 -c 2,3 -o min,max'
    ps = subprocess.Popen(cmd ,shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    output, error = ps.communicate()
    mbg = pd.read_csv(StringIO(output.decode('utf-8')), sep='\t', header=None)[[0,2,3,1]]
        
    return mbg