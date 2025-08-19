#!/usr/bin/env python
# coding: utf-8

"""
@author      CS Lim
@create date 2020-10-10 16:49:00
@modify date 2025-08-05 11:01:02
@desc        RIBOSS module for binary wrappers
"""



import subprocess, os, sys, logging.config
from io import StringIO
import pandas as pd
import os
from Bio import SeqIO


# DEFAULT_LOGGING = {
#     'version': 1,
#     'disable_existing_loggers': False,
#     'formatters': {
#         'detailed': {
#             'format': '%(asctime)s - %(levelname)s - %(name)s - %(message)s',
#             'datefmt': '%Y-%m-%d %H:%M:%S'
#         },
#     },
#     'handlers': {
#         'console': {
#             'class': 'logging.StreamHandler',
#             'formatter': 'detailed'
#         },
#     },
#     'loggers': {
#         '': {
#             'level': 'INFO',
#             'handlers': ['console'],
#         },
#         'another.module': {
#             'level': 'ERROR',
#             'handlers': ['console'],
#         },
#     }
# }

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




def filename(infname, outfname=None, outdir=None, pathbase=False):
    if outdir:
        os.makedirs(outdir, exist_ok=True)
        if pathbase==False:
            fname = os.path.basename(infname)
            fname = outdir + '/' + os.path.splitext(fname)[0]
            if outfname:
                fname = fname.replace('//','/') + '.' + outfname
            else:
                fname = fname.replace('//','/')
        else:
            path, ext = os.path.split(infname)
            fname = path + '/' + ext.split(os.extsep)[0]
            if outfname:
                fname = fname.replace('//','/') + outfname
            else:
                fname = fname.replace('//','/')
    else:
        if outfname:
            fname = './' + outfname
        else:
            fname = os.path.splitext(infname)[0]

    return fname



def transcriptome_assembly(superkingdom, genome, long_reads, short_reads=None, strandness=None, annotation=None, num_threads=4, trim=True,
                           min_length=100, coverage=1, single_exon_coverage=1.5, fraction=0.1, outdir=None):
    """
    Assembly transcriptome using long reads or a mix of short and long reads.
    
    Input:
        * superkingdom: Archaea, Bacteria or Eukaryota (required)
        * genome: genomic fasta file (required)
        * long_reads: BAM file for PacBio or Oxford Nanopore long reads (required)
        * short_reads: BAM file for Illumina long reads (default: None)
        * annotation: reference annotation in GTF (default: None)
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
    else:
        cmd = cmd + ['-L', long_reads]
    
    if annotation!=None:
        basename, ext = os.path.splitext(annotation)
        if (ext=='.gz') & (os.path.isfile(annotation)):
            subprocess.run(['gunzip', annotation], check=True)
            annotation = basename
        elif os.path.isfile(basename):
            annotation = basename
        else:
            pass
        cmd = cmd + ['-G', annotation]
        
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
    
    if annotation!=None:
        subprocess.run(['gzip', annotation], check=True)
    
    return fname + '.transcripts.fa', fname + '.gtf'



def fasta_to_dataframe(seq):
    fasta_df = pd.read_csv(seq, sep='>', lineterminator='>', header=None)
    df = fasta_df[0].str.split('\n', n=1, expand=True)
    df[1] = df[1].replace('\n','', regex=True)
    df = df[df[1] != '']
    df = df.dropna()
    df.columns = ['tid','seq']
    df['tid'] = df['tid'].str.split('|').str[0]
    return df

    

def fasta_record_generator(fasta_path):
    """
    Reads a FASTA file and yields each record as a dictionary.
    """
    handle = partial(gzip.open, mode='rt') if fasta_path.endswith('.gz') else open
    with handle(fasta_path, 'r') as f:
        for record in SeqIO.parse(f, 'fasta'):
            yield {'Name': record.id, 'Sequence': str(record.seq)}
            
            
def build_star_index(
        fasta_path,
        index,
        star='STAR',
        mode='genomeGenerate',
        num_threads=4,
        nbases=7,
        nbits=11,
        delim=None,
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

    fname = filename(fasta_path)
    fasta = fasta_to_dataframe(fasta_path)
    fasta['tid'] = '>' + fasta.tid.str.split().str[0].str.split(delim).str[0]
    
    if '.gz' in fasta_path:
        fasta.to_csv(fname, sep='\n', header=None, index=None)
        logging.info('cleaned up fasta headers and saved as ' + fname)
        fastafile = fname
    else:
        os.rename(fasta_path, fname + '.original.fasta')
        fasta.to_csv(fasta_path, sep='\n', header=None, index=None) 
        logging.info('renamed ' + fasta_path + ' as ' + fname + '.original.fasta')
        logging.info('cleaned up fasta headers and saved as ' + fasta_path)
        fastafile = fasta_path
        
    program = [star]
    options = [
        '--runMode', mode,
        '--runThreadN', str(num_threads),
        '--genomeDir', index,
        '--genomeFastaFiles', fastafile,
        '--genomeSAindexNbases', str(nbases), 
        '--genomeChrBinNbits', str(nbits)
        ]
    
    cmd = program + options
    subprocess.run(cmd, check=True)
      
    if os.path.isdir(index):
        logging.info('saved index to ' + index)
    else:
        logging.error('No index generated!')


def get_magic_number(filepath):
    """Safely reads the first two bytes of a file."""
    if not os.path.isfile(filepath) or os.path.getsize(filepath) < 2:
        raise TypeError("Invalid Read file path!")
    with open(filepath, 'rb') as f:
        return f.read(2)
    

def align_short_reads(
        reads,
        prefix,
        index,
        star='STAR',
        mode='alignReads',
        num_threads=4,
        read_files_command='zcat',
        seed_search_lmax=10,
        filter_multimap_score_range=0,
        filter_multimap_nmax=255,
        filter_mismatch_nmax=1,
        filter_intron_motifs='RemoveNoncanonical',
        sam_type='BAM SortedByCoordinate',
        # sam_mode='NoQS',
        # sam_attributes='NH NM',
        clip_3p_adapter_seq=None     
        ):
    """
    Wrapper to run STAR aligner.

    Input:
          * reads: str if single-end, a list of paired-end (required)
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
          # * sam_mode: (default: NoQS)
          # * sam_attributes: (default: NH NM)           

    Output:
          * output: BAM filename with prefix + Aligned.out.bam
    """
    
    program = [star]

    if type(reads)==str:
        options = [
            '--runMode', mode,
            '--readFilesIn', reads,# read2,
            '--outFileNamePrefix', prefix,
            '--genomeDir', index,
            '--runThreadN', str(num_threads), 
            '--seedSearchLmax', str(seed_search_lmax),
            '--outFilterMultimapScoreRange', str(filter_multimap_score_range),
            '--outFilterMultimapNmax', str(filter_multimap_nmax),
            '--outFilterMismatchNmax', str(filter_mismatch_nmax),
            '--outFilterIntronMotifs', filter_intron_motifs,
            '--outSAMtype', sam_type.split()[0], sam_type.split()[1],
            # '--outSAMattributes', sam_attributes.split()[0], sam_attributes.split()[1]
            ]
        
        magic_number = get_magic_number(reads)
        if magic_number == b'\x1f\x8b':
            print("Read file in gzip")
            options = options + ['--readFilesCommand', read_files_command]
        else:
            logging.warning(f'{reads} has not been compressed as gzip. Processing with mapping')
            
        if clip_3p_adapter_seq!=None:
            cmd = program + options + ['--clip3pAdapterSeq', str(clip_3p_adapter_seq)]
        else:
            cmd = program + options
            logging.warning('Adapter sequence is not provided! Please use the clip_3p_adapter_seq argument if it has not been trimmed.')
        
    elif type(reads)==list:
        options = [
            '--runMode', mode,
            '--readFilesIn', reads[0], reads[1],
            '--outFileNamePrefix', prefix,
            '--genomeDir', index,
            '--runThreadN', str(num_threads),
            '--seedSearchLmax', str(seed_search_lmax),
            '--outFilterMultimapScoreRange', str(filter_multimap_score_range),
            '--outFilterMultimapNmax', str(filter_multimap_nmax),
            '--outFilterMismatchNmax', str(filter_mismatch_nmax),
            '--outFilterIntronMotifs', filter_intron_motifs,
            '--outSAMtype', sam_type.split()[0], sam_type.split()[1],
            # '--outSAMattributes', sam_attributes.split()[0], sam_attributes.split()[1]
            ]
        
        magic_numbers = [get_magic_number(filepath) for filepath in reads]
        is_all_gzip = all(mn == b'\x1f\x8b' for mn in magic_numbers)

        if is_all_gzip:
            options = options + ['--readFilesCommand', read_files_command]
        else:
            logging.warning(f'{reads} have not been compressed as gzip. Processing with mapping')
            
        if clip_3p_adapter_seq!=None:
            cmd = program + options + ['--clip3pAdapterSeq', str(clip_3p_adapter_seq), str(clip_3p_adapter_seq), 
                                       '--clip3pAdapterMMp', str(0.1), str(0.1)]
        else:
            cmd = program + options
            logging.warning('Adapter sequence is not provided! Please use the clip_3p_adapter_seq argument if it has not been trimmed.')
    else:
        logging.error('Read input must be a string or a list!')

        
    subprocess.run(cmd, check=True)

    if sam_type.split()[1]=='SortedByCoordinate':
        bam = prefix + 'Aligned.sortedByCoord.out.bam'    
    elif sam_type.split()[1]=='Unsorted':
        bam = prefix + 'Aligned.out.bam'
    
    if os.path.exists(bam):
        logging.info('saved alignment to ' + bam)
    else:
        logging.error('No alignment file generated!')
        
    return bam




def quantify_transcripts(reads, fasta_path, adapter=None, index_prefix=None, outdir=None, overwrite=False):
    """
    Build pufferfish index and count reads by mRNA isoforms using Salmon.

    Input:
        * reads: str if single-end, a list if paired-end
        * fasta_path: path to FASTA file
        * index_prefix: prefix for index path
        * outdir: output directory
        * overwrite: whether or not to overwrite index (Default: False)

    Output:
        * [prefix]_puff and [prefix]_salmont_quant
    """

    if not index_prefix:
        logging.error('No index path provided!')
        sys.exit(1)

    index = index_prefix + '_puff'
    build_index = ['salmon', 'index', '-t', fasta_path, '-i', index]

    index_exists = os.path.exists(index)
    index_is_dir = os.path.isdir(index)

    if index_exists:
        if index_is_dir:
            logging.info(f'Index directory {index} exists.')
        else:  # index is a file
            logging.warning(f'{index} file exists. Consider renaming or moving it.')
            # sys.exit(1)  # Consider if you want to exit here
    elif overwrite:
        try:
            subprocess.run(build_index, check=True)
            logging.info(f'Saved index to {index}')
        except OSError as e:
            logging.error(f'Error creating index: {e}. Check permissions and Salmon installation.')
            sys.exit(1)
    else:  # index doesn't exist and not overwriting
        try:
            subprocess.run(build_index, check=True)
            logging.info(f'Saved index to {index}')
        except OSError as e:
            logging.error(f'Error creating index: {e}. Check permissions and Salmon installation.')
            sys.exit(1)


    if isinstance(reads, str):  # Single-end reads
        quant_dir = filename(reads, outfname='_salmon_quant/', outdir=outdir, pathbase=True)
        quant_args = ['salmon', 'quant', '-l', 'A', '-i', index, '-o', quant_dir, '-q']

        if adapter:
            trim_fastq = filename(reads, outfname='_trimmed.fastq.gz', outdir=outdir, pathbase=True)
            trim = ['fastp', '-i', reads, '-o', trim_fastq, '-q', '10', '-w', '8', '-a', adapter]
            subprocess.run(trim, check=True)
            quant_args.extend(['-r', trim_fastq])
        else:
            quant_args.extend(['-r', reads])
            logging.warning('No adapter trimming! Consider trimming if low number of transcripts quantified.\n')

        try:
            subprocess.run(quant_args, check=True)
        except subprocess.CalledProcessError as e:
            logging.error(f"Salmon quantification failed: {e}")
            sys.exit(1)


    elif isinstance(reads, list):  # Paired-end reads
        if len(reads) != 2:
            logging.error("For paired-end reads, provide a list of two file paths.")
            sys.exit(1)

        read1, read2 = reads  # Unpack read files
        quant_dir = filename(read1, outfname='_salmon_quant/', outdir=outdir, pathbase=True)
        quant_args = ['salmon', 'quant', '-l', 'A', '-i', index, '-o', quant_dir, '-q']

        if adapter:
            trim_fastqs = [filename(i, outdir=outdir, pathbase=True) for i in reads]
            trim = ['fastp', '-i', reads[0], '-I', reads[1],
                    '-o', trim_fastqs[0] + '_trimmed_1.fastq.gz',
                    '-O', trim_fastqs[1] + '_trimmed_2.fastq.gz',
                    '-q', '10', '-w', '8', '-a', adapter]
            subprocess.run(trim, check=True)
            quant_args.extend(['-1', trim_fastqs[0] + '_trimmed_1.fastq.gz',
                               '-2', trim_fastqs[1] + '_trimmed_2.fastq.gz'])
        else:
            quant_args.extend(['-1', reads[0], '-2', reads[1]])
            logging.warning('No adapter trimming! Consider trimming if low number of transcripts quantified.\n')

        try:
            subprocess.run(quant_args, check=True)
        except subprocess.CalledProcessError as e:
            logging.error(f"Salmon quantification failed: {e}")
            sys.exit(1)

    else:
        logging.error('Read input must be a string or a list!')
        sys.exit(1)

    if os.path.isdir(quant_dir):
        logging.info(f'Saved read counts to {quant_dir}')
        sf = pd.read_csv(os.path.join(quant_dir, 'quant.sf'), sep='\t') # Use os.path.join
        total_tx = sf.shape[0]
        quant_tx = sf[sf['NumReads'] != 0].shape[0] # Use string indexing
        logging.info(f'Quantified {quant_tx} of {total_tx} transcripts')
        return quant_dir # Return the quantification directory

    else:
        logging.error('No quant.sf generated!')
        return None  # Return None if quant.sf is not generated



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
        * min_fplen: minimum footprint size for analysis (default=25)
        * max_fplen: minimum footprint size for analysis (default=35)
    
    Output:
        * file path for .base. Output include .base., .codon, _abundant.list, and _scarce.list
    """

    cmd = [
        'riboprof',
        '--ribobam', '"{}"'.format(ribobam),
        '--mrnabam', '"{}"'.format(mrnabam),
        '--fasta', '"{}"'.format(fasta),
        '--cds_range', '"{}"'.format(cds_range_file),
        '--sf', '"{}"'.format(sf),
        '--tabd_cutoff', '"{}"'.format(str(tabd_cutoff)),
        '--offset', '"{}"'.format(offset),
        '--out', '"{}"'.format(out),
        '--useSecondary',
        '--min_fplen', '"{}"'.format(str(min_fplen)),
        '--max_fplen', '"{}"'.format(str(max_fplen))
        ]
    
    if os.path.isdir(out)==True:
        out = out + 'riboprof'
    elif (os.path.isdir(out)!=True) & (os.path.isdir(os.path.split(out)[0])!=True):
        os.makedirs(os.path.split(out)[0])
    
    if os.path.isfile(ribobam) & os.path.isfile(mrnabam) & os.path.isfile(fasta) & os.path.isfile(cds_range_file) & os.path.isfile(sf):
        subprocess.run(' '.join(cmd), check=True, shell=True)
    else:
        logging.error('Missing input file(s)!')
        
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