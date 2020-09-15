"""
@author      CS Lim
@create date 2020-09-15 17:40:16
@modify date 2020-09-16 06:34:51
@desc        RIBOSS pipline
"""

import os
work_dir = os.getcwd()

# IO
data_dir = work_dir + '/data/'
adapter = 'TCGTATGCCGTCTTCTGCTTG' # change me 'TGGAATTCTCGGGTGCCAAGG'
ribo_id = 'GSM546920' # change me
rna_id = 'GSM546921' # change me
ext = '_filtered_sequence.txt.gz' # change me
ribo_raw = data_dir + ribo_id + ext
rna_raw = data_dir + rna_id + ext
ribo_name = data_dir + ribo_id
rna_name = data_dir + rna_id
ribo_filt = ribo_name + '_Unmapped.out.mate1.gz'
rna_filt = rna_name + '_Unmapped.out.mate1.gz'
ribo_bam = ribo_name + '_Aligned.out.bam'
rna_bam = rna_name + '_Aligned.out.bam'
sm_dir = data_dir + 'sm_quant/'
sm_quant = sm_dir + 'quant.sf'
offset = ribo_name + '_Aligned.out_offset.txt'
riboprof_base = data_dir + ribo_id + '.base'
phase = data_dir + ribo_id + '_periodicity.txt'
stat = data_dir + ribo_id + '_fisher_test.txt'

# Reference sequences
# change me
species = 'Homo_sapiens'
database = 'gencode'
ref_dir = work_dir + '/ref/' + species + '/'
transcripts = ref_dir + 'gencode.v35.pc_transcripts.fa.gz',
pc_transcripts = ref_dir + 'gencode.v35.pc_transcripts_filtered.fasta'
cds_range = ref_dir + 'gencode.v35.pc_transcripts_filtered.txt'
pickle = ref_dir + 'gencode.v35.pc_transcripts_filtered.pkl.gz'
ncrna = ref_dir + 'Homo_sapiens.GRCh38.ncrna.fa.gz',
trna = ref_dir + 'eukaryotic-tRNAs.fa.gz'
contam = ref_dir + 'contaminant.fa'
# Index
star_contam_index = ref_dir + 'star_contam_index'
star_riboseq_index = ref_dir + 'star_riboseq_index'

# Tools
scripts_dir = work_dir + '/scripts/'
prep_ref = scripts_dir + 'prep_reference_files.py'
sel_fp = scripts_dir + 'select_footprint_size.py'
riboss = scripts_dir + 'riboss.py'
bin_dir = work_dir + '/bin/'
riboprof = bin_dir + 'riboprof'

# Parameters
thread = 20
contam_params = '--readFilesCommand zcat \
                --seedSearchLmax 10 \
                --outStd SAM \
                --outReadsUnmapped Fastx \
                --outFilterMultimapScoreRange 0 \
                --outFilterMultimapNmax 255 \
                --outFilterMismatchNmax 1 \
                --outFilterIntronMotifs RemoveNoncanonical > /dev/null'
map_params =    '--readFilesCommand zcat \
                --seedSearchLmax 10 \
                --outFilterMultimapScoreRange 0 \
                --outFilterMultimapNmax 255 \
                --outFilterMismatchNmax 1 \
                --outFilterIntronMotifs RemoveNoncanonical \
                --outSAMtype BAM Unsorted --outSAMmode NoQS \
                --outSAMattributes NH NM'
salmon_params = '-l A --seqBias --posBias'
riboprof_params = '--min_fplen 25 --max_fplen 35 --tabd_cutoff 0'

                
# 1
rule all:
    input:
        transcripts, 
        ncrna,
        trna,
        contam, 
        pc_transcripts,
        cds_range,
        pickle,
        star_contam_index,
        star_riboseq_index,
        ribo_raw,
        rna_raw,
        ribo_filt, 
        rna_filt,
        ribo_bam,
        rna_bam,
        offset,
        sm_quant,
        riboprof_base,
        phase,
        stat

# 2
rule fastq:
    output:
        ribor=ribo_raw,
        rnar=rna_raw
    params:
        data=data_dir
    run:
        shell(r"""
        mkdir -p {params.data};
        wget -P {params.data} -N ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM546nnn/GSM546920/suppl/GSM546920_filtered_sequence.txt.gz;
        wget -P {params.data} -N ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM546nnn/GSM546921/suppl/GSM546921_filtered_sequence.txt.gz;
        """)

# 3
rule reference:
    output:
        tx=transcripts, 
        nc=ncrna, 
        t=trna
    params:
        ref=ref_dir
    run:
        shell(r"""
        mkdir -p {params.ref};
        wget -P {params.ref} -N ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_35/gencode.v35.pc_transcripts.fa.gz;
        wget -P {params.ref} -N ftp://ftp.ensembl.org/pub/release-101/fasta/homo_sapiens/ncrna/Homo_sapiens.GRCh38.ncrna.fa.gz;
        wget -P {params.ref} -N http://gtrnadb.ucsc.edu/download/tRNAs/eukaryotic-tRNAs.fa.gz;
        """)

# 4
rule contamination:
    input:
        nc=ncrna, 
        t=trna
    output:
        c=contam
    params:
        sp=species
    run:
        shell(r"""
        zcat {input.nc} |
        awk 'BEGIN{{RS=">"}} NR>1 {{sub("\n","\t"); gsub("\n",""); print RS$0}}' |
        awk '/:Mt_rRNA|:Mt_tRNA|:rRNA|:scRNA|:snlRNA|:snoRNA|:snRNA|:tRNA|:tRNA_pseudogene/ && !seen[$NF]++ {{print $1 "\n" $NF}}' > {output.c};
        zcat {input.t} |
        awk 'BEGIN{{RS=">"}} NR>1 {{sub("\n","\t"); gsub("\n",""); print RS$0}}' |
        grep {params.sp} |
        awk '!seen[$NF]++ {{print $1 "\n" $NF}}' >> {output.c};
        """)

# 5     
rule prep_reference:
    input:
        tx=transcripts
    output:
        pc=pc_transcripts,
        cdsr=cds_range,
        pkl=pickle
    params:
        prep=prep_ref, 
        db=database, 
        ref=ref_dir
    run:
        shell(r"""
        python {params.prep} -d {params.db} -s {input.tx} -o {params.ref}
        """)

# 6
rule index:
    input:
        c=contam, 
        pc=pc_transcripts
    output:
        ci=directory(star_contam_index),
        pci=directory(star_riboseq_index)
    params:
        t=thread,
    run:
        shell(r"""
        STAR \
        --runThreadN {params.t} \
        --runMode genomeGenerate \
        --genomeDir {output.ci} \
        --genomeFastaFiles {input.c} \
        --genomeSAindexNbases 7 \
        --genomeChrBinNbits 11;

        STAR \
        --runThreadN {params.t} \
        --runMode genomeGenerate \
        --genomeDir {output.pci} \
        --genomeFastaFiles {input.pc} \
        --genomeSAindexNbases 11 \
        --genomeChrBinNbits 12;
        """)

# 7
rule filter:
    input:
        ribor=ribo_raw, 
        rnar=rna_raw,
        ci=star_contam_index
    output:
        ribof=ribo_filt, 
        rnaf=rna_filt,
    params:
        t=thread, 
        a=adapter,
        ribon=ribo_name + '_', 
        rnan=rna_name + '_',
        cp=contam_params,
        data=data_dir
    run:
        shell(r"""
        STAR \
        --runThreadN {params.t} \
        --genomeDir {input.ci} \
        --readFilesIn {input.ribor} \
        --clip3pAdapterSeq {params.a} \
        --outFileNamePrefix {params.ribon} \
        {params.cp};

        STAR \
        --runThreadN {params.t} \
        --genomeDir {input.ci} \
        --readFilesIn {input.rnar} \
        --clip3pAdapterSeq {params.a} \
        --outFileNamePrefix {params.rnan} \
        {params.cp};

        gzip {params.data}*.mate1;
        """)

# 8
rule alignment:
    input:
        ribof=ribo_filt, 
        rnaf=rna_filt,
        pci=star_riboseq_index
    output:
        ribob=ribo_bam,
        rnab=rna_bam,
    params:
        t=thread, 
        a=adapter, 
        ribon=ribo_name + '_', 
        rnan=rna_name + '_',
        mp=map_params
    run:
        shell(r"""
        STAR \
        --runThreadN {params.t} \
        --genomeDir {input.pci} \
        --readFilesIn {input.ribof} \
        --outFileNamePrefix {params.ribon} \
        --clip3pAdapterSeq {params.a} \
        {params.mp};

        STAR \
        --runThreadN {params.t} \
        --genomeDir {input.pci} \
        --readFilesIn {input.rnaf} \
        --outFileNamePrefix {params.rnan} \
        --clip3pAdapterSeq {params.a} \
        {params.mp};
        """)

# 9     
rule select_footprint:
    input:
        ribob=ribo_bam,
        cdsr=cds_range
    output:
        ost=offset
    params:
        sel=sel_fp, 
        data=data_dir
    run:
        shell(r"""
        python {params.sel} -b {input.ribob} -p {input.cdsr} -o {params.data}
        """)

# 10
rule quant:
    input:
        rnab=rna_bam,
        pc=pc_transcripts
    output:
        smq=sm_quant
    params:
        sp=salmon_params,
        t=thread,
        smd=sm_dir
    run:
        shell(r"""
        salmon quant -p {params.t} -t {input.pc} -a {input.rnab} -o {params.smd} {params.sp};
        """)

# 11
rule riboprof:
    input:
        ribob=ribo_bam,
        rnab=rna_bam,
        pc=pc_transcripts,
        cdsr=cds_range,
        ost=offset,
        smq=sm_quant,
    output:
        base=riboprof_base
    params:
        rbp=riboprof_params,
        riboprof=riboprof,
        ribon=ribo_name
    run:
        shell(r"""     
        {params.riboprof} \
        --ribobam {input.ribob} \
        --mrnabam {input.rnab} \
        --fasta {input.pc} \
        --cds_range {input.cdsr} \
        --offset {input.ost} \
        --sf {input.smq} \
        --out {params.ribon} \
        {params.rbp};
        """)

# 12          
rule riboss:
    input:
        base=riboprof_base,
        pkl=pickle
    output:
        phase=phase,
        stat=stat
    params:
        riboss=riboss,
        data=data_dir
    run:
        shell(r"""
        python {params.riboss} -i {input.base} -r {input.pkl} -o {params.data}
        """)
