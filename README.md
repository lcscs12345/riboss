![logo](doc/riboss_logo.svg)

## Comparing translation of open reading frames within individual transcripts

RIBOSS consists of Python modules for analysis of ribosome profiling data for prokaryotes and eukaryotes. See [`styphimurium.ipynb`](https://github.com/lcscs12345/riboss/blob/master/styphimurium.ipynb) where RIBOSS detects new ORFs in _S_. Typhimurium operons. This example starts from transcriptome assembly using long- and short-read RNA-seq data. It leverages the newly assembled transcriptome and highly phased ribosome profiling data to discover novel translation events.

![Flow Chart](doc/flow_chart.svg)

### User guide

#### Install Miniforge3 and create a conda environment

```
wget https://github.com/conda-forge/miniforge/releases/download/24.7.1-2/Miniforge3-24.7.1-2-Linux-x86_64.sh
bash Miniforge3-24.7.1-2-Linux-x86_64.sh -b -p $HOME/miniforge3
eval "$(/$HOME/miniforge3/bin/conda shell.bash hook)" # your terminal prompt will show (base) bash-5.1$
```

<!-- conda create -n riboss -y
conda activate riboss
conda install -y \
    -c conda-forge -c bioconda \
    boost-cpp seqan-library=1.4.2 \
    jupyter pandas \
    pysam seaborn matplotlib \
    stringtie=2.2.3 salmon \
    biopython htslib samtools bedtools pyranges minimap2 star tqdm jupyter \
    ucsc-gtftogenepred ucsc-bedtogenepred ucsc-genepredtobed ucsc-bedsort ucsc-bedtobigbed \
    pyfaidx rseqc
conda activate riboss
conda install bioconda::bowtie2 -y
conda env export > environment.yml -->

#### Install RIBOSS

```
git clone https://github.com/lcscs12345/riboss.git
cd riboss
conda env create -f environment.yml
conda activate riboss # your terminal prompt will show (riboss) bash-5.1$
DIRNAME=`which python | xargs dirname`
cp bin/riboprof $DIRNAME
chmod +x $DIRNAME/riboprof
pip install -e . # editable mode
```

<!-- pip install git+git://github.com/lcscs12345/riboss.git#egg=riboss -->

#### Activate the conda environment for next time

```
eval "$(/$HOME/miniforge3/bin/conda shell.bash hook)"
conda activate riboss
```

#### Basic usage

```
from riboss.orfs import translate
na='ATGGTCTGA'
translate(na)
```
You should see `'MV'`. For detail usage, see [`styphimurium.ipynb`](https://github.com/lcscs12345/riboss/blob/master/styphimurium.ipynb)

#### To reproduce plots in the manuscript

Create new directories `mkdir -p doc/ doc/metatranscriptome doc/styphimurium/ doc/styphimurium/rnaseq doc/styphimurium/riboseq`.

The alignment files for transcriptome assembly and ribosome profiling are available at [Zenodo](https://doi.org/10.5281/zenodo.13997374).

- SRR11215003.bam and SRR11215004.bam: Nanopore long-read direct RNA-seq. Download and `mv SRR24781620.bam doc/metatranscriptome`.
- SRR11215663.bam and SRR11215664.bam: Illumina short-read RNA-seq. Download and `mv SRR24781620.bam doc/metatranscriptome`.
- SRR24781620.bam: Nanopore long-read cDNA sequencing. Download and `mv SRR24781620.bam doc/styphimurium/rnaseq`.

Download the Ribosome profiling alignment files and `mv ERR913094*.out.bam doc/styphimurium/riboseq`.

- ERR9130942Aligned.out.bam: RNase I, 1000 U.
- ERR9130943Aligned.out.bam: RNase I, 500 U.
- ERR9130946Aligned.out.bam: matched RNA-seq.

### References:

- Lim, C. S., & Brown, C. M. (2024). RIBOSS detects novel translational events by combining long- and short-read transcriptome and translatome profiling. _BioRxiv_, DOI: [10.1101/2024.11.07.622529](https://doi.org/10.1101/2024.11.07.622529)
- Lim, C.S., Wardell, S.J.T., Kleffmann, T. & Brown, C.M. (2018) The exon-intron gene structure upstream of the initiation codon predicts translation efficiency. _Nucleic Acids Res_, 46:4575-4591. DOI: [10.1093/nar/gky282](https://doi.org/10.1093/nar/gky282)
- Bryant, O.J., Lastovka, F., Powell, J. et al. (2023) The distinct translational landscapes of gram-negative _Salmonella_ and gram-positive _Listeria_. _Nat Commun_, 14:8167. DOI: [10.1038/s41467-023-43759-1](https://doi.org/10.1038/s41467-023-43759-1)
- Yang, M., Cousineau, A., Liu, X., Luo, Y., Sun, D., Li, S., Gu, T., Sun, L., Dillow, H., Lepine, J., Xu, M., Zhang, B. (2020) Direct Metatranscriptome RNA-seq and Multiplex RT-PCR Amplicon Sequencing on Nanopore MinION - Promising Strategies for Multiplex Identification of Viable Pathogens in Food. _Front Microbiol_, 11:514. DOI: [10.3389/fmicb.2020.00514](https://doi.org/10.3389/fmicb.2020.00514)
