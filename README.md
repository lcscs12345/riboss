![logo](doc/riboss_logo.svg)

## Comparing the translatability of open reading frames within individual transcripts

RIBOSS consists of Python modules for analysis of ribosome profiling data for prokaryotes and eukaryotes. See [`styphimurium.ipynb`](https://github.com/lcscs12345/riboss/blob/master/styphimurium.ipynb) where RIBOSS detects new ORFs in _S_. Typhimurium operons. This example starts from transcriptome assembly using long- and short-read RNA-seq data. It leverages the newly assembled transcriptome and highly phased ribosome profiling data to discover novel translation events.

![Flow Chart](doc/flow_chart.svg)

### User guide

#### Install Miniforge3 and create a conda environment

```
wget https://github.com/conda-forge/miniforge/releases/download/24.7.1-2/Miniforge3-24.7.1-2-Linux-x86_64.sh
bash Miniforge3-24.7.1-2-Linux-x86_64.sh -b -p $HOME/miniforge3
eval "$(/$HOME/miniforge3/bin/conda shell.bash hook)"
conda env create -f environment.yml
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
conda env export > environment.yml -->

#### Install RIBOSS and dependencies

```
conda activate riboss
git clone https://github.com/lcscs12345/riboss.git

DIRNAME=`which python | xargs dirname`
cp riboss/bin/riboprof $DIRNAME
chmod +x $DIRNAME/riboprof
```

#### Activate the conda environment for next time

```
eval "$(/$HOME/miniforge3/bin/conda shell.bash hook)"
conda activate riboss
```

#### Test instructions

The alignment files for transcriptome assembly and ribosome profiling are available at [Zenodo](https://doi.org/10.5281/zenodo.13997374).
Download the RNA-seq alignment files and `mv D23005*.bam doc/styphimurium/rnaseq`.

- D23005.sorted.bam: Illumina short-read alignment file.
- D23005-sc-1962750.sorted.bam: PacBio long-read alignment file.

Download the Ribosome profiling alignment files and `mkdir doc/styphimurium/riboseq; mv ERR913094*.out.bam doc/styphimurium/riboseq`.

- ERR9130942Aligned.out.bam: RNase I, 1000 U.
- ERR9130943Aligned.out.bam: RNase I, 500 U.
- ERR9130946Aligned.out.bam: matched RNA-seq.

Follow the steps in `test.ipynb`.

### References:

- Lim, C.S., Wardell, S.J.T., Kleffmann, T. & Brown, C.M. (2018) The exon-intron gene structure upstream of the initiation codon predicts translation efficiency. _Nucleic Acids Res_, 46:4575-4591. DOI: [10.1093/nar/gky282](https://doi.org/10.1093/nar/gky282)
- Bryant, O.J., Lastovka, F., Powell, J. et al. (2023) The distinct translational landscapes of gram-negative _Salmonella_ and gram-positive _Listeria_. _Nat Commun_, 14:8167. DOI: [10.1038/s41467-023-43759-1](https://doi.org/10.1038/s41467-023-43759-1)
- Ondari, E.M., Klemm, E.J., Msefula, C.L. et al. (2019) Rapid transcriptional responses to serum exposure are associated with sensitivity and resistance to antibody-mediated complement killing in invasive _Salmonella_ Typhimurium ST313 [version 1; peer review: 2 approved]. _Wellcome Open Res_, 4:74. DOI: [10.12688/wellcomeopenres.15059.1](https://doi.org/10.12688/wellcomeopenres.15059.1)
