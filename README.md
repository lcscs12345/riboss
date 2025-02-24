![logo](doc/riboss_logo.svg)

## Comparing translation of open reading frames within individual transcripts

Ribosome profiling is routinely used for discovering actively translated ORFs. Standard ribosome profiling involves RNase digestion of ribosome-protected mRNA fragments, followed by sucrose gradient fractionation, and RNA sequencing (RNA-seq).

As ribosomes progress along the mRNA codon-by-codon, they generate a characteristic triplet periodicity profile in the footprint data. Triplet periodicity can be used to determine the correct reading frame for the translated ORFs and distinguish true translation events from background noise.

RIBOSS consists of Python modules for analysis of ribosome profiling data for prokaryotes and eukaryotes. See the use cases and benchmarking results [here](https://github.com/lcscs12345/riboss_paper).

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
which python | awk 'sub(/python/,"pip3") {print $1, "install -e ."}' | sh # editable mode
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

You should see `'MV'`.

### References:
- Ingolia, N. T., Hussmann, J. A., & Weissman, J. S. (2019) Ribosome Profiling: Global Views of Translation. Cold Spring Harb. Perspect. Biol., 11. DOI: [10.1101/cshperspect.a032698](https://doi.org/10.1101/cshperspect.a032698)
- Ingolia, N. T., Ghaemmaghami, S., Newman, J. R. S., & Weissman, J. S. (2009) Genome-wide analysis in vivo of translation with nucleotide resolution using ribosome profiling. Science, 324: 218–223. DOI: [10.1126/science.1168978](https://doi.org/10.1126/science.1168978)
- Lim, C. S., & Brown, C. M. (2024). RIBOSS detects novel translational events by combining long- and short-read transcriptome and translatome profiling. _BioRxiv_. DOI: [10.1101/2024.11.07.622529](https://doi.org/10.1101/2024.11.07.622529)
- Lim, C.S., Wardell, S.J.T., Kleffmann, T. & Brown, C.M. (2018) The exon-intron gene structure upstream of the initiation codon predicts translation efficiency. _Nucleic Acids Res_, 46:4575-4591. DOI: [10.1093/nar/gky282](https://doi.org/10.1093/nar/gky282)
