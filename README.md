![logo](./riboss_logo.svg)

## Comparing the translatability of open reading frames within individual transcripts
See styphimurium.ipynb for analysis of _Salmonella enterica_ serovar Typhimurium.

### Install dependencies

#### install Anaconda

```
wget https://repo.anaconda.com/archive/Anaconda3-2024.06-1-Linux-x86_64.sh -O ~/Anaconda3.sh
bash ~/Anaconda3.sh -b -p $HOME/Anaconda3
```

#### Create a conda environment and install packages

```
conda create -n riboss -y
conda activate riboss
pip install tables cython pysam quicksect cgatcore pandarallel rseqc
conda install -c bioconda -c conda-forge biopython pysam htslib bedtools minimap2 star fastp tqdm
git clone https://github.com/cgat-developers/cgat-apps.git
cd cgat-apps
python setup.py develop
```

#### Install samtools

```
wget https://github.com/samtools/samtools/releases/download/1.21/samtools-1.21.tar.bz2
tar jxvf samtools-1.21.tar.bz2
cd samtools-1.21
autoheader
autoconf -Wno-syntax
make
pwd | awk '{print "export PATH=\"" $1 ":$PATH\""}' >> ~/.bashrc
```

#### Download precompiled tools

```
git clone https://github.com/lcscs12345/riboss.git
DIRNAME=`which python | xargs dirname`
cp riboss/bin/riboprof $DIRNAME
chmod +x $DIRNAME/riboprof

wget https://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/gtfToGenePred
wget https://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/genePredToBed
mv gtfToGenePred genePredToBed $DIRNAME
# chmod +x $DIRNAME/gtfToGenePred
# chmod +x $DIRNAME/genePredToBed

wget https://github.com/gpertea/stringtie/releases/download/v2.2.3/stringtie-2.2.3.Linux_x86_64.tar.gz
tar zxvf stringtie-2.2.3.Linux_x86_64.tar.gz
mv stringtie-2.2.3.Linux_x86_64/stringtie $DIRNAME
# chmod +x $DIRNAME/stringtie

cd DIRNAME
wget https://github.com/COMBINE-lab/salmon/releases/download/v1.10.0/salmon-1.10.0_linux_x86_64.tar.gz
tar zxvf salmon-1.10.0_linux_x86_64.tar.gz
pwd | awk '{print "export PATH=\"" $1 "/salmon-latest_linux_x86_64/bin:$PATH\""}' >> ~/.bashrc
```

#### Related article:

- Lim, C.S., Wardell, S.J.T., Kleffmann, T. & Brown, C.M. (2018) The exon-intron gene structure upstream of the initiation codon predicts translation efficiency. Nucleic Acids Res. 46: 4575-4591.
