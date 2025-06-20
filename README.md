# RNA-seq workflow

Authors: Laura Barreales & Sara Lévano  
Date: 18/6/25

## What is this repository for?

This repository contains a series of directories to follow a RNA-seq analysis workflow. In each directory, there are two main folders: one for logs and other for results. Miniforge3 has to be installed in order to perform the RNAseq.

Due to the nature of the bioinformatic tools used, the stdout of some tools will be redirected to stderr.log by default.

This scripts were used in order to perform an analysis of 6 samples obtained from this article:

Panigrahi G, Candia J, Dorsey TH, et al. (2023). Diabetes-associated breast cancer is molecularly distinct and shows a DNA damage repair deficiency. JCI Insight, 8(23): e170105. Publicado el 8 de diciembre de 2023\. doi: 10.1172/jci.insight.170105 [https://pubmed.ncbi.nlm.nih.gov/37906280/](https://pubmed.ncbi.nlm.nih.gov/37906280/)

Then, this workflow is specifically created to obtain the biological processes that are differentially regulated comparing high glucose conditions and low glucose conditions in human breast cancer samples.

Regarding this, it is creted to perform a pair-end RNA-seq analysis with forward and reverse lectures and single-end analysis cannot be perform with these scripts.

## Repository structure

The repository is structured in this way:  
SRV0

- [README.md](http://README.md)  
- 00-raw\_reads  
- 01-pre\_fastq  
- 02-trimming\_filtering  
- 03-reads\_alignment  
- 04-counts\_alignments  
- 05-Differential\_analysis

Each folder’s content is explained later

## Environment installation

When executing the analysis this RNAseq environment must be activated  
In order to export the environment to your computer, run this command on the terminal (requires conda or mamba):  
	```mamba env export \> RNAseq.yml```  
Then activate it: 

	```conda activate RNAseq or source /absolute/path/to/environmentbinaries/activate RNAseq```

## About the environment

The environment has the following packages:

- Sra-toolkit (conda install bioconda::sra-tools) \# We have the symlinks to samples, but first we used it to download the samples  
- Samtools (conda install bioconda::samtools)  
- Fastp (conda install bioconda::fastp)  
- Fastqc (conda install bioconda::fastqc)  
- Multiqc (conda install bioconda::multiqc)  
- FeatureCounts (conda install bioconda::subread)

The tools used were chosen because they have been frequently used in previous analyses (tools such as FastQC are considered the standard to perform analysis like this one), for their robustness and ease of use.

More information about the environment (some of them like trimmomatic were used at first, but may not appear in the script. However, they may be helpful for future analysis or test for other RNAseq analysis):

\# Name                    Version                   Build  Channel  
\_libgcc\_mutex             0.1                 conda\_forge    conda-forge  
\_openmp\_mutex             4.5                       2\_gnu    conda-forge  
alsa-lib                  1.2.13               hb9d3cd8\_0    conda-forge  
annotated-types           0.7.0              pyhd8ed1ab\_1    conda-forge  
brotli-python             1.1.0           py313h46c70d0\_2    conda-forge  
bzip2                     1.0.8                h4bc722e\_7    conda-forge  
c-ares                    1.34.5               hb9d3cd8\_0    conda-forge  
ca-certificates           2025.4.26            hbd8a1cb\_0    conda-forge  
cairo                     1.18.4               h3394656\_0    conda-forge  
certifi                   2025.4.26          pyhd8ed1ab\_0    conda-forge  
cffi                      1.17.1          py313hfab6e84\_0    conda-forge  
charset-normalizer        3.4.1              pyhd8ed1ab\_0    conda-forge  
click                     8.1.8              pyh707e725\_0    conda-forge  
colorama                  0.4.6              pyhd8ed1ab\_1    conda-forge  
coloredlogs               15.0.1             pyhd8ed1ab\_4    conda-forge  
colormath                 3.0.0              pyhd8ed1ab\_4    conda-forge  
curl                      8.13.0               h332b0f4\_0    conda-forge  
expat                     2.6.4                h5888daf\_0    conda-forge  
**fastp                     0.24.1               heae3180\_0    bioconda**  
**fastqc                    0.12.1               hdfd78af\_0    bioconda**  
font-ttf-dejavu-sans-mono 2.37                 hab24e00\_0    conda-forge  
font-ttf-inconsolata      3.000                h77eed37\_0    conda-forge  
font-ttf-source-code-pro  2.038                h77eed37\_0    conda-forge  
font-ttf-ubuntu           0.83                 h77eed37\_3    conda-forge  
fontconfig                2.15.0               h7e30c49\_1    conda-forge  
fonts-conda-ecosystem     1                             0    conda-forge  
fonts-conda-forge         1                             0    conda-forge  
freetype                  2.12.1               h267a509\_2    conda-forge  
giflib                    5.2.2                hd590300\_0    conda-forge  
git                       2.49.0          pl5321h59d505e\_0    conda-forge  
graphite2                 1.3.13            h59595ed\_1003    conda-forge  
h2                        4.2.0              pyhd8ed1ab\_0    conda-forge  
harfbuzz                  10.4.0               h76408a6\_0    conda-forge  
**hisat2                    2.2.1                h503566f\_8    bioconda**  
hpack                     4.1.0              pyhd8ed1ab\_0    conda-forge  
htslib                    1.21                 h566b1c6\_1    bioconda  
humanfriendly             10.0               pyh707e725\_8    conda-forge  
humanize                  4.12.1             pyhd8ed1ab\_0    conda-forge  
hyperframe                6.1.0              pyhd8ed1ab\_0    conda-forge  
icu                       75.1                 he02047a\_0    conda-forge  
idna                      3.10               pyhd8ed1ab\_1    conda-forge  
importlib-metadata        8.6.1              pyha770c72\_0    conda-forge  
isa-l                     2.31.1               hb9d3cd8\_1    conda-forge  
jinja2                    3.1.6              pyhd8ed1ab\_0    conda-forge  
kaleido-core              0.2.1                h3644ca4\_0    conda-forge  
keyutils                  1.6.1                h166bdaf\_0    conda-forge  
krb5                      1.21.3               h659f571\_0    conda-forge  
lcms2                     2.17                 h717163a\_0    conda-forge  
ld\_impl\_linux-64          2.43                 h712a8e2\_4    conda-forge  
lerc                      4.0.0                h27087fc\_0    conda-forge  
libblas                   3.9.0           31\_h59b9bed\_openblas    conda-forge  
libcblas                  3.9.0           31\_he106b2a\_openblas    conda-forge  
libcups                   2.3.3                h4637d8d\_4    conda-forge  
libcurl                   8.13.0               h332b0f4\_0    conda-forge  
libdeflate                1.23                 h4ddbbb0\_0    conda-forge  
libedit                   3.1.20250104    pl5321h7949ede\_0    conda-forge  
libev                     4.33                 hd590300\_2    conda-forge  
libexpat                  2.6.4                h5888daf\_0    conda-forge  
libffi                    3.4.6                h2dba641\_0    conda-forge  
libgcc                    14.2.0               h767d61c\_2    conda-forge  
libgcc-ng                 14.2.0               h69a702a\_2    conda-forge  
libgfortran               14.2.0               h69a702a\_2    conda-forge  
libgfortran5              14.2.0               hf1ad2bd\_2    conda-forge  
libglib                   2.82.2               h2ff4ddf\_1    conda-forge  
libgomp                   14.2.0               h767d61c\_2    conda-forge  
libiconv                  1.18                 h4ce23a2\_1    conda-forge  
libjpeg-turbo             3.0.0                hd590300\_1    conda-forge  
liblapack                 3.9.0           31\_h7ac8fdf\_openblas    conda-forge  
liblzma                   5.6.4                hb9d3cd8\_0    conda-forge  
libmpdec                  4.0.0                h4bc722e\_0    conda-forge  
libnghttp2                1.64.0               h161d5f1\_0    conda-forge  
libopenblas               0.3.29          pthreads\_h94d23a6\_0    conda-forge  
libpng                    1.6.47               h943b412\_0    conda-forge  
libsqlite                 3.49.1               hee588c1\_1    conda-forge  
libssh2                   1.11.1               hcf80075\_0    conda-forge  
libstdcxx                 14.2.0               h8f9b012\_2    conda-forge  
libstdcxx-ng              14.2.0               h4852527\_2    conda-forge  
libtiff                   4.7.0                hd9ff511\_3    conda-forge  
libuuid                   2.38.1               h0b41bf4\_0    conda-forge  
libwebp-base              1.5.0                h851e524\_0    conda-forge  
libxcb                    1.17.0               h8a09558\_0    conda-forge  
libxcrypt                 4.4.36               hd590300\_1    conda-forge  
libxml2                   2.14.0               h8d12d68\_0    conda-forge  
libzlib                   1.3.1                hb9d3cd8\_2    conda-forge  
markdown                  3.6                pyhd8ed1ab\_0    conda-forge  
markdown-it-py            3.0.0              pyhd8ed1ab\_1    conda-forge  
markupsafe                3.0.2           py313h8060acc\_1    conda-forge  
mathjax                   2.7.7                ha770c72\_3    conda-forge  
mdurl                     0.1.2              pyhd8ed1ab\_1    conda-forge  
**multiqc                   1.27.1             pyhdfd78af\_0    bioconda**  
narwhals                  1.30.0             pyhd8ed1ab\_0    conda-forge  
natsort                   8.4.0              pyh29332c3\_1    conda-forge  
ncbi-vdb                  3.2.1                h9948957\_0    bioconda  
ncurses                   6.5                  h2d0b736\_3    conda-forge  
networkx                  3.4.2              pyh267e887\_2    conda-forge  
nspr                      4.36                 h5888daf\_0    conda-forge  
nss                       3.108                h159eef7\_0    conda-forge  
numpy                     2.2.3           py313h17eae1a\_0    conda-forge  
openjdk                   23.0.2               h68779a4\_1    conda-forge  
openjpeg                  2.5.3                h5fbd93e\_0    conda-forge  
openssl                   3.5.0                h7b32b05\_1    conda-forge  
ossuuid                   1.6.2             h5888daf\_1001    conda-forge  
packaging                 24.2               pyhd8ed1ab\_2    conda-forge  
pcre2                     10.44                hba22ea6\_2    conda-forge  
perl                      5.32.1          7\_hd590300\_perl5    conda-forge  
perl-alien-build          2.84            pl5321h7b50bb2\_1    bioconda  
perl-alien-libxml2        0.17            pl5321h577a1d6\_1    bioconda  
perl-business-isbn        3.007           pl5321hdfd78af\_0    bioconda  
perl-business-isbn-data   20210112.006    pl5321hdfd78af\_0    bioconda  
perl-capture-tiny         0.48            pl5321hdfd78af\_2    bioconda  
perl-carp                 1.38            pl5321hdfd78af\_4    bioconda  
perl-constant             1.33            pl5321hdfd78af\_2    bioconda  
perl-data-dumper          2.183           pl5321hec16e2b\_1    bioconda  
perl-encode               3.19            pl5321hec16e2b\_1    bioconda  
perl-exporter             5.72            pl5321hdfd78af\_2    bioconda  
perl-extutils-makemaker   7.70            pl5321hd8ed1ab\_0    conda-forge  
perl-ffi-checklib         0.28            pl5321hdfd78af\_0    bioconda  
perl-file-chdir           0.1010          pl5321hdfd78af\_3    bioconda  
perl-file-path            2.18            pl5321hd8ed1ab\_0    conda-forge  
perl-file-temp            0.2304          pl5321hd8ed1ab\_0    conda-forge  
perl-file-which           1.24            pl5321hd8ed1ab\_0    conda-forge  
perl-importer             0.026           pl5321hdfd78af\_0    bioconda  
perl-mime-base64          3.16            pl5321hec16e2b\_2    bioconda  
perl-parent               0.236           pl5321hdfd78af\_2    bioconda  
perl-path-tiny            0.122           pl5321hdfd78af\_0    bioconda  
perl-pathtools            3.75            pl5321hec16e2b\_3    bioconda  
perl-scope-guard          0.21            pl5321hdfd78af\_3    bioconda  
perl-sub-info             0.002           pl5321hdfd78af\_1    bioconda  
perl-term-table           0.024           pl5321hdfd78af\_0    bioconda  
perl-test2-suite          0.000163        pl5321hdfd78af\_0    bioconda  
perl-uri                  5.12            pl5321hdfd78af\_0    bioconda  
perl-xml-libxml           2.0210          pl5321hf886d80\_0    bioconda  
perl-xml-namespacesupport 1.12            pl5321hdfd78af\_1    bioconda  
perl-xml-sax              1.02            pl5321hdfd78af\_1    bioconda  
perl-xml-sax-base         1.09            pl5321hdfd78af\_1    bioconda  
pillow                    11.1.0          py313h8db990d\_0    conda-forge  
pip                       25.0.1             pyh145f28c\_0    conda-forge  
pixman                    0.44.2               h29eaf8c\_0    conda-forge  
plotly                    6.0.0              pyhd8ed1ab\_0    conda-forge  
pthread-stubs             0.4               hb9d3cd8\_1002    conda-forge  
pyaml-env                 1.2.2              pyhd8ed1ab\_0    conda-forge  
pycparser                 2.22               pyh29332c3\_1    conda-forge  
pydantic                  2.10.6             pyh3cfb1c2\_0    conda-forge  
pydantic-core             2.27.2          py313h920b4c0\_0    conda-forge  
pygments                  2.19.1             pyhd8ed1ab\_0    conda-forge  
pysocks                   1.7.1              pyha55dd90\_7    conda-forge  
python                    3.13.2          hf636f53\_101\_cp313    conda-forge  
python-dotenv             1.0.1              pyhd8ed1ab\_1    conda-forge  
python-kaleido            0.2.1              pyhd8ed1ab\_0    conda-forge  
python\_abi                3.13                    5\_cp313    conda-forge  
pyyaml                    6.0.2           py313h8060acc\_2    conda-forge  
readline                  8.2                  h8c095d6\_2    conda-forge  
regex                     2024.11.6       py313h536fd9c\_0    conda-forge  
requests                  2.32.3             pyhd8ed1ab\_1    conda-forge  
rich                      13.9.4             pyhd8ed1ab\_1    conda-forge  
rich-click                1.8.8              pyhd8ed1ab\_0    conda-forge  
samtools                  1.21                 h96c455f\_1    bioconda  
seqtk                     1.4                  h577a1d6\_3    bioconda  
spectra                   0.0.11             pyhd8ed1ab\_2    conda-forge  
sqlite                    3.49.1               h9eae976\_1    conda-forge  
**sra-tools                 3.2.1                h4304569\_0    bioconda**  
**subread                   2.1.1                h577a1d6\_0    bioconda**  
tiktoken                  0.9.0           py313h3f234b3\_0    conda-forge  
tk                        8.6.13          noxft\_h4845f30\_101    conda-forge  
tqdm                      4.67.1             pyhd8ed1ab\_1    conda-forge  
trimmomatic               0.39                 hdfd78af\_2    bioconda  
typeguard                 4.4.2              pyhd8ed1ab\_0    conda-forge  
typing-extensions         4.12.2               hd8ed1ab\_1    conda-forge  
typing\_extensions         4.12.2             pyha770c72\_1    conda-forge  
tzdata                    2025a                h78e105d\_0    conda-forge  
urllib3                   2.3.0              pyhd8ed1ab\_0    conda-forge  
xorg-libice               1.1.2                hb9d3cd8\_0    conda-forge  
xorg-libsm                1.2.5                he73a12e\_0    conda-forge  
xorg-libx11               1.8.11               h4f16b4b\_0    conda-forge  
xorg-libxau               1.0.12               hb9d3cd8\_0    conda-forge  
xorg-libxdmcp             1.1.5                hb9d3cd8\_0    conda-forge  
xorg-libxext              1.3.6                hb9d3cd8\_0    conda-forge  
xorg-libxfixes            6.0.1                hb9d3cd8\_0    conda-forge  
xorg-libxi                1.8.2                hb9d3cd8\_0    conda-forge  
xorg-libxrandr            1.5.4                hb9d3cd8\_0    conda-forge  
xorg-libxrender           0.9.12               hb9d3cd8\_0    conda-forge  
xorg-libxt                1.3.1                hb9d3cd8\_0    conda-forge  
xorg-libxtst              1.2.5                hb9d3cd8\_3    conda-forge  
yaml                      0.2.5                h7f98852\_2    conda-forge  
zipp                      3.21.0             pyhd8ed1ab\_1    conda-forge  
zlib                      1.3.1                hb9d3cd8\_2    conda-forge  
zstandard                 0.23.0          py313h536fd9c\_1    conda-forge  
zstd                      1.5.7                hb8e6e7a\_1    conda-forge

# Folder specifications

Below, the contents of each folder and what each script does will be explained in detail. However, each script has a more detailed explanation of each command and tool option, as well as the exact commands used to run each one, in the header.
log folders are created in order to make it simpler to evalute whether the script has run correctly and to easily check errors and warning messages.

# 00-raw\_reads

This folder contains two scripts. The first one: [00-data.sh](http://00-data.sh) accepts a .txt file with the accession numbers of the raw sequences of the samples from a sequencing experiment.  
The second script [symlinks.sh](http://symlinks.sh) is the one used in the project which creates symbolic links for .[fastqc.gz](http://fastqc.gz) files. One of this two scripts must be ran in order to star the RNAseq data analysis.

Also it contains two folders:

- logs: will contain log files created by the script  
- results: will contain the symlinks

## Input data:

- Directory where target files are located: ./results

## Output data:

- Symbolic links for [fastq.gz](http://fastq.gz) files named after their accession number  
- A logs folder will be created (unless it’s already created) where output and error logs will be stored. There will be two types of logs:  
  - All samples\_log containing stdout an stderr from the first lines in the code that make sure every file and directory given have the appropriate characteristics  
  - Individual sample logs for each sample

## About the script:

The [symlinks.sh](http://symlinks.sh) script will create symbolic links for those files in this folder.  
**Usage example: ```symlinks.sh -d <absolute/path/to/directory>```**  
Additionally you can write ```symlinks.sh \-h``` or ```symlinks.sh \-v``` to ask for usage help or to know which version the script is respectively.  
First the scripts checks whether there are the correct number of arguments provided  
These arguments are managed using a while-getopts-case loop. Afterwards it will initiate a loop which will read every file in the directory to search for [fastq.gz](http://fastq.gz) files that are not empty and are readable. Finally it will create the symbolic link in the results folder.

# 01-pre\_fastq

This folder contains the script that checks raw sequences’ quality  
Also it contains two folders:

- logs: will contain log files created by the script  
- results: will contain results for the FastQC analysis

## Input data:

- input directory where raw sequences are stored  
- output directory where fastqc and multiqc analysis reports will be stored

## Output data:

- Files of the FastQC analysis: and .html report and a .zip report for each sample  
- multiqc\_analysis.html   
- multiqc\_analysis\_data folder  
- all\_files logs (.err and .out) containing stdout and stderr from the first lines in the code that make sure every file and directory given have the appropriate characteristics  
- Individual sample logs for each sample from the analysis

## About the script:

**Usage example: ```[pre-fastqc.sh](http://pre-fastqc.sh) -i <input_dir> -o <output_dir>```**  
Additionally you can write ```[pre-fastqc.sh](http://pre-fastqc.sh) \-h``` or ```[pre-fastqc.sh](http://pre-fastqc.sh) \-v``` to ask for usage help or to know which version the script is respectively.  
First the script checks whether there are the correct number of arguments provided.  
These arguments are managed using a while-getopts-case loop. Afterwards it will check if the directories provided exist, if not they will be created. Then it will initiate a loop which will read every file in the input directory to search for.[fastq.gz](http://fastq.gz) files that are not empty and are readable. Then it will run a FastQC analysis of each sample and store its report in the results folder. If this analysis has gone correctly it will print a message on screen.  
After all the sequences are analyzed, MultiQC will create an .html report to have all the analysis in one file. It will be stored in the results folder  
This script will create these log files containing stdout and stderr of the general analysis of folder and files at the begining:

- all\_files.out  
- all\_files.err

And these files containing stdout and stderr of the specific analysis for each sample:

- sample.err   
- sample.out

Manual de FastQC: [https://github.com/s-andrews/FastQC/blob/master/fastqc](https://github.com/s-andrews/FastQC/blob/master/fastqc)   
Manual de MultiQC:[https://docs.seqera.io/multiqc/getting\_started/running\_multiqc](https://docs.seqera.io/multiqc/getting_started/running_multiqc) and [https://github.com/MultiQC/MultiQC](https://github.com/MultiQC/MultiQC) 

# 02-trimming\_filtering

This folder contains the script that trims and filters raw sequences. Then it analyzes filtered sequence quality.  
Also it contains two folders:

- logs: will contain log files created by the script  
- results: will contain results for the Fastp, FastQC and MultiQC analysis

## Input data:

Minimum input data (these are the directories we used):

- input directory="../00-raw\_data/results"  
- output directory="./results"  
- fastqc directory="./results/fastqc"  
- multiqc directory="./results/multiqc"

Additional data:

- minimum quality that each base sequencing should have  
- minimum length each sequence should have to be  aligned

if these two parameters are not provided, this script will run with a minimum quality of 20 and a minimum length of 40 bp

## Output data:

- Two \_filtered[.fastq.gz](http://.fastq.gz) files for each sample (paired-end)  
- .html and a .json report for each sample containing the analysis of the filtered sequences created with fastp  
- fastqc folder: contains .html and a [fastq.gz](http://fastq.gz) file for each sequence of the pair of sequences created by fastqc (analysis of filtered sequences quality)  
- multiqc folder: containing the .html report that groups all fastqc reports  
- all\_files logs (.err and .out) containing stdout and stderr from the first lines in the code that make sure every file and directory given have the appropriate characteristics  
- Individual sample logs for each sample from the analysis

## About the script:

**Usage: ```[fastp.sh](http://fastp.sh)  -m <MIN_QUAL> -l <MIN_LEN> -i <INPUT> -o <OUTPUT> -f <FASTQC_DIR> -c <MULTIQC_DIR>```**  
Additionally you can write ```[fastp.sh](http://fastp.sh) \-h``` or ```[fastp.sh](http://fastp.sh) \-v``` to ask for usage help or to know which version the script is respectively.

These options are managed using a while-getopts-case loop. Afterwards it will check if the directories provided exists, if not they will be created. Then it will initiate a loop which will read every file in the input directory to search for.[fastq.gz](http://fastq.gz) files that are not empty and are readable. Then it will filter and trim using the tool fastp to each pair of paired samples.  
This are the options used for fastp:

- i: input files from de input\_dir  
- o: outputfiles that will be stored in the output directory  
- q: minimum quality   
- l: minimum length after trimming  
- detect\_adapter\_for\_pe: Detects adapters and trims them  
- thread: Uses 4 threads  
- html and json: Generates .html and .json reports  
- trim\_poly\_x \--poly\_x\_min\_len 10: Detects polyA and polyT that are minimum 10 bases long and trims them   
- report\_title: Each sample has its report 

If user does not introduce any minimum quality or minimum length, minimum quality will be 20, and minimum length will be 40bp  
Finally a FastQC and MultiQC analysis will be executed

This script will create these log files containing stdout and stderr of the general analysis of folder and files at the begining:

- all\_files.out  
- all\_files.err

And these files containing stdout and stderr of the specific analysis for each sample:

- sample.err   
- sample.out

Link to fastp manual: [https://open.bioqueue.org/home/knowledge/showKnowledge/sig/fastp](https://open.bioqueue.org/home/knowledge/showKnowledge/sig/fastp)  and [https://github.com/OpenGene/fastp/blob/master/README.md](https://github.com/OpenGene/fastp/blob/master/README.md) 

# 03[\-reads\_alignment](https://github.com/SailaraNR/rna-seq-project/tree/main/04-reads_alignment)

This directory contains the script to index and perform the sample’s alignment with the genome.  
Also it contains two folders:

- logs: will contain log files created by the script  
- results: will contain results for Hisat2 and MultiQC analysis

## Input data:

- genome reference in FASTA file  
- genome annotation in GTF file  
- fastq forward and reverse sequences  
- output directory

## Output data:

- For each sample there will be created a folder that contains .bam file, a summary.txt and a sorted.bam file

## About the script:

This scripts index genome reference and align reads to the genome reference.  
**Usage example: ```./alignment_script.sh -d <input_reads_dir> -G <genome/genome.fa> -A <genome/annotation.gtf> -o <out_dir">```**

The bioinformatic tools used are:  
Hisat2 (alignment tool): [https://github.com/DaehwanKimLab/hisat2](https://github.com/DaehwanKimLab/hisat2)   
Samtools (convert .sam to .bam): [https://www.htslib.org/](https://www.htslib.org/)   
MultiQC (alignment analysis): [https://github.com/MultiQC/MultiQC](https://github.com/MultiQC/MultiQC) 

First, the script initializes the logs and checks if .out and .err files exist, if not they will be created. Then checks for the correct number of arguments and uses a getopts for managing options. The script validates the directories and files given and changes permissions if possible, or creates them if they don’t exist. If the user does not introduce the data correctly, the script will end and return a “help” message. Also, this help can be accessed by using \-h option ```./alignment_script.sh -h”```

Indexing genome:  
```
hisat2-build <genome_file> <output_dir>  
Samples alignment:  
-x <genoma_ref> -1 <FW> -2 <RV> -S <output_sam_files> -p <threads> --summary-file <name.txt>  
bam and sorted bam files:  
samtools view -bS <sam_file>  
samtools sort <bam_file> -o <name_sorted_bam>  
Finally, the script performs a MultiQC analysis.  
multiqc <samples_dir> -n <name_multiqc.html> -o <output_dir> --force
```

This script will create these log files containing stdout and stderr of the general analysis of folder and files at the begining:

- all\_files.out  
- all\_files.err

And these files containing stdout and stderr of the specific analysis for each sample:

- sample.err   
- sample.out

# 04[\-counts\_alignments](https://github.com/SailaraNR/rna-seq-project/tree/main/05-counts_alignments)

This directory contains the script to count the numbers of reads of each alignment.   
Also it contains two folders:

- logs: will contain log files created by the script  
- results: will contain results for the Featurecounts

## Input data:

- Genome annotation in GTF file  
- Directory with the sorted bam files

## Output data:

- txt file for each sample with the counts results  
- txt summary file for each sample

## About the script:

This script counts the number of reads of an alignment.  
**Usage example: ```./counts_alignment.sh -A <genome_annotation.gtf> -d <dir_sorted_bam>```**

The bioinformatic tools used is:

- FeatureCounts (reads counts): [https://rnnh.github.io/bioinfo-notebook/docs/featureCounts.html](https://rnnh.github.io/bioinfo-notebook/docs/featureCounts.html) 

First, the script initializes the logs and checks if .out and .err files exist, if not they will be created. Then checks for the correct number of arguments and uses a getopts for managing options. The script validates the directories and files given and changes permissions if possible, or creates them if they don’t exist. If the user does not introduce the data correctly, the script will end and return a “help” message. Also, this help can be accessed by using \-h option ```./counts\_alignments.sh \-h```.

```featureCounts -p -O -T <threads> -a <genome_annotations> -o <name_output_file> <alignment_file>```

This script will create these log files containing stdout and stderr of the general analysis of folder and files at the begining:

- all\_files.out  
- all\_files.err

And these files containing stdout and stderr of the specific analysis for each sample:

- sample.err   
- sample.out

# 05-Differential\_analysis

This directory contains an R script that will execute a differential expression analysis using the table of counts generated by FeatureCounts. Then it will run a Gene Ontology enrichment analysis.  
You will need to run this script using Rstudio so that all of the packages will be installed correctly and the analysis runs correctly.  
It also contains the file metadata.csv that will be needed to run R script. This file is specific for the analysis we have done and possibly not useful for other analysis.

## Input data:

- This script accepts  .txt files containing the table of results created by the analysis of FeatureCounts. For that the user must introduce the absolute path to the folder containing the folder where the FeatureCounts matrixes are stored and the absolute path where the results will be downloaded. These paths should be specified at the beginning of the script. 

## Output data:

It will create two folders:

- DEseq2/Condition\_HighGlucose\_vs\_LowGlucose: it will contain a table containing the analysis of results of DEseq2 analysis comparing High and Low glucose conditions  
- Functional\_Analysis\_GO\_BP.csv: Functional analysis results will be stored in this file   
- plots folder: This folder contains plots generated by the script in a .png format. Should include.  
  - Counts histogram (filtered and unfiltered)  
  - Dendrogram  
  - Heatmap using Pearson Correlation  
  - PCA  
  - MA plot  
  - Volcano plot  
  - GO barplot

Additionally you can download an interactive .html file  and a static pdf file with the results of the analysis by running this command on RStudio console:

```{r download-output, eval=FALSE}  
rmarkdown::render("DEseq2.Rmd", output_format = c("html_document", "pdf_document"))
```

## About the script:

This script installs every package that will be needed to execute the analysis if they are not already installed in the users’ computer. To continue it will read the directory where results will be downloaded and the directory where FeatureCounts matrices are stored (These variables should be filled manually). Then it will process FeatureCount files by creating a matrix especially prepared to be read by DEseq2.  
Later, it will perform a pre-filter of genes before running DEseq2. It will remove genes with low number of counts (\< 10 counts) to reduce noise and improve statistical quality. Then it performs the differential expression analysis using DESeq2 (including normalization, hierarchical clustering, Pearson correlations and PCA plots to help visualize data quality). The comparison was executed comparing HighGlucose and Low glucose conditions and stored in a csv file. Additionally an MA plot and a volcano plot will be plotted to help visualize upregulated, downregulated and significant differences among gene expression. Finally, a GO enrichment analysis will be carried out.

