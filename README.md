# SiteGA
# Requirement
To compile exetubable codes from source codes, you need:
in Linux system GCC compuler https://gcc.gnu.org/
in Windiws system Microsoft Visual C++, e.g. https://visualstudio.microsoft.com/vs/express/

# Description
Previous SiteGA version required the alignment of binding sites (Levitsky et al. (2007) Effective transcription factor binding site prediction using a combination of optimization, a genetic algorithm and discriminant analysis to capture distant interactions. BMC Bioinformatics 8, 481 https://doi.org/10.1186/1471-2105-8-481)

Current SiteGA version represented the algorithm of previous version (2007) adopted for de novo search in ChIP-seq dataset, i.e. the alignment of binding sites is not required

# Source code
Folder **src** contains five files with sitega source codes in c++ language, they respect to separate modules of pipeline analysis. 
## preparation
monte0dg.cpp prepare parameter file to train a model (andy02.cpp) or perform the bootsrap cross validation test (andy0bsn2.cpp)
## train a model
andy02.cpp  trains a model with a given train ChIP-seq dataset (peaks)
## estimate accuracy for a model
andy0bsn2.cpp performs the bootsrap cross validation test to estimate the performance of a model with a given train ChIP-seq dataset
## set threshold for a model
sitega_thr_dist_mat.cpp creates table of thresholds for the scaner (andy1_mat.cpp) based on score distribution for a background dataset
## scan test seauences with a model
andy1_mat.cpp scans a fasta file with DNA sequences with a given model

Among all modules pair of modules {preparation, train a model} and {preparation, estimate accuracy for a model} are tightly related, i.e. the second module each time require the result of the first module. 

{set threshold for a model} and {scan test seauences with a model} modules and require the sitega model which should be previosly computed by {train a model} module

{set threshold for a model} module only helps to select a correct threshold, since scan module {test seauences with a model} takes the threshold from command line

# How to run separate modules
command line arguments below described for each module

## preparation

monte0dg.cpp 

1int reg = length of region of one locally positioned dinucleotide (default value 6)

2file seq = peaks (fasta file) https://github.com/parthian-sterlet/sitega/blob/master/examples/peaks.fa each peak should consists of only four types of letters respecting to nucleotides ('a', 'c', 'g' and 't'), i.e. 'n' is forbidden

3file out = output file = parameter file from monte0dg.cpp https://github.com/parthian-sterlet/sitega/blob/master/examples/diagonal_cov.mnt

## train a model

andy02.cpp

1char file_cor = parameter file from monte0dg.cpp https://github.com/parthian-sterlet/sitega/blob/master/examples/diagonal_cov.mnt

2int motif_len = length of motif (default value 30)

4int size_end = end value for the number of LPDs (default value 90)

5int size_dif = variation value for the number of LPDs (default value 10)

## estimate accuracy for a model

andy0bsn2.cpp

1char file_cor = parameter file from monte0dg.cpp

2int motif_len = length of motif (integer value above 30 is recommended, default value 30)

3int size_start = start value for the number of locally positioned dinucleotides (LPDs) (default value 10)

4int size_end = end value for the number of LPDs (default value 90)

5int size_dif = variation value for the number of LPDs (default value 10)

6double ratio_cnt_of_all  = ratio of the number of peaks to the number of control peaks (default value 10)

7int num_iterations = number of iteration in bootatrap (default 1), but it is recomended to run bootstrap several times (at least 5 runs) to get reliable results

## set threshold for a model

sitega_thr_dist_mat.cpp

1sitega_matrix_file = file with sitega model https://github.com/parthian-sterlet/sitega/blob/master/examples/model.mat

2file_profile_fasta = background dataset (unzip files from folder genomes, use hs* & mm* files for human & mouse, respectively)

3file out_dist = output file, table SiteGA model threshold vs. False Positive Rate (FPR)

4double pvalue_large = low bound for FPR (default value 0.0005)

5double score_min = low bound for tested threshold of SiteGA model (default value 0.997)

6double dpvalue = granulation value for FPR compaction in table (threshold vs. FPR), default value 0.0000000005 implies the absence of compaction

## scan test seauences with a model

andy1_mat.cpp

1file.seq = test file

2sitega_matrix_file = file with sitega model https://github.com/parthian-sterlet/sitega/blob/master/examples/model.mat

3file_train = facultative file (default value train.fa)

4thr = threshold for sitega model

5cmpl = DNA strand parameter, default value 2 

6file.ipr = name for ouput files (default value chipseq)

7seq_head = (default value = 1)

8print_pos = (default value = 1)

9site_desc = (default value = 0)

10bit = (default value = 300)

# Interpretation of results

## preparation

monte0dg.cpp creates file with **mnt** extention https://github.com/parthian-sterlet/sitega/blob/master/examples/diagonal_cov.mnt that may used for training (andy02.cpp) or performance evaluation andy0bsn2.cpp)

## train a model

andy02.cpp gradually constructs several sitega models, with the numbers of locally positioned dinucleotides (LPDs) assigned in 3rd, 4th and 5th parameters of command line (size_start, size_end and size_dif), their default values 10, 90 and 10 define the search of nine SiteGA models - with 10, 20, 30, etc. up to 90 LPDs. Selection of the final best model among these {10, 20, 30, .., 90} models is performed according to FPR estimated (see file with {train.txt} extension). The final sitega model with the minimal FPR at true positive rate (TPR) 0.5 is written in file with **mat** extention, https://github.com/parthian-sterlet/sitega/blob/master/examples/model.mat

## estimate accuracy for a model

andy0bsn2.cpp may several tomes gradually construct several sitega models (parameter 7th num_iterations), but each time use only a part of dataset for training, the rest part of dataset is used to estimate FPR). Results represent the table of FPRs for TPR 0.1, 0.2, .. up to 0.9. The stored in file with extentsion **bs1.txt** https://github.com/parthian-sterlet/sitega/blob/master/examples/model_bs1.txt

## scan test seauences with a model

andy1_mat.cpp takes ready sitega model and threshold and construct the profile of hits for tested file in fasta format, main ouput file respect to 6th parameter of command line, the format of output **profile** file is following https://github.com/parthian-sterlet/sitega/blob/master/examples/hit_profile
I.e. after the header of each peak with first '>' symbol from zero to several lines respect to separate hits, for each hit are printed the start position, score, strand and whole sequence 

## set threshold for a model

sitega_thr_dist_mat.cpp computes the distribition of sitega scores, output **table** file (3rd parameter of command line) represents two columns with thresholds and respective FPRs, e.g. https://github.com/parthian-sterlet/sitega/blob/master/examples/thr_fpr
