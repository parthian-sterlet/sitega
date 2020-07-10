# SiteGA
# Requirement
To compile exetubable codes from source codes, you need:

in Linux system [GCC](https://gcc.gnu.org/) compuler 

in Windiws system [Microsoft Visual C++](https://visualstudio.microsoft.com/vs/express/)

# Description
Previous SiteGA version required the alignment of binding sites, [Levitsky et al., 2007](https://doi.org/10.1186/1471-2105-8-481)

Current SiteGA version represented the algorithm of previous version (2007) adopted for de novo search in a ChIP-seq dataset, i.e. the alignment of binding sites is not required

# Source code
Folder **src** contains five files with sitega source codes in c++ language, they respect to separate modules of pipeline analysis. 
## Module 1: Preparation
monte0dg.cpp prepares parameter file to train a model (andy02.cpp) or perform the bootsrap cross validation test (andy0bsn2.cpp)
## Module 2: Train a model
andy02.cpp  trains a model with a given train ChIP-seq dataset (peaks)
## Module 3: Estimate accuracy for a model
andy0bsn2.cpp performs the bootsrap cross validation test to estimate the performance of a model with a given train ChIP-seq dataset
## Module 4: Set threshold for a model
sitega_thr_dist_mat.cpp creates table of thresholds for the scaner (andy1_mat.cpp) based on score distribution for a background dataset
## Module 5: Scan sequences with a model
andy1_mat.cpp scans a fasta file with DNA sequences with a given model

## Interaction of modules
Pairs of modules **Estimate accuracy for a model** and **Train a model** modules must run with [parameter file](https://github.com/parthian-sterlet/sitega/blob/master/examples/diagonal_cov.mnt) from **Preparation** module

Modules **Set threshold for a model** and **Scan test seauences with a model** require the sitega model which should be previosly computed by **Train a model** module

Module **Set threshold for a model** only helps to select a correct threshold, since scan module **Scan test sequences with a model** module takes the threshold from command line

# How to run separate modules
List of command line arguments for all modules are described below.

## Preparation

monte0dg.cpp 
1. int reg = length of region of one locally positioned dinucleotide (default value 6)
2. file seq = [fasta file of peaks](https://github.com/parthian-sterlet/sitega/blob/master/examples/peaks.fa) each peak should consists of only four types of letters respecting to nucleotides ('a', 'c', 'g' and 't'), i.e. 'n' is forbidden
3. file out = output file = [parameter file](https://github.com/parthian-sterlet/sitega/blob/master/examples/diagonal_cov.mnt) from monte0dg.cpp 

## Train a model

andy02.cpp
1. char file_cor = [parameter file](https://github.com/parthian-sterlet/sitega/blob/master/examples/diagonal_cov.mnt) from monte0dg.cpp 
2. int motif_len = length of motif (default value 30)
3. int size_start = start value for the number of locally positioned dinucleotides (LPDs) (default value 10)
4. int size_end = end value for the number of LPDs (default value 90)
5. int size_dif = variation value for the number of LPDs (default value 10)

## Estimate accuracy for a model

andy0bsn2.cpp
1. char file_cor = [parameter file](https://github.com/parthian-sterlet/sitega/blob/master/examples/diagonal_cov.mnt) from monte0dg.cpp 
2. int motif_len = length of motif (integer value respecting to the optimal length of a traditional position weight matrix is recommended, default value 30 usually brought good results)
3. int size_start = start value for the number of locally positioned dinucleotides (LPDs) (default value 10)
4. int size_end = end value for the number of LPDs (default value 90)
5. int size_dif = variation value for the number of LPDs (default value 10)
6. double ratio_cnt_of_all  = ratio of the number of peaks to the number of control peaks (default value 10)
7. int num_iterations = number of iteration in bootatrap (default 1), but it is recomended to run bootstrap several times (at least 5 runs) to get reliable results

## Set threshold for a model

sitega_thr_dist_mat.cpp
1. sitega_matrix_file = [sitega model file](https://github.com/parthian-sterlet/sitega/blob/master/examples/model.mat)
2. file_profile_fasta = background dataset (unzip files from folder [genomes](https://github.com/parthian-sterlet/sitega/tree/master/genomes), use hs* & mm* files for human & mouse data, respectively)
3. file out_dist = output [Thr vs FPR table file](https://github.com/parthian-sterlet/sitega/blob/master/examples/thr_fpr), table SiteGA model threshold vs. False Positive Rate (FPR)
4. double pvalue_large = maximal FPR (default value 0.0005)
5. double score_min = lowest threshold of SiteGA model (default value 0.997)
6. double dpvalue = granulation value for FPR compaction in table (threshold vs. FPR), default value 0.0000000005 implies the absence of compaction

## Scan test sequences with a model

andy1_mat.cpp
1. file.seq = test file has the same format as [fasta file of peaks], non ('a', 'c', 'g' and 't') nucleotides are ignored
2. sitega_matrix_file = [sitega model file](https://github.com/parthian-sterlet/sitega/blob/master/examples/model.mat)
3. file_train = facultative file (default value train.fa)
4. thr = threshold for sitega model, it is recommended to use [Thr vs FPR table file](https://github.com/parthian-sterlet/sitega/blob/master/examples/thr_fpr) from **set threshold for a model** module
5. cmpl = DNA strand parameter, default value 2 
6. file.ipr = [profile file](https://github.com/parthian-sterlet/sitega/blob/master/examples/hit_profile) with recognized hits
7. seq_head = (default value = 1)
8. print_pos = (default value = 1)
9. site_desc = (default value = 0)
10. bit = (default value = 300)

# Interpretation of results

## Preparation

monte0dg.cpp creates [parameter file](https://github.com/parthian-sterlet/sitega/blob/master/examples/diagonal_cov.mnt) that may be used for training (andy02.cpp) or performance evaluation (andy0bsn2.cpp)

## Train a model

andy02.cpp gradually constructs several sitega models, with the numbers of locally positioned dinucleotides (LPDs) assigned according to 3rd, 4th and 5th parameters of the command line (size_start, size_end and size_dif), their default values 10, 90 and 10 define the search of nine SiteGA models - with 10, 20, 30, etc. up to 90 LPDs. Selection of the final best model among these {10, 20, 30, .., 90} models is performed by FPR estimated (see file with *{train.txt}* extension, it has the same format as output file of bootsrap procedure [FPR_vs TPR table file](https://github.com/parthian-sterlet/sitega/blob/master/examples/model_bs1.txt). The final sitega model with the minimal FPR at true positive rate (TPR) 0.5 is written in [sitega model file](https://github.com/parthian-sterlet/sitega/blob/master/examples/model.mat)

## Estimate accuracy for a model

andy0bsn2.cpp may several tomes gradually construct several sitega models (parameter 7th num_iterations), but each time use only a part of dataset for training, the rest part of dataset is used to estimate FPR). Results represent the table of FPRs for TPR 0.1, 0.2, .. up to 0.9. The stored in [FPR_vs TPR table file](https://github.com/parthian-sterlet/sitega/blob/master/examples/model_bs1.txt)

## Scan test sequences with a model

andy1_mat.cpp takes ready sitega model and threshold and construct the profile of hits for tested file in fasta format, main ouput file respect to 6th parameter of command line, the format of [profile file](https://github.com/parthian-sterlet/sitega/blob/master/examples/hit_profile)
I.e. after the header of each peak with first '>' symbol from zero to several lines respect to separate hits, for each hit a start position, score, strand and whole sequence are printed

## Set threshold for a model

sitega_thr_dist_mat.cpp computes the distribition of sitega scores, output [Thr vs FPR table file](https://github.com/parthian-sterlet/sitega/blob/master/examples/thr_fpr) (3rd parameter of command line) represents two columns with thresholds and respective FPRs
