# SiteGA
# Requirement
To compile exetubable codes from source code, that was written in C++ language, you need:

in Linux system [GCC](https://gcc.gnu.org/) compuler 

in Windiws system [Microsoft Visual C++](https://visualstudio.microsoft.com/vs/express/)

# Description
Previous SiteGA version required the alignment of binding sites, [Levitsky et al., 2007](https://doi.org/10.1186/1471-2105-8-481)

Current SiteGA version represented the algorithm of previous version (2007) adopted for de novo search in a ChIP-seq dataset, i.e. the alignment of binding sites is not required

# Source code
Folder [**src**](https://github.com/parthian-sterlet/sitega/tree/master/src) contains five files with SiteGA source codes, they respect to five separate modules of pipeline: 
## 1. Preparation
monte0dg.cpp prepares parameter file to train a model (andy02.cpp) or perform the bootsrap cross validation test (andy0bsn2.cpp)
## 2. Train a model
andy02.cpp  trains a model with a given train ChIP-seq dataset (peaks)
## 3. Estimate accuracy for a model
andy0bsn2.cpp performs the bootsrap cross-validation test to estimate the performance of a model with a given train ChIP-seq dataset
## 4. Set threshold for a model
sitega_thr_dist_mat.cpp creates table of thresholds for the scaner (andy1_mat.cpp) based on score distribution for a background dataset
## 5. Scan test sequences with a model
andy1_mat.cpp scans a fasta file with DNA sequences with a given model

## Integration of modules
![scheme](https://github.com/parthian-sterlet/sitega/blob/master/examples/scheme_github_sitega.jpg)

Pairs of modules **Estimate accuracy for a model** and **Train a model** modules must run with [parameter file](https://github.com/parthian-sterlet/sitega/blob/master/examples/diagonal_cov.mnt) which previously computed by **Preparation** module. See examples scripts for [only training a model](https://github.com/parthian-sterlet/sitega/blob/master/scripts/train) and [bootstrap test](https://github.com/parthian-sterlet/sitega/blob/master/scripts/bootstrap)

Modules **Set threshold for a model** and **Scan test seauences with a model** require [sitega model file](https://github.com/parthian-sterlet/sitega/blob/master/examples/model.mat) which should be previosly computed by **Train a model** module

Module **Set threshold for a model** is required to select a correct threshold for **Scan test sequences with a model** module, see example script [for training and scanning](https://github.com/parthian-sterlet/sitega/blob/master/scripts/scan)

Module **Estimate accuracy for a model** is not required for functionality of **Set threshold for a model** and **Scan test seauences with a model** modules. Though only bootstrap procedure correctly evaluates the accuracy (see block **ROC curve, control data**), estimates of False Positive Rate for a [sitega model file](https://github.com/parthian-sterlet/sitega/blob/master/examples/model.mat) may be retrieved from results of testing with training data (**Train a model**), i.e. see output file with *{train.txt}* extension [FPR_vs_TPR table file](https://github.com/parthian-sterlet/sitega/blob/master/examples/model_bs1.txt), it respects to **ROC curve, training data** block on the scheme

# How to run separate modules
List of command line arguments for all modules are described below.

## Preparation

monte0dg.cpp 
1. int reg = length of region of one locally positioned dinucleotide (default value 6)
2. file seq = [fasta file of peaks](https://github.com/parthian-sterlet/sitega/blob/master/examples/peaks.fa) each peak should consists of only four types of letters respecting to nucleotides ('a', 'c', 'g' and 't'), i.e. 'n' is forbidden
3. file out = output [parameter file](https://github.com/parthian-sterlet/sitega/blob/master/examples/diagonal_cov.mnt)

## Train a model

andy02.cpp
1. char file_cor = input [parameter file](https://github.com/parthian-sterlet/sitega/blob/master/examples/diagonal_cov.mnt) from monte0dg.cpp 
2. int motif_len = length of motif (default value 30)
3. int size_start = start value for the number of locally positioned dinucleotides (LPDs) (default value 10)
4. int size_end = end value for the number of LPDs (default value 90)
5. int size_dif = variation value for the number of LPDs (default value 10)

## Estimate accuracy for a model

andy0bsn2.cpp
1. char file_cor = input [parameter file](https://github.com/parthian-sterlet/sitega/blob/master/examples/diagonal_cov.mnt) from monte0dg.cpp 
2. int motif_len = length of motif (integer value respecting to the optimal length of a traditional position weight matrix is recommended, default value 30 usually brought good results)
3. int size_start = start value for the number of locally positioned dinucleotides (LPDs) (default value 10)
4. int size_end = end value for the number of LPDs (default value 90)
5. int size_dif = variation value for the number of LPDs (default value 10)
6. double ratio_cnt_of_all  = ratio of the number of peaks to the number of control peaks (default value 10)
7. int num_iterations = number of iteration in bootatrap (default 1), but it is recomended to run bootstrap several times (at least 5 runs) to get reliable results

## Set threshold for a model

sitega_thr_dist_mat.cpp
1. sitega_matrix_file = input [sitega model file](https://github.com/parthian-sterlet/sitega/blob/master/examples/model.mat) from andy02.cpp
2. file_profile_fasta = background dataset (unzip files from folder [genomes](https://github.com/parthian-sterlet/sitega/tree/master/genomes), use hs* & mm* files for human & mouse data, respectively)
3. output [Thr_vs_FPR table file](https://github.com/parthian-sterlet/sitega/blob/master/examples/thr_fpr), table SiteGA model threshold vs. False Positive Rate (FPR)
4. double pvalue_large = maximal FPR (default value 0.0005)
5. double score_min = lowest threshold of SiteGA model (default value 0.997)
6. double dpvalue = granulation value for FPR compaction in table (Threshold vs. FPR), default value 0.0000000005 implies the absence of compaction

## Scan test sequences with a model

andy1_mat.cpp
1. file.seq = input test file has the same format as [fasta file of peaks], non ('a', 'c', 'g' and 't') nucleotides are ignored
2. sitega_matrix_file = input [sitega model file](https://github.com/parthian-sterlet/sitega/blob/master/examples/model.mat) from andy02.cpp
4. FPR threshold = threshold for FPR of sitega model is used to select the sitega threshold according to [Thr vs FPR table file](https://github.com/parthian-sterlet/sitega/blob/master/examples/thr_fpr) from **Set threshold for a model** module
5. input [Thr vs FPR table file](https://github.com/parthian-sterlet/sitega/blob/master/examples/thr_fpr), table SiteGA model threshold vs. False Positive Rate (FPR)
6. output [profile file](https://github.com/parthian-sterlet/sitega/blob/master/examples/hit_profile) with recognized hits

# Interpretation of results

## Preparation

monte0dg.cpp creates [parameter file](https://github.com/parthian-sterlet/sitega/blob/master/examples/diagonal_cov.mnt) that may be used for training (andy02.cpp) or performance evaluation (andy0bsn2.cpp)

## Train a model

andy02.cpp gradually constructs several sitega models, with the numbers of locally positioned dinucleotides (LPDs) assigned according to 3rd, 4th and 5th parameters of the command line (size_start, size_end and size_dif), their default values 10, 90 and 10 define the search of nine SiteGA models - with 10, 20, 30, etc. up to 90 LPDs. Selection of the final best model among these {10, 20, 30, .., 90} models is performed by FPR estimated (see file with *{train.txt}* extension, it has the same format as output file of accuracy estimation procedure [FPR_vs TPR table file](https://github.com/parthian-sterlet/sitega/blob/master/examples/model_bs1.txt). The final sitega model with the minimal FPR at true positive rate (TPR) 0.5 is written in [sitega model file](https://github.com/parthian-sterlet/sitega/blob/master/examples/model.mat)

## Estimate accuracy for a model

andy0bsn2.cpp may several times gradually construct several sitega models (parameter 7th num_iterations), but each time use only a part of dataset for training, the rest part of dataset is used to estimate FPR). Output [FPR_vs TPR table file](https://github.com/parthian-sterlet/sitega/blob/master/examples/model_bs1.txt) represents the table of FPRs for TPR 0.1, 0.2, .. up to 0.9. 

## Set threshold for a model

sitega_thr_dist_mat.cpp computes the distribition of SiteGA scores, output [Thr vs FPR table file](https://github.com/parthian-sterlet/sitega/blob/master/examples/thr_fpr) represents two columns with thresholds and respective FPRs

## Scan test sequences with a model

andy1_mat.cpp takes ready sitega model and threshold  from input  [Thr vs FPR table file](https://github.com/parthian-sterlet/sitega/blob/master/examples/thr_fpr) and constructs the profile of hits for tested file in fasta format, main output file is [profile file](https://github.com/parthian-sterlet/sitega/blob/master/examples/hit_profile)
I.e. after the header of each peak with first '>' symbol from zero to several lines respect to separate hits, for each hit a start position, score, strand and whole sequence are printed
