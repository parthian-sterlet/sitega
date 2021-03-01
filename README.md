# SiteGA - binding Sites recognition by Genetic Algorithm

# Description
The current SiteGA version (Levitsky et al., in preparation) represented the algorithm of [Levitsky et al. (2007)](https://doi.org/10.1186/1471-2105-8-481) adopted for *de novo* motif search in a ChIP-seq dataset. SiteGA is stochastic algorithm that searches in input sequences the pattern of mutually dependent locally positioned dinucleotides (LPDs) that simulates short- and long-range intractions of nucleotide context within transcription factor binding sites. Hence, SiteGA approach is drastically different from that of traditional Position Weight Matrix (PWM), which searches for the most conserved motifs based on additive impacts of nucleotide frequencies from various site positions.

# Requirements
SiteGA source code was written in C++ language. Hence, to compile exetubables from source code you need:

* In Linux system, C++ compiler, e.g. [GCC](https://gcc.gnu.org/) compiler 
* In Windiws system any VC++ package, e.g. [Microsoft Visual Studio Express](https://visualstudio.microsoft.com/vs/express/)

# Source code
Folder [**src**](https://github.com/parthian-sterlet/sitega/tree/master/src) contains files with SiteGA source codes, they respect to separate modules of pipeline: 
## 1. Preparation
[monte0dg.cpp](https://github.com/parthian-sterlet/sitega/blob/master/src/monte0dg.cpp) prepares [Common settings of models](https://github.com/parthian-sterlet/sitega/blob/master/examples/diagonal_cov.mnt) that are required to perform the bootsrap cross validation test (**Estimate accuracy for a model** module) and to train a model (**Train a model** module). 'Common settings of models' are diaganal elements of the covariation matrix for LPDs of all dinucleotide types and all allowed lengths for the background dataset [(Levitsky et al. 2007)](https://doi.org/10.1186/1471-2105-8-481)
## 2. Set parameters of a model through accuracy estimation
[andy0bsn2.cpp](https://github.com/parthian-sterlet/sitega/blob/master/src/andy0bsn2.cpp) performs the bootsrap cross-validation test to select parameters of model: the optimial length and the number of LPDs providing the best performance. The maximal partial area under curve (pAUC) is used to estimate the accuracy of a model, i.e. the receiver operating characteristic (ROC) curve with dependence of True Positive Rate (TPR) from False Positive Rate (FPR) for control data. The term partial means that the criterion FPR < 0.001 is impied for pAUC computation. The maximal pAUC estimates the performance of a model with a given train ChIP-seq dataset.
## 3. Train a model
[andy02.cpp](https://github.com/parthian-sterlet/sitega/blob/master/src/andy02.cpp) trains a [SiteGA model](https://github.com/parthian-sterlet/sitega/blob/master/examples/model.mat) with selected parameters of length and number of LPDs for a given train ChIP-seq dataset (peaks)
## 4. Set threshold for a model
[sitega_thr_dist_mat.cpp](https://github.com/parthian-sterlet/sitega/blob/master/src/sitega_thr_dist_mat.cpp) creates table of thresholds for the scaner (**Scan test sequences with a model** module) based on score distribution for the background set of whole whole-genome promoters
## 5. Scan test sequences with a model
[andy1_mat.cpp](https://github.com/parthian-sterlet/sitega/blob/master/src/andy1_mat.cpp) scans test sequences with a constructed model and a selected for it threshold

# How to compile
* In Linux system: 

git clone https://github.com/parthian-sterlet/sitega \
cd sitega\src\
chmod a+x build.sh\
./build.sh

* In Windiws system:

separate compilation of five modules in VC++

## Integration of modules

Scheme of modules fucntioning is given below

![scheme](https://github.com/parthian-sterlet/sitega/blob/master/examples/scheme_github_sitega3.jpg)

Modules **Estimate accuracy for a model** and **Train a model** must run with file of [Common settings of models](https://github.com/parthian-sterlet/sitega/blob/master/examples/diagonal_cov.mnt) which previously computed by **Preparation** module

Modules **Set threshold for a model** and **Scan test seauences with a model** require file with [SiteGA model](https://github.com/parthian-sterlet/sitega/blob/master/examples/model.mat) which should be previosly computed by **Train a model** module

Module **Set threshold for a model** is required to select a correct threshold for **Scan test sequences with a model** module

Module **Estimate accuracy for a model** is not required for functionality of **Set threshold for a model** and **Scan test seauences with a model** modules. Though only bootstrap procedure correctly evaluates the accuracy (see output data block **ROC curve, control data**), estimates of False Positive Rate for a [SiteGA model](https://github.com/parthian-sterlet/sitega/blob/master/examples/model.mat) may be retrieved from results of testing with training data (**Train a model**), i.e. see output file with *{train.txt}* extension [FPR_vs_TPR table file](https://github.com/parthian-sterlet/sitega/blob/master/examples/model_bs1.txt), it respects to output data block **ROC curve, training data** on the scheme

# How to run separate modules
Lists of command line arguments for all modules are described below

## Preparation

[monte0dg.cpp](https://github.com/parthian-sterlet/sitega/blob/master/src/monte0dg.cpp)
1. int reg = maximal length of region of one locally positioned dinucleotide (default value 6)
2. file seq = input [Fasta file of peaks](https://github.com/parthian-sterlet/sitega/blob/master/examples/peaks.fa), each peak should consist of only four types of letters respecting to nucleotides ('a', 'c', 'g' and 't'), i.e. 'n' is forbidden
3. file out = output file of [Common settings of models](https://github.com/parthian-sterlet/sitega/blob/master/examples/diagonal_cov.mnt)

## Set parameters of a model through accuracy estimation

[andy0bsn2.cpp](https://github.com/parthian-sterlet/sitega/blob/master/src/andy0bsn2.cpp)
1. char file_cor = input file of [Common settings of models](https://github.com/parthian-sterlet/sitega/blob/master/examples/diagonal_cov.mnt) from **Preparation** module
2. int motif_len = length of motif (integer value respecting to the optimal length L, default value is optimized in the range from 40 to 100 bp)
3. int size_start = start value for the number of LPDs (default value the optimal length of motif, L)
4. int size_end = end value for the number of LPDs (default value twice larger than the optimal length of motif, 2L)
5. int size_dif = variation value for the number of LPDs (default value a quarter of the optimal length, L/4)
6. double ratio_cnt_of_all  = ratio of the number of peaks in training dataset to that control dataset (default value -1 means odd/even peaks for training/control sequence sets)
7. int num_iterations = number of iterations in bootatrap (default 1)

## Train a model

[andy02.cpp](https://github.com/parthian-sterlet/sitega/blob/master/src/andy02.cpp)
1. char file_cor = input file of [Common settings of models](https://github.com/parthian-sterlet/sitega/blob/master/examples/diagonal_cov.mnt) from **Preparation** module
2. int motif_len = length of motif (integer value L respecting to the optimal length of a respective traditional PWM model is recommended, default values from 50 to 100 are recommended, the best length is selected according the accuracy estimate pAUC, the partial area under curve, which is estimated in the bootstrap crossvaliation test, see the next paragraph)
3. int size_start = start value for the number of LPDs (a value is estimated in the bootstrap crossvaliation test, see the next paragraph)
4. int size_end = end value for the number of LPDs (default value is equal to the previous parameter size_start)
5. int size_dif = variation value for the number of LPDs (default value 10)

## Set threshold for a model

[sitega_thr_dist_mat.cpp](https://github.com/parthian-sterlet/sitega/blob/master/src/sitega_thr_dist_mat.cpp)
1. sitega_matrix_file = input file [SiteGA model](https://github.com/parthian-sterlet/sitega/blob/master/examples/model.mat) from **Train a model** module
2. file_profile_fasta = input background dataset in fasta format (unzip files from folder [genomes](https://github.com/parthian-sterlet/sitega/tree/master/genomes), use hs* & mm* files for human & mouse data, respectively)
3. char output file [Table Threshold_vs_FPR](https://github.com/parthian-sterlet/sitega/blob/master/examples/thr_fpr), table SiteGA model threshold vs. False Positive Rate (FPR)
4. double pvalue_large = maximal FPR (default value 0.0005)
5. double score_min = lowest threshold of SiteGA model (default value 0.995)
6. double dpvalue = granulation value for FPR compaction in [Table Threshold_vs_FPR](https://github.com/parthian-sterlet/sitega/blob/master/examples/thr_fpr), the default value 0.0000000005 implies the absence of compaction

## Scan test sequences with a model

[andy1_mat.cpp](https://github.com/parthian-sterlet/sitega/blob/master/src/andy1_mat.cpp)
1. file.seq = input file of test sequences, it has the same format as [Fasta file of peaks](https://github.com/parthian-sterlet/sitega/blob/master/examples/peaks.fa), non ('a', 'c', 'g' and 't') nucleotides are ignored
2. char sitega_matrix_file = input file [SiteGA model](https://github.com/parthian-sterlet/sitega/blob/master/examples/model.mat) from **Train a model** module
3. int site_description_mode = 0 or 1. 0 means default mode, 1 means computation of frequencies of all LPDs for all tested sequences (option is used for the train fasta file to describe a [SiteGA model](https://github.com/parthian-sterlet/sitega/blob/master/examples/model.mat))
4. double FPR threshold = threshold for FPR of [SiteGA model](https://github.com/parthian-sterlet/sitega/blob/master/examples/model.mat) is used to select the SiteGA threshold according to input file [Table Threshold vs FPR](https://github.com/parthian-sterlet/sitega/blob/master/examples/thr_fpr) from **Set threshold for a model** module
5. input file [Threshold vs FPR table](https://github.com/parthian-sterlet/sitega/blob/master/examples/thr_fpr) for [SiteGA model](https://github.com/parthian-sterlet/sitega/blob/master/examples/model.mat) threshold selection by FPR threshold
6. output file [Profile with recognized hits](https://github.com/parthian-sterlet/sitega/blob/master/examples/hit_profile)

# Examples scripts:

These scripts implement various pipelines for Linux:
1. [only training a model](https://github.com/parthian-sterlet/sitega/blob/master/src/train) - **Preparation** and  **Train a model** modules
2. [bootstrap test for a model](https://github.com/parthian-sterlet/sitega/blob/master/src/bootstrap) - **Preparation** and **Estimate accuracy for a model** modules
3. [training and scanning with a model](https://github.com/parthian-sterlet/sitega/blob/master/src/scan) - **Preparation**, **Train a model**, **Set threshold for a model** and **Scan test seauences with a model** modules

# Interpretation of results

## Preparation

[monte0dg.cpp](https://github.com/parthian-sterlet/sitega/blob/master/src/monte0dg.cpp) creates file of [Common settings of models](https://github.com/parthian-sterlet/sitega/blob/master/examples/diagonal_cov.mnt) that may be used for training (**Train a model** module) or performance evaluation (**Estimate accuracy for a model** module)

## Set parameters of a model through accuracy estimation

[andy0bsn2.cpp](https://github.com/parthian-sterlet/sitega/blob/master/src/andy0bsn2.cpp) may several times gradually construct several distinct [SiteGA models](https://github.com/parthian-sterlet/sitega/blob/master/examples/model.mat) (parameter 7th num_iterations), but each time it uses only a part of ChIP-seq dataset for training, the rest (control) part of dataset is used to estimate FPR). The output file [Table FPR_vs TPR](https://github.com/parthian-sterlet/sitega/blob/master/examples/model_bs1.txt) represents the table of FPRs for TPR 0.01, 0.02, .. up to 0.99. Selection of the one model among several ones with different LPDs and lengths L is performed by respective estimated pAUC values computed for the receiver operating characteristic (ROC) curve (see file with *{train.txt}* extension, [FPR_vs TPR table file](https://github.com/parthian-sterlet/sitega/blob/master/examples/model_bs1.txt). 

## Train a model

[andy02.cpp](https://github.com/parthian-sterlet/sitega/blob/master/src/andy02.cpp) constructs one sitega model, with the numbers of locally positioned dinucleotides (LPDs) assigned according to 3rd, 4th and 5th parameters of the command line (size_start, size_end and size_dif), their values deduced from the bootstrap cross validation test (see the next paragraph). The selected model respecting the maximal pAUC in the bootstrap cross validation test, this model is written in output file [SiteGA model](https://github.com/parthian-sterlet/sitega/blob/master/examples/model.mat) with *.mat* extension

## Set threshold for a model

[sitega_thr_dist_mat.cpp](https://github.com/parthian-sterlet/sitega/blob/master/src/sitega_thr_dist_mat.cpp) computes the distribition of the recognition scores of [SiteGA model](https://github.com/parthian-sterlet/sitega/blob/master/examples/model.mat), output file [Table Threshold vs FPR](https://github.com/parthian-sterlet/sitega/blob/master/examples/thr_fpr) represents two columns with respective thresholds and FPRs

## Scan test sequences with a model

[andy1_mat.cpp](https://github.com/parthian-sterlet/sitega/blob/master/src/andy1_mat.cpp) takes ready sitega model and threshold  from input file [Table Threshold vs FPR](https://github.com/parthian-sterlet/sitega/blob/master/examples/thr_fpr) and constructs the profile of hits for tested file in fasta format, main output file is [Profile file](https://github.com/parthian-sterlet/sitega/blob/master/examples/hit_profile) I.e. after the header of each peak with first '>' symbol from zero to several lines respect to separate hits, for each hit a start position, score, strand and whole sequence are printed. Interpretation of [SiteGA model](https://github.com/parthian-sterlet/sitega/blob/master/examples/model.mat) is performed through computation of LPD frequencies for the training fasta dataset (see fourth argument of the command line). Hence, a computed matrix Number_of_sequences vs. Number_of_LPDs can be used for the correlation analysis, e.g. for revealing the most correlated LPDs in a [SiteGA model](https://github.com/parthian-sterlet/sitega/blob/master/examples/model.mat).
