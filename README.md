# SiteGA
# Requirement
To compile exetubable codes from source codes, you need:
Linux - GCC compuler https://gcc.gnu.org/
Windiws - Microsoft Visual C++, e.g. https://visualstudio.microsoft.com/vs/express/

# Description
Previous SiteGA version required the alignment of binding sites, see of the description of the algorithm here
Levitsky et al. Effective transcription factor binding site prediction using a combination of optimization, a genetic algorithm and discriminant analysis to capture distant interactions. BMC Bioinformatics 8, 481 (2007). https://doi.org/10.1186/1471-2105-8-481
Current SiteGA version represented the algorithm of previous version (2007) adopted for de novo search in ChIP-seq dataset, i.e. the alignment of binding sites is not required

# Source code
Folder src contains five files with sitega source codes in c++ language, they respect to separate modeules of pipeline analysis:
monte0dg.cpp prepare parameter file to train a model (andy02.cpp) or perform the bootsrap cross validation test (andy0bsn2.cpp)
andy02.cpp trains a model with a given train ChIP-seq dataset (peaks)
andy0bsn2.cpp performs the bootsrap cross validation test to estimate the performance of a model with a given train ChIP-seq dataset
andy1_mat.cpp scans a fasta file with DNA sequences with a given model
sitega_thr_dist_mat.cpp creates table of thresholds for the scaner (andy1_mat.cpp) based on score distribution for a background dataset

# How to run separate modules
command line arguments below described for each module

monte0dg.cpp (preparation step to train a model):
1int reg = length of region (default value 6)
2file seq = peaks (fasta file) 
3file out = output file = parameter file from monte0dg.cpp

andy02.cpp (train a model):
1char file_cor = parameter file from monte0dg.cpp
2int motif_len = length of motif (integer value above 30 is recommended, default value 30)
3int size_start = start value for the number of locally positioned dinucleotides (default value 10)
4int size_end = end value for the number of locally positioned dinucleotides (default value 90)
5int size_dif = variation value for the number of locally positioned dinucleotides (default value 10)

andy0bsn2.cpp (performace estimation by cross-validation):
1char file_cor = parameter file from monte0dg.cpp
2int motif_len = length of motif (integer value above 30 is recommended, default value 30)
3int size_start = start value for the number of locally positioned dinucleotides (default value 10)
4int size_end = end value for the number of locally positioned dinucleotides (default value 90)
5int size_dif = variation value for the number of locally positioned dinucleotides (default value 10)
6double ratio_cnt_of_all(0=jk)  = ratio of the number of peaks to the number of control peaks (default value 10)
7int num_iterations = bumber of iteration in bootatrap (default 1), but it is recomended to run bootstrap several times to get reliable results

andy1_mat.cpp (scaner to apply a trained model for a test file)
1file.seq = test file
2sitega_matrix_file = file with sitega model
3file_train = facultative file (default value train.fa)
4thr = threshold for  sitega model
5cmpl = default value 2 
6file.ipr = name for ouput files (default value chipseq)
7seq_head = (default value = 1)
8print_pos = (default value = 1)
9site_desc = (default value = 0)
10bit = (default value = 300)

sitega_thr_dist_mat.cpp (threshold selection for a scaner by false positive rate):
1sitega_matrix_file = file with sitega model
2file_profile_fasta = background dataset (unzip files from folder genomes, use hs* & mm* files for human & mouse, respectively)
3file out_dist = output file, table SiteGA model threshold vs. False Positive Rate (FPR)
4double pvalue_large = low bound for FPR (default value 0.0005)
5double score_min = low bound for tested threshold of SiteGA model (default value 0.997)
6double dpvalue = granulation value for FPR compaction in table (threshold vs. FPR), default value 0.0000000005 implies the absence of compaction

# Interpretation of results

monte0dg.cpp creates file with {mnt} extention https://github.com/parthian-sterlet/sitega/blob/master/examples/diagonal_cov.mnt that may used for training (andy02.cpp) or performance evaluation andy0bsn2.cpp)

andy02.cpp gradually constructs several sitega models, with the numbers of locally positioned dinucleotides (LPDs) assigned in 3rd, 4th and 5th parameters of command line (size_start, size_end and size_dif), their default values 10, 90 and 10 define the search of nine SiteGA models - with 10, 20, 30, etc. up to 90 LPDs. Selection of the final best model among these {10, 20, 30, .., 90} models is performed according to FPR estimated (see file with {train.txt} extension). The final sitega model with the minimal FPR at true positive rate (TPR) 0.5 is written in file with {mat} extention, https://github.com/parthian-sterlet/sitega/blob/master/examples/model

andy0bsn2.cpp may several tomes gradually construct several sitega models (parameter 7th num_iterations), but each time use only a part of dataset for training, the rest part of dataset is used to estimate FPR). Results represent the table of FPRs for TPR 0.1, 0.2, .. up to 0.9. The stored in file with extentsion {bs1.txt} https://github.com/parthian-sterlet/sitega/blob/master/examples/bootstrap.txt

andy1_mat.cpp takes ready sitega model and threshold and construct the profile of hits for tested file in fasta format, main ouput file  respect to 6th parameter of command line, the format of output file is following

https://github.com/parthian-sterlet/sitega/blob/master/examples/hit_profile

after the header of each peak with first '>' symbol from 0 to several line respect to separate hits, for each hit are printed the start position, score, strand and sequence

sitega_thr_dist_mat.cpp compute the distribition of sitega scores, output file (3rd parameter of command line) represents two columns with thresholds and respective FPRs, e.g.

https://github.com/parthian-sterlet/sitega/blob/master/examples/thr_fpr
