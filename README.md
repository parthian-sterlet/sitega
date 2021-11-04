# SiteGA - binding Sites recognition by Genetic Algorithm

# Description
The current SiteGA version (Levitsky et al., in preparation) represented the algorithm of [Levitsky et al. (2007)](https://doi.org/10.1186/1471-2105-8-481) adopted for *de novo* motif search in a single ChIP-seq dataset. SiteGA is stochastic algorithm that searches in input sequences the pattern of mutually dependent locally positioned dinucleotides (LPDs) that simulates short- and long-range intractions of nucleotide context within transcription factor binding sites. Hence, SiteGA approach is drastically different from that of traditional Position Weight Matrix (PWM), which searches for the most conserved motifs based on additive impacts of nucleotide frequencies from various site positions.

# Requirements
SiteGA source code is written in C++ language. To compile exetubables from the source code you need:

* In Linux system, C++ compiler, e.g. [GCC](https://gcc.gnu.org/) compiler 
* In Windiws system any VC++ package, e.g. [Microsoft Visual Studio Express](https://visualstudio.microsoft.com/vs/express/)

# Source code
Folder [**src**](https://github.com/parthian-sterlet/sitega/tree/master/src) contains files with SiteGA source codes, they respect to decribed below separate modules of pipeline.
## 1. Background dataset generation
The background (negative) dataset is required as a complement to the foreground (positive, ChIP-seq peaks) dataset to select the two parameters of the SiteGA model: the length of motif and the number of LPDs. The main purpose of the background dataset is to exclude from the results of de novo motif search artifact motifs related to a genome-specific sequence content bias, e.g. polyA. The background sequences notably influence the excpected false positive rate for various motifs in the consequent bootstrap cross-validation procedure (see next parargraph). Any custom sequence set may be used as a background dataset, but it is recommended that this dataset (1) reflected the bias of a genome, i.e. the expected distribution of frequencies of short oligonucleotide (respecting the motif lentghs and (2) was at least several times larger that the foreground dataset. We propose the propgram [background_genome.cpp](https://github.com/parthian-sterlet/sitega/blo/master/src/background_genome.cpp) to generate the background dataset for de novo motif search for certain genome (hg38, mm10, at10). The program generates four output fasta files respecting datasets of genome sequences adopted by mononucleotide content, or additionally applies the relative abundance measures of di-, tri- and tetranucleotides, see [Karlin & Campbell (1994)](https://doi.org/10.1073/pnas.91.26.12842). The analysis of longer oligonucleotide measures implies application of addition criteria on the shorter one, i.e. the dinucleotide measure is applied only to a genomic sequence meeting the criterion on the content of mononicleotides,  the trinucleotide measure is applied only to a sequence meeting the criteria on mono- and dinucleotides, etc. Currently, we recommend the application of only mononucleotide content option for the background dataset generation.
## 2. Set parameters of a model through accuracy estimation
[andy0bsn2.cpp](https://github.com/parthian-sterlet/sitega/blob/master/src/andy0bsn2.cpp) performs the bootsrap cross-validation test to select parameters of a model providing the best performance for given foreground and background datasets. A model has parameters of the optimal length of motif and the number of LPDs. The program consecutively checks various combinations of the motif length and the number of LPDs. The range and step of the motif length are parameters of command line, so that for each length the fixed numbers of LPDs equal to 20, 40, 60 and 80 are checked. The bootsrap cross-validation test denotes the partitioning of the foreground dataset into subsets of training and control sequence sets, the former is used to train a model, while the latter is applied to measure its accuracy. The maximal partial area under curve (pAUC) is used to estimate the accuracy of a model for certain motif length and the number of LPDs. The receiver operating characteristic (ROC) curve with dependence of True Positive Rate (TPR) from False Positive Rate (FPR) allows to compute pAUC. The term partial means that only the part of a ROC curve respecting the criterion FPR < 0.001 is impied for pAUC computation. The background dataset is used to compute FPR for each sequence from the positive control set.
## 3. Train a model
[andy02.cpp](https://github.com/parthian-sterlet/sitega/blob/master/src/andy02.cpp) trains a SiteGA model with selected by the bootstrap cross-validation procedure parameters of the motif length and the number of LPDs for given foreground and background datasets. The resulting model is written in a [special matrix file](https://github.com/parthian-sterlet/sitega/blob/master/examples/model.mat)
## 4. Set threshold for a model
[sitega_thr_dist_mat.cpp](https://github.com/parthian-sterlet/sitega/blob/master/src/sitega_thr_dist_mat.cpp) creates table of thresholds for model's hits search in test sequences (**Scan test sequences with a model** module) based on score distribution for the dataset of whole-genome promoter sequences selected for respective species. Prepeared datasets are stored in [genomes](https://github.com/parthian-sterlet/sitega/tree/master/genomes) folder. The threshold selection implies the estimation of FPR of a model for promoter sequences of whole genome. The dependence of the threshold from FPR is stored in a [special file](https://github.com/parthian-sterlet/sitega/blob/master/examples/thr_fpr)
## 5. Scan test sequences with a model
[andy1_mat.cpp](https://github.com/parthian-sterlet/sitega/blob/master/src/andy1_mat.cpp) scans test sequences with a constructed model and the selected threshold of a model.

# How to compile
* In Linux system: 

git clone https://github.com/parthian-sterlet/sitega \
cd sitega\src\
chmod a+x build.sh\
./build.sh

* In Windiws system:

separate compilation of all source files in VC++

## Integration of modules

Scheme of modules functioning is given below

![scheme](https://github.com/parthian-sterlet/sitega/blob/master/examples/scheme_github_sitega8.jpg)

Modules **Set parameters of a model through accuracy estimation** and **Train a model** should run with file of the background sequence dataset, e.g. it was previously computed by **Background dataset generation** module

Module **Set parameters of a model through accuracy estimation** is required for functionality of **Train a model** and all consequent modules since only the bootstrap procedure correctly selects parameters of a model (see output data block **Table FPR vs. TPR, ROC curve & pAUC**)

Modules **Set threshold for a model** and **Scan test seauences with a model** require file with [SiteGA model](https://github.com/parthian-sterlet/sitega/blob/master/examples/model.mat) which should be previously computed by **Train a model** module

Module **Set threshold for a model** is required to select a correct threshold for **Scan test sequences with a model** module

# How to run separate modules
Lists of command line arguments for all modules are described below

## Background dataset generation

[background_genome.cpp](https://github.com/parthian-sterlet/sitega/blob/master/src/background_genome.cpp)
1. path to whole genome sequences of chromosomes in plain format (see the paragraph below, the last symbol of path must be '/' and '\\' for Linux and Windows OS, respectively)
2. input fasta file
3. output fasta file with genome sequences adopted by mononucleotide content
4. output fasta file with genome sequences adopted by mononucleotide content and dinucleotide measure
5. output fasta file with genome sequences adopted by mononucleotide content, di- and trinucleotide measures
6. output fasta file with genome sequences adopted by mononucleotide content, di-, tri- and tetranucleotide measures
7. maximal number of background sequences per one peak (default value 10)
8. deviation of mononucleitide content of a background sequence from that for a foreground sequence (default value 0.01)
9. percentile threshold for deviation between foreground and background sequences by the dinucleotide measure (default value 1, selection is absent)
10 percentile threshold for deviation between foreground and background sequences by the trinucleotide measure (default value 1, selection is absent)
11. percentile threshold for deviation between foreground and background sequences by the tetranucleotide measure (default value 1, selection is absent)
12. average total number of attemtps to get a background sequence from genome per one foreground sequence (default value 500)
13. genome release (default values are at10, mm10 and hg38 for Arabidopsis, human and mouse genomes)

Whole chromosome sequences in plain format are required to run the program, i.e. headers lines >... should be deleted from the whole chromosome files in fasta format. These plain files should contain only nucleotide letters, IUPAC nucleotides codes N,W,S etc. are ignored by program, all other symbols like ' ', '\t' etc. should deleted, e.g. for Arabidopsis five files are required: chr1.plain, chr2.plain, chr3.plain, chr4.plain, chr5.plain, for human/mouse respective files refer to whole chromosomes 1-22,X,Y / 1-19,X,Y. To see example unzip chr4.plain file from folder [genomes](https://github.com/parthian-sterlet/sitega/tree/master/genomes). Any one of the four output fasta files can be used as the background dataset in consequent analysis (see parameters 2, 3, 4 and 5 of the command line, genome sequences adopted by mononucleotide content, di-, tri-, or tetranucleotide measures, respectively).

## Set parameters of a model through accuracy estimation

[andy0bsn2.cpp](https://github.com/parthian-sterlet/sitega/blob/master/src/andy0bsn2.cpp)
1. path to fasta files with datasets of foreground and background sequences (the last symbol of path must be '/' and '\\' for Linux and Windows OS, respectively)
2. fasta file with dataset of foreground sequences
3. fasta file with dataset of background sequences
4. maximal length of one LPD (default value 6)
5. minimal length of motif (integer value respecting to a tested length L, default value is 8)
6. maximal length of motif (default value is 28)
7. step of length of motif (default value is 4, i.e. lengths 8, 12, 16 etc. are considered)
8. cross-validation type specification: positive value below 1 means the ratio of the training subset size to that of control subset for repeated random subsampling validation, default value -1 means equal sizes of training and control subsets, odd/even peaks are used either for training and control subsets)
9. number of iterations in bootatrap (default 2)

## Train a model

[andy02.cpp](https://github.com/parthian-sterlet/sitega/blob/master/src/andy02.cpp)
1. path to fasta files with datasets of foreground and background sequences (the last symbol of path must be '/' and '\\' for Linux and Windows OS, respectively)
2. fasta file with dataset of foreground sequences
3. fasta file with dataset of background sequences
4. maximal length of one LPD (default value 6)
5. length of motif (integer value respecting to a tested length L, this value is selected by the bootstrap cross-valiation test, see the previous paragraph)
6. size, the number of LPDs (a value is estimated in the bootstrap cross-valiation test, see the previous paragraph)

## Set threshold for a model

[sitega_thr_dist_mat.cpp](https://github.com/parthian-sterlet/sitega/blob/master/src/sitega_thr_dist_mat.cpp)
1. sitega_matrix_file = input file [SiteGA model](https://github.com/parthian-sterlet/sitega/blob/master/examples/model.mat) from **Train a model** module
2. file_profile_fasta = input Whole-genome promoters dataset in fasta format (unzip files from folder [genomes](https://github.com/parthian-sterlet/sitega/tree/master/genomes), use hs*, mm* and at* files for human, mouse and Arabidopsis data, respectively)
3. output file [Table Threshold_vs_FPR](https://github.com/parthian-sterlet/sitega/blob/master/examples/thr_fpr), table SiteGA model threshold vs. False Positive Rate (FPR)
4. pvalue_large = maximal FPR (default value 0.0005)
5. score_min = lowest threshold of SiteGA model (default value 0.995)
6. dpvalue = granulation value for FPR compaction in [Table Threshold_vs_FPR](https://github.com/parthian-sterlet/sitega/blob/master/examples/thr_fpr), the default value 0.0000000005 implies the absence of compaction

## Scan test sequences with a model

[andy1_mat.cpp](https://github.com/parthian-sterlet/sitega/blob/master/src/andy1_mat.cpp)
1. file.seq = input file of test sequences, it has the same format as [Fasta file of peaks](https://github.com/parthian-sterlet/sitega/blob/master/examples/peaks.fa), non ('a', 'c', 'g' and 't') nucleotides are ignored
2. sitega_matrix_file = input file [SiteGA model](https://github.com/parthian-sterlet/sitega/blob/master/examples/model.mat) from **Train a model** module
3. site_description_mode = 0 or 1. 0 means default mode, 1 means computation of frequencies of all LPDs for all tested sequences (option is used for the train fasta file to describe a [SiteGA model](https://github.com/parthian-sterlet/sitega/blob/master/examples/model.mat))
4. FPR threshold = threshold for FPR of [SiteGA model](https://github.com/parthian-sterlet/sitega/blob/master/examples/model.mat) is used to select the SiteGA threshold according to input file [Table Threshold vs FPR](https://github.com/parthian-sterlet/sitega/blob/master/examples/thr_fpr) from **Set threshold for a model** module
5. input file [Threshold vs FPR table](https://github.com/parthian-sterlet/sitega/blob/master/examples/thr_fpr) for [SiteGA model](https://github.com/parthian-sterlet/sitega/blob/master/examples/model.mat) threshold selection by FPR threshold
6. output file [Profile with recognized hits](https://github.com/parthian-sterlet/sitega/blob/master/examples/hit_profile)

# Examples scripts:

These scripts implement various pipelines for Linux:
1. [only training a model](https://github.com/parthian-sterlet/sitega/blob/master/src/train) - **Preparation** and  **Train a model** modules
2. [bootstrap test for a model](https://github.com/parthian-sterlet/sitega/blob/master/src/bootstrap) - **Preparation** and **Estimate accuracy for a model** modules
3. [training and scanning with a model](https://github.com/parthian-sterlet/sitega/blob/master/src/scan) - **Preparation**, **Train a model**, **Set threshold for a model** and **Scan test seauences with a model** modules

# Interpretation of results

## Background dataset generation

[background_genome.cpp](https://github.com/parthian-sterlet/sitega/blob/master/src/background_genome.cpp) prepares the background dataset that may be used for training (**Train a model** module) or performance evaluation (**Set parameters of a model through accuracy estimation** module)

## Set parameters of a model through accuracy estimation

[andy0bsn2.cpp](https://github.com/parthian-sterlet/sitega/blob/master/src/andy0bsn2.cpp) may several times gradually construct several distinct [SiteGA models](https://github.com/parthian-sterlet/sitega/blob/master/examples/model.mat) (parameter of command line 'number of iterations'), but each time it uses only a part of the foreground dataset for training, the rest (control) part of dataset is used to estimate FPR). The output file [Table FPR_vs TPR](https://github.com/parthian-sterlet/sitega/blob/master/examples/model_bs1.txt) represents the table of FPRs for TPR 0.01, 0.02, etc. up to 0.99. Selection of the one model among several ones with different numbers of LPDs and lengths L is performed by respective estimated pAUC values computed for the receiver operating characteristic (ROC) curve (see file with *{train.txt}* extension, [FPR_vs TPR table file](https://github.com/parthian-sterlet/sitega/blob/master/examples/model_bs1.txt). 

## Train a model

[andy02.cpp](https://github.com/parthian-sterlet/sitega/blob/master/src/andy02.cpp) constructs one sitega model, with the numbers of locally positioned dinucleotides (LPDs) assigned according the parameter of the command line, this value deduced from the bootstrap cross validation test (see the previous paragraph). The selected model respecting the maximal pAUC in the bootstrap cross validation test, this model is written in output file [SiteGA model](https://github.com/parthian-sterlet/sitega/blob/master/examples/model.mat) with *.mat* extension

## Set threshold for a model

[sitega_thr_dist_mat.cpp](https://github.com/parthian-sterlet/sitega/blob/master/src/sitega_thr_dist_mat.cpp) computes the distribition of the recognition scores of [SiteGA model](https://github.com/parthian-sterlet/sitega/blob/master/examples/model.mat), output file [Table Threshold vs FPR](https://github.com/parthian-sterlet/sitega/blob/master/examples/thr_fpr) represents two columns with respective thresholds and FPRs

## Scan test sequences with a model

[andy1_mat.cpp](https://github.com/parthian-sterlet/sitega/blob/master/src/andy1_mat.cpp) takes ready sitega model and threshold  from input file [Table Threshold vs FPR](https://github.com/parthian-sterlet/sitega/blob/master/examples/thr_fpr) and constructs the profile of hits for tested file in fasta format, main output file is [Profile file](https://github.com/parthian-sterlet/sitega/blob/master/examples/hit_profile) I.e. after the header of each peak with first '>' symbol from zero to several lines respect to separate hits, for each hit a start position, score, strand and whole sequence are printed. Interpretation of [SiteGA model](https://github.com/parthian-sterlet/sitega/blob/master/examples/model.mat) is performed through computation of LPD frequencies for the training fasta dataset (see fourth argument of the command line). Hence, a computed matrix Number_of_sequences vs. Number_of_LPDs can be used for the correlation analysis, e.g. for revealing the most correlated LPDs in a [SiteGA model](https://github.com/parthian-sterlet/sitega/blob/master/examples/model.mat).
