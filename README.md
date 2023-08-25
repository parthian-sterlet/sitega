# SiteGA - binding Sites recognition by Genetic Algorithm

# Description
SiteGA implements the earlier algorithm [Levitsky et al. (2007)](https://doi.org/10.1186/1471-2105-8-481) for *de novo* motif search in a single ChIP-seq dataset (ChIP-seq peaks) [(Tsukanov et al., 2022)](https://doi.org/10.3389/fpls.2022.938545). SiteGA applies Genetic Algorithm that searches in input sequences the motif as a pattern of mutually dependent locally positioned dinucleotides (LPDs) that simulates short- and long-range interactions of nucleotide context within transcription factor binding sites. Hence, SiteGA approach is drastically different from that of the traditional model of Position Weight Matrix (PWM), which searches for the motif consisting of additive impacts of nucleotide frequencies from various site positions. 

SiteGA estimates the alignment X = {x(1), x(2), ..., x(N)} of N best predicted sites in N peaks with the fitness function F(X) = D(X) * E(X). The first factor D(X) reflects dependencies between positions within an alignment through the linear discriminant analysis approach as in the original paper [(Levitsky et al. 2007)](https://doi.org/10.1186/1471-2105-8-481), while the second factor E(X) implies the average enrichment of k-mers within an alignment of predicted sites [(Tsukanov et al., 2022)](https://doi.org/10.3389/fpls.2022.938545). The program indexed all possible positions of sites in peaks (foreground set) with the log-ratio enrichment of the inherent k-mers in the comparison of foreground and background sequences sets. The average log-ratio enrichment for k-mers located within a whole alignment defines the second factor E(X). 
DepLogo tool [(Grau et al. 2019)](https://doi.org/10.1093/bioinformatics/btz507) can be used to visualize the nucleotide context pattern explaining the recognition of binding sites by the SiteGA model, the required alignment of predicted sites can be deduced from [one](https://github.com/parthian-sterlet/sitega/blob/master/examples/hit_profile) of output files.

# Requirements
SiteGA source code is written in C++ language. To compile exetubables from the source code you need:

* In Linux system, C++ compiler, e.g. [GCC](https://gcc.gnu.org/) compiler 
* In Windiws system any VC++ package, e.g. [Microsoft Visual Studio Community](https://visualstudio.microsoft.com/vs/community/)

# Input data
The main part of input data represent ChIP-seq peaks in a Fasta format. To optimize the computation time, the restriction 3000 bp for the length of any peak is used, though a length of peak is not restricted by an algorithm. The same reason requires the application in analysis of at least moderate number of peaks, typically ~1000 peaks are enough to derive a motif.

# Source code
Folder [**src**](https://github.com/parthian-sterlet/sitega/tree/master/src) contains files with SiteGA source codes, they respect to decribed below separate modules of pipeline.

## 1. Background sequence set generation
The background (negative) set of sequences is required as a complement to the foreground (positive) set to select the two parameters of the SiteGA model: the length of motif and the number of LPDs. The main purpose of the background set is to exclude from the results of *de novo* motif search artifact motifs related to a genome-specific sequence content bias, e.g. polyA. The background sequences notably influence the excpected false positive rate for various motifs in the consequent bootstrap cross-validation procedure (see the next parargraph). Any custom sequence set may be used as a background set, but it is recommended that this set (1) reflects the bias of a genome, i.e. the expected distribution of frequencies of nucleotides and (2) is at least several times larger that the foreground set. We propose the propgram [background_genome_mono.cpp](https://github.com/parthian-sterlet/antinoise/blob/main/src/background_genome_mono.cpp) to generate the background set for de novo motif search for certain genome (hg38, mm10, at10). This program generates an output fasta file respecting the foreground set adopted by the content of A/T nucleotides. The command line tool [AntiNoise](https://github.com/parthian-sterlet/antinoise) and web service [AntiNoise](https://denovosea.icgbio.ru/antinoise/) are specially devoted to this genomic approach of background sequences generation.
## 2. Set parameters of a model through accuracy estimation
[andy0bsn5.cpp](https://github.com/parthian-sterlet/sitega/blob/master/src/andy0bsn5.cpp) performs the bootsrap cross-validation test to select parameters of a model providing the best performance for given foreground and background sets. A model has parameters of the optimal length of motif and the number of LPDs. The program consecutively checks various combinations of the motif length and the number of LPDs. The range and step of the motif length are parameters of command line, so that for each length the fixed numbers of LPDs equal to 40, 60, 80 and 100 are checked. The bootsrap cross-validation test denotes the partitioning of the foreground set into subsets of training and control sequence sets, the former is used to train a model, while the latter is applied to measure its accuracy. The maximal partial area under curve (pAUC) and area under precision-recall curve (AUC-PR) are used to estimate the accuracy of a model for certain motif length and the number of LPDs. The receiver operating characteristic (ROC) curve with dependence of True Positive Rate (TPR) from False Positive Rate (FPR) allows to compute pAUC, while the dependence of Precision (or Positive Predictive Value, PPV) TP/(TP+FP) from Recall (TPR) TP/(TP+FN) represents AUC-PR (FP & FN are counts of false positives & false negatives). The term partial means that only the part of a ROC curve respecting the criterion FPR < 0.001 is impied for pAUC computation. We estimate FPR for ROC curve as the frequency of motif occurrences, but not as a fraction of background sequences. TPR (Recall) for ROC (AUC-PR), PPV for AUC-PR imply a fraction of sequences.
## 3. Train a model
[andy05.cpp](https://github.com/parthian-sterlet/sitega/blob/master/src/andy05.cpp) trains a SiteGA model with selected by the bootstrap cross-validation procedure parameters of the motif length and the number of LPDs for given foreground and background sets. The resulting model is written in a [special matrix file](https://github.com/parthian-sterlet/sitega/blob/master/examples/model.mat)
## 4. Set threshold for a model
[sitega_thr_dist_mat.cpp](https://github.com/parthian-sterlet/sitega/blob/master/src/sitega_thr_dist_mat.cpp) creates table of thresholds for model's hits search in test sequences (**Scan test sequences with a model** module) based on score distribution for the set of whole-genome promoter sequences selected for respective species. Prepeared sets are stored in [genomes](https://github.com/parthian-sterlet/sitega/tree/master/genomes) folder. The threshold selection implies the estimation of Expected Recognition Rate (ERR) of a model for promoter sequences of whole genome. The dependence of the threshold from ERR is stored in a [special file](https://github.com/parthian-sterlet/sitega/blob/master/examples/thr_err)
## 5. Scan test sequences with a model
[andy1_mat.cpp](https://github.com/parthian-sterlet/sitega/blob/master/src/andy1_mat.cpp) scans test sequences with a constructed model and the selected threshold of a model.
## 6. Scan whole genome with a model
[andy1_mat_long.cpp](https://github.com/parthian-sterlet/sitega/blob/master/src/andy1_mat_long.cpp) scans whole genome as a set of full-sized chromosomes with a constructed model and the selected threshold of a model.

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

![scheme](https://github.com/parthian-sterlet/sitega/blob/master/examples/sitega_github.png)

Modules **Set parameters of a model through accuracy estimation** and **Train a model** should run with file of the background sequence set, e.g. it was previously computed by **Background set generation** module

Module **Set parameters of a model through accuracy estimation** is required for functionality of **Train a model** and all consequent modules since only the bootstrap procedure correctly selects parameters of a model (see output data block **Table FPR vs. TPR, ROC curve & pAUC**)

Modules **Set threshold for a model** and **Scan test seauences with a model** require file with [SiteGA model](https://github.com/parthian-sterlet/sitega/blob/master/examples/model.mat) which should be previously computed by **Train a model** module

Module **Set threshold for a model** is required to select a correct threshold for **Scan test sequences with a model** module

# How to run separate modules
Lists of command line arguments for all modules are described below

## Background set generation

Use [background_genome_mono.cpp](https://github.com/parthian-sterlet/antinoise/blob/main/src/background_genome_mono.cpp), 
see the github repositiory [AntiNoise](https://github.com/parthian-sterlet/antinoise)

## Set parameters of a model through accuracy estimation

[andy0bsn5.cpp](https://github.com/parthian-sterlet/sitega/blob/master/src/andy0bsn5.cpp)
1. path to fasta files with sets of foreground and background sequences (the last symbol of path must be '/' and '\\' for Linux and Windows OS, respectively)
2. fasta file with set of foreground sequences
3. fasta file with set of background sequences
4. maximal length of one LPD (default value 6)
5. minimal length of motif (integer value respecting to a tested length L, default value is 8)
6. maximal length of motif (default value is 28)
7. step of length of motif (default value is 4, i.e. lengths 8, 12, 16 etc. are considered)
8. cross-validation type specification: positive value below 1 means the ratio of the training subset size to that of control subset for repeated random subsampling validation, default value -1 means equal sizes of training and control subsets, odd/even peaks are used either for training and control subsets)
9. number of iterations in bootatrap (default 2)
10. k-mer length to take into account the sequence bias between foreground and background sequences (default 6, i.e. hexamer frequencies are involved)
11. path to output files (the last symbol of path must be '/' and '\\' for Linux and Windows OS, respectively)
12. maximal peak length (default value is 3000)

## Train a model

[andy05.cpp](https://github.com/parthian-sterlet/sitega/blob/master/src/andy05.cpp)
1. path to fasta files with sets of foreground and background sequences (the last symbol of path must be '/' and '\\' for Linux and Windows OS, respectively)
2. fasta file with set of foreground sequences
3. fasta file with set of background sequences
4. maximal length of one LPD (default value 6)
5. length of motif (integer value respecting to a tested length L, this value is selected by the bootstrap cross-valiation test, see the previous paragraph)
6. size, the number of LPDs (a value is estimated in the bootstrap cross-valiation test, see the previous paragraph)
7. k-mer length to take into account the sequence bias between foreground and background sequences (default 6, i.e. hexamer frequencies are involved)
8. path to output files (the last symbol of path must be '/' and '\\' for Linux and Windows OS, respectively)
9. maximal peak length (default value is 3000)

## Set threshold for a model

[sitega_thr_dist_mat.cpp](https://github.com/parthian-sterlet/sitega/blob/master/src/sitega_thr_dist_mat.cpp)
1. path to file_profile_fasta (see argument #3 below, the last symbol of path must be '/' and '\\' for Linux and Windows OS, respectively)
2. sitega_matrix_file = input file [SiteGA model](https://github.com/parthian-sterlet/sitega/blob/master/examples/model.mat) from **Train a model** module
3. file_profile_fasta = input Whole-genome promoters set in fasta format (unzip files from folder [genomes](https://github.com/parthian-sterlet/sitega/tree/master/genomes), use hs*, mm* and at* files for human, mouse and Arabidopsis data, respectively; [additional file](https://github.com/parthian-sterlet/sitega/blob/master/genomes/prom_all_bed.zip) provides promoters in bed format for other species)
4. output file [Table Threshold_vs_ERR](https://github.com/parthian-sterlet/sitega/blob/master/examples/thr_err), table SiteGA model threshold vs. Expected Recognition Rate (ERR)
5. pvalue_large = maximal ERR (default value 0.001)
6. score_min = lowest threshold of SiteGA model (default value 0.75)
7. dpvalue = granulation value for ERR compaction in [Table Threshold_vs_ERR](https://github.com/parthian-sterlet/sitega/blob/master/examples/thr_err), the default value 0.0000000005 implies the absence of compaction

## Scan test sequences with a model

[andy1_mat.cpp](https://github.com/parthian-sterlet/sitega/blob/master/src/andy1_mat.cpp)
1. path to whole genome sequences of chromosomes in plain format (see above, the last symbol of path must be '/' and '\\' for Linux and Windows OS, respectively)
2. sitega_matrix_file = input file [SiteGA model](https://github.com/parthian-sterlet/sitega/blob/master/examples/model.mat) from **Train a model** module
3. site_description_mode = 0 or 1. 0 means default mode, 1 means computation of frequencies of all LPDs for all tested sequences (option is used for the train fasta file to describe a [SiteGA model](https://github.com/parthian-sterlet/sitega/blob/master/examples/model.mat))
4. ERR threshold = threshold for ERR of [SiteGA model](https://github.com/parthian-sterlet/sitega/blob/master/examples/model.mat) is used to select the SiteGA threshold according to input file [Table Threshold vs ERR](https://github.com/parthian-sterlet/sitega/blob/master/examples/thr_err) from **Set threshold for a model** module
5. input file [Threshold vs ERR table](https://github.com/parthian-sterlet/sitega/blob/master/examples/thr_err) for [SiteGA model](https://github.com/parthian-sterlet/sitega/blob/master/examples/model.mat) threshold selection by ERR threshold
6. output file [Profile with recognized hits](https://github.com/parthian-sterlet/sitega/blob/master/examples/hit_profile)

## Scan whole genome with a model

[andy1_mat_long.cpp](https://github.com/parthian-sterlet/sitega/blob/master/src/andy1_mat_long.cpp)
1. file.seq = input file of test sequences, it has the same format as [Fasta file of peaks](https://github.com/parthian-sterlet/sitega/blob/master/examples/peaks.fa), non ('a', 'c', 'g' and 't') nucleotides are ignored
2. genome release (hg38, mm10, dm6, ce235, sc64 and at10 for H.sapiens, M.musculus, D.,melanogaster, C.elehans, S.cerevisiae and A.thaliana  genomes, respectively)
3. sitega_matrix_file = input file [SiteGA model](https://github.com/parthian-sterlet/sitega/blob/master/examples/model.mat) from **Train a model** module
4. input file [Threshold vs ERR table](https://github.com/parthian-sterlet/sitega/blob/master/examples/thr_err) for [SiteGA model](https://github.com/parthian-sterlet/sitega/blob/master/examples/model.mat) threshold selection by ERR threshold
5. ERR threshold = threshold for ERR of [SiteGA model](https://github.com/parthian-sterlet/sitega/blob/master/examples/model.mat) is used to select the SiteGA threshold according to input file [Table Threshold vs ERR](https://github.com/parthian-sterlet/sitega/blob/master/examples/thr_err) from **Set threshold for a model** module
6. output file [Profile with recognized hits](https://github.com/parthian-sterlet/sitega/blob/master/examples/hit_profile)

# Examples scripts:

These scripts implement various pipelines for Linux:
1. [only training a model](https://github.com/parthian-sterlet/sitega/blob/master/src/train) - **Preparation** and  **Train a model** modules
2. [bootstrap test for a model](https://github.com/parthian-sterlet/sitega/blob/master/src/bootstrap) - **Preparation** and **Estimate accuracy for a model** modules
3. [training and scanning a fasta file with a model](https://github.com/parthian-sterlet/sitega/blob/master/src/scan_fasta) - **Preparation**, **Train a model**, **Set threshold for a model** and **Scan test seauences with a model** modules
4. [training and scanning whole genome with a model](https://github.com/parthian-sterlet/sitega/blob/master/src/scan_genome) - **Preparation**, **Train a model**, **Set threshold for a model** and **Scan whole genome with a model** modules
# Interpretation of results

## Background set generation

[background_genome.cpp](https://github.com/parthian-sterlet/sitega/blob/master/src/background_genome.cpp) prepares the background set that may be used for training (**Train a model** module) or performance evaluation (**Set parameters of a model through accuracy estimation** module)

## Set parameters of a model through accuracy estimation

[andy0bsn5.cpp](https://github.com/parthian-sterlet/sitega/blob/master/src/andy0bsn5.cpp) may several times gradually construct several distinct [SiteGA models](https://github.com/parthian-sterlet/sitega/blob/master/examples/model.mat) (parameter of command line 'number of iterations'), but each time it uses only a part of the foreground set for training, the rest (control) part of set is used to estimate FPR). The output file [Table FPR_vs TPR](https://github.com/parthian-sterlet/sitega/blob/master/examples/model_bs1.txt) represents the table of FPRs for TPR 0.01, 0.02, etc. up to 0.99. Selection of the one model among several ones with different numbers of LPDs and lengths L is performed by respective estimated pAUC values computed for the receiver operating characteristic (ROC) curve (see file with *{train.txt}* extension, [FPR_vs TPR table file](https://github.com/parthian-sterlet/sitega/blob/master/examples/model_bs1.txt). 

## Train a model

[andy05.cpp](https://github.com/parthian-sterlet/sitega/blob/master/src/andy05.cpp) constructs one sitega model, with the numbers of locally positioned dinucleotides (LPDs) assigned according the parameter of the command line, this value deduced from the bootstrap cross validation test (see the previous paragraph). The selected model respecting the maximal pAUC in the bootstrap cross validation test, this model is written in output file [SiteGA model](https://github.com/parthian-sterlet/sitega/blob/master/examples/model.mat) with *.mat* extension

## Set threshold for a model

[sitega_thr_dist_mat.cpp](https://github.com/parthian-sterlet/sitega/blob/master/src/sitega_thr_dist_mat.cpp) computes the distribition of the recognition scores of [SiteGA model](https://github.com/parthian-sterlet/sitega/blob/master/examples/model.mat), output file [Table Threshold vs ERR](https://github.com/parthian-sterlet/sitega/blob/master/examples/thr_err) represents two columns with respective thresholds and ERRs

## Scan test sequences with a model

[andy1_mat.cpp](https://github.com/parthian-sterlet/sitega/blob/master/src/andy1_mat.cpp) takes ready sitega model and threshold  from input file [Table Threshold vs FPR](https://github.com/parthian-sterlet/sitega/blob/master/examples/thr_fpr) and constructs the profile of hits for tested file in fasta format, main output file is [Profile file](https://github.com/parthian-sterlet/sitega/blob/master/examples/hit_profile) I.e. after the header of each peak with first '>' symbol from zero to several lines respect to separate hits, for each hit a start position, score, strand and whole sequence are printed. Interpretation of [SiteGA model](https://github.com/parthian-sterlet/sitega/blob/master/examples/model.mat) is performed through computation of LPD frequencies for the training fasta set (see fourth argument of the command line). Hence, a computed matrix Number_of_sequences vs. Number_of_LPDs can be used for the correlation analysis, e.g. for revealing the most correlated LPDs in a [SiteGA model](https://github.com/parthian-sterlet/sitega/blob/master/examples/model.mat).

## Scan whole genome with a model

[andy1_mat_long.cpp](https://github.com/parthian-sterlet/sitega/blob/master/src/andy1_mat.cpp) takes ready sitega model and threshold  from input file [Table Threshold vs FPR](https://github.com/parthian-sterlet/sitega/blob/master/examples/thr_fpr) and constructs the profile of hits for a genome in plain format, main output file is [Profile file](https://github.com/parthian-sterlet/sitega/blob/master/examples/hit_profile) I.e. after the header of each peak with first '>' symbol from zero to several lines respect to separate hits, for each hit a start position, score, strand and whole sequence are printed. Interpretation of [SiteGA model](https://github.com/parthian-sterlet/sitega/blob/master/examples/model.mat) is performed through computation of LPD frequencies for the training fasta set (see fourth argument of the command line). Hence, a computed matrix Number_of_sequences vs. Number_of_LPDs can be used for the correlation analysis, e.g. for revealing the most correlated LPDs in a [SiteGA model](https://github.com/parthian-sterlet/sitega/blob/master/examples/model.mat).
