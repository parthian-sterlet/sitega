# SiteGA - binding Sites recognition by Genetic Algorithm

# Description
SiteGA implements the algorithm of [Levitsky et al. (2007)](https://doi.org/10.1186/1471-2105-8-481) for *de novo* motif search in a single ChIP-seq dataset (ChIP-seq peaks) [(Tsukanov et al., 2022)](https://doi.org/10.3389/fpls.2022.938545). SiteGA applies a genetic algorithm that evaluates motifs as patterns of frequencies of interdependent locally located dinucleotides (LPDs). This approach models short- and long-range nucleotide context interactions within transcription factor binding sites (TFBS). Thus, the SiteGA approach is radically different from that of the traditional position weight matrix (PWM) model, which estimatessites by sums of additive contributions of nucleotide frequencies from all positions. The PWM approach ignores any dependencies between different positions within the TFBS. 

For a given set of N foreground sequences (peaks), the fitness function F(X) = D(X) * E(X) estimates the alignment X = {x(1), x(2), ..., x(N)} of N best predicted sites in the the SiteGA model. The first factor D(X) reflects dependencies of positions within an alignment through the linear discriminant analysis approach [(Levitsky et al. 2007)](https://doi.org/10.1186/1471-2105-8-481). The second factor E(X) implies the average enrichment of k-mers within the alignment of the predicted sites [(Tsukanov et al., 2022)](https://doi.org/10.3389/fpls.2022.938545). The program indexes all possible positions of sites in the foreground set bwith weights equal to the logarithms of the enrichment of their intrinsic k-mers when comparing foreground and background sequence sets. 

# Requirements
SiteGA source code is written in C++ language. To compile exetubables from the source code you need:

* In Linux system, C++ compiler, e.g. [GCC](https://gcc.gnu.org/) 
* In Windows system any VC++ package, e.g. [Microsoft Visual Studio Community](https://visualstudio.microsoft.com/vs/community/)

# Input data
The bulk of input data are ChIP-seq peaks in FASTA format. To optimize computation time, the default limit of 3000 bp for the length of any peak is used, although the length of a peak is not restricted by the algorithm. This same reason requires at least a moderate number of peaks to be used in the analysis, usually ~1000 peaks are sufficient to obtain a motif.

# Source code
Folder [**src**](https://github.com/parthian-sterlet/sitega/tree/master/src) contains SiteGA source code files that relate to individual pipeline modules and are described below.

## 1. Background sequence set generation
The background set of sequences is required as a complement to the foreground set to select two parameters of the SiteGA model: motif length and the number of LPDs. The main purpose of the background set is to exclude artifact motifs related to a genome-specific sequence content bias, e.g. polyA, from the results of *de novo* motif search. To prepare genome-specific sets of background sequences, it is recommended to use the [AntiNoise](https://github.com/parthian-sterlet/antinoise) package. 
## 2. Set parameters of a model through accuracy estimation
[andy0bsn5.cpp](https://github.com/parthian-sterlet/sitega/blob/master/src/andy0bsn5.cpp) performs the bootsrap cross-validation test to select parameters of a model providing the best performance for given foreground and background sets. A model has parameters of motif length and the number of LPDs. The program sequentially checks all possible combinations of motif length and number of LPDs. The minimal, maximal and step of the motif length, and those values for the number of LPDs are the parameters of command line, default lengths are 8, 12, 16 and 20 nt, default numbers of LPDs are 40, 60, 80 and 100. The bootsrap cross-validation test denotes the partitioning of the foreground set into subsets of training and control sequence sets, the former is used to train a model and the latter to measure its performance. The maximal partial area under curve (pAUC) and area under precision-recall curve (AUCPR) are used to estimate the accuracy of a model for certain motif length and the number of LPDs. The term partial means that only the part of a ROC curve respecting the criterion FPR < 0.001 is impied for pAUC computation. The receiver operating characteristic (ROC) curve means the dependence between True Positive Rate (TPR) and False Positive Rate (FPR). TPR/FPR (axes Y/X of the ROC curve) are defined as the fraction of sequences from the foreground set containing predicted sites and the frequency of predicted sites in the background set, respectively. FPR = Nb/Wb, here Nb - number of predicted sites in the background set, Wb - total number of checked posiitons for predicted sites in the background set. The dependence between Precision and Recall represents AUCPR. Recall = TPR = TP/(TP+FN), Precision = TP/(TP+FP), here TP & FP denote counts of predicted sequences from foreground & background sets, FN denote counts of not predicted sequences from the foreground set. 
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

git clone https://github.com/parthian-sterlet/sitega

cd sitega/src

chmod a+x build.sh

./build.sh

* In Windiws system:

separate compilation of all source files in VC++

## Integration of modules

Scheme of modules functioning is given below

![scheme](https://github.com/parthian-sterlet/sitega/blob/master/examples/scheme_github_sitega9.png)

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
2. fasta file, set of foreground sequences
3. fasta file, set of background sequences
4. integer value, maximal length of one LPD (default value 6)
5. integer value, minimal length of motif (integer value respecting to a tested length L, default value is 8)
6. integer value, maximal length of motif (default value is 28)
7. integer value, step of length of motif (default value is 4, i.e. lengths 8, 12, 16 etc. are considered)
8. cross-validation type specification: positive value below 1 means the ratio of the training subset size to that of control subset for repeated random subsampling validation, default value -1 means equal sizes of training and control subsets, odd/even peaks are used either for training and control subsets)
9. integer value, number of iterations in bootatrap (default 2)
10. integer value, k-mer length to take into account the sequence bias between foreground and background sequences (default 6, i.e. hexamer frequencies are involved)
11. path to output files (the last symbol of path must be '/' and '\\' for Linux and Windows OS, respectively)
12. integer value, maximal peak length (default value is 3000)
13. output log file
    
## Train a model

[andy05.cpp](https://github.com/parthian-sterlet/sitega/blob/master/src/andy05.cpp)
1. path to fasta files with sets of foreground and background sequences (the last symbol of path must be '/' and '\\' for Linux and Windows OS, respectively)
2. fasta file with set of foreground sequences
3. fasta file with set of background sequences
4. integer value, maximal length of one LPD (default value 6)
5. integer value, length of motif (integer value respecting to a tested length L, this value is selected by the bootstrap cross-valiation test, see the previous paragraph)
6. integer value, size, the number of LPDs (a value is estimated in the bootstrap cross-valiation test, see the previous paragraph)
7. integer value, k-mer length to take into account the sequence bias between foreground and background sequences (default 6, i.e. hexamer frequencies are involved)
8. path to output files (the last symbol of path must be '/' and '\\' for Linux and Windows OS, respectively)
9. integer value, maximal peak length (default value is 3000)
10. output log file

## Set threshold for a model

[sitega_thr_dist_mat.cpp](https://github.com/parthian-sterlet/sitega/blob/master/src/sitega_thr_dist_mat.cpp)
1. path to file_profile_fasta (see argument #3 below, the last symbol of path must be '/' and '\\' for Linux and Windows OS, respectively)
2. sitega_matrix_file = input file [SiteGA model](https://github.com/parthian-sterlet/sitega/blob/master/examples/model.mat) from **Train a model** module
3. file_profile_fasta = input Whole-genome promoters set in fasta format ([example](https://github.com/parthian-sterlet/sitega/blob/master/genomes/sc64/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa) provides one for _S. cerevisiae_)
4. output file [Table Threshold_vs_ERR](https://github.com/parthian-sterlet/sitega/blob/master/examples/thr_err), table SiteGA model threshold vs. Expected Recognition Rate (ERR)
5. double pvalue_large = maximal ERR (default value 0.001)
6. double score_min = lowest threshold of SiteGA model (default value 0.75)
7. double dpvalue = granulation value for ERR compaction in [Table Threshold_vs_ERR](https://github.com/parthian-sterlet/sitega/blob/master/examples/thr_err), the default value 0.0000000005 implies the absence of compaction

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
1. path to chromosome sequences in [Fasta format](https://github.com/parthian-sterlet/sitega/blob/master/examples/peaks.fa), these are files for individual chromosomes in separate files, they are prepared by splitting  [the whole-genome fasta file](https://github.com/parthian-sterlet/sitega/blob/master/genomes/sc64/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa) by the [https://github.com/parthian-sterlet/antinoise/blob/main/src/fasta_muliplefiles.cpp](https://github.com/parthian-sterlet/antinoise/blob/main/src/fasta_muliplefiles.cpp) from the [AntiNoise](https://github.com/parthian-sterlet/antinoise) package.
2. species and genome release (values hg38, mm10, rn6, zf11, dm6, and ce235; at10, gm21, zm73, and mp61; sc64 and sch294). The animals inludes human _Homo sapiens_ hg38, mouse _Mus musculus_ mm10, rat _Rattus norvegicus_ Rnor_6.0, zebrafish _Danio rerio_ GRCz11, fly _Drosophila melanogaster_ dm6, and roundworm _Caenorhabditis elegans_ WBcel235; the plants are arabidopsis _Arabidopsis thaliana_ TAIR10, soybean _Glycine max_ v2.1, maize _Zea mays_ B73, and liverwort _Marchantia polymorpha_ MpTak v6.1; the fungi are baker's yeast _Saccharomyces cerevisiae_ R64-1-1 and fission yeast _Schizosaccharomyces pombe_ ASM294v2.
3. sitega_matrix_file = input file [SiteGA model](https://github.com/parthian-sterlet/sitega/blob/master/examples/model.mat) from **Train a model** module
4. input file [Threshold vs ERR table](https://github.com/parthian-sterlet/sitega/blob/master/examples/thr_err) for [SiteGA model](https://github.com/parthian-sterlet/sitega/blob/master/examples/model.mat) threshold selection by ERR threshold
5. double value, ERR threshold = threshold for ERR of [SiteGA model](https://github.com/parthian-sterlet/sitega/blob/master/examples/model.mat) is used to select the SiteGA threshold according to input file [Table Threshold vs ERR](https://github.com/parthian-sterlet/sitega/blob/master/examples/thr_err) from **Set threshold for a model** module
6. output file [Profile with recognized hits](https://github.com/parthian-sterlet/sitega/blob/master/examples/hit_profile)

# Examples scripts:

These scripts implement various pipelines for Linux:
1. [only training a model](https://github.com/parthian-sterlet/sitega/blob/master/src/train) - **Preparation** and  **Train a model** modules
2. [bootstrap test for a model](https://github.com/parthian-sterlet/sitega/blob/master/src/bootstrap) - **Preparation** and **Estimate accuracy for a model** modules
3. [training and scanning a fasta file with a model](https://github.com/parthian-sterlet/sitega/blob/master/src/scan_fasta) - **Preparation**, **Train a model**, **Set threshold for a model** and **Scan test seauences with a model** modules
4. [training and scanning whole genome with a model](https://github.com/parthian-sterlet/sitega/blob/master/src/scan_genome) - **Preparation**, **Train a model**, **Set threshold for a model** and **Scan whole genome with a model** modules
   
# Interpretation of results

## Background set generation

[background_genome.cpp](https://github.com/parthian-sterlet/antinoise/blob/master/src/background_genome.cpp) from the [AntiNoise](https://github.com/parthian-sterlet/antinoise) package prepares the background set that may be used for training (**Train a model** module) or performance evaluation (**Set parameters of a model through accuracy estimation** module)

## Set parameters of a model through accuracy estimation

[andy0bsn5.cpp](https://github.com/parthian-sterlet/sitega/blob/master/src/andy0bsn5.cpp) several times construct distinct [SiteGA models](https://github.com/parthian-sterlet/sitega/blob/master/examples/model.mat) (parameter of command line 'number of iterations'), but each time it uses only a part of the foreground set for training, the rest (control) part of set is used to estimate FPR). The output file [Table pAUC file](https://github.com/parthian-sterlet/sitega/blob/master/examples/PEAKS035427_ATOH1_P48985_MACS2.fa_roc_bs1.txt) reoresenting pAUC values computed by [ROC curves](https://github.com/parthian-sterlet/sitega/blob/master/examples/PEAKS035427_ATOH1_P48985_MACS2.fa_roc_bs1.txt) for various values of parameters of a model allows selection of the one model among several ones. These several ones come from different numbers of LPDs and lengths L of a motif. 

## Train a model

[andy05.cpp](https://github.com/parthian-sterlet/sitega/blob/master/src/andy05.cpp) constructs one sitega model, with the numbers of locally positioned dinucleotides (LPDs) assigned according the parameter of the command line, this value deduced from the bootstrap cross validation test (see the previous paragraph). The selected model respecting the maximal pAUC in the bootstrap cross validation test, this model is written in output file [SiteGA model](https://github.com/parthian-sterlet/sitega/blob/master/examples/model.mat) with *.mat* extension

## Set threshold for a model

[sitega_thr_dist_mat.cpp](https://github.com/parthian-sterlet/sitega/blob/master/src/sitega_thr_dist_mat.cpp) computes the distribition of the recognition scores of [SiteGA model](https://github.com/parthian-sterlet/sitega/blob/master/examples/model.mat), output file [Table Threshold vs ERR](https://github.com/parthian-sterlet/sitega/blob/master/examples/thr_err) represents two columns with respective thresholds and ERRs

## Scan test sequences with a model

[andy1_mat.cpp](https://github.com/parthian-sterlet/sitega/blob/master/src/andy1_mat.cpp) takes ready sitega model and threshold  from input file [Table Threshold vs FPR](https://github.com/parthian-sterlet/sitega/blob/master/examples/thr_fpr) and constructs the profile of hits for tested file in fasta format, main output file is [Profile file](https://github.com/parthian-sterlet/sitega/blob/master/examples/hit_profile) I.e. after the header of each peak with first '>' symbol from zero to several lines respect to separate hits, for each hit a start position, score, strand and whole sequence are printed. Interpretation of [SiteGA model](https://github.com/parthian-sterlet/sitega/blob/master/examples/model.mat) is performed through computation of LPD frequencies for the training fasta set (see fourth argument of the command line). Hence, a computed matrix Number_of_sequences vs. Number_of_LPDs can be used for the correlation analysis, e.g. for revealing the most correlated LPDs in a [SiteGA model](https://github.com/parthian-sterlet/sitega/blob/master/examples/model.mat).

## Scan whole genome with a model

[andy1_mat_long.cpp](https://github.com/parthian-sterlet/sitega/blob/master/src/andy1_mat.cpp) takes ready sitega model and threshold  from input file [Table Threshold vs FPR](https://github.com/parthian-sterlet/sitega/blob/master/examples/thr_fpr) and constructs the profile of hits for a genome in plain format, main output file is [Profile file](https://github.com/parthian-sterlet/sitega/blob/master/examples/hit_profile) I.e. after the header of each peak with first '>' symbol from zero to several lines respect to separate hits, for each hit a start position, score, strand and whole sequence are printed. Interpretation of [SiteGA model](https://github.com/parthian-sterlet/sitega/blob/master/examples/model.mat) is performed through computation of LPD frequencies for the training fasta set (see fourth argument of the command line). Hence, a computed matrix Number_of_sequences vs. Number_of_LPDs can be used for the correlation analysis, e.g. for revealing the most correlated LPDs in a [SiteGA model](https://github.com/parthian-sterlet/sitega/blob/master/examples/model.mat).
