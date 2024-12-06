# SiteGA - binding Sites recognition by Genetic Algorithm

# Description
SiteGA implements the algorithm of [Levitsky et al. (2007)](https://doi.org/10.1186/1471-2105-8-481) for *de novo* motif search in a single ChIP-seq dataset (ChIP-seq peaks) [(Tsukanov et al., 2022)](https://doi.org/10.3389/fpls.2022.938545). SiteGA applies a genetic algorithm that evaluates motifs as patterns of frequencies of interdependent locally located dinucleotides (LPDs). This approach models short- and long-range interactions of nucleotide context within transcription factor binding sites (TFBS). Thus, the SiteGA approach is quite different from that of the traditional position weight matrix (PWM) model, which estimates sites by sums of additive contributions of nucleotide frequencies from all positions. The PWM approach ignores any dependencies between different positions within the TFBS. 

For a given set of N foreground sequences (peaks), the fitness function F(X) = D(X) * E(X) estimates the alignment X = {x(1), x(2), ..., x(N)} of N best predicted sites in the the SiteGA model. The first factor D(X) reflects dependencies of positions within an alignment through the linear discriminant analysis approach [(Levitsky et al. 2007)](https://doi.org/10.1186/1471-2105-8-481). The second factor E(X) implies the average enrichment of k-mers within the alignment of the predicted sites [(Tsukanov et al., 2022)](https://doi.org/10.3389/fpls.2022.938545). The program indexes all possible positions of sites in the foreground set bwith weights equal to the logarithms of the enrichment of their intrinsic k-mers when comparing foreground and background sequence sets. 

# Requirements
SiteGA source code is written in C++ language. To compile exetubables from the source code you need:

* In Linux system, C++ compiler, the minimal one is [GCC](https://gcc.gnu.org/) compiler, more prominent is [Clang](https://clang.llvm.org/) under the management of [Spack](https://spack.readthedocs.io/en/latest/getting_started.html) tool
* In Windows system any VC++ package, e.g. [Microsoft Visual Studio Community](https://visualstudio.microsoft.com/vs/community/)

# Input data
The bulk of input data are ChIP-seq peaks in FASTA format. To optimize computation time, the default limit of 3000 bp for the length of any peak is used, although the length of a peak is not restricted by the algorithm. This same reason requires at least a moderate number of peaks to be used in the analysis, usually ~1000 peaks are sufficient to obtain a motif.

# Source code
Folder [**src**](https://github.com/parthian-sterlet/sitega/tree/master/src) contains SiteGA source code files that relate to individual pipeline modules and are described below.

# Binary files
Folder [**bin**](https://github.com/parthian-sterlet/sitega/tree/master/bin) contains SiteGA source code files compiled by [Clang](https://clang.llvm.org/) under the management of [Spack](https://spack.readthedocs.io/en/latest/getting_started.html) tool. It is recomended to use these binary files instead of those compiled by default g++ compiler due to computation speed. Just copy these binary files to the [**src**](https://github.com/parthian-sterlet/sitega/tree/master/src) folder and ensure that they are functional. For illustrative purpose only all scripts are started from binary files compiled with default g++ compiler, see [build.sh](https://github.com/parthian-sterlet/sitega/blob/master/src/build.sh). 

## 1. Background sequence set generation
The background set of sequences is required as a complement to the foreground set to select the values of two parameters of the SiteGA model: motif length and the number of LPDs. The main purpose of the background set is to exclude artifact motifs related to a genome-specific sequence content bias, e.g. polyA, from the results of *de novo* motif search. To prepare genome-specific sets of background sequences, it is recommended to use the [AntiNoise](https://github.com/parthian-sterlet/antinoise) package [(Raditsa et al., 2024)](https://doi.org/10.1093/nargab/lqae090). 
## 2. Set parameters of a model through accuracy estimation
[andy0bsn5cell.cpp](https://github.com/parthian-sterlet/sitega/blob/master/src/andy0bsn5cell.cpp) performs the bootsrap cross-validation test to select parameters of a model providing the best performance for given foreground and background sets. A model has parameters of the motif length and the range for the number of LPDs. The motif length, the minimal value for the number of LPDs and it step are the parameters of command line. Recommended motif lengths are from 8 to 12 nt. The genetic algorithm simultaneously tests multiple values of the number of LPDs. The internal parameters of genetic algorithm are the range of the number of LPDs and its increment, default range consists of sixteen numbers of LPDs, 18, 20, etc. up to 48. The bootsrap cross-validation test denotes the partitioning of the foreground set into subsets of training and control sequence sets, the former is used to train a model and the latter to measure its performance. The maximal partial area under curve (pAUC) and area under precision-recall curve (AUCPR) are used to estimate the accuracy of a model for certain motif length and the number of LPDs. The term partial means that only the parts of the ROC or PR curve is impied for pAUC computation. In both cureves the criterion of motif frequency (Expected Recognition Rate, ERR) below 0.01 is applied. 
The receiver operating characteristic (ROC) curve means the dependence between True Positive Rate (TPR) and False Positive Rate (FPR). TPR/FPR (axes Y/X of the ROC curve) are defined as the fraction of sequences from the foreground set containing predicted sites and the frequency of predicted sites in the background set, respectively. For ROC curve FPR = Nb/Wb, here Nb - number of predicted sites in the background set, Wb - total number of checked positions for predicted sites in the background set. 
The PR curve is the dependence between Precision and Recall, the matrics of PR curve is the partial area under PR curve, pAUPRC. Recall = TPR = TP/(TP+FN), Precision = TPR/(TPR+FPR) = [TP/NF]/[TP/NF+FP/NB] = TP/[TP+FP*(NF/NF)], here TPR & FPR are fraction of all NF foreground & NB background sequences that are predicted, TP & FP denote counts of predicted sequences from foreground & background sets. The calculation of AUPRC also implies the precision values corrected for their expectation value 0.5, see [MetArea](https://github.com/parthian-sterlet/metarea). The two-fold cross-validation procedure is performed, its two iterations mean odd/even peaks either for training and control subsets. The sizes of training and control subsets are approximately equal.
## 3. Train a model
[andy05cell.cpp](https://github.com/parthian-sterlet/sitega/blob/master/src/andy05cell.cpp) trains a SiteGA model with selected by the bootstrap cross-validation procedure parameters of the motif length and the number of LPDs for given foreground and background sets. The resulting model is written in a [special matrix file](https://github.com/parthian-sterlet/sitega/blob/master/examples/PEAKS035427_ATOH1_P48985_MACS2.mat)
## 4. Set threshold for a model
[sitega_thr_dist_mat.cpp](https://github.com/parthian-sterlet/sitega/blob/master/src/sitega_thr_dist_mat.cpp) creates table of thresholds for model's hits search in test sequences (**Scan test sequences with a model** module) based on the distribution of recognition score for the set of whole-genome promoter sequences selected for the respective species. Prepeared sets are stored in [genomes](https://github.com/parthian-sterlet/sitega/tree/master/genomes) folder. The threshold selection implies the estimation of Expected Recognition Rate (ERR) of a model for promoter sequences of whole genome. The dependence of the threshold from ERR is stored in a [special file](https://github.com/parthian-sterlet/sitega/blob/master/examples/PEAKS035427_ATOH1_P48985_MACS2_sga.dist) containing a list of two columns: recognition thresholds in descending order and corresponding -Log10(ERR) values
## 5. Scan test sequences with a model
[andy1_mat.cpp](https://github.com/parthian-sterlet/sitega/blob/master/src/1_mat.cpp) scans test sequences with a constructed model and the selected threshold of a model.
## 6. Scan whole genome with a model
[andy1_mat_long.cpp](https://github.com/parthian-sterlet/sitega/blob/master/src/1_mat_long.cpp) scans whole genome as a set of full-sized chromosomes with a constructed model and the selected threshold of a model.

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

![scheme](https://github.com/parthian-sterlet/sitega/blob/master/examples/github_sitega.jpg)

Modules **Set parameters of a model through accuracy estimation** and **Train a model** take the input data of the foreground and background sequence sets, the background set is prepared by the **Background set generation** module

The module **Set parameters of a model through accuracy estimation** is required for functionality of the module **Train a model** and all subsequent modules since the bootstrap procedure correctly selects parameters of a model (see output data block **ROC & PR curves, pAUC ROC & pAUPRC performance metrics**)

Modules **Set threshold for a model** and **Scan test seauences with a model** require file with [SiteGA model](https://github.com/parthian-sterlet/sitega/blob/master/examples/PEAKS035427_ATOH1_P48985_MACS2.mat) which should be previously prepared by the module **Train a model** 

The module **Set threshold for a model** is required to select a correct threshold for **Scan test sequences with a model** module

# How to run separate modules
Lists of command line arguments for all modules are described below

## Background set generation

Use [background_genome_mono.cpp](https://github.com/parthian-sterlet/antinoise/blob/main/src/background_genome_mono.cpp) from the github repositiory [AntiNoise](https://github.com/parthian-sterlet/antinoise)

## Set parameters of a model through accuracy estimation

[andy0bsn5cell.cpp](https://github.com/parthian-sterlet/sitega/blob/master/src/andy0bsn5cell.cpp)
1. path to fasta files with sets of foreground and background sequences (the last symbol of path must be '/' and '\\' for Linux and Windows OS, respectively)
2. fasta file, set of foreground sequences
3. fasta file, set of background sequences
4. integer value, maximal length of one LPD (default value 6)
5. integer value, minimal number of LPDs (LPDmin, default value in the range from 20 to 50)
6. integer value, step of the number of LPDs (default value 1, four numbers are {LPDmin, LPDmin + 1, LPDmin + 2, LPDmin + 3} since the default number of distinct numbers of LPDs is four)
7. integer value, length of motif in nucleotides, the motif length from 8 to 12 nt is recommended. The computation time increases with the growth of motif length. 
8. double value, Expected Recognition Rate - maximum frequency of motif in the background set. The value is used to restrict the X axes in ROC and PR curves to define partial areas under curves.
9. integer value, k-mer length to take into account the sequence bias between foreground and background sequences (default 6, i.e. hexamer frequencies are involved)
10. path to output files (the last symbol of path must be '/' and '\\' for Linux and Windows OS, respectively)
11. integer value, maximal peak length (default value is 3000)
12. output log file
    
## Train a model

[andy05cell.cpp](https://github.com/parthian-sterlet/sitega/blob/master/src/andy05cell.cpp)
1. path to fasta files with sets of foreground and background sequences (the last symbol of path must be '/' and '\\' for Linux and Windows OS, respectively)
2. fasta file with set of foreground sequences
3. fasta file with set of background sequences
4. integer value, maximal length of one LPD (default value 6)
5. integer value, length of motif (this value is selected in the bootstrap cross-valiation test, [see above](https://github.com/parthian-sterlet/sitega/tree/master#set-parameters-of-a-model-through-accuracy-estimation))
6. integer value, size, the number of LPDs (a value is selected in the bootstrap cross-valiation test, [see above](https://github.com/parthian-sterlet/sitega/tree/master#set-parameters-of-a-model-through-accuracy-estimation))
7. integer value, k-mer length to take into account the sequence bias between foreground and background sequences (default 6, i.e. hexamer frequencies are involved)
8. path to output files (the last symbol of path must be '/' and '\\' for Linux and Windows OS, respectively)
9. integer value, maximal peak length (default value is 3000)
10. output log file

## Set threshold for a model

[sitega_thr_dist_mat.cpp](https://github.com/parthian-sterlet/sitega/blob/master/src/sitega_thr_dist_mat.cpp)
1. path to file_profile_fasta (see argument #3 below, the last symbol of path must be '/' and '\\' for Linux and Windows OS, respectively)
2. sitega_matrix_file = input file [SiteGA model](https://github.com/parthian-sterlet/sitega/blob/master/examples/PEAKS035427_ATOH1_P48985_MACS2.mat) from **Train a model** module
3. file_profile_fasta = input Whole-genome promoters set in fasta format ([archive with example](https://github.com/parthian-sterlet/sitega/blob/master/genomes/mm10/ups2kb_mm10.seq.tar.gz) provides one for _M. musculus_)
4. output file, text mode, [Table Threshold_vs_ERR](https://github.com/parthian-sterlet/sitega/blob/master/examples/PEAKS035427_ATOH1_P48985_MACS2_sga.dist), table SiteGA model threshold vs. Expected Recognition Rate (ERR)
5. output file, binary mode, [SiteGA model](https://github.com/parthian-sterlet/sitega/blob/master/examples/PEAKS035427_ATOH1_P48985_MACS2.mat) as a matrix and [Table Threshold_vs_ERR](https://github.com/parthian-sterlet/sitega/blob/master/examples/PEAKS035427_ATOH1_P48985_MACS2_sga.dist), table SiteGA model threshold vs. ERR)
6. double pvalue_large = maximal ERR, default value 0.001, value below 0.01 is recomended
7. double dpvalue = granulation value for ERR compaction in [Table Threshold_vs_ERR](https://github.com/parthian-sterlet/sitega/blob/master/examples/PEAKS035427_ATOH1_P48985_MACS2.mat), the default value 0.0000005 implies a moderate compaction

## Scan test sequences with a model

[andy1_mat.cpp](https://github.com/parthian-sterlet/sitega/blob/master/src/andy1_mat.cpp)
1. test file in fasta format with its path
2. sitega_matrix_file = input file [SiteGA model](https://github.com/parthian-sterlet/sitega/blob/master/examples/PEAKS035427_ATOH1_P48985_MACS2.mat) from **Train a model** module
3. site_description_mode = 0 or 1. the value 0 means default mode, the value 1 means computation of frequencies of all LPDs for all tested sequences (option is used for the train fasta file to describe a [SiteGA model](https://github.com/parthian-sterlet/sitega/blob/master/examples/PEAKS035427_ATOH1_P48985_MACS2.mat))
4. ERR threshold = threshold for ERR of [SiteGA model](https://github.com/parthian-sterlet/sitega/blob/master/examples/PEAKS035427_ATOH1_P48985_MACS2.mat) is used to select the SiteGA threshold according to input file [Threshold vs ERR table](https://github.com/parthian-sterlet/sitega/blob/master/examples/PEAKS035427_ATOH1_P48985_MACS2_sga.dist) from **Set threshold for a model** module
5. input file [Threshold vs ERR table](https://github.com/parthian-sterlet/sitega/blob/master/examples/PEAKS035427_ATOH1_P48985_MACS2_sga.dist) for [SiteGA model](https://github.com/parthian-sterlet/sitega/blob/master/examples/PEAKS035427_ATOH1_P48985_MACS2.mat) threshold selection by ERR threshold
6. output file [Profile with recognized hits](https://github.com/parthian-sterlet/sitega/blob/master/examples/sga_out)

## Scan whole genome with a model

[andy1_mat_long.cpp](https://github.com/parthian-sterlet/sitega/blob/master/src/andy1_mat_long.cpp)
1. path to chromosome sequences in [Plain format](https://github.com/parthian-sterlet/sitega/blob/master/examples/chrI.plain), individual chromosomes are in separate files, they are prepared by splitting  [the whole-genome fasta file](https://github.com/parthian-sterlet/sitega/blob/master/genomes/sc64/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa) by the [https://github.com/parthian-sterlet/antinoise/blob/main/src/fasta_muliplefiles.cpp](https://github.com/parthian-sterlet/antinoise/blob/main/src/fasta_muliplefiles.cpp) and [https://github.com/parthian-sterlet/antinoise/blob/main/src/fasta_to_plain0.cpp](https://github.com/parthian-sterlet/antinoise/blob/main/src/fasta_to_plain0.cpp) from the [AntiNoise](https://github.com/parthian-sterlet/antinoise) package.
2. species and genome release (values hg38, mm10, rn6, zf11, dm6, and ce235; at10, gm21, zm73, and mp61; sc64 and sch294). The animals inludes human _Homo sapiens_ hg38, mouse _Mus musculus_ mm10, rat _Rattus norvegicus_ Rnor_6.0, zebrafish _Danio rerio_ GRCz11, fly _Drosophila melanogaster_ dm6, and roundworm _Caenorhabditis elegans_ WBcel235; the plants are arabidopsis _Arabidopsis thaliana_ TAIR10, soybean _Glycine max_ v2.1, maize _Zea mays_ B73, and liverwort _Marchantia polymorpha_ MpTak v6.1; the fungi are baker's yeast _Saccharomyces cerevisiae_ R64-1-1 and fission yeast _Schizosaccharomyces pombe_ ASM294v2.
3. sitega_matrix_file = input file [SiteGA model](https://github.com/parthian-sterlet/sitega/blob/master/examples/PEAKS035427_ATOH1_P48985_MACS2.mat) from **Train a model** module
4. input file [Threshold vs ERR table](https://github.com/parthian-sterlet/sitega/blob/master/examples/PEAKS035427_ATOH1_P48985_MACS2_sga.dist) for [SiteGA model](https://github.com/parthian-sterlet/sitega/blob/master/examples/PEAKS035427_ATOH1_P48985_MACS2.mat) threshold selection by ERR threshold
5. double value, ERR threshold = threshold for ERR of [SiteGA model](https://github.com/parthian-sterlet/sitega/blob/master/examples/PEAKS035427_ATOH1_P48985_MACS2.mat) is used to select the SiteGA threshold according to input file [Table Threshold vs ERR](https://github.com/parthian-sterlet/sitega/blob/master/examples/PEAKS035427_ATOH1_P48985_MACS2_sga.dist) from **Set threshold for a model** module
6. output file [Profile with recognized hits](https://github.com/parthian-sterlet/sitega/blob/master/examples/sga_long_out)

# Examples scripts:

These scripts implement various pipelines for Linux:
1. [bootstrap test for a model](https://github.com/parthian-sterlet/sitega/blob/master/src/bootstrap) - **Estimate accuracy for a model** module
2. [training a model](https://github.com/parthian-sterlet/sitega/blob/master/src/train) - **Train a model** module
3. [thresholds selection by ERRs](https://github.com/parthian-sterlet/sitega/blob/master/src/thr_err) - **Set threshold for a model** module
4. [sites recognition for a test fasta file with a model](https://github.com/parthian-sterlet/sitega/blob/master/src/scan_fasta) - **Scan test seauences with a model** module
5. [sites recognition for a whole genome with a model](https://github.com/parthian-sterlet/sitega/blob/master/src/scan_genome) - **Scan whole genome with a model** module
   
# Interpretation of results

## Background set generation

[background_genome.cpp](https://github.com/parthian-sterlet/antinoise/blob/master/src/background_genome.cpp) from the [AntiNoise](https://github.com/parthian-sterlet/antinoise) package prepares the background set that may be used for training (**Train a model** module) or performance evaluation (**Set parameters of a model through accuracy estimation** module)

## Set parameters of a model through accuracy estimation

[andy0bsn5cell.cpp](https://github.com/parthian-sterlet/sitega/blob/master/src/andy0bsn5cell.cpp) several times construct distinct [SiteGA models](https://github.com/parthian-sterlet/sitega/blob/master/examples/PEAKS035427_ATOH1_P48985_MACS2.mat) (parameter of command line 'number of iterations'), but each time it uses only a part of the foreground set for training, the rest (control) part of set is used to estimate FPR). The output files [ROC curve (pAUC ROC)](https://github.com/parthian-sterlet/sitega/blob/master/examples/PEAKS035427_ATOH1_P48985_MACS2_roc_bs1.txt) and [PR curve (pAUPRC)](https://github.com/parthian-sterlet/sitega/blob/master/examples/PEAKS035427_ATOH1_P48985_MACS2_prc_bs1.txt) for various values of parameters of a model allows selection of the one model among several ones. These several ones come from different numbers of LPDs and lengths L of a motif. [Summary file](https://github.com/parthian-sterlet/sitega/blob/master/examples/PEAKS035427_ATOH1_P48985_MACS2_auc_bs.txt) lists pAUC ROC and pAUPRC values for various numbers of LPDs.

## Train a model

[andy05cell.cpp](https://github.com/parthian-sterlet/sitega/blob/master/src/andy05cell.cpp) constructs one sitega model, with the numbers of locally positioned dinucleotides (LPDs) assigned according the parameter of the command line, this value deduced from the bootstrap cross validation test (see the previous paragraph). The selected model respecting the maximal pAUC in the bootstrap cross validation test, this model is written in output file [SiteGA model](https://github.com/parthian-sterlet/sitega/blob/master/examples/PEAKS035427_ATOH1_P48985_MACS2.mat) with *.mat* extension

## Set threshold for a model

[sitega_thr_dist_mat.cpp](https://github.com/parthian-sterlet/sitega/blob/master/src/sitega_thr_dist_mat.cpp) computes the distribition of the recognition scores of [SiteGA model](https://github.com/parthian-sterlet/sitega/blob/master/examples/PEAKS035427_ATOH1_P48985_MACS2.mat), output file [Table Threshold vs ERR](https://github.com/parthian-sterlet/sitega/blob/master/examples/PEAKS035427_ATOH1_P48985_MACS2_sga.dist) represents two columns with respective thresholds and ERRs

## Scan test sequences with a model

[andy1_mat.cpp](https://github.com/parthian-sterlet/sitega/blob/master/src/andy1_mat.cpp) takes ready sitega model and threshold  from input file [Table Threshold vs FPR](https://github.com/parthian-sterlet/sitega/blob/master/examples/thr_fpr) and constructs the profile of hits for tested file in fasta format, main output file is [Profile file](https://github.com/parthian-sterlet/sitega/blob/master/examples/sga_out) I.e. after the header of each peak with first '>' symbol from zero to several lines respect to separate hits, for each hit a start position, score, strand and whole sequence are printed. Interpretation of [SiteGA model](https://github.com/parthian-sterlet/sitega/blob/master/examples/PEAKS035427_ATOH1_P48985_MACS2.mat) is performed through computation of LPD frequencies for the training fasta set (see fourth argument of the command line). Hence, a computed matrix Number_of_sequences vs. Number_of_LPDs can be used for the correlation analysis, e.g. for revealing the most correlated LPDs in a [SiteGA model](https://github.com/parthian-sterlet/sitega/blob/master/examples/PEAKS035427_ATOH1_P48985_MACS2.mat).

## Scan whole genome with a model

[andy1_mat_long.cpp](https://github.com/parthian-sterlet/sitega/blob/master/src/andy1_mat.cpp) takes ready sitega model and threshold  from input file [Table Threshold vs FPR](https://github.com/parthian-sterlet/sitega/blob/master/examples/thr_fpr) and constructs the profile of hits for a genome in plain format, main output file is [Profile file](https://github.com/parthian-sterlet/sitega/blob/master/examples/sga_long_out) I.e. after the header of each peak with first '>' symbol from zero to several lines respect to separate hits, for each hit a start position, score, strand and whole sequence are printed. Interpretation of [SiteGA model](https://github.com/parthian-sterlet/sitega/blob/master/examples/PEAKS035427_ATOH1_P48985_MACS2.mat) is performed through computation of LPD frequencies for the training fasta set (see fourth argument of the command line). Hence, a computed matrix Number_of_sequences vs. Number_of_LPDs can be used for the correlation analysis, e.g. for revealing the most correlated LPDs in a [SiteGA model](https://github.com/parthian-sterlet/sitega/blob/master/examples/PEAKS035427_ATOH1_P48985_MACS2.mat).
