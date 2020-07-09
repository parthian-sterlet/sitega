# sitega
folder src contains four files with sitega source codes, c++ language

monte0dg.cpp prepare parameter file to train a model (andy02.cpp) or perform the bootsrap cross validation test (andy0bsn2.cpp)
command line arguments:
1int reg -2file seq  -3file out

1int reg = length of region (default value 6)
2file seq = peaks (fasta file) 
3file out = output file = parameter file from monte0dg.cpp

andy02.cpp trains a model with a given train ChIP-seq dataset (peaks)
command line arguments:
1char file_cor 2int motif_len 3int size_start 4int size_end 5int size_dif

1char file_cor = parameter file from monte0dg.cpp
2int motif_len = length of motif (integer value above 30 is recommended, default value 30)
3int size_start = start value for the number of locally positioned dinucleotides (default value 10)
4int size_end = end value for the number of locally positioned dinucleotides (default value 90)
5int size_dif = variation value for the number of locally positioned dinucleotides (default value 10)

andy0bsn2.cpp performs the bootsrap cross validation test to estimate the performance of a model with a given train ChIP-seq dataset
command line arguments:
1char file_cor 2int motif_len 3int size_start 4int size_end 5int size_dif 6double ratio_cnt_of_all(0=jk) 7int num_iterations 

1char file_cor = parameter file from monte0dg.cpp
2int motif_len = length of motif (integer value above 30 is recommended, default value 30)
3int size_start = start value for the number of locally positioned dinucleotides (default value 10)
4int size_end = end value for the number of locally positioned dinucleotides (default value 90)
5int size_dif = variation value for the number of locally positioned dinucleotides (default value 10)
6double ratio_cnt_of_all(0=jk)  = ratio of the number of peaks to the number of control peaks (default value 10)
7int num_iterations = bumber of iteration in bootatrap (default 1), but it is recomended to run bootstrap several times to get reliable results

andy1_mat.cpp scans a fasta file with DNA sequences with a given model
command line arguments:
1file.seq  2sitename 3file_train 4thr 5cmpl 6file.ipr 7seq_head 8print_pos 9site_desc 10bit

1file.seq = test file
2sitename = file with sitega model
3file_train = faculatative file (default value train.fa)
4thr = threshold for  sitega model
5cmpl = default value 2 
6file.ipr = name for ouput files (default value chipseq)
7seq_head = (default value = 1)
8print_pos = (default value = 1)
9site_desc = (default value = 0)
10bit = (default value = 300)
