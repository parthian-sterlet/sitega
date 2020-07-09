# sitega
folder src contains four files with sitega source codes, c++ language

andy02.cpp train a model with a given train ChIP-seq dataset (peaks)
1char file_cor 2int motif_len 3int size_start 4int size_end 5int size_dif

andy1_mat.cpp scan a fasta file with DNA sequences with a given model

andy0bsn2.cpp perform bootsrap cross validation test to estimate the performance of a model with a given train ChIP-seq dataset
1char file_cor 2int motif_len 3int size_start 4int size_end 5int size_dif 6double ratio_cnt_of_all(0=jk) 7int num_iterations 
