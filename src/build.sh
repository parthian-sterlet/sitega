#!/bin/sh

g++ -o sitega_thr_dist_mat.exe sitega_thr_dist_mat.cpp
g++ -o andy0bsn5cell.exe andy0bsn5cell.cpp
g++ -o andy1_mat.exe andy1_mat.cpp
g++ -o andy05cell.exe andy05cell.cpp
g++ -o fasta_to_plain0.exe fasta_to_plain0.cpp
g++ -o fasta_muliplefiles.exe fasta_muliplefiles.cpp

chmod a+x andy0bsn5cell.exe
chmod a+x andy1_mat.exe
chmod a+x andy05cell.exe
chmod a+x sitega_thr_dist_mat.exe
chmod a+x fasta_to_plain0.exe 
chmod a+x fasta_muliplefiles.exe

chmod a+x bootstrap
chmod a+x scan_fasta
chmod a+x scan_genome
chmod a+x thr_err
chmod a+x train

cd ../genomes/mm10
tar -xzvf ups2kb_mm10.seq.tar.gz
cd ..
cd ..

cd bin
chmod a+x andy0bsn5cell.exe
chmod a+x andy1_mat.exe
chmod a+x andy05cell.exe
chmod a+x sitega_thr_dist_mat.exe
chmod a+x fasta_to_plain0.exe 
chmod a+x fasta_muliplefiles.exe
cd ..
