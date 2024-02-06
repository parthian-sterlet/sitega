#!/bin/sh
cd src
g++ -o sitega_thr_dist_mat.exe sitega_thr_dist_mat.cpp
g++ -o andy0bsn5.exe andy0bsn5.cpp
g++ -o andy1_mat.exe andy1_mat.cpp
g++ -o andy05.exe andy05.cpp
chmod a+x andy0bsn5.exe
chmod a+x andy1_mat.exe
chmod a+x andy05.exe
chmod a+x sitega_thr_dist_mat.exe
chmod a+x bootstrap
chmod a+x scan 
chmod a+x train
cd ..
cd genomes
cd mm10
tar -cvzf ups2kb_mm10.seq.tar.gz ups2kb_mm10.seq
cd ..
cd ..
