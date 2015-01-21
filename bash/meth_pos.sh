#!/bin/bash

#pos_file='/home/guoxing/storage1/lima/12878_chr8_meth/meth_pos_all.txt'
pos_file='/home/guoxing/storage1/lima/12878_chr8_meth/pos_all_jan17.txt'
ngs_path="/home/guoxing/disk2/ngs/morehouse/python/"
geno_file=$1
echo $geno_file

#$ngs_path"snpPick_all.py" -d 'chr8_c_indel_sorted_qs'$qs'.txt' -c chr8 -m output

gzip -c $geno_file > $geno_file'.gz'
zfgrep -w --file=$pos_file $geno_file'.gz' > $geno_file'_cnv'
cut -f1 $geno_file'_cnv' | diff $pos_file - > $geno_file'_diff'
rm $geno_file'.gz'