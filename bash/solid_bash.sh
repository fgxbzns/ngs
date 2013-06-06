#!/bin/bash

ngs_path="/home/guoxing/disk2/ngs/morehouse/python/"
solid="/home/guoxing/disk2/solid/"
cd $solid
pwd
declare -a chr=('5' 'X' '9' '1' '11' '7' '17' '13' 'X' '3');
for((i=$1; i<$2; i++))
#for i in $1
do
#current_folder=$(pwd)"/song_$i/"
current_folder="song_$i"
echo $current_folder
#echo "current folder $current_folder"
cd $current_folder 
#rm primary.2010101* -r
#mkdir prem
#mkdir rmsk 
#mkdir prem_indel
#mkdir prem_rmsk_indel
chr_name="chr${chr[$i-1]}"
echo $chr_name

prem_file_name="song_$i_prem_$chr_name_sorted.sam"
#echo "grep $chr_name song_"$i"_prem_chr_sorted.sam > song_"$i"_prem_"$chr_name"_sorted.sam &"
#grep $chr_name "song_"$i"_prem_chr_sorted.sam > song_"$i"_prem_"$chr_name"_sorted.sam"
#$ngs_path"sam_process.py" -s $current_folder/"$prem_file_name" &"

#mv $prem_file_name prem/ &

prem_rmsk_file_name="song_$i_prem_$chr_name_sorted_indel.sam"
echo $prem_rmsk_file_name
mv $prem_rmsk_file_name prem_rmsk/ &


ls *indel.sam
cd $solid



done
exit 0
