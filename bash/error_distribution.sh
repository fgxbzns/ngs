#!/bin/bash
#total_reads=$(($(wc song_5_2_noPri_chr.sam)-24))
ngs_path="/home/guoxing/disk2/ngs/morehouse/python"
solid="/home/guoxing/disk2/solid/"
cd $solid
declare -a chr=('5' 'X' '9' '1' '11' '7' '17' '13' 'X' '3' '11');
for((i=$1; i<$2; i++))
do

current_folder=$solid"song_$i"
cd $current_folder 
mkdir -p error_dirstribution

#chr_name="chr${chr[$i-1]}"
#echo $chr_name

#seed_file="song_"$i"_prem_"$chr_name"_sorted_rmsk_indel_3_hifi.txt"
error_position_file="song_"$i"_prem_"$chr_name"_sorted_rmsk_indel_3_B.txt"
#cp prem_rnsk_indel/$seed_file error_dirstribution/ &
cp prem_rnsk_indel/$error_position_file error_dirstribution/ &
cp haplotype.txt genotype.txt refHaplos.txt error_dirstribution/ &

cd $solid



done
exit 0
