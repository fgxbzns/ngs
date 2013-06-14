#!/bin/bash
#total_reads=$(($(wc song_5_2_noPri_chr.sam)-24))
ngs_path="/home/guoxing/disk2/ngs/morehouse/python"
solid="/home/guoxing/disk2/solid/"
cd $solid
#pwd
declare -a chr=('5' 'X' '9' '1' '11' '7' '17' '13' 'X' '3' '11');
for((i=$1; i<$2; i++))
#for i in $1
do
#current_folder=$(pwd)"/song_$i/"
current_folder=$solid"song_$i"
cd $current_folder 
#rm primary.2010101* -r

chr_name="chr${chr[$i-1]}"
echo $chr_name

prem_file_name="song_"$i"_prem_"$chr_name"_sorted.sam"
#echo "grep $chr_name song_"$i"_prem_chr_sorted.sam > song_"$i"_prem_"$chr_name"_sorted.sam &"
#grep $chr_name "song_"$i"_prem_chr_sorted.sam > song_"$i"_prem_"$chr_name"_sorted.sam"
#$ngs_path"sam_process.py" -s $current_folder/"$prem_file_name" &"

prem_rmsk_file_name="song_"$i"_prem_"$chr_name"_sorted_rmsk.sam"
#echo $prem_rmsk_file_name
#mv hifi/$prem_rmsk_file_name prem_rmsk/ &

prem_indel_file_name="song_"$i"_prem_"$chr_name"_sorted_indel.sam"
#echo $prem_indel_file_name
#mv prem_indel/prem_rmsk/$prem_indel_file_name prem_indel/ &

prem_rmsk_indel_file_name="song_"$i"_prem_"$chr_name"_sorted_rmsk_indel.sam"
#echo $prem_rmsk_indel_file_name
#cp hifi/$prem_rmsk_indel_file_name prem_rmsk_indel/ &

#mkdir -p rmsk_ori
#mv *chr_sorted* rmsk_ori
#mv prem_rmsk rmsk_ori/
#mv -f prem_rmsk_indel/ rmsk_ori/


#ls

#cd prem_rmsk
#$ngs_path/"sam_process.py" -s "$prem_rmsk_file_name" &
#wait
#mv $prem_rmsk_indel_file_name ../ &
#wait
#cd ..

#cp prem/$prem_file_name ./ &
#ls prem_rmsk_indel/haplotype.txt || rm -r prem_rmsk_indel
#cd prem_rmsk_indel
#$ngs_path/solid_process_4.py -s $prem_rmsk_indel_file_name -c $chr_name -d 3 &
#$ngs_path/hifi &
#$ngs_path/hifiAccuCheck_v2.py -c $chr_name -n $i &

# for filter combination
#for dir_name in prem prem_rmsk prem_indel prem_rmsk_indel
for dir_name in prem_rmsk_indel
do
	#mkdir -p $dir_name
	#mv $prem_rmsk_indel_file_name $dir_name &
	#wait
	#cd $dir_name
	#ls
	fname=$(eval echo \${"$dir_name"_file_name})
	#$ngs_path/repeatRemove_sorted.py -s $prem_file_name &
	#wait
	#mv $prem_file_name ../
	#wait
	
	


	#$ngs_path/solid_process_4.py -s $fname -c $chr_name -d 3 &
	#wait
	#$ngs_path/hifi &
	#$ngs_path/hifiAccuCheck_v2.py -c $chr_name -n $1 &
	#wait
	#cd $current_folder
done


# for depth
#mkdir -p depth
#cp prem_rmsk_indel/$prem_rmsk_indel_file_name depth/ &
cd depth

for depth in {6..10}
do
	echo $depth
#	$ngs_path/solid_process_4.py -s "$prem_rmsk_indel_file_name" -c $chr_name -d $depth &
	wait
#	mkdir -p $depth
#	mv song_"$i"_prem_"$chr_name"_sorted_rmsk_indel_"$depth"_* $depth
	#mv hifi_*.txt imput* $depth &
	#mv infosum* match_* preimputehaplotype* quality_score* runingtime* $depth &
#	mv genotype.txt haplotype.txt refHaplos.txt $depth &
	
	cd $depth
	$ngs_path/hifi &
	wait
	#$ngs_path/hifiAccuCheck_v2.py -c $chr_name &
	wait
	cd ..
	#mv genotype.txt haplotype.txt refHaplos.txt $depth &
	#mv hifi_*.txt imput* $depth &
	#mv infosum* match_* preimputehaplotype* quality_score* runingtime* $depth &
done



#cd prem 
#$ngs_path/solid_process_4.py -s "$prem_file_name" -c $chr_name &
#wait
#$ngs_path/hifi &
#$ngs_path/hifiAccuCheck_v2.py -c $chr_name &
#wait
#cd $current_folder

#cd prem_rmsk 
#$ngs_path/solid_process_4.py -s "$prem_rmsk_file_name" -c $chr_name &
#wait
#$ngs_path/hifi &
#$ngs_path/hifiAccuCheck_v2.py -c $chr_name &
#wait
#cd $current_folder


#cd prem_indel 
#$ngs_path/solid_process_4.py -s "$prem_indel_file_name" -c $chr_name &
#wait
#$ngs_path/hifi &
#$ngs_path/hifiAccuCheck_v2.py -c $chr_name &
#wait
#cd $current_folder


#cd prem_rmsk_indel 
#$ngs_path/solid_process_4.py -s "$prem_rmsk_indel_file_name" -c $chr_name &
#wait
#$ngs_path/hifi &
#$ngs_path/hifiAccuCheck_v2.py -c $chr_name &
#wait
#cd $current_folder

cd $solid



done
exit 0
