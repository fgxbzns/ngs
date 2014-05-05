#!/bin/bash
ngs_path="/home/guoxing/disk2/ngs/morehouse/python/"

# for hg18 solid mimi data
project_path='/home/guoxing/disk2/wli/'
rmsk_path=$project_path'rmsk'

sam_path="/home/guoxing/disk2/lima/mimi_solid/mimi_solid_sam/"
db_path="/home/guoxing/disk2/lima/mimi_solid/mimi_solid_snpPick_db/"
chr_9_length=140273252

#cd $sam_path

# to sort chr in rmsk file by starting position
cd $rmsk_path

for((i=1; i<=25; i++))
	do

	file_name='zebrafish_rmsk_danRer7.txt'
	chr_name='chr'$i
	grep_name='rmsk_'$chr_name'.txt'
	echo 'greping'
	#grep -w $chr_name $file_name > $grep_name
	echo 'sorting'
	sorted_name='rmsk_'$chr_name'_sorted.txt'
	#sort -k 7 -n $grep_name > $sorted_name
	echo 'cat '$chr_name
	cat rmsk_sorted.txt $sorted_name > rmsk_sorted_tmp.txt
	mv rmsk_sorted_tmp.txt rmsk_sorted.txt
	#rm $grep_name
	done

: '
for((i=3; i<=3; i++))

	do
	chr_name="chr${chr[$i-1]}"
	sam_file="song_"$i"_prem_"$chr_name"_sorted.sam"
	folder_name="song_"$i

	#mkdir $folder_name
	#mv $sam_file $folder_name
	#cd $folder_name

	# repeat remove
	#$ngs_path"sam_process.py" -s $sam_file -c $chr_name -m solid_lima &

	# variation call
	processed_sam_file="song_"$i"_prem_"$chr_name"_sorted_rmsk_single_indel_base_cleaned.sam"
	echo $processed_sam_file
	db_file="song_"$i"_"$chr_name"_qs30"
	#$ngs_path"snpPick_fish_sql_lima.py" -s $processed_sam_file -c $chr_name -m update -d $db_file &

	# output variation call
	#snpPick_fish -c chr -m mf -b startLine -e endLine -d db_name
	$ngs_path"snpPick_fish_sql_lima.py" -c $chr_name -m mf -b 0 -e $chr_9_length -d $db_file &


	#wait
	cd ..
	done
'