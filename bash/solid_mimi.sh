#!/bin/bash
ngs_path="/home/guoxing/disk2/ngs/morehouse/python/"

# for hg18 solid mimi data
sam_path="/home/guoxing/disk2/lima/mimi_solid/mimi_solid_sam/"
db_path="/home/guoxing/disk2/lima/mimi_solid/mimi_solid_snpPick_db/"
chr_9_length=140273252

declare -a chr_length=('180857866' '154913754' '140273252' '247249719' '134452384' '158821424' '78774742' '114142980' '154913754' '199501827' '134452384');

#cd $sam_path
cd $db_path

declare -a chr=('5' 'X' '9' '1' '11' '7' '17' '13' 'X' '3' '11');
: '
data_path="/home/guoxing/disk2/solid/"
cd $data_path

for((i=4; i<=11; i++))
	do
	folder_name="song_"$i
	cd $folder_name
	chr_name="chr${chr[$i-1]}"
	echo $chr_name
	file_name="song_"$i"_prem_"$chr_name"_sorted.sam"
	echo $file_name
	#cp $file_name $sam_path &
	cd ..
	done
'
for((i=5; i<=11; i++))

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
	target_chr_length=${chr_length[$i-1]}
	echo $chr_name $target_chr_length
	#snpPick_fish -c chr -m mf -b startLine -e endLine -d db_name
	$ngs_path"snpPick_fish_sql_lima.py" -c $chr_name -m mf -b 0 -e $target_chr_length -d $db_file &


	wait
	cd ..
	done