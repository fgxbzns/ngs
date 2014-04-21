#!/bin/bash
ngs_path="/home/guoxing/disk2/ngs/morehouse/python/"

# for hg18 solid mimi data
sam_path="/home/guoxing/disk2/lima/solid_mimi/mimi_solid_sam/"
db_path="/home/guoxing/disk2/lima/solid_mimi/mimi_solid_snpPick_db/"
chr_9_length=140273252

#cd $sam_path
#cd $db_path

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
for((i=78; i<=79; i++))

	do
	#sam_file="NA128"$i"_S1_ChrXnew.sam"
	sam_file="NA128"$i"_S1_chrX.sam"
	#mv $sam_file $sam_file_new
	folder_name="NA128"$i

	#mkdir $folder_name
	#mv $sam_file $folder_name
	#cd $folder_name

	# pre-process
	#$ngs_path"sam_process.py" -s $sam_file -c chrX -m mimi &

	# variation call
	processed_sam_file="NA128"$i"_S1_chrX_pairend_XA_sorted_rmsk_combined_indel.sam"
	db_file="NA128"$i"_S1_chrX"
	#$ngs_path"snpPick_fish_sql_lima.py" -s $processed_sam_file -c chrX -m update -d $db_file &

	# output variation call
	#snpPick_fish -c chr -m mf -b startLine -e endLine -d db_name
	#$ngs_path"snpPick_fish_sql_lima.py" -c chrX -m mf -b 0 -e $chrX_length -d $db_file &


	#wait
	#cd ..
	done