#!/bin/bash
ngs_path="/home/guoxing/disk2/ngs/morehouse/python/"
sam_path="/home/guoxing/disk2/lima/mimi_chrx_sam/"
db_path="/home/guoxing/disk2/lima/mimi_snpPick_db/"
chrX_length=158375978

#cd $sam_path
cd $db_path
for((i=84; i<=92; i++))

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
	$ngs_path"snpPick_fish_sql_lima.py" -c chrX -m mf -b 0 -e $chrX_length -d $db_file &


	wait
	#cd ..
	done