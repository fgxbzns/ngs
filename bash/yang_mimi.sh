#!/bin/bash
ngs_path="/home/guoxing/disk2/ngs/morehouse/python/"

# for hg18 solid mimi data
#project_path="/home/guoxing/disk2/lima/yang/"
#sam_path="/home/guoxing/disk2/lima/yang/mimi_yang_sam/"
#db_path="/home/guoxing/disk2/lima/yang/mimi_yang_snpPick_db/"


project_path="/home/guoxing/storage1/lima/yang/"
sam_path=$project_path"mimi_yang_sam/"
db_path=$project_path"mimi_yang_snpPick/"

chr19_length='63811651'
#cd $project_path
cd $sam_path
#cd $db_path


chr_name="chr19"

#for((i=$1; i<=$2; i++))
for((i=1; i<=24; i++))
	do
	fastq_file="tag_"$i".fastq"
	folder_name="tag_"$i
	#mkdir $folder_name
	#mv $fastq_file $folder_name
	#cd $folder_name
	#$ngs_path"primerRemove.py" -i $fastq_file

	seqRem_file_name="tag_"$i"_seqRem.fastq"
	echo $seqRem_file_name
	#$ngs_path"bwaRun.py" -i $seqRem_file_name &

	sam_file_name="tag_"$i"_seqRem_all.sam"
	sam_chr_file_name="tag_"$i"_seqRem_all_chr.sam"
	#grep "chr" $sam_file_name > sam_chr_file_name

	#$ngs_path"chrPercentage.py" -i $sam_chr_file_name
	#ls -lh
	#chr_percentage_file="tag_"$i"_seqRem_all_chr.sam_percentage"
	#head $chr_percentage_file

	#sam_process
	#$ngs_path"sam_process.py" -s $sam_chr_file_name -c $chr_name -m yang_mimi

	# variation call
	processed_sam_file="tag_"$i"_seqRem_all_"$chr_name"_XA_sorted_rmsk_single_indel.sam"
	echo $processed_sam_file
	#db_file="tag_"$i"_"$chr_name"_qs20"
	$ngs_path"snpPick_all.py" -s $processed_sam_file -c $chr_name -m yang

	# output variation call
	#snpPick_fish -c chr -m mf -b startLine -e endLine -d db_name
	#$ngs_path"snpPick_fish_sql_lima.py" -c $chr_name -m mf -b 0 -e $chr19_length -d $db_file &

	wait
	#cd ..
	done

mv *.txt $db_path

: '
# to sort chr in rmsk file by starting position
rmsk_path="/home/guoxing/disk2/solid/common_files/hg18_rmsh_chr"
cd $rmsk_path

for((i=1; i<=25; i++))
	do

	file_name="hg18_rmsk.txt_original"
	chr_name="chr"$i
	grep_name="rmsk_"$chr_name".txt"
	echo "greping"
	grep -w $chr_name $file_name > $grep_name
	#echo "sorting"
	#sorted_name="rmsk_"$chr_name"_sorted.txt"
	#sort -k 7 -n $grep_name > $sorted_name
	#echo "cat "$chr_name
	#cat rmsk_sorted.txt $sorted_name > rmsk_sorted_tmp.txt
	#mv rmsk_sorted_tmp.txt rmsk_sorted.txt
	#rm $grep_name
	done
'

: '
for((i=11; i<=11; i++))

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
	#$ngs_path"snpPick_fish_sql_lima.py" -c $chr_name -m mf -b 0 -e $target_chr_length -d $db_file &

	# map rmsk file
	#repeatRemove_sorted -s song_8_chr13_qs30_0_114142980_filtered_2nddepth_5.txt -r a -c chr13 -m mimi_map_rmsk
	pos_file="song_"$i"_"$chr_name"_qs30_0_"$target_chr_length"_filtered_2nddepth_5.txt"
	$ngs_path"repeatRemove_sorted.py" -s $pos_file -c $chr_name -r a -m mimi_map_rmsk &

	wait
	#cd ..
	done
'
