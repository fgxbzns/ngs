#!/bin/bash
ngs_path="/home/guoxing/disk2/ngs/morehouse/python/"

# for hg18 solid mimi data
project_path='/home/guoxing/storage1/zebrafish/'
sample_name="A2"
sample_path=$project_path'A2'
rmsk_path=$project_path'rmsk'

#sam_path="/home/guoxing/disk2/lima/mimi_solid/mimi_solid_sam/"
#db_path="/home/guoxing/disk2/lima/mimi_solid/mimi_solid_snpPick_db/"

chr_3_length=63268876

declare -a chr_length=('xxx' '60300536' '63268876' 'xxx' 'xxx' 'xxx' 'xxx' 'xxx' 'xxx' 'xxx' 'xxx');

for((i=17; i<=20; i++))

	do
	cd $sample_path
	chr_name="chr"$i
	sam_file=$chr_name".sam"
	folder_name="chr"$i

	mkdir -p $folder_name
	#mv $sam_file $folder_name
	cd $folder_name
	echo "processing " $chr_name

	# repeat remove
	#$ngs_path"sam_process.py" -s $sam_file -c $chr_name -m fish_wli > sam_process_record.txt

	# variation call
	processed_sam_file=$chr_name"_pairend_XA_sorted_rmsk_combined_indel.sam"
	echo $processed_sam_file
	db_file=$sample_name"_"$chr_name"_qs30"
	$ngs_path"snpPick_fish_sql_mimi_wli.py" -s $processed_sam_file -c $chr_name -m update -d $db_file > variation_call_record.txt &

	# output variation call
	#snpPick_fish -c chr -m mf -b startLine -e endLine -d db_name
	#$ngs_path"snpPick_fish_sql_mimi_wli.py" -c $chr_name -m mf -b 0 -e $chr_3_length -d $db_file > mimi_output_record.txt


	#wait
	#echo "processing finished " $chr_name
	cd ..
	done



: '
# to sort chr in rmsk file by starting position
cd $rmsk_path

for((i=1; i<=25; i++))
	do

	file_name="zebrafish_rmsk_danRer7.txt"
	chr_name="chr"$i
	grep_name="rmsk_"$chr_name".txt"
	echo "greping"
	#grep -w $chr_name $file_name > $grep_name
	echo "sorting"
	sorted_name="rmsk_"$chr_name"_sorted.txt"
	#sort -k 7 -n $grep_name > $sorted_name
	echo "cat "$chr_name
	cat rmsk_sorted.txt $sorted_name > rmsk_sorted_tmp.txt
	mv rmsk_sorted_tmp.txt rmsk_sorted.txt
	#rm $grep_name
	done
'
: '
# to concanaet sorted chr into one file.
sorted_sam_path=$project_path"sorted_chr"
cd $sorted_sam_path
# need to remove empty line in rmsk_sorted file. Maybe combine 1 and 2 first. or rename chr1 to tmp
for((i=1; i<=25; i++))
	do

	chr_name="chr"$i

	sorted_name=$chr_name".sam"
	echo "cat "$chr_name
	cat rmsk_sorted.txt $sorted_name > rmsk_sorted_tmp.txt
	mv rmsk_sorted_tmp.txt rmsk_sorted.txt
	done
'