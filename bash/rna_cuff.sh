#!/bin/bash
ngs_path="/home/guoxing/disk2/ngs/morehouse/python"

project_path="/home/guoxing/disk3/rna_seq"
cd $project_path
echo running cufflinks
tophat_output=$project_path"/tophat"
: '
for ((i=1; i<=1; i++))

	do
	echo "processing sample "$i
	output_folder=$project_path"/cufflinks/"$i
	combined_sam=$output_folder"/"$i"_tophat_accepted.sam"
	sorted_sam=$output_folder"/"$i"_tophat_accepted_sorted.sam"
	sample_folder=$tophat_output"/"$i

	for ((j=2; j<=8; j++))
		do
		bam_file=$sample_folder"/tophat_out_L00"$j"/accepted_hits.bam"
		echo "samtools processing "$bam_file
		samtools view $bam_file >> $combined_sam
		wait
		done
	echo "sorting started"
	sort -k 3,3 -k 4,4n $combined_sam > $sorted_sam
	wait
	echo "sorting finished"
	rm $combined_sam
	wait
	done
echo "All done"
'

for ((i=1; i<=24; i++))

	do
	echo "processing sample "$i
	output_folder=$project_path"/cufflinks/"$i
	cd $output_folder
	#sorted_sam=$output_folder"/"$i"_tophat_accepted_sorted.sam"
	#echo processing $sorted_sam
	#cufflinks $sorted_sam &
	$ngs_path/rnaseq_analyze.py -f transcripts.gtf
	wait
	mv transcripts_statistics.txt $project_path/cufflinks/statistics/"transcripts_statistics_"$i".txt"

	#folder=$tophat_output"/"$i
	#cd $folder
	#rm -r $folder/logs


	echo finished processing $i
	#cp align_summary.txt $tophat_output/summary/"align_summary_"$i".txt"

	done
echo "All done"