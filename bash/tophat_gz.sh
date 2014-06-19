#!/b{i}n/bash

project_path="/home/guoxing/disk3/rna_seq"
cd $project_path
echo tophat mapping start...
fastq_gz="/home/guoxing/disk3/bcl2fastq_06122014"
for ((i=15; i<=15; i++))
 do 
  folder=$fastq_gz"/Project_BCM-3/Sample_15"
  cd $folder
  for ((j=2; j<=8; j++))
   do
    echo "sample${i} lane${j} processing"
    tag="TCCGGAGA-CAGGACGT"
    r1=$i"_"$tag"_L00"$j"_R1_00"
    r2=$i"_"$tag"_L00"$j"_R2_00"
    output_folder=$project_path"/tophat_output/"$i
    tophat2 -r -20 -o $output_folder/tophat_out_L00${j} $project_path/sugar_beet_ref/sugar_beet_ref $r1"1.fastq.gz",$r1"2.fastq.gz",$r1"3.fastq.gz",$r1"4.fastq.gz" $r2"1.fastq.gz",$r2"2.fastq.gz",$r2"3.fastq.gz",$r2"4.fastq.gz" &
    echo "sample${i} lane${j} finished"
   done
  wait
 done 
echo "All done"
