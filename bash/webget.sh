#!/bin/bash

#link='http://hapmap.ncbi.nlm.nih.gov/downloads/frequencies/latest_phaseIII_ncbi_b36/fwd_strand/non-redundant/'
link='http://hapmap.ncbi.nlm.nih.gov/downloads/ld_data/latest/'


population='CHB'
#declare -a pos_list=('10316582' '117223133' '172751522' '78130319')
for i in {2..22} X
	do
		#filename='allele_freqs_chr'$i'_'$population'_phase3.2_nr.b36_fwd.txt.gz'
		filename='ld_chr'$i'_'$population'.txt.gz'
		echo "downloading " $link$filename
		wget $link$filename
	done


#http://hapmap.ncbi.nlm.nih.gov/downloads/frequencies/latest_phaseIII_ncbi_b36/fwd_strand/non-redundant/allele_freqs_chr1_ASW_phase3.2_nr.b36_fwd.txt.gz
#http://hapmap.ncbi.nlm.nih.gov/downloads/ld_data/latest/ld_chr1_CHB.txt.gz
