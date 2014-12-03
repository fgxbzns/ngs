#!/bin/bash

#link='http://hapmap.ncbi.nlm.nih.gov/downloads/frequencies/latest_phaseIII_ncbi_b36/fwd_strand/non-redundant/'
#link='http://hapmap.ncbi.nlm.nih.gov/downloads/ld_data/latest/'
link='ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20130502/'

population='CHB'
#declare -a pos_list=('10316582' '117223133' '172751522' '78130319')
for i in {1..22}
#for i in 1
	do
		#filename='allele_freqs_chr'$i'_'$population'_phase3.2_nr.b36_fwd.txt.gz'
		#filename='ld_chr'$i'_'$population'.txt.gz'
		#filename='ALL.chr'$i'.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz'
		filename='chb_chs_chr'$i'.vcf.gz'
		#echo "downloading " $link$filename
		#wget $link$filename
		echo "getting columns " $filename
		#time gzip -cd $filename | cut -f1-9,1784-1886 | gzip > 'chb_chr'$i'.vcf.gz'
		#time gzip -cd $filename | cut -f1-9,195-306,1784-1886 | gzip > 'chb_chs_chr'$i'.vcf.gz'
		#mv 'chb_chs_chr'$i'.vcf.gz' chb_chs_1000g
		#filename='chb_chs_chr'$i'.vcf.gz'
		zfgrep -w --file=rsID.txt $filename >> chr_rsID.txt






	done


#frequency   http://hapmap.ncbi.nlm.nih.gov/downloads/frequencies/latest_phaseIII_ncbi_b36/fwd_strand/non-redundant/allele_freqs_chr1_ASW_phase3.2_nr.b36_fwd.txt.gz
#ld data   http://hapmap.ncbi.nlm.nih.gov/downloads/ld_data/latest/ld_chr1_CHB.txt.gz
#genotype   ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20130502/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz