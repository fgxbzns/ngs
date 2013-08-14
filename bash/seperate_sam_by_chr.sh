#!/bin/bash
##############################################################
# give the sam file name in parameters and the content in the 
# sam file will be seperated by chr names. chr_name.sam files 
# will be generated
##############################################################
for i in "$@";do
	echo "processing $i"
	awk '{if(substr($1,0,1)!="@" && $3!="*" && $7=="="){print $0 >> $3 ".sam"};if(substr($1,0,1)!="@" && $3!="*" && $7!="*" &&$7!="="){print $0 >> "ctx.sam"}}' $i
done

