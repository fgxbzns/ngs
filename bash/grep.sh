#!/bin/bash
file_name=$1
#declare -a pos_list=('10316582' '117223133' '172751522' '78130319')
for pos in '10316582' '117223133' '172751522' '78130319'
do
grep $pos $file_name
done
