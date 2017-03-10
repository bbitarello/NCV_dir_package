#!/bin/sh

#  process_T1_T2.sh
#  
#
#  Created by Dutty on 4/14/16.
#
#Copied this from my github...bash_script_ballet.sh


cd /mnt/sequencedb/PopGen/barbara/BALLET/tmp_Eur/

for i in {1..1000};
do


awk '{print $2}' OutFile_T1_neu_$i |sort -rn |head -n 1  >> tmp_neu_T1.txt


#awk '{print $2}' OutFile_T1_f0.1_bs_$i |sort -rn |head -n 1  >> tmp_f0.1_bs_T1.txt


#awk '{print $2}' OutFile_T1_f0.2_bs_$i |sort -rn |head -n 1  >> tmp_f0.2_bs_T1.txt

awk '{print $2}' OutFile_T1_f0.3_bs_$i |sort -rn |head -n 1  >> tmp_f0.3_bs_T1.txt

awk '{print $2}' OutFile_T1_f0.4_bs_$i |sort -rn |head -n 1  >> tmp_f0.4_bs_T1.txt

awk '{print $2}' OutFile_T1_f0.5_bs_$i |sort -rn |head -n 1  >> tmp_f0.5_bs_T1.txt






awk '{print $2}' OutFile_T2_neu_$i |sort -rn|head -n 1  >> tmp_neu_T2.txt


#awk '{print $2}' OutFile_T2_f0.1_bs_$i |sort -rn|head -n 1  >> tmp_f0.1_bs_T2.txt

#awk '{print $2}' OutFile_T2_f0.2_bs_$i |sort -rn|head -n 1  >> tmp_f0.2_bs_T2.txt

awk '{print $2}' OutFile_T2_f0.3_bs_$i |sort -rn|head -n 1  >> tmp_f0.3_bs_T2.txt

awk '{print $2}' OutFile_T2_f0.4_bs_$i |sort -rn|head -n 1  >> tmp_f0.4_bs_T2.txt

awk '{print $2}' OutFile_T2_f0.5_bs_$i |sort -rn|head -n 1  >> tmp_f0.5_bs_T2.txt


done
