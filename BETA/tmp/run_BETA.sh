##############################
#	Barbara Bitarello
#
#	11.Sept.2017
##############################


#for i in {1..1000};
#do

#sed -i '1d' /mnt/sequencedb/PopGen/barbara/NCV_dir_package/BETA/tmp/neuSNPFile$i 
#sed -i '1d' /mnt/sequencedb/PopGen/barbara/NCV_dir_package/BETA/tmp/bs_f0.5_SNPFile$i

#sed -i '1d' /mnt/sequencedb/PopGen/barbara/NCV_dir_package/BETA/tmp/bs_f0.4_SNPFile$i
#sed -i '1d' /mnt/sequencedb/PopGen/barbara/NCV_dir_package/BETA/tmp/bs_f0.3_SNPFile$i

#done
#each time I change sometihing like window size, I save the R object, but the tmp files are overwritten

for i in {1..1000};
do

python /mnt/sequencedb/PopGen/barbara/NCV_dir_package/BETA/BetaScan.py -i /mnt/sequencedb/PopGen/barbara/NCV_dir_package/BETA/tmp/neuSNPFile$i -w 3000 -p 20 -fold

python /mnt/sequencedb/PopGen/barbara/NCV_dir_package/BETA/BetaScan.py -i /mnt/sequencedb/PopGen/barbara/NCV_dir_package/BETA/tmp/bs_f0.5_SNPFile$i -w 3000 -p 20 -fold

python /mnt/sequencedb/PopGen/barbara/NCV_dir_package/BETA/BetaScan.py -i /mnt/sequencedb/PopGen/barbara/NCV_dir_package/BETA/tmp/bs_f0.4_SNPFile$i -w 3000 -p 20 -fold

python /mnt/sequencedb/PopGen/barbara/NCV_dir_package/BETA/BetaScan.py -i /mnt/sequencedb/PopGen/barbara/NCV_dir_package/BETA/tmp/bs_f0.3_SNPFile$i -w 3000 -p 20 -fold

python /mnt/sequencedb/PopGen/barbara/NCV_dir_package/BETA/BetaScan.py -i /mnt/sequencedb/PopGen/barbara/NCV_dir_package/BETA/tmp/Eu_neuSNPFile$i  -w 3000 -p 20 -fold

python /mnt/sequencedb/PopGen/barbara/NCV_dir_package/BETA/BetaScan.py -i /mnt/sequencedb/PopGen/barbara/NCV_dir_package/BETA/tmp/Eu_bs_f0.5_SNPFile$i -w 3000 -p 20 -fold

python /mnt/sequencedb/PopGen/barbara/NCV_dir_package/BETA/BetaScan.py -i /mnt/sequencedb/PopGen/barbara/NCV_dir_package/BETA/tmp/Eu_bs_f0.4_SNPFile$i -w 3000 -p 20 -fold

python /mnt/sequencedb/PopGen/barbara/NCV_dir_package/BETA/BetaScan.py -i /mnt/sequencedb/PopGen/barbara/NCV_dir_package/BETA/tmp/Eu_bs_f0.3_SNPFile$i -w 3000 -p 20 -fold

python /mnt/sequencedb/PopGen/barbara/NCV_dir_package/BETA/BetaScan.py -i /mnt/sequencedb/PopGen/barbara/NCV_dir_package/BETA/tmp/As_neuSNPFile$i -w 3000 -p 20 -fold

python /mnt/sequencedb/PopGen/barbara/NCV_dir_package/BETA/BetaScan.py -i /mnt/sequencedb/PopGen/barbara/NCV_dir_package/BETA/tmp/As_bs_f0.5_SNPFile$i -w 3000 -p 20 -fold

python /mnt/sequencedb/PopGen/barbara/NCV_dir_package/BETA/BetaScan.py -i /mnt/sequencedb/PopGen/barbara/NCV_dir_package/BETA/tmp/As_bs_f0.4_SNPFile$i -w 3000 -p 20 -fold

python /mnt/sequencedb/PopGen/barbara/NCV_dir_package/BETA/BetaScan.py -i  /mnt/sequencedb/PopGen/barbara/NCV_dir_package/BETA/tmp/As_bs_f0.3_SNPFile$i -w 3000 -p 20 -fold

echo $i
done
