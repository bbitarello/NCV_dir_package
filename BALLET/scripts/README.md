as modified: 26.11.2014
##Barbara Bitarello
#for simulations


#open R in this directory and run the script 


barbara_sims_T1_T2.R

This generates the input files

In this directory, run the commands below as a test



#DivFile
./BALLET -inter_coal_time CombinedSNPFile_neu 15000000 0.000731 DivFile_neu


#coalescence time is greater for BS, so so far it looks better

#PolySub
./BALLET -poly_sub CombinedSNPFile_neu PolySubFile_neu



#Spect
./BALLET -spect CombinedSNPFile_neu SpectFile_neu




# run 4 jobs (t1 /t2 , neu, bs)

qsub -o ballet_T1_neu.log run_ballet_T1_neu.sge 
qsub -o ballet_T1_bs.log run_ballet_T1_bs.sge 
qsub -o ballet_T2_bs.log run_ballet_T2_bs.sge 
qsub -o ballet_T2_neu.log run_ballet_T2_neu.sge 



#after running everything, run 

bash_script_ballet.sh


#to get one output file for each of the sets, and compare the distributions of T/T2 for both.

#after runnign everything I moved these files to the 'tmp' directory.

