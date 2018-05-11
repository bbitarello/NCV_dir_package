# run 4 jobs (t1 /t2 , neu, bs)
	
qsub -o ballet_T1_neu.log run_ballet_T1_neu.sge 
qsub -o ballet_T1_bs.log run_ballet_T1_bs.sge 
qsub -o ballet_T2_bs.log run_ballet_T2_bs.sge 
qsub -o ballet_T2_neu.log run_ballet_T2_neu.sge 
