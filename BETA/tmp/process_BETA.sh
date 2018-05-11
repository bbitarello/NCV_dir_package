rm tmp_*


for i in {1..1000};
do


awk '{print $2}' Betas_neuSNPFile$i |sort -rn |head -n 1    >> tmp_neu.txt

awk '{print $2}' Betas_Eu_neuSNPFile$i |sort -rn |head -n 1    >> tmp_Eu_neu.txt

awk '{print $2}' Betas_As_neuSNPFile$i |sort -rn |head -n 1    >> tmp_As_neu.txt

awk '{print $2}' Betas_bs_f0.5_SNPFile$i | sort -rn |head -n 1    >> tmp_bs_f0.5.txt 

awk '{print $2}' Betas_Eu_bs_f0.5_SNPFile$i | sort -rn |head -n 1    >> tmp_Eu_bs_f0.5.txt 

awk '{print $2}' Betas_As_bs_f0.5_SNPFile$i |sort -rn |head -n 1    >> tmp_As_bs_f0.5.txt 

awk '{print $2}' Betas_bs_f0.4_SNPFile$i |sort -rn |head -n 1    >> tmp_bs_f0.4.txt

awk '{print $2}' Betas_Eu_bs_f0.4_SNPFile$i |sort -rn |head -n 1    >> tmp_Eu_bs_f0.4.txt

awk '{print $2}' Betas_As_bs_f0.4_SNPFile$i |sort -rn |head -n 1    >> tmp_As_bs_f0.4.txt 

awk '{print $2}' Betas_bs_f0.3_SNPFile$i |sort -rn |head -n 1    >> tmp_bs_f0.3.txt

awk '{print $2}' Betas_Eu_bs_f0.3_SNPFile$i |sort -rn |head -n 1    >> tmp_Eu_bs_f0.3.txt

awk '{print $2}' Betas_As_bs_f0.3_SNPFile$i |sort -rn |head -n 1    >> tmp_As_bs_f0.3.txt 

echo $i
done
