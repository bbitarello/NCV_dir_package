./make_beta_input.R  #make input files for beta

cd tmp/

./run_BETA.sh  #run beta. creates a lot of files


./process_BETA.sh #takes hoghest beta value from each sismulation

./power_analyses_BETA.R #check power
