#This script minimizes Ala10. 

#----------------
#Create res dir : 
res_dir=../results/minimization-$(date '+%Y-%b-%d-%Hh-%Mmin')
mkdir $res_dir

#----------------
#Get .tpr file, a portable binary run input file containing the starting structure (.pdb) of the simulation, the molecular topology (.top) and all the simulation parameters (.tpr):
gmx grompp -f ../data/mini2.mdp -c ../data/ala10_md000.pdb -p ../data/ala10.top -po ../data/grompp_mdout.mdp -o ../data/minimi.tpr


