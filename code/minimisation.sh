#This script minimizes Ala10. 

#----------------
#Get args : 
res_dir=$1
echo $res_dir

#----------------
#Get .tpr file, a portable binary run input file containing the starting structure (.pdb) of the simulation, the molecular topology (.top) and all the simulation parameters (.tpr):

#To ignore these warnings : 
#"
#Warning: atom name 1 in ala10.top and ala10_md000.pdb does not match (CH3 - HH31)
#Warning: atom name 2 in ala10.top and ala10_md000.pdb does not match (HH31 - CH3)
#If you are sure all warnings are harmless, use the -maxwarn option.
#"                                                  --------
#XXX TODO : should be ignored?
gmx grompp -f ../data/mini2.mdp -c ../data/ala10_md000.pdb -p ../data/ala10.top -po ../data/grompp_mdout.mdp -o ../data/minimi.tpr -maxwarn 1

#read new .tpr file : 
#gmx dump -s topol.tpr | more

#----------------
#Run minimi : 
gmx mdrun -v -s ../data/minimi.tpr -o $res_dir/traj_minimi.trr -c ../data/ala10_minimized.gro -e $res_dir/energy_minimi.edr -g $res_dir/minimi.log

#----------------
#Check energy diminution on the plot: 

echo "10" | gmx energy -f $res_dir/energy_minimi.edr -o $res_dir/energy_minimi.xvg
#(choose Potential)

xmgrace $res_dir/energy_minimi.xvg
