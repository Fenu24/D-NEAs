#/bin/bash

############################################################
#      AUTHOR: Marco Fenucci                               #
#        DATE: Agust 2022                                  #
# INSTITUTION: University of Belgrade                      #
#                                                          #
# This script lanches the thermal inertia estimation in    #
# parallel for a list of objects, each of them with 1 CPU  #
#                                                          #
############################################################

# TODO: This script must be improved:
#  1) Input distributions of the objects must be placed in a folder
#     called data
#  2) Instead of running all the objects in data, we should provide 
#     a list of objects to run, maybe including their semimajor axis and eccentricity
#  3) Instead of using fixed input files for the main program, I think we should give
#     in input other arguments to choose the options:
#        -nproc=        (mandatory)
#        -C=            (optional)
#        -Kmin=         (optional)
#        -Kmax=         (optional)
#        -method=       (optional)
#        -max_iter=     (optional)
#        -expo=         (optional)
#        -epsi=         (optional)
#        -alpha=        (optional)
# For the optional arguments, we should assign some default value


function usage () {
   echo " Usage: gamma_est_mc_par.x <nproc>"
   exit 1
}

if [ $# -lt 1 ]; then
   echo 'No number of CPUs supplied'
   usage
fi

nproc=$1

cd data
# Delete possible old links
rm * 2> /dev/null
# Take the directory list
ls -d */ | sed -e "s/\///g" > directory_list
cd ..

while IFS=" " read -r dir ; do
   cd data

   rm *.txt 2> /dev/null
   rm ../input/gamma_est_mc.nml 2> /dev/null
   
   ln -s $dir/diam_mc.txt
   ln -s $dir/rho_mc.txt
   ln -s $dir/period_mc.txt
   ln -s $dir/gamma_mc.txt
   ln -s $dir/dadt_mc.txt

   cd ../input 
   
   ln -s ../data/$dir/gamma_est_mc.nml 

   cd ..

   echo " $dir "

   ./gamma_est_mc.x & 2> /dev/null 

   sleep 1

   # Take the number of processes running. If nproc processes
   # are running, check every 10 seconds if some of them has finished.
   # In this manner we almost always keep occupied nproc processes
   currRunning=$(ps -el | grep gamma_est_mc.x | wc -l)
   while [ $currRunning -ge $nproc ]; do
      sleep 10
      currRunning=$(ps -el | grep gamma_est_mc.x | wc -l)
   done

done < data/directory_list 

wait

cd ..
