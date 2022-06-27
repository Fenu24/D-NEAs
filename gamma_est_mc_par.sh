#/bin/bash

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

usedproc=0
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

   let usedproc=usedproc+1
   if [ $usedproc -ge $nproc ]; then
     wait
     usedproc=0
   fi
done < data/directory_list 

wait

cd ..

