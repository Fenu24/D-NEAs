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


function usage () {
   RED="\033[0;31m"
   NC="\033[0m"
   echo " "
   echo -e "${RED}SYNTAX${NC}"
   echo "   gamma_est_mc_batches.x -nproc <nproc> -list <list> [-C|Kmin|Kmax|expo|alpha|epsi|method|max_iter]"
   echo " "
   echo -e "${RED}ARGUMENTS${NC}"
   echo "   MANDATORY ARGUMENTS: "
   echo "      nproc      Number of processors to use"
   echo "      list       Name of the list of asteroids to process"
   echo " "
   echo "   OPTIONAL ARGUMENTS:  "
   echo "      C          Value of the heat capacity                         [DEFAULT = 800.d0]"
   echo "      Kmin       Minimum conductivity solution                      [DEFAULT = 1.d-6]"
   echo "      Kmax       Maximum conductivity soluion                       [DEFAULT = 100.d0]"
   echo "      expo       Exponent for conductivity variation                [DEFAULT = 0.d0]"
   echo "      alpha      Absorprion coefficien                              [DEFAULT = 1.d0]"
   echo "      epsi       Emissivity                                         [DEFAULT = 1.d0]"
   echo "      method     Flag for the Yarkovsky model                       [DEFAULT = 2]"
   echo "      max_iter   Max. num. of iterations for the Monte Carlo method [DEFAULT = 1000000]"
   exit 1
}

function help () {
   RED="\033[0;31m"
   NC="\033[0m"
   echo -e "${RED}DESCRIPTION${NC}"
   echo "   This script runs the thermal inertia estimation of a list of objects in parall,"
   echo "   each of them on a single core. "
   echo ""
   echo "   The list of objects is chosen by the user, and on each line it must contain:"
   echo "      1) the name of the object "
   echo "      3) the semi-major axis "
   echo "      3) the eccentricity "
   echo ""
   echo "   The input for the object obj must be places in the folder data/obj. "
   echo ""
   echo "   The number of objects to run at the same time is chosen by the user."
   usage
   exit 1
}


# Read the input flags
while test $# -gt 0; do
   case "$1" in 
        -h)
          help
          ;;
      -nproc)
          shift
          nproc=$1
          shift
          ;;
       -list)
          shift
          list=$1
          shift
          ;;
       -C)
          shift
          C=$1
          shift
          ;;
       -Kmin)
          shift
          Kmin=$1
          shift
          ;;
       -Kmax)
          shift
          Kmax=$1
          shift
          ;;
      -method)
         shift
         method=$1
         shift
         ;;
      -max_iter)
         shift
         max_iter=$1
         shift
         ;;
      -expo)
         shift
         expo=$1
         shift
         ;;
      -epsi)
         shift
         epsi=$1
         shift
         ;;
      -alpha)
         shift
         alpha=$1
         shift
         ;;
       *)
         echo "ERROR: $1 is not a valid flag."
         usage
         ;;
   esac
done

# Check if nproc and list are both provided
if [ -z "$nproc" ]; then
   echo "ERROR: nproc must be provided."
   usage
fi

if [ -z "$list" ]; then
   echo "ERROR: list must be provided."
   usage
fi

# Set the default values if the other options are not provided
if [ -z "$C" ]; then
   C=800.d0 
fi
if [ -z "$Kmin" ]; then
   Kmin=1.d-6
fi
if [ -z "$Kmax" ]; then
   Kmax=100.d0
fi
if [ -z "$method" ]; then
   method=2
fi
if [ -z "$max_iter" ]; then
   max_iter=1000000
fi
if [ -z "$expo" ]; then
   expo=0.d0
fi
if [ -z "$epsi" ]; then
   epsi=1.d0
fi
if [ -z "$aplha" ]; then
   alpha=1.d0
fi

# Print the options on the screen
echo "-------------------------------    "
echo "--------- BATCHES RUN ---------    "
echo "-------------------------------    "
echo "                                   "
echo " Options of the run:               "
echo "   Number of processors: $nproc    "
echo "      List of asteroids: $list     "
echo "                                   "
echo " Parameters for main program:      "
echo "          Heat Capacity: $C        "
echo "              K minimum: $Kmin     "
echo "              K maximum: $Kmax     "
echo "        K var. exponent: $expo     "
echo "      Absorprion coeff.: $alpha    "
echo "             Emissivity: $epsi     "
echo "          Yarko. method: $method   "
echo "        Max. iterations: $max_iter "

exit 1

cd data
# Delete possible old links
#rm * 2> /dev/null
# Take the directory list
ls -d */ | sed -e "s/\///g" > directory_list
cd ..

while IFS=" " read -r dir ; do
   cd data

   #rm *.txt 2> /dev/null
   #rm ../input/gamma_est_mc.nml 2> /dev/null
   
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


