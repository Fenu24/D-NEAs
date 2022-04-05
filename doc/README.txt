Author:      Marco Fenucci
Institution: University of Belgrade
Email:       marco_fenucci@matf.bg.ac.rs
Date:        September 2020
Version:     1.0

 C O N T E N T S
 ===============

 (1) Introduction

 (2) Input parameters and data

 (3) Output

----------------------------------------------------

 (1) I N T R O D U C T I O N 
 ===========================

 This software is used to produce the results contained in the paper:
   Fenucci, Novakovic, Vokrouhlicky: The low thermal conductivity of the fast rotator (499998) 2011 PT

 The main idea is that the thermal conductivity of an asteroid can be deduced
 from the comparison of the measured Yarkovsky effect and the estimated one. 
 The estimation of the Yarkovsky effect used in this software is taken from
 Vokrouhlicky 1999. The parameters needed for the inversion are
 - the semimajor axis
 - the diameter 
 - the density 
 - the heat capacity
 - the rotation period
 - the obliquity
 - the absorption coefficient
 - the emissivity
 Some of the parameters are modeled according to some distributions, while other ones
 are simply taken as constant, see the reference paper for details.

 The software is structured as follows:
   -    bin: contains the binary files of the executables
   -   data: contains the distributions of the input parameters
   -  input: contains the input files with constant parameters
   - matlab: contains scripts and functions to produce input and output
   -    mod: contains the .mod files needed at compilation time
   -    obj: contains the .o files needed at compilation time
   - output: contains the output files
   -    src: contains the source codes
   -   test: contains executables for tests of the code
 
 To build the program simply type
   make
 and to unbuild
   make clean
 To clean the executables, the data folder, and the output folder:
   make veryclean
 The executables produced are:
   yarko_est_mc.x:       the main program for the estimation of the thermal conductivity of a
                         NEO.
   eccentric_yarko_mc.x: compute the difference in the Yarkovsky effect produced by an
                         eccentric orbit.
   test/xu2020_test.x:   test for the implementation of the Yarkovsky drift of Xu et al.
                         2020.
   test/yarko_compare.x: test program to estimate the difference in the computation of the
                         Yarkovsky drift, using the implementations of Vokrouhlicky 1999
                         and Xu et al. 2020.

 (2) I N P U T   P A R A M E T E R S   A N D   D A T A
 =====================================================
 
 ----------------------------
 | Program | yarko_est_mc.x |
 ----------------------------

 For the main program yarko_est_mc.x there are two kinds of input: one for the fixed 
 parameters and another one for the parameters with a distribution.

 2.1 FIXED PARAMETERS
 
 The input file for the fixed parameters is
   input/yarko_est_mc.nml
 Here you have to specify: the heat capacity, the bounds for K in which you want to 
 search for solutions, the semimajor axis of the object, the eccentricity of the orbit,
 the absorption coefficient, the emissivity, the method used to compute the Yarkovskky
 effect, the name of the output file, and a flag to indicate if you want to compare the
 output with an eccentric model.
 The output file will be placed in the folder
   output

 2.2 PARAMETERS WITH A DISTRIBUTION

 To generate the distributions of 
   - the diameter
   - the density
   - the obliquity
   - the measured Yarkovsky drift
   - the rotation period
 you have to use the Matlab script 
   output/gen_distrib.m
 Here you have to provide: 
   - the absolute magnitude and its uncertainty
   - the measured Yarkovsky drift and its uncertainty
   - the rotation period and its uncertainty
   - the probabilities for the object to come from main belt
     source region, as provided by the NEO model of Granvik et. al. 2018
 For the last point, the vector with the probabilities have to be ordered 
 as follows:
   p(1) : nu6 secular resonance
   p(2) : 3:1 JMMR
   p(3) : 5:2 JMMR
   p(4) : Hungaria region
   p(5) : Phocaea region
   p(6) : 2:1 JMMR
   p(7) : Jupiter Family Comets
 
 Once you set these parameters, you can run the script, and it will save the distributions
 in different files, placed in the folder
   data

 ----------------------------------
 | Program | eccentric_yarko_mc.x |
 ----------------------------------

 This program computes the difference of the estimated Yarkovsky drift using two different 
 models: a circular one and onther that takes into account the eccentricity. Also here
 there are fixed parameters and parameters with a distribution. The fixed parameters can
 be setted in the file
   input/eccentric_yarko_eval.txt
 
 The parameters with a distribution are still in the directory
   data
 For the density, the diameter, the obliquity, and the rotation period, they can be
 produced as before, using the matlab script matlab/gen_distrib.m
 The sample for the thermal conductivity K is produced using the output distriution
 obtained with the program yarko_est_mc.x, and it is produced inside the matlab script
   output/plotHist_mc.m
 used to analyze the output.


 (3) O U T P U T 
 ================

 ---------------------------
 | Program: yarko_est_mc.x |
 ---------------------------

 The output file will be placed in the folder
   output
 and it consists of all the values of the thermal conductivity K which satisfy the 
 measured Yarkovsky drift. To plot the histogram, and to compute the probabilities
   P ( K < 0.001 )
   P ( 0.001 < K < 0.1 )
 as well as the local maxima and the corresponding values of the probability
 density function, you can use the Matlab script
   output/plotHist_mc.m
 
 ----------------------------------
 | Program | eccentric_yarko_mc.x |
 ----------------------------------

 The output of this program is the file
   output/diffEst.txt
 The results can be analyzed using the matlab script
   output/plotPercentageError.m

