# Manual for the D-NEAs software

## Table of contents

* [Introduction](#introduction)
* [Compilation](#compilation)
* [Thermal inertia estimation](#thermal-inertia-estimation)
* [License and authors](#license-and-authors)
* [References](#references)

## Introduction

## Compilation

Before using the package for the first time, the code needs to be compiled. To facilitate the user, the distribution
comes with a configuration script and a Makefile that automatically does the job. To compile the source code, 
please follow these steps
1. Choose the compiler and the compilation options by running the config.sh script. By running the script without further options, you will receive an help message. The script permits to choose between two different compiler: GNU gfortran, and Intel ifort. An additional option defines the compilation flags, and the final user can select the optimization flags "-O". For instance, if you want to use the GNU gfortran compiler, you can run the script as
           
           ./config.sh -O gfortran
            
2. Run the Makefile to compile the source code. To this purpose, you can type

            make

   and the executable binary files will be placed in the *bin* directory. Symbolic links will be created in the main directory.


**NOTE 1.** For the compilation to work correctly, you need two hidden directories called *.mod* and *.obj*. Please be sure that these two directories are correctly contained on your local distribution.


## Thermal inertia estimation 

### Generating parameters distributions

### Running the thermal inertia estimation

## License and authors

This software is released under the XXX license, and it developed and maintained by

- [Marco Fenucci](http://adams.dm.unipi.it/~fenucci/index.html), Department of Astronomy, Faculty of Mathematics, University of Belgrade (<marco_fenucci@matf.bg.ac.rs>) 
- [Bojan Novaković](http://poincare.matf.bg.ac.rs/~bojan/index_e.html), Department of Astronomy, Faculty of Mathematics, University of Belgrade (<bojan@matf.bg.ac.rs>) 
- [Dušan Marčeta](http://poincare.matf.bg.ac.rs/~dmarceta/), Department of Astronomy, Faculty of Mathematics, University of Belgrade (<dmarceta@matf.bg.ac.rs>) 

## References
