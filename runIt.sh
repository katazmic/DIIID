#!/bin/bash

module load gsl
#/1.13

make DIIID 
./DIIID 


gnuplot -persist plot.gnu

