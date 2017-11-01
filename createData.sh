#!/bin/bash
#PBS -N createData
#PBS -l select=1:ncpus=16:mem=60gb,walltime=72:00:00
#PBS -o /dev/null
#PBS -e /dev/null

cd $PBS_O_WORKDIR

HOME=$PWD

PATH=$PATH:$PWD/bin

module load matlab

createDataMatrix results/dayData results/simData

