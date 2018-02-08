#!/bin/bash
#PBS -N runBridge
#PBS -l select=1:ncpus=16:mem=60gb,walltime=24:00:00
#PBS -o /dev/null
#PBS -e /dev/null
#PBS -J 0-24


module load matlab
module load gnu-parallel

if [ -z "$PBS_O_WORKDIR" ]; then
  echo "Must be submitted as a PBS job"
  exit 1
fi

cd $PBS_O_WORKDIR
if [ -z "$PBS_ARRAY_INDEX" ]; then
  echo "Must be submitted as array job"
  exit 1
fi

NUM_DAYS=361
JOBS_PER_NODE=$(($NUM_DAYS / 24))
START=$(($JOBS_PER_NODE * $PBS_ARRAY_INDEX))
END=$(($START + $JOBS_PER_NODE - 1))
if [ $PBS_ARRAY_INDEX == 24 ]; then
  END=$NUM_DAYS
fi

if [ $START == 0 ];then
  START=1
fi

HOME=$(pwd)
RES=$HOME/results
DATA=$RES/intermediate_data
OUT_DIR=$RES/dayData

mkdir -p $RES
mkdir -p $OUT_DIR

if [ -f "$DATA.mat" ]; then
  echo "Using existing data file: $DATA"
else
  ./bin/Bridge_Pre_Allocation $DATA
fi

parallel "./bin/SimulateDay $DATA {} $OUT_DIR/{}" ::: $(seq $START $END)
