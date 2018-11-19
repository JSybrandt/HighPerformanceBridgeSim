#!/bin/bash
#PBS -N btrain
#PBS -l select=1:ncpus=24:mem=120gb,walltime=72:00:00
#PBS -o /dev/null
#PBS -e /dev/null

MODEL_ROOT=/scratch2/jsybran/old_bridge_stuff/original
DATA_ROOT=/scratch2/jsybran/old_bridge_stuff/rebuilt_bridge
# DATA_ROOT=/scratch2/jsybran/old_bridge_stuff/retrain
# DATA_ROOT=/scratch2/jsybran/old_bridge_stuff/original

for MAT_FILE in `find $DATA_ROOT/ -iname "*scaled.mat"`; do

  echo $MAT_FILE

  #example: /scratch2/jsybran/old_bridge_stuff/rebuilt_bridge/M1_D1_E0
  DATA_EXP_ROOT=`dirname $MAT_FILE`

  #example M1_D1_E0
  EXP=`basename $DATA_EXP_ROOT`

  # gets the name without .mat
  # example: no_other.scaled
  BASE=`basename -s .mat $MAT_FILE`
  MODEL_FILE=$MODEL_ROOT/$EXP/model_$BASE

  IMG_RES=$DATA_ROOT/$EXP/$EXP-$BASE.png

#  ./bin/evaluate.py $MAT_FILE $ROOT_DIR/model_$BASE $ROOT_DIR/$BASE.png &

  ./code/train_model/confidence_plot.py \
    $MAT_FILE $MODEL_FILE $IMG_RES &


done


wait

IMG_DIR=`mktemp -d`

find $DATA_ROOT/ -iname "*.png" -exec mv {} $IMG_DIR/ \;
rm img_results.zip
zip img_results.zip $IMG_DIR/*.png
post img_results.zip
rm -rf $IMG_DIR
