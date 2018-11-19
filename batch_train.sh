#!/bin/bash
#PBS -N btrain
#PBS -l select=1:ncpus=24:mem=120gb,walltime=72:00:00
#PBS -o /dev/null
#PBS -e /dev/null
#PBS -J 0-11

cd $PBS_O_WORKDIR
echo $PWD

NUM_MACHINES=12

MODEL_ROOT=/scratch2/jsybran/old_bridge_stuff/rebuilt_bridge
DATA_ROOT_A=/scratch2/jsybran/old_bridge_stuff/rebuilt_bridge
DATA_ROOT_B=/scratch2/jsybran/old_bridge_stuff/original

MY_LIST=`mktemp`
find $DATA_ROOT_A/ -iname "*scaled.mat" | sort > $MY_LIST

TOTAL=`wc -l $MY_LIST | awk '{print $1}'`
echo $TOTAL

Q_PER=$(($TOTAL / $NUM_MACHINES))
START=$(($PBS_ARRAY_INDEX * $Q_PER))
END=$(($START + $Q_PER))
if [ $PBS_ARRAY_INDEX -eq $(($NUM_MACHINES-1)) ];then
  END=$TOTAL
fi

# offset for line numbers
START=$(($START + 1))
END=$(($END + 1))

for ((i = $START ; i < $END ; i++)); do
  MAT_A=`sed -n "$i"p $MY_LIST`
  ROOT_DIR=`dirname $MAT_A`
  EXP=`basename $ROOT_DIR`
  BASE=`basename -s .mat $MAT_A`
  MAT_B=$DATA_ROOT_B/$EXP/$(basename $MAT_A)
	MAT_B="NONE"
  echo $MAT_A $MAT_B

  LOG=$ROOT_DIR/$BASE.log
  rm $LOG

  ./bin/train.py $MAT_A $MAT_B $MODEL_ROOT/$EXP/model_$BASE v # &>> $LOG
  # ./bin/evaluate.py $MAT_A $ROOT_DIR/model_$BASE $ROOT_DIR/$BASE.png v &>> $LOG
  echo "DONE" >> $LOG

done

rm $MY_LIST
