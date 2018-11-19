#!/bin/bash
#PBS -N runBridge
#PBS -l select=10:ncpus=24:mem=120gb,walltime=72:00:00
#PBS -o /dev/null
#PBS -e /dev/null


module load matlab gnu-parallel

if [ -z "$PBS_O_WORKDIR" ]; then
  echo "Not Moving"
else
  cd $PBS_O_WORKDIR
  echo $PBS_O_WORKDIR
fi


TOTAL_DAYS=719

function setup(){
  FILE="code/simulate/Bridge_Pre_Allocation.m"
  MULT=$1
  DAM=$2
  ENV_EFF=$3
  sed -i -E "s/^Multiple_Vehicles=[01]/Multiple_Vehicles=$MULT/g" $FILE
  sed -i -E "s/^Damage_Case=[12]/Damage_Case=$DAM/g" $FILE
  sed -i -E "s/^Temp=[01]/Temp=$ENV_EFF/g" $FILE
  sed -i -E "s/^RainEffects=[01]/RainEffects=0/g" $FILE
  sed -i -E "s/^Surface=[01]/Surface=$ENV_EFF/g" $FILE

  make bin/Bridge_Pre_Allocation
}

function makeResDir(){
  D="/scratch2/jsybran/old_bridge_stuff/rebuilt_bridge/M$1_D$2_E$3"
  mkdir -p $D
  mkdir -p $D/dayData
  echo $D
}

for MULT in 0 1; do
  for DAM in 1 2; do
    for ENV_EFF in 0 1; do

      RES_DIR=$(makeResDir $MULT $DAM $ENV_EFF)
      echo $RES_DIR
      if [ ! -f $RES_DIR/intermediate_data.mat ]; then
        # compile preallocation and setup res dir
        setup $MULT $DAM $ENV_EFF
        # actually preallocate
        ./bin/Bridge_Pre_Allocation $RES_DIR/intermediate_data
      fi

      if [ ! -f $RES_DIR/done.flag ]; then
        echo "WRITING DAYS FOR $RES_DIR"
        parallel --workdir $PBS_O_WORKDIR --sshloginfile $PBS_NODEFILE \
          "module load matlab; \
           ./bin/SimulateDay $RES_DIR/intermediate_data {} $RES_DIR/dayData/{}" \
        ::: $(seq 1 $TOTAL_DAYS)
        touch $RES_DIR/done.flag
      fi

      #COMPILING

      if [ ! -f $RES_DIR/no_other.mat ]; then
        echo "COMPILING $RES_DIR WITH ONLY PRIMARY VEHICLES"
        ./bin/createDataMatrix $RES_DIR/dayData $RES_DIR/no_other 0
      fi

      if [ $MULT -eq 1 ] && [ ! -f $RES_DIR/with_other.mat ]; then
        echo "COMPILING $RES_DIR WITH ALL VEHICLES"
        ./bin/createDataMatrix $RES_DIR/dayData $RES_DIR/with_other.mat 1
      fi

      # SCALING

      if [ ! -f $RES_DIR/no_other.scaled.mat ]; then
        echo "Scaling $RES_DIR/no_other.mat"
        ./bin/zero_one_scale $RES_DIR/no_other.mat $RES_DIR/no_other.scaled.mat
      fi

      if [ $MULT -eq 1 ] && [ ! -f $RES_DIR/with_other.scaled.mat ]; then
        echo "Scaling $RES_DIR/with_other.mat"
        ./bin/zero_one_scale $RES_DIR/with_other.mat $RES_DIR/with_other.scaled.mat
      fi

    done
  done
done

echo "Time to run batch job training"
