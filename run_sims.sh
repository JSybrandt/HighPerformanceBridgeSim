#!/bin/bash

FLAG=$1
DAYS_TO_SIM=${DAYS_TO_SIM:-710}
JOBS_PER_NODE=${JOBS_PER_NODE:-12}
SCRATCH_DIR="/scratch4/$USER/bridge_sim"
mkdir -p $SCRATCH_DIR

# Support conda deactivate
MODULE_REF_FILE="/home/jsybran/module_groups/bridge_proj"
CONDA_PREP_FILE=~/anaconda3/etc/profile.d/conda.sh
source $CONDA_PREP_FILE
source $MODULE_REF_FILE

#Gets the directory containing this file
PROJ_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

PREP_BIN=$PROJ_DIR/bin/Bridge_Pre_Allocation
SIM_BIN=$PROJ_DIR/bin/SimulateDay

WEATHER_FILES=(\
  $PROJ_DIR/data/chicago_linear_interp_10-11.mat \
  $PROJ_DIR/data/USW00003870.hourly.mat \
)

WEATHER_SETTINGS=( 0 1 )
VEHICLE_SETTINGS=( 0 1 )
DAMAGE_SETTINGS=( 1 2 )
BRIDGE_DATA=$PROJ_DIR/data/BridgeVariables.dat
VEHICLE_DATA=$PROJ_DIR/data/VehicleData.dat

# these files describe each parallel task
SCENARIO_DESCRIPTIONS=$(mktemp)
DAILY_SIM_DESCRIPTIONS=$(mktemp)
for VEHICLE_SETTING in ${VEHICLE_SETTINGS[@]}; do
for DAMAGE_SETTING in ${DAMAGE_SETTINGS[@]}; do
for WEATHER_SETTING in ${WEATHER_SETTINGS[@]}; do
for WEATHER_DATA in ${WEATHER_FILES[@]}; do
  WORK_DIR=$(mktemp -dp $SCRATCH_DIR)
  DAY_DIR=$WORK_DIR/daily_sim
  mkdir -p $DAY_DIR
  OUT_REF_FILE=$WORK_DIR/reference_data

  echo $OUT_REF_FILE \
       $VEHICLE_SETTING \
       $DAMAGE_SETTING \
       $WEATHER_SETTING \
       $WEATHER_DATA \
  >> $SCENARIO_DESCRIPTIONS
  echo """
  This file was created on : $(date)
  In this directory we simulate:
    Multiple_Vehicles=$VEHICLE_SETTING
    Damage_Case=$DAMAGE_SETTING
    Environmental_Effects=$WEATHER_SETTING
    num_days=$DAYS_TO_SIM
    weather_data_path=$WEATHER_DATA
    bridge_data_path=$BRIDGE_DATA
    vehicle_data_path=$VEHICLE_DATA
  The results of the scenario prep are stored in $OUT_REF_FILE
  The results of each daily trial are stored in $DAY_DIR
  """ > $WORK_DIR/description.txt

  for DAY in $(seq $DAYS_TO_SIM); do
    echo $OUT_REF_FILE $DAY $DAY_DIR/$DAY >> $DAILY_SIM_DESCRIPTIONS
  done

done done done done

if [[ $FLAG == '--dry-run' ]]; then
  cat $SCENARIO_DESCRIPTIONS
  cat $DAILY_SIM_DESCRIPTIONS
elif [[ $FLAG == '--clean' ]]; then
  rm -rf $SCRATCH_DIR
elif [[ -z $FLAG ]]; then
  echo "Running $PREP_BIN on each scenario"
  # Run prep
  parallel \
    --sshloginfile "$PBS_NODEFILE" \
    --workdir $PROJ_DIR \
    -j $JOBS_PER_NODE \
    --colsep ' ' \
    """
    OUT_REF_FILE={1}
    VEHICLE_SETTING={2}
    DAMAGE_SETTING={3}
    WEATHER_SETTING={4}
    WEATHER_DATA={5}
    source $CONDA_PREP_FILE
    source $MODULE_REF_FILE
    $PREP_BIN \
      \$OUT_REF_FILE \
      \$VEHICLE_SETTING \
      \$DAMAGE_SETTING \
      \$WEATHER_SETTING \
      $DAYS_TO_SIM \
      $WEATHER_DATA \
      $BRIDGE_DATA \
      $VEHICLE_DATA
    echo Created \$OUT_REF_FILE
    """ < $SCENARIO_DESCRIPTIONS

  echo "Simulating each day"
  parallel \
    --sshloginfile "$PBS_NODEFILE" \
    --workdir $PROJ_DIR \
    -j $JOBS_PER_NODE \
    --colsep ' ' \
    """
    REFERENCE_DATA={1}
    DAY={2}
    OUTPUT={3}
    source $CONDA_PREP_FILE
    source $MODULE_REF_FILE
    echo Running $SIM_BIN \
                \$REFERENCE_DATA \
                \$DAY \
                \$OUTPUT
    $SIM_BIN \
      \$REFERENCE_DATA \
      \$DAY \
      \$OUTPUT
    echo Simulated \$OUTPUT
    """ < $DAILY_SIM_DESCRIPTIONS
else
  echo "Did not understand $FLAG"
  echo "Usage: ./run_sims.sh [--dry-run|--clean]"
  exit 1
fi
rm $SCENARIO_DESCRIPTIONS
rm $DAILY_SIM_DESCRIPTIONS
