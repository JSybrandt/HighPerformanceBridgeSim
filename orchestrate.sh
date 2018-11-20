#!/bin/bash
#PBS -N runBridge
#PBS -l select=10:ncpus=24:mem=120gb,walltime=72:00:00

# Doop Doop
#  - - PBS -o /dev/null
#  - - PBS -e /dev/null

# Includes
module load matlab
module load gnu-parallel

# Move to location of qsub
cd $PBS_O_WORKDIR

# Constants
SCRATCH="/scratch2/$USER"
PROJECT="bridge"
DAYS_TO_SIMULATE=719
RECORD_BOTH_CARS=0

PROJECT_DIR="$SCRATCH/$PROJECT"
mkdir -p $PROJECT_DIR

PRE_ALLOCATION="bin/Bridge_Pre_Allocation"
if [ ! -f "$PRE_ALLOCATION" ]; then
  echo "Failed to find $PRE_ALLOCATION. Run make"
  exit 1
fi

SIMULATE="bin/SimulateDay"
if [ ! -f "$SIMULATE" ]; then
  echo "Failed to find $SIMULATE. Run make"
  exit 1
fi

COMBINE_DATA="bin/createDataMatrix"
if [ ! -f "$COMBINE_DATA" ]; then
  echo "Failed to find $COMBINE_DATA. Run make"
  exit 1
fi

SCALE_DATA="bin/zero_one_scale"
if [ ! -f "$SCALE_DATA" ]; then
  echo "Failed to find $SCALE_DATA. Run make"
  exit 1
fi

# For each of our bridge cases
for MULT_VEHICLES in 0 1; do
  for DAMAGE_CASE in 1 2; do
    for ENVIRONMENTAL_EFFECTS in 0 1; do

      # Get project name
      if [ $MULT_VEHICLES -eq 0 ]; then
        MV_TAG="one_car"
      else
        MV_TAG="two_cars"
      fi
      if [ $DAMAGE_CASE -eq 1 ]; then
        DC_TAG="continuous_damage"
      else
        DC_TAG="discrete_damage"
      fi
      if [ $ENVIRONMENTAL_EFFECTS -eq 0 ]; then
        EN_TAG="env_off"
      else
        EN_TAG="env_on"
      fi

      CASE_DIR="$PROJECT_DIR/$MV_TAG.$DC_TAG.$EN_TAG"
      mkdir -p $CASE_DIR

      # If the pre_allocation file hasn't been made
      if [ ! -f "$CASE_DIR/pre_allocation.mat" ]; then
        # Run this case
        $PRE_ALLOCATION \
          $CASE_DIR/pre_allocation \
          $MULT_VEHICLES \
          $DAMAGE_CASE \
          $ENVIRONMENTAL_EFFECTS
      fi

      DAY_DIR="$CASE_DIR/daily_sim"
      mkdir -p $DAY_DIR
      # If the number of days simulated is less than the days we would like
      if [ $( ls $DAY_DIR | wc -l ) -lt $DAYS_TO_SIMULATE ]; then
        # For each day (1 to $DAYS_TO_SIMULATE) run SimulateDay
        parallel --workdir $PBS_O_WORKDIR --sshloginfile $PBS_NODEFILE \
          "module load matlab; \
           $SIMULATE $CASE_DIR/pre_allocation {} $DAY_DIR/{}" \
        ::: $(seq 1 $DAYS_TO_SIMULATE)
      fi

      if [ ! -f "$CASE_DIR/combined_raw_data.mat" ]; then
        $COMBINE_DATA \
          $CASE_DIR/dayData \
          $CASE_DIR/combined_raw_data \
          $RECORD_BOTH_CARS
      fi

      # Scale data and prep for machine learning!
      if [ ! -f "$CASE_DIR/training_data.mat" ]; then
        $SCALE_DATA \
          $CASE_DIR/combined_raw_data.mat \
          $CASE_DIR/training_data.mat
      fi

    done # end ENVIRONMENTAL_EFFECTS
  done # end DAMAGE_CASE
done # end MULT_VEHICLES
