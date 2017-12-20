module load matlab
DATA_FILE=results/simData.mat
INTERM_DATA=results/svrData
mkdir -p $INTERM_DATA
TRAINING=$INTERM_DATA/training
VALIDATION=$INTERM_DATA/validation
bin/writeForLibSVM $DATA_FILE $TRAINING $VALIDATION "2 3" 14
