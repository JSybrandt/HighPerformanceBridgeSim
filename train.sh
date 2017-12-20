module load matlab
DATA_FILE=results/simData.mat
MODEL_FILE=results/subsetmodel
bin/trainSubset $DATA_FILE $MODEL_FILE "2 3" 14
