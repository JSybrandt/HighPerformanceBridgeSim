# Bridge simulation

In order to prep your palmetto environment, I recommend
running the following commands:

``` bash
# Remove other modules that may conflict
module purge
# Add matlab for simulation
module load matlab
# add the most up-to-date version of python
module load anaconda3/5.01
# add nvidia modules for GPU acceleration of the ML model
module load cuda-toolkit/9.0.176
module load cuDNN/9.0v7
```

## Running the simulation

From there, `make -j` will compile the matlab scripts from
`code/analysis` and `code/simulate` into `bin` so that they
may be run.

`orchestrate.sh` is a script that can launch your simulations.
Its setup to run the gamut, and will write results to your
scratch2 directory.

Once the results are finished. You can then run train.py, 
pointing it at the finished `.../training_data.mat` file.
This will handle the ML pipeline, and spits out confidence
numbers as well as the plot.

Finally, with `qsub orchestrate.sh` you will launch the bridge
simulation on the cluster.

## Running the ML Model

In order to prep the machine learning model, you will need
to install tensorflow for palmetto following
[this guide](https://www.palmetto.clemson.edu/palmetto/software_tensorflow.html).
(The script above will have already prepped your modules, 
this guide will setup a python environment.)

Additionally, there are some python dependencies for train.py
beyond tensorflow. These can be installed by running the following
AFTER you have completed the tensorflow setup.

```
pip3 install -r code/train_model/requirements.txt
```

From there on, you can run the training by first loading
modules using the first script, and then secondly running
`source activate tf_env` (or whatever you named your python
environment from the TF tutorial.)
