#!/bin/bash

HOME=$PWD

PATH=$PATH:$PWD/bin

module load matlab

cd data

Bridge_Pre_Allocation ../results/intermediate_data

