#!/bin/bash
protoc --python_out ../convert_mat_to_proto result.proto
protoc --python_out ../train_model result.proto
