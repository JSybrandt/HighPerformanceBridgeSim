#!/usr/bin/env python3

import argparse
from pathlib import Path
import logging
import os
import scipy.io as sio
import numpy as np
import sys

# print confidence plot
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

# Needed for classification
import keras
from keras.layers import Dense
from keras.layers import Input
from keras.models import Model
from sklearn.preprocessing import OneHotEncoder
from keras.callbacks import EarlyStopping
from sklearn.metrics import classification_report
from sklearn.model_selection import train_test_split

def config_args():
  parser = argparse.ArgumentParser()
  parser.add_argument("--log-to-stderr", action='store_true')
  parser.add_argument("--matrix", type=str)
  parser.add_argument("--model", type=str)
  parser.add_argument("--confidence-plot", type=str)
  parser.add_argument("--test-ratio", type=float, default=0.2)
  parser.add_argument("--epochs", type=int, default=100)
  return parser.parse_args()

col = {
    'feature_start': 0,
    'peak_freq_1': 0,
    'peak_freq_2': 1,
    'peak_freq_3': 2,
    'peak_freq_4': 3,
    'peak_freq_5': 4,
    'peak_freq_6': 5,
    'peak_freq_7': 6,
    'peak_freq_8': 7,
    'peak_freq_9': 8,
    'peak_freq_10': 9,
    'peak_amp_1': 10,
    'peak_amp_2': 11,
    'peak_amp_3': 12,
    'peak_amp_4': 13,
    'peak_amp_5': 14,
    'peak_amp_6': 15,
    'peak_amp_7': 16,
    'peak_amp_8': 17,
    'peak_amp_9': 18,
    'peak_amp_10': 19,
    'mass': 20,
    'speed': 21,
    'temp': 22,
    'rain': 23,
    'feature_end': 23,
    'day': 24,
    'damage_class': 25
}

def get_count_per_damage_class(matrix):
    dc_obs = matrix[:, col['damage_class']]
    damage_classes, counts = np.unique(np.asarray(dc_obs), return_counts=True)
    damage_classes = [int(d) for d in damage_classes]
    return damage_classes, counts

def sample_from_damage_classes(matrix, dc, count):
    target_rows = np.where(matrix[:, col['damage_class']] == dc)[0]
    sample_rows = np.random.choice(target_rows, count)
    return matrix[sample_rows, :]

def get_classifier(num_training_features):
  input_layer = Input((num_training_features,), name="input", dtype=np.float32)
  hidden_1 = Dense(num_training_features, activation='relu')(input_layer)
  hidden_2 = Dense(num_training_features*2, activation='relu')(hidden_1)
  hidden_3 = Dense(num_training_features*3, activation='relu')(hidden_2)
  hidden_4 = Dense(num_training_features*3, activation='relu')(hidden_3)
  hidden_5 = Dense(num_training_features*2, activation='relu')(hidden_4)
  hidden_6 = Dense(num_training_features, activation='relu')(hidden_5)
  damage_class = Dense(len(damage_classes), activation='softmax', name='damage_class')(hidden_6)
  classifier = Model(inputs=[input_layer], outputs=[damage_class])
  classifier.compile(optimizer='adagrad', loss='mse')
  return classifier

def log_to_stderr():
  logging.basicConfig(level=logging.INFO)

def feature_columns(data):
  return data[:, col['feature_start']:col['feature_end']].astype(np.float32)

def label_columns(data, one_hot_encoder, fit=False):
  labels = data[:, col['damage_class']].astype(int)
  labels = labels.reshape(len(labels), -1)
  if fit:
    return one_hot_encoder.fit_transform(labels)
  else:
    return one_hot_encoder.transform(labels)

def print_confidence_plot(plot_data, one_hot_encoder):
  

if __name__ == "__main__":
  args = config_args()
  if args.log_to_stderr:
    log_to_stderr()

  logging.info("Checking input path")
  input_matrix_path = Path(args.matrix)
  assert input_matrix_path.is_file()

  logging.info("Checking output path")
  output_model_path = Path(args.model)
  assert output_model_path.parent.is_dir()
  assert not output_model_path.exists()

  if args.confidence_plot is not None:
    logging.info("Checking that we can write confidence plot")
    confidence_plot_path = Path(args.confidence_plot)
    assert confidence_plot_path.parent.is_dir()
    assert not confidence_plot_path.exists()
    assert confidence_plot_path.suffix == ".png"

  logging.info("Checking for valid test-ratio")
  assert args.test_ratio > 0
  assert args.test_ratio < 1

  logging.info("Checking for valid # epochs")
  assert args.epochs > 0

  logging.info("Loading %s", input_matrix_path)
  input_data = sio.loadmat(input_matrix_path)
  training_data, testing_data = train_test_split(input_data,
                                                 test_size=args.test_ratio)

  num_observations = training_data.shape[0]

  # Balance damage classes for better performance
  damage_classes, dc_counts = get_count_per_damage_class(training_data)
  damaged_obs = sum(dc_counts[1:])
  logging.info("Counts per dc: %s", str(dc_counts))

  dc2target_count = {}
  dc2target_count[0] = int(num_observations*0.3)
  intended_per_dam = 1/len(damage_classes[1:])
  for d in damage_classes[1:]:
    dc_per_dam = (dc_counts[d] / damaged_obs)
    total_dam = num_observations * 0.7
    if dc_per_dam < intended_per_dam:
      # oversample underrepresented classes
      dc_per_dam = intended_per_dam
    dc2target_count[d] = int(dc_per_dam * total_dam)

  balanced_class_mat = np.concatenate([
    sample_from_damage_classes(training_data, dc,i dc2target_count[dc])
    for dc in damage_classes])

  damage_classes, dc_counts = get_count_per_damage_class(balanced_class_mat)
  logging.info("New Counts per dc: %s", str(dc_counts))

  logging.info("Shuffling training_data")
  shuffled_matrix = np.random.permutation(balanced_class_mat)

  logging.info("Setting up training data")
  # Turns (0-5) damage class into [0, 0, 1, 0, 0, 0]
  one_hot_encoder = OneHotEncoder(len(damage_classes))

  X_train = feature_columns(shuffled_matrix)
  Y_train = label_columns(shuffled_matrix, one_hot_encoder, fit=True)
  num_training_examples, num_training_features = X_train.shape

  #Model
  logging.info("Training classifier")
  classifier = get_classifier(num_training_features)
  stopper = EarlyStopping(monitor="loss")
  classifier.fit(X_train, Y_train, epochs=args.epochs, batch_size=500, callbacks=[stopper])

  X_test = feature_columns(testing_data)
  Y_test = label_columns(testing_data, one_hot_encoder)



def write_confidence_plot(matrix, classifier, out_path):

  damage_classes, dc_counts = get_count_per_damage_class(matrix)


  X = feature_columns(matrix)

  days = matrix[:, col['day']].tolist()
  # Note that here, damage class is an int
  ground_truth_damage_class = matrix[:, col['damage_class']].tolist()
  # and here, damage class is a probability distribution
  predicted_labels = [np.argmax(per_dc_conf)
                      for per_dc_conf
                      in classifier.predict(X)]

  day2classcounts = {d: [0 for d in damage_classes] for d in np.unique(days)}
  day2correct = {}

  vprint("Prepared counts for", len(day2classcounts), "days")

  for i in range(X.shape[0]):
      day = days[i]
      real = ground_truth_damage_class[i]
      pred = predicted_labels[i]
      day2classcounts[day][pred] += 1
      if day not in day2correct:
          day2correct[day] = int(real)

  vprint("Starting plot")
  fig, (ax1, ax3, ax2) = plt.subplots(3, sharex=True,
                                      figsize=(8, 4),
                                      gridspec_kw={
                                          'height_ratios': [5, 1, 1]})

  vprint("Plotting days:")
  uniq_days = np.unique(days)

  dc2color = {
          0: "darkgreen",
          1: "gold",
          2: "darkorange",
          3: "sienna",
          4: "coral",
          5: "firebrick"
  }
  grey_color = (1, 1, 1, 0.6)
  for dc in reversed(damage_classes):
      fill_data = [
          sum([x for i, x in enumerate(day2classcounts[day]) if i <= dc])
          / sum(day2classcounts[day])
          for day in uniq_days
      ]
      is_max = [
          day2classcounts[day][dc] == max(day2classcounts[day])
          for day in uniq_days
      ]

      ax1.fill_between(uniq_days, 0, fill_data, label="dc: {}".format(dc),
                       facecolor=dc2color[dc], edgecolor=dc2color[dc])
      ax3.fill_between(uniq_days, 0, 1, where=is_max,
                       facecolor=dc2color[dc], edgecolor=dc2color[dc])

  vprint("Plotting DC separators")
  last_dc = 0
  last_day = 0
  for day in uniq_days:
      dc = day2correct[day]
      if dc != last_dc:
          ax2.fill_between([last_day, day], 0, 1,
                           facecolor=dc2color[last_dc],
                           edgecolor=dc2color[last_dc])
          ax2.text((day+last_day)/2, 0.5, "D{}".format(last_dc),
                   horizontalalignment="center",
                   verticalalignment="center")
          last_dc = dc
          last_day = day
          ax1.axvline(x=day, color='black')
          ax2.axvline(x=day, color='black')
          ax3.axvline(x=day, color='black')
  ax2.fill_between([last_day, uniq_days[-1]], 0, 1,
                   facecolor=dc2color[damage_classes[-1]], edgecolor=dc2color[damage_classes[-1]])
  ax2.text((uniq_days[-1]+last_day)/2, 0.5, "D{}".format(damage_classes[-1]),
           horizontalalignment="center",
           verticalalignment="center")

  ax1.set_ylim(0, 1)
  ax1.set_xlim(0, 720)
  ax2.set_xlim(0, 720)
  ax3.set_xlim(0, 720)
  ax2.yaxis.set_visible(False)
  ax3.yaxis.set_visible(False)

  ax1.set_ylabel("Percentage of Cars")
  ax1.set_title("Car Predictions Per Day")

  ax2.set_xlabel("Time (Days)")

  ax2.text(0, 0.5, "Truth", fontsize=10,
           horizontalalignment='right',
           verticalalignment='center',
           rotation="vertical", transform=ax2.transAxes)
  ax3.text(0, 0.5, "Max", fontsize=10,
           horizontalalignment='right',
           verticalalignment='center',
           rotation="vertical", transform=ax3.transAxes)

  fig.tight_layout()

  plt.savefig(out_path, format="png", bbox_inches="tight")
  plt.close()

