#!/usr/bin/env python3

"""
This file is responsible for training a machine learning model provided the
input from the matlab simulation.
Input from the --matrix argument is parsed and used to train out model.
"""

import argparse
from pathlib import Path
import logging
import os
import scipy.io as sio
import numpy as np
import sys
from copy import copy

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


def config_args():
  parser = argparse.ArgumentParser()
  parser.add_argument("--log-to-stderr", action='store_true')
  parser.add_argument("--matrix", type=str)
  parser.add_argument("--model", type=str)
  parser.add_argument("--confidence-plot", type=str)
  parser.add_argument("--test-ratio", type=float, default=0.2)
  parser.add_argument("--epochs", type=int, default=100)
  return parser.parse_args()


def get_count_per_damage_class(matrix):
  dc_obs = matrix[:, col['damage_class']]
  damage_classes, counts = np.unique(np.asarray(dc_obs), return_counts=True)
  damage_classes = [int(d) for d in damage_classes]
  return damage_classes, counts


def sample_from_damage_classes(matrix, damage_class, count):
  target_rows = np.where(matrix[:, col['damage_class']] == damage_class)[0]
  sample_rows = np.random.choice(target_rows, count)
  return matrix[sample_rows, :]


def get_classifier(num_training_features):
  input_layer = Input((num_training_features,), name="input", dtype=np.float32)
  hidden_1 = Dense(num_training_features, activation='relu')(input_layer)
  hidden_2 = Dense(num_training_features * 2, activation='relu')(hidden_1)
  hidden_3 = Dense(num_training_features * 3, activation='relu')(hidden_2)
  hidden_4 = Dense(num_training_features * 3, activation='relu')(hidden_3)
  hidden_5 = Dense(num_training_features * 2, activation='relu')(hidden_4)
  hidden_6 = Dense(num_training_features, activation='relu')(hidden_5)
  damage_class = Dense(
      len(damage_classes), activation='softmax', name='damage_class')(
          hidden_6)
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

def real_and_est_damage_classes_per_day(matrix, classifer, damage_classes):
  X = feature_columns(matrix)
  days = matrix[:, col['day']].tolist()
  actual_damage_classes = matrix[:, col['damage_class']].tolist()
  predicted_damage_classes = [np.argmax(p) for p in classifier.predict(X)]
  day2prediction_counts = {d: [0 for d in damage_classes] for d in np.unique(days)}
  day2actual = {}

  logging.info("Prepared counts for %s days", len(day2prediction_counts))
  for i in range(X.shape[0]):
    day = days[i]
    actual = actual_damage_classes[i]
    pred = predicted_damage_classes[i]
    day2prediction_counts[day][pred] += 1
    if day not in day2actual:
      day2actual[day] = int(actual)
  return day2actual, day2prediction_counts

def write_confidence_plot(day2correct, day2classcounts, out_path):
  days = list(day2correct.keys())

  fig, (ax2, ax1) = plt.subplots(2, sharex=True,
                               figsize=(8, 2.5),
                               gridspec_kw={'height_ratios': [1, 4]})

  logging.info("Plotting days:")
  uniq_days = np.unique(days)
  damage2color = {
      0: "darkgreen",
      1: "gold",
      2: "darkorange",
      3: "sienna",
      4: "coral",
      5: "firebrick"
  }

  largest_incorrect_each_day = [
      max([
          x for i, x in enumerate(day2classcounts[day]) if i != day2correct[day]
      ]) / sum(day2classcounts[day]) for day in uniq_days
  ]
  correct_amt_each_day = [
      day2classcounts[day][day2correct[day]] / sum(day2classcounts[day])
      for day in uniq_days
  ]

  confidence_data = [
      i - j for i, j in zip(correct_amt_each_day, largest_incorrect_each_day)
  ]

  # reversed to improve coloring
  for damage_class in reversed(damage_classes):
    days_when_largest = [
        day2classcounts[day][damage_class] == max(day2classcounts[day])
        for day in uniq_days
    ]

    ax1.fill_between(
        uniq_days,
        0,
        confidence_data,
        where=days_when_largest,
        facecolor=damage2color[damage_class],
        edgecolor=damage2color[damage_class])

  logging.info("Plotting DC separators")
  last_dc = 0
  last_day = 0
  for day in uniq_days:
    damage_class = day2correct[day]
    if damage_class != last_dc:
      ax2.fill_between([last_day, day],
                       0,
                       1,
                       facecolor=damage2color[last_dc],
                       edgecolor=damage2color[last_dc])
      ax2.text(
          (day + last_day) / 2,
          0.5,
          "D{}".format(last_dc),
          horizontalalignment="center",
          verticalalignment="center")
      last_dc = damage_class
      last_day = day
      ax1.axvline(x=day, color='black')
      ax2.axvline(x=day, color='black')

  ax2.fill_between([last_day, uniq_days[-1]],
                   0,
                   1,
                   facecolor=damage2color[damage_classes[-1]],
                   edgecolor=damage2color[damage_classes[-1]])

  ax2.text(
      (uniq_days[-1] + last_day) / 2,
      0.5,
      "D{}".format(damage_classes[-1]),
      horizontalalignment="center",
      verticalalignment="center")

  ax1.set_ylim(-1, 1)
  ax1.set_xlim(0, 720)
  ax2.set_xlim(0, 720)
  ax2.yaxis.set_visible(False)

  ax1.set_ylabel("Confidence")
  ax1.set_xlabel("Time (Days)")

  ax2.text(
      0,
      0.5,
      "Truth",
      fontsize=10,
      horizontalalignment='right',
      verticalalignment='center',
      rotation="vertical",
      transform=ax2.transAxes)

  ax2.set_title("Car Predictions Per Day")
  fig.tight_layout()

  plt.savefig(out_path, format="png", bbox_inches="tight")
  plt.close()

def summarize_classifer_performance(day2actual, day2counts):
  num_exactly_correct = 0 # True if exactly correct
  num_boolean_correct = 0 # Healthy vs unhealthy
  for day, actual in day2actual.items():
    counts = day2counts[day]
    actual_count = counts[actual]
    tmp_counts = copy(counts)
    tmp_counts[actual] = 0
    largest_incorrect_class = np.argmax(tmp_counts)
    largest_incorrect_class_count = tmp_counts[largest_incorrect_class]
    # If we are exactly correct
    if actual_count > largest_incorrect_class_count:
      num_exactly_correct += 1
      num_boolean_correct += 1
    # if we are instead, close to correct (both healthy or both not)
    elif (actual == 0) == (largest_incorrect_class == 0):
      num_boolean_correct += 1
  return (num_exactly_correct / len(day2actual),
          num_boolean_correct / len(day2actual))

if __name__ == "__main__":
  args = config_args()
  if args.log_to_stderr:
    log_to_stderr()

  ##################
  # Checking input #
  ##################

  logging.info("Checking that ANY output flag is set")
  logging.info("Otherwise, whats the point?")
  assert args.confidence_plot is not None or args.model is not None

  logging.info("Checking input path")
  input_matrix_path = Path(args.matrix)
  assert input_matrix_path.is_file()

  if args.model is not None:
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

  #####################################
  # Parsing matrix, balancing classes #
  #####################################

  # Load
  logging.info("Loading %s", input_matrix_path)
  input_data = sio.loadmat(str(input_matrix_path))['data']

  # Split
  logging.info("Train-Test Split")
  training_data, testing_data = train_test_split(
      input_data, test_size=args.test_ratio)

  num_observations = training_data.shape[0]
  intended_total_unhealthy_examples = num_observations * 0.7

  # Balance damage classes for better performance
  damage_classes, damage_class2count = get_count_per_damage_class(training_data)
  num_unhealthy_examples = sum(damage_class2count[1:])
  logging.info("Counts per damage_class: %s", str(damage_class2count))

  damage_class2intended_count = {}
  damage_class2intended_count[0] = int(num_observations * 0.3)
  # we need to make some intermediate damage classes healthier
  min_unhealthy_ratio = 1 / len(damage_classes[1:])
  for damage_class in damage_classes[1:]:
    current_ratio = (damage_class2count[damage_class] / num_unhealthy_examples)
    if current_ratio < min_unhealthy_ratio:
      # oversample underrepresented classes
      current_ratio = min_unhealthy_ratio
    damage_class2intended_count[damage_class] = \
        int(current_ratio * intended_total_unhealthy_examples)

  balanced_class_mat = np.concatenate([
      sample_from_damage_classes(training_data, damage_class,
                                 damage_class2intended_count[damage_class])
      for damage_class in damage_classes
  ])

  damage_classes, damage_class2count = get_count_per_damage_class(balanced_class_mat)
  logging.info("New Counts per damage_class: %s", str(damage_class2count))

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
  classifier.fit(
      X_train, Y_train, epochs=args.epochs, batch_size=500, callbacks=[stopper], verbose=(1 if args.log_to_stderr else 0))

  train_day2actual, train_day2precited_counts = real_and_est_damage_classes_per_day(training_data, classifier, damage_classes)
  test_day2actual, test_day2precited_counts = real_and_est_damage_classes_per_day(testing_data, classifier, damage_classes)

  logging.info("Printing training data performance summary")
  exact, boolean = summarize_classifer_performance(train_day2actual, train_day2precited_counts)
  print("Training Set Exact Confidence:  ", exact)
  print("Training Set Boolean Confidence:", boolean)
  logging.info("Printing testing data performance summary")
  exact, boolean = summarize_classifer_performance(test_day2actual, test_day2precited_counts)
  print("Test Set Exact Confidence:  ", exact)
  print("Test Set Boolean Confidence:", boolean)

  if args.model is not None:
    logging.info("Saving model")
    classifier.save(str(output_model_path))

  if args.confidence_plot is not None:
    logging.info("Creating confidence plot")
    write_confidence_plot(test_day2actual, test_day2precited_counts,
                          confidence_plot_path)
