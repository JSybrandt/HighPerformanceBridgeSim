#!/usr/bin/env python3

import os
import scipy.io as sio
import numpy as np
import tensorflow.contrib.learn as skflow
import tensorflow as tf
import math
import os.path
import sys
from util import col, getObCountPerDc, sampleFromDc, getClassifierShape
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

def write_confidence_plot(matrix, out_path):

  num_training_features = col['feature_end'] - col['feature_start']
  dcs, dc_counts = getObCountPerDc(matrix)

  feat_cols = [tf.contrib.layers.real_valued_column("",
               dimension=num_training_features)]

  vprint("Loading / initializing classifier")
  classifier = skflow.DNNClassifier(
      getClassifierShape(),
      feat_cols,
      n_classes=len(dcs),
      model_dir=model_dir
  )

  X = matrix[:, col['feature_start']:col['feature_end']].astype(np.float32)

  days = [x[0] for x in matrix[:, col['day']].tolist()]
  realL = [x[0] for x in matrix[:, col['damage_class']].tolist()]
  predL = [p for p in classifier.predict(X)]

  day2classcounts = {d: [0 for d in dcs] for d in np.unique(days)}
  day2correct = {}

  vprint("Prepared counts for", len(day2classcounts), "days")

  for i in range(X.shape[0]):
      day = days[i]
      real = realL[i]
      pred = predL[i]
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
  for dc in reversed(dcs):
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
                   facecolor=dc2color[dcs[-1]], edgecolor=dc2color[dcs[-1]])
  ax2.text((uniq_days[-1]+last_day)/2, 0.5, "D{}".format(dcs[-1]),
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

