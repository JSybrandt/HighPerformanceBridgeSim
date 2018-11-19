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

global VERBOSE
VERBOSE = False


def vprint(*args, **kargs):
    global VERBOSE
    if VERBOSE:
        print(*args, **kargs)


# PARSE ARGS
if(len(sys.argv) >= 2):
    mat_file = sys.argv[1]
    if not os.path.isfile(mat_file) or not mat_file.endswith(".mat"):
        raise ValueError("Didn't give me a valid mat file")
else:
    raise ValueError("Must supply mat file as arg1")

if(len(sys.argv) >= 3):
    model_dir = sys.argv[2]
    if not os.path.isdir(model_dir):
        raise ValueError("Didn't give me a valid model dir")
else:
    raise ValueError("Must supply model file as arg2")

if(len(sys.argv) >= 4):
    out_path = sys.argv[3]
    if not out_path.endswith(".png"):
        raise ValueError("out file must be a png")
else:
    raise ValueError("Must supply result file as arg1")

if(len(sys.argv) >= 5):
    VERBOSE = True
    tf.logging.set_verbosity(tf.logging.INFO)


vprint("STARTED IN VERBOSE MODE")
vprint("Loading ", mat_file)
mat_data = sio.loadmat(mat_file)
matrix = np.asmatrix(mat_data["data"])
pre_scale_mins = np.asmatrix(mat_data["mins"])
pre_scale_disp = np.asmatrix(mat_data["dists"])

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
fig, (ax2, ax1) = plt.subplots(2, sharex=True,
                               figsize=(8, 2.5),
                               gridspec_kw={'height_ratios': [1, 4]})

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


max_wrong = [
    max([x for i, x in enumerate(day2classcounts[day])
         if i != day2correct[day]])
    / sum(day2classcounts[day])
    for day in uniq_days
]

correct_data = [
    day2classcounts[day][day2correct[day]]
    / sum(day2classcounts[day])
    for day in uniq_days
]

confidence_data = [i-j for i, j in zip(correct_data, max_wrong)]


for dc in reversed(dcs):
    is_max = [
        day2classcounts[day][dc] == max(day2classcounts[day])
        for day in uniq_days
    ]

    ax1.fill_between(uniq_days, 0, confidence_data, where=is_max,
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
ax2.fill_between([last_day, uniq_days[-1]], 0, 1,
                 facecolor=dc2color[dcs[-1]], edgecolor=dc2color[dcs[-1]])
ax2.text((uniq_days[-1]+last_day)/2, 0.5, "D{}".format(dcs[-1]),
         horizontalalignment="center",
         verticalalignment="center")


# ax1.axhline(0.5, linestyle="--", color="lightgray")
# ax1.axhline(0.2, linestyle=":", color="grey")

# ax1.set_ylim(-1, 1)
ax1.set_ylim(-1, 1)
ax1.set_xlim(0, 720)
ax2.set_xlim(0, 720)
ax2.yaxis.set_visible(False)

ax1.set_ylabel("Confidence")

ax1.set_xlabel("Time (Days)")
# ax2.set_title("Correct Damage Classes", fontsize=12)
# ax3.set_title("Majority Vote", fontsize=12)

ax2.text(0, 0.5, "Truth", fontsize=10,
         horizontalalignment='right',
         verticalalignment='center',
         rotation="vertical", transform=ax2.transAxes)

ax2.set_title("Car Predictions Per Day")
fig.tight_layout()

# ax1.legend()
plt.savefig(out_path, format="png", bbox_inches="tight")
plt.close()

