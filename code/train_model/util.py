# This function holds shared functions needed for both train and evaluate
import numpy as np

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


def getClassifierShape():
    num_training_features = col['feature_end'] - col['feature_start']
    return [num_training_features,
            2*num_training_features,
            3*num_training_features,
            3*num_training_features,
            2*num_training_features,
            num_training_features]  # CNN Layers


def getObCountPerDc(matrix):
    dc_obs = matrix[:, col['damage_class']]
    vals, counts = np.unique(np.asarray(dc_obs), return_counts=True)
    vals = [int(d) for d in vals]
    return vals, counts


def sampleFromDc(matrix, dc, count):
    target_rows = np.where(matrix[:, col['damage_class']] == dc)[0]
    sample_rows = np.random.choice(target_rows, count)
    return matrix[sample_rows, :]
