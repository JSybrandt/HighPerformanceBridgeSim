from result_pb2 import Trial
import stream  # proto stream
from tensorflow.keras.utils import Sequence
from glob import glob
from multiprocessing import Pool
from random import shuffle
import numpy as np

class BaseTrialSequence(Sequence):
  def __init__(self,
               proto_dir_path,
               batch_size,
               fft_type="original",
               value_type="both"):
    """
    proto_dir_path : Location containing '<day#>.proto_stream' files
    fft_type : Either "original" or "filtered"
    value_type : Either "upper", "lower", or "both"
    """
    proto_dir_path = str(proto_dir_path)
    self.stream_paths = glob(proto_dir_path + "/*.proto_stream")
    if len(self.stream_paths) == 0:
      raise ValueError(proto_dir_path + " does not contain any *.proto_stream")

    self.fft_type = fft_type
    assert self.fft_type in ["original", "filtered"]

    self.value_type = value_type
    assert self.value_type in ["upper", "lower", "both"]

    self.batch_size = batch_size
    assert self.batch_size > 0

    # collection of (x, y) pairs
    self.data_labels = []
    with Pool() as pool:
      for res in pool.map(self.load_stream, self.stream_paths):
        self.data_labels += res

    self.on_epoch_end()

  def load_stream(self, stream_path):
    res = []
    with stream.open(stream_path, 'rb') as proto_stream:
      for data in proto_stream:
        trial = Trial()
        trial.ParseFromString(data)
        assert trial.HasField("damage_class")
        res.append((self.to_vec(trial), trial.damage_class))
    return res

  def get_select_fft(self, trial):
    if self.fft_type == "original":
      fft = trial.fft_original
    elif self.fft_type == "filtered":
      fft = trial.fft_filtered
    else:
      raise Exception("Invalid fft_type")
    return fft

  def get_select_amps(self, fft_pt):
    if self.value_type == "upper":
      return [fft_pt.upper_amp]
    elif self.value_type == "lower":
      return [fft_pt.lower_amp]
    elif self.value_type == "both":
      return [fft_pt.upper_amp, fft_pt.lower_amp]
    else:
      raise Exception("Invalid fft_type")

  def to_vec(self, trial):
    raise NotImplementedError()

  def __getitem__(self, batch_idx):
    batch_start_idx = self.batch_size * batch_idx
    batch_end_idx = batch_start_idx + self.batch_size

    data = []
    labels = []
    for d, l in self.data_labels[batch_start_idx:batch_end_idx]:
      data.append(d)
      labels.append(l)

    return data, labels

  def on_epoch_end(self):
    shuffle(self.data_labels)

  def __len__(self):
    return int(np.ceil(len(self.data_labels)/self.batch_size))

class HighestPeakSequence(BaseTrialSequence):
  def __init__(self, proto_dir_path, batch_size, num_peaks=10, **kwargs):
    self.num_peaks = num_peaks
    assert self.num_peaks > 0
    BaseTrialSequence.__init__(self, proto_dir_path, batch_size, **kwargs)

  def to_vec(self, trial):
    freq_amp = [(pt.frequency, np.mean(self.get_select_amps(pt)))
                for pt in self.get_select_fft(trial)]
    freq_amp.sort(key=lambda f_a: f_a[1], reverse=True)
    return np.array(
        [f for f, a in freq_amp[:self.num_peaks]] + \
        [a for f, a in freq_amp[:self.num_peaks]] + \
        [trial.vehicle.mass, trial.vehicle.speed] + \
        ([trial.temperature] if trial.HasField("temperature") else []))
