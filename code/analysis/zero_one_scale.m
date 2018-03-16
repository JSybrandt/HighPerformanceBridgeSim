function out = createDataMatrix(inPath, outPath)

NUM_PEAKS = 10;
NUM_FEATURES = 22;

MASS_COL=21;
SPEED_COL=22;
RECORD_DAY_COL = 23;
MIN_BRIDGE_FREQ_COL = 24;
MAX_BRIDGE_FREQ_COL = 25;
DAMAGE_CLASS_COL = 26;

SAMPLE_RATIO = 1;
SELECT_MASS = 8000;

load(inPath, 'data');

%subset
if SELECT_MASS > 0
  data = data(data(:, MASS_COL)>=SELECT_MASS, :);
end

% zero-one scale
mins = min(data);
maxs = max(data);

% put all freq in the same scale
mins(1:NUM_PEAKS) = min(mins(1:NUM_PEAKS));
maxs(1:NUM_PEAKS) = max(maxs(1:NUM_PEAKS));

% put all amplitude in the same scale
mins(NUM_PEAKS+1:2*NUM_PEAKS) = min(mins(NUM_PEAKS+1:2*NUM_PEAKS));
maxs(NUM_PEAKS+1:2*NUM_PEAKS) = max(maxs(NUM_PEAKS+1:2*NUM_PEAKS));

% put bridge freqs in same scale
mins(MIN_BRIDGE_FREQ_COL:MAX_BRIDGE_FREQ_COL) = ...
  min(mins(MIN_BRIDGE_FREQ_COL:MAX_BRIDGE_FREQ_COL));
maxs(MIN_BRIDGE_FREQ_COL:MAX_BRIDGE_FREQ_COL) = ...
  max(maxs(MIN_BRIDGE_FREQ_COL:MAX_BRIDGE_FREQ_COL));

dists = maxs - mins;

% don't scale class
dists(DAMAGE_CLASS_COL) = 1;
mins(DAMAGE_CLASS_COL) = 0;

% don't scale day
dists(RECORD_DAY_COL) = 1;
mins(RECORD_DAY_COL) = 0;

data = (data - mins) ./ dists;

% shuffle
[num_examples, num_features] = size(data);
shuffleRows = randperm(num_examples);
data = data(shuffleRows, :);

% subsample
if SAMPLE_RATIO < 1
  num_examples = floor(num_examples * SAMPLE_RATIO);
  data = data(1:num_examples, :);
end

save(outPath, 'data', 'mins', 'dists')

out = 0
