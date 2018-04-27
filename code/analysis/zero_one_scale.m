function out = createDataMatrix(inPath, outPath)

NUM_PEAKS = 10;
NUM_FEATURES = 26;

PEAK_DATA_COLS=20
MASS_COL=21;
SPEED_COL=22;
TEMP_COL=23;
RAIN_COL=24
DAY_COL = 25;
DAMAGE_CLASS_COL = 26;

% select everything
SAMPLE_RATIO = 1;
% don't subset by mass
SELECT_MASS = 0;

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

dists = maxs - mins;

% If any columns are all the same number, don't scale
% Fixes the div by zero
dists(dists==0) = 1

% don't scale class
dists(DAMAGE_CLASS_COL) = 1;
mins(DAMAGE_CLASS_COL) = 0;

% don't scale day
dists(DAY_COL) = 1;
mins(DAY_COL) = 0;

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
