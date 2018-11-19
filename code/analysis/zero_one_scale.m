% In order for matlab to run as a compiled linux utility, we need a function
% signature. The parameters (input_dir_path...) are supplied as strings from the
% command line

function error_code = createDataMatrix(input_path, output_path)

  % Note that these columns should equal the python code
  % and the order of paramters in dat2mat
  NUM_PEAKS = 10;
  NUM_FEATURES = 26;
  PEAK_DATA_COLS=20
  MASS_COL=21;
  SPEED_COL=22;
  TEMP_COL=23;
  RAIN_COL=24
  DAY_COL = 25;
  DAMAGE_CLASS_COL = 26;

  load(input_path, 'data');

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

  save(output_path, 'data', 'mins', 'dists')

  error_code = 0

end % function
