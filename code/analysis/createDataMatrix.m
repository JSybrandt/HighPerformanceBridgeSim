% In order for matlab to run as a compiled linux utility, we need a function
% signature. The parameters (input_dir_path...) are supplied as strings from the
% command line

function error_code = createDataMatrix(input_dir_path, out_path, record_both_cars)

  % record_both_cars is going to be 0 or 1
  record_both_cars = str2num(record_both_cars);

  % Get all files in the input directory (output of simulation)
  % We expect one file per day
  files = dir(strcat(input_dir_path, '/*.mat'));
  num_days = size(files, 1);
  paths = @(i) strcat(files(i).folder, '/',  files(i).name);

  % We are going to use the first day to prep some info
  first_path = paths(1);
  load(first_path, 'Multiple_Vehicles', 'n')
  cars_per_day = n;

  total_cars = cars_per_day * num_days;
  if Multiple_Vehicles && record_both_cars
    % multiple vechicles doubles the number of cars per day
    % We're okay if this is an over-estimate
    total_cars = total_cars * 2;
  end

  % Load the first day,set matrix dimensions
  first_day_data = day2Mat(first_path, record_both_cars);
  [num_rows, num_features] = size(first_day_data);

  % "data" is the overall result
  data = zeros(total_cars, num_features);

  data(1:num_rows, :) = first_day_data;
  start_idx = num_rows + 1;

  for file_num = 2:num_days
      % dat2Mat extracts only the relevant feature information from a trial
      day_data = day2Mat(paths(file_num), record_both_cars);
      rows = size(day_data, 1);
      data(start_idx:start_idx+rows-1, :) = day_data;
      start_idx = start_idx + rows;
  end

  % only save rows we use
  data = data(1:start_idx-1, :);

  save(out_path, 'data')

  % no error
  error_code = 0;

end % end function
