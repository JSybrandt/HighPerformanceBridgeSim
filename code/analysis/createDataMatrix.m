function out = createDataMatrix(dir_path, out_path, use_other)

use_other = str2num(use_other);

files = dir(strcat(dir_path, '/*.mat'));

get_path = @(i) strcat(files(i).folder, '/',  files(i).name);

num_days = size(files, 1);

prepFile = get_path(1);

load(prepFile, 'Multiple_Vehicles', 'n')
cars_per_day = n;

total_cars = cars_per_day * num_days;

if Multiple_Vehicles && use_other
    % multiple vechicles doubles the number of cars per day
     total_cars = total_cars * 2;
end

% test build first day
pre_compiled_mat = day2Mat(prepFile, use_other);
[one_day_rows, num_features] = size(pre_compiled_mat);

% setup the overall matrix
data = zeros(total_cars, num_features);

data(1:one_day_rows, :) = pre_compiled_mat;
start_idx = one_day_rows + 1;

for file_num = 2:num_days
    tmp_mat = day2Mat(get_path(file_num), use_other);
    rows = size(tmp_mat, 1);
    data(start_idx:start_idx+rows-1, :) = tmp_mat;
    start_idx = start_idx + rows;
end

% only save rows we use
data = data(1:start_idx-1, :);

save(out_path, 'data')

out = 0;
