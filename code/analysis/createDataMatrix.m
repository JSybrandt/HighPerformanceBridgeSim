function out = createDataMatrix(dirPath, outPath)

NUM_DAYS = 360;

% I am lazy
CARS_PER_DAY = 1440;
NUM_FEATURES = 26;
NUM_OBSERVATIONS = CARS_PER_DAY * NUM_DAYS;

data = zeros(NUM_OBSERVATIONS, NUM_FEATURES);

for day = 1:NUM_DAYS
    disp(day)
    startIdx = (day-1) * CARS_PER_DAY + 1;
    endIdx = startIdx + CARS_PER_DAY - 1;
    disp(size(data(startIdx:endIdx, :)))
    data(startIdx:endIdx, :) = day2Mat(day, dirPath);
end

save(outPath, 'data')

out = 0;
